#!/usr/bin/env python
"""Detect inter-plasmid recombination events.

This module identifies reads with alignments to multiple different reference plasmids,
indicating potential recombination events between transgene, helper, and rep-cap plasmids.

Key features:
- Uses ref_ids.json for accurate reference categorization (not keyword matching)
- Filters transgene alignments to ITR-ITR cassette region only
- Applies minimum alignment length thresholds to reduce false positives
- Uses stricter threshold for host genome due to large size
- Classifies recombination by biological severity
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


# Default minimum alignment lengths (can be overridden via CLI)
DEFAULT_MIN_ALIGNMENT_LENGTH = 200  # For helper/rep-cap plasmids
DEFAULT_MIN_HOST_ALIGNMENT_LENGTH = 500  # Stricter for large host genomes


def argparser():
    """Create argument parser."""
    parser = wf_parser("recombination")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam containing alignment information",
        type=Path)
    parser.add_argument(
        '--ref_ids_json',
        help="JSON file with reference ID mappings (category -> ref names)",
        type=Path)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
    parser.add_argument(
        '--itr_range',
        help="ITR1 start and ITR2 end positions [itr1_start, itr2_end] for transgene filtering",
        nargs='*',
        type=int,
        default=None)
    parser.add_argument(
        '--min_alignment_length',
        type=int,
        default=DEFAULT_MIN_ALIGNMENT_LENGTH,
        help=f"Minimum alignment length (bp) for helper/rep-cap plasmids to count as "
             f"evidence of reference presence. Shorter alignments may be spurious. "
             f"Default: {DEFAULT_MIN_ALIGNMENT_LENGTH}bp. Adjust based on sequencing "
             f"technology (lower for high-accuracy reads like PacBio HiFi).")
    parser.add_argument(
        '--min_host_alignment_length',
        type=int,
        default=DEFAULT_MIN_HOST_ALIGNMENT_LENGTH,
        help=f"Minimum alignment length (bp) for host genome alignments. Higher threshold "
             f"than plasmids because large genomes (e.g., 3GB human) produce many short "
             f"random matches by chance. Default: {DEFAULT_MIN_HOST_ALIGNMENT_LENGTH}bp.")
    parser.add_argument(
        '--sample_id',
        help="Sample ID for output labeling")
    parser.add_argument(
        '--outfile',
        help="Path to output recombination events TSV",
        type=Path)
    parser.add_argument(
        '--summary_outfile',
        help="Path to output recombination summary TSV",
        type=Path,
        default=None)

    return parser


def classify_recombination_type(
    ref_categories_present: dict,
    has_valid_transgene: bool
) -> str:
    """Classify recombination type based on which reference categories are present.
    
    Uses priority-based classification - returns the most severe classification:
    1. Host integration (highest risk)
    2. Transgene + plasmid recombination
    3. Helper/RepCap only (likely backbone, lowest priority)
    
    :param ref_categories_present: Dict of category -> bool for each ref type
    :param has_valid_transgene: Whether there's a valid transgene alignment (in ITR-ITR)
    :return: Recombination type classification string
    """
    has_host = ref_categories_present.get('host', False)
    has_helper = ref_categories_present.get('helper', False)
    has_repcap = ref_categories_present.get('repcap', False)
    
    # Priority 1: Host genome integration (most severe)
    if has_valid_transgene and has_host:
        if has_helper or has_repcap:
            return "transgene_host_plasmid_recombination"
        else:
            return "transgene_host_recombination"
    
    # Priority 2: Transgene + other plasmid recombination
    if has_valid_transgene:
        if has_helper and has_repcap:
            return "transgene_helper_repcap_recombination"
        elif has_helper:
            return "transgene_helper_recombination"
        elif has_repcap:
            return "transgene_repcap_recombination"
        # Has valid transgene but no other valid alignments
        return "transgene_only"
    
    # Priority 3: No valid transgene alignment
    # Check if it's helper + repcap only (likely backbone artifact)
    if has_helper and has_repcap:
        return "helper_repcap_only"
    elif has_helper:
        return "helper_only"
    elif has_repcap:
        return "repcap_only"
    elif has_host:
        return "host_only"
    
    return "no_valid_alignments"


def main(args):
    """Run main entry point."""
    import json
    
    # Load reference ID mappings from JSON
    # This gives us authoritative category assignments (not keyword guessing)
    ref_to_category = {}
    if args.ref_ids_json and args.ref_ids_json.exists():
        with open(args.ref_ids_json) as f:
            ref_ids = json.load(f)
            # Build reverse mapping: ref_name -> category
            for category, ref_names in ref_ids.items():
                if isinstance(ref_names, list):
                    for name in ref_names:
                        ref_to_category[name] = category.lower()
                else:
                    ref_to_category[ref_names] = category.lower()
    
    # Add transgene to mapping if not already present
    if args.transgene_plasmid_name and args.transgene_plasmid_name not in ref_to_category:
        ref_to_category[args.transgene_plasmid_name] = 'transgene'
    
    # Parse ITR range for transgene filtering
    itr1_start, itr2_end = None, None
    if args.itr_range and len(args.itr_range) >= 2:
        itr1_start, itr2_end = args.itr_range[0], args.itr_range[1]
    
    # Read the per-alignment read summaries
    loaded_dataframes = []
    
    with pd.read_csv(
        args.bam_info,
        sep='\t',
        usecols=['Read', 'Ref', 'Pos', 'EndPos'],
        dtype={
            'Read': str,
            'Ref': str,
            'Pos': np.uint32,
            'EndPos': np.uint32
        },
        chunksize=50000
    ) as reader:
        for df_chunk in reader:
            # Calculate alignment length
            df_chunk['AlnLen'] = df_chunk['EndPos'] - df_chunk['Pos']
            
            # Map reference names to categories
            df_chunk['Category'] = df_chunk['Ref'].map(ref_to_category).fillna('unknown')
            
            loaded_dataframes.append(df_chunk)
    
    if not loaded_dataframes:
        _write_empty_outputs(args)
        return

    # Combine all chunks
    df_all = pd.concat(loaded_dataframes, ignore_index=True)
    
    # Apply filtering based on reference category
    def is_valid_alignment(row) -> bool:
        """Check if alignment passes quality filters for its category."""
        category = row['Category']
        aln_len = row['AlnLen']
        pos = row['Pos']
        end_pos = row['EndPos']
        
        if category == 'transgene':
            # For transgene: must be within ITR-ITR cassette region
            # Using 3bp tolerance as in the original truncations.py
            if itr1_start is not None and itr2_end is not None:
                in_cassette = (pos > itr1_start - 3) and (end_pos < itr2_end + 3)
                return in_cassette
            # If no ITR range provided, accept any transgene alignment
            return True
        
        elif category == 'host':
            # Stricter threshold for host genome
            return aln_len >= args.min_host_alignment_length
        
        elif category in ('helper', 'repcap'):
            # Standard threshold for plasmids
            return aln_len >= args.min_alignment_length
        
        else:
            # Unknown category - apply standard threshold
            return aln_len >= args.min_alignment_length
    
    # Apply filtering
    df_all['IsValid'] = df_all.apply(is_valid_alignment, axis=1)
    df_valid = df_all[df_all['IsValid']].copy()
    
    # Group by read and aggregate categories
    def aggregate_read_categories(group):
        """Aggregate valid alignments per read."""
        categories_present = set(group['Category'].unique())
        refs_involved = set(group['Ref'].unique())
        
        # Check which categories have valid alignments
        has_transgene = 'transgene' in categories_present
        has_host = 'host' in categories_present
        has_helper = 'helper' in categories_present
        has_repcap = 'repcap' in categories_present
        
        return pd.Series({
            'categories_present': categories_present,
            'refs_involved': refs_involved,
            'has_transgene': has_transgene,
            'has_host': has_host,
            'has_helper': has_helper,
            'has_repcap': has_repcap,
            'num_categories': len(categories_present)
        })
    
    df_per_read = (
        df_valid.groupby('Read')
        .apply(aggregate_read_categories, include_groups=False)
        .reset_index()
    )
    
    # Classify recombination type for each read
    def classify_read(row):
        ref_cats = {
            'host': row['has_host'],
            'helper': row['has_helper'],
            'repcap': row['has_repcap']
        }
        return classify_recombination_type(ref_cats, row['has_transgene'])
    
    df_per_read['recombination_type'] = df_per_read.apply(classify_read, axis=1)
    
    # Filter for multi-category reads (potential recombination)
    df_multi = df_per_read[df_per_read['num_categories'] > 1].copy()
    
    # Identify recombination events (exclude single-category and "only" types)
    excluded_types = [
        'transgene_only', 'host_only', 'helper_only', 
        'repcap_only', 'no_valid_alignments'
    ]
    df_recomb = df_multi[~df_multi['recombination_type'].isin(excluded_types)].copy()
    
    # Prepare output events dataframe
    df_events = pd.DataFrame({
        'Read': df_recomb['Read'],
        'refs_involved': df_recomb['refs_involved'].apply(lambda x: ','.join(sorted(x))),
        'categories_involved': df_recomb['categories_present'].apply(lambda x: ','.join(sorted(x))),
        'recombination_type': df_recomb['recombination_type'],
        'sample_id': args.sample_id
    })
    df_events.to_csv(args.outfile, sep='\t', index=False)
    
    # Generate summary if output path provided
    if args.summary_outfile:
        total_reads = len(df_per_read)
        
        # Define recombination type order (by severity - most severe first)
        type_order = [
            "transgene_host_recombination",
            "transgene_host_plasmid_recombination",
            "transgene_helper_recombination",
            "transgene_repcap_recombination",
            "transgene_helper_repcap_recombination",
            "helper_repcap_only",  # Likely backbone artifact
        ]
        
        # Count by recombination type
        type_counts = df_recomb['recombination_type'].value_counts()
        
        # Count backbone events separately
        backbone_count = type_counts.get('helper_repcap_only', 0)
        true_recomb_count = len(df_recomb) - backbone_count
        
        # Create summary dataframe
        summary_records = []
        for rtype in type_order:
            count = type_counts.get(rtype, 0)
            pct = (count / total_reads * 100) if total_reads > 0 else 0
            summary_records.append({
                'recombination_type': rtype,
                'count': count,
                'percentage': round(pct, 4),
                'sample_id': args.sample_id
            })
        
        # Add summary rows
        summary_records.append({
            'recombination_type': 'total_recombination_events',
            'count': true_recomb_count,
            'percentage': round((true_recomb_count / total_reads * 100) if total_reads > 0 else 0, 4),
            'sample_id': args.sample_id
        })
        
        summary_records.append({
            'recombination_type': 'total_backbone_only_events',
            'count': backbone_count,
            'percentage': round((backbone_count / total_reads * 100) if total_reads > 0 else 0, 4),
            'sample_id': args.sample_id
        })
        
        no_recomb = total_reads - len(df_recomb)
        summary_records.append({
            'recombination_type': 'no_inter_plasmid_recombination',
            'count': no_recomb,
            'percentage': round((no_recomb / total_reads * 100) if total_reads > 0 else 0, 4),
            'sample_id': args.sample_id
        })
        
        df_summary = pd.DataFrame(summary_records)
        df_summary.to_csv(args.summary_outfile, sep='\t', index=False)


def _write_empty_outputs(args):
    """Write empty output files when no data is available."""
    df_events = pd.DataFrame({
        'Read': [],
        'refs_involved': [],
        'categories_involved': [],
        'recombination_type': [],
        'sample_id': []
    })
    df_events.to_csv(args.outfile, sep='\t', index=False)
    
    if args.summary_outfile:
        df_summary = pd.DataFrame({
            'recombination_type': [],
            'count': [],
            'percentage': [],
            'sample_id': []
        })
        df_summary.to_csv(args.summary_outfile, sep='\t', index=False)
