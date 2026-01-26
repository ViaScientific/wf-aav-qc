"""Calculate transgene integrity scores.

Calculates per-read mapped and continuously mapped scores from BAM file.
- Mapped score: Total matching bases within ITR region
- Continuously mapped score: Full score if all bases match, else 0
- %Intact = sum(continuously_mapped) / sum(mapped) * 100
"""

from pathlib import Path

import pandas as pd
import pysam

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("integrity")

    parser.add_argument(
        '--bam',
        help="Path to input BAM file (must be indexed)",
        type=Path)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
    parser.add_argument(
        '--transgene_plasmid_fasta',
        help="Path to transgene plasmid reference FASTA",
        type=Path)
    parser.add_argument(
        '--itr_range', help="[itr1_start, itr_2_end]",
        nargs='*', type=int)
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--per_read_outfile',
        help="Path to output per-read scores TSV",
        type=Path)
    parser.add_argument(
        '--summary_outfile',
        help="Path to output summary TSV",
        type=Path)

    return parser


def calculate_read_integrity(alignment, itr_start, itr_end, ref_fasta, ref_name):
    """Calculate integrity scores for a single alignment.
    
    :param alignment: pysam AlignedSegment
    :param itr_start: ITR region start position
    :param itr_end: ITR region end position
    :param ref_fasta: pysam.FastaFile object
    :param ref_name: Reference sequence name
    :return: tuple (mapped_score, continuously_mapped_score)
    """
    if alignment.is_unmapped:
        return 0, 0
    
    # Get aligned pairs (query_pos, ref_pos)
    # with_seq=False beucase MD tag might be missing
    aligned_pairs = alignment.get_aligned_pairs(with_seq=False)
    
    if not aligned_pairs:
        return 0, 0
    
    matches = 0
    has_mismatch_or_indel = False
    query_seq = alignment.query_sequence
    
    if query_seq is None:
        return 0, 0
        
    try:
        ref_seq_str = ref_fasta.fetch(ref_name)
    except KeyError:
        # Reference name not found in FASTA
        return 0, 0
    
    for query_pos, ref_pos in aligned_pairs:
        # Skip positions outside ITR region
        if ref_pos is None:
            # Insertion in query (no reference position)
            # Indels count against continuously mapped score
            if itr_start <= alignment.reference_start < itr_end:
                 # Approximating location check for insertion
                 # Ideally check if insertion is within ITR bounds based on adjacent matches
                 # But simplistic check: if we are processing this read which overlaps ITR...
                 # More accurate: check previous/next ref_pos.
                 # For now, simplistic: if alignment overlaps ITR, indels matter.
                 pass
            has_mismatch_or_indel = True
            continue
        
        if ref_pos < itr_start or ref_pos >= itr_end:
            continue
        
        if query_pos is None:
            # Deletion in query (no query position)
            has_mismatch_or_indel = True
            continue
        
        # Both positions exist - check if it's a match
        query_base = query_seq[query_pos]
        # ref_pos is 0-based
        if ref_pos < len(ref_seq_str):
            ref_base = ref_seq_str[ref_pos]
            if ref_base.upper() == query_base.upper():
                matches += 1
            else:
                has_mismatch_or_indel = True
        else:
            has_mismatch_or_indel = True
    
    mapped_score = matches
    continuously_mapped_score = matches if not has_mismatch_or_indel else 0
    
    return mapped_score, continuously_mapped_score


def main(args):
    """Run main entry point."""
    itr_start, itr_end = args.itr_range
    
    # Open BAM file
    bam = pysam.AlignmentFile(args.bam, "rb")
    ref_fasta = pysam.FastaFile(args.transgene_plasmid_fasta)
    
    # Track per-read scores (aggregate multiple alignments per read)
    read_scores = {}  # read_id -> {'mapped': int, 'continuously_mapped': int}
    
    # Fetch alignments to transgene plasmid in ITR region
    try:
        alignments = bam.fetch(args.transgene_plasmid_name, itr_start, itr_end)
    except ValueError:
        # Reference not in BAM or region out of bounds
        alignments = []
    
    for aln in alignments:
        read_id = aln.query_name
        mapped, continuously_mapped = calculate_read_integrity(
            aln, itr_start, itr_end, ref_fasta, args.transgene_plasmid_name)
        
        if read_id not in read_scores:
            read_scores[read_id] = {
                'mapped': 0,
                'continuously_mapped': 0,
                'has_mismatch': False
            }
        
        read_scores[read_id]['mapped'] += mapped
        # If any alignment for this read has mismatch, mark it
        if mapped > 0 and continuously_mapped == 0:
            read_scores[read_id]['has_mismatch'] = True
        elif continuously_mapped > 0:
            read_scores[read_id]['continuously_mapped'] += continuously_mapped
    
    bam.close()
    ref_fasta.close()
    
    # Create per-read dataframe
    per_read_data = []
    for read_id, scores in read_scores.items():
        # If any alignment had mismatch, continuously_mapped should be 0
        cont_mapped = 0 if scores['has_mismatch'] else scores['continuously_mapped']
        per_read_data.append({
            'read_id': read_id,
            'mapped_score': scores['mapped'],
            'continuously_mapped_score': cont_mapped,
            'sample_id': args.sample_id
        })
    
    if per_read_data:
        df_per_read = pd.DataFrame(per_read_data)
    else:
        df_per_read = pd.DataFrame({
            'read_id': [],
            'mapped_score': [],
            'continuously_mapped_score': [],
            'sample_id': []
        })
    
    # Write per-read output
    df_per_read.to_csv(args.per_read_outfile, sep='\t', index=False)
    
    # Calculate summary
    total_mapped = df_per_read['mapped_score'].sum()
    total_continuously_mapped = df_per_read['continuously_mapped_score'].sum()
    
    if total_mapped > 0:
        percent_intact = round(total_continuously_mapped / total_mapped * 100, 2)
    else:
        percent_intact = 0.0
    
    df_summary = pd.DataFrame({
        'total_mapped': [int(total_mapped)],
        'total_continuously_mapped': [int(total_continuously_mapped)],
        'percent_intact': [percent_intact],
        'total_reads': [len(df_per_read)],
        'sample_id': [args.sample_id]
    })
    
    df_summary.to_csv(args.summary_outfile, sep='\t', index=False)
    
    # Also create distribution data for bar plots
    # Mapped score distribution
    if not df_per_read.empty:
        mapped_counts = df_per_read['mapped_score'].value_counts().sort_index()
        total = len(df_per_read)
        
        df_mapped_dist = pd.DataFrame({
            'score': mapped_counts.index.astype(int),
            'count': mapped_counts.values,
            'percentage': (mapped_counts.values / total * 100).round(2),
            'score_type': 'mapped',
            'sample_id': args.sample_id
        })
        
        cont_mapped_counts = df_per_read['continuously_mapped_score'].value_counts().sort_index()
        df_cont_dist = pd.DataFrame({
            'score': cont_mapped_counts.index.astype(int),
            'count': cont_mapped_counts.values,
            'percentage': (cont_mapped_counts.values / total * 100).round(2),
            'score_type': 'continuously_mapped',
            'sample_id': args.sample_id
        })
        
        df_dist = pd.concat([df_mapped_dist, df_cont_dist], ignore_index=True)
    else:
        df_dist = pd.DataFrame({
            'score': [],
            'count': [],
            'percentage': [],
            'score_type': [],
            'sample_id': []
        })
    
    # Append distribution to summary file (or output separately)
    dist_outfile = args.summary_outfile.parent / (
        args.summary_outfile.stem + '_distribution.tsv')
    df_dist.to_csv(dist_outfile, sep='\t', index=False)
