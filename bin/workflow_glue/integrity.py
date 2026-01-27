"""Calculate transgene integrity scores.

Calculates per-read mapped and continuously mapped scores from BAM file.
- Full-length: Read starts within/before ITR1 end AND ends within/after ITR2 start
- Mapped score: Inner region length for full-length; calculated length for partial
- Continuously mapped score: Inner region length for full-length; 0 for partial
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
        '--itr_locations',
        help="ITR locations: itr1_start itr1_end itr2_start itr2_end",
        nargs=4, type=int)
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


def get_alignment_score(alignment):
    """Get alignment score for ranking alignments.
    
    Uses AS tag if available, otherwise mapping quality.
    
    :param alignment: pysam AlignedSegment
    :return: alignment score (higher is better)
    """
    try:
        return alignment.get_tag('AS')
    except KeyError:
        return alignment.mapping_quality


def calculate_read_integrity(read_start, read_end, itr1_end, itr2_start, inner_length):
    """Calculate integrity scores for a single read based on position.
    
    Full-length: read starts at or before ITR1 end AND ends at or after ITR2 start.
    
    :param read_start: Alignment reference start position
    :param read_end: Alignment reference end position
    :param itr1_end: End position of ITR1
    :param itr2_start: Start position of ITR2
    :param inner_length: Pre-calculated inner region length (itr2_start - itr1_end)
    :return: tuple (mapped_score, continuously_mapped_score, is_full_length)
    """
    # Check if full-length: starts within/before ITR1 and ends within/after ITR2
    is_full_length = (read_start <= itr1_end) and (read_end >= itr2_start)
    
    if is_full_length:
        # Full-length reads get the constant inner_length for both scores
        return inner_length, inner_length, True
    else:
        # Partial reads: calculate length excluding ITR regions
        # Clamp positions to inner region (itr1_end to itr2_start)
        clamped_start = max(read_start, itr1_end)
        clamped_end = min(read_end, itr2_start)
        
        # If the read doesn't overlap the inner region at all
        if clamped_start >= clamped_end:
            mapped_score = 0
        else:
            mapped_score = clamped_end - clamped_start
        
        # Partial reads get 0 for continuously mapped
        return mapped_score, 0, False


def main(args):
    """Run main entry point."""
    itr1_start, itr1_end, itr2_start, itr2_end = args.itr_locations
    
    # Calculate the inner region length (excluding ITRs)
    # This is the constant score for full-length reads
    inner_length = itr2_start - itr1_end
    
    # Open BAM file
    bam = pysam.AlignmentFile(args.bam, "rb")
    
    # Track best alignment per read (by alignment score)
    # read_id -> {'start': int, 'end': int, 'score': int}
    best_alignments = {}
    
    # Fetch alignments to transgene plasmid in the full ITR-ITR region
    try:
        alignments = bam.fetch(args.transgene_plasmid_name, itr1_start, itr2_end)
    except ValueError:
        # Reference not in BAM or region out of bounds
        alignments = []
    
    for aln in alignments:
        if aln.is_unmapped:
            continue
            
        read_id = aln.query_name
        aln_score = get_alignment_score(aln)
        read_start = aln.reference_start
        read_end = aln.reference_end
        
        # Filter: only include reads CONTAINED within ITR-ITR region
        # (same filter as truncations.py - exclude reads extending beyond ITRs)
        if read_start <= itr1_start - 3 or read_end >= itr2_end + 3:
            continue
        
        # Keep only the best alignment per read
        if read_id not in best_alignments:
            best_alignments[read_id] = {
                'start': read_start,
                'end': read_end,
                'score': aln_score
            }
        elif aln_score > best_alignments[read_id]['score']:
            best_alignments[read_id] = {
                'start': read_start,
                'end': read_end,
                'score': aln_score
            }
    
    bam.close()
    
    # Calculate integrity scores for each read's best alignment
    per_read_data = []
    for read_id, aln_info in best_alignments.items():
        mapped, continuously_mapped, is_full_length = calculate_read_integrity(
            aln_info['start'], aln_info['end'],
            itr1_end, itr2_start, inner_length
        )
        
        per_read_data.append({
            'read_id': read_id,
            'mapped_score': mapped,
            'continuously_mapped_score': continuously_mapped,
            'is_full_length': is_full_length,
            'sample_id': args.sample_id
        })
    
    if per_read_data:
        df_per_read = pd.DataFrame(per_read_data)
    else:
        df_per_read = pd.DataFrame({
            'read_id': [],
            'mapped_score': [],
            'continuously_mapped_score': [],
            'is_full_length': [],
            'sample_id': []
        })
    
    # Write per-read output
    df_per_read.to_csv(args.per_read_outfile, sep='\t', index=False)
    
    # Calculate summary
    total_mapped = df_per_read['mapped_score'].sum()
    total_continuously_mapped = df_per_read['continuously_mapped_score'].sum()
    full_length_count = df_per_read['is_full_length'].sum() if not df_per_read.empty else 0
    
    if total_mapped > 0:
        percent_intact = round(total_continuously_mapped / total_mapped * 100, 2)
    else:
        percent_intact = 0.0
    
    df_summary = pd.DataFrame({
        'total_mapped': [int(total_mapped)],
        'total_continuously_mapped': [int(total_continuously_mapped)],
        'percent_intact': [percent_intact],
        'total_reads': [len(df_per_read)],
        'full_length_reads': [int(full_length_count)],
        'inner_region_length': [inner_length],
        'sample_id': [args.sample_id]
    })
    
    df_summary.to_csv(args.summary_outfile, sep='\t', index=False)
    
    # Create distribution data for line plots
    # Mapped score distribution (frequency per score value)
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
    
    # Output distribution file
    dist_outfile = args.summary_outfile.parent / (
        args.summary_outfile.stem + '_distribution.tsv')
    df_dist.to_csv(dist_outfile, sep='\t', index=False)
