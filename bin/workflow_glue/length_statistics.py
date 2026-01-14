#!/usr/bin/env python
"""Calculate read alignment length statistics for transgene-mapped reads.

This module bins reads by their alignment length (in base pairs) to the transgene
plasmid, providing a distribution of read sizes independent of completeness metrics.

This differs from truncation_severity which measures completeness percentage relative
to the ITR-ITR cassette region. This module shows the raw alignment length distribution.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


# Default bin size in base pairs
DEFAULT_BIN_SIZE = 500


def argparser():
    """Create argument parser."""
    parser = wf_parser("length_statistics")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--bin_size',
        type=int,
        default=DEFAULT_BIN_SIZE,
        help=f"Bin size in base pairs for alignment length histogram. Default: {DEFAULT_BIN_SIZE}bp")
    parser.add_argument(
        '--max_length',
        type=int,
        default=5000,
        help="Maximum length for binning (reads longer than this go in final bin). Default: 5000bp")
    parser.add_argument(
        '--outfile',
        help="Path to output TSV",
        type=Path)

    return parser


def merge_intervals(intervals):
    """Merge overlapping intervals to avoid double-counting coverage.
    
    :param intervals: List of (start, end) tuples
    :return: List of merged non-overlapping (start, end) tuples
    """
    if not intervals:
        return []
    
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    
    return merged


def calculate_merged_coverage(intervals):
    """Calculate total coverage from merged intervals."""
    merged = merge_intervals(intervals)
    return sum(end - start for start, end in merged)


def assign_length_bin(aligned_length: int, bin_size: int, max_length: int) -> str:
    """Assign an alignment length to a bin.
    
    :param aligned_length: Alignment length in base pairs
    :param bin_size: Size of each bin in base pairs
    :param max_length: Maximum length for explicit bins
    :return: Bin label string (e.g., "0-500bp", "500-1000bp", ">5000bp")
    """
    if aligned_length >= max_length:
        return f"≥{max_length}bp"
    
    lower = (aligned_length // bin_size) * bin_size
    upper = lower + bin_size
    return f"{lower}-{upper}bp"


def generate_bin_order(bin_size: int, max_length: int) -> list:
    """Generate ordered list of bin labels.
    
    :param bin_size: Size of each bin in base pairs
    :param max_length: Maximum length for explicit bins
    :return: List of bin labels in order
    """
    bins = []
    for lower in range(0, max_length, bin_size):
        upper = lower + bin_size
        bins.append(f"{lower}-{upper}bp")
    bins.append(f"≥{max_length}bp")
    return bins


def main(args):
    """Run main entry point."""
    bin_size = args.bin_size
    max_length = args.max_length
    
    # Read the per-alignment read summaries in chunks
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
            # Filter for reads mapping to transgene plasmid
            # Note: NO ITR-ITR filtering - we want ALL transgene-mapped reads
            # This gives the raw read length distribution
            df_transgene = df_chunk.loc[
                df_chunk.Ref == args.transgene_plasmid_name
            ].copy()
            
            if not df_transgene.empty:
                loaded_dataframes.append(
                    df_transgene[['Read', 'Pos', 'EndPos']]
                )
    
    if not loaded_dataframes:
        # No reads mapped to transgene - create empty output
        df_summary = pd.DataFrame({
            'bin': [],
            'count': [],
            'percentage': [],
            'sample_id': []
        })
        df_summary.to_csv(args.outfile, sep='\t', index=False)
        return
    
    # Combine all chunks
    df_all = pd.concat(loaded_dataframes, ignore_index=True)
    
    # Aggregate per read using interval merging for accurate length calculation
    def aggregate_read_intervals(group):
        """Aggregate multiple alignments of a read using interval merging."""
        intervals = list(zip(group['Pos'].values, group['EndPos'].values))
        merged_coverage = calculate_merged_coverage(intervals)
        return pd.Series({'aligned_length': merged_coverage})
    
    df_per_read = (
        df_all.groupby('Read')
        .apply(aggregate_read_intervals, include_groups=False)
        .reset_index()
    )
    
    # Assign bins based on merged alignment length
    df_per_read['bin'] = df_per_read['aligned_length'].apply(
        lambda x: assign_length_bin(x, bin_size, max_length)
    )
    
    # Generate bin order for proper sorting
    bin_order = generate_bin_order(bin_size, max_length)
    
    # Count reads per bin
    bin_counts = df_per_read['bin'].value_counts()
    total_reads = len(df_per_read)
    
    # Create summary dataframe with all bins (including zeros)
    df_summary = pd.DataFrame({'bin': bin_order})
    df_summary['count'] = df_summary['bin'].map(bin_counts).fillna(0).astype(int)
    df_summary['percentage'] = (df_summary['count'] / total_reads * 100).round(2)
    df_summary['sample_id'] = args.sample_id
    
    # Add summary statistics
    mean_length = df_per_read['aligned_length'].mean()
    median_length = df_per_read['aligned_length'].median()
    min_length = df_per_read['aligned_length'].min()
    max_len = df_per_read['aligned_length'].max()
    
    stats_rows = pd.DataFrame({
        'bin': ['SUMMARY', 'mean_length_bp', 'median_length_bp', 'min_length_bp', 'max_length_bp'],
        'count': [total_reads, 0, 0, 0, 0],
        'percentage': [100.0, round(mean_length, 0), round(median_length, 0), min_length, max_len],
        'sample_id': [args.sample_id] * 5
    })
    
    df_output = pd.concat([df_summary, stats_rows], ignore_index=True)
    df_output.to_csv(args.outfile, sep='\t', index=False)
