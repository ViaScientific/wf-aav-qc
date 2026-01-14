"""
Identify truncation hotspots.

Create a CSV file containing start and end locations for alignments fully contained
within the ITR-ITR region.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("truncations")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--itr_range', help="[itr1_start, itr_2_end]",
        nargs='*', type=int)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--outfile',
        help="Path to output",
        type=Path)
    parser.add_argument(
        '--summary_outfile',
        help="Path to output truncation severity summary TSV",
        type=Path,
        default=None)

    return parser


def assign_severity(completeness_pct: float) -> str:
    """Assign a severity category based on completeness percentage.
    
    Categories:
    - Full length: >= 90%
    - Minor truncation: 75-90%
    - Moderate truncation: 50-75%
    - Severe truncation: 25-50%
    - Very severe truncation: < 25%
    
    :param completeness_pct: Completeness percentage
    :return: Severity category string
    """
    if completeness_pct >= 90:
        return "Full length (≥90%)"
    elif completeness_pct >= 75:
        return "Minor truncation (75-90%)"
    elif completeness_pct >= 50:
        return "Moderate truncation (50-75%)"
    elif completeness_pct >= 25:
        return "Severe truncation (25-50%)"
    else:
        return "Very severe truncation (<25%)"


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


def main(args):
    """Run main entry point."""
    # Get the ITR locations
    itr1_start_pos, itr2_end_pos = args.itr_range
    expected_length = itr2_end_pos - itr1_start_pos

    loaded_dataframes = []
    # Get all the reads that map to one of the trans plasmids
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
        for df_bam in reader:
            df_bam = df_bam.loc[df_bam.Ref.isin([args.transgene_plasmid_name])].copy()

            # Filter for alignments that start and end within the ITR-ITR region
            df_bam = df_bam.loc[
                (df_bam.Pos > itr1_start_pos - 3) &
                (df_bam.EndPos < itr2_end_pos + 3)
            ].copy()

            if not df_bam.empty:
                # Keep Read, Pos, EndPos for interval merging later
                loaded_dataframes.append(df_bam[['Read', 'Pos', 'EndPos']])

    if not loaded_dataframes:
        # No reads - create empty outputs
        df_bam = pd.DataFrame({
            'Ref Start': [],
            'Ref End': [],
            'sample_id': []
        })

        df_bam.to_csv(args.outfile, sep='\t', index=False)
        
        if args.summary_outfile:
            df_summary = pd.DataFrame({
                'severity': [],
                'count': [],
                'percentage': [],
                'sample_id': []
            })
            df_summary.to_csv(args.summary_outfile, sep='\t', index=False)
        return

    # Combine all chunks
    df_all = pd.concat(loaded_dataframes, ignore_index=True)

    # Write csv of start and end positions for plotting (original output)
    df_positions = (
        df_all[['Pos', 'EndPos']]
        .rename(columns={'Pos': 'Ref Start', 'EndPos': 'Ref End'})
    )
    df_positions['sample_id'] = args.sample_id
    df_positions.to_csv(args.outfile, sep='\t', index=False)


    # Generate severity summary if output path provided
    if args.summary_outfile:
        # Define severity order for display
        severity_order = [
            "Full length (≥90%)",
            "Minor truncation (75-90%)",
            "Moderate truncation (50-75%)",
            "Severe truncation (25-50%)",
            "Very severe truncation (<25%)"
        ]
        
        # Aggregate per read using interval merging for accurate coverage
        def aggregate_read_intervals(group):
            """Aggregate multiple alignments using interval merging."""
            intervals = list(zip(group['Pos'].values, group['EndPos'].values))
            merged_coverage = calculate_merged_coverage(intervals)
            completeness_pct = (merged_coverage / expected_length * 100)
            return pd.Series({'completeness_pct': completeness_pct})
        
        df_per_read = (
            df_all.groupby('Read')
            .apply(aggregate_read_intervals, include_groups=False)
            .reset_index()
        )
        df_per_read['severity'] = df_per_read['completeness_pct'].apply(assign_severity)

        
        # Count by severity
        severity_counts = df_per_read['severity'].value_counts()
        total_reads = len(df_per_read)
        
        # Create summary dataframe with all severities (including zeros)
        df_summary = pd.DataFrame({'severity': severity_order})
        df_summary['count'] = df_summary['severity'].map(severity_counts).fillna(0).astype(int)
        df_summary['percentage'] = (df_summary['count'] / total_reads * 100).round(2)
        df_summary['sample_id'] = args.sample_id
        
        df_summary.to_csv(args.summary_outfile, sep='\t', index=False)

