"""Create workflow report."""
import json
import math

from dominate.tags import p
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
import numpy as np
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def plot_trucations(report, truncations_file):
    """Make a report section with start and end position histograms.

    The truncations_file contains start and end positions of alignments that are fully
    contained within the ITR-ITR regions.
    """
    df = pd.read_csv(
        truncations_file, sep='\t',
        dtype={
            'Read start': str,
            'Read end': np.uint32,
            'sample_id': str
        }
    )

    with report.add_section("Truncations", "Truncations"):
        p(
            "This plot illustrates the frequency of start and end positions of "
            "alignments that map completely within the transgene plasmid ITR-ITR "
            "region, helping to identify potential truncation hotspots."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    df_sample.drop(columns=['sample_id'], inplace=True)
                    plt = ezc.histplot(data=df_sample, binwidth=5)
                    plt._fig.xaxis.axis_label = 'Transgene cassette genome position'
                    plt._fig.yaxis.axis_label = 'Number of alignments'
                    EZChart(plt, theme='epi2melabs')


def plot_truncation_severity(report, severity_file):
    """Make report section with truncation severity breakdown.
    
    Shows categorization of reads by truncation severity.
    """
    df = pd.read_csv(
        severity_file,
        sep='\t',
        dtype={
            'severity': str,
            'count': np.int64,
            'percentage': np.float64,
            'sample_id': str
        }
    )

    # Define severity order for display
    severity_order = [
        "Full length (≥90%)",
        "Minor truncation (75-90%)",
        "Moderate truncation (50-75%)",
        "Severe truncation (25-50%)",
        "Very severe truncation (<25%)"
    ]

    with report.add_section("Truncation Severity", "Severity"):
        p(
            "This shows the breakdown of reads by truncation severity. "
            "'Full length' reads cover ≥90% of the expected AAV genome size."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    # Sort by severity order
                    df_sample = df_sample.copy()
                    df_sample['severity'] = pd.Categorical(
                        df_sample['severity'], categories=severity_order, ordered=True)
                    df_sample = df_sample.sort_values('severity')

                    with Grid(columns=2):
                        # Pie chart of severity distribution
                        if not df_sample.empty and df_sample['count'].sum() > 0:
                            # Filter out zero counts for pie chart
                            df_nonzero = df_sample[df_sample['count'] > 0]
                            plt = ezc.barplot(
                                data=df_nonzero, x='severity', y='percentage')
                            plt.title = dict(text='Truncation Severity Distribution')
                            plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                            plt._fig.yaxis.axis_label = 'Percentage of reads'
                            EZChart(plt, theme='epi2melabs', height='400px')

                        # Paginated table with counts
                        df_display = df_sample[['severity', 'count', 'percentage']].copy()
                        df_display.columns = ['Severity Category', 'Count', 'Percentage (%)']
                        DataTable.from_pandas(df_display, use_index=False)


def plot_truncation_length_granular(report, granular_file):
    """Make report section with granular truncation length distribution.
    
    Shows bar plot with each unique aligned length as its own bar.
    Uses truncation data (reads within ITR-ITR region).
    """
    df = pd.read_csv(
        granular_file,
        sep='\t',
        dtype={
            'aligned_length': np.int64,
            'count': np.int64,
            'percentage': np.float64,
            'sample_id': str
        }
    )

    with report.add_section("Truncation Length Distribution", "Trunc Length"):
        p(
            "This shows the distribution of aligned lengths (in base pairs) for reads "
            "within the ITR-ITR region. Each bar represents a unique aligned length, "
            "providing granular visibility into truncation patterns."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    # Sort by aligned length
                    df_sample = df_sample.sort_values('aligned_length')
                    
                    # Display total reads
                    total_reads = df_sample['count'].sum()
                    p(f"Total reads: {total_reads:,}")

                    with Grid(columns=1):
                        if not df_sample.empty:
                            plt = ezc.barplot(
                                data=df_sample, x='aligned_length', y='percentage')
                            plt.title = dict(text='Truncation Length Distribution')
                            plt._fig.xaxis.axis_label = 'Aligned length (bp)'
                            plt._fig.yaxis.axis_label = 'Percentage of reads'
                            EZChart(plt, theme='epi2melabs', height='400px')


def plot_integrity(report, summary_file, distribution_file):
    """Make report section with transgene integrity metrics.
    
    Shows %Intact summary and bar plots of mapped/continuously mapped score distributions.
    """
    df_summary = pd.read_csv(
        summary_file,
        sep='\t',
        dtype={
            'total_mapped': np.int64,
            'total_continuously_mapped': np.int64,
            'percent_intact': np.float64,
            'total_reads': np.int64,
            'sample_id': str
        }
    )
    
    df_dist = pd.read_csv(
        distribution_file,
        sep='\t',
        dtype={
            'score': np.int64,
            'count': np.int64,
            'percentage': np.float64,
            'score_type': str,
            'sample_id': str
        }
    )

    with report.add_section("Transgene Integrity", "Integrity"):
        p(
            "This shows the transgene integrity metrics based on base-level matching. "
            "%Intact = (continuously mapped bases) / (total mapped bases) × 100. "
            "Continuously mapped reads have all bases matching with no mismatches or indels."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df_summary.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    # Display summary metrics
                    row = df_sample.iloc[0]
                    p(f"**%Intact: {row['percent_intact']:.2f}%**")
                    p(f"Total mapped bases: {row['total_mapped']:,}")
                    p(f"Continuously mapped bases: {row['total_continuously_mapped']:,}")
                    p(f"Total reads analyzed: {row['total_reads']:,}")

                    # Get distribution data for this sample
                    df_dist_sample = df_dist[df_dist['sample_id'] == sample]
                    
                    with Grid(columns=2):
                        # Mapped score distribution
                        df_mapped = df_dist_sample[
                            df_dist_sample['score_type'] == 'mapped'].copy()
                        if not df_mapped.empty:
                            df_mapped = df_mapped.sort_values('score')
                            plt = ezc.barplot(
                                data=df_mapped, x='score', y='percentage')
                            plt.title = dict(text='Mapped Score Distribution')
                            plt._fig.xaxis.axis_label = 'Mapped bases per read'
                            plt._fig.yaxis.axis_label = 'Percentage of reads'
                            EZChart(plt, theme='epi2melabs', height='400px')
                        
                        # Continuously mapped score distribution
                        df_cont = df_dist_sample[
                            df_dist_sample['score_type'] == 'continuously_mapped'].copy()
                        if not df_cont.empty:
                            df_cont = df_cont.sort_values('score')
                            plt = ezc.barplot(
                                data=df_cont, x='score', y='percentage')
                            plt.title = dict(text='Continuously Mapped Score Distribution')
                            plt._fig.xaxis.axis_label = 'Continuously mapped bases per read'
                            plt._fig.yaxis.axis_label = 'Percentage of reads'
                            EZChart(plt, theme='epi2melabs', height='400px')


def plot_length_statistics(report, length_stats_file):
    """Make report section with read alignment length distribution.
    
    Shows histogram of alignment lengths in base pairs for all transgene-mapped reads.
    """
    df = pd.read_csv(
        length_stats_file,
        sep='\t',
        dtype={
            'bin': str,
            'count': np.int64,
            'percentage': np.float64,
            'sample_id': str
        }
    )

    with report.add_section("Read Length Distribution", "Length Stats"):
        p(
            "This shows the distribution of alignment lengths (in base pairs) for all reads "
            "mapping to the transgene plasmid. This provides the raw size distribution of "
            "sequenced molecules, independent of completeness calculations."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    # Separate summary stats from bin data
                    stats_rows = ['SUMMARY', 'mean_length_bp', 'median_length_bp', 
                                  'min_length_bp', 'max_length_bp']
                    df_bins = df_sample[~df_sample['bin'].isin(stats_rows)].copy()
                    df_stats = df_sample[df_sample['bin'].isin(stats_rows)].copy()

                    # Extract and display key statistics
                    if not df_stats.empty:
                        summary_row = df_stats[df_stats['bin'] == 'SUMMARY']
                        mean_row = df_stats[df_stats['bin'] == 'mean_length_bp']
                        median_row = df_stats[df_stats['bin'] == 'median_length_bp']
                        
                        if not summary_row.empty:
                            total = int(summary_row['count'].values[0])
                            p(f"Total reads: {total:,}")
                        if not mean_row.empty:
                            mean_len = int(mean_row['percentage'].values[0])
                            p(f"Mean alignment length: {mean_len:,} bp")
                        if not median_row.empty:
                            median_len = int(median_row['percentage'].values[0])
                            p(f"Median alignment length: {median_len:,} bp")

                    # Plot histogram - only bins with at least 1 read
                    with Grid(columns=2):
                        # Filter to non-empty bins only
                        df_bins_nonzero = df_bins[df_bins['count'] > 0].copy()
                        
                        if not df_bins_nonzero.empty:
                            plt = ezc.barplot(
                                data=df_bins_nonzero, x='bin', y='percentage')
                            plt.title = dict(text='Alignment Length Distribution')
                            plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                            plt._fig.xaxis.axis_label = 'Alignment length (bp)'
                            plt._fig.yaxis.axis_label = 'Percentage of reads'
                            EZChart(plt, theme='epi2melabs', height='400px')

                        # Display bin counts table (only non-empty bins)
                        df_bins_display = df_bins_nonzero[['bin', 'count', 'percentage']].copy()
                        df_bins_display.columns = ['Length Bin', 'Count', 'Percentage (%)']
                        DataTable.from_pandas(df_bins_display, use_index=False)


def plot_itr_coverage(report, coverage_file):
    """Make report section with ITR-ITR coverage of transgene cassette region."""
    df = pd.read_csv(
        coverage_file,
        sep=r"\s+",
        dtype={
            'ref': str,
            'pos': np.uint32,
            'depth': np.uint32,
            'strand': str,
            'sample_id': str
        })

    with report.add_section("ITR-ITR coverage", "Coverage"):
        p(
            "For each transgene reference, sequencing depth is calculated "
            "for both forward and reverse mapping alignments."

        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    with Grid(columns=1):
                        for ref, df_ref in df_sample.groupby('ref'):

                            plt = ezc.lineplot(
                                data=df_ref, title=ref,
                                x='pos', y='depth', hue='strand',
                                marker=False, s=2
                            )
                            EZChart(plt, theme='epi2melabs', height='300px')


def plot_contamination(report, class_counts):
    """Make report section with contamination plots.

    Two plots: (1) mapped/unmapped; (2) mapped reads per reference
    """
    df_class_counts = pd.read_csv(
        class_counts,
        sep='\t',
        dtype={
            'Reference': str,
            'Number of alignments': np.uint32,
            'Percentage of alignments': np.float32,
            'sample_id': str
        }
    )

    with report.add_section("Contamination", "Contamination"):
        p(
            "These two plots show mapping summaries that can highlight "
            "potential contamination issues."
        )
        p(
            "The first plot shows the percentage of reads that either map to any "
            "combined reference sequence or are unmapped."
        )
        p(
            "The second plot breaks down the the alignment numbers into the "
            "specific references (host, helper plasmid, Rep-Cap plasmid, and transgene "
            "plasmid)."
        )

        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df_class_counts.groupby('sample_id'):
                # Reference x-axis labels are taken from the refrence filenames.
                # Trucate these to a maximum of 25 characters.
                df_sample['Reference'] = df_sample['Reference'].apply(
                    lambda x: x if len(x) < 25 else x[:22] + '...')
                with tabs.add_dropdown_tab(sample):
                    with Grid(columns=2):
                        df_reads = df_sample[
                            df_sample.Reference.isin(['Mapped', 'Unmapped'])]
                        df_reads = df_reads.rename(columns={
                            'Percentage of alignments': 'Percentage of Reads'})
                        plt = ezc.barplot(
                            data=df_reads, x='Reference', y='Percentage of Reads')
                        plt.title = dict(text='Reads mapped/unmapped')
                        plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                        EZChart(plt, theme='epi2melabs', height='400px')

                        df_alns = df_sample[
                            ~df_sample.Reference.isin(['Mapped', 'Unmapped'])]
                        plt = ezc.barplot(
                            data=df_alns, x='Reference', y='Percentage of alignments')
                        plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                        plt.title = dict(text='Alignment counts per target')
                        EZChart(plt, theme='epi2melabs', height='400px')

def plot_recombination(report, recomb_summary_file):
    """Make report section showing inter-plasmid recombination events.
    
    Shows counts of reads with alignments to multiple different reference plasmids.
    Includes explanatory text about filtering and backbone artifacts.
    """
    df = pd.read_csv(
        recomb_summary_file,
        sep='\t',
        dtype={
            'recombination_type': str,
            'count': np.int64,
            'percentage': np.float64,
            'sample_id': str
        }
    )
    
    with report.add_section("Recombination Events", "Recombination"):
        p(
            "This section identifies potential inter-plasmid recombination events - "
            "reads with valid alignments to multiple different reference plasmids."
        )
        p(
            "Filtering applied: Transgene alignments are restricted to the ITR-ITR "
            "cassette region. Helper/RepCap plasmids require ≥200bp alignments. "
            "Host genome requires ≥500bp alignments to reduce false positives from "
            "short spurious matches."
        )
        p(
            "Note: 'helper_repcap_only' events (without transgene involvement) are "
            "typically artifacts from shared plasmid backbone sequences (e.g., AmpR, "
            "pUC origin) and usually do not indicate quality issues with the transgene "
            "product. These are reported separately from true recombination events."
        )
        
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    # Get summary stats
                    total_row = df_sample[
                        df_sample['recombination_type'] == 'total_recombination_events'
                    ]
                    backbone_row = df_sample[
                        df_sample['recombination_type'] == 'total_backbone_only_events'
                    ]
                    
                    if not total_row.empty:
                        total_pct = total_row['percentage'].values[0]
                        total_count = total_row['count'].values[0]
                        p(f"Transgene recombination events: {total_count:,} reads ({total_pct:.2f}%)")
                    
                    if not backbone_row.empty:
                        bb_pct = backbone_row['percentage'].values[0]
                        bb_count = backbone_row['count'].values[0]
                        if bb_count > 0:
                            p(f"Backbone-only events (helper+repcap, no transgene): "
                              f"{bb_count:,} reads ({bb_pct:.2f}%) - likely artifacts")
                    
                    # Filter to show recombination types only (not summary rows)
                    df_types = df_sample[
                        ~df_sample['recombination_type'].isin([
                            'total_recombination_events',
                            'total_backbone_only_events',
                            'no_inter_plasmid_recombination'
                        ])
                    ].copy()
                    
                    with Grid(columns=2):
                        # Bar plot of recombination types
                        if not df_types.empty and df_types['count'].sum() > 0:
                            df_nonzero = df_types[df_types['count'] > 0]
                            if not df_nonzero.empty:
                                # Shorten labels for display
                                df_nonzero = df_nonzero.copy()
                                df_nonzero['type_short'] = df_nonzero['recombination_type'].str.replace(
                                    '_recombination', '').str.replace('_', '+')
                                plt = ezc.barplot(
                                    data=df_nonzero, x='type_short', y='count')
                                plt.title = dict(text='Recombination Types')
                                plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                                plt._fig.yaxis.axis_label = 'Number of reads'
                                EZChart(plt, theme='epi2melabs', height='400px')
                        else:
                            # No recombination events - show the "no recombination" bar
                            no_recomb_row = df_sample[
                                df_sample['recombination_type'] == 'no_inter_plasmid_recombination'
                            ]
                            if not no_recomb_row.empty:
                                df_plot = no_recomb_row.copy()
                                df_plot['type_short'] = 'no_recombination'
                                plt = ezc.barplot(
                                    data=df_plot, x='type_short', y='count')
                                plt.title = dict(text='Recombination Types')
                                plt._fig.yaxis.axis_label = 'Number of reads'
                                EZChart(plt, theme='epi2melabs', height='400px')
                        
                        # Paginated table
                        df_display = df_sample[['recombination_type', 'count', 'percentage']].copy()
                        df_display.columns = ['Type', 'Count', 'Percentage (%)']
                        DataTable.from_pandas(df_display, use_index=False)


def plot_aav_structures(report, structures_file):
    """Make report section barplots detailing the AAV structures found."""
    df = pd.read_csv(
        structures_file,
        sep='\t',
        dtype={
            'Assigned_genome_type': str,
            'count': np.uint32,
            'percentage': np.float32,
            'sample_id': str

        })

    with report.add_section("AAV Structures", "Structures"):
        p(
            "The numbers of different of the various AAV transgene genome types  "
            "identified in the sample(s) are summarised here."
        )

        p(
            "A detailed report containing more granular genome type assignments "
            "per read can be found at:"
            " `output/<sample_id>/<sample_id>_aav_per_read_info.tsv` "
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):

                with tabs.add_dropdown_tab(sample):
                    df_sample = df_sample.sort_values('percentage', ascending=False)
                    # Plot of main genome type counts
                    plt = ezc.barplot(
                        df_sample,
                        x='Assigned_genome_type',
                        y='percentage')
                    plt.title = dict(text='Genome types')
                    plt._fig.xaxis.major_label_orientation = 45 * (math.pi / 180)
                    EZChart(plt, theme='epi2melabs')

                    # Table with counts and percentages
                    # (in lieu of being able to annotate bar plots in ezcharts)
                    df_sample = (
                        df_sample.astype({'percentage': 'float64'})
                        .round({'count': 2, 'percentage': 2})
                    )
                    DataTable.from_pandas(df_sample, use_index=False)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "AAV QC workflow report", "wf-aav-qc",
        args.params, args.versions, args.wf_version)

    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["sample"])

    with report.add_section('Read summary', 'Read summary'):
        names = tuple(d['sample'] for d in sample_details)
        stats = tuple(args.stats)
        if len(stats) == 1:
            stats = stats[0]
            names = names[0]
        fastcat.SeqSummary(stats, sample_names=names)

    plot_contamination(
        report,
        args.contam_class_counts)
    plot_trucations(report, args.truncations)
    plot_truncation_severity(report, args.truncation_severity)
    plot_truncation_length_granular(report, args.truncation_length_granular)
    plot_integrity(report, args.integrity_summary, args.integrity_distribution)
    plot_length_statistics(report, args.length_statistics)
    plot_itr_coverage(report, args.itr_coverage)
    plot_recombination(report, args.recombination_summary)
    plot_aav_structures(report, args.aav_structures)

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in sample_details:
                with tabs.add_dropdown_tab(d["sample"]):
                    df = pd.DataFrame.from_dict(
                        d, orient="index", columns=["Value"])
                    df.index.name = "Key"
                    DataTable.from_pandas(df)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='+',
        help="Fastcat per-read stats, ordered as per entries in --metadata.")
    parser.add_argument(
        "--truncations", help="TSV with start and end columns for.")
    parser.add_argument(
        "--truncation_severity", help="TSV with truncation severity summary.")
    parser.add_argument(
        "--truncation_length_granular", help="TSV with per-length truncation data.")
    parser.add_argument(
        "--length_statistics", help="TSV with read length completeness statistics.")
    parser.add_argument(
        "--itr_coverage", help="TSV with alignment Pos and EndPos columns.")


    parser.add_argument(
        "--contam_class_counts", help="TSV of reference mapping counts.")
    parser.add_argument(
        "--recombination_summary", help="TSV with recombination event summary.")
    parser.add_argument(
        "--integrity_summary", help="TSV with integrity summary metrics.")
    parser.add_argument(
        "--integrity_distribution", help="TSV with integrity score distribution.")
    parser.add_argument(
        "--aav_structures", help="TSV of reads with AAV structure assignment.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    return parser
