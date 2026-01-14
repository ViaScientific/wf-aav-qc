"""Test the length statistics module."""
from pathlib import Path
from types import SimpleNamespace
import tempfile

import pandas as pd
import pytest

from workflow_glue.length_statistics import (
    calculate_completeness,
    assign_bin,
    main
)


class TestCalculateCompleteness:
    """Test the completeness calculation function."""

    def test_full_length_read(self):
        """Test a read that covers 100% of expected length."""
        assert calculate_completeness(1000, 1000) == 100.0

    def test_half_length_read(self):
        """Test a read that covers 50% of expected length."""
        assert calculate_completeness(500, 1000) == 50.0

    def test_longer_than_expected(self):
        """Test a read longer than expected (>100%)."""
        assert calculate_completeness(1200, 1000) == 120.0

    def test_zero_aligned_length(self):
        """Test zero aligned length."""
        assert calculate_completeness(0, 1000) == 0.0

    def test_zero_expected_length(self):
        """Test zero expected length (edge case)."""
        assert calculate_completeness(1000, 0) == 0.0


class TestAssignBin:
    """Test the bin assignment function."""

    def test_bin_0_10(self):
        """Test 0-10% bin assignment."""
        assert assign_bin(5.0) == "0-10%"
        assert assign_bin(0.0) == "0-10%"
        assert assign_bin(9.9) == "0-10%"

    def test_bin_50_60(self):
        """Test 50-60% bin assignment."""
        assert assign_bin(50.0) == "50-60%"
        assert assign_bin(55.5) == "50-60%"
        assert assign_bin(59.9) == "50-60%"

    def test_bin_90_100(self):
        """Test 90-100% bin assignment."""
        assert assign_bin(90.0) == "90-100%"
        assert assign_bin(99.9) == "90-100%"

    def test_bin_100_plus(self):
        """Test 100%+ bin assignment."""
        assert assign_bin(100.0) == "100%+"
        assert assign_bin(150.0) == "100%+"

    def test_boundary_values(self):
        """Test boundary values between bins."""
        assert assign_bin(10.0) == "10-20%"
        assert assign_bin(20.0) == "20-30%"
        assert assign_bin(30.0) == "30-40%"


class TestMain:
    """Test the main function end-to-end."""

    def test_basic_calculation(self):
        """Test basic length statistics calculation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create test bam_info file
            bam_info = tmpdir / "bam_info.tsv"
            test_data = pd.DataFrame({
                'Read': ['read1', 'read2', 'read3', 'read4', 'read5'],
                'Ref': ['transgene'] * 5,
                'Pos': [0, 0, 0, 0, 0],
                'EndPos': [1000, 500, 250, 100, 1100],  # 100%, 50%, 25%, 10%, 110%
                'ReadLen': [1000, 500, 250, 100, 1100]
            })
            test_data.to_csv(bam_info, sep='\t', index=False)
            
            outfile = tmpdir / "length_stats.tsv"
            
            args = SimpleNamespace(
                bam_info=bam_info,
                itr_range=[0, 1000],  # Expected length = 1000
                transgene_plasmid_name='transgene',
                sample_id='test_sample',
                outfile=outfile
            )
            
            main(args)
            
            # Read and verify output
            result = pd.read_csv(outfile, sep='\t')
            
            # Check that we have the expected bins
            assert '0-10%' in result['bin'].values
            assert '100%+' in result['bin'].values
            
            # Check mean_completeness is present
            assert 'mean_completeness_%' in result['bin'].values
            
            # Verify count for 100%+ bin (reads 1 and 5)
            count_100_plus = result[result['bin'] == '100%+']['count'].values[0]
            assert count_100_plus == 2

    def test_no_transgene_reads(self):
        """Test handling when no reads map to transgene."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create test bam_info file with reads to different ref
            bam_info = tmpdir / "bam_info.tsv"
            test_data = pd.DataFrame({
                'Read': ['read1', 'read2'],
                'Ref': ['other_ref', 'other_ref'],
                'Pos': [0, 0],
                'EndPos': [500, 500],
                'ReadLen': [500, 500]
            })
            test_data.to_csv(bam_info, sep='\t', index=False)
            
            outfile = tmpdir / "length_stats.tsv"
            
            args = SimpleNamespace(
                bam_info=bam_info,
                itr_range=[0, 1000],
                transgene_plasmid_name='transgene',
                sample_id='test_sample',
                outfile=outfile
            )
            
            main(args)
            
            # Should create empty output
            result = pd.read_csv(outfile, sep='\t')
            assert result.empty or result['count'].sum() == 0
