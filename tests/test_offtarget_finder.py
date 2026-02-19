import pandas as pd
import pytest

from plexus.blast.offtarget_finder import AmpliconFinder, Position


@pytest.fixture
def sample_bound_df():
    """
    Simulate predicted_bound hits for primers.
    qseqid sseqid sstart sstrand
    """
    cols = ["qseqid", "sseqid", "sstart", "sstrand"]
    data = [
        # Target A: Intended amplicons on chr7
        ["A_F", "chr7", 1000, "plus"],
        ["A_R", "chr7", 1200, "minus"],
        # Target B: Intended amplicons on chr22
        ["B_F", "chr22", 5000, "plus"],
        ["B_R", "chr22", 5150, "minus"],
        # Off-target: A_F also binds to chr22 near B_R
        ["A_F", "chr22", 4950, "plus"],
        # Very distant hit on chr7
        ["A_R", "chr7", 10000, "minus"],
    ]
    return pd.DataFrame(data, columns=cols)


def test_find_amplicons_basic(sample_bound_df):
    finder = AmpliconFinder(sample_bound_df)
    finder.find_amplicons(max_size_bp=1000)

    df = finder.amplicon_df
    # Expected amplicons:
    # 1. A_F (1000) + A_R (1200) on chr7
    # 2. B_F (5000) + B_R (5150) on chr22
    # 3. A_F (4950) + B_R (5150) on chr22 (Cross-target/Off-target)

    assert len(df) == 3

    chr7_amps = df[df["chrom"] == "chr7"]
    assert len(chr7_amps) == 1
    assert chr7_amps.iloc[0]["product_bp"] == 201

    chr22_amps = df[df["chrom"] == "chr22"]
    assert len(chr22_amps) == 2


def test_find_amplicons_max_size(sample_bound_df):
    finder = AmpliconFinder(sample_bound_df)
    # Distance between A_F (1000) and A_R (10000) is 9000
    finder.find_amplicons(max_size_bp=5000)

    df = finder.amplicon_df
    # Still 3, because the 9000bp one is excluded
    assert len(df) == 3
    assert not any(df["product_bp"] > 5000)


def test_annotate_expected_binding(sample_bound_df):
    finder = AmpliconFinder(sample_bound_df)
    finder.find_amplicons()

    # Define expected binding positions for primers
    # primer_name, chrom, start (0-based in finder.expected_dt, but finder adds 1)
    # wait, AmpliconFinder._get_expected_dict does row["start"] + 1
    # and _annotate_expected_binding compares Position(chrom, F_start) == expected_dt[F_primer]
    # In find_amplicons: F_start = int(F_df["sstart"])

    primer_data = [
        {"primer_name": "A_F", "chrom": "chr7", "start": 999},  # 999 + 1 = 1000
        {"primer_name": "A_R", "chrom": "chr7", "start": 1199},
        {"primer_name": "B_F", "chrom": "chr22", "start": 4999},
        {"primer_name": "B_R", "chrom": "chr22", "start": 5149},
    ]
    primer_df = pd.DataFrame(primer_data)

    finder.create_ontarget_dataframe(primer_df)
    finder.create_offtarget_dataframe()

    # On-target should be A_F/A_R on chr7 and B_F/B_R on chr22
    assert len(finder.ontarget_df) == 2

    # Off-target should be A_F (4950) / B_R (5150) on chr22
    assert len(finder.offtarget_df) == 1
    assert finder.offtarget_df.iloc[0]["F_primer"] == "A_F"
    assert finder.offtarget_df.iloc[0]["chrom"] == "chr22"


def test_amplicon_size_inclusive_coordinates():
    """
    Test for BUG-02 fix: product_bp should account for inclusive coordinates.

    When forward primer is at position 1000 and reverse primer at 1200,
    the amplicon spans 1000-1200 inclusive = 201 bp, not 200 bp.
    """
    cols = ["qseqid", "sseqid", "sstart", "sstrand"]
    data = [
        ["FWD", "chr1", 1000, "plus"],
        ["REV", "chr1", 1200, "minus"],
    ]
    bound_df = pd.DataFrame(data, columns=cols)

    finder = AmpliconFinder(bound_df)
    finder.find_amplicons(max_size_bp=1000)

    amplicon = finder.amplicon_df.iloc[0]

    # Should be 201 (inclusive), not 200
    assert amplicon.product_bp == 201
    assert amplicon.F_start == 1000
    assert amplicon.R_start == 1200

    # Verify the calculation: R_start - F_start + 1
    expected = amplicon.R_start - amplicon.F_start + 1
    assert amplicon.product_bp == expected


@pytest.mark.parametrize(
    "fwd_pos,rev_pos,expected_size",
    [
        (1000, 1050, 51),  # Small amplicon
        (2000, 2100, 101),  # Medium amplicon
        (5000, 5200, 201),  # Original test case
        (10000, 10001, 2),  # Tiny amplicon (edge case)
    ],
)
def test_amplicon_sizes_various(fwd_pos, rev_pos, expected_size):
    """Test various amplicon sizes to ensure inclusive coordinate calculation works generally."""
    cols = ["qseqid", "sseqid", "sstart", "sstrand"]
    data = [
        [f"FWD_{fwd_pos}", "chr1", fwd_pos, "plus"],
        [f"REV_{rev_pos}", "chr1", rev_pos, "minus"],
    ]
    bound_df = pd.DataFrame(data, columns=cols)

    finder = AmpliconFinder(bound_df)
    finder.find_amplicons(max_size_bp=2000)

    amplicon = finder.amplicon_df.iloc[0]
    assert amplicon.product_bp == expected_size
    assert amplicon.product_bp == rev_pos - fwd_pos + 1  # Inclusive calculation


def test_find_amplicons_uses_target_map(sample_bound_df):
    target_map = {"A_F": "GENE_A", "A_R": "GENE_A", "B_F": "GENE_B", "B_R": "GENE_B"}
    finder = AmpliconFinder(sample_bound_df, target_map=target_map)
    finder.find_amplicons(max_size_bp=1000)
    df = finder.amplicon_df

    # On-target chr7 amplicon: both primers belong to GENE_A
    chr7 = df[df["chrom"] == "chr7"].iloc[0]
    assert chr7["F_target"] == "GENE_A"
    assert chr7["R_target"] == "GENE_A"

    # Cross-target chr22 amplicon: A_F (4950) + B_R (5150)
    cross = df[(df["chrom"] == "chr22") & (df["F_primer"] == "A_F")]
    assert len(cross) == 1
    assert cross.iloc[0]["F_target"] == "GENE_A"
    assert cross.iloc[0]["R_target"] == "GENE_B"


def test_position_equality():
    p1 = Position("chr1", 100)
    p2 = Position("chr1", 100)
    p3 = Position("chr1", 101)
    assert p1 == p2
    assert p1 != p3
