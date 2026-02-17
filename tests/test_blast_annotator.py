import pandas as pd
import pytest

from plexus.blast.annotator import BlastResultsAnnotator


@pytest.fixture
def sample_blast_df():
    """
    Create a sample BLAST dataframe matching -outfmt 6 with columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen
    """
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen".split()
    data = [
        # 1. Perfect 3' match (from_3prime=True, length=20, qend=20)
        ["P1_F", "chr7", 100.0, 20, 0, 0, 1, 20, 100, 119, 0.0, 40.0, "plus", 20],
        # 2. Internal match (from_3prime=False, qend=15, length=10)
        ["P1_F", "chr7", 100.0, 10, 0, 0, 6, 15, 200, 209, 0.0, 30.0, "plus", 20],
        # 3. Mismatching 3' match (pident=90)
        ["P1_R", "chr7", 90.0, 20, 2, 0, 1, 20, 300, 319, 0.0, 25.0, "minus", 20],
        # 4. High E-value match (evalue=10.0)
        ["P2_F", "chr22", 100.0, 12, 0, 0, 9, 12, 500, 511, 10.0, 20.0, "plus", 20],
        # 5. Low E-value 3' match (evalue=0.1, from_3prime=True)
        ["P2_R", "chr22", 100.0, 12, 0, 0, 9, 12, 600, 611, 0.1, 20.0, "minus", 20],
    ]
    return pd.DataFrame(data, columns=cols)


def test_build_annotation_dict(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=1.0)

    assert "from_3prime" in annotator.annotations
    assert "predicted_bound" in annotator.annotations


def test_add_annotations(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict(length_threshold=15, evalue_threshold=1.0)
    annotator.add_annotations()

    df = annotator.blast_df

    # Hit 0: Perfect 3' match
    assert bool(df.iloc[0]["from_3prime"]) is True
    assert bool(df.iloc[0]["length_pass_3prime"]) is True
    assert bool(df.iloc[0]["predicted_bound"]) is True

    # Hit 1: Internal match (qend=15, qlen=20)
    assert bool(df.iloc[1]["from_3prime"]) is False
    assert bool(df.iloc[1]["predicted_bound"]) is False

    # Hit 2: Mismatching (pident=90)
    assert bool(df.iloc[2]["length_pass_3prime"]) is False

    # Hit 4: Low E-value 3' match should pass via evalue_pass_3prime
    assert bool(df.iloc[4]["from_3prime"]) is True
    assert bool(df.iloc[4]["evalue_pass_3prime"]) is True
    assert bool(df.iloc[4]["predicted_bound"]) is True


def test_get_predicted_bound(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict()
    annotator.add_annotations()

    bound_df = annotator.get_predicted_bound()
    assert len(bound_df) == 4  # Hits 0, 2, 4, and 5
    assert all(bound_df["predicted_bound"])


def test_summarise_by_primer(sample_blast_df):
    annotator = BlastResultsAnnotator(sample_blast_df)
    annotator.build_annotation_dict()
    annotator.add_annotations()
    annotator.summarise_by_primer()

    summary = annotator.blast_primer_df
    assert "primer_name" in summary.columns
    assert "total_alignments" in summary.columns

    # P1_F has 2 hits, both are predicted_bound with default thresholds
    p1f = summary[summary["primer_name"] == "P1_F"].iloc[0]
    assert p1f["total_alignments"] == 2
    assert p1f["predicted_bound"] == 1  # Wait, Hit 0 is bound, Hit 1 is not. Correct.
