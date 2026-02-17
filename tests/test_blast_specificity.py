from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from plexus.blast.specificity import run_specificity_check
from plexus.designer.multiplexpanel import MultiplexPanel


@pytest.fixture
def mock_panel():
    panel = MagicMock(spec=MultiplexPanel)
    panel.unique_primer_map = {"ACTG": "P1_F", "CAGT": "P1_R"}

    # Mock Junction
    junction = MagicMock()
    junction.chrom = "chr7"
    junction.name = "J1"

    # Mock PrimerPair
    pair = MagicMock()
    pair.forward.seq = "ACTG"
    pair.reverse.seq = "CAGT"
    pair.amplicon_length = 200
    pair.pair_id = "J1_P1"

    junction.primer_pairs = [pair]
    panel.junctions = [junction]

    return panel


def test_run_specificity_check_integration(mock_panel, tmp_path):
    # Setup mocks for internal classes
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("plexus.blast.specificity.BlastResultsAnnotator") as MockAnnotator,
        patch("plexus.blast.specificity.AmpliconFinder") as MockFinder,
        patch("os.makedirs"),
    ):
        # 1. Mock Runner behavior
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame({"dummy": [1]})

        # 2. Mock Annotator behavior
        annotator_instance = MockAnnotator.return_value
        annotator_instance.get_predicted_bound.return_value = pd.DataFrame(
            {"dummy_bound": [1]}
        )

        # 3. Mock Finder behavior
        finder_instance = MockFinder.return_value
        # One off-target amplicon
        finder_instance.amplicon_df = pd.DataFrame(
            [
                {
                    "chrom": "chr7",
                    "F_primer": "P1_F",
                    "R_primer": "P1_R",
                    "product_bp": 500,  # Not 200, so off-target
                    "F_start": 1000,
                    "R_start": 1500,
                }
            ]
        )

        # Run
        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        # Verify
        pair = mock_panel.junctions[0].primer_pairs[0]
        assert pair.specificity_checked is True
        assert len(pair.off_target_products) == 1
        assert pair.off_target_products[0]["product_bp"] == 500


def test_run_specificity_check_no_hits(mock_panel, tmp_path):
    with (
        patch("plexus.blast.specificity.BlastRunner") as MockRunner,
        patch("os.makedirs"),
    ):
        runner_instance = MockRunner.return_value
        runner_instance.get_dataframe.return_value = pd.DataFrame()  # Empty

        # Ensure it's False or unset initially
        pair = mock_panel.junctions[0].primer_pairs[0]
        pair.specificity_checked = False

        # Should return early without error
        run_specificity_check(mock_panel, str(tmp_path), "fake_genome.fa")

        assert pair.specificity_checked is False
