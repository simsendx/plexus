# ================================================================================
# Tests for SNP overlap checking module
# ================================================================================

from unittest.mock import MagicMock, patch

import pytest

from plexus.config import DesignerConfig, SnpCheckParameters
from plexus.designer.primer import Primer, PrimerPair
from plexus.snpcheck.checker import (
    _count_snps_in_region,
    _get_allele_frequency,
    _primer_genomic_coords,
    run_snp_check,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_primer(name="test_fwd", direction="forward", start=10, length=20):
    return Primer(name=name, seq="A" * length, direction=direction, start=start, length=length)


def _make_pair(fwd_start=10, fwd_len=20, rev_start=50, rev_len=20):
    fwd = _make_primer("fwd", "forward", fwd_start, fwd_len)
    rev = _make_primer("rev", "reverse", rev_start, rev_len)
    return PrimerPair(
        forward=fwd,
        reverse=rev,
        insert_size=rev_start - (fwd_start + fwd_len),
        amplicon_sequence="A" * 80,
        amplicon_length=80,
        pair_penalty=1.0,
        pair_id="test_pair_0",
    )


def _make_junction(design_start=1000, chrom="chr7"):
    junction = MagicMock()
    junction.name = "TEST_JUNCTION"
    junction.chrom = chrom
    junction.start = design_start + 200
    junction.end = design_start + 200
    junction.design_start = design_start
    junction.primer_pairs = [_make_pair()]
    return junction


def _make_panel(junctions=None):
    panel = MagicMock()
    panel.junctions = junctions or [_make_junction()]
    return panel


# ---------------------------------------------------------------------------
# Coordinate mapping
# ---------------------------------------------------------------------------


class TestPrimerGenomicCoords:
    def test_basic_mapping(self):
        """Primer at relative start=10 with design_start=1000 → genomic 1010."""
        primer = _make_primer(start=10, length=20)
        junction = _make_junction(design_start=1000)
        chrom, gstart, gend = _primer_genomic_coords(primer, junction)
        assert chrom == "chr7"
        assert gstart == 1010
        assert gend == 1030

    def test_zero_offset(self):
        """Primer at the start of the design region."""
        primer = _make_primer(start=0, length=25)
        junction = _make_junction(design_start=5000)
        _, gstart, gend = _primer_genomic_coords(primer, junction)
        assert gstart == 5000
        assert gend == 5025


# ---------------------------------------------------------------------------
# Allele frequency extraction
# ---------------------------------------------------------------------------


class TestGetAlleleFrequency:
    def test_af_field_single(self):
        record = MagicMock()
        record.info = {"AF": 0.05}
        assert _get_allele_frequency(record) == 0.05

    def test_af_field_list(self):
        record = MagicMock()
        record.info = {"AF": [0.01, 0.03]}
        assert _get_allele_frequency(record) == 0.03

    def test_caf_field(self):
        record = MagicMock()
        record.info = {"CAF": ["0.95", "0.03", "0.02"]}
        assert _get_allele_frequency(record) == 0.03

    def test_caf_with_missing(self):
        record = MagicMock()
        record.info = {"CAF": ["0.95", ".", "0.05"]}
        assert _get_allele_frequency(record) == 0.05

    def test_no_frequency_info(self):
        record = MagicMock()
        record.info = {}
        assert _get_allele_frequency(record) is None


# ---------------------------------------------------------------------------
# SNP counting
# ---------------------------------------------------------------------------


class TestCountSnpsInRegion:
    def test_no_variants(self):
        vcf = MagicMock()
        vcf.fetch.return_value = []
        count, positions = _count_snps_in_region(vcf, "chr7", 1000, 1020, 0.01)
        assert count == 0
        assert positions == []

    def test_variants_above_threshold(self):
        record1 = MagicMock()
        record1.info = {"AF": 0.05}
        record1.pos = 1005
        record2 = MagicMock()
        record2.info = {"AF": 0.001}  # Below threshold
        record2.pos = 1010

        vcf = MagicMock()
        vcf.fetch.return_value = [record1, record2]

        count, positions = _count_snps_in_region(vcf, "chr7", 1000, 1020, 0.01)
        assert count == 1
        assert positions == [1005]

    def test_contig_not_in_vcf(self):
        """ValueError from pysam when contig doesn't exist should be handled."""
        vcf = MagicMock()
        vcf.fetch.side_effect = ValueError("contig not found")
        count, positions = _count_snps_in_region(vcf, "chrUn", 100, 200, 0.01)
        assert count == 0
        assert positions == []


# ---------------------------------------------------------------------------
# run_snp_check integration
# ---------------------------------------------------------------------------


class TestRunSnpCheck:
    def test_applies_penalty(self):
        """SNP overlaps should increase snp_count and pair_penalty."""
        panel = _make_panel()
        pair = panel.junctions[0].primer_pairs[0]
        original_penalty = pair.pair_penalty

        # Mock pysam.VariantFile
        mock_record = MagicMock()
        mock_record.info = {"AF": 0.05}
        mock_record.pos = 1015

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_vcf.fetch.return_value = [mock_record]
            mock_variant_file.return_value = mock_vcf

            run_snp_check(
                panel=panel,
                vcf_path="/fake/dbsnp.vcf.gz",
                af_threshold=0.01,
                snp_penalty_weight=5.0,
            )

        # Each primer gets 1 SNP → pair gets 2 SNPs
        assert pair.snp_count == 2
        assert pair.snp_penalty == 10.0
        assert pair.pair_penalty == original_penalty + 10.0

    def test_no_snps_no_penalty(self):
        """No SNPs means penalty unchanged."""
        panel = _make_panel()
        pair = panel.junctions[0].primer_pairs[0]
        original_penalty = pair.pair_penalty

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_vcf.fetch.return_value = []
            mock_variant_file.return_value = mock_vcf

            run_snp_check(
                panel=panel,
                vcf_path="/fake/dbsnp.vcf.gz",
                af_threshold=0.01,
                snp_penalty_weight=5.0,
            )

        assert pair.snp_count == 0
        assert pair.snp_penalty == 0.0
        assert pair.pair_penalty == original_penalty

    def test_skips_junction_without_design_start(self):
        """Junctions without design_start are skipped gracefully."""
        junction = _make_junction()
        junction.design_start = None
        panel = _make_panel([junction])

        with patch("pysam.VariantFile") as mock_variant_file:
            mock_vcf = MagicMock()
            mock_variant_file.return_value = mock_vcf

            # Should not raise
            run_snp_check(panel=panel, vcf_path="/fake/dbsnp.vcf.gz")

        # fetch should never be called
        mock_vcf.fetch.assert_not_called()


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


class TestSnpCheckConfig:
    def test_default_values(self):
        params = SnpCheckParameters()
        assert params.af_threshold == 0.01
        assert params.snp_penalty_weight == 5.0

    def test_custom_values(self):
        params = SnpCheckParameters(af_threshold=0.05, snp_penalty_weight=10.0)
        assert params.af_threshold == 0.05
        assert params.snp_penalty_weight == 10.0

    def test_included_in_designer_config(self):
        config = DesignerConfig()
        assert hasattr(config, "snp_check_parameters")
        assert config.snp_check_parameters.af_threshold == 0.01

    def test_af_threshold_validation(self):
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            SnpCheckParameters(af_threshold=-0.1)
        with pytest.raises(ValidationError):
            SnpCheckParameters(af_threshold=1.5)


# ---------------------------------------------------------------------------
# PrimerPair SNP fields
# ---------------------------------------------------------------------------


class TestPrimerPairSnpFields:
    def test_default_zero(self):
        pair = _make_pair()
        assert pair.snp_count == 0
        assert pair.snp_penalty == 0.0

    def test_primer_default_zero(self):
        primer = _make_primer()
        assert primer.snp_count == 0
