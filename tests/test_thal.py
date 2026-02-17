# ================================================================================
# Tests for Thermodynamic alignment model
# ================================================================================

import pytest

from plexus.designer.thal import (
    InvalidConcentrationError,
    InvalidSequenceError,
    TmResult,
    divalent_to_monovalent,
    end_oligodg,
    oligodg,
    oligotm,
    seqtm,
    symmetry,
)


class TestThal:
    def test_seqtm_short_sequence(self):
        """Test seqtm function with a very short sequence."""
        seq = "ATGC"
        result = seqtm(seq)
        # A 4bp sequence has a very low Tm (primer3-py also gives -51.8)
        assert result.Tm == pytest.approx(-51.8, abs=1.0)

    def test_seqtm_primer_sequence_at_53C(self):
        """Test seqtm function with a realistic primer sequence at 53°C annealing.

        At 53°C (below Tm of 58.1°C), expect high binding (~90%).
        """
        seq = "CAGTGGCTCTATTGAATTTCTGTG"
        result = seqtm(
            seq,
            dna_conc=50.0,
            salt_conc=50.0,
            divalent_conc=1.5,
            dntp_conc=0.6,
            annealing_temp=53.0,
        )

        # Tm should match primer3-py calculation
        assert result.Tm == pytest.approx(58.1, abs=0.5)
        # Amount bound should be high (~90%) since we're below Tm
        assert result.bound == pytest.approx(90.2, abs=2.0)

    def test_seqtm_primer_sequence_at_60C(self):
        """Test seqtm function with a realistic primer sequence at 60°C annealing.

        At 60°C (above Tm of 58.1°C), expect lower binding (~31%).
        """
        seq = "CAGTGGCTCTATTGAATTTCTGTG"
        result = seqtm(
            seq,
            dna_conc=50.0,
            salt_conc=50.0,
            divalent_conc=1.5,
            dntp_conc=0.6,
            annealing_temp=60.0,
        )

        # Tm should match primer3-py calculation
        assert result.Tm == pytest.approx(58.1, abs=0.5)
        # Amount bound should be lower (~31%) since we're above Tm
        assert result.bound == pytest.approx(31.1, abs=2.0)

    def test_seqtm_high_gc_primer(self):
        """Test seqtm with a high GC content primer."""
        seq = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC, 20bp
        result = seqtm(seq)
        # High GC primers have higher Tm
        assert result.Tm > 70.0

    def test_oligotm_basic(self):
        """Test oligotm function directly."""
        seq = "ACGTACGTACGTACGTACGT"  # 20bp, 50% GC
        result = oligotm(seq)
        # Should return a reasonable Tm for a 20-mer
        assert 50.0 < result.Tm < 70.0
        assert 0.0 < result.bound < 100.0

    def test_symmetry(self):
        """Test symmetry detection for self-complementary sequences."""
        # Self-complementary sequences
        assert symmetry("ATAT") is True
        assert symmetry("GCGC") is True
        assert symmetry("ACGT") is True

        # Non-self-complementary sequences
        assert symmetry("AAAA") is False
        assert symmetry("ACGTA") is False  # Odd length
        assert symmetry("AATTGG") is False

    def test_seqtm_long_sequence_uses_gc_formula(self):
        """Test that seqtm uses the GC% formula for sequences > 60bp."""
        seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        assert len(seq) > 60
        result = seqtm(seq)
        assert isinstance(result, TmResult)
        assert isinstance(result.Tm, float)
        assert result.bound == 0.0
        assert result.ddG == 0.0

    def test_formamide_correction_without_dmso(self):
        """Test that formamide correction applies even when DMSO is 0."""
        seq = "CAGTGGCTCTATTGAATTTCTGTG"
        baseline = oligotm(seq, dmso_conc=0.0, formamide_conc=0.0)
        with_formamide = oligotm(seq, dmso_conc=0.0, formamide_conc=1.5)
        # Formamide lowers Tm (the correction term is negative for typical GC%)
        assert with_formamide.Tm < baseline.Tm

    def test_annealing_temperature_effect(self):
        """Test that lower annealing temp gives higher bound fraction."""
        seq = "CAGTGGCTCTATTGAATTTCTGTG"

        result_low_temp = seqtm(seq, annealing_temp=50.0)
        result_high_temp = seqtm(seq, annealing_temp=70.0)

        # Lower annealing temp should give higher bound fraction
        assert result_low_temp.bound > result_high_temp.bound
        # Tm should be the same regardless of annealing temp
        assert result_low_temp.Tm == pytest.approx(result_high_temp.Tm, abs=0.1)


class TestOligodg:
    def test_known_sequence(self):
        dg = oligodg("ACGT")
        assert isinstance(dg, float)
        assert dg == oligodg("ACGT")

    def test_too_short(self):
        with pytest.raises(InvalidSequenceError):
            oligodg("A")


class TestEndOligodg:
    def test_shorter_than_length(self):
        """Seq shorter than length returns oligodg(full seq)."""
        assert end_oligodg("ACGT", 10) == oligodg("ACGT")

    def test_longer_than_length(self):
        """Returns oligodg of last 3 bases."""
        assert end_oligodg("AACGT", 3) == oligodg("CGT")


class TestDivalentToMonovalent:
    def test_zero_divalent(self):
        assert divalent_to_monovalent(0, 0.6) == 0.0

    def test_negative_divalent(self):
        with pytest.raises(InvalidConcentrationError):
            divalent_to_monovalent(-1.0, 0.0)


class TestOligotmErrors:
    def test_too_short(self):
        with pytest.raises(InvalidSequenceError):
            oligotm("A")

    def test_invalid_chars(self):
        with pytest.raises(InvalidSequenceError):
            oligotm("ACGTXYZ")
