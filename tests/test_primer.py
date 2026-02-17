# ================================================================================
# Tests for Primer and PrimerPair classes
# ================================================================================

import pytest

from plexus.designer.primer import Primer, PrimerPair


class TestPrimerAddTail:
    def test_add_tail_five_prime(self):
        p = Primer(name="p", seq="ACGT", direction="forward", start=0, length=4)
        p.add_tail("TTT", "five_prime")
        assert p.seq == "TTTACGT"

    def test_add_tail_three_prime(self):
        p = Primer(name="p", seq="ACGT", direction="forward", start=0, length=4)
        p.add_tail("TTT", "three_prime")
        assert p.seq == "ACGTTTT"


class TestPrimerPairPenalty:
    def test_pair_penalty_no_product_opt(self):
        """product_opt_size=0 -> only primer penalties + Tm diff."""
        penalty = PrimerPair.calculate_primer_pair_penalty_th(
            primer_left_penalty=1.0,
            primer_right_penalty=2.0,
            primer_left_tm=60.0,
            primer_right_tm=62.0,
            product_size=100,
            product_opt_size=0,
            wt_diff_tm=1.0,
        )
        assert penalty == pytest.approx(5.0)  # 1*(1+2) + 1*|60-62|

    def test_pair_penalty_with_product_size_gt(self):
        """Product larger than optimal -> wt_product_size_gt penalty."""
        penalty = PrimerPair.calculate_primer_pair_penalty_th(
            primer_left_penalty=0.0,
            primer_right_penalty=0.0,
            primer_left_tm=60.0,
            primer_right_tm=60.0,
            product_size=150,
            product_opt_size=100,
            wt_product_size_gt=0.5,
        )
        assert penalty == pytest.approx(25.0)  # 0.5*(150-100)
