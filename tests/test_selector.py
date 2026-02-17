# ================================================================================
# Tests for selector modules (Multiplex, GreedySearch, RandomSearch, BruteForce)
# ================================================================================

from unittest.mock import MagicMock

import pandas as pd
import pytest

from plexus.selector.multiplex import Multiplex
from plexus.selector.selectors import BruteForce, GreedySearch, RandomSearch


class TestMultiplex:
    def test_defaults(self):
        m = Multiplex()
        assert m.primer_pairs == []
        assert m.cost == 0.0


@pytest.fixture
def selector_inputs():
    df = pd.DataFrame(
        {
            "target_id": ["T1", "T1", "T2", "T2"],
            "pair_name": ["P1a", "P1b", "P2a", "P2b"],
        }
    )
    cost_fn = MagicMock()
    cost_fn.calc_cost = MagicMock(return_value=1.0)
    return df, cost_fn


class TestGreedySearch:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = GreedySearch(df, cost_fn).run(N=5)
        assert len(results) == 5
        assert all(isinstance(m, Multiplex) for m in results)
        assert all(len(m.primer_pairs) == 2 for m in results)


class TestRandomSearch:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = RandomSearch(df, cost_fn).run(N=5)
        assert len(results) == 5
        assert all(isinstance(m, Multiplex) for m in results)


class TestBruteForce:
    def test_returns_multiplexes(self, selector_inputs):
        df, cost_fn = selector_inputs
        results = BruteForce(df, cost_fn).run(store_maximum=10)
        assert len(results) == 4  # 2 pairs x 2 pairs = 4 combos
        assert all(isinstance(m, Multiplex) for m in results)
