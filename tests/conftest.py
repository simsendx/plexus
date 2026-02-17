"""Shared pytest fixtures for integration tests."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def fixture_fasta() -> Path:
    """Path to the small 2-contig FASTA fixture."""
    p = DATA_DIR / "fixtures.fa"
    assert p.exists(), f"Missing fixture: {p}"
    return p


@pytest.fixture(scope="session")
def fixture_vcf() -> Path:
    """Path to the remapped gnomAD VCF fixture."""
    p = DATA_DIR / "snps.vcf.gz"
    assert p.exists(), f"Missing fixture: {p}"
    return p


@pytest.fixture(scope="session")
def fixture_csv() -> Path:
    """Path to the 3-junction CSV fixture."""
    p = DATA_DIR / "junctions.csv"
    assert p.exists(), f"Missing fixture: {p}"
    return p
