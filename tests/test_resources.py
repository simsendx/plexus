"""Tests for plexus.resources — SHA-256, checksums, registry, and global config."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from plexus.resources import (
    _compute_sha256,
    _parse_checksums_file,
    _register_genome_resources,
    get_operational_mode,
    set_operational_mode,
    verify_resource_checksums,
)

# ---------------------------------------------------------------------------
# SHA-256
# ---------------------------------------------------------------------------


def test_compute_sha256(tmp_path: Path):
    f = tmp_path / "test.txt"
    content = b"hello world\n"
    f.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert _compute_sha256(f) == expected


def test_compute_sha256_large_file(tmp_path: Path):
    """Verify chunked reading produces the same result as a single read."""
    f = tmp_path / "big.bin"
    data = b"x" * (2 * 1024 * 1024 + 7)  # just over 2 chunks
    f.write_bytes(data)
    assert _compute_sha256(f) == hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Checksums file parsing
# ---------------------------------------------------------------------------


def test_parse_checksums_file_standard(tmp_path: Path):
    """Standard sha256sum format with two-space separator."""
    f = tmp_path / "sha256sums.txt"
    f.write_text("abc123  genome.fa\ndef456  snps.vcf.gz\n")
    result = _parse_checksums_file(f)
    assert result == {"genome.fa": "abc123", "snps.vcf.gz": "def456"}


def test_parse_checksums_file_with_path(tmp_path: Path):
    """Filenames with directory paths are reduced to basename."""
    f = tmp_path / "sums.txt"
    f.write_text("abc123  /data/genomes/genome.fa\n")
    result = _parse_checksums_file(f)
    assert result == {"genome.fa": "abc123"}


def test_parse_checksums_file_skips_comments_and_blanks(tmp_path: Path):
    f = tmp_path / "sums.txt"
    f.write_text("# comment\n\nabc123  file.fa\n  \n")
    result = _parse_checksums_file(f)
    assert result == {"file.fa": "abc123"}


def test_parse_checksums_file_binary_mode_asterisk(tmp_path: Path):
    """Some tools use '*' prefix for binary mode."""
    f = tmp_path / "sums.txt"
    f.write_text("abc123 *genome.fa\n")
    result = _parse_checksums_file(f)
    assert result == {"genome.fa": "abc123"}


def test_parse_checksums_file_malformed(tmp_path: Path):
    f = tmp_path / "bad.txt"
    f.write_text("just_a_hash_no_filename\n")
    with pytest.raises(ValueError, match="Malformed"):
        _parse_checksums_file(f)


# ---------------------------------------------------------------------------
# Registry with checksums
# ---------------------------------------------------------------------------


def test_register_genome_resources_computes_checksums(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    vcf = tmp_path / "snps.vcf.gz"
    vcf.write_bytes(b"fake vcf data")

    registry_file = tmp_path / "registry.json"

    with patch("plexus.resources._registry_path", return_value=registry_file):
        _register_genome_resources("hg38", fasta, snp_vcf_path=vcf)

    reg = json.loads(registry_file.read_text())
    assert "hg38" in reg
    assert reg["hg38"]["fasta"] == str(fasta)
    assert reg["hg38"]["fasta_sha256"] == _compute_sha256(fasta)
    assert reg["hg38"]["snp_vcf"] == str(vcf)
    assert reg["hg38"]["snp_vcf_sha256"] == _compute_sha256(vcf)


def test_register_genome_resources_without_vcf(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    registry_file = tmp_path / "registry.json"

    with patch("plexus.resources._registry_path", return_value=registry_file):
        _register_genome_resources("hg38", fasta)

    reg = json.loads(registry_file.read_text())
    assert "fasta_sha256" in reg["hg38"]
    assert "snp_vcf" not in reg["hg38"]


def test_register_with_expected_checksums_match(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    expected = {fasta.name: _compute_sha256(fasta)}
    registry_file = tmp_path / "registry.json"

    with patch("plexus.resources._registry_path", return_value=registry_file):
        _register_genome_resources("hg38", fasta, expected_checksums=expected)

    reg = json.loads(registry_file.read_text())
    assert reg["hg38"]["fasta_sha256"] == expected[fasta.name]


def test_register_with_expected_checksums_mismatch(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    bad_checksums = {
        fasta.name: "0000000000000000000000000000000000000000000000000000000000000000"
    }
    registry_file = tmp_path / "registry.json"

    with patch("plexus.resources._registry_path", return_value=registry_file):
        with pytest.raises(ValueError, match="checksum mismatch"):
            _register_genome_resources("hg38", fasta, expected_checksums=bad_checksums)

    # Registry should NOT have been updated
    assert not registry_file.exists()


# ---------------------------------------------------------------------------
# Checksum verification
# ---------------------------------------------------------------------------


def test_verify_checksums_ok(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    sha = _compute_sha256(fasta)

    registry_file = tmp_path / "registry.json"
    registry_file.write_text(
        json.dumps({"hg38": {"fasta": str(fasta), "fasta_sha256": sha}})
    )

    with patch("plexus.resources._registry_path", return_value=registry_file):
        result = verify_resource_checksums("hg38")
    assert result["fasta"] is True
    assert result["snp_vcf"] is None  # not registered


def test_verify_checksums_mismatch(tmp_path: Path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")

    registry_file = tmp_path / "registry.json"
    registry_file.write_text(
        json.dumps({"hg38": {"fasta": str(fasta), "fasta_sha256": "wrong_hash"}})
    )

    with patch("plexus.resources._registry_path", return_value=registry_file):
        result = verify_resource_checksums("hg38")
    assert result["fasta"] is False


def test_verify_checksums_missing_file(tmp_path: Path):
    registry_file = tmp_path / "registry.json"
    registry_file.write_text(
        json.dumps(
            {"hg38": {"fasta": str(tmp_path / "gone.fa"), "fasta_sha256": "abc"}}
        )
    )

    with patch("plexus.resources._registry_path", return_value=registry_file):
        result = verify_resource_checksums("hg38")
    assert result["fasta"] is False


def test_verify_checksums_no_entry(tmp_path: Path):
    registry_file = tmp_path / "registry.json"
    registry_file.write_text("{}")

    with patch("plexus.resources._registry_path", return_value=registry_file):
        result = verify_resource_checksums("hg38")
    assert result["fasta"] is None
    assert result["snp_vcf"] is None


# ---------------------------------------------------------------------------
# Global config / operational mode
# ---------------------------------------------------------------------------


def test_get_operational_mode_default(tmp_path: Path):
    config_path = tmp_path / "config.json"
    with patch("plexus.resources._global_config_path", return_value=config_path):
        assert get_operational_mode() == "research"


def test_set_and_get_operational_mode(tmp_path: Path):
    config_path = tmp_path / "config.json"
    with patch("plexus.resources._global_config_path", return_value=config_path):
        set_operational_mode("compliance")
        assert get_operational_mode() == "compliance"
        set_operational_mode("research")
        assert get_operational_mode() == "research"


def test_set_invalid_mode(tmp_path: Path):
    config_path = tmp_path / "config.json"
    with patch("plexus.resources._global_config_path", return_value=config_path):
        with pytest.raises(ValueError, match="Invalid mode"):
            set_operational_mode("invalid")
