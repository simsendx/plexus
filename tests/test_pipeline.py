"""Tests for AUDT-03: provenance.json lifecycle — status, completed_at, steps_completed.

These tests focus solely on the provenance update behaviour introduced in AUDT-03.
They use run_blast=False and skip_snpcheck=True to avoid requiring BLAST+/bcftools
so they run in any environment.  Full pipeline correctness is covered by test_integration.py.
"""

import json
from unittest.mock import patch

import pytest

from plexus.pipeline import run_pipeline


def test_provenance_completed_status(fixture_fasta, fixture_csv, tmp_path):
    """provenance.json should have status=completed after a successful run."""
    out = tmp_path / "output"
    run_pipeline(
        input_file=fixture_csv,
        fasta_file=fixture_fasta,
        output_dir=out,
        run_blast=False,
        skip_snpcheck=True,
    )

    prov = json.loads((out / "provenance.json").read_text())
    assert prov["status"] == "completed"
    assert prov["completed_at"] is not None
    assert isinstance(prov["steps_completed"], list)
    assert isinstance(prov["errors"], list)
    # Original fields are preserved
    assert "plexus_version" in prov
    assert "run_timestamp" in prov


def test_provenance_failed_status_on_exception(fixture_fasta, fixture_csv, tmp_path):
    """provenance.json must have status=failed and a non-null completed_at when the
    pipeline raises an unhandled exception mid-run."""
    out = tmp_path / "output"

    with patch(
        "plexus.pipeline.design_primers",
        side_effect=RuntimeError("simulated design crash"),
    ):
        with pytest.raises(RuntimeError, match="simulated design crash"):
            run_pipeline(
                input_file=fixture_csv,
                fasta_file=fixture_fasta,
                output_dir=out,
                run_blast=False,
                skip_snpcheck=True,
            )

    assert (out / "provenance.json").exists()
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["status"] == "failed"
    assert prov["completed_at"] is not None
    assert isinstance(prov["steps_completed"], list)
    assert isinstance(prov["errors"], list)


def test_provenance_started_before_pipeline_work(fixture_fasta, fixture_csv, tmp_path):
    """provenance.json should have status=started at the point design_primers is called,
    and status=completed after the pipeline finishes."""
    out = tmp_path / "output"
    statuses_during_run: list[str] = []

    import plexus.designer.design as _design_mod

    original_design = _design_mod.design_multiplex_primers

    def capture_and_design(*args, **kwargs):
        try:
            prov = json.loads((out / "provenance.json").read_text())
            statuses_during_run.append(prov.get("status", "missing"))
        except Exception:
            statuses_during_run.append("unreadable")
        return original_design(*args, **kwargs)

    with patch(
        "plexus.pipeline.design_primers",
        wraps=lambda panel, method="simsen", **kw: capture_and_design(panel, **kw),
    ):
        run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=out,
            run_blast=False,
            skip_snpcheck=True,
        )

    # Provenance was in "started" state when design_primers was called
    assert statuses_during_run, "design_primers was never called"
    assert statuses_during_run[0] == "started"

    # After the run, provenance should be "completed"
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["status"] == "completed"


def test_provenance_steps_completed_recorded(fixture_fasta, fixture_csv, tmp_path):
    """steps_completed in provenance.json should list the pipeline stages that ran."""
    out = tmp_path / "output"
    run_pipeline(
        input_file=fixture_csv,
        fasta_file=fixture_fasta,
        output_dir=out,
        run_blast=False,
        skip_snpcheck=True,
    )

    prov = json.loads((out / "provenance.json").read_text())
    assert "primers_designed" in prov["steps_completed"]
    assert "multiplex_optimized" in prov["steps_completed"]


# ── REPR-01 tests ─────────────────────────────────────────────────────────────


def test_chrom_naming_check_skipped_recorded(fixture_fasta, fixture_csv, tmp_path):
    """chrom_naming_check=skipped when skip_snpcheck=True."""
    out = tmp_path / "output"
    run_pipeline(
        input_file=fixture_csv,
        fasta_file=fixture_fasta,
        output_dir=out,
        run_blast=False,
        skip_snpcheck=True,
    )
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["chrom_naming_check"] == "skipped"


def test_chrom_naming_mismatch_raises_in_compliance(
    fixture_fasta, fixture_csv, fixture_vcf, tmp_path
):
    """Compliance mode + mismatch → ValueError before any pipeline work."""
    out = tmp_path / "output"
    with (
        patch("plexus.resources.get_operational_mode", return_value="compliance"),
        patch("plexus.utils.env.validate_environment", return_value={}),
        patch("plexus.pipeline.get_missing_tools", return_value=[]),
        patch("plexus.utils.utils.get_fasta_contigs", return_value={"chr1", "chr2"}),
        patch("plexus.utils.utils.get_vcf_contigs", return_value={"1", "2"}),
    ):
        with pytest.raises(ValueError, match="[Cc]hromosome naming mismatch"):
            run_pipeline(
                input_file=fixture_csv,
                fasta_file=fixture_fasta,
                output_dir=out,
                snp_vcf=str(fixture_vcf),
                run_blast=False,
                skip_snpcheck=False,
            )
    # provenance should NOT exist — pipeline aborted before the write
    assert not (out / "provenance.json").exists()


def test_chrom_naming_mismatch_warns_in_research(
    fixture_fasta, fixture_csv, fixture_vcf, tmp_path
):
    """Research mode + mismatch → warning, pipeline continues, provenance records 'mismatch'."""
    out = tmp_path / "output"
    with (
        patch("plexus.resources.get_operational_mode", return_value="research"),
        patch("plexus.pipeline.get_missing_tools", return_value=[]),
        patch("plexus.utils.utils.get_fasta_contigs", return_value={"chr1", "chr2"}),
        patch("plexus.utils.utils.get_vcf_contigs", return_value={"1", "2"}),
        patch("plexus.snpcheck.snp_data.get_snp_vcf", return_value=fixture_vcf),
        patch("plexus.snpcheck.checker.run_snp_check"),
    ):
        run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=out,
            snp_vcf=str(fixture_vcf),
            run_blast=False,
            skip_snpcheck=False,
        )
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["chrom_naming_check"] == "mismatch"
    assert prov["status"] == "completed"


def test_chrom_naming_check_pass_recorded(
    fixture_fasta, fixture_csv, fixture_vcf, tmp_path
):
    """When contigs match, chrom_naming_check=pass is recorded in provenance."""
    out = tmp_path / "output"
    matching = {"chr1", "chr2", "chr3"}
    with (
        patch("plexus.pipeline.get_missing_tools", return_value=[]),
        patch("plexus.utils.utils.get_fasta_contigs", return_value=matching),
        patch("plexus.utils.utils.get_vcf_contigs", return_value=matching),
        patch("plexus.snpcheck.snp_data.get_snp_vcf", return_value=fixture_vcf),
        patch("plexus.snpcheck.checker.run_snp_check"),
    ):
        run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=out,
            snp_vcf=str(fixture_vcf),
            run_blast=False,
            skip_snpcheck=False,
        )
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["chrom_naming_check"] == "pass"


# ── Registry-resolution tests ──────────────────────────────────────────────


def test_chrom_naming_check_runs_via_registry(
    fixture_fasta, fixture_csv, fixture_vcf, tmp_path
):
    """chrom_naming_check should not be 'skipped' when registry resolves a VCF."""
    out = tmp_path / "output"
    matching = {"chr1", "chr2", "chr3"}
    with (
        patch("plexus.pipeline.get_missing_tools", return_value=[]),
        patch("plexus.resources.get_registered_snp_vcf", return_value=fixture_vcf),
        patch("plexus.utils.utils.get_fasta_contigs", return_value=matching),
        patch("plexus.utils.utils.get_vcf_contigs", return_value=matching),
        patch("plexus.snpcheck.snp_data.get_snp_vcf", return_value=fixture_vcf),
        patch("plexus.snpcheck.checker.run_snp_check"),
    ):
        run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=out,
            snp_vcf=None,  # no explicit VCF — registry should be used
            run_blast=False,
            skip_snpcheck=False,
        )
    prov = json.loads((out / "provenance.json").read_text())
    assert (
        prov["chrom_naming_check"] != "skipped"
    ), "chrom_naming_check should run when a registry VCF is available"


def test_snp_vcf_provenance_backfilled_from_registry(
    fixture_fasta, fixture_csv, fixture_vcf, tmp_path
):
    """snp_vcf_path and snp_vcf_sha256 should be non-null when the registry resolves a VCF."""
    out = tmp_path / "output"
    fake_sha = "abc123fake"

    def fake_registry():
        return {
            "hg38": {
                "snp_vcf": str(fixture_vcf),
                "snp_vcf_sha256": fake_sha,
                "fasta": "",
                "fasta_sha256": None,
            }
        }

    with (
        patch("plexus.pipeline.get_missing_tools", return_value=[]),
        patch("plexus.resources.get_registered_snp_vcf", return_value=fixture_vcf),
        patch("plexus.pipeline._load_registry", fake_registry, create=True),
        patch("plexus.resources._load_registry", fake_registry),
        patch("plexus.utils.utils.get_fasta_contigs", return_value={"chr1"}),
        patch("plexus.utils.utils.get_vcf_contigs", return_value={"chr1"}),
        patch("plexus.snpcheck.snp_data.get_snp_vcf", return_value=fixture_vcf),
        patch("plexus.snpcheck.checker.run_snp_check"),
    ):
        run_pipeline(
            input_file=fixture_csv,
            fasta_file=fixture_fasta,
            output_dir=out,
            snp_vcf=None,
            run_blast=False,
            skip_snpcheck=False,
        )
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["snp_vcf_path"] == str(
        fixture_vcf
    ), "snp_vcf_path should be the registry-resolved path"
    assert (
        prov["snp_vcf_sha256"] == fake_sha
    ), "snp_vcf_sha256 should be populated from registry"
