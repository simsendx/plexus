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
