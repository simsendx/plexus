# Changelog

All notable changes to plexus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 13-02-2026

### Added

- SNP overlap checking for primer binding sites with two backends:
  - Ensembl REST API (default, no setup required)
  - Local tabix-indexed VCF/BCF via pysam
- New CLI options: `--snp-vcf`, `--skip-snpcheck`, `--snp-af-threshold`.
- `SnpCheckParameters` configuration model with `af_threshold` and
  `snp_penalty_weight` settings.
- SNP columns (`SNP_Count`, `SNP_Penalty`, `Forward_SNP_Count`,
  `Reverse_SNP_Count`) in output CSVs.
- `snp_count` field on `Primer`; `snp_count` and `snp_penalty` fields on
  `PrimerPair`.
- 19 new tests for the snpcheck module.
- README badges for CI, Python version, and license.

### Changed

- Updated Dockerfile to Python 3.13 with optimized layer caching.
- Capped Python version to `<3.14` due to pysam compatibility.
- Added `requests` to dependencies.

### Removed

- Legacy `snpcheck` module (multiply/bedtools dependencies).

## [0.1.1] - 13-02-2026

### Added

- Changes project name to `plexus`
- Multi-panel support: input CSV can now include an optional `Panel` column to
  define multiple panels in a single file. Each panel is designed and optimized
  independently with output saved to `output_dir/<panel_id>/`.
- New `orchestrator.py` module with `run_multi_panel()` entry point that detects
  panels, dispatches per-panel pipeline runs, and aggregates results.
- `MultiPanelResult` dataclass for aggregated results across panels.
- `--parallel` CLI flag to run panels concurrently via `ProcessPoolExecutor`.
- `--max-workers` CLI flag to control parallel worker count.
- `multi_panel_summary.json` output file for multi-panel runs.
- Tests for orchestrator (20 new tests) and `JunctionInput` Panel field.

### Changed

- CLI now routes through `run_multi_panel()` which delegates to `run_pipeline()`
  per panel. Single-panel CSVs (no Panel column) behave identically to before.
- Marked `primer3py` design method as experimental with warnings; it does not yet
  produce `PrimerPair` objects for downstream multiplex optimization.
- Marked `primer3` (local binary) design method as not implemented with a clear
  error message pointing users to the `simsen` method.
- `PipelineResult` is now constructed after panel creation succeeds, eliminating
  the `panel=None` / `type: ignore` workaround. The result always contains a
  valid `MultiplexPanel`.

### Removed

- Removed `run_designer()` convenience function and associated hardcoded paths
  from `designer/design.py`. Use `run_pipeline()` or the CLI instead.

## [0.1.0] - 13-02-2026

Initial release.
