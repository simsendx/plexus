# ================================================================================
# Environment and dependency verification utilities
# ================================================================================

from __future__ import annotations

import importlib.metadata
import shutil
import subprocess


def check_executable(name: str) -> bool:
    """Check if an executable exists on the PATH."""
    return shutil.which(name) is not None


def get_missing_tools(need_blast: bool = False, need_snp: bool = False) -> list[str]:
    """
    Identify missing system dependencies based on run requirements.

    Args:
        need_blast: Whether BLAST+ tools are required.
        need_snp: Whether bcftools is required for SNP checking.

    Returns:
        List of missing executable names.
    """
    missing = []
    if need_blast:
        # We need both for a full run (db creation, searching, and formatting)
        for tool in ("blastn", "makeblastdb", "blast_formatter"):
            if not check_executable(tool):
                missing.append(tool)
    if need_snp:
        if not check_executable("bcftools"):
            missing.append("bcftools")
    return missing


# ---------------------------------------------------------------------------
# Tool version capture
# ---------------------------------------------------------------------------


def get_tool_version(name: str) -> str | None:
    """Return version string for a system tool, or ``None`` if absent.

    Tries ``<name> --version`` first, then ``<name> -version``.
    Returns the first non-empty line of stdout/stderr.
    """
    if not check_executable(name):
        return None
    for flag in ("--version", "-version"):
        try:
            result = subprocess.run(
                [name, flag],
                capture_output=True,
                text=True,
                timeout=5,
            )
            output = (result.stdout or result.stderr or "").strip()
            if output:
                return output.splitlines()[0]
        except (subprocess.TimeoutExpired, OSError):
            continue
    return "unknown"


def get_tool_versions(tools: list[str]) -> dict[str, str | None]:
    """Return ``{tool_name: version_string}`` for each tool (``None`` if missing)."""
    return {t: get_tool_version(t) for t in tools}


def get_plexus_version() -> str:
    """Return the installed plexus version."""
    try:
        return importlib.metadata.version("plexus")
    except importlib.metadata.PackageNotFoundError:
        from plexus.version import __version__

        return __version__


def get_primer3_version() -> str | None:
    """Return the installed primer3-py version, or ``None``."""
    try:
        return importlib.metadata.version("primer3-py")
    except importlib.metadata.PackageNotFoundError:
        return None
