# ================================================================================
# SNP overlap checking using the Ensembl REST API
#
# Queries the Ensembl overlap/region endpoint for variants overlapping primer
# binding sites. No local data download required — only the primer regions
# (~20-30bp each) are queried. Rate-limited to 14 req/sec.
#
# Author: Stefan Filges (stefan@simsendiagnostics.com)
# Copyright (c) 2025-2026 Simsen Diagnostics AB
# ================================================================================

from __future__ import annotations

import time

import requests
from loguru import logger

from plexus.designer.multiplexpanel import MultiplexPanel
from plexus.snpcheck.checker import _primer_genomic_coords

ENSEMBL_REST_URL = "https://rest.ensembl.org"
REQUEST_INTERVAL = 0.07  # ~14 req/sec (Ensembl allows 15/sec)


def _query_variants_in_region(chrom: str, start: int, end: int) -> list[dict]:
    """Query Ensembl REST API for variants overlapping a genomic region.

    Parameters
    ----------
    chrom : str
        Chromosome (e.g. "chr7" or "7"). Ensembl expects no "chr" prefix.
    start : int
        1-based inclusive start.
    end : int
        1-based inclusive end.

    Returns
    -------
    list[dict]
        List of variant records from Ensembl.
    """
    # Ensembl uses chromosome names without "chr" prefix
    chrom_clean = chrom.removeprefix("chr")

    url = f"{ENSEMBL_REST_URL}/overlap/region/human/{chrom_clean}:{start}-{end}"
    params = {"feature": "variation", "content-type": "application/json"}

    response = requests.get(url, params=params, timeout=30)
    if response.status_code == 429:
        # Rate limited — wait and retry once
        retry_after = float(response.headers.get("Retry-After", 1))
        logger.debug(f"Rate limited, waiting {retry_after}s...")
        time.sleep(retry_after)
        response = requests.get(url, params=params, timeout=30)

    if not response.ok:
        logger.warning(
            f"Ensembl API error for {chrom_clean}:{start}-{end}: "
            f"{response.status_code} {response.text[:200]}"
        )
        return []

    return response.json()


def _get_max_maf(variant: dict) -> float | None:
    """Extract the minor allele frequency from an Ensembl variant record.

    Returns the MAF as reported by Ensembl, or None if not available.
    """
    # Ensembl overlap endpoint returns "minor_allele_freq" directly
    maf = variant.get("minor_allele_freq")
    if maf is not None:
        return float(maf)
    return None


def run_snp_check_api(
    panel: MultiplexPanel,
    af_threshold: float = 0.01,
    snp_penalty_weight: float = 5.0,
    padding: int = 200,
) -> None:
    """Check all primer pairs for SNP overlaps using the Ensembl REST API.

    Same interface and behavior as run_snp_check (VCF version) but queries
    the Ensembl API instead of a local file.

    Parameters
    ----------
    panel : MultiplexPanel
    af_threshold : float
        Minimum allele frequency (default 0.01 = 1%).
    snp_penalty_weight : float
        Penalty per SNP (default 5.0).
    padding : int
        Unused, kept for API compatibility.
    """
    logger.info("Running SNP check via Ensembl REST API")
    logger.info(f"AF threshold: {af_threshold}, penalty weight: {snp_penalty_weight}")

    total_snps = 0
    primers_with_snps = 0
    queries = 0

    for junction in panel.junctions:
        if not junction.primer_pairs:
            continue
        if not hasattr(junction, "design_start") or junction.design_start is None:
            logger.warning(
                f"No design_start for junction {junction.name}, skipping SNP check"
            )
            continue

        for pair in junction.primer_pairs:
            pair_snp_count = 0

            for primer in (pair.forward, pair.reverse):
                chrom, gstart, gend = _primer_genomic_coords(primer, junction)

                # Query Ensembl (1-based inclusive on both ends)
                variants = _query_variants_in_region(chrom, gstart, gend - 1)
                queries += 1
                time.sleep(REQUEST_INTERVAL)

                # Filter by allele frequency
                count = 0
                positions = []
                for var in variants:
                    maf = _get_max_maf(var)
                    if maf is not None and maf >= af_threshold:
                        count += 1
                        positions.append(var.get("start", 0))

                primer.snp_count = count
                pair_snp_count += count

                if count > 0:
                    primers_with_snps += 1
                    total_snps += count
                    logger.debug(
                        f"Primer {primer.name}: {count} SNPs at positions {positions}"
                    )

            pair.snp_count = pair_snp_count
            pair.snp_penalty = pair_snp_count * snp_penalty_weight

            if pair.snp_penalty > 0:
                if pair.pair_penalty is not None:
                    pair.pair_penalty += pair.snp_penalty
                else:
                    pair.pair_penalty = pair.snp_penalty

    logger.info(
        f"SNP check complete: {total_snps} SNPs found across {primers_with_snps} primers "
        f"({queries} API queries)"
    )
