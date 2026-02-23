"""Panel QC report generation (REPT-01)."""

from __future__ import annotations

import re
import statistics
from itertools import combinations
from typing import TYPE_CHECKING

from plexus.aligner import PrimerDimerPredictor

if TYPE_CHECKING:
    from plexus.designer.multiplexpanel import Junction


def generate_panel_qc(
    junctions: list[Junction],
    *,
    gc_high_threshold: float = 70.0,
    gc_low_threshold: float = 30.0,
    homopolymer_min_run: int = 4,
    dimer_threshold: float = 0.0,
) -> dict:
    """Generate panel QC metrics for the selected primer pairs."""
    # Build working lists
    selected = []  # list of (junction_name, PrimerPair)
    for j in junctions:
        for pair in j.primer_pairs:
            if pair.selected:
                selected.append((j.name, pair))
                break

    all_primers = []  # (junction_name, direction, Primer)
    for jname, pair in selected:
        all_primers.append((jname, "forward", pair.forward))
        all_primers.append((jname, "reverse", pair.reverse))

    # Tm distribution
    tms = [p.tm for _, _, p in all_primers]
    tm_distribution = {
        "mean": round(statistics.mean(tms), 2) if tms else None,
        "std": round(statistics.stdev(tms), 2) if len(tms) >= 2 else None,
        "min": round(min(tms), 2) if tms else None,
        "max": round(max(tms), 2) if tms else None,
        "per_primer": [
            {"junction": jname, "direction": d, "name": p.name, "tm": round(p.tm, 2)}
            for jname, d, p in all_primers
        ],
    }

    # Sequence flags
    _hp_re = re.compile(
        r"([ACGT])\1{" + str(homopolymer_min_run - 1) + r",}", re.IGNORECASE
    )
    flagged_primers = []
    high_gc_count = low_gc_count = homopolymer_count = 0
    for jname, direction, p in all_primers:
        flags = []
        if p.gc > gc_high_threshold:
            flags.append("high_gc")
            high_gc_count += 1
        if p.gc < gc_low_threshold:
            flags.append("low_gc")
            low_gc_count += 1
        if _hp_re.search(p.seq):
            flags.append("homopolymer")
            homopolymer_count += 1
        if flags:
            flagged_primers.append(
                {
                    "junction": jname,
                    "direction": direction,
                    "name": p.name,
                    "gc": round(p.gc, 1),
                    "sequence": p.seq,
                    "flags": flags,
                }
            )
    sequence_flags = {
        "gc_high_threshold": gc_high_threshold,
        "gc_low_threshold": gc_low_threshold,
        "homopolymer_min_run": homopolymer_min_run,
        "high_gc_count": high_gc_count,
        "low_gc_count": low_gc_count,
        "homopolymer_count": homopolymer_count,
        "flagged_primers": flagged_primers,
    }

    # Cross-reactivity matrix
    predictor = PrimerDimerPredictor()
    matrix: dict[str, dict] = {}
    for (jname_a, pair_a), (jname_b, pair_b) in combinations(selected, 2):
        scores = []
        for seq_a, name_a, seq_b, name_b in [
            (
                pair_a.forward.seq,
                pair_a.forward.name,
                pair_b.forward.seq,
                pair_b.forward.name,
            ),
            (
                pair_a.forward.seq,
                pair_a.forward.name,
                pair_b.reverse.seq,
                pair_b.reverse.name,
            ),
            (
                pair_a.reverse.seq,
                pair_a.reverse.name,
                pair_b.forward.seq,
                pair_b.forward.name,
            ),
            (
                pair_a.reverse.seq,
                pair_a.reverse.name,
                pair_b.reverse.seq,
                pair_b.reverse.name,
            ),
        ]:
            predictor.set_primers(seq_a, seq_b, name_a, name_b)
            predictor.align()
            scores.append(predictor.score or 0.0)

        cell = {
            "min_dimer_score": round(min(scores), 4),
            "interaction_count": sum(1 for s in scores if s < dimer_threshold),
        }
        matrix.setdefault(jname_a, {})[jname_b] = cell
        matrix.setdefault(jname_b, {})[jname_a] = cell

    cross_reactivity_matrix = {
        "dimer_threshold": dimer_threshold,
        "matrix": matrix,
    }

    return {
        "tm_distribution": tm_distribution,
        "sequence_flags": sequence_flags,
        "cross_reactivity_matrix": cross_reactivity_matrix,
    }
