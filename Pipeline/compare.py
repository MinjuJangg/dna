from __future__ import annotations
from pathlib import Path
from typing import List, Tuple
from datetime import datetime


def sequence_identity(seq1: str, seq2: str) -> float:
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 0.0
    min_len = min(len(seq1), len(seq2))
    matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
    return matches / max_len


def compare(orphan_genes: List[str], final_genes: List[str], threshold: float,
    results_dir: str | Path, filename: str, *, epoch_idx: int) -> float:
    if not orphan_genes or not final_genes:
        return 0.0

    matched_pairs: List[Tuple[int, int, float]] = []

    max_ident = 0.0
    max_i = -1
    max_j = -1

    for i, og in enumerate(orphan_genes):
        for j, fg in enumerate(final_genes):
            ident = sequence_identity(og, fg)

            if ident > max_ident:
                max_ident = ident
                max_i, max_j = i, j

            if ident >= threshold:
                matched_pairs.append((i, j, ident))

    if not matched_pairs:
        return max_ident

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / filename

    with out_path.open("a", encoding="utf-8") as w:
        w.write("\n" + "=" * 80 + "\n")
        w.write(f"EPOCH: {epoch_idx}\n")
        w.write(f"RUN RESULT | threshold={threshold}\n")
        w.write(f"Time: {datetime.now().isoformat(timespec='seconds')}\n")
        w.write(f"Orphan genes: {len(orphan_genes)}\n")
        w.write(f"Final genes: {len(final_genes)}\n")
        w.write(f"Max identity (any pair): {max_ident*100:.2f}% (orphan_idx={max_i}, final_idx={max_j})\n")
        w.write(f"Matched pairs (>=threshold): {len(matched_pairs)}\n\n")
        w.write("[MAX PAIR]\n")
        w.write(f"identity={max_ident*100:.2f}% orphan_idx={max_i} final_idx={max_j}\n")
        w.write("[ORPHAN]\n")
        w.write(orphan_genes[max_i] + "\n")
        w.write("[FINAL]\n")
        w.write(final_genes[max_j] + "\n")
        w.write("-" * 80 + "\n\n")
        for i, j, ident in matched_pairs:
            w.write(f"[PAIR] orphan_idx={i}, final_idx={j}, identity={ident*100:.2f}%\n")
            w.write("[ORPHAN]\n")
            w.write(orphan_genes[i] + "\n")
            w.write("[FINAL]\n")
            w.write(final_genes[j] + "\n")
            w.write("-" * 80 + "\n")

    return max_ident
