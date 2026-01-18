from pathlib import Path
from typing import Set


def setting_genelist(faa_path: str | Path, orphan_path: str | Path) -> Set[str]:
    faa_path = Path(faa_path)
    orphan_path = Path(orphan_path)

    remove_set: Set[str] = set()
    with orphan_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            acc = line.strip()
            if acc:
                remove_set.add(acc)

    genelist: Set[str] = set()

    with faa_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            acc = line[1:].split()[0]
            if not acc.startswith("NP_"):
                continue
            if acc in remove_set:
                continue
            genelist.add(acc)

    return genelist
