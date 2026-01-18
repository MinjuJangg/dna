from pathlib import Path
from typing import Set, Dict, List


def load_faa_sequences(faa_path: str | Path) -> Dict[str, str]:
    faa_path = Path(faa_path)
    seq_dict: Dict[str, str] = {}
    current_acc = None
    current_seq = []
    with faa_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_acc is not None:
                    seq_dict[current_acc] = "".join(current_seq)
                current_acc = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_acc is not None:
            seq_dict[current_acc] = "".join(current_seq)
    return seq_dict


def compute_average_gene_length(faa_path: str | Path) -> float:
    faa_path = Path(faa_path)
    lengths = []
    current_len = 0

    with faa_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    if not lengths:
        raise ValueError("No gene sequences found in FAA file.")
    avg_len = sum(lengths) / len(lengths)
    return avg_len


def get_remaining_genes(genelist: Set[str], seq_dict: Dict[str, str]) -> List[str]:
    remaining_genes: List[str] = []

    for acc in genelist:
        if acc not in seq_dict:
            raise KeyError(f"Sequence not found for accession: {acc}")
        remaining_genes.append(seq_dict[acc])

    return remaining_genes


def load_orphan_genes_as_sequences(orphan_path: str | Path, seq_dict: Dict[str, str]) -> List[str]:
    orphan_path = Path(orphan_path)
    orphan_genes: List[str] = []

    with orphan_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            acc = line.strip()
            if not acc:
                continue

            if acc not in seq_dict:
                raise KeyError(f"Orphan gene accession not found in FAA: {acc}")

            orphan_genes.append(seq_dict[acc])

    return orphan_genes
