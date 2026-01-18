import random
from typing import Set, Dict, List


def duplication(genelist: Set[str], seq_dict: Dict[str, str],
    duplication_num: int, rng: random.Random | None = None) -> List[str]:
    
    rng = rng or random
    sampled_genes = set(rng.sample(tuple(genelist), duplication_num))
    genelist.difference_update(sampled_genes)

    sequences: List[str] = []
    for acc in sampled_genes:
        if acc not in seq_dict:
            raise KeyError(f"Sequence not found in FAA for {acc}")
        sequences.append(seq_dict[acc])

    return sequences
