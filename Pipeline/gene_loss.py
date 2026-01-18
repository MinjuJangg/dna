import random
from typing import List, Optional


def gene_loss(duplicated_genes: List[str], loss_ratio: float = 0.3, seed: Optional[int] = None) -> List[str]:
    n_total = len(duplicated_genes)
    if n_total == 0:
        return []
    
    rng = random.Random(seed)
    n_remove = int(round(n_total * loss_ratio))
    n_remove = min(n_remove, n_total)
    remove_indices = set(rng.sample(range(n_total), n_remove))

    return [gene for i, gene in enumerate(duplicated_genes) if i not in remove_indices]
