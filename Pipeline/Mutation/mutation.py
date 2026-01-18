import math
from typing import List

from .recombination import recombination
from .point_mutation import point_mutation


def run_mutation(duplicated_genes: List[str], remaining_genes: List[str], point_mutation_num: int) -> List[str]:
    if point_mutation_num == 0:
        return duplicated_genes

    rounds = math.ceil(point_mutation_num / 10)

    for r in range(rounds):
        duplicated_genes = recombination(target_genes=duplicated_genes, donor_genes=remaining_genes)
        remaining = point_mutation_num - (r * 10)
        k = 10 if remaining >= 10 else remaining
        for _ in range(k):
            duplicated_genes = point_mutation(genes=duplicated_genes)

    return duplicated_genes
