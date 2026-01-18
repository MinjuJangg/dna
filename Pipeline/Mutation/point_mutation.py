from __future__ import annotations
import random
from typing import List, Dict, Optional


def point_mutation(
    genes: List[str],
    n: int = 1,
    *,
    p_ins: float = 0.01,
    p_del: float = 0.06,
    p_sub: float = 0.93,
    insert_base_probs: Optional[Dict[str, float]] = None,
    sub_matrix: Optional[Dict[str, Dict[str, float]]] = None,
    seed: Optional[int] = None,
) -> List[str]:
    """
    genes: 유전자 문자열 리스트 (예: ["ATGC", "GGTA", ...])
    n: 변이 이벤트를 총 몇 번 수행할지 (매 이벤트마다 유전자를 다시 랜덤 선택)
    반환: 변이가 반영된 유전자 문자열 리스트 (원본 genes는 변경하지 않음)
    """

    # 0) 원본 유지: 복사본을 만들어 이걸 수정
    mutated_genes = list(genes)
    if not mutated_genes or n <= 0:
        return mutated_genes

    rng = random.Random(seed)

    # 1) INS 삽입 염기 확률 기본값
    if insert_base_probs is None:
        insert_base_probs = {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}

    # 2) SUB 치환 행렬 기본값 (질문 표)
    if sub_matrix is None:
        sub_matrix = {
            "A": {"T": 0.303, "G": 0.455, "C": 0.242},
            "T": {"A": 0.268, "G": 0.293, "C": 0.439},
            "G": {"A": 0.479, "T": 0.333, "C": 0.188},
            "C": {"A": 0.212, "T": 0.577, "G": 0.212},
        }

    # 3) INS/DEL/SUB 확률 정규화 (합이 1이 아니어도 상대비로 동작)
    total_p = p_ins + p_del + p_sub
    if total_p <= 0:
        raise ValueError("p_ins + p_del + p_sub 는 0보다 커야 합니다.")
    p_ins /= total_p
    p_del /= total_p
    p_sub /= total_p

    # 4) n번 변이 수행 (매번 유전자 선택도 다시 함)
    for _ in range(n):
        # 4-1) 변이할 유전자 하나를 랜덤 선택
        gene_idx = rng.randrange(len(mutated_genes))
        gene = mutated_genes[gene_idx]

        # 4-2) 변이 타입 선택(확률)
        r = rng.random()
        if r < p_ins:
            mut_type = "INS"
        elif r < p_ins + p_del:
            mut_type = "DEL"
        else:
            mut_type = "SUB"

        # 4-3) 선택된 유전자를 list로 바꿔서 수정(직관적으로 다루기 쉬움)
        seq = list(gene)
        L = len(seq)

        # 4-4) 변이 타입에 따라 처리
        if mut_type == "INS":
            # INS: "사이사이" 위치는 모두 동일 확률
            pos = rng.randrange(L + 1)

            # 삽입될 염기는 insert_base_probs 확률대로 선택
            bases = list(insert_base_probs.keys())
            weights = list(insert_base_probs.values())
            base_to_insert = rng.choices(bases, weights=weights, k=1)[0]

            seq.insert(pos, base_to_insert)

        elif mut_type == "DEL":
            # DEL: 염기 1개를 동일 확률로 선택 후 제거
            if L == 0:
                continue
            pos = rng.randrange(L)
            seq.pop(pos)

        else:  # SUB
            # SUB: 염기 1개를 동일 확률로 선택 후 표 확률대로 치환
            if L == 0:
                continue
            pos = rng.randrange(L)
            from_base = seq[pos]

            # 표에 없는 염기면 치환 불가 → 이번 변이는 스킵
            if from_base not in sub_matrix:
                continue

            to_dict = sub_matrix[from_base]
            to_bases = list(to_dict.keys())
            to_weights = list(to_dict.values())
            to_base = rng.choices(to_bases, weights=to_weights, k=1)[0]

            seq[pos] = to_base

        # 4-5) 수정된 seq를 문자열로 되돌려 리스트에 반영
        mutated_genes[gene_idx] = "".join(seq)

    return mutated_genes
