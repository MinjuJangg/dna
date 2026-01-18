import random

def recombination(
    target_genes: list, 
    donor_genes: list, 
    mutation_rate: float = 0.03,
    recomb_ratio: float = 0.1, 
    min_copy_len: int = 10, 
    max_copy_len: int = 100
) -> list:
    """
    [Recombination Module]
    문서의 '7.1 Recombination Process' 구현
    기증자(Donor) 유전자의 일부를 복사하여 타겟(Target)에 덮어씁니다.
    """
    
    current_genes = list(target_genes)
    
    if not current_genes or not donor_genes:
        print("Warning: 유전자 리스트가 비어있습니다.")
        return current_genes

    for idx, _ in enumerate(current_genes):
        target_seq = current_genes[idx]
        
        # 유전자 길이에 따른 이벤트 횟수 자동 계산 (반올림 적용)
        # 공식: 유전자 길이 * 변이율(0.03) * 재조합비율(0.1)
        expected_events = len(target_seq) * mutation_rate * recomb_ratio
        num_events = int(round(expected_events))
        
        for _ in range(num_events):
            target_seq = current_genes[idx] # 최신 서열 갱신

            # 1. 기증자 유전자 선택 [cite: 235]
            donor_seq = random.choice(donor_genes)
            if len(donor_seq) <= min_copy_len:
                continue

            # 2. 복사할 위치와 길이 선택 [cite: 236]
            actual_max_len = min(len(donor_seq), max_copy_len)
            if actual_max_len < min_copy_len:
                continue
                
            segment_len = random.randint(min_copy_len, actual_max_len)
            src_start = random.randint(0, len(donor_seq) - segment_len)
            copied_segment = donor_seq[src_start : src_start + segment_len]
            
            # 3. 붙여넣기 수행 (Overwrite) [cite: 238]
            if len(target_seq) > 0:
                tgt_start = random.randint(0, len(target_seq) - 1)
                
                new_gene_seq = (
                    target_seq[:tgt_start] + 
                    copied_segment + 
                    target_seq[tgt_start + segment_len:]
                )
                current_genes[idx] = new_gene_seq
                
            else:
                current_genes[idx] = copied_segment

    return current_genes