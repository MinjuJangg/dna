import argparse
from pathlib import Path
from typing import Dict, Any, List
import yaml

from util import load_faa_sequences, compute_average_gene_length, get_remaining_genes, load_orphan_genes_as_sequences
from setting_genelist import setting_genelist
from duplication import duplication
from Mutation.mutation import run_mutation
from gene_loss import gene_loss
from compare import compare


def load_yaml(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"YAML config not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("YAML must be a key-value mapping")

    return cfg


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_yaml(args.config)

    if "faa_path" not in cfg:
        raise KeyError("Missing required key: faa_path")
    if "orphan_path" not in cfg:
        raise KeyError("Missing required key: orphan_path")
    if "year_num" not in cfg:
        raise KeyError("Missing required key: year_num")
    if "epoch_num" not in cfg:
        raise KeyError("Missing required key: epoch_num")

    faa_path = cfg["faa_path"]
    orphan_path = cfg["orphan_path"]
    year_num = cfg['year_num']
    epoch_num = cfg['epoch_num']

    yaml_name = Path(args.config).stem 
    result_filename = f"{yaml_name}_gene_overlap.txt"

    seq_dict = load_faa_sequences(faa_path)
    genelist: List[str] = setting_genelist(faa_path=faa_path, orphan_path=orphan_path)
    orphan_genelist = load_orphan_genes_as_sequences(orphan_path=orphan_path, seq_dict=seq_dict)
    print(f"[OK] genelist size = {len(genelist)}")

    num_genes = len(genelist)
    avg_gene_length = round(compute_average_gene_length(faa_path))
    duplication_num = round(num_genes * 10**(-9) * year_num)

    duplicated_genes = duplication(genelist=genelist, seq_dict=seq_dict, duplication_num=duplication_num)
    remaining_genes = get_remaining_genes(genelist=genelist, seq_dict=seq_dict)
    print(f"[OK] genelist size = {len(genelist)}")
    print(f"[OK] duplicated_genes size = {len(duplicated_genes)}")

    point_mutation_num = round(0.03 * avg_gene_length)
    duplicated_genes = run_mutation(duplicated_genes, remaining_genes, point_mutation_num)
    ratio = compare(orphan_genes=orphan_genelist, final_genes=duplicated_genes,
        threshold=0.8, results_dir="Results",filename=result_filename,epoch_idx=0)
    if ratio >= 0.999999:
        print("[TERMINATE] 100% overlap achieved at epoch 0.")
        return
    left_genes = gene_loss(duplicated_genes)
    
    for r in range(1, epoch_num+1):
        if len(genelist) <= 0:
            print(f"[STOP] genelist is empty at epoch {r}.")
            break
        if duplication_num <= 0:
            print(f"[STOP] duplication_num={duplication_num} (<=0) at epoch {r}.")
            break
        if duplication_num > len(genelist):
            print(f"[STOP] Not enough genes to sample: need {duplication_num}, have {len(genelist)} at epoch {r}.")
            break
        duplicated_genes = duplication(genelist=genelist, seq_dict=seq_dict, duplication_num=duplication_num)
        duplicated_genes += left_genes
        print(f"[OK] genelist size = {len(genelist)}")
        duplicated_genes = run_mutation(duplicated_genes, remaining_genes, point_mutation_num)
        ratio = compare(orphan_genes=orphan_genelist, final_genes=duplicated_genes,
            threshold=0.8, results_dir="Results",filename=result_filename, epoch_idx=r)
        if ratio >= 0.999999:
            print(f"[TERMINATE] 100% overlap achieved at epoch {r}.")
            break
        left_genes = gene_loss(duplicated_genes)



if __name__ == "__main__":
    main()
