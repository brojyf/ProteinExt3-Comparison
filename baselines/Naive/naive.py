from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path


def extract_entry_id(header: str) -> str:
    title = header[1:] if header.startswith(">") else header
    parts = title.split("|")
    if len(parts) >= 2:
        return parts[1]
    return title.split()[0]


def read_query_ids(fasta_path: Path) -> list[str]:
    protein_ids = []
    with fasta_path.open("r", encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                protein_ids.append(extract_entry_id(line.strip()))
    return protein_ids


def read_truth_by_pid(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    truth_by_pid: dict[str, set[str]] = defaultdict(set)
    go_to_aspect: dict[str, str] = {}

    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line:
                continue

            parts = [part.strip() for part in line.split("\t")]
            if len(parts) < 2:
                continue

            if parts[0] in {"AUTHOR", "MODEL", "KEYWORDS", "END"}:
                continue
            if parts[0] == "EntryID" and len(parts) >= 2 and parts[1] == "term":
                continue

            protein_id = parts[0]
            go_term = parts[1]
            if not protein_id or not go_term.startswith("GO:"):
                continue

            truth_by_pid[protein_id].add(go_term)
            if len(parts) >= 3 and parts[2].lower() in {"c", "f", "p"}:
                go_to_aspect[go_term] = parts[2].lower()

    if not truth_by_pid:
        raise ValueError(f"无法从 {path} 解析任何 GO 注释")

    return truth_by_pid, go_to_aspect


def build_io_paths(base_prefix: str, input_prefix: str, output_path: str) -> dict[str, Path]:
    return {
        "base_truth": Path(f"{base_prefix}.tsv"),
        "input_fasta": Path(f"{input_prefix}.fasta"),
        "input_truth": Path(f"{input_prefix}.tsv"),
        "output_tsv": Path(output_path),
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="fasta/cafa/cafa", help="训练集前缀，会自动使用 .tsv")
    parser.add_argument("--in", dest="input_prefix", default="fasta/test", help="测试集前缀，会自动使用 .fasta")
    parser.add_argument("--out", default="baselines/Naive/predictions.tsv")
    parser.add_argument("--min-prob", type=float, default=0.01)
    return parser.parse_args()


def main():
    args = parse_args()
    io_paths = build_io_paths(args.base, args.input_prefix, args.out)

    base_truth = io_paths["base_truth"]
    input_fasta = io_paths["input_fasta"]
    output_tsv = io_paths["output_tsv"]

    for path in [base_truth, input_fasta]:
        if not path.exists():
            raise FileNotFoundError(f"找不到文件: {path}")

    truth_by_pid, go_to_aspect = read_truth_by_pid(base_truth)
    query_ids = read_query_ids(input_fasta)

    go_counts = Counter()
    for go_terms in truth_by_pid.values():
        for go_term in go_terms:
            go_counts[go_term] += 1

    total_proteins = len(truth_by_pid)
    go_probs = {
        go_term: count / total_proteins
        for go_term, count in go_counts.items()
        if count / total_proteins >= args.min_prob
    }

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8") as file:
        for protein_id in query_ids:
            for go_term, prob in sorted(go_probs.items(), key=lambda item: item[1], reverse=True):
                aspect = go_to_aspect.get(go_term, "")
                file.write(f"{protein_id}\t{go_term}\t{prob:.6f}\t{aspect}\n")

    print(f"base_truth: {base_truth}")
    print(f"input_fasta: {input_fasta}")
    print(f"input_truth: {io_paths['input_truth']}")
    print(f"output_tsv: {output_tsv}")
    print(f"base_proteins: {total_proteins}")
    print(f"query_proteins: {len(query_ids)}")
    print(f"predicted_go_terms: {len(go_probs)}")
    print(f"prediction_rows: {len(query_ids) * len(go_probs)}")


if __name__ == "__main__":
    main()
