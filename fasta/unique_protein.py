from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


ASPECT_BY_NAMESPACE = {
    "biological_process": "p",
    "molecular_function": "f",
    "cellular_component": "c",
}


def read_fasta(path):
    records = []
    header = None
    sequence_lines = []

    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_lines).upper()
                    records.append((header, sequence))
                header = line
                sequence_lines = []
                continue

            sequence_lines.append(line)

    if header is not None:
        sequence = "".join(sequence_lines).upper()
        records.append((header, sequence))

    return records


def write_fasta(records, path, width=60):
    with open(path, "w") as file:
        for header, sequence in records:
            file.write(header + "\n")
            for start in range(0, len(sequence), width):
                file.write(sequence[start : start + width] + "\n")


def extract_entry_id(header):
    title = header[1:] if header.startswith(">") else header
    parts = title.split("|")
    if len(parts) >= 2:
        return parts[1]
    return title.split()[0]


def _is_truth_header(parts):
    left = parts[0].strip().lower()
    right = parts[1].strip().lower()
    return left in {"entryid", "pid", "protein_id"} and right in {
        "term",
        "goterm",
        "go_term",
    }


def read_truth_tsv(path):
    rows = []
    truth_by_pid = defaultdict(set)

    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue
            if _is_truth_header(parts):
                continue

            protein_id = parts[0].strip()
            go_term = parts[1].strip()
            if not protein_id or not go_term:
                continue

            rows.append((protein_id, go_term))
            truth_by_pid[protein_id].add(go_term)

    return rows, {protein_id: sorted(go_terms) for protein_id, go_terms in truth_by_pid.items()}


def load_go_aspect_map(obo_path):
    go_aspect = {}
    current_id = None
    current_namespace = None
    in_term = False

    def flush():
        if current_id and current_namespace in ASPECT_BY_NAMESPACE:
            go_aspect[current_id] = ASPECT_BY_NAMESPACE[current_namespace]

    with open(obo_path, "r") as file:
        for raw_line in file:
            line = raw_line.strip()

            if line == "[Term]":
                flush()
                current_id = None
                current_namespace = None
                in_term = True
                continue

            if line.startswith("[") and line != "[Term]":
                flush()
                current_id = None
                current_namespace = None
                in_term = False
                continue

            if not in_term:
                continue

            if line.startswith("id: GO:"):
                current_id = line.split("id:", 1)[1].strip()
            elif line.startswith("namespace:"):
                current_namespace = line.split("namespace:", 1)[1].strip()

    flush()
    return go_aspect


def write_truth_tsv(rows, path):
    with open(path, "w") as file:
        for protein_id, go_term, aspect in rows:
            file.write(f"{protein_id}\t{go_term}\t{aspect}\n")


def _unique_in_order(values):
    seen = set()
    unique_values = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        unique_values.append(value)
    return unique_values


def _build_truth_rows(output_pids, truth_source_by_pid, go_aspect):
    truth_rows = []
    missing_aspect_terms = set()
    missing_truth_pids = []

    for protein_id in _unique_in_order(output_pids):
        go_terms = truth_source_by_pid.get(protein_id)
        if not go_terms:
            missing_truth_pids.append(protein_id)
            continue

        for go_term in go_terms:
            aspect = go_aspect.get(go_term)
            if aspect is None:
                missing_aspect_terms.add(go_term)
                continue
            truth_rows.append((protein_id, go_term, aspect))

    if missing_aspect_terms:
        preview = ", ".join(sorted(missing_aspect_terms)[:10])
        raise KeyError(f"以下 GO term 没有在 OBO 中找到 aspect: {preview}")

    return truth_rows, missing_truth_pids


def merge(infasta1, infasta2, infastatruth1, infastatruth2, output_fasta, output_truth, obo_path):
    _, truth_by_pid1 = read_truth_tsv(infastatruth1)
    _, truth_by_pid2 = read_truth_tsv(infastatruth2)

    go_aspect = load_go_aspect_map(obo_path)
    records1 = read_fasta(infasta1)
    records2 = read_fasta(infasta2)

    sequence_to_candidates = defaultdict(list)
    for source_order, (records, truth_by_pid) in enumerate(
        [(records1, truth_by_pid1), (records2, truth_by_pid2)]
    ):
        for index, (header, sequence) in enumerate(records):
            protein_id = extract_entry_id(header)
            sequence_to_candidates[sequence].append(
                {
                    "header": header,
                    "sequence": sequence,
                    "protein_id": protein_id,
                    "truth_by_pid": truth_by_pid,
                    "has_truth": protein_id in truth_by_pid,
                    "source_order": source_order,
                    "index": index,
                }
            )

    output_records = []
    output_pids = []
    truth_source_by_pid = {}
    skipped_without_truth = 0

    for candidates in sequence_to_candidates.values():
        selected = min(
            candidates,
            key=lambda item: (
                0 if item["has_truth"] else 1,
                item["source_order"],
                item["index"],
            ),
        )

        if not selected["has_truth"]:
            skipped_without_truth += 1
            continue

        output_records.append((selected["header"], selected["sequence"]))
        output_pids.append(selected["protein_id"])
        truth_source_by_pid[selected["protein_id"]] = selected["truth_by_pid"][
            selected["protein_id"]
        ]

    truth_rows, missing_truth_pids = _build_truth_rows(output_pids, truth_source_by_pid, go_aspect)

    write_fasta(output_records, output_fasta)
    write_truth_tsv(truth_rows, output_truth)

    return {
        "mode": "merge",
        "input_records_1": len(records1),
        "input_records_2": len(records2),
        "unique_sequences_after_merge": len(sequence_to_candidates),
        "output_fasta_records": len(output_records),
        "output_truth_rows": len(truth_rows),
        "output_truth_proteins": len(_unique_in_order(output_pids)),
        "skipped_without_truth": skipped_without_truth,
        "missing_truth_pids": missing_truth_pids,
        "output_fasta": output_fasta,
        "output_truth": output_truth,
        "obo_path": obo_path,
    }


def set_difference(
    infasta1, infasta2, infastatruth1, infastatruth2, output_fasta, output_truth, obo_path
):
    del infastatruth2

    _, truth_by_pid1 = read_truth_tsv(infastatruth1)
    go_aspect = load_go_aspect_map(obo_path)
    records1 = read_fasta(infasta1)
    records2 = read_fasta(infasta2)

    sequences_b = {sequence for _, sequence in records2}
    output_records = []
    output_pids = []
    truth_source_by_pid = {}
    kept_sequences = set()
    skipped_overlap = 0
    skipped_without_truth = 0
    skipped_duplicate_in_a = 0

    for header, sequence in records1:
        protein_id = extract_entry_id(header)

        if sequence in sequences_b:
            skipped_overlap += 1
            continue

        if sequence in kept_sequences:
            skipped_duplicate_in_a += 1
            continue

        if protein_id not in truth_by_pid1:
            skipped_without_truth += 1
            continue

        kept_sequences.add(sequence)
        output_records.append((header, sequence))
        output_pids.append(protein_id)
        truth_source_by_pid[protein_id] = truth_by_pid1[protein_id]

    truth_rows, missing_truth_pids = _build_truth_rows(output_pids, truth_source_by_pid, go_aspect)

    write_fasta(output_records, output_fasta)
    write_truth_tsv(truth_rows, output_truth)

    return {
        "mode": "set_difference",
        "input_records_a": len(records1),
        "input_records_b": len(records2),
        "output_fasta_records": len(output_records),
        "output_truth_rows": len(truth_rows),
        "output_truth_proteins": len(_unique_in_order(output_pids)),
        "skipped_overlap": skipped_overlap,
        "skipped_without_truth": skipped_without_truth,
        "skipped_duplicate_in_a": skipped_duplicate_in_a,
        "missing_truth_pids": missing_truth_pids,
        "output_fasta": output_fasta,
        "output_truth": output_truth,
        "obo_path": obo_path,
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices=["set-difference", "merge"])
    parser.add_argument("--in1", required=True, help="输入 1 的公共前缀，脚本会自动补 .fasta 和 .tsv")
    parser.add_argument("--in2", required=True, help="输入 2 的公共前缀，脚本会自动补 .fasta 和 .tsv")
    parser.add_argument("--out", required=True, help="输出公共前缀，脚本会自动补 .fasta 和 .tsv")
    parser.add_argument(
        "--obo",
        default=str(Path.cwd() / "go-basic.obo"),
        help="GO OBO 文件路径。默认使用当前目录下的 go-basic.obo",
    )
    return parser.parse_args()


def build_io_paths(input_prefix_1, input_prefix_2, output_prefix):
    return {
        "infasta1": f"{input_prefix_1}.fasta",
        "infasta2": f"{input_prefix_2}.fasta",
        "infastatruth1": f"{input_prefix_1}.tsv",
        "infastatruth2": f"{input_prefix_2}.tsv",
        "output_fasta": f"{output_prefix}.fasta",
        "output_truth": f"{output_prefix}.tsv",
    }


def main():
    args = parse_args()
    current_dir = str(Path.cwd())
    io_paths = build_io_paths(args.in1, args.in2, args.out)

    if args.method == "merge":
        summary = merge(
            infasta1=io_paths["infasta1"],
            infasta2=io_paths["infasta2"],
            infastatruth1=io_paths["infastatruth1"],
            infastatruth2=io_paths["infastatruth2"],
            output_fasta=io_paths["output_fasta"],
            output_truth=io_paths["output_truth"],
            obo_path=args.obo,
        )
    else:
        summary = set_difference(
            infasta1=io_paths["infasta1"],
            infasta2=io_paths["infasta2"],
            infastatruth1=io_paths["infastatruth1"],
            infastatruth2=io_paths["infastatruth2"],
            output_fasta=io_paths["output_fasta"],
            output_truth=io_paths["output_truth"],
            obo_path=args.obo,
        )

    print(f"current_dir: {current_dir}")
    for key, value in io_paths.items():
        print(f"{key}: {value}")
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
