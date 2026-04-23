from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


ASPECTS = {"p", "f", "c"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    sequence_lines = []

    with path.open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))
                header = line
                sequence_lines = []
                continue

            sequence_lines.append(line)

    if header is not None:
        records.append((header, "".join(sequence_lines)))

    return records


def write_fasta(records: list[tuple[str, str]], path: Path, width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as file:
        for header, sequence in records:
            file.write(f"{header}\n")
            for start in range(0, len(sequence), width):
                file.write(f"{sequence[start:start + width]}\n")


def fasta_record_id(header: str) -> str:
    title = header[1:] if header.startswith(">") else header
    parts = title.split("|")
    if len(parts) >= 2:
        return parts[1]
    return title.split()[0]


def is_tsv_header(parts: list[str]) -> bool:
    if len(parts) < 2:
        return False
    left = parts[0].strip().lower()
    right = parts[1].strip().lower()
    return left in {"entryid", "pid", "protein_id"} and right in {
        "term",
        "goterm",
        "go_term",
    }


def read_label_space(path: Path) -> set[str]:
    labels = set()
    with path.open("r", encoding="utf-8") as file:
        for line in file:
            raw_line = line.rstrip("\n")
            if not raw_line:
                continue

            parts = raw_line.split("\t")
            if len(parts) < 2 or is_tsv_header(parts):
                continue

            go_term = parts[1].strip()
            if go_term.startswith("GO:"):
                labels.add(go_term)

    return labels


def load_obo_parents(path: Path) -> dict[str, set[str]]:
    parents = defaultdict(set)
    current_term = None

    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
            elif line.startswith("id: "):
                if current_term is None:
                    current_term = line[4:].strip()
            elif line.startswith("is_a: "):
                if current_term:
                    parent_term = line[6:16]
                    parents[current_term].add(parent_term)
            elif line.startswith("relationship: part_of "):
                if current_term:
                    parts = line.split()
                    if len(parts) >= 3 and parts[2].startswith("GO:"):
                        parents[current_term].add(parts[2][:10])

    return dict(parents)


def read_test_labels(
    test_labels: Path,
    excluded_protein_ids: set[str],
) -> tuple[list[str], dict[str, list[str]], dict[str, int | str]]:
    header_lines = []
    rows_by_protein_id = defaultdict(list)
    input_rows = 0
    skipped_homologous = 0
    skipped_invalid = 0

    with test_labels.open("r", encoding="utf-8") as file:
        for line in file:
            raw_line = line.rstrip("\n")
            if not raw_line:
                continue

            parts = raw_line.split("\t")
            if is_tsv_header(parts):
                header_lines.append(raw_line)
                continue

            input_rows += 1
            if len(parts) < 3:
                skipped_invalid += 1
                continue

            protein_id = parts[0].strip()
            go_term = parts[1].strip()
            aspect = parts[2].strip().lower()
            if not protein_id or not go_term.startswith("GO:") or aspect not in ASPECTS:
                skipped_invalid += 1
                continue

            if protein_id in excluded_protein_ids:
                skipped_homologous += 1
                continue

            rows_by_protein_id[protein_id].append(raw_line)

    summary = {
        "test_label_input_rows": input_rows,
        "test_label_rows_after_homology_filter": sum(len(rows) for rows in rows_by_protein_id.values()),
        "skipped_label_rows_homologous": skipped_homologous,
        "skipped_label_rows_invalid": skipped_invalid,
        "test_proteins_with_labels_after_homology_filter": len(rows_by_protein_id),
        "test_label_header_rows": len(header_lines),
    }
    return header_lines, dict(rows_by_protein_id), summary


def parse_label_row(row: str) -> tuple[str, str, str] | None:
    parts = row.split("\t")
    if len(parts) < 3:
        return None

    protein_id = parts[0].strip()
    go_term = parts[1].strip()
    aspect = parts[2].strip().lower()
    if not protein_id or not go_term.startswith("GO:") or aspect not in ASPECTS:
        return None

    return protein_id, go_term, aspect


def propagate_labels(
    rows_by_protein_id: dict[str, list[str]],
    parents: dict[str, set[str]],
) -> tuple[dict[str, list[str]], dict[str, int]]:
    propagated_rows_by_protein_id = {}
    original_annotations = 0
    propagated_annotations = 0

    for protein_id, rows in rows_by_protein_id.items():
        labels_by_aspect = defaultdict(set)
        for row in rows:
            parsed = parse_label_row(row)
            if parsed is None:
                continue
            _, go_term, aspect = parsed
            labels_by_aspect[aspect].add(go_term)

        propagated_labels = set()
        for aspect, terms in labels_by_aspect.items():
            queue = list(terms)
            seen = set()
            while queue:
                term = queue.pop(0)
                if term in seen:
                    continue
                seen.add(term)

                propagated_labels.add((term, aspect))
                queue.extend(parents.get(term, []))

        original_annotations += sum(len(terms) for terms in labels_by_aspect.values())
        propagated_annotations += len(propagated_labels)
        propagated_rows_by_protein_id[protein_id] = [
            f"{protein_id}\t{go_term}\t{aspect}"
            for go_term, aspect in sorted(propagated_labels, key=lambda item: (item[1], item[0]))
        ]

    return propagated_rows_by_protein_id, {
        "labels_before_propagation": original_annotations,
        "labels_after_propagation": propagated_annotations,
        "propagated_label_rows_added": propagated_annotations - original_annotations,
    }


def filter_labels_to_train_space(
    rows_by_protein_id: dict[str, list[str]],
    train_label_space: set[str],
) -> tuple[dict[str, list[str]], dict[str, int]]:
    filtered_rows_by_protein_id = {}
    input_rows = 0
    output_rows = 0
    skipped_out_of_label_space = 0

    for protein_id, rows in rows_by_protein_id.items():
        kept_rows = []
        for row in rows:
            input_rows += 1
            parsed = parse_label_row(row)
            if parsed is None:
                continue

            _, go_term, _ = parsed
            if go_term not in train_label_space:
                skipped_out_of_label_space += 1
                continue

            kept_rows.append(row)
            output_rows += 1

        if kept_rows:
            filtered_rows_by_protein_id[protein_id] = kept_rows

    return filtered_rows_by_protein_id, {
        "label_rows_before_label_space_filter": input_rows,
        "label_rows_after_label_space_filter": output_rows,
        "skipped_label_rows_out_of_label_space": skipped_out_of_label_space,
        "proteins_with_labels_after_label_space_filter": len(filtered_rows_by_protein_id),
    }


def write_labels(header_lines: list[str], rows_by_protein_id: dict[str, list[str]], path: Path) -> int:
    output_rows = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as file:
        for header_line in header_lines:
            file.write(f"{header_line}\n")
        for protein_id in sorted(rows_by_protein_id):
            for row in rows_by_protein_id[protein_id]:
                file.write(f"{row}\n")
                output_rows += 1
    return output_rows


def run_mmseqs_search(
    mmseqs_path: str,
    query_fasta: Path,
    target_fasta: Path,
    threshold: float,
    coverage: float,
    cov_mode: int,
) -> set[str]:
    with tempfile.TemporaryDirectory(prefix="mmseqs_filter_") as work_dir:
        work_path = Path(work_dir)
        result_path = work_path / "matches.tsv"
        tmp_path = work_path / "tmp"
        command = [
            mmseqs_path,
            "easy-search",
            str(query_fasta),
            str(target_fasta),
            str(result_path),
            str(tmp_path),
            "--min-seq-id",
            str(threshold),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
            "--format-output",
            "query,target,pident,alnlen,evalue,bits",
        ]
        subprocess.run(command, check=True)

        if not result_path.exists():
            return set()

        homologous_query_ids = set()
        with result_path.open("r", encoding="utf-8") as file:
            for line in file:
                if not line.strip():
                    continue
                homologous_query_ids.add(fasta_record_id(line.split("\t", 1)[0]))

    return homologous_query_ids


def filter_fasta_by_ids(
    test_fasta: Path,
    output_fasta: Path,
    retained_protein_ids: set[str],
) -> dict[str, int | str]:
    input_records = read_fasta(test_fasta)
    output_records = [
        (header, sequence)
        for header, sequence in input_records
        if fasta_record_id(header) in retained_protein_ids
    ]
    write_fasta(output_records, output_fasta)

    input_protein_ids = {fasta_record_id(header) for header, _ in input_records}
    missing_fasta_ids = sorted(retained_protein_ids - input_protein_ids)

    return {
        "test_fasta_input_records": len(input_records),
        "test_fasta_output_records": len(output_records),
        "removed_fasta_records": len(input_records) - len(output_records),
        "label_proteins_missing_from_fasta": len(missing_fasta_ids),
        "output_fasta": str(output_fasta),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a protein function prediction test set by train homology and train label space."
        )
    )
    parser.add_argument("--train-fasta", required=True, help="training FASTA path")
    parser.add_argument("--train-labels", required=True, help="training TSV label path")
    parser.add_argument("--threshold", type=float, required=True, help="MMseqs2 --min-seq-id, range [0, 1]")
    parser.add_argument("--test-fasta", required=True, help="test FASTA path")
    parser.add_argument("--test-labels", required=True, help="test TSV label path")
    parser.add_argument("--out-fasta", required=True, help="filtered FASTA output path")
    parser.add_argument("--out-labels", required=True, help="filtered TSV label output path")
    parser.add_argument("--restricted", action="store_true", help="Restrict test labels to the training label space")
    parser.add_argument(
        "--obo",
        default=str(Path(__file__).with_name("go-basic.obo")),
        help="GO OBO path used for label propagation, default: fasta/go-basic.obo",
    )
    parser.add_argument("--coverage", type=float, default=0.0, help="MMseqs2 coverage threshold, default 0.0")
    parser.add_argument("--cov-mode", type=int, default=0, help="MMseqs2 --cov-mode, default 0")
    parser.add_argument("--mmseqs", default="mmseqs", help="mmseqs executable path")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not 0.0 <= args.threshold <= 1.0:
        raise ValueError("--threshold must be between 0 and 1")
    if not 0.0 <= args.coverage <= 1.0:
        raise ValueError("--coverage must be between 0 and 1")

    for path in [
        Path(args.train_fasta),
        Path(args.train_labels),
        Path(args.test_fasta),
        Path(args.test_labels),
        Path(args.obo),
    ]:
        if not path.exists():
            raise FileNotFoundError(path)


def main() -> None:
    args = parse_args()
    validate_args(args)

    train_fasta = Path(args.train_fasta)
    train_labels = Path(args.train_labels)
    test_fasta = Path(args.test_fasta)
    test_labels = Path(args.test_labels)
    out_fasta = Path(args.out_fasta)
    out_labels = Path(args.out_labels)
    obo_path = Path(args.obo)

    mmseqs_path = shutil.which(args.mmseqs)
    if mmseqs_path is None:
        raise FileNotFoundError(
            f"Cannot find MMseqs2 executable: {args.mmseqs}. Install MMseqs2 or pass --mmseqs."
        )

    train_label_space = read_label_space(train_labels)
    if not train_label_space:
        raise ValueError(f"No GO labels found in train labels: {train_labels}")
    parents = load_obo_parents(obo_path)

    homologous_test_ids = run_mmseqs_search(
        mmseqs_path=mmseqs_path,
        query_fasta=test_fasta,
        target_fasta=train_fasta,
        threshold=args.threshold,
        coverage=args.coverage,
        cov_mode=args.cov_mode,
    )

    header_lines, rows_by_protein_id, label_summary = read_test_labels(
        test_labels=test_labels,
        excluded_protein_ids=homologous_test_ids,
    )
    rows_by_protein_id, propagation_summary = propagate_labels(
        rows_by_protein_id=rows_by_protein_id,
        parents=parents,
    )
    if args.restricted:
        rows_by_protein_id, label_space_summary = filter_labels_to_train_space(
            rows_by_protein_id=rows_by_protein_id,
            train_label_space=train_label_space,
        )
    retained_protein_ids = set(rows_by_protein_id)

    output_label_rows = write_labels(header_lines, rows_by_protein_id, out_labels)
    fasta_summary = filter_fasta_by_ids(test_fasta, out_fasta, retained_protein_ids)

    print(f"train_fasta: {train_fasta}")
    print(f"train_labels: {train_labels}")
    print(f"test_fasta: {test_fasta}")
    print(f"test_labels: {test_labels}")
    print(f"obo: {obo_path}")
    print(f"threshold: {args.threshold:g}")
    print(f"coverage: {args.coverage:g}")
    print(f"train_label_space: {len(train_label_space)}")
    print(f"homologous_test_proteins: {len(homologous_test_ids)}")
    for key, value in label_summary.items():
        print(f"{key}: {value}")
    for key, value in propagation_summary.items():
        print(f"{key}: {value}")
    if args.restricted:
        for key, value in label_space_summary.items():
            print(f"{key}: {value}")
    print(f"removed_unlabeled_test_proteins: {label_summary['test_proteins_with_labels_after_homology_filter'] - len(retained_protein_ids)}")
    print(f"output_label_rows: {output_label_rows}")
    print(f"output_labels: {out_labels}")
    for key, value in fasta_summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
