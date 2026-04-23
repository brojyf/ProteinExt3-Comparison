from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path


ASPECTS = ["p", "f", "c"]
ASPECT_NAMES = {
    "p": "BP",
    "f": "MF",
    "c": "CC",
}


def read_fasta_ids(path: Path) -> set[str]:
    protein_ids = set()
    with path.open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                continue
            protein_ids.add(fasta_record_id(line))
    return protein_ids


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


def read_labels(path: Path) -> tuple[dict[str, set[tuple[str, str]]], int, int]:
    labels_by_protein = defaultdict(set)
    total_rows = 0
    skipped_rows = 0

    with path.open("r", encoding="utf-8") as file:
        for line in file:
            raw_line = line.rstrip("\n")
            if not raw_line:
                continue

            parts = raw_line.split("\t")
            if is_tsv_header(parts):
                continue

            total_rows += 1
            if len(parts) < 3:
                skipped_rows += 1
                continue

            protein_id = parts[0].strip()
            go_term = parts[1].strip()
            aspect = parts[2].strip().lower()
            if not protein_id or not go_term.startswith("GO:") or aspect not in ASPECTS:
                skipped_rows += 1
                continue

            labels_by_protein[protein_id].add((go_term, aspect))

    return dict(labels_by_protein), total_rows, skipped_rows


def safe_divide(numerator: int | float, denominator: int | float) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def analyze(fasta_path: Path, labels_path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    fasta_protein_ids = read_fasta_ids(fasta_path)
    labels_by_protein, total_label_rows, skipped_label_rows = read_labels(labels_path)
    labeled_protein_ids = set(labels_by_protein)

    labels_by_aspect = {aspect: [] for aspect in ASPECTS}
    proteins_by_aspect = {aspect: set() for aspect in ASPECTS}
    terms_by_aspect = {aspect: set() for aspect in ASPECTS}
    term_frequency = Counter()
    term_aspect = {}

    for protein_id, labels in labels_by_protein.items():
        for go_term, aspect in labels:
            labels_by_aspect[aspect].append((protein_id, go_term))
            proteins_by_aspect[aspect].add(protein_id)
            terms_by_aspect[aspect].add(go_term)
            term_frequency[go_term] += 1
            term_aspect[go_term] = aspect

    total_annotations = sum(len(labels) for labels in labels_by_protein.values())
    total_label_space = len(term_frequency)
    labeled_proteins_in_fasta = labeled_protein_ids & fasta_protein_ids

    summary = {
        "fasta": str(fasta_path),
        "labels": str(labels_path),
        "fasta_proteins": str(len(fasta_protein_ids)),
        "labeled_proteins": str(len(labeled_protein_ids)),
        "labeled_proteins_in_fasta": str(len(labeled_proteins_in_fasta)),
        "labeled_proteins_missing_from_fasta": str(len(labeled_protein_ids - fasta_protein_ids)),
        "fasta_proteins_without_labels": str(len(fasta_protein_ids - labeled_protein_ids)),
        "input_label_rows": str(total_label_rows),
        "valid_unique_annotations": str(total_annotations),
        "skipped_label_rows": str(skipped_label_rows),
        "annotation_density": f"{safe_divide(total_annotations, len(labeled_protein_ids)):.6f}",
        "label_space_size": str(total_label_space),
    }

    for aspect in ASPECTS:
        aspect_name = ASPECT_NAMES[aspect]
        annotation_count = len(labels_by_aspect[aspect])
        protein_count = len(proteins_by_aspect[aspect])
        term_count = len(terms_by_aspect[aspect])
        summary[f"{aspect_name}_annotations"] = str(annotation_count)
        summary[f"{aspect_name}_annotation_ratio"] = f"{safe_divide(annotation_count, total_annotations):.6f}"
        summary[f"{aspect_name}_proteins"] = str(protein_count)
        summary[f"{aspect_name}_protein_ratio"] = f"{safe_divide(protein_count, len(labeled_protein_ids)):.6f}"
        summary[f"{aspect_name}_label_space_size"] = str(term_count)
        summary[f"{aspect_name}_annotation_density"] = f"{safe_divide(annotation_count, protein_count):.6f}"

    frequency_rows = []
    for rank, (go_term, count) in enumerate(term_frequency.most_common(), start=1):
        aspect = term_aspect[go_term]
        frequency_rows.append(
            {
                "rank": str(rank),
                "go_term": go_term,
                "aspect": ASPECT_NAMES[aspect],
                "count": str(count),
                "protein_ratio": f"{safe_divide(count, len(labeled_protein_ids)):.6f}",
                "annotation_ratio": f"{safe_divide(count, total_annotations):.6f}",
            }
        )

    return summary, frequency_rows


def write_report(path: Path, summary: dict[str, str], frequency_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as file:
        file.write("# Dataset Analysis\n\n")
        file.write("## Inputs\n\n")
        file.write("| Field | Value |\n")
        file.write("|---|---|\n")
        file.write(f"| FASTA | `{summary['fasta']}` |\n")
        file.write(f"| Labels | `{summary['labels']}` |\n")

        file.write("\n## Dataset Summary\n\n")
        file.write("| Metric | Value |\n")
        file.write("|---|---:|\n")
        for key in [
            "fasta_proteins",
            "labeled_proteins",
            "labeled_proteins_in_fasta",
            "labeled_proteins_missing_from_fasta",
            "fasta_proteins_without_labels",
            "input_label_rows",
            "valid_unique_annotations",
            "skipped_label_rows",
            "annotation_density",
            "label_space_size",
        ]:
            file.write(f"| {key} | {summary[key]} |\n")

        file.write("\n## Aspect Distribution\n\n")
        file.write(
            "| Aspect | Annotations | Annotation Ratio | Proteins | Protein Ratio | "
            "Label Space Size | Annotation Density |\n"
        )
        file.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for aspect_name in ["BP", "MF", "CC"]:
            file.write(
                f"| {aspect_name} | "
                f"{summary[f'{aspect_name}_annotations']} | "
                f"{summary[f'{aspect_name}_annotation_ratio']} | "
                f"{summary[f'{aspect_name}_proteins']} | "
                f"{summary[f'{aspect_name}_protein_ratio']} | "
                f"{summary[f'{aspect_name}_label_space_size']} | "
                f"{summary[f'{aspect_name}_annotation_density']} |\n"
            )

        file.write("\n## GO Term Frequency Long Tail\n\n")
        file.write("| Rank | GO Term | Aspect | Count | Protein Ratio | Annotation Ratio |\n")
        file.write("|---:|---|---|---:|---:|---:|\n")
        for row in frequency_rows:
            file.write(
                f"| {row['rank']} | {row['go_term']} | {row['aspect']} | {row['count']} | "
                f"{row['protein_ratio']} | {row['annotation_ratio']} |\n"
            )


def print_summary(summary: dict[str, str]) -> None:
    for key, value in summary.items():
        print(f"{key}: {value}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze FASTA/label TSV statistics for protein function prediction datasets."
    )
    parser.add_argument("--fasta", required=True, help="input FASTA path")
    parser.add_argument("--labels", required=True, help="input label TSV path: protein_id, GO term, aspect")
    parser.add_argument(
        "--out",
        help="optional Markdown report path containing summary and full GO term frequency output",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fasta_path = Path(args.fasta)
    labels_path = Path(args.labels)

    for path in [fasta_path, labels_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    summary, frequency_rows = analyze(fasta_path, labels_path)
    print_summary(summary)

    if args.out:
        output_path = Path(args.out)
        write_report(output_path, summary, frequency_rows)
        print(f"report_output: {output_path}")


if __name__ == "__main__":
    main()
