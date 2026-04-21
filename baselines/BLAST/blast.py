from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


def extract_entry_id(header: str) -> str:
    title = header[1:] if header.startswith(">") else header
    parts = title.split("|")
    if len(parts) >= 2:
        return parts[1]
    return title.split()[0]


def read_fasta_headers(path: Path) -> list[str]:
    headers = []
    with path.open("r", encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                headers.append(line.strip())
    return headers


def read_base_truth(path: Path) -> dict[str, set[tuple[str, str]]]:
    truth_by_pid: dict[str, set[tuple[str, str]]] = defaultdict(set)

    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line:
                continue

            parts = [part.strip() for part in line.split("\t")]
            if len(parts) < 3:
                continue

            if parts[0] in {"AUTHOR", "MODEL", "KEYWORDS", "END"}:
                continue
            if parts[0] == "EntryID" and len(parts) >= 3 and parts[1] == "term":
                continue

            protein_id, go_term, aspect = parts[:3]
            if not protein_id or not go_term.startswith("GO:") or aspect not in {"c", "f", "p", "C", "F", "P"}:
                continue

            truth_by_pid[protein_id].add((go_term, aspect.lower()))

    if not truth_by_pid:
        raise ValueError(f"无法从 {path} 解析任何 pid/go_term/aspect 记录")

    return truth_by_pid


def blast_db_files(prefix: Path) -> list[Path]:
    return [prefix.with_suffix(suffix) for suffix in [".pin", ".phr", ".psq"]]


def blast_db_ready(prefix: Path, fasta_path: Path) -> bool:
    db_files = blast_db_files(prefix)
    if not all(path.exists() for path in db_files):
        return False
    fasta_mtime = fasta_path.stat().st_mtime
    return all(path.stat().st_mtime >= fasta_mtime for path in db_files)


def ensure_blast_tools() -> tuple[str, str]:
    makeblastdb_path = shutil.which("makeblastdb")
    blastp_path = shutil.which("blastp")
    if makeblastdb_path is None or blastp_path is None:
        raise FileNotFoundError("需要先安装 NCBI BLAST+，并确保 makeblastdb 和 blastp 在 PATH 中。")
    return makeblastdb_path, blastp_path


def build_blast_db(makeblastdb_path: str, base_fasta: Path, db_prefix: Path, reuse_existing_db: bool) -> None:
    if reuse_existing_db and blast_db_ready(db_prefix, base_fasta):
        return

    cmd = [
        makeblastdb_path,
        "-in",
        str(base_fasta),
        "-dbtype",
        "prot",
        "-out",
        str(db_prefix),
    ]
    subprocess.run(cmd, check=True)


def run_blastp(
    blastp_path: str,
    query_fasta: Path,
    db_prefix: Path,
    blast_output_tsv: Path,
    num_threads: int,
    evalue: float,
    max_target_seqs: int,
) -> None:
    cmd = [
        blastp_path,
        "-query",
        str(query_fasta),
        "-db",
        str(db_prefix),
        "-out",
        str(blast_output_tsv),
        "-evalue",
        str(evalue),
        "-num_threads",
        str(num_threads),
        "-max_target_seqs",
        str(max_target_seqs),
        "-max_hsps",
        "1",
        "-outfmt",
        "6 qseqid sseqid bitscore pident evalue",
    ]
    subprocess.run(cmd, check=True)


def parse_blast_hits(path: Path) -> list[dict[str, object]]:
    if not path.exists() or path.stat().st_size == 0:
        return []

    hits = []
    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != 5:
                continue

            query_header, subject_header, bitscore, pident, evalue = parts
            hits.append(
                {
                    "query_header": query_header,
                    "subject_header": subject_header,
                    "query_id": extract_entry_id(query_header),
                    "subject_id": extract_entry_id(subject_header),
                    "bitscore": float(bitscore),
                    "pident": float(pident),
                    "evalue": float(evalue),
                }
            )

    return hits


def keep_best_subject_hit(hits: list[dict[str, object]]) -> list[dict[str, object]]:
    best_hits: dict[tuple[str, str], dict[str, object]] = {}

    for hit in hits:
        key = (str(hit["query_id"]), str(hit["subject_id"]))
        previous = best_hits.get(key)
        if previous is None:
            best_hits[key] = hit
            continue

        rank = (
            float(hit["bitscore"]),
            float(hit["pident"]),
            -float(hit["evalue"]),
        )
        previous_rank = (
            float(previous["bitscore"]),
            float(previous["pident"]),
            -float(previous["evalue"]),
        )
        if rank > previous_rank:
            best_hits[key] = hit

    return list(best_hits.values())


def select_hits_for_prediction(
    hits: list[dict[str, object]], min_bitscore: float, top_hits_for_terms: int
) -> list[dict[str, object]]:
    hits_by_query: dict[str, list[dict[str, object]]] = defaultdict(list)

    for hit in hits:
        if float(hit["bitscore"]) < min_bitscore:
            continue
        hits_by_query[str(hit["query_id"])].append(hit)

    selected_hits = []
    for query_id, query_hits in hits_by_query.items():
        del query_id
        ranked_hits = sorted(query_hits, key=lambda item: float(item["bitscore"]), reverse=True)
        selected_hits.extend(ranked_hits[:top_hits_for_terms])

    return selected_hits


def propagate_terms(
    hits: list[dict[str, object]], base_truth: dict[str, set[tuple[str, str]]]
) -> list[tuple[str, str, float, str]]:
    best_score: dict[tuple[str, str], float] = {}
    aspect_map: dict[tuple[str, str], str] = {}

    for hit in hits:
        query_id = str(hit["query_id"])
        subject_id = str(hit["subject_id"])
        score = float(hit["pident"]) / 100.0
        for go_term, aspect in base_truth.get(subject_id, set()):
            key = (query_id, go_term)
            if key not in best_score or score > best_score[key]:
                best_score[key] = score
                aspect_map[key] = aspect

    return sorted(
        (query_id, go_term, best_score[(query_id, go_term)], aspect_map[(query_id, go_term)])
        for query_id, go_term in best_score
    )


def write_predictions(path: Path, predictions: list[tuple[str, str, float, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as file:
        for protein_id, go_term, score, aspect in predictions:
            file.write(f"{protein_id}\t{go_term}\t{score:.6f}\t{aspect}\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-fasta", required=True)
    parser.add_argument("--base-truth", required=True)
    parser.add_argument("--in", dest="input_fasta", required=True)
    parser.add_argument("--out", default="baselines/BLAST/prediction.tsv")
    parser.add_argument("--num-threads", type=int, default=max(1, (os.cpu_count() or 1) - 1))
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--max-target-seqs", type=int, default=20)
    parser.add_argument("--min-bitscore", type=float, default=50.0)
    parser.add_argument("--top-hits-for-terms", type=int, default=10)
    parser.add_argument("--reuse-existing-db", action="store_true", default=False)
    return parser.parse_args()


def main():
    args = parse_args()

    base_fasta = Path(args.base_fasta)
    base_truth = Path(args.base_truth)
    input_fasta = Path(args.input_fasta)
    output_tsv = Path(args.out)

    for path in [base_fasta, base_truth, input_fasta]:
        if not path.exists():
            raise FileNotFoundError(f"找不到文件: {path}")

    makeblastdb_path, blastp_path = ensure_blast_tools()

    work_dir = output_tsv.parent / f"{output_tsv.stem}_blast_work"
    db_prefix = work_dir / "base_db"
    blast_hits_tsv = work_dir / "blast_hits.tsv"
    work_dir.mkdir(parents=True, exist_ok=True)

    build_blast_db(makeblastdb_path, base_fasta, db_prefix, args.reuse_existing_db)
    run_blastp(
        blastp_path=blastp_path,
        query_fasta=input_fasta,
        db_prefix=db_prefix,
        blast_output_tsv=blast_hits_tsv,
        num_threads=max(1, int(args.num_threads)),
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
    )

    base_truth_map = read_base_truth(base_truth)
    raw_hits = parse_blast_hits(blast_hits_tsv)
    best_hits = keep_best_subject_hit(raw_hits)
    selected_hits = select_hits_for_prediction(
        best_hits,
        min_bitscore=args.min_bitscore,
        top_hits_for_terms=args.top_hits_for_terms,
    )
    predictions = propagate_terms(selected_hits, base_truth_map)
    write_predictions(output_tsv, predictions)

    query_ids = [extract_entry_id(header) for header in read_fasta_headers(input_fasta)]
    predicted_query_ids = {protein_id for protein_id, _, _, _ in predictions}
    missing_queries = [query_id for query_id in query_ids if query_id not in predicted_query_ids]

    print(f"base_fasta: {base_fasta}")
    print(f"base_truth: {base_truth}")
    print(f"input_fasta: {input_fasta}")
    print(f"output_tsv: {output_tsv}")
    print(f"work_dir: {work_dir}")
    print(f"raw_hits: {len(raw_hits)}")
    print(f"best_hits: {len(best_hits)}")
    print(f"selected_hits: {len(selected_hits)}")
    print(f"prediction_rows: {len(predictions)}")
    print(f"query_proteins: {len(query_ids)}")
    print(f"queries_with_predictions: {len(predicted_query_ids)}")
    print(f"queries_without_predictions: {len(missing_queries)}")
    if missing_queries:
        print(f"first_missing_query: {missing_queries[0]}")


if __name__ == "__main__":
    main()
