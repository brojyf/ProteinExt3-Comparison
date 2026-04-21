from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


ASPECTS = ["c", "f", "p"]
ASPECT_NAMES = {
    "c": "cc",
    "f": "mf",
    "p": "bp",
}

def load_obo(path: Path):
    parents = defaultdict(set)
    with path.open("r", encoding="utf-8") as f:
        current_term = None
        for raw_line in f:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
            elif line.startswith("id: "):
                if not current_term:
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
    return parents


def propagate_truth(truth_by_aspect, go_to_aspect, parents):
    for aspect in ASPECTS:
        for protein_id, terms in list(truth_by_aspect[aspect].items()):
            propagated_terms = set()
            queue = list(terms)
            while queue:
                term = queue.pop(0)
                if term not in propagated_terms:
                    propagated_terms.add(term)
                    queue.extend(parents.get(term, []))
            truth_by_aspect[aspect][protein_id] = propagated_terms
            for term in propagated_terms:
                go_to_aspect[term] = aspect


def propagate_predictions(predictions_by_aspect, parents):
    for aspect in ASPECTS:
        for protein_id, terms_dict in predictions_by_aspect[aspect].items():
            propagated_dict = dict(terms_dict)
            queue = list(propagated_dict.keys())
            while queue:
                term = queue.pop(0)
                score = propagated_dict[term]
                for p in parents.get(term, []):
                    if p not in propagated_dict or propagated_dict[p] < score:
                        propagated_dict[p] = score
                        queue.append(p)
            predictions_by_aspect[aspect][protein_id] = propagated_dict



METHOD_NAME_BY_STEM = {
    "predictions": None,
    "blast": "BLAST",
    "naive": "Naive",
    "dgs": "DeepGO-SE",
    "esm2_mlp": "ESM2+MLP",
    "pe3": "ProteinExt3",
    "pe3_esm2": "ProteinExt3 (ESM2)",
    "pe3_t5": "ProteinExt3 (ProtT5)",
    "pe3_cnn": "ProteinExt3 (CNN)",
    "predictions_esm2": "ProteinExt3 (ESM2)",
    "predictions_t5": "ProteinExt3 (ProtT5)",
    "predictions_cnn": "ProteinExt3 (CNN)",
    "predictions_blast": "ProteinExt3 (BLAST)",
    "predictions_blastboosted": "ProteinExt3 (Boosted)",
}


def normalize_protein_id(text: str) -> str:
    text = text.strip()
    parts = text.split("|")
    if len(parts) >= 3 and parts[1]:
        return parts[1]
    return text.split()[0]


def detect_method_name(path: Path) -> str:
    stem = path.stem.lower()
    parent_name = path.parent.name
    mapped = METHOD_NAME_BY_STEM.get(stem)
    if mapped:
        return mapped
    if stem == "predictions":
        return parent_name
    return stem


def iter_tsv_rows(path: Path):
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
            if parts[0].lower() in {"entryid", "entry", "protein_id", "protein", "query_id", "pid"}:
                continue

            yield parts


def load_truth(path: Path):
    truth_by_aspect = {aspect: defaultdict(set) for aspect in ASPECTS}
    go_to_aspect = {}

    for parts in iter_tsv_rows(path):
        if len(parts) < 3:
            raise ValueError(f"truth 文件需要三列: pid go_term aspect, 出错文件: {path}")

        protein_id = normalize_protein_id(parts[0])
        go_term = parts[1]
        aspect = parts[2].lower()

        if aspect not in ASPECTS or not go_term.startswith("GO:"):
            continue

        truth_by_aspect[aspect][protein_id].add(go_term)
        go_to_aspect[go_term] = aspect

    return truth_by_aspect, go_to_aspect


def load_predictions(path: Path, go_to_aspect: dict[str, str]):
    predictions_by_aspect = {aspect: defaultdict(dict) for aspect in ASPECTS}

    for parts in iter_tsv_rows(path):
        protein_id = normalize_protein_id(parts[0])
        go_term = parts[1]
        if not go_term.startswith("GO:"):
            continue

        aspect = None
        score = 1.0

        if len(parts) >= 3:
            token = parts[2].strip().lower()
            if token in ASPECTS:
                aspect = token
                if len(parts) >= 4:
                    try:
                        score = float(parts[3])
                    except ValueError:
                        score = 1.0
            else:
                try:
                    score = float(parts[2])
                except ValueError:
                    score = 1.0
                if len(parts) >= 4:
                    token4 = parts[3].strip().lower()
                    if token4 in ASPECTS:
                        aspect = token4

        if aspect is None:
            aspect = go_to_aspect.get(go_term)
        if aspect not in ASPECTS:
            continue
        if go_term not in go_to_aspect:
            continue

        old_score = predictions_by_aspect[aspect][protein_id].get(go_term)
        if old_score is None or score > old_score:
            predictions_by_aspect[aspect][protein_id][go_term] = score

    return predictions_by_aspect


def compute_ic_by_aspect(truth_by_aspect):
    ic_by_aspect = {aspect: {} for aspect in ASPECTS}

    for aspect in ASPECTS:
        truth_by_protein = truth_by_aspect[aspect]
        protein_count = len(truth_by_protein)
        if protein_count == 0:
            continue

        proteins_by_term = defaultdict(set)
        for protein_id, terms in truth_by_protein.items():
            for term in terms:
                proteins_by_term[term].add(protein_id)

        for term, proteins in proteins_by_term.items():
            frequency = len(proteins)
            ic_by_aspect[aspect][term] = -math.log2(frequency / protein_count)

    return ic_by_aspect


def compute_binary_curve_metrics(labels: list[int], scores: list[float]) -> tuple[float, float]:
    total_pos = sum(labels)
    total_neg = len(labels) - total_pos
    if total_pos == 0 or total_neg == 0:
        return 0.0, 0.0

    items = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    tp = 0
    fp = 0
    prev_recall = 0.0
    prev_fpr = 0.0
    prev_tpr = 0.0
    aupr = 0.0
    auc = 0.0

    index = 0
    while index < len(items):
        score = items[index][0]
        group_tp = 0
        group_fp = 0

        while index < len(items) and items[index][0] == score:
            group_tp += items[index][1]
            group_fp += 1 - items[index][1]
            index += 1

        tp += group_tp
        fp += group_fp

        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / total_pos
        aupr += (recall - prev_recall) * precision
        prev_recall = recall

        fpr = fp / total_neg
        tpr = recall
        auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2
        prev_fpr = fpr
        prev_tpr = tpr

    return aupr, auc


def compute_curve_metrics(truth_by_protein, predictions_by_protein, min_positives):
    proteins = sorted(truth_by_protein)
    proteins_by_term = defaultdict(set)

    for protein_id, true_terms in truth_by_protein.items():
        for term in true_terms:
            proteins_by_term[term].add(protein_id)

    evaluated_terms = [
        term
        for term, positive_proteins in proteins_by_term.items()
        if len(positive_proteins) >= min_positives and len(positive_proteins) < len(proteins)
    ]

    if not evaluated_terms:
        return 0.0, 0.0

    aupr_values = []
    auc_values = []
    for term in sorted(evaluated_terms):
        positive_proteins = proteins_by_term[term]
        labels = [1 if protein_id in positive_proteins else 0 for protein_id in proteins]
        scores = [predictions_by_protein.get(protein_id, {}).get(term, 0.0) for protein_id in proteins]
        aupr, auc = compute_binary_curve_metrics(labels, scores)
        aupr_values.append(aupr)
        auc_values.append(auc)

    return sum(aupr_values) / len(aupr_values), sum(auc_values) / len(auc_values)


def sweep_threshold_metrics(truth_by_protein, predictions_by_protein, ic_by_term):
    proteins = sorted(truth_by_protein)
    protein_count = len(proteins)
    if protein_count == 0:
        return 0.0, 0.0, 0.0

    prediction_events = []
    for protein_id in proteins:
        for term, score in predictions_by_protein.get(protein_id, {}).items():
            prediction_events.append((score, protein_id, term))
    prediction_events.sort(key=lambda item: item[0], reverse=True)

    if not prediction_events:
        ru_total = sum(sum(ic_by_term.get(term, 0.0) for term in truth_by_protein[p]) for p in proteins)
        return 0.0, ru_total / protein_count, 0.0

    predicted_terms_by_protein = {protein_id: set() for protein_id in proteins}
    tp_count = {protein_id: 0 for protein_id in proteins}
    pred_count = {protein_id: 0 for protein_id in proteins}

    predicted_protein_count = 0
    precision_sum = 0.0
    recall_sum = 0.0
    ru_total = sum(sum(ic_by_term.get(term, 0.0) for term in truth_by_protein[p]) for p in proteins)
    mi_total = 0.0

    best_fmax = 0.0
    best_fmax_threshold = 0.0
    best_smin = math.sqrt((ru_total / protein_count) ** 2 + (mi_total / protein_count) ** 2)

    index = 0
    while index < len(prediction_events):
        current_score = prediction_events[index][0]

        while index < len(prediction_events) and prediction_events[index][0] == current_score:
            _, protein_id, term = prediction_events[index]
            index += 1

            if term in predicted_terms_by_protein[protein_id]:
                continue

            old_tp = tp_count[protein_id]
            old_pred = pred_count[protein_id]
            old_precision = old_tp / old_pred if old_pred else 0.0
            truth_count = len(truth_by_protein[protein_id])
            old_recall = old_tp / truth_count if truth_count else 0.0

            predicted_terms_by_protein[protein_id].add(term)
            pred_count[protein_id] = old_pred + 1
            if old_pred == 0:
                predicted_protein_count += 1

            if term in truth_by_protein[protein_id]:
                tp_count[protein_id] = old_tp + 1
                ru_total -= ic_by_term.get(term, 0.0)
            else:
                mi_total += ic_by_term.get(term, 0.0)

            new_tp = tp_count[protein_id]
            new_pred = pred_count[protein_id]
            new_precision = new_tp / new_pred
            new_recall = new_tp / truth_count if truth_count else 0.0

            precision_sum += new_precision - old_precision
            recall_sum += new_recall - old_recall

        avg_precision = precision_sum / predicted_protein_count if predicted_protein_count else 0.0
        avg_recall = recall_sum / protein_count
        if avg_precision + avg_recall > 0:
            f_score = 2 * avg_precision * avg_recall / (avg_precision + avg_recall)
            if f_score > best_fmax:
                best_fmax = f_score
                best_fmax_threshold = current_score

        s_value = math.sqrt((ru_total / protein_count) ** 2 + (mi_total / protein_count) ** 2)
        best_smin = min(best_smin, s_value)

    return best_fmax, best_smin, best_fmax_threshold


def evaluate_aspect(truth_by_protein, predictions_by_protein, ic_by_term, min_positives):
    truth_by_protein = {protein_id: terms for protein_id, terms in truth_by_protein.items() if terms}

    fmax, smin, threshold = sweep_threshold_metrics(truth_by_protein, predictions_by_protein, ic_by_term)
    aupr, auc = compute_curve_metrics(truth_by_protein, predictions_by_protein, min_positives)
    return fmax, smin, aupr, auc, threshold


def discover_prediction_files(baselines_dir: Path, prediction_name: str) -> list[tuple[str, Path]]:
    discovered = []

    for baseline_dir in sorted(p for p in baselines_dir.iterdir() if p.is_dir()):
        tsv_files = sorted(p for p in baseline_dir.glob("*.tsv"))
        for tsv_file in tsv_files:
            discovered.append((detect_method_name(tsv_file), tsv_file))

    return discovered


def write_result_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["Method", "Fmax", "Threshold", "Smin", "AUPR", "AUC"])
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default="fasta/test.tsv", help="truth tsv，格式为 pid go_term aspect")
    parser.add_argument("--obo", default="fasta/go-basic.obo", help="go-basic.obo 路径")
    parser.add_argument("--baselines-dir", default="baselines")
    parser.add_argument("--prediction-name", default="predictions.tsv")
    parser.add_argument("--out-dir", default="comparison/results")
    parser.add_argument("--min-positives", type=int, default=10)
    return parser.parse_args()


def main():
    args = parse_args()

    truth_path = Path(args.truth)
    obo_path = Path(args.obo)
    baselines_dir = Path(args.baselines_dir)
    out_dir = Path(args.out_dir)

    if not truth_path.exists():
        raise FileNotFoundError(f"找不到 truth 文件: {truth_path}")
    if not obo_path.exists():
        raise FileNotFoundError(f"找不到 obo 文件: {obo_path}")
    if not baselines_dir.exists():
        raise FileNotFoundError(f"找不到 baselines 目录: {baselines_dir}")

    parents = load_obo(obo_path)
    truth_by_aspect, go_to_aspect = load_truth(truth_path)
    propagate_truth(truth_by_aspect, go_to_aspect, parents)
    
    ic_by_aspect = compute_ic_by_aspect(truth_by_aspect)
    prediction_files = discover_prediction_files(baselines_dir, args.prediction_name)

    if not prediction_files:
        raise FileNotFoundError(f"没有在 {baselines_dir} 下发现任何预测文件")

    rows_by_aspect = {aspect: [] for aspect in ASPECTS}

    for method_name, prediction_path in prediction_files:
        predictions_by_aspect = load_predictions(prediction_path, go_to_aspect)
        propagate_predictions(predictions_by_aspect, parents)
        print(f"evaluating: {method_name} <- {prediction_path}")

        for aspect in ASPECTS:
            fmax, smin, aupr, auc, threshold = evaluate_aspect(
                truth_by_aspect[aspect],
                predictions_by_aspect[aspect],
                ic_by_aspect[aspect],
                args.min_positives,
            )
            rows_by_aspect[aspect].append(
                {
                    "Method": method_name,
                    "Fmax": f"{fmax:.3f}",
                    "Threshold": f"{threshold:.4f}",
                    "Smin": f"{smin:.3f}",
                    "AUPR": f"{aupr:.3f}",
                    "AUC": f"{auc:.3f}",
                }
            )

    for aspect in ASPECTS:
        output_path = out_dir / f"{ASPECT_NAMES[aspect]}.csv"
        write_result_csv(output_path, rows_by_aspect[aspect])
        print(f"wrote: {output_path}")


if __name__ == "__main__":
    main()
