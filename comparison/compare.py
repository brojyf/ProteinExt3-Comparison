from __future__ import annotations

import argparse
import csv
import math
import pickle
import random
from collections import defaultdict, deque
from functools import lru_cache
from pathlib import Path

import numpy as np


ASPECTS = ["c", "f", "p"]
ASPECT_NAMES = {
    "c": "c",
    "f": "f",
    "p": "p",
}
NAMESPACE_TO_ASPECT = {
    "cellular_component": "c",
    "molecular_function": "f",
    "biological_process": "p",
}
THRESHOLDS = np.arange(0, 101, dtype=np.int16) / 100.0
IC_PKL_KEYS = {
    "c": ("C", "c", "CC", "cc"),
    "f": ("F", "f", "MF", "mf"),
    "p": ("P", "p", "BP", "bp"),
}


def parse_go_id(raw_text: str) -> str | None:
    token = raw_text.strip().split()[0]
    if token.startswith("GO:"):
        return token
    return None


def load_obo(path: Path):
    parents = defaultdict(set)
    go_to_aspect = {}

    with path.open("r", encoding="utf-8") as handle:
        current_term = None
        current_aspect = None
        in_term = False

        for raw_line in handle:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
                current_aspect = None
                in_term = True
            elif line.startswith("["):
                if current_term and current_aspect:
                    go_to_aspect[current_term] = current_aspect
                current_term = None
                current_aspect = None
                in_term = False
            elif not in_term:
                continue
            elif line.startswith("id: "):
                current_term = line[4:].strip()
            elif line.startswith("namespace: "):
                namespace = line[len("namespace: "):].strip()
                current_aspect = NAMESPACE_TO_ASPECT.get(namespace)
            elif line.startswith("is_a: "):
                if current_term:
                    parent_term = parse_go_id(line[6:])
                    if parent_term:
                        parents[current_term].add(parent_term)
            elif line.startswith("relationship: part_of "):
                if current_term:
                    parent_term = parse_go_id(line[len("relationship: part_of "):])
                    if parent_term:
                        parents[current_term].add(parent_term)

        if current_term and current_aspect:
            go_to_aspect[current_term] = current_aspect

    return parents, go_to_aspect


def normalize_protein_id(text: str) -> str:
    text = text.strip()
    parts = text.split("|")
    if len(parts) >= 3 and parts[1]:
        return parts[1]
    return text.split()[0]


def iter_tsv_rows(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
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


def build_ancestor_getter(parents, go_to_aspect):
    @lru_cache(maxsize=None)
    def get_ancestors(term: str, aspect: str) -> tuple[str, ...]:
        seen = set()
        queue = deque([term])

        while queue:
            cur = queue.popleft()
            for parent in parents.get(cur, ()):
                if parent in seen:
                    continue
                if go_to_aspect.get(parent) != aspect:
                    continue
                seen.add(parent)
                queue.append(parent)

        return tuple(seen)

    return get_ancestors


def propagate_truth(truth_by_aspect, get_ancestors):
    for aspect in ASPECTS:
        for protein_id, terms in list(truth_by_aspect[aspect].items()):
            propagated_terms = set(terms)
            for term in terms:
                propagated_terms.update(get_ancestors(term, aspect))
            truth_by_aspect[aspect][protein_id] = propagated_terms


def load_raw_predictions(path: Path, go_to_aspect: dict[str, str]):
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

        old_score = predictions_by_aspect[aspect][protein_id].get(go_term)
        if old_score is None or score > old_score:
            predictions_by_aspect[aspect][protein_id][go_term] = score

    return predictions_by_aspect


def propagate_predictions(predictions_by_aspect, get_ancestors):
    for aspect in ASPECTS:
        for protein_id, terms_dict in predictions_by_aspect[aspect].items():
            propagated_dict = dict(terms_dict)
            for term, score in terms_dict.items():
                for parent_term in get_ancestors(term, aspect):
                    old_score = propagated_dict.get(parent_term, 0.0)
                    if score > old_score:
                        propagated_dict[parent_term] = score
            predictions_by_aspect[aspect][protein_id] = propagated_dict


def load_eval_label_space(path: Path, get_ancestors, go_to_aspect):
    eval_truth_by_aspect, eval_go_to_aspect = load_truth(path)
    merged_go_to_aspect = dict(go_to_aspect)
    merged_go_to_aspect.update(eval_go_to_aspect)
    _ = merged_go_to_aspect  # 兼容原逻辑，实际 ancestor aspect 已由 get_ancestors 固定

    propagate_truth(eval_truth_by_aspect, get_ancestors)

    eval_terms_by_aspect = {aspect: set() for aspect in ASPECTS}
    for aspect in ASPECTS:
        for terms in eval_truth_by_aspect[aspect].values():
            eval_terms_by_aspect[aspect].update(terms)
    return eval_terms_by_aspect


def restrict_predictions_to_eval_space(predictions_by_aspect, eval_terms_by_aspect):
    restricted = {aspect: defaultdict(dict) for aspect in ASPECTS}

    for aspect in ASPECTS:
        allowed_terms = eval_terms_by_aspect[aspect]
        for protein_id, terms_dict in predictions_by_aspect[aspect].items():
            kept = {}
            for term, score in terms_dict.items():
                if term in allowed_terms:
                    old_score = kept.get(term)
                    if old_score is None or score > old_score:
                        kept[term] = score
            restricted[aspect][protein_id] = kept

    return restricted


def compute_ic_by_aspect(truth_by_aspect, go_to_aspect, parents):
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
            parent_terms = [
                parent_term
                for parent_term in parents.get(term, set())
                if go_to_aspect.get(parent_term) == aspect
            ]

            if not parent_terms:
                parent_count = protein_count
            else:
                parent_sets = [proteins_by_term.get(parent_term, set()) for parent_term in parent_terms]
                parent_count = len(set.intersection(*parent_sets)) if parent_sets else protein_count

            if parent_count == 0:
                ic_by_aspect[aspect][term] = 0.0
                continue

            conditional_probability = len(proteins) / parent_count
            if conditional_probability <= 0.0:
                ic_by_aspect[aspect][term] = 0.0
                continue

            ic_by_aspect[aspect][term] = -math.log2(conditional_probability)

    return ic_by_aspect


def load_ic_by_aspect_from_pickle(path: Path):
    with path.open("rb") as handle:
        payload = pickle.load(handle)

    if isinstance(payload, tuple):
        dict_payload = None
        for item in reversed(payload):
            if isinstance(item, dict):
                dict_payload = item
                break
        if dict_payload is None:
            raise ValueError(f"IC pickle 不包含可识别的 dict payload: {path}")
        payload = dict_payload

    if not isinstance(payload, dict):
        raise ValueError(f"IC pickle 顶层对象必须是 dict 或 tuple[... , dict]: {path}")

    ic_by_aspect = {aspect: {} for aspect in ASPECTS}
    for aspect in ASPECTS:
        aspect_payload = None
        for key in IC_PKL_KEYS[aspect]:
            if key in payload:
                aspect_payload = payload[key]
                break

        if aspect_payload is None:
            continue
        if not isinstance(aspect_payload, dict):
            raise ValueError(f"IC pickle 中 {aspect!r} 对应的数据不是 dict: {path}")

        normalized = {}
        for term, score in aspect_payload.items():
            if isinstance(term, str) and term.startswith("GO:"):
                normalized[term] = float(score)
        ic_by_aspect[aspect] = normalized

    return ic_by_aspect


def build_label_space(aspect, truth_by_protein, predictions_by_protein, eval_terms_by_aspect):
    if eval_terms_by_aspect is not None:
        return sorted(eval_terms_by_aspect[aspect])

    label_space = {term for terms in truth_by_protein.values() for term in terms}
    label_space.update(term for terms in predictions_by_protein.values() for term in terms)
    return sorted(label_space)


def build_matrices(truth_by_protein, predictions_by_protein, ic_by_term, aspect, eval_terms_by_aspect):
    truth_by_protein = {protein_id: terms for protein_id, terms in truth_by_protein.items() if terms}
    proteins = sorted(truth_by_protein)
    if not proteins:
        return None

    label_space = build_label_space(aspect, truth_by_protein, predictions_by_protein, eval_terms_by_aspect)
    term_to_idx = {term: i for i, term in enumerate(label_space)}

    n_proteins = len(proteins)
    n_terms = len(label_space)

    truth = np.zeros((n_proteins, n_terms), dtype=bool)
    pred = np.zeros((n_proteins, n_terms), dtype=np.float32)
    ic = np.zeros(n_terms, dtype=np.float32)
    explicit_pred_count = 0

    for term, idx in term_to_idx.items():
        ic[idx] = float(ic_by_term.get(term, 0.0))

    for i, protein_id in enumerate(proteins):
        for term in truth_by_protein[protein_id]:
            idx = term_to_idx.get(term)
            if idx is not None:
                truth[i, idx] = True

        predicted_scores = predictions_by_protein.get(protein_id, {})
        for term, score in predicted_scores.items():
            idx = term_to_idx.get(term)
            if idx is not None and score > pred[i, idx]:
                pred[i, idx] = float(score)
                explicit_pred_count += 1

    return {
        "proteins": proteins,
        "label_space": label_space,
        "truth": truth,
        "pred": pred,
        "ic": ic,
        "explicit_pred_count": explicit_pred_count,
    }


def grouped_pr_area(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = labels.astype(np.int8, copy=False)
    total_pos = int(labels.sum())
    total_neg = labels.size - total_pos
    if total_pos == 0 or total_neg == 0:
        return 0.0

    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    y = labels[order]

    tp_cum = np.cumsum(y, dtype=np.int64)
    fp_cum = np.cumsum(1 - y, dtype=np.int64)

    group_end = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s) - 1]
    tp = tp_cum[group_end].astype(np.float64)
    fp = fp_cum[group_end].astype(np.float64)

    recall = tp / total_pos
    precision = tp / (tp + fp)

    prev_recall = np.r_[0.0, recall[:-1]]
    prev_precision = np.r_[1.0, precision[:-1]]

    area = np.sum((recall - prev_recall) * (precision + prev_precision) / 2.0)
    return float(area)


def grouped_roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = labels.astype(np.int8, copy=False)
    total_pos = int(labels.sum())
    total_neg = labels.size - total_pos
    if total_pos == 0 or total_neg == 0:
        return 0.0

    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    y = labels[order]

    tp_cum = np.cumsum(y, dtype=np.int64)
    fp_cum = np.cumsum(1 - y, dtype=np.int64)

    group_end = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s) - 1]
    tp = tp_cum[group_end].astype(np.float64)
    fp = fp_cum[group_end].astype(np.float64)

    tpr = tp / total_pos
    fpr = fp / total_neg

    prev_tpr = np.r_[0.0, tpr[:-1]]
    prev_fpr = np.r_[0.0, fpr[:-1]]

    auc = np.sum((fpr - prev_fpr) * (tpr + prev_tpr) / 2.0)
    return float(auc)


def quantize_scores_0_100(pred: np.ndarray) -> np.ndarray:
    # 与原逻辑对齐：threshold 为 step/100，score >= threshold 算命中
    # 用 floor(score*100) 存到 0..100
    clipped = np.clip(pred, 0.0, 1.0)
    return np.floor(clipped * 100.0 + 1e-12).astype(np.int16)


def compute_threshold_profiles(truth: np.ndarray, pred: np.ndarray, ic: np.ndarray):
    n_proteins, _ = truth.shape

    truth_count = truth.sum(axis=1).astype(np.int32)
    truth_ic_sum = (truth * ic[None, :]).sum(axis=1, dtype=np.float64)

    q = quantize_scores_0_100(pred)

    pred_count_ge = np.zeros((n_proteins, 101), dtype=np.int32)
    inter_count_ge = np.zeros((n_proteins, 101), dtype=np.int32)
    pred_ic_ge = np.zeros((n_proteins, 101), dtype=np.float64)
    inter_ic_ge = np.zeros((n_proteins, 101), dtype=np.float64)

    ic_row = ic.astype(np.float64, copy=False)

    for i in range(n_proteins):
        steps = q[i]

        bucket_pred_count = np.bincount(steps, minlength=101)
        bucket_pred_ic = np.bincount(steps, weights=ic_row, minlength=101)

        true_mask = truth[i]
        true_steps = steps[true_mask]
        true_ic = ic_row[true_mask]

        if true_steps.size > 0:
            bucket_inter_count = np.bincount(true_steps, minlength=101)
            bucket_inter_ic = np.bincount(true_steps, weights=true_ic, minlength=101)
        else:
            bucket_inter_count = np.zeros(101, dtype=np.int64)
            bucket_inter_ic = np.zeros(101, dtype=np.float64)

        pred_count_ge[i] = np.cumsum(bucket_pred_count[::-1])[::-1]
        inter_count_ge[i] = np.cumsum(bucket_inter_count[::-1])[::-1]
        pred_ic_ge[i] = np.cumsum(bucket_pred_ic[::-1])[::-1]
        inter_ic_ge[i] = np.cumsum(bucket_inter_ic[::-1])[::-1]

    return {
        "truth_count": truth_count,
        "truth_ic_sum": truth_ic_sum,
        "pred_count_ge": pred_count_ge,
        "inter_count_ge": inter_count_ge,
        "pred_ic_ge": pred_ic_ge,
        "inter_ic_ge": inter_ic_ge,
    }


def compute_protein_centric_pr_from_profiles(profiles):
    truth_count = profiles["truth_count"]
    pred_count_ge = profiles["pred_count_ge"]
    inter_count_ge = profiles["inter_count_ge"]

    n_proteins = truth_count.shape[0]
    if n_proteins == 0:
        zeros = np.zeros(101, dtype=np.float64)
        return zeros, zeros

    pred_any = pred_count_ge > 0
    precision_per_protein = np.divide(
        inter_count_ge,
        pred_count_ge,
        out=np.zeros_like(inter_count_ge, dtype=np.float64),
        where=pred_count_ge > 0,
    )
    predicted_protein_count = pred_any.sum(axis=0)
    precision_sum = precision_per_protein.sum(axis=0)

    avg_precision = np.divide(
        precision_sum,
        predicted_protein_count,
        out=np.zeros(101, dtype=np.float64),
        where=predicted_protein_count > 0,
    )

    recall_per_protein = inter_count_ge / truth_count[:, None]
    avg_recall = recall_per_protein.mean(axis=0)
    return avg_precision, avg_recall


def compute_protein_centric_aupr(avg_precision: np.ndarray, avg_recall: np.ndarray) -> float:
    recall = avg_recall[::-1]
    precision = avg_precision[::-1]

    recall = np.r_[0.0, recall]
    precision = np.r_[1.0, precision]
    recall = np.maximum.accumulate(recall)

    area = np.sum((recall[1:] - recall[:-1]) * (precision[1:] + precision[:-1]) / 2.0)
    return float(area)


def compute_fmax_smin_from_profiles(profiles):
    truth_count = profiles["truth_count"]
    truth_ic_sum = profiles["truth_ic_sum"]
    pred_ic_ge = profiles["pred_ic_ge"]
    inter_ic_ge = profiles["inter_ic_ge"]

    n_proteins = truth_count.shape[0]
    if n_proteins == 0:
        return 0.0, 0.0, 0.0, 0.0

    avg_precision, avg_recall = compute_protein_centric_pr_from_profiles(profiles)

    denom = avg_precision + avg_recall
    f_all = np.divide(
        2.0 * avg_precision * avg_recall,
        denom,
        out=np.zeros_like(denom, dtype=np.float64),
        where=denom > 0,
    )

    ru = truth_ic_sum[:, None] - inter_ic_ge
    mi = pred_ic_ge - inter_ic_ge
    s_all = np.sqrt((ru.mean(axis=0) ** 2) + (mi.mean(axis=0) ** 2))

    best_f = float(np.max(f_all))
    f_candidates = np.flatnonzero(np.isclose(f_all, best_f))
    best_f_step = int(f_candidates[-1]) if f_candidates.size else 0

    best_s = float(np.min(s_all))
    s_candidates = np.flatnonzero(np.isclose(s_all, best_s))
    best_s_step = int(s_candidates[-1]) if s_candidates.size else 0

    return best_f, best_s, best_f_step / 100.0, best_s_step / 100.0


def compute_class_centric_mean_auc(truth: np.ndarray, pred: np.ndarray) -> float:
    n_proteins, n_terms = truth.shape
    if n_proteins == 0 or n_terms == 0:
        return 0.0

    positive_counts = truth.sum(axis=0).astype(np.int32)
    valid_mask = (positive_counts > 0) & (positive_counts < n_proteins)
    valid_indices = np.flatnonzero(valid_mask)
    if valid_indices.size == 0:
        return 0.0

    auc_values = []
    for idx in valid_indices:
        auc_values.append(
            grouped_roc_auc(
                truth[:, idx].astype(np.int8, copy=False),
                pred[:, idx],
            )
        )
    return float(np.mean(auc_values))


def compute_metrics_from_matrices(matrix_pack):
    truth = matrix_pack["truth"]
    pred = matrix_pack["pred"]
    ic = matrix_pack["ic"]

    n_proteins, n_terms = truth.shape
    if n_proteins == 0 or n_terms == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    profiles = compute_threshold_profiles(truth, pred, ic)
    avg_precision, avg_recall = compute_protein_centric_pr_from_profiles(profiles)
    fmax, smin, f_thr, s_thr = compute_fmax_smin_from_profiles(profiles)
    aupr = compute_protein_centric_aupr(avg_precision, avg_recall)
    auc = compute_class_centric_mean_auc(truth, pred)

    return fmax, smin, aupr, auc, f_thr, s_thr


def make_bootstrap_indices(n_proteins: int, n_boot: int, seed: int) -> np.ndarray:
    rng = random.Random(seed)
    out = np.empty((n_boot, n_proteins), dtype=np.int32)
    for b in range(n_boot):
        out[b] = [rng.randrange(n_proteins) for _ in range(n_proteins)]
    return out


def bootstrap_ci_from_matrices(matrix_pack, n_boot=1000, seed=42, alpha=0.05):
    truth = matrix_pack["truth"]
    pred = matrix_pack["pred"]
    ic = matrix_pack["ic"]

    n_proteins = truth.shape[0]
    if n_proteins == 0:
        return {}

    metric_names = ["Fmax", "Smin", "AUPR", "AUC"]
    boot_results = {m: [] for m in metric_names}

    boot_indices = make_bootstrap_indices(n_proteins, n_boot, seed)

    for idx in boot_indices:
        sampled_pack = {
            "truth": truth[idx],
            "pred": pred[idx],
            "ic": ic,
        }
        fmax, smin, aupr, auc, _, _ = compute_metrics_from_matrices(sampled_pack)
        boot_results["Fmax"].append(fmax)
        boot_results["Smin"].append(smin)
        boot_results["AUPR"].append(aupr)
        boot_results["AUC"].append(auc)

    ci = {}
    lo_idx = int(n_boot * alpha / 2)
    hi_idx = int(n_boot * (1 - alpha / 2))
    hi_idx = min(hi_idx, n_boot - 1)

    for m in metric_names:
        vals = np.sort(np.asarray(boot_results[m], dtype=np.float64))
        ci[f"{m}_lo"] = float(vals[lo_idx])
        ci[f"{m}_hi"] = float(vals[hi_idx])

    return ci


def discover_prediction_files(baselines_dir: Path, prediction_name: str) -> list[tuple[str, Path]]:
    discovered = []
    for baseline_dir in sorted(path for path in baselines_dir.iterdir() if path.is_dir()):
        for tsv_file in sorted(baseline_dir.glob(prediction_name)):
            discovered.append((tsv_file.stem, tsv_file))
    return discovered

RESULT_FIELDS_BASE = ["Method", "Fmax", "Smin", "AUPR", "AUC"]
RESULT_FIELDS_CI = [
    "Fmax_lo", "Fmax_hi", "Smin_lo", "Smin_hi",
    "AUPR_lo", "AUPR_hi", "AUC_lo", "AUC_hi",
]


def write_result_csv(path: Path, rows: list[dict[str, object]]) -> None:
    has_ci = any(k for row in rows for k in row if k.endswith("_lo"))
    fieldnames = RESULT_FIELDS_BASE + (RESULT_FIELDS_CI if has_ci else [])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_threshold_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Aspect", "Method", "FmaxThreshold", "SminThreshold"])
        writer.writeheader()
        writer.writerows(rows)


def warn_if_sparse_predictions(method_name: str, aspect: str, matrix_pack) -> None:
    protein_count = len(matrix_pack["proteins"])
    term_count = len(matrix_pack["label_space"])
    total_slots = protein_count * term_count
    if total_slots == 0:
        return
    density = matrix_pack["explicit_pred_count"] / total_slots
    if density < 0.05:
        print(
            f"WARNING: {method_name} aspect={ASPECT_NAMES[aspect]} looks sparse "
            f"(explicit score density={density:.4%}). "
            "If this file was thresholded before evaluation, Fmax/AUPR/AUC will be underestimated."
        )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default="fasta/test.tsv", help="truth tsv，格式为 pid go_term aspect")
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--restricted", action="store_true", help="在 restricted label space 上评估")
    mode_group.add_argument("--open", action="store_true", help="在 open label space 上评估")
    parser.add_argument("--obo", default="fasta/go-basic.obo", help="go-basic.obo 路径")
    parser.add_argument("--baselines-dir", default="baselines")
    parser.add_argument("--prediction-name", default="*.tsv", help="预测文件 glob，默认读取每个 baseline 目录下所有 .tsv")
    parser.add_argument("--out-dir", default="comparison/results")
    parser.add_argument(
        "--ic-pkl",
        default=None,
        help="可选的 IA/IC pickle 文件路径；若提供，则 Smin 使用该文件中的 term 信息量",
    )
    parser.add_argument("--min-positives", type=int, default=10, help="保留该参数以兼容旧命令；标准评估中不再使用")
    parser.add_argument("--n-boot", type=int, default=0, help="bootstrap 重采样次数，0 则跳过 CI 计算")
    return parser.parse_args()


def main():
    args = parse_args()

    truth_path = Path(args.truth)
    obo_path = Path(args.obo)
    baselines_dir = Path(args.baselines_dir)
    out_dir = Path(args.out_dir)
    ic_pkl_path = Path(args.ic_pkl) if args.ic_pkl else None

    if not truth_path.exists():
        raise FileNotFoundError(f"找不到 truth 文件: {truth_path}")
    if not obo_path.exists():
        raise FileNotFoundError(f"找不到 obo 文件: {obo_path}")
    if not baselines_dir.exists():
        raise FileNotFoundError(f"找不到 baselines 目录: {baselines_dir}")
    if ic_pkl_path is not None and not ic_pkl_path.exists():
        raise FileNotFoundError(f"找不到 IC pickle 文件: {ic_pkl_path}")

    parents, obo_go_to_aspect = load_obo(obo_path)

    truth_by_aspect, truth_go_to_aspect = load_truth(truth_path)
    go_to_aspect = dict(obo_go_to_aspect)
    go_to_aspect.update(truth_go_to_aspect)

    get_ancestors = build_ancestor_getter(parents, go_to_aspect)
    propagate_truth(truth_by_aspect, get_ancestors)

    if ic_pkl_path is not None:
        ic_by_aspect = load_ic_by_aspect_from_pickle(ic_pkl_path)
        print(f"using IC pickle for Smin: {ic_pkl_path}")
    else:
        ic_by_aspect = compute_ic_by_aspect(truth_by_aspect, go_to_aspect, parents)
        print("WARNING: IC estimated from test set truth — not suitable for publication. "
              "Use --ic-pkl with IC computed from training annotations instead.")
        print("using IC estimated from propagated truth")

    eval_terms_by_aspect = None
    if args.restricted:
        eval_terms_by_aspect = load_eval_label_space(truth_path, get_ancestors, go_to_aspect)

    prediction_files = discover_prediction_files(baselines_dir, args.prediction_name)
    if not prediction_files:
        raise FileNotFoundError(f"没有在 {baselines_dir} 下发现任何预测文件")

    rows_by_aspect = {aspect: [] for aspect in ASPECTS}
    threshold_rows = []

    for method_name, prediction_path in prediction_files:
        predictions_by_aspect = load_raw_predictions(prediction_path, go_to_aspect)
        propagate_predictions(predictions_by_aspect, get_ancestors)

        if args.restricted:
            predictions_by_aspect = restrict_predictions_to_eval_space(
                predictions_by_aspect,
                eval_terms_by_aspect,
            )

        print(f"evaluating: {method_name} <- {prediction_path}")

        for aspect in ASPECTS:
            matrix_pack = build_matrices(
                truth_by_aspect[aspect],
                predictions_by_aspect[aspect],
                ic_by_aspect[aspect],
                aspect,
                eval_terms_by_aspect,
            )

            if matrix_pack is None:
                row = {
                    "Method": method_name,
                    "Fmax": "0.000",
                    "Smin": "0.000",
                    "AUPR": "0.000",
                    "AUC": "0.000",
                }
                rows_by_aspect[aspect].append(row)
                threshold_rows.append(
                    {
                        "Aspect": ASPECT_NAMES[aspect],
                        "Method": method_name,
                        "FmaxThreshold": "0.0000",
                        "SminThreshold": "0.0000",
                    }
                )
                continue

            warn_if_sparse_predictions(method_name, aspect, matrix_pack)
            fmax, smin, aupr, auc, fmax_threshold, smin_threshold = compute_metrics_from_matrices(matrix_pack)

            row = {
                "Method": method_name,
                "Fmax": f"{fmax:.3f}",
                "Smin": f"{smin:.3f}",
                "AUPR": f"{aupr:.3f}",
                "AUC": f"{auc:.3f}",
            }

            if args.n_boot > 0:
                ci = bootstrap_ci_from_matrices(
                    matrix_pack,
                    n_boot=args.n_boot,
                )
                for key, val in ci.items():
                    row[key] = f"{val:.3f}"

            rows_by_aspect[aspect].append(row)
            threshold_rows.append(
                {
                    "Aspect": ASPECT_NAMES[aspect],
                    "Method": method_name,
                    "FmaxThreshold": f"{fmax_threshold:.4f}",
                    "SminThreshold": f"{smin_threshold:.4f}",
                }
            )

    for aspect in ASPECTS:
        output_path = out_dir / f"{ASPECT_NAMES[aspect]}.csv"
        write_result_csv(output_path, rows_by_aspect[aspect])
        print(f"wrote: {output_path}")

    threshold_path = out_dir / "th.csv"
    write_threshold_csv(threshold_path, threshold_rows)
    print(f"wrote: {threshold_path}")


if __name__ == "__main__":
    main()
