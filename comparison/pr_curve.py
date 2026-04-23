from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

try:
    from comparison.compare import (
        ASPECT_NAMES,
        ASPECTS,
        compute_ic_by_aspect,
        discover_prediction_files,
        load_obo,
        load_raw_predictions,
        load_truth,
        load_ic_by_aspect_from_pickle,
        load_eval_label_space,
        propagate_predictions,
        propagate_truth,
        restrict_predictions_to_eval_space,
    )
except ModuleNotFoundError:
    from compare import (
        ASPECT_NAMES,
        ASPECTS,
        compute_ic_by_aspect,
        discover_prediction_files,
        load_obo,
        load_raw_predictions,
        load_truth,
        load_ic_by_aspect_from_pickle,
        load_eval_label_space,
        propagate_predictions,
        propagate_truth,
        restrict_predictions_to_eval_space,
    )


def collect_score_thresholds(predictions_by_protein: dict[str, dict[str, float]]) -> list[float]:
    scores = {score for terms_dict in predictions_by_protein.values() for score in terms_dict.values()}
    return sorted(scores)


def compute_curve(
    truth_by_protein: dict[str, set[str]],
    predictions_by_protein: dict[str, dict[str, float]],
    ic_by_term: dict[str, float],
):
    proteins = sorted(truth_by_protein)
    protein_count = len(proteins)
    if protein_count == 0:
        return [], 0.0, 0.0, 0.0, 0.0

    events = [
        (score, protein_id, term)
        for protein_id, terms_dict in predictions_by_protein.items()
        if protein_id in truth_by_protein
        for term, score in terms_dict.items()
    ]
    if not events:
        return [], 0.0, 0.0, 0.0, 0.0

    events.sort(key=lambda item: item[0], reverse=True)

    best_fmax = 0.0
    best_fmax_threshold = 0.0
    best_smin = math.inf
    best_smin_threshold = 0.0

    truth_sizes = {protein_id: len(truth_by_protein[protein_id]) for protein_id in proteins}
    truth_ic_sums = {
        protein_id: sum(ic_by_term.get(term, 0.0) for term in truth_by_protein[protein_id])
        for protein_id in proteins
    }
    predicted_term_count = {protein_id: 0 for protein_id in proteins}
    true_positive_count = {protein_id: 0 for protein_id in proteins}

    precision_sum = 0.0
    recall_sum = 0.0
    predicted_protein_count = 0
    ru_total = sum(truth_ic_sums.values())
    mi_total = 0.0

    rows = []

    index = 0
    while index < len(events):
        threshold = events[index][0]
        while index < len(events) and events[index][0] == threshold:
            _, protein_id, term = events[index]
            truth_terms = truth_by_protein[protein_id]
            truth_count = truth_sizes[protein_id]

            old_predicted_count = predicted_term_count[protein_id]
            old_true_positive_count = true_positive_count[protein_id]
            old_precision = old_true_positive_count / old_predicted_count if old_predicted_count else 0.0
            old_recall = old_true_positive_count / truth_count if truth_count else 0.0

            new_predicted_count = old_predicted_count + 1
            new_true_positive_count = old_true_positive_count + (1 if term in truth_terms else 0)
            predicted_term_count[protein_id] = new_predicted_count
            true_positive_count[protein_id] = new_true_positive_count

            new_precision = new_true_positive_count / new_predicted_count
            new_recall = new_true_positive_count / truth_count if truth_count else 0.0

            precision_sum += new_precision - old_precision
            recall_sum += new_recall - old_recall

            if old_predicted_count == 0:
                predicted_protein_count += 1

            term_ic = ic_by_term.get(term, 0.0)
            if term in truth_terms:
                ru_total -= term_ic
            else:
                mi_total += term_ic

            index += 1

        precision = precision_sum / predicted_protein_count if predicted_protein_count else 0.0
        recall = recall_sum / protein_count
        fmax = 0.0
        if precision + recall > 0.0:
            fmax = 2 * precision * recall / (precision + recall)
        smin = math.sqrt((ru_total / protein_count) ** 2 + (mi_total / protein_count) ** 2)
        rows.append(
            {
                "Threshold": f"{threshold:.6f}",
                "Precision": f"{precision:.6f}",
                "Recall": f"{recall:.6f}",
                "Fmax": f"{fmax:.6f}",
                "Smin": f"{smin:.6f}",
                "PredictedProteins": predicted_protein_count,
            }
        )

        if fmax > best_fmax or (fmax == best_fmax and threshold > best_fmax_threshold):
            best_fmax = fmax
            best_fmax_threshold = threshold
        if smin < best_smin or (smin == best_smin and threshold > best_smin_threshold):
            best_smin = smin
            best_smin_threshold = threshold

    return rows, best_fmax, best_smin, best_fmax_threshold, best_smin_threshold


def plot_pr_curves(path: Path, curve_rows: list[dict[str, object]], threshold_rows: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / ".mplconfig"))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    path.parent.mkdir(parents=True, exist_ok=True)

    methods = sorted({row["Method"] for row in curve_rows})
    method_colors = {method: plt.get_cmap("tab20")(index % 20) for index, method in enumerate(methods)}

    fig, axes = plt.subplots(1, len(ASPECTS), figsize=(18, 5), sharex=True, sharey=True)
    if len(ASPECTS) == 1:
        axes = [axes]

    for axis, aspect in zip(axes, ASPECTS):
        aspect_name = ASPECT_NAMES[aspect]
        aspect_curve_rows = [row for row in curve_rows if row["Aspect"] == aspect_name]
        for method in methods:
            method_rows = [row for row in aspect_curve_rows if row["Method"] == method]
            if not method_rows:
                continue
            recalls = [float(row["Recall"]) for row in method_rows]
            precisions = [float(row["Precision"]) for row in method_rows]
            fmax_values = [float(row["Fmax"]) for row in method_rows]
            best_index = max(range(len(method_rows)), key=lambda idx: fmax_values[idx])
            axis.plot(
                recalls,
                precisions,
                color=method_colors[method],
                linewidth=1.8,
                alpha=0.9,
            )
            axis.scatter(
                [recalls[best_index]],
                [precisions[best_index]],
                color=method_colors[method],
                s=18,
                zorder=3,
            )

        axis.set_title(aspect_name.upper())
        axis.set_xlabel("Recall")
        axis.grid(True, alpha=0.2)
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(0.0, 1.0)

    axes[0].set_ylabel("Precision")
    handles = [Line2D([0], [0], color=method_colors[method], lw=2, label=method) for method in methods]
    labels = methods
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle("Precision-Recall Curves", y=1.02)
    fig.tight_layout(rect=(0, 0.08, 1, 0.98))
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", default="fasta/open/test.tsv", help="truth tsv，格式为 pid go_term aspect")
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--restricted", action="store_true", help="在 restricted label space 上评估")
    mode_group.add_argument("--open", action="store_true", help="在 open label space 上评估")
    parser.add_argument("--obo", default="fasta/go-basic.obo", help="go-basic.obo 路径")
    parser.add_argument("--baselines-dir", default="baselines")
    parser.add_argument("--prediction-name", default="*.tsv", help="预测文件 glob，默认读取每个 baseline 目录下所有 .tsv")
    parser.add_argument("--ic-pkl", default=None, help="可选的 IC pickle 路径；不提供则从 truth 估计")
    parser.add_argument("--min-positives", type=int, default=0, help="保留至少多少正例的 term")
    parser.add_argument("--out-dir", default="comparison/pr_curve", help="输出目录")
    return parser.parse_args()


def normalize_output_dir(raw_path: str) -> Path:
    out_dir = Path(raw_path)
    if out_dir.is_absolute():
        out_dir = Path.cwd() / out_dir.as_posix().lstrip("/")
    return out_dir


def main():
    args = parse_args()
    truth_path = Path(args.truth)
    obo_path = Path(args.obo)
    baselines_dir = Path(args.baselines_dir)
    out_dir = normalize_output_dir(args.out_dir)

    if not truth_path.exists():
        raise FileNotFoundError(f"找不到 truth 文件: {truth_path}")
    if not obo_path.exists():
        raise FileNotFoundError(f"找不到 OBO 文件: {obo_path}")
    if not baselines_dir.exists():
        raise FileNotFoundError(f"找不到 baselines 目录: {baselines_dir}")

    truth_by_aspect, go_to_aspect = load_truth(truth_path)
    parents, go_to_aspect_from_obo = load_obo(obo_path)
    go_to_aspect.update(go_to_aspect_from_obo)
    propagate_truth(truth_by_aspect, parents, go_to_aspect)

    if args.ic_pkl is not None:
        ic_by_aspect = load_ic_by_aspect_from_pickle(Path(args.ic_pkl))
        print(f"using IC pickle for Smin: {args.ic_pkl}")
    else:
        ic_by_aspect = compute_ic_by_aspect(truth_by_aspect, go_to_aspect, parents)
        print("using IC estimated from propagated truth")

    eval_terms_by_aspect = None
    if args.restricted:
        eval_terms_by_aspect = load_eval_label_space(truth_path, parents, go_to_aspect)

    prediction_files = discover_prediction_files(baselines_dir, args.prediction_name)
    if not prediction_files:
        raise FileNotFoundError(f"没有在 {baselines_dir} 下发现任何预测文件")

    curve_rows = []
    threshold_rows = []

    for method_name, prediction_path in prediction_files:
        predictions_by_aspect = load_raw_predictions(prediction_path, go_to_aspect)
        propagate_predictions(predictions_by_aspect, parents, go_to_aspect)
        if args.restricted:
            predictions_by_aspect = restrict_predictions_to_eval_space(predictions_by_aspect, eval_terms_by_aspect)

        print(f"evaluating: {method_name} <- {prediction_path}")

        for aspect in ASPECTS:
            rows, best_fmax, best_smin, best_fmax_threshold, best_smin_threshold = compute_curve(
                truth_by_aspect[aspect],
                predictions_by_aspect[aspect],
                ic_by_aspect[aspect],
            )
            for row in rows:
                curve_rows.append(
                    {
                        "Aspect": ASPECT_NAMES[aspect],
                        "Method": method_name,
                        **row,
                    }
                )

            threshold_rows.append(
                {
                    "Aspect": ASPECT_NAMES[aspect],
                    "Method": method_name,
                    "FmaxThreshold": f"{best_fmax_threshold:.4f}",
                    "SminThreshold": f"{best_smin_threshold:.4f}",
                }
            )
            print(
                f"  {ASPECT_NAMES[aspect]}: best_fmax={best_fmax:.4f} @ {best_fmax_threshold:.4f}, "
                f"best_smin={best_smin:.4f} @ {best_smin_threshold:.4f}"
            )

    curve_path = out_dir / "pr_curve.png"
    plot_pr_curves(curve_path, curve_rows, threshold_rows)
    print(f"wrote: {curve_path}")


if __name__ == "__main__":
    main()
