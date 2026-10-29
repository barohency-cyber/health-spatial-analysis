#!/usr/bin/env python3
"""GWR 带宽敏感性分析：换不同 k，局部系数方向是否稳定？

输出：
    output/gwr_sensitivity.csv
    output/gwr_sensitivity_summary.json
    output/gwr_sensitivity.png / .svg

用法：
    python code/gwr_sensitivity.py
    python code/gwr_sensitivity.py --min-k 20 --max-k 60 --step 5
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACCESS = ROOT / "output" / "accessibility.csv"
DEFAULT_OUTCOMES = ROOT / "data" / "outcomes.csv"
DEFAULT_OUTPUT = ROOT / "output"
DEFAULT_PREFIX = "gwr_sensitivity"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gwr  # noqa: E402
from spatial_regression import build_data, read_access, read_outcomes  # noqa: E402


def rank_correlation(left: np.ndarray, right: np.ndarray) -> float:
    left_rank = np.argsort(np.argsort(left))
    right_rank = np.argsort(np.argsort(right))
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def analyze(
    observations,
    y: np.ndarray,
    x: np.ndarray,
    min_k: int,
    max_k: int,
    step: int,
) -> tuple[list[dict[str, float]], dict[str, object]]:
    distances = gwr.distance_matrix(observations)
    candidates = list(range(min_k, max_k + 1, step))
    if candidates[-1] != max_k:
        candidates.append(max_k)

    table: list[dict[str, float]] = []
    coefficients: dict[int, np.ndarray] = {}
    for k in candidates:
        rows, score, _ = gwr.fit_gwr_at_k(observations, y, x, k, distances=distances)
        values = np.array([float(row["access_coef"]) for row in rows])
        coefficients[k] = values
        table.append(
            {
                "k": float(k),
                "aicc": float(score["aicc"]),
                "access_coef_min": float(values.min()),
                "access_coef_median": float(np.median(values)),
                "access_coef_max": float(values.max()),
                "negative_share": float((values < 0).mean()),
                "local_r2_median": float(
                    np.median([float(row["local_r2"]) for row in rows])
                ),
                "condition_number_max": float(
                    max(float(row["condition_number"]) for row in rows)
                ),
            }
        )

    selected_k = int(min(table, key=lambda row: row["aicc"])["k"])
    baseline = coefficients[selected_k]
    stability = []
    for row in table:
        k = int(row["k"])
        values = coefficients[k]
        stability.append(
            {
                "k": k,
                "pearson_with_selected": float(np.corrcoef(baseline, values)[0, 1]),
                "rank_corr_with_selected": rank_correlation(baseline, values),
                "same_sign_share": float((np.sign(baseline) == np.sign(values)).mean()),
            }
        )
    by_k = {item["k"]: item for item in stability}
    for row in table:
        extra = by_k[int(row["k"])]
        row["pearson_with_selected"] = extra["pearson_with_selected"]
        row["rank_corr_with_selected"] = extra["rank_corr_with_selected"]
        row["same_sign_share"] = extra["same_sign_share"]

    summary = {
        "method": "GWR bandwidth sensitivity",
        "n": len(y),
        "selected_k": selected_k,
        "candidate_ks": candidates,
        "table": table,
        "stability": {
            "min_pearson_with_selected": float(
                min(item["pearson_with_selected"] for item in stability)
            ),
            "min_rank_corr_with_selected": float(
                min(item["rank_corr_with_selected"] for item in stability)
            ),
            "min_same_sign_share": float(
                min(item["same_sign_share"] for item in stability)
            ),
        },
        "note": (
            "如果 k 变化时局部系数方向或排序大幅变化，GWR 结果只能作为探索性证据；"
            "不要只报告 AICc 最优的单一带宽。"
        ),
    }
    return table, summary


def write_outputs(
    table: Sequence[dict[str, float]],
    summary: dict[str, object],
    output_dir: Path,
    prefix: str,
    no_plot: bool,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{prefix}.csv"
    summary_path = output_dir / f"{prefix}_summary.json"
    png_path = output_dir / f"{prefix}.png"
    svg_path = output_dir / f"{prefix}.svg"
    fields = [
        "k", "aicc", "access_coef_min", "access_coef_median", "access_coef_max",
        "negative_share", "local_r2_median", "condition_number_max",
        "pearson_with_selected", "rank_corr_with_selected", "same_sign_share",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(table)
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)
    paths = {"csv": csv_path, "summary": summary_path}
    if not no_plot:
        plot_sensitivity(table, summary, png_path, svg_path)
        paths["plot_png"] = png_path
        paths["plot_svg"] = svg_path
    return paths


def plot_sensitivity(
    table: Sequence[dict[str, float]],
    summary: dict[str, object],
    png_path: Path,
    svg_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("提示：未安装 matplotlib，已跳过敏感性图。", file=sys.stderr)
        return

    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS", "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False

    ks = [int(row["k"]) for row in table]
    aicc = [row["aicc"] for row in table]
    low = [row["access_coef_min"] for row in table]
    median = [row["access_coef_median"] for row in table]
    high = [row["access_coef_max"] for row in table]
    selected = int(summary["selected_k"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=160)
    fig.patch.set_facecolor("#f7f9fb")
    for ax in axes:
        ax.set_facecolor("#fbfcfd")
        ax.grid(color="#c9d3dd", alpha=0.35, linewidth=0.6, linestyle="--")
        for spine in ax.spines.values():
            spine.set_color("#c9d3dd")

    axes[0].plot(ks, aicc, marker="o", color="#2c7fb8", linewidth=1.8, markersize=4)
    axes[0].axvline(selected, color="#d73027", linestyle="--", linewidth=1, label=f"AICc 最优 k={selected}")
    axes[0].set_xlabel("自适应带宽 k（最近邻数）", fontsize=9)
    axes[0].set_ylabel("AICc", fontsize=9)
    axes[0].set_title("带宽选择", fontsize=11, loc="left")
    axes[0].legend(fontsize=8, frameon=False)

    axes[1].fill_between(ks, low, high, color="#9ecae1", alpha=0.35, label="局部系数范围")
    axes[1].plot(ks, median, marker="o", color="#08519c", linewidth=1.8, markersize=4, label="中位数")
    axes[1].axhline(0, color="#555", linewidth=0.8)
    axes[1].set_xlabel("自适应带宽 k（最近邻数）", fontsize=9)
    axes[1].set_ylabel("局部可达性系数", fontsize=9)
    axes[1].set_title("局部系数随带宽变化", fontsize=11, loc="left")
    axes[1].legend(fontsize=8, frameon=False)

    fig.suptitle("GWR 带宽敏感性：局部系数是否稳定？", fontsize=15, fontweight="bold", x=0.06, ha="left")
    fig.text(
        0.06, 0.01,
        "模拟数据；若系数范围或排序随 k 大幅变化，GWR 只能作为探索性证据。",
        fontsize=8, color="#777",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GWR 带宽敏感性分析",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--access", type=Path, default=DEFAULT_ACCESS)
    parser.add_argument("--outcomes", type=Path, default=DEFAULT_OUTCOMES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--min-k", type=int, default=20)
    parser.add_argument("--max-k", type=int, default=60)
    parser.add_argument("--step", type=int, default=5)
    parser.add_argument("--no-plot", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        observations = read_access(args.access)
        outcomes = read_outcomes(args.outcomes)
        y, x, _, matched = build_data(observations, outcomes)
        table, summary = analyze(matched, y, x, args.min_k, args.max_k, args.step)
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2

    paths = write_outputs(table, summary, args.output_dir, args.prefix, args.no_plot)
    print(f"候选 k：{summary['candidate_ks']}")
    print(f"AICc 最优 k：{summary['selected_k']}")
    print("k    AICc     系数范围             与最优 k 的相关系数")
    for row in table:
        print(
            f"{int(row['k']):>3}  {row['aicc']:>8.2f}  "
            f"{row['access_coef_min']:.4f}–{row['access_coef_max']:.4f}  "
            f"r={row['pearson_with_selected']:.3f}"
        )
    print("输出文件：")
    for name, path in paths.items():
        print(f"  {name:10s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
