#!/usr/bin/env python3
"""地理加权回归（GWR）教学示例：看“可达性效应”是否随空间变化。

模型：
    y_i = β0(u_i, v_i) + β1(u_i, v_i) access_i
          + β2(u_i, v_i) elderly_i + β3(u_i, v_i) income_i + ε_i

每个社区单独做一次加权最小二乘，权重用自适应 bisquare 核：
    w_ij = (1 - (d_ij / b_i)^2)^2,  d_ij < b_i
其中 b_i 是第 k 个最近邻的距离。k 用 AICc 在候选带宽中自动选择。

输入：
    output/accessibility.csv
    data/outcomes.csv

输出：
    output/gwr_coefficients.csv
    output/gwr_summary.json
    output/gwr_map.png / .svg

注意：GWR 是探索性方法，局部系数不是因果效应；局部共线性、带宽选择和
多重比较都要在论文里讨论。
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
DEFAULT_PREFIX = "gwr"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spatial_autocorr import Observation, haversine_km  # noqa: E402
from spatial_regression import build_data, read_access, read_outcomes  # noqa: E402


def distance_matrix(observations: Sequence[Observation]) -> np.ndarray:
    n = len(observations)
    distances = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            value = haversine_km(
                observations[i].lon, observations[i].lat,
                observations[j].lon, observations[j].lat,
            )
            distances[i, j] = value
            distances[j, i] = value
    return distances


def adaptive_weights(distance_row: np.ndarray, k: int) -> np.ndarray:
    n = len(distance_row)
    k = min(max(k, 2), n - 1)
    order = np.argsort(distance_row, kind="stable")[:k]
    bandwidth = float(distance_row[order[-1]])
    weights = np.zeros(n, dtype=float)
    if bandwidth <= 1e-12:
        weights[order] = 1.0
        return weights
    scaled = distance_row[order] / bandwidth
    weights[order] = (1 - scaled ** 2) ** 2
    weights[order[-1]] = 0.0
    return weights


def local_fit(
    y: np.ndarray, x: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    xtwx = x.T @ (x * weights[:, None])
    xtwy = x.T @ (weights * y)
    try:
        beta = np.linalg.solve(xtwx, xtwy)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(xtwx) @ xtwy
    residual = y - x @ beta
    rss = float((weights * residual * residual).sum())
    return beta, xtwx, residual, rss


def local_r2(y: np.ndarray, weights: np.ndarray, residual: np.ndarray) -> float:
    weight_sum = float(weights.sum())
    if weight_sum <= 0:
        return float("nan")
    mean = float((weights * y).sum() / weight_sum)
    tss = float((weights * (y - mean) ** 2).sum())
    rss = float((weights * residual ** 2).sum())
    return 1 - rss / tss if tss > 0 else float("nan")


def gwr_aicc(
    y: np.ndarray,
    x: np.ndarray,
    distances: np.ndarray,
    k: int,
) -> dict[str, float]:
    n = len(y)
    fitted = np.zeros(n)
    hat_trace = 0.0
    hat_squared_trace = 0.0
    for i in range(n):
        weights = adaptive_weights(distances[i], k)
        beta, xtwx, _, _ = local_fit(y, x, weights)
        fitted[i] = float(x[i] @ beta)
        # GWR hat row: x_i' (X'WX)^-1 X'W
        row = x[i] @ np.linalg.pinv(xtwx) @ (x.T * weights)
        hat_trace += float(row[i])
        hat_squared_trace += float(row @ row)
    residual = y - fitted
    rss = float(residual @ residual)
    sigma = math.sqrt(rss / n)
    denominator = n - 2 - hat_trace
    if sigma <= 0 or denominator <= 0:
        aicc = math.inf
    else:
        aicc = (
            2 * n * math.log(sigma)
            + n * math.log(2 * math.pi)
            + n * (n + hat_trace) / denominator
        )
    return {
        "k": k,
        "aicc": aicc,
        "trace_hat": hat_trace,
        "trace_hat_squared": hat_squared_trace,
        "rss": rss,
        "sigma": sigma,
    }


def ols(y: np.ndarray, x: np.ndarray) -> dict[str, object]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    rss = float(residual @ residual)
    tss = float((y - y.mean()) @ (y - y.mean()))
    n, p = x.shape
    sigma2 = rss / n
    log_likelihood = -0.5 * n * (math.log(2 * math.pi) + math.log(sigma2) + 1)
    aic = -2 * log_likelihood + 2 * (p + 1)
    return {
        "beta": [float(value) for value in beta],
        "rss": rss,
        "r2": 1 - rss / tss,
        "aic": aic,
        "sigma2": sigma2,
    }


def fit_gwr(
    observations: Sequence[Observation],
    y: np.ndarray,
    x: np.ndarray,
    terms: Sequence[str],
    min_k: int,
    max_k: int,
    step: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    n = len(y)
    distances = distance_matrix(observations)
    if not 5 <= min_k < max_k <= n - 1:
        raise ValueError(f"带宽范围不合法：min_k={min_k}, max_k={max_k}, n={n}")
    candidates = list(range(min_k, max_k + 1, step))
    if candidates[-1] != max_k:
        candidates.append(max_k)
    scores = [gwr_aicc(y, x, distances, k) for k in candidates]
    best = min(scores, key=lambda item: item["aicc"])
    k = int(best["k"])

    # 用最优 k 重新计算每个社区的局部结果与标准误
    sigma2_global = float(best["rss"]) / (
        n - 2 * float(best["trace_hat"]) + float(best["trace_hat_squared"])
    )
    rows = []
    for i, observation in enumerate(observations):
        weights = adaptive_weights(distances[i], k)
        beta, xtwx, residual, _ = local_fit(y, x, weights)
        try:
            covariance = sigma2_global * np.linalg.pinv(xtwx)
            standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0))
        except np.linalg.LinAlgError:
            standard_errors = np.full(len(beta), float("nan"))
        t_values = beta / standard_errors
        rows.append(
            {
                "community_id": observation.community_id,
                "name": observation.name,
                "lon": observation.lon,
                "lat": observation.lat,
                "intercept": float(beta[0]),
                "access_coef": float(beta[1]),
                "elderly_coef": float(beta[2]),
                "income_coef": float(beta[3]),
                "access_t": float(t_values[1]),
                "local_r2": local_r2(y, weights, residual),
                "condition_number": float(np.linalg.cond(xtwx)),
                "bandwidth_k": k,
            }
        )

    access_values = [row["access_coef"] for row in rows]
    local_r2_values = [row["local_r2"] for row in rows]
    conditions = [row["condition_number"] for row in rows]
    summary = {
        "method": "Geographically Weighted Regression (adaptive bisquare)",
        "dependent_variable": "outpatient_utilization",
        "covariates": list(terms[1:]),
        "n": n,
        "bandwidth_k": k,
        "candidate_bandwidths": [
            {"k": int(item["k"]), "aicc": float(item["aicc"])} for item in scores
        ],
        "gwr_aicc": float(best["aicc"]),
        "trace_hat": float(best["trace_hat"]),
        "effective_parameters": float(best["trace_hat"]),
        "sigma2": sigma2_global,
        "global_ols": ols(y, x),
        "access_coef_summary": {
            "min": float(min(access_values)),
            "median": float(np.median(access_values)),
            "max": float(max(access_values)),
            "negative_share": float(sum(value < 0 for value in access_values) / n),
        },
        "local_r2_summary": {
            "min": float(min(local_r2_values)),
            "median": float(np.median(local_r2_values)),
            "max": float(max(local_r2_values)),
        },
        "condition_number_max": float(max(conditions)),
        "note": (
            "GWR 是探索性空间异质性方法，局部系数不是因果效应；"
            "需报告带宽选择、局部共线性和多重比较问题。"
        ),
    }
    return rows, summary


def write_outputs(
    rows: Sequence[dict[str, object]],
    summary: dict[str, object],
    observations: Sequence[Observation],
    output_dir: Path,
    prefix: str,
    facilities_path: Path | None,
    no_map: bool,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{prefix}_coefficients.csv"
    summary_path = output_dir / f"{prefix}_summary.json"
    map_png = output_dir / f"{prefix}_map.png"
    map_svg = output_dir / f"{prefix}_map.svg"

    fields = [
        "community_id", "name", "lon", "lat", "intercept", "access_coef",
        "elderly_coef", "income_coef", "access_t", "local_r2",
        "condition_number", "bandwidth_k",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)

    paths = {"coefficients": csv_path, "summary": summary_path}
    if not no_map:
        plot_gwr(rows, summary, facilities_path, map_png, map_svg)
        paths["map_png"] = map_png
        paths["map_svg"] = map_svg
    return paths


def _read_facilities(path: Path | None) -> list[dict[str, float | str]]:
    if path is None or not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [
            {"lon": float(row["lon"]), "lat": float(row["lat"])}
            for row in csv.DictReader(handle)
        ]


def plot_gwr(
    rows: Sequence[dict[str, object]],
    summary: dict[str, object],
    facilities_path: Path | None,
    png_path: Path,
    svg_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize, TwoSlopeNorm
        from matplotlib.lines import Line2D
    except ImportError:
        print("提示：未安装 matplotlib，已跳过 GWR 地图。", file=sys.stderr)
        return

    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS", "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False

    lons = [float(row["lon"]) for row in rows]
    lats = [float(row["lat"]) for row in rows]
    values = [float(row["access_coef"]) for row in rows]
    facilities = _read_facilities(facilities_path)
    minimum, maximum = min(values), max(values)
    if minimum < 0 < maximum:
        norm = TwoSlopeNorm(vmin=minimum, vcenter=0, vmax=maximum)
        cmap = "RdBu_r"
    else:
        norm = Normalize(vmin=minimum, vmax=maximum)
        cmap = "YlOrRd" if maximum > 0 else "Blues"

    fig, ax = plt.subplots(figsize=(11.2, 8.4), dpi=160)
    fig.patch.set_facecolor("#f7f9fb")
    ax.set_facecolor("#fbfcfd")
    scatter = ax.scatter(
        lons, lats, c=values, cmap=cmap, norm=norm,
        s=60, edgecolors="white", linewidths=0.6, zorder=3,
    )
    if facilities:
        ax.scatter(
            [item["lon"] for item in facilities],
            [item["lat"] for item in facilities],
            marker="^", s=95, c="#1a1a1a", edgecolors="white", linewidths=0.8, zorder=5,
        )

    span_x = max(lons) - min(lons)
    span_y = max(lats) - min(lats)
    ax.set_xlim(min(lons) - span_x * 0.06, max(lons) + span_x * 0.06)
    ax.set_ylim(min(lats) - span_y * 0.07, max(lats) + span_y * 0.07)
    ax.set_aspect(1 / math.cos(math.radians(float(np.mean(lats)))))
    ax.grid(color="#c9d3dd", alpha=0.35, linewidth=0.6, linestyle="--")
    ax.set_xlabel("经度", fontsize=9, color="#444")
    ax.set_ylabel("纬度", fontsize=9, color="#444")
    ax.tick_params(labelsize=8, colors="#555")
    for spine in ax.spines.values():
        spine.set_color("#c9d3dd")

    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    colorbar.set_label("局部可达性系数", fontsize=9)
    colorbar.ax.tick_params(labelsize=8)

    if facilities:
        ax.legend(
            handles=[Line2D([0], [0], marker="^", linestyle="none", markersize=8,
                            markerfacecolor="#1a1a1a", markeredgecolor="white", label="医疗机构")],
            loc="lower left", frameon=True, framealpha=0.95, fontsize=8.5,
        )

    global_beta = float(summary["global_ols"]["beta"][1])
    acc = summary["access_coef_summary"]
    fig.text(0.075, 0.955, "可达性效应的空间异质性（GWR）", fontsize=18, fontweight="bold", color="#17212b", ha="left")
    fig.text(
        0.075, 0.917,
        f"自适应 bisquare · k = {summary['bandwidth_k']} · "
        f"局部系数 {acc['min']:.4f} ~ {acc['max']:.4f} · 全局 OLS β = {global_beta:.4f}",
        fontsize=10, color="#5b6773", ha="left",
    )
    fig.text(
        0.075, 0.025,
        "模拟数据；GWR 局部系数是探索性结果，不是因果效应；负值表示该社区周边可达性—就医利用关系方向反转。",
        fontsize=8, color="#777", ha="left",
    )
    fig.subplots_adjust(left=0.075, right=0.88, top=0.875, bottom=0.09)
    fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="地理加权回归（GWR）：可达性效应的空间异质性",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--access", type=Path, default=DEFAULT_ACCESS, help="可达性 CSV")
    parser.add_argument("--outcomes", type=Path, default=DEFAULT_OUTCOMES, help="结果变量 CSV")
    parser.add_argument("--facilities", type=Path, default=ROOT / "data" / "facilities.csv", help="机构 CSV")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT, help="输出目录")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help="输出文件前缀")
    parser.add_argument("--min-k", type=int, default=20, help="自适应带宽候选最小值")
    parser.add_argument("--max-k", type=int, default=100, help="自适应带宽候选最大值")
    parser.add_argument("--step", type=int, default=5, help="候选带宽步长")
    parser.add_argument("--no-map", action="store_true", help="只输出表格，不画图")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        observations = read_access(args.access)
        outcomes = read_outcomes(args.outcomes)
        y, x, terms, matched = build_data(observations, outcomes)
        rows, summary = fit_gwr(
            matched, y, x, terms, args.min_k, args.max_k, args.step
        )
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2

    paths = write_outputs(
        rows, summary, matched, args.output_dir, args.prefix, args.facilities, args.no_map
    )
    acc = summary["access_coef_summary"]
    print(
        f"GWR 带宽 k={summary['bandwidth_k']}，AICc={float(summary['gwr_aicc']):.2f}；"
        f"局部可达性系数 {acc['min']:.4f} ~ {acc['max']:.4f}，"
        f"中位数 {acc['median']:.4f}。"
    )
    print("输出文件：")
    for name, path in paths.items():
        print(f"  {name:12s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
