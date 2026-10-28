#!/usr/bin/env python3
"""医疗资源可达性的空间自相关分析：全局 Moran's I + 局部 LISA。

输入：`output/accessibility.csv`（先运行 two_step_fca.py）
输出：
    output/lisa.csv            每个社区的局部 Moran's I、p 值、FDR q 值、聚类类型
    output/lisa.geojson        QGIS 可直接打开的 LISA 结果
    output/lisa_summary.json   全局 Moran's I、置换检验、聚类计数
    output/lisa_map.png/svg    LISA 聚类图

方法要点：
- 空间权重：k 个最近邻（默认 k=6），行标准化；
- 全局 Moran's I 使用 999 次随机置换做伪 p 值；
- 局部 Moran's I 使用条件置换，并用 Benjamini-Hochberg 控制 FDR；
- 显著且 z_i>0、空间滞后>0 为高-高聚集，余下类推。

注意：这里的“显著”是统计上的空间聚集，不是因果结论。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACCESS = ROOT / "output" / "accessibility.csv"
DEFAULT_FACILITIES = ROOT / "data" / "facilities.csv"
DEFAULT_OUTPUT = ROOT / "output"
DEFAULT_PREFIX = "lisa"
EARTH_RADIUS_KM = 6371.0088

CLUSTER_COLORS = {
    "High-High": "#d73027",
    "Low-Low": "#4575b4",
    "High-Low": "#fdae61",
    "Low-High": "#abd9e9",
    "Not Significant": "#bdbdbd",
}
CLUSTER_LABELS = {
    "High-High": "高-高聚集",
    "Low-Low": "低-低聚集",
    "High-Low": "高-低异常",
    "Low-High": "低-高异常",
    "Not Significant": "不显著",
}


@dataclass(frozen=True)
class Observation:
    community_id: str
    name: str
    lon: float
    lat: float
    population: float
    access_per_1000: float


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    )
    a = min(1.0, max(0.0, a))
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


def read_observations(path: Path) -> list[Observation]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
    except FileNotFoundError as exc:
        raise ValueError(f"找不到文件：{path}；请先运行 python code/two_step_fca.py") from exc

    required = {"community_id", "name", "lon", "lat", "population", "access_per_1000"}
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"{path} 缺少字段：{', '.join(missing)}")
    if not rows:
        raise ValueError(f"{path} 没有数据行")

    observations: list[Observation] = []
    seen: set[str] = set()
    for line_no, row in enumerate(rows, start=2):
        community_id = (row.get("community_id") or "").strip()
        if not community_id:
            raise ValueError(f"{path} 第 {line_no} 行：community_id 不能为空")
        if community_id in seen:
            raise ValueError(f"{path} 第 {line_no} 行：community_id 重复")
        seen.add(community_id)
        try:
            observation = Observation(
                community_id=community_id,
                name=(row.get("name") or community_id).strip(),
                lon=float(row["lon"]),
                lat=float(row["lat"]),
                population=float(row["population"]),
                access_per_1000=float(row["access_per_1000"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path} 第 {line_no} 行存在无法解析的数值") from exc
        observations.append(observation)

    values = [item.access_per_1000 for item in observations]
    if max(values) - min(values) <= 1e-12:
        raise ValueError("access_per_1000 没有变异，无法计算 Moran's I")
    return observations


def build_knn_weights(
    observations: Sequence[Observation], k: int
) -> tuple[list[list[float]], list[list[int]]]:
    n = len(observations)
    if n < 3:
        raise ValueError("至少需要 3 个需求点")
    if not 1 <= k < n:
        raise ValueError(f"k 必须在 1 到 {n - 1} 之间")

    weights = [[0.0] * n for _ in range(n)]
    neighbors: list[list[int]] = []
    for i, origin in enumerate(observations):
        distances = []
        for j, target in enumerate(observations):
            if i == j:
                continue
            distance = haversine_km(origin.lon, origin.lat, target.lon, target.lat)
            distances.append((distance, j))
        distances.sort(key=lambda item: (item[0], item[1]))
        nearest = [j for _, j in distances[:k]]
        neighbors.append(nearest)
        for j in nearest:
            weights[i][j] = 1.0 / k
    return weights, neighbors


def _moran_i(values: Sequence[float], weights: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    n = len(values)
    mean = statistics.fmean(values)
    deviations = [value - mean for value in values]
    s0 = sum(sum(row) for row in weights)
    denominator = sum(value * value for value in deviations)
    if denominator <= 0 or s0 <= 0:
        raise ValueError("无法计算 Moran's I：数值没有变异或权重矩阵为空")
    numerator = n * sum(
        weights[i][j] * deviations[i] * deviations[j]
        for i in range(n)
        for j in range(n)
    )
    observed = numerator / (s0 * denominator)
    expected = -1.0 / (n - 1)

    s1 = 0.5 * sum(
        (weights[i][j] + weights[j][i]) ** 2
        for i in range(n)
        for j in range(n)
    )
    s2 = sum(
        (sum(weights[i]) + sum(weights[j][i] for j in range(n))) ** 2
        for i in range(n)
    )
    variance = (
        (n * n * s1 - n * s2 + 3 * s0 * s0) / ((n * n - 1) * s0 * s0)
        - expected * expected
    )
    z_score = (observed - expected) / math.sqrt(variance) if variance > 0 else 0.0
    return observed, expected, z_score


def global_moran(
    values: Sequence[float],
    weights: Sequence[Sequence[float]],
    permutations: int,
    seed: int,
) -> dict[str, float]:
    observed, expected, z_score = _moran_i(values, weights)
    n = len(values)
    rng = random.Random(seed)
    extreme = 0
    threshold = abs(observed - expected) - 1e-12
    for _ in range(permutations):
        shuffled = list(values)
        rng.shuffle(shuffled)
        permuted, _, _ = _moran_i(shuffled, weights)
        if abs(permuted - expected) >= threshold:
            extreme += 1
    p_permutation = (extreme + 1) / (permutations + 1)
    p_normality = math.erfc(abs(z_score) / math.sqrt(2))
    return {
        "global_i": observed,
        "expected_i": expected,
        "z_score": z_score,
        "p_normality": p_normality,
        "p_permutation": p_permutation,
    }


def bh_fdr(p_values: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg FDR；返回与输入同序的 q 值。"""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    q_values = [1.0] * m
    previous = 1.0
    for rank in range(m, 0, -1):
        index = order[rank - 1]
        value = min(previous, p_values[index] * m / rank)
        q_values[index] = value
        previous = value
    return q_values


def classify_cluster(z_value: float, spatial_lag: float, significant: bool) -> str:
    if not significant:
        return "Not Significant"
    if z_value > 0 and spatial_lag > 0:
        return "High-High"
    if z_value < 0 and spatial_lag < 0:
        return "Low-Low"
    if z_value > 0 and spatial_lag < 0:
        return "High-Low"
    return "Low-High"


def local_moran(
    values: Sequence[float],
    weights: Sequence[Sequence[float]],
    neighbors: Sequence[Sequence[int]],
    permutations: int,
    seed: int,
    alpha: float,
) -> list[dict[str, float | str]]:
    n = len(values)
    mean = statistics.fmean(values)
    deviations = [value - mean for value in values]
    m2 = sum(value * value for value in deviations) / n
    if m2 <= 0:
        raise ValueError("无法计算局部 Moran's I：数值没有变异")

    rows: list[dict[str, float | str]] = []
    rng = random.Random(seed + 1)
    for i in range(n):
        neighbor_indices = list(neighbors[i])
        observed_lag = sum(weights[i][j] * deviations[j] for j in neighbor_indices)
        observed_local = deviations[i] * observed_lag / m2

        others = deviations[:i] + deviations[i + 1:]
        extreme = 0
        threshold = abs(observed_local) - 1e-12
        for _ in range(permutations):
            rng.shuffle(others)
            permuted_lag = 0.0
            for j in neighbor_indices:
                position = j if j < i else j - 1
                permuted_lag += weights[i][j] * others[position]
            permuted_local = deviations[i] * permuted_lag / m2
            if abs(permuted_local) >= threshold:
                extreme += 1
        p_sim = (extreme + 1) / (permutations + 1)
        rows.append(
            {
                "z": deviations[i],
                "spatial_lag": observed_lag,
                "local_i": observed_local,
                "p_sim": p_sim,
            }
        )

    q_values = bh_fdr([float(row["p_sim"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        significant = q_value <= alpha
        row["q_value"] = q_value
        row["cluster"] = classify_cluster(
            float(row["z"]), float(row["spatial_lag"]), significant
        )
    return rows


def write_outputs(
    observations: Sequence[Observation],
    rows: Sequence[dict[str, float | str]],
    summary: dict[str, object],
    output_dir: Path,
    prefix: str,
    facilities_path: Path | None,
    no_map: bool,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{prefix}.csv"
    geojson_path = output_dir / f"{prefix}.geojson"
    summary_path = output_dir / f"{prefix}_summary.json"
    map_png = output_dir / f"{prefix}_map.png"
    map_svg = output_dir / f"{prefix}_map.svg"

    fieldnames = [
        "community_id",
        "name",
        "lon",
        "lat",
        "population",
        "access_per_1000",
        "z",
        "spatial_lag",
        "local_i",
        "p_sim",
        "q_value",
        "cluster",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for observation, row in zip(observations, rows):
            writer.writerow(
                {
                    "community_id": observation.community_id,
                    "name": observation.name,
                    "lon": observation.lon,
                    "lat": observation.lat,
                    "population": observation.population,
                    "access_per_1000": f"{observation.access_per_1000:.6f}",
                    "z": f"{float(row['z']):.6f}",
                    "spatial_lag": f"{float(row['spatial_lag']):.6f}",
                    "local_i": f"{float(row['local_i']):.6f}",
                    "p_sim": f"{float(row['p_sim']):.6f}",
                    "q_value": f"{float(row['q_value']):.6f}",
                    "cluster": row["cluster"],
                }
            )

    features = []
    for observation, row in zip(observations, rows):
        features.append(
            {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [observation.lon, observation.lat],
                },
                "properties": {
                    "community_id": observation.community_id,
                    "name": observation.name,
                    "population": observation.population,
                    "access_per_1000": round(observation.access_per_1000, 6),
                    "local_i": round(float(row["local_i"]), 6),
                    "q_value": round(float(row["q_value"]), 6),
                    "cluster": row["cluster"],
                },
            }
        )
    with geojson_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {"type": "FeatureCollection", "features": features},
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)

    paths = {"csv": csv_path, "geojson": geojson_path, "summary": summary_path}
    if not no_map:
        plot_lisa(observations, rows, summary, facilities_path, map_png, map_svg)
        paths["map_png"] = map_png
        paths["map_svg"] = map_svg
    return paths


def _read_facilities(path: Path | None) -> list[dict[str, float | str]]:
    if path is None or not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [
            {
                "name": row["name"],
                "lon": float(row["lon"]),
                "lat": float(row["lat"]),
            }
            for row in csv.DictReader(handle)
        ]


def plot_lisa(
    observations: Sequence[Observation],
    rows: Sequence[dict[str, float | str]],
    summary: dict[str, object],
    facilities_path: Path | None,
    png_path: Path,
    svg_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        print("提示：未安装 matplotlib，已跳过 LISA 地图。", file=sys.stderr)
        return

    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False

    lons = [item.lon for item in observations]
    lats = [item.lat for item in observations]
    clusters = [str(row["cluster"]) for row in rows]
    colors = [CLUSTER_COLORS.get(cluster, CLUSTER_COLORS["Not Significant"]) for cluster in clusters]
    facilities = _read_facilities(facilities_path)

    fig, ax = plt.subplots(figsize=(11.2, 8.4), dpi=160)
    fig.patch.set_facecolor("#f7f9fb")
    ax.set_facecolor("#fbfcfd")
    ax.scatter(lons, lats, s=58, c=colors, edgecolors="white", linewidths=0.6, zorder=3)
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
    ax.set_aspect(1 / math.cos(math.radians(statistics.fmean(lats))))
    ax.grid(color="#c9d3dd", alpha=0.35, linewidth=0.6, linestyle="--")
    ax.set_xlabel("经度", fontsize=9, color="#444")
    ax.set_ylabel("纬度", fontsize=9, color="#444")
    ax.tick_params(labelsize=8, colors="#555")
    for spine in ax.spines.values():
        spine.set_color("#c9d3dd")

    handles = []
    for cluster in ("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"):
        count = int(summary["cluster_counts"].get(cluster, 0))
        handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=7,
                markerfacecolor=CLUSTER_COLORS[cluster], markeredgecolor="white",
                label=f"{CLUSTER_LABELS[cluster]}（{count}）",
            )
        )
    if facilities:
        handles.append(
            Line2D(
                [0], [0], marker="^", linestyle="none", markersize=8,
                markerfacecolor="#1a1a1a", markeredgecolor="white", label="医疗机构",
            )
        )
    legend = ax.legend(
        handles=handles, title="LISA 聚类", loc="center left",
        bbox_to_anchor=(1.015, 0.5), frameon=True, framealpha=0.96,
        fontsize=8.5, title_fontsize=9.5, borderpad=0.8, labelspacing=0.65,
    )
    legend.get_frame().set_edgecolor("#d0d7de")

    alpha = float(summary["alpha"])
    fig.text(0.075, 0.955, "医疗资源可达性的空间聚集（LISA）", fontsize=18, fontweight="bold", color="#17212b", ha="left")
    fig.text(
        0.075, 0.917,
        f"全局 Moran's I = {float(summary['global_i']):.3f} · "
        f"置换 p = {float(summary['p_permutation']):.3f} · "
        f"k = {summary['k']} 最近邻 · FDR q ≤ {alpha:g}",
        fontsize=10, color="#5b6773", ha="left",
    )
    fig.text(
        0.075, 0.025,
        "模拟数据；LISA 的“显著”仅表示空间聚集，不代表因果；正式研究需做权重与多重检验敏感性分析。",
        fontsize=8, color="#777", ha="left",
    )
    fig.subplots_adjust(left=0.075, right=0.79, top=0.875, bottom=0.09)
    fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def analyze(
    observations: Sequence[Observation],
    k: int,
    permutations: int,
    seed: int,
    alpha: float,
) -> tuple[list[dict[str, float | str]], dict[str, object]]:
    values = [item.access_per_1000 for item in observations]
    weights, neighbors = build_knn_weights(observations, k)
    global_stats = global_moran(values, weights, permutations, seed)
    local_rows = local_moran(values, weights, neighbors, permutations, seed, alpha)

    counts = {name: 0 for name in CLUSTER_LABELS}
    for row in local_rows:
        counts[str(row["cluster"])] = counts.get(str(row["cluster"]), 0) + 1
    summary: dict[str, object] = {
        "method": "Global and local Moran's I (LISA)",
        "variable": "access_per_1000",
        "n": len(observations),
        "k": k,
        "weights": {"type": "k nearest neighbors", "row_standardized": True},
        "permutations": permutations,
        "seed": seed,
        "alpha": alpha,
        "global_i": global_stats["global_i"],
        "expected_i": global_stats["expected_i"],
        "z_score": global_stats["z_score"],
        "p_normality": global_stats["p_normality"],
        "p_permutation": global_stats["p_permutation"],
        "cluster_counts": counts,
        "note": (
            "LISA 聚类基于局部 Moran's I 的条件置换 p 值，经 Benjamini-Hochberg FDR 校正；"
            "k 近邻权重用直线距离构造，正式研究应做不同 k 和距离阈值的敏感性分析。"
        ),
    }
    return local_rows, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="医疗资源可达性的全局 Moran's I 与局部 LISA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--accessibility", type=Path, default=DEFAULT_ACCESS, help="可达性 CSV")
    parser.add_argument("--facilities", type=Path, default=DEFAULT_FACILITIES, help="机构 CSV（仅用于地图）")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT, help="输出目录")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help="输出文件前缀")
    parser.add_argument("--k", type=int, default=6, help="k 近邻空间权重")
    parser.add_argument("--permutations", type=int, default=999, help="置换次数")
    parser.add_argument("--seed", type=int, default=20261028, help="随机种子，保证可复现")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR 显著性水平")
    parser.add_argument("--no-map", action="store_true", help="只输出表格，不画 LISA 图")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.permutations < 99:
        print("错误：permutations 至少为 99", file=sys.stderr)
        return 2
    if not 0 < args.alpha < 1:
        print("错误：alpha 必须在 (0, 1) 内", file=sys.stderr)
        return 2
    try:
        observations = read_observations(args.accessibility)
        rows, summary = analyze(
            observations, args.k, args.permutations, args.seed, args.alpha
        )
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2

    paths = write_outputs(
        observations, rows, summary, args.output_dir, args.prefix,
        args.facilities, args.no_map,
    )
    counts = summary["cluster_counts"]
    print(
        f"全局 Moran's I = {float(summary['global_i']):.3f}，"
        f"置换 p = {float(summary['p_permutation']):.3f}；"
        f"高-高 {counts['High-High']} 个，低-低 {counts['Low-Low']} 个，"
        f"高-低 {counts['High-Low']} 个，低-高 {counts['Low-High']} 个。"
    )
    print("输出文件：")
    for name, path in paths.items():
        print(f"  {name:10s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
