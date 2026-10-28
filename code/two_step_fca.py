#!/usr/bin/env python3
"""两步移动搜索法（2SFCA）极简实现：CSV 进，可达性表 + 地图出。

方法（两步）：
    第一步 对每个供给点 j：
        R_j = S_j / Σ_i [P_i * W_ij]
    第二步 对每个需求点 i：
        A_i = Σ_j [R_j * W_ij]
其中 S_j 为机构容量（如床位数），P_i 为需求点人口，W_ij 为距离衰减权重，
只保留距离 d_ij <= d0 的供需对。A_i 的单位是“容量/人”，乘 1000 后
解释为“每千人可获得的床位数”。

注意：这是直线距离版 2SFCA，适合课程作业、方法入门和复试展示；
正式论文若声称“路网可达性”，应把 haversine 距离替换为路网出行时间/距离。

运行：
    python code/two_step_fca.py
    python code/two_step_fca.py --threshold-km 10 --decay gaussian
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DEMAND = ROOT / "data" / "demand.csv"
DEFAULT_FACILITIES = ROOT / "data" / "facilities.csv"
DEFAULT_OUTPUT = ROOT / "output"
EARTH_RADIUS_KM = 6371.0088
ZERO_COLOR = "#d73027"
BLUE_STEPS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]


@dataclass(frozen=True)
class DemandPoint:
    community_id: str
    name: str
    lon: float
    lat: float
    population: float


@dataclass(frozen=True)
class Facility:
    facility_id: str
    name: str
    lon: float
    lat: float
    capacity: float


@dataclass(frozen=True)
class AccessRecord:
    demand: DemandPoint
    access_per_1000: float
    facilities_in_catchment: int
    nearest_facility_km: float


@dataclass(frozen=True)
class TwoStepResult:
    records: list[AccessRecord]
    facility_ratios: list[float]
    facilities_without_demand: int


# ---------------------------------------------------------------------------
# 读数据与校验：报错时直接指出文件、行号和字段，方便替换成自己的数据
# ---------------------------------------------------------------------------

def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
    except FileNotFoundError as exc:
        raise ValueError(f"找不到数据文件：{path}") from exc
    return fields, rows


def _require_columns(fields: Sequence[str], required: Sequence[str], path: Path) -> None:
    missing = [name for name in required if name not in fields]
    if missing:
        raise ValueError(
            f"{path} 缺少字段：{', '.join(missing)}；"
            f"当前字段为：{', '.join(fields) or '（空表头）'}"
        )


def _text(row: dict[str, str], key: str, path: Path, line_no: int) -> str:
    value = (row.get(key) or "").strip()
    if not value:
        raise ValueError(f"{path} 第 {line_no} 行：{key} 不能为空")
    return value


def _number(row: dict[str, str], key: str, path: Path, line_no: int) -> float:
    raw = (row.get(key) or "").strip()
    if not raw:
        raise ValueError(f"{path} 第 {line_no} 行：{key} 不能为空")
    try:
        return float(raw)
    except ValueError as exc:
        raise ValueError(f"{path} 第 {line_no} 行：{key} 不是数字（{raw!r}）") from exc


def load_demand(path: Path) -> list[DemandPoint]:
    fields, rows = _read_table(path)
    _require_columns(fields, ["community_id", "name", "lon", "lat", "population"], path)
    if not rows:
        raise ValueError(f"{path} 没有数据行")
    points: list[DemandPoint] = []
    seen: set[str] = set()
    for line_no, row in enumerate(rows, start=2):
        community_id = _text(row, "community_id", path, line_no)
        if community_id in seen:
            raise ValueError(f"{path} 第 {line_no} 行：community_id 重复（{community_id}）")
        seen.add(community_id)
        lon = _number(row, "lon", path, line_no)
        lat = _number(row, "lat", path, line_no)
        population = _number(row, "population", path, line_no)
        if not -180 <= lon <= 180 or not -90 <= lat <= 90:
            raise ValueError(f"{path} 第 {line_no} 行：经纬度超出合法范围")
        if population < 0:
            raise ValueError(f"{path} 第 {line_no} 行：population 不能为负数")
        points.append(
            DemandPoint(
                community_id=community_id,
                name=_text(row, "name", path, line_no),
                lon=lon,
                lat=lat,
                population=population,
            )
        )
    if sum(point.population for point in points) <= 0:
        raise ValueError(f"{path}：总人口必须大于 0")
    return points


def load_facilities(path: Path) -> list[Facility]:
    fields, rows = _read_table(path)
    _require_columns(fields, ["facility_id", "name", "lon", "lat", "capacity"], path)
    if not rows:
        raise ValueError(f"{path} 没有数据行")
    facilities: list[Facility] = []
    seen: set[str] = set()
    for line_no, row in enumerate(rows, start=2):
        facility_id = _text(row, "facility_id", path, line_no)
        if facility_id in seen:
            raise ValueError(f"{path} 第 {line_no} 行：facility_id 重复（{facility_id}）")
        seen.add(facility_id)
        lon = _number(row, "lon", path, line_no)
        lat = _number(row, "lat", path, line_no)
        capacity = _number(row, "capacity", path, line_no)
        if not -180 <= lon <= 180 or not -90 <= lat <= 90:
            raise ValueError(f"{path} 第 {line_no} 行：经纬度超出合法范围")
        if capacity <= 0:
            raise ValueError(
                f"{path} 第 {line_no} 行：capacity 必须大于 0；"
                "没有供给能力的机构请先从供给点表中剔除"
            )
        facilities.append(
            Facility(
                facility_id=facility_id,
                name=_text(row, "name", path, line_no),
                lon=lon,
                lat=lat,
                capacity=capacity,
            )
        )
    return facilities


# ---------------------------------------------------------------------------
# 2SFCA 核心
# ---------------------------------------------------------------------------

def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """球面大圆距离（km）；在城市尺度上足够接近平面距离。"""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    )
    a = min(1.0, max(0.0, a))
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


def distance_weight(distance_km: float, threshold_km: float, decay: str) -> float:
    """距离衰减权重；超过阈值的供需对不参与计算。"""
    if distance_km > threshold_km:
        return 0.0
    if decay == "binary":
        return 1.0
    if decay == "gaussian":
        # 截断高斯核：sigma = d0 / 2，权重在阈值处约为 0.135。
        # 这是教学示例的默认核，正式研究应对核函数与带宽做敏感性分析。
        sigma = threshold_km / 2
        return math.exp(-0.5 * (distance_km / sigma) ** 2)
    raise ValueError(f"未知衰减方式：{decay}")


def compute_2sfca(
    demand: Sequence[DemandPoint],
    facilities: Sequence[Facility],
    threshold_km: float,
    decay: str = "binary",
) -> TwoStepResult:
    if threshold_km <= 0:
        raise ValueError("threshold_km 必须大于 0")
    if not demand:
        raise ValueError("至少需要一个需求点")
    if not facilities:
        raise ValueError("至少需要一家医疗机构")

    distances = [
        [haversine_km(point.lon, point.lat, facility.lon, facility.lat) for facility in facilities]
        for point in demand
    ]
    weights = [
        [distance_weight(distance, threshold_km, decay) for distance in row]
        for row in distances
    ]

    ratios: list[float] = []
    facilities_without_demand = 0
    for j, facility in enumerate(facilities):
        potential_population = math.fsum(
            point.population * weights[i][j] for i, point in enumerate(demand)
        )
        if potential_population > 0:
            ratios.append(facility.capacity / potential_population)
        else:
            # 该机构在阈值内没有任何需求点：按 2SFCA 定义，比值为 0。
            ratios.append(0.0)
            facilities_without_demand += 1

    records: list[AccessRecord] = []
    for i, point in enumerate(demand):
        access = math.fsum(ratios[j] * weights[i][j] for j in range(len(facilities)))
        in_catchment = sum(1 for weight in weights[i] if weight > 0)
        records.append(
            AccessRecord(
                demand=point,
                access_per_1000=access * 1000,
                facilities_in_catchment=in_catchment,
                nearest_facility_km=min(distances[i]),
            )
        )
    return TwoStepResult(
        records=records,
        facility_ratios=ratios,
        facilities_without_demand=facilities_without_demand,
    )


# ---------------------------------------------------------------------------
# 统计与导出
# ---------------------------------------------------------------------------

def _quantile(values: Sequence[float], q: float) -> float:
    """线性插值分位数；用于图例分段。"""
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def weighted_quantile(values: Sequence[float], weights: Sequence[float], q: float) -> float:
    """按人口加权的分位数；人口为 0 的点会被自然忽略。"""
    pairs = sorted(zip(values, weights))
    total = math.fsum(weight for _, weight in pairs)
    if total <= 0:
        raise ValueError("加权分位数的权重总和必须大于 0")
    target = total * q
    cumulative = 0.0
    for value, weight in pairs:
        cumulative += weight
        if cumulative >= target:
            return value
    return pairs[-1][0]


def summarize(
    result: TwoStepResult,
    threshold_km: float,
    decay: str,
) -> dict[str, object]:
    records = result.records
    population = math.fsum(record.demand.population for record in records)
    values = [record.access_per_1000 for record in records]
    populations = [record.demand.population for record in records]
    covered_population = math.fsum(
        record.demand.population for record in records if record.facilities_in_catchment > 0
    )
    zero_access_population = math.fsum(
        record.demand.population for record in records if record.access_per_1000 <= 0
    )
    mean_access = math.fsum(
        value * weight for value, weight in zip(values, populations)
    ) / population
    p10 = weighted_quantile(values, populations, 0.10)
    p50 = weighted_quantile(values, populations, 0.50)
    p90 = weighted_quantile(values, populations, 0.90)
    return {
        "method": "2SFCA",
        "distance": "haversine straight-line",
        "decay": decay,
        "threshold_km": threshold_km,
        "demand_count": len(records),
        "facility_count": len(result.facility_ratios),
        "total_population": population,
        "covered_population": covered_population,
        "coverage_rate": covered_population / population,
        "zero_access_population": zero_access_population,
        "zero_access_rate": zero_access_population / population,
        "mean_access_per_1000": mean_access,
        "median_access_per_1000": p50,
        "p10_access_per_1000": p10,
        "p90_access_per_1000": p90,
        "p90_p10_ratio": (p90 / p10) if p10 > 0 else None,
        "min_access_per_1000": min(values),
        "max_access_per_1000": max(values),
        "facilities_without_demand": result.facilities_without_demand,
        "note": (
            "access_per_1000 为每人可获得容量乘 1000；"
            "本例 capacity 为床位数，因此单位为“每千人床位”。"
            "直线距离结果不能直接等同于路网可达性。"
        ),
    }


def write_outputs(
    result: TwoStepResult,
    summary: dict[str, object],
    output_dir: Path,
    prefix: str = "accessibility",
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{prefix}.csv"
    geojson_path = output_dir / f"{prefix}.geojson"
    summary_path = output_dir / f"{prefix}_summary.json"

    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fieldnames = [
            "community_id",
            "name",
            "lon",
            "lat",
            "population",
            "access_per_1000",
            "facilities_in_catchment",
            "nearest_facility_km",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in result.records:
            point = record.demand
            writer.writerow(
                {
                    "community_id": point.community_id,
                    "name": point.name,
                    "lon": point.lon,
                    "lat": point.lat,
                    "population": point.population,
                    "access_per_1000": f"{record.access_per_1000:.6f}",
                    "facilities_in_catchment": record.facilities_in_catchment,
                    "nearest_facility_km": f"{record.nearest_facility_km:.6f}",
                }
            )

    features = []
    for record in result.records:
        point = record.demand
        features.append(
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [point.lon, point.lat]},
                "properties": {
                    "community_id": point.community_id,
                    "name": point.name,
                    "population": point.population,
                    "access_per_1000": round(record.access_per_1000, 6),
                    "facilities_in_catchment": record.facilities_in_catchment,
                    "nearest_facility_km": round(record.nearest_facility_km, 6),
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

    return {"csv": csv_path, "geojson": geojson_path, "summary": summary_path}


# ---------------------------------------------------------------------------
# 快速出图（matplotlib 可选；没有安装时只跳过地图）
# ---------------------------------------------------------------------------

def _class_edges(values: Sequence[float]) -> list[float]:
    positive = sorted(value for value in values if value > 0)
    if not positive:
        return [0.0, 1.0]
    edges = [0.0]
    for q in (0.2, 0.4, 0.6, 0.8):
        value = _quantile(positive, q)
        if value > edges[-1] + 1e-9:
            edges.append(value)
    if edges[-1] >= positive[-1] - 1e-9:
        edges[-1] = positive[-1]
    else:
        edges.append(positive[-1])
    return edges


def _class_colors(count: int) -> list[str]:
    if count <= 0:
        return []
    if count == 1:
        return ["#3182bd"]
    return [
        BLUE_STEPS[round(index * (len(BLUE_STEPS) - 1) / (count - 1))]
        for index in range(count)
    ]


def plot_map(
    result: TwoStepResult,
    facilities: Sequence[Facility],
    summary: dict[str, object],
    png_path: Path,
    svg_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        print("提示：未安装 matplotlib，已跳过地图输出；pip install matplotlib 后重跑即可。", file=sys.stderr)
        return

    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False

    records = result.records
    lons = [record.demand.lon for record in records]
    lats = [record.demand.lat for record in records]
    values = [record.access_per_1000 for record in records]
    edges = _class_edges(values)
    colors = _class_colors(len(edges) - 1)
    point_colors = []
    for value in values:
        if value <= 0:
            point_colors.append(ZERO_COLOR)
            continue
        index = bisect.bisect_right(edges, value) - 1
        index = min(max(index, 0), len(colors) - 1)
        point_colors.append(colors[index])

    fig, ax = plt.subplots(figsize=(11.2, 8.4), dpi=160)
    fig.patch.set_facecolor("#f7f9fb")
    ax.set_facecolor("#fbfcfd")

    ax.scatter(lons, lats, s=48, c=point_colors, edgecolors="white", linewidths=0.55, zorder=3)
    ax.scatter(
        [facility.lon for facility in facilities],
        [facility.lat for facility in facilities],
        marker="^",
        s=105,
        c="#1a1a1a",
        edgecolors="white",
        linewidths=0.9,
        zorder=5,
    )

    span_x = max(lons) - min(lons)
    span_y = max(lats) - min(lats)
    pad_x, pad_y = span_x * 0.06, span_y * 0.07
    ax.set_xlim(min(lons) - pad_x, max(lons) + pad_x)
    ax.set_ylim(min(lats) - pad_y, max(lats) + pad_y)
    ax.set_aspect(1 / math.cos(math.radians(sum(lats) / len(lats))))
    ax.grid(color="#c9d3dd", alpha=0.35, linewidth=0.6, linestyle="--")
    ax.set_xlabel("经度", fontsize=9, color="#444")
    ax.set_ylabel("纬度", fontsize=9, color="#444")
    ax.tick_params(labelsize=8, colors="#555")
    for spine in ax.spines.values():
        spine.set_color("#c9d3dd")

    # 图例：0 值单独用红色，正数按分位数分五档
    handles = []
    if any(value <= 0 for value in values):
        handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=7,
                markerfacecolor=ZERO_COLOR, markeredgecolor="white", label="0（阈值内无机构）",
            )
        )
    for index, color in enumerate(colors):
        low, high = edges[index], edges[index + 1]
        handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=7,
                markerfacecolor=color, markeredgecolor="white",
                label=f"{low:.1f}–{high:.1f}",
            )
        )
    handles.append(
        Line2D(
            [0], [0], marker="^", linestyle="none", markersize=8,
            markerfacecolor="#1a1a1a", markeredgecolor="white", label="医疗机构",
        )
    )
    legend = ax.legend(
        handles=handles,
        title="每千人床位",
        loc="center left",
        bbox_to_anchor=(1.015, 0.5),
        frameon=True,
        framealpha=0.96,
        fontsize=8.5,
        title_fontsize=9.5,
        borderpad=0.8,
        labelspacing=0.65,
    )
    legend.get_frame().set_edgecolor("#d0d7de")

    # 指北针放在绘图区外，避免遮挡边缘的需求点
    ax.annotate(
        "N", xy=(1.035, 0.985), xytext=(1.035, 0.905),
        xycoords="axes fraction", ha="center", va="center", fontsize=10, fontweight="bold",
        annotation_clip=False,
        arrowprops={"arrowstyle": "-|>", "color": "#333", "lw": 1.2},
    )

    # 比例尺（按研究区平均纬度换算）
    km_per_degree_lon = 111.32 * math.cos(math.radians(sum(lats) / len(lats)))
    bar_km = 5 if span_x * km_per_degree_lon >= 10 else 2
    bar_degree = bar_km / km_per_degree_lon
    x0 = min(lons) + span_x * 0.035
    y0 = min(lats) + span_y * 0.035
    ax.plot([x0, x0 + bar_degree], [y0, y0], color="#222", linewidth=3.5, solid_capstyle="butt", zorder=6)
    ax.text(x0 + bar_degree / 2, y0 + span_y * 0.012, f"{bar_km} km", ha="center", va="bottom", fontsize=8, color="#333")

    fig.text(0.075, 0.955, "医疗资源空间可达性：两步移动搜索法（2SFCA）",
             fontsize=18, fontweight="bold", color="#17212b", ha="left")
    decay_label = "二值衰减" if summary["decay"] == "binary" else "高斯衰减"
    fig.text(
        0.075, 0.917,
        f"直线距离 {summary['threshold_km']:g} km · {decay_label} · 颜色越深，人均床位资源越高",
        fontsize=10, color="#5b6773", ha="left",
    )
    fig.text(
        0.075, 0.025,
        "模拟数据，仅用于方法演示；正式研究请用路网出行时间、真实人口和机构容量替换。",
        fontsize=8, color="#777", ha="left",
    )
    ratio = summary.get("p90_p10_ratio")
    ratio_text = f"{ratio:.1f}" if isinstance(ratio, (int, float)) else "—"
    ax.text(
        1.015, 0.015,
        f"需求点  {summary['demand_count']} 个\n"
        f"医疗机构  {summary['facility_count']} 家\n"
        f"覆盖人口  {summary['coverage_rate'] * 100:.1f}%\n"
        f"中位数  {summary['median_access_per_1000']:.2f}/千人\n"
        f"P90/P10  {ratio_text}",
        transform=ax.transAxes, va="bottom", ha="left", fontsize=8.2, color="#333",
        linespacing=1.55,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "#d0d7de", "alpha": 0.95},
    )
    fig.subplots_adjust(left=0.075, right=0.79, top=0.875, bottom=0.09)

    fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="两步移动搜索法（2SFCA）：计算医疗资源空间可达性并导出地图。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--demand", type=Path, default=DEFAULT_DEMAND, help="需求点 CSV")
    parser.add_argument("--facilities", type=Path, default=DEFAULT_FACILITIES, help="医疗机构 CSV")
    parser.add_argument("--threshold-km", type=float, default=5.0, help="服务阈值 d0（km）")
    parser.add_argument(
        "--decay", choices=("binary", "gaussian"), default="binary",
        help="距离衰减：binary 为阈值内等权，gaussian 为截断高斯",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT, help="输出目录")
    parser.add_argument("--prefix", default="accessibility", help="输出文件前缀")
    parser.add_argument("--no-map", action="store_true", help="只算表格，不画地图")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        demand = load_demand(args.demand)
        facilities = load_facilities(args.facilities)
        result = compute_2sfca(demand, facilities, args.threshold_km, args.decay)
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2

    summary = summarize(result, args.threshold_km, args.decay)
    summary["total_capacity"] = math.fsum(facility.capacity for facility in facilities)
    summary["capacity_per_1000_residents"] = (
        summary["total_capacity"] / summary["total_population"] * 1000
    )
    paths = write_outputs(result, summary, args.output_dir, args.prefix)
    if not args.no_map:
        paths["map_png"] = args.output_dir / f"{args.prefix}_map.png"
        paths["map_svg"] = args.output_dir / f"{args.prefix}_map.svg"
        plot_map(result, facilities, summary, paths["map_png"], paths["map_svg"])

    print(
        f"完成：{len(demand)} 个需求点、{len(facilities)} 家机构；"
        f"覆盖人口 {summary['coverage_rate'] * 100:.1f}%，"
        f"每千人床位中位数 {summary['median_access_per_1000']:.2f}。"
    )
    if result.facilities_without_demand:
        print(
            f"提示：{result.facilities_without_demand} 家机构在 {args.threshold_km:g} km 内没有需求点，"
            "其两步比值为 0；请检查研究区边界或机构是否落在区外。",
            file=sys.stderr,
        )
    print("输出文件：")
    for name, path in paths.items():
        print(f"  {name:10s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
