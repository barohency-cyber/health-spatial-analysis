#!/usr/bin/env python3
"""生成演示用的“就医利用”结果变量（模拟数据）。

依赖 `output/accessibility.csv` 中的可达性结果，生成：
    data/outcomes.csv
字段：
    community_id, outpatient_utilization, elderly_share, income_index

公式（仅用于演示）：
    可达性效应 = 0.10 + 0.35 × 市中心衰减权重
    两周就诊率 = 0.15 + 可达性效应 × 可达性标准化 + 0.30 × 老年人口占比
                 + 0.10 × 收入指数 + 空间相关扰动
所有数值均为模拟，不对应真实社区。
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACCESS = ROOT / "output" / "accessibility.csv"
OUTCOME = ROOT / "data" / "outcomes.csv"
EARTH_RADIUS_KM = 6371.0088


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    )
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(min(1.0, max(0.0, a))))


def main() -> None:
    with ACCESS.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"{ACCESS} 没有数据；请先运行 python code/two_step_fca.py")

    points = [
        {
            "id": row["community_id"],
            "name": row["name"],
            "lon": float(row["lon"]),
            "lat": float(row["lat"]),
            "access": float(row["access_per_1000"]),
        }
        for row in rows
    ]
    center_lon = sum(point["lon"] for point in points) / len(points)
    center_lat = sum(point["lat"] for point in points) / len(points)
    max_distance = max(
        haversine_km(center_lon, center_lat, point["lon"], point["lat"])
        for point in points
    )
    max_access = max(point["access"] for point in points) or 1.0

    # 6 个最近邻，用于生成空间相关的扰动项
    neighbors: list[list[int]] = []
    for i, origin in enumerate(points):
        distances = []
        for j, target in enumerate(points):
            if i == j:
                continue
            distances.append(
                (haversine_km(origin["lon"], origin["lat"], target["lon"], target["lat"]), j)
            )
        distances.sort(key=lambda item: (item[0], item[1]))
        neighbors.append([j for _, j in distances[:6]])

    rng = random.Random(20261028)
    raw_noise = [rng.gauss(0, 0.035) for _ in points]
    elderly = [0.16 + rng.gauss(0, 0.035) for _ in points]
    income = [1.0 + rng.gauss(0, 0.16) for _ in points]

    output_rows = []
    for i, point in enumerate(points):
        distance = haversine_km(center_lon, center_lat, point["lon"], point["lat"])
        spatial_noise = 0.55 * raw_noise[i] + 0.45 * (
            sum(raw_noise[j] for j in neighbors[i]) / len(neighbors[i])
        )
        access_norm = point["access"] / max_access
        # 让可达性效应随位置变化：市中心更强，边缘更弱，用于演示 GWR
        centrality = math.exp(-((distance / (max_distance * 0.55)) ** 2))
        access_effect = 0.10 + 0.35 * centrality
        utilization = (
            0.15
            + access_effect * access_norm
            + 0.30 * elderly[i]
            + 0.10 * income[i]
            + spatial_noise
            + 0.015 * (distance / max_distance)
        )
        utilization = min(0.88, max(0.05, utilization))
        output_rows.append(
            {
                "community_id": point["id"],
                "outpatient_utilization": f"{utilization:.4f}",
                "elderly_share": f"{min(0.35, max(0.08, elderly[i])):.4f}",
                "income_index": f"{min(1.5, max(0.6, income[i])):.4f}",
            }
        )

    OUTCOME.parent.mkdir(parents=True, exist_ok=True)
    with OUTCOME.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "community_id",
                "outpatient_utilization",
                "elderly_share",
                "income_index",
            ],
        )
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"已生成 {len(output_rows)} 条模拟就医利用数据：{OUTCOME}")


if __name__ == "__main__":
    main()
