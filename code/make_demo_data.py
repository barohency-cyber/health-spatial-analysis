#!/usr/bin/env python3
"""生成一套虚构的“临江市”示例数据，方便零门槛跑通 2SFCA。

数据全部为模拟值，机构名称、位置、人口、床位均不对应真实地点。
运行：python code/make_demo_data.py
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"

# 虚构研究区范围（约 17 km × 13 km）
WEST, EAST = 120.100, 120.276
SOUTH, NORTH = 30.230, 30.350
CENTER_LON, CENTER_LAT = 120.188, 30.290

STREETS = [
    "青禾街道", "望江街道", "云岭街道", "南塘街道",
    "北辰街道", "东港街道", "西林街道", "长桥街道",
    "新安街道", "江湾街道", "凤凰街道", "石湖街道",
]

# 15 家虚构医疗机构：市中心 5 家 + 外围 8 家 + 远郊 2 家
FACILITIES = [
    ("F01", "临江市第一人民医院", 120.188, 30.290, 1500),
    ("F02", "临江市中医院", 120.176, 30.282, 850),
    ("F03", "临江市妇幼保健院", 120.200, 30.298, 600),
    ("F04", "临江市第二人民医院", 120.182, 30.304, 500),
    ("F05", "临江市康复医院", 120.196, 30.280, 400),
    ("F06", "青禾社区卫生服务中心", 120.143, 30.290, 260),
    ("F07", "望江社区卫生服务中心", 120.233, 30.290, 280),
    ("F08", "云岭社区卫生服务中心", 120.188, 30.332, 240),
    ("F09", "南塘社区卫生服务中心", 120.188, 30.248, 220),
    ("F10", "北辰社区卫生服务中心", 120.153, 30.320, 180),
    ("F11", "东港社区卫生服务中心", 120.223, 30.320, 200),
    ("F12", "西林社区卫生服务中心", 120.153, 30.260, 170),
    ("F13", "长桥社区卫生服务中心", 120.223, 30.260, 190),
    ("F14", "江湾镇卫生院", 120.120, 30.340, 120),
    ("F15", "石湖乡卫生院", 120.256, 30.240, 100),
]


def population_at(lon: float, lat: float, rng: random.Random) -> int:
    """市中心人口密度高、外围低，并加入小幅随机波动。"""
    dx_km = (lon - CENTER_LON) * 111.32 * math.cos(math.radians(CENTER_LAT))
    dy_km = (lat - CENTER_LAT) * 110.57
    distance_km = math.hypot(dx_km, dy_km)
    base = 1_200 + 18_000 * math.exp(-((distance_km / 3.4) ** 2))
    value = base * rng.uniform(0.86, 1.14)
    return max(300, int(round(value, -1)))


def make_demand() -> list[dict[str, object]]:
    rng = random.Random(20261028)
    rows: list[dict[str, object]] = []
    n_cols, n_rows = 12, 10
    for col in range(n_cols):
        lon0 = WEST + (EAST - WEST) * col / (n_cols - 1)
        for row in range(n_rows):
            lat0 = SOUTH + (NORTH - SOUTH) * row / (n_rows - 1)
            lon = lon0 + rng.uniform(-0.0018, 0.0018)
            lat = lat0 + rng.uniform(-0.0018, 0.0018)
            idx = col * n_rows + row + 1
            street = STREETS[col % len(STREETS)]
            rows.append(
                {
                    "community_id": f"C{idx:03d}",
                    "name": f"{street}{row + 1:02d}社区",
                    "lon": round(lon, 6),
                    "lat": round(lat, 6),
                    "population": population_at(lon, lat, rng),
                }
            )
    return rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    demand = make_demand()
    facilities = [
        {
            "facility_id": fid,
            "name": name,
            "lon": lon,
            "lat": lat,
            "capacity": capacity,
        }
        for fid, name, lon, lat, capacity in FACILITIES
    ]
    write_csv(
        DATA_DIR / "demand.csv",
        ["community_id", "name", "lon", "lat", "population"],
        demand,
    )
    write_csv(
        DATA_DIR / "facilities.csv",
        ["facility_id", "name", "lon", "lat", "capacity"],
        facilities,
    )
    print(f"已生成 {len(demand)} 个需求点、{len(facilities)} 家医疗机构：{DATA_DIR}")


if __name__ == "__main__":
    main()
