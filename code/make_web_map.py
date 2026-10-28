#!/usr/bin/env python3
"""生成 GitHub Pages 用的交互式可达性地图（Leaflet，无额外 Python 依赖）。

前置：先运行 `python code/two_step_fca.py` 生成 output/accessibility.csv，
然后运行：

    python code/make_web_map.py

输出默认写入 docs/index.html，可直接用浏览器打开，也可以部署到 GitHub Pages。
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import sys
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACCESS = ROOT / "output" / "accessibility.csv"
DEFAULT_FACILITIES = ROOT / "data" / "facilities.csv"
DEFAULT_SUMMARY = ROOT / "output" / "accessibility_summary.json"
DEFAULT_OUTPUT = ROOT / "docs" / "index.html"
ZERO_COLOR = "#d73027"
BLUE_STEPS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]


def read_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except FileNotFoundError as exc:
        raise ValueError(f"找不到文件：{path}；请先运行 python code/two_step_fca.py") from exc


def quantile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * q
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def class_edges(values: Sequence[float]) -> list[float]:
    positive = sorted(value for value in values if value > 0)
    if not positive:
        return [0.0, 1.0]
    edges = [0.0]
    for q in (0.2, 0.4, 0.6, 0.8):
        value = quantile(positive, q)
        if value > edges[-1] + 1e-9:
            edges.append(value)
    if edges[-1] >= positive[-1] - 1e-9:
        edges[-1] = positive[-1]
    else:
        edges.append(positive[-1])
    return edges


def class_colors(count: int) -> list[str]:
    if count <= 0:
        return []
    if count == 1:
        return ["#3182bd"]
    return [
        BLUE_STEPS[round(index * (len(BLUE_STEPS) - 1) / (count - 1))]
        for index in range(count)
    ]


def legend_html(edges: Sequence[float], colors: Sequence[str], has_zero: bool) -> str:
    rows = []
    if has_zero:
        rows.append(
            '<div class="legend-row"><span class="dot" style="background:'
            + ZERO_COLOR
            + '"></span><span>0（阈值内无机构）</span></div>'
        )
    for index, color in enumerate(colors):
        low, high = edges[index], edges[index + 1]
        rows.append(
            f'<div class="legend-row"><span class="dot" style="background:{color}"></span>'
            f"<span>{low:.1f}–{high:.1f}</span></div>"
        )
    rows.append('<div class="legend-row"><span class="hospital-symbol">✚</span><span>医疗机构</span></div>')
    return (
        '<div class="panel legend-panel"><div class="legend-title">每千人床位</div>'
        + "".join(rows)
        + '<div class="legend-note">圆点大小 ≈ 需求人口</div></div>'
    )


def title_html(summary: dict[str, object], n_records: int, n_facilities: int) -> str:
    threshold = summary.get("threshold_km", "—")
    decay = "二值衰减" if summary.get("decay") == "binary" else "高斯衰减"
    coverage = summary.get("coverage_rate")
    median = summary.get("median_access_per_1000")
    coverage_text = f"{coverage * 100:.1f}%" if isinstance(coverage, (int, float)) else "—"
    median_text = f"{median:.2f}" if isinstance(median, (int, float)) else "—"
    return f"""
    <div class="panel title-panel">
      <div class="title">医疗资源空间可达性 · 2SFCA</div>
      <div class="subtitle">直线距离 {threshold} km · {decay} · 模拟数据</div>
      <div class="subtitle">覆盖人口 {coverage_text} · 每千人床位中位数 {median_text}</div>
      <div class="hint">{n_records} 个社区 · {n_facilities} 家机构 · 点击圆点查看指标</div>
    </div>
    """


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>医疗资源空间可达性 · 2SFCA</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<style>
  html, body { height: 100%; margin: 0; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Microsoft YaHei", "Noto Sans CJK SC", sans-serif; color: #17212b; }
  #map { height: 100%; width: 100%; background: #eef2f6; }
  .panel { position: fixed; z-index: 1000; background: rgba(255,255,255,.96);
           border: 1px solid #d0d7de; border-radius: 10px;
           box-shadow: 0 3px 14px rgba(0,0,0,.10); }
  .title-panel { top: 14px; left: 52px; max-width: min(390px, calc(100vw - 76px)); padding: 12px 15px; }
  .title { font-size: 17px; font-weight: 800; line-height: 1.35; }
  .subtitle { font-size: 12px; color: #5b6773; margin-top: 4px; }
  .hint { font-size: 11px; color: #8a949e; margin-top: 6px; }
  .legend-panel { bottom: 24px; left: 24px; padding: 10px 12px; font-size: 12px; line-height: 1.5; }
  .legend-title { font-weight: 700; margin-bottom: 5px; }
  .legend-row { display: flex; align-items: center; gap: 7px; margin: 3px 0; }
  .dot { width: 11px; height: 11px; border-radius: 50%; display: inline-block; flex: 0 0 auto; }
  .hospital-symbol { color: #8b1a1a; width: 11px; text-align: center; font-size: 14px; line-height: 11px; }
  .legend-note { color: #8a949e; font-size: 11px; margin-top: 6px; }
  .popup { font-family: -apple-system, BlinkMacSystemFont, "Microsoft YaHei", sans-serif; min-width: 190px; }
  .popup-title { font-size: 14px; font-weight: 700; color: #17212b; }
  .popup-id { font-size: 11px; color: #8a949e; margin-bottom: 5px; }
  .popup-body { font-size: 12px; line-height: 1.7; }
  .hospital-pin { width: 22px; height: 22px; border-radius: 50%; background: #8b1a1a; color: #fff;
                  border: 2px solid #fff; box-shadow: 0 1px 5px rgba(0,0,0,.35);
                  display: flex; align-items: center; justify-content: center; font-size: 13px; font-weight: 700; }
  .leaflet-control-layers { border-radius: 8px; }
  .grid-label { background: transparent; color: #7b8794; font-size: 10px; white-space: nowrap;
                text-shadow: 0 0 3px #fff, 0 0 3px #fff; }
  @media (max-width: 640px) {
    .title-panel { left: 52px; right: 10px; max-width: none; padding: 9px 11px; }
    .title { font-size: 14px; }
    .subtitle, .hint { font-size: 10px; }
    .legend-panel { bottom: 12px; left: 12px; padding: 8px 10px; font-size: 11px; }
  }
</style>
</head>
<body>
<div id="map"></div>
__TITLE__
__LEGEND__
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
const DATA = __DATA_JSON__;
const map = L.map("map", { zoomControl: true, preferCanvas: true });

const osm = L.tileLayer(
  "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
  { maxZoom: 19, attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors' }
);
const imagery = L.tileLayer(
  "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
  { maxZoom: 19, attribution: "Tiles &copy; Esri" }
);

// 离线示意底图：不依赖外部瓦片，断网也能看清点位、图例和比例尺
const allLats = DATA.records.map(r => r.lat);
const allLons = DATA.records.map(r => r.lon);
const dataBounds = L.latLngBounds(
  [Math.min(...allLats), Math.min(...allLons)],
  [Math.max(...allLats), Math.max(...allLons)]
);
const padded = dataBounds.pad(0.10);
const south = padded.getSouth(), west = padded.getWest();
const north = padded.getNorth(), east = padded.getEast();
const schematic = L.layerGroup();
L.rectangle(
  [[dataBounds.getSouth(), dataBounds.getWest()], [dataBounds.getNorth(), dataBounds.getEast()]],
  { color: "#c9d3dd", weight: 1.2, fillColor: "#f8fafc", fillOpacity: 1 }
).addTo(schematic);
for (let lon = Math.ceil(west / 0.04) * 0.04; lon <= east + 1e-9; lon += 0.04) {
  L.polyline([[south, lon], [north, lon]], {
    color: "#dde5ed", weight: 1, interactive: false
  }).addTo(schematic);
  L.marker([south + 0.008, lon], {
    icon: L.divIcon({
      className: "grid-label", html: lon.toFixed(2) + "°E",
      iconSize: [64, 14], iconAnchor: [32, 0]
    }),
    interactive: false
  }).addTo(schematic);
}
for (let lat = Math.ceil(south / 0.04) * 0.04; lat <= north + 1e-9; lat += 0.04) {
  L.polyline([[lat, west], [lat, east]], {
    color: "#dde5ed", weight: 1, interactive: false
  }).addTo(schematic);
  L.marker([lat, west + 0.006], {
    icon: L.divIcon({
      className: "grid-label", html: lat.toFixed(2) + "°N",
      iconSize: [54, 14], iconAnchor: [0, 7]
    }),
    interactive: false
  }).addTo(schematic);
}
schematic.addTo(map);

const accessLayer = L.layerGroup().addTo(map);
const facilityLayer = L.layerGroup().addTo(map);
const catchmentLayer = L.layerGroup();

L.control.layers(
  {
    "离线示意底图": schematic,
    "OpenStreetMap（需联网）": osm,
    "影像底图（需联网）": imagery
  },
  { "社区可达性": accessLayer, "医疗机构": facilityLayer, "5 km 服务范围": catchmentLayer },
  { collapsed: true, position: "topright" }
).addTo(map);
L.control.scale({ imperial: false, position: "bottomright" }).addTo(map);

const maxPopulation = Math.max(...DATA.records.map(r => r.population), 1);
const hospitalIcon = L.divIcon({
  className: "",
  html: '<div class="hospital-pin">✚</div>',
  iconSize: [24, 24],
  iconAnchor: [12, 12]
});

function esc(value) {
  return String(value).replace(/[&<>"']/g, ch => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
  }[ch]));
}

DATA.records.forEach(r => {
  const radius = 4 + 7 * Math.sqrt(r.population / maxPopulation);
  const marker = L.circleMarker([r.lat, r.lon], {
    radius: radius, color: "#ffffff", weight: 0.8,
    fillColor: r.color, fillOpacity: 0.88
  });
  marker.bindPopup(
    '<div class="popup"><div class="popup-title">' + esc(r.name) + '</div>' +
    '<div class="popup-id">' + esc(r.id) + '</div>' +
    '<div class="popup-body">每千人床位：<b>' + r.access.toFixed(2) + '</b><br>' +
    '需求人口：' + r.population.toLocaleString() + ' 人<br>' +
    '阈值内机构：' + r.facilities + ' 家<br>' +
    '最近机构直线距离：' + r.nearest.toFixed(2) + ' km</div></div>',
    { maxWidth: 300 }
  );
  marker.bindTooltip(r.name + " · " + r.access.toFixed(2) + "/千人", { direction: "top" });
  marker.addTo(accessLayer);
});

DATA.facilities.forEach(f => {
  const marker = L.marker([f.lat, f.lon], { icon: hospitalIcon, title: f.name });
  marker.bindPopup(
    '<div class="popup"><div class="popup-title">' + esc(f.name) + '</div>' +
    '<div class="popup-id">' + esc(f.id) + '</div>' +
    '<div class="popup-body">容量：<b>' + f.capacity.toLocaleString() + '</b> 床位</div></div>',
    { maxWidth: 280 }
  );
  marker.bindTooltip(f.name, { direction: "top" });
  marker.addTo(facilityLayer);

  L.circle([f.lat, f.lon], {
    radius: DATA.threshold_km * 1000, color: "#7b8794", weight: 1,
    opacity: 0.45, fill: false, dashArray: "4 4"
  }).addTo(catchmentLayer);
});

const bounds = L.latLngBounds(DATA.records.map(r => [r.lat, r.lon]));
map.fitBounds(bounds.pad(0.08));
</script>
</body>
</html>
"""


def build_map(
    access_rows: list[dict[str, str]],
    facility_rows: list[dict[str, str]],
    summary: dict[str, object],
    output: Path,
    threshold_km: float,
) -> None:
    records = []
    values = []
    for row in access_rows:
        value = float(row["access_per_1000"])
        values.append(value)
        records.append(
            {
                "id": row["community_id"],
                "name": row["name"],
                "lon": float(row["lon"]),
                "lat": float(row["lat"]),
                "population": float(row["population"]),
                "access": value,
                "facilities": int(float(row["facilities_in_catchment"])),
                "nearest": float(row["nearest_facility_km"]),
            }
        )
    facilities = [
        {
            "id": row["facility_id"],
            "name": row["name"],
            "lon": float(row["lon"]),
            "lat": float(row["lat"]),
            "capacity": float(row["capacity"]),
        }
        for row in facility_rows
    ]
    if not records:
        raise ValueError("accessibility.csv 没有数据行")
    if not facilities:
        raise ValueError("facilities.csv 没有数据行")

    edges = class_edges(values)
    colors = class_colors(len(edges) - 1)
    for record in records:
        if record["access"] <= 0:
            record["color"] = ZERO_COLOR
        else:
            index = bisect.bisect_right(edges, record["access"]) - 1
            record["color"] = colors[min(max(index, 0), len(colors) - 1)]

    payload = {
        "records": records,
        "facilities": facilities,
        "threshold_km": threshold_km,
    }
    data_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

    page = (
        HTML_TEMPLATE.replace("__TITLE__", title_html(summary, len(records), len(facilities)))
        .replace("__LEGEND__", legend_html(edges, colors, any(v <= 0 for v in values)))
        .replace("__DATA_JSON__", data_json)
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    (output.parent / ".nojekyll").write_text("", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="生成 Leaflet 交互地图（GitHub Pages / 本地浏览器）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--accessibility", type=Path, default=DEFAULT_ACCESS, help="可达性 CSV")
    parser.add_argument("--facilities", type=Path, default=DEFAULT_FACILITIES, help="机构 CSV")
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY, help="汇总 JSON")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="输出 HTML")
    parser.add_argument("--threshold-km", type=float, default=None, help="覆盖画布上的阈值标注")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        access_rows = read_rows(args.accessibility)
        facility_rows = read_rows(args.facilities)
        summary: dict[str, object] = {}
        if args.summary.exists():
            summary = json.loads(args.summary.read_text(encoding="utf-8"))
        threshold = args.threshold_km
        if threshold is None:
            threshold = float(summary.get("threshold_km", 5.0))
        build_map(access_rows, facility_rows, summary, args.output, threshold)
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2
    print(f"已生成交互地图：{args.output}")
    print("本地预览：python -m http.server 8000 -d docs，然后打开 http://localhost:8000/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
