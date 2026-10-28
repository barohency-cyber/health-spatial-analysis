# 医疗资源空间可达性：2SFCA 极简示例

面向 **卫生事业管理 / 公共卫生 / 数字治理** 方向的 2SFCA（两步移动搜索法）最小项目：
不用装 ArcGIS，不依赖 GeoPandas，准备好两张 CSV 就能算；算完直接出 CSV、GeoJSON 和论文汇报用地图。

> 仓库内全部为**模拟数据**：虚构的“临江市”、120 个社区、15 家医疗机构。
> 机构名、位置、人口、床位均不对应真实地点，可直接公开上传 GitHub。

![医疗资源空间可达性地图](output/accessibility_map.png)

---

## 1. 30 秒跑通

需要 Python 3.10+。核心计算只用标准库，画图需要 matplotlib。

```powershell
python -m pip install -r requirements.txt
python code/two_step_fca.py
```

看到类似输出即成功：

```text
完成：120 个需求点、15 家机构；覆盖人口 98.0%，每千人床位中位数 17.86。
```

不想画图、只要表格：

```powershell
python code/two_step_fca.py --no-map
```

跑单元测试：

```powershell
python -m unittest discover -s tests -v
```

---

## 2. 它算了什么

2SFCA 分两步，核心是同时考虑两件事：

1. **供给端竞争**：一家医院不只服务门口的人，周边居民越多，每个人能分到的床位越少；
2. **需求端多机构叠加**：一个社区可能同时处在几家医院的服务范围内，可达性是这些医院供给能力的总和。

### 第一步：算每家机构的“供需比”

```text
R_j = S_j / Σ(P_i × W_ij)      （只对 d_ij ≤ d0 的需求点求和）
```

- `S_j`：机构 j 的容量，例如床位数、医生数；
- `P_i`：需求点 i 的人口；
- `d_ij`：需求点 i 到机构 j 的距离；
- `d0`：服务阈值，例如 5 km；
- `W_ij`：距离衰减权重。`binary` 表示阈值内权重都为 1；`gaussian` 为截断高斯核。

`R_j` 可以理解为：**这家机构服务范围内，平均每个人能分到多少资源**。

### 第二步：算每个需求点的可达性

```text
A_i = Σ(R_j × W_ij)           （只对 d_ij ≤ d0 的机构求和）
```

`A_i` 的单位是“容量 / 人”。代码乘以 1000，输出 `access_per_1000`，若 `capacity` 用床位数，单位就是**每千人床位**。

> 这不是“居民实际利用了多少医疗资源”，而是“在给定供给、人口和距离阈值下，居民理论上能获得的资源”。论文里要区分 **potential access** 与 **realized utilization**。

经典方法表述见：

- Luo W, Wang F. Measures of spatial accessibility to health care in a GIS environment: synthesis and a case study in the Chicago region. *Environment and Planning B: Planning and Design*, 2003, 30(6): 865–884.

---

## 3. 输入数据格式

只需要两张 CSV，坐标用 **WGS84 经纬度（EPSG:4326）**。

### `data/demand.csv`：需求点（社区 / 街道 / 网格）

| community_id | name | lon | lat | population |
|---|---|---|---|---|
| C001 | 青禾街道01社区 | 120.101701 | 30.230211 | 1280 |
| C002 | 青禾街道02社区 | 120.101575 | 30.243615 | 1110 |

- `community_id`：唯一编号，后面和边界图层连接要用；
- `population`：该需求点的常住人口（或你研究口径下的人口）；
- 一个社区一行，不要重复。

### `data/facilities.csv`：供给点（医院 / 卫生院 / 社区卫生服务中心）

| facility_id | name | lon | lat | capacity |
|---|---|---|---|---|
| F01 | 临江市第一人民医院 | 120.188 | 30.290 | 1500 |
| F06 | 青禾社区卫生服务中心 | 120.143 | 30.290 | 260 |

- `capacity`：机构容量。**整列必须同一单位**，不能一半填床位数、一半填门诊量；
- `capacity` 必须大于 0；没有供给能力的机构请先从供给点表中剔除；
- 如果研究“医疗资源可及性”，常用床位数或执业（助理）医师数；
- 机构可以落在研究区外，但要注意边界效应（见第 7 节）。

### 从 QGIS 导出

1. 把需求点/机构点图层的坐标系设为 `EPSG:4326`（WGS 84）；
2. 右键图层 → **导出 → 要素另存为**；
3. 格式选 **CSV**，在“几何”选项中选择 `AS_XY`；
4. 把导出的字段重命名为上表字段名（可在 QGIS 字段计算器里新建 `lon = $x`、`lat = $y`）；
5. 文件编码选 UTF-8（脚本也兼容带 BOM 的 UTF-8）。

如果你只有行政区面图层：

- 用 **矢量 → 几何工具 → 质心** 生成点，再把人口属性带过去；
- 或直接用面图层字段计算器取 `$x`、`$y`；
- 面转点会引入 MAUP 问题，论文里要说明用的是“社区/街道质心”还是“人口加权中心”。

---

## 4. 换成自己的研究

```powershell
python code/two_step_fca.py `
  --demand data/my_demand.csv `
  --facilities data/my_facilities.csv `
  --threshold-km 5 `
  --output-dir output/my_city
```

常用参数：

| 参数 | 说明 |
|---|---|
| `--threshold-km` | 服务阈值，单位 km；示例默认 5 |
| `--decay` | `binary`（阈值内等权）或 `gaussian`（截断高斯） |
| `--output-dir` | 输出目录 |
| `--prefix` | 输出文件前缀，默认 `accessibility` |
| `--no-map` | 只输出表格，不画图 |

做阈值敏感性分析（论文里很加分）：

```powershell
python code/two_step_fca.py --threshold-km 3 --output-dir output/d3
python code/two_step_fca.py --threshold-km 5 --output-dir output/d5
python code/two_step_fca.py --threshold-km 10 --output-dir output/d10
```

然后比较三套结果的排序是否稳定（例如 Spearman 相关系数）。**不要只挑一个看起来最显著的结果写进论文。**

---

## 5. 输出怎么用

| 文件 | 用途 |
|---|---|
| `output/accessibility.csv` | 主结果表；按 `community_id` 连回 QGIS 面图层做分级设色 |
| `output/accessibility.geojson` | 可直接拖进 QGIS / kepler.gl / GitHub 预览 |
| `output/accessibility_summary.json` | 覆盖人口、均值、中位数、P10/P90 等汇总指标 |
| `output/accessibility_map.png` | 快速检查用图；正式论文图建议在 QGIS 里重做 |
| `output/accessibility_map.svg` | 矢量图，可进 Illustrator / Inkscape 调整 |

字段解释：

- `access_per_1000`：2SFCA 可达性，数值越高表示人均资源越多；
- `facilities_in_catchment`：阈值内机构数量；为 0 表示“医疗资源空白区”；
- `nearest_facility_km`：最近机构的直线距离，只作辅助，不等于可达性。

### 在 QGIS 里做论文图

1. 把 `accessibility.geojson` 拖进 QGIS；
2. 用 **按字段连接**（Join attributes by field value）把 `access_per_1000` 连到社区/街道面图层，连接字段为 `community_id`；
3. 图层属性 → 符号化 → **分级** → 字段选 `access_per_1000`，分 5 类，选单色系色带；
4. 新建打印布局，加比例尺、指北针、图例、数据来源；
5. 导出 300 dpi PNG 或 PDF。

### 结果解读示例

```text
研究区 120 个社区中，98.0% 的人口在 5 km 直线距离内至少能到达一家医疗机构；
但可达性分布并不均衡：每千人床位中位数为 17.86，P10 仅为 2.65，P90 为 20.18，
P90/P10 约为 7.6。高值区集中在市中心，低值区和 0 值区主要出现在东北、西南边缘，
提示这些区域存在医疗资源覆盖不足的风险。
```

注意：`P90/P10` 在 `P10 = 0` 时没有定义，脚本输出 `null`，不要硬写成无穷大。

### 一个容易被问倒的点

在封闭研究区、且每家机构都至少服务到一个需求点的假设下，**人口加权平均可达性会等于“总容量 ÷ 总人口”**。这是 2SFCA 的守恒性质，不是bug。所以：

- 不要只报“平均每千人床位”；
- 要同时看中位数、P10/P90、覆盖率和空间分布；
- 平均值相同，不代表资源分配公平。

---

## 6. 论文与复试表述模板

### 方法段

```text
采用两步移动搜索法（2SFCA）测度医疗资源空间可达性。以社区为需求单元，
以各级医疗机构为供给点，以机构床位数衡量供给容量，以社区常住人口衡量需求规模。
参考已有研究将服务阈值设为 5 km，采用二值距离衰减函数，分别计算每家机构的
供需比 R_j 和每个社区的可达性指数 A_i；A_i 越高，表示该社区在竞争后可获得的
人均医疗资源越多。进一步以人口加权分位数和 P90/P10 描述可达性不平等程度。
```

### 图注

```text
图 X 研究区医疗资源空间可达性分布（2SFCA，服务阈值 5 km）
注：颜色越深表示每千人床位资源越高；红色表示 5 km 内无医疗机构。
基于模拟数据绘制，仅用于方法示例。
```

### 复试 30 秒口述

```text
我用 2SFCA 做医疗资源可达性。第一步，以每家医院为中心，把 5 公里内的社区人口
加总，得到这家医院的供需比；第二步，对每个社区，把它能覆盖到的所有医院的供需比
加起来，就是该社区的人均资源可得性。这样既考虑了居民之间的竞争，也考虑了
跨社区就医。结果用覆盖人口、P10/P90 和空间分布来解读，而不是只看平均值。
这个方法可以用 QGIS 出图，也可以用 Python 复现。
```

---

## 7. 常见坑（医管学生高频踩雷）

1. **把直线距离说成路网可达性。**  
   本代码是 haversine 直线距离。论文若写“路网可达性”，必须换成路网出行时间/距离，否则方法表述不成立。

2. **阈值拍脑袋。**  
   5 km 只是示例。应按研究问题、城市尺度、居民就医行为或文献设定，并做 3/5/10 km 敏感性分析。

3. **边界效应。**  
   研究区外的医院没纳入，边缘社区会被低估；研究区外的居民没纳入，边缘医院会被高估。解决：研究区外扩缓冲带、纳入周边机构，或明确说明假设。

4. **容量单位混用。**  
   床位数、医师数、门诊量不能混在一列。`capacity` 整列必须同单位；换成医师数后，结论要相应改写。

5. **人口口径不一致。**  
   常住人口、户籍人口、老年人口不要混用；需求点人口要与研究问题匹配，并注明来源和年份。

6. **坐标系搞错。**  
   脚本要的是经纬度（WGS84）。把投影坐标（米）当经纬度传进去，距离会完全错误。QGIS 导出前先确认 CRS。

7. **零值不等于“没有需求”。**  
   `facilities_in_catchment = 0` 是“阈值内没有机构”，不是“没有人”。两者要分开解释。

8. **把可达性当成利用率。**  
   可达性高不代表居民真的去看病；就医行为还受医保、价格、信任、交通、服务质量影响。

9. **MAUP 与生态谬误。**  
   社区/街道是人为划分的，换一种划分方式结果可能变化；群体层面的关系不能直接推到个人。

10. **数据合规。**  
    真实边界 shp、地址、问卷个体数据不要随便传公开仓库；公开仓库只放模拟数据、聚合结果或脱敏数据。敏感数据放本地私有仓库，并在 `.gitignore` 中排除。

11. **公式混用。**  
    2SFCA、E2SFCA、3SFCA、Huff 模型的公式不同，不要把一种方法的权重套到另一种方法上还沿用原名。

---

## 8. 仓库结构

```text
.
├── code/
│   ├── two_step_fca.py       # 2SFCA 核心 + CSV/GeoJSON/地图输出
│   └── make_demo_data.py     # 重新生成模拟数据
├── data/
│   ├── demand.csv            # 需求点（模拟）
│   └── facilities.csv        # 医疗机构（模拟）
├── output/
│   ├── accessibility.csv
│   ├── accessibility.geojson
│   ├── accessibility_summary.json
│   ├── accessibility_map.png
│   └── accessibility_map.svg
├── tests/
│   └── test_two_step_fca.py  # 手算例子校验公式
├── README.md
└── requirements.txt
```

重新生成模拟数据（会覆盖 `data/` 下两个 CSV）：

```powershell
python code/make_demo_data.py
python code/two_step_fca.py
```

---

## 9. 下一步可以怎么升级

- **换路网**：用 QGIS 的 ORS Tools、OSRM，或高德/百度路径规划 API 生成 OD 时间矩阵，把 `haversine_km` 换成出行时间；
- **换面数据**：用 GeoPandas 把社区面质心化，或直接做面到面 OD；
- **空间自相关**：用 PySAL/esda 对 `access_per_1000` 做全局/局部 Moran's I，识别高-高、低-低聚集；
- **空间回归**：以可达性为自变量，做 SLM/SEM，解释医保利用或健康结果；
- **方法扩展**：E2SFCA、3SFCA、Huff 模型；但先把手里的 2SFCA 假设、阈值和边界效应讲清楚，比堆模型更加分。

GitHub 仓库描述和 topics 可以写：

```text
description: 2SFCA-based spatial accessibility analysis for healthcare resources
topics: health-gis, 2sfca, healthcare-accessibility, public-health, qgis
```

---

## 10. 数据与许可

- 仓库内数据为程序生成的模拟数据，不涉及真实人口、机构或边界；
- 示例代码可自由修改用于课程作业、论文复现和复试展示；
- 正式研究请替换为合规数据，并在论文中写明数据来源、年份、口径和伦理审批（如适用）。
