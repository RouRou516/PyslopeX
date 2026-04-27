# pyslopex

2D 边坡平面滑动法和折线滑动法（传递系数法）边坡稳定性分析 Python 库。

## 安装

```bash
pip install pyslopex
```

依赖：`numpy >= 1.20`，`matplotlib >= 3.3`

## 快速开始

```python
from pyslopex import Slope, Material

s = Slope(height=10, angle=45)
s.set_materials(Material(unit_weight=20, friction_angle=30,
                         cohesion=10, depth_to_bottom=10))

# 平面滑动法
r = s.analyse_planar()
print(f"Fs = {r.fos:.3f}, θ = {r.critical_angle:.1f}°")

# 折线滑动法（自动搜索最危险滑动面）
r = s.analyse_polyline(auto_search=True, num_internal_points=2)
print(f"Fs = {r.fos:.3f}")
```

## 坐标系

```
y ↑
  │        ___________  y = H（坡顶平台）
  │       /
  │      /  坡面
  │     /
  │    /    坡脚(toe): (0, 0)
  │───/───────────────────→ x
  │  (0,0)     坡顶(crest): (L, H)
```

- x 轴水平向右，y 轴竖直向上
- 坡脚固定在原点 `(0, 0)`
- 坡顶位于 `(L, H)`，L 为坡面水平投影长度
- 滑动方向：从坡顶（右上）向坡脚（左下）

## 核心 CODE

### 1. 创建边坡

```python
# 方式一：指定坡角（度）
s = Slope(height=10, angle=45)

# 方式二：指定坡面水平投影长度
# 坡比 1:1.5 → 水平距离 = 1.5 × 高度
s = Slope(height=8, length=12)
```

| 参数 | 说明 |
|---|---|
| `height` | 坡高（m） |
| `angle` | 坡角（度），与 `length` 二选一 |
| `length` | 坡面水平投影长度（m），与 `angle` 二选一 |

### 2. 定义土层

```python
from pyslopex import Material

# 单层土
m = Material(unit_weight=20, friction_angle=30,
             cohesion=10, depth_to_bottom=10)
s.set_materials(m)

# 多层土：按 depth_to_bottom 自动排序，最后一层向下无限延伸
m1 = Material(unit_weight=20, friction_angle=35, cohesion=15,
              depth_to_bottom=4, name="密实砂")
m2 = Material(unit_weight=19, friction_angle=28, cohesion=12,
              depth_to_bottom=7, name="粉质粘土")
m3 = Material(unit_weight=18, friction_angle=22, cohesion=8,
              depth_to_bottom=12, name="软粘土")
s.set_materials(m1, m2, m3)
```

| 参数 | 说明 |
|---|---|
| `unit_weight` | 重度（kN/m³） |
| `friction_angle` | 有效内摩擦角 φ'（度） |
| `cohesion` | 有效粘聚力 c'（kPa） |
| `depth_to_bottom` | 从坡顶到该层底面的深度（m） |
| `name` | 可选，名称（用于图例） |
| `color` | 可选，绘图颜色 |

### 3. 设置水条件

```python
# 浸润线距坡顶深度 5m → 即 y = H - 5 处的水平线
s.set_water_table(depth=5)

# 移除水条件
s.set_water_table(depth=None)
```

孔隙水压力按静水压力计算：`u = γ_w × (y_water - y)`，γ_w = 9.81 kN/m³。

### 4. 设置外荷载

```python
from pyslopex import Udl, LineLoad

# 均布荷载：q=50kPa，从坡顶偏移 2m 处开始，宽 3m
s.set_udls(Udl(magnitude=50, offset=2, length=3))

# 均布荷载：无限宽（从坡顶开始）
s.set_udls(Udl(magnitude=20))

# 集中线荷载：P=30kN/m，距坡顶 1m
s.set_line_loads(LineLoad(magnitude=30, offset=1))
```

| Udl 参数 | 说明 |
|---|---|
| `magnitude` | 荷载强度（kPa） |
| `offset` | 距坡顶的水平偏移（m），默认 0 |
| `length` | 荷载宽度（m），None 表示无限宽 |

| LineLoad 参数 | 说明 |
|---|---|
| `magnitude` | 线荷载强度（kN/m） |
| `offset` | 距坡顶的水平偏移（m），默认 0 |

---

## 平面滑动法

假设滑动面为通过坡脚的平面，搜索不同倾角找到最小安全系数。

```python
result = s.analyse_planar(
    min_angle=5.0,       # 最小搜索角度（度）
    max_angle=None,      # 最大搜索角度，默认 = 坡角 - 0.1°
    num_angles=200,      # 搜索角度数
    num_slices=50        # 竖向条分数
)
```

### 返回结果（PlanarResult）

| 属性 | 说明 |
|---|---|
| `fos` | 最小安全系数 |
| `critical_angle` | 临界滑动面角度（度） |
| `failure_plane` | `[(x1,y1), (x2,y2)]` 滑动面两端坐标 |
| `sliding_wedge` | 滑动楔体顶点坐标 |
| `all_results` | 所有试算角度的 `(角度, Fs)` 列表 |

### 原理

对过坡脚、倾角为 θ 的平面：

$$Fs = \frac{c \cdot L + (W\cos\theta - U)\tan\phi}{W\sin\theta}$$

---

## 折线滑动法（传递系数法）

采用国标 GB 50330 传递系数法，适用于任意折线形滑动面。

### 方式一：手动指定滑动面

```python
# 滑动面控制点：从坡顶（右上）到坡脚（左下）排列
points = [(15, 8), (8, 5), (0, 0)]
result = s.analyse_polyline(failure_surface_points=points)
```

### 方式二：自动搜索最危险滑动面

```python
result = s.analyse_polyline(
    auto_search=True,
    num_internal_points=2,    # 内部折点数（1-3）
    grid_resolution=20,       # 搜索网格精度
    method='implicit'         # 'implicit' 或 'explicit'
)
```

### 参数

| 参数 | 说明 |
|---|---|
| `failure_surface_points` | 滑动面控制点列表，从入口到出口排列 |
| `method` | `'implicit'`（隐式，二分法求解）或 `'explicit'`（显式，直接计算） |
| `auto_search` | 是否自动搜索最危险面 |
| `num_internal_points` | 自动搜索的内部折点数 |
| `grid_resolution` | 网格搜索精度 |
| `entry_x_range` | 入口点搜索范围 `(x_min, x_max)` |

### 返回结果（PolylineResult）

| 属性 | 说明 |
|---|---|
| `fos` | 安全系数 |
| `failure_surface` | 滑动面控制点 `[(x,y), ...]` |
| `slice_data` | 各分块详细信息（角度、重量、力等） |
| `thrust_distribution` | 各分块边界的剩余推力（kN/m） |
| `thrust_at_toe` | 坡脚处剩余推力 |
| `method` | 使用的求解方法 |

### 原理

传递系数法核心公式：

$$P_i = F_s \cdot T_i + P_{i-1} \cdot \psi_{i-1} - R_i$$

其中：
- $T_i = W_i \sin\alpha_i$（下滑力）
- $R_i = c_i L_i + (W_i \cos\alpha_i - U_i) \tan\phi_i$（抗滑力）
- $\psi_{i-1} = \cos(\alpha_{i-1}-\alpha_i) - \frac{\sin(\alpha_{i-1}-\alpha_i) \tan\phi_i}{F_s}$（传递系数）

隐式法通过二分法迭代求 $F_s$，使坡脚推力 $P_n = 0$。中间分块推力小于 0 时置 0（国标做法）。

---

## 绘图

```python
# 边坡边界图（含土层、水位线、荷载）
fig = s.plot_boundary()

# 平面滑动法结果图
result = s.analyse_planar()
fig = s.plot_critical_planar(result)

# 折线滑动法结果图
fig = s.plot_polyline(failure_surface_points=points)

# 保存图片
fig.savefig("result.png", dpi=150)
fig.savefig("result.pdf")
```

---

## 完整示例

### 示例 1：基本分析

```python
from pyslopex import Slope, Material

s = Slope(height=10, angle=45)
s.set_materials(Material(unit_weight=20, friction_angle=30,
                         cohesion=10, depth_to_bottom=10))

# 平面滑动法
r_plan = s.analyse_planar()
print(f"平面法: Fs = {r_plan.fos:.3f}, θ = {r_plan.critical_angle:.1f}°")
# → 平面法: Fs = 1.513, θ = 33.1°

# 折线滑动法（自动搜索）
r_poly = s.analyse_polyline(auto_search=True, num_internal_points=2)
print(f"折线法: Fs = {r_poly.fos:.3f}")
# → 折线法: Fs = 1.258
```

### 示例 2：考虑水压力和荷载

```python
from pyslopex import Slope, Material, Udl

# 8m 高，坡比 1:1.5
s = Slope(height=8, length=12)
s.set_materials(Material(unit_weight=19, friction_angle=30,
                         cohesion=5, depth_to_bottom=12))

# 浸润线 y=3m（距坡顶 5m）
s.set_water_table(depth=5)

# 坡顶均布荷载 20kPa，宽 3m
s.set_udls(Udl(magnitude=20, offset=0, length=3))

# 自动搜索
r = s.analyse_polyline(auto_search=True, num_internal_points=2,
                        grid_resolution=20)
print(f"Fs = {r.fos:.3f}")
```

### 示例 3：指定滑动面

```python
s = Slope(height=10, angle=45)
s.set_materials(Material(unit_weight=20, friction_angle=30,
                         cohesion=10, depth_to_bottom=10))

# 手动指定 3 段折线滑动面
points = [(16, 10), (11, 7), (5, 3), (0, 0)]

# 隐式法
r_imp = s.analyse_polyline(failure_surface_points=points, method='implicit')
print(f"隐式法 Fs = {r_imp.fos:.3f}")

# 显式法
r_exp = s.analyse_polyline(failure_surface_points=points, method='explicit')
print(f"显式法 Fs = {r_exp.fos:.3f}")

# 查看各分块详情
for i, sd in enumerate(r_imp.slice_data):
    print(f"段{i+1}: α={sd.alpha_deg:.1f}°, W={sd.weight:.1f}kN/m, "
          f"U={sd.U:.1f}kN/m, T={sd.driving:.1f}, R={sd.resisting:.1f}")
```

### 示例 4：多层土 + 完整分析

```python
from pyslopex import Slope, Material, Udl, LineLoad

s = Slope(height=10, angle=45)

# 三层土
s.set_materials(
    Material(unit_weight=21, friction_angle=35, cohesion=18,
             depth_to_bottom=3, name="密实砂", color="#D2B48C"),
    Material(unit_weight=20, friction_angle=30, cohesion=12,
             depth_to_bottom=7, name="粉质粘土", color="#A0522D"),
    Material(unit_weight=19, friction_angle=25, cohesion=8,
             depth_to_bottom=12, name="软粘土", color="#8B7355"),
)

s.set_water_table(depth=3)
s.set_udls(Udl(magnitude=50, offset=1, length=2))
s.set_line_loads(LineLoad(magnitude=20, offset=5))

# 平面法
r_plan = s.analyse_planar()
print(f"平面法 Fs = {r_plan.f_plan:.3f}")

# 折线法（自动搜索）
r_poly = s.analyse_polyline(auto_search=True, num_internal_points=2,
                             grid_resolution=20)
print(f"折线法 Fs = {r_poly.fos:.3f}")
print(f"最危险面: {[(f'{p[0]:.1f}',f'{p[1]:.1f}') for p in r_poly.failure_surface]}")

# 绘图
s.plot_critical_planar(r_plan)
s.plot_polyline(result=r_poly)
```

---

## 项目结构

```
pyslopex/
├── __init__.py        # 包入口，导出 Slope, Material, Udl, LineLoad
├── slope.py           # Slope 边坡模型类
├── material.py        # Material 土层类
├── loads.py           # Udl / LineLoad 荷载类
├── planar.py          # 平面滑动法分析
├── polyline.py        # 折线滑动法（传递系数法）分析
├── plotting.py        # matplotlib 绘图
└── utils.py           # 几何工具函数
tests/
├── test_planar.py     # 平面法测试（13 项）
└── test_polyline.py   # 折线法测试（12 项）
examples/
└── example_usage.py   # 使用示例
```

## 运行测试

```bash
pip install pytest
pytest tests/ -v
```

## 算法说明

### 平面滑动法

过坡脚搜索不同倾角的平面，用竖向条分法计算每个试算平面的 Fs，取最小值。支持分层土、孔隙水压力、外荷载。

### 折线滑动法（传递系数法）

基于国标 GB 50330 。将滑体按折线顶点划分为竖直分块，自上而下逐块计算剩余推力，通过传递系数 ψ 将上方分块的推力传递到下方。

**隐式法**（推荐）：F_s 出现在传递系数中，用二分法迭代至坡脚推力为零。中间分块推力小于 0 时置零（符合国标做法）。

**显式法**：F_s 不出现在传递系数中，可一次计算得出。当传递系数 ψ < 0 或分块下滑力为负时自动回退为隐式法。
