# PD-PWB-SMF 参数说明

本文说明 `pwb_core.py` 中 PD-PWB-SMF 几何生成相关参数，重点解释 `_1/_2` 的圆弧半径 `R` 与 `_3` 复杂中心线中的“弯曲半径”概念差异。

## 0. 路径配置

脚本中的工程文件、材料库和结果目录路径统一从根目录 `sim_config.py` 获取。常用项包括：

- `PD_DIR`：当前场景目录。
- `PD_RESULTS_DIR`：结果输出目录。
- `PD_SMF_FSP`：基础 `SMF.fsp` 工程文件。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

如需迁移项目或调整 Lumerical 安装位置，优先修改环境变量 `SIM_PROJECT_ROOT`、`LUMERICAL_API_PATH`、`SIM_MATERIAL_DB`，或调整 `sim_config.py`，不要在各个脚本中分散修改绝对路径。

## 1. `_1/_2` 的简单圆弧逻辑

`generate_pwb_path_1` 和 `generate_pwb_path_2` 使用固定 90 度圆弧：

```python
center_x = l1 + l
center_z = -R
x2 = center_x + R * np.sin(t)
z2 = center_z + R * np.cos(t)
```

因此在 `_1/_2` 中，`R` 是非常直观的几何参数：

- `R` 就是 90 度 bend 的圆弧半径。
- 圆弧段的局部曲率半径处处等于 `R`。
- 增大 `R` 会降低弯曲曲率，一般有利于降低弯曲损耗，但会占用更大仿真空间。

这类结构的弯曲半径是一个固定值。

## 2. `_3` 的复杂中心线逻辑

`generate_pwb_path_3` 使用二维 `x-z` 平面中心线，`y = 0`。当前结构由三部分组成：

```text
水平 taper/直段 -> 上拱 Bezier -> 下弯 Bezier
```

水平直段：

```text
x: 0 -> l1
z: 0
```

弯曲段使用两段 cubic Bezier。关键点如下：

```python
dx = L - l1

start = [l1, 0, 0]
peak  = [l1 + arch_position * dx, 0, bend_lift * h]
end   = [L, 0, -h]
```

第一段 Bezier 从 `start` 平滑上拱到 `peak`。第二段 Bezier 从 `peak` 平滑下弯到 `end`。

当前实现中：

- 入口切线保持水平，避免水平 taper 到 bend 之间出现突兀折角。
- `peak` 是实际拱顶，拱高为 `bend_lift * h`。
- 末端切线趋向竖直方向，便于接近 taper2 的方向。

## 3. `_3` 中的“弯曲半径”是什么

在 `_3` 中，中心线不是圆弧，而是 Bezier 曲线。因此不存在一个像 `_1/_2` 那样全局固定的圆弧半径。

更准确的概念是**局部曲率半径**：

```text
局部曲率半径 rho(s) = 1 / kappa(s)
```

其中 `s` 是沿中心线的位置，`kappa(s)` 是该点曲率。对于二维 `x-z` 曲线，曲率可写为：

```text
kappa = |x' z'' - z' x''| / (x'^2 + z'^2)^(3/2)
rho = 1 / kappa
```

所以 `_3` 中的弯曲半径是随位置变化的：

- 上拱段有一组局部曲率半径。
- 下弯段有另一组局部曲率半径。
- 曲率最大的地方对应最小局部曲率半径。
- 对光学损耗而言，通常最需要关注的是 `min(rho)`，即最小曲率半径，而不是某一个输入参数 `R`。

结论：

```text
_1/_2: R 是固定圆弧半径。
_3: 没有固定弯曲半径，只有沿路径变化的局部曲率半径 rho(s)。
```

## 4. `_3` 中 `R` 当前的作用

当前 `_3` 算法中，`R` 不直接控制中心线的弯曲半径。

它只在半径函数中用于估算末端输出 taper 的长度：

```python
output_len = max(h - R, 0.0)
```

也就是说，`R` 在 `_3` 中目前更像一个“输出 taper 长度参考参数”，而不是几何曲率半径。真正控制 `_3` 中心线形状和局部曲率半径的是：

- `L`
- `h`
- `l1`
- `bend_shape`
- `bend_lift`
- `arch_position`
- `drop_shape`

如果需要把 `_3` 也改造成“指定最小弯曲半径”的形式，需要额外加入曲率约束或曲率优化逻辑。目前实现还没有自动保证 `min(rho) >= R_target`。

## 5. `_3` 关键参数

### `L`

总水平长度。复杂路径终点为：

```text
x = L
z = -h
```

增大 `L` 会给弯曲段更多水平距离，通常会增大局部曲率半径，使弯曲更缓。

### `h`

总高度差。复杂路径终点为：

```text
z = -h
```

增大 `h` 会增加下弯幅度，若 `L` 不变，通常会减小局部曲率半径。

### `l1`

入口水平直段长度。

- 增大 `l1`：入口直段更长，但弯曲段可用水平距离 `L - l1` 更短，可能导致曲率增大。
- 减小 `l1`：弯曲段可用水平距离更长，通常更容易获得缓弯。

### `bend_shape`

控制入口上拱段的形状，主要影响从水平 taper 到拱顶之间的过渡。

- 增大：入口水平趋势保持更久，上拱更舒展。
- 减小：更早开始明显上拱。

当前代码会把该值限制在 `0.05 ~ 0.95`。

### `bend_lift`

控制实际上拱高度：

```text
arch_height = bend_lift * h
```

- 增大：上拱更明显，更接近“先拱起再下弯”的结构。
- 减小：上拱减弱，路径更直接地下弯。
- 设为 `0`：取消上拱。

### `arch_position`

控制拱顶在弯曲段水平长度中的位置：

```text
peak_x = l1 + arch_position * (L - l1)
```

- 较小：拱顶更靠近入口，前段上拱更快。
- 较大：拱顶更靠近出口，前段更长、更平缓。

当前代码会把该值限制在 `0.05 ~ 0.90`。

### `drop_shape`

控制从拱顶下弯到终点的形状。

- 增大：下弯段前半部分更舒展，末端接近竖直方向的过渡更明显。
- 减小：下弯更早发生，局部曲率可能更集中。

当前代码会把该值限制在 `0.05 ~ 0.95`。

### `complex_segments`

Bezier 弯曲段采样点数。

- 增大：FDTD 中分段圆柱更多，几何更平滑，但对象数量和仿真设置开销增加。
- 减小：生成更快，但圆柱拼接更粗。

建议先用 `200 ~ 500` 调形状，最终仿真前再按需要提高。

## 6. 半径函数

`generate_radius_profile_3` 的“半径”指 PWB 横截面半径，不是弯曲半径。

它按中心线归一化弧长 `s` 变化：

```text
s = 0   路径起点
s = 1   路径终点
```

当前分三段：

```text
入口 taper: r1 -> r
中间 bend: r
出口 taper: r -> r2
```

相关参数：

- `r1`：入口半径。
- `r`：bend 主体半径。
- `r2`：出口半径。

注意区分：

```text
r / r1 / r2: 波导横截面半径
R 或 rho: 路径弯曲半径/曲率半径
```

## 7. 如何估算 `_3` 的局部曲率半径

可以用中心线点列做数值估算。示例：

```python
import numpy as np
from pwb_core import PWBParameters, generate_pwb_path_3

params = PWBParameters()
path = generate_pwb_path_3(params)

x = path[:, 0]
z = path[:, 2]

dx = np.gradient(x)
dz = np.gradient(z)
ddx = np.gradient(dx)
ddz = np.gradient(dz)

kappa = np.abs(dx * ddz - dz * ddx) / np.clip((dx**2 + dz**2) ** 1.5, 1e-30, None)
rho = 1 / np.clip(kappa, 1e-30, None)

print("min curvature radius (um):", np.nanmin(rho) * 1e6)
```

这个 `min curvature radius` 更接近 `_3` 中需要关注的“最小弯曲半径”指标。

## 8. 调参建议

如果希望整体弯曲半径变大、损耗降低：

- 增大 `L`。
- 减小 `l1`，给 bend 留更多水平空间。
- 减小 `h` 或 `bend_lift`。
- 增大 `arch_position`，让上拱更平缓。
- 调整 `drop_shape`，避免曲率集中在某一小段。

如果希望上拱更明显：

- 增大 `bend_lift`。
- 将 `arch_position` 调到合适位置，例如 `0.3 ~ 0.6`。

如果希望入口更平滑：

- 增大 `bend_shape`。
- 提高 `complex_segments`。

如果希望末端接 taper2 更自然：

- 调整 `drop_shape`。
- 观察末端附近曲率是否过大。

## 9. 运行方式

可以在 `test_setup_complex.py` 中设置参数，例如：

```python
params = PWBParameters()
params.L = 250e-6
params.h = 150e-6
params.l1 = 100e-6

params.bend_lift = 0.12
params.arch_position = 0.60
params.bend_shape = 0.45
params.drop_shape = 0.45
params.complex_segments = 400
```

然后运行：

```powershell
python PD-PWB-SMF\test_setup_complex.py
```

