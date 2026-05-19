# LD-PWB-SMF 参数说明

本目录用于 LD -> PWB -> SMF 的耦合仿真。当前已将 notebook 中的主要逻辑整理到 `pwb_core.py`，并提供：

- `test_setup.py`：只生成并保存结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `simulation-LD-PWB-SMF.ipynb`：保留原始探索和历史扫描记录。

## 0. 路径配置

脚本中的工程文件、材料库、数据源和结果目录路径统一从根目录 `sim_config.py` 获取。常用项包括：

- `LD_DIR`：当前场景目录。
- `LD_RESULTS_DIR`：结果输出目录。
- `LD_BASE_FSP`：基础 `LD.fsp` 工程文件。
- `LD_SOURCE_DATASET`：输入场数据 `Section2 output.mat`。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

如需迁移项目或调整 Lumerical 安装位置，优先修改环境变量 `SIM_PROJECT_ROOT`、`LUMERICAL_API_PATH`、`SIM_MATERIAL_DB`，或调整 `sim_config.py`，不要在各个脚本中分散修改绝对路径。

## 1. 中心线生成

中心线由 `generate_pwb_path(params)` 生成，位于 `pwb_core.py`。

路径是二维 `x-z` 平面结构，`y = 0`，共 7 段：

```text
1. 输入 taper
2. 输入水平直段
3. 向上弯曲
4. 中间弯曲
5. 向下弯曲
6. 输出水平直段
7. 输出 taper
```

核心几何参数：

```python
r              # 弯曲圆弧半径
l1             # 输入 taper 长度
l2             # 输出 taper 长度
total_length   # 总目标水平长度，默认 170 um
vertical_offset # 弯曲上拱/偏移高度，默认 20 um
curve_points   # 每段中心线采样点数
```

`r` 在 LD 结构中是实际圆弧构造半径。代码中通过：

```python
theta = arcsin(vertical_offset / r)
```

计算弯曲角度。因此必须满足：

```text
vertical_offset <= r
```

否则无法生成实数角度。

## 2. 水平直段长度

输入和输出之间的中间水平补偿长度为：

```python
l = (total_length - (l1 + l2)) / 2
```

因此需要满足：

```text
l1 + l2 <= total_length
```

如果增大 `l1` 或 `l2`，中间两个水平直段会变短。若二者之和超过 `total_length`，路径生成会报错。

## 3. 截面半径模型

LD-PWB-SMF 中 PWB 输入 taper 使用椭圆截面，当前有两套半径函数：

```python
get_radius_at_position_1(t, params)
get_radius_at_position_2(t, params)
```

相关参数：

```python
r_11  # 输入 taper 的第一轴起始半径
r_12  # 输入 taper 的第二轴起始半径
r_2   # 中间 PWB 主体半径
r_3   # 输出 taper 末端半径，通常用于接 SMF
```

半径变化逻辑：

```text
路径前 1/7: 输入 taper
    r_11 -> r_2
    r_12 -> r_2

路径中间 5/7: 主体 bend
    r_2

路径最后 1/7: 输出 taper
    r_2 -> r_3
```

注意：

```text
r, r_11, r_12, r_2, r_3 含义不同
```

- `r`：路径弯曲圆弧半径。
- `r_11/r_12/r_2/r_3`：PWB 横截面半径参数。

## 4. 仿真设置

`setup_fdtd_simulation(fdtd, params, path)` 负责添加 FDTD 区域、光源、监视器和 mode expansion。

当前整理版本保留 notebook 的主要设置：

- FDTD x 范围：`140 um -> 270 um`
- z 范围：`-7 um -> 7 um`
- source 位置：`x = 150 um`
- input monitor：`x = 160 um`
- output monitor：`x = 260 um`
- transmission monitor：`x = 205 um`, `x span = 130 um`

## 5. 导入光源

原 notebook 使用：

```python
fdtd.addimportedsource()
fdtd.importdataset("Section2 output.mat")
```

当前整理版本默认仍尝试使用：

```python
params.use_imported_source = True
params.source_dataset = LD-PWB-SMF/Section2 output.mat
```

如果 `Section2 output.mat` 不存在，代码会自动回退到普通 mode source，避免 setup 脚本直接失败。若要强制使用普通 mode source，可设置：

```python
params.use_imported_source = False
```

## 6. 常用调参方向

降低弯曲损耗：

- 增大 `r`。
- 减小 `vertical_offset`。
- 适当增大 `total_length`，给结构更多水平距离。

优化输入 taper：

- 调整 `l1`。
- 调整 `r_11`、`r_12`、`r_2`。

优化输出到 SMF：

- 调整 `l2`。
- 调整 `r_3`。

提高结构光滑度：

- 增大 `curve_points`。
- 但注意 FDTD 中对象数量会随之增加。

## 7. 运行方式

只生成结构：

```powershell
python LD-PWB-SMF\test_setup.py
```

单次仿真：

```powershell
python LD-PWB-SMF\run_single.py
```

完整参数扫描仍保留在 notebook 中。后续如需复用，可按 PD-PWB-SMF 的方式拆分成独立 `sweep_*.py` 脚本。
