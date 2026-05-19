# LNOI-PWB-SMF 参数说明

本目录用于 LNOI -> PWB -> SMF 的耦合仿真。当前已将 notebook 中的主要逻辑整理到 `pwb_core.py`，并提供：

- `test_setup.py`：只生成并保存结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `sweep_h1_h2_w1_w2.py`：执行 `h1/h2/w1/w2` 四参数扫描。
- `simulation-LNOI-PWB-SMF.ipynb`：保留原始探索和历史扫描记录。

## 0. 路径配置

脚本中的工程文件、材料库和结果目录路径统一从根目录 `sim_config.py` 获取。常用项包括：

- `LNOI_DIR`：当前场景目录。
- `LNOI_RESULTS_DIR`：结果输出目录。
- `LNOI_BASE_FSP`：基础 `LNOI.fsp` 工程文件。
- `MATERIAL_DB`：材料库文件。
- `add_lumerical_api_path()`：将 Lumerical Python API 路径加入 `sys.path`。

如需迁移项目或调整 Lumerical 安装位置，优先修改环境变量 `SIM_PROJECT_ROOT`、`LUMERICAL_API_PATH`、`SIM_MATERIAL_DB`，或调整 `sim_config.py`，不要在各个脚本中分散修改绝对路径。

## 1. 结构生成

核心函数：

```python
create_pwb_structure_in_fdtd(fdtd, params)
```

该函数会加载：

```text
LNOI-PWB-SMF/LNOI.fsp
database.mdf
```

然后添加两部分 PWB 结构：

```text
1. PWB taper：addpyramid()
2. PWB-Straight：addrect()
```

## 2. PWB taper 参数

主要参数位于 `PWBParameters`：

```python
h1 = 1.3e-6
h2 = 1.5e-6
w1 = 3.5e-6
w2 = 3.0e-6
l  = 80e-6
```

在 `addpyramid()` 中的对应关系：

```python
x = l / 2
x span bottom = 2 * h1
x span top    = 2 * h2
y span bottom = w1
y span top    = w2
z span         = l
rotation 1     = 90 around y-axis
```

因此：

- `l`：taper 长度。
- `h1`：taper bottom 的半厚度参数，对应 `x span bottom = 2*h1`。
- `h2`：taper top 的半厚度参数，对应 `x span top = 2*h2`。
- `w1`：taper bottom 的宽度。
- `w2`：taper top 的宽度。

注意这里的 `h1/h2` 名称来自原 notebook。由于结构经过 `rotation 1 = 90`，它们最终体现为截面尺寸参数，而不是简单的全局 z 高度。

## 3. PWB-Straight 参数

整理版本新增了两个命名参数：

```python
straight_length = 30e-6
straight_offset = 15e-6
```

对应原 notebook 中：

```python
x = 15e-6 + l
x span = 30e-6
y span = w2
z = h2 / 2
z span = h2
```

即：

- `straight_offset`：直波导中心相对 taper 末端的偏移。
- `straight_length`：直波导长度。
- `w2/h2`：直波导截面尺寸。

## 4. 仿真设置

核心函数：

```python
setup_fdtd_simulation(fdtd, params)
```

当前保留 notebook 中的主要监视器设置：

- source：`x = -5 um`
- monitor_0：`x = 0`
- monitor_60：`x = 60 um`
- monitor_80：`x = 80 um`
- monitor_90：`x = 90 um`
- mode expansion：绑定 `monitor_90`
- transmission_monitor：`2D Y-normal`，覆盖 `x = -10 um -> 95 um`

FDTD 区域：

```text
x: -10 um -> 95 um
z: -2 um -> h2 + 2 um
y: -w2/2 - 2 um -> w2/2 + 2 um
```

## 5. 结果指标

`get_data(fdtd, params)` 返回：

```python
source_E
transmission_E
monitor_0_E
monitor_60_E
monitor_80_E
monitor_90_E
T_forward
```

其中 `T_forward` 来自：

```python
mode_data = fdtd.getresult("mode_expansion", "expansion for output")
T_forward = mode_data["T_forward"][0][0]
```

## 6. 参数扫描

`sweep_h1_h2_w1_w2.py` 扫描：

```python
h1_values = np.arange(0.7, 1.7, 0.2) * 1e-6
h2_values = np.arange(0.7, 1.7, 0.2) * 1e-6
w1_values = np.arange(1.5, 4.0, 0.5) * 1e-6
w2_values = np.arange(1.5, 4.0, 0.5) * 1e-6
```

输出位置：

```text
results/h1_h2_w1_w2_scan/T_forward_results.csv
results/h1_h2_w1_w2_scan/*.jpg
```

该扫描组合数量较多，会启动大量 FDTD 运行。运行前建议先缩小参数范围做测试。

## 7. 运行方式

只生成结构：

```powershell
python LNOI-PWB-SMF\test_setup.py
```

单次仿真：

```powershell
python LNOI-PWB-SMF\run_single.py
```

四参数扫描：

```powershell
python LNOI-PWB-SMF\sweep_h1_h2_w1_w2.py
```

运行扫描前请确认 Lumerical 可用，并预留足够时间。
