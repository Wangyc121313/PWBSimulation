# AGENTS.md

## 项目背景

这个工作区是一个 PWB（Photonic Wire Bonding，光子引线键合）仿真项目。核心流程是通过 Python 调用 Lumerical/Ansys FDTD Python API（`lumapi`）来构建仿真结构、运行 FDTD 仿真，并对导出的数据、图片等结果进行分析。

项目按耦合场景组织：

- `PD-PWB-SMF/`：当前脚本化程度最高的目录。包含可复用 Python 代码 `pwb_core.py`、单次运行脚本、参数扫描脚本、`.fsp` 工程文件、notebook 和已生成结果。
- `LD-PWB-SMF/`：LD 到 PWB 到 SMF 的仿真资源，包含 notebook、`.fsp` 工程、LSF 脚本片段和扫描结果。
- `LNOI-PWB-SMF/`：LNOI 到 PWB 到 SMF 的仿真资源，以及参数扫描图片和 CSV 结果。
- `SMF-PWB-SMF/`：SMF 到 PWB 到 SMF 的仿真工程和导出结果。
- `Taper/`：用于构建 taper 的 Lumerical 脚本片段。
- `U-Turn/`：用途需先检查后再判断，不要直接假设其角色。
- 根目录下的 `data analysis.ipynb` 和 `results_analysis.ipynb` 主要用于数据分析和绘图。

大型二进制文件和生成结果是研究流程的一部分：

- `.fsp` 文件是 Lumerical FDTD 工程文件。
- `results/` 下的 `.h5`、`.png`、`.jpg`、`.txt`、`.csv` 文件通常是仿真输出或分析产物。
- 除非用户明确要求，不要删除、重新生成或批量改动这些结果文件。

## 运行前提

运行仿真脚本需要本机安装 Lumerical/Ansys，并能访问 FDTD Python API。当前统一通过根目录 `sim_config.py` 配置并引入 `lumapi`：

```python
from sim_config import add_lumerical_api_path

add_lumerical_api_path()
import lumapi
```

默认 API 路径为 `D:\Program Files\Lumerical\v241\api\python`，可通过环境变量 `LUMERICAL_API_PATH` 覆盖。

脚本还依赖常见科学计算包，例如 `numpy`、`matplotlib` 和 `pandas`。

FDTD 仿真通常耗时较长，并依赖商业软件许可、GPU 和本机资源。不要随意运行完整参数扫描；除非用户明确要求，优先做代码检查、轻量分析或后处理。

## 路径注意事项

路径配置集中在根目录 `sim_config.py`。新增或修改脚本时，应优先从该文件导入路径常量，而不是在各个 `.py` 文件中硬编码绝对路径。

常用配置：

- `PROJECT_ROOT`：项目根目录，默认等于 `sim_config.py` 所在目录，可通过 `SIM_PROJECT_ROOT` 覆盖。
- `LUMERICAL_API_PATH`：Lumerical Python API 路径，可通过 `LUMERICAL_API_PATH` 覆盖。
- `MATERIAL_DB`：材料库路径，默认是根目录 `database.mdf`，可通过 `SIM_MATERIAL_DB` 覆盖。
- `PD_DIR`、`LD_DIR`、`LNOI_DIR`：三个主要场景目录。
- `PD_RESULTS_DIR`、`LD_RESULTS_DIR`、`LNOI_RESULTS_DIR`：三个主要场景的结果目录。
- `PD_SMF_FSP`、`LD_BASE_FSP`、`LNOI_BASE_FSP`：常用基础 `.fsp` 工程文件。

当前工作区根目录是：

```text
D:/simulation/Simulation Project
```

如果发现历史脚本里还有类似 `D:/simulation/Simulation Project/simulation/...` 的旧路径，应优先改为 `sim_config.py` 中的路径常量，并注意额外的 `simulation/` 目录层级通常来自旧目录布局。

## 关键 Python 入口

`PD-PWB-SMF/` 中的主要脚本：

- `pwb_core.py`：共享参数类、PWB 路径生成、结构构建、FDTD 设置、结果读取和绘图辅助函数。
- `run_single.py`：构建 Section 1，保存 `.fsp`，运行 FDTD，绘制电场图片，并打印 `T_total`。
- `sweep_r_R.py`：对波导半径 `r` 和弯曲半径 `R` 做二维扫描；输出 `results/r_R_scan/T_total_sweep_results.txt` 和对应电场图。
- `sweep_h_R.py`：对总高度 `h` 和弯曲半径 `R` 做二维扫描；输出 `results/h_R_scan/T_total_sweep_results.txt`。
- `analyze_results.py`：读取已有扫描结果并绘制 loss 热力图。它不启动 FDTD，是最适合轻量运行的脚本。
- `test_setup.py`：只构建并保存 FDTD 结构，不执行完整仿真；但仍然依赖 `lumapi`。

`LD-PWB-SMF/` 中的主要脚本：

- `pwb_core.py`：从 notebook 抽出的 LD-PWB-SMF 核心逻辑，包含中心线生成、椭圆/圆形截面半径函数、FDTD 结构生成、仿真设置、结果读取和绘图。
- `test_setup.py`：只生成并保存 LD-PWB-SMF 结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `parameter.md`：说明 LD-PWB-SMF 的 `r/l1/l2/r_11/r_12/r_2/r_3` 等参数。
- `simulation-LD-PWB-SMF.ipynb`：保留原始探索、历史扫描和人工绘图记录。

`LNOI-PWB-SMF/` 中的主要脚本：

- `pwb_core.py`：从 notebook 抽出的 LNOI-PWB-SMF 核心逻辑，包含 pyramid/straight 结构生成、仿真设置、结果读取和绘图。
- `test_setup.py`：只生成并保存 LNOI-PWB-SMF 结构，不运行 FDTD。
- `run_single.py`：生成结构、运行一次仿真并绘图。
- `sweep_h1_h2_w1_w2.py`：执行 `h1/h2/w1/w2` 四参数扫描，输出 CSV 和场图。
- `parameter.md`：说明 LNOI-PWB-SMF 的 `h1/h2/w1/w2/l` 等参数。
- `simulation-LNOI-PWB-SMF.ipynb`：保留原始探索和历史扫描记录。

## 开发约定

- 修改应尽量限制在用户提到的具体场景目录和脚本内。
- 若要增强可复现性，优先使用 `sim_config.py` 参数化路径；不要在未经明确同意的情况下重组目录结构。
- 保持物理单位一致。代码中的几何参数大多使用 SI 单位，通常以微米数值乘以 `1e-6` 表示。
- `PD-PWB-SMF` 中 `_1/_2` 与 `_3` 的弯曲半径含义不同：`_1/_2` 的 `R` 是固定 90 度圆弧半径；`_3` 使用分段 Bezier 复杂中心线，不存在全局固定弯曲半径，应关注沿中心线变化的局部曲率半径 `rho(s)=1/kappa(s)`，尤其是最小曲率半径 `min(rho)`。
- 在 `_3` 当前实现中，`R` 不直接控制中心线曲率，只用于 `generate_radius_profile_3` 中估算输出 taper 区域长度；复杂路径形状主要由 `L`、`h`、`l1`、`bend_lift`、`arch_position`、`bend_shape`、`drop_shape` 和 `complex_segments` 控制。
- 谨慎处理 `T_total` 的符号和 loss 换算。现有分析使用 `loss = -10 * log10(abs(-1.0 * T_total))`；修改公式前需确认物理约定。
- 新增分析代码时，输出文件应放入对应场景目录的 `results/` 子目录。
- 对 notebook 不做大规模重写。若逻辑需要复用，优先抽取到 `.py` 文件。
- 生成图片的文件名应清楚包含扫描参数值。

## 验证方式

轻量验证：

- 仅运行后处理：`python PD-PWB-SMF/analyze_results.py`
- 对 Python 文件做语法检查：`python -m py_compile sim_config.py PD-PWB-SMF/pwb_core.py`

仿真相关验证仅在用户要求且本机 Lumerical 可用时执行：

- 只构建结构、不完整运行：`python PD-PWB-SMF/test_setup.py`
- 单次 FDTD 仿真：`python PD-PWB-SMF/run_single.py`
- 参数扫描：`python PD-PWB-SMF/sweep_r_R.py` 或 `python PD-PWB-SMF/sweep_h_R.py`

LD/LNOI 验证入口：

- LD 只构建结构：`python LD-PWB-SMF/test_setup.py`
- LD 单次仿真：`python LD-PWB-SMF/run_single.py`
- LNOI 只构建结构：`python LNOI-PWB-SMF/test_setup.py`
- LNOI 单次仿真：`python LNOI-PWB-SMF/run_single.py`
- LNOI 参数扫描：`python LNOI-PWB-SMF/sweep_h1_h2_w1_w2.py`
