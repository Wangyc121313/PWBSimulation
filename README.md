# PWB 仿真项目 (Photonic Wire Bonding Simulation)

基于 Lumerical/Ansys FDTD Python API (`lumapi`) 的光子引线键合仿真平台，用于评估不同光电芯片间通过 PWB 聚合物波导互连的耦合效率。

## 仿真场景

| 场景 | 路径 | 描述 |
|------|------|------|
| PD-PWB-SMF | 光电探测器 → PWB → 单模光纤 | 含 Bezier 复杂中心线、参数扫描 |
| LD-PWB-SMF | 激光器 → PWB → 单模光纤 | 含导入光源、椭圆截面 taper |
| LNOI-PWB-SMF | 铌酸锂薄膜 → PWB → 单模光纤 | pyramid taper + 直波导 |
| SOA-PWB-SOA | 片上 SOA → PWB → 外置 SOA | 直波导 + 两端独立 taper + FDE 模式分析 |

## 环境要求

- **Lumerical FDTD Solutions** (安装于 `D:\Program Files\Lumerical\`)
- **Python 3.9.9** (使用 Lumerical 自带的 Python)
- Python 包：`numpy`, `matplotlib` (Lumerical Python 已自带)
- **MODE Solutions** (可选，用于 FDE 模式分析)

## 快速开始

```powershell
# 设置 Lumerical Python 路径
$env:PY = "D:\Program Files\Lumerical\python\python.exe"

# 语法检查（不需要 Lumerical）
python -m py_compile SOA-PWB-SOA\pwb_core.py

# 仅构建结构（需要 Lumerical）
& $env:PY SOA-PWB-SOA\test_setup.py

# 运行单次 FDTD 仿真（需要 Lumerical，耗时较长）
& $env:PY SOA-PWB-SOA\run_single.py
```

## 项目结构

```
PWBSimulation/
├── config/
│   ├── sim_config.py        # 统一路径配置
│   └── database.mdf         # 材料库
├── PD-PWB-SMF/              # PD → PWB → SMF 场景
│   ├── pwb_core.py          # 核心模块
│   ├── scripts/             # 仿真脚本（单次 + 参数扫描）
│   ├── tests/               # 测试（仅构建结构）
│   ├── analysis/            # 后处理分析
│   └── results/             # 仿真输出
├── LD-PWB-SMF/              # LD → PWB → SMF 场景
├── LNOI-PWB-SMF/            # LNOI → PWB → SMF 场景
├── SOA-PWB-SOA/             # SOA → PWB → SOA 场景
├── notebooks/               # Jupyter notebooks（数据分析）
├── docs/
│   └── superpowers/
│       ├── specs/           # 设计文档
│       └── plans/           # 实现计划
├── AGENTS.md                # AI Agent 使用说明
├── CLAUDE.md                # Claude Code 使用说明
└── README.md                # 本文件
```

## 开发约定

- **单位**：几何参数统一用 SI（米），以 `X * 1e-6` 表示微米
- **路径**：所有路径从 `config/sim_config.py` 导入，不硬编码
- **FDTD 代价**：优先做代码检查，未经明确要求不运行完整仿真或参数扫描
- **结果输出**：放入对应场景的 `results/` 子目录

## 参考

- 设计文档：`docs/superpowers/specs/`
- 实现计划：`docs/superpowers/plans/`
- 物理参数说明：各场景 `parameter.md`
