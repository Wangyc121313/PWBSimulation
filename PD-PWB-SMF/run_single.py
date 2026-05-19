"""
run_single.py — 单次完整仿真
搭建 Section 1 结构，运行仿真，可视化电场分布并打印 T_total。
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import PD_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_1,
    setup_fdtd_simulation_1,
    visualize_and_save_results_1,
)

SAVE_PATH = PD_DIR / "test" / "test.fsp"

params = PWBParameters()
fdtd = lumapi.FDTD()

generate_pwb_structure_1(fdtd, params)
setup_fdtd_simulation_1(fdtd, params)
SAVE_PATH.parent.mkdir(parents=True, exist_ok=True)
fdtd.save(str(SAVE_PATH))

fdtd.run()

T_total = visualize_and_save_results_1(fdtd, params)
print(f"T_total = {T_total}")

fdtd.close()
