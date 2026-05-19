"""
sweep_h_R.py — Section 2 的 h / R 参数化扫描
对 PWB 总高度 h 和弯曲半径 R 进行二维网格扫描，
保存每组参数的 T_total 到 results/h_R_scan/。
"""

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import PD_DIR, PD_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

import numpy as np

from pwb_core import (
    PWBParameters,
    generate_pwb_path_2,
    generate_pwb_structure_2,
    setup_fdtd_simulation_2,
    get_data_2,
)

RESULTS_DIR = PD_RESULTS_DIR / "h_R_scan"
SAVE_PATH = PD_DIR / "temp.fsp"

os.makedirs(RESULTS_DIR, exist_ok=True)

# 定义参数范围
R_values = np.arange(100, 120, 10) * 1e-6
h_values = np.arange(120, 180, 10) * 1e-6

all_results = []

for i, R in enumerate(R_values):
    for j, h in enumerate(h_values):
        total = len(R_values) * len(h_values)
        idx   = i * len(h_values) + j + 1
        print(f"[{idx}/{total}] R={R*1e6:.1f} μm, h={h*1e6:.1f} μm")

        params = PWBParameters()
        params.R = R
        params.h = h

        fdtd = lumapi.FDTD()
        generate_pwb_path_2(params)   # 预计算路径（内部调用）
        generate_pwb_structure_2(fdtd, params)
        setup_fdtd_simulation_2(fdtd, params)
        fdtd.save(str(SAVE_PATH))
        fdtd.run()

        results = get_data_2(fdtd, params)
        T_total = results['T_total']
        all_results.append({'R': R * 1e6, 'h': h * 1e6, 'T_total': T_total})
        print(f"  T_total = {T_total}")

        fdtd.close()

# 保存汇总结果
results_file = os.path.join(RESULTS_DIR, "T_total_sweep_results.txt")
with open(results_file, 'w', encoding='utf-8') as f:
    f.write("R(μm)\th(μm)\tT_total\n")
    for res in all_results:
        f.write(f"{res['R']:.1f}\t{res['h']:.1f}\t{res['T_total']:.6f}\n")

print(f"\n扫描完成，结果已保存至: {results_file}")
