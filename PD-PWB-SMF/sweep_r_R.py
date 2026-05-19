"""
sweep_r_R.py — Section 1 的 r / R 参数化扫描
对纤芯半径 r 和弯曲半径 R 进行二维网格扫描，
保存每组参数的 T_total 及电场分布图到 results/r_R_scan/。
"""

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import PD_DIR, PD_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

import numpy as np
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from pwb_core import (
    PWBParameters,
    generate_pwb_path_1,
    generate_pwb_structure_1,
    setup_fdtd_simulation_1,
    get_data_1,
)

RESULTS_DIR = PD_RESULTS_DIR / "r_R_scan"
SAVE_PATH = PD_DIR / "test" / "test.fsp"

os.makedirs(RESULTS_DIR, exist_ok=True)
SAVE_PATH.parent.mkdir(parents=True, exist_ok=True)

# 定义参数范围
r_values = np.arange(0.4, 2.4, 0.2) * 1e-6
R_values = np.arange(20, 120, 10) * 1e-6

all_results = []

for i, r in enumerate(r_values):
    for j, R in enumerate(R_values):
        total = len(r_values) * len(R_values)
        idx   = i * len(R_values) + j + 1
        print(f"[{idx}/{total}] r={r*1e6:.2f} μm, R={R*1e6:.1f} μm")

        params = PWBParameters()
        params.r = r
        params.R = R

        fdtd = lumapi.FDTD()
        generate_pwb_path_1(params)   # 预计算路径（内部调用）
        generate_pwb_structure_1(fdtd, params)
        setup_fdtd_simulation_1(fdtd, params)
        fdtd.save(str(SAVE_PATH))
        fdtd.run()

        results  = get_data_1(fdtd, params)
        T_total  = results['T_total']
        all_results.append({'r': r * 1e6, 'R': R * 1e6, 'T_total': T_total})

        # 保存电场分布图
        try:
            plt.rcParams['font.family'] = 'Times New Roman'
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12
            fig, ax = plt.subplots(figsize=(6, 5))

            trans_E = results['transmission_E']
            E_field = trans_E['E']
            Ex = E_field[:, 0, :, 0, 0]
            Ey = E_field[:, 0, :, 0, 1]
            Ez = E_field[:, 0, :, 0, 2]
            E_intensity = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2

            x_trans = np.squeeze(trans_E['x']) * 1e6
            z_trans = np.squeeze(trans_E['z']) * 1e6
            X_trans, Z_trans = np.meshgrid(x_trans, z_trans, indexing='ij')

            norm = LogNorm(vmin=1e-5, vmax=np.nanmax(E_intensity))
            pcm  = ax.pcolormesh(X_trans, Z_trans, E_intensity, cmap="inferno", shading='auto', norm=norm)
            fig.colorbar(pcm, ax=ax, pad=0.02, format=ticker.LogFormatterSciNotation())
            ax.set_xlabel('X (μm)')
            ax.set_ylabel('Z (μm)')
            ax.set_title('Transmission Monitor')
            plt.tight_layout()

            fig_path = os.path.join(RESULTS_DIR, f"r_{r*1e6:.2f}_R_{R*1e6:.1f}_field_plot.png")
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  T_total = {T_total:.4f}  图片: {fig_path}")
        except Exception as e:
            print(f"  绘图时出错: {e}")
            print(f"  T_total = {T_total:.4f}")

        fdtd.close()

# 保存汇总结果
results_file = os.path.join(RESULTS_DIR, "T_total_sweep_results.txt")
with open(results_file, 'w', encoding='utf-8') as f:
    f.write("r(μm)\tR(μm)\tT_total\n")
    for res in all_results:
        f.write(f"{res['r']:.2f}\t{res['R']:.1f}\t{res['T_total']:.6f}\n")

print(f"\n扫描完成，结果已保存至: {results_file}")
