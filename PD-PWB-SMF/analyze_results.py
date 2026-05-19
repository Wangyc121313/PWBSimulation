"""Post-process existing PD-PWB-SMF sweep results without launching FDTD."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import PD_RESULTS_DIR
from pwb_core import plot_Ttotal_loss_heatmap_rR


data_file_rR = PD_RESULTS_DIR / "r_R_scan" / "T_total_sweep_results.txt"
loss_matrix, r_vals, R_vals = plot_Ttotal_loss_heatmap_rR(data_file_rR)
print("r values (um):", r_vals)
print("R values (um):", R_vals)

# To analyze h/R sweep results, use:
# data_file_hR = PD_RESULTS_DIR / "h_R_scan" / "T_total_sweep_results.txt"
# loss_matrix_hR, R_vals_hR, h_vals_hR = plot_Ttotal_loss_heatmap_rR(data_file_hR)
