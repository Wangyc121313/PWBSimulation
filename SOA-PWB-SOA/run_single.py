"""Run one SOA-PWB-SOA simulation and plot the results."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import SOA_DIR, add_lumerical_api_path

# Ensure local pwb_core.py is found before pwb_core.py in other directories
# (sim_config inserts PD-PWB-SMF at sys.path[0], so we must insert after it)
sys.path.insert(0, str(Path(__file__).resolve().parent))

add_lumerical_api_path()
import lumapi

from pwb_core import (
    SOAPWBParams,
    create_pwb_structure_in_fdtd,
    setup_fdtd_simulation,
    visualize_and_save_results,
)

SAVE_PATH = SOA_DIR / "SOA_PWB_SOA_temp.fsp"

params = SOAPWBParams()
fdtd = lumapi.FDTD()

try:
    path = create_pwb_structure_in_fdtd(fdtd, params)
    setup_fdtd_simulation(fdtd, params, path)
    fdtd.save(str(SAVE_PATH))
    print(f"Structure saved. Total PWB length: {params.total_length * 1e6:.1f} um")
    print("Running FDTD simulation...")
    fdtd.run()
    print("Simulation complete. Extracting results...")
    T_forward = visualize_and_save_results(fdtd, params)
    print(f"T_forward = {T_forward}")
finally:
    fdtd.close()
