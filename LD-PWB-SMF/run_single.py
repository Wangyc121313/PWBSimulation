"""Run one LD-PWB-SMF simulation and plot the monitor fields."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LD_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PWBParameters,
    create_pwb_structure_in_fdtd,
    setup_fdtd_simulation,
    visualize_and_save_results,
)


SAVE_PATH = LD_DIR / "LD_with_PWB_temp.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()
path = create_pwb_structure_in_fdtd(fdtd, params)
setup_fdtd_simulation(fdtd, params, path)
fdtd.save(str(SAVE_PATH))
fdtd.run()

T_forward = visualize_and_save_results(fdtd, params)
print(f"T_forward = {T_forward}")
fdtd.close()
