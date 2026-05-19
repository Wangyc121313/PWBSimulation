"""Run one LNOI-PWB-SMF simulation and plot monitor fields."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LNOI_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PICTURES_DIR,
    PWBParameters,
    create_pwb_structure_in_fdtd,
    setup_fdtd_simulation,
    visualize_and_save_results,
)


SAVE_PATH = LNOI_DIR / "temp_single.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()
create_pwb_structure_in_fdtd(fdtd, params)
setup_fdtd_simulation(fdtd, params)
fdtd.save(str(SAVE_PATH))
fdtd.run()

PICTURES_DIR.mkdir(parents=True, exist_ok=True)
T_forward = visualize_and_save_results(fdtd, params, PICTURES_DIR / "single_run.jpg")
print(f"T_forward = {T_forward}")
fdtd.close()
