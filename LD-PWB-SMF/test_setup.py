"""Build and save the LD-PWB-SMF structure without running FDTD."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LD_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import PWBParameters, create_pwb_structure_in_fdtd, setup_fdtd_simulation


SAVE_PATH = LD_DIR / "LD_with_PWB_temp.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()
path = create_pwb_structure_in_fdtd(fdtd, params)
setup_fdtd_simulation(fdtd, params, path)
fdtd.save(str(SAVE_PATH))
print("LD-PWB-SMF structure saved:", SAVE_PATH)
fdtd.close()
