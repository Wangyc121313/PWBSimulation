"""Build and save the LNOI-PWB-SMF structure without running FDTD."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LNOI_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import PWBParameters, create_pwb_structure_in_fdtd, setup_fdtd_simulation


SAVE_PATH = LNOI_DIR / "test_setup.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()
create_pwb_structure_in_fdtd(fdtd, params)
setup_fdtd_simulation(fdtd, params)
fdtd.save(str(SAVE_PATH))
print("LNOI-PWB-SMF structure saved:", SAVE_PATH)
fdtd.close()
