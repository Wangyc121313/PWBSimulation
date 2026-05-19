"""
Build and save the complex planar PWB structure without running FDTD.

Use this script to inspect the generated geometry in Lumerical before running
full simulations or sweeps.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import PD_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_3,
    setup_fdtd_simulation_3,
)


SAVE_PATH = PD_DIR / "temp_complex.fsp"


params = PWBParameters()
fdtd = lumapi.FDTD()

generate_pwb_structure_3(fdtd, params)
setup_fdtd_simulation_3(fdtd, params)
fdtd.save(str(SAVE_PATH))
print("Complex planar PWB structure saved:", SAVE_PATH)
fdtd.close()
