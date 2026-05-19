"""Run the LNOI-PWB-SMF h1/h2/w1/w2 parameter sweep."""

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LNOI_DIR, LNOI_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PWBParameters,
    create_pwb_structure_in_fdtd,
    get_data,
    setup_fdtd_simulation,
    visualize_and_save_results,
)


RESULTS_DIR = LNOI_RESULTS_DIR / "h1_h2_w1_w2_scan"
RESULTS_FILE = RESULTS_DIR / "T_forward_results.csv"
SAVE_PATH = LNOI_DIR / "temp_sweep.fsp"


RESULTS_DIR.mkdir(parents=True, exist_ok=True)
file_exists = RESULTS_FILE.exists()
with open(RESULTS_FILE, "a", newline="") as csvfile:
    fieldnames = ["h1", "h2", "w1", "w2", "T_forward"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    if not file_exists:
        writer.writeheader()


h1_values = np.arange(0.7, 1.7, 0.2) * 1e-6
h2_values = np.arange(0.7, 1.7, 0.2) * 1e-6
w1_values = np.arange(1.5, 4.0, 0.5) * 1e-6
w2_values = np.arange(1.5, 4.0, 0.5) * 1e-6


for h1 in h1_values:
    for h2 in h2_values:
        for w1 in w1_values:
            for w2 in w2_values:
                print(
                    f"Processing h1={h1*1e6:.1f} um, h2={h2*1e6:.1f} um, "
                    f"w1={w1*1e6:.1f} um, w2={w2*1e6:.1f} um"
                )

                params = PWBParameters()
                params.h1 = h1
                params.h2 = h2
                params.w1 = w1
                params.w2 = w2

                fdtd = lumapi.FDTD()
                create_pwb_structure_in_fdtd(fdtd, params)
                setup_fdtd_simulation(fdtd, params)
                fdtd.save(str(SAVE_PATH))
                fdtd.run()

                results = get_data(fdtd, params)
                T_forward = results["T_forward"]

                with open(RESULTS_FILE, "a", newline="") as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=["h1", "h2", "w1", "w2", "T_forward"])
                    writer.writerow({
                        "h1": h1 * 1e6,
                        "h2": h2 * 1e6,
                        "w1": w1 * 1e6,
                        "w2": w2 * 1e6,
                        "T_forward": T_forward,
                    })

                fig_path = RESULTS_DIR / (
                    f"h1_{h1*1e6:.1f}_h2_{h2*1e6:.1f}_"
                    f"w1_{w1*1e6:.1f}_w2_{w2*1e6:.1f}.jpg"
                )
                try:
                    visualize_and_save_results(
                        fdtd,
                        params,
                        fig_path,
                        title=(
                            f"h1={h1*1e6:.1f} um, h2={h2*1e6:.1f} um, "
                            f"w1={w1*1e6:.1f} um, w2={w2*1e6:.1f} um"
                        ),
                    )
                except Exception as exc:
                    print(f"Error while saving plot: {exc}")

                print("T_forward:", T_forward)
                fdtd.close()


print("Parameter scan completed. Results saved:", RESULTS_FILE)
