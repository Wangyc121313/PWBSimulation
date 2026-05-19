"""
Core helpers for LNOI-PWB-SMF simulations.

This module extracts reusable logic from simulation-LNOI-PWB-SMF.ipynb.
Callers create and pass the lumapi FDTD handle.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LNOI_BASE_FSP, LNOI_RESULTS_DIR, MATERIAL_DB

PICTURES_DIR = LNOI_RESULTS_DIR / "Pictures"


class PWBParameters:
    def __init__(self):
        self.h1 = 1.3e-6
        self.h2 = 1.5e-6
        self.w1 = 3.5e-6
        self.w2 = 3.0e-6
        self.l = 80e-6
        self.straight_length = 30e-6
        self.straight_offset = 15e-6
        self.wavelength = 1.55e-6
        self.mesh_accuracy = 1
        self.simulation_time = 2000e-15


def create_pwb_structure_in_fdtd(fdtd, params):
    """Create the LNOI-PWB-SMF taper structure in FDTD."""
    h1, h2, w1, w2, l = params.h1, params.h2, params.w1, params.w2, params.l

    fdtd.deleteall()
    fdtd.load(str(LNOI_BASE_FSP))
    fdtd.importmaterialdb(str(MATERIAL_DB))

    fdtd.addpyramid()
    fdtd.set("name", "PWB")
    fdtd.set("material", "Vancore B")
    fdtd.set("override mesh order from material database", 1)
    fdtd.set("mesh order", 3)
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.5)
    fdtd.set("x", l / 2)
    fdtd.set("x span bottom", 2 * h1)
    fdtd.set("x span top", 2 * h2)
    fdtd.set("y", 0)
    fdtd.set("y span bottom", w1)
    fdtd.set("y span top", w2)
    fdtd.set("z", 0)
    fdtd.set("z span", l)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", 90)

    fdtd.addrect()
    fdtd.set("name", "PWB-Straight")
    fdtd.set("material", "Vancore B")
    fdtd.set("override mesh order from material database", 1)
    fdtd.set("mesh order", 3)
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.5)
    fdtd.set("x", params.straight_offset + l)
    fdtd.set("x span", params.straight_length)
    fdtd.set("y", 0)
    fdtd.set("y span", w2)
    fdtd.set("z", h2 / 2)
    fdtd.set("z span", h2)


def setup_fdtd_simulation(fdtd, params):
    """Set up FDTD region, source, monitors, and mode expansion."""
    h2, w2 = params.h2, params.w2
    margin = 2e-6

    fdtd.addfdtd()
    fdtd.set("dimension", "3D")
    fdtd.set("x min", -10e-6)
    fdtd.set("x max", 95e-6)
    fdtd.set("z min", -2e-6)
    fdtd.set("z max", h2 + margin)
    fdtd.set("y min", -w2 / 2 - margin)
    fdtd.set("y max", w2 / 2 + margin)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("x", -5e-6)
    fdtd.set("y min", -w2 / 2 - margin)
    fdtd.set("y max", w2 / 2 + margin)
    fdtd.set("z min", -2e-6)
    fdtd.set("z max", h2 + margin)
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    for name, x in (
        ("monitor_0", 0),
        ("monitor_60", 60e-6),
        ("monitor_80", 80e-6),
        ("monitor_90", 90e-6),
    ):
        fdtd.addpower()
        fdtd.set("name", name)
        fdtd.set("monitor type", "2D X-normal")
        fdtd.set("x", x)
        fdtd.set("y min", -w2 / 2 - margin)
        fdtd.set("y max", w2 / 2 + margin)
        fdtd.set("z min", -2e-6)
        fdtd.set("z max", h2 + margin)

    fdtd.addmodeexpansion()
    fdtd.set("name", "mode_expansion")
    fdtd.setexpansion("output", "monitor_90")
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("x", 90e-6)
    fdtd.set("y min", -w2 / 2 - margin)
    fdtd.set("y max", w2 / 2 + margin)
    fdtd.set("z min", -2e-6)
    fdtd.set("z max", h2 + margin)

    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x min", -10e-6)
    fdtd.set("x max", 95e-6)
    fdtd.set("z min", -2e-6)
    fdtd.set("z max", h2 + margin)


def get_data(fdtd, params):
    """Read LNOI-PWB-SMF simulation results."""
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    monitor_0_E = fdtd.getresult("monitor_0", "E")
    monitor_60_E = fdtd.getresult("monitor_60", "E")
    monitor_80_E = fdtd.getresult("monitor_80", "E")
    monitor_90_E = fdtd.getresult("monitor_90", "E")
    mode_data = fdtd.getresult("mode_expansion", "expansion for output")
    T_forward = mode_data["T_forward"][0][0]
    return {
        "source_E": source_E,
        "transmission_E": transmission_E,
        "monitor_0_E": monitor_0_E,
        "monitor_60_E": monitor_60_E,
        "monitor_80_E": monitor_80_E,
        "monitor_90_E": monitor_90_E,
        "T_forward": T_forward,
    }


def _plot_x_monitor(ax, monitor_data, title):
    e_field = monitor_data["E"]
    ex = e_field[0, :, :, 0, 0]
    ey = e_field[0, :, :, 0, 1]
    ez = e_field[0, :, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    y = np.squeeze(monitor_data["y"]) * 1e6
    z = np.squeeze(monitor_data["z"]) * 1e6
    y_grid, z_grid = np.meshgrid(y, z, indexing="ij")
    ax.pcolormesh(y_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax.set_xlabel("Y (um)")
    ax.set_ylabel("Z (um)")
    ax.set_title(title)


def visualize_and_save_results(fdtd, params, output_path=None, title="LNOI-PWB-SMF Field"):
    """Visualize monitor fields and return T_forward."""
    results = get_data(fdtd, params)
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    ax1, ax2, ax3 = axes[0]
    ax4, ax5, ax6 = axes[1]
    _plot_x_monitor(ax1, results["source_E"], "Source")
    _plot_x_monitor(ax2, results["monitor_0_E"], "Monitor_0")
    _plot_x_monitor(ax3, results["monitor_60_E"], "Monitor_60")
    _plot_x_monitor(ax4, results["monitor_80_E"], "Monitor_80")
    _plot_x_monitor(ax5, results["monitor_90_E"], "Monitor_90")

    trans_e = results["transmission_E"]
    e_field = trans_e["E"]
    ex = e_field[:, 0, :, 0, 0]
    ey = e_field[:, 0, :, 0, 1]
    ez = e_field[:, 0, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    x = np.squeeze(trans_e["x"]) * 1e6
    z = np.squeeze(trans_e["z"]) * 1e6
    x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
    ax6.pcolormesh(x_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax6.set_xlabel("X (um)")
    ax6.set_ylabel("Z (um)")
    ax6.set_title("Transmission Monitor")

    plt.suptitle(title)
    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    return results["T_forward"]
