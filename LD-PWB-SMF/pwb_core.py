"""
Core helpers for LD-PWB-SMF simulations.

This module extracts the reusable notebook logic into script-friendly functions.
It does not import lumapi directly; callers create and pass an FDTD handle.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.colors import LogNorm

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import LD_BASE_FSP, LD_DIR, LD_RESULTS_DIR, LD_SOURCE_DATASET, MATERIAL_DB

PICTURES_DIR = LD_RESULTS_DIR / "Pictures"


class PWBParameters:
    def __init__(self):
        self.r = 100e-6
        self.l1 = 16e-6
        self.l2 = 100e-6
        self.curve_points = 200
        self.r_11 = 0.6e-6
        self.r_12 = 1.5e-6
        self.r_2 = 1.1e-6
        self.r_3 = 6.2e-6
        self.wavelength = 1.55e-6
        self.mesh_accuracy = 1
        self.simulation_time = 1500e-15
        self.total_length = 170e-6
        self.vertical_offset = 20e-6
        self.use_imported_source = True
        self.source_dataset = LD_SOURCE_DATASET


def generate_pwb_path(params):
    """Generate the LD-PWB-SMF centerline."""
    r = params.r
    l1 = params.l1
    l2 = params.l2
    l = (params.total_length - (l1 + l2)) / 2
    if l < 0:
        raise ValueError("params.l1 + params.l2 must not exceed params.total_length")
    if params.vertical_offset > r:
        raise ValueError("params.vertical_offset must not exceed params.r")

    theta = np.arcsin(params.vertical_offset / r)
    curve_points = params.curve_points

    x1 = np.linspace(0, l1, curve_points)
    z1 = np.zeros(curve_points)

    x2 = np.linspace(x1[-1], x1[-1] + l, curve_points)
    z2 = np.zeros(curve_points)

    t1 = np.linspace(0, theta, curve_points)
    center3_x = x2[-1]
    center3_z = r
    x3 = center3_x + r * np.sin(t1)
    z3 = center3_z - r * np.cos(t1)

    tangent_angle = theta
    center4_x = x3[-1] + r * np.sin(tangent_angle)
    center4_z = z3[-1] - r * np.cos(tangent_angle)
    t3 = np.linspace(np.pi / 2 - tangent_angle, np.pi / 2 + tangent_angle, curve_points)
    x4 = center4_x - r * np.cos(t3)
    z4 = center4_z + r * np.sin(t3)

    t4 = np.linspace(-tangent_angle, 0, curve_points)
    center5_x = x4[-1] + r * np.sin(abs(tangent_angle))
    center5_z = r
    x5 = center5_x + r * np.sin(t4)
    z5 = center5_z - r * np.cos(t4)

    x6 = np.linspace(x5[-1], x5[-1] + l, curve_points)
    z6 = np.zeros(curve_points)

    x7 = np.linspace(x6[-1], x6[-1] + l2, curve_points)
    z7 = np.zeros(curve_points)

    x = np.concatenate((x1, x2, x3, x4, x5, x6, x7))
    y = np.zeros_like(x)
    z = np.concatenate((z1, z2, z3, z4, z5, z6, z7))
    return np.column_stack((x, y, z))


def get_radius_at_position_1(t, params):
    """Radius profile for the first elliptical-axis radius."""
    point_idx = int(t * (7 * params.curve_points - 1))
    if point_idx < params.curve_points:
        local_t = point_idx / params.curve_points
        return params.r_11 - (params.r_11 - params.r_2) * local_t
    if point_idx < 6 * params.curve_points:
        return params.r_2
    local_t = (point_idx - 6 * params.curve_points) / params.curve_points
    return params.r_2 + (params.r_3 - params.r_2) * local_t


def get_radius_at_position_2(t, params):
    """Radius profile for the second elliptical-axis radius."""
    point_idx = int(t * (7 * params.curve_points - 1))
    if point_idx < params.curve_points:
        local_t = point_idx / params.curve_points
        return params.r_12 - (params.r_12 - params.r_2) * local_t
    if point_idx < 6 * params.curve_points:
        return params.r_2
    local_t = (point_idx - 6 * params.curve_points) / params.curve_points
    return params.r_2 + (params.r_3 - params.r_2) * local_t


def _set_segment_rotation(fdtd, direction):
    dz, dy, dx = direction[2], direction[1], direction[0]
    theta = np.arctan2(dx, dz) * 180 / np.pi
    phi = np.arcsin(np.clip(dy, -1.0, 1.0)) * 180 / np.pi
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", theta)
    fdtd.set("second axis", "x")
    fdtd.set("rotation 2", phi)


def _add_path_segments(
    fdtd,
    path,
    params,
    start_idx,
    end_idx,
    count,
    name_prefix,
    ellipsoid=False,
    radius_func_1=get_radius_at_position_1,
    radius_func_2=get_radius_at_position_2,
):
    path_length = len(path)
    span = end_idx - start_idx
    segment_length = max(1, span // count)

    for i in range(count):
        idx0 = start_idx + i * segment_length
        idx1 = min(start_idx + (i + 1) * segment_length, end_idx - 1)
        if idx0 >= idx1:
            continue

        start_pos = path[idx0]
        end_pos = path[idx1]
        center = (start_pos + end_pos) / 2
        direction = end_pos - start_pos
        segment_len = np.linalg.norm(direction)
        if segment_len <= 0:
            continue
        direction = direction / segment_len

        t0 = idx0 / (path_length - 1)
        t1 = idx1 / (path_length - 1)
        radius_1 = (radius_func_1(t0, params) + radius_func_1(t1, params)) / 2
        radius_2 = (radius_func_2(t0, params) + radius_func_2(t1, params)) / 2

        fdtd.addcircle()
        fdtd.set("name", f"{name_prefix}_{i}")
        fdtd.set("material", "Vancore B")
        fdtd.set("make ellipsoid", 1 if ellipsoid else 0)
        fdtd.set("x", center[0])
        fdtd.set("y", center[1])
        fdtd.set("z", center[2])
        fdtd.set("radius", radius_1)
        if ellipsoid:
            fdtd.set("radius 2", radius_2)
        fdtd.set("z span", segment_len)
        _set_segment_rotation(fdtd, direction)


def create_pwb_structure_in_fdtd(fdtd, params):
    """Create the LD-PWB-SMF PWB structure in FDTD and return its path."""
    path = generate_pwb_path(params)

    fdtd.deleteall()
    fdtd.load(str(LD_BASE_FSP))
    fdtd.importmaterialdb(str(MATERIAL_DB))

    fdtd.addcircle()
    fdtd.set("name", "Fiber Cladding")
    fdtd.set("material", "Fiber Cladding")
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.3)
    fdtd.set("x", 300e-6)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("make ellipsoid", 0)
    fdtd.set("radius", 50e-6)
    fdtd.set("z span", 100e-6)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", 90)

    fdtd.addcircle()
    fdtd.set("name", "Fiber Core")
    fdtd.set("material", "Fiber Core")
    fdtd.set("x", 300e-6)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("make ellipsoid", 0)
    fdtd.set("radius", 4.1e-6)
    fdtd.set("z span", 100e-6)
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", 90)

    fdtd.addrect()
    fdtd.set("name", "PWB Cladding")
    fdtd.set("material", "PWB Cladding")
    fdtd.set("override color opacity from material database", 1)
    fdtd.set("alpha", 0.3)
    fdtd.set("x", 125e-6)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("x span", 250e-6)
    fdtd.set("y span", 50e-6)
    fdtd.set("z span", 50e-6)

    cp = params.curve_points
    _add_path_segments(fdtd, path, params, 0, cp, 100, "PWB_taper_in", ellipsoid=True)
    _add_path_segments(fdtd, path, params, cp, 6 * cp, 100, "PWB_bend")
    _add_path_segments(fdtd, path, params, 6 * cp, 7 * cp, 100, "PWB_taper_out")
    return path


def setup_fdtd_simulation(fdtd, params, path):
    """Set up the LD-PWB-SMF FDTD region, source, monitors, and mode expansion."""
    y_min, y_max = np.min(path[:, 1]), np.max(path[:, 1])
    margin = 7e-6

    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("dimension", "3D")
    fdtd.set("x min", 140e-6)
    fdtd.set("x max", 270e-6)
    fdtd.set("y min", y_min - margin)
    fdtd.set("y max", y_max + margin)
    fdtd.set("z min", -margin)
    fdtd.set("z max", margin)
    fdtd.set("express mode", True)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    if params.use_imported_source and Path(params.source_dataset).exists():
        fdtd.addimportedsource()
        fdtd.set("name", "source")
        fdtd.importdataset(str(params.source_dataset))
        fdtd.set("x", 150e-6)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
    else:
        fdtd.addmode()
        fdtd.set("name", "source")
        fdtd.set("injection axis", "x-axis")
        fdtd.set("direction", "forward")
        fdtd.set("x", 150e-6)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
        fdtd.set("y span", 2 * margin)
        fdtd.set("z span", 2 * margin)
        fdtd.set("wavelength start", params.wavelength)
        fdtd.set("wavelength stop", params.wavelength)

    for name, x in (("input_monitor", 160e-6), ("output_monitor", 260e-6), ("output_monitor_2", 240e-6)):
        fdtd.addpower()
        fdtd.set("name", name)
        fdtd.set("monitor type", "2D X-normal")
        fdtd.set("x", x)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
        fdtd.set("y span", 2 * margin)
        fdtd.set("z span", 2 * margin)

    fdtd.addmodeexpansion()
    fdtd.set("name", "mode_expansion")
    fdtd.setexpansion("input", "output_monitor")
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("x", 260e-6)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)

    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x", 205e-6)
    fdtd.set("x span", 130e-6)
    fdtd.set("z", 0)
    fdtd.set("z span", 2 * margin)


def get_data(fdtd, params):
    """Read LD-PWB-SMF simulation results."""
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    input_E = fdtd.getresult("input_monitor", "E")
    output_E = fdtd.getresult("output_monitor", "E")
    mode_data = fdtd.getresult("mode_expansion", "expansion for input")
    T_forward = mode_data["T_forward"]
    return {
        "source_E": source_E,
        "transmission_E": transmission_E,
        "input_E": input_E,
        "output_E": output_E,
        "T_forward": T_forward,
    }


def _plot_monitor_x(ax, monitor_data, title):
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


def visualize_and_save_results(fdtd, params, output_path=None):
    """Visualize monitor fields and return T_forward."""
    results = get_data(fdtd, params)
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    _plot_monitor_x(ax1, results["source_E"], "Source")
    _plot_monitor_x(ax2, results["input_E"], "Input Monitor")
    _plot_monitor_x(ax3, results["output_E"], "Output Monitor")

    trans_e = results["transmission_E"]
    e_field = trans_e["E"]
    ex = e_field[:, 0, :, 0, 0]
    ey = e_field[:, 0, :, 0, 1]
    ez = e_field[:, 0, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    x = np.squeeze(trans_e["x"]) * 1e6
    z = np.squeeze(trans_e["z"]) * 1e6
    x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
    ax4.pcolormesh(x_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax4.set_xlabel("X (um)")
    ax4.set_ylabel("Z (um)")
    ax4.set_title("Transmission Monitor")

    plt.suptitle("LD-PWB-SMF Field")
    plt.tight_layout()
    if output_path is None:
        PICTURES_DIR.mkdir(parents=True, exist_ok=True)
        output_path = PICTURES_DIR / "single_run.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    return results["T_forward"]


def plot_transmission_monitor(fdtd):
    """Plot the Y-normal transmission monitor on a logarithmic color scale."""
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    fig, ax = plt.subplots(figsize=(15, 5))
    e_field = transmission_E["E"]
    ex = e_field[:, 0, :, 0, 0]
    ey = e_field[:, 0, :, 0, 1]
    ez = e_field[:, 0, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    x = np.squeeze(transmission_E["x"]) * 1e6
    z = np.squeeze(transmission_E["z"]) * 1e6
    x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
    norm = LogNorm(vmin=1e-5, vmax=np.nanmax(intensity))
    pcm = ax.pcolormesh(x_grid, z_grid, intensity, cmap="inferno", shading="auto", norm=norm)
    fig.colorbar(pcm, ax=ax, pad=0.02, format=ticker.LogFormatterSciNotation())
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Z (um)")
    ax.set_title("Transmission Monitor")
    plt.tight_layout()
    plt.show()
