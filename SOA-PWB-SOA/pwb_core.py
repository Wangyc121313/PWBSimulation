"""
Core helpers for SOA-PWB-SOA simulations.

Models the on-chip SOA output → PWB → external SOA input coupling path.
Both SOA waveguides are loaded from a base .fsp file (SOA_base.fsp).
The PWB is built as circular/elliptical segments using addcircle().

Structure (straight line along +x, total length = 250 um):
  [SOA_out, from .fsp] → Taper-1 (radius expansion) → PWB straight
  → Taper-2 (radius compression) → [SOA_in, from .fsp]

This module does not import lumapi directly; callers create and pass an FDTD handle.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import MATERIAL_DB, SOA_BASE_FSP, SOA_RESULTS_DIR

PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"


class SOAPWBParams:
    """Parameters for SOA-PWB-SOA simulation.

    All geometric dimensions in SI units (meters). Use X * 1e-6 for microns.

    PWB cross-section is circular/elliptical, defined by radius (first axis)
    and optionally radius_2 (second axis for ellipsoid).  Tapers transition
    between the SOA facet mode size and the PWB cross-section.
    """

    def __init__(self):
        # ---- Wavelength & simulation ----
        self.wavelength = 1.55e-6          # center wavelength [m]
        self.mesh_accuracy =1
        self.simulation_time = 2000e-15    # [s]

        # ---- Total PWB length constraint ----
        self.total_length = 250e-6         # total PWB length (tapers + straight) [m]

        # ---- Taper-1: SOA output → PWB (radius expansion) ----
        self.taper1_length = 30e-6         # [m]

        # ---- Taper-2: PWB → external SOA input (radius compression) ----
        self.taper2_length = 30e-6         # [m]

        # PWB straight section length is derived:
        #   pwb_length = total_length - taper1_length - taper2_length

        # ---- PWB radius profile ----
        # First axis radius along the path:
        #   taper1_start (r_in) → taper1_end (r_pwb) → taper2_start → taper2_end (r_out)
        self.r_in = 1.25e-6                 # start radius at SOA output facet [m]
        self.r_pwb = 1.1e-6                # middle PWB section radius [m]
        self.r_out = 1.0e-6                # end radius at external SOA input facet [m]

        # ---- Elliptical cross-section support ----
        # If the PWB needs an elliptical cross-section, set use_ellipsoid = True
        # and provide second-axis radius values.  Otherwise the cross-section
        # is circular (radius_2 defaults to the first-axis value).
        self.use_ellipsoid = True
        self.r_in_2 =3.0e-6               # second-axis radius at taper1 start [m]
        self.r_pwb_2 = 1.1e-6              # second-axis radius in PWB middle [m]
        self.r_out_2 = 2.0e-6              # second-axis radius at taper2 end [m]

        # ---- Path discretisation ----
        self.curve_points = 200            # number of centreline sample points

        # ---- Materials ----
        self.material_pwb = "Vancore B"    # PWB polymer material


# ---------------------------------------------------------------------------
# Centreline generation (straight path, no bend)
# ---------------------------------------------------------------------------

def generate_pwb_path(params):
    """Generate a straight centreline along +x for SOA-PWB-SOA.

    Returns (N, 3) array of [x, y, z] coordinates.
    """
    n = params.curve_points
    x = np.linspace(0, params.total_length, n)
    y = np.zeros(n)
    z = np.zeros(n)
    return np.column_stack((x, y, z))


# ---------------------------------------------------------------------------
# Radius profiles along the path
# ---------------------------------------------------------------------------

def get_radius_at_position_1(t, params):
    """First-axis radius at normalised path position t ∈ [0, 1].

    Profile: r_in → r_pwb (taper-1) | r_pwb constant (PWB) | r_pwb → r_out (taper-2)
    """
    frac1 = params.taper1_length / params.total_length
    frac2 = 1.0 - params.taper2_length / params.total_length

    if t < frac1:
        local_t = t / frac1
        return params.r_in + (params.r_pwb - params.r_in) * local_t
    elif t < frac2:
        return params.r_pwb
    else:
        local_t = (t - frac2) / (1.0 - frac2)
        return params.r_pwb + (params.r_out - params.r_pwb) * local_t


def get_radius_at_position_2(t, params):
    """Second-axis radius at normalised path position t ∈ [0, 1].

    Profile: r_in_2 → r_pwb_2 (taper-1) | r_pwb_2 constant (PWB) | r_pwb_2 → r_out_2 (taper-2)
    """
    if not params.use_ellipsoid:
        return get_radius_at_position_1(t, params)

    frac1 = params.taper1_length / params.total_length
    frac2 = 1.0 - params.taper2_length / params.total_length

    if t < frac1:
        local_t = t / frac1
        return params.r_in_2 + (params.r_pwb_2 - params.r_in_2) * local_t
    elif t < frac2:
        return params.r_pwb_2
    else:
        local_t = (t - frac2) / (1.0 - frac2)
        return params.r_pwb_2 + (params.r_out_2 - params.r_pwb_2) * local_t


# ---------------------------------------------------------------------------
# Segment rotation helper (same as LD-PWB-SMF)
# ---------------------------------------------------------------------------

def _set_segment_rotation(fdtd, direction):
    """Orient a circle primitive so its normal aligns with `direction`."""
    dz, dy, dx = direction[2], direction[1], direction[0]
    theta = np.arctan2(dx, dz) * 180 / np.pi
    phi = np.arcsin(np.clip(dy, -1.0, 1.0)) * 180 / np.pi
    fdtd.set("first axis", "y")
    fdtd.set("rotation 1", theta)
    fdtd.set("second axis", "x")
    fdtd.set("rotation 2", phi)


# ---------------------------------------------------------------------------
# Path-segment construction (adapted from LD-PWB-SMF)
# ---------------------------------------------------------------------------

def _add_path_segments(
    fdtd,
    path,
    params,
    start_idx,
    end_idx,
    count,
    name_prefix,
):
    """Place `count` circle/ellipsoid primitives along a subsection of `path`.

    Each segment is an addcircle() primitive.  When `params.use_ellipsoid` is
    True the primitive has two independent radii (elliptical cross-section).
    """
    path_length = len(path)
    span = end_idx - start_idx
    if span <= 0:
        return
    segment_length = max(1, span // count)

    for i in range(count):
        idx0 = start_idx + i * segment_length
        idx1 = min(start_idx + (i + 1) * segment_length, end_idx)
        if idx0 >= idx1 or idx1 >= path_length:
            continue

        start_pos = path[idx0]
        end_pos = path[idx1]
        center = (start_pos + end_pos) / 2
        direction = end_pos - start_pos
        seg_len = np.linalg.norm(direction)
        if seg_len <= 0:
            continue
        direction = direction / seg_len

        t0 = idx0 / (path_length - 1)
        t1 = idx1 / (path_length - 1)
        r1 = (get_radius_at_position_1(t0, params)
              + get_radius_at_position_1(t1, params)) / 2
        r2 = (get_radius_at_position_2(t0, params)
              + get_radius_at_position_2(t1, params)) / 2

        fdtd.addcircle()
        fdtd.set("name", f"{name_prefix}_{i}")
        fdtd.set("material", params.material_pwb)
        fdtd.set("make ellipsoid", 1 if params.use_ellipsoid else 0)
        fdtd.set("x", center[0])
        fdtd.set("y", center[1])
        fdtd.set("z", center[2])
        fdtd.set("radius", r1)
        if params.use_ellipsoid:
            fdtd.set("radius 2", r2)
        fdtd.set("z span", seg_len)
        _set_segment_rotation(fdtd, direction)


# ---------------------------------------------------------------------------
# Main structure creation
# ---------------------------------------------------------------------------

def create_pwb_structure_in_fdtd(fdtd, params):
    """Build the PWB structure in FDTD and return its centreline path.

    Loads SOA_base.fsp (which contains both on-chip and external SOA
    waveguides at their correct positions), then builds the PWB with
    tapers as circular/elliptical segments along the +x axis.

    Args:
        fdtd: lumapi.FDTD handle (already opened by caller).
        params: SOAPWBParams instance.

    Returns:
        np.ndarray: (N, 3) centreline path.
    """
    path = generate_pwb_path(params)

    fdtd.deleteall()
    fdtd.load(str(SOA_BASE_FSP))
    fdtd.importmaterialdb(str(MATERIAL_DB))

    cp = params.curve_points

    # Segment boundaries derived from actual taper fractions
    frac1 = params.taper1_length / params.total_length
    frac2 = 1.0 - params.taper2_length / params.total_length
    idx_taper1_end = max(1, int(frac1 * cp))
    idx_taper2_start = min(cp - 1, int(frac2 * cp))

    # Each section uses count = span so every path step gets a segment.
    # end_idx for the last section is cp-1 (last valid path index);
    # for other sections it is the neighbour's start index (inclusive boundary).

    # Taper-1: radius expansion  (r_in → r_pwb)
    _add_path_segments(fdtd, path, params,
                       0, idx_taper1_end, idx_taper1_end, "PWB_taper1")

    # PWB straight: constant radius (r_pwb)
    _add_path_segments(fdtd, path, params,
                       idx_taper1_end, idx_taper2_start,
                       idx_taper2_start - idx_taper1_end, "PWB_straight")

    # Taper-2: radius compression (r_pwb → r_out)
    _add_path_segments(fdtd, path, params,
                       idx_taper2_start, cp - 1,
                       cp - 1 - idx_taper2_start, "PWB_taper2")

    return path


# ---------------------------------------------------------------------------
# FDTD simulation setup
# ---------------------------------------------------------------------------

def setup_fdtd_simulation(fdtd, params, path):
    """Configure the FDTD simulation region, source, monitors, and mode expansion.

    Args:
        fdtd: lumapi.FDTD handle.
        params: SOAPWBParams instance.
        path: (N, 3) centreline from generate_pwb_path().
    """
    # ---- Lateral margin ----
    max_radius = max(params.r_pwb, params.r_in, params.r_out,
                     params.r_in_2, params.r_pwb_2, params.r_out_2)
    margin = max_radius  + 2e-6  # extra margin for PML and field decay

    # FDTD region: include a short SOA waveguide lead-in (from .fsp) before
    # the PWB taper starts at x=0, so the mode source can be placed inside
    # the SOA output waveguide.
    soa_lead_in = 10e-6               # length of SOA output waveguide in simulation
    x_min = -soa_lead_in
    x_max = params.total_length + params.taper2_length * 0.3

    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("express mode", True)
    fdtd.set("dimension", "3D")
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("y min", -margin)
    fdtd.set("y max", margin)
    fdtd.set("z min", -margin)
    fdtd.set("z max", margin)
    fdtd.set("x min bc", "PML")
    fdtd.set("x max bc", "PML")
    fdtd.set("y min bc", "PML")
    fdtd.set("y max bc", "PML")
    fdtd.set("z min bc", "PML")
    fdtd.set("z max bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    # ---- Mode source: placed inside the SOA output waveguide (from .fsp),
    #      before its output facet at x=0, to correctly excite the SOA
    #      waveguide's fundamental mode.
    source_x = -soa_lead_in * 0.5      # centre of SOA lead-in segment
    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("x", source_x)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    # ---- X-normal power monitors at key positions ----
    # Positions (in x, relative to PWB path: x=0 is SOA output facet):
    #   input  — at SOA output facet / PWB taper-1 start
    #   pwb_in — at end of taper-1 (start of PWB straight section)
    #   pwb_out — at end of PWB straight section (start of taper-2)
    #   output — at end of taper-2 (external SOA input facet)
    mon_input = 5e-6  # just inside SOA output facet
    mon_pwb_in = params.taper1_length
    mon_pwb_out = params.total_length - params.taper2_length
    mon_output = params.total_length + 5e-6  # just inside external SOA input facet

    monitor_positions = [
        ("input_monitor", mon_input),
        ("pwb_in_monitor", mon_pwb_in),
        ("pwb_out_monitor", mon_pwb_out),
        ("output_monitor", mon_output),
    ]
    for name, mx in monitor_positions:
        fdtd.addpower()
        fdtd.set("name", name)
        fdtd.set("monitor type", "2D X-normal")
        fdtd.set("x", mx)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
        fdtd.set("y span", 2 * margin)
        fdtd.set("z span", 2 * margin)

    # ---- Mode expansion at several positions to track fundamental-mode
    #      transmission through each section ----
    expansion_monitors = [
        ("mode_exp_pwb_in", mon_input),    # after taper-1
        ("mode_exp_output", mon_output),    # after taper-2
    ]
    for name, mx in expansion_monitors:
        fdtd.addmodeexpansion()
        fdtd.set("name", name)
        fdtd.setexpansion("input", name.replace("mode_exp_", "") + "_monitor")
        fdtd.set("mode selection", "fundamental TE mode")
        fdtd.set("x", mx)
        fdtd.set("y", 0)
        fdtd.set("z", 0)
        fdtd.set("y span", 2 * margin)
        fdtd.set("z span", 2 * margin)

    # ---- Y-normal transmission monitor (side view for field propagation) ----
    fdtd.addpower()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("z", 0)
    fdtd.set("z span", 2 * margin)


# ---------------------------------------------------------------------------
# Result extraction
# ---------------------------------------------------------------------------

def get_data(fdtd, params):
    """Read simulation results from all monitors.

    Returns:
        dict with keys: source_E, transmission_E, monitors dict, T_forward.
    """
    source_E = fdtd.getresult("source", "mode profile")
    transmission_E = fdtd.getresult("transmission_monitor", "E")
    input_E = fdtd.getresult("input_monitor", "E")
    pwb_in_E = fdtd.getresult("pwb_in_monitor", "E")
    pwb_out_E = fdtd.getresult("pwb_out_monitor", "E")
    output_E = fdtd.getresult("output_monitor", "E")

    # T_forward after taper-2 (overall transmission into fundamental mode)
    try:
        mode_data_out = fdtd.getresult("mode_exp_output", "expansion for input")
        T_forward = mode_data_out.get("T_forward", None)
    except Exception:
        T_forward = None

    # T_forward after taper-1 (transmission into PWB fundamental mode)
    try:
        mode_data_pwb = fdtd.getresult("mode_exp_pwb_in", "expansion for input")
        T_pwb_in = mode_data_pwb.get("T_forward", None)
    except Exception:
        T_pwb_in = None

    return {
        "source_E": source_E,
        "transmission_E": transmission_E,
        "input_E": input_E,
        "pwb_in_E": pwb_in_E,
        "pwb_out_E": pwb_out_E,
        "output_E": output_E,
        "T_forward": T_forward,
        "T_pwb_in": T_pwb_in,
    }


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def _plot_x_normal_monitor(ax, monitor_data, title):
    """Plot |E|^2 intensity on an X-normal (Y-Z) monitor plane."""
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
    """Generate field visualisation figure and return T_forward.

    Creates a 2×3 subplot: source mode + 4 cross-section monitors
    + transmission side view.
    """
    results = get_data(fdtd, params)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    (ax_src, ax_in, ax_pwb_in), (ax_pwb_out, ax_out, ax_trans) = axes

    _plot_x_normal_monitor(ax_src, results["source_E"], "Source Mode")
    _plot_x_normal_monitor(ax_in, results["input_E"], "Input (SOA facet / Taper-1 start)")
    _plot_x_normal_monitor(ax_pwb_in, results["pwb_in_E"], "PWB In (Taper-1 end)")
    _plot_x_normal_monitor(ax_pwb_out, results["pwb_out_E"], "PWB Out (Taper-2 start)")
    _plot_x_normal_monitor(ax_out, results["output_E"], "Output (Taper-2 end)")

    # Transmission monitor (side view: X-Z plane, Y=0)
    trans_e = results["transmission_E"]
    e_field = trans_e["E"]
    ex = e_field[:, 0, :, 0, 0]
    ey = e_field[:, 0, :, 0, 1]
    ez = e_field[:, 0, :, 0, 2]
    intensity = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2

    x = np.squeeze(trans_e["x"]) * 1e6
    z = np.squeeze(trans_e["z"]) * 1e6
    x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
    ax_trans.pcolormesh(x_grid, z_grid, intensity, cmap="jet", shading="auto")
    ax_trans.set_xlabel("X (um)")
    ax_trans.set_ylabel("Z (um)")
    ax_trans.set_title("Field Propagation (X-Z)")

    plt.suptitle("SOA-PWB-SOA Field Distribution", fontsize=14)
    plt.tight_layout()

    if output_path is None:
        PICTURES_DIR.mkdir(parents=True, exist_ok=True)
        output_path = PICTURES_DIR / "single_run.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

    return results["T_forward"]


# ============================================================
# FDE Mode Analysis (requires Lumerical MODE Solutions)
# ============================================================
# NOTE: FDE is secondary for this project.  The primary optimisation
# workflow uses FDTD mode-expansion monitors placed just after each
# taper to directly read fundamental-mode transmission efficiency.
# FDE is kept for sanity checks (e.g. verifying that chosen PWB
# cross-sections support a guided mode).


def run_fde_mode_analysis(params, mode_obj=None):
    """Compute fundamental mode field distributions at key cross-sections.

    Uses Lumerical MODE Solutions FDE (Finite-Difference Eigenmode) solver.

    Sections analyzed:
      1. SOA output waveguide cross-section
      2. PWB polymer waveguide cross-section
      3. External SOA input waveguide cross-section

    Args:
        params: SOAPWBParams instance.
        mode_obj: lumapi.MODE handle. If None, creates a new one.

    Returns:
        dict with keys: soa_mode, pwb_mode, ext_soa_mode.
        Each contains: E (3D array), x, y, neff.
    """
    close_on_exit = mode_obj is None
    if mode_obj is None:
        try:
            import lumapi
            mode_obj = lumapi.MODE()
        except AttributeError:
            print("WARNING: MODE Solutions not available. FDE analysis skipped.")
            return None

    try:
        results = {}

        # --- 1. SOA output waveguide mode ---
        # Dimensions from the SOA_base.fsp — these are the ridge waveguide
        # cross-section at the output facet (5 um × 2 um, with 200 nm active).
        results["soa_mode"] = _compute_single_mode(
            mode_obj, params,
            width=5.0e-6,
            height=2.0e-6,
            material="InP",
            label="SOA_output",
        )

        # --- 2. PWB waveguide mode ---
        # Effective circular/elliptical cross-section at the PWB middle.
        # Modeled as a rectangular approximation for FDE purposes.
        results["pwb_mode"] = _compute_single_mode(
            mode_obj, params,
            width=params.r_pwb * 2,
            height=params.r_pwb_2 * 2 if params.use_ellipsoid else params.r_pwb * 2,
            material=params.material_pwb,
            label="PWB",
        )

        # --- 3. External SOA input waveguide mode ---
        results["ext_soa_mode"] = _compute_single_mode(
            mode_obj, params,
            width=5.0e-6,
            height=2.0e-6,
            material="InP",
            label="Ext_SOA_input",
        )

        return results

    finally:
        if close_on_exit and mode_obj is not None:
            mode_obj.close()


def _compute_single_mode(mode, params, width, height, material, label):
    """Compute fundamental TE mode for a single rectangular waveguide cross-section.

    Args:
        mode: lumapi.MODE handle.
        params: SOAPWBParams instance.
        width, height: waveguide cross-section dimensions [m].
        material: material name string.
        label: name for the solver region.

    Returns:
        dict with E, x, y, neff, or None on failure.
    """
    margin = max(width, height) * 2

    mode.addfde()
    mode.set("name", label)
    mode.set("solver type", "2D X normal")
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z", 0)
    mode.set("y span", width + margin)
    mode.set("z span", height + margin)
    mode.set("wavelength", params.wavelength)
    mode.set("number of trial modes", 5)
    mode.set("mesh refinement", "conformal variant 1")

    # Add waveguide structure in MODE
    mode.addrect()
    mode.set("name", f"{label}_core")
    mode.set("material", material)
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z", 0)
    mode.set("x span", 1e-9)  # thin in propagation direction
    mode.set("y span", width)
    mode.set("z span", height)

    # Add substrate
    mode.addrect()
    mode.set("name", f"{label}_substrate")
    mode.set("material", "SiO2 (Glass) - Palik")
    mode.set("x", 0)
    mode.set("y", 0)
    mode.set("z min", -margin * 2)
    mode.set("z max", 0)
    mode.set("x span", 1e-9)
    mode.set("y span", width + margin * 2)

    mode.run()
    mode.findmodes()

    try:
        neff = mode.getdata("FDE::data::mode1", "neff")
        E = mode.getresult("FDE::data::mode1", "mode profile")
        x = mode.getdata("FDE::data::mode1", "x")
        y = mode.getdata("FDE::data::mode1", "y")
        mode.deleteall()
        return {"E": E, "x": x, "y": y, "neff": neff}
    except Exception:
        mode.deleteall()
        return None


def compute_mode_overlap(mode1, mode2):
    """Compute the mode overlap factor eta between two modes.

    eta = |integral E1(x,y) . E2*(x,y) dxdy|^2 /
          [integral|E1|^2 dxdy * integral|E2|^2 dxdy]

    Args:
        mode1, mode2: dicts from _compute_single_mode() with keys E, x, y.

    Returns:
        float: overlap factor eta (0 to 1), or None on failure.
    """
    if mode1 is None or mode2 is None:
        return None

    E1 = mode1["E"]  # shape: (3, ny, nx, ...)
    E2 = mode2["E"]

    # Extract transverse components: Ex, Ey, Ez at the central slice
    e1_x = E1[0, :, :, 0, 0].flatten()
    e1_y = E1[1, :, :, 0, 0].flatten()
    e1_z = E1[2, :, :, 0, 0].flatten()

    e2_x = E2[0, :, :, 0, 0].flatten()
    e2_y = E2[1, :, :, 0, 0].flatten()
    e2_z = E2[2, :, :, 0, 0].flatten()

    # Full vector overlap
    numerator = np.abs(
        np.dot(e1_x, np.conj(e2_x))
        + np.dot(e1_y, np.conj(e2_y))
        + np.dot(e1_z, np.conj(e2_z))
    ) ** 2

    denom1 = (
        np.dot(e1_x, np.conj(e1_x))
        + np.dot(e1_y, np.conj(e1_y))
        + np.dot(e1_z, np.conj(e1_z))
    )
    denom2 = (
        np.dot(e2_x, np.conj(e2_x))
        + np.dot(e2_y, np.conj(e2_y))
        + np.dot(e2_z, np.conj(e2_z))
    )

    if denom1 == 0 or denom2 == 0:
        return 0.0

    return float(numerator / (denom1 * denom2))


def plot_mode_profiles(fde_results, output_path=None):
    """Plot normalized mode field distributions from FDE analysis.

    Args:
        fde_results: dict from run_fde_mode_analysis().
        output_path: optional path for saving the figure.
    """
    if fde_results is None:
        print("No FDE results to plot.")
        return

    plt.rcParams["font.family"] = "Times New Roman"
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    labels = ["SOA Output Mode", "PWB Mode", "Ext SOA Input Mode"]
    keys = ["soa_mode", "pwb_mode", "ext_soa_mode"]

    for ax, label, key in zip(axes, labels, keys):
        mode = fde_results.get(key)
        if mode is None:
            ax.set_title(f"{label}\n(not available)")
            continue
        E = mode["E"]
        intensity = (
            np.abs(E[0, :, :, 0, 0]) ** 2
            + np.abs(E[1, :, :, 0, 0]) ** 2
            + np.abs(E[2, :, :, 0, 0]) ** 2
        )
        y = np.squeeze(mode["y"]) * 1e6
        z = np.squeeze(mode["x"]) * 1e6  # MODE may swap x/y axes
        y_grid, z_grid = np.meshgrid(y, z, indexing="ij")
        ax.pcolormesh(y_grid, z_grid, intensity / np.max(intensity),
                       cmap="hot", shading="auto")
        ax.set_xlabel("Y (um)")
        ax.set_ylabel("Z (um)")
        neff = mode.get("neff", "?")
        ax.set_title(f"{label}\nn_eff = {neff}")

    plt.suptitle("Fundamental Mode Profiles (|E|^2 normalized)", fontsize=14)
    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
