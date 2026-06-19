"""
CFD Manifold Post-Processing
=============================
Reads Ansys Fluent surface-report exports and produces:
  - Bar chart of outlet mass flow rates (coarse mesh, Geometry 1)
  - Bar chart of outlet mass flow rates (refined mesh, Geometry 1)
  - Side-by-side coarse vs. refined mesh comparison (Geometry 1)
  - Bar chart for Geometry 2 refined mesh  ← only when data is present
  - Side-by-side Geometry 1 vs. Geometry 2 comparison  ← only when Geom2 data is present

Expected folder layout (relative to repository root):

    Geometry_1/
        coarse/   flow_rates  outlets_y  outlets_z
        refined/  mass_flow_rates  y_coord  z_coord
                  avg_tot_pressure_out  pressure_in
    Geometry_2/
        refined/  mass_flow_rates  y_coord  z_coord
                  avg_tot_pressure_out  pressure_in
                  (folder is empty until Geometry 2 simulation is complete)

Plots are saved next to the data files; the geometry comparison plot goes
to the repository root.

Run from any working directory – all paths are resolved from this file's location.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from utils import THRESHOLD, parse_fluent_report, sort_outlets_by_angle

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Repository root is one level above the scripts/ folder
ROOT = Path(__file__).parent.parent

GEOM1_DIR    = ROOT / "Geometry_1"
GEOM1_COARSE = GEOM1_DIR / "coarse"
GEOM1_FINE   = GEOM1_DIR / "refined"

GEOM2_DIR  = ROOT / "Geometry_2"
GEOM2_FINE = GEOM2_DIR / "refined"


# ---------------------------------------------------------------------------
# Single-variant analysis
# ---------------------------------------------------------------------------

def analyse_variant(data_dir, label, mfr_file, y_file, z_file,
                    pin_file=None, out_file=None):
    """
    Load Fluent exports from *data_dir*, compute distribution statistics,
    print a summary table to stdout, and save a bar-chart PNG.

    Parameters
    ----------
    data_dir  : Path – folder containing the Fluent export files
    label     : str  – human-readable name shown in plot title and console
    mfr_file  : str  – filename of the mass-flow-rate report
    y_file    : str  – filename of the outlet Y-coordinate report
    z_file    : str  – filename of the outlet Z-coordinate report
    pin_file  : str or None – filename of the inlet total pressure report;
                              if provided, avg_tot_pressure_out must also exist
    out_file  : Path or None – PNG output path; defaults to
                              data_dir/<slug>_analysis.png

    Returns
    -------
    dict with keys: entries, mdot, mean, cv, ui, delta_p
      delta_p is None when no pressure data is available (coarse mesh).
    """
    # Load reports
    mfr_raw = parse_fluent_report(data_dir / mfr_file)
    y_raw   = parse_fluent_report(data_dir / y_file)
    z_raw   = parse_fluent_report(data_dir / z_file)
    entries = sort_outlets_by_angle(mfr_raw, y_raw, z_raw)

    mdot  = np.array([e["mfr"] for e in entries])
    mean  = np.mean(mdot)
    cv    = np.std(mdot)  / mean * 100   # Coefficient of Variation
    max_d = (np.max(mdot) - np.min(mdot)) / mean * 100
    ui    = 1 - np.sum(np.abs(mdot - mean)) / (2 * len(mdot) * mean)  # Uniformity Index

    # Optional pressure drop (inlet total pressure − mean outlet total pressure)
    delta_p = None
    if pin_file and (data_dir / pin_file).exists():
        pin_raw  = parse_fluent_report(data_dir / pin_file)
        pout_raw = parse_fluent_report(data_dir / "avg_tot_pressure_out")
        p_in     = pin_raw["inlet"]
        p_out    = np.mean([pout_raw[e["orig_name"]]
                            for e in entries if e["orig_name"] in pout_raw])
        delta_p  = p_in - p_out

    # --- Console output ---
    print("\n=== {} ===".format(label))
    print("{:>4}  {:>14}  {:>8}  {:>8}  {:>12}".format(
        "No.", "Zone", "Y [mm]", "Z [mm]", "mdot [g/s]"))
    print("-" * 60)
    for i, e in enumerate(entries):
        print("{:>4}  {:>14}  {:>8.2f}  {:>8.2f}  {:>12.4f}".format(
            i + 1, e["orig_name"],
            e["y"] * 1000, e["z"] * 1000, e["mfr"] * 1000))
    print()
    print("  Outlets:          {}".format(len(entries)))
    print("  Mean:             {:.4f} g/s".format(mean * 1000))
    print("  CV:               {:.2f} %".format(cv))
    print("  Max deviation:    {:.2f} %".format(max_d))
    print("  Uniformity Index: {:.4f}".format(ui))
    if delta_p is not None:
        print("  Pressure drop:    {:.3f} kPa".format(delta_p / 1000))

    # --- Bar chart ---
    # Bars outside the ±THRESHOLD band are highlighted in red
    colors = ["tab:red" if abs(v - mean) / mean > THRESHOLD else "tab:blue"
              for v in mdot]

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(range(len(mdot)), mdot * 1000, color=colors, width=0.6)
    ax.axhline(mean * 1000, color="gray", linewidth=1, linestyle="--",
               label="Mean = {:.4f} g/s".format(mean * 1000))
    ax.axhline((1 + THRESHOLD) * mean * 1000, color="red",
               linewidth=0.8, linestyle=":",
               label="±{:g}% band".format(THRESHOLD * 100))
    ax.axhline((1 - THRESHOLD) * mean * 1000, color="red",
               linewidth=0.8, linestyle=":")
    ax.set_xticks(range(len(mdot)))
    ax.set_xticklabels([str(i + 1) for i in range(len(mdot))], fontsize=7)
    ax.set_xlabel("Outlet number – counter-clockwise starting from inlet")
    ax.set_ylabel("Mass flow rate [g/s]")
    ax.set_title("{} | CV = {:.2f}%  UI = {:.4f}".format(label, cv, ui))
    ax.legend()
    plt.tight_layout()

    if out_file is None:
        slug     = label.lower().replace(" ", "_").replace("–", "").replace("-", "")
        out_file = data_dir / "{}_analysis.png".format(slug)
    plt.savefig(str(out_file), dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Plot saved: {}".format(out_file))

    return {"entries": entries, "mdot": mdot, "mean": mean,
            "cv": cv, "ui": ui, "delta_p": delta_p}


# ---------------------------------------------------------------------------
# Coarse vs. refined mesh comparison
# ---------------------------------------------------------------------------

def compare_meshes(fine_result, coarse_result, out_file, geom_label):
    """
    Side-by-side bar chart comparing the coarse and refined mesh results
    for the same geometry.

    Uses up to n = min(len(fine), len(coarse)) channels so both datasets
    fit on the same x-axis even if outlet counts differ between meshes.
    """
    n = min(len(fine_result["mdot"]), len(coarse_result["mdot"]))
    mf = fine_result["mdot"][:n]
    mc = coarse_result["mdot"][:n]

    mean_f = np.mean(mf);  mean_c = np.mean(mc)
    cv_f   = np.std(mf) / mean_f * 100
    cv_c   = np.std(mc) / mean_c * 100
    ui_f   = 1 - np.sum(np.abs(mf - mean_f)) / (2 * n * mean_f)
    ui_c   = 1 - np.sum(np.abs(mc - mean_c)) / (2 * n * mean_c)

    print("\n=== Mesh comparison – {} ===".format(geom_label))
    print("  Fine mesh:   CV = {:.2f}%  UI = {:.4f}".format(cv_f, ui_f))
    print("  Coarse mesh: CV = {:.2f}%  UI = {:.4f}".format(cv_c, ui_c))

    x, w = np.arange(n), 0.4
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(x - w / 2, mf * 1000, w, label="Fine mesh",
           color="steelblue",  edgecolor="k", linewidth=0.2)
    ax.bar(x + w / 2, mc * 1000, w, label="Coarse mesh",
           color="darkorange", edgecolor="k", linewidth=0.2)
    ax.axhline(mean_f * 1000, color="steelblue",  linewidth=1, linestyle="--",
               label="Mean fine   = {:.4f} g/s".format(mean_f * 1000))
    ax.axhline(mean_c * 1000, color="darkorange", linewidth=1, linestyle="--",
               label="Mean coarse = {:.4f} g/s".format(mean_c * 1000))
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(n)], fontsize=7)
    ax.set_xlabel("Outlet number – counter-clockwise starting from inlet")
    ax.set_ylabel("Mass flow rate [g/s]")
    ax.set_title("{} – Mesh comparison | "
                 "Fine: CV={:.2f}% UI={:.4f} | "
                 "Coarse: CV={:.2f}% UI={:.4f}".format(
                     geom_label, cv_f, ui_f, cv_c, ui_c))
    ax.legend()
    plt.tight_layout()
    plt.savefig(str(out_file), dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Plot saved: {}".format(out_file))


# ---------------------------------------------------------------------------
# Geometry 1 vs. Geometry 2 comparison (refined meshes only)
# ---------------------------------------------------------------------------

def compare_geometries(geom1_result, geom2_result, out_file):
    """
    Side-by-side comparison of the two refined geometry variants.
    Called only when Geometry 2 data has been exported from Fluent.
    """
    n = min(len(geom1_result["mdot"]), len(geom2_result["mdot"]))
    m1 = geom1_result["mdot"][:n]
    m2 = geom2_result["mdot"][:n]

    mean1 = np.mean(m1);  mean2 = np.mean(m2)
    cv1   = np.std(m1) / mean1 * 100
    cv2   = np.std(m2) / mean2 * 100
    ui1   = 1 - np.sum(np.abs(m1 - mean1)) / (2 * n * mean1)
    ui2   = 1 - np.sum(np.abs(m2 - mean2)) / (2 * n * mean2)

    print("\n=== Geometry comparison (refined meshes) ===")
    print("  Geometry 1: CV = {:.2f}%  UI = {:.4f}".format(cv1, ui1))
    print("  Geometry 2: CV = {:.2f}%  UI = {:.4f}".format(cv2, ui2))

    x, w = np.arange(n), 0.4
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(x - w / 2, m1 * 1000, w, label="Geometry 1 (refined)",
           color="steelblue", edgecolor="k", linewidth=0.2)
    ax.bar(x + w / 2, m2 * 1000, w, label="Geometry 2 (refined)",
           color="seagreen",  edgecolor="k", linewidth=0.2)
    ax.axhline(mean1 * 1000, color="steelblue", linewidth=1, linestyle="--",
               label="Mean Geom 1 = {:.4f} g/s".format(mean1 * 1000))
    ax.axhline(mean2 * 1000, color="seagreen",  linewidth=1, linestyle="--",
               label="Mean Geom 2 = {:.4f} g/s".format(mean2 * 1000))
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(n)], fontsize=7)
    ax.set_xlabel("Outlet number – counter-clockwise starting from inlet")
    ax.set_ylabel("Mass flow rate [g/s]")
    ax.set_title("Geometry comparison (refined) | "
                 "Geom 1: CV={:.2f}% UI={:.4f} | "
                 "Geom 2: CV={:.2f}% UI={:.4f}".format(cv1, ui1, cv2, ui2))
    ax.legend()
    plt.tight_layout()
    plt.savefig(str(out_file), dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Plot saved: {}".format(out_file))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # --- Geometry 1: coarse mesh ---
    # Data files live in CFD_Export/, plots are written to plots/
    coarse_result = analyse_variant(
        data_dir = GEOM1_COARSE / "CFD_Export",
        label    = "Geometry 1 - Coarse mesh",
        mfr_file = "flow_rates",
        y_file   = "outlets_y",
        z_file   = "outlets_z",
        out_file = GEOM1_COARSE / "plots" / "geom1_coarse_analysis.png",
    )

    # --- Geometry 1: refined mesh ---
    fine_result = analyse_variant(
        data_dir = GEOM1_FINE / "CFD_Export",
        label    = "Geometry 1 - Refined mesh",
        mfr_file = "mass_flow_rates",
        y_file   = "y_coord",
        z_file   = "z_coord",
        pin_file = "pressure_in",
        out_file = GEOM1_FINE / "plots" / "geom1_fine_analysis.png",
    )

    # --- Geometry 1: coarse vs. refined comparison ---
    compare_meshes(
        fine_result   = fine_result,
        coarse_result = coarse_result,
        out_file      = GEOM1_DIR / "plots" / "geom1_mesh_comparison.png",
        geom_label    = "Geometry 1",
    )

    # --- Geometry 2: only processed when Fluent exports are present ---
    if (GEOM2_FINE / "CFD_Export" / "mass_flow_rates").exists():
        geom2_result = analyse_variant(
            data_dir = GEOM2_FINE / "CFD_Export",
            label    = "Geometry 2 - Refined mesh",
            mfr_file = "mass_flow_rates",
            y_file   = "y_coord",
            z_file   = "z_coord",
            pin_file = "pressure_in",
            out_file = GEOM2_FINE / "plots" / "geom2_fine_analysis.png",
        )
        compare_geometries(
            geom1_result = fine_result,
            geom2_result = geom2_result,
            out_file     = ROOT / "geometry_comparison.png",
        )
    else:
        print("\n[INFO] Geometry 2 data not found - Geometry 2 analysis skipped.")
        print("       Add Fluent exports to: {}".format(GEOM2_FINE / "CFD_Export"))

    print("\nAll done.")
