"""
Manifold Circle Plot
=====================
Renders a top-down view of all manifold outlets arranged on their actual
ring positions, colour-coded by deviation from the mean mass flow rate:

    red   -> more than THRESHOLD above mean
    green -> within +-THRESHOLD of mean
    blue  -> more than THRESHOLD below mean

The script automatically discovers every Geometry_* folder in the repository
root and generates one circle plot per geometry whose refined/ data is present.
When Geometry 2 (or later geometries) are not yet simulated, they are silently
skipped with an informational message.

Run from any working directory - all paths are resolved from this file's location.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.transforms import Affine2D
from pathlib import Path
from utils import THRESHOLD, parse_fluent_report, sort_outlets_by_angle

# Repository root is one level above the scripts/ folder
ROOT = Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# Circle plot
# ---------------------------------------------------------------------------

def make_circle_plot(refined_dir, geom_name):
    """
    Create and save the circle plot for one geometry's refined simulation.

    Each outlet is drawn as a small rectangle at its actual YZ position,
    rotated so the long axis points tangentially around the ring.
    Rectangle size is scaled to the ring circumference so outlets tile neatly.

    Parameters
    ----------
    refined_dir : Path - folder containing mass_flow_rates, y_coord, z_coord
    geom_name   : str  - used in the plot title and output filename
    """
    mfr_raw = parse_fluent_report(refined_dir / "CFD_Export" / "mass_flow_rates")
    y_raw   = parse_fluent_report(refined_dir / "CFD_Export" / "y_coord")
    z_raw   = parse_fluent_report(refined_dir / "CFD_Export" / "z_coord")

    data = sort_outlets_by_angle(mfr_raw, y_raw, z_raw)

    mdot = np.array([d["mfr"] for d in data])
    mean = np.mean(mdot)

    # Outlet rectangle dimensions derived from the outlet-ring geometry
    r_mean      = np.mean([math.sqrt(d["y"]**2 + d["z"]**2) for d in data])
    arc_per_out = 2 * math.pi * r_mean / len(data)   # arc length per slot
    rect_w      = arc_per_out * 0.75   # tangential width (75 % of slot arc)
    rect_h      = r_mean * 0.13        # radial depth (13 % of ring radius)

    fig, ax = plt.subplots(figsize=(11, 11))
    ax.set_aspect("equal")
    ax.axis("off")

    # Dashed reference circle at the mean outlet radius
    ax.add_patch(plt.Circle((0, 0), r_mean, fill=False,
                             color="lightgray", linewidth=0.8, linestyle="--"))

    for i, d in enumerate(data):
        # Map CFD coordinates to plot axes: CFD Y -> plot X, CFD Z -> plot Y
        cx = d["y"]
        cy = d["z"]

        deviation = (d["mfr"] - mean) / mean
        color = ("tab:red"   if deviation >  THRESHOLD else
                 "tab:blue"  if deviation < -THRESHOLD else
                 "tab:green")

        # Angle from positive plot-X axis (used for rotation)
        radial_angle_deg = math.degrees(math.atan2(cy, cx))

        # Rectangle centred at (0, 0), then rotated so its long side is
        # tangential, then translated to the outlet position
        rect = mpatches.Rectangle(
            (-rect_w / 2, -rect_h / 2), rect_w, rect_h,
            facecolor=color, edgecolor="black", linewidth=0.4,
        )
        t = (Affine2D()
             .rotate_deg(radial_angle_deg + 90)   # +90 deg makes the long axis tangential
             .translate(cx, cy)
             + ax.transData)
        rect.set_transform(t)
        ax.add_patch(rect)

        # Channel number printed just outside the rectangle
        label_r   = math.sqrt(cx**2 + cy**2) + rect_h * 0.7 + r_mean * 0.06
        angle_rad = math.atan2(cy, cx)
        ax.text(label_r * math.cos(angle_rad), label_r * math.sin(angle_rad),
                str(i + 1), ha="center", va="center", fontsize=6.5)

    # Inlet arrow - points at the gap at the negative-Z position
    ax.annotate(
        "Inlet",
        xy     = (0, -r_mean),
        xytext = (0, -r_mean - r_mean * 0.35),
        ha="center", va="top", fontsize=10, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="black", linewidth=0.8),
    )

    # Legend
    ax.legend(handles=[
        mpatches.Patch(color="tab:red",
                       label="> +{:.0f}% above mean".format(THRESHOLD * 100)),
        mpatches.Patch(color="tab:green",
                       label="within +/-{:.0f}%".format(THRESHOLD * 100)),
        mpatches.Patch(color="tab:blue",
                       label="< -{:.0f}% below mean".format(THRESHOLD * 100)),
    ], loc="upper right", fontsize=9)

    ax.set_title(
        "{} (refined) - Channel Distribution\n"
        "Mean: {:.4f} g/s  |  Threshold: +/-{:.0f}%".format(
            geom_name, mean * 1000, THRESHOLD * 100),
        fontsize=11, fontweight="bold",
    )

    margin = r_mean * 1.8
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin * 1.15, margin)
    plt.tight_layout()

    out = refined_dir / "plots" / "{}_circle.png".format(geom_name.lower().replace(" ", "_"))
    plt.savefig(str(out), dpi=150, bbox_inches="tight")
    plt.close()
    print("  -> Circle plot saved: {}".format(out))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Discover all Geometry_* folders sorted by name
    geom_dirs = sorted(ROOT.glob("Geometry_*"))

    if not geom_dirs:
        print("No Geometry_* folders found in {}.".format(ROOT))
    else:
        for geom_dir in geom_dirs:
            refined_dir = geom_dir / "refined"
            required_files = [
                refined_dir / "CFD_Export" / "mass_flow_rates",
                refined_dir / "CFD_Export" / "y_coord",
                refined_dir / "CFD_Export" / "z_coord",
            ]
            if all(p.exists() for p in required_files):
                print("Processing {}...".format(geom_dir.name))
                make_circle_plot(refined_dir, geom_dir.name)
            else:
                print("[SKIP] {} - refined data not yet available.".format(
                    geom_dir.name))

    print("\nAll done.")
