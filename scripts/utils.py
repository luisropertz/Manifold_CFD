"""
Shared utilities for Manifold CFD post-processing scripts.
"""

import math
from pathlib import Path

# Deviation threshold for colour-coding outlets (fraction of mean mass flow).
# Applied consistently in both the bar-chart and the circle plot.
THRESHOLD = 0.03   # 3 %


def parse_fluent_report(filepath):
    """
    Parse a single Ansys Fluent surface-report text file.

    Fluent writes data lines in the form:
        zone-name    <numeric value>
    interspersed with header/summary lines that contain keywords we skip.
    Returns {zone_name: float_value}.
    """
    SKIP_KEYWORDS = ["Net", "Flow Rate", "Coordinate", "Integral",
                     "Average", "---", '"']
    result = {}
    with open(filepath) as f:
        for line in f:
            if any(kw in line for kw in SKIP_KEYWORDS):
                continue
            parts = line.rsplit(None, 1)
            if len(parts) == 2:
                try:
                    result[parts[0].strip()] = float(parts[1])
                except ValueError:
                    pass
    return result


def sort_outlets_by_angle(mfr_raw, y_raw, z_raw):
    """
    Return outlets sorted counter-clockwise starting from the inlet (negative Z-axis).

    The manifold outlets lie in the YZ-plane.  Channel 1 is at the negative-Z
    position (closest to the inlet); subsequent channels go counter-clockwise
    when viewed from the positive X-direction.

    Only zones that appear in all three report dictionaries are included.

    Returns a list of dicts with keys:
        orig_name, y, z, angle_deg, sort_key, mfr
    """
    common = sorted(set(mfr_raw) & set(y_raw) & set(z_raw))
    entries = []
    for zone in common:
        y     = y_raw[zone]
        z     = z_raw[zone]
        angle = math.atan2(y, z)
        entries.append({
            "orig_name": zone,
            "y":         y,
            "z":         z,
            "angle_deg": math.degrees(angle),
            "sort_key":  (math.pi - angle) % (2 * math.pi),
            "mfr":       abs(mfr_raw[zone]),
        })
    entries.sort(key=lambda e: e["sort_key"])
    return entries
