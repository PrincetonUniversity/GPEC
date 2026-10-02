"""Write a bounded synthetic circular source in native VACIN5 format."""
import math
from pathlib import Path
import sys

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
nodes, low, high, n = 128, -2, 2, 3
center, radius, zcenter = 6.2, 0.62, 0.0
theta = [i / nodes for i in range(nodes + 1)]
R = [center + radius * math.cos(2 * math.pi * t) for t in theta[:-1]]
# The native vacuum contour convention is clockwise in the (R,Z) plane.
Z = [zcenter - radius * math.sin(2 * math.pi * t) for t in theta[:-1]]
R.append(R[0])
Z.append(Z[0])
coefficients = {-1: 0.04 - 0.03j, 0: 0.02 + 0.07j, 2: -0.015 + 0.01j}
lines = ["scalars", "", "2", "2", f"{center-radius:.17e}",
         f"{center+radius:.17e}", f"{zcenter-radius:.17e}",
         f"{zcenter+radius:.17e}", str(nodes), str(low), str(high), str(n), "2.5"]


def block(label, values):
    lines.extend([label, ""])
    for i in range(0, len(values), 4):
        lines.append(" ".join(f"{v:.17e}" for v in values[i:i+4]))


block("Poloidal Coordinate Theta", theta)
block("Normalized Arclength", theta)
block("Major Radius", R)
block("Vertical Coordinate", Z)
block("Toroidal Angle Shift", [0.0] * (nodes + 1))
block("Normal Field Real", [coefficients.get(m, 0j).real for m in range(low, high + 1)])
block("Normal Field Imaginary", [coefficients.get(m, 0j).imag for m in range(low, high + 1)])
(out / "vacin5").write_text("\n".join(lines) + "\n")
(out / "query_contract.in").write_text(f"{nodes} {high-low+1}\n")
