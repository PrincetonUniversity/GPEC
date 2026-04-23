#!/usr/bin/env python3
"""
Plot curvature diagnostic quantities from GPEC reconstruction output.

Reads spatial domain data:
- gpec_recon_curvature_sol*.out: Spatial domain (r-z plane)

Produces single r-z plane plot with curvature values colored.
"""

import glob
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import griddata


def find_curvature_file():
    """Find curvature output file with wildcard matching."""
    curv_files = glob.glob("gpec_recon_curvature_sol*.out")

    if not curv_files:
        print("Error: Could not find gpec_recon_curvature_sol*.out file")
        sys.exit(1)

    return curv_files[0]


def read_curvature_data():
    """Read spatial curvature data from output file."""
    curv_file = find_curvature_file()

    print(f"Reading spatial data from: {curv_file}")
    curv = np.loadtxt(curv_file, skiprows=1)

    psi = curv[:, 0]
    theta = curv[:, 1]
    r = curv[:, 2]
    z = curv[:, 3]
    curv_re = curv[:, 4]
    curv_im = curv[:, 5]

    return psi, theta, r, z, curv_re, curv_im


def main():
    """Main plotting routine."""
    print("=" * 70)
    print("GPEC Curvature Diagnostic Plotter (r-z plane only)")
    print("=" * 70)

    psi, theta, r, z, curv_re, curv_im = read_curvature_data()

    unique_psi = np.unique(psi)

    print(f"Spatial grid: {len(unique_psi)} psi levels")
    print(f"Total grid points: {len(r)}")
    print(
        f"Curvature range: min={np.min(curv_re):.4e}, max={np.max(curv_re):.4e}, "
        f"imag_max={np.max(np.abs(curv_im)):.4e}"
    )

    fig, ax = plt.subplots(figsize=(10, 8))

    # Use the same compression used by the shear plot so sign structure is visible.
    curv_re_asinh = np.arcsinh(curv_re)

    r_min, r_max = np.min(r), np.max(r)
    z_min, z_max = np.min(z), np.max(z)

    grid_r = np.linspace(r_min, r_max, 150)
    grid_z = np.linspace(z_min, z_max, 150)
    grid_r_mesh, grid_z_mesh = np.meshgrid(grid_r, grid_z)

    points = np.c_[r, z]
    grid_curv = griddata(
        points, curv_re_asinh, (grid_r_mesh, grid_z_mesh),
        method="cubic", fill_value=np.nan
    )

    vmin, vmax = np.nanmin(grid_curv), np.nanmax(grid_curv)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    ax.contourf(
        grid_r_mesh, grid_z_mesh, grid_curv, levels=30,
        cmap="RdBu_r", norm=norm
    )

    scatter = ax.scatter(
        r, z, c=curv_re_asinh, cmap="RdBu_r",
        norm=norm, s=15, alpha=0.6, edgecolors="none"
    )

    ax.contour(
        grid_r_mesh, grid_z_mesh, grid_curv, levels=[0],
        colors="black", linewidths=2
    )

    ax.set_xlabel("R", fontsize=22, fontweight="bold")
    ax.set_ylabel("Z", fontsize=22, fontweight="bold")
    ax.set_title("Curvature (spatial) - r-z plane", fontsize=24,
                 fontweight="bold")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.tick_params(labelsize=16)

    r_ticks = np.linspace(r_min, r_max, 6)
    ax.set_xticks(r_ticks)
    ax.set_xticklabels([f"{rv:.3f}" for rv in r_ticks], fontsize=14)

    z_ticks = np.linspace(z_min, z_max, 6)
    ax.set_yticks(z_ticks)
    ax.set_yticklabels([f"{zv:.3f}" for zv in z_ticks], fontsize=14)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("asinh(Curvature)", fontsize=18, fontweight="bold")
    cbar.ax.tick_params(labelsize=14)

    cbar_ticks = np.linspace(vmin, vmax, 7)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{val:.1e}" for val in cbar_ticks], fontsize=12)

    plt.tight_layout()

    output_file = "gpec_recon_curvature_rz.png"
    print(f"\nSaving plot to: {output_file}")
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print("Plot complete!")


if __name__ == "__main__":
    main()
