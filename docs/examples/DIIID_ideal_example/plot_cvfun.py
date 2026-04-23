#!/usr/bin/env python
"""
Plot cvfun vs cv2fun on the R-Z plane and report fun-space errors.

Outputs:
  - gpec_recon_cvfun_solX.png
  - gpec_recon_cv2fun_solX.png
  - gpec_recon_cvdiff_solX.png
"""

import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


def find_required_file(pattern):
    matches = sorted(glob.glob(pattern))
    if not matches:
        print(f"Error: Could not find {pattern}")
        sys.exit(1)
    return matches[0]


def extract_suffix(path):
    name = os.path.basename(path)
    match = re.search(r"(sol\d+)\.out$", name)
    if not match:
        print(f"Error: Could not extract solution suffix from {name}")
        sys.exit(1)
    return match.group(1)


def load_data(path):
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data


def make_diverging_norm(values):
    vmin_temp, vmax_temp = np.percentile(values, [2, 98])
    vmax = max(abs(vmin_temp), abs(vmax_temp), 1.0e-14)
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def scatter_component(ax, r, z, values, title, cbar_label):
    scatter = ax.scatter(
        r,
        z,
        c=values,
        s=10,
        cmap="RdBu_r",
        alpha=0.75,
        edgecolors="none",
        norm=make_diverging_norm(values),
    )
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("R (m)", fontsize=10)
    ax.set_ylabel("Z (m)", fontsize=10)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(cbar_label, fontsize=9)
    return scatter


def make_six_panel_figure(r, z, values_by_comp, figure_title, output):
    fig, axes = plt.subplots(3, 2, figsize=(13, 15))
    panel_order = [
        ("p_re", r"$C_\psi$ Re", "Re"),
        ("p_im", r"$C_\psi$ Im", "Im"),
        ("t_re", r"$C_\theta$ Re", "Re"),
        ("t_im", r"$C_\theta$ Im", "Im"),
        ("z_re", r"$C_\zeta$ Re", "Re"),
        ("z_im", r"$C_\zeta$ Im", "Im"),
    ]
    for ax, (key, title, label) in zip(axes.ravel(), panel_order):
        scatter_component(
            ax,
            r,
            z,
            values_by_comp[key],
            f"{figure_title} {title}",
            label,
        )
    fig.suptitle(figure_title, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def symmetric_relative_l2(lhs, rhs):
    denom = np.linalg.norm(lhs) + np.linalg.norm(rhs)
    if denom == 0.0:
        return 0.0
    return 2.0 * np.linalg.norm(lhs - rhs) / denom


def component_arrays(data):
    return {
        "p_re": data[:, 4],
        "p_im": data[:, 5],
        "t_re": data[:, 6],
        "t_im": data[:, 7],
        "z_re": data[:, 8],
        "z_im": data[:, 9],
    }


def main():
    cv_file = find_required_file("gpec_recon_cvfun_sol*.out")
    cv2_file = find_required_file("gpec_recon_cv2fun_sol*.out")
    suffix = extract_suffix(cv_file)

    print(f"Reading {os.path.basename(cv_file)}")
    print(f"Reading {os.path.basename(cv2_file)}")

    cv_data = load_data(cv_file)
    cv2_data = load_data(cv2_file)

    if cv_data.shape != cv2_data.shape or not np.allclose(
        cv_data[:, :4], cv2_data[:, :4]
    ):
        print("Error: cvfun and cv2fun grids do not match.")
        sys.exit(1)

    r = cv_data[:, 2]
    z = cv_data[:, 3]
    cv = component_arrays(cv_data)
    cv2 = component_arrays(cv2_data)
    panel_keys = ["p_re", "p_im", "t_re", "t_im", "z_re", "z_im"]
    comp_labels = [("p", r"$C_\psi$"), ("t", r"$C_\theta$"), ("z", r"$C_\zeta$")]

    direct_values = {}
    metric_values = {}
    diff_values = {}

    print("")
    print("Fun-space error summary:")

    for comp, latex_name in comp_labels:
        for part in ["re", "im"]:
            key = f"{comp}_{part}"
            direct = cv[key]
            metric = cv2[key]
            diff = direct - metric
            direct_values[key] = direct
            metric_values[key] = metric
            diff_values[key] = diff

            rel_l2 = symmetric_relative_l2(direct, metric)
            peak = np.max(np.maximum(np.abs(direct), np.abs(metric)))
            mask = np.maximum(np.abs(direct), np.abs(metric)) > 0.05 * peak
            masked_rel_l2 = symmetric_relative_l2(direct[mask], metric[mask])
            abs_energy_frac = np.sum(diff[mask] ** 2) / max(np.sum(diff**2), 1.0e-30)
            idx = np.argmax(np.abs(diff))

            print(
                f"  cv{comp}_{part}: rel_L2 = {100*rel_l2:6.3f}%"
                f"   rel_L2(>5% peak) = {100*masked_rel_l2:6.3f}%"
                f"   mismatch energy in >5% peak region = {100*abs_energy_frac:6.2f}%"
            )
            print(
                f"        max |Δ| = {abs(diff[idx]):.6e}"
                f" at psi={cv_data[idx,0]:.6f}, theta={cv_data[idx,1]:.6f},"
                f" R={cv_data[idx,2]:.6f}, Z={cv_data[idx,3]:.6f}"
            )
            print(
                f"        direct = {direct[idx]:.6e}, metric = {metric[idx]:.6e}"
            )

    cv_png = f"gpec_recon_cvfun_{suffix}.png"
    cv2_png = f"gpec_recon_cv2fun_{suffix}.png"
    diff_png = f"gpec_recon_cvdiff_{suffix}.png"

    make_six_panel_figure(
        r, z, direct_values, "Direct cvfun components", cv_png
    )
    make_six_panel_figure(
        r, z, metric_values, "Metric cv2fun components", cv2_png
    )
    make_six_panel_figure(
        r, z, diff_values, "Difference components (direct - metric)", diff_png
    )

    print("")
    print("Saved:")
    print(f"  {cv_png}")
    print(f"  {cv2_png}")
    print(f"  {diff_png}")


if __name__ == "__main__":
    main()
