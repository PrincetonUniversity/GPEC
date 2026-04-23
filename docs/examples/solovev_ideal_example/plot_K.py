#!/usr/bin/env python3
"""
Plot K diagnostics on the r-z plane.

Running this script produces two figures:
1. gpec_recon_k_sol1.png
2. gpec_recon_k-term_sol1.png
"""

import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


def find_required_file(pattern):
    """Find a reconstruction output file by wildcard pattern."""
    matches = sorted(glob.glob(pattern))
    if not matches:
        print(f"Error: Could not find {pattern}")
        sys.exit(1)
    return matches[0]


def extract_suffix(path):
    """Extract the solution suffix such as sol1 from the filename."""
    name = os.path.basename(path)
    match = re.search(r"(sol\d+)\.out$", name)
    if not match:
        print(f"Error: Could not extract solution suffix from {name}")
        sys.exit(1)
    return match.group(1)


def make_norm(values):
    """Build a diverging normalization centered at zero."""
    vmin_temp, vmax_temp = np.percentile(values, [2, 98])
    vmax = max(abs(vmin_temp), abs(vmax_temp))
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def style_axis(ax, title):
    """Apply common axis styling."""
    ax.set_xlabel("R (m)", fontsize=11)
    ax.set_ylabel("Z (m)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")


def scatter_panel(ax, r, z, values, title, cbar_label):
    """Create a single scatter panel with symmetric color scaling."""
    scatter = ax.scatter(
        r,
        z,
        c=values,
        cmap="RdBu_r",
        s=20,
        alpha=0.6,
        norm=make_norm(values),
    )
    style_axis(ax, title)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(cbar_label, fontsize=10)


def main():
    """Main plotting routine."""
    print("=" * 70)
    print("GPEC K Diagnostic Plotter")
    print("=" * 70)

    k_file = find_required_file("gpec_recon_k_sol*.out")
    suffix = extract_suffix(k_file)

    print(f"Reading K data from: {k_file}")

    k_data = np.loadtxt(k_file, skiprows=1)

    r = k_data[:, 2]
    z = k_data[:, 3]

    k_total = k_data[:, 4]
    t1 = k_data[:, 6]
    t2 = k_data[:, 7]
    t3 = k_data[:, 8]

    # Figure 1: total K
    fig_k, ax_k = plt.subplots(figsize=(7, 8))
    scatter_panel(ax_k, r, z, k_total, "K (Total)", "K")
    plt.tight_layout()
    k_png = f"gpec_recon_k_{suffix}.png"
    fig_k.savefig(k_png, dpi=150, bbox_inches="tight")

    # Figure 2: three terms in one row
    fig_terms, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    term_specs = [
        (t1, r"Term1: $|\nabla\psi_{\mathrm{dcon}}|^2 \sigma S_{\mathrm{dcon}}$", "T1"),
        (t2, r"Term2: $B^2 \sigma^2$", "T2"),
        (t3, r"Term3: $2 P' \kappa^\psi$", "T3"),
    ]

    for ax, (values, title, label) in zip(axes, term_specs):
        scatter_panel(ax, r, z, values, title, label)

    plt.tight_layout()
    terms_png = f"gpec_recon_k-term_{suffix}.png"
    fig_terms.savefig(terms_png, dpi=150, bbox_inches="tight")

    print("Saved:")
    print(f"  {k_png}")
    print(f"  {terms_png}")


if __name__ == "__main__":
    main()
