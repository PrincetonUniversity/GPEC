#!/usr/bin/env python
"""
Plot reconstructed covariant C components in Fourier space.

This script reads:
  - gpec_recon_cvmn_sol*.out
  - gpec_recon_cv2mn_sol*.out

and writes:
  - gpec_recon_cvmn.png
  - gpec_recon_cv2mn.png

Each PNG contains 6 subplots:
  cvp_re, cvp_im, cvt_re, cvt_im, cvz_re, cvz_im
with psi on the x-axis and one line per m value.
"""

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def find_required_file(pattern):
    matches = sorted(glob.glob(pattern))
    if not matches:
        print(f"Error: Could not find {pattern}")
        sys.exit(1)
    return matches[0]


def load_mn_data(path):
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data


def build_series(data):
    psi_vals = np.unique(data[:, 0])
    m_vals = np.unique(data[:, 1].astype(int))
    series = {}
    for m in m_vals:
        mask = data[:, 1].astype(int) == m
        block = data[mask]
        order = np.argsort(block[:, 0])
        series[m] = block[order]
    return psi_vals, m_vals, series


def style_axis(ax, title):
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("psi", fontsize=10)
    ax.grid(True, alpha=0.25, color="0.75", linewidth=0.6)
    ax.set_facecolor("#fbfbfb")


def plot_family(data, titles, png_name, figure_title):
    _, m_vals, series = build_series(data)
    cmap = plt.colormaps["nipy_spectral"]
    linestyles = ["-", "--", "-.", ":"]

    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
    axes = axes.ravel()

    for idx, ax in enumerate(axes):
        col = idx + 2
        for i, m in enumerate(m_vals):
            block = series[m]
            color = cmap((i + 0.5) / max(len(m_vals), 1))
            ax.plot(
                block[:, 0],
                block[:, col],
                color=color,
                linewidth=1.6 if i % 5 else 2.1,
                linestyle=linestyles[i % len(linestyles)],
                alpha=0.95,
            )
        style_axis(ax, titles[idx])

    axes[0].set_ylabel("value", fontsize=10)
    axes[1].set_ylabel("value", fontsize=10)
    axes[2].set_ylabel("value", fontsize=10)
    axes[3].set_ylabel("value", fontsize=10)
    axes[4].set_ylabel("value", fontsize=10)
    axes[5].set_ylabel("value", fontsize=10)

    # Add a compact colorbar-like legend using m-index colors.
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=m_vals.min(), vmax=m_vals.max())
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.02, pad=0.02)
    cbar.set_label("m", fontsize=10)

    fig.suptitle(figure_title, fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.96, 0.97])
    fig.savefig(png_name, dpi=160, bbox_inches="tight")
    plt.close(fig)


def sort_rows(data):
    return data[np.lexsort((data[:, 1], data[:, 0]))]


def symmetric_relative_l2(lhs, rhs):
    denom = np.linalg.norm(lhs) + np.linalg.norm(rhs)
    if denom == 0.0:
        return 0.0
    return 2.0 * np.linalg.norm(lhs - rhs) / denom


def max_pointwise_relative(lhs, rhs):
    scale = np.maximum(np.maximum(np.abs(lhs), np.abs(rhs)), 1.0e-14)
    return np.max(np.abs(lhs - rhs) / scale)


def compare_families(cvmn_data, cv2mn_data):
    lhs = sort_rows(cvmn_data)
    rhs = sort_rows(cv2mn_data)

    if lhs.shape != rhs.shape:
        raise ValueError(
            f"Shape mismatch: cvmn {lhs.shape} vs cv2mn {rhs.shape}"
        )

    if not np.allclose(lhs[:, :2], rhs[:, :2], atol=0.0, rtol=0.0):
        raise ValueError("psi/m grids do not match between cvmn and cv2mn.")

    labels = ["p", "t", "z"]
    results = []
    all_lhs = []
    all_rhs = []

    for i, label in enumerate(labels):
        col = 2 + 2 * i
        lhs_comp = lhs[:, col] + 1j * lhs[:, col + 1]
        rhs_comp = rhs[:, col] + 1j * rhs[:, col + 1]
        all_lhs.append(lhs_comp)
        all_rhs.append(rhs_comp)
        results.append(
            (
                label,
                symmetric_relative_l2(lhs_comp, rhs_comp),
                max_pointwise_relative(lhs_comp, rhs_comp),
            )
        )

    all_lhs = np.concatenate(all_lhs)
    all_rhs = np.concatenate(all_rhs)
    overall = (
        symmetric_relative_l2(all_lhs, all_rhs),
        max_pointwise_relative(all_lhs, all_rhs),
    )
    return results, overall


def main():
    cvmn_file = find_required_file("gpec_recon_cvmn_sol*.out")
    cv2mn_file = find_required_file("gpec_recon_cv2mn_sol*.out")

    print(f"Reading {os.path.basename(cvmn_file)}")
    print(f"Reading {os.path.basename(cv2mn_file)}")

    cvmn_data = load_mn_data(cvmn_file)
    cv2mn_data = load_mn_data(cv2mn_file)
    component_results, overall = compare_families(cvmn_data, cv2mn_data)

    plot_family(
        cvmn_data,
        [
            "cvp_re",
            "cvp_im",
            "cvt_re",
            "cvt_im",
            "cvz_re",
            "cvz_im",
        ],
        "gpec_recon_cvmn.png",
        "Covariant C from direct route (mn)",
    )

    plot_family(
        cv2mn_data,
        [
            "cv2p_re",
            "cv2p_im",
            "cv2t_re",
            "cv2t_im",
            "cv2z_re",
            "cv2z_im",
        ],
        "gpec_recon_cv2mn.png",
        "Covariant C from metric route (mn)",
    )

    print("Saved:")
    print("  gpec_recon_cvmn.png")
    print("  gpec_recon_cv2mn.png")
    print("")
    print("Symmetric relative error between cvmn and cv2mn:")
    for label, rel_l2, rel_max in component_results:
        print(
            f"  cv{label}: rel_L2 = {100.0 * rel_l2:8.3f}%"
            f"   max_pointwise = {100.0 * rel_max:8.3f}%"
        )
    print(
        f"  overall: rel_L2 = {100.0 * overall[0]:8.3f}%"
        f"   max_pointwise = {100.0 * overall[1]:8.3f}%"
    )


if __name__ == "__main__":
    main()
