#!/usr/bin/env python3
"""
Plot overlaid equilibria and kinetic profiles for benchmark scans.

Produces multi-panel figures showing:
  - Flux surface contours (R,Z) colored by scan parameter
  - 1D profiles: q(psi_N), p(psi_N), f(psi_N)
  - Kinetic profiles from .gpeckf files: ne(psi_N), Te(psi_N)

Inspired by GPEC OMFIT module's plot_equilibrium.py and plot_equil_summary.py.

Usage:
    python plot_scan_equilibria.py --scan beta
    python plot_scan_equilibria.py --scan epsilon
    python plot_scan_equilibria.py --scan beta --stride 3  # plot every 3rd
    python plot_scan_equilibria.py --scan beta --save
"""

import argparse
import glob
import os
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

SCRIPT_DIR = Path(__file__).resolve().parent
SCANS_DIR = SCRIPT_DIR.parent
EQUILIBRIA_DIR = SCANS_DIR / 'equilibria'

# ============================================================================
# Lightweight geqdsk reader (no external dependencies)
# ============================================================================


def _read_line_values(f, n):
    """Read n values from a geqdsk file (5 per line, 16 chars each)."""
    vals = []
    while len(vals) < n:
        line = f.readline()
        # Each value is 16 characters wide, 5 per line
        for i in range(0, len(line.rstrip('\n')), 16):
            chunk = line[i:i + 16].strip()
            if chunk:
                vals.append(float(chunk))
            if len(vals) >= n:
                break
    return np.array(vals[:n])


def read_geqdsk(filepath):
    """
    Read a G-EQDSK file and return key equilibrium quantities.

    Returns dict with: fpol, pres, q, psi_norm, psi_grid, R, Z, rbnd, zbnd, etc.
    """
    with open(filepath, 'r') as f:
        # First line: identification string and grid dimensions
        line0 = f.readline()
        parts = line0.split()
        nw = int(parts[-2])  # radial grid points
        nh = int(parts[-1])  # vertical grid points

        # Line 2-5: 4 values per read, 5 lines total = 20 scalar values
        scalars = _read_line_values(f, 20)
        rdim, zdim, rcentr, rleft, zmid = scalars[0:5]
        rmag, zmag, simag, sibry, bcentr = scalars[5:10]
        current, simag2, xdum, rmag2, xdum2 = scalars[10:15]
        zmag2, xdum3, sibry2, xdum4, xdum5 = scalars[15:20]

        # 1D profiles on uniform psi grid (nw points)
        fpol = _read_line_values(f, nw)
        pres = _read_line_values(f, nw)
        ffprim = _read_line_values(f, nw)
        pprime = _read_line_values(f, nw)

        # 2D poloidal flux (nw x nh)
        psi_grid = _read_line_values(f, nw * nh).reshape(nh, nw)

        # Safety factor on uniform psi grid
        q = _read_line_values(f, nw)

        # Boundary and limiter
        line = f.readline()
        parts = line.split()
        nbbbs = int(parts[0])  # boundary points
        limitr = int(parts[1])  # limiter points

        rbnd = np.zeros(0)
        zbnd = np.zeros(0)
        if nbbbs > 0:
            bnd = _read_line_values(f, 2 * nbbbs)
            rbnd = bnd[0::2]
            zbnd = bnd[1::2]

    # Derived quantities
    psi_norm = np.linspace(0, 1, nw)
    R = np.linspace(rleft, rleft + rdim, nw)
    Z = np.linspace(zmid - zdim / 2, zmid + zdim / 2, nh)

    return {
        'fpol': fpol,
        'pres': pres,
        'q': q,
        'psi_norm': psi_norm,
        'psi_grid': psi_grid,
        'R': R,
        'Z': Z,
        'rbnd': rbnd,
        'zbnd': zbnd,
        'R0': rcentr,
        'Rmag': rmag,
        'Zmag': zmag,
        'Bt0': bcentr,
        'psi_axis': simag,
        'psi_bnd': sibry,
    }


# ============================================================================
# Kinetic file reader
# ============================================================================


def read_gpeckf(filepath):
    """Read a .gpeckf kinetic profiles file (6-column ASCII)."""
    data = np.loadtxt(filepath, skiprows=1)
    return {
        'psi_N': data[:, 0],
        'ni': data[:, 1],
        'ne': data[:, 2],
        'Ti': data[:, 3],
        'Te': data[:, 4],
        'wexb': data[:, 5],
    }


# ============================================================================
# File discovery
# ============================================================================


def discover_files(scan_type, stride=1):
    """
    Discover geqdsk and gpeckf files for a given scan type.

    Returns list of (param_value, geqdsk_path, gpeckf_path) tuples.
    """
    if scan_type == 'beta':
        geqdsk_pattern = 'TJ_betascan_*.geqdsk'
        prefix = 'TJ_betascan_'
    elif scan_type == 'epsilon':
        geqdsk_pattern = 'TJ_epsilon_scan_*.geqdsk'
        prefix = 'TJ_epsilon_scan_'
    else:
        raise ValueError(f"Unknown scan type: {scan_type}")

    files = sorted(glob.glob(str(EQUILIBRIA_DIR / geqdsk_pattern)))
    entries = []

    for geq_path in files:
        base = os.path.basename(geq_path)
        val_str = base.replace(prefix, '').replace('.geqdsk', '')
        try:
            val = float(val_str)
        except ValueError:
            continue

        gpeckf_path = geq_path.replace('.geqdsk', '.gpeckf')
        if not os.path.exists(gpeckf_path):
            gpeckf_path = None

        entries.append((val, geq_path, gpeckf_path))

    # Apply stride
    return entries[::stride]


# ============================================================================
# Plotting
# ============================================================================


def plot_equilibria(entries, scan_type, save=False):
    """
    Create multi-panel equilibrium overview figure.

    Panels:
    - [0,0] Overlaid flux surface boundaries
    - [0,1] q(psi_N) profiles
    - [1,0] p(psi_N) profiles
    - [1,1] f(psi_N) profiles
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ax_flux = axes[0, 0]
    ax_q = axes[0, 1]
    ax_p = axes[1, 0]
    ax_f = axes[1, 1]

    param_label = r'$\beta$' if scan_type == 'beta' else r'$\epsilon$'

    values = [e[0] for e in entries]
    norm = Normalize(vmin=min(values), vmax=max(values))
    cmap = cm.viridis

    for val, geq_path, _ in entries:
        color = cmap(norm(val))
        try:
            geq = read_geqdsk(geq_path)
        except Exception as e:
            print(f'  WARNING: Could not read {geq_path}: {e}')
            continue

        # Flux surface boundary
        if len(geq['rbnd']) > 0:
            ax_flux.plot(geq['rbnd'], geq['zbnd'], color=color, alpha=0.6, lw=0.8)

        # 1D profiles
        psi_N = geq['psi_norm']
        ax_q.plot(psi_N, np.abs(geq['q']), color=color, alpha=0.6, lw=0.8)
        ax_p.plot(psi_N, geq['pres'], color=color, alpha=0.6, lw=0.8)
        ax_f.plot(psi_N, np.abs(geq['fpol']), color=color, alpha=0.6, lw=0.8)

    # Formatting
    ax_flux.set_aspect('equal')
    ax_flux.set_xlabel('R (m)')
    ax_flux.set_ylabel('Z (m)')
    ax_flux.set_title('Flux Surface Boundaries')

    ax_q.set_xlabel(r'$\psi_N$')
    ax_q.set_ylabel('q')
    ax_q.set_title('Safety Factor')

    ax_p.set_xlabel(r'$\psi_N$')
    ax_p.set_ylabel('p (Pa)')
    ax_p.set_title('Pressure')

    ax_f.set_xlabel(r'$\psi_N$')
    ax_f.set_ylabel('f = RB_t')
    ax_f.set_title('Poloidal Current Function')

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02)
    cbar.set_label(param_label)

    fig.suptitle(f'Equilibrium Overview - {scan_type.capitalize()} Scan',
                 fontsize=14, fontweight='bold')
    fig.subplots_adjust(right=0.88, top=0.93)

    if save:
        outfile = SCRIPT_DIR / f'{scan_type}_scan_equilibria.png'
        fig.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f'Saved: {outfile}')


def plot_kinetic_profiles(entries, scan_type, save=False):
    """
    Plot overlaid kinetic profiles (ne, Te) from .gpeckf files.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax_ne = axes[0]
    ax_Te = axes[1]

    param_label = r'$\beta$' if scan_type == 'beta' else r'$\epsilon$'

    values = [e[0] for e in entries if e[2] is not None]
    if not values:
        print('No .gpeckf files found for kinetic profile plot.')
        return

    norm = Normalize(vmin=min(values), vmax=max(values))
    cmap = cm.viridis

    for val, _, gpeckf_path in entries:
        if gpeckf_path is None:
            continue
        color = cmap(norm(val))
        try:
            kin = read_gpeckf(gpeckf_path)
        except Exception as e:
            print(f'  WARNING: Could not read {gpeckf_path}: {e}')
            continue

        ax_ne.plot(kin['psi_N'], kin['ne'], color=color, alpha=0.6, lw=0.8)
        ax_Te.plot(kin['psi_N'], kin['Te'], color=color, alpha=0.6, lw=0.8)

    ax_ne.set_xlabel(r'$\psi_N$')
    ax_ne.set_ylabel(r'$n_e$ (m$^{-3}$)')
    ax_ne.set_title('Electron Density')

    ax_Te.set_xlabel(r'$\psi_N$')
    ax_Te.set_ylabel(r'$T_e$ (eV)')
    ax_Te.set_title('Electron Temperature')

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.8, pad=0.02)
    cbar.set_label(param_label)

    fig.suptitle(f'Kinetic Profiles - {scan_type.capitalize()} Scan',
                 fontsize=14, fontweight='bold')
    fig.subplots_adjust(right=0.88, top=0.90)

    if save:
        outfile = SCRIPT_DIR / f'{scan_type}_scan_kinetic.png'
        fig.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f'Saved: {outfile}')


# ============================================================================
# Main
# ============================================================================


def main():
    parser = argparse.ArgumentParser(
        description='Plot overlaid equilibria and kinetic profiles for benchmark scans')
    parser.add_argument(
        '--scan', choices=['beta', 'epsilon'], required=True,
        help='Which scan to plot')
    parser.add_argument(
        '--stride', type=int, default=1,
        help='Plot every Nth equilibrium (default: all)')
    parser.add_argument(
        '--save', action='store_true',
        help='Save figures as PNG files')
    parser.add_argument(
        '--no-show', action='store_true',
        help='Do not display figures interactively')
    args = parser.parse_args()

    entries = discover_files(args.scan, stride=args.stride)
    if not entries:
        print(f'No {args.scan} scan files found in {EQUILIBRIA_DIR}')
        sys.exit(1)

    print(f'Found {len(entries)} equilibria for {args.scan} scan')

    plot_equilibria(entries, args.scan, save=args.save)
    plot_kinetic_profiles(entries, args.scan, save=args.save)

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
