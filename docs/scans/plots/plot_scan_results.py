#!/usr/bin/env python3
"""
Plot STRIDE scan results: Delta' and delta-W vs scan parameter.

Reads the summary CSVs produced by run_stride_beta_scan.py and
run_stride_epsilon_scan.py.

Usage:
    python plot_scan_results.py --scan beta
    python plot_scan_results.py --scan epsilon
    python plot_scan_results.py --scan both
    python plot_scan_results.py --scan both --save
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SCANS_DIR = SCRIPT_DIR.parent
OUTPUT_DIR = SCANS_DIR / 'outputs'

CSV_PATHS = {
    'beta': OUTPUT_DIR / 'beta_scan' / 'pressure_factor_scan_summary.csv',
    'epsilon': OUTPUT_DIR / 'epsilon_scan' / 'epsilon_scan_summary.csv',
}

PARAM_LABELS = {
    'beta': r'Pressure factor',
    'epsilon': r'$\epsilon = a/R_0$',
}


def load_csv(scan_type):
    """Load a scan summary CSV, returning a dict of arrays."""
    path = CSV_PATHS[scan_type]
    if not path.exists():
        print(f'ERROR: CSV not found: {path}')
        print(f'  Run the scan first: python run_stride_{scan_type}_scan.py')
        return None
    data = np.genfromtxt(str(path), delimiter=',', names=True)
    return data


def plot_single_scan(data, scan_type, save=False):
    """
    Plot Delta' and delta-W for a single scan.

    3-panel figure:
      - Delta' (real, log scale) for 2/1 and 3/1 modes
      - delta-W components vs parameter
      - q_edge vs parameter
    """
    if scan_type == 'beta':
        x = data['pressure_factor']
    else:
        x = data['epsilon']

    xlabel = PARAM_LABELS[scan_type]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: Delta' (real part, log scale) for both modes
    ax = axes[0]
    dp21 = data['delta_prime_21_real']
    dp31 = data['delta_prime_31_real']
    ax.semilogy(x, np.abs(dp21), 'o-', color='C0', markersize=3, lw=1.2, label='2/1')
    ax.semilogy(x, np.abs(dp31), 's-', color='C1', markersize=3, lw=1.2, label='3/1')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$|\Delta'|$")
    ax.set_title(r"$\Delta'$ (real)")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3, which='both')

    # Panel 2: delta-W components
    ax = axes[1]
    ax.plot(x, data['delta_W_plasma'], 's-', color='C1', markersize=3, lw=1.2, label=r'$\delta W_p$')
    ax.plot(x, data['delta_W_vacuum'], '^-', color='C2', markersize=3, lw=1.2, label=r'$\delta W_v$')
    ax.plot(x, data['delta_W_total'], 'o-', color='C3', markersize=3, lw=1.2, label=r'$\delta W_t$')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\delta W$')
    ax.set_title('Energy components')
    ax.axhline(0, color='k', lw=0.5, ls='--')
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)

    # Panel 3: edge q
    ax = axes[2]
    ax.plot(x, data['q_edge'], 'D-', color='C4', markersize=3, lw=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$q_{\mathrm{edge}}$')
    ax.set_title('Edge safety factor')
    ax.grid(True, alpha=0.3)

    scan_name = 'Pressure Factor' if scan_type == 'beta' else 'Epsilon'
    fig.suptitle(f'STRIDE Benchmark — {scan_name} Scan', fontsize=14, fontweight='bold')
    fig.tight_layout()

    if save:
        outfile = SCRIPT_DIR / f'{scan_type}_scan_results.png'
        fig.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f'Saved: {outfile}')


def plot_both(data_beta, data_epsilon, save=False):
    """
    Side-by-side comparison of both scans: Delta' and delta-W.

    2x2 figure: top row = Delta' (log), bottom row = delta-W.
    Left column = pressure_factor scan, right column = epsilon scan.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    scans = [
        ('beta', data_beta, data_beta['pressure_factor']),
        ('epsilon', data_epsilon, data_epsilon['epsilon']),
    ]

    for col, (scan_type, data, x) in enumerate(scans):
        xlabel = PARAM_LABELS[scan_type]

        # Top: Delta' (log scale, both modes)
        ax = axes[0, col]
        dp21 = data['delta_prime_21_real']
        dp31 = data['delta_prime_31_real']
        ax.semilogy(x, np.abs(dp21), 'o-', color='C0', markersize=3, lw=1.2, label='2/1')
        ax.semilogy(x, np.abs(dp31), 's-', color='C1', markersize=3, lw=1.2, label='3/1')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$|\Delta'|$")
        ax.set_title(r"$\Delta'$ (real)")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3, which='both')

        # Bottom: delta-W
        ax = axes[1, col]
        ax.plot(x, data['delta_W_plasma'], 's-', color='C1', markersize=3, lw=1.2, label=r'$\delta W_p$')
        ax.plot(x, data['delta_W_vacuum'], '^-', color='C2', markersize=3, lw=1.2, label=r'$\delta W_v$')
        ax.plot(x, data['delta_W_total'], 'o-', color='C3', markersize=3, lw=1.2, label=r'$\delta W_t$')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\delta W$')
        ax.set_title('Energy components')
        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3)

    fig.suptitle('STRIDE Benchmark Scans', fontsize=14, fontweight='bold')
    fig.tight_layout()

    if save:
        outfile = SCRIPT_DIR / 'both_scans_results.png'
        fig.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f'Saved: {outfile}')


def main():
    parser = argparse.ArgumentParser(
        description='Plot STRIDE scan results (Delta\' and delta-W)')
    parser.add_argument(
        '--scan', choices=['beta', 'epsilon', 'both'], default='both',
        help='Which scan to plot (default: both)')
    parser.add_argument(
        '--save', action='store_true',
        help='Save figures as PNG files')
    parser.add_argument(
        '--no-show', action='store_true',
        help='Do not display figures interactively')
    args = parser.parse_args()

    if args.scan in ('beta', 'both'):
        data_beta = load_csv('beta')
        if data_beta is None and args.scan == 'beta':
            sys.exit(1)
    else:
        data_beta = None

    if args.scan in ('epsilon', 'both'):
        data_epsilon = load_csv('epsilon')
        if data_epsilon is None and args.scan == 'epsilon':
            sys.exit(1)
    else:
        data_epsilon = None

    if args.scan == 'beta' and data_beta is not None:
        plot_single_scan(data_beta, 'beta', save=args.save)
    elif args.scan == 'epsilon' and data_epsilon is not None:
        plot_single_scan(data_epsilon, 'epsilon', save=args.save)
    elif args.scan == 'both':
        if data_beta is not None:
            plot_single_scan(data_beta, 'beta', save=args.save)
        if data_epsilon is not None:
            plot_single_scan(data_epsilon, 'epsilon', save=args.save)
        if data_beta is not None and data_epsilon is not None:
            plot_both(data_beta, data_epsilon, save=args.save)

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
