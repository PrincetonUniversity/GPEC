#!/usr/bin/env python3
"""
Convert TJ Equilibrium .nc files to GPEC kinetic format (.gpeckf).

Reads ne (electron density) and Te (electron temperature) profiles from
TJ equilibrium netCDF files and writes them in the 6-column .kin format
expected by GPEC/PENTRC:

    psi_N  ni(m^-3)  ne(m^-3)  Ti(eV)  Te(eV)  omega_ExB(rad/s)

Assumptions for these analytic LAR equilibria:
  - ni = ne  (quasineutrality, single ion species)
  - Ti = Te  (equal temperatures)
  - omega_ExB = 0  (no rotation)

Usage:
    python convert_tj_to_gpeckf.py --tj-dir /path/to/TJ_equils --output-dir equilibria/
    python convert_tj_to_gpeckf.py  # uses defaults relative to this script
"""

import argparse
import glob
import os
import sys

import numpy as np

try:
    from netCDF4 import Dataset as NCDataset
    NC_BACKEND = 'netCDF4'
except ImportError:
    try:
        import h5py
        NC_BACKEND = 'h5py'
    except ImportError:
        print("ERROR: Requires netCDF4 or h5py to read TJ equilibrium files.")
        sys.exit(1)


def read_tj_equilibrium(nc_path):
    """Read PsiN, ne, Te from a TJ Equilibrium netCDF file."""
    if NC_BACKEND == 'netCDF4':
        ds = NCDataset(nc_path, 'r')
        psiN = ds.variables['PsiN'][:].copy()
        ne = ds.variables['ne'][:].copy()
        Te = ds.variables['Te'][:].copy()
        ds.close()
    else:
        f = h5py.File(nc_path, 'r')
        psiN = f['PsiN'][:].copy()
        ne = f['ne'][:].copy()
        Te = f['Te'][:].copy()
        f.close()
    return psiN, ne, Te


def write_gpeckf(output_path, psiN, ne, Te):
    """
    Write a .gpeckf kinetic profiles file.

    Format matches GPEC/PENTRC .kin reader (pentrc/inputs.f90 read_kin):
    6 columns: psi_N, ni(m^-3), ne(m^-3), Ti(eV), Te(eV), omega_ExB(rad/s)
    """
    ni = ne    # quasineutrality
    Ti = Te    # equal temperatures
    wexb = np.zeros_like(psiN)  # no rotation

    header = ("             psi         ni(m^-3)         ne(m^-3)"
              "           ti(eV)           te(eV)      wexb(rad/s)")

    with open(output_path, 'w') as f:
        f.write(header + '\n')
        for i in range(len(psiN)):
            f.write(f"  {psiN[i]:14.8e}   {ni[i]:14.8e}   {ne[i]:14.8e}"
                    f"   {Ti[i]:14.8e}   {Te[i]:14.8e}   {wexb[i]:14.8e}\n")


def nc_basename_to_gpeckf(nc_filename):
    """
    Convert Equilibrium NC filename to .gpeckf filename.

    Examples:
        Equilibrium_betascan_0.001.nc -> TJ_betascan_0.001.gpeckf
        Equilibrium_epsilon_scan_0.125.nc -> TJ_epsilon_scan_0.125.gpeckf
        Equilibrium_0.001.nc -> TJ_0.001.gpeckf
    """
    base = os.path.basename(nc_filename)
    # Remove 'Equilibrium_' prefix and '.nc' suffix
    stem = base.replace('Equilibrium_', '').replace('.nc', '')
    return f'TJ_{stem}.gpeckf'


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description='Convert TJ Equilibrium .nc files to GPEC .gpeckf format')
    parser.add_argument(
        '--tj-dir',
        default=os.path.join(os.path.dirname(script_dir),
                             '..', '..', '..',
                             'CTM-processing', 'SLAYER_benchmark_paper', 'TJ_equils'),
        help='Directory containing TJ Equilibrium_*.nc files')
    parser.add_argument(
        '--output-dir',
        default=os.path.join(script_dir, 'equilibria'),
        help='Output directory for .gpeckf files')
    parser.add_argument(
        '--pattern', default='Equilibrium_*.nc',
        help='Glob pattern for NC files (default: Equilibrium_*.nc)')
    parser.add_argument(
        '--dry-run', action='store_true',
        help='List files that would be converted without writing')
    args = parser.parse_args()

    tj_dir = os.path.abspath(args.tj_dir)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.isdir(tj_dir):
        print(f"ERROR: TJ directory not found: {tj_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    nc_files = sorted(glob.glob(os.path.join(tj_dir, args.pattern)))
    if not nc_files:
        print(f"No files matching '{args.pattern}' found in {tj_dir}")
        sys.exit(1)

    print(f"Found {len(nc_files)} Equilibrium NC files in {tj_dir}")
    print(f"Output directory: {output_dir}")

    n_converted = 0
    n_skipped = 0

    for nc_path in nc_files:
        gpeckf_name = nc_basename_to_gpeckf(nc_path)
        gpeckf_path = os.path.join(output_dir, gpeckf_name)

        if args.dry_run:
            print(f"  Would convert: {os.path.basename(nc_path)} -> {gpeckf_name}")
            n_converted += 1
            continue

        try:
            psiN, ne, Te = read_tj_equilibrium(nc_path)
            write_gpeckf(gpeckf_path, psiN, ne, Te)
            n_converted += 1
        except Exception as e:
            print(f"  WARNING: Failed to convert {os.path.basename(nc_path)}: {e}")
            n_skipped += 1

    action = "Would convert" if args.dry_run else "Converted"
    print(f"\n{action} {n_converted} files, skipped {n_skipped}")


if __name__ == '__main__':
    main()
