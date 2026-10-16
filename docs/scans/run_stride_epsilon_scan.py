#!/usr/bin/env python3
"""
Run STRIDE epsilon (inverse aspect ratio) scan for LAR tokamak benchmark.

Iterates over a range of epsilon = a/R0 values, running STRIDE for each
equilibrium to compute delta-prime and delta-W stability quantities.

Usage:
    python run_stride_epsilon_scan.py
    python run_stride_epsilon_scan.py --gpechome /path/to/GPEC --test
    python run_stride_epsilon_scan.py --dry-run
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# ============================================================================
# Scan parameters
# ============================================================================

# Downsampled from 109 original points: NaNs removed, cropped before the
# Delta' pole at eps~0.662, thinned in flat regions, full resolution kept
# in the steep gradient region (eps > 0.62).
EPSILONS = [
    # Flat region (every ~3rd point)
    0.125, 0.1499, 0.1748, 0.1997, 0.2246, 0.2495, 0.2744, 0.2993,
    0.3242, 0.3491,
    # Moderate decline (every ~2nd point)
    0.3574, 0.374, 0.3906, 0.4072, 0.4238, 0.4404, 0.457, 0.4736,
    0.4902, 0.5005, 0.5151, 0.5317, 0.54275,
    # Slow approach (every ~2nd point)
    0.551, 0.55475, 0.55925, 0.56475, 0.57025, 0.57575, 0.58125,
    0.58675, 0.59225, 0.59775, 0.60325, 0.60875, 0.61425, 0.61975,
    # Steep rise to pole (all points kept)
    0.6225, 0.62525, 0.628, 0.63075, 0.6335, 0.63625, 0.639,
    0.64175, 0.6445, 0.64725, 0.65, 0.65125, 0.65375, 0.655,
    0.65625, 0.6575, 0.65875, 0.66, 0.66125,
]

EPSILONS_TEST = [0.2495, 0.4072, 0.595]

SCRIPT_DIR = Path(__file__).resolve().parent
EQUILIBRIA_DIR = SCRIPT_DIR / 'equilibria'
NAMELISTS_DIR = SCRIPT_DIR / 'inputs' / 'epsilon_scan'
OUTPUT_DIR = SCRIPT_DIR / 'outputs' / 'epsilon_scan'


# ============================================================================
# Namelist utilities (same as beta scan)
# ============================================================================


def _format_nml_value(val):
    if isinstance(val, bool):
        return '.TRUE.' if val else '.FALSE.'
    elif isinstance(val, str):
        return f'"{val}"'
    elif isinstance(val, float):
        if val != 0 and (abs(val) >= 1e6 or abs(val) < 1e-3):
            return f'{val:.6E}'
        return str(val)
    elif isinstance(val, (list, tuple)):
        return ','.join(_format_nml_value(v) for v in val)
    return str(val)


def _parse_nml_value(val_str):
    val_str = val_str.strip()
    if val_str.lower() in ('.true.', 't', '.t.'):
        return True
    if val_str.lower() in ('.false.', 'f', '.f.'):
        return False
    if (val_str.startswith('"') and val_str.endswith('"')) or \
       (val_str.startswith("'") and val_str.endswith("'")):
        return val_str[1:-1]
    rmatch = re.match(r'^(\d+)\*(.+)$', val_str)
    if rmatch:
        return [_parse_nml_value(rmatch.group(2))] * int(rmatch.group(1))
    try:
        return int(val_str)
    except ValueError:
        pass
    try:
        return float(val_str)
    except ValueError:
        pass
    return val_str


def read_namelist(filepath):
    text = Path(filepath).read_text()
    lines = []
    for line in text.split('\n'):
        in_str = False
        qc = None
        result = []
        for ch in line:
            if in_str:
                result.append(ch)
                if ch == qc:
                    in_str = False
            elif ch in ('"', "'"):
                in_str = True
                qc = ch
                result.append(ch)
            elif ch == '!':
                break
            else:
                result.append(ch)
        lines.append(''.join(result))
    text = '\n'.join(lines)

    groups = {}
    for match in re.finditer(r'&(\w+)(.*?)/', text, re.DOTALL):
        gname = match.group(1)
        body = re.sub(r'\s+', ' ', match.group(2)).strip()
        variables = {}
        tokens = re.split(r'(\w+(?:\([^)]*\))?\s*=)', body)
        i = 1
        while i < len(tokens) - 1:
            key = re.sub(r'\(.*\)', '', tokens[i].strip().rstrip('=').strip())
            val = tokens[i + 1].strip().rstrip(',').strip()
            variables[key] = _parse_nml_value(val)
            i += 2
        groups[gname] = variables
    return groups


def write_namelist(groups, filepath):
    lines = []
    for gname, variables in groups.items():
        lines.append(f'&{gname}')
        for key, val in variables.items():
            lines.append(f'    {key}={_format_nml_value(val)}')
        lines.append('/')
        lines.append('')
    Path(filepath).write_text('\n'.join(lines))


def update_namelist(filepath, overrides):
    groups = read_namelist(filepath)
    for gname, variables in overrides.items():
        if gname not in groups:
            groups[gname] = {}
        groups[gname].update(variables)
    write_namelist(groups, filepath)


# ============================================================================
# STRIDE runner
# ============================================================================


def run_stride_single(gpec_home, epsilon, run_dir, dry_run=False):
    """
    Run STRIDE for a single epsilon value.

    Returns path to stride_output_n1.nc on success, None on failure.
    """
    gpec_home = Path(gpec_home)
    run_dir = Path(run_dir)

    J = round(epsilon, 6)
    geqdsk = EQUILIBRIA_DIR / f'TJ_epsilon_scan_{J}.geqdsk'
    if not geqdsk.exists():
        print(f'  WARNING: geqdsk not found for epsilon={J}: {geqdsk}')
        return None

    run_dir.mkdir(parents=True, exist_ok=True)
    for fname in ['equil.in', 'dcon.in', 'stride.in', 'vac.in']:
        shutil.copy2(NAMELISTS_DIR / fname, run_dir / fname)

    update_namelist(
        run_dir / 'equil.in',
        {'EQUIL_CONTROL': {'eq_filename': str(geqdsk.resolve())}}
    )

    if dry_run:
        print(f'  [DRY RUN] epsilon={J}: would run STRIDE in {run_dir}')
        return None

    stride_bin = gpec_home / 'bin' / 'stride'
    if not stride_bin.exists():
        print(f'ERROR: STRIDE binary not found: {stride_bin}')
        sys.exit(1)

    print(f'  Running STRIDE for epsilon={J}...')
    result = subprocess.run(
        [str(stride_bin)],
        cwd=str(run_dir),
        capture_output=True,
    )

    if result.returncode != 0:
        print(f'  WARNING: STRIDE failed for epsilon={J} (rc={result.returncode})')
        stderr = result.stderr.decode('utf-8', errors='replace')
        if stderr:
            print(f'    stderr: {stderr[:500]}')
        return None

    output_nc = run_dir / 'stride_output_n1.nc'
    if not output_nc.exists():
        print(f'  WARNING: STRIDE output not found for epsilon={J}')
        return None

    return output_nc


def collect_results(output_files, epsilons):
    """Read delta-prime and delta-W from STRIDE output NC files."""
    try:
        from netCDF4 import Dataset as NCDataset
    except ImportError:
        import h5py

    n = len(epsilons)
    results = {
        'epsilon': np.array(epsilons),
        'delta_prime_21_real': np.full(n, np.nan),
        'delta_prime_21_imag': np.full(n, np.nan),
        'delta_prime_31_real': np.full(n, np.nan),
        'delta_prime_31_imag': np.full(n, np.nan),
        'delta_W_plasma': np.full(n, np.nan),
        'delta_W_vacuum': np.full(n, np.nan),
        'delta_W_total': np.full(n, np.nan),
        'q_edge': np.full(n, np.nan),
    }

    for i, (eps, nc_path) in enumerate(zip(epsilons, output_files)):
        if nc_path is None or not Path(nc_path).exists():
            continue
        try:
            if 'NCDataset' in dir():
                ds = NCDataset(str(nc_path), 'r')
                dp = ds.variables['Delta_prime'][:]
                results['delta_prime_21_real'][i] = dp[0, 0, 0]
                results['delta_prime_21_imag'][i] = dp[1, 0, 0]
                if dp.shape[1] > 1:
                    results['delta_prime_31_real'][i] = dp[0, 1, 1]
                    results['delta_prime_31_imag'][i] = dp[1, 1, 1]
                if 'plasma1' in ds.ncattrs():
                    results['delta_W_plasma'][i] = ds.getncattr('plasma1')
                if 'vacuum1' in ds.ncattrs():
                    results['delta_W_vacuum'][i] = ds.getncattr('vacuum1')
                if 'total1' in ds.ncattrs():
                    results['delta_W_total'][i] = ds.getncattr('total1')
                q = ds.variables['q'][:]
                results['q_edge'][i] = q[-1]
                ds.close()
            else:
                f = h5py.File(str(nc_path), 'r')
                dp = f['Delta_prime'][:]
                results['delta_prime_21_real'][i] = dp[0, 0, 0]
                results['delta_prime_21_imag'][i] = dp[1, 0, 0]
                if dp.shape[1] > 1:
                    results['delta_prime_31_real'][i] = dp[0, 1, 1]
                    results['delta_prime_31_imag'][i] = dp[1, 1, 1]
                if 'plasma1' in f.attrs:
                    results['delta_W_plasma'][i] = f.attrs['plasma1']
                if 'vacuum1' in f.attrs:
                    results['delta_W_vacuum'][i] = f.attrs['vacuum1']
                if 'total1' in f.attrs:
                    results['delta_W_total'][i] = f.attrs['total1']
                q = f['q'][:]
                results['q_edge'][i] = q[-1]
                f.close()
        except Exception as e:
            print(f'  WARNING: Could not read results for epsilon={eps}: {e}')

    return results


def save_summary(results, output_path):
    """Save scan results as a CSV summary file."""
    output_path = Path(output_path)
    header = ('epsilon,delta_prime_21_real,delta_prime_21_imag,'
              'delta_prime_31_real,delta_prime_31_imag,'
              'delta_W_plasma,delta_W_vacuum,delta_W_total,q_edge')
    data = np.column_stack([
        results['epsilon'],
        results['delta_prime_21_real'],
        results['delta_prime_21_imag'],
        results['delta_prime_31_real'],
        results['delta_prime_31_imag'],
        results['delta_W_plasma'],
        results['delta_W_vacuum'],
        results['delta_W_total'],
        results['q_edge'],
    ])
    np.savetxt(str(output_path), data, header=header, delimiter=',',
               fmt='%.10e', comments='')
    print(f'Summary saved to {output_path}')


# ============================================================================
# Main
# ============================================================================


def main():
    parser = argparse.ArgumentParser(
        description='Run STRIDE epsilon scan for LAR tokamak benchmark')
    parser.add_argument(
        '--gpechome', default=os.environ.get('GPECHOME', str(SCRIPT_DIR.parents[1])),
        help='Path to GPEC installation (default: $GPECHOME or inferred)')
    parser.add_argument(
        '--test', action='store_true',
        help='Run only a small subset of epsilon values for testing')
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Validate setup without actually running STRIDE')
    parser.add_argument(
        '--keep-workdirs', action='store_true',
        help='Keep intermediate working directories')
    args = parser.parse_args()

    epsilons = EPSILONS_TEST if args.test else EPSILONS
    gpec_home = Path(args.gpechome).resolve()

    print(f'GPEC home: {gpec_home}')
    print(f'Running epsilon scan with {len(epsilons)} points')
    print(f'Output directory: {OUTPUT_DIR}')

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    output_files = []

    for i, eps in enumerate(epsilons):
        J = round(eps, 6)
        print(f'[{i + 1}/{len(epsilons)}] epsilon = {J}')

        if args.keep_workdirs:
            run_dir = OUTPUT_DIR / f'run_eps_{J}'
        else:
            run_dir = Path(tempfile.mkdtemp(prefix=f'stride_eps_{J}_'))

        nc_path = run_stride_single(gpec_home, eps, run_dir, dry_run=args.dry_run)

        if nc_path is not None:
            dest = OUTPUT_DIR / f'stride_output_eps_{J}.nc'
            shutil.copy2(nc_path, dest)
            output_files.append(dest)
        else:
            output_files.append(None)

        if not args.keep_workdirs and not args.dry_run:
            shutil.rmtree(run_dir, ignore_errors=True)

    if not args.dry_run:
        results = collect_results(output_files, epsilons)
        save_summary(results, OUTPUT_DIR / 'epsilon_scan_summary.csv')

    print('Done.')


if __name__ == '__main__':
    main()
