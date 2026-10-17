#!/usr/bin/env python3
"""Native DCON regression for an equilibrium with no n=3 rational surface.

Build with gfortran -O0 -fcheck=all, then execute this module. As in the other
regression/test modules, run(repo, loc) returns True on success.
"""
import argparse
import math
import os
from pathlib import Path
import struct
import subprocess


def run(repo=None, loc=None):
    repo = Path(repo or Path(__file__).resolve().parents[2]).resolve()
    loc = Path(loc or repo / 'regression/runs/dcon_no_rational_surface').resolve()
    loc.mkdir(parents=True, exist_ok=True)
    (loc / 'lar.in').write_text('''&lar_input
 ma=256, mtau=128, lar_r0=6.2, lar_a=.62, beta0=1e-6,
 q0=2.1666666666666665, p_pres=2, p_sig=100, sigma_type='default',
 zeroth=.true.
/
''')
    (loc / 'equil.in').write_text('''&equil_control
 eq_type='lar', eq_filename='lar.in', jac_type='boozer',
 grid_type='ldp', psilow=.005, psihigh=.99,
 mpsi=32, mtheta=64, newq0=0, etol=1e-7
/
&equil_output
 gse_flag=.false., bin_eq_2d=.true.
/
''')
    (loc / 'dcon.in').write_text('''&dcon_control
 bal_flag=.false., mat_flag=.true., ode_flag=.true.,
 vac_flag=.false., mer_flag=.false., sas_flag=.false.,
 psiedge=.99, qlow=1.02, qhigh=30.2, sing_start=0,
 reform_eq_with_psilim=.false., nn=3, delta_mlow=0,
 delta_mhigh=0, delta_mband=0, parallel_threads=1
/
&dcon_output
 bin_euler=.true., euler_stride=3, netcdf_out=.false.
/
''')
    euler = loc / 'euler.bin'
    if euler.exists():
        euler.unlink()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = subprocess.run([str(repo / 'bin/dcon')], cwd=loc, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=60, text=True)
    (loc / 'dcon.log').write_text(result.stdout)
    if result.returncode or 'Normal termination.' not in result.stdout:
        raise AssertionError('DCON did not finish; see ' + str(loc / 'dcon.log'))
    # Independent oracle: n*q spans no integer, read from the native carrier.
    # Fortran sequential records use four-byte markers with GNU defaults.
    with euler.open('rb') as stream:
        records = []
        for _ in range(7):
            length = struct.unpack('=i', stream.read(4))[0]
            records.append(stream.read(length))
            assert struct.unpack('=i', stream.read(4))[0] == length
    nn = struct.unpack('=5i2d', records[0])[2]
    qmin, qmax = struct.unpack_from('=21d', records[4])[16:18]
    assert nn == 3 and 2.16 < qmin <= qmax < 2.3, (nn, qmin, qmax)
    assert math.ceil(nn * qmin) > math.floor(nn * qmax), (nn, qmin, qmax)
    assert 'psi = 9.900E-01' in result.stdout, 'ODE did not reach the edge'
    print(f'PASS: n={nn}, q=[{qmin:.9f}, {qmax:.9f}], no rational surface; ODE reaches .99')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo')
    parser.add_argument('--loc')
    args = parser.parse_args()
    run(args.repo, args.loc)
