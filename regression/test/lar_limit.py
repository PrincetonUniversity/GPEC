#!/usr/bin/env python
"""

Regression Test - Test Large Aspect Ratio Limit Non-Ambipolar Transport Calculations
====================================================================================

:Date: 08/18/2016

Run notes:
    -

Result notes:
    -


"""

import os
import numpy as np
import xarray
import gpec, data

plt = data.plt

repo = '/'.join(__file__.split('/')[:-3])
loc = repo + '/regression/runs/lar_lim'


def run(repo=repo, loc=None):
    """
    Run the test cases.
    By default it tests this version in this repository, but you can use this to
    run the test for any repository (i.e. older versions that don't have this in pypec)
    and run it in any location (e.g. /p/gpec/benchmark on portal).

    Args:
        repo: str. Path to GPEC repository. Default is this file's repository.
        loc: str. Path where runs will be made. Default is repo/regression/runs.

    Returns:
        True.

    """
    # set up paths
    repo = os.path.abspath(repo)
    if loc is None:
        loc = repo + '/regression/runs'
    loc = os.path.abspath(loc)
    exedir = repo + '/bin'
    exmdir = '{:}/docs/examples/solovev_ideal_example/'.format(repo)
    # backwards comparability
    if not os.path.isdir(exmdir):
        exmdir = exmdir.replace('_ideal_example', '_example')

    # start from example input files so we don't need to maintain too many different sets
    inputs = {}
    for f in os.listdir(exmdir):
        if f.endswith('.in'):
            inputs[f[:-3]] = gpec.namelist.read(exmdir + f)

    # make any necessary modifications to input files

    # make input paths explicit if they aren't already
    pf = inputs['pentrc']['PENT_INPUT']['kinetic_file']
    if not pf.startswith('/'):
        inputs['pentrc']['PENT_INPUT']['kinetic_file'] = os.path.abspath(exmdir + pf)
    inputs['pentrc']['PENT_INPUT']['data_dir'] = repo + '/pentrc'

    # SOL - no rationals, no elongation, order unity lengths
    inputs['sol']['SOL_INPUT']['e'] = 1.0
    inputs['sol']['SOL_INPUT']['q0'] = 1.2
    inputs['sol']['SOL_INPUT']['r0'] = 1.0
    inputs['sol']['SOL_INPUT']['a'] = 1.0
    # VAC - close, stabilizing wall
    inputs['vac']['SHAPE']['a'] = 0.06
    # EQUIL - use example settings
    # DCON - use all example settings
    # GPEC - apply simple m=2 field and minimize output
    inputs['gpec']['GPEC_INPUT']['mode_flag'] = False
    inputs['gpec']['GPEC_INPUT']['harmonic_flag'] = True
    inputs['gpec']['GPEC_INPUT']['cosmn(2)'] = 1e-4
    inputs['gpec']['GPEC_INPUT']['sinmn(2)'] = -1e-4
    #inputs['gpec']['GPEC_CONTROL']['atol_psi'] = 1e-2
    # turn everything off
    for key, value in inputs['gpec']['GPEC_OUTPUT'].iteritems():
        if type(value) is bool:
            inputs['gpec']['GPEC_OUTPUT'][key] = False
    for key, value in inputs['gpec']['GPEC_DIAGNOSE'].iteritems():
        if type(value) is bool:
            inputs['gpec']['GPEC_DIAGNOSE'][key] = False
    # turn minimal output on
    inputs['gpec']['GPEC_OUTPUT']['xclebsch_flag'] = True
    inputs['gpec']['GPEC_OUTPUT']['verbose'] = True
    # PENTRC - compare RLAR, CLAR, and TGAR. Don't bother running anything else.
    for key, value in inputs['pentrc']['PENT_OUTPUT'].iteritems():
        if key.endswith('_flag'):
            inputs['pentrc']['PENT_OUTPUT'][key] = key in ['rlar_flag', 'clar_flag', 'tgar_flag']

    # scan inverse aspect ratio (solovev defined between 0 and 0.5)
    for a in np.logspace(-3, 0, 6, endpoint=False)[::-1]:
        qsub = not (a == 1e-3)
        inputs['sol']['SOL_INPUT']['a'] = a
        gpec.run(loc='{l}/lar_lim/a{a:.2e}'.format(l=loc, a=a), rundir=exedir, qsub=qsub,
                 mem=2e3, mailon='', rundcon=True, rungpec=True, runpentrc=True, **inputs)

    return True


def check(loc=loc):
    """
    Plot results of run for a quick check of the limit.

    Args:
        loc: str. Path to directory containing aspect ratio scan
            done by run function.

    Returns:
        bool. True if the large aspect ratio limit torque method agree within 1 percent.

    """
    # setup
    loc = os.path.abspath(loc)
    results = None

    # collect results
    for d in os.listdir(loc):
        filename = '{l}/{d}/pentrc_output_n1.nc'.format(l=loc, d=d)
        if os.path.isfile(filename):
            pentrc = data.open_dataset(filename)
            sol = gpec.namelist.read('{l}/{d}/sol.in'.format(l=loc, d=d))
            if results is None:
                # assign first results
                results = pentrc.assign_coords(a=[sol['SOL_INPUT']['a']])
            else:
                # grow the a dimension
                #results = xarray.concat((results, pentrc.assign_coords(a=[sol['SOL_INPUT']['a']])), dim='a')
                results = xarray.auto_combine((results, pentrc.assign_coords(a=[sol['SOL_INPUT']['a']])), concat_dim='a')
    # make sure order is monotonic
    results = results.isel(a=results['a'].argsort())

    # prep test of success, simple figure and simple table
    success = True
    fig, ax = plt.subplots(2, 1, sharex=True)
    print('{:>24} {:>24} {:>24}'.format('Method', '% Error', 'Success'))

    # test each case
    tgar = results['T_tgar'].sel(psi_tgar=1,method='nearest').sum(dim='ell')
    for key in ['tgar', 'clar', 'rlar']:
        tphi = results['T_' + key].sel(method='nearest', **{'psi_'+key: 1}).sum(dim='ell')
        ax[0].plot(tphi['a'], np.real(tphi), marker='o', mec='none', lw=2, label=key)
        ax[1].plot(tphi['a'], np.real(tphi) / np.real(tgar), marker='o', mec='none', lw=2, label=key + '/tgar')
        perr = 100 * np.abs(1 - tphi / tgar).sel(a=0, method='nearest')
        if perr > 1:
            success = False
        print('{:>24} {:>24.2} {:>24}'.format(key, perr.values, success))
    for a in fig.axes:
        a.legend(numpoints=1)
        a.set_xscale('log')
        a.set_yscale('log') #'('symlog',ylinthresh=1e-24)

    return success


# this is what will happen if you run the script as an exe
if __name__ == "__main__":
    if len(sys.argv[1:]) > 0 and sys.argv[1] == '-c':
        check(*sys.argv[2:])
    elif len(sys.argv[1:]) > 0 and sys.argv[1] == '-h':
        print("Regression test requires two inputs.\n" +
              "Use <loc> <rundir> to run test.\n" +
              "Use -c <loc> to check results.")
    else:
        run(*sys.argv[1:])
