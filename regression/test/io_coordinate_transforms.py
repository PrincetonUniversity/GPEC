#!/usr/bin/env python
"""

Regression Test - Test Coordinate Transformed IO between GPEC and PENTRC
=================================================================================

:Date: 08/02/2016

Run notes:
    -

Result notes:
    -


"""

import os
import numpy as np
import gpec, data

plt = data.plt


repo = '/'.join(__file__.split('/')[:-3])
loc = repo + '/regression/runs'


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
    rundir = repo + '/bin'
    exdir = '{:}/docs/examples/DIIID_ideal_example/'.format(repo)
    if not os.path.isdir(exdir):
        exdir = exdir.replace('_ideal_example','_example')

    # gather example input files
    inputs = {}
    for f in os.listdir(exdir):
        if f.endswith('.in'):
            inputs[f[:-3]] = gpec.namelist.read(exdir + f)
    # fill out explicit paths
    eqf = inputs['equil']['EQUIL_CONTROL']['eq_filename']
    if not eqf.startswith('/'):
        inputs['equil']['EQUIL_CONTROL']['eq_filename'] = os.path.abspath(exdir+eqf)
    pf = inputs['pentrc']['PENT_INPUT']['kinetic_file']
    if not pf.startswith('/'):
        inputs['pentrc']['PENT_INPUT']['kinetic_file'] = os.path.abspath(exdir+pf)
    inputs['pentrc']['PENT_INPUT']['data_dir'] = repo + '/pentrc'
    inputs['coil']['COIL_CONTROL']['data_dir'] = repo + '/coil'

    # output transformed coordinate perturbation profiles
    inputs['gpec']['GPEC_OUTPUT']['pmodb_flag'] = True
    # output diagnostic working coordinate perturbation profiles from gpec
    inputs['gpec']['GPEC_DIAGNOSE']['pmodbmn_flag'] = True
    # output diagnostic working coordinate perturbation profiles from pentrc
    if 'indebug' in inputs['pentrc']['PENT_ADMIN']:
        inputs['pentrc']['PENT_ADMIN']['indebug'] = True
    else:
        inputs['pentrc']['PENT_ADMIN']['tdebug'] = True

    # run mix of cases
    for coord in ['pest', 'hamada']:
        for t in [0, 1]:
            inputs['gpec']['GPEC_OUTPUT']['jac_out'] = coord
            inputs['gpec']['GPEC_OUTPUT']['tmag_out'] = t
            inputs['pentrc']['PENT_INPUT']['jac_in'] = coord
            inputs['pentrc']['PENT_INPUT']['tmag_in'] = t
            gpec.run(loc='{l}/{c}{t}'.format(l=loc, c=coord, t=t), rundir=rundir, qsub=True, mem=2e3, mailon='',
                     rundcon=True, rungpec=True, runpentrc=True, **inputs)

    return True

def check(loc=loc,m=2):
    """
    Plot results of run for a quick check of the perturbations and torque result.

    Args:
        loc:

    Returns:

    """
    loc = loc.rstrip('/')
    files = ['ipdiag_pmodb_n1.out', 'ipdiag_pmodbmn_n1.out', 'pentrc_peq_n1.out']
    results = {}

    fx, ax = plt.subplots(2, 2, sharey=True, sharex=True)
    fb, ab = plt.subplots(2, 2, sharey=True, sharex=True)
    fc, ac = plt.subplots(3, 1, sharey=True, sharex=True)
    ft, at = plt.subplots(2, 2, subplot_kw={'yscale': 'symlog'})

    for i, coord in enumerate(['hamada', 'pest']):
        for j, t in enumerate([1, 0]):
            dirname = '{c}{t}'.format(c=coord, t=t)
            dirpath = loc + '/' + dirname
            results[dirname] = {}
            print(dirpath)

            # check the displacements
            f = 'pentrc_peq_n1.out'
            ds = data.open_dataset(dirpath + '/' + f, quiet=True)
            for k,a in [('xi^psi',ac[0]),('xi^psi1',ac[1]),('xi^alpha',ac[2])]:
                l, = np.abs(ds[k].sel(m=m, method='nearest')).plot(ax=a)
                l.set_label(dirname)

            # check the perturbation profiles
            for f in files:
                filename = f.split('.')[0]
                if os.path.isfile(dirpath + '/' + f):
                    ds = data.open_dataset(dirpath + '/' + f, quiet=True)
                    xkey = ([''] + [k for k in ds.data_vars if k.startswith('JBBdivx')])[-1]
                    bkey = (['JBdeltaB'] + [k for k in ds.data_vars if k.startswith('JBdeltaB')])[-1]
                    if bkey not in ds and 'JBBdivx' in ds and 'JBBkx' in ds:
                        ds[bkey] = -(ds['JBBdivx']+ds['JBBkx'])
                    for k,a in [(xkey,ax),(bkey,ab)]:
                        if k in ds:
                            l, = np.abs(ds[k]).sel(m=m).plot(ax=a[i, j])
                            l.set_label(filename)
                        a[i, j].set_title(a[i, j].get_title()+', '+dirname)
                    results[dirname][filename] = ds.copy()
                else:
                    print('** No file ' + dirpath + '/' + f)

            # check final results
            f = 'pentrc_fgar_ell_n1.out'
            ds = data.open_dataset(dirpath + '/' + f, quiet=True)
            l, = ds['intT_phi'].sel(psi_n=1, method='nearest').plot(ax=at[0,0])
            l, = ds['intT_phi'].sum(dim='ell').plot(ax=at[1,0])
            l, = ds['int2ndeltaW'].sel(psi_n=1, method='nearest').plot(ax=at[0,1])
            l, = ds['int2ndeltaW'].sum(dim='ell').plot(ax=at[1,1])
            for a in ft.axes:
                a.lines[-1].set_label(dirname)
            results[dirname]['pentrc_fgar_ell_n1'] = ds.copy()

    # label lines
    for a in fx.axes + fb.axes + fc.axes + ft.axes:
        a.legend()

    # print table
    totals = []
    print('{:24} {:24} {:24}'.format('directory','Re(dW_k)','Im(dW_k)'))
    for k,v in results.iteritems():
        dwk = 0.5*(v['pentrc_fgar_ell_n1']['int2ndeltaW'] + 1j*v['pentrc_fgar_ell_n1']['intT_phi'])
        totals.append(dwk.sel(psi_n=1, method='nearest').sum(dim='ell').values[0])
        print('{:24} {:<+24.4e} {<+24.4e}'.format(k,np.real(totals[-1]),np.imag(totals[-1])))
    totals = np.array(totals)
    werrs = 100*np.abs(1 - np.real(totals)/np.real(totals[0]))
    terrs = 100*np.abs(1 - np.imag(totals)/np.imag(totals[0]))
    success = np.all(werrs<1) and np.all(terrs<1)

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
