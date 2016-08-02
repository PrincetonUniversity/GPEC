#!/usr/bin/env python
"""

Regression Running - Test Coordinate Transformed IO between IPEC and PENTRC
=================================================================================

:Date: 08/02/2016

Run notes:
    -

Result notes:
    -


"""

import os
import gpec


def coordinate_tranformed_io(loc, repo):
    """
    Run as a regression test for any repository under development.

    Args:
        loc: str. Path where runs will be made.
        repo: str. Path to GPEC repository.

    Returns:

    """
    loc = os.path.abspath(loc)
    repo = os.path.abspath(repo)
    rundir = repo + '/bin'
    exdir = '{:}/docs/examples/DIIID_ideal_example/'.format(repo)

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

    # output diagnostic perturbation profiles
    inputs['ipec']['IPEC_DIAGNOSE']['pmodbmn_flag'] = True
    inputs['pentrc']['PENT_ADMIN']['indebug'] = True

    for coord in ['pest', 'hamada']:
        for t in [0, 1]:
            inputs['ipec']['IPEC_OUTPUT']['jac_out'] = coord
            inputs['ipec']['IPEC_OUTPUT']['tmag_out'] = t
            inputs['pentrc']['PENT_INPUT']['jac_in'] = coord
            inputs['pentrc']['PENT_INPUT']['tmag_in'] = t
            gpec.run(loc='{l}/{c}{t}'.format(l=loc, c=coord, t=t), rundir=rundir, qsub=True, mem=2e3, mailon='',
                     rundcon=True, runipec=True, runpentrc=True, **inputs)
