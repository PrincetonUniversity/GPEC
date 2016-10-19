#!/usr/bin/env python
"""

Regression Test: Test Self Consistency of GPEC Torque Calculations
==================================================================



"""

import os
import numpy as np
import gpec, data

plt = data.plt

repo = '/'.join(__file__.split('/')[:-3])
loc = repo + '/docs/examples/DIIID_kinetic_example'


def check(loc=loc, frac_error=0.05):
    """

    :param loc: str. Path to a kinetic run of GPEC.
    :return: bool. All methods agree within frac_error.

    """
    # find out toroidal mode number
    files = [f for f in os.listdir(loc) if '.nc' in f]
    nn = files[0].split('.')[0].split('n')[-1]
    # read results
    con = data.open_dataset(loc + '/gpec_control_output_n{:}.nc'.format(nn))
    prof = data.open_dataset(loc + '/gpec_profile_output_n{:}.nc'.format(nn))
    pent = data.open_dataset(loc + '/pentrc_output_n{:}.nc'.format(nn))

    # get the three methods
    torques = data.xarray.Dataset()
    torques['PENTRC'] = pent['T_fgar'].sum('ell')
    torques['GPEC_profile'] = prof['T']
    M = prof['T_xe'].transpose('psi_n', 'm', 'm_prime')
    T = np.dot(np.dot(con['Phi_xe'].T.conj(), M), con['Phi_xe']) / 2 # 1/2 for the quadratic use of linear exp
    torques['GPEC_matrix'] = data.xarray.DataArray(T, coords={'psi_n': prof['psi_n']})

    def xreal(darray):
        return darray * 0 + np.real(darray)

    def ximag(darray):
        return darray * 0 + np.imag(darray)

    def isclose(val1, val2):
        rerr = np.abs(np.real(val1 - val2) / np.real(val1))
        ierr = np.abs(np.imag(val1 - val2) / np.imag(val1))
        return (rerr < frac_error) & (ierr < frac_error)

    # summary
    tt = torques['PENTRC'].values[-1]
    test = True
    f, ax = plt.subplots(1, 2)
    print('{:>24} {:>24} {:>24}\n'.format('Method', 'Total', 'Fractional-Error') + '-' * 25 * 3)
    for key, val in torques.data_vars.iteritems():
        tk = val.values[-1]
        test = test & isclose(tt, tk)
        rerr = np.abs(np.real(tt - tk) / np.real(tt))
        ierr = np.abs(np.imag(tt - tk) / np.imag(tt))
        print'{:>24} {:24.3e} {:24.3f}'.format(key, tk, rerr+1j*ierr)
        xreal(val).plot(ax=ax[0], label=key)
        ximag(val).plot(ax=ax[1], label=key)

    ax[0].set_ylabel('REAL')
    ax[1].set_ylabel('IMAG')
    for a in ax:
        a.legend()

    return test
