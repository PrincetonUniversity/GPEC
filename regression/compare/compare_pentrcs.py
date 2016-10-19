#!/usr/bin/env python
"""

:mod:`compare.compare_pentrcs` -- Regression Testing Standard PENTRC Output
===========================================================================

This module contains a library of functions for visualizing standard PENTRC
output in a way that is conducive to comparing multiple PENTRC runs.

"""

import numpy as np
import os, sys, inspect, itertools
from matplotlib import pyplot
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-", "--", "-.", ":"]
linecycler = itertools.cycle(linestyles)
markerstyles = ["x", "+", "|", "."]
markercycler = itertools.cycle(markerstyles)


def _reset_lines(linecycler=linecycler, markercycler=markercycler):
    '''Reset the line cycle'''
    i = 0  # just to be safe and never have an infinite loop
    while linecycler.next() != linestyles[-1]:
        i += 1
        if i > 1000: break
    i = 0  # just to be safe and never have an infinite loop
    while markercycler.next() != markerstyles[-1]:
        i += 1
        if i > 1000: break


############################################ Functions

def check_profiles(*directories):
    """

    Check the ntv profile.

    Args:
        *args: strings. Any number of directory paths containing gpec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile gpec netcdfs

    """
    methods = ['fgar']#, 'tgar']
    nm = 0
    all_datasets = []

    for method in methods:
        print('\n' + method.upper())
        _reset_lines()
        datasets = []
        fig, axes = plt.subplots(3, 2, sharex=True)
        for i, d in enumerate(directories):
            # default is to use netcdf
            ascii = d + '/pentrc_{:}_ell_n1.out'.format(method)
            netcdf = d + '/pentrc_output_n1.nc'
            if os.path.isfile(netcdf):
                ds = data.open_dataset(netcdf)
            elif os.path.isfile(ascii):
                ds = data.open_dataset(ascii, quiet=True, auto_complex=False).rename({'psi_n': 'psi_' + method})
                ds['T_{:}'.format(method)] = ds['T_phi'] + 1j * ds['int2ndeltaW']
                ds = ds.rename({'Gamma': 'Gamma_' + method, 'chi': 'chi_' + method})
            else:
                continue
            nm +=1
            ls = next(linecycler)
            ms = next(markercycler)
            datasets.append(ds)
            values = []
            print(ls + ' '+ ms + ' ' + d)
            for key in ['T_' + method, 'Gamma_' + method, 'chi_' + method]:
                values.append(ds[key].real)
                values.append(ds[key].imag)
            for ax, value in zip(axes.flat, values):
                lines = []
                for ell in ds['ell'].values:
                    if np.abs(ell)>1: continue
                    l, = value.sel(ell=ell).plot(ax=ax)
                    if i == 0:
                        l.set_label('ell = {:}'.format(ell))
                    l.set_linestyle(ls)
                    l.set_marker(ms)
                    lines.append(l)
                sm = plt.set_linearray(lines, cmap='viridis')
        # clean up axes
        for ax in axes.flatten():
            if len(ax.lines) > 0:
                ax.set_title('REAL' * (ax is axes[0, 0]) + 'IMAG' * (ax is axes[0, 1]))
                leg1 = ax.legend(loc=2, ncol=1)
                ax.set_yscale('symlog', linthreshy=1e-4)
        all_datasets += datasets
    plt.show()

    return all_datasets


def check_integrals(*directories):
    """

    Check the ell dependence of the integrated torques.

    Args:
        *args: strings. Any number of directory paths containing gpec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile gpec netcdfs

    """
    methods = ['fgar']#, 'tgar']
    short_directories = [d[:18] + '...' if len(d) > 18 else d for d in directories]
    all_datasets = []

    fig, ax2d = plt.subplots(len(methods), 2, sharex=True, squeeze=False)
    for m,method in enumerate(methods):
        _reset_lines()
        datasets = []
        #fig, axes = plt.subplots(2, sharex=True)
        axes = ax2d[m]
        nd = 0
        short_directories = []
        for i, d in enumerate(directories):
            # default is to use netcdf
            ascii = d + '/pentrc_{:}_ell_n1.out'.format(method)
            netcdf = d + '/pentrc_output_n1.nc'
            if os.path.isfile(netcdf):
                ds = data.open_dataset(netcdf)
            elif os.path.isfile(ascii):
                ds = data.open_dataset(ascii, quiet=True).rename({'psi_n': 'psi_' + method})
                ds['T_{:}'.format(method)] = ds['T_phi'] + 1j * ds['int2ndeltaW']
            else:
                continue
            nd += 1
            if len(d)>18:
                shortd = d[:18] + '...'
            else:
                shortd = d
            short_directories.append(shortd)
            ls = next(linecycler)
            ms = next(markercycler)
            datasets.append(ds)
            key = 'T_{:}'.format(method)
            values = [ds[key].real, ds[key].imag]
            for ax, value in zip(axes, values):
                l, = np.abs(value.sel(method='nearest', **{'psi_' + method: 1})).plot(ax=ax)
                l.set_linestyle(ls)
                l.set_marker(ms)
                l.set_label(short_directories[i])
        # clean up axes
        for ax in axes.flat:
            ri = 'REAL' * (ax in ax2d[:,0]) + 'IMAG' * (ax in ax2d[:,1])
            ax.set_title(ri + ' ' + method.upper() + ', psi_n = 1')
            ax.legend()
            ax.set_yscale('log')
        # print totals
        print('\n{:<24}'.format('Value') + ('{:<24}' * nd).format(*short_directories))
        print('-' * (24 * (nd + 1)))
        for key in ['T_' + method]:
            values = map(lambda dataset: sum(dataset[key].sel(method='nearest', **{'psi_' + method: 1}).values),
                         datasets)
            print('{:<24}'.format(key) + ('{:<24.4e}' * nd).format(*values))
        all_datasets += datasets
    plt.show()

    return all_datasets


############################################ Run


def check_all(*directories):
    # generic way to collect all the functions
    this_module = sys.modules[__name__]
    all_functions = inspect.getmembers(this_module, inspect.isfunction)
    # call each function
    for key, function in all_functions:
        if key[0] != '_' and key != 'check_all':
            print("\nCalling {:}".format(key))
            datasets = function(*directories)
    plt.show()


if __name__ == "__main__":
    check_all(*sys.argv[1:])
