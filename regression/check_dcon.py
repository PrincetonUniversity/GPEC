#!/usr/bin/env python
"""

Regression Testing - Testing Standard DCON Output
======================================================

:Date: 07/02/2016

Run notes:
    -

Result notes:
    -


"""

import numpy as np
import os, sys, inspect, itertools
from matplotlib import pyplot
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-","--","-.",":"]
linecycler = itertools.cycle(linestyles)
markerstyles = ["x","+","|","."]
markercycler = itertools.cycle(markerstyles)

def _reset_lines(linecycler=linecycler,markercycler=markercycler):
    '''Reset the line cycle'''
    i=0 # just to be safe and never have an infinite loop
    while linecycler.next()!=linestyles[-1]:
        i+=1
        if i>1000: break
    i=0 # just to be safe and never have an infinite loop
    while markercycler.next()!=markerstyles[-1]:
        i+=1
        if i>1000: break

############################################ Functions

def check_control_matrices(*directories):
    """

    Check the W matrix

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    datasets = []

    return datasets

def check_energies(*directories):
    """

    Check the eigenvalues

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    datasets = []
    nd = len(directories)
    short_directories = [d[:18]+'...' if len(d)>18 else d for d in directories]
    fig,axes = plt.subplots(4,sharex=True)
    _reset_lines()

    for i,d in enumerate(directories):
        ls = next(linecycler)
        ms = next(markercycler)
        lbl = short_directories[i]
        ds = data.read(d+'/dcon.out',quiet=True)
        datasets.append(ds)
        axes[0].plot(ds[2].x[0],np.real(ds[2].y['total']),label=lbl,ls=ls,marker=ms)
        axes[1].plot(ds[2].x[0],np.imag(ds[2].y['total']),label=lbl,ls=ls,marker=ms)
        axes[2].plot(ds[2].x[0],ds[2].y['plasma'],label=lbl,ls=ls,marker=ms)
        axes[3].plot(ds[2].x[0],ds[2].y['vacuum'],label=lbl,ls=ls,marker=ms)

    # label axes
    axes[0].set_ylabel('Re Total Eigenvalue')
    axes[1].set_ylabel('Im Total Eigenvalue')
    axes[2].set_ylabel('Plasma Eigenvalue')
    axes[3].set_ylabel('Vacuum Eigenvalue')
    axes[3].set_xlabel('Eigenmode Index')
    for ax in axes:
        ax.set_yscale('symlog',linthreshy=1e-3)

    # print totals
    print('\n{:<24}'.format('')+('{:<24}'*nd).format(*short_directories))
    print('-'*(24*(nd+1)))
    for key in ['total','plasma','vacuum']:
        values = map(lambda dataset: dataset[2].y[key][0], datasets)
        print('{:<24}'.format(key)+('{:<+24.4e}'*nd).format(*values))
    return datasets

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


if __name__ == "__main__":
    check_all(*sys.argv[1:])