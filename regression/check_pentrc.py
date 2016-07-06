#!/usr/bin/env python
"""

Regression Testing - Testing Standard PENTRC Output
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

def check_ascii_profiles(*directories):
    """

    Check the ntv profile.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    methods=['fgar','tgar']
    nm = len(methods)
    all_datasets = []

    for method in methods:
        print('\n'+method.upper())
        _reset_lines()
        datasets = []
        fig,axes = plt.subplots(4,sharex=True)
        for i,d in enumerate(directories):
            filename = d+'/pentrc_{:}_ell_n1.out'.format(method)
            if os.path.isfile(filename):
                ls = next(linecycler)
                ms = next(markercycler)
                ds = data.open_dataset(filename,quiet=True)
                datasets.append(ds)
                keys = ['T_phi','2ndeltaW','Gamma','chi']
                for ax,key in zip(axes,keys):
                    value = ds[key]
                    lines = []
                    for ell in ds['ell'].values:
                        l, = np.abs(value.sel(ell=ell)).plot(ax=ax)
                        if i==0: l.set_label('ell = {:}'.format(ell))
                        l.set_linestyle(ls)
                        l.set_marker(ms)
                        lines.append(l)
                    sm = plt.set_linearray(lines,ds['ell'].values,cmap='viridis')
        # clean up axes
        for ax in axes.flatten():
            if len(ax.lines)>0:
                ax.set_title('')
                leg1 = ax.legend(loc=2,ncol=2)
        axes[0].set_title(method.upper())
        all_datasets += datasets
    plt.show()

    return all_datasets

def check_ascii_integrals(*directories):
    """

    Check the ell dependence of the integrated torques.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    methods=['fgar','tgar']
    nd = len(directories)
    short_directories = [d[:18]+'...' if len(d)>18 else d for d in directories]
    all_datasets = []

    for method in methods:
        print('\n'+method.upper())
        _reset_lines()
        datasets = []
        fig,axes = plt.subplots(2,sharex=True)
        for i,d in enumerate(directories):
            filename = d+'/pentrc_{:}_ell_n1.out'.format(method)
            if os.path.isfile(filename):
                ls = next(linecycler)
                ms = next(markercycler)
                ds = data.open_dataset(filename,quiet=True)
                datasets.append(ds)
                keys = ['intT_phi','int2ndeltaW']
                for ax,key in zip(axes,keys):
                    value = ds[key]
                    l, = np.abs(value.sel(psi_n=1,method='nearest')).plot(ax=ax)
                    l.set_linestyle(ls)
                    l.set_marker(ms)
                    l.set_label(short_directories[i])
        # clean up axes
        for ax in axes.flat:
            ax.set_title(method.upper()+', '+ax.get_title())
        # print totals
        print('\n{:<24}'.format('')+('{:<24}'*nd).format(*short_directories))
        print('-'*(24*(nd+1)))
        for key in ['intT_phi','int2ndeltaW']:
            values = map(lambda dataset: sum(dataset[key].sel(psi_n=1,method='nearest').values), datasets)
            print('{:<24}'.format(key)+('{:<24.4e}'*nd).format(*values))
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


if __name__ == "__main__":
    check_all(*sys.argv[1:])