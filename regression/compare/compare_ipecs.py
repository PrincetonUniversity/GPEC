#!/usr/bin/env python
"""

Regression Testing - Testing Standard IPEC Output
======================================================

:Date: 07/02/2016

Run notes:
    -

Result notes:
    -


"""

import numpy as np
import sys, inspect, itertools
from matplotlib import pyplot
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-","--","-.",":"]
linecycler = itertools.cycle(linestyles)

def _reset_lines(linecycler):
    '''Reset the line cycle'''
    i=0 # just to be safe and never have an infinite loop
    while linecycler.next()!=linestyles[-1]:
        i+=1
        if i>1000: break

############################################ Functions

def check_control_attributes(*directories):
    """

    Check basic global ipec results.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    nd = len(directories)
    datasets = []

    for i,d in enumerate(directories):
        ds = data.open_dataset(d+'/ipec_control_output_n1.nc')
        datasets.append(ds)
        print('{:} from {:}'.format(ds.attrs['version'],d))
    # check that the basic parameters match
    for key in ['machine','shot','time','n','jacobian','helicity','q_lim']:
        values = map(lambda dataset: dataset.attrs[key], datasets)
        if values.count(values[0]) != len(values):
            raise ValueError("The {k} not the same for all cases.".format(k=key))
    # show the global results side by side
    key = 'version'
    values = map(lambda dataset: dataset.attrs[key], datasets)
    print('\n{:<24}'.format(key.upper())+('{:<24}'*nd).format(*values))
    print('-'*(24*(nd+1)))
    for key in ['energy_vacuum','energy_surface','energy_plasma']:
        values = map(lambda dataset: dataset.attrs[key], datasets)
        print('{:<24}'.format(key)+('{:<24.4e}'*nd).format(*values))
    print('')

    return datasets

def check_control_perturbations(*directories):
    """

    Check control surface xi and b outputs.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    _reset_lines(linecycler)
    fig,axes = plt.subplots(2,2)
    #ffun,afuns = plt.subplots(2,2)
    datasets = []

    for i,d in enumerate(directories):
        ls = next(linecycler)
        ds = data.open_dataset(d+'/ipec_control_output_n1.nc')
        datasets.append(ds)

        # spectral plots
        keys = ['b_xnm','xi_xnm','b_nm','xi_nm']
        altkeys = ['b_n_x','xi_n_x','b_n','xi_n']
        maxes = list(axes[:,0])*2
        for ax,key,altkey in zip(maxes,keys,altkeys):
            value = ds.get(key,ds.get(altkey,None))
            if value is None:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_control_output_n1.nc'))
            else:
                l, = np.abs(value).plot(ax=ax)
                l.set_linestyle(ls)
                l.set_label('{:} {:}'.format(value.name,ds.attrs['version']))
        # function plots
        keys = ['b_n_x_fun','xi_n_x_fun','b_n_fun','xi_n_fun']
        altkeys = ['b_xn','xi_xn','b_n','xi_n']
        faxes = list(axes[:,1])*2
        for ax,key,altkey in zip(faxes,keys,altkeys):
            value = ds.get(key,ds.get(altkey,None))
            if value is None or value.dims[0]!='theta':
                print("KeyError: {:} function not in {:}".format(key,d+'/ipec_control_output_n1.nc'))
            else:
                l, = np.abs(value).plot(ax=ax)
                l.set_linestyle(ls)
                l.set_label('{:} {:}'.format(value.name,ds.attrs['version']))

    # clean up the figure labeling
    for ax in fig.axes:
        if len(ax.lines)>0:
            ax.set_title('')
            ax.legend()

    return datasets

def check_control_matrices(*directories):
    """

    Check control surface L, Lambda, P, and rho matrices.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    _reset_lines(linecycler)
    nd = len(directories)
    datasets = []

    for d1,d2 in zip(directories,np.roll(directories,1)):
        # no need to repeat comparison
        if nd==2 and len(datasets)>0: continue
        # compare two directories per figure
        fig,axes = plt.subplots(3,4)
        ds1 = data.open_dataset(d1+'/ipec_control_output_n1.nc')
        ds2 = data.open_dataset(d2+'/ipec_control_output_n1.nc')
        datasets.append(ds1)
        datasets.append(ds2)

        keys = ['L','Lambda','P','rho']
        for i,key in enumerate(keys):
            if key not in ds1:
                print("KeyError: {:} not in {:}".format(key,d1+'/ipec_control_output_n1.nc'))
            elif key not in ds2:
                print("KeyError: {:} not in {:}".format(key,d2+'/ipec_control_output_n1.nc'))
            else:
                ds1[key].plot(ax=axes[0,i])
                axes[0,i].set_title(ds2.attrs.get('version','Unknown'))
                ds2[key].plot(ax=axes[1,i])
                axes[1,i].set_title(ds2.attrs.get('version','Unknown'))
                (ds1[key]-ds2[key]).plot(ax=axes[2,i])
                axes[2,i].set_title('Difference')
        for ax in axes.ravel():
            ax.set_aspect('equal')

    return datasets

def check_xbnormal(*directories):
    """

    Check xbnormal netcdf outputs.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    _reset_lines(linecycler)
    fig,axes = plt.subplots(2,sharex=True)
    ffun,afuns = plt.subplots(2,sharex=True)
    datasets = []

    for i,d in enumerate(directories):
        ls = next(linecycler)
        ds = data.open_dataset(d+'/ipec_profile_output_n1.nc')
        datasets.append(ds)

        # 1D plots
        ms = range(1,6)
        keys = ['b_n','xi_n']
        for ax,key in zip(axes,keys):
            if key not in ds:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_profile_output_n1.nc'))
            else:
                lines = []
                for m in ms:
                    l, = np.abs(ds[key]).sel(m=m).plot(ax=ax)
                    l.set_linestyle(ls)
                    # label the first set of ms
                    if i==0:
                        l.set_label('m = {:}'.format(m))
                    lines.append(l)
                # color the lines based on m
                sm = plt.set_linearray(lines,ms,cmap='viridis')

        # plot 1D slices of functions
        thetas = np.linspace(0,1,4,False)
        keys = ['b_n_fun','xi_n_fun']
        for ax,key in zip(afuns,keys):
            if key not in ds:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_profile_output_n1.nc'))
            else:
                lines = []
                for theta in thetas:
                    l, = np.abs(ds[key]).sel(theta=theta).plot(ax=ax)
                    l.set_linestyle(ls)
                    # label the first set of thetas
                    if i==0:
                        l.set_label('theta = {:.2f}'.format(theta))
                    lines.append(l)
                # color the lines based on theta
                sm = plt.set_linearray(lines,thetas,cmap='viridis')

    # clean up the figure labeling
    versions = map(lambda dataset: dataset.attrs['version'],datasets)
    for ax in fig.axes+ffun.axes:
        if len(ax.lines)>0:
            ax.set_title('')
            leg1 = ax.legend(loc=2)
            leg2 = ax.legend(ax.lines[::len(ms)],versions,loc=1)
            ax.add_artist(leg1)
    plt.show()

    return datasets

def check_pmodb(*directories):
    """

    Check xbnormal netcdf outputs.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    _reset_lines(linecycler)
    fig,axes = plt.subplots(2,sharex=True)
    ffun,afuns = plt.subplots(2,sharex=True)
    datasets = []

    for i,d in enumerate(directories):
        ls = next(linecycler)
        ds = data.open_dataset(d+'/ipec_profile_output_n1.nc')
        datasets.append(ds)

        # 1D plots
        ms = range(1,6)
        keys = ['b_lag','Bdivxi_perp']
        for ax,key in zip(axes,keys):
            if key not in ds:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_profile_output_n1.nc'))
            else:
                lines = []
                for m in ms:
                    l, = np.abs(ds[key]).sel(m=m).plot(ax=ax)
                    l.set_linestyle(ls)
                    # label the first set of ms
                    if i==0:
                        l.set_label('m = {:}'.format(m))
                    lines.append(l)
                # color the lines based on m
                sm = plt.set_linearray(lines,ms,cmap='viridis')

        # plot 1D slices of functions
        thetas = np.linspace(0,1,4,False)
        keys = ['b_lag_fun','Bdivxi_perp_fun']
        for ax,key in zip(afuns,keys):
            if key not in ds:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_profile_output_n1.nc'))
            else:
                lines = []
                for theta in thetas:
                    l, = np.abs(ds[key]).sel(theta=theta).plot(ax=ax)
                    l.set_linestyle(ls)
                    # label the first set of thetas
                    if i==0:
                        l.set_label('theta = {:.2f}'.format(theta))
                    lines.append(l)
                # color the lines based on theta
                sm = plt.set_linearray(lines,thetas,cmap='viridis')

    # clean up the figure labeling
    versions = map(lambda dataset: dataset.attrs['version'],datasets)
    for f in [fig,ffun]:
        nlines = 0
        for ax in f.axes:
            nlines+=len(ax.lines)
            if len(ax.lines)>0:
                ax.set_title('')
                leg1 = ax.legend(loc=2)
                leg2 = ax.legend(ax.lines[::len(ms)],versions,loc=1)
                ax.add_artist(leg1)
        if nlines==0:
            pyplot.close(f)
    plt.show()

    return datasets

def check_xclebsch(*directories):
    """

    Check xclebsch netcdf outputs.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    _reset_lines(linecycler)
    fig,axes = plt.subplots(3,sharex=True)
    datasets = []

    for i,d in enumerate(directories):
        ls = next(linecycler)
        ds = data.open_dataset(d+'/ipec_profile_output_n1.nc')
        datasets.append(ds)

        # 1D plots
        ms = range(1,6)
        keys = ['xi_psi1','xi_psi','xi_alpha']
        altkeys = ['derxi_m_contrapsi','xi_m_contrapsi','xi_m_contraalpha']
        for ax,key,altkey in zip(axes,keys,altkeys):
            value = ds.get(key,ds.get(altkey,None))
            if value is None:
                print("KeyError: {:} not in {:}".format(key,d+'/ipec_profile_output_n1.nc'))
            else:
                lines = []
                for m in ms:
                    l, = np.abs(value).sel(m=m).plot(ax=ax)
                    l.set_linestyle(ls)
                    # label the first set of m's
                    if i==0:
                        l.set_label('m = {:}'.format(m))
                    lines.append(l)
                # color the lines based on m
                sm = plt.set_linearray(lines,ms,cmap='viridis')

    # clean up the figure labeling
    versions = map(lambda dataset: dataset.attrs['version'],datasets)
    for ax in fig.axes:
        if len(ax.lines)>0:
            ax.set_title('')
            leg1 = ax.legend(loc=2)
            leg2 = ax.legend(ax.lines[::len(ms)],versions,loc=1)
            ax.add_artist(leg1)
    plt.show()

    return datasets

############################################ Run


def check_all(*directories):
    # generic way to collect all the functions
    this_module = sys.modules[__name__]
    all_functions = inspect.getmembers(this_module, inspect.isfunction)
    # call each function
    for key, function in all_functions:
        if key[0] != '_' and key != 'check_all':
            print("Calling {:}".format(key))
            datasets = function(*directories)


if __name__ == "__main__":
    check_all(*sys.argv[1:])