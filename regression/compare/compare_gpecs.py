#!/usr/bin/env python
"""

Regression Testing - Testing Standard GPEC Output
======================================================

:Date: 07/02/2016

Run notes:
    -

Result notes:
    -


"""

import numpy as np
import sys, inspect, itertools, os
from matplotlib import pyplot
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-", "--", "-.", ":"]
linecycler = itertools.cycle(linestyles)


############################################ Functions


def _reset_lines(cycler):
    """
    Reset the line cycle.

    :param cycler: cycle. The global linecycler.

    """
    i = 0  # just to be safe and never have an infinite loop
    while cycler.next() != linestyles[-1]:
        i += 1
        if i > 1000: break


def _update_name_conventions(dataset, version=None, inplace=False):
    """
    Update the naming conventions in a dataset to the latest conventions.

    Use this to avoid lots of if statements in cross-checking versions.

    :param dataset:
    :param version: str. Override version attribute.
    :return:

    """
    # set to operate on
    if inplace:
        newset = dataset
    else:
        newset = dataset.copy(deep=True)
    # version number
    if version is None:
        version = dataset.attrs.get('version', None)
        if version is None:
            raise ValueError('Must specify version')
    version = version.split()[-1]  # just the major.minor.patch numbers
    version = np.sum(np.array([1, 0.1, 0.001]) * map(np.float, version.split('.')))  # float representation

    if version < 0.3:
        # profile output changes
        old += ['xi_psi1', 'xi_psi', 'xi_alpha']
        new += ['xigradpsi_dpsi', 'xigradpsi', 'xigradalpha']
    if version < 0.4:
        # control output changes
        old = ['b_xnm', 'xi_xnm', 'b_nm', 'xi_nm', 'b_xn', 'xi_xn', 'b_n', 'xi_n']
        new = ['b_n_x', 'xi_n_x', 'b_n', 'xi_n', 'b_n_x_fun', 'xi_n_x_fun', 'b_n_fun', 'xi_n_fun']
        # profile output changes
        old += ['derxi_m_contrapsi', 'xi_m_contrapsi', 'xi_m_contraalpha']
        new += ['xigradpsi_dpsi', 'xigradpsi', 'xigradalpha']
    else:
        old = []
        new = []

    # swap out any old names that are in the dataset
    name_dict = dict(oldnew for oldnew in zip(old, new) if oldnew[0] in dataset)
    newset = dataset.rename(name_dict, inplace=inplace)

    return newset


def check_control_attributes(*directories):
    """

    Check basic global gpec results.

    Args:
        *args: strings. Any number of directory paths containing gpec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile gpec netcdfs

    """
    nd = len(directories)
    datasets = []

    for i, d in enumerate(directories):
        gpath = d + '/gpec_control_output_n1.nc'
        ipath = d + '/ipec_control_output_n1.nc'
        if os.path.isfile(gpath):
            ds = data.open_dataset(gpath)
        elif os.path.isfile(ipath):
            ds = data.open_dataset(ipath)
        else:
            print('WARNING: No output in ' + d)
            continue
        datasets.append(ds)
        print('{:} from {:}'.format(ds.attrs['version'], d))
    # check that the basic parameters match
    for key in ['machine', 'shot', 'time', 'n', 'jacobian', 'helicity', 'q_lim']:
        values = map(lambda dataset: dataset.attrs[key], datasets)
        if values.count(values[0]) != len(values):
            raise ValueError("The {k} not the same for all cases.".format(k=key))
    # show the global results side by side
    key = 'version'
    values = map(lambda dataset: dataset.attrs[key], datasets)
    print('\n{:<24}'.format(key.upper()) + ('{:<24}' * nd).format(*values))
    print('-' * (24 * (nd + 1)))
    for key in ['energy_vacuum', 'energy_surface', 'energy_plasma']:
        values = map(lambda dataset: dataset.attrs[key], datasets)
        print('{:<24}'.format(key) + ('{:<24.4e}' * nd).format(*values))
    print('')

    return datasets


def check_control_1d(*directories):
    """

    Check control surface xi and b outputs.

    Args:
        *args: strings. Any number of directory paths containing gpec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile gpec netcdfs

    """
    _reset_lines(linecycler)
    axes_dict = {}
    datasets = []

    for dnum, d in enumerate(directories):
        # new linestyle for each directory
        ls = next(linecycler)
        print((ls, d))

        # read and store the control output
        ds = None
        for f in os.listdir(d):
            if 'pec_control_output_n' in f:
                ds = data.open_dataset(d + '/' + f)
                break
        if ds is None:
            print('WARNING: No output in ' + d)
            continue
        ds = _update_name_conventions(ds, inplace=True)
        datasets.append(ds)

        # 1D plots
        for dim in ds.dims.keys():
            allkeys = [k for k, v in ds.data_vars.iteritems() if v.dims == (dim,) and len(v.data) > 1]
            groups_of_4 = [allkeys[i:i + 4] for i in xrange(0, len(allkeys), 4)]
            for keys in groups_of_4:
                if ls == '-':  # first in cycle
                    f, axes = plt.subplots(len(keys), 2, squeeze=False)
                    axes_dict.update(dict(zip(keys, axes)))
                for key in keys:
                    print(key)
                    tmp = ds[key].copy()
                    label = '{:} {:}'.format(key, ds.attrs['version'].split()[-1])
                    ax = axes_dict.get(key, None)
                    if ax is not None:
                        np.abs(tmp).plot(ax=ax[0], linestyle=ls, label=label)
                        tmp.data = np.angle(tmp)
                        tmp.plot(ax=ax[1], linestyle=ls, label=label.replace(key, key + ' phase'))

    # clean up the figure labeling
    for key, axes in axes_dict.iteritems():
        for ax in axes:
            ax.set_title('')
            ax.legend()

    return datasets


def check_control_2d(*directories):
    """

    Check control surface xi and b outputs.

    Args:
        *args: strings. Any number of directory paths containing gpec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile gpec netcdfs.

    """
    _reset_lines(linecycler)
    axes_dict = {}
    matrix_keys = []
    datasets = []

    for dnum, d in enumerate(directories):
        # new linestyle for each directory
        ls = next(linecycler)
        print((ls, d))

        # read and store the control output
        ds = None
        for f in os.listdir(d):
            if 'pec_control_output_n' in f:
                ds = data.open_dataset(d + '/' + f)
                break
        if ds is None:
            print('WARNING: No output in ' + d)
            continue
        ds = _update_name_conventions(ds, inplace=True)
        datasets.append(ds)

        # 2D plots
        allkeys = [k for k, v in ds.data_vars.iteritems() if len(v.dims) == 2 and min(np.shape(v)) > 1]
        groups_of_4 = [allkeys[i:i + 4] for i in xrange(0, len(allkeys), 4)]
        for keys in groups_of_4:
            if ls == '-':  # first in cycle
                matrix_keys = allkeys
                f, axes = plt.subplots(len(keys), len(directories) + 1, squeeze=False)
                axes_dict.update(dict(zip(keys, axes)))
            for key in keys:
                print(key)
                tmp = ds[key].copy()
                if 'theta' in tmp.dims:
                    tmp = tmp.sel(theta=np.linspace(0, 1, 60), method='nearest')
                ax = axes_dict.get(key, None)
                if ax is not None:
                    np.abs(tmp).plot.imshow(ax=ax[dnum])

    # standard deviation of 2D plots
    for key in matrix_keys:
        std = np.apply_over_axes(np.std, np.array([np.abs(tmp[key]) for tmp in datasets]), 0)[0]
        std = data.xarray.DataArray(std, coords=tmp[key].coords, name='std')
        if 'theta' in std.dims:
            std = std.sel(theta=np.linspace(0, 1, 60), method='nearest')
        std.plot.imshow(ax=axes_dict[key][-1])

    # clean up the figure labeling
    for key, axes in axes_dict.iteritems():
        for ax in axes:
            ax.set_title('')

    return datasets


def check_profile_output(*directories):
    """

    :param directories: strings. Any number of directory paths containing gpec outputs.
    :return: List of Datasets from the corresponding profile gpec netcdfs.

    """
    _reset_lines(linecycler)
    axes_dict = {}
    datasets = []

    for dnum, d in enumerate(directories):
        # new linestyle for each directory
        ls = next(linecycler)

        # read and store the control output
        ds = None
        for f in os.listdir(d):
            if 'pec_profile_output_n' in f:
                ds = data.open_dataset(d + '/' + f)
                break
        if ds is None:
            print('WARNING: No output in ' + d)
            continue
        ds = _update_name_conventions(ds, inplace=True)
        datasets.append(ds)

        # 1D plots
        for dim, dsel in [('m_out', np.arange(5)), ('theta_dcon', np.linspace(0, 1, 4, False))]:
            allkeys = [k for k, v in ds.data_vars.iteritems() if v.dims == (dim, 'psi_n')]
            groups_of_4 = [allkeys[i:i + 4] for i in xrange(0, len(allkeys), 4)]
            for keys in groups_of_4:
                if ls == '-':  # first in cycle
                    f, axes = plt.subplots(len(keys), 2, squeeze=False)
                    axes_dict.update(dict(zip(keys, axes)))
                for key in keys:
                    ax = axes_dict.get(key, None)
                    if ax is not None:
                        for sel in dsel:
                            label = '{:} = {:}'.format(dim, sel) * (dnum == 0)
                            ds[key].sel(**{dim: sel}).real.plot(ax=ax[0], linestyle=ls, label=label)
                            ds[key].sel(**{dim: sel}).imag.plot(ax=ax[1], linestyle=ls,
                                                                label=label.replace('real ', 'imag '))
                        sm = plt.set_linearray(ax[0].lines[-len(dsel):], dsel, cmap='viridis')
                        sm = plt.set_linearray(ax[1].lines[-len(dsel):], dsel, cmap='viridis')

    # clean up the figure labeling
    versions = map(lambda dataset: dataset.attrs['version'], datasets)
    for key, axes in axes_dict.iteritems():
        axes[0].set_title('Real')
        axes[1].set_title('Imaginary')
        for ax in axes:
            if len(ax.lines) > 0:
                nsel = len(ax.lines) / len(versions)
                leg1 = ax.legend(loc=2)
                leg2 = ax.legend(ax.lines[::nsel], versions, loc=1)
                ax.add_artist(leg1)
    plt.show()

    return datasets


def check_cylindircal_output(*directories):
    """

    Plot the reference and std of cylindircal output

    :param directories: strings. Any number of directory paths containing gpec outputs.

    :returns datasets: List of datasets from the corresponding profile gpec netcdfs.

    """
    _reset_lines(linecycler)
    axes_dict = {}
    datasets = []

    for dnum, d in enumerate(directories):
        # new linestyle for each directory
        ls = next(linecycler)
        print((ls, d))

        # read and store the control output
        ds = None
        for f in os.listdir(d):
            if 'pec_cylindrical_output_n' in f:
                ds = data.open_dataset(d + '/' + f)
                break
        if ds is None:
            print('WARNING: No output in ' + d)
            continue
        ds = _update_name_conventions(ds, inplace=True)
        datasets.append(ds)

        # 2D plots
        allkeys = [k for k, v in ds.data_vars.iteritems() if v.dims == ('z', 'R')]
        groups_of_4 = [allkeys[i:i + 4] for i in xrange(0, len(allkeys), 4)]
        for keys in groups_of_4:
            if dnum == 0:  # first in cycle
                f, axes = plt.subplots(len(keys), 2, squeeze=False)
                axes_dict.update(dict(zip(keys, axes)))
            for key in keys:
                ax = axes_dict.get(key, None)
                if axes_dict is not None:
                    for ndim, dim in enumerate(['z', 'R']):
                        sel = np.float(ds[dim].mean())
                        label = '{:} = {:.3}'.format(dim, sel) * (dnum == 0)
                        c = None
                        if dnum > 0: c = ax[0].lines[ndim].get_color()
                        ds[key].sel(method='nearest', **{dim: sel}).real.plot(ax=ax[0], label=label, ls=ls, color=c)
                        ds[key].sel(method='nearest', **{dim: sel}).imag.plot(ax=ax[1], label=label, ls=ls, color=c)

    # clean up the figure labeling
    versions = map(lambda dataset: dataset.attrs['version'], datasets)
    for key, axes in axes_dict.iteritems():
        axes[0].set_title('Real')
        axes[1].set_title('Imaginary')
        for ax in axes:
            if len(ax.lines) > 0:
                nsel = len(ax.lines) / len(versions)
                leg1 = ax.legend(loc=2)
                leg2 = ax.legend(ax.lines[::nsel], versions, loc=1)
                ax.add_artist(leg1)
                ax.set_xlabel('z, R')
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
