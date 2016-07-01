#!/usr/bin/env python
"""

Regression Testing - Testing IPEC's xbnormal profiles
======================================================

:Date: 07/02/2016

Run notes:
    -

Result notes:
    -


"""

import numpy as np
import sys, itertools
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-","--","-.",":"]
linecycler = itertools.cycle(linestyles)

############################################ Functions

def plot_xbnormal(*dirs):
    """

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    fig,axes = plt.subplots(2)
    datasets = []

    for d in dirs:
        ls = next(linecycler)
        ds = data.open_dataset(d+'/ipec_profile_output_n1.nc')
        datasets.append(ds)

        # 1D plots
        lines = []
        for m in range(1,6):
            for ax,key in zip(axes,['b_n','xi_n']):
                l, = np.abs(ds[key]).sel(m=m).plot(ax=ax)
                l.set_linestyle(ls)
                if d==dirs[0]:
                    l.set_label('m = {:}'.format(m))
                lines.append(l)
        sm = plt.set_linearray(lines,cmap='viridis')

    for ax in axes:
        ax.set_title('')
        leg = ax.legend()
    plt.show()

    return datasets

############################################ Run

if __name__ == "__main__":
    plot_xbnormal(*sys.argv[1:])