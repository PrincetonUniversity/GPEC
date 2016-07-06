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
import sys, inspect, itertools
from pypec import data
from pypec import modplot as plt

############################################ Variables

linestyles = ["-","--","-.",":"]
linecycler = itertools.cycle(linestyles)

############################################ Functions

def check_profile(*directories,method='fgar'):
    """

    Check the ntv profile.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    datasets = []

    return datasets

def check_ell(*directories):
    """

    Check the ell dependence of the integrated torques.

    Args:
        *args: strings. Any number of directory paths containing ipec outputs.

    Returns:
        datasets: List of datasets from the corresponding profile ipec netcdfs

    """
    datasets = []

    return datasets

############################################ Run


def check_all(*directories):
    # generic way to collect all the functions
    this_module = sys.modules[__name__]
    all_functions = inspect.getmembers(this_module, inspect.isfunction)
    # call each function
    for key, function in all_functions.iteritems():
        print("Calling {:}".format(key))
        function(*directories)


if __name__ == "__main__":
    check_all(*sys.argv[1:])