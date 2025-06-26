#!/usr/local/env python
"""
:mod:`regression.compare` -- For Comparing Restults
=========================================================

Python modules for regression testing and benchmarking GPEC
code developments.

:Author:
  N.C. Logan
:Email:
  nikolas.logan@columbia.edu


Interface
-----------

The compare submodules can be run as an executable from a linux terminal, and take one or more directory paths as arguments. Each check function will be called for the directories and the plots displayed. Be warned, this may create a large number of plots.

For example,

  $ ./compare_dcons.py /p/gpec/GPEC-0.4/docs/examples/DIIID_example /p/gpec/GPEC-1.0/docs/examples/DIIID_ideal_example

Individual checks can be called from an interactive python environment (see :py:mod:`pypec`).

For example,

  $ ipython --pylab

>>> from regression.compare import compare_dcons
>>> data = compare_dcons.check_energies('/p/gpec/GPEC-0.4/docs/examples/DIIID_example','/p/gpec/GPEC-1.0/docs/examples/DIIID_ideal_example')

               /p/gpec/GPEC-0.4/d...              /p/gpec/GPEC-1.0/d...
========= ================================== ============================
total                   +7.3090e-01             +9.0210e-01+7.0180e-05j
plasma                  -1.7540e+00             -1.6000e+00
vacuum                  +2.4850e+00             +2.5020e+00
========= ================================== ============================

>>> pyplot.gcf().savefig('examples/figures/example_compare_dcons.png')

.. image:: examples/figures/example_compare_dcons.png
   :width: 600px

"""

# this file tells python to treat the folder as a package

# for "from package import *" to work use __all__ = ["file1","file2",...]
__all__ = ['compare_dcons','compare_gpecs','compare_pentrcs']

# calling "imprort package" will run this code
#import compare_dcons, compare_gpecs, compare_pentrcs
