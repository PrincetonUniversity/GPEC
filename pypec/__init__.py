#!/usr/local/env python
"""
:mod:`pypec` -- Main Package
============================

Python modules for data visualization and postprocessing, 
as well as python wrappers for GPEC and PENT.

:Author: 
  N.C. Logan
:Location:
  Princeton Plasma Physics Laboratory
:Email:
  nlogan@pppl.gov


Python at PPPL 
==============

To set up python on portal at PPPL, insert the 
following lines into your .bashrc file:

  export PYTHONPATH=$PYTHONPATH:/path/to/gpec
  module load anaconda/2.3.0
  
..note: Anaconda is the most complete and current python distribution
available at pppl. Users can also use locally built python/2.7.2, but
will loose 3D plotting capabilities and may experience problems reading
large data files.

Obviously, "/path/to/gpec" needs to be replaced with the path
to your installation of the GPEC package (containing the pypec directory).

For tcsh shell users, replace export with setenv syntax
and insert in .cshrc.

..note: Mayavi's default gui api is not consistent with the anaconda
interactive python default. To enable 3D plotting, you must add
"export QT_API=pyqt" to your .bashrc or a similar line to your .cshrc.

To put the changes into effect:

  $ source ~/.bashrc

Users are futher encouraged to create a matplotlib rc file in
user ~/.config/matplotlib similar to /u/nlogan/.config/matplotlib/matplotlibrc.
This file sets a number of plotting defaults.


Python Tutorials
-----------------

The 3 workhorse modules for scientific programming
in python are numpy, scipy and matplotlib (although
the third is really a personal preference). Tutorials
for each can be found at

`numpy tutorial <http://wiki.scipy.org/Tentative_NumPy_Tutorial>`_

`scipy tutorial <http://docs.scipy.org/doc/scipy/reference/tutorial/>`_

`matplotlib tutorial <http://matplotlib.org/users/pyplot_tutorial.html>`_


Setting Your Personal Defaults
------------------------------

To set personal defualts for runs, edit the _defaults_.py file in
this directory. Some of these (like email) should be obvious, others
you may not understand until becoming more familiar with the package.
Do not worry if these are confusing, and feel free to skip them for now.


Using the Best Environement
---------------------------

The ipython (interactive) environment is recomended for commandline
data analysis. To start, it is easiest to auto-import a number of 
basic mathematic and visualization modules and use an enhanced 
interactivity for displaying figures. To do this enter:

  $ ipython --pylab

Look into online tutorial on numpy and matplotlib. Advanced users 
are recomended to edit ~/.matplotlib/matplotlibrc and 
~/.ipython/profile_default/startup/autoimports.ipy files.

Now in the ipython envirnment, type

>>> from pypec import gpec,data # doctest:+ELLIPSIS

There may be a few warnings/hints... take them or leave them.


"""

# This file tells python to treat the folder as a package

# for "from package import *" to work use __all__ = ["file1","file2",...]
__all__ = ['data','gpec']#,'namelist','modplot']

import gpec,data,synthetics,modplot
import numpy as np

