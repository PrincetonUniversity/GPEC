#!/usr/local/env python
"""
:mod:`pypec.data` -- ASCII Table Data Turned Pythonic
=====================================================

The data module is a very generalized tool for visualizing and manipulating scientific data from netcdf and ascii files.

GPEC data is migrating to netcdf, and this module uses the xarray module for its python netcdf interface. The xarray module is very powerful and well documented `here <http://xarray.readthedocs.org/>`_.

Previous GPEC results and some current outputs are still writen to ascii. This module contains a custom I/O for converting ascii data to python objects. A custom Data object was created for this purpose, and conversion to the xarray Dataset is automated in the open_dataset function to facilitate migration to netcdf. The only requirement for creating a data object from an ascii file is that the file have one or more tables of data, with appropriate labels in a header line directly above the data.

Here, we show some examples using common outputs from both IPEC and PENT.

Beginners
---------

First, add your release of GPEC to your PYTHONPATH environment variable.

In an open ipython session, import the data module

>>> from pypec import data

Explore the data module by typing 'data.' and hitting tab. For information on any particular object in the module, use the '?' syntax or 'help' function. For example, enter 'data.read?'. You will see that it has a lot of metadata on this function, but the important things to look at for now are the arguments and key word arguments. These are passed to the function as shown in the Definition line... you can either pass a single string, or a string and specification for squeeze. If this confuses you, consult any of the many great python tutorials online.

Lets use it! To get a typical GPEC output into python, use the open_dataset function.

>>> con = data.open_dataset('examples/DIIID_example/ipec_control_output_n1.nc')

Printing to the terminal is intelegent, giving a conviniently human readable summary of the data.

>>> con
<xarray.Dataset>
Dimensions:        (i: 2, m: 129, mode: 34, mode_PX: 34, mode_RX: 34, mode_SC: 4, mode_WX: 34, mode_XT: 34, theta: 513)
Coordinates:
  * i              (i) int32 0 1
  * m              (m) int32 -64 -63 -62 -61 -60 -59 -58 -57 -56 -55 -54 -53 ...
  * mode           (mode) int32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 ...
  * mode_SC        (mode_SC) int32 1 2 3 4
  * theta          (theta) float64 0.0 0.001953 0.003906 0.005859 0.007812 ...
  * mode_XT        (mode_XT) int32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
  * mode_WX        (mode_WX) int32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
  * mode_RX        (mode_RX) int32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
  * mode_PX        (mode_PX) int32 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
Data variables:
    L              (mode, m) complex128 (-0.478146242223-0.0251176845733j) ...
    Lambda         (mode, m) complex128 (1.4032106893e-07-3.56355823737e-08j) ...
    P              (mode, m) complex128 (0.000212307471204-9.52839135614e-05j) ...
    rho            (mode, m) complex128 (9.96920996839e+36+9.96920996839e+36j) ...
    X_EDT          (mode_XT, m) complex128 (3.88116613338e-05+1.90243722796e-05j) ...
    X_EVT          (mode_XT) float64 2.517e+04 1.334e+05 5.081e+05 1.119e+06 ...
    W_EDX          (mode_WX, m) complex128 (2.0673520088e-06-4.52521668335e-07j) ...
    W_EVX          (mode_WX) float64 8.543 2.882 2.516 1.137 0.8107 0.6662 ...
    W_EVX_A        (mode_WX) float64 3.922e+05 7.805e+04 2.466e+05 3.042e+05 ...
    W_EVX_energyv  (mode_WX) float64 6.701e+06 4.499e+05 1.241e+06 6.917e+05 ...
    W_EVX_energys  (mode_WX) float64 6.383e+06 8.481e+05 1.724e+06 6.867e+05 ...
    W_EVX_energyp  (mode_WX) float64 6.817e+06 5.582e+05 8.381e+05 5.807e+05 ...
    R_EDX          (mode_RX, m) complex128 (-9.04352944467e-07-5.23229281269e-07j) ...
    R_EVX          (mode_RX) float64 1.881e+06 8.165e+05 2.477e+05 1.895e+05 ...
    R_EVX_RL       (mode_RX) float64 5.689 0.7149 0.3795 0.8092 0.2417 ...
    R_EVX_energyv  (mode_RX) float64 4.298e+05 1.373e+06 1.01e+06 7.832e+05 ...
    R_EVX_energys  (mode_RX) float64 8.994e+05 1.745e+06 1.096e+06 8.803e+05 ...
    R_EVX_energyp  (mode_RX) float64 3.785e+05 8.975e+05 8.31e+05 6.814e+05 ...
    P_EDX          (mode_PX, m) complex128 (-5.0711563448e-06-2.51239391938e-06j) ...
    P_EVX          (mode_PX) float64 493.5 161.8 146.4 128.7 127.8 110.2 ...
    O_WRX          (mode_RX, mode_WX) float64 0.082 0.99 0.08839 0.01132 ...
    O_WPX          (mode_PX, mode_WX) float64 0.01414 0.9607 0.08296 0.01058 ...
    O_RPX          (mode_PX, mode_RX) float64 0.9571 0.01583 0.02807 0.1558 ...
    O_XT           (mode_XT) complex128 (0.00296215168999+0.00489536349583j) ...
    O_WX           (mode_XT) complex128 (0.0101206436883-0.00397946238035j) ...
    O_RX           (mode_RX) complex128 (-0.00967090555989-0.000112969521773j) ...
    O_PX           (mode_PX) complex128 (-0.00598278876631+0.0105074013581j) ...
    Phi_EX         (m) complex128 (-1.42853908779e-10-4.20811059626e-11j) ...
    Phi_ET         (m) complex128 (-6.89310992668e-10-1.01978717443e-09j) ...
    W_EDX_FUN      (mode_WX, theta) complex128 (-1.45421344314-0.0472501327701j) ...
    R_EDX_FUN      (mode_RX, theta) complex128 (1.45665015126-0.521253285202j) ...
    b_nm           (m) complex128 (-9.46126704696e-10+1.48737377032e-09j) ...
    b_xnm          (m) complex128 (-5.48365827795e-10+5.89720875704e-10j) ...
    xi_nm          (m) complex128 (-3.3871234611e-09-1.86371281018e-09j) ...
    xi_xnm         (m) complex128 (-6.60388933097e-10+2.7635104984e-09j) ...
    xi_n           (theta) complex128 (0.000402066428626+0.00260098774146j) ...
    xi_xn          (theta) complex128 (-0.000601114809848-0.00197978522356j) ...
    b_n            (theta) complex128 (-0.00157772177652+0.00052698051739j) ...
    b_xn           (theta) complex128 (-0.000920387597747+0.000315158824532j) ...
    R              (theta) float64 2.289 2.289 2.288 2.287 2.286 2.284 2.282 ...
    z              (theta) float64 0.02442 0.03972 0.05506 0.07046 0.08596 ...
    R_n            (theta) float64 1.0 0.9995 0.9981 0.9959 0.9928 0.9889 ...
    z_n            (theta) float64 0.003568 0.03241 0.06134 0.09036 0.1195 ...
    dphi           (theta) float64 0.0 -0.039 -0.078 -0.117 -0.1559 -0.1948 ...
Attributes:
    title: IPEC outputs in Fourier or alternate modal bases
    Jacobian: pest
    q_lim: 5.2
    psi_n_lim: 0.991536136179
    shot: 147131
    time: 2300
    machine: DIII-D
    n: 1
    version: GPEC version 0.3.1
    long_name: Surface Inductance
    energy_vacuum: 3.28697635252
    energy_surface: 7.8714975121
    energy_plasma: 3.7597335618

We see the data contains dimensions, variables, and global attributes. Some individual data arrays may have attributes as well.

For quick visualization, we can use xarray's built in plotting routines,

>>> import xarray.plot as xplot
>>> f1,a1 = plt.subplots()
>>> l1 = xplot.plot(con['W_EVX']) # 1D data, returns a list of lines
>>> f1.savefig('examples/example_eigenvalues.png')

.. image:: examples/example_eigenvalues.png
   :width: 600px

>>> f2,a2 = plt.subplots()
>>> l2 = xplot.plot(con['W_EDX']) # This is 2D, returns Quadmesh
>>> f2.savefig('examples/example_eigenvectors.png')

.. image:: examples/example_eigenvectors.png
   :width: 600px

Note this uses matplotlib, and the plots are fully interactive. Navigate and zoom using the toolbar below plots.

This is not the recommended way to display detailed studies using GPEC data, however. It is
recommended that the user familiarise themselves with matplotlib (there are many great online
tutorials) and write their own scripts for plotting the data of interest, as details will
inevitably need tweaking.

For a common example of more customized plots, lets look at some 2d data

>>> cyl = data.open_dataset('examples/DIIID_example/ipec_cylindrical_output_n1.nc')

We could use the built in plotting methods like ds['b_z'].plot), but these are (R,z) figures
should be done with a little care to look right. Lets look at b_z for example,

>>> plt = data.plt
>>> f,ax = plt.subplots(figsize=(8,8))
>>> v = cyl['b_z']
>>> im = ax.imshow(np.real(v),aspect='equal',origin='lower',cmap='RdBu',vmin=-2e-3,vmax=2e-3,
...          extent=[v['R'].min(),v['R'].max(),v['z'].min(),v['z'].max()],
...          interpolation='gaussian')
>>> cb = f.colorbar(im)
>>> cb.set_label('$\delta B_z$')
>>> xtks = ax.set_xticks([1,1.7,2.4])
>>> xtxt = ax.set_xlabel('R (m)')
>>> ytxt = ax.set_ylabel('z (m)')
>>> f.savefig('examples/example_pbrzphi.png')

.. image:: examples/example_pbrzphi.png
   :width: 600px

There we go. That looks publication ready!

Note, one of the many cool additions to matplotlib in the pypec package is the ability to
write line plots to tables. Lets look at a 1D mode spectrum plot for example,

>>> prof = data.open_dataset('examples/DIIID_example/ipec_profile_output_n1.nc')
>>> f,ax = plt.subplots()
>>> for m in range(1,9):
...     lines = ax.plot(prof['psi_n'],np.abs(prof['xi_m_contrapsi'].sel(m=m)),label=r'm={:}'.format(m))
>>> mappable = plt.set_linearray(ax.lines,cmap='viridis')
>>> xtxt = ax.set_ylabel(r'$\\xi^\psi_m$')
>>> ytxt = ax.set_xlabel(r'$\psi$')
>>> leg = ax.legend(ncol=2)
>>> f.savefig('examples/example_xim.png')
>>> f.printlines('examples/example_xim.dat',squeeze=True)
Wrote lines to examples/example_xim.dat
True

.. image:: examples/example_xim.png
   :width: 600px

That's it! Read the table with your favorite editor. It will probably need a little
cleaning at the top since it tries to use the lengthy legend labels as column headers.


Advanced Users
--------------

In this section, we will assume you have a pretty good hold on both python and the basic ploting functions of these modules.

Now that you use this module regularly, save yourself some time by placing the line:

  from pypec import *

into your autoimport file located at
~/.ipython/profile_default/startup/autoimports.ipy


Deprecated Ascii Interface
---------------------------

As mentioned above, GPEC has moved the vast majority of its I/O to netcdf, and open_dataset should
be used to create xarray Dataset objects whenever possible. This section describes the deprecated
ascii interface for context, and motivates the move to netcdf.

Data in ascii tables can be read into a custom Data object using the read function,

>>> mydata = read('examples/DIIID_example/g147131.02300_DIIID_KEFIT.kin')
Casting table 1 into Data object.
>>> type(mydata),len(mydata)
(<type 'list'>, 1)

The read function creates a list of data objects corresponding to the tables in the text file. In this case, there is only one table. If we knew this ahead of time, we could have gotten the data object directly by splitting the list.

>>> mydata = read('examples/DIIID_example/g147131.02300_DIIID_KEFIT.kin')[0]
Casting table 1 into Data object.

At its heart, the data object consists of independent data and dependent data. This are stored differently, namely as a list in the pts attribute and as a dictionary in y. The x attribute is a dictionary if the data is 1 dimensional or a regular grid is found in the first n columns of the data (n=2,3), and left empty if the n>1 dependent data is irregular.

Explore the data using the corresponding python tools.

>>> mydata.xnames
['psi']
>>> mydata.y.keys()
['nim^3', 'tieV', 'teeV', 'wexbrad/s', 'nem^3']

You would access the ion density data using mydata.y['nim3'] for example.

The data object, however, is more that just a place holder for parced text files. It contains visualization and data manipulation tools. Plotting functions return interactive matplotlib Figure objects.

>>> fig = mydata.plot1d()
>>> fig.savefig('examples/example_profiles.png')

.. image:: examples/example_profiles.png
   :width: 600px

Interpolation returns a new Data object of the desired data, which can be one, multiple, or all of the dependent variables found in y.

>>> myinterp = mydata.interp(np.arange(0.1,1.0,0.1),['nim^3','tieV'])
Forming interpolator for nim^3.
Interpolating values for nim^3.
Forming interpolator for tieV.
Interpolating values for tieV.
>>> f = mydata.plot1d('nim^3')
>>> f = myinterp.plot1d('nim^3',marker='s',lw=0,figure=f)
>>> f.savefig('examples/example_profiles2.png')

.. image:: examples/example_profiles2.png
   :width: 400px

Note that the interpolator initialized the interpolation function on the first call, this is saved internally for subsequent calls.

One common need is to look at spectra. For this we want to utilize
the full functionality of the data instances' plot1d function.

>>> xc, = read('examples/DIIID_example/ipec_xclebsch_n1.out')
Casting table 1 into Data object.
>>> f = xc.plot1d('xi^psi','psi',x2rng=(1,3))
>>> f.savefig('examples/example_spectrum.png')

.. image:: examples/example_spectrum.png
   :width: 600px

We knew it was one table, so we used the "," to automatically
unpack returned list. We were then so familiar with the args and
kwargs of the plotting function that we did not bother typing the names.
If you do this, be sure you get the order consistent with the order
given in the documentation "Definition" line (watch out for inconsistent
"Docstring" listings).


.. note::
  These examples can be tested by developers using ipython 
  in this directory as follows:

    In [1]: import data,doctest
   
    In [2]: doctest.testmod(data,verbose=True)

"""
"""
    @package pypec
    @author NC Logan
    @email nlogan@pppl.gov
"""

import os,copy,time

from string import join,capitalize,count         # faster than '+' loop
from StringIO import StringIO
from types import MethodType

import numpy as np                               # math
from scipy.interpolate import interp1d,interp2d,LinearNDInterpolator,RegularGridInterpolator,griddata
from scipy.special import ellipk,ellipe

#in this package
import modplot as plt
import namelist 

try:
    import mayavi.mlab as mmlab
except ImportError:
    mmlab = False
    print('WARNING: Mayavi not in python path')
try:
    import xarray
except ImportError:
    xarray = False
    print('WARNING: xarray not in python path')
    print(' -> We recomend loading anaconda/2.3.0 on portal')
try:
    import seaborn as sns
except ImportError:
    sns = False
    print('WARNING: seaborn not in python path')
    print(' -> We recomend loading anaconda/2.3.0 on portal')

# better label recognition using genfromtxt
for c in '^|<>/':
    if c in np.lib._iotools.NameValidator.defaultdeletechars:
        np.lib._iotools.NameValidator.defaultdeletechars.remove(c)

######################################################## Global Variables

default_quiet = False
cmap_div = plt.pyplot.get_cmap('RdBu_r')
cmap_seq = plt.pyplot.get_cmap('viridis') #_load_default_cmap()

######################################################## Change default keywords

def interp1d_unbound(x, y, kind='linear', axis=-1, copy=True, bounds_error=False,
                     fill_value=np.nan, assume_sorted=False):
    return interp1d(x, y, kind=kind, axis=axis, copy=copy, bounds_error=bounds_error,
                     fill_value=fill_value, assume_sorted=assume_sorted)
interp1d_unbound.__doc__ = interp1d.__doc__

######################################################## Helper functions

def _set_color_defaults(calc_data,center=None,**kwargs):
    """
    Stolen from xarray.plot. Sets vmin, vmax, cmap.
    """
    vmin = kwargs.get('vmin',None)
    vmax = kwargs.get('vmax',None)
    cmap = kwargs.get('cmap',None)

    if vmin is None:
        vmin = np.percentile(calc_data, 2)
    if vmax is None:
        vmax = np.percentile(calc_data, 98)
    # Simple heuristics for whether these data should  have a divergent map
    divergent = ((vmin < 0) and (vmax > 0))
    # Now set center to 0 so math below makes sense
    if center is None:
        center = 0
    # A divergent map should be symmetric around the center value
    if divergent:
        vlim = max(abs(vmin - center), abs(vmax - center))
        vmin, vmax = -vlim, vlim
    # Now add in the centering value and set the limits
    vmin += center
    vmax += center
    # Choose default colormaps if not provided
    if cmap is None:
        if divergent:
            cmap = "RdBu_r"
        else:
            cmap = "viridis"

    kwargs['vmin'] = vmin
    kwargs['vmax'] = vmax
    kwargs['cmap'] = cmap

    return kwargs

######################################################## IO FOR DATA OBJECTs

def open_dataset(filename_or_obj,complex_dim='i',**kwargs):
    """
    Wrapper for xarray.open_dataset that allows automated reduction of
    a dimension destinguishing real and imaginary components.

    New Parameter
    -------------
    complex_dim : str, Dimension designating real/imaginary (0,1)

    """

    try:
        ds = xarray.open_dataset(filename_or_obj,**kwargs)
    except:
        dat = read(filename_or_obj,**kwargs)
        ds = xarray.Dataset()
        for i,d in enumerate(dat):
            if i==0:
                for k,v in d.params.iteritems():
                    ds.attrs[k]=v
            if not d.x:
                raise ValueError("No regular grid for dataset")
            for yk,yv in d.y.iteritems():
                ds[yk] = xarray.DataArray(yv.reshape(d.shape),coords=d.x,dims=d.xnames,attrs=d.params)

    if complex_dim in ds.dims:
        for k,v in ds.data_vars.iteritems():
            if complex_dim in v.dims:
                tmp = v.attrs
                ds[k] = v.loc[{complex_dim:0}]+1j*v.loc[{complex_dim:1}]
                ds[k].attrs = tmp

    return ds
open_dataset.__doc__+= xarray.open_dataset.__doc__

def read(fname,squeeze=False,forcex=[],forcedim=[],maxnumber=999,maxlength=1e6,
         auto_complex=True,quiet=default_quiet):
    """
    Get the data from any ipec output as a list of python 
    class-type objects using numpy.genfromtxt.

    Arguments
      fname      : str. 
        Path of ipec output file.
    
    Key Word Arguments:
      squeeze : bool. 
        Sets all attributes to single data object.
      forcex  : list. 
        Set independent data labels.
      forcedim : list.
        Set dimensionality of each x.
      maxnumber : int.
        Reads only first maxnumber of tables.
      maxlength : int.
        Tables with more rows than this are downsampled for speed.
      auto_complex : bool.
        Include amplitude and phase of complex data as variables.
      quiet   : bool.
        Prevent non-warning messages printed to terminal.
    
    Returns:
      list. 
        Class objects for each block of data in the file.
    
    """
    #debug start_time = time.time()

    collection = []
    dcollection = []
    pcollection = []
    
    with open(fname) as f:
        #Get preamble
        preamble = ''
        lastline = ''
        for lstart,line in enumerate(f):
            try:
                test = float(line.split()[0])
                break
            except:
                preamble+= lastline
                lastline = line
        f.seek(0)
        try: #Simple file : single table
            #data = np.genfromtxt(f,skip_header=lstart-1,names=True,dtype=np.float)
            raise ValueError # genfromtxt read pfile (multiple 2-column tables) as one... '\r' bug
            pcollection.append(preamble)
            dcollection.append(data)
        #look through file manually
        except ValueError as e:
            data = []
            f.seek(0)
            # clean up messy dcon output formats
            if 'dcon.out' in fname:
                lines = f.read()
                lines=lines.replace('mu0 p','mu0p')
                lines=lines.replace(' abs ',' abs')
                lines=lines.replace(' re ',' real')
                lines=lines.replace(' im ',' imag')
                lines=lines.replace('*','')
                lines = StringIO(lines).readlines()
            else:
                lines = f.readlines()
            top,bot = 0,0
            count = 0
            length= len(lines)
            while bot<length and top<(length-2) and count<maxnumber:
                preamble=''
                lastline=''
                for i,line in enumerate(lines[bot:]):
                    try:
                        # Find top defined as start of numbers
                        test = float(line.split()[0])
                        top = bot+i
                        # Throw out single lines of numbers 
                        #if not lines[top+1].translate(None,' \n\t\r'):#=='\n': 
                        #    raise ValueError
                        # Find bottom defined by end of numbers
                        for j,bline in enumerate(lines[top:]):
                            try:
                                test = float(bline.split()[0])
                            except:
                                break # end of numbers
                        # Throw out single lines of numbers
                        if j==1:
                            # but include them as preamble (DCON one-liners)
                            vals = lines[top]
                            keys = lines[top-1]
                            if not keys.translate(None,' \n\t'): #empty line
                                keys = lines[top-2]
                            if '=' not in keys and len(keys.split())==len(vals.split()):
                                for k,v in zip(keys.split(),vals.split()):
                                    preamble+='{:} = {:}\n'.format(k,v)
                            raise ValueError
                        else:
                            bot = top+j+1
                        break
                    except:
                        preamble+= lastline
                        lastline = line
                if line==lines[-1] and line==lastline:
                    break #end of file without another table
                """
                try:
                    bot = top+lines[top:].index(' \n')
                except ValueError:
                    try:
                        bot = top+lines[top:].index('\n')
                    except ValueError:
                        try:
                            bot = top+lines[top:].index('\r\n')
                        except:
                            try:
                                bot = top+lines[top:].index(' \r\n')
                            except:
                                bot = length
                """
                
                # include headers
                top-=1
                if not lines[top].translate(None,' \n\t'): #empty line
                    top-=1
                skipfoot = length-bot
                f.seek(0)
                table = lines[top:bot]
                if '\n' in table: #empty space
                    table.remove('\n')
                data = np.genfromtxt(StringIO(join(table)),names=True,
                                     deletechars='?',dtype=np.float)
                pcollection.append(preamble)
                dcollection.append(data)
                count+=1
        #debug print("Finished parsing file in "
        #debug       +"{} seconds".format(time.time()-start_time))
        #turn arrays into classes
        for i,(data,preamble) in enumerate(zip(dcollection,pcollection)): 
            if not quiet: 
                print("Casting table "+str(i+1)+" into Data object.")
            collection.append(DataBase(data,preamble,forcex=forcex,forcedim=forcedim,maxlength=maxlength,quiet=quiet))
            
    # force all attributes to single object
    if squeeze:
        dc1 = collection[0]
        if len(collection)>1:
            for dc in collection[1:]:
                dc1+=dc
        return dc1

    return collection


def readall(base='.',filetype='pent_n1.out',**kwargs):
    """
    Recursively searches base directory for files with name
    filetype, reads them into python objects, and stores 
    objects in a dictionary with directory names as keys.
    
    Key Word Arguments:
      base     : str. 
        Top level directory in which to begin search.
      filetype : str. 
        Files to read.
      kwargs :dict. 
        Passed to appropriate reader function.
    
    Returns:
      dict.
        Key-Value pairs are directories and their contained data.
    
    """

    results = {}
    if('.in' in filetype):
        reader = namelist.read
    else:
        reader = read
    base = os.path.abspath(base)+'/'
    #if(base[-1]!='/'): base+='/'
        
    subs = [name for name in os.listdir(base) if os.path.isdir(base+name)]
    if subs:
        results = {}
        for subd in subs:
            results[subd] = readall(base=base+subd,filetype=filetype,**kwargs)
    if os.path.isfile(base+filetype):
        results[filetype] = reader(base+filetype,**kwargs)
    _cleandict(results)

    # not a dict if only result in the directory
    if results.keys()==[filetype]:
        results = results[filetype]

    return results

def _cleandict(d):
    """ 
    Helper function to clear empty keys from dictionary.
    
    """
    ks = d.keys()
    for k in ks:
        if type(d[k])==dict:
            _cleandict(d[k])
        if not d[k]:
            del d[k]
    #return d


def write(dataobj,fname='',ynames=[],**kwargs):
    """
    Write data object to file. The aim is to replicate the original file
    structure... the header text will have been lost in the read however,
    as will some column label symbols in the use of genfromtxt.
    
    Arguments:
      dataobj : obj. 
        Data_Base type object.
      fname   : str. 
        File to write to.
    
      kwargs  : dict. 
      Passed to numpy.savetxt

    Returns:
      bool. 
        True.
    
    """
    if not fname: return False
    
    if 'fmt' not in kwargs: kwargs['fmt'] = '%16.8E'
    if 'delimiter' not in kwargs: kwargs['delimiter']=' '
    if 'comments' not in kwargs: kwargs['comments']=''
    if 'header' in kwargs:
        kwargs['header']+= '\n\n'
    else:
        kwargs['header']=''

    for k,v in dataobj.params.iteritems():
        kwargs['header']+= k+" = "+str(v)+"\n"

    nums = []
    if not ynames: ynames = dataobj.y.keys()    
    for i,name in enumerate(dataobj.xnames):
        nums.append(dataobj.pts[:,i])
        kwargs['header']+=kwargs['delimiter']*(i!=0)+'{:>16}'.format(name)
    for i,name in enumerate(ynames):
        if np.iscomplex(dataobj.y[name][0]):
            nums.append(np.real(dataobj.y[name]))
            kwargs['header']+=kwargs['delimiter']+'{:>16}'.format('real('+name.replace(' ','')+')')
            nums.append(np.imag(dataobj.y[name]))
            kwargs['header']+=kwargs['delimiter']+'{:>16}'.format('imag('+name.replace(' ','')+')')
        else:
            nums.append(dataobj.y[name])
            kwargs['header']+=kwargs['delimiter']+'{:>16}'.format(name.replace(' ',''))
    np.savetxt(fname,np.array(nums).T,**kwargs)

    return True




def plotall(results,xfun,yfun,label='',axes=None,**kwargs):
    """
    Line plot of data gethered recursively from a dictionary 
    containing multiple namelists or data objects.

    Arguments:
      results : dict. 
        A readall result.
      xfun    : func. 
        A function that takes args (key,val) for each dictionary item 
        and whose result is appended to the xaxis.
      yfun    : func. 
        A function that takes args (key,val) for each dictionary item 
        and whose result is appended to the yaxis.

    Key Word Arguments:
      axes    : obj. 
        matplotlib axes object.
      kwargs  : dict. 
        Passed to matplotlib plot function
    
    Returns:
      figure.
    
    """
    x = []
    y = []
    if not axes:
        f,ax = plt.subplots()
    else:
        ax = axes

    def newpoint(key,val):
        if type(val) == dict:
            for k,v in val.iteritems():
                newpoint(k,v)
        else:
            try:
                x.append(xfun(key,val))
                y.append(yfun(key,val))
            except:
                pass
    for key,val in results.iteritems():
        newpoint(key,val)

    z = zip(x,y)
    z.sort()
    x,y = zip(*z)
    ax.plot(x,y,label=label,**kwargs)
    ax.legend()
    f = ax.get_figure()
    f.show()
    return f


######################################################## THE BASE IPEC DATA OBJECTS

class DataBase(object):
    """
    An emtpy class object with the following operations overloaded::

      +,-,*,/,**

    such that the operation is only performed on the attributes 
    starting with .real_ and .imag_ and ending with r,z, or phi. 
    If the two factors in the operation have different r or z attributes
    data is interpolated to the wider spaced axis.
    
    """
    
    def __init__(self,fromtxt,preamble,forcex=[],forcedim=[],maxlength=1e6,
                 auto_complex=True,quiet=default_quiet):
        """
        Takes structured array from numpy.genfromtxt, and creates
        dictionaries of dependent and independent data. Also breaks
        down preamble to find any global parameters.
        
        **Arguments:**
          fromtxt : structured array. 
          preamble: str.
          
        **Key Word Arguments:**
          forcex  : list.
            x-axis labels.
          maxlength : int.
            Tables with more rows than this are downsampled for speed.
          auto_complex : bool.
            Include amplitude and phase of complex data as variables.
          quiet : bool.
            Supress printed information and warnings.

        """
        #debug start_time = time.time()
        names = list(fromtxt.dtype.names)
        l = len(fromtxt)
        potxnames = names[:3]
        # Set independent/dependent variables (max 3D,min 3pt grid)
        if not np.all([name in names for name in forcex]):
            print('WARNING: Requested axes not available. Using default left column(s).')
        if forcex and np.all([name in names for name in forcex]):
            nd = len(forcex)
            potxnames = forcex
        elif len(names)==2:
            nd = 1
            potxnames = names[0:1]
        else:
            nd = 1
            # hack for psi,q pairs
            if 'q' in potxnames:
                potxnames = names[:4]
                potxnames.remove('q')
            # hack for rzphi
            if 'l' in potxnames:
                potxnames = names[:4]
                potxnames.remove('l')
            if len(set(fromtxt[potxnames[1]]))<l/3: 
                nd+=1
                if len(set(fromtxt[potxnames[2]]))<l/3:
                    nd+=1
            # hack for pent
            if nd==1 and potxnames[1]=='Lambda':
                nd+=1
            if nd==2 and potxnames[2]=='x':
                nd+=1
        #debug print("Determined dimensionality in {} seconds".format(time.time()-start_time))
        #debug start_time = time.time()
        self.xnames= potxnames[:nd]    
        x = fromtxt[potxnames[:nd]]
        if nd==1:
            self.x = [x[self.xnames[0]]] # can be any order
        else:    # Must be written in ascending order
            self.x = [np.sort(list(set(x[name]))) for name in self.xnames]
        self.shape = [len(ax) for ax in self.x]
        #debug print("Set x in {} seconds".format(time.time()-start_time))
        #debug print("shape = "+str(self.shape))
        #debug start_time = time.time()        
        # reduce exesively large files, assuming 1st column is x
        if len(fromtxt)>maxlength:
            name = self.xnames[0]
            step = l/int(maxlength)+1
            if np.product(self.shape) != l:
                fromtxt = fromtxt[::step] # no axes to muck up
            else:
                mask = np.ones(self.shape,dtype=bool)
                mask[::step,...] = False
                #debug print('set mask in {} seconds'.format(time.time()-start_time))
                fromtxt = fromtxt[mask.ravel()]
            #debug print('masked fromtxt in {} seconds'.format(time.time()-start_time))
            if not quiet:
                print("WARNING: Reducing length of "
                +"table from {:} to {:}".format(l,len(fromtxt))
                +" by reducing x[0] by a factor of {}".format(step))
            l = len(fromtxt)
            x = fromtxt[potxnames[:nd]]
            self.x = [np.sort(list(set(x[name]))) for name in self.xnames]
            self.shape = [len(ax) for ax in self.x]
            #debug print("Set x in {} seconds".format(time.time()-start_time))
            #debug start_time = time.time()
        self.pts  = x.view(np.float).reshape(l,-1)
        if forcedim:
            #debug print(len(self.pts),self.pts.shape)
            forcedim = list(forcedim)+[-1]
            newx = self.pts.reshape(*forcedim)
            if nd==2: self.x = [newx[:,0,0],newx[0,:,1]]
            if nd==3: self.x = [newx[:,0,0,0],newx[0,:,0,1],newx[0,0,:,2]]
            self.shape = [len(ax) for ax in self.x]
        #debug print("Set pts in {} seconds".format(time.time()-start_time))
        #debug start_time = time.time()
        self.nd    = np.shape(self.pts)[-1]
        #debug print(self.shape,l)
        if np.product(self.shape) != l: #probably psi to closely packed
            # maybe first axis is packed too tight, try trusting others
            if self.nd ==1:
                self.shape = [l]
                self.x[0] = self.pts[:,0]
            if self.nd ==2:
                self.shape[0] = l/self.shape[1] 
                self.x[0]  = self.pts[:,0][::self.shape[1]]
            if self.nd ==3:
                self.shape[0] = l/np.product(self.shape[1:])
                self.x[0] = self.pts[:,0][::np.product(self.shape[1:])]
            # test result
            if np.any(list(set(self.pts[:,0]).difference(self.x[0]))):
                if not quiet: print("WARNING: Irregular dependent data: can not form axes for interpolation.")
                self.x=None
            else:
                if not quiet: print("Warning: Axes contain repeated values. Forced axes may be incorrect.")
        #debug print("Checked axes in {} seconds".format(time.time()-start_time))
        #debug start_time = time.time()

        self.y={}
        ynames = [name for name in names if name not in self.xnames]
        for name in ynames:
            if 'real' in name:
                newname = name.replace('real','')
                self.y[newname] = fromtxt[name]+1j*fromtxt[name.replace('real','imag')]
                if auto_complex:
                    self.y['|'+newname+'|'] = np.abs(self.y[newname]).real
                    self.y['angle '+newname] = np.angle(self.y[newname])
            elif not 'imag' in name:
                self.y[name]=fromtxt[name]
        #debug print("Set y in {} seconds".format(time.time()-start_time))
        #debug start_time = time.time()

        #scipy interpolation function for ND (stored so user can override)
        self._interpdict={}

        # Set any variables from preamble
        preamble=preamble.split()
        params=[]
        names = []
        while preamble.count('='):    #paramter specifications
            idx = preamble.index('=')
            param = preamble[idx+1]
            name = preamble[idx-1]
            name.translate(None,'()[]{}\/*-+^%.,:!@#&')
            if name in names and idx>1:
                name = join(preamble[idx-2:idx],' ')
                name.translate(None,'()[]{}\/*-+^%.,:!@#&')
            elif name in names:
                for sfx in range(1,100):
                    if name+str(sfx) not in names:
                        name+= str(sfx)
                        break
            try: 
                params.append((name,float(param))) #if its a number
            except: 
                params.append((name,param))        #if its not a number
            names.append(name)
            preamble.remove(preamble[idx-1])
            preamble.remove(str(param))
            preamble.remove('=')
        self.params=dict(params)
        self.__doc__ = join(preamble,'')    #whatever is left

        #debug print("Set preamble in {} seconds".format(time.time()-start_time))
        #debug start_time = time.time()


    # Operations overrides
    def __add__(self,other): return _data_op(self,other,np.add,'+')
    def __sub__(self,other): return _data_op(self,other,np.subtract,'-')
    def __mul__(self,other): return _data_op(self,other,np.multiply,'*')
    def __div__(self,other): return _data_op(self,other,np.divide,'/')
    def __pow__(self,other): return _data_op(self,other,pow,'^')
    def __radd__(self,other): return self.__add__(other)
    def __rmul__(self,other): return self.__mul__(other)
    def __rsub__(self,other): return self.__sub__(other)
    def __rdiv__(self,other): return self.__div__(other)
    

    # Built in interpolation
    def interp(self,x,ynames=None,quiet=False,**kwargs):
        """
        Interpolate data to point(s).
        
        Arguments:
          x  : ndarray shape (npts,ndim).
            Point(s) on dependent axis(axes).
        
        Key Word Arguments:
          ynames : list. 
            Any number of keys from y attribute
          quiet  : bool. 
            Suppres status messages

        Returns:
          dict. 
            Each name contains array of interpolated data.  
        
        """
        #housekeeping
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=(ynames,)
        if not type(x) in [list,tuple,np.ndarray]: x=[x,]
        NewData = copy.deepcopy(self)
        if self.nd==1:
            x = np.atleast_1d(x)
            NewData.pts = x
            NewData.x = [x]
            NewData.shape = [len(x)]
        else:
            if self.nd==2:
                x = np.atleast_2d(x)
            elif self.nd==3:
                x = np.atleast_3d(x)
            NewData.pts = x
            NewData.x = [np.array(sorted(set(xi))) for xi in x.T]
            NewData.shape = [len(xi) for xi in NewData.x]
            if np.product(NewData.shape)!=len(NewData.pts):
                NewData.x = None
        
        #function to form interpolator if needed
        def new_interp(name):
            if not quiet: print("Forming interpolator for "+name+".")
            args = []
            if self.x: # regular grid is fast
                for n in range(self.nd):
                    args.append(self.x[n])
                return RegularGridInterpolator(tuple(args),self.y[name].reshape(self.shape),bounds_error=False)
            else: # irregular grid is slower but ok
                step = max(1,len(self.pts[:,0])/1e6)
                for n in range(self.nd):
                    args.append(self.pts[::step,n])
                return LinearNDInterpolator(zip(*args),self.y[name][::step])

        # for each name check if interpolator is up to date and get values
        values={}
        for name in ynames:
            if name in self._interpdict.keys():
                if not all(self._interpdict[name].values.ravel() == self.y[name]):
                    self._interpdict[name] = new_interp(name)
            else:
                self._interpdict[name] = new_interp(name)
            if not quiet: print("Interpolating values for "+name+".")
            values[name]=self._interpdict[name](x)
        NewData.y = values
        return NewData
        
    
    # Built-in visualizations
    def plot1d(self,ynames=None,xname=None,x1rng=None,x2rng=None,x3rng=None,squeeze=False,
               **kwargs):
        """
        Line plots of from 1D or 2D data.
        
        Key Word Arguments:

          ynames : list
            Strings specifying dependent data displayed.
          xname  : string
            The independent axis if 2D data.
          x1rng  : optional
            Valid formats are:

              - x               Plot single nearest point x-axis.
              - (xmin,xmax)     Plots all data in range (min,max). 
              - (xmin,xmax,n)   Plots n evenly spaced points in range.

          x2rng  : . 
            Same as x1rng. Applicable to 2D or 3D data.
          x3rng  : . 
            Same as x2rng for 3D data only. Order of x is 
            xnames with xname placed in front if given.
          squeeze : bool
            Plots all ynames in same axis.
            
        kwargs : dict
            Valid matplotlib pyplot.plot keyword arguments.

        Returns:

          f : figure
            poloidal cross sections of plasma response

        """
        # housekeeping
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=(ynames,)
        if not xname: xname=self.xnames[0]

        indx = self.xnames.index(xname)
        if self.x:
            x=self.x[indx]
        else:
            x = self.pts[:,indx]

        # helper to sort through axes of ND data
        def maskmaker(ax,rng):
            if type(rng) in [float,int]: 
                mask = ax==ax[abs(ax-rng).argmin()]
            elif type(rng) in [tuple,list] and len(rng)==2:
                mask = (ax>=rng[0]) * (ax<=rng[1])
            elif type(rng) in [tuple,list] and len(rng)>2:
                print("Warning: mapping all axes to regular grid.")
                if self.nd==2:
                    x,ax = np.mgrid[x.min():x.max():x.size,rng[0]:rng[1]:rng[2]]
                    x = x[0]
                    mask = range(rng[2])
                if self.nd==3:
                    raise ValueError("Griding for 3D data is under developement.")
            else:
                mask = np.ones_like(ax,dtype=bool)
            return ax,mask
        
        # display data
        if 'figure' in kwargs:
            f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
            del kwargs['figure']
        elif 'axes' in kwargs:
            if type(kwargs['axes'])==np.ndarray:
                f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
            else:
                f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
            del kwargs['axes']
        else:
            if squeeze:
                nrow,ncol = 1,1
            else:
                nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
            f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=plt.rcp_size*(ncol,nrow))
        if squeeze: ax = np.array([ax.ravel()[0]]*len(ynames))

        if 'label' in kwargs:
            ulbl = kwargs['label']+', '
            del kwargs['label']
        else:
            ulbl = ''
        
        for a,name in zip(ax.ravel(),ynames):
            lbl = ulbl + _mathtext(name)
            if self.nd==1:
                a.plot(x,self.y[name],label=lbl,**kwargs)
            elif self.nd==2:
                if self.x:
                    x2,rng2 = maskmaker(self.x[indx-1],x2rng)
                    ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
                    xlbl = _mathtext(self.xnames[indx-1])
                    for y,x2 in zip(ys[rng2],x2[rng2]):
                        x,rng1 = maskmaker(x,x1rng)
                        a.plot(x[rng1],y[rng1],label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
                else:
                    x2s,rng2 = maskmaker(self.pts[:,indx-1],x2rng)
                    xlbl = _mathtext(self.xnames[indx-1])
                    x2s = np.sort(list(set(x2s[rng2])))
                    for x2 in x2s:
                        x,rng2 = maskmaker(self.pts[:,indx-1],float(x2))
                        x = self.pts[:,indx][rng2]
                        y = self.y[name][rng2]
                        x,y = np.array(zip(*sorted(zip(x,y))))
                        x,rng1 = maskmaker(x,x1rng)
                        a.plot(x[rng1],y[rng1],label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
            elif self.nd==3:
                indx23 = range(3)
                indx23.remove(indx)
                if self.x:
                    x2,rng2 = maskmaker(self.x[indx23[0]],x2rng)
                    x3,rng3 = maskmaker(self.x[indx23[1]],x3rng)
                    ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
                    xlbls = [_mathtext(name) for name in self.xnames if name!=xname]
                    for y2,x2 in zip(ys[rng2],x2[rng2]):
                        x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
                        for y,x3 in zip(y2[rng3],x3[rng3]):
                            x,rng1 = maskmaker(x,x1rng)
                            a.plot(x[rng1],y[rng1],label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
                else:
                    x2s,rng2 = maskmaker(self.pts[:,indx23[0]],x2rng)
                    x2s = np.sort(list(set(x2s[rng2])))
                    xs = self.pts[:,indx][rng2]
                    ys = self.y[name][rng2]
                    xlbls = [_mathtext(xlbl) for xlbl in self.xnames if xlbl!=xname]
                    for x2 in x2s:
                        x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
                        x3s,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],x3rng)
                        x3s = np.sort(list(set(x3s[rng3])))
                        for x3 in x3s:
                            x3d,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],float(x3))
                            x = xs[rng3]
                            y = ys[rng3]
                            np.array(zip(*sorted(zip(x,y))))
                            x,rng1 = maskmaker(x,x1rng)
                            a.plot(x[rng1],y[rng1],label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
            a.set_xlabel(_mathtext(self.xnames[indx]))
            a.legend(ncol=max(1,len(a.lines)/6))
        f.show()
        return f
                              
    def plot2d(self,ynames=None,aspect='auto',plot_type='imshow',cbar=True,
               center=None,grid_size=(256,256),swap=False,**kwargs):
        """
        Matplotlib 2D plots of the data.
        
        Key Word Arguments:
          ynames : list. 
            Strings specifying dependent data displayed.
          aspect : str/float. 
            Choose from 'auto' or 'equal' (see matplotlib axes.set_aspect).
            A float stretches a circle to hight=aspect*width. 
          plot_type: str. 
            Choose one of:
            
            - "imshow" Use matplotlib imshow to plot on grid.
            - "pcolormesh" Use matplotlib pcolormesh to plot on grid.
            - "contour"    Use matplotlib contour to plot on grid.
            - "contourf"   Use matplotlib contourf to plot on grid.
            - "tripcolor"  Use matplotlib tripcolor to plot raw points.
            - "tripcontour"  Use matplotlib tripcontour to plot raw points.
            - "tripcontourf" Use matplotlib tripcontourf to plot raw points.

          cbar   : bool. 
            Show colorbar.
          grid_size : tuple. 
            Size of grid used if no regular axes x (in order of xnames).
          swap   : bool. 
            Swap x/y axes.
          kwargs : dict. 
            Valid pcolormesh/contour/contourf key arguments.
        
        Returns:
          figure. 
        
        """
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=[ynames]
        if self.nd != 2: raise IndexError("Data not 2D.")
            
        if self.x:
            x1,x2 = map(np.reshape,self.pts.T,[self.shape]*self.nd)
        elif plot_type in ['imshow','pcolormesh','contour','contourf']:
            x1,x2 = np.mgrid[self.pts[:,0].min():self.pts[:,0].max():1j*grid_size[0],
                             self.pts[:,1].min():self.pts[:,1].max():1j*grid_size[1]]
        else:
            step = max(len(self.pts[:,0])/np.product(grid_size),1)
            if step>1: print('Reducing data by factor of {:}'.format(step))
            x1,x2 = self.pts[::step].T

        for name in ynames:
            try:
                if np.any(np.iscomplex(self.y[name])):
                    idx = ynames.index(name)
                    ynames[idx] = 'imag '+name
                    ynames.insert(idx,'real '+name)
            except: # looks at new names -> KeyError
                pass

        # display data
        if 'figure' in kwargs:
            f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
        elif 'axes' in kwargs:
            if type(kwargs['axes'])==np.ndarray:
                f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
            else:
                f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
        else:
            nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
            f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=plt.rcp_size*(ncol,nrow))
        
        for a,name in zip(ax.ravel(),ynames):
            print("Plotting "+name+".")
            plotter = getattr(a,plot_type)
            #if(plot_type=='line'):
            #    #plotter = getattr(a,'trip'*bool(self.x)+'contour')
            #    plotter=a.contour
            #elif(plot_type=='fill'):
            #    #plotter = getattr(a,'trip'*bool(self.x)+'contourf')
            #    plotter=a.contourf
            #else:
            #    if self.x:
            #        plotter = a.pcolormesh
            #    else:
            #        plotter = a.tripcolor

            # convert to reals and grid if necessary
            reducedname = name.replace('real ','').replace('imag ','')
            if 'imag ' in name:
                raw = np.imag(self.y[reducedname])
            else:
                raw = np.real(self.y[reducedname])
            if self.x:
                y = raw.reshape(self.shape)
            elif np.all(x1==self.pts[:,0]) or np.all(x1==self.pts[:,1]):
                y = raw
            else:
                y = griddata(self.pts,np.nan_to_num(raw),(x1,x2),method='linear')
            if swap:
                x1,x2,y = x2.T,x1.T,y.T

            kwargs = _set_color_defaults(y,center=center,**kwargs)

            # plot type specifics
            if plotter==a.imshow:
                kwargs.setdefault('origin','lower')
                kwargs['aspect']=aspect
                kwargs.setdefault('extent',[x1.min(),x1.max(),x2.min(),x2.max()])
                kwargs.setdefault('interpolation','gaussian')
                args = [y.T]
            else:
                args = [x1,x2,y]
                if plotter in [a.pcolormesh,a.tripcolor]:
                    kwargs.setdefault('edgecolor','None')
                    kwargs.setdefault('shading','gouraud')
            pcm = plotter(*args,**kwargs)

            a.set_xlabel(_mathtext(self.xnames[0+swap]))
            a.set_ylabel(_mathtext(self.xnames[1-swap]))
            a.set_xlim(x1.min(),x1.max())
            a.set_title(_mathtext(name))
            a.set_aspect(aspect)
            if cbar: f.colorbar(pcm,ax=a)
        f.show()
        return f

    def plot3d(self,ynames=None,filter={'psi':1},cbar=False,size=(600,600),
               plot_type='',phi=None,center=None,**kwargs):
        """
        Three dimensional plots. Data with xnames r,z will be plotted
        with third dimension phi where Y = Re[y]cos(n*phi)+Im[y]sin(n*phi).

            1D data > lines in a plane in 3-space.
            2D data > z axis becomes y value
            3D data > chosen display type

        Must have mayavi.mlab module available in pythonpath.
        
        *Key Word Arguments:*
          ynames   : list. 
            Strings specifying dependent data displayed.
          filter   : dict. 
            Filter points as only those closest to this y value.
          cbar : bool.
            Display a vertical colorbar in the Mayavi Scene.
          size : tuple.
            Size of figure in pixels.
          plot_type : str. 
            Valid mayavi.mlab function name. 
            Overrides default function choices (plot3d,mesh,quiver3d).
          phi : ndarray.
            Toroidal angle grid for complex 2D data. Default is 360 degrees, 180 pts.
          center : float.
            Center colormap on this value.

        *Returns:*
          figure.
          
          
        *Example:*
        
        xb, = data.read('examples/DIIID_example/ipec_xbnormal_fun_n1.out',forcex=['r','z'])
        xb.plot3d('bno')
        
        ..note: I think what we usually want out of IPEC is
        the 3D surface perturbation. If we run with fun_flag
        xn, = data.read('ipec_xbnormal_fun_n2.out')
        psi = 1.0
        n = 2
        fltr= xb.y['psi']==xb.y['psi'][np.abs(xb.y['psi']-psi).argmin()]
        r,z,p = xb.pts[:,0][fltr],xb.pts[:,1][fltr],np.linspace(0,2*np.pi,60)
        s = xb.y['bno'][fltr]
        cx = r.reshape(-1,1)*np.exp(1j*p.reshape(1,-1))
        X,Y = cx.real,cx.imag
        Z = z.reshape(-1,1)*np.ones_like(X)
        S = (s.reshape(-1,1)*np.exp(-n*1j*p.reshape(1,-1))).real
        mmlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
        mmlab.mesh(X, Y, Z, scalars=S) #, colormap='YlGnBu'
        
        """
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=[ynames]
        if not mmlab: raise ImportError("Unable to import mayavi.mlab.")

        # the correct 3d plotting function
        plotfunc = getattr(mmlab,['plot3d','mesh','quiver3d'][self.nd-1])
        if plot_type:
            if plot_type=='volume':
                def plotfunc(x,y,z,s,name='',**kwargs):
                    mmlab.pipeline.sc
                    src = mmlab.pipeline.scalar_field(xx,yy,zz,S,name=name)
                    mmlab.pipeline.volume(src,**kwargs)
                    mmlab.pipeline.scalar_cut_plane(thr,plane_orientation='x_axes')
                    mmlab.pipeline.scalar_cut_plane(thr,plane_orientation='y_axes')
            else:
                plotfunc = getattr(mmlab,plot_type)
                

        # filter data by other y variable
        fltr = self.y[ynames[0]]==self.y[ynames[0]]
        for k in filter:
            if k in self.y.keys():
                fltr*= self.y[k]==self.y[k][np.abs(self.y[k]-filter[k]).argmin()]
        
        # general multidimensional axes
        X = np.array(self.pts[fltr,:]).T.tolist()
        if self.nd==1: X += np.zeros_like(X.ravel())
        # Extra work to convert the cylindrical axes
        if self.xnames==['r','z']:
            r,z,p = self.pts[:,0][fltr],self.pts[:,1][fltr],np.linspace(0,2*np.pi,180)
            if phi!=None: p = np.array(phi) # Added for more custamization
            XY = r.reshape(-1,1)*np.exp(1j*p.reshape(1,-1))
            Z = z.reshape(-1,1)*np.ones_like(XY.real)
            X = [XY.real,XY.imag,Z]
        
        # display data
        for name in ynames:
            if 'figure' not in kwargs: f = mmlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=size)
            S = self.y[name][fltr]
            if self.xnames==['r','z']:
                S = np.real((S.reshape(-1,1)*np.exp(-self.params['n']*1j*p.reshape(1,-1))))
            if plotfunc in [mmlab.mesh,mmlab.points3d]:
                XYZ = X
                kwargs['scalars'] = S
            else:
                XYZ = X+np.real(S).tolist()
            #print(np.array(XYZ).shape)
            kwargs = _set_color_defaults(S,center=center,**kwargs)
            mapper = {'viridis':'YlGnBu'}
            cmap = kwargs.pop('cmap','')
            cmap = mapper.get(cmap,cmap)
            reverse = cmap.endswith('_r')
            kwargs['colormap'] = cmap.rstrip('_r')
            plotobj = plotfunc(*XYZ,name=name,**kwargs)
            if reverse: plotobj.module_manager.scalar_lut_manager.reverse_lut = True
            if cbar: mmlab.colorbar(title=name,orientation='vertical')
            f = mmlab.gcf()
        return f #X,S
    

    def scatter(self,ynames=None,xname=None,x1=None,x2=None,x3=None,**kwargs):
        """
        Scatter plot display for non-regular data. Display mimics plot1d in t
        hat it displays a single independent axis. If is another independent 
        axis exists, and is not sliced using one of the optional keyword arguments, 
        it s represented as the color of the plotted symbols. If the data is 3D 
        and unsliced, preference if given in order of xnames.

        Arguments:     
          ynames : list.
            Dependent data to be displayed. Default is all y data.
          xname  : str.
            Independent axis used in scatter plots. Default is first in xnames.
          x1     : float.
            Use only data with first axis clossest to value.
          x2     : float.
            Use only data with second axis clossest to value.
          x3     : float.
            Use only data with thrid axis clossest to value.
        
        Returns:
          figure
        
        """
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=[ynames]
        if not xname: xname=self.xnames[-1]
        indx = self.xnames.index(xname)

        rng = np.ones_like(self.pts[:,0],dtype=bool)
        ylbl=''
        indx2 = range(self.nd)
        indx2.remove(indx)
        if x1: 
            indx2.remove(0)
            x1 = self.pts[:,0][abs(self.pts[:,0]-x1).argmin()]
            ylbl+=', '+self.xnames[0]+'={0:.3}'.format(x1)
            rng*=self.pts[:,0]==x1
        if x2: 
            indx2.remove(1)
            x2 = self.pts[:,1][abs(self.pts[:,1]-x2).argmin()]
            ylbl+=', '+self.xnames[1]+'={0:.3}'.format(x2)
            rng*=self.pts[:,1]==x2
        if x3: 
            indx2.remove(2)
            x3 = self.pts[:,2][abs(self.pts[:,2]-x3).argmin()]
            ylbl+=', '+self.xnames[2]+'={0:.3}'.format(x3)
            rng*=self.pts[:,2]==x3
        if indx2: 
            c = self.pts[:,indx2[0]][rng]
            cset = np.array(list(set(c))).round(4)
            clbl = self.xnames[indx2[0]]
        else:
            c = np.ones_like(rng[rng])
            clbl = ''

        #double number of axes for complex data
        for name in ynames:
            try:
                if np.any(np.iscomplex(self.y[name])):
                    idx = ynames.index(name)
                    ynames[idx] = 'imag '+name
                    ynames.insert(idx,'real '+name)
            except: # looks at new names -> KeyError
                pass
        # actual plotting
        if 'figure' in kwargs:
            f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
        elif 'axes' in kwargs:
            if type(kwargs['axes'])==np.ndarray:
                f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
            else:
                f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
        else:
            nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
            f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=plt.rcp_size*(ncol,nrow))

        for a,name in zip(ax.ravel(),ynames):
            """if 'real' in name:
                sc=a.scatter(self.pts[:,indx][rng],np.real(self.y[name.replace('real ','')][rng]),
                         c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
            elif 'imag' in name:
                sc=a.scatter(self.pts[:,indx][rng],np.imag(self.y[name.replace('imag ','')][rng]),
                         c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
            else:
                sc=a.scatter(self.pts[:,indx][rng],np.real(self.y[name][rng]),
                         c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
        """
            y = self.y[name.replace('real ','').replace('imag ','')][rng]
            if 'real' in name:
                y = np.real(y)
            elif 'imag' in name:
                y = np.imag(y)
            sc=a.scatter(self.pts[:,indx][rng],y,c=c,
                         cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
            # for printing purposes
            line = plt.matplotlib.lines.Line2D(
                    self.pts[:,indx][rng],y,visible=False)
            if 'label' in kwargs: line.set_label(kwargs['label'])
            a.add_line(line)
            
            a.set_xlabel(_mathtext(self.xnames[indx]))
            a.set_ylabel(_mathtext(name+ylbl))
            if clbl:
                f.colorbar(sc,ax=a,ticks=cset)
                cb=f.get_axes()[-1]
                cb.set_ylabel(_mathtext(clbl),fontsize='large')
        f.show()
        return f

    def slice(self,x,xname=None,ynames=[]):#npts=100,**kwargs):
        """
        UNDER CONSTRUCTION

        Return 1D data object from 2D. Method grids new data, 
        so beware of loosing accuracy near rational surfaces.
        
        Arguments:    
          x       : float. 
            Lower bound of the axis.
            
        Key Word Arguments:
          xname   : str. 
            Axis allong which slice is performed.
        
        Returns: 
          obj. 
            New data nd-1 dimensional data object.
        
        """

        # housekeeping
        if not ynames: ynames=np.sort(self.y.keys()).tolist()
        if not type(ynames) in (list,tuple): ynames=(ynames,)
        if not xname: xname=self.xnames[0]

        indx = self.xnames.index(xname)
        x = self.pts[:,indx]

        # helper to sort through axes of ND data
        def maskmaker(ax,rng):
            if type(rng) in [float,int]: 
                mask = ax==ax[abs(ax-rng).argmin()]
            elif type(rng) in [tuple,list] and len(rng)==2:
                mask = (ax>=rng[0]) * (ax<=rng[1])
            elif type(rng) in [tuple,list] and len(rng)>2:
                print("Warning: mapping all axes to regular grid.")
                if self.nd==2:
                    x,ax = np.mgrid[x.min():x.max():x.size,rng[0]:rng[1]:rng[2]]
                    x = x[0]
                    mask = range(rng[2])
                if self.nd==3:
                    raise ValueError("Griding for 3D data is under developement.")
            else:
                mask = np.ones_like(ax,dtype=bool)
            return ax,mask
        
        for a,name in zip(ax.ravel(),ynames):
            if self.nd==1:
                a.plot(x,self.y[name],label=lbl,**kwargs)
            elif self.nd==2:
                if self.x:
                    x2,rng2 = maskmaker(self.x[indx-1],x2rng)
                    ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
                    xlbl = _mathtext(self.xnames[indx-1])
                    for y,x2 in zip(ys[rng2],x2[rng2]):
                        a.plot(x,y,label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
                else:
                    x2s,rng2 = maskmaker(self.pts[:,indx-1],x2rng)
                    xlbl = _mathtext(self.xnames[indx-1])
                    x2s = np.sort(list(set(x2s[rng2])))
                    for x2 in x2s:
                        x,rng2 = maskmaker(self.pts[:,indx-1],float(x2))
                        x = self.pts[:,indx][rng2]
                        y = self.y[name][rng2]
                        x,y = zip(*sorted(zip(x,y)))
                        a.plot(x,y,label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
            elif self.nd==3:
                indx23 = range(3)
                indx23.remove(indx)
                if self.x:
                    x2,rng2 = maskmaker(self.x[indx23[0]],x2rng)
                    x3,rng3 = maskmaker(self.x[indx23[1]],x3rng)
                    ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
                    xlbls = [_mathtext(name) for name in self.xnames if name!=xname]
                    for y2,x2 in zip(ys[rng2],x2[rng2]):
                        x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
                        for y,x3 in zip(y2[rng3],x3[rng3]):
                            a.plot(x,y,label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
                else:
                    x2s,rng2 = maskmaker(self.pts[:,indx23[0]],x2rng)
                    x2s = np.sort(list(set(x2s[rng2])))
                    xs = self.pts[:,indx][rng2]
                    ys = self.y[name][rng2]
                    xlbls = [_mathtext(xlbl) for xlbl in self.xnames if xlbl!=xname]
                    for x2 in x2s:
                        x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
                        x3s,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],x3rng)
                        x3s = np.sort(list(set(x3s[rng3])))
                        for x3 in x3s:
                            x3d,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],float(x3))
                            x = xs[rng3]
                            y = ys[rng3]
                            x,y = zip(*sorted(zip(x,y)))
                            a.plot(x,y,label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
        
    
    





######################################################## helper functions

##class operations

def _data_op(self,other,fun,op,quiet=default_quiet):
    """
    This function generalizes the process of creating a new class
    and filling it with the appropriate dependent data, info, and
    indepenedent data from an operation on one or more data classes.
    
    """
    NewData = copy.deepcopy(self) #intialize new data object.

    if type(other)==DataBase:
        # check if same dependent variables
        if self.xnames != other.xnames:
            raise ArithmeticError('Two data objects do not have same dependent variables.')
        if not np.all(self.pts==other.pts):#[np.all(sx==ox) for sx,ox in zip(self.x,other.x)]):
            if not quiet: print('interpolating to same axes.')
            otherinterp = other.interp(NewData.pts).y
        else:
            otherinterp = other.y
        # problem from unbound interpolation?? Fixed when switched to returning objects.
        #if self.nd==1:
        #    for k in otherinterp:
        #        otherinterp[k] = otherinterp[k].ravel()
            
    
        for key in other.y:
            if key in NewData.y:
                #if not quiet: print(key+op+key)
                NewData.y[key] = fun(NewData.y[key],otherinterp[key])
            else: 
                if not quiet: print('0'+op+key)
                NewData.y[key] = fun(np.zeros_like(otherinterp[key]),otherinterp[key])
        for key in NewData.y:
            if key not in other.y and not quiet: print(key+op+'0')
        
    else: # apply operation to y values
        for k,v in NewData.y.iteritems():
            #if not quiet: print(k+op+str(other))
            NewData.y[k] = fun(v,other)
            
    return NewData




def _mathtext(text):
    """
    Helper function to convert lazy ascii conventions to 
    nicer latex strings.
    
    Arguments:
      text : str.
    
    """

    simplemaps = ['epsilon','psi','omega','nabla','cdot','kappa','vartheta','theta',
                  'nu','times','perp','parallel','lambda','ln',
                  'Lambda','xi','chi','varphi','phi','Phi','mu','int', 'Re','Im',
                  'delta','angle','ell','alpha','beta','zeta','rho']
    
    trickymaps = [('par','parallel'), ('divxprp','nablacdotxi_perp'),
                  ('exb','EtimesB'), ('lmda','lambda'),#('lam','lambda'),
                  ('wdian','omega_*n'), ('wdiat','omega_{*T}'),
                  ('wrot','omega_phi'), ('real','Re'),
                  ('imag','Im'), ('log','ln'), ('ln10','log_{10}'),
                  ('kxprp','kappacdotxi_perp'),
                  ('prp','perp'), ('int','int '),('eps_','epsilon_')
                 ]

    # replace some lazy data labeling with proper math
    for tm in trickymaps:
        text=text.replace(*tm)
    # tex up symbols
    for sm in simplemaps:
        text=text.replace(sm,"$\\"+sm+'$')
    # clean up known bugs
    text = text.replace('e$\\psi$lon','epsilon')
    text = text.replace('var$\\theta','vartheta')
    text = text.replace('var$\\phi','varphi')
    # deliniate start/stop of latexing
    for word in text.split(' '):
        #if '_' in word: # got too fancy (required spaces after subscripts)
        #    newword='$'+word.replace('_','_{').replace('$','')+'}'*word.count('_')+'$'
        if '^' in word or '_' in word:
            text=text.replace(word,'$'+word.replace('$',' ')+'$')
    #text.replace('$$','')
    
    return text#.encode('string-escape') #raw string

def getshot(path='.',full_name=False):
    """
    Find 6 digit shot number in path to directory
    or file.
    
    **Arguments:**
      
    **Key Word Arguments:**
      d : str.
        Path.
    **Returns:**
      str.
        Shot number if found, '' if not.
    
    """
    pth = os.path.abspath(path)
    shot = ''
    for i in range(len(pth)-5):
        try:
            shot = str(int(pth[i:i+6]))
            break
        except:
            pass
    if full_name:
        dirs = pth.split('/')
        for d in dirs:
            if shot in d:
                shot = d
                break
    return shot

def add_control_geometry(ds,overwrite=False):
    """
    Add geometric dimensions to dataset from ipec_control_output_n#.nc.

    **Arguments:**
        ds : Dataset. xarray Dataset opened from ipec_control_output_n#.nc

    **Key Word Arguments:**
        overwrite : bool. Overwrite geometric quantities if they already exist.

    **Examples:**

    After opening, and adding geometry,
    >>> ds = open_dataset('examples/DIIID_example/ipec_control_output_n1.nc')
    >>> ds = add_control_geometry(ds)
    
    it is easy to make 3D surface plots using mayavi,
    >>> mesh = mmlab.mesh(ds['X'].values,ds['Y'].values,ds['Z'].values,scalars=np.real(ds['b_n']*ds['expn']))

    Make 2D perturbed last-closed-flux-surface plots quickly using,
    >>> fac = 1e-1 # good for DIII-D unit eigenmodes
    >>> f,ax = plt.subplots()
    >>> for i in range(4):
    ...     v = np.real(ds['W_EDX_FUN'][0,i]*fac)
    ...     l, = ax.plot(ds['R']+v*ds['R_n'],ds['z']+v*ds['z_n'])
    >>> sm = plt.set_linearray(ax.lines[-4:],range(1,5))
    >>> cb = f.colorbar(sm)
    >>> cb.set_label('Mode Index')
    >>> txt = ax.set_xlabel('R (m)')
    >>> txt = ax.set_ylabel('z (m)')


    
    """
    # Toroidal angle
    if 'phi' not in ds or overwrite:
        phi = np.linspace(0,2*np.pi,180)
        phi = xarray.DataArray(phi,coords={'phi':phi})
        ds['phi'] = phi
        ds['expn'] = np.exp(-1j*ds.attrs['n']*ds['phi'])

    # 3D cartesian coords for 3D plots
    if 'X' not in ds or overwrite:
        xy = (ds['R']*np.exp(1j*ds['phi'])).to_dataset(name='xy')
        ds['X'] = xy.apply(np.real)['xy']
        ds['Y'] = xy.apply(np.imag)['xy']
        ds['Z'] = ds['z']*(1+0*ds['phi'])

    # Normal vectors
    if 'z_n' not in ds or 'R_n' not in ds or overwrite:
        dr = np.roll(ds['R'],1)-np.roll(ds['R'],-1)
        dz = np.roll(ds['z'],1)-np.roll(ds['z'],-1)
        norm = np.sqrt(dr**2+dz**2)
        ds['z_n'] = xarray.DataArray(dr/norm,coords=ds['theta'].to_dataset())
        ds['R_n'] =-xarray.DataArray(dz/norm,coords=ds['theta'].to_dataset())

    return ds

######################################################## Developer functions


def _write_ellipk(npts=3e2,fname='ellipk01.dat'):
    """
    Writes complete elliptic integral of the first 
    kind to file.

    .. note:: This is not the original routine. Could not find it.

    Key Word Arguments: 
      npts : int.
      fname: str.
    
    """
    # ellipk(1-1e-17)=inf
    x = np.log(1./1e-15)
    kappas = 1.0-np.exp(np.linspace(0,x,npts))/np.exp(x)
    with open(fname,'w') as f:
        f.write('{:.0f} \n'.format(npts))
        for k in kappas[::-1]:
            f.write('{:.21} '.format(k))
        f.write('\n')
        for k in kappas[::-1]:
            f.write('{} '.format(ellipk(k)))
    return True

def _write_ellipe(npts=3e2,fname='ellipe01.dat'):
    """
    Writes complete elliptic integral of the first 
    kind to file.
    
    .. note:: This is not the original routine. Could not find it.
    
    Key Word Arguments:
      npts : int.
      fname: str.
    
    """
    kappas = np.linspace(0,1,npts)
    with open(fname,'w') as f:
        f.write('{:.0f} \n'.format(npts))
        for k in kappas:
            f.write('{} '.format(k))
        f.write('\n')
        for k in kappas:
            f.write('{} '.format(ellipe(k)))
    return True
        
def _factors(n):
    return list(sorted(set(reduce(list.__add__, 
                      ([i, n//i] for i in range(1, int(n**0.5) + 1) 
                       if n % i == 0)))))

def _plot3dvolume(self,ynames=None,npts=124,cbar=True,plot_type='',**kwargs):
    """
    Under construction volume plotting.
    """
    
    # general multidimensional axes
    X = np.array(self.pts).T.tolist()
    if self.nd==1: X += np.zeros_like(X.ravel())
    # Extra work to convert the cylindrical axes
    if self.xnames==['r','z']:
        # regularly spaced grid in x,y,z
        rlim = self.pts[:,0].max()
        zlim = np.abs(self.pts[:,1]).max()
        xx,yy,zz = np.mgrid[-rlim:rlim:npts*1j,-rlim:rlim:npts*1j,-zlim:zlim:npts*1j]
        rr 	 = np.sqrt(xx**2.0+xx**2.0)
        pp = np.arctan2(yy,xx)
        X = [xx,yy,zz]
    
        # To index coordinates
        cordpts = [(npts-1)*(rr-rlim)/(2*rlim),(npts-1)*(zz-zlim)/(2*zlim)]
    
    # display data
    for name in ynames:
        if 'figure' not in kwargs: f = mmlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=size)
        S = self.y[name]
        if self.xnames==['r','z']:
            
            # 2D interpolation
            print('interpolating real component...')
            S = griddata(self.pts,np.real(self.y[name]),(rr,zz),method='linear')
            print('interpolating imag component...')
            S+= griddata(self.pts,np.imag(self.y[name]),(rr,zz),method='linear')*1j
            S = np.real(S*np.exp(-self.params['n']*1j*pp))
            print(np.shape(S),np.shape(rr))
            # Take out wedge for better view
            wedge = (xx>=0.0)*(yy>=0.0)
            #m[wedge==False] = 0.0
            S[wedge==False] = -1e99
            S[(S==S)==False]= -1e99
            src = mmlab.pipeline.scalar_field(xx,yy,zz,S)
            thr = mmlab.pipeline.threshold(src,low=-1e98)
            #mmlab.pipeline.volume(src)
            mmlab.pipeline.scalar_cut_plane(src,plane_orientation='x_axes')
            mmlab.pipeline.scalar_cut_plane(src,plane_orientation='y_axes')
            if cbar: mmlab.colorbar(title=name,orientation='vertical')
        
    f = mmlab.gcf()
    return xx,yy,zz,S

def _plot3d_dev(self,ynames=None,filter={'psi':1},cbar=False,plot_type='',**kwargs):
    """
    Three dimensional plots. Data with xnames r,z will be plotted
    with third dimension phi where Y = Re[y]cos(n*phi)+Im[y]sin(n*phi).

        1D data > lines in a plane in 3-space.
        2D data > z axis becomes y value
        3D data > chosen display type

    Must have mayavi.mlab module available in pythonpath.
    
    *Key Word Arguments:*
      ynames   : list. 
        Strings specifying dependent data displayed.
      filter   : dict. 
        Filter points as only those closest to this y value.
      plot_type : str. 
        Valid mayavi.mlab function name. 
        Overrides default function choices (plot3d,mesh,quiver3d).

    *Returns:*
      figure.
      
      
    *Example:*
    
    xb, = data.read('examples/DIIID_example/ipec_xbnormal_fun_n1.out',forcex=['r','z'])
    xb.plot3d('bno')
    
    ..note: I think what we usually want out of IPEC is
    the 3D surface perturbation. If we run with fun_flag
    xn, = data.read('ipec_xbnormal_fun_n2.out')
    psi = 1.0
    n = 2
    fltr= xb.y['psi']==xb.y['psi'][np.abs(xb.y['psi']-psi).argmin()]
    r,z,p = xb.pts[:,0][fltr],xb.pts[:,1][fltr],np.linspace(0,2*np.pi,60)
    s = xb.y['bno'][fltr]
    cx = r.reshape(-1,1)*np.exp(1j*p.reshape(1,-1))
    X,Y = cx.real,cx.imag
    Z = z.reshape(-1,1)*np.ones_like(X)
    S = (s.reshape(-1,1)*np.exp(-n*1j*p.reshape(1,-1))).real
    mmlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    mmlab.mesh(X, Y, Z, scalars=S) #, colormap='YlGnBu'
    
    """
    if not ynames: ynames=np.sort(self.y.keys()).tolist()
    if not type(ynames) in (list,tuple): ynames=[ynames]
    if not mmlab: raise ImportError("Unable to import mayavi.mlab.")

    # the correct 3d plotting function
    plotfunc = getattr(mmlab,['plot3d','mesh','quiver3d'][self.nd-1])
    if plot_type:
        if plot_type=='volume':
            def plotfunc(X,Y,Z,S,name='',**kwargs):
                src = mmlab.pipeline.scalar_field(X,Y,Z,S,name=name)
                thr = mmlab.pipeline.threshold(src,low=1e-9)
                mmlab.pipeline.volume(src,**kwargs)
                mmlab.pipeline.scalar_cut_plane(thr,plane_orientation='x_axes')
                mmlab.pipeline.scalar_cut_plane(thr,plane_orientation='y_axes')
        else:
            plotfunc = getattr(mmlab,plot_type)
            

    # filter data by other y variable
    fltr = self.y[ynames[0]]==self.y[ynames[0]]
    for k in filter:
        if k in self.y.keys():
            fltr*= self.y[k]==self.y[k][np.abs(self.y[k]-filter[k]).argmin()]
    
    # general multidimensional axes
    X = np.array(self.pts[fltr,:]).T.tolist()
    if self.nd==1: X += np.zeros_like(X.ravel())
    # Extra work to convert the cylindrical axes
    if self.xnames==['r','z']:
        r,z,p = self.pts[:,0][fltr],self.pts[:,1][fltr],np.linspace(0,2*np.pi,180)
        XY = r.reshape(-1,1)*np.exp(1j*p.reshape(1,-1))
        Z = z.reshape(-1,1)*np.ones_like(XY.real)
        X = [XY.real,XY.imag,Z]
        if plot_type=='volume':
            X = np.meshgrid(XY.real,XY.imag,Z)
            x = [x[0],X[1],X[2]]
            p = np.arctan2(X[1],X[0])
    
    # display data
    for name in ynames:
        if 'figure' not in kwargs: f = mmlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=size)
        S = self.y[name][fltr]
        if self.xnames==['r','z']:
            S = np.real((S.reshape(-1,1)*np.exp(-self.params['n']*1j*p.reshape(1,-1))))
        if plotfunc in [mmlab.mesh,mmlab.points3d]:
            XYZ = X
            kwargs['scalars'] = S
        else:
            XYZ = X+[np.real(S)]
        #print(np.array(XYZ).shape)
        plotfunc(*XYZ,name=name,**kwargs)
        if cbar: mmlab.colorbar(title=name,orientation='vertical')
        f = mmlab.gcf()
    return XYZ
