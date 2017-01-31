#!/usr/local/env python
"""
:mod:`pypec.post` -- Post-Processing Fundamental GPEC Outputs
==============================================================

The post module is a repository of basic post-processing functions for expanding,
manipulating, and visualizing fundamental gpec output. This is not an exhaustive
list of post-processing procedures. The focus is on a) transforming between
coordinate systems and b) matrix optimizations.

Here, we show some examples using common outputs.

Examples
---------

First, add your release of GPEC to your PYTHONPATH environment variable.

In an open ipython session, import the data module

>>> from pypec import data, post

To get a typical GPEC output into python, use the open_dataset function.

>>> con = data.open_dataset('examples/DIIID_ideal_example/gpec_control_output_n1.nc')
>>> prof = data.open_dataset('examples/DIIID_ideal_example/gpec_profile_output_n1.nc')

First, 3D equilibrium calculations lend themselves naturally to 3D figures. This module
has a helpful function for forming the necessary x,y,z meshes from the 2D boundary.
After adding the requisite geometry to the control surface output,

>>> con = add_3dsurface(con, phi_lim=pi * 1.5)

it is easy to make 3D surface plots using mayavi,

>>> x,y,z = con['x_surf'].values,con['y_surf'].values,con['z_surf'].values
>>> b_n_3d = np.real(con['b_n_fun']*con['phase_surf'])
>>> fig = mmlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=(600,600))
>>> mesh = mmlab.mesh(x,y,z,scalars=b_n_3d,colormap='RdBu',vmin=-2e-3,vmax=2e-3,figure=fig)

Note that mayavi explicitly requires numpy.ndarray objects so we had to
use the values of each xarray DataArray.

Now, lets cap the ends with our full profile information. This is a little
more involved. For the toroidal angle of each end, we need to rotate the
2D arrays for the mesh and the field values.

>>> for phi in [0,np.pi*3/2]:
...     b = np.real(prof['b_n_fun']*np.exp(-1j*1*phi)) # poloidal cross section at phi
...     xy = prof['R']*np.exp(1j*phi) # rotate the plane
...     x,y,z = np.real(xy),np.imag(xy),np.real(prof['z'])
...     x,y,z,b = np.array([x,y,z,b])[:,:,::10] # downsample psi by 10 to speed up mayavi
...     cap = mmlab.mesh(x,y,z,scalars=b,colormap='RdBu',vmin=-2e-3,vmax=2e-3,figure=fig)
>>> mmlab.savefig(filename='examples/figures/example_bn_3d_capped.png',figure=fig)

.. image:: examples/figures/example_bn_3d_capped.png
   :width: 600px

That looks awesome!

There are, of course, plenty more outputs and important ways of visualizing
them intuitively. See, for example, the 2D plots in the add_3dsurface
function's documentation. At this point we have outlined the basics. Using
this module and similar commands to those used here, the user should be able
to do all sorts of post-processing tasks.


"""
"""
    @package pypec
    @author NC Logan
    @email nlogan@pppl.gov
"""

import numpy as np  # math
from collections import OrderedDict

# in this package
import modplot as plt
import data
from data import xarray, mmlab

######################################################## Global Variables

pi = np.pi
dtor = np.pi / 180


######################################################## Helper functions


def update_name_conventions(dataset, version=None, inplace=False):
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

    translator = OrderedDict()
    if version < 0.3:
        # profile output changes
        translator['xi_psi1'] = 'xigradpsi_dpsi'
        translator['xi_psi'] = 'xigradpsi'
        translator['xi_alpha'] = 'xigradalpha'
    if version < 0.5:
        # control output changes
        translator['b_xn'] = 'b_n_x_fun'
        translator['b_n'] = 'b_n_fun'
        translator['xi_xn'] = 'xi_n_x_fun'
        translator['xi_n'] = 'xi_n_fun'
        translator['dphi'] = 'delta_phi'
        translator['b_xnm'] = 'b_n_x'
        translator['b_nm'] = 'b_n'
        translator['xi_xnm'] = 'xi_n_x'
        translator['xi_nm'] = 'xi_n'
        translator['Phi_X'] = 'Phi_x'
        translator['Phi_EX'] = 'Phi_xe'
        translator['Phi_T'] = 'Phi'
        translator['Phi_ET'] = 'Phi_e'
        translator['X_EVT'] = 'X_eigenvalue'
        translator['X_EDT'] = 'X_eigenvector'
        translator['W_EVX'] = 'W_xe_eigenvalue'
        translator['R_EVX'] = 'R_xe_eigenvalue'
        translator['P_EVX'] = 'P_xe_eigenvalue'
        translator['C_EVX'] = 'C_xe_eigenvalue'
        translator['W_EDX'] = 'W_xe_eigenvector'
        translator['R_EDX'] = 'R_xe_eigenvector'
        translator['P_EDX'] = 'P_xe_eigenvector'
        translator['C_EDX'] = 'C_xe_eigenvector'
        translator['W_EVX_energyv'] = 'W_xe_energyv'
        translator['W_EVX_energys'] = 'W_xe_energys'
        translator['W_EVX_energyp'] = 'W_xe_energyp'
        translator['R_EVX_energyv'] = 'R_xe_energyv'
        translator['R_EVX_energys'] = 'R_xe_energys'
        translator['R_EVX_energyp'] = 'R_xe_energyp'
        translator['W_EVX_A'] = 'W_xe_amp'
        translator['R_EVX_RL'] = 'R_xe_RL'
        translator['O_XT'] = 'O_Xxi_n'
        translator['O_WX'] = 'O_WPhi_xe'
        translator['O_PX'] = 'O_PPhi_xe'
        translator['O_RX'] = 'O_RPhi_xe'
        # profile output changes
        translator['derxi_m_contrapsi'] = 'xigradpsi_dpsi'
        translator['xi_m_contrapsi'] = 'xigradpsi'
        translator['xi_m_contraalpha'] = 'xigradalpha'
        # cylindrical output changes
        translator['b_r_plas'] = 'b_r_plasma'
        translator['b_z_plas'] = 'b_z_plasma'
        translator['b_t_plas'] = 'b_t_plasma'

    # swap out any old names that are in the dataset
    #name_dict = OrderedDict(oldnew for oldnew in zip(old, new) if oldnew[0] in dataset)
    #newset = dataset.rename(name_dict, inplace=inplace)
    # do this in an explicit loop because some names already exist and need to get replaced in order
    for okey, nkey in translator.iteritems():
        if okey in newset:
            newset = newset.rename({okey:nkey}, inplace=True)

    return newset

######################################################## Post Processing Control Output


def add_3dsurface(control_output, phi_lim=2 * pi, overwrite=False, inplace=False):
    """
    Add 3D geometric dimensions to dataset from gpec_control_output_n#.nc.

    Note, this surface uses the machine toroidal angle. It should only by used
    with functions that do the same (tmag_out=0).

    :param control_output: Dataset. xarray Dataset opened from gpec_control_output_n#.nc
    :param phi_lim: float. Toroidal angle extended from 0 to phi_lim radians.
    :param overwrite: bool. Overwrite geometric quantities if they already exist.
    :param inplace: bool. Modify dataset inplace. Otherwise, return a new dataset.

    :returns: Dataset. New dataset with additional dimensions.

    :Examples:

    After opening, adding geometry and confirming the functional forms
    are using the machine angle,

    >>> con = data.open_dataset('examples/DIIID_ideal_example/gpec_control_output_n1.nc')
    >>> con = add_3dsurface(con)
    >>> con = add_fun(con,tmag=False)

    it is easy to make 3D surface plots using mayavi,

    >>> fig = mmlab.figure(bgcolor=(1, 1, 1), fgcolor=(0,0,0), size=(600, 600))
    >>> mesh = mmlab.mesh(con['x_surf'].values,con['y_surf'].values,con['z_surf'].values,
    ...                   scalars=np.real(con['b_n_fun'] * con['phase_surf']), colormap='RdBu', figure=fig)
    >>> mmlab.savefig(filename='examples/figures/example_bn_3d.png', figure=fig)

    .. image:: examples/figures/example_bn_3d.png
       :width: 600px

    Make 2D perturbed last-closed-flux-surface plots quickly using,

    >>> fac = 1e-1 # good for DIII-D unit eigenmodes
    >>> f,ax = plt.subplots()
    >>> ax.set_aspect('equal')
    >>> l, = ax.plot(con['R'], con['z'], color='k')
    >>> for i in range(4): # most and least stable modes
    ...     v = np.real(con['R_xe_eigenvector_fun'][i] * fac)
    ...     l, = ax.plot(con['R'] + v * con['R_n'], con['z'] + v * con['z_n'])
    >>> sm = plt.set_linearray(ax.lines[-4:], cmap='viridis')
    >>> cb = f.colorbar(sm, ticks=range(4))
    >>> cb.set_label('Mode Index')
    >>> xtl = ax.set_xticks([1, 1.7, 2.4])
    >>> txt = ax.set_xlabel('R (m)')
    >>> txt = ax.set_ylabel('z (m)')
    >>> f.savefig('examples/figures/example_boundary_wiggle.png')

    .. image:: examples/figures/example_boundary_wiggle.png
       :width: 600px

    Combining these two concepts, we can even plot a 3D perturbed surface distorted
    by the normal displacements and colored to show the normal field.

    >>> fig = mmlab.figure(bgcolor=(1, 1, 1), fgcolor=(0,0,0), size=(600, 600))
    >>> xinorm = 10 * con['phase_surf'] # good for ~cm displacements
    >>> x = np.real(con['x_surf'].values + con['xi_n_fun'] * xinorm * con['x_surf_n'])
    >>> y = np.real(con['y_surf'].values + con['xi_n_fun'] * xinorm * con['y_surf_n'])
    >>> z = np.real(con['z_surf'].values + con['xi_n_fun'] * xinorm * con['z_surf_n'])
    >>> s = np.real(con['b_n_fun'] * con['phase_surf'])
    >>> mesh = mmlab.mesh(x,y,z,scalars=s, colormap='RdBu', figure=fig)
    >>> mmlab.savefig(filename='examples/figures/example_bn_xin_3d.png', figure=fig)

    .. image:: examples/figures/example_bn_xin_3d.png
       :width: 600px

    Fancy!

    """
    if inplace:
        ds = control_output
    else:
        ds = control_output.copy(deep=True)

    # machine toroidal angle
    if 'phi_surf' not in ds or overwrite:
        phi = np.linspace(0, phi_lim, 180)
        phi = xarray.DataArray(phi, coords={'phi_surf': phi}, name='phi_surf')
        ds['phi_surf'] = phi
        ds['phase_surf'] = np.exp(-1j * ds.attrs['n'] * ds['phi_surf'])

    # 3D cartesian coords for 3D plots
    if 'x_surf' not in ds or overwrite:
        xy = (ds['R'] * np.exp(1j * ds['phi_surf'])).to_dataset(name='xy')
        ds['x_surf'] = xy.apply(np.real)['xy']
        ds['y_surf'] = xy.apply(np.imag)['xy']
        ds['z_surf'] = ds['z'] * (1 + 0 * ds['phi_surf'])

        # normal vectors
        xy = (ds['R_n'] * np.exp(1j * ds['phi_surf'])).to_dataset(name='xy')
        ds['x_surf_n'] = xy.apply(np.real)['xy']
        ds['y_surf_n'] = xy.apply(np.imag)['xy']
        ds['z_surf_n'] = ds['z_n'] * (1 + 0 * ds['phi_surf'])

    return ds


def add_fun(control_output, keys=None, tmag=False, inplace=True):
    """
    Add real space functions of theta calculated from the spectral outputs.
    This is helpful when analyzing a GPEC run for which fun_flag was off
    for speed/efficiency.

    :param control_output: Dataset. xarray Dataset opened from gpec_control_output_n#.nc
    :param keys: iterable. List of keys for which <key>_fun will be added. Default is all spectral variables.
    :param tmag: bool. Keep the magnetic coordinate angle (default is machine angle).
    :param inplace: bool. Add to dataset. Otherwise, make a new dataset object.

    :returns: Dataset. New dataset with real space functions of theta.

    :Example:

    You can check this against fun_flag outputs,

    >>> con = data.open_dataset('examples/DIIID_kinetic_example/gpec_control_output_n1.nc')
    >>> con_mag = post.add_fun(con, tmag=True, inplace=False)
    >>> con_mac = post.add_fun(con, tmag=False, inplace=False)
    >>> f,ax = plt.subplots(2)
    >>> l, = ax[0].plot(con['theta'], np.abs(con['xi_n_fun']), label='fun_flag Output')
    >>> l, = ax[0].plot(con['theta'], np.abs(con_mag['xi_n_fun']), label='Magnetic Angle')
    >>> l, = ax[0].plot(con['theta'], np.abs(con_mac['xi_n_fun']), ls='--', label='Machine Angle')
    >>> leg = ax[0].legend()
    >>> txt = ax[0].set_title('Normal Displacement')
    >>> txt = ax[0].set_ylabel('Amplitude')
    >>> l, = ax[1].plot(con['theta'], np.angle(con['xi_n_fun']), label='fun_flag Output')
    >>> l, = ax[1].plot(con['theta'], np.angle(con_mag['xi_n_fun']), label='Magnetic Angle')
    >>> l, = ax[1].plot(con['theta'], np.angle(con_mac['xi_n_fun']), ls='--', label='Machine Angle')
    >>> txt = ax[1].set_ylabel('Phase')
    >>> txt = ax[1].set_xlabel('theta')
    >>> f.savefig('examples/figures/example_add_fun.png')

    .. image:: examples/figures/example_add_fun.png
       :width: 600px

    To see the function is properly recovered. Note that the fun_flag outputs are
    on the machine coordinate.

    """
    if inplace:
        ds = control_output
    else:
        ds = control_output.copy(deep=True)

    # default to all spectral variables
    if keys is None:
        keys = [k for k, v in ds.data_vars.iteritems() if v.dims == ('m',)]

    # calculate functions
    for k in keys:
        # inverse fourier transform
        fun = (ds['xi_n'] * np.exp(1j * (ds['m'] * ds['theta'] * 2 * pi))).sum('m')
        # convert to machine toroidal angle
        if not tmag:
            fun *= np.exp(1j * ds.attrs['n'] * ds['delta_phi'])
        # add to the dataset
        ds[k + '_fun'] = fun

    return ds


######################################################## Post Processing Control Output

def optimize_torque(matrixprofile, psilow=0, psihigh=1, normalize=False, energy=False, minimize=False):
    """
    Calculate the eigenvalue and eigenvector corresponding to maximum the torque
    within the specified range of normalized poloidal flux.

    :param matrixprofile: DataArray. T_x or T_coil from gpec_profile_output_n#.nc.
    :param psilow: float. Lower limit of normalized poloidal flux.
    :param psihigh: float. Upper limit of normalized poloidal flux.
    :param normalize: bool. Optimize value within psi range normalized to the total.
    :param energy: bool. Optimize 2ndW (default is torque).
    :param minimize: bool. FInd the minimum (default is maximum).

    :returns: w,v. Eigenvalue, eigenvector of optimum.


    :Examples:

    With the output of a kinetic GPEC run in hand,

    >>> con = data.open_dataset('examples/DIIID_kinetic_example/gpec_control_output_n1.nc')
    >>> prof = data.open_dataset('examples/DIIID_kinetic_example/gpec_profile_output_n1.nc')

    This function provides the spectrum that optimizes the integral torque
    and the corresponding profile is easily obtained.

    >>> w,v = optimize_torque(prof['T_xe']) # handles any order of dims
    >>> T_mat = prof['T_xe'].transpose('psi_n', 'm_prime', 'm') # order matters in  dot products
    >>> T_opt = np.real(np.dot(np.dot(v.conj(), T_mat), v))/2 # 0.5 for quadratic fourier analysis

    We can then compare the spectrum and torque from our GPEC perturbed equilibrium to
    the theoretical optima,

    >>> f,ax = plt.subplots(2)
    >>> input_power = np.dot(con['Phi_xe'].conj(), con['Phi_xe'])
    >>> l = np.abs(con['Phi_xe']).plot(ax=ax[0], label='Input')
    >>> l, = ax[0].plot(con['m'], np.abs(v) * np.sqrt(input_power), label='Optimum Torque')
    >>> leg = ax[0].legend()
    >>> l = prof['T'].real.plot(ax=ax[1], label='Input')
    >>> l, = ax[1].plot(prof['psi_n'], T_opt*input_power, label='Optimized Torque')
    >>> leg = ax[1].legend()

    >>> f.savefig('examples/figures/example_optimize_torque.png')

    .. image:: examples/figures/example_optimize_torque.png
       :width: 600px

    """
    # enforce correct order (default)
    if 'm' in matrixprofile.dims:
        mprof = matrixprofile.transpose('psi_n', 'm', 'm_prime')
    elif 'coil_index' in matrixprofile.dims:
        mprof = matrixprofile.transpose('psi_n', 'coil_index', 'coil_index_prime')
    else:
        raise ValueError("Expected m by m_prime or coil_index by coil_index_prime matrix")

    # total torque matrix
    T1 = mprof.sel(method='nearest', psi_n=1)
    T1 = (T1 + T1.conj().T) / 2
    T1inv = np.linalg.inv(T1)

    # hermitian torque matrices on either end of the range
    Th = mprof.sel(method='nearest', psi_n=psihigh)
    Th = (Th + Th.conj().T) / 2
    Tl = mprof.sel(method='nearest', psi_n=psilow)
    Tl = (Tl + Tl.conj().T) / 2
    if psilow == 0:
        Tl *= 0  # explicitly remove lower range

    # optimal torque within range
    Topt = Th - Tl
    # optimal torque normalized to total torque
    if normalize:
        Topt = np.dot(T1inv, Topt)

    # calculate eigenvalues
    w, v = np.linalg.eig(Topt)

    # find relevant "optimum"
    argopt = np.argmax
    if minimize:
        argopt = np.argmin
    if energy:
        i = argopt(np.imag(w))
    else:
        i = argopt(np.real(w))

    return w[i], v[:, i]
