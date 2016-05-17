******************
Accessing Releases
******************


Official releases are available for download from the `github release page <https://github.com/PrincetonUniversity/GPEC/releases>`_.

On the PPPL portal computers, releases are available at /p/gpec. To use a release,::

    module use /p/gpec/modules

will give you access to the official gpec modules. Use,::

    module avail

to see the available releases, or simply,::

    module load gpec

to load the current default. This will add the appropriate directories to your path, such that you can simply
call the executables (i.e. dcon, etc.).

On the General Atomics iris computer, releases are available at /fusion/projects/codes/gpec. To use a release, simply load the publicly available gpec module.


******************
Release Notes
******************

The `github release notes <https://github.com/PrincetonUniversity/GPEC/releases>`_ are reproduced below.


GPEC v0.3.3
===========

This release features a critical bug fix for control surface netcdf output and pmodb/xbnormal outputs

 - All area normalized or energy normalized quantities were incorrectly converted to the users specified jac_out coordinates. All quantities are now in the DCON jac_type coordinate system unless specifically noted otherwise.
 - The jacobian and surface area have been added to the control netcdf as global attributes
 - The filtering of singular coupling modes is now done entirely within the DCON coordinate system, for which a new singular coupling matrix is formed and SVD'd.

 - Bugs in the use of bcoordsout for pmodb and xbnormal profile quantities that wrote the first variable to multiple variables (i.e. eulb to lagb) were fixed.
 - A Bug in the weighting of the bwp profile was fixed


GPEC v0.3.2
===========

- This release features a critical bug fix for control surface ascii output Phi^x.
  - If the jac_out was not the working jac_type Phi^x outputs in the jac_out table were mistakenly in the jac_in coordinate system.
- The external and total flux have been added to the control netcdf alongside their previously stored energy normalized values.


GPEC v0.3.1
===========

Bugfixes
------------
 - IPEC fixed mistaken use of Hermitian lapack subroutines for permeability matrix
 - PYPEC synthetics properly closes synthetic surfaces that cover the full poloidal angle (vessel wall, etc.)
 - PYPEC coil plotting bug fixes for axes and color key words
 - PYPEC updated to reflect move from xray to xarray

Features
-------------
 - IPEC netcdf additions, including control surface matrices, profile quantities, shot/time/machine, and more
 - IPEC netcdf names conform to netcdf conventions
 - IPEC all netcdf outputs converted to jac_out
 - IPEC filter decomposition modes are now all in ascending order (SVD convention)
 - IPEC added amplification to filter modes
 - COIL added MAST coils
 - PENTRC added new grid options, which now include equil_grid or input_grid (i.e. the DCON grid)
 - PYPEC improved ascii/netcdf interface using data.open_dataset
 - PYPEC synthetics now includes magnetic sensors
 - PYPEC add_control_geometry function expands control surface geometry for 2D and 3D plots
 - PYPEC improved colormaps and automatic colormap choices
 - PYPEC now uses seaborn for context/palettes, has custom set_context function
 - PYPEC custom subplots automatically re-size figure to keep axes size
 - PYPEC now has png_to_gif function for making movies

Performance
------------------
 - Improved speed of ipeq_bcoordsout/ipeq_fcoordsout by checking for unnecessary calls to ipeq_fcoords/ipeq_bcoords


GPEC v0.3.0
===========

Bugfixes
------------
 - DCON qhigh is enforced independent of sas_flag
 - IPEC longstanding bug that caused crashes when brzphi was requested without eqbrzphi is fixed
 - IPEC mthsurf bug fixed
   - Benchmarks show perfect recovery of excessively high DCON mthsurf results using mthsurf=1

Features
-------------
 - COIL now includes 4x48 RFX-mod coils
 - IPEC netCDF output is now available for major output flags (more will be transitioned soon)
   - Currently netCDF files include: filter_flag, x/brzphi_flag, xbnormal_flag, pmodb_flag, and control surface fun_flag outputs
 - IPEC output subroutines can now be individually timed using the timeit flag
 - IPEC mode filtering has a new filter_type, filter_modes interface in IPEC_INPUT
 - IPEC reduced terminal printing - no longer is every eigenmode printed to the terminal

Performance
------------------
 - IPEC speed was increased by saving coordinate transformation information on a surface when performing multiple transformations on one surface
 - IPEC brzphi speed was increase by 1-2 orders of magnitude by calculating (r,z,phi) quantities on the requested grid points instead of across surfaces
 - IPEC speed can now be confidently increased by a large factor using the mthsurf flag (see bug fix)
 - IPEC compiler optimizability increased with the switch from pointers to allocatable arrays