.. _releases:

This is an active, developing scientific research code. The release structure attempt to follow the semantic versioning guidelines. Given a version number MAJOR.MINOR.PATCH,

- MAJOR version is incremented when incompatible API changes are made,
- MINOR version is incremented when functionality is added in a backwards-compatible manner, and
- PATCH version is incremented when backwards-compatible bug fixes are made.

Note, backwards compatibility is defined as the ability to return to a previous run and re-run any previously successful executable. This does not guarantee all output file formats are exactly maintained.

The `github release notes <https://github.com/PrincetonUniversity/GPEC/releases>`_ are reproduced below.

GPEC v1.4.3
===========

Fixes
------
- VACUUM - Respects the namelist a and lspark when determining whether to write a wall geometry file output
- DCON - Properly spaces the vacuum matrices used when finding the peak dW with psiedge < psilim

Changes
------
- DCON - dW_edge in the netcdf output now includes the imaginary component of the least stable eigenvalue


GPEC v1.4.2
===========

Fixes
------
- PENTRC - Respects the namelist verbose setting throughout the program

Changes
------
- COIL - Updates ASDEX Upgrade coil geometry to have more accurate, rounded cornering


GPEC v1.4.1
===========

Fixes
------
- GPEC - Corrects helicity for brzphi calculations in machine angle
- PENTRC - Fixes reading of xclebsch displacements from other codes that use psi as the outer loop


GPEC v1.4.0
===========

Fixes
------
- INSTALL - Support added for compiling on MacOS as well as with gfortran on linux
- VACUUM - Makes ishape 41 use on manual wall robust to a>10 infinite wall tests
- DCON - Corrects for wall time wrap-around when run spans multiple days
- GPEC - Removes (incorrect) units from normalized minor radius netcdf output

Changes
--------
- DCON - Improves parallel assembly of kinetic matrices using collapsed psi and ell loops and dynamic scheduling. This enables efficient use of large parallel_threads.

Adds
------
- EQUIL - Adds new pow1 and pow2 grid_type options. These grids are linear in core psi_n and approach zero spacing at psihigh linearly/quadratically (packed in the edge).
- DCON - Adds new psiedge variable to dcon_control. When this is below psihigh, DCON records dW(psi) between psiedge and psihigh and re-runs the integration with truncation adjusted to max(dW).
- PENTRC - Adds new option (force_xialpha) that triggers a calculation of tangential from radial displacement (useful for radial displacement inputs from nonlinear codes like JOREK and M3DC1).


GPEC v1.3.8
===========

Fixes
------
- GPEC - Fixes bug that caused GPEC netcdf to crash when there are no rational surfaces (common when running DCON with kin_flag)
- DCON - Defaults the reform_eq_with_psilim flag to false for better consistency with earlier versions


GPEC v1.3.7
===========

Fixes
------
- PENTRC - Corrects sqdivxi_perp_SA output, which was being overwritten with sqdBoB_L_SA


GPEC v1.3.6
===========

Fixes
------
- DCON - Fixes openmp reduction of kinetic matrices when kin_flag=t

GPEC v1.3.5
===========

Fixes
------
- INSTALL - Minor makefile robustness improvements

GPEC v1.3.4
===========

Fixes
------
- STRIDE - Fixes faulty writing of 2Nx2N Delta matrix into a NxN Delta_prime in netcdf

Adds
------
- STRIDE - Writes full 2Nx2N Delta matrix as well as NxN A', B', Gamma', & Delta' sub-matrices to netcdf

GPEC v1.3.3
===========

Fixes
------
- DCON - Re-enables inverse equilibria, which raised errors in the more rigorous bicube checking introduced in 1.3.2

GPEC v1.3.2
===========

Adds
------
- INSTALL - Generalizes makefiles for better portability (intel/pgi/gfortran)
- INSTALL - Enables parallel or serial compilation (controlled by OMPFLAG)
- EXAMPLES - Adds a new example showing how to run rdcon and rmatch to determine tearing instability
- EQUIL - Adds spline root finding subroutines
- DCON - Adds the reform_eq_with_psilim flag to avoid reforming the equilibrium splines to psilim
- GPEC - Adds explicit signfld overlap quantities to netcdf, identical to ascii values
- PENTRC - Adds new dynamic_grid flag to enable skipping dynamic integration if only a static grid output is desired
- PENTRC - Parallelizes pentrc bounce harmonic integrations

Changes
--------
- INSTALL - Updates instructions for new AUG cluster
- GPEC - Moves overlap outputs to control netcdf, consistent with correct dimensions (mode, not psi_rational)
- GPEC - Renames netcdf dimensions *_rat to *_rational for clarity

Fixes
------
- DCON - Makes dcon robust to peak_flag when truncation is already very near the seperatrix
- GPEC - Fixes calculation of vacuum fields when using mode_flag
- GPEC - Fixes dimension of local coupling matrices in netcdf
- GPEC - Fixes netcdf output of full complex vsbrzphi
- GPEC - Enables GPEC when DCON used peak_flag
- PENTRC - Makes pitch integration spline evaluations external to avoid conflicts that caused LSODE to fail

GPEC v1.3.1
===========

Fixes
------
- INSTALL - Fixes a makefile error

GPEC v1.3.0
===========

Adds
------
- STRIDE - Adds the State Transition Rapid Integration with DCON (Asymptotic) Expansions code
  * Citations included in docs (website)
  * Parallel calculations of plasma stability (delta-W and Delta-prime)
- DCON - Adds new peak_flag that checks the dW every step after the last rational and ends integration at the first local maximum
  * This provides a useful for consistent, physics based truncation choice that can be more robust to small EFIT changes

Changes
--------
- DOCS - Improves online documentation
- INSTALL - Updates flags for intel 2018 and adds instructions for IPP Max-Planck Garching

Fixes
------
- DCON - Fixes a formatting error in sing_find.out
- DOCS - Fixes rmatch eta and massdens inputs for DIIID_resisitive_example
- DCON - Fixes inappropriate uses of psihigh, which may not be the end of integration psilim if sas_flag, qhigh, or peak_flag are used

GPEC v1.2.3
===========

Changes
-------
- DOCS - Modernizes online documentation

Fixes
------
- GPEC - Fixes incorrect w_isl unit labeling in netcdf and makes it the full width for consistency with w_isl_v


GPEC v1.2.2
===========

Adds
------
- GPEC - Adds new coils for DIII-D, COMPASS, ASDEX Upgrade, and NSTX-U

Changes
--------
- DOCS - Updates documentation of authors, public installations, and compilation at each institution
- DCON - Makes formatting of surface quantities table in dcon.out consistent between DCON and RDCON
- DCON - Makes sum1.bin (previously sum1.dat) outputs consistent with CALTRANS DCON.
- GPEC - Singular coupling routines are skipped entirely (avoiding possible issues) if no rationals exist in the domain

Fixes
------
- GPEC - Fixes a ss_flag bug so the penetrated vacuum field is calculated on the correct surfaces


GPEC v1.2.1
===========

Fixes
------
- DOCS - Fixes a mismatch in the jacobian's between GPEC and PENTRC (the torques profiles now match)


GPEC v1.2.0
===========

Adds
------
- RDCON & RMATCH - Adds resistive DCON packages
  * Calculates inner layer model and performs matching with ideal outer layer
  * Cite [A.H. Glasser, Z.R. Wang, and J.-K. Park, Physics of Plasmas 23, 112506 (2016)]
  * Includes new resistive examples
- DCON - Adds ability to calculate bounce harmonics in parallel when forming kinetic matrices for Euler-Lagrange Equation
- DCON - Adds ability to start the Euler-Lagrange ODE integration at or above a minimum safety factor qlow
- DCON - Adds ability to include electron kinetic terms in Euler-Lagrange equation (controlled by the new electron_flag)
- VACUUM - Adds ability to read prescribed wall geometry file
- GPEC - Adds new singfld and signcoup calculations and includes all singcoup and singfld outputs in netcdf
  * Delta: is the published resonant coupling metric [Park, Phys. Plasmas 2007] normalized by B.grad(theta) instead of B.grad(phi).
  + This is similar to what is sometimes called external Delta' in tearing stability theory **need to divide by vacuum**
  * Penetrated resonant flux: interpolated across singspot, and physically meaningful for kinetic MHD equilibria only
- DCON, GPEC, & PENTRC - Updates the version based on the compile-time git commit
- DCON, GPEC & PENTRC - Can use a classical spline coefficient solution for "extrap" boundary condition splines, avoiding  a suspected (minor) bug in the original tri-diagonal solution that resulted in large grad-shafranov errors in poorer quality equilibrium (especially inverse or modified equilibrium).
  * Previous tri-diagonal spline solutions can be recovered by setting use_classic_splines to false

Changes
--------
- DCON - Improves clarity of singular surface search messages
- GPEC - Improves clarity and consistency of singular coupling outputs
  * Uses iszinv to invert hermitian fldflxmat
  * Uses area normalization of penetrated flux for consistency with effective flux
  * Adds unique names for the singcoup mat and svd ascii outputs (enables python reading)
- PYPEC - Improves automatic selection of partitions and threads in job submission and adds rdcon to exe options

Fixes
------
- DCON - Fixes only the the plasma energy matrix written to dcon.out to include full matrix (previously only 2 columns)
- GPEC - Improves clarity and consistency of singular coupling outputs
  * Corrects units of Phi_res in netcdf (area normalized, so T not Wb)
  * Corrects units and calculation of island width in netcdf (unitless width in psi_n, required a sqrt)
- GPEC - Fixes bug in iszinv for m/=mpert matrices (no impact on previous results, which all used m=mpert)
- GPEC - Fixes bug in the normalization of singular coupling islandwidths (singdfld unchanged)
- GPEC - Fixes poor formatting in response file header
- PENTRC - Corrects the sign of the charge when calculating NTV torque and kinetic delta-W for electrons
- VACUUM - Makes vacuum code robust to namelists without a header line

Removes
--------
- ALL - Removes official support for all compilers other than intel
  * Parallel openmpi calls unique to intel
  * Move is consistent with RDCON development path


GPEC v1.1.7
===========

Features
---------
- DCON - A new, explicit ion flag toggles whether the ion kinetic energy is included in the kinetic Euler-Lagrange equation


GPEC v1.1.6
===========

This release corrects a bug that may have made previous GPEC electron NTV have the incorrect sign.

Fixes
----------
- PENTRC - Corrected the sign of the charge (diamagnetic frequencies, etc) for electron calculations.


GPEC v1.1.5
===========

This version includes a minor but important change to make the ideal GPEC eigenfunctions almost identical to those from DCON in IPEC. A power extraction essential for numerical stability when forming the fundamental H and G matrices in the kinetic solutions has been removed from the ideal calculations for consistency with the previous calculations in the ideal case.

Adds
---------
- COIL - New coils are available for JET, NSTX, and COMPASS. The number of coils usable in a run increased.
- GPEC - The q, rho, and volume profiles are included in the netcdf output if any profile output is requested.
- GPEC - The local coupling matrix between opsi1 and opsi2 and corresponding svd vectors are available. **needs netcdf output??**

Fixes
----------
- DCON - Fundamental matrices only use power extraction technique when kin_flag is true.
- PENTRC - Progressbars are now called at the end of do loops for more precise reporting.
- PENTRC - Torque estimation from surface currents is now recorded in harvest and netcdf.

Documentation
--------------
- EXAMPLES - Examples now include "run" examples with J.-K. Park's typical workflow and settings.
- INPUT - Annotations and settings of default input namelists include minor changes.
- PYPEC - Mayavi instructions are updated for latest portal python installations.


GPEC v1.1.4
===========

Fixes
----------
- COIL - Fixed faulty 1.1.3 implementation of increasing the east coil windings.


GPEC v1.1.3
===========

Fixes
----------
- COIL - Increased the number of windings for the up and down EAST coil arrays


GPEC v1.1.2
===========

Fixes
--------------
- PENTRC - Now successfully writes kinetic profiles on the equilibrium grid to netcdf files


GPEC v1.1.1
===========

Fixes
------------
- PYPEC - A bug was fixed in the python processing tools' optimize_torque function


GPEC v1.1.0
===========

This release includes a new DCON netcdf output file and SLURM job submission interface in PYPEC for compatibility with the new portal and iris computing standards. Details are below.

Adds
---------
- DCON - A clean, efficient netcdf file replicates the information in the complicated dcon.out ascii.
- DCON - The new namelist variable, out_fund, toggles fundamental matrix output (ABCDEH in imats.out fs.bin, ks.bin and gs.bin).
- COIL - KSTAR and EAST coils are available.
- COIL - A NSTX-U error field model is available.
- GPEC - Control netcdf outputs include the external flux applied from each coil and coil names.
- GPEC - Profile netcdf outputs include rational surface quantities, coil names, and vsbrzphi, xbrzphifun, and arzphifun outputs.
- GPEC - Code is robust to singfld_flag with con_flag.
- GPEC - The new namelist variables, ascii_flag and netcdf_flag, toggle all ascii and netcdf outputs respectively.
- PYPEC - SLURM job submission.
- PYPEC - Post processing includes a function that updates netcdf naming conventions to be consistent with the latest version.
- PYPEC - Backwards compatibility for running ipec is available.
- REGRESSION - Tools for comparing versions are available.

Fixes
----------
- DCON, GPEC, PENTRC - Timers were fixed to correctly handle multi-day runs.
- DCON - Ascii formatting is updated for complex eigenvalue energies.
- GPEC - An indexing offset in calculation in dw_flag torque matrix output was fixed.
- GPEC - Appropriate ascii closing was added.

Documentation
--------------
- DOCS - Documentation includes compare module.
- INPUT - Annotations and settings of default input namelists include minor changes.


GPEC v1.0.6
===========

This patch features fixes to a number of deeply embedded indexing and memory allocation bugs. This is necessary for compiler robustness. The regression examples show essentially no change in the results to machine precision on portal.

Fixes
----------
- VACUUM & LSODE - This patch fixes the misallocation of memory for input arrays in a number of old subroutines.
- EQUIL - This patch fixes the misallocation of memory for temporary arrays in Fourier spline fitting.
- GPEC - This patch fixes an index offset in the matrices forming the torque matrix profile.


GPEC v1.0.5
===========

Fixes
-----------
- Fixed normalization of filter_flag energy normalized field decomposition.

This bug was introduced with the new normalized field (T) convention in 1.0.2. To correct the decomposed energy normalized flux O_*Phi_xe in versions 1.0.2-1.0.4, multiply by 1/sqrt(A).


GPEC v1.0.4
===========

Avoids repetition of dimensions in control netcdf J_surf_2.
Note this is not critical for the netcdf, but necessary for the way pypec and xarray treat dimensions.

GPEC v1.0.3
===========

This patch fixes a mis-labeling of the control netcdf Phi_fun and Phi_x_fun units. The units are Wb.


GPEC v1.0.2
===========

This patch features one bug fix and one addition to the netcdf output.

Adds
--------------
- A transform matrix J_surf_2 has been added to the control netcdf. This matrix applies a dimensionless half-area weighting.

Fixes
-------------
- The netcdf output Phi_xe has been changed from "energy-normalized flux" with units Wb/m to "energy-normalized field" with units of Tesla. The related \*_xe matrices have been similarly normalized. No physics is changed, only the scalar area normalization.


GPEC v1.0.1
===========

This patch cleans up the input directory, removing deprecated files.


GPEC v1.0.0
===========

This major release marks the true transition from individual ideal perturbed equilibrium calculations to a fully generalized perturbed equilibrium package.

The Perturbed Equilibrium Nonambipolar TRansport Code (PENTRC) is used to calculate the neoclassical drift kinetic pressure matrixes required to minimize the hybrid kinetic-MHD perturbed energy and find a set of force balance states. The computational structure of the ideal DCON code is largely maintained in finding these states, although generalizations and modifications have been made to account for new mathematical properties. Foremost among these are 1) the absence of hermitian properties and 2) the integrable nature of singularities near the rational surfaces. Generalization of the linear algebra and new decomposition / recomposition of the matrices required by these changes are now used for both the ideal and kinetic calculations.

The Ideal Perturbed Equilibrium Code (IPEC) has officially been deprecated and is now the package namesake: the Generalized Perturbed Equilibrium Code (GPEC). The foundational computational changes are much less than in the above case however, with only a few minor generalizations of hermitian linear algebra assumptions.

Adds
-------------
 - DCON inclusion of kinetic terms is now determined by the kin_flag input.

    - Additional dcon_control namelist inputs can be used to control the kinetic calculations

 - IPEC now calculates generalized perturbed equilibrium (no assumption that the force balance states form a hermitian matrix)
 - IPEC netcdf output is nearly complete and naming conventions are official
 - PENTRC now has fully netcdf output unless ascii is specifically requested by the user

    - Output is now separated from calculations, setting the stage for parallelization

Documentation
----------------------
 - Example runs have been split into ideal and kinetic examples to show the kinetic effects
 - An "a10" example has been added for simple circular-large-aspect-ratio intuition


GPEC v0.4.0
===========

This release includes a number of minor I/O changes and convenient default input features as well as a few minor bug fixes.

Fixes
--------------

- MATCH updated interface for changes DCON file formats
- IPEC fixed alignment of columns in xclebsch_fun output

Features
-------------

- DCON, IPEC, PENTRC all accept the additional Jacobian type 'park'

  + Sets the power of (b,bp,r) to (1,0,0)

- IPEC includes (r,z) in xclebsch_fun output
- COIL, PENTRC the data_dir used to look up hardcoded data now accepts defaults to $GPECHOME/pentrc

  + This option is used when set to 'default' or ''

- PENTRC now includes a valid circular large-aspect-ratio calculation

  + Calculates Eq. (19) from [Logan, Phys. Plasmas, 2013] using Eqs. (10-12) from [Park, Phys. Rev. Lett. 2009] with the kappa dependence
  + Previous versions included this flag as a placeholder only and should not be used

Documentation
----------------------

- Example namelists updated to use native coordinates throughout for increased speed and clarity


GPEC v0.3.5
===========

This release includes critical bug fixes for the nonambipolar transport calculations in PENTRC.

Fixes
--------------

- PENTRC a correction factor of 1/2 has been applied to the fcgl, \*gar, and \*mm methods to correctly represent quadratic terms using complex analysis
- PENTRC xclebsch is now correctly transformed back to DCON working coordinates when output on more m than the DCON mpert.


GPEC v0.3.4
===========

This release includes a number of critical bug fixes found and fixed in a general review of the ideal MHD package in preparation of the move to kinetic MHD version 0.4.0 under development. It also includes a few (re-)standardizations of features.

Fixes
--------------

- PENTRC +/- omega_b included for passing and not trapped particles, removing unphysical symmetry in ell of trapped particle torques

  + **All previous 0.3 version torques should be considered incorrect**

- PENTRC fixed bug in inverse Fourier transformation of perturbed quantities and fixed (removed) JBB normalization of perturbed quantity splines for consistent treatment in GAR, LAR, and GCL methods (now benchmarked with PENT).

  + **All previous 0.3 version LAR and CGL torques should be considered incorrect**

- PENTRC returned factor of 2 to all GAR methods (now benchmarked against PENT for MDC2 cases)
- PENTRC fixed radial grid outputs from (over)writing sum and individual ell profiles to same file
- PENTRC enforce psi limits on grid outputs
- IPEC fixed bug in writing O_CX, b_nm, b_xnm, xi_nm, and xi_xnm to control netcdf file

  + **All previous 0.3 version values should be considered incorrect**

- IPEC working jacobian power factors are explicitly enforced when jac_in or jac_out re not specified
- IPEC fixed bug using wrong jacobian and angle in ipeq_fcoordsout conversions (not used in any previous version)
- IPEC ipeq_fcoordsout and ipeq_bcoordsout always perform transformation on larger of the working/output m grids (not expected to be an issue for previous versions)

Features
-------------

- IPEC output coordinate m range is now determined by a new IPEC_OUTPUT variable mlim_out
- IPEC the control surface theta-space function values are now always calculated and output
- IPEC bwp_pest_flag is now true by default and produces pest ouputs for both xbnormal and vbnormal
- IPEC xclebsch outputs are now converted to output coordinates and theta-space outputs are available
- PENTRC now accepts jsurf_in, tmag_in and all individual powers of the jac_in, allowing it to interface with IPEC's new xclebsch outputs that are transformed from the working to ipec output coordinates

  + Coordinate transformation back to the DCON working coordinates is done on the large of the working/input m grids

- IPEC added helicity to control and profile netcdf outputs
- PENTRC now has the option to override the perturbed quantities calculated using the xclebsch interface with a direct ipec_pmodb ascii interface (when the user specifies a pmodb_file)
- PENTRC now enforces that a substring of the form 'n#' where # is the DCON toroidal mode number be in the peq_file file name

Speed and Stability
---------------------------

- PENTRC only runs the psi_out surfaces if detailed outputs are actually requested
- PENTRC exclude trapped/passing boundary from pitch-space splines using power-grids approaching from either side
- INSTALL and all individual makefiles have updated from the develop branch, reorganizing the linking order and allowing diverse machine/compiler options.

Documentation
----------------------

- Updated input and example namelists and their annotation


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

  + If the jac_out was not the working jac_type Phi^x outputs in the jac_out table were mistakenly in the jac_in coordinate system.

- The external and total flux have been added to the control netcdf alongside their previously stored energy normalized values.


GPEC v0.3.1
===========

Fixes
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

Fixes
------------

- DCON qhigh is enforced independent of sas_flag
- IPEC longstanding bug that caused crashes when brzphi was requested without eqbrzphi is fixed
- IPEC mthsurf bug fixed

  + Benchmarks show perfect recovery of excessively high DCON mthsurf results using mthsurf=1

Features
-------------

- COIL now includes 4x48 RFX-mod coils
- IPEC netCDF output is now available for major output flags (more will be transitioned soon)

  + Currently netCDF files include: filter_flag, x/brzphi_flag, xbnormal_flag, pmodb_flag, and control surface fun_flag outputs

- IPEC output subroutines can now be individually timed using the timeit flag
- IPEC mode filtering has a new filter_type, filter_modes interface in IPEC_INPUT
- IPEC reduced terminal printing - no longer is every eigenmode printed to the terminal

Performance
------------------

- IPEC speed was increased by saving coordinate transformation information on a surface when performing multiple transformations on one surface
- IPEC brzphi speed was increase by 1-2 orders of magnitude by calculating (r,z,phi) quantities on the requested grid points instead of across surfaces
- IPEC speed can now be confidently increased by a large factor using the mthsurf flag (see bug fix)
- IPEC compiler optimizability increased with the switch from pointers to allocatable arrays
