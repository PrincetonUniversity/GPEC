.. _releases:

This is an active, developing scientific research code. The release structure attempt to follow the semantic versioning guidelines. Given a version number MAJOR.MINOR.PATCH,

- MAJOR version is incremented when incompatible API changes are made,
- MINOR version is incremented when functionality is added in a backwards-compatible manner, and
- PATCH version is incremented when backwards-compatible bug fixes are made.

Note, backwards compatibility is defined as the ability to return to a previous run and re-run any previously successful executable. This does not guarantee all output file formats are exactly maintained.

The `github release notes <https://github.com/PrincetonUniversity/GPEC/releases>`_ are reproduced below.

GPEC v1.6.0
===========

This is planned to be the last minor release of the Fortran GPEC 1.x series.
New feature development continues in the Julia GPEC 2.x code; the 1.6.x line will receive bug-fix patch releases only.
GPEC is now distributed under the MIT license.

Read the Upgrade notes first: the build system, several binary interface files, and a few namelist and output names changed, and many numerical results change because long-standing bugs were fixed (see Fixes).

Upgrade notes
-------------
- INSTALL - The build system was rewritten. Build only from ``install/`` with ``make`` (``make -j N`` is supported); the per-directory makefiles are gone. ``FC`` and ``CC`` must be set, git is required at build time, and missing netCDF/HDF5/OpenBLAS are built from bundled git submodules (``git submodule update --init``). ``DEPSINSTALLDIR`` installs them to a shared location; ``make realclean`` is required before switching compilers.
- ALL - Rerun the full chain after upgrading. DCON ``euler.bin`` has a new record (``wv_farwall_flag``), and the resistive interface files ``delta.bin``, ``gal_solution.bin``, ``globalsol.bin`` and ``vmat.bin`` changed format, so RDCON, RMATCH and GPEC must be rerun together.
- DCON, RDCON - ``vacuum.bin`` is no longer written by default (set ``bin_vac=t``). GPEC always recomputes the vacuum matrices and ignores ``ivacuumfile``. ``ahg2msc.out`` is no longer written; vacuum data are passed in memory (``out_ahg2msc`` in ``equil.in`` restores the files).
- RDCON - ``coil%rpec_flag`` moved from ``&GAL_OUTPUT`` to ``&RDCON_CONTROL`` in ``rdcon.in``. Setting it in ``&GAL_OUTPUT`` still works but prints a deprecation warning.
- RMATCH - ``ntor`` is no longer an input (n is read from ``delta.bin``); remove it from any ``rmatch.in`` that sets it. The default inner-layer solver is now ``inps_type="inps"`` with ``nx=512``, and the eigenvalue search uses the ``2*pi*n*Omega`` rotation convention of RPEC.
- DCON, PENTRC, COIL, STRIDE - Thread counts default to ``OMP_NUM_THREADS`` (inputs default to 0). ``parallel_threads`` is renamed ``dcon_kin_threads`` and ``openmp_threads`` is renamed ``pentrc_threads``; the old names still work with a deprecation warning.
- EQUIL - ``newq0 < 0`` (including the undocumented ``newq0=-1`` sign flip) now stops with an error. ``eq_type`` ``marklin_direct``, ``marklin_inverse`` and ``hansen_inverse`` were removed; ``chease2`` is renamed ``chease_ascii`` (the old name still works with a warning).
- GPEC - Renamed outputs: netCDF ``w_isl_crit`` is now ``w_isl_v_crit``; cylindrical ``vcbr``, ``vcbz``, ``vcbp`` are now ``b_r_vacuum``, ``b_z_vacuum``, ``b_p_vacuum``. The ``gpec_singfld_n#.out`` ASCII columns changed (new threshold, island-width and pressure columns). ``I_overlap_norm`` and the DCON/RDCON ``W_v_eigenvalue`` are now 1-D real arrays.
- RDCON - The netCDF output is ``rdcon_output_n#.nc`` (it no longer overwrites DCON's ``dcon_output_n#.nc``).
- PENTRC - With ``pmodb_file``, PENTRC now reads the supplied ``Bdivxi_perp`` column. Regenerate ``gpec_pmodb_n#.out`` with GPEC 1.6 before using it.
- DOCS - All DIII-D examples now use a public TokaMaker DIII-D-like H-mode equilibrium and kinetic profiles in place of the KEFIT 147131 files.

Adds
------
- SLAYER - New slab-layer linear drift-MHD module (Park, PoP 2022) with a standalone ``slayer`` executable, ``input/slayer.in`` and an a10 example.
- SLAYER - Tearing-mode growth rates at the STRIDE rational surfaces by matching the inner-layer and outer Delta' (optional critical Delta, ``dc_type``), per surface or coupled, with the E x B Doppler shift (``slayer_output_n#.nc``; see the DIII-D ideal example).
- GPEC - SLAYER resonant-field penetration threshold at each rational surface (``singthresh_slayer_flag``; ``singthresh_callen_flag`` selects the Callen threshold and ``singthresh_flag`` enables both). New netCDF ``Phi_res_crit``.
- GPEC - Prandtl number for the SLAYER threshold set by ``singthresh_slayer_inpr`` (default 5.0) or per rational surface by ``singthresh_slayer_inpr_prof``.
- GPEC - Callen saturated and minimum island widths (``w_isl_sat``, ``w_isl_min``) and the Callen critical resonant field ``Phi_res_crit_callen``.
- GPEC - Estimate of the pressure lost to Callen-saturated islands (``dP_res``, ``P_res`` and the modified pressure profile ``p_mod``).
- GPEC - Bernstein/Chance decomposition of the eigenmode potential energy (``effective_field_flag``, ``standard_form_flag``, ``effective_field_terms_flag``, ``curl_check_flag``, ``recon_out``, ``recon_int``) with ``gpec_recon*`` outputs and a derivation note in ``docs/tex/recon``.
- GPEC - Mutual inductance between the plasma surface and a conformal wall and its inverse (``M_w``, ``M_w_inv``); ``wv_farwall_flag`` in ``dcon.in`` keeps no-wall eigenmodes for GPEC while the wall inductance is computed; ``mutual_test_flag`` diagnostics.
- GPEC - Resistive perturbed equilibria (RPEC) carry the analytic Galerkin derivative d(xi)/dpsi and the RMATCH field through to GPEC. Per-surface resistive sweet spots from the Lundquist number (``use_res_spot``, ``sing_resspot``, ``sing_interpspot``) and new outputs ``sing_spots``, ``Sfac``, ``tauA``, ``tauR``, ``Jbgradpsi_rmatch`` and ``b_pen``.
- GPEC - Dominant-mode overlap outputs ``dominant_mode_overlap`` and ``dominant_mode_local_overlap``.
- GPEC - Full set of DCON plasma attributes and local SLAYER inputs in the profile and control netCDF outputs.
- COIL - Shift and tilt coils from ``coil.in`` (``coil_shift[xyz]``, ``coil_tilt[xyz]``, ``coil_[xyz]nom``, ``cnpert``), with optional netCDF output of the modified geometry (``coil_out``).
- COIL - OpenMP-parallel Biot-Savart fields with a new ``coil_threads`` input; the maximum number of coils is raised from 50 to 500; MAST-U, JET ex-vessel and extra ITER coil sets.
- RDCON - netCDF output ``rdcon_output_n#.nc`` (``netcdf_out``, default true) with the ``dcon.out`` content, the Delta, Delta', Delta_gw and Delta_coil matrices and the Galerkin cell layout.
- RDCON - Optional modified Rutherford equation terms (``MRE_flag``, ``geom_flag``) with a radial Zeff profile (``Zeff``, ``psi_N_Zeff``).
- RDCON - STRIDE-style truncation controls ``qlow``, ``qhigh`` and ``reform_eq_with_psilim``; ``gnstep`` sets the maximum LSODE steps for the Galerkin basis; ``interp_np_res`` refines output in resonant cells.
- MATCH - ``sol_num`` selects which DCON eigenfunction MATCH reconstructs.
- PENTRC, DCON - ``kinetic_format="rates-v1"`` reads a 9-column kinetic file that supplies logLambda, nu_i and nu_e directly (default ``"legacy"`` unchanged).
- PENTRC - Separate energy-integration tolerances ``atol_x``, ``rtol_x`` (tighter by default, so runs take longer); Zeff and ``omega_p`` in the netCDF profile output.
- EQUIL - ``profile_source="integrate"`` builds F and p from the g-file's FF' and p' (recommended only for single-precision g-files such as TokaMaker's); default ``"values"`` is unchanged.
- EQUIL - ``ldp_i`` (TokaMaker i-file) reader accepts F < 0 and an optional ``[psi_axis, psi_lcfs]`` record.
- EQUIL - ``enstep`` sets the maximum field-line integration steps in ``direct_run`` (default 16384, previously fixed at 2048).
- EQUIL - LAR ``sigma_type='wesson'`` current profile; ``chease_ascii`` eq_type.
- VACUUM - New Green's function evaluation (Chance 2007), more accurate at high n; ``use_legacy_greens_function`` restores the old one.
- INSTALL - gfortran (with OpenMP), ifx and ARM macOS support; ``OPENBLASHOME``; ``make deps`` / ``make checkdeps``; build instructions for NERSC, SNU 3DST, PLARE and WSL.
- CI - GitHub Actions builds with gfortran and ifort and runs Solov'ev ideal, kinetic, resistive and energy-decomposition examples on every pull request.
- DOCS - Installation page, sign-conventions reference, resonant coupling, energy-decomposition and pmodb gauge output docs, and the ``docs/tex`` notes (moved from ``tex/``).

Changes
--------
- ALL - GPEC is released under the MIT license (``install/LICENSE``), replacing the previous user agreement.
- ALL - Floating-point kinds are standardized on ``REAL(r8)``; LSODE updated to Fortran 2018; duplicate utility, grid, spline and jacobi code consolidated; obsolete files removed.
- ALL - Consistent ``!! WARNING`` terminal messages; GPEC reports CPU and wall time separately for OpenMP runs.
- DCON - Stops with a clear error when ``nn = 0``.
- GPEC - Singular-coupling overlap outputs have correct units (Tesla, or % for ``*_overlap_norm``) and consistent descriptions.
- GPEC - ``singthresh_callen_flag`` forces ``vsingfld_flag`` when coils are present; ``gal_flag`` runs always use the resistive sweet spot.
- RDCON - Free-boundary energies use the DCON vacuum-call pattern and a complex eigen-solve; ``euler.bin`` stores the complex eigenvalues.
- RDCON - ``solver="cholesky"`` switches to LU with a warning when ``rpec_flag=t``.
- RDCON, STRIDE - Modules renamed with ``rdcon_``/``stride_`` prefixes; STRIDE builds without OpenMP and uses ``OMP_SET_MAX_ACTIVE_LEVELS``.
- STRIDE - Banded solve of the Delta' boundary-value problem (faster at high n).
- RDCON - Lower peak memory in the LU Galerkin solve.
- PENTRC - Exact-zero ``omega_E`` knots are replaced by 1e-9 rad/s with a warning.
- PENTRC, STRIDE - LSODE failures in the PENTRC energy and radial integrations and ZVODE failures in STRIDE stop with a nonzero exit status.
- EQUIL - Tabulated Zeff and kinetic profiles use monotone (pchip) fits.
- DOCS - Template and example decks document the new inputs (including ``enstep`` in ``equil.in``), ``rmatch.in`` is fully commented, and the a10 example writes the PENTRC ASCII tables.

Fixes
------
- PENTRC - ``fcgl`` recovers all of kappa.xi_perp from dB/B and div(xi_perp), not half.
- PENTRC - Removed a spurious factor of the major radius from the general-aspect-ratio ``omega_D`` (``fgar``, ``tgar``, ``pgar``, ``*mm``); recompute results from these methods.
- PENTRC - Removed spurious factors of ``sqrt(B0)`` from the ``clar`` ``omega_b`` and ``B0`` from its ``omega_D``.
- PENTRC - Fixed OpenMP data races (energy contour flag, ``tpsi``, shared splines), so results no longer depend on the thread count.
- PENTRC - Fixed an out-of-bounds ``bpts`` array and a kinetic-spline length mismatch in ``tpsi``.
- PENTRC - ``rlar`` uses the evaluated ``fnml_f`` in the pitch integrand, so it no longer returns zero torque.
- PENTRC, DCON - Large automatic arrays in the pitch and energy integrations are heap allocated, so OpenMP kinetic runs no longer overflow the worker-thread stacks.
- PENTRC, DCON - The kinetic Euler-Lagrange matrices (``kin_flag``) conjugate the bounce-averaged W vectors, so they act on the actual perturbation instead of its up-down mirror image (up-down asymmetric equilibria only).
- PENTRC - Collision frequencies use a monotone (pchip) fit with temperature and logLambda floors, so they stay positive near a cold edge; non-positive densities and negative temperatures stop with an error.
- PENTRC - Parallel-velocity bounce-root splines use periodic end conditions.
- PENTRC - The ``omega_D`` integrand is no longer padded with ``omega_b`` values near the trapped/passing boundary.
- PENTRC - Fixed anti-damping with ``nutype="small"`` for negative ``omega_E``.
- PENTRC - Radial torque-density export writes the last surface's own value.
- DCON - Fixed the kinetic matrix ``kinfac`` for ``method > 0`` and race conditions in the kinetic matrix loop.
- DCON - Plasmas with q0 >= qa are no longer truncated in the deep core.
- DCON - Fixed out-of-bounds singular-surface access when the q range has no rational surface.
- DCON, RDCON, STRIDE - With ``sas_flag``, ``dmlim`` no longer carries ``qlim`` past ``qhigh``.
- DCON, RDCON, STRIDE - The free-boundary vacuum matrix and eigenvector norm are put in the plasma's Fourier frame, so delta-W no longer depends on the theta=0 origin.
- DCON, RDCON - The singular reduction at a rational surface is kept when a norm reset falls just before it.
- DCON - ``W_v_eigenvalue`` holds the correct real eigenvalues.
- DCON, RDCON, STRIDE - Frobenius exponents at rational surfaces include the trace of the resonant block (no spurious imaginary parts on the Delta' diagonal); Mercier alpha is computed in double precision.
- DCON, RDCON - Fixed Taylor-series slot typos in ``sing_mmat`` (F 6th order; DCON G 3rd order for ``sing_order >= 3``).
- DCON, RDCON, STRIDE, GPEC - Fixed LAPACK and ZVODE argument types, sizes and work arrays (``W_xe_*`` used a stale matrix; ``s_P`` precision; uninitialized ``out_evals`` columns).
- RDCON - Added missing zero initializations, safe deallocations and the correct ``restore_us`` condition in the Galerkin solver; ``vmat.bin`` is written correctly.
- RDCON - Galerkin F is built from the same factored spline as the asymptotics.
- RDCON - Galerkin vacuum matrix keeps its imaginary part; the LU solve is diagonally scaled; a failed Cholesky falls back to LU.
- RDCON - Galerkin elements below a singular surface are sized with ``|n q1|`` instead of ``|q1|``.
- RDCON - Fixed the ``delta_gw.out`` crash with 25 or more singular surfaces and the netCDF failure with ``rpec_flag=f``; ODE integration checks that it reaches psilim.
- RDCON - The netCDF Galerkin cell grid is initialized when the Galerkin output flags are off.
- RDCON, GPEC - Removed spurious shielding at the RPEC control surface; fixed the ``globalsol.bin`` point-count mismatch that blocked field penetration at rational surfaces.
- RMATCH - Removed spurious "imap error" messages with the default ``inps_type="inps"``.
- RMATCH - Fixed an out-of-bounds read in ``match_solution`` when ``msing`` is less than RDCON's surface count, an uninitialized ``delta`` and an undefined-behavior logical in ``inps_ua``.
- GPEC - Control-surface normals ``R_n`` and ``z_n`` use the current geometry derivatives.
- GPEC - ``xbnormal`` normalizes each surface with its own major radius.
- GPEC - pmodb netCDF functions (``b_eul_fun``, ``b_lag_fun``, ``Bdivxi_perp_fun``, ``Bkappaxi_perp_fun``) have correct imaginary parts and include the theta=1 row.
- GPEC - pmodb divergence output honours ``reg_flag=f``; the pmodb diagnostic projection uses the current Jacobian.
- GPEC, PENTRC - pmodb ``Bdivxi_perp`` and ``Bkappaxi_perp`` are formed directly from xi_perp instead of subtracting a Fourier-truncated parallel displacement; PENTRC's ``pmodb_file`` reader uses them instead of ``b_lag``.
- GPEC - Fixed a use-after-free of the local coupling vectors; exported singular vectors use the adjoint area weighting; ``C_l<tag>_x_out`` holds the local coupling matrix.
- GPEC - ``q1`` at resonant surfaces is computed instead of set to 0; fixed out-of-bounds access with a single rational surface, a netCDF buffer error when ``coil_flag`` is set without coils, the ``gpout_xbnormal`` array length and long_name attributes written to the wrong files.
- STRIDE - Fixed an OpenMP data race in ZLANGE, stack overflow at large ``mpert``, an out-of-bounds write and a double allocation.
- STRIDE - The interval from the axis to the first rational surface is always split into at least 8 pieces.
- VACUUM - Fixed an off-by-one between ``mthsurf`` and ``mthvac``, a divide by zero in ``defglo`` and out-of-bounds writes; observer-grid coordinates are double precision.
- COIL - Fixed undefined spline memory and an out-of-range write in ``field.F``; zero-amplitude shifts no longer distort the coils.
- EQUIL - Fixed spline "extrap" end derivatives on clustered grids and real-valued classic complex spline coefficients; ``spline_roots`` finds roots near cell endpoints.
- EQUIL - Fixed uninitialized ``ro``/``zo`` for EFIT input, ``eq_type="dump"`` for direct equilibria, an out-of-bounds write for g-files without a q column, a crash re-splining inverse equilibria, and the CHEASE ASCII reader; the field-line integration failure message is printed in full.
- EQUIL - Direct equilibria are now sampled exactly at the common straight-fieldline angles instead of interpolated from each surface's own field-line trace.
- EQUIL - ``grid_type="original"`` honours ``psilow`` and ``psihigh``.
- INSTALL - The embedded version string uses release tags only; bundled dependencies honour ``FC``/``CC`` and embed an rpath; xdraw compiles with strict C99/C11 compilers; added missing build dependencies.


GPEC v1.5.7
===========

Fixes
------
- ALL - BUMP - Update version to 1.5.7 and document changes in release notes
- INSTALL - Remove harvest submodule reference from git tree
- ALL - REMOVE - Remove all harvest dependencies and capabilities


GPEC v1.5.6
===========

Fixes
------
- EQUIL - Fixes bug in the reading of the f profile from chease2 equilibrium files


GPEC v1.5.5
===========

Fixes
------
- DCON - Fixes bug in the calculation of the kinetic torque matrix when using either the electron terms or method==1.


GPEC v1.5.4
===========

Fixes
------
- STRIDE - Fixes reformation of splines with psilow, psilim from truncation logic in STRIDE


GPEC v1.5.3
===========

Fixes
------
- INSTALL - Respects user's env when choosing NETCDFINC


GPEC v1.5.2
===========

Fixes
------
- GPEC - Fixes bug in flux-field conversion of P_xe
- GPEC - Fixes orientation of L, Lambda, P, rho, W_e, W_xe, R_xe, P_xe matrices in netcdf outputs
- GPEC - Fixes bug that had scrambled C_*_x_out netcdf output


GPEC v1.5.1
===========

Fixes
------
- GPEC - Fixes bug in jacfac diagnostic, now benchmarked with MEPHIT
- GPEC - Fixes bug in terminal output of q during vsingfld


GPEC v1.5.0
===========

Adds
------
- GPEC - Outputs new J. Callen threshold (w_isl_crit) for edge island penetration, controlled by the singtresh_flag. This requires PENTRC kinetic file inputs.
- GPEC - Outputs new coil to resonant field coupling matrix
- GPEC - Outputs vector potential in cylindrical coordinates
- GPEC - Outputs Delta and Jbgradpsi (bwp in ascii) to netcdf outputs, as well as PEST coordinate Jbgradpsi_pest
- GPEC - Outputs jafac and delpsi profiles via an optional diagnostic routine
- COIL - Adds leading DIII-D M-coil designs
- RMATCH - Adds variable inps_eps=1e-4,1e-7,1e-8 for automatic determination of computational elements, removing xmax_method and replacing xfac with inps_xfac
- RMATCH - Includes production code used to form initialization scans published by Froese & Brennan for General Fusion plasmas
- RMATCH - Adds TRANSP interface module

Changes
--------
- VACUUM - Implicit uses of integer variables made explicit along with initialization of variables at zero (towards static compilation)
- GPEC - Improves terminal output with table-like formatting and a change from reported averages to vector magnitudes
- COIL - Updates coil dat file reading to use arbitrary x,y,z numerical formatting for higher accuracy
- COIL - Accepts unformatted header lines in coil dat files (still needs to be 3 integers and a float)
- PENTRC - Reduces number of LSODE warnings in energy and pitch integrals
- PYPEC - Makes necessary changes to switch to python3 and updates portal partitions

Fixes
------
- EQUIL - Adds proper zmid displacement to z grid from g-files (previously, z of equilibria were off by this value - effectively moving the magnetic axis to z=0 in all cases)
- DCON - Forces euler.bin outputs at critical steps on either side of singular surfaces, enabling safe use of euler_stride with GPEC
- DCON - Scans initial guess for inboard and outboard separatrix positions for robustness to deep magnetic wells within g-file domains
- STRIDE - Forces reforming of equil splines with user modified bounds
- GPEC - Fixed units of the perturbed field and displacement in netcdf output (previously, the former's "m" units were overwritting the field units)
- GPEC - Fixes units of the vacuum resonant flux in the netcdf output
- GPEC - Initializes the coil_indmat values to zero to prevent erroneous values beyond cmlow-cmhigh
- GPEC - Truncates rzphi extrapolations to the core when psilow or qlow sets the start far from the magnetic axis, preventing non-physical values from swamping the output
- VACUUM - Eliminates warning due to searchx declaring incorrect x-array size
- VACUUM - Adds ishape=42 to residual calculation in kernel
- INSTALL - Makes run_examples bash script explicitly use installation's exes


GPEC v1.4.4
===========

Fixes
------
- DCON - Fixes output of plasma energy eigenvectors in netcdf


GPEC v1.4.3
===========

Fixes
------
- VACUUM - Respects the namelist a and lspark when determining whether to write a wall geometry file output
- DCON - Properly spaces the vacuum matrices used when finding the peak dW with psiedge < psilim

Changes
--------
- DCON - dW_edge in the netcdf output now includes the imaginary component of the least stable eigenvalue


GPEC v1.4.2
===========

Fixes
------
- PENTRC - Respects the namelist verbose setting throughout the program

Changes
--------
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
- GPEC - Renames netcdf dimensions \*_rat to \*_rational for clarity

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
