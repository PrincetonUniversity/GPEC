c=======================================================================
c     MODULE layerinputs_mod
c
c     Reads STRIDE NetCDF output and constructs the per-surface
c     input arrays required by the SLAYER layer-physics solver.
c
c     Subprograms contained:
c       1. read_stride_netcdf_diagonal -- read STRIDE NetCDF, extract
c              Deltaprime diagonal, geometry, and equilibrium scalars.
c       2. issurfint                   -- surface integral by simple
c              quadrature (adapted from EQUIL).
c       3. build_inputs                -- master driver: reads kinetic
c              profiles, evaluates params(), and populates
c              slayer_inputs_type for every rational surface.
c=======================================================================
      MODULE layerinputs_mod

      USE inputs,   ONLY : read_kin, read_equil, kin, chi1
      USE spline_mod, ONLY : spline_alloc, spline_eval, spline_type,
     $                       spline_dealloc, spline_int, spline_fit
      USE sglobal_mod          ! SLAYER global scalars, types, constants
      USE params_mod           ! params() -- compute derived layer params
      USE netcdf               ! NetCDF Fortran bindings
      USE equil_mod, ONLY : equil_read, rzphi, twopi, ro, zo, sq
      USE bicube_mod, ONLY : bicube_eval_external, bicube_type
      USE slayer_netcdf_mod    ! sl_check(), SLAYER NetCDF output

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. read_stride_netcdf_diagonal.
c     Read the STRIDE NetCDF file and extract the Deltaprime matrix,
c     rational-surface geometry (q, psi_n, shear, dgeo, dr), and
c     equilibrium scalars (R0, B_t0, psi0, m_psi, n, resm).
c
c     BUG FLAG 1 -- NetCDF variable IDs (dp_id, qr_id, ...) are
c       declared but never initialised via nf90_inq_varid before
c       being passed to nf90_inquire_attribute.  The attribute
c       reads for ro, bt0, psio, mpsi, n use *_id variables that
c       are still zero at that point; the calls appear to succeed
c       only because the NetCDF library treats the id argument as
c       a global-attribute flag when it is NF90_GLOBAL (=0).
c       Suggested fix: either replace the id arguments with
c       NF90_GLOBAL explicitly, or move the nf90_inq_varid calls
c       above the nf90_inquire_attribute calls.
c-----------------------------------------------------------------------
      SUBROUTINE read_stride_netcdf_diagonal(ncfile,msing,dp_mat,
     $   Re_dp_diagonal,Im_dp_diagonal,q_rational,psi_n_rational,dgeo,
     $   shear,r_o,my_bt0,my_psio,dr_vals,mpsi,nn,resm)

c --- arguments (all INTENT(OUT) except ncfile)
      CHARACTER(512), INTENT(IN) :: ncfile           ! path to STRIDE NetCDF
      INTEGER, INTENT(OUT)       :: msing            ! number of singular surfaces
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: dp_mat
                                                      ! Deltaprime matrix (msing x msing x 2)
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $      Re_dp_diagonal,    ! Re diag(Deltaprime)
     $      Im_dp_diagonal     ! Im diag(Deltaprime)
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $      q_rational,        ! rational-surface q values
     $      psi_n_rational,    ! normalised psi at each surface
     $      shear,             ! magnetic shear
     $      dgeo               ! geometric delta (Shafranov shift)
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $      r_o,               ! major radius R0
     $      my_bt0,            ! toroidal field B_t0
     $      my_psio,           ! poloidal flux psi_0
     $      mpsi,              ! poloidal mode number array
     $      dr_vals            ! radial width dr at each surface
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $      nn,                ! toroidal mode number(s)
     $      resm               ! resonant poloidal mode numbers

c --- locals
      REAL(r8), DIMENSION(:), ALLOCATABLE :: msing_arr  ! temp for reading msing attribute
      INTEGER(kind=nf90_int) :: ncid, stat             ! NetCDF file id / return status
      INTEGER(kind=nf90_int) :: r_dim_id, r_dim        ! dimension id / length (unused)
      INTEGER(kind=nf90_int) :: dp_id, qr_id, pr_id    ! variable ids
      INTEGER(kind=nf90_int) :: dgeo_id, shear_id      ! variable ids
      INTEGER(kind=nf90_int) :: ro_id, bt0_id, psio_id ! attribute ids (see BUG FLAG 1)
      INTEGER(kind=nf90_int) :: mpsi_id, msing_id      ! attribute ids
      INTEGER(kind=nf90_int) :: nn_id, resm_id, drr_id ! variable / attribute ids
      INTEGER(kind=nf90_int), DIMENSION(1) :: start, count  ! NetCDF hyperslab
      INTEGER :: i                                      ! loop index
      INTEGER :: bt0_len, ro_len, psio_len              ! attribute lengths
      INTEGER :: mpsi_len, msing_len, nn_len, dr_len    ! attribute lengths

c-----------------------------------------------------------------------
c     open the STRIDE NetCDF file and read dimension / attribute data.
c-----------------------------------------------------------------------
      WRITE(*,*) '$^$ opening netcdf file', ncfile

      stat = nf90_open(path=ncfile, mode=NF90_WRITE, ncid=ncid)
      CALL sl_check(stat)

c --- read msing (number of singular surfaces) from global attribute
      stat = nf90_inquire_attribute(ncid, msing_id, 'msing',
     $        len = msing_len)
      CALL sl_check(stat)
      ALLOCATE(msing_arr(msing_len))
      stat = nf90_get_att(ncid, msing_id, 'msing', msing_arr)
      CALL sl_check(stat)
      msing = INT(msing_arr(1))

c --- allocate output arrays sized by msing
      ALLOCATE(Re_dp_diagonal(msing), q_rational(msing),
     $         psi_n_rational(msing), shear(msing), dgeo(msing),
     $         resm(msing), Im_dp_diagonal(msing), dr_vals(msing))
      ALLOCATE(dp_mat(msing, msing, 2))

c --- read lengths of scalar / small-array global attributes
      stat = nf90_inquire_attribute(ncid, ro_id,   'ro',   len=ro_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, bt0_id,  'bt0',  len=bt0_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, psio_id, 'psio', len=psio_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, mpsi_id, 'mpsi', len=mpsi_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, nn_id,   'n',    len=nn_len)
      CALL sl_check(stat)

      ALLOCATE(my_bt0(INT(bt0_len)), r_o(INT(ro_len)),
     $         my_psio(INT(psio_len)),
     $         mpsi(INT(mpsi_len)), nn(INT(nn_len)))

c --- obtain NetCDF variable IDs
      stat = nf90_inq_varid(ncid, 'Delta_prime',   dp_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'q_rational',    qr_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'psi_n_rational', pr_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'Delta_geo',     dgeo_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'shear',         shear_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'resm',          resm_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, 'dr_rational',   drr_id)
      CALL sl_check(stat)

c --- read global attributes (equilibrium scalars)
      stat = nf90_get_att(ncid, ro_id,   'ro',   r_o)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, bt0_id,  'bt0',  my_bt0)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, psio_id, 'psio', my_psio)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, mpsi_id, 'mpsi', mpsi)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, nn_id,   'n',    nn)
      CALL sl_check(stat)

c --- read variable data: Deltaprime matrix and 1-D surface arrays
      stat = nf90_get_var(ncid, dp_id, dp_mat, start=(/1,1,1/))
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, qr_id, q_rational)
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, pr_id, psi_n_rational)
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, dgeo_id, dgeo)
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, shear_id, shear)
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, resm_id, resm)
      CALL sl_check(stat)
      stat = nf90_get_var(ncid, drr_id, dr_vals)
      CALL sl_check(stat)

c --- extract diagonal of the complex Deltaprime matrix
      DO i = 1, msing
        Re_dp_diagonal(i) = dp_mat(i, i, 1)
        Im_dp_diagonal(i) = dp_mat(i, i, 2)
      END DO

c --- close file
      stat = nf90_close(ncid)
      CALL sl_check(stat)

      END SUBROUTINE read_stride_netcdf_diagonal
c-----------------------------------------------------------------------
c     subprogram 2. issurfint.
c     Surface integral by simple quadrature, adapted from EQUIL.
c     Computes  int f(theta) * W(theta) d(theta)  where W depends
c     on the weight flag `wegt`:
c       0 = jac * |grad psi|
c       1 = R * jac * |grad psi|      (R-weighted)
c       2 = jac * |grad psi| / R      (1/R-weighted)
c       3 = a * jac * |grad psi|       (minor-radius weighted)
c     If ave==1, the result is divided by the unweighted surface area.
c
c     Geometry is cached via first/fsave/psave to avoid recomputing
c     Jacobians when the same surface is queried repeatedly.
c
c     BUG FLAG 2 -- `first = .FALSE.` is commented out (line after
c       `fsave = fs`).  The caching logic therefore never sets
c       first=.FALSE., so geometry is recomputed on every call
c       even when psi and fs are unchanged.  Uncomment the line
c       or remove the caching branch entirely.
c
c     BUG FLAG 3 -- `z` (toroidal-Z coordinate) is computed inside
c       the geometry loop but never used outside it.  Remove or
c       document its intended purpose.
c-----------------------------------------------------------------------
      FUNCTION issurfint(func,fs,inpsi,wegt,ave,
     $     fsave,psave,jacs,delpsi,inr,ina,first)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
c --- arguments
      INTEGER, INTENT(IN)    :: fs      ! number of poloidal segments
      INTEGER, INTENT(IN)    :: wegt    ! weight flag (0..3)
      INTEGER, INTENT(IN)    :: ave     ! 1 = return surface average
      REAL(r8), INTENT(IN)   :: inpsi   ! normalised psi of surface
      REAL(r8), DIMENSION(0:fs), INTENT(IN) :: func  ! integrand array

      LOGICAL, INTENT(INOUT) :: first   ! .TRUE. on first call (see BUG FLAG 2)
      INTEGER, INTENT(INOUT) :: fsave   ! cached fs
      REAL(r8), INTENT(INOUT) :: psave  ! cached psi
      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: jacs   ! Jacobian cache
      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: delpsi ! |grad psi| cache
      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: inr    ! R(theta) cache
      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: ina    ! a(theta) cache
c --- return value
      REAL(r8) :: issurfint
c --- locals
      INTEGER  :: itheta                            ! poloidal loop index
      INTEGER  :: ix, iy                            ! bicube grid hints
      REAL(r8) :: rfac, ineta, injac, inarea        ! geometry intermediates
      REAL(r8), DIMENSION(1,2)  :: w                ! gradient components
      REAL(r8), DIMENSION(0:fs) :: z                ! Z coords (UNUSED -- BUG FLAG 3)
      REAL(r8), DIMENSION(0:fs) :: thetas           ! normalised theta grid
      REAL(r8), DIMENSION(4) :: rzphi_f, rzphi_fx, rzphi_fy
                                                     ! bicube_eval_external outputs
c-----------------------------------------------------------------------
c     compute / cache geometry if surface changed.
c     [bicube_eval_external]: external bicubic interpolation from EQUIL.
c-----------------------------------------------------------------------
      issurfint = 0
      inarea = 0
      ix = 0
      iy = 0

      IF (first .OR. inpsi /= psave .OR. fs /= fsave) THEN
         psave = inpsi
         fsave = fs
         !first = .FALSE.                ! BUG FLAG 2: should be uncommented
         DO itheta = 0, fs
            thetas(itheta) = REAL(itheta, r8) / REAL(fs, r8)
         ENDDO
         DO itheta = 0, fs-1
            CALL bicube_eval_external(rzphi, inpsi, thetas(itheta), 1,
     $           ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
            rfac  = SQRT(rzphi_f(1))
            ineta = twopi * (thetas(itheta) + rzphi_f(2))
            ina(itheta) = rfac
            inr(itheta) = ro + rfac * COS(ineta)
            z(itheta)   = zo + rfac * SIN(ineta)
            injac = rzphi_f(4)
            jacs(itheta) = injac
c           gradient magnitude: |grad psi| from metric components
            w(1,1) = (1 + rzphi_fy(2)) * twopi**2
     $                * rfac * inr(itheta) / injac
            w(1,2) = -rzphi_fy(1) * pi * inr(itheta)
     $                / (rfac * injac)
            delpsi(itheta) = SQRT(w(1,1)**2 + w(1,2)**2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     perform weighted surface integral based on wegt flag.
c-----------------------------------------------------------------------
      IF (wegt == 0) THEN
         DO itheta = 0, fs-1
            issurfint = issurfint
     $           + jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt == 1) THEN
         DO itheta = 0, fs-1
            issurfint = issurfint
     $         + inr(itheta)*jacs(itheta)*delpsi(itheta)*
     $           func(itheta)/fs
         ENDDO
      ELSE IF (wegt == 2) THEN
         DO itheta = 0, fs-1
            issurfint = issurfint
     $           + jacs(itheta)*delpsi(itheta)*
     $             func(itheta)/inr(itheta)/fs
         ENDDO
      ELSE IF (wegt == 3) THEN
         DO itheta = 0, fs-1
            issurfint = issurfint
     $        + ina(itheta)*jacs(itheta)*delpsi(itheta)*
     $          func(itheta)/fs
         ENDDO
      ELSE
         STOP 'ERROR: issurfint wegt must be in [0,1,2,3]'
      ENDIF
c-----------------------------------------------------------------------
c     optionally normalise by unweighted surface area.
c-----------------------------------------------------------------------
      IF (ave == 1) THEN
         DO itheta = 0, fs-1
            inarea = inarea + jacs(itheta)*delpsi(itheta)/fs
         ENDDO
         issurfint = issurfint / inarea
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfint
c-----------------------------------------------------------------------
c     subprogram 3. build_inputs.
c     Master driver: reads the STRIDE NetCDF file, sets up kinetic
c     profiles from the input file, reads the equilibrium, and then
c     loops over every singular surface to compute derived layer
c     parameters via params() and populate slayer_inputs_type.
c
c     BUG FLAG 4 -- many local scalars declared here (lpsi, rpsi,
c       hdist, sbnosurf, spl, sr, my_inpe, tau_i, b_l, v_a, tau_h,
c       rho, tau_v, Qconv, lbeta, qintb, tau_ee_num..chi_par,
c       psitor, rhotor, my_rhotor, my_psitor, rfac, jac, wit)
c       are never used.  They appear to be left over from an earlier
c       version.  Remove to reduce confusion.
c
c     BUG FLAG 5 -- `ising` is declared REAL(r8) but is used as a
c       DO-loop index (integer context).  This works in Fortran but
c       is non-standard in F90+ and may fail with strict compilers.
c       Declare as INTEGER.
c
c     BUG FLAG 6 -- `zeff = 2.0` is hardcoded on every surface
c       (the kin%f(9) alternative is commented out).  If the
c       kinetic file provides Z_eff, this should use it.
c
c     BUG FLAG 7 -- `mrs` and `nrs` are assigned
c       `real(mms,4)` / `real(nns,4)` (single precision) but are
c       declared INTEGER, so the float is silently truncated.
c       They are also never used afterwards.  Remove or fix.
c
c     BUG FLAG 8 -- the local arrays ne_arr, te_arr, ... mu_i_arr,
c       nns_arr, dr_arr, omegas_e_arr, omegas_i_arr are populated
c       inside the loop but never used outside it (the old
c       slayer_netcdf_inputs call is commented out).  Remove or
c       gate behind a diagnostic flag.
c-----------------------------------------------------------------------
      SUBROUTINE build_inputs(infile,ncfile,sl_in)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
c --- arguments
      CHARACTER(512), INTENT(IN) :: infile  ! kinetic input file path
      CHARACTER(512), INTENT(IN) :: ncfile  ! STRIDE NetCDF file path
      TYPE(slayer_inputs_type), INTENT(INOUT) :: sl_in
c --- surface-loop control
      LOGICAL  :: firstsurf                 ! first-call flag for issurfint
      REAL(r8) :: respsi                    ! normalised psi at current surface
      REAL(r8) :: ising                     ! loop index (BUG FLAG 5: should be INTEGER)
c --- unused scalars (BUG FLAG 4 -- remove)
      REAL(r8) :: lpsi, rpsi, hdist, sbnosurf
c --- kinetic profile parameters (for read_kin)
      INTEGER  :: zi, zimp, mi, mimp        ! charge/mass species ids
      REAL(r8) :: nfac, tfac, wefac, wpfac  ! profile scale factors
      REAL(r8) :: e                         ! elementary charge [C]
c --- unused spline temporaries (BUG FLAG 4)
      TYPE(spline_type) :: spl
      TYPE(spline_type) :: sr
c --- mode number workspace
      INTEGER  :: mms, nns                  ! poloidal / toroidal mode nums
      INTEGER  :: mrs, nrs                  ! UNUSED, wrongly typed (BUG FLAG 7)
      INTEGER  :: mpsi                      ! poloidal-flux index from attr
c --- local plasma quantities at current surface
      REAL(r8) :: n_e, t_e, n_i, t_i       ! densities [m^-3], temperatures [eV]
      REAL(r8) :: omega, omega_e, omega_i   ! toroidal & diamagnetic freqs [rad/s]
      REAL(r8) :: my_qval, my_sval          ! safety factor, magnetic shear
      REAL(r8) :: my_bt, my_rs, R_0         ! toroidal field, minor radius, major radius
      REAL(r8) :: my_inpe                   ! UNUSED (BUG FLAG 4)
      REAL(r8) :: zeff                      ! effective charge (hardcoded -- BUG FLAG 6)
      REAL(r8) :: dgeo_val                  ! geometric delta (Shafranov)
      REAL(r8) :: mu_i                      ! ion mass ratio to proton
      REAL(r8) :: dr_val                    ! radial width dr at surface
      REAL(r8) :: l_n, l_t                  ! density / temperature gradient lengths
      REAL(r8) :: gammafac                  ! growth-rate conversion factor
      REAL(r8), DIMENSION(3) :: chi_s       ! chi_perp, chi_tor, kappa
c --- unused derived quantities (BUG FLAG 4 -- left over from params duplication)
      REAL(r8) :: tau_i, b_l, v_a, tau_h, rho, tau_v
      REAL(r8) :: Qconv, lbeta, qintb
      REAL(r8) :: tau_ee_num, tau_ee_denom, tau_ee
      REAL(r8) :: sigma_par_1, sigma_par_2, sigma_par
      REAL(r8) :: tau_perp, Wd, vte
      REAL(r8) :: chi_par_smfp, chi_par_lmfp, chi_par
      INTEGER  :: wit
c --- unused flux-coordinate arrays (BUG FLAG 4)
      REAL(r8), DIMENSION(0:128) :: psitor, rhotor
      REAL(r8), DIMENSION(:), ALLOCATABLE :: my_rhotor, my_psitor
c --- STRIDE data (from read_stride_netcdf_diagonal)
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dp_mat
      REAL(r8), DIMENSION(:), ALLOCATABLE :: Re_dp_diagonal,
     $           Im_dp_diagonal, q_rational, psi_n_rational,
     $           shear, dgeo, r_o, my_bt0, my_psio, mpsi_arr,
     $           dr_vals, dr_arr
c --- local per-surface kinetic arrays (UNUSED outside loop -- BUG FLAG 8)
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ne_arr, te_arr,
     $    ni_arr, ti_arr, zeff_arr, bt_arr, rs_arr, R0_arr,
     $    mu_i_arr, omegas_e_arr, omegas_i_arr
      INTEGER, DIMENSION(:), ALLOCATABLE :: nn, resm, nns_arr
c --- surface-integral workspace
      INTEGER  :: msing, i, mthsurf         ! surface count, loop idx, theta pts
      REAL(r8), DIMENSION(0:512) :: unitfun ! unit function for area integral
      INTEGER  :: fsave                     ! cached fs for issurfint
      REAL(r8) :: psave                     ! cached psi for issurfint
      REAL(r8), DIMENSION(:), ALLOCATABLE :: jacs, delpsi, rsurf, asurf
      REAL(r8) :: rfac, jac, a_surf         ! rfac/jac UNUSED (BUG FLAG 4)
c-----------------------------------------------------------------------
c     read STRIDE NetCDF: Deltaprime matrix, geometry, equilibrium scalars.
c-----------------------------------------------------------------------
      CALL read_stride_netcdf_diagonal(ncfile,msing,dp_mat,
     $           Re_dp_diagonal,Im_dp_diagonal,q_rational,
     $           psi_n_rational,dgeo,shear,r_o,my_bt0,my_psio,dr_vals,
     $           mpsi_arr,nn,resm)

      mpsi = INT(mpsi_arr(1))
      mthsurf = 512             ! poloidal segments for surface integral

c --- allocate slayer_inputs_type arrays (one entry per surface)
      ALLOCATE(sl_in%qval_arr(msing),  sl_in%omegas_arr(msing),
     $  sl_in%omegas_e_arr(msing),     sl_in%dp_matrix(msing,msing),
     $  sl_in%omegas_i_arr(msing),
     $  sl_in%Q_e_arr(msing),          sl_in%Q_i_arr(msing),
     $  sl_in%psi_n_arr(msing),
     $  sl_in%Re_dp_arr(msing),        sl_in%Im_dp_arr(msing),
     $  sl_in%d_crit_arr(msing),       sl_in%P_tor_arr(msing),
     $  sl_in%P_perp_arr(msing),       sl_in%tau_arr(msing),
     $  sl_in%D_norm_arr(msing),
     $  sl_in%d_beta_arr(msing),       sl_in%gammafac_arr(msing),
     $  sl_in%c_beta_arr(msing),       sl_in%lu_arr(msing),
     $  sl_in%Qconv_arr(msing))

c --- allocate local kinetic arrays (diagnostic, BUG FLAG 8)
      ALLOCATE(ne_arr(msing), te_arr(msing), ni_arr(msing),
     $    ti_arr(msing), zeff_arr(msing), bt_arr(msing),
     $    rs_arr(msing), R0_arr(msing), mu_i_arr(msing),
     $    nns_arr(msing), dr_arr(msing),
     $    omegas_e_arr(msing), omegas_i_arr(msing))

      ALLOCATE(jacs(0:mthsurf), delpsi(0:mthsurf),
     $         rsurf(0:mthsurf), asurf(0:mthsurf))
c-----------------------------------------------------------------------
c     set up kinetic profiles and equilibrium.
c     [read_kin]: external, reads kinetic input file into kin spline.
c     [equil_read]: external, reads equilibrium into EQUIL module.
c-----------------------------------------------------------------------
      zi   = 1                  ! main-ion charge
      zimp = 6                  ! impurity charge (carbon)
      mi   = 2                  ! main-ion mass (deuterium)
      mimp = 12                 ! impurity mass
      nfac  = 1.0
      tfac  = 1.0
      wefac = 1.0
      wpfac = 1.0
      e = 1.6021917e-19         ! elementary charge [C]
      chi1 = twopi * my_psio(1) ! total poloidal flux (module-level)

      CALL read_kin(infile,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,.false.)

      CALL equil_read(out_unit)

c     store full complex Deltaprime matrix in sl_in
      sl_in%dp_matrix(:,:) = CMPLX(dp_mat(:,:,1), dp_mat(:,:,2))

c-----------------------------------------------------------------------
c     loop over singular surfaces: evaluate kinetic/equilibrium
c     quantities via spline interpolation, compute derived layer
c     parameters via params(), and populate sl_in arrays.
c-----------------------------------------------------------------------
      DO ising = 1, msing

         respsi = psi_n_rational(ising)    ! normalised psi
         firstsurf = .TRUE.
         unitfun = 1

c        compute flux-surface-averaged minor radius
         a_surf = issurfint(unitfun,mthsurf,respsi,3,1,
     $           fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)

c-----------------------------------------------------------------------
c        evaluate kinetic splines at this surface.
c        [spline_eval]: external, evaluates kin spline at respsi.
c        kin%f(1..5) = n_i, n_e, t_i, t_e, omega  (SI units)
c        kin%f1(1..5) = d/d(psi_n) of the above
c-----------------------------------------------------------------------
         CALL spline_eval(kin,respsi,1)

c        diamagnetic frequencies (rad/s)
         omega_i = -twopi*kin%f(3)*kin%f1(1)/(e*zi*chi1*kin%f(1))
     $             -twopi*kin%f1(3)/(e*zi*chi1)
         omega_e =  twopi*kin%f(4)*kin%f1(2)/(e*chi1*kin%f(2))
     $             +twopi*kin%f1(4)/(e*chi1)

         sl_in%omegas_e_arr(ising) = omega_e
         sl_in%omegas_i_arr(ising) = omega_i

c        extract local plasma quantities from spline
         n_e = kin%f(2)
         t_e = kin%f(4) / e          ! convert J -> eV
         n_i = kin%f(1)
         t_i = kin%f(3) / e

         zeff = 2.0                   ! hardcoded (BUG FLAG 6)

         omega    = kin%f(5)
         my_qval  = q_rational(ising)
         my_sval  = shear(ising)
         dgeo_val = dgeo(ising)
         my_bt    = my_bt0(1)
         my_rs    = a_surf
         R_0      = r_o(1)
         mu_i     = 2.0               ! deuterium
         dr_val   = dr_vals(ising)

c        transport coefficients from caller-provided arrays
         chi_s(1) = sl_in%chi_p_arr(ising) ! chi_perp
         chi_s(2) = sl_in%chi_t_arr(ising) ! chi_tor
         chi_s(3) = sl_in%kappa_arr(ising) ! kappa (thermal cond.)

c        store local kinetic arrays (BUG FLAG 8: unused)
         ne_arr(ising)   = n_e
         te_arr(ising)   = t_e
         ni_arr(ising)   = n_i
         ti_arr(ising)   = t_i
         zeff_arr(ising) = zeff
         bt_arr(ising)   = my_bt
         rs_arr(ising)   = my_rs
         R0_arr(ising)   = R_0
         mu_i_arr(ising) = mu_i

         mms = resm(ising)
         nns = nn(1)
         mrs = real(mms, 4)           ! BUG FLAG 7: float -> integer truncation
         nrs = real(nns, 4)
         nns_arr(ising) = nn(1)
         nr = nn(1)                   ! module-level toroidal mode number

         l_n = 0.0
         l_t = 0.0

c-----------------------------------------------------------------------
c        compute derived layer parameters.
c        [params]: external (params_mod), sets module-level globals
c          tau, tau_r, tauk, lu, c_beta, d_beta, D_norm, P_perp,
c          P_tor, dc_tmp, etc. in sglobal_mod.
c-----------------------------------------------------------------------
         CALL params(n_e,t_e,t_i,omega,chi_s,dr_val,dgeo_val,
     $        l_n,l_t,my_qval,my_sval,my_bt,my_rs,R_0,mu_i,
     $        zeff,.false.)

c        growth-rate conversion factor: Deltaprime -> gamma
         gammafac = (my_rs * Re_dp_diagonal(ising)) / tau_r

c-----------------------------------------------------------------------
c        populate sl_in for this surface from params() globals.
c-----------------------------------------------------------------------
         sl_in%qval_arr(ising)    = INT(my_qval)
         sl_in%lu_arr(ising)      = lu
         sl_in%Q_e_arr(ising)     = -tauk * omega_e
         sl_in%Q_i_arr(ising)     = -tauk * omega_i
         sl_in%c_beta_arr(ising)  = c_beta
         sl_in%d_beta_arr(ising)  = d_beta
         sl_in%D_norm_arr(ising)  = D_norm
         sl_in%tau_arr(ising)     = tau
         sl_in%omegas_arr(ising)  = omega
         sl_in%psi_n_arr(ising)   = respsi
         sl_in%gammafac_arr(ising)= gammafac
         sl_in%Re_dp_arr(ising)   = Re_dp_diagonal(ising)
         sl_in%Im_dp_arr(ising)   = Im_dp_diagonal(ising)
         sl_in%d_crit_arr(ising)  = dc_tmp
         sl_in%P_perp_arr(ising)  = P_perp
         sl_in%P_tor_arr(ising)   = P_tor
         sl_in%Qconv_arr(ising)   = tauk
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod