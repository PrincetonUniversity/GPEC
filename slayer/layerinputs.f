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
c     All global-attribute reads use NF90_GLOBAL explicitly;
c     variable IDs are obtained via nf90_inq_varid before use.
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
      INTEGER(kind=nf90_int) :: resm_id, drr_id         ! variable ids
      INTEGER(kind=nf90_int), DIMENSION(1) :: start, count  ! NetCDF hyperslab
      INTEGER :: i                                      ! loop index
      INTEGER :: bt0_len, ro_len, psio_len              ! attribute lengths
      INTEGER :: mpsi_len, msing_len, nn_len, dr_len    ! attribute lengths

c-----------------------------------------------------------------------
c     open the STRIDE NetCDF file and read dimension / attribute data.
c-----------------------------------------------------------------------
      stat = nf90_open(path=ncfile, mode=NF90_NOWRITE,
     $                  ncid=ncid)
      CALL sl_check(stat)

c --- read msing (number of singular surfaces) from global attribute
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'msing',
     $        len = msing_len)
      CALL sl_check(stat)
      ALLOCATE(msing_arr(msing_len))
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'msing', msing_arr)
      CALL sl_check(stat)
      msing = INT(msing_arr(1))

c --- allocate output arrays sized by msing
      ALLOCATE(Re_dp_diagonal(msing), q_rational(msing),
     $         psi_n_rational(msing), shear(msing), dgeo(msing),
     $         resm(msing), Im_dp_diagonal(msing), dr_vals(msing))
      ALLOCATE(dp_mat(msing, msing, 2))

c --- read lengths of scalar / small-array global attributes
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'ro',
     $        len=ro_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'bt0',
     $        len=bt0_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'psio',
     $        len=psio_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'mpsi',
     $        len=mpsi_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'n',
     $        len=nn_len)
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
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'ro',   r_o)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'bt0',  my_bt0)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'psio', my_psio)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'mpsi', mpsi)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, NF90_GLOBAL, 'n',    nn)
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

      LOGICAL, INTENT(INOUT) :: first   ! .TRUE. on first call for caching
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
         first = .FALSE.
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
      INTEGER  :: ising                     ! loop index over singular surfaces
c --- kinetic profile parameters (for read_kin)
      INTEGER  :: zi, zimp, mi, mimp        ! charge/mass species ids
      REAL(r8) :: nfac, tfac, wefac, wpfac  ! profile scale factors
      REAL(r8) :: e                         ! elementary charge [C]
c --- mode number workspace
      INTEGER  :: mms, nns                  ! poloidal / toroidal mode nums
      INTEGER  :: mpsi                      ! poloidal-flux index from attr
c --- local plasma quantities at current surface
      REAL(r8) :: n_e, t_e, n_i, t_i       ! densities [m^-3], temperatures [eV]
      REAL(r8) :: omega, omega_e, omega_i   ! toroidal & diamagnetic freqs [rad/s]
      REAL(r8) :: my_qval, my_sval          ! safety factor, magnetic shear
      REAL(r8) :: my_bt, my_rs, R_0         ! toroidal field, minor radius, major radius
      REAL(r8) :: zeff                      ! effective charge from kin%f(9)
      REAL(r8) :: dgeo_val                  ! geometric delta (Shafranov)
      REAL(r8) :: mu_i                      ! ion mass ratio to proton
      REAL(r8) :: dr_val                    ! radial width dr at surface
      REAL(r8) :: l_n, l_t                  ! density / temperature gradient lengths
      REAL(r8) :: gammafac                  ! growth-rate conversion factor
      REAL(r8), DIMENSION(3) :: chi_s       ! chi_perp, chi_tor, kappa
c --- STRIDE data (from read_stride_netcdf_diagonal)
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dp_mat
      REAL(r8), DIMENSION(:), ALLOCATABLE :: Re_dp_diagonal,
     $           Im_dp_diagonal, q_rational, psi_n_rational,
     $           shear, dgeo, r_o, my_bt0, my_psio, mpsi_arr,
     $           dr_vals, dr_arr
c --- local per-surface kinetic arrays (for future NetCDF diagnostic output)
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
      REAL(r8) :: a_surf                     ! flux-surface-averaged minor radius
c --- Fitzpatrick (r-based) shear workspace
      REAL(r8) :: a_surf_p, a_surf_m         ! a_surf at psiN +/- h
      REAL(r8) :: da_dpsiN                   ! da_surf/dpsiN (Jacobian)
      REAL(r8) :: dpsi_h                     ! finite-diff step for Jacobian
      REAL(r8) :: s_fitz                     ! Fitzpatrick shear r*dq/dr/q
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

c --- allocate local kinetic arrays (for future NetCDF diagnostic output)
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

      dpsi_h = 0.002               ! finite-difference step for Jacobian

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

c        compute Jacobian da_surf/dpsiN by central difference
         a_surf_p = issurfint(unitfun,mthsurf,
     $        MIN(respsi+dpsi_h, REAL(1.0,r8)),3,1,
     $        fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)
         a_surf_m = issurfint(unitfun,mthsurf,
     $        MAX(respsi-dpsi_h, REAL(0.001,r8)),3,1,
     $        fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)
         da_dpsiN = (a_surf_p - a_surf_m)
     $        / (MIN(respsi+dpsi_h, REAL(1.0,r8))
     $         - MAX(respsi-dpsi_h, REAL(0.001,r8)))

c-----------------------------------------------------------------------
c        evaluate kinetic splines at this surface.
c        [spline_eval]: external, evaluates kin spline at respsi.
c        kin%f(1..5) = n_i, n_e, t_i, t_e, omega  (SI units)
c        kin%f1(1..5) = d/d(psi_n) of the above
c-----------------------------------------------------------------------
         CALL spline_eval(kin,respsi,1)

c        diamagnetic frequencies (rad/s) from GPEC kinetic splines.
c        These compute the ELECTRON diamagnetic frequency directly.
         omega_e =  twopi*kin%f(4)*kin%f1(2)/(e*chi1*kin%f(2))
     $             +twopi*kin%f1(4)/(e*chi1)
         omega_i = -twopi*kin%f(3)*kin%f1(1)/(e*zi*chi1*kin%f(1))
     $             -twopi*kin%f1(3)/(e*zi*chi1)
         sl_in%omegas_e_arr(ising) = omega_e
         sl_in%omegas_i_arr(ising) = omega_i

c        extract local plasma quantities from spline
         n_e = kin%f(2)
         t_e = kin%f(4) / e          ! convert J -> eV
         n_i = kin%f(1)
         t_i = kin%f(3) / e

c        Z_eff: kinetic spline value may be incorrect if ni=ne in gpeckf
c        (quasi-neutrality assumption makes Zeff=1). Override to 2.0
c        for deuterium plasma with carbon impurities (matching TJ).
c        TODO: fix gpeckf generation to include proper ni for Zeff,
c        or read Zeff from a namelist parameter.
         zeff = 2.0

         omega    = kin%f(5)
         my_qval  = q_rational(ising)
         my_sval  = shear(ising)
         dgeo_val = dgeo(ising)
         my_bt    = my_bt0(1)
         my_rs    = a_surf
         R_0      = r_o(1)
         mu_i     = 2.0               ! deuterium
         dr_val   = dr_vals(ising)

c        convert STRIDE shear (psiN-based) to Fitzpatrick shear (r-based).
c        s_Fitz = s_psiN * r_s / (psiN * da_surf/dpsiN)
         s_fitz = my_sval * my_rs / (respsi * da_dpsiN)

c        transport coefficients from caller-provided arrays.
c        guard: arrays may be smaller than msing (e.g. from
c        fixed-size namelist); reuse last element if exceeded.
         i = MIN(ising, SIZE(sl_in%chi_p_arr))
         chi_s(1) = sl_in%chi_p_arr(i)
         i = MIN(ising, SIZE(sl_in%chi_t_arr))
         chi_s(2) = sl_in%chi_t_arr(i)
         i = MIN(ising, SIZE(sl_in%kappa_arr))
         chi_s(3) = sl_in%kappa_arr(i)

c        store local kinetic arrays (for future NetCDF diagnostic output)
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
         nns_arr(ising) = nn(1)
         nr = nn(1)                   ! module-level toroidal mode number

         l_n = 0.0
         l_t = 0.0

c-----------------------------------------------------------------------
c        compute derived layer parameters using Fitzpatrick (r-based)
c        shear. params() sees s_fitz as `sval`, so lu, tauk, D_norm,
c        dc_tmp, etc. are all Fitzpatrick-consistent. Gradient lengths
c        are zero here; params() skips its own omega_e/omega_i and we
c        set Q_e/Q_i below from the spline-derived frequencies.
c-----------------------------------------------------------------------
         CALL params(n_e,t_e,t_i,omega,chi_s,dr_val,dgeo_val,
     $        l_n,l_t,my_qval,s_fitz,my_bt,my_rs,R_0,mu_i,
     $        zeff,.false.)

c        growth-rate conversion factor: Deltaprime -> gamma
         gammafac = (my_rs * Re_dp_diagonal(ising)) / tau_r

c-----------------------------------------------------------------------
c        populate sl_in for this surface from params() globals.
c        Q_e/Q_i use spline-derived omega_e/omega_i normalised by
c        tauk (= Fitzpatrick S^(1/3) * tau_H, from params()).
c-----------------------------------------------------------------------
         sl_in%qval_arr(ising)    = INT(my_qval)
         sl_in%lu_arr(ising)      = lu
         sl_in%Q_e_arr(ising)     = -tauk * omega_e
         sl_in%Q_i_arr(ising)     =  tauk * omega_i
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