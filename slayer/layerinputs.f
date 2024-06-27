      MODULE layerinputs_mod

      USE inputs, ONLY : read_kin,kin,chi1
      !USE params
      !USE dcon_interface_mod, ONLY : respsi
      !USE direct_mod, ONLY : direct_run
      USE spline_mod, ONLY : spline_alloc,spline_eval,spline_type,
     $                       spline_dealloc
      USE sglobal_mod, ONLY: m_p, chag, lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,mu0,r8, ! NOT out_unit
     $   eta,visc,rho_s,lu,omega_e,omega_i,delta_n,Q
      USE netcdf


      ! STILL NEED ro AND bt0, GLOBAL GPEC VARIABLES

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE check(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         !STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c-----------------------------------------------------------------------
c     Read and build equilibrium inputs
c-----------------------------------------------------------------------
        SUBROUTINE read_stride_netcdf_diagonal(ncfile, r_dim,
     $   dp_diagonal, q_rational, psi_n_rational, shear,
     $   r0, bt0)

        !USE netcdf   ! NetCDF module for Fortran
        !USE stride_netcdf_mod ! For the 'check' subroutine (error handling)

        ! Input/Output Arguments
        CHARACTER(128), INTENT(IN) :: ncfile
        REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: dp_diagonal
        REAL(r8), DIMENSION(:), INTENT(OUT) :: q_rational,
     $  psi_n_rational, shear
        REAL(r8), INTENT(OUT) :: r0, bt0

        ! Internal Variables
        INTEGER(kind=nf90_int) :: ncid, stat, r_dim_id, r_dim,
     $  dp_id, qr_id, pr_id, shear_id, r0_id, bt0_id  ! Explicit kind for NetCDF variables
        INTEGER(kind=nf90_int), DIMENSION(1) :: start, count ! Explicit kind for NetCDF variables
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: delta_prime
        INTEGER :: i

        ! Open the NetCDF file
        stat = nf90_open(path=ncfile,mode=NF90_WRITE,ncid=ncid)
        WRITE(*,*)"ncfile=",ncfile
        CALL check(stat)  ! Error handling

        ! Get Dimension Information
        stat = nf90_inq_varid(ncid, "r", r_dim_id)
        CALL check(stat)
        stat = nf90_get_var(ncid, r_dim_id, r_dim)  ! Use len= for dimension
        CALL check(stat)

        ! Allocate Arrays (based on dimension)
        ALLOCATE(dp_diagonal(r_dim))
        !ALLOCATE(delta_prime(r_dim, r_dim))

        ! Get Variable IDs
        stat = nf90_inq_varid(ncid, "Delta_prime", dp_id)
        CALL check(stat)
        stat = nf90_inq_varid(ncid, "q_rational", qr_id)
        CALL check(stat)
        stat = nf90_inq_varid(ncid, "psi_n_rational", pr_id)
        CALL check(stat)
        stat = nf90_inq_varid(ncid, "shear", shear_id)
        CALL check(stat)
        stat = nf90_inq_varid(ncid, "bt0", bt0_id)
        CALL check(stat)
        stat = nf90_inq_varid(ncid, "r0", r0_id)
        CALL check(stat)

        ! Read Data from NetCDF File
        ! Set up start and count for reading only the diagonal
        start(1) = 1
        count(1) = 1

        ! Read the diagonal of delta_prime. The results will be put on a 1D temporary array.
        stat = nf90_get_var(ncid, dp_id, delta_prime)
        CALL check(stat)
        ! Read 1D variables
        stat = nf90_get_var(ncid, qr_id, q_rational)
        CALL check(stat)
        stat = nf90_get_var(ncid, pr_id, psi_n_rational)
        CALL check(stat)
        stat = nf90_get_var(ncid, shear_id, shear)
        CALL check(stat)
        stat = nf90_get_var(ncid, bt0_id, bt0)
        CALL check(stat)
        stat = nf90_get_var(ncid, r0_id, r0)
        CALL check(stat)

        ! Extract Diagonal, with 3rd index signifying REAL part
        DO i = 1, r_dim
            dp_diagonal(i) = REAL(delta_prime(i, i, 1))
        END DO

        ! Clean Up
        DEALLOCATE(delta_prime)
        stat = nf90_close(ncid)
        CALL check(stat)

        END SUBROUTINE read_stride_netcdf_diagonal
c-----------------------------------------------------------------------
c     subprogram 1. build_inputs.
c     compute
c-----------------------------------------------------------------------
      SUBROUTINE build_inputs(ncfile,slayer_inpr,growthrate_flag,
     $               qval_arr,inQ_arr,inQ_e_arr,inQ_i_arr,inc_beta_arr,
     $               inds_arr,intau_arr,Q0_arr,inpr_arr,inpe_arr,
     $               omegas_arr,outer_delta_arr)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: ncfile
      LOGICAL, INTENT(IN) :: growthrate_flag
      REAL(r8), INTENT(IN) ::slayer_inpr

      REAL(r8) :: respsi,lpsi,rpsi,hdist,sbnosurf,
     $ ising

      INTEGER :: zi, zimp, mi, mimp, msing
      REAL(r8) :: nfac,tfac,wefac,wpfac,e,twopi

      TYPE(spline_type) :: spl

      TYPE(spline_type) :: sr

      INTEGER :: mms,nns,mrs,nrs

      REAL(r8) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     qval,sval,bt,rs,zeff,inpe,r0,bt0
      REAL(r8) :: mu_i,tau_i,b_l,v_a,tau_r,tau_h,
     $            rho,tau_v,inpr,Qconv,lbeta

      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: qval_arr,
     $          inQ_arr,inQ_e_arr,
     $          inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,Q0_arr,
     $          inpr_arr,inpe_arr,omegas_arr,outer_delta_arr

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dp_diagonal, q_rational,
     $                                  psi_n_rational, shear
      INTEGER :: r_dim

c-----------------------------------------------------------------------
c     Read in STRIDE netcdf
c-----------------------------------------------------------------------
      character(512) ::
     $  kinetic_file = '/fusion/projects/codes/gpec/GPEC-1.5/
     $                   docs/examples/a10_ideal_example/a10_prof1.txt'
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: q_rational, q_rational_coords
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: Delta_prime, Delta_prime_coords
      WRITE(*,*)"ncfile=",ncfile
      CALL read_stride_netcdf_diagonal(ncfile,
     $              r_dim, dp_diagonal, q_rational, psi_n_rational,
     $              shear, r0, bt0)

      !CALL READ_Q_RATIONAL_AND_DELTA_PRIME_FROM_NETCDF(
      !&      filename, q_rational, q_rational_coords,
      !&      Delta_prime, Delta_prime_coords)

      msing = r_dim

      ALLOCATE(qval_arr(msing),inQ_arr(msing),inQ_e_arr(msing),
     $         inQ_i_arr(msing),
     $         inc_beta_arr(msing),inds_arr(msing),intau_arr(msing),
     $         Q0_arr(msing),outer_delta_arr(msing))
      !DEALLOCATE(q_rational, q_rational_coords, Delta_prime, Delta_prime_coords)
c-----------------------------------------------------------------------
c     Read and build equilibrium inputs
c-----------------------------------------------------------------------

    !=======================================================================
    !subroutine read_kin(file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,write_log)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Read ascii file containing table of kinetic profiles ni, ne, ti,
    !   te, and omegaE, then form kin spline containing some additional
    !   information (krook nui,nue).
    !
    !   Assumes table consists of 6 columns: psi_n, n_i(m^-3), n_e(m^-3),
    !   T_i(eV), T_e(eV), omega_E(rad/s). File can (nearly) arbitrary
    !   header and/or footer, with the exception that no lines start with
    !   a number.
    !
    !*ARGUMENTS:
    !    file : character(256) (in)
    !       File path.
    !   zi : integer
    !       Ion charge in fundamental units
    !   zimp : integer
    !       Impurity ion charge in fundamental units
    !   mi : integer
    !       Ion mass in fundamental units (mass proton)
    !   mimp : integer
    !       Impurity ion mass in fundamental units
    !   wefac : real
    !       Direct multiplier for omegaE profiles
    !   wpfac : real
    !       Scaling of rotation profile, done via manipulation of omegaE
    !   write_log : bool
    !       Writes kinetic spline to log file
    !
    !-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     set up kin
c-----------------------------------------------------------------------
        ! manually set the kinetic profiles
      zi = 1
      zimp = 1
      mi = 1
      mimp = 1
      nfac = 1.0
      tfac = 1.0
      wefac = 1.0
      wpfac = 1.0
      e=1.6021917e-19
      twopi = 6.28318530718

      CALL read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,.false.)
            ! manually set the perturbed equilibrium displacements
            ! use false flat xi and xi' for equal weighting

c-----------------------------------------------------------------------
c     loop across singular surfaces, evaluate spline quantities.
c-----------------------------------------------------------------------
      ! j_c is j_c/(chi1*sq%f(4))
      DO ising=1,msing
         !resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         !respsi=singtype(ising)%psifac
         !CALL spline_eval(sq,respsi,1)
         respsi = psi_n_rational(ising)
c-----------------------------------------------------------------------
c     prepare layer analysis.
c-----------------------------------------------------------------------
         !resm = mfac(resnum(ising))
         CALL spline_eval(sr,respsi,1)
         CALL spline_eval(kin,respsi,1)
c-----------------------------------------------------------------------
c     SLAYER inputs for sing surface
c-----------------------------------------------------------------------
         omega_i=-twopi*kin%f(3)*kin%f1(1)/(e*zi*chi1*kin%f(1))
     $           -twopi*kin%f1(3)/(e*zi*chi1)
         omega_e=twopi*kin%f(4)*kin%f1(2)/(e*chi1*kin%f(2))
     $           +twopi*kin%f1(4)/(e*chi1)

      ! Here's where I'm getting these from
         !CALL gpec_slayer(kin%f(2),kin%f(4)/e,kin%f(1),kin%f(3)/e,
      !$           kin%f(5),kin%f(9),omega_e,omega_i,sq%f(4),sq%f1(4),
      !$           bt0,sr%f1(1),ro,mi,slayer_inpr,resm,nn,ascii_flag,
      !$           delta_s,psi0,jxb,omega_sol,br_th)

         n_e = kin%f(2)
         t_e = kin%f(4)/e
         n_i = kin%f(1)
         t_i = kin%f(3)/e
         zeff = kin%f(5)
         omega = kin%f(9)
         qval = q_rational(ising)!sq%f(4)
         sval = shear(ising)
         bt = bt0
         rs = sr%f1(1)
         R0 = r0
         mu_i = mi
         inpr = slayer_inpr
         !mms = resm
         !nns = nn
         mrs = 2.0 ! FROM NAMELIST ???? real(mms,4)
         nrs = 1.0 ! FROM NAMELIST ???? real(nns,4)

          ! String representations of the m and n mode numbers
          !IF (nns<10) THEN
          !   WRITE(UNIT=sn,FMT='(I1)') nns
          !   sn=ADJUSTL(sn)
          !ELSE
          !   WRITE(UNIT=sn,FMT='(I2)') nns
          !ENDIF
          !IF (mms<10) THEN
          !   WRITE(UNIT=sm,FMT='(I1)') mms
          !   sm=ADJUSTL(sm)
          !ELSEIF (mms<100) THEN
          !   WRITE(UNIT=sm,FMT='(I2)') mms
          !   sm=ADJUSTL(sm)
          !ELSE
          !   WRITE(UNIT=sm,FMT='(I3)') mms
          !ENDIF

         inpe=0.0                         ! Waybright added this

         tau= t_i/t_e                     ! ratio of ion to electron temperature
         tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
         eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
         rho=(mu_i*m_p)*n_e               ! mass density

         b_l=(nrs/mrs)*nrs*sval*bt/R0     ! characteristic magnetic field
         v_a=b_l/(mu0*rho)**0.5           ! alfven velocity
         rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Lamour by elec. Temp.

         tau_h=R0*(mu0*rho)**0.5/(nns*sval*bt) ! alfven time across surface
         tau_r=mu0*rs**2.0/eta            ! resistive time scale
         tau_v=tau_r/inpr                   ! rho*rs**2.0/visc ! viscous time scale

          ! this one must be anomalous. calculated back from pr.
         visc= rho*rs**2.0/tau_v

         lu=tau_r/tau_h                   ! Lundquist number

         Qconv=lu**(1.0/3.0)*tau_h        ! conversion to Qs based on Cole

          ! note Q depends on Qconv even if omega is fixed.
         Q=Qconv*omega
         Q_e=-Qconv*omega_e
         Q_i=-Qconv*omega_i

          ! This is the most critical parameter
         ds=lu**(1.0/3.0)*rho_s/rs        ! conversion based on Cole.

         lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
         c_beta=(lbeta/(1.0+lbeta))**0.5

         delta_n=lu**(1.0/3.0)/rs         ! norm factor for delta primes

         qval_arr(ising) = qval
         inQ_arr(ising)=Q
         inQ_e_arr(ising)=Q_e
         inQ_i_arr(ising)=Q_i
         inc_beta_arr(ising)=c_beta
         inds_arr(ising)=ds
         intau_arr(ising)=tau
         Q0_arr(ising)=Q
         inpr_arr(ising) = inpr
         inpe_arr(ising) = 0.0 !!! TEMPORARY?
         omegas_arr(ising) = omega
         outer_delta_arr(ising) = dp_diagonal(ising)

      ENDDO
      CALL spline_dealloc(sr)
      !CALL cspline_dealloc(fsp_sol)
      !CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------

      RETURN
      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod