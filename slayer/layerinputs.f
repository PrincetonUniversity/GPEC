      MODULE layerinputs_mod

      USE inputs, ONLY : read_kin,read_equil,kin,chi1
      !USE params
      USE dcon_interface, ONLY : set_geom,geom
      !USE direct_mod, ONLY : direct_run
      USE spline_mod, ONLY : spline_alloc,spline_eval,spline_type,
     $                       spline_dealloc,spline_int,spline_fit
      USE sglobal_mod, ONLY: m_p, chag, lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,mu0,r8, ! NOT out_unit
     $   eta,visc,rho_s,lu,omega_e,omega_i,delta_n,Q
      USE netcdf
      USE pentrc_interface, ONLY : pentrc_timer=>timer
      !USE equil_mod
      !USE equil_out_mod

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
      SUBROUTINE read_stride_netcdf_diagonal(ncfile, msing,
     $   dp_diagonal, q_rational, psi_n_rational, shear,
     $   ro,bt0,psio,mpsi,nn,resm)

        !USE netcdf   ! NetCDF module for Fortran
        !USE stride_netcdf_mod ! For the 'check' subroutine (error handling)

        ! Input/Output Arguments
      CHARACTER(512), INTENT(IN) :: ncfile
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: dp_diagonal
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: q_rational,
     $  psi_n_rational, shear
      REAL(r8), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ro,bt0,psio,
     $            mpsi
      INTEGER, DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: nn,resm
      INTEGER, INTENT(OUT) :: msing

      REAL(r8), DIMENSION(:), ALLOCATABLE :: msing_arr

        ! Internal Variables
      INTEGER(kind=nf90_int) :: ncid, stat, r_dim_id, r_dim,
     $  dp_id, qr_id,pr_id,shear_id,ro_id,bt0_id,psio_id,mpsi_id,
     $  msing_id,nn_id,resm_id  ! Explicit kind for NetCDF variables
      INTEGER(kind=nf90_int), DIMENSION(1) :: start, count ! Explicit kind for NetCDF variables
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: delta_prime
      INTEGER :: i
      INTEGER :: bt0_len,ro_len,psio_len,mpsi_len,
     $             msing_len,nn_len    ! Attribute lengths


        ! Open the NetCDF file
      stat = nf90_open(path=ncfile,mode=NF90_WRITE,ncid=ncid)
      WRITE(*,*)"ncfile=",ncfile
      CALL check(stat)  ! Error handling

      stat = nf90_inquire_attribute(ncid,msing_id,"msing",
     $        len = msing_len)
      CALL check(stat)
      ALLOCATE(msing_arr(msing_len))
      stat = nf90_get_att(ncid,msing_id,"msing",msing_arr)
      CALL check(stat)

      msing=INT(msing_arr(1))

        ! Allocate Arrays (based on dimension)
      ALLOCATE(dp_diagonal(msing),q_rational(msing),
     $           psi_n_rational(msing),shear(msing),
     $           resm(msing))
      ALLOCATE(delta_prime(msing, msing,2))

      stat = nf90_inquire_attribute(ncid,ro_id,"ro",len = ro_len)
      CALL check(stat)
      stat = nf90_inquire_attribute(ncid,bt0_id,"bt0",len=bt0_len)
      CALL check(stat)

      bt0_id=0 !!!!! THIS COULD BE A PROBLEM

      stat = nf90_inquire_attribute(ncid,psio_id,"psio",len = psio_len)
      CALL check(stat)

      stat = nf90_inquire_attribute(ncid,mpsi_id,"mpsi",len = mpsi_len)
      CALL check(stat)
      stat = nf90_inquire_attribute(ncid,nn_id,"n",len = nn_len)
      CALL check(stat)

      ALLOCATE(bt0(INT(bt0_len)),ro(INT(ro_len)),psio(INT(psio_len)),
     $mpsi(INT(mpsi_len)),nn(INT(nn_len)))

        ! Get Variable IDs
      stat = nf90_inq_varid(ncid, "Delta_prime", dp_id)
      CALL check(stat)
      stat = nf90_inq_varid(ncid, "q_rational", qr_id)
      CALL check(stat)
      stat = nf90_inq_varid(ncid, "psi_n_rational", pr_id)
      CALL check(stat)
      stat = nf90_inq_varid(ncid, "shear", shear_id)
      CALL check(stat)
      stat = nf90_inq_varid(ncid, "resm", resm_id)
      CALL check(stat)

        ! Read Data from NetCDF File
        ! Set up start and count for reading only the diagonal
        !start(1) = 1
        !count(1) = 1

      ! Get attributes
      stat = nf90_get_att(ncid, ro_id, "ro", ro)
      CALL check(stat)

      stat = nf90_get_att(ncid, bt0_id, "bt0", bt0)
      CALL check(stat)
      stat = nf90_get_att(ncid, psio_id, "psio", psio)
      CALL check(stat)
      stat = nf90_get_att(ncid, mpsi_id, "mpsi", mpsi)
      CALL check(stat)
      stat = nf90_get_att(ncid, nn_id, "n", nn)
      CALL check(stat)

        ! Read the diagonal of delta_prime. The results will be put on a 1D temporary array.
      stat = nf90_get_var(ncid, dp_id, delta_prime,start=(/ 1,1,1 /))
      CALL check(stat)
      ! Read 1D variables
      stat = nf90_get_var(ncid, qr_id, q_rational)
      CALL check(stat)
      stat = nf90_get_var(ncid, pr_id, psi_n_rational)
      CALL check(stat)
      stat = nf90_get_var(ncid, shear_id, shear)
      CALL check(stat)
      stat = nf90_get_var(ncid, resm_id, resm)
      CALL check(stat)

      ! Extract Diagonal, with 3rd index signifying REAL part
      DO i = 1, msing
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
      SUBROUTINE build_inputs(infile,ncfile,slayer_inpr,
     $               growthrate_flag,
     $               qval_arr,psi_n_rational,inQ_arr,inQ_e_arr,
     $               inQ_i_arr,inc_beta_arr,
     $               inds_arr,intau_arr,Q0_arr,inpr_arr,inpe_arr,
     $               omegas_arr,outer_delta_arr)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      LOGICAL, INTENT(IN) :: growthrate_flag
      REAL(r8), INTENT(IN) ::slayer_inpr

      REAL(r8) :: respsi,lpsi,rpsi,hdist,sbnosurf,
     $ ising

      INTEGER :: zi, zimp, mi, mimp
      REAL(r8) :: nfac,tfac,wefac,wpfac,e,pi,twopi

      TYPE(spline_type) :: spl
      TYPE(spline_type) :: sr

      INTEGER :: mms,nns,mrs,nrs,mpsi

      REAL(r8) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     qval,sval,bt,rs,zeff,inpe,R0
      REAL(r8) :: mu_i,tau_i,b_l,v_a,tau_r,tau_h,
     $            rho,tau_v,inpr,Qconv,lbeta,qintb

      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: qval_arr,
     $          inQ_arr,inQ_e_arr,psi_n_rational,
     $          inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,Q0_arr,
     $          inpr_arr,inpe_arr,omegas_arr,outer_delta_arr
      REAL(r8), DIMENSION(0:128) :: psitor, rhotor
      REAL(r8), DIMENSION(:), ALLOCATABLE :: my_rhotor,my_psitor

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dp_diagonal, q_rational,
     $                      shear,ro,bt0,psio,mpsi_arr
      INTEGER,DIMENSION(:),ALLOCATABLE :: nn,resm
      INTEGER :: msing,i


c-----------------------------------------------------------------------
c     Read in STRIDE netcdf
c-----------------------------------------------------------------------
      CHARACTER(512), INTENT(IN) :: infile,ncfile

      !character(512) ::
      !$ kinetic_file = '/fusion/projects/codes/gpec/users/burgessd/GPEC/
      !$bin/g147131.02300.txt'
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: q_rational, q_rational_coords
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: Delta_prime, Delta_prime_coords

      CALL read_stride_netcdf_diagonal(ncfile,
     $              msing, dp_diagonal, q_rational, psi_n_rational,
     $              shear, ro, bt0, psio, mpsi_arr,nn,resm)
      WRITE(*,*)"msing_out=",msing
      WRITE(*,*)"dp_diagonal=",dp_diagonal
      WRITE(*,*)"q_rational=",q_rational
      WRITE(*,*)"psi_n_rational=",psi_n_rational
      WRITE(*,*)"shear=",shear
      WRITE(*,*)"ro=",ro
      WRITE(*,*)"bt0=",bt0
      WRITE(*,*)"psio=",psio
      WRITE(*,*)"nn=",nn
      WRITE(*,*)"resm=",resm

      mpsi = INT(mpsi_arr(1))
      WRITE(*,*)"mpsi=",mpsi

      !CALL READ_Q_RATIONAL_AND_DELTA_PRIME_FROM_NETCDF(
      !&      filename, q_rational, q_rational_coords,
      !&      Delta_prime, Delta_prime_coords)
      ALLOCATE(qval_arr(msing),inQ_arr(msing),inQ_e_arr(msing),
     $         inQ_i_arr(msing),
     $         inc_beta_arr(msing),inds_arr(msing),intau_arr(msing),
     $         Q0_arr(msing),inpr_arr(msing),inpe_arr(msing),
     $         omegas_arr(msing),outer_delta_arr(msing))
      !ALLOCATE(my_rhotor(mpsi),my_psitor(1))
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
      zimp = 6
      mi = 2
      mimp = 12
      nfac = 1.0
      tfac = 1.0
      wefac = 1.0
      wpfac = 1.0
      e=1.6021917e-19
      pi = 3.1415926535897932
      twopi = 2*pi
      chi1 = twopi*psio(1)

      CALL read_kin(infile,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,.false.)
            ! manually set the perturbed equilibrium displacements
            ! use false flat xi and xi' for equal weighting

      !CALL spline_int(sq)
      !WRITE(*,*)"sq=",sq
      !WRITE(*,*)"mpsi=",mpsi
      !CALL spline_alloc(sq,mpsi,4)

      !CALL spline_int(sq)

      !CALL timer(0,out_unit)
      !CALL equil_read(out_unit)
      !CALL equil_out_global
      !CALL equil_out_qfind

      !if(.not. sq%allocated)
      !$      stop 'ERROR: Cannot define geometric splines without sq'


      !CALL set_geom
      !WRITE(*,*)"geom%f(2)=",geom%fs(:,2)




      ! minor radius defined using toroidal flux. Used for threshold
      !CALL spline_int(sq)
      !qintb = sq%fsi(mpsi, 4)
      !psitor(:) = sq%fsi(:, 4) / qintb  ! normalized toroidal flux
      !rhotor(:) = SQRT(sq%fsi(:, 4)*twopi*psio / (pi * bt0))  ! effective minor radius in Callen
      !CALL spline_alloc(sr,mpsi,1)
      !sr%xs = sq%xs
      !sr%fs(:, 1) = rhotor(:)
      !CALL spline_fit(sr,"extrap")


      !WRITE(*,*)"twopi*psio / (pi * bt0)=",twopi*psio / (pi * bt0)
      !WRITE(*,*)"sq%fsi(:, 4)=",sq%fsi(:, 4)
      !WRITE(*,*)"sq%fsi(:, 4)*twopi*psio / (pi * bt0)="
      !$(twopi*psio / (pi * bt0))*
      !DO i=1,mpsi
      !  WRITE(*,*)"spline(i)=",sq%fsi(i, 4)
      !  WRITE(*,*)"(i)=",SIZE(SQRT(sq%fsi(i, 4)*(twopi*psio /
      !$(pi * bt0))))
      !  my_psitor = REAL(sq%fsi(i, 4))
        !WRITE(*,*)"my_psitor(1)=",my_psitor(1)
        !WRITE(*,*)"(my_psitor(1))=",my_psitor(1)
        !WRITE(*,*)"(my_rhotor(i))=",my_rhotor(i)
        !my_rhotor(i)=SQRT(my_psitor(1)*(twopi*psio /
      !$(pi * bt0)))

      !ENDDO
      !WRITE(*,*)"rhotor=",rhotor

c-----------------------------------------------------------------------
c     loop across singular surfaces, evaluate spline quantities.
c-----------------------------------------------------------------------
      ! j_c is j_c/(chi1*sq%f(4))
      DO ising=1,msing
         !resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         !respsi=singtype(ising)%psifac
         !CALL spline_eval(sq,respsi,1)
         respsi = psi_n_rational(ising)
         WRITE(*,*)"respsi=",respsi
         WRITE(*,*)"chi1=",chi1

c-----------------------------------------------------------------------
c     prepare layer analysis.
c-----------------------------------------------------------------------
         !resm = mfac(resnum(ising))
         !CALL spline_eval(sr,respsi,1)
         CALL spline_eval(kin,respsi,1)
         !CALL spline_eval(geom,respsi,1)
c-----------------------------------------------------------------------
c     SLAYER inputs for sing surface
c-----------------------------------------------------------------------



         omega_i=-twopi*kin%f(3)*kin%f1(1)/(e*zi*chi1*kin%f(1))
     $           -twopi*kin%f1(3)/(e*zi*chi1)
         omega_e=twopi*kin%f(4)*kin%f1(2)/(e*chi1*kin%f(2))
     $           +twopi*kin%f1(4)/(e*chi1)
         WRITE(*,*)"omega_i=",omega_i
         WRITE(*,*)"omega_e=",omega_e

      ! Here's where I'm getting these from
         !CALL gpec_slayer(kin%f(2),kin%f(4)/e,kin%f(1),kin%f(3)/e,
      !$           kin%f(5),kin%f(9),omega_e,omega_i,sq%f(4),sq%f1(4),
      !$           bt0,sr%f1(1),ro,mi,slayer_inpr,resm,nn,ascii_flag,
      !$           delta_s,psi0,jxb,omega_sol,br_th)

         !singtype(ising)%q = q(ising)?
         !resnum=NINT(singtype(ising)%q*nn)-mlow+1
         !shear(ising)=mfac(resnum)*sq%f1(4)/sq%f(4)**2

         n_e = kin%f(2)
         t_e = kin%f(4)/e
         n_i = kin%f(1)
         t_i = kin%f(3)/e
         zeff = kin%f(9)
         omega = kin%f(5)
         qval = q_rational(ising)!sq%f(4)
         sval = shear(ising) ! SHEAR SO MUCH SMALLER THAN BEFORE???? 500 puts it on the root
         bt = bt0(1)
         rs = 0.167
         R0 = ro(1)
         mu_i = mi
         inpr = slayer_inpr

         !!! TEMPORARY
         !mm=2
         !nn=1
         !mr = real(mm,4)
         !nr = real(nn,4)
         mms = resm(ising)
         nns = nn(1)
         mrs = real(mms,4)
         nrs = real(nns,4)

         WRITE(*,*)"nns=",nns
         WRITE(*,*)"mms=",mms

         WRITE(*,*)"n_e=",n_e
         WRITE(*,*)"t_e=",t_e
         WRITE(*,*)"n_i=",n_i
         WRITE(*,*)"t_i=",t_i
         WRITE(*,*)"zeff=",zeff
         WRITE(*,*)"omega=",omega
         WRITE(*,*)"qval=",qval
         WRITE(*,*)"sval=",sval
         WRITE(*,*)"bt=",bt
         WRITE(*,*)"rs=",rs
         WRITE(*,*)"R0=",R0
         WRITE(*,*)"mu_i=",mu_i
         WRITE(*,*)"inpr=",inpr
         WRITE(*,*)"dp_diagonal(ising)=",dp_diagonal
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

         WRITE(*,*)"lnLamb=",lnLamb
         WRITE(*,*)"mu_i=",mu_i
         WRITE(*,*)"m_p=",m_p

         tau= t_i/t_e                     ! ratio of ion to electron temperature
         tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
         eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
         rho=(mu_i*m_p)*n_e               ! mass density




         b_l=(nrs/mrs)*nrs*sval*bt/R0     ! characteristic magnetic field
         WRITE(*,*)"b_l=",b_l

         v_a=b_l/(mu0*rho)**0.5           ! alfven velocity, B_L IS BROKEN
         rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Lamour by elec. Temp.

         tau_h=R0*(mu0*rho)**0.5/(nns*sval*bt) ! alfven time across surface
         tau_r=mu0*rs**2.0/eta            ! resistive time scale
         tau_v=tau_r/inpr                   ! rho*rs**2.0/visc ! viscous time scale
         WRITE(*,*)"tau_v=",tau_v

          ! this one must be anomalous. calculated back from pr.
         visc= rho*rs**2.0/tau_v

         lu=tau_r/tau_h                   ! Lundquist number
         WRITE(*,*)"lu=",lu

         Qconv=lu**(1.0/3.0)*tau_h        ! conversion to Qs based on Cole

          ! note Q depends on Qconv even if omega is fixed.
         Q=Qconv*omega
         Q_e=-Qconv*omega_e
         Q_i=-Qconv*omega_i
         WRITE(*,*)"Q=",Q

          ! This is the most critical parameter
         ds=lu**(1.0/3.0)*rho_s/rs        ! conversion based on Cole.

         lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
         c_beta=(lbeta/(1.0+lbeta))**0.5
         WRITE(*,*)"c_beta=",c_beta

         delta_n=lu**(1.0/3.0)/rs         ! norm factor for delta primes
         WRITE(*,*)"delta_n=",delta_n

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

      !WRITE(*,*)"qval_arr=",qval_arr
      !WRITE(*,*)"inQ_arr=",inQ_arr
      !WRITE(*,*)"inQ_e_arr=",inQ_e_arr
      !WRITE(*,*)"inQ_i_arr=",inQ_i_arr
      !WRITE(*,*)"inc_beta_arr=",inc_beta_arr
      !WRITE(*,*)"inds_arr=",inds_arr
      !WRITE(*,*)"intau_arr=",intau_arr
      !WRITE(*,*)"Q0_arr=",Q0_arr
      !!WRITE(*,*)"inpr_arr=",inpr_arr
      !WRITE(*,*)"inpe_arr=",inpe_arr
      !WRITE(*,*)"omegas_arr=",omegas_arr
      !WRITE(*,*)"outer_delta_arr=",outer_delta_arr

      !CALL spline_dealloc(sr)
      !CALL cspline_dealloc(fsp_sol)
      !CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------

      RETURN
      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod