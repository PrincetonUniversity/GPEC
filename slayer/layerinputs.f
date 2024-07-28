      MODULE layerinputs_mod

      USE inputs, ONLY : read_kin,read_equil,kin,chi1
      USE spline_mod, ONLY : spline_alloc,spline_eval,spline_type,
     $                       spline_dealloc,spline_int,spline_fit
      USE sglobal_mod, ONLY: m_p, chag, lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,r8,mu0,pi,out_unit, ! NOT out_unit
     $   eta,visc,rho_s,lu,omega_e,omega_i,delta_n,Q
      USE netcdf
      USE equil_mod, ONLY: equil_read,rzphi,twopi,ro,zo,sq
      USE bicube_mod, ONLY: bicube_eval_external,bicube_type
      USE slayer_netcdf_mod, ONLY: slayer_netcdf_inputs

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
c     subprogram 2. read_stride_netcdf_diagonal.
c     Read STRIDE netcdf file for SLAYER inputs only.
c-----------------------------------------------------------------------
      SUBROUTINE read_stride_netcdf_diagonal(ncfile, msing,
     $   dp_diagonal, q_rational, psi_n_rational, shear,
     $   r_o,my_bt0,my_psio,mpsi,nn,resm,prandtl)

        ! Input/Output Arguments
      CHARACTER(512), INTENT(IN) :: ncfile
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: dp_diagonal
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: q_rational,
     $  psi_n_rational, shear, prandtl
      REAL(r8), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: r_o,my_bt0,
     $ my_psio,mpsi
      INTEGER, DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: nn,resm
      INTEGER, INTENT(OUT) :: msing

      REAL(r8), DIMENSION(:), ALLOCATABLE :: msing_arr

        ! Internal Variables
      INTEGER(kind=nf90_int) :: ncid, stat, r_dim_id, r_dim,
     $  dp_id, qr_id,pr_id,shear_id,ro_id,bt0_id,psio_id,mpsi_id,
     $  msing_id,nn_id,resm_id,prandtl_id  ! Explicit kind for NetCDF variables
      INTEGER(kind=nf90_int), DIMENSION(1) :: start, count ! Explicit kind for NetCDF variables
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: delta_prime
      INTEGER :: i
      INTEGER :: bt0_len,ro_len,psio_len,mpsi_len,
     $             msing_len,nn_len    ! Attribute lengths

        ! Open the NetCDF file
      stat = nf90_open(path=ncfile,mode=NF90_WRITE,ncid=ncid)
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
     $           resm(msing),prandtl(msing))
      ALLOCATE(delta_prime(msing, msing,2))

      stat = nf90_inquire_attribute(ncid,ro_id,"ro",len = ro_len)
      CALL check(stat)
      stat = nf90_inquire_attribute(ncid,bt0_id,"bt0",len=bt0_len)
      CALL check(stat)


      stat = nf90_inquire_attribute(ncid,psio_id,"psio",len=psio_len)
      CALL check(stat)

      stat = nf90_inquire_attribute(ncid,mpsi_id,"mpsi",len=mpsi_len)
      CALL check(stat)
      stat = nf90_inquire_attribute(ncid,nn_id,"n",len = nn_len)
      CALL check(stat)

      bt0_id=0 !!!!! THIS COULD BE A PROBLEM
      nn_id=0
      mpsi_id=0
      psio_id=0
      ro_id=0

      ALLOCATE(my_bt0(INT(bt0_len)),r_o(INT(ro_len)),
     $ my_psio(INT(psio_len)),
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
      stat = nf90_inq_varid(ncid, "prandtl", prandtl_id)
      CALL check(stat)
      stat = nf90_inq_varid(ncid, "resm", resm_id)
      CALL check(stat)
      ! Get attributes
      stat = nf90_get_att(ncid, ro_id, "ro", r_o)
      CALL check(stat)

      stat = nf90_get_att(ncid, bt0_id, "bt0", my_bt0)
      CALL check(stat)
      stat = nf90_get_att(ncid, psio_id, "psio", my_psio)
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
      stat = nf90_get_var(ncid, prandtl_id, prandtl)
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
c     subprogram 2. issurfint.
c     surface integration by simple method.
c-----------------------------------------------------------------------
      FUNCTION issurfint(func,fs,inpsi,wegt,ave,
     $     fsave,psave,jacs,delpsi,inr,ina,first)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      !IMPLICIT NONE
      INTEGER, INTENT(IN) :: fs,wegt,ave
      REAL(r8), INTENT(IN) :: inpsi
      REAL(r8), DIMENSION(0:fs), INTENT(IN) :: func

      LOGICAL, INTENT(INOUT) :: first
      INTEGER, INTENT(INOUT)  :: fsave
      REAL(r8), INTENT(INOUT) :: psave
      REAL(r8),DIMENSION(0:),INTENT(INOUT) :: jacs,delpsi,inr,ina
      INTEGER  :: itheta, ix, iy
      REAL(r8) :: issurfint
      REAL(r8) :: rfac,ineta,injac,inarea
      REAL(r8), DIMENSION(1,2) :: w
      REAL(r8), DIMENSION(0:fs) :: z,thetas
      REAL(r8), dimension(4) :: rzphi_f, rzphi_fx, rzphi_fy

      issurfint=0
      inarea=0
      ix = 0
      iy = 0
      IF(first .OR. inpsi/=psave .OR. fs/=fsave)THEN
         psave = inpsi
         fsave = fs
         !first = .FALSE.
         DO itheta=0,fs
            thetas(itheta) = REAL(itheta,r8)/REAL(fs,r8)
         ENDDO
         DO itheta=0,fs-1
            CALL bicube_eval_external(rzphi, inpsi, thetas(itheta), 1,
     $           ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
            rfac=SQRT(rzphi_f(1))
            ineta=twopi*(thetas(itheta)+rzphi_f(2))
            ina(itheta)=rfac
            inr(itheta)=ro+rfac*COS(ineta)
            z(itheta)=zo+rfac*SIN(ineta)
            injac=rzphi_f(4)
            jacs(itheta)=injac
            w(1,1)=(1+rzphi_fy(2))*twopi**2*rfac*inr(itheta)/injac
            w(1,2)=-rzphi_fy(1)*pi*inr(itheta)/(rfac*injac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         ENDDO
      ENDIF

      IF (wegt==0) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==1) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $         inr(itheta)*jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==2) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/inr(itheta)/fs
         ENDDO
      ELSE IF (wegt==3) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $        ina(itheta)*jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE
         STOP "ERROR: issurfint wegt must be in [0,1,2,3]"
      ENDIF

      IF (ave==1) THEN
         DO itheta=0,fs-1
            inarea=inarea+jacs(itheta)*delpsi(itheta)/fs
         ENDDO
         issurfint=issurfint/inarea
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfint
c-----------------------------------------------------------------------
c     subprogram 3. build_inputs.
c     build input arrays for SLAYER
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
      LOGICAL :: firstsurf
      REAL(r8) :: respsi,lpsi,rpsi,hdist,sbnosurf,
     $ ising
      INTEGER :: zi, zimp, mi, mimp
      REAL(r8) :: nfac,tfac,wefac,wpfac,e!,pi,twopi

      TYPE(spline_type) :: spl
      TYPE(spline_type) :: sr

      INTEGER :: mms,nns,mrs,nrs,mpsi

      REAL(r8) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     my_qval,my_sval,my_bt,my_rs,zeff,inpe,R_0
      REAL(r8) :: mu_i,tau_i,b_l,v_a,tau_r,tau_h,
     $            rho,tau_v,inpr,Qconv,lbeta,qintb

      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $          inQ_arr,inQ_e_arr,psi_n_rational,
     $          inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,Q0_arr,
     $          inpr_arr,inpe_arr,omegas_arr,
     $          outer_delta_arr
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: qval_arr
      REAL(r8), DIMENSION(0:128) :: psitor, rhotor
      REAL(r8), DIMENSION(:), ALLOCATABLE :: my_rhotor,my_psitor

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dp_diagonal, q_rational,
     $                      shear,r_o,my_bt0,my_psio,mpsi_arr,
     $                      omegas_e_arr,omegas_i_arr,prandtl
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ne_arr,te_arr,ni_arr,
     $    ti_arr,zeff_arr,bt_arr,rs_arr,
     $    R0_arr,mu_i_arr
      INTEGER,DIMENSION(:),ALLOCATABLE :: nn,resm,nns_arr
      INTEGER :: msing,i,mthsurf
      REAL(r8), DIMENSION(0:512) :: unitfun
      INTEGER :: fsave
      REAL(r8) :: psave
      REAL(r8), DIMENSION(:), ALLOCATABLE :: jacs,delpsi,rsurf,asurf
      REAL(r8) :: rfac,jac,a_surf

      CHARACTER(512), INTENT(IN) :: infile,ncfile

c-----------------------------------------------------------------------
c     Read in STRIDE netcdf
c-----------------------------------------------------------------------
      CALL read_stride_netcdf_diagonal(ncfile,
     $           msing,dp_diagonal,q_rational,psi_n_rational,
     $           shear,r_o,my_bt0,my_psio,mpsi_arr,nn,resm,prandtl)
      WRITE(*,*)"msing_out=",msing
      WRITE(*,*)"dp_diagonal=",dp_diagonal
      WRITE(*,*)"q_rational=",q_rational
      WRITE(*,*)"psi_n_rational=",psi_n_rational
      WRITE(*,*)"shear=",shear
      WRITE(*,*)"r_o=",r_o
      WRITE(*,*)"my_bt0=",my_bt0
      WRITE(*,*)"my_psio=",my_psio
      WRITE(*,*)"nn=",nn
      WRITE(*,*)"resm=",resm
      WRITE(*,*)"resm=",prandtl

      mpsi = INT(mpsi_arr(1))
      mthsurf = 512

      ALLOCATE(qval_arr(msing),inQ_arr(msing),inQ_e_arr(msing),
     $     inQ_i_arr(msing),
     $     inc_beta_arr(msing),inds_arr(msing),intau_arr(msing),
     $     Q0_arr(msing),inpr_arr(msing),inpe_arr(msing),
     $     omegas_arr(msing),omegas_e_arr(msing),omegas_i_arr(msing),
     $     outer_delta_arr(msing))
      ALLOCATE(ne_arr(msing),te_arr(msing),ni_arr(msing),
     $    ti_arr(msing),zeff_arr(msing),bt_arr(msing),rs_arr(msing),
     $    R0_arr(msing),mu_i_arr(msing),nns_arr(msing))

      ALLOCATE(jacs(0:mthsurf),delpsi(0:mthsurf),
     $                 rsurf(0:mthsurf),asurf(0:mthsurf))
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
      chi1 = twopi*my_psio(1)

      CALL read_kin(infile,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,.false.)

      CALL equil_read(out_unit)

c-----------------------------------------------------------------------
c     loop across singular surfaces, evaluate spline quantities.
c-----------------------------------------------------------------------

      DO ising=1,msing

         respsi = psi_n_rational(ising)

         firstsurf = .TRUE.
         unitfun = 1

         ! Minor radius!
         a_surf = issurfint(unitfun,mthsurf,respsi,3,1,
     $           fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)

c-----------------------------------------------------------------------
c     SLAYER inputs for sing surface
c-----------------------------------------------------------------------
         CALL spline_eval(kin,respsi,1)

         omega_i=-twopi*kin%f(3)*kin%f1(1)/(e*zi*chi1*kin%f(1))
     $           -twopi*kin%f1(3)/(e*zi*chi1)
         omega_e=twopi*kin%f(4)*kin%f1(2)/(e*chi1*kin%f(2))
     $           +twopi*kin%f1(4)/(e*chi1)

         n_e = kin%f(2)
         t_e = kin%f(4)/e
         n_i = kin%f(1)
         t_i = kin%f(3)/e
         zeff = kin%f(9)
         omega = kin%f(5)
         my_qval = q_rational(ising)!sq%f(4)
         my_sval = shear(ising) ! SHEAR SO MUCH SMALLER THAN BEFORE???? 500 puts it on the root
         my_bt = my_bt0(1)
         my_rs = a_surf
         R_0 = r_o(1)
         mu_i = 2.0

         eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)

         inpr = prandtl(ising)

         inpe=0.0!0.0165*inpr                        ! Waybright added this

         ne_arr(ising) = n_e
         te_arr(ising) = t_e
         ni_arr(ising) = n_i
         ti_arr(ising) = t_i
         zeff_arr(ising) = zeff
         bt_arr(ising) = my_bt
         rs_arr(ising) = my_rs
         R0_arr(ising) = R_0
         mu_i_arr(ising) = mu_i

         mms = resm(ising)
         nns = nn(1)
         mrs = real(mms,4)
         nrs = real(nns,4)

         nns_arr = nn(1)

         tau= t_i/t_e                     ! ratio of ion to electron temperature
         tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
         rho=(mu_i*m_p)*n_e               ! mass density

         b_l=(nrs/mrs)*nrs*my_sval*my_bt/R_0     ! characteristic magnetic field

         v_a=b_l/(mu0*rho)**0.5           ! alfven velocity, B_L IS BROKEN
         rho_s=1.02e-4*(mu_i*t_e)**0.5/my_bt ! ion Lamour by elec. Temp.

         tau_h=R_0*(mu0*rho)**0.5/(nns*my_sval*my_bt) ! alfven time across surface
         tau_r=mu0*my_rs**2.0/eta            ! resistive time scale
         tau_v=tau_r/inpr                   ! rho*rs**2.0/visc ! viscous time scale

          ! this one must be anomalous. calculated back from pr.
         visc= rho*my_rs**2.0/tau_v

         lu=tau_r/tau_h                   ! Lundquist number

         Qconv=lu**(1.0/3.0)*tau_h        ! conversion to Qs based on Cole

          ! note Q depends on Qconv even if omega is fixed.
         Q=Qconv*omega
         Q_e=-Qconv*omega_e
         Q_i=-Qconv*omega_i

          ! This is the most critical parameter
         ds=lu**(1.0/3.0)*rho_s/my_rs        ! conversion based on Cole.

         lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/my_bt**2.0
         c_beta=(lbeta/(1.0+lbeta))**0.5

         delta_n=lu**(1.0/3.0)/my_rs         ! norm factor for delta primes

         qval_arr(ising) = INT(my_qval)
         inQ_arr(ising)=Q
         inQ_e_arr(ising)=Q_e
         inQ_i_arr(ising)=Q_i
         inc_beta_arr(ising)=c_beta
         inds_arr(ising)=ds
         intau_arr(ising)=tau
         Q0_arr(ising)=Q
         inpr_arr(ising) = inpr
         inpe_arr(ising) = inpe!0.0!0.0165*inpr !!! TEMPORARY?
         omegas_arr(ising) = omega
         omegas_e_arr(ising) = omega_e
         omegas_i_arr(ising) = omega_i
         outer_delta_arr(ising) = dp_diagonal(ising)

      ENDDO

      !WRITE(*,*)"zeff_arr=",zeff_arr

      CALL slayer_netcdf_inputs(msing,ne_arr,te_arr,ni_arr,ti_arr,
     $  zeff_arr,shear,bt_arr,rs_arr,R0_arr,mu_i_arr,resm,nns_arr,
     $  qval_arr,inQ_arr,inQ_e_arr,inQ_i_arr,inc_beta_arr,inds_arr,
     $  intau_arr,inpr_arr,inpe_arr,omegas_arr,
     $  outer_delta_arr)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------

      RETURN


      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod