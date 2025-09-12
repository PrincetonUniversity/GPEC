      MODULE layerinputs_mod

      USE inputs, ONLY : read_kin,read_equil,kin,chi1
      USE spline_mod, ONLY : spline_alloc,spline_eval,spline_type,
     $                       spline_dealloc,spline_int,spline_fit
      USE sglobal_mod
      USE params_mod
      USE netcdf
      USE equil_mod, ONLY: equil_read,rzphi,twopi,ro,zo,sq
      USE bicube_mod, ONLY: bicube_eval_external,bicube_type
      USE slayer_netcdf_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. read_stride_netcdf_diagonal.
c     Read STRIDE netcdf file for SLAYER inputs only.
c-----------------------------------------------------------------------
      SUBROUTINE read_stride_netcdf_diagonal(ncfile,msing,dp_mat,
     $   Re_dp_diagonal,Im_dp_diagonal,q_rational,psi_n_rational,dgeo,
     $   shear,r_o,my_bt0,my_psio,dr_vals,mpsi,nn,resm)

        ! Input/Output Arguments
      CHARACTER(512), INTENT(IN) :: ncfile
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $                                 Re_dp_diagonal,Im_dp_diagonal
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT):: dp_mat
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: q_rational,
     $                                     psi_n_rational, shear,dgeo
      REAL(r8), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: r_o,my_bt0,
     $                                         my_psio,mpsi,dr_vals
      INTEGER, DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: nn,resm
      INTEGER, INTENT(OUT) :: msing

      REAL(r8), DIMENSION(:), ALLOCATABLE :: msing_arr

        ! Internal Variables
      INTEGER(kind=nf90_int) :: ncid, stat, r_dim_id, r_dim,
     $  dp_id, qr_id,pr_id,dgeo_id,shear_id,ro_id,bt0_id,psio_id,
     $  mpsi_id,msing_id,nn_id,resm_id,drr_id ! Explicit kind for NetCDF variables
      INTEGER(kind=nf90_int), DIMENSION(1) :: start, count ! Explicit kind for NetCDF variables
      INTEGER :: i
      INTEGER :: bt0_len,ro_len,psio_len,mpsi_len,
     $             msing_len,nn_len,dr_len    ! Attribute lengths

      WRITE(*,*)"$^$ opening netcdf file",ncfile

        ! Open the NetCDF file
      stat = nf90_open(path=ncfile,mode=NF90_WRITE,ncid=ncid)
      CALL sl_check(stat)  ! Error handling

      stat = nf90_inquire_attribute(ncid,msing_id,"msing",
     $        len = msing_len)
      CALL sl_check(stat)
      ALLOCATE(msing_arr(msing_len))
      stat = nf90_get_att(ncid,msing_id,"msing",msing_arr)
      CALL sl_check(stat)

      msing=INT(msing_arr(1))

      ! Allocate Arrays (based on dimension)
      ALLOCATE(Re_dp_diagonal(msing),q_rational(msing),
     $           psi_n_rational(msing),shear(msing),dgeo(msing),
     $           resm(msing),Im_dp_diagonal(msing),dr_vals(msing))
      ALLOCATE(dp_mat(msing, msing,2))

      stat = nf90_inquire_attribute(ncid,ro_id,"ro",len = ro_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid,bt0_id,"bt0",len=bt0_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid,psio_id,"psio",len=psio_len)
      CALL sl_check(stat)

      stat = nf90_inquire_attribute(ncid,mpsi_id,"mpsi",len=mpsi_len)
      CALL sl_check(stat)
      stat = nf90_inquire_attribute(ncid,nn_id,"n",len = nn_len)
      CALL sl_check(stat)

      !bt0_id=0 !!!!! THIS COULD BE A PROBLEM
      !nn_id=0
      !mpsi_id=0
      !psio_id=0
      !ro_id=0

      ALLOCATE(my_bt0(INT(bt0_len)),r_o(INT(ro_len)),
     $         my_psio(INT(psio_len)),
     $         mpsi(INT(mpsi_len)),nn(INT(nn_len)))

      ! Get Variable IDs
      stat = nf90_inq_varid(ncid, "Delta_prime", dp_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "q_rational", qr_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "psi_n_rational", pr_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "Delta_geo", dgeo_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "shear", shear_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "resm", resm_id)
      CALL sl_check(stat)
      stat = nf90_inq_varid(ncid, "dr_rational", drr_id)
      CALL sl_check(stat)
      ! Get attributes
      stat = nf90_get_att(ncid, ro_id, "ro", r_o)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, bt0_id, "bt0", my_bt0)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, psio_id, "psio", my_psio)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, mpsi_id, "mpsi", mpsi)
      CALL sl_check(stat)
      stat = nf90_get_att(ncid, nn_id, "n", nn)
      CALL sl_check(stat)

      ! Read the diagonal of delta prime. The results will be put on a 1D temporary array.
      stat = nf90_get_var(ncid, dp_id, dp_mat,start=(/ 1,1,1 /))
      CALL sl_check(stat)
      ! Read 1D variables
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

      ! Extract Diagonal, with 3rd index signifying REAL part
      DO i = 1, msing
        Re_dp_diagonal(i) = dp_mat(i, i, 1)
        Im_dp_diagonal(i) = dp_mat(i, i, 2)
      END DO
      ! Clean Up
      stat = nf90_close(ncid)
      CALL sl_check(stat)

      END SUBROUTINE read_stride_netcdf_diagonal
c-----------------------------------------------------------------------
c     subprogram 2. issurfint.
c     surface integration by simple method. copied from EQUIL
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
      SUBROUTINE build_inputs(infile,ncfile,sl_in)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      ! Inputs
      CHARACTER(512), INTENT(IN) :: infile,ncfile
      TYPE(slayer_inputs_type), INTENT(INOUT) :: sl_in
      ! Internals
      LOGICAL :: firstsurf
      REAL(r8) :: respsi,lpsi,rpsi,hdist,sbnosurf,ising
      INTEGER :: zi, zimp, mi, mimp
      REAL(r8) :: nfac,tfac,wefac,wpfac,e

      TYPE(spline_type) :: spl
      TYPE(spline_type) :: sr

      INTEGER :: mms,nns,mrs,nrs,mpsi

      REAL(r8) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     my_qval,my_sval,my_bt,my_rs,my_inpe,zeff,R_0,dgeo_val
      REAL(r8) :: mu_i,tau_i,b_l,v_a,tau_h,l_n,l_t,
     $            rho,tau_v,Qconv,lbeta,qintb,gammafac
      REAL(r8) :: tau_ee_num,tau_ee_denom,tau_ee,sigma_par_1,
     $            sigma_par_2,sigma_par,tau_perp,Wd,vte,
     $            dr_val,chi_par_smfp,chi_par_lmfp,chi_par
      REAL(r8), DIMENSION(3) :: chi_s
      INTEGER :: wit

      REAL(r8), DIMENSION(0:128) :: psitor, rhotor
      REAL(r8), DIMENSION(:), ALLOCATABLE :: my_rhotor,my_psitor
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dp_mat
      REAL(r8), DIMENSION(:), ALLOCATABLE :: Re_dp_diagonal,dr_arr,
     $           q_rational,shear,r_o,my_bt0,my_psio,mpsi_arr,
     $           omegas_e_arr,omegas_i_arr,Im_dp_diagonal,dr_vals,
     $           psi_n_rational,dgeo
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ne_arr,te_arr,ni_arr,
     $    ti_arr,zeff_arr,bt_arr,rs_arr,R0_arr,mu_i_arr
      INTEGER,DIMENSION(:),ALLOCATABLE :: nn,resm,nns_arr
      INTEGER :: msing,i,mthsurf
      REAL(r8), DIMENSION(0:512) :: unitfun
      INTEGER :: fsave
      REAL(r8) :: psave
      REAL(r8), DIMENSION(:), ALLOCATABLE :: jacs,delpsi,rsurf,asurf
      REAL(r8) :: rfac,jac,a_surf
c-----------------------------------------------------------------------
c     Read in STRIDE netcdf
c-----------------------------------------------------------------------

      CALL read_stride_netcdf_diagonal(ncfile,msing,dp_mat,
     $           Re_dp_diagonal,Im_dp_diagonal,q_rational,
     $           psi_n_rational,dgeo,shear,r_o,my_bt0,my_psio,dr_vals,
     $           mpsi_arr,nn,resm)
      WRITE(*,*)"msing_out=",msing
      WRITE(*,*)"Re_dp_diagonal=",Re_dp_diagonal
      WRITE(*,*)"Im_dp_diagonal=",Im_dp_diagonal
      WRITE(*,*)"q_rational=",q_rational
      WRITE(*,*)"psi_n_rational=",psi_n_rational
      WRITE(*,*)"dgeo=",dgeo
      WRITE(*,*)"shear=",shear
      WRITE(*,*)"r_o=",r_o
      WRITE(*,*)"my_bt0=",my_bt0
      WRITE(*,*)"my_psio=",my_psio
      WRITE(*,*)"nn=",nn
      WRITE(*,*)"resm=",resm

      mpsi = INT(mpsi_arr(1))
      mthsurf = 512 ! Hardcoded, but this is a default value

c     Allocate SLAYER input type arrays
      ALLOCATE(sl_in%qval_arr(msing),sl_in%omegas_arr(msing),
     $  sl_in%omegas_e_arr(msing),sl_in%dp_matrix(msing,msing),
     $  sl_in%omegas_i_arr(msing),!sl_in%chi_prof_arr(msing),
     $  sl_in%Q_e_arr(msing),sl_in%Q_i_arr(msing),
     $  sl_in%psi_n_arr(msing),
     $  sl_in%Re_dp_arr(msing),sl_in%Im_dp_arr(msing),
     $  sl_in%d_crit_arr(msing),sl_in%P_tor_arr(msing),
     $  sl_in%P_perp_arr(msing),sl_in%tau_arr(msing),
     $  sl_in%D_norm_arr(msing),
     $  sl_in%d_beta_arr(msing),sl_in%gammafac_arr(msing),
     $  sl_in%c_beta_arr(msing),sl_in%lu_arr(msing),
     $  sl_in%Qconv_arr(msing))

c     Allocate local kinetic arrays
      ALLOCATE(ne_arr(msing),te_arr(msing),ni_arr(msing),
     $    ti_arr(msing),zeff_arr(msing),bt_arr(msing),rs_arr(msing),
     $    R0_arr(msing),mu_i_arr(msing),nns_arr(msing),dr_arr(msing),
     $    omegas_e_arr(msing),omegas_i_arr(msing))

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

      ! Input Delta' matrix
      sl_in%dp_matrix(:,:) = CMPLX(dp_mat(:,:,1),dp_mat(:,:,2))

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

         sl_in%omegas_e_arr(ising) = omega_e
         sl_in%omegas_i_arr(ising) = omega_i

         n_e = kin%f(2)
         t_e = kin%f(4)/e
         n_i = kin%f(1)
         t_i = kin%f(3)/e

         zeff = 2.0!kin%f(9)
         
         omega = kin%f(5)
         my_qval = q_rational(ising)!sq%f(4)
         my_sval = shear(ising)
         dgeo_val = dgeo(ising)
         my_bt = my_bt0(1)
         my_rs = a_surf
         R_0 = r_o(1)
         mu_i = 2.0
         dr_val = dr_vals(ising)
         
         chi_s(1) = sl_in%chi_p_arr(ising) ! chi_perp
         chi_s(2) = sl_in%chi_t_arr(ising) ! chi_tor
         chi_s(3) = sl_in%kappa_arr(ising) ! kappa (thermal cond.)

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

         nns_arr(ising) = nn(1)
         nr = nn(1)

         l_n = 0.0
         l_t = 0.0
         WRITE(*,*)"$^$ calling params()"

         CALL params(n_e,t_e,t_i,omega,chi_s,dr_val,dgeo_val,
     $        l_n,l_t,my_qval,my_sval,my_bt,my_rs,R_0,mu_i,zeff,.false.)

!!!!!!!!!!!
         gammafac = (my_rs*Re_dp_diagonal(ising))/tau_r ! scalar to convert thickness into growth rate

         sl_in%qval_arr(ising) = INT(my_qval)
         sl_in%lu_arr(ising)=lu
         sl_in%Q_e_arr(ising)=-tauk*omega_e ! skipping params() calculation
         sl_in%Q_i_arr(ising)=-tauk*omega_i ! skipping params() calculation
         sl_in%c_beta_arr(ising)=c_beta
         sl_in%d_beta_arr(ising)=d_beta
         sl_in%D_norm_arr(ising)=D_norm
         sl_in%tau_arr(ising)=tau
         sl_in%omegas_arr(ising) = omega
         sl_in%psi_n_arr(ising) = respsi
         sl_in%gammafac_arr(ising) = gammafac
         sl_in%Re_dp_arr(ising) = Re_dp_diagonal(ising)
         sl_in%Im_dp_arr(ising) = Im_dp_diagonal(ising)
         sl_in%d_crit_arr(ising) = dc_tmp
         sl_in%P_perp_arr(ising) = P_perp
         sl_in%P_tor_arr(ising) = P_tor
         sl_in%Qconv_arr(ising) = tauk
      ENDDO

      !WRITE(*,*)"msing=",msing
      !WRITE(*,*)"qval_arr=",qval_arr
      !WRITE(*,*)"ne_arr=",ne_arr
      !WRITE(*,*)"te_arr=",te_arr
      !WRITE(*,*)"ni_arr=",ni_arr
      !WRITE(*,*)"ti_arr=",ti_arr
      !WRITE(*,*)"zeff_arr=",zeff_arr
      !WRITE(*,*)"shear=",shear
      !WRITE(*,*)"bt_arr=",bt_arr
      !WRITE(*,*)"rs_arr=",rs_arr
      !WRITE(*,*)"R0_arr=",R0_arr
      !WRITE(*,*)"resm=",resm
      !WRITE(*,*)"nns_arr=",nns_arr
      !WRITE(*,*)"inc_beta_arr=",inc_beta_arr
      !WRITE(*,*)"inds_arr=",inds_arr
      !WRITE(*,*)"intau_arr=",intau_arr
      !WRITE(*,*)"inpr_arr=",inpr_arr
      !WRITE(*,*)"inpe_arr=",inpe_arr
      !WRITE(*,*)"omegas_arr=",omegas_arr
      !WRITE(*,*)"omegas_e_arr=",omegas_e_arr
      !WRITE(*,*)"omegas_i_arr=",omegas_i_arr
      !WRITE(*,*)"Re_deltaprime_arr=",Re_deltaprime_arr
      !WRITE(*,*)"Im_deltaprime_arr=",Im_deltaprime_arr
      !stop
      !CALL slayer_netcdf_inputs(msing,qval_arr,ne_arr,te_arr,ni_arr,
      !$           ti_arr,zeff_arr,shear,bt_arr,rs_arr,R0_arr,
      !$           resm,nns_arr,inc_beta_arr,inds_arr,
      !$           intau_arr,inpr_arr,inpe_arr,inQ_arr,omegas_arr,
      !$           omegas_e_arr,omegas_i_arr,
      !$           Re_deltaprime_arr,Im_deltaprime_arr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod