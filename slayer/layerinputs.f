      MODULE layerinputs_mod

      USE inputs, ONLY : kin
      USE pentrc_interface, ONLY : zi,mi,wefac,wpfac,initialize_pentrc
      USE read_eq_mode, ONLY : read_eq_efit
      USE direct_mod, ONLY : direct_run
      !resm = mfac(resnum(ising))
      !CALL spline_eval(kin,respsi,1)

      ! STILL NEED ro AND bt0, GLOBAL GPEC VARIABLES

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. build_inputs.
c     compute
c-----------------------------------------------------------------------
      SUBROUTINE build_inputs(egnum,xspmn,spot,nspot,
     $             growthrate_flag,
     $             slayer_inpr)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      LOGICAL, INTENT(IN) :: growthrate_flag
      INTEGER, INTENT(IN) :: egnum,nspot
      REAL(r8), INTENT(IN) :: spot, slayer_inpr

      !INTEGER :: i_id,q_id,m_id,p_id,c_id,w_id,k_id,n_id,d_id,a_id,
      !$           pp_id,cp_id,wp_id,np_id,dp_id,wc_id,bc_id,
      !$           astat

      REAL(r8) :: respsi,lpsi,rpsi,shear,hdist,sbnosurf

      INTEGER, DIMENSION(msing) :: resnum

      TYPE(spline_type) :: spl

      INTEGER :: resm
      TYPE(spline_type) :: sr

      REAL(r8), DIMENSION(msing) :: inQ_arr,inQ_e_arr,inQ_i_arr,
     $                       inc_beta_arr,inds_arr,intau_arr,Q0_arr,
     $                       outer_delta_arr

c-----------------------------------------------------------------------
c     Find # of rational surfaces msing.
c-----------------------------------------------------------------------
      CHARACTER(*) :: vmat_filename
      INTEGER :: ising

      OPEN(UNIT=debug_unit,FILE=TRIM(vmat_filename),STATUS="OLD",
     $     FORM="UNFORMATTED")
      READ(debug_unit)mpert,msing

c-----------------------------------------------------------------------
c     Read and build equilibrium inputs
c-----------------------------------------------------------------------
      CALL fourfit_action_matrix
            ! call the automatic reading and distributing of inputs
      CALL initialize_pentrc(op_kin=.FALSE.,op_deq=.FALSE.,
     $          op_peq=.FALSE.)
            ! manually set the pentrc equilibrium description
      CALL set_eq(eqfun,sq,rzphi,smats,tmats,xmats,ymats,zmats,
     $           twopi*psio,ro,nn,jac_type,mlow,mhigh,mpert,mthvac)
            ! manually set the kinetic profiles
      CALL read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,indebug)
        ! manually set the perturbed equilibrium displacements
        ! use false flat xi and xi' for equal weighting
      ALLOCATE(psitmp(sq%mx+1),mtmp(mpert),xtmp(sq%mx+1,mpert))
      psitmp(:) = sq%xs(0:)
      mtmp = (/(m,m=mlow,mhigh)/)
      xtmp = 1e-4
      CALL set_peq(psitmp,mtmp,xtmp,xtmp,xtmp,.false.,tdebug)
      DEALLOCATE(xtmp,mtmp,psitmp)
      IF(verbose) WRITE(*,*)"Computing Kinetic Matrices"
      CALL fourfit_kinetic_matrix(kingridtype,out_fund)

      CALL sing_scan
      DO ising=1,msing
         CALL resist_eval(sing(ising))
      ENDDO

      CALL ksing_find


      CALL gpeq_alloc
      CALL idcon_build(egnum,xspmn)

      CALL gpeq_interp_singsurf(fsp_sol,spot,nspot)

      IF (vsbrzphi_flag) ALLOCATE(singbno_mn(mpert,msing))

      ! minor radius defined using toroidal flux. Used for threshold
      CALL spline_int(sq)
      qintb = sq%fsi(mpsi, 4)
      psitor(:) = sq%fsi(:, 4) / qintb  ! normalized toroidal flux
      rhotor(:) = SQRT(sq%fsi(:, 4)*twopi*psio / (pi * bt0))  ! effective minor radius in Callen
      CALL spline_alloc(sr,mpsi,1)
      sr%xs = sq%xs
      sr%fs(:, 1) = rhotor(:)
      CALL spline_fit(sr,"extrap")

c-----------------------------------------------------------------------
c     HERE I NEED TO READ IN STRIDE NETCDF ??
c-----------------------------------------------------------------------
      IF growthrate_flag THEN
        outer_delta_arr = 0.0
       ! TBD
      END IF
c-----------------------------------------------------------------------
c     loop across singular surfaces, evaluate spline quantities.
c-----------------------------------------------------------------------
      ! j_c is j_c/(chi1*sq%f(4))
      DO ising=1,msing
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)

c-----------------------------------------------------------------------
c     prepare layer analysis.
c-----------------------------------------------------------------------
         resm = mfac(resnum(ising))
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
         qval = sq%f(4)
         sval = sq%f1(4)
         bt = bt0
         rs = sr%f1(1)
         R0 = ro
         mu_i = mi
         inpr = slayer_inpr
         mms = resm
         nns = nn

         mrs = real(mms,4)
         nrs = real(nns,4)

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

         qvals(ising) = qval
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
         outer_delta_arr(ising) = outer_layer_delta

      ENDDO
      CALL spline_dealloc(sr)
      CALL cspline_dealloc(fsp_sol)
      CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE build_inputs

      END MODULE layerinputs_mod