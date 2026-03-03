      MODULE gslayer_mod
c-----------------------------------------------------------------------
c     gslayer_mod: Single-surface SLAYER driver called by GPEC.
c
c     This module contains the gpec_slayer subroutine, which computes
c     Delta, torque, and critical field threshold for one rational
c     surface.  It is the primary interface used by gpec/gpout.f.
c
c     Growth-rate scanning, AMR dispersion solvers, and I/O helpers
c     have been split out into growthrates_mod (growthrates.f).
c
c     Subprograms:
c       1. gpec_slayer           - single-surface delta/torque driver
c-----------------------------------------------------------------------

      USE omp_lib

      USE sglobal_mod, ONLY: out_unit,r8,mu0,m_p,chag,lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,
     $   eta,visc,rho_s,lu,omega_e,omega_i,
     $   delta_n,Q,
     $   ifac,g_tmp,pi,                              ! used by AMR
     $   tauk,iota_e,D_norm,P_perp,P_tor,delta_eff,  ! used by det
     $   amr_cell_type,amr_cells,n_amr_cells,         ! v2 types
     $   Q_store,D_store,n_pts,                       ! output arrays
     $   MAX_PTS,HASH_SZ,HASH_SCALE,                  ! v1 constants
     $   hash_head,hash_next,                          ! v1 hash
     $   MAX_CELLS,                                    ! v2 limit
     $   slayer_inputs_type,slayer_outputs_type,
     $   deltas_outputs_type,
     $   tau_r,dc_tmp,dc_type,
     $   sn_str,sm_str
      USE delta_mod, ONLY: riccati,riccati_f,riccati_out,
     $   parflow_flag,PeOhmOnly_flag
      USE params_mod
      USE layerinputs_mod
      USE slayer_netcdf_mod

      IMPLICIT NONE

c --- reconnection regulariser used in psi0 / JxB expressions;
c     expose via namelist to make user-configurable.
      REAL(r8), PARAMETER :: DELTA_N_PERT = 1.0e-2_r8

      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. gpec_slayer.
c     Single-surface SLAYER driver: compute Delta, torque, and b_crit
c     for one rational surface characterised by its plasma profiles.
c
c     Steps:
c       1. Derive normalised SLAYER parameters (lu, ds, Q, Q_e, etc.)
c          from dimensional inputs.
c       2. Compute baseline Delta via riccati().
c       3. Scan over a range of Q (rotation) to build a torque
c          balance curve, then identify the critical threshold br_th.
c
c     TODO: `zeff` and `qval` are INTENT(IN)
c       arguments accepted in the interface but never referenced in
c       the subroutine body.  Remove from signature (breaking API
c       change) or add a comment documenting their reserved intent.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_slayer(n_e,t_e,n_i,t_i,zeff,omega,omega_e,
     $   omega_i,qval,sval,bt,rs,R0,mu_i,inpr,mms,nns,ascii_flag,
     $     delta,psi0,jxb,omega_sol,br_th)

c --- input arguments: dimensional plasma profiles for this surface
      REAL(r8),INTENT(IN) :: n_e       ! electron density [m^-3]
      REAL(r8),INTENT(IN) :: t_e       ! electron temperature [eV]
      REAL(r8),INTENT(IN) :: n_i       ! ion density [m^-3]
      REAL(r8),INTENT(IN) :: t_i       ! ion temperature [eV]
      REAL(r8),INTENT(IN) :: zeff      ! (TODO: unused; see header)
      REAL(r8),INTENT(IN) :: omega     ! plasma rotation frequency
      REAL(r8),INTENT(IN) :: omega_e   ! electron diamagnetic freq
      REAL(r8),INTENT(IN) :: omega_i   ! ion diamagnetic freq
      REAL(r8),INTENT(IN) :: qval      ! (TODO: unused; see header)
      REAL(r8),INTENT(IN) :: sval      ! magnetic shear
      REAL(r8),INTENT(IN) :: bt        ! toroidal field [T]
      REAL(r8),INTENT(IN) :: rs        ! minor radius of surface [m]
      REAL(r8),INTENT(IN) :: R0        ! major radius [m]
      REAL(r8),INTENT(IN) :: inpr      ! Prandtl number
      INTEGER, INTENT(IN) :: mms       ! poloidal mode number
      INTEGER, INTENT(IN) :: nns       ! toroidal mode number
      INTEGER, INTENT(IN) :: mu_i      ! ion mass number (AMU)
      LOGICAL, INTENT(IN) :: ascii_flag ! write ASCII torque-balance file
c --- output arguments
      COMPLEX(r8),INTENT(OUT) :: delta  ! complex Delta
      COMPLEX(r8),INTENT(OUT) :: psi0   ! reconnected flux (a.u.)
      REAL(r8),INTENT(OUT) :: jxb       ! JxB torque (a.u.)
      REAL(r8),INTENT(OUT) :: omega_sol ! rotation at torque-balance
      REAL(r8),INTENT(OUT) :: br_th     ! critical radial field threshold

c --- loop / scan control
      INTEGER :: i               ! loop index
      INTEGER :: inum             ! number of scan points
      INTEGER, DIMENSION(1) :: max_idx ! index of max torque balance
c --- local copies of normalised parameters for riccati()
      REAL(r8) :: inQ,inQ_e,inQ_i  ! normalised frequencies
      REAL(r8) :: inpe             ! electron Prandtl (set to 0)
      REAL(r8) :: inc_beta,inds    ! c_beta, ds copies
      REAL(r8) :: intau            ! tau copy
c --- derived dimensional quantities
      REAL(r8) :: mrs,nrs          ! real-valued mode numbers
      REAL(r8) :: rho              ! mass density [kg/m^3]
      REAL(r8) :: b_l              ! characteristic magnetic field
      REAL(r8) :: Qconv            ! Q normalisation factor
      REAL(r8) :: Q0               ! initial Q (before scan)
      REAL(r8) :: delta_n_p        ! Delta_n perturbation
      REAL(r8) :: lbeta            ! local beta
      REAL(r8) :: tau_i            ! ion collision time
      REAL(r8) :: tau_h            ! Alfven time across surface
      REAL(r8) :: tau_v            ! viscous diffusion time
c --- scan workspace
      REAL(r8) :: inQ_min,inQ_max  ! scan bounds in Q
      REAL(r8) :: Q_sol            ! Q at torque-balance threshold
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs  ! Q scan array
      REAL(r8), DIMENSION(:), ALLOCATABLE :: jxbl   ! JxB scan array
      REAL(r8), DIMENSION(:), ALLOCATABLE :: bal    ! balance scan array
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal ! delta scan array
      CHARACTER(3) :: l_sn,l_sm   ! local mode number strings

c --- set delta_mod control flags for this single-surface call
      parflow_flag=.FALSE.
      PeOhmOnly_flag=.TRUE.
      riccati_out=.FALSE.

      mrs = real(mms,4)
      nrs = real(nns,4)
c-----------------------------------------------------------------------
c     build string representations of m and n for file names.
c-----------------------------------------------------------------------
      IF (nns<10) THEN
         WRITE(UNIT=l_sn,FMT='(I1)') nns
         l_sn=ADJUSTL(l_sn)
      ELSE
         WRITE(UNIT=l_sn,FMT='(I2)') nns
      ENDIF
      IF (mms<10) THEN
         WRITE(UNIT=l_sm,FMT='(I1)') mms
         l_sm=ADJUSTL(l_sm)
      ELSEIF (mms<100) THEN
         WRITE(UNIT=l_sm,FMT='(I2)') mms
         l_sm=ADJUSTL(l_sm)
      ELSE
         WRITE(UNIT=l_sm,FMT='(I3)') mms
      ENDIF

      inpe=0.0                  ! electron Prandtl not used here

c-----------------------------------------------------------------------
c     derive normalised SLAYER parameters from dimensional inputs.
c     tau, eta, rho, b_l, v_a are intermediate dimensional quantities;
c     lu (Lundquist), ds, Qconv, Q, Q_e, Q_i, c_beta are the
c     normalised parameters used by riccati().
c-----------------------------------------------------------------------
      tau= t_i/t_e                     ! ion-to-electron temp ratio
      tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion coll. time
      eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! Spitzer resistivity (Wesson)
      rho=(mu_i*m_p)*n_e               ! mass density

      b_l=(nrs/mrs)*nrs*sval*bt/R0     ! characteristic magnetic field
      rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Larmor radius at T_e

      tau_h=R0*(mu0*rho)**0.5/(nns*sval*bt) ! Alfven transit time
      tau_r=mu0*rs**2.0/eta            ! resistive diffusion time
      tau_v=tau_r/inpr                  ! viscous diffusion time
      visc= rho*rs**2.0/tau_v          ! back-calculated viscosity

      lu=tau_r/tau_h                    ! Lundquist number
      Qconv=lu**(1.0/3.0)*tau_h        ! Q normalisation factor

c --- normalised frequencies
      Q=Qconv*omega
      Q_e=-Qconv*omega_e
      Q_i=-Qconv*omega_i

      ds=lu**(1.0/3.0)*rho_s/rs        ! normalised ion sound radius

      lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
      c_beta=(lbeta/(1.0+lbeta))**0.5  ! compressibility parameter

      delta_n=lu**(1.0/3.0)/rs         ! Delta normalisation factor

c --- copy normalised values into local variables for riccati() call
      inQ=Q
      inQ_e=Q_e
      inQ_i=Q_i
      inc_beta=c_beta
      inds=ds
      intau=tau
      Q0=Q
c-----------------------------------------------------------------------
c     compute baseline Delta, reconnected flux, and JxB torque.
c     delta_n_p uses the module-level DELTA_N_PERT constant; promote
c     that constant to a namelist input to make it user-configurable.
c-----------------------------------------------------------------------
      delta_n_p = DELTA_N_PERT
      delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
      psi0=1.0/ABS(delta+delta_n_p)     ! reconnected flux (a.u.)
      jxb=-AIMAG(1.0/(delta+delta_n_p)) ! JxB torque (a.u.)
c-----------------------------------------------------------------------
c     torque-balance scan: sweep Q over [inQ_min, inQ_max] and
c     compute delta(Q), JxB(Q), and the balance parameter.
c     The threshold br_th is sqrt(max(bal)/lu * s^2/2).
c
c     TODO: physics bounds from Q0/Q_e (computed
c       in the IF/ELSE above) are immediately overridden by the two
c       fixed assignments below.  To use physics bounds, remove the
c       inQ_max=10 / inQ_min=-10 lines.  To keep fixed bounds, remove
c       the dead IF/ELSE block above.
c-----------------------------------------------------------------------
      IF (Q0>inQ_e) THEN
         inQ_max=2.0*Q0
         inQ_min=1.05*inQ_e
      ELSE
         inQ_max=0.95*inQ_e
         IF (Q0>0) THEN
            inQ_min=0.8*inQ_i
         ELSE
            inQ_min=1.5*MINVAL((/Q0,inQ_i/))
         ENDIF
      ENDIF

      inQ_max=10.0              ! TODO: fixed bound; see header
      inQ_min=-10.0
      inum=200
      ALLOCATE(inQs(0:inum),deltal(0:inum),jxbl(0:inum),bal(0:inum))
      DO i=0,inum
         inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
         deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $        inpr,inc_beta,inds,intau,inpe)
         jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
         bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
      ENDDO

c --- optionally write torque balance curve to ASCII file
      IF(ascii_flag)THEN
         OPEN(UNIT=out_unit,FILE="gpec_slayer_torque_balance_m"//
     $        TRIM(l_sm)//"_n"//TRIM(l_sn)//".OUT",
     $        STATUS="UNKNOWN")
         WRITE(out_unit,'(1x,5(a17))') "inQ","RE(delta)",
     $        "IM(delta)","jxb","bal"
         DO i=0,inum
            WRITE(out_unit,'(1x,5(es17.8e3))')
     $           inQs(i),REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
         ENDDO
         CLOSE(out_unit)
      ENDIF

c --- identify the critical threshold from the maximum balance value
      max_idx=MAXLOC(bal)
      Q_sol=inQs(max_idx(1))
      omega_sol=inQs(max_idx(1))/Qconv
      br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))
      DEALLOCATE(inQs,deltal,jxbl,bal)

      RETURN
      END SUBROUTINE gpec_slayer
      END MODULE gslayer_mod