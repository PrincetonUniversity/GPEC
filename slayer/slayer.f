c-----------------------------------------------------------------------
c     SLAYER: Slab LAYER linear drift-MHD code.
c     Main driver program.
c
c     Computes tearing-mode layer quantities (inner layer Delta, 
c     growth rates, torque balance, field thresholds) from 
c     slab-geometry drift-MHD matching via a Riccati integration method.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM slayer

      USE sglobal_mod
      USE delta_mod, ONLY: riccati,riccati_f,riccati_del_s,
     $                     riccati_out,parflow_flag,PeOhmOnly_flag
      USE gslayer_mod
      USE layerinputs_mod

      IMPLICIT NONE
c-----------------------------------------------------------------------
c     local scalars — loop indices and counters.
c-----------------------------------------------------------------------
      INTEGER :: i,j,k              ! general loop indices
      INTEGER :: inn                 ! number of input-file surfaces
      INTEGER :: count               ! generic counter
c-----------------------------------------------------------------------
c     local scalars — numerical resolution.
c-----------------------------------------------------------------------
      INTEGER :: inum                ! resolution for 1-D scans
      INTEGER :: jnum                ! resolution for 2-D scan axis 1
      INTEGER :: knum                ! resolution for 2-D scan axis 2
      INTEGER :: Q_num               ! resolution for stab. scan Re(Q)
      INTEGER :: msing_max           ! max number of singular surfaces
      INTEGER :: n_k                 ! number of rational surfaces
c-----------------------------------------------------------------------
c     local scalars — MAXLOC result holder.
c-----------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: iloc  ! result of MAXLOC
c-----------------------------------------------------------------------
c     control flags — workflow.
c-----------------------------------------------------------------------
      LOGICAL :: params_flag         ! compute params from kinetic data
      LOGICAL :: input_flag          ! read multi-surface input file
      LOGICAL :: read_eq             ! read equilibrium files
      LOGICAL :: verbose             ! enable progress messages
      LOGICAL :: params_check        ! print diagnostic output in params
c-----------------------------------------------------------------------
c     control flags — physics modes.
c-----------------------------------------------------------------------
      LOGICAL :: est_gamma_flag      ! estimate growth rate
      LOGICAL :: match_gamma_flag    ! asymptotically matched gamma
      LOGICAL :: fitz_flag           ! use Fitzpatrick layer model
      LOGICAL :: coupling_flag       ! coupled rational surfaces
      LOGICAL :: br_th_flag          ! Br threshold test scan
      LOGICAL :: bal_flag            ! torque balance scan
      LOGICAL :: stability_flag      ! complex-Q delta scan
      LOGICAL :: Pe_flag             ! include electron pressure
c-----------------------------------------------------------------------
c     control flags — parameter-space scans.
c-----------------------------------------------------------------------
      LOGICAL :: QPscan_flag         ! (Q,P) scan
      LOGICAL :: QPescan_flag        ! (Q,Pe) scan
      LOGICAL :: QPscan2_flag        ! (Q,P) scan variant 2
      LOGICAL :: QDscan2_flag        ! (Q,D) scan variant 2
      LOGICAL :: Qbscan_flag         ! (Q,beta) scan
      LOGICAL :: Qscan_flag          ! 1-D Q scan
      LOGICAL :: onscan_flag         ! (omega,n) scan
      LOGICAL :: otscan_flag         ! (omega,T) scan
      LOGICAL :: ntscan_flag         ! (n,T) scan
      LOGICAL :: nbtscan_flag        ! (n,Bt) scan
      LOGICAL :: riccatiscan_flag    ! Riccati-variable scan
      LOGICAL :: stabscan_flag       ! stability scan (single surface)
      LOGICAL :: coupled_stabscan_flag ! stability scan (coupled)
      LOGICAL :: amr_flag            ! adaptive mesh refinement scan
c-----------------------------------------------------------------------
c     control flags — output format.
c-----------------------------------------------------------------------
      LOGICAL :: ascii_flag          ! write ASCII output files
      LOGICAL :: bin_flag            ! write binary output files
      LOGICAL :: netcdf_flag         ! write NetCDF output files
c-----------------------------------------------------------------------
c     local scalars — physical input quantities.
c-----------------------------------------------------------------------
      REAL(r8) :: n_e                ! electron density [m^-3]
      REAL(r8) :: t_e                ! electron temperature [eV]
      REAL(r8) :: t_i                ! ion temperature [eV]
      REAL(r8) :: omega              ! toroidal rotation [rad/s]
      REAL(r8) :: omega0             ! unused
      REAL(r8) :: l_n                ! density gradient scale length
      REAL(r8) :: l_t                ! temperature gradient scale length
      REAL(r8) :: qval               ! safety factor at surface
      REAL(r8) :: sval               ! magnetic shear at surface
      REAL(r8) :: bt                 ! toroidal field [T]
      REAL(r8) :: rs                 ! minor radius of surface [m]
      REAL(r8) :: R0                 ! major radius [m]
      REAL(r8) :: mu_i               ! ion mass number
      REAL(r8) :: zeff               ! effective charge
      REAL(r8) :: dr_val             ! radial derivative parameter
      REAL(r8) :: dgeo_val           ! geometric factor parameter
      REAL(r8) :: scan_width         ! half-width of complex-Q scan
c-----------------------------------------------------------------------
c     local scalars — normalized layer parameters (namelist overrides).
c-----------------------------------------------------------------------
      REAL(r8) :: inQ                ! normalized ExB rotation freq.
      REAL(r8) :: inQ_e              ! normalized electron diamagnetic
      REAL(r8) :: inQ_i              ! normalized ion diamagnetic
      REAL(r8) :: inpr               ! normalized pressure gradient
      REAL(r8) :: inpe               ! normalized electron pressure
      REAL(r8) :: inc_beta           ! normalized beta
      REAL(r8) :: inds               ! normalized D (magnetic diffusion)
      REAL(r8) :: intau              ! normalized tau = T_i/T_e
      REAL(r8) :: inlu               ! normalized Lundquist number
c-----------------------------------------------------------------------
c     local scalars — derived / scratch quantities.
c-----------------------------------------------------------------------
      REAL(r8) :: psi0               ! reconnected flux (a.u.)
      REAL(r8) :: jxb                ! j x B torque (a.u.)
      REAL(r8) :: Q0                 ! unperturbed rotation frequency
      REAL(r8) :: Q_sol              ! solved rotation frequency
      REAL(r8) :: br_th              ! radial field threshold
      REAL(r8) :: Qratio             ! Q_e/Q ratio for scan2 variants
c-----------------------------------------------------------------------
c     local scalars — scan grid helpers.
c-----------------------------------------------------------------------
      REAL(r8) :: inQ_min,inQ_max    ! rotation scan bounds
      REAL(r8) :: j_min,j_max,jpower ! 2-D scan axis 1
      REAL(r8) :: k_min,k_max,kpower ! 2-D scan axis 2
      REAL(r8) :: ing_step           ! growth-rate grid step
      REAL(r8) :: ing_coarse         ! Re(gamma) grid value
      REAL(r8) :: iing_coarse        ! Im(gamma) grid value
c-----------------------------------------------------------------------
c     local scalars — complex quantities.
c-----------------------------------------------------------------------
      COMPLEX(r8) :: delta           ! layer Delta (tearing index)
      COMPLEX(r8) :: delta_n_p       ! Deltaprime scale factor
      COMPLEX(r8) :: dels_db         ! delta_s / d_beta
      COMPLEX(r8) :: del_s           ! delta_s
      COMPLEX(r8) :: ingamma         ! initial gamma guess (namelist)
      COMPLEX(r8) :: delta_prime     ! external Deltaprime (namelist)
c-----------------------------------------------------------------------
c     local arrays — transport profile coefficients.
c-----------------------------------------------------------------------
      REAL(r8) :: chis(3)            ! chi_perp, chi_tor, kappa
      REAL(r8), DIMENSION(8) :: chi_p_prof   ! chi_perp radial profile
      REAL(r8), DIMENSION(8) :: chi_t_prof   ! chi_tor  radial profile
      REAL(r8), DIMENSION(8) :: kappa_prof   ! kappa    radial profile
c-----------------------------------------------------------------------
c     local arrays — multi-surface input-file storage.
c-----------------------------------------------------------------------
      INTEGER,  DIMENSION(:), ALLOCATABLE :: mms,nns
      REAL(r8), DIMENSION(:), ALLOCATABLE :: prs,n_es,t_es,t_is,
     $     omegas,l_ns,l_ts,svals,qvals,bts,rss,R0s,mu_is,zeffs,
     $     Q_soll,br_thl,pes
c-----------------------------------------------------------------------
c     local arrays — scan workspace.
c-----------------------------------------------------------------------
      REAL(r8), DIMENSION(:),   ALLOCATABLE :: inQs,iinQs
      REAL(r8), DIMENSION(:),   ALLOCATABLE :: jxbl,bal
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: js,ks,psis,jxbs,
     $     Q_sols,br_ths
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Q_solss,br_thss
      COMPLEX(r8), DIMENSION(:),   ALLOCATABLE :: deltal,outer_deltas
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas
c-----------------------------------------------------------------------
c     AMR (adaptive mesh refinement) growth-rate scan.
c-----------------------------------------------------------------------
      INTEGER :: AMR_passes          ! number of AMR refinement levels
      INTEGER :: m_AMR               ! effective number of surfaces
c-----------------------------------------------------------------------
c     structured input/output types (defined in sglobal_mod).
c-----------------------------------------------------------------------
      TYPE(slayer_inputs_type)  :: sl_in
      TYPE(slayer_outputs_type) :: sl_out
      TYPE(deltas_outputs_type), ALLOCATABLE :: all_deltas_out(:)
c-----------------------------------------------------------------------
c     file-path strings.
c-----------------------------------------------------------------------
      CHARACTER(512) :: infile       ! multi-surface input file path
      CHARACTER(512) :: ncfile       ! NetCDF equilibrium file path

c-----------------------------------------------------------------------
c     namelist groups.
c-----------------------------------------------------------------------
      NAMELIST/slayer_input/input_flag,infile,
     $    ncfile,params_flag,mm,nn,n_e,t_e,t_i,sval,bt,rs,R0,omega,
     $    l_t,l_n,qval,mu_i,zeff,dr_val,dgeo_val,chi_p_prof,
     $    chi_t_prof,kappa_prof,inpr,inpe,inQ,inQ_e,inQ_i,inc_beta,
     $    inds,intau,Q0,delta_prime,delta_n_p,ingamma
      NAMELIST/slayer_control/inum,jnum,knum,Q_num,scan_width,
     $    AMR_passes,msing_max,dc_type,read_eq,fitz_flag,coupling_flag,
     $    QPscan_flag,Qscan_flag,QPescan_flag,Qbscan_flag,onscan_flag,
     $    otscan_flag,ntscan_flag,nbtscan_flag,parflow_flag,
     $    peohmonly_flag,Pe_flag,layfac
      NAMELIST/slayer_output/verbose,ascii_flag,bin_flag,netcdf_flag,
     $    est_gamma_flag,match_gamma_flag,stability_flag,
     $    stabscan_flag,coupled_stabscan_flag,amr_flag,br_th_flag,
     $    bal_flag
      NAMELIST/slayer_diagnose/riccati_out,riccatiscan_flag,
     $    params_check
c-----------------------------------------------------------------------
c     set initial values.
c     defaults are overridden by namelist reads below.
c-----------------------------------------------------------------------

      ! mode numbers (sglobal_mod: INTEGER mm,nn; REAL mr,nr)
      mm   = 0
      nn   = 0
      mr   = 0.0
      nr   = 0.0

      ! kinetic / equilibrium inputs
      n_e     = 0.0
      t_e     = 0.0
      t_i     = 0.0
      omega   = 0.0
      l_n     = 0.0
      l_t     = 0.0
      qval    = 0.0
      sval    = 0.0
      bt      = 0.0
      rs      = 0.0
      R0      = 0.0
      mu_i    = 0.0
      zeff    = 0.0
      dr_val  = 0.0
      dgeo_val= 0.0

      ! normalized layer-parameter overrides
      inQ      = 0.0
      inQ_e    = 0.0
      inQ_i    = 0.0
      inpr     = 0.0
      inpe     = 0.0
      inc_beta = 0.0
      inds     = 0.0
      intau    = 0.0
      inlu     = 0.0
      Q0       = 0.0

      ! transport profile coefficients
      chi_p_prof = 0.0
      chi_t_prof = 0.0
      kappa_prof = 0.0
      chis       = 0.0

      ! complex namelist inputs
      delta_prime = (0.0,0.0)
      delta_n_p   = (0.0,0.0)
      ingamma     = (0.0,0.0)

      ! global module scalars (sglobal_mod)
      gamma_fac = 0.0
      dc_type   = ""

      ! scan resolution defaults
      inum       = 400   ! 1-D resolution (error-field threshold scans)
      jnum       = 500   ! 2-D scan axis-1 resolution
      knum       = 100   ! 2-D scan axis-2 resolution
      Q_num      = 100   ! stability scan Re(Q) resolution
      scan_width = 2.0   ! half-width for complex-Q scans
      AMR_passes = 4     ! AMR refinement levels
      msing_max  = 2     ! max singular surfaces to process

      ! I/O unit numbers (sglobal_mod)
      in_unit    = 1
      out_unit   = 2
      out2_unit  = 3
      out3_unit  = 4
      bin_unit   = 5
      bin_2d_unit= 6
      input_unit = 7

      ! workflow flags
      read_eq              = .FALSE.
      est_gamma_flag       = .FALSE.
      match_gamma_flag     = .FALSE.
      fitz_flag            = .FALSE.
      coupling_flag        = .FALSE.
      params_flag          = .TRUE.
      input_flag           = .FALSE.

      ! parameter-space scan flags
      QPscan_flag          = .FALSE.
      QPescan_flag         = .FALSE.
      Qbscan_flag          = .FALSE.
      onscan_flag          = .FALSE.
      otscan_flag          = .FALSE.
      ntscan_flag          = .FALSE.
      nbtscan_flag         = .FALSE.

      ! physics / model flags
      layfac               = 0.02
      Qratio               = 0.5
      parflow_flag         = .FALSE.
      PeOhmOnly_flag       = .TRUE.
      Pe_flag              = .FALSE.

      ! file paths
      infile               = ""
      ncfile               = ""

      ! output control
      verbose              = .TRUE.
      ascii_flag           = .TRUE.
      bin_flag             = .TRUE.
      netcdf_flag          = .FALSE.

      ! diagnostic flags
      riccati_out          = .FALSE.
      riccatiscan_flag     = .FALSE.
      params_check         = .FALSE.

      ! remaining physics-mode flags
      bal_flag             = .FALSE.
      stability_flag       = .FALSE.
      stabscan_flag        = .FALSE.
      coupled_stabscan_flag= .FALSE.
      amr_flag             = .FALSE.
      br_th_flag           = .FALSE.
c-----------------------------------------------------------------------
c     read slayer.in.
c     four namelist groups: input, control, output, diagnose.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)""
      IF(verbose) WRITE(*,*)"SLAYER START"
      IF(verbose) WRITE(*,*)"__________________________________________"
      OPEN(UNIT=in_unit,FILE="slayer.in",STATUS="OLD")
      READ(in_unit,NML=slayer_input)
      READ(in_unit,NML=slayer_control)
      READ(in_unit,NML=slayer_output)
      READ(in_unit,NML=slayer_diagnose)
      CLOSE(UNIT=in_unit)

      ! Build toroidal-mode-number string for output filenames.
      IF (nn<10) THEN
         WRITE(UNIT=sn_str,FMT='(I1)') nn
         sn_str=ADJUSTL(sn_str)
      ELSE
         WRITE(UNIT=sn_str,FMT='(I2)') nn
      ENDIF
c-----------------------------------------------------------------------
c     compute normalized layer parameters from kinetic inputs.
c     params() (params_mod) converts dimensional plasma profiles into
c     the normalized quantities (Q, Q_e, Q_i, c_beta, ds, tau, lu)
c     used by the Riccati solver.
c-----------------------------------------------------------------------
      IF (params_flag) THEN
         CALL params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)
         ! Copy module-level results into local working variables.
         inQ=Q
         inQ_e=Q_e
         inQ_i=Q_i
         inc_beta=c_beta
         inds=ds
         intau=tau
         Q0=Q
      ELSE
         lu=inlu   ! manual Lundquist number when params not computed
      ENDIF
c-----------------------------------------------------------------------
c     baseline single-surface delta, reconnected flux, & torque.
c     skipped when the matched-gamma path is active (it computes
c     its own delta internally).
c-----------------------------------------------------------------------
      IF (.NOT. (match_gamma_flag)) THEN
         delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
         psi0=1.0/ABS(delta+delta_n_p)     ! reconnected flux  [a.u.]
         jxb=-AIMAG(1.0/(delta+delta_n_p)) ! j x B torque      [a.u.]
         IF (verbose) THEN
            WRITE(*,*)"delta=",delta
            WRITE(*,*)"psi0=",psi0
            WRITE(*,*)"jxb=",jxb
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     multi-surface input-file mode.
c     reads an external file of (m,n) surfaces with per-surface kinetic
c     profiles, computes delta and the error-field threshold for each.
c-----------------------------------------------------------------------
      IF (input_flag) THEN
         OPEN(UNIT=input_unit,FILE=infile,STATUS="old")
         READ(input_unit,*)inn
         ALLOCATE(mms(0:inn-1),nns(0:inn-1),prs(0:inn-1),
     $        n_es(0:inn-1),t_es(0:inn-1),t_is(0:inn-1),
     $        omegas(0:inn-1),
     $        l_ns(0:inn-1),l_ts(0:inn-1),qvals(0:inn-1),
     $        svals(0:inn-1),
     $        bts(0:inn-1),rss(0:inn-1),R0s(0:inn-1),
     $        mu_is(0:inn-1),zeffs(0:inn-1),
     $        Q_soll(0:inn-1),br_thl(0:inn-1))
         DO k=0,inn-1
            READ(input_unit,'(2(1x,I2),14(1x,e12.4))')
     $           mms(k),nns(k),prs(k),
     $           n_es(k),t_es(k),t_is(k),omegas(k),
     $           l_ns(k),l_ts(k),qvals(k),svals(k),
     $           bts(k),rss(k),R0s(k),mu_is(k),zeffs(k)
         ENDDO
         CLOSE(input_unit)
         
         DO k=0,inn-1
            WRITE(*,*)k    ! surface index
            mr=REAL(mms(k))
            nr=REAL(nns(k))
            inpr=prs(k)
            CALL params(n_es(k),t_es(k),t_is(k),omegas(k),chis,dr_val,
     $           dgeo_val,l_ns(k),l_ts(k),qvals(k),svals(k),bts(k),
     $           rss(k),R0s(k),mu_is(k),zeffs(k),params_check)
            inQ=Q
            inQ_e=Q_e
            inQ_i=Q_i
            inc_beta=c_beta
            inds=ds
            intau=tau
            Q0=Q
 
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
            riccati_out=.FALSE.
         
            ALLOCATE(inQs(0:inum),deltal(0:inum),
     $           jxbl(0:inum),bal(0:inum)) 
         
            DO i=0,inum
               inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
               deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $              inpr,inc_beta,inds,intau,inpe)
               jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
               bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
            ENDDO
        
            iloc=MAXLOC(bal)
            Q_soll(k)=inQs(iloc(1))
            br_thl(k)=sqrt(MAXVAL(bal)/lu*(svals(k)**2.0/2.0))*1e4

            IF (verbose) WRITE(*,*)"Q_sol=",Q_soll(k)
            IF (verbose) WRITE(*,*)"br_th=",br_thl(k)
            DEALLOCATE(inQs,deltal,jxbl,bal)         
         ENDDO
         OPEN(UNIT=out_unit,FILE="slayer_input_bal_n"//
     $      TRIM(sn_str)//".out",STATUS="UNKNOWN")
         WRITE(out_unit,'(1x,(2a17))') "Q_sol","br_th"
         
         DO k=0,inn-1
            WRITE(out_unit,'(1x,2(es17.8e3))')
     $           Q_soll(k),br_thl(k)
         ENDDO
         CLOSE(out_unit)
         DEALLOCATE(prs,n_es,t_es,t_is,omegas,l_ns,l_ts,qvals,svals,
     $        bts,rss,R0s,mu_is,zeffs,Q_soll,br_thl,mms,nns)
      ENDIF
c-----------------------------------------------------------------------
c     estimate growth rate via resistive-layer thickness.
c     Uses riccati_del_s to get delta_s/d_beta, then scales by
c     d_beta to obtain the layer thickness delta_s and the estimated 
c     gamma. Inputs may come from equilibrium files (read_eq) or 
c     namelist.Subroutines: build_inputs, allocate_inputs, 
c     allocate_outputsare defined in gslayer_mod / layerinputs_mod.
c-----------------------------------------------------------------------
      IF (est_gamma_flag) THEN
      WRITE(*,*)"------------------------------------------"
      WRITE(*,*)">>> Estimating growth rate"

         IF (read_eq) THEN
            ! Read equilibrium files for multi-surface inputs.
            ! build_inputs (layerinputs_mod) reads STRIDE NetCDF data.
            sl_in%chi_p_arr = chi_p_prof
            sl_in%chi_t_arr = chi_t_prof
            sl_in%kappa_arr = kappa_prof

            CALL build_inputs(infile,ncfile,sl_in)

            n_k = SIZE(sl_in%qval_arr)
            CALL allocate_outputs(n_k,sl_out)

         ELSE
            ! Single-surface mode: build inputs from namelist.
            n_k = 1
            mr = mm
            nr = nn

            chis(1) = chi_p_prof(1) ! chi_perp
            chis(2) = chi_t_prof(1) ! chi_tor
            chis(3) = kappa_prof(1) ! kappa (thermal cond.)

            CALL params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

            ! Override computed parameters with nonzero namelist values.
            IF (ABS(inQ) > 0.0) THEN
               Q = inQ ! NAMELIST
            END IF
            IF (ABS(inQ_e) > 0.0) THEN
               Q_e = inQ_e ! NAMELIST
            END IF
            IF (ABS(inQ_i) > 0.0) THEN
               Q_i = inQ_i ! NAMELIST
            END IF
            IF (inpr > 0.0) THEN
               pr = inpr ! NAMELIST
            END IF
            IF (intau > 0.0) THEN
               tau = intau ! NAMELIST
            END IF
            IF (inds > 0.0) THEN
               D_norm = inds ! NAMELIST
            END IF

            CALL allocate_inputs(n_k,sl_in)  ! gslayer_mod
            CALL allocate_outputs(n_k,sl_out) ! gslayer_mod

            sl_in%qval_arr = (/ qval /)
            sl_in%omegas_arr = (/ omega /)
            !sl_in%Q_arr = (/ Q /)
            sl_in%Q_e_arr = (/ Q_e /)
            sl_in%Q_i_arr = (/ Q_i /)
            sl_in%psi_n_arr = (/ 0.0 /)
            sl_in%Re_dp_arr = (/ REAL(delta_prime) /)
            sl_in%Im_dp_arr = (/ AIMAG(delta_prime) /)
            sl_in%d_crit_arr = (/ dc_tmp /)
            sl_in%P_perp_arr = (/ P_perp /)
            sl_in%P_tor_arr = (/ P_tor /)
            sl_in%tau_arr = (/ tau /)
            sl_in%D_norm_arr = (/ D_norm /)
            sl_in%d_beta_arr = (/ d_beta /)
            sl_in%gammafac_arr = (/ gamma_fac /)
            sl_in%c_beta_arr = (/ c_beta /)
            sl_in%lu_arr = (/ lu /)
            sl_in%Qconv_arr = (/ tauk /)
         END IF 

         ! Loop over rational surfaces to estimate growth rates.
         DO k=1,n_k
            WRITE(*,*)
            WRITE(*,'(A,I0,A)') 'Calculating growth rate '//
     $             'estimate on q = ',
     $       sl_in%qval_arr(k),' rational surface'

            D_norm = sl_in%D_norm_arr(k)
            ! First arg is Q_e (electron diamagnetic freq), not Q
            ! (ExB freq).  This is intentional per riccati_del_s API.
            dels_db=riccati_del_s(sl_in%Q_e_arr(k),
     $                   sl_in%Q_i_arr(k),sl_in%P_perp_arr(k),
     $                   5.0*sl_in%D_norm_arr(k))

            del_s = dels_db * sl_in%d_beta_arr(k)

            sl_out%gamma_est_arr(k) = sl_in%gammafac_arr(k)/del_s
            sl_out%dels_db_arr(k) = dels_db
            WRITE(*,*)
            WRITE(*,'(A,F0.3,A)')'Growth rate estimate = ',
     $                 REAL(sl_out%gamma_est_arr(k)),' [Hz]'
            
         ENDDO

         IF (.NOT. (match_gamma_flag)) THEN
            sl_out%gamma_sol_arr = (/0./)
            CALL output_gamma(est_gamma_flag,m_AMR,sl_in,sl_out,
     $      all_deltas_out)
         END IF
      ENDIF
c-----------------------------------------------------------------------
c     asymptotically matched growth rate.
c     Matches the inner-layer Delta to the outer-region Delta' to
c     find the self-consistent complex growth rate.  Supports both
c     single-surface and coupled multi-surface (AMR) modes.
c     Subroutines: dispersion_AMR_v2, dispersion_det (gslayer_mod),
c                  riccati_f (delta_mod).
c-----------------------------------------------------------------------
      IF (match_gamma_flag) THEN
         WRITE(*,*)"------------------------------------------"
         WRITE(*,*)">>> Calculating asymptotically matched growth rate"

         IF (read_eq) THEN
         
            IF (.NOT. est_gamma_flag) THEN
               sl_in%chi_p_arr = chi_p_prof
               sl_in%chi_t_arr = chi_t_prof
               sl_in%kappa_arr = kappa_prof

               CALL build_inputs(infile,ncfile,sl_in)

               n_k = SIZE(sl_in%qval_arr)
               CALL allocate_outputs(n_k,sl_out)
            END IF
         ELSE
            n_k = 1

            IF (.NOT. est_gamma_flag) THEN

            chis(1) = chi_p_prof(1) ! chi_perp
            chis(2) = chi_t_prof(1) ! chi_tor
            chis(3) = kappa_prof(1) ! kappa (thermal cond.)

            WRITE(*,*)"chis(1) (chi_perp): ",chis(1)
            WRITE(*,*)"chis(2) (chi_tor): ",chis(2)
            WRITE(*,*)"chis(3) (kappa): ",chis(3)

            ! Use namelist kinetic inputs instead of equilibrium files
            CALL params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

            ! Override desired normalized parameters
            IF (ABS(inQ) > 0.0) THEN
               Q = inQ ! NAMELIST
            END IF
            IF (ABS(inQ_e) > 0.0) THEN
               Q_e = inQ_e ! NAMELIST
            END IF
            IF (ABS(inQ_i) > 0.0) THEN
               Q_i = inQ_i ! NAMELIST
            END IF
            IF (inpr > 0.0) THEN
               P_perp = inpr ! NAMELIST
            END IF
            IF (intau > 0.0) THEN
               tau = intau ! NAMELIST
            END IF
            IF (inds > 0.0) THEN
               D_norm = inds ! NAMELIST
            END IF

            IF (.NOT. est_gamma_flag) THEN
               CALL allocate_inputs(n_k,sl_in)
               CALL allocate_outputs(n_k,sl_out)
            END IF

            sl_in%qval_arr = (/ qval /)
            sl_in%omegas_arr = (/ omega /)
            !sl_in%Q_arr = (/ Q /)
            sl_in%Q_e_arr = (/ Q_e /)
            sl_in%Q_i_arr = (/ Q_i /)
            sl_in%psi_n_arr = (/ 0.0 /)
            sl_in%Re_dp_arr = (/ REAL(delta_prime) /)
            sl_in%Im_dp_arr = (/ AIMAG(delta_prime) /)
            sl_in%d_crit_arr = (/ dc_tmp /)
            sl_in%P_perp_arr = (/ P_perp /)
            sl_in%P_tor_arr = (/ P_tor /)
            sl_in%tau_arr = (/ tau /)
            sl_in%D_norm_arr = (/ D_norm /)
            sl_in%d_beta_arr = (/ d_beta /)
            sl_in%gammafac_arr = (/ gamma_fac /)
            sl_in%c_beta_arr = (/ c_beta /)
            sl_in%lu_arr = (/ lu /)
            sl_in%Qconv_arr = (/ tauk /)

            END IF
         END IF 

c-----------------------------------------------------------------------
c     allocate output arrays for AMR delta storage.
c-----------------------------------------------------------------------
         IF (AMR_flag .AND. .NOT. coupling_flag) THEN
            ALLOCATE(all_deltas_out(n_k))
         ELSEIF (AMR_flag .AND. coupling_flag) THEN
            ALLOCATE(all_deltas_out(1))
         END IF

         IF (AMR_flag) THEN
            IF (coupling_flag) THEN
               m_AMR = 1
            ELSE
               m_AMR = MIN(n_k,msing_max)
            END IF
         END IF

         WRITE(*,*),"Rational q domain: ",sl_in%qval_arr
c-----------------------------------------------------------------------
c     loop over rational surfaces to find matched growth rates.
c-----------------------------------------------------------------------
         DO k=1,MIN(n_k,msing_max)
            WRITE(*,*)
            WRITE(*,'(A,I0,A)') 'Calculating growth rate on q = ',
     $       sl_in%qval_arr(k),' rational surface:'

            ! Load per-surface parameters into module-level scalars.
            Q_e = sl_in%Q_e_arr(k)
            Q_i = sl_in%Q_i_arr(k)
            P_perp = sl_in%P_perp_arr(k)
            P_tor = sl_in%P_tor_arr(k)
            tau = sl_in%tau_arr(k)
            D_norm = sl_in%D_norm_arr(k)
            c_beta = sl_in%c_beta_arr(k)
            tauk = sl_in%Qconv_arr(k)
            iota_e = Q_e / (Q_e - Q_i)

            WRITE(*,*)"Q_e: ",Q_e
            WRITE(*,*)"Q_i: ",Q_i
            WRITE(*,*)"P_perp: ",P_perp
            WRITE(*,*)"P_tor: ",P_tor
            WRITE(*,*)"tau: ",tau
            WRITE(*,*)"D_norm: ",D_norm
            WRITE(*,*)"tauk: ",tauk
            WRITE(*,*)"iota_e: ",iota_e
            WRITE(*,*)"Delta_prime: ",sl_in%Re_dp_arr(k)
            WRITE(*,*)"Delta_crit: ",sl_in%d_crit_arr(k)

            ! Calculate (Deltaprime - D_crit)/S^1/3
            delta_eff = (sl_in%Re_dp_arr(k) - 
     $          sl_in%d_crit_arr(k))/(sl_in%lu_arr(k)**(1.0/3.0))
            pe = 0.0

            ! Placeholder: gamma_sol_arr filled for external root-finder.
            sl_out%gamma_sol_arr(k) = 0.0

c-----------------------------------------------------------------------
c     uncoupled AMR scan (one surface at a time).
c     dispersion_AMR_v2 (gslayer_mod) populates Q_store, D_store.
c-----------------------------------------------------------------------
            IF (AMR_flag .AND. .NOT. coupling_flag) THEN

               WRITE(*,'(A,I0,A)') 'Calling uncoupled AMR scan on q = ',
     $       sl_in%qval_arr(k),' rational surface:'
               CALL dispersion_AMR_v2(n_k,sl_in,msing_max,scan_width,
     $                  Q_num,AMR_passes,coupling_flag)

               ! Re-allocate output arrays for this surface.
               IF (ALLOCATED(all_deltas_out(k)%inQs)) 
     $             DEALLOCATE(all_deltas_out(k)%inQs)
               IF (ALLOCATED(all_deltas_out(k)%iinQs)) 
     $             DEALLOCATE(all_deltas_out(k)%iinQs)
               IF (ALLOCATED(all_deltas_out(k)%real_deltas)) 
     $             DEALLOCATE(all_deltas_out(k)%real_deltas)
               IF (ALLOCATED(all_deltas_out(k)%imag_deltas)) 
     $             DEALLOCATE(all_deltas_out(k)%imag_deltas)

               ALLOCATE(all_deltas_out(k)%inQs(n_pts),
     $                  all_deltas_out(k)%iinQs(n_pts))
               ALLOCATE(all_deltas_out(k)%real_deltas(n_pts), 
     $         all_deltas_out(k)%imag_deltas(n_pts)) 

               ! Flatten unique AMR points into 1-D output arrays.
               DO i = 1, n_pts
                  all_deltas_out(k)%inQs(i) = REAL(Q_store(i))
                  all_deltas_out(k)%iinQs(i) = -AIMAG(Q_store(i))

                  all_deltas_out(k)%real_deltas(i) = REAL(D_store(i))
                  all_deltas_out(k)%imag_deltas(i) = AIMAG(D_store(i))
               END DO
         
               ! Clean up temporary AMR memory.
               DEALLOCATE(Q_store, D_store)

               WRITE(*,'(A,I2,A,I7,A,2ES14.6)')
     $          '   Surface', k, ': n_pts=', 
     $          SIZE(all_deltas_out(k)%real_deltas),
     $          ' out_chksum=',
     $          SUM(all_deltas_out(k)%real_deltas),
     $          SUM(all_deltas_out(k)%imag_deltas)

            END IF

c-----------------------------------------------------------------------
c     single-surface stability scan on [Re(Q), Im(Q)] grid.
c     Uses riccati_f() (new Fitzpatrick TJ-like formalism) or 
c.    riccati() (orig. SLAYER/Waelbroeck).
c-----------------------------------------------------------------------
            IF ((stabscan_flag)) THEN
               WRITE(*,*)"------------------------------------------"
               WRITE(*,'(A,F0.1)')' >>> Running [Re(Q),'//
     $            'Im(Q)] scan with Q width = ',
     $                scan_width

               ing_step = (2.0*scan_width) / (Q_num - 1)
               count = 0
         
               ALLOCATE(inQs(1:(Q_num+1)),iinQs(1:Q_num))
               ALLOCATE(deltas(1:(Q_num+1),1:Q_num))

               DO i = 1, (Q_num+1)
                  DO j = 1, Q_num
                     ing_coarse = -scan_width + (i - 1) * ing_step ! added 0.5*
                     iing_coarse = -scan_width + (j - 1) * ing_step
                     ! Evaluate riccati function
                     g_tmp = CMPLX(ing_coarse,iing_coarse)
                     IF (fitz_flag) THEN
                        delta=riccati_f()
                     ELSE
                        delta=riccati(iing_coarse,Q_e,Q_i,P_perp,
     $                             c_beta,D_norm,tau,pe,
     $                             iinQ=ing_coarse)
                     END IF
                     inQs(i) = ing_coarse
                     IF (fitz_flag) THEN
                        iinQs(j) = iing_coarse
                     ELSE
                        iinQs(j) = -iing_coarse
                     END IF
                     deltas(i,j) = delta
                  ENDDO
               ENDDO


               IF (k<10) THEN
                  WRITE(UNIT=sm_str,FMT='(I1)') sl_in%qval_arr(k)
                  sm_str=ADJUSTL(sm_str)
               ELSE
                  WRITE(UNIT=sm_str,FMT='(I2)') sl_in%qval_arr(k)
               ENDIF

               OPEN(UNIT=out_unit,FILE="slayer_stability_n"//
     $         TRIM(sn_str)//"m"//TRIM(sm_str)//".out", STATUS="UNKNOWN")
               WRITE(out_unit,'(1x,4(a17))') "RE(Q)",
     $           "IM(Q)","RE(delta)","IM(delta)"
               DO i=1,Q_num+1
                  DO j=1,Q_num
                     WRITE(out_unit,'(1x,4(es17.8e3))')
     $                 inQs(i),iinQs(j),
     $                 REAL(deltas(i,j)),AIMAG(deltas(i,j))
                  ENDDO
               ENDDO
               CLOSE(out_unit)

            DEALLOCATE(inQs,iinQs,deltas)
            ENDIF 
         ENDDO 

         IF (.NOT. (est_gamma_flag)) THEN
            sl_in%d_beta_arr = (/ 0. /)
            sl_out%dels_db_arr = (/ 0. /)
         END IF

c-----------------------------------------------------------------------
c     coupled AMR scan (all surfaces simultaneously).
c     dispersion_AMR_v2 (gslayer_mod) with coupling_flag = .TRUE.
c-----------------------------------------------------------------------
         IF (AMR_flag .AND. coupling_flag) THEN

            CALL dispersion_AMR_v2(n_k,sl_in,msing_max,scan_width,
     $                  Q_num,AMR_passes,coupling_flag)

            ! Re-allocate output arrays.
            IF (ALLOCATED(all_deltas_out(1)%inQs)) 
     $             DEALLOCATE(all_deltas_out(1)%inQs)
            IF (ALLOCATED(all_deltas_out(1)%iinQs)) 
     $             DEALLOCATE(all_deltas_out(1)%iinQs)
            IF (ALLOCATED(all_deltas_out(1)%real_deltas)) 
     $             DEALLOCATE(all_deltas_out(1)%real_deltas)
            IF (ALLOCATED(all_deltas_out(1)%imag_deltas)) 
     $             DEALLOCATE(all_deltas_out(1)%imag_deltas)

            ALLOCATE(all_deltas_out(1)%inQs(n_pts),
     $                  all_deltas_out(1)%iinQs(n_pts))
            ALLOCATE(all_deltas_out(1)%real_deltas(n_pts), 
     $         all_deltas_out(1)%imag_deltas(n_pts)) 

            ! Flatten unique AMR points into 1-D output arrays.
            DO i = 1, n_pts
               all_deltas_out(1)%inQs(i) = REAL(Q_store(i))
               
               all_deltas_out(1)%iinQs(i) = -AIMAG(Q_store(i)) ! verify this sign convention

               all_deltas_out(1)%real_deltas(i) = REAL(D_store(i))
               all_deltas_out(1)%imag_deltas(i) = AIMAG(D_store(i))
            END DO
      
            ! Clean up temporary AMR memory.
            DEALLOCATE(Q_store, D_store)

         END IF

c-----------------------------------------------------------------------
c     coupled-surface stability scan on [Re(Q), Im(Q)] grid.
c     Uses dispersion_det (gslayer_mod) for the full dispersion
c     determinant including inter-surface coupling.
c-----------------------------------------------------------------------
         IF (coupled_stabscan_flag) THEN
            WRITE(*,*)"------------------------------------------"
            WRITE(*,'(A,F0.1)')' >>> Running [Re(Q),'//
     $            'Im(Q)] determinant scan with radius = ',
     $                scan_width

            ing_step = (2.0 * scan_width) / (Q_num - 1)
            count = 0

            !IF (.NOT. stabscan_flag) THEN
            ALLOCATE(inQs(1:(Q_num+1)),iinQs(1:Q_num))
            ALLOCATE(deltas(1:(Q_num+1),1:Q_num))
            !END IF

            inQs=0.0; iinQs=0.0; deltas=(0.0,0.0)

            DO i = 1, (Q_num+1)
               DO j = 1, Q_num
                  ing_coarse = -scan_width + (i - 1) * ing_step
                  iing_coarse = -scan_width + (j - 1) * ing_step
                  inQs(i) = ing_coarse
                  iinQs(j) = iing_coarse

                  ! Evaluate determinant
                  g_tmp = CMPLX(ing_coarse,iing_coarse)
                  deltas(i,j)=dispersion_det(g_tmp,n_k,sl_in,msing_max)
               ENDDO
            ENDDO

            OPEN(UNIT=out_unit,FILE="slayer_determinants_n"//
     $         TRIM(sn_str)//".out", STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,4(a17))') "RE(Q)",
     $           "IM(Q)","RE(det)","IM(det)"
            DO i=1,Q_num+1
               DO j=1,Q_num
                  WRITE(out_unit,'(1x,4(es17.8e3))')
     $                 inQs(i),iinQs(j),
     $                 REAL(deltas(i,j)),AIMAG(deltas(i,j))
               ENDDO
            ENDDO
            CLOSE(out_unit)

         DEALLOCATE(inQs,iinQs,deltas)
         END IF

         CALL output_gamma(est_gamma_flag,m_AMR,sl_in,sl_out,
     $                     all_deltas_out)
         stop
      ENDIF  ! match_gamma_flag
c-----------------------------------------------------------------------
c     Br threshold test scan (analytic).
c     For testing & verification only.  Scans rotation to find the
c     critical radial-field threshold from a simple torque balance.
c-----------------------------------------------------------------------
      IF (br_th_flag) THEN
         WRITE(*,*)"------------------------------------------"
         WRITE(*,*)">>> Computing Br threshold"

         IF (read_eq) THEN
            sl_in%chi_p_arr = chi_p_prof
            sl_in%chi_t_arr = chi_t_prof
            sl_in%kappa_arr = kappa_prof
            CALL build_inputs(infile,ncfile,sl_in)
            n_k = SIZE(sl_in%qval_arr)
            CALL allocate_outputs(n_k,sl_out)
         ELSE
            n_k = 1

            chis(1) = chi_p_prof(1)
            chis(2) = chi_t_prof(1)
            chis(3) = kappa_prof(1)

            CALL params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $           l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

            inQ=Q
            inQ_e=Q_e
            inQ_i=Q_i
            inc_beta=c_beta
            inds=ds
            intau=tau

            CALL allocate_inputs(n_k,sl_in)
            CALL allocate_outputs(n_k,sl_out)

            sl_in%qval_arr    = (/ qval /)
            sl_in%omegas_arr  = (/ omega /)
            sl_in%Q_e_arr     = (/ Q_e /)
            sl_in%Q_i_arr     = (/ Q_i /)
            sl_in%psi_n_arr   = (/ 0.0 /)
            sl_in%Re_dp_arr   = (/ 0.0 /)
            sl_in%Im_dp_arr   = (/ 0.0 /)
            sl_in%d_crit_arr  = (/ 0.0 /)
            sl_in%P_perp_arr  = (/ P_perp /)
            sl_in%P_tor_arr   = (/ P_tor /)
            sl_in%tau_arr     = (/ tau /)
            sl_in%D_norm_arr  = (/ D_norm /)
            sl_in%d_beta_arr  = (/ d_beta /)
            sl_in%gammafac_arr = (/ gamma_fac /)
            sl_in%c_beta_arr  = (/ c_beta /)
            sl_in%lu_arr      = (/ lu /)
            sl_in%Qconv_arr   = (/ tauk /)
         END IF
c-----------------------------------------------------------------------
c     loop over rational surfaces to compute Br threshold.
c-----------------------------------------------------------------------
         delta_n_p = 1e-2
         inum = 200
         inQ_max = 10.0
         inQ_min = -10.0

         DO k=1,n_k
            WRITE(*,*)
            WRITE(*,'(A,I0,A)') 'Computing Br threshold on q = ',
     $         sl_in%qval_arr(k),' rational surface'

            Q_e   = sl_in%Q_e_arr(k)
            Q_i   = sl_in%Q_i_arr(k)
            inQ_e = Q_e
            inQ_i = Q_i
            inpr  = sl_in%P_perp_arr(k)
            inc_beta = sl_in%c_beta_arr(k)
            inds  = sl_in%D_norm_arr(k)
            intau = sl_in%tau_arr(k)
            Q0    = sl_in%Q_e_arr(k)

            ALLOCATE(inQs(0:inum),deltal(0:inum),
     $               jxbl(0:inum),bal(0:inum))
            DO i=0,inum
               inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
               deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $                      inpr,inc_beta,inds,intau,inpe)
               jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
               bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
            ENDDO

            iloc=MAXLOC(bal)
            Q_sol=inQs(iloc(1))
            br_th=SQRT(MAXVAL(bal)/sl_in%lu_arr(k)
     $           *(sval**2.0/2.0))

            sl_out%br_th_arr(k) = br_th
            sl_out%gamma_sol_arr(k) = 0.0
            sl_out%gamma_est_arr(k) = 0.0
            sl_out%dels_db_arr(k)   = 0.0

            WRITE(*,'(A,ES12.4)') '  br_th = ', br_th
            DEALLOCATE(inQs,deltal,jxbl,bal)
         ENDDO
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
         CALL output_gamma(est_gamma_flag,m_AMR,sl_in,sl_out,
     $                     all_deltas_out)
         STOP
      ENDIF  ! br_th_flag
c-----------------------------------------------------------------------
c     find solutions based on simple torque balance.
c-----------------------------------------------------------------------
      IF (bal_flag)THEN
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
         ! just for diagnostics
         inQ_max=10.0
         inQ_min=-10.0
         
         riccati_out=.FALSE.
         
         ALLOCATE(inQs(0:inum),deltal(0:inum),jxbl(0:inum),bal(0:inum)) 
         
         DO i=0,inum
            inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
            deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $           inpr,inc_beta,inds,intau,inpe)
            jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
            bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
         ENDDO

         ! write components of torque balance
         IF(ascii_flag)THEN
            OPEN(UNIT=out_unit,FILE="slayer_bal_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,5(a17))') "inQ","RE(delta)",
     $           "IM(delta)","jxb","bal"

            DO i=0,inum
               WRITE(out_unit,'(1x,5(es17.8e3))') inQs(i),
     $              REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
            ENDDO
            CLOSE(out_unit)
         ENDIF

         iloc=MAXLOC(bal)
         Q_sol=inQs(iloc(1))
         br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))*1e4
         WRITE(*,*)"Q_sol=",Q_sol
         WRITE(*,*)"br_th=",br_th
         DEALLOCATE(inQs,deltal,jxbl,bal)         
      ENDIF
c-----------------------------------------------------------------------
c     examine delta dependencies on complex Q for stability.
c-----------------------------------------------------------------------
      IF (stability_flag) THEN
         ALLOCATE(inQs(0:inum),iinQs(0:200))
         ALLOCATE(deltas(0:inum,0:200))

         inQ_max=10.0
         inQ_min=-10.0

         DO i=0,inum
            DO j=0,200
               inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
               iinQs(j)=inQ_min+(REAL(j)/200)*(inQ_max-inQ_min)
               deltas(i,j)=riccati(inQs(i),inQ_e,inQ_i,inpr,inc_beta,
     $              inds,intau,inpe,iinQ=iinQs(j))
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_stability_n"//
     $         TRIM(sn_str)//".out", STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,4(a17))') "RE(Q)",
     $           "IM(Q)","RE(delta)","IM(delta)"
            DO i=0,inum
               DO j=0,200
                  WRITE(out_unit,'(1x,4(es17.8e3))')
     $                 inQs(i),iinQs(j),
     $                 REAL(deltas(i,j)),AIMAG(deltas(i,j))
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF
         DEALLOCATE(inQs,iinQs,deltas)
      ENDIF
c-----------------------------------------------------------------------
c     riccati scan.
c-----------------------------------------------------------------------
      IF (riccatiscan_flag) THEN
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        deltas(0:jnum,0:knum))

         j_min=0.0  ! in log scale
         j_max=2.0  ! in log scale
         DO j=0,jnum
            jpower=j_min+(j_max-j_min)/jnum*REAL(j)
            js(j,:)=10.0**jpower
            DO k=0,knum
               ks(j,k)=inc_beta/sqrt((1+intau)*inds)*js(j,k)**2.0
               deltas(j,k)=riccati(inQ,inQ_e,inQ_i,inpr,
     $              inc_beta,inds,intau,inpe,inx=js(j,k),
     $              iny=ks(j,k)*EXP(ifac*2*pi*REAL(k)/knum))
               WRITE(*,*)"deltas",deltas(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_riccatiscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,5(a17))') "x","yphs","yamp",
     $           "RE(delta)","IM(delta)"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,5(es17.8e3))')
     $                 js(j,k),2*pi*REAL(k)/knum,ks(j,k),
     $                 REAL(deltas(j,k)),AIMAG(deltas(j,k))
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

      DEALLOCATE(js,ks,deltas)
      ENDIF
c-----------------------------------------------------------------------
c     (Q,Pe) scan.
c-----------------------------------------------------------------------
      IF (QPescan_flag) THEN
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        deltas(0:jnum,0:knum),psis(0:jnum,0:knum),
     $        jxbs(0:jnum,0:knum))

         j_min=0.05 ! extended from 20.0
         j_max=50.0 ! extended from 20.0
         k_min=-3.0 ! in log scale
         k_max=0.0 ! in log scale
         DO j=0,jnum            
            js(j,:)=j_min+(j_max-j_min)/jnum*REAL(j)
            DO k=0,knum
               kpower=k_min+(k_max-k_min)/knum*REAL(k)
               ks(j,k)=10.0**kpower
               deltas(j,k)=riccati(js(j,k),inQ_e,inQ_i,inpr,
     $              inc_beta,inds,intau,ks(j,k))
               psis(j,k)=1.0/ABS(deltas(j,k)+delta_n_p)
               jxbs(j,k)=-AIMAG(1.0/(deltas(j,k)+delta_n_p))
               WRITE(*,*)"deltas",deltas(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_QPescan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Q","Pe","RE(delta)",
     $           "IM(delta)","psi","jxb"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 js(j,k),ks(j,k),REAL(deltas(j,k)),
     $                 AIMAG(deltas(j,k)),psis(j,k),jxbs(j,k)
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_QPescan_n'
     $         //TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(js,4),REAL(ks,4)
            WRITE(bin_2d_unit)REAL(REAL(deltas),4)
            WRITE(bin_2d_unit)REAL(AIMAG(deltas),4)
            WRITE(bin_2d_unit)REAL(psis,4)
            WRITE(bin_2d_unit)REAL(jxbs,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(js,ks,deltas,psis,jxbs)
      ENDIF
c-----------------------------------------------------------------------
c     (Q,P) scan.
c-----------------------------------------------------------------------
      IF (QPscan_flag) THEN
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        deltas(0:jnum,0:knum),psis(0:jnum,0:knum),
     $        jxbs(0:jnum,0:knum))

         j_min=0.05 ! extended from 20.0
         j_max=50.0 ! extended from 20.0
         k_min=-3.0 ! in log scale
         k_max=0.0 ! in log scale
         DO j=0,jnum            
            js(j,:)=j_min+(j_max-j_min)/jnum*REAL(j)
            DO k=0,knum
               kpower=k_min+(k_max-k_min)/knum*REAL(k)
               ks(j,k)=10.0**kpower
               deltas(j,k)=riccati(js(j,k),inQ_e,inQ_i,60.6*ks(j,k),
     $              inc_beta,inds,intau,inpe)
               psis(j,k)=1.0/ABS(deltas(j,k)+delta_n_p)
               jxbs(j,k)=-AIMAG(1.0/(deltas(j,k)+delta_n_p))
               WRITE(*,*)"deltas",deltas(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_QPscan"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Q","Pr","RE(delta)",
     $           "IM(delta)","psi","jxb"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 js(j,k),ks(j,k),REAL(deltas(j,k)),
     $                 AIMAG(deltas(j,k)),psis(j,k),jxbs(j,k)
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE="slayer_QPscan_"
     $         //TRIM(sn_str)//".bin",
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(js,4),REAL(ks,4)
            WRITE(bin_2d_unit)REAL(REAL(deltas),4)
            WRITE(bin_2d_unit)REAL(AIMAG(deltas),4)
            WRITE(bin_2d_unit)REAL(psis,4)
            WRITE(bin_2d_unit)REAL(jxbs,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(js,ks,deltas,psis,jxbs)
      ENDIF
c-----------------------------------------------------------------------
c     (Q) scan.
c-----------------------------------------------------------------------
      IF (Qscan_flag) THEN
         ALLOCATE(js(0:jnum,0:1),
     $        deltas(0:jnum,0:1),psis(0:jnum,0:1),
     $        jxbs(0:jnum,0:1))

         j_min=0.05 ! extended from 20.0
         j_max=50.0 ! extended from 20.0
         DO j=0,jnum            
            js(j,:)=j_min+(j_max-j_min)/jnum*REAL(j)
            deltas(j,0)=riccati(js(j,0),inQ_e,inQ_i,inpr,
     $         inc_beta,inds,intau,inpe)
            psis(j,0)=1.0/ABS(deltas(j,0)+delta_n_p)
            jxbs(j,0)=-AIMAG(1.0/(deltas(j,0)+delta_n_p))
            WRITE(*,*)"deltas",deltas(j,0)
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_Qscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Q","Pr","RE(delta)",
     $           "IM(delta)","psi","jxb"
            DO j=0,jnum
               WRITE(out_unit,'(1x,6(es17.8e3))')
     $              js(j,0),inpr,REAL(deltas(j,0)),
     $              AIMAG(deltas(j,0)),psis(j,0),jxbs(j,0)
            ENDDO
            CLOSE(out_unit)
         ENDIF

      DEALLOCATE(js,deltas,psis,jxbs)
      ENDIF
c-----------------------------------------------------------------------
c     (Q,P) scan 2.
c-----------------------------------------------------------------------
      IF (QPscan2_flag) THEN
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        deltas(0:jnum,0:knum),psis(0:jnum,0:knum),
     $        jxbs(0:jnum,0:knum))

         j_min=-3.0 ! in log scale
         j_max=1.7 ! in log scale for cb02_ds005_Qr05, for cb01_ds02_Qr05
         !j_max=1.4 ! in log scale for cb03_ds40_Qr05
         k_min=-4.0 ! in log scale
         k_max=3.0 ! in log scale
         DO j=0,jnum
            jpower=j_min+(j_max-j_min)/jnum*REAL(j)
            js(j,:)=10.0**jpower
            DO k=0,knum
               kpower=k_min+(k_max-k_min)/knum*REAL(k)
               ks(j,k)=10.0**kpower
               deltas(j,k)=riccati(js(j,k),js(j,k)*Qratio,
     $              -js(j,k)*Qratio,ks(j,k),inc_beta,inds,intau,inpe)
               psis(j,k)=1.0/ABS(deltas(j,k)+delta_n_p)
               jxbs(j,k)=-AIMAG(1.0/(deltas(j,k)+delta_n_p))
               WRITE(*,*)"deltas",deltas(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_QPscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Q","Pr","RE(delta)",
     $           "IM(delta)","psi","jxb"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 js(j,k),ks(j,k),REAL(deltas(j,k)),
     $                 AIMAG(deltas(j,k)),psis(j,k),jxbs(j,k)
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE="slayer_QPscan_n"//
     $         TRIM(sn_str)//".bin",
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(js,4),REAL(ks,4)
            WRITE(bin_2d_unit)REAL(REAL(deltas),4)
            WRITE(bin_2d_unit)REAL(AIMAG(deltas),4)
            WRITE(bin_2d_unit)REAL(psis,4)
            WRITE(bin_2d_unit)REAL(jxbs,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(js,ks,deltas,psis,jxbs)
      ENDIF
c-----------------------------------------------------------------------
c     (Q,D) scan 2.
c-----------------------------------------------------------------------
      IF (QDscan2_flag) THEN
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        deltas(0:jnum,0:knum),psis(0:jnum,0:knum),
     $        jxbs(0:jnum,0:knum))

         j_min=1.0 
         j_max=100.0 
         k_min=1.0 ! in log scale
         k_max=100.0 ! in log scale
         DO j=0,jnum
            js(j,:)=j_min+(j_max-j_min)/jnum*REAL(j)
            DO k=0,knum
               ks(j,k)=k_min+(k_max-k_min)/knum*REAL(k)
               deltas(j,k)=riccati(js(j,k),js(j,k)*Qratio,
     $              -js(j,k)*Qratio,inpr,inc_beta,ks(j,k),intau,inpe)
               psis(j,k)=1.0/ABS(deltas(j,k)+delta_n_p)
               jxbs(j,k)=-AIMAG(1.0/(deltas(j,k)+delta_n_p))
               WRITE(*,*)"deltas",deltas(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="QDscan.out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Q","D","RE(delta)",
     $           "IM(delta)","psi","jxb"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 js(j,k),ks(j,k),REAL(deltas(j,k)),
     $                 AIMAG(deltas(j,k)),psis(j,k),jxbs(j,k)
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_QDscan_n'//
     $         TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(js,4),REAL(ks,4)
            WRITE(bin_2d_unit)REAL(REAL(deltas),4)
            WRITE(bin_2d_unit)REAL(AIMAG(deltas),4)
            WRITE(bin_2d_unit)REAL(psis,4)
            WRITE(bin_2d_unit)REAL(jxbs,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(js,ks,deltas,psis,jxbs)
      ENDIF
c-----------------------------------------------------------------------
c     (o,n) scan.
c-----------------------------------------------------------------------
      IF (onscan_flag) THEN
         ALLOCATE(inQs(0:inum),bal(0:inum)) 
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        Q_sols(0:jnum,0:knum),br_ths(0:jnum,0:knum))

         j_min=-2.0 ! -2.0
         j_max=5.0 ! 5.0
         k_min=0.2
         k_max=8.0

         DO j=0,jnum
            js(j,:)=j_min+(j_max-j_min)*(REAL(j)/jnum)
            DO k=0,knum
               ks(j,k)=k_min+(k_max-k_min)*(REAL(k)/knum)
               
               CALL params(n_e*ks(j,k),t_e,t_i,omega*js(j,k),chis,
     $                    dr_val,dgeo_val,l_n,l_t,qval,sval,bt,rs,
     $                    R0,mu_i,zeff,params_check)
               inQ=Q
               inQ_e=Q_e
               inQ_i=Q_i
               inc_beta=c_beta
               inds=ds
               intau=tau
               Q0=Q
               

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


               
               DO i=0,inum
                  inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
                  delta=riccati(inQs(i),inQ_e,inQ_i,
     $                 inpr,inc_beta,inds,intau,inpe)
                  jxb=-AIMAG(1.0/(delta+delta_n_p))
                  bal(i)=2.0*inpr*(Q0-inQs(i))/jxb
               ENDDO
               iloc=MAXLOC(bal)
               Q_sols(j,k)=inQs(iloc(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_onscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Omega","Density",
     $           "Omega_i","Omega_e","Omega_sol","Field_Threshold"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 omega*js(j,k),n_e*ks(j,k),
     $                 omega_i,omega_e,
     $                 Q_sols(j,k),br_ths(j,k)
                  
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_onscan_n'//
     $         TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(omega*js,4),REAL(n_e*ks,4)
            WRITE(bin_2d_unit)REAL(Q_sols,4)
            WRITE(bin_2d_unit)REAL(br_ths,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(inQs,bal,js,ks,Q_sols,br_ths)
      ENDIF
c-----------------------------------------------------------------------
c     (o,t) scan.
c-----------------------------------------------------------------------
      IF (otscan_flag) THEN
         ALLOCATE(inQs(0:inum),bal(0:inum)) 
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        Q_sols(0:jnum,0:knum),br_ths(0:jnum,0:knum))

         j_min=-2.0 ! -2.0
         j_max=5.0 ! 5.0
         k_min=0.2
         k_max=8.0

         DO j=0,jnum
            js(j,:)=j_min+(j_max-j_min)*(REAL(j)/jnum)
            DO k=0,knum
               ks(j,k)=k_min+(k_max-k_min)*(REAL(k)/knum)
               
               CALL params(n_e,t_e*ks(j,k),t_i*ks(j,k),
     $              omega*js(j,k),chis,dr_val,dgeo_val,l_n,l_t,qval,
     $              sval,bt,rs,R0,mu_i,zeff,params_check)
               inQ=Q
               inQ_e=Q_e
               inQ_i=Q_i
               inc_beta=c_beta
               inds=ds
               intau=tau
               Q0=Q
               
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

               
               DO i=0,inum
                  inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
                  delta=riccati(inQs(i),inQ_e,inQ_i,
     $                 inpr,inc_beta,inds,intau,inpe)
                  jxb=-AIMAG(1.0/(delta+delta_n_p))
                  bal(i)=2.0*inpr*(Q0-inQs(i))/jxb
               ENDDO
               iloc=MAXLOC(bal)
               Q_sols(j,k)=inQs(iloc(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"t_e=",t_e*ks(j,k),"br_ths=",br_ths(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_otscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Omega","Temperature",
     $           "Omega_i","Omega_e","Omega_sol","Field_Threshold"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 omega*js(j,k),t_e*ks(j,k),
     $                 omega_i,omega_e,
     $                 Q_sols(j,k),br_ths(j,k)
                  
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_otscan_n'//
     $         TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(omega*js,4),REAL(t_e*ks,4)
            WRITE(bin_2d_unit)REAL(Q_sols,4)
            WRITE(bin_2d_unit)REAL(br_ths,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(inQs,bal,js,ks,Q_sols,br_ths)
      ENDIF
c-----------------------------------------------------------------------
c     (n,t) scan.
c-----------------------------------------------------------------------
      IF (ntscan_flag) THEN
         ALLOCATE(inQs(0:inum),bal(0:inum)) 
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        Q_sols(0:jnum,0:knum),br_ths(0:jnum,0:knum))

         j_min=0.2
         j_max=8.0
         k_min=0.2
         k_max=8.0

         DO j=0,jnum
            js(j,:)=j_min+(j_max-j_min)*(REAL(j)/jnum)
            DO k=0,knum
               ks(j,k)=k_min+(k_max-k_min)*(REAL(k)/knum)
               
               CALL params(n_e*ks(j,k),t_e*js(j,k),t_i*js(j,k),omega,
     $              chis,dr_val,dgeo_val,l_n,l_t,qval,sval,bt,rs,R0,
     $              mu_i,zeff,params_check)
               inQ=Q
               inQ_e=Q_e
               inQ_i=Q_i
               inc_beta=c_beta
               inds=ds
               intau=tau
               Q0=Q
               
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

               
               DO i=0,inum
                  inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
                  delta=riccati(inQs(i),inQ_e,inQ_i,
     $                 inpr,inc_beta,inds,intau,inpe)
                  jxb=-AIMAG(1.0/(delta+delta_n_p))
                  bal(i)=2.0*inpr*(Q0-inQs(i))/jxb
               ENDDO
               iloc=MAXLOC(bal)
               Q_sols(j,k)=inQs(iloc(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_ntscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))') "Temperature","Density",
     $           "Omega_i","Omega_e","Omega_sol","Field_Threshold"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,6(es17.8e3))')
     $                 t_e*js(j,k),n_e*ks(j,k),
     $                 omega_i,omega_e,
     $                 Q_sols(j,k),br_ths(j,k)
               ENDDO
            ENDDO                  
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_ntscan_n'//
     $         TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(t_e*js,4),REAL(n_e*ks,4)
            WRITE(bin_2d_unit)REAL(Q_sols,4)
            WRITE(bin_2d_unit)REAL(br_ths,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(inQs,bal,js,ks,Q_sols,br_ths)
      ENDIF
c-----------------------------------------------------------------------
c     (n,bt) scan.
c-----------------------------------------------------------------------
      IF (nbtscan_flag) THEN
         ALLOCATE(inQs(0:inum),bal(0:inum)) 
         ALLOCATE(js(0:jnum,0:knum),ks(0:jnum,0:knum),
     $        Q_sols(0:jnum,0:knum),br_ths(0:jnum,0:knum))

         j_min=0.3
         j_max=8.0
         k_min=0.2
         k_max=8.0

         DO j=0,jnum
            js(j,:)=j_min+(j_max-j_min)*(REAL(j)/jnum)
            DO k=0,knum
               ks(j,k)=k_min+(k_max-k_min)*(REAL(k)/knum)

             
               CALL params(n_e*ks(j,k),t_e,t_i,omega,chis,dr_val,
     $              dgeo_val,l_n,l_t,qval,sval,bt*js(j,k),rs,R0,mu_i,
     $              zeff,params_check)
               inQ=Q
               inQ_e=Q_e
               inQ_i=Q_i
               inc_beta=c_beta
               inds=ds
               intau=tau
               Q0=Q
           

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

        
               DO i=0,inum
                  inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
                  delta=riccati(inQs(i),inQ_e,inQ_i,
     $                 inpr,inc_beta,inds,intau,inpe)
                  jxb=-AIMAG(1.0/(delta+delta_n_p))
                  bal(i)=2.0*inpr*(Q0-inQs(i))/jxb
               ENDDO
               iloc=MAXLOC(bal)
               Q_sols(j,k)=inQs(iloc(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_nbtscan_n"//
     $         TRIM(sn_str)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,4(a17))') "Bt","Density",
     $           "Omega_sol","Field_Threshold"
            DO j=0,jnum
               DO k=0,knum
                  WRITE(out_unit,'(1x,4(es17.8e3))')
     $                 bt*js(j,k),n_e*ks(j,k),
     $                 Q_sols(j,k),br_ths(j,k)
                  
               ENDDO
            ENDDO
            CLOSE(out_unit)
         ENDIF

         IF (bin_flag) THEN
            OPEN(UNIT=bin_2d_unit,FILE='slayer_nbtscan_n'//
     $         TRIM(sn_str)//'.bin',
     $         STATUS='UNKNOWN',POSITION='REWIND',FORM='UNFORMATTED')
            WRITE(bin_2d_unit)1,0
            WRITE(bin_2d_unit)jnum,knum
            WRITE(bin_2d_unit)REAL(bt*js,4),REAL(n_e*ks,4)
            WRITE(bin_2d_unit)REAL(Q_sols,4)
            WRITE(bin_2d_unit)REAL(br_ths,4)
            CLOSE(bin_2d_unit)
         ENDIF
      DEALLOCATE(inQs,bal,js,ks,Q_sols,br_ths)
      ENDIF

      END PROGRAM slayer
     
