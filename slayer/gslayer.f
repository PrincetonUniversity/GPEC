      MODULE gslayer_mod
c-----------------------------------------------------------------------
c     gslayer_mod: SLAYER driver routines and AMR dispersion scanners.
c
c     Subprograms:
c       1. gpec_slayer           - single-surface delta/torque driver
c       2. output_gamma          - write results to netCDF via
c                                  slayer_netcdf_mod
c       3. allocate_inputs       - allocate slayer_inputs_type arrays
c       4. allocate_outputs      - allocate slayer_outputs_type arrays
c       5. shrink_array          - trim over-allocated scan arrays
c       6. grow_array            - expand scan arrays dynamically
c       7. calc_determinant      - 2x2 / 3x3 complex determinant
c       8. dispersion_det        - coupled dispersion determinant
c       9. dispersion_AMR        - AMR scan v1 (hash-based dedup)
c      10. dispersion_AMR_v2     - AMR scan v2 (cell-based storage)
c
c     Helper subroutines (v1): get_or_compute
c     Helper subroutines (v2): get_or_compute_v2,
c       compute_delta_sub, check_cell_crossing_sub,
c       subdivide_cell_sub, flatten_cells_to_points_sub
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
c-----------------------------------------------------------------------
c     subprogram 2. output_gamma.
c     Write growth-rate results to netCDF via slayer_netcdf_out.
c     Passes the input/output structured types and AMR results
c     through to the netCDF writer (slayer_netcdf_mod).
c-----------------------------------------------------------------------
      SUBROUTINE output_gamma(est_gamma_flag,m_AMR,sl_in,sl_out,
     $                        all_deltas_out)

      LOGICAL, INTENT(IN) :: est_gamma_flag  ! single-surface mode?
      INTEGER, INTENT(IN) :: m_AMR           ! AMR pass count
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      TYPE(slayer_outputs_type), INTENT(IN) :: sl_out
      TYPE(deltas_outputs_type), INTENT(IN) ::
     $                            all_deltas_out(SIZE(sl_in%qval_arr))

      CALL slayer_netcdf_out(SIZE(sl_in%qval_arr),m_AMR,est_gamma_flag,
     $                       sl_in,sl_out,all_deltas_out)

      END SUBROUTINE output_gamma
c-----------------------------------------------------------------------
c     subprogram 3. allocate_inputs.
c     Allocate all per-surface arrays inside slayer_inputs_type.
c-----------------------------------------------------------------------
      SUBROUTINE allocate_inputs(n_k,sl_in)

      INTEGER, INTENT(IN) :: n_k                           ! number of surfaces
      TYPE(slayer_inputs_type), INTENT(INOUT) :: sl_in

      ALLOCATE(sl_in%qval_arr(n_k),sl_in%omegas_arr(n_k),
     $  sl_in%Q_e_arr(n_k),sl_in%Q_i_arr(n_k),sl_in%psi_n_arr(n_k),
     $  sl_in%Re_dp_arr(n_k),sl_in%Im_dp_arr(n_k),
     $  sl_in%d_crit_arr(n_k),sl_in%P_tor_arr(n_k),
     $  sl_in%P_perp_arr(n_k),sl_in%tau_arr(n_k),
     $  sl_in%D_norm_arr(n_k),
     $  sl_in%d_beta_arr(n_k),sl_in%gammafac_arr(n_k),
     $  sl_in%c_beta_arr(n_k),sl_in%lu_arr(n_k),sl_in%Qconv_arr(n_k))
      RETURN
      END SUBROUTINE allocate_inputs
c-----------------------------------------------------------------------
c     subprogram 4. allocate_outputs.
c     Allocate per-surface arrays inside slayer_outputs_type.
c-----------------------------------------------------------------------
      SUBROUTINE allocate_outputs(n_k,sl_out)

      INTEGER, INTENT(IN) :: n_k                            ! number of surfaces
      TYPE(slayer_outputs_type), INTENT(INOUT) :: sl_out

      ALLOCATE(sl_out%dels_db_arr(n_k),sl_out%gamma_sol_arr(n_k),
     $         sl_out%gamma_est_arr(n_k),sl_out%br_th_arr(n_k)  )
      RETURN
      END SUBROUTINE allocate_outputs
c-----------------------------------------------------------------------
c     subprogram 5. shrink_array.
c     Trim an over-allocated REAL(r8) array down to new_size using
c     MOVE_ALLOC (no copy of trailing elements).
c-----------------------------------------------------------------------
      SUBROUTINE shrink_array(arr, new_size)

          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: new_size  ! target size
          REAL(r8), ALLOCATABLE :: temp(:)  ! temporary buffer

          ALLOCATE(temp(new_size))
          temp(1:new_size) = arr(1:new_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE shrink_array
c-----------------------------------------------------------------------
c     subprogram 6. grow_array.
c     Expand a REAL(r8) array from old_size to new_size, preserving
c     existing data via MOVE_ALLOC.
c-----------------------------------------------------------------------
      SUBROUTINE grow_array(arr, old_size, new_size)

          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: old_size   ! current valid element count
          INTEGER, INTENT(IN) :: new_size   ! target allocation size
          REAL(r8), ALLOCATABLE :: temp(:)   ! temporary buffer

          ALLOCATE(temp(new_size))
          temp(1:old_size) = arr(1:old_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE grow_array
c-----------------------------------------------------------------------
c     subprogram 7. calc_determinant.
c     Compute the determinant of a 2x2 or 3x3 complex matrix.
c     Returns (0,0) and sets status=-1 for unsupported sizes.
c     status=0 on success, -1 when nk is neither 2 nor 3.
c-----------------------------------------------------------------------
      SUBROUTINE calc_determinant(matk, nk, detk, status)

      IMPLICIT NONE

c --- arguments
      INTEGER, INTENT(IN) :: nk                            ! matrix rank (2 or 3)
      COMPLEX(r8), DIMENSION(nk,nk), INTENT(IN) :: matk    ! input matrix
      COMPLEX(r8), INTENT(OUT) :: detk                     ! determinant result
      INTEGER, INTENT(OUT) :: status        ! 0=success, -1=unsupported rank
                        
      status = 0  ! Initialize status as success
            
      SELECT CASE (nk)
        CASE (2)
        ! 2x2 determinant: ad - bc
        detk = matk(1,1) * matk(2,2) - matk(1,2) * matk(2,1)
                    
        CASE (3)
        ! 3x3 determinant using cofactor expansion along first row
        detk = matk(1,1)*(matk(2,2)*matk(3,3)-matk(2,3)
     $      *matk(3,2))-matk(1,2)*(matk(2,1)*matk(3,3)
     $      -matk(2,3)*matk(3,1))+matk(1,3)*(matk(2,1)
     $      *matk(3,2)-matk(2,2)*matk(3,1))
                    
        CASE default
        ! Unsupported matrix size
            detk = (0.0, 0.0)
            status = -1
                    
      END SELECT
      RETURN
      END SUBROUTINE calc_determinant
c-----------------------------------------------------------------------
c     subprogram 8. dispersion_det.
c     Compute the coupled dispersion determinant for n_k surfaces.
c
c     For n_k = 1 (single surface):
c       Evaluate riccati_f() (uses module g_tmp), de-normalise by lu^(1/3), and
c       return Deltaprime - delta(Q).
c
c     For n_k = 2 or 3 (coupled surfaces):
c       Build the diagonal delta(Q) matrix, subtract from dp_matrix,
c       and return det(dp_matrix - delta_Q).
c
c-----------------------------------------------------------------------
      FUNCTION dispersion_det(g_in,n_k,sl_in,msing_max)

c --- arguments
      COMPLEX(r8), INTENT(IN) :: g_in       ! complex growth rate
      INTEGER, INTENT(IN) :: n_k            ! number of surfaces
      INTEGER, INTENT(IN) :: msing_max      ! max surfaces to include
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
c --- function result and locals
      COMPLEX(r8) :: dispersion_det         ! returned determinant
      COMPLEX(r8) :: det_val                ! intermediate determinant
      COMPLEX(r8) :: tmp_delta              ! single-surface delta
      COMPLEX(r8), ALLOCATABLE :: delta_Q(:,:)      ! diagonal delta matrix
      COMPLEX(r8), ALLOCATABLE :: result_matrix(:,:) ! dp - delta_Q
      INTEGER :: k                          ! surface loop index
      INTEGER :: det_status                 ! calc_determinant error status

c --- single-surface branch
      IF (msing_max < 2) THEN
c        set module-level variables for riccati_f
         Q_e = sl_in%Q_e_arr(1)
         Q_i = sl_in%Q_i_arr(1)
         P_perp = sl_in%P_perp_arr(1)
         P_tor = sl_in%P_tor_arr(1)
         tau = sl_in%tau_arr(1)
         D_norm = sl_in%D_norm_arr(1)
         c_beta = sl_in%c_beta_arr(1)
         tauk = sl_in%Qconv_arr(1)
         iota_e = Q_e / (Q_e - Q_i)

         g_tmp = g_in
         tmp_delta=riccati_f()
c        de-normalise delta by lu^(1/3)
         det_val=tmp_delta*(sl_in%lu_arr(1)**(1.0/3.0))

c        return Deltaprime - delta(Q)
         dispersion_det = sl_in%Re_dp_arr(1) - det_val

c --- coupled-surface branch (2 or 3 surfaces)
      ELSEIF ((msing_max == 2) .OR. (msing_max == 3)) THEN
         ALLOCATE(delta_Q(msing_max,msing_max))
         delta_Q=(0.0,0.0)
         DO k=1,msing_max
c           set module-level variables for this surface
            Q_e = sl_in%Q_e_arr(k)
            Q_i = sl_in%Q_i_arr(k)
            P_perp = sl_in%P_perp_arr(k)
            P_tor = sl_in%P_tor_arr(k)
            tau = sl_in%tau_arr(k)
            D_norm = sl_in%D_norm_arr(k)
            c_beta = sl_in%c_beta_arr(k)
            tauk = sl_in%Qconv_arr(k)
            iota_e = Q_e / (Q_e - Q_i)

c           evaluate riccati_f at rescaled growth rate, de-normalise
            g_tmp = (g_in*sl_in%Qconv_arr(1))/tauk ! sets module-level g_tmp to SCALED value
            delta_Q(k,k)=riccati_f()
            delta_Q(k,k)=delta_Q(k,k)*sl_in%lu_arr(k)**(1.0/3.0)
         END DO

c        compute det(dp_matrix - delta_Q)
         result_matrix = sl_in%dp_matrix - delta_Q
         CALL calc_determinant(result_matrix, msing_max, det_val,
     $        det_status)
         IF (det_status /= 0) THEN
            WRITE(*,*) 'ERROR: calc_determinant unsupported rank=',
     $           msing_max
            STOP 'dispersion_det: calc_determinant failed'
         END IF
         dispersion_det = det_val
      ELSE
         WRITE(*,*) "Error: no support for msing > 3"
         STOP "dispersion_det: unsupported n_k"
      END IF
      END FUNCTION dispersion_det

c-----------------------------------------------------------------------
c     dispersion_AMR: hash-based adaptive mesh refinement scanner.\nc     
c     Scans a 2D complex-Q grid for zeros of the dispersion relation
c     D(Q) using adaptive refinement.  A coarse grid is evaluated
c     first (two-pass: nodes then cells), then cells that span a zero
c     crossing in Re(D) or Im(D) are subdivided.
c
c     Point deduplication uses a spatial hash table (HASH_SZ buckets,
c     chained) so that midpoints shared between neighbouring cells
c     are evaluated only once.
c-----------------------------------------------------------------------
      SUBROUTINE dispersion_AMR(n_k,sl_in,msing_max,
     $                          scan_width,Q_num,AMR_passes,
     $                          coupling_flag)
c     DEPRECATED: use dispersion_AMR_v2 instead.
c     This v1 hash-based scanner is retained for
c     backwards compatibility only.

c --- arguments
      INTEGER, INTENT(IN)  :: n_k           ! number of rational surfaces
      INTEGER, INTENT(IN)  :: msing_max     ! max surfaces for coupling
      INTEGER, INTENT(IN)  :: Q_num         ! grid points per axis
      INTEGER, INTENT(IN)  :: AMR_passes    ! refinement passes
      REAL(r8), INTENT(IN) :: scan_width    ! half-width of Re/Im scan window
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)  :: coupling_flag ! coupled dispersion_det?
c --- cell storage
      INTEGER, ALLOCATABLE :: cells(:,:)       ! (4, max) corner indices per cell
      INTEGER, ALLOCATABLE :: new_cells(:,:)   ! scratch for next level
c --- loop / index variables
      INTEGER :: n_cells, n_new_cells, i, j
      INTEGER :: c_idx, pass
      INTEGER :: idx_TL, idx_TR, idx_BL, idx_BR   ! corner indices
      INTEGER :: idx_TM, idx_BM, idx_LM, idx_RM, idx_MM  ! midpoint indices
c --- scan workspace
      REAL(r8) :: r_min, r_max, i_min, i_max   ! min/max Re/Im across corners
      REAL(r8) :: ing_step                      ! coarse grid spacing
      REAL(r8) :: ing_coarse, iing_coarse       ! Re/Im coords for coarse node
      LOGICAL  :: cross_real, cross_imag        ! zero-crossing flags
      COMPLEX(r8) :: q_curr                     ! current evaluation point
      INTEGER, ALLOCATABLE :: coarse_indices(:,:)  ! (Q_num,Q_num) node index map

c --- 1. initialise hash table and point/cell storage
      IF (ALLOCATED(Q_store)) DEALLOCATE(Q_store)
      IF (ALLOCATED(D_store)) DEALLOCATE(D_store)
      IF (ALLOCATED(hash_head)) DEALLOCATE(hash_head)
      IF (ALLOCATED(hash_next)) DEALLOCATE(hash_next)

      ALLOCATE(Q_store(MAX_PTS), D_store(MAX_PTS))
      ALLOCATE(hash_head(HASH_SZ), hash_next(MAX_PTS))
      ALLOCATE(cells(4, 200000), new_cells(4, 200000))

      hash_head = 0
      hash_next = 0
      n_pts = 0
      n_cells = 0

c --- 2. build initial coarse grid (two-pass method)
c     Pass 1 computes and hashes every node; Pass 2 stitches cells
c     from the stored indices -- no floating-point comparison needed.
      ALLOCATE(coarse_indices(Q_num, Q_num))
      ing_step = (2.0*scan_width) / (Q_num - 1)

c     Pass 1: compute all grid nodes and store their hash indices
      DO i = 1, Q_num
         DO j = 1, Q_num
             ing_coarse = -scan_width + (i - 1) * ing_step
             iing_coarse = -scan_width + (j - 1) * ing_step
             q_curr = CMPLX(ing_coarse, iing_coarse)

             CALL get_or_compute(q_curr, coarse_indices(i,j), n_k,
     $                        sl_in, msing_max, coupling_flag)
         END DO
      END DO

c     Pass 2: stitch cells from the stored integer indices
      DO i = 1, Q_num - 1
         DO j = 1, Q_num - 1
             n_cells = n_cells + 1
             cells(1, n_cells) = coarse_indices(i, j)     ! TL
             cells(2, n_cells) = coarse_indices(i+1, j)   ! TR
             cells(3, n_cells) = coarse_indices(i, j+1)   ! BL
             cells(4, n_cells) = coarse_indices(i+1, j+1) ! BR
         END DO
      END DO
      DEALLOCATE(coarse_indices)
  
c --- 3. refinement passes: subdivide cells with zero crossings
      DO pass = 1, AMR_PASSES
          WRITE(*,'(A,I2,A,I6,A)') '   > Pass ', pass,
     $         ': Checking ', n_cells, ' cells...'
          n_new_cells = 0

          DO c_idx = 1, n_cells
              idx_TL = cells(1, c_idx)
              idx_TR = cells(2, c_idx)
              idx_BL = cells(3, c_idx)
              idx_BR = cells(4, c_idx)

c             check for sign change in Re(D) across cell corners
              r_min = MIN(REAL(D_store(idx_TL)),
     $                         REAL(D_store(idx_TR)),
     $                         REAL(D_store(idx_BL)),
     $                         REAL(D_store(idx_BR)))
              r_max = MAX(REAL(D_store(idx_TL)),
     $                         REAL(D_store(idx_TR)),
     $                         REAL(D_store(idx_BL)),
     $                         REAL(D_store(idx_BR)))
              cross_real = (r_min * r_max <= 0.0d0)

c             check for sign change in Im(D) across cell corners
              i_min = MIN(AIMAG(D_store(idx_TL)),
     $                         AIMAG(D_store(idx_TR)),
     $                         AIMAG(D_store(idx_BL)),
     $                         AIMAG(D_store(idx_BR)))
              i_max = MAX(AIMAG(D_store(idx_TL)),
     $                         AIMAG(D_store(idx_TR)),
     $                         AIMAG(D_store(idx_BL)),
     $                         AIMAG(D_store(idx_BR)))
              cross_imag = (i_min * i_max <= 0.0d0)
  
              IF (cross_real .OR. cross_imag) THEN
c                 refine: compute 5 midpoints, create 4 sub-cells
                  
                  ! Top-Mid
                  q_curr = 0.5d0*(Q_store(idx_TL)+Q_store(idx_TR))
                  CALL get_or_compute(q_curr, idx_TM,n_k,
     $              sl_in,msing_max,coupling_flag)

                  ! Bot-Mid
                  q_curr = 0.5d0*(Q_store(idx_BL)+Q_store(idx_BR))
                  CALL get_or_compute(q_curr, idx_BM,n_k,
     $              sl_in,msing_max,coupling_flag)

                  ! Left-Mid
                  q_curr = 0.5d0*(Q_store(idx_TL)+Q_store(idx_BL))
                  CALL get_or_compute(q_curr, idx_LM,n_k,
     $              sl_in,msing_max,coupling_flag)

                  ! Right-Mid
                  q_curr = 0.5d0*(Q_store(idx_TR)+Q_store(idx_BR))
                  CALL get_or_compute(q_curr, idx_RM,n_k,
     $              sl_in,msing_max,coupling_flag)

                  ! Center
                  q_curr = 0.5d0*(Q_store(idx_TL)+Q_store(idx_BR))
                  CALL get_or_compute(q_curr, idx_MM,n_k,
     $              sl_in,msing_max,coupling_flag)
                  
                  ! Create 4 sub-cells (TL, TR, BL, BR quadrants)
                  n_new_cells = n_new_cells + 1
                  new_cells(1, n_new_cells) = idx_TL
                  new_cells(2, n_new_cells) = idx_TM
                  new_cells(3, n_new_cells) = idx_LM
                  new_cells(4, n_new_cells) = idx_MM
                  
                  ! Sub 2 (Top-Right)
                  n_new_cells = n_new_cells + 1
                  new_cells(1, n_new_cells) = idx_TM
                  new_cells(2, n_new_cells) = idx_TR
                  new_cells(3, n_new_cells) = idx_MM
                  new_cells(4, n_new_cells) = idx_RM
                  
                  ! Sub 3 (Bot-Left)
                  n_new_cells = n_new_cells + 1
                  new_cells(1, n_new_cells) = idx_LM
                  new_cells(2, n_new_cells) = idx_MM
                  new_cells(3, n_new_cells) = idx_BL
                  new_cells(4, n_new_cells) = idx_BM
                  
                  ! Sub 4 (Bot-Right)
                  n_new_cells = n_new_cells + 1
                  new_cells(1, n_new_cells) = idx_MM
                  new_cells(2, n_new_cells) = idx_RM
                  new_cells(3, n_new_cells) = idx_BM
                  new_cells(4, n_new_cells) = idx_BR
                  
              ELSE
                  ! No refinement needed, keep original cell
                  n_new_cells = n_new_cells + 1
                  new_cells(:, n_new_cells) = cells(:, c_idx)
              END IF
          END DO
          
c --- swap arrays for next refinement pass
          n_cells = n_new_cells
          cells(:, 1:n_cells) = new_cells(:, 1:n_cells)

      END DO
      DEALLOCATE(cells, new_cells)
      WRITE(*,*) "AMR Scan Complete. Total Points:", n_pts
      RETURN
      END SUBROUTINE dispersion_AMR
          
c-----------------------------------------------------------------------
c     get_or_compute: hash-based point lookup for dispersion_AMR v1.
c     If q_in is already in the hash table, return its index.
c     Otherwise, evaluate the dispersion relation at q_in, store
c     the result, and insert into the hash table.
c
c     Uses 64-bit arithmetic internally to avoid integer overflow
c     in the hash function.
c-----------------------------------------------------------------------
      SUBROUTINE get_or_compute(q_in,idx_out,n_k,sl_in,msing_max,
     $                          coupling_flag)

c --- arguments
      COMPLEX(r8), INTENT(IN) :: q_in        ! complex-Q evaluation point
      INTEGER, INTENT(OUT) :: idx_out        ! returned point index
      INTEGER, INTENT(IN) :: n_k             ! number of surfaces
      INTEGER, INTENT(IN) :: msing_max       ! max surfaces to include
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN) :: coupling_flag   ! use coupled dispersion_det?
c --- locals
      INTEGER :: h               ! hash bucket index
      INTEGER :: curr             ! linked-list traversal index
      COMPLEX(r8) :: delta_val    ! computed dispersion result
      INTEGER(8) :: ix8, iy8, h8  ! 64-bit intermediates for hash

c --- 1. compute hash from quantised Re/Im coordinates (64-bit safe)
      ix8 = NINT(REAL(q_in) * HASH_SCALE, KIND=8)
      iy8 = NINT(AIMAG(q_in) * HASH_SCALE, KIND=8)
      h8 = MOD(ABS(ix8 * 73856093_8 + iy8 * 19349663_8),
     $         INT(HASH_SZ, 8)) + 1_8
      h = INT(h8)

      IF (h < 1 .OR. h > HASH_SZ) THEN
         WRITE(*,*) "HASH ERROR: h=", h, " q_in=", q_in
         STOP "get_or_compute: hash out of bounds"
      END IF

c --- 2. search hash chain for existing point
      curr = hash_head(h)
      DO WHILE (curr /= 0)
          IF (ABS(Q_store(curr) - q_in) < 1.0d-8) THEN
              idx_out = curr
              RETURN
          END IF
          curr = hash_next(curr)
      END DO

c --- 3. point not found: evaluate dispersion relation and store
      n_pts = n_pts + 1
      IF (n_pts > MAX_PTS) THEN
          WRITE(*,*) "ERROR: AMR exceeded MAX_PTS"
          STOP "get_or_compute: MAX_PTS exceeded"
      END IF

      idx_out = n_pts
      Q_store(idx_out) = q_in

      IF (coupling_flag) THEN
           g_tmp = q_in
           delta_val = dispersion_det(g_tmp, n_k, sl_in, msing_max)
      ELSE
           g_tmp = q_in
           delta_val = riccati_f()
           delta_val = delta_val - delta_eff
      END IF
      D_store(idx_out) = delta_val

c --- 4. insert into hash chain (prepend)
      hash_next(idx_out) = hash_head(h)
      hash_head(h) = idx_out

      END SUBROUTINE get_or_compute

c-----------------------------------------------------------------------
c     get_or_compute_v2: hash-cached dispersion evaluation for AMR v2.
c     Identical to get_or_compute but applies the ifac (imaginary-unit)
c     Wick rotation that compute_delta_sub uses:  g_tmp = q_in * ifac.
c-----------------------------------------------------------------------
      SUBROUTINE get_or_compute_v2(q_in, idx_out, n_k, sl_in,
     $                              msing_max, coupling_flag)

      IMPLICIT NONE

c --- arguments
      COMPLEX(r8), INTENT(IN)  :: q_in
      INTEGER, INTENT(OUT)     :: idx_out
      INTEGER, INTENT(IN)      :: n_k, msing_max
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)      :: coupling_flag

c --- locals
      INTEGER     :: h, curr
      COMPLEX(r8) :: delta_val
      INTEGER(8)  :: ix8, iy8, h8

c --- 1. compute hash bucket
      ix8 = NINT(REAL(q_in) * HASH_SCALE, KIND=8)
      iy8 = NINT(AIMAG(q_in) * HASH_SCALE, KIND=8)
      h8  = MOD(ABS(ix8 * 73856093_8 + iy8 * 19349663_8),
     $          INT(HASH_SZ, 8)) + 1_8
      h   = INT(h8)

c --- 2. search hash chain for existing point
      curr = hash_head(h)
      DO WHILE (curr /= 0)
          IF (ABS(Q_store(curr) - q_in) < 1.0d-8) THEN
              idx_out = curr
              RETURN
          END IF
          curr = hash_next(curr)
      END DO

c --- 3. not found: evaluate with ifac rotation and store
      n_pts = n_pts + 1
      IF (n_pts > MAX_PTS) THEN
          WRITE(*,*) 'ERROR: AMR v2 cache exceeded MAX_PTS'
          STOP 'get_or_compute_v2: MAX_PTS exceeded'
      END IF

      idx_out = n_pts
      Q_store(idx_out) = q_in

      IF (coupling_flag) THEN
          g_tmp = q_in * ifac
          delta_val = dispersion_det(g_tmp, n_k, sl_in,
     $                               msing_max)
      ELSE
          g_tmp = q_in * ifac
          delta_val = riccati_f()
          delta_val = delta_val - delta_eff
      END IF
      D_store(idx_out) = delta_val

c --- 4. insert into hash chain (prepend)
      hash_next(idx_out) = hash_head(h)
      hash_head(h) = idx_out

      END SUBROUTINE get_or_compute_v2

c-----------------------------------------------------------------------
c     dispersion_AMR_v2: cell-based adaptive mesh refinement scanner.
c     Unlike v1 (hash-based point deduplication), v2 stores complete
c     cells (TYPE amr_cell_type), each carrying 4 corner Q- and D-
c     values.  Refinement subdivides cells that contain a zero in
c     Re(D) or Im(D) and re-evaluates the dispersion relation at the
c     5 new midpoints.
c
c     All dispersion evaluations go through get_or_compute_v2, which
c     caches results in Q_store / D_store via a hash table.  This
c     eliminates redundant evaluations for shared corners (initial grid)
c     and shared edge-midpoints (refinement).  At completion, Q_store
c     and D_store are trimmed to n_pts unique output points.
c
c-----------------------------------------------------------------------
      SUBROUTINE dispersion_AMR_v2(n_k, sl_in, msing_max,
     $                             scan_width, Q_num, AMR_passes,
     $                             coupling_flag)

      IMPLICIT NONE

c --- arguments
      INTEGER, INTENT(IN)  :: n_k           ! number of rational surfaces
      INTEGER, INTENT(IN)  :: msing_max     ! max surfaces for coupling
      INTEGER, INTENT(IN)  :: Q_num         ! grid points per axis
      INTEGER, INTENT(IN)  :: AMR_passes    ! refinement passes
      REAL(r8), INTENT(IN) :: scan_width    ! half-width of Re/Im scan window
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)  :: coupling_flag ! coupled dispersion_det?
c --- locals
      TYPE(amr_cell_type), ALLOCATABLE :: new_cells(:)
      TYPE(amr_cell_type), ALLOCATABLE :: swap_tmp(:)  ! for pointer swap
      INTEGER     :: i, j, c, corner, pass   ! loop counters
      REAL(r8)    :: step                     ! grid spacing
      REAL(r8)    :: x, y                     ! real / imag grid coords
      LOGICAL     :: cross_real, cross_imag   ! zero-crossing flags
      INTEGER     :: n_new_cells              ! count during refinement
      INTEGER     :: cells_to_refine          ! cells flagged per pass
      INTEGER     :: cells_kept               ! cells kept per pass
      INTEGER     :: idx_tmp                  ! hash-cache index
      COMPLEX(r8), ALLOCATABLE :: temp_Q(:)  ! for trimming output
      COMPLEX(r8), ALLOCATABLE :: temp_D(:)  ! for trimming output
      
c --- 1. initialise cell storage and hash cache

      IF (ALLOCATED(amr_cells)) DEALLOCATE(amr_cells)
      IF (ALLOCATED(Q_store))   DEALLOCATE(Q_store)
      IF (ALLOCATED(D_store))   DEALLOCATE(D_store)
      IF (ALLOCATED(hash_head)) DEALLOCATE(hash_head)
      IF (ALLOCATED(hash_next)) DEALLOCATE(hash_next)

      ALLOCATE(amr_cells(MAX_CELLS))
      ALLOCATE(new_cells(MAX_CELLS))
      ALLOCATE(Q_store(MAX_PTS))
      ALLOCATE(D_store(MAX_PTS))
      ALLOCATE(hash_head(HASH_SZ))
      ALLOCATE(hash_next(MAX_PTS))

      hash_head = 0
      hash_next = 0
      n_pts     = 0
      n_amr_cells = 0
      step = (2.0d0 * scan_width) / DBLE(Q_num - 1)
      
c --- 2. build initial coarse grid of (Q_num-1)^2 cells

      DO i = 1, Q_num - 1
          DO j = 1, Q_num - 1
              x = -scan_width + DBLE(i-1) * step
              y = -scan_width + DBLE(j-1) * step

              n_amr_cells = n_amr_cells + 1

              IF (n_amr_cells > MAX_CELLS) THEN
                  WRITE(*,*) 'ERROR: Exceeded MAX_CELLS in init'
                  STOP 'dispersion_AMR_v2: MAX_CELLS in init'
              END IF

c             corner order: BL=1, BR=2, TL=3, TR=4
              amr_cells(n_amr_cells)%Q(1) = CMPLX(x, y, KIND=r8)
              amr_cells(n_amr_cells)%Q(2) = CMPLX(x+step, y, KIND=r8)
              amr_cells(n_amr_cells)%Q(3) = CMPLX(x, y+step, KIND=r8)
              amr_cells(n_amr_cells)%Q(4) = CMPLX(x+step, y+step,
     $                                            KIND=r8)

c             evaluate dispersion at each corner (hash-cached)
              DO corner = 1, 4
                  CALL get_or_compute_v2(
     $                amr_cells(n_amr_cells)%Q(corner),
     $                idx_tmp, n_k, sl_in, msing_max,
     $                coupling_flag)
                  amr_cells(n_amr_cells)%D(corner) =
     $                D_store(idx_tmp)
              END DO

              amr_cells(n_amr_cells)%needs_refine = .FALSE.
          END DO
      END DO

c --- 3. refinement passes: subdivide cells with zero crossings
      DO pass = 1, AMR_passes
          WRITE(*,'(A,I2,A,I7,A)') '   Pass ', pass,
     $         ': Processing ', n_amr_cells, ' cells'

c         flag cells that span a zero in Re(D) or Im(D)
          cells_to_refine = 0
          DO c = 1, n_amr_cells
              CALL check_cell_crossing_sub(amr_cells(c),
     $                                     cross_real, cross_imag)
              amr_cells(c)%needs_refine = (cross_real .OR. cross_imag)
              IF (amr_cells(c)%needs_refine) THEN
                  cells_to_refine = cells_to_refine + 1
              END IF
          END DO

c         build new cell list: subdivide flagged, keep the rest
          n_new_cells = 0
          cells_kept = 0

          DO c = 1, n_amr_cells
              IF (amr_cells(c)%needs_refine) THEN
                  CALL subdivide_cell_sub(amr_cells(c),
     $                 new_cells, n_new_cells, MAX_CELLS,
     $                 n_k, sl_in, msing_max, coupling_flag)
              ELSE
                  n_new_cells = n_new_cells + 1
                  IF (n_new_cells > MAX_CELLS) THEN
                      WRITE(*,*) 'ERROR: Exceeded MAX_CELLS in refine'
                      STOP 'dispersion_AMR_v2: MAX_CELLS in refine'
                  END IF
                  new_cells(n_new_cells) = amr_cells(c)
                  cells_kept = cells_kept + 1
              END IF
          END DO

c         swap arrays for next pass (pointer swap, no element copy)
          CALL MOVE_ALLOC(new_cells, swap_tmp)
          CALL MOVE_ALLOC(amr_cells, new_cells)   ! old amr_cells becomes new_cells
          CALL MOVE_ALLOC(swap_tmp, amr_cells)     ! filled array becomes amr_cells
          n_amr_cells = n_new_cells

      END DO

c --- 4. output: Q_store/D_store already populated by hash cache.
c     Trim to exact size n_pts and deallocate hash infrastructure.

      ALLOCATE(temp_Q(n_pts))
      ALLOCATE(temp_D(n_pts))
      temp_Q(1:n_pts) = Q_store(1:n_pts)
      temp_D(1:n_pts) = D_store(1:n_pts)
      CALL MOVE_ALLOC(temp_Q, Q_store)
      CALL MOVE_ALLOC(temp_D, D_store)

      IF (ALLOCATED(hash_head)) DEALLOCATE(hash_head)
      IF (ALLOCATED(hash_next)) DEALLOCATE(hash_next)
      DEALLOCATE(new_cells)
c     keep amr_cells allocated for potential post-run inspection

      WRITE(*,*) 'AMR v2 Complete. Unique output points:', n_pts
      WRITE(*,'(A,2ES14.6)') '   D_store checksum (Re,Im):',
     $   SUM(REAL(D_store(1:n_pts))),
     $   SUM(AIMAG(D_store(1:n_pts)))
      WRITE(*,'(A,2ES14.6)') '   D_store(1) sample:',
     $   REAL(D_store(1)), AIMAG(D_store(1))

      RETURN
      END SUBROUTINE dispersion_AMR_v2

c-----------------------------------------------------------------------
c     compute_delta_sub: evaluate the dispersion relation at a single
c     complex-Q point for dispersion_AMR_v2.  Multiplies q_in by ifac
c     before passing to the Riccati solver or coupled-surface
c     determinant routine.
c-----------------------------------------------------------------------
      SUBROUTINE compute_delta_sub(q_in, n_k, sl_in, msing_max,
     $                             coupling_flag, delta_out)

      IMPLICIT NONE

c --- arguments
      COMPLEX(r8), INTENT(IN)  :: q_in          ! complex-Q evaluation point
      INTEGER, INTENT(IN)      :: n_k           ! number of surfaces
      INTEGER, INTENT(IN)      :: msing_max     ! max surfaces for coupling
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)      :: coupling_flag ! use coupled det?
      COMPLEX(r8), INTENT(OUT) :: delta_out     ! dispersion result

      IF (coupling_flag) THEN
          g_tmp = q_in*ifac
          delta_out = dispersion_det(g_tmp, n_k, sl_in, msing_max)
      ELSE
          g_tmp = q_in*ifac
          delta_out = riccati_f()
          delta_out = delta_out - delta_eff
      END IF

      RETURN
      END SUBROUTINE compute_delta_sub


c-----------------------------------------------------------------------
c     check_cell_crossing_sub: test whether a cell’s 4 corner D-values
c     span a zero crossing in Re(D) and/or Im(D).  Used by
c     dispersion_AMR_v2 to decide which cells to refine.
c-----------------------------------------------------------------------
      SUBROUTINE check_cell_crossing_sub(cell, cross_real, cross_imag)

      IMPLICIT NONE

      TYPE(amr_cell_type), INTENT(IN) :: cell
      LOGICAL, INTENT(OUT) :: cross_real, cross_imag

      REAL(r8) :: r_vals(4), i_vals(4)  ! corner Re/Im values
      REAL(r8) :: r_min, r_max, i_min, i_max
      INTEGER  :: k

      DO k = 1, 4
          r_vals(k) = REAL(cell%D(k), KIND=r8)
          i_vals(k) = AIMAG(cell%D(k))
      END DO

      r_min = MINVAL(r_vals)
      r_max = MAXVAL(r_vals)
      cross_real = (r_min * r_max <= 0.0d0)

      i_min = MINVAL(i_vals)
      i_max = MAXVAL(i_vals)
      cross_imag = (i_min * i_max <= 0.0d0)

      RETURN
      END SUBROUTINE check_cell_crossing_sub


c-----------------------------------------------------------------------
c     subdivide_cell_sub: split a parent cell into 4 child cells by
c     computing 5 midpoints (bottom-mid, top-mid, left-mid, right-mid,
c     centre) and evaluating the dispersion relation at each.  The 4
c     resulting child cells are appended to new_cells.
c-----------------------------------------------------------------------
      SUBROUTINE subdivide_cell_sub(parent, new_cells, n_new,
     $                              max_cells, n_k, sl_in,
     $                              msing_max, coupling_flag)

      IMPLICIT NONE

c --- arguments
      TYPE(amr_cell_type), INTENT(IN)    :: parent
      TYPE(amr_cell_type), INTENT(INOUT) :: new_cells(*)
      INTEGER, INTENT(INOUT) :: n_new       ! running count of new cells
      INTEGER, INTENT(IN)    :: max_cells   ! capacity of new_cells
      INTEGER, INTENT(IN)    :: n_k
      INTEGER, INTENT(IN)    :: msing_max
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)    :: coupling_flag
c --- corner coordinates and D-values from parent
      COMPLEX(r8) :: q_bl, q_br, q_tl, q_tr
      COMPLEX(r8) :: d_bl, d_br, d_tl, d_tr
c --- midpoint coordinates and D-values (cached via hash)
      COMPLEX(r8) :: q_bm, q_tm, q_lm, q_rm, q_mm
      COMPLEX(r8) :: d_bm, d_tm, d_lm, d_rm, d_mm
      INTEGER     :: idx_tmp                  ! hash-cache index
      
c --- extract parent corners (BL=1, BR=2, TL=3, TR=4)
      q_bl = parent%Q(1)
      q_br = parent%Q(2)
      q_tl = parent%Q(3)
      q_tr = parent%Q(4)

      d_bl = parent%D(1)
      d_br = parent%D(2)
      d_tl = parent%D(3)
      d_tr = parent%D(4)

c --- compute 5 midpoint coordinates
      q_bm = 0.5d0 * (q_bl + q_br)
      q_tm = 0.5d0 * (q_tl + q_tr)
      q_lm = 0.5d0 * (q_bl + q_tl)
      q_rm = 0.5d0 * (q_br + q_tr)
      q_mm = 0.25d0 * (q_bl + q_br + q_tl + q_tr)

c --- evaluate dispersion at new midpoints (hash-cached)
      CALL get_or_compute_v2(q_bm, idx_tmp,
     $     n_k, sl_in, msing_max, coupling_flag)
      d_bm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_tm, idx_tmp,
     $     n_k, sl_in, msing_max, coupling_flag)
      d_tm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_lm, idx_tmp,
     $     n_k, sl_in, msing_max, coupling_flag)
      d_lm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_rm, idx_tmp,
     $     n_k, sl_in, msing_max, coupling_flag)
      d_rm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_mm, idx_tmp,
     $     n_k, sl_in, msing_max, coupling_flag)
      d_mm = D_store(idx_tmp)

c --- check space for 4 new cells
      IF (n_new + 4 > max_cells) THEN
          WRITE(*,*) 'ERROR: Would exceed MAX_CELLS in subdivide'
          STOP 'subdivide_cell_sub: MAX_CELLS exceeded'
      END IF

c --- child 1: bottom-left quadrant (BL, BM, LM, MM)
      n_new = n_new + 1
      new_cells(n_new)%Q(1) = q_bl
      new_cells(n_new)%Q(2) = q_bm
      new_cells(n_new)%Q(3) = q_lm
      new_cells(n_new)%Q(4) = q_mm
      new_cells(n_new)%D(1) = d_bl
      new_cells(n_new)%D(2) = d_bm
      new_cells(n_new)%D(3) = d_lm
      new_cells(n_new)%D(4) = d_mm
      new_cells(n_new)%needs_refine = .FALSE.
      
c --- child 2: bottom-right quadrant (BM, BR, MM, RM)
      n_new = n_new + 1
      new_cells(n_new)%Q(1) = q_bm
      new_cells(n_new)%Q(2) = q_br
      new_cells(n_new)%Q(3) = q_mm
      new_cells(n_new)%Q(4) = q_rm
      new_cells(n_new)%D(1) = d_bm
      new_cells(n_new)%D(2) = d_br
      new_cells(n_new)%D(3) = d_mm
      new_cells(n_new)%D(4) = d_rm
      new_cells(n_new)%needs_refine = .FALSE.
      
c --- child 3: top-left quadrant (LM, MM, TL, TM)
      n_new = n_new + 1
      new_cells(n_new)%Q(1) = q_lm
      new_cells(n_new)%Q(2) = q_mm
      new_cells(n_new)%Q(3) = q_tl
      new_cells(n_new)%Q(4) = q_tm
      new_cells(n_new)%D(1) = d_lm
      new_cells(n_new)%D(2) = d_mm
      new_cells(n_new)%D(3) = d_tl
      new_cells(n_new)%D(4) = d_tm
      new_cells(n_new)%needs_refine = .FALSE.
      
c --- child 4: top-right quadrant (MM, RM, TM, TR)
      n_new = n_new + 1
      new_cells(n_new)%Q(1) = q_mm
      new_cells(n_new)%Q(2) = q_rm
      new_cells(n_new)%Q(3) = q_tm
      new_cells(n_new)%Q(4) = q_tr
      new_cells(n_new)%D(1) = d_mm
      new_cells(n_new)%D(2) = d_rm
      new_cells(n_new)%D(3) = d_tm
      new_cells(n_new)%D(4) = d_tr
      new_cells(n_new)%needs_refine = .FALSE.
      
      RETURN
      END SUBROUTINE subdivide_cell_sub


c-----------------------------------------------------------------------
c     flatten_cells_to_points_sub: extract unique (Q, D) points from
c     the cell array into the module-level Q_store / D_store arrays.
c     Uses a brute-force O(n^2) duplicate check which is acceptable
c     for moderate cell counts; could be replaced by a hash set for
c     very large scans.
c
c     n_total_corners = num_cells*4; guarded against MAX_PTS overflow.
c-----------------------------------------------------------------------
      SUBROUTINE flatten_cells_to_points_sub(num_cells)

      IMPLICIT NONE

c --- arguments
      INTEGER, INTENT(IN) :: num_cells       ! number of cells to flatten
c --- locals
      INTEGER :: c, corner, i, idx
      INTEGER :: n_total_corners              ! = num_cells * 4
      COMPLEX(r8), ALLOCATABLE :: temp_Q(:)  ! all corner Q-values
      COMPLEX(r8), ALLOCATABLE :: temp_D(:)  ! all corner D-values
      INTEGER, ALLOCATABLE :: sort_idx(:)    ! sort permutation
      REAL(r8) :: tol                        ! duplicate tolerance

      tol = 1.0d-10
      n_total_corners = num_cells * 4

      IF (n_total_corners > MAX_PTS) THEN
          WRITE(*,*) 'ERROR: n_total_corners=', n_total_corners,
     $               ' exceeds MAX_PTS=', MAX_PTS
          STOP 'flatten_cells_to_points_sub: MAX_PTS exceeded'
      END IF

      IF (num_cells <= 0) THEN
          WRITE(*,*) 'ERROR: No cells to flatten'
          n_pts = 0
          RETURN
      END IF

c --- gather all corners from cells
      ALLOCATE(temp_Q(n_total_corners))
      ALLOCATE(temp_D(n_total_corners))
      ALLOCATE(sort_idx(n_total_corners))

      idx = 0
      DO c = 1, num_cells
          DO corner = 1, 4
              idx = idx + 1
              temp_Q(idx) = amr_cells(c)%Q(corner)
              temp_D(idx) = amr_cells(c)%D(corner)
              sort_idx(idx) = idx
          END DO
      END DO

c --- sort by (Re(Q), Im(Q)) via quicksort on the index array
      CALL qsort_complex_idx(temp_Q, sort_idx, 1, n_total_corners)

c --- linear scan to count unique points (sorted order)
      n_pts = 1
      DO i = 2, n_total_corners
          IF (ABS(temp_Q(sort_idx(i)) - temp_Q(sort_idx(i-1)))
     $        >= tol) THEN
              n_pts = n_pts + 1
          END IF
      END DO

c --- copy unique points to module-level output arrays
      IF (ALLOCATED(Q_store)) DEALLOCATE(Q_store)
      IF (ALLOCATED(D_store)) DEALLOCATE(D_store)
      ALLOCATE(Q_store(n_pts))
      ALLOCATE(D_store(n_pts))

      idx = 1
      Q_store(1) = temp_Q(sort_idx(1))
      D_store(1) = temp_D(sort_idx(1))
      DO i = 2, n_total_corners
          IF (ABS(temp_Q(sort_idx(i)) - temp_Q(sort_idx(i-1)))
     $        >= tol) THEN
              idx = idx + 1
              Q_store(idx) = temp_Q(sort_idx(i))
              D_store(idx) = temp_D(sort_idx(i))
          END IF
      END DO

      DEALLOCATE(temp_Q, temp_D, sort_idx)

      RETURN
      END SUBROUTINE flatten_cells_to_points_sub

c-----------------------------------------------------------------------
c     qsort_complex_idx: in-place quicksort of an index array by
c     the complex keys (Re then Im).  Operates on sort_idx so the
c     Q/D data arrays remain untouched.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE qsort_complex_idx(keys, idx, lo, hi)

      IMPLICIT NONE

      COMPLEX(r8), INTENT(IN)    :: keys(:)
      INTEGER,     INTENT(INOUT) :: idx(:)
      INTEGER,     INTENT(IN)    :: lo, hi

      INTEGER :: i, j, pivot_idx, tmp
      REAL(r8) :: p_re, p_im, k_re, k_im

      IF (lo >= hi) RETURN

c     median-of-three pivot selection
      pivot_idx = idx((lo + hi) / 2)
      p_re = REAL(keys(pivot_idx), KIND=r8)
      p_im = AIMAG(keys(pivot_idx))

      i = lo
      j = hi
      DO WHILE (i <= j)
c         advance i while keys(idx(i)) < pivot
          k_re = REAL(keys(idx(i)), KIND=r8)
          k_im = AIMAG(keys(idx(i)))
          DO WHILE (k_re < p_re .OR.
     $             (k_re == p_re .AND. k_im < p_im))
              i = i + 1
              k_re = REAL(keys(idx(i)), KIND=r8)
              k_im = AIMAG(keys(idx(i)))
          END DO
c         retreat j while keys(idx(j)) > pivot
          k_re = REAL(keys(idx(j)), KIND=r8)
          k_im = AIMAG(keys(idx(j)))
          DO WHILE (k_re > p_re .OR.
     $             (k_re == p_re .AND. k_im > p_im))
              j = j - 1
              k_re = REAL(keys(idx(j)), KIND=r8)
              k_im = AIMAG(keys(idx(j)))
          END DO
c         swap if pointers haven't crossed
          IF (i <= j) THEN
              tmp    = idx(i)
              idx(i) = idx(j)
              idx(j) = tmp
              i = i + 1
              j = j - 1
          END IF
      END DO

c     recurse on partitions
      IF (lo < j) CALL qsort_complex_idx(keys, idx, lo, j)
      IF (i < hi) CALL qsort_complex_idx(keys, idx, i, hi)

      RETURN
      END SUBROUTINE qsort_complex_idx
      END MODULE gslayer_mod