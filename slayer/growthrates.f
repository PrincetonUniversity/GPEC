      MODULE growthrates_mod
c-----------------------------------------------------------------------
c     growthrates_mod: Growth-rate scanning and AMR dispersion solvers.
c
c     Split from gslayer_mod, this module contains every subroutine
c     after gpec_slayer that was formerly inside gslayer.f.  It
c     provides I/O helpers, array utilities, the dispersion
c     determinant, and both AMR scanner variants used by the main
c     SLAYER driver (slayer.f).
c
c     Subprograms:
c       1. output_gamma          - write results to netCDF
c       2. allocate_inputs       - allocate slayer_inputs_type
c       3. allocate_outputs      - allocate slayer_outputs_type
c       4. shrink_array          - trim over-allocated scan arrays
c       5. grow_array            - expand scan arrays dynamically
c       6. calc_determinant      - 2x2 / 3x3 complex determinant
c       7. dispersion_det        - coupled dispersion determinant
c       8. get_or_compute_v2     - hash-cached dispersion eval
c       9. dispersion_AMR_v2     - AMR scan (cell-based storage)
c      10. check_cell_crossing_sub - zero-crossing test
c      11. subdivide_cell_sub    - cell refinement
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
c     subprogram 1. output_gamma.
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
c     subprogram 2. allocate_inputs.
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
c     subprogram 3. allocate_outputs.
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
c     subprogram 4. shrink_array.
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
c     subprogram 5. grow_array.
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
c     subprogram 6. calc_determinant.
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
            detk = CMPLX(0.0_r8, 0.0_r8, KIND=r8)
            status = -1
                    
      END SELECT
      RETURN
      END SUBROUTINE calc_determinant
c-----------------------------------------------------------------------
c     subprogram 7. dispersion_det.
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
         det_val=tmp_delta*
     $      (sl_in%lu_arr(1)**(1.0_r8/3.0_r8))

c        return Deltaprime - delta(Q)
         dispersion_det = sl_in%Re_dp_arr(1) - det_val

c --- coupled-surface branch (2 or 3 surfaces)
      ELSEIF ((msing_max == 2) .OR. (msing_max == 3)) THEN
         ALLOCATE(delta_Q(msing_max,msing_max))
         delta_Q=CMPLX(0.0_r8, 0.0_r8, KIND=r8)
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
c           rescale g_in to this surface's normalisation
            g_tmp = (g_in*sl_in%Qconv_arr(1))/tauk
            delta_Q(k,k)=riccati_f()
            delta_Q(k,k)=delta_Q(k,k)*
     $         sl_in%lu_arr(k)**(1.0_r8/3.0_r8)
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
         DEALLOCATE(delta_Q, result_matrix)
         dispersion_det = det_val
      ELSE
         WRITE(*,*) "Error: no support for msing > 3"
         STOP "dispersion_det: unsupported n_k"
      END IF
      END FUNCTION dispersion_det

c-----------------------------------------------------------------------
c     get_or_compute_v2: hash-cached dispersion evaluation for AMR v2.
c     Applies ifac Wick rotation: g_tmp = q_in * ifac.
c-----------------------------------------------------------------------
      SUBROUTINE get_or_compute_v2(q_in, idx_out, n_k,
     $                              sl_in, msing_max,
     $                              coupling_flag, full)

      IMPLICIT NONE

c --- arguments
      COMPLEX(r8), INTENT(IN)  :: q_in
      INTEGER, INTENT(OUT)     :: idx_out
      INTEGER, INTENT(IN)      :: n_k, msing_max
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN)      :: coupling_flag
      LOGICAL, INTENT(OUT)     :: full

c --- locals
      INTEGER     :: h, curr
      COMPLEX(r8) :: delta_val
      INTEGER(8)  :: ix8, iy8, h8

      full = .FALSE.
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
          n_pts = n_pts - 1
          full = .TRUE.
          idx_out = -1
          RETURN
      END IF

      idx_out = n_pts
      Q_store(idx_out) = q_in

      IF (coupling_flag) THEN
c          dispersion_det sets g_tmp per-surface internally;
c          pass q_in*ifac directly as g_in argument.
          delta_val = dispersion_det(q_in * ifac, n_k, sl_in,
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
      LOGICAL     :: pts_full                 ! MAX_PTS reached flag
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
     $                coupling_flag, pts_full)
                  IF (pts_full) GOTO 800
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
                  CALL subdivide_cell_sub(
     $                 amr_cells(c),
     $                 new_cells, n_new_cells,
     $                 MAX_CELLS, n_k, sl_in,
     $                 msing_max, coupling_flag,
     $                 pts_full)
                  IF (pts_full) GOTO 800
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
      GOTO 810

 800  CONTINUE
      WRITE(*,'(A,I8,A)')
     $   ' WARNING: MAX_PTS (', MAX_PTS,
     $   ') reached during AMR v2.'
      WRITE(*,'(A)')
     $   '   Saving existing results.'

 810  CONTINUE
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
      IF (ALLOCATED(new_cells)) DEALLOCATE(new_cells)
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
      SUBROUTINE subdivide_cell_sub(parent,
     $      new_cells, n_new, max_cells, n_k,
     $      sl_in, msing_max, coupling_flag,
     $      pts_full)

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
      LOGICAL, INTENT(OUT)   :: pts_full    ! MAX_PTS flag
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
      pts_full = .FALSE.
      CALL get_or_compute_v2(q_bm, idx_tmp,
     $     n_k, sl_in, msing_max,
     $     coupling_flag, pts_full)
      IF (pts_full) RETURN
      d_bm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_tm, idx_tmp,
     $     n_k, sl_in, msing_max,
     $     coupling_flag, pts_full)
      IF (pts_full) RETURN
      d_tm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_lm, idx_tmp,
     $     n_k, sl_in, msing_max,
     $     coupling_flag, pts_full)
      IF (pts_full) RETURN
      d_lm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_rm, idx_tmp,
     $     n_k, sl_in, msing_max,
     $     coupling_flag, pts_full)
      IF (pts_full) RETURN
      d_rm = D_store(idx_tmp)
      CALL get_or_compute_v2(q_mm, idx_tmp,
     $     n_k, sl_in, msing_max,
     $     coupling_flag, pts_full)
      IF (pts_full) RETURN
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

      END MODULE growthrates_mod
