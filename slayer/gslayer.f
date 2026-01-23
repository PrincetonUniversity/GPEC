      MODULE gslayer_mod

      USE omp_lib

      USE sglobal_mod, ONLY: out_unit,r8, mu0, m_p, chag, lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,
     $   eta,visc,rho_s,lu,omega_e,omega_i,
     $   delta_n,
     $   Q
      USE delta_mod     
      !, ONLY: riccati,riccati_out,
      !$   parflow_flag,PeOhmOnly_flag

      USE params_mod

      USE layerinputs_mod

      USE slayer_netcdf_mod

      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. gpec_slayer.
c     run slayer to provide b_crit(ising).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_slayer(n_e,t_e,n_i,t_i,zeff,omega,omega_e,
     $   omega_i,qval,sval,bt,rs,R0,mu_i,inpr,mms,nns,ascii_flag,
     $     delta,psi0,jxb,omega_sol,br_th)

      REAL(r8),INTENT(IN) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     qval,sval,bt,rs,R0,zeff,inpr
      INTEGER, INTENT(IN) :: mms,nns,mu_i
      LOGICAL, INTENT(IN) :: ascii_flag
      COMPLEX(r8),INTENT(OUT) :: delta,psi0
      REAL(r8),INTENT(OUT) :: jxb,omega_sol,br_th

      INTEGER :: i,inum
      INTEGER, DIMENSION(1) :: index

      REAL(r8) :: inQ,inQ_e,inQ_i,inpe,inc_beta,inds,intau,inlu
      REAL(r8) :: mrs,nrs,rho,b_l,v_a,Qconv,Q0,delta_n_p,
     $            lbeta,tau_i,tau_h,tau_v
      REAL(r8) :: inQ_min,inQ_max,Q_sol

      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs,jxbl,bal
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal
      CHARACTER(3) :: sn,sm

      parflow_flag=.FALSE.
      PeOhmOnly_flag=.TRUE.
      riccati_out=.FALSE.

      mrs = real(mms,4)
      nrs = real(nns,4)

      ! String representations of the m and n mode numbers
      IF (nns<10) THEN
         WRITE(UNIT=sn,FMT='(I1)') nns
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)') nns
      ENDIF
      IF (mms<10) THEN
         WRITE(UNIT=sm,FMT='(I1)') mms
         sm=ADJUSTL(sm)
      ELSEIF (mms<100) THEN
         WRITE(UNIT=sm,FMT='(I2)') mms
         sm=ADJUSTL(sm)
      ELSE
         WRITE(UNIT=sm,FMT='(I3)') mms
      ENDIF

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

      ! note Q depends on Qconv even IF omega is fixed.
      Q=Qconv*omega
      Q_e=-Qconv*omega_e
      Q_i=-Qconv*omega_i

      ! This is the most critical PARAMETER
      ds=lu**(1.0/3.0)*rho_s/rs        ! conversion based on Cole.

      lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
      c_beta=(lbeta/(1.0+lbeta))**0.5

      delta_n=lu**(1.0/3.0)/rs         ! norm factor for delta primes

      inQ=Q
      inQ_e=Q_e
      inQ_i=Q_i
      inc_beta=c_beta
      inds=ds
      intau=tau
      Q0=Q
c-----------------------------------------------------------------------
c     calculate basic delta, torque, balance, error fields.
c-----------------------------------------------------------------------
      delta_n_p=1e-2
      delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
      psi0=1.0/ABS(delta+delta_n_p)     ! a.u.
      jxb=-AIMAG(1.0/(delta+delta_n_p)) ! a.u.
c-----------------------------------------------------------------------
c     find solutions based on simple torque balance.
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

      ! Scan of rotation
      inQ_max=10.0
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

      ! Write torque balance curves to file for diagnostic purposes
      IF(ascii_flag)THEN
         OPEN(UNIT=out_unit,FILE="gpec_slayer_torque_balance_m"//
     $        TRIM(sm)//"_n"//TRIM(sn)//".OUT",
     $        STATUS="UNKNOWN")
         WRITE(out_unit,'(1x,5(a17))') "inQ","RE(delta)",
     $        "IM(delta)","jxb","bal"
         DO i=0,inum
            WRITE(out_unit,'(1x,5(es17.8e3))')
     $           inQs(i),REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
         ENDDO
         CLOSE(out_unit)
      ENDIF

      ! Identify the threshold from the maximum of the balance PARAMETER
      index=MAXLOC(bal)
      Q_sol=inQs(index(1))
      omega_sol=inQs(index(1))/Qconv
      br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))
      DEALLOCATE(inQs,deltal,jxbl,bal)

      RETURN
      END SUBROUTINE gpec_slayer
c-----------------------------------------------------------------------
c     Subprogram 2. output_gamma
c     Take SLAYER input and output dicts, send to netCDF subroutine
c-----------------------------------------------------------------------
      SUBROUTINE output_gamma(est_gamma_flag,sl_in,sl_out,
     $                        all_deltas_out)

      ! Declarations (include necessary type declarations from original code)
      LOGICAL, INTENT(IN) :: est_gamma_flag
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      TYPE(slayer_outputs_type), INTENT(IN) :: sl_out
      TYPE(deltas_outputs_type), INTENT(IN) :: 
     $                            all_deltas_out(SIZE(sl_in%qval_arr))

      CALL slayer_netcdf_out(SIZE(sl_in%qval_arr),est_gamma_flag,
     $                       sl_in,sl_out,all_deltas_out)

      END SUBROUTINE output_gamma
c-----------------------------------------------------------------------
c     Subprogram 3. allocate_inputs
c     Allocate arrays inside SLAYER inputs type (dictionary-esque)
c-----------------------------------------------------------------------
      SUBROUTINE allocate_inputs(n_k,sl_in)
      INTEGER, INTENT(IN) :: n_k
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
c     Subprogram 4. allocate_outputs
c     Allocate arrays inside SLAYER outputs type (dictionary-esque)
c-----------------------------------------------------------------------
      SUBROUTINE allocate_outputs(n_k,sl_out)
      INTEGER, INTENT(IN) :: n_k
      TYPE(slayer_outputs_type), INTENT(INOUT) :: sl_out

      ALLOCATE(sl_out%dels_db_arr(n_k),sl_out%gamma_sol_arr(n_k),
     $         sl_out%gamma_est_arr(n_k)  )
      RETURN
      END SUBROUTINE allocate_outputs
c-----------------------------------------------------------------------
c     Subprogram 5. shrink_array
c     Remove excess scan array size from memory
c-----------------------------------------------------------------------
      SUBROUTINE shrink_array(arr, new_size)
          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: new_size
          REAL(r8), ALLOCATABLE :: temp(:)

          ALLOCATE(temp(new_size))
          temp(1:new_size) = arr(1:new_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE shrink_array
c-----------------------------------------------------------------------
c     Subprogram 6. grow_array
c     Increase scan array size IF necessary
c-----------------------------------------------------------------------
      SUBROUTINE grow_array(arr, old_size, new_size)
          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: old_size,new_size
          REAL(r8), ALLOCATABLE :: temp(:)

          ALLOCATE(temp(new_size))
          temp(1:old_size) = arr(1:old_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE grow_array
c-----------------------------------------------------------------------
c     Subprogram 7. calc_determinant
c     Calculate determinant of 2x2 and 3x3 matrices
c-----------------------------------------------------------------------
      SUBROUTINE calc_determinant(matk, nk, detk)
      IMPLICIT NONE
            
      ! Arguments
      INTEGER, INTENT(IN) :: nk                           ! Matrix size (2 or 3)
      COMPLEX(r8), DIMENSION(nk,nk), INTENT(IN) :: matk      ! Input matrix
      COMPLEX(r8), INTENT(OUT) :: detk                        ! Determinant result
      INTEGER :: status                     ! Status (0=success, -1=error)
                        
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
c     Subprogram 8. dispersion_det
c     Calculate determinant of coupling matrix e-value problem
c-----------------------------------------------------------------------
      FUNCTION dispersion_det(g_tmp,n_k,sl_in,msing_max)
      
      COMPLEX(r8), INTENT(IN) :: g_tmp
      INTEGER, INTENT(IN) :: n_k,msing_max
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      COMPLEX(r8) :: dispersion_det,det_val,tmp_delta
      COMPLEX(r8), ALLOCATABLE :: delta_Q(:,:),result_matrix(:,:)
      INTEGER :: k

      IF (n_k < 2) THEN
         Q_e = sl_in%Q_e_arr(1)
         Q_i = sl_in%Q_i_arr(1)
         P_perp = sl_in%P_perp_arr(1)
         P_tor = sl_in%P_tor_arr(1)
         tau = sl_in%tau_arr(1)
         D_norm = sl_in%D_norm_arr(1)
         c_beta = sl_in%c_beta_arr(1)
         tauk = sl_in%Qconv_arr(1)
         iota_e = Q_e / (Q_e - Q_i)

         tmp_delta=riccati_f(g_tmp)
         det_val=tmp_delta*(sl_in%lu_arr(1)**(1.0/3.0)) ! DE-NORMALIZE          

         ! Calculate Deltaprime - Delta(Q)
         dispersion_det = sl_in%Re_dp_arr(1) - det_val

      ELSEIF ((n_k == 2) .OR. (n_k == 3)) THEN
         ALLOCATE(delta_Q(n_k,n_k))
         delta_Q=(0.0,0.0)
         DO k=1,msing_max ! maxing out at msing_max
            Q_e = sl_in%Q_e_arr(k)
            Q_i = sl_in%Q_i_arr(k)
            P_perp = sl_in%P_perp_arr(k)
            P_tor = sl_in%P_tor_arr(k)
            tau = sl_in%tau_arr(k)
            D_norm = sl_in%D_norm_arr(k)
            c_beta = sl_in%c_beta_arr(k)
            tauk = sl_in%Qconv_arr(k)
            iota_e = Q_e / (Q_e - Q_i)
     
            delta_Q(k,k)=riccati_f(((g_tmp*sl_in%Qconv_arr(1))
     $           /tauk))
            delta_Q(k,k)=delta_Q(k,k)*sl_in%lu_arr(k)**(1.0/3.0) ! DE-NORMALIZE          
         END DO
   
         ! Calculate Deltaprime - Delta(Q)
         result_matrix = sl_in%dp_matrix - delta_Q

         ! Calculate determinant
         CALL calc_determinant(result_matrix, n_k, det_val)
         dispersion_det = det_val
      ELSE
         WRITE(*,*)"Error: no support for msing > 3"
         stop
      END IF
      END FUNCTION dispersion_det

      SUBROUTINE dispersion_AMR(n_k,sl_in,msing_max,
     $                          scan_width,Q_num,AMR_passes,
     $                          coupling_flag)
      !WRITE(*,*)"------------------------------------------"
      !WRITE(*,'(A,F0.1,A,I2,A)')' >>> Running Adaptive AMR Scan [Width=',scan_width, &
      !                          ', Passes=', AMR_PASSES, ']'
      INTEGER, INTENT(IN) :: n_k,msing_max,Q_num,AMR_passes
      REAL(r8), INTENT(IN) :: scan_width
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      LOGICAL, INTENT(IN) :: coupling_flag

      !COMPLEX(r8), INTENT(OUT), ALLOCATABLE :: Q_store(:)    ! Stores Q coordinates
      !COMPLEX(r8), INTENT(OUT), ALLOCATABLE :: D_store(:)    ! Stores Result Delta
      INTEGER, ALLOCATABLE :: cells(:,:)       ! 4 corners (indices) per cell
      INTEGER, ALLOCATABLE :: new_cells(:,:)   ! Temp array for next level
      INTEGER :: n_cells, n_new_cells, i, j
      INTEGER :: h_idx, pt_idx, c_idx, pass
      INTEGER :: idx_TL, idx_TR, idx_BL, idx_BR ! Corner indices
      INTEGER :: idx_TM, idx_BM, idx_LM, idx_RM, idx_MM ! Midpoint indices
      REAL(r8) :: r_min, r_max, i_min, i_max, ing_step,
     $            ing_coarse,iing_coarse
      LOGICAL :: cross_real, cross_imag
      !INTEGER, INTENT(OUT) :: n_pts

      COMPLEX(r8) :: q_curr

      ! (Re-allocate or just reset counters. Re-allocating ensures clean slate)
      !IF (ALLOCATED(Q_store)) DEALLOCATE(Q_store, D_store)
      !IF (ALLOCATED(hash_head)) DEALLOCATE(hash_head, hash_next)
      !IF (ALLOCATED(cells)) DEALLOCATE(cells, new_cells)

      ! --- 1. Initialize Memory ---
      ALLOCATE(Q_store(MAX_PTS), D_store(MAX_PTS))
      ALLOCATE(hash_head(HASH_SZ), hash_next(MAX_PTS))
      ALLOCATE(cells(4, 200000), new_cells(4, 200000)) ! Estimate cell count
      
      hash_head = 0
      hash_next = 0
      n_pts = 0
      n_cells = 0
  
      ! --- 2. Build Initial Coarse Grid (100x100) ---
      ! We treat the grid as a collection of quadrilateral cells
      ing_step = (2.0*scan_width) / (Q_num - 1)
      
      ! A. Generate Points & Evaluate
      DO i = 1, Q_num
          DO j = 1, Q_num
             ing_coarse = -scan_width + (i - 1) * ing_step
             iing_coarse = -scan_width + (j - 1) * ing_step
             q_curr = CMPLX(ing_coarse, iing_coarse)
             
             ! Check/Compute (Using inline logic to simulate a function call)
             CALL get_or_compute(q_curr, pt_idx,n_k,sl_in,msing_max,
     $                          coupling_flag)
             
             ! If we are not at the right/bottom edge, form a cell with neighbors
             IF (i < Q_num .AND. j < Q_num) THEN
                 n_cells = n_cells + 1
                 ! Store indices of corners: TL, TR, BL, BR (row-major logic)
                 ! Note: This indexing assumes we inserted in order, but for AMR 
                 ! we must rely on the returned pt_idx, not loop counters.
                 ! To simplify, we just store the TL index and calculate others? 
                 ! No, AMR breaks structure. We must look up all 4 corners.
                 
                 ! Top-Left (current)
                 cells(1, n_cells) = pt_idx 
                 
                 ! Top-Right (i+1, j)
                 q_curr = CMPLX(ing_coarse + ing_step, iing_coarse)
                 CALL get_or_compute(q_curr, cells(2, n_cells),n_k,
     $              sl_in,msing_max,coupling_flag)
                 
                 ! Bottom-Left (i, j+1)
                 q_curr = CMPLX(ing_coarse, iing_coarse+ing_step)
                 CALL get_or_compute(q_curr, cells(3, n_cells),n_k,
     $              sl_in,msing_max,coupling_flag)
                 
                 ! Bottom-Right (i+1, j+1)
                 q_curr = CMPLX(ing_coarse + ing_step, iing_coarse+
     $                              ing_step)
                 CALL get_or_compute(q_curr, cells(4, n_cells),n_k,
     $              sl_in,msing_max,coupling_flag)
             END IF
          END DO
      END DO
  
      ! --- 3. Refinement Loops ---
      DO pass = 1, AMR_PASSES
          WRITE(*,'(A,I2,A,I6,A)') '   > Pass ', pass, 
     $         ': Checking ', n_cells, ' cells...'
          n_new_cells = 0
          
          DO c_idx = 1, n_cells
              idx_TL = cells(1, c_idx)
              idx_TR = cells(2, c_idx)
              idx_BL = cells(3, c_idx)
              idx_BR = cells(4, c_idx)
              
              ! Check for contours (Sign changes across the cell)
              ! Real Part Check
              r_min = MIN(REAL(D_store(idx_TL)), 
     $                         REAL(D_store(idx_TR)),
     $                         REAL(D_store(idx_BL)), 
     $                         REAL(D_store(idx_BR)))
              r_max = MAX(REAL(D_store(idx_TL)), 
     $                         REAL(D_store(idx_TR)),
     $                         REAL(D_store(idx_BL)), 
     $                         REAL(D_store(idx_BR)))
              cross_real = (r_min * r_max <= 0.0d0)
              
              ! Imag Part Check
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
                  ! --- REFINE THIS CELL ---
                  ! We need 5 new points: Top-Mid, Bot-Mid, Left-Mid, Right-Mid, Center
                  
                  ! Calculate coords from corners
                  ! TL: Q_store(idx_TL), BR: Q_store(idx_BR)
                  
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
                  
                  ! Create 4 new sub-cells (Top-Left, Top-Right, Bot-Left, Bot-Right)
                  ! Sub 1 (Top-Left)
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
          
          ! Swap arrays for next iteration
          n_cells = n_new_cells
          cells(:, 1:n_cells) = new_cells(:, 1:n_cells)
          
      END DO
      DEALLOCATE(cells, new_cells)
      WRITE(*,*) "AMR Scan Complete. Total Points:", n_pts
      WRITE(*,*)"gslayer.f Q_store(10) = ",Q_store(10)
      RETURN
      END SUBROUTINE dispersion_AMR
          
      SUBROUTINE get_or_compute(q_in,idx_out,n_k,sl_in,msing_max,
     $                          coupling_flag)
      COMPLEX(r8), INTENT(IN) :: q_in
      !INTEGER, INTENT(INOUT) :: n_pts
      !COMPLEX(r8), INTENT(INOUT) :: Q_store(:),D_store(:)
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      INTEGER, INTENT(IN) :: n_k, msing_max
      INTEGER, INTENT(OUT) :: idx_out
      LOGICAL, INTENT(IN) :: coupling_flag
      INTEGER :: h, curr, ix, iy
      COMPLEX(r8) :: delta_val

      ! 1. Calculate Hash
      ix = NINT(REAL(q_in) * HASH_SCALE)
      iy = NINT(AIMAG(q_in) * HASH_SCALE)
      ! Simple hash mix
      h = MOD(ABS(ix * 73856093 + iy * 19349663), HASH_SZ) + 1
      
      ! 2. Check collisions
      curr = hash_head(h)
      DO WHILE (curr /= 0)
          ! Check if point matches (with small tolerance)
          IF (ABS(Q_store(curr) - q_in) < 1.0d-8) THEN
              idx_out = curr
              RETURN ! Found it, return existing index
          END IF
          curr = hash_next(curr)
      END DO
      
      ! 3. Not found: Compute and Store
      n_pts = n_pts + 1

      IF (n_pts > MAX_PTS) THEN
          WRITE(*,*) "ERROR: AMR exceeded MAX_PTS"
          STOP
      END IF
      
      idx_out = n_pts
      Q_store(idx_out) = q_in
      
      ! --- PHYSICS EVALUATION ---
      IF (coupling_flag) THEN
           g_tmp = q_in
           delta_val = dispersion_det(g_tmp, n_k, sl_in, msing_max)
      ELSE
           g_tmp = q_in
           delta_val = riccati_f(g_tmp)
           !delta_val = delta_val - delta_eff
      END IF
      D_store(idx_out) = delta_val
      ! --------------------------
      
      ! 4. Add to Hash Table
      hash_next(idx_out) = hash_head(h)
      hash_head(h) = idx_out
      
      END SUBROUTINE get_or_compute
      END MODULE gslayer_mod