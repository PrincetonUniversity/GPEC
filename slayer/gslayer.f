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
c     Subprogram 3. scan_grid
c     Run stability scan on real and imaginary rotation axes
c-----------------------------------------------------------------------
      SUBROUTINE output_gamma(est_gamma_flag,sl_in,sl_out)

      ! Declarations (include necessary type declarations from original code)
      LOGICAL, INTENT(IN) :: est_gamma_flag
      TYPE(slayer_inputs_type), INTENT(IN) :: sl_in
      TYPE(slayer_outputs_type), INTENT(IN) :: sl_out

      CALL slayer_netcdf_out(SIZE(sl_in%qval_arr),est_gamma_flag,
     $                       sl_in,sl_out)

      END SUBROUTINE output_gamma

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

      SUBROUTINE allocate_outputs(n_k,sl_out)
      INTEGER, INTENT(IN) :: n_k
      TYPE(slayer_outputs_type), INTENT(INOUT) :: sl_out

      ALLOCATE(sl_out%dels_db_arr(n_k),sl_out%gamma_sol_arr(n_k),
     $         sl_out%gamma_est_arr(n_k)  )
      RETURN
      END SUBROUTINE allocate_outputs
c-----------------------------------------------------------------------
c     Subprogram 2. growthrate_scan
c     Set up and iterate stability scans IF no match is found
c-----------------------------------------------------------------------
      SUBROUTINE growthrate_scan(qval,my_lu,inQ,inQ_e,inQ_i,inc_beta,
     $         inds,intau,inQ0,inpr,inpe,scan_radius,ncoarse,
     $         compress_deltas,deltaprime,results)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe,my_lu
      INTEGER, INTENT(IN) :: qval,scan_radius,ncoarse
      REAL(r8), INTENT(IN) :: deltaprime
      LOGICAL, INTENT(IN) :: compress_deltas
      TYPE(result_type), INTENT(INOUT) :: results

      COMPLEX(r8) :: delta
      INTEGER :: new_scan_radius,new_ncoarse
      INTEGER :: nfine, new_nfine
      REAL(r8), PARAMETER :: tolerance = 1.0E-6
      REAL(r8) :: delta_real, delta_imag, threshold
      INTEGER :: i, j, k, l, m, count, match_count
      LOGICAL :: repeat
      REAL(r8) :: inQ_step, iinQ_step, inQ_fine, iinQ_fine,
     $            inQ_coarse, iinQ_coarse
      INTEGER :: max_points, new_max_points
      INTEGER :: ci, cj, nx, ny
      REAL(r8) :: dx, dy, overlap_factor
      INTEGER :: fi, fj
      REAL(r8) :: fine_dx, fine_dy, overlap_x, overlap_y
      REAL(r8) :: x_start, x_end, y_start, y_end, x, y
      !!!!!!!!!!!!!!!!
      repeat = .FALSE.
      dx = 1.0
      dy = 1.0
      nfine = 6
      overlap_factor = 0.5
      max_points = ncoarse**2 * ((nfine)**2 - 1)

      ! Allocate arrays with maximum possible size
      ALLOCATE(results%inQs(max_points), results%iinQs(max_points))
      ALLOCATE(results%Re_deltas(max_points),
     $ results%Im_deltas(max_points))

      results%inQs=0.0
      results%iinQs=0.0
      results%Re_deltas=0.0
      results%Im_deltas=0.0
      ! Initialize counter
      count = 0

      ! Calculate step sizes
      inQ_step = (2.0 * scan_radius) / (ncoarse - 1)
      iinQ_step = (2.0 * scan_radius) / (ncoarse - 1)
      dx = inQ_step
      dy = iinQ_step

      match_count = 0
      ! Run scan
      CALL scan_grid(inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,my_lu,
     $    scan_radius,ncoarse,nfine,deltaprime,compress_deltas,
     $    results,count,match_count,dx,dy)

      ! Set the actual count of points
      results%count = count

      IF (count < max_points) THEN
        ! Resize arrays to actual number of points
        CALL shrink_array(results%inQs, count)
        CALL shrink_array(results%iinQs, count)
        CALL shrink_array(results%Re_deltas, count)
        CALL shrink_array(results%Im_deltas, count)
      END IF

      RETURN
      END SUBROUTINE growthrate_scan
c-----------------------------------------------------------------------
c     Subprogram 3. scan_grid
c     Run stability scan on real and imaginary rotation axes
c-----------------------------------------------------------------------
      SUBROUTINE scan_grid(inQ_e,inQ_i,inpr,inc_beta,inds,intau, 
     $     inpe,my_lu,scan_radius,ncoarse,nfine,deltaprime,
     $     compress_deltas,results,count,match_count,dx,dy)
      
      ! Declarations (include necessary type declarations from original code)
      REAL(r8), INTENT(IN) :: inQ_e,inQ_i,inpr,inc_beta,inds,
     $     intau,inpe,my_lu,deltaprime
      INTEGER, INTENT(IN) :: scan_radius, ncoarse, nfine
      LOGICAL, INTENT(IN) :: compress_deltas
      TYPE(result_type), INTENT(INOUT) :: results
      INTEGER, INTENT(INOUT) :: count, match_count
      REAL(r8), INTENT(INOUT) :: dx, dy
      
      ! Local variables
      REAL(r8) :: inQ_step, iinQ_step, inQ_fine, iinQ_fine,
     $     inQ_coarse, iinQ_coarse
      REAL(r8) :: delta_real, delta_imag, threshold
      COMPLEX(r8) :: delta
      REAL(r8) :: fine_dx, fine_dy, overlap_x, overlap_y
      REAL(r8) :: x_start, x_end, y_start, y_end
      INTEGER :: i, j, fi, fj
      REAL(r8), PARAMETER :: tolerance = 1.0E-6
      REAL(r8) :: overlap_factor = 0.5

      ! Calculate step sizes
      inQ_step = (2.0 * scan_radius) / (ncoarse - 1)
      iinQ_step = (2.0 * scan_radius) / (ncoarse - 1)
      dx = inQ_step
      dy = iinQ_step
      count = 0
      
      DO i = 1, ncoarse
        DO j = 1, ncoarse
          inQ_coarse = -scan_radius + (i - 1) * inQ_step
          iinQ_coarse = -scan_radius + (j - 1) * iinQ_step
          ! Evaluate riccati FUNCTION
          delta = riccati(inQ_coarse,inQ_e,inQ_i,inpr,inc_beta,
     $                        inds,intau,inpe,iinQ=iinQ_coarse)
          delta_real = REAL(delta)*(my_lu**(1.0/3.0)) ! Critical normalization
          delta_imag = AIMAG(delta)*(my_lu**(1.0/3.0)) ! Critical normalization

          count = count + 1
          results%inQs(count) = inQ_coarse
          results%iinQs(count) = iinQ_coarse
          results%Re_deltas(count) = delta_real
          results%Im_deltas(count) = delta_imag

        END DO
      END DO
      END SUBROUTINE scan_grid
c-----------------------------------------------------------------------
c     Subprogram 4. shrink_array
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
c     Subprogram 5. grow_array
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
c
c
c
c-----------------------------------------------------------------------
c     Subprogram 6. determinant
c     Increase scan array size IF necessary
c-----------------------------------------------------------------------
      SUBROUTINE calc_determinant(matk, nk, detk)
      IMPLICIT NONE
            
      ! Arguments
      INTEGER, INTENT(IN) :: nk                           ! Matrix size (2 or 3)
      COMPLEX(r8), DIMENSION(nk,nk), INTENT(IN) :: matk      ! Input matrix
      COMPLEX(r8), INTENT(OUT) :: detk                        ! Determinant result
      INTEGER :: status                     ! Status (0=success, -1=error)
            
      ! Local variables
            
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
      
      !===========================================================================
      ! Initialize the adaptive grid
      !===========================================================================
      subroutine init_grid(grid, omega_min, omega_max, gamma_min, gamma_max, initial_capacity)
          type(adaptive_grid), intent(out) :: grid
          real(dp), intent(in) :: omega_min, omega_max, gamma_min, gamma_max
          integer, intent(in) :: initial_capacity
          
          grid%omega_min = omega_min
          grid%omega_max = omega_max
          grid%gamma_min = gamma_min
          grid%gamma_max = gamma_max
          grid%capacity = initial_capacity
          grid%npoints = 0
          
          allocate(grid%points(initial_capacity))
          
      end subroutine init_grid
      
      !===========================================================================
      ! Add a point to the grid (with automatic reallocation if needed)
      !===========================================================================
      subroutine add_point(grid, omega, gamma, delta_val)
          type(adaptive_grid), intent(inout) :: grid
          real(dp), intent(in) :: omega, gamma
          complex(dp), intent(in) :: delta_val
          
          type(grid_point), allocatable :: temp(:)
          integer :: new_capacity
          
          ! Check if we need to reallocate
          if (grid%npoints >= grid%capacity) then
              new_capacity = grid%capacity * 2
              allocate(temp(new_capacity))
              temp(1:grid%npoints) = grid%points(1:grid%npoints)
              call move_alloc(temp, grid%points)
              grid%capacity = new_capacity
          end if
          
          ! Add the new point
          grid%npoints = grid%npoints + 1
          grid%points(grid%npoints)%omega = omega
          grid%points(grid%npoints)%gamma = gamma
          grid%points(grid%npoints)%delta = delta_val
          grid%points(grid%npoints)%computed = .true.
          grid%points(grid%npoints)%distance_to_contour = huge(1.0_dp)
          
      end subroutine add_point
      
      !===========================================================================
      ! Check if a point already exists in the grid (within tolerance)
      !===========================================================================
      logical function point_exists(grid, omega, gamma, tol)
          type(adaptive_grid), intent(in) :: grid
          real(dp), intent(in) :: omega, gamma, tol
          integer :: i
          
          point_exists = .false.
          do i = 1, grid%npoints
              if (abs(grid%points(i)%omega - omega) < tol .and. &
                  abs(grid%points(i)%gamma - gamma) < tol) then
                  point_exists = .true.
                  return
              end if
          end do
          
      end function point_exists
      
      !===========================================================================
      ! Perform initial coarse scan
      !===========================================================================
      subroutine coarse_scan(grid, delta_function, n_omega, n_gamma)
          type(adaptive_grid), intent(inout) :: grid
          interface
              function delta_function(omega, gamma) result(delta)
                  import :: dp
                  real(dp), intent(in) :: omega, gamma
                  complex(dp) :: delta
              end function delta_function
          end interface
          integer, intent(in) :: n_omega, n_gamma
          
          real(dp) :: omega, gamma, domega, dgamma
          complex(dp) :: delta_val
          integer :: i, j
          
          domega = (grid%omega_max - grid%omega_min) / real(n_omega - 1, dp)
          dgamma = (grid%gamma_max - grid%gamma_min) / real(n_gamma - 1, dp)
          
          do i = 1, n_omega
              omega = grid%omega_min + real(i-1, dp) * domega
              do j = 1, n_gamma
                  gamma = grid%gamma_min + real(j-1, dp) * dgamma
                  
                  ! Compute delta at this point
                  delta_val = delta_function(omega, gamma)
                  call add_point(grid, omega, gamma, delta_val)
                  
              end do
          end do
          
          print *, 'Coarse scan complete. Points computed:', grid%npoints
          
      end subroutine coarse_scan
      
      !===========================================================================
      ! Identify contour regions based on sign changes
      !===========================================================================
      subroutine identify_contour_regions(grid, deltaprime, contour_tol)
          type(adaptive_grid), intent(inout) :: grid
          complex(dp), intent(in) :: deltaprime
          real(dp), intent(in) :: contour_tol
          
          integer :: i
          real(dp) :: dist_real, dist_imag, min_dist
          
          ! For each point, compute distance to contours
          do i = 1, grid%npoints
              dist_real = abs(real(grid%points(i)%delta) - real(deltaprime))
              dist_imag = abs(aimag(grid%points(i)%delta) - aimag(deltaprime))
              min_dist = min(dist_real, dist_imag)
              grid%points(i)%distance_to_contour = min_dist
          end do
          
      end subroutine identify_contour_regions
      
      !===========================================================================
      ! Adaptive refinement around contours
      !===========================================================================
      subroutine adaptive_refine(grid, delta_function, deltaprime, &
                                refinement_width, refinement_levels, min_spacing)
          type(adaptive_grid), intent(inout) :: grid
          interface
              function delta_function(omega, gamma) result(delta)
                  import :: dp
                  real(dp), intent(in) :: omega, gamma
                  complex(dp) :: delta
              end function delta_function
          end interface
          complex(dp), intent(in) :: deltaprime
          real(dp), intent(in) :: refinement_width
          integer, intent(in) :: refinement_levels
          real(dp), intent(in) :: min_spacing
          
          integer :: level, i, j, n_original
          real(dp) :: omega, gamma, spacing
          complex(dp) :: delta_val
          logical, allocatable :: needs_refinement(:)
          real(dp) :: omega_new, gamma_new
          integer :: n_refined
          
          do level = 1, refinement_levels
              n_original = grid%npoints
              allocate(needs_refinement(n_original))
              
              ! Identify points that need refinement
              do i = 1, n_original
                  needs_refinement(i) = grid%points(i)%distance_to_contour < refinement_width
              end do
              
              n_refined = 0
              spacing = refinement_width / (2.0_dp**level)
              
              ! Add refined points around identified regions
              do i = 1, n_original
                  if (needs_refinement(i)) then
                      omega = grid%points(i)%omega
                      gamma = grid%points(i)%gamma
                      
                      ! Add points in a 3x3 grid around this point
                      do j = -1, 1
                          omega_new = omega + real(j, dp) * spacing
                          if (omega_new < grid%omega_min .or. omega_new > grid%omega_max) cycle
                          
                          gamma_new = gamma - spacing
                          if (gamma_new >= grid%gamma_min .and. gamma_new <= grid%gamma_max) then
                              if (.not. point_exists(grid, omega_new, gamma_new, min_spacing)) then
                                  delta_val = delta_function(omega_new, gamma_new)
                                  call add_point(grid, omega_new, gamma_new, delta_val)
                                  n_refined = n_refined + 1
                              end if
                          end if
                          
                          if (j /= 0) then
                              gamma_new = gamma
                              if (.not. point_exists(grid, omega_new, gamma_new, min_spacing)) then
                                  delta_val = delta_function(omega_new, gamma_new)
                                  call add_point(grid, omega_new, gamma_new, delta_val)
                                  n_refined = n_refined + 1
                              end if
                          end if
                          
                          gamma_new = gamma + spacing
                          if (gamma_new >= grid%gamma_min .and. gamma_new <= grid%gamma_max) then
                              if (.not. point_exists(grid, omega_new, gamma_new, min_spacing)) then
                                  delta_val = delta_function(omega_new, gamma_new)
                                  call add_point(grid, omega_new, gamma_new, delta_val)
                                  n_refined = n_refined + 1
                              end if
                          end if
                      end do
                  end if
              end do
              
              print *, 'Refinement level', level, ': Added', n_refined, 'points'
              
              ! Update distances for new points
              call identify_contour_regions(grid, deltaprime, refinement_width)
              
              deallocate(needs_refinement)
              
              ! Stop if no new points were added
              if (n_refined == 0) exit
          end do
          
      end subroutine adaptive_refine
      
      !===========================================================================
      ! Marching squares helper for more accurate contour following
      !===========================================================================
      subroutine refine_with_marching_squares(grid, delta_function, deltaprime, &
                                             contour_tol, max_new_points)
          type(adaptive_grid), intent(inout) :: grid
          interface
              function delta_function(omega, gamma) result(delta)
                  import :: dp
                  real(dp), intent(in) :: omega, gamma
                  complex(dp) :: delta
              end function delta_function
          end interface
          complex(dp), intent(in) :: deltaprime
          real(dp), intent(in) :: contour_tol
          integer, intent(in) :: max_new_points
          
          ! Implementation of marching squares refinement
          ! This would trace along detected contours for extra precision
          ! Left as a stub for brevity, but can be expanded if needed
          
      end subroutine refine_with_marching_squares
      
      !===========================================================================
      ! Export grid to file
      !===========================================================================
      subroutine export_grid(grid, filename)
          type(adaptive_grid), intent(in) :: grid
          character(len=*), intent(in) :: filename
          
          integer :: i, unit_num
          
          open(newunit=unit_num, file=filename, status='replace', action='write')
          
          ! Write header
          write(unit_num, '(A)') '# omega, gamma, Re(delta), Im(delta), distance_to_contour'
          write(unit_num, '(A,I0)') '# Number of points: ', grid%npoints
          
          ! Write data
          do i = 1, grid%npoints
              write(unit_num, '(5E16.8)') grid%points(i)%omega, grid%points(i)%gamma, &
                                          real(grid%points(i)%delta), &
                                          aimag(grid%points(i)%delta), &
                                          grid%points(i)%distance_to_contour
          end do
          
          close(unit_num)
          
          print *, 'Grid exported to ', trim(filename)
          print *, 'Total points: ', grid%npoints
          
      end subroutine export_grid
      
      !===========================================================================
      ! Main driver routine
      !===========================================================================
      subroutine find_contours(delta_function, omega_min, omega_max, &
                              gamma_min, gamma_max, deltaprime, &
                              n_coarse_omega, n_coarse_gamma, &
                              refinement_width, refinement_levels, &
                              min_spacing, output_file)
          interface
              function delta_function(omega, gamma) result(delta)
                  import :: dp
                  real(dp), intent(in) :: omega, gamma
                  complex(dp) :: delta
              end function delta_function
          end interface
          real(dp), intent(in) :: omega_min, omega_max, gamma_min, gamma_max
          complex(dp), intent(in) :: deltaprime
          integer, intent(in) :: n_coarse_omega, n_coarse_gamma
          real(dp), intent(in) :: refinement_width
          integer, intent(in) :: refinement_levels
          real(dp), intent(in) :: min_spacing
          character(len=*), intent(in) :: output_file
          
          type(adaptive_grid) :: grid
          integer :: initial_capacity
          
          print *, '========================================'
          print *, 'Starting adaptive contour finding'
          print *, '========================================'
          
          ! Initialize grid
          initial_capacity = n_coarse_omega * n_coarse_gamma * 4
          call init_grid(grid, omega_min, omega_max, gamma_min, gamma_max, initial_capacity)
          
          ! Perform coarse scan
          print *, 'Step 1: Coarse scanning...'
          call coarse_scan(grid, delta_function, n_coarse_omega, n_coarse_gamma)
          
          ! Identify contour regions
          print *, 'Step 2: Identifying contour regions...'
          call identify_contour_regions(grid, deltaprime, refinement_width)
          
          ! Adaptive refinement
          print *, 'Step 3: Adaptive refinement...'
          call adaptive_refine(grid, delta_function, deltaprime, &
                             refinement_width, refinement_levels, min_spacing)
          
          ! Export results
          print *, 'Step 4: Exporting results...'
          call export_grid(grid, output_file)
          
          ! Clean up
          deallocate(grid%points)
          
      end subroutine find_contours
  
      END MODULE gslayer_mod
