      MODULE gslayer_mod

      USE omp_lib

      USE sglobal_mod, ONLY: out_unit,r8, mu0, m_p, chag, lnLamb,
     $   Q_e,Q_i,pr,pe,c_beta,ds,tau,
     $   eta,visc,rho_s,lu,omega_e,omega_i,
     $   delta_n,
     $   Q
      USE delta_mod, ONLY: riccati,riccati_out,
     $   parflow_flag,PeOhmOnly_flag

      USE params_mod

      USE layerinputs_mod

      USE slayer_netcdf_mod

      USE grid, ONLY : powspace,linspace

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
     $            lbeta,tau_i,tau_h,tau_r,tau_v
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

      ! note Q depends on Qconv even if omega is fixed.
      Q=Qconv*omega
      Q_e=-Qconv*omega_e
      Q_i=-Qconv*omega_i

      ! This is the most critical parameter
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
     $        TRIM(sm)//"_n"//TRIM(sn)//".out",
     $        STATUS="UNKNOWN")
         WRITE(out_unit,'(1x,5(a17))'),"inQ","RE(delta)",
     $        "IM(delta)","jxb","bal"
         DO i=0,inum
            WRITE(out_unit,'(1x,5(es17.8e3))')
     $           inQs(i),REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
         ENDDO
         CLOSE(out_unit)
      ENDIF

      ! Identify the threshold from the maximum of the balance parameter
      index=MAXLOC(bal)
      Q_sol=inQs(index(1))
      omega_sol=inQs(index(1))/Qconv
      br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))
      DEALLOCATE(inQs,deltal,jxbl,bal)

      RETURN
      END SUBROUTINE gpec_slayer
c-----------------------------------------------------------------------
c     Subprogram 2. growthrate_scan
c     Run stability scan on real and imaginary rotation axes
c-----------------------------------------------------------------------
      SUBROUTINE old_growthrate_scan(qval,inQ,inQ_e,inQ_i,inc_beta,
     $           inds,intau,inQ0,inpr,inpe,scan_radius,reQ_num,
     $           Re_deltas,Im_deltas,inQs,iinQs)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe
      INTEGER, INTENT(IN) :: qval,scan_radius,reQ_num
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: inQs,iinQs
      INTEGER :: n_fine
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x_s,y_s
      REAL(r8) :: x, y, dx, dy, ric_val,start_time,end_time
      REAL(r8), DIMENSION(2) :: start, end, step
      REAL(r8) :: coarse_threshold,fine_threshold,grad_tol,wr
      COMPLEX(r8) :: delta
      INTEGER :: k, w, n_steps, i2, i,j
      logical :: found_zero

      INTEGER :: n_points, idx_1, idx_2
      REAL(r8) :: x_1, y_1, x_2, y_2, dist, min_dist1, min_dist2

      REAL(r8), DIMENSION(:,:), ALLOCATABLE,
     $             INTENT(OUT) :: Re_deltas,Im_deltas
      wr=0.0
      ! Parameters
      grad_tol = 1.0!abs(0.1*deltaprime)     ! Tolerance for steep gradient
      n_steps = reQ_num       ! Initial number of steps
!     n_fine = 10
!     if (abs(deltaprime)>100) then
!       coarse_threshold = abs(0.95*deltaprime)    ! Tolerance for zero value
!       fine_threshold = abs(0.1*deltaprime)
!       n_fine = n_fine * 2
!     else if (abs(deltaprime)>50) then
!       coarse_threshold = abs(0.8*deltaprime)    ! Tolerance for zero value
!       fine_threshold = abs(0.1*deltaprime)
!       n_fine = n_fine * 2
!     else if (abs(deltaprime)>20) then
!       coarse_threshold = abs(0.66*deltaprime)    ! Tolerance for zero value
!       fine_threshold = abs(0.05*deltaprime)
!       n_fine = n_fine * 2
!     else if (abs(deltaprime)>10) then
!       coarse_threshold = abs(0.3*deltaprime)    ! Tolerance for zero value
!       fine_threshold = abs(0.03*deltaprime)
!     else
!       coarse_threshold = abs(0.1*deltaprime)    ! Tolerance for zero value
!       fine_threshold = abs(0.01*deltaprime)
!     end if
      ALLOCATE(inQs(n_steps))
      ALLOCATE(iinQs(n_steps))
      ALLOCATE(Re_deltas(n_steps,n_steps))
      ALLOCATE(Im_deltas(n_steps,n_steps))
      inQs = linspace((-scan_radius+wr),(scan_radius+wr),n_steps)
      iinQs = linspace((-scan_radius+wr),(scan_radius+wr),n_steps)
      i=1
      j=1
      !CALL OMP_SET_NUM_THREADS(4)
      !PRINT *, "Max threads: ",OMP_GET_MAX_THREADS()
      !CALL cpu_time(start_time)
      !!$OMP PARALLEL DO PRIVATE(j,y,delta) !BIG REWRITE NEEDED
      DO i = 1, n_steps
          DO j = 1, n_steps
              !PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
              y=iinQs(j)
              delta = riccati(inQs(i),inQ_e,inQ_i,inpr,inc_beta,
     $                        inds,intau,inpe,iinQ=y)
              Re_deltas(i,j) = REAL(delta)
              Im_deltas(i,j) = AIMAG(delta)
              !PRINT *,'inQs(i):',inQs(i)
              !PRINT *,'y:',y
              !PRINT *,'delta:',delta
          END DO
          !stop
      END DO
      !!$OMP END PARALLEL DO
      !CALL cpu_time(end_time)
      !PRINT *,'Layer scan time:',end_time-start_time,'seconds'
      !test comment
      RETURN
      END SUBROUTINE old_growthrate_scan
c-----------------------------------------------------------------------
c     Subprogram 2. growthrate_scan
c     Run stability scan on real and imaginary rotation axes
c-----------------------------------------------------------------------
      SUBROUTINE growthrate_scan(qval,inQ,inQ_e,inQ_i,inc_beta,
     $         inds,intau,inQ0,inpr,inpe,scan_radius,coarse_grid_size,
     $         deltaprime,results)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe
      INTEGER, INTENT(IN) :: qval,scan_radius,coarse_grid_size
      INTEGER :: new_scan_radius,new_coarse_grid_size
      !REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: inQs,iinQs
      COMPLEX(r8) :: delta
      REAL(r8), INTENT(IN) :: deltaprime
      INTEGER :: fine_grid_size, new_fine_grid_size
      REAL(r8), PARAMETER :: tolerance = 1.0E-6
      REAL(r8) :: delta_real, delta_imag, threshold
      INTEGER :: i, j, k, l, m, count, match_count
      LOGICAL :: repeat
      REAL(r8) :: inQ_step, iinQ_step, inQ_fine, iinQ_fine,
     $            inQ_coarse, iinQ_coarse, inQ_coarse_shifted
      REAL(r8), DIMENSION(2) :: both_coarse_inQs
      TYPE(result_type), INTENT(INOUT) :: results
      INTEGER :: max_points, new_max_points

      !!!!!!!!!!!!!!!!
      repeat = .FALSE.

      fine_grid_size = 6
      max_points = coarse_grid_size**2 * (1 + (fine_grid_size-1)**2)

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
      inQ_step = (2.0 * scan_radius) / (coarse_grid_size - 1)
      iinQ_step = (2.0 * scan_radius) / (coarse_grid_size - 1)

      match_count = 0
      ! Coarse grid loop
      DO i = 1, coarse_grid_size
          DO j = 1, coarse_grid_size
              inQ_coarse = -scan_radius + (i - 1) * inQ_step
              inQ_coarse_shifted = -scan_radius + (i - 2) * inQ_step
              both_coarse_inQs(1) = inQ_coarse
              both_coarse_inQs(2) = inQ_coarse_shifted
              iinQ_coarse = -scan_radius + (j - 1) * iinQ_step

              ! Evaluate riccati function
              delta = riccati(inQ_coarse,inQ_e,inQ_i,inpr,inc_beta,
     $                        inds,intau,inpe,iinQ=iinQ_coarse)
              delta_real = REAL(delta)
              delta_imag = AIMAG(delta)

              ! Store coarse grid point
              count = count + 1
              results%inQs(count) = inQ_coarse
              results%iinQs(count) = iinQ_coarse
              results%Re_deltas(count) = delta_real
              results%Im_deltas(count) = delta_imag

              IF (ABS(deltaprime) > 8) THEN
                threshold = ABS(deltaprime)**(1./3.)
              ELSE
                threshold = 0.25 * ABS(deltaprime)
              END IF

              ! Check if refinement is needed
              IF ((ABS(delta_real) > threshold)) THEN
      !        IF ((ABS(delta_real) > threshold) .AND. (SIGN(1.0,
      !$                   delta_real) == SIGN(1.0, deltaprime))) THEN

                  ! Fine grid loop
                  DO m = 1, 2
                    DO k = 2, fine_grid_size
                      DO l = 2, fine_grid_size
                        inQ_fine = both_coarse_inQs(m)+(k-1)*inQ_step/
     $                   (fine_grid_size - 1)
                        iinQ_fine = iinQ_coarse + (l-1) * iinQ_step /
     $                   (fine_grid_size - 1)

                       IF ((ABS(both_coarse_inQs(m) - inQ_fine) <
     $                   tolerance) .AND. (ABS(iinQ_coarse -
     $                   iinQ_fine) < tolerance)) CYCLE

                        ! Evaluate riccati function
                        delta = riccati(inQ_fine,inQ_e,inQ_i,inpr,
     $                     inc_beta,inds,intau,inpe,iinQ=iinQ_fine)
                        delta_real = REAL(delta)
                        delta_imag = AIMAG(delta)

                        IF (ABS(delta_real) > ABS(deltaprime)) THEN
                          match_count = match_count + 1
                        END IF

                        ! Store fine grid point
                        count = count + 1
                        results%inQs(count) = inQ_fine
                        results%iinQs(count) = iinQ_fine
                        results%Re_deltas(count) = delta_real
                        results%Im_deltas(count) = delta_imag
                      END DO
                    END DO
                  END DO
              END IF
          END DO
      END DO

      !!!!!

      IF (match_count == 0) THEN
        WRITE(*,*)"No match found, increasing scan radius"
        repeat = .TRUE.
        new_scan_radius = scan_radius + 2
        new_coarse_grid_size = coarse_grid_size + 100
        new_fine_grid_size = 8
      ELSE IF (match_count > 0 .AND. match_count < 3) THEN
        WRITE(*,*)"Match not definitive, increasing scan resolution"
        repeat = .TRUE.
        new_scan_radius = scan_radius
        new_coarse_grid_size = coarse_grid_size
        new_fine_grid_size = 10
      ELSE
        repeat = .FALSE.
        WRITE(*,*)"Match found"

      END IF

      IF (repeat) THEN
        WRITE(*,*)"Rerunning growth rate scan"

        !DEALLOCATE(results%inQs(max_points), results%iinQs(max_points))
        !DEALLOCATE(results%Re_deltas(max_points),
      !$  results%Im_deltas(max_points))

        new_max_points = new_coarse_grid_size**2 * (1 +
     $       (new_fine_grid_size-1)**2)

        ! Resize arrays to new max number of points
        CALL grow_array(results%inQs, max_points, new_max_points)
        CALL grow_array(results%iinQs, max_points, new_max_points)
        CALL grow_array(results%Re_deltas, max_points, new_max_points)
        CALL grow_array(results%Im_deltas, max_points, new_max_points)

        results%inQs=0.0
        results%iinQs=0.0
        results%Re_deltas=0.0
        results%Im_deltas=0.0
        ! Initialize counter
        count = 0

        ! Calculate step sizes
        inQ_step = (2.0 * new_scan_radius)/(new_coarse_grid_size - 1)
        iinQ_step = (2.0 * new_scan_radius)/(new_coarse_grid_size - 1)

        match_count = 0
        ! Coarse grid loop
        DO i = 1, new_coarse_grid_size
            DO j = 1, new_coarse_grid_size
                inQ_coarse = -new_scan_radius + (i - 1) * inQ_step
                inQ_coarse_shifted = -scan_radius + (i - 2) * inQ_step
                both_coarse_inQs(1) = inQ_coarse
                both_coarse_inQs(2) = inQ_coarse_shifted
                iinQ_coarse = -new_scan_radius + (j - 1) * iinQ_step

              ! Evaluate riccati function
                delta = riccati(inQ_coarse,inQ_e,inQ_i,inpr,inc_beta,
     $                          inds,intau,inpe,iinQ=iinQ_coarse)
                delta_real = REAL(delta)
                delta_imag = AIMAG(delta)

                ! Store coarse grid point
                count = count + 1
                results%inQs(count) = inQ_coarse
                results%iinQs(count) = iinQ_coarse
                results%Re_deltas(count) = delta_real
                results%Im_deltas(count) = delta_imag

                IF (ABS(deltaprime) > 8) THEN
                  threshold = ABS(deltaprime)**(1./3.)
                ELSE
                  threshold = 0.25 * ABS(deltaprime)
                END IF

                ! Check if refinement is needed
                IF ((ABS(delta_real) > threshold)) THEN
                    ! Fine grid loop
                    DO m = 1,2
                      DO k = 2, new_fine_grid_size
                        DO l = 2, new_fine_grid_size
                          inQ_fine=both_coarse_inQs(m)+(k-1)*inQ_step/
     $                     (new_fine_grid_size - 1)
                          iinQ_fine = iinQ_coarse + (l-1) * iinQ_step/
     $                     (new_fine_grid_size - 1)

                          IF ((ABS(both_coarse_inQs(m) - inQ_fine) <
     $                      tolerance) .AND. (ABS(iinQ_coarse -
     $                      iinQ_fine) < tolerance)) CYCLE

                            ! Evaluate riccati function
                          delta = riccati(inQ_fine,inQ_e,inQ_i,inpr,
     $                       inc_beta,inds,intau,inpe,iinQ=iinQ_fine)
                          delta_real = REAL(delta)
                          delta_imag = AIMAG(delta)

                          IF (ABS(delta_real)>ABS(deltaprime)) THEN
                              match_count = match_count + 1
                          END IF

                            ! Store fine grid point
                          count = count + 1
                          results%inQs(count) = inQ_fine
                          results%iinQs(count) = iinQ_fine
                          results%Re_deltas(count) = delta_real
                          results%Im_deltas(count) = delta_imag
                        END DO
                      END DO
                    END DO
                END IF
            END DO
        END DO

      END IF

      ! Set the actual count of points
      results%count = count

      ! Resize arrays to actual number of points
      CALL shrink_array(results%inQs, count)
      CALL shrink_array(results%iinQs, count)
      CALL shrink_array(results%Re_deltas, count)
      CALL shrink_array(results%Im_deltas, count)

      RETURN
      END SUBROUTINE growthrate_scan

      SUBROUTINE shrink_array(arr, new_size)
          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: new_size
          REAL(r8), ALLOCATABLE :: temp(:)

          ALLOCATE(temp(new_size))
          temp(1:new_size) = arr(1:new_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE shrink_array

      SUBROUTINE grow_array(arr, old_size, new_size)
          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: old_size,new_size
          REAL(r8), ALLOCATABLE :: temp(:)

          ALLOCATE(temp(new_size))
          temp(1:old_size) = arr(1:old_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE grow_array
      END MODULE gslayer_mod
