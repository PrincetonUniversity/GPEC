      MODULE gslayer_mod
      
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

      FUNCTION log_prob(params,inQ,inQ_e,inQ_i,inc_beta,inds,
     $       intau,inQ0,inpr,inpe,deltaprime,sigma) RESULT(lp)
          real(r8), dimension(2), intent(in) :: params
          real(r8), intent(in) :: deltaprime, sigma
          real(r8) :: lp, distance
          REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inc_beta,inds,
     $                      intau,inQ0,inpr,inpe
          ! Calculate distance from the target value
          distance = abs((riccati(params(1),
     $                               inQ_e,
     $                               inQ_i,
     $                               inpr,
     $                               inc_beta,
     $                               inds,
     $                               intau,
     $                               inpe,
     $                               iinQ=params(2)) -
     $              deltaprime))

          ! Calculate the log probability (negative log-likelihood of a Gaussian)
          lp = -0.5 * (distance / sigma)**2
          RETURN
      END FUNCTION log_prob
c-----------------------------------------------------------------------
c     Subprogram 2. gamma_stability_scan
c     Run grid packed slayer stab. scan around omega_ExB and gamma axes
c-----------------------------------------------------------------------
      SUBROUTINE growthrate_search(qval,inQ,inQ_e,inQ_i,inc_beta,
     $           inds,intau,inQ0,inpr,inpe,deltaprime,scan_radius,
     $           deltas,inQs,iinQs,roots,growthrate,growthrate_loc)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe,deltaprime
      INTEGER, INTENT(IN) :: qval,scan_radius
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      real(r8), intent(out), dimension(:,:), ALLOCATABLE :: roots ! Array to store (x, y) pairs of zeros
      integer :: num_zeros             ! Number of zeros found
      real(r8), dimension(:), ALLOCATABLE, INTENT(OUT) :: inQs,iinQs
      integer :: n_fine
      real(r8), dimension(:), ALLOCATABLE :: x_s,y_s
      real(r8) :: x, y, dx, dy, ric_val, grad_x, grad_y
      real(r8), dimension(2) :: start, end, step
      real(r8) :: zero_tol, grad_tol, wr, delta
      integer :: i, j, k, w, n_steps
      logical :: found_zero
      real(r8), intent(out) :: growthrate_loc, growthrate    ! Output (x, y) pair

      integer :: n_points, idx_1, idx_2
      real(r8) :: x_1, y_1, x_2, y_2, dist, min_dist1, min_dist2

      REAL(r8), DIMENSION(:,:), ALLOCATABLE,
     $             INTENT(OUT) :: deltas
      wr = inQ+inQ_e
      ! Parameters
      zero_tol = abs(0.1*deltaprime)    ! Tolerance for zero value
      grad_tol = 1.0!abs(0.1*deltaprime)     ! Tolerance for steep gradient
      n_steps = 100       ! Initial number of steps
      n_fine = 10
      if (abs(deltaprime)>50) then
        n_fine = n_fine * 2
      end if
      ALLOCATE(x_s(n_fine),y_s(n_fine))
      ALLOCATE(inQs(n_steps))
      ALLOCATE(roots(10,2))
      ALLOCATE(deltas(n_steps,n_steps))

      ! Initial coarse grid search

      dx = 2 * dble(scan_radius) / dble(n_steps)
      dy = 2 * dble(scan_radius) / dble(n_steps)
      num_zeros = 0

      inQs = linspace((-scan_radius+wr),(scan_radius+wr),n_steps)
      ALLOCATE(iinQs(n_steps))

      do i = 1, n_steps
          x = inQs(i)!wr - scan_radius + (i - 0.5) * dx
          do j = 1, n_steps
              !y = -scan_radius + (j - 0.5) * dy
              y=-scan_radius+(REAL(j)/n_steps)*(scan_radius+
     $                 scan_radius)
              iinQs(j) = y
              delta = REAL(riccati(x,inQ_e,inQ_i,inpr,inc_beta,inds,
     $                            intau,inpe,iinQ=y))
              ric_val = abs(delta - deltaprime)
              deltas(i,j) = delta
              grad_x = 10!(REAL(riccati(x + dx,inQ_e,inQ_i,inpr,inc_beta,
      !$                           inds,intau,inpe,iinQ=y)) -
      !$                  REAL(riccati(x - dx,inQ_e,inQ_i,inpr,inc_beta,
      !$                            inds,intau,inpe,iinQ=y))) / (2 * dx)
              grad_y = 10!(REAL(riccati(x,inQ_e,inQ_i,inpr,inc_beta,inds,
      !$                            intau,inpe,iinQ=y + dy)) -
      !$                  REAL(riccati(x,inQ_e,inQ_i,inpr,inc_beta,inds,
      !$                            intau,inpe,iinQ=y - dy))) / (2 * dy)

              !WRITE(*,*) "Point", [x,y]
              !WRITE(*,*) "ric_val", ric_val
              !WRITE(*,*) "grad_x", grad_x
              !WRITE(*,*) "grad_y", grad_y
              !WRITE(*,*) "[dx,dy]", [dx,dy]


              if ((abs(ric_val) < zero_tol) .and. ((abs(grad_x) >
     $                grad_tol) .or. (abs(grad_y) > grad_tol))) then
                  ! Potential zero found, refine with a finer grid
                  !start = [x - 0.1, y - 0.1]
                  !end = [x + 0.1, y + 0.1]
                  !step = (end - start) / 20
                  !WRITE(*,*) "Potential zero found", [x,y]
                 ! WRITE(*,*) "dx", dx
                 ! WRITE(*,*) "dy", dy
                 ! WRITE(*,*) "step", step
                 ! WRITE(*,*) "scan_radius", scan_radius

                  x_s = linspace((x-0.1),(x+0.1),n_fine)
                  y_s = linspace((y-0.1),(y+0.1),n_fine)
                  found_zero = .false.
                  do k=1,n_fine
                    do w=1,n_fine
                          !x_s = start(1) + (k - 0.5) * step(1)
                          !y_s = start(2) + (k - 0.5) * step(2)

                          !WRITE(*,*) "fine point", [x_s(k),y_s(w)]
                          !WRITE(*,*) "y_s", y_s
                          !stop
                          !WRITE(*,*) "fine scan", k

                          ric_val = REAL(riccati(x_s(k),inQ_e,inQ_i,
     $                     inpr,inc_beta,inds,intau,inpe,
     $                     iinQ=y_s(w))) - deltaprime
                          !WRITE(*,*) "ric_val", ric_val

                          if (abs(ric_val) < abs(0.01*deltaprime)) then
                              ! Found a zero within tolerance, store it
                              num_zeros = num_zeros + 1
                              roots(num_zeros, 1) = x
                              roots(num_zeros, 2) = y
                              WRITE(*,*) "FOUND ROOT #", num_zeros
                              found_zero=.true.
                              !exit ! Exit the inner loop since zero is found
                          end if
                          if (found_zero==.true.) then
                            exit
                          endif
                    end do
                    if (found_zero==.true.) then
                        exit
                    endif
                  end do
              end if
          end do
      end do

      ! 1. Find the index of the (x, y) pair with x closest to wr
      n_points = size(roots, 1)
      min_dist1 = abs(roots(1, 1) - wr)
      idx_1 = 1

      do i = 2, n_points
          dist = abs(roots(i, 1) - wr)
          if (dist < min_dist1) then
              min_dist1 = dist
              idx_1 = i
          end if
      end do

      x_1 = roots(idx_1, 1)
      y_1 = roots(idx_1, 2)

      ! 2. Find the index of the next closest (x, y) pair with y more than 0.5 away from y_1
      min_dist2 = huge(min_dist2) ! Initialize to a very large value
      idx_2 = 0  ! Initialize to an invalid index

      do i = 1, n_points
          if (i /= idx_1 .and. abs(roots(i, 2) - y_1) > 0.5) then  ! Check if it's a different point and far enough in y
              dist = abs(roots(i, 1) - wr)
              if (dist < min_dist2) then
                  min_dist2 = dist
                  idx_2 = i
              end if
          end if
      end do

      if (idx_2 > 0) then  ! Check if a valid second point was found
          x_2 = roots(idx_2, 1)
          y_2 = roots(idx_2, 2)
      end if

      ! 3. Apply your specific conditions to select the final (x, y) pair
      if (idx_2 > 0 .and. abs(x_1 - x_2) < 0.1) then
          if (y_1 < 0 .or. y_2 < 0) then  ! Choose the pair with negative y
              if (y_1 < 0) then
                  growthrate_loc = x_1
                  growthrate = y_1
              else
                  growthrate_loc = x_2
                  growthrate = y_2
              end if
          end if
      else  ! If the x values are far apart, choose the first pair
          growthrate_loc = x_1
          growthrate = y_1
      end if
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      ! Allocate grid packing arrays and 2D complex deltas array
      !ALLOCATE(deltas(0:3+ReQ_num,0:ImQ_num))

      !DO i=0,ReQ_num+1
      !   DO j=0,ImQ_num
      !      iinQs(j)=Im_inQ_min+(REAL(j)/ImQ_num)*(Im_inQ_max-
      !$       Im_inQ_min)
      !      ! Run riccati() at each Q index to give delta
      !      deltas(i,j)=riccati(inQs(i),inQ_e,inQ_i,inpr,
      !$           inc_beta,inds,intau,inpe,iinQ=iinQs(j)) ! NOT USING GRID PACKING
      !   ENDDO
      !ENDDO

      RETURN
      END SUBROUTINE growthrate_search
c-----------------------------------------------------------------------
c     Subprogram 2. gamma_stability_scan
c     Run grid packed slayer stab. scan around omega_ExB and gamma axes
c-----------------------------------------------------------------------
      SUBROUTINE gamma_stability_scan(qval,inQ,inQ_e,inQ_i,inc_beta,
     $           inds,intau,inQ0,inpr,inpe,ReQ_num,ImQ_num,scan_radius,
     $           deltas,inQs,iinQs)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe
      REAL(r8), INTENT(IN) :: inQ ! REAL???
      INTEGER, INTENT(IN) :: qval,ReQ_num,ImQ_num,scan_radius

      ! Outputs
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE,
     $             INTENT(OUT) :: deltas
      REAL(r8), DIMENSION(:), ALLOCATABLE,
     $                        INTENT(OUT) :: inQs,iinQs

      ! Local variables
      INTEGER :: i,j,k
      REAL(r8) :: Re_inQ_min,Re_inQ_max,Im_inQ_min,Im_inQ_max
      CHARACTER(3) :: q_str
      CHARACTER(len=8) :: fmt ! format descriptor for stab file naming
      CHARACTER(len=8) :: x1 ! string for stab file naming
      INTEGER :: i1 ! integer for stab file naming
      !REAL(r8), DIMENSION(:,:), ALLOCATABLE :: inQs_left,inQs_right
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs_log
c-----------------------------------------------------------------------
c     Build exponential grid packing for stability scan, then run
c-----------------------------------------------------------------------
      ! Allocate grid packing arrays and 2D complex deltas array
      ALLOCATE(inQs(0:ReQ_num+1),iinQs(0:ImQ_num))
      !ALLOCATE(inQs_left(0:2+ReQ_num/2,0:2+ReQ_num/2))
      !ALLOCATE(inQs_right(0:2+ReQ_num/2,0:2+ReQ_num/2))
      !ALLOCATE(inQs_log(0:3+ReQ_num))
      ALLOCATE(deltas(0:3+ReQ_num,0:ImQ_num))

      Im_inQ_max=scan_radius ! max growth rate in scan, OPEN TO USER?
      Im_inQ_min=-scan_radius ! min growth rate in scan, OPEN TO USER?
      Re_inQ_max=scan_radius ! max growth rate in scan, OPEN TO USER?
      Re_inQ_min=-scan_radius ! min growth rate in scan, OPEN TO USER?

      ! Grid packing - right now going to Q +/- 0.2 -- OPEN TO USER?
      !inQs_left = powspace(inQ-Re_inQ_max,inQ,1, ! omega-3.0
      !$                        2+ReQ_num/2,"upper")
      !inQs_right = powspace(inQ,inQ+Re_inQ_max,1, ! omega+3.0
      !$                         2+ReQ_num/2,"lower")
      !inQs_log = (/inQs_left(1,1:2+ReQ_num/2),
      !$               inQs_right(1,2:1+ReQ_num/2)/)

      inQs = linspace((Re_inQ_min+inQ),(Re_inQ_max+inQ),ReQ_num+1)

      DO i=0,ReQ_num+1
         DO j=0,ImQ_num
            iinQs(j)=Im_inQ_min+(REAL(j)/ImQ_num)*(Im_inQ_max-
     $       Im_inQ_min)
            ! Run riccati() at each Q index to give delta
            deltas(i,j)=riccati(inQs(i),inQ_e,inQ_i,inpr,
     $           inc_beta,inds,intau,inpe,iinQ=iinQs(j)) ! NOT USING GRID PACKING
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     Alter qval to string to add to output filename
c-----------------------------------------------------------------------
      fmt = '(I5.5)' ! an integer of width 2 for q surface
      i1 = qval
      write (x1,fmt) i1 ! integer to string using a 'internal file'
      ! Write stability scan output file
      OPEN(UNIT=out_unit,FILE="slayer_stability_q"//
     $   TRIM(x1)//".out", STATUS="UNKNOWN")
      WRITE(out_unit,'(1x,4(a17))'),"RE(Q)",
     $     "IM(Q)","RE(delta)","IM(delta)"
      DO i=0,ReQ_num+1
         DO j=0,ImQ_num
            WRITE(out_unit,'(1x,4(es17.8e3))')
     $           inQs(i),iinQs(j),
     $           REAL(deltas(i,j)),AIMAG(deltas(i,j))
         ENDDO
      ENDDO
          CLOSE(out_unit)

      RETURN
      END SUBROUTINE gamma_stability_scan
c-----------------------------------------------------------------------
c     Subprogram 3. gamma_match
c     Loop stability across k rational surfaces
c-----------------------------------------------------------------------
      SUBROUTINE gamma_match(qval_arr,psi_n_rational,inQ_arr,inQ_e_arr,
     $                    inQ_i_arr,
     $                    inc_beta_arr,inds_arr,intau_arr,inQ0_arr,
     $                    inpr_arr,inpe_arr,omegas_arr,outer_delta_arr,
     $                    ReQ_num,ImQ_num,scan_radius,inQs,iinQs,
     $                    all_Re_deltas,all_inQs)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8), DIMENSION(:), INTENT(IN) :: inQ_e_arr,
     $                       inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,
     $                       inQ0_arr,inpr_arr,inpe_arr,omegas_arr,
     $                       inQ_arr,psi_n_rational
      INTEGER, DIMENSION(:), INTENT(IN) :: qval_arr
      REAL(r8), DIMENSION(:), INTENT(IN) :: outer_delta_arr
      INTEGER, INTENT(IN) :: ReQ_num,ImQ_num,scan_radius
      ! Outputs
      REAL(r8), DIMENSION(:), ALLOCATABLE :: growthrates
      ! Local variables
      INTEGER :: n_k ! Number of rational surfaces
      INTEGER :: k,w
      ! Local variables received from internal subroutines
      REAL(r8), DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: inQs,iinQs
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: all_RE_deltas
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: all_slices,all_inQs
      REAL(r8), DIMENSION(:), ALLOCATABLE :: slice
      REAL(r8) :: layer_Q, ImQ_gamma

      n_k = SIZE(qval_arr)
      ! Allocate growthrates arrays
      ALLOCATE(growthrates(n_k))

c-----------------------------------------------------------------------
c     Loop across rational surfaces
c-----------------------------------------------------------------------
      ! Summary: for each rational surface, run narrow stability scan
      ! for analysis, then slice out 1D array of growth rates at
      ! given omega_ExB (Q), then find growth rate corresponding
      ! to delta-deltaprime match
      DO k=1,n_k
         WRITE(*,*) "Scanning q=", qval_arr(k), " rational surface"
         ! Run stability scan
         CALL gamma_stability_scan(qval_arr(k),inQ_arr(k),inQ_e_arr(k),
     $            inQ_i_arr(k),inc_beta_arr(k),inds_arr(k),
     $            intau_arr(k),inQ0_arr(k),inpr_arr(k),inpe_arr(k),
     $            ReQ_num,ImQ_num,scan_radius,deltas,inQs,iinQs)

         layer_Q = inQ_arr(k) ! REAL???

         IF (k==1) THEN
            ALLOCATE(all_RE_deltas(SIZE(inQs),SIZE(iinQs),n_k))
            ALLOCATE(all_slices(SIZE(iinQs),n_k),
     $                all_inQs(SIZE(inQs),n_k))
         ENDIF

         all_RE_deltas(:,:,k) = REAL(deltas)
         all_inQs(:,k) = inQs
      ENDDO

      DEALLOCATE(deltas)

      RETURN
      END SUBROUTINE gamma_match

      END MODULE gslayer_mod
