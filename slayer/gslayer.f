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
         WRITE(out_unit,'(1x,5(a17))') "inQ","RE(delta)",
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
c     Subprogram 3. scan_grid
c     Run stability scan on real and imaginary rotation axes
c-----------------------------------------------------------------------
      SUBROUTINE output_lar_gamma(lar_gamma_eq_flag,lar_gamma_flag,
     $       stabscan_eq_flag,stabscan_flag,br_th_flag,qval_arr,
     $       omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,ind_beta_arr,
     $       D_beta_norm_arr,inpr_arr,psi_n_rational,Re_deltaprime_arr,
     $       Im_deltaprime_arr,dels_db_arr,lu_arr,lar_gamma_arr)

      ! Declarations (include necessary type declarations from original code)
      LOGICAL, INTENT(IN) :: lar_gamma_eq_flag,lar_gamma_flag,
     $         stabscan_eq_flag,stabscan_flag,br_th_flag

      INTEGER, INTENT(IN), DIMENSION(:), ALLOCATABLE :: qval_arr
      REAL(r8), INTENT(IN), DIMENSION(:), ALLOCATABLE :: omegas_arr,
     $      inQ_arr,inQ_e_arr,inQ_i_arr,psi_n_rational,
     $      Re_deltaprime_arr,Im_deltaprime_arr,inpr_arr,ind_beta_arr,
     $      D_beta_norm_arr,lu_arr
  
      COMPLEX(r8), INTENT(IN), DIMENSION(:), ALLOCATABLE :: dels_db_arr,
     $                                         lar_gamma_arr

      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs
      TYPE(result_type) :: results(8)

      REAL(r8) :: br_th = 0

      WRITE(*,*)"Successfully entered output_lar_gamma()"

      inQs = (/0.0/)
      iinQs = (/0.0/)

      CALL slayer_netcdf_out(SIZE(qval_arr),lar_gamma_eq_flag,
     $    lar_gamma_flag,stabscan_eq_flag,stabscan_flag,br_th_flag,
     $            qval_arr,omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,
     $            psi_n_rational,inpr_arr,br_th,Re_deltaprime_arr,
     $            Im_deltaprime_arr,dels_db_arr,lu_arr,ind_beta_arr,
     $            D_beta_norm_arr,lar_gamma_arr,inQs,iinQs,results)
   
      END SUBROUTINE output_lar_gamma
c-----------------------------------------------------------------------
c     Subprogram 2. growthrate_scan
c     Set up and iterate stability scans if no match is found
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
          ! Evaluate riccati function
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
c     Increase scan array size if necessary
c-----------------------------------------------------------------------
      SUBROUTINE grow_array(arr, old_size, new_size)
          REAL(r8), ALLOCATABLE, INTENT(INOUT) :: arr(:)
          INTEGER, INTENT(IN) :: old_size,new_size
          REAL(r8), ALLOCATABLE :: temp(:)

          ALLOCATE(temp(new_size))
          temp(1:old_size) = arr(1:old_size)
          CALL move_alloc(temp, arr)
      END SUBROUTINE grow_array
      END MODULE gslayer_mod
