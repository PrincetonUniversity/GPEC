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
c
c
c     Adapted from
      FUNCTION dispersion_det(g_tmp,n_k,sl_in)
      
      COMPLEX(r8), INTENT(IN) :: g_tmp
      INTEGER, INTENT(IN) :: n_k
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
         DO k=1,2 !!! MAXING OUT AT 2X2
            Q_e = sl_in%Q_e_arr(k)
            Q_i = sl_in%Q_i_arr(k)
            P_perp = sl_in%P_perp_arr(k)
            P_tor = sl_in%P_tor_arr(k)
            tau = sl_in%tau_arr(k)
            D_norm = sl_in%D_norm_arr(k)
            c_beta = sl_in%c_beta_arr(k)
            tauk = sl_in%Qconv_arr(k)
            iota_e = Q_e / (Q_e - Q_i)

            delta_eff = (sl_in%Re_dp_arr(k) - 
     $          sl_in%d_crit_arr(k))/(sl_in%lu_arr(k)**(1.0/3.0))
     
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

      END MODULE gslayer_mod