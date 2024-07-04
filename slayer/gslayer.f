      MODULE gslayer_mod
      
      USE sglobal_mod, ONLY: out_unit, r8, mu0, m_p, chag, lnLamb,
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
c     Subprogram 2. interpolate_slice_at_Q
c     Either extract or interpolate a 1D deltas slice at given Q
c-----------------------------------------------------------------------
      SUBROUTINE interpolate_slice_at_Q(deltas, Q, inQs_log,
     $                                      slice)
      ! Input
      REAL(r8), DIMENSION(:, :), INTENT(IN) :: deltas ! 2D delta array
      REAL(r8), DIMENSION(:), INTENT(IN) :: inQs_log ! Re(Q) axis
      REAL(r8), INTENT(IN) :: Q ! slice value (Q)
      ! Output
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: slice
      ! Local variables
      INTEGER :: NROWS, NCOLS, i_lower, i_upper, j
      REAL(r8) :: Q_lower,Q_upper,lower_slice,upper_slice

      slice = 0.0
      ! Extract array dimensions
      NROWS = SIZE(deltas, 1)
      NCOLS = SIZE(deltas, 2)
      ! Allocate slice array (same size as second dimension of deltas)
      ALLOCATE(slice(NCOLS))
      ! Find bracketing rows based on Q (assuming sorted x-axis)
      i_lower = 1
      DO WHILE (i_lower < NROWS .AND. inQs_log(i_lower+1) < Q)
        i_lower = i_lower + 1
      END DO

      ! Check if Q is included in inQs_log (it should be)
      IF (inQs_log(i_lower) - Q < 1.0E-6_r8) THEN !!! this threshold could be an issue?
            WRITE(*,*)"Indexing slice at Q!"
            DO j = 1, NCOLS
                slice(j) = deltas(i_lower, j)
            END DO
      ELSE

      !! Handle cases where Q is outside the range
      !    IF (i_lower == 1 .AND. Q < deltas(1, 1)) THEN
      !      i_lower = 1
      !      i_upper = 2
      !    ELSE IF (i_lower == NROWS) THEN
      !      i_lower = NROWS - 1
      !      i_upper = NROWS
      !    ELSE
      !      i_upper = i_lower + 1
      !    END IF

      ! Interpolate each element of the slice (along the columns)
          DO j = 1, NCOLS
            i_upper = i_lower+1

            lower_slice = deltas(i_lower, j)
            upper_slice = deltas(i_upper, j)
            Q_lower = deltas(i_lower, 1)
            Q_upper = deltas(i_upper, 1)

            slice(j) = (lower_slice + upper_slice )/2.0

      ! This linear interpolation isn't correct ???
            !slice(j) = lower_slice +
      !$             (Q - Q_lower) * (upper_slice -
      !$             lower_slice) / (Q_upper - Q_lower)
          END DO
      END IF

      RETURN
      END SUBROUTINE interpolate_slice_at_Q
c-----------------------------------------------------------------------
c     Subprogram 3. gamma_from_delta_match
c     Interpolate gamma corresponding to delta-deltaprime match
c-----------------------------------------------------------------------
      SUBROUTINE gamma_from_delta_match(slice, iinQs, deltap,
     $         ImQ_gamma)
      ! Inputs
      REAL(r8), DIMENSION(:), INTENT(IN) :: slice ! 1D gammas array, i.e. deltas(Q,:)
      REAL(r8), DIMENSION(:), INTENT(IN) :: iinQs ! 1D array of gammas, i.e. Im(Q)
      REAL(r8), INTENT(IN) :: deltap ! Target outer layer delta prime
      ! Output
      REAL(r8), INTENT(OUT) :: ImQ_gamma ! Matched gamma(delta = deltaprime), giving growth rate
      ! Local variables
      INTEGER :: n, i_lower, i_upper, i, j, i_mid
      REAL(r8) :: deltap_lower, deltap_upper, slope,
     $            lower_val, upper_val,n_poles,dx
      REAL(r8), DIMENSION(:), ALLOCATABLE :: temp_slice, temp_iinQs, ! Temporary arrays
     $                   grad_slice, pospole_gamma, negpole_gamma

      ! Array size check
      n = SIZE(slice)
      IF (SIZE(iinQs) /= n) THEN
        WRITE(*, *) 'ERROR: slice and iinQs arrays
     $       must have the same size.'
        RETURN
      END IF

      ! Allocate temporary arrays
      ALLOCATE(temp_slice(n), temp_iinQs(n))
      ALLOCATE(grad_slice(n-1)) ! gradient

      ImQ_gamma = 0.0

      ! Copy input arrays to temporary arrays
      temp_slice = slice
      temp_iinQs = iinQs

      ! GRADIENT THRESHOLD IS 50???
      n_poles = 1.0
      ALLOCATE(pospole_gamma(2),negpole_gamma(2))
      pospole_gamma = 1e+20
      negpole_gamma = 1e+20

      ! Calculate the gradient
      DO i = 1, n - 1
        dx = iinQs(i + 1) - iinQs(i)  ! Change in x
        grad_slice(i) = (slice(i + 1) - slice(i)) / dx  ! Slope (gradient)
      END DO

      ! Find first match
      DO i = 1, n-1
        ! is this a pole point?
        IF (((SIGN(1.0,slice(i)) /= SIGN(1.0,slice(i+1))) .AND.
     $     ABS(grad_slice(i) > 50))) THEN!.AND.
      !$         (SIGN(1.0,grad_slice(i)) /=
      !$          SIGN(1.0,grad_slice(i+1)))) THEN
          ! it is a pole
            IF (slice(i) > 0) THEN
                WRITE(*,*)"FOUND A POS POLE ",0.0
                pospole_gamma(n_poles) = iinQs(i)
                negpole_gamma(n_poles) = iinQs(i+1)
                n_poles = n_poles + 1
            ELSE
                WRITE(*,*)"FOUND A NEG POLE ",i
                negpole_gamma(n_poles) = iinQs(i)
                pospole_gamma(n_poles) = iinQs(i+1)
                n_poles = n_poles+1
            END IF

            WRITE(*,*)"n_poles=",n_poles

            IF ((SIGN(1.0,slice(i))<0.0 .AND. deltap<slice(i))
     $             .OR.
     $         (SIGN(1.0,slice(i))>0.0 .AND. deltap>slice(i))) THEN
                ImQ_gamma = iinQs(i)  ! stride value is huge and not captured by our grid, but will be close to the pole
            END IF
        ELSE
          ! Determine which element of the slice is lower and which is upper
            IF (slice(i) <= slice(i + 1)) THEN
                lower_val = slice(i)
                upper_val = slice(i + 1)
            ELSE
                lower_val = slice(i + 1)
                upper_val = slice(i)
            END IF

          ! Check if deltap is within the interval
            IF ((lower_val <= deltap) .AND.
     $          (deltap <= upper_val)) THEN
                !ImQ_gamma = (iinQs(i) + iinQs(i+1)) / 2  ! do a better linear interpolation here

                ! Linear interpolation (revised for robustness)
                IF (upper_val - lower_val < 1.0E-10_r8) THEN  ! Handle nearly equal values
                    ImQ_gamma = (iinQs(i) +
     $                 iinQs(i+1)) / 2
                ELSE
                  slope = (iinQs(i+1) -
     $                  iinQs(i)) / (upper_val - lower_val)
                  ImQ_gamma = iinQs(i) +
     $                  slope * (deltap - lower_val)
                END IF
            END IF
        END IF
      END DO

      IF (n_poles > 2) THEN
        WRITE(*,*)"Alert: more than two poles detected"
      END IF


      !WRITE(*,*)"temp_iinQs(i_lower)=",temp_iinQs(i_lower)
      !WRITE(*,*)"temp_iinQs(i_upper)=",temp_iinQs(i_upper)
      !WRITE(*,*)"NEWImQ_gamma=",ImQ_gamma
      !WRITE(*,*)"pospole_gamma=",pospole_gamma
      !WRITE(*,*)"negpole_gamma=",negpole_gamma

      !WRITE(*,*)"slice=",slice

      ! Deallocate temporary arrays
      DEALLOCATE(temp_slice, temp_iinQs)

      RETURN
      END SUBROUTINE gamma_from_delta_match
c-----------------------------------------------------------------------
c     Subprogram 4. gamma_stability_scan
c     Run grid packed slayer stab. scan around omega_ExB and gamma axes
c-----------------------------------------------------------------------
      SUBROUTINE gamma_stability_scan(qval,inQ,inQ_e,inQ_i,inc_beta,
     $            inds,intau,inQ0,inpr,inpe,ReQ_num,ImQ_num,deltas,
     $            inQs_log,iinQs)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8),INTENT(IN) :: qval,inQ_e,inQ_i,inc_beta,inds,
     $     intau,inQ0,inpr,inpe
      REAL(r8), INTENT(IN) :: inQ ! REAL???
      INTEGER, INTENT(IN) :: ReQ_num,ImQ_num

      ! Outputs
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE,
     $             INTENT(OUT) :: deltas
      REAL(r8), DIMENSION(:), ALLOCATABLE,
     $                        INTENT(OUT) :: inQs_log,iinQs

      ! Local variables
      INTEGER :: i,j,k
      REAL(r8) :: inQ_min,inQ_max
      CHARACTER(3) :: q_str
      CHARACTER(len=8) :: fmt ! format descriptor for stab file naming
      CHARACTER(len=8) :: x1 ! string for stab file naming
      INTEGER :: i1 ! integer for stab file naming
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: inQs_left,inQs_right
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs
c-----------------------------------------------------------------------
c     Build exponential grid packing for stability scan, then run
c-----------------------------------------------------------------------
      ! Allocate grid packing arrays and 2D complex deltas array
      ALLOCATE(inQs(0:ReQ_num),iinQs(0:ImQ_num))
      ALLOCATE(inQs_left(0:2+ReQ_num/2,0:2+ReQ_num/2))
      ALLOCATE(inQs_right(0:2+ReQ_num/2,0:2+ReQ_num/2))
      ALLOCATE(inQs_log(0:3+ReQ_num))
      ALLOCATE(deltas(0:3+ReQ_num,0:ImQ_num))

      inQ_max=3.0 ! max growth rate in scan, OPEN TO USER?
      inQ_min=-3.0 ! min growth rate in scan, OPEN TO USER?

      ! Grid packing - right now going to Q +/- 0.2 -- OPEN TO USER?
      inQs_left = powspace(inQ-0.5,inQ,1, ! omega-0.5
     $                        2+ReQ_num/2,"upper")
      inQs_right = powspace(inQ,inQ+0.5,1, ! omega+0.5
     $                         2+ReQ_num/2,"lower")
      inQs_log = (/inQs_left(1,1:2+ReQ_num/2),
     $               inQs_right(1,2:1+ReQ_num/2)/)
      !WRITE(*,*)"inQs_log=",inQs_log

      DO i=0,ReQ_num+1
         DO j=0,ImQ_num
            ! Getting rid of "inQs" weirdly broke things??
            inQs(i)=inQ_min+(REAL(i)/ReQ_num)*(inQ_max-inQ_min)
            iinQs(j)=inQ_min+(REAL(j)/ImQ_num)*(inQ_max-inQ_min)
            ! Run riccati() at each Q index to give delta
            deltas(i,j)=riccati(inQs_log(i),inQ_e,inQ_i,inpr,
     $           inc_beta,inds,intau,inpe,iinQ=iinQs(j))
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
     $           inQs_log(i),iinQs(j),
     $           REAL(deltas(i,j)),AIMAG(deltas(i,j))
         ENDDO
      ENDDO
          CLOSE(out_unit)

      DEALLOCATE(inQs_left,inQs_right)

      RETURN
      END SUBROUTINE gamma_stability_scan
c-----------------------------------------------------------------------
c     Subprogram 5. gamma_match
c     Loop stability scans and gamma matches across k rational surfaces
c-----------------------------------------------------------------------
      SUBROUTINE gamma_match(qval_arr,psi_n_rational,inQ_arr,inQ_e_arr,
     $ inQ_i_arr,
     $                    inc_beta_arr,inds_arr,intau_arr,inQ0_arr,
     $                    inpr_arr,inpe_arr,omegas_arr,outer_delta_arr,
     $                    ReQ_num,ImQ_num,growthrates,growthrate_err)
c-----------------------------------------------------------------------
c     Declarations
c-----------------------------------------------------------------------
      ! Inputs
      REAL(r8), DIMENSION(:), INTENT(IN) :: qval_arr,inQ_e_arr,
     $                       inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,
     $                       inQ0_arr,inpr_arr,inpe_arr,omegas_arr,
     $                       inQ_arr,psi_n_rational
      REAL(r8), DIMENSION(:), INTENT(IN) :: outer_delta_arr
      INTEGER, INTENT(IN) :: ReQ_num,ImQ_num
      ! Outputs
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::
     $                                   growthrates, growthrate_err
      ! Local variables
      INTEGER :: n_k ! Number of rational surfaces
      INTEGER :: k,w
      ! Local variables received from internal subroutines
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs_log,iinQs
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: all_RE_deltas
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: all_slices
      REAL(r8), DIMENSION(:), ALLOCATABLE :: slice
      REAL(r8), DIMENSION(3) :: Q_range ! (-Q_err, Q, +Q_err )
      REAL(r8), DIMENSION(3) :: growthrate_range ! (-gamma, +gamma )
      REAL(r8) :: layer_Q, ImQ_gamma

      n_k = SIZE(qval_arr)
      ! Allocate growthrates arrays
      ALLOCATE(growthrates(n_k))
      ALLOCATE(growthrate_err(n_k))

c-----------------------------------------------------------------------
c     Loop across rational surfaces
c-----------------------------------------------------------------------
      ! Summary: for each rational surface, run narrow stability scan
      ! for analysis, then slice out 1D array of growth rates at
      ! given omega_ExB (Q), then find growth rate corresponding
      ! to delta-deltaprime match
      DO k=1,n_k
         WRITE(*,*)"layer #: ",k
         ! Run stability scan
         CALL gamma_stability_scan(qval_arr(k),inQ_arr(k),inQ_e_arr(k),
     $            inQ_i_arr(k),inc_beta_arr(k),inds_arr(k),
     $            intau_arr(k),inQ0_arr(k),inpr_arr(k),inpe_arr(k),
     $            ReQ_num,ImQ_num,deltas,inQs_log,iinQs)

         layer_Q = inQ_arr(k) ! REAL???

         ! Hardcoding rudimentary +/- 10% omega_ExB errorbars
         Q_range = (/0.9*layer_Q, 1.1*layer_Q, layer_Q/)

         ! Calculate growth rate +/- omega_ExB = Q errorbars
         DO w=1,3
            IF (w==3) THEN

                ! Slice out growth rates at layer_Q (Re(Q))
                CALL interpolate_slice_at_Q(REAL(deltas),
     $                      Q_range(w), inQs_log, slice)
                ! Match delta to delta prime to obtain growth rate
                CALL gamma_from_delta_match(slice, iinQs,
     $                             outer_delta_arr(k),
     $                             ImQ_gamma)
                CALL gamma_from_delta_match(slice, iinQs,
     $                             outer_delta_arr(k),
     $                             ImQ_gamma)
             ! Qconv = Q / omega_ExB
             ! gamma = Im(Q) / Qconv
                growthrates(k) = ImQ_gamma / (layer_Q / omegas_arr(k))
            !ELSE
            !   growthrate_range(w) = ImQ_gamma / (layer_Q / omegas(k))
            ENDIF

         ENDDO
         growthrate_err(k) = ABS(growthrate_range(2) -
     $                            growthrate_range(1))

         ALLOCATE(all_RE_deltas(SIZE(inQs_log),SIZE(iinQs),n_k))
         ALLOCATE(all_slices(SIZE(iinQs),n_k))

         all_RE_deltas(:,:,k) = REAL(deltas)
         all_slices(:,k) = slice

         DEALLOCATE(slice) ! Free memory after use for each layer

      ENDDO

      CALL slayer_netcdf_out(n_k,SIZE(inQs_log),SIZE(iinQs),qval_arr,
     $ inQs_log,iinQs,growthrates,omegas_arr,inQ_arr,psi_n_rational,
     $ all_Re_deltas,all_slices)

      DEALLOCATE(deltas)

      WRITE(*,*)"inQs_log len ",SIZE(inQs_log)
      WRITE(*,*)"iinQs len ",SIZE(iinQs)

      RETURN
      END SUBROUTINE gamma_match

      END MODULE gslayer_mod
