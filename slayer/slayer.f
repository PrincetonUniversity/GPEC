c-----------------------------------------------------------------------
c     Slab LAYER based on linear drift MHD
c     SLAYER: main program
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     slayer.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM slayer

      USE sglobal_mod
      USE delta_mod, ONLY: riccati,riccati_del_s,riccati_out,
     $                     parflow_flag,PeOhmOnly_flag
      USE gslayer_mod
      USE layerinputs_mod

      IMPLICIT NONE

      CHARACTER(512) :: infile,ncfile
      INTEGER :: i,j,k,inum,jnum,knum,inn,count,
     $           Q_num,n_k,scan_radius
      INTEGER, DIMENSION(1) :: index

      LOGICAL :: params_flag,QPscan_flag,QPescan_flag,QPscan2_flag,
     $     QDscan2_flag,Qbscan_flag,Qscan_flag,
     $     onscan_flag,otscan_flag,ntscan_flag,nbtscan_flag,
     $     Pe_flag,verbose,ascii_flag,bin_flag,netcdf_flag,
     $     bal_flag,stability_flag,riccatiscan_flag,input_flag,
     $     params_check,stabscan_flag,read_eq,est_gamma_flag,
     $     match_gamma_flag,fitz_flag,br_th_flag
      REAL(r8) :: n_e,t_e,t_i,omega,omega0,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,dr_val
      REAL(r8) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds,intau,inlu
      REAL(r8) :: psi0,jxb,Q0,Q_sol,br_th,d_b,Residual
      COMPLEX(r8) :: delta,delta_n_p,dels_db,del_s,lar_gamma,
     $               tmp_gamma,ingamma,delta_prime

      REAL(r8) :: inQ_min,inQ_max,j_min,j_max,jpower,k_min,k_max,
     $     kpower,ing_step,ing_coarse,iing_coarse,delta_real,
     $     delta_imag,Qratio,chi

      INTEGER, DIMENSION(:), ALLOCATABLE :: mms,nns

      REAL(r8), DIMENSION(:), ALLOCATABLE :: jxbl,bal,
     $        prs,n_es,t_es,t_is,omegas,l_ns,l_ts,svals,qvals,
     $        bts,rss,R0s,mu_is,zeffs,Q_soll,br_thl,pes
      REAL(r8), DIMENSION(:), ALLOCATABLE :: Q_e_arr,Q_i_arr,
     $                   c_beta_arr,D_norm_arr,tau_arr,P_perp_arr,
     $                   omegas_arr,Q_arr,lu_arr,psi_n_rational,
     $                   d_beta_arr,Qconv_arr
      REAL(r8), DIMENSION(8) :: chi_prof
      INTEGER, DIMENSION(:), ALLOCATABLE :: qval_arr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs,
     $                       Re_deltaprime_arr,Im_deltaprime_arr,
     $                       gammafac_arr,delta_crit_arr
      REAL(r8), DIMENSION(:,:), ALLOCATABLE ::
     $     js,ks,psis,jxbs,Q_sols,br_ths
      REAL(r8) :: spot, slayer_inpr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Q_solss,br_thss
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal,outer_deltas,
     $                                     dels_db_arr,lar_gamma_arr,
     $                                     gamma_sol_arr,gamma_est_arr
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas

      NAMELIST/slayer_input/input_flag,infile,
     $     ncfile,params_flag,mm,nn,n_e,t_e,t_i,sval,bt,rs,R0,omega,
     $     l_t,l_n,qval,mu_i,zeff,dr_val,chi_prof,inpr,inpe,inQ,
     $     inQ_e,inQ_i,inc_beta,inds,intau,Q0,delta_prime,delta_n_p,
     $     ingamma
      NAMELIST/slayer_control/inum,jnum,knum,Q_num,scan_radius,
     $     read_eq,fitz_flag,QPscan_flag,Qscan_flag,QPescan_flag,
     $     Qbscan_flag,onscan_flag,otscan_flag,ntscan_flag,
     $     nbtscan_flag,parflow_flag,peohmonly_flag,Pe_flag,layfac
      NAMELIST/slayer_output/verbose,ascii_flag,bin_flag,netcdf_flag,
     $     est_gamma_flag,match_gamma_flag,stability_flag,
     $     stabscan_flag,br_th_flag,bal_flag
      NAMELIST/slayer_diagnose/riccati_out,riccatiscan_flag,
     $     params_check
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      mm=0.0
      nn=0.0
      mr = 0.0!real(mm,4)
      nr = 0.0!real(nn,4)
      n_e=0.0
      t_e=0.0
      t_i=0.0
      omega=0.0
      l_n=0.0
      l_t=0.0
      qval=0.0
      sval=0.0
      bt=0.0
      rs=0.0
      R0=0.0
      mu_i=0.0
      zeff=0.0
      dr_val=0.0
      inQ=0.0
      inQ_e=0.0
      inQ_i=0.0
      inpr=0.0
      chi_prof=0.0
      inpe=0.0
      inc_beta=0.0
      inds=0.0
      intau=0.0
      inlu=0.0
      Q0=0.0
      chi=0.0
      dr_val=0.0
      gamma_fac=0.0
      delta_prime=(0.0,0.0)
      delta_n_p=(0.0,0.0)
      ingamma=(0.0,0.0)
      inum=400 ! resolution to find error field thresholds.
      jnum=500 ! resolution for 2d scan along with Q,omega.
      knum=100 ! resolution for 2d scan alont with the other.
      Q_num=100 ! resolution for stab. scan along Re(Q) axis
      scan_radius = 2.0
      in_unit=1
      out_unit=2
      out2_unit=3
      out3_unit=4
      bin_unit=5
      bin_2d_unit=6
      input_unit=7
      read_eq=.FALSE.
      est_gamma_flag=.FALSE.
      match_gamma_flag=.FALSE.
      fitz_flag=.FALSE.
      QPscan_flag=.FALSE. ! scan (Q,P) space for delta and torque.
      QPescan_flag=.FALSE. ! scan (Q,Pe) space for delta and torque.
      Qbscan_flag=.FALSE. ! scan (Q,beta) space for delta and torque.
      onscan_flag=.FALSE. ! scan (omega,n) space for error fields.
      otscan_flag=.FALSE. ! scan (omega,t) space for error fields.
      ntscan_flag=.FALSE. ! scan (n,te) space for error fields.
      nbtscan_flag=.FALSE. ! scan (n,bt) space for error fields.
      layfac=0.02
      Qratio=0.5
      parflow_flag=.FALSE.
      PeOhmOnly_flag=.TRUE.
      Pe_flag=.FALSE.
      params_flag=.TRUE.
      input_flag=.FALSE.
      infile=""
      ncfile=""
      verbose=.TRUE.
      ascii_flag=.TRUE.
      bin_flag=.TRUE.
      netcdf_flag=.FALSE.
      riccati_out=.FALSE.
      riccatiscan_flag=.FALSE.
      params_check=.FALSE.
      bal_flag=.FALSE.
      stability_flag=.FALSE.
      stabscan_flag=.FALSE.
      br_th_flag=.FALSE.
c-----------------------------------------------------------------------
c     read slayer.in.
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

      IF (nn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)') nn
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)') nn
      ENDIF
c-----------------------------------------------------------------------
c     calculate parameters as needed.
c-----------------------------------------------------------------------
      IF (params_flag) THEN
         CALL params(n_e,t_e,t_i,omega,chi,dr_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)
         inQ=Q
         inQ_e=Q_e
         inQ_i=Q_i
         inc_beta=c_beta
         inds=ds
         intau=tau
         Q0=Q
      ELSE
         lu=inlu
      ENDIF
c-----------------------------------------------------------------------
c     calculate basic delta, torque, balance, error fields.
c-----------------------------------------------------------------------
      IF (.NOT. (match_gamma_flag)) THEN
      delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
      psi0=1.0/ABS(delta+delta_n_p) ! a.u.
      jxb=-AIMAG(1.0/(delta+delta_n_p)) ! a.u.
      WRITE(*,*)"delta=",delta
      WRITE(*,*)"psi0=",psi0
      WRITE(*,*)"jxb=",jxb
      ENDIF
c-----------------------------------------------------------------------
c     calculate parameters as needed.
c-----------------------------------------------------------------------
      IF (input_flag) THEN
         OPEN(UNIT=input_unit,FILE=infile,STATUS="old")
         READ(input_unit,*)inn
         ALLOCATE(mms(inn),nns(inn),prs(inn),
     $        n_es(inn),t_es(inn),t_is(inn),omegas(inn),
     $        l_ns(inn),l_ts(inn),qvals(inn),svals(inn),
     $        bts(inn),rss(inn),R0s(inn),mu_is(inn),zeffs(inn),
     $        Q_soll(inn),br_thl(inn))
         DO k=0,inn-1
            READ(input_unit,'(2(1x,I2),14(1x,e12.4))')
     $           mms(k),nns(k),prs(k),
     $           n_es(k),t_es(k),t_is(k),omegas(k),
     $           l_ns(k),l_ts(k),qvals(k),svals(k),
     $           bts(k),rss(k),R0s(k),mu_is(k),zeffs(k)
         ENDDO
         CLOSE(input_unit)
         
         DO k=0,inn-1
            WRITE(*,*)k
            mr=REAL(mms(k))
            nr=REAL(nns(k))
            inpr=prs(k)
            CALL params(n_es(k),t_es(k),t_is(k),omegas(k),chi,dr_val,
     $           l_ns(k),l_ts(k),qvals(k),svals(k),bts(k),rss(k),R0s(k),
     $           mu_is(k),zeffs(k),params_check)
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
        
            index=MAXLOC(bal)
            Q_soll(k)=inQs(index(1))
            br_thl(k)=sqrt(MAXVAL(bal)/lu*(svals(k)**2.0/2.0))*1e4

            WRITE(*,*)"Q_sol=",Q_soll(k)
            WRITE(*,*)"br_th=",br_thl(k)
            DEALLOCATE(inQs,deltal,jxbl,bal)         
         ENDDO
         OPEN(UNIT=out_unit,FILE="slayer_input_bal_n"//
     $      TRIM(sn)//".out",STATUS="UNKNOWN")
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
c     LAR (cylindrical) growthrates via restive layer thickness
c-----------------------------------------------------------------------
      IF (est_gamma_flag) THEN
      WRITE(*,*)"------------------------------------------"
      WRITE(*,*)">>> Estimating growth rate"

         IF (read_eq) THEN

            CALL build_inputs(infile,ncfile,chi_prof,qval_arr,
     $           psi_n_rational,lu_arr,Qconv_arr,Q_arr,Q_e_arr,Q_i_arr,
     $           c_beta_arr,d_beta_arr,D_norm_arr,tau_arr,
     $           P_perp_arr,omegas_arr,gammafac_arr,
     $           Re_deltaprime_arr,Im_deltaprime_arr,delta_crit_arr)

            n_k = SIZE(qval_arr)
         ELSE
            n_k = 1
            mr = mm
            nr = nn

            CALL params(n_e,t_e,t_i,omega,chi_prof(1),dr_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

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

            qval_arr = (/ qval /)
            omegas_arr = (/ omega /)
            Q_arr = (/ Q /)
            Q_e_arr = (/ Q_e /)
            Q_i_arr = (/ Q_i /)
            psi_n_rational = (/ 0.0 /)
            Re_deltaprime_arr = (/ REAL(delta_prime) /)
            Im_deltaprime_arr = (/ AIMAG(delta_prime) /)
            P_perp_arr = (/ P_perp /)
            tau_arr = (/ tau /)
            D_norm_arr = (/ D_norm /)
            d_beta_arr = (/ d_beta /)
            delta_eff = REAL(delta_prime)-delta_crit
            gamma_fac = (rs*delta_eff)/tau_r
            gammafac_arr = (/ gamma_fac /)
         END IF 

         ALLOCATE(gamma_est_arr(n_k),dels_db_arr(n_k))
         DO k=1,n_k
            WRITE(*,*) "Calculating growth rate estimate on q=",
     $       qval_arr(k)," rational surface"

            dels_db=riccati_del_s(Q_arr(k),Q_e_arr(k),
     $                   Q_i_arr(k),P_perp_arr(k),c_beta_arr(k),
     $                   D_norm_arr(k),tau_arr(k),
     $                   5.0*D_norm_arr(k))

            del_s = dels_db * d_beta_arr(k)

            gamma_est_arr(k) = gammafac_arr(k)/del_s
            dels_db_arr(k) = dels_db
            WRITE(*,*)"Growth rate estimate=",
     $                 REAL(gamma_est_arr(k))," [Hz]"

         ENDDO

         IF (.NOT. (match_gamma_flag)) THEN
            gamma_sol_arr = (/0./)
            CALL output_gamma(est_gamma_flag,qval_arr,
     $         omegas_arr,Q_arr,Q_e_arr,Q_i_arr,d_beta_arr,
     $         c_beta_arr,D_norm_arr,P_perp_arr,lu_arr,psi_n_rational,
     $         Re_deltaprime_arr,Im_deltaprime_arr,delta_crit_arr,
     $         dels_db_arr,gamma_sol_arr,gamma_est_arr,Qconv_arr,
     $         re_trace,im_trace)
         END IF
      ENDIF
c-----------------------------------------------------------------------
c     LAR (cylindrical) growthrates via restive layer thickness
c-----------------------------------------------------------------------
      IF (match_gamma_flag) THEN
         WRITE(*,*)"------------------------------------------"
         WRITE(*,*)">>> Calculating asymptotically matched growth rate"

         IF (read_eq) THEN

            CALL build_inputs(infile,ncfile,chi_prof,qval_arr,
     $          psi_n_rational,lu_arr,Qconv_arr,Q_arr,Q_e_arr,Q_i_arr,
     $          c_beta_arr,d_beta_arr,D_norm_arr,tau_arr,
     $          P_perp_arr,omegas_arr,gammafac_arr,
     $          Re_deltaprime_arr,Im_deltaprime_arr,delta_crit_arr)

            n_k = SIZE(qval_arr)

         ELSE
            n_k = 1

            ! Use namelist kinetic inputs instead of equilibrium files
            CALL params(n_e,t_e,t_i,omega,chi_prof(1),dr_val,
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

            qval_arr = (/ qval /)
            omegas_arr = (/ omega /)
            Q_arr = (/ Q /)
            Q_e_arr = (/ Q_e /)
            Q_i_arr = (/ Q_i /)
            psi_n_rational = (/ 0.0 /)
            Re_deltaprime_arr = (/ REAL(delta_prime) /)
            Im_deltaprime_arr = (/ AIMAG(delta_prime) /)
            delta_crit_arr = (/ delta_crit /)
            P_perp_arr = (/ P_perp /)
            tau_arr = (/ tau /)
            D_norm_arr = (/ D_norm /)
            d_beta_arr = (/ d_beta /)
            gammafac_arr = (/ gamma_fac /)
            c_beta_arr = (/ c_beta /)
            lu_arr = (/ lu /)
            Qconv_arr = (/ tauk /)
         END IF 

         ALLOCATE(gamma_sol_arr(n_k)) 

         DO k=1,n_k
            WRITE(*,*) "Calculating growth rates on q=",
     $       qval_arr(k)," rational surface"
   
            Q_e = Q_e_arr(k)
            Q_i = Q_i_arr(k)
            P_perp = P_perp_arr(k)
            tau = tau_arr(k)
            D_norm = D_norm_arr(k)
            c_beta = c_beta_arr(k)
            tauk = Qconv_arr(k)

            ! (Deltaprime - delta_crit)/S^1/3
            delta_eff = (Re_deltaprime_arr(k) - 
     $                   delta_crit_arr(k))/(lu_arr(k)**(1.0/3.0))
            pe = 0.0

            ALLOCATE(re_trace(100),im_trace(100))
            re_trace = 0.0
            im_trace = 0.0
            n_trace = 1

            IF (fitz_flag) THEN
               g_r = 0.0
               g_i = -Q_e
               re_trace(1) = 0.0
               im_trace(1) = -Q_e
            ELSE
               g_i = 0.0
               g_r = -Q_e
               im_trace(1) = 0.0
               re_trace(1) = -Q_e
            END IF

            CALL newton_root(g_r,g_i,1,fitz_flag)

            WRITE(*,*)"Success! growth rate = ",
     $                g_r/tauk," [Hz]"

            CALL shrink_array(re_trace, n_trace)
            CALL shrink_array(im_trace, n_trace)

            re_trace = re_trace/tauk
            im_trace = im_trace/tauk

            gamma_sol_arr(k) = g_r/tauk

         ENDDO 

         IF (.NOT. (est_gamma_flag)) THEN
            d_beta_arr = (/ 0. /)
            dels_db_arr = (/ 0. /)
         END IF

         IF (stabscan_flag) THEN
            WRITE(*,*)"------------------------------------------"
            WRITE(*,*)">>> Running Re(Q),Im(Q) scan, radius=",
     $                scan_radius

            ing_step = (2.0 * scan_radius) / (Q_num - 1)
            count = 0
      
            ALLOCATE(inQs(1:(Q_num+1)),iinQs(1:Q_num))
            ALLOCATE(deltas(1:(Q_num+1),1:Q_num))

            DO i = 1, (Q_num+1)
               DO j = 1, Q_num
                  ing_coarse = -scan_radius + (i - 1) * ing_step
                  iing_coarse = -scan_radius + (j - 1) * ing_step
                  ! Evaluate riccati function
                  g_tmp = CMPLX(ing_coarse,iing_coarse)
                  IF (fitz_flag) THEN
                     delta=riccati_f(g_tmp)
                  ELSE
                     delta=riccati(iing_coarse,Q_e,Q_i,P_perp,
     $                             c_beta,D_norm,tau,pe,
     $                             iinQ=ing_coarse)
                  END IF
                  inQs(i) = ing_coarse
                  iinQs(j) = iing_coarse
                  deltas(i,j) = delta
               ENDDO
            ENDDO

            OPEN(UNIT=out_unit,FILE="slayer_stability_n"//
     $         TRIM(sn)//".out", STATUS="UNKNOWN")
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

         ENDIF 

         CALL output_gamma(est_gamma_flag,qval_arr,
     $         omegas_arr,Q_arr,Q_e_arr,Q_i_arr,d_beta_arr,
     $         c_beta_arr,D_norm_arr,P_perp_arr,lu_arr,psi_n_rational,
     $         Re_deltaprime_arr,Im_deltaprime_arr,delta_crit_arr,
     $         dels_db_arr,gamma_sol_arr,gamma_est_arr,Qconv_arr,
     $         re_trace,im_trace)
         stop
      ENDIF
c-----------------------------------------------------------------------
c     TEST ANALYTIC SCAN IN GSLAYER.F, FOR TESTING ONLY
c-----------------------------------------------------------------------
      IF (br_th_flag) THEN

         WRITE(*,*)"running br_th scan"

         CALL params(n_e,t_e,t_i,omega,chi,dr_val,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)
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
         WRITE(*,*)"basic delta successful"
         delta_n_p=1e-2
         delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,
     $                   inpe)
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
         WRITE(*,*)"rotation scan"
         ! Scan of rotation
         inQ_max=10.0
         inQ_min=-10.0
         inum=200
         ALLOCATE(inQs(0:inum),deltal(0:inum),jxbl(0:inum),
     $                bal(0:inum))
         DO i=0,inum
            inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
            deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $                   inpr,inc_beta,inds,intau,inpe)
            jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
            bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
         ENDDO

         ! Identify the threshold from the maximum of the balance parameter
         index=MAXLOC(bal)
         Q_sol=inQs(index(1))
         !omega_sol=inQs(index(1))/Qconv
         br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))
         DEALLOCATE(inQs,deltal,jxbl,bal)

         WRITE(*,*)"allocating"
         qval_arr = (/ 3 /)
         inQs = (/ 1.0 /)
         inQs = (/ 1.0 /)

         n_k = SIZE(qval_arr)

         qval_arr = (/ 3 /)
         inQs = (/ 1.0 /)
         inQs = (/ 1.0 /)

         omegas_arr = (/ 0.0 /)
         Q_arr = (/ inQ /)
         Q_e_arr = (/ inQ_e /)
         Q_i_arr = (/ inQ_i /)
         psi_n_rational = (/ 0.0 /)
         Re_deltaprime_arr = (/ 0.0 /)
         Im_deltaprime_arr = (/ 0.0 /)
         P_perp_arr = (/ inpr /)

         WRITE(*,*)"allocations successful"

c         CALL slayer_netcdf_out(n_k,lar_gamma_eq_flag,lar_gamma_flag,
c     $                     stabscan_eq_flag,stabscan_flag,br_th_flag)
         stop
      ENDIF
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,5(a17))') "inQ","RE(delta)",
     $           "IM(delta)","jxb","bal"

            DO i=0,inum
               WRITE(out_unit,'(1x,5(es17.8e3))') inQs(i),
     $              REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
            ENDDO
            CLOSE(out_unit)
         ENDIF

         index=MAXLOC(bal)
         Q_sol=inQs(index(1))
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
     $         TRIM(sn)//".out", STATUS="UNKNOWN")
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         //TRIM(sn)//'.bin',
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         //TRIM(sn)//".bin",
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//".bin",
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
     $         TRIM(sn)//'.bin',
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
               
               CALL params(n_e*ks(j,k),t_e,t_i,omega*js(j,k),chi,dr_val,
     $              l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)
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
               index=MAXLOC(bal)
               Q_sols(j,k)=inQs(index(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_onscan_n"//
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//'.bin',
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
     $              omega*js(j,k),chi,dr_val,l_n,l_t,qval,sval,bt,
     $              rs,R0,mu_i,zeff,params_check)
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
               index=MAXLOC(bal)
               Q_sols(j,k)=inQs(index(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"t_e=",t_e*ks(j,k),"br_ths=",br_ths(j,k)
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_otscan_n"//
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//'.bin',
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
     $              chi,dr_val,l_n,l_t,qval,sval,bt,rs,R0,mu_i,
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
               index=MAXLOC(bal)
               Q_sols(j,k)=inQs(index(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_ntscan_n"//
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//'.bin',
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

             
               CALL params(n_e*ks(j,k),t_e,t_i,omega,chi,dr_val,
     $              l_n,l_t,qval,sval,bt*js(j,k),rs,R0,mu_i,zeff,
     $              params_check)
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
               index=MAXLOC(bal)
               Q_sols(j,k)=inQs(index(1))
               br_ths(j,k)=sqrt(MAXVAL(bal)/lu)*1e4
               WRITE(*,*)"br_ths=",br_ths(j,k)               
            ENDDO
         ENDDO

         IF (ascii_flag) THEN
            OPEN(UNIT=out_unit,FILE="slayer_nbtscan_n"//
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
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
     $         TRIM(sn)//'.bin',
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
     
