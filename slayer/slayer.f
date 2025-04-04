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
      !USE params_mod
      USE delta_mod, ONLY: riccati,riccati_del_s,riccati_out,
     $                     parflow_flag,PeOhmOnly_flag

      ! FOR TESTING:
      USE gslayer_mod
      USE layerinputs_mod

      IMPLICIT NONE

      CHARACTER(512) :: infile,ncfile
      INTEGER :: i,j,k,inum,jnum,knum,inn,
     $           ReQ_num,ImQ_num,n_k,scan_radius
      INTEGER, DIMENSION(1) :: index

      LOGICAL :: params_flag,QPscan_flag,QPescan_flag,QPscan2_flag,
     $     QDscan2_flag,Qbscan_flag,Qscan_flag,
     $     onscan_flag,otscan_flag,ntscan_flag,nbtscan_flag,
     $     verbose,ascii_flag,bin_flag,netcdf_flag,
     $     bal_flag,stability_flag,riccatiscan_flag,input_flag,
     $     params_check,stabscan_eq_flag,stabscan_flag,
     $     lar_gamma_eq_flag,lar_gamma_flag,Pe_flag,
     $     br_th_flag,compress_deltas

      REAL(r8) :: n_e,t_e,t_i,omega,omega0,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff
      REAL(r8) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds,intau,inlu
      REAL(r8) :: psi0,jxb,Q0,Q_sol,br_th,tau_r,d_b
      COMPLEX(r8) :: delta,delta_n_p,dels_db,del_s,lar_gamma

      REAL(r8) :: inQ_min,inQ_max,j_min,j_max,jpower,k_min,k_max,
     $     kpower,
     $     Qratio

      INTEGER, DIMENSION(:), ALLOCATABLE :: mms,nns

      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs,jxbl,
     $        bal,
     $        prs,n_es,t_es,t_is,omegas,l_ns,l_ts,svals,qvals,
     $        bts,rss,R0s,mu_is,zeffs,Q_soll,br_thl,pes
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQ_e_arr,
     $                   inQ_i_arr,inc_beta_arr,inds_arr,
     $                   intau_arr,inQ0_arr,inpr_arr,
     $                   inpe_arr,omegas_arr,inQ_arr,
     $                   psi_n_rational,ind_beta_arr,D_beta_norm_arr
      REAL(r8), DIMENSION(8) :: inpr_prof
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: all_Re_deltas,
     $                           all_Im_deltas,all_roots
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: all_inQs
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: Re_deltas,Im_deltas

      INTEGER, DIMENSION(:), ALLOCATABLE :: qval_arr
      REAL(r8), DIMENSION(:), ALLOCATABLE ::
     $                       Re_deltaprime_arr,Im_deltaprime_arr,
     $                       all_growthrates,all_growthrate_locs
      REAL(r8), DIMENSION(:,:), ALLOCATABLE ::
     $     js,ks,psis,jxbs,Q_sols,br_ths,
     $     inQs_left,inQs_right,coarse_deltas
      REAL(r8) :: spot, slayer_inpr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Q_solss,br_thss
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal,outer_deltas,
     $                                     dels_db_arr,lar_gamma_arr
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas

      TYPE(result_type) :: results(8) ! Assuming max 8 rational surfaces

      REAL(r8), DIMENSION(:), ALLOCATABLE :: gammafac_arr,growthrates,
     $                                        growthrate_err

      NAMELIST/slayer_input/input_flag,infile,
     $     ncfile,params_flag,mm,nn,n_e,t_e,t_i,omega,l_n,l_t,
     $     qval,sval,bt,rs,R0,zeff,mu_i,inQ,inQ_e,
     $     inQ_i,inpr,inpr_prof,inpe,inc_beta,inds,
     $     intau,inlu,Q0,delta_n_p
      NAMELIST/slayer_control/inum,jnum,knum,ReQ_num,ImQ_num,
     $     scan_radius,QPscan_flag,QPscan2_flag,
     $     QPescan_flag,QDscan2_flag,Qbscan_flag,Qscan_flag,
     $     onscan_flag,otscan_flag,ntscan_flag,nbtscan_flag,
     $     layfac,Qratio,parflow_flag,peohmonly_flag,Pe_flag
      NAMELIST/slayer_output/verbose,ascii_flag,bin_flag,netcdf_flag,
     $     stability_flag,lar_gamma_eq_flag,lar_gamma_flag,
     $     stabscan_eq_flag,stabscan_flag,br_th_flag,compress_deltas,
     $     bal_flag
      NAMELIST/slayer_diagnose/riccati_out,riccatiscan_flag,
     $     params_check
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      mm=2
      nn=1
      mr = real(mm,4)
      nr = real(nn,4)
      n_e=1e19
      t_e=1e3
      t_i=1e3
      omega=1e4
      l_n=2e-1
      l_t=2e-1
      qval=2.0
      sval=2.0
      bt=1.0
      rs=0.5
      R0=1.0
      mu_i=2.0
      zeff=2.0
      inQ=20.0  ! Q=23.0 for DIII-D example.
      inQ_e=2.0 ! Q_e=2.0 for DIII-D example.
      inQ_i=-2.6 ! Q_i=-2.6 for DIII-D example.
      inpr=0.5 ! 0.5 for DIII-D example.
      inpr_prof=1.0
      inpe=0.1 !I added this
      inc_beta=0.7 ! c_beta=0.7 for DIII-D example.
      inds=6.0 ! 6.0 for DIII-D example.
      intau=1.0 ! 1.0 for DIII-D example.     
      inlu=1e8
      Q0=4.0     
      delta_n_p=(1e-2,1e-2)
      inum=400 ! resolution to find error field thresholds.
      jnum=500 ! resolution for 2d scan along with Q,omega.
      knum=100 ! resolution for 2d scan alont with the other.
      ReQ_num=350 ! resolution for stab. scan along Re(Q) axis
      ImQ_num=350 ! resolution for stab. scan along Im(Q) axis
      scan_radius = 3
      in_unit=1
      out_unit=2
      out2_unit=3
      out3_unit=4
      bin_unit=5
      bin_2d_unit=6
      input_unit=7
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
      stabscan_eq_flag=.FALSE.
      lar_gamma_flag=.FALSE.
      lar_gamma_eq_flag=.FALSE.
      br_th_flag=.FALSE.
      compress_deltas=.FALSE.
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
         CALL params(n_e,t_e,t_i,omega,
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
      delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
      psi0=1.0/ABS(delta+delta_n_p) ! a.u.
      jxb=-AIMAG(1.0/(delta+delta_n_p)) ! a.u.
      WRITE(*,*)"delta=",delta
      WRITE(*,*)"psi0=",psi0
      WRITE(*,*)"jxb=",jxb
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
            CALL params(n_es(k),t_es(k),t_is(k),omegas(k),
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
         WRITE(out_unit,'(1x,(2a17))'),"Q_sol","br_th"   
         
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
      IF (lar_gamma_eq_flag) THEN

         ! propagate inpr value to inpr_prof if inpr_prof is turned off
         IF (inpr_prof(1) < 0) THEN 
            inpr_prof = inpr
         END IF

         CALL build_inputs(infile,ncfile,inpr_prof,
     $               inpe,Pe_flag,qval_arr,psi_n_rational,
     $               inQ_arr,inQ_e_arr,inQ_i_arr,inc_beta_arr,
     $               inds_arr,ind_beta_arr,D_beta_norm_arr,
     $               intau_arr,inQ0_arr,inpr_arr,inpe_arr,
     $               omegas_arr,gammafac_arr,
     $               Re_deltaprime_arr,Im_deltaprime_arr)

         WRITE(*,*)"Safety factor values=",qval_arr
         WRITE(*,*)"inQ values=",inQ_arr
         WRITE(*,*)"Prantdl numbers=",inpr_arr
         WRITE(*,*)"Electron viscosities=",inpe_arr
         WRITE(*,*)"Omega ExB values=",omegas_arr
         WRITE(*,*)"outer region real deltaprimes=",Re_deltaprime_arr
         WRITE(*,*)"outer region imag deltaprimes=",Im_deltaprime_arr
         WRITE(*,*)"inQ_e_arr=",inQ_e_arr
         WRITE(*,*)"ind_beta_arr=",ind_beta_arr
         WRITE(*,*)"D_beta_norm_arr=",D_beta_norm_arr
         WRITE(*,*)"intau_arr=",intau_arr
         WRITE(*,*)"gammafac_arr=",gammafac_arr

         n_k = SIZE(qval_arr)

         ALLOCATE(lar_gamma_arr(n_k),dels_db_arr(n_k))

         DO k=1,n_k
            WRITE(*,*) "Calculating growth rate on q=", qval_arr(k),
     $       " rational surface"

            dels_db=riccati_del_s(inQ_arr(k),inQ_e_arr(k),
     $                   inQ_i_arr(k),inpr_arr(k),inc_beta_arr(k),
     $                   D_beta_norm_arr(k),intau_arr(k))

            del_s = dels_db * ind_beta_arr(k)
            lar_gamma = gammafac_arr(k)/del_s

            lar_gamma_arr(k) = lar_gamma
            dels_db_arr(k) = dels_db

         ENDDO
         WRITE(*,*)"Calling slayer_netcdf_out"
         WRITE(*,*)"qval_arr=",qval_arr
         WRITE(*,*)"omegas_arr=",omegas_arr
         WRITE(*,*)"inQ_arr=",inQ_arr
         WRITE(*,*)"inQ_e_arr=",inQ_e_arr
         WRITE(*,*)"inQ_i_arr=",inQ_i_arr
         WRITE(*,*)"ind_beta_arr=",ind_beta_arr
         WRITE(*,*)"D_beta_norm_arr=",D_beta_norm_arr
         WRITE(*,*)"inpr_arr=",inpr_arr
         WRITE(*,*)"psi_n_rational=",psi_n_rational
         WRITE(*,*)"lu_arr=",lu_arr
         WRITE(*,*)"Re_deltaprime_arr=",Re_deltaprime_arr
         WRITE(*,*)"Im_deltaprime_arr=",Im_deltaprime_arr
         WRITE(*,*)"dels_db_arr=",dels_db_arr
         WRITE(*,*)"lar_gamma_arr=",lar_gamma_arr

         br_th = 0.0

         CALL output_lar_gamma(lar_gamma_eq_flag,lar_gamma_flag,
     $         stabscan_eq_flag,stabscan_flag,br_th_flag,qval_arr,
     $         omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,ind_beta_arr,
     $         D_beta_norm_arr,inpr_arr,psi_n_rational,
     $         Re_deltaprime_arr,Im_deltaprime_arr,dels_db_arr,
     $         lar_gamma_arr)

         stop
      ENDIF
c-----------------------------------------------------------------------
c     LAR (cylindrical) growthrates via restive layer thickness
c-----------------------------------------------------------------------
      IF (lar_gamma_flag) THEN
         WRITE(*,*)"intau=",intau
         WRITE(*,*)"inQ=",inQ
         WRITE(*,*)"Prantdl numbers=",inpr
         WRITE(*,*)"Electron viscosities=",inpe
         WRITE(*,*)"inQ_e=",inQ_e
         WRITE(*,*)"inds=",inds
         WRITE(*,*)"inc_beta=",inc_beta

         WRITE(*,*)"Calculating LAR growth rate"

         CALL params(n_e,t_e,t_i,omega,
     $        l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

         !inQ=Q ! TAKEN FROM NAMELIST
         !inQ_e=Q_e ! TAKEN FROM NAMELIST
         !inQ_i=Q_i ! TAKEN FROM NAMELIST
         !inc_beta=c_beta
         !inds=ds
         !intau=tau
         !Q0=Q

         qval_arr = (/ qval /)
         omegas_arr = (/ omega /)
         inQ_arr = (/ inQ /)
         inQ_e_arr = (/ inQ_e /)
         inQ_i_arr = (/ inQ_i /)
         psi_n_rational = (/ 0.0 /)
         Re_deltaprime_arr = (/ REAL(delta_n_p) /)
         Im_deltaprime_arr = (/ AIMAG(delta_n_p) /)
         inpr_arr = (/ inpr /)
         D_beta_norm_arr = (/ D_beta_norm /)

         !D_beta_norm = inds ! NAMELIST

         WRITE(*,*)"D_beta_norm = ",D_beta_norm
         WRITE(*,*)"inds = ",inds

         WRITE(*,*)"lar_gamma inQ = ",inQ

         dels_db=riccati_del_s(inQ,inQ_e,inQ_i,inpr,inc_beta,d_beta,
     $                        intau)

         WRITE(*,*)"dels_db() call successful"

         eta = 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
         tau_r = mu0*rs**2.0/eta ! resistive time scale

         del_s = dels_db * d_beta
         lar_gamma = (REAL(delta_n_p)/tau_r) * (rs/del_s)

         ind_beta_arr = (/ d_beta /)
         dels_db_arr = (/ dels_db /)
         lar_gamma_arr = (/ lar_gamma /)

         WRITE(*,*)"dels_db=",dels_db
         WRITE(*,*)"lar_gamma=",lar_gamma

         WRITE(*,*)"slayer.f lar_gamma=",lar_gamma

         br_th = 0.0

         CALL output_lar_gamma(lar_gamma_eq_flag,lar_gamma_flag,
     $         stabscan_eq_flag,stabscan_flag,br_th_flag,qval_arr,
     $         omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,ind_beta_arr,
     $         D_beta_norm_arr,inpr_arr,psi_n_rational,
     $         Re_deltaprime_arr,Im_deltaprime_arr,dels_db_arr,
     $         lar_gamma_arr)

         stop
      ENDIF
c-----------------------------------------------------------------------
c     TEST GAMMA_MATCH IN GSLAYER.F, FOR TESTING ONLY
c-----------------------------------------------------------------------
      IF (stabscan_eq_flag) THEN
         !WRITE(*,*)"infile=",infile
         !WRITE(*,*)"ncfile=",ncfile

         ! propagate inpr value to inpr_prof if inpr_prof is turned off
         IF (inpr_prof(1) < 0) THEN 
            inpr_prof = inpr
         END IF

         CALL build_inputs(infile,ncfile,inpr_prof,
     $               inpe,Pe_flag,qval_arr,psi_n_rational,
     $               inQ_arr,inQ_e_arr,inQ_i_arr,inc_beta_arr,
     $               inds_arr,ind_beta_arr,D_beta_norm_arr,
     $               intau_arr,inQ0_arr,inpr_arr,inpe_arr,
     $               omegas_arr,gammafac_arr,
     $               Re_deltaprime_arr,Im_deltaprime_arr)

         WRITE(*,*)"Safety factor values=",qval_arr
         WRITE(*,*)"inQ values=",inQ_arr
         WRITE(*,*)"Prantdl numbers=",inpr_arr
         WRITE(*,*)"Electron viscosities=",inpe_arr
         WRITE(*,*)"Omega ExB values=",omegas_arr
         WRITE(*,*)"outer region real deltaprimes=",Re_deltaprime_arr
         WRITE(*,*)"outer region imag deltaprimes=",Im_deltaprime_arr
         WRITE(*,*)"inQ_e_arr=",inQ_e_arr
         WRITE(*,*)"inds_arr=",inds_arr
         WRITE(*,*)"intau_arr=",intau_arr
         n_k = SIZE(qval_arr)

         DO k=1,n_k
            WRITE(*,*) "Finding roots on q=", qval_arr(k),
     $       " rational surface"

            CALL growthrate_scan(qval_arr(k),lu_arr(k),inQ_arr(k),
     $         inQ_e_arr(k),inQ_i_arr(k),inc_beta_arr(k),inds_arr(k),
     $         intau_arr(k),inQ0_arr(k),inpr_arr(k),inpe_arr(k),
     $         scan_radius,reQ_num,compress_deltas,
     $         Re_deltaprime_arr(k),results(k))
            WRITE(*,*)"Exited growthrate_scan"

         ENDDO
         WRITE(*,*)"Calling slayer_netcdf_out"

         br_th = 0.0
         WRITE(*,*)"Successfully entered output_lar_gamma()"
         WRITE(*,*)"qval_arr = ",qval_arr

         CALL slayer_netcdf_out(SIZE(qval_arr),lar_gamma_eq_flag,
     $    lar_gamma_flag,stabscan_eq_flag,stabscan_flag,br_th_flag,
     $            qval_arr,omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,
     $            psi_n_rational,inpr_arr,br_th,Re_deltaprime_arr,
     $            Im_deltaprime_arr,dels_db_arr,lu_arr,ind_beta_arr,
     $            D_beta_norm_arr,lar_gamma_arr,inQs,iinQs,results)
c         CALL slayer_netcdf_out(n_k,lar_gamma_eq_flag,lar_gamma_flag,
c     $                     stabscan_eq_flag,stabscan_flag,br_th_flag)
        stop
      ENDIF
c-----------------------------------------------------------------------
c     TEST ANALYTIC SCAN IN GSLAYER.F, FOR TESTING ONLY
c-----------------------------------------------------------------------
      IF (stabscan_flag) THEN
         WRITE(*,*)"intau=",intau
         WRITE(*,*)"inQ=",inQ
         WRITE(*,*)"Prantdl numbers=",inpr
         WRITE(*,*)"Electron viscosities=",inpe
         WRITE(*,*)"inQ_e=",inQ_e
         WRITE(*,*)"inds=",inds
         WRITE(*,*)"inc_beta=",inc_beta

         WRITE(*,*)"running analytic scan"

         qval_arr = (/ qval /)
         omegas_arr = (/ omega /)
         inQ_arr = (/ inQ /)
         inQ_e_arr = (/ inQ_e /)
         inQ_i_arr = (/ inQ_i /)
         psi_n_rational = (/ 0.0 /)
         Re_deltaprime_arr = (/ REAL(delta_n_p) /)
         Im_deltaprime_arr = (/ AIMAG(delta_n_p) /)
         inpr_arr = (/ inpr /)

         CALL growthrate_scan(qval_arr(1),lu,inQ,inQ_e,
     $             inQ_i,inc_beta,inds,
     $             intau,inQ,inpr,inpe,
     $             scan_radius,reQ_num,compress_deltas,
     $             Re_deltaprime_arr(1),results(1))

         WRITE(*,*)"allocations successful"

         br_th = 0.0
c         CALL slayer_netcdf_out(1,lar_gamma_eq_flag,lar_gamma_flag,
c     $                     stabscan_eq_flag,stabscan_flag,br_th_flag)
         stop
      ENDIF
c-----------------------------------------------------------------------
c     TEST ANALYTIC SCAN IN GSLAYER.F, FOR TESTING ONLY
c-----------------------------------------------------------------------
      IF (br_th_flag) THEN

         WRITE(*,*)"running br_th scan"

         CALL params(n_e,t_e,t_i,omega,
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
         k=1
         IF (k==1) THEN
            ALLOCATE(all_RE_deltas(1,1,n_k))
            ALLOCATE(all_Im_deltas(1,1,n_k))
            ALLOCATE(all_inQs(1,n_k))
            ALLOCATE(all_growthrates(n_k))
            ALLOCATE(all_growthrate_locs(n_k))
         ENDIF
         all_Re_deltas(:,:,k) = 0.0
         all_Im_deltas(:,:,k) = 0.0
         all_inQs(:,k) = 0.0

         qval_arr = (/ 3 /)
         inQs = (/ 1.0 /)
         inQs = (/ 1.0 /)

         omegas_arr = (/ 0.0 /)
         inQ_arr = (/ inQ /)
         inQ_e_arr = (/ inQ_e /)
         inQ_i_arr = (/ inQ_i /)
         psi_n_rational = (/ 0.0 /)
         Re_deltaprime_arr = (/ 0.0 /)
         Im_deltaprime_arr = (/ 0.0 /)
         inpr_arr = (/ inpr /)

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
            WRITE(out_unit,'(1x,5(a17))'),"inQ","RE(delta)",
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
            WRITE(out_unit,'(1x,4(a17))'),"RE(Q)",
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
            WRITE(out_unit,'(1x,5(a17))'),"x","yphs","yamp",
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
            WRITE(out_unit,'(1x,6(a17))'),"Q","Pe","RE(delta)",
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
            WRITE(out_unit,'(1x,6(a17))'),"Q","Pr","RE(delta)",
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
            WRITE(out_unit,'(1x,6(a17))'),"Q","Pr","RE(delta)",
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
            WRITE(out_unit,'(1x,6(a17))'),"Q","Pr","RE(delta)",
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
            WRITE(out_unit,'(1x,6(a17))'),"Q","D","RE(delta)",
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
               
               CALL params(n_e*ks(j,k),t_e,t_i,omega*js(j,k),
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
            WRITE(out_unit,'(1x,6(a17))'),"Omega","Density",
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
     $              omega*js(j,k),l_n,l_t,qval,sval,bt,
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
            WRITE(out_unit,'(1x,6(a17))'),"Omega","Temperature",
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
            OPEN(UNIT=out_unit,FILE="slayer_ntscan_n"//
     $         TRIM(sn)//".out",STATUS="UNKNOWN")
            WRITE(out_unit,'(1x,6(a17))'),"Temperature","Density",
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

             
               CALL params(n_e*ks(j,k),t_e,t_i,omega,
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
            WRITE(out_unit,'(1x,4(a17))'),"Bt","Density",
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
     
