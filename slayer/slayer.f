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
      USE params_mod
      USE delta_mod, ONLY: riccati,riccati_out,
     $                     parflow_flag,PeOhmOnly_flag

      ! FOR TESTING:
      USE gslayer_mod

      IMPLICIT NONE

      CHARACTER(128) :: infile
      INTEGER :: i,j,k,inum,jnum,knum,inn,ReQ_num,ImQ_num,n_k
      INTEGER, DIMENSION(1) :: index

      LOGICAL :: params_flag,QPscan_flag,QPescan_flag,QPscan2_flag,
     $     QDscan2_flag,Qbscan_flag,Qscan_flag,
     $     onscan_flag,otscan_flag,ntscan_flag,nbtscan_flag,
     $     verbose,ascii_flag,bin_flag,netcdf_flag,
     $     bal_flag,stability_flag,riccatiscan_flag,input_flag,
     $     params_check,gamma_match_flag !FOR TESTING

      REAL(r8) :: n_e,t_e,t_i,omega,omega0,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff
      REAL(r8) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds,intau,inlu
      REAL(r8) :: psi0,jxb,Q0,Q_sol,br_th
      COMPLEX(r8) :: delta,delta_n_p

      REAL(r8) :: inQ_min,inQ_max,j_min,j_max,jpower,k_min,k_max,kpower,
     $     Qratio

      INTEGER, DIMENSION(:), ALLOCATABLE :: mms,nns

      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs,jxbl,bal,
     $        prs,n_es,t_es,t_is,omegas,l_ns,l_ts,svals,
     $        bts,rss,R0s,mu_is,zeffs,Q_soll,br_thl,pes,
     $        inQs_log
      REAL(r8), DIMENSION(:), ALLOCATABLE:: qvals,inQ_e_arr,
     $                       inQ_i_arr,inc_beta_arr,inds_arr,intau_arr,
     $                       inQ0_arr,inpr_arr,inpe_arr,omegas_arr
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: inQ_arr,outer_delta_arr
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: 
     $     js,ks,psis,jxbs,Q_sols,br_ths,
     $     inQs_left,inQs_right
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: Q_solss,br_thss
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal,outer_deltas
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: deltas

      REAL(r8), DIMENSION(:), ALLOCATABLE :: growthrates,
     $                                        growthrate_err

      NAMELIST/slayer_input/params_flag,input_flag,infile,
     $     mm,nn,n_e,t_e,t_i,omega,l_n,l_t,
     $     qval,sval,bt,rs,R0,zeff,mu_i,inQ,inQ_e,
     $     inQ_i,inpr,inpe,inc_beta,inds,intau,inlu,Q0,delta_n_p
      NAMELIST/slayer_control/inum,jnum,knum,QPscan_flag,QPscan2_flag,
     $     QPescan_flag,QDscan2_flag,Qbscan_flag,Qscan_flag,
     $     onscan_flag,otscan_flag,ntscan_flag,nbtscan_flag,
     $     layfac,Qratio,parflow_flag,peohmonly_flag
      NAMELIST/slayer_output/verbose,ascii_flag,bin_flag,netcdf_flag,
     $     stability_flag
      NAMELIST/slayer_diagnose/riccati_out,riccatiscan_flag,
     $     params_check,bal_flag,gamma_match_flag
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
      ReQ_num=100 ! resolution for stab. scan along Re(Q) axis
      ImQ_num=200 ! resolution for stab. scan along Im(Q) axis
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
      params_flag=.TRUE.
      input_flag=.FALSE.
      infile="input_params.dat"
      verbose=.TRUE.
      ascii_flag=.TRUE.
      bin_flag=.TRUE.
      netcdf_flag=.FALSE.
      riccati_out=.FALSE.
      riccatiscan_flag=.FALSE.
      params_check=.FALSE.
      bal_flag=.FALSE.
      stability_flag=.FALSE.
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
c     TEST GAMMA_MATCH IN GSLAYER.F, FOR TESTING ONLY
c-----------------------------------------------------------------------
      IF (gamma_match_flag) THEN
         n_k=2
         ALLOCATE(mms(n_k),nns(n_k),prs(n_k),pes(n_k),
     $        n_es(n_k),t_es(n_k),t_is(n_k),omegas(n_k),
     $        l_ns(n_k),l_ts(n_k),qvals(n_k),svals(n_k),
     $        bts(n_k),rss(n_k),R0s(n_k),mu_is(n_k),zeffs(n_k),
     $        outer_deltas(n_k))
         ! Approximate data for q=2 and q=3, for testing only
         ! These are of the form needed to output from
         ! upcoming build_inputs.f script, interfacing with equil.f
         ! Kinetic inputs (old)
         !mms = (/2, 2 /)
         !nns = (/1, 1 /)
         !prs = (/0.1, 0.005 /)
         !pes = (/0.0, 0.0 /)
         !n_es = (/7.765e+17, 9.999e+16 /)
         !t_es = (/26.0, 2.0 /)
         !t_is = (/26.0, 2.0 /)
         !l_ns = (/0.2, 0.2 /)
         !l_ts = (/0.1, 0.1 /)
         !svals = (/729.5, 500.5 /)
         !bts = (/1.0, 1.0 /)
         !rss = (/0.16, 0.2 /)
         !R0s = (/2.0, 2.0 /)
         !mu_is = (/2.0, 2.0 /)
         !zeffs = (/2.0, 2.0 /)
         qvals = (/2.0, 3.0 /)
         inQ_arr = (/ (2.6744E-003,0.0), (5.5874E-004,0.0) /)
         inQ_e_arr =  (/ 1.9127E-005, 4.6109E-007 /)
         inQ_i_arr = (/ -1.9127E-005, -4.6109E-007 /)
         inc_beta_arr= (/ 3.6809E-003, 3.6635E-004 /)
         inds_arr = (/ 1.4394, 0.1276 /)
         intau_arr = (/ 1.0, 1.0 /)
         inQ0_arr = (/ (2.6744E-003,0.0), (5.5874E-004,0.0) /)
         !eta_s =  (/ 6.6907E-006, 3.1361E-004 /)
         !S_s =   (/ 30696266.0, 1956419.0 /)
         inpr_arr = (/0.1, 0.005 /)
         inpe_arr = (/0.0, 0.0 /)
         omegas_arr = (/54530.0, 54530.0 /)
         outer_delta_arr = (/(27.15,0.1), (20.15,0.1) /)

         CALL gamma_match(qvals,inQ_arr,inQ_e_arr,inQ_i_arr,
     $                    inc_beta_arr,inds_arr,intau_arr,
     $                    inQ0_arr,inpr_arr,inpe_arr,
     $                    outer_delta_arr,omegas_arr,
     $                    ReQ_num,ImQ_num,growthrates,growthrate_err)

         WRITE(*,*)"growthrates=",growthrates
         WRITE(*,*)"growthrate error=",growthrate_err

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
     
