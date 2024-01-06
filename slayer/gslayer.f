      MODULE gslayer_mod
      
      USE sglobal_mod
      USE delta_mod

      IMPLICIT NONE
      
      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. gpec_slayer.
c     run slayer to provide b_crit(ising).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_slayer(n_e,t_e,n_i,t_i,omega,omega_e,
     $   omega_i,qval,sval,bt,rs,R0,mms,nns,
     $     delta,psi0,jxb,omega_sol,br_th)

      REAL(r8),INTENT(IN) :: n_e,t_e,n_i,t_i,omega,omega_e,omega_i,
     $     qval,sval,bt,rs,R0
      INTEGER, INTENT(IN) :: mms,nns
      COMPLEX(r8),INTENT(OUT) :: delta,psi0
      REAL(r8),INTENT(OUT) :: jxb,omega_sol,br_th
   
      INTEGER :: i,inum
      INTEGER, DIMENSION(1) :: index

      REAL(r8) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds,intau,inlu
      REAL(r8) :: mrs,nrs,rho,b_l,v_a,Qconv,lbeta,Q0,delta_n_p,mu_i,zeff
      REAL(r8) :: inQ_min,inQ_max,Q_sol
      
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inQs,iinQs,jxbl,bal
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: deltal

      parflow_flag=.FALSE.
      PeOhmOnly_flag=.TRUE.

      mrs = real(mms,4)
      nrs = real(nns,4)

      ! initial parameters
      inQ=20.0  ! Q=23.0 for DIII-D example.
      inQ_e=2.0 ! Q_e=2.0 for DIII-D example.
      inQ_i=-2.6 ! Q_i=-2.6 for DIII-D example.
      inc_beta=0.7 ! c_beta=0.7 for DIII-D example.
      inds=6.0 ! 6.0 for DIII-D example.
      intau=1.0 ! 1.0 for DIII-D example.     

      mu_i=2.0
      zeff=2.0

      tau= t_i/t_e ! ratio of ion to electron temperature 
      tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
      eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
      rho=(mu_i*m_p)*n_e ! mass density

      b_l=(nrs/mrs)*nrs*sval*bt/R0 ! characteristic magnetic field
      v_a=b_l/(mu0*rho)**0.5 ! alfven velocity
      rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Lamour by elec. Temp.

      tau_h=R0*(mu0*rho)**0.5/(nns*sval*bt) ! alfven time across surface
      tau_r=mu0*rs**2.0/eta ! resistive time scale
      tau_v=tau_r/pr !rho*rs**2.0/visc ! viscous time scale
      
      ! this one must be anomalous. calculated back from pr.
      visc= rho*rs**2.0/tau_v 
      
      lu=tau_r/tau_h ! Lundquist number 

      Qconv=lu**(1.0/3.0)*tau_h ! conversion to Qs based on Cole
      
      ! note Q depends on Qconv even if omega is fixed.     
      Q=Qconv*omega
      Q_e=-Qconv*omega_e
      Q_i=-Qconv*omega_i

      ! This is the most critical parameter
      ds=lu**(1.0/3.0)*rho_s/rs ! conversion based on Cole.

      lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
      c_beta=(lbeta/(1.0+lbeta))**0.5

      delta_n=lu**(1.0/3.0)/rs ! norm factor for delta primes

      inQ=Q
      inQ_e=Q_e
      inQ_i=Q_i
      inc_beta=c_beta
      inds=ds
      intau=tau
      Q0=Q
      inpr=0.5 ! 0.5 for DIII-D example.
      inpe=0.0 !I added this
c-----------------------------------------------------------------------
c     calculate basic delta, torque, balance, error fields.
c-----------------------------------------------------------------------
      delta_n_p=1e-2
      delta=riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)
      psi0=1.0/ABS(delta+delta_n_p) ! a.u.
      jxb=-AIMAG(1.0/(delta+delta_n_p)) ! a.u.
!      WRITE(*,*)"delta=",delta
!      WRITE(*,*)"psi0=",psi0
!      WRITE(*,*)"jxb=",jxb
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
      ! just for diagnostics
      inQ_max=10.0
      inQ_min=-10.0
      
      riccati_flag=.FALSE.
      
      inum=200

      ALLOCATE(inQs(0:inum),deltal(0:inum),jxbl(0:inum),bal(0:inum)) 
      
      DO i=0,inum
         inQs(i)=inQ_min+(REAL(i)/inum)*(inQ_max-inQ_min)
         deltal(i)=riccati(inQs(i),inQ_e,inQ_i,
     $        inpr,inc_beta,inds,intau,inpe)
         jxbl(i)=-AIMAG(1.0/(deltal(i)+delta_n_p))
         bal(i)=2.0*inpr*(Q0-inQs(i))/jxbl(i)
      ENDDO
      OPEN(UNIT=out_unit,FILE="bal.out",STATUS="UNKNOWN")
      WRITE(out_unit,'(1x,5(a17))'),"inQ","RE(delta)",
     $     "IM(delta)","jxb","bal"    
      
      DO i=0,inum
         WRITE(out_unit,'(1x,5(es17.8e3))')
     $        inQs(i),REAL(deltal(i)),AIMAG(deltal(i)),jxbl(i),bal(i)
      ENDDO
      CLOSE(out_unit)
      
      index=MAXLOC(bal)
      Q_sol=inQs(index(1))
      omega_sol=inQs(index(1))/Qconv
      br_th=sqrt(MAXVAL(bal)/lu*(sval**2.0/2.0))
!      WRITE(*,*)"Q_sol=",Q_sol
!      WRITE(*,*)"br_th=",br_th
      DEALLOCATE(inQs,deltal,jxbl,bal)  

      RETURN
      END SUBROUTINE gpec_slayer

      END MODULE gslayer_mod


