      MODULE params_mod

      USE sglobal_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     calculate parameters.
c-----------------------------------------------------------------------
      SUBROUTINE params(n_e,t_e,t_i,omega,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

      REAL(r8), INTENT(IN) :: n_e,t_e,t_i,omega,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff

      LOGICAL, INTENT(IN) :: params_check

      REAL(r8) :: rho,b_l,v_a,Qconv,
     $            lbeta,tau_i,tau_h,tau_r,tau_v

      ! mu_i: ion mass ratio to proton
      tau= t_i/t_e ! ratio of ion to electron temperature 
      tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
      eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
      rho=(mu_i*m_p)*n_e ! mass density

      b_l=(nr/mr)*rs*sval*bt/R0 ! characteristic magnetic field
      v_a=b_l/(mu0*rho)**0.5 ! alfven velocity
      rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Lamour by elec. Temp.

      tau_h=R0*(mu0*rho)**0.5/(nn*sval*bt) ! alfven time across surface
      tau_r=mu0*rs**2.0/eta ! resistive time scale
      tau_v=tau_r/pr !rho*rs**2.0/visc ! viscous time scale
      
      ! this one must be anomalous. calculated back from pr.
      visc= rho*rs**2.0/tau_v 
      
      lu=tau_r/tau_h ! Lundquist number 
!      pr=tau_r/tau_v ! Prandtl number. Only place needs viscosity.
      
      omega_e=-t_e/(bt*R0)*(1.0/l_n+1.0/l_t)*qval ! elec. diamag
      omega_i=t_i/(bt*R0)*(1.0/l_n+1.0/l_t)*qval ! ion diamag

      ! now calculate the main 7 normalized parameters.

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

      ! quick diagnostics.
      IF (params_check) THEN
         WRITE(*,*)"eta=",eta
         WRITE(*,*)"S=",lu
         WRITE(*,*)"Q=",Q
         WRITE(*,*)"Q_e=",Q_e
         WRITE(*,*)"Q_i=",Q_i
         WRITE(*,*)"ds=",ds
         WRITE(*,*)"c_beta=",c_beta     
      ENDIF

      RETURN
      END SUBROUTINE params  
      
      END MODULE params_mod
