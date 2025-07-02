      MODULE params_mod

      USE sglobal_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     calculate parameters.
c-----------------------------------------------------------------------
      SUBROUTINE params(n_e,t_e,t_i,omega,chi,dr_val,dgeo_val,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

      REAL(r8), INTENT(IN) :: n_e,t_e,t_i,omega,chi,dr_val,dgeo_val,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff

      LOGICAL, INTENT(IN) :: params_check

      REAL(r8) :: rho,b_l,v_a,Qconv,
     $            lbeta,tau_i,tau_h,tau_v
      REAL(r8) :: tau_ee_num,tau_ee_denom,tau_ee,sigma_par_1,
     $            sigma_par_2,sigma_par,tau_perp,Wd,vte,
     $            chi_par_smfp,chi_par_lmfp,chi_par
      INTEGER :: wit

      lnLamb = 24 + 3.0*LOG(10.0) - 0.5*LOG(n_e) + LOG(t_e)

      ! mu_i: ion mass ratio to proton
      tau= t_i/t_e ! ratio of ion to electron temperature 
      tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5/(n_e*lnLamb) ! ion colls.
      eta= 1.65e-9*lnLamb/(t_e/1e3)**1.5 ! spitzer resistivity (wesson)
      rho=(mu_i*m_p)*n_e ! mass density

      tau_ee_num = 6.0*SQRT(2.0)*(pi**1.5)*
     $             (eps0**2.0)*(m_e**0.5)*(t_e**1.5)
      tau_ee_denom = lnLamb*(chag**2.5)*n_e
      tau_ee = tau_ee_num / tau_ee_denom
      !tau_ee = 1.09e16 * (t_e/1.0e3)**1.5 / (n_e * lnLamb)


      sigma_par_1 = ( SQRT(2.0) + 13.0*(Zeff/4.0) ) / 
     $              (Zeff*(SQRT(2.0) + Zeff))
      sigma_par_2 = (n_e * (chag**2.0) * tau_ee) / m_e
      sigma_par = sigma_par_1*sigma_par_2

      b_l=(nr/mr)*rs*sval*bt/R0 ! characteristic magnetic field
      v_a=b_l/(mu0*rho)**0.5 ! alfven velocity
      rho_s=1.02e-4*(mu_i*t_e)**0.5/bt ! ion Lamour by elec. Temp.
      d_i = ( (mu_i*m_p)/(n_e * (chag**2) * mu0) )**0.5 ! collisionless ion skin depth

      tau_h=R0*(mu0*rho)**0.5/(nn*sval*bt) ! alfven time across surface
      !tau_r=mu0*rs**2.0/eta ! resistive time scale
      tau_r=mu0*(rs**2.0)*(sigma_par) ! R. Fitzpatrick resistive time scale
      tau_v=tau_r/pr !rho*rs**2.0/visc ! viscous time scale
      IF (ABS(chi) > 0.0) THEN
          tau_perp = ( rs**2.0 ) / chi
      ELSE
          tau_perp = 0.0
      END IF

      ! this one must be anomalous. calculated back from pr.
      visc= rho*rs**2.0/tau_v 
      
      lu=tau_r/tau_h ! Lundquist number 
      
      omega_e=-t_e/(bt*R0)*(1.0/l_n+1.0/l_t)*qval ! elec. diamag
      omega_i=t_i/(bt*R0)*(1.0/l_n+1.0/l_t)*qval ! ion diamag

      ! now calculate the main 7 normalized parameters.

      Qconv=lu**(1.0/3.0)*tau_h ! conversion to Qs based on Cole
      tauk = Qconv

      ! note Q depends on Qconv even if omega is fixed.     
      Q=Qconv*omega
      Q_e=-Qconv*omega_e
      Q_i=-Qconv*omega_i
      
      ! This is the most critical parameter
      ds=lu**(1.0/3.0)*rho_s/rs ! conversion based on Cole.

      lbeta=(5.0/3.0)*mu0*n_e*chag*(t_e+t_i)/bt**2.0
      c_beta=(lbeta/(1.0+lbeta))**0.5

      IF (ABS(tau_perp) > 0.0) THEN
          P_perp = tau_r / tau_perp ! perpendicular magnetic Prandtl number
      ELSE
          P_perp = 0.0
      END IF

      ! this is using Fitzpatrick's tau', we need tau eventually
      d_beta = c_beta*d_i  
      D_norm = (d_beta/rs)*(lu**(1.0/3.0))*(tau/(1+tau))**(0.5)

      delta_n=lu**(1.0/3.0)/rs ! norm factor for delta primes

      ! Calculate Delta_crit
      IF (ABS(dr_val) > 0.0) THEN
      vte = SQRT((2.0*(t_e*chag))/m_e)
      chi_par_smfp = (1.581*tau_ee*(vte**2.0))/
     $               (1.0+0.2535*Zeff)

      Wd = 0.1
      DO wit = 1,10

          chi_par_lmfp = (2.0*R0*vte)/(SQRT(pi)*nr*sval*Wd)
          chi_par = (chi_par_smfp*chi_par_lmfp)/
     $              (chi_par_smfp+chi_par_lmfp)
          Wd = SQRT(8.0)*((chi/chi_par)**0.25)*
     $         (1.0/SQRT((rs/R0)*sval*nr))
      END DO

      SELECT CASE(dc_type)
            CASE("lar")
              dc_tmp = 0.5*(-dr_val)*(pi**1.5)*((chi_par/chi)**0.25)*
     $                 ( (nr*sval)/(R0*rs) )**0.5
            CASE("rfitzp")
               dc_tmp = -(SQRT(2.0)*(pi**(1.5))*dr_val)/Wd
            CASE("toroidal")
               dc_tmp = 0.5*(-dr_val)*(pi**1.5)*
     $                  ((chi_par/chi)**0.25)*dgeo_val
            CASE default
               dc_tmp = 0.5*(-dr_val)*(pi**1.5)*
     $                  ((chi_par/chi)**0.25)*dgeo_val
      END SELECT

      ELSE
      dc_tmp = 0.0
      END IF

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
