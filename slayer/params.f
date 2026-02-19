c=======================================================================
c     MODULE params_mod
c
c     Computes the normalised layer-physics parameters needed by the
c     SLAYER dispersion solver from dimensional equilibrium and
c     kinetic-profile inputs.
c
c     Subprograms contained:
c       1. params -- derive all normalised quantities (Q, ds, c_beta,
c              D_norm, P_perp, P_tor, lu, delta_crit, ...) for a
c              single rational surface and store them in sglobal_mod.
c=======================================================================
      MODULE params_mod

      USE sglobal_mod            ! SLAYER global scalars, types, constants

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. params.
c     Compute all derived layer-physics parameters for a single
c     rational surface from dimensional equilibrium / kinetic inputs.
c     Results are written to module-level variables in sglobal_mod
c     (tau, tau_r, tauk, lu, Q, Q_e, Q_i, ds, c_beta, d_beta,
c     D_norm, P_perp, P_tor, delta_n, dc_tmp, eta, visc, rho_s, ...).
c
c     BUG FLAG 1 -- `pr` (magnetic Prandtl number) and `pe` are read
c       from sglobal_mod but never set within this routine.  They must
c       be initialised elsewhere before calling params(), otherwise
c       `tau_v = tau_r / pr` will divide by zero or garbage.
c       Suggested fix: add pr/pe as INTENT(IN) arguments, or
c       document the required initialisation order.
c
c     BUG FLAG 2 -- Several debug WRITE statements print to stdout
c       unconditionally on every call.  For a public release these
c       should either be removed or guarded behind `params_check`.
c-----------------------------------------------------------------------
      SUBROUTINE params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

c --- arguments
      REAL(r8), INTENT(IN) :: n_e       ! electron density   [m^-3]
      REAL(r8), INTENT(IN) :: t_e       ! electron temperature [eV]
      REAL(r8), INTENT(IN) :: t_i       ! ion temperature      [eV]
      REAL(r8), INTENT(IN) :: omega     ! toroidal rotation    [rad/s]
      REAL(r8), INTENT(IN) :: dr_val    ! radial width dr at surface
      REAL(r8), INTENT(IN) :: dgeo_val  ! geometric delta (Shafranov shift factor)
      REAL(r8), INTENT(IN) :: l_n       ! density gradient length
      REAL(r8), INTENT(IN) :: l_t       ! temperature gradient length
      REAL(r8), INTENT(IN) :: qval      ! safety factor
      REAL(r8), INTENT(IN) :: sval      ! magnetic shear
      REAL(r8), INTENT(IN) :: bt        ! toroidal field  [T]
      REAL(r8), INTENT(IN) :: rs        ! minor radius    [m]
      REAL(r8), INTENT(IN) :: R0        ! major radius    [m]
      REAL(r8), INTENT(IN) :: mu_i      ! ion mass ratio to proton
      REAL(r8), INTENT(IN) :: zeff      ! effective charge
      REAL(r8), DIMENSION(3), INTENT(IN) :: chis
                                         ! (1) chi_perp [m^2/s]
                                         ! (2) chi_tor  [m^2/s]
                                         ! (3) kappa    [m^2/s]
      LOGICAL, INTENT(IN)  :: params_check  ! .TRUE. = print diagnostics

c --- local variables: basic plasma
      REAL(r8) :: rho              ! mass density [kg/m^3]
      REAL(r8) :: b_l              ! characteristic magnetic field [T]
      REAL(r8) :: v_a              ! Alfvén velocity [m/s]
      REAL(r8) :: lbeta            ! local beta-related quantity
c --- local variables: electron collision time (Braginskii)
      REAL(r8) :: tau_ee_num       ! numerator of tau_ee formula
      REAL(r8) :: tau_ee_denom     ! denominator of tau_ee formula
      REAL(r8) :: tau_ee           ! electron-electron collision time [s]
c --- local variables: parallel conductivity (Spitzer-Härm)
      REAL(r8) :: sigma_par_1      ! neoclassical correction factor
      REAL(r8) :: sigma_par_2      ! classical conductivity [1/(Ohm*m)]
      REAL(r8) :: sigma_par        ! parallel conductivity  [1/(Ohm*m)]
c --- local variables: timescales
      REAL(r8) :: tau_i            ! ion collision time      [s]
      REAL(r8) :: tau_h            ! Alfvén transit time     [s]
      REAL(r8) :: tau_v            ! viscous time            [s]
      REAL(r8) :: tau_tor          ! toroidal diffusion time [s]
      REAL(r8) :: tau_perp         ! perp. diffusion time    [s]
c --- local variables: delta_crit iteration
      REAL(r8) :: vte              ! thermal electron speed  [m/s]
      REAL(r8) :: chi_par_smfp     ! chi_par in short mfp limit
      REAL(r8) :: chi_par_lmfp     ! chi_par in long  mfp limit
      REAL(r8) :: chi_par          ! effective parallel thermal cond.
      REAL(r8) :: Wd               ! magnetic island width proxy
      INTEGER  :: wit              ! iteration counter
c --- local variables: unused intermediates
      REAL(r8) :: Qconv            ! (shadowed by module-level Qconv)
      REAL(r8) :: K_val, Csq       ! kappa/eta and composite quantity

c-----------------------------------------------------------------------
c     Coulomb logarithm, basic plasma quantities, and Spitzer
c     resistivity.
c-----------------------------------------------------------------------
      lnLamb = 24 + 3.0*LOG(10.0) - 0.5*LOG(n_e) + LOG(t_e)

      tau   = t_i / t_e                             ! T_i / T_e
      tau_i = 6.6e17*mu_i**0.5*(t_i/1e3)**1.5
     $        / (n_e*lnLamb)                         ! ion collision time [s]
      eta   = 1.65e-9*lnLamb / (t_e/1e3)**1.5       ! Spitzer resistivity (Wesson)
      rho   = (mu_i*m_p)*n_e                         ! mass density [kg/m^3]

c-----------------------------------------------------------------------
c     Electron-electron collision time (Braginskii) and Spitzer-Härm
c     parallel conductivity.
c-----------------------------------------------------------------------
      tau_ee_num   = 6.0*SQRT(2.0)*(pi**1.5)
     $              *(eps0**2.0)*(m_e**0.5)*(t_e**1.5)
      tau_ee_denom = lnLamb*(chag**2.5)*n_e
      tau_ee       = tau_ee_num / tau_ee_denom

      sigma_par_1 = ( SQRT(2.0) + 13.0*(Zeff/4.0) )
     $              / (Zeff*(SQRT(2.0) + Zeff))
      sigma_par_2 = (n_e * (chag**2.0) * tau_ee) / m_e
      sigma_par   = sigma_par_1 * sigma_par_2

c-----------------------------------------------------------------------
c     Characteristic field, Alfvén speed, length scales, and
c     fundamental timescales.
c     Note: mr, nr (rational m/n) and nn (toroidal mode number) are
c     module-level variables set before calling params().
c-----------------------------------------------------------------------
      b_l   = (nr/mr)*rs*sval*bt/R0               ! characteristic B [T]
      v_a   = b_l / (mu0*rho)**0.5                 ! Alfvén velocity [m/s]
      rho_s = 1.02e-4*(mu_i*t_e)**0.5 / bt        ! ion Larmor at T_e [m]
      d_i   = ( (mu_i*m_p) / (n_e*(chag**2)*mu0) )**0.5
                                                    ! ion skin depth [m]

      tau_h = R0*(mu0*rho)**0.5 / (nn*sval*bt)     ! Alfvén time [s]
      tau_r = mu0*(rs**2.0)*sigma_par               ! resistive time [s] (Fitzpatrick)
      tau_v = tau_r / pr                            ! viscous time [s] (BUG FLAG 1)

c     back-compute anomalous viscosity from tau_v
      visc = rho*rs**2.0 / tau_v

c     Lundquist number
      lu = tau_r / tau_h
      
c-----------------------------------------------------------------------
c     Diamagnetic frequencies and normalised Q parameters.
c     Qconv converts dimensional frequencies to the normalised Q
c     used in the SLAYER dispersion relation (Cole scaling).
c-----------------------------------------------------------------------
      omega_e = -t_e/(bt*R0)*(1.0/l_n + 1.0/l_t)*qval  ! electron diamagnetic [rad/s]
      omega_i =  t_i/(bt*R0)*(1.0/l_n + 1.0/l_t)*qval  ! ion diamagnetic     [rad/s]

      Qconv = lu**(1.0/3.0) * tau_h     ! frequency normalisation (Cole)
      tauk  = Qconv                      ! stored in sglobal_mod

      Q   = Qconv * omega               ! normalised rotation frequency
      Q_e = -Qconv * omega_e            ! normalised electron diamagnetic
      Q_i = -Qconv * omega_i            ! normalised ion diamagnetic

c     normalised ion Larmor radius (critical stability parameter)
      ds = lu**(1.0/3.0) * rho_s / rs

c-----------------------------------------------------------------------
c     Plasma beta and Prandtl-number-like transport ratios.
c-----------------------------------------------------------------------
      lbeta  = (5.0/3.0)*mu0*n_e*chag*(t_e+t_i) / bt**2.0
      c_beta = (lbeta / (1.0+lbeta))**0.5

c      kappa-based P_perp (JKP's definition, to implement later if desired)
c      K_val = chis(3) / eta
c      Csq   = c_beta**2.0 + (1.0 - c_beta**2.0)*K_val

c      IF (ABS(Csq) > 0.0) THEN
c        P_perp = Csq
c      ELSE
c        tau_perp = (rs**2.0) / chis(1)
c      END IF

c     effective perpendicular and toroidal Prandtl numbers
      tau_perp = (rs**2.0) / chis(1)
      P_perp   = tau_r / tau_perp           ! perp magnetic Prandtl number

      tau_tor = (rs**2.0) / chis(2)
      P_tor   = tau_r / tau_tor             ! toroidal magnetic Prandtl number

c-----------------------------------------------------------------------
c     Normalised beta-related width and Delta norm factor.
c-----------------------------------------------------------------------
c     d_beta uses Fitzpatrick's definition (tau' form)
      d_beta = c_beta * d_i
      D_norm = (d_beta/rs) * lu**(1.0/3.0)
     $         * (tau/(1+tau))**(0.5)

      delta_n = lu**(1.0/3.0) / rs  ! normalisation for Delta values

c-----------------------------------------------------------------------
c     Critical Deltaprime (dc_tmp) via iterative chi_parallel calculation.
c     The island-width Wd is iterated 10 times to converge the
c     short-mfp / long-mfp interpolation for chi_parallel.
c     dc_type (from sglobal_mod) selects the formula:
c       'lar'      -- cylindrical (Lutjens)
c       'rfitzp'   -- R. Fitzpatrick
c       'toroidal' -- toroidal geometry using dgeo_val
c       default    -- dc_tmp = 0
c
c     BUG FLAG 4 -- iteration count 10 is hardcoded with no
c       convergence check.  Consider adding a tolerance test.
c-----------------------------------------------------------------------
      IF (ABS(dr_val) > 0.0) THEN

          vte = SQRT((2.0*(t_e*chag)) / m_e)
          chi_par_smfp = (1.581*tau_ee*(vte**2.0))
     $                   / (1.0 + 0.2535*Zeff)

          Wd = 0.1              ! initial guess
          DO wit = 1, 10
              chi_par_lmfp = (2.0*R0*vte)
     $                       / (SQRT(pi)*nr*sval*Wd)
              chi_par = (chi_par_smfp*chi_par_lmfp)
     $                  / (chi_par_smfp + chi_par_lmfp)
              Wd = SQRT(8.0)*((chis(1)/chi_par)**0.25)
     $           * (1.0/SQRT((rs/R0)*sval*nr))
          END DO

          SELECT CASE(dc_type)
              CASE('lar')
                  dc_tmp = 0.5*(-dr_val)*(pi**1.5)
     $                     *((chi_par/chis(1))**0.25)
     $                     *( (nr*sval)/(R0*rs) )**0.5
              CASE('rfitzp')
                  dc_tmp = -(SQRT(2.0)*(pi**(1.5))*dr_val) / Wd
              CASE('toroidal')
                  dc_tmp = 0.5*(-dr_val)*(pi**1.5)
     $                     *((chi_par/chis(1))**0.25)*dgeo_val
              CASE default
                  dc_tmp = 0.0
          END SELECT

      ELSE
          dc_tmp = 0.0
      END IF

c-----------------------------------------------------------------------
c     optional diagnostics (guarded by params_check flag).
c-----------------------------------------------------------------------
      IF (params_check) THEN
         WRITE(*,*) 'eta    = ', eta
         WRITE(*,*) 'S      = ', lu
         WRITE(*,*) 'Q      = ', Q
         WRITE(*,*) 'Q_e    = ', Q_e
         WRITE(*,*) 'Q_i    = ', Q_i
         WRITE(*,*) 'ds     = ', ds
         WRITE(*,*) 'c_beta = ', c_beta
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE params

      END MODULE params_mod