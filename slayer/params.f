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
c       2. slayer_d_norm -- D_norm from D_fac and the diamagnetic
c              frequencies.
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
c     Note: `pr` (magnetic Prandtl number) is read from sglobal_mod
c     and must be set by the caller (e.g. pr = inpr from namelist).
c     If pr == 0, tau_v and visc are set to zero (no viscosity).
c-----------------------------------------------------------------------
      SUBROUTINE params(n_e,t_e,t_i,omega,chis,dr_val,dgeo_val,
     $     l_n,l_t,qval,sval,bt,rs,R0,mu_i,zeff,params_check)

c --- arguments
      REAL(r8), INTENT(IN) :: n_e       ! electron density   [m^-3]
      REAL(r8), INTENT(IN) :: t_e       ! electron temperature [eV]
      REAL(r8), INTENT(IN) :: t_i       ! ion temperature      [eV]
      REAL(r8), INTENT(IN) :: omega     ! toroidal rotation    [rad/s]
      REAL(r8), INTENT(IN) :: dr_val    ! radial width dr at surface
      REAL(r8), INTENT(IN) :: dgeo_val  ! unused; kept for the interface
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
      REAL(r8) :: Wd_new            ! updated Wd for convergence check
      INTEGER  :: wit              ! iteration counter
c --- local variables: unused intermediates
      REAL(r8) :: K_val, Csq       ! kappa/eta and composite quantity

c-----------------------------------------------------------------------
c     stop on inputs that would give NaN or Inf below.
c-----------------------------------------------------------------------
      IF (.NOT. (n_e > 0.0_r8 .AND. t_e > 0.0_r8 .AND.
     $     t_i > 0.0_r8)) THEN
         WRITE(*,*) 'ERROR: params needs positive n_e, t_e and t_i,'
         WRITE(*,*) '  got n_e =',n_e,' [m^-3], t_e =',t_e,
     $        ' [eV], t_i =',t_i,' [eV]'
         WRITE(*,*) '  at q =',qval,'. Check the kinetic profiles.'
         STOP 'params: non-positive density or temperature'
      END IF
      IF (sval == 0.0_r8) THEN
         WRITE(*,*) 'ERROR: params needs a nonzero magnetic shear,'
         WRITE(*,*) '  got sval = 0 at q =',qval,'. The Alfven time'
         WRITE(*,*) '  and the layer widths are undefined.'
         STOP 'params: zero magnetic shear'
      END IF
c-----------------------------------------------------------------------
c     Coulomb logarithm, basic plasma quantities, and Spitzer
c     resistivity.
c-----------------------------------------------------------------------
      lnLamb = 24 + 3.0_r8*LOG(10.0_r8) - 0.5_r8*LOG(n_e) + LOG(t_e)

      tau   = t_i / t_e                             ! T_i / T_e
      tau_i = 6.6e17_r8*mu_i**0.5_r8*(t_i/1e3_r8)**1.5_r8
     $        / (n_e*lnLamb)                         ! ion collision time [s]
      eta   = 1.65e-9_r8*lnLamb / (t_e/1e3_r8)**1.5_r8       ! Spitzer resistivity (Wesson)
      rho   = (mu_i*m_p)*n_e                         ! mass density [kg/m^3]

c-----------------------------------------------------------------------
c     Electron-electron collision time (Braginskii) and Spitzer-Härm
c     parallel conductivity.
c-----------------------------------------------------------------------
      tau_ee_num   = 6.0_r8*SQRT(2.0_r8)*(pi**1.5_r8)
     $              *(eps0**2.0_r8)*(m_e**0.5_r8)*(t_e**1.5_r8)
      tau_ee_denom = lnLamb*(chag**2.5_r8)*n_e
      tau_ee       = tau_ee_num / tau_ee_denom

      sigma_par_1 = ( SQRT(2.0_r8) + 13.0_r8*(Zeff/4.0_r8) )
     $              / (Zeff*(SQRT(2.0_r8) + Zeff))
      sigma_par_2 = (n_e * (chag**2.0_r8) * tau_ee) / m_e
      sigma_par   = sigma_par_1 * sigma_par_2

c-----------------------------------------------------------------------
c     Characteristic field, Alfvén speed, length scales, and
c     fundamental timescales.
c     Note: mr, nr (rational m/n) and nn (toroidal mode number) are
c     module-level variables set before calling params().
c-----------------------------------------------------------------------
      b_l   = (nr/mr)*rs*sval*bt/R0               ! characteristic B [T]
      v_a   = b_l / (mu0*rho)**0.5_r8                 ! Alfvén velocity [m/s]
      rho_s = 1.02e-4_r8*(mu_i*t_e)**0.5_r8 / bt        ! ion Larmor at T_e [m]
      d_i   = ( (mu_i*m_p) / (n_e*(chag**2)*mu0) )**0.5_r8
                                                    ! ion skin depth [m]

c     The layer timescales and widths below depend on the magnitude of
c     the shear, not its sign: ABS(sval) keeps them positive on a
c     reverse-shear surface.
      tau_h = R0*(mu0*rho)**0.5_r8 / (nn*ABS(sval)*bt) ! Alfvén time [s]
      tau_r = mu0*(rs**2.0_r8)*sigma_par               ! resistive time [s] (Fitzpatrick)
      IF (pr > 0.0_r8) THEN
         tau_v = tau_r / pr                         ! viscous time [s]
         visc  = rho*rs**2.0_r8 / tau_v
      ELSE
         tau_v = 0.0_r8
         visc  = 0.0_r8
      END IF

c     Lundquist number
      lu = tau_r / tau_h

c-----------------------------------------------------------------------
c     Diamagnetic frequencies and normalised Q parameters.
c     Qconv converts dimensional frequencies to the normalised Q
c     used in the SLAYER dispersion relation (Cole scaling).
c-----------------------------------------------------------------------
      Qconv = lu**(1.0_r8/3.0_r8) * tau_h     ! frequency normalisation (Cole)
      tauk  = Qconv                      ! stored in sglobal_mod
c     the layer sees the mode at the Doppler frequency n*omega_E, as in
c     the coupled determinant's shift (slayer.f) and gpec_slayer
      Q     = Qconv * REAL(nn,r8) * omega   ! normalised ExB frequency

c     Diamagnetic frequencies require gradient lengths; skip when the
c     caller does not provide them (l_n <= 0 or l_t <= 0). Callers that
c     compute omega_e/omega_i by other means must set Q_e and Q_i
c     themselves afterwards.
c     omega_* = k_theta*v_* carries the toroidal mode number through
c     m = n*q: the bare gradient-length expression is multiplied by
c     n, as in build_inputs (layerinputs.f).
      IF (l_n > 0.0_r8 .AND. l_t > 0.0_r8) THEN
         omega_e = -REAL(nn,r8)*t_e/(bt*R0)
     $        *(1.0_r8/l_n + 1.0_r8/l_t)*qval
         omega_i =  REAL(nn,r8)*t_i/(bt*R0)
     $        *(1.0_r8/l_n + 1.0_r8/l_t)*qval
         Q_e = -Qconv * omega_e
         Q_i = -Qconv * omega_i
      END IF

c     normalised ion Larmor radius (critical stability parameter)
      ds = lu**(1.0_r8/3.0_r8) * rho_s / rs

c-----------------------------------------------------------------------
c     Plasma beta and Prandtl-number-like transport ratios.
c-----------------------------------------------------------------------
      lbeta  = (5.0_r8/3.0_r8)*mu0*n_e*chag*(t_e+t_i) / bt**2.0_r8
      c_beta = (lbeta / (1.0_r8+lbeta))**0.5_r8

c      kappa-based P_perp (JKP's definition, to implement later if desired)
c      K_val = chis(3) / eta
c      Csq   = c_beta**2.0 + (1.0 - c_beta**2.0)*K_val

c      IF (ABS(Csq) > 0.0) THEN
c        P_perp = Csq
c      ELSE
c        tau_perp = (rs**2.0) / chis(1)
c      END IF

c     effective perpendicular and toroidal Prandtl numbers
      tau_perp = (rs**2.0_r8) / chis(1)
      P_perp   = tau_r / tau_perp           ! perp magnetic Prandtl number

      tau_tor = (rs**2.0_r8) / chis(2)
      P_tor   = tau_r / tau_tor             ! toroidal magnetic Prandtl number

c-----------------------------------------------------------------------
c     Normalised beta-related width and Delta norm factor.
c-----------------------------------------------------------------------
c     d_beta uses Fitzpatrick's definition (tau' form).
c     D = S^(1/3) * iota_e^(1/2) * d_beta/r_s with
c     iota_e = Q_e/(Q_e - Q_i), the electron share of the diamagnetic
c     frequency (see slayer_d_norm). When this call does not set Q_e
c     and Q_i, the equal-gradient-ratio value iota_e = T_e/(T_e + T_i)
c     = 1/(1 + tau) is used provisionally; callers that set Q_e and Q_i
c     afterwards must reset D_norm = slayer_d_norm(D_fac,Q_e,Q_i).
      d_beta = c_beta * d_i
      D_fac  = (d_beta/rs) * lu**(1.0_r8/3.0_r8)
      IF (l_n > 0.0_r8 .AND. l_t > 0.0_r8) THEN
         D_norm = slayer_d_norm(D_fac,Q_e,Q_i)
      ELSE
         D_norm = D_fac * SQRT(1.0_r8/(1.0_r8+tau))
      END IF

      delta_n = lu**(1.0_r8/3.0_r8) / rs  ! normalisation for Delta values

c-----------------------------------------------------------------------
c     Critical Deltaprime (dc_tmp) via iterative chi_parallel calculation.
c     The island-width Wd is iterated to converge the
c     short-mfp / long-mfp interpolation for chi_parallel.
c     dc_type (from sglobal_mod) selects the formula:
c       'rfitzp'   -- R. Fitzpatrick
c       default    -- dc_tmp = 0
c-----------------------------------------------------------------------
      IF (ABS(dr_val) > 0.0_r8) THEN

          IF (.NOT. (chis(1) > 0.0_r8)) THEN
             WRITE(*,*) 'ERROR: the critical Delta-prime needs'
             WRITE(*,*) '  chi_perp > 0, got chi_perp =',chis(1),
     $            ' [m^2/s] at q =',qval,'.'
             WRITE(*,*) '  Set chi_p_prof in slayer.in.'
             STOP 'params: chi_perp must be positive'
          END IF
          vte = SQRT((2.0_r8*(t_e*chag)) / m_e)
          chi_par_smfp = (1.581_r8*tau_ee*(vte**2.0_r8))
     $                   / (1.0_r8 + 0.2535_r8*Zeff)

          Wd = 0.1_r8              ! initial guess
          DO wit = 1, 100
              chi_par_lmfp = (2.0_r8*R0*vte)
     $                       / (SQRT(pi)*nr*ABS(sval)*Wd)
              chi_par = (chi_par_smfp*chi_par_lmfp)
     $                  / (chi_par_smfp + chi_par_lmfp)
              Wd_new = SQRT(8.0_r8)*((chis(1)/chi_par)**0.25_r8)
     $           * (1.0_r8/SQRT((rs/R0)*ABS(sval)*nr))
              IF (ABS(Wd_new - Wd) / MAX(ABS(Wd), 1.0e-30_r8)
     $            < 1.0e-10_r8) THEN
                  Wd = Wd_new
                  EXIT
              END IF
              Wd = Wd_new
          END DO
          IF (wit > 100) THEN
              WRITE(*,*) 'params: Wd iteration failed to converge'
              STOP 'params: Wd iteration did not converge'
          END IF

          SELECT CASE(dc_type)
              CASE('rfitzp')
                  dc_tmp = -(SQRT(2.0_r8)*(pi**(1.5_r8))*dr_val) / Wd
              CASE default
                  dc_tmp = 0.0_r8
          END SELECT

      ELSE
          dc_tmp = 0.0_r8
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
c-----------------------------------------------------------------------
c     subprogram 2. slayer_d_norm.
c     Normalised ion-sound-radius parameter of the Fitzpatrick layer
c     equations, D = S^(1/3) * iota_e^(1/2) * d_beta/r_s, from
c     D_fac = S^(1/3)*d_beta/r_s and iota_e = Q_e/(Q_e - Q_i).
c     iota_e must be positive, which holds when the electron and ion
c     diamagnetic frequencies have opposite signs.
c-----------------------------------------------------------------------
      FUNCTION slayer_d_norm(D_fac_in,Q_e_in,Q_i_in)

      REAL(r8), INTENT(IN) :: D_fac_in  ! S^(1/3)*d_beta/r_s
      REAL(r8), INTENT(IN) :: Q_e_in    ! normalised electron diamag.
      REAL(r8), INTENT(IN) :: Q_i_in    ! normalised ion diamag.
      REAL(r8) :: slayer_d_norm
      REAL(r8) :: iota                  ! Q_e/(Q_e - Q_i)

      IF (Q_e_in == Q_i_in) THEN
         WRITE(*,*) 'ERROR: Q_e = Q_i =',Q_e_in,
     $        ', so iota_e = Q_e/(Q_e - Q_i) is undefined.'
         STOP 'slayer_d_norm: Q_e equals Q_i'
      END IF
      iota = Q_e_in / (Q_e_in - Q_i_in)
      IF (.NOT. (iota > 0.0_r8)) THEN
         WRITE(*,*) 'ERROR: iota_e = Q_e/(Q_e - Q_i) =',iota,
     $        ' is not positive (Q_e =',Q_e_in,', Q_i =',Q_i_in,').'
         WRITE(*,*) '  D_norm = S^(1/3)*sqrt(iota_e)*d_beta/r_s is'
         WRITE(*,*) '  undefined. The electron and ion diamagnetic'
         WRITE(*,*) '  frequencies are expected to have opposite signs.'
         STOP 'slayer_d_norm: iota_e is not positive'
      END IF
      slayer_d_norm = D_fac_in * SQRT(iota)

      RETURN
      END FUNCTION slayer_d_norm

      END MODULE params_mod