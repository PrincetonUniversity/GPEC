c-----------------------------------------------------------------------
c     file bernstein.f.
c     bernstein equilibrium calculations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. bernstein_mod.
c     1. metric_calculation.
c     2. shear_calculation.
c     3. curvature_calculation.
c     4. bernstein_calculation.
c     7. bernstein_out.
c-----------------------------------------------------------------------
c     subprogram 0. bernstein_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE bernstein_mod
      USE global_mod
      USE bicube_mod
      USE spline_mod
      USE direct_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: w,v
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: jacobian

      TYPE(bicube_type):: shear, curvature, B_magnitude

      LOGICAL :: shear_flag=.TRUE.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. metric_calculation.
c     compute the metric (covariant componentsv_ij and contravariant 
c     components w_ij).
c-----------------------------------------------------------------------
      SUBROUTINE metric_calculation

      REAL(r8), PARAMETER :: pi = 3.141592653589793_r8
      REAL(r8), PARAMETER :: twopi = 2.0_r8 * pi
      REAL(r8), PARAMETER :: r_eps = 1e-10_r8  ! r = 0 
      INTEGER :: ipsi, itheta, i, j

      REAL(r8) :: theta, psi
      REAL(r8) :: r_minor_sq, eta_periodic_part, nu_func
      REAL(r8) :: dr2_dpsi, dr2_dtheta,
     $              deta_p_dpsi, deta_p_dtheta,
     $              dnu_dpsi, dnu_dtheta, jacobian_ij
      REAL(r8) :: eta_norm, deta_dpsi, deta_dtheta
      REAL(r8) :: r_minor, R_major, eta_rad, cos_eta, sin_eta
      REAL(r8) :: inv_J, inv_r, inv_2rJ, inv_rJ, inv_2pi

      REAL(r8), PARAMETER :: small_r = 1.0E-9_r8

      ALLOCATE(w(0:rzphi%mx, 0:rzphi%my, 3, 3))
      ALLOCATE(v(0:rzphi%mx, 0:rzphi%my, 3, 3))
      ALLOCATE(jacobian(0:rzphi%mx, 0:rzphi%my))

      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)
            jacobian(ipsi, itheta) = 0.0
            DO i = 1, 3
               w(ipsi, itheta, i, 1) = 0.0
               w(ipsi, itheta, i, 2) = 0.0
               w(ipsi, itheta, i, 3) = 0.0

               v(ipsi, itheta, i, 1) = 0.0
               v(ipsi, itheta, i, 2) = 0.0
               v(ipsi, itheta, i, 3) = 0.0

            END DO
         END DO
      END DO

      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)
            
            CALL bicube_eval(rzphi, psi, theta, 1)

            r_minor_sq        = rzphi%f(1)
            eta_periodic_part = rzphi%f(2)
            nu_func           = rzphi%f(3)
            jacobian_ij = rzphi%f(4)
            dr2_dpsi    = rzphi%fx(1)
            deta_p_dpsi = rzphi%fx(2)
            dnu_dpsi    = rzphi%fx(3)
            dr2_dtheta    = rzphi%fy(1)
            deta_p_dtheta = rzphi%fy(2)
            dnu_dtheta    = rzphi%fy(3)

            eta_norm    = eta_periodic_part + theta
            deta_dpsi   = deta_p_dpsi
            deta_dtheta = deta_p_dtheta + 1.0_r8

            r_minor = SQRT(MAX(0.0_r8, r_minor_sq))
            eta_rad = eta_norm * twopi
            cos_eta = COS(eta_rad)
            sin_eta = SIN(eta_rad)
            R_major = ro + r_minor * cos_eta

            inv_J = 1.0_r8 / jacobian_ij
            inv_2pi = 1.0_r8 / twopi
            IF (r_minor > small_r) THEN
                  inv_r = 1.0_r8 / r_minor
            ELSE
                  inv_r = 0.0_r8
            ENDIF
            inv_2rJ = 0.5_r8 * inv_r * inv_J
            inv_rJ  = inv_r * inv_J

            jacobian(ipsi, itheta) = jacobian_ij
c-----------------------------------------------------------------------
            v(ipsi, itheta, 1, 1) = inv_2rJ * dr2_dpsi
            v(ipsi, itheta, 1, 2) = inv_J * (r_minor*twopi) * deta_dpsi
            v(ipsi, itheta, 1, 3) = inv_J * R_major * dnu_dpsi
            v(ipsi, itheta, 2, 1) = inv_2rJ * dr2_dtheta
            v(ipsi, itheta, 2, 2) = inv_J 
     $                          * (r_minor*twopi) * deta_dtheta
            v(ipsi, itheta, 2, 3) = inv_J * R_major * dnu_dtheta
            v(ipsi, itheta, 3, 3) = inv_J * R_major * twopi

            w(ipsi, itheta, 1, 1) = (twopi**2*r_minor*R_major*inv_J)
     $                      * deta_dtheta
            w(ipsi, itheta, 1, 2) = (-pi * R_major * inv_rJ) 
     $                           * dr2_dtheta
            w(ipsi, itheta, 2, 1) = (-twopi**2*r_minor*R_major*inv_J)
     $                      * deta_dpsi
            w(ipsi, itheta, 2, 2) = (pi * R_major * inv_rJ) * dr2_dpsi
            w(ipsi, itheta, 3, 1) = inv_2pi 
     $                        * (-dnu_dpsi*w(ipsi, itheta, 1, 1)
     $                      - dnu_dtheta*w(ipsi, itheta, 2, 1))
            w(ipsi, itheta, 3, 2) = inv_2pi 
     $                        * (-dnu_dpsi*w(ipsi, itheta, 1, 2)
     $                      - dnu_dtheta*w(ipsi, itheta, 2, 2))
            w(ipsi, itheta, 3, 3) = inv_2pi * ((1.0_r8 / R_major)
     $                      - dnu_dpsi*w(ipsi, itheta, 1, 3)
     $                      - dnu_dtheta*w(ipsi, itheta, 2, 3))
         END DO
      END DO

      PRINT *, ' > bernstein_compute: metric calculation finshed'
      PRINT *, ' > v(0,0,1,1) = ', v(0, 0, 1, 1)
      PRINT *, ' > w(0,0,1,1) = ', w(0, 0, 1, 1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE metric_calculation

c-----------------------------------------------------------------------
c     subprogram 2. shear_calculation.
c     computes local magnetic shear S(ψ, θ) from Eq. (29) of GPEC Notes.
c-----------------------------------------------------------------------
      SUBROUTINE shear_calculation

      INTEGER :: ipsi, itheta, shear_unit
      REAL(r8) :: psi, theta, g_psi_g_theta, g_psi_g_zeta, g_psi_g_psi
      REAL(r8) :: r2, r, eta, R_val, Z, q_val,f_psi, r_minor, eta_norm
      REAL(r8) :: R_major_local, chi_prime, q_prime, dTerm_dtheta
      INTEGER, PARAMETER :: out_unit = 98
      REAL(r8), PARAMETER :: pi = 3.141592653589793_r8
      REAL(r8), PARAMETER :: twopi = 2.0_r8 * pi

      TYPE(bicube_type) :: temp 
      !this is for q(gpsipsi)-(gpsitheta)/(gpsizeta)

      CALL bicube_alloc(shear, rzphi%mx, rzphi%my, 1)

      shear%xs = rzphi%xs
      shear%ys = rzphi%ys
      shear%name = "shear"
      shear%xtitle = "psi"
      shear%ytitle = "theta"

      CALL bicube_alloc(temp, rzphi%mx, rzphi%my, 1)

      temp%name="temporary_term"
      temp%xs=rzphi%xs
      temp%ys=rzphi%ys
      temp%xtitle = "psi"
      temp%ytitle = "theta"
      PRINT *, ' > shear calcuation'
      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)

            g_psi_g_psi = w(ipsi,itheta,1,1)**2
     $                    +w(ipsi,itheta,1,2)**2
     $                    +w(ipsi,itheta,1,3)**2   

            g_psi_g_theta = w(ipsi,itheta,1,1)*w(ipsi,itheta,2,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,2,2)

            g_psi_g_zeta = w(ipsi,itheta,1,1)*w(ipsi,itheta,3,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,3,2)

            CALL spline_eval(sq, psi, 0)
            q_val = sq%f(4)
            f_psi = sq%f(1) / twopi

            IF (g_psi_g_psi > 1.0E-12_r8) THEN
              temp%fs(ipsi, itheta, 1) = 
     $             (q_val * g_psi_g_theta - g_psi_g_zeta) / g_psi_g_psi
            ELSE
              temp%fs(ipsi, itheta, 1) = 0.0_r8
            ENDIF

            r_minor = SQRT(MAX(0.0_r8, rzphi%fs(ipsi,itheta,1)))
            eta_norm = rzphi%fs(ipsi,itheta,2) + theta
            R_major_local = ro + r_minor * COS(eta_norm * twopi)
         
         END DO
      END DO
      PRINT *, ' > shear calcuation'
      CALL bicube_fit(temp, "extrap", "periodic")

      chi_prime = psio * twopi
      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)
            CALL spline_eval(sq, psi, 1)
            q_prime = sq%f1(4)

            CALL bicube_eval(temp, psi, theta, 1)
            dTerm_dtheta = temp%fy(1)

            shear%fs(ipsi, itheta, 1) = (chi_prime**2 
     $                                 / jacobian(ipsi,itheta)) 
     $                               * (q_prime + dTerm_dtheta)
         END DO
      END DO
      PRINT *, ' > shear calcuation'
      CALL bicube_fit(shear, "extrap", "periodic")

      shear_unit = 101
      IF (shear_flag) THEN
         CALL ascii_open(shear_unit, "shear.out", "UNKNOWN")
         CALL bicube_write_xy(shear, .TRUE., .FALSE.
     $       , shear_unit, 0, .FALSE.)
         CALL ascii_close(shear_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE shear_calculation
c-----------------------------------------------------------------------
c     subprogram 3. curvature_calculation
c     CAUTION - this does not calculate curvature it only cal kappa dot (delpsi)
c-----------------------------------------------------------------------

      SUBROUTINE curvature_calculation

      INTEGER :: ipsi, itheta, curvature_unit
      REAL(r8) :: psi, theta
      REAL(r8) :: kappa_dot_grad_psi
      REAL(r8), PARAMETER :: pi = 3.141592653589793_r8
      REAL(r8), PARAMETER :: twopi = 2.0_r8 * pi

      CALL bicube_alloc(curvature, rzphi%mx, rzphi%my, 1)

      curvature%xs = rzphi%xs
      curvature%ys = rzphi%ys
      curvature%name = "curvature"
      curvature%xtitle = "psi"
      curvature%ytitle = "theta"

      PRINT *, ' > curvature calculation'
      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)
            
            ! 여기서 κ·∇ψ 계산
            ! 구체적인 공식이 필요해요!
            
            curvature%fs(ipsi, itheta, 1) = kappa_dot_grad_psi
         END DO
      END DO

      CALL bicube_fit(curvature, "extrap", "periodic")
      PRINT *, ' > curvature calculation finished'

      curvature_unit = 102
      CALL ascii_open(curvature_unit, "curvature.out", "UNKNOWN")
      CALL bicube_write_xy(PgradB_spline, .TRUE., .FALSE.
   $       , curvature_unit, 0, .FALSE.)
      CALL ascii_close(curvature_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curvature_calculation

      END MODULE bernstein_mod