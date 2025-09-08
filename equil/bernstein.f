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

      TYPE(bicube_type):: shear, curvature, B_squared
      TYPE(bicube_type):: bernstein_k, sigma
      LOGICAL :: shear_flag=.TRUE.
      LOGICAL :: curvature_flag=.TRUE.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. metric_calculation.
c     compute the metric (covariant componentsv_ij and contravariant 
c     components w_ij).
c-----------------------------------------------------------------------
      SUBROUTINE metric_calculation

      INTEGER :: ipsi, itheta, i, j

      REAL(r8) :: theta, psi
      REAL(r8) :: rfac,eta,r,jacfac

      ALLOCATE(w(0:rzphi%mx, 0:rzphi%my, 3, 3))
      ALLOCATE(v(0:rzphi%mx, 0:rzphi%my, 3, 3))

      DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
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

      DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(twopi*eta) ! major radius
c-----------------------------------------------------------------------
            v(ipsi, itheta, 1, 1) = rzphi%fx(1)/(2*rfac*jacfac)
            v(ipsi, itheta, 1, 2) = rzphi%fx(2)*twopi*rfac/jacfac
            v(ipsi, itheta, 1, 3) = rzphi%fx(3)*r/jacfac

            v(ipsi, itheta, 2, 1) = rzphi%fy(1)/(2*rfac*jacfac)
            v(ipsi, itheta, 2, 2) = (1+rzphi%fy(2))*twopi*rfac/jacfac
            v(ipsi, itheta, 2, 3) = rzphi%fy(3)*r/jacfac

            v(ipsi, itheta, 3, 3) = twopi*r/jacfac
c-----------------------------------------------------------------------
            w(ipsi, itheta, 1, 1) = (1+rzphi%fy(2))*
     $                              twopi**2*rfac*r/jacfac
            w(ipsi, itheta, 1, 2) = -rzphi%fy(1)*pi*r/(rfac*jacfac)

            w(ipsi, itheta, 2, 1) = -twopi**2*rfac*r*rzphi%fx(2)/jacfac
            w(ipsi, itheta, 2, 2) = pi*r*rzphi%fx(1)/(rfac*jacfac)

            w(ipsi, itheta, 3, 1) = (twopi*r*rfac/jacfac)*
     $            (rzphi%fx(2)*rzphi%fy(3)-rzphi%fx(3)*(1+rzphi%fy(2)))
            w(ipsi, itheta, 3, 2) = (r/(2*rfac*jacfac))*
     $             (rzphi%fx(3)*rzphi%fy(1)-rzphi%fx(1)*rzphi%fy(3))
            w(ipsi, itheta, 3, 3) = 1/(twopi*r)
         END DO
      END DO

      PRINT *, ' > bernstein_compute: metric calculation finshed'
      
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
      REAL(r8) :: psi, theta, rfac, eta, r, z
      REAL(r8) :: chi_prime,q, q_prime, jacfac
      REAL(r8) :: g_psi_g_psi, g_psi_g_theta, g_psi_g_zeta
      REAL(r8) :: dterm_dtheta
      REAL(r8) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
      REAL(r8) :: w11, w12, w13, w21, w22, w23, w31, w32, w33

      TYPE(bicube_type) :: temp 

      CALL bicube_alloc(shear, rzphi%mx, rzphi%my, 3)
      shear%xs = rzphi%xs
      shear%ys = rzphi%ys
      shear%name = "shear"
      shear%title = (/"  S  "," r ", "z"/)
      shear%xtitle = "psi"
      shear%ytitle = "theta"

      CALL bicube_alloc(temp, rzphi%mx, rzphi%my, 5)
      temp%name="temporary_term"
      temp%xs=rzphi%xs
      temp%ys=rzphi%ys
      temp%xtitle = "psi"
      temp%ytitle = "theta"
      PRINT *, ' > shear calcuation'

      chi_prime = psio * twopi

      PRINT *, 'psio, chi_prime = ', psio, chi_prime

      DO ipsi = 0, mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         q = sq%f(4)
         DO itheta = 0, mtheta
            psi = rzphi%xs(ipsi)
            theta = rzphi%ys(itheta)

            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(twopi*eta) ! major radius

            v11 = v(ipsi,itheta,1,1)
            v12 = v(ipsi,itheta,1,2)
            v13 = v(ipsi,itheta,1,3)
            v21 = v(ipsi,itheta,2,1)
            v22 = v(ipsi,itheta,2,2)
            v23 = v(ipsi,itheta,2,3)
            v31 = v(ipsi,itheta,3,1)
            v32 = v(ipsi,itheta,3,2)
            v33 = v(ipsi,itheta,3,3)

            g_psi_g_psi = w(ipsi,itheta,1,1)**2
     $                    +w(ipsi,itheta,1,2)**2

            g_psi_g_theta = w(ipsi,itheta,1,1)*w(ipsi,itheta,2,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,2,2)

            g_psi_g_zeta = w(ipsi,itheta,1,1)*w(ipsi,itheta,3,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,3,2)

            temp%fs(ipsi, itheta, 1) = twopi*r*jacfac
     $           * (-(v11*v21+v12*v22)*(v23+q*v33)+v13*(v21**2+v22**2))
            temp%fs(ipsi, itheta, 2) = (q* g_psi_g_theta - g_psi_g_zeta) 
            temp%fs(ipsi, itheta, 3) = g_psi_g_psi
            temp%fs(ipsi, itheta, 4) = twopi**2 * r**2 
     $           * (v21**2 + v22**2)

            temp%fs(ipsi, itheta, 5) = 
     $               temp%fs(ipsi,itheta,1)/temp%fs(ipsi,itheta,4)
         END DO
      END DO

      CALL bicube_fit(temp, "extrap", "periodic")

      DO ipsi = 0, mpsi
         CALL spline_eval(sq,sq%xs(ipsi),1)
         q_prime = sq%f1(4)
         DO itheta = 0, mtheta
            psi = rzphi%xs(ipsi)
            theta = rzphi%ys(itheta)

            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(twopi*eta) ! major radius
            z=zo+rfac*SIN(twopi*eta)

            CALL bicube_eval(temp, psi, theta, 1)
            dterm_dtheta = temp%fy(5)

            shear%fs(ipsi, itheta, 1) = (twopi**2 / jacfac)
     $                               * (q_prime + dterm_dtheta)
            shear%fs(ipsi, itheta, 2) = r
            shear%fs(ipsi, itheta, 3) = z
         END DO
      END DO
      PRINT *, ' > shear calcuation finished'
      CALL bicube_fit(shear, "extrap", "periodic")

      shear_unit=58
      CALL ascii_open(shear_unit, "shear.out", "UNKNOWN")
      CALL bicube_write_xy(shear, .TRUE., .FALSE.,
     $                            shear_unit, 0, .FALSE.) 
      CALL ascii_close(shear_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE shear_calculation
c-----------------------------------------------------------------------
c     subprogram 3. curvature_calculation
c     CAUTION - this does not calculate curvature it only cal kappa dot (delpsi)
c     κ·∇ψ = (|∇ψ|²/B²)[μ₀p' + (1/2)(∂B²/∂ψ) + (1/2)(∂B²/∂θ)(∇ψ·∇ψ)/(∇ψ·∇θ)]
c-----------------------------------------------------------------------
      SUBROUTINE curvature_calculation

      INTEGER :: ipsi, itheta, curvature_unit
      REAL(r8) :: psi, theta
      REAL(r8) :: q,f,chi_prime, rfac, eta, r, jacfac, p_prime
      REAL(r8) :: kappa_dot_grad_psi, grad_psi_squared
      REAL(r8) ::delpsi, B_squared2,grad_psi_dot_grad_theta

      CALL bicube_alloc(curvature, rzphi%mx, rzphi%my, 1)
      curvature%xs = rzphi%xs
      curvature%ys = rzphi%ys
      curvature%name = "curvature"
      curvature%xtitle = "psi"
      curvature%ytitle = "theta"

      CALL bicube_alloc(B_squared, rzphi%mx, rzphi%my, 2)
      B_squared%xs = rzphi%xs
      B_squared%ys = rzphi%ys
      B_squared%name = "B^2"
      B_squared%xtitle = "psi"
      B_squared%ytitle = "theta"

      PRINT *, ' > curvature calculation: computing B²'
      chi_prime = psio * twopi

       DO ipsi = 0, mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         q= sq%f(4)
         DO itheta = 0, mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)

            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=twopi*(itheta/REAL(mtheta,r8)+rzphi%f(2))
            r=ro+rfac*COS(eta) ! major radius
            jacfac=rzphi%f(4) ! Jacobian
            
            delpsi = SQRT(w(ipsi,itheta,1,1)**2 +
     $                      w(ipsi,itheta,1,2)**2 +
     $                      w(ipsi,itheta,1,3)**2)
            B_squared%fs(ipsi,itheta,1) = (((chi_prime*delpsi)**2+
     $           sq%f(1)**2)/(twopi*r)**2)
            B_squared%fs(ipsi,itheta,2) = eqfun%fs(ipsi,itheta,1)**2

         END DO
      END DO

      CALL bicube_fit(B_squared, "extrap", "periodic")

      PRINT *, ' > curvature calculation: computing κ·∇ψ'

       DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
            psi = rzphi%xs(ipsi)
            theta = rzphi%ys(itheta)

            g_psi_g_psi = w(ipsi,itheta,1,1)**2
     $                    +w(ipsi,itheta,1,2)**2

            g_psi_g_theta = w(ipsi,itheta,1,1)*w(ipsi,itheta,2,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,2,2)
            
            CALL spline_eval(sq, psi, 1)
            p_prime = sq%f1(2)

            CALL bicube_eval(B_squared, psi, theta, 1)

            kappa_dot_grad_psi = (g_psi_g_psi* twopi**4) * 
     $          (p_prime + 0.5_r8 * B_squared%fx(1) +
     $             0.5_r8 * B_squared%fy(1)*  
     $             _psi_g_theta / g_psi_g_psi)
     $            / B_squared%f(1)
            
            curvature%fs(ipsi, itheta, 1) = kappa_dot_grad_psi
         END DO
      END DO

c-----------------------------------------------------------------------
c     curvature spline fitting
c-----------------------------------------------------------------------
      CALL bicube_fit(curvature, "extrap", "periodic")

c-----------------------------------------------------------------------
c     file print
c-----------------------------------------------------------------------
      curvature_unit = 59
      IF (curvature_flag) THEN
         CALL ascii_open(curvature_unit, "curvature.out", "UNKNOWN")
         CALL bicube_write_xy(curvature, .TRUE., .FALSE.,
     $                     curvature_unit, 0, .FALSE.)
         CALL ascii_close(curvature_unit)
      ENDIF

      PRINT *, ' > curvature calculation finished'
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curvature_calculation

c-----------------------------------------------------------------------
c     subprogram 3. bernstein_calculation.
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_calculation

      INTEGER :: ipsi, itheta
      REAL(r8) :: psi, theta

      REAL(r8) :: J, chi_p, q, f_p, p_p
      REAL(r8) :: e1(3), e2(3), e3(3)                         ! e_psi, e_theta, e_zeta
      REAL(r8) :: G11, G22, G33, G12, G13, G23                ! metric G_ij = e_i·e_j
      REAL(r8) :: Bth, Bze, jth, jze, B_squared, B_squared2
      REAL(r8), parameter :: eps_chi = 1.0e-12_r8
      REAL(r8), parameter :: eps_B2  = 1.0e-20_r8

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: jdotb

      ALLOCATE(jdotb(0:rzphi%mx, 0:rzphi%my))

      chi_p = psio * twopi

      CALL metric_calculation
      CALL shear_calculation
      CALL curvature_calculation

      CALL bicube_alloc(sigma, rzphi%mx, rzphi%my, 1)

      sigma%xs = rzphi%xs
      sigma%ys = rzphi%ys
      sigma%name = "sigma"
      sigma%xtitle = "psi"
      sigma%ytitle = "theta"

      CALL bicube_alloc(bernstein_k, rzphi%mx, rzphi%my, 1)

      bernstein_k%name="bernstein_k"
      bernstein_k%xs=rzphi%xs
      bernstein_k%ys=rzphi%ys
      bernstein_k%xtitle = "psi"
      bernstein_k%ytitle = "theta"
      bernstein_k%title = (/"  K  "," term1 "," term2 ","  term3 "/)

      DO itheta = 0, rzphi%my
         theta = rzphi%ys(itheta)
         DO ipsi = 0, rzphi%mx
            psi = rzphi%xs(ipsi)

            CALL spline_eval(sq, psi, 0)
            f_p   = sq%f1(1) / twopi
            p_p   = sq%f1(2)
            q    = sq%f(4)

            ! ---- covariant basis e_i components from w ----
            e1(:) = w(ipsi, itheta, 1, 1:3)   ! e_psi
            e2(:) = w(ipsi, itheta, 2, 1:3)   ! e_theta
            e3(:) = w(ipsi, itheta, 3, 1:3)   ! e_zeta

            G11 = e1(1)*e1(1) + e1(2)*e1(2) + e1(3)*e1(3)
            G22 = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
            G33 = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
            G12 = e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)
            G13 = e1(1)*e3(1) + e1(2)*e3(2) + e1(3)*e3(3)
            G23 = e2(1)*e3(1) + e2(2)*e3(2) + e2(3)*e3(3)

            Bth  = chi_p / J
            Bze  = q * chi_p / J

         END DO
      END DO
      

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bernstein_calculation
c-----------------------------------------------------------------------
c     module end
c-----------------------------------------------------------------------
      END MODULE bernstein_mod