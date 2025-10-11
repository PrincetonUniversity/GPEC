c-----------------------------------------------------------------------
c     file bernstein.f.
c     bernstein equilibrium calculations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. bernstein_mod.
c     1. bernstein_metric.
c     2. bernstein_shear.
c     3. bernstein_curvature.
c     4. bernstein_calculation.
c     5. bernstein_xi.
c     6. bernstein_k.
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
      USE free_mod
      USE dcon_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: w,v 

      TYPE(bicube_type):: shear, curvature, B_squared, temp
      TYPE(bicube_type):: g
      TYPE(bicube_type):: bernstein_k, sigma
      LOGICAL :: shear_flag=.TRUE.
      LOGICAL :: curvature_flag=.TRUE.
      REAL(r8) :: phi_level=0
      COMPLEX(r8), ALLOCATABLE :: xi_psi(:,:), b_psi(:,:)
      COMPLEX(r8), ALLOCATABLE :: xi_norm(:,:), b_norm(:,:)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. bernstein_metric.
c     compute the metric (covariant componentsv_ij and contravariant 
c     components w_ij).
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_metric

      INTEGER :: ipsi, itheta, i, j

      REAL(r8) :: theta, psi
      REAL(r8) :: rfac,eta,r,jacfac

      ALLOCATE(w(0:rzphi%mx, 0:rzphi%my, 3, 3))
      ALLOCATE(v(0:rzphi%mx, 0:rzphi%my, 3, 3))

      CALL bicube_alloc(g, rzphi%mx, rzphi%my, 6)
      g%xs = rzphi%xs
      g%ys = rzphi%ys
      g%name = "metric"
      g%xtitle = "psi"
      g%ytitle = "theta"
      g%title = (/" g_12 "," g_22 "," g_33 "," g_23 "," g_31 ","g_12"/)

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

      DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
c-----------------------------------------------------------------------
            g%fs(ipsi,itheta,1) = (v(ipsi,itheta,1,1)**2+
     $                            v(ipsi,itheta,1,2)**2+
     $                            v(ipsi,itheta,1,3)**2)*jacfac*jacfac
            g%fs(ipsi,itheta,2) = (v(ipsi,itheta,2,1)**2+
     $                            v(ipsi,itheta,2,2)**2+
     $                            v(ipsi,itheta,2,3)**2)*jacfac*jacfac
            g%fs(ipsi,itheta,3) = v(ipsi,itheta,3,1)**2*+
     $                            v(ipsi,itheta,3,2)**2*+
     $                            v(ipsi,itheta,3,3)**2*jacfac*jacfac
            g%fs(ipsi,itheta,4) = (v(ipsi,itheta,2,1)*v(ipsi,itheta,3,1)
     $                           +v(ipsi,itheta,2,2)*v(ipsi,itheta,3,2)
     $                           +v(ipsi,itheta,2,3)*v(ipsi,itheta,3,3))
     $                           *jacfac*jacfac 
            g%fs(ipsi,itheta,5) = (v(ipsi,itheta,3,1)*v(ipsi,itheta,1,1)
     $                           +v(ipsi,itheta,3,2)*v(ipsi,itheta,1,2)
     $                           +v(ipsi,itheta,3,3)*v(ipsi,itheta,1,3))
     $                           *jacfac*jacfac 
            g%fs(ipsi,itheta,6) = (v(ipsi,itheta,1,1)*v(ipsi,itheta,2,1)
     $                           +v(ipsi,itheta,1,2)*v(ipsi,itheta,2,2)
     $                           +v(ipsi,itheta,1,3)*v(ipsi,itheta,2,3))
     $                           *jacfac*jacfac 
c-----------------------------------------------------------------------
         END DO
      END DO

      CALL bicube_fit(g, "extrap", "periodic")
      PRINT *, ' > bernstein_compute: metric calculation finshed'
      
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bernstein_metric

c-----------------------------------------------------------------------
c     subprogram 2. bernstein_shear.
c     computes local magnetic shear S(ψ, θ) from Eq. (29) of GPEC Notes.
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_shear

      INTEGER :: ipsi, itheta, shear_unit
      REAL(r8) :: psi, theta, rfac, eta, r, z
      REAL(r8) :: chi1, q, q1, jacfac
      REAL(r8) :: g_psi_g_psi, g_psi_g_theta, g_psi_g_zeta
      REAL(r8) :: dterm_dtheta
      REAL(r8) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
      REAL(r8) :: w11, w12, w13, w21, w22, w23, w31, w32, w33

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

      chi1 = psio * twopi

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
         q1 = sq%f1(4)
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
     $                               * (q1 + dterm_dtheta)
            shear%fs(ipsi, itheta, 2) = r
            shear%fs(ipsi, itheta, 3) = z
         END DO
      END DO
      CALL bicube_fit(shear, "extrap", "periodic")

      shear_unit=58
      CALL ascii_open(shear_unit, "shear.out", "UNKNOWN")
      CALL bicube_write_xy(shear, .TRUE., .FALSE.,
     $                            shear_unit, 0, .FALSE.) 
      CALL ascii_close(shear_unit)

      CALL bicube_dealloc(temp)

      PRINT *, ' > shear calcuation finished'

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bernstein_shear
c-----------------------------------------------------------------------
c     subprogram 3. bernstein_curvature
c     CAUTION - this does not calculate curvature it only cal kappa dot (delpsi)
c     κ·∇ψ = (|∇ψ|²/B²)[μ₀p' + (1/2)(∂B²/∂ψ) + (1/2)(∂B²/∂θ)(∇ψ·∇ψ)/(∇ψ·∇θ)]
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_curvature

      INTEGER :: ipsi, itheta, curvature_unit
      REAL(r8) :: psi, theta
      REAL(r8) :: q,f,chi1, rfac, eta, r, jacfac, p1
      REAL(r8) :: kappa_psi, delpsi, g_psi_g_theta

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
      chi1 = psio * twopi

      DO ipsi = 0, mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         q= sq%f(4)
         DO itheta = 0, mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(twopi*eta) ! major radius
            
            delpsi = SQRT(w(ipsi,itheta,1,1)**2 +
     $                      w(ipsi,itheta,1,2)**2 +
     $                      w(ipsi,itheta,1,3)**2)
            B_squared%fs(ipsi,itheta,1) = (((chi1*delpsi)**2+
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

            delpsi = SQRT(w(ipsi,itheta,1,1)**2 +
     $                      w(ipsi,itheta,1,2)**2 +
     $                      w(ipsi,itheta,1,3)**2)

            g_psi_g_theta = w(ipsi,itheta,1,1)*w(ipsi,itheta,2,1)
     $                    +w(ipsi,itheta,1,2)*w(ipsi,itheta,2,2)
            
            CALL spline_eval(sq, psi, 1)
            p1 = sq%f1(2)

            CALL bicube_eval(B_squared, psi, theta, 1)
            kappa_psi = (delpsi*twopi*chi1)**2 * 
     $          (p1 + 0.5_r8 * B_squared%fx(1) +
     $             0.5_r8 * B_squared%fy(1)*  
     $             g_psi_g_theta / (delpsi**2))
     $            / B_squared%f(1)
            
            curvature%fs(ipsi, itheta, 1) = kappa_psi
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
      END SUBROUTINE bernstein_curvature

c-----------------------------------------------------------------------
c     subprogram 3. bernstein_calculation.
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_calculation

      INTEGER :: ipsi, itheta, k_unit
      REAL(r8) :: psi, theta, curvature_psi

      REAL(r8) :: jacfac, chi1, q, f1, p1
      REAL(r8) :: Bth, Bze, jth, jze, bsq, delpsi
      REAL(r8) :: g22, g23, g33, rfac, eta, r, z

      chi1 = psio * twopi

      CALL metric_calculation
      CALL shear_calculation
      CALL curvature_calculation

      CALL bicube_alloc(sigma, rzphi%mx, rzphi%my, 1)
      sigma%xs = rzphi%xs
      sigma%ys = rzphi%ys
      sigma%name = "sigma"
      sigma%xtitle = "psi"
      sigma%ytitle = "theta"

      CALL bicube_alloc(bernstein_k, rzphi%mx, rzphi%my, 6)
      bernstein_k%name="bernstein_k"
      bernstein_k%xs=rzphi%xs
      bernstein_k%ys=rzphi%ys
      bernstein_k%xtitle = "psi"
      bernstein_k%ytitle = "theta"
      ! term1 = sigma * Shear * (delChi)^2
      ! term2 = j dot B
      ! term3 = kappa dot ∇ψ * 2 * p'
      bernstein_k%title = (/" K ","term1","term2","term3","r","z"/)

      PRINT *, ' > K value calculation started'

      DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
            psi = rzphi%xs(ipsi)
            theta = rzphi%ys(itheta)

            CALL spline_eval(sq, psi, 1)
            f1   = sq%f1(1) / twopi
            p1   = sq%f1(2)
            q    = sq%f(4)

            CAll bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),0)
            jacfac = rzphi%f(4)

            jth = -f1*twopi/jacfac
            jze = q*jth - p1/chi1

            Bth  = chi1/ jacfac
            Bze  = q * chi1 / jacfac

            CALL bicube_eval(g,rzphi%xs(ipsi),rzphi%ys(itheta),0)
            g22 = g%f(2)
            g33 = g%f(3)
            g23 = g%f(4)

            CALL bicube_eval(B_squared,psi,theta,0)
            bsq = B_squared%f(1)
            sigma%fs(ipsi,itheta,1)=(g22*Bth*jth + g33*Bze*jze 
     $                             + g23*(Bth*jze + Bze*jth))/bsq

         END DO
      END DO

      PRINT *, ' > K value calculation started2'

      CALL bicube_fit(sigma, "extrap", "periodic")

       DO ipsi = 0, mpsi
         DO itheta = 0, mtheta
            psi = rzphi%xs(ipsi)
            theta = rzphi%ys(itheta)

            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jacfac = rzphi%f(4)
            rfac=SQRT(rzphi%f(1)) ! minor radius
            eta=(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(twopi*eta) ! major radius
            z=zo+rfac*SIN(twopi*eta)

            delpsi = SQRT(w(ipsi,itheta,1,1)**2 +
     $                      w(ipsi,itheta,1,2)**2 +
     $                      w(ipsi,itheta,1,3)**2)

            bernstein_k%fs(ipsi,itheta,2)= delpsi**2 * chi1**2 *
     $           sigma%fs(ipsi,itheta,1) * shear%fs(ipsi,itheta,1)

            CALL bicube_eval(B_squared, psi,theta,0)
            bernstein_k%fs(ipsi,itheta,3)= B_squared%f(1) *
     $           sigma%fs(ipsi,itheta,1)**2   

            CALL bicube_eval(curvature, psi, theta, 0)
            CALL spline_eval(sq, psi, 1)
            p1   = sq%f1(2)
            curvature_psi = curvature%f(1)
            bernstein_k%fs(ipsi,itheta,4)= curvature_psi * 2.0_r8 * p1

            bernstein_k%fs(ipsi,itheta,1)=
     $            bernstein_k%fs(ipsi,itheta,2) + 
     $            bernstein_k%fs(ipsi,itheta,3) +
     $            bernstein_k%fs(ipsi,itheta,4)

            bernstein_k%fs(ipsi, itheta, 5) = r
            bernstein_k%fs(ipsi, itheta, 6) = z

         END DO
       END DO

       CALL bicube_fit(bernstein_k, "extrap", "periodic")

      PRINT *, ' > K value calculation started3'
c-----------------------------------------------------------------------
c     file print
c-----------------------------------------------------------------------
      k_unit = 60
      CALL ascii_open(k_unit, "bernstein_k.out", "UNKNOWN")
      CALL bicube_write_xy(bernstein_k, .TRUE., .FALSE.,
     $                     k_unit, 0, .FALSE.)
      CALL ascii_close(k_unit)

      PRINT *, ' > K calculation finished'
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bernstein_calculation

c-----------------------------------------------------------------------
c     subprogram 4.
c     for displacement calculation - 100% same with match
c-----------------------------------------------------------------------
       SUBROUTINE bernstein_xi
       COMPLEX(r8) :: xipsi,bpsi
       REAL(r8):: r, z, psi, theta
       
       INTEGER :: sol_num=1
       INTEGER :: istep,ifix,jfix,kfix,ieq,info
       INTEGER :: ipert
       REAL(r8) :: r2,deta,dphi,rfac,eta,jac,v21,v22,dpsisq,norm,singfac
       COMPLEX(r8) :: expfac,expfac1, mstep

       INTEGER :: ipsi,itheta,nqty=4

       COMPLEX(r8), DIMENSION(mpert):: xi_vec
       INTEGER, DIMENSION(mpert) :: ipiv
       COMPLEX(r8), DIMENSION(mpert) :: uedge,temp1
       COMPLEX(r8), DIMENSION(mpert,mpert) :: temp2

       ALLOCATE(xi_psi(0:rzphi%mx, 0:rzphi%my))
       ALLOCATE(xi_norm(0:rzphi%mx, 0:rzphi%my))
       ALLOCATE(b_psi(0:rzphi%mx, 0:rzphi%my))
       ALLOCATE(b_norm(0:rzphi%mx, 0:rzphi%my))
c-----------------------------------------------------------------------
c     construct uedge.
c-----------------------------------------------------------------------
c       uedge=wt(:,sol_num)
c       temp2=soltype(mstep)%u(:,1:mpert,1)
c       CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
c       CALL zgetrs('N',mpert,1,temp2,mpert,ipiv,uedge,mpert,info)
c-----------------------------------------------------------------------
c     construct eigenfunctions.
c-----------------------------------------------------------------------
c        jfix=0
c        DO ifix=0,mfix
c           temp1=MATMUL(fixtype(ifix)%transform,uedge)
c           kfix=fixstep(ifix+1)
c           DO ieq=1,4
c              DO istep=jfix,kfix
c                 xi_vec(:,ieq,istep)
c      $              =MATMUL(soltype(istep)%u(:,1:mpert,ieq),temp1)
c              ENDDO
c           ENDDO
c           jfix=kfix+1
c        ENDDO
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
c       DO ifix=0,mfix
c          DEALLOCATE(fixtype(ifix)%transform)
c       ENDDO
c-----------------------------------------------------------------------
c     compute contour values.
c-----------------------------------------------------------------------
c       WRITE(*,*)"Compute values for contour plots"
c       DO ipsi=0,mpsi
c         DO itheta=0,mtheta
c              psi = rzphi%xs(ipsi)
c              theta = rzphi%ys(itheta)
c-----------------------------------------------------------------------
c     evaluate radial position.
c-----------------------------------------------------------------------
c              CALL bicube_eval(rzphi,psi,theta,1)
c              r2=rzphi%f(1)
c              deta=rzphi%f(2)
c              dphi=rzphi%f(3)
c              rfac=SQRT(r2)
c              eta=theta+deta
c              r=ro+rfac*COS(twopi*eta)
c              z=zo+rfac*SIN(twopi*eta)
c-----------------------------------------------------------------------
c     evaluate normalizing factors.
c-----------------------------------------------------------------------
c              jac=rzphi%f(4)
c              v21=v(ipsi,itheta, 2,1)
c              v22=v(ipsi, itheta, 2, 2)
c              dpsisq=(twopi*r)**2*(v21**2+v22**2)
c              CALL spline_eval(sq,psi,0)
c              expfac=EXP(ifac*(mlow*twopi*theta-nn*(phi_level-dphi)))
c              expfac1=EXP(ifac*twopi*theta)
c              singfac=mlow-nn*sq%f(4)
c              xipsi=0
c              bpsi=0
c              DO ipert=1,mpert
c                     xipsi=xipsi+xi_vec(ipert)*expfac
c                     bpsi=bpsi+xi_vec(ipert)*expfac*singfac*ifac
c                     singfac=singfac+1
c                     expfac=expfac*expfac1
c              ENDDO
c              xi_psi(ipsi,itheta) = xipsi
c-----------------------------------------------------------------------
c     compute xi_norm and b_norm.
c-----------------------------------------------------------------------
c              norm=1/SQRT(dpsisq)
c              xi_norm(ipsi,itheta)=xipsi*norm
c              b_psi(ipsi,itheta)=bpsi/jac
c              b_norm(ipsi,itheta)=bpsi*norm
c         ENDDO
c       ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
       RETURN
       END SUBROUTINE bernstein_xi
c-----------------------------------------------------------------------
c     subprogram 4. K integral
c-----------------------------------------------------------------------
      SUBROUTINE bernstein_k

       INTEGER:: mstep, mthet

       INTEGER :: ipsi, itheta
       REAL(r8) :: psi, theta
       REAL(r8) :: jacfac, k, delpsi, k_1, k_2, k_3
       COMPLEX(r8) :: xipsi
       TYPE(bicube_type):: temp

       CALL bicube_alloc(temp, rzphi%mx, rzphi%my, 3)
       temp%xs = rzphi%xs
       temp%ys = rzphi%ys
       temp%name = "temp"
       temp%xtitle = "psi"
       temp%ytitle = "theta"
       temp%title = (/"temp"," temp1 "," temp2 ","  temp3"/)

       DO ipsi = 0, mpsi
          DO itheta = 0, mtheta
              psi = rzphi%xs(ipsi)
              theta = rzphi%ys(itheta)

              delpsi = SQRT(w(ipsi,itheta,1,1)**2 +
     $                      w(ipsi,itheta,1,2)**2 +
     $                      w(ipsi,itheta,1,3)**2)
              CALL bicube_eval(rzphi, psi, theta, 1)
              CALL bicube_eval(bernstein_k, psi, theta, 1)
              jacfac = rzphi%f(4)
              k = bernstein_k%f(1)
              k_1 = bernstein_k%f(2)
              k_2 = bernstein_k%f(3)
              k_3 = bernstein_k%f(4)

              xipsi=xi_psi(ipsi, itheta)

              temp%fs(ipsi, itheta, 1) = jacfac*k*ABS(xipsi)**2
     $                     /(delpsi**2)

              temp%fs(ipsi, itheta, 2) = jacfac*k_1*ABS(xipsi)**2
     $                     /(delpsi**2)
              temp%fs(ipsi, itheta, 3) = jacfac*k_2*ABS(xipsi)**2
     $                     /(delpsi**2)
              temp%fs(ipsi, itheta, 4) = jacfac*k_3*ABS(xipsi)**2
     $                     /(delpsi**2)

          END DO
       END DO

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bernstein_k

c-----------------------------------------------------------------------
c     module end
c-----------------------------------------------------------------------
      END MODULE bernstein_mod