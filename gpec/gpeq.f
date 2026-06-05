c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     calculate functions for perturbed equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. gpeq_mod
c      1. gpeq_sol
c      2. gpeq_contra
c      3. gpeq_cova
c      4. gpeq_normal
c      5. gpeq_tangent
c      6. gpeq_parallel
c      7. gpeq_rzphi
c      8. gpeq_surface
c      9. gpeq_epf          (reconstruction: C vector for EPF)
c     10. gpeq_dst          (reconstruction: C vector for DST)
c     11. gpeq_shear        (reconstruction: magnetic shear)
c     12. gpeq_curvature    (reconstruction: curvature)
c     13. gpeq_K            (reconstruction: Bernstein K quantity)
c    13b. gpeq_recon3       (reconstruction: |C|^2 decomposition)
c     14. gpeq_fcoords
c     15. gpeq_fcoordsout
c     16. gpeq_bcoords
c     17. gpeq_bcoordsout
c     18. gpeq_weight
c     19. gpeq_rzpgrid
c     20. gpeq_rzpdiv
c     21. gpeq_alloc
c     22. gpeq_dealloc
c     23. gpeq_interp_singsurf
c     24. gpeq_interp_sol
c-----------------------------------------------------------------------
c     subprogram 0. gpeq_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gpeq_mod
      USE idcon_mod
 
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gpeq_sol.
c     obtain solutions of perturbed quantities.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_sol(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER, DIMENSION(mpert) :: ipiva
      COMPLEX(r8), DIMENSION(mpert) :: xspfac

      COMPLEX(r8), DIMENSION(mpert*mpert) :: work

      IF(debug_flag) PRINT *, "Entering gpeq_sol"
c-----------------------------------------------------------------------
c     evaluate matrices and solutions.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,1)
      q=sq%f(4)
      q1=sq%f1(4)
      singfac=mfac-nn*q

      ALLOCATE(amat(mpert,mpert),bmat(mpert,mpert),cmat(mpert,mpert))
      ALLOCATE(fmats(mband+1,mpert),gmats(mband+1,mpert))
      ALLOCATE(kmats(2*mband+1,mpert))
      CALL idcon_matrix(psi)
c-----------------------------------------------------------------------
c     compute preliminary quantities.
c-----------------------------------------------------------------------
      CALL cspline_eval(u1,psi,0)
      xsp_mn=u1%f
      IF (kin_flag) THEN
         CALL cspline_eval(u3,psi,0)
         CALL cspline_eval(u4,psi,0)
         xsp1_mn=u3%f
         xss_mn=u4%f
      ELSE
         IF (galsol%gal_flag) THEN
            CALL cspline_eval(u1,psi,1)
            xsp1_mn=u1%f1
         ELSE
            CALL cspline_eval(u2,psi,0)
            xspfac=u2%f/singfac
            CALL zgbmv('N',mpert,mpert,mband,mband,-ione,kmats,
     $         2*mband+1,u1%f,1,ione,xspfac,1)
            CALL zpbtrs('L',mpert,mband,1,fmats,mband+1,xspfac,mpert,
     $         info)
            xsp1_mn=xspfac/singfac
         ENDIF
         CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,bmat,mpert,info)
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,cmat,mpert,info)
         xss_mn=-MATMUL(bmat,xsp1_mn)-MATMUL(cmat,xsp_mn)
      ENDIF
c-----------------------------------------------------------------------
c     compute Jacobian-weighted contravariant b fields: J b^i.
c-----------------------------------------------------------------------
      bwp_mn=(chi1*singfac*twopi*ifac*xsp_mn)
      bwt_mn=-(chi1*xsp1_mn+twopi*ifac*nn*xss_mn)
      bwz_mn=-(chi1*(q1*xsp_mn+sq%f(4)*xsp1_mn)+twopi*ifac*mfac*xss_mn)
c-----------------------------------------------------------------------
c     compute derivative of b fields.
c-----------------------------------------------------------------------
      bwp1_mn=(twopi*ifac*chi1*singfac)*xsp1_mn-
     $     twopi*ifac*chi1*nn*q1*xsp_mn
c-----------------------------------------------------------------------
c     compute modified quantities.
c-----------------------------------------------------------------------
      IF (reg_flag) THEN
         xmp1_mn=xsp1_mn*(singfac**2/(singfac**2+reg_spot**2))
         IF (kin_flag) THEN
            xms_mn=xss_mn*(singfac**2/(singfac**2+reg_spot**2))
         ELSE
            xms_mn=-MATMUL(bmat,xmp1_mn)-MATMUL(cmat,xsp_mn)
         ENDIF
         bmt_mn=-(chi1*xmp1_mn+twopi*ifac*nn*xms_mn)
         bmz_mn=-(chi1*(q1*xsp_mn+sq%f(4)*xmp1_mn)+
     $        twopi*ifac*mfac*xms_mn)
      ELSE
         xmp1_mn=xsp1_mn
         xms_mn=xss_mn
         bmt_mn=bwt_mn
         bmz_mn=bwz_mn
      ENDIF

      DEALLOCATE(amat,bmat,cmat,fmats,gmats,kmats)
      IF(debug_flag) PRINT *, "->Leaving gpeq_sol"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_sol
c-----------------------------------------------------------------------
c     subprogram 2. gpeq_contra.
c     compute contravariant components of perturbed quantities.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_contra(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      
      INTEGER :: ipert,jpert,m1,dm
      COMPLEX(r8), DIMENSION(-mband:mband) ::jmat,jmat1

      IF(debug_flag) PRINT *, "Entering gpeq_contra"

      CALL spline_eval(sq,psi,1)
      CALL cspline_eval(metric%cs,psi,0)

      q=sq%f(4)
      q1=sq%f1(4)
      singfac=mfac-nn*q
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
      jmat(0:-mband:-1)=metric%cs%f(6*mband+7:7*mband+7)
      jmat1(0:-mband:-1)=metric%cs%f(7*mband+8:8*mband+8)
c-----------------------------------------------------------------------
c     compute upper half of matrices.
c-----------------------------------------------------------------------
      jmat(1:mband)=CONJG(jmat(-1:-mband:-1))
      jmat1(1:mband)=CONJG(jmat1(-1:-mband:-1))
c-----------------------------------------------------------------------
c     compute contravariant and modified quantities.
c-----------------------------------------------------------------------
      IF (reg_flag) THEN
         ipert=0
         xwp_mn=0
         xwt_mn=0
         xwz_mn=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               jpert=ipert+dm
               xwp_mn(ipert)=xwp_mn(ipert)+jmat(dm)*xsp_mn(jpert)
               xwt_mn(ipert)=xwt_mn(ipert)-(jmat(dm)*xmp1_mn(jpert)+
     $              jmat1(dm)*xsp_mn(jpert)+
     $              twopi*ifac*nn/chi1*jmat(dm)*xms_mn(jpert))/
     $              (twopi*ifac*(m1-nn*q))
               xwz_mn(ipert)=xwz_mn(ipert)-(q*jmat(dm)*xmp1_mn(jpert)+
     $              q*jmat1(dm)*xsp_mn(jpert)+
     $              twopi*ifac*m1/chi1*jmat(dm)*xms_mn(jpert))/
     $              (twopi*ifac*(m1-nn*q))
            ENDDO
         ENDDO
         xmt_mn=xwt_mn*(singfac**2/(singfac**2+reg_spot**2))
         xmz_mn=xwz_mn*(singfac**2/(singfac**2+reg_spot**2))
      ELSE
         ipert=0
         xwp_mn=0
         xwt_mn=0
         xwz_mn=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               jpert=ipert+dm
               xwp_mn(ipert)=xwp_mn(ipert)+jmat(dm)*xsp_mn(jpert)
               xwt_mn(ipert)=xwt_mn(ipert)-(jmat(dm)*xsp1_mn(jpert)+
     $              jmat1(dm)*xsp_mn(jpert)+
     $              twopi*ifac*nn/chi1*jmat(dm)*xss_mn(jpert))/
     $              (twopi*ifac*(m1-nn*q))
               xwz_mn(ipert)=xwz_mn(ipert)-(q*jmat(dm)*xsp1_mn(jpert)+
     $              q*jmat1(dm)*xsp_mn(jpert)+
     $              twopi*ifac*m1/chi1*jmat(dm)*xss_mn(jpert))/
     $              (twopi*ifac*(m1-nn*q))
            ENDDO
         ENDDO
         xmt_mn=xwt_mn
         xmz_mn=xwz_mn
      ENDIF
      IF(debug_flag) PRINT *, "->Leaving gpeq_contra"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_contra
c-----------------------------------------------------------------------
c     subprogram 3. gpeq_cova.
c     compute Jacobian-weighted covariant components from
c     Jacobian-weighted contravariant components.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_cova(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8) ,INTENT(IN) :: psi

      INTEGER :: ipert,jpert,m1,dm
      COMPLEX(r8), DIMENSION(-mband:mband) :: g11,g22,g33,g23,g31,g12

      IF(debug_flag) PRINT *, "Entering gpeq_cova"
      
      CALL spline_eval(sq,psi,1)
      CALL cspline_eval(metric%cs,psi,0)
      q=sq%f(4)
      singfac=mfac-nn*q
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
      g11(0:-mband:-1)=metric%cs%f(1:mband+1)
      g22(0:-mband:-1)=metric%cs%f(mband+2:2*mband+2)
      g33(0:-mband:-1)=metric%cs%f(2*mband+3:3*mband+3)
      g23(0:-mband:-1)=metric%cs%f(3*mband+4:4*mband+4)
      g31(0:-mband:-1)=metric%cs%f(4*mband+5:5*mband+5)
      g12(0:-mband:-1)=metric%cs%f(5*mband+6:6*mband+6)      
c-----------------------------------------------------------------------
c     compute upper half of matrices.
c-----------------------------------------------------------------------
      g11(1:mband)=CONJG(g11(-1:-mband:-1))
      g22(1:mband)=CONJG(g22(-1:-mband:-1))
      g33(1:mband)=CONJG(g33(-1:-mband:-1))
      g23(1:mband)=CONJG(g23(-1:-mband:-1))
      g31(1:mband)=CONJG(g31(-1:-mband:-1))
      g12(1:mband)=CONJG(g12(-1:-mband:-1))
c-----------------------------------------------------------------------
c     compute Jacobian-weighted covariant components with metric tensors.
c-----------------------------------------------------------------------
      ipert=0
      xvp_mn=0
      xvt_mn=0
      xvz_mn=0
      bvp_mn=0
      bvt_mn=0
      bvz_mn=0 
      DO m1=mlow,mhigh
         ipert=ipert+1
         DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
            jpert=ipert+dm
            xvp_mn(ipert)=xvp_mn(ipert)+g11(dm)*xwp_mn(jpert)+
     $           g12(dm)*xmt_mn(jpert)+g31(dm)*xmz_mn(jpert)
            xvt_mn(ipert)=xvt_mn(ipert)+g12(dm)*xwp_mn(jpert)+
     $           g22(dm)*xmt_mn(jpert)+g23(dm)*xmz_mn(jpert)
            xvz_mn(ipert)=xvz_mn(ipert)+g31(dm)*xwp_mn(jpert)+
     $           g23(dm)*xmt_mn(jpert)+g33(dm)*xmz_mn(jpert)
            bvp_mn(ipert)=bvp_mn(ipert)+g11(dm)*bwp_mn(jpert)+
     $           g12(dm)*bmt_mn(jpert)+g31(dm)*bmz_mn(jpert)
            bvt_mn(ipert)=bvt_mn(ipert)+g12(dm)*bwp_mn(jpert)+
     $           g22(dm)*bmt_mn(jpert)+g23(dm)*bmz_mn(jpert)
            bvz_mn(ipert)=bvz_mn(ipert)+g31(dm)*bwp_mn(jpert)+
     $           g23(dm)*bmt_mn(jpert)+g33(dm)*bmz_mn(jpert)
         ENDDO
      ENDDO
      IF(debug_flag) PRINT *, "->Leaving gpeq_cova"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_cova
c-----------------------------------------------------------------------
c     subprogram 4. gpeq_normal.
c     compute normal components.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_normal(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xno_fun,bno_fun
      IF(debug_flag) PRINT *, "Entering gpeq_normal"
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         jacs(itheta)=jac
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
      ENDDO
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
      xno_fun=xwp_fun/(jacs*delpsi)
      bno_fun=bwp_fun/(jacs*delpsi)
      CALL iscdftf(mfac,mpert,xno_fun,mthsurf,xno_mn)
      CALL iscdftf(mfac,mpert,bno_fun,mthsurf,bno_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_normal"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_normal
c-----------------------------------------------------------------------
c     subprogram 5. gpeq_tangent.
c     compute tangent components.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_tangent(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: jacs,bs,rfun,zfun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwt_fun,xvt_fun,xvz_fun,
     $     bwt_fun,bvt_fun,bvz_fun,xta_fun,bta_fun
      IF(debug_flag) PRINT *, "Entering gpeq_tangent"
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      q=sq%f(4)
      DO itheta=0,mthsurf
         CALL bicube_eval(eqfun,psi,theta(itheta),0)
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jacs(itheta)=rzphi%f(4)
         bs(itheta)=eqfun%f(1)
         v(2,1)=rzphi%fy(1)/(2*rfac)
         v(2,2)=(1+rzphi%fy(2))*twopi*rfac
         rfun=v(2,1)*cos(eta)-v(2,2)*sin(eta)
         zfun=v(2,1)*sin(eta)+v(2,2)*cos(eta)
      ENDDO
c-----------------------------------------------------------------------
c     compute tangential components, b times delpsi.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
      CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
      CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
      CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
      CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
      CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
      xta_fun=xwt_fun/jacs-(chi1/(jacs*bs))**2*(xvt_fun+q*xvz_fun)
      bta_fun=bwt_fun/jacs-(chi1/(jacs*bs))**2*(bvt_fun+q*bvz_fun)
      xta_fun=xta_fun*sqrt(rfun**2+zfun**2)
      bta_fun=bta_fun*sqrt(rfun**2+zfun**2)
      CALL iscdftf(mfac,mpert,xta_fun,mthsurf,xta_mn)
      CALL iscdftf(mfac,mpert,bta_fun,mthsurf,bta_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_tangent"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_tangent
c-----------------------------------------------------------------------
c     subprogram 6. gpeq_parallel.
c     compute parallel components for xi and b.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_parallel(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: eqb
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xvt_fun,bvt_fun,
     $     xvz_fun,bvz_fun,xpa_fun,bpa_fun
      IF(debug_flag) PRINT *, "Entering gpeq_parallel"
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      q=sq%f(4)
      DO itheta=0,mthsurf
         CALL bicube_eval(eqfun,psi,theta(itheta),0)
         eqb(itheta)=eqfun%f(1)
      ENDDO
c-----------------------------------------------------------------------
c     compute tangential components, b times delpsi.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
      CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
      CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
      CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
      xpa_fun=(xvt_fun+q*xvz_fun)/eqb
      bpa_fun=(bvt_fun+q*bvz_fun)/eqb
      CALL iscdftf(mfac,mpert,xpa_fun,mthsurf,xpa_mn)
      CALL iscdftf(mfac,mpert,bpa_fun,mthsurf,bpa_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_parallel"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_parallel
c-----------------------------------------------------------------------
c     subprogram 7. gpeq_rzphi.
c     compute rzphi components.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_rzphi(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: t11,t12,t21,t22,t33,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun,xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun

      IF(debug_flag) PRINT *, "Entering gpeq_rzphi"
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         jacs(itheta)=jac
         v(1,1)=rzphi%fx(1)/(2*rfac)
         v(1,2)=rzphi%fx(2)*twopi*rfac
         v(2,1)=rzphi%fy(1)/(2*rfac)
         v(2,2)=(1+rzphi%fy(2))*twopi*rfac
         v(3,3)=twopi*r(itheta)
         t11(itheta)=cos(eta)*v(1,1)-sin(eta)*v(1,2)
         t12(itheta)=cos(eta)*v(2,1)-sin(eta)*v(2,2)
         t21(itheta)=sin(eta)*v(1,1)+cos(eta)*v(1,2)
         t22(itheta)=sin(eta)*v(2,1)+cos(eta)*v(2,2)
         t33(itheta)=-1.0/v(3,3)
      ENDDO
c-----------------------------------------------------------------------
c     three vector components.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
      CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
      CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
      CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
      CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
      xrr_fun=(t11*xwp_fun+t12*xwt_fun)/jacs
      brr_fun=(t11*bwp_fun+t12*bwt_fun)/jacs
      xrz_fun=(t21*xwp_fun+t22*xwt_fun)/jacs
      brz_fun=(t21*bwp_fun+t22*bwt_fun)/jacs
      xrp_fun=t33*xvz_fun
      brp_fun=t33*bvz_fun
      CALL iscdftf(mfac,mpert,xrr_fun,mthsurf,xrr_mn)
      CALL iscdftf(mfac,mpert,brr_fun,mthsurf,brr_mn)
      CALL iscdftf(mfac,mpert,xrz_fun,mthsurf,xrz_mn)
      CALL iscdftf(mfac,mpert,brz_fun,mthsurf,brz_mn)
      CALL iscdftf(mfac,mpert,xrp_fun,mthsurf,xrp_mn)
      CALL iscdftf(mfac,mpert,brp_fun,mthsurf,brp_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_rzphi"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_rzphi
c-----------------------------------------------------------------------
c     subprogram 8. gpeq_surface.
c     compute surface currents and potentials.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_surface(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER :: i,j,itheta,rtheta

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: chi_fun,che_fun,kax_fun,
     $     xwp_fun,xwt_fun,bvt_fun,bvz_fun,bwp_fun,
     $     chy_fun,chx_fun,chw_fun,kaw_fun
      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
      COMPLEX(r8), DIMENSION(4,0:mthsurf) :: chp_fun,kap_fun

      REAL(r8), DIMENSION(:), POINTER :: 
     $     grri_real,grri_imag,grre_real,grre_imag,
     $     griw_real,griw_imag,grrw_real,grrw_imag

      ALLOCATE(grri_real(nths2),grri_imag(nths2),
     $     grre_real(nths2),grre_imag(nths2),
     $     griw_real(nths2),griw_imag(nths2),
     $     grrw_real(nths2),grrw_imag(nths2))
      IF(debug_flag) PRINT *, "Entering gpeq_surface"
c-----------------------------------------------------------------------
c     take into account reverse-theta in vacuum code.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
      CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
      CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
      CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
      
      CALL spline_eval(sq,psi,1)
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         dphi(itheta)=rzphi%f(3)
         jac=rzphi%f(4)
         jacs(itheta)=jac
      ENDDO
c-----------------------------------------------------------------------
c     compute vacuum magnetic potentials with reverse normal vector.
c-----------------------------------------------------------------------
      rbwp_mn=CONJG(bwp_mn)
      grri_real=MATMUL(grri,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      grri_imag=MATMUL(grri,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
      grre_real=MATMUL(grre,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      grre_imag=MATMUL(grre,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
      griw_real=MATMUL(griw,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      griw_imag=MATMUL(griw,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
      grrw_real=MATMUL(grrw,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      grrw_imag=MATMUL(grrw,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
c-----------------------------------------------------------------------
c     return into original coordinates.
c-----------------------------------------------------------------------
      DO itheta=0,mthsurf-1
         rtheta=mthsurf-itheta
         chi_fun(itheta+1)=(grri_real(rtheta)-ifac*grri_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
         che_fun(itheta+1)=(grre_real(rtheta)-ifac*grre_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
         chy_fun(itheta+1)=(griw_real(rtheta)-ifac*griw_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
         chx_fun(itheta+1)=(grrw_real(rtheta)-ifac*grrw_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
      ENDDO
      chi_fun(0)=chi_fun(mthsurf)
      che_fun(0)=che_fun(mthsurf)
      chy_fun(0)=chy_fun(mthsurf)
      chx_fun(0)=chx_fun(mthsurf)
      ! mutual inductance
      DO itheta=0,mthsurf-1
         rtheta=2*mthsurf-itheta
         chw_fun(itheta+1)=grrw_real(rtheta)-ifac*grrw_imag(rtheta)
      ENDDO
      chw_fun(0)=chw_fun(mthsurf)
c-----------------------------------------------------------------------
c     normalize chi functions of vacuum.
c-----------------------------------------------------------------------
      chi_fun=chi_fun/(twopi**2)
      che_fun=-che_fun/(twopi**2)
      chy_fun=chy_fun/(twopi**2)
      chx_fun=-chx_fun/(twopi**2)
      chw_fun=-chw_fun/(twopi**2)
      CALL iscdftf(mfac,mpert,chi_fun,mthsurf,chi_mn)
      CALL iscdftf(mfac,mpert,che_fun,mthsurf,che_mn)
      CALL iscdftf(mfac,mpert,chy_fun,mthsurf,chy_mn)
      CALL iscdftf(mfac,mpert,chx_fun,mthsurf,chx_mn)
      CALL iscdftf(mfac,mpert,chw_fun,mthsurf,chw_mn)
c-----------------------------------------------------------------------
c     compute plasma magnetic potential on the surface.
c-----------------------------------------------------------------------
      chp_fun(3,:)=-bvz_fun/(twopi*ifac*nn)
      chp_fun(1,:)=chp_fun(3,:)-sq%f1(1)*xwp_fun/(twopi*ifac*nn*jacs)
      chp_fun(4,:)=bvt_fun
      chp_fun(2,:)=chp_fun(4,:)-(sq%f1(1)*sq%f(4)+jacs*sq%f1(2)/chi1)*
     $           xwp_fun/jacs
      DO j=1,4
         CALL iscdftf(mfac,mpert,chp_fun(j,:),mthsurf,chp_mn(j,:))
      ENDDO
      DO i=1,mpert
         IF ((i+mlow-1).EQ.0) THEN
            chp_mn(2,i)=chp_mn(1,i)
            chp_mn(4,i)=chp_mn(3,i)
         ELSE
            chp_mn(2,i)=chp_mn(2,i)/(twopi*ifac*mfac(i))
            chp_mn(4,i)=chp_mn(4,i)/(twopi*ifac*mfac(i))
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     fourier transform and compute necessary matrices.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,chp_fun(2,:),mthsurf,chp_mn(2,:))
      CALL iscdftb(mfac,mpert,chp_fun(4,:),mthsurf,chp_mn(4,:))
      DO j=1,4
         kap_fun(j,:)=(chp_fun(j,:)-che_fun(:))/mu0
         CALL iscdftf(mfac,mpert,kap_fun(j,:),mthsurf,kap_mn(j,:))
      ENDDO
      kax_fun=(chi_fun-che_fun)/mu0
      kaw_fun=(chy_fun-chx_fun)/mu0
      CALL iscdftf(mfac,mpert,kax_fun,mthsurf,kax_mn)
      CALL iscdftf(mfac,mpert,kaw_fun,mthsurf,kaw_mn)
     
      DEALLOCATE(grri_real,grri_imag,grre_real,grre_imag)
      DEALLOCATE(griw_real,griw_imag,grrw_real,grrw_imag)
      IF(debug_flag) PRINT *, "->Leaving gpeq_surface"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_surface
c-----------------------------------------------------------------------
c     subprogram 9. gpeq_c.
c     compute C vector (covariant + contravariant) for EPF calculation.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_c(psi, ipsi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      INTEGER, INTENT(IN) :: ipsi

      INTEGER :: itheta, ipert, jpert, m1, dm
      
      REAL(r8) :: q, q1, f1raw, p1, chi1, jac
      REAL(r8) :: eta, rfac, v21, v22, v23, v33
      REAL(r8), DIMENSION(0:mthsurf) :: jacs, dphi, r_vec, z_vec

      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun, xno_fun, bno_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun, bmt_fun, bmz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bvp_fun, bvt_fun, bvz_fun

c     Local temporal arrays for computation
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cwp_fun, cwt_fun, cwz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cvp_fun, cvt_fun, cvz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cvp_funp, cvt_funp, cvz_funp

c     Supporting variables for Fourier reconstruction
      REAL(r8) :: jwt, jwz
      
c     Supporting functions
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi, dpdt, dpdz
      REAL(r8), DIMENSION(0:mthsurf) :: g_22, g_23, g_33
      COMPLEX(r8), DIMENSION(-mband:mband) :: g11,g22,g33,g23,g31,g12

      IF(debug_flag) PRINT *, "Entering gpeq_epf at ipsi=", ipsi
c-----------------------------------------------------------------------
c     1) Setup: equilibrium and metric at this psi.
c-----------------------------------------------------------------------
      CALL spline_eval(sq, psi, 1)
      CALL cspline_eval(metric%cs, psi, 0)
      q = sq%f(4)
      q1 = sq%f1(4)
      f1raw = sq%f1(1)
      p1 = sq%f1(2) / mu0
      chi1 = psio * twopi
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
c-----------------------------------------------------------------------
c     3) Reconstruct eigenfunctions in spatial representation.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac, mpert, bwp_fun, mthsurf, bwp_mn)
      CALL iscdftb(mfac, mpert, bmt_fun, mthsurf, bmt_mn)
      CALL iscdftb(mfac, mpert, bmz_fun, mthsurf, bmz_mn)
      CALL iscdftb(mfac, mpert, bvp_fun, mthsurf, bvp_mn)
      CALL iscdftb(mfac, mpert, bvt_fun, mthsurf, bvt_mn)
      CALL iscdftb(mfac, mpert, bvz_fun, mthsurf, bvz_mn)

c-----------------------------------------------------------------------
c     4) Extract contravariant components: divide by jacobian.
c     
c     bwp_mn, bwt_mn, bwz_mn are complex Fourier mode coefficients
c     (derived from complex eigenfunctions xsp_mn, xss_mn).
c     After IFFT via iscdftb(), they remain complex-valued functions
c     of theta. Dividing by real jacobian preserves complex nature.
c     This is essential because perturbations are inherently complex
c     functions in MHD stability analysis.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 1)
         jac = rzphi%f(4)
         rfac = SQRT(rzphi%f(1))
         eta = twopi*(theta(itheta) + rzphi%f(2))
         r_vec(itheta) = ro + rfac*COS(eta)
         z_vec(itheta) = zo + rfac*SIN(eta)
         jacs(itheta) = jac
         dphi(itheta) = rzphi%f(3)

         w(1,1) = (1.0+ rzphi%fy(2))*twopi**2*rfac*r_vec(itheta)/jac
         w(1,2) = -rzphi%fy(1)*pi*r_vec(itheta)/(rfac*jac)
         
         w(2,1) = -twopi**2*rfac*r_vec(itheta)*rzphi%fx(2)/jac
         w(2,2) = pi*r_vec(itheta)*rzphi%fx(1)/(rfac*jac)

         w(3,1) = (twopi*r_vec(itheta)*rfac/jac)*
     $          (rzphi%fx(2)*rzphi%fy(3)-rzphi%fx(3)*(1+rzphi%fy(2)))
         w(3,2) = (r_vec(itheta)/(2*rfac*jac))*
     $             (rzphi%fx(3)*rzphi%fy(1)-rzphi%fx(1)*rzphi%fy(3))

         delpsi(itheta) = SQRT(w(1,1)**2 + w(1,2)**2)
         dpdt(itheta) = w(1,1)*w(2,1) + w(1,2)*w(2,2)
         dpdz(itheta) = w(1,1)*w(3,1) + w(1,2)*w(3,2)

c     contravariant basis vectors (idcon_metric style - NO jac)
         v21 = rzphi%fy(1)/(2*rfac)
         v22 = (1+rzphi%fy(2))*twopi*rfac
         v23 = rzphi%fy(3)*r_vec(itheta)
         v33 = twopi*r_vec(itheta)

c     metric tensor: g_ij = sum(v_i * v_j)
         g_22(itheta) = (v21**2 + v22**2 + v23**2) 
         g_33(itheta) = (v33**2) 
         g_23(itheta) = (v23*v33) 
      ENDDO

      xno_fun=xwp_fun/(jacs*delpsi)
      bno_fun=bwp_fun/(jacs*delpsi)
c-----------------------------------------------------------------------
c     6) Compute C components.
c     bwp stores J Q^psi, while bmt/bmz store the modified
c     upper-family J Q^theta/J Q^zeta used consistently by gpeq_cova.
c
c     The metric stored in metric%cs is g_ij / J, so lowering via gpeq_cova
c     maps J Q^i -> Q_i. Therefore bv* and cv2* are lower/covariant
c     components without an extra Jacobian factor.
c
c     xwp_fun stores J xi^psi. Therefore:
c       - upper J C^i corrections use xwp_fun/(J |grad psi|^2),
c       - lower C_i corrections use xwp_fun/|grad psi|^2 or xwp_fun*j^i.
c     This matches Eqs. (148)-(151) together with the metric convention
c     used in idcon_metric.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         jac = jacs(itheta)
         jwt = - f1raw / jac
         jwz = - p1 * mu0 / chi1 - f1raw * q / jac

c     J C^psi = J Q^psi
         cwp_fun(itheta)= bwp_fun(itheta)

c     J C^theta = J Q^theta + xi^psi/|grad psi|^2
c     * [mu0 j^theta g_theta_zeta + mu0 j^zeta g_zeta_zeta]
c     = J Q^theta + xwp_fun/(J |grad psi|^2) * [...]
         cwt_fun(itheta) = bmt_fun(itheta) + xwp_fun(itheta)
     $      / ((delpsi(itheta)**2) * jac )*
     $     (jwt * g_23(itheta) + jwz * g_33(itheta))

c     J C^zeta = J Q^zeta - xi^psi/|grad psi|^2
c     * [mu0 j^theta g_theta_theta + mu0 j^zeta g_theta_zeta]
c     = J Q^zeta - xwp_fun/(J |grad psi|^2) * [...]
         cwz_fun(itheta) = bmz_fun(itheta) - xwp_fun(itheta)
     $     / ((delpsi(itheta)**2) * jac )*
     $     (jwt * g_22(itheta) + jwz * g_23(itheta))

c     C_psi = Q_psi + (J xi^psi)/|grad psi|^2
c     * [mu0 j^theta (grad psi.grad zeta) - mu0 j^zeta (grad psi.grad theta)]
         cvp_fun(itheta) = bvp_fun(itheta)  + 
     $        xwp_fun(itheta) / (delpsi(itheta)**2) *
     $        (jwt * dpdz(itheta) - jwz * dpdt(itheta))

c     C_theta = Q_theta + (J xi^psi) mu0 j^zeta
         cvt_fun(itheta) = bvt_fun(itheta) + 
     $        xwp_fun(itheta) * jwz

c     C_zeta = Q_zeta - (J xi^psi) mu0 j^theta
         cvz_fun(itheta) = bvz_fun(itheta) - 
     $        xwp_fun(itheta) * jwt

      ENDDO
      cwp_mn = 0
      cwt_mn = 0
      cwz_mn = 0
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
      g11(0:-mband:-1)=metric%cs%f(1:mband+1)
      g22(0:-mband:-1)=metric%cs%f(mband+2:2*mband+2)
      g33(0:-mband:-1)=metric%cs%f(2*mband+3:3*mband+3)
      g23(0:-mband:-1)=metric%cs%f(3*mband+4:4*mband+4)
      g31(0:-mband:-1)=metric%cs%f(4*mband+5:5*mband+5)
      g12(0:-mband:-1)=metric%cs%f(5*mband+6:6*mband+6)      
c-----------------------------------------------------------------------
c     compute upper half of matrices.
c-----------------------------------------------------------------------
      g11(1:mband)=CONJG(g11(-1:-mband:-1))
      g22(1:mband)=CONJG(g22(-1:-mband:-1))
      g33(1:mband)=CONJG(g33(-1:-mband:-1))
      g23(1:mband)=CONJG(g23(-1:-mband:-1))
      g31(1:mband)=CONJG(g31(-1:-mband:-1))
      g12(1:mband)=CONJG(g12(-1:-mband:-1))

      CALL iscdftf(mfac, mpert, cwp_fun, mthsurf,cwp_mn)
      CALL iscdftf(mfac, mpert, cwt_fun, mthsurf,cwt_mn)
      CALL iscdftf(mfac, mpert, cwz_fun, mthsurf,cwz_mn)
      CALL iscdftf(mfac, mpert, cvp_fun, mthsurf,cvp_mn)
      CALL iscdftf(mfac, mpert, cvt_fun, mthsurf,cvt_mn)
      CALL iscdftf(mfac, mpert, cvz_fun, mthsurf,cvz_mn)

      ipert = 0
      c2vp_mn = 0
      c2vt_mn = 0
      c2vz_mn = 0
      DO m1=mlow,mhigh
         ipert=ipert+1
         DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
            jpert=ipert+dm
            c2vp_mn(ipert)=c2vp_mn(ipert)+g11(dm)*cwp_mn(jpert)+
     $                    g12(dm)*cwt_mn(jpert)+g31(dm)*cwz_mn(jpert)
            c2vt_mn(ipert)=c2vt_mn(ipert)+g12(dm)*cwp_mn(jpert)+
     $                    g22(dm)*cwt_mn(jpert)+g23(dm)*cwz_mn(jpert)
            c2vz_mn(ipert)=c2vz_mn(ipert)+g31(dm)*cwp_mn(jpert)+
     $                    g23(dm)*cwt_mn(jpert)+g33(dm)*cwz_mn(jpert)
         ENDDO
      ENDDO

      CALL iscdftb(mfac, mpert, cvp_funp, mthsurf, c2vp_mn)
      CALL iscdftb(mfac, mpert, cvt_funp, mthsurf, c2vt_mn)
      CALL iscdftb(mfac, mpert, cvz_funp, mthsurf, c2vz_mn)

      IF(debug_flag) PRINT *, "->Leaving gpeq_c at ipsi=", ipsi
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_c
c-----------------------------------------------------------------------
c     subprogram 9c. gpeq_cveri.
c     verify the C-vector identity
c
c        (curl C) · grad(psi)
c        = (d C_zeta / d theta - d C_theta / d zeta) / J
c
c     using the single-n Fourier convention
c
c        exp[2*pi*i*(m*theta - n*zeta)].
c
c     Therefore
c
c        d/dtheta ->  2*pi*i*m
c        d/dzeta -> -2*pi*i*n
c
c     and
c
c        (curl C) · grad(psi)
c        = (d_theta C_zeta + 2*pi*i*n*C_theta)/J.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_cveri(psi, cveri_fun)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(0:mthsurf), INTENT(OUT) :: cveri_fun

      INTEGER :: itheta, ipert
      COMPLEX(r8), DIMENSION(mpert) :: dth_cvz_mn, dzt_cvt_mn
      COMPLEX(r8), DIMENSION(mpert) :: curlpsi_mn

      IF(debug_flag) PRINT *, "Entering gpeq_cveri"
c-----------------------------------------------------------------------
c     prepare psi-local equilibrium and reconstruct C.
c-----------------------------------------------------------------------
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)
      CALL gpeq_normal(psi)
      CALL gpeq_c(psi, 0)

c-----------------------------------------------------------------------
c     exact theta/zeta derivatives in mode space.
c-----------------------------------------------------------------------
      DO ipert = 1, mpert
         dth_cvz_mn(ipert) = twopi * ifac * mfac(ipert) * cvz_mn(ipert)
         dzt_cvt_mn(ipert) = -twopi * ifac * nn * cvt_mn(ipert)
         curlpsi_mn(ipert) = dth_cvz_mn(ipert) - dzt_cvt_mn(ipert)
      ENDDO

      CALL iscdftb(mfac, mpert, cveri_fun, mthsurf, curlpsi_mn)

c-----------------------------------------------------------------------
c     divide by Jacobian pointwise:
c        (curl C) · grad(psi) = (d_theta C_zeta - d_zeta C_theta)/J
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 0)
         jac = rzphi%f(4)
         cveri_fun(itheta) = cveri_fun(itheta) / jac
      ENDDO

      IF(debug_flag) PRINT *, "->Leaving gpeq_cveri"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_cveri
c-----------------------------------------------------------------------
c     subprogram 9d. gpeq_firstform.
c     evaluate the first-form plasma energy density components on one
c     flux surface with gamma term omitted:
c
c        |Q|^2 / mu0
c        j . (Q x xi*) / mu0
c        (div xi)* (xi . grad p)
c
c     Here xi.grad p is evaluated directly as p'(psi) xi^psi, and
c     div xi is evaluated directly from the Jacobian-weighted
c     contravariant displacement:
c
c        div xi = 1/J [ d(J xi^psi)/dpsi
c                      + d(J xi^theta)/dtheta
c                      + d(J xi^zeta)/dzeta ].
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_firstform(psi, q2_int, jqx_int, pdiv_int,
     $     total_int, q2_fun, jqx_fun, pdiv_fun, total_fun)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), INTENT(OUT) :: q2_int
      COMPLEX(r8), INTENT(OUT) :: jqx_int, pdiv_int, total_int
      COMPLEX(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     q2_fun, jqx_fun, pdiv_fun, total_fun

      INTEGER :: itheta
      REAL(r8) :: f1raw, p1
      REAL(r8) :: q2_density
      REAL(r8) :: mu0jtheta, mu0jzeta
      COMPLEX(r8) :: jqx_density, total_density, pdiv_density
      COMPLEX(r8) :: det_term, divxi_density, xigradp_density
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun, xmp1_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun, xmt_fun, xmz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun, bmt_fun, bmz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bvp_fun, bvt_fun, bvz_fun
      TYPE(cspline_type) :: divspl

      IF(debug_flag) PRINT *, "Entering gpeq_firstform"
c-----------------------------------------------------------------------
c     prepare psi-local perturbed equilibrium quantities.
c-----------------------------------------------------------------------
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)

      CALL spline_eval(sq, psi, 1)
      f1raw = sq%f1(1)
      p1 = sq%f1(2) / mu0

      CALL iscdftb(mfac, mpert, xsp_fun, mthsurf, xsp_mn)
      CALL iscdftb(mfac, mpert, xmp1_fun, mthsurf, xmp1_mn)
      CALL iscdftb(mfac, mpert, xwp_fun, mthsurf, xwp_mn)
      CALL iscdftb(mfac, mpert, xmt_fun, mthsurf, xmt_mn)
      CALL iscdftb(mfac, mpert, xmz_fun, mthsurf, xmz_mn)
      CALL iscdftb(mfac, mpert, bwp_fun, mthsurf, bwp_mn)
      CALL iscdftb(mfac, mpert, bmt_fun, mthsurf, bmt_mn)
      CALL iscdftb(mfac, mpert, bmz_fun, mthsurf, bmz_mn)
      CALL iscdftb(mfac, mpert, bvp_fun, mthsurf, bvp_mn)
      CALL iscdftb(mfac, mpert, bvt_fun, mthsurf, bvt_mn)
      CALL iscdftb(mfac, mpert, bvz_fun, mthsurf, bvz_mn)

      CALL cspline_alloc(divspl, mthsurf, 2)
      divspl%xs = theta
      DO itheta = 0, mthsurf
         divspl%fs(itheta,1) = xmt_fun(itheta)
         divspl%fs(itheta,2) = xmz_fun(itheta)
      ENDDO
      CALL cspline_fit(divspl, "periodic")

      q2_int = 0.0_r8
      jqx_int = CMPLX(0.0_r8, 0.0_r8, r8)
      pdiv_int = CMPLX(0.0_r8, 0.0_r8, r8)
      total_int = CMPLX(0.0_r8, 0.0_r8, r8)

      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 1)
         CALL cspline_eval(divspl, theta(itheta), 1)
         jac = rzphi%f(4)
         jac1 = rzphi%fx(4)

c        pointwise contravariant current components: mu0 j^theta,
c        mu0 j^zeta.  These match the definitions already used in gpeq_c.
         mu0jtheta = -f1raw / jac
         mu0jzeta = -mu0 * p1 / chi1 - sq%f(4) * f1raw / jac

         q2_density = REAL(CONJG(bwp_fun(itheta)) * bvp_fun(itheta) +
     $        CONJG(bmt_fun(itheta)) * bvt_fun(itheta) +
     $        CONJG(bmz_fun(itheta)) * bvz_fun(itheta), r8)
     $        / (mu0 * jac)

         det_term = -mu0jtheta *
     $        (bwp_fun(itheta) * CONJG(xmz_fun(itheta)) -
     $        bmz_fun(itheta) * CONJG(xwp_fun(itheta))) +
     $        mu0jzeta *
     $        (bwp_fun(itheta) * CONJG(xmt_fun(itheta)) -
     $        bmt_fun(itheta) * CONJG(xwp_fun(itheta)))
         jqx_density = det_term / (mu0 * jac)

         divxi_density = xmp1_fun(itheta) +
     $        (jac1 / jac) * xsp_fun(itheta) +
     $        divspl%f1(1) / jac -
     $        (twopi * ifac * nn) * divspl%f(2) / jac
         xigradp_density = p1 * xsp_fun(itheta)
         pdiv_density = CONJG(divxi_density) * xigradp_density
         total_density = CMPLX(q2_density, 0.0_r8, r8) -
     $        jqx_density + pdiv_density

         IF (PRESENT(q2_fun)) q2_fun(itheta) =
     $        CMPLX(q2_density, 0.0_r8, r8)
         IF (PRESENT(jqx_fun)) jqx_fun(itheta) = jqx_density
         IF (PRESENT(pdiv_fun)) pdiv_fun(itheta) = pdiv_density
         IF (PRESENT(total_fun)) total_fun(itheta) = total_density

         IF (itheta < mthsurf) THEN
            q2_int = q2_int + q2_density * jac / REAL(mthsurf, r8)
            jqx_int = jqx_int + jqx_density * jac / REAL(mthsurf, r8)
            pdiv_int = pdiv_int + pdiv_density * jac /
     $           REAL(mthsurf, r8)
            total_int = total_int + total_density * jac /
     $           REAL(mthsurf, r8)
         ENDIF
      ENDDO

      CALL cspline_dealloc(divspl)

      IF(debug_flag) PRINT *, "->Leaving gpeq_firstform"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_firstform
c-----------------------------------------------------------------------
c     subprogram 9b. gpeq_epf.
c     compute the flux-surface integral of the first EPF kernel,
c
c        \int dtheta dzeta J |C|^2 / mu0
c
c     and return the theta/zeta integrated value for one psi.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_epf(psi, epf_int, epf_p, epf_t, epf_z,
     $     epf_den_fun, epf_p_fun, epf_t_fun, epf_z_fun)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), INTENT(OUT) :: epf_int
      REAL(r8), OPTIONAL, INTENT(OUT) :: epf_p, epf_t, epf_z
c     optional per-theta C^2/mu0 densities for R-Z heatmaps. Same
c     convention as gpeq_recon3 epf_density: integral = sum den*jac/mth.
      REAL(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     epf_den_fun, epf_p_fun, epf_t_fun, epf_z_fun
      INTEGER :: itheta
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cwp_fun, cwt_fun, cwz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cvp_fun, cvt_fun, cvz_fun
      REAL(r8) :: epf_psi_int, epf_theta_int, epf_zeta_int
      REAL(r8) :: dpsi_d, dthe_d, dzet_d
      
      IF(debug_flag) PRINT *, "Entering gpeq_epf"
c-----------------------------------------------------------------------
c     Prepare psi-local perturbed equilibrium and compute C components.
c-----------------------------------------------------------------------
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)
      CALL gpeq_normal(psi)
      CALL gpeq_c(psi, 0)

      CALL iscdftb(mfac, mpert, cwp_fun,  mthsurf, cwp_mn)
      CALL iscdftb(mfac, mpert, cwt_fun,  mthsurf, cwt_mn)
      CALL iscdftb(mfac, mpert, cwz_fun,  mthsurf, cwz_mn)
      CALL iscdftb(mfac, mpert, cvp_fun, mthsurf, cvp_mn)
      CALL iscdftb(mfac, mpert, cvt_fun, mthsurf, cvt_mn)
      CALL iscdftb(mfac, mpert, cvz_fun, mthsurf, cvz_mn)

      epf_int = 0.0_r8
      epf_psi_int = 0.0_r8
      epf_theta_int = 0.0_r8
      epf_zeta_int = 0.0_r8
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 0)
         jac = rzphi%f(4)
c        per-theta densities (Re(conjg(cw_i) cv_i)/(mu0 jac)).
         dpsi_d = REAL(CONJG(cwp_fun(itheta)) * cvp_fun(itheta), r8)
     $        / (mu0 * jac)
         dthe_d = REAL(CONJG(cwt_fun(itheta)) * cvt_fun(itheta), r8)
     $        / (mu0 * jac)
         dzet_d = REAL(CONJG(cwz_fun(itheta)) * cvz_fun(itheta), r8)
     $        / (mu0 * jac)
         IF (PRESENT(epf_p_fun)) epf_p_fun(itheta) = dpsi_d
         IF (PRESENT(epf_t_fun)) epf_t_fun(itheta) = dthe_d
         IF (PRESENT(epf_z_fun)) epf_z_fun(itheta) = dzet_d
         IF (PRESENT(epf_den_fun)) epf_den_fun(itheta) =
     $        dpsi_d + dthe_d + dzet_d
         IF (itheta < mthsurf) THEN
            epf_psi_int = epf_psi_int + dpsi_d*jac/REAL(mthsurf, r8)
            epf_theta_int = epf_theta_int + dthe_d*jac/REAL(mthsurf, r8)
            epf_zeta_int = epf_zeta_int + dzet_d*jac/REAL(mthsurf, r8)
         ENDIF
      ENDDO
      epf_int = epf_psi_int + epf_theta_int + epf_zeta_int
      IF (PRESENT(epf_p)) epf_p = epf_psi_int
      IF (PRESENT(epf_t)) epf_t = epf_theta_int
      IF (PRESENT(epf_z)) epf_z = epf_zeta_int

      IF(debug_flag) PRINT *, "->Leaving gpeq_epf"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_epf
c-----------------------------------------------------------------------
c     subprogram 10. gpeq_dst.
c     compute the flux-surface integral of the destabilizing kernel,
c
c        \int dtheta dzeta J * K * xi_n^2
c
c     and return the theta/zeta integrated value for one psi.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_dst(psi, dst_int, dst_t1, dst_t2, dst_t3,
     $     dst1_den_fun, dst2_den_fun, dst3_den_fun)
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), INTENT(OUT) :: dst_int
      COMPLEX(r8), OPTIONAL, INTENT(OUT) :: dst_t1, dst_t2, dst_t3
c     optional per-theta K_i xi_n^2 densities for R-Z heatmaps. Same
c     convention as gpeq_recon3: integral = sum den*jac/mthsurf.
      REAL(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     dst1_den_fun, dst2_den_fun, dst3_den_fun
      INTEGER :: itheta
      COMPLEX(r8), DIMENSION(0:mthsurf) :: K_fun, xno_fun
      REAL(r8), DIMENSION(0:mthsurf) :: K_term1, K_term2, K_term3
      REAL(r8) :: xin2_fac, xin2_loc, d1, d2, d3
      COMPLEX(r8) :: dst1_int, dst2_int, dst3_int

c-----------------------------------------------------------------------
c     DST(psi) = \int dtheta dzeta J * K * xi_n^2
c-----------------------------------------------------------------------
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)
      CALL gpeq_normal(psi)
      CALL gpeq_K(psi, K_fun, K_term1, K_term2, K_term3)
      CALL iscdftb(mfac, mpert, xno_fun, mthsurf, xno_mn)
      dst_int = CMPLX(0.0_r8, 0.0_r8, r8)
      dst1_int = CMPLX(0.0_r8, 0.0_r8, r8)
      dst2_int = CMPLX(0.0_r8, 0.0_r8, r8)
      dst3_int = CMPLX(0.0_r8, 0.0_r8, r8)
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 0)
         jac = rzphi%f(4)
         xin2_loc = ABS(xno_fun(itheta))**2
         d1 = K_term1(itheta) * xin2_loc
         d2 = K_term2(itheta) * xin2_loc
         d3 = K_term3(itheta) * xin2_loc
         IF (PRESENT(dst1_den_fun)) dst1_den_fun(itheta) = d1
         IF (PRESENT(dst2_den_fun)) dst2_den_fun(itheta) = d2
         IF (PRESENT(dst3_den_fun)) dst3_den_fun(itheta) = d3
         IF (itheta < mthsurf) THEN
            xin2_fac = jac / REAL(mthsurf, r8)
            dst1_int = dst1_int + CMPLX(d1 * xin2_fac, 0.0_r8, r8)
            dst2_int = dst2_int + CMPLX(d2 * xin2_fac, 0.0_r8, r8)
            dst3_int = dst3_int + CMPLX(d3 * xin2_fac, 0.0_r8, r8)
         ENDIF
      ENDDO
      dst_int = dst1_int + dst2_int + dst3_int
      IF (PRESENT(dst_t1)) dst_t1 = dst1_int
      IF (PRESENT(dst_t2)) dst_t2 = dst2_int
      IF (PRESENT(dst_t3)) dst_t3 = dst3_int
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_dst
c-----------------------------------------------------------------------
c     subprogram 11. gpeq_shear.
c     compute magnetic shear at a single psi level.
c     approach: compute spatial domain first (theta), then FFT to mode space.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_shear(psi, shear_fun)
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(0:mthsurf), INTENT(OUT) :: shear_fun

      INTEGER :: itheta
      REAL(r8) :: r_val, theta_val
      REAL(r8) :: dpdp, dpdt, dpdz, shear_deriv
      REAL(r8) :: shear_contra
      TYPE(spline_type) :: shear_temp

c-----------------------------------------------------------------------
c     Setup: get equilibrium at this psi.
c-----------------------------------------------------------------------
      CALL spline_eval(sq, psi, 1)
      q = sq%f(4)
      q1 = sq%f1(4)

c-----------------------------------------------------------------------
c     Allocate spline for shear_temp with 5 components:
c     1: legacy covariant numerator term
c     2: (q*dpdt - dpdz) numerator
c     3: g_psi_g_psi denominator
c     4: legacy covariant denominator
c     5: DCON shear geometry term (component 2 / component 3)
c-----------------------------------------------------------------------
      CALL spline_alloc(shear_temp, mthsurf, 5)
      shear_temp%xs = theta(0:mthsurf)
      shear_temp%name = "shear_"

c-----------------------------------------------------------------------
c     Step 1: Compute shear components at all theta points.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         theta_val = theta(itheta)
         CALL bicube_eval(rzphi, psi, theta(itheta), 1)
         
         jac = rzphi%f(4)
         rfac = SQRT(MAX(rzphi%f(1), 0.0_r8))
         eta = twopi*(theta_val + rzphi%f(2))
         r_val = ro + rfac*COS(eta)

         
c        Compute w_ij contravariant metric in Cartesian components
         w(1,1) = (1.0_r8 + rzphi%fy(2))*twopi**2*rfac*r_val/jac
         w(1,2) = -rzphi%fy(1)*pi*r_val/(rfac*jac)
         
         w(2,1) = -twopi**2*rfac*r_val*rzphi%fx(2)/jac
         w(2,2) = pi*r_val*rzphi%fx(1)/(rfac*jac)
         
         w(3,1) = (twopi*r_val*rfac/jac)*
     $        (rzphi%fx(2)*rzphi%fy(3) - rzphi%fx(3)*(1.0_r8 + 
     $        rzphi%fy(2)))
         w(3,2) = (r_val/(2.0_r8*rfac*jac))*
     $        (rzphi%fx(3)*rzphi%fy(1) - rzphi%fx(1)*rzphi%fy(3))

c        Compute covariant basis v_ij = ∂r/∂ξ_j (Cartesian components)
c        From recon_metric:
         v(1,1) = rzphi%fx(1)/(2.0_r8*rfac*jac)
         v(1,2) = rzphi%fx(2)*twopi*rfac/jac
         v(1,3) = rzphi%fx(3)*r_val/jac

         v(2,1) = rzphi%fy(1)/(2.0_r8*rfac*jac)
         v(2,2) = (1.0_r8 + rzphi%fy(2))*twopi*rfac/jac
         v(2,3) = rzphi%fy(3)*r_val/jac
         
c        Component 33 of covariant basis
         v(3,3) = twopi*r_val/jac

c        Contravariant dot products: dpdp = ∇ψ·∇ψ, etc.
         dpdp = w(1,1)**2 + w(1,2)**2
         dpdt = w(1,1)*w(2,1) + w(1,2)*w(2,2)
         dpdz = w(1,1)*w(3,1) + w(1,2)*w(3,2)

         shear_temp%fs(itheta, 1) = twopi*r_val*jac*
     $        (-(v(1,1)*v(2,1)+v(1,2)*v(2,2))*(v(2,3)+q*v(3,3))
     $                +v(1,3)*(v(2,1)**2+v(2,2)**2))
         shear_temp%fs(itheta, 2) = (q*dpdt - dpdz)
         shear_temp%fs(itheta, 3) = dpdp
         shear_temp%fs(itheta, 4) = 
     $        twopi**2 * r_val**2 * (v(2,1)**2 + v(2,2)**2)


         shear_temp%fs(itheta, 5) = 
     $           shear_temp%fs(itheta, 2) / shear_temp%fs(itheta, 3)

      ENDDO
c-----------------------------------------------------------------------
c     Step 2: Fit spline to compute derivatives.
c-----------------------------------------------------------------------
      CALL spline_fit(shear_temp, "periodic")
c-----------------------------------------------------------------------
c     Step 3: Compute shear_fun = (chi1^2/J)*(q' + d shear_temp/d theta)
c     with shear_temp = (q*g^psi_theta - g^psi_zeta)/g^psi_psi.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         theta_val = theta(itheta)

         CALL spline_eval(shear_temp, theta_val, 1)
         CALL bicube_eval(rzphi, psi, theta_val, 1)
         jac = rzphi%f(4)
         
c        Contravariant approach: derivative of component 5
         shear_deriv = shear_temp%f1(5)
         shear_contra = (chi1**2/jac)*(q1 + shear_deriv)
         shear_fun(itheta) = CMPLX(shear_contra, 0.0_r8, r8)

      ENDDO

c-----------------------------------------------------------------------
c     Step 4: Transform to mode space via FFT.
c-----------------------------------------------------------------------
      CALL spline_dealloc(shear_temp)
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_shear
c-----------------------------------------------------------------------
c     subprogram 12. gpeq_curvature.
c     compute curvature in spatial domain.
c     κ·∇ψ = (|∇ψ|²/B²)[μ₀p' + (1/2)(∂B²/∂ψ) + (1/2)(∂B²/∂θ)(∇ψ·∇θ)/(∇ψ·∇ψ)]
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_curvature(psi, curv_fun)
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(0:mthsurf), INTENT(OUT) :: curv_fun

      INTEGER :: itheta
      REAL(r8) :: r_val, theta_val
      REAL(r8) :: delpsi
      REAL(r8) :: bsq_val, bsq_psi, bsq_theta
      REAL(r8) :: dpdt , kappa_psi

c-----------------------------------------------------------------------
c     Setup: get equilibrium at this psi.
c-----------------------------------------------------------------------
      CALL spline_eval(sq, psi, 1)
      p1 = sq%f1(2) / mu0

c-----------------------------------------------------------------------
c     Step 3: Compute curvature at all theta points.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         theta_val = theta(itheta)
         CALL bicube_eval(rzphi, psi, theta_val, 1)
         
         jac = rzphi%f(4)
         rfac = SQRT(rzphi%f(1))
         eta = twopi*(theta_val + rzphi%f(2))
         r_val = ro + rfac*COS(eta)

c        Compute contravariant metrics
         w(1,1) = (1.0+ rzphi%fy(2))*(twopi**2)*rfac*r_val/jac
         w(1,2) = -rzphi%fy(1)*pi*r_val/(rfac*jac)
         
         w(2,1) = -(twopi**2)*rfac*r_val*rzphi%fx(2)/jac
         w(2,2) = pi*r_val*rzphi%fx(1)/(rfac*jac)
         
         w(3,1) = (twopi*r_val*rfac/jac)*
     $        (rzphi%fx(2)*rzphi%fy(3) - rzphi%fx(3)*
     $        (1.0_r8 + rzphi%fy(2)))
         w(3,2) = (r_val/(2.0_r8*rfac*jac))*
     $        (rzphi%fx(3)*rzphi%fy(1) - rzphi%fx(1)*rzphi%fy(3))

c        Compute |∇ψ|² and ∇ψ·∇θ
         delpsi = SQRT(w(1,1)**2 + w(1,2)**2)
         dpdt = w(1,1)*w(2,1) + w(1,2)*w(2,2)

         CALL bicube_eval(eqfun, psi, theta_val, 1)
c        Retrieve B² and derivatives
         bsq_val = eqfun%f(1) ** 2
         bsq_psi = eqfun%fx(1) * 2 * eqfun%f(1)
         bsq_theta = eqfun%fy(1) * 2 * eqfun%f(1)

c        Compute κ·∇ψ safely

         kappa_psi = (delpsi**2 / bsq_val) *
     $           (p1*mu0 + 0.5_r8*bsq_psi +
     $           0.5_r8*bsq_theta*dpdt/(delpsi**2))

         curv_fun(itheta) = CMPLX(kappa_psi, 0.0_r8, r8)

      ENDDO
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_curvature
c-----------------------------------------------------------------------
c     subprogram 13. gpeq_K.
c     compute Bernstein K quantity (stability indicator).
c     K = |∇ψ_dcon|^2 * σ * S_dcon
c       + mu0 * B^2 * σ^2 + 2 p' * κ^psi
c     Uses pre-computed metric from idcon_metric
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_K(psi, K_fun, K_term1, K_term2, K_term3,
     $     sigma_fun, jdotb_fun, shear_out, curv_out)
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(0:mthsurf), INTENT(OUT) :: K_fun
      REAL(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     K_term1, K_term2, K_term3, sigma_fun, jdotb_fun
      COMPLEX(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     shear_out, curv_out

      INTEGER :: itheta
      REAL(r8) :: f1raw, delpsi
      REAL(r8) :: jwt, jwz, bth, bze, bsq_val, sigma, jdotb_val
      REAL(r8) :: g22, g23, g33, r_val
      REAL(r8), DIMENSION(0:mthsurf) :: K_t1, K_t2, K_t3
      REAL(r8), DIMENSION(0:mthsurf) :: sigma_vals, jdotb_vals
      COMPLEX(r8), DIMENSION(0:mthsurf) :: shear_fun, curv_fun

      chi1 = psio * twopi

      IF(debug_flag) PRINT *, "Entering gpeq_K"
c-----------------------------------------------------------------------
c     load equilibrium parameters
c-----------------------------------------------------------------------
      CALL spline_eval(sq, psi, 1)
      f1raw = sq%f1(1)
      p1 = sq%f1(2) / mu0
      q = sq%f(4)
c-----------------------------------------------------------------------
c     compute shear and curvature
c-----------------------------------------------------------------------
      CALL gpeq_shear(psi, shear_fun)
      CALL gpeq_curvature(psi, curv_fun)
c-----------------------------------------------------------------------
c     compute K in spatial domain (idcon_metric convention)
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 1)
         jac = rzphi%f(4)
         rfac = SQRT(rzphi%f(1))
         eta = twopi*(theta(itheta) + rzphi%f(2))
         r_val = ro + rfac*COS(eta)

c     |∇ψ| magnitude
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r_val/jac
         w(1,2)=-rzphi%fy(1)*pi*r_val/(rfac*jac)
         delpsi=SQRT(w(1,1)**2+w(1,2)**2)

c     current and magnetic field components
c     note: -f1raw/jac is the mu0*j^theta component (DCON convention),
c     so divide by mu0 to get the physical current j for K = mu0*sigma^2*B^2.
         jwt = -f1raw/(jac*mu0)
         jwz = q*jwt - p1/chi1
         bth = chi1 / jac
         bze = q * chi1 / jac

c     contravariant basis vectors (idcon_metric style - NO jac)
         v(2,1) = rzphi%fy(1)/(2*rfac)
         v(2,2) = (1+rzphi%fy(2))*twopi*rfac
         v(2,3) = rzphi%fy(3)*r_val
         
         v(3,3) = twopi*r_val

c     metric tensor: g_ij = sum(v_i * v_j)  (idcon_metric convention)
         g22 = (v(2,1)**2 + v(2,2)**2 + v(2,3)**2) 
         g33 = (v(3,3)**2) 
         g23 = (v(2,3)*v(3,3)) 

c     j·B and σ = (j·B)/B²
         CALL bicube_eval(eqfun, psi, theta(itheta), 0)
         bsq_val = eqfun%f(1)**2

         jdotb_val = g22*bth*jwt + g33*bze*jwz +
     $        g23*(bth*jwz + bze*jwt)
         sigma = jdotb_val / bsq_val
         sigma_vals(itheta) = sigma
         jdotb_vals(itheta) = jdotb_val

c     Term1: |∇ψ_dcon|^2 * σ * S_dcon
         K_t1(itheta) = (delpsi**2) * sigma * REAL(shear_fun(itheta))

c     Term2: mu0 * B^2 * σ^2
         K_t2(itheta) = bsq_val * sigma**2 * mu0

c     Term3: 2 p' * κ^psi
         K_t3(itheta) = 2.0_r8 * REAL(curv_fun(itheta)) * p1
      ENDDO

c     total K = T1 + T2 + T3
      DO itheta = 0, mthsurf
         K_fun(itheta) = CMPLX(K_t1(itheta) + K_t2(itheta) +
     $                         K_t3(itheta), 0.0_r8, r8)
      ENDDO

c     optionally return individual term values
      IF (PRESENT(K_term1)) K_term1 = K_t1
      IF (PRESENT(K_term2)) K_term2 = K_t2
      IF (PRESENT(K_term3)) K_term3 = K_t3
      IF (PRESENT(sigma_fun)) sigma_fun = sigma_vals
      IF (PRESENT(jdotb_fun)) jdotb_fun = jdotb_vals
      IF (PRESENT(shear_out)) shear_out = shear_fun
      IF (PRESENT(curv_out)) curv_out = curv_fun

      IF(debug_flag) PRINT *, "Exiting gpeq_K"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_K
c-----------------------------------------------------------------------
c     subprogram 13b. gpeq_recon3.
c     Bernstein-form |C|^2 decomposition for recon3 diagnostic.
c     Vector identity: C = Q + V with V = xi_n * mu0 * j x n_hat.
c     |C|^2 = |Q|^2 + 2 Re(V*.Q) + |V|^2
c     j is purely tangent to flux surface (j.n=0), so
c     |V|^2 = mu0^2 xi_n^2 |j|^2 = mu0^2 xi_n^2 (sigma^2 B^2
c           + p'^2 |grad psi|^2 / B^2).
c     Returns the surface integrals (theta only, jac measure):
c        A_int    = int |Q|^2 / mu0 * J dtheta / mthsurf
c        B_int    = int 2 Re(V*.Q) / mu0 * J dtheta / mthsurf
c                   directly via V = C - Q (covariant pairing)
c        Cpar_int = int mu0 sigma^2 B^2 xi_n^2 * J dtheta / mthsurf
c                   (= K_2 xi_n^2 integral; matches gpeq_dst's dst2)
c        Iperp_int= int mu0 p'^2 |grad psi|^2 xi_n^2 / B^2 *
c                   J dtheta / mthsurf
c        epf_int  = int |C|^2 / mu0 * J dtheta / mthsurf
c     Identity check: epf_int = A_int + B_int + Cpar_int + Iperp_int
c     should hold within numerical roundoff.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_recon3(psi, A_int, B_int, Iperp_int, Cpar_int,
     $     epf_int, A_fun, B_fun, Iperp_fun, Cpar_fun, surf_extra)
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), INTENT(OUT) :: A_int, B_int, Iperp_int, Cpar_int,
     $     epf_int
      REAL(r8), DIMENSION(0:mthsurf), OPTIONAL, INTENT(OUT) ::
     $     A_fun, B_fun, Iperp_fun, Cpar_fun
c     surf_extra: packed per-theta fields for R-Z heatmaps (recon_out).
c     columns: 1 VdotQ_re, 2 VdotQ_im, 3 Bcur_den, 4 Bpre_den,
c       5/6 Qp_re/im, 7/8 Qt_re/im, 9/10 Qz_re/im,
c       11/12 Vt_re/im, 13/14 Vz_re/im, 15 jpar, 16 jperp,
c       17/18 xin_re/im, 19 delpsi, 20 Bmod,
c       21 dst1_den (K1 xin2), 22 dst3_den (K3 xin2).
      REAL(r8), DIMENSION(0:mthsurf,22), OPTIONAL, INTENT(OUT) ::
     $     surf_extra

      INTEGER :: itheta
      REAL(r8) :: f1raw, p1_local
      REAL(r8) :: A_density, B_density, Iperp_density, Cpar_density,
     $     epf_density
      REAL(r8) :: delpsi_val, bsq_val, sigma_local, jdotb_val
      REAL(r8) :: jwt, jwz, bth, bze, g22, g33, g23, r_val
      REAL(r8) :: xin2, mu0jwt, mu0jwz
      REAL(r8) :: mu0jpar_t, mu0jpar_z, mu0jper_t, mu0jper_z
      REAL(r8) :: jsq_val, jpar_val, bmod_val
      COMPLEX(r8) :: jvt_loc, jvz_loc, vdotq_loc
      COMPLEX(r8) :: jvt_cur, jvz_cur, jvt_pre, jvz_pre
      COMPLEX(r8), DIMENSION(0:mthsurf) :: kfun_loc
      REAL(r8), DIMENSION(0:mthsurf) :: kt1_loc, kt2_loc, kt3_loc
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun, bmt_fun, bmz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bvp_fun, bvt_fun, bvz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cwp_fun, cwt_fun, cwz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: cvp_fun, cvt_fun, cvz_fun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xno_fun, xwp_fun

      IF(debug_flag) PRINT *, "Entering gpeq_recon3"
c-----------------------------------------------------------------------
c     prepare psi-local perturbed equilibrium and C, Q.
c-----------------------------------------------------------------------
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)
      CALL gpeq_normal(psi)
      CALL gpeq_c(psi, 0)

      CALL spline_eval(sq, psi, 1)
      f1raw = sq%f1(1)
      p1_local = sq%f1(2) / mu0
      chi1 = psio * twopi

      CALL iscdftb(mfac, mpert, xno_fun, mthsurf, xno_mn)
      CALL iscdftb(mfac, mpert, bwp_fun, mthsurf, bwp_mn)
      CALL iscdftb(mfac, mpert, bmt_fun, mthsurf, bmt_mn)
      CALL iscdftb(mfac, mpert, bmz_fun, mthsurf, bmz_mn)
      CALL iscdftb(mfac, mpert, bvp_fun, mthsurf, bvp_mn)
      CALL iscdftb(mfac, mpert, bvt_fun, mthsurf, bvt_mn)
      CALL iscdftb(mfac, mpert, bvz_fun, mthsurf, bvz_mn)
      CALL iscdftb(mfac, mpert, cwp_fun, mthsurf, cwp_mn)
      CALL iscdftb(mfac, mpert, cwt_fun, mthsurf, cwt_mn)
      CALL iscdftb(mfac, mpert, cwz_fun, mthsurf, cwz_mn)
      CALL iscdftb(mfac, mpert, cvp_fun, mthsurf, cvp_mn)
      CALL iscdftb(mfac, mpert, cvt_fun, mthsurf, cvt_mn)
      CALL iscdftb(mfac, mpert, cvz_fun, mthsurf, cvz_mn)
      CALL iscdftb(mfac, mpert, xwp_fun, mthsurf, xwp_mn)

c     K term breakdown (K1 shear, K3 curvature) for dst1/dst3 densities;
c     only needed for the R-Z heatmap output.
      IF (PRESENT(surf_extra)) THEN
         CALL gpeq_K(psi, kfun_loc, kt1_loc, kt2_loc, kt3_loc)
      ENDIF

      A_int = 0.0_r8
      B_int = 0.0_r8
      Iperp_int = 0.0_r8
      Cpar_int = 0.0_r8
      epf_int = 0.0_r8
c-----------------------------------------------------------------------
c     accumulate densities in theta with jac measure.
c-----------------------------------------------------------------------
      DO itheta = 0, mthsurf
         CALL bicube_eval(rzphi, psi, theta(itheta), 1)
         CALL bicube_eval(eqfun, psi, theta(itheta), 0)
         jac = rzphi%f(4)
         rfac = SQRT(rzphi%f(1))
         eta = twopi*(theta(itheta) + rzphi%f(2))
         r_val = ro + rfac*COS(eta)
         bsq_val = eqfun%f(1)**2

c        |grad psi|
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r_val/jac
         w(1,2)=-rzphi%fy(1)*pi*r_val/(rfac*jac)
         delpsi_val=SQRT(w(1,1)**2+w(1,2)**2)

c        currents and B (idcon_metric convention, same as gpeq_K).
         jwt = -f1raw/(jac*mu0)
         jwz = sq%f(4)*jwt - p1_local/chi1
         bth = chi1 / jac
         bze = sq%f(4) * chi1 / jac

         v(2,1) = rzphi%fy(1)/(2*rfac)
         v(2,2) = (1+rzphi%fy(2))*twopi*rfac
         v(2,3) = rzphi%fy(3)*r_val
         v(3,3) = twopi*r_val
         g22 = (v(2,1)**2 + v(2,2)**2 + v(2,3)**2)
         g33 = (v(3,3)**2)
         g23 = (v(2,3)*v(3,3))

         jdotb_val = g22*bth*jwt + g33*bze*jwz +
     $        g23*(bth*jwz + bze*jwt)
         sigma_local = jdotb_val / bsq_val

         xin2 = REAL(CONJG(xno_fun(itheta))*xno_fun(itheta), r8)

c        |Q|^2 / mu0 density (matches gpeq_firstform q2_density)
         A_density = REAL(CONJG(bwp_fun(itheta)) * bvp_fun(itheta) +
     $        CONJG(bmt_fun(itheta)) * bvt_fun(itheta) +
     $        CONJG(bmz_fun(itheta)) * bvz_fun(itheta), r8)
     $        / (mu0 * jac)

c        |C|^2 / mu0 density (matches gpeq_epf integrand)
         epf_density = REAL(
     $        CONJG(cwp_fun(itheta)) * cvp_fun(itheta) +
     $        CONJG(cwt_fun(itheta)) * cvt_fun(itheta) +
     $        CONJG(cwz_fun(itheta)) * cvz_fun(itheta), r8)
     $        / (mu0 * jac)

c        B density = 2 Re(V*.Q) / mu0, with V = xi_n (mu0 j x n_hat)
c        computed DIRECTLY from the equilibrium current, NOT from C - Q.
c        These are the same C^i - Q^i terms derived in main.tex (Sec. 5,
c        the C component box): J V^psi = 0, and with xwp_fun = J xi^psi,
c          J V^theta =  xwp/(|grad psi|^2 J) (mu0 j^theta g23 + mu0 j^zeta g33)
c          J V^zeta  = -xwp/(|grad psi|^2 J) (mu0 j^theta g22 + mu0 j^zeta g23).
c        mu0 j^i = mu0 * (physical j^i); jwt,jwz already hold physical j^i.
         mu0jwt = mu0 * jwt
         mu0jwz = mu0 * jwz
         jvt_loc = xwp_fun(itheta) / (delpsi_val**2 * jac) *
     $        (mu0jwt * g23 + mu0jwz * g33)
         jvz_loc = -xwp_fun(itheta) / (delpsi_val**2 * jac) *
     $        (mu0jwt * g22 + mu0jwz * g23)
c        V is contravariant (J V^i); pair with covariant Q_i (bvt,bvz).
         B_density = 2.0_r8 * REAL(
     $        CONJG(jvt_loc) * bvt_fun(itheta) +
     $        CONJG(jvz_loc) * bvz_fun(itheta), r8) / (mu0 * jac)

c        K_2 xi_n^2 density (parallel current; should match dst2)
         Cpar_density = mu0 * sigma_local**2 * bsq_val * xin2

c        I_perp density (perpendicular Pfirsch-Schlueter current^2)
         Iperp_density = mu0 * p1_local**2 * delpsi_val**2 *
     $        xin2 / bsq_val

c        Extra per-theta fields for R-Z heatmaps (only if requested).
         IF (PRESENT(surf_extra)) THEN
            bmod_val = SQRT(bsq_val)
c           V*.Q (complex); B_density = 2 Re(vdotq_loc)/mu0.
            vdotq_loc = (CONJG(jvt_loc) * bvt_fun(itheta) +
     $           CONJG(jvz_loc) * bvz_fun(itheta)) / jac
c           Split V = V_par + V_perp via j = j_par + j_perp.
c           Parallel current: mu0 j_par^i = mu0 sigma B^i (B^t=bth,B^z=bze).
            mu0jpar_t = mu0 * sigma_local * bth
            mu0jpar_z = mu0 * sigma_local * bze
            mu0jper_t = mu0jwt - mu0jpar_t
            mu0jper_z = mu0jwz - mu0jpar_z
            jvt_cur = xwp_fun(itheta) / (delpsi_val**2 * jac) *
     $           (mu0jpar_t * g23 + mu0jpar_z * g33)
            jvz_cur = -xwp_fun(itheta) / (delpsi_val**2 * jac) *
     $           (mu0jpar_t * g22 + mu0jpar_z * g23)
            jvt_pre = jvt_loc - jvt_cur
            jvz_pre = jvz_loc - jvz_cur
c           j_par = sigma B; |j_perp|^2 = |j|^2 - j_par^2 (j physical).
            jpar_val = sigma_local * bmod_val
            jsq_val = g22 * jwt**2 + g33 * jwz**2 +
     $           2.0_r8 * g23 * jwt * jwz
            surf_extra(itheta,1) = REAL(vdotq_loc, r8)
            surf_extra(itheta,2) = AIMAG(vdotq_loc)
            surf_extra(itheta,3) = 2.0_r8 * REAL(
     $           CONJG(jvt_cur) * bvt_fun(itheta) +
     $           CONJG(jvz_cur) * bvz_fun(itheta), r8) / (mu0 * jac)
            surf_extra(itheta,4) = 2.0_r8 * REAL(
     $           CONJG(jvt_pre) * bvt_fun(itheta) +
     $           CONJG(jvz_pre) * bvz_fun(itheta), r8) / (mu0 * jac)
            surf_extra(itheta,5) = REAL(bvp_fun(itheta), r8)
            surf_extra(itheta,6) = AIMAG(bvp_fun(itheta))
            surf_extra(itheta,7) = REAL(bvt_fun(itheta), r8)
            surf_extra(itheta,8) = AIMAG(bvt_fun(itheta))
            surf_extra(itheta,9) = REAL(bvz_fun(itheta), r8)
            surf_extra(itheta,10) = AIMAG(bvz_fun(itheta))
            surf_extra(itheta,11) = REAL(jvt_loc, r8)
            surf_extra(itheta,12) = AIMAG(jvt_loc)
            surf_extra(itheta,13) = REAL(jvz_loc, r8)
            surf_extra(itheta,14) = AIMAG(jvz_loc)
            surf_extra(itheta,15) = jpar_val
            surf_extra(itheta,16) = SQRT(MAX(0.0_r8,
     $           jsq_val - jpar_val**2))
            surf_extra(itheta,17) = REAL(xno_fun(itheta), r8)
            surf_extra(itheta,18) = AIMAG(xno_fun(itheta))
            surf_extra(itheta,19) = delpsi_val
            surf_extra(itheta,20) = bmod_val
c           dst1 = K1 xin2 (shear), dst3 = K3 xin2 (curvature).
            surf_extra(itheta,21) = kt1_loc(itheta) * xin2
            surf_extra(itheta,22) = kt3_loc(itheta) * xin2
         ENDIF

         IF (PRESENT(A_fun)) A_fun(itheta) = A_density
         IF (PRESENT(B_fun)) B_fun(itheta) = B_density
         IF (PRESENT(Iperp_fun)) Iperp_fun(itheta) = Iperp_density
         IF (PRESENT(Cpar_fun)) Cpar_fun(itheta) = Cpar_density

         IF (itheta < mthsurf) THEN
            A_int = A_int + A_density * jac / REAL(mthsurf, r8)
            B_int = B_int + B_density * jac / REAL(mthsurf, r8)
            epf_int = epf_int + epf_density * jac / REAL(mthsurf, r8)
            Cpar_int = Cpar_int + Cpar_density * jac /
     $           REAL(mthsurf, r8)
            Iperp_int = Iperp_int + Iperp_density * jac /
     $           REAL(mthsurf, r8)
         ENDIF
      ENDDO

      IF(debug_flag) PRINT *, "->Leaving gpeq_recon3"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_recon3
c-----------------------------------------------------------------------
c     subprogram 14. gpeq_fcoords.
c     transform coordinates to dcon coordinates.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_fcoords(psi,ftnmn,amf,amp,ri,bpi,bi,rci,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,ri,bpi,bi,rci,ti,ji
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      LOGICAL :: first = .TRUE.
      INTEGER :: i,itheta
      REAL(r8) :: thetai,jarea,psave=0
      INTEGER, DIMENSION(6) :: isave=0,itmp=0

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dphi,thetas,jacfac
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      TYPE(spline_type) :: spl

      ! note automatic arrays are allocated and deallocated on entry/exit
      ! instead, we use allocatables and just allocate once for all
      SAVE :: first,psave,isave,jarea,dphi,thetas,jacfac
      IF(first) ALLOCATE(dphi(0:mthsurf),thetas(0:mthsurf),
     $      jacfac(0:mthsurf))
      first  = .FALSE.

      IF(debug_flag) PRINT *, "Entering gpeq_fcoords"

      ! global sq may have been eval'd elsewhere inbetween bcoords calls
      CALL spline_eval(sq,psi,0)
      ! expensive spline formation, do only if asking for new coords or surface
      itmp = (/ri,bpi,bi,rci,ti,ji/)
      IF(.NOT.ALL(itmp==isave).OR. (psave/=psi))THEN
         isave = (/ri,bpi,bi,rci,ti,ji/)
         psave = psi
         dphi   = 0
         thetas = 0
         jacfac = 0

         CALL spline_alloc(spl,mthsurf,2)
         spl%xs=theta

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            bpfac=psio*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            spl%fs(itheta,1)=fac/(r(itheta)**ri*rfac**rci)*
     $           bpfac**bpi*bfac**bi
            ! jacobian for coordinate angle at dcon angle
            spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**ri*rfac**rci/
     $           (bpfac**bpi*bfac**bi)
            IF (ti .EQ. 0) THEN
               dphi(itheta)=rzphi%f(3)
            ENDIF
         ENDDO

         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)

         ! coordinate angle at dcon angle
         thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
         IF (ji /= 0) THEN
            DO itheta=0,mthsurf
               ! jacobian at coordinate angle
               thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
               CALL spline_eval(spl,thetai,0)
               jacfac(itheta) = spl%f(2)
            ENDDO
            SELECT CASE(ji)
               CASE(-2)
                  jacfac = 1/sqrt(jacfac)
               CASE(-1)
                  jacfac = 1/jacfac
               CASE(1)
                  jacfac = jacfac
               CASE(2)
                  jacfac = sqrt(jacfac)
            END SELECT
            ! surface area
            jarea=0
            DO itheta=0,mthsurf-1
               jarea=jarea+jacfac(itheta)/mthsurf
            ENDDO
         ENDIF
         CALL spline_dealloc(spl)
      ENDIF
      ! convert to unweighted spectrum
      IF (ji /= 0) THEN
         CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)
         ftnfun=ftnfun/jacfac*jarea
         CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------      
      ! compute given function in dcon angle
      DO itheta=0,mthsurf
         ftnfun(itheta)=0
         DO i=1,amp
            ftnfun(itheta)=ftnfun(itheta)+
     $           ftnmn(i)*EXP(ifac*twopi*amf(i)*thetas(itheta))
         ENDDO
      ENDDO

      ! multiply toroidal factor for dcon angle
      IF (ti .EQ. 0) THEN
         ftnfun(:)=ftnfun(:)*EXP(-ifac*nn*dphi(:))
      ELSE
         ftnfun(:)=ftnfun(:)*
     $        EXP(-twopi*ifac*nn*sq%f(4)*(thetas(:)-theta(:)))
      ENDIF

      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_fcoords"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_fcoords
c-----------------------------------------------------------------------
c     subprogram 14. gpeq_fcoordsout.
c     transform to dcon coordinates. Assumes mpert,lmpert,jac_out
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_fcoordsout(fmo,fmi,psi,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN), OPTIONAL :: ti,ji
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(lmpert), INTENT(IN) :: fmi
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: fmo

      COMPLEX(r8), DIMENSION(lmpert) :: tmp
      INTEGER :: i,tout,jout
c-----------------------------------------------------------------------
c     Defaults.
c-----------------------------------------------------------------------
      IF(PRESENT(ti))THEN
         tout = ti
      ELSE
         tout = tmag_out
      ENDIF
      IF(PRESENT(ji))THEN
         jout = ji
      ELSE
         jout = 0
      ENDIF
      fmo = 0
c-----------------------------------------------------------------------
c     Transform new vector in larger m space and then transfer to mo
c-----------------------------------------------------------------------
      IF(mpert>lmpert)THEN
         ! transfer to mo space
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               fmo(i) = fmi(mlow-lmlow+i)
            ENDIF
         ENDDO
         ! transform new vector
         IF((jac_out/=jac_type).OR.(tout==0).OR.(jout/=0))THEN
            CALL gpeq_fcoords(psi,fmo,mfac,mpert,power_rout,power_bpout,
     $         power_bout,power_rcout,tout,jout)
         ENDIF
      ELSE
         ! transform in mi space
         tmp = fmi
         IF((jac_out/=jac_type).OR.(tout==0).OR.(jout/=0))THEN
            CALL gpeq_fcoords(psi,tmp,lmfac,lmpert,power_rout,
     $         power_bpout,power_bout,power_rcout,tout,jout)
         ENDIF
         ! transfer to mo space
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               fmo(i) = tmp(mlow-lmlow+i)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_fcoordsout
c-----------------------------------------------------------------------
c     subprogram 15. gpeq_bcoords.
c     transform dcon coordinates to other coordinates.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_bcoords(psi,ftnmn,amf,amp,ri,bpi,bi,rci,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: amp,ri,bpi,bi,rci,ti,ji
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      LOGICAL :: first = .TRUE.
      INTEGER :: i,itheta
      INTEGER, DIMENSION(6) :: isave=0,itmp=0
      REAL(r8) :: thetai,jarea,psave=0

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dphi,thetas,jacfac
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      TYPE(spline_type) :: spl       

      ! note automatic arrays are allocated and deallocated on entry/exit
      ! instead, we use allocatables and just allocate once for all
      SAVE :: first,psave,isave,jarea,spl,dphi,thetas,jacfac
      IF(first) ALLOCATE(dphi(0:mthsurf),thetas(0:mthsurf),
     $   jacfac(0:mthsurf))
      first = .FALSE.

      IF(debug_flag) PRINT *, "Entering gpeq_bcoords"
      
      ! global sq may have been eval'd elsewhere inbetween bcoords calls
      CALL spline_eval(sq,psi,0)
      ! expensive spline formation, do only if asking for new bcoords
      itmp = (/ri,bpi,bi,rci,ti,ji/)
      IF(.NOT.ALL(itmp==isave).OR. psave/=psi)THEN
         isave = (/ri,bpi,bi,rci,ti,ji/)
         psave = psi
         dphi   = 0
         thetas = 0
         jacfac = 0

         CALL spline_alloc(spl,mthsurf,2)
         spl%xs=theta

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            bpfac=psio*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            spl%fs(itheta,1)=fac/(r(itheta)**ri*rfac**rci)*
     $           bpfac**bpi*bfac**bi
            spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**ri*rfac**rci/
     $           (bpfac**bpi*bfac**bi)
            IF (ti .EQ. 0) THEN
               dphi(itheta)=rzphi%f(3)
            ENDIF
         ENDDO

         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)

         ! coordinate angle at dcon angle
         thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
         DO itheta=0,mthsurf
            ! dcon angle at coordinate angle
            thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
            ! jacobian at coordinate angle
            CALL spline_eval(spl,thetai,0)
            jacfac(itheta) = spl%f(2)
         ENDDO
         SELECT CASE(ji)
            CASE(-2)
               jacfac = 1/sqrt(jacfac)
            CASE(-1)
               jacfac = 1/jacfac
            CASE(1)
               jacfac = jacfac
            CASE(2)
               jacfac = sqrt(jacfac)
         END SELECT
         ! surface area
         jarea=0
         DO itheta=0,mthsurf-1
            jarea=jarea+jacfac(itheta)/mthsurf
         ENDDO
         CALL spline_dealloc(spl)
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      ! functional form of dcon fourier vector
      CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)
      ! take toroidal factor from dcon angle
      IF (ti .EQ. 0) THEN
         ftnfun=ftnfun*EXP(ifac*nn*dphi)
      ENDIF
      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)

      ! compute given function in coordinate angle
      DO itheta=0,mthsurf
         ftnfun(itheta)=0
         ! dcon angle at coordinate angle
         thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
         DO i=1,amp
            ftnfun(itheta)=ftnfun(itheta)+
     $           ftnmn(i)*EXP(ifac*twopi*amf(i)*thetai)
         ENDDO
         ! take toroidal factor back for coordinate angle
         IF (ti .NE. 0) THEN
            ftnfun(itheta)=ftnfun(itheta)*
     $           EXP(-twopi*ifac*nn*sq%f(4)*(thetai-theta(itheta)))
         ENDIF
      ENDDO

      ! optional jacobian wieghting
      IF (ji /= 0) THEN
         ftnfun=ftnfun*jacfac/jarea
      ENDIF
      ! forward transform function to coordinate fourier spectrum
      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_bcoords"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_bcoords
c-----------------------------------------------------------------------
c     subprogram 16. gpeq_bcoordsout.
c     transform dcon to other coordinates. Assumes mpert,lmpert,jac_out
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_bcoordsout(fmo,fmi,psi,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN), OPTIONAL :: ti,ji
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: fmi
      COMPLEX(r8), DIMENSION(lmpert), INTENT(OUT) :: fmo

      COMPLEX(r8), DIMENSION(mpert) :: tmp
      INTEGER :: i,tout,jout
c-----------------------------------------------------------------------
c     Defaults.
c-----------------------------------------------------------------------
      IF(PRESENT(ti))THEN
         tout = ti
      ELSE
         tout = tmag_out
      ENDIF
      IF(PRESENT(ji))THEN
         jout = ji
      ELSE
         jout = 0
      ENDIF
      fmo = 0
c-----------------------------------------------------------------------
c     Transform vector in larger m space and then transfer to mo
c-----------------------------------------------------------------------
      IF(mpert<lmpert)THEN
         ! transfer to mo space
         DO i=1,mpert
           IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
              fmo(mlow-lmlow+i) = fmi(i)
           ENDIF
         ENDDO
         ! transform new vector
         IF((jac_out/=jac_type).OR.(tout==0).OR.(jout/=0))THEN
            CALL gpeq_bcoords(psi,fmo,lmfac,lmpert,power_rout,
     $         power_bpout,power_bout,power_rcout,tout,jout)
         ENDIF
      ELSE
         ! transform in mi space
         tmp = fmi
         IF((jac_out/=jac_type).OR.(tout==0).OR.(jout/=0))THEN
            CALL gpeq_bcoords(psi,tmp,mfac,mpert,power_rout,
     $         power_bpout,power_bout,power_rcout,tout,jout)
         ENDIF
         ! transfer to mo space
         DO i=1,mpert
           IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
              fmo(mlow-lmlow+i) = tmp(i)
           ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_bcoordsout
c-----------------------------------------------------------------------
c     subprogram 17. gpeq_weight.
c     switch between a function and a weighted function.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_weight(psi,ftnmn,amf,amp,wegt)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,wegt
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      LOGICAL :: first = .TRUE.
      INTEGER :: itheta
      REAL(r8) :: psave = 0

      REAL(r8), DIMENSION(:), ALLOCATABLE :: delpsi,wgtfun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      ! note automatic arrays are allocated and deallocated on entry/exit
      ! instead, we use allocatables and just allocate once for all
      SAVE :: psave,wgtfun,delpsi

      IF(debug_flag) PRINT *, "Entering gpeq_weight"

      ! form expensive geometry factors if called on new surface
      IF(psave/=psi) THEN
         IF(first) ALLOCATE(delpsi(0:mthsurf),wgtfun(0:mthsurf))
         first = .FALSE.
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            wgtfun(itheta)=1.0/(jac*delpsi(itheta))
         ENDDO
      ENDIF
      ! convert from fourier to real space
      CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     weight function.
c-----------------------------------------------------------------------
      SELECT CASE(wegt)
      CASE(0) ! flux to field
         ftnfun=ftnfun*wgtfun
      CASE(1) ! field to flux
         ftnfun=ftnfun/wgtfun
      CASE(2) ! field to sqrt(A)b
         ftnfun=ftnfun/sqrt(wgtfun)
      CASE(3)
         ftnfun=ftnfun/sqrt(wgtfun*r)
      CASE(4) ! x^psi to x^norm
         ftnfun=ftnfun/delpsi
      CASE(5) ! x^norm to x^psi
         ftnfun=ftnfun*delpsi
      CASE(6) ! flux to sqrt(A)b
         ftnfun=ftnfun*sqrt(wgtfun)
      END SELECT
      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_weight"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_weight
c-----------------------------------------------------------------------
c     subprogram 18. gpeq_rzpgrid.
c     find magnetic coordinates for given rz coords.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_rzpgrid(nr,nz,psixy)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz,psixy

      INTEGER :: i,j,itheta
      REAL(r8) :: xint,zint,ttheta,ptheta

      REAL(r8), DIMENSION(0:mthsurf) :: thetas,rbarr,etarr

      TYPE(spline_type) :: rbeta

      ALLOCATE(gdr(0:nr,0:nz),gdz(0:nr,0:nz),gdl(0:nr,0:nz),
     $     gdpsi(0:nr,0:nz),gdthe(0:nr,0:nz),gdphi(0:nr,0:nz))

      IF(debug_flag) WRITE(*,*) "Entering gpeq_rzgrid"
      IF(debug_flag) WRITE(*,*) "  ",nr,nz
c-----------------------------------------------------------------------
c     invert given rzphi to magnetic coordinates.
c-----------------------------------------------------------------------
      gdr=0
      gdz=0
      gdl=0
      gdpsi=0
      gdthe=0
      gdphi=0

      IF(psixy<1)THEN
         DO i=0,nr
            DO j=0,nz
               gdr(i,j) = rmin+i*(rmax-rmin)/nr
               gdz(i,j) = -zlim+j*2.0*zlim/nz
            ENDDO
         ENDDO
         RETURN
      ENDIF
      xint=(psi_in%xs(mr)-psi_in%xs(0))/nr
      zint=(psi_in%ys(mz)-psi_in%ys(0))/nz
      IF(debug_flag) WRITE(*,*) "Used psi_in"

      ! information of plasma boundary
      CALL spline_alloc(rbeta,mthsurf,1)
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psilim,theta(itheta),0)
         rbarr(itheta)=SQRT(rzphi%f(1))
         etarr(itheta)=theta(itheta)+rzphi%f(2)
      ENDDO
      rbeta%xs=etarr
      rbeta%fs(:,1)=rbarr
      CALL spline_fit(rbeta,"periodic")

      DO i=0,nr 
         DO j=0,nz
            gdr(i,j)=psi_in%xs(0)+i*xint
            gdz(i,j)=psi_in%ys(0)+j*zint

            ! compare grid with equilibrium input
            IF ((nr==mr) .AND. (nz==mz)) THEN
               gdpsi(i,j)=psi_in%fs(i,j,1)
            ELSE
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),0)
               gdpsi(i,j)=psi_in%f(1)
            ENDIF
            ! avoid o-point 
            IF (gdpsi(i,j)<psilow) gdpsi(i,j)=psilow
            ttheta=ATAN2((gdz(i,j)-zo),(gdr(i,j)-ro))
            IF (ttheta >= 0) THEN 
               ptheta=ttheta/twopi
            ELSE
               ptheta=1+ttheta/twopi
            ENDIF

            IF (gdpsi(i,j)<psilim) THEN
               CALL spline_eval(rbeta,ptheta,0)
            ! recheck whether it is inside the boundary 
               IF (SQRT((gdr(i,j)-ro)**2+(gdz(i,j)-zo)**2)<rbeta%f(1))
     $              THEN
                  gdl(i,j)=1
                  DO itheta=0,mthsurf
                     CALL bicube_eval(rzphi,gdpsi(i,j),theta(itheta),0)
                     thetas(itheta)=(theta(itheta)+rzphi%f(2))
                  ENDDO
                  gdthe(i,j)=issect(mthsurf,theta(:),thetas(:),ptheta)
                  CALL bicube_eval(rzphi,gdpsi(i,j),gdthe(i,j),0)
                  gdphi(i,j)=-rzphi%f(3)/twopi
               ENDIF
            ENDIF
            ! mark inside the domain of calculations
            IF (gdpsi(i,j)<psifac(1)) gdl(i,j)=2
         ENDDO
      ENDDO
      CALL spline_dealloc(rbeta)
      IF(debug_flag) WRITE(*,*) "Leaving gpeq_rzgrid"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_rzpgrid
c-----------------------------------------------------------------------
c     subprogram 19. gpeq_rzpdiv.
c     make zero divergence of rzphi functions.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_rzpdiv(nr,nz,rval,zval,fr,fz,fp)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz
      REAL(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: rval,zval
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: fr,fz
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: fp
 
      INTEGER :: i,j

      TYPE(bicube_type) :: rfr,ifr,rfz,ifz

      IF(verbose) WRITE(*,*)"Modifying bphi to make zero divergence"

      CALL bicube_alloc(rfr,nr,nz,1)
      CALL bicube_alloc(ifr,nr,nz,1)
      CALL bicube_alloc(rfz,nr,nz,1)
      CALL bicube_alloc(ifz,nr,nz,1)

      rfr%xs=rval(:,0)
      ifr%xs=rval(:,0)
      rfz%xs=rval(:,0)
      ifz%xs=rval(:,0) 

      rfr%ys=zval(0,:)
      ifr%ys=zval(0,:)
      rfz%ys=zval(0,:)
      ifz%ys=zval(0,:) 

      rfr%fs(:,:,1)=rval*REAL(fr)
      ifr%fs(:,:,1)=rval*AIMAG(fr)
      rfz%fs(:,:,1)=REAL(fz)
      ifz%fs(:,:,1)=AIMAG(fz)

      CALL bicube_fit(rfr,"extrap","extrap")
      CALL bicube_fit(ifr,"extrap","extrap")
      CALL bicube_fit(rfz,"extrap","extrap")
      CALL bicube_fit(ifz,"extrap","extrap")      

      DO i=0,nr
         DO j=0,nz
            CALL bicube_eval(rfr,rval(i,j),zval(i,j),1)
            CALL bicube_eval(ifr,rval(i,j),zval(i,j),1)
            CALL bicube_eval(rfz,rval(i,j),zval(i,j),1)
            CALL bicube_eval(ifz,rval(i,j),zval(i,j),1)

            fp(i,j)=((ifr%fx(1)+ifz%fy(1)*rval(i,j))-
     $          ifac*(rfr%fx(1)+rfz%fy(1)*rval(i,j)))/nn

         ENDDO
      ENDDO
       
      CALL bicube_dealloc(rfr)
      CALL bicube_dealloc(ifr)
      CALL bicube_dealloc(rfz)
      CALL bicube_dealloc(ifz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_rzpdiv
c-----------------------------------------------------------------------
c     subprogram 20. gpeq_alloc.
c     allocate essential vectors in fourier space 
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_alloc
      IF(debug_flag) PRINT *, "Entering gpeq_alloc"

      ALLOCATE(xsp_mn(mpert),xsp1_mn(mpert),xss_mn(mpert),xms_mn(mpert),
     $     xwp_mn(mpert),xwt_mn(mpert),xwz_mn(mpert),xmt_mn(mpert),
     $     bwp_mn(mpert),bwt_mn(mpert),bwz_mn(mpert),bmt_mn(mpert),
     $     bwp1_mn(mpert),xmp1_mn(mpert),
     $     xvp_mn(mpert),xvt_mn(mpert),xvz_mn(mpert),xmz_mn(mpert),
     $     bvp_mn(mpert),bvt_mn(mpert),bvz_mn(mpert),bmz_mn(mpert),
     $     xno_mn(mpert),xta_mn(mpert),xpa_mn(mpert),
     $     bno_mn(mpert),bta_mn(mpert),bpa_mn(mpert),
     $     xrr_mn(mpert),xrz_mn(mpert),xrp_mn(mpert),
     $     brr_mn(mpert),brz_mn(mpert),brp_mn(mpert),
     $     c2vp_mn(mpert),c2vt_mn(mpert),c2vz_mn(mpert),
     $     cvp_mn(mpert),cvt_mn(mpert),cvz_mn(mpert),
     $     cwp_mn(mpert),cwt_mn(mpert),cwz_mn(mpert))
      IF(debug_flag) PRINT *, "->Leaving gpeq_alloc"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_alloc
c-----------------------------------------------------------------------
c     subprogram 21. gpeq_dealloc.
c     deallocate essential vectors in fourier space 
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_dealloc
      IF(debug_flag) PRINT *, "Entering gpeq_dealloc"

      DEALLOCATE(xsp_mn,xsp1_mn,xss_mn,xms_mn,bwp1_mn,xmp1_mn,
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,xmt_mn,bmt_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,xmz_mn,bmz_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn,
     $     c2vp_mn,c2vt_mn,c2vz_mn,cvp_mn,cvt_mn,cvz_mn,
     $     cwp_mn,cwt_mn,cwz_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_dealloc"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_dealloc
c-----------------------------------------------------------------------
c     subprogram 22. gpeq_interp_singsurf.
c     create spline for interpretation of solution near singular surface.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_interp_singsurf(fsp_sol,spot,npsi)
      TYPE(cspline_type), INTENT(INOUT)::fsp_sol    ! spline of bwp smoothly crossing rationals
      REAL(r8), INTENT(IN) :: spot                  ! roughly the span in  m-nq to cross
      INTEGER, INTENT(IN) :: npsi                   ! number of points between rationals in the spline

      INTEGER::psisize,ising,ix,icount
      INTEGER,PARAMETER:: method=1
      REAL(r8)::nq1,x,x0,x1
      REAL(r8),DIMENSION(msing)::respsi,dxl,dxr
c-----------------------------------------------------------------------
c     detemine spline allocation.
c-----------------------------------------------------------------------
      DO ising=1,msing
         respsi(ising)=singtype(ising)%psifac
      ENDDO
      DO ising=1,msing
         nq1=singtype(ising)%q1*nn
         SELECT CASE(method)
         CASE(1)
         IF (ising==1) THEN
            dxl(ising)=respsi(ising)
     $                 -spot*(respsi(ising)-psilow)
            dxr(ising)=respsi(ising)
     $                 +spot*(respsi(ising+1)-respsi(ising))
         ELSEIF (ising==msing) THEN
            dxl(ising)=respsi(ising)
     $                 -spot*(respsi(ising)-respsi(ising-1))
            dxr(ising)=respsi(ising)
     $                 +spot*(psilim-respsi(ising))
         ELSE
            dxl(ising)=respsi(ising)
     $                 -spot*(respsi(ising)-respsi(ising-1))
            dxr(ising)=respsi(ising)
     $                 +spot*(respsi(ising+1)-respsi(ising))
         ENDIF
         CASE(2)
            dxl(ising)=respsi(ising)-spot/nq1
            dxr(ising)=respsi(ising)+spot/nq1
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     construct solution spline.
c-----------------------------------------------------------------------
      icount=npsi*(msing+1)-1
      CALL cspline_alloc(fsp_sol,icount,mpert)
      icount=0
      DO ising=1,msing+1
         IF (ising==1) THEN
            x0=psilow
            x1=dxl(ising)
         ELSEIF (ising==msing+1) THEN
            x0=dxr(ising-1)
            x1=psilim
         ELSE
            x0=dxr(ising-1)
            x1=dxl(ising)
         ENDIF
         DO ix=1,npsi
            x=x0+(ix-1)*(x1-x0)/(npsi-1)
            CALL gpeq_sol(x)
            fsp_sol%xs(icount)=x
            fsp_sol%fs(icount,:)=bwp_mn
            icount=icount+1
         ENDDO
      ENDDO
      CALL cspline_fit(fsp_sol,"not-a-knot")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_interp_singsurf
c-----------------------------------------------------------------------
c     subprogram 23. gpeq_interp_sol.
c     get bwn at psi after calling gpeq_interp_singsurf.
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_interp_sol(fsp_sol,psi,interpbwn)
      TYPE(cspline_type), INTENT(INOUT)::fsp_sol
      REAL(r8), INTENT(IN):: psi
      COMPLEX(r8),DIMENSION(mpert), INTENT(OUT)::interpbwn

      CALL cspline_eval(fsp_sol,psi,0)
      interpbwn=fsp_sol%f
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_interp_sol

      END MODULE gpeq_mod
