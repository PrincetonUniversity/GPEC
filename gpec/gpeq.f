c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     GPEQ: calculate functions for perturbed equilibrium.
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
c      9. gpeq_fcoords
c     10. gpeq_fcoordsout
c     11. gpeq_bcoords
c     12. gpeq_bcoordsout
c     13. gpeq_weight
c     14. gpeq_rzpgrid
c     15. gpeq_rzpdiv
c     16. gpeq_alloc
c     17. gpeq_dealloc
c     18. gpeq_interp_singsurf
c     19. gpeq_interp_sol
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
            CALL cspline_eval(u2,psi,0)
            xsp1_mn=u1%f1
         ELSE
            CALL cspline_eval(u1,psi,0)
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
c     compute contravariant b fields.
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
c     compute covariant components.
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
c     compute covariant components with metric tensors.
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
c-----------------------------------------------------------------------
c     subprogram 3-2. gpeq_cur
c     compute perturbed current components.
c     gpeq_sol, gpeq_contra and gpeq_cova should be called prior to this
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_cur(psi)

      REAL(r8), INTENT(IN) :: psi

      INTEGER :: istep,itheta,ipert,m1,dm,jpert,loc(1)
      REAL(r8), DIMENSION(0:mthsurf) :: eqb
      COMPLEX(r8), DIMENSION(mpert) :: bvtl_mn,bvzl_mn,bvt1_mn,bvz1_mn
      COMPLEX(r8), DIMENSION(-mband:mband) :: g11,g22,g33,g23,g31,g12
      COMPLEX(r8), DIMENSION(0:mthsurf) :: jvt_fun, jvz_fun, jpa_fun

      IF (debug_flag) PRINT * ,"Entering gpeq_cur"
c-----------------------------------------------------------------------
c     compute contravarient component.
c-----------------------------------------------------------------------
      ! psi derivative of covarent field
      loc = MINLOC(ABS(psifac(1:mstep-1) - psi))
      istep = loc(1)
      CALL gpeq_sol(psifac(istep-1))
      CALL gpeq_contra(psifac(istep-1))
      CALL gpeq_cova(psifac(istep-1))
      bvtl_mn = bvt_mn
      bvzl_mn = bvz_mn
      CALL gpeq_sol(psifac(istep+1))
      CALL gpeq_contra(psifac(istep+1))
      CALL gpeq_cova(psifac(istep+1))
      bvt1_mn = (bvt_mn - bvtl_mn) / (psifac(istep+1) - psifac(istep-1))
      bvz1_mn = (bvz_mn - bvzl_mn) / (psifac(istep+1) - psifac(istep-1))
      CALL gpeq_sol(psi)
      CALL gpeq_contra(psi)
      CALL gpeq_cova(psi)
      ! contravarient current
      jwp_mn=twopi*ifac*mfac*bvz_mn-(-twopi*ifac*nn*bvt_mn)
      jwt_mn=twopi*ifac*nn*bvp_mn-(bvz1_mn)
      jwz_mn=bvt1_mn-(twopi*ifac*mfac*bvp_mn)
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
      CALL cspline_eval(metric%cs,psi,0)
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
c     compute covariant components with metric tensors.
c-----------------------------------------------------------------------
      ipert=0
      jvp_mn=0
      jvt_mn=0
      jvz_mn=0 
      DO m1=mlow,mhigh
         ipert=ipert+1
         DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
            jpert=ipert+dm
            jvp_mn(ipert)=jvp_mn(ipert)+g11(dm)*jwp_mn(jpert)+
     $           g12(dm)*jwt_mn(jpert)+g31(dm)*jwz_mn(jpert)
            jvt_mn(ipert)=jvt_mn(ipert)+g12(dm)*jwp_mn(jpert)+
     $           g22(dm)*jwt_mn(jpert)+g23(dm)*jwz_mn(jpert)
            jvz_mn(ipert)=jvz_mn(ipert)+g31(dm)*jwp_mn(jpert)+
     $           g23(dm)*jwt_mn(jpert)+g33(dm)*jwz_mn(jpert)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute parallel component.
c-----------------------------------------------------------------------
      IF (psifac(istep) == psi) THEN
          eqb(:)=eqfun%fs(istep,:,1)
          q = sq%fs(istep,4)
      ELSE
          CALL spline_eval(sq,psi,0)
          q=sq%f(4)
          DO itheta=0,mthsurf
             CALL bicube_eval(eqfun,psi,theta(itheta),0)
             eqb(itheta)=eqfun%f(1)
          ENDDO
      ENDIF
      CALL iscdftb(mfac,mpert,jvt_fun,mthsurf,jvt_mn)
      CALL iscdftb(mfac,mpert,jvz_fun,mthsurf,jvz_mn)
      jpa_fun=(jvt_fun+q*jvz_fun)/eqb!*chi1
      CALL iscdftf(mfac,mpert,jpa_fun,mthsurf,jpa_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_cur"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_cur
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
     $     xwp_fun,xwt_fun,bvt_fun,bvz_fun,bwp_fun
      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
      COMPLEX(r8), DIMENSION(4,0:mthsurf) :: chp_fun,kap_fun

      REAL(r8), DIMENSION(:), POINTER :: grri_real,grri_imag,
     $     grre_real,grre_imag

      ALLOCATE(grri_real(nths2),grri_imag(nths2),
     $     grre_real(nths2),grre_imag(nths2))
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
c-----------------------------------------------------------------------
c     return into original coordinates.
c-----------------------------------------------------------------------
      DO itheta=0,mthsurf-1
         rtheta=mthsurf-itheta
         chi_fun(itheta+1)=(grri_real(rtheta)-ifac*grri_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
         che_fun(itheta+1)=(grre_real(rtheta)-ifac*grre_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
      ENDDO
      chi_fun(0)=chi_fun(mthsurf)
      che_fun(0)=che_fun(mthsurf)
c-----------------------------------------------------------------------
c     normalize chi functions of vacuum.
c-----------------------------------------------------------------------
      chi_fun=chi_fun/(twopi**2)
      che_fun=-che_fun/(twopi**2)
      CALL iscdftf(mfac,mpert,chi_fun,mthsurf,chi_mn)
      CALL iscdftf(mfac,mpert,che_fun,mthsurf,che_mn)
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
      CALL iscdftf(mfac,mpert,kax_fun,mthsurf,kax_mn)
     
      DEALLOCATE(grri_real,grri_imag,grre_real,grre_imag)
      IF(debug_flag) PRINT *, "->Leaving gpeq_surface"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_surface
c-----------------------------------------------------------------------
c     subprogram 9. gpeq_fcoords.
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
c     subprogram 10. gpeq_fcoordsout.
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
c     subprogram 11. gpeq_bcoords.
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
c     subprogram 12. gpeq_bcoordsout.
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
c     subprogram 13. gpeq_weight.
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
c     subprogram 14. gpeq_rzpgrid.
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
c     subprogram 15. gpeq_rzpdiv.
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
c     subprogram 16. gpeq_alloc.
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
     $     jwp_mn(mpert),jwt_mn(mpert),jwz_mn(mpert),
     $     jvp_mn(mpert),jvt_mn(mpert),jvz_mn(mpert),jpa_mn(mpert))
      IF(debug_flag) PRINT *, "->Leaving gpeq_alloc"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_alloc
c-----------------------------------------------------------------------
c     subprogram 17. gpeq_dealloc.
c     deallocate essential vectors in fourier space 
c-----------------------------------------------------------------------
      SUBROUTINE gpeq_dealloc
      IF(debug_flag) PRINT *, "Entering gpeq_dealloc"

      DEALLOCATE(xsp_mn,xsp1_mn,xss_mn,xms_mn,bwp1_mn,xmp1_mn,
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,xmt_mn,bmt_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,xmz_mn,bmz_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn,
     $     jwp_mn,jwt_mn,jwz_mn,jvp_mn,jvt_mn,jvz_mn,jpa_mn)
      IF(debug_flag) PRINT *, "->Leaving gpeq_dealloc"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpeq_dealloc
c-----------------------------------------------------------------------
c     subprogram 18. gpeq_interp_singsurf.
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
c     subprogram 19. gpeq_interp_sol.
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
