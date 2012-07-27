c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     IPEQ: calculate functions for perturbed equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipeq_mod
c      1. ipeq_sol
c      2. ipeq_contra
c      3. ipeq_cova
c      4. ipeq_normal
c      5. ipeq_tangent
c      6. ipeq_parallel
c      7. ipeq_rzphi
c      8. ipeq_surface
c      9. ipeq_fcoords
c     10. ipeq_bcoords
c     11. ipeq_bntoxp
c     12. ipeq_xptobn
c     13. ipeq_weight
c     14. ipeq_rzpgrid
c     15. ipeq_rzpdiv
c     16. ipeq_alloc
c     17. ipeq_dealloc
c-----------------------------------------------------------------------
c     subprogram 0. ipeq_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipeq_mod
      USE idcon_mod
 
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipeq_sol.
c     obtain solutions of perturbed quantities.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_sol(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER :: ipert,jpert,iqty

      INTEGER, DIMENSION(mpert) :: ipiva
      COMPLEX(r8), DIMENSION(mpert*mpert) :: work

      COMPLEX(r8), DIMENSION(mpert) :: xspfac
c-----------------------------------------------------------------------
c     evaluate matrices and solutions.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,1)
      CALL cspline_eval(u1,psi,0)
      CALL cspline_eval(u2,psi,0)

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
      xsp_mn=u1%f
      xspfac=u2%f/singfac
      CALL zgbmv('N',mpert,mpert,mband,mband,-ione,kmats,
     $     2*mband+1,u1%f,1,ione,xspfac,1)
      CALL zpbtrs('L',mpert,mband,1,fmats,mband+1,xspfac,mpert,info)
      xsp1_mn=xspfac/singfac
      CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,bmat,mpert,info)
      CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,cmat,mpert,info)
      xss_mn=-MATMUL(bmat,xsp1_mn)-MATMUL(cmat,xsp_mn)
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
         xms_mn=xss_mn*(singfac**2/(singfac**2+reg_spot**2))
         bmt_mn=bwt_mn*(singfac**2/(singfac**2+reg_spot**2))
      ELSE
         xms_mn=xss_mn
         bmt_mn=bwt_mn
      ENDIF

      DEALLOCATE(amat,bmat,cmat,fmats,gmats,kmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_sol
c-----------------------------------------------------------------------
c     subprogram 2. ipeq_contra.
c     compute contravariant components of perturbed quantities.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_contra(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      
      INTEGER :: ipert,jpert,m1,m2,dm
      COMPLEX(r8), DIMENSION(-mband:mband) ::jmat,jmat1

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
c     compute contravariant quantities.
c-----------------------------------------------------------------------
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
     $           jmat1(dm)*xsp_mn(jpert)+
     $           twopi*ifac*nn/chi1*jmat(dm)*xss_mn(jpert))/
     $           (twopi*ifac*(m1-nn*q))
            xwz_mn(ipert)=xwz_mn(ipert)-(q*jmat(dm)*xsp1_mn(jpert)+
     $           q*jmat1(dm)*xsp_mn(jpert)+
     $           twopi*ifac*m1/chi1*jmat(dm)*xss_mn(jpert))/
     $           (twopi*ifac*(m1-nn*q))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute modified quantities.
c-----------------------------------------------------------------------
      IF (reg_flag) THEN
         xmt_mn=xwt_mn*(singfac**2/(singfac**2+reg_spot**2))
         xmz_mn=xwz_mn*(singfac**2/(singfac**2+reg_spot**2))
      ELSE
         xmt_mn=xwt_mn
         xmz_mn=xwz_mn
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_contra
c-----------------------------------------------------------------------
c     subprogram 3. ipeq_cova.
c     compute covariant components.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_cova(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8) ,INTENT(IN) :: psi

      INTEGER :: ipert,jpert,m1,m2,dm
      COMPLEX(r8), DIMENSION(-mband:mband) :: g11,g22,g33,g23,g31,g12

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
     $           g12(dm)*bwt_mn(jpert)+g31(dm)*bwz_mn(jpert)
            bvt_mn(ipert)=bvt_mn(ipert)+g12(dm)*bwp_mn(jpert)+
     $           g22(dm)*bwt_mn(jpert)+g23(dm)*bwz_mn(jpert)
            bvz_mn(ipert)=bvz_mn(ipert)+g31(dm)*bwp_mn(jpert)+
     $           g23(dm)*bwt_mn(jpert)+g33(dm)*bwz_mn(jpert)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_cova
c-----------------------------------------------------------------------
c     subprogram 4. ipeq_normal.
c     compute normal components.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_normal(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xno_fun,bno_fun
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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_normal
c-----------------------------------------------------------------------
c     subprogram 5. ipeq_tangent.
c     compute tangential components.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_tangent(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqrpsi,
     $     psithe,psiphi,eqb,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,xwt_fun,
     $     bwt_fun,xwz_fun,bwz_fun,xta_fun,bta_fun
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      q=sq%f(4)
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
         w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
         w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
         w(3,1)=rzphi%fx(3)*2*rfac
         w(3,2)=rzphi%fy(3)/(twopi*rfac)
         sqrpsi(itheta)=w(1,1)**2+w(1,2)**2
         delpsi(itheta)=SQRT(sqrpsi(itheta))
         psithe(itheta)=w(1,1)*w(2,1)+w(1,2)*w(2,2)
         psiphi(itheta)=w(1,1)*w(3,1)+w(1,2)*w(3,2)
         eqb(itheta)=SQRT((sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $           /(twopi*rfac)**2)
      ENDDO
c-----------------------------------------------------------------------
c     compute tangential components, b times delpsi.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
      CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
      CALL iscdftb(mfac,mpert,xwz_fun,mthsurf,xwz_mn)
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
      CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
      CALL iscdftb(mfac,mpert,bwz_fun,mthsurf,bwz_mn)
      xta_fun=chi1*((q*xwt_fun-xwz_fun)*sqrpsi-
     $     xwp_fun*(q*psithe-psiphi))/(eqb*jacs*delpsi)
      bta_fun=chi1*((q*bwt_fun-bwz_fun)*sqrpsi-
     $     bwp_fun*(q*psithe-psiphi))/(eqb*jacs*delpsi)
      CALL iscdftf(mfac,mpert,xta_fun,mthsurf,xta_mn)
      CALL iscdftf(mfac,mpert,bta_fun,mthsurf,bta_mn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_tangent
c-----------------------------------------------------------------------
c     subprogram 6. ipeq_parallel.
c     compute parallel components for xi and b.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_parallel(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: eqb
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xvt_fun,bvt_fun,
     $     xvz_fun,bvz_fun,xpa_fun,bpa_fun
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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_parallel
c-----------------------------------------------------------------------
c     subprogram 7. ipeq_rzphi.
c     compute rzphi components.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_rzphi(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: t11,t12,t21,t22,t33,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun,xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun
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
      CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
      CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_rzphi
c-----------------------------------------------------------------------
c     subprogram 8. ipeq_surface.
c     compute surface currents and potentials.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_surface(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER :: i,j,itheta,rtheta

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: rbwp_fun,bwp_fun,
     $     chi_fun,che_fun,kax_fun,xwp_fun,xwt_fun,
     $     bvt_fun,bvz_fun
      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
      COMPLEX(r8), DIMENSION(4,0:mthsurf) :: chp_fun,kap_fun

      REAL(r8), DIMENSION(:), POINTER :: grri_real,grri_imag,
     $     grre_real,grre_imag

      ALLOCATE(grri_real(nths2),grri_imag(nths2),
     $     grre_real(nths2),grre_imag(nths2))
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
c     reverse poloidal and torodial coordinates.
c-----------------------------------------------------------------------
      DO itheta=0,mthvac
         rtheta=mthvac-itheta
         rbwp_fun(itheta)=CONJG(bwp_fun(rtheta))
      ENDDO
      CALL iscdftf(mfac,mpert,rbwp_fun,mthvac,rbwp_mn)
c-----------------------------------------------------------------------
c     compute vacuum magnetic potentials with reverse normal vector.
c-----------------------------------------------------------------------
      grri_real=MATMUL(grri,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      grri_imag=MATMUL(grri,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
      grre_real=MATMUL(grre,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
      grre_imag=MATMUL(grre,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
c-----------------------------------------------------------------------
c     return into original coordinates.
c-----------------------------------------------------------------------
      DO itheta=0,mthvac-1
         rtheta=mthvac-itheta
         chi_fun(itheta+1)=(grri_real(rtheta)-ifac*grri_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
         che_fun(itheta+1)=(grre_real(rtheta)-ifac*grre_imag(rtheta))
     $        *EXP(-ifac*nn*dphi(itheta+1))
      ENDDO
      chi_fun(0)=chi_fun(mthvac)
      che_fun(0)=che_fun(mthvac)
c-----------------------------------------------------------------------
c     normalize chi functions of vacuum.
c-----------------------------------------------------------------------
      chi_fun=chi_fun/twopi**2
      che_fun=-che_fun/twopi**2
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
      CALL iscdftb(mfac,mpert,chp_fun(2,:),mthvac,chp_mn(2,:))
      CALL iscdftb(mfac,mpert,chp_fun(4,:),mthvac,chp_mn(4,:))
      DO j=1,4
         kap_fun(j,:)=(chp_fun(j,:)-che_fun(:))/mu0
         CALL iscdftf(mfac,mpert,kap_fun(j,:),mthsurf,kap_mn(j,:))
      ENDDO
      kax_fun=(chi_fun-che_fun)/mu0
      CALL iscdftf(mfac,mpert,kax_fun,mthsurf,kax_mn)
     
      DEALLOCATE(grri_real,grri_imag,grre_real,grre_imag)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_surface
c-----------------------------------------------------------------------
c     subprogram 9. ipeq_fcoords.
c     transform coordinates to dcon coordinates. 
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_fcoords(psi,ftnmn,amf,amp,ri,bpi,bi,rci,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,ri,bpi,bi,rci,ti,ji
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      INTEGER :: i,ising,itheta
      REAL(r8) :: thetai,jarea

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacfac
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      TYPE(spline_type) :: spl       
      
      CALL spline_alloc(spl,mthsurf,2)
      spl%xs=theta

      CALL spline_eval(sq,psi,0)
      dphi=0
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
     $        bpfac**bpi*bfac**bi
         ! jacobian for coordinate angle at dcon angle
         spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**ri*rfac**rci/
     $        (bpfac**bpi*bfac**bi)  
         IF (ti .EQ. 0) THEN
            dphi(itheta)=rzphi%f(3)
         ENDIF
      ENDDO      

      CALL spline_fit(spl,"periodic")
      CALL spline_int(spl)

      ! coordinate angle at dcon angle
      thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
      IF (ji .EQ. 1) THEN
         DO itheta=0,mthsurf
            ! jacobian at coordinate angle
            thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
            CALL spline_eval(spl,thetai,0)
            jacfac(itheta)=spl%f(2)
         ENDDO
         jarea=0
         DO itheta=0,mthsurf-1
            jarea=jarea+jacfac(itheta)/mthsurf
         ENDDO
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

      CALL spline_dealloc(spl)

      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_fcoords
c-----------------------------------------------------------------------
c     subprogram 10. ipeq_bcoords.
c     transform dcon coordinates to other coordinates.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_bcoords(psi,ftnmn,amf,amp,ri,bpi,bi,rci,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,ri,bpi,bi,rci,ti,ji
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      INTEGER :: i,ising,itheta
      REAL(r8) :: thetai,jarea

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacfac
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      TYPE(spline_type) :: spl       
      
      CALL spline_alloc(spl,mthsurf,2)
      spl%xs=theta

      CALL spline_eval(sq,psi,0)
      dphi=0
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
     $        bpfac**bpi*bfac**bi
         spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**ri*rfac**rci/
     $        (bpfac**bpi*bfac**bi)
         IF (ti .EQ. 0) THEN
            dphi(itheta)=rzphi%f(3)
         ENDIF
      ENDDO

      CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)

      CALL spline_fit(spl,"periodic")
      CALL spline_int(spl)

      ! coordinate angle at dcon angle
      thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
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
     $           ftnmn(i)*EXP(ifac*twopi*(i+amf(1)-1)*thetai)
         ENDDO

         ! take toroidal factor back for coordinate angle
         IF (ti .NE. 0) THEN
            ftnfun(itheta)=ftnfun(itheta)*
     $           EXP(-twopi*ifac*nn*sq%f(4)*(thetai-theta(itheta)))
         ENDIF

         ! jacobian at coordinate angle
         CALL spline_eval(spl,thetai,0)
         jacfac(itheta)=spl%f(2)
      ENDDO
         
      IF (ji .EQ. 1) THEN
         jarea=0
         DO itheta=0,mthsurf-1
            jarea=jarea+jacfac(itheta)/mthsurf
         ENDDO
         ftnfun=ftnfun*jacfac/jarea
      ENDIF
   
      CALL spline_dealloc(spl)
      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_bcoords
c-----------------------------------------------------------------------
c     subprogram 11. ipeq_bntoxp.
c     transform b normal to xi normal solution.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_bntoxp(psi,finmn,foutmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: finmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: foutmn

      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs

      COMPLEX(r8), DIMENSION(0:mthsurf) :: finfun,foutfun

      CALL spline_eval(sq,psi,0)
      singfac=mfac-nn*sq%f(4)
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

      CALL iscdftb(mfac,mpert,finfun,mthsurf,finmn)      
      foutfun=finfun*jacs*delpsi
      CALL iscdftf(mfac,mpert,foutfun,mthsurf,foutmn)
  
      foutmn=foutmn/(chi1*singfac*twopi*ifac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_bntoxp
c-----------------------------------------------------------------------
c     subprogram 12. ipeq_xptobn.
c     transforms xi normal solution to b normal.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_xptobn(psi,finmn,foutmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: finmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: foutmn      

      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs

      COMPLEX(r8), DIMENSION(0:mthsurf) :: finfun,foutfun

      CALL spline_eval(sq,psi,0)
      singfac=mfac-nn*sq%f(4)
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

      foutmn=finmn*chi1*singfac*twopi*ifac

      CALL iscdftb(mfac,mpert,foutfun,mthsurf,foutmn)      
      foutfun=foutfun/(jacs*delpsi)
      CALL iscdftf(mfac,mpert,foutfun,mthsurf,foutmn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_xptobn
c-----------------------------------------------------------------------
c     subprogram 13. ipeq_weight.
c     switch between a function and a weighted function.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_weight(psi,ftnmn,amf,amp,wegt)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,wegt
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,wgtfun
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

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
      CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     weight function.
c-----------------------------------------------------------------------
      SELECT CASE(wegt)
      CASE(0)
         ftnfun=ftnfun*wgtfun
      CASE(1)
         ftnfun=ftnfun/wgtfun
      CASE(2)
         ftnfun=ftnfun/sqrt(wgtfun)
      CASE(3)
         ftnfun=ftnfun/sqrt(wgtfun*r)
      CASE(4)
         ftnfun=ftnfun/delpsi
      CASE(5)
         ftnfun=ftnfun*delpsi
      END SELECT
      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_weight
c-----------------------------------------------------------------------
c     subprogram 14. ipeq_rzpgrid.
c     find magnetic coordinates for given rz coords.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_rzpgrid(nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz

      INTEGER :: i,j,itheta
      REAL(r8) :: xint,zint,ttheta,ptheta

      REAL(r8), DIMENSION(0:mthsurf) :: thetas,rbarr,etarr

      TYPE(spline_type) :: rbeta

      ALLOCATE(gdr(0:nr,0:nz),gdz(0:nr,0:nz),gdl(0:nr,0:nz),
     $     gdpsi(0:nr,0:nz),gdthe(0:nr,0:nz),gdphi(0:nr,0:nz))
c-----------------------------------------------------------------------
c     invert given rzphi to magnetic coordinates.
c-----------------------------------------------------------------------
      gdr=0
      gdz=0
      gdl=0
      gdpsi=0
      gdthe=0
      gdphi=0

      xint=(psi_in%xs(mr)-psi_in%xs(0))/nr
      zint=(psi_in%ys(mz)-psi_in%ys(0))/nz

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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_rzpgrid
c-----------------------------------------------------------------------
c     subprogram 15. ipeq_rzpdiv.
c     make zero divergence of rzphi functions.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_rzpdiv(nr,nz,lval,rval,zval,fr,fz,fp)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz
      INTEGER, DIMENSION(0:nr,0:nz), INTENT(IN) :: lval
      REAL(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: rval,zval
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: fr,fz
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: fp
 
      INTEGER :: i,j

      TYPE(bicube_type) :: rfr,ifr,rfz,ifz

      WRITE(*,*)"Modifying bphi to make zero divergence"

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
      END SUBROUTINE ipeq_rzpdiv     
c-----------------------------------------------------------------------
c     subprogram 16. ipeq_alloc.
c     allocate essential vectors in fourier space 
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_alloc

      ALLOCATE(xsp_mn(mpert),xsp1_mn(mpert),xss_mn(mpert),xms_mn(mpert),
     $     xwp_mn(mpert),xwt_mn(mpert),xwz_mn(mpert),xmt_mn(mpert),
     $     bwp_mn(mpert),bwt_mn(mpert),bwz_mn(mpert),bmt_mn(mpert),
     $     bwp1_mn(mpert),
     $     xvp_mn(mpert),xvt_mn(mpert),xvz_mn(mpert),xmz_mn(mpert),
     $     bvp_mn(mpert),bvt_mn(mpert),bvz_mn(mpert),bmz_mn(mpert),
     $     xno_mn(mpert),xta_mn(mpert),xpa_mn(mpert),
     $     bno_mn(mpert),bta_mn(mpert),bpa_mn(mpert),
     $     xrr_mn(mpert),xrz_mn(mpert),xrp_mn(mpert),
     $     brr_mn(mpert),brz_mn(mpert),brp_mn(mpert))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_alloc
c-----------------------------------------------------------------------
c     subprogram 17. ipeq_dealloc.
c     deallocate essential vectors in fourier space 
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_dealloc

      DEALLOCATE(xsp_mn,xsp1_mn,xss_mn,xms_mn,bwp1_mn,
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,xmt_mn,bmt_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,xmz_mn,bmz_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_dealloc

      END MODULE ipeq_mod
