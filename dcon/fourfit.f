c-----------------------------------------------------------------------
c     file fourfit.f
c     fits equilibrium quantities to Fourier series, evaluates matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fourfit_mod.
c     1. fourfit_make_metric.
c     2. fourfit_make_matrix.
c     3. fourfit_action_matrix.
c     4. fourfit_write_metric.
c     5. fourfit_write_matrix.
c     6. fourfit_evals.
c     7. fourfit_diagnose_1.
c-----------------------------------------------------------------------
c     subprogram 0. fourfit_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fourfit_mod
      USE fspline_mod
      USE dcon_mod
      USE torque, only : tpsi
      IMPLICIT NONE

      TYPE(fspline_type), PRIVATE :: metric,fmodb
      TYPE(cspline_type) :: dmats,emats,hmats,
     $     baats,caats,eaats,kaats,gaats,fmats,gmats,kmats,fbats,kbats
      TYPE(spline_type) :: k0s
      INTEGER, DIMENSION(:), POINTER :: ipiva
      COMPLEX(r8), DIMENSION(:,:), POINTER :: asmat,bsmat,csmat,
     $     fsmat,ksmat
      COMPLEX(r8), DIMENSION(:), POINTER :: jmat

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fourfit_make_metric.
c     computes fourier series of metric tensor components.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_make_metric

      INTEGER :: ipsi,itheta
      REAL(r8) :: theta,rfac,eta,r,jac,jac1,psifac,p1,q,
     $     g12,g22,g13,g23,g33,b2h,b2hp,b2ht,chi1
      REAL(r8), DIMENSION(3,3) :: v
c-----------------------------------------------------------------------
c     set up Fourier-spline type.
c-----------------------------------------------------------------------
      CALL fspline_alloc(metric,mpsi,mtheta,mband,8)
      metric%xs=rzphi%xs
      metric%ys=rzphi%ys*twopi
      metric%name="metric"
      metric%xtitle=" psi  "
      metric%ytitle="theta "
      metric%title=(/" g11  "," g22  "," g33  "," g23  "," g31  ",
     $     " g12  "," jmat ","jmat1 "/)
c-----------------------------------------------------------------------
c     set up kinetic Fourier-spline type.
c-----------------------------------------------------------------------
      CALL fspline_alloc(fmodb,mpsi,mtheta,mband,8)
      fmodb%xs=rzphi%xs
      fmodb%ys=rzphi%ys*twopi
      fmodb%name="fmodb"
      fmodb%xtitle=" psi  "
      fmodb%ytitle="theta "
      fmodb%title=(/" smat  "," tmat  "," xmat  ",
     $     " ymat1 "," ymat2 "," zmat1 ", " zmat2 "," zmat3 "/)
c-----------------------------------------------------------------------
c     begin loop over nodes.
c-----------------------------------------------------------------------
      chi1=twopi*psio      
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         p1=sq%fs1(ipsi,2)
         q=sq%fs(ipsi,4)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            CALL bicube_eval(eqfun,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            theta=rzphi%ys(itheta)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta+rzphi%f(2))
            r=ro+rfac*COS(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
            b2h=eqfun%f(1)**2/2
            b2hp=eqfun%f(1)*eqfun%fx(1)
            b2ht=eqfun%f(1)*eqfun%fy(1)  
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
            v(1,1)=rzphi%fx(1)/(2*rfac*jac)
            v(1,2)=rzphi%fx(2)*twopi*rfac/jac
            v(1,3)=rzphi%fx(3)*r/jac
            v(2,1)=rzphi%fy(1)/(2*rfac*jac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac/jac
            v(2,3)=rzphi%fy(3)*r/jac
            v(3,3)=twopi*r/jac
            g12=SUM(v(1,:)*v(2,:))*jac**2
            g13=v(3,3)*v(1,3)*jac**2
            g22=SUM(v(2,:)**2)*jac**2
            g23=v(2,3)*v(3,3)*jac**2
            g33=v(3,3)*v(3,3)*jac**2
c-----------------------------------------------------------------------
c     compute metric tensor components.
c-----------------------------------------------------------------------
            metric%fs(ipsi,itheta,1)=SUM(v(1,:)**2)*jac
            metric%fs(ipsi,itheta,2)=SUM(v(2,:)**2)*jac
            metric%fs(ipsi,itheta,3)=v(3,3)*v(3,3)*jac
            metric%fs(ipsi,itheta,4)=v(2,3)*v(3,3)*jac
            metric%fs(ipsi,itheta,5)=v(3,3)*v(1,3)*jac
            metric%fs(ipsi,itheta,6)=SUM(v(1,:)*v(2,:))*jac
            metric%fs(ipsi,itheta,7)=jac
            metric%fs(ipsi,itheta,8)=jac1
c-----------------------------------------------------------------------
c     compute kinetic metric tensor components.
c-----------------------------------------------------------------------
            fmodb%fs(ipsi,itheta,1)=jac*(p1+b2hp)
     $           -chi1**2*b2ht*(g12+q*g13)/(jac*b2h*2)
            fmodb%fs(ipsi,itheta,2)=
     $           chi1**2*b2ht*(g23+q*g33)/(jac*b2h*2)
            fmodb%fs(ipsi,itheta,3)=jac*b2h*2
            fmodb%fs(ipsi,itheta,4)=jac1*b2h*2-chi1**2*b2h*2*eqfun%fy(2)
            fmodb%fs(ipsi,itheta,5)=-twopi*chi1**2/jac*(g12+q*g13)
            fmodb%fs(ipsi,itheta,6)=chi1**2*b2h*2*eqfun%fy(3)
            fmodb%fs(ipsi,itheta,7)=twopi*chi1**2/jac*(g23+q*g33)
            fmodb%fs(ipsi,itheta,8)=twopi*chi1**2/jac*(g22+q*g23)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit Fourier-spline type.
c-----------------------------------------------------------------------
      IF(fft_flag)THEN
         CALL fspline_fit_2(metric,"extrap",.FALSE.)
         CALL fspline_fit_2(fmodb,"extrap",.FALSE.)
      ELSE
         CALL fspline_fit_1(metric,"extrap",.FALSE.)
         CALL fspline_fit_1(fmodb,"extrap",.FALSE.)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_make_metric
c-----------------------------------------------------------------------
c     subprogram 2. fourfit_action_matrix.
c     equilibrium matrices necessary to calc perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_action_matrix
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: ipsi,istep,ipert,jpert,itheta,dm,m1,m2
      REAL(r8) :: q,singfac2
      COMPLEX(r8), DIMENSION(-mband:mband) :: 
     $     sband,tband,xband,yband1,yband2,zband1,zband2,zband3
      COMPLEX(r8), DIMENSION(mpert,mpert) :: smat,tmat,xmat,ymat,zmat

      WRITE(*,*)"Computing action matrices"
c-----------------------------------------------------------------------
c     set up fourier-spline type.
c-----------------------------------------------------------------------
      CALL cspline_alloc(smats,mpsi,mpert**2)
      CALL cspline_alloc(tmats,mpsi,mpert**2)
      CALL cspline_alloc(xmats,mpsi,mpert**2)
      CALL cspline_alloc(ymats,mpsi,mpert**2)
      CALL cspline_alloc(zmats,mpsi,mpert**2)
      smats%xs=sq%xs
      tmats%xs=sq%xs
      xmats%xs=sq%xs
      ymats%xs=sq%xs
      zmats%xs=sq%xs
      DO ipsi=0,mpsi
         q=sq%fs(ipsi,4)
         sband(0:-mband:-1)=fmodb%cs%fs(ipsi,1:mband+1)
         tband(0:-mband:-1)=fmodb%cs%fs(ipsi,mband+2:2*mband+2)
         xband(0:-mband:-1)=fmodb%cs%fs(ipsi,2*mband+3:3*mband+3)
         yband1(0:-mband:-1)=fmodb%cs%fs(ipsi,3*mband+4:4*mband+4)
         yband2(0:-mband:-1)=fmodb%cs%fs(ipsi,4*mband+5:5*mband+5)
         zband1(0:-mband:-1)=fmodb%cs%fs(ipsi,5*mband+6:6*mband+6)
         zband2(0:-mband:-1)=fmodb%cs%fs(ipsi,6*mband+7:7*mband+7)
         zband3(0:-mband:-1)=fmodb%cs%fs(ipsi,7*mband+8:8*mband+8)

         sband(1:mband)=CONJG(sband(-1:-mband:-1))
         tband(1:mband)=CONJG(tband(-1:-mband:-1))
         xband(1:mband)=CONJG(xband(-1:-mband:-1))
         yband1(1:mband)=CONJG(yband1(-1:-mband:-1))
         yband2(1:mband)=CONJG(yband2(-1:-mband:-1))
         zband1(1:mband)=CONJG(zband1(-1:-mband:-1))
         zband2(1:mband)=CONJG(zband2(-1:-mband:-1))
         zband3(1:mband)=CONJG(zband3(-1:-mband:-1))

         ipert=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               m2=m1+dm
               singfac2=m2-nn*q
               jpert=ipert+dm
               smat(ipert,jpert)=sband(dm)
               tmat(ipert,jpert)=tband(dm)
               xmat(ipert,jpert)=xband(dm)
               ymat(ipert,jpert)=yband1(dm)+ifac*singfac2*yband2(dm)
               zmat(ipert,jpert)=zband1(dm)+
     $              ifac*(m2*zband2(dm)+nn*zband3(dm))
            ENDDO
         ENDDO
         smats%fs(ipsi,:)=RESHAPE(smat,(/mpert**2/))
         tmats%fs(ipsi,:)=RESHAPE(tmat,(/mpert**2/))
         xmats%fs(ipsi,:)=RESHAPE(xmat,(/mpert**2/))
         ymats%fs(ipsi,:)=RESHAPE(ymat,(/mpert**2/))
         zmats%fs(ipsi,:)=RESHAPE(zmat,(/mpert**2/))
      ENDDO

      CALL cspline_fit(smats,"extrap")
      CALL cspline_fit(tmats,"extrap")
      CALL cspline_fit(xmats,"extrap")
      CALL cspline_fit(ymats,"extrap")
      CALL cspline_fit(zmats,"extrap")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_action_matrix
c-----------------------------------------------------------------------
c     subprogram 3. fourfit_make_matrix.
c     constructs the coefficient matrices and fits them to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_make_matrix

      CHARACTER(128) :: message
      INTEGER :: ipsi,ipert,jpert,m1,m2,m,dm,info,iqty,l
      REAL(r8) :: chi1,jtheta,nq,p1,psifac,q,q1,singfac1,singfac2
      COMPLEX(r8) :: tphi

      INTEGER, DIMENSION(mpert) :: mfac

      COMPLEX(r8), DIMENSION(mpert*mpert) :: work
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
      COMPLEX(r8), DIMENSION(-mband:mband) ::
     $     g11,g22,g33,g23,g31,g12,jmat1,imat
      COMPLEX(r8), DIMENSION(mpert,mpert) :: amat,bmat,cmat,dmat,emat,
     $     fmat,gmat,hmat,kmat,temp0,temp1,temp2,baat,caat,eaat,
     $     fbat,kbat,gaat

      COMPLEX(r8), DIMENSION(3*mband+1,mpert) :: amatlu,fmatlu
      COMPLEX(r8), DIMENSION(mpert,mpert,6) :: kinmats,kinmats_l,
     $     kinmata,kinmata_l

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER, PARAMETER :: unit=99

      COMPLEX(r8), DIMENSION(0:mpsi,mpert**2) :: ksa,ksb,ksc,ksd,
     $        kse,ksh,kaa,kab,kac,kad,kae,kah

      mfac =(/(m,m=mlow,mhigh)/)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",5x,"j",4x,"re fb",6x,"im fb"/)
 20   FORMAT(2i6,1p,2e11.3)
c-----------------------------------------------------------------------
c     set up complex cubic splines for matrices.
c-----------------------------------------------------------------------
      IF(diagnose)CALL bin_open(unit,"coefs.bin","UNKNOWN","REWIND",
     $     "none")
      ALLOCATE(asmat(mpert,mpert),bsmat(mpert,mpert),csmat(mpert,mpert),
     $     fsmat(3*mband+1,mpert),ksmat(2*mband+1,mpert),
     $     ipiva(mpert),jmat(-mband:mband))
      CALL cspline_alloc(fmats,mpsi,(mband+1)*(2*mpert-mband)/2)
      CALL cspline_alloc(gmats,mpsi,(mband+1)*(2*mpert-mband)/2)
      CALL cspline_alloc(kmats,mpsi,(2*mband+1)*mpert)
      CALL cspline_alloc(kaats,mpsi,(2*mband+1)*mpert)
      CALL cspline_alloc(gaats,mpsi,(2*mband+1)*mpert)
      fmats%xs=rzphi%xs
      gmats%xs=rzphi%xs
      kmats%xs=rzphi%xs
      kaats%xs=rzphi%xs
      gaats%xs=rzphi%xs
      fmats%fs=0
      gmats%fs=0
      kmats%fs=0
      kaats%fs=0
      gaats%fs=0
      imat=0
      imat(0)=1
c-----------------------------------------------------------------------
c     set up cubic splines for interpolation to psilim.
c-----------------------------------------------------------------------
      CALL cspline_alloc(amats,mpsi,mpert**2)
      CALL cspline_alloc(bmats,mpsi,mpert**2)
      CALL cspline_alloc(cmats,mpsi,mpert**2)
      CALL cspline_alloc(dmats,mpsi,mpert**2)
      CALL cspline_alloc(emats,mpsi,mpert**2)
      CALL cspline_alloc(hmats,mpsi,mpert**2) 
      CALL cspline_alloc(baats,mpsi,mpert**2)
      CALL cspline_alloc(caats,mpsi,mpert**2)
      CALL cspline_alloc(eaats,mpsi,mpert**2)
      CALL cspline_alloc(fbats,mpsi,mpert**2)
      CALL cspline_alloc(kbats,mpsi,mpert**2)
  
      amats%xs=rzphi%xs
      bmats%xs=rzphi%xs
      cmats%xs=rzphi%xs
      dmats%xs=rzphi%xs
      emats%xs=rzphi%xs
      hmats%xs=rzphi%xs
      baats%xs=rzphi%xs
      caats%xs=rzphi%xs
      eaats%xs=rzphi%xs
      fbats%xs=rzphi%xs
      kbats%xs=rzphi%xs
      amats%fs=0
      bmats%fs=0
      cmats%fs=0
      dmats%fs=0
      emats%fs=0
      hmats%fs=0
      baats%fs=0
      caats%fs=0
      eaats%fs=0
      fbats%fs=0
      kbats%fs=0
c-----------------------------------------------------------------------
c     define flux surface quantities.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         p1=sq%fs1(ipsi,2)
         q=sq%fs(ipsi,4)
         q1=sq%fs1(ipsi,4)
         chi1=twopi*psio
         nq=nn*q
         jtheta=-sq%fs1(ipsi,1)
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
         g11(0:-mband:-1)=metric%cs%fs(ipsi,1:mband+1)
         g22(0:-mband:-1)=metric%cs%fs(ipsi,mband+2:2*mband+2)
         g33(0:-mband:-1)=metric%cs%fs(ipsi,2*mband+3:3*mband+3)
         g23(0:-mband:-1)=metric%cs%fs(ipsi,3*mband+4:4*mband+4)
         g31(0:-mband:-1)=metric%cs%fs(ipsi,4*mband+5:5*mband+5)
         g12(0:-mband:-1)=metric%cs%fs(ipsi,5*mband+6:6*mband+6)
         jmat(0:-mband:-1)=metric%cs%fs(ipsi,6*mband+7:7*mband+7)
         jmat1(0:-mband:-1)=metric%cs%fs(ipsi,7*mband+8:8*mband+8)
c-----------------------------------------------------------------------
c     compute upper half of matrices.
c-----------------------------------------------------------------------
         g11(1:mband)=CONJG(g11(-1:-mband:-1))
         g22(1:mband)=CONJG(g22(-1:-mband:-1))
         g33(1:mband)=CONJG(g33(-1:-mband:-1))
         g23(1:mband)=CONJG(g23(-1:-mband:-1))
         g31(1:mband)=CONJG(g31(-1:-mband:-1))
         g12(1:mband)=CONJG(g12(-1:-mband:-1))
         jmat(1:mband)=CONJG(jmat(-1:-mband:-1))
         jmat1(1:mband)=CONJG(jmat1(-1:-mband:-1))
c-----------------------------------------------------------------------
c     begin loops over perturbed fourier components.
c-----------------------------------------------------------------------
         ipert=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            singfac1=m1-nq
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               m2=m1+dm
               singfac2=m2-nq
               jpert=ipert+dm
c-----------------------------------------------------------------------
c     construct primitive matrices.
c-----------------------------------------------------------------------
               amat(ipert,jpert)=twopi**2*(nn*nn*g22(dm)
     $              +nn*(m1+m2)*g23(dm)+m1*m2*g33(dm))
               bmat(ipert,jpert)=-twopi*ifac*chi1
     $              *(nn*g22(dm)+(m1+nq)*g23(dm)+m1*q*g33(dm))
               cmat(ipert,jpert)=twopi*ifac*(
     $              twopi*ifac*chi1*singfac2*(nn*g12(dm)+m1*g31(dm))
     $              -q1*chi1*(nn*g23(dm)+m1*g33(dm)))
     $              -twopi*ifac*(jtheta*singfac1*imat(dm)
     $              +nn*p1/chi1*jmat(dm))
               dmat(ipert,jpert)=twopi*chi1*(g23(dm)+g33(dm)*m1/nn)
               emat(ipert,jpert)=-chi1/nn*(q1*chi1*g33(dm)
     $              -twopi*ifac*chi1*g31(dm)*singfac2
     $              +jtheta*imat(dm))
               hmat(ipert,jpert)=(q1*chi1)**2*g33(dm)
     $              +(twopi*chi1)**2*singfac1*singfac2*g11(dm)
     $              -twopi*ifac*chi1*dm*q1*chi1*g31(dm)
     $              +jtheta*q1*chi1*imat(dm)+p1*jmat1(dm)
               fmat(ipert,jpert)=(chi1/nn)**2*g33(dm)
               kmat(ipert,jpert)=twopi*ifac*chi1*(g23(dm)+g33(dm)*m1/nn)
            ENDDO
         ENDDO
         temp0=amat
         fbat=fmat
         kbat=kmat  
c-----------------------------------------------------------------------
c     factor A.
c-----------------------------------------------------------------------
         CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
         IF(info /= 0)THEN
            WRITE(message,'(2(a,i2))')
     $           "zhetrf: amat singular at ipsi = ",ipsi,
     $           ", ipert = ",info,", increase delta_mband"
            CALL program_stop(message)
         ENDIF
c-----------------------------------------------------------------------
c     compute composite matrices F, G, and K.
c-----------------------------------------------------------------------
         temp1=dmat
         temp2=cmat
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,temp1,mpert,info)
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,temp2,mpert,info)
         fmat=fmat-MATMUL(CONJG(TRANSPOSE(dmat)),temp1)
         kmat=emat-MATMUL(CONJG(TRANSPOSE(kmat)),temp2)
         gmat=hmat-MATMUL(CONJG(TRANSPOSE(cmat)),temp2)
c-----------------------------------------------------------------------
c     kinetic matrices
c-----------------------------------------------------------------------
         amat=temp0
         IF (kin_flag) THEN
            ipert=0
            DO m1=mlow,mhigh
               ipert=ipert+1
               singfac1=m1-nq
               DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
                  m2=m1+dm
                  singfac2=m2-nq
                  jpert=ipert+dm
                  dmat(ipert,jpert)=chi1**2*(g22(dm)+q*g23(dm)+
     $                 q*(g23(dm)+q*g33(dm)))
                  emat(ipert,jpert)=chi1**2*(q1*(g23(dm)+q*g33(dm))-
     $                 twopi*ifac*(g12(dm)+q*g31(dm))*singfac2)+
     $                 p1*jmat(dm)
               ENDDO
            ENDDO
            kinmats = 0
            kinmata = 0
            print *,psifac," - kinetic mats"
            DO l=-nl,nl
               print *,"l=",l
               kinmats_l = 0
               kinmata_l = 0
               tphi = tpsi(psifac,nn,l,zi,mi,wdfac,divxfac,electron,
     $              "twmm",keq_out,theta_out,xlmda_out,kinmats_l)
               kinmats = kinmats+kinmats_l
               tphi = tpsi(psifac,nn,l,zi,mi,wdfac,divxfac,electron,
     $              "ttmm",keq_out,theta_out,xlmda_out,kinmata_l)
               kinmata = kinmata+kinmata_l
            ENDDO

            kinmats=kinfac1*kinmats
            kinmata=kinfac2*kinmata

            kinmats(:,:,1)=kinmats(:,:,1)*2*mu0/chi1**2
            kinmata(:,:,1)=kinmata(:,:,1)*2*mu0/chi1**2
            kinmats(:,:,2)=kinmats(:,:,2)*2*mu0/chi1
            kinmata(:,:,2)=kinmata(:,:,2)*2*mu0/chi1
            kinmats(:,:,3)=kinmats(:,:,3)*2*mu0/chi1
            kinmata(:,:,3)=kinmata(:,:,3)*2*mu0/chi1
            kinmats(:,:,4)=kinmats(:,:,4)*2*mu0
            kinmata(:,:,4)=kinmata(:,:,4)*2*mu0
            kinmats(:,:,5)=kinmats(:,:,5)*2*mu0
            kinmata(:,:,5)=kinmata(:,:,5)*2*mu0
            kinmats(:,:,6)=kinmats(:,:,6)*2*mu0
            kinmata(:,:,6)=kinmata(:,:,6)*2*mu0

            amat=amat+kinmats(:,:,1)+kinmata(:,:,1)
            bmat=bmat+kinmats(:,:,2)+kinmata(:,:,2)
            cmat=cmat+kinmats(:,:,3)+kinmata(:,:,3)
            dmat=dmat+kinmats(:,:,4)+kinmata(:,:,4)
            emat=emat+kinmats(:,:,5)+kinmata(:,:,5)
            hmat=hmat+kinmats(:,:,6)+kinmata(:,:,6)
            baat=bmat-2*kinmata(:,:,2)
            caat=cmat-2*kinmata(:,:,3)
            eaat=emat-2*kinmata(:,:,5)
c-----------------------------------------------------------------------
c     only for diagnostic purposes.
c-----------------------------------------------------------------------
            IF(.TRUE.)THEN
               ksa(ipsi,:)=RESHAPE(kinmats(:,:,1),(/mpert**2/))
               ksb(ipsi,:)=RESHAPE(kinmats(:,:,2),(/mpert**2/))
               ksc(ipsi,:)=RESHAPE(kinmats(:,:,3),(/mpert**2/))
               ksd(ipsi,:)=RESHAPE(kinmats(:,:,4),(/mpert**2/))
               kse(ipsi,:)=RESHAPE(kinmats(:,:,5),(/mpert**2/))
               ksh(ipsi,:)=RESHAPE(kinmats(:,:,6),(/mpert**2/))
               kaa(ipsi,:)=RESHAPE(kinmata(:,:,1),(/mpert**2/))
               kab(ipsi,:)=RESHAPE(kinmata(:,:,2),(/mpert**2/))
               kac(ipsi,:)=RESHAPE(kinmata(:,:,3),(/mpert**2/))
               kad(ipsi,:)=RESHAPE(kinmata(:,:,4),(/mpert**2/))
               kae(ipsi,:)=RESHAPE(kinmata(:,:,5),(/mpert**2/))
               kah(ipsi,:)=RESHAPE(kinmata(:,:,6),(/mpert**2/))
            ENDIF
            amatlu=0
            DO jpert=1,mpert
               DO ipert=1,mpert
                  amatlu(2*mband+1+ipert-jpert,jpert)=amat(ipert,jpert)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     factor non-hermitian matrix A.
c-----------------------------------------------------------------------
            CALL zgbtrf(mpert,mpert,mband,mband,amatlu,2*mband+mband+1,
     $           ipiva,info)
            IF(info /= 0)THEN
               WRITE(message,'(2(a,i2))')
     $              "zgbtrf: amat singular at ipsi = ",ipsi,
     $              ", ipert = ",info,", increase delta_mband"
               CALL program_stop(message)
            ENDIF
c-----------------------------------------------------------------------
c     modify matrix G.
c-----------------------------------------------------------------------
            temp2=cmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           2*mband+mband+1,ipiva,temp2,mpert,info)
            gaat=hmat-MATMUL(CONJG(TRANSPOSE(caat)),temp2) 

            iqty=1
            DO jpert=1,mpert
               DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
                  gaats%fs(ipsi,iqty)=gaat(ipert,jpert)
                  iqty=iqty+1
               ENDDO
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c     store matrices for interpolation.
c-----------------------------------------------------------------------
         amats%fs(ipsi,:)=RESHAPE(amat,(/mpert**2/))
         bmats%fs(ipsi,:)=RESHAPE(bmat,(/mpert**2/))
         cmats%fs(ipsi,:)=RESHAPE(cmat,(/mpert**2/))
         dmats%fs(ipsi,:)=RESHAPE(dmat,(/mpert**2/))
         emats%fs(ipsi,:)=RESHAPE(emat,(/mpert**2/))
         hmats%fs(ipsi,:)=RESHAPE(hmat,(/mpert**2/))
         baats%fs(ipsi,:)=RESHAPE(baat,(/mpert**2/))
         caats%fs(ipsi,:)=RESHAPE(caat,(/mpert**2/))
         eaats%fs(ipsi,:)=RESHAPE(eaat,(/mpert**2/))
         fbats%fs(ipsi,:)=RESHAPE(fbat,(/mpert**2/))
         kbats%fs(ipsi,:)=RESHAPE(kbat,(/mpert**2/))
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
         IF(feval_flag)CALL fourfit_evals(ipsi,psifac,fmat)
         IF(diagnose)WRITE(unit)REAL(psifac,4),
     $        REAL(fmat(1,1),4),REAL(kmat(1,1),4),
     $        REAL(gmat(1,1)*psifac,4),REAL(g11(0)*psifac,4),
     $        REAL(g22(0),4),REAL(g33(0),4)
c-----------------------------------------------------------------------
c     transfer F to banded matrix.
c-----------------------------------------------------------------------
         DO jpert=1,mpert
            DO ipert=jpert,MIN(mpert,jpert+mband)
               fmatb(1+ipert-jpert,jpert)=fmat(ipert,jpert)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     factor F.
c-----------------------------------------------------------------------
         CALL zpbtrf('L',mpert,mband,fmatb,mband+1,info)
         IF(info /= 0)THEN
            WRITE(message,'(2(a,i3),a)')
     $           "zpbtrf: fmat singular at ipsi = ",ipsi,
     $           ", ipert = ",info,", reduce delta_mband"
            CALL program_stop(message)
         ENDIF
c-----------------------------------------------------------------------
c     store Hermitian matrices F and G.
c-----------------------------------------------------------------------
         iqty=1
         DO jpert=1,mpert
            DO ipert=jpert,MIN(mpert,jpert+mband)
               fmats%fs(ipsi,iqty)=fmatb(1+ipert-jpert,jpert)
               gmats%fs(ipsi,iqty)=gmat(ipert,jpert)
               iqty=iqty+1
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     store non-Hermitian matrix K.
c-----------------------------------------------------------------------
         iqty=1
         DO jpert=1,mpert
            DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
               kmats%fs(ipsi,iqty)=kmat(ipert,jpert)
               iqty=iqty+1
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     close diagnostic file.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(unit)
         CALL bin_close(unit)
      ENDIF
c-----------------------------------------------------------------------
c     set powers.
c-----------------------------------------------------------------------
      gmats%xpower(1,:)=-1
      IF(power_flag)THEN
         m=mlow
         iqty=1
         DO jpert=1,mpert
            DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
               dm=ipert-jpert
               IF(m == 1 .AND. dm == -1 .OR. m == -1 .AND. dm == 1)
     $              kmats%xpower(1,iqty)=-1
               iqty=iqty+1
            ENDDO
            m=m+1
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     cubic spline fit banded matrices.
c-----------------------------------------------------------------------
      CALL cspline_fit(amats,"extrap")
      CALL cspline_fit(bmats,"extrap")
      CALL cspline_fit(cmats,"extrap")
      CALL cspline_fit(dmats,"extrap")
      CALL cspline_fit(emats,"extrap")
      CALL cspline_fit(hmats,"extrap")
      CALL cspline_fit(baats,"extrap")
      CALL cspline_fit(caats,"extrap")
      CALL cspline_fit(eaats,"extrap")
      CALL cspline_fit(fmats,"extrap")
      CALL cspline_fit(kmats,"extrap")
      CALL cspline_fit(gmats,"extrap")
      CALL cspline_fit(gaats,"extrap")
      CALL cspline_fit(fbats,"extrap")
      CALL cspline_fit(kbats,"extrap")
c-----------------------------------------------------------------------
c     write binary output for diagnosis.
c-----------------------------------------------------------------------  
      IF(.TRUE.)THEN
         WRITE(*,*)"Write binary output for graphs."
         CALL bin_open(bin_unit,"fs.bin","UNKNOWN","REWIND","none")
         DO ipert=1,(mband+1)*(2*mpert-mband)/2
            DO ipsi=0,mpsi
               WRITE(bin_unit)REAL(sq%xs(ipsi),4),
     $              REAL(sq%fs(ipsi,4),4),
     $              REAL(REAL(fmats%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(fmats%fs(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"ks.bin","UNKNOWN","REWIND","none")
         DO ipert=1,(2*mband+1)*mpert
            DO ipsi=0,mpsi
               WRITE(bin_unit)REAL(sq%xs(ipsi),4),
     $              REAL(sq%fs(ipsi,4),4),
     $              REAL(REAL(kmats%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kmats%fs(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"gs.bin","UNKNOWN","REWIND","none")
         DO ipert=1,(mband+1)*(2*mpert-mband)/2
            DO ipsi=0,mpsi
               WRITE(bin_unit)REAL(sq%xs(ipsi),4),
     $              REAL(sq%fs(ipsi,4),4),
     $              REAL(REAL(gmats%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(gmats%fs(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"kss.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert**2
            DO ipsi=0,mpsi
               WRITE(bin_unit)REAL(sq%xs(ipsi),4),
     $              REAL(REAL(ksa(ipsi,ipert)),4),
     $              REAL(AIMAG(ksa(ipsi,ipert)),4),
     $              REAL(REAL(ksb(ipsi,ipert)),4),
     $              REAL(AIMAG(ksb(ipsi,ipert)),4),
     $              REAL(REAL(ksc(ipsi,ipert)),4),
     $              REAL(AIMAG(ksc(ipsi,ipert)),4),
     $              REAL(REAL(ksd(ipsi,ipert)),4),
     $              REAL(AIMAG(ksd(ipsi,ipert)),4),
     $              REAL(REAL(kse(ipsi,ipert)),4),
     $              REAL(AIMAG(kse(ipsi,ipert)),4),
     $              REAL(REAL(ksh(ipsi,ipert)),4),
     $              REAL(AIMAG(ksh(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"kas.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert**2
            DO ipsi=0,mpsi
               WRITE(bin_unit)REAL(sq%xs(ipsi),4),
     $              REAL(REAL(kaa(ipsi,ipert)),4),
     $              REAL(AIMAG(kaa(ipsi,ipert)),4),
     $              REAL(REAL(kab(ipsi,ipert)),4),
     $              REAL(AIMAG(kab(ipsi,ipert)),4),
     $              REAL(REAL(kac(ipsi,ipert)),4),
     $              REAL(AIMAG(kac(ipsi,ipert)),4),
     $              REAL(REAL(kad(ipsi,ipert)),4),
     $              REAL(AIMAG(kad(ipsi,ipert)),4),
     $              REAL(REAL(kae(ipsi,ipert)),4),
     $              REAL(AIMAG(kae(ipsi,ipert)),4),
     $              REAL(REAL(kah(ipsi,ipert)),4),
     $              REAL(AIMAG(kah(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     store matrices for plasma inductance and permeability.
c-----------------------------------------------------------------------
      CALL cspline_eval(amats,psilim,0)
      CALL cspline_eval(bmats,psilim,0)
      CALL cspline_eval(cmats,psilim,0)
      CALL cspline_eval(dmats,psilim,0)
      CALL cspline_eval(emats,psilim,0)
      CALL cspline_eval(baats,psilim,0)
      amat=RESHAPE(amats%f,(/mpert,mpert/))
      bmat=RESHAPE(bmats%f,(/mpert,mpert/))
      cmat=RESHAPE(cmats%f,(/mpert,mpert/))
      dmat=RESHAPE(dmats%f,(/mpert,mpert/))
      emat=RESHAPE(emats%f,(/mpert,mpert/))
      baat=RESHAPE(baats%f,(/mpert,mpert/))

      amatlu=0
      DO jpert=1,mpert
         DO ipert=1,mpert
            amatlu(2*mband+1+ipert-jpert,jpert)=amat(ipert,jpert)
         ENDDO
      ENDDO
         
      CALL zgbtrf(mpert,mpert,mband,mband,amatlu,3*mband+1,
     $     ipiva,info)
      IF(info /= 0)THEN
         WRITE(message,'(a,e16.9,a,i2)')
     $        "zgbtrf: amat singular at psifac = ",psilim,
     $        ", ipert = ",info,", reduce delta_mband"
         CALL program_stop(message)
      ENDIF
      temp1=bmat
      temp2=cmat
      CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $     3*mband+1,ipiva,temp1,mpert,info)
      CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $     3*mband+1,ipiva,temp2,mpert,info)
      fmat=dmat-MATMUL(CONJG(TRANSPOSE(baat)),temp1)
      kmat=emat-MATMUL(CONJG(TRANSPOSE(baat)),temp2)

      fsmat=0
      DO jpert=1,mpert
         DO ipert=1,mpert
            fsmat(2*mband+1+ipert-jpert,jpert)=fmat(ipert,jpert)
         ENDDO
      ENDDO

      ksmat=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            ksmat(1+mband+ipert-jpert,jpert)=kmat(ipert,jpert)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     interpolate matrices to psilim (need modification for valen3d).
c-----------------------------------------------------------------------
      IF(sas_flag)THEN
         asmat=RESHAPE(amats%f,(/mpert,mpert/))
         bsmat=RESHAPE(bmats%f,(/mpert,mpert/))
         csmat=RESHAPE(cmats%f,(/mpert,mpert/))
         CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose and deallocate.
c-----------------------------------------------------------------------
      IF(bin_metric)CALL fourfit_write_metric
      IF(bin_fmat)CALL fourfit_write_matrix(fmats,"fmat",.TRUE.)
      IF(bin_gmat)CALL fourfit_write_matrix(gmats,"gmat",.TRUE.)
      IF(bin_kmat)CALL fourfit_write_matrix(kmats,"kmat",.FALSE.)
      CALL fspline_dealloc(metric)
      CALL fspline_dealloc(fmodb)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_make_matrix
c-----------------------------------------------------------------------
c     subprogram 4. fourfit_write_metric.
c     uses cspline_write_log to diagnose cspline_types.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_write_metric

      REAL(r8), DIMENSION(2) :: xend=(/zero,one/)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"metric.bin","UNKNOWN","REWIND","none")
      CALL cspline_write_log(metric%cs,.FALSE.,.TRUE.,out_unit,
     $     bin_unit,.FALSE.,mband+1,xend)
      CALL bin_close(fourfit_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_write_metric
c-----------------------------------------------------------------------
c     subprogram 5. fourfit_write_matrix
c     produces ascii and binary output of logs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_write_matrix(matrix,name,sym_flag)

      TYPE(cspline_type) :: matrix
      CHARACTER(*), INTENT(IN) :: name
      LOGICAL, INTENT(IN) :: sym_flag

      INTEGER :: iqty,ix,jx,j,ipert,jpert,iband,iband0,mband0=8,imin
      REAL(r8) :: dx
      REAL(r8), DIMENSION(2) :: xend=(/zero,one/)
      REAL(r8), DIMENSION(0:4*mpsi) :: x
      REAL(r8), DIMENSION(0:4*mpsi,2) :: xlog

      REAL(r8), DIMENSION(:), POINTER :: epsilon
      REAL(r8), DIMENSION(:,:,:), POINTER :: flog
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: f
c-----------------------------------------------------------------------
c     set limits and allocate arrays.
c-----------------------------------------------------------------------
      mband0=MIN(mband0,mband)
      IF(sym_flag)THEN
         iband0=0
      ELSE
         iband0=-mband0
      ENDIF
      ALLOCATE(f(0:4*mpsi,mpert,iband0:mband0),
     $     flog(0:4*mpsi,mpert,iband0:mband0),
     $     epsilon(iband0:mband0))
c-----------------------------------------------------------------------
c     compute values.
c-----------------------------------------------------------------------
      f=0
      jx=0
      DO ix=1,mpsi
         dx=(matrix%xs(ix)-matrix%xs(ix-1))/4
         DO j=0,4
            IF(j == 4 .AND. ix < mpsi)CYCLE
            x(jx)=matrix%xs(ix-1)+j*dx
            xlog(jx,:)=LOG10(ABS(x(jx)-xend))
            CALL cspline_eval(matrix,x(jx),0)
            iqty=1
            DO jpert=1,mpert
               IF(sym_flag)THEN
                  imin=jpert
               ELSE
                  imin=MAX(1,jpert-mband)
               ENDIF
               DO ipert=imin,MIN(mpert,jpert+mband)
                  IF(ipert-jpert <= mband0
     $                 .AND. ipert-jpert >= -mband0)
     $                 f(jx,jpert,ipert-jpert)=matrix%f(iqty)
                  iqty=iqty+1
               ENDDO
            ENDDO
            jx=jx+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute logs.
c-----------------------------------------------------------------------
      flog=HUGE(flog)
      WHERE(f /= 0)
         flog=LOG(f)
      ELSEWHERE
         flog=HUGE(flog)
      ENDWHERE
      DO iband=iband0,mband0
         epsilon(iband)=MINVAL(flog(:,:,iband))
         WHERE(f(:,:,iband) == 0)
            flog(:,:,iband)=epsilon(iband)
         ENDWHERE
      ENDDO
      flog=flog/alog10
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,TRIM(name)//".bin","UNKNOWN","REWIND",
     $     "none")
c-----------------------------------------------------------------------
c     print node values.
c-----------------------------------------------------------------------
      DO ipert=1,mpert
         DO ix=0,4*mpsi,4
            WRITE(bin_unit)REAL(x(ix),4),REAL(xlog(ix,:),4),
     $           REAL(flog(ix,ipert,:),4)
         ENDDO
         WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO ipert=1,mpert
         DO ix=0,4*mpsi
            WRITE(bin_unit)REAL(x(ix),4),REAL(xlog(ix,:),4),
     $           REAL(flog(ix,ipert,:),4)
         ENDDO
         WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c      close output file and deallocate arrays.
c-----------------------------------------------------------------------
      CALL bin_close(bin_unit)
      DEALLOCATE(f,flog,epsilon)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_write_matrix
c-----------------------------------------------------------------------
c     subprogram 6. fourfit_evals
c     computes and prints eigenvalues Hermitian matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_evals(ipsi,psifac,matrix)
      USE global_mod
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ipsi
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: matrix

      INTEGER :: info,lwork
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp
      COMPLEX(r8), DIMENSION(2*(mpert+1)*mpert) :: work
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      REAL(r8), DIMENSION(mpert) :: evals
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"ipsi",3x,"psifac",5x,"eval1",6x,"eval2"/)
 20   FORMAT(i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(ipsi == 0)THEN
         CALL ascii_open(evals_out_unit,"feval.out","UNKNOWN")
         WRITE(evals_out_unit,10)
         CALL bin_open(evals_bin_unit,"feval.bin","UNKNOWN","REWIND",
     $        "none")
      ENDIF
c-----------------------------------------------------------------------
c     compute eigenvalues.
c-----------------------------------------------------------------------
      lwork=SIZE(work)
      temp=matrix
      CALL zheev('N','U',mpert,temp,mpert,evals,work,lwork,rwork,info)
c-----------------------------------------------------------------------
c     print results.
c-----------------------------------------------------------------------
      WRITE(evals_out_unit,20)ipsi,psifac,evals(1:2)
      WRITE(evals_bin_unit)REAL(psifac,4),REAL(evals(1:2),4)
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(ipsi == mpsi)THEN
         WRITE(evals_out_unit,10)
         CALL ascii_close(evals_out_unit)
         CALL bin_close(evals_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_evals
c-----------------------------------------------------------------------
c     subprogram 7. fourfit_diagnose_1.
c     diagnoses coefficient matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_diagnose_1(g11,g22,g33,g23,g31,g12)

      COMPLEX(r8), DIMENSION(-mband:mband), INTENT(IN) ::
     $     g11,g22,g33,g23,g31,g12

      INTEGER :: dm,unit=98
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"dm",5x,"g11",8x,"g22",8x,"g33",8x,"g23",8x,"g31",8x,
     $     "g12"/)
 20   FORMAT(i6,1p,8e11.3)
c-----------------------------------------------------------------------
c     write binary output.
c-----------------------------------------------------------------------
      CALL ascii_open(unit,"metric.out","UNKNOWN")
      WRITE(unit,10)
      DO dm=0,mband
         WRITE(unit,20)dm,
     $        REAL(g11(dm),4),REAL(g22(dm),4),REAL(g33(dm),4),
     $        REAL(g23(dm),4),REAL(g31(dm),4),REAL(g12(dm),4)
      ENDDO
      WRITE(unit,10)
      CALL ascii_close(unit)
      CALL program_stop("Termination by fourfit_diagnose_1") 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_diagnose_1
      END MODULE fourfit_mod
