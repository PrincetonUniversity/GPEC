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
c     8. fourfit_kinetic_matrix.
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
      USE pentrc_interface,         ! rename overlapping names
     $    pentrc_verbose=>verbose,  ! should get a more fundamental fix
     $    pentrc_mpert=>mpert,
     $    pentrc_nn=>nn,
     $    pentrc_r8=>r8,
     $    pentrc_timer=>timer
      USE utilities, only : progressbar
      USE torque, only : kelmm      ! cspline Euler-Lagrange mats for local use
      USE inputs, only : dbob_m,divx_m,kin,xs_m,fnml
      USE energy_integration
      USE pitch_integration
      USE dcon_interface, only: geom,
     $     dcon_int_rzphi=>rzphi,
     $     dcon_int_eqfun=>eqfun,
     $     dcon_int_sq=>sq,
     $     dcon_int_smats=>smats,
     $     dcon_int_tmats=>tmats,
     $     dcon_int_xmats=>xmats,
     $     dcon_int_ymats=>ymats,
     $     dcon_int_zmats=>zmats
      IMPLICIT NONE

      TYPE(fspline_type), PRIVATE :: metric,fmodb
      TYPE(cspline_type) :: dmats,emats,hmats,dbats,ebats,fbats,
     $     fmats,kmats,gmats,kaats,gaats,
     $     f0mats,pmats,paats,kkmats,kkaats,r1mats,r2mats,r3mats,
     $     akmats,bkmats,ckmats
      TYPE(spline_type) :: k0s
      INTEGER, DIMENSION(:), POINTER :: ipiva
      COMPLEX(r8), DIMENSION(:,:), POINTER :: asmat,bsmat,csmat
      COMPLEX(r8), DIMENSION(:), POINTER :: jmat

      INTEGER :: parallel_threads = 0, dcon_kin_threads = 0

      ! kientic ABCDEH mats for sing_mod
      TYPE(cspline_type) :: kwmats(6),ktmats(6)

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
c     subprogram 2. fourfit_make_matrix.
c     constructs the coefficient matrices and fits them to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_make_matrix(op_diagnose)

      LOGICAL, OPTIONAL, INTENT(IN) :: op_diagnose

      CHARACTER(128) :: message
      INTEGER :: ipsi,ipert,jpert,m1,m2,m,dm,info,iqty,l,i,j,iindex
      REAL(r8) :: chi1,jtheta,nq,p1,psifac,q,q1,singfac1,singfac2,ileft
      COMPLEX(r8) :: tphi

      INTEGER, DIMENSION(mpert) :: mfac,ipiv

      COMPLEX(r8), DIMENSION(mpert*mpert) :: work
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
      COMPLEX(r8), DIMENSION(-mband:mband) ::
     $     g11,g22,g33,g23,g31,g12,jmat1,imat
      COMPLEX(r8), DIMENSION(mpert,mpert) :: amat,bmat,cmat,dmat,emat,
     $     fmat,gmat,hmat,kmat,temp0,temp1,temp2,dbat,ebat,fbat
      COMPLEX(r8), DIMENSION(3*mband+1,mpert) :: amatlu,fmatlu

      LOGICAL :: diagnose=.FALSE.
      INTEGER, PARAMETER :: unit=99

      mfac =(/(m,m=mlow,mhigh)/)
      IF(PRESENT(op_diagnose)) diagnose = op_diagnose
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
     $     jmat(-mband:mband),ipiva(mpert))
      CALL cspline_alloc(fmats,mpsi,(mband+1)*(2*mpert-mband)/2)
      CALL cspline_alloc(gmats,mpsi,(mband+1)*(2*mpert-mband)/2)
      CALL cspline_alloc(kmats,mpsi,(2*mband+1)*mpert)
      fmats%xs=rzphi%xs
      gmats%xs=rzphi%xs
      kmats%xs=rzphi%xs
      fmats%fs=0
      gmats%fs=0
      kmats%fs=0
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
      CALL cspline_alloc(dbats,mpsi,mpert**2) 
      CALL cspline_alloc(ebats,mpsi,mpert**2) 
      CALL cspline_alloc(fbats,mpsi,mpert**2) 
  
      amats%xs=rzphi%xs
      bmats%xs=rzphi%xs
      cmats%xs=rzphi%xs
      dmats%xs=rzphi%xs
      emats%xs=rzphi%xs
      hmats%xs=rzphi%xs
      dbats%xs=rzphi%xs
      ebats%xs=rzphi%xs
      fbats%xs=rzphi%xs

      amats%fs=0
      bmats%fs=0
      cmats%fs=0
      dmats%fs=0
      emats%fs=0
      hmats%fs=0
      dbats%fs=0
      ebats%fs=0
      fbats%fs=0
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
         dbat=dmat
         ebat=emat
         fbat=fmat
c-----------------------------------------------------------------------
c     factor A.
c-----------------------------------------------------------------------
         temp0=amat
         CALL zhetrf('L',mpert,amat,mpert,ipiv,work,mpert*mpert,info)
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
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiv,temp1,mpert,info)
         CALL zhetrs('L',mpert,mpert,amat,mpert,ipiv,temp2,mpert,info)
         fmat=fmat-MATMUL(CONJG(TRANSPOSE(dmat)),temp1)
         kmat=emat-MATMUL(CONJG(TRANSPOSE(kmat)),temp2)
         gmat=hmat-MATMUL(CONJG(TRANSPOSE(cmat)),temp2)
         amat=temp0
c-----------------------------------------------------------------------
c     kinetic matrices.
c-----------------------------------------------------------------------
         ipert=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            singfac1=m1-nq
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               m2=m1+dm
               singfac2=m2-nq
               jpert=ipert+dm
               dmat(ipert,jpert)=chi1**2*(g22(dm)+q*g23(dm)+
     $              q*(g23(dm)+q*g33(dm)))
               emat(ipert,jpert)=chi1**2*(q1*(g23(dm)+q*g33(dm))-
     $              twopi*ifac*(g12(dm)+q*g31(dm))*singfac2)+
     $              p1*jmat(dm)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     store matrices for interpolation.
c-----------------------------------------------------------------------
         amats%fs(ipsi,:)=RESHAPE(amat,(/mpert**2/))
         bmats%fs(ipsi,:)=RESHAPE(bmat,(/mpert**2/))
         cmats%fs(ipsi,:)=RESHAPE(cmat,(/mpert**2/))
         dmats%fs(ipsi,:)=RESHAPE(dmat,(/mpert**2/))
         emats%fs(ipsi,:)=RESHAPE(emat,(/mpert**2/))
         hmats%fs(ipsi,:)=RESHAPE(hmat,(/mpert**2/))
         dbats%fs(ipsi,:)=RESHAPE(dbat,(/mpert**2/))
         ebats%fs(ipsi,:)=RESHAPE(ebat,(/mpert**2/))
         fbats%fs(ipsi,:)=RESHAPE(fbat,(/mpert**2/))
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
         IF (kin_flag)THEN
            hmats%x0(2)=1.0
            hmats%xpower(1,:)=-1
            hmats%xpower(2,:)=-1
            gmats%x0(2)=1.0
            gmats%xpower(2,:)=-1
         ENDIF
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
      CALL cspline_fit(fmats,"extrap")
      CALL cspline_fit(kmats,"extrap")
      CALL cspline_fit(gmats,"extrap")
      CALL cspline_fit(dbats,"extrap")
      CALL cspline_fit(ebats,"extrap")
      CALL cspline_fit(fbats,"extrap")

c-----------------------------------------------------------------------
c     write binary output for diagnosis.
c-----------------------------------------------------------------------  
      IF(diagnose)THEN
         WRITE(*,*)"Write binary output for graphs."
         mfac =(/(i,i=mlow,mhigh)/)
         CALL ascii_open(fourfit_out_unit,"imats.out","UNKNOWN")
         WRITE(fourfit_out_unit,*)"DCON ideal energy matrices"
         WRITE(fourfit_out_unit,'(1/,1x,a12,1x,I6,1x,1(a12,I4),1/)')
     $        "mpsi =",mpsi,"mpert =",mpert
         WRITE(fourfit_out_unit,'(1x,a16,2(1x,a4),12(1x,a16))')
     $        "psi","m1","m2",
     $        "real(Ai)","imag(Ai)","real(Bi)","imag(Bi)",
     $        "real(Ci)","imag(Ci)","real(Di)","imag(Di)",
     $        "real(Ei)","imag(Ei)","real(Hi)","imag(Hi)"
         DO ipsi=0,mpsi
            DO i=1,mpert
               DO j=1,mpert
                  ipert = (i-1)*mpert + j
                  WRITE(fourfit_out_unit,'(1x,es16.8,2(1x,I4),'//
     $                 '12(1x,es16.8))')
     $                 REAL(amats%xs(ipsi),4),mfac(i),mfac(j),
     $                 REAL(REAL(amats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(amats%fs(ipsi,ipert)),4),
     $                 REAL(REAL(bmats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(bmats%fs(ipsi,ipert)),4),
     $                 REAL(REAL(cmats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(cmats%fs(ipsi,ipert)),4),
     $                 REAL(REAL(dmats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(dmats%fs(ipsi,ipert)),4),
     $                 REAL(REAL(emats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(emats%fs(ipsi,ipert)),4),
     $                 REAL(REAL(hmats%fs(ipsi,ipert)),4),
     $                 REAL(AIMAG(hmats%fs(ipsi,ipert)),4)
               ENDDO
            ENDDO
         ENDDO
         WRITE(fourfit_out_unit,*)
         CALL ascii_close(fourfit_out_unit)

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
      ENDIF
c-----------------------------------------------------------------------
c     interpolate matrices to psilim (need modification for valen3d).
c-----------------------------------------------------------------------
      IF(sas_flag)THEN
         CALL cspline_eval(amats,psilim,0)
         CALL cspline_eval(bmats,psilim,0)
         CALL cspline_eval(cmats,psilim,0)
         asmat=RESHAPE(amats%f,(/mpert,mpert/))
         bsmat=RESHAPE(bmats%f,(/mpert,mpert/))
         csmat=RESHAPE(cmats%f,(/mpert,mpert/))
         CALL zhetrf('L',mpert,asmat,mpert,ipiva,work,mpert*mpert,info)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose and deallocate.
c-----------------------------------------------------------------------
      IF(bin_metric)CALL fourfit_write_metric
      IF(bin_fmat)CALL fourfit_write_matrix(fmats,"fmat",.TRUE.)
      IF(bin_gmat)CALL fourfit_write_matrix(gmats,"gmat",.TRUE.)
      IF(bin_kmat)CALL fourfit_write_matrix(kmats,"kmat",.FALSE.)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_make_matrix
c-----------------------------------------------------------------------
c     subprogram 3. fourfit_action_matrix.
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

      CALL fspline_dealloc(metric)
      CALL fspline_dealloc(fmodb)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_action_matrix
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
c-----------------------------------------------------------------------
c     subprogram 8. fourfit_kinetic_matrix.
c     Use PENTRC to calculated the coefficient matrices on a dynamic
c     grid and fit them to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_kinetic_matrix(methodin,writein)
      IMPLICIT NONE
      
      LOGICAL, OPTIONAL :: writein
      INTEGER, OPTIONAL :: methodin

      LOGICAL :: output
      INTEGER :: ipsi,ipert,l,i,j,iindex,method = 0
      CHARACTER(1) :: ft
      INTEGER, DIMENSION(mpert) :: mfac
      REAL(r8) :: ileft,psifac,chi1,plim(2)
      COMPLEX(r8) :: tphi
      COMPLEX(r8), DIMENSION(mpert,mpert,6) :: kwmat,kwmat_l,
     $     ktmat,ktmat_l
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: kwmatls, ktmatls
c-----------------------------------------------------------------------
c     declarations for diagnostics.
c-----------------------------------------------------------------------
      INTEGER :: info,iqty,jpert,ix,idone
      INTEGER, DIMENSION(mpert) :: ipiv
      CHARACTER(128) :: message
      COMPLEX(r8), DIMENSION(mpert*mpert) :: work,amats_f,bmats_f,
     $     cmats_f,dmats_f,emats_f,hmats_f,dbats_f,ebats_f,fbats_f
      COMPLEX(r8), DIMENSION(mpert,mpert) :: amat,bmat,cmat,dmat,emat,
     $     fmat,gmat,hmat,kmat,temp0,temp1,temp2,caat,gaat,
     $     f0mat,pmat,paat,kkmat,kkaat,r1mat,r2mat,r3mat,
     $     dbat,ebat,umat,aamat,bkmat,bkaat,b1mat
      COMPLEX(r8), DIMENSION(3*mband+1,mpert) :: amatlu
c-----------------------------------------------------------------------
c     declarations for parallelization.
c-----------------------------------------------------------------------
      LOGICAL :: debug_omp = .FALSE.
      INTEGER :: t_start, t_end, cr, num_threads, ithread
      REAL(r8) :: t_sec
      INTEGER :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
c-----------------------------------------------------------------------
c     output formats
c-----------------------------------------------------------------------
      IF(debug_omp)THEN
         PRINT *,"In serial region..."
#ifdef _OPENMP
         num_threads = OMP_GET_NUM_THREADS()
         PRINT *,"# of OMP threads = ",num_threads
#else
         num_threads = 1
         PRINT *,"# NO OPENMP MULTITHREADING"
#endif
      ENDIF
c-----------------------------------------------------------------------
c     some basic variables
c-----------------------------------------------------------------------
      IF(PRESENT(methodin)) method = methodin
      output=.FALSE.
      IF(PRESENT(writein)) output = writein
      chi1=twopi*psio
      plim = (/0.0,1.0/)
      IF (passing_flag .AND. trapped_flag) THEN
         ft="f"
      ELSE IF (trapped_flag) THEN
         ft="t"
      ELSE IF (passing_flag) THEN
         ft="p"
      ELSE
         CALL program_stop("Kinetic calculations require "//
     $                "passing_flag and/or trapped_flag")
      ENDIF
      ALLOCATE(kwmatls(mpert,mpert,6,0:mpsi,-nl:nl),
     $         ktmatls(mpert,mpert,6,0:mpsi,-nl:nl))
      kwmatls(:,:,:,:,:)=0
      ktmatls(:,:,:,:,:)=0
c-----------------------------------------------------------------------
c     Original approach using eqgrid loop to calculate kinetic matrices
c-----------------------------------------------------------------------
      CALL cspline_alloc(kaats,mpsi,(2*mband+1)*mpert)
      CALL cspline_alloc(gaats,mpsi,(2*mband+1)*mpert)
      kaats%xs=rzphi%xs
      gaats%xs=rzphi%xs
      kaats%fs=0
      gaats%fs=0

      IF(method==-1)THEN
         output = .FALSE.
      ELSEIF(method==0)THEN
         fkg_kmats_flag=.TRUE.

         CALL cspline_alloc(akmats,mpsi,mpert**2)
         CALL cspline_alloc(bkmats,mpsi,mpert**2)
         CALL cspline_alloc(ckmats,mpsi,mpert**2)
         CALL cspline_alloc(f0mats,mpsi,mpert**2)
         CALL cspline_alloc(pmats,mpsi,mpert**2)
         CALL cspline_alloc(paats,mpsi,mpert**2)
         CALL cspline_alloc(kkmats,mpsi,mpert**2)
         CALL cspline_alloc(kkaats,mpsi,mpert**2)
         CALL cspline_alloc(r1mats,mpsi,mpert**2)
         CALL cspline_alloc(r2mats,mpsi,mpert**2)
         CALL cspline_alloc(r3mats,mpsi,mpert**2)

         akmats%xs=rzphi%xs
         bkmats%xs=rzphi%xs
         ckmats%xs=rzphi%xs
         f0mats%xs=rzphi%xs
         pmats%xs=rzphi%xs
         paats%xs=rzphi%xs
         kkmats%xs=rzphi%xs
         kkaats%xs=rzphi%xs
         r1mats%xs=rzphi%xs
         r2mats%xs=rzphi%xs
         r3mats%xs=rzphi%xs

         akmats%fs=0
         bkmats%fs=0
         ckmats%fs=0
         f0mats%fs=0
         pmats%fs=0
         paats%fs=0
         kkmats%fs=0
         kkaats%fs=0
         r1mats%fs=0
         r2mats%fs=0
         r3mats%fs=0

         DO i=1,6
            CALL cspline_alloc(kwmats(i),mpsi,mpert**2)
            CALL cspline_alloc(ktmats(i),mpsi,mpert**2)
            kwmats(i)%xs=rzphi%xs
            ktmats(i)%xs=rzphi%xs
         ENDDO

c-----------------------------------------------------------------------
c     Parallel calls of expensive kinetic calculations for each psi,ell
c-----------------------------------------------------------------------
         idone = 0
         CALL SYSTEM_CLOCK(COUNT_RATE=cr)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,8) DEFAULT(NONE) COLLAPSE(2)
!$OMP& PROC_BIND(CLOSE)
!$OMP& PRIVATE(
!$OMP& ipsi,ix,psifac,l,num_threads,t_start,t_end,t_sec,
!$OMP& tphi,kwmat_l,ktmat_l,kwmat,ktmat)
!$OMP& SHARED(
!$OMP& rzphi,ion_flag,electron_flag,debug_omp,verbose,cr,idone,
!$OMP& nn,nl,ft,zi,mi,wdfac,divxfac,mpsi,kwmatls,ktmatls,
c!!!!!!... these can be SHARED because all calls use spline_*_external
c!!!!!!...    to avoid conflicts
!$OMP& dcon_int_smats,dcon_int_tmats,dcon_int_xmats,dcon_int_ymats,
!$OMP& dcon_int_zmats)

!$OMP& COPYIN(
c!!!!!!...from inputs.f90...
!$OMP& dbob_m,divx_m,kin,xs_m,fnml,
c!!!!!!...from dcon_interface.f90
!$OMP& geom, dcon_int_sq,dcon_int_eqfun,dcon_int_rzphi)

         DO ipsi=0,mpsi
            DO l=-nl,nl
               CALL SYSTEM_CLOCK(COUNT=t_start)
#ifdef _OPENMP
               IF(ipsi==0 .AND. l==-nl)THEN
                  num_threads = OMP_GET_NUM_THREADS()
                  WRITE(*,'(1x,a,i3,a)') "Running in parallel with ",
     $                 num_threads," OMP threads"
               ENDIF
#endif
               psifac=rzphi%xs(ipsi)
               IF (ion_flag) THEN
                  kwmat_l = 0
                  ktmat_l = 0
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,
     $                 .FALSE.,ft//"wmm",op_wmats=kwmat_l)
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,
     $                 .FALSE.,ft//"tmm",op_wmats=ktmat_l)
                  kwmatls(:,:,:,ipsi,l) = kwmatls(:,:,:,ipsi,l)+kwmat_l
                  ktmatls(:,:,:,ipsi,l) = ktmatls(:,:,:,ipsi,l)+ktmat_l
               ENDIF
               IF (electron_flag) THEN
                  kwmat_l = 0
                  ktmat_l = 0
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,.TRUE.,
     $                 ft//"wmm",op_wmats=kwmat_l)
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,.TRUE.,
     $                 ft//"tmm",op_wmats=ktmat_l)
                  kwmatls(:,:,:,ipsi,l) = kwmatls(:,:,:,ipsi,l)+kwmat_l
                  ktmatls(:,:,:,ipsi,l) = ktmatls(:,:,:,ipsi,l)+ktmat_l
               ENDIF
            ! helpful status messages
!$OMP CRITICAL
            idone = idone + 1
            IF(verbose) CALL progressbar(idone,1,(mpsi+1)*(2*nl+1),
     $         op_percent=10)
!$OMP END CRITICAL
            IF(debug_omp)THEN
               CALL SYSTEM_CLOCK(COUNT=t_end)
               t_sec = REAL(t_end-t_start,r8)/REAL(cr,r8)
               WRITE(*,'(1x,a5,i3,a6,es16.8,a8)')
     $            "ipsi=",ipsi," took ",t_sec," seconds"
            ENDIF
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

c-----------------------------------------------------------------------
c     Now go back totypical serial psi operations
c-----------------------------------------------------------------------
         DO ipsi=0,mpsi
            psifac=rzphi%xs(ipsi)

            ! collect all the ell for a given psi from the collapsed parallel calcs
            kwmat = 0
            ktmat = 0
            DO l=-nl,nl
               kwmat = kwmat + kwmatls(:,:,:,ipsi,l)
               ktmat = ktmat + ktmatls(:,:,:,ipsi,l)
            ENDDO

            ! apply normalizations and hypertangent smoothing for core
            IF (ktanh_flag) THEN
               kwmat=kinfac1*kwmat*(1+tanh((psifac-ktc)*ktw))
               ktmat=kinfac2*ktmat*(1+tanh((psifac-ktc)*ktw))
            ELSE
               kwmat=kinfac1*kwmat
               ktmat=kinfac2*ktmat
            ENDIF
            ! store to splines
            DO i=1,6
               kwmats(i)%fs(ipsi,:)=RESHAPE(kwmat(:,:,i),(/mpert**2/))
               ktmats(i)%fs(ipsi,:)=RESHAPE(ktmat(:,:,i),(/mpert**2/))
            ENDDO
c-----------------------------------------------------------------------
c     Pass only essential matrices to splines
c-----------------------------------------------------------------------
            ix = ipsi
            CALL cspline_eval_external(amats,psifac,ix,amats_f)
            CALL cspline_eval_external(bmats,psifac,ix,bmats_f)
            CALL cspline_eval_external(cmats,psifac,ix,cmats_f)
            CALL cspline_eval_external(dmats,psifac,ix,dmats_f)
            CALL cspline_eval_external(emats,psifac,ix,emats_f)
            CALL cspline_eval_external(hmats,psifac,ix,hmats_f)
            CALL cspline_eval_external(dbats,psifac,ix,dbats_f)
            CALL cspline_eval_external(ebats,psifac,ix,ebats_f)
            CALL cspline_eval_external(fbats,psifac,ix,fbats_f)
            
            amat=RESHAPE(amats_f,(/mpert,mpert/))
            bmat=RESHAPE(bmats_f,(/mpert,mpert/))
            cmat=RESHAPE(cmats_f,(/mpert,mpert/))
            dmat=RESHAPE(dmats_f,(/mpert,mpert/))
            emat=RESHAPE(emats_f,(/mpert,mpert/))
            hmat=RESHAPE(hmats_f,(/mpert,mpert/))
            dbat=RESHAPE(dbats_f,(/mpert,mpert/))
            ebat=RESHAPE(ebats_f,(/mpert,mpert/))
            fmat=RESHAPE(fbats_f,(/mpert,mpert/))
               
            amat=amat+kwmat(:,:,1)+ktmat(:,:,1)
            bmat=bmat+kwmat(:,:,2)+ktmat(:,:,2)
            cmat=cmat+kwmat(:,:,3)+ktmat(:,:,3)
            dmat=dmat+kwmat(:,:,4)+ktmat(:,:,4)
            emat=emat+kwmat(:,:,5)+ktmat(:,:,5)
            hmat=hmat+kwmat(:,:,6)+ktmat(:,:,6)
            caat=cmat-2*ktmat(:,:,3)
            b1mat=ifac*dbat

            ! invert non-hermitian a matrix
            amatlu=0
            umat=0
            DO jpert=1,mpert
               DO ipert=1,mpert
                  amatlu(2*mband+1+ipert-jpert,jpert)=
     $                 amat(ipert,jpert)
                  IF(ipert==jpert)umat(ipert,jpert)=1
               ENDDO
            ENDDO
            CALL zgbtrf(mpert,mpert,mband,mband,amatlu,3*mband+1,
     $           ipiv,info)
            IF(info /= 0)THEN
               WRITE(message,'(2(a,i2))')
     $              "zgbtrf: amat singular at ipsi = ",ipsi,
     $              ", ipert = ",info,", increase delta_mband"
               CALL program_stop(message)
            ENDIF
            
            temp1=dbat 
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp1,mpert,info)
            f0mat=fmat-MATMUL(CONJG(TRANSPOSE(dbat)),temp1)

            ! calculate 3 submatrices for kinetic f matrix
            temp2=amat
            CALL zgbtrs("C",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info) ! close to unit matrix.
            aamat=CONJG(TRANSPOSE(temp2))
            umat=umat-aamat
            
            bkmat=kwmat(:,:,2)+ktmat(:,:,2)+ifac*chi1/(twopi*nn)*
     $           (kwmat(:,:,1)+ktmat(:,:,1))
            bkaat=kwmat(:,:,2)-ktmat(:,:,2)+ifac*chi1/(twopi*nn)*
     $           (kwmat(:,:,1)+ktmat(:,:,1))
            temp2=bkmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            pmat=MATMUL(CONJG(TRANSPOSE(b1mat)),temp2)
               
            temp2=b1mat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            paat=MATMUL(CONJG(TRANSPOSE(bkaat)),temp2)
     $           -ifac*chi1/(twopi*nn)*MATMUL(umat,b1mat)
            paat=CONJG(TRANSPOSE(paat))
            
            temp1=kwmat(:,:,1)+ktmat(:,:,1)
            temp2=bkmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            r1mat=kwmat(:,:,4)+ktmat(:,:,4)-
     $           (chi1/(twopi*nn))**2*CONJG(TRANSPOSE(temp1))+
     $           ifac*chi1/(twopi*nn)*CONJG(TRANSPOSE(bkaat))-
     $           ifac*chi1/(twopi*nn)*MATMUL(aamat,bkmat)-
     $           MATMUL(CONJG(TRANSPOSE(bkaat)),temp2) 
            
            ! calculate 4 submatrices for kinetic k matrix
            temp1=cmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp1,mpert,info)
            kkmat=ebat-MATMUL(CONJG(TRANSPOSE(b1mat)),temp1)
            
            temp1=b1mat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp1,mpert,info)
            kkaat=CONJG(TRANSPOSE(ebat))-
     $           MATMUL(CONJG(TRANSPOSE(caat)),temp1)
            
            temp1=kwmat(:,:,5)+ktmat(:,:,5)-ifac*chi1/(twopi*nn)*
     $           (kwmat(:,:,3)+ktmat(:,:,3))
            temp2=cmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            r2mat=temp1+ifac*chi1/(twopi*nn)*MATMUL(umat,cmat)-
     $           MATMUL(CONJG(TRANSPOSE(bkaat)),temp2)
            
            temp1=kwmat(:,:,5)-ktmat(:,:,5)-ifac*chi1/(twopi*nn)*
     $           (kwmat(:,:,3)-ktmat(:,:,3))
            temp2=bkmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            r3mat=CONJG(TRANSPOSE(temp1))-
     $           MATMUL(CONJG(TRANSPOSE(caat)),temp2)

            ! calculate kinetic g matrix
            temp2=cmat
            CALL zgbtrs("N",mpert,mband,mband,mpert,amatlu,
     $           3*mband+1,ipiv,temp2,mpert,info)
            gaat=hmat-MATMUL(CONJG(TRANSPOSE(caat)),temp2) 

            ! pass only essential kinetic matrices
            akmats%fs(ipsi,:)=RESHAPE(amat,(/mpert**2/))
            bkmats%fs(ipsi,:)=RESHAPE(bmat,(/mpert**2/))
            ckmats%fs(ipsi,:)=RESHAPE(cmat,(/mpert**2/))
            f0mats%fs(ipsi,:)=RESHAPE(f0mat,(/mpert**2/))
            pmats%fs(ipsi,:)=RESHAPE(pmat,(/mpert**2/))
            paats%fs(ipsi,:)=RESHAPE(paat,(/mpert**2/))
            kkmats%fs(ipsi,:)=RESHAPE(kkmat,(/mpert**2/))
            kkaats%fs(ipsi,:)=RESHAPE(kkaat,(/mpert**2/))
            r1mats%fs(ipsi,:)=RESHAPE(r1mat,(/mpert**2/))
            r2mats%fs(ipsi,:)=RESHAPE(r2mat,(/mpert**2/))
            r3mats%fs(ipsi,:)=RESHAPE(r3mat,(/mpert**2/))

            ! pass banded g matrix directly
            iqty=1
            DO jpert=1,mpert
               DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
                  gaats%fs(ipsi,iqty)=gaat(ipert,jpert)
                  iqty=iqty+1
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(kwmatls, ktmatls)  ! these are really big, get rid of them asap
c-----------------------------------------------------------------------
c     fit splines
c-----------------------------------------------------------------------
         DO i=1,6
            CALL cspline_fit(kwmats(i),"extrap")
            CALL cspline_fit(ktmats(i),"extrap")
         ENDDO

         CALL cspline_fit(akmats,"extrap")
         CALL cspline_fit(bkmats,"extrap")
         CALL cspline_fit(ckmats,"extrap")
         CALL cspline_fit(f0mats,"extrap")
         CALL cspline_fit(pmats,"extrap")
         CALL cspline_fit(paats,"extrap")
         CALL cspline_fit(kkmats,"extrap")
         CALL cspline_fit(kkaats,"extrap")
         CALL cspline_fit(r1mats,"extrap")
         CALL cspline_fit(r2mats,"extrap")
         CALL cspline_fit(r3mats,"extrap")         
c-----------------------------------------------------------------------
c     Use built in PENTRC spline integration options to form matrixes
c-----------------------------------------------------------------------
      ELSEIF(method==1)THEN
         DO i=1,6
            CALL cspline_alloc(kwmats(i),mpsi,mpert**2)
            CALL cspline_alloc(ktmats(i),mpsi,mpert**2)
            kwmats(i)%xs=rzphi%xs
            ktmats(i)%xs=rzphi%xs
         ENDDO         
         DO ipsi=0,mpsi
            iindex = FLOOR(REAL(ipsi+1,r8)/FLOOR((mpsi+1)/10.0))*10
            ileft = REAL(ipsi+1,r8)/FLOOR((mpsi+1)/10.0)*10-iindex
            IF ((ipsi /= 0) .AND. (ileft == 0) .AND. verbose)
     $           WRITE(*,*)"  ...",iindex,"% of kinetic computations"
            kwmat = 0
            ktmat = 0
            psifac=rzphi%xs(ipsi)
            ! get matrices for all ell at this one psi
            IF (ion_flag) THEN
               DO l=-nl,nl
                  kwmat_l = 0
                  ktmat_l = 0
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,
     $                 .FALSE.,ft//"wmm",op_wmats=kwmat_l)
                  kwmat = kwmat+kwmat_l
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,
     $                 .FALSE.,ft//"tmm",op_wmats=ktmat_l)
                  ktmat = ktmat+ktmat_l
               ENDDO
            ENDIF
            IF (electron_flag) THEN
               DO l=-nl,nl
                  kwmat_l = 0
                  ktmat_l = 0
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,.TRUE.,
     $                 ft//"wmm",op_wmats=kwmat_l)
                  kwmat = kwmat+kwmat_l
                  call tpsi(tphi,psifac,nn,l,zi,mi,wdfac,divxfac,.TRUE.,
     $                 ft//"tmm",op_wmats=ktmat_l)
                  ktmat = ktmat+ktmat_l
               ENDDO
            ENDIF
            ! apply normalizations and hypertangent smoothing for core
            IF (ktanh_flag) THEN
               kwmat=kinfac1*kwmat*(1+tanh((psifac-ktc)*ktw))
               ktmat=kinfac2*ktmat*(1+tanh((psifac-ktc)*ktw))
            ELSE
               kwmat=kinfac1*kwmat
               ktmat=kinfac2*ktmat
            ENDIF
            ! store to splines
            DO i=1,6
               kwmats(i)%fs(ipsi,:)=RESHAPE(kwmat(:,:,i),(/mpert**2/))
               ktmats(i)%fs(ipsi,:)=RESHAPE(ktmat(:,:,i),(/mpert**2/))
            ENDDO
         ENDDO
         DO i=1,6
            CALL cspline_fit(kwmats(i),"extrap")
            CALL cspline_fit(ktmats(i),"extrap")
         ENDDO
c-----------------------------------------------------------------------
c     Use built in PENTRC spline integration options to form matrixes
c-----------------------------------------------------------------------
      ELSEIF(method==2)THEN
         WRITE(*,*) " Kinetic energy calculation using MXM euler "//
     $        "lagrange matrix on equilibrium grid"
         tphi = tintgrl_grid('equil',plim,nn,nl,zi,mi,wdfac,divxfac,
     $        .FALSE.,ft//"wmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),kwmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  kwmats(i)%fs(ipsi,:) = kinfac1*kwmats(i)%fs(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
                  kwmats(i)%fs1(ipsi,:) = kinfac1*kwmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               kwmats(i)%fs = kinfac1*kwmats(i)%fs
               kwmats(i)%fs1 = kinfac1*kwmats(i)%fs1
            ENDIF
         ENDDO
         WRITE(*,*) " Kinetic torque calculation using MXM euler "//
     $      "lagrange matrix on equilibrium grid"
         tphi = tintgrl_grid('equil',plim,nn,nl,zi,mi,wdfac,divxfac,
     $        .FALSE.,ft//"tmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),ktmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  ktmats(i)%fs(ipsi,:) = kinfac2*ktmats(i)%fs(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
                  ktmats(i)%fs1(ipsi,:) = kinfac2*ktmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               ktmats(i)%fs = kinfac2*ktmats(i)%fs
               ktmats(i)%fs1 = kinfac2*ktmats(i)%fs1
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     Use built in PENTRC LSODE integration options to form matrixes
c      -> Grid determined by T & dW from flat xi spectrum
c-----------------------------------------------------------------------
      ELSEIF(method==3)THEN
         WRITE(*,*) " Kinetic energy calculation using MXM euler "//
     $      "lagrange matrix"
         tphi = tintgrl_lsode(plim,nn,nl,zi,mi,wdfac,divxfac,.FALSE.,
     $        ft//"wmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),kwmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  kwmats(i)%fs(ipsi,:) = kinfac1*kwmats(i)%fs(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
                  kwmats(i)%fs1(ipsi,:) = kinfac1*kwmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               kwmats(i)%fs = kinfac1*kwmats(i)%fs
               kwmats(i)%fs1 = kinfac1*kwmats(i)%fs1
            ENDIF
         ENDDO
         WRITE(*,*) " Kinetic torque calculation using MXM euler "//
     $      "lagrange matrix"
         tphi = tintgrl_lsode(plim,nn,nl,zi,mi,wdfac,divxfac,.FALSE.,
     $        ft//"tmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),ktmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  ktmats(i)%fs(ipsi,:) = kinfac2*ktmats(i)%fs(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
                  ktmats(i)%fs1(ipsi,:) = kinfac2*ktmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               ktmats(i)%fs = kinfac2*ktmats(i)%fs
               ktmats(i)%fs1 = kinfac2*ktmats(i)%fs1
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     Use built in PENTRC LSODE integration options to form matrixes
c      -> Grid determined by norm of EL matrices
c-----------------------------------------------------------------------
      ELSEIF(method==4)THEN
         WRITE(*,*) " Kinetic MXM euler lagrange energy matrix norm "
     $      //"calculation"
         tphi = tintgrl_lsode(plim,nn,nl,zi,mi,wdfac,divxfac,.FALSE.,
     $        ft//"kmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),kwmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  kwmats(i)%fs(ipsi,:) = kinfac1*kwmats(i)%fs(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
                  kwmats(i)%fs1(ipsi,:) = kinfac1*kwmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((kwmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               kwmats(i)%fs = kinfac1*kwmats(i)%fs
               kwmats(i)%fs1 = kinfac1*kwmats(i)%fs1
            ENDIF
         ENDDO
         WRITE(*,*) " Kinetic MXM euler lagrange torque matrix norm "
     $      //"calculation"
         tphi = tintgrl_lsode(plim,nn,nl,zi,mi,wdfac,divxfac,.FALSE.,
     $        ft//"rmm")
         ! copy and apply factor to splines
         DO i=1,6
            CALL cspline_copy(kelmm(i),ktmats(i))
            IF (ktanh_flag) THEN
               DO ipsi=0,kwmats(1)%mx
                  ktmats(i)%fs(ipsi,:) = kinfac2*ktmats(i)%fs(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
                  ktmats(i)%fs1(ipsi,:) = kinfac2*ktmats(i)%fs1(ipsi,:)*
     $                 (1+tanh((ktmats(i)%xs(ipsi)-ktc)*ktw))
               ENDDO
            ELSE
               ktmats(i)%fs = kinfac2*ktmats(i)%fs
               ktmats(i)%fs1 = kinfac2*ktmats(i)%fs1
            ENDIF
         ENDDO
      ELSE
         CALL program_stop("ERROR: Valid kingridtypes are 0,1,2,3,4")
      ENDIF
      gaats%x0(2)=1.0
      gaats%xpower(1,:)=-1
      gaats%xpower(2,:)=-1
      CALL cspline_fit(kaats,"extrap")
      CALL cspline_fit(gaats,"extrap")
c-----------------------------------------------------------------------
c     Optionally write matrices to binary files for diagnostics
c-----------------------------------------------------------------------
      IF(output)THEN
         ! binary output
         CALL bin_open(bin_unit,"kwmats.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert**2
            DO ipsi=0,kwmats(1)%mx
               WRITE(bin_unit)   REAL(kwmats(1)%xs(ipsi),4),
     $              REAL(REAL(kwmats(1)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(1)%fs(ipsi,ipert)),4),
     $              REAL(REAL(kwmats(2)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(2)%fs(ipsi,ipert)),4),
     $              REAL(REAL(kwmats(3)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(3)%fs(ipsi,ipert)),4),
     $              REAL(REAL(kwmats(4)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(4)%fs(ipsi,ipert)),4),
     $              REAL(REAL(kwmats(5)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(5)%fs(ipsi,ipert)),4),
     $              REAL(REAL(kwmats(6)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(kwmats(6)%fs(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"ktmats.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert**2
            DO ipsi=0,ktmats(1)%mx
               WRITE(bin_unit)   REAL(ktmats(1)%xs(ipsi),4),
     $              REAL(REAL(ktmats(1)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(1)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(2)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(2)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(3)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(3)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(4)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(4)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(5)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(5)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(6)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(6)%fs(ipsi,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         
         mfac =(/(i,i=mlow,mhigh)/)
         CALL ascii_open(fourfit_out_unit,"kwmats.out","UNKNOWN")
         WRITE(fourfit_out_unit,*)"DCON Kinetic energy matrices"
         WRITE(fourfit_out_unit,'(1/,1x,a12,1x,I6,1x,1(a12,I4),1/)')
     $        "mpsi =",mpsi,"mpert =",mpert
         WRITE(fourfit_out_unit,'(1x,a16,2(1x,a4),12(1x,a16))')
     $      "psi","m1","m2",
     $      "real(Ak)","imag(Ak)","real(Bk)","imag(Bk)",
     $      "real(Ck)","imag(Ck)","real(Dk)","imag(Dk)",
     $      "real(Ek)","imag(Ek)","real(Hk)","imag(Hk)"
         DO ipsi=0,kwmats(1)%mx 
            DO i=1,mpert
               DO j=1,mpert
                  ipert = (i-1)*mpert + j
                  WRITE(fourfit_out_unit,'(1x,es16.8,2(1x,I4),'//
     $               '12(1x,es16.8))')
     $               REAL(kwmats(1)%xs(ipsi),4),mfac(i),mfac(j),
     $               REAL(REAL(kwmats(1)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(1)%fs(ipsi,ipert)),4),
     $               REAL(REAL(kwmats(2)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(2)%fs(ipsi,ipert)),4),
     $               REAL(REAL(kwmats(3)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(3)%fs(ipsi,ipert)),4),
     $               REAL(REAL(kwmats(4)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(4)%fs(ipsi,ipert)),4),
     $               REAL(REAL(kwmats(5)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(5)%fs(ipsi,ipert)),4),
     $               REAL(REAL(kwmats(6)%fs(ipsi,ipert)),4),
     $               REAL(AIMAG(kwmats(6)%fs(ipsi,ipert)),4)
               ENDDO
            ENDDO
         ENDDO
         WRITE(fourfit_out_unit,*)
         CALL ascii_close(fourfit_out_unit)
         CALL ascii_open(fourfit_out_unit,"ktmats.out","UNKNOWN")
         WRITE(fourfit_out_unit,*) "DCON Kinetic torque matrices"
         WRITE(fourfit_out_unit,'(1/,1x,a12,1x,I6,1x,1(a12,I4),1/)')
     $      "mpsi =",mpsi,"mpert =",mpert
         WRITE(fourfit_out_unit,'(1x,a16,2(1x,a4),12(1x,a16))')
     $      "psi","m1","m2",
     $      "real(Ak)","imag(Ak)","real(Bk)","imag(Bk)",
     $      "real(Ck)","imag(Ck)","real(Dk)","imag(Dk)",
     $      "real(Ek)","imag(Ek)","real(Hk)","imag(Hk)"
         DO ipsi=0,kwmats(1)%mx 
            DO i=1,mpert
               DO j=1,mpert
                  ipert = (i-1)*mpert + j
                  WRITE(fourfit_out_unit,'(1x,es16.8,2(1x,I4),'//
     $              '12(1x,es16.8))')
     $              REAL(kwmats(1)%xs(ipsi),4),mfac(i),mfac(j),
     $              REAL(REAL(ktmats(1)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(1)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(2)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(2)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(3)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(3)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(4)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(4)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(5)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(5)%fs(ipsi,ipert)),4),
     $              REAL(REAL(ktmats(6)%fs(ipsi,ipert)),4),
     $              REAL(AIMAG(ktmats(6)%fs(ipsi,ipert)),4)
               ENDDO
            ENDDO
         ENDDO
         WRITE(fourfit_out_unit,*)
         CALL ascii_close(fourfit_out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_kinetic_matrix
      END MODULE fourfit_mod
