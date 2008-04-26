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
c     3. fourfit_write_metric.
c     4. fourfit_write_matrix.
c     5. fourfit_evals.
c     6. fourfit_diagnose_1.
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
      IMPLICIT NONE

      TYPE(fspline_type), PRIVATE :: metric
      TYPE(cspline_type) :: fmats,gmats,kmats
      TYPE(spline_type) :: k0s
      INTEGER, DIMENSION(:), POINTER :: ipiva
      COMPLEX(r8), DIMENSION(:,:), POINTER :: amat,bmat,cmat
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
      REAL(r8) :: theta,rfac,eta,r,jac,jac1,psifac
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
c     begin loop over nodes.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            theta=rzphi%ys(itheta)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta+rzphi%f(2))
            r=ro+rfac*COS(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
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
c-----------------------------------------------------------------------
c     compute metric tensor components.
c-----------------------------------------------------------------------
            jac=jac
            metric%fs(ipsi,itheta,1)=SUM(v(1,:)**2)*jac
            metric%fs(ipsi,itheta,2)=SUM(v(2,:)**2)*jac
            metric%fs(ipsi,itheta,3)=v(3,3)*v(3,3)*jac
            metric%fs(ipsi,itheta,4)=v(2,3)*v(3,3)*jac
            metric%fs(ipsi,itheta,5)=v(3,3)*v(1,3)*jac
            metric%fs(ipsi,itheta,6)=SUM(v(1,:)*v(2,:))*jac
            metric%fs(ipsi,itheta,7)=jac
            metric%fs(ipsi,itheta,8)=jac1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit Fourier-spline type.
c-----------------------------------------------------------------------
      IF(fft_flag)THEN
         CALL fspline_fit_2(metric,"extrap",.FALSE.)
      ELSE
         CALL fspline_fit_1(metric,"extrap",.FALSE.)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_make_metric
c-----------------------------------------------------------------------
c     subprogram 2. fourfit_make_matrix.
c     constructs the coefficient matrices and fits them to cubic splines
c     .
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfit_make_matrix

      CHARACTER(128) :: message
      INTEGER :: ipsi,ipert,jpert,m1,m2,m,dm,info,iqty
      REAL(r8) :: chi1,jtheta,nq,p1,psifac,q,q1,singfac1,singfac2
      COMPLEX(r8), DIMENSION(mpert*mpert) :: work
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb
      COMPLEX(r8), DIMENSION(-mband:mband) ::
     $     g11,g22,g33,g23,g31,g12,jmat1,imat
      COMPLEX(r8), DIMENSION(mpert,mpert) :: dmat,emat,
     $     fmat,gmat,hmat,kmat,temp1,temp2
      TYPE(cspline_type) :: amats,bmats,cmats

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER, PARAMETER :: unit=99
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
      ALLOCATE(amat(mpert,mpert),bmat(mpert,mpert),cmat(mpert,mpert),
     $     ipiva(mpert),jmat(-mband:mband))
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
      IF(sas_flag)THEN
         CALL cspline_alloc(amats,mpsi,mpert**2)
         CALL cspline_alloc(bmats,mpsi,mpert**2)
         CALL cspline_alloc(cmats,mpsi,mpert**2)
         amats%xs=sq%xs
         bmats%xs=sq%xs
         cmats%xs=sq%xs
      ENDIF
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
c     CALL fourfit_diagnose_1(g11,g22,g33,g23,g31,g12)
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
     $              +jtheta*imat(dm)
     $              )
               hmat(ipert,jpert)=(q1*chi1)**2*g33(dm)
     $              +(twopi*chi1)**2*singfac1*singfac2*g11(dm)
     $              -twopi*ifac*chi1*dm*q1*chi1*g31(dm)
     $              +jtheta*q1*chi1*imat(dm)+p1*jmat1(dm)
               fmat(ipert,jpert)=(chi1/nn)**2*g33(dm)
               kmat(ipert,jpert)=twopi*ifac*chi1*(g23(dm)+g33(dm)*m1/nn)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     store matrices for interpolation to psilim.
c-----------------------------------------------------------------------
         IF(sas_flag)THEN
            amats%fs(ipsi,:)=RESHAPE(amat,(/mpert**2/))
            bmats%fs(ipsi,:)=RESHAPE(bmat,(/mpert**2/))
            cmats%fs(ipsi,:)=RESHAPE(cmat,(/mpert**2/))
         ENDIF
c-----------------------------------------------------------------------
c     factor A.
c-----------------------------------------------------------------------
         CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
         IF(info /= 0)THEN
            WRITE(message,'(2(a,i2))')
     $           "zhetrf: amat singular at ipsi = ",ipsi,
     $           ", ipert = ",info,", reduce delta_mband"
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
      CALL cspline_fit(fmats,"extrap")
      CALL cspline_fit(kmats,"extrap")
      CALL cspline_fit(gmats,"extrap")
c-----------------------------------------------------------------------
c     interpolate matrices to psilim.
c-----------------------------------------------------------------------
      IF(sas_flag)THEN
         CALL cspline_fit(amats,"extrap")
         CALL cspline_fit(bmats,"extrap")
         CALL cspline_fit(cmats,"extrap")
         CALL cspline_eval(amats,psilim,0)
         CALL cspline_eval(bmats,psilim,0)
         CALL cspline_eval(cmats,psilim,0)
         amat=RESHAPE(amats%f,(/mpert,mpert/))
         bmat=RESHAPE(bmats%f,(/mpert,mpert/))
         cmat=RESHAPE(cmats%f,(/mpert,mpert/))
         CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
         CALL cspline_dealloc(amats)
         CALL cspline_dealloc(bmats)
         CALL cspline_dealloc(cmats)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose and deallocate.
c-----------------------------------------------------------------------
      IF(bin_metric)CALL fourfit_write_metric
      IF(bin_fmat)CALL fourfit_write_matrix(fmats,"fmat",.TRUE.)
      IF(bin_gmat)CALL fourfit_write_matrix(gmats,"gmat",.TRUE.)
      IF(bin_kmat)CALL fourfit_write_matrix(kmats,"kmat",.FALSE.)
      CALL fspline_dealloc(metric)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfit_make_matrix
c-----------------------------------------------------------------------
c     subprogram 3. fourfit_write_metric.
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
c     subprogram 4. fourfit_write_matrix
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
c     subprogram 5. fourfit_evals
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
c     subprogram 6. fourfit_diagnose_1.
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
