c-----------------------------------------------------------------------
c     file fspline.f.
c     fits functions to cubic spline in x and Fourier series in y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fspline_mod.
c     1. fspline_alloc.
c     2. fspline_dealloc.
c     3. fspline_fit_1.
c     4. fspline_fit_2.
c     5. fspline_eval.
c     6. fspline_all_eval.
c     7. fspline_write_xy.
c     8. fspline_write_yx.
c     9. fspline_copy.
c-----------------------------------------------------------------------
c     subprogram 0. fspline_type definition.
c     defines fspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fspline_mod
      USE spline_mod
      USE cspline_mod
      USE utils_mod
      USE fft_mod
      IMPLICIT NONE

      TYPE :: fspline_type
      INTEGER :: mx,my,mband,nqty
      REAL(r8), DIMENSION(:), POINTER :: xs,ys
      REAL(r8), DIMENSION(:,:,:), POINTER :: fs
      TYPE(cspline_type) :: cs
      REAL(r8), DIMENSION(:), POINTER :: f,fx,fy,fxx,fxy,fyy
      CHARACTER(6) :: xtitle,ytitle
      CHARACTER(6), DIMENSION(:), POINTER :: title
      CHARACTER(6) :: name
      END TYPE fspline_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fspline_alloc.
c     allocates space for fspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_alloc(fst,mx,my,mband,nqty)

      INTEGER, INTENT(IN) :: mx,my,mband,nqty
      TYPE(fspline_type), INTENT(OUT) :: fst
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      fst%mx=mx
      fst%my=my
      fst%mband=mband
      fst%nqty=nqty
      ALLOCATE(fst%xs(0:mx))
      ALLOCATE(fst%ys(0:my))
      ALLOCATE(fst%fs(0:mx,0:my,nqty))
      ALLOCATE(fst%title(nqty))
      ALLOCATE(fst%f(nqty))
      ALLOCATE(fst%fx(nqty))
      ALLOCATE(fst%fy(nqty))
      ALLOCATE(fst%fxx(nqty))
      ALLOCATE(fst%fxy(nqty))
      ALLOCATE(fst%fyy(nqty))
      CALL cspline_alloc(fst%cs,mx,(mband+1)*nqty)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_alloc
c-----------------------------------------------------------------------
c     subprogram 2. fspline_dealloc.
c     deallocates space for fspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_dealloc(fst)

      TYPE(fspline_type), INTENT(INOUT) :: fst
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(fst%xs)
      DEALLOCATE(fst%ys)
      DEALLOCATE(fst%fs)
      DEALLOCATE(fst%title)
      DEALLOCATE(fst%f)
      DEALLOCATE(fst%fx)
      DEALLOCATE(fst%fy)
      DEALLOCATE(fst%fxx)
      DEALLOCATE(fst%fxy)
      DEALLOCATE(fst%fyy)
      CALL cspline_dealloc(fst%cs)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. fspline_fit_1.
c     fits functions to fsplines by integrating periodic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_fit_1(fst,endmode,fit_flag)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      CHARACTER(*), INTENT(IN) :: endmode
      LOGICAL, INTENT(IN) :: fit_flag

      INTEGER :: m,ix,iq,mx,my,mband,nqty,j
      REAL(r8), PARAMETER :: eps=1e-3
      REAL(r8), DIMENSION(fst%my) :: delta,d4fac
      REAL(r8), DIMENSION(0:fst%mx,0:fst%my,fst%nqty) :: fs,fsy
      COMPLEX(r8), DIMENSION(0:fst%my) :: expfac0,expfac
      COMPLEX(r8), DIMENSION(fst%my) :: dexpfac0,dexpfac
      COMPLEX(r8), DIMENSION(fst%my) :: alpha1,alpha2,beta1,beta2
      COMPLEX(r8), DIMENSION(0:fst%mx,0:fst%mband,fst%nqty) :: coef
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     copy sizes and zero Fourier coefficients.
c-----------------------------------------------------------------------
      mx=fst%mx
      my=fst%my
      mband=fst%mband
      nqty=fst%nqty
      coef=0
c-----------------------------------------------------------------------
c     prepare principal periodic functions.
c-----------------------------------------------------------------------
      delta=fst%ys(1:my)-fst%ys(0:my-1)
      d4fac=1/delta**4
      expfac0=EXP(-ifac*fst%ys)
      dexpfac0=expfac0(0:my-1)/expfac0(1:my)
      expfac=1
      dexpfac=1
c-----------------------------------------------------------------------
c     compute y derivatives.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,my,mx+1)
      spl%xs=fst%ys
      DO iq=1,nqty
         spl%fs=TRANSPOSE(fst%fs(:,:,iq))
         CALL spline_fit(spl,"periodic")
         fsy(:,:,iq)=TRANSPOSE(spl%fs1)
      ENDDO
      fs=fst%fs
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     compute alpha's and beta's.
c-----------------------------------------------------------------------
      DO m=0,mband
         WHERE(ABS(m*delta) > eps)
            alpha1=(dexpfac*(12._r8-ifac*m*delta*(6._r8+(m*delta)**2))
     $           -6._r8*(2._r8+ifac*m*delta))*d4fac/m**4
            beta1=(dexpfac*(6._r8-m*delta*(4*ifac+m*delta))
     $           -2._r8*(3._r8+ifac*m*delta))*d4fac/m**4
         ELSEWHERE
            alpha1=.5_r8+7._r8*ifac*m*delta/20._r8-2._r8*
     $           (m*delta)**2/15._r8
            beta1=1._r8/12._r8+ifac*m*delta/20._r8-(m*delta)**2/60._r8
         ENDWHERE
         alpha1=alpha1*delta/twopi
         beta1=beta1*delta**2/twopi
         alpha2=CONJG(alpha1)*expfac(0:mx-1)
         beta2=CONJG(beta1)*expfac(0:mx-1)
         alpha1=alpha1*expfac(1:mx)
         beta1=beta1*expfac(1:mx)
c-----------------------------------------------------------------------
c     compute Fourier coefficients.
c-----------------------------------------------------------------------
         DO iq=1,nqty
            DO ix=0,mx
               coef(ix,m,iq)=coef(ix,m,iq)+SUM(
     $              alpha1*fs(ix,0:my-1,iq)+alpha2*fs(ix,1:my,iq)
     $              +beta1*fsy(ix,0:my-1,iq)-beta2*fsy(ix,1:my,iq))
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     advance to next Fourier component.
c-----------------------------------------------------------------------
         expfac=expfac*expfac0
         dexpfac=dexpfac*dexpfac0
      ENDDO
c-----------------------------------------------------------------------
c     fit Fourier coefficients to cubic splines as functions of x.
c-----------------------------------------------------------------------
      fst%cs%xs=fst%xs
      fst%cs%fs=RESHAPE(coef,(/mx+1,(mband+1)*nqty/))
      IF(fit_flag)CALL cspline_fit(fst%cs,endmode)
      fst%cs%name=fst%name
      fst%cs%title(0)="  x   "
      j=0
      DO iq=1,nqty
         DO m=0,mband
            j=j+1
            WRITE(fst%cs%title(j),'("cs",i1,"_",i2.2)')iq,m
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_fit_1
c-----------------------------------------------------------------------
c     subprogram 4. fspline_fit_2.
c     fits functions to fsplines using Fast Fourier Transform.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_fit_2(fst,endmode,fit_flag)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      CHARACTER(*), INTENT(IN) :: endmode
      LOGICAL, INTENT(IN) :: fit_flag

      CHARACTER(64) :: message
      INTEGER :: iqty,j,m,mband,mx,my,nqty,p2
      COMPLEX(r8), DIMENSION(0:fst%my-1,0:fst%mx,fst%nqty) :: f,g
c-----------------------------------------------------------------------
c     abort if my-1 is not a power of 2.
c-----------------------------------------------------------------------
      p2=1
      my=fst%my
      DO
         my=my/2
         IF(my == 0)EXIT
         p2=p2*2
      ENDDO
      IF(fst%my /= p2)THEN
         WRITE(message,'(a,i3,a)')
     $     "fft_fit_2: my = ",fst%my," is not a power of 2"
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     abort if 2*mband > my-1.
c-----------------------------------------------------------------------
      IF(2*fst%mband > fst%my-1)THEN
         WRITE(message,'(a,i3,a,i3)')
     $        "fft_fit_2: 2*mband = ",2*fst%mband," > my-1 = ",fst%my-1
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     copy sizes and zero Fourier coefficients.
c-----------------------------------------------------------------------
      mx=fst%mx
      my=fst%my
      mband=fst%mband
      nqty=fst%nqty
c-----------------------------------------------------------------------
c     set up and solve Fast Fourier Transform.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         f(:,:,iqty)=TRANSPOSE(fst%fs(:,0:my-1,iqty))
      ENDDO
      g=RESHAPE(fft_run(RESHAPE(f,(/my,(mx+1)*nqty/)),-1),
     $     (/my,mx+1,nqty/))
c-----------------------------------------------------------------------
c     copy Fourier coefficients.
c-----------------------------------------------------------------------
      j=1
      DO iqty=1,nqty
         DO m=0,mband
            fst%cs%fs(:,j)=g(m,:,iqty)
            j=j+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit Fourier coefficients to cubic splines as functions of x.
c-----------------------------------------------------------------------
      fst%cs%xs=fst%xs
      IF(fit_flag)CALL cspline_fit(fst%cs,endmode)
      fst%cs%name=fst%name
      fst%cs%title(0)="  x   "
      j=0
      DO iqty=1,nqty
         DO m=0,mband
            j=j+1
            WRITE(fst%cs%title(j),'("cs",i1,"_",i2.2)')iqty,m
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_fit_2
c-----------------------------------------------------------------------
c     subprogram 5. fspline_eval.
c     evaluates fspline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_eval(fst,x,y,mode)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      REAL(r8), INTENT(IN) :: x,y
      INTEGER, INTENT(IN) :: mode

      INTEGER :: m
      COMPLEX(r8) :: expfac,expfac0
      COMPLEX(r8), DIMENSION(fst%nqty) :: term,termx,termxx
      COMPLEX(r8), DIMENSION(0:fst%mband,fst%nqty) :: c,cx,cxx
c-----------------------------------------------------------------------
c     evaluate cubic splines and m = 0 terms for functions.
c-----------------------------------------------------------------------
      CALL cspline_eval(fst%cs,x,mode)
      c=RESHAPE(fst%cs%f,(/fst%mband+1,fst%nqty/))
      fst%f=c(0,:)
c-----------------------------------------------------------------------
c     evaluate m = 0 terms for first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         cx=RESHAPE(fst%cs%f1,(/fst%mband+1,fst%nqty/))
         fst%fx=cx(0,:)
         fst%fy=0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate m = 0 terms for second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         cxx=RESHAPE(fst%cs%f2,(/fst%mband+1,fst%nqty/))
         fst%fxx=cxx(0,:)
         fst%fyy=0
         fst%fxy=0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate m > 1 for functions.
c-----------------------------------------------------------------------
      expfac0=EXP(ifac*y)
      expfac=2
      DO m=1,fst%mband
         expfac=expfac*expfac0
         term=c(m,:)*expfac
         fst%f=fst%f+term
         IF(mode < 1)CYCLE
c-----------------------------------------------------------------------
c     evaluate m > 1 for first derivatives.
c-----------------------------------------------------------------------
         termx=cx(m,:)*expfac
         fst%fx=fst%fx+termx
         fst%fy=fst%fy+term*ifac*m
         IF(mode < 2)CYCLE
c-----------------------------------------------------------------------
c     evaluate m > 1 for second derivatives.
c-----------------------------------------------------------------------
         termxx=cxx(m,:)*expfac
         fst%fxx=fst%fxx+termxx
         fst%fxy=fst%fxy+termx*ifac*m
         fst%fyy=fst%fyy-term*m*m
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_eval
c-----------------------------------------------------------------------
c     subprogram 6. fspline_all_eval.
c     evaluates fsplines in all intervals.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_all_eval(fst,dx,dy,f,fx,fy,fxx,fyy,fxy,mode)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      REAL(r8), INTENT(IN) :: dx,dy
      REAL(r8), INTENT(OUT), DIMENSION(fst%mx,fst%my,fst%nqty) ::
     $     f,fx,fy,fxx,fyy,fxy
      INTEGER, INTENT(IN) :: mode

      INTEGER :: m,iy
      REAL(r8), DIMENSION(fst%my) :: y
      COMPLEX(r8), DIMENSION(fst%my) :: expfac,expfac0
      COMPLEX(r8), DIMENSION(fst%my,fst%nqty*(fst%mx+1)) ::
     $     f0,f1,f2,f3
      COMPLEX(r8), DIMENSION(fst%mx,fst%my,fst%nqty) ::
     $     term,termx,termxx
      COMPLEX(r8), DIMENSION(fst%mx,0:fst%mband,fst%nqty) :: c,cx,cxx
c-----------------------------------------------------------------------
c     evaluate cubic splines and m = 0 terms for functions.
c-----------------------------------------------------------------------
      CALL cspline_all_eval(fst%cs,dx,f0,f1,f2,f3,mode)
      c=RESHAPE(f0,(/fst%mx,fst%mband+1,fst%nqty/))
      DO iy=1,fst%my
         f(:,iy,:)=c(:,0,:)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate m = 0 terms for first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         cx=RESHAPE(f1,(/fst%mx,fst%mband+1,fst%nqty/))
         DO iy=1,fst%my
            fx(:,iy,:)=cx(:,0,:)
         ENDDO
         fy=0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate m = 0 terms for second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         cxx=RESHAPE(f2,(/fst%mx,fst%mband+1,fst%nqty/))
         DO iy=1,fst%my
            fxx(:,iy,:)=cxx(:,0,:)
         ENDDO
         fxy=0
         fyy=0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate m > 0 terms for functions.
c-----------------------------------------------------------------------
      y=fst%ys(0:fst%my-1)+dy*(fst%ys(1:fst%my)-fst%ys(0:fst%my-1))
      expfac0=EXP(ifac*y)
      expfac=2
      DO m=1,fst%mband
         expfac=expfac*expfac0
         DO iy=1,fst%my
            term(:,iy,:)=c(:,m,:)*expfac(iy)
         ENDDO
         f=f+term
         IF(mode < 1)CYCLE
c-----------------------------------------------------------------------
c     evaluate m > 0 terms for first derivatives.
c-----------------------------------------------------------------------
         DO iy=1,fst%my
            termx(:,iy,:)=cx(:,m,:)*expfac(iy)
         ENDDO
         fx=fx+termx
         fy=fy+term*ifac*m
         IF(mode < 2)CYCLE
c-----------------------------------------------------------------------
c     evaluate m > 0 terms for second derivatives.
c-----------------------------------------------------------------------
         DO iy=1,fst%my
            termxx(:,iy,:)=cxx(:,m,:)*expfac(iy)
         ENDDO
         fxx=fxx+termxx
         fxy=fxy+termx*ifac*m
         fyy=fyy-term*m*m
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_all_eval
c-----------------------------------------------------------------------
c     subprogram 7. fspline_write_xy.
c     produces ascii and binary output for fspline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_write_xy(fst,out,bin,interpolate,filename)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      LOGICAL, INTENT(IN) :: out,bin,interpolate
      CHARACTER(*) :: filename

      INTEGER :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      CHARACTER(80) :: format2,format1
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/1x,'ix = ',i3,', x = ',1p,e11.3)
 20   FORMAT(/1x,'ix = ',i3,', jx = ',i1,', x = ',1p,e11.3)
 30   FORMAT('(/4x,"iy",6x,"y",4x,',i1,'(6x,"f",i1,3x)/)')
 40   FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     open binary output file.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
      IF(bin)CALL bin_open(bin_unit,TRIM(filename),
     $     "UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(out_unit,'(1x,a)')"input data"
         WRITE(format1,30)fst%nqty
         WRITE(format2,40)fst%nqty
      ENDIF
      DO iy=0,fst%my
         y=fst%ys(iy)
         IF(out)then
            WRITE(out_unit,10)iy,fst%ys(iy)
            WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
         ENDIF
         DO ix=0,fst%mx
            x=fst%xs(ix)
            fst%f=fst%fs(ix,iy,:)
            IF(out)WRITE(out_unit,format2)ix,x,fst%f
            IF(bin)WRITE(bin_unit)REAL(x,4),REAL(fst%f,4)
         ENDDO
         IF(out)WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
         IF(bin)WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
      IF(interpolate)THEN
         IF(out)WRITE(out_unit,'(1x,a)')"interpolated data"
         DO iy=0,fst%my-1
            dy=(fst%ys(iy+1)-fst%ys(iy))/4
            DO jy=0,4
               y=fst%ys(iy)+dy*jy
               IF(out)then
                  WRITE(out_unit,20)iy,jy,y
                  WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
               DO ix=0,fst%mx-1
                  dx=(fst%xs(ix+1)-fst%xs(ix))/4
                  DO jx=0,4
                     x=fst%xs(ix)+dx*jx
                     CALL fspline_eval(fst,x,y,0)
                     IF(out)WRITE(out_unit,format2)ix,x,fst%f
                     IF(bin)WRITE(bin_unit)REAL(x,4),REAL(fst%f,4)
                  ENDDO
                  IF(out)WRITE(out_unit,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over y.
c-----------------------------------------------------------------------
               IF(out)WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
               IF(bin)WRITE(bin_unit)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close binary output file.
c-----------------------------------------------------------------------
      IF(bin)CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_write_xy
c-----------------------------------------------------------------------
c     subprogram 8. fspline_write_yx.
c     produces ascii and binary output for fspline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_write_yx(fst,out,bin,interpolate,filename)

      TYPE(fspline_type), INTENT(INOUT) :: fst
      LOGICAL, INTENT(IN) :: out,bin,interpolate
      CHARACTER(*) :: filename

      INTEGER :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      CHARACTER(80) :: format2,format1
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/1x,'ix = ',i3,', x = ',1p,e11.3)
 20   FORMAT(/1x,'ix = ',i3,', jx = ',i1,', x = ',1p,e11.3)
 30   FORMAT('(/4x,"iy",6x,"y",4x,',i1,'(6x,"f",i1,3x)/)')
 40   FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     open binary output file.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
      IF(bin)CALL bin_open(bin_unit,TRIM(filename),
     $     "UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(out_unit,'(1x,a)')"input data"
         WRITE(format1,30)fst%nqty
         WRITE(format2,40)fst%nqty
      ENDIF
      DO ix=0,fst%mx
         x=fst%xs(ix)
         IF(out)then
            WRITE(out_unit,10)ix,fst%xs(ix)
            WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
         ENDIF
         DO iy=0,fst%my
            y=fst%ys(iy)
            fst%f=fst%fs(ix,iy,:)
            IF(out)WRITE(out_unit,format2)iy,y,fst%f
            IF(bin)WRITE(bin_unit)REAL(y,4),REAL(fst%f,4)
         ENDDO
         IF(out)WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
         IF(bin)WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
      IF(interpolate)THEN
         IF(out)WRITE(out_unit,'(1x,a)')"interpolated data"
         DO ix=0,fst%mx-1
            dx=(fst%xs(ix+1)-fst%xs(ix))/4
            DO jx=0,4
               x=fst%xs(ix)+dx*jx
               IF(out)then
                  WRITE(out_unit,20)ix,jx,x
                  WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
               DO iy=0,fst%my-1
                  dy=(fst%ys(iy+1)-fst%ys(iy))/4
                  DO jy=0,4
                     y=fst%ys(iy)+dy*jy
                     CALL fspline_eval(fst,x,y,0)
                     IF(out)WRITE(out_unit,format2)iy,y,fst%f
                     IF(bin)WRITE(bin_unit)REAL(y,4),REAL(fst%f,4)
                  ENDDO
                  IF(out)WRITE(out_unit,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over x.
c-----------------------------------------------------------------------
               IF(out)WRITE(out_unit,format1)(iqty,iqty=1,fst%nqty)
               IF(bin)WRITE(bin_unit)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close binary output file.
c-----------------------------------------------------------------------
      IF(bin)CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_write_yx
c-----------------------------------------------------------------------
c     subprogram 9. fspline_copy.
c     copies one fspline type to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fspline_copy(fst1,fst2)

      TYPE(fspline_type), INTENT(IN) :: fst1
      TYPE(fspline_type), INTENT(INOUT) :: fst2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(fst2%xs))CALL fspline_dealloc(fst2)
      CALL fspline_alloc(fst2,fst1%mx,fst1%my,fst1%mband,fst1%nqty)
      CALL cspline_copy(fst1%cs,fst2%cs)
      fst2%xs=fst1%xs
      fst2%ys=fst1%ys
      fst2%fs=fst1%fs
      fst2%name=fst1%name
      fst2%title=fst1%title
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fspline_copy
      END MODULE fspline_mod
