c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     ISMATH: basic math routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ismath_mod.
c     1. iscdftf.
c     2. iscdftb.
c-----------------------------------------------------------------------
c     subprogram 0. ismath_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ismath_mod
      USE local_mod
      USE bicube_mod
      USE fspline_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iscdftf.
c     compute 1D discrete Fourier forward transform by simple way.
c-----------------------------------------------------------------------
      SUBROUTINE iscdftf(m,ms,func,fs,funcm)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,ms,fs
      INTEGER, DIMENSION(ms), INTENT(IN) :: m
      COMPLEX(r8), DIMENSION(0:fs), INTENT(IN) :: func
      COMPLEX(r8), DIMENSION(ms), INTENT(OUT) :: funcm

      funcm=0

      DO i=1,ms
         DO j=0,fs-1
            funcm(i) = funcm(i)+(1.0/fs)*func(j)*
     $           EXP(-(twopi*ifac*m(i)*j)/fs)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iscdftf
c-----------------------------------------------------------------------
c     subprogram 2. iscdfb.
c     compute 1D discrete Fourier backward transform by simple way.
c-----------------------------------------------------------------------
      SUBROUTINE iscdftb(m,ms,func,fs,funcm)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,ms,fs
      INTEGER, DIMENSION(ms), INTENT(IN) :: m
      COMPLEX(r8), DIMENSION(ms), INTENT(IN) :: funcm
      COMPLEX(r8), DIMENSION(0:fs), INTENT(OUT) :: func

      func=0

      DO i=0,fs-1
         DO j=1,ms
            func(i) = func(i)+funcm(j)*
     $           EXP((twopi*ifac*m(j)*i)/fs)
         ENDDO
      ENDDO
      func(fs)=func(0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iscdftb
c-----------------------------------------------------------------------
c     subprogram 3. isintg.
c     compute root of discrete function by the secant method.
c-----------------------------------------------------------------------
      FUNCTION issect(gridnum,xvec,yvec,yval)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: gridnum
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: xvec
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: yvec
      REAL(r8), INTENT(IN) :: yval

      REAL(r8) :: xval1,xval2,xval3,yval1,yval2,yval3,issect
      TYPE(spline_type) :: yfun

      CALL spline_alloc(yfun,gridnum,1)
      yfun%xs=xvec
      yfun%fs(:,1)=yvec
      CALL spline_fit(yfun,"extrap")

      xval1=xvec(0)
      yval1=yvec(0)
      xval2=xvec(gridnum)
      yval2=yvec(gridnum)

      xval3=yval
      CALL spline_eval(yfun,xval3,0)
      yval3=yfun%f(1)

      DO
         IF (ABS(yval3-yval) .LT. 1e-4) EXIT
         IF (yval3 .GT. yval) THEN
            xval2=xval3
            yval2=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ELSE
            xval1=xval3
            yval1=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ENDIF
         CALL spline_eval(yfun,xval3,0)
         yval3=yfun%f(1)
      ENDDO
      issect=xval3
      CALL spline_dealloc(yfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION issect

      END MODULE ismath_mod
