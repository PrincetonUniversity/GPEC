c-----------------------------------------------------------------------
c     file fft.f.
c     fast fourier transform routines.
c     Numerical Recipes in Fortran, Second Edition, p. 501, four1.
c     Fortran 90 complex arithmetic version by Alan H. Glasser.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fft_mod.
c     1. fft_run.
c     2. fft_write.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. fft_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE fft_mod
      USE local_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fft_run.
c     performs fast fourier transform.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION fft_run(f,sign) RESULT(g)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: f
      COMPLEX(r8), DIMENSION(SIZE(f,1),SIZE(f,2)) :: g
      INTEGER, INTENT(IN) :: sign

      INTEGER :: i,istep,j,m,mmax,n
      COMPLEX(r8) :: w,wp
      COMPLEX(r8), DIMENSION(SIZE(f,2)) :: temp
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      n=SIZE(f,1)
      SELECT CASE(sign)
      CASE(-1)
         g=f/n
      CASE(1)
         g=f
      CASE DEFAULT
         WRITE(*,'(a,i2,a)')"sign = ",sign," is illegal, must be +/- 1"
         STOP
      END SELECT
c-----------------------------------------------------------------------
c     bit reversal.
c-----------------------------------------------------------------------
      j=1
      DO i=1,n
         IF(j > i)THEN
            temp(:)=g(j,:)
            g(j,:)=g(i,:)
            g(i,:)=temp
         ENDIF
         m=n/2
         DO
            IF(m < 1 .OR. j <= m)EXIT
            j=j-m
            m=m/2
         ENDDO
         j=j+m
      ENDDO
c-----------------------------------------------------------------------
c     Danielson-Lanczos loop.
c-----------------------------------------------------------------------
      mmax=1
      DO
         IF(n <= mmax)EXIT
         istep=2*mmax
         wp=EXP(pi*ifac/(sign*mmax))
         w=1
         DO m=1,mmax
            DO i=m,n,istep
               j=i+mmax
               temp=w*g(j,:)
               g(j,:)=g(i,:)-temp
               g(i,:)=g(i,:)+temp
            ENDDO
            w=w*wp
         ENDDO
         mmax=istep
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION fft_run
c-----------------------------------------------------------------------
c     subprogram 2. fft_write.
c     ascii and binary output of fast fourier transform.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fft_write(t,f,g,h,out,bin,out_unit,bin_unit)

      REAL(r8), DIMENSION(:), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: f,g,h
      LOGICAL :: out,bin
      INTEGER :: out_unit,bin_unit

      INTEGER :: i,k,n
      REAL(r8) :: dt,tperiod,dfreq,freqmax

      INTEGER, DIMENSION(SIZE(t)) :: m
      REAL(r8), DIMENSION(SIZE(t)) :: freq
      COMPLEX(r8), DIMENSION(SIZE(g,1),SIZE(g,2)) :: gg
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/6x,"n",7x,"dt",6x,"tperiod",5x,"dfreq",5x,"freqmax"
     $     //i8,1p,4e11.3/)
 20   FORMAT(/7x,"i",6x,"m",7x,"t",9x,"freq",7x,"re f",7x,"im f",7x,
     $     "re g",7x,"im g",7x,"re h",7x,"im h"/)
 30   FORMAT(2i8,1p,8e11.3)
c-----------------------------------------------------------------------
c     compute scalars.
c-----------------------------------------------------------------------
      n=SIZE(t)
      dt=t(2)-t(1)
      tperiod=n*dt
      dfreq=1/(n*dt)
      freqmax=dfreq*n/2
c-----------------------------------------------------------------------
c     compute auxiliary arrays.
c-----------------------------------------------------------------------
      m(1:n/2-1)=(/(i,i=1-n/2,-1)/)
      m(n/2:n)=(/(i,i=0,n/2)/)
      freq=m*dfreq
      gg(1:n/2-1,:)=g(n/2+2:n,:)
      gg(n/2:n,:)=g(1:n/2+1,:)
c-----------------------------------------------------------------------
c     diagnose fft.
c-----------------------------------------------------------------------
      WRITE(out_unit,10)n,dt,tperiod,dfreq,freqmax
      DO k=1,SIZE(f,2)
         IF(out)THEN
            WRITE(out_unit,'(/a,i3)')"k = ",k
            WRITE(out_unit,20)
         ENDIF
         DO i=1,n
            IF(out)WRITE(out_unit,30)i-1,m(i),t(i),freq(i),
     $           f(i,k),gg(i,k),h(i,k)
            IF(bin)WRITE(bin_unit)
     $           REAL(i-1,4),REAL(m(i),4),REAL(t(i),4),REAL(freq(i),4),
     $           REAL(REAL(f(i,k)),4),REAL(AIMAG(f(i,k)),4),
     $           REAL(REAL(gg(i,k)),4),REAL(AIMAG(gg(i,k)),4),
     $           REAL(REAL(h(i,k)),4),REAL(AIMAG(h(i,k)),4)
         ENDDO
         IF(out)WRITE(out_unit,20)
         IF(bin)WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fft_write
      END MODULE fft_mod
