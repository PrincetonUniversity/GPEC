c-----------------------------------------------------------------------
c     file zroot.f.
c     uses secant method to find root of scalar complex function f(z).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. zroot_mod.
c     1. zroot_newton.
c     2. zroot_mueller.
c-----------------------------------------------------------------------
c     subprogram 0. zroot_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE zroot_mod
      USE local_mod
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: itmax=50
      REAL(r8), PARAMETER, PRIVATE :: tol=1e-4,dzfac=.01

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. zroot_newton.
c     uses Newton's method to find root of scalar complex f(z).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE zroot_newton(ff,z,err,it)

      COMPLEX(r8) :: ff
      COMPLEX(r8), INTENT(INOUT) :: z
      REAL(r8), INTENT(OUT) :: err
      INTEGER, INTENT(OUT) :: it

      COMPLEX(r8) :: z_old,f_old,f,dz
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(a,i2,a,1p,e9.3,2(a,e10.3),a)
c-----------------------------------------------------------------------
c     find initial guess.
c-----------------------------------------------------------------------
      f=ff(z)
      dz=z*dzfac
      it=0
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO
         it=it+1
         err=ABS(dz/z)
         WRITE(*,10)"it = ",it,", err = ",err,
     $        ", z = (",REAL(z),",",AIMAG(z),")"
         IF(err < tol)EXIT
         IF(it > itmax)
     $        CALL program_stop("zroot_newton can't find root.")
         z_old=z
         z=z+dz
         f_old=f
         f=ff(z)
         dz=-f*(z-z_old)/(f-f_old)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE zroot_newton
c-----------------------------------------------------------------------
c     subprogram 2. zroot_mueller.
c     uses Mueller's method to find root of scalar complex f(z).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE zroot_mueller(ff,zz,err,it)

      COMPLEX(r8) :: ff
      COMPLEX(r8), INTENT(INOUT) :: zz
      REAL(r8), INTENT(OUT) :: err
      INTEGER, INTENT(OUT) :: it

      INTEGER :: i
      COMPLEX(r8) :: q,a,b,c,d,d1,s
      COMPLEX(r8), DIMENSION(0:2) :: f,z
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(a,i2,a,1p,e9.3,2(a,e10.3),a)
c-----------------------------------------------------------------------
c     find initial guess.
c-----------------------------------------------------------------------
      z=zz*(1+(/0,-1,1/)*dzfac)
      DO i=0,2
         f(i)=ff(z(i))
      ENDDO
      it=0
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO
         it=it+1
         err=ABS((z(0)-z(1))/z(0))
         WRITE(*,10)"it = ",it,", err = ",err,
     $        ", z = (",REAL(zz),",",AIMAG(zz),")"
         IF(err < tol)EXIT
         IF(it > itmax)
     $        CALL program_stop("zroot_mueller can't find root.")
         q=(z(0)-z(1))/(z(1)-z(2))
         a=q*(f(0)-(1+q)*f(1)+q*f(2))
         b=(2*q+1)*f(0)-(1+q)**2*f(1)+q**2*f(2)
         c=(1+q)*f(0)
         s=SQRT(b*b-4*a*c)
         d=b+s
         d1=b-s
         IF(ABS(d) < ABS(d1))d=d1
         zz=z(0)-(z(0)-z(1))*2*c/d
         z(1:2)=z(0:1)
         f(1:2)=f(0:1)
         z(0)=zz
         f(0)=ff(zz)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE zroot_mueller
      END MODULE zroot_mod
