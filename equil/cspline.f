c-----------------------------------------------------------------------
c     file cspline.f
c     fits complex functions to cubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. cspline_mod.
c     1. cspline_alloc.
c     2. cspline_dealloc.
c     3. cspline_fit.
c     4. cspline_fac.
c     5. cspline_eval.
c     6. cspline_all_eval.
c     7. cspline_write.
c     8. cspline_write_log.
c     9. cspline_int.
c     10. cspline_triluf.
c     11. cspline_trilus.
c     12. cspline_sherman.
c     13. cspline_morrison.
c     14. cspline_copy.
c-----------------------------------------------------------------------
c     subprogram 0. cspline_type definition.
c     defines cspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE cspline_mod
      USE local_mod
      IMPLICIT NONE
      
      TYPE :: cspline_type
      INTEGER :: mx,nqty,ix
      REAL(r8), DIMENSION(:), POINTER :: xs
      REAL(r8), DIMENSION(2) :: x0
      REAL(r8), DIMENSION(:,:), POINTER :: xpower
      COMPLEX(r8), DIMENSION(:), POINTER :: f,f1,f2,f3
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fs,fs1,fsi
      CHARACTER(6), DIMENSION(:), POINTER :: title
      CHARACTER(6) :: name
      LOGICAL :: periodic
      END TYPE cspline_type
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. cspline_alloc.
c     allocates space for cspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_alloc(spl,mx,nqty)
      
      INTEGER, INTENT(IN) :: mx,nqty
      TYPE(cspline_type), INTENT(INOUT) :: spl
c-----------------------------------------------------------------------
c     set scalars.
c-----------------------------------------------------------------------
      spl%mx=mx
      spl%nqty=nqty
      spl%ix=0
      spl%periodic=.FALSE.
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(spl%xs(0:mx))
      ALLOCATE(spl%f(nqty))
      ALLOCATE(spl%f1(nqty))
      ALLOCATE(spl%f2(nqty))
      ALLOCATE(spl%f3(nqty))
      ALLOCATE(spl%title(0:nqty))
      ALLOCATE(spl%fs(0:mx,nqty))
      ALLOCATE(spl%fs1(0:mx,nqty))
      ALLOCATE(spl%xpower(2,nqty))
      spl%xpower=0
      spl%x0=0
      NULLIFY(spl%fsi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_alloc
c-----------------------------------------------------------------------
c     subprogram 2. cspline_dealloc.
c     deallocates space for cspline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_dealloc(spl)
      
      TYPE(cspline_type), INTENT(INOUT) :: spl
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(spl%xs)
      DEALLOCATE(spl%f)
      DEALLOCATE(spl%f1)
      DEALLOCATE(spl%f2)
      DEALLOCATE(spl%f3)
      DEALLOCATE(spl%title)
      DEALLOCATE(spl%fs)
      DEALLOCATE(spl%fs1)
      DEALLOCATE(spl%xpower)
      IF(ASSOCIATED(spl%fsi))DEALLOCATE(spl%fsi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. cspline_fit.
c     fits complex functions to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_fit(spl,endmode)
      
      TYPE(cspline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode
      
      INTEGER :: iqty,iside
      REAL(r8), DIMENSION(-1:1,0:spl%mx) :: a
      REAL(r8), DIMENSION(spl%mx) :: b
      REAL(r8), DIMENSION(4) :: cl,cr
      REAL(r8), DIMENSION(0:spl%mx) :: xfac
c-----------------------------------------------------------------------
c     extract powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(spl%xs-spl%x0(iside))**spl%xpower(iside,iqty)
               spl%fs(:,iqty)=spl%fs(:,iqty)*xfac
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set up grid matrix.
c-----------------------------------------------------------------------
      CALL cspline_fac(spl,a,b,cl,cr,endmode)
c-----------------------------------------------------------------------
c     compute first derivatives, interior.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         spl%fs1(1:spl%mx-1,iqty)=
     $        3*((spl%fs(2:spl%mx,iqty)-spl%fs(1:spl%mx-1,iqty))
     $        *b(2:spl%mx)
     $        +(spl%fs(1:spl%mx-1,iqty)-spl%fs(0:spl%mx-2,iqty))
     $        *b(1:spl%mx-1))
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         DO iqty=1,spl%nqty
            spl%fs1(0,iqty)=SUM(cl(1:4)*spl%fs(0:3,iqty))
            spl%fs1(spl%mx,iqty)=SUM(cr(1:4)
     $           *spl%fs(spl%mx:spl%mx-3:-1,iqty))
            spl%fs1(1,iqty)=spl%fs1(1,iqty)-spl%fs1(0,iqty)
     $           /(spl%xs(1)-spl%xs(0))
            spl%fs1(spl%mx-1,iqty)=
     $           spl%fs1(spl%mx-1,iqty)-spl%fs1(spl%mx,iqty)
     $           /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
         ENDDO
         CALL cspline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         spl%fs1(1,:)=spl%fs1(1,:)-(2*spl%fs(1,:)
     $        -spl%fs(0,:)-spl%fs(2,:))*2*b(1)
         spl%fs1(spl%mx-1,:)=spl%fs1(spl%mx-1,:)
     $        +(2*spl%fs(spl%mx-1,:)-spl%fs(spl%mx,:)
     $        -spl%fs(spl%mx-2,:))*2*b(spl%mx)
         CALL cspline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
         spl%fs1(0,:)=(2*(2*spl%fs(1,:)-spl%fs(0,:)-spl%fs(2,:))
     $        +(spl%fs1(1,:)+spl%fs1(2,:))*(spl%xs(2)-spl%xs(1))
     $        -spl%fs1(1,:)*(spl%xs(1)-spl%xs(0)))/(spl%xs(1)-spl%xs(0))
         spl%fs1(spl%mx,:)=
     $        (2*(spl%fs(spl%mx-2,:)+spl%fs(spl%mx,:)
     $        -2*spl%fs(spl%mx-1,:))
     $        +(spl%fs1(spl%mx-1,:)+spl%fs1(spl%mx-2,:))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        -spl%fs1(spl%mx-1,:)
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-1)))
     $        /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boudary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         spl%periodic=.TRUE.
         spl%fs1(0,:)=3*((spl%fs(1,:)-spl%fs(0,:))*b(1)
     $        +(spl%fs(0,:)-spl%fs(spl%mx-1,:))*b(spl%mx))
         CALL cspline_morrison(a(:,0:spl%mx-1),spl%fs1(0:spl%mx-1,:))
         spl%fs1(spl%mx,:)=spl%fs1(0,:)
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_fit
c-----------------------------------------------------------------------
c     subprogram 4. cspline_fac.
c     sets up matrix for cubic spline fitting.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_fac(spl,a,b,cl,cr,endmode)
      
      TYPE(cspline_type), INTENT(IN) :: spl
      REAL(r8), DIMENSION(-1:1,0:spl%mx), INTENT(OUT) :: a
      REAL(r8), DIMENSION(spl%mx), INTENT(OUT) :: b
      REAL(r8), DIMENSION(4), INTENT(OUT) :: cl,cr
      CHARACTER(*), INTENT(IN) :: endmode
      
      INTEGER :: j
c-----------------------------------------------------------------------
c     compute interior matrix.
c-----------------------------------------------------------------------
      b=1/(spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1))
      DO j=1,spl%mx-1
         a(-1,j)=b(j)
         a(0,j)=2*(b(j)+b(j+1))
         a(1,j)=b(j+1)
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         b=b*b
         cl(1)=(spl%xs(0)*(3*spl%xs(0)
     $        -2*(spl%xs(1)+spl%xs(2)+spl%xs(3)))
     $        +spl%xs(1)*spl%xs(2)+spl%xs(1)*spl%xs(3)
     $        +spl%xs(2)*spl%xs(3))
     $        /((spl%xs(0)-spl%xs(1))*(spl%xs(0)-spl%xs(2))
     $        *(spl%xs(0)-spl%xs(3)))
         cl(2)=((spl%xs(2)-spl%xs(0))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(1)-spl%xs(0))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(1)-spl%xs(3)))
         cl(3)=((spl%xs(0)-spl%xs(1))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(0)-spl%xs(2))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(3)-spl%xs(2)))
         cl(4)=((spl%xs(1)-spl%xs(0))*(spl%xs(2)-spl%xs(0)))
     $        /((spl%xs(3)-spl%xs(0))*(spl%xs(3)-spl%xs(1))
     $        *(spl%xs(3)-spl%xs(2)))
         cr(1)=(spl%xs(spl%mx)*(3*spl%xs(spl%mx)
     $        -2*(spl%xs(spl%mx-1)+spl%xs(spl%mx-2)
     $        +spl%xs(spl%mx-3)))+spl%xs(spl%mx-1)
     $        *spl%xs(spl%mx-2)+spl%xs(spl%mx-1)*spl%xs(spl%mx
     $        -3)+spl%xs(spl%mx-2)*spl%xs(spl%mx-3))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-2))*(spl%xs(spl%mx
     $        )-spl%xs(spl%mx-3)))
         cr(2)=((spl%xs(spl%mx-2)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-3)))
         cr(3)=((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         cr(4)=((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-2)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-3)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         CALL cspline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         b=b*b
         a(0,1)=a(0,1)+(spl%xs(2)+spl%xs(0)-2*spl%xs(1))*b(1)
         a(1,1)=a(1,1)+(spl%xs(2)-spl%xs(1))*b(1)
         a(0,spl%mx-1)=a(0,spl%mx-1)
     $        +(2*spl%xs(spl%mx-1)-spl%xs(spl%mx-2)
     $        -spl%xs(spl%mx))*b(spl%mx)
         a(-1,spl%mx-1)=a(-1,spl%mx-1)
     $        +(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))*b(spl%mx)
         CALL cspline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boundary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         a(0,0:spl%mx:spl%mx)=2*(b(spl%mx)+b(1))
         a(1,0)=b(1)
         a(-1,0)=b(spl%mx)
         b=b*b
         CALL cspline_sherman(a(:,0:spl%mx-1))
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_fac
c-----------------------------------------------------------------------
c     subprogram 5. cspline_eval.
c     evaluates complex cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_eval(spl,x,mode)
      
      TYPE(cspline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: mode
      
      INTEGER :: iqty,iside
      REAL(r8) :: xx,d,z,z1,xfac,dx
      COMPLEX(r8) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      xx=x
      spl%ix=MAX(spl%ix,0)
      spl%ix=MIN(spl%ix,spl%mx-1)
c-----------------------------------------------------------------------
c     normalize interval for periodic splines.
c-----------------------------------------------------------------------
      IF(spl%periodic)THEN
         DO
            IF(xx < spl%xs(spl%mx))EXIT
            xx=xx-spl%xs(spl%mx)
         ENDDO
         DO
            IF(xx >= spl%xs(0))EXIT
            xx=xx+spl%xs(spl%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find cubic spline interval.
c-----------------------------------------------------------------------
      DO
         IF(xx >= spl%xs(spl%ix).OR.spl%ix <= 0)EXIT
         spl%ix=spl%ix-1
      ENDDO
      DO
         IF(xx < spl%xs(spl%ix+1).OR.spl%ix >= spl%mx-1)EXIT
         spl%ix=spl%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      d=spl%xs(spl%ix+1)-spl%xs(spl%ix)
      z=(xx-spl%xs(spl%ix))/d
      z1=1-z
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      spl%f=spl%fs(spl%ix,:)*z1*z1*(3-2*z1)
     $     +spl%fs(spl%ix+1,:)*z*z*(3-2*z)
     $     +d*z*z1*(spl%fs1(spl%ix,:)*z1
     $     -spl%fs1(spl%ix+1,:)*z)
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         spl%f1=6*(spl%fs(spl%ix+1,:)
     $        -spl%fs(spl%ix,:))*z*z1/d
     $        +spl%fs1(spl%ix,:)*z1*(3*z1-2)
     $        +spl%fs1(spl%ix+1,:)*z*(3*z-2)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         spl%f2=(6*(spl%fs(spl%ix+1,:)
     $        -spl%fs(spl%ix,:))*(z1-z)/d
     $        -spl%fs1(spl%ix,:)*(6*z1-2)
     $        +spl%fs1(spl%ix+1,:)*(6*z-2))/d
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         spl%f3=(12*(spl%fs(spl%ix,:)
     $        -spl%fs(spl%ix+1,:))/d
     $        +6*(spl%fs1(spl%ix,:)
     $        +spl%fs1(spl%ix+1,:)))/(d*d)
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=spl%f(iqty)*xfac
            IF(mode > 0)g1=(spl%f1(iqty)+spl%f(iqty)
     $           *spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(spl%f2(iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*spl%f1(iqty)+(spl%xpower(iside,iqty)-1)
     $           *spl%f(iqty)/dx))*xfac
            IF(mode > 2)g3=(spl%f3(iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*spl%f2(iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*spl%f1(iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *spl%f(iqty))))*xfac
            spl%f(iqty)=g
            IF(mode > 0)spl%f1(iqty)=g1
            IF(mode > 1)spl%f2(iqty)=g2
            IF(mode > 2)spl%f3(iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_eval
c-----------------------------------------------------------------------
c     subprogram 6. cspline_all_eval.
c     evaluates cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_all_eval(spl,z,f,f1,f2,f3,mode)

      TYPE(cspline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: z
      COMPLEX(r8), DIMENSION(spl%mx,spl%nqty), INTENT(OUT) ::
     $     f,f1,f2,f3
      INTEGER, INTENT(IN) :: mode

      INTEGER :: iqty,nqty,n,iside
      REAL(r8) :: z1
      REAL(r8), DIMENSION(spl%mx) :: d,xfac,dx
      COMPLEX(r8), DIMENSION(spl%mx) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      n=spl%mx
      nqty=spl%nqty
      z1=1-z
      d=spl%xs(1:n)-spl%xs(0:n-1)
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         f(:,iqty)=spl%fs(0:n-1,iqty)*z1*z1*(3-2*z1)
     $        +spl%fs(1:n,iqty)*z*z*(3-2*z)
     $        +d*z*z1*(spl%fs1(0:n-1,iqty)*z1-spl%fs1(1:n,iqty)*z)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         DO iqty=1,nqty
            f1(:,iqty)=6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*z*z1/d
     $           +spl%fs1(0:n-1,iqty)*z1*(3*z1-2)
     $           +spl%fs1(1:n,iqty)*z*(3*z-2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         DO iqty=1,nqty
            f2(:,iqty)=(6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*(z1-z)/d
     $           -spl%fs1(0:n-1,iqty)*(6*z1-2)
     $           +spl%fs1(1:n,iqty)*(6*z-2))/d
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         DO iqty=1,nqty
            f3(:,iqty)=(12*(spl%fs(0:n-1,iqty)-spl%fs(1:n,iqty))/d
     $           +6*(spl%fs1(0:n-1,iqty)+spl%fs1(1:n,iqty)))/(d*d)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=(spl%xs(0:spl%mx-1)+z*d(1:spl%mx))-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=f(:,iqty)*xfac
            IF(mode > 0)g1=(f1(:,iqty)
     $           +f(:,iqty)*spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(f2(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*f1(:,iqty)+(spl%xpower(iside,iqty)-1)
     $           *f(:,iqty)/dx))*xfac
     $           
            IF(mode > 2)g3=(f3(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*f2(:,iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*f1(:,iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *f(:,iqty))))*xfac
            f(:,iqty)=g
            IF(mode > 0)f1(:,iqty)=g1
            IF(mode > 1)f2(:,iqty)=g2
            IF(mode > 2)f3(:,iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_all_eval
c-----------------------------------------------------------------------
c     subprogram 7. cspline_write.
c     produces ascii and binary output.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_write(spl,out,bin,iua,iub,interp)

      TYPE(cspline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      CHARACTER(80) :: format1,format2
      INTEGER :: i,j
      REAL(r8) :: x,dx
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"i",4x,a6,1x,',i2.2,'(2x,"re ",a6,2x,"im ",a6)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
 30   FORMAT('(/4x,"i",2x,"j",',i2.2,'(4x,a6,1x)/)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN
c-----------------------------------------------------------------------
c     print ascii tables of node values and derivatives.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(format1,10)spl%nqty
         WRITE(format2,20)2*spl%nqty+1
         WRITE(iua,'(/1x,a)')'node values:'
         WRITE(iua,format1)spl%title(0),
     $        (spl%title(i),spl%title(i),i=1,spl%nqty)
      ENDIF
      DO i=0,spl%mx
         CALL cspline_eval(spl,spl%xs(i),0)
         IF(out)WRITE(iua,format2)spl%xs(i),spl%f
         IF(bin)WRITE(iub)REAL(spl%xs(i),4),
     $        REAL(spl%f,4),REAL(AIMAG(spl%f),4)
      ENDDO
      IF(out)WRITE(iua,format1)spl%title(0),
     $     (spl%title(i),spl%title(i),i=1,spl%nqty)
      IF(bin)WRITE(iub)
      IF(.NOT. interp)RETURN
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,'(/1x,a)')'interpolated values:'
         WRITE(iua,format1)spl%title(0),
     $        (spl%title(i),spl%title(i),i=1,spl%nqty)
      ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         dx=(spl%xs(i+1)-spl%xs(i))/4
         DO j=0,4
            x=spl%xs(i)+j*dx
            CALL cspline_eval(spl,x,0)
            IF(out)WRITE(iua,format2)i,x,spl%f
            IF(bin)WRITE(iub)REAL(x,4),
     $           REAL(spl%f,4),REAL(AIMAG(spl%f),4)
         ENDDO
         IF(out)WRITE(iua,'(1x)')
      ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
      x=spl%xs(spl%mx)
      CALL cspline_eval(spl,x,0)
      IF(out)THEN
         WRITE(iua,format2)i,x,spl%f
         WRITE(iua,format1)spl%title(0),
     $        (spl%title(i),spl%title(i),i=1,spl%nqty)
      ENDIF
      IF(bin)THEN
         WRITE(iub)REAL(x,4),
     $        REAL(spl%f,4),REAL(AIMAG(spl%f),4)
         WRITE(iub)
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_write
c-----------------------------------------------------------------------
c     subprogram 8. cspline_write_log
c     produces ascii and binary output of logs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_write_log(spl,out,bin,iua,iub,interp,stride,
     $     xend)

      TYPE(cspline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      INTEGER, INTENT(IN) :: stride
      LOGICAL, INTENT(IN) :: interp
      REAL(r8), DIMENSION(2) :: xend

      CHARACTER(50) :: format1,format2
      INTEGER :: iqty,ix,j,offset
      REAL(r8), PARAMETER :: epsilon=-20*alog10
      REAL(r8) :: x,dx
      REAL(r8), DIMENSION(2) :: xlog
      COMPLEX(r8), DIMENSION(spl%nqty/stride) :: flog
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"i",4x,a6,1x,',i2.2,'(2x,"re ",a6,2x,"im ",a6)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. out .AND. .NOT. bin)RETURN
      IF(MOD(spl%nqty,stride) /= 0)THEN
         WRITE(*,'(2(a,i3))')"Cspline_write_log: nqty = ",spl%nqty,
     $        " is not an integral multiple of stride = ",stride
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     write title and start loop over offsets.
c-----------------------------------------------------------------------
      IF(OUT)WRITE(iua,'(1x,a)')
     $     "Output from cspline_write_log for "//TRIM(spl%name)//":"
      DO offset=1,stride
c-----------------------------------------------------------------------
c     print ascii table of node values.
c-----------------------------------------------------------------------
         IF(out)THEN
            WRITE(format1,10)spl%nqty/stride
            WRITE(format2,20)2*spl%nqty/stride+1
            WRITE(iua,'(/1x,a,i3,a)')
     $           "input values for offset = ",offset-1,":"
            WRITE(iua,format1)spl%title(0),(spl%title(iqty),
     $           spl%title(iqty),iqty=offset,spl%nqty,stride)
            DO ix=0,spl%mx
               CALL cspline_eval(spl,spl%xs(ix),0)
               WRITE(iua,format2)ix,spl%xs(ix),
     $              spl%f(offset:spl%nqty:stride)
            ENDDO
            WRITE(iua,format1)spl%title(0),(spl%title(iqty),
     $           spl%title(iqty),iqty=offset,spl%nqty,stride)
         ENDIF
c-----------------------------------------------------------------------
c     compute logs.
c-----------------------------------------------------------------------
         IF(bin)THEN
            DO ix=0,spl%mx
               xlog=LOG10(ABS(spl%xs(ix)-xend))
               CALL cspline_eval(spl,spl%xs(ix),0)
               WHERE(spl%f(offset:spl%nqty:stride) /= 0)
                  flog=LOG(spl%f(offset:spl%nqty:stride))
               ELSEWHERE
                  flog=epsilon
               ENDWHERE
c-----------------------------------------------------------------------
c     print binary table of node values.
c-----------------------------------------------------------------------
               WRITE(iub)REAL(spl%xs(ix),4),REAL(xlog,4),
     $              (REAL(REAL(flog(iqty))/alog10,4),
     $              REAL(AIMAG(flog(iqty))*rtod,4),
     $              iqty=1,SIZE(flog))
            ENDDO
            WRITE(iub)
         ENDIF
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
         IF(interp)THEN
            IF(out)THEN
               WRITE(iua,'(/1x,a,i3)')
     $              "interpolated values for offset = ",offset,":"
               WRITE(iua,format1)spl%title(0),
     $              (spl%title(iqty),spl%title(iqty),
     $              iqty=offset,spl%nqty,stride)
            ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
            DO ix=0,spl%mx-1
               dx=(spl%xs(ix+1)-spl%xs(ix))/4
               DO j=0,3
                  x=spl%xs(ix)+j*dx
                  xlog=LOG10(ABS(x-xend))
                  CALL cspline_eval(spl,x,0)
                  IF(out)WRITE(iua,format2)ix,x,(spl%f(iqty),
     $                 iqty=offset,spl%nqty,stride)
                  IF(bin)THEN
                     WHERE(spl%f(offset:spl%nqty:stride) /= 0)
                        flog=LOG(spl%f(offset:spl%nqty:stride))
                     ELSEWHERE
                        flog=epsilon
                     ENDWHERE
                     WRITE(iub)REAL(x,4),REAL(xlog,4),
     $                    (REAL(DREAL(flog(iqty))/alog10,4),
     $                    REAL(AIMAG(flog(iqty))*rtod,4),
     $                    iqty=1,SIZE(flog))
                  ENDIF
               ENDDO
               IF(out)WRITE(iua,'(1x)')
            ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
            x=spl%xs(spl%mx)
            xlog=LOG10(ABS(x-xend))
            CALL cspline_eval(spl,x,0)
            IF(out)THEN
               WRITE(iua,format2)ix,x,(spl%f(iqty),
     $              iqty=offset,spl%nqty,stride)
               WRITE(iua,format1)spl%title(0),(spl%title(iqty),
     $              spl%title(iqty),iqty=offset,spl%nqty,stride)
            ENDIF
            IF(bin)THEN
               WHERE(spl%f(offset:spl%nqty:stride) /= 0)
                  flog=LOG(spl%f(offset:spl%nqty:stride))
               ELSEWHERE
                  flog=epsilon
               ENDWHERE
               WRITE(iub)REAL(x,4),REAL(xlog,4),
     $              (REAL(DREAL(flog(iqty))/alog10,4),
     $              REAL(AIMAG(flog(iqty))*rtod,4),
     $              iqty=1,SIZE(flog))
               WRITE(iub)
            ENDIF
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_write_log
c-----------------------------------------------------------------------
c     subprogram 9. cspline_int.
c     integrates complex cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_int(spl)

      TYPE(cspline_type), INTENT(INOUT) :: spl

      INTEGER :: ix,iqty,ig
      REAL(r8), DIMENSION(spl%mx) :: dx
      COMPLEX(r8), DIMENSION(spl%mx,spl%nqty) :: term,f

      INTEGER, PARAMETER :: mg=4
      REAL(r8), DIMENSION(mg) :: xg=(1+(/-0.861136311594053_r8,
     $     -0.339981043584856_r8,0.339981043584856_r8,
     $     0.861136311594053_r8/))/2
      REAL(r8), DIMENSION(mg) :: wg=(/0.347854845137454_r8,
     $     0.652145154862546_r8,0.652145154862546_r8,
     $     0.347854845137454_r8/)/2
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(spl%fsi))ALLOCATE(spl%fsi(0:spl%mx,spl%nqty))
      dx=spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1)
      term=0
c-----------------------------------------------------------------------
c     compute integrals over intervals.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         IF(spl%xpower(1,iqty) == 0 .AND. spl%xpower(2,iqty) == 0)THEN
            term(:,iqty)=dx/12
     $           *(6*(spl%fs(0:spl%mx-1,iqty)+spl%fs(1:spl%mx,iqty))
     $           +dx*(spl%fs1(0:spl%mx-1,iqty)-spl%fs1(1:spl%mx,iqty)))
         ELSE
            DO ig=1,mg
               CALL cspline_all_eval(spl,xg(ig),f,f,f,f,0)
               term(:,iqty)=term(:,iqty)+dx*wg(ig)*f(:,iqty)
            ENDDO
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     accumulate over intervals.
c-----------------------------------------------------------------------
      spl%fsi(0,:)=0
      DO ix=1,spl%mx
         spl%fsi(ix,:)=spl%fsi(ix-1,:)+term(ix,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_int
c-----------------------------------------------------------------------
c     subprogram 10. cspline_triluf.
c     performs tridiagonal LU factorization.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_triluf(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: i,j,k,jmin,jmax,n
c-----------------------------------------------------------------------
c     begin loop over rows and define limits.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         jmin=MAX(1-i,-1)
         jmax=MIN(n-i,1)
c-----------------------------------------------------------------------
c     compute lower elements.
c-----------------------------------------------------------------------
         DO j=jmin,-1
            DO k=MAX(jmin,j-1),j-1
               a(j,i)=a(j,i)-a(k,i)*a(j-k,i+k)
            ENDDO
            a(j,i)=a(j,i)*a(0,i+j)
         ENDDO
c-----------------------------------------------------------------------
c     compute diagonal element
c-----------------------------------------------------------------------
         DO k=MAX(jmin,-1),-1
            a(0,i)=a(0,i)-a(k,i)*a(-k,i+k)
         ENDDO
         a(0,i)=1/a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_triluf
c-----------------------------------------------------------------------
c     subprogram 11. cspline_trilus.
c     performs tridiagonal LU solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_trilus(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      COMPLEX(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: i,j,n
c-----------------------------------------------------------------------
c     down sweep.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         DO j=MAX(1-i,-1),-1
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     up sweep.
c-----------------------------------------------------------------------
      DO i=n,1,-1
         DO j=1,MIN(n-i,1)
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
         x(i,:)=x(i,:)*a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_trilus
c-----------------------------------------------------------------------
c     subprogram 12. cspline_sherman.
c     uses Sherman-Morrison formula to factor periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_sherman(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: j,n
      COMPLEX(r8), DIMENSION(SIZE(a,2),1) :: u
c-----------------------------------------------------------------------
c     prepare matrices.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      a(0,1)=a(0,1)-a(-1,1)
      a(0,n)=a(0,n)-a(-1,1)
      u=RESHAPE((/one,(zero,j=2,n-1),one/),SHAPE(u))
      CALL cspline_triluf(a)
      CALL cspline_trilus(a,u)
      a(-1,1)=a(-1,1)/(1+a(-1,1)*(u(1,1)+u(n,1)))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_sherman
c-----------------------------------------------------------------------
c     subprogram 13. cspline_morrison.
c     uses Sherman-Morrison formula to solve periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_morrison(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      COMPLEX(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: n
      COMPLEX(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: y
c-----------------------------------------------------------------------
c     solve for x.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      y=x
      CALL cspline_trilus(a,y)
      x(1,:)=x(1,:)-a(-1,1)*(y(1,:)+y(n,:))
      x(n,:)=x(n,:)-a(-1,1)*(y(1,:)+y(n,:))
      CALL cspline_trilus(a,x)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_morrison
c-----------------------------------------------------------------------
c     subprogram 14. cspline_copy.
c     copies one cspline_type to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cspline_copy(spl1,spl2)

      TYPE(cspline_type), INTENT(IN) :: spl1
      TYPE(cspline_type), INTENT(INOUT) :: spl2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(spl2%xs))CALL cspline_dealloc(spl2)
      CALL cspline_alloc(spl2,spl1%mx,spl1%nqty)
      spl2%xs=spl1%xs
      spl2%fs=spl1%fs
      spl2%fs1=spl1%fs1
      spl2%name=spl1%name
      spl2%title=spl1%title
      spl2%periodic=spl1%periodic
      spl2%xpower=spl1%xpower
      spl2%x0=spl1%x0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cspline_copy
      END MODULE cspline_mod
