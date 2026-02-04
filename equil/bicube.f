c-----------------------------------------------------------------------
c     file bicube.f.
c     fits functions to bicubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. bicube_type definition.
c      1. bicube_alloc.
c      2. bicube_dealloc.
c      3. bicube_fit.
c      4. bicube_lsfit.
c      5. bicube_eval.
c      5a. bicube_eval_external.
c      6. bicube_getco.
c      6a. bicube_getco_external.
c      7. bicube_all_eval.
c      8. bicube_all_getco.
c      9. bicube_write_xy.
c     10. bicube_write_yx.
c     11. bicube_write_arrays.
c     12. bicube_copy.
c     13. bicube_extrema.
c-----------------------------------------------------------------------
c     subprogram 0. bicube_type definition.
c     defines bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE bicube_mod
      USE spline_mod
      IMPLICIT NONE

      TYPE :: bicube_type
      INTEGER :: mx,my,nqty,ix,iy
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xext,yext,fext
      REAL(r8), DIMENSION(2) :: x0,y0
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xs,ys
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xpower,ypower
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs,fsx,fsy,fsxy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy,fxx,fxy,fyy
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: cmats
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gs,gsx,gsy,gsxy,
     $     gsxx,gsyy
      CHARACTER(6) :: xtitle,ytitle
      CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
      CHARACTER(6) :: name
      LOGICAL, DIMENSION(2) :: periodic
      END TYPE bicube_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. bicube_alloc.
c     allocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_alloc(bcs,mx,my,nqty)

      INTEGER, INTENT(IN) :: mx,my,nqty
      TYPE(bicube_type), INTENT(OUT) :: bcs
c-----------------------------------------------------------------------
c     set scalars.
c-----------------------------------------------------------------------
      bcs%mx=mx
      bcs%my=my
      bcs%ix=0
      bcs%iy=0
      bcs%nqty=nqty
      bcs%periodic=.FALSE.
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(bcs%xs(0:mx))
      ALLOCATE(bcs%ys(0:my))
      ALLOCATE(bcs%fs(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsx(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsy(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsxy(0:mx,0:my,nqty))
      ALLOCATE(bcs%title(nqty))
      ALLOCATE(bcs%f(nqty))
      ALLOCATE(bcs%fx(nqty))
      ALLOCATE(bcs%fy(nqty))
      ALLOCATE(bcs%fxx(nqty))
      ALLOCATE(bcs%fxy(nqty))
      ALLOCATE(bcs%fyy(nqty))
      ALLOCATE(bcs%xpower(2,nqty),bcs%ypower(2,nqty))
      ALLOCATE(bcs%xext(2,nqty),bcs%yext(2,nqty),bcs%fext(2,nqty))
      bcs%xpower=0
      bcs%ypower=0
      bcs%x0=0
      bcs%y0=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_alloc
c-----------------------------------------------------------------------
c     subprogram 2. bicube_dealloc.
c     deallocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_dealloc(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(bcs%xs)
      DEALLOCATE(bcs%ys)
      DEALLOCATE(bcs%fs)
      DEALLOCATE(bcs%fsx)
      DEALLOCATE(bcs%fsy)
      DEALLOCATE(bcs%fsxy)
      DEALLOCATE(bcs%title)
      DEALLOCATE(bcs%f)
      DEALLOCATE(bcs%fx)
      DEALLOCATE(bcs%fy)
      DEALLOCATE(bcs%fxx)
      DEALLOCATE(bcs%fxy)
      DEALLOCATE(bcs%fyy)
      DEALLOCATE(bcs%xpower,bcs%ypower)
      DEALLOCATE(bcs%xext,bcs%yext,bcs%fext)
      IF(ALLOCATED(bcs%cmats))DEALLOCATE(bcs%cmats)
      IF(ALLOCATED(bcs%gs))DEALLOCATE(bcs%gs)
      IF(ALLOCATED(bcs%gsx))DEALLOCATE(bcs%gsx)
      IF(ALLOCATED(bcs%gsy))DEALLOCATE(bcs%gsy)
      IF(ALLOCATED(bcs%gsxx))DEALLOCATE(bcs%gsxx)
      IF(ALLOCATED(bcs%gsxy))DEALLOCATE(bcs%gsxy)
      IF(ALLOCATED(bcs%gsyy))DEALLOCATE(bcs%gsyy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. bicube_fit.
c     fits functions to bicubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_fit(bcs,endmode1,endmode2)

      TYPE(bicube_type), INTENT(INOUT), TARGET :: bcs
      CHARACTER(*), INTENT(IN) :: endmode1,endmode2

      INTEGER :: iqty,iside,ix,iy
      REAL(r8), DIMENSION(0:bcs%mx) :: xfac
      REAL(r8), DIMENSION(0:bcs%my) :: yfac
      TYPE(spline_type) :: spl

      REAL(r8), DIMENSION(:,:,:), POINTER :: fs,fsx,fsy,fsxy
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      fs => bcs%fs
      fsx => bcs%fsx
      fsy => bcs%fsy
      fsxy => bcs%fsxy
c-----------------------------------------------------------------------
c     extract x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(bcs%xs-bcs%x0(iside))**bcs%xpower(iside,iqty)
               DO iy=0,bcs%my
                  bcs%fs(:,iy,iqty)=bcs%fs(:,iy,iqty)*xfac
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     extract y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) /= 0)THEN
               yfac=1/ABS(bcs%ys-bcs%y0(iside))**bcs%ypower(iside,iqty)
               DO ix=0,bcs%mx
                  bcs%fs(ix,:,iqty)=bcs%fs(ix,:,iqty)*yfac
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set periodicity.
c-----------------------------------------------------------------------
      bcs%periodic=(/endmode1 == "periodic",endmode2 == "periodic"/)
      IF(bcs%periodic(1))bcs%fs(bcs%mx,:,:)=bcs%fs(0,:,:)
      IF(bcs%periodic(2))bcs%fs(:,bcs%my,:)=bcs%fs(:,0,:)
c-----------------------------------------------------------------------
c     evaluate y derivatives.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%my,bcs%mx+1)
      spl%xs=bcs%ys
      DO iqty=1,bcs%nqty
         spl%fs=TRANSPOSE(bcs%fs(:,:,iqty))
         CALL spline_fit(spl,endmode2)
         bcs%fsy(:,:,iqty)=TRANSPOSE(spl%fs1)
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     evaluate x derivatives.
c-----------------------------------------------------------------------
      spl%mx=bcs%mx
      spl%nqty=bcs%my+1
      CALL spline_alloc(spl,bcs%mx,bcs%my+1)
      spl%xs=bcs%xs
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fs(:,:,iqty)
         CALL spline_fit(spl,endmode1)
         bcs%fsx(:,:,iqty)=spl%fs1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate mixed derivatives.
c-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fsy(:,:,iqty)
         CALL spline_fit(spl,endmode1)
         bcs%fsxy(:,:,iqty)=spl%fs1
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_fit
c-----------------------------------------------------------------------
c     subprogram 4. bicube_lsfit.
c     least-square fit to cubic splines of piecewise-constant functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_lsfit(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ix,iy,nx,ny,mx,my,nrhs,n,info,i,j,k,l,kd,ldab
      INTEGER, DIMENSION(4*(bcs%mx+1)*(bcs%my+1)) :: ipiv
      REAL(r8), PARAMETER ::
     $     a0=13/35._r8,a1=9/70._r8,b0=11/210._r8,b1=13/420._r8
      REAL(r8), DIMENSION(0:bcs%mx+1) :: dx
      REAL(r8), DIMENSION(0:bcs%my+1) :: dy
      REAL(r8), DIMENSION(4,0:bcs%mx,0:bcs%my,bcs%nqty) :: rhs
      REAL(r8), DIMENSION(0:bcs%mx+1,0:bcs%my+1,bcs%nqty) :: g
      REAL(r8), DIMENSION(4,4,-1:1,-1:1,0:bcs%mx,0:bcs%my) :: amat
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ab
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nx=bcs%mx
      ny=bcs%my
      nrhs=bcs%nqty
      kd=4*(nx+1)+7
      ldab=3*kd+1
      n=4*(nx+1)*(ny+1)
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      amat=0
      rhs=0
      dx=0
      dy=0
      g=0
      dx(1:nx)=bcs%xs(1:nx)-bcs%xs(0:nx-1)
      dy(1:ny)=bcs%ys(1:ny)-bcs%ys(0:ny-1)
      g(1:nx,1:ny,:)=bcs%fs(1:nx,1:ny,:)
c-----------------------------------------------------------------------
c     least squares fit, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,1,0,0,0:nx,iy)=a0**2
     $        *(dx(0:nx)+dx(1:nx+1))*(dy(iy)+dy(iy+1))
         amat(1,1,-1,0,0:nx,iy)=a0*a1*dx(0:nx)*(dy(iy)+dy(iy+1))
         amat(1,1,+1,0,0:nx,iy)=a0*a1*dx(1:nx+1)*(dy(iy)+dy(iy+1))
         amat(1,1,0,-1,0:nx,iy)=a0*a1*(dx(0:nx)+dx(1:nx+1))*dy(iy)
         amat(1,1,0,+1,0:nx,iy)=a0*a1*(dx(0:nx)+dx(1:nx+1))*dy(iy+1)
         amat(1,1,-1,-1,0:nx,iy)=a1**2*dx(0:nx)*dy(iy)
         amat(1,1,-1,+1,0:nx,iy)=a1**2*dx(0:nx)*dy(iy+1)
         amat(1,1,+1,-1,0:nx,iy)=a1**2*dx(1:nx+1)*dy(iy)
         amat(1,1,+1,+1,0:nx,iy)=a1**2*dx(1:nx+1)*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     least squares fit, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,2,0,0,0:nx,iy)=a0*b0
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy)+dy(iy+1))
         amat(1,2,-1,0,0:nx,iy)=+a0*b1*(dy(iy)+dy(iy+1))*dx(0:nx)**2
         amat(1,2,+1,0,0:nx,iy)=-a0*b1*(dy(iy)+dy(iy+1))*dx(1:nx+1)**2
         amat(1,2,-1,-1,0:nx,iy)=+a1*b1*dx(0:nx)**2*dy(iy)
         amat(1,2,+1,-1,0:nx,iy)=-a1*b1*dx(1:nx+1)**2*dy(iy)
         amat(1,2,-1,+1,0:nx,iy)=+a1*b1*dx(0:nx)**2*dy(iy+1)
         amat(1,2,+1,+1,0:nx,iy)=-a1*b1*dx(1:nx+1)**2*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     least squares fit, y derivatiives.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(1,3,0,0,ix,0:ny)=a0*b0
     $        *(dy(1:ny+1)**2-dy(0:ny)**2)*(dx(ix)+dx(ix+1))
         amat(1,3,-1,0,ix,0:ny)=+a0*b1*(dx(ix)+dx(ix+1))*dy(0:ny)**2
         amat(1,3,+1,0,ix,0:ny)=-a0*b1*(dx(ix)+dx(ix+1))*dy(1:ny+1)**2
         amat(1,3,-1,-1,ix,0:ny)=+a1*b1*dy(0:ny)**2*dx(ix)
         amat(1,3,+1,-1,ix,0:ny)=-a1*b1*dy(1:ny+1)**2*dx(ix)
         amat(1,3,-1,+1,ix,0:ny)=+a1*b1*dy(0:ny)**2*dx(ix+1)
         amat(1,3,+1,+1,ix,0:ny)=-a1*b1*dy(1:ny+1)**2*dx(ix+1)
      ENDDO
c-----------------------------------------------------------------------
c     least squares fit, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,4,0,0,0:nx,iy)=b0**2
     $        *(dx(0:nx)**2-dx(1:nx+1)**2)*(dy(iy)**2-dy(iy+1)**2)
         amat(1,4,-1,0,0:nx,iy)=+b0*b1*dx(0:nx)**2
     $        *(dy(iy)**2+dy(iy+1)**2)
         amat(1,4,+1,0,0:nx,iy)=-b0*b1*dx(1:nx+1)**2
     $        *(dy(iy)**2+dy(iy+1)**2)
         amat(1,4,0,-1,0:nx,iy)=+b0*b1*(dx(0:nx)**2+dx(1:nx+1)**2)**2
     $        *dy(iy)**2
         amat(1,4,0,+1,0:nx,iy)=-b0*b1*(dx(0:nx)**2+dx(1:nx+1)**2)**2
     $        *dy(iy+1)
         amat(1,4,-1,-1,0:nx,iy)=b1**2*dx(0:nx)**2*dy(iy)**2
         amat(1,4,-1,+1,0:nx,iy)=-b1**2*dx(0:nx)**2*dy(iy+1)**2
         amat(1,4,+1,-1,0:nx,iy)=-b1**2*dx(1:nx+1)**2*dy(iy)**2
         amat(1,4,+1,+1,0:nx,iy)=b1**2*dx(1:nx+1)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     least squares fit, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx
            rhs(1,ix,iy,:)
     $           =(dx(ix)*g(ix,iy,:)+dx(ix+1)*g(ix+1,iy,:))*dy(iy)
     $           +(dx(ix)*g(ix,iy+1,:)+dx(ix+1)*g(ix+1,iy+1,:))*dy(iy+1)
         ENDDO
      ENDDO
      rhs=rhs/4
c-----------------------------------------------------------------------
c     continuity of second x-derivatives.
c-----------------------------------------------------------------------
      DO ix=1,nx-1
         amat(2,1,-1,0,ix,0:ny)=3/dx(ix)**2
         amat(2,1,0,0,ix,0:ny)=3/dx(ix+1)**2-3/dx(ix)**2
         amat(2,1,1,0,ix,0:ny)=-3/dx(ix+1)**2
         amat(2,2,-1,0,ix,0:ny)=1/dx(ix)
         amat(2,2,0,0,ix,0:ny)=2/dx(ix)+2/dx(ix+1)
         amat(2,2,1,0,ix,0:ny)=1/dx(ix+1)
      ENDDO
      amat(2,2,0,0,0:nx:nx,0:ny)=1
c-----------------------------------------------------------------------
c     continuity of second y-derivatives.
c-----------------------------------------------------------------------
      DO iy=1,ny-1
         amat(3,1,0,-1,0:nx,iy)=3/dy(iy)**2
         amat(3,1,0,0,0:nx,iy)=3/dy(iy+1)**2-3/dy(iy)**2
         amat(3,1,0,1,0:nx,iy)=-3/dy(iy+1)**2
         amat(3,3,0,-1,0:nx,iy)=1/dy(iy)
         amat(3,3,0,0,0:nx,iy)=2/dy(iy)+2/dy(iy+1)
         amat(3,3,0,1,0:nx,iy)=1/dy(iy+1)
      ENDDO
      amat(3,3,0,0,0:nx,0:ny:ny)=1
c-----------------------------------------------------------------------
c     continuity of mixed second derivatives.
c-----------------------------------------------------------------------
      DO ix=1,nx-1
         amat(4,3,-1,0,ix,0:ny)=3/dx(ix)**2
         amat(4,3,0,0,ix,0:ny)=3/dx(ix+1)**2-3/dx(ix)**2
         amat(4,3,1,0,ix,0:ny)=-3/dx(ix+1)**2
         amat(4,4,-1,0,ix,0:ny)=1/dx(ix)
         amat(4,4,0,0,ix,0:ny)=2/dx(ix)+2/dx(ix+1)
         amat(4,4,1,0,ix,0:ny)=1/dx(ix+1)
      ENDDO
      amat(4,4,0,0,0:nx:nx,0:ny)=1
c-----------------------------------------------------------------------
c     transfer matrix to lapack band storage.
c-----------------------------------------------------------------------
      ALLOCATE(ab(ldab,n))
      ab=0
      DO ix=0,nx
         DO mx=MAX(-ix,-1),MIN(nx-ix,1)
            DO iy=0,ny
               DO my=MAX(-iy,-1),MIN(ny-iy,1)
                  DO k=1,4
                     DO l=1,4
                        i=4*(iy*(nx+1)+ix)+k
                        j=4*((iy+my)*(nx+1)+ix+mx)+l
                        ab(2*kd+1+i-j,j)=amat(k,l,mx,my,ix,iy)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     factor and solve.
c-----------------------------------------------------------------------
      CALL dgbtrf(n,n,kd,kd,ab,ldab,ipiv,info)
      CALL dgbtrs('N',n,kd,kd,nrhs,ab,ldab,ipiv,rhs,n,info)
      DEALLOCATE(ab)
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      bcs%fs=rhs(1,:,:,:)
      bcs%fsx=rhs(2,:,:,:)
      bcs%fsy=rhs(3,:,:,:)
      bcs%fsxy=rhs(4,:,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_lsfit
c-----------------------------------------------------------------------
c     subprogram 5. bicube_eval.
c     evaluates bicubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_eval(bcs,x,y,mode)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      REAL(r8), INTENT(IN) :: x,y
      INTEGER, INTENT(IN) :: mode

      INTEGER :: i,iqty,iside
      REAL(r8) :: dx,dy,xx,yy,g,gx,gy,gxx,gyy,gxy,xfac,yfac
      REAL(r8), DIMENSION (4,4,bcs%nqty) :: c
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      IF(bcs%mx==0 .OR. bcs%my==0)THEN
         WRITE(*, *) 'ERROR: Bicubic spline entitled "',bcs%title,'"'
         IF(bcs%mx==0)THEN
            WRITE(*, *) 'has 0 elements specified for x axis '
     $                  //TRIM(bcs%xtitle)
         ELSE
            WRITE(*, *) 'has 0 elements specified for y axis '
     $                  //TRIM(bcs%ytitle)
         ENDIF
         STOP
      ENDIF
      bcs%ix=max(bcs%ix,0)
      bcs%ix=min(bcs%ix,bcs%mx-1)
      bcs%iy=max(bcs%iy,0)
      bcs%iy=min(bcs%iy,bcs%my-1)
      xx=x
      yy=y
c-----------------------------------------------------------------------
c     normalize x interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(1))THEN
         DO
            IF(xx < bcs%xs(bcs%mx))EXIT
            xx=xx-bcs%xs(bcs%mx)
         ENDDO
         DO
            IF(xx >= bcs%xs(0))EXIT
            xx=xx+bcs%xs(bcs%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find x interval.
c-----------------------------------------------------------------------
      DO
         IF(bcs%ix <= 0)EXIT
         IF(xx >= bcs%xs(bcs%ix))EXIT
         bcs%ix=bcs%ix-1
      ENDDO
      DO
         IF(bcs%ix >= bcs%mx-1)EXIT
         IF(xx < bcs%xs(bcs%ix+1))EXIT
         bcs%ix=bcs%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     normalize y interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(2))THEN
         DO
            IF(yy < bcs%ys(bcs%my))EXIT
            yy=yy-bcs%ys(bcs%my)
         ENDDO
         DO
            IF(yy >= bcs%ys(0))EXIT
            yy=yy+bcs%ys(bcs%my)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find y interval.
c-----------------------------------------------------------------------
      DO
         IF(bcs%iy <= 0)EXIT
         IF(yy >= bcs%ys(bcs%iy))EXIT
         bcs%iy=bcs%iy-1
      ENDDO
      DO
         IF(bcs%iy >= bcs%my-1)EXIT
         IF(yy < bcs%ys(bcs%iy+1))EXIT
         bcs%iy=bcs%iy+1
      ENDDO
c-----------------------------------------------------------------------
c     find offsets and compute local coefficients.
c-----------------------------------------------------------------------
      dx=xx-bcs%xs(bcs%ix)
      dy=yy-bcs%ys(bcs%iy)
      IF(ALLOCATED(bcs%cmats))THEN
         c=bcs%cmats(:,:,bcs%ix+1,bcs%iy+1,:)
      ELSE
         c=bicube_getco(bcs)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      bcs%f=0
      DO i=4,1,-1
         bcs%f=bcs%f*dx
     $        +((c(i,4,:)*dy
     $        +c(i,3,:))*dy
     $        +c(i,2,:))*dy
     $        +c(i,1,:)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate first derivatives of f
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         bcs%fx=0
         bcs%fy=0
         DO i=4,1,-1
            bcs%fy=bcs%fy*dx
     $           +(c(i,4,:)*3*dy
     $           +c(i,3,:)*2)*dy
     $           +c(i,2,:)
            bcs%fx=bcs%fx*dy
     $           +(c(4,i,:)*3*dx
     $           +c(3,i,:)*2)*dx
     $           +c(2,i,:)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives of f
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         bcs%fxx=0
         bcs%fyy=0
         bcs%fxy=0
         DO i=4,1,-1
            bcs%fyy=bcs%fyy*dx
     $           +(c(i,4,:)*3*dy
     $           +c(i,3,:))*2
            bcs%fxx=bcs%fxx*dy
     $           +(c(4,i,:)*3*dx
     $           +c(3,i,:))*2
         ENDDO
         DO i=4,2,-1
            bcs%fxy=bcs%fxy*dx
     $           +((c(i,4,:)*3*dy
     $           +c(i,3,:)*2)*dy
     $           +c(i,2,:))*(i-1)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-bcs%x0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**bcs%xpower(iside,iqty)
            g=bcs%f(iqty)*xfac
            IF(mode > 0)THEN
               gx=(bcs%fx(iqty)+bcs%f(iqty)
     $              *bcs%xpower(iside,iqty)/dx)*xfac
               gy=bcs%fy(iqty)*xfac
            ENDIF
            IF(mode > 1)THEN
               gxx=(bcs%fxx(iqty)+bcs%xpower(iside,iqty)/dx
     $              *(2*bcs%fx(iqty)+(bcs%xpower(iside,iqty)-1)
     $              *bcs%f(iqty)/dx))*xfac
               gxy=(bcs%fxy(iqty)+bcs%fy(iqty)
     $              *bcs%xpower(iside,iqty)/dx)*xfac
               gyy=bcs%fyy(iqty)*xfac
            ENDIF
            bcs%f(iqty)=g
            IF(mode > 0)THEN
               bcs%fx(iqty)=gx
               bcs%fy(iqty)=gy
            ENDIF
            IF(mode > 1)THEN
               bcs%fxx(iqty)=gxx
               bcs%fxy(iqty)=gxy
               bcs%fyy(iqty)=gyy
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     restore y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dy=y-bcs%y0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) == 0)CYCLE
            yfac=ABS(dy)**bcs%ypower(iside,iqty)
            g=bcs%f(iqty)*yfac
            IF(mode > 0)THEN
               gx=bcs%fx(iqty)*yfac
               gy=(bcs%fy(iqty)+bcs%f(iqty)
     $              *bcs%ypower(iside,iqty)/dy)*yfac
            ENDIF
            IF(mode > 1)THEN
               gxx=bcs%fxx(iqty)*yfac
               gxy=(bcs%fxy(iqty)+bcs%fy(iqty)
     $              *bcs%ypower(iside,iqty)/dy)*yfac
               gyy=(bcs%fyy(iqty)+bcs%ypower(iside,iqty)/dy
     $              *(2*bcs%fy(iqty)+(bcs%ypower(iside,iqty)-1)
     $              *bcs%f(iqty)/dy))*yfac
            ENDIF
            bcs%f(iqty)=g
            IF(mode > 0)THEN
               bcs%fx(iqty)=gx
               bcs%fy(iqty)=gy
            ENDIF
            IF(mode > 1)THEN
               bcs%fxx(iqty)=gxx
               bcs%fxy(iqty)=gxy
               bcs%fyy(iqty)=gyy
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_eval
c-----------------------------------------------------------------------
c     subprogram 5a. bicube_eval_external.
c     evaluates bicubic spline function with external arrays (parallel).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_eval_external(bcs, x, y, mode,
     $     b_ix, b_iy, b_f, b_fx, b_fy)

      TYPE(bicube_type), INTENT(IN) :: bcs
      REAL(r8), INTENT(IN) :: x,y
      INTEGER, INTENT(IN) :: mode

      INTEGER :: i,iqty,iside
      REAL(r8) :: dx,dy,xx,yy,g,gx,gy,gxx,gyy,gxy,xfac,yfac
      REAL(r8), DIMENSION (4,4,bcs%nqty) :: c

      INTEGER, INTENT(INOUT) :: b_ix,b_iy
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: b_f,b_fx,b_fy
c-----------------------------------------------------------------------
c     error-check for mode number--external array is limited.
c-----------------------------------------------------------------------
      IF (mode > 1) THEN
         CALL program_stop("Set bicube_eval_external mode <=1 !")
      ENDIF
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      b_ix=max(b_ix,0)
      b_ix=min(b_ix,bcs%mx-1)
      b_iy=max(b_iy,0)
      b_iy=min(b_iy,bcs%my-1)
      xx=x
      yy=y
c-----------------------------------------------------------------------
c     normalize x interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(1))THEN
         DO
            IF(xx < bcs%xs(bcs%mx))EXIT
            xx=xx-bcs%xs(bcs%mx)
         ENDDO
         DO
            IF(xx >= bcs%xs(0))EXIT
            xx=xx+bcs%xs(bcs%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find x interval.
c-----------------------------------------------------------------------
      DO
         IF(b_ix <= 0)EXIT
         IF(xx >= bcs%xs(b_ix))EXIT
         b_ix=b_ix-1
      ENDDO
      DO
         IF(b_ix >= bcs%mx-1)EXIT
         IF(xx < bcs%xs(b_ix+1))EXIT
         b_ix=b_ix+1
      ENDDO
c-----------------------------------------------------------------------
c     normalize y interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(2))THEN
         DO
            IF(yy < bcs%ys(bcs%my))EXIT
            yy=yy-bcs%ys(bcs%my)
         ENDDO
         DO
            IF(yy >= bcs%ys(0))EXIT
            yy=yy+bcs%ys(bcs%my)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find y interval.
c-----------------------------------------------------------------------
      DO
         IF(b_iy <= 0)EXIT
         IF(yy >= bcs%ys(b_iy))EXIT
         b_iy=b_iy-1
      ENDDO
      DO
         IF(b_iy >= bcs%my-1)EXIT
         IF(yy < bcs%ys(b_iy+1))EXIT
         b_iy=b_iy+1
      ENDDO
c-----------------------------------------------------------------------
c     find offsets and compute local coefficients.
c-----------------------------------------------------------------------
      dx=xx-bcs%xs(b_ix)
      dy=yy-bcs%ys(b_iy)
      call bicube_getco_external_sub(bcs,b_ix,b_iy,c)
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      b_f=0
      DO i=4,1,-1
         b_f=b_f*dx
     $        +((c(i,4,:)*dy
     $        +c(i,3,:))*dy
     $        +c(i,2,:))*dy
     $        +c(i,1,:)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate first derivatives of f
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         b_fx=0
         b_fy=0
         DO i=4,1,-1
            b_fy=b_fy*dx
     $           +(c(i,4,:)*3*dy
     $           +c(i,3,:)*2)*dy
     $           +c(i,2,:)
            b_fx=b_fx*dy
     $           +(c(4,i,:)*3*dx
     $           +c(3,i,:)*2)*dx
     $           +c(2,i,:)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-bcs%x0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**bcs%xpower(iside,iqty)
            g=b_f(iqty)*xfac
            IF(mode > 0)THEN
               gx=(b_fx(iqty)+b_f(iqty)
     $              *bcs%xpower(iside,iqty)/dx)*xfac
               gy=b_fy(iqty)*xfac
            ENDIF
            b_f(iqty)=g
            IF(mode > 0)THEN
               b_fx(iqty)=gx
               b_fy(iqty)=gy
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     restore y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dy=y-bcs%y0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) == 0)CYCLE
            yfac=ABS(dy)**bcs%ypower(iside,iqty)
            g=b_f(iqty)*yfac
            IF(mode > 0)THEN
               gx=b_fx(iqty)*yfac
               gy=(b_fy(iqty)+b_f(iqty)
     $              *bcs%ypower(iside,iqty)/dy)*yfac
            ENDIF
            b_f(iqty)=g
            IF(mode > 0)THEN
               b_fx(iqty)=gx
               b_fy(iqty)=gy
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_eval_external
c-----------------------------------------------------------------------
c     subprogram 6. bicube_getco.
c     computes coefficient matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bicube_getco(bcs) RESULT(cmat)

      TYPE(bicube_type), INTENT(IN) :: bcs
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: cmat

      REAL(r8) :: hxfac,hxfac2,hxfac3
      REAL(r8) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(bcs%ix+1)-bcs%xs(bcs%ix))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1)=-3*hxfac2
      gxmat(3,2)=-2*hxfac
      gxmat(3,3)=3*hxfac2
      gxmat(3,4)=-hxfac
      gxmat(4,1)=2*hxfac3
      gxmat(4,2)=hxfac2
      gxmat(4,3)=-2*hxfac3
      gxmat(4,4)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(bcs%iy+1)-bcs%ys(bcs%iy))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1)=-3*hyfac2
      gymat(3,2)=-2*hyfac
      gymat(3,3)=3*hyfac2
      gymat(3,4)=-hyfac
      gymat(4,1)=2*hyfac3
      gymat(4,2)=hyfac2
      gymat(4,3)=-2*hyfac3
      gymat(4,4)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(bcs%ix,bcs%iy,:)
      cmat(1,2,:)=bcs%fsy(bcs%ix,bcs%iy,:)
      cmat(1,3,:)=bcs%fs(bcs%ix,bcs%iy+1,:)
      cmat(1,4,:)=bcs%fsy(bcs%ix,bcs%iy+1,:)
      cmat(2,1,:)=bcs%fsx(bcs%ix,bcs%iy,:)
      cmat(2,2,:)=bcs%fsxy(bcs%ix,bcs%iy,:)
      cmat(2,3,:)=bcs%fsx(bcs%ix,bcs%iy+1,:)
      cmat(2,4,:)=bcs%fsxy(bcs%ix,bcs%iy+1,:)
      cmat(3,1,:)=bcs%fs(bcs%ix+1,bcs%iy,:)
      cmat(3,2,:)=bcs%fsy(bcs%ix+1,bcs%iy,:)
      cmat(3,3,:)=bcs%fs(bcs%ix+1,bcs%iy+1,:)
      cmat(3,4,:)=bcs%fsy(bcs%ix+1,bcs%iy+1,:)
      cmat(4,1,:)=bcs%fsx(bcs%ix+1,bcs%iy,:)
      cmat(4,2,:)=bcs%fsxy(bcs%ix+1,bcs%iy,:)
      cmat(4,3,:)=bcs%fsx(bcs%ix+1,bcs%iy+1,:)
      cmat(4,4,:)=bcs%fsxy(bcs%ix+1,bcs%iy+1,:)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:)
      temp(:,3,:)
     $     =cmat(:,1,:)*gymat(3,1)
     $     +cmat(:,2,:)*gymat(3,2)
     $     +cmat(:,3,:)*gymat(3,3)
     $     +cmat(:,4,:)*gymat(3,4)
      temp(:,4,:)
     $     =cmat(:,1,:)*gymat(4,1)
     $     +cmat(:,2,:)*gymat(4,2)
     $     +cmat(:,3,:)*gymat(4,3)
     $     +cmat(:,4,:)*gymat(4,4)
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      cmat(1:2,:,:)=temp(1:2,:,:)
      cmat(3,:,:)
     $     =gxmat(3,1)*temp(1,:,:)
     $     +gxmat(3,2)*temp(2,:,:)
     $     +gxmat(3,3)*temp(3,:,:)
     $     +gxmat(3,4)*temp(4,:,:)
      cmat(4,:,:)
     $     =gxmat(4,1)*temp(1,:,:)
     $     +gxmat(4,2)*temp(2,:,:)
     $     +gxmat(4,3)*temp(3,:,:)
     $     +gxmat(4,4)*temp(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bicube_getco
c-----------------------------------------------------------------------
c     subprogram 6a. bicube_getco_external.
c     computes coefficient matrices for external arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bicube_getco_external(bcs,b_ix,b_iy) RESULT(cmat)

      TYPE(bicube_type), INTENT(IN) :: bcs
      INTEGER, INTENT(IN) :: b_ix, b_iy
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: cmat

      REAL(r8) :: hxfac,hxfac2,hxfac3
      REAL(r8) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(b_ix+1)-bcs%xs(b_ix))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1)=-3*hxfac2
      gxmat(3,2)=-2*hxfac
      gxmat(3,3)=3*hxfac2
      gxmat(3,4)=-hxfac
      gxmat(4,1)=2*hxfac3
      gxmat(4,2)=hxfac2
      gxmat(4,3)=-2*hxfac3
      gxmat(4,4)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(b_iy+1)-bcs%ys(b_iy))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1)=-3*hyfac2
      gymat(3,2)=-2*hyfac
      gymat(3,3)=3*hyfac2
      gymat(3,4)=-hyfac
      gymat(4,1)=2*hyfac3
      gymat(4,2)=hyfac2
      gymat(4,3)=-2*hyfac3
      gymat(4,4)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(b_ix,b_iy,:)
      cmat(1,2,:)=bcs%fsy(b_ix,b_iy,:)
      cmat(1,3,:)=bcs%fs(b_ix,b_iy+1,:)
      cmat(1,4,:)=bcs%fsy(b_ix,b_iy+1,:)
      cmat(2,1,:)=bcs%fsx(b_ix,b_iy,:)
      cmat(2,2,:)=bcs%fsxy(b_ix,b_iy,:)
      cmat(2,3,:)=bcs%fsx(b_ix,b_iy+1,:)
      cmat(2,4,:)=bcs%fsxy(b_ix,b_iy+1,:)
      cmat(3,1,:)=bcs%fs(b_ix+1,b_iy,:)
      cmat(3,2,:)=bcs%fsy(b_ix+1,b_iy,:)
      cmat(3,3,:)=bcs%fs(b_ix+1,b_iy+1,:)
      cmat(3,4,:)=bcs%fsy(b_ix+1,b_iy+1,:)
      cmat(4,1,:)=bcs%fsx(b_ix+1,b_iy,:)
      cmat(4,2,:)=bcs%fsxy(b_ix+1,b_iy,:)
      cmat(4,3,:)=bcs%fsx(b_ix+1,b_iy+1,:)
      cmat(4,4,:)=bcs%fsxy(b_ix+1,b_iy+1,:)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:)
      temp(:,3,:)
     $     =cmat(:,1,:)*gymat(3,1)
     $     +cmat(:,2,:)*gymat(3,2)
     $     +cmat(:,3,:)*gymat(3,3)
     $     +cmat(:,4,:)*gymat(3,4)
      temp(:,4,:)
     $     =cmat(:,1,:)*gymat(4,1)
     $     +cmat(:,2,:)*gymat(4,2)
     $     +cmat(:,3,:)*gymat(4,3)
     $     +cmat(:,4,:)*gymat(4,4)
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      cmat(1:2,:,:)=temp(1:2,:,:)
      cmat(3,:,:)
     $     =gxmat(3,1)*temp(1,:,:)
     $     +gxmat(3,2)*temp(2,:,:)
     $     +gxmat(3,3)*temp(3,:,:)
     $     +gxmat(3,4)*temp(4,:,:)
      cmat(4,:,:)
     $     =gxmat(4,1)*temp(1,:,:)
     $     +gxmat(4,2)*temp(2,:,:)
     $     +gxmat(4,3)*temp(3,:,:)
     $     +gxmat(4,4)*temp(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bicube_getco_external
c-----------------------------------------------------------------------
c     subprogram 6b. bicube_getco_external_sub.
c     computes coefficient matrices for external arrays.
c     Subroutine version.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_getco_external_sub(bcs,b_ix,b_iy, cmat)

      TYPE(bicube_type), INTENT(IN) :: bcs
      INTEGER, INTENT(IN) :: b_ix, b_iy
      REAL(r8), DIMENSION(4,4,bcs%nqty), INTENT(OUT) :: cmat

      REAL(r8) :: hxfac,hxfac2,hxfac3
      REAL(r8) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(b_ix+1)-bcs%xs(b_ix))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1)=-3*hxfac2
      gxmat(3,2)=-2*hxfac
      gxmat(3,3)=3*hxfac2
      gxmat(3,4)=-hxfac
      gxmat(4,1)=2*hxfac3
      gxmat(4,2)=hxfac2
      gxmat(4,3)=-2*hxfac3
      gxmat(4,4)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(b_iy+1)-bcs%ys(b_iy))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1)=-3*hyfac2
      gymat(3,2)=-2*hyfac
      gymat(3,3)=3*hyfac2
      gymat(3,4)=-hyfac
      gymat(4,1)=2*hyfac3
      gymat(4,2)=hyfac2
      gymat(4,3)=-2*hyfac3
      gymat(4,4)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(b_ix,b_iy,:)
      cmat(1,2,:)=bcs%fsy(b_ix,b_iy,:)
      cmat(1,3,:)=bcs%fs(b_ix,b_iy+1,:)
      cmat(1,4,:)=bcs%fsy(b_ix,b_iy+1,:)
      cmat(2,1,:)=bcs%fsx(b_ix,b_iy,:)
      cmat(2,2,:)=bcs%fsxy(b_ix,b_iy,:)
      cmat(2,3,:)=bcs%fsx(b_ix,b_iy+1,:)
      cmat(2,4,:)=bcs%fsxy(b_ix,b_iy+1,:)
      cmat(3,1,:)=bcs%fs(b_ix+1,b_iy,:)
      cmat(3,2,:)=bcs%fsy(b_ix+1,b_iy,:)
      cmat(3,3,:)=bcs%fs(b_ix+1,b_iy+1,:)
      cmat(3,4,:)=bcs%fsy(b_ix+1,b_iy+1,:)
      cmat(4,1,:)=bcs%fsx(b_ix+1,b_iy,:)
      cmat(4,2,:)=bcs%fsxy(b_ix+1,b_iy,:)
      cmat(4,3,:)=bcs%fsx(b_ix+1,b_iy+1,:)
      cmat(4,4,:)=bcs%fsxy(b_ix+1,b_iy+1,:)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:)
      temp(:,3,:)
     $     =cmat(:,1,:)*gymat(3,1)
     $     +cmat(:,2,:)*gymat(3,2)
     $     +cmat(:,3,:)*gymat(3,3)
     $     +cmat(:,4,:)*gymat(3,4)
        temp(:,4,:)
     $     =cmat(:,1,:)*gymat(4,1)
     $     +cmat(:,2,:)*gymat(4,2)
     $     +cmat(:,3,:)*gymat(4,3)
     $     +cmat(:,4,:)*gymat(4,4)
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      cmat(1:2,:,:)=temp(1:2,:,:)
      cmat(3,:,:)
     $     =gxmat(3,1)*temp(1,:,:)
     $     +gxmat(3,2)*temp(2,:,:)
     $     +gxmat(3,3)*temp(3,:,:)
     $     +gxmat(3,4)*temp(4,:,:)
         cmat(4,:,:)
     $     =gxmat(4,1)*temp(1,:,:)
     $     +gxmat(4,2)*temp(2,:,:)
     $     +gxmat(4,3)*temp(3,:,:)
     $     +gxmat(4,4)*temp(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_getco_external_sub
c-----------------------------------------------------------------------
c     subprogram 7. bicube_all_eval.
c     evaluates bicubic splines in all intervals.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_all_eval(bcs,dx,dy,f,fx,fy,fxx,fyy,fxy,mode)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      REAL(r8), INTENT(IN) :: dx,dy
      REAL(r8), INTENT(OUT), DIMENSION(bcs%mx,bcs%my,bcs%nqty) ::
     $     f,fx,fy,fxx,fyy,fxy
      INTEGER, INTENT(IN) :: mode

      INTEGER :: i,ix,iy,iqty,iside
      REAL(r8), DIMENSION(bcs%mx) :: dxv
      REAL(r8), DIMENSION(bcs%my) :: dyv

      REAL(R8), DIMENSION(bcs%mx) :: dxx,xfac
      REAL(R8), DIMENSION(bcs%my) :: dyy,yfac
      REAL(R8), DIMENSION(bcs%mx,bcs%my) :: g,gx,gy,gxx,gxy,gyy
c-----------------------------------------------------------------------
c     compute local displacements and coefficients.
c-----------------------------------------------------------------------
      dxv=(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1))*dx
      dyv=(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1))*dy
      CALL bicube_all_getco(bcs)
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      f=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO ix=1,bcs%mx
               f(ix,:,:)=f(ix,:,:)*dxv(ix)
            ENDDO
         ENDIF
         DO iy=1,bcs%my
            f(:,iy,:)=f(:,iy,:)
     $           +((bcs%cmats(i,4,:,iy,:)*dyv(iy)
     $           +bcs%cmats(i,3,:,iy,:))*dyv(iy)
     $           +bcs%cmats(i,2,:,iy,:))*dy
     $           +bcs%cmats(i,1,:,iy,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     evaluate fx.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         fx=0
         DO i=4,1,-1
            IF(i /= 4)THEN
               DO iy=1,bcs%my
                  fx(:,iy,:)=fx(:,iy,:)*dyv(iy)
               ENDDO
            ENDIF
            DO ix=1,bcs%mx
               fx(ix,:,:)=fx(ix,:,:)
     $              +(bcs%cmats(4,i,ix,:,:)*3*dxv(ix)
     $              +bcs%cmats(3,i,ix,:,:)*2)*dxv(ix)
     $              +bcs%cmats(2,i,ix,:,:)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate fy.
c-----------------------------------------------------------------------
         fy=0
         DO i=4,1,-1
            IF(i /= 4)THEN
               DO ix=1,bcs%mx
                  fy(ix,:,:)=fy(ix,:,:)*dxv(ix)
               ENDDO
            ENDIF
            DO iy=1,bcs%my
               fy(:,iy,:)=fy(:,iy,:)
     $              +(bcs%cmats(i,4,:,iy,:)*3*dyv(iy)
     $              +bcs%cmats(i,3,:,iy,:)*2)*dyv(iy)
     $              +bcs%cmats(i,2,:,iy,:)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate fxx.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         fxx=0
         DO i=4,1,-1
            IF(i /= 4)THEN
               DO iy=1,bcs%my
                  fxx(:,iy,:)=fxx(:,iy,:)*dyv(iy)
               ENDDO
            ENDIF
            DO ix=1,bcs%mx
               fxx(ix,:,:)=fxx(ix,:,:)
     $              +(bcs%cmats(4,i,ix,:,:)*3*dxv(ix)
     $              +bcs%cmats(3,i,ix,:,:))*2
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate fyy and fxy
c-----------------------------------------------------------------------
         fyy=0
         DO i=4,1,-1
            IF(i /= 4)THEN
               DO ix=1,bcs%mx
                  fyy(ix,:,:)=fyy(ix,:,:)*dxv(ix)
                  fxy(ix,:,:)=fxy(ix,:,:)*dxv(ix)
               ENDDO
            ENDIF
            DO iy=1,bcs%my
               fyy(:,iy,:)=fyy(:,iy,:)
     $              +(bcs%cmats(i,4,:,iy,:)*3*dyv(iy)
     $              +bcs%cmats(i,3,:,iy,:))*2
               fxy(:,iy,:)=fxy(:,iy,:)
     $              +((bcs%cmats(i,4,:,iy,:)*3*dyv(iy)
     $              +bcs%cmats(i,3,:,iy,:)*2)*dyv(iy)
     $              +bcs%cmats(i,2,:,iy,:))*(i-1)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dxx=(bcs%xs(0:bcs%mx-1)+dxv(1:bcs%mx))-bcs%x0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) == 0)CYCLE
            xfac=dxx**bcs%xpower(iside,iqty)
            DO iy=1,bcs%my
               g(:,iy)=f(:,iy,iqty)*xfac
               IF(mode > 0)THEN
                  gx(:,iy)=(fx(:,iy,iqty)+f(:,iy,iqty)
     $                 *bcs%xpower(iside,iqty)/dxx)*xfac
                  gy(:,iy)=fy(:,iy,iqty)*xfac
               ENDIF
               IF(mode > 1)THEN
                  gxy(:,iy)=(fxy(:,iy,iqty)+fy(:,iy,iqty)
     $                 *bcs%xpower(iside,iqty)/dxx)*xfac
                  gxx(:,iy)=(fxx(:,iy,iqty)+bcs%xpower(iside,iqty)/dxx
     $                 *(2*fx(:,iy,iqty)+(bcs%xpower(iside,iqty)-1)
     $                 *f(:,iy,iqty)/dxx))*xfac
                  gyy(:,iy)=fyy(:,iy,iqty)*xfac
               ENDIF
               f(:,iy,iqty)=g(:,iy)
               IF(mode > 0)THEN
                  fx(:,iy,iqty)=gx(:,iy)
                  fy(:,iy,iqty)=gy(:,iy)
               ENDIF
               IF(mode > 1)THEN
                  fxx(:,iy,iqty)=gxx(:,iy)
                  fxy(:,iy,iqty)=gxy(:,iy)
                  fyy(:,iy,iqty)=gyy(:,iy)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     restore y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dyy=(bcs%ys(0:bcs%my-1)+dyv(1:bcs%my))-bcs%y0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) == 0)CYCLE
            yfac=dyy**bcs%ypower(iside,iqty)
            DO ix=1,bcs%mx
               g(ix,:)=f(ix,:,iqty)*yfac
               IF(mode > 0)THEN
                  gy(ix,:)=(fy(ix,:,iqty)+f(ix,:,iqty)
     $                 *bcs%ypower(iside,iqty)/dyy)*yfac
                  gx(ix,:)=fx(ix,:,iqty)*yfac
               ENDIF
               IF(mode > 1)THEN
                  gxx(ix,:)=fxx(ix,:,iqty)*yfac
                  gxy(ix,:)=(fxy(ix,:,iqty)+fx(ix,:,iqty)
     $                 *bcs%ypower(iside,iqty)/dyy)*yfac
                  gyy(ix,:)=(fyy(ix,:,iqty)+bcs%ypower(iside,iqty)/dyy
     $                 *(2*fy(ix,:,iqty)+(bcs%ypower(iside,iqty)-1)
     $                 *f(ix,:,iqty)/dyy))*yfac
               ENDIF
               f(ix,:,iqty)=g(ix,:)
               IF(mode > 0)THEN
                  fx(ix,:,iqty)=gx(ix,:)
                  fy(ix,:,iqty)=gy(ix,:)
               ENDIF
               IF(mode > 1)THEN
                  fxx(ix,:,iqty)=gxx(ix,:)
                  fxy(ix,:,iqty)=gxy(ix,:)
                  fyy(ix,:,iqty)=gyy(ix,:)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_all_eval
c-----------------------------------------------------------------------
c     subprogram 8. bicube_all_getco.
c     computes coefficient matrices in all intervals.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_all_getco(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      INTEGER :: ix,iy
      REAL(r8), DIMENSION(bcs%mx) :: hxfac,hxfac2,hxfac3
      REAL(r8), DIMENSION(bcs%my) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4,bcs%mx) :: gxmat
      REAL(r8), DIMENSION(3:4,4,bcs%my) :: gymat
      REAL(r8), DIMENSION(4,4,bcs%mx,bcs%my,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      IF(ALLOCATED(bcs%cmats))THEN
         RETURN
      ELSE
         ALLOCATE(bcs%cmats(4,4,bcs%mx,bcs%my,bcs%nqty))
      ENDIF
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1,:)=-3*hxfac2
      gxmat(3,2,:)=-2*hxfac
      gxmat(3,3,:)=3*hxfac2
      gxmat(3,4,:)=-hxfac
      gxmat(4,1,:)=2*hxfac3
      gxmat(4,2,:)=hxfac2
      gxmat(4,3,:)=-2*hxfac3
      gxmat(4,4,:)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1,:)=-3*hyfac2
      gymat(3,2,:)=-2*hyfac
      gymat(3,3,:)=3*hyfac2
      gymat(3,4,:)=-hyfac
      gymat(4,1,:)=2*hyfac3
      gymat(4,2,:)=hyfac2
      gymat(4,3,:)=-2*hyfac3
      gymat(4,4,:)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      bcs%cmats(1,1,:,:,:)=bcs%fs(0:bcs%mx-1,0:bcs%my-1,:)
      bcs%cmats(1,2,:,:,:)=bcs%fsy(0:bcs%mx-1,0:bcs%my-1,:)
      bcs%cmats(1,3,:,:,:)=bcs%fs(0:bcs%mx-1,1:bcs%my,:)
      bcs%cmats(1,4,:,:,:)=bcs%fsy(0:bcs%mx-1,1:bcs%my,:)
      bcs%cmats(2,1,:,:,:)=bcs%fsx(0:bcs%mx-1,0:bcs%my-1,:)
      bcs%cmats(2,2,:,:,:)=bcs%fsxy(0:bcs%mx-1,0:bcs%my-1,:)
      bcs%cmats(2,3,:,:,:)=bcs%fsx(0:bcs%mx-1,1:bcs%my,:)
      bcs%cmats(2,4,:,:,:)=bcs%fsxy(0:bcs%mx-1,1:bcs%my,:)
      bcs%cmats(3,1,:,:,:)=bcs%fs(1:bcs%mx,0:bcs%my-1,:)
      bcs%cmats(3,2,:,:,:)=bcs%fsy(1:bcs%mx,0:bcs%my-1,:)
      bcs%cmats(3,3,:,:,:)=bcs%fs(1:bcs%mx,1:bcs%my,:)
      bcs%cmats(3,4,:,:,:)=bcs%fsy(1:bcs%mx,1:bcs%my,:)
      bcs%cmats(4,1,:,:,:)=bcs%fsx(1:bcs%mx,0:bcs%my-1,:)
      bcs%cmats(4,2,:,:,:)=bcs%fsxy(1:bcs%mx,0:bcs%my-1,:)
      bcs%cmats(4,3,:,:,:)=bcs%fsx(1:bcs%mx,1:bcs%my,:)
      bcs%cmats(4,4,:,:,:)=bcs%fsxy(1:bcs%mx,1:bcs%my,:)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:,:,:)=bcs%cmats(:,1:2,:,:,:)
      DO iy=1,bcs%my
         temp(:,3,:,iy,:)
     $        =bcs%cmats(:,1,:,iy,:)*gymat(3,1,iy)
     $        +bcs%cmats(:,2,:,iy,:)*gymat(3,2,iy)
     $        +bcs%cmats(:,3,:,iy,:)*gymat(3,3,iy)
     $        +bcs%cmats(:,4,:,iy,:)*gymat(3,4,iy)
         temp(:,4,:,iy,:)
     $        =bcs%cmats(:,1,:,iy,:)*gymat(4,1,iy)
     $        +bcs%cmats(:,2,:,iy,:)*gymat(4,2,iy)
     $        +bcs%cmats(:,3,:,iy,:)*gymat(4,3,iy)
     $        +bcs%cmats(:,4,:,iy,:)*gymat(4,4,iy)
      ENDDO
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      bcs%cmats(1:2,:,:,:,:)=temp(1:2,:,:,:,:)
      DO ix=1,bcs%mx
         bcs%cmats(3,:,ix,:,:)
     $        =gxmat(3,1,ix)*temp(1,:,ix,:,:)
     $        +gxmat(3,2,ix)*temp(2,:,ix,:,:)
     $        +gxmat(3,3,ix)*temp(3,:,ix,:,:)
     $        +gxmat(3,4,ix)*temp(4,:,ix,:,:)
         bcs%cmats(4,:,ix,:,:)
     $        =gxmat(4,1,ix)*temp(1,:,ix,:,:)
     $        +gxmat(4,2,ix)*temp(2,:,ix,:,:)
     $        +gxmat(4,3,ix)*temp(3,:,ix,:,:)
     $        +gxmat(4,4,ix)*temp(4,:,ix,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_all_getco
c-----------------------------------------------------------------------
c     subprogram 9. bicube_write_xy.
c     produces ascii and binary output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_xy(bcs,out,bin,iua,iub,interp)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      INTEGER :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      CHARACTER(80) :: format1,format2
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"iy = ",i3,", ",a6," = ",1p,e11.3)
 20   FORMAT(1x,"iy = ",i3,", jy = ",i1,", ",a6," = ",1p,e11.3)
 30   FORMAT('(/3x,"ix",4x,a,1x,',i3.3,'(4x,a6,1x)/)')
 40   FORMAT('(i5,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(format1,30)bcs%nqty
         WRITE(format2,40)bcs%nqty
      ENDIF
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out)WRITE(iua,'(1x,a/)')"input data"
      DO iy=0,bcs%my
         y=bcs%ys(iy)
         IF(out)then
            WRITE(iua,10)iy,bcs%ytitle,bcs%ys(iy)
            WRITE(iua,format1)bcs%xtitle,
     $           (bcs%title(iqty),iqty=1,bcs%nqty)
         ENDIF
         DO ix=0,bcs%mx
            x=bcs%xs(ix)
            CALL bicube_eval(bcs,x,y,0)
            IF(out)WRITE(iua,format2)ix,x,bcs%f
            IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4)
         ENDDO
         IF(out)WRITE(iua,format1)bcs%xtitle,
     $        (bcs%title(iqty),iqty=1,bcs%nqty)
         IF(bin)WRITE(iub)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
      IF(interp)THEN
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data"
         DO iy=0,bcs%my-1
            dy=(bcs%ys(iy+1)-bcs%ys(iy))/4
            DO jy=0,4
               y=bcs%ys(iy)+dy*jy
               IF(out)then
                  WRITE(iua,20)iy,jy,bcs%ytitle,y
                  WRITE(iua,format1)bcs%xtitle,
     $                 (bcs%title(iqty),iqty=1,bcs%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
               DO ix=0,bcs%mx-1
                  dx=(bcs%xs(ix+1)-bcs%xs(ix))/4
                  DO jx=0,4
                     x=bcs%xs(ix)+dx*jx
                     CALL bicube_eval(bcs,x,y,0)
                     IF(out)WRITE(iua,format2)ix,x,bcs%f
                     IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4)
                  ENDDO
                  IF(out)WRITE(iua,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over y.
c-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)bcs%xtitle,
     $                 (bcs%title(iqty),iqty=1,bcs%nqty)
               IF(bin)WRITE(iub)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_xy
c-----------------------------------------------------------------------
c     subprogram 10. bicube_write_yx.
c     produces ascii and binary output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_yx(bcs,out,bin,iua,iub,interp)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      INTEGER :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      CHARACTER(80) :: format1,format2
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"ix = ",i3,", ",a6," = ",1p,e11.3)
 20   FORMAT(1x,"ix = ",i3,", jx = ",i1,", ",a6," = ",1p,e11.3)
 30   FORMAT('(/4x,"iy",4x,a6,1x,',i3.3,'(4x,a6,1x)/)')
 40   FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(format1,30)bcs%nqty
         WRITE(format2,40)bcs%nqty
         WRITE(iua,'(1x,a/)')"input data"
      ENDIF
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      DO ix=0,bcs%mx
         x=bcs%xs(ix)
         IF(out)then
            WRITE(iua,10)ix,bcs%xtitle,bcs%xs(ix)
            WRITE(iua,format1)bcs%ytitle,
     $           (bcs%title(iqty),iqty=1,bcs%nqty)
         ENDIF
         DO iy=0,bcs%my
            y=bcs%ys(iy)
            CALL bicube_eval(bcs,x,y,0)
            IF(out)WRITE(iua,format2)iy,y,bcs%f
            IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4)
         ENDDO
         IF(out)WRITE(iua,format1)bcs%ytitle,
     $        (bcs%title(iqty),iqty=1,bcs%nqty)
         IF(bin)WRITE(iub)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
      IF(interp)THEN
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data"
         DO ix=0,bcs%mx-1
            dx=(bcs%xs(ix+1)-bcs%xs(ix))/4
            DO jx=0,4
               x=bcs%xs(ix)+dx*jx
               IF(out)then
                  WRITE(iua,20)ix,jx,bcs%xtitle,x
                  WRITE(iua,format1)bcs%ytitle,
     $                 (bcs%title(iqty),iqty=1,bcs%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
               DO iy=0,bcs%my-1
                  dy=(bcs%ys(iy+1)-bcs%ys(iy))/4
                  DO jy=0,4
                     y=bcs%ys(iy)+dy*jy
                     CALL bicube_eval(bcs,x,y,0)
                     IF(out)WRITE(iua,format2)iy,y,bcs%f
                     IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4)
                  ENDDO
                  IF(out)WRITE(iua,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over x.
c-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)bcs%ytitle,
     $                 (bcs%title(iqty),iqty=1,bcs%nqty)
               IF(bin)WRITE(iub)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_yx
c-----------------------------------------------------------------------
c     subprogram 11. bicube_write_arrays.
c     produces ascii and binary output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_arrays(bcs,out,iua,iqty)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out
      INTEGER, INTENT(IN) :: iua,iqty

      CHARACTER(80) :: format1,format2
      INTEGER :: ix,iy,my
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/2x,"ix/iy",',i3.3,'(3x,i3.3,5x)/)')
 20   FORMAT('(i5,1p,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. out)RETURN
      my=MIN(bcs%my,32)
      WRITE(iua,'(a,i2/)')"iqty = ",iqty
c-----------------------------------------------------------------------
c     write fs.
c-----------------------------------------------------------------------
      WRITE(format1,10)my+1
      WRITE(format2,20)my+1
      WRITE(iua,"(a)")"fs:"
      WRITE(iua,format1)(iy,iy=0,my)
      WRITE(iua,format2)(ix,(bcs%fs(ix,iy,iqty),iy=0,my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,my)
c-----------------------------------------------------------------------
c     write fsx.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsx:"
      WRITE(iua,format1)(iy,iy=0,my)
      WRITE(iua,format2)(ix,(bcs%fsx(ix,iy,iqty),iy=0,my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,my)
c-----------------------------------------------------------------------
c     write fsy.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsy:"
      WRITE(iua,format1)(iy,iy=0,my)
      WRITE(iua,format2)(ix,(bcs%fsy(ix,iy,iqty),iy=0,my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,my)
c-----------------------------------------------------------------------
c     write fsxy.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsxy:"
      WRITE(iua,format1)(iy,iy=0,my)
      WRITE(iua,format2)(ix,(bcs%fsxy(ix,iy,iqty),iy=0,my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,my)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_arrays
c-----------------------------------------------------------------------
c     subprogram 12. bicube_copy.
c     copies one bicube type to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_copy(bcs1,bcs2)

      TYPE(bicube_type), INTENT(IN) :: bcs1
      TYPE(bicube_type), INTENT(INOUT) :: bcs2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ALLOCATED(bcs2%xs))CALL bicube_dealloc(bcs2)
      CALL bicube_alloc(bcs2,bcs1%mx,bcs1%my,bcs1%nqty)
      bcs2%xs=bcs1%xs
      bcs2%ys=bcs1%ys
      bcs2%fs=bcs1%fs
      bcs2%fsx=bcs1%fsx
      bcs2%fsy=bcs1%fsy
      bcs2%fsxy=bcs1%fsxy
      bcs2%name=bcs1%name
      bcs2%xtitle=bcs1%xtitle
      bcs2%ytitle=bcs1%ytitle
      bcs2%title=bcs1%title
      bcs2%periodic=bcs1%periodic
      bcs2%xpower=bcs1%xpower
      bcs2%ypower=bcs1%ypower
      bcs2%x0=bcs1%x0
      bcs2%y0=bcs1%y0
      bcs2%xext=bcs1%xext
      bcs2%yext=bcs1%yext
      bcs2%fext=bcs1%fext
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_copy
c-----------------------------------------------------------------------
c     subprogram 13. bicube_extrema.
c     finds extrema.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_extrema(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      CHARACTER(80) :: message
      INTEGER :: iext,iqty,it
      INTEGER, PARAMETER :: itmax=20
      INTEGER, DIMENSION(2,2) :: jext
      REAL(r8), PARAMETER :: eps=1e-10
      REAL(r8) :: x,y,dx,dy,lx,ly,lf,f,df,adet
      REAL(r8), DIMENSION(2,2) :: amat,ainv
c-----------------------------------------------------------------------
c     compute lengths.
c-----------------------------------------------------------------------
      lx=bcs%xs(bcs%mx)-bcs%xs(0)
      ly=bcs%ys(bcs%my)-bcs%ys(0)
c-----------------------------------------------------------------------
c     start loops over iqty.
c-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty
         jext(:,1)=MINLOC(bcs%fs(:,:,iqty))-1
         jext(:,2)=MAXLOC(bcs%fs(:,:,iqty))-1
         bcs%fext(1,iqty)=bcs%fs(jext(1,1),jext(2,1),iqty)
         bcs%fext(2,iqty)=bcs%fs(jext(1,2),jext(2,2),iqty)
         lf=bcs%fext(2,iqty)-bcs%fext(1,iqty)
c-----------------------------------------------------------------------
c     start loops over extrema.
c-----------------------------------------------------------------------
         DO iext=1,2
            x=bcs%xs(jext(1,iext))
            y=bcs%ys(jext(2,iext))
            f=HUGE(f)
            dx=lx
            dy=ly
            it=0
c-----------------------------------------------------------------------
c     locate extema by newton iteration.
c-----------------------------------------------------------------------
            DO
               CALL bicube_eval(bcs,x,y,2)
               df=bcs%f(iqty)-f
               IF(ABS(dx) < eps*lx .OR. ABS(dy) < eps*ly
     $              .OR. ABS(df) < eps*lf .OR. it >= itmax)EXIT
               it=it+1
               f=bcs%f(iqty)
               amat(1,1)=bcs%fxx(iqty)
               amat(2,2)=bcs%fyy(iqty)
               amat(1,2)=bcs%fxy(iqty)
               amat(2,1)=bcs%fxy(iqty)
               adet=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
               ainv(1,1)=amat(2,2)
               ainv(2,2)=amat(1,1)
               ainv(1,2)=-amat(1,2)
               ainv(2,1)=-amat(2,1)
               ainv=ainv/adet
               dx=-ainv(1,1)*bcs%fx(iqty)-ainv(1,2)*bcs%fy(iqty)
               dy=-ainv(2,1)*bcs%fx(iqty)-ainv(2,2)*bcs%fy(iqty)
               x=x+dx
               y=y+dy
            ENDDO
c-----------------------------------------------------------------------
c     abort on failure.
c-----------------------------------------------------------------------
            IF(it >= itmax)THEN
               WRITE(message,'(a,i3,a)')
     $              "bicube_extrema: convergence failure for iqty = ",
     $              iqty,"."
               CALL program_stop(message)
            ENDIF
c-----------------------------------------------------------------------
c     finish loops over iext and iqty.
c-----------------------------------------------------------------------
            bcs%xext(iext,iqty)=x
            bcs%yext(iext,iqty)=y
            bcs%fext(iext,iqty)=bcs%f(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_extrema
      END MODULE bicube_mod
