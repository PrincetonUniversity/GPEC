c-----------------------------------------------------------------------
c     file debug2.f.
c     prints large us.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. debug_mod.
c     1. debug1.
c     2. debug2.
c     3. debug3.
c-----------------------------------------------------------------------
c     subprogram 0. debug_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE debug_mod
      USE local_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. debug1.
c     prints 1d arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE debug1(u,name,unit)

      REAL(r8), DIMENSION(:,0:), INTENT(IN) :: u
      CHARACTER(*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: unit

      CHARACTER(80) :: format1,format2
      INTEGER :: nqty,nx,iqty,i
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"i",',i2,'(6x,i2,4x)/)')
 20   FORMAT('(i4,1p,',i2,'e11.3)')
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nqty=SIZE(u,1)
      nx=SIZE(u,2)
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      WRITE(format1,10)nqty
      WRITE(format2,20)nqty
c-----------------------------------------------------------------------
c     write array.
c-----------------------------------------------------------------------
      WRITE(unit,'(a)')TRIM(name)//":"
      WRITE(unit,format1)(iqty,iqty=1,nqty)
      DO i=0,nx-1
         WRITE(unit,format2)i,u(:,i)
      ENDDO
      WRITE(unit,format1)(iqty,iqty=1,nqty)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE debug1
c-----------------------------------------------------------------------
c     subprogram 2. debug2.
c     prints 2d arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE debug2(u,name,unit)

      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: u
      CHARACTER(*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: unit

      CHARACTER(80) :: format1,format2
      INTEGER :: nqty,nx,ny,iqty,ix,iy
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/2x,"iy",2x,"ix",',i2,'(5x,i2,4x)/)')
 20   FORMAT('(2i4,1p,',i2,'e11.3)')
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nqty=SIZE(u,1)
      nx=SIZE(u,2)
      ny=SIZE(u,3)
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      WRITE(format1,10)nqty
      WRITE(format2,20)nqty
c-----------------------------------------------------------------------
c     write array.
c-----------------------------------------------------------------------
      WRITE(unit,'(a)')TRIM(name)//":"
      WRITE(unit,format1)(iqty,iqty=1,nqty)
      DO iy=0,ny-1
         DO ix=0,nx-1
            WRITE(unit,format2)iy,ix,u(:,ix,iy)
         ENDDO
         WRITE(unit,format1)(iqty,iqty=1,nqty)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE debug2
c-----------------------------------------------------------------------
c     subprogram 3. debug3.
c     prints 2d arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE debug3(u,name,unit)

      REAL(r8), DIMENSION(:,:,0:,0:), INTENT(IN) :: u
      CHARACTER(*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: unit

      CHARACTER(80) :: format1,format2
      INTEGER :: nqty1,nqty2,nx,ny,iqty,jqty,ix,iy
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/2y,"iy",2x,"ix",',i3,'(5x,i2,",",i2,3x)/)')
 20   FORMAT('(2i4,1p,',i3,'e11.3)')
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nqty1=SIZE(u,1)
      nqty2=SIZE(u,2)
      nx=SIZE(u,3)
      ny=SIZE(u,4)
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      WRITE(format1,10)nqty1*nqty2
      WRITE(format2,20)nqty1*nqty2
c-----------------------------------------------------------------------
c     write array.
c-----------------------------------------------------------------------
      WRITE(unit,'(a)')TRIM(name)//":"
      WRITE(unit,format1)((iqty,jqty,iqty=1,nqty1),jqty=1,nqty2)
      DO iy=0,ny-1
         DO ix=0,nx-1
            WRITE(unit,format2)iy,ix,u(:,:,ix,iy)
         ENDDO
         WRITE(unit,format1)((iqty,jqty,iqty=1,nqty1),jqty=1,nqty2)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE debug3
      END MODULE debug_mod
