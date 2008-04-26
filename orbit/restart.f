c-----------------------------------------------------------------------
c     file restart.f.
c     reads and writes data for restarting run.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. restart_mod.
c     1. restart_write.
c     2. restart_read.
c-----------------------------------------------------------------------
c     subprogram 0. restart_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE restart_mod
      USE orbit_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. restart_write.
c     writes restart data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE restart_write(istep,t,y)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(IN) :: y
c-----------------------------------------------------------------------
c     write restart data.
c-----------------------------------------------------------------------
      CALL ascii_open(field_out_unit,"restart.out","UNKNOWN")
      WRITE(field_out_unit,'(a9,2i4)')particle,ai,zi
      WRITE(field_out_unit,'(a9,i7,1p,e24.16)')orbit_type,istep,t
      WRITE(field_out_unit,'(1p,3e24.16)')energy0,mu_0,tau_report
      WRITE(field_out_unit,'(1p,3e24.16)')y
      CALL ascii_close(field_out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE restart_write
c-----------------------------------------------------------------------
c     subprogram 2. restart_read.
c     reads restart data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE restart_read(istep,t,y)

      INTEGER, INTENT(OUT) :: istep
      REAL(r8), INTENT(OUT) :: t
      REAL(r8), DIMENSION(3,2), INTENT(OUT) :: y

      INTEGER :: ios
c-----------------------------------------------------------------------
c     read restart data.
c-----------------------------------------------------------------------
      OPEN(UNIT=field_out_unit,FILE="restart.out",STATUS="OLD",
     $     IOSTAT=ios)
      IF(ios == 0)THEN
         READ(field_out_unit,'(a9,2i4)')particle,ai,zi
         READ(field_out_unit,'(a9,i7,1p,e24.16)')orbit_type,istep,t
         READ(field_out_unit,'(1p,3e24.16)')energy0,mu_0,tau_report
         READ(field_out_unit,'(1p,3e24.16)')y
         CLOSE(UNIT=field_out_unit)
      ELSE
         istep=0
         t=0
         restart_flag=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE restart_read
      END MODULE restart_mod
