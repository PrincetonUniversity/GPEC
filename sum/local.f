c-----------------------------------------------------------------------
c     file local.f
c     local defintions for most computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. local.
c     1. timer.
c     2. bin_open.
c     3. bin_close.
c     4. ascii_open.
c     5. ascii_close.
c     6. program_stop.
c-----------------------------------------------------------------------
c     subprogram 0. local.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE local_mod
      USE io_mod
      IMPLICIT NONE

      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.false.
      INTEGER, PARAMETER ::
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi,mu0=4e-7_r8*pi,
     $     rtod=180/pi,dtor=pi/180,alog10=2.302585093_r8
      REAL(r8), PARAMETER :: e=1.6021917e-19,
     $     mp=1.672614e-27,me=9.1091e-31,
     $     c=2.997925e10,jzero=2.4048255577_r8,ev=e
      REAL(r8), PARAMETER :: zero=0,one=1
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(mode,unit)
      
      INTEGER, INTENT(IN) :: mode,unit

      CHARACTER(10) :: date,time,zone
      INTEGER, DIMENSION(8) :: values
      REAL(4), SAVE :: seconds
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
      ELSE
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
     $        -seconds
         WRITE(unit,10)"Total cpu time = ",seconds," seconds"
         WRITE(*,10)"Total cpu time = ",seconds," seconds"//CHAR(7)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     subprogram 2. bin_open.
c     opens a binary input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bin_open(unit,name,stat,pos)

      CHARACTER(*), INTENT(IN) :: name,stat,pos
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $     FORM="UNFORMATTED")  !,CONVERT="BIG_ENDIAN")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bin_open
c-----------------------------------------------------------------------
c     subprogram 3. bin_close.
c     close a binary input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bin_close(unit)

      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      CLOSE(UNIT=unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bin_close
c-----------------------------------------------------------------------
c     subprogram 4. ascii_open.
c     opens a ascii input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ascii_open(unit,name,stat)

      CHARACTER(*), INTENT(IN) :: name,stat
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=unit,FILE=name,STATUS=stat)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ascii_open
c-----------------------------------------------------------------------
c     subprogram 5. ascii_close.
c     close a ascii input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ascii_close(unit)

      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     work.
c-----------------------------------------------------------------------
      CLOSE(UNIT=unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ascii_close
c-----------------------------------------------------------------------
c     subprogram 6. program_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE program_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      CALL timer(1,out_unit)
      CALL ascii_close(out_unit)
      WRITE(*,'(1x,2a)') 'PROGRAM_STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop
      END MODULE local_mod
