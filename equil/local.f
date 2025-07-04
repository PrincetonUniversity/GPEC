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
c     7. floored_log.
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

      CHARACTER(6) :: convert_type="none"
      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.false.
      INTEGER, PARAMETER ::
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi,mu0=4e-7_r8*pi,
     $     rtod=180/pi,dtor=pi/180,alog10=2.3025850929940459_r8
      REAL(r8), PARAMETER :: e=1.6021917e-19,
     $     mp=1.672614e-27,me=9.1091e-31,
     $     c=2.997925e10,jzero=2.4048255577_r8,ev=e
      REAL(r8), PARAMETER :: zero=0,one=1,two=2,half=.5
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      COMPLEX(r8), PARAMETER :: one_c=(1,0)
      COMPLEX(r8), PARAMETER :: negone_c=(-1,0)
      COMPLEX(r8), PARAMETER :: zero_c=(0,0)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(mode,unit,op_cpuseconds,op_wallseconds)
      
      INTEGER, INTENT(IN) :: mode,unit
      REAL(r4), INTENT(OUT), OPTIONAL :: op_cpuseconds,op_wallseconds

      INTEGER(8), SAVE :: count_rate, wall_start
      REAL(r4), SAVE :: start
      REAL(r4) :: seconds
      INTEGER(8) :: hrs,mins,secs, wall_seconds, count_max
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL CPU_TIME(start)
         CALL SYSTEM_CLOCK(wall_start)
      ELSE
         ! report cpu time
         CALL CPU_time(seconds)
         seconds=seconds-start
         secs = int(seconds)
         hrs = secs/(60*60)
         mins = (secs-hrs*60*60)/60
         secs = secs-hrs*60*60-mins*60
         IF(PRESENT(op_cpuseconds))THEN
            ! simply provide the time to the caller
            op_cpuseconds = seconds
         ELSE
            ! write the time to terminal and file
            IF(hrs>0)THEN
               WRITE(*,'(1x,a,i3,a,i2,a,i2,a)') "Total cpu time = ",
     $            hrs," hours, ",mins," minutes, ",secs," seconds"
            ELSEIF(mins>0)THEN
               WRITE(*,'(1x,a,i2,a,i2,a)') "Total cpu time = ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(secs>0)THEN
               WRITE(*,'(1x,a,i2,a)') "Total cpu time = ",secs,
     $            " seconds"
            ENDIF
            WRITE(unit,10) "Total cpu time = ",seconds," seconds"
         ENDIF
         ! report wall time
         CALL SYSTEM_CLOCK(wall_seconds, count_rate, count_max)
         seconds=(wall_seconds-wall_start)/REAL(count_rate, 8)
         secs = int(seconds)
         hrs = secs/(60*60)
         mins = (secs-hrs*60*60)/60
         secs = secs-hrs*60*60-mins*60
         IF(PRESENT(op_wallseconds))THEN
            op_wallseconds = seconds
         ELSE
            IF(hrs>0)THEN
               WRITE(*,'(1x,a,i3,a,i2,a,i2,a)') "Total wall time = ",
     $            hrs," hours, ",mins," minutes, ",secs," seconds"
            ELSEIF(mins>0)THEN
               WRITE(*,'(1x,a,i2,a,i2,a)') "Total wall time = ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(secs>0)THEN
               WRITE(*,'(1x,a,i2,a)') "Total wall time = ",secs,
     $            " seconds"
            ENDIF
            WRITE(unit,10) "Total wall time = ",seconds," seconds"
         ENDIF
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
      SUBROUTINE bin_open(unit,name,stat,pos,convert_type)

      CHARACTER(*), INTENT(IN) :: name,stat,pos,convert_type
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      SELECT CASE(convert_type)
      CASE("none")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED")
      CASE("big")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED",CONVERT="BIG_ENDIAN")
      CASE("little")
         OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $        FORM="UNFORMATTED",CONVERT="LITTLE_ENDIAN")
      CASE DEFAULT
         CALL program_stop
     $        ("Cannot recognize convert_type = "//TRIM(convert_type))
      END SELECT
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
      WRITE(*,'(1x,2a)') 'PROGRAM STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop
c-----------------------------------------------------------------------
c     subprogram 7. floored_log.
c     returns bounded log of specific component.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION floored_log(u) RESULT(ulog)

      COMPLEX(r8), INTENT(IN) :: u
      REAL(r4) :: ulog

      REAL, PARAMETER :: minlog=-15
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(u == 0)THEN
         ulog=minlog
      ELSE
         ulog=LOG10(ABS(u))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION floored_log
      END MODULE local_mod
