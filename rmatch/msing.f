c-----------------------------------------------------------------------
c     file sing.f.
c     computations relating to singular surfaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. msing_mod.
c     1. msing_read.
c     2. msing_get_ua.
c     3. msing_write.
c     4. msing_run.
c     5. msing_estimate_zo.
c     6. msing_init.
c-----------------------------------------------------------------------
c     subprogram 0. msing_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE msing_mod
      USE local_mod
      IMPLICIT NONE

      TYPE :: sing_type
      INTEGER :: order
      INTEGER, DIMENSION(1) :: r1
      INTEGER, DIMENSION(2) :: r2
      REAL(r8) :: psifac,q
      COMPLEX(r8) :: alpha
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmatl,vmatr
      END TYPE sing_type
      
      LOGICAL :: msing_diagnose=.FALSE.
      LOGICAL :: us_flag=.TRUE.,ul_flag=.TRUE.
      CHARACTER(80) :: vmat_filename
      INTEGER, PRIVATE :: mpert,msing
      TYPE(sing_type), DIMENSION(:), POINTER, PRIVATE :: sing
      TYPE(sing_type), POINTER, PRIVATE :: singp

      INTEGER :: sing_nout=100
      REAL(r8) :: sing_frac=.1

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. msing_read.
c     reads vmat.bin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE msing_read(vmat_filename)
      CHARACTER(*) :: vmat_filename
      INTEGER :: ising
c-----------------------------------------------------------------------
c     read power series data to binary file.
c-----------------------------------------------------------------------
      OPEN(UNIT=debug_unit,FILE=TRIM(vmat_filename),STATUS="OLD",
     $     FORM="UNFORMATTED")
      READ(debug_unit)mpert,msing
      ALLOCATE(sing(msing))
      DO ising=1,msing
         singp => sing(ising)
         READ(debug_unit)singp%order,singp%alpha,singp%psifac,
     $        singp%q,singp%r1,singp%r2
         ALLOCATE(singp%vmatl(mpert,2*mpert,2,0:2*singp%order))
         ALLOCATE(singp%vmatr(mpert,2*mpert,2,0:2*singp%order))
         READ(debug_unit)singp%vmatl,singp%vmatr
      ENDDO
      CLOSE(UNIT=debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE msing_read
c-----------------------------------------------------------------------
c     subprogram 2. msing_get_ua.
c     computes asymptotic series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE msing_get_ua(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ua

      INTEGER :: iorder
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: dpsi,sqrtfac,pfac
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
c-----------------------------------------------------------------------
c     set pointers and compute distance from singular surface.
c-----------------------------------------------------------------------
      singp => sing(ising)
      IF(psifac < singp%psifac)THEN
         vmat => singp%vmatl
         dpsi=singp%psifac-psifac
      ELSE
         vmat => singp%vmatr
         dpsi=psifac-singp%psifac         
      ENDIF
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      sqrtfac=SQRT(dpsi)
      pfac=dpsi**singp%alpha
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      ua=vmat(:,:,:,2*singp%order)
      DO iorder=2*singp%order-1,0,-1
         ua=ua*sqrtfac+vmat(:,:,:,iorder)
      ENDDO
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      ua(r1,:,1)=ua(r1,:,1)/sqrtfac
      ua(r1,:,2)=ua(r1,:,2)*sqrtfac
      ua(:,r2(1),:)=ua(:,r2(1),:)/pfac
      ua(:,r2(2),:)=ua(:,r2(2),:)*pfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE msing_get_ua
c-----------------------------------------------------------------------
c     subprogram 3. msing_write.
c     computes and writes resonant power-series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE msing_write(ising,cofout,delta,unit)

      INTEGER, INTENT(IN) :: ising,unit
      COMPLEX(r8), INTENT(IN), DIMENSION(:) :: cofout
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: delta

      INTEGER :: ipsi,jsing,ksing
      REAL(r8) :: psifac,psil,psir,dpsil,dpsir,dpsi
      COMPLEX(r8) :: sol,ul,us
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
c-----------------------------------------------------------------------
c     compute range of psi values on left of singular surface.
c----------------------------------------------------------------------- 
      singp => sing(ising)
      IF(ising == 1)THEN
         psil=(1-sing_frac)*singp%psifac
      ELSE
         psil=sing(ising-1)%psifac
     $        +(1-sing_frac)*(singp%psifac-sing(ising-1)%psifac)
      ENDIF
      dpsil=(singp%psifac-psil)/sing_nout
c-----------------------------------------------------------------------
c     compute range of psi values on right of singular surface.
c----------------------------------------------------------------------- 
      IF(ising == msing)THEN
         psir=singp%psifac+sing_frac*(1-singp%psifac)
      ELSE
         psir=singp%psifac+sing_frac*(sing(ising+1)%psifac-singp%psifac)
      ENDIF
      dpsir=(psir-singp%psifac)/sing_nout
c-----------------------------------------------------------------------
c     compute and write solutions.
c-----------------------------------------------------------------------
      psifac=psil
      dpsi=dpsil
      ksing=2*ising-1
      DO ipsi=-sing_nout,sing_nout
         IF (ipsi==0) THEN
            psifac=singp%psifac+dpsir
            dpsi=dpsir
            ksing=2*ising
            WRITE(unit)
            CYCLE
         ENDIF
         sol=0
         CALL msing_get_ua(ising,psifac,ua)
         ul=ua(singp%r1(1),singp%r2(1),1)
         us=ua(singp%r1(1),singp%r2(2),1)
         IF (ul_flag) THEN
            sol=sol+cofout(ksing)*ul
         ENDIF
         IF (us_flag) THEN
            DO jsing=1,2*msing 
               sol=sol+cofout(jsing)*delta(jsing,ksing)*us
            ENDDO
         ENDIF
         WRITE(unit)REAL(psifac,4),REAL(sol,4),
     $             REAL(IMAG(sol),4),floored_log(sol)
         psifac=psifac+dpsi
      ENDDO
      WRITE(unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE msing_write
c-----------------------------------------------------------------------
c     subprogram 4. msing_run.
c     run msing diagnostic.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE msing_run(cofout,delta)
      COMPLEX(r8), INTENT(IN), DIMENSION(:) :: cofout
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: delta
      
      CHARACTER(100):: filename
      INTEGER :: ising
      IF (.NOT.msing_diagnose) RETURN
c-----------------------------------------------------------------------
c     do msing diagnostic.
c-----------------------------------------------------------------------     
      DO ising=1,msing
         WRITE (filename,"(I4)")ising
         filename=ADJUSTL(filename)
         WRITE(filename,*) 'msing_'
     $                       //TRIM(filename)//'.bin'

         CALL bin_open(bin_unit,filename,"REPLACE",
     $                 "REWIND","none")
         CALL msing_write(ising,cofout,delta,bin_unit)
         CALL bin_close(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE msing_run
c-----------------------------------------------------------------------
c     subprogram 5. msing_estimate_zo.
c     estimate zo value for zeroth order dominant 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE msing_estimate_zo(zo,ising)
      REAL(r8), INTENT(OUT) :: zo
      INTEGER, INTENT(IN) :: ising 
      REAL(r8) :: tmp1,tmp2,tmp3,tmp4
      singp => sing(ising)
      tmp1=MAXVAL(ABS(singp%vmatr(:,singp%r2(1),1,0)))
      tmp2=MAXVAL(ABS(singp%vmatr(:,singp%r2(1),1,2)))
      tmp3=tmp1/tmp2
      tmp1=MAXVAL(ABS(singp%vmatr(:,singp%r2(2),1,0)))
      tmp2=MAXVAL(ABS(singp%vmatr(:,singp%r2(2),1,2)))
      tmp4=tmp1/tmp2
      IF (tmp4 < tmp3) tmp3=tmp4
      zo=tmp3
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------      
      END SUBROUTINE msing_estimate_zo
c-----------------------------------------------------------------------
c     subprogram 6. msing_init.
c     initialize msing module.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE msing_init
      NAMELIST/msing_list/ msing_diagnose,vmat_filename,sing_frac,
     $                     sing_nout,us_flag,ul_flag
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE='rmatch.in',STATUS="OLD")
      READ(in_unit,NML=msing_list)
      CLOSE(UNIT=in_unit)
      CALL msing_read(vmat_filename)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE msing_init
      END MODULE msing_mod
