c-----------------------------------------------------------------------
c     file match_mod.f.
c     global definitions for match code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. match_mod.
c     1. match_dealloc.
c-----------------------------------------------------------------------
c     subprogram 0. match_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE match_mod
      USE deltar_mod
      USE root_mod
      USE utils_mod
      IMPLICIT NONE

      TYPE :: solution_type
      INTEGER :: msol
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u
      END TYPE solution_type

      TYPE :: fixfac_type
      INTEGER :: msol
      INTEGER, DIMENSION(:), POINTER :: index
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fixfac,transform,gauss
      END TYPE fixfac_type

      TYPE :: sing_type
      INTEGER :: msol_l,msol_r,jfix,jpert
      REAL(r8) :: psifac,q,q1
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: ca_l,ca_r
      TYPE(resist_type) :: restype
      END TYPE sing_type

      LOGICAL :: ideal_flag=.FALSE.,res_flag=.FALSE.
      LOGICAL, DIMENSION(:), POINTER :: sing_flag
      INTEGER :: mfix,mhigh,mlow,mpert,mstep,nn,msing,mmatch,mbit,mterm
      INTEGER, DIMENSION(:), POINTER :: fixstep
      REAL(r8) :: sfac0=1
      REAL(r8), DIMENSION(:), POINTER :: psifac,rho,q,et
      COMPLEX(r8), DIMENSION(:), POINTER :: deltap1,deltap2,ff
      COMPLEX(r8), DIMENSION(:,:), POINTER :: wt,match,bmatch
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: v

      TYPE(solution_type), DIMENSION(:), POINTER :: soltype
      TYPE(fixfac_type), DIMENSION(:), POINTER :: fixtype
      TYPE(sing_type), DIMENSION(:), POINTER :: singtype

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. match_dealloc.
c     deallocate storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_dealloc

      INTEGER :: istep,ifix,ising
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DO istep=0,mstep
         DEALLOCATE(soltype(istep)%u)
      ENDDO
      IF(ideal_flag)THEN
         DO ifix=1,mfix
            DEALLOCATE(fixtype(ifix)%fixfac,fixtype(ifix)%index)
         ENDDO
      ENDIF
      IF(res_flag)THEN
         DEALLOCATE(match,ff,deltap1,deltap2)
      ENDIF
      DO ising=1,msing
         DEALLOCATE(singtype(ising)%ca_l,singtype(ising)%ca_r)
      ENDDO
      DEALLOCATE(psifac,rho,q,soltype,v,singtype,fixstep,fixtype,
     $     sing_flag,et,wt)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_dealloc
      END MODULE match_mod
