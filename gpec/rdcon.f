c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     read resistive DCON data and control jump conditions
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. rdcon_mod.
c     1. rdcon_read.
c     2. rdcon_read_solution.
c-----------------------------------------------------------------------
c     subprogram 0. rdcon_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rdcon_mod
      USE gpglobal_mod
      USE idcon_mod
      USE ismath_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. rdcon_read.
c     reads gal_solution.bin
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE rdcon_read

      CHARACTER(128) :: message
      INTEGER :: istep,ising,ipert
c-----------------------------------------------------------------------
c     allocate arrays and prepare to read data.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,rdconfile,"OLD","REWIND","none")
      READ(in_unit)rsing,rnqty,rnx
      IF (rsing /= msing) THEN
         WRITE(message,'(a)')"GPEC needs the same msing number"
         CALL gpec_stop(message)
      ENDIF
      ALLOCATE(rsoltype(rsing),rpsifac(0:rnx)) 
      rsoltype(:)%msol=rsing
      READ(in_unit)rpsifac
      DO ising=1,rsing
          ALLOCATE(rsoltype(ising)%u(rnqty,0:rnx,2))
          READ(in_unit)rsoltype(ising)%u(:,:,1)
          READ(in_unit)rsoltype(ising)%u(:,:,2)
      ENDDO
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF (bin_flag) THEN
          CALL bin_open(bin_unit,
     $           "rsol1.bin","UNKNOWN","REWIND","none")
          DO ipert=1,rnqty
            DO istep=0,rnx
               IF (.NOT.(ABS(rsoltype(3)%u(ipert,istep,1)) .LT. 1e4)) 
     $            cycle
               IF (.NOT.(ABS(rsoltype(3)%u(ipert,istep,2)) .LT. 1e4)) 
     $            cycle
               CALL spline_eval(sq,rpsifac(istep),0)
               singfac(ipert)=mfac(ipert)-nn*sq%f(4)
               WRITE(bin_unit)REAL(rpsifac(istep),4),
     $              REAL((REAL(rsoltype(3)%u(ipert,istep,1))-
     $              REAL(rsoltype(3)%u(ipert,istep,2)))*
     $              singfac(ipert),4),
     $              REAL((AIMAG(rsoltype(3)%u(ipert,istep,1))-
     $              AIMAG(rsoltype(3)%u(ipert,istep,2)))*
     $              singfac(ipert),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rdcon_read
c-----------------------------------------------------------------------
c     subprogram 2. rdcon_read_solution.
c     reads gal_solution.bin
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE rdcon_read_solution
      INTEGER isol,ip,ipert
c-----------------------------------------------------------------------
c     force the coil matching with option 1.
c-----------------------------------------------------------------------
      IF (.NOT.galsol%gal_flag) RETURN
      resp_index=1
c-----------------------------------------------------------------------
c     read solutions for coupling to the coil (final rpec run).
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,rdconfile,"OLD","REWIND","none")
      READ(bin_unit) galsol%mpert,galsol%tot_grids,galsol%mtot,
     $                 galsol%mlow,galsol%mhigh
      ALLOCATE (galsol%psifac(0:galsol%tot_grids),
     $          galsol%u(galsol%mpert,0:galsol%tot_grids,galsol%mtot))
      READ(bin_unit) galsol%psifac
      DO isol=1,galsol%mtot
         DO ip=0,galsol%tot_grids
            DO ipert=1,galsol%mpert
               READ(bin_unit) galsol%u(ipert,ip,isol)
            ENDDO
         ENDDO
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rdcon_read_solution
      END MODULE rdcon_mod
