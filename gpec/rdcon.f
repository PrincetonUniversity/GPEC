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
c     reads globalsol.bin
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE rdcon_read_solution
      INTEGER isol,ip,ipert,ifix
c-----------------------------------------------------------------------
c     Only read if flag is on
c-----------------------------------------------------------------------
      IF (.NOT.galsol%gal_flag) RETURN
c-----------------------------------------------------------------------
c     read solutions for coupling to the coil (final rpec run).
c-----------------------------------------------------------------------
      IF(verbose) PRINT *,"Reading RDCON solutions"
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
c     force global variables
c-----------------------------------------------------------------------
      ! deallocate soltype so it forces an error if we ever try to use it
      IF(ALLOCATED(soltype))THEN
         DO ip=0,mstep
            DEALLOCATE(soltype(ip)%u)
         ENDDO
         DEALLOCATE(soltype)
      ENDIF
      IF(ALLOCATED(fixtype))THEN
         DO ifix=0,mfix
            IF(ALLOCATED(fixtype(ifix)%fixfac))
     $         DEALLOCATE(fixtype(ifix)%fixfac)
            IF(ALLOCATED(fixtype(ifix)%index))
     $         DEALLOCATE(fixtype(ifix)%index)
            IF(ALLOCATED(fixtype(ifix)%transform))
     $         DEALLOCATE(fixtype(ifix)%transform)
            IF(ALLOCATED(fixtype(ifix)%gauss))
     $         DEALLOCATE(fixtype(ifix)%gauss)
         ENDDO
         DEALLOCATE(fixtype)
      ENDIF
      ! force coil matching to method 1 since we don't have energies (method 0)
      IF(verbose)THEN
        PRINT *,"  > Forcing resp_index to 1 for RDCON interfacing"
      ENDIF
      resp_index=1
      ! force new mstep and all the corresponding radial grids
      IF(verbose)THEN
        PRINT '(1x,a38,a7,I4,a4,I4)',
     $    "  > Redefining radial grid with RDCON.",
     $    " mstep ",mstep," -> ",galsol%tot_grids
      ENDIF
      mstep = galsol%tot_grids
      DEALLOCATE(psifac,rhofac,qfac)
      ALLOCATE(psifac(0:mstep),rhofac(0:mstep),qfac(0:mstep))
      psifac = galsol%psifac
      rhofac = SQRT(psifac)
      DO ip=0,mstep
         CALL spline_eval(sq,psifac(ip),0)
         qfac(ip) = sq%f(4)
      ENDDO
      ! m's must already match
      IF(galsol%mpert/=mpert) THEN
         PRINT *,'Galerkin mpert = ',galsol%mpert
         PRINT *,'Ideal    mpert = ',mpert
         stop "ERROR: Galerkin mpert not ideal mpert."
      ENDIF
      ! These are the solution splines that can actually be used by idcon_build
      CALL cspline_dealloc(u1)
      CALL cspline_dealloc(u2)
      CALL cspline_dealloc(u3)
      CALL cspline_dealloc(u4)
      CALL cspline_alloc(u1,mstep,mpert)
      CALL cspline_alloc(u2,mstep,mpert)
      CALL cspline_alloc(u3,mstep,mpert)
      CALL cspline_alloc(u4,mstep,mpert)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rdcon_read_solution
      END MODULE rdcon_mod
