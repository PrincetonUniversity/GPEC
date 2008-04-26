c-----------------------------------------------------------------------
c     program orbit.
c     solves for the motion of a charged particle in a combined trap.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. orbit_main.
c-----------------------------------------------------------------------
c     subprogram 1. orbit_main.
c     controls program flow and i/o.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM orbit_main
      USE inner_mod
      USE outer_mod
      USE inner_gc_mod
      USE outer_gc_mod
      IMPLICIT NONE

      LOGICAL :: done,cross=.FALSE.
      INTEGER :: istep
      REAL(r8) :: t,tau,errstep
      REAL(r8), DIMENSION(3,2) :: y
      REAL(r8), DIMENSION(4) :: y_gc

      NAMELIST/orbit_input/orbit_type,particle,ai,zi,psi0,theta0,zeta0,
     $     energy0,alpha0,phase0,taumax,report,stride,tol,nstep,
     $     orbit_out,orbit_bin,break_flag,errmax,out_flag,r0,z0,phi0,
     $     cross,field_out,field_bin,interp_flag,restart_flag,ode_bin,
     $     gc_flag
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(//1x,"particle",2x,"ai",2x,"zi",2x,"orbit",3x,"cross"
     $     //1x,a8,2i4,a7,l6/)
 20   FORMAT(/5x,"psi0",6x,"theta0",5x,"zeta0",8x,"r0",9x,"z0",8x,"phi0"
     $     //1p,6e11.3/)
 30   FORMAT(/3x,"energy0",5x,"alpha0",5x,"phase0"//1p,3e11.3/)
 40   FORMAT(/4x,"taumax",5x,"report",3x,"stride",3x,"tol",7x,"errmax"
     $     //1p,2e11.3,i6,2e11.3/)
 50   FORMAT(/5x,"tau",7x,"taumax",7x,"t",9x,"tmax",6x,
     $     "mstep"//1p,4e11.3,i10/)
 60   FORMAT(/4x,"error",6x,"errmax",4x,"errstep",6x,"tol"
     $     //1p,4e11.3//)
c-----------------------------------------------------------------------
c     open output files, read, process, and diagnose equilibrium.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"orbit.out","UNKNOWN")
      CALL timer(0,out_unit)
      CALL equil_read(out_unit)
      CALL equil_out_global
      CALL equil_out_qfind
      CALL equil_out_diagnose(.TRUE.,out_unit)
      CALL equil_out_write_2d
c-----------------------------------------------------------------------
c     read orbit input and open output files.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"orbit.in","OLD")
      READ(in_unit,NML=orbit_input)
      CALL ascii_close(in_unit)
      IF(orbit_bin)CALL bin_open(bin_unit,"orbit.bin","UNKNOWN",
     $     "REWIND","none")
c-----------------------------------------------------------------------
c     print input parameters.
c-----------------------------------------------------------------------
      WRITE(out_unit,10)particle,ai,zi,TRIM(orbit_type),cross
      WRITE(out_unit,20)psi0,theta0,zeta0,r0,z0,phi0
      WRITE(out_unit,30)energy0,alpha0,phase0
      WRITE(out_unit,40)taumax,report,stride,tol,errmax
c-----------------------------------------------------------------------
c     start or restart.
c-----------------------------------------------------------------------
      IF(restart_flag)THEN
         CALL restart_read(istep,t,y)
      ELSE
         istep=0
         t=0
      ENDIF
      done=.FALSE.
c-----------------------------------------------------------------------
c     set particle constants.
c-----------------------------------------------------------------------
      SELECT CASE(TRIM(particle))
      CASE("ion")
         q=zi*e
         m=ai*mp
      CASE("electron")
         q=-e
         m=me
      CASE DEFAULT
         CALL program_stop
     $        ("Cannot recognize particle = "//TRIM(particle))
      END SELECT
c-----------------------------------------------------------------------
c     prepare fields.
c-----------------------------------------------------------------------
      WRITE(*,'(1x,a)')"Preparing fields"
      IF(gc_flag)THEN
         IF(orbit_type == "inner" .OR. cross)CALL field_gc_in
         IF(direct_flag .AND. (orbit_type == "outer" .OR. cross))
     $        CALL field_gc_out
      ELSE
         CALL field_1d
         CALL field_2d
      ENDIF
c-----------------------------------------------------------------------
c     compute orbit inside plasma region.
c-----------------------------------------------------------------------
      WRITE(*,'(1x,a)')"Computing orbit"
      DO
         SELECT CASE(TRIM(orbit_type))
         CASE("inner")
            IF(orbit_out)WRITE(out_unit,'(a)')
     $           "Computing orbit inside plasma region:"
            IF(gc_flag)THEN
               CALL inner_gc_run(istep,t,y_gc,done)
            ELSE
               CALL inner_run(istep,t,y,done)
            ENDIF
            IF(.NOT. done)orbit_type="outer"
c-----------------------------------------------------------------------
c     compute orbit outside plasma region.
c-----------------------------------------------------------------------
         CASE("outer")
            IF(direct_flag)THEN
               IF(orbit_out)WRITE(out_unit,'(a)')
     $              "Computing orbit outside plasma region:"
               IF(gc_flag)THEN
                  CALL outer_gc_run(istep,t,y_gc,done)
               ELSE
                  CALL outer_run(istep,t,y,done)
               ENDIF
               IF(.NOT. done)orbit_type="inner"
            ELSE
               done=.TRUE.
            ENDIF
c-----------------------------------------------------------------------
c     null and default cases.
c-----------------------------------------------------------------------
         CASE("none")
            done=.TRUE.
         CASE DEFAULT
            CALL program_stop
     $           ("Can't recognize orbit_type "//TRIM(orbit_type))
         END SELECT
c-----------------------------------------------------------------------
c     finish loop over orbit segments.
c-----------------------------------------------------------------------
         IF(done .OR. .NOT. cross)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     write final output.
c-----------------------------------------------------------------------
      tau=ABS(omega0)*t/twopi
      IF(istep == 0)THEN
         errstep=0
      ELSE
         errstep=error/istep
      ENDIF
      WRITE(out_unit,50)tau,taumax,t,tmax,istep
      WRITE(out_unit,60)error,errmax,errstep,tol
c-----------------------------------------------------------------------
c     deallocate arrays and close binary file.
c-----------------------------------------------------------------------
      DEALLOCATE(rwork)
      CALL bicube_dealloc(rzphi)
      CALL spline_dealloc(sq)
      IF(gc_flag)THEN
         IF(orbit_type == "inner" .OR. cross)
     $        CALL bicube_dealloc(gc_infield)
         IF(direct_flag .AND. (orbit_type == "outer" .OR. cross))
     $        CALL bicube_dealloc(gc_outfield)
      ELSE
         CALL bicube_dealloc(metric)
         CALL spline_dealloc(avec)
      ENDIF
      IF(direct_flag)CALL bicube_dealloc(psi_in)
      IF(orbit_bin)CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END PROGRAM orbit_main
