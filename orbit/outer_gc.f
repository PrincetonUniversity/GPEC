c-----------------------------------------------------------------------
c     file outer_gc.f
c     sets up and integrates the orbit equations inside the plasma.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. outer_gc_mod.
c     1. outer_gc_run.
c     2. outer_gc_der.
c     3. outer_gc_energy.
c     4. outer_gc_init.
c     5. outer_gc_reinit.
c     6. outer_gc_output.
c     7. outer_gc_bound_check.
c-----------------------------------------------------------------------
c     subprogram 0. outer_gc_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE outer_gc_mod
      USE field_mod
      USE restart_mod
      IMPLICIT NONE

      CHARACTER(7), PRIVATE :: region
      INTEGER, PARAMETER, PRIVATE :: neq=4,lrw=20+neq*16

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. outer_gc_run.
c     sets up and integrates differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_run(istep,t,y,done)
      
      INTEGER, INTENT(INOUT) :: istep
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y
      LOGICAL, INTENT(INOUT) :: done

      LOGICAL :: bound
      REAL(r8), DIMENSION(4) :: atol
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/6x,"i",5x,"tau",9x,"r",10x,"z",9x,"phi",8x,"vpar",
     $     7x,"error"/)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      IF(t == 0 .OR. restart_flag)THEN
         CALL outer_gc_init(t,y)
      ELSE
         CALL outer_gc_reinit(y)
      ENDIF
      istate=1
      bound=.FALSE.
      IF(orbit_out)WRITE(out_unit,10)
c-----------------------------------------------------------------------
c     print data for each time step
c-----------------------------------------------------------------------
      DO
         CALL outer_gc_bound_check(istep,t,y,bound)
         IF(mod(istep,stride) == 0 .OR. t >= tmax)
     $        CALL outer_gc_output(istep,t,y)
c-----------------------------------------------------------------------
c     compute absolute tolerances.
c-----------------------------------------------------------------------
         WHERE(y == 0)
            atol=tol*tol
         ELSEWHERE
            atol=0
         ENDWHERE
c-----------------------------------------------------------------------
c     advance differential equations.
c-----------------------------------------------------------------------
         IF(istep >= nstep .OR. t >= tmax .OR. istate < 0 .OR. bound
     $        .OR. y(1) > rmax .OR. y(1) < rmin
     $        .OR. y(2) > zmax .OR. y(2) < zmin
     $        .OR. ABS(error) > errmax)THEN
            done=.TRUE.
            EXIT
         ENDIF
         istep=istep+1
         CALL lsode(outer_gc_der,neq,y,t,tmax,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     finish up.
c-----------------------------------------------------------------------
      IF(orbit_bin)WRITE(bin_unit)
      IF(orbit_out)WRITE(out_unit,10)
      IF(bound)done=.FALSE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_run
c-----------------------------------------------------------------------
c     subprogram 2. outer_gc_der.
c     contains the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_der(neq,t,y,dy)
      
      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(neq), INTENT(IN) :: y
      REAL(r8), DIMENSION(neq), INTENT(OUT) :: dy

      REAL(r8) :: r,bmod,bgradb,curlbgradb,qbparfac,vpar
      REAL(r8), DIMENSION(3) :: curlb,gradb,bvec
c-----------------------------------------------------------------------
c     evaluate cubic splines.
c-----------------------------------------------------------------------
      r=y(1)
      vpar=y(4)
      CALL bicube_eval(gc_outfield,y(1),y(2),1)
      bmod=gc_outfield%f(4)
      bvec=gc_outfield%f(1:3)*bmod
      bvec(3)=bvec(3)/r**2
c-----------------------------------------------------------------------
c     compute curl b.
c-----------------------------------------------------------------------
      curlb(1)=gc_outfield%fy(3)
      curlb(2)=-gc_outfield%fx(3)
      curlb(3)=gc_outfield%fx(2)-gc_outfield%fy(1)
      curlb=curlb/r
c-----------------------------------------------------------------------
c     compute curl b.
c-----------------------------------------------------------------------
      gradb(1)=-gc_outfield%f(3)*gc_outfield%fy(4)
      gradb(2)=gc_outfield%f(3)*gc_outfield%fx(4)
      gradb(3)=gc_outfield%f(1)*gc_outfield%fy(4)
     $     -gc_outfield%f(2)*gc_outfield%fx(4)
      gradb=gradb/r
c-----------------------------------------------------------------------
c     composite factors.
c-----------------------------------------------------------------------
      bgradb=bvec(1)*gc_outfield%fx(4)+bvec(2)*gc_outfield%fy(4)
      curlbgradb=curlb(1)*gc_outfield%fx(4)
     $     +curlb(2)*gc_outfield%fy(4)
      qbparfac=q*bmod+m*vpar*curlbgradb
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1:3)=(vpar*(q*bvec+m*vpar*curlb)+mu*gradb)/qbparfac
      dy(4)=-(q*bgradb+m*vpar*curlbgradb)*mu/(m*qbparfac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_der
c-----------------------------------------------------------------------
c     subprogram 3. outer_gc_energy.
c     computes particle energy.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_energy(y,energy)

      REAL(r8), DIMENSION(:), INTENT(IN) :: y
      REAL(r8), INTENT(OUT) :: energy

      REAL(r8) :: bmod
c-----------------------------------------------------------------------
c     compute energy.
c-----------------------------------------------------------------------
      CALL bicube_eval(gc_outfield,y(1),y(2),0)
      bmod=gc_outfield%f(4)
      energy=mu*bmod+m*y(4)**2/2
      energy=energy/ev
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_energy
c-----------------------------------------------------------------------
c     subprogram 4. outer_gc_init.
c     initializes dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_init(t,y)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y

      REAL(r8) :: v,vr,vphi,vz,psi,r,z,br,bz,bt,bmod,omegac
      REAL(r8), DIMENSION(3) :: vvec,bvec
c-----------------------------------------------------------------------
c     determine region.
c-----------------------------------------------------------------------
      IF(restart_flag)THEN
         r=y(1)
         z=y(2)
      ELSE
         r=rs2+r0*(rmax-rs2)
         z=(zmax+zmin)/2+z0*(zmax-zmin)/2
      ENDIF
      CALL bicube_eval(psi_in,r,z,1)
      psi=1-psi_in%f(1)/psio
      IF(r > rsep(1) .OR. r < rsep(2)
     $     .OR. z > zsep(1) .OR. z < zsep(2))THEN
         region="outside"
      ELSE IF(psi >= psihigh)THEN
         region="outside"
      ELSE
         region="inside"
      ENDIF
c-----------------------------------------------------------------------
c     compute velocity.
c-----------------------------------------------------------------------
      IF(restart_flag)GOTO 10
      v=SQRT(2*energy0*ev/m)
      vr=v*SIN(dtor*alpha0)*COS(dtor*phase0)
      vz=v*COS(dtor*alpha0)
      vphi=v*SIN(dtor*alpha0)*SIN(dtor*phase0)
      vvec=(/vr,vz,vphi/)
c-----------------------------------------------------------------------
c     compute magnetic field.
c-----------------------------------------------------------------------
      br=psi_in%fy(1)/r
      bz=-psi_in%fx(1)/r
      bt=sq%fs(mpsi,1)/(twopi*r)
      bvec=(/br,bz,bt/)
      bmod=SQRT(SUM(bvec*bvec))
      bvec=bvec/bmod
      omegac=q*bmod/m
c-----------------------------------------------------------------------
c     compute dependent variables.
c-----------------------------------------------------------------------
      y(1:3)=(/r,z,phi0/)
     $     +(/vz*bt-vphi*bz,vphi*br-vr*bt,vr*bz-vz*br/)/omegac
      y(4)=SUM(vvec*bvec)
c-----------------------------------------------------------------------
c     compute initial energy and zero tau_report.
c-----------------------------------------------------------------------
      CALL bicube_eval(gc_outfield,y(1),y(2),0)
      bmod=gc_outfield%f(4)
      mu=(SUM(vvec**2)-y(4)**2)*m/(2*bmod)
      CALL outer_gc_energy(y,energy0)
      tau_report=0
c-----------------------------------------------------------------------
c     initialize remaining variables.
c-----------------------------------------------------------------------
 10   CONTINUE
      omega0=SQRT(2*energy0*ev/m)/ro
      tmax=taumax*twopi/ABS(omega0)+t
      iopt=1
      mf=10
      itol=2
      rtol=tol
      itask=5
      iwork=0
      ALLOCATE(rwork(lrw))
      rwork=0
      rwork(1)=tmax
      rwork(11)=1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_init
c-----------------------------------------------------------------------
c     subprogram 5. outer_gc_reinit.
c     reinitializes dependent variables from data inside separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_reinit(y)

      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y

      REAL(r8) :: rfac,eta,r,z,phi,psi0,psi,dpsi,det,dr,dz,tol=1e-14
c-----------------------------------------------------------------------
c     compute cylindrical coordinates.
c-----------------------------------------------------------------------
      CALL bicube_eval(rzphi,y(1),y(2),0)
      rfac=SQRT(rzphi%f(1))
      eta=twopi*(y(2)+rzphi%f(2))
      r=ro+rfac*COS(eta)
      z=zo+rfac*SIN(eta)
      phi=twopi*y(3)+rzphi%f(3)
c-----------------------------------------------------------------------
c     use Newton iteration to refine flux surface.
c-----------------------------------------------------------------------
      psi0=psio*(1-y(1))
      DO
         CALL bicube_eval(psi_in,r,z,1)
         dpsi=psi0-psi_in%f(1)
         det=psi_in%fx(1)**2+psi_in%fy(1)**2
         dr=dpsi*psi_in%fx(1)/det
         dz=dpsi*psi_in%fy(1)/det
         r=r+dr
         z=z+dz
         error=ABS(dpsi/psio)
         IF(error < tol)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     initialize dependent variables.
c-----------------------------------------------------------------------
      y(1:3)=(/r,z,phi/)
c-----------------------------------------------------------------------
c     determine region.
c-----------------------------------------------------------------------
      psi=psio*(1-y(1))
      IF(1-psi/psio >= psihigh)THEN
         region="outside"
      ELSE
         region="inside"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_reinit
c-----------------------------------------------------------------------
c     subprogram 6. outer_gc_output.
c     produces ascii and binary output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_output(istep,t,y)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(4), INTENT(IN) :: y

      REAL(r8) :: bmod,phi,r,rl,rlfac,vpar,vprpsq,z,dmufac,psi,
     $     xx,yy,tau
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(i7,1p,6e11.3)
c-----------------------------------------------------------------------
c     evaluate splines and energy.
c-----------------------------------------------------------------------
      CALL bicube_eval(psi_in,y(1),y(2),0)
      CALL bicube_eval(gc_outfield,y(1),y(2),0)
      CALL outer_gc_energy(y,energy)
c-----------------------------------------------------------------------
c     compute output quantities.
c-----------------------------------------------------------------------
      bmod=gc_outfield%f(4)
      tau=ABS(omega0)*t/twopi
      error=ABS(energy/energy0-1)
      r=y(1)
      z=y(2)
      phi=y(3)
      xx=r*COS(phi)
      yy=r*SIN(phi)
      vpar=y(4)
      vprpsq=2*mu*bmod/m
      rl=SQRT(vprpsq)*m/(q*bmod)
      rlfac=rl/r
      psi=1-psi_in%f(1)/psio
      dmufac=0
c-----------------------------------------------------------------------
c     write output to files.
c-----------------------------------------------------------------------
      IF(orbit_out)WRITE(out_unit,10)istep,tau,y,error
      IF(orbit_bin)WRITE(bin_unit)REAL(tau,4),REAL(r,4),REAL(z,4),
     $     REAL(phi,4),REAL(error,4),REAL(bmod,4),REAL(rl,4),
     $     REAL(rlfac,4),REAL(dmufac,4),REAL(psi,4),REAL(SQRT(psi),4),
     $     REAL(vpar,4),REAL(xx,4),REAL(yy,4)
c-----------------------------------------------------------------------
c     report progress to screen.
c-----------------------------------------------------------------------
      IF(report > 0 .AND. tau >= tau_report .OR. restart_flag
     $     .OR. istep >= nstep .OR. t >= tmax .OR. istate < 0
     $     .OR. error > errmax)THEN
         WRITE(*,'(1x,a,i7,1p,2(a,e10.3))')
     $        "istep = ",istep,", tau = ",tau,", error = ",error
c         CALL restart_write(istep,t,y)
         tau_report=tau_report+report
         IF(tau > 0 .AND. orbit_bin .AND. break_flag)WRITE(bin_unit)
         restart_flag=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_output
c-----------------------------------------------------------------------
c     subprogram 7. outer_gc_bound_check.
c     checks for separatrix crossing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_gc_bound_check(istep,t,y,bound)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y
      LOGICAL, INTENT(OUT) :: bound

      CHARACTER(64) :: message1,message2,message3
      CHARACTER(128) :: message
      INTEGER :: iflag
      REAL(r8), SAVE :: psiold,told
      REAL(r8) :: psi,r,z
      REAL(r8), PARAMETER :: eps=1e-3
c-----------------------------------------------------------------------
c     compute local value of psi.
c-----------------------------------------------------------------------
      r=y(1)
      z=y(2)
      CALL bicube_eval(psi_in,r,z,0)
      psi=1-psi_in%f(1)/psio
      iflag=0
      bound=.FALSE.
c-----------------------------------------------------------------------
c     interpolate to outer boundary.
c-----------------------------------------------------------------------
      IF(psi < psihigh .AND. istate /= 1 .AND. region == "outside"
     $     .AND. r < rsep(1)+eps .AND. r > rsep(2)-eps
     $     .AND. z < zsep(1)+eps .AND. z > zsep(2)-eps)THEN
         t=told+(t-told)*(psihigh-psiold)/(psi-psiold)
         CALL dintdy(t,0,rwork(21),neq,y,iflag)
         bound=.TRUE.
      ENDIF
c-----------------------------------------------------------------------
c     interpolate to inner boundary.
c-----------------------------------------------------------------------
      IF(psi > psilow .AND. istate /= 1 .AND. region == "inside")THEN
         t=told+(t-told)*(psilow-psiold)/(psi-psiold)
         CALL dintdy(t,0,rwork(21),neq,y,iflag)
         bound=.TRUE.
      ENDIF
c-----------------------------------------------------------------------
c     check for interpolation error.
c-----------------------------------------------------------------------
      IF(iflag < 0)THEN
         WRITE(message1,'(a,i2,a,i6)')
     $        "Termination by dintdy, iflag = ",iflag,", istep = ",istep
         WRITE(message2,'(1p,2(a,e13.6))')"t = ",t,", told = ",told
         WRITE(message3,'(1p,2(a,e13.6))')
     $        "psi = ",psi,", psiold = ",psiold
         message=TRIM(message1)
     $        //CHAR(10)//TRIM(message2)//CHAR(10)//TRIM(message3)
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     finish up.
c-----------------------------------------------------------------------
      told=t
      psiold=psi
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_gc_bound_check
      END MODULE outer_gc_mod
