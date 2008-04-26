c-----------------------------------------------------------------------
c     file inner_gc.f
c     sets up and integrates the orbit equations inside the plasma.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inner_gc_mod.
c     1. inner_gc_run.
c     2. inner_gc_der.
c     3. inner_gc_energy.
c     4. inner_gc_init.
c     5. inner_gc_reinit.
c     6. inner_gc_output.
c     7. inner_gc_triple.
c-----------------------------------------------------------------------
c     subprogram 0. inner_gc_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inner_gc_mod
      USE field_mod
      USE restart_mod
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: neq=4,lrw=20+neq*16

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. inner_gc_run.
c     sets up and integrates differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_run(istep,t,y,done)
      
      INTEGER, INTENT(INOUT) :: istep
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y
      LOGICAL, INTENT(INOUT) :: done

      REAL(r8) :: told,psiold
      REAL(r8), DIMENSION(4) :: atol
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/6x,"i",5x,"tau",8x,"psi",7x,"theta",7x,"zeta",7x,"vpar",
     $     7x,"error"/)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      IF(t == 0 .OR. restart_flag)THEN
         CALL inner_gc_init(t,y)
         IF(y(1) > psihigh)THEN
            orbit_type="outer"
            RETURN
         ENDIF
      ELSE
         CALL inner_gc_reinit(y)
      ENDIF
      istate=1
      IF(orbit_out)WRITE(out_unit,10)
      IF(ode_bin)CALL bin_open(ode_bin_unit,"ode.bin",
     $     "UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     print data for each time step
c-----------------------------------------------------------------------
      DO
         IF(ode_bin)WRITE(ode_bin_unit)REAL(t,4),REAL(y,4)
         IF(mod(istep,stride) == 0 .OR. t >= tmax)
     $        CALL inner_gc_output(istep,t,y)
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
         IF(istep >= nstep .OR. t >= tmax .OR. istate < 0
     $        .OR. y(1) < psilow .OR. y(1) > psihigh+tol
     $        .OR. error > errmax)EXIT
         istep=istep+1
         told=t
         psiold=y(1)
         CALL lsode(inner_gc_der,neq,y,t,tmax,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     finish up.
c-----------------------------------------------------------------------
      IF(ode_bin)CALL bin_close(ode_bin_unit)
      IF(orbit_out)WRITE(out_unit,10)
      IF(y(1) < psihigh .AND. y(1) > psilow
     $     .OR. error > errmax)done=.TRUE.
      IF(orbit_bin)WRITE(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_gc_run
c-----------------------------------------------------------------------
c     subprogram 2. inner_gc_der.
c     contains the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_der(neq,t,y,dy)
      
      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(neq), INTENT(IN) :: y
      REAL(r8), DIMENSION(neq), INTENT(OUT) :: dy

      REAL(r8) :: jac,bmod,bgradb,curlbgradb,qbparfac,vpar
      REAL(r8), DIMENSION(3) :: curlb,gradb,bvec
c-----------------------------------------------------------------------
c     evaluate cubic splines.
c-----------------------------------------------------------------------
      CALL bicube_eval(gc_infield,y(1),y(2),1)
      CALL spline_eval(sq,y(1),0)
      bmod=gc_infield%f(4)
      jac=gc_infield%f(5)
      bvec=(/0._r8,1._r8,sq%f(4)/)*twopi*psio/jac
      vpar=y(4)
c-----------------------------------------------------------------------
c     compute curl b.
c-----------------------------------------------------------------------
      curlb(1)=gc_infield%fy(3)
      curlb(2)=-gc_infield%fx(3)
      curlb(3)=gc_infield%fx(2)-gc_infield%fy(1)
      curlb=curlb/jac
c-----------------------------------------------------------------------
c     compute curl b.
c-----------------------------------------------------------------------
      gradb(1)=-gc_infield%f(3)*gc_infield%fy(4)
      gradb(2)=gc_infield%f(3)*gc_infield%fx(4)
      gradb(3)=gc_infield%f(1)*gc_infield%fy(4)
     $     -gc_infield%f(2)*gc_infield%fx(4)
      gradb=gradb/jac
c-----------------------------------------------------------------------
c     composite factors.
c-----------------------------------------------------------------------
      bgradb=gc_infield%fy(4)*twopi*psio/jac
      curlbgradb=curlb(1)*gc_infield%fx(4)
     $     +curlb(2)*gc_infield%fy(4)
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
      END SUBROUTINE inner_gc_der
c-----------------------------------------------------------------------
c     subprogram 3. inner_gc_energy.
c     computes particle energy.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_energy(y,energy)

      REAL(r8), DIMENSION(:), INTENT(IN) :: y
      REAL(r8), INTENT(OUT) :: energy

      REAL(r8) :: bmod
c-----------------------------------------------------------------------
c     compute energy.
c-----------------------------------------------------------------------
      CALL bicube_eval(gc_infield,y(1),y(2),0)
      bmod=gc_infield%f(4)
      energy=mu*bmod+m*y(4)**2/2
      energy=energy/ev
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_gc_energy
c-----------------------------------------------------------------------
c     subprogram 4. inner_gc_init.
c     initializes dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_init(t,y)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(4), INTENT(OUT) :: y

      REAL(r8) :: vmod,vr,vz,vphi,bmod,chi1,eta,omegac,r,rfac,
     $     cosfac,sinfac
      REAL(r8), DIMENSION(3) :: vvec,bvec
      REAL(r8), DIMENSION(3,3) :: v,w
c-----------------------------------------------------------------------
c     compute coordinates.
c-----------------------------------------------------------------------
      IF(restart_flag)GOTO 10
      psi0=MIN(psi0,psihigh)
      CALL bicube_eval(rzphi,psi0,theta0,1)
      CALL spline_eval(sq,psi0,0)
      rfac=SQRT(rzphi%f(1))
      eta=twopi*(theta0+rzphi%f(2))
      cosfac=COS(eta)
      sinfac=SIN(eta)
      r=ro+rfac*cosfac
      jac=rzphi%f(4)
      chi1=twopi*psio
c-----------------------------------------------------------------------
c     compute velocity vector.
c-----------------------------------------------------------------------
      vmod=SQRT(2*energy0*ev/m)
      vr=vmod*SIN(dtor*alpha0)*COS(dtor*phase0)
      vz=vmod*COS(dtor*alpha0)
      vphi=vmod*SIN(dtor*alpha0)*SIN(dtor*phase0)
      vvec=(/vr*cosfac+vz*sinfac,-vr*sinfac+vz*cosfac,vphi/)
c-----------------------------------------------------------------------
c     compute magnetic field vector.
c-----------------------------------------------------------------------
      v=0
      v(2,1)=rzphi%fy(1)/(2*rfac*jac)
      v(2,2)=(1+rzphi%fy(2))*twopi*rfac/jac
      v(2,3)=rzphi%fy(3)*r/jac
      v(3,3)=twopi*r/jac
      bvec=chi1*(v(2,:)+sq%f(4)*v(3,:))
      bmod=SQRT(SUM(bvec**2))
      bvec=bvec/bmod
      omegac=q*bmod/m
c-----------------------------------------------------------------------
c     compute covariant basis vectors.
c-----------------------------------------------------------------------
      w=0
      w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r/jac
      w(1,2)=-rzphi%fy(1)*pi*r/(rfac*jac)
      w(2,1)=-rzphi%fx(2)*twopi**2*r*rfac/jac
      w(2,2)=rzphi%fx(1)*pi*r/(rfac*jac)
      w(3,1)=(rzphi%fx(2)*rzphi%fy(3)-rzphi%fx(3)*(1+rzphi%fy(2)))
     $     *twopi*r*rfac/jac
      w(3,2)=(rzphi%fx(3)*rzphi%fy(1)-rzphi%fx(1)*rzphi%fy(3))
     $     *r/(2*rfac*jac)
      w(3,3)=1/(twopi*r)
c-----------------------------------------------------------------------
c     compute initial conditions.
c-----------------------------------------------------------------------
      y(1:3)=(/psi0,theta0,zeta0/)+inner_gc_triple(vvec,bvec,w)/omegac
      y(4)=SUM(vvec*bvec)
c-----------------------------------------------------------------------
c     compute initial energy and zero tau_report.
c-----------------------------------------------------------------------
      CALL bicube_eval(gc_infield,y(1),y(2),0)
      bmod=gc_infield%f(4)
      mu=(SUM(vvec**2)-y(4)**2)*m/(2*bmod)
      CALL inner_gc_energy(y,energy0)
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
      END SUBROUTINE inner_gc_init
c-----------------------------------------------------------------------
c     subprogram 5. inner_gc_reinit.
c     reinitializes dependent variables from data outside separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_reinit(y)

      REAL(r8), DIMENSION(4), INTENT(INOUT) :: y

      INTEGER :: ipsi,itheta
      INTEGER, DIMENSION(1) :: jtheta
      REAL(r8) :: dtheta,eta0,f,f1,phi,r,z,psi
      REAL(r8), PARAMETER :: tol=1e-14
      REAL(r8), DIMENSION(0:mtheta) :: deta
c-----------------------------------------------------------------------
c     find coordinates.
c-----------------------------------------------------------------------
      r=y(1)
      z=y(2)
      phi=y(3)
      CALL bicube_eval(psi_in,r,z,0)
      psi=1-psi_in%f(1)/psio
      IF(ABS(psi-psihigh) < ABS(psi-psilow))THEN
         ipsi=mpsi
      ELSE
         ipsi=0
      ENDIF
      psi0=sq%xs(ipsi)
      eta0=ATAN2(z-zo,r-ro)/twopi
      IF(eta0 < 0)eta0=eta0+1
      deta=rzphi%ys+rzphi%fs(ipsi,:,2)-eta0
      jtheta=MINLOC(ABS(deta))
      itheta=jtheta(1)
      theta0=rzphi%ys(itheta)
c-----------------------------------------------------------------------
c     refine by Newton iteration.
c-----------------------------------------------------------------------
      DO
         CALL bicube_eval(rzphi,psi0,theta0,1)
         f=theta0+rzphi%f(2)-eta0
         f1=1+rzphi%fy(2)
         IF(ABS(f) < tol)EXIT
         dtheta=-f/f1
         theta0=theta0+dtheta
      ENDDO
      zeta0=(phi-rzphi%f(3))/twopi
c-----------------------------------------------------------------------
c     initialize dependent variables.
c-----------------------------------------------------------------------
      y(1:3)=(/psi0,theta0,zeta0/)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_gc_reinit
c-----------------------------------------------------------------------
c     subprogram 6. inner_gc_output.
c     produces ascii and binary output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_gc_output(istep,t,y)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(4), INTENT(IN) :: y

      REAL(r8) :: bmod,eta,phi,r,rfac,rl,rlfac,vpar,vprpsq,
     $     z,dmufac,psi,xx,yy,tau
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(i7,1p,6e11.3)
c-----------------------------------------------------------------------
c     evaluate splines and energy.
c-----------------------------------------------------------------------
      CALL bicube_eval(rzphi,y(1),y(2),0)
      CALL bicube_eval(gc_infield,y(1),y(2),0)
      CALL inner_gc_energy(y,energy)
c-----------------------------------------------------------------------
c     compute output quantities.
c-----------------------------------------------------------------------
      bmod=gc_infield%f(4)
      tau=ABS(omega0)*t/twopi
      error=ABS(energy/energy0-1)
      rfac=SQRT(rzphi%f(1))
      eta=twopi*(y(2)+rzphi%f(2))
      r=ro+rfac*COS(eta)
      z=zo+rfac*SIN(eta)
      phi=twopi*y(3)+rzphi%f(3)
      xx=r*COS(phi)
      yy=r*SIN(phi)
      vpar=y(4)
      vprpsq=2*mu*bmod/m
      rl=SQRT(vprpsq)*m/(q*bmod)
      rlfac=rl/r
      psi=y(1)
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
      END SUBROUTINE inner_gc_output
c-----------------------------------------------------------------------
c     subprogram 7. inner_gc_triple.
c     produces a triple product of 3-vectors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION inner_gc_triple(a,b,c) RESULT(triple)

      REAL(r8), DIMENSION(3), INTENT(IN) :: a,b
      REAL(r8), DIMENSION(3,3), INTENT(IN) :: c
      REAL(r8), DIMENSION(3) :: triple
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      triple
     $     =a(1)*(b(2)*c(:,3)-b(3)*c(:,2))
     $     +a(2)*(b(3)*c(:,1)-b(1)*c(:,3))
     $     +a(3)*(b(1)*c(:,2)-b(2)*c(:,1))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION inner_gc_triple
      END MODULE inner_gc_mod
