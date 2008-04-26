c-----------------------------------------------------------------------
c     file inner.f
c     sets up and integrates the orbit equations inside the plasma.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inner_mod.
c     1. inner_run.
c     2. inner_der.
c     3. inner_energy.
c     4. inner_init.
c     5. inner_reinit.
c     6. inner_output.
c-----------------------------------------------------------------------
c     subprogram 0. inner_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inner_mod
      USE field_mod
      USE restart_mod
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: neq=6,lrw=20+neq*16

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. inner_run.
c     sets up and integrates differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_run(istep,t,y,done)
      
      INTEGER, INTENT(INOUT) :: istep
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y
      LOGICAL, INTENT(INOUT) :: done

      REAL(r8) :: told,psiold
      REAL(r8), DIMENSION(3,2) :: atol
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/6x,"i",5x,"tau",8x,"psi",7x,"theta",7x,"zeta",6x,"p_psi",
     $     5x,"p_theta",4x,"p_zeta",7x,"error"/)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      IF(t == 0 .OR. restart_flag)THEN
         CALL inner_init(t,y)
      ELSE
         CALL inner_reinit(y)
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
     $        CALL inner_output(istep,t,y)
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
     $        .OR. y(1,1) < psilow .OR. y(1,1) > psihigh+tol
     $        .OR. error > errmax)EXIT
         istep=istep+1
         told=t
         psiold=y(1,1)
         CALL lsode(inner_der,neq,y,t,tmax,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     finish up.
c-----------------------------------------------------------------------
      IF(ode_bin)CALL bin_close(ode_bin_unit)
      IF(orbit_out)WRITE(out_unit,10)
      IF(y(1,1) < psihigh .AND. y(1,1) > psilow
     $     .OR. error > errmax)done=.TRUE.
      IF(orbit_bin)WRITE(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_run
c-----------------------------------------------------------------------
c     subprogram 2. inner_der.
c     contains the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_der(neq,t,y,dy)
      
      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(IN) :: y
      REAL(r8), DIMENSION(3,2), INTENT(OUT) :: dy

      INTEGER :: info
      INTEGER, DIMENSION(3) :: ipiv
      REAL(r8), DIMENSION(3) :: vvec,vvecx,vvecy
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     evaluate splines.
c-----------------------------------------------------------------------
      CALL field_get_jac(y(1,1),y(2,1),1,jac,jacx,jacy)
      CALL spline_eval(avec,y(1,1),1)
c-----------------------------------------------------------------------
c     compute vvec and its derivatives.
c-----------------------------------------------------------------------
      CALL dgetrf(3,3,jac,3,ipiv,info)
      vvec=(y(:,2)-(/zero,avec%f(1),avec%f(2)/)*q)/m
      CALL dgetrs("T",3,1,jac,3,ipiv,vvec,3,info)
      vvecx=-MATMUL(TRANSPOSE(jacx),vvec)
      vvecx=vvecx-(/zero,avec%f1(1),avec%f1(2)/)*q/m
      CALL dgetrs("T",3,1,jac,3,ipiv,vvecx,3,info)
      vvecy=-MATMUL(TRANSPOSE(jacy),vvec)
      CALL dgetrs("T",3,1,jac,3,ipiv,vvecy,3,info)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(:,1)=vvec
      CALL dgetrs("N",3,1,jac,3,ipiv,dy(:,1),3,info)
      dy(1,2)=-SUM(vvec*vvecx)*m
      dy(2,2)=-SUM(vvec*vvecy)*m
      dy(3,2)=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_der
c-----------------------------------------------------------------------
c     subprogram 3. inner_energy.
c     computes particle energy.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_energy(y,energy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: y
      REAL(r8), INTENT(OUT) :: energy

      INTEGER :: info
      INTEGER, DIMENSION(3) :: ipiv
      REAL(r8), DIMENSION(3) :: vvec
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     evaluate splines.
c-----------------------------------------------------------------------
      CALL field_get_jac(y(1,1),y(2,1),0,jac,jacx,jacy)
      CALL spline_eval(avec,y(1,1),0)
c-----------------------------------------------------------------------
c     compute energy.
c-----------------------------------------------------------------------
      vvec=(y(:,2)-q*(/zero,avec%f(1),avec%f(2)/))/m
      CALL dgetrf(3,3,jac,3,ipiv,info)
      CALL dgetrs("T",3,1,jac,3,ipiv,vvec,3,info)
      energy=SUM(vvec**2)*m/(2*ev)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_energy
c-----------------------------------------------------------------------
c     subprogram 4. inner_init.
c     initializes dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_init(t,y)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(OUT) :: y

      REAL(r8) :: v,vr,vz,vphi
      REAL(r8), DIMENSION(3) :: vvec
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      IF(.NOT. restart_flag)THEN
         psi0=MIN(psi0,psihigh)
         CALL field_get_jac(psi0,theta0,0,jac,jacx,jacy)
         CALL spline_eval(avec,psi0,0)
         v=SQRT(2*energy0*ev/m)
         vr=v*SIN(dtor*alpha0)*COS(dtor*phase0)
         vz=v*COS(dtor*alpha0)
         vphi=v*SIN(dtor*alpha0)*SIN(dtor*phase0)
         vvec=(/vr,vz,vphi/)
         y(:,1)=(/psi0,theta0,zeta0/)
         y(:,2)=m*MATMUL(TRANSPOSE(jac),vvec)
     $        +q*(/zero,avec%f(1),avec%f(2)/)
         CALL inner_energy(y,energy0)
         tau_report=0
      ENDIF
c-----------------------------------------------------------------------
c     initialize remaining variables.
c-----------------------------------------------------------------------
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
      END SUBROUTINE inner_init
c-----------------------------------------------------------------------
c     subprogram 5. inner_reinit.
c     reinitializes dependent variables from data outside separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_reinit(y)

      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y

      INTEGER :: ipsi,itheta
      INTEGER, DIMENSION(1) :: jtheta
      REAL(r8) :: az,dtheta,eta0,f,f1,phi,r,vphi,vr,vz,z,psi
      REAL(r8), PARAMETER :: tol=1e-14
      REAL(r8), DIMENSION(3) :: vvec
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
      REAL(r8), DIMENSION(0:mtheta) :: deta
c-----------------------------------------------------------------------
c     find coordinates.
c-----------------------------------------------------------------------
      r=y(1,1)
      z=y(2,1)
      phi=y(3,1)
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
c     initialize velocities.
c-----------------------------------------------------------------------
      az=f0*LOG(r)
      vr=y(1,2)/m
      vz=(y(2,2)-q*az)/m
      vphi=(y(3,2)-q*psi_in%f(1))/(m*r)
      vvec=(/vr,vz,vphi/)
c-----------------------------------------------------------------------
c     evaluate bicube splines.
c-----------------------------------------------------------------------
      CALL field_get_jac(psi0,theta0,0,jac,jacx,jacy)
      CALL spline_eval(avec,psi0,0)
c-----------------------------------------------------------------------
c     initialize dependent variables.
c-----------------------------------------------------------------------
      y(:,1)=(/psi0,theta0,zeta0/)
      y(:,2)=m*MATMUL(TRANSPOSE(jac),vvec)
     $     +q*(/zero,avec%f(1),avec%f(2)/)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_reinit
c-----------------------------------------------------------------------
c     subprogram 6. inner_output.
c     produces ascii and binary output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inner_output(istep,t,y)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(IN) :: y

      INTEGER :: info
      INTEGER, DIMENSION(3) :: ipiv
      REAL(r8) :: bmod,eta,phi,r,rfac,rl,rlfac,vpar,vprpsq,vsq,
     $     z,dmufac,psi,xx,yy,tau
      REAL(r8), DIMENSION(3) :: vvec,bvec
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(i7,1p,8e11.3)
c-----------------------------------------------------------------------
c     compute velocity and energy.
c-----------------------------------------------------------------------
      tau=ABS(omega0)*t/twopi
      CALL inner_energy(y,energy)
      error=ABS(energy/energy0-1)
c-----------------------------------------------------------------------
c     evaluate splines.
c-----------------------------------------------------------------------
      CALL bicube_eval(rzphi,y(1,1),y(2,1),0)
      CALL spline_eval(avec,y(1,1),1)
      CALL bicube_eval(rzphi,y(1,1),y(2,1),0)
      CALL field_get_jac(y(1,1),y(2,1),0,jac,jacx,jacy)
c-----------------------------------------------------------------------
c     compute cylindrical coordinates.
c-----------------------------------------------------------------------
      rfac=SQRT(rzphi%f(1))
      eta=twopi*(y(2,1)+rzphi%f(2))
      r=ro+rfac*COS(eta)
      z=zo+rfac*SIN(eta)
      phi=twopi*y(3,1)+rzphi%f(3)
      xx=r*COS(phi)
      yy=r*SIN(phi)
c-----------------------------------------------------------------------
c     compute magnetic moment.
c-----------------------------------------------------------------------
      bvec=(/zero,-avec%f1(2),avec%f1(1)/)/rzphi%f(4)
      bvec=MATMUL(jac,bvec)
      bmod=SQRT(SUM(bvec*bvec))
      CALL dgetrf(3,3,jac,3,ipiv,info)
      vvec=(y(:,2)-(/zero,avec%f(1),avec%f(2)/)*q)/m
      CALL dgetrs("T",3,1,jac,3,ipiv,vvec,3,info)
      vsq=SUM(vvec**2)
      vpar=SUM(vvec*bvec)/bmod
      vprpsq=vsq-vpar**2
      mu=vprpsq/bmod
      IF(istep == 0)mu_0=mu
      dmufac=mu/mu_0-1
c-----------------------------------------------------------------------
c     compute Larmor radius and poloidal flux.
c-----------------------------------------------------------------------
      rl=SQRT(vprpsq)*m/(q*bmod)
      rlfac=rl/r
      psi=y(1,1)
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
         CALL restart_write(istep,t,y)
         tau_report=tau_report+report
         IF(tau > 0 .AND. orbit_bin .AND. break_flag)WRITE(bin_unit)
         restart_flag=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inner_output
      END MODULE inner_mod
