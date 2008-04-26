c-----------------------------------------------------------------------
c     file outer.f
c     sets up and integrates the orbit equations outside the plasma.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. outer_mod.
c     1. outer_run.
c     2. outer_der.
c     3. outer_energy.
c     4. outer_init.
c     5. outer_reinit.
c     6. outer_output.
c     7. outer_bound_check.
c-----------------------------------------------------------------------
c     subprogram 0. outer_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE outer_mod
      USE field_mod
      USE restart_mod
      IMPLICIT NONE

      CHARACTER(7), PRIVATE :: region
      INTEGER, PARAMETER, PRIVATE :: neq=6,lrw=20+neq*16

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. outer_run.
c     sets up and integrates differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_run(istep,t,y,done)

      INTEGER, INTENT(INOUT) :: istep
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y
      LOGICAL, INTENT(INOUT) :: done

      LOGICAL :: bound
      REAL(r8), DIMENSION(3,2) :: atol
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/6x,"i",5x,"tau",9x,"r",10x,"z",9x,"phi",8x,"p_r",8x,
     $     "p_z",8x,"p_phi",6x,"error"/)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      IF(t == 0 .OR. restart_flag)THEN
         CALL outer_init(t,y)
      ELSE
         CALL outer_reinit(y)
      ENDIF
      istate=1
      bound=.FALSE.
      IF(orbit_out)WRITE(out_unit,10)
c-----------------------------------------------------------------------
c     print data for each time step
c-----------------------------------------------------------------------
      DO
         CALL outer_bound_check(istep,t,y,bound)
         IF(mod(istep,stride) == 0 .OR. t >= tmax)
     $        CALL outer_output(istep,t,y)
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
     $        .OR. y(1,1) > rmax .OR. y(1,1) < rmin
     $        .OR. y(2,1) > zmax .OR. y(2,1) < zmin
     $        .OR. ABS(error) > errmax)THEN
            done=.TRUE.
            EXIT
         ENDIF
         istep=istep+1
         CALL lsode(outer_der,neq,y,t,tmax,itol,rtol,atol,
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
      END SUBROUTINE outer_run
c-----------------------------------------------------------------------
c     subprogram 2. outer_der.
c     contains the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_der(neq,t,y,dy)
      
      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(IN) :: y
      REAL(r8), DIMENSION(3,2), INTENT(OUT) :: dy

      REAL(r8) :: r,az,az_r,psi,psi_r,psi_z
c-----------------------------------------------------------------------
c     compute fields.
c-----------------------------------------------------------------------
      CALL bicube_eval(psi_in,y(1,1),y(2,1),1)
      psi=psi_in%f(1)
      psi_r=psi_in%fx(1)
      psi_z=psi_in%fy(1)
      r=y(1,1)
      az=f0*LOG(r)
      az_r=f0/r
c-----------------------------------------------------------------------
c     compute velocities.
c-----------------------------------------------------------------------
      dy(1,1)=y(1,2)/m
      dy(2,1)=(y(2,2)-q*az)/m
      dy(3,1)=(y(3,2)-q*psi)/(m*r**2)
c-----------------------------------------------------------------------
c     compute forces.
c-----------------------------------------------------------------------
      dy(1,2)=q*(dy(2,1)*az_r+dy(3,1)*psi_r)+m*r*dy(3,1)**2
      dy(2,2)=q*dy(3,1)*psi_z
      dy(3,2)=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_der
c-----------------------------------------------------------------------
c     subprogram 3. outer_energy.
c     computes particle energy.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_energy(y,energy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: y
      REAL(r8), INTENT(OUT) :: energy

      REAL(R8) :: t=0
      REAL(r8), DIMENSION(SIZE(y,1),SIZE(y,2)) :: dy
c-----------------------------------------------------------------------
c     compute momentum vector and energy.
c-----------------------------------------------------------------------
      CALL outer_der(neq,t,y,dy)
      energy=(dy(1,1)**2+dy(2,1)**2+(dy(3,1)*y(1,1))**2)*m/(2*ev)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_energy
c-----------------------------------------------------------------------
c     subprogram 4. outer_init.
c     initializes dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_init(t,y)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y

      REAL(r8) :: v,vr,vphi,vz,psi,r,z,az
c-----------------------------------------------------------------------
c     determine region.
c-----------------------------------------------------------------------
      IF(restart_flag)THEN
         r=y(1,1)
         z=y(2,1)
      ELSE
         r=rs2+r0*(rmax-rs2)
         z=(zmax+zmin)/2+z0*(zmax-zmin)/2
      ENDIF
      CALL bicube_eval(psi_in,r,z,0)
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
c     initialize dependent variables.
c-----------------------------------------------------------------------
      f0=sq%fs(mpsi,1)/twopi
      IF(.NOT. restart_flag)THEN
         psi=psi_in%f(1)
         az=f0*LOG(r)
         v=SQRT(2*energy0*ev/m)
         vr=v*SIN(dtor*alpha0)*COS(dtor*phase0)
         vz=v*COS(dtor*alpha0)
         vphi=v*SIN(dtor*alpha0)*SIN(dtor*phase0)
         y(:,1)=(/r,z,phi0/)
         y(:,2)=(/m*vr,m*vz+q*az,m*r*vphi+q*psi/)
         CALL outer_energy(y,energy0)
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
      END SUBROUTINE outer_init
c-----------------------------------------------------------------------
c     subprogram 5. outer_reinit.
c     reinitializes dependent variables from data inside separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_reinit(y)

      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y

      INTEGER :: info
      INTEGER, DIMENSION(3) :: ipiv
      REAL(r8) :: rfac,eta,r,z,phi,psi0,psi,az,dpsi,det,dr,dz,tol=1e-14
      REAL(r8), DIMENSION(3) :: vvec
      REAL(r8), DIMENSION(3,3) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     evaluate splines.
c-----------------------------------------------------------------------
      CALL bicube_eval(rzphi,y(1,1),y(2,1),0)
      CALL spline_eval(avec,y(1,1),0)
      CALL spline_eval(sq,y(1,1),0)
c-----------------------------------------------------------------------
c     compute cylindrical coordinates.
c-----------------------------------------------------------------------
      rfac=SQRT(rzphi%f(1))
      eta=twopi*(y(2,1)+rzphi%f(2))
      r=ro+rfac*COS(eta)
      z=zo+rfac*SIN(eta)
      phi=twopi*y(3,1)+rzphi%f(3)
      f0=sq%f(1)/twopi
c-----------------------------------------------------------------------
c     use Newton iteration to refine flux surface.
c-----------------------------------------------------------------------
      psi0=psio*(1-y(1,1))
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
c     compute velocity vvec.
c-----------------------------------------------------------------------
      CALL field_get_jac(y(1,1),y(2,1),0,jac,jacx,jacy)
      CALL dgetrf(3,3,jac,3,ipiv,info)
      vvec=(y(:,2)-(/0._r8,avec%f(1),avec%f(2)/)*q)/m
      CALL dgetrs("T",3,1,jac,3,ipiv,vvec,3,info)
      vvec(3)=vvec(3)*r
c-----------------------------------------------------------------------
c     initialize dependent variables.
c-----------------------------------------------------------------------
      psi=psio*(1-y(1,1))
      az=f0*LOG(r)
      y(:,1)=(/r,z,phi/)
      y(:,2)=m*vvec+q*(/0._r8,az,psi/)
c-----------------------------------------------------------------------
c     determine region.
c-----------------------------------------------------------------------
      IF(1-psi/psio >= psihigh)THEN
         region="outside"
      ELSE
         region="inside"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE outer_reinit
c-----------------------------------------------------------------------
c     subprogram 6. outer_output.
c     produces ascii and binary output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_output(istep,t,y)

      INTEGER, INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(3,2), INTENT(IN) :: y

      REAL(r8) :: bmod,phi,r,rl,rlfac,vpar,vprpsq,vsq,z,dmufac,
     $     br,bz,bt,psi,az,xx,yy,tau
      REAL(r8), DIMENSION(3) :: vvec,bvec
      REAL(r8), DIMENSION(3,2) :: dy
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(i7,1p,8e11.3)
c-----------------------------------------------------------------------
c     compute velocity and energy.
c-----------------------------------------------------------------------
      tau=ABS(omega0)*t/twopi
      CALL outer_energy(y,energy)
      error=ABS(energy/energy0-1)
c-----------------------------------------------------------------------
c     compute cylindrical coordinates.
c-----------------------------------------------------------------------
      r=y(1,1)
      z=y(2,1)
      phi=y(3,1)
      xx=r*COS(phi)
      yy=r*SIN(phi)
c-----------------------------------------------------------------------
c     compute fields.
c-----------------------------------------------------------------------
      CALL bicube_eval(psi_in,r,z,0)
      br=psi_in%fy(1)/r
      bz=-psi_in%fx(1)/r
      bt=f0/r
      az=f0*LOG(r)
c-----------------------------------------------------------------------
c     compute magnetic moment.
c-----------------------------------------------------------------------
      bvec=(/br,bz,bt/)
      bmod=SQRT(SUM(bvec**2))
      CALL outer_der(neq,t,y,dy)
      vvec=dy(:,1)
      vvec(3)=vvec(3)*r
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
      psi=1-psi_in%f(1)/psio
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
     $     .OR. y(1,1) < psilow .OR. error > errmax)THEN
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
      END SUBROUTINE outer_output
c-----------------------------------------------------------------------
c     subprogram 7. outer_bound_check.
c     checks for separatrix crossing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE outer_bound_check(istep,t,y,bound)

      INTEGER, INTENT(IN) :: istep

      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(3,2), INTENT(INOUT) :: y
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
      r=y(1,1)
      z=y(2,1)
      CALL bicube_eval(psi_in,y(1,1),y(2,1),0)
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
      END SUBROUTINE outer_bound_check
      END MODULE outer_mod
