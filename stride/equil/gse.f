c-----------------------------------------------------------------------
c     file gse.f.
c     solves cylindrical Grad_Shafranov equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. gse_mod.
c     1. gse_run.
c     2. gse_int.
c     3. gse_der.
c     4. gse_prof.
c-----------------------------------------------------------------------
c     subprogram 0. gse_mod.
c     profile declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gse_mod
      USE inverse_mod
      IMPLICIT NONE

      LOGICAL, PARAMETER :: diagnose=.TRUE.
      INTEGER :: ma
      REAL(r8), PRIVATE :: chi1,a,r0,b0=1,beta0,p_pres,p_sig,tol=1e-12
      TYPE(spline_type) :: eq

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gse_run.
c     controls 1D Grad-Shafranov solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gse_run

      INTEGER :: it,itmax=10,mstep,mtau,ia,itau
      REAL(r8) :: psi,rho,drho,u1max,eps=1e-10,theta,r
      REAL(r8) :: a_in,r0_in,b0_in,beta0_in,p_pres_in,p_sig_in,tol_in

      NAMELIST/gse_input/a_in,r0_in,b0_in,q0,beta0_in,p_pres_in,
     $     p_sig_in,tol_in,ma,mtau
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,NML=gse_input)
      CALL ascii_close(in_unit)
      a=a_in
      r0=r0_in
      b0=b0_in
      beta0=beta0_in
      p_pres=p_pres_in
      p_sig=p_sig_in
      tol=tol_in
c-----------------------------------------------------------------------
c     allocate spline and bicube types.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
c-----------------------------------------------------------------------
c     iterate to find chi1.
c-----------------------------------------------------------------------
      chi1=a*a/(2*q0)
      it=0
      DO
         CALL gse_int(u1max,mstep,.FALSE.,eq)
         IF((u1max-1._r8) < eps .OR. it > itmax)EXIT
         chi1=chi1/u1max
         it=it+1
      ENDDO
      IF(it > itmax)CALL program_stop
     $     ("Termination by gse_run: unable to find chi1.")
c-----------------------------------------------------------------------
c     compute and store solution, fit to cubice splines.
c-----------------------------------------------------------------------
      CALL spline_alloc(eq,mstep,4)
      CALL gse_int(u1max,mstep,.TRUE.,eq)
      CALL spline_fit(eq,"extrap")
c-----------------------------------------------------------------------
c     diagnose eq.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         eq%name="eq"
         eq%title=(/" psi  ","   f  ","   p  ","   q  ","  r   "/)
         OPEN(UNIT=debug_unit,FILE="gse.out",STATUS="UNKNOWN")
         OPEN(UNIT=bin_unit,FILE="gse.bin",STATUS="UNKNOWN",
     $        FORM="UNFORMATTED")
         CALL spline_write(eq,.TRUE.,.TRUE.,debug_unit,bin_unit,.FALSE.)
         CLOSE(UNIT=bin_unit)
         CLOSE(UNIT=debug_unit)
      ENDIF
c-----------------------------------------------------------------------
c     fill scalars.
c-----------------------------------------------------------------------
      ro=r0
      zo=0
      psio=chi1*b0
c-----------------------------------------------------------------------
c     fill surface quantities and 2D arrays.
c-----------------------------------------------------------------------
      drho=SQRT(psihigh)/(ma+1)
      rho=0
      DO ia=0,ma
         rho=rho+drho
         psi=rho**2
         CALL spline_eval(eq,psi,0)
         sq_in%xs(ia)=psi
         sq_in%fs(ia,1:3)=eq%f(1:3)
         sq_in%fs(ia,1)=sq_in%fs(ia,1)*b0
         sq_in%fs(ia,2)=sq_in%fs(ia,2)*b0**2
         r=eq%f(4)
         DO itau=0,mtau
            theta=twopi*itau/mtau
            rz_in%fs(ia,itau,1)=ro+r*COS(theta)
            rz_in%fs(ia,itau,2)=zo+r*SIN(theta)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate local arrays and process equilibrium.
c-----------------------------------------------------------------------
      CALL spline_dealloc(eq)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gse_run
c-----------------------------------------------------------------------
c     subprogram 2. gse_int.
c     integrates cylindrical Grad-Shafranov equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gse_int(u1max,mstep,final,eq)

      REAL(r8), INTENT(OUT) :: u1max
      INTEGER :: mstep
      LOGICAL, INTENT(IN) :: final
      TYPE(spline_type), INTENT(OUT) :: eq

      INTEGER :: iopt,istate,itask,itol,jac,mf
      INTEGER :: istep
      REAL(r8) :: psi,psiout,atol,rtol,f,p,q,r,pfac,p1fac,jfac,rho

      INTEGER, PARAMETER :: neq=3,liw=20,lrw=22+neq*16,nstep=1000
      INTEGER, DIMENSION(liw) :: iwork
      REAL(r8), DIMENSION(lrw) :: rwork
      REAL(r8), DIMENSION(neq) :: u
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      rho=SQRT(psihigh)/(ma+1)
      psi=rho**2
      psiout=1
      u(1)=psi*4*chi1/a**2*q0/2
      u(2)=psi
      u(3)=r0
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istate=1
      itask=5
      mf=10
      iopt=0
      itol=1
      rtol=tol
      atol=0
      iwork=0
      rwork=0
      rwork(1)=1
      istep=0
c-----------------------------------------------------------------------
c     store.
c-----------------------------------------------------------------------
      DO
         IF(final)THEN
            CALL gse_prof(psi,u,pfac,p1fac,jfac)
            eq%xs(istep)=psi
            f=u(3)
            p=pfac*beta0/2
            q=(a*a*f/(2*r0*chi1))*(u(1)/u(2))
            r=SQRT(u(1))*a
            eq%fs(istep,:)=(/f,p,q,r/)
         ENDIF
c-----------------------------------------------------------------------
c     advance.
c-----------------------------------------------------------------------
         IF(psi >= 1 .OR. istep > nstep .OR. istate < 0)EXIT
         istep=istep+1
         CALL lsode(gse_der,neq,u,psi,psiout,itol,rtol,atol,itask,
     $        istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     compute output values.
c-----------------------------------------------------------------------
      u1max=u(1)
      mstep=istep
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gse_int
c-----------------------------------------------------------------------
c     subprogram 3. gse_der.
c     contains newcomb's differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gse_der(neq,psi,u,du)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), DIMENSION(neq), INTENT(IN) :: u
      REAL(r8), DIMENSION(neq), INTENT(OUT) :: du

      REAL(r8) :: pfac,p1fac,jfac,det
      REAL(r8), DIMENSION(2) :: rhs
      REAL(r8), DIMENSION(2,2) :: mat
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      du(1)=u(1)/u(2)
      CALL gse_prof(psi,u,pfac,p1fac,jfac)
      mat(1,1)=(2*chi1/a)**2/du(1)
      mat(1,2)=u(3)
      mat(2,1)=u(3)
      mat(2,2)=-u(2)
      det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      rhs(1)=-r0**2*p1fac*beta0/2
      rhs(2)=(a/2)**2*du(1)*u(3)/chi1*2/q0*jfac
      du(2)=(mat(2,2)*rhs(1)-mat(1,2)*rhs(2))/det
      du(3)=(mat(1,1)*rhs(2)-mat(2,1)*rhs(1))/det
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gse_der
c-----------------------------------------------------------------------
c     subprogram 4. gse_prof.
c     profiles for cylindrical Grad-Shafranov equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gse_prof(psi,u,pfac,p1fac,jfac)

      REAL(r8), INTENT(IN) :: psi
      REAL(r8), DIMENSION(:), INTENT(IN) :: u
      REAL(r8), INTENT(OUT) :: pfac,p1fac,jfac

      REAL(r8) :: psifac
c-----------------------------------------------------------------------
c     internal.
c-----------------------------------------------------------------------
      psifac=1-psi*psi
      IF(psifac > 0)THEN
         jfac=(1-psi)**p_sig
         pfac=psifac**p_pres
         p1fac=-2*p_pres*psi*psifac**(p_pres-1)
c-----------------------------------------------------------------------
c     external.
c-----------------------------------------------------------------------
      ELSE
         jfac=0
         p1fac=0
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gse_prof
      END MODULE gse_mod
