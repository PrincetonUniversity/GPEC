c-----------------------------------------------------------------------
c     file lar.f.
c     constructs large-aspect-ratio equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. lar_mod.
c     1. lar_der.
c     2. lar_run.
c     3. lar_coefs.
c-----------------------------------------------------------------------
c     subprogram 0. lar_mod.
c     profile declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE lar_mod
      USE inverse_mod
      IMPLICIT NONE

      INTEGER :: lar_m=2,lar_n=1
      REAL(r8) :: lar_a,lar_r0,p_pres,p_sig
      REAL(r8), PRIVATE :: beta0,p,p00,q,sigma,sigma0
      TYPE(spline_type), PRIVATE :: spl

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. lar_der.
c     contains newcomb's differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lar_der(neq,r,y,dy)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: r
      REAL(r8), DIMENSION(neq), INTENT(IN) :: y
      REAL(r8), DIMENSION(neq), INTENT(OUT) :: dy

      REAL(r8) :: x,xfac,bsq,pp
c-----------------------------------------------------------------------
c     compute equilibrium profiles.
c-----------------------------------------------------------------------
      x=r/lar_a
      xfac=1-x*x
      p=p00*xfac**p_pres
      pp=-2*p_pres*p00*x*xfac**(p_pres-1)
      sigma=sigma0/(1+x**(2*p_sig))**(1+1/p_sig)
      bsq=(y(1)/r)**2+y(2)**2
      q=r**2*y(2)/(lar_r0*y(1))
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1)=-pp/bsq*y(1)+sigma*y(2)*r
      dy(2)=-pp/bsq*y(2)-sigma*y(1)/r
      dy(3)=y(1)*lar_r0/r
      dy(4)=y(5)*(q/r)**2/r
      dy(5)=y(2)*(r/q)**2*(r/lar_r0)*(1-2*(lar_r0*q/y(2))**2*pp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lar_der
c-----------------------------------------------------------------------
c     subprogram 2. lar_run.
c     constructs large-aspect-ratio equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lar_run

      LOGICAL :: zeroth=.FALSE.,out=.FALSE.,bin=.FALSE.
      INTEGER, PARAMETER :: neq=5,liw=20,lrw=22+neq*16,nstep=1000
      INTEGER :: iopt,istate,itask,itol,jac,mf
      INTEGER, DIMENSION(liw) :: iwork
      INTEGER :: i,itau,ia,mtau,ma,istep,mstep
      REAL(r8) :: r,dr,atol,rtol,theta,cosfac,sinfac,rfac,beta0,q0
      REAL(r8), PARAMETER :: rmin=1e-4,tol=1e-6
      REAL(r8), DIMENSION(lrw) :: rwork
      REAL(r8), DIMENSION(neq) :: y,dy
      REAL(r8), DIMENSION(:), POINTER :: r2
      REAL(r8), DIMENSION(:,:), POINTER :: temp

      NAMELIST/lar_input/mtau,ma,lar_r0,lar_a,beta0,q0,p_pres,p_sig,
     $     zeroth,out,bin,lar_m,lar_n

      REAL(r8), DIMENSION(:), POINTER :: xptr,fsptr,fs1ptr
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",6x,"r",10x,"dr",3x,5(5x,"y(",i1,")",2x),
     $     6x,"p",8x,"sigma"/)
 20   FORMAT(i6,1p,9e11.3)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,NML=lar_input)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
      ALLOCATE(r2(0:ma),temp(0:nstep,0:8))
      r2=0
c-----------------------------------------------------------------------
c     initialize constants.
c-----------------------------------------------------------------------
      p00=beta0/2
      p_pres=MAX(p_pres,1.001_r8)
      sigma0=2/(q0*lar_r0)
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      r=rmin*lar_a
      y=0
      y(1)=r**2/(lar_r0*q0)
      y(2)=1
      y(3)=y(1)*lar_r0/2
      CALL lar_der(neq,r,y,dy)
      y(5)=dy(5)*r/4
      y(4)=(q/r)**2*y(5)/2
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
      rwork(1)=lar_a
      istep=0
c-----------------------------------------------------------------------
c     open output files and create formats.
c-----------------------------------------------------------------------
      IF(out)THEN
         CALL ascii_open(lar_out_unit,"lar.out","UNKNOWN")
         WRITE(lar_out_unit,10)(i,i=1,neq)
      ENDIF
      IF(bin)CALL bin_open(lar_bin_unit,"lar.bin","UNKNOWN","REWIND",
     $     "none")
c-----------------------------------------------------------------------
c     store, diagnose, and advance.
c-----------------------------------------------------------------------
      DO
         temp(istep,:)=(/r,y,p,sigma,q/)
         IF(out)WRITE(lar_out_unit,20)istep,r,rwork(11),y,p,sigma
         IF(bin)WRITE(lar_bin_unit)REAL(r,4),REAL(y,4),REAL(p,4),
     $        REAL(sigma/sigma0,4),REAL(q,4),REAL(q0/q,4)
         IF(r >= lar_a .OR. istep > nstep .OR. istate < 0)EXIT
         istep=istep+1
         CALL lsode(lar_der,neq,y,r,lar_a,itol,rtol,atol,itask,
     $        istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(lar_out_unit,10)(i,i=1,neq)
         CALL ascii_close(lar_out_unit)
      ENDIF
      IF(bin)CALL bin_close(lar_bin_unit)
c-----------------------------------------------------------------------
c     fit to cubic splines.
c-----------------------------------------------------------------------
      mstep=istep
      CALL spline_alloc(spl,mstep,8)
      spl%xs=temp(0:mstep,0)
      spl%fs=temp(0:mstep,1:8)
      CALL spline_fit(spl,"extrap")
      xptr => spl%xs
      fsptr => spl%fs(:,1)
      fs1ptr => spl%fs1(:,1)
c-----------------------------------------------------------------------
c     fill scalars.
c-----------------------------------------------------------------------
      ro=lar_r0
      zo=0
      psio=temp(mstep,3)
c-----------------------------------------------------------------------
c     interpolate and fill surface quantities.
c-----------------------------------------------------------------------
      dr=lar_a/(ma+1)
      r=0
      DO ia=0,ma
         r=r+dr
         CALL spline_eval(spl,r,1)
         sq_in%xs(ia)=spl%f(3)/psio
         sq_in%fs(ia,1)=lar_r0*spl%f(2)
         sq_in%fs(ia,2)=spl%f(6)
         sq_in%fs(ia,3)=spl%f(8)
         r2(ia)=-(spl%f(4)*r/spl%f(8))/spl%f1(3)
      ENDDO
      IF(zeroth)r2=0
c-----------------------------------------------------------------------
c     fill 2D arrays.
c-----------------------------------------------------------------------
      DO itau=0,mtau
         theta=twopi*itau/mtau
         cosfac=COS(theta)
         sinfac=SIN(theta)
         DO ia=0,ma
            r=lar_a*(ia+1)/(ma+1)
            rfac=r+r2(ia)*cosfac
            rz_in%fs(ia,itau,1)=ro+rfac*cosfac
            rz_in%fs(ia,itau,2)=zo+rfac*sinfac
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate local arrays and process equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(r2,temp)
c      CALL lar_coefs
      CALL spline_dealloc(spl)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lar_run
c-----------------------------------------------------------------------
c     subprogram 3. lar_coefs.
c     diagnoses cylindrical limit of diagonal matrix elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lar_coefs

      INTEGER :: ir,mr=100
      REAL(r8) :: chi1,chi2,f,g,g11,g22,g33,jtheta,k,p1,psi,q1,r,v1,
     $     singfac,dr
c-----------------------------------------------------------------------
c     diagnose matrix elements.
c-----------------------------------------------------------------------
      CALL bin_open(lar_bin_unit,"lar1.bin","UNKNOWN","REWIND","none")
      dr=lar_a/mr
      DO ir=1,mr
         r=dr*ir
         CALL spline_eval(spl,r,1)
         psi=spl%f(3)/psio
         v1=(twopi)**2*lar_r0*r
         p=spl%f(6)
         p1=spl%f1(6)/v1
         q=spl%f(8)
         q1=spl%f1(8)/v1
         chi1=twopi*spl%f(1)*lar_r0/(r*v1)
         chi2=chi1*(spl%f1(1)/spl%f(1)-2/r)
         jtheta=-twopi*lar_r0*spl%f1(2)/v1
         singfac=lar_m-lar_n*q
         g11=1/(twopi**2*lar_r0*r)**2
         g22=(twopi*r)**2
         g33=(twopi*lar_r0)**2
         f=(chi1*singfac)**2*g22*g33/(lar_m**2*g33+lar_n**2*g22)
         k=chi1*chi2*g22+q*chi1*(q*chi2+q1*chi1)+2*p1
     $        -chi1*(lar_n*g22+lar_m*q*g33)/(lar_m**2*g33+lar_n**2*g22)
     $        *(lar_n*chi2*g22+lar_m*(q*chi2+q1*chi1)*g33
     $        -2*(singfac*jtheta+lar_n*p1/chi1))
         g=(twopi*chi1*singfac)**2*g11
     $        +chi2*g22+(q*chi2+q1*chi1)**2*g33
     $        +p1*chi2/chi1+jtheta*q1*chi1
     $        -(lar_n*chi2*g22+lar_m*(q*chi2+q1*chi1)*g33
     $        -2*(singfac*jtheta+lar_n*p1/chi1))**2
     $        /(lar_m**2*g33+lar_n**2*g22)
         WRITE(lar_bin_unit)REAL(psi,4),REAL(SQRT(psi),4),REAL(q,4),
     $        REAL(g11,4),REAL(g22,4),REAL(g33,4),
     $        REAL(f,4),REAL(g,4),REAL(k,4)
      ENDDO
      CALL bin_close(lar_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lar_coefs
      END MODULE lar_mod
