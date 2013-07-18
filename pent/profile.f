c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     calculate ntv torque and kinetic perturbed energy
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. natorq_mod
c     a. nu
c     b. xcoeffs
c     c. xintegral
c     1. natorq_gen
c     2. natorq_lar
c-----------------------------------------------------------------------
c     subprogram 0. natorq_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE profile_mod
      USE pentio_mod      
      USE grid_mod
      USE energy_mod

      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     Function C. xintegral
c     Energy integral spline formation
c-----------------------------------------------------------------------
      FUNCTION xintegral(wn,wt,we,wd,wb,l,nutype,s,ximag,xmax,offset)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      type(cspline_type) :: xintegral
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax
      INTEGER, INTENT(IN) :: s
      CHARACTER(16), INTENT(IN) :: nutype
      LOGICAL, INTENT(IN) :: offset

      REAL(r8) :: quadb,quadc,xmin,ccntr
      REAL(r8), DIMENSION(1:2) :: sqrtxres,xres
      INTEGER  :: i,nres,neo
      REAL(r8), DIMENSION(0:nx) :: nux
      REAL(r8), DIMENSION(0:1,0:nx) :: xdx
      COMPLEX(r8), DIMENSION(0:nx) :: cx,fx,gx

      neo = 0
      ccntr = 0
      IF(neorot_flag) neo=1
      IF(nutype=='zero') ccntr = ximag
c-----------------------------------------------------------------------
c     Solve quadratic eq to find bounce harmonic resonances
c-----------------------------------------------------------------------
      quadb =-l*wb/(nn*wd)
      quadc = we/wd
      sqrtxres = 0
      nres = 1  ! Default at zero if no resonances
      IF(quadb**2-4*quadc>=0)THEN
         sqrtxres = (/ 0.5*(-quadb+SQRT(quadb**2-4*quadc)), 
     $        0.5*(-quadb-SQRT(quadb**2-4*quadc)) /)
         DO i=1,2
            IF(sqrtxres(i)<0 .OR. sqrtxres(i)>SQRT(xmax)) sqrtxres(i)=0
         ENDDO
         IF(sqrtxres(1)==sqrtxres(2)) sqrtxres(2)=0
         IF(sqrtxres(2)>0)THEN
            sqrtxres = (/ sqrtxres(2),sqrtxres(1) /)
            IF(sqrtxres(2)>0) nres = 2
         ENDIF
      ENDIF

c-----------------------------------------------------------------------
c     Virtual grid concentrated at bounce harmonic resonances
c-----------------------------------------------------------------------
      CALL cspline_alloc(xintegral,nx,1)
      xintegral%xs(:) = xnorm
      xres = sqrtxres**2
      xmin = 0.0

      IF(xres(1)==0)THEN
         CALL powergrid(xnorm,xdx,(/xmin,xmax/),5,"lower")
      ELSE
         CALL powergrid(xnorm,xdx,(/xmin,xres(1:nres),xmax/),5,"none")
      ENDIF
      cx = xdx(0,:)+ifac*ccntr
      nux = nu_all(nutype,s,l,xdx(0,:))
      fx = -(we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
      IF(offset) fx = -(wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
      gx = ifac*(l*wb*SQRT(cx)-nn*(we+wd*cx))-nux
      xintegral%fs(:,1) = xdx(1,:)*fx/gx

      IF(kolim_flag)THEN
         !ignore all the above complications
         CALL powergrid(xnorm,xdx,(/xmin,xmax/),3,"lower")
         cx = xdx(0,:)!+ifac*ccntr
         xintegral%fs(:,1) = xdx(1,:)*cx**2.5*EXP(-cx)
      ENDIF

      CALL cspline_fit(xintegral,"extrap")
      CALL cspline_int(xintegral)

c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION xintegral



c-----------------------------------------------------------------------
c     function d. spl_roots.
c     finds roots of a spline.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spl_roots(roots,spl,iqty)
      INTEGER, INTENT(IN) :: iqty
      TYPE(spline_type) :: spl
      
      REAL(r8), DIMENSION(:), ALLOCATABLE :: roots
      CHARACTER(80) :: message
      INTEGER :: iroot,it,ix,nroots
      INTEGER, PARAMETER :: itmax=20
      REAL(r8), PARAMETER :: eps=1e-10
      REAL(r8) :: x,dx,lx,lf,f,df

c-----------------------------------------------------------------------
c     compute length.
c-----------------------------------------------------------------------
      lx=spl%xs(spl%mx)-spl%xs(0)
      iroot = 1
      lf = MAXVAL(spl%fs(:,iqty))-MINVAL(spl%fs(:,iqty))
      nroots = 0
      DO ix=1,spl%mx
         IF(spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .LE. 0.0) nroots=nroots+1
         !zeros counted twice
         IF(spl%fs(ix,iqty)==0.0 .AND. ix<spl%mx) nroots=nroots-1 
      ENDDO
      IF(spl%fs(0,iqty)==0) nroots=nroots-1
      IF(spl%fs(spl%mx,iqty)==0) nroots=nroots-1
      ALLOCATE(roots(0:nroots+1))
      roots(0) = spl%xs(0)
      roots(nroots+1) = spl%xs(spl%mx)
c-----------------------------------------------------------------------
c     find all zero passings, intialize at larger gradient.
c-----------------------------------------------------------------------
      DO ix=1,spl%mx
         IF(spl%fs(ix,iqty)==0.0 .AND. ix<spl%mx) CYCLE !don't do zeros twice
         IF (spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .LE. 0.0) THEN
            x=spl%xs(ix-1)-spl%fs(ix-1,iqty)*(spl%xs(ix)
     $           -spl%xs(ix-1))/(spl%fs(ix,iqty)-spl%fs(ix-1,iqty))
            f=HUGE(f)
            dx=lx
            it=0
c-----------------------------------------------------------------------
c     locate roots by newton iteration.
c-----------------------------------------------------------------------
            DO
               CALL spline_eval(spl,x,1)
               df=spl%f(iqty)-f
               IF(ABS(dx) < eps*lx 
     $                 .OR. ABS(df) < eps*lf 
     $              .OR. it >= itmax)EXIT
               it=it+1
               f=spl%f(iqty)
               dx=-spl%f(iqty)/spl%f1(iqty)
               x=x+dx
            ENDDO
c-----------------------------------------------------------------------
c     abort on failure.
c-----------------------------------------------------------------------
            IF(it >= itmax)THEN
               !x=spl%xs(ix)
               !IF(spl%fs(ix,iqty) .LT. 0.0) x = spl%xs(ix-1)
               WRITE(*,*) "WARNING: roots convergence failure!"
               WRITE(*,*) "-- failure accured at index ",ix
               WRITE(*,*) "-- indexed x value ",spl%xs(ix)
               WRITE(*,*) "-- estimated root ",x
            ENDIF
c-----------------------------------------------------------------------
c     store each root, and allocate it to the spline.
c-----------------------------------------------------------------------
            roots(iroot) = x
            iroot = iroot+1
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spl_roots












c-----------------------------------------------------------------------
c     subprogram 1. profile_gen.
c     Calculate NTV on a magnetic surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE profile_gen(lagbpar,divxprp,nl,maskpsi,ptfac,ximag,
     $     xmax,wefac,wdfac,collision,outsurfs)

      INTEGER, INTENT(IN) :: nl,maskpsi
      REAL(r8), INTENT(IN) :: ptfac,ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10), INTENT(IN) :: outsurfs
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf), INTENT(IN) :: 
     $     lagbpar,divxprp
      CHARACTER(16) :: collision

      LOGICAL :: outsurf
      INTEGER :: istep,itheta,ibmin,i,ntheta,ilmda,ll,sigma,
     $     iters,iter,indx,s,iout,ix,pass,xsigma
      REAL(r8):: chrg,mass,epsa,epsr,lmda,welec,wdian,wdiat,wbhat,wdhat,
     $     wbbar,wdbar,djdjbar,bo,ls,dpsi,bmax,bmin,dt
      REAL(r8), DIMENSION(:), ALLOCATABLE :: bpts,wd0
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: htheta,tdt,ldl
      COMPLEX(r8) :: totals
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: pl

      REAL(r8), DIMENSION(100,5,10,1,8) :: x_out
      REAL(r8), DIMENSION(1000,5,10,1,15) :: b_out
      REAL(r8), DIMENSION(1,1000,10,1:2*nl+1,16) :: p_out
      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
      CHARACTER(32) :: outfile
      CHARACTER(512) :: outtext,outlbls
      
      TYPE(spline_type) :: omega,tspl,vpar,lspl
      TYPE(cspline_type) o1spl,action,xspl,ntvdl,offdl,rot,rotp

c-----------------------------------------------------------------------
c     Set species.
c-----------------------------------------------------------------------
      chrg = -1*echarge*icharge
      IF(chrg >= 0)THEN
         s = 1
         mass = pmass*imass
      ELSE IF(icharge == -1)THEN
         s = 2
         mass = emass
      ELSE
         CALL pentio_stop('Negative charge must be electrons (-1).')
      ENDIF
      pass = 0
      xsigma=0
      IF(passing_flag) pass=1
      IF(divxprp_flag) xsigma=1
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------  
      PRINT *, "PENT - general calculation v1.0"

      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
      x_out = 0
      b_out = 0
      p_out = 0

      epsa = SQRT(SUM(rzphi%fs(rzphi%mx,:,1))/rzphi%my)/ro

      iters = (mstep-1)/maskpsi +1
      CALL cspline_alloc(ntv,iters-1,2*nl+1)
      CALL cspline_alloc(ntvp,iters-1,2*nl+1)
      CALL cspline_alloc(rot,iters-1,2*nl+1+1)
      CALL cspline_alloc(rotp,iters-1,2*nl+1+1)

      CALL spline_alloc(tspl,mthsurf,5)
      CALL cspline_alloc(o1spl,mthsurf,2)
      CALL spline_alloc(vpar,mthsurf,1)
      tspl%xs(:) = theta(:)
      o1spl%xs(:) = theta(:)
      vpar%xs(:) = theta(:)

      CALL spline_alloc(omega,nt,2)
      CALL cspline_alloc(action,nt,3)
      ALLOCATE(htheta(0:nt,2),pl(0:nt),tdt(0:1,0:nt))
      omega%xs(:) = tnorm(:)
      action%xs(:) = tnorm(:)

      CALL spline_alloc(lspl,nlmda,2*nl+1+3)
      CALL cspline_alloc(ntvdl,nlmda,2*nl+1)
      ALLOCATE(ldl(0:1,0:nlmda))

      DO istep=1,mstep,maskpsi
         iter = (istep-1)/maskpsi+1
         ! print details of specified surfaces
         outsurf = .FALSE.
         IF (MOD(iter,iters/10)==0)THEN 
            WRITE(*,'(1x,a9,i3,a19)') "volume = ",iter/(iters/10)*10,
     $           " % NTV computations"
         ENDIF
         IF(istep==1)THEN !First point
            DO i=1,10
               dpsi = ABS(psifac(istep)-outsurfs(i))
               IF(dpsi<=ABS(psifac(istep+maskpsi)-outsurfs(i)))THEN 
                  outsurf = .TRUE.
                  iout = i
               ENDIF
            ENDDO
         ELSEIF(istep+maskpsi>mstep)THEN ! last point
            DO i=1,10
               dpsi = ABS(psifac(istep)-outsurfs(i))
               IF(dpsi<=ABS(psifac(istep-maskpsi)-outsurfs(i)))THEN
                  outsurf = .TRUE.
                  iout = i
               ENDIF
            ENDDO
         ELSE
            DO i = 1,10
               dpsi = ABS(psifac(istep)-outsurfs(i))
               IF(dpsi<=ABS(psifac(istep-maskpsi)-outsurfs(i)).AND.
     $               dpsi<=ABS(psifac(istep+maskpsi)-outsurfs(i)))THEN 
                  outsurf = .TRUE.
                  iout = i
               ENDIF
            ENDDO            
         ENDIF

         ntv%xs(iter-1) = psifac(istep)
         ntvp%xs(iter-1)= psifac(istep)
         rot%xs(iter-1) = psifac(istep)
         rotp%xs(iter-1)= psifac(istep)
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL spline_eval(kin,psifac(istep),1)
         epsr = epsa*SQRT(psifac(istep))
         !Rotation frequencies: electric and density and temp gradient
         welec = wefac*kin%f(5)
         wdian =-twopi*kin%f(s+2)*kin%f1(s)/(chrg*chi1*kin%f(s))
         wdiat =-twopi*kin%f1(s+2)/(chrg*chi1)
         rot%fs(iter-1,1) = welec+wdian+wdiat
         rotp%fs(iter-1,1) = welec+wdian+wdiat
         
         !Poloidal functions - note ABS(A*clebsch) = ABS(A)
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            tspl%fs(itheta,1)= eqfun%f(1)            !b
            tspl%fs(itheta,2)= eqfun%fx(1)/chi1      !db/dpsi
            tspl%fs(itheta,3)= eqfun%fy(1)           !db/dtheta
            tspl%fs(itheta,4)= rzphi%f(4)/chi1       !jac
            tspl%fs(itheta,5)= rzphi%fx(4)/chi1**2   !dj/dpsi
            o1spl%fs(itheta,1)= lagbpar(istep,itheta) !lagb
            o1spl%fs(itheta,2)= divxprp(istep,itheta) !divx
         ENDDO
         ibmin = MINLOC(tspl%fs(:,1),DIM=1)-1 !not used anymore?
         ! clebsch conversion * now in djdt
         !o1spl%fs(:,1)=o1spl%fs(:,1)*EXP(-twopi*ifac*nn*sq%f(4)*
     $   !        (theta(:)-theta(ibmin)))
         !o1spl%fs(:,2)=o1spl%fs(:,2)*EXP(-twopi*ifac*nn*sq%f(4)*
     $   !        (theta(:)-theta(ibmin)))
         CALL spline_fit(tspl,"periodic")
         CALL cspline_fit(o1spl,"periodic")
         bmax = MAXVAL(tspl%fs(:,1),DIM=1)
         bmin = MINVAL(tspl%fs(:,1),DIM=1) 
         DO i=2,8 !8th smallest so spline has more than 1 pt
            bmin = MINVAL(tspl%fs(:,1),MASK=tspl%fs(:,1)>bmin,DIM=1)
         ENDDO         
         
c-----------------------------------------------------------------------
c     Energy integration for trapped and pasing particles
c-----------------------------------------------------------------------
         DO sigma = 0,pass
            ldl = powspace(sigma*1e-6/bmax+
     $           (1-sigma)*(1+(1/bmin-1/bmax)*ptfac)/bmax
     $           ,sigma*(1-ptfac)/bmax+(1-sigma)/bmin,1,nlmda+1,"both")
            DO ilmda=0,nlmda
               lmda = ldl(0,ilmda)
               vpar%fs(:,1) = 1-lmda*tspl%fs(:,1)
               CALL spline_fit(vpar,"extrap")
               CALL spl_roots(bpts,vpar,1)
               CALL spline_eval(vpar,SUM(bpts(1:2))/2,0)
               ! the usual case with bpts centered around 0 requires shift
               IF(vpar%f(1)<=0.0 .AND. sigma==0) 
     $              bpts(1:2)=(/bpts(2)-1.0,bpts(1)/)
c               print *,"lmda_n=",lmda/ldl(0,nlmda),"bpts = ",
c     $              bpts(1-sigma:2-sigma)
               ! hack - shouldn't happen (had repeated zeros in spl_roots)
               IF(bpts(1)>=bpts(2) .AND. sigma==0)THEN
                  lspl%fs(ilmda,:) = lspl%fs(ilmda-1,:)
                  lspl%xs(ilmda) = lmda
                  PRINT *, "WARNING: Bounce point crossing"
                  PRINT *, "Normalized (psi,lmda) = ",psifac(istep),
     $                 lmda/ldl(0,nlmda)
                  PRINT *,"Bounce pts = ",bpts
                  DEALLOCATE(bpts)
                  CYCLE
               ENDIF
               CALL powergrid(tnorm,tdt,bpts(1-sigma:2-sigma),
     $              4,"both")
               DEALLOCATE(bpts)
c-----------------------------------------------------------------------
c     Bounce averages. Bounce & precession freqs, variation in action
c-----------------------------------------------------------------------
               omega%fs(0,:) = 0
               omega%fs(omega%mx,:) = 0.0
               action%fs(0,:) = 0
               action%fs(action%mx,:) = 0.0
               DO itheta=1,nt-1
                  dt = tdt(1,itheta)
                  CALL spline_eval(tspl,MOD(tdt(0,itheta)+1,1.0),0)
                  CALL cspline_eval(o1spl,MOD(tdt(0,itheta)+1,1.0),0)
                  CALL spline_eval(vpar,MOD(tdt(0,itheta)+1,1.0),0)
                  omega%fs(itheta,1) = dt*tspl%f(4)*tspl%f(1)
     $                 /SQRT(vpar%f(1))
                  omega%fs(itheta,2) = dt*tspl%f(4)*tspl%f(2)*
     $                 (1-1.5*lmda*tspl%f(1))/SQRT(vpar%f(1))+
     $                 dt*tspl%f(5)*tspl%f(1)*SQRT(vpar%f(1))
                  action%fs(itheta,1) = dt*tspl%f(4)*tspl%f(1)*(
     $                 xsigma*o1spl%f(2)*SQRT(vpar%f(1))+
     $                 o1spl%f(1)*(1-1.5*lmda*tspl%f(1))/
     $                 SQRT(vpar%f(1)))
     $                 *EXP(-twopi*ifac*nn*sq%f(4)*(tdt(0,itheta)-
     $                 tdt(0,1))) ! doesn't matter since |dj|^2
               ENDDO
               call spline_fit(omega,"extrap")
               CALL spline_int(omega)
               htheta(:,1)= omega%fsi(:,1)
               htheta(:,2)= omega%fsi(:,2)
 !!!!!              !missing some analytic piece of h(:,2)?
               htheta(:,1)=htheta(:,1)/((2-sigma)*omega%fsi(omega%mx,1))
               htheta(:,2)=htheta(:,2)/((2-sigma)*omega%fsi(omega%mx,2))
               !Bounce averaged frequencies without energy dependencies
               wbbar = twopi/((2-sigma)*omega%fsi(omega%mx,1))
               wdbar = wdfac*wbbar*2*(2-sigma)*omega%fsi(omega%mx,2)
               lspl%xs(ilmda)   = lmda
               lspl%fs(ilmda,1) = wbbar
               lspl%fs(ilmda,2) = wdbar
               wbhat = SQRT(2*kin%f(s+2)/mass)*wbbar
               wdhat = wdbar*kin%f(s+2)/chrg
               DO ll=-nl,nl
                  ! added in MISK benchmark
                  ls = REAL(ll+sigma*nn*sq%f(4))
                  !SHOULD HAVE MORE SUFFISTICATED PHASE FACTOR HERE
                  IF(welec<wdhat*3)THEN !mag precession dominates drift
                     pl(:) = EXP(-twopi*ifac*ls*htheta(:,2))
                  ELSE !electric precession dominates drift (normal)
                     pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
                  ENDIF
                  pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
                  action%fs(:,2) =action%fs(:,1)*(pl(:)+1/pl(:))
                  action%fs(:,3) =omega%fs(:,1)*(1/pl(:)+pl(:))**2.0
                  CALL cspline_fit(action,"extrap")
                  CALL cspline_int(action)
                  djdjbar=REAL(action%fsi(action%mx,2))*
     $                 REAL(action%fsi(action%mx,2))+
     $                 AIMAG(action%fsi(action%mx,2))*
     $                 AIMAG(action%fsi(action%mx,2))
                  lspl%fs(ilmda,4+ll+nl) = djdjbar
c-----------------------------------------------------------------------
c     Record a sampling of intermediate steps for debuging
c-----------------------------------------------------------------------
                  IF(outsurf)THEN
                     !iout = MAX(MIN(iter/(iters/10),10),1)
                     p_out(1,ilmda+1,iout,ll+nl+1,14:16) = (/
     $                    wbbar*REAL(action%fsi(action%mx,2)),
     $                    wbbar*AIMAG(action%fsi(action%mx,2)),lmda /)
                     IF(ll==0 .AND. (MOD(ilmda,nlmda/4)==0))THEN
                        DO i = 0,omega%mx
                           b_out(i+1,ilmda/(nlmda/4)+1,iout,1,:) = (/ 
     $                          psifac(istep),lmda,omega%xs(i),tdt(:,i),
     $                          REAL(action%fs(i,2)),
     $                          AIMAG(action%fs(i,2)),
     $                          REAL(action%fsi(i,2)),
     $                          AIMAG(action%fsi(i,2))
     $                          ,omega%fs(i,:),omega%fsi(i,:),
     $                          htheta(i,:)/)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO ! l loop
            ENDDO ! lmda loop
            
c-----------------------------------------------------------------------
c     Pitch angle integrand consentrated near magnetic precession nulls
c-----------------------------------------------------------------------
            ! watch out for resonance in lmda (esp. low nu, low welec)
            ! assume peak is when BHR in energy integrand peak, which is
            ! 1.5<~x<~2.5 for nu=0. Not exact -> low power concentration
            lspl%fs(:,3) = 2.5*lspl%fs(:,2)+welec/(kin%f(s+2)/chrg)
            CALL spline_fit(lspl,'extrap')
            CALL spl_roots(wd0,lspl,3)
            CALL powergrid(lnorm,ldl,wd0,1,"both") !includes ends
            DEALLOCATE(wd0)
            DO ilmda=0,nlmda
               lmda = ldl(0,ilmda)
               CALL spline_eval(lspl,lmda,0)
               wbhat = SQRT(2*kin%f(s+2)/mass)*lspl%f(1)
               wdhat = (kin%f(s+2)/chrg)*lspl%f(2)
               DO ll=-nl,nl
                  ls = REAL(ll+sigma*nn*sq%f(4)) 
                  IF(offset_flag)THEN
                     xspl = xintegral(wdian,wdiat,welec,wdhat,wbhat,ls,
     $                    collision,s,ximag,xmax,.TRUE.)
                     offdl%fs(ilmda,nl+ll+1)=ldl(1,ilmda)*lspl%f(1)*
     $                 lspl%f(4+nl+ll)*xspl%fsi(xspl%mx,1)
                  ENDIF
                  xspl = xintegral(wdian,wdiat,welec,wdhat,wbhat,ls,
     $                 collision,s,ximag,xmax,.FALSE.)
                  ntvdl%fs(ilmda,nl+ll+1) = !ldl(1,ilmda)*
     $                 lspl%f(1)*lspl%f(4+nl+ll)*xspl%fsi(xspl%mx,1)
                  ntvdl%xs(ilmda) = lmda!lnorm(ilmda)
c-----------------------------------------------------------------------
c     Record a sampling of intermediate steps for debuging
c-----------------------------------------------------------------------
                  IF(outsurf)THEN
                     !iout = MAX(MIN(iter/(iters/10),10),1)
                     p_out(1,ilmda+1,iout,ll+nl+1,1:11) = 
     $                    (/psifac(istep),lnorm(ilmda),ldl(:,ilmda),
     $                    REAL(xspl%fsi(nx,1)),
     $                    AIMAG(xspl%fsi(nx,1)),
     $                    REAL(ntvdl%fs(ilmda,nl+ll+1)),
     $                    AIMAG(ntvdl%fs(ilmda,nl+ll+1)),lspl%f(1)*ro,
     $                    ro**2*bo*(lspl%f(2)),
     $                    lspl%f(4+nl+ll)*lspl%f(1)**2/ro**2/)
                     IF(ll==0 .AND. (MOD(ilmda,nlmda/4)==0))THEN
                        DO i = 0,nx,nx/99
                           ix = MIN(i/(nx/99)+1,100)
                           x_out(ix,ilmda/(nlmda/4)+1,iout,1,:)=(/ 
     $                          psifac(istep),lmda,xspl%xs(i),
     $                          REAL(xspl%fs(i,1)),AIMAG(xspl%fs(i,1)),
     $                          REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1))
     $                          ,lspl%f(1)*lspl%f(4+nl+ll) /)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO ! l loop
            ENDDO ! lmda loop
c-----------------------------------------------------------------------
c     Final integration over Lambda and normalization to Nm
c-----------------------------------------------------------------------
            CALL cspline_fit(ntvdl,"extrap")
            CALL cspline_int(ntvdl)
            
            IF(outsurf)THEN
               DO ll=-nl,nl
                  DO i=0,lspl%mx
                     p_out(1,i+1,iout,ll+nl+1,12:13) = 
     $                    (/ REAL(ntvdl%fsi(i,ll+nl+1)),
     $                    AIMAG(ntvdl%fsi(i,ll+nl+1)) /)
                  ENDDO
               ENDDO
            ENDIF

            IF(sigma==0) ntv%fs(iter-1,:)=ntvdl%fsi(ntvdl%mx,:)*
     $           ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
            IF(sigma==1) ntvp%fs(iter-1,:)=ntvdl%fsi(ntvdl%mx,:)*
     $           ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
            IF(offset_flag)THEN
               CALL cspline_fit(offdl,"extrap")
               CALL cspline_int(offdl)
               IF(sigma==0)THEN
                  rot%fs(iter-1,2:) = rot%fs(iter-1,1)*
     $                 offdl%fsi(offdl%mx,:)/(ntvdl%fsi(ntvdl%mx,:)-
     $                 offdl%fsi(offdl%mx,:))
               ELSEIF(sigma==1)THEN
                  rotp%fs(iter-1,2:) = rotp%fs(iter-1,1)*
     $                 offdl%fsi(offdl%mx,:)/(ntvdl%fsi(ntvdl%mx,:)-
     $                 offdl%fsi(offdl%mx,:))
               ENDIF
            ENDIF
         ENDDO ! sigma loop
      ENDDO ! psi loop
c-----------------------------------------------------------------------
c     Final integration over normalized flux
c-----------------------------------------------------------------------
      !Torque in Nm, Energy*R0 in Jm
      CALL cspline_fit(ntv,"extrap")
      CALL cspline_int(ntv)
      totals = SUM(ntv%fsi(ntv%mx,:))
      PRINT *, "Total NTV torque = ",REAL(totals)
      PRINT *, "Total Kinetic Energy = ",AIMAG(totals)/(2*nn)
      PRINT *, "alpha/s  = ",REAL(totals)/(-1*AIMAG(totals))

      CALL spline_dealloc(tspl)
      CALL cspline_dealloc(o1spl)
      CALL spline_dealloc(vpar)

      CALL spline_dealloc(omega)
      CALL cspline_dealloc(action)

      CALL spline_dealloc(lspl)
      CALL cspline_dealloc(ntvdl)

c-----------------------------------------------------------------------
c     Write output files
c-----------------------------------------------------------------------
      outfile ="pent_n"//TRIM(sn)
      outtext ="Trapped particle torque (Nm) and energy (J)"
      IF(offset_flag)THEN
         CALL cspline_fit(rot,"extrap")
         CALL pentio_write_profile(outfile,outtext,ntv,rot)
      ELSE
         CALL pentio_write_profile(outfile,outtext,ntv)
      ENDIF
      CALL cspline_dealloc(rot)

      IF(passing_flag)THEN
         outfile ="pent_passing_n"//TRIM(sn)
         outtext ="Passing particle torque (Nm) and energy (J)"
         CALL cspline_fit(ntvp,"extrap")
         CALL cspline_int(ntvp)
         IF(offset_flag)THEN
            CALL cspline_fit(rotp,"extrap")
            CALL pentio_write_profile(outfile,outtext,ntvp,rotp)
         ELSE
            CALL pentio_write_profile(outfile,outtext,ntvp)
         ENDIF
      ENDIF
      CALL cspline_dealloc(rotp)

      IF(pitch_flag)THEN
         outfile = "pent_pitch_n"//TRIM(sn)
         outtext = "Pitch angle integrand and functions"
     $        //newl//"Normalizations are as follows:"//newl//
     $        "Lambda: m*v_perp^2/(2B)"//newl//
     $        "omega_D: omega_D*B0*R0^2/(E/Ze)"
     $        //newl//"deltaJ^2: (dJ_l*dJ_-l)*/R0^2"//newl//
     $        "T_phi,2ndeltaW: T_phi,2ndeltaW/(-2n^2TN*chi'/sqrt(pi))"
     $        //newl//"Ix = unitless"//newl//"H/x: Lagrangian/(E/T)"
         WRITE(outlbls,'(16(1x,a16))') "psi_n","Lambda_n","Lambda",
     $        "dLambda_n","real(I_x)",
     $        "imag(I_x)","T_phi","2ndeltaW","omega_b","omega_D",
     $        "deltaJ^2","int(T_phi)","int(2ndeltaW)",
     $        "real(H/x)","imag(H/x)","Lambda_H"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(p_out),p_out,6)
      ENDIF

      IF(energy_flag)THEN
         outfile = "pent_energy_n"//TRIM(sn)
         outtext = "Energy integral samples"//newl//
     $        "Normalizations are as follows:"//newl//
     $        "Lambda: m*v_perp^2/(2B)"//newl//"x: E/T"//newl//
     $        "omega_b*deltaJ^2: wb*dJ^2/R0"//newl//
     $        "omega_b*deltaJ^2*T_phi: unitless"//newl//
     $        "T_phi: T_phi/(-2n^2TN*chi'/sqrt(pi))"//newl//
     $        "2ndeltaW: 2*n*deltaW /(-2n^2TN*chi'/sqrt(pi))"
         WRITE(outlbls,'(8(1x,a16))') "psi","Lambda","x_n","T_phi",
     $        "2ndeltaW","int(T_phi)","int(2ndeltaW)",
     $        "deltaJ^2omegab"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(x_out),x_out,7)
      ENDIF

      IF(bounce_flag)THEN
         outfile = "pent_bounce_n"//TRIM(sn)
         outtext = "Bounce averaging integrands"
         WRITE(outlbls,'(15(1x,a16))') "psi","Lambda","theta_n",
     $        "vartheta","dvartheta/dtheta_n",
     $        "real(deltaJ)","imag(deltaJ)","realint(deltaJ)",
     $        "imagint(deltaJ)","omega_b","omega_D","int(omega_b)",
     $        "int(omega_D)","h_E","h_B"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(b_out),b_out,7)

      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE profile_gen


c-----------------------------------------------------------------------
c     subprogram 5. profile_lar.
c     Calculate NTV using large aspect ratio approximations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE profile_lar(lagbpar,divxprp,nl,maskpsi,ximag,xmax,
     $     wefac,wdfac,collision,outsurfs)
      
      INTEGER, INTENT(IN) :: nl,maskpsi
      REAL(r8), INTENT(IN) :: ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10), INTENT(IN) :: outsurfs
      COMPLEX(r8), DIMENSION(mstep,mpert), INTENT(IN) :: lagbpar,divxprp
      CHARACTER(16) :: collision

      LOGICAL :: outsurf
      INTEGER, PARAMETER :: nk = 50
      INTEGER :: istep,k,i,ipert,jpert,ll,sigma,
     $     pass,iters,iter,indx,s,nsplt,iout,ix
      REAL(r8) :: chrg,mass,kappa,wbhat,wdhat,welec,wdian,wdiat,
     $     epsa,epsr,wt,wg,bo,ls,dpsi
      REAL(r8), DIMENSION(mpert,2):: fm
      COMPLEX(r8) :: totals

      REAL(r8), DIMENSION(1,100,10,1:2*nl+1,6)  :: x_out
      REAL(r8), DIMENSION(1,nk+1,10,1:2*nl+1,4) :: p_out
      REAL(r8), DIMENSION(1,1,mstep,1,10) :: b_out
      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
      CHARACTER(32) :: outfile
      CHARACTER(512) :: outtext,outlbls
      
      TYPE(spline_type) :: kspl
      TYPE(cspline_type):: oxspl,xspl,rot,rotp

c-----------------------------------------------------------------------
c     Set species.
c-----------------------------------------------------------------------
      PRINT *, "PENT - Large Aspect Ratio"

      chrg = -1*echarge*icharge
      IF(chrg >= 0)THEN
         s = 1
         mass = pmass*imass
      ELSE IF(icharge == -1)THEN
         s = 2
         mass = emass
      ELSE
         CALL pentio_stop('Negative charge must be electrons (-1).')
      ENDIF
      pass = 0
      IF(passing_flag) pass=1

c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------  
      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
      x_out = 0
      p_out = 0
      b_out = 0

      epsa = SQRT(SUM(rzphi%fs(rzphi%mx,:,1))/rzphi%my)/ro

      iters = (mstep-1)/maskpsi +1
      CALL cspline_alloc(lar,iters-1,2*nl+1)
      CALL cspline_alloc(larp,iters-1,2*nl+1)
      CALL cspline_alloc(rot,iters-1,2*nl+1+1)
      CALL cspline_alloc(rotp,iters-1,2*nl+1+1)

      DO istep=1,mstep,maskpsi
         iter = (istep-1)/maskpsi+1
         outsurf = .FALSE.
         IF (MOD(iter,iters/10)==0) THEN 
            WRITE(*,'(1x,a9,i3,a19)') "volume = ",iter/(iters/10)*10,
     $           " % LAR computations"
            outsurf = .TRUE.
         ELSEIF(istep-maskpsi>=0 .AND. (istep+maskpsi)<=mstep)THEN
            !Print details at specified surfaces
            DO i = 1,10
               dpsi = ABS(psifac(istep)-outsurfs(i))
               IF(dpsi<=ABS(psifac(istep-maskpsi)-outsurfs(i)).AND.
     $               dpsi<=ABS(psifac(istep+maskpsi)-outsurfs(i))) 
     $              outsurf = .TRUE.
            ENDDO            
         ENDIF
         lar%xs(iter-1) = psifac(istep)
         larp%xs(iter-1)= psifac(istep)
         rot%xs(iter-1) = psifac(istep)
         rotp%xs(iter-1)= psifac(istep)
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL spline_eval(kin,psifac(istep),1)
         !Rotation frequencies: electric and density and temp gradient
         welec = wefac*kin%f(5)
         wdian =-twopi*kin%f(s+2)*kin%f1(s)/(chrg*chi1*kin%f(s))
         wdiat =-twopi*kin%f1(s+2)/(chrg*chi1)
         ! Bounce and precession frequencies
         q     = sq%f(4)
         epsr  = epsa*SQRT(psifac(istep)) 
         wt    = SQRT(2*kin%f(s+2)/mass)/(q*ro)
         wg    = chrg*bo/mass
         wbhat = (pi/4)*SQRT(epsr/2)*wt
         wdhat = q**3*wt**2/(4*epsr*wg)*wdfac

         rot%fs(iter-1,1) = welec+wdian+wdiat
         rotp%fs(iter-1,1)= welec+wdian+wdiat
         b_out(1,1,iter,1,:) = (/psifac(istep),q,epsr,wt,wg,wbhat,wdhat,
     $        welec,wdian,wdiat/)

         DO sigma = 0,pass
            DO ll = -nl,nl
               ls    = ll-sigma*nn*sq%f(4)
c-----------------------------------------------------------------------
c     Energy integration
c-----------------------------------------------------------------------
               xspl  = xintegral(wdian,wdiat,welec,wdhat,wbhat,ls,
     $              collision,s,ximag,xmax,.FALSE.)
               IF(offset_flag) oxspl  = xintegral(wdian,wdiat,welec,
     $              wdhat,wbhat,ls,collision,s,ximag,xmax,.TRUE.)
c-----------------------------------------------------------------------
c     Kappa integration and sum over perturbations
c-----------------------------------------------------------------------
               CALL spline_alloc(kspl,nk,1)
               kspl%fs(:,1) = 0
               DO k = 0,nk
                  kappa = (k+1)/(nk+2.0)
                  kspl%xs(k) = kappa
                  CALL spline_eval(ellipk,kappa,0)
                  DO ipert = 1,mpert
                     ! Don't let things go beyond precalculated spline
                     IF(ABS(mfac(ipert)-nn*q-ll) > fmnl%ys(fmnl%my))THEN 
                        PRINT *,mfac(ipert),nn*q,ll
                        CALL pentio_stop("m-nq-l out of fkmnql range")
                     ENDIF
                     ! Calculate for forward and backward banana
                     CALL bicube_eval(fmnl,kappa,mfac(ipert)-nn*q-ll,0)
                     fm(ipert,1) = fmnl%f(1)
                     CALL bicube_eval(fmnl,kappa,mfac(ipert)-nn*q+ll,0)
                     fm(ipert,2) = fmnl%f(1)
                  ENDDO
                  DO ipert = 1,mpert
                     DO jpert = 1,mpert
                        kspl%fs(k,1) = kspl%fs(k,1)+
     $                       2*kappa*0.25*(fm(ipert,1)*fm(jpert,1)+
     $                       fm(ipert,2)*fm(jpert,1)+
     $                       fm(ipert,1)*fm(jpert,2)+
     $                       fm(ipert,2)*fm(jpert,2))*
     $                       REAL(lagbpar(istep,ipert)*CONJG(
     $                       lagbpar(istep,jpert)))/(4*ellipk%f(1))
                     ENDDO
                  ENDDO
               ENDDO
               CALL spline_fit(kspl,'extrap')
               CALL spline_int(kspl)
c-----------------------------------------------------------------------
c     Record a sampling of intermediate steps for debuging
c-----------------------------------------------------------------------
               IF(outsurf)THEN
                  iout = MAX(MIN(iter/(iters/10),10),1)
                  DO k = 0,nk
                     p_out(1,k+1,iout,ll+nl+1,:)=(/psifac(istep)
     $                    ,kspl%xs(k),kspl%fs(k,1),kspl%fsi(k,1)/)
                  ENDDO
                  DO i = 0,nx,nx/99
                     ix = MIN(i/(nx/99)+1,100) 
                     x_out(1,ix,iout,ll+nl+1,:)=(/psifac(istep),
     $                    xspl%xs(i),REAL(xspl%fs(i,1)),
     $                    AIMAG(xspl%fs(i,1)),REAL(xspl%fsi(i,1)),
     $                    AIMAG(xspl%fsi(i,1)) /)
                  ENDDO
               ENDIF
               IF(sigma==0) lar%fs(iter-1,ll+nl+1) = sq%f(3)*
     $              kspl%fsi(kspl%mx,1)*0.5*xspl%fsi(xspl%mx,1)*
     $              SQRT(epsr/(2*pi**3))*nn*nn*kin%f(s)*kin%f(s+2)
               IF(sigma==1) larp%fs(iter-1,ll+nl+1) = sq%f(3)*
     $              kspl%fsi(kspl%mx,1)*0.5*xspl%fsi(xspl%mx,1)*
     $              SQRT(epsr/(2*pi**3))*nn*nn*kin%f(s)*kin%f(s+2)
               IF(offset_flag .AND. sigma==0) 
     $              rot%fs(iter-1,ll+nl+1+1)=rot%fs(iter-1,1)
     $              *oxspl%fsi(xspl%mx,1)/xspl%fsi(xspl%mx,1)
               IF(offset_flag .AND. sigma==1) 
     $              rotp%fs(iter-1,ll+nl+1+1)=rotp%fs(iter-1,1)
     $              *oxspl%fsi(xspl%mx,1)/xspl%fsi(xspl%mx,1)
               CALL spline_dealloc(kspl)
               CALL cspline_dealloc(xspl)
            ENDDO               ! ll
         ENDDO                  ! trapped/passing
      ENDDO                     ! psi

c-----------------------------------------------------------------------
c     Final integration over normalized flux
c-----------------------------------------------------------------------
      CALL cspline_fit(lar,"extrap")
      CALL cspline_int(lar)
      totals = SUM(lar%fsi(lar%mx,:))
      PRINT *, "Total LAR torque = ",REAL(totals)
      PRINT *, "alpha/s  = ",REAL(totals)/(-1*AIMAG(totals))

c-----------------------------------------------------------------------
c     Write output files
c-----------------------------------------------------------------------
      outfile ="pent_lar_n"//TRIM(sn)
      outtext ="Large aspect ratio, trapped particle "//
     $     "torque (Nm) and energy (J)"
      IF(offset_flag)THEN
         CALL cspline_fit(rot,"extrap")
         CALL pentio_write_profile(outfile,outtext,lar,rot)
      ELSE
         CALL pentio_write_profile(outfile,outtext,lar)
      ENDIF
      CALL cspline_dealloc(rot)

      IF(passing_flag)THEN
         outfile ="pent_lar_passing_n"//TRIM(sn)
         outtext ="Large aspect ratio, passing particle "//
     $        "torque (Nm) and energy (J)"
         IF(offset_flag)THEN
            CALL cspline_fit(rotp,"extrap")
            CALL pentio_write_profile(outfile,outtext,larp,rotp)
         ELSE
            CALL pentio_write_profile(outfile,outtext,larp)
         ENDIF
      ENDIF
      CALL cspline_dealloc(rotp)
      
      IF(energy_flag) THEN
         outfile = "pent_lar_energy_n"//TRIM(sn)
         outtext = "Sampling of large aspect ratio energy integrands"
     $        //newl//"Normalizations are as follows:"
     $        //newl//"x: E/T"//newl//
     $        "T_phi: ??T_phi/(-2n^2TN*chi'/sqrt(pi))"//newl//
     $        "2ndeltaW: ??2*n*deltaW /(-2n^2TN*chi'/sqrt(pi))"
         WRITE(outlbls,'(6(1x,a16))') "psi","x","T_phi","2ndeltaW",
     $        "int(T_phi)","int(2ndeltaW)"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(x_out),x_out,6)
      ENDIF

      IF(pitch_flag) THEN
         outfile ="pent_lar_pitch_n"//TRIM(sn)
         outtext ="Sampling of large aspect ratio pitch angle functions"
         WRITE(outlbls,'(4(1x,a16))') "psi","kappa","deltawl",
     $        "int(deltawl)"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(p_out),p_out) 
      ENDIF

      outfile ="pent_lar_freq_n"//TRIM(sn)
      outtext ="Radial profile of various frequencies"//newl//
     $     "Frequencies are normalized as follows:"//newl//
     $     "omega_b = omega/sqrt(x)"//newl//"omega_D = omega/x"
      WRITE(outlbls,'(10(1x,a16))') "psi","q","eps_r","omega_t",
     $     "omega_g","omega_b","omega_D","omega_E","omega_*N","omega_*T"
      CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $     SHAPE(b_out),b_out) 

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE profile_lar



      END MODULE profile_mod

