c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     calculate ntv torque and kinetic perturbed energy
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. profile_mod
c     1. profile_gen
c     2. profile_rlar
c-----------------------------------------------------------------------
c     subprogram 0. profile_mod.
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
c     subprogram 1. profile_gen.
c     Calculate NTV on a magnetic surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE profile_gen(lagbpar,divxprp,nl,maskpsi,ptfac,ximag,
     $     xmax,wefac,wdfac,collision,outsurfs,clar)

      LOGICAL, INTENT(IN) :: clar
      INTEGER, INTENT(IN) :: nl,maskpsi
      REAL(r8), INTENT(IN) :: ptfac,ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10), INTENT(IN) :: outsurfs
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf), INTENT(IN) :: 
     $     lagbpar,divxprp
      CHARACTER(32), INTENT(IN) :: collision

      LOGICAL :: outsurf
      INTEGER :: istep,ilast,inext,itheta,i,ilmda,ll,sigma,ibmax,
     $     iters,iter,indx,s,iout,ix,pass,xsigma
      REAL(r8):: chrg,mass,epsa,epsr,lmda,welec,wdian,wdiat,bhat,dhat,
     $     wbbar,wdbar,bo,ls,dpsi,bmax,bmin,dt,kk,k,lmdamin,lmdamax,dq
      REAL(r8), DIMENSION(:), ALLOCATABLE :: bpts
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: htheta,tdt,ldl
      COMPLEX(r8) :: totals
      COMPLEX(r8), DIMENSION(-nl:nl) :: lmdaint,lmdaoff
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: pl,act0

      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
      CHARACTER(32) :: outfile,lbl
      CHARACTER(512) :: outtext,outlbls
      
      TYPE(spline_type) :: omega,tspl,vpar,omegabar,djdjbar
      TYPE(cspline_type) o1spl,action,xspl,ntvdl,offdl,rot,rotp

c-----------------------------------------------------------------------
c     Set species.
c-----------------------------------------------------------------------
      chrg = e*icharge
      IF(chrg >= 0)THEN
         s = 1
         mass = mp*imass
      ELSE IF(icharge == -1)THEN
         s = 2
         mass = me
      ELSE
         CALL pentio_stop('Negative charge must be electrons (-1).')
      ENDIF
      pass = 0
      xsigma=0
      IF(passing_flag) pass=1
      IF(divxprp_flag .AND. .NOT. clar) xsigma=1 !LAR->no fieldline bend
c-----------------------------------------------------------------------
c     Calculations.
c----------------------------------------------------------------------- 
      IF(clar)THEN
         PRINT *, "PENT - circular large-aspect-ratio calculation v1.0"
         lbl = 'clar_'
      ELSE
         PRINT *, "PENT - general calculation v1.0"
         lbl = ''
      ENDIF

      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)

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
      CALL cspline_alloc(action,nt,1)
      ALLOCATE(htheta(0:nt,2),pl(0:nt),act0(0:nt),tdt(0:1,0:nt))
      omega%xs(:) = tnorm(:)
      action%xs(:) = tnorm(:)

      CALL spline_alloc(omegabar,nlmda,2)
      CALL spline_alloc(djdjbar,nlmda,2*nl+1)
      ALLOCATE(ldl(0:1,0:nlmda))

      DO istep=1,mstep,maskpsi
         iter = (istep-1)/maskpsi+1
         ! print details of specified surfaces
         IF (MOD(iter,iters/10)==0 .OR. iter==1)THEN 
            !WRITE(*,'(1x,a9,i3,a19)') "volume = ",iter/(iters/10)*10,
     $      !     " % NTV computations"
            CALL pentio_progress(iter/(iters/10))
         ENDIF
         outsurf = .FALSE.
         ilast = MAX(1,istep-maskpsi)
         inext = MIN(mstep,istep+maskpsi)
         DO i = 1,10
            dpsi = ABS(psifac(istep)-outsurfs(i))
            IF(dpsi<=ABS(psifac(ilast)-outsurfs(i)) .AND.
     $           dpsi<=ABS(psifac(inext)-outsurfs(i)))THEN 
               outsurf = .TRUE.
               iout = i
            ENDIF
         ENDDO

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
         ! clebsch conversion now in djdt o1*EXP(-twopi*ifac*nn*q*theta)
         CALL spline_fit(tspl,"periodic")
         CALL cspline_fit(o1spl,"periodic")
         bmax = MAXVAL(tspl%fs(:,1),DIM=1)
         bmin = MINVAL(tspl%fs(:,1),DIM=1)
         ibmax= 0-1+MAXLOC(tspl%fs(:,1),DIM=1)
         IF(bmax/=tspl%fs(ibmax,1)) CALL pentio_stop("ibmax off")
         DO i=2,8 !8th smallest so spline has more than 1 pt
            bmin = MINVAL(tspl%fs(:,1),MASK=tspl%fs(:,1)>bmin,DIM=1)
         ENDDO         
         DO sigma = 0,pass
c-----------------------------------------------------------------------
c     Calculate pitch dependent variables
c      Concentrate near trapped/passing boundary and deeply trapped lim
c-----------------------------------------------------------------------
            lmdamin = sigma*1e-6/bmax+
     $           (1-sigma)*(1+(1/bmin-1/bmax)*ptfac)/bmax
            lmdamax = sigma*(1-ptfac)/bmax+(1-sigma)/bmin
            IF(clar)THEN
               IF(sigma==1) EXIT ! No passing particle contrib. in clar
               lmdamin = MAX(lmdamin,
     $              1.0/(bo*(1+2*epsr*(1.0-1e-6)-epsr)) ) !kk=1
               lmdamax = MIN(lmdamax,1./(bo*(1+2*epsr*0-epsr))) !kk=0
            ENDIF
            ldl = powspace(lmdamin,lmdamax,1,nlmda+1,"both")
            DO ilmda=0,nlmda
               lmda = ldl(0,ilmda)
c-----------------------------------------------------------------------
c     Determine bounce points
c-----------------------------------------------------------------------
               vpar%fs(:,1) = 1-lmda*tspl%fs(:,1)
               CALL spline_fit(vpar,"extrap")
               IF(sigma==0)THEN
                  CALL spl_roots(bpts,vpar,1)
                  CALL spline_eval(vpar,SUM(bpts(1:2))/2,0)
                  ! the usual case with bpts centered around 0 requires shift
                  IF(vpar%f(1)<=0.0 .AND. sigma==0) 
     $                 bpts(1:2)=(/bpts(2)-1.0,bpts(1)/)
                  ! Debugger
                  IF(bpts(1)>=bpts(2) .AND. sigma==0)THEN
                     PRINT *," Error at (psi,lmda) = ",psifac(istep),
     $                    lmda/ldl(0,nlmda)
                     PRINT *,"Bounce pts = ",bpts
                     CALL pentio_stop("Bounce point crossing."
     $                    //" Check if repeated zeros in spl_roots.")
                  ENDIF
                  ! was bpts(1-sigma,2-sigma) before IF
                  CALL powergrid(tnorm,tdt,bpts(1:2),4,"both")
                  DEALLOCATE(bpts)
               ELSE
                  ALLOCATE(bpts(2))
                  bpts(1) = tspl%xs(ibmax)
                  bpts(2) = tspl%xs(ibmax)+1
                  CALL powergrid(tnorm,tdt,bpts,2,"both")
                  DEALLOCATE(bpts)
               ENDIF
c-----------------------------------------------------------------------
c     Bounce averages. Bounce & precession freqs, variation in action
c-----------------------------------------------------------------------
               omega%fs(0,:) = 0
               omega%fs(omega%mx,:) = 0.0
               action%fs(0,:) = 0
               action%fs(action%mx,:) = 0.0
               act0(:) = 0.0
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
                  act0(itheta) = dt*tspl%f(4)*tspl%f(1)*(
     $                 xsigma*o1spl%f(2)*SQRT(vpar%f(1))+
     $                 o1spl%f(1)*(1-1.5*lmda*tspl%f(1))/
     $                 SQRT(vpar%f(1)))
     $                 *EXP(-twopi*ifac*nn*sq%f(4)*(tdt(0,itheta)
     $                 -tdt(0,1))) ! doesn't matter since |dJ|^2
               ENDDO
               call spline_fit(omega,"extrap")
               CALL spline_int(omega)
               htheta(:,1)= omega%fsi(:,1)
               htheta(:,2)= omega%fsi(:,2) ! missing analytic piece?
               htheta(:,1)=htheta(:,1)/((2-sigma)*omega%fsi(omega%mx,1))
               htheta(:,2)=htheta(:,2)/((2-sigma)*omega%fsi(omega%mx,2))
               ! normalized bounce averaged frequencies
               IF(clar)THEN
                  !kk = MIN((1/bo - lmda +epsr*lmda)/(2*epsr*lmda),
     $            !     1.0-1e-6)   
                  kk = (1/bo - lmda +epsr*lmda)/(2*epsr*lmda)
                  k = SQRT(kk)
                  CALL spline_eval(ellipk,kk,0)
                  CALL spline_eval(ellipe,kk,0)
                  wbbar = pi*SQRT(2*epsr*lmda*bo)
     $                 /(4*sq%f(4)*ro*ellipk%f(1))
                  wdbar = (2*sq%f(4)*lmda*(ellipe%f(1)/ellipk%f(1)-0.5)
     $                 /(ro**2*epsr))*wdfac
                  htheta(:,1) = (tdt(0,:)-tdt(0,0))/2.0
                  ! Null torque outside clar trapped phase space
                  !IF(kk==1.0-1e-6 .OR. sigma==1) wbbar=0
               ELSE
                  wbbar = twopi/((2-sigma)*omega%fsi(omega%mx,1))
                  wdbar = wdfac*wbbar*2*(2-sigma)*omega%fsi(omega%mx,2)
               ENDIF
               omegabar%xs(ilmda)   = lmda
               omegabar%fs(ilmda,1) = wbbar
               omegabar%fs(ilmda,2) = wdbar
               djdjbar%xs(ilmda)   = lmda
               DO ll=-nl,nl
                  ls = REAL(ll+sigma*nn*sq%f(4))
                  IF(welec<SQRT(2*kin%f(s+2)/mass)*wdbar*3)THEN 
                     ! mag precession dominates drift
                     pl(:) = EXP(-twopi*ifac*ls*htheta(:,2))
                  ELSE 
                     ! electric precession dominates drift (normal)
                     pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
                  ENDIF
                  pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
                  action%fs(:,1) =CONJG(act0(:))*(pl(:)+(1-sigma)/pl(:))
                  CALL cspline_fit(action,"extrap")
                  CALL cspline_int(action)
                  djdjbar%xs(ilmda) = lmda
                  djdjbar%fs(ilmda,ll+nl+1) = 
     $                 REAL(action%fsi(action%mx,1))**2+
     $                 AIMAG(action%fsi(action%mx,1))**2
c-----------------------------------------------------------------------
c     Optionally write bounce integration to file
c-----------------------------------------------------------------------
                  IF(outsurf .AND. bounce_flag .AND. 
     $                 MOD(ilmda,nlmda/4)==0)THEN
                     WRITE(UNIT=sl,FMT='(I2)') ll
                     sl = ADJUSTL(sl)
                     OPEN(UNIT=out_unit,FILE="pent_bounce_n"//TRIM(sn)
     $                    //"_l"//TRIM(sl)//".out",STATUS="UNKNOWN",
     $                    POSITION="APPEND")
                     DO i=0,nt
                        WRITE(out_unit,'(1x,15(1x,es16.8E3))') kin%x,
     $                       lmda,omega%xs(i),tdt(:,i),
     $                       REAL(action%fs(i,1)),AIMAG(action%fs(i,1)),
     $                       REAL(action%fsi(i,1)),
     $                       AIMAG(action%fsi(i,1)),
     $                       omega%fs(i,:),omega%fsi(i,:),htheta(i,:)
                     ENDDO
                     CALL ascii_close(out_unit)
                  ENDIF
               ENDDO ! l loop
            ENDDO ! lmda loop
            CALL spline_fit(omegabar,'extrap')
            CALL spline_fit(djdjbar,'extrap')
            
c-----------------------------------------------------------------------
c     Lambda and E integrations
c-----------------------------------------------------------------------
            ! normalisation factors without energy dependence
            bhat = SQRT(2*kin%f(s+2)/mass)
            dhat = (kin%f(s+2)/chrg)
            DO ll=-nl,nl
               ! hack for mdc2 benchmark for mars-k analytic x integrals
               IF(mdc2_flag .AND. ll/=0)THEN
                  dhat = 0.0
               ELSE
                  dhat = (kin%f(s+2)/chrg)
               ENDIF
               IF(lmdalsode_flag)THEN
                  lmdaint(ll) =lsode_lambda(wdian,wdiat,welec,bhat,dhat,
     $                 omegabar,djdjbar,ll,ximag,xmax,collision,s,sigma,
     $                 .FALSE.,outsurf,lbl)
                  IF(offset_flag) lmdaoff(ll) =lsode_lambda(wdian,
     $                 wdiat,welec,bhat,dhat,omegabar,djdjbar,ll,
     $                 ximag,xmax,collision,s,sigma,.TRUE.,outsurf,lbl)
               ELSE
                  lmdaint(ll)=intspl_lambda(wdian,wdiat,welec,bhat,dhat,
     $                 omegabar,djdjbar,ll,ximag,xmax,collision,s,sigma,
     $                 .FALSE.,outsurf,lbl)
                  IF(offset_flag) lmdaoff(ll) =lsode_lambda(wdian,
     $                 wdiat,welec,bhat,dhat,omegabar,djdjbar,ll,
     $                 ximag,xmax,collision,s,sigma,.TRUE.,outsurf,lbl)
               ENDIF
            ENDDO               
c-----------------------------------------------------------------------
c     Final profile in normalized flux (Nm)
c-----------------------------------------------------------------------
            IF(sigma==0) ntv%fs(iter-1,:)=-lmdaint(:)*
     $           (chi1/twopi)*2*nn**2*kin%f(s)*kin%f(s+2)/SQRT(pi)
            IF(sigma==1) ntvp%fs(iter-1,:)=-lmdaint(:)*
     $           (chi1/twopi)*2*nn**2*kin%f(s)*kin%f(s+2)/SQRT(pi)
            IF(offset_flag)THEN
               IF(sigma==0) rot%fs(iter-1,2:) = rot%fs(iter-1,1)*
     $                 lmdaoff/(lmdaint-lmdaoff)
               IF(sigma==1) rotp%fs(iter-1,2:) = rotp%fs(iter-1,1)*
     $                 lmdaoff/(lmdaint-lmdaoff)
            ENDIF

            ! Remove singularities near resonant surfaces
            IF(treg_flag)THEN
               dq = sq%f(4) - NINT(sq%f(4))
               if(sigma==0) ntv%fs(iter-1,:)=ntv%fs(iter-1,:)
     $              *dq**4/(dq**4+0.01**4)
               if(sigma==1) ntvp%fs(iter-1,:)=ntvp%fs(iter-1,:)
     $              *dq**4/(dq**4+0.01**4)
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
      PRINT *, "-----------------------------------"
      PRINT "(a24,es11.3E3)", "Total torque = ",
     $     REAL(totals)
      PRINT "(a24,es11.3E3)", "Total Kinetic Energy = ",
     $     AIMAG(totals)/(2*nn)
      PRINT "(a24,es11.3E3)", "alpha/s  = ",
     $     REAL(totals)/(-1*AIMAG(totals))
      PRINT *, "-----------------------------------"

      CALL spline_dealloc(tspl)
      CALL cspline_dealloc(o1spl)
      CALL spline_dealloc(vpar)

      CALL spline_dealloc(omega)
      CALL cspline_dealloc(action)

c-----------------------------------------------------------------------
c     Write output files
c-----------------------------------------------------------------------
      outfile ="pent_"//TRIM(lbl)//"n"//TRIM(sn)
      outtext ="Trapped particle torque (Nm) and energy (J)"
      IF(offset_flag)THEN
         CALL cspline_fit(rot,"extrap")
         CALL pentio_write_profile(outfile,outtext,ntv,rot)
      ELSE
         CALL pentio_write_profile(outfile,outtext,ntv)
      ENDIF
      CALL cspline_dealloc(rot)

      IF(passing_flag)THEN
         outfile ="pent_"//TRIM(lbl)//"passing_n"//TRIM(sn)
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

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE profile_gen


c-----------------------------------------------------------------------
c     subprogram 5. profile_rlar.
c     Calculate NTV using reduced large aspect ratio approximations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE profile_rlar(lagbpar,divxprp,nl,maskpsi,ximag,xmax,
     $     wefac,wdfac,collision,outsurfs)
      
      INTEGER, INTENT(IN) :: nl,maskpsi
      REAL(r8), INTENT(IN) :: ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10), INTENT(IN) :: outsurfs
      COMPLEX(r8), DIMENSION(mstep,mpert), INTENT(IN) :: lagbpar,divxprp
      CHARACTER(32), INTENT(IN) :: collision

      LOGICAL :: outsurf
      INTEGER, PARAMETER :: nk = 50
      INTEGER :: istep,k,i,ipert,jpert,ll,sigma,
     $     pass,iters,iter,indx,s,nsplt,iout,ix
      REAL(r8) :: chrg,mass,kappa,wbhat,wdhat,welec,wdian,wdiat,
     $     epsa,epsr,wt,wg,bo,ls,dpsi
      REAL(r8), DIMENSION(mpert,2):: fm
      COMPLEX(r8) :: totals,xint,oint

      REAL(r8), DIMENSION(1,100,10,1:2*nl+1,6)  :: x_out
      REAL(r8), DIMENSION(1,nk+1,10,1:2*nl+1,4) :: p_out
      REAL(r8), DIMENSION(1,1,mstep,1,10) :: b_out
      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
      CHARACTER(32) :: outfile,lbl = 'rlar_'
      CHARACTER(512) :: outtext,outlbls
      
      TYPE(spline_type) :: kspl
      TYPE(cspline_type):: oxspl,xspl,rot,rotp

c-----------------------------------------------------------------------
c     Set species.
c-----------------------------------------------------------------------
      PRINT *, "PENT - reduced large-aspect-ratio calculation v1.0"

      chrg = e*icharge
      IF(chrg >= 0)THEN
         s = 1
         mass = mp*imass
      ELSE IF(icharge == -1)THEN
         s = 2
         mass = me
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
      CALL cspline_alloc(ntv,iters-1,2*nl+1)
      CALL cspline_alloc(ntvp,iters-1,2*nl+1)
      CALL cspline_alloc(rot,iters-1,2*nl+1+1)
      CALL cspline_alloc(rotp,iters-1,2*nl+1+1)

      DO istep=1,mstep,maskpsi
         iter = (istep-1)/maskpsi+1
         outsurf = .FALSE.
         IF (MOD(iter,iters/10)==0 .OR. iter==1) THEN 
            !WRITE(*,'(1x,a9,i3,a19)') "volume = ",iter/(iters/10)*10,
     $      !     " % NTV computations"
            CALL pentio_progress(iter/(iters/10))
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
         ntv%xs(iter-1) = psifac(istep)
         ntvp%xs(iter-1)= psifac(istep)
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
            IF(sigma==1) EXIT ! wb,d not valid for passing
            DO ll = -nl,nl
               ls    = ll-sigma*nn*sq%f(4)
               WRITE(UNIT=sl,FMT='(I2)') ll
               sl = ADJUSTL(sl)
c-----------------------------------------------------------------------
c     Energy integration
c-----------------------------------------------------------------------
               IF(xlsode_flag)THEN
                  xint = lsode_x(wdian,wdiat,welec,wdhat,wbhat,ls,
     $                 ximag,xmax,collision,s,sigma,.FALSE.,
     $                 (outsurf .AND. energy_flag),
     $                 lbl)
                  IF(offset_flag) oint = lsode_x(wdian,wdiat,welec,
     $                 wdhat,wbhat,ls,ximag,xmax,collision,s,sigma,
     $                 .TRUE.,.FALSE.,'')
               ELSE
                  xint = intspl_x(wdian,wdiat,welec,wdhat,wbhat,ls,
     $                 ximag,xmax,collision,s,sigma,.FALSE.,outsurf,
     $                 lbl)
                  IF(offset_flag) oint = intspl_x(wdian,wdiat,welec,
     $                 wdhat,wbhat,ls,ximag,xmax,collision,s,sigma,
     $                 .TRUE.,.FALSE.,'')
               ENDIF   
c               IF(offset_flag) oxspl  = xintegral(wdian,wdiat,welec,
c     $              wdhat,wbhat,ls,collision,s,ximag,xmax,.TRUE.)
c-----------------------------------------------------------------------
c     Kappa integration and sum over perturbations
c-----------------------------------------------------------------------
               CALL spline_alloc(kspl,nk,1)
               kspl%fs(:,1) = 0
               DO k = 0,nk
                  kappa = (k+1)/(nk+2.0) !0<=kappa<1 (spline fails at 1)
                  kspl%xs(k) = kappa
                  ! use K(k) not K(k^2) since integrating over kappa
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
               ENDIF
c-----------------------------------------------------------------------
c     Complete torque expression
c-----------------------------------------------------------------------
               IF(sigma==0) ntv%fs(iter-1,ll+nl+1) = sq%f(3)*
     $              kspl%fsi(kspl%mx,1)*0.5*(-xint)*
     $              SQRT(epsr/(2*pi**3))*nn*nn*kin%f(s)*kin%f(s+2)
               IF(sigma==1) ntvp%fs(iter-1,ll+nl+1) = sq%f(3)*
     $              kspl%fsi(kspl%mx,1)*0.5*(-xint)*
     $              SQRT(epsr/(2*pi**3))*nn*nn*kin%f(s)*kin%f(s+2)
               IF(offset_flag .AND. sigma==0) 
     $              rot%fs(iter-1,ll+nl+1+1)=rot%fs(iter-1,1)
     $              *oint/xint
               IF(offset_flag .AND. sigma==1) 
     $              rotp%fs(iter-1,ll+nl+1+1)=rotp%fs(iter-1,1)
     $              *oint/xint
               CALL spline_dealloc(kspl)
c               CALL cspline_dealloc(xspl)
            ENDDO               ! ll
         ENDDO                  ! trapped/passing
      ENDDO                     ! psi

c-----------------------------------------------------------------------
c     Final integration over normalized flux
c-----------------------------------------------------------------------
      CALL cspline_fit(ntv,"extrap")
      CALL cspline_int(ntv)
      totals = SUM(ntv%fsi(ntv%mx,:))
      PRINT *, "-----------------------------------"
      PRINT "(a24,es11.3E3)", "Total torque = ",
     $     REAL(totals)
      PRINT "(a24,es11.3E3)", "Total Kinetic Energy = ",
     $     AIMAG(totals)/(2*nn)
      PRINT "(a24,es11.3E3)", "alpha/s  = ",
     $     REAL(totals)/(-1*AIMAG(totals))
      PRINT *, "-----------------------------------"

c-----------------------------------------------------------------------
c     Write output files
c-----------------------------------------------------------------------
      outfile ="pent_rlar_n"//TRIM(sn)
      outtext ="Reduced large-aspect-ratio, trapped particle "//
     $     "torque (Nm) and energy (J)"
      IF(offset_flag)THEN
         CALL cspline_fit(rot,"extrap")
         CALL pentio_write_profile(outfile,outtext,ntv,rot)
      ELSE
         CALL pentio_write_profile(outfile,outtext,ntv)
      ENDIF
      CALL cspline_dealloc(rot)

      IF(passing_flag)THEN
         CALL cspline_fit(ntvp,"extrap")
         CALL cspline_int(ntvp)
         outfile ="pent_rlar_passing_n"//TRIM(sn)
         outtext ="Reduced large-aspect-ratio, passing particle "//
     $        "torque (Nm) and energy (J)"
         IF(offset_flag)THEN
            CALL cspline_fit(rotp,"extrap")
            CALL pentio_write_profile(outfile,outtext,ntvp,rotp)
         ELSE
            CALL pentio_write_profile(outfile,outtext,ntvp)
         ENDIF
      ENDIF
      CALL cspline_dealloc(rotp)

      IF(pitch_flag) THEN
         outfile ="pent_rlar_pitch_n"//TRIM(sn)
         outtext ="Reduced large-aspect-ratio pitch angle functions"
         WRITE(outlbls,'(4(1x,a16))') "psi","kappa","deltawl",
     $        "int(deltawl)"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(p_out),p_out) 
      ENDIF

      IF(kinetic_flag)THEN
         outfile ="pent_rlar_freq_n"//TRIM(sn)
         outtext ="Radial profile of various frequencies"//newl//
     $        "Frequencies are normalized as follows:"//newl//
     $        "omega_b = omega/sqrt(x)"//newl//"omega_D = omega/x"
         WRITE(outlbls,'(10(1x,a16))') "psi","q","eps_r","omega_t",
     $        "omega_g","omega_b","omega_D","omega_E","omega_*N",
     $        "omega_*T"
         CALL pentio_write_bxpflag(outfile,outtext,outlbls,
     $        SHAPE(b_out),b_out) 
      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE profile_rlar



      END MODULE profile_mod

