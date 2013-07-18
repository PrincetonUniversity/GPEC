c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     various diagnostics of individual components
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.... is a mess
c-----------------------------------------------------------------------
c     0. diagnostic_mod
c     ?
c     ?
c-----------------------------------------------------------------------
c     subprogram 0. diag_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE diagnostic_mod
      USE profile_mod
      USE energy_mod

      IMPLICIT NONE

      CONTAINS


      SUBROUTINE diag_unpack()
      REAL(r8), DIMENSION(3) :: a
      REAL(r8) :: a1,a2,a3
      a = (/ 1,2,3 /)
      !a1,a2,a3 = a
      a1 = a(1)
      a2 = a(2)
      a3 = a(3)
      PRINT *, a
      PRINT *, a1,a2,a3
      RETURN
      END SUBROUTINE diag_unpack

c-----------------------------------------------------------------------
c     function a. diag_inner.
c     Test function namespaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION diag_inner(i)
      !INTEGER, INTENT(INOUT) :: i
      IMPLICIT NONE
      INTEGER diag_inner
      INTEGER i,j
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      j = 1
      i = i+1
      PRINT *,"inner ",j,i
      diag_inner = i
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION diag_inner

c-----------------------------------------------------------------------
c     function a. diag_inner.
c     Test function namespaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION diag_outer(i)
      IMPLICIT NONE
      INTEGER diag_outer
      INTEGER i,j
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      diag_outer = 0
      PRINT *,"outer start ",i,diag_outer
      DO j = 0,3
         diag_outer = diag_inner(i+j)
         PRINT *,'outer ',j,i," -->",diag_outer
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION diag_outer




c-----------------------------------------------------------------------
c     subprogram 1. diag_grid.
c     Test grid making module.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_grid()
      REAL(r8), DIMENSION(0:1,0:nx) :: ydylow1,ydylow2,ydyup,ydyboth
      REAL(r8) :: pow = 3.0
      REAL(r8), DIMENSION(0:1) :: set
      REAL(r8), DIMENSION(4) :: set4
      INTEGER :: i
      LOGICAL :: test
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      
c$$$      set = (/0.0,0.0,10.0/)
c$$$      CALL  grid_power(xnorm,ydylow1,set,pow,"lower")
c$$$      set4 = (/0.0,3.0,6.0,10.0/)
c$$$      CALL  grid_power(xnorm,ydylow2,set4,pow,"lower")
c$$$      CALL  grid_power(xnorm,ydyup,set4,pow,"upper")
c$$$      set4 = (/0.0,5.0,7.0,10.0/)
c$$$      CALL  grid_power(xnorm,ydyboth,set4,pow,"both")

      set = (/0.0,10.0/)
      CALL  powergrid(xnorm,ydylow1,set,3,"lower")
      set4 = (/0.0,3.0,6.0,10.0/)
      CALL  powergrid(xnorm,ydylow2,set4,3,"lower")
      CALL  powergrid(xnorm,ydyup,set4,3,"upper")
      set4 = (/0.0,5.0,7.0,10.0/)
      CALL  powergrid(xnorm,ydyboth,set4,3,"none")
      
      print *,"*************** - grid diagnostic - ***************"
      test = .TRUE.
      DO i=1,nx
         IF(ydylow2(0,i-1)>=ydylow2(0,i)) test = .FALSE.
      ENDDO
      print *,"y always increasing is ",test
      print *,"---------------------------------------------------"

      CALL ascii_open(out_unit,"diag_grid.out","UNKNOWN")
      WRITE(out_unit,*)"NON-AMBIPOLAR TRANSPORT and TOROIDAL TORQUE: "//
     $     "Hardcoded grid diagnostic for developer use only"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,9(1x,a16))')"x","y1","y1prime","y2","y2prime",
     $     "y3","y3prime","y4","y4prime"
      DO i = 0,nx
         WRITE(out_unit,'(1x,9(1x,es16.8E3))') xnorm(i),ydylow1(0,i),
     $        ydylow1(1,i),ydylow2(0,i),ydylow2(1,i),ydyup(0,i),
     $        ydyup(1,i),ydyboth(0,i),ydyboth(1,i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE diag_grid


c-----------------------------------------------------------------------
c     subprogram 2. diag_lsode_example.
c     Test lsode module.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_lsode_example
c      EXTERNAL  FEX, JEX
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
     $     MF, NEQ
      REAL*8  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
      NEQ = 3
      Y(1) = 1.
      Y(2) = 0.
      Y(3) = 0.
      T = 0.
      TOUT = .4
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-6
      ATOL(2) = 1.D-10
      ATOL(3) = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = 58
      LIW = 23
      MF = 21
      PRINT *,"*************** - lsode example - ***************"
      DO 40 IOUT = 1,12
        CALL LSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
        WRITE(6,20)  T, Y(1), Y(2), Y(3)
  20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
        IF (ISTATE .LT. 0)  GO TO 80
  40    TOUT = TOUT*10.
      WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
  60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
      STOP
  80  WRITE(6,90)  ISTATE
  90  FORMAT(///' Error halt.. ISTATE =',I3)
      STOP
      END SUBROUTINE diag_lsode_example

      SUBROUTINE  FEX (NEQ, T, Y, YDOT)
      INTEGER  NEQ
      REAL*8  T, Y(3), YDOT(3)
      YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
      YDOT(3) = 3.D7*Y(2)*Y(2)
      YDOT(2) = -YDOT(1) - YDOT(3)
      RETURN
      END SUBROUTINE FEX

      SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y(3), PD(NRPD,3)
      PD(1,1) = -.04
      PD(1,2) = 1.D4*Y(3)
      PD(1,3) = 1.D4*Y(2)
      PD(2,1) = .04
      PD(2,3) = -PD(1,3)
      PD(3,2) = 6.D7*Y(2)
      PD(2,2) = -PD(1,2) - PD(3,2)
      RETURN
      END SUBROUTINE JEX

c-----------------------------------------------------------------------
c     subprogram 2. diag_lsode_custom.
c     Test lsode module.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_lsode_custom
c      EXTERNAL  FEX, JEX
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     $     MF, NEQ
      REAL*8  ATOL(1), RTOL, RWORK(36), T, TOUT, Y(1)
      NEQ = 1
      Y(:) = (/0./)
      T = 0
      TOUT = 0.5
      ITOL = 2
      RTOL = 1.D-8
      ATOL(:) = (/1.D-9/)
      ITASK = 4 !don't overshoot tcrit defined in rwork(1)
      ISTATE = 1
      IOPT = 0
      LRW = 58
      LIW = 23
      MF = 10
      RWORK(1) = 1
      PRINT *,"*************** - lsode custom - ***************"
      DO 40 IOUT = 1,10
        IF(IOUT==10) TOUT=1 
        CALL LSODE (FCUSTOM, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JCUSTOM, MF)
        WRITE(6,20)  T, (/ Y(1), -2.*SQRT(1-T)+2 /)
  20    FORMAT(' At t =',D17.9,'   numeric,analytic =',2D14.6)
        IF (ISTATE .LT. 0)  GO TO 80
  40    TOUT = 1.-1./(10.**IOUT)
      WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
  60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
      STOP
  80  WRITE(6,90)  ISTATE
  90  FORMAT(///' Error halt.. ISTATE =',I3)
      STOP
      END SUBROUTINE diag_lsode_custom

      SUBROUTINE  FCUSTOM (NEQ, T, Y, YDOT)
      INTEGER  NEQ
      REAL*8  T, Y(1), YDOT(1)
      YDOT(1) = 1.0/SQRT(1-T)
      RETURN
      END SUBROUTINE FCUSTOM

      SUBROUTINE  JCUSTOM (NEQ, T, Y, ML, MU, PD, NRPD)
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y, PD(NRPD,1)
      PD(1,1) = 0
      RETURN
      END SUBROUTINE JCUSTOM





      FUNCTION xfun(wn,wt,we,wd,wb,l,ximag,xmax,nutype,s,offset,quiet)
c     Wrapper for lsode energy integration
      IMPLICIT NONE
      COMPLEX(r8) :: xfun
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax
      INTEGER, INTENT(IN) :: nutype,s,offset,quiet
c     grabbed from lsode example
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     $     MF, NEQ(4),i
      REAL*8  ATOL(2), RTOL, RWORK(52), x, xout, energy(9)
      CHARACTER(16) :: collision = "zero"
      Logical :: outsurf
      !REAL(r8), DIMENSION(:,:,:,:,:) :: x_out
      TYPE(cspline_type) xspl
      NEQ(1) = 2
      NEQ(2:4) = (/ nutype,s,offset /)
      energy(1:2) = (/ 0,0 /)
      energy(3:9) = (/ wn,wt,we,wd,wb,l,ximag /)
      x=1e-9
      xout = xmax
      ITOL = 2
      RTOL = 1.D-6
      ATOL(:) = (/1.D-3,1.D-3/)
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = 58
      LIW = 23
      MF = 10
      IF(quiet<=0) THEN
         ITASK = 2 !one step
         PRINT *,"Starting x integral lsode"
         CALL ascii_open(out_unit,"diag_lsodeIx.out","UNKNOWN")
         WRITE(out_unit,*)"NON-AMBIPOLAR TRANSPORT and TOROIDAL "//
     $        "TORQUE: LSODE x integral diagnostic for developer only"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,3(1x,a16))')"x","ReIx","ImIx"
         DO WHILE (x<xout)
            CALL LSODE(Fxintegral, NEQ, energy, x, xout, ITOL, RTOL,  
     $           ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, 
     $           JCUSTOM, MF)
            WRITE(out_unit,'(1x,3(1x,es16.8E3))') x,energy(1:2)
         ENDDO
         CALL ascii_close(out_unit)
         CALL ascii_open(out_unit,"diag_csplIx.out","UNKNOWN")
         xspl = xintegral(wn,wt,we,wd,wb,l,
     $                    collision,s,ximag,xmax,.FALSE.)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,5(1x,a16))')"x","ReIx","ImIx","Reix","Imix"
         DO i=0,xspl%mx
            WRITE(out_unit,'(1x,5(1x,es16.8E3))') xspl%xs(i),
     $           REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1)),
     $           REAL(xspl%fs(i,1)),AIMAG(xspl%fs(i,1))
         ENDDO
         CALL ascii_close(out_unit)
         PRINT *,' No. steps =',IWORK(11),',  No. f-s =',IWORK(12)
         PRINT *,'Ix = ',energy(1),' xspl = ',REAL(xspl%fsi(nx,1))
         !IF(REAL(xspl%fsi(nx,1)) .GE. 0.0) 
         IF(xspl%mx/=nx) CALL pentio_stop("check files")
      ELSE
         CALL LSODE(Fxintegral, NEQ, energy, x, xout, ITOL, RTOL, ATOL, 
     $        ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JCUSTOM, MF)
      ENDIF
      xfun = energy(1)+ifac*energy(2)
      END FUNCTION xfun

      SUBROUTINE Fxintegral(neq, x, y, ydot)
      INTEGER,DIMENSION(4) ::  NEQ
      REAL*8 x, Y(2+7), YDOT(2+7)
      INTEGER :: nutype,s,offset
      REAL(r8) :: wn,wt,we,wd,wb,l,ximag,nux,neo
      COMPLEX(r8) fx,gx,cx
      wn = Y(3)
      wt = Y(4)
      we = Y(5)
      wd = Y(6)
      wb = Y(7)
      l  = Y(8)
      ximag  = Y(9)
      nutype = NEQ(2)
      s      = NEQ(3)
      offset = NEQ(4)
      cx = x
      SELECT CASE (nutype)
         CASE (0)
            nux = 0.0
            cx = cx+ximag*ifac
         CASE (1)
            nux = 1E-5*kin%f(5)
         CASE (2)
            nux = kin%f(s+6)
         CASE (3)
            nux = kin%f(s+6)*(1+0.25*l*l)*x**(-1.5)/(2.0*kin%f(9))
         CASE DEFAULT
            nux = kin%f(s+6)*(1+0.25*l*l)*x**(-1.5)/(2.0*kin%f(9))
      END SELECT
      neo=0
      IF(neorot_flag) neo=1
      fx = -(we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
      IF(offset>0) fx = -(wt*((1-neo)*(x-2.5)+neo))*x**2.5*EXP(-x)
      gx = ifac*(l*wb*SQRT(cx)-nn*(we+wd*cx))-nux
      YDOT(1) = REAL(fx/gx)
      YDOT(2) = AIMAG(fx/gx)
      RETURN
      END SUBROUTINE Fxintegral

      SUBROUTINE  Jxintegral (NEQ, T, Y, ML, MU, PD, NRPD)
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y, PD(NRPD,2)
      PD(:,:) = 0
      RETURN
      END SUBROUTINE Jxintegral




c-----------------------------------------------------------------------
c     subprogram 1. natorq_gen.
c     Calculate NTV on a magnetic surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_lsode_gen(lagbpar,divxprp,nl,maskpsi,ptfac,
     $     ximag,xmax,wefac,wdfac,collision,outsurfs)

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
     $     wbbar,wdbar,djdjbar,bo,ls,dpsi,bmax,bmin,dt,bhat,dhat
      REAL(r8), DIMENSION(:), ALLOCATABLE :: bpts,wd0
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: htheta,tdt,ldl
      COMPLEX(r8) :: totals,xint
      COMPLEX(r8), DIMENSION(-nl:nl) :: lint
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
      PRINT *, "PENT - General calculation diagnostic"

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
         ! clebsch conversion *Exp(thetanq) now in djdt
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
               !action%xs(:) = tdt(0,:)
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
     $                 *EXP(-twopi*ifac*nn*sq%f(4)*tdt(0,itheta))
               ENDDO
               call spline_fit(omega,"extrap")
               CALL spline_int(omega)
               htheta(:,1)= omega%fsi(:,1)!-omega%fsi(omega%mx/2,1)
               htheta(:,2)= omega%fsi(:,2)!-omega%fsi(omega%mx/2,2)
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
                     iout = MAX(MIN(iter/(iters/10),10),1)
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
c$$$            lspl%fs(:,3) = 2.5*lspl%fs(:,2)+welec/(kin%f(s+2)/chrg)
c$$$            CALL spline_fit(lspl,'extrap')
c$$$            CALL spl_roots(wd0,lspl,3)
c$$$            CALL powergrid(lnorm,ldl,wd0,1,"both") !includes ends
c$$$            DEALLOCATE(wd0)
c$$$            DO ilmda=0,nlmda
c$$$               lmda = ldl(0,ilmda)
c$$$               CALL spline_eval(lspl,lmda,0)
c$$$               wbhat = SQRT(2*kin%f(s+2)/mass)*lspl%f(1)
c$$$               wdhat = (kin%f(s+2)/chrg)*lspl%f(2)
c$$$               DO ll=-nl,nl
c$$$                  ls = REAL(ll+sigma*nn*sq%f(4)) 
c$$$                  IF(offset_flag)THEN
c$$$                     xspl = xintegral(wdian,wdiat,welec,wdhat,wbhat,ls,
c$$$     $                    collision,s,ximag,xmax,.TRUE.)
c$$$                     offdl%fs(ilmda,nl+ll+1)=ldl(1,ilmda)*lspl%f(1)*
c$$$     $                 lspl%f(4+nl+ll)*xspl%fsi(xspl%mx,1)
c$$$                  ENDIF
c$$$                  xint = lsode_x(wdian,wdiat,welec,wdhat,wbhat,ls,
c$$$     $                    ximag,xmax,3,s,0,0)
c$$$                  xspl = xintegral(wdian,wdiat,welec,wdhat,wbhat,ls,
c$$$     $                 collision,s,ximag,xmax,.FALSE.)
c$$$                  ntvdl%fs(ilmda,nl+ll+1) = !ldl(1,ilmda)*
c$$$     $                 lspl%f(1)*lspl%f(4+nl+ll)*xint
c$$$                  ntvdl%xs(ilmda) = lmda!lnorm(ilmda)
c$$$c-----------------------------------------------------------------------
c$$$c     Record a sampling of intermediate steps for debuging
c$$$c-----------------------------------------------------------------------
c$$$                  IF(outsurf)THEN
c$$$                     iout = MAX(MIN(iter/(iters/10),10),1)
c$$$                     p_out(1,ilmda+1,iout,ll+nl+1,1:11) = 
c$$$     $                    (/psifac(istep),lnorm(ilmda),ldl(:,ilmda),
c$$$     $                    REAL(xint),
c$$$     $                    AIMAG(xint),
c$$$     $                    REAL(xspl%fsi(xspl%mx,1)),
c$$$     $                    AIMAG(xspl%fsi(xspl%mx,1)),lspl%f(1)*ro,
c$$$     $                    ro**2*bo*(lspl%f(2)),
c$$$     $                    lspl%f(4+nl+ll)*lspl%f(1)**2/ro**2/)
c$$$                     IF(ll==0 .AND. (MOD(ilmda,nlmda/4)==0))THEN
c$$$                        DO i = 0,nx,nx/99
c$$$                           ix = MIN(i/(nx/99)+1,100)
c$$$                           x_out(ix,ilmda/(nlmda/4)+1,iout,1,:)=(/ 
c$$$     $                          psifac(istep),lmda,xspl%xs(i),
c$$$     $                          REAL(xspl%fs(i,1)),AIMAG(xspl%fs(i,1)),
c$$$     $                          REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1))
c$$$     $                          ,lspl%f(1)*lspl%f(4+nl+ll) /)
c$$$                        ENDDO
c$$$                     ENDIF
c$$$                  ENDIF
c$$$               ENDDO ! l loop
c$$$            ENDDO ! lmda loop
c$$$c-----------------------------------------------------------------------
c$$$c     Final integration over Lambda and normalization to Nm
c$$$c-----------------------------------------------------------------------
c$$$            CALL cspline_fit(ntvdl,"extrap")
c$$$            CALL cspline_int(ntvdl)
c$$$            
c$$$            IF(outsurf)THEN
c$$$               DO ll=-nl,nl
c$$$                  DO i=0,lspl%mx
c$$$                     p_out(1,i+1,iout,ll+nl+1,12:13) = 
c$$$     $                    (/ REAL(ntvdl%fsi(i,ll+nl+1)),
c$$$     $                    AIMAG(ntvdl%fsi(i,ll+nl+1)) /)
c$$$                  ENDDO
c$$$               ENDDO
c$$$            ENDIF
c$$$            IF(sigma==0) ntv%fs(iter-1,:)=ntvdl%fsi(ntvdl%mx,:)*
c$$$     $           ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
c$$$            IF(sigma==1) ntvp%fs(iter-1,:)=ntvdl%fsi(ntvdl%mx,:)*
c$$$     $           ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
c-----------------------------------------------------------------------
c     lsode lambda integral to resolve wd~0 resonance in Ix
c-----------------------------------------------------------------------
            bhat = SQRT(2*kin%f(s+2)/mass)
            dhat = (kin%f(s+2)/chrg)
            CALL spline_fit(lspl,'extrap')
            DO ll=-nl,nl
               ls = REAL(ll+sigma*nn*sq%f(4))
               IF(outsurf) print *, psifac(istep)
               lint(ll) = lsode_lambda(wdian,wdiat,welec,bhat,dhat,
     $              lspl,4+nl+ll,ls,ximag,xmax,collision,s,
     $              .FALSE.,outsurf)
            ENDDO
            IF(sigma==0) ntv%fs(iter-1,:)=lint(:)*
     $           ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
            IF(sigma==1) ntvp%fs(iter-1,:)=lint(:)*
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
     $        "dLambda_n","real(lsode)",
     $        "imag(lsode)","real(I_x)","imag(I_x)","omega_b","omega_D",
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
      END SUBROUTINE diag_lsode_gen









c-----------------------------------------------------------------------
c     terminate module.
c-----------------------------------------------------------------------
      END MODULE diagnostic_mod
