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
      SUBROUTINE diag_complexpow()
      COMPLEX(r8) :: x
      REAL(r8) :: y
      LOGICAL :: test = .FALSE.
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      x = 1+ifac
      PRINT *, "(1+i)**3 = ",x**3
      PRINT *, "(1+i)**-3 = ",x**(-3)
      PRINT *, "(1+i)**-3/2 = ",x**(-1.5)
      x = x*0
      y = 0
      if(x==0) test=.TRUE.
      PRINT *, "0**-3/2 = ",y**(-3/2)
      PRINT *, "(0+0i)**-3/2 = ",x**(-3/2)
      PRINT *, "x = 0 is ",test
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE diag_complexpow

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
      IMPLICIT NONE 
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
c      PRINT *,"*************** - lsode example - ***************"
      DO 40 IOUT = 1,1!12
         CALL LSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, 
     *        ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
c        WRITE(6,20)  T, Y(1), Y(2), Y(3)
c  20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
c        IF (ISTATE .LT. 0)  GO TO 80
  40    TOUT = TOUT*10.
c      WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
c  60  FORMAT(/'EX No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
c      STOP
c  80  WRITE(6,90)  ISTATE
c  90  FORMAT(///' Error halt.. ISTATE =',I3)
c      STOP
      END SUBROUTINE diag_lsode_example

      SUBROUTINE  FEX (NEQ, T, Y, YDOT)
      IMPLICIT NONE 
      INTEGER  NEQ
      REAL*8  T, Y(3), YDOT(3)
      !PRINT *, "FEX"
      YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
      YDOT(3) = 3.D7*Y(2)*Y(2)
      YDOT(2) = -YDOT(1) - YDOT(3)
      RETURN
      END SUBROUTINE FEX

      SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
      IMPLICIT NONE 
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
      IMPLICIT NONE 
c      EXTERNAL  FEX, JEX
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     $     MF
      INTEGER, PARAMETER :: NEQ(1) =1
      REAL*8  ATOL(1), RTOL(1), RWORK(36), T, TOUT, Y(1)

c      common /CALLBLOCK/ 
c     $     NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
c     $     ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF
c     $     neq,itol,itask,istate,iopt,lrw,iwork,liw,mf, ! 28 integers
c     $     y,t,tout,rtol,atol,rwork ! 41 reals

c      NEQ = 1
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
      DO 40 IOUT = 1,1!0
        IF(IOUT==10) TOUT=1 
        PRINT *,rwork(11:14)!neq,itol,itask,istate,iopt,lrw,iwork,liw,mf
        CALL LSODE2(FCUSTOM, NEQ, Y, T, TOUT, ITOL, RTOL,ATOL,ITASK,
     *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JCUSTOM, MF)
        PRINT *,rwork(11:14)
        WRITE(6,20)  T, (/ Y(1), -2.*SQRT(1-T)+2 /)
  20    FORMAT(' At t =',D17.9,'   numeric,analytic =',2D14.6)
c        IF (ISTATE .LT. 0)  GO TO 80
  40    TOUT = 1.-1./(10.**IOUT)
      WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
  60  FORMAT(/'COSTUM No. steps =',i4,',  No. f-s =',i4,
     $     ',  No. J-s =',i4)
c      STOP
c  80  WRITE(6,90)  ISTATE
c  90  FORMAT(///' Error halt.. ISTATE =',I3)
      STOP
      END SUBROUTINE diag_lsode_custom

      RECURSIVE SUBROUTINE  FCUSTOM (NEQc, Tc, Yc, YDOTc)
      IMPLICIT NONE 
      INTEGER  NEQc
      REAL*8  Tc, Yc(1), YDOTc(1)
      INTEGER, DIMENSION(100) :: ISAV !37
      REAL(r8), DIMENSION(1000) :: RSAV  !218 
      INTEGER, DIMENSION(28) :: ISAVCALL
      REAL(r8), DIMENSION(41) :: RSAVCALL

      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     $     MF, NEQ
      REAL*8  ATOL(1), RTOL, RWORK(36), T, TOUT, Y(1)

      REAL*8 ATOL2(1),RTOL2,T2,TOUT2,Y2(1)
      common /CALLBLOCK/ 
     $     NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     $     ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF

      !CALL CSRCOM(RSAVCALL,ISAVCALL,1)
      PRINT *, "start ", Tc,Yc
c      CALL DSRCOM(RSAV,ISAV,1)  ! save lsode common block
      CALL diag_lsode_example()
      ! dummy call to re-intate non-dls001 vars
c      PRINT *,NEQ, Y, TOUT, TOUT, ITOL, RTOL, ATOL, ITASK
c      PRINT *,ISTATE, IOPT, RWORK(11:13), LRW, IWORK(11:13), LIW
c      CALL DSRCOM(RSAV,ISAV,2)
      !ISTATE = 2
c      TOUT2 = Tc
c      Y2 = Yc
c      CALL LSODE (FDUM, NEQ, Y2,TOUT,TOUT, ITOL, RTOL, ATOL, ITASK,
c     *     ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JCUSTOM, MF)
      PRINT *,"finish ", Tc,Yc
c      PRINT *,NEQ, Y, TOUT, TOUT, ITOL, RTOL, ATOL, ITASK
c      PRINT *,ISTATE, IOPT, RWORK(11:13), LRW, IWORK(11:13), LIW
      ! restore lsode common block
c      CALL DSRCOM(RSAV,ISAV,2) 
      !CALL CSRCOM(RSAVCALL,ISAVCALL,2)
      YDOTc(1) = 1.0/SQRT(1-Tc)
      RETURN
      END SUBROUTINE FCUSTOM

      SUBROUTINE  JCUSTOM (NEQ, T, Y, ML, MU, PD, NRPD)
      IMPLICIT NONE 
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y, PD(NRPD,1)
      PD(1,1) = 0
      RETURN
      END SUBROUTINE JCUSTOM

      SUBROUTINE  FDUM (NEQ, T, Y, YDOT)
      IMPLICIT NONE 
      INTEGER  NEQ
      REAL*8  T, Y(1), YDOT(1)
      YDOT(1) = 0
      RETURN
      END SUBROUTINE FDUM

c$$$      SUBROUTINE CSRCOM(RSAV,ISAV,JOB)
c$$$      INTEGER JOB
c$$$      INTEGER, PARAMETER :: LR = 36, LI = 20
c$$$      REAL(r8), DIMENSION(*) :: RSAV
c$$$      INTEGER, DIMENSION(*) :: ISAV
c$$$      REAL(r8), DIMENSION(LR) :: RCALL
c$$$      INTEGER, DIMENSION(LI) :: ICALL
c$$$      COMMON /CALLBLOCK/ ICALL,RCALL
c$$$
c$$$      IF(JOB==1)THEN
c$$$         PRINT *,RCALL(11:14),ICALL(11:13)
c$$$         ISAV(:LI) = ICALL(:LI)
c$$$         RSAV(:LR) = RCALL(:LR)
c$$$      ELSEIF(JOB==2)THEN
c$$$         ICALL(:LI) = ISAV(:LI)
c$$$         RCALL(:LR) = RSAV(:LR)
c$$$      ENDIF
c$$$      RETURN
c$$$      END SUBROUTINE CSRCOM


c     This test shows that y of the first loop
c     is changed by each call... eventually printing
c     a post value of 2 rather than 3 if it was unique
c     to the first call (this is a bug in embedded lsode's)
      SUBROUTINE DIAG_SUB1(x)
      integer x,y,x2
      y = 2
      x2 = x
      if(x<5)then
         x = x+y
         y = y+1
         print *,"pre ", x,y
         call diag_sub2(x)
         print *,"post ",x,y
      endif
      return
      end subroutine diag_sub1

      subroutine diag_sub2(x)
      integer x
      CALL diag_sub1(x)
      return
      end subroutine diag_sub2
      



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
      CHARACTER(32) :: collision = "zero"
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
      REAL(r8) :: wn,wt,we,wd,wb,l,ximag,neo
      COMPLEX(r8) fx,gx,cx,nux
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
      CHARACTER(32), INTENT(IN) :: nutype
      LOGICAL, INTENT(IN) :: offset

      REAL(r8) :: quadb,quadc,xmin,ccntr
      REAL(r8), DIMENSION(1:2) :: sqrtxres,xres
      INTEGER  :: i,nres,neo
      REAL(r8), DIMENSION(0:1,0:nx) :: xdx
      COMPLEX(r8), DIMENSION(0:nx) :: cx,fx,gx,nux

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
      nux = nu_all(nutype,s,1,l,cx)
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
c     subprogram 1. profile_gen.
c     Calculate NTV on a magnetic surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_singlesurf(lagbpar,divxprp,nl,psis,ptfac,ximag,
     $     xmax,wefac,wdfac,collision,outsurfs)

      INTEGER, INTENT(IN) :: nl
      REAL(r8), INTENT(IN) :: ptfac,ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10), INTENT(IN) :: outsurfs
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf), INTENT(IN) :: 
     $     lagbpar,divxprp
      CHARACTER(32) :: collision

      LOGICAL :: outsurf
      INTEGER :: istep,ilast,inext,itheta,i,ilmda,ll,sigma,
     $     iters,iter,indx,s,iout,ix,pass,xsigma
      REAL(r8):: chrg,mass,epsa,epsr,lmda,welec,wdian,wdiat,bhat,dhat,
     $     wbbar,wdbar,bo,ls,dpsi,bmax,bmin,dt,psis
      REAL(r8), DIMENSION(:), ALLOCATABLE :: bpts
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: htheta,tdt,ldl
      COMPLEX(r8) :: totals
      COMPLEX(r8), DIMENSION(-nl:nl) :: lmdaint,lmdaoff
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: pl

      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
      CHARACTER(32) :: outfile
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
      IF(divxprp_flag) xsigma=1
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------  
      PRINT *, "PENT - general calculation v1.0"

      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)

      epsa = SQRT(SUM(rzphi%fs(rzphi%mx,:,1))/rzphi%my)/ro

      CALL spline_alloc(tspl,mthsurf,5)
      CALL cspline_alloc(o1spl,mthsurf,2)
      CALL spline_alloc(vpar,mthsurf,1)
      tspl%xs(:) = theta(:)
      o1spl%xs(:) = theta(:)
      vpar%xs(:) = theta(:)

      CALL spline_alloc(omega,nt,2)
      CALL cspline_alloc(action,nt,2)
      ALLOCATE(htheta(0:nt,2),pl(0:nt),tdt(0:1,0:nt))
      omega%xs(:) = tnorm(:)
      action%xs(:) = tnorm(:)

      CALL spline_alloc(omegabar,nlmda,2)
      CALL spline_alloc(djdjbar,nlmda,2*nl+1)
      ALLOCATE(ldl(0:1,0:nlmda))

      outsurf = .TRUE.
! Find closest surface
      PRINT *,psis
      DO i=1,mstep
         ilast = MAX(1,i-1)
         inext = MIN(mstep,i+1)
         dpsi = ABS(psifac(i)-psis)
         IF(dpsi<=ABS(psifac(ilast)-psis) .AND.
     $        dpsi<=ABS(psifac(inext)-psis))THEN 
            PRINT *,"yes ", i,psifac(i)
            istep = i
            PRINT *,"istep ",istep
            psis = psifac(i)
            PRINT *, "psi = ",psis
         ENDIF
      ENDDO

      PRINT *,"psi,istep = ",psis,istep
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psis,1)
      CALL spline_eval(kin,psis,1)
      epsr = epsa*SQRT(psis)
!Rotation frequencies: electric and density and temp gradient
      welec = wefac*kin%f(5)
      wdian =-twopi*kin%f(s+2)*kin%f1(s)/(chrg*chi1*kin%f(s))
      wdiat =-twopi*kin%f1(s+2)/(chrg*chi1)
      
!Poloidal functions - note ABS(A*clebsch) = ABS(A)
      DO itheta=0,mthsurf
         CALL bicube_eval(eqfun,psis,theta(itheta),1)
         CALL bicube_eval(rzphi,psis,theta(itheta),1)
         tspl%fs(itheta,1)= eqfun%f(1) !b
         tspl%fs(itheta,2)= eqfun%fx(1)/chi1 !db/dpsi
         tspl%fs(itheta,3)= eqfun%fy(1) !db/dtheta
         tspl%fs(itheta,4)= rzphi%f(4)/chi1 !jac
         tspl%fs(itheta,5)= rzphi%fx(4)/chi1**2 !dj/dpsi
         o1spl%fs(itheta,1)= lagbpar(istep,itheta) !lagb
         o1spl%fs(itheta,2)= divxprp(istep,itheta) !divx
      ENDDO
! clebsch conversion now in djdt o1*EXP(-twopi*ifac*nn*q*theta)
      CALL spline_fit(tspl,"periodic")
      CALL cspline_fit(o1spl,"periodic")
      bmax = MAXVAL(tspl%fs(:,1),DIM=1)
      bmin = MINVAL(tspl%fs(:,1),DIM=1) 
      DO i=2,8                  !8th smallest so spline has more than 1 pt
         bmin = MINVAL(tspl%fs(:,1),MASK=tspl%fs(:,1)>bmin,DIM=1)
      ENDDO         
      DO sigma = 0,pass
c-----------------------------------------------------------------------
c     Calculate pitch dependent variables
c      Concentrate near trapped/passing boundary and deeply trapped lim
c-----------------------------------------------------------------------
         ldl = powspace(sigma*1e-6/bmax+
     $        (1-sigma)*(1+(1/bmin-1/bmax)*ptfac)/bmax
     $        ,sigma*(1-ptfac)/bmax+(1-sigma)/bmin,1,nlmda+1,"both")
         DO ilmda=0,nlmda
            lmda = ldl(0,ilmda)
c-----------------------------------------------------------------------
c     Determine bounce points
c-----------------------------------------------------------------------
            vpar%fs(:,1) = 1-lmda*tspl%fs(:,1)
            CALL spline_fit(vpar,"extrap")
            CALL spl_roots(bpts,vpar,1)
            CALL spline_eval(vpar,SUM(bpts(1:2))/2,0)
! the usual case with bpts centered around 0 requires shift
            IF(vpar%f(1)<=0.0 .AND. sigma==0) 
     $           bpts(1:2)=(/bpts(2)-1.0,bpts(1)/)
! hack - shouldn't happen (had repeated zeros in spl_roots)
            IF(bpts(1)>=bpts(2) .AND. sigma==0)THEN
               PRINT *," Error at (psi,lmda) = ",psis,
     $              lmda/ldl(0,nlmda)
               PRINT *,"Bounce pts = ",bpts
               CALL pentio_stop("Bounce point crossing")
            ENDIF
            CALL powergrid(tnorm,tdt,bpts(1-sigma:2-sigma),4,"both")
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
     $              /SQRT(vpar%f(1))
               omega%fs(itheta,2) = dt*tspl%f(4)*tspl%f(2)*
     $              (1-1.5*lmda*tspl%f(1))/SQRT(vpar%f(1))+
     $              dt*tspl%f(5)*tspl%f(1)*SQRT(vpar%f(1))
               action%fs(itheta,1) = dt*tspl%f(4)*tspl%f(1)*(
     $              xsigma*o1spl%f(2)*SQRT(vpar%f(1))+
     $              o1spl%f(1)*(1-1.5*lmda*tspl%f(1))/
     $              SQRT(vpar%f(1)))
     $              *EXP(-twopi*ifac*nn*sq%f(4)*(tdt(0,itheta)-
     $              tdt(0,1)))  ! doesn't matter since |dJ|^2
            ENDDO
            call spline_fit(omega,"extrap")
            CALL spline_int(omega)
            htheta(:,1)= omega%fsi(:,1)
            htheta(:,2)= omega%fsi(:,2) ! missing analytic piece?
            htheta(:,1)=htheta(:,1)/((2-sigma)*omega%fsi(omega%mx,1))
            htheta(:,2)=htheta(:,2)/((2-sigma)*omega%fsi(omega%mx,2))
! normalized bounce averaged frequencies
            wbbar = twopi/((2-sigma)*omega%fsi(omega%mx,1))
            wdbar = wdfac*wbbar*2*(2-sigma)*omega%fsi(omega%mx,2)
            omegabar%xs(ilmda)   = lmda
            omegabar%fs(ilmda,1) = wbbar
            omegabar%fs(ilmda,2) = wdbar
            djdjbar%xs(ilmda)   = lmda
            DO ll=-nl,nl
               ls = REAL(ll+sigma*nn*sq%f(4))
               IF(welec<SQRT(2*kin%f(s+2)/mass)*wdbar*3)THEN !mag precession dominates drift
                  pl(:) = EXP(-twopi*ifac*ls*htheta(:,2))
               ELSE             !electric precession dominates drift (normal)
                  pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
               ENDIF
               pl(:) = EXP(-twopi*ifac*ls*htheta(:,1))
               action%fs(:,2) = action%fs(:,1)*(pl(:)+1/pl(:))
               CALL cspline_fit(action,"extrap")
               CALL cspline_int(action)
               djdjbar%xs(ilmda) = lmda
               djdjbar%fs(ilmda,ll+nl+1) = 
     $              REAL(action%fsi(action%mx,2))**2 +
     $              AIMAG(action%fsi(action%mx,2))**2
c-----------------------------------------------------------------------
c     Optionally write bounce integration to file
c-----------------------------------------------------------------------
               IF(outsurf .AND. bounce_flag)THEN
                  OPEN(UNIT=out_unit,FILE="pent_bounce_n"//TRIM(sn)
     $                 //"_l"//TRIM(sl)//".out",STATUS="UNKNOWN",
     $                 POSITION="APPEND")
                  DO i=0,nx
                     WRITE(out_unit,'(1x,15(1x,es16.8E3))') kin%x,
     $                    lmda,omega%xs(i),tdt(:,i),
     $                    REAL(action%fs(i,2)),AIMAG(action%fs(i,2)),
     $                    REAL(action%fsi(i,2)),
     $                    AIMAG(action%fsi(i,2)),
     $                    omega%fs(i,:),omega%fsi(i,:),htheta(i,:)
                  ENDDO
                  CALL ascii_close(out_unit)
               ENDIF
            ENDDO               ! l loop
         ENDDO                  ! lmda loop
         CALL spline_fit(omegabar,'extrap')
         CALL spline_fit(djdjbar,'extrap')
            
c-----------------------------------------------------------------------
c     Lambda and E integrations
c-----------------------------------------------------------------------
! normalisation factors without energy dependence
         bhat = SQRT(2*kin%f(s+2)/mass)
         dhat = (kin%f(s+2)/chrg)
         DO ll=-nl,nl
            IF(lmdalsode_flag)THEN
               lmdaint(ll) =lsode_lambda(wdian,wdiat,welec,bhat,dhat,
     $              omegabar,djdjbar,ll,ximag,xmax,collision,0.0_r8,s,
     $              sigma,.FALSE.,outsurf,'diag_')
               IF(offset_flag) lmdaoff(ll) =lsode_lambda(wdian,
     $             wdiat,welec,bhat,dhat,omegabar,djdjbar,ll,ximag,xmax,
     $             collision,0.0_r8,s,sigma,.TRUE.,outsurf,'diag_')
            ELSE
               lmdaint(ll)=intspl_lambda(wdian,wdiat,welec,bhat,dhat,
     $              omegabar,djdjbar,ll,ximag,xmax,collision,0.0_r8,s,
     $              sigma,.FALSE.,outsurf,'diag_')
               IF(offset_flag) lmdaoff(ll) =lsode_lambda(wdian,
     $             wdiat,welec,bhat,dhat,omegabar,djdjbar,ll,ximag,xmax,
     $             collision,0.0_r8,s,sigma,.TRUE.,outsurf,'diag_')
            ENDIF
         ENDDO               
c-----------------------------------------------------------------------
c     Final profile in normalized flux (Nm)
c-----------------------------------------------------------------------
      ENDDO                     ! sigma loop
c-----------------------------------------------------------------------
c     Final integration over normalized flux
c-----------------------------------------------------------------------
      !Torque in Nm, Energy*R0 in Jm
      totals = SUM(lmdaint(:))*
     $        ro*chi1*nn**2*kin%f(s)*kin%f(s+2)/(twopi*SQRT(pi))
      PRINT *, "Total NTV torque = ",REAL(totals)
      PRINT *, "Total Kinetic Energy = ",AIMAG(totals)/(2*nn)
      PRINT *, "alpha/s  = ",REAL(totals)/(-1*AIMAG(totals))

      CALL spline_dealloc(tspl)
      CALL cspline_dealloc(o1spl)
      CALL spline_dealloc(vpar)

      CALL spline_dealloc(omega)
      CALL cspline_dealloc(action)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE diag_singlesurf





c-----------------------------------------------------------------------
c     subprogram 1. diag_elliptics.
c     Test NAG library elliptic functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE diag_elliptics()

      INTEGER :: istep
      REAL(r8):: k

c-----------------------------------------------------------------------
c     Write elliptics to a table
c-----------------------------------------------------------------------
      CALL pentio_read_special_functions
      PRINT *, "Writting Elliptic integrals to diag_elliptics.out"
      CALL ascii_open(out_unit,"diag_elliptics.out","UNKNOWN")
      WRITE(out_unit,*)"COMPLETE ELLIPTIC FUNCTIONS:"
      WRITE(out_unit,'(1/,1x,3(1x,a16))') "k","K(k)","E(k)"
      DO istep=1,5000
         k = istep/5001.0
         CALL spline_eval(ellipk,k,0)
         CALL spline_eval(ellipe,k,0)
         WRITE(out_unit,'(1x,15(1x,es16.8E3))')k,ellipk%f(1),ellipe%f(1)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE diag_elliptics


















c-----------------------------------------------------------------------
c     terminate module.
c-----------------------------------------------------------------------
      END MODULE diagnostic_mod
