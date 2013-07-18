c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     energy integration
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. energy_mod
c     A. lsode_x
c     B. lsode_lambda
c     C. intspl_x
c     1. noj
c     2. lsode_x_integrand
c     3. lsode_lambda_integrand
c-----------------------------------------------------------------------
c     subprogram 0. energy_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE energy_mod
      USE collision_mod
      USE pentglobal_mod
      USE grid_mod

      IMPLICIT NONE

      type(spline_type), PRIVATE :: lambda_spline

      CONTAINS

c-----------------------------------------------------------------------
c     Function A. lsode_x.
c     Wrapper for lsode energy integration
c     
c     args - wn     : density gradient diamagnetic drift frequency 
c          - wt     : temperature gradient diamagnetic drift frequency
c          - we     : electric precession frequency
c          - wd     : magnetic precession frequency
c          - wb     : bounce frequency
c          - l      : precession mode number (with -nq if passing)
c          - ximag  : step off of real axes (used if nutype "zero")
c          - xmax   : upper limit of integration
c          - nutype : collision operator
c          - s      : species (1 for ions 2 for electrons)
c          - offset : >0 to calculate offset rotation 
c          - out    : integrate step by step, recording as you go
c-----------------------------------------------------------------------
      FUNCTION lsode_x(wn,wt,we,wd,wb,l,ximag,xmax,nutype,s,offset,out)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: lsode_x
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax
      INTEGER, INTENT(IN) ::s
      LOGICAL, INTENT(IN) :: offset,out
      CHARACTER(16), INTENT(IN) :: nutype

      REAL*8  RBLOCK(7)
      INTEGER IBLOCK(1)
      LOGICAL LBLOCK(1)
      CHARACTER(16) CBLOCK
      COMMON /x_block/ RBLOCK,IBLOCK,LBLOCK,CBLOCK

      !EXTERNAL lsode_x_integrand
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, MF, i
      INTEGER, PARAMETER :: 
     $     NEQ = 2,               ! true number of equations
     $     LIW  = 20 + NEQ,       ! for MF 22 ! only uses 20 if MF 10
     $     LRW  = 20 + 16*NEQ     ! for MF 10 
c    $     LRW = 22+9*NEQ+NEQ**2 ! for MF 22
      INTEGER IWORKX(LIW)
      REAL*8  ATOL, RTOL, RWORKX(LRW), X, XOUT, Y(NEQ)
c-----------------------------------------------------------------------
c     set lsode options - see lsode package for documentation
c-----------------------------------------------------------------------
      Y(1:2) = (/ 0,0 /)
      X = 1e-9
      XOUT = xmax
      RWORKX(1) = xmax
      ITOL = 1                  ! RTOL and ATOL are scalars
      RTOL = 1.D-9              ! 14
      ATOL = 1.D-6              ! 15
      ISTATE = 1                ! first step
      IOPT = 0                  ! no optional inputs
      MF = 10                   ! not stiff with unknown J 
c      MF = 22                   ! stiff with unknown J 
c-----------------------------------------------------------------------
c     Share all the relavent variables in a common block with integrand.
c-----------------------------------------------------------------------
      RBLOCK = (/ wn,wt,we,wd,wb,l,ximag /)
      IBLOCK = s
      LBLOCK = offset
      CBLOCK = nutype
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      IF(out) THEN
         ITASK = 2              ! single step
         OPEN(UNIT=out_unit,FILE="pent_energy_lsode_n"//
     $        TRIM(sn)//".out",STATUS="UNKNOWN",POSITION="APPEND")
         DO WHILE (x<xout)
            CALL LSODE(lsode_x_integrand, NEQ, Y, X, XOUT,   
     $           ITOL, RTOL, ATOL,ITASK,ISTATE, IOPT, RWORKX, LRW,  
     $           IWORKX, LIW, NOJ, MF)
            WRITE(out_unit,'(1x,3(1x,es16.8E3))') X,Y(1:2)
         ENDDO
         CALL ascii_close(out_unit)
      ELSE
         ITASK = 1              ! full integral
         CALL LSODE(lsode_x_integrand, NEQ, Y, X, XOUT, ITOL,  
     $        RTOL,ATOL,ITASK,ISTATE, IOPT, RWORKX, LRW, IWORKX, LIW,  
     $        NOJ, MF)
      ENDIF
c-----------------------------------------------------------------------
c     Convert to complex space solution
c-----------------------------------------------------------------------
      lsode_x = Y(1)+ifac*Y(2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lsode_x

c-----------------------------------------------------------------------
c     Function A. lsode_lambda.
c     Wrapper for lsode pitch integration
c     
c     args - wn     : density gradient diamagnetic drift frequency 
c          - wt     : temperature gradient diamagnetic drift frequency
c          - we     : electric precession frequency
c          - lspl   : spline of wb,wd,and djdj as functions of lambda
c          - wb     : bounce frequency
c          - l      : precession mode number (with -nq if passing)
c          - ximag  : step off of real axes (used if nutype "zero")
c          - xmax   : upper limit of integration
c          - nutype : collision operator
c          - s      : species (1 for ions 2 for electrons)
c          - offset : >0 to calculate offset rotation 
c          - out    : integrate step by step, recording as you go
c-----------------------------------------------------------------------
      FUNCTION lsode_lambda(wn,wt,we,bhat,dhat,lspl,lh,ls,ximag,
     $     xmax,nutype,s,offset,out)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: lsode_lambda
      REAL(r8), INTENT(IN) :: wn,wt,we,bhat,dhat,ls,ximag,xmax
      INTEGER, INTENT(IN) :: s,lh
      LOGICAL, INTENT(IN) :: offset,out
      CHARACTER(16), INTENT(IN) :: nutype
      TYPE(spline_type), INTENT(INOUT) :: lspl

      REAL*8  RBLOCK(8)
      INTEGER IBLOCK(2)
      CHARACTER(16) CBLOCK
      LOGICAL LBLOCK(2)
      COMMON /lmda_block/ RBLOCK,IBLOCK,LBLOCK,CBLOCK

      !EXTERNAL lsode_lambda_integrand
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, MF, i
      INTEGER, PARAMETER :: 
     $     NEQ = 2,               ! true number of equations
     $     LIW  = 20 + NEQ,       ! for MF 22 ! only uses 20 if MF 10
     $     LRW  = 20 + 16*NEQ     ! for MF 10 
c    $     LRW = 22+9*NEQ+NEQ**2 ! for MF 22
      INTEGER IWORK(LIW)
      REAL*8  ATOL, RTOL, RWORK(LRW), LMDA, LMDAOUT, Y(NEQ),CRIT

      REAL*8 RLS
      INTEGER ILS
      COMMON /DLS001/ RLS(218),ILS(37)
c-----------------------------------------------------------------------
c     set lsode options - see lsode package for documentation
c-----------------------------------------------------------------------
      Y(1:2) = (/ 0,0 /)
      LMDA    = lspl%xs(0)
      LMDAOUT = lspl%xs(lspl%mx)
      CRIT    = lspl%xs(lspl%mx)
      RWORK(1) = CRIT
      ITOL = 1                  ! RTOL and ATOL are scalars
      RTOL = 1.D-9              ! 14
      ATOL = 1.D-9              ! 15
      ISTATE = 1                ! first step
      IOPT = 0                  ! no optional inputs
      MF = 10                   ! not stiff with unknown J 
c      MF = 22                   ! stiff with unknown J 
c-----------------------------------------------------------------------
c     Share all the relavent variables in a common block with integrand.
c-----------------------------------------------------------------------
      RBLOCK = (/ wn,wt,we,dhat,bhat,ls,ximag,xmax /)
      IBLOCK = (/ s,lh /)
      CBLOCK = nutype
      LBLOCK = (/ offset,out /)
      lambda_spline = lspl
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      IF(out) THEN
         ITASK = 5              ! single step without passing rwork(1)
         !PRINT *,"Starting lmda integral lsode"
         print *,"Lambda range ",LMDA,CRIT
         PRINT *,"lmda ",lmda,lmdaout
         PRINT *,"y    ",Y
         PRINT *,"inputs ",itol,itask,iopt,lrw,liw,mf,rtol,atol
         PRINT *,"step ",IWORK(11), IWORK(12)
         WRITE(out_unit,'(1x,3(1x,es16.8E3))') LMDA,Y(1:2)
         OPEN(UNIT=out_unit,FILE="pent_pitch_lsode_n"//TRIM(sn)//".out",
     $        STATUS="UNKNOWN",POSITION="APPEND")
         DO WHILE (LMDA<LMDAOUT)
            CALL LSODE(lsode_lambda_integrand, NEQ, Y, LMDA, LMDAOUT,   
     $           ITOL,RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK,  
     $           LIW, NOJ, MF)
            PRINT *,"lmda ",lmda,lmdaout
            PRINT *,"y    ",Y
            PRINT *,"inputs ",itol,itask,iopt,lrw,liw,mf,rtol,atol
            PRINT *,"step ",IWORK(11), IWORK(12)
            WRITE(out_unit,'(1x,3(1x,es16.8E3))') LMDA,Y(1:2)
         ENDDO
         CALL ascii_close(out_unit)
      ELSE
         ITASK = 4              ! full int without overshooting crit
         CALL LSODE(lsode_lambda_integrand, NEQ, Y, LMDA, LMDAOUT, ITOL,  
     $        RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, 
     $        NOJ, MF)
      ENDIF
      lsode_lambda = Y(1)+ifac*Y(2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lsode_lambda


c-----------------------------------------------------------------------
c     Function C. intspl_x
c     Energy integral spline formation
c-----------------------------------------------------------------------
      FUNCTION intspl_x(wn,wt,we,wd,wb,l,ximag,xmax,nutype,s,offset,out)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) intspl_x
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax
      INTEGER, INTENT(IN) :: s
      LOGICAL, INTENT(IN) :: offset,out
      CHARACTER(16), INTENT(IN) :: nutype

      REAL(r8) :: quadb,quadc,xmin,ccntr
      REAL(r8), DIMENSION(1:2) :: sqrtxres,xres
      INTEGER  :: i,nres,neo
      REAL(r8), DIMENSION(0:nx) :: nux
      REAL(r8), DIMENSION(0:1,0:nx) :: xdx
      COMPLEX(r8), DIMENSION(0:nx) :: cx,fx,gx
      type(cspline_type) :: xspl
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
      xmin = 0.0
      xres = sqrtxres**2
      IF(kolim_flag)THEN
         CALL powergrid(xnorm,xdx,(/xmin,xmax/),3,"lower")
      ELSEIF(xres(1)==0)THEN
         CALL powergrid(xnorm,xdx,(/xmin,xmax/),5,"lower")
      ELSE
         CALL powergrid(xnorm,xdx,(/xmin,xres(1:nres),xmax/),5,"none")
      ENDIF
c-----------------------------------------------------------------------
c     get collisionality according to chosen operator 
c     -> integrate over complex contour if collisionality is zero
c-----------------------------------------------------------------------
      IF(nutype=="zero")THEN
         cx = xdx(0,:)+ifac*ximag
      ELSE
         cx = xdx(0,:)
      ENDIF
      nux = nu_all(nutype,s,l,xdx(0,:))
c-----------------------------------------------------------------------
c     form spline on virtual grid
c-----------------------------------------------------------------------
      CALL cspline_alloc(xspl,nx,1)
      xspl%xs(:) = xnorm

      IF(neorot_flag)THEN
         neo = 1
      ELSE
         neo = 0
      ENDIF
      IF(kolim_flag)THEN           ! unity resonant operator (we->inf)
         fx = cx**2.5*EXP(-cx)
         gx = 1.0
      ELSE
         IF(offset)THEN
            fx = -(wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
         ELSE
            fx = -(we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
         ENDIF
         gx = ifac*(l*wb*SQRT(cx)-nn*(we+wd*cx))-nux
      ENDIF
      xspl%fs(:,1) = xdx(1,:)*fx/gx
      CALL cspline_fit(xspl,"extrap")
      CALL cspline_int(xspl)
c-----------------------------------------------------------------------
c     optionally write intermediate steps to file
c-----------------------------------------------------------------------
      IF(out)THEN
         OPEN(UNIT=out_unit,FILE="pent_energy_n"//TRIM(sn)//".out",
     $        STATUS="UNKNOWN",POSITION="APPEND")
         DO i=0,nx
            WRITE(out_unit,'(1x,5(1x,es16.8E3))') xspl%xs(i),
     $           REAL(xspl%fs(i,1)),AIMAG(xspl%fs(i,1)),
     $           REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1))
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     complex solution
c-----------------------------------------------------------------------
      intspl_x = xspl%fsi(xspl%mx,1)
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION intspl_x























c-----------------------------------------------------------------------
c     subprogram 1. noj
c     Dummy jacobian for use in lsode when true value unknown
c     
c     See lsode module for details.
c-----------------------------------------------------------------------
      SUBROUTINE  noj (NEQ, T, Y, ML, MU, PD, NRPD)
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y, PD(NRPD,2)
      PD(:,:) = 0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE noj

c-----------------------------------------------------------------------
c     subprogram 2. lsode_x_integrand
c     Integrand of energy integral.
c
c     see lsode module for details.
c-----------------------------------------------------------------------
      SUBROUTINE lsode_x_integrand(NEQ, X, Y, YDOT)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER  NEQ
      REAL*8 X, Y(NEQ), YDOT(NEQ)

      COMPLEX(r8) fx,gx,cx
      REAL(r8) :: wn,wt,we,wd,wb,l,ximag,nux,neo
      INTEGER :: s
      LOGICAL :: offset
      CHARACTER(16) :: nutype

      COMMON /x_block/
     $     wn,wt,we,wd,wb,l,ximag, ! 7 reals
     $     s,                   ! 1 integers
     $     offset,              ! 1 logicals
     $     nutype               ! 1 16 character
c      PRINT *,wn,wt,we,wd,wb,l,ximag
c      PRINT *,s,offset,nutype
c-----------------------------------------------------------------------
c     get collisionality according to chosen operator 
c     -> integrate over complex contour if collisionality is zero
c-----------------------------------------------------------------------
      IF(nutype=="zero")THEN
         cx = x+ifac*ximag
      ELSE
         cx = x
      ENDIF
      nux = nu_x(nutype,s,l,x)
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      IF(neorot_flag)THEN
         neo = 1
      ELSE
         neo = 0
      ENDIF

      IF(kolim_flag)THEN           ! unity resonant operator (we->inf)
         fx = cx**2.5*EXP(-cx)
         gx = 1.0
      ELSE
         IF(offset)THEN
            fx = -(wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
         ELSE
            fx = -(we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
         ENDIF
         gx = ifac*(l*wb*SQRT(cx)-nn*(we+wd*cx))-nux
      ENDIF
c-----------------------------------------------------------------------
c     Convert to 2 real space solutions
c-----------------------------------------------------------------------
      YDOT(1) = REAL(fx/gx)
      YDOT(2) = AIMAG(fx/gx)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lsode_x_integrand


c-----------------------------------------------------------------------
c     subprogram 2. lsode_lambda_integrand
c     Integrand of pitch integral.
c
c     see lsode module for input details.
c
c     Additional inputs passed through common block with lsode_lambda
c     wrapper function rather than lsode's capability to store variables
c     in extra neq/y dimensions. This provides clarity and flexibility
c     in the type of variables shared.
c-----------------------------------------------------------------------
      SUBROUTINE lsode_lambda_integrand(NEQ, LMDA, Y, YDOT)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER  NEQ
      REAL*8 LMDA, Y(NEQ), YDOT(NEQ)

      INTEGER, DIMENSION(37) :: ISAV
      REAL(r8), DIMENSION(218) :: RSAV

      COMPLEX(r8) :: xint,xint2
      REAL(r8) :: wn,wt,we,bhat,dhat,wd,wb,ls,ximag,xmax,djdj
      INTEGER :: s,lindx
      LOGICAL :: off,out
      CHARACTER(16) :: nutype
      type(cspline_type) :: xspl

      COMMON /lmda_block/ 
     $     wn,wt,we,dhat,bhat,ls,ximag,xmax, ! 8 reals
     $     s,lindx,             ! 2 integers
     $     off,out,             ! 2 logicals
     $     nutype               ! 1 16 character

      REAL*8 RLS
      INTEGER ILS
      COMMON /DLS001/ RLS(218),ILS(37)
c-----------------------------------------------------------------------
c     Evaluate well behaved lambda functions.
c-----------------------------------------------------------------------
      CALL spline_eval(lambda_spline,lmda,0)
      wb = bhat*lambda_spline%f(1)
      wd = dhat*lambda_spline%f(2)
      djdj = lambda_spline%f(lindx)
c-----------------------------------------------------------------------
c     Evaluate Energy integral (possible wd~0 resonance)
c-----------------------------------------------------------------------
      IF(xlsode_flag)THEN
         CALL DSRCOM(RSAV,ISAV,1) ! save lsode common block
         xint = lsode_x(wn,wt,we,wd,wb,ls,ximag,xmax,nutype,s,off,out)!3,s,0,0)!nutype,s,off,out)
         CALL DSRCOM(RSAV,ISAV,2) ! restore lsode common block
         IF(out)THEN
            xint2 = intspl_x(wn,wt,we,wd,wb,ls,ximag,xmax,
     $           nutype,s,off,out)
            PRINT *,"xint ",xint,xint2
         ENDIF
      ELSE
         xint = intspl_x(wn,wt,we,wd,wb,ls,ximag,xmax,nutype,s,off,out)
      ENDIF
c-----------------------------------------------------------------------
c     Convert to 2 real space solutions
c-----------------------------------------------------------------------
      YDOT(1) = REAL(lambda_spline%f(1)*djdj*xint)
      YDOT(2) = AIMAG(lambda_spline%f(1)*djdj*xint)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lsode_lambda_integrand



c-----------------------------------------------------------------------
c     terminate module.
c-----------------------------------------------------------------------
      END MODULE energy_mod




