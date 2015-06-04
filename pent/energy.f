c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     energy integration
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. energy_mod
c     The components of this module are organized into 4 blocks. Blocks
c     1 & 2 are lsode and spline methods of Lambda integration. Either
c     of these blocks use within them calls to either block 3 or 4,
c     where 3 & 4 represent lsode and spline energy integration 
c     respectively. 
c
c     BLOCK 1
c     1A. lsode_lambda               - Function wrapper for lsode
c     1B. lsode_lambda_integrand     - Subproutine integrand for lsode
c     1C. noj                        - Dummy subroutine for lsode
c     BLOCK 2
c     2A. intspl_lambda              - Function for spline integration
c
c     BLOCK 3
c     3A. lsode_x                    - Function wrapper for lsode
c     3B. lsode_x_integrand          - Subroutine integrand for lsode
c     BLOCK 4
c     3. intspl_x                    - Function for spline integration
c
c     BLOCK 5
c     5A. spl_roots                  - Subproutine root finding
c-----------------------------------------------------------------------
c     MODULE 0. energy_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE energy_mod
      USE collision_mod
      USE pentglobal_mod
      USE grid_mod
      USE lsode1_mod
      USE lsode2_mod

      IMPLICIT NONE

      type(spline_type), PRIVATE :: wbar_spline,jbar_spline

      CONTAINS


c-----------------------------------------------------------------------
c     BLOCK 1. Lambda integration using lsode
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     FUNCTION 1A. lsode_lambda.
c     Wrapper for lsode pitch integration
c     
c     args - wn     : density gradient diamagnetic drift frequency 
c          - wt     : temperature gradient diamagnetic drift frequency
c          - we     : electric precession frequency
c          - lspl   : spline of wb,wd,and djdj as functions of lambda
c          - wb     : bounce frequency
c          - l      : precession mode number (with -nq if passing)
c          - ximag  : step off of real axes
c          - xmax   : upper limit of integration
c          - nutype : collision operator
c          - s      : species (1 for ions 2 for electrons)
c          - offset : >0 to calculate offset rotation 
c          - out    : integrate step by step, recording as you go
c-----------------------------------------------------------------------
      FUNCTION lsode_lambda(wn,wt,we,bhat,dhat,wbar,jbar,ll,ximag,xmax,
     $     nutype,gamma,s,sigma,off,out,lbl)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: lsode_lambda
      REAL(r8), INTENT(IN) :: wn,wt,we,bhat,dhat,ximag,xmax,gamma
      INTEGER, INTENT(IN) :: s,sigma,ll
      LOGICAL, INTENT(IN) :: off,out
      CHARACTER(32), INTENT(IN) :: nutype,lbl

      TYPE(spline_type), INTENT(INOUT) :: wbar,jbar

      REAL*8  RBLOCK(8)
      INTEGER IBLOCK(3)
      CHARACTER(32) CBLOCK1,CBLOCK2
      LOGICAL LBLOCK(2)
      COMMON /lmda_block/ RBLOCK,IBLOCK,LBLOCK,CBLOCK1,CBLOCK2

      !EXTERNAL lsode_lambda_integrand
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, MF, IFLAG,i
      INTEGER, PARAMETER :: 
     $     NEQ = 2,               ! true number of equations
     $     LIW  = 20 + NEQ,       ! for MF 22 ! only uses 20 if MF 10
     $     LRW  = 20 + 16*NEQ     ! for MF 10 
c    $     LRW = 22+9*NEQ+NEQ**2 ! for MF 22
      INTEGER IWORK(LIW)
      REAL*8  ATOL, RTOL, RWORK(LRW), LMDA, LMDAOUT, Y(NEQ), DKY(NEQ),
     $     CRIT
c-----------------------------------------------------------------------
c     set lsode options - see lsode package for documentation
c-----------------------------------------------------------------------
      Y(1:2) = (/ 0,0 /)
      LMDA    = wbar%xs(0)
      LMDAOUT = wbar%xs(wbar%mx)
      CRIT    = wbar%xs(wbar%mx)
      ITOL = 1                  ! RTOL and ATOL are scalars
      RTOL = lsode_rtol         !1.D-12              ! 14
      ATOL = lsode_atol         !1.D-12              ! 15
      ISTATE = 1                ! first step
      IOPT = 1                  ! optional inputs
      IWORK(:) = 0
      RWORK(:) = 0
      RWORK(1) = CRIT
      IWORK(6) = 50000          ! max number of steps
      MF = 10                   ! not stiff with unknown J 
c      MF = 22                   ! stiff with unknown J 
c-----------------------------------------------------------------------
c     Share all the relavent variables in a common block with integrand.
c-----------------------------------------------------------------------
      RBLOCK = (/ wn,wt,we,bhat,dhat,ximag,xmax,gamma /)
      IBLOCK = (/ ll,s,sigma /)
      CBLOCK1 = nutype
      CBLOCK2 = lbl
      LBLOCK = (/ off,.FALSE. /) ! don't write function calls
      wbar_spline = wbar
      jbar_spline = jbar
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      IF(out .AND. pitch_flag) THEN
         ITASK = 5              ! single step without overshooting crit
         WRITE(UNIT=sl,FMT='(I2)') ll
         sl = ADJUSTL(sl)
         OPEN(UNIT=pitch_unit,FILE="pent_"//TRIM(lbl)//"pitch_n"//
     $        TRIM(sn)//"_l"//TRIM(sl)//".out",STATUS="UNKNOWN",
     $        POSITION="APPEND")
         DO WHILE (LMDA<LMDAOUT)
            CALL LSODE(lsode_lambda_integrand, NEQ, Y, LMDA, LMDAOUT,   
     $           ITOL,RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK,  
     $           LIW, NOJ, MF)
            CALL DINTDY (LMDA, 1, RWORK(21), NEQ, DKY, IFLAG)
            WRITE(pitch_unit,'(1x,9(1x,es16.8E3))') 
     $           kin%x,LMDA,DKY(1:2),Y(1:2),
     $           wbar_spline%f(1:2),jbar_spline%f(ll+nl+1)
         ENDDO
         CALL ascii_close(pitch_unit)
      ELSE
         ITASK = 4              ! full int without overshooting crit
         CALL LSODE(lsode_lambda_integrand, NEQ, Y, LMDA, LMDAOUT, ITOL,  
     $        RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, 
     $        NOJ, MF)
      ENDIF
      IF(IWORK(11)>IWORK(6)/2) PRINT *, "WARNING: ",IWORK(11)," of ",
     $     IWORK(6)," maximum steps used in Lambda integration."
      IF(ISTATE==-1) CALL pentio_stop("Too many steps in Lambda "
     $     //"required.")
c-----------------------------------------------------------------------
c     Repeat for optional output purposes only.
c-----------------------------------------------------------------------
      IF(out .AND. energy_flag)THEN
         WRITE(UNIT=sl,FMT='(I2)') ll
         sl = ADJUSTL(sl)
         DO i = 0,4
            lmda = wbar%xs(0) + (i/4.0)*(lmdaout-wbar%xs(0))
            LBLOCK(2) = .TRUE. ! re-run purely for file write
            CALL lsode_lambda_integrand(NEQ,LMDA,Y,DKY)
            LBLOCK(2) = .FALSE.
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complex solution.
c-----------------------------------------------------------------------
      lsode_lambda = Y(1)+ifac*Y(2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lsode_lambda

c-----------------------------------------------------------------------
c     SUBROUTINE 1B. lsode_lambda_integrand
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
      IMPLICIT NONE
      INTEGER  NEQ
      REAL*8 LMDA, Y(NEQ), YDOT(NEQ)

      INTEGER, DIMENSION(37) :: ISAV
      REAL(r8), DIMENSION(218) :: RSAV

      COMPLEX(r8) :: xint,test
      REAL(r8) :: wn,wt,we,bhat,dhat,wd,wb,ls,ximag,xmax,djdj,gamma
      INTEGER :: s,ll,sigma
      LOGICAL :: off,out
      CHARACTER(32) :: nutype,lbl
      type(cspline_type) :: xspl

      COMMON /lmda_block/ 
     $     wn,wt,we,bhat,dhat,ximag,xmax,gamma, ! 8 reals
     $     ll,s,sigma,          ! 3 integers
     $     off,out,             ! 2 logicals
     $     nutype,lbl           ! 2 32 character
c-----------------------------------------------------------------------
c     Evaluate well behaved lambda functions.
c-----------------------------------------------------------------------
      CALL spline_eval(wbar_spline,lmda,0)
      CALL spline_eval(jbar_spline,lmda,0)
      wb = bhat*wbar_spline%f(1)
      wd = dhat*wbar_spline%f(2)
      ls = REAL(ll+sigma*nn*sq%f(4))
      djdj = jbar_spline%f(ll+nl+1)
c-----------------------------------------------------------------------
c     Evaluate Energy integral (possible wd~0 resonance)
c-----------------------------------------------------------------------
      IF(xlsode_flag)THEN
         xint = lsode_x(wn,wt,we,wd,wb,ls,ximag,xmax,nutype,gamma,s,
     $        sigma,off,out,lbl)
         IF(sigma/=0) xint = xint 
     $        +lsode_x(wn,wt,we,wd,-wb,ls,ximag,xmax,nutype,gamma,s,
     $        sigma,off,out,lbl)
      ELSE
         xint = intspl_x(wn,wt,we,wd,wb,ls,ximag,xmax,nutype,gamma,s,
     $        sigma,off,out,lbl)
         IF(sigma/=0) xint = xint 
     $        +intspl_x(wn,wt,we,wd,-wb,ls,ximag,xmax,nutype,gamma,s,
     $        sigma,off,out,lbl)
      ENDIF
c-----------------------------------------------------------------------
c     Convert to 2 real space solutions
c-----------------------------------------------------------------------
      YDOT(1) = REAL(wbar_spline%f(1)*djdj*xint)
      YDOT(2) = AIMAG(wbar_spline%f(1)*djdj*xint)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lsode_lambda_integrand

c-----------------------------------------------------------------------
c     SUBROUTINE 1C. noj
c     Dummy jacobian for use in lsode when true value unknown
c     
c     See lsode module for details.
c-----------------------------------------------------------------------
      SUBROUTINE  noj (NEQ, T, Y, ML, MU, PD, NRPD)
      IMPLICIT NONE
      INTEGER  NEQ, ML, MU, NRPD
      REAL*8  T, Y, PD(NRPD,2)
      PD(:,:) = 0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE noj


c-----------------------------------------------------------------------
c     BLOCK 2. Lambda integration using splines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     FUNCTION 2A. intspl_lamda
c     Lambda (mu/E) integral spline formation
c
c     watch out for resonance in lmda (esp. low nu, low welec)
c     assume peak is when BHR in energy integrand peak, which is
c     1.5<~x<~2.5 for nu=0. Not exact -> low power concentration
c-----------------------------------------------------------------------
      FUNCTION intspl_lambda(wn,wt,we,bhat,dhat,wbar,jbar,ll,ximag,xmax,
     $     nutype,gamma,s,sigma,off,out,lbl)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: intspl_lambda
      INTEGER, INTENT(IN) :: s,ll,sigma
      REAL(r8), INTENT(IN) :: wn,wt,we,bhat,dhat,ximag,xmax,gamma
      LOGICAL, INTENT(IN) :: off,out
      CHARACTER(32), INTENT(IN) :: nutype,lbl
      type(spline_type) :: wbar,jbar,resspl
      
      type(cspline_type) :: lspl

      LOGICAL :: xout = .FALSE.
      INTEGER :: ilmda,i
      REAL(r8) :: wbhat,wdhat,lmda,ls
      REAL(r8), DIMENSION(:), ALLOCATABLE :: wd0,gridpts
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ldl
      COMPLEX(r8) :: xint

      wbar_spline = wbar
      jbar_spline = jbar

      IF(out)THEN
         WRITE(UNIT=sl,FMT='(I2)') ll
         sl = ADJUSTL(sl)
      ENDIF
c-----------------------------------------------------------------------
c     grid resolving approx location of lmda resonance.
c-----------------------------------------------------------------------
      ALLOCATE(ldl(0:1,0:nlmda))
      CALL spline_alloc(resspl,nlmda,1)
      resspl%xs(:) = wbar_spline%xs(:)
      resspl%fs(:,1) = 2.5*dhat*wbar_spline%fs(:,2)+we
      CALL spline_fit(resspl,'extrap')
      CALL spl_roots(wd0,resspl,1)
      ALLOCATE(gridpts(0:SIZE(wd0+1)))
      gridpts(1:SIZE(wd0)) = wd0(:)
      gridpts(0) = wbar_spline%xs(0)
      gridpts(SIZE(wd0)+1) = wbar_spline%xs(nlmda)
      CALL powergrid(lnorm,ldl,gridpts,2,"both") !includes ends
      DEALLOCATE(wd0,gridpts)
      CALL spline_dealloc(resspl)
c-----------------------------------------------------------------------
c     form lamda spline on new grid.
c-----------------------------------------------------------------------
      CALL cspline_alloc(lspl,nlmda,1)
      DO ilmda=0,nlmda
         lmda = ldl(0,ilmda)
         CALL spline_eval(wbar_spline,lmda,0)
         CALL spline_eval(jbar_spline,lmda,0)
         wbhat = bhat*wbar_spline%f(1)
         wdhat = dhat*wbar_spline%f(2)
         ls = REAL(ll+sigma*nn*sq%f(4))
c-----------------------------------------------------------------------
c     Evaluate Energy integral (possible wd~0 resonance)
c-----------------------------------------------------------------------
         IF(xlsode_flag)THEN
            xint = lsode_x(wn,wt,we,wdhat,wbhat,ls,ximag,xmax,
     $           nutype,gamma,s,sigma,off,xout,lbl)
            IF(sigma/=0) xint = xint 
     $           + lsode_x(wn,wt,we,wdhat,-wbhat,ls,ximag,xmax,
     $           nutype,gamma,s,sigma,off,xout,lbl)
         ELSE
            xint = intspl_x(wn,wt,we,wdhat,wbhat,ls,ximag,xmax,
     $           nutype,gamma,s,sigma,off,xout,lbl)
            IF(sigma/=0) xint = xint 
     $           + intspl_x(wn,wt,we,wdhat,-wbhat,ls,ximag,xmax,
     $           nutype,gamma,s,sigma,off,xout,lbl)
         ENDIF    
         lspl%xs(ilmda) = lmda
         lspl%fs(ilmda,1) = wbar_spline%f(1)*jbar_spline%f(nl+ll+1)*xint
      ENDDO
      CALL cspline_fit(lspl,'extrap')
      CALL cspline_int(lspl)
c-----------------------------------------------------------------------
c     optionally write intermediate steps to file
c-----------------------------------------------------------------------
      IF(out .AND. pitch_flag)THEN
         OPEN(UNIT=pitch_unit,FILE="pent_"//TRIM(lbl)//"pitch_n"//
     $        TRIM(sn)//"_l"//TRIM(sl)//".out",STATUS="UNKNOWN",
     $        POSITION="APPEND")
         DO ilmda=0,nlmda
            lmda = ldl(0,ilmda)
            CALL spline_eval(wbar_spline,lmda,0)
            CALL spline_eval(jbar_spline,lmda,0)
            WRITE(pitch_unit,'(1x,9(1x,es16.8E3))') kin%x,lmda,
     $           REAL(lspl%fs(ilmda,1)),AIMAG(lspl%fs(ilmda,1)),
     $           REAL(lspl%fsi(ilmda,1)),AIMAG(lspl%fsi(ilmda,1)),
     $           wbar_spline%f(1:2),jbar_spline%f(ll+nl+1)
         ENDDO
         CALL ascii_close(pitch_unit)
      ENDIF
c-----------------------------------------------------------------------
c     Repeat energy integrals for optional output purposes only.
c-----------------------------------------------------------------------
      IF(out .AND. energy_flag)THEN
         WRITE(UNIT=sl,FMT='(I2)') ll
         sl = ADJUSTL(sl)
         xout = .TRUE.
         DO i = 0,4
            lmda = wbar_spline%xs(0)+
     $           0.25*FLOAT(i)*(wbar_spline%xs(nlmda)-wbar_spline%xs(0))
            IF(xlsode_flag)THEN
               xint = lsode_x(wn,wt,we,wdhat,wbhat,ls,ximag,xmax,
     $              nutype,gamma,s,sigma,off,xout,lbl)
               IF(sigma/=0) xint = xint 
     $              + lsode_x(wn,wt,we,wdhat,-wbhat,ls,ximag,xmax,
     $              nutype,gamma,s,sigma,off,xout,lbl)
            ELSE
               xint = intspl_x(wn,wt,we,wdhat,wbhat,ls,ximag,xmax,
     $              nutype,gamma,s,sigma,off,xout,lbl)
               IF(sigma/=0) xint = xint 
     $              + intspl_x(wn,wt,we,wdhat,-wbhat,ls,ximag,xmax,
     $              nutype,gamma,s,sigma,off,xout,lbl)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complex solution
c-----------------------------------------------------------------------
      intspl_lambda = lspl%fsi(lspl%mx,1)
      CALL cspline_dealloc(lspl)
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION intspl_lambda





c-----------------------------------------------------------------------
c     BLOCK 3. x integration using lsode
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     FUNCTION 3A. lsode_x.
c     Wrapper for lsode energy integration
c
c     !!! NOTE this wrapper uses the local lsode2 mode. Due to the f90
c     compilation, (*) quantities may no longer be single int/floats, but
c     must explicitely have length (1) or higher. The use of lsode2 is
c     necessary to avoid memory overwrite bug in recursive lsode calls.
c
c     
c     args - wn     : density gradient diamagnetic drift frequency 
c          - wt     : temperature gradient diamagnetic drift frequency
c          - we     : electric precession frequency
c          - wd     : magnetic precession frequency
c          - wb     : bounce frequency
c          - l      : precession mode number (with -nq if passing)
c          - ximag  : step off of real axes
c          - xmax   : upper limit of integration
c          - nutype : collision operator
c          - s      : species (1 for ions 2 for electrons)
c          - offset : >0 to calculate offset rotation 
c          - out    : integrate step by step, recording as you go
c-----------------------------------------------------------------------
      FUNCTION lsode_x(wn,wt,we,wd,wb,l,ximag,xmax,nutype,gamma,
     $     s,sigma,off,out,lbl)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: lsode_x
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax,gamma
      INTEGER, INTENT(IN) ::s,sigma
      LOGICAL, INTENT(IN) :: off,out
      CHARACTER(32), INTENT(IN) :: nutype,lbl

      REAL*8  RBLOCK(8)
      INTEGER IBLOCK(2)
      LOGICAL LBLOCK(2)
      CHARACTER(32) CBLOCK
      COMMON /x_block/ RBLOCK,IBLOCK,LBLOCK,CBLOCK

      !EXTERNAL lsode_x_integrand
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, MF, IFLAG, i
      INTEGER, PARAMETER :: 
     $     NEQ(1) = 2,          ! true number of equations
     $     LIW  = 20 + NEQ(1),  ! for MF 22 ! only uses 20 if MF 10
     $     LRW  = 20 + 16*NEQ(1) ! for MF 10 
c    $     LRW = 22+9*NEQ(1)+NEQ(1)**2 ! for MF 22
      INTEGER IWORK(LIW)
      REAL*8  ATOL(1),RTOL(1),RWORK(LRW),X,XOUT,Y(NEQ(1)),DKY(NEQ(1))
c-----------------------------------------------------------------------
c     set lsode options - see lsode package for documentation
c-----------------------------------------------------------------------
      Y(1:2) = (/ 0,0 /)
      X = 1e-15
      XOUT = xmax
      ITOL = 2                  ! RTOL and ATOL are arrays
      RTOL = lsode_rtol         !1.D-7!9              ! 14
      ATOL = lsode_atol         !1.D-7!9              ! 15
      ISTATE = 1                ! first step
      IOPT = 1                  ! optional inputs
      IWORK(:) = 0              ! defaults
      RWORK(:) = 0              ! defaults
      RWORK(1) = xmax           ! only used if itask 4,5
      IWORK(6) = 10000          ! max number steps
      MF = 10                   ! not stiff with unknown J 
c      MF = 22                   ! stiff with unknown J 
c-----------------------------------------------------------------------
c     Share all the relavent variables in a common block with integrand.
c-----------------------------------------------------------------------
      RBLOCK = (/ wn,wt,we,wd,wb,l,ximag,gamma /)
      IBLOCK = (/ s,sigma /)
      LBLOCK = (/ off,.FALSE. /)
      CBLOCK = nutype
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      IF(out) THEN
         ITASK = 2              ! single step
         OPEN(UNIT=energy_unit,FILE="pent_"//TRIM(lbl)//
     $        "energy_n"//TRIM(sn)//"_l"//TRIM(sl)//".out",
     $        STATUS="UNKNOWN",POSITION="APPEND")
         DO WHILE (x<xout)
            CALL LSODE2(lsode_x_integrand, NEQ, Y, X, XOUT,   
     $           ITOL, RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW,  
     $           IWORK, LIW, NOJ, MF)
            CALL DINTDY2(X, 1, RWORK(21), NEQ(1), DKY, IFLAG)
            WRITE(energy_unit,'(1x,7(1x,es16.8E3))') kin%x,
     $           wbar_spline%x,X,DKY(1:2),Y(1:2)
         ENDDO
         CALL ascii_close(energy_unit)
      ELSE
         ITASK = 1              ! full integral
         CALL LSODE2(lsode_x_integrand, NEQ, Y, X, XOUT, ITOL,  
     $        RTOL,ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW,  
     $        NOJ, MF)
      ENDIF
c-----------------------------------------------------------------------
c     Additional output if integration failing
c-----------------------------------------------------------------------
      IF(IWORK(11)>IWORK(6)/2 .AND. ISTATE/=-1) PRINT *, "WARNING: ",
     $     IWORK(11)," of maximum ",IWORK(6)," steps in x integration"
      IF(ISTATE==-1) THEN
         OPEN(UNIT=energy_unit,FILE="pent_energy_stop_n"//TRIM(sn)//
     $        "_l"//TRIM(sl)//".out",STATUS="UNKNOWN")
         WRITE(out_unit,*) "psi = ",kin%x," Lambda = ",wbar_spline%x
         WRITE(energy_unit,'(5(1x,a16))') "x","T_phi",
     $           "2ndeltaW","int(T_phi)","int(2ndeltaW)"
         ITASK = 2
         Y(1:2) = (/ 0,0 /)
         X = 1e-15              !1e-12
         ISTATE = 1             ! first step
         IWORK(:) = 0           ! defaults
         RWORK(:) = 0           ! defaults
         RWORK(1) = xmax        ! only used if itask 4,5
         IWORK(6) = 10000       ! max number steps
         DO WHILE (x<xout .AND. IWORK(11)<IWORK(6))
            CALL LSODE2(lsode_x_integrand, NEQ, Y, X, XOUT,   
     $           ITOL, RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW,  
     $           IWORK, LIW, NOJ, MF)
            CALL DINTDY2(X, 1, RWORK(21), NEQ(1), DKY, IFLAG)
            WRITE(energy_unit,'(1x,5(1x,es16.8E3))') X,DKY(1:2),Y(1:2)
         ENDDO
         CALL ascii_close(energy_unit)
         CALL pentio_stop("Too many steps in x required."
     $     //" Consider complex contour (ximag>0).")
      ENDIF
c-----------------------------------------------------------------------
c     Convert to complex space solution
c-----------------------------------------------------------------------
      lsode_x = Y(1)+ifac*Y(2)
c-----------------------------------------------------------------------
c     Imaginary axis contribution
c-----------------------------------------------------------------------
      IF(ximag/=0.0)THEN
         LBLOCK(2) = .TRUE.
         Y(:) = 0
         X = 1e-15
         XOUT = ximag
         IWORK(:) = 0 
         IWORK(6) = 5000        ! max number of steps
         RWORK(:) = 0
         RWORK(1) = ximag       ! only relavent if using crit task
         ISTATE = 1
         IF(out) THEN
            ITASK = 2           ! single step
            OPEN(UNIT=ienergy_unit,FILE="pent_"//TRIM(lbl)//
     $           "ienergy_n"//TRIM(sn)//"_l"//TRIM(sl)//".out",
     $           STATUS="UNKNOWN",POSITION="APPEND")
            DO WHILE (x<xout)
               CALL LSODE2(lsode_x_integrand, NEQ, Y, X, XOUT,   
     $              ITOL, RTOL, ATOL,ITASK,ISTATE, IOPT, RWORK, LRW,  
     $              IWORK, LIW, NOJ, MF)
               CALL DINTDY2(X, 1, RWORK(21), NEQ(1), DKY, IFLAG)
               WRITE(ienergy_unit,'(1x,7(1x,es16.8E3))') kin%x,
     $              wbar_spline%x,X,DKY(1:2),Y(1:2)  
            ENDDO
         ELSE
            ITASK = 1           ! full integral
            CALL LSODE2(lsode_x_integrand, NEQ, Y, X, XOUT, ITOL,  
     $           RTOL,ATOL,ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW,  
     $           NOJ, MF)
         ENDIF
         IF(IWORK(11)>IWORK(6)/2) PRINT *, "WARNING: ",IWORK(11),
     $        " of ",IWORK(6)," maximum steps used in imaginary x "
     $        //"integration."
         IF(ISTATE==-1) CALL pentio_stop("Too many steps in " 
     $        //"imaginary x required. Change complex offset ximag.")
         lsode_x = lsode_x + Y(1)+ifac*Y(2)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lsode_x

c-----------------------------------------------------------------------
c     SUBROUTINE 3B. lsode_x_integrand
c     Integrand of energy integral.
c
c     see lsode module for details.
c-----------------------------------------------------------------------
      SUBROUTINE lsode_x_integrand(NEQ, X, Y, YDOT)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  NEQ
      REAL*8 X, Y(NEQ), YDOT(NEQ)

      COMPLEX(r8) :: fx,gx,cx,nux
      REAL(r8) :: wn,wt,we,wd,wb,l,ximag,neo,gamma
      INTEGER :: s,sigma
      LOGICAL :: offset,imaxis
      CHARACTER(32) :: nutype

      COMMON /x_block/
     $     wn,wt,we,wd,wb,l,ximag,gamma, ! 8 reals
     $     s,sigma,             ! 2 integers
     $     offset,imaxis,       ! 2 logicals
     $     nutype               ! 1 16 character
c-----------------------------------------------------------------------
c     get collisionality according to chosen operator 
c     -> integrate over complex contour if collisionality is zero
c-----------------------------------------------------------------------
      IF(imaxis)THEN
         cx = ifac*x
      ELSE
         cx = x + ifac*ximag
      ENDIF
      nux = nu_x(nutype,s,sigma,l,cx)
c-----------------------------------------------------------------------
c     calculations.
c-----------------------------------------------------------------------
      IF(neorot_flag)THEN
         neo = 1
      ELSE
         neo = 0
      ENDIF

      IF(kolim_flag)THEN           ! unity resonant operator (we->inf)
         fx = cx**2.5*EXP(-cx)/(ifac*nn)
         IF(offset) fx = 0
         gx = 1.0
      ELSE
         IF(offset)THEN
            fx = (wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
         ELSE
            fx = (we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
         ENDIF
         gx = ifac*(l*wb*SQRT(cx)+nn*(we+wd*cx))-nux-gamma
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
c     BLOCK 4. x integration using splines
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     FUNCTION 4A. intspl_x
c     Energy integral spline formation
c-----------------------------------------------------------------------
      FUNCTION intspl_x(wn,wt,we,wd,wb,l,ximag,xmax,nutype,gamma,
     $     s,sigma,off,out,lbl)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) intspl_x
      REAL(r8), INTENT(IN) :: wn,wt,we,wd,wb,l,ximag,xmax,gamma
      INTEGER, INTENT(IN) :: s,sigma
      LOGICAL, INTENT(IN) :: off,out
      CHARACTER(32), INTENT(IN) :: nutype,lbl

      REAL(r8) :: quadb,quadc,xmin,ccntr
      REAL(r8), DIMENSION(1:2) :: sqrtxres,xres
      INTEGER  :: i,nres,neo
c      REAL(r8), DIMENSION(0:nx) :: nux
      REAL(r8), DIMENSION(0:1,0:nx) :: xdx
      COMPLEX(r8), DIMENSION(0:nx) :: cx,fx,gx,nux
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
      cx = xdx(0,:) + ifac*ximag
c-----------------------------------------------------------------------
c     get collisionality according to chosen operator 
c     -> integrate over complex contour if collisionality is zero
c-----------------------------------------------------------------------
      nux = nu_all(nutype,s,sigma,l,cx)
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
         fx = cx**2.5*EXP(-cx)/(ifac*nn)
         IF(off) fx=0
         gx = 1.0
      ELSE
         IF(off)THEN
            fx = (wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
         ELSE
            fx = (we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
         ENDIF
         gx = ifac*(l*wb*SQRT(cx)+nn*(we+wd*cx))-nux-gamma
      ENDIF
      xspl%fs(:,1) = xdx(1,:)*fx/gx
      CALL cspline_fit(xspl,"extrap")
      CALL cspline_int(xspl)
c-----------------------------------------------------------------------
c     optionally write intermediate steps to file
c-----------------------------------------------------------------------
      IF(out)THEN
         OPEN(UNIT=energy_unit,FILE="pent_"//TRIM(lbl)//
     $        "energy_n"//TRIM(sn)//"_l"//TRIM(sl)//".out",
     $        STATUS="UNKNOWN",POSITION="APPEND")
         DO i=0,nx
            WRITE(energy_unit,'(1x,7(1x,es16.8E3))') kin%x,
     $           wbar_spline%x,xdx(0,i),REAL(xspl%fs(i,1)/xdx(1,i)),
     $           AIMAG(xspl%fs(i,1)/xdx(1,i)),
     $           REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1))
         ENDDO
         CALL ascii_close(energy_unit)
      ENDIF
c-----------------------------------------------------------------------
c     complex solution
c-----------------------------------------------------------------------
      intspl_x = xspl%fsi(xspl%mx,1)

c-----------------------------------------------------------------------
c     Imaginary axis integration
c-----------------------------------------------------------------------
      IF(ximag/=0)THEN
         CALL powergrid(xnorm,xdx,(/0.0_r8,ximag/),1,"both")
         cx = ifac*xdx(0,:)
         nux = nu_all(nutype,s,sigma,l,cx)
         IF(off)THEN
            fx = (wt*((1-neo)*(cx-2.5)+neo))*cx**2.5*EXP(-cx)
         ELSE
            fx = (we+wn+wt*((1-neo)*(cx-1.5)+neo*2))*cx**2.5*EXP(-cx)
         ENDIF
         gx = ifac*(l*wb*SQRT(cx)+nn*(we+wd*cx))-nux-gamma
         xspl%fs(:,1) = xdx(1,:)*fx/gx
         CALL cspline_fit(xspl,"extrap")
         CALL cspline_int(xspl)
c-----------------------------------------------------------------------
c     optionally write intermediate steps to file
c-----------------------------------------------------------------------
         IF(out)THEN
            OPEN(UNIT=ienergy_unit,FILE="pent_"//TRIM(lbl)//
     $           "ienergy_n"//TRIM(sn)//"_l"//TRIM(sl)//".out",
     $           STATUS="UNKNOWN",POSITION="APPEND")
            DO i=0,nx
               WRITE(ienergy_unit,'(1x,7(1x,es16.8E3))') kin%x,
     $              wbar_spline%x,xspl%xs(i),
     $              REAL(xspl%fs(i,1)),AIMAG(xspl%fs(i,1)),
     $              REAL(xspl%fsi(i,1)),AIMAG(xspl%fsi(i,1))
            ENDDO
            CALL ascii_close(ienergy_unit)
         ENDIF
c-----------------------------------------------------------------------
c     complex solution
c-----------------------------------------------------------------------
         intspl_x = intspl_x + xspl%fsi(xspl%mx,1)
      ENDIF
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION intspl_x












c-----------------------------------------------------------------------
c     BLOCK 5. root finding
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     SUBROUTINE 5A. spl_roots.
c     finds roots of a spline.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spl_roots(roots,spl,iqty)
      IMPLICIT NONE
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
c      IF(spl%fs(0,iqty)==0) nroots=nroots-1
      IF(spl%periodic .AND. spl%fs(spl%mx,iqty)==0) nroots=nroots-1
      IF(ALLOCATED(roots)) DEALLOCATE(roots)
      ALLOCATE(roots(1:nroots))
c      ALLOCATE(roots(0:nroots+1))
c      roots(0) = spl%xs(0)
c      roots(nroots+1) = spl%xs(spl%mx)
c-----------------------------------------------------------------------
c     find all zero passings, intialize at larger gradient.
c-----------------------------------------------------------------------
      DO ix=1,spl%mx
         ! dont calculate exact zeros twice
         IF(spl%fs(ix,iqty)==0.0 .AND. ix<spl%mx) CYCLE
         IF(spl%fs(ix,iqty)==0.0 .AND. spl%periodic) CYCLE
         ! find crossing window
         IF (spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .LE. 0.0) THEN
c            x=spl%xs(ix-1)-spl%fs(ix-1,iqty)*(spl%xs(ix)
c     $           -spl%xs(ix-1))/(spl%fs(ix,iqty)-spl%fs(ix-1,iqty))
            x = SUM(spl%xs(ix-1:ix))/2.0
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
c     terminate module.
c-----------------------------------------------------------------------
      END MODULE energy_mod




