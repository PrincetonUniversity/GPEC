c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     lsode integration copy of all subroutines that involve common
c     block to avoid recursive overwrites. Note that calls to 
c     DSRCOM were insufficient to avoid this.

c     !!! THIS MODULE STILL RELIES ON LSODE LIBRARY  for all subroutines
c     /functions that don't use the common block !!!!!
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. lsode_mod
c     1. lsode1
c-----------------------------------------------------------------------
c     MODULE 0. lsode_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE LSODE1_MOD
      USE local_mod, ONLY: r8
     

      IMPLICIT NONE

      REAL(r8) REAL_LSODE1
      INTEGER INT_LSODE1
      COMMON /DLS011/ REAL_LSODE1(218),INT_LSODE1(37)
!$OMP THREADPRIVATE(/DLS011/)

      CONTAINS

      SUBROUTINE LSODE1(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
     $     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
C***BEGIN PROLOGUE  LSODE
C***PURPOSE  Livermore solver for ordinary differential equations.
C            LSODE solves the initial-value problem for stiff or
C            nonstiff systems of first-order ODE's,
C               dy/dt = f(t,y),   or, in component form,
C               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
C***LIBRARY   MATHLIB (ODEPACK)
C***CATEGORY  I1A
C***TYPE      REAL(r8) (SLSODE-S, LSODE-D)
C***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
C             STIFF, NONSTIFF
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C             Computing and Mathematics Research Div., L-316
C             Lawrence Livermore National Laboratory
C             Livermore, CA 94550.
C***DESCRIPTION
C
C     NOTE: The LSODE solver is not re-entrant, and so is usable on
C           the Cray multi-processor machines only if it is not used
C           in a multi-tasking environment.
C           If re-entrancy is required, use NLSODE instead.
C
C           The formats of the LSODE and NLSODE writeups differ from
C           those of the other MATHLIB routines.
C
C           The "Usage" and "Arguments" sections treat only a subset of
C           available options, in condensed fashion.  The options
C           covered and the information supplied will support most
C           standard uses of LSODE.
C
C           For more sophisticated uses, full details on all options are
C           given in the concluding section, headed "Long Description."
C           A synopsis of the LSODE Long Description is provided at the
C           beginning of that section; general topics covered are:
C           - Elements of the call sequence; optional input and output
C           - Optional supplemental routines in the LSODE package
C           - internal COMMON block
C
C *Usage:
C     Communication between the user and the LSODE package, for normal
C     situations, is summarized here.  This summary describes a subset
C     of the available options.  See "Long Description" for complete
C     details, including optional communication, nonstandard options,
C     and instructions for special situations.
C
C     A sample program is given in the "Examples" section.
C
C     Refer to the argument descriptions for the definitions of the
C     quantities that appear in the following sample declarations.
C
C     For MF = 10,
C        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20)
C     For MF = 21 or 22,
C        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ)
C     For MF = 24 or 25,
C        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ,
C       *                                         LIW = 20 + NEQ)
C
C        EXTERNAL F, JAC
C        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW),
C       *         LIW, MF
C        REAL(r8) Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW)
C
C        CALL LSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
C       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
C
C *Arguments:
C     F     :EXT    Name of subroutine for right-hand-side vector f.
C                   This name must be declared EXTERNAL in calling
C                   program.  The form of F must be:
C
C                   SUBROUTINE  F (NEQ, T, Y, YDOT)
C                   INTEGER  NEQ
C                   REAL(r8)  T, Y(NEQ), YDOT(NEQ)
C
C                   The inputs are NEQ, T, Y.  F is to set
C
C                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)),
C                                                     i = 1, ..., NEQ .
C
C     NEQ   :IN     Number of first-order ODE's.
C
C     Y     :INOUT  Array of values of the y(t) vector, of length NEQ.
C                   Input:  For the first call, Y should contain the
C                           values of y(t) at t = T. (Y is an input
C                           variable only if ISTATE = 1.)
C                   Output: On return, Y will contain the values at the
C                           new t-value.
C
C     T     :INOUT  Value of the independent variable.  On return it
C                   will be the current value of t (normally TOUT).
C
C     TOUT  :IN     Next point where output is desired (.NE. T).
C
C     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or
C                   an array.
C
C     RTOL  :IN     Relative tolerance parameter (scalar).
C
C     ATOL  :IN     Absolute tolerance parameter (scalar or array).
C                   If ITOL = 1, ATOL need not be dimensioned.
C                   If ITOL = 2, ATOL must be dimensioned at least NEQ.
C
C                   The estimated local error in Y(i) will be controlled
C                   so as to be roughly less (in magnitude) than
C
C                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C
C                   Thus the local error test passes if, in each
C                   component, either the absolute error is less than
C                   ATOL (or ATOL(i)), or the relative error is less
C                   than RTOL.
C
C                   Use RTOL = 0.0 for pure absolute error control, and
C                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative
C                   error control.  Caution:  Actual (global) errors may
C                   exceed these local tolerances, so choose them
C                   conservatively.
C
C     ITASK :IN     Flag indicating the task LSODE is to perform.
C                   Use ITASK = 1 for normal computation of output
C                   values of y at t = TOUT.
C
C     ISTATE:INOUT  Index used for input and output to specify the state
C                   of the calculation.
C                   Input:
C                    1   This is the first call for a problem.
C                    2   This is a subsequent call.
C                   Output:
C                    2   LSODE was successful (otherwise, negative).
C                        Note that ISTATE need not be modified after a
C                        successful return.
C                   -1   Excess work done on this call (perhaps wrong
C                        MF).
C                   -2   Excess accuracy requested (tolerances too
C                        small).
C                   -3   Illegal input detected (see printed message).
C                   -4   Repeated error test failures (check all
C                        inputs).
C                   -5   Repeated convergence failures (perhaps bad
C                        Jacobian supplied or wrong choice of MF or
C                        tolerances).
C                   -6   Error weight became zero during problem
C                        (solution component i vanished, and ATOL or
C                        ATOL(i) = 0.).
C
C     IOPT  :IN     Flag indicating whether optional inputs are used:
C                   0   No.
C                   1   Yes.  (See "Optional inputs" under "Long
C                       Description," Part 1.)
C
C     RWORK :WORK   Real work array of length at least:
C                   20 + 16*NEQ                    for MF = 10,
C                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C
C     LRW   :IN     Declared length of RWORK (in user's DIMENSION
C                   statement).
C
C     IWORK :WORK   Integer work array of length at least:
C                   20        for MF = 10,
C                   20 + NEQ  for MF = 21, 22, 24, or 25.
C
C                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the
C                   lower and upper Jacobian half-bandwidths ML,MU.
C
C                   On return, IWORK contains information that may be
C                   of interest to the user:
C
C            Name   Location   Meaning
C            -----  ---------  -----------------------------------------
C            NST    IWORK(11)  Number of steps taken for the problem so
C                              far.
C            NFE    IWORK(12)  Number of f evaluations for the problem
C                              so far.
C            NJE    IWORK(13)  Number of Jacobian evaluations (and of
C                              matrix LU decompositions) for the problem
C                              so far.
C            NQU    IWORK(14)  Method order last used (successfully).
C            LENRW  IWORK(17)  Length of RWORK actually required.  This
C                              is defined on normal returns and on an
C                              illegal input return for insufficient
C                              storage.
C            LENIW  IWORK(18)  Length of IWORK actually required.  This
C                              is defined on normal returns and on an
C                              illegal input return for insufficient
C                              storage.
C
C     LIW   :IN     Declared length of IWORK (in user's DIMENSION
C                   statement).
C
C     JAC   :EXT    Name of subroutine for Jacobian matrix (MF =
C                   21 or 24).  If used, this name must be declared
C                   EXTERNAL in calling program.  If not used, pass a
C                   dummy name.  The form of JAC must be:
C
C                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C                   INTEGER  NEQ, ML, MU, NROWPD
C                   REAL(r8)  T, Y(NEQ), PD(NROWPD,NEQ)
C
C                   See item c, under "Description" below for more
C                   information about JAC.
C
C     MF    :IN     Method flag.  Standard values are:
C                   10  Nonstiff (Adams) method, no Jacobian used.
C                   21  Stiff (BDF) method, user-supplied full Jacobian.
C                   22  Stiff method, internally generated full
C                       Jacobian.
C                   24  Stiff method, user-supplied banded Jacobian.
C                   25  Stiff method, internally generated banded
C                       Jacobian.
C
C *Description:
C     LSODE solves the initial value problem for stiff or nonstiff
C     systems of first-order ODE's,
C
C        dy/dt = f(t,y) ,
C
C     or, in component form,
C
C        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ))
C                                                  (i = 1, ..., NEQ) .
C
C     LSODE is a package based on the GEAR and GEARB packages, and on
C     the October 23, 1978, version of the tentative ODEPACK user
C     interface standard, with minor modifications.
C
C     The steps in solving such a problem are as follows.
C
C     a. First write a subroutine of the form
C
C           SUBROUTINE  F (NEQ, T, Y, YDOT)
C           INTEGER  NEQ
C           REAL(r8)  T, Y(NEQ), YDOT(NEQ)
C
C        which supplies the vector function f by loading YDOT(i) with
C        f(i).
C
C     b. Next determine (or guess) whether or not the problem is stiff.
C        Stiffness occurs when the Jacobian matrix df/dy has an
C        eigenvalue whose real part is negative and large in magnitude
C        compared to the reciprocal of the t span of interest.  If the
C        problem is nonstiff, use method flag MF = 10.  If it is stiff,
C        there are four standard choices for MF, and LSODE requires the
C        Jacobian matrix in some form.  This matrix is regarded either
C        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the
C        banded case, LSODE requires two half-bandwidth parameters ML
C        and MU. These are, respectively, the widths of the lower and
C        upper parts of the band, excluding the main diagonal.  Thus the
C        band consists of the locations (i,j) with
C
C           i - ML <= j <= i + MU ,
C
C        and the full bandwidth is ML + MU + 1 .
C
C     c. If the problem is stiff, you are encouraged to supply the
C        Jacobian directly (MF = 21 or 24), but if this is not feasible,
C        LSODE will compute it internally by difference quotients (MF =
C        22 or 25).  If you are supplying the Jacobian, write a
C        subroutine of the form
C
C           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C           INTEGER  NEQ, ML, MU, NRWOPD
C           REAL(r8)  Y, Y(NEQ), PD(NROWPD,NEQ)
C
C        which provides df/dy by loading PD as follows:
C        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
C          the partial derivative of f(i) with respect to y(j).  (Ignore
C          the ML and MU arguments in this case.)
C        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
C          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the
C          rows of PD from the top down.
C        - In either case, only nonzero elements need be loaded.
C
C     d. Write a main program that calls subroutine LSODE once for each
C        point at which answers are desired.  This should also provide
C        for possible use of logical unit 6 for output of error messages
C        by LSODE.
C
C        Before the first call to LSODE, set ISTATE = 1, set Y and T to
C        the initial values, and set TOUT to the first output point.  To
C        continue the integration after a successful return, simply
C        reset TOUT and call LSODE again.  No other parameters need be
C        reset.
C
C *Examples:
C     The following is a simple example problem, with the coding needed
C     for its solution by LSODE. The problem is from chemical kinetics,
C     and consists of the following three rate equations:
C
C        dy1/dt = -.04*y1 + 1.E4*y2*y3
C        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
C        dy3/dt = 3.E7*y2**2
C
C     on the interval from t = 0.0 to t = 4.E10, with initial conditions
C     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
C
C     The following coding solves this problem with LSODE, using 
C     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses 
C     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2 
C     has much smaller values.  At the end of the run, statistical 
C     quantities of interest are printed.
C
C        EXTERNAL  FEX, JEX
C        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
C       *         MF, NEQ
C        REAL(r8)  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
C        NEQ = 3
C        Y(1) = 1.
C        Y(2) = 0.
C        Y(3) = 0.
C        T = 0.
C        TOUT = .4
C        ITOL = 2
C        RTOL = 1.D-4
C        ATOL(1) = 1.D-6
C        ATOL(2) = 1.D-10
C        ATOL(3) = 1.D-6
C        ITASK = 1
C        ISTATE = 1
C        IOPT = 0
C        LRW = 58
C        LIW = 23
C        MF = 21
C        DO 40 IOUT = 1,12
C          CALL LSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
C       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
C          WRITE(6,20)  T, Y(1), Y(2), Y(3)
C    20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
C          IF (ISTATE .LT. 0)  GO TO 80
C    40    TOUT = TOUT*10.
C        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
C    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
C        STOP
C    80  WRITE(6,90)  ISTATE
C    90  FORMAT(///' Error halt.. ISTATE =',I3)
C        STOP
C        END
C
C        SUBROUTINE  FEX (NEQ, T, Y, YDOT)
C        INTEGER  NEQ
C        REAL(r8)  T, Y(3), YDOT(3)
C        YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
C        YDOT(3) = 3.D7*Y(2)*Y(2)
C        YDOT(2) = -YDOT(1) - YDOT(3)
C        RETURN
C        END
C
C        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
C        INTEGER  NEQ, ML, MU, NRPD
C        REAL(r8)  T, Y(3), PD(NRPD,3)
C        PD(1,1) = -.04
C        PD(1,2) = 1.D4*Y(3)
C        PD(1,3) = 1.D4*Y(2)
C        PD(2,1) = .04
C        PD(2,3) = -PD(1,3)
C        PD(3,2) = 6.D7*Y(2)
C        PD(2,2) = -PD(1,2) - PD(3,2)
C        RETURN
C        END
C
C     The output from this program (on a Cray-1 in single precision)
C     is as follows.
C
C     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
C     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
C     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
C     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
C     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
C     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
C     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
C     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
C     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
C     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01
C     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
C     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00
C
C     No. steps = 330,  No. f-s = 405,  No. J-s = 69
C
C *Accuracy:
C     The accuracy of the solution depends on the choice of tolerances
C     RTOL and ATOL.  Actual (global) errors may exceed these local
C     tolerances, so choose them conservatively.
C
C *Cautions:
C     The work arrays should not be altered between calls to LSODE for
C     the same problem, except possibly for the conditional and optional
C     inputs.
C
C *Portability:
C     Since NEQ is dimensioned inside LSODE, some compilers may object
C     to a call to LSODE with NEQ a scalar variable.  In this event, 
C     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL.
C
C     Note to Cray users:
C     For maximum efficiency, use the CFT77 compiler.  Appropriate
C     compiler optimization directives have been inserted for CFT77 
C     (but not CIVIC).
C
C     NOTICE:  If moving the LSODE source code to other systems,
C     contact the author for notes on nonstandard Fortran usage,
C     COMMON block, and other installation details.
C
C *Reference:
C     Alan C. Hindmarsh, "ODEPACK, a systematized collection of ODE
C     solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
C     (North-Holland, Amsterdam, 1983), pp. 55-64.
C
C *Long Description:
C     The following complete description of the user interface to
C     LSODE consists of four parts:
C
C     1.  The call sequence to subroutine LSODE, which is a driver
C         routine for the solver.  This includes descriptions of both
C         the call sequence arguments and user-supplied routines.
C         Following these descriptions is a description of optional
C         inputs available through the call sequence, and then a
C         description of optional outputs in the work arrays.
C
C     2.  Descriptions of other routines in the LSODE package that may
C         be (optionally) called by the user.  These provide the ability
C         to alter error message handling, save and restore the internal
C         COMMON, and obtain specified derivatives of the solution y(t).
C
C     3.  Descriptions of COMMON block to be declared in overlay or
C         similar environments, or to be saved when doing an interrupt
C         of the problem and continued solution later.
C
C     4.  Description of two routines in the LSODE package, either of
C         which the user may replace with his own version, if desired.
C         These relate to the measurement of errors.
C
C
C                         Part 1.  Call Sequence
C                         ----------------------
C
C     Arguments
C     ---------
C     The call sequence parameters used for input only are
C
C        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
C
C     and those used for both input and output are
C
C        Y, T, ISTATE.
C
C     The work arrays RWORK and IWORK are also used for conditional and
C     optional inputs and optional outputs.  (The term output here
C     refers to the return from subroutine LSODE to the user's calling
C     program.)
C
C     The legality of input parameters will be thoroughly checked on the
C     initial call for the problem, but not checked thereafter unless a
C     change in input parameters is flagged by ISTATE = 3 on input.
C
C     The descriptions of the call arguments are as follows.
C
C     F        The name of the user-supplied subroutine defining the ODE
C              system.  The system must be put in the first-order form
C              dy/dt = f(t,y), where f is a vector-valued function of
C              the scalar t and the vector y. Subroutine F is to compute
C              the function f. It is to have the form
C
C                 SUBROUTINE F (NEQ, T, Y, YDOT)
C                 REAL(r8)  Y(NEQ), YDOT(NEQ)
C
C              where NEQ, T, and Y are input, and the array YDOT =
C              f(T,Y) is output.  Y and YDOT are arrays of length NEQ.
C              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be
C              declared EXTERNAL in the calling program.
C
C              Subroutine F may access user-defined quantities in
C              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array
C              (dimensioned in F) and/or Y has length exceeding NEQ(1).
C              See the descriptions of NEQ and Y below.
C
C              If quantities computed in the F routine are needed
C              externally to LSODE, an extra call to F should be made
C              for this purpose, for consistent and accurate results.
C              If only the derivative dy/dt is needed, use DINTDY1
C              instead.
C
C     NEQ      The size of the ODE system (number of first-order
C              ordinary differential equations).  Used only for input.
C              NEQ may be decreased, but not increased, during the
C              problem.  If NEQ is decreased (with ISTATE = 3 on input),
C              the remaining components of Y should be left undisturbed,
C              if these are to be accessed in F and/or JAC.
C
C              Normally, NEQ is a scalar, and it is generally referred
C              to as a scalar in this user interface description.
C              However, NEQ may be an array, with NEQ(1) set to the
C              system size.  (The LSODE package accesses only NEQ(1).)
C              In either case, this parameter is passed as the NEQ
C              argument in all calls to F and JAC.  Hence, if it is an
C              array, locations NEQ(2),... may be used to store other
C              integer data and pass it to F and/or JAC.  Subroutines
C              F and/or JAC must include NEQ in a DIMENSION statement
C              in that case.
C
C     Y        A real array for the vector of dependent variables, of
C              length NEQ or more.  Used for both input and output on
C              the first call (ISTATE = 1), and only for output on
C              other calls.  On the first call, Y must contain the
C              vector of initial values.  On output, Y contains the
C              computed solution vector, evaluated at T. If desired,
C              the Y array may be used for other purposes between
C              calls to the solver.
C
C              This array is passed as the Y argument in all calls to F
C              and JAC.  Hence its length may exceed NEQ, and locations
C              Y(NEQ+1),... may be used to store other real data and
C              pass it to F and/or JAC.  (The LSODE package accesses
C              only Y(1),...,Y(NEQ).)
C
C     T        The independent variable.  On input, T is used only on
C              the first call, as the initial point of the integration.
C              On output, after each call, T is the value at which a
C              computed solution Y is evaluated (usually the same as
C              TOUT).  On an error return, T is the farthest point
C              reached.
C
C     TOUT     The next value of T at which a computed solution is
C              desired.  Used only for input.
C
C              When starting the problem (ISTATE = 1), TOUT may be equal
C              to T for one call, then should not equal T for the next
C              call.  For the initial T, an input value of TOUT .NE. T
C              is used in order to determine the direction of the
C              integration (i.e., the algebraic sign of the step sizes)
C              and the rough scale of the problem.  Integration in
C              either direction (forward or backward in T) is permitted.
C
C              If ITASK = 2 or 5 (one-step modes), TOUT is ignored
C              after the first call (i.e., the first call with
C              TOUT .NE. T).  Otherwise, TOUT is required on every call.
C
C              If ITASK = 1, 3, or 4, the values of TOUT need not be
C              monotone, but a value of TOUT which backs up is limited
C              to the current internal T interval, whose endpoints are
C              TCUR - HU and TCUR.  (See "Optional Outputs" below for
C              TCUR and HU.)
C
C
C     ITOL     An indicator for the type of error control.  See
C              description below under ATOL.  Used only for input.
C
C     RTOL     A relative error tolerance parameter, either a scalar or
C              an array of length NEQ.  See description below under
C              ATOL.  Input only.
C
C     ATOL     An absolute error tolerance parameter, either a scalar or
C              an array of length NEQ.  Input only.
C
C              The input parameters ITOL, RTOL, and ATOL determine the
C              error control performed by the solver.  The solver will
C              control the vector e = (e(i)) of estimated local errors
C              in Y, according to an inequality of the form
C
C                 rms-norm of ( e(i)/EWT(i) ) <= 1,
C
C              where
C
C                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C
C              and the rms-norm (root-mean-square norm) here is
C
C                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ).
C
C              Here EWT = (EWT(i)) is a vector of weights which must
C              always be positive, and the values of RTOL and ATOL
C              should all be nonnegative.  The following table gives the
C              types (scalar/array) of RTOL and ATOL, and the
C              corresponding form of EWT(i).
C
C              ITOL    RTOL      ATOL      EWT(i)
C              ----    ------    ------    -----------------------------
C              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
C              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
C              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
C              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C              When either of these parameters is a scalar, it need not
C              be dimensioned in the user's calling program.
C
C              If none of the above choices (with ITOL, RTOL, and ATOL
C              fixed throughout the problem) is suitable, more general
C              error controls can be obtained by substituting
C              user-supplied routines for the setting of EWT and/or for
C              the norm calculation.  See Part 4 below.
C
C              If global errors are to be estimated by making a repeated
C              run on the same problem with smaller tolerances, then all
C              components of RTOL and ATOL (i.e., of EWT) should be
C              scaled down uniformly.
C
C     ITASK    An index specifying the task to be performed.  Input
C              only.  ITASK has the following values and meanings:
C              1   Normal computation of output values of y(t) at
C                  t = TOUT (by overshooting and interpolating).
C              2   Take one step only and return.
C              3   Stop at the first internal mesh point at or beyond
C                  t = TOUT and return.
C              4   Normal computation of output values of y(t) at
C                  t = TOUT but without overshooting t = TCRIT.  TCRIT
C                  must be input as RWORK(1).  TCRIT may be equal to or
C                  beyond TOUT, but not behind it in the direction of
C                  integration.  This option is useful if the problem
C                  has a singularity at or beyond t = TCRIT.
C              5   Take one step, without passing TCRIT, and return.
C                  TCRIT must be input as RWORK(1).
C
C              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C              (within roundoff), it will return T = TCRIT (exactly) to
C              indicate this (unless ITASK = 4 and TOUT comes before
C              TCRIT, in which case answers at T = TOUT are returned
C              first).
C
C     ISTATE   An index used for input and output to specify the state
C              of the calculation.
C
C              On input, the values of ISTATE are as follows:
C              1   This is the first call for the problem
C                  (initializations will be done).  See "Note" below.
C              2   This is not the first call, and the calculation is to
C                  continue normally, with no change in any input
C                  parameters except possibly TOUT and ITASK.  (If ITOL,
C                  RTOL, and/or ATOL are changed between calls with
C                  ISTATE = 2, the new values will be used but not
C                  tested for legality.)
C              3   This is not the first call, and the calculation is to
C                  continue normally, but with a change in input
C                  parameters other than TOUT and ITASK.  Changes are
C                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C                  ML, MU, and any of the optional inputs except H0.
C                  (See IWORK description for ML and MU.)
C
C              Note:  A preliminary call with TOUT = T is not counted as
C              a first call here, as no initialization or checking of
C              input is done.  (Such a call is sometimes useful for the
C              purpose of outputting the initial conditions.)  Thus the
C              first call for which TOUT .NE. T requires ISTATE = 1 on
C              input.
C
C              On output, ISTATE has the following values and meanings:
C               1  Nothing was done, as TOUT was equal to T with
C                  ISTATE = 1 on input.
C               2  The integration was performed successfully.
C              -1  An excessive amount of work (more than MXSTEP steps)
C                  was done on this call, before completing the
C                  requested task, but the integration was otherwise
C                  successful as far as T. (MXSTEP is an optional input
C                  and is normally 500.)  To continue, the user may
C                  simply reset ISTATE to a value >1 and call again (the
C                  excess work step counter will be reset to 0).  In
C                  addition, the user may increase MXSTEP to avoid this
C                  error return; see "Optional Inputs" below.
C              -2  Too much accuracy was requested for the precision of
C                  the machine being used.  This was detected before
C                  completing the requested task, but the integration
C                  was successful as far as T. To continue, the
C                  tolerance parameters must be reset, and ISTATE must
C                  be set to 3. The optional output TOLSF may be used
C                  for this purpose.  (Note:  If this condition is
C                  detected before taking any steps, then an illegal
C                  input return (ISTATE = -3) occurs instead.)
C              -3  Illegal input was detected, before taking any
C                  integration steps.  See written message for details.
C                  (Note:  If the solver detects an infinite loop of
C                  calls to the solver with illegal input, it will cause
C                  the run to stop.)
C              -4  There were repeated error-test failures on one
C                  attempted step, before completing the requested task,
C                  but the integration was successful as far as T.  The
C                  problem may have a singularity, or the input may be
C                  inappropriate.
C              -5  There were repeated convergence-test failures on one
C                  attempted step, before completing the requested task,
C                  but the integration was successful as far as T. This
C                  may be caused by an inaccurate Jacobian matrix, if
C                  one is being used.
C              -6  EWT(i) became zero for some i during the integration.
C                  Pure relative error control (ATOL(i)=0.0) was
C                  requested on a variable which has now vanished.  The
C                  integration was successful as far as T.
C
C              Note:  Since the normal output value of ISTATE is 2, it
C              does not need to be reset for normal continuation.  Also,
C              since a negative input value of ISTATE will be regarded
C              as illegal, a negative output value requires the user to
C              change it, and possibly other inputs, before calling the
C              solver again.
C
C     IOPT     An integer flag to specify whether any optional inputs
C              are being used on this call.  Input only.  The optional
C              inputs are listed under a separate heading below.
C              0   No optional inputs are being used.  Default values
C                  will be used in all cases.
C              1   One or more optional inputs are being used.
C
C     RWORK    A real working array (real(r8)).  The length of
C              RWORK must be at least
C
C                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
C
C              where
C                 NYH = the initial value of NEQ,
C              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                       smaller value is given as an optional input),
C                 LWM = 0           if MITER = 0,
C                 LWM = NEQ**2 + 2  if MITER = 1 or 2,
C                 LWM = NEQ + 2     if MITER = 3, and
C                 LWM = (2*ML + MU + 1)*NEQ + 2
C                                   if MITER = 4 or 5.
C              (See the MF description below for METH and MITER.)
C
C              Thus if MAXORD has its default value and NEQ is constant,
C              this length is:
C              20 + 16*NEQ                    for MF = 10,
C              22 + 16*NEQ + NEQ**2           for MF = 11 or 12,
C              22 + 17*NEQ                    for MF = 13,
C              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15,
C              20 +  9*NEQ                    for MF = 20,
C              22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C              22 + 10*NEQ                    for MF = 23,
C              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C
C              The first 20 words of RWORK are reserved for conditional
C              and optional inputs and optional outputs.
C
C              The following word in RWORK is a conditional input:
C              RWORK(1) = TCRIT, the critical value of t which the
C                         solver is not to overshoot.  Required if ITASK
C                         is 4 or 5, and ignored otherwise.  See ITASK.
C
C     LRW      The length of the array RWORK, as declared by the user.
C              (This will be checked by the solver.)
C
C     IWORK    An integer work array.  Its length must be at least
C              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
C              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25).
C              (See the MF description below for MITER.)  The first few
C              words of IWORK are used for conditional and optional
C              inputs and optional outputs.
C
C              The following two words in IWORK are conditional inputs:
C              IWORK(1) = ML   These are the lower and upper half-
C              IWORK(2) = MU   bandwidths, respectively, of the banded
C                              Jacobian, excluding the main diagonal.
C                         The band is defined by the matrix locations
C                         (i,j) with i - ML <= j <= i + MU. ML and MU
C                         must satisfy 0 <= ML,MU <= NEQ - 1. These are
C                         required if MITER is 4 or 5, and ignored
C                         otherwise.  ML and MU may in fact be the band
C                         parameters for a matrix to which df/dy is only
C                         approximately equal.
C
C     LIW      The length of the array IWORK, as declared by the user.
C              (This will be checked by the solver.)
C
C     Note:  The work arrays must not be altered between calls to LSODE
C     for the same problem, except possibly for the conditional and
C     optional inputs, and except for the last 3*NEQ words of RWORK.
C     The latter space is used for internal scratch space, and so is
C     available for use by the user outside LSODE between calls, if
C     desired (but not for use by F or JAC).
C
C     JAC      The name of the user-supplied routine (MITER = 1 or 4) to
C              compute the Jacobian matrix, df/dy, as a function of the
C              scalar t and the vector y.  (See the MF description below
C              for MITER.)  It is to have the form
C
C                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C                 REAL(r8)  Y(NEQ), PD(NROWPD,NEQ)
C
C              where NEQ, T, Y, ML, MU, and NROWPD are input and the
C              array PD is to be loaded with partial derivatives
C              (elements of the Jacobian matrix) on output.  PD must be
C              given a first dimension of NROWPD.  T and Y have the same
C              meaning as in subroutine F.
C
C              In the full matrix case (MITER = 1), ML and MU are
C              ignored, and the Jacobian is to be loaded into PD in
C              columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C
C              In the band matrix case (MITER = 4), the elements within
C              the band are to be loaded into PD in columnwise manner,
C              with diagonal lines of df/dy loaded into the rows of PD.
C              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML
C              and MU are the half-bandwidth parameters (see IWORK).
C              The locations in PD in the two triangular areas which
C              correspond to nonexistent matrix elements can be ignored
C              or loaded arbitrarily, as they are overwritten by LSODE.
C
C              JAC need not provide df/dy exactly. A crude approximation
C              (possibly with a smaller bandwidth) will do.
C
C              In either case, PD is preset to zero by the solver, so
C              that only the nonzero elements need be loaded by JAC.
C              Each call to JAC is preceded by a call to F with the same
C              arguments NEQ, T, and Y. Thus to gain some efficiency,
C              intermediate quantities shared by both calculations may
C              be saved in a user COMMON block by F and not recomputed
C              by JAC, if desired.  Also, JAC may alter the Y array, if
C              desired.  JAC must be declared EXTERNAL in the calling
C              program.
C
C              Subroutine JAC may access user-defined quantities in
C              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C              (dimensioned in JAC) and/or Y has length exceeding
C              NEQ(1).  See the descriptions of NEQ and Y above.
C
C     MF       The method flag.  Used only for input.  The legal values
C              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24,
C              and 25.  MF has decimal digits METH and MITER:
C                 MF = 10*METH + MITER .
C
C              METH indicates the basic linear multistep method:
C              1   Implicit Adams method.
C              2   Method based on backward differentiation formulas
C                  (BDF's).
C
C              MITER indicates the corrector iteration method:
C              0   Functional iteration (no Jacobian matrix is
C                  involved).
C              1   Chord iteration with a user-supplied full (NEQ by
C                  NEQ) Jacobian.
C              2   Chord iteration with an internally generated
C                  (difference quotient) full Jacobian (using NEQ
C                  extra calls to F per df/dy value).
C              3   Chord iteration with an internally generated
C                  diagonal Jacobian approximation (using one extra call
C                  to F per df/dy evaluation).
C              4   Chord iteration with a user-supplied banded Jacobian.
C              5   Chord iteration with an internally generated banded
C                  Jacobian (using ML + MU + 1 extra calls to F per
C                  df/dy evaluation).
C
C              If MITER = 1 or 4, the user must supply a subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For other values of MITER, a dummy argument can be used.
C
C     Optional Inputs
C     ---------------
C     The following is a list of the optional inputs provided for in the
C     call sequence.  (See also Part 2.)  For each such input variable,
C     this table lists its name as used in this documentation, its
C     location in the call sequence, its meaning, and the default value.
C     The use of any of these inputs requires IOPT = 1, and in that case
C     all of these inputs are examined.  A value of zero for any of
C     these optional inputs will cause the default value to be used.
C     Thus to use a subset of the optional inputs, simply preload
C     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively,
C     and then set those of interest to nonzero values.
C
C     Name    Location   Meaning and default value
C     ------  ---------  -----------------------------------------------
C     H0      RWORK(5)   Step size to be attempted on the first step.
C                        The default value is determined by the solver.
C     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
C                        default value is infinite.
C     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
C                        default value is 0.  (This lower bound is not
C                        enforced on the final step before reaching
C                        TCRIT when ITASK = 4 or 5.)
C     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
C                        is 12 if METH = 1, and 5 if METH = 2. (See the
C                        MF description above for METH.)  If MAXORD
C                        exceeds the default value, it will be reduced
C                        to the default value.  If MAXORD is changed
C                        during the problem, it may cause the current
C                        order to be reduced.
C     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
C                        allowed during one call to the solver.  The
C                        default value is 500.
C     MXHNIL  IWORK(7)   Maximum number of messages printed (per
C                        problem) warning that T + H = T on a step
C                        (H = step size).  This must be positive to
C                        result in a nondefault value.  The default
C                        value is 10.
C
C     Optional Outputs
C     ----------------
C     As optional additional output from LSODE, the variables listed
C     below are quantities related to the performance of LSODE which 
C     are available to the user.  These are communicated by way of the
C     work arrays, but also have internal mnemonic names as shown. 
C     Except where stated otherwise, all of these outputs are defined on
C     any successful return from LSODE, and on any return with ISTATE =
C     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3),
C     they will be unchanged from their existing values (if any), except
C     possibly for TOLSF, LENRW, and LENIW.  On any error return,
C     outputs relevant to the error will be defined, as noted below.
C
C     Name   Location   Meaning
C     -----  ---------  ------------------------------------------------
C     HU     RWORK(11)  Step size in t last used (successfully).
C     HCUR   RWORK(12)  Step size to be attempted on the next step.
C     TCUR   RWORK(13)  Current value of the independent variable which
C                       the solver has actually reached, i.e., the
C                       current internal mesh point in t. On output,
C                       TCUR will always be at least as far as the
C                       argument T, but may be farther (if interpolation
C                       was done).
C     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0,
C                       computed when a request for too much accuracy
C                       was detected (ISTATE = -3 if detected at the
C                       start of the problem, ISTATE = -2 otherwise).
C                       If ITOL is left unaltered but RTOL and ATOL are
C                       uniformly scaled up by a factor of TOLSF for the
C                       next call, then the solver is deemed likely to
C                       succeed.  (The user may also ignore TOLSF and
C                       alter the tolerance parameters in any other way
C                       appropriate.)
C     NST    IWORK(11)  Number of steps taken for the problem so far.
C     NFE    IWORK(12)  Number of F evaluations for the problem so far.
C     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU
C                       decompositions) for the problem so far.
C     NQU    IWORK(14)  Method order last used (successfully).
C     NQCUR  IWORK(15)  Order to be attempted on the next step.
C     IMXER  IWORK(16)  Index of the component of largest magnitude in
C                       the weighted local error vector ( e(i)/EWT(i) ),
C                       on an error return with ISTATE = -4 or -5.
C     LENRW  IWORK(17)  Length of RWORK actually required.  This is
C                       defined on normal returns and on an illegal
C                       input return for insufficient storage.
C     LENIW  IWORK(18)  Length of IWORK actually required.  This is
C                       defined on normal returns and on an illegal
C                       input return for insufficient storage.
C
C     The following two arrays are segments of the RWORK array which may
C     also be of interest to the user as optional outputs.  For each
C     array, the table below gives its internal name, its base address
C     in RWORK, and its description.
C
C     Name  Base address  Description
C     ----  ------------  ----------------------------------------------
C     YH    21            The Nordsieck history array, of size NYH by
C                         (NQCUR + 1), where NYH is the initial value of
C                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of
C                         YH contains HCUR**j/factorial(j) times the jth
C                         derivative of the interpolating polynomial
C                         currently representing the solution, evaluated
C                         at t = TCUR.
C     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated
C                         corrections on each step, scaled on output to
C                         represent the estimated local error in Y on
C                         the last step.  This is the vector e in the
C                         description of the error control.  It is
C                         defined only on successful return from LSODE.
C
C
C                    Part 2.  Other Callable Routines
C                    --------------------------------
C
C     The following are optional calls which the user may make to gain
C     additional capabilities in conjunction with LSODE.
C
C     Form of call              Function
C     ------------------------  ----------------------------------------
C     CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                               output of messages from LSODE, if the
C                               default is not desired.  The default
C                               value of LUN is 6. This call may be made
C                               at any time and will take effect
C                               immediately.
C     CALL XSETF(MFLAG)         Set a flag to control the printing of
C                               messages by LSODE.  MFLAG = 0 means do
C                               not print.  (Danger:  this risks losing
C                               valuable information.)  MFLAG = 1 means
C                               print (the default).  This call may be
C                               made at any time and will take effect
C                               immediately.
C     CALL DSRCOM1(RSAV,ISAV,JOB)  Saves and restores the contents of the
C                               internal COMMON blocks used by LSODE
C                               (see Part 3 below).  RSAV must be a
C                               real array of length 218 or more, and
C                               ISAV must be an integer array of length
C                               37 or more.  JOB = 1 means save COMMON
C                               into RSAV/ISAV.  JOB = 2 means restore
C                               COMMON from same.  DSRCOM1 is useful if
C                               one is interrupting a run and restarting
C                               later, or alternating between two or
C                               more problems solved with LSODE.
C     CALL DINTDY1(,,,,,)        Provide derivatives of y, of various
C     (see below)               orders, at a specified point t, if
C                               desired.  It may be called only after a
C                               successful return from LSODE.  Detailed
C                               instructions follow.
C
C     Detailed instructions for using DINTDY1
C     --------------------------------------
C     The form of the CALL is:
C
C           CALL DINTDY1 (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C     The input parameters are:
C
C     T          Value of independent variable where answers are
C                desired (normally the same as the T last returned by
C                LSODE).  For valid results, T must lie between
C                TCUR - HU and TCUR.  (See "Optional Outputs" above
C                for TCUR and HU.)
C     K          Integer order of the derivative desired.  K must
C                satisfy 0 <= K <= NQCUR, where NQCUR is the current
C                order (see "Optional Outputs").  The capability
C                corresponding to K = 0, i.e., computing y(t), is
C                already provided by LSODE directly.  Since
C                NQCUR >= 1, the first derivative dy/dt is always
C                available with DINTDY1.
C     RWORK(21)  The base address of the history array YH.
C     NYH        Column length of YH, equal to the initial value of NEQ.
C
C     The output parameters are:
C
C     DKY        Real array of length NEQ containing the computed value
C                of the Kth derivative of y(t).
C     IFLAG      Integer flag, returned as 0 if K and T were legal,
C                -1 if K was illegal, and -2 if T was illegal.
C                On an error return, a message is also written.
C
C
C                          Part 3.  Common Blocks
C                          ----------------------
C
C     If LSODE is to be used in an overlay situation, the user must
C     declare, in the primary overlay, the variables in:
C     (1) the call sequence to LSODE,
C     (2) the internal COMMON block /DLS011/, of length 255 
C         (218 real(r8) words followed by 37 integer words).
C
C     If LSODE is used on a system in which the contents of internal
C     COMMON blocks are not preserved between calls, the user should
C     declare the above COMMON block in his main program to insure that
C     its contents are preserved.
C
C     If the solution of a given problem by LSODE is to be interrupted
C     and then later continued, as when restarting an interrupted run or
C     alternating between two or more problems, the user should save,
C     following the return from the last LSODE call prior to the
C     interruption, the contents of the call sequence variables and the
C     internal COMMON block, and later restore these values before the
C     next LSODE call for that problem.   In addition, if XSETUN and/or
C     XSETF was called for non-default handling of error messages, then
C     these calls must be repeated.  To save and restore the COMMON
C     block, use subroutine DSRCOM1 (see Part 2 above).
C
C
C              Part 4.  Optionally Replaceable Solver Routines
C              -----------------------------------------------
C
C     Below are descriptions of two routines in the LSODE package which
C     relate to the measurement of errors.  Either routine can be
C     replaced by a user-supplied version, if desired.  However, since
C     such a replacement may have a major impact on performance, it
C     should be done only when absolutely necessary, and only with great
C     caution.  (Note:  The means by which the package version of a
C     routine is superseded by the user's version may be system-
C     dependent.)
C
C     DEWSET
C     ------
C     The following subroutine is called just before each internal
C     integration step, and sets the array of error weights, EWT, as
C     described under ITOL/RTOL/ATOL above:
C
C           SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C
C     where NEQ, ITOL, RTOL, and ATOL are as in the LSODE call
C     sequence, YCUR contains the current dependent variable vector,
C     and EWT is the array of weights set by DEWSET.
C
C     If the user supplies this subroutine, it must return in EWT(i)
C     (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C     in Y(i) to.  The EWT array returned by DEWSET is passed to the
C     DVNORM routine (see below), and also used by LSODE in the
C     computation of the optional output IMXER, the diagonal Jacobian
C     approximation, and the increments for difference quotient
C     Jacobians.
C
C     In the user-supplied version of DEWSET, it may be desirable to use
C     the current values of derivatives of y. Derivatives up to order NQ
C     are available from the history array YH, described above under
C     optional outputs.  In DEWSET, YH is identical to the YCUR array,
C     extended to NQ + 1 columns with a column length of NYH and scale
C     factors of H**j/factorial(j).  On the first call for the problem,
C     given by NST = 0, NQ is 1 and H is temporarily set to 1.0.  The
C     quantities NQ, NYH, H, and NST can be obtained by including in
C     DEWSET the statements:
C
C           REAL(r8)  RLS
C           COMMON /DLS011/ RLS(218),ILS(37)
C           NQ = ILS(33)
C           NYH = ILS(12)
C           NST = ILS(34)
C           H = RLS(212)
C
C     Thus, for example, the current value of dy/dt can be obtained as
C     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary
C     when NST = 0).
C
C     DVNORM
C     ------
C     DVNORM is a real function routine which computes the weighted
C     root-mean-square norm of a vector v:
C
C        d = DVNORM (n, v, w)
C
C     where:
C     n = the length of the vector,
C     v = real array of length n containing the vector,
C     w = real array of length n containing weights,
C     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ).
C
C     DVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where
C     EWT is as set by subroutine DEWSET.
C
C     If the user supplies this function, it should return a nonnegative
C     value of DVNORM suitable for use in the error control in LSODE.
C     None of the arguments should be altered by DVNORM.  For example, a
C     user-supplied DVNORM routine might:
C     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
C     - Ignore some components of v in the norm, with the effect of
C       suppressing the error control on those components of Y.
C  ---------------------------------------------------------------------
C***REFERENCES  Alan C. Hindmarsh, "ODEPACK, a systematized collection
C                 of ODE solvers", in Scientific Computing, R. S.
C                 Stepleman, et al. (Eds.), (North-Holland, Amsterdam,
C                 1983), pp. 55-64.
C***ROUTINES CALLED  DEWSET, DINTDY1, DUMACH, DSTODE1, DVNORM, XERRWD
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   870330  Major update by ACH.
C   890426  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890501  Many improvements to prologue.  (FNF)
C   890503  A few final corrections to prologue.  (FNF)
C   890504  Minor cosmetic changes.  (FNF)
C   890510  Corrected description of Y in Arguments section.  (FNF)
C   890517  Minor corrections to prologue.  (FNF)
C   920514  Updated with prologue edited 891025 by G. Shaw for manual.
C   920515  Converted source lines to upper case.  (FNF)
C   920603  Revised XERRWV calls using mixed upper-lower case.  (ACH)
C   920616  Revised prologue comment regarding CFT.  (ACH)
C   921116  Revised prologue comments regarding Common.  (ACH).
C   930326  Added comment about non-reentrancy.  (FNF)
C   930723  Changed R1MACH to RUMACH. (FNF)
C   930801  Removed Common variables ILLIN and NTREP (affects driver
C           logic and Common references); minor changes to prologue and
C           internal comments; changed Hollerith strings to quoted 
C           strings; changed internal comments to mixed case; changed
C           dummy dimensions from 1 to *. (ACH)
C   930809  Changed to generic intrinsic names; changed names of
C           subprograms and Common blocks to SLSODE etc. (ACH)
C   930929  Eliminated use of REAL intrinsic; other minor changes. (ACH)
C   931005  Generated real(r8) version. (ACH)
C***END PROLOGUE  LSODE
C
C*Internal Notes:
C
C Other Routines in the LSODE Package.
C
C In addition to Subroutine LSODE, the LSODE package includes the
C following subroutines and function routines:
C  DINTDY1   computes an interpolated value of the y vector at t = TOUT.
C  DSTODE1   is the core integrator, which does one step of the
C           integration and the associated error control.
C  DCFODE   sets all method coefficients and test constants.
C  DPREPJ1   computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  DSOLSY1   manages solution of linear system in chord iteration.
C  DEWSET   sets the error weight vector EWT before each step.
C  DVNORM   computes the weighted R.M.S. norm of a vector.
C  DSRCOM1   is a user-callable routine to save and restore
C           the contents of the internal Common block.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  DGBFA and DGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  DUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWD, XSETUN, and XSETF   handle the printing of all error
C           messages and warnings.  XERRWD is machine-dependent.
C Note.. DVNORM and DUMACH are function routines.  All the others 
C are subroutines.
C
C The intrinsic routines used by LSODE are..
C ABS, MAX, MIN, MOD, SIGN, and SQRT.
C
C**End
C
C  Declare arguments.
C
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL(r8) Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C
C  Declare externals.
C
c      EXTERNAL DPREPJ1, DSOLSY1
      REAL(r8) DUMACH, DVNORM
C
C  Declare all other variables.
C
      INTEGER INIT, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     1   MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     1   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      REAL(r8) ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL(r8) ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*80 MSG
C-----------------------------------------------------------------------
C The following internal common block contains
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine (own variables), and
C (b) variables which are communicated between subroutines.
C The structure of the block is as follows..  All real variables are
C listed first, followed by all integers.  Within each type, the
C variables are grouped with those local to Subroutine LSODE first,
C then those local to Subroutine DSTODE1, and finally those used
C for communication.  The block is declared in subroutines
C LSODE, DINTDY1, DSTODE1, DPREPJ1, and DSOLSY1.  Groups of variables are
C replaced by dummy arrays in the common declarations in routines
C where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /DLS011/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     3   MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!$OMP THREADPRIVATE(/DLS011/)

C
      DATA  MORD/12,5/
      DATA  MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .GT. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  LSODE
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0
      HMXI = 0.0
      HMIN = 0.0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0) GO TO 615
      HMXI = 0.0
      IF (HMAX .GT. 0.0) HMXI = 1.0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENWM = N*N + 2
      IF (MITER .EQ. 3) LENWM = N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 20
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0) GO TO 619
        IF (ATOLI .LT. 0.0) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODE1. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
 80     RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0) GO TO 625
      IF (H0 .NE. 0.0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0
      NQU = 0
      CCMAX = 0.3
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
C if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0*UROUND)
      TOL = MIN(TOL,0.001d0)
      SUM = DVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0) GO TO 250
      CALL DINTDY1 (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0 + 100.0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0) GO TO 245
      CALL DINTDY1 (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0 + 4.0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0) GO TO 250
      H = (TCRIT - TN)*(1.0 - 4.0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODE1.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0) GO TO 280
      TOLSF = TOLSF*2.0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'LSODE- Pitch -  Warning..internal T (=R1) and H (=R2) are '
      CALL XERRWD (MSG, 58, 101, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      (H = step size). Solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'LSODE- Pitch -  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 58, 102, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      It will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0, 0.0)
 290  CONTINUE
C-----------------------------------------------------------------------
C  CALL DSTODE1(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPREPJ1,DSOLSY1)
C-----------------------------------------------------------------------
      CALL DSTODE1 (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, DPREPJ1, DSOLSY1)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0) GO TO 250
      CALL DINTDY1 (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0) GO TO 345
      CALL DINTDY1 (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0 + 4.0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0) GO TO 250
      H = (TCRIT - TN)*(1.0 - 4.0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from LSODE.
C If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, the illegal input counter is zeroed, and the
C optional outputs are loaded into the work arrays before returning.
C If ISTATE = 1 and TOUT = T, there is a return with no action taken.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.  The optional outputs
C are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'LSODE- Pitch -  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 58, 201, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0)
      ISTATE = -1
      GO TO 580
C EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'LSODE- Pitch -  At T (=R1), EWT(I1) has become R2 .LE. 0.'
      CALL XERRWD (MSG, 58, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'LSODE- Pitch -  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 58, 203, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'LSODE- Pitch -  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 58, 204, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'LSODE- Pitch -  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 58, 205, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'LSODE- Pitch -  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 38, 1, 0, 1, ISTATE, 0, 0, 0.0, 0.0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'LSODE- Pitch -  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 38, 2, 0, 1, ITASK, 0, 0, 0.0, 0.0)
      GO TO 700
 603  MSG = 'LSODE- Pitch -  ISTATE .GT. 1 but LSODE not initialized '
      CALL XERRWD (MSG, 58, 3, 0, 0, 0, 0, 0, 0.0, 0.0)
      GO TO 700
 604  MSG = 'LSODE- Pitch -  NEQ (=I1) .LT. 1     '
      CALL XERRWD (MSG, 38, 4, 0, 1, NEQ(1), 0, 0, 0.0, 0.0)
      GO TO 700
 605  MSG = 'LSODE- Pitch -  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 58, 5, 0, 2, N, NEQ(1), 0, 0.0, 0.0)
      GO TO 700
 606  MSG = 'LSODE- Pitch -  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 38, 6, 0, 1, ITOL, 0, 0, 0.0, 0.0)
      GO TO 700
 607  MSG = 'LSODE- Pitch -  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 38, 7, 0, 1, IOPT, 0, 0, 0.0, 0.0)
      GO TO 700
 608  MSG = 'LSODE- Pitch -  MF (=I1) illegal     '
      CALL XERRWD (MSG, 38, 8, 0, 1, MF, 0, 0, 0.0, 0.0)
      GO TO 700
 609  MSG = 'LSODE- Pitch -  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 58, 9, 0, 2, ML, NEQ(1), 0, 0.0, 0.0)
      GO TO 700
 610  MSG = 'LSODE- Pitch -  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 58, 10, 0, 2, MU, NEQ(1), 0, 0.0, 0.0)
      GO TO 700
 611  MSG = 'LSODE- Pitch -  MAXORD (=I1) .LT. 0  '
      CALL XERRWD (MSG, 38, 11, 0, 1, MAXORD, 0, 0, 0.0, 0.0)
      GO TO 700
 612  MSG = 'LSODE- Pitch -  MXSTEP (=I1) .LT. 0  '
      CALL XERRWD (MSG, 38, 12, 0, 1, MXSTEP, 0, 0, 0.0, 0.0)
      GO TO 700
 613  MSG = 'LSODE- Pitch -  MXHNIL (=I1) .LT. 0  '
      CALL XERRWD (MSG, 38, 13, 0, 1, MXHNIL, 0, 0, 0.0, 0.0)
      GO TO 700
 614  MSG = 'LSODE- Pitch -  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 48, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0)
      GO TO 700
 615  MSG = 'LSODE- Pitch -  HMAX (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 38, 15, 0, 0, 0, 0, 1, HMAX, 0.0)
      GO TO 700
 616  MSG = 'LSODE- Pitch -  HMIN (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 38, 16, 0, 0, 0, 0, 1, HMIN, 0.0)
      GO TO 700
 617  CONTINUE
      MSG='LSODE- Pitch -  RWORK length needed, LENRW (=I1), '//
     $    'exceeds LRW (=I2)'
      CALL XERRWD (MSG, 68, 17, 0, 2, LENRW, LRW, 0, 0.0, 0.0)
      GO TO 700
 618  CONTINUE
      MSG='LSODE- Pitch -  IWORK length needed, LENIW (=I1), '//
     $    'exceeds LIW (=I2)'
      CALL XERRWD (MSG, 68, 18, 0, 2, LENIW, LIW, 0, 0.0, 0.0)
      GO TO 700
 619  MSG = 'LSODE- Pitch -  RTOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 48, 19, 0, 1, I, 0, 1, RTOLI, 0.0)
      GO TO 700
 620  MSG = 'LSODE- Pitch -  ATOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 48, 20, 0, 1, I, 0, 1, ATOLI, 0.0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'LSODE- Pitch -  EWT(I1) is R1 .LE. 0.0         '
      CALL XERRWD (MSG, 48, 21, 0, 1, I, 0, 1, EWTI, 0.0)
      GO TO 700
 622  CONTINUE
      MSG='LSODE- Pitch -  TOUT (=R1) too close to T(=R2) '//
     $    'to start integration'
      CALL XERRWD (MSG, 68, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='LSODE- Pitch -  ITASK = I1 and TOUT (=R1) behind '//
     $    'TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 68, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='LSODE- Pitch -  ITASK = 4 OR 5 and TCRIT (=R1) '//
     $    'behind TCUR (=R2)   '
      CALL XERRWD (MSG, 68, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='LSODE- Pitch -  ITASK = 4 or 5 and TCRIT (=R1) '//
     $    'behind TOUT (=R2)   '
      CALL XERRWD (MSG, 68, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'LSODE- Pitch -  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 58, 26, 0, 0, 0, 0, 0, 0.0, 0.0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'LSODE- Pitch -  Trouble in DINTDY2.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 58, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'LSODE- Pitch -  Run aborted.. apparent infinite loop     '
      CALL XERRWD (MSG, 58, 303, 2, 0, 0, 0, 0, 0.0, 0.0)
      RETURN
C----------------------- END OF SUBROUTINE LSODE ----------------------
      END SUBROUTINE LSODE1












      SUBROUTINE DSRCOM1 (RSAV, ISAV, JOB)
C***BEGIN PROLOGUE  DSRCOM1
C***SUBSIDIARY
C***PURPOSE  Save/restore ODEPACK COMMON blocks.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SSRCOM-S, DSRCOM1-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine saves or restores (depending on JOB) the contents of
C  the COMMON block DLS011, which is used internally
C  by one or more ODEPACK solvers.
C
C  RSAV = real array of length 218 or more.
C  ISAV = integer array of length 37 or more.
C  JOB  = flag indicating to save or restore the COMMON blocks:
C         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
C         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
C         A call with JOB = 2 presumes a prior call with JOB = 1.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   921116  Deleted treatment of block /EH0001/.  (ACH)
C   930801  Reduced Common block length by 2.  (ACH)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DSRCOM1
C**End
      INTEGER ISAV, JOB
      INTEGER ILS
      INTEGER I, LENILS, LENRLS
      REAL(r8) RSAV,   RLS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS011/ RLS(218), ILS(37)
!$OMP THREADPRIVATE(/DLS011/)
      DATA LENRLS/218/, LENILS/37/
C
C***FIRST EXECUTABLE STATEMENT  DSRCOM1
      IF (JOB .EQ. 2) GO TO 100
C
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      RETURN
C----------------------- END OF SUBROUTINE DSRCOM1 ----------------------
      END SUBROUTINE DSRCOM1










      SUBROUTINE DINTDY1 (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTDY1
C***SUBSIDIARY
C***PURPOSE  Interpolate solution derivatives.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SINTDY-S, DINTDY1-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DINTDY1 computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  XERRWD
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DINTDY1
C**End
      INTEGER K, NYH, IFLAG
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      REAL(r8) T, YH, DKY
      REAL(r8) ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL(r8) C, R, S, TP
      CHARACTER*80 MSG
      DIMENSION YH(NYH,*), DKY(*)
      COMMON /DLS011/ ROWNS(209),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(12), IOWNS(6),
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!$OMP THREADPRIVATE(/DLS011/)
C
C***FIRST EXECUTABLE STATEMENT  DINTDY1
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0*UROUND*(TN + HU)
      IF ((T-TP)*(T-TN) .GT. 0.0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   MSG = 'DINTDY1-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0, 0.0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY1-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY1 ----------------------
      END SUBROUTINE DINTDY1









      SUBROUTINE DPREPJ1 (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC)
C***BEGIN PROLOGUE  DPREPJ1
C***SUBSIDIARY
C***PURPOSE  Compute and process Newton iteration matrix.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SPREPJ-S, DPREPJ1-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DPREPJ1 is called by DSTODE1 to compute and process the matrix
C  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
C  Here J is computed by the user-supplied routine JAC if
C  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
C  If MITER = 3, a diagonal approximation to J is used.
C  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
C  subjected to LU decomposition in preparation for later solution
C  of linear systems with P as coefficient matrix.  This is done
C  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C  In addition to variables described in DSTODE1 and LSODE prologues,
C  communication with DPREPJ1 uses the following:
C  Y     = array containing predicted values on entry.
C  FTEM  = work array of length N (ACOR in DSTODE1).
C  SAVF  = array containing f evaluated at predicted y.
C  WM    = real work space for matrices.  On output it contains the
C          inverse diagonal matrix if MITER = 3 and the LU decomposition
C          of P if MITER is 1, 2 , 4, or 5.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C          WM(2) = H*EL0, saved for later use if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  EL0   = EL(1) (input).
C  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C          P matrix found to be singular.
C  JCUR  = output flag = 1 to indicate that the Jacobian matrix
C          (or approximation) is now current.
C  This routine also uses the COMMON variables EL0, H, TN, UROUND,
C  MITER, N, NFE, and NJE.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  DGBFA, DGEFA, DVNORM
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890504  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DPREPJ1
C**End
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      REAL(r8) Y, YH, EWT, FTEM, SAVF, WM
      REAL(r8) ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL(r8) CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DVNORM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      COMMON /DLS011/ ROWNS(209),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(12), IOWNS(6),
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!$OMP THREADPRIVATE(/DLS011/)
C
C***FIRST EXECUTABLE STATEMENT  DPREPJ1
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0) R0 = 1.0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
C Add identity matrix. -------------------------------------------------
 240  J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C If MITER = 3, construct a diagonal approximation to J and P. ---------
 300  WM(2) = HL0
      R = EL0*0.1
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
      CALL F (NEQ, TN, Y, WM(3))
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = 0.1*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = 1.0
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. 0.0) GO TO 330
        WM(I+2) = 0.1*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 1
      RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0) R0 = 1.0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
C Add identity matrix. -------------------------------------------------
 570  II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- END OF SUBROUTINE DPREPJ1 ----------------------
      END SUBROUTINE DPREPJ1







      SUBROUTINE DSOLSY1 (WM, IWM, X, TEM)
C***BEGIN PROLOGUE  DSOLSY1
C***SUBSIDIARY
C***PURPOSE  ODEPACK linear system solver.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SSOLSY-S, DSOLSY1-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls DGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls DGBSL.
C  Communication with DSOLSY1 uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DSOLSY1
C**End
      INTEGER IWM
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      REAL(r8) WM, X, TEM
      REAL(r8) ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL(r8) DI, HL0, PHL0, R
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      COMMON /DLS011/ ROWNS(209),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(12), IOWNS(6),
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!$OMP THREADPRIVATE(/DLS011/)
C
C***FIRST EXECUTABLE STATEMENT  DSOLSY1
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0 - R*(1.0 - 1.0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0) GO TO 390
 320    WM(I+2) = 1.0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
C----------------------- END OF SUBROUTINE DSOLSY1 ----------------------
      END SUBROUTINE DSOLSY1








      SUBROUTINE DSTODE1 (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
C***BEGIN PROLOGUE  DSTODE1
C***SUBSIDIARY
C***PURPOSE  Performs one step of an ODEPACK integration.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SSTODE-S, DSTODE1-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DSTODE1 performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Note:  DSTODE1 is independent of the value of the iteration method
C  indicator MITER, when this is .ne. 0, and hence is independent
C  of the type of chord method used, or the Jacobian structure.
C  Communication with DSTODE1 is done with the following variables:
C
C  NEQ    = integer array containing problem size in NEQ(1), and
C           passed as the NEQ argument in all calls to F and JAC.
C  Y      = an array of length .ge. N used as the Y argument in
C           all calls to F and JAC.
C  YH     = an NYH by LMAX array containing the dependent variables
C           and their approximate scaled derivatives, where
C           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C           j-th derivative of y(i), scaled by h**j/factorial(j)
C           (j = 0,1,...,NQ).  on entry for the first step, the first
C           two columns of YH must be set from the initial values.
C  NYH    = a constant integer .ge. N, the first dimension of YH.
C  YH1    = a one-dimensional array occupying the same space as YH.
C  EWT    = an array of length N containing multiplicative weights
C           for local error measurements.  Local errors in Y(i) are
C           compared to 1.0/EWT(i) in various error tests.
C  SAVF   = an array of working storage, of length N.
C           Also used for input of YH(*,MAXORD+2) when JSTART = -1
C           and MAXORD .lt. the current order NQ.
C  ACOR   = a work array of length N, used for the accumulated
C           corrections.  On a successful return, ACOR(i) contains
C           the estimated one-step local error in Y(i).
C  WM,IWM = real and integer work arrays associated with matrix
C           operations in chord iteration (MITER .ne. 0).
C  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C           and P = I - h*el0*JAC, if a chord method is being used.
C  SLVS   = name of routine to solve linear system in chord iteration.
C  CCMAX  = maximum relative change in h*el0 before PJAC is called.
C  H      = the step size to be attempted on the next step.
C           H is altered by the error control algorithm during the
C           problem.  H can be either positive or negative, but its
C           sign must remain constant throughout the problem.
C  HMIN   = the minimum absolute value of the step size h to be used.
C  HMXI   = inverse of the maximum absolute value of h to be used.
C           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
C           HMIN and HMXI may be changed at any time, but will not
C           take effect until the next change of h is considered.
C  TN     = the independent variable. TN is updated on each step taken.
C  JSTART = an integer used for input only, with the following
C           values and meanings:
C                0  perform the first step.
C            .gt.0  take a new step continuing from the last.
C               -1  take the next step with a new value of H, MAXORD,
C                     N, METH, MITER, and/or matrix parameters.
C               -2  take the next step with a new value of H,
C                     but with other inputs unchanged.
C           On return, JSTART is set to 1 to facilitate continuation.
C  KFLAG  = a completion code with the following meanings:
C                0  the step was succesful.
C               -1  the requested error could not be achieved.
C               -2  corrector convergence could not be achieved.
C               -3  fatal error in PJAC or SLVS.
C           A return with KFLAG = -1 or -2 means either
C           abs(H) = HMIN or 10 consecutive failures occurred.
C           On a return with KFLAG negative, the values of TN and
C           the YH array are as of the beginning of the last
C           step, and H is the last step size attempted.
C  MAXORD = the maximum order of integration method to be allowed.
C  MAXCOR = the maximum number of corrector iterations allowed.
C  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C  MXNCF  = maximum number of convergence failures allowed.
C  METH/MITER = the method flags.  See description in driver.
C  N      = the number of first-order differential equations.
C  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
C  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  DCFODE, DVNORM
C***COMMON BLOCKS    DLS011
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DSTODE1
C**End
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      REAL(r8) Y, YH, YH1, EWT, SAVF, ACOR, WM
      REAL(r8) CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL(r8) DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      COMMON /DLS011/ CONIT, CRATE, EL(13), ELCO(13,12),
     1   HOLD, RMAX, TESCO(3,12),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, IOWND(12),
     3   IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
!$OMP THREADPRIVATE(/DLS011/)
C
C***FIRST EXECUTABLE STATEMENT  DSTODE1
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set to 2
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0
      RC = 0.0
      EL0 = 1.0
      CRATE = 0.7
      HOLD = H
      MEO = METH
      NSLP = 0
      IPUP = MITER
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If the caller has changed MAXORD to a value less than the current
C order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5/(NQ+2)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0/L
      RHDN = 1.0/(1.3*DDN**EXDN + 0.0000013)
      RH = MIN(RHDN,1.0d0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C DCFODE is called to get all the integration coefficients for the
C current METH.  Then the EL vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0d0,ABS(H)*HMXI*RH)
      R = 1.0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal Triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the R.M.S. norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - h*el(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0
      NSLP = NST
      CRATE = 0.7
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M.gt.0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(0.2*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0d0,1.5d0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .LE. 1.0) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0
      IF (ABS(H) .LE. HMIN*1.00001) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0/(L+1)
      RHUP = 1.0/(1.4*DUP**EXUP + 0.0000014)
 540  EXSM = 1.0/L
      RHSM = 1.0/(1.2*DSM**EXSM + 0.0000012)
      RHDN = 0.0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0/NQ
      RHDN = 1.0/(1.3*DDN**EXDN + 0.0000013)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0) RH = 1.0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2d0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, l, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0
 700  R = 1.0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- END OF SUBROUTINE DSTODE1 ----------------------
      END SUBROUTINE DSTODE1





      END MODULE LSODE1_MOD
