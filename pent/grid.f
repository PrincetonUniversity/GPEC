c-----------------------------------------------------------------------
c     NON AMBIPOLAR TRANSPORT and TOROIDAL TORQUE
c     form irregular grids for splines with resonances
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. grid_mod
c     1. grid_power
c     2. grid_cos
c     3. grid_logistic
c-----------------------------------------------------------------------
c     subprogram 0. grid_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE grid_mod
      USE pentio_mod

      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     function a. linspace.
c     Create linear grid.
c     args - xmin : grid minimum
c          - xmax : grid maximum
c          - num  : ubound of points in grid (NOTE size = num+1)
c     returns     : linear grid [xmin,xmax] indexed from zero
c-----------------------------------------------------------------------
      FUNCTION linspace(xmin,xmax,num)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: xmin, xmax
      INTEGER :: num,i
      REAL(r8), DIMENSION(1:num) :: linspace

      linspace = (/(i/(num-1.0),i=0,num-1)/)*(xmax-xmin)+xmin
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION linspace

c-----------------------------------------------------------------------
c     function b. powspace.
c     Create power function grid consentrated at x0.
c     args - xmin : grid minimum
c          - xmax : grid maximum
c          - pow  : power of grid consentration
c          - num  : power with which endpoint derivative approaches zero
c     returns     : 2D array containing x and dx/dnorm
c-----------------------------------------------------------------------
      FUNCTION powspace(xmin,xmax,pow,num,endpoints)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      CHARACTER(8), INTENT(IN) :: endpoints
      REAL(r8), INTENT(IN) :: xmin, xmax
      INTEGER, INTENT(IN) :: num,pow
      REAL(r8), DIMENSION(1:num) :: x
      REAL(r8), DIMENSION(0:1,1:num) :: powspace

      REAL(r8) :: deltay,deltax
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     Check inputs.
c-----------------------------------------------------------------------
      IF(xmax<=xmin)THEN
         PRINT *,"xmin,xmax = ",xmin,xmax
         CALL pentio_stop(
     $     "Error in limits for creating power law space")
      ENDIF
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      SELECT CASE (endpoints)
         CASE ("lower")
            x = linspace(-1.0_r8,0.0_r8,num)
         CASE ("upper")
            x = linspace(0.0_r8,1.0_r8,num)
         CASE ("both")
            x = linspace(-1.0_r8,1.0_r8,num)
         CASE DEFAULT
            x = linspace(-1.0_r8,1.0_r8,num)
      END SELECT
      powspace(1,:) = ( (x-1)*(x+1) )**pow
      SELECT CASE (pow)
         CASE (1)
            powspace(0,:) = -x+x**3/3
         CASE (2)
            powspace(0,:) = x-(2*x**3)/3+x**5/5
         CASE (3)
            powspace(0,:) =-x + x**3 - (3*x**5)/5 + x**7/7
         CASE (4)
            powspace(0,:) =x -(4*x**3)/3+(6*x**5)/5-(4*x**7)/7+x**9/9
         CASE (5)
            powspace(0,:) =-x + (5*x**3)/3 - 2*x**5 + (10*x**7)/7 - 
     $           (5*x**9)/9 + x**11/11
         CASE (6)
            powspace(0,:) =x - 2*x**3 + 3*x**5 - (20*x**7)/7 +(5*x**9)/3 
     $           - (6*x**11)/11 + x**13/13
         CASE (7)
            powspace(0,:) =-x + (7*x**3)/3 - (21*x**5)/5 + 5*x**7 - 
     $           (35*x**9)/9 + (21*x**11)/11 - (7*x**13)/13 + x**15/15
         CASE (8)
            powspace(0,:) =x-(8*x**3)/3+(28*x**5)/5-8*x**7+(70*x**9)/9 
     $           -(56*x**11)/11 + (28*x**13)/13 - (8*x**15)/15+x**17/17
         CASE (9)
            powspace(0,:) = -x + 3*x**3 - (36*x**5)/5 + 12*x**7 - 
     $           14*x**9 + (126*x**11)/11 - (84*x**13)/13 + 
     $           (12*x**15)/5 - (9*x**17)/17 + x**19/19
         CASE DEFAULT
            PRINT *, "WARNING: Power not in analytic database,"//
     $           "attempting numeric integration"
            CALL spline_alloc(spl,num-1,1)
            spl%xs = x
            spl%fs(:,1) = powspace(0,:)
            CALL spline_fit(spl,"extrap")
            CALL spline_int(spl)
            powspace(0,:) = spl%fsi(:,1)
      END SELECT

      !stretch to the desired range
      deltay = powspace(0,num) - powspace(0,1)
      deltax = xmax-xmin
      powspace(:,:) = powspace(:,:)*deltax/deltay
c      print *,pow,endpoints,deltax,deltay
      !move by offset
      powspace(0,:) = powspace(0,:) - powspace(0,1) + xmin
      !dirivative with respect to [0,1]
      powspace(1,:) = powspace(1,:)*(x(num)-x(1))

c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION powspace

c-----------------------------------------------------------------------
c     function c. pow2space.
c     Create power function grid consentrated at single resonance x0.
c     args - norm : Normalized array [0,1]
c          - xmin : grid minimum
c          - x0   : grid consentrated point
c          - xmax : grid maximum
c          - pow  : power of grid consentration
c     returns     : 2D array containing x and dx/dnorm
c-----------------------------------------------------------------------
      FUNCTION pow2space(xmin,x0,xmax,pow,num)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: xmin, xmax, x0,pow
      INTEGER, INTENT(IN) :: num
      REAL(r8), DIMENSION(1:num) :: norm
      REAL(r8), DIMENSION(0:1,1:num) :: pow2space

      REAL(r8) :: c0,c1,c2
c-----------------------------------------------------------------------
c     Check inputs.
c-----------------------------------------------------------------------
      IF(xmax<=xmin .OR. x0<xmin .OR. x0>xmax) CALL pentio_stop(
     $     "Error in limits for creating power law space")
      IF(MOD(pow,2)/=1) CALL pentio_stop(
     $     "Power law grid consentration must use odd power.")
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      norm = linspace(0.0_r8,1.0_r8,num)
      c0 = x0
      IF(x0>xmin)THEN
         c1 = 1/(1+((xmax-x0)/(x0-xmin))**(1/pow))
         c2 = c1**pow/(x0-xmin)
      ELSE
         c1 = 0.0
         c2 = 1.0/(xmax-c0)
      ENDIF
      
      pow2space(0,:) = (norm-c1)**pow/c2+c0
      pow2space(1,:) = pow*(norm-c1)**(pow-1)/c2
      !print *, "pow2space xn=",norm(:5)
      !print *, "powc = ", num,pow, c0, c1, c2
      !print *, "powmin,max=",MINVAL(pow2space(0,:)),MAXVAL(pow2space(0,:))
      !print *, "pow2space = ",pow2space(:,num-5:num)
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION pow2space

c-----------------------------------------------------------------------
c     function c. cosspace.
c     Create cosine function grid consentrated at enpoints.
c     args - norm : Normalized array [0,1]
c          - xmin : grid minimum
c          - xmax : grid maximum
c     returns     : 2D array containing x and dx/dnorm
c-----------------------------------------------------------------------
      FUNCTION cosspace(xmin,xmax,num)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: xmin, xmax
      INTEGER, INTENT(IN) :: num
      REAL(r8), DIMENSION(1:num) :: norm
      REAL(r8), DIMENSION(0:1,1:num) :: cosspace

      REAL(r8) :: c0,c1

c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      norm = linspace(0.0_r8,1.0_r8,num)
      c0 = xmin
      c1 = (xmax-xmin)/2
      IF(c1<=0) 
     $     CALL pentio_stop('Error ordering points for cosine grid.')
      cosspace(0,:) = c1*(1.0-COS(norm*pi))+c0
      cosspace(1,:) = c1*pi*SIN(norm*pi)

c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION cosspace

c-----------------------------------------------------------------------
c     function d. logspace.
c     Create cosine function grid consentrated at enpoints.
c     args - norm : Normalized array [0,1]
c          - xmin : grid minimum
c          - xmax : grid maximum
c          - pow  : exponential power with which xmax is approached
c     returns     : 2D array containing x and dx/dnorm
c-----------------------------------------------------------------------
      FUNCTION logspace(xmin,xmax,pow,num)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: xmin, xmax, pow
      INTEGER, INTENT(IN) :: num
      REAL(r8), DIMENSION(1:num) :: norm
      REAL(r8), DIMENSION(0:1,1:num) :: logspace

      REAL(r8) :: c0,c1

c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      norm = linspace(0.0_r8,1.0_r8,num)
      c0 = xmin
      c1 = xmax-xmin
      IF(c1<=0) 
     $     CALL pentio_stop('Error ordering points in logistic grid')
      logspace(0,:) = c1/(1-EXP(-pow*2*(norm-0.5)))+c0
      logspace(1,:) = c1*pow*EXP(-pow*2*(norm-0.5))/
     $        (1-EXP(-pow*2*norm))**2.0
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION logspace


c-----------------------------------------------------------------------
c     subprogram 1. powergrid.
c     Create grid with power law concentration at multiple points.
c     args - x    : A pre-formed linear grid.
c          - ydy  : 2D array of the new grid, and its derivative
c          - yset : Array of grid locations such that dy/dx(y)=0
c                   should include ymin,ymax on either end
c          - pow  : Power with which the grid spacing approaches 0 
c                   at each y in yset
c          - endpoints : Specify which points should have dy=0
c-----------------------------------------------------------------------
      SUBROUTINE powergrid(x,ydy,yset,pow,endpoints)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      CHARACTER(8), INTENT(IN) :: endpoints
      INTEGER, INTENT(IN) :: pow
      REAL(r8), DIMENSION(:), INTENT(IN) :: yset
      REAL(r8), DIMENSION(:), INTENT(IN) :: x
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: ydy

      INTEGER :: i,block,nblocks,nset,l,r

c-----------------------------------------------------------------------
c     Doubel check inputs.
c-----------------------------------------------------------------------
      nset = SIZE(yset,DIM=1)
      IF(nset<2)
     $     CALL pentio_stop('Must specify at least 2 pts for powergrid')
      IF(nset==2 .AND. endpoints=="none")
     $     CALL pentio_stop('Must specify at least 1 concentration '//
     $     'point for powergrid')
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      nblocks = SIZE(yset,DIM=1)-1
      block = SIZE(x,DIM=1)/nblocks
      ! include remainder in first block
      r = MOD(SIZE(x,DIM=1),block)
      IF(endpoints=="none")THEN
         ydy(:,:block+r) = powspace(yset(1),yset(2),pow,block+r,"upper")
      ELSEIF(endpoints=="upper")THEN
         ydy(:,:block+r) = powspace(yset(1),yset(2),pow,block+r,"upper")
      ELSEIF(nblocks==1 .AND. endpoints=="lower")THEN
         ydy(:,:block+r) = powspace(yset(1),yset(2),pow,block+r,"lower")
      ELSE
         ydy(:,:block+r) = powspace(yset(1),yset(2),pow,block+r,"both")
      ENDIF
      ydy(2,:block+r) = ydy(2,:block+r)/(x(block+r)-x(1))
      IF(nblocks>1)THEN
         DO i = 2,nblocks
c            print *,SIZE(x,DIM=1),r+(i-1)*block,r+i*block,block
            IF(i==nblocks .AND. (endpoints=="lower" .OR. 
     $           endpoints=="none"))THEN
               ydy(:,r+(i-1)*block:r+i*block)=powspace(
     $              yset(i),yset(i+1),pow,block+1,"lower")
            ELSE
               ydy(:,r+(i-1)*block:r+i*block)=powspace(
     $           yset(i),yset(i+1),pow,block+1,"both")
            ENDIF
            ydy(2,r+(i-1)*block:r+i*block) = 
     $           ydy(2,r+(i-1)*block:r+i*block)/
     $           (x(i*block+r)-x(r+(i-1)*block))
         ENDDO
      ENDIF
    
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE powergrid

c-----------------------------------------------------------------------
c     subprogram 1. grid_power.
c     Create grid with power law concentration multiple points.
c     args - x    : A pre-formed linear grid [0,1]. Indexed from 0.
c          - ydy  : 2D array of the new grid, and its derivative
c          - yset : Array of (ymin,(yres),ymax) such that y(y'=0) = yres
c          - pow  : Power with which the grid approaches the resonance
c          - endpoints : Specify which points should have dy=0
c-----------------------------------------------------------------------
      SUBROUTINE grid_power(x,ydy,yset,pow,endpoints)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      CHARACTER(8), INTENT(IN) :: endpoints
      REAL(r8), INTENT(IN) :: pow
      REAL(r8), DIMENSION(:), INTENT(IN) :: yset
      REAL(r8), DIMENSION(:), INTENT(IN) :: x
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: ydy

      INTEGER :: i,block,nblocks,nyset,nset,nend,l,r
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norm,set

c-----------------------------------------------------------------------
c     Set up endpoints.
c-----------------------------------------------------------------------
      nyset = SIZE(yset,DIM=1)
      nend=0
      l=0
      SELECT  CASE (endpoints)
         CASE ("none")
            nend = 0
         CASE ("lower")
            IF(yset(1)/=yset(2))THEN
               nend = 1
               l = 1
            ENDIF
         CASE ("upper")
            IF(yset(nyset-1)/=yset(nyset)) nend = 1
         CASE ("both")
            IF(yset(1)/=yset(2))THEN
               nend = nend+1
               l = 1
            ENDIF
            IF(yset(nyset-1)/=yset(nyset)) nend = nend+1
         CASE DEFAULT
            nend = 0
      END SELECT
      nset = nyset+nend+(nyset+nend-3)
      !IF(yset(1)==yset(2)) nset = nset-1
      !IF(yset(nyset-1)==yset(nyset)) nset=nset-1
      ALLOCATE(set(nset))
      set(1:2) = yset(1)
      set(nset-1:nset) = yset(nyset)
      DO i=2,nyset-1
         set((i-1+l)*2) = yset(i)
      ENDDO
      ! fill in the half points
      !print *,"sparse set =",set
      nset = SIZE(set,DIM=1)
      IF(nset>3)THEN
         DO i = 3,nset-2,2
            set(i) = (set(i-1)+set(i+1))/2
         ENDDO
      ENDIF
      !print *, "yset =",yset
      !print *, "nend =",nend,"lower=",l,"set = ",set
c-----------------------------------------------------------------------
c     Doubel check inputs.
c-----------------------------------------------------------------------
      IF(nset<3)
     $     CALL pentio_stop('Error setting limits for power grid.')
      DO i =2,nset
         IF(set(i)<set(i-1)) CALL pentio_stop(
     $     'Error ordering points for power grid.')
      ENDDO
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      nblocks = nset/2
      block = SIZE(x,DIM=1)/nblocks
      r = MOD(SIZE(x,DIM=1),block)
      ALLOCATE(norm(block+r))
      ydy(:,:block+r) = pow2space(set(1),set(2),set(3),pow,block+r)
      ydy(2,:block+r) = ydy(2,:block+r)/x(block+r)
      !print *, block+r
      DEALLOCATE(norm)
      IF(nblocks>1)THEN
         ALLOCATE(norm(block))
         IF(nblocks>2)THEN
            DO i=2,nblocks-1
               !print *, r+(i-1)*block
               ydy(:,r+(i-1)*block:r+i*block) = pow2space(set(2*i-1),
     $              set(2*i),set(2*i+1),pow,block+1)
               ydy(2,r+(i-1)*block:r+i*block)=
     $              ydy(2,r+(i-1)*block:r+i*block)/x(block+1)
            ENDDO
         ENDIF
         ! "upper" could have made nset even
         ydy(:,r+(nblocks-1)*block:r+nblocks*block) = 
     $        pow2space(set(nset-2),set(nset-1),set(nset),pow,block+1)
            ydy(2,r+(nblocks-1)*block:r+nblocks*block)=
     $           ydy(2,r+(nblocks-1)*block:r+nblocks*block)/x(block+1)
         !print *,r+(nblocks-1)*block,r+nblocks*block
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE grid_power


      END MODULE grid_mod
