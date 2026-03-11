!@GENERAL PERTURBED EQUILIBRIUM CODE

    module grid_mod
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Functional grid creation.
    !
    !*DEPENDENCIES:
    !   
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !   All functions are public
    !
    !*PUBLIC DATA MEMBERS:
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nikolas.logan@columbia.edu
    !-----------------------------------------------------------------------
    
    use spline_mod
    
    implicit none
    
    contains
    
    !=======================================================================
    function linspace(xmin,xmax,num)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Create linear grid.
    !
    !*ARGUMENTS:
    !    xmin : real (in)
    !       Grid minimum.
    !    xmax : real (in)
    !       Grid maximum.
    !    num : integer (in)
    !       Number of points.
    !
    !*RETURNS:
    !     real array.
    !        Linear grid from xmin to xmax.
    !-----------------------------------------------------------------------
        implicit none
        !declare function
        real(r8), dimension(1:num) :: linspace
        ! declare arguments
        real(r8), intent(in) :: xmin, xmax
        integer :: num
        ! declare variables
        integer :: i

        linspace = (/(i/(num-1.0),i=0,num-1)/)*(xmax-xmin)+xmin
        return
    end function linspace

    subroutine linspace_sub(xmin,xmax,num,linspacevar)
      !----------------------------------------------------------------------- 
      !*DESCRIPTION: 
      !   Create linear grid.
      !
      !*ARGUMENTS:
      !    xmin : real (in)
      !       Grid minimum.
      !    xmax : real (in)
      !       Grid maximum.
      !    num : integer (in)
      !       Number of points.
      !
      !*RETURNS:
      !     real array.
      !        Linear grid from xmin to xmax.
      !-----------------------------------------------------------------------
          implicit none
          !declare function
          real(r8), dimension(1:num), intent(out) :: linspacevar
          ! declare arguments
          real(r8), intent(in) :: xmin, xmax
          integer :: num
          ! declare variables
          integer :: i
  
          linspacevar = (/(i/(num-1.0),i=0,num-1)/)*(xmax-xmin)+xmin
          return
      end subroutine linspace_sub

    !=======================================================================
    function powspace(xmin,xmax,pow,num,endpoints)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Create grid spaced such that the spacing approaches zero at the rate
    !   specified by power.
    !   -> Derivatives are with respect to a unit linear space.
    !
    !*ARGUMENTS:
    !    xmin : real (in)
    !       Grid minimum.
    !    xmax : real (in)
    !       Grid maximum.
    !   pow : real.
    !       power of grid consentration
    !    num : integer (in)
    !       Number of points.
    !   endpoints : character
    !       Concentration at edges of grid (left,right,both) 
    !
    !*RETURNS:
    !     powspace : real 2D array.
    !        x,dx/dnorm where x is the grid and norm is a linear space.
    !-----------------------------------------------------------------------
        ! declare function and arguments
        real(r8), dimension(2,num) :: powspace
        character(*), intent(in) :: endpoints
        real(r8), intent(in) :: xmin, xmax
        integer, intent(in) :: num,pow
        ! declare variables
        real(r8), dimension(1:num) :: x
        real(r8) :: deltay,deltax
        type(spline_type) :: spl
    
        ! check input validity
        if(xmax<=xmin)then
           print *,"xmin,xmax = ",xmin,xmax
           stop 'ERROR: powspace - xmax must be less than xmin.'
        endif
    
        ! endpoint characterization
        select case (endpoints)
           case ("lower")
              x = linspace(-1.0_r8,0.0_r8,num)
           case ("upper")
              x = linspace(0.0_r8,1.0_r8,num)
           case ("both")
              x = linspace(-1.0_r8,1.0_r8,num)
           case default
              stop "ERROR: powspace - not a valid endpoint"
        end select
        
        ! concentrate grid points
        powspace(2,:) = abs( (x-1)*(x+1) )**pow
        if(any(powspace(2,2:num-1)<=0))then
           print *,x
           print *,''
           print *,( (x-1)*(x+1) )**pow
        endif
        select case (pow)
           case (1)
                powspace(1,:) = -x+x**3/3
           case (2)
                powspace(1,:) = x-(2*x**3)/3+x**5/5
           case (3)
                powspace(1,:) =-x + x**3 - (3*x**5)/5 + x**7/7
           case (4)
                powspace(1,:) =x -(4*x**3)/3+(6*x**5)/5-(4*x**7)/7+x**9/9
           case (5)
                powspace(1,:) =-x + (5*x**3)/3 - 2*x**5 + (10*x**7)/7 &
                    - (5*x**9)/9 + x**11/11
           case (6)
                powspace(1,:) =x - 2*x**3 + 3*x**5 - (20*x**7)/7 +(5*x**9)/3 &
                    - (6*x**11)/11 + x**13/13
           case (7)
                powspace(1,:) =-x + (7*x**3)/3 - (21*x**5)/5 + 5*x**7 & 
                    - (35*x**9)/9 + (21*x**11)/11 - (7*x**13)/13 + x**15/15
           case (8)
                powspace(1,:) =x-(8*x**3)/3+(28*x**5)/5-8*x**7+(70*x**9)/9 &
                    - (56*x**11)/11 + (28*x**13)/13 - (8*x**15)/15+x**17/17
           case (9)
                powspace(1,:) = -x + 3*x**3 - (36*x**5)/5 + 12*x**7 & 
                    - 14*x**9 + (126*x**11)/11 - (84*x**13)/13 & 
                    + (12*x**15)/5 - (9*x**17)/17 + x**19/19
           case default
              print *, "!! WARNING: power not in analytic database, &
                    &attempting numeric integration"
              call spline_alloc(spl,num-1,1)
              spl%xs = x
              spl%fs(:,1) = powspace(1,:)
              call spline_fit(spl,"extrap")
              call spline_int(spl)
              powspace(1,:) = spl%fsi(:,1)
        end select
    
        !stretch to the desired range
        deltay = powspace(1,num) - powspace(1,1)
        deltax = xmax-xmin
        powspace(:,:) = powspace(:,:)*deltax/deltay
        !move by offset
        powspace(1,:) = powspace(1,:) - powspace(1,1) + xmin
        !dirivative with respect to [0,1]???
        powspace(2,:) = powspace(2,:)*(x(num)-x(1))
        return
    end function powspace

    !=======================================================================
    subroutine powspace_sub(xmin,xmax,pow,num,endpoints,powspacevar)
      !----------------------------------------------------------------------- 
      !*DESCRIPTION: 
      !   Create grid spaced such that the spacing approaches zero at the rate
      !   specified by power.
      !   -> Derivatives are with respect to a unit linear space.
      !
      !*ARGUMENTS:
      !    xmin : real (in)
      !       Grid minimum.
      !    xmax : real (in)
      !       Grid maximum.
      !   pow : real.
      !       power of grid consentration
      !    num : integer (in)
      !       Number of points.
      !   endpoints : character
      !       Concentration at edges of grid (left,right,both) 
      !
      !*RETURNS:
      !     powspace : real 2D array.
      !        x,dx/dnorm where x is the grid and norm is a linear space.
      !-----------------------------------------------------------------------
          ! declare function and arguments
          real(r8), dimension(2,num), intent(out) :: powspacevar
          character(*), intent(in) :: endpoints
          real(r8), intent(in) :: xmin, xmax
          integer, intent(in) :: num,pow
          ! declare variables
          real(r8), dimension(1:num) :: x
          real(r8) :: deltay,deltax
          type(spline_type) :: spl
      
          ! check input validity
          if(xmax<=xmin)then
             print *,"xmin,xmax = ",xmin,xmax
             stop 'ERROR: powspace - xmax must be less than xmin.'
          endif
      
          ! endpoint characterization
          select case (endpoints)
             case ("lower")
                call linspace_sub(-1.0_r8,0.0_r8,num,x)
             case ("upper")
                call linspace_sub(0.0_r8,1.0_r8,num,x)
             case ("both")
                call linspace_sub(-1.0_r8,1.0_r8,num,x)
             case default
                stop "ERROR: powspace - not a valid endpoint"
          end select
          
          ! concentrate grid points
          powspacevar(2,:) = abs( (x-1)*(x+1) )**pow
          if(any(powspacevar(2,2:num-1)<=0))then
             print *,x
             print *,''
             print *,( (x-1)*(x+1) )**pow
          endif
          select case (pow)
             case (1)
                  powspacevar(1,:) = -x+x**3/3
             case (2)
                  powspacevar(1,:) = x-(2*x**3)/3+x**5/5
             case (3)
                  powspacevar(1,:) =-x + x**3 - (3*x**5)/5 + x**7/7
             case (4)
                  powspacevar(1,:) =x -(4*x**3)/3+(6*x**5)/5-(4*x**7)/7+x**9/9
             case (5)
                  powspacevar(1,:) =-x + (5*x**3)/3 - 2*x**5 + (10*x**7)/7 &
                      - (5*x**9)/9 + x**11/11
             case (6)
                  powspacevar(1,:) =x - 2*x**3 + 3*x**5 - (20*x**7)/7 +(5*x**9)/3 &
                      - (6*x**11)/11 + x**13/13
             case (7)
                  powspacevar(1,:) =-x + (7*x**3)/3 - (21*x**5)/5 + 5*x**7 & 
                      - (35*x**9)/9 + (21*x**11)/11 - (7*x**13)/13 + x**15/15
             case (8)
                  powspacevar(1,:) =x-(8*x**3)/3+(28*x**5)/5-8*x**7+(70*x**9)/9 &
                      - (56*x**11)/11 + (28*x**13)/13 - (8*x**15)/15+x**17/17
             case (9)
                  powspacevar(1,:) = -x + 3*x**3 - (36*x**5)/5 + 12*x**7 & 
                      - 14*x**9 + (126*x**11)/11 - (84*x**13)/13 & 
                      + (12*x**15)/5 - (9*x**17)/17 + x**19/19
             case default
                print *, "WARNING: power not in analytic database, &
                      &attempting numeric integration"
                call spline_alloc(spl,num-1,1)
                spl%xs = x
                spl%fs(:,1) = powspacevar(1,:)
                call spline_fit(spl,"extrap")
                call spline_int(spl)
                powspacevar(1,:) = spl%fsi(:,1)
          end select
      
          !stretch to the desired range
          deltay = powspacevar(1,num) - powspacevar(1,1)
          deltax = xmax-xmin
          powspacevar(:,:) = powspacevar(:,:)*deltax/deltay
          !move by offset
          powspacevar(1,:) = powspacevar(1,:) - powspacevar(1,1) + xmin
          !dirivative with respect to [0,1]???
          powspacevar(2,:) = powspacevar(2,:)*(x(num)-x(1))
          return
      end subroutine powspace_sub

    !=======================================================================
    function cosspace(xmin,xmax,num)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Create power grid.
    !
    !*ARGUMENTS:
    !    xmin : real (in)
    !       Grid minimum.
    !    xmax : real (in)
    !       Grid maximum.
    !    num : integer (in)
    !       Number of points.
    !
    !*RETURNS:
    !     cosspace : real 2D array.
    !        x,dx/dnorm where x is the grid and norm is a linear space.
    !-----------------------------------------------------------------------
        real(r8), intent(in) :: xmin, xmax
        integer, intent(in) :: num
        real(r8), dimension(1:num) :: norm
        real(r8), dimension(2,num) :: cosspace
        ! variables
        real(r8) :: c0,c1
    
        ! calculations
        norm = linspace(0.0_r8,1.0_r8,num)
        c0 = xmin
        c1 = (xmax-xmin)/2
        if(c1<=0) stop 'ERROR: cosspace - xmax must be less than xmin.'
        cosspace(1,:) = c1*(1.0-cos(norm*pi))+c0
        cosspace(2,:) = c1*pi*sin(norm*pi)
        return
    end function cosspace

    !=======================================================================
    function logspace(xmin,xmax,pow,num)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Create power grid.
    !
    !*ARGUMENTS:
    !    xmin : real (in)
    !       Grid minimum.
    !    xmax : real (in)
    !       Grid maximum.
    !   pow : real
    !       Exponential power with which xmax is approached.
    !    num : integer (in)
    !       Number of points.
    !
    !*RETURNS:
    !     cosspace : real 2D array.
    !        x,dx/dnorm where x is the grid and norm is a linear space.
    !-----------------------------------------------------------------------
        ! declare function and arguments
        real(r8), intent(in) :: xmin, xmax, pow
        integer, intent(in) :: num
        real(r8), dimension(2,num) :: logspace
        ! declare variables
        real(r8), dimension(num) :: norm
        real(r8) :: c0,c1
    
        ! calculations
        norm = linspace(0.0_r8,1.0_r8,num)
        c0 = xmin
        c1 = xmax-xmin
        if(c1<=0) stop 'ERROR: logspace - xmax must be less than xmin.'
        logspace(1,:) = c1/(1-exp(-pow*2*(norm-0.5)))+c0
        logspace(2,:) = c1*pow*exp(-pow*2*(norm-0.5))/(1-exp(-pow*2*norm))**2.0
       return
    end function logspace
    
    
    end module grid_mod
