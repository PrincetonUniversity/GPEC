c-----------------------------------------------------------------------
c     file utils.f
c     a few useful utilities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. utils_mod.
c     1. asinh.
c     2. interpolate.
c     3. cinterpolate.
c     4. put_seed.
c     5. bubble.
c     6. my_transpose.
c     7. my_ctranspose.
c     8. mypack.
c     9. ident.
c    10. powspace
c-----------------------------------------------------------------------
c     subprogram 0. utils.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      MODULE utils_mod
      USE local_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. asinh.
c     computes inverse hyperbolic sin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION asinh(x) RESULT(y)

      REAL(r8), INTENT(IN) :: x
      REAL(r8) :: y

      REAL(r8) :: x2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      x2=x*x
      IF(x < -1e3)THEN
         y=-LOG(-2.*x)-.25/x2
      ELSE IF(ABS(x) < 1e-2)THEN
         y=x*(1-x2/6*(1-x2*9/20))
      ELSE IF(x > 1e3)THEN
         y=LOG(2.*x)+.25/x2
      ELSE
         y=LOG(x+SQRT(1.+x2))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION asinh
c-----------------------------------------------------------------------
c     subprogram 2. interpolate.
c     polynomial interpolation of function and first two derivative 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION interpolate(xs,fs,x) RESULT(f)
      
      REAL(r8), DIMENSION(:), INTENT(IN) :: xs
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: fs
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(0:2,SIZE(fs,2)) :: f

      INTEGER :: i,j,k,n
      REAL(r8) :: term0,term1,term2,dx,dxmin,xx
      REAL(r8), PARAMETER :: eps=1e-10
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      n=SIZE(xs)
      f=0
      xx=x
      dx=MINVAL(ABS(xx-xs))
      dxmin=eps*MINVAL(ABS(xs(2:n)-xs(1:n-1)))
      IF(dx < dxmin)xx=xx+dxmin
c-----------------------------------------------------------------------
c     start loop over nodes.
c-----------------------------------------------------------------------
      DO i=1,n
         term0=1
         term1=0
         term2=0
c-----------------------------------------------------------------------
c     compute term0.
c-----------------------------------------------------------------------
         DO j=1,n
            IF(j /= i)term0=term0*(xx-xs(j))/(xs(i)-xs(j))
         ENDDO
c-----------------------------------------------------------------------
c     compute term1 and term2.
c-----------------------------------------------------------------------
         DO j=1,n
            IF(j /= i)THEN
               term1=term1+term0/(xx-xs(j))
               DO k=1,n
                  IF(k /= i .AND. k /= j)THEN
                     term2=term2
     $                    +term0/((xx-xs(j))*(xx-xs(k)))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     finish loop over nodes.
c-----------------------------------------------------------------------
         f(0,:)=f(0,:)+fs(i,:)*term0
         f(1,:)=f(1,:)+fs(i,:)*term1
         f(2,:)=f(2,:)+fs(i,:)*term2
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION interpolate
c-----------------------------------------------------------------------
c     subprogram 3. cinterpolate.
c     polynomial interpolation of complex function and first two derivative 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION cinterpolate(xs,fs,x) RESULT(f)
      
      REAL(r8), DIMENSION(:), INTENT(IN) :: xs
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: fs
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(0:2,SIZE(fs,2)) :: f

      INTEGER :: i,j,k,n
      REAL(r8) :: dx,dxmin,xx
      COMPLEX(r8) :: term0,term1,term2
      REAL(r8), PARAMETER :: eps=1e-10
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      n=SIZE(xs)
      f=0
      xx=x
      dx=MINVAL(ABS(xx-xs))
      dxmin=eps*MINVAL(ABS(xs(2:n)-xs(1:n-1)))
      IF(dx < dxmin)xx=xx+dxmin
c-----------------------------------------------------------------------
c     start loop over nodes.
c-----------------------------------------------------------------------
      DO i=1,n
         term0=1
         term1=0
         term2=0
c-----------------------------------------------------------------------
c     compute term0.
c-----------------------------------------------------------------------
         DO j=1,n
            IF(j /= i)term0=term0*(xx-xs(j))/(xs(i)-xs(j))
         ENDDO
c-----------------------------------------------------------------------
c     compute term1 and term2.
c-----------------------------------------------------------------------
         DO j=1,n
            IF(j /= i)THEN
               term1=term1+term0/(xx-xs(j))
               DO k=1,n
                  IF(k /= i .AND. k /= j)THEN
                     term2=term2
     $                    +term0/((xx-xs(j))*(xx-xs(k)))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     finish loop over nodes.
c-----------------------------------------------------------------------
         f(0,:)=f(0,:)+fs(i,:)*term0
         f(1,:)=f(1,:)+fs(i,:)*term1
         f(2,:)=f(2,:)+fs(i,:)*term2
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION cinterpolate
c-----------------------------------------------------------------------
c     subprogram 4. put_seed.
c     deposits random number seed.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE put_seed(seed)

      INTEGER, INTENT(IN) :: seed

      INTEGER :: seed_size
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed_array
c-----------------------------------------------------------------------
c     deposit seed.
c-----------------------------------------------------------------------
      CALL RANDOM_SEED(SIZE=seed_size)
      ALLOCATE(seed_array(seed_size))
      seed_array=0
      seed_array(1)=seed
      CALL RANDOM_SEED(PUT=seed_array)
      DEALLOCATE(seed_array)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE put_seed
c-----------------------------------------------------------------------
c     subprogram 5. bubble.
c     performs a bubble sort in decreasing order of value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bubble(key,index,mmin,mmax)

      REAL(r8), DIMENSION(:), INTENT(IN) :: key
      INTEGER, DIMENSION(:), INTENT(INOUT) :: index
      INTEGER :: mmin,mmax

      LOGICAL :: switch
      INTEGER :: i,temp
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      switch= .TRUE.
      DO while(switch)
         switch= .FALSE.
         DO i=mmin,mmax-1
            IF(key(index(i)) < key(index(i+1)))THEN
               temp=index(i)
               index(i)=index(i+1)
               index(i+1)=temp
               switch= .TRUE.
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bubble
c-----------------------------------------------------------------------
c     subprogram 6. my_transpose.
c     transposes a rank-2 array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION my_transpose(a) RESULT(b)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(SIZE(a,2),SIZE(a,1)) :: b

      INTEGER(4) :: i,j
c-----------------------------------------------------------------------
c     compute transpose.
c-----------------------------------------------------------------------
      DO i=1,SIZE(a,1)
         DO j=1,SIZE(a,2)
            b(j,i)=a(i,j)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION my_transpose
c-----------------------------------------------------------------------
c     subprogram 7. my_ctranspose.
c     transposes a rank-2 array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION my_ctranspose(a) RESULT(b)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: a
      COMPLEX(r8), DIMENSION(SIZE(a,2),SIZE(a,1)) :: b

      INTEGER(4) :: i,j
c-----------------------------------------------------------------------
c     compute transpose.
c-----------------------------------------------------------------------
      DO i=1,SIZE(a,1)
         DO j=1,SIZE(a,2)
            b(j,i)=a(i,j)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION my_ctranspose
c-----------------------------------------------------------------------
c     subprogram 8. direct_mypack.
c     computes packed grid on (0,1).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION mypack(nx,pfac,side) RESULT(x)

      INTEGER, INTENT(IN) :: nx
      REAL(r8), INTENT(IN) :: pfac
      CHARACTER(*), INTENT(IN) :: side
      REAL(r8), DIMENSION(-nx:nx) :: x

      INTEGER :: ix
      REAL(r8) :: lambda,a
      REAL(r8), DIMENSION(-nx:nx) :: xi,num,den
c-----------------------------------------------------------------------
c     fill logical grid.
c-----------------------------------------------------------------------
      SELECT CASE(side)
      CASE("left")
         xi=(/(ix,ix=-2*nx,0)/)/REAL(2*nx,r8)
      CASE("right")
         xi=(/(ix,ix=0,2*nx)/)/REAL(2*nx,r8)
      CASE("both")
         xi=(/(ix,ix=-nx,nx)/)/REAL(nx,r8)
      CASE DEFAULT
         CALL program_stop
     $        ("gal_pack: cannot recognize side = "//TRIM(side))
      END SELECT
c-----------------------------------------------------------------------
c     compute grid.
c-----------------------------------------------------------------------
      IF(pfac > 1)THEN
         lambda=SQRT(1-1/pfac)
         num=LOG((1+lambda*xi)/(1-lambda*xi))
         den=LOG((1+lambda)/(1-lambda))
         x=num/den
      ELSEIF(pfac == 1)THEN
         x=xi
      ELSEIF(pfac > 0)THEN
         lambda=SQRT(1-pfac)
         a=LOG((1+lambda)/(1-lambda))
         num=EXP(a*xi)-1
         den=EXP(a*xi)+1
         x=num/den/lambda
      ELSE
         x=xi**(-pfac)
      ENDIF
c-----------------------------------------------------------------------
c     shift output to (-1,1)
c-----------------------------------------------------------------------
      SELECT CASE(side)
      CASE("left")
         x=2*x+1
      CASE("right")
         x=2*x-1
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION mypack
c-----------------------------------------------------------------------
c     subprogram 9. ident.
c     computes and returns integer identity matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION ident(n) RESULT(mat)

      INTEGER, INTENT(IN) :: n
      INTEGER, DIMENSION(n,n) :: mat

      INTEGER :: i,j
c-----------------------------------------------------------------------
c     compute identity matrix.
c-----------------------------------------------------------------------
      mat=RESHAPE((/(1,(0,j=1,n),i=1,n-1),1/),(/n,n/))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION ident

c-----------------------------------------------------------------------
c     subprogram 10. powspace.
c     Taken from pentrc/grid.f
c-----------------------------------------------------------------------
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
        integer :: i
        real(r8), dimension(1:num) :: x
        real(r8) :: deltay,deltax

        ! check input validity
        if(xmax<=xmin)then
           print *,"xmin,xmax = ",xmin,xmax
           stop 'ERROR: powspace - xmax must be less than xmin.'
        endif

        ! endpoint characterization
        select case (endpoints)
           case ("lower")
              x = -1 + (/(i,i=0,num-1)/)/(num-1.0) ! -1 to 0
           case ("upper")
              x = (/(i,i=0,num-1)/)/(num-1.0)  ! 0 to 1
           case ("both")
              x = -1 + 2*(/(i,i=0,num-1)/)/(num-1.0)  ! -1 to 1
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
                powspace(1,:)=x-(4*x**3)/3+(6*x**5)/5-(4*x**7)/7+x**9/9
           case (5)
                powspace(1,:) =-x + (5*x**3)/3 - 2*x**5 +
     $              (10*x**7)/7 - (5*x**9)/9 + x**11/11
           case (6)
                powspace(1,:) =x - 2*x**3 + 3*x**5 - (20*x**7)/7 +
     $              (5*x**9)/3- (6*x**11)/11 + x**13/13
           case (7)
                powspace(1,:) =-x + (7*x**3)/3 - (21*x**5)/5 + 5*x**7 -
     $              (35*x**9)/9 +(21*x**11)/11 -(7*x**13)/13 + x**15/15
           case (8)
                powspace(1,:) =x-(8*x**3)/3+(28*x**5)/5-8*x**7+
     $              (70*x**9)/9- (56*x**11)/11 + (28*x**13)/13 -
     $              (8*x**15)/15+x**17/17
           case (9)
                powspace(1,:) = -x + 3*x**3 - (36*x**5)/5 + 12*x**7
     $              - 14*x**9 + (126*x**11)/11 - (84*x**13)/13
     $              + (12*x**15)/5 - (9*x**17)/17 + x**19/19
           case default
              call program_stop("Grid power not in analytic database")
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
      end function

      END MODULE utils_mod
