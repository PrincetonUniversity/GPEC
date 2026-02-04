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

      END MODULE utils_mod