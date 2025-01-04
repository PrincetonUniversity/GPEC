  subroutine inverse(a,c,n)
!========================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!----------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!========================================
implicit none
integer n
double complex a(n,n), acopy(n,n), c(n,n)
double complex L(n,n), U(n,n), b(n), d(n), x(n), Ld(n), Ux(n)
double complex coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0d0
U=0.0d0
b=0.0d0

do i=1,n
  do j=1,n
    acopy(i,j)=a(i,j)
  end do
end do

! step 1: forward elimination
do k=1, n-1
  do i=k+1,n
    coeff=acopy(i,k)/acopy(k,k)
    L(i,k) = coeff
    do j=k+1,n
      acopy(i,j) = acopy(i,j)-coeff*acopy(k,j)
    end do
  end do
end do

! Step 2: prepare L and U matrices
! L matrix is a mtrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = acopy(i,j)
  end do
end do

!write (*,401)
!401 format (/,' Lower Triangular Matrix')
!do i=1,n
!  write (*,301) (L(i,j),j=1,n)
!end do
!301 format (6f12.6)

!write (*,402)
!402 format (/,' Upper Triangular Matrix')
!do i=1,n
!  write (*,301) (U(i,j),j=1,n)
!end do
!302 format (6f12.6)

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0d0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i)=d(i)-L(i,j)*d(j)
    end do
  end do

!  write (*,501)
!  501 format (/,' Ld')
  do i=1,n
    Ld(i)=0.0d0
    do j=1,n
      Ld(i)=Ld(i)+L(i,j)*d(j)
    end do
!    write (*,502) Ld(i)
  end do
!  502 format (6f12.6)

! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n) 
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do

!  write (*,601)
!  601 format (/,' d')
!  do i=1,n
!    write (*,602) d(i)
!  end do
!  602 format (6f12.6)

!  write (*,701)
!  701 format (/,' Ux')
  do i=1,n
    Ux(i)=0.0d0
    do j=1,n
      Ux(i)=Ux(i)+U(i,j)*x(j)
    end do
!    write (*,702) Ux(i)
  end do
!  702 format (6f12.6)

! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0d0
end do
end subroutine inverse

