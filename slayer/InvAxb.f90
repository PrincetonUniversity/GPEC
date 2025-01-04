  subroutine InvAxb(a,x,b,n)
implicit none
integer n
!complex*8, dimension(n,n) :: a(n,n), acopy(n,n)
!complex*8, dimension(n) :: b(n), d(n), x(n), InvLb(n)
!complex*8 coeff
double complex :: a(n,n), acopy(n,n)
double complex :: b(n), d(n), x(n), InvLb(n)
double complex :: coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
x = 0.0
InvLb = b

do i=1,n
  do j=1,n
    acopy(i,j)=a(i,j)
  end do
end do

! step 1: forward elimination
do k=1, n-1
  do i=k+1,n
    coeff=acopy(i,k)/acopy(k,k)
    InvLb(i) = InvLb(i)-coeff*InvLb(k)
    !L(i,k) = coeff
    do j=k+1,n
      acopy(i,j) = acopy(i,j)-coeff*acopy(k,j)
    end do
  end do
end do

! Step 2: prepare L and U matrices
!do i=1,n
!  L(i,i) = 1.0
!end do
do j=n,1,-1
  x(j) = InvLb(j)/acopy(j,j)
  do i=j+1,n
    !U(i,j) = acopy(i,j)
    x(j) = x(j)-acopy(j,i)*x(i)/acopy(j,j)
  end do
end do
end subroutine InvAxb
