      PROGRAM test_ismath
      USE ismath_mod

      IMPLICIT NONE

      INTEGER :: i,itheta,m,n,mlimit,nlimit
      REAL(8) :: xval,yval,eta,omega
      REAL(8), DIMENSION(0:480) :: theta,thetas
      REAL(8), DIMENSION(-1000:1000) :: x
      REAL(8), DIMENSION(-1000:1000,3) :: y

      REAL(8), DIMENSION(:,:), POINTER :: surfcosmn,surfsinmn

      theta=(/(itheta,itheta=0,480)/)/REAL(480,8)
      thetas=theta**2
      
      DO i=1,11

         yval=0.1*(i-1)
         xval=issect(480,theta,thetas,yval)
         WRITE(*,*) xval,yval,xval**2-yval
         
      ENDDO

      OPEN(UNIT=3, FILE="surfmn.out.cs", STATUS="old")
 24   FORMAT(1x,25f12.6)

      mlimit=64
      nlimit=64
      ALLOCATE(surfcosmn(-mlimit:mlimit,0:nlimit),
     $     surfsinmn(-mlimit:mlimit,0:nlimit))
      DO m=-mlimit,mlimit
         READ(3,24) (surfcosmn(m,n),n=0,nlimit)
         READ(3,24) (surfsinmn(m,n),n=0,nlimit)
      ENDDO
      CLOSE(3)

      DO m=1,30
         WRITE(*,*)surfcosmn(m,1)
      ENDDO
c-----------------------------------------------------------------------
c     quick generation of data sets.
c-----------------------------------------------------------------------
      OPEN(UNIT=2, FILE="fig1_solution1.out", STATUS="unknown")
      omega=1.0
      DO i=-1000,1000
         x(i)=REAL(i,8)/10.0
         eta=0.5
         y(i,1)=eta/(eta**2+(x(i)+omega)**2)+
     $        eta/(eta**2+(x(i)-omega)**2)
         eta=1.0
         y(i,2)=eta/(eta**2+(x(i)+omega)**2)+
     $        eta/(eta**2+(x(i)-omega)**2)
         eta=2.0
         y(i,3)=eta/(eta**2+(x(i)+omega)**2)+
     $        eta/(eta**2+(x(i)-omega)**2) 
         WRITE(2,'(4(2x,es12.3))')x(i),y(i,1),y(i,2),y(i,3)
      ENDDO

      CLOSE(2)
c-----------------------------------------------------------------------
c     end.
c-----------------------------------------------------------------------

      END PROGRAM test_ismath

      
