      program test_simple
      
      implicit none

      integer :: i,j,k

      character(3) :: si,sii,sj,sk

      i=2
      j=12
      k=-3

      IF (i<10) THEN 
         WRITE(unit=si,fmt='(I1)')i
      ELSE
         WRITE(unit=si,fmt='(I2)')i  
      ENDIF

      WRITE(unit=sj,fmt='(I2)')j
      WRITE(unit=sk,fmt='(I2)')k

      WRITE(*,*)si//'>'
      WRITE(*,*)TRIM(si)//'>'
      WRITE(*,*)sj//'>'
      WRITE(*,*)sk//'>'

      end program test_simple

      
