      PROGRAM read_equil
      
      IMPLICIT NONE

      INTEGER :: i,j,nw,nh,ia,mr,mz,ma,in_unit
      INTEGER :: ios
      REAL(8) :: bcentr,cpasma,rgrid,rmaxis,rzero,ssibry1,
     $     ssibry2,ssimag1,ssimag2,xdim,xdum,zdim,zmaxis,zmid
      REAL(8), DIMENSION(:,:), POINTER :: sq,psi
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      in_unit=3
      OPEN(UNIT=in_unit,FILE='equil3_conv.txt',STATUS='old')
      READ(in_unit,'(52x,2i4)')nw,nh
      READ(in_unit,'(5e16.9)')xdim,zdim,rzero,rgrid,zmid
      READ(in_unit,'(5e16.9)')rmaxis,zmaxis,ssimag1,ssibry1,bcentr
      READ(in_unit,'(5e16.9)')cpasma,ssimag2,xdum,rmaxis,xdum
      READ(in_unit,'(5e16.9)')zmaxis,xdum,ssibry2,xdum,xdum
      ALLOCATE(sq(0:nw-1,3),psi(0:nw-1,0:nh-1))
      READ(in_unit,'(5e16.9)')(sq(i,1),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq(i,2),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq(i,3),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq(i,3),i=0,nw-1)
      READ(in_unit,'(5e16.9)')((psi(i,j),i=0,nw-1),j=0,nh-1)
      READ(in_unit,'(5e16.9)',iostat=ios)(sq(i,3),i=0,nw-1)
      IF(ios /= 0)sq(i,3)=0
      CLOSE(in_unit)
      WRITE(*,*) "DONE!"
c-----------------------------------------------------------------------
c     translate to internal quantities.
c-----------------------------------------------------------------------      
      END PROGRAM read_equil
