!-----------------------------------------------------------------------
!     Manage term coefficients used by general FEM solver.
!-----------------------------------------------------------------------
MODULE DeltagCoeffcientMod
   USE fspline_mod
   IMPLICIT NONE

   INTEGER, PRIVATE, PARAMETER :: maxQueueSize = 500	    !define max size of term queue
   INTEGER, PRIVATE, PARAMETER :: nameSize = 10            !max size of variable name
   LOGICAL, PRIVATE, PARAMETER :: fft_flag = .FALSE.       !choose type of fft
   !element coefficient spline 
   TYPE :: CoefElement
      CHARACTER (len = nameSize) :: coefName      !coefficient name
      INTEGER ::                    mband         !band of fourier spectrum
      TYPE(fspline_type)         :: coefSpline    !fourier spline of coefficient
   END TYPE CoefElement

   !Queue to maintain equation terms
   TYPE, PUBLIC :: CoefManage

      INTEGER :: queueSize = 0                 !number of coefficients maintained by queue 
      TYPE(CoefElement), DIMENSION (maxQueueSize) :: coefElement

   CONTAINS
      PROCEDURE :: addCoef => deltagAddCoef                !add coefficient into coefficient manage
      PROCEDURE :: findCoefIndex => deltagFindCoefIndex    !find index by coefficient name
      PROCEDURE :: getCoefMat => deltagGetCoefMat          !return coeefficient matrix with mlow and mhigh
   END TYPE CoefManage

   CONTAINS

   SUBROUTINE deltagAddCoef (self,coefName,psi,theta,values,mband)
      CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: coefName
      REAL(r8), DIMENSION(:), INTENT(IN) :: psi,theta  ! psi [0,1] and theta [0,2pi]
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: values
      INTEGER, INTENT(IN) :: mband
      INTEGER :: mpsi,mtheta

      self%queueSize=self%queueSize+1
      self%coefElement(self%queueSize)%coefName=coefName
      mpsi=SIZE(psi)
      mtheta=SIZE(theta)
      !allocate fourier spline
      CALL fspline_alloc( self%coefElement(self%queueSize)%coefSpline, mpsi, mtheta, mband, 1)
      self%coefElement(self%queueSize)%coefSpline%xtitle='psi'
      self%coefElement(self%queueSize)%coefSpline%ytitle='theta'
      self%coefElement(self%queueSize)%coefSpline%xs=psi
      self%coefElement(self%queueSize)%coefSpline%ys=theta
      self%coefElement(self%queueSize)%coefSpline%fs(:,:,1)=values
      self%coefElement(self%queueSize)%mband=mband

      IF (fft_flag) THEN
         CALL fspline_fit_2(self%coefElement(self%queueSize)%coefSpline, "extrap",.TRUE.)
      ELSE
         CALL fspline_fit_1(self%coefElement(self%queueSize)%coefSpline, "extrap",.TRUE.)
      ENDIF
   END SUBROUTINE deltagAddCoef

   FUNCTION deltagFindCoefIndex(self,coefName) RESULT(idx)
	CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: coefName
      INTEGER :: idx
      INTEGER :: i

      DO i = 1, self%queueSize, 1
         IF (coefName == self%coefElement(i)%coefName ) THEN
            idx=i
            return
         ENDIF
      ENDDO
      idx=-1
   END FUNCTION deltagFindCoefIndex

   SUBROUTINE deltagGetCoefMat(self, coefName, psi, mlow,mhigh,mat)
      CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: coefName
      REAL(r8) :: psi
      INTEGER, INTENT(IN) :: mlow, mhigh
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: mat

      INTEGER :: cofIdx,ipert,jpert,mpert,mband,m1,dm
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: harmo,spf

      cofIdx = self%findCoefIndex(coefName)
      mband=self%coefElement(cofIdx)%mband

      ALLOCATE(harmo(-mband:mband),spf(mband))

      CALL cspline_eval( self%coefElement(cofIdx)%coefSpline%cs,psi,0 )
      spf = self%coefElement(cofIdx)%coefSpline%f
      !from mband to 0
      harmo(0:-mband:-1)=spf
      !from 0 to -mband
      harmo(1:mband)=CONJG(spf)
      !create matrix for convolution
      mat=0.
      mpert=mhigh-mlow+1
      ipert=0
      DO m1=mlow,mhigh
         ipert=ipert+1
         DO dm=MAX(1-ipert,-mband), MIN(mpert-ipert,mband)
            jpert=ipert+dm
            mat(ipert,jpert) = harmo(dm)
         ENDDO
      ENDDO
      
      DEALLOCATE(harmo,mband) 

   END SUBROUTINE deltagGetCoefMat
END MODULE DeltagCoeffcientMod

