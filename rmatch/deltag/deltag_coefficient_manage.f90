!-----------------------------------------------------------------------
!     Manage term coefficients used by general FEM solver.
!-----------------------------------------------------------------------
MODULE DeltagCoeffcientMod
   USE fspline_mod
   IMPLICIT NONE

   INTEGER, PARAMETER :: maxQueueSize = 500		!define max size of term queue
   INTEGER, PARAMETER :: nameSize = 10            !max size of variable name

   !element coefficient spline 
   TYPE :: CoefElement
      CHARACTER (len = nameSize) :: coefName      !coefficient name
      TYPE(fspline_type)         :: coefSpline    !fourier spline of coefficient
   END TYPE CoefElement

   !Queue to maintain equation terms
   TYPE, PUBLIC :: CoefManage

      INTEGER :: queueSize = 0                 !number of coefficients maintained by queue 
      TYPE(CoefElement), DIMENSION (maxQueueSize) :: coefElement

   CONTAINS
      PROCEDURE :: addCoef => deltagAddCoef                !add coefficient into coefficient manage
      PROCEDURE :: findCoefIndex => deltagFindCoefIndex    !find index by coefficient spline

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

   END SUBROUTINE deltagAddCoef

   FUNCTION deltagFindCoefIndex(self,coefName) RESULT(idx)
	CLASS(CoefManage), INTENT(INOUT) :: self
     CHARACTER(len=*), INTENT(IN) :: coefName
     INTEGER :: idx
     INTEGER :: i

     DO i = 1, self.queueSize, 1
        IF (coefName == self%coefElement(i)%coefName ) THEN
           idx=i
           return
        ENDIF
     ENDDO
     idx=-1
   END FUNCTION deltagFindCoefIndex
END MODULE DeltagCoeffcientMod

