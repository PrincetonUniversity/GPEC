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
      PROCEDURE :: addCoef => deltagAddCoef    !add coefficient into coefficient manage

   END TYPE CoefManage

   CONTAINS

   SUBROUTINE deltagAddCoef (self,coefName)
      CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: coefName
      self%queueSize=self%queueSize+1
      self%coefElement(self%queueSize)%coefName=coefName
   END SUBROUTINE deltagAddCoef

END MODULE DeltagCoeffcientMod

