!-----------------------------------------------------------------------
!     Queue to maintain equation terms in general FEM solver.
!-----------------------------------------------------------------------
MODULE DeltagTermsQueueMod
   USE local_mod
   IMPLICIT NONE

   INTEGER, PARAMETER :: maxQueueSize = 500		!define max size of term queue
   INTEGER, PARAMETER :: nameSize = 10            !max size of variable name

   !define equation term (scalCoef * coefName * varName in equName)
   TYPE :: QueueElement
      CHARACTER (len = nameSize) :: equName       !equation name
      CHARACTER (len = nameSize) :: varName       !variable name in corresponding equation
      CHARACTER (len = nameSize) :: coefName      !Coefficient matrix name of term
      INTEGER                    :: varOrder      !variable: order of derivative (0: zeroth order derivative, 1: first order derivative)
      COMPLEX(r8)                :: scalCoef=1.0  !scalar coefficient of term
   END TYPE QueueElement

   !Queue to maintain equation terms
   TYPE, PUBLIC :: DeltagTermsQueue

      INTEGER :: queueSize = 0                 !number of terms maintained by queue 
      TYPE(QueueElement), DIMENSION (maxQueueSize) :: qelement

   CONTAINS
      PROCEDURE :: addTerm => deltagAddTerm    !add terms of equation in queue

   END TYPE DeltagTermsQueue

   CONTAINS

   SUBROUTINE deltagAddTerm (self,equName,varName,coefName,varOrder,scalCoef)
      CLASS(DeltagTermsQueue), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: equName, varName, coefName
      INTEGER, INTENT(IN) :: varOrder
      COMPLEX(r8), OPTIONAL, INTENT(IN) :: scalCoef
      self%queueSize=self%queueSize+1
      self%qelement(self%queueSize)%equName=equName
      self%qelement(self%queueSize)%varName=varName
      self%qelement(self%queueSize)%coefName=coefName
      self%qelement(self%queueSize)%varOrder=varOrder
      IF (PRESENT(scalCoef)) THEN
         self%qelement(self%queueSize)%scalCoef=scalCoef         
      ENDIF

   END SUBROUTINE deltagAddTerm

END MODULE DeltagTermsQueueMod

PROGRAM test
   USE DeltagTermsQueueMod
   IMPLICIT NONE
   TYPE(DeltagTermsQueue) :: queue
   REAL(r8) :: res
   COMPLEX(r8) :: cof
   cof=CMPLX(-1,0)
   CALL queue%addTerm('Xi1','B1','G11',0)
   CALL queue%addTerm('Xi1','B2','G22',1,cof)

END PROGRAM test
