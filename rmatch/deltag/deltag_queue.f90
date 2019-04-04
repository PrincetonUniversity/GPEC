!-----------------------------------------------------------------------
!     Queue to maintain equation terms in general FEM solver.
!-----------------------------------------------------------------------
MODULE DeltagQueueMod
   IMPLICIT NONE

   INTEGER, PARAMETER :: maxQueueSize = 500		!define max size of term queue
   INTEGER, PARAMETER :: nameSize = 10            !max size of variable name

   !define equation term (scalCoef * matCoef * varName in equName)
   TYPE :: QueueElement
      CHARACTER (len = nameSize) :: equName       !equation name
      CHARACTER (len = nameSize) :: varName       !variable name in corresponding equation
      CHARACTER (len = nameSize) :: matCoef       !Coefficient matrix name of term
      INTEGER                    :: varOrder      !variable: order of derivative (0: zeroth order derivative, 1: first order derivative)
      COMPLEX                    :: scalCoef=1.0  !scalar coefficient of term
   END TYPE QueueElement

   !Queue to maintain equation terms
   TYPE, PUBLIC :: DeltagTermsQueue

      INTEGER :: queueSize = 0                 !number of terms maintained by queue 
      TYPE(QueueElement), DIMENSION (maxQueueSize) :: qelement

   CONTAINS
      PROCEDURE :: addTerm => deltagAddTerm    !add terms of equation in queue

   END TYPE DeltagTermsQueue

   CONTAINS

   SUBROUTINE deltagAddTerm (self,equName,varName,matCoef,varOrder,scalCoef)
      CLASS(DeltagTermsQueue), INTENT(INOUT) :: self
      CHARACTER(len=*), INTENT(IN) :: equName, varName, matCoef
      INTEGER, INTENT(IN) :: varOrder
      COMPLEX, OPTIONAL, INTENT(IN) :: scalCoef
      PRINT *, equName
      PRINT *, varName
      self%queueSize=self%queueSize+1
      self%qelement(self%queueSize)%equName=equName
      self%qelement(self%queueSize)%varName=varName
      self%qelement(self%queueSize)%matCoef=matCoef
      self%qelement(self%queueSize)%varOrder=varOrder
      IF (PRESENT(scalCoef)) THEN
         self%qelement(self%queueSize)%scalCoef=scalCoef         
      ENDIF

   END SUBROUTINE deltagAddTerm

END MODULE DeltagQueueMod

PROGRAM test
   USE DeltagQueueMod
   IMPLICIT NONE
   TYPE(DeltagTermsQueue) :: queue
   REAL:: res
   COMPLEX :: cof
   cof=CMPLX(-1,0)
   CALL queue%addTerm('Xi1','B1','G11',0)
   CALL queue%addTerm('Xi1','B2','G22',1,cof)

END PROGRAM test
