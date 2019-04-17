!-----------------------------------------------------------------------
!     Finite Element Method matrix management.
!-----------------------------------------------------------------------
MODULE DeltagFEMMod
   USE local_mod
   IMPLICIT NONE

   !define equation term (scalCoef * coefName * varName in equName)
   TYPE :: LocalElementMatrix
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat  ! local FEM matrix mat(basis_i, basis_j, equIndex*mpert, varIndex*mpert)
   END TYPE LocalElementMatrix

   !Manage FEM matrix
   TYPE, PUBLIC :: FiniteElementMethodMatrixManage

      INTEGER :: elementNum = 0                 !number of elements
      INTEGER :: mpert

      TYPE(LocalElementMatrix), DIMENSION (:), ALLOCATABLE :: locMatArray

   CONTAINS
      PROCEDURE :: initFEMManage => deltagInitFEMManage                !Init FEM matrix Manage
      PROCEDURE :: deallocate => deltagDeallocate                   !Deallocate FEM matrix Manage
   END TYPE FiniteElementMethodMatrixManage

   CONTAINS

   SUBROUTINE deltagInitFEMManage (self,equNameList,varNameList,mpert,elementNum)
      CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER,DIMENSION(:,:),ALLOCATABLE, INTENT(IN) :: equNameList, varNameList
      INTEGER, INTENT(IN) :: mpert, elementNum
      INTEGER :: i, equNum,varNum,compNum
      equNum = size(equNameList,1)
      varNum = size(varNameList,1)
      compNum = equNum * varNum * mpert
      self%elementNum = elementNum
      ALLOCATE(self%locMatArray(elementNum))
      DO i = 1, self%elementNum, 1
         ALLOCATE(self%locMatArray(i)%mat(4,4,compNum,compNum))
      ENDDO
   END SUBROUTINE deltagInitFEMManage

   SUBROUTINE deltagDeallocate(self)
      CLASS(FiniteElementMethodMatrixManage), INTENT(INOUT) :: self
      INTEGER :: i

      DO i = 1, self%elementNum, 1
         DEALLOCATE(self%locMatArray(i)%mat)
      ENDDO

      DEALLOCATE(self%locMatArray)
   END SUBROUTINE deltagDeallocate
END MODULE DeltagFEMMod

