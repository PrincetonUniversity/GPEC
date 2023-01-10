!-----------------------------------------------------------------------
!     Finite Element Method matrix management.
!-----------------------------------------------------------------------
MODULE DeltagFEMMod
   USE local_mod
   IMPLICIT NONE

   !define equation term (scalCoef * coefName * varName in equName)
   TYPE :: LocalElementMatrix
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat  ! local FEM matrix mat(basis_i, basis_j, equIndex*mpert, varIndex*mpert)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: map  ! maping between local and global matrix, location index of variables, map(basis_i, varIndex * mpert)
      CHARACTER, DIMENSION(:,:), ALLOCATABLE :: varNameList  !save variable list
      INTEGER, DIMENSION(:), ALLOCATABLE :: varMpert         !save mpert for each vairable
      REAL(r8), DIMENSION(2) :: x                            !start and end point
   END TYPE LocalElementMatrix

   !Manage FEM matrix
   TYPE, PUBLIC :: FiniteElementMethodMatrixManage

      INTEGER :: elementNum = 0                 !number of elements
      INTEGER :: mpert

      TYPE(LocalElementMatrix), DIMENSION (:), ALLOCATABLE :: locMats

   CONTAINS
      PROCEDURE :: initFEMManage => deltagInitFEMManage                !Init FEM matrix Manage
      PROCEDURE :: deallocate => deltagDeallocate                   !Deallocate FEM matrix Manage
   END TYPE FiniteElementMethodMatrixManage

   CONTAINS

   SUBROUTINE deltagInitFEMManage (self,elementNum)
      CLASS(CoefManage), INTENT(INOUT) :: self
      INTEGER, INTENT(IN) :: mpert, elementNum
      self%elementNum = elementNum
      ALLOCATE(self%locMats(elementNum))
   END SUBROUTINE deltagInitFEMManage

   SUBROUTINE deltagInitNormalElement (self,eleidx,varNameList,mpert)
      CLASS(CoefManage), INTENT(INOUT) :: self
      CHARACTER,DIMENSION(:,:),ALLOCATABLE, INTENT(IN) :: varNameList
      INTEGER, INTENT(IN) :: eleidx, mpert
      INTEGER :: equNum,varNum,compNum
      equNum = size(equNameList,1)
      varNum = size(varNameList,1)
      compNum = varNum * mpert
      ALLOCATE(self%locMatArray(i)%mat(4,4,compNum,compNum))
      ALLOCATE(self%locMatArray(i)%map(4,4,compNum,compNum))
   END SUBROUTINE deltagInitElement

   SUBROUTINE deltagDeallocate(self)
      CLASS(FiniteElementMethodMatrixManage), INTENT(INOUT) :: self
      INTEGER :: i

      DO i = 1, self%elementNum, 1
         DEALLOCATE(self%locMats(i)%mat)
         DEALLOCATE(self%locMats(i)%map)
      ENDDO

      DEALLOCATE(self%locMatArray)
   END SUBROUTINE deltagDeallocate
END MODULE DeltagFEMMod

