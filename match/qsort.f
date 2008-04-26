c-----------------------------------------------------------------------
c     subprogram 1. qsort.
c     quick sort.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. qsort_mod.
c     1. qsort
c     2. compare
c-----------------------------------------------------------------------
c     subprogram 0. qsort_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE qsort_mod
      USE local_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. qsort.
c     quick sort.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE qsort(data,index)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: data
      INTEGER, DIMENSION(:), INTENT(OUT) :: index

      INTEGER, PARAMETER :: m=9
      INTEGER :: l,r,i,j,p,indexp,indext,istk,n
      INTEGER, DIMENSION(31) :: lstk,rstk
      REAL(r8), DIMENSION(SIZE(data,1)) :: datap
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      n=SIZE(data,2)
      index=(/(i,i=1,n)/)
      IF(n <= m)GOTO 900
      istk=0
      l=1
      r=n
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
 200  CONTINUE
      i=l
      j=r
      p=(l+r)/2
      indexp=index(p)
      datap=data(:,indexp)
      IF(compare(data(:,index(l)),datap) > 0)THEN
         index(p)=index(l)
         index(l)=indexp
         indexp=index(p)
         datap=data(:,indexp)
      ENDIF
      IF(compare(datap,data(:,index(r))) > 0)THEN
         IF(compare(data(:,index(l)),data(:,index(r))) > 0)THEN
            index(p)=index(l)
            index(l)=index(r)
         else
            index(p)=index(r)
         ENDIF
         index(r)=indexp
         indexp=index(p)
         datap=data(:,indexp)
      ENDIF
 300  CONTINUE
      i=i+1
      IF(compare(data(:,index(i)),datap) < 0)GOTO 300

      DO
         j=j-1
         IF(compare(data(:,index(j)),datap) <= 0)EXIT
      ENDDO

      IF(i < j)THEN
         indext=index(i)
         index(i)=index(j)
         index(j)=indext
         GOTO 300
      else
         IF(r-j >= i-l .AND. i-l  >  m)THEN
            istk=istk+1
            lstk(istk)=j+1
            rstk(istk)=r
            r=i-1
         else IF(i-l  >  r-j .AND. r-j  >  m)THEN
            istk=istk+1
            lstk(istk)=l
            rstk(istk)=i-1
            l=j+1
         else IF(r-j  >  m)THEN
            l=j+1
         else IF(i-l  >  m)THEN
            r=i-1
         else
            IF(istk < 1)GOTO 900
            l=lstk(istk)
            r=rstk(istk)
            istk=istk-1
         ENDIF
         GOTO 200
      ENDIF
c-----------------------------------------------------------------------
c     primitive sort.
c-----------------------------------------------------------------------
 900  CONTINUE
      DO i=2,n
         IF(compare(data(:,index(i-1)),data(:,index(i))) <= 0)CYCLE
         indexp=index(i)
         datap=data(:,indexp)
         p=i-1
         DO
            index(p+1)=index(p)
            p=p-1
            IF(p <= 0 .OR. compare(data(:,index(p)),datap) <= 0)EXIT
         ENDDO
         index(p+1)=indexp
      ENDDO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qsort
c-----------------------------------------------------------------------
c     subprogram 2. compare.
c     comparison routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      FUNCTION compare(x,y) RESULT(z)

      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8) :: z

      INTEGER :: i,n
c-----------------------------------------------------------------------
c     comparison.
c-----------------------------------------------------------------------
      n=SIZE(x)
      z=0
      i=1
      DO
         IF(x(i) > y(i))THEN
            z=1
            EXIT
         ELSE IF(x(i) < y(i))THEN
            z=-1
            EXIT
         ELSE
            i=i+1
            IF(i > n)EXIT
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION compare
      END MODULE qsort_mod
