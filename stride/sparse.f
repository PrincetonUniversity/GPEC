c-----------------------------------------------------------------------
c     program sparse.
c     defines sparse algebra types and performs algebraic manipulations
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  sparse_mod.
c     1.  sparse_form_array.
c     2.  sparse_form_diag.
c     3.  sparse_copy_array.
c     4.  sparse_solve_ADiag_x_equals_b.
c     5.  sparse_solve_A_x_equals_b.
c     6.  sparse_mult_ADiag_x.
c     7.  sparse_mult_A_x.
c     8.  sparse_heapsort.
c     9.  sparse_maxheapify.
c     10. sparse_whittle_array.
c     11. sparse_whittle_diag.
c     12. sparse_print_array.
c     13. sparse_i2str
c-----------------------------------------------------------------------
c     subprogram 0. sparse_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sparse_mod
      USE local_mod
      IMPLICIT NONE

      TYPE sparse_array
         INTEGER :: n !# of sparse elements
         INTEGER :: nr, nc !# of rows, columns in the dense array
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: a
         INTEGER, DIMENSION(:), ALLOCATABLE :: i
         INTEGER, DIMENSION(:), ALLOCATABLE :: ic
      END TYPE sparse_array

      TYPE sparse_diag
         INTEGER :: n !# of sparse elements
         INTEGER :: nr, nc !# of rows, columns in the dense array
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: a
         INTEGER, DIMENSION(:), ALLOCATABLE :: i
      END TYPE sparse_diag

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1.  sparse_form_array.
c     defines a real 2D sparse matrix--in compressed col and row formats
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_form_array(nsp, n, asp, isp, jsp, b)
         INTEGER, INTENT(IN) :: nsp !# of sparse elements
         INTEGER, INTENT(IN) :: n !size of dense matrix (n x n)
         COMPLEX(r8), DIMENSION(nsp), INTENT(IN) :: asp !values
         INTEGER, DIMENSION(nsp), INTENT(IN) :: isp, jsp !rows, cols
         TYPE(sparse_array), INTENT(OUT) :: b

         LOGICAL, DIMENSION(nsp) :: tanz
         LOGICAL, DIMENSION(:), ALLOCATABLE :: anz
         INTEGER, DIMENSION(:), ALLOCATABLE :: ti,tj,ts,ttj,tti,tempJ,tO
         COMPLEX, DIMENSION(:), ALLOCATABLE :: ta
         INTEGER :: i, j, k, m, nnz
c-----------------------------------------------------------------------
c     check inputs for errors.
c-----------------------------------------------------------------------
         IF(n < MAXVAL(isp) .OR. n < MAXVAL(jsp))THEN
            CALL program_stop("Indices do not fit in sparse matrix.")
         ENDIF
c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
         b%nr = n
         b%nc = n
         tanz = (asp /= 0.0_r8)
         nnz = COUNT(tanz) !nnz = the number of nonzero elements
c-----------------------------------------------------------------------
c     sum over multiply specified elements and remove zeros
c-----------------------------------------------------------------------
         ALLOCATE(ti(nnz), tj(nnz), ta(nnz), ts(nnz), anz(nnz))
         ti = PACK(isp, tanz)
         tj = PACK(jsp, tanz)
ccccccccccccccccccccccccccccccccccccccccccccccccCCCCCCCCCCCCCCCCCCCC
c         ta = PACK(asp, tanz)
c-------------------------------------------------------------------
c     A FORTRAN BUG MADE THIS LOOP MORE RELIABLE THAN THE LINE ABOVE
c-------------------------------------------------------------------
         k=0
         DO i = 1,nsp
            IF(asp(i) /= 0.0_r8)THEN
               k=k+1
               ta(k) = asp(i)
            ENDIF
         ENDDO
         IF (k/=nnz) THEN
            CALL program_stop("PACK proxy loop misfired.")
         ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccCCCCCCCCCCCCCCCCCCCC
c-----------------------------------------------------------------------
c     now that zeros are removed, sort by rows & sum up repeat elements.
c-----------------------------------------------------------------------
         CALL sparse_heapsort(nnz, ti, ts)
         tj = tj(ts)
         ta = ta(ts)
         ALLOCATE(ttj(1),tti(1))
         m = 1
         DO
            k = m
            DO WHILE (ti(k) == ti(m))
               k=k+1
               IF (k > nnz) EXIT
            ENDDO !k is now one ahead of the last index with row=ti(m)
c-----------------------------------------------------------------------
c     sort elements of the same row by columns.
c-----------------------------------------------------------------------
            DEALLOCATE(ttj,tti)
            ALLOCATE(ttj(k-m),tti(k-m))
            ttj = tj(m:k-1)
            CALL sparse_heapsort(k-m,ttj,tti)
            tj(m:k-1) = ttj
            ta(m:k-1) = ta(m+tti-1)
c-----------------------------------------------------------------------
c     now all sorted, so sum elements occurring with equal (i,j).
c-----------------------------------------------------------------------
            DO j = k-1,m+1,-1 !run the loop backward (recall k > m)
               IF ( ti(j-1)==ti(j) .AND. tj(j-1)==tj(j)) THEN
                  ta(j-1) = ta(j-1)+ta(j)
                  ta(j) = 0.0_r8
               ENDIF
            ENDDO
            m = k
            IF (m > nnz) EXIT
         ENDDO
c-----------------------------------------------------------------------
c     re-remove zeros, since zeroed excess (i,j) & form sparse matrix.
c-----------------------------------------------------------------------
         anz = (ta /= 0.0_r8)
         nnz = COUNT(anz)
         b%n = nnz
         ALLOCATE(b%a(nnz),b%i(nnz),b%ic(n+1), tempJ(nnz), tO(nnz))
         b%i = PACK(ti, anz)
         b%a = PACK(ta, anz)
         tempJ = PACK(tj, anz)
         CALL sparse_heapsort(nnz,tempJ,tO) !sort sparse elems by col
         b%i = b%i(tO)
         b%a = b%a(tO)
         b%ic = 0 !b%ic gives the 1st index holding an elem of column j
         b%ic(1) = 1
         DO i = 1,nnz !count the # of sparse elems in each column
            b%ic(tempJ(i)+1) = b%ic(tempJ(i)+1) + 1
         ENDDO
         DO i = 1,n !cumulatively sum these counts-by-column
            b%ic(i+1) =  b%ic(i+1) + b%ic(i)
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_form_array
c-----------------------------------------------------------------------
c     subprogram 2.  sparse_form_diag.
c     defines a real 2D sparse diagonal matrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_form_diag(nsp, n, asp, isp, b)
         INTEGER, INTENT(IN) :: nsp !# of sparse elements
         INTEGER, INTENT(IN) :: n !size of dense matrix (n x n)
         COMPLEX(r8), DIMENSION(nsp), INTENT(IN) :: asp
         INTEGER, DIMENSION(nsp), INTENT(IN) :: isp
         TYPE(sparse_diag), INTENT(OUT) :: b

         LOGICAL, DIMENSION(nsp) :: tanz
         INTEGER, DIMENSION(:), ALLOCATABLE :: ti, ts, tti
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: ta
         LOGICAL, DIMENSION(:), ALLOCATABLE :: anz
         INTEGER :: nnz, m, k
c-----------------------------------------------------------------------
c     check inputs for errors.
c-----------------------------------------------------------------------
         IF(n < MAXVAL(isp))THEN
            CALL program_stop("Indices do not fit sparse diag matrix.")
         ENDIF
c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
         b%nr = n
         b%nc = n
         tanz = (asp /= 0.0_r8)
         nnz = COUNT(tanz)
c-----------------------------------------------------------------------
c     sum over multiply specified elements and remove zeros
c-----------------------------------------------------------------------
         ALLOCATE(ti(nnz), ta(nnz), ts(nnz), anz(nnz))
         ti = PACK(isp, tanz)
         ta = PACK(asp, tanz)
         CALL sparse_heapsort(nnz, ti, ts)
         ta = ta(ts)
         ALLOCATE(tti(1))
         m = 1
         DO WHILE ( m < nnz )
            k = m+1
            DO WHILE (ti(k) == ti(m))
               ta(m) = ta(m)+ta(k)
               ta(k) = 0.0_r8
               k = k+1
               IF (k > nnz) EXIT
            ENDDO
            m = k
         ENDDO
         anz = (ta /= 0.0_r8)
         nnz = COUNT(anz)
         b%n = nnz
         ALLOCATE(b%a(nnz),b%i(nnz))
         b%i = PACK(ti, anz)
         b%a = PACK(ta, anz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_form_diag
c-----------------------------------------------------------------------
c     subprogram 3. sparse_copy_array.
c     copy one sparse real 2D array into another
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_copy_array(B, X)
         TYPE(sparse_array), INTENT(IN) :: B
         TYPE(sparse_array), INTENT(INOUT) :: X
         INTEGER :: m
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
         X%n = B%n
         X%nr = B%nr
         X%nc = B%nc
         m = SIZE(B%i)
         ALLOCATE(X%i(m), X%a(m), X%ic(X%nc+1))
         X%i = B%i
         X%a = B%a
         X%ic = B%ic
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_copy_array
c-----------------------------------------------------------------------
c     subprogram 4. sparse_solve_ADiag_x_equals_b.
c     Solve Ax=B with A real diagonal sparse and B real sparse
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_solve_ADiag_x_equals_b(A,B,X)
         TYPE(sparse_diag), INTENT(IN) :: A
         TYPE(sparse_array), INTENT(IN) :: B
         TYPE(sparse_array), INTENT(OUT) :: X
         INTEGER :: n, k
c-----------------------------------------------------------------------
c     check inputs for errors.
c-----------------------------------------------------------------------
         n = SIZE(A%i)
         IF(A%nr /= A%nc)THEN
            CALL program_stop("SolveARealDiagxEqualsB intended for"//
     $           " square A.")
         ELSEIF(ANY(A%a == 0.0_r8))THEN
            CALL program_stop("Cannot invert a sparse diagonal matrix"//
     $           " with diagonal elements equal to zero.")
         ELSEIF(ANY(A%i /= (/(k,k=1,n)/) ))THEN
            CALL program_stop("Cannot invert a sparse diagonal matrix"//
     $           " with diagonal elements implicitly equal to zero. #1")
         ENDIF
         IF(MAXVAL(B%i) > n)THEN
            CALL program_stop("Cannot invert a sparse diagonal matrix"//
     $           " with diagonal elements implicitly equal to zero. #2")
         ENDIF
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
         CALL sparse_copy_array(B, X)
         X%a = X%a / A%a(X%i) !recall only the diagonal is stored in A!
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_solve_ADiag_x_equals_b
c-----------------------------------------------------------------------
c     subprogram 5.  sparse_solve_A_x_equals_b.
c     SuperLU sparse matrix direct solver with dense rhs vector b.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_solve_A_x_equals_b(a, b, x, runMode)
         TYPE(sparse_array), INTENT(IN) :: a
         COMPLEX(r8), DIMENSION(1:), INTENT(IN) :: b
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: x
         INTEGER, INTENT(IN) :: runMode

         INTEGER, DIMENSION(8), SAVE :: factors
         INTEGER :: info, iopt, ldb, nrhs
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
         ALLOCATE(x(SIZE(b)))
         x = b
         nrhs = 1
c-----------------------------------------------------------------------
c     factor the matrix: runMode=1 has this addt'l step for LU-decomp.
c-----------------------------------------------------------------------
         IF (runMode==1) THEN
            iopt = 1
            CALL c_fortran_zgssv(iopt, a%nc, a%n, nrhs, a%a, a%i, a%ic,
     $           x, a%nr, factors, info)
            IF (info .NE. 0) THEN
               CALL program_stop("Sparse factorize failed: "//
     $              sparse_int2str(info))
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     solve the factored system.
c-----------------------------------------------------------------------
         IF (runMode==1 .OR. runMode==2) THEN
            iopt = 2
            CALL c_fortran_zgssv(iopt, a%nc, a%n, nrhs, a%a, a%i, a%ic,
     $           x, a%nr, factors, info) !x now contains the solution X
            IF (info .NE. 0) THEN
               CALL program_stop("Sparse backsolve failed: "//
     $              sparse_int2str(info))
            ENDIF
c-----------------------------------------------------------------------
c     free memory: runMode=3 frees up this memory.
c-----------------------------------------------------------------------
         ELSEIF (runMode==3) THEN
            iopt = 3
            CALL c_fortran_zgssv(iopt, a%nc, a%n, nrhs, a%a, a%i, a%ic,
     $           x, a%nr, factors, info)
         ELSE
            CALL program_stop("Unknown runMode!")
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_solve_A_x_equals_b
c-----------------------------------------------------------------------
c     subprogram 6.  sparse_mult_ADiag_x.
c     multiplies (sparse matrix a)*(dense vector x) = (dense vector b)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_mult_ADiag_x(a, x, b)
         TYPE(sparse_diag), INTENT(IN) :: a
         COMPLEX(r8), DIMENSION(1:), INTENT(IN) :: x
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::  b

         INTEGER :: j ,k, n
c-----------------------------------------------------------------------
c     initialize and check inputs for errors.
c-----------------------------------------------------------------------
         n = a%nc
         IF (n /= SIZE(x)) THEN
            CALL program_stop("Mismatched size in sparse_mult_ADiag_x!")
         ENDIF
         ALLOCATE(b(n))
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
         b = 0.0_r8
         DO k = 1, a%n
            b(a%i(k)) = b(a%i(k)) + a%a(k) * x(a%i(k))
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_mult_ADiag_x
c-----------------------------------------------------------------------
c     subprogram 7.  sparse_mult_A_x.
c     multiplies (sparse matrix a)*(dense vector x) = (dense vector b)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_mult_A_x(a, x, b)
         TYPE(sparse_array), INTENT(IN) :: a
         COMPLEX(r8), DIMENSION(1:), INTENT(IN) :: x
         COMPLEX(r8), DIMENSION(:), INTENT(OUT) ::  b

         INTEGER :: j ,k, n
c-----------------------------------------------------------------------
c     initialize and check inputs for errors.
c-----------------------------------------------------------------------
         n = a%nc
         IF (n /= SIZE(x)) THEN
            CALL program_stop("Mismatched arrays in sparse_mult_A_x!")
         ENDIF
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
         b = 0.0_r8
         DO j = 1, n
            DO k = a%ic(j), a%ic(j+1)-1
               b(a%i(k)) = b(a%i(k)) + a%a(k) * x(j)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_mult_A_x
c-----------------------------------------------------------------------
c     subprogram 8.  sparse_heapsort.
c     sorts and returns an integer vector in ascending order,
c     (and its original index)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_heapsort(n, a, ai)
         INTEGER, INTENT(IN) :: n  !length of array a
         INTEGER, DIMENSION(n), INTENT(INOUT) :: a
         INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ai
         INTEGER :: i,l,heapsize,vali
         REAL(r8) :: val

         ALLOCATE(ai(n))
         ai = (/(i,i=1,n)/)

         l = n/2
         DO i = l,1,-1
            CALL sparse_maxheapify(a,ai,i,n)
         ENDDO

         heapsize = n
         DO i = n,2,-1
            val = a(1)
            vali = ai(1)
            a(1) = a(i)
            ai(1) = ai(i)
            a(i) = val
            ai(i) = vali

            heapsize = heapsize - 1
            CALL sparse_maxheapify(a,ai,1,heapsize)
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_heapsort
c-----------------------------------------------------------------------
c     subprogram 9.  sparse_maxheapify.
c     a recursive subroutine of the heapsort algorithm
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE sparse_maxheapify(x,xi,i,n)
         INTEGER, DIMENSION(:), INTENT(INOUT) :: x, xi
         INTEGER, INTENT(IN) :: i,n

         INTEGER :: ll, rr, largest, vali
         REAL(r8) :: val

         ll = 2*i
         rr = ll+1
         largest = i
         IF (ll .LE. n) THEN
            IF (x(ll) > x(i)) THEN
               largest = ll
            ENDIF
         ENDIF
         IF (rr .LE. n) THEN
            IF (x(rr) > x(largest)) THEN
               largest = rr
            ENDIF
         ENDIF
         IF (largest /= i) THEN
            val = x(i)
            vali = xi(i)
            x(i) = x(largest)
            xi(i) = xi(largest)
            x(largest) = val
            xi(largest) = vali

            CALL sparse_maxheapify(x,xi,largest,n)
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_maxheapify
c-----------------------------------------------------------------------
c     subprogram 10. sparse_whittle_array.
c     masks parts of a sparse array according to the keep logical vector
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_whittle_array(a, keep, b, idxmap)
         TYPE(sparse_array), INTENT(IN) :: a
         LOGICAL, DIMENSION(1:), INTENT(IN) :: keep
         TYPE(sparse_array), INTENT(OUT) :: b
         INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: idxmap

         INTEGER, DIMENSION(:), ALLOCATABLE :: ai,aj,ai0,aj0,idx1,idx2
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: aa, aa0
         LOGICAL, DIMENSION(:), ALLOCATABLE :: ak0
         INTEGER :: i, j, k, n, na, nnz, keepcount
c-----------------------------------------------------------------------
c     initialize and check for errors.
c-----------------------------------------------------------------------
         n = SIZE(keep)
         IF( n/=a%nr .OR. n/=a%nc ) THEN
            CALL program_stop("sparse_whittle only works on square"//
     $           " matrices, and the masking vector must match.")
         ENDIF
         na = a%n
         ALLOCATE(ai0(na),aj0(na),aa0(na),ak0(na),idx1(na),idx2(na))
         ai0 = a%i
         aa0 = a%a
         DO j=1,a%nc
            aj0(a%ic(j):a%ic(j+1)-1) = j
         ENDDO
         ak0 = .TRUE.
c-----------------------------------------------------------------------
c     generate mapping from originals rows to whittled rows.
c-----------------------------------------------------------------------
         ALLOCATE(idxmap(n,2))
         idxmap = 0
         k = 0
         DO i = 1,n
            idxmap(i,1) = i
            IF (keep(i)) THEN
               k = k+1
               idxmap(i,2) = k
            ENDIF
         ENDDO
         keepcount = k
c-----------------------------------------------------------------------
c     whittle.
c-----------------------------------------------------------------------
         CALL sparse_heapsort(na, aj0, idx1)
         ai0 = ai0(idx1)
         aa0 = aa0(idx1)
         i = 1
         DO k = 1,n
            IF (.NOT. keep(k)) THEN
               DO WHILE (aj0(i) .LE. k)
                  IF (aj0(i)==k) ak0(i) = .FALSE.
                  i = i+1
                  IF (i > na) EXIT
               ENDDO
            ENDIF
            IF (i > na) EXIT
         ENDDO
         CALL sparse_heapsort(na, ai0, idx2)
         aj0 = aj0(idx2)
         aa0 = aa0(idx2)
         ak0 = ak0(idx2)
         i = 1
         DO k = 1,n
            IF (.NOT. keep(k)) THEN
               DO WHILE (ai0(i) .LE. k)
                  IF (ai0(i)==k) ak0(i) = .FALSE.
                  i = i+1
                  IF (i > na) EXIT
               ENDDO
            ENDIF
            IF (i > na) EXIT
         ENDDO
         nnz = COUNT(ak0)
         ALLOCATE(ai(nnz),aj(nnz),aa(nnz))
         ai = PACK(ai0, ak0)
         aj = PACK(aj0, ak0)
         aa = PACK(aa0, ak0)
c-----------------------------------------------------------------------
c     relabel original rows as the whittled rows
c-----------------------------------------------------------------------
         DO i = 1,nnz
            ai(i) = idxmap(ai(i),2)
            aj(i) = idxmap(aj(i),2)
         ENDDO
         CALL sparse_form_array(nnz, keepcount, aa, ai, aj, b)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_whittle_array
c-----------------------------------------------------------------------
c     subprogram 11. sparse_whittle_diag.
c     masks parts of a sparse array according to the keep logical vector
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_whittle_diag(a, keep, b, idxmap)
         TYPE(sparse_diag), INTENT(IN) :: a
         LOGICAL, DIMENSION(1:), INTENT(IN) :: keep
         TYPE(sparse_diag), INTENT(OUT) :: b
         INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: idxmap

         INTEGER, DIMENSION(:), ALLOCATABLE :: ai,ai0,idx1
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: aa, aa0
         LOGICAL, DIMENSION(:), ALLOCATABLE :: ak0
         INTEGER :: i, k, n, na, nnz, keepcount
c-----------------------------------------------------------------------
c     initialize and check for errors.
c-----------------------------------------------------------------------
         n = SIZE(keep)
         IF( n/=a%nr .OR. n/=a%nc ) THEN
            CALL program_stop("sparse_whittle only works on square"//
     $           " diag matrices, and the masking vector must match.")
         ENDIF
         na = a%n
         ALLOCATE(ai0(na),aa0(na),ak0(na),idx1(na))
         ai0 = a%i
         aa0 = a%a
         ak0 = .TRUE.
c-----------------------------------------------------------------------
c     generate mapping from originals rows to whittled rows.
c-----------------------------------------------------------------------
         ALLOCATE(idxmap(n,2))
         idxmap = 0
         k = 0
         DO i = 1,n
            idxmap(i,1) = i
            IF (keep(i)) THEN
               k = k+1
               idxmap(i,2) = k
            ENDIF
         ENDDO
         keepcount = k
c-----------------------------------------------------------------------
c     whittle.
c-----------------------------------------------------------------------
         CALL sparse_heapsort(na, ai0, idx1)
         aa0 = aa0(idx1)
         i = 1
         DO k = 1,n
            IF (.NOT. keep(k)) THEN
               DO WHILE (ai0(i) .LE. k)
                  IF (ai0(i)==k) ak0(i) = .FALSE.
                  i = i+1
                  IF (i > na) EXIT
               ENDDO
            ENDIF
            IF (i > na) EXIT
         ENDDO
         nnz = COUNT(ak0)
         ALLOCATE(ai(nnz),aa(nnz))
         ai = PACK(ai0, ak0)
         aa = PACK(aa0, ak0)
c-----------------------------------------------------------------------
c     relabel original rows as the whittled rows
c-----------------------------------------------------------------------
         DO i = 1,nnz
            ai(i) = idxmap(ai(i),2)
         ENDDO
         CALL sparse_form_diag(nnz, keepcount, aa, ai, b)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sparse_whittle_diag
c-----------------------------------------------------------------------
c     subprogram 12. sparse_print_array.
c     print the sparse array
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sparse_print_array(x, aname, printstyle)
         TYPE(sparse_array), INTENT(IN) :: x
         CHARACTER(LEN=*), INTENT(IN) :: aname
         INTEGER, INTENT(IN) :: printstyle
         INTEGER :: i,kcol
 10      FORMAT(I5,A,I5,A,I5,A,E14.7)
 20      FORMAT(I5,A,I3,A,I3,A,E14.7)
c-----------------------------------------------------------------------
c     print.
c-----------------------------------------------------------------------
         IF (printstyle == 1) THEN
            print *,aname,"%n=",x%n
            print *,aname,"%nr=",x%nr
            print *,aname,"%nc=",x%nc
            print *,aname,"%ic=",x%ic
            print *,aname,"%i=",x%i
            print *,aname,"%ac=",x%a
         ELSEIF (printstyle == 2) THEN
            kcol = 1
            DO i = 1, x%n
               IF (i .GE. x%ic(kcol+1)) THEN
                  kcol=kcol+1
               ENDIF

               IF (x%nr > 999 .OR. x%nc > 999) THEN
                  WRITE(*, 10) i,"  (",x%i(i),",",kcol,")   ",x%a(i)
               ELSE
                  WRITE(*, 20) i,"  (",x%i(i),",",kcol,")   ",x%a(i)
               ENDIF
            ENDDO
         ELSE
            CALL program_stop("Unknown printstyle in sparse_print!")
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE sparse_print_array
c-----------------------------------------------------------------------
c     function 13. sparse_int2str
c     converts integers to strings
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sparse_int2str(num) RESULT(res)
         CHARACTER(:),ALLOCATABLE :: res
         INTEGER, INTENT(in) :: num
         CHARACTER(RANGE(num)+2) :: tmp

         WRITE(tmp,'(i0)') num
         res = trim(tmp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION sparse_int2str

      END MODULE sparse_mod
