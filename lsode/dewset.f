      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
      USE local_mod, ONLY: r8

C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      REAL(r8) RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
 15     CONTINUE
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
 25     CONTINUE
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
 35     CONTINUE
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
 45     CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
