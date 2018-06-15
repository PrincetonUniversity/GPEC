      REAL*8 FUNCTION DVNORM (N, V, W)
C***BEGIN PROLOGUE  DVNORM
C***SUBSIDIARY
C***PURPOSE  Weighted root-mean-square vector norm.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL*8 (SVNORM-S, DVNORM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This function routine computes the weighted root-mean-square norm
C  of the vector of length N contained in the array V, with weights
C  contained in the array W of length N:
C    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/real*8 versions. (ACH)
C***END PROLOGUE  DVNORM
C**End
      INTEGER N,   I
      REAL*8 V, W,   SUM
      DIMENSION V(N), W(N)
C
C***FIRST EXECUTABLE STATEMENT  DVNORM
      SUM = 0.0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION DVNORM ------------------------
      END
