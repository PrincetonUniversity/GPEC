      REAL*8 FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***LIBRARY   MATHLIB
C***CATEGORY  R1
C***TYPE      REAL*8 (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        REAL*8  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930216  DATE WRITTEN
C   930818  Added SLATEC-format prologue.  (FNF)
C***END PROLOGUE  DUMACH
C
C*Internal Notes:
C-----------------------------------------------------------------------
C Subroutines/functions called by DUMACH.. None
C-----------------------------------------------------------------------
C**End
C
      REAL*8 U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0
 10   U = U*0.5
      COMP = 1.0 + U
      IF (COMP .NE. 1.0) GO TO 10
      DUMACH = U*2.0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
