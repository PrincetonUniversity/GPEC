      SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
      USE local_mod, ONLY: r8

C***BEGIN PROLOGUE  DSRCOM
C***SUBSIDIARY
C***PURPOSE  Save/restore ODEPACK COMMON blocks.
C***LIBRARY   MATHLIB (ODEPACK)
C***TYPE      REAL(r8) (SSRCOM-S, DSRCOM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine saves or restores (depending on JOB) the contents of
C  the COMMON block DLS001, which is used internally
C  by one or more ODEPACK solvers.
C
C  RSAV = real array of length 218 or more.
C  ISAV = integer array of length 37 or more.
C  JOB  = flag indicating to save or restore the COMMON blocks:
C         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
C         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
C         A call with JOB = 2 presumes a prior call with JOB = 1.
C
C***SEE ALSO  LSODE
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   921116  Deleted treatment of block /EH0001/.  (ACH)
C   930801  Reduced Common block length by 2.  (ACH)
C   930809  Renamed to allow single/real(r8) versions. (ACH)
C***END PROLOGUE  DSRCOM
C**End
      INTEGER ISAV, JOB
      INTEGER ILS
      INTEGER I, LENILS, LENRLS
      REAL(r8) RSAV,   RLS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(218), ILS(37)
      DATA LENRLS/218/, LENILS/37/
C
C***FIRST EXECUTABLE STATEMENT  DSRCOM
      IF (JOB .EQ. 2) GO TO 100
C
      DO 10 I = 1,LENRLS
        RSAV(I) = RLS(I)
 10     CONTINUE
      DO 20 I = 1,LENILS
        ISAV(I) = ILS(I)
 20     CONTINUE
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
         RLS(I) = RSAV(I)
 110     CONTINUE
      DO 120 I = 1,LENILS
         ILS(I) = ISAV(I)
 120     CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DSRCOM ----------------------
      END
