!-----------------------------------------------------------------------
! Unit test for spline_fit_hermite: flagged columns keep the supplied
! knot derivatives; other columns are those of spline_fit.
! Exits with status 1 on failure.
!-----------------------------------------------------------------------
PROGRAM test_spline_hermite
  USE spline_mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: mx=32
  INTEGER :: ix
  REAL(r8), DIMENSION(0:mx) :: slopes
  REAL(r8), DIMENSION(0:mx,2) :: fs1_ref
  TYPE(spline_type) :: spl
  LOGICAL :: ok=.TRUE.

  CALL spline_alloc(spl,mx,2)
  spl%xs=(/(ix,ix=0,mx)/)/REAL(mx,r8)
  spl%fs(:,1)=SIN(3*spl%xs)
  spl%fs(:,2)=spl%fs(:,1)
  CALL spline_fit(spl,"extrap")
  fs1_ref=spl%fs1
  slopes=3*COS(3*spl%xs)
  spl%fs1(:,1)=slopes
  CALL spline_fit_hermite(spl,"extrap",(/.TRUE.,.FALSE./))
  IF(ANY(spl%fs1(:,1) /= slopes))THEN
     WRITE(*,*)"FAIL: flagged column lost the supplied derivatives"
     ok=.FALSE.
  ENDIF
  IF(ANY(spl%fs1(:,2) /= fs1_ref(:,2)))THEN
     WRITE(*,*)"FAIL: column not flagged differs from spline_fit"
     ok=.FALSE.
  ENDIF
  CALL spline_dealloc(spl)
  IF(.NOT.ok)STOP 1
  WRITE(*,*)"PASS"

END PROGRAM test_spline_hermite
