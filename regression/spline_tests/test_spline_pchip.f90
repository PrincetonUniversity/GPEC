!-----------------------------------------------------------------------
! Unit tests for spline_fit_pchip.
!   1. matches scipy.interpolate.PchipInterpolator on a pedestal-like table.
!   2. a 2-point table gives the straight line (flat for equal values).
! Exits with status 1 on failure.
!-----------------------------------------------------------------------
PROGRAM test_spline_pchip
  USE spline_mod
  IMPLICIT NONE

  INTEGER :: ix
  REAL(r8) :: x
  REAL(r8), DIMENSION(10), PARAMETER :: px=(/0._r8,.1_r8,.25_r8,.3_r8,.5_r8, &
       .7_r8,.85_r8,.9_r8,.95_r8,1._r8/), py=(/1._r8,1.05_r8,1.1_r8,1.1_r8, &
       1.3_r8,1.6_r8,2.4_r8,3.5_r8,3.6_r8,3.62_r8/)
  REAL(r8), DIMENSION(5), PARAMETER :: pxe=(/.05_r8,.2_r8,.4_r8,.875_r8,.97_r8/), &
       ! scipy 1.x PchipInterpolator(px,py)(pxe) and its derivative
       pye=(/1.02701577_r8,1.09154154_r8,1.17_r8,2.98681184_r8,3.61184_r8/), &
       pde=(/0.50698198_r8,0.30930931_r8,1.2_r8,29.69419306_r8,0.496_r8/)
  TYPE(spline_type) :: spl
  LOGICAL :: ok=.TRUE.

  CALL spline_alloc(spl,9,1)
  spl%xs=px
  spl%fs(:,1)=py
  CALL spline_fit_pchip(spl)
  DO ix=1,5
     CALL spline_eval(spl,pxe(ix),1)
     IF(ABS(spl%f(1)-pye(ix)) > 1e-7_r8 .OR. ABS(spl%f1(1)-pde(ix)) > 1e-6_r8)THEN
        WRITE(*,*)"FAIL: pchip differs from scipy at x = ",pxe(ix),spl%f(1),spl%f1(1)
        ok=.FALSE.
     ENDIF
  ENDDO
  CALL spline_dealloc(spl)
  WRITE(*,*)"pchip vs scipy checked at 5 points"

  CALL spline_alloc(spl,1,2)
  spl%xs=(/0._r8,1._r8/)
  spl%fs(:,1)=1.52_r8
  spl%fs(:,2)=(/1._r8,3._r8/)
  CALL spline_fit_pchip(spl)
  DO ix=0,4
     x=ix/4._r8
     CALL spline_eval(spl,x,1)
     IF(ABS(spl%f(1)-1.52_r8) > 1e-14_r8 .OR. ABS(spl%f1(1)) > 1e-14_r8 &
          .OR. ABS(spl%f(2)-1-2*x) > 1e-14_r8 .OR. ABS(spl%f1(2)-2) > 1e-14_r8)THEN
        WRITE(*,*)"FAIL: 2-point pchip is not the straight line at x = ",x,spl%f,spl%f1
        ok=.FALSE.
     ENDIF
  ENDDO
  CALL spline_dealloc(spl)
  WRITE(*,*)"2-point tables checked at 5 points"
  IF(.NOT.ok)STOP 1
  WRITE(*,*)"PASS"

END PROGRAM test_spline_pchip
