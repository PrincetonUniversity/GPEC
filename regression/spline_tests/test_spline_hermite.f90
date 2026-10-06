!-----------------------------------------------------------------------
! Unit tests for spline_fit_hermite, spline_fit_pchip and the "extrap" end slopes.
!   1. cubic Hermite with exact knot slopes converges at 4th order.
!   2. "extrap" end slopes converge at 3rd order or better.
!   3. columns not kept are bit-for-bit those of spline_fit.
!   4. spline_fit_pchip matches scipy.interpolate.PchipInterpolator.
! Exits with status 1 on failure.
!-----------------------------------------------------------------------
PROGRAM test_spline_hermite
  USE spline_mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: nres=4, mfine=2001
  INTEGER, DIMENSION(nres), PARAMETER :: mxs=(/16,32,64,128/)
  INTEGER :: ires,ix,mx
  REAL(r8), DIMENSION(nres) :: err_val,err_end
  REAL(r8) :: x,order_val,order_end
  REAL(r8), DIMENSION(:,:), ALLOCATABLE :: fs1_ref
  TYPE(spline_type) :: spl
  LOGICAL :: ok=.TRUE.
  REAL(r8), DIMENSION(10), PARAMETER :: px=(/0._r8,.1_r8,.25_r8,.3_r8,.5_r8, &
       .7_r8,.85_r8,.9_r8,.95_r8,1._r8/), py=(/1._r8,1.05_r8,1.1_r8,1.1_r8, &
       1.3_r8,1.6_r8,2.4_r8,3.5_r8,3.6_r8,3.62_r8/)
  REAL(r8), DIMENSION(5), PARAMETER :: pxe=(/.05_r8,.2_r8,.4_r8,.875_r8,.97_r8/), &
       ! scipy 1.x PchipInterpolator(px,py)(pxe) and its derivative
       pye=(/1.02701577_r8,1.09154154_r8,1.17_r8,2.98681184_r8,3.61184_r8/), &
       pde=(/0.50698198_r8,0.30930931_r8,1.2_r8,29.69419306_r8,0.496_r8/)

  DO ires=1,nres
     mx=mxs(ires)
     CALL spline_alloc(spl,mx,2)
     spl%xs=(/(ix,ix=0,mx)/)/REAL(mx,r8)
     spl%fs(:,1)=f(spl%xs)
     spl%fs(:,2)=f(spl%xs)
     spl%fs1(:,1)=df(spl%xs)
     CALL spline_fit_hermite(spl,"extrap",(/.TRUE.,.FALSE./))
     err_val(ires)=0
     DO ix=0,mfine
        x=ix/REAL(mfine,r8)
        CALL spline_eval(spl,x,0)
        err_val(ires)=MAX(err_val(ires),ABS(spl%f(1)-f1(x)))
     ENDDO
     err_end(ires)=MAX(ABS(spl%fs1(0,2)-df1(0._r8)),ABS(spl%fs1(mx,2)-df1(1._r8)))
     ALLOCATE(fs1_ref(0:mx,2))
     fs1_ref=spl%fs1
     spl%fs(:,2)=f(spl%xs)
     CALL spline_fit(spl,"extrap")
     IF(ANY(spl%fs1(:,2) /= fs1_ref(:,2)))THEN
        WRITE(*,*)"FAIL: column not kept differs from spline_fit, mx = ",mx
        ok=.FALSE.
     ENDIF
     DEALLOCATE(fs1_ref)
     CALL spline_dealloc(spl)
  ENDDO
  order_val=LOG(err_val(nres-1)/err_val(nres))/LOG(2._r8)
  order_end=LOG(err_end(nres-1)/err_end(nres))/LOG(2._r8)
  WRITE(*,'(a,4es10.2,a,f5.2)')" hermite max error:",err_val,", order",order_val
  WRITE(*,'(a,4es10.2,a,f5.2)')" extrap end slope error:",err_end,", order",order_end
  IF(order_val < 3.8_r8)THEN
     WRITE(*,*)"FAIL: Hermite with exact slopes is not 4th order"
     ok=.FALSE.
  ENDIF
  IF(order_end < 2.8_r8)THEN
     WRITE(*,*)"FAIL: extrap end slopes are not 3rd order"
     ok=.FALSE.
  ENDIF
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
  IF(.NOT.ok)STOP 1
  WRITE(*,*)"PASS"

CONTAINS

  ELEMENTAL REAL(r8) FUNCTION f1(x)
    REAL(r8), INTENT(IN) :: x
    f1=SIN(3*x)+EXP(-x)*x**2
  END FUNCTION f1

  ELEMENTAL REAL(r8) FUNCTION df1(x)
    REAL(r8), INTENT(IN) :: x
    df1=3*COS(3*x)+EXP(-x)*(2*x-x**2)
  END FUNCTION df1

  FUNCTION f(x)
    REAL(r8), DIMENSION(:), INTENT(IN) :: x
    REAL(r8), DIMENSION(SIZE(x)) :: f
    f=f1(x)
  END FUNCTION f

  FUNCTION df(x)
    REAL(r8), DIMENSION(:), INTENT(IN) :: x
    REAL(r8), DIMENSION(SIZE(x)) :: df
    df=df1(x)
  END FUNCTION df

END PROGRAM test_spline_hermite
