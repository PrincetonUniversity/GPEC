! Module   : crmath_r_sp
! Purpose  : Fortran interface to crlibm (real, single precision)
!
! Copyright 2020 Rich Townsend
!
! This file is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!

submodule (crmath) crmath_r_sp

  ! No implicit typing

  implicit none

  ! Procedures

contains

  elemental real(SP) function log_r_sp_ (x) result (log_x)

    real(SP), intent(in) :: x

    log_x = REAL(log(REAL(x, DP)), SP)

  end function log_r_sp_

  !****

  elemental real(SP) function log1p_r_sp_ (x) result (log1p_x)

    real(SP), intent(in) :: x

    log1p_x = REAL(log1p(REAL(x, DP)), SP)

  end function log1p_r_sp_

  !****

  elemental real(SP) function log2_r_sp_ (x) result (log2_x)

    real(SP), intent(in) :: x

    log2_x = REAL(log2(REAL(x, DP)), SP)

  end function log2_r_sp_

  !****

  elemental real(SP) function log10_r_sp_ (x) result (log10_x)

    real(SP), intent(in) :: x

    log10_x = REAL(log10(REAL(x, DP)), SP)

  end function log10_r_sp_

  !****

  elemental real(SP) function exp_r_sp_ (x) result (exp_x)

    real(SP), intent(in) :: x

    exp_x = REAL(exp(REAL(x, DP)), SP)

  end function exp_r_sp_

  !****

  elemental real(SP) function expm1_r_sp_ (x) result (expm1_x)

    real(SP), intent(in) :: x

    expm1_x = REAL(expm1(REAL(x, DP)), SP)

  end function expm1_r_sp_

  !****

  elemental real(SP) function hypot_r_r_sp_ (x, y) result (hypot_xy)

    real(SP), intent(in) :: x
    real(SP), intent(in) :: y

    if (ABS(x) > ABS(y)) then
       hypot_xy = ABS(x)*SQRT(1._SP + (y/x)**2)
    else
       if (ABS(y) == 0._SP) then
          hypot_xy = 0._SP
       else
          hypot_xy = ABS(y)*SQRT(1._SP + (x/y)**2)
       endif
    endif

  end function hypot_r_r_sp_

  !****

  elemental real(SP) function cos_r_sp_ (x) result (cos_x)

    real(SP), intent(in) :: x

    cos_x = REAL(cos(REAL(x, DP)), SP)

  end function cos_r_sp_

  !****

  elemental real(SP) function sin_r_sp_ (x) result (sin_x)

    real(SP), intent(in) :: x

    sin_x = REAL(sin(REAL(x, DP)), SP)

  end function sin_r_sp_

  !****

  elemental real(SP) function tan_r_sp_ (x) result (tan_x)

    real(SP), intent(in) :: x

    tan_x = REAL(tan(REAL(x, DP)), SP)

  end function tan_r_sp_

  !****

  elemental real(SP) function cospi_r_sp_ (x) result (cospi_x)

    real(SP), intent(in) :: x

    cospi_x = REAL(cospi(REAL(x, DP)), SP)

  end function cospi_r_sp_

  !****

  elemental real(SP) function sinpi_r_sp_ (x) result (sinpi_x)

    real(SP), intent(in) :: x

    sinpi_x = REAL(sinpi(REAL(x, DP)), SP)

  end function sinpi_r_sp_

  !****

  elemental real(SP) function tanpi_r_sp_ (x) result (tanpi_x)

    real(SP), intent(in) :: x

    tanpi_x = REAL(tanpi(REAL(x, DP)), SP)

  end function tanpi_r_sp_

  !****

  elemental real(SP) function acos_r_sp_ (x) result (acos_x)

    real(SP), intent(in) :: x

    acos_x = REAL(acos(REAL(x, DP)), SP)

  end function acos_r_sp_

  !****

  elemental real(SP) function asin_r_sp_ (x) result (asin_x)

    real(SP), intent(in) :: x

    asin_x = REAL(asin(REAL(x, DP)), SP)

  end function asin_r_sp_

  !****

  elemental real(SP) function atan_r_sp_ (x) result (atan_x)

    real(SP), intent(in) :: x

    atan_x = REAL(atan(REAL(x, DP)), SP)

  end function atan_r_sp_

  !****

  elemental real(SP) function atan2_r_r_sp_ (y, x) result (atan2_yx)

    real(SP), intent(in) :: y
    real(SP), intent(in) :: x

    atan2_yx = REAL(atan2(REAL(y, DP), REAL(x, DP)), SP)

  end function atan2_r_r_sp_

  !****

  elemental real(SP) function acospi_r_sp_ (x) result (acospi_x)

    real(SP), intent(in) :: x

    acospi_x = REAL(acospi(REAL(x, DP)), SP)

  end function acospi_r_sp_

  !****

  elemental real(SP) function asinpi_r_sp_ (x) result (asinpi_x)

    real(SP), intent(in) :: x

    asinpi_x = REAL(asinpi(REAL(x, DP)), SP)

  end function asinpi_r_sp_

  !****

  elemental real(SP) function atanpi_r_sp_ (x) result (atanpi_x)

    real(SP), intent(in) :: x

    atanpi_x = REAL(atanpi(REAL(x, DP)), SP)

  end function atanpi_r_sp_

  !****

  elemental real(SP) function cosh_r_sp_ (x) result (cosh_x)

    real(SP), intent(in) :: x

    cosh_x = REAL(cosh(REAL(x, DP)), SP)

  end function cosh_r_sp_

  !****

  elemental real(SP) function sinh_r_sp_ (x) result (sinh_x)

    real(SP), intent(in) :: x

    sinh_x = REAL(sinh(REAL(x, DP)), SP)

  end function sinh_r_sp_

  !****

  elemental real(SP) function tanh_r_sp_ (x) result (tanh_x)

    real(SP), intent(in) :: x

    tanh_x = REAL(tanh(REAL(x, DP)), SP)

  end function tanh_r_sp_

  !****

  elemental real(SP) function acosh_r_sp_ (x) result (acosh_x)

    real(SP), intent(in) :: x

    acosh_x = REAL(acosh(REAL(x, DP)), SP)

  end function acosh_r_sp_

  !****

  elemental real(SP) function asinh_r_sp_ (x) result (asinh_x)

    real(SP), intent(in) :: x

    asinh_x = REAL(asinh(REAL(x, DP)), SP)

  end function asinh_r_sp_

  !****

  elemental real(SP) function atanh_r_sp_ (x) result (atanh_x)

    real(SP), intent(in) :: x

    atanh_x = REAL(atanh(REAL(x, DP)), SP)

  end function atanh_r_sp_

  !****

  elemental complex(DP) function dcmplx_r_sp_ (x) result (dcmplx_x)

    real(SP), intent(in) :: x

    dcmplx_x = CMPLX(x, KIND=DP)

  end function dcmplx_r_sp_
  
  !****

  elemental complex(DP) function dcmplx_r_r_sp_ (x, y) result (dcmplx_x)

    real(SP), intent(in) :: x
    real(SP), intent(in) :: y

    dcmplx_x = CMPLX(x, y, KIND=DP)

  end function dcmplx_r_r_sp_

end submodule crmath_r_sp
