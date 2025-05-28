! Module   : crmath_r_dp
! Purpose  : Fortran interfacr to crlibm (real, double precision)
!
! Copyright 2020 Rich Townsend
!
! This file is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!

submodule (crmath) crmath_r_dp

  ! No implicit typing

  implicit none

  ! C bindings

  interface

     pure real(C_DOUBLE) function log_rz (x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function log_rz

     pure real(C_DOUBLE) function log1p_rz (x) result (log1p_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function log1p_rz

     pure real(C_DOUBLE) function log2_rz (x) result (log2_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in)  :: x
     end function log2_rz

     pure real(C_DOUBLE) function log10_rz (x) result (log10_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function log10_rz

     pure real(C_DOUBLE) function exp_rd (x) result (exp_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function exp_rd

     pure real(C_DOUBLE) function expm1_rz (x) result (expm1_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function expm1_rz

     pure real(C_DOUBLE) function cos_rz (x) result (cos_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function cos_rz

     pure real(C_DOUBLE) function sin_rz (x) result (sin_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function sin_rz

     pure real(C_DOUBLE) function tan_rz (x) result (tan_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function tan_rz

     pure real(C_DOUBLE) function cospi_rz (x) result (cospi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function cospi_rz

     pure real(C_DOUBLE) function sinpi_rz (x) result (sinpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function sinpi_rz

     pure real(C_DOUBLE) function tanpi_rz (x) result (tanpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function tanpi_rz

     pure real(C_DOUBLE) function acos_rd (x) result (acos_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function acos_rd

     pure real(C_DOUBLE) function asin_rz (x) result (asin_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function asin_rz

     pure real(C_DOUBLE) function atan_rz (x) result (atan_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function atan_rz

     pure real(C_DOUBLE) function acospi_rd (x) result (acospi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function acospi_rd

     pure real(C_DOUBLE) function asinpi_rz (x) result (asinpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function asinpi_rz

     pure real(C_DOUBLE) function atanpi_rz (x) result (atanpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function atanpi_rz

     pure real(C_DOUBLE) function cosh_rz (x) result (cosh_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function cosh_rz

     pure real(C_DOUBLE) function sinh_rz (x) result (sinh_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
     end function sinh_rz

  end interface

  ! Procedures

contains

  elemental real(DP) function log_r_dp_ (x) result (log_x)

    real(DP), intent(in) :: x

    log_x = log_rz(x)

  end function log_r_dp_

  !****

  elemental real(DP) function log1p_r_dp_ (x) result (log1p_x)

    real(DP), intent(in) :: x

    log1p_x = log1p_rz(x)

  end function log1p_r_dp_

  !****

  elemental real(DP) function log2_r_dp_ (x) result (log2_x)

    real(DP), intent(in) :: x

    log2_x = log2_rz(x)

  end function log2_r_dp_

  !****

  elemental real(DP) function log10_r_dp_ (x) result (log10_x)

    real(DP), intent(in) :: x

    log10_x = log10_rz(x)

  end function log10_r_dp_

  !****

  elemental real(DP) function exp_r_dp_ (x) result (exp_x)

    real(DP), intent(in) :: x

    exp_x = exp_rd(x)

  end function exp_r_dp_

  !****

  elemental real(DP) function expm1_r_dp_ (x) result (expm1_x)

    real(DP), intent(in) :: x

    expm1_x = expm1_rz(x)

  end function expm1_r_dp_

  !****

  elemental real(DP) function hypot_r_r_dp_ (x, y) result (hypot_xy)

    real(DP), intent(in) :: x
    real(DP), intent(in) :: y

    if (ABS(x) > ABS(y)) then
       hypot_xy = ABS(x)*SQRT(1._DP + (y/x)**2)
    else
       if (ABS(y) == 0._DP) then
          hypot_xy = 0._DP
       else
          hypot_xy = ABS(y)*SQRT(1._DP + (x/y)**2)
       endif
    endif

  end function hypot_r_r_dp_

  !****

  elemental real(DP) function cos_r_dp_ (x) result (cos_x)

    real(DP), intent(in) :: x

    cos_x = cos_rz(x)

  end function cos_r_dp_

  !****

  elemental real(DP) function sin_r_dp_ (x) result (sin_x)

    real(DP), intent(in) :: x

    sin_x = sin_rz(x)

  end function sin_r_dp_

  !****

  elemental real(DP) function tan_r_dp_ (x) result (tan_x)

    real(DP), intent(in) :: x

    tan_x = tan_rz(x)

  end function tan_r_dp_

  !****

  elemental real(DP) function cospi_r_dp_ (x) result (cospi_x)

    real(DP), intent(in) :: x

    cospi_x = cospi_rz(x)

  end function cospi_r_dp_

  !****

  elemental real(DP) function sinpi_r_dp_ (x) result (sinpi_x)

    real(DP), intent(in) :: x

    sinpi_x = sinpi_rz(x)

  end function sinpi_r_dp_

  !****

  elemental real(DP) function tanpi_r_dp_ (x) result (tanpi_x)

    real(DP), intent(in) :: x

    tanpi_x = tanpi_rz(x)

  end function tanpi_r_dp_

  !****

  elemental real(DP) function acos_r_dp_ (x) result (acos_x)

    real(DP), intent(in) :: x

    acos_x = acos_rd(x)

  end function acos_r_dp_

  !****

  elemental real(DP) function asin_r_dp_ (x) result (asin_x)

    real(DP), intent(in) :: x

    asin_x = asin_rz(x)

  end function asin_r_dp_

  !****

  elemental real(DP) function atan_r_dp_ (x) result (atan_x)

    real(DP), intent(in) :: x

    atan_x = atan_rz(x)

  end function atan_r_dp_

  !****

  elemental real(DP) function atan2_r_r_dp_ (y, x) result (atan2_yx)

    real(DP), intent(in) :: y
    real(DP), intent(in) :: x

    if (x > 0._DP) then
       atan2_yx = atan(y/x)
    elseif (y > 0._DP) then
       atan2_yx = HALFPI_DP - atan(x/y)
    elseif (y < 0._DP) then
       atan2_yx = -HALFPI_DP - atan(x/y)
    elseif (x < 0._DP) then
       atan2_yx = atan(y/x) + PI_DP
    endif

  end function atan2_r_r_dp_

  !****

  elemental real(DP) function acospi_r_dp_ (x) result (acospi_x)

    real(DP), intent(in) :: x

    acospi_x = acospi_rd(x)

  end function acospi_r_dp_

  !****

  elemental real(DP) function asinpi_r_dp_ (x) result (asinpi_x)

    real(DP), intent(in) :: x

    asinpi_x = asinpi_rz(x)

  end function asinpi_r_dp_

  !****

  elemental real(DP) function atanpi_r_dp_ (x) result (atanpi_x)

    real(DP), intent(in) :: x

    atanpi_x = atanpi_rz(x)

  end function atanpi_r_dp_

  !****

  elemental real(DP) function cosh_r_dp_ (x) result (cosh_x)

    real(DP), intent(in) :: x

    cosh_x = cosh_rz(x)

  end function cosh_r_dp_

  !****

  elemental real(DP) function sinh_r_dp_ (x) result (sinh_x)

    real(DP), intent(in) :: x

    sinh_x = sinh_rz(x)

  end function sinh_r_dp_

  !****

  elemental real(DP) function tanh_r_dp_ (x) result (tanh_x)

    real(DP), intent(in) :: x

    tanh_x = sinh(x)/cosh(x)

  end function tanh_r_dp_

  !****

  elemental real(DP) function acosh_r_dp_ (x) result (acosh_x)

    real(DP), intent(in) :: x

    acosh_x = log(x + SQRT(x - 1._DP)*SQRT(x + 1._DP))

  end function acosh_r_dp_

  !****

  elemental real(DP) function asinh_r_dp_ (x) result (asinh_x)

    real(DP), intent(in) :: x

    asinh_x = log(x + SQRT(1._DP + x**2))

  end function asinh_r_dp_

  !****

  elemental real(DP) function atanh_r_dp_ (x) result (atanh_x)

    real(DP), intent(in) :: x

    atanh_x = log((1._DP + x)/(1._DP - x))/2._DP

  end function atanh_r_dp_

  !****

  elemental complex(DP) function dcmplx_r_dp_ (x) result (dcmplx_x)

    real(DP), intent(in) :: x

    dcmplx_x = CMPLX(x, KIND=DP)

  end function dcmplx_r_dp_
  
  !****

  elemental complex(DP) function dcmplx_r_r_dp_ (x, y) result (dcmplx_x)

    real(DP), intent(in) :: x
    real(DP), intent(in) :: y

    dcmplx_x = CMPLX(x, y, KIND=DP)

  end function dcmplx_r_r_dp_

end submodule crmath_r_dp
