! Module   : crmath_r_dp
! Purpose  : Fortran interface to crlibm (complex, double precision)
!
! Copyright 2020 Rich Townsend
!
! This file is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!

submodule (crmath) crmath_c_dp

  ! No implicit typing

  implicit none

  ! Procedures

contains

  elemental complex(DP) function log_c_dp_ (x) result (log_x)

    complex(DP), intent(in) :: x

    log_x = CMPLX(log(abs(x)), atan2(x%im, x%re), DP)

  end function log_c_dp_

  !****

  elemental complex(DP) function log1p_c_dp_ (x) result (log1p_x)

    complex(DP), intent(in) :: x

    complex(DP) :: x1p

    x1p = x + 1._DP

    log1p_x = CMPLX(log(abs(x1p)), atan2(x1p%im, x1p%re), DP)

  end function log1p_c_dp_

  !****

  elemental complex(DP) function exp_c_dp_ (x) result (exp_x)

    complex(DP), intent(in) :: x

    exp_x = exp(x%re)*CMPLX(cos(x%im), sin(x%im), DP)

  end function exp_c_dp_

  !****

  elemental complex(DP) function expm1_c_dp_ (x) result (expm1_x)

    complex(DP), intent(in) :: x

    expm1_x = exp(x) - 1._DP

  end function expm1_c_dp_

  !****

  elemental complex(DP) function sqrt_c_dp_ (x) result (sqrt_x)

    complex(DP), intent(in) :: x

    real(DP) :: d
    real(DP) :: r
    real(DP) :: s

    if (x%im == 0._DP) then

       if (x%re < 0._DP) then

          sqrt_x = CMPLX(0._DP, SIGN(SQRT(-x%re), x%im), DP)

       else

          sqrt_x = CMPLX(SQRT(x%re), SIGN(0._DP, x%im), DP)

       endif

    elseif (x%re == 0._DP) then

       r = SQRT(0.5_DP*ABS(x%im))

       sqrt_x = CMPLX(r, SIGN(r, x%im), DP)

    else

       d = hypot(x%re, x%im)

       if (x%re > 0._DP) then
          r = SQRT(0.5_DP*d + 0.5_DP*x%re)
          s = 0.5_DP*x%im/r
       else
          s = SQRT(0.5_DP*d - 0.5_DP*x%re)
          r = ABS(0.5_DP*x%im/s)
       endif

       sqrt_x = CMPLX(r, SIGN(s, x%im), DP)

    endif

  end function sqrt_c_dp_

  !****

  elemental real(DP) function abs_c_dp_ (x) result (abs_x)

    complex(DP), intent(in) :: x

    abs_x = hypot(x%re, x%im)

  end function abs_c_dp_

  !****

  elemental complex(DP) function cos_c_dp_ (x) result (cos_x)

    complex(DP), intent(in) :: x

    cos_x = CMPLX(cos(x%re)*cosh(x%im), -sin(x%re)*sinh(x%im), DP)

  end function cos_c_dp_

  !****

  elemental complex(DP) function sin_c_dp_ (x) result (sin_x)

    complex(DP), intent(in) :: x

    sin_x = CMPLX(sin(x%re)*cosh(x%im), cos(x%re)*sinh(x%im), DP)

  end function sin_c_dp_

  !****

  elemental complex(DP) function tan_c_dp_ (x) result (tan_x)

    complex(DP), intent(in) :: x

    real(DP) :: rt
    real(DP) :: it

    rt = tan(x%re)
    it = tanh(x%im)

    tan_x = CMPLX(rt, it, DP)/CMPLX(1._DP, -rt*it, DP)

  end function tan_c_dp_

  !****

  elemental complex(DP) function acos_c_dp_ (x) result (acos_x)

    complex(DP), intent(in) :: x

    acos_x = -I_DP*log(x + I_DP*SQRT(1._dp - x**2))

  end function acos_c_dp_

  !****

  elemental complex(DP) function asin_c_dp_ (x) result (asin_x)

    complex(DP), intent(in) :: x

    asin_x = -I_DP*log(I_DP*x + SQRT(1._DP - x**2))

  end function asin_c_dp_

  !****

  elemental complex(DP) function atan_c_dp_ (x) result (atan_x)

    complex(DP), intent(in) :: x

    atan_x = I_DP*log((I_DP + x)/(I_DP - x))/2._DP

  end function atan_c_dp_

  !****

  elemental complex(DP) function cosh_c_dp_ (x) result (cosh_x)

    complex(DP), intent(in) :: x

    cosh_x = CMPLX(cosh(x%re)*cos(x%im), sinh(x%re)*sin(x%im), DP)

  end function cosh_c_dp_

  !****

  elemental complex(DP) function sinh_c_dp_ (x) result (sinh_x)

    complex(DP), intent(in) :: x

    sinh_x = CMPLX(sinh(x%re)*cos(x%im), cosh(x%re)*sin(x%im), DP)

  end function sinh_c_dp_

  !****

  elemental complex(DP) function tanh_c_dp_ (x) result (tanh_x)

    complex(DP), intent(in) :: x

    real(DP) :: rt
    real(DP) :: it

    rt = tanh(x%re)
    it = tan(x%im)

    tanh_x = CMPLX(rt, it, DP)/CMPLX(1._DP, rt*it, DP)

  end function tanh_c_dp_
  
  !****

  elemental complex(DP) function acosh_c_dp_ (x) result (acosh_x)

    complex(DP), intent(in) :: x

    acosh_x = log(x + sqrt(x - 1._DP)*sqrt(x + 1._DP))

  end function acosh_c_dp_

  !****

  elemental complex(DP) function asinh_c_dp_ (x) result (asinh_x)

    complex(DP), intent(in) :: x

    asinh_x = log(x + sqrt(1._DP + x**2))

  end function asinh_c_dp_

  !****

  elemental complex(DP) function atanh_c_dp_ (x) result (atanh_x)

    complex(DP), intent(in) :: x

    real(DP) :: rt
    real(DP) :: it

    atanh_x = log((1._DP + x)/(1._DP - x))/2._DP

  end function atanh_c_dp_

end submodule crmath_c_dp
