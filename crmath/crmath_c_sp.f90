! Module   : crmath_c_sp
! Purpose  : Fortran interface to crlibm (complex, single precision)
!
! Copyright 2020 Rich Townsend
!
! This file is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!

submodule (crmath) crmath_c_sp

  ! No implicit typing

  implicit none

  ! Procedures

contains

  elemental complex(SP) function log_c_sp_ (x) result (log_x)

    complex(SP), intent(in) :: x

    log_x = CMPLX(log(CMPLX(x, KIND=DP)), KIND=SP)

  end function log_c_sp_

  !****

  elemental complex(SP) function log1p_c_sp_ (x) result (log1p_x)

    complex(SP), intent(in) :: x

    log1p_x = CMPLX(log1p(CMPLX(x, KIND=DP)), KIND=SP)

  end function log1p_c_sp_

  !****

  elemental complex(SP) function exp_c_sp_ (x) result (exp_x)

    complex(SP), intent(in) :: x

    exp_x = CMPLX(exp(CMPLX(x, KIND=DP)), KIND=SP)

  end function exp_c_sp_

  !****

  elemental complex(SP) function expm1_c_sp_ (x) result (expm1_x)

    complex(SP), intent(in) :: x

    expm1_x = CMPLX(expm1(CMPLX(x, KIND=DP)), KIND=SP)

  end function expm1_c_sp_

  !****

  elemental complex(SP) function sqrt_c_sp_ (x) result (sqrt_x)

    complex(SP), intent(in) :: x

    sqrt_x = CMPLX(sqrt(CMPLX(x, KIND=DP)), KIND=SP)

  end function sqrt_c_sp_

  !****

  elemental real(SP) function abs_c_sp_ (x) result (abs_x)

    complex(SP), intent(in) :: x

    abs_x = real(abs(CMPLX(x, KIND=DP)), KIND=SP)

  end function abs_c_sp_

  !****

  elemental complex(SP) function cos_c_sp_ (x) result (cos_x)

    complex(SP), intent(in) :: x

    cos_x = CMPLX(cos(CMPLX(x, KIND=DP)), KIND=SP)

  end function cos_c_sp_

  !****

  elemental complex(SP) function sin_c_sp_ (x) result (sin_x)

    complex(SP), intent(in) :: x

    sin_x = CMPLX(sin(CMPLX(x, KIND=DP)), KIND=SP)

  end function sin_c_sp_

  !****

  elemental complex(SP) function tan_c_sp_ (x) result (tan_x)

    complex(SP), intent(in) :: x

    tan_x = CMPLX(tan(CMPLX(x, KIND=DP)), KIND=SP)

  end function tan_c_sp_

  !****

  elemental complex(SP) function acos_c_sp_ (x) result (acos_x)

    complex(SP), intent(in) :: x

    acos_x = CMPLX(acos(CMPLX(x, KIND=DP)), KIND=SP)

  end function acos_c_sp_

  !****

  elemental complex(SP) function asin_c_sp_ (x) result (asin_x)

    complex(SP), intent(in) :: x

    asin_x = CMPLX(asin(CMPLX(x, KIND=DP)), KIND=SP)

  end function asin_c_sp_

  !****

  elemental complex(SP) function atan_c_sp_ (x) result (atan_x)

    complex(SP), intent(in) :: x

    atan_x = CMPLX(atan(CMPLX(x, KIND=DP)), KIND=SP)

  end function atan_c_sp_

  !****

  elemental complex(SP) function cosh_c_sp_ (x) result (cosh_x)

    complex(SP), intent(in) :: x

    cosh_x = CMPLX(cosh(CMPLX(x, KIND=DP)), KIND=SP)

  end function cosh_c_sp_

  !****

  elemental complex(SP) function sinh_c_sp_ (x) result (sinh_x)

    complex(SP), intent(in) :: x

    sinh_x = CMPLX(sinh(CMPLX(x, KIND=DP)), KIND=SP)

  end function sinh_c_sp_

  !****

  elemental complex(SP) function tanh_c_sp_ (x) result (tanh_x)

    complex(SP), intent(in) :: x

    tanh_x = CMPLX(tanh(CMPLX(x, KIND=DP)), KIND=SP)

  end function tanh_c_sp_

  !****

  elemental complex(SP) function acosh_c_sp_ (x) result (acosh_x)

    complex(SP), intent(in) :: x

    acosh_x = CMPLX(acosh(CMPLX(x, KIND=DP)), KIND=SP)

  end function acosh_c_sp_

  !****

  elemental complex(SP) function asinh_c_sp_ (x) result (asinh_x)

    complex(SP), intent(in) :: x

    asinh_x = CMPLX(asinh(CMPLX(x, KIND=DP)), KIND=SP)

  end function asinh_c_sp_

  !****

  elemental complex(SP) function atanh_c_sp_ (x) result (atanh_x)

    complex(SP), intent(in) :: x

    atanh_x = CMPLX(atanh(CMPLX(x, KIND=DP)), KIND=SP)

  end function atanh_c_sp_

  !****

  elemental complex(DP) function dcmplx_c_sp_ (x) result (dcmplx_x)

    complex(SP), intent(in) :: x

    dcmplx_x = CMPLX(x, KIND=DP)

  end function dcmplx_c_sp_
  
end submodule crmath_c_sp
