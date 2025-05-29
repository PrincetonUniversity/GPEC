! Module   : crmath
! Purpose  : Fortran interface to crlibm
!
! Copyright 2020 Rich Townsend
!
! This file is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
module crmath

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: SP = KIND(0.)
  integer, parameter :: DP = KIND(0.D0)

  complex(SP), parameter :: I_SP = CMPLX(0._SP, 1._SP, SP)
  complex(DP), parameter :: I_DP = CMPLX(0._DP, 1._DP, DP)

  real(SP), parameter :: PI_SP = 3.1415926535897932385_SP
  real(SP), parameter :: TWOPI_SP = 6.2831853071795864769_SP
  real(SP), parameter :: HALFPI_SP = 1.5707963267948966192_SP

  real(DP), parameter :: PI_DP = 3.1415926535897932385_DP
  real(DP), parameter :: TWOPI_DP = 6.2831853071795864769_DP
  real(DP), parameter :: HALFPI_DP = 1.5707963267948966192_DP

  ! Private interfaces (all but the first, implemented in submodules)

  interface

     pure subroutine crlibm_init () bind (C)
       use ISO_C_BINDING
     end subroutine crlibm_init

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

  ! Public interfaces

  interface log
     module procedure log_r_sp_
     module procedure log_c_sp_
     module procedure log_r_dp_
     module procedure log_c_dp_
  end interface log

  interface alog
     module procedure log_r_sp_
  end interface alog

  interface clog
     module procedure log_c_sp_
  end interface clog

  interface dlog
     module procedure log_r_dp_
  end interface dlog

  interface cdlog
     module procedure log_c_dp_
  end interface cdlog

  interface log1p
     module procedure log1p_r_sp_
     module procedure log1p_c_sp_
     module procedure log1p_r_dp_
     module procedure log1p_c_dp_
  end interface log1p

  interface log2
     module procedure log2_r_sp_
     module procedure log2_r_dp_
  end interface log2

  interface log10
     module procedure log10_r_sp_
     module procedure log10_r_dp_
  end interface log10

  interface alog10
     module procedure log10_r_sp_
  end interface alog10

  interface dlog10
     module procedure log10_r_dp_
  end interface dlog10

  interface exp
     module procedure exp_r_sp_
     module procedure exp_c_sp_
     module procedure exp_r_dp_
     module procedure exp_c_dp_
  end interface exp

  interface cexp
     module procedure exp_c_sp_
  end interface cexp

  interface dexp
     module procedure exp_r_dp_
  end interface dexp

  interface cdexp
     module procedure exp_c_dp_
  end interface cdexp

  interface expm1
     module procedure expm1_r_sp_
     module procedure expm1_c_sp_
     module procedure expm1_r_dp_
     module procedure expm1_c_dp_
  end interface expm1

  interface sqrt
     module procedure sqrt_c_sp_
     module procedure sqrt_c_dp_
  end interface sqrt

  interface csqrt
     module procedure sqrt_c_sp_
  end interface csqrt

  interface cdsqrt
     module procedure sqrt_c_dp_
  end interface cdsqrt

  interface abs
     module procedure abs_c_sp_
     module procedure abs_c_dp_
  end interface abs

  interface cabs
     module procedure abs_c_sp_
  end interface cabs
     
  interface cdabs
     module procedure abs_c_dp_
  end interface cdabs

  interface hypot
     module procedure hypot_r_r_sp_
     module procedure hypot_r_r_dp_
  end interface hypot

  interface cos
     module procedure cos_r_sp_
     module procedure cos_c_sp_
     module procedure cos_r_dp_
     module procedure cos_c_dp_
  end interface cos

  interface ccos
     module procedure cos_c_sp_
  end interface ccos

  interface dcos
     module procedure cos_r_dp_
  end interface dcos

  interface cdcos
     module procedure cos_c_dp_
  end interface cdcos

  interface sin
     module procedure sin_r_sp_
     module procedure sin_c_sp_
     module procedure sin_r_dp_
     module procedure sin_c_dp_
  end interface sin

  interface csin
     module procedure sin_c_sp_
  end interface csin

  interface dsin
     module procedure sin_r_dp_
  end interface dsin

  interface cdsin
     module procedure sin_c_dp_
  end interface cdsin

  interface tan
     module procedure tan_r_sp_
     module procedure tan_c_sp_
     module procedure tan_r_dp_
     module procedure tan_c_dp_
  end interface tan

  interface ctan
     module procedure tan_c_sp_
  end interface ctan

  interface dtan
     module procedure tan_r_dp_
  end interface dtan

  interface cdtan
     module procedure tan_c_dp_
  end interface cdtan

  interface cospi
     module procedure cospi_r_sp_
     module procedure cospi_r_dp_
  end interface cospi

  interface sinpi
     module procedure sinpi_r_sp_
     module procedure sinpi_r_dp_
  end interface sinpi

  interface tanpi
     module procedure tanpi_r_sp_
     module procedure tanpi_r_dp_
  end interface tanpi

  interface acos
     module procedure acos_r_sp_
     module procedure acos_c_sp_
     module procedure acos_r_dp_
     module procedure acos_c_dp_
  end interface acos

  interface dacos
     module procedure acos_r_dp_
  end interface dacos

  interface asin
     module procedure asin_r_sp_
     module procedure asin_c_sp_
     module procedure asin_r_dp_
     module procedure asin_c_dp_
  end interface asin

  interface dasin
     module procedure asin_r_dp_
  end interface dasin

  interface atan
     module procedure atan_r_sp_
     module procedure atan_c_sp_
     module procedure atan_r_dp_
     module procedure atan_c_dp_
  end interface atan

  interface datan
     module procedure atan_r_dp_
  end interface datan

  interface atan2
     module procedure atan2_r_r_sp_
     module procedure atan2_r_r_dp_
  end interface atan2

  interface datan2
     module procedure atan2_r_r_dp_
  end interface datan2

  interface acospi
     module procedure acospi_r_sp_
     module procedure acospi_r_dp_
  end interface acospi

  interface asinpi
     module procedure asinpi_r_sp_
     module procedure asinpi_r_dp_
  end interface asinpi

  interface atanpi
     module procedure atanpi_r_sp_
     module procedure atanpi_r_dp_
  end interface atanpi

  interface cosh
     module procedure cosh_r_sp_
     module procedure cosh_c_sp_
     module procedure cosh_r_dp_
     module procedure cosh_c_dp_
  end interface cosh

  interface dcosh
     module procedure cosh_r_dp_
  end interface dcosh

  interface sinh
     module procedure sinh_r_sp_
     module procedure sinh_c_sp_
     module procedure sinh_r_dp_
     module procedure sinh_c_dp_
  end interface sinh

  interface dsinh
     module procedure sinh_r_dp_
  end interface dsinh

  interface tanh
     module procedure tanh_r_sp_
     module procedure tanh_c_sp_
     module procedure tanh_r_dp_
     module procedure tanh_c_dp_
  end interface tanh

  interface dtanh
     module procedure tanh_r_dp_
  end interface dtanh

  interface acosh
     module procedure acosh_r_sp_
     module procedure acosh_c_sp_
     module procedure acosh_r_dp_
     module procedure acosh_c_dp_
  end interface acosh

  interface dacosh
     module procedure acosh_r_dp_
  end interface dacosh

  interface asinh
     module procedure asinh_r_sp_
     module procedure asinh_c_sp_
     module procedure asinh_r_dp_
     module procedure asinh_c_dp_
  end interface asinh

  interface dasinh
     module procedure asinh_r_dp_
  end interface dasinh

  interface atanh
     module procedure atanh_r_sp_
     module procedure atanh_c_sp_
     module procedure atanh_r_dp_
     module procedure atanh_c_dp_
  end interface atanh

  interface datanh
     module procedure atanh_r_dp_
  end interface datanh

  interface dcmplx
     module procedure dcmplx_r_sp_
     module procedure dcmplx_r_r_sp_
     module procedure dcmplx_c_sp_
     module procedure dcmplx_r_dp_
     module procedure dcmplx_r_r_dp_
  end interface dcmplx

  ! Access specifiers

  private

  public :: crmath_init
  public :: log
  public :: clog
  public :: dlog
  public :: cdlog
  public :: log1p
  public :: log2
  public :: log10
  public :: dlog10
  public :: exp
  public :: expm1
  public :: cexp
  public :: dexp
  public :: cdexp
  public :: sqrt
  public :: csqrt
  public :: cdsqrt
  public :: abs
  public :: cabs
  public :: cdabs
  public :: hypot
  public :: cos
  public :: ccos
  public :: dcos
  public :: cdcos
  public :: sin
  public :: csin
  public :: dsin
  public :: cdsin
  public :: tan
  public :: ctan
  public :: dtan
  public :: cdtan
  public :: cospi
  public :: sinpi
  public :: tanpi
  public :: acos
  public :: dacos
  public :: asin
  public :: dasin
  public :: atan
  public :: datan
  public :: atan2
  public :: datan2
  public :: acospi
  public :: asinpi
  public :: atanpi
  public :: cosh
  public :: dcosh
  public :: sinh
  public :: dsinh
  public :: tanh
  public :: dtanh
  public :: acosh
  public :: dacosh
  public :: asinh
  public :: dasinh
  public :: atanh
  public :: datanh
  public :: dcmplx

contains

  subroutine crmath_init ()

    call crlibm_init()

  end subroutine crmath_init

  ! Real, single precision functions

  !****

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

  ! Complex, single precision functions

  !****

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

  ! Real, double precision functions

  !****

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


  ! Complex, double precision functions

  !****

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


end module crmath
