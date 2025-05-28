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

     ! log

     elemental real(SP) module function log_r_sp_ (x)
       real(SP), intent(in) :: x
     end function log_r_sp_

     elemental complex(SP) module function log_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function log_c_sp_

     elemental real(DP) module function log_r_dp_ (x)
       real(DP), intent(in) :: x
     end function log_r_dp_

     elemental complex(DP) module function log_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function log_c_dp_

     ! log1p

     elemental real(SP) module function log1p_r_sp_ (x)
       real(SP), intent(in) :: x
     end function log1p_r_sp_

     elemental complex(SP) module function log1p_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function log1p_c_sp_

     elemental real(DP) module function log1p_r_dp_ (x)
       real(DP), intent(in) :: x
     end function log1p_r_dp_

     elemental complex(DP) module function log1p_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function log1p_c_dp_

     ! log2

     elemental real(SP) module function log2_r_sp_ (x)
       real(SP), intent(in) :: x
     end function log2_r_sp_

     elemental real(DP) module function log2_r_dp_ (x)
       real(DP), intent(in) :: x
     end function log2_r_dp_

     ! log10

     elemental real(SP) module function log10_r_sp_ (x)
       real(SP), intent(in) :: x
     end function log10_r_sp_

     elemental real(DP) module function log10_r_dp_ (x)
       real(DP), intent(in) :: x
     end function log10_r_dp_

     ! exp

     elemental real(SP) module function exp_r_sp_ (x)
       real(SP), intent(in) :: x
     end function exp_r_sp_

     elemental complex(SP) module function exp_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function exp_c_sp_

     elemental real(DP) module function exp_r_dp_ (x)
       real(DP), intent(in) :: x
     end function exp_r_dp_

     elemental complex(DP) module function exp_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function exp_c_dp_

     ! expm1

     elemental real(SP) module function expm1_r_sp_ (x)
       real(SP), intent(in) :: x
     end function expm1_r_sp_

     elemental complex(SP) module function expm1_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function expm1_c_sp_

     elemental real(DP) module function expm1_r_dp_ (x)
       real(DP), intent(in) :: x
     end function expm1_r_dp_

     elemental complex(DP) module function expm1_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function expm1_c_dp_

     ! sqrt

     elemental complex(SP) module function sqrt_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function sqrt_c_sp_

     elemental complex(DP) module function sqrt_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function sqrt_c_dp_

     ! abs

     elemental real(SP) module function abs_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function abs_c_sp_

     elemental real(DP) module function abs_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function abs_c_dp_

     ! hypot

     elemental real(SP) module function hypot_r_r_sp_ (x, y)
       real(SP), intent(in) :: x
       real(SP), intent(in) :: y
     end function hypot_r_r_sp_
 
     elemental real(DP) module function hypot_r_r_dp_ (x, y)
       real(DP), intent(in) :: x
       real(DP), intent(in) :: y
     end function hypot_r_r_dp_

     ! cos

     elemental real(SP) module function cos_r_sp_ (x)
       real(SP), intent(in) :: x
     end function cos_r_sp_
    
     elemental complex(SP) module function cos_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function cos_c_sp_
    
     elemental real(DP) module function cos_r_dp_ (x)
       real(DP), intent(in) :: x
     end function cos_r_dp_
    
     elemental complex(DP) module function cos_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function cos_c_dp_
    
     ! sin

     elemental real(SP) module function sin_r_sp_ (x)
       real(SP), intent(in) :: x
     end function sin_r_sp_
    
     elemental complex(SP) module function sin_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function sin_c_sp_
    
     elemental real(DP) module function sin_r_dp_ (x)
       real(DP), intent(in) :: x
     end function sin_r_dp_
    
     elemental complex(DP) module function sin_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function sin_c_dp_
    
     ! tan

     elemental real(SP) module function tan_r_sp_ (x)
       real(SP), intent(in) :: x
     end function tan_r_sp_
    
     elemental complex(SP) module function tan_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function tan_c_sp_
    
     elemental real(DP) module function tan_r_dp_ (x)
       real(DP), intent(in) :: x
     end function tan_r_dp_
    
     elemental complex(DP) module function tan_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function tan_c_dp_
    
     ! cospi

     elemental real(SP) module function cospi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function cospi_r_sp_
    
     elemental real(DP) module function cospi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function cospi_r_dp_
    
     ! sinpi

     elemental real(SP) module function sinpi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function sinpi_r_sp_
    
     elemental real(DP) module function sinpi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function sinpi_r_dp_
    
     ! tanpi

     elemental real(SP) module function tanpi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function tanpi_r_sp_
    
     elemental real(DP) module function tanpi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function tanpi_r_dp_
    
     ! acos

     elemental real(SP) module function acos_r_sp_ (x)
       real(SP), intent(in) :: x
     end function acos_r_sp_
    
     elemental complex(SP) module function acos_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function acos_c_sp_
    
     elemental real(DP) module function acos_r_dp_ (x)
       real(DP), intent(in) :: x
     end function acos_r_dp_
    
     elemental complex(DP) module function acos_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function acos_c_dp_
    
     ! asin

     elemental real(SP) module function asin_r_sp_ (x)
       real(SP), intent(in) :: x
     end function asin_r_sp_
    
     elemental complex(SP) module function asin_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function asin_c_sp_
    
     elemental real(DP) module function asin_r_dp_ (x)
       real(DP), intent(in) :: x
     end function asin_r_dp_
    
     elemental complex(DP) module function asin_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function asin_c_dp_
    
     ! atan

     elemental real(SP) module function atan_r_sp_ (x)
       real(SP), intent(in) :: x
     end function atan_r_sp_
    
     elemental complex(SP) module function atan_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function atan_c_sp_
    
     elemental real(DP) module function atan_r_dp_ (x)
       real(DP), intent(in) :: x
     end function atan_r_dp_
    
     elemental complex(DP) module function atan_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function atan_c_dp_

     ! atan2

     elemental real(SP) module function atan2_r_r_sp_ (y, x)
       real(SP), intent(in) :: x
       real(SP), intent(in) :: y
     end function atan2_r_r_sp_
    
     elemental real(DP) module function atan2_r_r_dp_ (y, x)
       real(DP), intent(in) :: x
       real(DP), intent(in) :: y
     end function atan2_r_r_dp_
     
     ! acospi

     elemental real(SP) module function acospi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function acospi_r_sp_
    
     elemental real(DP) module function acospi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function acospi_r_dp_
    
     ! asinpi

     elemental real(SP) module function asinpi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function asinpi_r_sp_
    
     elemental real(DP) module function asinpi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function asinpi_r_dp_
    
     ! atanpi

     elemental real(SP) module function atanpi_r_sp_ (x)
       real(SP), intent(in) :: x
     end function atanpi_r_sp_
    
     elemental real(DP) module function atanpi_r_dp_ (x)
       real(DP), intent(in) :: x
     end function atanpi_r_dp_
    
     ! cosh

     elemental real(SP) module function cosh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function cosh_r_sp_
    
     elemental complex(SP) module function cosh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function cosh_c_sp_
    
     elemental real(DP) module function cosh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function cosh_r_dp_
    
     elemental complex(DP) module function cosh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function cosh_c_dp_
    
     ! sinh

     elemental real(SP) module function sinh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function sinh_r_sp_
    
     elemental complex(SP) module function sinh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function sinh_c_sp_
    
     elemental real(DP) module function sinh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function sinh_r_dp_
    
     elemental complex(DP) module function sinh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function sinh_c_dp_
    
     ! tanh

     elemental real(SP) module function tanh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function tanh_r_sp_
    
     elemental complex(SP) module function tanh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function tanh_c_sp_
    
     elemental real(DP) module function tanh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function tanh_r_dp_
    
     elemental complex(DP) module function tanh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function tanh_c_dp_
    
     ! acosh

     elemental real(SP) module function acosh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function acosh_r_sp_
    
     elemental complex(SP) module function acosh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function acosh_c_sp_
    
     elemental real(DP) module function acosh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function acosh_r_dp_
    
     elemental complex(DP) module function acosh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function acosh_c_dp_
    
     ! asinh

     elemental real(SP) module function asinh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function asinh_r_sp_
    
     elemental complex(SP) module function asinh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function asinh_c_sp_
    
     elemental real(DP) module function asinh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function asinh_r_dp_
    
     elemental complex(DP) module function asinh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function asinh_c_dp_
    
     ! atanh

     elemental real(SP) module function atanh_r_sp_ (x)
       real(SP), intent(in) :: x
     end function atanh_r_sp_
    
     elemental complex(SP) module function atanh_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function atanh_c_sp_
    
     elemental real(DP) module function atanh_r_dp_ (x)
       real(DP), intent(in) :: x
     end function atanh_r_dp_
    
     elemental complex(DP) module function atanh_c_dp_ (x)
       complex(DP), intent(in) :: x
     end function atanh_c_dp_

     ! dcmplx

     elemental complex(DP) module function dcmplx_r_sp_ (x)
       real(SP), intent(in) :: x
     end function dcmplx_r_sp_
    
     elemental complex(DP) module function dcmplx_r_r_sp_ (x, y)
       real(SP), intent(in) :: x
       real(SP), intent(in) :: y
     end function dcmplx_r_r_sp_
    
     elemental complex(DP) module function dcmplx_c_sp_ (x)
       complex(SP), intent(in) :: x
     end function dcmplx_c_sp_
    
     elemental complex(DP) module function dcmplx_r_dp_ (x)
       real(DP), intent(in) :: x
     end function dcmplx_r_dp_
     
     elemental complex(DP) module function dcmplx_r_r_dp_ (x, y)
       real(DP), intent(in) :: x
       real(DP), intent(in) :: y
     end function dcmplx_r_r_dp_
     
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

end module crmath
