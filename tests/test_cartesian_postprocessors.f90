! Exercise the production guard with manufactured finite and unavailable fields.
program test_cartesian_postprocessors
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use vacuum_mod, only: check_cartesian_postprocessors
    implicit none
    complex(dp) :: br(3,4), bz(3,4), bp(3,4)
    real(dp) :: unavailable
    character(24) :: mode
    integer :: component

    call get_command_argument(1, mode)
    unavailable = ieee_value(0.0_dp, ieee_quiet_nan)
    br = cmplx(1.0_dp, 2.0_dp, dp)
    bz = cmplx(-3.0_dp, 4.0_dp, dp)
    bp = cmplx(5.0_dp, -6.0_dp, dp)
    call check_cartesian_postprocessors(br, bz, bp, .true., .true.)
    select case (trim(mode))
    case ('chebyshev')
        bz(2,3) = cmplx(-3.0_dp, unavailable, dp)
        call check_cartesian_postprocessors(br, bz, bp, .false., .true.)
        print *, 'UNEXPECTED: Chebyshev accepted an unavailable total field'
    case ('divzero')
        bp(2,3) = cmplx(unavailable, -6.0_dp, dp)
        call check_cartesian_postprocessors(br, bz, bp, .true., .false.)
        print *, 'UNEXPECTED: divzero accepted an unavailable total field'
    case default
        do component = 1, 6
            br = cmplx(1.0_dp, 2.0_dp, dp)
            bz = cmplx(-3.0_dp, 4.0_dp, dp)
            bp = cmplx(5.0_dp, -6.0_dp, dp)
            select case (component)
            case (1)
                br(2,3) = cmplx(unavailable, 2.0_dp, dp)
            case (2)
                br(2,3) = cmplx(1.0_dp, unavailable, dp)
            case (3)
                bz(2,3) = cmplx(unavailable, 4.0_dp, dp)
            case (4)
                bz(2,3) = cmplx(-3.0_dp, unavailable, dp)
            case (5)
                bp(2,3) = cmplx(unavailable, -6.0_dp, dp)
            case (6)
                bp(2,3) = cmplx(5.0_dp, unavailable, dp)
            end select
            call check_cartesian_postprocessors(br, bz, bp, .false., .false.)
        end do
        print *, 'PASS finite postprocessor inputs and raw masked exports'
    end select
end program test_cartesian_postprocessors
