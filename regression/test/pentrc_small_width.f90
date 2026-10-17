! A causal Krook width stays positive when electric precession reverses.
! With wb = wd = wt = 0 and wn = 1 - we, the native energy integral is
! Gamma(7/2)/(i*we - nu), where nu = 1e-5*abs(we).
program pentrc_small_width
    use params, only: r8
    use energy_integration, only: xintgrl_lsode, xatol, xrtol, xnutype, &
        xf0type, ximag, xnufac
    implicit none

    complex(r8) :: integral(2), expected
    real(r8) :: we, width, gamma_seven_halves, relative_error
    integer :: i

    xnutype = "small"
    xf0type = "maxwellian"
    ximag = 0.0_r8
    xnufac = 1.0_r8
    xatol = 1.0e-14_r8
    xrtol = 1.0e-11_r8
    gamma_seven_halves = 15.0_r8*sqrt(acos(-1.0_r8))/8.0_r8

    do i = 1, 2
        we = real(3 - 2*i, r8)
        width = 1.0e-5_r8*abs(we)
        expected = gamma_seven_halves/cmplx(-width, we, kind=r8)
        integral(i) = xintgrl_lsode(1.0_r8-we, 0.0_r8, we, 0.0_r8, &
            0.0_r8, 1.0_r8, 0, 0.0_r8, 1, 0.5_r8, 0.5_r8, "fgar")
        relative_error = abs(real(integral(i), r8)/real(expected, r8)-1.0_r8)
        print '(a,f5.1,a,2es23.14,a,es12.3)', "we = ", we, &
            " integral = ", integral(i), " dissipative error = ", relative_error
        if (real(integral(i), r8) >= 0.0_r8) &
            error stop "small collision width gives anti-damping"
        if (relative_error > 1.0e-7_r8) &
            error stop "small-width dissipative integral differs from analytic oracle"
        if (abs(aimag(integral(i))/aimag(expected)-1.0_r8) > 1.0e-9_r8) &
            error stop "small-width reactive integral differs from analytic oracle"
    end do

    if (abs(integral(2)-conjg(integral(1))) > 1.0e-9_r8*abs(integral(1))) &
        error stop "reversing precession violates the causal reflection symmetry"
end program pentrc_small_width
