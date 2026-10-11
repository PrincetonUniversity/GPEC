program spline_endpoint_polynomials
    use spline_mod, only: r8, spline_type, spline_alloc, spline_fit, &
        spline_eval, spline_int, spline_dealloc
    use cspline_mod, only: cspline_type, cspline_alloc, cspline_fit, &
        cspline_dealloc
    implicit none
    integer :: resolution, degree, m, i
    real(r8) :: value_error, derivative_error, integral_error, s, pi
    real(r8) :: q0, q1, expected_integral
    type(spline_type) :: profile, extended
    type(cspline_type) :: cprofile

    pi = acos(-1.0_r8)
    do resolution = 1, 2
        m = 256 * 2**resolution
        do degree = 0, 3, 3
            call spline_alloc(profile, m, 1)
            profile%xs = 1e-4_r8 + (1.0_r8 - 2e-4_r8) * &
                sin([(real(i, r8) / m, i = 0, m)] * pi / 2)**2
            profile%fs(:, 1) = polynomial(profile%xs, degree)
            call spline_fit(profile, 'extrap')
            derivative_error = maxval(abs(profile%fs1(:, 1) - &
                derivative(profile%xs, degree)))
            call spline_eval(profile, 0.0_r8, 0)
            q0 = profile%f(1)
            call spline_eval(profile, 1.0_r8, 0)
            q1 = profile%f(1)
            value_error = max(abs(q0 - polynomial(0.0_r8, degree)), &
                abs(q1 - polynomial(1.0_r8, degree)))
            call spline_alloc(extended, m + 2, 1)
            extended%xs = [0.0_r8, profile%xs, 1.0_r8]
            extended%fs(:, 1) = [q0, profile%fs(:, 1), q1]
            call spline_fit(extended, 'extrap')
            call spline_int(extended)
            expected_integral = 1.0_r8
            if (degree == 3) expected_integral = 2.0_r8 + 1.5_r8 + &
                4.0_r8 / 3 + 1.25_r8
            integral_error = abs(extended%fsi(m + 2, 1) - expected_integral)
            write(*, '(2i6,3es25.16)') m, degree, derivative_error, &
                value_error, integral_error
            call spline_dealloc(extended)
            call spline_dealloc(profile)
        end do
    end do
    ! complex profile with a large constant offset, as the F, G, K matrices have
    do resolution = 1, 2
        m = 256 * 2**resolution
        call cspline_alloc(cprofile, m, 1)
        cprofile%xs = 1e-4_r8 + (1.0_r8 - 2e-4_r8) * &
            sin([(real(i, r8) / m, i = 0, m)] * pi / 2)**2
        cprofile%fs(:, 1) = cmplx(300.0_r8, -2.0_r8, r8)
        call cspline_fit(cprofile, 'extrap')
        derivative_error = maxval(abs(cprofile%fs1(:, 1)))
        write(*, '(a,i6,es25.16)') 'complex', m, derivative_error
        call cspline_dealloc(cprofile)
    end do
contains
    elemental real(r8) function polynomial(x, degree) result(q)
        real(r8), intent(in) :: x
        integer, intent(in) :: degree
        q = 1.0_r8
        if (degree == 3) q = 2.0_r8 + x * (3.0_r8 + x * (4.0_r8 + 5.0_r8*x))
    end function polynomial

    elemental real(r8) function derivative(x, degree) result(dq)
        real(r8), intent(in) :: x
        integer, intent(in) :: degree
        dq = 0.0_r8
        if (degree == 3) dq = 3.0_r8 + x * (8.0_r8 + 15.0_r8*x)
    end function derivative
end program spline_endpoint_polynomials
