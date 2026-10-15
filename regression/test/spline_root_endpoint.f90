! A monotone Hermite cell from a PENTRC trapped orbit. Its unique root lies
! just inside the right endpoint; a rounded Cardano exponent can move the
! initial candidate outside the cell and silently discard the bounce point.
program spline_root_endpoint
    use local_mod, only: r8
    use spline_mod, only: spline_type, spline_alloc, spline_dealloc, spline_roots
    implicit none
    type(spline_type) :: spl
    real(r8), parameter :: left = .767578125_r8, right = .76953125_r8
    real(r8), parameter :: fleft = -1.99139065907738555e-3_r8
    real(r8), parameter :: fright = 1.32122017770797129e-8_r8
    real(r8), parameter :: dleft = 1.01883408452285917_r8
    real(r8), parameter :: dright = 1.02042520765343880_r8
    real(r8) :: roots(3), lo, hi, mid, expected, z, value, scale
    integer :: i, j, nroots

    ! Independent bisection of the Hermite polynomial with opposite endpoint
    ! signs and positive derivative throughout this cell gives its unique root.
    lo = 0.0_r8
    hi = 1.0_r8
    do i = 1, 60
        mid = (lo + hi)/2
        if (hermite(mid) < 0) then
            lo = mid
        else
            hi = mid
        end if
    end do
    expected = left + (right - left)*(lo + hi)/2
    do j = -3, 3, 3
        scale = 10.0_r8**j
        call spline_alloc(spl, 1, 1)
        spl%xs = [left, right]
        spl%fs(:, 1) = scale*[fleft, fright]
        spl%fs1(:, 1) = scale*[dleft, dright]
        call spline_roots(spl, 1, nroots, roots)
        print '(a,i0,a,es24.16)', 'roots=', nroots, ' oracle=', expected
        if (nroots /= 1) error stop 'Spline lost a bracketed simple root'
        if (abs(roots(1) - expected) > 5.0e-14_r8) &
            error stop 'Spline root disagrees with independent bracket oracle'
        z = (roots(1) - left)/(right - left)
        value = hermite(z)
        if (z <= 0 .or. z >= 1 .or. abs(value) > 1.0e-13_r8) &
            error stop 'Root is outside its bracket or has a large residual'
        call spline_dealloc(spl)
    end do
    print *, 'PASS: bracketed endpoint root survives amplitude scaling'
contains
    function hermite(z) result(value)
        real(r8), intent(in) :: z
        real(r8) :: value
        value = (2*z**3 - 3*z**2 + 1)*fleft + &
            (-2*z**3 + 3*z**2)*fright + &
            (z**3 - 2*z**2 + z)*(right - left)*dleft + &
            (z**3 - z**2)*(right - left)*dright
    end function hermite
end program spline_root_endpoint
