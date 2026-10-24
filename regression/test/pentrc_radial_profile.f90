! Run in a PENTRC case directory (pentrc.in, euler.bin, profiles and peq).
! The oracle evaluates the surface response directly, independently of the
! adaptive integral's dense-output polynomial.
program pentrc_radial_profile
    use pentrc_interface, only: initialize_pentrc, r8, nn, zi, mi, &
                                wdfac, divxfac, electron, nl
    use torque, only: tintgrl_lsode, trans, tpsi, atol_psi, rtol_psi, &
                      output_ascii, output_netcdf, torque_record
    implicit none
    integer :: i
    real(r8) :: bounds(2), relative_error, max_error
    complex(r8) :: integral, expected, actual

    call initialize_pentrc
    output_ascii = .false.
    output_netcdf = .false.
    atol_psi = 1.0e-4_r8
    rtol_psi = 1.0e-4_r8
    bounds = [0.02_r8, 0.6_r8]
    integral = tintgrl_lsode(bounds, nn, nl, zi, mi, wdfac, divxfac, &
                             electron, 'fgar')
    max_error = 0.0_r8
    do i = 0, trans%mx
        call tpsi(expected, trans%xs(i), nn, 0, zi, mi, wdfac, &
                  divxfac, electron, 'fgar')
        actual = cmplx(torque_record(1, 1)%ell_record(7, nl + 1, i + 1), &
                       torque_record(1, 1)%ell_record(8, nl + 1, i + 1), r8)
        relative_error = abs(actual - expected)/ &
                         max(abs(expected), tiny(1.0_r8))
        max_error = max(max_error, relative_error)
    end do
    print '(a,i0,a,es12.4)', 'Surface evaluations: ', trans%mx + 1, &
        '; max relative error: ', max_error
    if (abs(trans%xs(trans%mx) - bounds(2)) > 1.0e-14_r8) &
        error stop 'Adaptive profile does not reach the requested endpoint'
    if (max_error > 1.0e-10_r8) &
        error stop 'Adaptive profile differs from direct surface response'
    print '(a)', 'PASS: adaptive profile equals direct surface response'
end program pentrc_radial_profile
