! Run in a normal PENTRC case directory, as for pentrc_radial_profile.
! Replace only the external radial LSODE entry with a deliberate failure.
! The production radial wrapper must reject it before evaluating another RHS
! or publishing a partial integral. Energy/pitch solvers use other entries.
program pentrc_radial_istate
    use pentrc_interface, only: initialize_pentrc, r8, nn, zi, mi, &
        wdfac, divxfac, electron, nl
    use torque, only: tintgrl_lsode, output_ascii, output_netcdf
    implicit none
    real(r8) :: bounds(2)
    complex(r8) :: integral

    call initialize_pentrc
    output_ascii = .false.
    output_netcdf = .false.
    bounds = [0.90_r8, 0.900001_r8]
    integral = tintgrl_lsode(bounds, nn, nl, zi, mi, wdfac, divxfac, &
        electron, 'fgar')
    error stop 'Radial wrapper returned after a failed LSODE call'
end program pentrc_radial_istate

subroutine lsode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
    iopt, rwork, lrw, iwork, liw, jac, mf)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    external :: f, jac
    integer, intent(in) :: neq(*), itol, itask, iopt, lrw, liw, mf
    integer, intent(inout) :: istate, iwork(liw)
    real(dp), intent(in) :: tout, rtol(*), atol(*)
    real(dp), intent(inout) :: y(*), t, rwork(lrw)
    integer, save :: calls = 0
    character(16) :: argument
    integer :: requested_status

    calls = calls + 1
    if (calls > 1) error stop 'Radial wrapper retried failed LSODE'
    call get_command_argument(1, argument)
    read(argument, *) requested_status
    if (requested_status >= 0) error stop 'Test needs a negative LSODE status'
    print '(a,i0)', 'Injected radial LSODE status ', requested_status
    istate = requested_status
end subroutine lsode
