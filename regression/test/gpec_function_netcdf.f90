! Run: make -C install test-gpec-function-netcdf
! Exercise the production pmodb packer and the native NetCDF interface.
program gpec_function_netcdf
    use gpout_mod, only: gpout_function_channels, r8
    use netcdf
    implicit none
    integer, parameter :: npsi = 3, ntheta = 5
    complex(r8) :: values(npsi, ntheta)
    real(r8) :: restored(npsi, ntheta, 2), expected_real, expected_imag
    real(r8) :: imag_sign
    integer :: i, j, sign_case, field, ncid, pdim, tdim, cdim, varid
    integer :: shape_read(3), dimension_ids(3), k

    do j = 1, ntheta - 1
        do i = 1, npsi
            values(i, j) = cmplx(10*i + j, -100*i - 3*j, r8)
        end do
    end do
    values(:, ntheta) = values(:, 1)

    ! Both helicities, and both signs of the four pmodb quantities.
    do sign_case = -1, 1, 2
        imag_sign = real(sign_case, r8)
        do field = -1, 1, 2
            call check(nf90_create('gpec_function_roundtrip.nc', nf90_clobber, ncid))
            call check(nf90_def_dim(ncid, 'psi', npsi, pdim))
            call check(nf90_def_dim(ncid, 'theta', ntheta, tdim))
            call check(nf90_def_dim(ncid, 'complex', 2, cdim))
            call check(nf90_def_var(ncid, 'function', nf90_double, &
                [pdim, tdim, cdim], varid))
            call check(nf90_enddef(ncid))
            call check(nf90_put_var(ncid, varid, &
                gpout_function_channels(field*values, imag_sign)))
            call check(nf90_close(ncid))

            call check(nf90_open('gpec_function_roundtrip.nc', nf90_nowrite, ncid))
            call check(nf90_inq_varid(ncid, 'function', varid))
            call check(nf90_inquire_variable(ncid, varid, dimids=dimension_ids))
            do k = 1, 3
                call check(nf90_inquire_dimension(ncid, dimension_ids(k), len=shape_read(k)))
            end do
            if (any(shape_read /= [npsi, ntheta, 2])) error stop 'Wrong NetCDF shape'
            call check(nf90_get_var(ncid, varid, restored))
            call check(nf90_close(ncid))
            do j = 1, ntheta
                do i = 1, npsi
                    k = j
                    if (j == ntheta) k = 1
                    ! Independent scalar oracle, including the periodic endpoint.
                    expected_real = real(field*(10*i + k), r8)
                    expected_imag = imag_sign*real(field*(-100*i - 3*k), r8)
                    if (restored(i, j, 1) /= expected_real) error stop 'Real channel corrupted'
                    if (restored(i, j, 2) /= expected_imag) error stop 'Imaginary channel corrupted'
                end do
            end do
        end do
    end do
    open(unit=99, file='gpec_function_roundtrip.nc', status='old')
    close(99, status='delete')
    print *, 'PASS: native pmodb NetCDF channels and periodic endpoint'
contains
    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            error stop 'Native NetCDF call failed'
        end if
    end subroutine check
end program gpec_function_netcdf
