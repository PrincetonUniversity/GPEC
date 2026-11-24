! Exercise the native pitch integral with recording on and off at two surfaces.
program pentrc_record_first_surface
    use params, only: r8
    use netcdf
    use utilities, only: check
    use cspline_mod, only: cspline_type, cspline_alloc, cspline_fit
    use spline_mod, only: spline_type, spline_alloc, spline_fit
    use pitch_integration, only: lambdaintgrl_lsode, output_pitch_netcdf
    use energy_integration, only: output_energy_netcdf
    implicit none

    type(cspline_type) :: bounce
    type(spline_type) :: turns
    complex(r8), allocatable :: recorded(:), ordinary(:)
    integer :: i, surface, ell
    real(r8) :: psi

    call cspline_alloc(bounce,4,3)
    call spline_alloc(turns,4,6)
    do i=0,4
        bounce%xs(i) = 0.875_r8 + 0.03125_r8*i
        bounce%fs(i,:) = cmplx([2.0_r8,0.25_r8,1.0_r8 + 0.1_r8*i],0.0_r8,r8)
        turns%xs(i) = bounce%xs(i)
        turns%fs(i,:) = 0.0_r8
    end do
    call cspline_fit(bounce,"extrap")
    call spline_fit(turns,"extrap")
    do surface=1,2
        psi = 0.25_r8*surface
        do ell=-1,1
            ordinary = integral(.false.)
            recorded = integral(.true.)
            if (abs(ordinary(1)-recorded(1)) > 1.0e-6_r8*abs(ordinary(1))) &
                error stop "recording changes the integral"
        end do
    end do
    call output_pitch_netcdf(1)
    call output_energy_netcdf(1)
    call check_labels("pentrc_pitch_output_n1.nc", energy=.false.)
    call check_labels("pentrc_energy_output_n1.nc", energy=.true.)
contains
    subroutine check_labels(path, energy)
        character(*), intent(in) :: path
        logical, intent(in) :: energy
        integer :: ncid, varid, harmonics(3), j
        real(r8) :: surfaces(2), pitches(5,3,2)

        call check(nf90_open(path, nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, "psi_n", varid))
        call check(nf90_get_var(ncid, varid, surfaces))
        if (any(surfaces /= [0.25_r8, 0.5_r8])) error stop "Incorrect surface labels"
        call check(nf90_inq_varid(ncid, "ell", varid))
        call check(nf90_get_var(ncid, varid, harmonics))
        if (any(harmonics /= [-1, 0, 1])) error stop "Incorrect harmonic labels"
        if (energy) then
            call check(nf90_inq_varid(ncid, "Lambda_tgar", varid))
            call check(nf90_get_var(ncid, varid, pitches))
            do j=1,5
                if (any(abs(pitches(j,:,:) - (0.875_r8+0.03125_r8*(j-1))) > 1.0e-14_r8)) &
                    error stop "Incorrect pitch labels"
            end do
        end if
        call check(nf90_close(ncid))
    end subroutine check_labels

    function integral(record_this) result(response)
        logical, intent(in) :: record_this
        complex(r8), allocatable :: response(:)
        response = lambdaintgrl_lsode(0.2_r8,0.3_r8,0.4_r8,0.1_r8,0.8_r8,0.2_r8, &
                                     2.5_r8,bounce,ell,1,1.0_r8,1.0_r8,psi,turns,"tgar",record_this)
    end function integral
end program pentrc_record_first_surface
