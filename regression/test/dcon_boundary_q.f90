! After make -C install dcon, link only the DCON boundary solver on Linux.
! gfortran -fopenmp -ffunction-sections -fdata-sections -Idcon -Iequil -Ipentrc -Ivacuum -Jdcon &
!    -c dcon/sing.f -o /path/to/test_sing.o
! gfortran -fopenmp -Wl,--gc-sections -Idcon -Iequil -Ipentrc -Ivacuum &
!    regression/test/dcon_boundary_q.f90 /path/to/test_sing.o dcon/dcon_mod.o &
!    -Ldeps/lib -lequil -llapack -lblas -lnetcdff -lnetcdf -o dcon_boundary_q
program dcon_boundary_q
    use sing_mod, only: nn, mpsi, sq, spline_alloc, spline_fit_hermite, sing_lim, &
                        psihigh, psiedge, psilim, sas_flag, qmax, qhigh, qlim, q1lim
    use params, only: r8
    implicit none
    real(r8), parameter :: tol = 1.0e-10_r8
    integer :: i

    mpsi = 8
    nn = 1
    sas_flag = .false.
    psihigh = 1.0_r8
    psiedge = psihigh
    call spline_alloc(sq, mpsi, 4)
    sq%xs = [(real(i, r8) / mpsi, i = 0, mpsi)]
    sq%fs = 0.0_r8
    sq%fs1 = 0.0_r8

    sq%fs(:, 4) = 2.6_r8 + 2.0_r8 * (sq%xs - 0.6_r8)**2
    sq%fs1(:, 4) = 4.0_r8 * (sq%xs - 0.6_r8)
    qmax = 3.32_r8
    call spline_fit_hermite(sq, "extrap", [.false., .false., .false., .true.])
    qhigh = 1000.0_r8
    call check_limit(1.0_r8, 2.92_r8, 1.6_r8)
    qhigh = 3.1_r8
    call check_limit(1.0_r8, 2.92_r8, 1.6_r8)

    sq%fs(:, 4) = 2.0_r8 + 2.0_r8 * sq%xs
    sq%fs1(:, 4) = 2.0_r8
    qmax = 4.0_r8
    call spline_fit_hermite(sq, "extrap", [.false., .false., .false., .true.])
    qhigh = 3.75_r8
    call check_limit(0.875_r8, 3.75_r8, 2.0_r8)

contains

    subroutine check_limit(expected_psi, expected_q, expected_qprime)
        real(r8), intent(in) :: expected_psi, expected_q, expected_qprime

        call sing_lim
        print '(a,3f12.8)', "psi, q, qprime at the boundary ", psilim, qlim, q1lim
        if (abs(psilim - expected_psi) > tol) error stop "wrong boundary flux"
        if (abs(qlim - expected_q) > tol) error stop "wrong boundary safety factor"
        if (abs(q1lim - expected_qprime) > tol) error stop "wrong boundary shear"
    end subroutine check_limit

end program dcon_boundary_q
