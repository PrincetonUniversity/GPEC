! Energy integral with an LSODE failure other than ISTATE = -1.
! "default": default tolerances, must return a finite integral.
! "tight": rtol = 1e-17, atol = 1e-30 makes LSODE return ISTATE = -2 (too much
! accuracy requested) near x = 8e-4; xintgrl_lsode must not return.
program pentrc_energy_istate
    use params, only: r8
    use energy_integration, only: xintgrl_lsode, xatol, xrtol
    implicit none
    character(16) :: mode
    complex(r8) :: integral

    call get_command_argument(1, mode)
    if (trim(mode) == "tight") then
        xatol = 1.0e-30_r8
        xrtol = 1.0e-17_r8
    end if
    integral = xintgrl_lsode(0.1_r8, 1.0_r8, 0.5_r8, 0.1_r8, 0.0_r8, 1.0e-2_r8, &
                             0, 0.0_r8, 1, 0.5_r8, 0.5_r8, "fgar")
    print '(a,2es24.15)', "integral ", integral
    if (trim(mode) == "tight") error stop "xintgrl_lsode returned after an LSODE failure"
    if (.not. (abs(integral) > 0.0_r8 .and. abs(integral) < huge(1.0_r8))) &
        error stop "default integral is not finite and nonzero"
end program pentrc_energy_istate
