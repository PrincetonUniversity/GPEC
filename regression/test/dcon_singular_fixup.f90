! Link the production DCON objects, replacing only the main program.
! Both scenarios require a native singular reduction record after new=.true.
! Run: make -C install test-dcon-singular-fixup (normal compiler/dependency flags).
program dcon_singular_fixup
    use ode_mod, only: r8, ode_unorm, new, u, unorm, unorm0, fixfac, index, &
        mpert, msol, mlow, mhigh, istep, psifac, q, ucrit, bin_euler, &
        crit_out_unit, euler_bin_unit
    implicit none

    complex(r8) :: before(2, 2, 2), saved_fix(2, 2), reconstructed(2, 2)
    integer :: scenario, expected_count, count, record_type, saved_msol, status, i
    integer :: saved_index(2)
    logical :: saved_singular
    character(16) :: argument

    call get_command_argument(1, argument)
    read(argument, *, iostat=status) scenario
    if (status /= 0) error stop 'Specify scenario 1 or 2'
    if (scenario < 1 .or. scenario > 2) error stop 'Unknown scenario'
    mpert = 2
    msol = 2
    mlow = 1
    mhigh = 2
    istep = 1
    psifac = 0.5_r8
    q = 3.0_r8
    ucrit = 10.0_r8
    bin_euler = .true.
    allocate(u(2, 2, 2), unorm(2), unorm0(2), fixfac(2, 2), index(2))
    u = cmplx(0.0_r8, 0.0_r8, r8)
    u(:, 1, 1) = [1.0_r8, 0.0_r8]
    u(:, 2, 1) = [1.0_r8, 100.0_r8]
    unorm0 = 1.0_r8
    open(unit=crit_out_unit, status='scratch', form='formatted')
    open(unit=euler_bin_unit, status='scratch', form='unformatted')

    expected_count = 1
    new = .true.
    if (scenario == 2) then
        new = .false.
        call ode_unorm(.false.)
        if (.not. new) error stop 'Manufactured ordinary reduction did not occur'
        expected_count = 2
    end if
    before = u
    call ode_unorm(.true.)
    rewind(euler_bin_unit)
    count = 0
    do
        read(euler_bin_unit, iostat=status) record_type
        if (status < 0) exit
        if (status /= 0) error stop 'Unreadable native reduction record'
        if (record_type /= 2) error stop 'Unexpected native record type'
        read(euler_bin_unit) saved_singular, saved_msol
        if (saved_msol /= 2) error stop 'Wrong reduction solution count'
        read(euler_bin_unit) saved_fix, saved_index
        count = count + 1
        if (count == expected_count) then
            if (.not. saved_singular) error stop 'Singular reduction record missing'
            if (any(saved_index < 1) .or. any(saved_index > 2)) &
                error stop 'Invalid reduction permutation'
            if (saved_index(1) == saved_index(2)) error stop 'Repeated reduction index'
            do i = 1, 2
                reconstructed = matmul(before(:, :, i), saved_fix)
                if (maxval(abs(reconstructed - u(:, :, i))) > 1.0e-13_r8) &
                    error stop 'Recorded transform does not reproduce native basis'
            end do
        else
            if (saved_singular) error stop 'Ordinary reduction incorrectly marked singular'
        end if
    end do
    if (count /= expected_count) error stop 'Forced singular reduction record missing'
    if (.not. new) error stop 'Forced reduction did not request norm reset'
    close(euler_bin_unit)
    close(crit_out_unit)
    print *, 'PASS: native forced singular reduction and serialized transform', scenario
end program dcon_singular_fixup
