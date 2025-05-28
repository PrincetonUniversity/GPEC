program compare_2
    ! use crmath

    integer, parameter :: r8 = selected_real_kind(13, 307)

    real(r8) :: eta = 2.1127504979480101_r8
    real(r8) :: coseta
    real(r8) :: sineta

    write(*,*) "Testing compiler functions with real kind r8"
    write(*,*) "eta = ", eta
    coseta = cos(eta)
    sineta = sin(eta)
    write(*,*) "cos(eta) = ", coseta
    write(*,*) "sin(eta) = ", sineta

end program compare_2