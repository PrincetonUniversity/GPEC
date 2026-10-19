program test_pmodb_target_gauge
    use local_mod, only: r8
    use pmodb_gauge_mod, only: pmodb_target_gauge
    implicit none
    integer, parameter :: nt=128, nm=25
    integer :: m(nm), i, power
    real(r8) :: th, q, pi, b(0:nt), bt(0:nt)
    complex(r8) :: u(0:nt), du(0:nt), delta(0:nt), d(0:nt), result(0:nt)
    character(16) :: mode
    call get_command_argument(1,mode)
    pi=acos(-1.0_r8)
    q=1.3_r8
    m=[(i,i=-12,12)]
    do i=0,nt
        th=2*pi*i/nt
        b(i)=3.0_r8*(1-0.2_r8*cos(th))
        bt(i)=2*pi*0.6_r8*sin(th)
        u(i)=exp(cmplx(0.0_r8,th,r8))
        du(i)=cmplx(0.0_r8,2*pi*(1-2*q),r8)*u(i)
    end do
    if (mode=='resonant') then
        d=exp(cmplx(0.0_r8,4*pi,r8)*[(real(i,r8)/nt,i=0,nt)])
        delta=0
        call pmodb_target_gauge(m,2,1.0_r8,b,bt,delta,d,0,result)
        error stop 'nonzero resonant source was accepted'
    end if
    do power=0,2,2
        ! A parallel displacement has delta=u*d_theta B and T=D u.
        delta=u*bt
        d=du-power*delta/b
        call pmodb_target_gauge(m,2,q,b,bt,delta,d,power,result)
        if (maxval(abs(result))>1.e-11_r8) error stop 'parallel displacement is not zero'
        ! With T=0 the target scalar must stay unchanged.
        delta=0.4_r8*u
        d=-power*delta/b
        call pmodb_target_gauge(m,2,q,b,bt,delta,d,power,result)
        if (maxval(abs(result-delta))>1.e-12_r8) error stop 'zero-source limit changed'
    end do
    ! Hamada scalar already relabeled to Boozer base angles; div xi_H=0.
    delta=0.2_r8*u
    d=0
    call pmodb_target_gauge(m,2,q,b,bt,delta,d,2,result)
    open(unit=31,file='pmodb_hamada_to_boozer.txt',status='replace')
    do i=0,nt-1
        write(31,'(5es24.15)') real(i,r8)/nt,b(i),real(result(i)),aimag(result(i)),bt(i)
    end do
    close(31)
    print '(a)', 'PASS: parallel displacement and zero-source limits in both gauges'
end program test_pmodb_target_gauge
