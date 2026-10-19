module pmodb_gauge_mod
    use local_mod, only: r8
    use ismath_mod, only: iscdftf, iscdftb
    use pmodb_gauge_kernel, only: pmodb_gauge_inverse
    implicit none
contains
    subroutine pmodb_target_gauge(m, n, q, b, bt, delta, divxi, power, target)
        integer, intent(in) :: m(:), n, power
        real(r8), intent(in) :: q, b(0:), bt(0:)
        complex(r8), intent(in) :: delta(0:), divxi(0:)
        complex(r8), intent(out) :: target(0:)
        complex(r8) :: source(size(m)), w(0:size(b)-1)
        real(r8) :: s, weight
        integer :: j, nt
        nt = size(b)-1
        call iscdftf(m, size(m), divxi+power*delta/b, nt, source)
        do j = 1, size(m)
            s = m(j)-n*q
            if (abs(s) < 1.0e-10_r8) then
                if (abs(source(j)) > 1.0e-12_r8*max(maxval(abs(source)), tiny(1.0_r8))) &
                    error stop 'Nonzero resonant pmodb gauge source'
                source(j) = 0
            else
                call pmodb_gauge_inverse(s, weight)
                source(j) = cmplx(0.0_r8, 1.0_r8, r8)*weight*source(j)
            end if
        end do
        call iscdftb(m, size(m), w, nt, source)
        target = delta+w*bt
    end subroutine pmodb_target_gauge
end module pmodb_gauge_mod
