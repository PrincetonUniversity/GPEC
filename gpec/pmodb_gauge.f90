module pmodb_gauge_mod
    use local_mod, only: r8
    use ismath_mod, only: iscdftf, iscdftb
    use pmodb_gauge_kernel, only: pmodb_gauge_inverse, pmodb_gauge_source, pmodb_gauge_apply
    implicit none
contains
    subroutine pmodb_target_gauge(m, n, q, b, bt, delta, divxi, power, target)
        integer, intent(in) :: m(:), n, power
        real(r8), intent(in) :: q, b(0:), bt(0:)
        complex(r8), intent(in) :: delta(0:), divxi(0:)
        complex(r8), intent(out) :: target(0:)
        complex(r8) :: source(size(m)), w(0:size(b)-1), t(0:size(b)-1)
        real(r8) :: s, out_r, out_i
        integer :: j, nt
        nt = size(b)-1
        do j = 0, nt
            call pmodb_gauge_source(real(power,r8), b(j), real(delta(j),r8), aimag(delta(j)), &
                real(divxi(j),r8), aimag(divxi(j)), out_r, out_i)
            t(j) = cmplx(out_r, out_i, r8)
        end do
        call iscdftf(m, size(m), t, nt, source)
        do j = 1, size(m)
            s = m(j)-n*q
            if (abs(s) < 1.0e-10_r8) then
                if (abs(source(j)) > 1.0e-12_r8*max(maxval(abs(source)), tiny(1.0_r8))) &
                    error stop 'Nonzero resonant pmodb gauge source'
                source(j) = 0
            else
                call pmodb_gauge_inverse(s, real(source(j),r8), aimag(source(j)), out_r, out_i)
                source(j) = cmplx(out_r, out_i, r8)
            end if
        end do
        call iscdftb(m, size(m), w, nt, source)
        do j = 0, nt
            call pmodb_gauge_apply(real(delta(j),r8), aimag(delta(j)), real(w(j),r8), aimag(w(j)), bt(j), out_r, out_i)
            target(j) = cmplx(out_r, out_i, r8)
        end do
    end subroutine pmodb_target_gauge
end module pmodb_gauge_mod
