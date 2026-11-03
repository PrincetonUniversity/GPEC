! Biot-Savart field of polygonal coils (coil/field.F).
! 1. field_bs_segment against the textbook closed form of a straight segment,
!    B = mu0 I/(4 pi rho) (cos(theta1) - cos(theta2)), from near to far points.
! 2. field_bs_rzphi for a rectangular coil passed as its 4 sides and with each
!    side split into 16 collinear pieces. The conductor is the same, so the
!    n = 1 field must agree (the former midpoint rule differs by 10 %).
! 3. With rmin > 0, points within rmin of the filament give a finite field
!    bounded by that of a wire of radius rmin.
! Build after `make` in install/ (gfortran):
!   gfortran -fopenmp -I../coil -I../equil ../regression/test/coil_segment_field.f90 \
!     -L../deps/lib -lcoil -lequil -lopenblas -lnetcdff -lnetcdf -o coil_segment_field
program coil_segment_field
    use coil_mod, only: coil, coil_num, cro, cnn, r8
    use field_mod, only: field_bs_rzphi, field_bs_segment
    implicit none
    integer, parameter :: nr = 4, nz = 4, np = 128
    real(r8), parameter :: a(3) = [0.1_r8, -0.2_r8, -0.5_r8]
    real(r8), parameter :: d(3) = [0.2_r8, 0.1_r8, 1.0_r8]
    real(r8) :: x(3), b(3), bref(3), err, rmin
    real(r8) :: gdr(0:nr, 0:nz), gdz(0:nr, 0:nz)
    complex(r8), dimension(0:nr, 0:nz) :: br1, bz1, bp1, br16, bz16, bp16
    integer :: i, ir, iz

    err = 0
    do i = 1, 8
        x = a + [0.02_r8*3**(i - 1), 0.01_r8, 0.3_r8*i - 0.9_r8]
        b = field_bs_segment(x(1) - a(1), x(2) - a(2), x(3) - a(3), &
                             d(1), d(2), d(3), 0.0_r8)*cross(d, x - a)
        bref = segment_closed_form(a, a + d, x)
        err = max(err, norm2(b - bref)/norm2(bref))
    end do
    print '(a,es10.2)', "segment vs closed form, max rel. error:", err
    if (err > 1.0e-12_r8) error stop "field_bs_segment differs from the closed form"

    do iz = 0, nz
        do ir = 0, nr
            gdr(ir, iz) = 0.95_r8 + 0.05_r8*ir
            gdz(ir, iz) = -0.3_r8 + 0.15_r8*iz
        end do
    end do
    cro = 1.0_r8
    cnn = 1
    coil_num = 1
    allocate(coil(1))
    call rectangle_field(1, br1, bz1, bp1)
    call rectangle_field(16, br16, bz16, bp16)
    err = max(maxval(abs(br1 - br16)), maxval(abs(bz1 - bz16)), &
              maxval(abs(bp1 - bp16))) &
          /max(maxval(abs(br16)), maxval(abs(bz16)), maxval(abs(bp16)))
    print '(a,es10.2)', "rzphi 4 vs 64 pieces, max rel. difference:", err
    if (err > 1.0e-12_r8) error stop "field_bs_rzphi depends on the polyline split"

    rmin = 1.0e-3_r8
    err = 0
    do i = 0, 2
        x = a + 0.3_r8*d + 0.5_r8*i*rmin*cross(d, [1.0_r8, 0.0_r8, 0.0_r8]) &
            /norm2(cross(d, [1.0_r8, 0.0_r8, 0.0_r8]))
        b = field_bs_segment(x(1) - a(1), x(2) - a(2), x(3) - a(3), &
                             d(1), d(2), d(3), rmin)*cross(d, x - a)
        if (.not. norm2(b) <= 2/rmin*(1 + 1.0e-6_r8)) &
            error stop "regularized field is not bounded by the rmin wire"
        err = max(err, norm2(b)*rmin/2)
    end do
    print '(a,f8.5)', "near-filament max |B|/(mu0 I/(2 pi rmin)):", err

contains

    pure function cross(u, v) result(w)
        real(r8), intent(in) :: u(3), v(3)
        real(r8) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]
    end function cross

    ! Field per mu0 I/(4 pi) from the angles subtended at x by the segment p--q.
    pure function segment_closed_form(p, q, x) result(bv)
        real(r8), intent(in) :: p(3), q(3), x(3)
        real(r8) :: bv(3), e(3), rho(3), r, z1, z2
        e = (q - p)/norm2(q - p)
        rho = (x - p) - dot_product(x - p, e)*e
        r = norm2(rho)
        z1 = -dot_product(x - p, e)
        z2 = z1 + norm2(q - p)
        bv = (z2/hypot(z2, r) - z1/hypot(z1, r))/r*cross(e, rho/r)
    end function segment_closed_form

    ! Rectangle 0.6 m x 0.4 m in the plane x = 1.25 m, I = 1 kA, each side in m pieces.
    subroutine rectangle_field(m, cbr, cbz, cbp)
        integer, intent(in) :: m
        complex(r8), dimension(0:nr, 0:nz), intent(out) :: cbr, cbz, cbp
        real(r8), parameter :: cy(5) = [-0.3_r8, 0.3_r8, 0.3_r8, -0.3_r8, -0.3_r8]
        real(r8), parameter :: cz(5) = [-0.2_r8, -0.2_r8, 0.2_r8, 0.2_r8, -0.2_r8]
        integer :: iside, ipiece, ipt
        real(r8) :: t

        if (associated(coil(1)%x)) &
            deallocate(coil(1)%x, coil(1)%y, coil(1)%z, coil(1)%cur)
        coil(1)%coil_name = "rectangle"
        coil(1)%ncoil = 1
        coil(1)%s = 1
        coil(1)%nsec = 4*m + 1
        coil(1)%nw = 1.0_r8
        allocate(coil(1)%x(1, 1, 4*m + 1), coil(1)%y(1, 1, 4*m + 1), &
                 coil(1)%z(1, 1, 4*m + 1), coil(1)%cur(1))
        coil(1)%cur(1) = 1.0e3_r8
        coil(1)%x = 1.25_r8
        ipt = 0
        do iside = 1, 4
            do ipiece = 0, m - 1
                t = real(ipiece, r8)/m
                ipt = ipt + 1
                coil(1)%y(1, 1, ipt) = cy(iside) + t*(cy(iside + 1) - cy(iside))
                coil(1)%z(1, 1, ipt) = cz(iside) + t*(cz(iside + 1) - cz(iside))
            end do
        end do
        coil(1)%y(1, 1, ipt + 1) = cy(5)
        coil(1)%z(1, 1, ipt + 1) = cz(5)
        cbr = 0
        cbz = 0
        cbp = 0
        call field_bs_rzphi(nr, nz, np, gdr, gdz, cbr, cbz, cbp)
    end subroutine rectangle_field

end program coil_segment_field
