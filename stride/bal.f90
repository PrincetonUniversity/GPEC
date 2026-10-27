module stride_bal_mod
  use stride_dcon_mod
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, private :: ipsi
  real(r8), private :: theta_max
  type(spline_type), private :: bf

contains

  subroutine bal_scan
    real(r8) :: di, h, psi, delta_prime
    real(r8), dimension(0:mpsi) :: volume
    type(spline_type) :: volume_spline

    volume = 0.0_r8
    if (bal_flag .and. bal_alpha_flag) then
      call spline_alloc(volume_spline,mpsi,1)
      volume_spline%xs = sq%xs
      volume_spline%fs(:,1) = sq%fs(:,3)
      volume_spline%fs1(:,1) = sq%fs1(:,3)
      call spline_int(volume_spline)
      volume = max(volume_spline%fsi(:,1),0.0_r8)
      call spline_dealloc(volume_spline)
    end if

    do ipsi = 0, mpsi
      psi = sq%xs(ipsi)
      call bal_prep(di,h,0.0_r8)
      locstab%fs(ipsi,1) = di * psi
      locstab%fs(ipsi,2) = (di + (h - 0.5_r8)**2) * psi
      locstab%fs(ipsi,3) = h
      locstab%fs(ipsi,4) = ieee_value(0.0_r8,ieee_quiet_nan)
      locstab%fs(ipsi,5) = ieee_value(0.0_r8,ieee_quiet_nan)
      if (bal_flag .and. psi <= 1.0_r8) then
        call bal_int(delta_prime)
        locstab%fs(ipsi,4) = delta_prime
      end if
      call spline_dealloc(bf)
      if (bal_flag .and. bal_alpha_flag) &
        call bal_alpha_first_boundary(volume(ipsi))
    end do
  end subroutine bal_scan

  subroutine bal_prep(di,h,corr_pprime)
    real(r8), intent(out) :: di, h
    real(r8), intent(in) :: corr_pprime

    integer :: itheta
    real(r8) :: psi, twopif, p1, v1, q, q1, chi1
    real(r8) :: theta, rho, eta, ceta, seta, r, jac
    real(r8) :: drfac_dp, eta_dp, dr_dp
    real(r8) :: dv21_dp, dv22_dp, dv23_dp, dv33_dp
    real(r8) :: iper_ref, local_scale, cscale
    real(r8) :: abar_w, bbar_w, response_denom
    real(r8) :: dbase(3), dB(3), bvec(3)
    real(r8) :: gradpsi(3), gradtheta(3), gradzeta(3)
    real(r8) :: v(3,3), w(3,3)
    real(r8) :: d0bar(2,2)
    real(r8), dimension(0:mtheta) :: bsq, gradpsi_sq, dpsisq, jacarr
    real(r8), dimension(0:mtheta) :: bdot_theta_phi, iper, htheta
    real(r8), dimension(0:mtheta) :: dbsq_dtheta, kappas
    real(r8), dimension(0:mtheta) :: jac_psi, bsq_psi
    real(r8), dimension(0:mtheta) :: spl1_vals, spl1_deriv
    real(r8), dimension(0:mtheta) :: a0, a1, a2
    real(r8), dimension(0:mtheta) :: kappaw0, kappaw1
    real(r8), dimension(0:mtheta) :: aw, bw, aperp
    real(r8), dimension(0:mtheta) :: diper_dtheta, iper_perturbed
    real(r8), dimension(0:mtheta) :: jacchi, pchi, sigma
    real(r8), dimension(0:mtheta) :: m012, m021
    real(r8), dimension(0:mtheta,3) :: hfs
    real(r8), dimension(0:mtheta,4) :: n0
    real(r8), dimension(0:mtheta,7) :: coeff
    real(r8), dimension(0:mtheta,4) :: f0, fx0, fy0, fxy0
    real(r8), dimension(0:mtheta) :: rho0, eta0, r0
    real(r8), dimension(0:mtheta,3,3) :: v0
    real(r8), dimension(0:mtheta,3) :: vzp0, vpt0
    type(spline_type) :: work

    psi = sq%xs(ipsi)
    twopif = sq%fs(ipsi,1)
    p1 = sq%fs1(ipsi,2)
    v1 = sq%fs(ipsi,3)
    q = sq%fs(ipsi,4)
    q1 = sq%fs1(ipsi,4)
    chi1 = twopi * psio

    do itheta = 0, mtheta
      theta = rzphi%ys(itheta)
      call bicube_eval(rzphi, rzphi%xs(ipsi), theta, 2)
      f0(itheta,:) = rzphi%f
      fx0(itheta,:) = rzphi%fx
      fy0(itheta,:) = rzphi%fy
      fxy0(itheta,:) = rzphi%fxy

      rho = sqrt(f0(itheta,1))
      eta = twopi * (theta + f0(itheta,2))
      ceta = cos(eta)
      seta = sin(eta)
      r = ro + rho * ceta
      jac = f0(itheta,4)
      rho0(itheta) = rho
      eta0(itheta) = eta
      r0(itheta) = r
      jacarr(itheta) = jac

      v = 0.0_r8
      v(1,1) = fx0(itheta,1) / (2.0_r8*rho*jac)
      v(1,2) = fx0(itheta,2) * twopi*rho / jac
      v(1,3) = fx0(itheta,3) * r / jac
      v(2,1) = fy0(itheta,1) / (2.0_r8*rho*jac)
      v(2,2) = (1.0_r8 + fy0(itheta,2)) * twopi*rho / jac
      v(2,3) = fy0(itheta,3) * r / jac
      v(3,3) = twopi*r / jac

      w = 0.0_r8
      w(1,1) = (1.0_r8 + fy0(itheta,2))*twopi**2*rho*r/jac
      w(1,2) = -fy0(itheta,1)*pi*r/(rho*jac)
      w(2,1) = -fx0(itheta,2)*twopi**2*r*rho/jac
      w(2,2) = fx0(itheta,1)*pi*r/(rho*jac)
      w(3,1) = (fx0(itheta,2)*fy0(itheta,3) - &
        fx0(itheta,3)*(1.0_r8 + fy0(itheta,2)))*twopi*r*rho/jac
      w(3,2) = (fx0(itheta,3)*fy0(itheta,1) - &
        fx0(itheta,1)*fy0(itheta,3))*r/(2.0_r8*rho*jac)
      w(3,3) = 1.0_r8/(twopi*r)

      gradpsi = w(1,:)
      gradtheta = w(2,:)
      gradzeta = w(3,:)
      bvec = (v(2,:) + q*v(3,:))*chi1
      bsq(itheta) = dot_product(bvec,bvec)
      gradpsi_sq(itheta) = dot_product(gradpsi,gradpsi)
      dpsisq(itheta) = (twopi*r)**2 * &
        (v(2,1)**2 + v(2,2)**2)
      iper(itheta) = dot_product(q*gradtheta-gradzeta,gradpsi) / &
        gradpsi_sq(itheta)
      htheta(itheta) = fy0(itheta,3)/twopi
      bdot_theta_phi(itheta) = dot_product(bvec,v(1,:))
      v0(itheta,:,:) = v
      vzp0(itheta,:) = v(2,:)
      vpt0(itheta,:) = v(3,:)
    end do

    call periodic_lowpass_derivative(bsq, rzphi%ys, 16, dbsq_dtheta)
    kappas = dbsq_dtheta*twopif/(2.0_r8*jacarr*bsq**2)

    do itheta = 0, mtheta
      jac_psi(itheta) = fx0(itheta,4)
      drfac_dp = fx0(itheta,1)/(2.0_r8*rho0(itheta))
      eta_dp = twopi*fx0(itheta,2)
      dr_dp = drfac_dp*cos(eta0(itheta)) - &
        rho0(itheta)*sin(eta0(itheta))*eta_dp

      dv21_dp = (fxy0(itheta,1)/(2.0_r8*rho0(itheta)) - &
        fx0(itheta,1)*fy0(itheta,1)/(4.0_r8*f0(itheta,1)* &
        rho0(itheta)))/jacarr(itheta) - &
        v0(itheta,2,1)*jac_psi(itheta)/jacarr(itheta)
      dv22_dp = (fxy0(itheta,2)*twopi*rho0(itheta) + &
        (1.0_r8+fy0(itheta,2))*twopi*drfac_dp)/jacarr(itheta) - &
        v0(itheta,2,2)*jac_psi(itheta)/jacarr(itheta)
      dv23_dp = (fxy0(itheta,3)*r0(itheta) + &
        fy0(itheta,3)*dr_dp)/jacarr(itheta) - &
        v0(itheta,2,3)*jac_psi(itheta)/jacarr(itheta)
      dv33_dp = twopi*dr_dp/jacarr(itheta) - &
        v0(itheta,3,3)*jac_psi(itheta)/jacarr(itheta)

      dbase = (/dv21_dp,dv22_dp,dv23_dp/) + &
        q*(/0.0_r8,0.0_r8,dv33_dp/) + q1*v0(itheta,3,:)
      dB = dbase*chi1
      bsq_psi(itheta) = 2.0_r8*dot_product( &
        vzp0(itheta,:)+q*vpt0(itheta,:),dB)*chi1
    end do

    spl1_vals = jacarr*bdot_theta_phi/sqrt(bsq)
    call periodic_lowpass_derivative(spl1_vals, rzphi%ys, 32, &
      spl1_deriv)
    kappaw0 = -(jac_psi/jacarr + 0.5_r8*bsq_psi/bsq)/chi1 + &
      spl1_deriv/(jacarr*sqrt(bsq)) + twopif*q1/(bsq*jacarr)

    ! corr_pprime also drives the I_per response.  This is the
    ! corr_qprime=0 specialization of Julia's response calculation.
    aw = twopi**2*r0**2/(chi1**2*gradpsi_sq)
    bw = twopif**2/(chi1**2*gradpsi_sq)
    hfs(:,1) = (q+htheta)*aw
    hfs(:,2) = (q+htheta)*bw
    hfs(mtheta,1:2) = hfs(0,1:2)
    call spline_alloc(work,mtheta,2)
    work%xs = rzphi%ys
    work%fs = hfs(:,1:2)
    call spline_fit(work,"periodic")
    call spline_int(work)
    abar_w = work%fsi(mtheta,1)/q
    bbar_w = work%fsi(mtheta,2)/q
    call spline_dealloc(work)
    response_denom = 1.0_r8+bbar_w
    aperp = (q+htheta)*((aw-abar_w) - &
      abar_w/response_denom*(bw-bbar_w))
    diper_dtheta = aperp*corr_pprime
    diper_dtheta(mtheta) = diper_dtheta(0)
    call spline_alloc(work,mtheta,1)
    work%xs = rzphi%ys
    work%fs(:,1) = diper_dtheta
    call spline_fit(work,"periodic")
    call spline_int(work)
    iper_perturbed = work%fsi(:,1)
    call spline_dealloc(work)

    iper_ref = iper(0)
    iper = iper - iper_ref + iper_perturbed
    p1 = p1+corr_pprime
    kappaw0 = kappaw0-kappas*iper_perturbed
    a0 = 1.0_r8/(chi1**2*gradpsi_sq) + &
      gradpsi_sq/bsq*iper**2
    a1 = gradpsi_sq/bsq*(2.0_r8*q1*iper)
    a2 = gradpsi_sq/bsq*q1**2
    jacchi = jacarr/chi1
    pchi = p1/chi1
    kappaw1 = q1*kappas

    coeff(:,1) = a0
    coeff(:,2) = a1
    coeff(:,3) = a2
    coeff(:,4) = kappaw0
    coeff(:,5) = kappaw1
    coeff(:,6) = jacchi
    coeff(:,7) = pchi

    theta_max = 0.0_r8
    do itheta = 0, mtheta
      if (abs(a2(itheta)) > epsilon(1.0_r8)) then
        cscale = a0(itheta)-a1(itheta)**2/(4.0_r8*a2(itheta))
        local_scale = sqrt(abs(cscale/a2(itheta)))
        if (ieee_is_finite(local_scale)) &
          theta_max = max(theta_max,local_scale)
      end if
    end do
    if (theta_max <= 0.0_r8) theta_max = 16.5_r8
    theta_max = min(16.5_r8,10.0_r8*theta_max)
    theta_max = max(theta_max,1.0e-2_r8)

    sigma = -(twopif*p1*q1)/(chi1**2*bsq)
    m012 = jacchi/a2
    m021 = -2.0_r8*jacchi*pchi*kappaw0
    n0(:,1) = 0.5_r8 + m012*sigma
    n0(:,2) = m012
    n0(:,3) = m021 - sigma - m012*sigma**2
    n0(:,4) = -0.5_r8 - m012*sigma
    n0(mtheta,:) = n0(0,:)

    call spline_alloc(work,mtheta,4)
    work%xs = rzphi%ys
    work%fs = n0
    call spline_fit(work,"periodic")
    call spline_int(work)
    d0bar(1,1) = work%fsi(mtheta,1)
    d0bar(1,2) = work%fsi(mtheta,2)
    d0bar(2,1) = work%fsi(mtheta,3)
    d0bar(2,2) = work%fsi(mtheta,4)
    call spline_dealloc(work)
    di = d0bar(1,1)*d0bar(2,2)-d0bar(1,2)*d0bar(2,1)

    hfs(:,1) = bsq/dpsisq*jacarr/v1
    hfs(:,2) = 1.0_r8/dpsisq*jacarr/v1
    hfs(:,3) = bsq*jacarr/v1
    hfs(mtheta,:) = hfs(0,:)
    call spline_alloc(work,mtheta,3)
    work%xs = rzphi%ys
    work%fs = hfs
    call spline_fit(work,"periodic")
    call spline_int(work)
    h = twopif*p1*v1/(q1*chi1**3)*(work%fsi(mtheta,2) - &
      work%fsi(mtheta,1)/work%fsi(mtheta,3))
    call spline_dealloc(work)

    coeff(mtheta,:) = coeff(0,:)
    call spline_alloc(bf,mtheta,7)
    bf%xs = rzphi%ys
    bf%fs = coeff
    call spline_fit(bf,"periodic")
  end subroutine bal_prep

  subroutine periodic_lowpass_derivative(vals,xs,nkeep,deriv)
    real(r8), intent(in) :: vals(0:), xs(0:)
    integer, intent(in) :: nkeep
    real(r8), intent(out) :: deriv(0:)
    integer :: n, kmax, j, l, k
    real(r8) :: period, omega, phase, ak, bk

    n = ubound(vals,1)
    period = xs(n)-xs(0)
    kmax = min(nkeep,(n-1)/2)
    deriv = 0.0_r8
    do k = 1, kmax
      ak = 0.0_r8
      bk = 0.0_r8
      do l = 0, n-1
        phase = twopi*real(k*l,r8)/real(n,r8)
        ak = ak + vals(l)*cos(phase)
        bk = bk + vals(l)*sin(phase)
      end do
      ak = 2.0_r8*ak/real(n,r8)
      bk = 2.0_r8*bk/real(n,r8)
      omega = twopi*real(k,r8)/period
      do j = 0, n-1
        phase = twopi*real(k*j,r8)/real(n,r8)
        deriv(j) = deriv(j) + omega*(-ak*sin(phase)+bk*cos(phase))
      end do
    end do
    deriv(n) = deriv(0)
  end subroutine periodic_lowpass_derivative

  subroutine bal_der(neq,theta,y,dy)
    integer, intent(in) :: neq
    real(r8), intent(in) :: theta
    real(r8), intent(in) :: y(neq)
    real(r8), intent(out) :: dy(neq)
    real(r8) :: f, kappaw, g, theta_wrapped

    theta_wrapped = modulo(theta-bf%xs(0), &
      bf%xs(bf%mx)-bf%xs(0))+bf%xs(0)
    call spline_eval(bf,theta_wrapped,0)
    f = (bf%f(1)+theta*bf%f(2)+theta**2*bf%f(3))/bf%f(6)
    kappaw = bf%f(4)-theta*bf%f(5)
    g = -2.0_r8*bf%f(6)*bf%f(7)*kappaw
    dy(1) = y(2)/f
    dy(2) = g*y(1)
  end subroutine bal_der

  subroutine bal_integrate_side(theta_start,theta_end,y,success)
    integer, parameter :: neq=2, liw=20, lrw=22+neq*16
    real(r8), intent(in) :: theta_start, theta_end
    real(r8), intent(inout) :: y(neq)
    logical, intent(out) :: success
    integer :: istate, itask, itol, iopt, jac, mf
    integer :: iwork(liw)
    real(r8) :: theta, rtol, atol, rwork(lrw)

    iwork = 0
    rwork = 0.0_r8
    theta = theta_start
    itol = 1
    rtol = 1.0e-8_r8
    atol = 1.0e-16_r8
    itask = 1
    istate = 1
    iopt = 1
    iwork(6) = 100000
    rwork(7) = 1.0e-10_r8
    jac = 0
    mf = 10
    call lsode(bal_der,neq,y,theta,theta_end,itol,rtol,atol,itask, &
      istate,iopt,rwork,lrw,iwork,liw,jac,mf)
    success = istate == 2
  end subroutine bal_integrate_side

  subroutine bal_int(delta_prime)
    real(r8), parameter :: match_point=1.0e-3_r8
    real(r8), intent(out) :: delta_prime
    real(r8) :: yl(2), yr(2)
    logical :: ok_left, ok_right

    yl = (/0.0_r8,1.0_r8/)
    yr = (/0.0_r8,-1.0_r8/)
    call bal_integrate_side(-theta_max,-match_point,yl,ok_left)
    call bal_integrate_side(theta_max,match_point,yr,ok_right)
    if (ok_left .and. ok_right .and. &
        abs(yl(1)) > tiny(1.0_r8) .and. &
        abs(yr(1)) > tiny(1.0_r8)) then
      delta_prime = yr(2)/yr(1)-yl(2)/yl(1)
    else
      delta_prime = ieee_value(0.0_r8,ieee_quiet_nan)
    end if
  end subroutine bal_int

  subroutine bal_delta_at_scale(scale,delta_prime)
    real(r8), intent(in) :: scale
    real(r8), intent(out) :: delta_prime
    real(r8) :: corr_pprime, di_dummy, h_dummy

    corr_pprime = sq%fs1(ipsi,2)*(scale-1.0_r8)
    call bal_prep(di_dummy,h_dummy,corr_pprime)
    call bal_int(delta_prime)
    call spline_dealloc(bf)
  end subroutine bal_delta_at_scale

  subroutine bal_alpha_first_boundary(volume)
    real(r8), intent(in) :: volume
    integer :: iscan
    real(r8), parameter :: pole_cap=3.0_r8
    real(r8) :: alpha_ref, scale, scale_prev, step
    real(r8) :: d, d_prev, d_anchor, lo, hi, mid, dm, sign_lo
    logical :: have_anchor

    if (sq%xs(ipsi) < 0.1_r8 .or. &
        sq%xs(ipsi) > min(0.99_r8,sq%xs(mpsi))) return
    if (bal_alpha_nscan < 1 .or. bal_alpha_max_scale <= 0.0_r8) return

    alpha_ref = -2.0_r8*sq%fs1(ipsi,2)*sq%fs(ipsi,3) / &
      psio**2*sqrt(volume/(twopi*max(ro,epsilon(1.0_r8))))
    step = bal_alpha_max_scale/real(bal_alpha_nscan,r8)
    have_anchor = .false.
    d_prev = ieee_value(0.0_r8,ieee_quiet_nan)
    scale_prev = 0.0_r8

    do iscan = 0, bal_alpha_nscan
      scale = step*real(iscan,r8)
      call bal_delta_at_scale(scale,d)
      if (.not. have_anchor) then
        if (ieee_is_finite(d) .and. d /= 0.0_r8) then
          have_anchor = .true.
          d_anchor = abs(d)
        end if
      else if (ieee_is_finite(d_prev) .and. ieee_is_finite(d) .and. &
          sign(1.0_r8,d) /= sign(1.0_r8,d_prev) .and. &
          max(abs(d_prev),abs(d)) <= pole_cap*d_anchor) then
        lo = scale_prev
        hi = scale
        sign_lo = sign(1.0_r8,d_prev)
        do while (abs(hi-lo) > bal_alpha_tol)
          mid = 0.5_r8*(lo+hi)
          call bal_delta_at_scale(mid,dm)
          if (ieee_is_finite(dm) .and. &
              sign(1.0_r8,dm) == sign_lo) then
            lo = mid
          else
            hi = mid
          end if
        end do
        locstab%fs(ipsi,5) = alpha_ref*0.5_r8*(lo+hi)
        return
      end if
      d_prev = d
      scale_prev = scale
    end do
  end subroutine bal_alpha_first_boundary

end module stride_bal_mod
