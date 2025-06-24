!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module torque
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Torque calculations.
    !
    !*DEPENDENCIES:
    !   - Must be called after inputs module read_dcon,read_kin establishes
    !   equilibrium (and perturbation) defining splines
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !   tpsi                    - complex NTV torque as a function of psi
    !
    !*PUBLIC DATA MEMBERS:
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------

    use params, only : r8, xj, mp, me, e, mu0, pi, twopi, &
        npsi_out, nell_out, nlambda_out, nmethods, methods, ngrids, grids, version
    use utilities, only : get_free_file_unit, check, median, append_2d, &
        ri, btoi, itob
    use special, only : ellipk,ellipe
    use grid_mod, only : powspace_sub,linspace_sub
    ! use lsode_mod just a subroutine in the lsode directory...
    use spline_mod, only :  spline_type,spline_eval,spline_alloc,spline_dealloc,&
                            spline_fit,spline_int,spline_write1,spline_eval_external,&
                            spline_roots
    use cspline_mod, only : cspline_type,cspline_eval,cspline_alloc,cspline_dealloc,&
                            cspline_fit,cspline_int,cspline_eval_external
    use bicube_mod, only : bicube_eval_external
    use pitch_integration, only : lambdaintgrl_lsode,kappaintgrl,kappadjsum
    use energy_integration, only : xintgrl_lsode,qt
    use dcon_interface, only : issurfint,mpert                  ! intel	< 2018 doesn't like mpert from inputs
    use inputs, only : eqfun,sq,geom,rzphi,smats,tmats,xmats,ymats,zmats,&
        kin,xs_m,dbob_m,divx_m,fnml, &                          ! equilib and pert. equilib splines
        chi1,ro,zo,bo,mfac,mthsurf,shotnum,shottime, &          ! reals or integers
        verbose, &                                              ! logical
        machine                                                 ! character
    use netcdf

    implicit none

    real(r8) :: atol_psi = 1e-3, rtol_psi= 1e-6
    real(r8) :: psi_warned = 0.0
    logical :: tdebug=.false., output_ascii =.false., output_netcdf=.true.
    integer :: nlmda=128, ntheta=128,nrecorded

    integer, parameter :: nfluxfuns = 19, nthetafuns = 12

    type record
        logical :: is_recorded
        real(r8), dimension(:,:), allocatable :: psi_record
        real(r8), dimension(:,:,:), allocatable :: ell_record
    endtype record
    type(record) :: torque_record(nmethods,ngrids)

    type diagnostic_record
        logical :: is_recorded
        integer :: psi_index, lambda_index, ell_index
        integer, dimension(:), allocatable :: ell
        real(r8), dimension(:), allocatable :: psi
        real(r8), dimension(:,:), allocatable :: fpsi
        real(r8), dimension(:,:,:), allocatable :: lambda
        real(r8), dimension(:,:,:,:,:), allocatable :: fs
    endtype diagnostic_record
    type(diagnostic_record) :: orbit_record(nmethods)

    complex(r8), dimension(:,:,:), allocatable :: elems
    TYPE(cspline_type) :: kelmm(6) ! kinetic euler lagrange matrix splines
    TYPE(cspline_type) :: trans ! nonambipolar transport and torque profile

    contains

    !=======================================================================
    subroutine tpsi(tpsi_var,psi,n,l,zi,mi,wdfac,divxfac,electron, &
                   method,op_erecord,op_ffuns,op_orecord,op_wmats)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Toroidal torque resulting from nonambipolar transport in perturbed
    !   equilibrium.
    !   Imaginary component is proportional to the kinetic energy Im(T) = 2*n*dW_k.
    !
    !*ARGUMENTS:
    !    psi : real (in)
    !       Normalized poloidal flux.
    !    n : integer (in)
    !       Toroidal mode number
    !    l : integer (in)
    !       Bounce harmonic number
    !    zi : integer (in)
    !       Ion charge in fundamental units (e).
    !    mi : integer (in)
    !       Ion mass (units of proton mass).
    !   electron : logical (in)
    !       Calculate quantities for electrons (zi,mi ignored)
    !   method : string (in)
    !       Choose from
    !           'RLAR' - Reduced Large-Aspect-Ratio,
    !           'CLAR' - Circular Large-Aspect-Ratio,
    !           '*GAR' - General-Aspect-Ratio
    !           '*TMM' - T_phi m-coupling matrix calculation.
    !           '*WMM' - dW_k m-coupling matrix calculation.
    !           '*KMM' - Euler-Lagrange Equation Matrix calculation.
    !                   -- Return value is sum of the 6 spectral norms (largest singular values)
    !       * = F,T,P for full,trapped,passing
    !       Note, these labels are also inserted in output file names.
    !
    !*OPTIONAL ARGUMENTS:
    !   op_erecord : logical
    !       Store energy space integrands in memory in the energy and pitch modules
    !   op_orecord : logical
    !       Store poloidal orbit functions in memory
    !   op_ffuns : real allocatable
    !       Store a record of relavent flux function quantities in this variable
    !   op_wmats : complex MxMx6
    !       Store a record of the 6 DCON matrice elements in this variable
    !       (requires extra computation)
    !
    !*RETURNS:
    !     complex.
    !        Toroidal torque do to nonambipolar transport.
    !-----------------------------------------------------------------------
        implicit none
        !declare function
        complex(r8), intent(out) :: tpsi_var
        ! declare arguments
        logical, intent(in) :: electron
        integer, intent(in) :: l,n,zi,mi
        real(r8), intent(in) :: psi,wdfac,divxfac
        character(*), intent(in) :: method
        logical, optional, intent(in) :: op_erecord, op_orecord
        real(r8), dimension(nfluxfuns), optional, intent(out) :: op_ffuns
        complex(r8), dimension(mpert,mpert,6), optional, intent(out) :: op_wmats
        ! declare local variables
        logical :: erecord, orecord
        integer :: i,j,k,s,ibmin,ibmax,sigma,ilmda,iqty,nbpts,nextrema
        integer, dimension(nlambda_out) :: ilambda_out
        real(r8) :: chrg,mass,welec,wdian,wdiat,wphi,wtran,wgyro,&
            wbhat,wdhat,nuk,nueff,q,epsr,bmin,bmax,lnq,theta,theta_bmin,theta_bmax,&
            lmdamin,lmdamax,lmdatpb,lmdatpe,lmda,t1,t2,&
            wbbar,wdbar,bhat,dhat,dbave,dxave,&
            vpar,kappaint,kappa,kk,djdj,jbb,&
            rex,imx,tnorm,he_t,hd_t,wb_t,wd_t
        real(r8), dimension(nthetafuns,ntheta) :: orbitfs
        real(r8), dimension(mthsurf*3) :: extrema,bpts
        real(r8), dimension(2,nlmda) :: ldl
        real(r8), dimension(2,2+nlmda) :: ldl_inc
        real(r8), dimension(2,2+nlmda/2) :: ldl_p
        real(r8), dimension(2,2+nlmda-nlmda/2) :: ldl_t
        real(r8), dimension(2,ntheta) :: tdt
        real(r8), dimension(:), allocatable :: dbfun,dxfun,fbnce_norm
        complex(r8) :: dbob,divx,kapx,xint,wtwnorm
        complex(r8), dimension(mpert) :: expm
        complex(r8), dimension(ntheta) :: jvtheta,pl
        complex(r8), dimension(1,1) :: t_zz,t_xx,t_yy,t_zx,t_zy,t_xy
        complex(r8), dimension(mpert,ntheta) :: wmu_mt,wen_mt
        complex(r8), dimension(1,mpert) :: wmmt,wemt,wxmt,wymt,wzmt
        complex(r8), dimension(mpert,mpert) :: smat,tmat,xmat,ymat,zmat
        complex(r8), dimension(mpert,1) :: wxmc,wymc,wzmc,xix,xiy,xiz
        complex(r8), dimension(:), allocatable :: lxint
        type(spline_type) :: tspl,vspl,bspl,cglspl,turns,dbdtspl
        type(cspline_type) :: bjspl,bwspl(2),fbnce
        ! for euclidean norm of wtw
        integer :: lwork,info
        real(r8), dimension(mpert) :: svals
        complex(r8), dimension(mpert,mpert) :: a,u,vt
        real(r8), dimension(5*mpert) :: rwork
        complex(r8), dimension(3*mpert) :: work

        ! for calling this in function parallel
        integer :: OMP_GET_THREAD_NUM
        integer :: ix, iy
        real(r8), dimension(3) :: geom_f, eqfun_f, eqfun_fx, eqfun_fy
        real(r8), dimension(4) :: rzphi_f, rzphi_fx, rzphi_fy, sq_s_f
        real(r8), dimension(9) :: kin_f, kin_f1
        complex(r8), dimension(mpert) :: xs_m1_f, xs_m2_f, xs_m3_f
        complex(r8), dimension(mpert) :: dbob_m_f, divx_m_f
        complex(r8), dimension(mpert**2) :: flatmat
        integer :: fsave
        real(r8) :: psave
        real(r8), dimension(0:mthsurf) :: jacs,delpsi,rsurf,asurf
        logical :: firstsurf

        ! debug initiation
        if(tdebug) print *,"torque - tpsi function, psi = ",psi
        if(tdebug) print *,"  electron ",electron
        if(tdebug) print *,"  ell ",l
        if(tdebug) print *,"  mpert ",mpert
        if(tdebug) print *,"  mfac ",mfac
        if(tdebug) print *,"  sq   psi ",sq%xs(0:sq%mx:sq%mx/5)
        if(tdebug) print *,"  dbob psi ",dbob_m%xs(0:dbob_m%mx:dbob_m%mx/5)
        if(tdebug) print *,"  db ",dbob_m%fs(0:dbob_m%mx:dbob_m%mx/5,1)
        if(tdebug) print *,"  xs ",xs_m(1)%fs(0:xs_m(1)%mx:xs_m(1)%mx/5,1)

        ! defaults of optional arguments
        if(present(op_erecord))then
            erecord = op_erecord
        else
            erecord = .false.
        endif
        if(present(op_orecord))then
            orecord = op_orecord
        else
            orecord = .false.
        endif

        ! enforce bounds
        if(psi>1) then
            tpsi_var = 0
            return
        endif

        ! set species
        if(electron)then
            chrg = -1*e
            mass = me
            s = 2
        else
            chrg = zi*e
            mass = mi*mp
            s = 1
        endif

        ! Get perturbations
        !cspline external evaluation
        CALL cspline_eval_external(dbob_m,psi,ix,dbob_m_f)
        CALL cspline_eval_external(divx_m,psi,ix,divx_m_f)

        !Poloidal functions - note ABS(A*clebsch) = ABS(A)
        allocate(dbfun(0:mthsurf),dxfun(0:mthsurf))
        call spline_alloc(tspl,mthsurf,5)
        tspl%xs(0:) = (/(i / real(mthsurf,r8), i=0, mthsurf)/) ! DCON theta normalized to unity
        do i=0,mthsurf
           theta = i/real(mthsurf,r8)
           call bicube_eval_external(eqfun, psi, theta, 1,&
                ix, iy, eqfun_f, eqfun_fx, eqfun_fy)
           call bicube_eval_external(rzphi, psi, theta, 1,&
                ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
           tspl%fs(i,1)= eqfun_f(1)            !b
           tspl%fs(i,2)= eqfun_fx(1)/chi1      !db/dpsi
           tspl%fs(i,3)= eqfun_fy(1)           !db/dtheta
           tspl%fs(i,4)= rzphi_f(4)/chi1       !jac
           tspl%fs(i,5)= rzphi_fx(4)/chi1**2   !dj/dpsi
           ! for flux fun outputs
           expm = exp(xj*twopi*mfac*theta)
           jbb  = rzphi_f(4)*eqfun_f(1)**2 ! chi1 cuz its the DCON working J
           dbfun(i) = ABS( sum(dbob_m_f(:)*expm) )**2
           dxfun(i) = ABS( sum(divx_m_f(:)*expm) * divxfac )**2
        enddo

        if (tdebug) then
#ifdef _OPENMP
        ! Compiled with OpenMP multithreading. Print out thread number.
           print *,"Pass1::: thread=",OMP_GET_THREAD_NUM(),"l=",l,"tspl%fs=",tspl%fs(:,1)
#endif
        endif

        ! clebsch conversion now in djdt o1*exp(-twopi*ifac*nn*q*theta)
        call spline_fit(tspl,"periodic")

        if (tdebug) then
#ifdef _OPENMP
        ! Compiled with OpenMP multithreading. Print out thread number.
           print *,"Pass2::: thread=",OMP_GET_THREAD_NUM(),"l=",l,"tspl%fs=",tspl%fs(:,1)
#endif
        endif

        ! rough estimate of bmin & bmax defining range of lambda
        bmax = maxval(tspl%fs(:,1),dim=1)
        bmin = minval(tspl%fs(:,1),dim=1)
        ibmax= 0-1+maxloc(tspl%fs(:,1),dim=1)
        if(bmax/=tspl%fs(ibmax,1)) stop "ERROR: tpsi_var - &
           &Equilibrium field maximum not consistent with index"
        do i=2,4 !4th smallest so spline has more than 1 pt
           bmin = minval(tspl%fs(:,1),mask=tspl%fs(:,1)>bmin,dim=1)
        enddo
        ibmin = 0-1+MINLOC(tspl%fs(:,1),MASK=tspl%fs(:,1)>=bmin,DIM=1)
        if(bmin/=tspl%fs(ibmin,1)) stop "ERROR: tpsi_var - &
           &Equilibrium field maximum not consistent with index"
        ! find precise bmin and bmax
        call spline_alloc(dbdtspl, mthsurf, 1)
        dbdtspl%xs(:) = tspl%xs(:)
        dbdtspl%fs(:, 1) = tspl%fs1(:, 1)
        call spline_fit(dbdtspl, "periodic")
        call spline_roots(dbdtspl, 1, nextrema, extrema)
        do i=1,nextrema
            call spline_eval(tspl, extrema(i), 0)
            if(tspl%f(1) < bmin)then
                bmin = tspl%f(1)
                theta_bmin = extrema(i)
             endif
             if (tspl%f(1) > bmax) then
                bmax = tspl%f(1)
                theta_bmax = extrema(i)
             end if
        end do
        call spline_dealloc(dbdtspl)
        if(tdebug) print *,"  bmin,bo,bmax = ", bmin, bo, bmax
        if(tdebug) print *,"  theta bmin,bmax = ", theta_bmin, theta_bmax  ! expect ~ 0, 0.5 (theta 0 is LFS)

        ! flux function variables  !!WARNING WHEN MODIFYING: THESE ARE CALCULATED SEPARATELY FOR I/O!!
        call spline_eval_external(sq,psi,ix,sq_s_f)
        call spline_eval_external(kin,psi,ix,kin_f,kin_f1)
        call spline_eval_external(geom,psi,ix,geom_f)
        q     = sq_s_f(4)
        welec = kin_f(5)                            ! electric precession
        wdian =-twopi*kin_f(s+2)*kin_f1(s)/(chrg*chi1*kin_f(s)) ! density gradient drive
        wdiat =-twopi*kin_f1(s+2)/(chrg*chi1)       ! temperature gradient drive
        wphi  = welec+wdian+wdiat                    ! toroidal rotation
        wtran = SQRT(2*kin_f(s+2)/mass)/(q*ro)      ! transit freq
        wgyro = chrg*bo/mass                        ! gyro frequency
        nuk = kin_f(s+6)                            ! krook collisionality

        call bicube_eval_external(rzphi, psi, theta_bmin, 0,&
             ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
        if(rzphi_f(1)<=0)then
           print *,"  psi = ",psi," -> r^2 at min(B) = ",rzphi_f(1)
           print *,"  -- theta at min(B) = ",theta_bmin
           do i=0,10
              theta = i/10.0
              call spline_eval(tspl,theta,0)
              call bicube_eval_external(rzphi, psi, theta, 0, ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
              print *,"  -- theta,B(theta),r^2(theta) = ",theta,tspl%f(1),rzphi_f(1)
           enddo
           stop "ERROR: torque - minor radius is negative"
        endif
        epsr = geom_f(2)/geom_f(3)
        wbhat = (pi/4)*SQRT(epsr/2)*wtran           ! RLAR normalized by x^1/2
        wdhat = q**3*wtran**2/(4*epsr*wgyro)*wdfac  ! RLAR normalized by x
        nueff = kin_f(s+6)/(2*epsr)

        firstsurf = .true.  ! first surface integral at this psi
        dbave = issurfint(dbfun,mthsurf,psi,0,1,fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)
        dxave = issurfint(dxfun,mthsurf,psi,0,1,fsave,psave,jacs,delpsi,rsurf,asurf,firstsurf)

        deallocate(dbfun,dxfun)
        if(tdebug) print('(a15,7(es10.1E2),i4)'), "   eq values = ",wdian,&
                        wdiat,welec,wdhat,wbhat,nueff,q

        ! optional record of flux functions
        if(present(op_ffuns)) then
            if(tdebug) print *, "  storing ",nfluxfuns," flux functions in ",size(op_ffuns,dim=1)
            op_ffuns = (/ epsr, kin_f(1:8), q, sq_s_f(2), &
                wdian, wdiat, wtran, wgyro, wbhat, wdhat, dbave, dxave /)
        endif

        ! optional orbit information at nlambda_out points in pitch angle
        if(orecord)then
            ilambda_out = (/(nint(ilmda*(nlmda-1)/real(nlambda_out-1))+1,ilmda=0,nlambda_out-1)/)
        endif



        if(tdebug) print *,'  method = '//method
        select case(method)

            case('fcgl')
                if(tdebug) print *,'  fcgl integrand. modes = ',n,l
                if(l==0)then
                    ! Poloidal functions (with jacobian!)
                    call spline_alloc(cglspl,mthsurf,2)
                    do i=0,mthsurf
                        theta = i/real(mthsurf,r8)
                        call bicube_eval_external(eqfun, psi, theta, 0,&
                             ix, iy, eqfun_f, eqfun_fx, eqfun_fy)
                        call bicube_eval_external(rzphi, psi, theta, 0,&
                             ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
                        expm = exp(xj*twopi*mfac*theta)
                        jbb = (rzphi_f(4) * eqfun_f(1)**2)
                        dbob = sum( dbob_m_f(:)*expm )   ! dB/B
                        divx = sum( divx_m_f(:)*expm ) ! nabla.xi_perp
                        kapx = -0.5*(dbob+divx)
                        cglspl%xs(i) = theta
                        ! include 0.5 correction for quadratic calculations in complex analysis
                        cglspl%fs(i,1) = rzphi_f(4)*0.5*abs(divx)**2
                        cglspl%fs(i,2) = rzphi_f(4)*0.5*abs(divx+3.0*kapx)**2
                    enddo
                    call spline_fit(cglspl,"periodic")
                    call spline_int(cglspl)
                    ! torque
                    tpsi_var = 2.0*n*xj*kin_f(s)*kin_f(s+2) &       ! T = 2nidW
                         *(0.5*(5.0/3.0)*cglspl%fsi(cglspl%mx,1)&
                         + 0.5*(1.0/3.0)*cglspl%fsi(cglspl%mx,2))
                    call spline_dealloc(cglspl)
               else
                    tpsi_var=0
               endif




            case("rlar")
                lnq = 1.0*l
                sigma = 0
                ! independent energy integration
                if(tdebug) print *,'  rlar integrations. modes = ',n,l
                ! energy space integrations
                lmdamax = min(1.0/(1-epsr),bo/bmin) ! just for record keeping - not rigorous
                xint = xintgrl_lsode(wdian,wdiat,welec,wdhat,wbhat,nueff,l,lnq,n,&
                    psi,lmdamax,method,op_record=erecord)
                ! kappa integration
                if(tdebug) print *,"  <|dB/B|> = ",sum(abs(dbob_m_f(:))/mpert)
                kappaint = kappaintgrl(n,l,q,mfac,dbob_m_f(:),fnml)
                ! dT/dpsi
                tpsi_var = sq_s_f(3)*kappaint*0.5*(-xint) &
                     *SQRT(epsr/(2*pi**3))*n*n*kin_f(s)*kin_f(s+2)
                if(tdebug) print *,'  ->  xint',xint,', kint',kappaint,', tpsi ',tpsi_var







            ! full general aspect ratio, trapped general aspect ratio
            case("clar")
                ! set up
                call cspline_alloc(fbnce,nlmda-1,3) ! <omegab,d> Lambda functions
                call spline_alloc(vspl,tspl%mx,1) ! vparallel(theta) -> roots are bounce pts
                vspl%xs(:) = tspl%xs(:)
                lmdatpb = bo/bmax
                lmdamin = max(1.0/(1+epsr),bo/bmax)
                lmdamax = min(1.0/(1-epsr),bo/bmin) ! kappa 0 to 1
                call powspace_sub(lmdamin,lmdamax,1,nlmda,"both",ldl) ! trapped space

                ! form smooth pitch angle functions
                call spline_alloc(turns,nlmda-1,6) ! (theta,r,z) of lower and upper turns
                do ilmda=1,nlmda
                    lmda = ldl(1,ilmda)
                    kk = (1 - lmda*(1- epsr))/(2*epsr*lmda)
                    if(kk<=0)then
                        kappa = 0
                    elseif(kk>=1)then
                        kappa = 1-1e-6
                    else
                        kappa = sqrt((1 - lmda*(1- epsr))/(2*epsr*lmda))
                    endif
                    lnq = 1.0*l
                    ! cylindrical bounce and precessions
                    wbbar = pi*SQRT(2*epsr*lmda*bo)/(4*q*ro*ellipk(kappa**2))
                    wdbar = (2*q*lmda*(ellipe(kappa**2)/ellipk(kappa**2)-0.5)&
                        /(ro**2*epsr))*wdfac
                    bhat = SQRT(2*kin_f(s+2)/mass)
                    dhat = (kin_f(s+2)/chrg)
                    ! perturbed action Eq. (12) [Park, Phys. Rev. Lett. 2009] divided by 2pi for DCON phi normalization
                    djdj = kappadjsum(kappa,n,l,q,mfac,dbob_m_f(:),fnml) * (0.5*pi/epsr) * (mass*bhat*q*ro)**2 / (2*pi)
                    ! bar normalization from [Logan, Phys. Plasmas 2013]
                    djdj = djdj / (mass*bhat*ro)**2
                    ! Lambda functions
                    fbnce%xs(ilmda-1) = lmda
                    fbnce%fs(ilmda-1,1) = wbbar*bhat
                    fbnce%fs(ilmda-1,2) = wdbar*dhat
                    fbnce%fs(ilmda-1,3) = wbbar*djdj
                    ! bounce locations recorded for optional output
                    vspl%fs(:,1) = 1.0-(lmda/bo)*tspl%fs(:,1)
                    call spline_fit(vspl,"extrap")
                    call spline_roots(vspl,1,nbpts,bpts)
                    if(nbpts<1)then
                        print *, "!! WARNING: Found passing particle in bounce particle integrals"
                        print *, "  >> Consider refining equilibrium or changing equil spline sizes"
                        t1 = 0.0
                        t2 = 1.0
                    else if(nbpts<2)then
                        t1 = bpts(1)
                        t1 = bpts(1) - 1.0
                    else
                        ! find deepest potential well
                        vpar = 0.0
                        do i=1,nbpts-1
                            call spline_eval(vspl,sum(bpts(i:i+1))/2,0)
                            if(i==1 .OR. vspl%f(1)>vpar)then
                                t1 = bpts(i)
                                t2 = bpts(i+1)
                                vpar = vspl%f(1)
                           endif
                        enddo
                        if (vpar==0) then
                            print *, "psi =",psi,", lambda =",lmda,", lambdamax =",lmdamax
                            stop "ERROR: Could not find potential well with positive vpar"
                        end if
                    end if
                    turns%fs(ilmda-1,1) = t1
                    call bicube_eval_external(rzphi, psi, t1, 0, ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
                    turns%fs(ilmda-1,2)=ro+SQRT(rzphi_f(1))*COS(twopi*(t1+rzphi_f(2)))
                    turns%fs(ilmda-1,3)=zo+SQRT(rzphi_f(1))*SIN(twopi*(t1+rzphi_f(2)))
                    turns%fs(ilmda-1,4) = t2
                    call bicube_eval_external(rzphi, psi, t2, 0, ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
                    turns%fs(ilmda-1,5)=ro+SQRT(rzphi_f(1))*COS(twopi*(t2+rzphi_f(2)))
                    turns%fs(ilmda-1,6)=zo+SQRT(rzphi_f(1))*SIN(twopi*(t2+rzphi_f(2)))
                    turns%xs(ilmda-1) = lmda
                enddo
                allocate(fbnce_norm(1))
                fbnce_norm(1) = 1/median(abs(fbnce%fs(:,3)))
                fbnce%fs(:,3) = fbnce%fs(:,3) * fbnce_norm(1)
                CALL cspline_fit(fbnce,'extrap')

                ! energy space integrations
                allocate(lxint(1))
                rex = 1.0
                imx = 1.0
                lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                    epsr,q,fbnce,l,n,rex,imx,psi,turns,method,op_record=erecord)

                ! dT/dpsi
                tpsi_var = (-2*n**2/sqrt(pi))*(ro/bo)*kin_f(s)*kin_f(s+2) &
                    *lxint(1)/fbnce_norm(1) &       ! lsode normalization
                    *(chi1/twopi)
                if(tdebug) print *,'  ->  lxint',lxint(1),', tpsi ',tpsi_var

                ! wrap up
                deallocate(fbnce_norm,lxint)
                call spline_dealloc(vspl)   ! vparallel(theta)
                call spline_dealloc(turns) ! bounce points
                call cspline_dealloc(fbnce) ! <omegab,d> Lambda functions







        ! full general aspect ratio, trapped general aspect ratio
            case("fgar","tgar","pgar", &
                 "fwmm","twmm","pwmm","ftmm","ttmm","ptmm",&
                 "fkmm","tkmm","pkmm","frmm","trmm","prmm")
                ! set up
                call spline_alloc(vspl,tspl%mx,1) ! vparallel(theta) -> roots are bounce pts
                call spline_alloc(bspl,ntheta-1,2) ! omegab,d bounce integration
                call cspline_alloc(bjspl,ntheta-1,1) ! action bounce integration
                vspl%xs(:) = tspl%xs(:)
                call linspace_sub(0.0_r8,1.0_r8,ntheta,bspl%xs(:))
                bjspl%xs(:)= bspl%xs(:)
                call spline_alloc(turns,nlmda-1,6) ! (theta,r,z) of lower and upper turns
                if(present(op_wmats))then
                    call cspline_alloc(fbnce,nlmda-1,3+mpert*mpert*6) ! <omegab,d>, <dJdJ>, <wtw> Lambda functions
                    do i=1,2
                        call cspline_alloc(bwspl(i),ntheta-1,mpert)
                        bwspl(i)%xs(:) = bjspl%xs(:)
                    enddo
                    ! build equilibrium geometric matrices
                    !Use the external array equivalent (for parallelization purposes, to avoid allocatable subcomponents)
                    ix = -1
                    call cspline_eval_external(smats,psi,ix,flatmat)
                    smat=reshape(flatmat,(/mpert,mpert/))
                    call cspline_eval_external(tmats,psi,ix,flatmat)
                    tmat=reshape(flatmat,(/mpert,mpert/))
                    call cspline_eval_external(xmats,psi,ix,flatmat)
                    xmat=reshape(flatmat,(/mpert,mpert/))
                    call cspline_eval_external(ymats,psi,ix,flatmat)
                    ymat=reshape(flatmat,(/mpert,mpert/))
                    call cspline_eval_external(zmats,psi,ix,flatmat)
                    zmat=reshape(flatmat,(/mpert,mpert/))
                else
                    call cspline_alloc(fbnce,nlmda-1,3) ! <omegab,d>, <dJdJ> Lambda functions
                endif
                fbnce%title(:) = "elemij"
                fbnce%title(0:3) = (/"Lambda","omegab","omegaD","wbdJdJ"/)
                lmdamin = 0.0
                lmdatpb = bo/bmax
                lmdatpe = min(bo/(tspl%fs(ibmax+1,1)),bo/(tspl%fs(ibmax-1,1)))
                lmdamax = bo/bmin
                if(method(1:1)=='t')then
                    call powspace_sub(lmdatpb,lmdamax,1,2+nlmda,"both",ldl_inc) ! trapped space including boundary
                    ldl = ldl_inc(:,2:1+nlmda) ! exclude boundary and max (both only have 1 bounce point)
                elseif(method(1:1)=='p')then
                    call powspace_sub(lmdamin,lmdatpb,1,2+nlmda,"both",ldl_inc) ! passing space including boundary
                    ldl = ldl_inc(:,2:1+nlmda) ! exclude boundary
                else
                    if(lmdatpb==lmdamax) then
                        print *, ""
                        print *,'!! WARNING: bmax = bmin @ psi',psi
                        print *, ""
                    end if
                    call powspace_sub(lmdamin,lmdatpb,2,2+nlmda/2,"upper",ldl_p) ! passing space including boundary
                    call powspace_sub(lmdatpb,lmdamax,2,2+nlmda-nlmda/2,"lower",ldl_t) ! trapped space including boundary
                    ldl(1,:) = (/ldl_p(1,2:1+nlmda/2),ldl_t(1,2:1+nlmda-nlmda/2)/) ! full space with no point on boundary or max
                    ldl(2,:) = (/ldl_p(2,2:1+nlmda/2),ldl_t(2,2:1+nlmda-nlmda/2)/)
                endif
                if(tdebug) print *," Lambda space ",ldl(1,1),ldl(1,nlmda),", t/p boundary = ",lmdatpb
                ! form smooth pitch angle functions
                do ilmda=1,nlmda
                    lmda = ldl(1,ilmda)
                    ! if(lmda==lmdatpb) lmda = lmda+1e-1*(ldl(1,2)-ldl(1,1)) ! tenth step off boundary
                    if(lmda>(bo/bmax)) then
                        sigma = 0 !trapped
                    else
                        sigma = 1 !passing
                    endif
                    lnq = l+sigma*n*q

                    ! determine bounce points
                    vspl%fs(:,1) = 1.0-(lmda/bo)*tspl%fs(:,1)
                    call spline_fit(vspl,"extrap")
                    if(sigma==0)then ! find trapped particle bounce pts
                        call spline_roots(vspl, 1, nbpts, bpts)
                        if(nbpts < 1)then
                            print *, "!!"
                            print *, "!! WARNING: Found passing particle in bounce particle integrals"
                            print *, "  > psi =",psi,", lambda=",lmda
                            print *, "  > lambda min, boundary, max =",lmdamin, lmdatpb, lmdamax
                            print *, "  > Consider refining equilibrium or changing equil spline sizes"
                            print *, "!!"
                            t1 = tspl%xs(ibmax)
                            t2 = tspl%xs(ibmax) + 1
                        else if(nbpts<2)then  ! marginally trapped
                            t1 = bpts(1)
                            t2 = bpts(1) + 1.0
                        else
                            ! find deepest potential well
                            vpar = 0.0
                            t1 = 0.0
                            t2 = 1.0
                            do i=1,nbpts
                                j = i + 1
                                if(j > nbpts) j = 1
                                if(bpts(i) > bpts(j)) then
                                    ! this actually the standard case since theta 0 is on the LFS
                                    ! so we expect vpar > 0 until bpts(1), negative to bpts(2), then positive back to 0
                                    ! in this case we actually want to integrate from bpts(2) to bpts(1)
                                    ! wrapping through the theta jump at theta = 0/1
                                    call spline_eval(vspl,modulo(0.5*(bpts(i) + bpts(j) + 1.0), 1.0_r8), 0)
                                else  ! normal case ordered and within 0-1
                                   call spline_eval(vspl,0.5*(bpts(i) + bpts(j)), 0)
                                end if
                                if(vspl%f(1) > vpar)then
                                    t1 = bpts(i)
                                    t2 = bpts(j)
                                    if(t2 < t1) t2 = t2 + 1.0
                                    vpar = vspl%f(1)
                               endif
                            enddo
                            if (vpar==0) print *, "ERROR: Could not find potential well with potive vpar"
                        end if
                        if(t1==t2)then
                            print *,"psi, lmda = ",psi,lmda
                            print *,"lmdamin, lambdatpb, lambdamax =",lmdamin,lmdatpb,lmdamax
                            print *,"t1 = t2 =",t1
                        end if
                        call powspace_sub(t1,t2,4,ntheta,"both",tdt)
                    else ! transit -> full theta integral
                        t1 = tspl%xs(ibmax)
                        t2 = tspl%xs(ibmax)+1
                        call powspace_sub(t1,t2,2,ntheta,"both",tdt)
                    endif
                    if(tdebug) then
                        if(mod(ilmda,ntheta/10)==0)then
                            print *,"(t1,t2) = ",t1,t2
                            print *,"tdt rng = ",tdt(1,1),tdt(1,ntheta)
                        endif
                    endif
                    ! bounce locations recorded for optional output
                    turns%fs(ilmda-1,1) = t1
                    call bicube_eval_external(rzphi, psi, t1, 0, ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
                    turns%fs(ilmda-1,2)=ro+SQRT(rzphi_f(1))*COS(twopi*(t1+rzphi_f(2)))
                    turns%fs(ilmda-1,3)=zo+SQRT(rzphi_f(1))*SIN(twopi*(t1+rzphi_f(2)))
                    turns%fs(ilmda-1,4) = t2
                    call bicube_eval_external(rzphi, psi, t2, 0, ix, iy, rzphi_f, rzphi_fx, rzphi_fy)
                    turns%fs(ilmda-1,5)=ro+SQRT(rzphi_f(1))*COS(twopi*(t2+rzphi_f(2)))
                    turns%fs(ilmda-1,6)=zo+SQRT(rzphi_f(1))*SIN(twopi*(t2+rzphi_f(2)))
                    turns%xs(ilmda-1) = lmda

                    ! bounce averages
                    wmu_mt = 0
                    wen_mt = 0
                    jvtheta(:) = 0
                    bspl%fs(:,:) = 0
                    do i=2,ntheta-1 ! edge dts are 0
                        call spline_eval(tspl,modulo(tdt(1,i),1.0_r8),0)
                        call spline_eval(vspl,modulo(tdt(1,i),1.0_r8),0)
                        vpar = vspl%f(1) ! more consistent w/ bnce pts than direct from tspl
                        if(vpar<=0)then ! local zero between nodes
                            if(lmda>lmdatpe .or. lmda<(2*lmdatpb-lmdatpe))then ! expected near t/p bounry
                                if(ABS(psi_warned-psi)>0.1)then !avoid flood of warnings
                                    print '(1x,a35,es10.3E2,a3,es10.3E2,a2,es10.3E2,a2,es10.3E2)', &
                                    "!! WARNING: vpar zero crossing at psi=",psi,",t=",t1,"<=",tdt(1,i),"<=",t2
                                    psi_warned = psi
                                endif
                            endif
                            if(i<ntheta/2)then
                                bspl%fs(:i-1,1) = 0
                                bspl%fs(:i-1,2) = 0
                                jvtheta(:i) = 0
                                cycle
                            else
                                bspl%fs(i-1:,1) = bspl%fs(i-2,1)
                                bspl%fs(i-1:,2) = bspl%fs(i-2,1)
                                jvtheta(i:) = jvtheta(i-1)
                                exit
                            endif
                        endif
                        bspl%fs(i-1,1) = tdt(2,i)*tspl%f(4)*tspl%f(1)/sqrt(vpar)
                        bspl%fs(i-1,2) = tdt(2,i)*tspl%f(4)*tspl%f(2)&
                            *(1-1.5*lmda*tspl%f(1)/bo)/sqrt(vpar)&
                            +tdt(2,i)*tspl%f(5)*tspl%f(1)*sqrt(vpar)
                        expm = exp(xj*twopi*mfac*tdt(1,i))
                        jbb = chi1*tspl%f(4)*tspl%f(1)**2 ! chi1 cuz its the DCON working J
                        dbob = sum(dbob_m_f(:)*expm)
                        divx = sum(divx_m_f(:)*expm) * divxfac
                        jvtheta(i) = tdt(2,i)*tspl%f(4)*tspl%f(1) &
                            *(divx*sqrt(vpar)+dbob*(1-1.5*lmda*tspl%f(1)/bo)/sqrt(vpar))&
                            *exp(-twopi*xj*n*q*(tdt(1,i)-tdt(1,1)))
                            ! theta0 doesn't really matter since |dj|^2
                        ! debuging
                        if(jvtheta(i)/=jvtheta(i))then
                            print *,'itheta,ntheta,theta',i,ntheta,tdt(1,i)
                            print *,'psi ',psi
                            print *,'dbob_m = ',dbob_m_f(:)
                            print *,'divx_m = ',divx_m_f(:)
                            print *,'vpar, action = ',vpar, jvtheta(i)
                            stop
                        endif
                        ! euler lagrange matrix vectors
                        if(present(op_wmats))then
                            !! chi1 comes from jac in s,t,x,y,zmats unconverted
                            wmu_mt(:,i) = tdt(2,i)*(lmda/bo)*expm/sqrt(vpar) &
                                *exp(-twopi*xj*n*q*(tdt(1,i)-tdt(1,1))) &
                                *1.0/(2*chi1) !! JKP s,t,x,y,zmats normalization factor??
                            wen_mt(:,i) = tdt(2,i)*expm/(tspl%f(1)*sqrt(vpar)) &
                                *exp(-twopi*xj*n*q*(tdt(1,i)-tdt(1,1))) &
                                *1.0/(2*chi1)  !! JKP s,t,x,y,zmats normalization factor??
                        endif

                        if(bspl%fs(i-2,1)==0)then ! smooth fill for pts beyond bounce
                            bspl%fs(2:i-2,1) = bspl%fs(i-1,1)
                            bspl%fs(2:i-2,2) = bspl%fs(i-1,2)
                            jvtheta(2:i-1) = jvtheta(i)
                        endif
                    enddo
                    call spline_fit(bspl,"extrap")
                    call spline_int(bspl)
                    ! Bounce averaged Lambda functions
                    fbnce%xs(ilmda-1) = lmda
                    wbbar = ro*twopi/((2-sigma)*bspl%fsi(bspl%mx,1))
                    wdbar = ro*ro*bo*wdfac*wbbar*2*(2-sigma)*bspl%fsi(bspl%mx,2)
                    bhat = sqrt(2*kin_f(s+2)/mass)/ro
                    dhat = (kin_f(s+2)/chrg)/(bo*ro*ro)
                    fbnce%fs(ilmda-1,1) = wbbar*bhat
                    fbnce%fs(ilmda-1,2) = wdbar*dhat
                    ! phase factor and action
                    !if(welec<fbnce%fs(ilmda-1,2))then ! mag prec. dominates h
                    !    pl=exp(-twopi*xj*lnq* bspl%fsi(0:,2)/((2-sigma)*bspl%fsi(bspl%mx,2)))
                    !else ! electric precession dominates drift (normal)
                    pl(:)=exp(-twopi*xj*lnq* bspl%fsi(0:,1)/((2-sigma)*bspl%fsi(bspl%mx,1)))
                    !endif
                    bjspl%fs(0:,1) =conjg(jvtheta(:))*(pl(:)+(1-sigma)/pl(:))
                    call cspline_fit(bjspl,"extrap")
                    call cspline_int(bjspl)
                    ! division by 2 corrects quadratic use of 2A_+n instead of proper A_+n + A_-n
                    fbnce%fs(ilmda-1,3) = wbbar * abs(bjspl%fsi(bjspl%mx, 1))**2 / 2 / ro**2

                    ! mxmx6 bounce averaged euler lagrange matrix elements
                    if(present(op_wmats))then
                        ! bounce integrate vectors W_mu,m and W_E,m
                        do i=1,mpert
                            bwspl(1)%fs(0:,i) = conjg(wmu_mt(i,:))*(pl(:)+(1-sigma)/pl(:))
                            bwspl(2)%fs(0:,i) = conjg(wen_mt(i,:))*(pl(:)+(1-sigma)/pl(:))
                        enddo
                        do i=1,2
                            call cspline_fit(bwspl(i),"extrap")
                            call cspline_int(bwspl(i))
                        enddo

                        ! build complete action mxm matrices W_X, W_Y, W_Z
                        wmmt(1,:) = bwspl(1)%fsi(bwspl(1)%mx,:)
                        wemt(1,:) = bwspl(2)%fsi(bwspl(2)%mx,:)
                        wxmt = matmul(wmmt,xmat)
                        wymt = matmul(wmmt,3*smat+ymat)-2.0*matmul(wemt,smat)
                        wzmt = matmul(wmmt,3*tmat+zmat)-2.0*matmul(wemt,tmat)
                        wxmc= conjg(transpose(wxmt))
                        wymc= conjg(transpose(wymt))
                        wzmc= conjg(transpose(wzmt))
                        op_wmats(:,:,1) = matmul(wzmc,wzmt)  !A
                        op_wmats(:,:,2) = matmul(wzmc,wxmt)  !B
                        op_wmats(:,:,3) = matmul(wzmc,wymt)  !C
                        op_wmats(:,:,4) = matmul(wxmc,wxmt)  !D
                        op_wmats(:,:,5) = matmul(wxmc,wymt)  !E
                        op_wmats(:,:,6) = matmul(wymc,wymt)  !H
                        do k=1,6
                            do j=1,mpert
                                do i=1,mpert
                                    iqty = ((k-1)*mpert+j-1)*mpert + i + 3
                                    fbnce%fs(ilmda-1,iqty) = wbbar*op_wmats(i,j,k)/ro**2
                                enddo
                            enddo
                        enddo
                    endif



                    ! optional deeply trapped bounce motion output for nlambda_out pitches
                    if(orecord .and. any(ilambda_out==ilmda))then
                        if(tdebug) print *, "  recording bounce functions"
                        do i=1,ntheta
                            expm = exp(xj*twopi*mfac*tdt(1,i))
                            dbob = sum(dbob_m_f(:)*expm)
                            divx = sum(divx_m_f(:)*expm) * divxfac
                            he_t = bspl%fsi(i-1,1)/((2-sigma)*bspl%fsi(bspl%mx,1))
                            hd_t = bspl%fsi(i-1,2)/((2-sigma)*bspl%fsi(bspl%mx,2))
                            wb_t = bspl%fs(i-1,1)
                            wd_t = bspl%fs(i-1,2)
                            orbitfs(:,i) = (/ tdt(1,i), tdt(2,i), &
                                real(bjspl%fs(i-1,1)), aimag((bjspl%fs(i-1,1))), &
                                real(dbob), aimag(dbob), real(divx), aimag(divx), &
                                wb_t, wd_t, he_t, hd_t /)
                        enddo
                        call record_orbit(method,psi,l,lmda,(/q,lmdatpb/),orbitfs)
                    endif


                enddo ! lambda loop
                ! normalize lambda functions for lsode (note: no resonant operator yet)
                allocate(fbnce_norm(fbnce%nqty-2))
                do i=1,fbnce%nqty-2
                    fbnce_norm(i) = 1.0/median(abs(fbnce%fs(:,i+2)))!maxval(abs(fbnce%fs(:,i+2)),1)!median(fbnce%fs(:,i+2))!
                    fbnce%fs(:,i+2) = fbnce%fs(:,i+2)*fbnce_norm(i)
                enddo
                call cspline_fit(fbnce,'extrap')
                call spline_fit(turns,'extrap')

                if(tdebug)then
                    print *,"lambda ~ ",ldl(1,::10),ldl(1,nlmda)
                    print *,"wb(lmda) ~ ",fbnce%fs(::10,1),fbnce%fs(nlmda-1,1)
                    print *,"wd(lmda) ~ ",fbnce%fs(::10,2),fbnce%fs(nlmda-1,2)
                    do i=fbnce%nqty/2,fbnce%nqty/2 !3,fbnce%nqty
                        print *,"dJdJ(lmda) ~ ",fbnce%fs(::10,i)
                    enddo
                endif


                ! energy space integrations
                allocate(lxint(fbnce%nqty-2))
                wtwnorm = 1.0 ! euler-lagrange real energy matrices (default)
                rex = 1.0 ! include real part of resonance operator (default)
                imx = 1.0 ! include imag part of resonance operator (default)
                if(method(2:4)=='wmm' .or. method(2:4)=='kmm')then
                    rex = 0.0
                elseif(method(2:4)=='tmm' .or. method(2:4)=='rmm')then
                    imx = 0.0
                    wtwnorm = -1.0
                endif
                lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                    epsr,q,fbnce,l,n,rex,imx,psi,turns,method,op_record=erecord)

                ! dT/dpsi
                tnorm = (-2 * n**2 / sqrt(pi)) * (ro / bo) * kin_f(s) * kin_f(s + 2) & ! Eq (19) [N.C. Logan, et al., Physics of Plasmas 20, (2013)]
                    * (chi1 / twopi) ! unit conversion from psi to psi_n, theta_n to theta
                tpsi_var = tnorm * ( lxint(1) / fbnce_norm(1) )       ! remove lsode normalization
                if(tdebug) print *,'  ->  lxint',lxint(1),', tpsi ',tpsi_var

                ! Euler lagrange matrices (wtw ~ dW ~ T/i)
                if(present(op_wmats))then
                    op_wmats(:,:,:) = 0
                    do k=1,6
                        do j=1,mpert
                            do i=1,mpert
                                iqty = ((k-1)*mpert+j-1)*mpert + i + 1
                                op_wmats(i,j,k) = lxint(iqty)/fbnce_norm(iqty) & ! remove lsode normalization
                                    * tnorm &               ! psi profile factors
                                    * (1 / (2*xj*n))        ! convert from torque to energy
                            enddo
                        enddo
                    enddo
                    call cspline_dealloc(bwspl(1))
                    call cspline_dealloc(bwspl(2))

                    ! DCON norms
                    op_wmats = 2*mu0*op_wmats
                    op_wmats(:,:,1:3) = op_wmats(:,:,1:3)/chi1
                    op_wmats(:,:,1) = op_wmats(:,:,1)/chi1

                    if(method(2:4)=='kmm' .or. method(2:4)=='rmm')then ! Euler-Lagrange Matrix norms indep xi
                        xix(:,1) = 1.0/sqrt(1.0*mpert) ! ||xix||=1
                        tpsi_var = 0
                        do i=1,6
                            !tpsi = tpsi + maxval(op_wmats(:,:,i))**2 ! Max m,m' couplings
                            !tpsi = tpsi + maxval(matmul(op_wmats(:,:,i),xix))**2 ! L_1 induced norms
                            !tpsi = tpsi + maxval(matmul(transpose(xix),op_wmats(:,:,i)))**2 ! L_inf induced norms
                            a=matmul(transpose(op_wmats(:,:,i)),op_wmats(:,:,i))
                            work=0
                            rwork=0
                            s=0
                            u=0
                            vt=0
                            lwork=3*mpert
                            call zgesvd('S','S',mpert,mpert,a,mpert,svals,u, &
                                mpert,vt,mpert,work,lwork,rwork,info)
                            if(tdebug) print *,i,maxval(svals)
                            tpsi_var = tpsi_var + maxval(svals) ! euclidean (spectral) norm
                        enddo
                        tpsi_var = (rex+imx*xj)*sqrt(tpsi_var) ! euclidean norm of the 6 norms
                    elseif(index(method,'mm')>0)then ! Mode-coupled dW of T_phi
                        call cspline_eval_external(xs_m(1),psi,ix,xs_m1_f)
                      call cspline_eval_external(xs_m(2),psi,ix,xs_m2_f)
                      call cspline_eval_external(xs_m(3),psi,ix,xs_m3_f)
                      xix(:,1) = xs_m1_f(:)
                      xiy(:,1) = xs_m2_f(:)
                      xiz(:,1) = xs_m3_f(:)
                        ! division by 2 corrects quadratic use of 2A_+n instead of proper A_+n + A_-n
                        t_zz = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,1),xiz))/2*chi1**2
                        t_zx = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,2),xix))/2*chi1
                        t_zy = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,3),xiy))/2*chi1
                        t_xx = matmul(conjg(transpose(xix)),matmul(op_wmats(:,:,4),xix))/2
                        t_xy = matmul(conjg(transpose(xix)),matmul(op_wmats(:,:,5),xiy))/2
                        t_yy = matmul(conjg(transpose(xiy)),matmul(op_wmats(:,:,6),xiy))/2
                        tpsi_var = (2*n*xj/(2*mu0))*(t_zz(1,1)+t_xx(1,1)+t_yy(1,1) &
                            +      t_zx(1,1)+t_zy(1,1)+t_xy(1,1) &
                            +wtwnorm*conjg(t_zx(1,1)+t_zy(1,1)+t_xy(1,1)))
                        if(tdebug)then
                            print *," -> WxWx ~ ",op_wmats(20:25,20,1)
                            print *,"expected to be all real (wmm) or all imag (tmm):"
                            print *,"  xx = ",t_xx
                            print *,"  yy = ",t_yy
                            print *,"  zz = ",t_zz
                        endif
                    endif

                endif

                ! wrap up
                deallocate(fbnce_norm,lxint)
                call spline_dealloc(vspl) ! vparallel(theta) -> roots are bounce pts
                call spline_dealloc(bspl) ! omegab,d bounce integration
                call spline_dealloc(turns) ! bounce points
                call cspline_dealloc(bjspl) ! action bounce integration
                call cspline_dealloc(fbnce) ! <omegab,d> Lambda functions


            case default
                stop "ERROR: torque - unknown method"
        end select

        if(tdebug) print *,"torque - end function, psi = ",psi

        return
    end subroutine tpsi

    !=======================================================================
    function tintgrl_grid(gtype,psilim,n,nl,zi,mi,wdfac,divxfac,electron,&
                          method)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Torque integratal over psi. This function forms a cubic spline of tpsi
    !   on the equilibrium grid taken from the DCON sq spline.
    !
    !*ARGUMENTS:
    !   gtype : string.
    !       Choose from
    !           - 'equil' for equilibrium spline grid from DCON
    !           - 'input' for input clebsch displacements grid (DCON Euler-Lagrange LSODE)
    !   psilim : real.
    !       (min,max) tuple specifying integration bounds.
    !   n : integer.
    !       mode number
    !   nl : integer.
    !       Number of bounce harmonics to include (-nl to nl)
    !    zi : integer (in)
    !       Ion charge in fundamental units (e).
    !    mi : integer (in)
    !       Ion mass (units of proton mass).
    !    wdfac : real (in)
    !       Add-hock multplied on magnetic precession.
    !    divxfac : real (in)
    !       Add-hock multiplier on divergence of perpendicular displacement.
    !   electron : logical
    !       Calculate quantities for electrons (zi,mi ignored)
    !   method : string
    !       Choose from 'RLAR', 'CLAR', 'FGAR', 'TGAR', 'PGAR', 'FWMM',
    !       'TWMM' or 'RWMM'.
    !       - also inserted in output file names.
    !
    !*RETURNS:
    !     complex. size of larray
    !        Integral int{dLambda f(Lambda)*int{dx Rln(x,Lambda)}, where
    !       f is input in the eq_spl and Rln is the resonance operator.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        complex(r8) :: tintgrl_grid
        ! declare arguments
        integer, intent(in) :: n,nl,zi,mi
        logical, intent(in) :: electron
        real(r8), intent(in) :: wdfac,divxfac
        real(r8), intent(inout) :: psilim(2)
        character(*) :: method,gtype
        ! declare variables
        integer :: i, j, l, s, mx, istrt, istop
        real(r8) :: x,xlast,chrg,drive,wdcom,dxcom
        real(r8), dimension(nfluxfuns) :: fcom
        real(r8), dimension(:), allocatable :: xs
        real(r8), dimension(:,:), allocatable :: profiles,ellprofiles
        complex(r8), dimension(:), allocatable :: gam,chi
        character(8) :: methcom
        ! lsode type variables
        integer  neqarray(6),neq
        real(r8), dimension(:), allocatable ::  y,dky
        ! declare new spline
        TYPE(cspline_type) :: tphi_spl

        common /tcom/ wdcom,dxcom,methcom,fcom

        ! set module variables
        psi_warned = 0.0
        ! set common variables
        wdcom = wdfac
        dxcom = divxfac
        methcom = method

        ! setup
        neq = 2*(1+2*nl)
        allocate(y(neq),dky(neq))
        neqarray(:) = (/neq,n,nl,zi,mi,btoi(electron)/)

        ! for flux calculations
        allocate(gam(2*nl+1),chi(2*nl+1))
        if(electron)then
            chrg = -1*e
            s = 2
        else
            chrg = zi*e
            s = 1
        endif

        ! allocate memory for kinetic DCON matrix profiles
        if(index(method,'mm')>0) then
            allocate(elems(mpert,mpert,6))
        endif

        ! set grid based on requested type
        if (gtype=='equil')then
            mx = sq%mx
            allocate(xs(0:mx))
            xs = sq%xs
        elseif (gtype=='input')then
            mx = xs_m(1)%mx
            allocate(xs(0:mx))
            xs = xs_m(1)%xs
        endif

        ! enforce integration bounds
        psilim(1) = max(psilim(1), sq%xs(0), xs_m(1)%xs(0))
        psilim(2) = min(psilim(2), sq%xs(sq%mx), xs_m(1)%xs(xs_m(1)%mx))
        istrt = -1
        istop = -1
        do i=0,mx
            if(xs(i)>=psilim(1)) exit
        enddo
        istrt = i
        do i=mx,0,-1
            if(xs(i)<=psilim(2)) exit
        enddo
        istop = i
        mx = istop-istrt

        ! prep allocations
        allocate(profiles(10+nfluxfuns,mx+1),ellprofiles(10,(1+2*nl)*(mx+1)))
        call cspline_alloc(tphi_spl,mx,2*nl+1)
        if(index(method,'mm')>0)then
            do j=1,6
                if(kelmm(j)%nqty/=0) call cspline_dealloc(kelmm(j))
                call cspline_alloc(kelmm(j),mx,mpert**2)
            enddo
        endif

        ! loop forming integrand
        xlast = 0
        do i=0,mx
            x = xs(i+istrt)
            ! torque profile
            call tintgrnd(neqarray,x,y,dky)
            tphi_spl%fs(i,:) = dky(1:neq:2)+dky(2:neq:2)*xj
            tphi_spl%xs(i) = x

            ! record of flux functions
            profiles(11:,i+1) = fcom
            ! save matrix of coefficients
            if(index(method,'mm')>0)then
                if(tdebug) print *, "Euler-Lagrange tmp vars, index=",i
                do j=1,6
                    kelmm(j)%xs(i) = x
                    kelmm(j)%fs(i,:) = reshape(elems(:,:,j),(/mpert**2/))
                enddo
            endif

            ! log progress
            if((mod(x,.1_r8)<mod(xlast,.1_r8) .or. xlast==xs(0)) .and. verbose)then
                print('(a7,f4.1,a13,2(es10.2e2),a1)'), " psi =",x,&
                    " -> dT/dpsi= ",sum(dky(1:neq:2),dim=1),sum(dky(2:neq:2),dim=1),'j'
            endif
            xlast = x
        enddo
        ! fit and integrate torque density spline
        call cspline_fit(tphi_spl,"extrap")
        call cspline_int(tphi_spl)

        ! flux/diffusivity profiles
        do i=0,mx
            x = xs(i+istrt)
            call spline_eval(sq,x,0)
            call spline_eval(kin,x,1)
            call spline_eval(geom,x,1)
            gam = tphi_spl%fs(i,:)*twopi/(chrg*chi1*geom%f(1))
            if(qt)then
                drive = (kin%f(s)/kin%f(s+2))*(kin%f1(s+2)/geom%f1(2)) ! (ndT/dr) /T
            else
                drive = kin%f1(s)/geom%f1(2)&   ! dn/dr
                       +(chrg*kin%f(s)/kin%f(s+2))*((kin%f(5)/twopi)/geom%f1(2)) ! (en/T)*dPhi/dr
            endif
            chi = -gam/(drive)

            ! save tables for ascii output
            profiles(:10,i+1) = (/ x, sq%f(3), &
                ri(sum(gam(:),dim=1)), &
                ri(sum(chi(:),dim=1)), &
                ri(sum(tphi_spl%fs(i,:),dim=1)), &
                ri(sum(tphi_spl%fsi(i,:),dim=1)) /)
            do l=-nl,nl
                ellprofiles(:,i*(1+2*nl) + (nl+l) + 1) = (/ x, l*1.0_r8, &
                    ri(gam(nl+l+1)), &
                    ri(chi(nl+l+1)), &
                    ri(tphi_spl%fs(i,nl+l+1)), &
                    ri(tphi_spl%fsi(i,nl+l+1)) /)
            enddo
        enddo

        ! (re-)set global transport and torque spline
        if(trans%nqty /= 0) call cspline_dealloc(trans)
        call cspline_alloc(trans,mx,2)
        trans%xs = profiles(1,:)
        trans%fs(:,1) = profiles(3,:)+xj*profiles(4,:) ! particle/heat flux gamma
        trans%fs(:,2) = profiles(7,:)+xj*profiles(8,:) ! torque & energy
        call cspline_fit(trans,"extrap")

        ! optionally fit global euler-lagrange coefficient splines
        if(index(method,'mm')>0)then
            do j=1,6
                call cspline_fit(kelmm(j),"extrap")
            enddo
        endif

        ! record results
        call record_method(trim(method),trim(gtype),profiles,ellprofiles)
        if(output_ascii)then
            call output_torque_ascii(n,zi,mi,electron,trim(method)//"_"//gtype//"_grid",profiles,ellprofiles)
        endif

        ! sum for ultimate result
        tintgrl_grid = sum(tphi_spl%fsi(mx,:),dim=1)

        deallocate(xs,y,dky,gam,chi,profiles,ellprofiles)
        call cspline_dealloc(tphi_spl)

        return
    end function tintgrl_grid

    !=======================================================================
    function tintgrl_lsode(psilim,n,nl,zi,mi,wdfac,divxfac,electron,method)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Torque integratal over psi. Integration boundes are set by the dcon
    !   equilibrium (taken from sq spline from inputs module's read_dcon).
    !
    !*ARGUMENTS:
    !   psilim : real.
    !       (min,max) tuple specifying integration bounds.
    !   n : integer.
    !       mode number
    !   nl : integer.
    !       Number of bounce harmonics to include (-nl to nl)
    !    zi : integer (in)
    !       Ion charge in fundamental units (e).
    !    mi : integer (in)
    !       Ion mass (units of proton mass).
    !    wdfac : real (in)
    !       Add-hock multplied on magnetic precession.
    !    divxfac : real (in)
    !       Add-hock multiplier on divergence of perpendicular displacement.
    !   electron : logical
    !       Calculate quantities for electrons (zi,mi ignored)
    !   method : string
    !       Choose from 'RLAR', 'CLAR', 'FGAR', 'TGAR', 'PGAR', 'FWMM',
    !       'TWMM' or 'RWMM'.
    !       - also inserted in output file names.
    !
    !*RETURNS:
    !     complex. size of larray
    !        Integral int{dLambda f(Lambda)*int{dx Rln(x,Lambda)}, where
    !       f is input in the eq_spl and Rln is the resonance operator.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        complex(r8) :: tintgrl_lsode
        ! declare arguments
        integer, intent(in) :: n,nl,zi,mi
        logical, intent(in) :: electron
        real(r8), intent(in) :: wdfac,divxfac
        real(r8), intent(inout) :: psilim(2)
        character(*) :: method
        ! declare variables
        integer, parameter :: maxsteps = 10000
        integer :: j,l,s
        real(r8) :: x0, xlast,wdcom,dxcom,chrg,drive
        real(r8), dimension(nfluxfuns) :: fcom
        real(r8), dimension(:), allocatable :: gam,chi
        real(r8), dimension(:,:), allocatable :: profiles,ellprofiles
        complex(r8), dimension(maxsteps,mpert**2,6) :: kel_flat_mats
        character(8) :: methcom
        ! declare lsode input variables
        integer  iopt, istate, itask, itol, mf, iflag,neqarray(6),&
            neq,liw,lrw
        integer, dimension(:), allocatable :: iwork
        real(r8) :: x,xout
        real(r8), dimension(:), allocatable ::  atol,rtol,rwork,y,dky

        common /tcom/ wdcom,dxcom,methcom,fcom

        ! set module variables
        psi_warned = 0.0
        ! set common variables
        wdcom = wdfac
        dxcom = divxfac
        methcom = method
        if(tdebug) print *, "torque - lsode wrapper method "//method//" -> "//methcom
        ! set lsode options - see lsode package for documentation
        neq = 2*(1+2*nl)        ! true number of equations
        liw  = 20 + neq         ! for mf 22 ! only uses 20 if mf 10
        lrw  = 20 + 16*neq      ! for mf 10
        allocate(iwork(liw))
        allocate(atol(neq),rtol(neq),rwork(lrw),y(neq),dky(neq))
        allocate(gam(neq),chi(neq))
        neqarray(:) = (/neq,n,nl,zi,mi,btoi(electron)/)
        y(:) = 0
        x    = max(psilim(1), sq%xs(0), xs_m(1)%xs(0)) ! xi grid and EQUIL grid might not match (DCON using qlow, etc.)
        x0 = x
        xout = min(psilim(2), sq%xs(sq%mx), xs_m(1)%xs(xs_m(1)%mx)) ! xi and equil grid might not match (psilim vs psihigh)
        itol = 2                  ! rtol and atol are arrays
        rtol(:) = rtol_psi              !1.d-7!9              ! 14
        atol(:) = atol_psi              !1.d-7!9              ! 15
        istate = 1                ! first step
        iopt = 1                  ! optional inputs
        iwork(:) = 0              ! defaults
        rwork(:) = 0              ! defaults
        rwork(1) = xout           ! only used if itask 4,5
        iwork(6) = maxsteps       ! max number steps
        mf = 10                   ! not stiff with unknown J
        if(tdebug) print *,"psilim = ",psilim
        if(tdebug) print *,"sq lim = ",sq%xs(0),sq%xs(sq%mx)
        if(tdebug) print *,"xs lim = ",xs_m(1)%xs(0),xs_m(1)%xs(xs_m(1)%mx)
        if(tdebug) print *,"x,xout = ",x,xout

        ! for flux calculation
        if(electron)then
            chrg = -1*e
            s = 2
        else
            chrg = zi*e
            s = 1
        endif

        ! integration
        itask = 5              ! single step without passing rwork(1)
        do while (x<xout)
            xlast = x
            call lsode(tintgrnd, neqarray, y, x, xout, itol, rtol,&
                atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
            call dintdy(x, 1, rwork(21), neq, dky, iflag)

            ! flux/diffusivity profiles
            call spline_eval(sq,x,0)
            call spline_eval(kin,x,1)
            call spline_eval(geom,x,1)
            gam = dky*twopi/(chrg*chi1*geom%f(1))
            if(qt)then
                drive = (kin%f(s)/kin%f(s+2))*(kin%f1(s+2)/geom%f1(2)) ! (ndT/dr) /T
            else
                drive = kin%f1(s)/geom%f1(2)&   ! dn/dr
                       +(chrg*kin%f(s)/kin%f(s+2))*((kin%f(5)/twopi)/geom%f1(2)) ! (en/T)*dPhi/dr
            endif
            chi = -gam/(drive)

            ! save tables for ascii output
            call append_2d(profiles, (/ x, sq%f(3), &
                sum(gam(1:neq:2),dim=1),sum(gam(2:neq:2),dim=1),&
                sum(chi(1:neq:2),dim=1),sum(chi(2:neq:2),dim=1),&
                sum(dky(1:neq:2),dim=1),sum(dky(2:neq:2),dim=1),&
                sum(  y(1:neq:2),dim=1),sum(  y(2:neq:2),dim=1),&
                fcom/) )
            do l=-nl,nl
                call append_2d(ellprofiles, (/ x, l*1.0_r8, &
                    gam(2*(nl+l)+1:2*(nl+l)+2), &
                    chi(2*(nl+l)+1:2*(nl+l)+2), &
                    dky(2*(nl+l)+1:2*(nl+l)+2), &
                      y(2*(nl+l)+1:2*(nl+l)+2) /) )
            enddo

            ! save matrix of coefficients
            if(index(method,'mm')>0)then
                if(tdebug) print *, "Euler-Lagrange tmp vars, iwork(11)=",iwork(11)
                do j=1,6
                    kel_flat_mats(iwork(11),:,j) = reshape(elems(:,:,j),(/mpert**2/))
                enddo
            endif

            ! print progress
            if((mod(x,.1_r8)<mod(xlast,.1_r8) .or. xlast==x0) .and. verbose)then
                print('(a7,f4.1,a13,2(es10.2e2),a1)'), " psi =",x,&
                    " -> T_phi = ",sum(y(1:neq:2),dim=1),sum(y(2:neq:2),dim=1),'j'
            endif
        enddo

        ! (re-)set global transport and torque spline
        if(trans%nqty /= 0) call cspline_dealloc(trans)
        call cspline_alloc(trans,iwork(11)-1,2)
        trans%xs = profiles(1,1:iwork(11))
        trans%fs(:,1) = profiles(3,:iwork(11))+xj*profiles(4,:iwork(11)) ! particle/heat flux gamma
        trans%fs(:,2) = profiles(7,:iwork(11))+xj*profiles(8,:iwork(11)) ! torque & energy
        call cspline_fit(trans,"extrap")

        ! optionally set global euler-lagrange coefficient splines
        if(index(method,'mm')>0)then
            do j=1,6
                if(kelmm(j)%nqty/=0) call cspline_dealloc(kelmm(j))
                call cspline_alloc(kelmm(j),iwork(11)-1,mpert**2)
                kelmm(j)%xs(:) = profiles(1,1:iwork(11))
                kelmm(j)%fs(:,:) = kel_flat_mats(1:iwork(11),:,j)
                call cspline_fit(kelmm(j),"extrap")
            enddo
        endif

        ! record results
        call record_method(trim(method),'lsode',profiles,ellprofiles)
        if(output_ascii)then
            call output_torque_ascii(n,zi,mi,electron,trim(method),profiles,ellprofiles)
        endif

        ! convert to complex space if integrations successful
        tintgrl_lsode = sum(y(1:neq:2),dim=1)+xj*sum(y(2:neq:2),dim=1)

        deallocate(iwork,rwork,atol,rtol,y,dky,gam,chi,profiles,ellprofiles)

        return
    end function tintgrl_lsode






    !=======================================================================
    subroutine tintgrnd(neq,x,y,ydot)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Wrapper routine for dynamic flux integration using lsode.
    !
    !*ARGUMENTS:
    !    neq : integer (in)
    !       Number of equations, must be 2 or greater.
    !    x : float (in)
    !       Normalized energy (E/T).
    !    y : real dimension(neq) (inout)
    !       First two elements are real and imaginary integral.
    !    ydot : real dimension(neq) (inout)
    !       First two elements are real and imaginary integrand.
    !
    !*RETURNS:
    !     complex.
    !        energy integral.
    !-----------------------------------------------------------------------
        implicit none
        integer ::  neq(*)
        real(r8) x, y(*), ydot(neq(1))

        real(r8) :: wdfac,xfac,psi
        real(r8), dimension(nfluxfuns) :: ffuns
        complex(r8) trq
        integer :: n,l,zi,mi,ee,nl
        logical :: electron=.false.,first=.true.
        character(8) :: method

        complex(r8), dimension (mpert,mpert,6) :: wtw_l

        common /tcom/ wdfac,xfac,method,ffuns

        !declarations for parallelization.
        integer :: omp_get_num_threads,omp_get_thread_num,lthreads


        if(tdebug .and. x<1e-2) print *, "torque - lsode subroutine wdfac ",wdfac
        if(tdebug .and. x<1e-2) print *, "torque - lsode subroutine method "//method

        n = neq(2)
        nl = neq(3)
        zi= neq(4)
        mi= neq(5)
        ee= neq(6)

        electron = itob(ee)

        if(first) then
            allocate(elems(mpert,mpert,6))
            first = .false.
        endif
        elems = 0

        psi = 1.0*x

        !$omp parallel default(shared) &
        !$omp& private(l,wtw_l,trq) &
        !$omp& reduction(+:elems) &
        !$omp& copyin(dbob_m,divx_m,kin,xs_m,fnml, &
        !$omp& geom, sq, eqfun, rzphi)
        
#ifdef _OPENMP
            IF(first .and. omp_get_thread_num() == 0)then
               lthreads = omp_get_num_threads()
               WRITE(*,'(1x,a,i3,a)') "Running in parallel with ",lthreads," OMP threads"
            ENDIF
#endif

        !$omp do
        do l=-nl,nl
            if(l==0)then
                if(index(method,'mm')>0)then
                    call tpsi(trq,psi,n,l,zi,mi,wdfac,xfac,electron,method,&
                               op_ffuns=ffuns,op_wmats=wtw_l)
                    elems = elems+wtw_l
                else
                    call tpsi(trq,psi,n,l,zi,mi,wdfac,xfac,electron,method,&
                               op_ffuns=ffuns)
                endif
            else
                if(index(method,'mm')>0)then
                    call tpsi(trq,psi,n,l,zi,mi,wdfac,xfac,electron,method,op_wmats=wtw_l)
                    elems = elems+wtw_l
                else
                    call tpsi(trq,psi,n,l,zi,mi,wdfac,xfac,electron,method)
                endif
            endif
            ! decouple two real space solutions
            ydot(2*(l+nl)+1) = real(trq)
            ydot(2*(l+nl)+2) = aimag(trq)
        enddo
        !$omp end do
        !$omp end parallel

        if(tdebug)then
            print *,'torque - intgrnd complete at psi ',x
            print *,'n,nl,zi,mi = ',n,nl,zi,mi
            print *,'x,psi,ydot = ',x,psi,ydot
        endif

        return
    end subroutine tintgrnd






    !=======================================================================
    subroutine noj(neq, t, y, ml, mu, pd, nrpd)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Dummy jacobian for use in lsode when true value unknown. See
    !   lsode module for details.
    !
    !*ARGUMENTS:
    !    neq : integer (in)
    !       Number of equations, must be 2 or greater.
    !    t : real (in)
    !       Normalized energy (E/T).
    !    y : real dimension(neq) (inout)
    !       First two elements are real and imaginary integral.
    !    ml : integer (in)
    !       Ignored.
    !    mu : integer.
    !       Ignored.
    !    pd : real dimension(nrpd,2).
    !       Set to 0.
    !    nrpd : integer.
    !       First dimension of pd.
    !
    !*OPTIONAL ARGUMENTS:
    !    wntedbk_ln : real, dimension(8) (in)
    !       wn : real.
    !           density gradient diamagnetic drift frequency
    !       wt : real.
    !           temperature gradient diamagnetic drift frequency
    !       we : real.
    !           electric precession frequency
    !       wd : real.
    !           magnetic precession frequency
    !       wb : real.
    !           bounce frequency divided by x (x=E/T)
    !       nuk : real.
    !           effective Krook collision frequency (i.e. /2eps for trapped)
    !       l : real.
    !           effective bounce harmonic (i.e. ell-nq for passing)
    !       n : real.
    !           toroidal mode number
    !
    !*RETURNS:
    !     complex.
    !        energy integral.
    !-----------------------------------------------------------------------
        implicit none
        integer  neq, ml, mu, nrpd
        real(r8)  t, y, pd(nrpd,2)
        ! null result
        pd(:,:) = 0
        return
    end subroutine noj


    !=======================================================================
    subroutine record_method(method,gridtype,prof,profl)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Store result of a integration to a module variable
    !
    !*ARGUMENTS:
    !    method : character
    !        Inserted into file name
    !    prof : real 2D
    !        Table of values
    !    profl : real 2D
    !        Table of (psi,ell) profiles
    !
    !*OPTIONAL ARGUMENTS:
    !
    !-----------------------------------------------------------------------
        implicit none
        character(*), intent(in) :: method,gridtype
        real(r8), dimension(:,:), intent(in) :: prof,profl

        integer :: i,j,nrow,ncol,nrowl,ncoll,nl,np

        if(tdebug) print *,"... Recording "//trim(method)//" on "//trim(gridtype)//" grid"

        ! learn sizes
        ncol = size(prof,dim=1)
        nrow = size(prof,dim=2)
        ncoll = size(profl,dim=1)
        nrowl = size(profl,dim=2)
        np = nrow
        if(np==0)then
            print *,"... failed attempt to record "//trim(method)//" on "//trim(gridtype)//" grid"
            return
        endif
        nl = nrowl/np
        if(nl*np/=nrowl) stop "ERROR: Unable to reshape torque output into P-by-L"

        ! allocate and store the method,grid records
        do i=1,ngrids
            if(gridtype==grids(i))then
                do j=1,nmethods
                    if(method==methods(j))then
                        torque_record(j,i)%is_recorded=.true.
                        allocate( torque_record(j,i)%psi_record(ncol,nrow) )
                        allocate( torque_record(j,i)%ell_record(ncoll,nl,np) )
                        torque_record(j,i)%psi_record = prof
                        torque_record(j,i)%ell_record = reshape(profl,(/ncoll,nl,np/))
                    endif
                enddo
            endif
        enddo

        return
    end subroutine record_method

    !=======================================================================
    subroutine record_orbit(method, psi, ell, lambda, fpsi, fs)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Save the bounce orbit theta functions for this method,
    !   at this surface and pitch.
    !
    !*ARGUMENTS:
    !   method : character
    !       Torque integrand method
    !   psi : real
    !       Normalized flux of integrand call
    !   ell : integer
    !       Bounce harmonic
    !   lambda : real
    !       Normalized pitch angle
    !   fpsi : real 2
    !       safety factor and trapped/Passing boundary for psi
    !   fs : real array
    !       Bounce orbit functions
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: ell
        real(r8), intent(in) :: psi,lambda,fpsi(2),fs(nthetafuns,ntheta)
        character(*), intent(in) :: method

        integer :: m, i, j, k
        logical, parameter :: debug = .false.

        if(debug) print *,"Recording method"

        ! find the right record
        do m=1,nmethods
            if(method==methods(m))then
                if(debug) print *,"  method "//trim(method)
                ! initialize the record of this method type if needed
                if(.not. orbit_record(m)%is_recorded)then
                    if(debug) print *,"   - is not recorded"
                    orbit_record(m)%is_recorded = .true.
                    orbit_record(m)%psi_index = 0
                    orbit_record(m)%lambda_index = 0
                    orbit_record(m)%ell_index = 0
                    allocate( orbit_record(m)%psi(npsi_out), &
                        orbit_record(m)%ell(nell_out), &
                        orbit_record(m)%lambda(nlambda_out,nell_out,npsi_out), &
                        orbit_record(m)%fpsi(2,npsi_out), &
                        orbit_record(m)%fs(nthetafuns,ntheta,nlambda_out,nell_out,npsi_out) )
                endif
                ! bump indexes
                if(debug) print *,"   - psi,ell,lambda = ",psi,ell,lambda
                if(debug) print *,"   - old i,j,k = ",i,j,k
                k = orbit_record(m)%psi_index
                j = orbit_record(m)%ell_index
                i = orbit_record(m)%lambda_index
                if(psi/=orbit_record(m)%psi(max(k, 1)) .or. k==0) k = k+1
                if(ell/=orbit_record(m)%ell(max(j, 1)) .or. j==0) j = mod(j,nell_out)+1
                if(lambda/=orbit_record(m)%lambda(max(i, 1),j,k) .or. i==0) i = mod(i,nlambda_out)+1
                if(debug) print *,"   - new i,j,k = ",i,j,k
                ! force fail if buggy
                if(k>npsi_out)then
                    print *,"ERROR: Too many psi energy records for method "//trim(method)
                    stop
                endif
                ! fill in new profiles
                orbit_record(m)%fpsi(:,k) = fpsi(:)
                orbit_record(m)%fs(:,:,i,j,k) = fs(:,:)
                orbit_record(m)%lambda(i,j,k) = lambda
                orbit_record(m)%ell(j) = ell
                orbit_record(m)%psi(k) = psi
                orbit_record(m)%psi_index = k
                orbit_record(m)%ell_index = j
                orbit_record(m)%lambda_index = i
            endif
        enddo
        if(debug) print *,"Done recording"
        return
    end subroutine record_orbit


    !=======================================================================
    subroutine reset_orbit_record()
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Deallocate records.
    !
    !*ARGUMENTS:
    !
    !
    !-----------------------------------------------------------------------
        implicit none

        integer :: m
        logical :: debug = .false.

        if(debug) print *,"Energy record resetting..."

        ! find the right record
        do m=1,nmethods
            ! initialize the record of this method type if needed
            if(orbit_record(m)%is_recorded)then
                orbit_record(m)%is_recorded = .false.
                orbit_record(m)%psi_index = 0
                orbit_record(m)%lambda_index = 0
                orbit_record(m)%ell_index = 0
                deallocate( orbit_record(m)%psi, &
                    orbit_record(m)%ell, &
                    orbit_record(m)%fpsi, &
                    orbit_record(m)%lambda, &
                    orbit_record(m)%fs )
            endif
        enddo
        if(debug) print *,"Energy record reset done"
        return
    end subroutine reset_orbit_record


    !=======================================================================
    subroutine output_torque_ascii(n,zi,mi,electron,method,prof,profl)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii torque profile files.
    !
    !*ARGUMENTS:
    !    n : integer
    !        Toroidal mode number for header
    !    zi : integer
    !        Ion charge for header
    !    mi : integer
    !        Ion mass for header
    !    electron : logical
    !        Modifies file name to indicate run was for electrons
    !    method : character
    !        Inserted into file name
    !    prof : real 2D
    !        Table of values
    !    op_profl : real 2D
    !        Table of (psi,ell) profiles
    !
    !*OPTIONAL ARGUMENTS:
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        character(*), intent(in) :: method
        logical :: electron
        real(r8), dimension(:,:), intent(in) :: prof,profl
        integer :: zi, mi

        integer :: i,s,nrow,ncol,nrowl,ncoll,out_unit,unit1,unit2
        real(r8) :: chrg,mass
        character(8) :: nstring
        character(32) :: fmt
        character(128) :: file1,file2

        ! set species
        if(electron)then
            chrg = -1*e
            mass = me
            s = 2
        else
            chrg = zi*e
            mass = mi*mp
            s = 2
        endif

        ! learn sizes
        ncol = size(prof,dim=1)
        nrow = size(prof,dim=2)
        ncoll = size(profl,dim=1)
        nrowl = size(profl,dim=2)

        ! open files
        write(nstring,'(I8)') n
        file1 = "pentrc_"//trim(method)//"_n"//trim(adjustl(nstring))//".out"
        file2 = "pentrc_"//trim(method)//"_ell_n"//trim(adjustl(nstring))//".out"
        if(qt)then
            file1 = file1(:7)//"heat_"//file1(8:)
            file2 = file2(:7)//"heat_"//file2(8:)
        endif
        unit1 = get_free_file_unit(-1)
        open(unit=unit1,file=trim(file1),status="unknown",action="write")
        unit2 = get_free_file_unit(-1)
        open(unit=unit2,file=trim(file2),status="unknown",action="write")

        ! write common header material
        do i=1,2
            if(i==1)then
                out_unit=unit1
            else
                out_unit=unit2
            endif
            write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
            write(out_unit,*)
            if(qt)then
                write(out_unit,*) " Gamma = Q/T = Heat flux in 1/(sm^2)"
                write(out_unit,*) " chi = Heat diffusivity in m^2/s"
                write(out_unit,*) " T_phi = Re[A*e*psi'*Gamma/2pi] in Nm (NOT NTV TORQUE!)"
                write(out_unit,*) " 2ndeltaW = Im[A*e*psi'*Gamma/2pi] in J (NOT dWk!)"
            else
                write(out_unit,*) " Gamma = Particle flux in 1/(sm^2)"
                write(out_unit,*) " chi = Particle diffusivity in m^2/s"
                write(out_unit,*) " T_phi = Torque in Nm"
                write(out_unit,*) " 2ndeltaW = Kinetic energy in J"
            endif
            write(out_unit,*) " * All energy dependent quantities are taken at E/T unity."
            write(out_unit,'(1/,1(a10,I4))') "n =",n
            write(out_unit,'(2(a10,es17.8E3))') "charge =",chrg,"mass = ",mass
            write(out_unit,'(3(a10,es17.8E3))') "R0 =",ro,"B0 =",bo,"chi1 = ",chi1
        enddo

        ! label columns
        if(ncol/=10+nfluxfuns) stop "ERROR: Torque ascii array dimensions do not match labels"
        write(fmt,*) '(1/,',10+nfluxfuns,'(a17))'
        write(unit1,fmt) "psi_n",  "dv/dpsi_n",  "real(Gamma)",  "imag(Gamma)", &
            "real(chi)",  "imag(chi)",  "T_phi",  "2ndeltaW",  "int(T_phi)",  "int(2ndeltaW)", &
            "eps_r","n_i","n_e","T_i","T_e","omega_E","logLambda","nu_i","nu_e","q",&
            "Pmu_0","omega_N","omega_T","omega_trans","omega_gyro","omega_b_rlar","omega_d_rlar",&
            "sqdBoB_L_SA","sqdivxi_perp_SA"!/J??

        if(ncoll/=10) stop "ERROR: Torque ascii ell array dimensions do not match labels"
        write(unit2,'(1/,10(a17))') "psi_n",  "ell",  "real(Gamma)",  "imag(Gamma)", &
            "real(chi)",  "imag(chi)",  "T_phi",  "2ndeltaW",  "int(T_phi)",  "int(2ndeltaW)"

        ! write tables
        write(fmt,*) '(',10+nfluxfuns,'(es17.8E3))'
        do i=1,nrow
            write(unit1,fmt) prof(:,i)
        enddo
        do i=1,nrowl
            write(unit2,'(10(es17.8E3))') profl(:,i)
        enddo

        close(unit1)
        close(unit2)
        return
    end subroutine output_torque_ascii

    !=======================================================================
    subroutine output_torque_netcdf(n,nl,zi,mi,electron,wdfac)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii torque profile files.
    !
    !*ARGUMENTS:
    !    n : integer
    !        Toroidal mode number for header
    !    zi : integer
    !        Ion charge for header
    !    mi : integer
    !        Ion mass for header
    !    electron : logical
    !        Modifies file name to indicate run was for electrons
    !    method : character
    !        Inserted into file name
    !    prof : real 2D
    !        Table of values
    !    op_profl : real 2D
    !        Table of (psi,ell) profiles
    !
    !*OPTIONAL ARGUMENTS:
    !
    !-----------------------------------------------------------------------
        implicit none
        real(r8), intent(in) :: wdfac
        integer, intent(in) :: n,nl,zi,mi
        logical, intent(in) :: electron

        integer :: i,j,s,np
        real(r8) :: chrg,mass,psi,q
        real(r8), dimension(:), allocatable :: &
            epsr,nuk,nueff,nui,nue,ni,ne,ti,te,llmda,zeff,&
            welec,wdian,wdiat,wphi,wtran,wgyro,wbhat,wdhat

        integer :: status, ncid,i_did,i_id,p_did,p_id,l_did,l_id, &
            v_id,g_id,c_id,d_id,t_id, b_id,x_id, er_id,q_id,mp_id, &
            ni_id,ne_id,ti_id,te_id,vi_id,ve_id,ze_id,ll_id,we_id, &
            wn_id,wt_id,ws_id,wg_id,wb_id,wd_id
        character(16) :: nstring,suffix
        character(128) :: ncfile

        if(verbose) print *,"Writing output to netcdf"

        ! set species
        if(electron)then
            chrg = -1*e
            mass = me
            s = 2
        else
            chrg = zi*e
            mass = mi*mp
            s = 1
        endif

        ! calculate the basic flux functions on the equilibrium grid
        allocate(epsr(sq%mx+1),nuk(sq%mx+1),nueff(sq%mx+1),nui(sq%mx+1),nue(sq%mx+1),&
            zeff(sq%mx+1),ni(sq%mx+1),ne(sq%mx+1),ti(sq%mx+1),te(sq%mx+1),llmda(sq%mx+1),&
            welec(sq%mx+1),wdian(sq%mx+1),wdiat(sq%mx+1),wphi(sq%mx+1),&
            wtran(sq%mx+1),wgyro(sq%mx+1),wbhat(sq%mx+1),wdhat(sq%mx+1))
        do i=1,sq%mx + 1
            psi = sq%xs(i-1)
            call spline_eval(sq,psi,1)
            call spline_eval(kin,psi,1)
            call spline_eval(geom,psi,0)
            q = sq%f(4)
            ni(i) = kin%f(1)
            ne(i) = kin%f(2)
            ti(i) = kin%f(3)
            te(i) = kin%f(4)
            welec(i) = kin%f(5)                            ! electric precession
            llmda(i) = kin%f(6)
            nui(i) = kin%f(7)
            nue(i) = kin%f(8)
            zeff(i) = kin%f(9)
            wdian(i) =-twopi*kin%f(s+2)*kin%f1(s)/(chrg*chi1*kin%f(s)) ! density gradient drive
            wdiat(i) =-twopi*kin%f1(s+2)/(chrg*chi1)       ! temperature gradient drive
            wphi(i)  = welec(i)+wdian(i)+wdiat(i)                    ! toroidal rotation
            wtran(i) = SQRT(2*kin%f(s+2)/mass)/(q*ro)      ! transit freq
            wgyro(i) = chrg*bo/mass                        ! gyro frequency
            nuk(i) = kin%f(s+6)                            ! krook collisionality
            epsr(i) = geom%f(2)/geom%f(3)
            wbhat(i) = (pi/4)*SQRT(epsr(i)/2)*wtran(i)           ! RLAR normalized by x^1/2
            wdhat(i) = q**3*wtran(i)**2/(4*epsr(i)*wgyro(i))*wdfac  ! RLAR normalized by x
            nueff(i) = kin%f(s+6)/(2*epsr(i))                 ! if trapped
        enddo

        ! create and open netcdf file
        write(nstring,'(I8)') n
        ncfile = "pentrc_output_n"//TRIM(ADJUSTL(nstring))//".nc"
        call check( nf90_create(ncfile, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
        ! store attributes
        call check( nf90_put_att(ncid, nf90_global, "title", "PENTRC fundamental outputs") )
        call check( nf90_put_att(ncid, nf90_global, "version", version))
        call check( nf90_put_att(ncid, nf90_global, "shot", INT(shotnum) ) )
        call check( nf90_put_att(ncid, nf90_global, "time", INT(shottime)) )
        call check( nf90_put_att(ncid, nf90_global, "machine", machine) )
        call check( nf90_put_att(ncid, nf90_global, "n", n) )
        call check( nf90_put_att(ncid, nf90_global, "charge", chrg) )
        call check( nf90_put_att(ncid, nf90_global, "mass", mass) )
        call check( nf90_put_att(ncid, nf90_global, "R0", ro) )
        call check( nf90_put_att(ncid, nf90_global, "B0", bo) )
        call check( nf90_put_att(ncid, nf90_global, "chi1", chi1) )
        ! define dimensions
        call check( nf90_def_dim(ncid, "i", 2, i_did) )
        call check( nf90_def_var(ncid, "i", nf90_int, i_did, i_id) )
        call check( nf90_def_dim(ncid, "ell", 2*nl+1, l_did) )
        call check( nf90_def_var(ncid, "ell", nf90_int, l_did, l_id) )
        call check( nf90_def_dim(ncid, "psi_n", sq%mx+1, p_did) )
        call check( nf90_def_var(ncid, "psi_n", nf90_double, p_did, p_id) )
        ! define variables
        call check( nf90_def_var(ncid, "dvdpsi", nf90_double, p_did, v_id) )
        call check( nf90_put_att(ncid, v_id, "long_name", "Differential Volume") )
        call check( nf90_put_att(ncid, v_id, "units", "m^3") )
        call check( nf90_def_var(ncid, "eps_r", nf90_double, p_did, er_id) )
        call check( nf90_put_att(ncid, er_id, "long_name", "Inverse Aspect Ratio") )
        call check( nf90_def_var(ncid, "q", nf90_double, p_did, q_id) )
        call check( nf90_put_att(ncid, q_id, "long_name", "Safety Factor") )
        call check( nf90_def_var(ncid, "mu0P", nf90_double, p_did, mp_id) )
        call check( nf90_put_att(ncid, mp_id, "long_name", "Equilibrium Pressure") )
        call check( nf90_def_var(ncid, "n_i", nf90_double, p_did, ni_id) )
        call check( nf90_put_att(ncid, ni_id, "long_name", "Ion Density") )
        call check( nf90_put_att(ncid, ni_id, "units", "m^-3") )
        call check( nf90_def_var(ncid, "n_e", nf90_double, p_did, ne_id) )
        call check( nf90_put_att(ncid, ne_id, "long_name", "Electron Density") )
        call check( nf90_put_att(ncid, ne_id, "units", "m^-3") )
        call check( nf90_def_var(ncid, "T_i", nf90_double, p_did, ti_id) )
        call check( nf90_put_att(ncid, ti_id, "long_name", "Ion Temperature") )
        call check( nf90_put_att(ncid, ti_id, "units", "eV") )
        call check( nf90_def_var(ncid, "T_e", nf90_double, p_did, te_id) )
        call check( nf90_put_att(ncid, te_id, "long_name", "Electron Temperature") )
        call check( nf90_put_att(ncid, te_id, "units", "eV") )
        call check( nf90_def_var(ncid, "logLambda", nf90_double, p_did, ll_id) )
        call check( nf90_put_att(ncid, ll_id, "long_name", "Logarithm of Plasma Parameter") )
        call check( nf90_def_var(ncid, "nu_i", nf90_double, p_did, vi_id) )
        call check( nf90_put_att(ncid, vi_id, "long_name", "Ion Collision Rate") )
        call check( nf90_put_att(ncid, vi_id, "units", "1/s") )
        call check( nf90_def_var(ncid, "nu_e", nf90_double, p_did, ve_id) )
        call check( nf90_put_att(ncid, ve_id, "long_name", "Electron Collision Rate") )
        call check( nf90_put_att(ncid, ve_id, "units", "1/s") )
        call check( nf90_def_var(ncid, "zeff", nf90_double, p_did, ze_id) )
        call check( nf90_put_att(ncid, ze_id, "long_name", "Effective Charge") )
        call check( nf90_def_var(ncid, "omega_E", nf90_double, p_did, we_id) )
        call check( nf90_put_att(ncid, we_id, "long_name", "Electric Precession Frequency") )
        call check( nf90_put_att(ncid, we_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_N", nf90_double, p_did, wn_id) )
        call check( nf90_put_att(ncid, wn_id, "long_name", "Density Gradient Diamagnetic Frequency") )
        call check( nf90_put_att(ncid, wn_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_T", nf90_double, p_did, wt_id) )
        call check( nf90_put_att(ncid, wt_id, "long_name", "Temperature Gradient Diamagnetic Frequency") )
        call check( nf90_put_att(ncid, wt_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_trans", nf90_double, p_did, ws_id) )
        call check( nf90_put_att(ncid, ws_id, "long_name", "Transit Frequency") )
        call check( nf90_put_att(ncid, ws_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_gyro", nf90_double, p_did, wg_id) )
        call check( nf90_put_att(ncid, wg_id, "long_name", "Gyro-Frequency") )
        call check( nf90_put_att(ncid, wg_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_b_rlar", nf90_double, p_did, wb_id) )
        call check( nf90_put_att(ncid, wb_id, "long_name", "Reduced Bounce Frequency") )
        call check( nf90_put_att(ncid, wb_id, "units", "rad/s") )
        call check( nf90_def_var(ncid, "omega_d_rlar", nf90_double, p_did, wd_id) )
        call check( nf90_put_att(ncid, wd_id, "long_name", "Reduced Magnetic Precession Frequency") )
        call check( nf90_put_att(ncid, wd_id, "units", "rad/s") )
        ! End definitions
        call check( nf90_enddef(ncid) )
        ! store variables
        call check( nf90_put_var(ncid, i_id, (/0,1/)) )
        call check( nf90_put_var(ncid, l_id, (/(i,i=-nl,nl)/)) )
        call check( nf90_put_var(ncid, p_id, sq%xs) )
        call check( nf90_put_var(ncid, v_id, sq%fs(:,3)) )
        call check( nf90_put_var(ncid, q_id, sq%fs(:,4)) )
        call check( nf90_put_var(ncid, mp_id, sq%fs(:,2)) )
        call check( nf90_put_var(ncid, er_id, epsr) )
        call check( nf90_put_var(ncid, ni_id, ni) )
        call check( nf90_put_var(ncid, ne_id, ne) )
        call check( nf90_put_var(ncid, ti_id, ti/1.602e-19) )
        call check( nf90_put_var(ncid, te_id, te/1.602e-19) )
        call check( nf90_put_var(ncid, we_id, welec) )
        call check( nf90_put_var(ncid, ll_id, llmda) )
        call check( nf90_put_var(ncid, vi_id, nui) )
        call check( nf90_put_var(ncid, ve_id, nue) )
        call check( nf90_put_var(ncid, ze_id, zeff) )
        call check( nf90_put_var(ncid, wn_id, wdian) )
        call check( nf90_put_var(ncid, wt_id, wdiat) )
        call check( nf90_put_var(ncid, ws_id, wtran) )
        call check( nf90_put_var(ncid, wg_id, wgyro) )
        call check( nf90_put_var(ncid, wb_id, wbhat) )
        call check( nf90_put_var(ncid, wd_id, wdhat) )
        deallocate(epsr,nuk,nueff,nui,nue,ni,ne,ti,te,llmda,&
            welec,wdian,wdiat,wphi,wtran,wgyro,wbhat,wdhat)

        ! store each method, grid combination that has been run
        do i=1,ngrids
            do j=1,nmethods
                if(torque_record(j,i)%is_recorded)then
                    ! create distinguishing labels
                    if(grids(i)=='lsode')then
                        suffix = '_'//trim(methods(j))
                    else
                        suffix = '_'//trim(methods(j))//'_'//trim(grids(i))
                    endif
                    ! check sizes
                    np = size(torque_record(j,i)%psi_record,dim=2)
                    if(np/=size(torque_record(j,i)%ell_record,dim=3))then
                        print *,SHAPE(torque_record(j,i)%psi_record)
                        print *,SHAPE(torque_record(j,i)%ell_record)
                        stop "Error: Record sizes are inconsistent"
                    endif
                    ! Re-open definitions
                    call check( nf90_redef(ncid) )
                    call check( nf90_put_att(ncid, nf90_global, "T_total"//trim(suffix), &
                                             torque_record(j,i)%psi_record(9,np)))
                    call check( nf90_put_att(ncid, nf90_global, "dW_total"//trim(suffix), &
                                             torque_record(j,i)%psi_record(10,np)/(2*n)))
                    if(trim(grids(i))=='lsode')then
                        call check( nf90_def_dim(ncid, "psi"//trim(suffix), np, p_did) )
                        call check( nf90_def_var(ncid, "psi"//trim(suffix), nf90_double, p_did, p_id) )
                        call check( nf90_def_var(ncid, "sqdBoB_L_SA"//trim(suffix), nf90_double, p_did, b_id) )
                        call check( nf90_def_var(ncid, "sqdivxi_perp_SA"//trim(suffix), nf90_double, p_did, x_id) )
                    else ! no need to store redundant flux functions on grids
                        status = nf90_inq_dimid(ncid, "psi_"//trim(grids(i)), p_did)
                        if(status == nf90_noerr )then ! re-use the psi for the repeat grid type
                            call check( nf90_inq_varid(ncid, "psi_"//trim(grids(i)), p_id) )
                            call check( nf90_inq_varid(ncid, "sqdBoB_L_SA_"//trim(grids(i)), b_id) )
                            call check( nf90_inq_varid(ncid, "sqdivxi_perp_SA_"//trim(grids(i)), x_id) )
                        else ! create the psi for the grid type
                            call check( nf90_def_dim(ncid, "psi_"//trim(grids(i)), np, p_did) )
                            call check( nf90_def_var(ncid, "psi_"//trim(grids(i)), nf90_double, p_did, p_id) )
                            call check( nf90_def_var(ncid, "sqdBoB_L_SA_"//trim(grids(i)), nf90_double, p_did, b_id) )
                            call check( nf90_def_var(ncid, "sqdivxi_perp_SA_"//trim(grids(i)), nf90_double, p_did, x_id) )
                        endif
                    endif
                    call check( nf90_def_var(ncid, "Gamma"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), g_id) )
                    call check( nf90_def_var(ncid, "chi"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), c_id) )
                    call check( nf90_def_var(ncid, "dTdpsi"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), d_id) )
                    call check( nf90_def_var(ncid, "T"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), t_id) )
                    ! Add metadata
                    if(qt)then
                        call check( nf90_put_att(ncid, g_id, "long_name", "Nonambipolar Heat Flux") )
                        call check( nf90_put_att(ncid, g_id, "units", "1/sm^2") )
                        call check( nf90_put_att(ncid, c_id, "long_name", "Nonambipolar Heat Diffusivity") )
                        call check( nf90_put_att(ncid, c_id, "units", "m^2/s") )
                        call check( nf90_put_att(ncid, d_id, "long_name", " A*e*psi'*Gamma/2pi profile") )
                        call check( nf90_put_att(ncid, d_id, "units", "Nm per psi") )
                        call check( nf90_put_att(ncid, t_id, "long_name", "Integrated A*e*psi'*Gamma/2pi") )
                        call check( nf90_put_att(ncid, t_id, "units", "Nm") )
                    else
                        call check( nf90_put_att(ncid, g_id, "long_name", "Nonambipolar Particle Flux") )
                        call check( nf90_put_att(ncid, g_id, "units", "1/sm^2") )
                        call check( nf90_put_att(ncid, c_id, "long_name", "Nonambipolar Particle Diffusivity") )
                        call check( nf90_put_att(ncid, c_id, "units", "m^2/s") )
                        call check( nf90_put_att(ncid, d_id, "long_name", "Toroidal Torque Profile") )
                        call check( nf90_put_att(ncid, d_id, "units", "Nm per unit normalized flux") )
                        call check( nf90_put_att(ncid, t_id, "long_name", "Integrated Toroidal Torque") )
                        call check( nf90_put_att(ncid, t_id, "units", "Nm") )
                    endif
                    call check( nf90_put_att(ncid, b_id, "long_name", "Surface Average Lagrangian Field Modulation") )
                    call check( nf90_put_att(ncid, x_id, "long_name", "Surface Average Perp. Displacement Divergence") )
                    call check( nf90_enddef(ncid) )
                    ! Put in variables
                    call check( nf90_put_var(ncid, p_id, torque_record(j,i)%ell_record(1,1,:)) )
                    call check( nf90_put_var(ncid, g_id, torque_record(j,i)%ell_record(3:4,:,:)) )
                    call check( nf90_put_var(ncid, c_id, torque_record(j,i)%ell_record(5:6,:,:)) )
                    call check( nf90_put_var(ncid, d_id, torque_record(j,i)%ell_record(7:8,:,:)) )
                    call check( nf90_put_var(ncid, t_id, torque_record(j,i)%ell_record(9:10,:,:)) )
                    call check( nf90_put_var(ncid, b_id, torque_record(j,i)%psi_record(28,:)) )
                    call check( nf90_put_var(ncid, x_id, torque_record(j,i)%psi_record(29,:)) )
                endif
            enddo
        enddo

        ! close file
        call check( nf90_close(ncid) )

        return
    end subroutine output_torque_netcdf


    !=======================================================================
    subroutine output_orbit_netcdf(n,op_label)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write recorded orbit functions to a netcdf file.
    !   Note, this is essentially a copy of output_energy_netcdf and the two
    !   should be kept consistent as much as possible.
    !
    !*ARGUMENTS:
    !    n : integer
    !        Toroidal mode number for filename
    !
    !*OPTIONAL ARGUMENTS:
    !    op_label : character
    !        Extra label inserted into filename
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        character(*), optional :: op_label

        integer :: i,m,npsi,nell,nlambda
        integer, dimension(:), allocatable :: ell_out
        real(r8), dimension(:), allocatable :: psi_out

        integer :: ncid,i_did,i_id,p_did,p_id,l_did,l_id, &
            a_did,a_id,x_did,x_id,aa_id,xx_id,dx_id, &
            dj_id,db_id,xi_id,wb_id,wd_id,he_id,hd_id,pq_id,at_id
        character(16) :: nstring,suffix,label
        character(128) :: ncfile

        logical :: debug = .false.

        if(verbose) print *,"Writing orbit record output to netcdf"

        ! optional labeling
        label = ''
        if(present(op_label)) label = "_"//trim(adjustl(op_label))

        ! assume all methods on the same psi's and ell's
        npsi = 0
        nell = 0
        nlambda = 0
        do m=1,nmethods
            if(orbit_record(m)%is_recorded)then
                if(debug) print *,"  Getting dims from "//trim(methods(m))
                npsi = orbit_record(m)%psi_index
                nell = orbit_record(m)%ell_index
                nlambda = orbit_record(m)%lambda_index
                allocate(psi_out(npsi),ell_out(nell))
                psi_out = orbit_record(m)%psi(:npsi)
                ell_out = orbit_record(m)%ell(:nell)
                exit
            endif
        enddo

        ! do nothing if no records
        if(npsi==0)then
            print *,"  WARNING: No energy records to output to netcdf"
            return
        endif

        ! create and open netcdf file
        write(nstring,'(I8)') n
        ncfile = "pentrc_orbit_output"//trim(label)//"_n"//trim(adjustl(nstring))//".nc"
        if(debug) print *, "  opening "//trim(ncfile)
        call check( nf90_create(ncfile, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid) )
        ! store attributes
        if(debug) print *, "  storing attributes"
        call check( nf90_put_att(ncid, nf90_global, "title", "PENTRC orbit outputs") )
        call check( nf90_put_att(ncid, nf90_global, "version", version))
        call check( nf90_put_att(ncid, nf90_global, "shot", INT(shotnum) ) )
        call check( nf90_put_att(ncid, nf90_global, "time", INT(shottime)) )
        call check( nf90_put_att(ncid, nf90_global, "machine", machine) )
        call check( nf90_put_att(ncid, nf90_global, "n", n) )
        ! define dimensions
        if(debug) print *, "  defining dimensions"
        if(debug) print *, "  npsi,nell,nlambda,nx = ",npsi,nell,nlambda,ntheta
        call check( nf90_def_dim(ncid, "i", 2, i_did) )
        call check( nf90_def_var(ncid, "i", nf90_int, i_did, i_id) )
        call check( nf90_def_dim(ncid, "ell", nell, l_did) )
        call check( nf90_def_var(ncid, "ell", nf90_int, l_did, l_id) )
        call check( nf90_def_dim(ncid, "psi_n", npsi, p_did) )
        call check( nf90_def_var(ncid, "psi_n", nf90_double, p_did, p_id) )
        call check( nf90_def_dim(ncid, "Lambda_index", nlambda, a_did) )
        call check( nf90_def_var(ncid, "Lambda_index", nf90_double, a_did, a_id) )
        call check( nf90_def_dim(ncid, "theta_index", ntheta, x_did) )
        call check( nf90_def_var(ncid, "theta_index", nf90_double, x_did, x_id) )
        ! End definitions
        call check( nf90_enddef(ncid) )
        ! store dimensions
        if(debug) print *, "  storing dimensions"
        call check( nf90_put_var(ncid, i_id, (/0,1/)) )
        call check( nf90_put_var(ncid, x_id, (/(i,i=1,ntheta)/)) )
        call check( nf90_put_var(ncid, a_id, (/(i,i=1,nlambda)/)) )
        call check( nf90_put_var(ncid, l_id, ell_out) )
        call check( nf90_put_var(ncid, p_id, psi_out) )

        ! store each method, grid combination that has been run
        do m=1,nmethods
            if(orbit_record(m)%is_recorded)then
                if(debug) print *, "  defining "//trim(methods(m))
                ! create distinguishing labels
                suffix = '_'//trim(methods(m))
                ! check sizes
                if(npsi/=orbit_record(m)%psi_index)then
                    print *,npsi,orbit_record(m)%psi_index
                    stop "Error: Record sizes are inconsistent"
                endif
                ! Re-open definitions
                call check( nf90_redef(ncid) )
                call check( nf90_def_var(ncid, "q"//trim(suffix), nf90_double, (/p_did/), pq_id) )
                call check( nf90_def_var(ncid, "Lambda_trap"//trim(suffix), nf90_double, (/p_did/), at_id) )
                call check( nf90_def_var(ncid, "Lambda"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), aa_id) )
                call check( nf90_def_var(ncid, "theta"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), xx_id) )
                call check( nf90_def_var(ncid, "dtheta"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), dx_id) )
                call check( nf90_def_var(ncid, "deltaJ"//trim(suffix), nf90_double, (/i_did,x_did,a_did,l_did,p_did/), dj_id) )
                call check( nf90_def_var(ncid, "b_lag_frac"//trim(suffix), nf90_double, (/i_did,x_did,a_did,l_did,p_did/), db_id) )
                call check( nf90_def_var(ncid, "divxi_perp"//trim(suffix), nf90_double, (/i_did,x_did,a_did,l_did,p_did/), xi_id) )
                call check( nf90_def_var(ncid, "omega_b"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), wb_id) )
                call check( nf90_def_var(ncid, "omega_D"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), wd_id) )
                call check( nf90_def_var(ncid, "h_E"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), he_id) )
                call check( nf90_def_var(ncid, "h_D"//trim(suffix), nf90_double, (/x_did,a_did,l_did,p_did/), hd_id) )
                call check( nf90_enddef(ncid) )
                ! Put in variables
                if(debug) print *, "  storing "//trim(methods(m))
                call check( nf90_put_var(ncid, pq_id, orbit_record(m)%fpsi(1,:npsi)) )
                call check( nf90_put_var(ncid, at_id, orbit_record(m)%fpsi(2,:npsi)) )
                call check( nf90_put_var(ncid, aa_id, orbit_record(m)%lambda(:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, aa_id, orbit_record(m)%lambda(:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, xx_id, orbit_record(m)%fs(1,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, dx_id, orbit_record(m)%fs(2,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, dj_id, orbit_record(m)%fs(3:4,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, db_id, orbit_record(m)%fs(5:6,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, xi_id, orbit_record(m)%fs(7:8,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, wb_id, orbit_record(m)%fs(9,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, wd_id, orbit_record(m)%fs(10,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, he_id, orbit_record(m)%fs(11,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, hd_id, orbit_record(m)%fs(12,:,:nlambda,:nell,:npsi)) )
            endif
        enddo

        ! close file
        call check( nf90_close(ncid) )
        ! clear the memory
        call reset_orbit_record( )
        if(debug) print *, "Finished energy netcdf output"

        return
    end subroutine output_orbit_netcdf


    !=======================================================================
    subroutine output_orbit_ascii(n, op_label)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii bounce function files.
    !
    !*ARGUMENTS:
    !    n : integer.
    !       mode number
    !*OPTIONAL ARGUMENTS:
    !    op_label : string
    !       Inserted into file name.
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        character(*), optional :: op_label

        integer :: m,i,j,k,istep,out_unit
        character(8) :: nstring
        character(16) :: label
        character(128) :: file

        ! optional labeling
        label = ''
        if(present(op_label)) label = "_"//trim(adjustl(op_label))

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstring,'(I8)') n
        file = "pentrc_orbit_output"//trim(label)//"_n"//trim(adjustl(nstring))//".out"
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Orbit outputs"
        write(out_unit,*) " - variables are:   lambda =  B0*m*v_perp^2/(2B)"
        write(out_unit,'(1/,1(a10,I4))') "n =",n

        ! write each method in a new table
        do m=1,nmethods
            if(orbit_record(m)%is_recorded)then
                write(out_unit,'(1/,2(a17))')"method =",trim(methods(m))
                write(out_unit,'(19(a17))') "psi_n","ell","Lambda_index","theta_index", &
                    "q","Lambda_trap","Lambda","theta","dtheta", &
                    "real(deltaJ)","imag(deltaJ)","real(b_lag_frac)","imag(b_lag_frac)", &
                    "real(divxi_perp)","imag(divxi_perp)","omega_b","omega_D","h_E","h_D"
                do k=1,orbit_record(m)%psi_index
                    do j=1,orbit_record(m)%ell_index
                        do i=1,orbit_record(m)%lambda_index
                            do istep=1,ntheta
                                write(out_unit,'(19(es17.8E3))') orbit_record(m)%psi(k), &
                                    real(orbit_record(m)%ell(j)),real(i),real(istep), &
                                    orbit_record(m)%fpsi(1,k),orbit_record(m)%fpsi(2,k), &
                                    orbit_record(m)%lambda(i,j,k),&
                                    orbit_record(m)%fs(:,istep,i,j,k)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        enddo
        close(out_unit)

        return
    end subroutine output_orbit_ascii


end module torque




