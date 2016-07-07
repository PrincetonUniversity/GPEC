!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE


module spline_help

    use params, only: r8
    use spline_mod, only: spline_type,spline_eval

    implicit none
    
    contains
    
    !=======================================================================
    subroutine spline_roots(roots,spl,iqty)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Find roots of a single spline quantity.
    !
    !*ARGUMENTS:
    !    roots : real, 1D allocatable (out)
    !       Zeros of the spline quantity.
    !    spl : spline_type (inout)
    !       Spline.
    !    iqty : integer (in)
    !       Spline quantity used.
    !
    !*RETURNS:
    !     real array.
    !        Zeros of spl%f(iqty).
    !-----------------------------------------------------------------------
        ! declarations.
        implicit none
        integer, intent(in) :: iqty
        type(spline_type) :: spl
        real(r8), dimension(:), allocatable :: roots
        ! declare variables
        integer :: iroot,it,ix,nroots
        integer, parameter :: itmax=20
        real(r8), parameter :: eps=1e-10
        real(r8) :: x,dx,lx,lf,f,df
    
        ! compute number of roots
        lx=spl%xs(spl%mx)-spl%xs(0)
        iroot = 1
        lf = maxval(spl%fs(:,iqty))-minval(spl%fs(:,iqty))
        nroots = 0
        do ix=1,spl%mx
           if(spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .le. 0.0) nroots=nroots+1
           !zeros counted twice
           if(spl%fs(ix,iqty)==0.0 .and. ix<spl%mx) nroots=nroots-1 
        enddo
        if(spl%periodic .and. spl%fs(spl%mx,iqty)==0) nroots=nroots-1
        if(allocated(roots)) deallocate(roots)
        allocate(roots(1:nroots))
    
        ! find all zero passings, intialize at larger gradient.
        do ix=1,spl%mx
           ! dont calculate exact zeros twice
           if(spl%fs(ix,iqty)==0.0 .and. ix<spl%mx) cycle
           if(spl%fs(ix,iqty)==0.0 .and. spl%periodic) cycle
           ! find crossing window
           if (spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .le. 0.0) then
              x = sum(spl%xs(ix-1:ix))/2.0 ! first estimate
              f=huge(f)
              dx=lx
              it=0
              ! locate roots by newton iteration.
              do
                 call spline_eval(spl,x,1)
                 df=spl%f(iqty)-f
                 if(abs(dx) < eps*lx .or. abs(df) < eps*lf .or. it >= itmax)exit
                 it=it+1
                 f=spl%f(iqty)
                 dx=-spl%f(iqty)/spl%f1(iqty)
                 x=x+dx
              enddo
              ! abort on failure.
              if(it >= itmax)then
                 !x=spl%xs(ix)
                 !if(spl%fs(ix,iqty) .lt. 0.0) x = spl%xs(ix-1)
                 write(*,*) "WARNING: roots convergence failure!"
                 write(*,*) "-- failure accured at index ",ix
                 write(*,*) "-- indexed x value ",spl%xs(ix)
                 write(*,*) "-- estimated root ",x
              endif
              !store each root, and allocate it to the spline.
              roots(iroot) = x
              iroot = iroot+1
           endif
        enddo
        return
        end subroutine spline_roots

end module spline_help


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
    
    use params, only : r8,xj,mp,me,e,mu0,pi,twopi
    use utilities, only : get_free_file_unit, check, median, append_2d, &
        ri, btoi, itob
    use special, only : ellipk,ellipe
    use grid, only : powspace,linspace
    ! use lsode_mod just a subroutine in the lsode directory...
    use spline_mod, only :  spline_type,spline_eval,spline_alloc,spline_dealloc,&
                            spline_fit,spline_int,spline_write1
    use cspline_mod, only : cspline_type,cspline_eval,cspline_alloc,cspline_dealloc,&
                            cspline_fit,cspline_int
    use fspline_mod, only : fspline_eval
    use bicube_mod, only : bicube_eval
    use spline_help, only: spline_roots
    use pitch_integration, only : lambdaintgrl_lsode,kappaintgrl,kappaintgnd
    use energy_integration, only : xintgrl_lsode,qt
    use dcon_interface, only : issurfint
    use inputs, only : eqfun,sq,geom,rzphi,smats,tmats,xmats,ymats,zmats,&
        kin,xs_m,dbob_m,divx_m,fnml, &                          ! equilib and pert. equilib splines
        chi1,ro,zo,bo,mfac,mpert,mthsurf,shotnum,shottime, &    ! reals or integers
        verbose, &                                              ! logical
        machine                                                 ! character
    use global_mod, only : version
    use netcdf

    implicit none
    
    real(r8) :: atol_psi = 1e-3, rtol_psi= 1e-6
    real(r8) :: psi_warned = 0.0
    logical :: tdebug=.false., output_ascii =.true., output_netcdf=.true.
    integer :: nlmda=128, ntheta=128,nrecorded
    
    integer, parameter :: nmethods = 5, ngrids = 3, nfluxfuns = 19, nthetafuns = 18
    character(4), dimension(nmethods) :: methods = (/'fcgl','rlar','clar','tgar','fgar'/)
    character(5), dimension(ngrids) :: grids = (/'lsode','equil','input'/)
    logical, dimension(nmethods,ngrids) :: isrecorded = .false.
    real(r8), dimension(:,:), allocatable, target :: &
        fcgl_lsode_psi_record, rlar_lsode_psi_record, clar_lsode_psi_record, &
        tgar_lsode_psi_record, fgar_lsode_psi_record, &
        fcgl_equil_psi_record, rlar_equil_psi_record, clar_equil_psi_record, &
        tgar_equil_psi_record, fgar_equil_psi_record, &
        fcgl_input_psi_record, rlar_input_psi_record, clar_input_psi_record, &
        tgar_input_psi_record, fgar_input_psi_record
    real(r8), dimension(:,:,:), allocatable, target :: &
        fcgl_lsode_ell_record, rlar_lsode_ell_record, clar_lsode_ell_record, &
        tgar_lsode_ell_record, fgar_lsode_ell_record, &
        fcgl_equil_ell_record, rlar_equil_ell_record, clar_equil_ell_record, &
        tgar_equil_ell_record, fgar_equil_ell_record, &
        fcgl_input_ell_record, rlar_input_ell_record, clar_input_ell_record, &
        tgar_input_ell_record, fgar_input_ell_record

    complex(r8), dimension(:,:,:), allocatable :: elems
    TYPE(cspline_type) :: kelmm(6) ! kinetic euler lagrange matrix splines
    TYPE(cspline_type) :: trans ! nonambipolar transport and torque profile

    contains
    
    !=======================================================================
    function tpsi(psi,n,l,zi,mi,wdfac,divxfac,electron,method,&
                    op_erecord,op_ffuns,op_tfuns,op_wmats)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Toroidal torque resulting from nonambipolar transport in perturbed
    !   equilibrium.
    !   Imaginary component isproportional to the kinetic energy Im(T) = 2*n*dW_k.
    !
    !*ARGUMENTS:
    !    psi : real (in)
    !       Normalized poloidal flux.
    !    n : integer (in)
    !       Toroidal mode number
    !    l : integer (in)
    !       Bounce harmonic number
    !    zi : integer (in)
    !       Ion charge in fundemental units (e).
    !    mi : integer (in)
    !       Ion mass (units of proton mass).
    !   electron : logical
    !       Calculate quantities for electrons (zi,mi ignored)
    !   method : string
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
    !   op_ffuns : real allocatable
    !       Store a record of relavent flux function quantities in this variable
    !   op_tfuns : real allocatable
    !       Store a record of bounce integrand quantities in this variable
    !       (requires extra computation)
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
        complex(r8) :: tpsi
        ! declare arguments
        logical, intent(in) :: electron
        integer, intent(in) :: l,n,zi,mi
        real(r8), intent(in) :: psi,wdfac,divxfac
        character(*) :: method
        logical, optional :: op_erecord
        real(r8), dimension(nfluxfuns), optional, intent(out) :: op_ffuns
        real(r8), dimension(nthetafuns,ntheta*3), optional, intent(out) :: op_tfuns
        complex(r8), dimension(mpert,mpert,6), optional, intent(out) :: op_wmats
        ! declare local variables
        logical :: erecord
        character(8) :: nstring,lstring
        character(32):: file
        integer :: i,j,k,s,ibmin,ibmax,out_unit,sigma,ilmda,iqty
        real(r8) :: chrg,mass,welec,wdian,wdiat,wphi,wtran,wgyro,&
            wbhat,wdhat,nuk,nueff,q,epsr,bmin,bmax,lnq,theta,&
            lmdamin,lmdamax,lmdatpb,lmdatpe,lmda,t1,t2,&
            wbbar,wdbar,bhat,dhat,dbave,dxave,&
            vpar,kappaint,kappa,kk,djdj,jbb,&
            rex,imx
        real(r8), dimension(2,nlmda) :: ldl
        real(r8), dimension(2,1+nlmda) :: ldl_inc
        real(r8), dimension(2,1+nlmda/2) :: ldl_p
        real(r8), dimension(2,1+nlmda-nlmda/2) :: ldl_t
        real(r8), dimension(2,ntheta) :: tdt
        real(r8), dimension(:), allocatable :: bpts,dbfun,dxfun
        complex(r8) :: dbob,divx,kapx,xint,wtwnorm
        complex(r8), dimension(mpert) :: expm
        complex(r8), dimension(ntheta) :: jvtheta,pl
        complex(r8), dimension(1,1) :: t_zz,t_xx,t_yy,t_zx,t_zy,t_xy
        complex(r8), dimension(mpert,ntheta) :: wmu_mt,wen_mt
        complex(r8), dimension(1,mpert) :: wmmt,wemt,wxmt,wymt,wzmt
        complex(r8), dimension(mpert,mpert) :: smat,tmat,xmat,ymat,zmat
        complex(r8), dimension(mpert,1) :: wxmc,wymc,wzmc,xix,xiy,xiz
        complex(r8), dimension(:), allocatable :: lxint,fbnce_norm
        type(spline_type) :: tspl,vspl,bspl,cglspl,turns
        type(cspline_type) :: bjspl,bwspl(2),fbnce
        ! for euclidean norm of wtw
        integer :: lwork,info
        real(r8), dimension(mpert) :: svals
        complex(r8), dimension(mpert,mpert) :: a,u,vt
        real(r8), dimension(5*mpert) :: rwork
        complex(r8), dimension(3*mpert) :: work
        
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

        ! enforce bounds
        if(psi>1) then
            tpsi = 0
            return
        endif

        ! set species
        if(electron)then
            chrg = e
            mass = me
            s = 2
        else
            chrg = zi*e
            mass = mi*mp
            s = 1
        endif
        
        ! Get perturbations
        call cspline_eval(dbob_m,psi,0)
        call cspline_eval(divx_m,psi,0)

        !Poloidal functions - note ABS(A*clebsch) = ABS(A)
        allocate(dbfun(0:eqfun%my),dxfun(0:eqfun%my))
        call spline_alloc(tspl,eqfun%my,5)
        tspl%xs(0:) = eqfun%ys(0:)
        do i=0,eqfun%my
           call bicube_eval(eqfun,psi,eqfun%ys(i),1)
           call bicube_eval(rzphi,psi,eqfun%ys(i),1)
           tspl%fs(i,1)= eqfun%f(1)            !b
           tspl%fs(i,2)= eqfun%fx(1)/chi1      !db/dpsi
           tspl%fs(i,3)= eqfun%fy(1)           !db/dtheta
           tspl%fs(i,4)= rzphi%f(4)/chi1       !jac
           tspl%fs(i,5)= rzphi%fx(4)/chi1**2   !dj/dpsi
           ! for flux fun outputs
           expm = exp(xj*twopi*mfac*eqfun%ys(i))
           jbb  = rzphi%f(4)*eqfun%f(1)**2 ! chi1 cuz its the DCON working J
           dbfun(i) = ABS( sum(dbob_m%f(:)*expm) )**2
           dxfun(i) = ABS( sum(divx_m%f(:)*expm) * divxfac )**2
        enddo
        ! clebsch conversion now in djdt o1*exp(-twopi*ifac*nn*q*theta)
        call spline_fit(tspl,"periodic")
        bmax = maxval(tspl%fs(:,1),dim=1)
        bmin = minval(tspl%fs(:,1),dim=1)
        ibmax= 0-1+maxloc(tspl%fs(:,1),dim=1)
        if(bmax/=tspl%fs(ibmax,1)) stop "ERROR: tpsi - &
           &Equilibirum field maximum not consistent with index"
        do i=2,4 !4th smallest so spline has more than 1 pt
           bmin = minval(tspl%fs(:,1),mask=tspl%fs(:,1)>bmin,dim=1)
        enddo
        ibmin = 0-1+MINLOC(tspl%fs(:,1),MASK=tspl%fs(:,1)>=bmin,DIM=1)
        if(bmin/=tspl%fs(ibmin,1)) stop "ERROR: tpsi - &
           &Equilibirum field maximum not consistent with index"
        if(tdebug) print *,"  bmin,bo,bmax = ",bmin,bo,bmax
        
        ! flux function variables
        call spline_eval(sq,psi,1)
        call spline_eval(kin,psi,1)
        call spline_eval(geom,psi,0)
        q     = sq%f(4)
        welec = kin%f(5)                            ! electric precession
        wdian =-twopi*kin%f(s+2)*kin%f1(s)/(chrg*chi1*kin%f(s)) ! density gradient drive
        wdiat =-twopi*kin%f1(s+2)/(chrg*chi1)       ! temperature gradient drive
        wphi  = welec+wdian+wdiat                    ! toroidal rotation
        wtran = SQRT(2*kin%f(s+2)/mass)/(q*ro)      ! transit freq
        wgyro = chrg*bo/mass                        ! gyro frequency
        nuk = kin%f(s+6)                            ! krook collisionality
        call bicube_eval(rzphi,psi,tspl%xs(ibmin),0)
        if(rzphi%f(1)<=0)then
            print *,"  psi = ",psi," -> r^2 at min(B) = ",rzphi%f(1)
            print *,"  -- theta at min(B) = ",tspl%xs(ibmin)
            do i=0,10
                theta = i/10.0
                call spline_eval(tspl,theta,0)
                call bicube_eval(rzphi,psi,theta,0)
                print *,"  -- theta,B(theta),r^2(theta) = ",theta,tspl%f(1),rzphi%f(1)
            enddo
            stop "ERROR: torque - minor radius is negative"
        endif
        epsr = geom%f(2)/geom%f(3)
        !epsr = sqrt(psi)*SQRT(SUM(rzphi%fs(rzphi%mx,:,1))/rzphi%my)/ro
        !epsr = sqrt(rzphi%f(1))/ro                  ! epsr at deep trapped limit
        wbhat = (pi/4)*SQRT(epsr/2)*wtran           ! RLAR normalized by x^1/2
        wdhat = q**3*wtran**2/(4*epsr*wgyro)*wdfac  ! RLAR normalized by x
        nueff = kin%f(s+6)/(2*epsr)                 ! if trapped
        dbave = issurfint(dbfun,eqfun%my,psi,0,1)
        dxave = issurfint(dxfun,eqfun%my,psi,0,1)
        deallocate(dbfun,dxfun)
        if(tdebug) print('(a14,7(es10.1E2),i4)'), "  eq values = ",wdian,&
                        wdiat,welec,wdhat,wbhat,nueff,q
        
        ! optional record of flux functions
        if(present(op_ffuns)) then
            if(tdebug) print *, "  storing ",nfluxfuns," flux functions in ",size(op_ffuns,dim=1)
            op_ffuns = (/ epsr, kin%f(1:8), q, sq%f(2), &
                wdian, wdiat, wtran, wgyro, wbhat, wdhat, dbave, dxave /)
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
                        call bicube_eval(eqfun,psi,theta,0)
                        call bicube_eval(rzphi,psi,theta,0)
                        expm = exp(xj*twopi*mfac*theta)
                        jbb = (rzphi%f(4) * eqfun%f(1)**2)
                        dbob = sum( dbob_m%f(:)*expm )   ! dB/B
                        divx = sum( divx_m%f(:)*expm ) ! nabla.xi_perp
                        kapx = -0.5*(dbob+divx)
                        cglspl%xs(i) = theta
                        cglspl%fs(i,1) = rzphi%f(4)*divx*CONJG(divx)
                        cglspl%fs(i,2) = rzphi%f(4)*(divx+3.0*kapx)*CONJG(divx+3.0*kapx)
                    enddo
                    CALL spline_fit(cglspl,"periodic")
                    CALL spline_int(cglspl)
                    ! torque
                    tpsi = 2.0*n*xj*kin%f(s)*kin%f(s+2) &       ! T = 2nidW
                        *(0.5*(5.0/3.0)*cglspl%fsi(cglspl%mx,1)&
                        + 0.5*(1.0/3.0)*cglspl%fsi(cglspl%mx,2))
                    call spline_dealloc(cglspl)
               else
                    tpsi=0
               endif
     
     
     
        
            case("rlar")
                lnq = 1.0*l
                ! independent energy integration
                if(tdebug) print *,'  rlar integrations. modes = ',n,l
                ! energy space integrations
                if(erecord)then
                    xint = xintgrl_lsode(wdian,wdiat,welec,wdhat,wbhat,nueff,lnq,n,&
                        (/psi,0.0_r8/))
                else
                    xint = xintgrl_lsode(wdian,wdiat,welec,wdhat,wbhat,nueff,lnq,n)
                endif
                ! kappa integration
                if(tdebug) print *,"  <|dB/B|> = ",sum(abs(dbob_m%f(:))/mpert)
                kappaint = kappaintgrl(n,l,q,mfac,dbob_m%f(:),fnml)
                ! dT/dpsi
                tpsi = sq%f(3)*kappaint*0.5*(-xint) &
                    *SQRT(epsr/(2*pi**3))*n*n*kin%f(s)*kin%f(s+2)
                if(tdebug) print *,'  ->  xint',xint,', kint',kappaint,', tpsi ',tpsi
            
            
            
            
            
            
            
            ! full general aspect ratio, trapped general aspect ratio
            case("clar")
                ! set up
                call cspline_alloc(fbnce,nlmda-1,3) ! <omegab,d> Lambda functions
                lmdatpb = bo/bmax
                lmdamin = max(1.0/(1+epsr),bo/bmax)
                lmdamax = min(1.0/(1-epsr),bo/bmin) ! kappa 0 to 1
                ldl = powspace(lmdamin,lmdamax,1,nlmda,"both") ! trapped space

                ! form smooth pitch angle functions
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
                    bhat = SQRT(2*kin%f(s+2)/mass)
                    dhat = (kin%f(s+2)/chrg)
                    ! JKP PRL 2009 perturbed action
                    djdj = kappaintgnd(kappa,n,l,q,mfac,dbob_m%f(:),fnml)
                    ! Lambda functions
                    fbnce%xs(ilmda-1) = lmda
                    fbnce%fs(ilmda-1,1) = wbbar*bhat
                    fbnce%fs(ilmda-1,2) = wdbar*dhat
                    fbnce%fs(ilmda-1,3) = djdj
                enddo
                CALL cspline_fit(fbnce,'extrap')
                
                ! energy space integrations
                allocate(lxint(1))
                if(erecord)then
                    lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                        epsr,q,fbnce,l,n,op_psi=psi)
                else
                    lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                        epsr,q,fbnce,l,n)
                endif

                ! dT/dpsi
                tpsi = sq%f(3)*(-1*lxint(1)) & ! may be missing some normalizations from djdj
                    *SQRT(epsr/(2*pi**3))*n*n*kin%f(s)*kin%f(s+2)
                if(tdebug) print *,'  ->  lxint',lxint(1),', tpsi ',tpsi
                
                ! wrap up
                deallocate(lxint)
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
                bspl%xs(:) = linspace(0.0_r8,1.0_r8,ntheta)
                bjspl%xs(:)= bspl%xs(:)
                call spline_alloc(turns,nlmda-1,6) ! (theta,r,z) of lower and upper turns               
                if(present(op_wmats))then
                    call cspline_alloc(fbnce,nlmda-1,3+mpert*mpert*6) ! <omegab,d>, <dJdJ>, <wtw> Lambda functions
                    do i=1,2
                        call cspline_alloc(bwspl(i),ntheta-1,mpert)
                        bwspl(i)%xs(:) = bjspl%xs(:)
                    enddo
                    ! build equilibrium geometric matrices 
                    CALL cspline_eval(smats,psi,0)
                    CALL cspline_eval(tmats,psi,0)
                    CALL cspline_eval(xmats,psi,0)
                    CALL cspline_eval(ymats,psi,0)
                    CALL cspline_eval(zmats,psi,0)
                    smat=RESHAPE(smats%f,(/mpert,mpert/))
                    tmat=RESHAPE(tmats%f,(/mpert,mpert/))
                    xmat=RESHAPE(xmats%f,(/mpert,mpert/))
                    ymat=RESHAPE(ymats%f,(/mpert,mpert/))
                    zmat=RESHAPE(zmats%f,(/mpert,mpert/))
                else                
                    call cspline_alloc(fbnce,nlmda-1,3) ! <omegab,d>, <dJdJ> Lambda functions
                endif
                fbnce%title(:) = "elemij"
                fbnce%title(0:3) = (/"Lambda","omegab","omegaD","wbdJdJ"/)
                lmdamin = 0.0
                lmdatpb = bo/bmax
                lmdatpe = min(bo/(tspl%fs(ibmax+1,1)),bo/(tspl%fs(ibmax-1,1)))
                lmdamax = bo/bmin ! has been reduced by 8 grid pts
                if(method(1:1)=='t')then
                    ldl_inc = powspace(lmdatpb,lmdamax,1,1+nlmda,"both") ! trapped space including boundary
                    ldl = ldl_inc(:,2:) ! exclude boundary
                elseif(method(1:1)=='p')then
                    ldl_inc = powspace(lmdamin,lmdatpb,1,1+nlmda,"both") ! passing space including boundary
                    ldl = ldl_inc(:,:nlmda) ! exclude boundary
                else
                    ldl_p = powspace(lmdamin,lmdatpb,2,1+nlmda/2,"upper") ! passing space including boundary
                    ldl_t = powspace(lmdatpb,lmdamax,2,1+nlmda-nlmda/2,"lower") ! trapped space including boundary
                    ldl(1,:) = (/ldl_p(1,:nlmda/2),ldl_t(1,2:)/) ! full space with no point on boundary
                    ldl(2,:) = (/ldl_p(2,:nlmda/2),ldl_t(2,2:)/)
                endif
                if(tdebug) print *," Lambda space ",ldl(1,1),ldl(1,nlmda),", t/p bounry = ",lmdatpb
                ! form smooth pitch angle functions
                do ilmda=1,nlmda
                    lmda = ldl(1,ilmda)
                    !if(lmda==lmdatpb) lmda = lmda+1e-1*(ldl(1,2)-ldl(1,1)) ! tenth step off boundry
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
                        call spline_roots(bpts,vspl,1)
                        if(size(bpts)>2 .and. ilmda==1) then
                            print *, "WARNING: using only deepest of &
                                &multiple magnetic wells at psi ",psi
                        endif
                        ! find deepest magnetic well ** not precise **
                        t1 = bpts(size(bpts))-1.0  
                        t2 = bpts(1)
                        call spline_eval(vspl,modulo((t1+t2)/2,1.0_r8),0)
                        vpar = vspl%f(1) ! bpts centered around 0 (standard)
                        do i=1,size(bpts)-1
                           call spline_eval(vspl,sum(bpts(i:i+1))/2,0)
                           if(vspl%f(1)>vpar)then
                              t1 = bpts(i)
                              t2 = bpts(i+1)
                              vpar = vspl%f(1)
                           endif
                        enddo
                        tdt = powspace(t1,t2,4,ntheta,"both")
                        deallocate(bpts)
                    else ! transit -> full theta integral
                        t1 = tspl%xs(ibmax)
                        t2 = tspl%xs(ibmax)+1
                        tdt = powspace(t1,t2,2,ntheta,"both")
                    endif
                    if(tdebug) then
                        if(mod(ilmda,ntheta/10)==0)then
                            print *,"(t1,t2) = ",t1,t2
                            print *,"tdt rng = ",tdt(1,1),tdt(1,ntheta)
                        endif
                    endif
                    ! bounce locations recorded for optional output
                    turns%fs(ilmda-1,1) = t1
                    CALL bicube_eval(rzphi,psi,t1,0)
                    turns%fs(ilmda-1,2)=ro+SQRT(rzphi%f(1))*COS(twopi*(t1+rzphi%f(2)))
                    turns%fs(ilmda-1,3)=zo+SQRT(rzphi%f(1))*SIN(twopi*(t1+rzphi%f(2)))
                    turns%fs(ilmda-1,4) = t2
                    CALL bicube_eval(rzphi,psi,t2,0)
                    turns%fs(ilmda-1,5)=ro+SQRT(rzphi%f(1))*COS(twopi*(t2+rzphi%f(2)))
                    turns%fs(ilmda-1,6)=zo+SQRT(rzphi%f(1))*SIN(twopi*(t2+rzphi%f(2)))
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
                                    print('(2x,a61,es10.3E2)'), "WARNING: vpar zero crossing &
                                            &internal to magnetic well at psi ",psi
                                    print('(2x,es10.3E2,a4,es10.3E2,a4,es10.3E2)'),&
                                        t1," <= ",tdt(1,i)," <= ",t2
                                    print *, "  -> Lambda, t/p boundry = ",lmda,lmdatpb
                                    print *, "  -> consider changing mtheta in equil.in"
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
                        dbob = sum(dbob_m%f(:)*expm)
                        divx = sum(divx_m%f(:)*expm) * divxfac
                        jvtheta(i) = tdt(2,i)*tspl%f(4)*tspl%f(1) &
                            *(divx*sqrt(vpar)+dbob*(1-1.5*lmda*tspl%f(1)/bo)/sqrt(vpar))&
                            *exp(-twopi*xj*n*q*(tdt(1,i)-tdt(1,1))) 
                            ! theta0 doesn't really matter since |dj|^2
                        ! debuging
                        if(jvtheta(i)/=jvtheta(i))then
                            print *,'itheta,ntheta,theta',i,ntheta,tdt(1,i)
                            print *,'psi ',psi
                            print *,'dbob_m = ',dbob_m%f(:)
                            print *,'divx_m = ',divx_m%f(:)
                            print *,jvtheta(i)
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
                    bhat = sqrt(2*kin%f(s+2)/mass)/ro
                    dhat = (kin%f(s+2)/chrg)/(bo*ro*ro)
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
                    fbnce%fs(ilmda-1,3) = wbbar*abs(bjspl%fsi(bjspl%mx,1))**2/ro**2
                    
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
                    
                
                
                    ! optional deeply trapped bounce motion output
                    if(present(op_tfuns) .and. (ilmda==1 .or. ilmda==nlmda/2 .or.ilmda==nlmda))then
                        if(tdebug) print *, "  recording bounce functions"
                        j = (ilmda*2)/nlmda ! 0,1,2
                        do i=1,ntheta
                            expm = exp(xj*twopi*mfac*tdt(1,i))
                            dbob = sum(dbob_m%f(:)*expm)
                            divx = sum(divx_m%f(:)*expm) * divxfac
                            op_tfuns(:,j*ntheta+i) = (/ &
                                real(i,r8), real(j+1,r8), real(l,r8), psi, lmda, &
                                bjspl%xs(i-1), tdt(:,i), bspl%fs(i-1,:), &
                                real(bjspl%fs(i-1,1)), aimag((bjspl%fs(i-1,1))), &
                                bspl%fsi(i-1,1)/((2-sigma)*bspl%fsi(bspl%mx,1)), &
                                bspl%fsi(i-1,2)/((2-sigma)*bspl%fsi(bspl%mx,2)), &
                                real(dbob), aimag(dbob), real(divx), aimag(divx) /)
                        enddo
                    endif
                
                
                enddo ! lambda loop
                ! normalize lambda functions for lsode (note: no resonant operator yet)
                allocate(fbnce_norm(fbnce%nqty-2))
                do i=1,fbnce%nqty-2
                    fbnce_norm(i) = 1.0/median(abs(fbnce%fs(:,i+2)))!maxval(abs(fbnce%fs(:,i+2)),1)!median(fbnce%fs(:,i+2))!
                    fbnce%fs(:,i+2) = fbnce%fs(:,i+2)*fbnce_norm(i)
                enddo
                CALL cspline_fit(fbnce,'extrap')
                CALL spline_fit(turns,'extrap')
                
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
                wtwnorm = 1.0 ! euler-lagrange real energy metrices (default)
                rex = 1.0 ! include real part of resonance operator (default)
                imx = 1.0 ! include imag part of resonance operator (default)
                if(method(2:4)=='wmm' .or. method(2:4)=='kmm')then
                    rex = 0.0
                endif
                if(method(2:4)=='tmm' .or. method(2:4)=='rmm')then
                    imx = 0.0
                    wtwnorm = -1.0
                endif
                if(erecord)then
                    lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                        epsr,q,fbnce,l,n,op_rex=rex,op_imx=imx,&
                        op_psi=psi,op_turns=turns)
                else
                    lxint = lambdaintgrl_lsode(wdian,wdiat,welec,nuk,bo/bmax,&
                        epsr,q,fbnce,l,n,op_rex=rex,op_imx=imx)
                endif
                
                ! dT/dpsi
                tpsi = (-2*n**2/sqrt(pi))*(ro/bo)*kin%f(s)*kin%f(s+2) &
                    *lxint(1)/fbnce_norm(1) &       ! lsode normalization
                    *(chi1/twopi) ! unit conversion from psi to psi_n, theta_n to theta
                if(tdebug) print *,'  ->  lxint',lxint(1),', tpsi ',tpsi

                ! Euler lagrange matrices (wtw ~ dW ~ T/i)
                if(present(op_wmats))then
                    op_wmats(:,:,:) = 0
                    do k=1,6
                        do j=1,mpert
                            do i=1,mpert
                                iqty = ((k-1)*mpert+j-1)*mpert + i + 1
                                op_wmats(i,j,k) = (1/(2*xj*n))*lxint(iqty)/fbnce_norm(iqty) &
                                    *kin%f(s)*kin%f(s+2) &
                                    *(-n**2/sqrt(pi))*(ro/bo)*(chi1/twopi)
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
                        tpsi = 0
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
                            CALL zgesvd('S','S',mpert,mpert,a,mpert,svals,u, &
                                mpert,vt,mpert,work,lwork,rwork,info)
                            if(tdebug) print *,i,maxval(svals)
                            tpsi = tpsi + maxval(svals) ! euclidean (spectral) norm
                        enddo
                        tpsi = (rex+imx*xj)*sqrt(tpsi) ! euclidean norm of the 6 norms
                    elseif(index(method,'mm')>0)then ! Mode-coupled dW of T_phi 
                        call cspline_eval(xs_m(1),psi,0)
                        call cspline_eval(xs_m(2),psi,0)
                        call cspline_eval(xs_m(3),psi,0)
                        xix(:,1) = xs_m(1)%f(:)
                        xiy(:,1) = xs_m(2)%f(:)
                        xiz(:,1) = xs_m(3)%f(:)
                        t_zz = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,1),xiz))*chi1**2
                        t_zx = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,2),xix))*chi1
                        t_zy = matmul(conjg(transpose(xiz)),matmul(op_wmats(:,:,3),xiy))*chi1
                        t_xx = matmul(conjg(transpose(xix)),matmul(op_wmats(:,:,4),xix))
                        t_xy = matmul(conjg(transpose(xix)),matmul(op_wmats(:,:,5),xiy))
                        t_yy = matmul(conjg(transpose(xiy)),matmul(op_wmats(:,:,6),xiy))
                        tpsi = (2*n*xj/(2*mu0))*(t_zz(1,1)+t_xx(1,1)+t_yy(1,1) &
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
    end function tpsi

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
    !       Ion charge in fundemental units (e).
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
        integer :: i, j, l, s, mx, istrt = -1, istop=-1
        real(r8) :: x,xlast,chrg,drive,wdcom,dxcom
        real(r8), dimension(nfluxfuns) :: fcom
        real(r8), dimension(:), allocatable :: xs
        real(r8), dimension(:,:), allocatable :: profiles,ellprofiles
        complex(r8), dimension(:), allocatable :: gam,chi
        character(8) :: methcom
        ! lsode type variables
        integer  neqarray(6),neq
        real*8, dimension(:), allocatable ::  y,dky
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
            chrg = e
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
        psilim(1) = max(psilim(1),xs(0))
        psilim(2) = min(psilim(2),xs(mx))
        do i=0,mx
            if(istrt<0 .and. xs(i)>=psilim(1)) exit
        enddo
        istrt = i
        do i=mx,0,-1
            if(sq%xs(i)<=xs(mx)) exit
        enddo
        istop = i
        mx = istop-istrt

        ! prep allocations
        allocate(profiles(mx+1,10+nfluxfuns),ellprofiles(neq*(mx+1),10))
        CALL cspline_alloc(tphi_spl,mx,2*nl+1)
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
            CALL tintgrnd(neqarray,x,y,dky)
            tphi_spl%fs(i,:) = dky(1:neq:2)+dky(2:neq:2)*xj
            tphi_spl%xs(i) = x

            ! save matrix of coefficients
            if(index(method,'mm')>0)then
                if(tdebug) print *, "Euler-Lagrange tmp vars, index=",i
                do j=1,6
                    kelmm(j)%xs(i) = x
                    kelmm(j)%fs(i,:) = RESHAPE(elems(:,:,j),(/mpert**2/))
                enddo
            endif

            ! log progress
            if(mod(x,.1_r8)<mod(xlast,.1_r8) .or. xlast==xs(0) .and. verbose)then
                print('(a7,f4.1,a13,2(es10.2e2),a1)'), " psi =",x,&
                    " -> dT/dpsi= ",sum(dky(1:neq:2),dim=1),sum(dky(2:neq:2),dim=1),'j'
            endif
            xlast = x
        enddo
        ! fit and integrate torque density spline
        CALL cspline_fit(tphi_spl,"extrap")
        CALL cspline_int(tphi_spl)

        ! flux/diffusivity profiles
        do i=0,mx
            x = xs(i+istrt)
            dky = tphi_spl%fs(i,:)
            y = tphi_spl%fsi(i,:)
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
            profiles(i+1,:) = (/ x, sq%f(3), &
                ri(sum(gam(:),dim=1)), &
                ri(sum(chi(:),dim=1)), &
                ri(sum(tphi_spl%fs(i,:),dim=1)), &
                ri(sum(tphi_spl%fsi(i,:),dim=1)), &
                fcom /)
            do l=-nl,nl
                ellprofiles(i*neq + (nl+l) + 1, :) = (/ x, l*1.0_r8, &
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
    !       Ion charge in fundemental units (e).
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
        integer :: i,j,l,s
        real(r8) :: xlast,wdcom,dxcom,chrg,drive
        real(r8), dimension(nfluxfuns) :: fcom
        real(r8), dimension(:), allocatable :: gam,chi
        real(r8), dimension(:,:), allocatable :: profiles,ellprofiles
        complex(r8), dimension(maxsteps,mpert**2,6) :: kel_flat_mats
        character(8) :: methcom
        ! declare lsode input variables
        integer  iopt, iout, istate, itask, itol, mf, iflag,neqarray(6),&
            neq,liw,lrw
        integer, dimension(:), allocatable :: iwork
        real*8 :: x,xout
        real*8, dimension(:), allocatable ::  atol,rtol,rwork,y,dky
        
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
        x = sq%xs(0)
        xout = min(xs_m(1)%xs(xs_m(1)%mx),sq%xs(sq%mx)) ! psilim, psihigh
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
        ! enforce integration bounds
        x = max(psilim(1),x)
        xout = min(psilim(2),xout)
        rwork(1) = xout
        if(tdebug) print *,"psilim = ",psilim
        if(tdebug) print *,"sq lim = ",sq%xs(0),sq%xs(sq%mx)
        if(tdebug) print *,"xs lim = ",xs_m(1)%xs(0),xs_m(1)%xs(xs_m(1)%mx)
        if(tdebug) print *,"x,xout = ",x,xout

        ! for flux calculation
        if(electron)then
            chrg = e
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
                    kel_flat_mats(iwork(11),:,j) = RESHAPE(elems(:,:,j),(/mpert**2/))
                enddo
            endif

            ! print progress
            if(mod(x,.1_r8)<mod(xlast,.1_r8) .or. xlast==sq%xs(0) .and. verbose)then
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
        !if(output_netcdf)then
        !    call output_torque_netcdf(n,zi,mi,electron,trim(method),profiles,ellprofiles)
        !endif

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
        real*8 x, y(*), ydot(neq(1))
    
        real(r8) :: wdfac,xfac,psi
        real(r8), dimension(nfluxfuns) :: ffuns
        complex(r8) trq
        integer :: n,l,zi,mi,ee,nl
        logical :: electron=.false.,first=.true.
        character(8) :: method
        
        integer :: i,j
        !complex(r8), dimension (mpert,mpert,6,-neq(3):neq(3)) :: wtw_l
        complex(r8), dimension (mpert,mpert,6) :: wtw_l
                
        common /tcom/ wdfac,xfac,method,ffuns
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
        do l=-nl,nl
            if(l==0)then
                if(index(method,'mm')>0)then
                    trq = tpsi(psi,n,l,zi,mi,wdfac,xfac,electron,method,&
                               op_ffuns=ffuns,op_wmats=wtw_l)
                    elems = elems+wtw_l
                else
                    trq = tpsi(psi,n,l,zi,mi,wdfac,xfac,electron,method,&
                               op_ffuns=ffuns)
                endif
            else
                if(index(method,'mm')>0)then
                    trq = tpsi(psi,n,l,zi,mi,wdfac,xfac,electron,method,op_wmats=wtw_l)
                    elems = elems+wtw_l
                else
                    trq = tpsi(psi,n,l,zi,mi,wdfac,xfac,electron,method)
                endif
            endif
            ! decouple two real space solutions
            ydot(2*(l+nl)+1) = real(trq)
            ydot(2*(l+nl)+2) = aimag(trq)
        enddo
        
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
        real*8  t, y, pd(nrpd,2)
        ! null result
        pd(:,:) = 0
        return
    end subroutine noj


    !=======================================================================
    subroutine ascii_table(file,list,labels,op_header,op_zi,op_mi)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii files.
    !
    !*ARGUMENTS:
    !    list : real 2D
    !       Table of values
    !    labels : character(17) 1D
    !       Column labels
    !
    !*OPTIONAL ARGUMENTS:
    !    header : character(5)
    !       "equil" -> Write flux function definition header
    !       "pentq" -> Write non-ambipolar heat transport header
    !       "pentt" -> Write non-ambipolar momentum transport header
    !       "theta" -> Write poloidal function definition header
    !
    !-----------------------------------------------------------------------
        implicit none
        real(r8), dimension(:,:), intent(in) :: list
        character(17), dimension(:), intent(in) :: labels
        character(128), intent(in) :: file
        character(5), optional, intent(in) :: op_header
        integer, optional :: op_zi, op_mi

        integer :: i,nrow,ncol,out_unit
        character(5) :: header = '     '
        character(32) :: label_fmt, table_fmt

        ! learn sizes
        ncol = size(list,dim=1)
        nrow = size(list,dim=2)

        ! open file
        out_unit = get_free_file_unit(-1)
        open(unit=out_unit,file=trim(file),status="unknown",action="write")
        ! default header
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"

        ! optional headers
        if (present(op_header)) header = op_header
        if(header=='equil')then
            write(out_unit,*) " Equilibrium profiles."
            write(out_unit,*) " - All energy dependent quantities are taken at E/T unity."
        elseif(header=='theta')then
            write(out_unit,*) " Bounce average functions"
        elseif(header=='pentq')then
            write(out_unit,*) " Gamma = Q/T = Heat flux in 1/(sm^2)"
            write(out_unit,*) " chi = Heat diffusivity in m^2/s"
            write(out_unit,*) " T_phi = Re[A*e*psi'*Gamma/2pi] in Nm (NOT NTV TORQUE!)"
            write(out_unit,*) " 2ndeltaW = Im[A*e*psi'*Gamma/2pi] in J (NOT dWk!)"
        elseif(header=='pentt')then
            write(out_unit,*) " Gamma = Particle flux in 1/(sm^2)"
            write(out_unit,*) " chi = Particle diffusivity in m^2/s"
            write(out_unit,*) " T_phi = Torque in Nm"
            write(out_unit,*) " 2ndeltaW = Kinetic energy in J"
        endif

        write(out_unit,*)
        if(present(op_zi)) write(out_unit,'((a10,es17.8E3))') "     Ze = ",op_zi*e
        if(present(op_mi)) write(out_unit,'((a10,es17.8E3))') "   mass = ",op_mi*mp
        write(out_unit,'(3(a10,es16.8E3))') "     R0 = ",ro,"     B0 = ",bo,"   chi1 = ",chi1

        ! write labels and table
        write(label_fmt,*) '(1/,',ncol,'(a17))'
        write(table_fmt,*) '(',ncol,'(es17.8E3))'
        write(out_unit,trim(label_fmt)) labels
        do i=1,nrow
            write(out_unit,trim(table_fmt)) list(:,i)
        enddo

        close(out_unit)
        return
    end subroutine ascii_table

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
        
        ! learn sizes
        ncol = size(prof,dim=1)
        nrow = size(prof,dim=2)
        ncoll = size(profl,dim=1)
        nrowl = size(profl,dim=2)
        np = nrow
        nl = nrowl/np
        if(nl*np/=nrowl) stop "ERROR: Unable to reshape torque output into P-by-L"

        select case (method)
            case ("fcgl")
                select case (gridtype)
                    case ('lsode')
                        allocate(fcgl_lsode_psi_record(ncol,nrow),fcgl_lsode_ell_record(ncoll,np,nl))
                        fcgl_lsode_psi_record = prof
                        fcgl_lsode_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('equil')
                        allocate(fcgl_equil_psi_record(ncol,nrow),fcgl_equil_ell_record(ncoll,np,nl))
                        fcgl_equil_psi_record = prof
                        fcgl_equil_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('input')
                        allocate(fcgl_input_psi_record(ncol,nrow),fcgl_input_ell_record(ncoll,np,nl))
                        fcgl_input_psi_record = prof
                        fcgl_input_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case default
                        Stop "ERROR: fcgl grid type does not have a global record"
                end select
            case ("rlar")
                select case (gridtype)
                    case ('lsode')
                        allocate(rlar_lsode_psi_record(ncol,nrow),rlar_lsode_ell_record(ncoll,np,nl))
                        rlar_lsode_psi_record = prof
                        rlar_lsode_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('equil')
                        allocate(rlar_equil_psi_record(ncol,nrow),rlar_equil_ell_record(ncoll,np,nl))
                        rlar_equil_psi_record = prof
                        rlar_equil_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('input')
                        allocate(rlar_input_psi_record(ncol,nrow),rlar_input_ell_record(ncoll,np,nl))
                        rlar_input_psi_record = prof
                        rlar_input_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case default
                        Stop "ERROR: rlar grid type does not have a global record"
                end select
            case ("clar")
                select case (gridtype)
                    case ('lsode')
                        allocate(clar_lsode_psi_record(ncol,nrow),clar_lsode_ell_record(ncoll,np,nl))
                        clar_lsode_psi_record = prof
                        clar_lsode_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('equil')
                        allocate(clar_equil_psi_record(ncol,nrow),clar_equil_ell_record(ncoll,np,nl))
                        clar_equil_psi_record = prof
                        clar_equil_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('input')
                        allocate(clar_input_psi_record(ncol,nrow),clar_input_ell_record(ncoll,np,nl))
                        clar_input_psi_record = prof
                        clar_input_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case default
                        Stop "ERROR: clar grid type does not have a global record"
                end select
            case ("tgar")
                select case (gridtype)
                    case ('lsode')
                        allocate(tgar_lsode_psi_record(ncol,nrow),tgar_lsode_ell_record(ncoll,np,nl))
                        tgar_lsode_psi_record = prof
                        tgar_lsode_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('equil')
                        allocate(tgar_equil_psi_record(ncol,nrow),tgar_equil_ell_record(ncoll,np,nl))
                        tgar_equil_psi_record = prof
                        tgar_equil_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('input')
                        allocate(tgar_input_psi_record(ncol,nrow),tgar_input_ell_record(ncoll,np,nl))
                        tgar_input_psi_record = prof
                        tgar_input_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case default
                        Stop "ERROR: tgar grid type does not have a global record"
                end select
            case ("fgar")
                select case (gridtype)
                    case ('lsode')
                        allocate(fgar_lsode_psi_record(ncol,nrow),fgar_lsode_ell_record(ncoll,np,nl))
                        fgar_lsode_psi_record = prof
                        fgar_lsode_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('equil')
                        allocate(fgar_equil_psi_record(ncol,nrow),fgar_equil_ell_record(ncoll,np,nl))
                        fgar_equil_psi_record = prof
                        fgar_equil_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case ('input')
                        allocate(fgar_input_psi_record(ncol,nrow),fgar_input_ell_record(ncoll,np,nl))
                        fgar_input_psi_record = prof
                        fgar_input_ell_record = RESHAPE(profl,(/ncoll,nl,np/))
                    case default
                        Stop "ERROR: fgar grid type does not have a global record"
                end select
            case default
                Stop "ERROR: Requested type does not have a global record"
        end select

        ! keep a key of what has been recorded
        do i=1,nmethods
            if (methods(i)==method)then
                do j=1,ngrids
                    if(grids(j)==gridtype) isrecorded(i,j) = .true.
                enddo
            endif
        enddo

        return
    end subroutine record_method

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

        integer :: i,nrow,ncol,nrowl,ncoll,out_unit,unit1,unit2
        real(r8) :: chrg,mass
        character(8) :: nstring
        character(32) :: fmt
        character(128) :: file1,file2

        ! set species
        if(electron)then
            chrg = e
            mass = me
        else
            chrg = zi*e
            mass = mi*mp
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
            "<(deltaB_L/B)^2>","<(divxprp)^2>"!/J??

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
    subroutine output_torque_netcdf(n,nl,zi,mi,electron)
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
        integer, intent(in) :: n,nl,zi,mi
        logical, intent(in) :: electron

        integer :: i,j,nrow,ncol,nrowl,ncoll,np
        real(r8) :: chrg,mass
        real(r8), dimension(:,:), pointer :: psi_record
        real(r8), dimension(:,:,:), pointer :: ell_record

        integer :: ncid,i_did,i_id,p_did,p_id,l_did,l_id, &
            v_id,g_id,c_id,d_id,t_id
        character(16) :: nstring,suffix,method,gridtype
        character(128) :: ncfile
        
        print *,"Writing output to netcdf"

        ! set species
        if(electron)then
            chrg = e
            mass = me
        else
            chrg = zi*e
            mass = mi*mp
        endif

        ! create and open netcdf file
        write(nstring,'(I8)') n
        ncfile = "pentrc_output_n"//TRIM(ADJUSTL(nstring))//".nc"
        CALL check( nf90_create(ncfile, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
        ! store attributes
        CALL check( nf90_put_att(ncid, nf90_global, "title", "PENTRC fundamental outputs") )
        CALL check( nf90_put_att(ncid, nf90_global, "version", version))
        CALL check( nf90_put_att(ncid, nf90_global, "shot", INT(shotnum) ) )
        CALL check( nf90_put_att(ncid, nf90_global, "time", INT(shottime)) )
        CALL check( nf90_put_att(ncid, nf90_global, "machine", machine) )
        CALL check( nf90_put_att(ncid, nf90_global, "n", n) )
        CALL check( nf90_put_att(ncid, nf90_global, "charge", chrg) )
        CALL check( nf90_put_att(ncid, nf90_global, "mass", mass) )
        CALL check( nf90_put_att(ncid, nf90_global, "R0", ro) )
        CALL check( nf90_put_att(ncid, nf90_global, "B0", bo) )
        CALL check( nf90_put_att(ncid, nf90_global, "chi1", chi1) )
        ! define dimensions
        CALL check( nf90_def_dim(ncid, "i", 2, i_did) )
        CALL check( nf90_def_var(ncid, "i", nf90_int, i_did, i_id) )
        CALL check( nf90_def_dim(ncid, "ell", nl, l_did) )
        CALL check( nf90_def_var(ncid, "ell", nf90_int, l_did, l_id) )
        ! End definitions
        CALL check( nf90_enddef(ncid) )
        ! store variables
        CALL check( nf90_put_var(ncid, i_id, (/0,1/)) )
        CALL check( nf90_put_var(ncid, l_id, (/(i,i=-nl,nl)/)) )

        do i=1,nmethods
            method = methods(i)
            do j=1,ngrids
                gridtype = grids(j)
                select case (method)
                    case ("fcgl")
                        select case (gridtype)
                            case ('lsode')
                                psi_record => fcgl_lsode_psi_record
                                ell_record => fcgl_lsode_ell_record
                            case ('equil')
                                psi_record => fcgl_equil_psi_record
                                ell_record => fcgl_equil_ell_record
                            case ('input')
                                psi_record => fcgl_input_psi_record
                                ell_record => fcgl_input_ell_record
                            case default
                                Stop "ERROR: fcgl grid type does not have a global record"
                        end select
                    case ("rlar")
                        select case (gridtype)
                            case ('lsode')
                                psi_record => rlar_lsode_psi_record
                                ell_record => rlar_lsode_ell_record
                            case ('equil')
                                psi_record => rlar_equil_psi_record
                                ell_record => rlar_equil_ell_record
                            case ('input')
                                psi_record => rlar_input_psi_record
                                ell_record => rlar_input_ell_record
                            case default
                                Stop "ERROR: rlar grid type does not have a global record"
                        end select
                    case ("clar")
                        select case (gridtype)
                            case ('lsode')
                                psi_record => clar_lsode_psi_record
                                ell_record => clar_lsode_ell_record
                            case ('equil')
                                psi_record => clar_equil_psi_record
                                ell_record => clar_equil_ell_record
                            case ('input')
                                psi_record => clar_input_psi_record
                                ell_record => clar_input_ell_record
                            case default
                                Stop "ERROR: clar grid type does not have a global record"
                        end select
                    case ("tgar")
                        select case (gridtype)
                            case ('lsode')
                                psi_record => tgar_lsode_psi_record
                                ell_record => tgar_lsode_ell_record
                            case ('equil')
                                psi_record => tgar_equil_psi_record
                                ell_record => tgar_equil_ell_record
                            case ('input')
                                psi_record => tgar_input_psi_record
                                ell_record => tgar_input_ell_record
                            case default
                                Stop "ERROR: tgar grid type does not have a global record"
                        end select
                    case ("fgar")
                        select case (gridtype)
                            case ('lsode')
                                psi_record => tgar_lsode_psi_record
                                ell_record => fgar_lsode_ell_record
                            case ('equil')
                                psi_record => fgar_equil_psi_record
                                ell_record => fgar_equil_ell_record
                            case ('input')
                                psi_record => fgar_input_psi_record
                                ell_record => fgar_input_ell_record
                            case default
                                Stop "ERROR: fgar grid type does not have a global record"
                        end select
                    case default
                        Stop "ERROR: Requested type does not have a global record"
                end select
                print *,i,j
                print *,method,gridtype
                print *,isrecorded(i,j)
                ! If this method and gridtype has been recorded, write it to the file
                if(isrecorded(i,j))then
                    ! learn sizes
                    np = size(psi_record,dim=2)
                    if(np/=size(ell_record,dim=3)) stop "Error: Record sizes are inconsistent"
                    ! create distinguishing labels
                    if(gridtype=='lsode')then
                        suffix = '_'//trim(method)
                    else
                        suffix = '_'//trim(method)//'_'//trim(gridtype)
                    endif
                    ! Re-open definitions
                    CALL check( nf90_redef(ncid) )
                    CALL check( nf90_put_att(ncid, nf90_global, "T_total"//trim(suffix), psi_record(9,np)))
                    CALL check( nf90_put_att(ncid, nf90_global, "dW_total"//trim(suffix), psi_record(10,np)/(2*n)))
                    CALL check( nf90_def_dim(ncid, "psi"//trim(suffix), np, p_did) )
                    CALL check( nf90_def_var(ncid, "psi"//trim(suffix), nf90_double, p_did, p_id) )
                    CALL check( nf90_def_var(ncid, "Gamma"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), g_id) )
                    CALL check( nf90_def_var(ncid, "chi"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), c_id) )
                    CALL check( nf90_def_var(ncid, "t_phi"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), d_id) )
                    CALL check( nf90_def_var(ncid, "T"//trim(suffix), nf90_double, (/i_did,l_did,p_did/), t_id) )
                    CALL check( nf90_def_var(ncid, "dvdpsi"//trim(suffix), nf90_double, (/i_did,p_did/), v_id) )
                    ! Add information
                    if(qt)then
                        CALL check( nf90_put_att(ncid, g_id, "long_name", "Nonambipolar Heat Flux") )
                        CALL check( nf90_put_att(ncid, g_id, "units", "1/sm^2") )
                        CALL check( nf90_put_att(ncid, c_id, "long_name", "Nonambipolar Heat Diffusivity") )
                        CALL check( nf90_put_att(ncid, c_id, "units", "m^2/s") )
                        CALL check( nf90_put_att(ncid, d_id, "long_name", " A*e*psi'*Gamma/2pi profile") )
                        CALL check( nf90_put_att(ncid, d_id, "units", "Nm per psi") )
                        CALL check( nf90_put_att(ncid, t_id, "long_name", "Integrated A*e*psi'*Gamma/2pi") )
                        CALL check( nf90_put_att(ncid, t_id, "units", "Nm") )
                        CALL check( nf90_put_att(ncid, v_id, "long_name", "Differential Volume") )
                        CALL check( nf90_put_att(ncid, v_id, "units", "m^3") )
                    else
                        CALL check( nf90_put_att(ncid, g_id, "long_name", "Nonambipolar Particle Flux") )
                        CALL check( nf90_put_att(ncid, g_id, "units", "1/sm^2") )
                        CALL check( nf90_put_att(ncid, c_id, "long_name", "Nonambipolar Particle Diffusivity") )
                        CALL check( nf90_put_att(ncid, c_id, "units", "m^2/s") )
                        CALL check( nf90_put_att(ncid, d_id, "long_name", "Toroidal Torque") )
                        CALL check( nf90_put_att(ncid, d_id, "units", "Nm per psi") )
                        CALL check( nf90_put_att(ncid, t_id, "long_name", "Integrated Toroidal Torque") )
                        CALL check( nf90_put_att(ncid, t_id, "units", "Nm") )
                        CALL check( nf90_put_att(ncid, v_id, "long_name", "Differential Volume") )
                        CALL check( nf90_put_att(ncid, v_id, "units", "m^3") )
                    endif
                    CALL check( nf90_enddef(ncid) )
                    ! Put in variables
                    CALL check( nf90_put_var(ncid, p_id, ell_record(1,1,:)) )
                    CALL check( nf90_put_var(ncid, g_id, ell_record(3:4,:,:)) )
                    CALL check( nf90_put_var(ncid, c_id, ell_record(5:6,:,:)) )
                    CALL check( nf90_put_var(ncid, d_id, ell_record(7:8,:,:)) )
                    CALL check( nf90_put_var(ncid, t_id, ell_record(9:10,:,:)) )
                    CALL check( nf90_put_var(ncid, v_id, psi_record(2,:)) )
                endif
            enddo
        enddo
        
        ! close file
        CALL check( nf90_close(ncid) )
        
            


        return
    end subroutine output_torque_netcdf

    !=======================================================================
    subroutine output_fluxfun_ascii(n,zi,mi,electron,method,table)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii flux function profile files.
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
    !    table : real 2D
    !       Table of values writen to file
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n, zi, mi
        character(*), intent(in) :: method
        logical :: electron
        real(r8), dimension(:,:), intent(in) :: table

        integer :: i,out_unit
        character(32) :: nstring,table_fmt,label_fmt
        character(128) :: file

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstring,'(I8)') n
        file = "pentrc_"//trim(method)//"_eqprofiles_n"//trim(adjustl(nstring))//".out"
        if(electron) file = file(:7)//"e_"//file(8:)
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Equilibrium profiles."
        write(out_unit,*) " - All energy dependent quantities are taken at E/T unity."
        write(out_unit,'(1/,1(a10,I4))') "n =",n
        write(out_unit,'(2(a10,es17.8E3))') "Ze =",zi*e,"mass =",mi*mp
        write(out_unit,'(3(a10,es17.8E3))') "R0 =",ro,"B0 =",bo,"chi1 =",chi1
        !write(out_unit,'(1/,a26,2(a10,es16.8E3),1/)') " Common normalizations:",&
        !    "    bbar = ",ro/sqrt(2.0/(mi*mp)),"    dbar = ",zi*e*bo*ro**2 ! needs T_is

        ! write column headers
        write(label_fmt,*) '(1/,',nfluxfuns,'(a17))'
        write(out_unit,label_fmt) "psi_n","eps_r","n_i","n_e",&
            "T_i","T_e","omega_E","logLambda","nu_i","nu_e","q","dconPmu_0","dvdpsi_n",&
            "omega_N","omega_T","omega_trans","omega_gyro","RLARomega_b","RLARomega_d",&
            "<(deltaB_L/B)^2>","<(divxprp)^2>"!/J??

        ! write tables
        write(table_fmt,*) '(',nfluxfuns,'(es17.8E3))'
        do i=1,size(table,dim=2)
            write(out_unit,table_fmt) table(:,i)
        enddo

        close(out_unit)
        return
    end subroutine output_fluxfun_ascii

    !=======================================================================
    subroutine output_bouncefun_ascii(n,zi,mi,electron,method,table)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii bounce function files.
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
    !    table : real 2D
    !       Table of values writen to file
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n, zi, mi
        character(*), intent(in) :: method
        logical :: electron
        real(r8), dimension(:,:), intent(in) :: table

        integer :: i,out_unit
        character(32) :: nstr,lstr,label_fmt,table_fmt
        character(128) :: file

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstr,'(I8)') n
        file = "pentrc_"//trim(method)//"_bounce_n"//trim(adjustl(nstr))//".out"
        if(electron) file = file(:7)//"e_"//file(8:)
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Bounce average functions"
        write(out_unit,'(1/,1(a10,I4))') "n =",n
        write(out_unit,'(2(a10,es17.8E3))') "Ze =",zi*e,"mass =",mi*mp
        write(out_unit,'(3(a10,es17.8E3))') "R0 =",ro,"B0 =",bo,"chi1 =",chi1

        ! write column headers
        write(label_fmt,*) '(1/,',nthetafuns,'(a17))'
        write(out_unit,label_fmt) "i_theta","i_Lambda","ell","psi_n", &
            "Lambda","theta_n","theta","dtheta", &
            "omega_b","omega_D","real(deltaJ)","imag(deltaJ)","h_E","h_D", &
            "real(deltaB)","imag(deltaB)","real(divxprp)","imag(divxprp)"

        ! write tables
        write(table_fmt,*) '(',nthetafuns,'(es17.8E3))'
        do i=1,size(table,dim=2)
            write(out_unit,table_fmt) table(:,i)
        enddo

        close(out_unit)
        return
    end subroutine output_bouncefun_ascii

end module torque




