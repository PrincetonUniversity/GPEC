!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module pitch_integration
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Functional form of the pitch integrand as found in
    !   [Logan, Park, et al., Phys. Plasma, 2013] and various integration
    !   methods.
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !   lambdaintgnd             - Resonance operator Eq. (20)
    !   lambdaintgrl_lsode       - Dynamic energy integration
    !
    !*PUBLIC DATA MEMBERS:
    !   latol               - absolute tolerance of lsode
    !   lrtol               - relative tolerance of lsode
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writting. 
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only : r8,xj
    use utilities, only : get_free_file_unit
    use special, only : ellipk
    use cspline_mod, only: cspline_type,cspline_eval
    use spline_mod, only: spline_type,spline_eval,spline_alloc,&
                            spline_dealloc,spline_fit,spline_int
    use bicube_mod, only: bicube_type,bicube_eval
    use energy_integration, only: xintgrl_lsode
    use lsode1_mod
    use inputs, only: rzphi !added for bounce locations
    
    implicit none
    
    private
    public &
        lambdaatol,lambdartol, &                            ! reals
        lambdadebug, &                                      ! logical
        lambdaintgrl_lsode,kappaintgrl,kappaintgnd          ! functions
    
    ! global variables with defaults
    logical :: lambdadebug = .false.
    real(r8) :: &
        lambdaatol = 1e-12, &
        lambdartol = 1e-9
    ! global variables for internal use
    !real(r8) :: wn_g,wt_g,we_g,nuk_g,l_g,n_g,&
    !    bobmax_g,epsr_g
    type(cspline_type) :: eqspl_g
    type(spline_type) :: turns_g
    
    contains

    !=======================================================================
    function lambdaintgrl_lsode(wn,wt,we,nuk,bobmax,epsr,q,eq_spl,l,n,&
                                op_rex,op_imx,op_psi,op_lbl,op_turns)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Dynamic pitch integration using lsode. Module global variabls
    !   lambdaatol, lambdartol are used for tolerances.
    !
    !   Integrand is defined by flux functions (wn,wt,we), mode numbers
    !   (l,n) and equilibrium pitch functions (eq_spl). See below for
    !   details.
    !
    !   Integral limits defined by range of eq_spl%xs.
    !
    !   Pro-tip: To calculate offset rotation use we=-wn-wt.
    !
    !*ARGUMENTS:
    !   wn : real.
    !       density gradient diamagnetic drift frequency
    !   wt : real.
    !       temperature gradient diamagnetic drift frequency
    !   we : real.
    !       electric precession frequency
    !   nuk : real.
    !       effective Krook collision frequency (i.e. /2eps for trapped)
    !   bmaxbo : real.
    !       Trapped passing boundary in Lambda (max[B(theta)|_psi]/B0))
    !   epsr : real.
    !       Inverse aspect ratio. Trapped particles have an effective
    !       collisionaity nuk/2epsr
    !   eq_spl : cspline_type.
    !       Equilibrium pitch angle functions. Must be following order
    !           wb(Lambda) : Bounce frequency
    !           wd(Lambda) : Bounce avaraged magnetic precession /x
    !           f(Lambda)  : Integrand function.
    !   l : real.
    !       effective bounce harmonic (i.e. ell-nq for passing)
    !   n : integer.
    !       toroidal mode number
    !*OPTIONAL ARGUMENTS:
    !   op_rex : real.
    !       Multiplier applied to real energy integral. The component
    !       gives the torque, but should be ignored when building Euler-Lagrange
    !       matrices for kinetic DCON.
    !   op_imx : real.
    !       Multiplier applied to imaginary energy integral. The component
    !       gives the energy, but should be ignored when building a complex
    !       mode-coupling torque matrix.
    !   op_psi : real.
    !       normalized flux. If this variable is used, integrand and
    !       integral are output to a file pentrc_<op_lbl>_lambda_n<n>_l<l>.out.
    !   op_lbl : character(32).
    !       output file modification
    !   op_turns : spline_type.
    !       If included with op_psi the bounce points are recorded.
    !           t1(Lambda) : Lower bounce theta
    !           t2(Lambda) : Upper bounce theta
    !
    !*RETURNS:
    !     complex.
    !        Integral int{dLambda f(Lambda)*int{dx Rln(x,Lambda)}, where
    !       f is input in the eq_spl and Rln is the resonance operator.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        complex(r8), dimension(:), allocatable :: lambdaintgrl_lsode        
        ! declare arguments
        integer, intent(in) :: n,l
        real(r8), intent(in) :: wn,wt,we,nuk,bobmax,epsr,q
        type(cspline_type) eq_spl      
        type(spline_type), optional :: op_turns      
        real(r8), intent(in), optional :: op_psi,op_rex,op_imx
        character(4), optional :: op_lbl
        ! declare variables
        integer :: out_unit,i,l_g,n_g
        real(r8) :: lmda,wb,wd,nueff,lnq,wn_g,wt_g,we_g,nuk_g,&
            bobmax_g,epsr_g,q_g,rex_g,imx_g
        complex(r8) :: xint
        character(16) :: flbl = ''
        character(64) ::file
        character(8) :: nstring,lstring
        logical :: fexists
        ! declare lsode input variables
        integer  iopt, iout, istate, itask, itol, mf, iflag,neqarray(1),&
            neq,liw,lrw
        integer, dimension(:), allocatable :: iwork
        real*8 :: x,xout
        real*8, dimension(:), allocatable ::  atol,rtol,rwork,y,dky
        !integer  iopt, iout, istate, itask, itol, mf, iflag
        !integer, parameter ::   &
        !    neq(1) = 2,         &   ! true number of equations
        !    liw  = 20 + neq(1), &   ! for mf 22 ! only uses 20 if mf 10
        !    lrw  = 20 + 16*neq(1)   ! for mf 10 
        !    !lrw = 22+9*neq(1)+neq(1)**2 ! for mf 22
        !integer iwork(liw)
        !real*8  atol(1),rtol(1),rwork(lrw),x,xout,y(neq(1)),dky(neq(1))
        
        common /lcom/ wn_g,wt_g,we_g,nuk_g,bobmax_g,epsr_g, &
                    q_g,l_g,n_g,rex_g,imx_g
        
        ! set lsode options - see lsode package for documentation
        neq = 2*(eq_spl%nqty-2)
        liw  = 20 + neq         ! for mf 22 ! only uses 20 if mf 10
        lrw  = 20 + 16*neq      ! for mf 10
        allocate(iwork(liw))
        allocate(atol(neq),rtol(neq),rwork(lrw),y(neq),dky(neq))
        neqarray(:) = (/ neq /)
        y(:) = 0
        x = eq_spl%xs(0)
        xout = eq_spl%xs(eq_spl%mx)
        itol = 2                  ! rtol and atol are arrays
        rtol = lambdartol              !1.d-7!9              ! 14
        atol = lambdaatol              !1.d-7!9              ! 15
        istate = 1                ! first step
        iopt = 1                  ! optional inputs
        iwork(:) = 0              ! defaults
        rwork(:) = 0              ! defaults
        rwork(1) = xout           ! only used if itask 4,5
        iwork(6) = 10000          ! max number steps
        mf = 10                   ! not stiff with unknown J 
    
        ! set glbl variables for access in integrand
        wn_g = wn
        wt_g = wt
        we_g = we
        nuk_g = nuk
        bobmax_g = bobmax
        epsr_g = epsr
        q_g = q
        l_g  = l
        n_g = n
        eqspl_g = eq_spl
        if(present(op_turns)) turns_g = op_turns
        if(present(op_rex))then
            rex_g = op_rex
        else
            rex_g = 1.0
        endif
        if(present(op_imx))then
            imx_g = op_imx
        else
            imx_g = 1.0
        endif
        
        if(present(op_psi)) then
            ! open and prepare file as needed
            out_unit = get_free_file_unit(-1)
            write(lstring,'(SPI3.2)') l
            write(nstring,'(I3)') n
            if(present(op_lbl)) flbl = trim(op_lbl)//"_l"//trim(adjustl(lstring))
            file = "pentrc_"//trim(flbl)//"_pitch_n"//trim(adjustl(nstring))//".out"
            inquire(file=trim(file),exist=fexists)
            if(fexists) then
                open(unit=out_unit,file=file,status="old",position="append")
            else
                open(unit=out_unit,file=file,status="new")!,action="write")
                write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
                write(out_unit,*) "Pitch angle integrand and functions"
                write(out_unit,*) "  variables are:   lambda =  B0*m*v_perp^2/(2B),  x = E/T"
                write(out_unit,*) "  frequencies are taken at x unity"
                write(out_unit,*) "n = ",int(n)," l = ",l
                write(out_unit,*)
                if(present(op_turns))then
                    write(out_unit,'(8(1x,a16),1x,a19,6(1x,a16))') "psi_n","Lambda","T_phi", &
                    "2ndeltaW","int(T_phi)","int(2ndeltaW)","omega_b","omega_D",&
                    "deltaJdeltaJomega_b","theta_l","r_l","z_l","theta_u","r_u","z_u"
                else
                    write(out_unit,'(8(1x,a16),1x,a19)') "psi_n","Lambda","T_phi", &
                    "2ndeltaW","int(T_phi)","int(2ndeltaW)","omega_b","omega_D",&
                    "deltaJdeltaJomega_b"
                endif
            endif
            itask = 5              ! single step
            do while (x<xout)
                call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol,&
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                call dintdy1(x, 1, rwork(21), neqarray(1), dky, iflag)
                if(present(op_turns))then
                    call spline_eval(turns_g,x,0)
                    write(out_unit,'(8(1x,es16.8e3),1x,es19.8e3,6(1x,es16.8e3))') &
                        op_psi,x,dky(1:2),y(1:2),real(eqspl_g%f(1:2)),real(eqspl_g%f(eqspl_g%nqty)),&
                        turns_g%f(:)
                else
                    write(out_unit,'(8(1x,es16.8e3),1x,es19.8e3)') &
                        op_psi,x,dky(1:2),y(1:2),real(eqspl_g%f(1:2)),real(eqspl_g%f(eqspl_g%nqty))
                endif
            enddo
            close(out_unit)
            ! write select energy integral output files
            do i = 0,4
                lmda = eq_spl%xs(0) + (i/4.0)*(xout-eq_spl%xs(0))
                call cspline_eval(eq_spl,lmda,0)
                wb = real(eq_spl%f(1))
                wd = real(eq_spl%f(2))
                if(lmda>bobmax)then
                    nueff = nuk/(2*epsr)
                    lnq = 1.0*l
                else
                    nueff = nuk
                    lnq = l+n*q
                endif
                ! note: currently ignores -wb case for trapped
                xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,lnq,n,(/op_psi,lmda/),flbl)
            enddo
        else
            itask = 4              ! full integral
            call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol,atol, &
                itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
        endif
        
        ! write error file and stop program if integration fials
        if(iwork(11)>iwork(6)/2 .and. istate/=-1) then
            print *, "WARNING: ",iwork(11)," of maximum ",iwork(6)," steps in lambda integration"
        endif
        if(istate==-1) then
            out_unit = get_free_file_unit(-1)
            open(unit=out_unit,file="pentrc_lambdaintrgl_lsode.err",status="unknown")
            if(present(op_psi)) then
                write(out_unit,*) "psi = ",op_psi
            endif
            write(out_unit,'(5(1x,a16))') "lambda","t_phi","2ndeltaw","int(t_phi)","int(2ndeltaw)"
            itask = 2
            y(1:2) = (/ 0,0 /)
            x = eq_spl%xs(0)       ! lower bound
            istate = 1             ! first step
            iwork(:) = 0           ! defaults
            rwork(:) = 0           ! defaults
            rwork(1) = xout        ! only used if itask 4,5
            iwork(6) = 10000       ! max number steps
            do while (x<xout .and. iwork(11)<iwork(6))
                call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol, &
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
               call dintdy1(x, 1, rwork(21), neqarray(1), dky, iflag)
               write(out_unit,'(1x,5(1x,es16.8e3))') x,dky(1:2),y(1:2)
            enddo
            close(out_unit)
            stop "ERROR: lambdaintgrl_lsode - too many steps in lambda required."
        endif

        ! convert to complex space if integrations successful
        if(lambdadebug) print *,"Finished Lambda lsode. Converting to complex..."
        if(allocated(lambdaintgrl_lsode)) deallocate(lambdaintgrl_lsode)
        allocate(lambdaintgrl_lsode(neq/2))
        if(lambdadebug) print *,"...allocated lambdaintgrl",neq/2,"..."        
        lambdaintgrl_lsode(:) = y(1::2)+xj*y(2::2)!(/(y(2*i-1)+xj*y(2*i),i=1,neq/2)/)

        deallocate(iwork,atol,rtol,rwork,y,dky)
        if(lambdadebug) print *,"Returning."        
        
        return
    end function lambdaintgrl_lsode

    
    
    
    
    

    !=======================================================================
    subroutine lintgrnd(neq,x,y,ydot)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Dynamic pitch integration using lsode. Module global variabls
    !   lambdaatol, lambdartol are used for tolerances.
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
    !*OPTIONAL ARGUMENTS: *** CURRENTLY UNAVAILABLE DUE TO LSODE ERROR ***
    !    wntedbk_ln : real, dimension(8) (in)
    !       wn : real.
    !           density gradient diamagnetic drift frequency
    !       wt : real.
    !           temperature gradient diamagnetic drift frequency
    !       we : real.
    !           electric precession frequency
    !       nuk : real.
    !           effective Krook collision frequency (i.e. /2eps for trapped)
    !       eq_spl : cspline_type.
    !           Equilibrium lambda functions
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
        integer ::  neq
        real*8 x, y(neq), ydot(neq)
    
        integer :: l,n,i
        real(r8) :: wn,wt,we,wd,wb,nuk,bobmax,epsr,lnq,q,fl,nueff,rex,imx
        complex(r8) :: xint,fres
    
        common /lcom/ wn,wt,we,nuk,bobmax,epsr,q,l,n,rex,imx
    
        ! use (input or) global variables
        if(lambdadebug) print *,'lintgrnd - lambda =',x
        call cspline_eval(eqspl_g,x,0)
        wb = real(eqspl_g%f(1))
        wd = real(eqspl_g%f(2))
        fl = real(eqspl_g%f(3))
        
        if(lambdadebug)then
            print *,'lambda,bobmax = ',x,bobmax
            print *,'nuk = ',nuk
            print *,'omegas = ',wn,wt,we,wd,wb
            print *,'f(lmda)= ',eqspl_g%f(3), &
                minval(abs(eqspl_g%f(4:))),maxval(abs(eqspl_g%f(4:)))
            print *,'n,l    = ',n,l
       endif
        
        ! energy integration of resonance opterator
        if(x<=bobmax)then      ! circulating particles
            nueff = nuk
            lnq = l+n*q
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,lnq,n)
        else                        ! trapped particles
            nueff = nuk/(2*epsr)    ! effective collisionality
            lnq = 1.0*l
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,lnq,n)
            xint = xint + xintgrl_lsode(wn,wt,we,wd,-wb,nueff,lnq,n)
        endif
        
        if(lambdadebug) print *,"lnq, nueff = ",lnq,nueff
        
        ! form full lambda function and decouple two real space solutions
        do i=3,eqspl_g%nqty
            fres = eqspl_g%f(i)*(rex*real(xint)+xj*imx*aimag(xint))
            ydot(1+2*(i-3)) = real(fres)
            ydot(2+2*(i-3)) = aimag(fres)
            !if(mod(i,100)==0 .or. i==3 .or. i==eqspl_g%nqty) print *,size(ydot),2+2*(i-3),aimag(fres)
        enddo        
        return
    end subroutine lintgrnd

    
    
    
    
    
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
    function kappaintgrl(n,l,q,marray,dbarray,fnml)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Kappa integral from Eq. (14) [Park, Phys. Rev. Let. 2009].
    !
    !*ARGUMENTS:
    !    n : integer (in)
    !       Toroidal mode.
    !    l : integer (in)
    !       Bounce harmonic.
    !    q : real (in)
    !       Safey factor.
    !    marray : reals
    !       Poloidal mode numbers.
    !    dbarray : real.
    !       Corresponding deltaB/B values.
    !    fnml : bicube_type.
    !       Spline of special function from [Park PRL 2009]
    !
    !*RETURNS:
    !     real.
    !        Kappa integral.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        real(r8) :: kappaintgrl
        ! declare arguments
        integer, intent(in) :: n,l
        real(r8), intent(in) :: q
        integer, dimension(:), intent(in) :: marray
        complex(r8), dimension(:) :: dbarray
        real(r8), dimension(:,:), allocatable :: fm
        type(bicube_type) :: fnml
        ! declare local variables
        integer, parameter :: nk = 50
        integer :: i,j,k,mstart,mpert
        real(r8) :: kappa
        type(spline_type) :: kspl
        
        mpert = size(marray,1)
        mstart = lbound(marray,1)
        allocate(fm(mpert,2))
        call spline_alloc(kspl,nk,1)
        kspl%fs(:,1) = 0
        do k=0,nk
            kappa = (k+1.0)/(nk+2.0) !0< kappa<1 !1 is singular
            kspl%xs(k) = kappa
            do i = mstart,mstart+mpert-1
               ! don't let things go beyond precalculated spline
               if(abs(marray(i)-n*q-l) > fnml%ys(fnml%my))then 
                  print *,marray(i),n*q,l
                  stop "ERROR: pitch - RLAR approximation m-nq-l out of Fkmnql range"
               endif
               ! calculate for forward and backward banana
               call bicube_eval(fnml,kappa,marray(i)-n*q-l,0)
               fm(i,1) = fnml%f(1)
               call bicube_eval(fnml,kappa,marray(i)-n*q+l,0)
               fm(i,2) = fnml%f(1)
            enddo
            do i = 1,mpert
               do j = 1,mpert
                  kspl%fs(k,1) = kspl%fs(k,1)&
                    +2*kappa*0.25&
                    *(fm(i,1)*fm(j,1)&
                    +fm(i,2)*fm(j,1)&
                    +fm(i,1)*fm(j,2)&
                    +fm(i,2)*fm(j,2))&
                    *real(dbarray(i)*conjg(dbarray(j)))&
                    /(4*ellipk(kappa))
               enddo
            enddo
        enddo
        call spline_fit(kspl,'extrap')
        call spline_int(kspl)
        kappaintgrl = kspl%fsi(nk,1)
        
        if(lambdadebug)then
            print *,"Kappa integration - marray = ",marray(:5),"...",marray(mpert-5:)
            print *,"Kappa integration - dbarray= ",abs(dbarray(:5)),"...",abs(dbarray(mpert-5:))
            print *,"Kappa integration - intgrnd= ",kspl%fs(:5,1),"...",kspl%fs(nk-5:,1)
        endif
        
        deallocate(fm)
        call spline_dealloc(kspl)
        return
    end function kappaintgrl
    
    
    !=======================================================================
    function kappaintgnd(kappa,n,l,q,marray,dbarray,fnml)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Kappa integrand from Eq. (14) [Park, Phys. Rev. Let. 2009].
    !
    !*ARGUMENTS:
    !    n : integer (in)
    !       Toroidal mode.
    !    l : integer (in)
    !       Bounce harmonic.
    !    q : real (in)
    !       Safey factor.
    !    marray : reals
    !       Poloidal mode numbers.
    !    dbarray : real.
    !       Corresponding deltaB/B values.
    !    fnml : bicube_type.
    !       Spline of special function from [Park PRL 2009]
    !
    !*RETURNS:
    !     real.
    !        Kappa integrand.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        real(r8) :: kappaintgnd
        ! declare arguments
        integer, intent(in) :: n,l
        real(r8), intent(in) :: q
        integer, dimension(:), intent(in) :: marray
        complex(r8), dimension(:) :: dbarray
        real(r8), dimension(:,:), allocatable :: fm
        type(bicube_type) :: fnml
        ! declare local variables
        integer, parameter :: nk = 50
        integer :: i,j,mstart,mpert
        real(r8) :: kappa
        type(spline_type) :: kspl
        
        kappaintgnd = 0
        mpert = size(marray,1)
        mstart = lbound(marray,1)
        allocate(fm(mpert,2))
        ! calculate forward and backwards banana for each m
        do i = mstart,mstart+mpert-1
           ! don't let things go beyond precalculated spline
           if(abs(marray(i)-n*q-l) > fnml%ys(fnml%my))then 
              print *,marray(i),n*q,l
              stop "ERROR: pitch - RLAR approximation m-nq-l out of Fkmnql range"
           endif
           ! calculate for forward and backward banana
           call bicube_eval(fnml,kappa,marray(i)-n*q-l,0)
           fm(i,1) = fnml%f(1)
           call bicube_eval(fnml,kappa,marray(i)-n*q+l,0)
           fm(i,2) = fnml%f(1)
        enddo
        ! sum over m,m'
        do i = 1,mpert
           do j = 1,mpert
              kappaintgnd = kappaintgnd&
                +2*kappa*0.25&
                *(fm(i,1)*fm(j,1)&
                +fm(i,2)*fm(j,1)&
                +fm(i,1)*fm(j,2)&
                +fm(i,2)*fm(j,2))&
                *real(dbarray(i)*conjg(dbarray(j)))&
                /(4*ellipk(kappa))
           enddo
        enddo
        
        deallocate(fm)
        return
    end function kappaintgnd
    

end module pitch_integration

