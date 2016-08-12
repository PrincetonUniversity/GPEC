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
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only : r8,xj,mp,e, npsi_out, nell_out, nlambda_out, &
        nmethods, methods
    use utilities, only : get_free_file_unit,append_2d,check
    use special, only : ellipk
    use cspline_mod, only: cspline_type,cspline_eval
    use spline_mod, only: spline_type,spline_eval,spline_alloc,&
                            spline_dealloc,spline_fit,spline_int
    use bicube_mod, only: bicube_type,bicube_eval
    use energy_integration, only: xintgrl_lsode
    use lsode1_mod
    use dcon_interface, only : shotnum, shottime, machine
    use global_mod, only : version
    use netcdf
    
    implicit none
    
    private
    public &
        lambdaatol,lambdartol, &                            ! real
        lambdadebug, &                                      ! logical
        lambdaintgrl_lsode,kappaintgrl,kappaintgnd, &       ! functions
        output_pitch_netcdf,output_pitch_ascii              ! subroutines
    
    ! global variables with defaults
    logical :: lambdadebug = .false.
    real(r8) :: &
        lambdaatol = 1e-12, &
        lambdartol = 1e-9

    ! global variables for internal use
    integer, parameter :: maxstep = 10000, nfs = 14
    real(r8) :: pitch_psi !wn_g,wt_g,we_g,nuk_g,l_g,n_g,bobmax_g,epsr_g
    character(4) :: pitch_method
    type(cspline_type) :: eqspl_g
    type(spline_type) :: turns_g

    type record
        logical :: is_recorded
        integer :: psi_index, ell_index
        integer, dimension(:), allocatable :: ell
        real(r8), dimension(:), allocatable :: psi
        real(r8), dimension(:,:,:,:), allocatable :: fs
    endtype record
    type(record) :: pitch_record(nmethods)
    
    contains

    !=======================================================================
    function lambdaintgrl_lsode(wn,wt,we,nuk,bobmax,epsr,q,eq_spl,ell,n,&
                                rex,imx,psi,turns,method,op_record)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Dynamic pitch integration using lsode. Module global variables
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
    !   ell : integer.
    !       bounce harmonic
    !   n : integer.
    !       toroidal mode number
    !   rex : real.
    !       Multiplier applied to real energy integral. The component
    !       gives the torque, but should be ignored when building Euler-Lagrange
    !       matrices for kinetic DCON.
    !   imx : real.
    !       Multiplier applied to imaginary energy integral. The component
    !       gives the energy, but should be ignored when building a complex
    !       mode-coupling torque matrix.
    !   psi : real.
    !       Normalized flux.
    !       **Not used in calculation, only in record keeping.**
    !   turns : spline_type.
    !       Pitch angle spline of lower (theta,r,z) and upper (theta,r,z) bounce points.
    !       **Not used in calculation, only in record keeping.**
    !   method : character.
    !       Torque psi profile method used.
    !       **Not used in calculation, only in record keeping.**
    !*OPTIONAL ARGUMENTS:
    !   op_record : logical, default .false.
    !       Record the pitch integral profiles in memory and
    !       record select energy space profiles in memory.
    !
    !*RETURNS:
    !     complex.
    !       Integral int{dLambda f(Lambda)*int{dx Rln(x,Lambda)}, where
    !       f is input in the eq_spl and Rln is the resonance operator.
    !-----------------------------------------------------------------------
        implicit none
        ! declare function
        complex(r8), dimension(:), allocatable :: lambdaintgrl_lsode        
        ! declare arguments
        integer, intent(in) :: n,ell
        real(r8), intent(in) :: wn,wt,we,nuk,bobmax,epsr,q,psi,rex,imx
        character(*), intent(in) :: method
        type(cspline_type) eq_spl      
        type(spline_type) :: turns
        logical, optional :: op_record
        ! declare variables
        logical :: record_this
        integer :: istep, ilambda, j, l_g, n_g, out_unit
        real(r8) :: lmda,wb,wd,nueff,lnq,wn_g,wt_g,we_g,nuk_g,&
            bobmax_g,epsr_g,q_g,rex_g,imx_g
        complex(r8) :: xint
        real(r8), dimension(nfs,maxstep) :: fs
        ! declare lsode input variables
        integer  iopt, iout, istate, itask, itol, mf, iflag, neqarray(1), &
            neq, liw, lrw
        integer, dimension(:), allocatable :: iwork
        real*8 :: x,xout
        real*8, dimension(:), allocatable ::  atol,rtol,rwork,y,dky

        common /lcom/ wn_g,wt_g,we_g,nuk_g,bobmax_g,epsr_g, &
                    q_g,l_g,n_g,rex_g,imx_g
        
        ! set lsode options - see lsode package for documentation
        neq = 2*(eq_spl%nqty-2)
        liw  = 20 + neq         ! for mf 22 ! only uses 20 if mf 10
        lrw  = 20 + 16*neq      ! for mf 10
        allocate(iwork(liw))
        allocate(atol(neq),rtol(neq),rwork(lrw),y(neq),dky(neq))
        neqarray(:) = (/ neq /)
        y(:) = 0                    ! initial value of integral
        x = eq_spl%xs(0)            ! lower bound of integration
        xout = eq_spl%xs(eq_spl%mx) ! upper bound of integration
        itol = 2                    ! rtol and atol are arrays
        rtol = lambdartol           !1.d-7!9              ! 14
        atol = lambdaatol           !1.d-7!9              ! 15
        istate = 1                  ! first step
        iopt = 1                    ! optional inputs
        iwork(:) = 0                ! defaults
        rwork(:) = 0                ! defaults
        rwork(1) = xout             ! only used if itask 4,5
        iwork(6) = maxstep          ! max number steps
        mf = 10                     ! not stiff with unknown J
    
        ! set module variables for access in integrand
        wn_g = wn
        wt_g = wt
        we_g = we
        nuk_g = nuk
        bobmax_g = bobmax
        epsr_g = epsr
        q_g = q
        l_g  = ell
        n_g = n
        eqspl_g = eq_spl
        turns_g = turns
        rex_g = rex
        imx_g = imx

        pitch_psi = psi
        pitch_method = method

        ! set default recording flag
        record_this = .false.
        if(present(op_record)) record_this = op_record

        if(record_this) then
            istep = 0
            itask = 5              ! single step
            do while (x<xout)
                call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol,&
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                call dintdy1(x, 1, rwork(21), neqarray(1), dky, iflag)
                call spline_eval(turns_g,x,0)
                istep = istep + 1
                fs(:,istep) = (/x, dky(1:2), y(1:2), real(eqspl_g%f(1:3)), turns_g%f(1:6)/)
            enddo
            ! write select energy integral output files
            do ilambda = 0,nlambda_out-1
                lmda = eq_spl%xs(0) + (ilambda/(nlambda_out-1.0))*(xout-eq_spl%xs(0))
                call cspline_eval(eq_spl,lmda,0)
                wb = real(eq_spl%f(1))
                wd = real(eq_spl%f(2))
                if(lmda>bobmax)then
                    nueff = nuk/(2*epsr)
                    lnq = real(ell,r8)
                else
                    nueff = nuk
                    lnq = ell+n*q
                endif
                ! note: currently ignores -wb case for trapped
                xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,ell,lnq,n,psi,lmda,method,record_this)
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
            write(out_unit,*) "psi = ",psi
            write(out_unit,'(5(1x,a16))') "lambda","t_phi","2ndeltaw","int(t_phi)","int(2ndeltaw)"
            itask = 2
            y(1:2) = (/ 0,0 /)
            x = eq_spl%xs(0)       ! lower bound
            istate = 1             ! first step
            iwork(:) = 0           ! defaults
            rwork(:) = 0           ! defaults
            rwork(1) = xout        ! only used if itask 4,5
            iwork(6) = maxstep     ! max number steps
            do while (x<xout .and. iwork(11)<iwork(6))
                call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol, &
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
               call dintdy1(x, 1, rwork(21), neqarray(1), dky, iflag)
               write(out_unit,'(1x,5(1x,es16.8e3))') x,dky(1:2),y(1:2)
            enddo
            close(out_unit)
            stop "ERROR: lambdaintgrl_lsode - too many steps in lambda required."
        endif

        ! save integration profile to memory
        if(record_this) then
            ! fill in unused x's with the endpoints
            do j=istep+1,maxstep
                fs(:,j) = fs(:,istep)
            enddo
            call record_method( method, psi, ell, fs )
        endif

        ! convert to complex space if integrations successful
        if(lambdadebug) print *,"Finished Lambda lsode. Converting to complex..."
        if(allocated(lambdaintgrl_lsode)) deallocate(lambdaintgrl_lsode)
        allocate(lambdaintgrl_lsode(neq/2))
        if(lambdadebug) print *,"...allocated lambdaintgrl",neq/2,"..."        
        lambdaintgrl_lsode(:) = y(1::2)+xj*y(2::2)

        deallocate(iwork,atol,rtol,rwork,y,dky)
        if(lambdadebug) print *,"Returning."        
        
        return
    end function lambdaintgrl_lsode

    
    
    
    
    

    !=======================================================================
    subroutine lintgrnd(neq,x,y,ydot)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Pitch integrand for us in LSODE.
    !   Note: Also uses a number of common variables.
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
        
        ! energy integration of resonance operator
        if(x<=bobmax)then      ! circulating particles
            nueff = nuk
            lnq = l+n*q
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,l,lnq,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
            xint = xint + xintgrl_lsode(wn,wt,we,wd,-wb,nueff,l,lnq,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
        else                        ! trapped particles
            nueff = nuk/(2*epsr)    ! effective collisionality
            lnq = real(l,r8)
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,l,lnq,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
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


    !=======================================================================
    subroutine record_method(method, psi, ell, fs)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Save the energy integrand and integral profiles for this method,
    !   at this surface and pitch.
    !
    !*ARGUMENTS:
    !
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: ell
        real(r8), intent(in) :: psi,fs(nfs,maxstep)
        character(*), intent(in) :: method

        integer :: m, j, k
        logical, parameter :: debug = .false.

        if(debug) print *,"Recording method"

        ! find the right record
        do m=1,nmethods
            if(method==methods(m))then
                if(debug) print *,"  method "//trim(method)
                ! initialize the record of this method type if needed
                if(.not. pitch_record(m)%is_recorded)then
                    if(debug) print *,"   - is not recorded"
                    pitch_record(m)%is_recorded = .true.
                    pitch_record(m)%psi_index = 0
                    pitch_record(m)%ell_index = 0
                    allocate( pitch_record(m)%psi(npsi_out), &
                        pitch_record(m)%ell(nell_out), &
                        pitch_record(m)%fs(nfs,maxstep,nell_out,npsi_out) )
                endif
                ! bump indexes
                if(debug) print *,"   - psi,ell = ",psi,ell
                if(debug) print *,"   - old j,k = ",j,k
                k = pitch_record(m)%psi_index
                j = pitch_record(m)%ell_index
                if(psi/=pitch_record(m)%psi(k) .or. k==0) k = k+1
                if(ell/=pitch_record(m)%ell(j) .or. j==0) j = mod(j,nell_out)+1
                if(debug) print *,"   - new j,k = ",j,k
                ! force fail if buggy
                if(k>npsi_out)then
                    print *,"ERROR: Too many psi energy records for method "//trim(method)
                    stop
                endif
                ! fill in new profiles
                pitch_record(m)%fs(:,:,j,k) = fs(:,:)
                pitch_record(m)%ell(j) = ell
                pitch_record(m)%psi(k) = psi
                pitch_record(m)%psi_index = k
                pitch_record(m)%ell_index = j
            endif
        enddo
        if(debug) print *,"Done recording"
        return
    end subroutine record_method

    !=======================================================================
    subroutine record_reset()
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Deallocate all energy records.
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
            if(pitch_record(m)%is_recorded)then
                pitch_record(m)%is_recorded = .false.
                pitch_record(m)%psi_index = 0
                pitch_record(m)%ell_index = 0
                deallocate( pitch_record(m)%psi, &
                    pitch_record(m)%ell, &
                    pitch_record(m)%fs )
            endif
        enddo
        if(debug) print *,"Energy record reset done"
        return
    end subroutine record_reset


    !=======================================================================
    subroutine output_pitch_netcdf(n,op_label)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write recorded pitch angle integrations to a netcdf file.
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

        integer :: i,m,npsi,nell
        integer, dimension(:), allocatable :: ell_out
        real(r8), dimension(:), allocatable :: psi_out

        integer :: status, ncid, i_did, i_id, &
            p_did, p_id, l_did, l_id, a_did, a_id, &
            aa_id, dt_id, it_id, wb_id, wd_id, jj_id, &
            tl_id, rl_id, zl_id, tu_id, ru_id, zu_id
        character(16) :: nstring,suffix,label
        character(128) :: ncfile

        logical :: debug = .false.

        print *,"Writing pitch record output to netcdf"

        ! optional labeling
        label = ''
        if(present(op_label)) label = "_"//trim(adjustl(op_label))

        ! assume all methods on the same psi's and ell's
        npsi = 0
        nell = 0
        do m=1,nmethods
            if(pitch_record(m)%is_recorded)then
                if(debug) print *,"  Getting dims from "//trim(methods(m))
                npsi = pitch_record(m)%psi_index
                nell = pitch_record(m)%ell_index
                allocate(psi_out(npsi),ell_out(nell))
                psi_out = pitch_record(m)%psi(:npsi)
                ell_out = pitch_record(m)%ell(:nell)
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
        ncfile = "pentrc_pitch_output"//trim(label)//"_n"//trim(adjustl(nstring))//".nc"
        if(debug) print *, "  opening "//trim(ncfile)
        call check( nf90_create(ncfile, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid) )
        ! store attributes
        if(debug) print *, "  storing attributes"
        call check( nf90_put_att(ncid, nf90_global, "title", "PENTRC pitch angle integration profiles") )
        call check( nf90_put_att(ncid, nf90_global, "version", version))
        call check( nf90_put_att(ncid, nf90_global, "shot", INT(shotnum) ) )
        call check( nf90_put_att(ncid, nf90_global, "time", INT(shottime)) )
        call check( nf90_put_att(ncid, nf90_global, "machine", machine) )
        call check( nf90_put_att(ncid, nf90_global, "n", n) )
        ! define dimensions
        if(debug) print *, "  defining dimensions"
        if(debug) print *, "  npsi,nell,nlambda = ",npsi,nell,maxstep
        call check( nf90_def_dim(ncid, "i", 2, i_did) )
        call check( nf90_def_var(ncid, "i", nf90_int, i_did, i_id) )
        call check( nf90_def_dim(ncid, "Lambda_index", maxstep, a_did) )
        call check( nf90_def_var(ncid, "Lambda_index", nf90_double, a_did, a_id) )
        call check( nf90_def_dim(ncid, "ell", nell, l_did) )
        call check( nf90_def_var(ncid, "ell", nf90_int, l_did, l_id) )
        call check( nf90_def_dim(ncid, "psi_n", npsi, p_did) )
        call check( nf90_def_var(ncid, "psi_n", nf90_double, p_did, p_id) )
        ! End definitions
        call check( nf90_enddef(ncid) )
        ! store dimensions
        if(debug) print *, "  storing dimensions"
        call check( nf90_put_var(ncid, i_id, (/0,1/)) )
        call check( nf90_put_var(ncid, a_id, (/(i,i=1,maxstep)/)) )
        call check( nf90_put_var(ncid, l_id, ell_out) )
        call check( nf90_put_var(ncid, p_id, psi_out) )

        ! store each method, grid combination that has been run
        do m=1,nmethods
            if(pitch_record(m)%is_recorded)then
                if(debug) print *, "  defining "//trim(methods(m))
                ! create distinguishing labels
                suffix = '_'//trim(methods(m))
                ! check sizes
                if(npsi/=pitch_record(m)%psi_index)then
                    print *,npsi,pitch_record(m)%psi_index
                    stop "Error: Record sizes are inconsistent"
                endif
                ! Re-open definitions
                call check( nf90_redef(ncid) )
                call check( nf90_def_var(ncid, "Lambda"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), aa_id) )
                call check( nf90_def_var(ncid, "T_psi_Lambda"//trim(suffix), nf90_double, (/i_did,a_did,l_did,p_did/), dt_id) )
                call check( nf90_def_var(ncid, "T_psi"//trim(suffix), nf90_double, (/i_did,a_did,l_did,p_did/), it_id) )
                call check( nf90_def_var(ncid, "omega_b"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), wb_id) )
                call check( nf90_def_var(ncid, "omega_D"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), wd_id) )
                call check( nf90_def_var(ncid, "dJdJomega_b"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), jj_id) )
                call check( nf90_def_var(ncid, "theta_l"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), tl_id) )
                call check( nf90_def_var(ncid, "R_l"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), rl_id) )
                call check( nf90_def_var(ncid, "z_l"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), zl_id) )
                call check( nf90_def_var(ncid, "theta_u"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), tu_id) )
                call check( nf90_def_var(ncid, "R_u"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), ru_id) )
                call check( nf90_def_var(ncid, "z_u"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), zu_id) )
                call check( nf90_enddef(ncid) )
                ! Put in variables
                if(debug) print *, "  storing "//trim(methods(m))
                call check( nf90_put_var(ncid, aa_id, pitch_record(m)%fs(1,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, dt_id, pitch_record(m)%fs(2:3,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, it_id, pitch_record(m)%fs(4:5,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, wb_id, pitch_record(m)%fs(6,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, wd_id, pitch_record(m)%fs(7,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, jj_id, pitch_record(m)%fs(8,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, tl_id, pitch_record(m)%fs(9,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, rl_id, pitch_record(m)%fs(10,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, zl_id, pitch_record(m)%fs(11,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, tu_id, pitch_record(m)%fs(12,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, ru_id, pitch_record(m)%fs(13,:,:nell,:npsi)) )
                call check( nf90_put_var(ncid, zu_id, pitch_record(m)%fs(14,:,:nell,:npsi)) )
            endif
        enddo

        ! close file
        call check( nf90_close(ncid) )
        ! clear the memory
        call record_reset( )
        if(debug) print *, "Finished energy netcdf output"

        return
    end subroutine output_pitch_netcdf


    !=======================================================================
    subroutine output_pitch_ascii(n,op_label)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii bounce function files.
    !
    !*ARGUMENTS:
    !    n : integer.
    !       toroidal mode number for filename
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n
        character(*), optional :: op_label

        integer :: m,i,j,k,out_unit
        character(8) :: nstring
        character(16) :: label
        character(128) :: file

        ! optional labeling
        label = ''
        if(present(op_label)) label = "_"//trim(adjustl(op_label))

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstring,'(I8)') n
        file = "pentrc_pitch_output"//trim(label)//"_n"//trim(adjustl(nstring))//".out"
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Pitch angle integrand and functions"
        write(out_unit,*) " - variables are:   lambda =  B0*m*v_perp^2/(2B),  x = E/T"
        write(out_unit,*) " - frequencies are taken at x unity"
        write(out_unit,'(1/,1(a10,I4))') "n =",n

        ! write each method in a new table
        do m=1,nmethods
            if(pitch_record(m)%is_recorded)then
                write(out_unit,'(1/,2(a17))')"method =",trim(methods(m))
                write(out_unit,'(1/,17(a17))') "psi_n","ell","Lambda_index","Lambda", &
                    "T_phi","2ndeltaW","int(T_phi)","int(2ndeltaW)","omega_b","omega_D",&
                    "dJdJomega_b","theta_l","r_l","z_l","theta_u","r_u","z_u"
                do k=1,pitch_record(m)%psi_index
                    do j=1,pitch_record(m)%ell_index
                        do i=1,maxstep
                            write(out_unit,'(17(es17.8E3))') pitch_record(m)%psi(k), &
                                1.0*pitch_record(m)%ell(j),1.0*i,pitch_record(m)%fs(:,i,j,k)
                        enddo
                    enddo
                enddo
            endif
        enddo
        close(out_unit)

        return
    end subroutine output_pitch_ascii

end module pitch_integration
