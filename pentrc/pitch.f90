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
    
    use params, only : r8,xj,mp,e, nlambda_out
    use utilities, only : get_free_file_unit,append_2d
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
        lambdaintgrl_lsode,kappaintgrl,kappaintgnd, &       ! functions
        output_pitch_record                                 ! subroutines
    
    ! global variables with defaults
    logical :: lambdadebug = .false.
    real(r8) :: &
        lambdaatol = 1e-12, &
        lambdartol = 1e-9
    real(r8), dimension(:,:), allocatable :: pitch_record
    ! global variables for internal use
    real(r8) :: pitch_psi !wn_g,wt_g,we_g,nuk_g,l_g,n_g,bobmax_g,epsr_g
    character(4) :: pitch_method
    type(cspline_type) :: eqspl_g
    type(spline_type) :: turns_g
    
    contains

    !=======================================================================
    function lambdaintgrl_lsode(wn,wt,we,nuk,bobmax,epsr,q,eq_spl,l,n,&
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
    !   l : real.
    !       effective bounce harmonic (i.e. ell-nq for passing)
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
        integer, intent(in) :: n,l
        real(r8), intent(in) :: wn,wt,we,nuk,bobmax,epsr,q,psi,rex,imx
        character(*), intent(in) :: method
        type(cspline_type) eq_spl      
        type(spline_type) :: turns
        logical, optional :: op_record
        ! declare variables
        logical :: record_this
        integer :: i,sigma,l_g,n_g,out_unit
        real(r8) :: lmda,wb,wd,nueff,lnq,wn_g,wt_g,we_g,nuk_g,&
            bobmax_g,epsr_g,q_g,rex_g,imx_g
        complex(r8) :: xint
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
        turns_g = turns
        rex_g = rex
        imx_g = imx

        pitch_psi = psi
        pitch_method = method

        ! set default recording flag
        record_this = .false.
        if(present(op_record)) record_this = op_record

        if(record_this) then
            itask = 5              ! single step
            do while (x<xout)
                call lsode1(lintgrnd, neqarray, y, x, xout, itol, rtol,&
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                call dintdy1(x, 1, rwork(21), neqarray(1), dky, iflag)
                call spline_eval(turns_g,x,0)
                call append_2d(pitch_record, (/psi,l*1.0_r8,x,dky(1:2),y(1:2), &
                    real(eqspl_g%f(1:2)),real(eqspl_g%f(eqspl_g%nqty)),turns_g%f(:)/) )
            enddo
            ! write select energy integral output files
            do i = 0,nlambda_out-1
                lmda = eq_spl%xs(0) + (i/(nlambda_out-1.0))*(xout-eq_spl%xs(0))
                call cspline_eval(eq_spl,lmda,0)
                wb = real(eq_spl%f(1))
                wd = real(eq_spl%f(2))
                if(lmda>bobmax)then
                    nueff = nuk/(2*epsr)
                    lnq = 1.0*l
                    sigma = 0
                else
                    nueff = nuk
                    lnq = l+n*q
                    sigma = 1
                endif
                ! note: currently ignores -wb case for trapped
                xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,l,sigma,n,psi,lmda,method,record_this)
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
        lambdaintgrl_lsode(:) = y(1::2)+xj*y(2::2)

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
    
        integer :: l,n,i,sigma
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
            sigma = 1
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,l,sigma,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
            xint = xint + xintgrl_lsode(wn,wt,we,wd,-wb,nueff,l,sigma,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
        else                        ! trapped particles
            nueff = nuk/(2*epsr)    ! effective collisionality
            lnq = 1.0*l
            sigma = 0
            xint = xintgrl_lsode(wn,wt,we,wd,wb,nueff,l,sigma,n,pitch_psi,real(x,r8),pitch_method,op_record=.false.)
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
    subroutine output_pitch_record(n,zi,mi,electron,method)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii bounce function files.
    !
    !*ARGUMENTS:
    !    n : integer.
    !       mode number
    !    zi : integer (in)
    !       Ion charge in fundamental units (e).
    !    mi : integer (in)
    !       Ion mass (units of proton mass).
    !    electron : logical
    !       Calculate quantities for electrons (zi,mi ignored)
    !    method : string
    !       Label inserted in output file names.
    !    table : real 2D
    !       Table of values writen to file
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n, zi, mi
        character(*), intent(in) :: method
        logical :: electron

        integer :: i,out_unit
        real(r8), dimension(:,:), allocatable :: table
        character(8) :: nstring
        character(128) :: file

        ! safety net
        if(.not. allocated(pitch_record))then
            print *,'WARNING: No pitch angle integrand record available'
            return
        endif

        ! some methods don't have the turn information, but I don't want to re-format file
        allocate(table(16,size(pitch_record,dim=2)))
        table = 0
        table(1:size(pitch_record,dim=1),:) = pitch_record(:,:)

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstring,'(I8)') n
        file = "pentrc_"//trim(method)//"_pitch_n"//trim(adjustl(nstring))//".out"
        if(electron) file = file(:7)//"e_"//file(8:)
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Pitch angle integrand and functions"
        write(out_unit,*) " - variables are:   lambda =  B0*m*v_perp^2/(2B),  x = E/T"
        write(out_unit,*) " - frequencies are taken at x unity"
        write(out_unit,'(1/,1(a10,I4))') "n =",n
        write(out_unit,'(2(a10,es17.8E3))') "Ze =",zi*e,"mass =",mi*mp

        ! write column headers
        write(out_unit,'(1/,16(a17))') "psi_n","ell","Lambda","T_phi", &
                    "2ndeltaW","int(T_phi)","int(2ndeltaW)","omega_b","omega_D",&
                    "dJdJomega_b","theta_l","r_l","z_l","theta_u","r_u","z_u"

        ! write tables
        do i=1,size(table,dim=2)
            write(out_unit,'(16(es17.8E3))') table(:,i)
        enddo

        close(out_unit)
        deallocate(table,pitch_record)
        return
    end subroutine output_pitch_record

end module pitch_integration
