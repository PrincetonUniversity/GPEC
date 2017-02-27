!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module energy_integration
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Functional form of the energy integrand as found in
    !   [Logan, Park, et al., Phys. Plasma, 2013] and various integration
    !   methods.
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !   xintgnd             - Resonance operator Eq. (20)
    !   xintgrl_lsode       - Dynamic energy integration
    !   xintrgl_spline      - Spline integration 
    !
    !*PUBLIC DATA MEMBERS:
    !   xatol               - absolute tolerance of lsode
    !   xrtol               - relative tolerance of lsode
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writting. 
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only : r8,mp,e
    use utilities, only : get_free_file_unit,append_2d
    
    use lsode2_mod
    
    implicit none
    
    private
    public &
        xatol,xrtol,xmax,ximag,xnufac, &    ! reals
        xnutype,xf0type, &                  ! characters
        qt,xdebug, &                        ! logical
        xintgrl_lsode, &                    ! functions
        output_energy_record, &             ! subroutines
        energy_record
    
    ! global variables with defaults
    real(r8) :: &
        xatol = 1e-12, &
        xrtol = 1e-9, &
        xmax  = 72.0, &
        ximag = 0.00, &
        xnufac = 1.00
    real(r8), dimension(:,:), allocatable :: energy_record
    character(32) :: &
        xnutype = "harmonic", &
        xf0type = "maxwellian"
    logical :: &
        qt     = .false., &
        xdebug = .false.
    ! global variables for internal use
    complex(r8), parameter :: xj = (0,1)
    !real(r8) :: wn_g,wt_g,we_g,wd_g,wb_g,nuk_g,l_g,n_g
    logical :: imaxis_g = .false.

    real(r8) xcom_real
    common /xcom/ xcom_real(8)
!$OMP THREADPRIVATE(energy_record)
!$OMP THREADPRIVATE(/xcom/)

    contains

    !=======================================================================
    function xintgrl_lsode(wn,wt,we,wd,wb,nuk,l,n,psilmda)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Dynamic energy integration using lsode. Module global variabls
    !   xatol, xrtol are used for tolerances and xmax, ximag are used
    !   to define integration path (0->i*ximag, i*ximag->xmax+i*ximag).
    !
    !   Module global variable xnutype is used to determine collision
    !   operator. Default is harmonic, and valid choices are:
    !       "zero" : collisionless
    !       "small" : 1e-6*(we+wd) (krook value ignored)
    !       "krook" : krook operator (unmodified argument)
    !       "harmoinc" : (1+(l/2)^2)*krook*x^-3/2
    !
    !   Module's global variable xf0type is used to determine
    !   denominator. Default is "maxwellian" and valid options ar:
    !       "maxwellian" : x^(5/2)*exp(-x)
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
    !   wd : real.
    !       magnetic precession frequency
    !   wb : real.
    !       bounce frequency divided by x (x=E/T)
    !   nuk : real.
    !       effective Krook collision frequency (i.e. /2eps for trapped)
    !   l : real.
    !       effective bounce harmonic (i.e. ell-nq for passing)
    !   n : integer.
    !       toroidal mode number
    !*OPTIONAL ARGUMENTS:
    !   psilmda : real(2).
    !       normalized flux and normalized pitch angle (muB0/E). If this
    !       variable is used, integrand and integral are recorded in memory
    !
    !*RETURNS:
    !     complex.
    !        energy integral.
    !-----------------------------------------------------------------------
        implicit none
        ! declare arguments
        complex(r8) :: xintgrl_lsode
        real(r8), intent(in) :: wn,wt,we,wd,wb,nuk,l
        integer, intent(in) :: n
        real(r8), dimension(2), intent(in), optional :: psilmda
        ! declare variables
        integer :: xout_unit
        real(r8) :: wn_g,wt_g,we_g,wd_g,wb_g,nuk_g,l_g,n_g
        ! declare lsode input variables
        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, MF, IFLAG, i
        INTEGER, PARAMETER ::   &
            NEQ(1) = 2,         &   ! true number of equations
            LIW  = 20 + NEQ(1), &   ! for MF 22 ! only uses 20 if MF 10
            LRW  = 20 + 16*NEQ(1)   ! for MF 10 
            !LRW = 22+9*NEQ(1)+NEQ(1)**2 ! for MF 22
        INTEGER IWORK(LIW)
        REAL*8  ATOL(NEQ(1)),RTOL(NEQ(1)),RWORK(LRW),X,XOUT,&
                Y(NEQ(1)),YI(NEQ(1)),DKY(NEQ(1))
        
        common /xcom/ wn_g,wt_g,we_g,wd_g,wb_g,nuk_g,l_g,n_g
        
!$OMP THREADPRIVATE(/xcom/)

        ! set lsode options - see lsode package for documentation
        Y(:) = 0
        YI(:) = 0
        X = 1e-15
        XOUT = xmax
        ITOL = 2                  ! RTOL and ATOL are arrays
        RTOL(:) = xrtol              !1.D-7!9              ! 14
        ATOL(:) = xatol              !1.D-7!9              ! 15
        ISTATE = 1                ! first step
        IOPT = 1                  ! optional inputs
        IWORK(:) = 0              ! defaults
        RWORK(:) = 0              ! defaults
        RWORK(1) = xmax           ! only used if itask 4,5
        IWORK(6) = 10000          ! max number steps
        MF = 10                   ! not stiff with unknown J 
        !MF = 22                   ! stiff with unknown J 
    
        ! set common variables for access in integrand
        wn_g = wn
        wt_g = wt
        we_g = we
        wd_g = wd
        wb_g = wb
        nuk_g = nuk
        l_g  = l
        n_g =1.0*n
        
        ! integration to xmax
        imaxis_g = .false.
        if(present(psilmda)) then
            itask = 2              ! single step
            do while (x<xout)
                call lsode2(xintgrnd, neq, y, x, xout, itol, rtol,&
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                call dintdy2(x, 1, rwork(21), neq(1), dky, iflag)
                call append_2d(energy_record, (/psilmda(:),x,0.0_r8,l,dky(1:2),y(1:2)/) )
            enddo
        else
            itask = 1              ! full integral
            call lsode2(xintgrnd, neq, y, x, xout, itol, rtol,atol, &
                itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
        endif

        ! write error file and stop program if integration fials
        if(iwork(11)>iwork(6)/2 .and. istate/=-1) then
            print *, "WARNING: ",iwork(11)," of maximum ",iwork(6)," steps in x integration"
        endif
        if(istate==-1) then
            xout_unit = get_free_file_unit(-1)
            open(unit=xout_unit,file="pentrc_xintrgl_lsode.err",status="unknown")
            if(present(psilmda)) then
                write(xout_unit,*) "psi = ",psilmda(1)," lambda = ",psilmda(1)
            endif
            write(xout_unit,'(5(1x,a16))') "x","T_phi","2ndeltaW","int(T_phi)","int(2ndeltaW)"
            itask = 2
            y(1:2) = (/ 0,0 /)
            x = 1e-15              !1e-12
            istate = 1             ! first step
            iwork(:) = 0           ! defaults
            rwork(:) = 0           ! defaults
            rwork(1) = xmax        ! only used if itask 4,5
            iwork(6) = 10000       ! max number steps
            do while (x<xout .and. iwork(11)<iwork(6))
                call lsode2(xintgrnd, neq, y, x, xout, itol, rtol, &
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
               call dintdy2(x, 1, rwork(21), neq(1), dky, iflag)
               write(xout_unit,'(1x,5(1x,es16.8e3))') x,dky(1:2),y(1:2)
            enddo
            close(xout_unit)
            stop "ERROR: xintgrl_lsode - too many steps in x required. &
                &consider complex contour (ximag>0)."
        endif


        ! imaginary axis contribution
        if(ximag/=0.0)then
            imaxis_g = .true.
            x = 1e-15
            xout = ximag
            iwork(:) = 0 
            iwork(6) = 5000        ! max number of steps
            rwork(:) = 0
            rwork(1) = ximag       ! only relavent if using crit task
            istate = 1
            if(present(psilmda)) then
                itask = 2              ! single step
                do while (x<xout)
                    call lsode2(xintgrnd, neq, yi, x, xout, itol, rtol,&
                        atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                    call dintdy2(x, 1, rwork(21), neq(1), dky, iflag)
                    call append_2d(energy_record, (/psilmda(:),0.0_r8,x,l,dky(1:2),y(1:2)/) )
                enddo
            else
                itask = 1              ! full integral
                call lsode2(xintgrnd, neq, yi, x, xout, itol, rtol,atol, &
                    itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
            endif
           
            if(iwork(11)>iwork(6)/2) print *, "WARNING: ",iwork(11)," of ",iwork(6),&
                " maximum steps used in imaginary x integration."
            if(istate==-1) then
                stop "ERROR: xintgrl_lsode - too many steps in x required &
                    &along imaginary axis. Consider changing ximag."
            endif
        endif

        ! convert to complex space if integrations successful
        xintgrl_lsode = y(1)+xj*y(2) + yi(1) + xj*yi(2)
        return
    end function xintgrl_lsode

    
    
    
    
    

    !=======================================================================
    subroutine xintgrnd(neq,x,y,ydot)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !    Energy integrand (y) as a function of normalize energy (x) as
    !    defined in Eq. (8) of [Logan, Park, et al., Phys. Plasma, 2013].
    !
    !    Module global variable xnutype is used to determine collision
    !    operator. Default is harmonic, and valid choices are:
    !       "zero" : collisionless
    !       "small" : 1e-6*(we+wd) (krook value ignored)
    !       "krook" : krook operator (unmodified argument)
    !       "harmoinc" : (1+(l/2)^2)*krook*x^-3/2
    !
    !    Module's global variable xf0type is used to determine
    !    denominator. Default is "maxwellian" and valid options ar:
    !       "maxwellian" : x^(5/2)*exp(-x)
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
        integer ::  neq
        real*8 x, y(neq), ydot(neq)
        !real(r8), dimension(8), optional :: wntedbk_ln
    
        complex(r8) :: denom,fx,cx,nux
        real(r8) :: wn,wt,we,wd,wb,nuk,l,n
    
        common /xcom/ wn,wt,we,wd,wb,nuk,l,n

!$OMP THREADPRIVATE(/xcom/)    

        ! use input or global variables
        !if(present(wntedbk_ln))then
        !    if(xdebug) print *,'xintgrnd -> using input freqs and modes'
        !    wn = wntedbk_ln(1)
        !    wt = wntedbk_ln(2)
        !    we = wntedbk_ln(3)
        !    wd = wntedbk_ln(4)
        !    wb = wntedbk_ln(5)
        !    nuk= wntedbk_ln(6)
        !    l  = wntedbk_ln(7)
        !    n  = wntedbk_ln(8)
        !else
        !if(xdebug) print *,'xintgrnd -> using global freqs and modes'
        !wn = wn_g 
        !wt = wt_g
        !we = we_g
        !wd = wd_g
        !wb = wb_g
        !nuk= nuk_g
        !l  = l_g
        !n  = n_g
        !endif
        
        ! complex contour determined by global variable for module
        if(imaxis_g)then
           cx = xj*x
        else
           cx = x + xj*ximag
        endif

        ! collisionality determined by global variable for module
        select case (xnutype)
            case ("zero")
                nux = 0.0
            case ("small")
                nux = 1e-5*we
            case ("krook")
                nux = nuk
            case ("harmonic")
                if(x==0)then        ! (0+0i)^-3 = (nan,nan) -> causes error
                    nux= HUGE(1.0_r8) !0.0**(-1.5) ! 0^-3 = inf -> smooth arithmatic
                else
                    nux = nuk*(1+0.25*l*l)*cx**(-1.5)
                endif
            case default
                Stop "ERROR: xintgrnd - nutype must be zero, small, krook, or harmonic"
        end select
        nux = xnufac*nux 
        
        ! zeroth order distribution behavior determined by global variable for module
        denom = xj*(l*wb*sqrt(cx)+n*(we+wd*cx))-nux
        select case (xf0type)
            ! Standard solution from [Logan, Park, et al., Phys. Plasma, 2013]
            case ("maxwellian")
                fx = (we+wn+wt*(cx-1.5))*cx**2.5*exp(-cx) /denom
            ! Jong-Kyu Park [Park,Boozer,Menard, PRL 2009] approx neoclassical offset
            case ("jkp")
                fx = (we+wn+wt*2)*cx**2.5*exp(-cx) /denom
            ! Chew-Goldberger-Low limit (we+wd -> inf)
            case ("cgl") 
                fx = cx**2.5*exp(-cx)/(xj*n) /1.0
            case default
                Stop "ERROR: xintgrnd - f0 type must be maxwellian, jkp, or cgl"
        end select
        
        ! Heat flux calculation
        if(qt) fx = (cx-2.5)*fx
        
        if(.false.)then
            print *,'nutype = ',xnutype
            print *,'f0type = ',xf0type
            print *,'ximag  = ',ximag
            print *,'imaxis = ',imaxis_g
            print *,'omegas = ',wn,wt,we,wd,wb,nuk
            print *,'n,l    = ',n,l
            print *,'x      = ',x
            print *,'fx      = ',fx
        endif
        
        ! decouple two real space solutions
        ydot(1) = real(fx)
        ydot(2) = aimag(fx)
        
        return
    end subroutine xintgrnd

    
    
    
    
    
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
    subroutine output_energy_record(n,zi,mi,electron,method)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write ascii bounce function files.
    !
    !*ARGUMENTS:
    !    n : integer.
    !       mode number
    !    zi : integer (in)
    !       Ion charge in fundemental units (e).
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
        character(8) :: nstring
        character(128) :: file

        ! safety net
        if(.not. allocated(energy_record))then
            print *,'WARNING: No energy integrand record available'
            return
        endif

        ! open and prepare file as needed
        out_unit = get_free_file_unit(-1)
        write(nstring,'(I8)') n
        file = "pentrc_"//trim(method)//"_energy_n"//trim(adjustl(nstring))//".out"
        if(electron) file = file(:7)//"e_"//file(8:)
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Energy integrand"
        write(out_unit,*) " - variables are:   lambda =  B0*m*v_perp^2/(2B),  x = E/T"
        write(out_unit,*) " - normalization is: 1/(-2n^2tn*chi'/sqrt(pi))"
        write(out_unit,'(1/,1(a10,I4))') "n =",n
        write(out_unit,'(2(a10,es17.8E3))') "Ze =",zi*e,"mass =",mi*mp

        ! write column headers
        write(out_unit,'(9(a17))') "psi_n","Lambda","real(x)","imag(x)","l_eff", &
            "T_phi","2ndeltaW","int(T_phi)","int(2ndeltaW)"

        ! write tables
        do i=1,size(energy_record,dim=2)
            write(out_unit,'(9(es17.8E3))') energy_record(:,i)
        enddo

        close(out_unit)
        deallocate(energy_record)
        return
    end subroutine output_energy_record

end module energy_integration

