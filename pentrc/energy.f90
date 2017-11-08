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
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only : r8, mp, me, e, npsi_out, nell_out, nlambda_out, &
        nmethods, methods
    use utilities, only : get_free_file_unit,append_2d,check
    use dcon_interface, only : shotnum, shottime, machine
    use global_mod, only : version
    use netcdf
    
    use lsode2_mod
    
    implicit none
    
    private
    public &
        maxstep, & ! integer
        xatol,xrtol,xmax,ximag,xnufac, &    ! real
        xnutype,xf0type,methods, &                  ! character
        qt,xdebug, &                        ! logical
        xintgrl_lsode, &                    ! function
        output_energy_netcdf,output_energy_ascii                ! subroutine
    
    ! global variables with defaults
    integer :: &
        maxstep = 10000
    real(r8) :: &
        xatol = 1e-12, &
        xrtol = 1e-9, &
        xmax  = 72.0, &
        ximag = 0.00, &
        xnufac = 1.00

    type record
        logical :: is_recorded
        integer :: psi_index, lambda_index, ell_index
        integer, dimension(:), allocatable :: ell
        real(r8), dimension(:), allocatable :: psi
        real(r8), dimension(:,:,:), allocatable :: lambda, leff
        real(r8), dimension(:,:,:,:,:), allocatable :: fs
    endtype record
    type(record) :: energy_record(nmethods)


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
!$OMP THREADPRIVATE(/xcom/)

    contains

    !=======================================================================
    function xintgrl_lsode(wn,wt,we,wd,wb,nuk,ell,leff,n,psi,lambda,method,op_record)
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
    !   ell : integer.
    !       Bounce harmonic
    !   leff : real.
    !       effective bounce harmonic ell-sigma*n*q where sigma=0(1) for trapped(passing)
    !   n : integer.
    !       toroidal mode number
    !   psi : real.
    !       Flux surface. Not used in calculation, only for record keeping.
    !   lmda : real.
    !       Normalized pitch angle muB0/E. Not used in calculation, only for record keeping.
    !*OPTIONAL ARGUMENTS:
    !   op_record : logical.
    !       Store integration energy-space profile in memory.
    !
    !*RETURNS:
    !     complex.
    !        energy integral.
    !-----------------------------------------------------------------------
        implicit none
        ! declare arguments
        complex(r8) :: xintgrl_lsode
        real(r8), intent(in) :: wn,wt,we,wd,wb,nuk,leff,psi,lambda
        integer, intent(in) :: ell,n
        character(4), intent(in) :: method
        logical, intent(in), optional :: op_record
        ! declare variables
        integer :: i, j, xout_unit
        real(r8) :: wn_g,wt_g,we_g,wd_g,wb_g,nuk_g,l_g,n_g
        real(r8), dimension(6,maxstep) :: xprofile
        logical :: record_this
        ! declare lsode input variables
        integer  iopt, istate, itask, itol, mf, iflag
        integer, parameter ::   &
            neq(1) = 2,         &   ! true number of equations
            liw  = 20 + neq(1), &   ! for mf 22 ! only uses 20 if mf 10
            lrw  = 20 + 16*neq(1)   ! for mf 10
            !lrw = 22+9*neq(1)+neq(1)**2 ! for mf 22
        integer iwork(liw)
        real*8  atol(neq(1)),rtol(neq(1)),rwork(lrw),x,xout,&
                y(neq(1)),yi(neq(1)),dky(neq(1))
        
        common /xcom/ wn_g,wt_g,we_g,wd_g,wb_g,nuk_g,l_g,n_g

!$OMP THREADPRIVATE(/xcom/)


        ! set default recording flag
        i = 0
        xprofile = 0.0
        record_this = .false.
        if(present(op_record)) record_this = op_record

        ! set lsode options - see lsode package for documentation
        y(:) = 0
        yi(:) = 0
        itol = 2                  ! rtol and atol are arrays
        rtol(:) = xrtol           ! relative tolerance
        atol(:) = xatol           ! absolute tolerance
        iopt = 1                  ! optional inputs
        mf = 10                   ! not stiff with unknown J

        ! set common variables for access in integrand
        wn_g = wn
        wt_g = wt
        we_g = we
        wd_g = wd
        wb_g = wb
        nuk_g = nuk
        l_g  = leff
        n_g =1.0*n
        
        ! optionally step off real axis to avoid poles
        if(ximag/=0.0)then
            imaxis_g = .true.
            x = 1e-15               ! lower bound of integration
            xout = ximag            ! upper bound of integration
            iwork(:) = 0            ! default
            iwork(6) = maxstep - i   ! max number of steps
            rwork(:) = 0            ! default
            rwork(1) = ximag        ! only relevant if using crit task (itask 4,5)
            istate = 1              ! (re)start initial step
            if(record_this) then
                itask = 2              ! single step
                do while (x<xout)
                    call lsode2(xintgrnd, neq, yi, x, xout, itol, rtol,&
                        atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                    call dintdy2(x, 1, rwork(21), neq(1), dky, iflag)
                    i = i+1
                    xprofile(:,i) = (/real(0,r8), x, dky(1:2), y(1:2)/)
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

        ! integration in real space to xmax
        imaxis_g = .false.
        x = 1e-15               ! lower bound of integration
        xout = xmax             ! upper bound of integration
        iwork(:) = 0            ! default
        iwork(6) = maxstep - i   ! max number of steps
        rwork(:) = 0            ! default
        rwork(1) = xmax         ! only relevant if using crit task (itask 4,5)
        istate = 1              ! (re)start initial step
        if(record_this) then
            itask = 2              ! single step
            do while (x<xout)
                call lsode2(xintgrnd, neq, y, x, xout, itol, rtol,&
                    atol,itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
                call dintdy2(x, 1, rwork(21), neq(1), dky, iflag)
                i = i+1
                xprofile(:,i) = (/x, real(0,r8), dky(1:2),y(1:2)/)
            enddo
        else
            itask = 1              ! full integral
            call lsode2(xintgrnd, neq, y, x, xout, itol, rtol,atol, &
                itask,istate, iopt, rwork, lrw, iwork, liw, noj, mf)
        endif

        ! write error file and stop program if integration fails
        if(iwork(11)>iwork(6)/2 .and. istate/=-1) then
            print *, "WARNING: ",iwork(11)," of maximum ",iwork(6)," steps in x integration"
        endif
        if(istate==-1) then
            xout_unit = get_free_file_unit(-1)
            open(unit=xout_unit,file="pentrc_xintrgl_lsode.err",status="unknown")
            write(xout_unit,*) "psi = ",psi," lambda = ",lambda
            write(xout_unit,'(5(1x,a16))') "x","T_phi","2ndeltaW","int(T_phi)","int(2ndeltaW)"
            itask = 2
            y(1:2) = (/ 0,0 /)
            x = 1e-15              !1e-12
            istate = 1             ! first step
            iwork(:) = 0           ! defaults
            rwork(:) = 0           ! defaults
            rwork(1) = xmax        ! only used if itask 4,5
            iwork(6) = maxstep       ! max number steps
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

        ! save integration profile to memory
        if(record_this) then
            ! fill in unused x's with the endpoints
            do j=i+1,maxstep
                xprofile(1,j) = xprofile(1,i) + 1e-9 ! last little step in real x
                xprofile(2,j) = 0 ! integral ends on real x axis
                xprofile(3:4,j) = 0 ! integrand goes to zero
                xprofile(5:6,j) = xprofile(5:6,i) ! final integral
            enddo
            call record_method( method, psi, ell, leff, lambda, xprofile )
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
    subroutine record_method(method, psi, ell, leff, lambda, fs)
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
        real(r8), intent(in) :: psi,lambda,leff,fs(6,maxstep)
        character(*), intent(in) :: method

        integer :: m, i, j, k
        logical, parameter :: debug = .false.

        if(debug) print *,"Recording method"

        ! find the right record
        do m=1,nmethods
            if(method==methods(m))then
                if(debug) print *,"  method "//trim(method)
                ! initialize the record of this method type if needed
                if(.not. energy_record(m)%is_recorded)then
                    if(debug) print *,"   - is not recorded"
                    energy_record(m)%is_recorded = .true.
                    energy_record(m)%psi_index = 0
                    energy_record(m)%lambda_index = 0
                    energy_record(m)%ell_index = 0
                    allocate( energy_record(m)%psi(npsi_out), &
                        energy_record(m)%ell(nell_out), &
                        energy_record(m)%leff(nlambda_out,nell_out,npsi_out), &
                        energy_record(m)%lambda(nlambda_out,nell_out,npsi_out), &
                        energy_record(m)%fs(6,maxstep,nlambda_out,nell_out,npsi_out) )
                endif
                ! bump indexes
                if(debug) print *,"   - psi,ell,leff,lambda = ",psi,ell,leff,lambda
                if(debug) print *,"   - old i,j,k = ",i,j,k
                k = energy_record(m)%psi_index
                j = energy_record(m)%ell_index
                i = energy_record(m)%lambda_index
                if(psi/=energy_record(m)%psi(k) .or. k==0) k = k+1
                if(ell/=energy_record(m)%ell(j) .or. j==0) j = mod(j,nell_out)+1
                if(lambda/=energy_record(m)%lambda(i,j,k) .or. i==0) i = mod(i,nlambda_out)+1
                if(debug) print *,"   - new i,j,k = ",i,j,k
                ! force fail if buggy
                if(k>npsi_out)then
                    print *,"ERROR: Too many psi energy records for method "//trim(method)
                    stop
                endif
                ! fill in new profiles
                energy_record(m)%fs(:,:,i,j,k) = fs(:,:)
                energy_record(m)%lambda(i,j,k) = lambda
                energy_record(m)%leff(i,j,k) = leff
                energy_record(m)%ell(j) = ell
                energy_record(m)%psi(k) = psi
                energy_record(m)%psi_index = k
                energy_record(m)%ell_index = j
                energy_record(m)%lambda_index = i
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
            if(energy_record(m)%is_recorded)then
                energy_record(m)%is_recorded = .false.
                energy_record(m)%psi_index = 0
                energy_record(m)%lambda_index = 0
                energy_record(m)%ell_index = 0
                deallocate( energy_record(m)%psi, &
                    energy_record(m)%ell, &
                    energy_record(m)%leff, &
                    energy_record(m)%lambda, &
                    energy_record(m)%fs )
            endif
        enddo
        if(debug) print *,"Energy record reset done"
        return
    end subroutine record_reset


    !=======================================================================
    subroutine output_energy_netcdf(n,op_label)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Write recorded energy space integrations to a netcdf file.
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
            le_id,a_did,a_id,x_did,x_id,aa_id, xx_id,dt_id, it_id
        character(16) :: nstring,suffix,label
        character(128) :: ncfile

        logical :: debug = .false.

        print *,"Writing energy record output to netcdf"

        ! optional labeling
        label = ''
        if(present(op_label)) label = "_"//trim(adjustl(op_label))

        ! assume all methods on the same psi's and ell's
        npsi = 0
        nell = 0
        nlambda = 0
        do m=1,nmethods
            if(energy_record(m)%is_recorded)then
                if(debug) print *,"  Getting dims from "//trim(methods(m))
                npsi = energy_record(m)%psi_index
                nell = energy_record(m)%ell_index
                nlambda = energy_record(m)%lambda_index
                allocate(psi_out(npsi),ell_out(nell))
                psi_out = energy_record(m)%psi(:npsi)
                ell_out = energy_record(m)%ell(:nell)
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
        ncfile = "pentrc_energy_output"//trim(label)//"_n"//trim(adjustl(nstring))//".nc"
        if(debug) print *, "  opening "//trim(ncfile)
        call check( nf90_create(ncfile, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid) )
        ! store attributes
        if(debug) print *, "  storing attributes"
        call check( nf90_put_att(ncid, nf90_global, "title", "PENTRC energy space outputs") )
        call check( nf90_put_att(ncid, nf90_global, "version", version))
        call check( nf90_put_att(ncid, nf90_global, "shot", INT(shotnum) ) )
        call check( nf90_put_att(ncid, nf90_global, "time", INT(shottime)) )
        call check( nf90_put_att(ncid, nf90_global, "machine", machine) )
        call check( nf90_put_att(ncid, nf90_global, "n", n) )
        ! define dimensions
        if(debug) print *, "  defining dimensions"
        if(debug) print *, "  npsi,nell,nlambda,nx = ",npsi,nell,nlambda,maxstep
        call check( nf90_def_dim(ncid, "i", 2, i_did) )
        call check( nf90_def_var(ncid, "i", nf90_int, i_did, i_id) )
        call check( nf90_def_dim(ncid, "ell", nell, l_did) )
        call check( nf90_def_var(ncid, "ell", nf90_int, l_did, l_id) )
        call check( nf90_def_dim(ncid, "psi_n", npsi, p_did) )
        call check( nf90_def_var(ncid, "psi_n", nf90_double, p_did, p_id) )
        call check( nf90_def_dim(ncid, "Lambda_index", nlambda, a_did) )
        call check( nf90_def_var(ncid, "Lambda_index", nf90_double, a_did, a_id) )
        call check( nf90_def_dim(ncid, "x_index", maxstep, x_did) )
        call check( nf90_def_var(ncid, "x_index", nf90_double, x_did, x_id) )
        ! End definitions
        call check( nf90_enddef(ncid) )
        ! store dimensions
        if(debug) print *, "  storing dimensions"
        call check( nf90_put_var(ncid, i_id, (/0,1/)) )
        call check( nf90_put_var(ncid, x_id, (/(i,i=1,maxstep)/)) )
        call check( nf90_put_var(ncid, a_id, (/(i,i=1,nlambda)/)) )
        call check( nf90_put_var(ncid, l_id, ell_out) )
        call check( nf90_put_var(ncid, p_id, psi_out) )

        ! store each method, grid combination that has been run
        do m=1,nmethods
            if(energy_record(m)%is_recorded)then
                if(debug) print *, "  defining "//trim(methods(m))
                ! create distinguishing labels
                suffix = '_'//trim(methods(m))
                ! check sizes
                if(npsi/=energy_record(m)%psi_index)then
                    print *,npsi,energy_record(m)%psi_index
                    stop "Error: Record sizes are inconsistent"
                endif
                ! Re-open definitions
                call check( nf90_redef(ncid) )
                call check( nf90_def_var(ncid, "Lambda"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), aa_id) )
                call check( nf90_def_var(ncid, "ell_eff"//trim(suffix), nf90_double, (/a_did,l_did,p_did/), le_id) )
                call check( nf90_def_var(ncid, "x"//trim(suffix), nf90_double, &
                    (/i_did,x_did,a_did,l_did,p_did/), xx_id) )
                call check( nf90_def_var(ncid, "T_psi_Lambda_x"//trim(suffix), nf90_double, &
                    (/i_did,x_did,a_did,l_did,p_did/), dt_id) )
                call check( nf90_def_var(ncid, "T_psi_Lambda"//trim(suffix), nf90_double, &
                    (/i_did,x_did,a_did,l_did,p_did/), it_id) )
                call check( nf90_enddef(ncid) )
                ! Put in variables
                if(debug) print *, "  storing "//trim(methods(m))
                call check( nf90_put_var(ncid, aa_id, energy_record(m)%lambda(:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, le_id, energy_record(m)%leff(:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, xx_id, energy_record(m)%fs(1:2,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, dt_id, energy_record(m)%fs(3:4,:,:nlambda,:nell,:npsi)) )
                call check( nf90_put_var(ncid, it_id, energy_record(m)%fs(5:6,:,:nlambda,:nell,:npsi)) )
            endif
        enddo

        ! close file
        call check( nf90_close(ncid) )
        ! clear the memory
        call record_reset( )
        if(debug) print *, "Finished energy netcdf output"

        return
    end subroutine output_energy_netcdf

    !=======================================================================
    subroutine output_energy_ascii(n, op_label)
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
        file = "pentrc_energy_output"//trim(label)//"_n"//trim(adjustl(nstring))//".out"
        open(unit=out_unit,file=file,status="unknown",action="write")

        ! write header material
        write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
        write(out_unit,*) " Energy integrand"
        write(out_unit,*) " - variables are:   lambda =  B0*m*v_perp^2/(2B),  x = E/T"
        write(out_unit,*) " - normalization is: 1/(-2n^2tn*chi'/sqrt(pi))"
        write(out_unit,'(1/,1(a10,I4))') "n =",n

        ! write each method in a new table
        do m=1,nmethods
            if(energy_record(m)%is_recorded)then
                write(out_unit,'(1/,2(a17))')"method =",trim(methods(m))
                write(out_unit,'(12(a17))') "psi_n","ell","Lambda_index","x_index", &
                    "Lambda","ell_eff","real(x)","imag(x)", &
                    "T_phi","2ndeltaW","int(T_phi)","int(2ndeltaW)"
                do k=1,energy_record(m)%psi_index
                    do j=1,energy_record(m)%ell_index
                        do i=1,energy_record(m)%lambda_index
                            do istep=1,maxstep
                                write(out_unit,'(12(es17.8E3))') energy_record(m)%psi(k), &
                                    real(energy_record(m)%ell(j)),real(i),real(istep), &
                                    energy_record(m)%lambda(i,j,k),energy_record(m)%leff(i,j,k),&
                                    energy_record(m)%fs(:,istep,i,j,k)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        enddo
        close(out_unit)

        return
    end subroutine output_energy_ascii

end module energy_integration

