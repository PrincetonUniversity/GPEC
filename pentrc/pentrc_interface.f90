!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE


module pentrc_interface
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Input interface for PENTRC. Specifically designed for use in kinetic
    !   DCON, the subprogram in this module essentially emulates a PENTRC
    !   program run without any of the calculations. It should be regularly
    !   synced with the PENTRC program for new input variables and defaults.
    !
    !*PUBLIC MEMBER SUBPROGRAMS:
    !   set_pentrc           - Read input namelists and distribute module variables
    !
    !*PUBLIC DATA MEMBERS:
    !   ro, real                - Major radius (idcon)
    !   bo, real                - Field on axis (determined from sq)
    !   kin, spline_type        - (psi) Kinetic profiles
    !   xs_m, fspline_type      - (psi,theta) First order perturbations
    !   eqfun, bicube_type      - (psi,theta) Equilibrium mod b
    !   sq, spline_type         - (psi,theta) Equilibrium flux functions
    !   rzphi, bicube_type      - (psi,theta) Cylindrical coordinate & Jacobian
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writting. 
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only: r8,xj
    use utilities, only : get_free_file_unit
    use inputs, only : read_kin,verbose
    use energy_integration, only: &
        xatol,xrtol,xmax,ximag,xnufac,& ! reals
        xnutype,xf0type,        &       ! character(32)
        xdebug                          ! logical
    use pitch_integration, only: &
        lambdaatol,lambdartol,&         ! reals
        lambdadebug                     ! logical
    use torque, only : tintgrl_lsode,tpsi,&
        ntheta,nlmda,&                  ! integers
        tatol,trtol,&                   ! reals
        tdebug,&                        ! logical
        mpert,mfac                      !! hacked for test writting

    implicit none
    private
    public &
        get_pentrc

    ! declarations and defaults
    integer, parameter :: nflags=18
    logical :: &
        fgar_flag=.true.,&
        tgar_flag=.false.,&
        pgar_flag=.false.,&
        rlar_flag=.false.,&
        clar_flag=.false.,&
        fcgl_flag=.false.,&
        wxyz_flag=.false.,&
        fkmm_flag=.false.,&
        tkmm_flag=.false.,&
        pkmm_flag=.false.,&
        frmm_flag=.false.,&
        trmm_flag=.false.,&
        prmm_flag=.false.,&
        fwmm_flag=.false.,&
        twmm_flag=.false.,&
        pwmm_flag=.false.,&
        ftmm_flag=.false.,&
        ttmm_flag=.false.,&
        ptmm_flag=.false.,&
        electron = .false.,&
        eq_out=.false.,&
        theta_out=.false.,&
        xlmda_out=.false.,&
        eqpsi_out=.false.,&
        equil_grid=.false.,&
        input_grid=.false.,&
        fnml_flag=.false.,&
        ellip_flag=.false.,&
        diag_flag=.false.,&
        term_flag=.false.,&
        clean=.true.,&
        flags(nflags)=.false.

    integer :: i,j,k,l,m, &
        mi=2, &
        zi=1, &
        zimp=6, &
        mimp=12, &
        nl=0, &
        tmag_in = 1,&
        jsurf_in = 0,&
        nout = 30

    real(r8) ::     &
        atol=1e-6, &
        rtol=1e-3,  &
        nfac=1.0,  &
        tfac=1.0,  &
        wefac=1.0,  &
        wdfac=1.0,  &
        wpfac=1.0,  &
        nufac=1.0,  &
        divxfac=1.0,&
        diag_psi = 0.7, &
        psilim(2) = (/0,1/),&
        psiout(30)= 0, &
        psi_out(30)= (/(i,i=1,30)/)/30.6

    complex(r8) :: tphi  = (0,0), tsurf = (0,0), teq = (0,0)
    complex(r8), dimension(:,:,:), allocatable :: wtw

    character(4) :: nstring,method,methods(nflags)
    character(512) :: &
        idconfile="euler.bin", &
        kinetic_file='kin.dat', &
        gpec_file  ="gpec_order1_n1.bin", &
        peq_file ="gpec_xclebsch_n1.out", &
        data_dir =".",&
        docs(nflags)=""
    character(32) :: &
        nutype = "harmonic",&
        f0type = "maxwellian",&
        jac_in = "default",&
        moment = "pressure"

    ! namelists
    namelist/pent_input/kinetic_file,gpec_file,peq_file,idconfile, &
        data_dir,zi,zimp,mi,mimp,nl,electron,nutype,f0type,&
        jac_in,jsurf_in,tmag_in

    namelist/pent_control/nfac,tfac,wefac,wdfac,wpfac,nufac,divxfac, &
        atol,rtol,tatol,trtol,nlmda,ntheta,ximag,xmax,psilim

    namelist/pent_output/moment,eq_out,theta_out,xlmda_out,eqpsi_out,equil_grid,input_grid,&
        fgar_flag,tgar_flag,pgar_flag,clar_flag,rlar_flag,fcgl_flag,&
        wxyz_flag,psiout,psi_out,fkmm_flag,tkmm_flag,pkmm_flag,frmm_flag,trmm_flag,prmm_flag,&
        fwmm_flag,twmm_flag,pwmm_flag,ftmm_flag,ttmm_flag,ptmm_flag,&
        term_flag,verbose,clean

    namelist/pent_admin/fnml_flag,ellip_flag,diag_flag,&
        tdebug,xdebug,lambdadebug
        
    contains    

    !=======================================================================
    subroutine get_pentrc(get_nl,get_zi,get_mi,get_wdfac,get_divxfac,&
                        get_electron,get_eq_out,get_theta_out,get_xlmda_out)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read the pentrc namelists and distribute the necessary variables 
    !   to the appropriate module circles.
    !
    !*ARGUMENTS:
    !    
    !
    !*RETURNS:
    !
    !-----------------------------------------------------------------------
        implicit none
    
        integer, intent(inout) :: get_nl,get_zi,get_mi
        real(r8), intent(inout) :: get_wdfac,get_divxfac
        logical, intent(inout) :: get_electron,get_eq_out,get_theta_out,&
            get_xlmda_out
    
        ! read interface and set modules
        print *,"Interfacing with PENTRC => v3.00"
        i = get_free_file_unit(-1)
        open(unit=i,file="pentrc.in",status="old")
        read(unit=i,nml=pent_input)
        read(unit=i,nml=pent_control)
        read(unit=i,nml=pent_output)
        read(unit=i,nml=pent_admin)
        close(i)
        
        ! distribute inputs to PENTRC module circles
        xatol = atol
        xrtol = rtol
        xnufac= nufac
        xnutype= nutype
        xf0type= f0type 
        lambdaatol = atol
        lambdartol = rtol
        
        ! set the kinetic spline
        if(tdebug)THEN
            print *,"read_kin args:"
            print *,"  ",trim(kinetic_file),zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,tdebug
        endif
        call read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,tdebug)

        ! transfer torque function variables to external program
        get_nl = nl
        get_zi = zi
        get_mi = mi
        get_wdfac = wdfac
        get_divxfac = divxfac
        get_electron = electron
        get_eq_out = eq_out
        get_theta_out = theta_out
        get_xlmda_out = xlmda_out
        
    end subroutine get_pentrc
    
end module pentrc_interface
    
    
