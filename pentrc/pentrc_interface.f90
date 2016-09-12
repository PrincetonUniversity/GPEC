!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE


module pentrc_interface
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Top level interface for PENTRC.
    !   This was essentially the PENTRC program without any of the calculations,
    !   and was moved to an independent module in order to allow kinetic DCON
    !   to inteface with and use the top-tier PENTRC functions / routines.
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
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only: r8,xj, npsi_out, nmethods, methods, docs
    use utilities, only: timer,to_upper,get_free_file_unit
    use special, only: set_fymnl,set_ellip
    use dcon_interface, only: set_eq, idcon_harvest
    use inputs, only : read_kin,read_equil,nn,read_peq,read_pmodb,&
                       set_peq,read_fnml,verbose
    use diagnostics, only: diagnose_all

    use energy_integration, only: &
        output_energy_netcdf,output_energy_ascii,&  ! subroutine
        xatol,xrtol,xmax,ximag,xnufac, &            ! real
        xnutype,xf0type,&                           ! character
        qt,xdebug                                   ! logical
    use pitch_integration, only: &
        output_pitch_netcdf,output_pitch_ascii,&                       ! subroutines
        lambdaatol,lambdartol,&                     ! real
        lambdadebug                                 ! logical
    use torque, only : &
        tintgrl_lsode,tintgrl_grid,tpsi,&           ! functions
        output_orbit_ascii,&
        output_orbit_netcdf,output_torque_netcdf,&  ! subroutines
        ntheta,nlmda,nthetafuns,&                   ! integers
        atol_psi,rtol_psi,&                         ! real
        tdebug,output_ascii,output_netcdf,&         ! logical
        mpert,mfac                                  !! hacked for test writing
    use global_mod, only: version                   ! GPEC package

    implicit none

    ! declarations and defaults
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
        flags(nmethods)=.false.,&
        indebug=.false.

    integer :: i, &
        mi=2, &
        zi=1, &
        zimp=6, &
        mimp=12, &
        nl=0, &
        tmag_in = 1,&
        jsurf_in = 0,&
        power_bin = -1,&
        power_bpin = -1,&
        power_rin = -1,&
        power_rcin = -1

    real(r8) ::    &
        atol_xlmda=1e-6, &
        rtol_xlmda=1e-3, &
        nfac=1.0,  &
        tfac=1.0,  &
        wefac=1.0, &
        wdfac=1.0, &
        wpfac=1.0, &
        nufac=1.0, &
        divxfac=1.0, &
        diag_psi = 0.7, &
        psi_out(npsi_out) = -1, &
        psilims(2) = (/0,1/)
    !real(r8), dimension(npsi_out) :: psi_out
    complex(r8) :: tphi  = (0,0), tsurf = (0,0), teq = (0,0)
    complex(r8), dimension(:,:,:), allocatable :: wtw

    character(4) :: nstring,method
    character(512) :: &
        idconfile="euler.bin", &
        kinetic_file='kin.dat', &
        gpec_file  ="gpec_order1_n1.bin", &
        peq_file ="gpec_xclebsch_n1.out", &
        pmodb_file ="none", &
        data_dir ="."
    character(32) :: &
        nutype = "harmonic",&
        f0type = "maxwellian",&
        jac_in = "default",&
        moment = "pressure"

    ! namelists
    namelist/pent_input/kinetic_file,gpec_file,peq_file,pmodb_file,idconfile, &
        data_dir,zi,zimp,mi,mimp,nl,electron,nutype,f0type,&
        jac_in,jsurf_in,tmag_in,power_bin,power_bpin,power_rin,power_rcin

    namelist/pent_control/nfac,tfac,wefac,wdfac,wpfac,nufac,divxfac, &
        atol_xlmda,rtol_xlmda,atol_psi,rtol_psi,nlmda,ntheta,ximag,xmax,psilims

    namelist/pent_output/moment,output_ascii,output_netcdf,&
        eq_out,theta_out,xlmda_out,eqpsi_out,equil_grid,input_grid,&
        fgar_flag,tgar_flag,pgar_flag,clar_flag,rlar_flag,fcgl_flag,&
        wxyz_flag,psi_out,fkmm_flag,tkmm_flag,pkmm_flag,frmm_flag,trmm_flag,prmm_flag,&
        fwmm_flag,twmm_flag,pwmm_flag,ftmm_flag,ttmm_flag,ptmm_flag,&
        term_flag,verbose,clean

    namelist/pent_admin/fnml_flag,ellip_flag,diag_flag,&
        tdebug,xdebug,lambdadebug,indebug

    private :: i ! these conflict with the dcon namespace

    contains    

    !=======================================================================
    subroutine initialize_pentrc(op_kin,op_deq,op_peq)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Read the pentrc namelists and distribute the necessary variables
    !   to the appropriate module circles.
    !   Optionally, read the other standard input files (otherwise they must be set)
    !
    !*ARGUMENTS:
    !
    !*OPTIONAL ARGUMENTS:
    !   op_kin : logical (default true)
    !       Read ascii table of kinetic profiles
    !   op_deq : logical (default true)
    !       Read dcon euler.bin description of the equilibrium
    !   op_peq : logical (default true)
    !       Read ascii table of clebsch displacement profiles
    !
    !*RETURNS:
    !
    !-----------------------------------------------------------------------
        implicit none
        logical, optional, intent(in) :: op_kin,op_deq,op_peq
        logical :: in_kin,in_deq,in_peq
        integer :: i

        ! defaults
        in_kin = .true.
        in_deq = .true.
        in_peq = .true.
        if(present(op_kin)) in_kin = op_kin
        if(present(op_deq)) in_deq = op_deq
        if(present(op_peq)) in_peq = op_peq

        ! read interface and set modules
        i = get_free_file_unit(-1)
        open(unit=i,file="pentrc.in",status="old")
        read(unit=i,nml=pent_input)
        read(unit=i,nml=pent_control)
        read(unit=i,nml=pent_output)
        read(unit=i,nml=pent_admin)
        close(i)

        ! defaults
        if(all(psi_out==-1)) psi_out = (/(i,i=1,30)/)/30.6 ! even spread if user doesn't pick any

        ! warnings if using deprecated inputs
        if(eq_out) print *, "WARNING: eq_out has been deprecated. Behavior is always true."
        if(eqpsi_out) print *, "WARNING: eqpsi_out has been deprecated. Use equil_grid."
        if(term_flag) print *, "WARNING: term_flag has been deprecated. Use verbose."
        !if(any(psiout/=0)) print *, "WARNING: psiout has been deprecated. Use psi_out." !! officially removed
        !if(any(psilim/=0)) print *, "WARNING: psilim has been deprecated. Use psilims."

        ! distribute some simplified inputs to module circles
        xatol = atol_xlmda
        xrtol = rtol_xlmda
        xnufac= nufac
        xnutype= nutype
        xf0type= f0type
        lambdaatol = atol_xlmda
        lambdartol = rtol_xlmda

        ! read (perturbed) equilibrium inputs
        if(indebug .and. in_kin) print *,"  read_kin args: ",trim(kinetic_file),zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,indebug
        if(in_deq) call read_equil(idconfile)
        if(in_kin) call read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,indebug)
        if(in_peq)then
            call read_peq(peq_file,jac_in,jsurf_in,tmag_in,indebug,&
                          op_powin=(/power_bin,power_bpin,power_rin,power_rcin/))
            !call read_gpec_peq(gpec_file,indebug)
            if(trim(pmodb_file)/="" .and. trim(pmodb_file)/="none")&
                call read_pmodb(pmodb_file,jac_in,jsurf_in,tmag_in,indebug,&
                                op_powin=(/power_bin,power_bpin,power_rin,power_rcin/))
        endif

    end subroutine initialize_pentrc

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
        integer :: i

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
        xatol = atol_xlmda
        xrtol = rtol_xlmda
        xnufac= nufac
        xnutype= nutype
        xf0type= f0type 
        lambdaatol = atol_xlmda
        lambdartol = rtol_xlmda
        
        ! set the kinetic spline
        if(indebug)THEN
            print *,"read_kin args:"
            print *,"  ",trim(kinetic_file),zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,tdebug
        endif
        call read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,indebug)

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
    

