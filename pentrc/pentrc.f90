!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE


program pentrc
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Main program. Interfaces with input file, sets corresponding
    !   global variables, and runs requested subproccesses.
    !
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writting. 
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only: r8,xj
    use utilities, only: timer
    use special, only: set_fymnl,set_ellip
    use inputs, only : read_kin,read_equil,nn,read_peq,&
                       read_ipec_peq,read_fnml,verbose
    use diagnostics, only: diagnose_all
    
    use energy_integration, only: &
        xatol,xrtol,xmax,ximag,xnufac,&  ! reals
        xnutype,xf0type,        &       ! character(32)
        qt,xdebug                          ! logical
    use pitch_integration, only: &
        lambdaatol,lambdartol,&         ! reals
        lambdadebug                     ! logical
    use torque, only : tintgrl_lsode,tintgrl_eqpsi,tpsi,&
        ntheta,nlmda,&                  ! integers
        tatol,trtol,&                   ! reals
        tdebug,&                          ! logical
        mpert,mfac                    !! hacked for test writting

    implicit none

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
        fnml_flag=.false.,&
        ellip_flag=.false.,&
        diag_flag=.false.,&
        term_flag=.true.,&
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
        psiout(30)= (/(i,i=1,30)/)/30.6
        
    complex(r8) :: tphi  = (0,0), tsurf = (0,0), teq = (0,0)
    complex(r8), dimension(:,:,:), allocatable :: wtw
        
    character(4) :: nstring,method,methods(nflags)
    character(512) :: &
        idconfile="euler.bin", &
        kinetic_file='kin.dat', &
        ipec_file  ="ipec_order1_n1.bin", &
        peq_file ="ipec_xclebsch_n1.out", &
        data_dir =".",&
        docs(nflags)=""
    character(32) :: &
        nutype = "harmonic",&
        f0type = "maxwellian",&
        jac_in = "",&
        moment = "pressure"
      
    namelist/pent_input/kinetic_file,ipec_file,peq_file,idconfile, &
        data_dir,zi,zimp,mi,mimp,nl,electron,nutype,f0type,&
        jac_in,jsurf_in,tmag_in
        
    namelist/pent_control/nfac,tfac,wefac,wdfac,wpfac,nufac,divxfac, &
        atol,rtol,tatol,trtol,nlmda,ntheta,ximag,xmax,psilim
        
    namelist/pent_output/moment,eq_out,theta_out,xlmda_out,eqpsi_out,&
        fgar_flag,tgar_flag,pgar_flag,clar_flag,rlar_flag,fcgl_flag,&
        wxyz_flag,psiout,fkmm_flag,tkmm_flag,pkmm_flag,frmm_flag,trmm_flag,prmm_flag,&
        fwmm_flag,twmm_flag,pwmm_flag,ftmm_flag,ttmm_flag,ptmm_flag,&
        term_flag,clean
        
    namelist/pent_admin/fnml_flag,ellip_flag,diag_flag,&
        tdebug,xdebug,lambdadebug

    
    ! read interface and set modules
    open(unit=1,file="pentrc.in",status="old")
    read(unit=1,nml=pent_input)
    read(unit=1,nml=pent_control)
    read(unit=1,nml=pent_output)
    read(unit=1,nml=pent_admin)
    close(1)
    
    ! distribute some simplified inputs to module circles
    xatol = atol
    xrtol = rtol
    xnufac= nufac
    xnutype= nutype
    xf0type= f0type 
    lambdaatol = atol
    lambdartol = rtol    
    verbose = term_flag
    if(verbose) print *,''
    if(verbose) print *,"PENTRC START => v3.00"
    if(verbose) print *,"______________________________"
    if(moment=="heat")then
        qt = .true.
        if(verbose) print *,"Heat transport calculation"
        if(verbose) print *,"------------------------------"
    elseif(moment=="pressure")then
        qt = .false.
        if(verbose) print *,"Particle transport and torque"
        if(verbose) print *,"------------------------------"
    else
        stop "ERROR: Input moment must be 'pressure' or 'heat'"
    endif
    
    ! start timer
    call timer(mode=0)
    ! clear working directory
    if(clean)then
        if(verbose) print *, "clearing working directory"
        call system ('rm pentrc_*.out')
    endif
    
    ! administrative setup/diagnostics/debugging.
    if(fnml_flag) call set_fymnl
    if(ellip_flag) call set_ellip
    if(diag_flag)then
        call diagnose_all
    else
        
    ! run models
        call read_equil(idconfile)
        call read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,tdebug)
        call read_peq(peq_file,jac_in,jsurf_in,tmag_in,tdebug)
        !call read_ipec_peq(ipec_file,tdebug) 
        if(wxyz_flag)then
            if(verbose) print *,"PENTRC - euler-lagrange matrix calculation"
            !! HACK - this should have its own flag
            allocate(wtw(mpert,mpert,6))
            write(nstring,'(I4)') nn
            open(unit=1,file="pentrc_tgar_elmat_n"//trim(adjustl(nstring))//".out",status="new")
            write(1,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
            write(1,*) "Kinetic additions to the ideal Euler-Lagrange matrices"
            write(1,'(1/,4x,2(a4,i4),1/)') "n = ",nn," l = ",0
            
            do i=1,nout
                if(verbose) print *," psi = ",psiout(i)
                write(1,'(1/,1x,a16,es16.8e3,1/)') "psi = ",psiout(i)
                write(1,'(1x,a4,a4,12(1x,a16))') "m_1","m_2","real(A_k)",&
                 "imag(A_k)","real(B_k)","imag(B_k)","real(C_k)","imag(C_k)",&
                 "real(D_k)","imag(D_k)","real(E_k)","imag(E_k)","real(H_k)",&
                 "imag(H_k)"
                do l=0,0!! should be all
                    if(psiout(i)>0 .and. psiout(i)<=1)then
                        tsurf = tpsi(psiout(i),nn,l,zi,mi,wdfac,divxfac,electron,'twmm',&
                            .false.,theta_out,xlmda_out,wtw)
                        do j=1,mpert
                            do k=1,mpert
                                write(1,'(1x,i4,i4,12(1x,es16.8e3))') &
                                    mfac(k),mfac(j),wtw(k,j,:)
                            enddo
                        enddo
                    endif
                enddo
            enddo
            close(1)
        endif
        flags  =(/&
                fgar_flag,tgar_flag,pgar_flag,rlar_flag,clar_flag,fcgl_flag,&
                fwmm_flag,twmm_flag,pwmm_flag,ftmm_flag,ttmm_flag,ptmm_flag,&
                fkmm_flag,tkmm_flag,pkmm_flag,frmm_flag,trmm_flag,prmm_flag &
                /)
        methods=(/&
                "fgar","tgar","pgar","rlar","clar","fcgl",&
                "fwmm","twmm","pwmm","ftmm","ttmm","ptmm",&
                "fkmm","tkmm","pkmm","frmm","trmm","prmm" &
                /)
        docs   =(/&
                "Full general-aspect-ratio calculation                       ",&
                "Trapped particle general-aspect-ratio calculation           ",&
                "Passing particle general-aspect-ratio calculation           ",&
                "Trapped particle large-aspect-ratio calculation             ",&
                "Trapped particle cylindrical large-aspect-ratio calculation ",&
                "Fluid Chew-Goldberger-Low calculation                       ",&
                "Full    energy calculation using MXM euler lagrange matrix  ",&
                "Trapped energy calculation using MXM euler lagrange matrix  ",&
                "Passing energy calculation using MXM euler lagrange matrix  ",&
                "Full    torque calculation using MXM euler lagrange matrix  ",&
                "Trapped torque calculation using MXM euler lagrange matrix  ",&
                "Passing torque calculation using MXM euler lagrange matrix  ",&
                "Full    MXM euler lagrange energy matrix norm calculation   ",&
                "Trapped MXM euler lagrange energy matrix norm calculation   ",&
                "Passing MXM euler lagrange energy matrix norm calculation   ",&
                "Full    MXM euler lagrange torque matrix norm calculation   ",&
                "Trapped MXM euler lagrange torque matrix norm calculation   ",&
                "Passing MXM euler lagrange torque matrix norm calculation   " &
                /)
        do m=1,nflags
            print *,m
            if(flags(m))then
                method = methods(m) !to_upper(methods(m))
                if(verbose)then
                    print *, "---------------------------------------------"
                    print *,method//" - "//TRIM(docs(m))
                ENDIF
                tphi = tintgrl_lsode(psilim,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m),eq_out)
                if(verbose) then
                    print "(a24,es11.3E3)", "Total torque = ", REAL(tphi)
                    print "(a24,es11.3E3)", "Total Kinetic Energy = ", AIMAG(tphi)/(2*nn)
                    print "(a24,es11.3E3)", "alpha/s  = ", REAL(tphi)/(-1*AIMAG(tphi))
                endif
                if(eqpsi_out)then
                    if(verbose) print *,method//" - "//"Recalculating on equilibrium grid"
                    teq = tintgrl_eqpsi(psilim,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m),eq_out)
                    if(verbose)then
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total torque = ", REAL(teq),&
                            ", % error = ",ABS(REAL(teq)-REAL(tphi))/REAL(tphi)
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn),&
                            ", % error = ",ABS(AIMAG(teq)-AIMAG(tphi))/AIMAG(tphi)
                        print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                    endif
                endif
                ! run select surfaces with detailed output
                if(verbose) print *,method//" - "//"Recalculating on psiout grid for detailed outputs"
                do i=1,nout
                    do l=-nl,nl,max(1,nl)
                        if(psiout(i)>0 .and. psiout(i)<=1)then
                            tsurf = tpsi(psiout(i),nn,l,zi,mi,wdfac,divxfac,electron,methods(m),&
                                .false.,theta_out,xlmda_out)
                        endif
                    enddo
                if(verbose)then
                    print *, method//" - Finished"
                    print *, "---------------------------------------------"
                endif
                enddo
            endif
        enddo
        
    endif
    
    ! display timer and stop
    call timer(mode=1)
    stop "PENTRC STOP=> normal termination."
end program pentrc
