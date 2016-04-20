!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE


program pentrc
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Main program. Interfaces with input file, sets corresponding
    !   global variables, and runs requested subprocess.
    !
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    use pentrc_interface

    ! local variables
    integer :: i,j,k,l,m

    ! harvest variables
    include 'harvest_lib.inc'
    integer :: ierr
    CHARACTER(LEN=50000) :: hnml
    character(len=65507) :: hlog
    character, parameter :: nul = char(0)

    call initialize_pentrc

    if(verbose) print *,''
    if(verbose) print *,"PENTRC START => "//trim(version)
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
        ! start log with harvest
        ierr=init_harvest('CODEDB_PENT'//nul,hlog,len(hlog))
        ierr=set_harvest_verbose(0)
        ! standard CODEDB records
        ierr=set_harvest_payload_str(hlog,'CODE'//nul,'PENT'//nul)
        ierr=set_harvest_payload_str(hlog,'VERSION'//nul,version//nul)
        ! PENT input records
        if(jac_in=="") jac_in = "default" ! harvest can't parse empty strings
        write(hnml,nml=pent_input)
        ierr=set_harvest_payload_nam(hlog,'PENT_INPUT'//nul,trim(hnml)//nul)
        write(hnml,nml=pent_control)
        ierr=set_harvest_payload_nam(hlog,'PENT_CONTROL'//nul,trim(hnml)//nul)
        write(hnml,nml=pent_output)
        ierr=set_harvest_payload_nam(hlog,'PENT_OUTPUT'//nul,trim(hnml)//nul)

        ! record dcon equilibrium basics
        call idcon_harvest(hlog)

        ! explicit matrix calculations
        if(wxyz_flag .and. output_ascii)then
            if(verbose) print *,"PENTRC - euler-lagrange matrix calculation"
            !! HACK - this should have its own flag
            allocate(wtw(mpert,mpert,6))
            write(nstring,'(I4)') nn
            open(unit=1,file="pentrc_tgar_elmat_n"//trim(adjustl(nstring))//".out",status="new")
            write(1,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
            write(1,*) "Kinetic additions to the ideal Euler-Lagrange matrices"
            write(1,'(1/,4x,2(a4,i4),1/)') "n = ",nn," l = ",0
            
            do i=1,nout
                if(verbose) print *," psi = ",psi_out(i)
                write(1,'(1/,1x,a16,es16.8e3,1/)') "psi = ",psi_out(i)
                write(1,'(1x,a4,a4,12(1x,a16))') "m_1","m_2","real(A_k)",&
                 "imag(A_k)","real(B_k)","imag(B_k)","real(C_k)","imag(C_k)",&
                 "real(D_k)","imag(D_k)","real(E_k)","imag(E_k)","real(H_k)",&
                 "imag(H_k)"
                do l=0,0!! should be all
                    if(psi_out(i)>0 .and. psi_out(i)<=1)then
                        tsurf = tpsi(psi_out(i),nn,l,zi,mi,wdfac,divxfac,electron,'twmm',&
                            op_wmats=wtw)
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
            if(flags(m))then
                method = methods(m) !to_upper(methods(m))
                if(verbose)then
                    print *, "---------------------------------------------"
                    print *,method//" - "//TRIM(docs(m))
                ENDIF
                if ((method=="clar" .or. method=="rlar")) then ! .and. fnml%nqty==0) then
                    call read_fnml(TRIM(data_dir)//'/fkmnl.dat')
                endif
                tphi = tintgrl_lsode(psilim,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                if(verbose) then
                    print "(a24,es11.3E3)", "Total torque = ", real(tphi)
                    print "(a24,es11.3E3)", "Total Kinetic Energy = ", aimag(tphi)/(2*nn)
                    print "(a24,es11.3E3)", "alpha/s  = ", real(tphi)/(-1*aimag(tphi))
                endif
                ierr=set_harvest_payload_dbl(hlog,'torque_'//method//nul,real(tphi))
                ierr=set_harvest_payload_dbl(hlog,'deltaW_'//method//nul,aimag(tphi)/(2*nn))
                if(equil_grid)then
                    if(verbose) print *,method//" - "//"Recalculating on equilibrium grid"
                    teq = tintgrl_grid('equil',psilim,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                    if(verbose)then
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total torque = ", REAL(teq),&
                            ", % error = ",ABS(REAL(teq)-REAL(tphi))/REAL(tphi)
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn),&
                            ", % error = ",ABS(AIMAG(teq)-AIMAG(tphi))/AIMAG(tphi)
                        print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                    endif
                endif
                if(input_grid)then
                    if(verbose) print *,method//" - "//"Recalculating on input displacements' grid"
                    teq = tintgrl_grid('input',psilim,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                    if(verbose)then
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total torque = ", REAL(teq),&
                            ", % error = ",ABS(REAL(teq)-REAL(tphi))/REAL(tphi)
                        print "(a24,es11.3E3,a12,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn),&
                            ", % error = ",ABS(AIMAG(teq)-AIMAG(tphi))/AIMAG(tphi)
                        print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                    endif
                endif
                ! run select surfaces with detailed output
                if(theta_out .or. xlmda_out)then
                    if(verbose) print *,method//" - "//"Recalculating on psi_out grid for detailed outputs"
                    allocate(tfuns(ntheta*3,3,nout,nthetafuns))
                    do i=1,nout
                        do l=-nl,nl,max(1,nl)
                            if(psi_out(i)>0 .and. psi_out(i)<=1)then
                                tsurf = tpsi(psi_out(i),nn,l,zi,mi,wdfac,divxfac,electron,methods(m),&
                                             op_tfuns=tfuns(:,1+nl/l,i,:))
                            endif
                        enddo
                    enddo
                    if(output_ascii)then
                        if(theta_out) call output_bouncefun_ascii(nn,zi,mi,electron,methods(m),&
                            reshape(tfuns,(/nout*3*ntheta*3,nthetafuns/)))
                        if(xlmda_out) call output_pitch_record(nn,zi,mi,electron,methods(m))
                        if(xlmda_out) call output_energy_record(nn,zi,mi,electron,methods(m))
                    endif
                    deallocate(tfuns)
                endif
                if(verbose)then
                    print *, method//" - Finished"
                    print *, "---------------------------------------------"
                endif
            endif
        enddo

    endif

    ! send harvest record
    ierr=harvest_send(hlog)
    
    ! display timer and stop
    call timer(mode=1)
    stop "PENTRC STOP=> normal termination."
end program pentrc
