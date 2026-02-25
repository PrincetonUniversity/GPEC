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
    ! EMAIL: nikolas.logan@columbia.edu
    !-----------------------------------------------------------------------
    use pentrc_interface
    use utilities, only: append_1d, append_2d, progressbar

    ! local variables
    integer :: i,j,k,l,m,nvalid
    real(r8), dimension(:), allocatable ::  psi_out_valid


    ! harvest variables
    include 'harvest_lib.inc'
    integer :: ierr
    character(len=65507) :: hlog
    character, parameter :: nul = char(0)

    if(verbose) print *,''
    if(verbose) print *,"PENTRC START => "//trim(version)
    if(verbose) print *,"_____________________________________________"

    call initialize_pentrc

    if(moment=="heat")then
        qt = .true.
        if(verbose) print *,"---------------------------------------------"
        if(verbose) print *,"Calculating heat transport"
        if(verbose) print *,"---------------------------------------------"
    elseif(moment=="pressure")then
        qt = .false.
        if(verbose) print *,"---------------------------------------------"
        if(verbose) print *,"Calculating particle transport and torque"
        if(verbose) print *,"---------------------------------------------"
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

    ! set the number of threads for OPENMP
#ifdef _OPENMP
    IF(pentrc_threads > 0)THEN
      CALL OMP_SET_NUM_THREADS(pentrc_threads)
    ENDIF
#else
    if(pentrc_threads > 1) print *,"!! WARNING: Not compiled with OPENMP. Forcing pentrc_threads = 1."
    pentrc_threads = 1
#endif

    ! run models
        ! start log with harvest
        ierr=init_harvest('CODEDB_PENT'//nul,hlog,len(hlog))
        ierr=set_harvest_verbose(0)
        ! standard CODEDB records
        ierr=set_harvest_payload_str(hlog,'CODE'//nul,'PENT'//nul)
        ierr=set_harvest_payload_str(hlog,'VERSION'//nul,version//nul)
        ! record PENT input
        ierr=set_harvest_payload_int(hlog,'zi'//nul,zi)
        ierr=set_harvest_payload_int(hlog,'zimp'//nul,zimp)
        ierr=set_harvest_payload_int(hlog,'mi'//nul,mi)
        ierr=set_harvest_payload_int(hlog,'mimp'//nul,mimp)
        ierr=set_harvest_payload_int(hlog,'mimp'//nul,mimp)
        ierr=set_harvest_payload_bol(hlog,'electron'//nul,electron)
        ierr=set_harvest_payload_str(hlog,'nutype'//nul,trim(nutype)//nul)
        ierr=set_harvest_payload_str(hlog,'f0type'//nul,trim(f0type)//nul)

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
            
            do i=1,npsi_out
                if(verbose) print *," psi = ",psi_out(i)
                write(1,'(1/,1x,a16,es16.8e3,1/)') "psi = ",psi_out(i)
                write(1,'(1x,a4,a4,12(1x,a16))') "m_1","m_2","real(A_k)",&
                 "imag(A_k)","real(B_k)","imag(B_k)","real(C_k)","imag(C_k)",&
                 "real(D_k)","imag(D_k)","real(E_k)","imag(E_k)","real(H_k)",&
                 "imag(H_k)"
                do l=0,0!! should be all
                    if(psi_out(i)>0 .and. psi_out(i)<=1)then
                        call tpsi(tsurf,psi_out(i),nn,l,zi,mi,wdfac,divxfac,electron,'twmm',&
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

        do m=1,nmethods
            if(flags(m))then
                method = methods(m) !to_upper(methods(m))
                if(verbose)then
                    print *, "---------------------------------------------"
                    print *,method//" - "//TRIM(docs(m))
                ENDIF
                if ((method=="clar" .or. method=="rlar")) then ! .and. fnml%nqty==0) then
                    if (trim(data_dir)=='' .or. trim(data_dir)=='default') then
                        call getenv('GPECHOME',data_dir)
                        if(len(trim(data_dir))==0) stop &
                            "ERROR: Use of fefault data directory requires GPECHOME environment variable"
                        data_dir = TRIM(data_dir)//'/pentrc'
                    endif
                    call read_fnml(trim(data_dir)//'/fkmnl.dat')
                endif
                if(dynamic_grid)then
                    if(verbose) print *,method//" - "//"Calculating using dynamic integration"
                    tphi = tintgrl_lsode(psilims,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                    if(verbose) then
                        print "(a24,es11.3E3)", "Total torque = ", real(tphi)
                        print "(a24,es11.3E3)", "Total Kinetic Energy = ", aimag(tphi)/(2*nn)
                        print "(a24,es11.3E3)", "alpha/s  = ", real(tphi)/(-1*aimag(tphi))
                    endif
                    ierr=set_harvest_payload_dbl(hlog,'torque_'//method//nul,real(tphi))
                    ierr=set_harvest_payload_dbl(hlog,'deltaW_'//method//nul,aimag(tphi)/(2*nn))
                end if
                if(equil_grid)then
                    if(verbose) print *,method//" - "//"Calculating on equilibrium grid"
                    teq = tintgrl_grid('equil',psilims,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                    if(verbose)then
                        if(dynamic_grid)then
                            print "(a24,es11.3E3,a12,es11.3E3)", "Total torque = ", REAL(teq),&
                                ", % error = ",ABS(REAL(teq)-REAL(tphi))/REAL(tphi)
                            print "(a24,es11.3E3,a12,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn),&
                                ", % error = ",ABS(AIMAG(teq)-AIMAG(tphi))/AIMAG(tphi)
                            print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                        else
                            print "(a24,es11.3E3)", "Total torque = ", REAL(teq)
                            print "(a24,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn)
                            print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                        end if
                    endif
                endif
                if(input_grid)then
                    if(verbose) print *,method//" - "//"Calculating on input displacements' grid"
                    teq = tintgrl_grid('input',psilims,nn,nl,zi,mi,wdfac,divxfac,electron,methods(m))
                    if(verbose)then
                        if(dynamic_grid)then
                            print "(a24,es11.3E3,a12,es11.3E3)", "Total torque = ", REAL(teq),&
                                ", % error = ",ABS(REAL(teq)-REAL(tphi))/REAL(tphi)
                            print "(a24,es11.3E3,a12,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn),&
                                ", % error = ",ABS(AIMAG(teq)-AIMAG(tphi))/AIMAG(tphi)
                            print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                        else
                            print "(a24,es11.3E3)", "Total torque = ", REAL(teq)
                            print "(a24,es11.3E3)", "Total Kinetic Energy = ", AIMAG(teq)/(2*nn)
                            print "(a24,es11.3E3)", "alpha/s  = ", REAL(teq)/(-1*AIMAG(teq))
                        end if
                    endif
                endif
                ! run select surfaces with detailed output
                if(theta_out .or. xlmda_out)then
                    if(verbose) print *,method//" - "//"Recalculating on psi_out grid for detailed outputs"
                    ! only use valid output surfaces
                    do i=1,npsi_out
                        if(psi_out(i)>0 .and. psi_out(i)<1) call append_1d(psi_out_valid,psi_out(i))
                    enddo
                    if(allocated(psi_out_valid))then
                        nvalid = size(psi_out_valid,dim=1)
                        do i=1,nvalid
                            if(nvalid<10) print '(a8,es11.3E3)',"  psi = ",psi_out_valid(i)
                            do l=-nl,nl,max(1,nl)
                                call tpsi(tsurf,psi_out_valid(i),nn,l,zi,mi,wdfac,divxfac,electron,methods(m),&
                                             op_erecord=xlmda_out,op_orecord=theta_out)
                            enddo
                            if(nvalid>10) call progressbar(i,1,nvalid,op_percent=20)
                        enddo
                        deallocate(psi_out_valid)
                    endif
                endif
                if(verbose)then
                    print *, method//" - Finished"
                    print *, "---------------------------------------------"
                endif
            endif
        enddo
        if(output_ascii)then
            if(theta_out) call output_orbit_ascii(nn)
            if(xlmda_out) call output_pitch_ascii(nn)
            if(xlmda_out) call output_energy_ascii(nn)
        endif
        if(output_netcdf)then
            call output_torque_netcdf(nn,nl,zi,mi,electron,wdfac)
            if(theta_out) call output_orbit_netcdf(nn)
            if(xlmda_out) call output_pitch_netcdf(nn)
            if(xlmda_out) call output_energy_netcdf(nn)
        endif
    endif

    ! send harvest record
    ierr=harvest_send(hlog)
    
    ! display timer and stop
    if(verbose)then
        call timer(mode=1)
        stop " PENTRC STOP => normal termination."
    else
        stop
    end if
end program pentrc
