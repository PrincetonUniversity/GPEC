!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module inputs
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Input interface for PENTRC. Includes DCON interface developed by
    !   J.-K. Park in IPEC package, as well as ascii interfaces for
    !   kinetic profiles and perturbed equilibrium files.
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !   read_profiles           - Read ascii kinetic profiles, form kin spline
    !   read_perteq             - Read ascii perturbartions, form xs_m spline
    !   read_dcon               - Read dcon bin files, form idcon splines
    !
    !*PUBLIC DATA MEMBERS:
    !   ro, real                - Major radius (idcon)
    !   bo, real                - Field on axis (determined from sq)
    !   kin, spline_type        - (psi) Kinetic profiles
    !   xs_m, fspline_type     - (psi,theta) First order perturbations
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
    
    use params, only : r8,xj,mp,me,e,mu0,twopi
    use utilities, only : get_free_file_unit,readtable,nunique,progressbar,iscdftf,iscdftb
    use spline_mod, only : spline_type,spline_alloc,spline_fit,spline_eval,spline_write1
    use cspline_mod, only : cspline_type,cspline_alloc,cspline_dealloc,&
                            cspline_fit,cspline_write,cspline_eval
    use fspline_mod, only : fspline_eval
    use bicube_mod, only : bicube_type,bicube_alloc,bicube_fit,bicube_eval
    
    use dcon_interface, only : idcon_read,idcon_transform,idcon_metric,&
        idcon_action_matrices,idcon_build,set_geom,idcon_harvest,&
        geom,eqfun,sq,rzphi,smats,tmats,xmats,ymats,zmats,&
        chi1,ro,zo,bo,nn,idconfile,jac_type,&
        mfac,psifac,mpert,mstep,mthsurf,theta,&
        idcon_coords

    implicit none

    private
    public &
         read_kin, &
         read_pmodb, &
         read_peq, &
         set_peq, &
         read_ipec_peq,&
         read_equil, &
         read_fnml, &
         kin, xs_m, dbob_m, divx_m, fnml, &
         chi1,ro,zo,bo,nn,mfac,mpert,mthsurf, &
         verbose

    ! global variables with defaults
    logical :: verbose=.TRUE.
    type(spline_type) :: kin
    type(cspline_type) :: dbob_m, divx_m, xs_m(3)
    type(bicube_type):: fnml

!$OMP THREADPRIVATE(dbob_m, divx_m, kin, xs_m, fnml)

    contains

    !=======================================================================
    function newm(len_i,m_i,fm_i,len_o,m_o)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Simply transfer the fm spectrum on m_i to the new modes m_o.
    !
    !*ARGUMENTS:
    !    len_i : integer (in)
    !       Length of initial spectrum array.
    !    m_i : integer array (in)
    !       Initial modes.
    !    fm_i : complex array (in)
    !       Initial spectrum array.
    !    len_o : integer (in)
    !       Length of output spectrum array.
    !    m_o : integer array (in)
    !       Output modes.
    !
    !-----------------------------------------------------------------------

        implicit none
        ! declare arguments
        logical :: debug = .false.
        integer :: len_i,len_o,i,ii,j,jj
        integer, dimension(len_i) :: m_i
        integer, dimension(len_o) :: m_o
        complex(r8), dimension(len_i) :: fm_i
        complex(r8), dimension(len_o) :: newm

        ! index the new range in the old range
        i = MAX(ABS(m_i(1))-ABS(m_o(1))+1,1)
        ii = len_i - MAX(ABS(m_i(len_i))-ABS(m_o(len_o)),0)
        j = MAX(ABS(m_o(1))-ABS(m_i(1))+1,1)
        jj = len_o - MAX(ABS(m_o(len_o))-ABS(m_i(len_i)),0)

        ! check alignment
        if(debug)then
            print *,'  -> Converting to DCON m-space',m_o(1),'to',m_o(len_o),&
                'filling',m_o(j),'to',m_o(jj),'from inputs'
        endif
        if(m_o(j)/=m_i(i) .or. m_o(jj)/=m_i(ii)) stop "ERROR: Misalignment in m-spaces"

        ! transfer
        newm = 0
        newm(j:jj) = fm_i(i:ii)

    end function newm


    !=======================================================================
    subroutine read_equil(file,hlog)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read dcon binary and form all the equilibrium splines.
    !
    !*ARGUMENTS:
    !    file : character(512) (in)
    !       File path.
    !
    !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        character(*), intent(in) :: file
        character(len=65507), optional :: hlog

        ! set idconfile
        idconfile = file
        ! prepare ideal solutions. (psixy=0)
        CALL idcon_read(0)
        if(present(hlog)) CALL idcon_harvest(hlog)
        CALL idcon_transform
        ! reconstruct metric tensors.
        CALL idcon_metric
        ! read vacuum data.
        !CALL idcon_vacuum
        ! form action matrices (this will move !!)
        call idcon_action_matrices!(egnum,xspmn)

        ! set additional geometry spline
        call set_geom

    end subroutine read_equil
    
    
    !=======================================================================
    subroutine read_kin(file,zi,zimp,mi,mimp,nfac,tfac,wefac,wpfac,write_log)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read ascii file containing table of kinetic profiles ni, ne, ti,
    !   te, and omegaE, then form kin spline containing some additional
    !   information (krook nui,nue).
    !
    !   Assumes table consists of 6 columns: psi_n, n_i(m^-3), n_e(m^-3),
    !   T_i(eV), T_e(eV), omega_E(rad/s). File can (nearly) arbitrary
    !   header and/or footer, with the exception that no lines start with
    !   a number.
    !
    !*ARGUMENTS:
    !    file : character(256) (in)
    !       File path.
    !   zi : integer
    !       Ion charge in fundemental units
    !   zimp : integer
    !       Impurity ion charge in fundemental units
    !   mi : integer
    !       Ion mass in fundemental units (mass proton)
    !   mimp : integer
    !       Impurity ion mass in fundemental units
    !   wefac : real
    !       Direct multiplier for omegaE profiles
    !   wpfac : real
    !       Scaling of rotation profile, done via manipulation of omegaE
    !   write_log : bool
    !       Writes kinetic spline to log file
    !   
    !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        character(512), intent(in) :: file
        integer, intent(in) :: zi,zimp,mi,mimp
        real(r8), intent(in) :: wefac,wpfac,nfac,tfac
        logical :: write_log
        ! declare local variables
        real(r8), dimension(:,:), allocatable :: table
        character(32), dimension(:), allocatable :: titles
        
        integer, parameter :: nkin = 100
        integer :: i,out_unit, tshape(2)
        real(r8) :: psi
        real(r8), dimension(0:nkin) :: zeff,zpitch,welec,wdian,wdiat,wphi,wpefac
        type(spline_type) :: tmp
        
        call readtable(file,table,titles,verbose,write_log)
        tshape = shape(table)
        if(write_log) print *,"Table shaped ",tshape
        
        ! temp spline on experimental grid
        call spline_alloc(tmp,tshape(1)-1,5)
        tmp%xs(0:) = table(1:,1)
        do i=1,5
            tmp%fs(0:,i) = table(1:,i+1)
        enddo
        call spline_fit(tmp,"extrap")
        if(write_log) print *,"Formed temporary spline"
        
        ! extrapolate to regular spline (helps smooth core??)
        call spline_alloc(kin,nkin,8)
        kin%title(0:) = (/"psi_n ","n_i   ","n_e   ","t_i   ","t_e   ","omegae",&
                    "loglam","nu_i  ","nu_e  " /)
        do i=0,kin%mx
              psi = (1.0*i)/kin%mx
              call spline_eval(tmp,psi,0)
              kin%xs(i) = psi
              kin%fs(i,1:5) = tmp%f(1:5)
              !print *,i,kin%fs(i,1:4)
        enddo
        
        ! convert temperatures to si units
        kin%fs(:,3:4) = kin%fs(:,3:4)*1.602e-19 !ev to j

        ! collisionalities **assumes si units for n,t**
        zeff = zimp-(kin%fs(:,1)/kin%fs(:,2))*zi*(zimp-zi)
        zpitch = 1.0+(1.0+mimp)/(2.0*mimp)*zimp*(zeff-1.0)/(zimp-zeff)   
        kin%fs(:,6) = 17.3-0.5*log(kin%fs(:,2)/1.0e20) &
            +1.5*log(kin%fs(:,4)/1.602e-16)
        kin%fs(:,7) = (zpitch/3.5e17)*kin%fs(:,1)*kin%fs(:,6) &
            /(sqrt(1.0*mi)*(kin%fs(:,3)/1.602e-16)**1.5)
        kin%fs(:,8) = (zpitch/3.5e17)*kin%fs(:,2)*kin%fs(:,6) &
            /(sqrt(me/mp)*(kin%fs(:,4)/1.602e-16)**1.5)
        call spline_fit(kin,"extrap")
        if(write_log) print *,"Formed kin spline"
    
        ! manipulation of N and T profiles
        kin%fs(:,1:2) = nfac*kin%fs(:,1:2)
        kin%fs(:,3:4) = tfac*kin%fs(:,3:4)
    
        ! manipulation of rotation variables
        welec(:) = wefac*kin%fs(:,5) ! direct manipulation of omegae
        if(wefac/=1.0 .and. verbose) print('(a55,es10.2e3)'),'  -> applying direct manipulation of omegaE by factor ',wefac
        wdian =-twopi*kin%fs(:,3)*kin%fs1(:,1)/(e*zi*chi1*kin%fs(:,1))
        wdiat =-twopi*kin%fs1(:,3)/(e*zi*chi1)
        wpefac= (wpfac*(welec+wdian+wdiat) - (wdian+wdiat))/welec
        welec = wpefac*welec   ! indirect manipulation of rotation
        if(wpfac/=1.0 .and. verbose) then
            print('(a40,es10.2e3)'),'  -> manipulating rotation by factor of ',wpfac
            print *,'     by indierct manipulation of omegae profile'
        endif
        kin%fs(:,5) = welec(:)
        call spline_fit(kin,"extrap")
        if(write_log) print *,"Reformed kin spline with rotation manipulatins"

        ! write log - designed as check of reading routines
        if(write_log)then
            out_unit = get_free_file_unit(-1)
            if(verbose) print *, "Writing kinetic spline to pentrc_kinetics.out"
            open(unit=out_unit,file="pentrc_kinetics.out",&
                status="unknown")!,action="write")
            call spline_write1(kin,.true.,.false.,out_unit,0,.true.)
            close(out_unit)
        endif
    
    end subroutine read_kin


    !=======================================================================
    subroutine read_pmodb(file,jac_in,jsurf_in,tmag_in,debug,op_powin)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Read psi,m matrix of displacements.
    !
    !*ARGUMENTS:
    !   file : character(512) (in)
    !       File path.
    !   jac_in : character
    !       Input file jacobian.
    !   jsurf_in : int.
    !       Surface weigted inputs should be 1
    !   tmag_in : int.
    !       Input toroidal angle specification: 1 = magnetic, 0 = cylindrical
    !   debug : logical
    !       Print intermidient messages to terminal.
    !
    !*OPTIONAL ARGUMENTS:
    !   op_powin : int(4)
    !       User specified powers of B, Bp, R, and Rc that define a Jacobian.
    !       Only used if jac_in is 'other'.
    !
    !-----------------------------------------------------------------------

        implicit none

        ! declare arguments
        logical, intent(in) :: debug
        integer, intent(in) :: jsurf_in,tmag_in
        integer, dimension(4), intent(in), optional :: op_powin
        character(32), intent(inout) :: jac_in
        character(512), intent(in) :: file
        ! declare local variables
        logical :: ncheck
        integer :: i,j,npsi,nm,ndigit,firstnm, powin(4), out_unit
        integer, dimension(:), allocatable :: ms
        real(r8), dimension(:), allocatable :: psi
        real(r8), dimension(:,:), allocatable :: table
        complex(r8), dimension(:,:), allocatable :: lagbmni,divxmni,kapxmni,&
            lagbmns,divxmns
        complex(r8), dimension(0:mthsurf) :: lagbfun,divxfun
        character(3) :: nstr
        character(32), dimension(:), allocatable :: titles

        ! file consistency check (requires naming convention)
        write(nstr,'(i3)') nn
        if(nn<-9)then
            ndigit = 3
        elseif(nn>=0 .and. nn<10)then
            ndigit = 1
        else ! assume only running with -99 to 99
            ndigit = 2
        endif
        ncheck = .false.
        DO i=1,512-ndigit
            if(file(i:i+ndigit)=='n'//trim(adjustl(nstr))) ncheck = .true.
        ENDDO
        if(.not. ncheck)then
            print *,"** Toroidal mode number determined from idconfile is",nn
            print *,"** Corresponding label 'n"//trim(adjustl(nstr))//"' must be in peq_file name"
            stop "ERROR: Inconsistent toroidal mode numbers"
        endif

        ! read file
        call readtable(file,table,titles,verbose,debug)
        ! should be npsi*nm by 8 (psi,m,realxi_1,imagxi_1,...)
        !npsi = nunique(table(:,1)) !! computationally expensive + ipec n=3's can have repeates
        nm = nunique(table(:,2),op_sorted=.True.)
        npsi = size(table,1)/nm
        if(npsi*nm/=size(table,1))then
            stop "ERROR - inputs - size of table not equal to product of unique psi & m"
        endif
        if(debug) print *,"  -> found ",npsi," steps in psi and ",nm," modes"

        allocate(ms(nm),psi(npsi))
        allocate(lagbmni(npsi,nm),divxmni(npsi,nm),kapxmni(npsi,nm))
        allocate(lagbmns(npsi,mpert),divxmns(npsi,mpert))
        firstnm = nunique(table(1:nm,2))
        if(firstnm==nm)then ! written with psi as outer loop
            ms = table(1:nm,2)
            psi = (/(table(j,1),j=1,npsi*nm,nm)/)
            if(debug) print *,"psi outerloop"
            lagbmni(:,:) = reshape(table(:,5),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,6),(/npsi,nm/),order=(/2,1/))
            divxmni(:,:) = reshape(table(:,7),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,8),(/npsi,nm/),order=(/2,1/))
            kapxmni(:,:) = reshape(table(:,9),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,10),(/npsi,nm/),order=(/2,1/))
        else ! written with m as outer loop
            ms = (/(table(i,2),i=1,npsi*nm,npsi)/)
            psi = table(1:npsi,1)
            if(debug) print *,"m outerloop"
            lagbmni(:,:) = reshape(table(:,5),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,6),(/npsi,nm/),order=(/1,2/))
            divxmni(:,:) = reshape(table(:,7),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,8),(/npsi,nm/),order=(/1,2/))
            kapxmni(:,:) = reshape(table(:,9),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,10),(/npsi,nm/),order=(/1,2/))
        endif

        ! For consistency with IPEC-0.3.0 and smaller values near rationals
        if(verbose) print *,"  -> Using div(xi_perp) = -(dB/B+kappa.xi_perp)"
        divxmni = -(lagbmni+kapxmni)

        ! convert to chebyshev coordinates
        if(jac_in=="" .or. jac_in=="default")then
            jac_in = jac_type
            if(verbose) print *,"  -> WARNING: Assuming DCON "//trim(jac_type)//" coordinates"
        endif
        powin = -1
        SELECT CASE(jac_in) ! set B,Bp,R, and Rc powers
            CASE("hamada")
                powin=(/0,0,0,0/)
            CASE("pest")
                powin=(/0,0,2,0/)
            CASE("equal_arc")
                powin=(/0,1,0,0/)
            CASE("boozer")
                powin=(/2,0,0,0/)
            CASE("polar")
                powin=(/0,1,0,1/)
            CASE("other")
                if(present(op_powin)) powin=op_powin
                if(powin(1)<0 .or. powin(2)<0 .or. powin(3)<0 .or. powin(4)<0)then
                    stop "ERROR: Did not receive all four powers for jac_in 'other'"
                endif
            CASE DEFAULT
                stop "ERROR: inputs - jac_in must be 'hamada','pest','equal_arc','boozer',&
                    & 'polar', or 'other'. Setting to 'default' uses idconfile jac_type."
        END SELECT
        if(verbose) print *,"  -> Displacements input in "//trim(jac_in)//" coordinates: b,bp,r,rc raised to",powin
        if(jac_in/=jac_type .or. tmag_in/=1)then
            if(tmag_in/=1 .and. verbose) print *,'     Displacements input in cylindrical toroidal angle'
            if(verbose) print *,'Converting to '//trim(jac_type)//' coordinates used by DCON'
            ! make sure to use the larger of the input and working spectra
            if(mpert>nm)then
                ! convert spectrum on each surface
                do i=1,npsi
                    if(verbose) call progressbar(i,1,npsi,op_percent=20)
                    lagbmns(i,:) = newm(nm,ms,lagbmni(i,:),mpert,mfac)
                    divxmns(i,:) = newm(nm,ms,divxmni(i,:),mpert,mfac)
                    CALL idcon_coords(psi(i),lagbmns(i,:),mfac,mpert,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),divxmns(i,:),mfac,mpert,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                enddo
            else
                ! convert spectrum on each surface
                do i=1,npsi
                    if(verbose) call progressbar(i,1,npsi,op_percent=20)
                    CALL idcon_coords(psi(i),lagbmni(i,:),ms,nm,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),divxmni(i,:),ms,nm,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    lagbmns(i,:) = newm(nm,ms,lagbmni(i,:),mpert,mfac)
                    divxmns(i,:) = newm(nm,ms,divxmni(i,:),mpert,mfac)
                enddo
            endif
        else
            do i=1,npsi
                lagbmns(i,:) = newm(nm,ms,lagbmni(i,:),mpert,mfac)
                divxmns(i,:) = newm(nm,ms,divxmni(i,:),mpert,mfac)
            enddo
        endif

        ! (re-)set global perturbation splines with 1/B weighting
        if(verbose) print *,"  -> Weighting by 1/B"
        if(associated(dbob_m%xs)) call cspline_dealloc(dbob_m)
        if(associated(divx_m%xs)) call cspline_dealloc(divx_m)
        call cspline_alloc(dbob_m,npsi-1,mpert)     ! (dB/B)
        call cspline_alloc(divx_m,npsi-1,mpert)     ! div(xi_prp)
        dbob_m%xs(0:) = psi(1:)
        divx_m%xs(0:) = psi(1:)
        do i=1,npsi
            if(verbose) call progressbar(i,1,npsi,op_step=1,op_percent=20)
            call iscdftb(mfac,mpert,lagbfun,mthsurf,lagbmns(i,:))
            call iscdftb(mfac,mpert,divxfun,mthsurf,divxmns(i,:))
            do j=0,mthsurf
                call bicube_eval(eqfun,psi(i),theta(j),0)
                lagbfun(j) = lagbfun(j) / eqfun%f(1) ! dB/B
                divxfun(j) = divxfun(j) / eqfun%f(1) ! nabla.xi_perp
            enddo
            call iscdftf(mfac,mpert,lagbfun,mthsurf,dbob_m%fs(i-1,:))
            call iscdftf(mfac,mpert,divxfun,mthsurf,divx_m%fs(i-1,:))
        enddo
        call cspline_fit(dbob_m,"extrap")
        call cspline_fit(divx_m,"extrap")

        ! write log to check reading/allocating routines
        if(debug)then
            print *,"  Writing log of pmodb splines"
            print *,"  inputs: mlow = ",mfac(1)," mhigh = ",mfac(mpert)," nm = ",nm
            print *,"  inputs: psilow = ",dbob_m%xs(0)," psihigh = ",dbob_m%xs(dbob_m%mx)," npsi = ",npsi
            out_unit = get_free_file_unit(-1)
            open(unit=out_unit,file="pentrc_pmodb_n"//trim(adjustl(nstr))//".out",&
                status="unknown")
            write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
            write(out_unit,'(a43,2/)') " Debuging log file for set_peq subroutine."
            write(out_unit,'(a8,es16.8e3,2/)') " chi1 = ",chi1
            write(out_unit,'(1x,a16,a5,8(1x,a16))')"psi_n","m",&
                'real(deltaB/B)','imag(deltaB/B)','real(divxprp)','imag(divxprp)',&
                'real(deltaB)','imag(deltaB)','real(Bdivxprp)','imag(Bdivxprp)'
            do i=0,npsi-1
                do j=1,mpert
                    write(out_unit,'(1x,es16.8e3,i5,14(1x,es16.8e3))') &
                        psi(i+1),mfac(j),&
                        dbob_m%fs(i,j),divx_m%fs(i,j),lagbmns(i,j),divxmns(i,j)
                enddo
            enddo
            close(out_unit)
        endif

        deallocate(ms,psi,lagbmni,divxmni,kapxmni,lagbmns,divxmns)

    end subroutine read_pmodb


    !=======================================================================
    subroutine read_peq(file,jac_in,jsurf_in,tmag_in,debug,op_powin)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Read psi,m matrix of displacements.
    !
    !*ARGUMENTS:
    !   file : character(512) (in)
    !       File path.
    !   jac_in : character
    !       Input file jacobian.
    !   jsurf_in : int.
    !       Surface weigted inputs should be 1
    !   tmag_in : int.
    !       Input toroidal angle specification: 1 = magnetic, 0 = cylindrical
    !   debug : logical
    !       Print intermidient messages to terminal.
    !
    !*OPTIONAL ARGUMENTS:
    !   op_powin : int(4)
    !       User specified powers of B, Bp, R, and Rc that define a Jacobian.
    !       Only used if jac_in is 'other'.
    !
    !-----------------------------------------------------------------------

        implicit none

        ! declare arguments
        logical, intent(in) :: debug
        integer, intent(in) :: jsurf_in,tmag_in
        integer, dimension(4), intent(in), optional :: op_powin
        character(32), intent(inout) :: jac_in
        character(512), intent(in) :: file
        ! declare local variables
        logical :: ncheck
        integer :: i,j,npsi,nm,ndigit,firstnm, powin(4)
        integer, dimension(:), allocatable :: ms
        real(r8), dimension(:), allocatable :: psi
        real(r8), dimension(:,:), allocatable :: table
        complex(r8), dimension(:,:), allocatable :: xmp1mns,xspmns,xmsmns,xmp1mni,xspmni,xmsmni
        character(3) :: nstr
        character(32), dimension(:), allocatable :: titles

        ! file consistency check (requires naming convention)
        write(nstr,'(i3)') nn
        if(nn<-9)then
            ndigit = 3
        elseif(nn>=0 .and. nn<10)then
            ndigit = 1
        else ! assume only running with -99 to 99
            ndigit = 2
        endif
        ncheck = .false.
        DO i=1,512-ndigit
            if(file(i:i+ndigit)=='n'//trim(adjustl(nstr))) ncheck = .true.
        ENDDO
        if(.not. ncheck)then
            print *,"** Toroidal mode number determined from idconfile is",nn
            print *,"** Corresponding label 'n"//trim(adjustl(nstr))//"' must be in peq_file name"
            stop "ERROR: Inconsistent toroidal mode numbers"
        endif

        ! read file
        call readtable(file,table,titles,verbose,debug)
        ! should be npsi*nm by 8 (psi,m,realxi_1,imagxi_1,...)
        !npsi = nunique(table(:,1)) !! computationally expensive + ipec n=3's can have repeats
        nm = nunique(table(:,2),op_sorted=.True.)
        npsi = size(table,1)/nm
        if(npsi*nm/=size(table,1))then
            stop "ERROR - inputs - size of table not equal to product of unique psi & m"
        endif
        if(debug) print *,"  -> found ",npsi," steps in psi and ",nm," modes"

        allocate(ms(nm),psi(npsi))
        allocate(xmp1mni(npsi,nm),xspmni(npsi,nm),xmsmni(npsi,nm))
        allocate(xmp1mns(npsi,mpert),xspmns(npsi,mpert),xmsmns(npsi,mpert))
        firstnm = nunique(table(1:nm,2))
        if(firstnm==nm)then ! written with psi as outer loop
            ms = table(1:nm,2)
            psi = (/(table(j,1),j=1,npsi*nm,nm)/)
            if(debug) print *,"psi outerloop"
            !if(debug) print *,"ms = ",ms
            !if(debug) print *,"psi range = ",psi(1),psi(npsi)
            xmp1mni(:,:) = reshape(table(:,3),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,4),(/npsi,nm/),order=(/2,1/))
            xspmni(:,:) = reshape(table(:,5),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,6),(/npsi,nm/),order=(/2,1/))
            xmsmni(:,:) = reshape(table(:,7),(/npsi,nm/),order=(/2,1/))&
                    +xj*reshape(table(:,8),(/npsi,nm/),order=(/2,1/))
        else ! written with m as outer loop
            ms = (/(table(i,2),i=1,npsi*nm,npsi)/)
            psi = table(1:npsi,1)
            if(debug) print *,"m outerloop"
            !if(debug) print *,"ms = ",ms
            !if(debug) print *,"psi range = ",psi(1),psi(npsi)
            xmp1mni(:,:) = reshape(table(:,3),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,4),(/npsi,nm/),order=(/1,2/))
            xspmni(:,:) = reshape(table(:,5),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,6),(/npsi,nm/),order=(/1,2/))
            xmsmni(:,:) = reshape(table(:,7),(/npsi,nm/),order=(/1,2/))&
                    +xj*reshape(table(:,8),(/npsi,nm/),order=(/1,2/))
        endif
        
        ! convert to chebyshev coordinates
        if(jac_in=="" .or. jac_in=="default")then
            jac_in = jac_type
            if(verbose) print *,"  -> WARNING: Assuming DCON "//trim(jac_type)//" coordinates"
        endif
        powin = -1
        SELECT CASE(jac_in) ! set B,Bp,R, and Rc powers
            CASE("hamada")
                powin=(/0,0,0,0/)
            CASE("pest")
                powin=(/0,0,2,0/)
            CASE("equal_arc")
                powin=(/0,1,0,0/)
            CASE("boozer")
                powin=(/2,0,0,0/)
            CASE("polar")
                powin=(/0,1,0,1/)
            CASE("other")
                if(present(op_powin)) powin=op_powin
                if(powin(1)<0 .or. powin(2)<0 .or. powin(3)<0 .or. powin(4)<0)then
                    stop "ERROR: Did not receive all four powers for jac_in 'other'"
                endif
            CASE DEFAULT
                stop "ERROR: inputs - jac_in must be 'hamada','pest','equal_arc','boozer',&
                    & 'polar', or 'other'. Setting to 'default' uses idconfile jac_type."
        END SELECT
        if(verbose) print *,"  -> Displacements input in "//trim(jac_in)//" coordinates: b,bp,r,rc raised to",powin
        if(jac_in/=jac_type .or. tmag_in/=1)then
            if(tmag_in/=1 .and. verbose) print *,'     Displacements input in cylindrical toroidal angle'
            if(verbose) print *,'Converting to '//trim(jac_type)//' coordinates used by DCON'
            ! make sure to use the larger of the input and working spectra
            if(mpert>nm)then
                ! convert spectrum on each surface
                do i=1,npsi
                    if(verbose) call progressbar(i,1,npsi,op_percent=20)
                    xmp1mns(i,:) = newm(nm,ms,xmp1mni(i,:),mpert,mfac)
                    xspmns(i,:) = newm(nm,ms,xspmni(i,:),mpert,mfac)
                    xmsmns(i,:) = newm(nm,ms,xmsmni(i,:),mpert,mfac)
                    CALL idcon_coords(psi(i),xmp1mns(i,:),mfac,mpert,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),xspmns(i,:),mfac,mpert,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),xmsmns(i,:),mfac,mpert,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                enddo
            else
                ! convert spectrum on each surface
                do i=1,npsi
                    if(verbose) call progressbar(i,1,npsi,op_percent=20)
                    CALL idcon_coords(psi(i),xmp1mns(i,:),ms,nm,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),xspmns(i,:),ms,nm,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    CALL idcon_coords(psi(i),xmsmns(i,:),ms,nm,&
                        powin(3),powin(2),powin(1),powin(4),tmag_in,jsurf_in)
                    xmp1mns(i,:) = newm(nm,ms,xmp1mni(i,:),mpert,mfac)
                    xspmns(i,:) = newm(nm,ms,xspmni(i,:),mpert,mfac)
                    xmsmns(i,:) = newm(nm,ms,xmsmni(i,:),mpert,mfac)
                enddo
            endif
        else
            do i=1,npsi
                xmp1mns(i,:) = newm(nm,ms,xmp1mni(i,:),mpert,mfac)
                xspmns(i,:) = newm(nm,ms,xspmni(i,:),mpert,mfac)
                xmsmns(i,:) = newm(nm,ms,xmsmni(i,:),mpert,mfac)
            enddo
        endif
        
        ! set global variables (perturbed quantity csplines)
        call set_peq(psi,mfac,xmp1mns,xspmns,xmsmns,.true.,debug)
        
        deallocate(ms,psi,xmp1mns,xspmns,xmsmns)
        
    end subroutine read_peq
    
    
    !=======================================================================
    subroutine set_peq(psi,ms,xmp1mns,xspmns,xmsmns,op_set_dbdx,op_debug)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Set m quantity xs_m complex splines in psi for the 2 Clebsch
    !   displacement components from 2D (psi,m) distributions.
    !   **This routine sets global variables xs_m dbob_m and divx_m**
    !
    !*ARGUMENTS:
    !    psi : real rank 1
    !       Normalized flux axis.
    !   ms : real rank 1
    !       Poloidal mode axis.
    !   xspmns : real rank 2
    !       The contravarient clebsch displacement xi^psi on (psi,m)
    !   xmsmns : real rank 2
    !       The contravarient clebsch displacement xi^alpha on (psi,m)
    !       *** NOTE this is NOT DCON CONVENTION -> xmsmns/chi1 in DCON ****
    !*OPTIONAL ARGUMENTS:
    !   op_set_dbdx : logical
    !       Calculate lagrangian mod-B and divergence of xi_perp and set
    !       them to global complex spline variable dbob_m,divx_m. (Defaut True)
    !   op_debug : logical
    !       Print intermidient messages to terminal. (Default False)
    !       
    !-----------------------------------------------------------------------
    
        implicit none
        
        ! declare arguments
        logical, intent(in), optional :: op_debug,op_set_dbdx
        integer,  dimension(:), intent(in) :: ms        
        real(r8), dimension(:), intent(in) :: psi
        complex(r8), dimension(:,:), intent(in) :: xmp1mns,xspmns,xmsmns
        ! declare local variables
        logical :: debug,set_dbdx
        integer :: i,j,ims,istrt_psi,istop_psi,npsi,nm, out_unit
        real(r8) :: r_mjr,r_mnr,jac,g12,g13,g22,g23,g33,gfac
        complex(r8), dimension(0:mthsurf) :: divxfun,dbobfun
        complex(r8), dimension(mpert,mpert) :: smat,tmat,xmat,ymat,zmat
        complex(r8), dimension(:,:), allocatable :: jbbkapxmns,jbbdivxmns,jbbdbobmns
        character(3) :: istring
        
        ! defaults for optional args
        debug = .false.
        set_dbdx = .true.
        if(present(op_debug)) debug = op_debug
        if(present(op_set_dbdx)) set_dbdx = op_set_dbdx
        ! get dimensions
        istrt_psi = LBOUND(psi,1)
        istop_psi = UBOUND(psi,1)
        npsi = size(psi)
        nm = size(ms)
        
        ! set global variables
        if(debug) print *,"  Setting global perturbed xi variables"
        !if(debug) print *,"  nm,mpert",nm,mpert
        if(debug) print *,"  mfac = ",mfac
        if(debug) print *,"  ms   = ",ms
        if(debug) print *,"  npsi = ",npsi
        if(debug) print *,"  psi lim = ",psi(1),psi(npsi)

        ! fill dcon m modes from available input
        do i=1,3
            if(associated(xs_m(i)%xs)) call cspline_dealloc(xs_m(i))
            call cspline_alloc(xs_m(i),npsi-1,mpert)
            xs_m(i)%xs(0:) = psi(:)
        enddo

        do i =1,mpert
            ims = i+(mfac(1)-ms(1))
            if(ims>0 .and. ims<=nm)then        
                xs_m(1)%fs(0:,i) =xmp1mns(1:,ims)
                xs_m(2)%fs(0:,i) = xspmns(1:,ims)
                xs_m(3)%fs(0:,i) = xmsmns(1:,ims)
            else
                print *,"WARNING: Not input for DCON m ",mfac(i)
                xs_m(1)%fs(0:,i) = 0
                xs_m(2)%fs(0:,i) = 0
                xs_m(3)%fs(0:,i) = 0
            endif
        enddo
        call cspline_fit(xs_m(1),"extrap")
        call cspline_fit(xs_m(2),"extrap")
        call cspline_fit(xs_m(3),"extrap")

        ! calculate db/b and divx to plug into old action integral formulation
        allocate(jbbkapxmns(npsi,mpert),jbbdivxmns(npsi,mpert),jbbdbobmns(npsi,mpert))
        if(associated(dbob_m%xs)) call cspline_dealloc(dbob_m)
        if(associated(divx_m%xs)) call cspline_dealloc(divx_m)
        call cspline_alloc(dbob_m,npsi-1,mpert)     ! dB/B
        call cspline_alloc(divx_m,npsi-1,mpert)     ! nabla.xi_perp

        dbob_m%xs(0:) = psi(1:)
        divx_m%xs(0:) = psi(1:)
        if(set_dbdx)then
            if(verbose) print *,'Calculating dB/B, div(xi_prp)'
            !call ipeq_alloc
            do i=istrt_psi,istop_psi
                j = i-istrt_psi+1
                if(verbose) call progressbar(j,1,npsi,op_percent=20)
                call cspline_eval(xs_m(1),psi(i),0)
                call cspline_eval(xs_m(2),psi(i),0)
                call cspline_eval(xs_m(3),psi(i),0)
                CALL cspline_eval(smats,psi(i),0)
                CALL cspline_eval(tmats,psi(i),0)
                CALL cspline_eval(xmats,psi(i),0)
                CALL cspline_eval(ymats,psi(i),0)
                CALL cspline_eval(zmats,psi(i),0)
                smat=RESHAPE(smats%f,(/mpert,mpert/))
                tmat=RESHAPE(tmats%f,(/mpert,mpert/))
                xmat=RESHAPE(xmats%f,(/mpert,mpert/))
                ymat=RESHAPE(ymats%f,(/mpert,mpert/))
                zmat=RESHAPE(zmats%f,(/mpert,mpert/))
                ! matrices give JBB weighted values
                jbbkapxmns(i,:) = MATMUL(smat,xs_m(2)%f(:))+MATMUL(tmat,xs_m(3)%f(:))
                jbbdivxmns(i,:) = MATMUL(xmat,xs_m(1)%f(:))+MATMUL(ymat,xs_m(2)%f(:)) &
                                +MATMUL(zmat,xs_m(3)%f(:))
                jbbdbobmns(i,:) = -(jbbdivxmns(i,:)+jbbkapxmns(i,:))
                ! remove weighting to get 1st order quantities
                call iscdftb(mfac,mpert,divxfun,mthsurf,jbbdivxmns(i,:))
                call iscdftb(mfac,mpert,dbobfun,mthsurf,jbbdbobmns(i,:))
                do j=0,mthsurf
                    call bicube_eval(eqfun,psi(i),theta(j),0)
                    call bicube_eval(rzphi,psi(i),theta(j),0)
                    divxfun(j) = divxfun(j) / (rzphi%f(4) * eqfun%f(1)**2) ! nabla.xi_perp
                    dbobfun(j) = dbobfun(j) / (rzphi%f(4) * eqfun%f(1)**2) ! dB/B
                enddo
                call iscdftf(mfac,mpert,divxfun,mthsurf,divx_m%fs(i-1,:))
                call iscdftf(mfac,mpert,dbobfun,mthsurf,dbob_m%fs(i-1,:))
            enddo
        else
            ! default dbdx is a flat spectrum
            if(verbose) print *,"Forming constant dB/B, div(xi_perp)"
            dbob_m%fs(:,:) = 1.0/sqrt(1.0*mpert)
            divx_m%fs(:,:) = 1.0/sqrt(1.0*mpert)
        endif
        call cspline_fit(dbob_m,"extrap")
        call cspline_fit(divx_m,"extrap")

        ! write log to check reading/allocating routines
        if(debug)then
            print *,"  inputs: mlow = ",mfac(1)," mhigh = ",mfac(mpert)," nm = ",nm
            print *,"  inputs: psilow = ",xs_m(1)%xs(0)," psihigh = ",xs_m(1)%xs(xs_m(1)%mx)," npsi = ",npsi
            out_unit = get_free_file_unit(-1)
            write(istring,'(i3)') nn
            open(unit=out_unit,file="pentrc_peq_n"//trim(adjustl(istring))//".out",&
                status="unknown")
            write(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE:"
            write(out_unit,'(a43,2/)') " Debugging log file for set_peq subroutine."
            write(out_unit,'(a8,es16.8e3,2/)') " chi1 = ",chi1
            write(out_unit,'(1x,a16,a5,14(1x,a16))')"psi_n","m",&
                "real(xi^psi1)","imag(xi^psi1)","real(xi^psi)","imag(xi^psi)","real(xi^alpha)","imag(xi^alpha)",&
                'real(deltaB/B)','imag(deltaB/B)','real(divxprp)','imag(divxprp)',&
                'real(JBdeltaB)','imag(JBdeltaB)','real(JBBdivxprp)','imag(JBBdivxprp)'
            do i=0,npsi-1
                do j=1,mpert
                    write(out_unit,'(1x,es16.8e3,i5,14(1x,es16.8e3))') &
                        xs_m(1)%xs(i),mfac(j),&
                        xs_m(1)%fs(i,j),xs_m(2)%fs(i,j),xs_m(3)%fs(i,j),&
                        dbob_m%fs(i,j),divx_m%fs(i,j),jbbdbobmns(i+1,j),jbbdivxmns(i+1,j)
                enddo
            enddo
            close(out_unit)
        endif

        deallocate(jbbkapxmns,jbbdivxmns,jbbdbobmns)
    end subroutine set_peq


    !=======================================================================
    subroutine read_fnml(file)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read matrix of pre-integrated functions, Eq. (13)
    !   [Park, Phys. Rev. Let 2009]
    !
    !*ARGUMENTS:
    !    file : character(512) (in)
    !       File path.
    !
    !-----------------------------------------------------------------------
        implicit none
        ! declare arguments
        character(*) :: file
        ! declare variables
        integer :: nfk,nft,in_unit
  
        if(verbose) print *, "Reading F^-1/2_mnl from file:"
        if(verbose) print *, "  ",trim(file)
        in_unit = get_free_file_unit(-1)
        open(unit=in_unit,file=trim(file),status="old",position="rewind",form="unformatted")
        read(in_unit) nfk,nft
        call bicube_alloc(fnml,nfk,nft,1)
        read(in_unit)fnml%xs(0:)
        read(in_unit)fnml%ys(0:)
        read(in_unit)fnml%fs(0:,0:,1)
        close(in_unit)
        call bicube_fit(fnml,'extrap','extrap')
        return
      end subroutine read_fnml

      
    !=======================================================================
    subroutine read_ipec_peq(file,write_log)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read pmodb and divxprp fourier components (psi,m) from input file.
    !   Holdover from PENT - used to test development.
    !
    !*ARGUMENTS:
    !    file : character (in)
    !       File path.
    !
    !-----------------------------------------------------------------------
        ! declare arguments
        logical :: write_log
        character(*), intent(in) :: file
        ! declare variables            
        complex(r8), dimension(mstep,mpert) :: lagbpar,divxprp
        integer :: i,ms,mp,iout,istep,in_unit,out_unit
        integer, dimension(:), allocatable :: mtemp
        type(cspline_type) :: outspl
    
          
        if(verbose) print *, "Reading pmodb and divxprp input file: "
        if(verbose) print *, '  ',trim(file)
        in_unit = get_free_file_unit(-1)
        open(unit=in_unit,file=file,status="old",position="rewind",form="unformatted")
        read(in_unit) ms,mp
        if((ms .ne. mstep) .or. (mp .ne. mpert)) then
            stop "ERROR: inputs - IPEC perturbations don't match equilibrium"
        endif
        read(in_unit)lagbpar
        read(in_unit)divxprp
        close(in_unit)
        
        ! form splines
        call cspline_alloc(dbob_m,mstep-1,mpert)
        call cspline_alloc(divx_m,mstep-1,mpert)
        dbob_m%xs(0:) = psifac(1:)
        divx_m%xs(0:) = psifac(1:)
        dbob_m%fs(0:,:) = lagbpar(:,:)
        divx_m%fs(0:,:) = divxprp(:,:)
        call cspline_fit(dbob_m,"extrap")
        call cspline_fit(divx_m,"extrap")
        
        ! write log - designed as check of reading routines
        if(write_log)then
            out_unit = get_free_file_unit(-1)
            print *, "Writing ipec lagb spline to pentrc_lagb.out"
            print *,"mpert = ",mpert
            print *,"mfac range = ",mfac(1),mfac(mpert)
            print *,"psifac range = ",psifac(1),psifac(mstep)
            iout = min(6,mpert)
            istep = mpert/iout
            call cspline_alloc(outspl,mstep-1,iout) ! reduced number for output
            outspl%xs = dbob_m%xs
            outspl%title(0) = "psi_n"
            do i=1,iout
                outspl%fs(:,i) = dbob_m%fs(:,i*istep)
                if(mfac(i*istep)<0)then
                    write(outspl%title(i),'(a1,i3.2)') 'm',mfac(i*istep)
                else
                    write(outspl%title(i),'(a1,i3.3)') 'm',mfac(i*istep)
                endif
                print *,"writing m ",mfac(i*istep)," to ",i,"th qty"
            enddo
            call cspline_fit(outspl,'extrap')
            open(unit=out_unit,file="pentrc_lagb.out",&
                status="unknown")
            call cspline_write(outspl,.true.,.false.,out_unit,0,.true.)
            close(out_unit)
        endif
    end subroutine read_ipec_peq

    
end module inputs
