! Native query oracle; supply query_contract.in, vac.in and vacin5 externally.
! The reference calls CHI directly at doubles, without pickup or loops.
program test_vacuum_query
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    use vacuum_mod, only: defglo, ent33, make_bltobp, diaplt, chi, pickup, mscfld
    use vglobal_mod, only: global_alloc, global_dealloc, kernelsign, ntsin0, &
        nths0, nfm, mtot, ndimlp, xobp, zobp, xloop, zloop, ldcon, lgpec, ieig, &
        iotty, outpest, inmode, outmod, farwal, mth, lmax, lmin, vacmat, vacmtiu, &
        xpla, zpla, plrad, deloop, epslp, grri, nths, nloop, nloopr, xplap, &
        zplap, bnlr, bnli, n, nxlpin, nzlpin
    implicit none
    integer :: surface_nodes, modes, iu, peak, side, l, k
    real(dp) :: qr(0:1,0:1), qz(0:1,0:1), rr(0:1,0:1), zz(0:1,0:1)
    real(dp) :: lr(0:1,0:1), lz(0:1,0:1), delta, band, center_r, center_z
    integer :: flags(0:1,0:1), legacy_flags(0:1,0:1)
    complex(dp) :: br(0:1,0:1), bz(0:1,0:1), bp(0:1,0:1)
    complex(dp) :: legacy_br(0:1,0:1), legacy_bz(0:1,0:1), legacy_bp(0:1,0:1)
    complex(dp) :: reference(3,2)
    complex(dp), allocatable :: matrix_ref(:,:), matrix_out(:,:)
    real(dp), allocatable :: re(:,:), im(:,:), region(:), cc(:,:), ss(:,:)
    real(dp) :: field_error, field_scale

    open(newunit=iu,file='query_contract.in',status='old')
    read(iu,*) surface_nodes,modes
    close(iu)
    if(surface_nodes<16 .or. surface_nodes>512) error stop 'Bounded surface size required'
    if(modes<1 .or. modes>40) error stop 'Bounded mode count required'
    kernelsign=1;ntsin0=surface_nodes+1;nths0=surface_nodes;nfm=modes;mtot=modes
    call global_alloc(nths0,nfm,mtot,ntsin0)
    ndimlp=4
    allocate(xobp(4),zobp(4),xloop(4),zloop(4))
    call defglo(surface_nodes)
    farwal=.true.;ldcon=0;lgpec=1;ieig=0
    open(iotty,file='mscvac.out',status='replace')
    open(outpest,file='pestotv',status='replace',form='formatted')
    open(inmode,file='vac.in',status='old',form='formatted')
    open(outmod,file='modovmc',status='replace',form='formatted')
    call ent33
    if(mth/=surface_nodes .or. lmax(1)-lmin(1)+1/=modes) error stop 'Source basis mismatch'
    call make_bltobp
    call diaplt
    allocate(matrix_ref(modes,modes),matrix_out(modes,modes))
    matrix_ref=cmplx(vacmat,vacmtiu,dp)
    open(newunit=iu,file='query_native_matrix.dat',status='replace')
    do k=1,modes
        do l=1,modes
            write(iu,'(2i6,2es26.17)') l,k,real(matrix_ref(l,k)),aimag(matrix_ref(l,k))
        end do
    end do
    close(iu)
    plrad=.5_dp*(maxval(xpla(:mth))-minval(xpla(:mth)))
    delta=plrad*deloop;band=plrad*epslp
    if(delta<=0 .or. band<=8*delta) error stop 'Oracle requires a resolvable legacy snap band'
    peak=maxloc(xpla(:mth),dim=1)
    center_r=.5_dp*(maxval(xpla(:mth))+minval(xpla(:mth)))
    center_z=.5_dp*(maxval(zpla(:mth))+minval(zpla(:mth)))
    qr(0,0)=xpla(peak)-.5_dp*band;qz(0,0)=zpla(peak)
    qr(0,1)=center_r+.13_dp*plrad;qz(0,1)=center_z+.07_dp*plrad
    qr(1,0)=xpla(peak);qz(1,0)=zpla(peak)
    qr(1,1)=xpla(peak)-.5_dp*delta;qz(1,1)=zpla(peak)

    ! The reference independently evaluates only the two resolved interior queries.
    allocate(re(5,4),im(5,4),region(4),cc(nths,modes),ss(nths,modes))
    cc=0;ss=0
    do l=1,modes
        cc(:mth,l)=grri(:mth,l)
        ss(:mth,l)=grri(:mth,modes+l)
    end do
    nloop=2;nloopr=0
    xloop(1:2)=qr(0,:);zloop(1:2)=qz(0,:)
    region=-1;re=0;im=0
    do side=1,5
        xobp(1:2)=xloop(1:2);zobp(1:2)=zloop(1:2)
        select case(side)
        case(1)
            zobp(1:2)=zobp(1:2)+delta
        case(2)
            zobp(1:2)=zobp(1:2)-delta
        case(3)
            xobp(1:2)=xobp(1:2)+delta
        case(4)
            xobp(1:2)=xobp(1:2)-delta
        end select
        call chi(xpla,zpla,xplap,zplap,-1,cc,ss,mth,1,re,im,side,bnlr,bnli,region)
    end do
    do k=1,2
        reference(1,k)=cmplx(re(3,k)-re(4,k),im(3,k)-im(4,k),dp)/(2*delta)
        reference(2,k)=cmplx(re(1,k)-re(2,k),im(1,k)-im(2,k),dp)/(2*delta)
        reference(3,k)=cmplx(n*im(5,k),-n*re(5,k),dp)/xloop(k)
    end do
    nxlpin=2;nzlpin=2
    call pickup(bnlr,bnli,1,1,flags,rr,zz,br,bz,bp, &
        preserve_query_coordinates=.true.,query_r=qr,query_z=qz)
    call check_exact_query()
    close(iotty);close(outpest);close(inmode);close(outmod)
    call global_dealloc
    call cleanup

    ! Exercise actual MSCFLD forwarding; its legacy GPEC path does not write WV.
    farwal=.true.
    call mscfld(matrix_out,modes,surface_nodes,surface_nodes,.true., &
        1,1,flags,rr,zz,br,bz,bp,preserve_query_coordinates=.true.,query_r=qr,query_z=qz)
    call check_exact_query()

    farwal=.true.
    call mscfld(matrix_out,modes,surface_nodes,surface_nodes,.true., &
        1,1,legacy_flags,lr,lz,legacy_br,legacy_bz,legacy_bp)
    farwal=.true.
    call mscfld(matrix_out,modes,surface_nodes,surface_nodes,.true., &
        1,1,flags,rr,zz,br,bz,bp,preserve_query_coordinates=.false.)
    if(any(flags/=legacy_flags) .or. any(rr/=lr) .or. any(zz/=lz)) &
        error stop 'Legacy optional contract changed coordinates or flags'
    if(any(br/=legacy_br) .or. any(bz/=legacy_bz) .or. any(bp/=legacy_bp)) &
        error stop 'Legacy optional contract changed field values'
    print *, 'PASS exact native query coordinates, direct CHI fields, masks and legacy contract'
    print *, 'DIRECT_CHI_RELATIVE_FIELD_ERROR',field_error
contains
    subroutine check_exact_query()
        integer :: j
        if(any(rr/=qr) .or. any(zz/=qz)) error stop 'Query coordinates relocated or downcast'
        if(any(flags(0,:)/=1)) error stop 'Resolved interior query misclassified'
        if(any(flags(1,:)/=-2)) error stop 'Source or crossing-stencil query not masked'
        if(.not.all(ieee_is_nan(real(br(1,:))))) error stop 'Invalid BR not unavailable'
        if(.not.all(ieee_is_nan(aimag(br(1,:))))) error stop 'Invalid imaginary BR not unavailable'
        if(.not.all(ieee_is_nan(real(bz(1,:))))) error stop 'Invalid BZ not unavailable'
        if(.not.all(ieee_is_nan(aimag(bz(1,:))))) error stop 'Invalid imaginary BZ not unavailable'
        if(.not.all(ieee_is_nan(real(bp(1,:))))) error stop 'Invalid Bphi not unavailable'
        if(.not.all(ieee_is_nan(aimag(bp(1,:))))) error stop 'Invalid imaginary Bphi not unavailable'
        field_error=0;field_scale=max(1e-20_dp,maxval(abs(reference)))
        do j=0,1
            field_error=max(field_error, &
                maxval(abs([br(0,j),bz(0,j),bp(0,j)]-reference(:,j+1)))/field_scale)
        end do
        if(.not.ieee_is_finite(field_error)) error stop 'Nonfinite resolved field'
        if(field_error>2e-11_dp) error stop 'Pickup field differs from independent direct CHI query'
    end subroutine check_exact_query
end program test_vacuum_query
