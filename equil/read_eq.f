c-----------------------------------------------------------------------
c     file read_eq.f.
c     reads equilibrium data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. read_eq_mod.
c     1. read_eq_get_flux.
c     2. read_eq_fluxgrid.
c     3. read_eq_miller.
c     4. read_eq_miller4.
c     5. read_eq_chease.
c     6. read_eq_chease2.
c     7 read_eq_chum.
c     8. read_eq_galkin.
c     9. read_eq_efit.
c     10. read_eq_rsteq.
c     11. read_eq_ldp_d.
c     12. read_eq_ldp_i.
c     13. read_eq_jsolver.
c     14. read_eq_lez.
c     15. read_eq_transp.
c     16. read_eq_sontag.
c     17. read_eq_tokamac.
c     18. read_eq_popov1.
c     19. read_eq_popov2.
c     20. read_eq_rtaylor.
c     21. read_eq_dump.
c     22. read_eq_wdn.
c     23. read_eq_lez_2
c-----------------------------------------------------------------------
c     subprogram 0. read_eq_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE read_eq_mod
      USE inverse_mod
      USE direct_mod
      USE utils_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. read_eq_get_flux.
c     computes poloidal flux from poloidal magnetic field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_get_flux(bfield,xs)

      REAL(r8), DIMENSION(0:,0:,:), INTENT(IN) :: bfield
      REAL(r8), DIMENSION(0:SIZE(bfield,1)-1), INTENT(OUT) :: xs

      INTEGER :: ix,mx,my
      TYPE(spline_type) :: tempx,tempy
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      mx=SIZE(bfield,1)-1
      my=SIZE(bfield,2)-1
      CALL spline_alloc(tempx,mx+1,1)
      CALL spline_alloc(tempy,my,1)
      tempx%xs=(/zero,rz_in%xs+one/)
      tempy%xs=rz_in%ys
      tempx%fs=0
c-----------------------------------------------------------------------
c     integrate over poloidal coordinate.
c-----------------------------------------------------------------------
      DO ix=0,mx
         tempy%fs(:,1)=rz_in%fs(ix,:,1)
     $        *(bfield(ix,:,2)*rz_in%fsx(ix,:,1)
     $        -bfield(ix,:,1)*rz_in%fsx(ix,:,2))
         CALL spline_fit(tempy,"periodic")
         CALL spline_int(tempy)
         tempx%fs(ix+1,1)=tempy%fsi(my,1)
      ENDDO
c-----------------------------------------------------------------------
c     integrate over radial coordinate.
c-----------------------------------------------------------------------
      CALL spline_fit(tempx,"extrap")
      CALL spline_int(tempx)
      xs=tempx%fsi(1:mx+1,1)
c-----------------------------------------------------------------------
c     deallocate spline types.
c-----------------------------------------------------------------------
      CALL spline_dealloc(tempx)
      CALL spline_dealloc(tempy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_get_flux
c-----------------------------------------------------------------------
c     subprogram 2. read_eq_fluxgrid.
c     reads equilibrium data from fluxgrid.dat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_fluxgrid

      LOGICAL :: file_stat
      INTEGER :: ix,iy,mxpie,mx,my
      REAL(r8) :: bo
      REAL(r8), DIMENSION(:,:,:), POINTER :: bfield
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
c-----------------------------------------------------------------------
c     read sizes.
c-----------------------------------------------------------------------
      READ(in_unit,*)mx
      READ(in_unit,*)my
      READ(in_unit,*)mxpie
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      ALLOCATE(bfield(0:mx,0:my,3))
      CALL spline_alloc(sq_in,mx,4)
      CALL bicube_alloc(rz_in,mx,my,2)
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit,*)sq_in%fs(0,1:3)
      READ(in_unit,*)ro,zo,bo
      DO ix=0,mx
         READ(in_unit,*)sq_in%fs(ix,1:3)
         READ(in_unit,*)(rz_in%fs(ix,iy,:),bfield(ix,iy,:),iy=0,my)
      ENDDO
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     fit input to bicubic splines.
c-----------------------------------------------------------------------
      rz_in%xs=(/(ix,ix=0,mx)/)
      rz_in%ys=(/(iy,iy=0,my)/)/REAL(my,r8)
      CALL bicube_fit(rz_in,"extrap","periodic")
c-----------------------------------------------------------------------
c     compute poloidal flux as new radial variable.
c-----------------------------------------------------------------------
      CALL read_eq_get_flux(bfield,sq_in%xs)
      DEALLOCATE(bfield)
      psio=sq_in%xs(mx)
      sq_in%xs=sq_in%xs/psio
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_fluxgrid
c-----------------------------------------------------------------------
c     subprogram 3. read_eq_miller.
c     reads equilibrium data from miller-type files.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_miller

      LOGICAL :: file_stat
      INTEGER :: ix,mx,my
      REAL(r8), DIMENSION(:), POINTER :: xs
      REAL(r8), DIMENSION(:,:), POINTER :: sq
      REAL(r8), DIMENSION(:,:,:), POINTER :: rz
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
c-----------------------------------------------------------------------
c     read sizes and allocate local arrays.
c-----------------------------------------------------------------------
      READ(in_unit)my,mx
      my=my-1
      mx=mx-1
      ALLOCATE(xs(0:mx))
      ALLOCATE(sq(0:mx,3))
      ALLOCATE(rz(0:mx,0:my,2))
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit)xs
      READ(in_unit)sq(:,1)
      READ(in_unit)sq(:,2)
      READ(in_unit)sq(:,3)
      READ(in_unit)(rz(ix,:,1),ix=0,mx)
      READ(in_unit)(rz(ix,:,2),ix=0,mx)
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     copy to eq_in.
c-----------------------------------------------------------------------
      ro=rz(0,0,1)
      zo=rz(0,0,2)
      psio=xs(mx)-xs(0)
      CALL spline_alloc(sq_in,mx,4)
      CALL bicube_alloc(rz_in,mx,my,2)
      sq_in%xs=(xs-xs(0))/psio
      sq_in%fs(:,1:3)=sq
      sq_in%fs(:,2)=sq_in%fs(:,2)*mu0
      rz_in%fs=rz
      DEALLOCATE(xs,sq,rz)
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_miller
c-----------------------------------------------------------------------
c     subprogram 4. read_eq_miller4.
c     reads equilibrium data from miller-type files.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_miller4

      LOGICAL :: file_stat
      INTEGER :: ix,mx,my
      REAL(r4), DIMENSION(:), POINTER :: xs
      REAL(r4), DIMENSION(:,:), POINTER :: sq
      REAL(r4), DIMENSION(:,:,:), POINTER :: rz
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
c-----------------------------------------------------------------------
c     read sizes and allocate local arrays.
c-----------------------------------------------------------------------
      READ(in_unit)my,mx
      my=my-1
      mx=mx-1
      ALLOCATE(xs(0:mx))
      ALLOCATE(sq(0:mx,3))
      ALLOCATE(rz(0:mx,0:my,2))
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit)xs
      READ(in_unit)sq(:,1)
      READ(in_unit)sq(:,2)
      READ(in_unit)sq(:,3)
      READ(in_unit)(rz(ix,:,1),ix=0,mx)
      READ(in_unit)(rz(ix,:,2),ix=0,mx)
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     copy to eq_in.
c-----------------------------------------------------------------------
      ro=rz(0,0,1)
      zo=rz(0,0,2)
      psio=xs(mx)-xs(0)
      CALL spline_alloc(sq_in,mx,4)
      CALL bicube_alloc(rz_in,mx,my,2)
      sq_in%xs=(xs-xs(0))/psio
      sq_in%fs(:,1:3)=sq
      sq_in%fs(:,2)=sq_in%fs(:,2)*mu0
      rz_in%fs=rz
      DEALLOCATE(xs,sq,rz)
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_miller4
c-----------------------------------------------------------------------
c     subprogram 5. read_eq_chease.
c     reads data from chease.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_chease

      INTEGER :: ntnova,npsi1,nsym,ma,mtau
      REAL(r8), DIMENSION(5) :: axx
      REAL(r8), DIMENSION(:), POINTER :: zcpr,zcppr,zq,zdq,ztmf,
     $     ztp,zfb,zfbp,zpsi,zpsim
      REAL(r8), DIMENSION(:,:), POINTER :: buffer
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)ntnova,npsi1,nsym
      READ(in_unit)axx
c-----------------------------------------------------------------------
c     allocate and read 1D arrays.
c-----------------------------------------------------------------------
      ALLOCATE(zcpr(npsi1-1),zcppr(npsi1),zq(npsi1),zdq(npsi1),
     $     ztmf(npsi1),ztp(npsi1),zfb(npsi1),zfbp(npsi1),zpsi(npsi1),
     $     zpsim(npsi1-1))
      READ(in_unit)zcpr
      READ(in_unit)zcppr
      READ(in_unit)zq
      READ(in_unit)zdq
      READ(in_unit)ztmf
      READ(in_unit)ztp
      READ(in_unit)zfb
      READ(in_unit)zfbp
      READ(in_unit)zpsi
      READ(in_unit)zpsim
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ma=npsi1-1
      CALL spline_alloc(sq_in,ma,4)
      psio=zpsi(npsi1)-zpsi(1)
      sq_in%xs=(zpsi(1:npsi1)-zpsi(1))/psio
      sq_in%fs(:,1)=ztmf(1:npsi1)
      sq_in%fs(:,2)=zcppr(1:npsi1)
      sq_in%fs(:,3)=zq(1:npsi1)
      CALL spline_fit(sq_in,"extrap")
      DEALLOCATE(zcpr,zcppr,zq,zdq,ztmf,ztp,zfb,zfbp,zpsi,zpsim)
c-----------------------------------------------------------------------
c     integrate pressure.
c-----------------------------------------------------------------------
      CALL spline_int(sq_in)
      sq_in%fs(:,2)=(sq_in%fsi(:,2)-sq_in%fsi(ma,2))*psio
      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     allocate, read, and copy 2D arrays.
c-----------------------------------------------------------------------
      mtau=ntnova
      CALL bicube_alloc(rz_in,ma,mtau,2)
      ALLOCATE(buffer(ntnova+3,npsi1))
      READ(in_unit)buffer
      ro=buffer(1,1)
      rz_in%fs(:,:,1)=TRANSPOSE(buffer(1:ntnova+1,1:))
      READ(in_unit)buffer
      zo=buffer(1,1)
      rz_in%fs(:,:,2)=TRANSPOSE(buffer(1:ntnova+1,1:))
      CALL bin_close(in_unit)
      DEALLOCATE(buffer)
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_chease
c-----------------------------------------------------------------------
c     subprogram 6. read_eq_chease2.
c     reads data from chease.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_chease2

      INTEGER :: ntnova,npsi1,nsym,ma,mtau
      REAL(r8), DIMENSION(5) :: axx
      REAL(r8), DIMENSION(:), POINTER :: zcpr,zcppr,zq,zdq,ztmf,
     $     ztp,zfb,zfbp,zpsi,zpsim
      REAL(r8), DIMENSION(:,:), POINTER :: zrcp,zzcp,zjacm,zjac
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(3i5)
 20   FORMAT(5e22.15)
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,10)ntnova,npsi1,nsym
      READ(in_unit,20)axx
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(zcpr(npsi1-1),zcppr(npsi1),zq(npsi1),zdq(npsi1),
     $     ztmf(npsi1),ztp(npsi1),zfb(npsi1),zfbp(npsi1),zpsi(npsi1),
     $     zpsim(npsi1-1))
      ALLOCATE(zrcp(ntnova+3,npsi1),zzcp(ntnova+3,npsi1),
     $     zjacm(ntnova+3,npsi1),zjac(ntnova+3,npsi1))
c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      READ(in_unit,20)zcpr
      READ(in_unit,20)zcppr
      READ(in_unit,20)zq
      READ(in_unit,20)zdq
      READ(in_unit,20)ztmf
      READ(in_unit,20)ztp
      READ(in_unit,20)zfb
      READ(in_unit,20)zfbp
      READ(in_unit,20)zpsi
      READ(in_unit,20)zpsim
      READ(in_unit,20)zrcp
      READ(in_unit,20)zzcp
      READ(in_unit,20)zjacm
      READ(in_unit,20)zjac
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ma=npsi1-1
      CALL spline_alloc(sq_in,ma,4)
      psio=zpsi(npsi1)-zpsi(1)
      sq_in%xs=(zpsi(1:npsi1)-zpsi(1))/psio
      sq_in%fs(:,1)=ztmf(1:npsi1)
      sq_in%fs(:,2)=zcppr(1:npsi1)
      sq_in%fs(:,3)=zq(1:npsi1)
      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     integrate pressure.
c-----------------------------------------------------------------------
      CALL spline_int(sq_in)
      sq_in%fs(:,2)=(sq_in%fsi(:,2)-sq_in%fsi(ma,2))*psio
      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      mtau=ntnova
      ro=zrcp(1,1)
      zo=zzcp(1,1)
      CALL bicube_alloc(rz_in,ma,mtau,2)
      rz_in%fs(:,:,1)=TRANSPOSE(zrcp(1:ntnova+1,1:))
      rz_in%fs(:,:,2)=TRANSPOSE(zzcp(1:ntnova+1,1:))
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(zcpr,zcppr,zq,zdq,ztmf,ztp,zfb,zfbp,zpsi,zpsim)
      DEALLOCATE(zrcp,zzcp,zjacm,zjac)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_chease2
c-----------------------------------------------------------------------
c     subprogram 7. read_eq_chum.
c     reads data from Ming Chu's equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_chum

      INTEGER :: jpsi,itht,ntor,ma,mtau
      REAL(4) :: rquot,omega0,growth
      REAL(4), DIMENSION(:), POINTER :: psi,f,p,q
      REAL(4), DIMENSION(:,:), POINTER :: rcc,zcc

      INTEGER, PARAMETER :: mpts=4
      INTEGER :: ipt
      REAL(r8), DIMENSION(mpts) :: xs
      REAL(r8), DIMENSION(mpts,2) :: fs
      REAL(r8), DIMENSION(0:2,2) :: fs0
c-----------------------------------------------------------------------
c     open file and read scalars.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)jpsi,itht,ntor
      READ(in_unit)rquot,omega0,growth
c-----------------------------------------------------------------------
c     allocate and read arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psi(jpsi),f(jpsi),p(jpsi),q(jpsi),
     $     rcc(itht,jpsi),zcc(itht,jpsi))
      READ(in_unit)psi
      READ(in_unit)f
      READ(in_unit)p
      READ(in_unit)q
      READ(in_unit)rcc
      READ(in_unit)zcc
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     copy and modify 1D arrays.
c-----------------------------------------------------------------------
      psio=psi(jpsi)
      ma=jpsi-1
      CALL spline_alloc(sq_in,ma,4)
      sq_in%xs(:)=psi/psio
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=mu0*p
      sq_in%fs(:,3)=q
c-----------------------------------------------------------------------
c     extrapolate coordinates to magnetic axis.
c-----------------------------------------------------------------------
      mtau=itht
      xs=psi(1:mpts)
      DO ipt=1,mpts
         fs(ipt,1)=SUM(rcc(1:mtau,ipt))/mtau
         fs(ipt,2)=SUM(zcc(1:mtau,ipt))/mtau
      ENDDO
      fs0=interpolate(xs,fs,zero)
      ro=fs0(0,1)
      zo=fs0(0,2)
c-----------------------------------------------------------------------
c     copy 2D equilibrium arrays.
c-----------------------------------------------------------------------
      CALL bicube_alloc(rz_in,ma,mtau,2)
      rz_in%fs(:,0:mtau-1,1)=TRANSPOSE(rcc)
      rz_in%fs(:,0:mtau-1,2)=TRANSPOSE(zcc)
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,f,p,q,rcc,zcc)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_chum
c-----------------------------------------------------------------------
c     subprogram 8. read_eq_galkin.
c     reads data from Sergei Galkin's equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_galkin

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: nfeq,n,m,i,j,ia,itau,ma,mtau
      REAL(r8) :: dadv,anun,anu0,um,vm
      REAL(r8), DIMENSION(:), POINTER :: uk,vk,psi,p,f,q,anu,pp,
     $     pff,a
      REAL(r8), DIMENSION(:,:), POINTER :: ro_in
c-----------------------------------------------------------------------
c     read data.  the reals and integers are likely a mix of 64 and 32
c     bit data, so using bin_open will not help porting.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)nfeq,n,m,dadv,anun,anu0,um,vm
      ALLOCATE(uk(2:m-1),vk(2:m-1),p(n),f(n),q(n),anu(n),pp(n),pff(n),
     $     a(n),ro_in(n,2:m-1))
      READ(in_unit)uk,vk,ro_in,p,f,q,anu,pp,pff,a
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     integrate psi profile.
c-----------------------------------------------------------------------
      ALLOCATE(psi(n))
      psi(n)=1
      DO i=n,2,-1
         psi(i-1)=psi(i)-anu(i)*(a(i)-a(i-1))
      ENDDO
c-----------------------------------------------------------------------
c     copy and modify 1D arrays.
c-----------------------------------------------------------------------
      ma=n-1
      CALL spline_alloc(sq_in,ma,4)
      psio=anun
      sq_in%xs=psi(1:n)
      sq_in%fs(:,1)=f(1:n)
      sq_in%fs(:,2)=p(1:n)
      sq_in%fs(:,3)=q(1:n)
c-----------------------------------------------------------------------
c     copy and modify 2D arrays.
c-----------------------------------------------------------------------
      mtau=m-2
      CALL bicube_alloc(rz_in,ma,mtau,2)
      DO ia=0,ma
         i=ia+1
         DO itau=0,mtau-1
            j=itau+2
            rz_in%fs(ia,itau,1)=um+ro_in(i,j)*(uk(j)-um)
            rz_in%fs(ia,itau,2)=vm+ro_in(i,j)*(vk(j)-vm)
         ENDDO
      ENDDO
      ro=um
      zo=vm
      psio=anun
c-----------------------------------------------------------------------
c     process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(uk,vk,p,f,q,anu,pp,pff,a,ro_in,psi)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_galkin
c-----------------------------------------------------------------------
c     subprogram 9. read_eq_efit.
c     reads data from General Atomic's EFIT equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_efit

      INTEGER :: i,j,nw,nh,ia,mr,mz,ma
      INTEGER :: ios
      REAL(r8) :: bcentr,cpasma,rgrid,rmaxis,rzero,ssibry1,
     $     ssibry2,ssimag1,ssimag2,xdim,xdum,zdim,zmaxis,zmid
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),'OLD')
      READ(in_unit,'(52x,2i4)')nw,nh
      READ(in_unit,'(5e16.9)')xdim,zdim,rzero,rgrid,zmid
      READ(in_unit,'(5e16.9)')rmaxis,zmaxis,ssimag1,ssibry1,bcentr
      READ(in_unit,'(5e16.9)')cpasma,ssimag2,xdum,rmaxis,xdum
      READ(in_unit,'(5e16.9)')zmaxis,xdum,ssibry2,xdum,xdum
      CALL spline_alloc(sq_in,nw-1,4)
      READ(in_unit,'(5e16.9)')(sq_in%fs(i,1),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq_in%fs(i,2),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq_in%fs(i,3),i=0,nw-1)
      READ(in_unit,'(5e16.9)')(sq_in%fs(i,3),i=0,nw-1)
      CALL bicube_alloc(psi_in,nw-1,nh-1,1)
      READ(in_unit,'(5e16.9)')((psi_in%fs(i,j,1),i=0,nw-1),j=0,nh-1)
      READ(in_unit,'(5e16.9)',iostat=ios)(sq_in%fs(i,3),i=0,nw-1)
      IF(ios /= 0)sq_in%fs(i,3)=0
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     translate to internal quantities.
c-----------------------------------------------------------------------
      mr=nw-1
      mz=nh-1
      ma=nw-1
      psi_in%mx=mr
      psi_in%my=mz
      rmin=rgrid
      rmax=rgrid+xdim
      zmin=-zdim/2
      zmax=zdim/2
      psio=ssibry1-ssimag1
      sq_in%xs=(/(ia,ia=0,ma)/)/dfloat(ma)
      sq_in%fs(:,1)=ABS(sq_in%fs(:,1))
      sq_in%fs(:,2)=MAX(sq_in%fs(:,2)*mu0,zero)
c-----------------------------------------------------------------------
c     copy and convert 2D quantities.
c-----------------------------------------------------------------------
      psi_in%fs=ssibry1-psi_in%fs
c-----------------------------------------------------------------------
c     normalize sign convention.
c-----------------------------------------------------------------------
      IF(psio < 0)THEN
         psio=-psio
         psi_in%fs=-psi_in%fs
      ENDIF
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_efit
c-----------------------------------------------------------------------
c     subprogram 10. read_eq_rsteq.
c     reads data from rsteq equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_rsteq
      
      INTEGER :: i,j,nw,nh,ia,mr,mz,ma
      INTEGER ::  nrst
c-----------------------------------------------------------------------
c     open equilibrium file and read sizes.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)nw,nh,nrst
      READ(in_unit)rmin,rmax,zmin,zmax
      READ(in_unit)ro,zo
      mr=nw-1
      mz=2*nh-3
      ma=nrst-1
      psi_in%mx=mr
      psi_in%my=mz
c-----------------------------------------------------------------------
c     allocate spline types.
c-----------------------------------------------------------------------
      CALL bicube_alloc(psi_in,mr,mz,1)
      CALL spline_alloc(sq_in,nrst-1,3)
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit)((psi_in%fs(i,j,1),i=0,nw-1),j=nh-2,mz)
      DO i=0,nw-1
        DO j=0,nh-3
          psi_in%fs(i,j,1)=psi_in%fs(i,mz-j,1)
        ENDDO
      ENDDO
      READ(in_unit)(sq_in%fs(i,1),i=0,nrst-1)   ! psi
      psio=sq_in%fs(0,1)-sq_in%fs(nrst-1,1)
      READ(in_unit)(sq_in%fs(i,1),i=0,nrst-1)   ! F
      READ(in_unit)(sq_in%fs(i,2),i=0,nrst-1)   ! pressure
      READ(in_unit)(sq_in%fs(i,3),i=0,nrst-1)   ! q
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     translate to internal quantities.
c-----------------------------------------------------------------------
      sq_in%xs=(/(ia,ia=0,ma)/)/dfloat(ma)
      sq_in%fs(:,1)=ABS(sq_in%fs(:,1))
      sq_in%fs(:,2)=MAX(sq_in%fs(:,2)*mu0,zero)
c-----------------------------------------------------------------------
c     normalize sign convention.
c-----------------------------------------------------------------------
      IF(psio < 0)THEN
         psio=-psio
         psi_in%fs=-psi_in%fs
      ENDIF
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_rsteq
c-----------------------------------------------------------------------
c     subprogram 11. read_eq_ldp_d.
c     reads data from L. Don Pearlstein's TEQ equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_ldp_d

      INTEGER :: mr,mz,ma
c-----------------------------------------------------------------------
c     open input, file, read scalar data, and allocate arrays.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)mr,mz,ma
      READ(in_unit)rmin,rmax,zmin,zmax
      ma=ma-1
      mr=mr-1
      mz=mz-1
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(psi_in,mr,mz,1)
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit)psi_in%fs
      READ(in_unit)sq_in%xs
      READ(in_unit)sq_in%fs(:,1)
      READ(in_unit)
      READ(in_unit)sq_in%fs(:,2)
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     revise 1D arrays.
c-----------------------------------------------------------------------
      psio=sq_in%xs(0)
      sq_in%xs=1-sq_in%xs/psio
      sq_in%fs(:,2)=sq_in%fs(:,2)*mu0
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_ldp_d
c-----------------------------------------------------------------------
c     subprogram 12. read_eq_ldp_i.
c     reads data from L. Don Pearlstein's inverse equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_ldp_i

      INTEGER :: mx,my
      REAL(r8), DIMENSION(:), POINTER :: psi,f,p,q
      REAL(r8), DIMENSION(:,:), POINTER :: r,z
c-----------------------------------------------------------------------
c     open data file, read scalar data.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)mx,my
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      my=my-1
      mx=mx-1
      ALLOCATE(psi(0:mx),f(0:mx),p(0:mx),q(0:mx),
     $     r(0:my,0:mx),z(0:my,0:mx))
c-----------------------------------------------------------------------
c     read binary data.
c-----------------------------------------------------------------------
      READ(in_unit)psi
      READ(in_unit)f
      READ(in_unit)p
      READ(in_unit)q
      READ(in_unit)r
      READ(in_unit)z
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     copy and revise 1D arrays.
c-----------------------------------------------------------------------
      psio=psi(mx)-psi(0)
      CALL spline_alloc(sq_in,mx,4)
      sq_in%xs=(psi-psi(0))/psio
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=p*mu0
      sq_in%fs(:,3)=q
      if(psio.lt.0)psio=-psio
c-----------------------------------------------------------------------
c     copy and revise 2D arrays.
c-----------------------------------------------------------------------
      ro=r(0,0)
      zo=z(0,0)
      CALL bicube_alloc(rz_in,mx,my,2)
      rz_in%fs(:,:,1)=TRANSPOSE(r)
      rz_in%fs(:,:,2)=TRANSPOSE(z)
c-----------------------------------------------------------------------
c     process inverse equilibrium and deallocate local arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,f,p,q,r,z)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_eq_ldp_i
c-----------------------------------------------------------------------
c     subprogram 13. read_eq_jsolver.
c     reads data from Steve Jardin's JSOLVER inverse equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_jsolver

      INTEGER :: nz1,nz2,nthe,npsi,isym
      REAL(r8) :: dr,dt,xzero,p0,plimpst,pminpst
      REAL(r8), DIMENSION(:), POINTER :: p,ppxx,q,qp,gxx,gpx,
     $     fb,fbp,f,fp,psival
      REAL(r8), DIMENSION(:,:), POINTER :: x,z,aj3,aj

      INTEGER :: ma,mtau
      REAL(r8), DIMENSION(:), POINTER :: debug
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(5e16.8)
c-----------------------------------------------------------------------
c     open input file, read scalars, and allocate arrays.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,*)nz1,nz2,nthe,npsi,isym,dr,dt
      READ(in_unit,*)xzero,p0,plimpst,pminpst
      ALLOCATE(p(nz2),ppxx(nz2),q(nz2),qp(nz2),gxx(nz2),gpx(nz2),
     $     fb(nz2),fbp(nz2),f(nz2),fp(nz2),psival(nz2),
     $     x(nz1,nz2),z(nz1,nz2),aj3(nz1,nz2),aj(nz1,nz2))
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit,10)p
      READ(in_unit,10)ppxx
      READ(in_unit,10)q
      READ(in_unit,10)qp
      READ(in_unit,10)gxx
      READ(in_unit,10)gpx
      READ(in_unit,10)fb
      READ(in_unit,10)fbp
      READ(in_unit,10)f
      READ(in_unit,10)fp
      READ(in_unit,10)psival
      READ(in_unit,10)x
      READ(in_unit,10)z
      READ(in_unit,10)aj3
c      READ(in_unit,10)aj
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ma=npsi-1
      CALL spline_alloc(sq_in,ma,4)
      psio=-psival(1)
      sq_in%xs=1+psival(1:npsi)/psio
      sq_in%fs(0:ma-1,1)=(gxx(1:npsi-1)+gxx(1:npsi))*xzero/2
      sq_in%fs(0:ma-1,2)=(p(1:npsi-1)+p(1:npsi))/2
      sq_in%fs(:,3)=q(1:npsi)
      sq_in%fs(ma,1)=xzero
      sq_in%fs(ma,2)=0
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      ro=x(3,1)
      zo=z(3,1)
      debug => z(1:nthe,2)
      IF(isym == 0)THEN
         mtau=nthe
         CALL bicube_alloc(rz_in,ma,mtau,2)
         rz_in%fs(:,:,1)=TRANSPOSE(x(3:nthe+2,1:npsi))
         rz_in%fs(:,:,2)=TRANSPOSE(z(3:nthe+2,1:npsi))
      ELSE
         mtau=nthe-1
         CALL bicube_alloc(rz_in,ma,2*mtau,2)
         rz_in%fs(:,0:mtau,1)=TRANSPOSE(x(3:nthe+2,1:npsi))
         rz_in%fs(:,0:mtau,2)=TRANSPOSE(z(3:nthe+2,1:npsi))
         rz_in%fs(:,mtau:2*mtau,1)=TRANSPOSE(x(nthe+1:3:-1,1:npsi))
         rz_in%fs(:,mtau:2*mtau,2)=-TRANSPOSE(z(nthe+1:3:-1,1:npsi))
      ENDIF
c-----------------------------------------------------------------------
c     deallocate local arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(p,ppxx,q,qp,gxx,gpx,fb,fbp,f,fp,psival,x,z,aj3,aj)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_jsolver
c-----------------------------------------------------------------------
c     subprogram 14. read_eq_lez.
c     reads data from Leonid E. Zakharov's equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_lez

      LOGICAL :: file_stat
      INTEGER :: ia,ma,mtau
      REAL(r4), DIMENSION(:), POINTER :: psis,qs,fs,ps
      REAL(r4), DIMENSION(:,:), POINTER :: rg,zg
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
c-----------------------------------------------------------------------
c     read sizes and allocate input arrays.
c-----------------------------------------------------------------------
      READ(in_unit)ma,mtau
      ma=ma-1
      mtau=mtau-1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psis(0:ma),qs(0:ma),fs(0:ma),ps(0:ma),
     $     rg(0:ma,0:mtau),zg(0:ma,0:mtau))
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
c-----------------------------------------------------------------------
c     read 2d data.
c-----------------------------------------------------------------------
      DO ia=0,ma
         READ(in_unit)rg(ia,:)
         READ(in_unit)zg(ia,:)
      ENDDO
c-----------------------------------------------------------------------
c     read 1d data.
c-----------------------------------------------------------------------
      READ(in_unit)psis
      READ(in_unit)qs
      READ(in_unit)fs
      READ(in_unit)
      READ(in_unit)ps
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     revise 1D data.
c-----------------------------------------------------------------------
      psio=-psis(ma)
      ma=ma-1
      sq_in%xs=-psis/psio
      sq_in%fs(:,1)=fs
      sq_in%fs(:,2)=ps*mu0*twopi*1e6
      sq_in%fs(:,3)=qs
      psio=psio/twopi
c-----------------------------------------------------------------------
c     integrate pressure profile.
c-----------------------------------------------------------------------
      CALL spline_fit(sq_in,"extrap")
      CALL spline_int(sq_in)
      sq_in%fs(:,2)=-(sq_in%fsi(:,2)-sq_in%fsi(ma,2))*psio
c-----------------------------------------------------------------------
c     revise 2D profiles.
c-----------------------------------------------------------------------
      ro=rg(0,0)
      zo=zg(0,0)
      rz_in%fs(:,:,1)=rg
      rz_in%fs(:,:,2)=zg
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psis,qs,fs,ps,rg,zg)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_lez
c-----------------------------------------------------------------------
c     subprogram 15. read_eq_transp.
c     reads data from transp.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_transp

      LOGICAL :: file_stat
      INTEGER :: ia,ma,mtau,nsurf,ntheta
      REAL(r4), DIMENSION(:), POINTER :: psis,qs,fs,ps,thetas,rhos
      REAL(r4), DIMENSION(:,:), POINTER :: rg,zg
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1p,5e14.6)
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
c-----------------------------------------------------------------------
c     read sizes and allocate input arrays.
c-----------------------------------------------------------------------
      READ(in_unit,'(//////)')
      READ(in_unit,*)ntheta
      READ(in_unit,'()')
      READ(in_unit,*)nsurf
      ma=nsurf-1
      mtau=ntheta-1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psis(0:ma),qs(0:ma),fs(0:ma),ps(0:ma),
     $     rg(0:ma,0:mtau),zg(0:ma,0:mtau),thetas(0:mtau),rhos(0:ma))
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
c-----------------------------------------------------------------------
c     read 1d data.
c-----------------------------------------------------------------------
      READ(in_unit,'()')
      READ(in_unit,10)thetas
      READ(in_unit,'()')
      READ(in_unit,10)rhos
      READ(in_unit,'()')
      READ(in_unit,10)psis
      READ(in_unit,'()')
      READ(in_unit,10)ps
      READ(in_unit,'()')
      READ(in_unit,10)qs
      READ(in_unit,'()')
      READ(in_unit,10)fs
c-----------------------------------------------------------------------
c     read 2d data.
c-----------------------------------------------------------------------
      READ(in_unit,'()')
      READ(in_unit,10)(rg(ia,:),ia=0,ma)
      READ(in_unit,'()')
      READ(in_unit,10)(zg(ia,:),ia=0,ma)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     revise 1D data.
c-----------------------------------------------------------------------
      psio=psis(ma)
      ma=ma
      sq_in%xs=psis/psio
      sq_in%fs(:,1)=fs
      sq_in%fs(:,2)=(ps-ps(ma))*mu0
      sq_in%fs(:,3)=qs
c-----------------------------------------------------------------------
c     revise 2D profiles.
c-----------------------------------------------------------------------
      ro=rg(0,0)
      zo=zg(0,0)
      rz_in%fs(:,:,1)=rg
      rz_in%fs(:,:,2)=zg
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psis,qs,fs,ps,rg,zg,thetas,rhos)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate program.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_transp
c-----------------------------------------------------------------------
c     subprogram 16. read_eq_sontag.
c     reads data from Aaron Sontag's direct solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_sontag

      INTEGER :: mr,mz,ma
      REAL(r8) :: psilim
      REAL(r8), DIMENSION(:), POINTER :: psi,f,p,q
      REAL(r8), DIMENSION(:,:), POINTER :: psig
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),'OLD')
      READ(in_unit,*)mr,mz,ma,rmin,rmax,zmin,zmax,psilim
      ALLOCATE(psi(ma),f(ma),p(ma),q(ma))
      READ(in_unit,*)psi,f,p,q
      ALLOCATE(psig(mz,mr))
      READ(in_unit,*)psig
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     store data.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,ma-2,4)
      psio=(psilim-psi(1))/twopi
      sq_in%xs=(psi(2:)-psi(1))/(psilim-psi(1))
      sq_in%fs(:,1)=f(2:)
      sq_in%fs(:,2)=p(2:)*mu0
      sq_in%fs(:,3)=q(2:)
      CALL bicube_alloc(psi_in,mr-1,mz-1,1)
      psi_in%fs(:,:,1)=(psilim-TRANSPOSE(psig))/twopi
c-----------------------------------------------------------------------
c     normalize sign convention.
c-----------------------------------------------------------------------
      IF(psio < 0)THEN
         psio=-psio
         psi_in%fs=-psi_in%fs
      ENDIF
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,f,p,q,psig)
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_sontag
c-----------------------------------------------------------------------
c     subprogram 17. read_eq_tokamac.
c     reads data from Mike Mauel's direct solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_tokamac

      INTEGER :: mr,mz,ma
      REAL(r8) :: beta_in,betap_in
      REAL(r8), DIMENSION(:), POINTER :: psi,f,p,q
      REAL(r8), DIMENSION(:,:), POINTER :: psig
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),'OLD')
      READ(in_unit,'(19(/))')
      READ(in_unit,*)mr,mz,ma
      ma=ma-1
      ALLOCATE(psi(0:ma),f(0:ma),p(0:ma),q(0:ma),psig(mr,mz))
      READ(in_unit,*)rmin,rmax,zmin,zmax
      READ(in_unit,*)beta_in,betap_in
      READ(in_unit,*)psig
      READ(in_unit,*)psi
      READ(in_unit,*)f
      READ(in_unit,*)p
      READ(in_unit,*)q
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     store data.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,ma,4)
      psio=psi(ma)-psi(0)
      sq_in%xs=(psi-psi(0))/psio
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=p*mu0
      sq_in%fs(:,3)=q
      psio=psio/twopi
      CALL bicube_alloc(psi_in,mr-1,mz-1,1)
      psi_in%fs(:,:,1)=(psi(ma)-psig(:,:))/twopi
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,f,p,q,psig)
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_tokamac
c-----------------------------------------------------------------------
c     subprogram 18. read_eq_popov1.
c     reads data from Sasha Popov's ascii file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_popov1

      INTEGER :: n,m1,ma,mtau
      REAL(r8) :: rum,psim,psip,dummy0
      REAL(r8), DIMENSION(:), POINTER :: psi,f,p,q,dummy
      REAL(r8), DIMENSION(:,:), POINTER :: r,z
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(4(e15.8,1x))
c-----------------------------------------------------------------------
c     open input file, read scalars, and allocate arrays.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,*)n,m1,rum,psim,psip,dummy0
c-----------------------------------------------------------------------
c     allocate and read arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psi(n),q(n),p(n),f(n),dummy(n),r(n,m1),z(n,m1))
      READ(in_unit,10)psi
      READ(in_unit,10)q
      READ(in_unit,10)p
      READ(in_unit,10)dummy
      READ(in_unit,10)f
      READ(in_unit,10)dummy
      READ(in_unit,10)dummy
      READ(in_unit,10)r(:,2:m1)
      READ(in_unit,'(/)')
      READ(in_unit,10)z(:,2:m1)
      call ascii_close(in_unit)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ma=n-1
      CALL spline_alloc(sq_in,ma,4)
      psio=ABS(psim-psip)
      sq_in%xs=psi
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=p
      sq_in%fs(:,3)=q
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      ro=r(1,m1)
      zo=z(1,m1)
      mtau=m1-1
      CALL bicube_alloc(rz_in,ma,mtau,2)
      rz_in%fs(:,0:mtau-1,1)=r(:,2:m1)
      rz_in%fs(:,0:mtau-1,2)=z(:,2:m1)
      rz_in%fs(:,mtau,:)=rz_in%fs(:,0,:)
c-----------------------------------------------------------------------
c     deallocate local arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,q,p,dummy,r,z)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_popov1
c-----------------------------------------------------------------------
c     subprogram 19. read_eq_popov2.
c     reads data from Sasha Popov's binary file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_popov2

      INTEGER :: nr1,nj,m1max,i,itheta,ma,m
      REAL(r8) :: r0,tmp,theta,dtheta
      REAL(r8), DIMENSION(:), POINTER :: qr,peqr,fr,psipr,rho
      REAL(r8), DIMENSION(:,:), POINTER :: rm,zm
c-----------------------------------------------------------------------
c     open input file, read scalars, and allocate arrays.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      read(1)nr1,nj,m1max,r0,i,tmp
c-----------------------------------------------------------------------
c     allocate and read arrays.
c-----------------------------------------------------------------------
      ALLOCATE(rho(nr1),qr(nr1),peqr(nr1),fr(nr1),psipr(nr1))
      ALLOCATE(rm(0:m1max,nr1),zm(m1max,nr1))
      READ(in_unit)rho
      READ(in_unit)qr
      READ(in_unit)peqr
      READ(in_unit)
      READ(in_unit)fr
      READ(in_unit)psipr
      READ(in_unit)
      READ(in_unit)
      READ(in_unit)
      READ(in_unit)rm
      READ(in_unit)zm
      call bin_close(in_unit)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ma=nr1-1
      CALL spline_alloc(sq_in,ma,4)
      psio=psipr(nr1)/2
      sq_in%xs=rho**2
      sq_in%fs(:,1)=fr
      sq_in%fs(:,2)=peqr
      sq_in%fs(:,3)=qr
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      ro=rm(0,1)
      zo=0
      dtheta=twopi/mtheta
      CALL bicube_alloc(rz_in,ma,mtheta,2)
      theta=0
      DO itheta=0,mtheta
         rz_in%fs(:,itheta,1)=rm(0,:)
         rz_in%fs(:,itheta,2)=0
         DO m=1,m1max
            rz_in%fs(:,itheta,1)=rz_in%fs(:,itheta,1)
     $           +rm(m,:)*COS(m*theta)
            rz_in%fs(:,itheta,2)=rz_in%fs(:,itheta,2)
     $           +zm(m,:)*SIN(m*theta)
         ENDDO
         theta=theta+dtheta
      ENDDO
c-----------------------------------------------------------------------
c     deallocate local arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(rho,qr,peqr,fr,psipr,rm,zm)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_popov2
c-----------------------------------------------------------------------
c     subprogram 20. read_eq_rtaylor.
c     reads data from Bob Taylor.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_rtaylor

      INTEGER :: mr,mz,ma
      REAL(8) :: ip
      REAL(8), DIMENSION(:), POINTER :: psi,f,p
      REAL(8), DIMENSION(:,:), POINTER :: psig
c-----------------------------------------------------------------------
c     open input file, read scalar data, allocate and read arrays.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND","none")
      READ(in_unit)mr,mz,ma,rmin,rmax,zmin,zmax,ip
      ALLOCATE(psi(ma),f(ma),p(ma),psig(mr,mz))
      READ(in_unit)psi,p,f,psig
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     store 1d data.
c-----------------------------------------------------------------------
      ma=ma-1
      psio=psi(1)
      CALL spline_alloc(sq_in,ma,4)
      sq_in%xs=1-psi/psi(1)
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=p
      sq_in%fs(:,3)=0
c-----------------------------------------------------------------------
c     store 2d data.
c-----------------------------------------------------------------------
      mr=mr-1
      mz=mz-1
      CALL bicube_alloc(psi_in,mr,mz,1)
      psi_in%fs(:,:,1)=psig
c-----------------------------------------------------------------------
c     deallocate local arrays and process equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,f,p,psig)
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_rtaylor
c-----------------------------------------------------------------------
c     subprogram 21. read_eq_dump.
c     reads dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_dump

      INTEGER :: nqty

      LOGICAL, PARAMETER :: diagnose=.TRUE.
      INTEGER :: ix,iy,iqty,iside,ipsi
      REAL(r8) :: f0,f0fac,ffac
      REAL(r8), DIMENSION(:), POINTER :: xfac
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(2x,"nqty",4x,"mx",4x,"my",
     $     5x,"x01",8x,"x02",8x,"y01",8x,"y02"//3i6,1p,4e11.3/)
 20   FORMAT(/2x,"iqty",5x,"xpr1",7x,"xpr2",7x,"ypr1",7x,"ypr2"/)
 30   FORMAT(i6,1p,4e11.3)
 40   FORMAT(/4x,"ix",6x,"xs"/)
 50   FORMAT(i6,1p,e11.3)
 60   FORMAT(/4x,"iy",6x,"ys"/)
 70   FORMAT(/4x,"iy",4x,"ix",4(5x,"fs",i1,3x)/)
 80   FORMAT(2i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     open dump file.
c-----------------------------------------------------------------------
      CALL bin_open(dump_unit,TRIM(eq_filename),
     $     "UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      READ(dump_unit)mpsi,mtheta,jac_type,power_bp,power_b,power_r,
     $     grid_type,psilow,psihigh,ro,zo,psio,q0,qa
c-----------------------------------------------------------------------
c     write 1D splines.
c-----------------------------------------------------------------------
      READ(dump_unit)mpsi,nqty
      CALL spline_alloc(sq,mpsi,nqty)
      READ(dump_unit)sq%xs,sq%fs,
     $     sq%xpower,sq%x0,sq%title,sq%name,sq%periodic
      CALL spline_fit(sq,"extrap")
c-----------------------------------------------------------------------
c     read 2D splines.
c-----------------------------------------------------------------------
      READ(dump_unit)mpsi,mtheta,nqty
      CALL bicube_alloc(rzphi,mpsi,mtheta,nqty)
      READ(dump_unit)rzphi%xs
      READ(dump_unit)rzphi%ys
      READ(dump_unit)rzphi%fs
      READ(dump_unit)rzphi%x0
      READ(dump_unit)rzphi%y0
      READ(dump_unit)rzphi%xpower
      READ(dump_unit)rzphi%ypower
      READ(dump_unit)rzphi%xtitle
      READ(dump_unit)rzphi%ytitle
      READ(dump_unit)rzphi%title
      READ(dump_unit)rzphi%name
      READ(dump_unit)rzphi%periodic
c-----------------------------------------------------------------------
c     restore x powers and fit.
c-----------------------------------------------------------------------
      ALLOCATE(xfac(0:rzphi%mx))
      DO iside=1,2
         DO iqty=1,rzphi%nqty
            IF(rzphi%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(rzphi%xs-rzphi%x0(iside))
     $              **rzphi%xpower(iside,iqty)
               DO iy=0,rzphi%my
                  rzphi%fs(:,iy,iqty)=rzphi%fs(:,iy,iqty)/xfac
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(xfac)
      CALL bicube_fit(rzphi,"extrap","periodic")
c-----------------------------------------------------------------------
c     revise q profile.
c-----------------------------------------------------------------------
      IF(newq0 /= 0)THEN
         f0=sq%fs(0,1)-sq%fs1(0,1)*sq%xs(0)
         f0fac=f0**2*((newq0/q0)**2-one)
         q0=newq0
         DO ipsi=0,mpsi
            ffac=SQRT(1+f0fac/sq%fs(ipsi,1)**2)
            sq%fs(ipsi,1)=sq%fs(ipsi,1)*ffac
            sq%fs(ipsi,4)=sq%fs(ipsi,4)*ffac
            rzphi%fs(ipsi,:,3)=rzphi%fs(ipsi,:,3)*ffac
         ENDDO
         CALL spline_fit(sq,"extrap")
      ENDIF
      qa=sq%fs(mpsi,4)+sq%fs1(mpsi,4)*(one-sq%xs(mpsi))
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL bin_close(dump_unit)
c-----------------------------------------------------------------------
c     diagnose rzphi.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(dump_unit,"rzphi1.out","UNKNOWN")
         WRITE(dump_unit,10)rzphi%nqty,rzphi%mx,rzphi%my,
     $        rzphi%x0,rzphi%y0
         WRITE(dump_unit,20)
         WRITE(dump_unit,30)(iqty,rzphi%xpower(:,iqty),
     $        rzphi%ypower(:,iqty),iqty=1,rzphi%nqty)
         WRITE(dump_unit,20)
         WRITE(dump_unit,40)
         WRITE(dump_unit,50)(ix,rzphi%xs(ix),ix=0,rzphi%mx)
         WRITE(dump_unit,40)
         WRITE(dump_unit,60)
         WRITE(dump_unit,50)(iy,rzphi%ys(iy),iy=0,rzphi%my)
         WRITE(dump_unit,60)
         WRITE(dump_unit,70)(iqty,iqty=1,rzphi%nqty)
         DO iy=0,rzphi%my
            DO ix=0,rzphi%mx
               WRITE(dump_unit,80)iy,ix,rzphi%fs(ix,iy,:)
            ENDDO
            WRITE(dump_unit,70)(iqty,iqty=1,rzphi%nqty)
         ENDDO
         CALL ascii_close(dump_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_dump
c-----------------------------------------------------------------------
c     subprogram 22. read_eq_wdn.
c     reads data from Dave Nystrom's equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_wdn

      INTEGER :: mr,mz,ma
c-----------------------------------------------------------------------
c     open input, file, read scalar data, and allocate arrays.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)mr,mz,ma
      READ(in_unit)rmin,rmax,zmin,zmax
      ma=ma-1
      mr=mr-1
      mz=mz-1
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(psi_in,mr,mz,1)
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(in_unit)psi_in%fs
      READ(in_unit)sq_in%xs
      READ(in_unit)sq_in%fs(:,1)
      READ(in_unit)sq_in%fs(:,2)
      READ(in_unit)sq_in%fs(:,3)
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     revise 1D arrays.
c-----------------------------------------------------------------------
      psio=sq_in%xs(0)
      sq_in%xs=1-sq_in%xs/psio
      sq_in%fs(:,1)=sq_in%fs(:,1)*twopi
      sq_in%fs(:,2)=sq_in%fs(:,2)*mu0
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL direct_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_wdn
c-----------------------------------------------------------------------
c     subprogram 23. read_eq_lez_2.
c     reads data from Leonid E. Zakharov's equilibrium code (by kkim).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_lez_2

      LOGICAL :: file_stat
      INTEGER :: ia,ja,ma,mtau
      REAL(r4), DIMENSION(:), POINTER :: psis,qs,fs,ps
      REAL(r4), DIMENSION(:,:), POINTER :: rg,zg
c-----------------------------------------------------------------------
c     open equilibrium file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eq_filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL program_stop
     $     ("Can't open input file "//TRIM(eq_filename))
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
c-----------------------------------------------------------------------
c     read sizes and allocate input arrays.
c-----------------------------------------------------------------------
      READ(in_unit,*)ma
      ma=ma-1
      READ(in_unit,*)mtau
      mtau=mtau-1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psis(0:ma),qs(0:ma),fs(0:ma),ps(0:ma),
     $     rg(0:ma,0:mtau),zg(0:ma,0:mtau))
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
c-----------------------------------------------------------------------
c     read 2d data.
c-----------------------------------------------------------------------
      DO ia=0,mtau
         DO ja=0,ma
           READ(in_unit,*)rg(ja,ia)
        ENDDO
      ENDDO
      DO ia=0,mtau
         DO ja=0,ma
           READ(in_unit,*)zg(ja,ia)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     read 1d data.
c-----------------------------------------------------------------------
      DO ja=0,ma
         READ(in_unit,*)psis(ja),qs(ja),fs(ja),ps(ja)
      ENDDO
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     revise 1D data.
c-----------------------------------------------------------------------
      psio=psis(ma)
      sq_in%xs=psis/psio
      sq_in%fs(:,1)=fs
      sq_in%fs(:,2)=ps
      sq_in%fs(:,3)=qs
c-----------------------------------------------------------------------
c     revise 2D profiles.
c-----------------------------------------------------------------------
      ro=rg(0,0)
      zo=zg(0,0)
      rz_in%fs(:,:,1)=rg
      rz_in%fs(:,:,2)=zg
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(psis,qs,fs,ps,rg,zg)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_lez_2
c-----------------------------------------------------------------------
c     subprogram 24. read_eq_chease3.
c     reads inverse equilibrium from CHEASE equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_chease3

      INTEGER :: npsi,ipsi,nchi
      REAL(r8) :: b0exp,r0exp,rmag,zmag,tpsio
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND",
     $     convert_type)
      READ(in_unit)npsi,nchi,r0exp,b0exp,rmag,zmag,tpsio
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(norpsi(npsi),psi(npsi),t(npsi),mu0p(npsi),q(npsi))
      ALLOCATE(req(npsi,nchi),zeq(npsi,nchi))
c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      DO ipsi=1,npsi
         READ(in_unit)norpsi(ipsi),psi(ipsi),t(ipsi),mu0p(ipsi),q(ipsi)
         READ(in_unit)req(ipsi,:)
         READ(in_unit)zeq(ipsi,:)
      ENDDO
      CALL ascii_close(in_unit)
      WRITE(*,*) "finish read"
c-----------------------------------------------------------------------
c     copy variables.
c-----------------------------------------------------------------------
      ro=rmag
      zo=zmag
      psio=abs(tpsio)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------      
      CALL spline_alloc(sq_in,npsi-1,4)
      sq_in%xs=norpsi
      sq_in%fs(:,1)=t
      sq_in%fs(:,2)=mu0p
      sq_in%fs(:,3)=q
c      WRITE(*,*) sq_in%fs(:,3)
      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      CALL bicube_alloc(rz_in,npsi-1,nchi-1,2)
      rz_in%fs(:,:,1)=req
      rz_in%fs(:,:,2)=zeq
c      WRITE(*,*) "rz_in2",rz_in%fs(:,:,2)
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(norpsi,psi,t,mu0p,q)
      DEALLOCATE(req,zeq)
      CALL bicube_fit(rz_in,"extrap","periodic")
c      WRITE(*,*) "finish bispline"
      CALL inverse_run
      WRITE(*,*) "finish inverse"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_chease3

      END MODULE read_eq_mod
