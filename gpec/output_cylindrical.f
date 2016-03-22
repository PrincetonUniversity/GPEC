c-----------------------------------------------------------------------
c     GENERALIZED PERTURBED EQUILIBRIUM CODE
c     write various output results of gpec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. output_cylindrical
c      1. xbrzphi
c      2. vsbrzphi
c      3. xbrzphifun
c      4. arzphifun
c-----------------------------------------------------------------------
c     subprogram 0. output_cylindrical.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE output_cylindrical
      USE gpec_response
      USE ipvacuum_mod
      USE gpec_diagnostic
      USE field_mod
      USE netcdf

      IMPLICIT NONE

      ! harvest variables
      INCLUDE 'harvest_lib.inc77'
      INTEGER  :: ierr
      INTEGER  :: hint = 0
      REAL(r8) :: hdbl = 0
      CHARACTER(LEN=16)    :: hkey
      CHARACTER(LEN=50000) :: hnml
      CHARACTER(LEN=65507) :: hlog
      CHARACTER, PARAMETER :: nul = char(0)

      ! netcdf ids
      INTEGER :: mncid,fncid,cncid
      CHARACTER(64) :: mncfile,fncfile,cncfile

      ! module wide output variables
      LOGICAL :: singcoup_set = .FALSE.
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: singcoup1mat,
     $   w1v,w2v,w3v,t1v,t2v,t3v,fldflxmat

      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
      SUBROUTINE xbrzphi(egnum,xspmn,nr,nz,bnimn,bnomn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnimn,bnomn

      INTEGER :: i,j,k,l,ipert,iindex,np
      REAL(r8) :: mid,btlim,rlim,ileft,delr,delz,cha,chb,chc,chd,
     $   rij,zij,tij,t11,t12,t21,t22,t33
      COMPLEX(r8) :: xwp,bwp,xwt,bwt,xvz,bvz

      INTEGER :: r_id,z_id,i_id,xr_id,xz_id,xp_id,br_id,bz_id,bp_id,
     $   bre_id,bze_id,bpe_id,brp_id,bzp_id,bpp_id

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz,ebr,ebz,ebp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: xrr,xrz,xrp,brr,brz,brp,
     $     bpr,bpz,bpp,vbr,vbz,vbp,vpbr,vpbz,vpbp,vvbr,vvbz,vvbp,
     $     btr,btz,btp,vcbr,vcbz,vcbp,xcr,xcz,xcp,bcr,bcz,bcp

      REAL(r8), DIMENSION(:), POINTER :: chex,chey
      COMPLEX(r8), DIMENSION(:,:), POINTER :: chear,cheaz,
     $     chxar,chxaz

c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      ! initialization
      ebr = 0
      ebz = 0
      ebp = 0
      xrr = 0
      xrz = 0
      xrp = 0
      xcr = 0
      xcz = 0
      xcp = 0
      brr = 0
      brz = 0
      brp = 0
      btr = 0
      btz = 0
      btp = 0
      bcr = 0
      bcz = 0
      bcp = 0
      bpr = 0
      bpz = 0
      bpp = 0
      vbr = 0
      vbz = 0
      vbp = 0
      vcbr = 0
      vcbz = 0
      vcbp = 0
      vpbr = 0
      vpbz = 0
      vpbp = 0
      vvbr = 0
      vvbz = 0
      vvbp = 0

      CALL idcon_build(egnum,xspmn)

      IF (eqbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing equilibrium magnetic fields"
         IF(timeit) CALL gpec_timer(-2)
         ! evaluate f value for vacuum
         mid = 0.0
         CALL spline_eval(sq,psilim,0)
         CALL bicube_eval(rzphi,psilim,mid,0)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*rzphi%f(2)
         rlim=ro+rfac*COS(eta)
         btlim = abs(sq%f(1))/(twopi*rlim)
         DO i=0,nr
            DO j=0,nz
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),1)
               ebr(i,j) = -psi_in%fy(1)/gdr(i,j)*psio
               ebz(i,j) = psi_in%fx(1)/gdr(i,j)*psio
               IF (gdl(i,j) == 1) THEN  
                  CALL spline_eval(sq,gdpsi(i,j),0)
                  ebp(i,j) = abs(sq%f(1))/(twopi*gdr(i,j))
               ELSE
                  ebp(i,j) = btlim*rlim/gdr(i,j)  
               ENDIF
            ENDDO   
         ENDDO

         IF(ipd>0)THEN
            ebr=-ebr
            ebz=-ebz
         ENDIF
         IF(btd<0)ebp=-ebp
         IF(timeit) CALL gpec_timer(2)
      ENDIF

      CALL peq_alloc

      IF(verbose) WRITE(*,*)"Mapping fields to cylindrical coordinates"
      IF(timeit) CALL gpec_timer(-2)

      IF (brzphi_flag .OR. xrzphi_flag) THEN
         DO i=0,nr
         iindex = FLOOR(REAL(i,8)/FLOOR((nr-1)/10.0))*10
         ileft = REAL(i,8)/FLOOR((nr-1)/10.0)*10-iindex
         IF ((i /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a10)')
     $        "volume = ",iindex,"% mappings"
            DO j=0,nz
               IF (gdl(i,j)==1) THEN
                  CALL peq_sol(gdpsi(i,j))
                  CALL peq_contra(gdpsi(i,j))
                  CALL peq_cova(gdpsi(i,j))
                  ! compute matric tensor components.
                  CALL bicube_eval(rzphi,gdpsi(i,j),gdthe(i,j),1)
                  rfac=SQRT(rzphi%f(1))
                  eta=twopi*(gdthe(i,j)+rzphi%f(2))
                  rij=ro+rfac*COS(eta)
                  zij=zo+rfac*SIN(eta)
                  jac=rzphi%f(4)
                  v(1,1)=rzphi%fx(1)/(2*rfac)
                  v(1,2)=rzphi%fx(2)*twopi*rfac
                  v(2,1)=rzphi%fy(1)/(2*rfac)
                  v(2,2)=(1+rzphi%fy(2))*twopi*rfac
                  v(3,3)=twopi*rij
                  t11=cos(eta)*v(1,1)-sin(eta)*v(1,2)
                  t12=cos(eta)*v(2,1)-sin(eta)*v(2,2)
                  t21=sin(eta)*v(1,1)+cos(eta)*v(1,2)
                  t22=sin(eta)*v(2,1)+cos(eta)*v(2,2)
                  t33=-1.0/v(3,3)
                  ! three vector components.
                  xwp = SUM(xwp_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bwp = SUM(bwp_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xwt = SUM(xwt_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bwt = SUM(bwt_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xvz = SUM(xvz_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bvz = SUM(bvz_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xrr(i,j) = (t11*xwp+t12*xwt)/jac
                  brr(i,j) = (t11*bwp+t12*bwt)/jac
                  xrz(i,j) = (t21*xwp+t22*xwt)/jac
                  brz(i,j) = (t21*bwp+t22*bwt)/jac
                  xrp(i,j) = t33*xvz
                  brp(i,j) = t33*bvz
                  ! Correction for helicity
                  IF(helicity>0)THEN
                     xrr(i,j) = CONJG(xrr(i,j))
                     brr(i,j) = CONJG(brr(i,j))
                     xrz(i,j) = CONJG(xrz(i,j))
                     brz(i,j) = CONJG(brz(i,j))
                     xrp(i,j) = CONJG(xrp(i,j))
                     brp(i,j) = CONJG(brp(i,j))
                  ENDIF
                  ! machine toroidal angle
                  xrr(i,j) = xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brr(i,j) = brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrz(i,j) = xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brz(i,j) = brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrp(i,j) = xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brp(i,j) = brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      IF(timeit) CALL gpec_timer(2)

      CALL peq_dealloc

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)
     $      "Computing vacuum fields by surface currents"
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
         IF (helicity<0) THEN
            vbr=CONJG(vbr)
            vbz=CONJG(vbz)
            vbp=CONJG(vbp)
         ENDIF
         DO i=0,nr
            DO j=0,nz
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  brr(i,j)=vbr(i,j)
                  brz(i,j)=vbz(i,j)
                  brp(i,j)=vbp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
         IF(timeit) CALL gpec_timer(2)
      ENDIF
      
      IF (divzero_flag) THEN
         CALL peq_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF
      IF (div_flag) THEN
         CALL diagnose_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF

      IF (brzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF

         IF (coil_flag) THEN
            IF(verbose) WRITE(*,*)"Computing vacuum fields by coils"
            np=nn*16
            CALL field_bs_rzphi(nr,nz,np,gdr,gdz,vcbr,vcbz,vcbp)
         ENDIF
         
         DO i=0,nr
            DO j=0,nz                  
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  bpr(i,j)=vpbr(i,j)
                  bpz(i,j)=vpbz(i,j)
                  bpp(i,j)=vpbp(i,j)
                  btr(i,j)=vpbr(i,j)+vcbr(i,j)
                  btz(i,j)=vpbz(i,j)+vcbz(i,j)
                  btp(i,j)=vpbp(i,j)+vcbp(i,j)
               ELSE
                  bpr(i,j)=brr(i,j)-vcbr(i,j)
                  bpz(i,j)=brz(i,j)-vcbz(i,j)
                  bpp(i,j)=brp(i,j)-vcbp(i,j)
                  btr(i,j)=brr(i,j)
                  btz(i,j)=brz(i,j)
                  btp(i,j)=brp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF

      IF (chebyshev_flag) THEN
         IF(verbose) WRITE(*,*)"Computing chebyshev for xbrzphi"
         ALLOCATE(chex(0:nr),chey(0:nz),
     $        chear(0:nchr,0:nchz),cheaz(0:nchr,0:nchz),
     $        chxar(0:nchr,0:nchz),chxaz(0:nchr,0:nchz))
         delr = gdr(1,0)-gdr(0,0)
         delz = gdz(0,1)-gdz(0,0)
         cha = (gdr(nr,0)-gdr(0,0))/2.0+delr
         chb = (gdr(nr,0)+gdr(0,0))/2.0
         chc = (gdz(0,nz)-gdz(0,0))/2.0+delz
         chd = (gdz(0,nz)+gdz(0,0))/2.0
         chex = (gdr(:,0)-chb)/cha
         chey = (gdz(0,:)-chd)/chc
         chear=0
         cheaz=0
         chxar=0
         chxaz=0

         DO i=0,nchr
            DO j=0,nchz
               DO k=0,nr
                  DO l=0,nz
                     chear(i,j)=chear(i,j)+btr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     cheaz(i,j)=cheaz(i,j)+btz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                     chxar(i,j)=chxar(i,j)+xrr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     chxaz(i,j)=chxaz(i,j)+xrz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                  ENDDO
               ENDDO
               chear(i,j)=chear(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               cheaz(i,j)=cheaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxar(i,j)=chxar(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxaz(i,j)=chxaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
            ENDDO
         ENDDO
         chear(0,:)=chear(0,:)/2.0
         cheaz(0,:)=cheaz(0,:)/2.0        
         chear(:,0)=chear(:,0)/2.0
         cheaz(:,0)=cheaz(:,0)/2.0  
         chxar(0,:)=chxar(0,:)/2.0
         chxaz(0,:)=chxaz(0,:)/2.0        
         chxar(:,0)=chxar(:,0)/2.0
         chxaz(:,0)=chxaz(:,0)/2.0   

         CALL ascii_open(out_unit,"gpec_brzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for brzphi field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,'(4(1x,a4,f12.3))')
     $        "a =",cha,"b =",chb,"c =",chc,"d =",chd
         WRITE(out_unit,*)   
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chear","cheaz",
     $        "real(chear)","imag(chear)","real(cheaz)","imag(cheaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $              REAL(chear(i,j)),AIMAG(chear(i,j)),
     $              REAL(cheaz(i,j)),AIMAG(cheaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

         CALL ascii_open(out_unit,"gpec_xrzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for displacement"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chxar","chxaz",
     $        "real(chxar)","imag(chxar)","real(chxaz)","imag(chxaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $              REAL(chxar(i,j)),AIMAG(chxar(i,j)),
     $              REAL(chxaz(i,j)),AIMAG(chxaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

         IF(verbose) WRITE(*,*)"Recontructing xbrzphi by chebyshev"

         DO i=0,nr
            DO j=0,nz
               DO k=0,nchr
                  DO l=0,nchz
                     bcr(i,j)=bcr(i,j)+chear(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcz(i,j)=bcz(i,j)+cheaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcp(i,j)=bcp(i,j)+chear(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+cheaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                     xcr(i,j)=xcr(i,j)+chxar(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcz(i,j)=xcz(i,j)+chxaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcp(i,j)=xcp(i,j)+chxar(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+chxaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                  ENDDO
               ENDDO
               bcp(i,j)=bcp(i,j)+bcr(i,j)
               xcp(i,j)=xcp(i,j)+xcr(i,j)
            ENDDO
         ENDDO
         bcp=bcp*(-ifac/nn)
         xcp=xcp*(-ifac/nn)

         DEALLOCATE(chex,chey,chear,cheaz,chxar,chxaz)

      ENDIF

      IF (pbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF
      ENDIF

      IF (vvbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)
     $      "Computing vacuum fields without plasma response"
         CALL ipvacuum_bnormal(psilim,bnimn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vvbr,vvbz,vvbp)
         IF (helicity<0) THEN
            vvbr=CONJG(vvbr)
            vvbz=CONJG(vvbz)
            vvbp=CONJG(vvbp)
         ENDIF
      ENDIF  
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      ! append to netcdf file once this is (mstep,mpert)
      IF(debug_flag) PRINT *,"Opening "//TRIM(cncfile)
      CALL check( nf90_open(cncfile,nf90_write,cncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(cncid,"i",i_id) )
      CALL check( nf90_inq_dimid(cncid,"R",r_id) )
      CALL check( nf90_inq_dimid(cncid,"z",z_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(cncid))
      CALL check( nf90_def_var(cncid, "b_r_equil", nf90_double,
     $            (/r_id,z_id/),bre_id) )
      CALL check( nf90_put_att(cncid,bre_id,"long_name",
     $            "Radial Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bre_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z_equil", nf90_double,
     $            (/r_id,z_id/),bze_id) )
      CALL check( nf90_put_att(cncid,bze_id,"long_name",
     $            "Vertical Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bze_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t_equil", nf90_double,
     $            (/r_id,z_id/),bpe_id) )
      CALL check( nf90_put_att(cncid,bpe_id,"long_name",
     $            "Toroidal Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bpe_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_r_plas", nf90_double,
     $            (/r_id,z_id,i_id/),brp_id) )
      CALL check( nf90_put_att(cncid,brp_id,"long_name",
     $            "Radial Plasma Field") )
      CALL check( nf90_put_att(cncid,brp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z_plas", nf90_double,
     $            (/r_id,z_id,i_id/),bzp_id) )
      CALL check( nf90_put_att(cncid,bzp_id,"long_name",
     $            "Vertical Plasma Field") )
      CALL check( nf90_put_att(cncid,bzp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t_plas", nf90_double,
     $            (/r_id,z_id,i_id/),bpp_id) )
      CALL check( nf90_put_att(cncid,bpp_id,"long_name",
     $            "Toroidal Plasma Field") )
      CALL check( nf90_put_att(cncid,bpp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_r", nf90_double,
     $            (/r_id,z_id,i_id/),br_id) )
      CALL check( nf90_put_att(cncid,br_id,"long_name",
     $            "Radial Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,br_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z", nf90_double,
     $            (/r_id,z_id,i_id/),bz_id) )
      CALL check( nf90_put_att(cncid,bz_id,"long_name",
     $            "Vertical Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,bz_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t", nf90_double,
     $            (/r_id,z_id,i_id/),bp_id) )
      CALL check( nf90_put_att(cncid,bp_id,"long_name",
     $            "Toroidal Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,bp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "xi_r", nf90_double,
     $            (/r_id,z_id,i_id/),xr_id) )
      CALL check( nf90_put_att(cncid,xr_id,"long_name",
     $            "Radial Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,xr_id,"units","m") )
      CALL check( nf90_def_var(cncid, "xi_z", nf90_double,
     $            (/r_id,z_id,i_id/),xz_id) )
      CALL check( nf90_put_att(cncid,xz_id,"long_name",
     $            "Vertical Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,bz_id,"units","m") )
      CALL check( nf90_def_var(cncid, "xi_t", nf90_double,
     $            (/r_id,z_id,i_id/),xp_id) )
      CALL check( nf90_put_att(cncid,xp_id,"long_name",
     $            "Toroidal Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,bp_id,"units","m") )
      CALL check( nf90_enddef(cncid) )
      IF(debug_flag) PRINT *,"  Writting variables"
      CALL check( nf90_put_var(cncid,bre_id,ebr) )
      CALL check( nf90_put_var(cncid,bze_id,ebz) )
      CALL check( nf90_put_var(cncid,bpe_id,ebp) )
      CALL check( nf90_put_var(cncid,br_id,RESHAPE((/REAL(btr),
     $             AIMAG(btr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bz_id,RESHAPE((/REAL(btz),
     $             AIMAG(btz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bp_id,RESHAPE((/REAL(btp),
     $             AIMAG(btp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,brp_id,RESHAPE((/REAL(bpr),
     $             AIMAG(bpr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bzp_id,RESHAPE((/REAL(bpz),
     $             AIMAG(bpz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bpp_id,RESHAPE((/REAL(bpp),
     $             AIMAG(bpp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xr_id,RESHAPE((/REAL(xrr),
     $             AIMAG(xrr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xz_id,RESHAPE((/REAL(xrz),
     $             AIMAG(xrz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xp_id,RESHAPE((/REAL(xrp),
     $             AIMAG(xrp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_close(cncid) )
      IF(debug_flag) PRINT *,"Closed "//TRIM(cncfile)


      IF (eqbrzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_eqbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_EQBRZPHI: Eq. b field in rzphi grid"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,5(a17))')"l","r","z","eb_r",
     $        "eb_z","eb_phi"
      
         DO i=0,nr
            DO j=0,nz 
               WRITE(out_unit,'(1x,I2,5(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              ebr(i,j),ebz(i,j),ebp(i,j)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (brzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_brzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BRZPHI: Total perturbed field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(btr(i,j)),AIMAG(btr(i,j)),
     $              REAL(btz(i,j)),AIMAG(btz(i,j)),
     $              REAL(btp(i,j)),AIMAG(btp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"gpec_chebrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CEHBRZPHI: Total perturbed field"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(bcr(i,j)),AIMAG(bcr(i,j)),
     $                 REAL(bcz(i,j)),AIMAG(bcz(i,j)),
     $                 REAL(bcp(i,j)),AIMAG(bcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         CALL ascii_open(out_unit,"gpec_pbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PBRZPHI: Perturbed field by plasma"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(bpr(i,j)),AIMAG(bpr(i,j)),
     $              REAL(bpz(i,j)),AIMAG(bpz(i,j)),
     $              REAL(bpp(i,j)),AIMAG(bpp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (coil_flag) THEN
            CALL ascii_open(out_unit,"gpec_cbrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CBRZPHI: External field by coils"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(vcbr(i,j)),AIMAG(vcbr(i,j)),
     $                 REAL(vcbz(i,j)),AIMAG(vcbz(i,j)),
     $                 REAL(vcbp(i,j)),AIMAG(vcbp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_vbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBRZPHI: Total perturbed field "//
     $        "by surface currents"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(brr(i,j)),AIMAG(brr(i,j)),
     $              REAL(brz(i,j)),AIMAG(brz(i,j)),
     $              REAL(brp(i,j)),AIMAG(brp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

      ENDIF

      IF (xrzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_xrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XRZPHI: Displacement"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(xi_r)","imag(xi_r)","real(xi_z)","imag(xi_z)",
     $        "real(xi_phi)","imag(xi_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(xrr(i,j)),AIMAG(xrr(i,j)),
     $              REAL(xrz(i,j)),AIMAG(xrz(i,j)),
     $              REAL(xrp(i,j)),AIMAG(xrp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"gpec_chexrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CEHXRZPHI: Displacement"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(x_r)","imag(x_r)","real(x_z)","imag(x_z)",
     $           "real(x_phi)","imag(x_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(xcr(i,j)),AIMAG(xcr(i,j)),
     $                 REAL(xcz(i,j)),AIMAG(xcz(i,j)),
     $                 REAL(xcp(i,j)),AIMAG(xcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

      IF (pbrzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_vpbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VPBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              vgdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(vpbr(i,j)),AIMAG(vpbr(i,j)),
     $              REAL(vpbz(i,j)),AIMAG(vpbz(i,j)),
     $              REAL(vpbp(i,j)),AIMAG(vpbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
      ENDIF

      IF (vvbrzphi_flag) THEN
         CALL ascii_open(out_unit,"gpec_vvbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VVBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              vgdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(vvbr(i,j)),AIMAG(vvbr(i,j)),
     $              REAL(vvbz(i,j)),AIMAG(vvbz(i,j)),
     $              REAL(vvbp(i,j)),AIMAG(vvbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
      ENDIF

      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag) DEALLOCATE(gdr,gdz,gdl,gdpsi,gdthe,gdphi)
      CALL peq_rzpgrid(nr,nz,psixy) ! reset the grid
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbrzphi
c-----------------------------------------------------------------------
c     subprogram 2. vsbrzphi.
c     write brzphi components restored by removing shielding currents.
c-----------------------------------------------------------------------
      SUBROUTINE vsbrzphi(snum,nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: snum,nr,nz

      INTEGER :: i,j,ipert,iindex
      REAL(r8) :: mid,bt0,ileft

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: vbr,vbz,vbp
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      vbr = 0
      vbz = 0
      vbp = 0

      IF (snum<10) THEN
         WRITE(UNIT=ss,FMT='(I1)')snum
         ss=TRIM(ADJUSTL(ss))
      ELSE
         WRITE(UNIT=ss,FMT='(I2)')snum
      ENDIF

      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vacuum fields by "//
     $     TRIM(ss)//"th resonant field"
      CALL ipvacuum_bnormal(singtype(snum)%psifac,
     $     singbno_mn(:,snum),nr,nz)
      CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $     nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
      IF (helicity<0) THEN
         vbr=CONJG(vbr)
         vbz=CONJG(vbz)
         vbp=CONJG(vbp)
      ENDIF
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"gpec_vsbrzphi_n"//
     $     TRIM(sn)//"_s"//TRIM(ss)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_VSBRZPHI: Vacuum field in rzphi grid by "//
     $     TRIM(ss)//"th resonant field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $     "real(vsb_r)","imag(vsb_r)","real(vsb_z)","imag(vsb_z)",
     $     "real(vsb_phi)","imag(vsb_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $           vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $           REAL(vbr(i,j)),AIMAG(vbr(i,j)),
     $           REAL(vbz(i,j)),AIMAG(vbz(i,j)),
     $           REAL(vbp(i,j)),AIMAG(vbp(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE vsbrzphi
c-----------------------------------------------------------------------
c     subprogram 3. xbrzphifun
c-----------------------------------------------------------------------
      SUBROUTINE xbrzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,dphi,t11,t12,t21,t22,t33
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      IF(verbose) WRITE(*,*)"Computing x and b rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL peq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a29)')
     $        "volume = ",iindex,"% xi and b rzphi computations"
         CALL peq_sol(psifac(istep))
         CALL peq_contra(psifac(istep))
         CALL peq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(3,3)=twopi*rs(istep,itheta)
            t11(itheta)=cos(eta)*v(1,1)-sin(eta)*v(1,2)
            t12(itheta)=cos(eta)*v(2,1)-sin(eta)*v(2,2)
            t21(itheta)=sin(eta)*v(1,1)+cos(eta)*v(1,2)
            t22(itheta)=sin(eta)*v(2,1)+cos(eta)*v(2,2)
            t33(itheta)=-1.0/v(3,3)
         ENDDO

         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         xrr_fun(istep,:)=(t11*xwp_fun+t12*xwt_fun)/jacs
         brr_fun(istep,:)=(t11*bwp_fun+t12*bwt_fun)/jacs
         xrz_fun(istep,:)=(t21*xwp_fun+t22*xwt_fun)/jacs
         brz_fun(istep,:)=(t21*bwp_fun+t22*bwt_fun)/jacs
         xrp_fun(istep,:)=t33*xvz_fun
         brp_fun(istep,:)=t33*bvz_fun

         xrr_fun(istep,:)=xrr_fun(istep,:)*EXP(ifac*nn*dphi)
         brr_fun(istep,:)=brr_fun(istep,:)*EXP(ifac*nn*dphi)
         xrz_fun(istep,:)=xrz_fun(istep,:)*EXP(ifac*nn*dphi)
         brz_fun(istep,:)=brz_fun(istep,:)*EXP(ifac*nn*dphi)
         xrp_fun(istep,:)=xrp_fun(istep,:)*EXP(ifac*nn*dphi)
         brp_fun(istep,:)=brp_fun(istep,:)*EXP(ifac*nn*dphi)         

      ENDDO

      IF(helicity>0)THEN
         xrr_fun=CONJG(xrr_fun)
         xrz_fun=CONJG(xrz_fun)
         xrp_fun=CONJG(xrp_fun)
         brr_fun=CONJG(brr_fun)
         brz_fun=CONJG(brz_fun)
         brp_fun=CONJG(brp_fun)
      ENDIF

      CALL peq_dealloc

      CALL ascii_open(out_unit,"gpec_xbrzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_XBRZPHI_FUN: "//
     $     "Rzphi components of displacement and field in functions"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $     "real(xp)","imag(xp)","real(br)","imag(br)",
     $     "real(bz)","imag(bz)","real(bp)","imag(bp)"

      DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(es17.8e3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(xrr_fun(istep,itheta)),
     $           AIMAG(xrr_fun(istep,itheta)),
     $           REAL(xrz_fun(istep,itheta)),
     $           AIMAG(xrz_fun(istep,itheta)),
     $           REAL(xrp_fun(istep,itheta)),
     $           AIMAG(xrp_fun(istep,itheta)),
     $           REAL(brr_fun(istep,itheta)),
     $           AIMAG(brr_fun(istep,itheta)),
     $           REAL(brz_fun(istep,itheta)),
     $           AIMAG(brz_fun(istep,itheta)),
     $           REAL(brp_fun(istep,itheta)),
     $           AIMAG(brp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbrzphifun
c-----------------------------------------------------------------------
c     subprogram 4. arzphifun
c-----------------------------------------------------------------------
      SUBROUTINE arzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xss_fun,
     $     ear,eat,eap,arr,art,arp
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: ear_fun,eaz_fun,
     $     eap_fun,arr_fun,arz_fun,arp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vector potential rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL peq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a37)')
     $        "volume = ",iindex,"% vector potential rzphi computations"
         CALL peq_sol(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xss_fun,mthsurf,xss_mn)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(istep,itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(2,1)=-rzphi%fx(2)*twopi**2*rs(istep,itheta)*rfac/jac
            w(2,2)=rzphi%fx(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(3,1)=rzphi%fx(3)*2*rfac
            w(3,2)=rzphi%fy(3)/(twopi*rfac)
            w(3,3)=1.0/(twopi*rs(istep,itheta))
            ear(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,1)-w(3,1))
            eat(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,2)-w(3,2))
            eap(itheta)=-chi1*psifac(istep)*w(3,3)
            ear_fun(istep,itheta)=ear(itheta)*cos(eta)-
     $           eat(itheta)*sin(eta)
            eaz_fun(istep,itheta)=ear(itheta)*sin(eta)+
     $           eat(itheta)*cos(eta)
            eap_fun(istep,itheta)=-eap(itheta)            
            arr(itheta)=xss_fun(itheta)*w(1,1)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,1)+xsp_fun(itheta)*w(3,1))
            art(itheta)=xss_fun(itheta)*w(1,2)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,2)+xsp_fun(itheta)*w(3,2))
            arp(itheta)=chi1*xsp_fun(itheta)*w(3,3)
            arr_fun(istep,itheta)=arr(itheta)*cos(eta)-
     $           art(itheta)*sin(eta)
            arz_fun(istep,itheta)=arr(itheta)*sin(eta)+
     $           art(itheta)*cos(eta)
            arp_fun(istep,itheta)=-arp(itheta)
         ENDDO
         ear_fun(istep,:)=ear_fun(istep,:)*EXP(ifac*nn*dphi)
         eaz_fun(istep,:)=eaz_fun(istep,:)*EXP(ifac*nn*dphi)
         eap_fun(istep,:)=eap_fun(istep,:)*EXP(ifac*nn*dphi)
         arr_fun(istep,:)=arr_fun(istep,:)*EXP(ifac*nn*dphi)
         arz_fun(istep,:)=arz_fun(istep,:)*EXP(ifac*nn*dphi)
         arp_fun(istep,:)=arp_fun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      IF(helicity>0)THEN
         ear_fun=CONJG(ear_fun)
         eaz_fun=CONJG(eaz_fun)
         eap_fun=CONJG(eap_fun)
         arr_fun=CONJG(arr_fun)
         arz_fun=CONJG(arz_fun)
         arp_fun=CONJG(arp_fun)
      ENDIF

      CALL peq_dealloc

      CALL ascii_open(out_unit,"gpec_arzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_ARZPHI_FUN: "//
     $     "Rzphi components of vector potential in functions"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(ear)","imag(ear)","real(eaz)","imag(eaz)",
     $     "real(eap)","imag(eap)","real(arr)","imag(arr)",
     $     "real(arz)","imag(arz)","real(arp)","imag(arp)"

      DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(es17.8e3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(ear_fun(istep,itheta)),
     $           AIMAG(ear_fun(istep,itheta)),
     $           REAL(eaz_fun(istep,itheta)),
     $           AIMAG(eaz_fun(istep,itheta)),
     $           REAL(eap_fun(istep,itheta)),
     $           AIMAG(eap_fun(istep,itheta)),
     $           REAL(arr_fun(istep,itheta)),
     $           AIMAG(arr_fun(istep,itheta)),
     $           REAL(arz_fun(istep,itheta)),
     $           AIMAG(arz_fun(istep,itheta)),
     $           REAL(arp_fun(istep,itheta)),
     $           AIMAG(arp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE arzphifun
c-----------------------------------------------------------------------
c     subprogram 5. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE check(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c-----------------------------------------------------------------------
c     subprogram 6. init_netcdf.
c     Initialize the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE init_netcdf
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER:: i,midid,mmdid,medid,mvdid,msdid,mtdid,
     $   fidid,fmdid,fpdid,ftdid,cidid,crdid,czdid,clvid,
     $   mivid,mmvid,mevid,mvvid,msvid,mtvid,fivid,fmvid,fpvid,ftvid,
     $   civid,crvid,czvid,id,fileids(3)
      INTEGER, DIMENSION(mpert) :: mmodes
c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *,"Initializing NETCDF files"
      mmodes = (/(i,i=1,mpert)/)
      mncfile = "gpec_control_output_n"//TRIM(sn)//".nc"
      fncfile = "gpec_profile_output_n"//TRIM(sn)//".nc"
      cncfile = "gpec_cylindrical_output_n"//TRIM(sn)//".nc"
c-----------------------------------------------------------------------
c     open files
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Creating netcdf files"
      CALL check( nf90_create(mncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=mncid) )
      CALL check( nf90_create(fncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=fncid) )
      CALL check( nf90_create(cncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=cncid) )
c-----------------------------------------------------------------------
c     define global file attributes
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining modal netcdf globals"
      CALL check( nf90_put_att(mncid,nf90_global,"title",
     $     "GPEC outputs in Fourier or alternate modal bases"))
      CALL check( nf90_def_dim(mncid,"i",2,       midid) )
      CALL check( nf90_def_var(mncid,"i",nf90_int,midid,mivid) )
      CALL check( nf90_def_dim(mncid,"m",lmpert,  mmdid) )
      CALL check( nf90_def_var(mncid,"m",nf90_int,mmdid,mmvid) )
      CALL check( nf90_def_dim(mncid,"mode",mpert,   medid) )
      CALL check( nf90_def_var(mncid,"mode",nf90_int,medid,mevid))
      CALL check( nf90_def_dim(mncid,"mode_SC",msing,   msdid) )
      CALL check( nf90_def_var(mncid,"mode_SC",nf90_int,msdid,msvid))
      CALL check( nf90_def_dim(mncid,"theta",mthsurf+1,  mtdid) )
      CALL check( nf90_def_var(mncid,"theta",nf90_double,mtdid,mtvid) )
      CALL check( nf90_put_att(mncid,nf90_global,"Jacobian",jac_out))
      CALL check( nf90_put_att(mncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(mncid,nf90_global,"psi_n_lim",psilim))

      IF(debug_flag) PRINT *," - Defining flux netcdf globals"
      CALL check( nf90_put_att(fncid,nf90_global,"title",
     $     "GPEC outputs in magnetic coordinate systems"))
      CALL check( nf90_def_dim(fncid,"i",2,       fidid) )
      CALL check( nf90_def_var(fncid,"i",nf90_int,fidid,fivid) )
      CALL check( nf90_def_dim(fncid,"m",lmpert,  fmdid) )
      CALL check( nf90_def_var(fncid,"m",NF90_INT,fmdid,fmvid) )
      CALL check( nf90_def_dim(fncid,"psi_n",mstep,    fpdid) )
      CALL check( nf90_def_var(fncid,"psi_n",nf90_double,fpdid,fpvid) )
      CALL check( nf90_def_dim(fncid,"theta",mthsurf+1,  ftdid) )
      CALL check( nf90_def_var(fncid,"theta",nf90_double,ftdid,ftvid) )
      CALL check( nf90_put_att(fncid,nf90_global,"Jacobian",jac_type))
      CALL check( nf90_put_att(fncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(fncid,nf90_global,"psi_n_lim",psilim))

      IF(debug_flag) PRINT *," - Defining cylindrical netcdf globals"
      CALL check( nf90_put_att(cncid,nf90_global,"title",
     $     "GPEC outputs in (R,z) coordinates"))
      CALL check( nf90_def_dim(cncid,"i",2,       cidid) )
      CALL check( nf90_def_var(cncid,"i",nf90_int,cidid,civid) )
      CALL check( nf90_def_dim(cncid,"R",nr+1,crdid) )
      CALL check( nf90_def_var(cncid,"R",nf90_double,crdid,crvid) )
      CALL check( nf90_put_att(cncid,crvid,"units","m") )
      CALL check( nf90_def_dim(cncid,"z",nz+1,czdid) )
      CALL check( nf90_def_var(cncid,"z",nf90_double,czdid,czvid) )
      CALL check( nf90_put_att(cncid,czvid,"units","m") )
      CALL check( nf90_def_var(cncid,"l",nf90_double,
     $                         (/crdid,czdid/),clvid) )

      ! Add common global attributes
      fileids = (/mncid,fncid,cncid/)
      DO i=1,3
         id = fileids(i)
         CALL check( nf90_put_att(id,nf90_global,"shot",INT(shotnum)) )
         CALL check( nf90_put_att(id,nf90_global,"time",INT(shottime)) )
         CALL check( nf90_put_att(id,nf90_global,"machine",machine) )
         CALL check( nf90_put_att(id,nf90_global,"n",nn) )
         CALL check( nf90_put_att(id,nf90_global,"version",version))
      ENDDO

      CALL check( nf90_enddef(mncid) )
      CALL check( nf90_enddef(fncid) )
      CALL check( nf90_enddef(cncid) )
c-----------------------------------------------------------------------
c     set dimensional variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Putting coordinates in control netcdfs"
      CALL check( nf90_put_var(mncid,mivid,(/0,1/)) )
      CALL check( nf90_put_var(mncid,mmvid,lmfac) )
      CALL check( nf90_put_var(mncid,mevid,mmodes) )
      CALL check( nf90_put_var(mncid,msvid,mmodes(:msing)) )
      CALL check( nf90_put_var(mncid,mtvid,theta) )

      IF(debug_flag) PRINT *," - Putting coordinates in flux netcdfs"
      CALL check( nf90_put_var(fncid,fivid,(/0,1/)) )
      CALL check( nf90_put_var(fncid,fmvid,lmfac) )
      CALL check( nf90_put_var(fncid,fpvid,psifac(1:)) )
      CALL check( nf90_put_var(fncid,ftvid,theta) )

      IF(debug_flag) PRINT *," - Putting coordinates in cyl. netcdfs"
      CALL check( nf90_put_var(cncid,civid,(/0,1/)) )
      CALL check( nf90_put_var(cncid,crvid,gdr(:,0)) )
      CALL check( nf90_put_var(cncid,czvid,gdz(0,:)) )
      CALL check( nf90_put_var(cncid,clvid,gdl(:,:)) )

      IF(debug_flag) PRINT *," - Closing netcdf files"
      CALL check( nf90_close(mncid) )
      CALL check( nf90_close(fncid) )
      CALL check( nf90_close(cncid) )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE init_netcdf
c-----------------------------------------------------------------------
c     subprogram 17. close_netcdf.
c     Close the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE close_netcdf
c-----------------------------------------------------------------------
c     close files
c-----------------------------------------------------------------------
      CALL check( nf90_close(mncid) )
      CALL check( nf90_close(fncid) )
      CALL check( nf90_close(cncid) )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE close_netcdf

      END MODULE output_cylindrical
