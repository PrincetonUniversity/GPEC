c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     manipulate vacuum code separately.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. gpvacuum_mod
c      1. gpvacuum_arbsurf
c      2. gpvacuum_flxsurf
c      3. gpvacuum_bnormal
c      4. gpvacuum_ideal_mutuals
c      5. gpvacuum_ahgwrite
c-----------------------------------------------------------------------
c     subprogram 0. gpvacuum_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gpvacuum_mod
      USE gpeq_mod
      USE vacuum_mod, ONLY: mscvac, mscfld

      IMPLICIT NONE

      REAL(r8), PRIVATE :: prad
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gpvacuum_arbsurf.
c     compute arbitrary surface inductance.
c-----------------------------------------------------------------------
      SUBROUTINE gpvacuum_arbsurf(majr,minr,psi,dist,surfile)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN), OPTIONAL :: majr,minr,psi,dist
      CHARACTER(128), INTENT(IN), OPTIONAL :: surfile

      INTEGER :: vmtheta,vmlow,vmhigh,vmpert,m,itheta,rtheta,i,lwork,vn
      REAL(r8) :: qa,kernelsignin,rplus,thetai
      CHARACTER(1), PARAMETER :: tab=CHAR(9)
      LOGICAL, PARAMETER :: complex_flag=.TRUE.,wall_flag=.FALSE.
      LOGICAL, PARAMETER :: farwal_flag=.TRUE.

      INTEGER, DIMENSION(:), POINTER :: ipiv,vmfac
      REAL(r8), DIMENSION(:), POINTER :: vtheta,vrfac,veta,vr,vz,
     $     dphi,delte,grri_real,grri_imag,grre_real,grre_imag,rwork
      COMPLEX(r8), DIMENSION(:), POINTER :: vbwp_mn,rbwp_mn,vbwp_fun,
     $     rbwp_fun,chi_fun,che_fun,flx_fun,kax_fun,work
      REAL(r8), DIMENSION(:,:), POINTER :: vgrri,vgrre,vxzpts
      COMPLEX(r8), DIMENSION(:,:), POINTER :: vflxmats,vkaxmats,
     $     temp1,temp2,vwv

c-----------------------------------------------------------------------
c     specify whatever boundary here with normal polar angles.
c-----------------------------------------------------------------------
      vmtheta=256
      vmlow=mlow
      vmhigh=mhigh
      vmpert=vmhigh-vmlow+1
      qlim=1.0
      vn=nn

      ALLOCATE(vmfac(vmpert))
      ALLOCATE(vtheta(0:vmtheta),vrfac(0:vmtheta),veta(0:vmtheta),
     $     vr(0:vmtheta),vz(0:vmtheta),
     $     dphi(0:vmtheta),delte(0:vmtheta))

      vmfac=(/(m,m=vmlow,vmhigh)/)      
      vtheta=(/(itheta,itheta=0,vmtheta)/)/REAL(vmtheta,r8)
      veta=twopi*vtheta
      rplus=0
      IF(present(dist))rplus=dist
      IF(present(majr).AND. present(minr)) THEN 
         vrfac=minr+rplus
         vr=majr+vrfac*COS(veta)
         vz=0.0+vrfac*SIN(veta)
      ENDIF

      IF(present(psi)) THEN
         DO itheta=0,vmtheta
            CALL bicube_eval(rzphi,psi,vtheta(itheta),0)
            veta(itheta)=twopi*(vtheta(itheta)+rzphi%f(2))
         ENDDO

         DO itheta=0,vmtheta
            thetai=issect(vmtheta,vtheta(:),veta(:),vtheta(itheta))
            CALL bicube_eval(rzphi,psi,thetai,0)
            vrfac(itheta)=SQRT(rzphi%f(1))+rplus
            vr(itheta)=ro+vrfac(itheta)*COS(veta(itheta))
            vz(itheta)=zo+vrfac(itheta)*SIN(veta(itheta))
            dphi(itheta)=rzphi%f(3)
         ENDDO
      ENDIF

c     IF(present(surfile)) THEN
c     ENDIF
      !surface jacobian for polar coordinates.
      dphi=0
      delte=-dphi/qlim
c-----------------------------------------------------------------------
c     invert values for vn < 0.
c-----------------------------------------------------------------------
      qa=qlim 
      IF(vn<0) THEN
         qa=-qa
         delte=-delte
         vn=-vn
      ENDIF  
c-----------------------------------------------------------------------
c     estimate prad for vacuum.
c-----------------------------------------------------------------------
      prad=(MAXVAL(vr)-MINVAL(vr))/2.0
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      CALL ascii_open(bin_unit,'ahg2msc_gpecarb.out',"UNKNOWN")
      WRITE(bin_unit,'(i4,a)')vmtheta,tab//tab//"mtheta"//tab//"mthin"
     $     //tab//"Number of poloidal nodes"
      WRITE(bin_unit,'(i4,a)')vmlow,tab//tab//"mlow"//tab//"lmin"//tab
     $     //"Lowest poloidal harmonic"
      WRITE(bin_unit,'(i4,a,a)')vmhigh,tab//tab//"mhigh"//tab//"lmax"
     $     //tab//"Highest poloidal harmonic"
      WRITE(bin_unit,'(i4,a)')vn,tab//tab//"nn"//tab//"nadj"//tab
     $     //"Toroidal harmonic"
      WRITE(bin_unit,'(f13.10,a)')qa,tab//"qa"//tab//"qa1"//tab
     $     //"Safety factor at plasma edge"
c-----------------------------------------------------------------------
c     write arrays.
c-----------------------------------------------------------------------
      WRITE(bin_unit,'(/a/)')"Poloidal Coordinate Theta:"
      WRITE(bin_unit,'(1p,4e18.10)')(1-vtheta(itheta),
     $     itheta=vmtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Polar Angle Eta:"
      WRITE(bin_unit,'(1p,4e18.10)')(twopi-veta(itheta),
     $     itheta=vmtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Radial Coordinate X:"
      WRITE(bin_unit,'(1p,4e18.10)')(vr(itheta),itheta=vmtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Axial Coordinate Z:"
      WRITE(bin_unit,'(1p,4e18.10)')(vz(itheta),itheta=vmtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Toroidal Angle Difference Delta:"
      WRITE(bin_unit,'(1p,4e18.10)')(delte(itheta),
     $     itheta=vmtheta,0,-1)
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     get grri and grre matrices by calling mscvac function.
c-----------------------------------------------------------------------
      ALLOCATE(vwv(vmpert,vmpert))
      nths=vmtheta+5
      nths2=nths*2
      nfm2=vmpert*2
      ALLOCATE(vgrri(nths2,nfm2),vgrre(nths2,nfm2),vxzpts(nths,4))
      kernelsignin=-1.0
      CALL mscvac(vwv,vmpert,vmtheta,vmtheta,complex_flag,
     $     kernelsignin,wall_flag,farwal_flag,vgrri,vxzpts)
      kernelsignin=1.0
      CALL mscvac(vwv,vmpert,vmtheta,vmtheta,complex_flag,
     $     kernelsignin,wall_flag,farwal_flag,vgrre,vxzpts)
      !vgrre(i+vmtheta,:) possibly provides coupling and mutuals!
c-----------------------------------------------------------------------
c     construct surface inductance matrix for specified boundary.
c-----------------------------------------------------------------------
      ALLOCATE(vbwp_mn(vmpert),rbwp_mn(vmpert))
      ALLOCATE(vbwp_fun(0:vmtheta),rbwp_fun(0:vmtheta),
     $     chi_fun(0:vmtheta),che_fun(0:vmtheta),flx_fun(0:vmtheta),
     $     kax_fun(0:vmtheta))
      ALLOCATE(grri_real(nths2),grri_imag(nths2),
     $     grre_real(nths2),grre_imag(nths2))
      ALLOCATE(vflxmats(vmpert,vmpert),vkaxmats(vmpert,vmpert))
      ALLOCATE(ipiv(vmpert),rwork(3*vmpert-2),work(2*vmpert-1))
      ALLOCATE(temp1(vmpert,vmpert),temp2(vmpert,vmpert))
      ALLOCATE(vsurf_indev(vmpert),vsurf_indmats(vmpert,vmpert))
      DO i=1,vmpert
         vbwp_mn=0
         vbwp_mn(i)=1.0
         CALL iscdftb(vmfac,vmpert,vbwp_fun,vmtheta,vbwp_mn)
         DO itheta=0,vmtheta
            rtheta=vmtheta-itheta
            rbwp_fun(itheta)=CONJG(vbwp_fun(rtheta))
         ENDDO
         CALL iscdftf(vmfac,vmpert,rbwp_fun,vmtheta,rbwp_mn)
         grri_real=MATMUL(vgrri,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grri_imag=MATMUL(vgrri,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
         grre_real=MATMUL(vgrre,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grre_imag=MATMUL(vgrre,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/)) 
                  
         DO itheta=0,vmtheta-1
            rtheta=vmtheta-itheta
            chi_fun(itheta+1)=(grri_real(rtheta)-
     $           ifac*grri_imag(rtheta))*EXP(-ifac*vn*dphi(itheta+1))
            che_fun(itheta+1)=(grre_real(rtheta)-
     $           ifac*grre_imag(rtheta))*EXP(-ifac*vn*dphi(itheta+1))
         ENDDO
         chi_fun(0)=chi_fun(vmtheta)
         che_fun(0)=che_fun(vmtheta)

         chi_fun=chi_fun/(twopi**2)
         che_fun=-che_fun/(twopi**2)
         IF (vmfac(i) == 0) THEN
            IF(verbose) WRITE(*,*) chi_fun(0)
            IF(verbose) WRITE(*,*) che_fun(0)
         ENDIF
         flx_fun=vbwp_fun
         kax_fun=(chi_fun-che_fun)/mu0

         CALL iscdftf(vmfac,vmpert,flx_fun,vmtheta,vflxmats(:,i))   
         CALL iscdftf(vmfac,vmpert,kax_fun,vmtheta,vkaxmats(:,i)) 
      ENDDO

      temp1=TRANSPOSE(vkaxmats)
      temp2=TRANSPOSE(vflxmats)
      CALL zgetrf(vmpert,vmpert,temp1,vmpert,ipiv,info)
      CALL zgetrs('N',vmpert,vmpert,temp1,vmpert,ipiv,temp2,vmpert,info)
      temp1=TRANSPOSE(temp2)
      temp1=0.5*(temp1+CONJG(TRANSPOSE(temp1)))
      vsurf_indmats=temp1
      lwork=2*vmpert-1
      CALL zheev('V','U',vmpert,temp1,vmpert,vsurf_indev,work,
     $     lwork,rwork,info)
      DEALLOCATE(vmfac,vtheta,vrfac,veta,vr,vz,dphi,delte)
      DEALLOCATE(vbwp_mn,rbwp_mn,vbwp_fun,rbwp_fun,chi_fun,che_fun,
     $     flx_fun,kax_fun,grri_real,grri_imag,grre_real,grre_imag,
     $     vflxmats,ipiv,rwork,work,temp1,temp2,
     $     vgrri,vgrre,vwv,vxzpts)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpvacuum_arbsurf
c-----------------------------------------------------------------------
c     subprogram 2. gpvacuum_flxsurf.
c     compute flux surface inductances.
c-----------------------------------------------------------------------
      SUBROUTINE gpvacuum_flxsurf(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER :: i,itheta,rtheta,lwork
      REAL(r8) :: kernelsignin
      CHARACTER(1), PARAMETER :: tab=CHAR(9)
      LOGICAL, PARAMETER :: complex_flag=.TRUE.,wall_flag=.FALSE.      
      LOGICAL, PARAMETER :: farwal_flag=.TRUE.

      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work

      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      
      COMPLEX(r8), DIMENSION(mpert) :: vbwp_mn,rbwp_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: chi_fun,che_fun,kax_fun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: fflxmats,fkaxmats,
     $     temp1,temp2,vwv

      REAL(r8), DIMENSION(:), POINTER :: grri_real,grri_imag,
     $     grre_real,grre_imag
      REAL(r8), DIMENSION(:,:), POINTER :: vgrri,vgrre,vxzpts
      CHARACTER(128) :: ahg_file
      CHARACTER(6) :: ahg_file_num

c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      WRITE(ahg_file,'(A16,F6.4,A4)') "ahg2msc_gpecflx_",psi,".out"

      CALL ahg_write(psi, ahg_file)
c-----------------------------------------------------------------------
c     get grri and grre matrices by calling mscvac functions.
c-----------------------------------------------------------------------
      nths=mthsurf+5
      nths2=nths*2
      nfm2=mpert*2
      ALLOCATE(vgrri(nths2,nfm2),vgrre(nths2,nfm2),vxzpts(nths,4))
      kernelsignin=-1.0
      CALL mscvac(vwv,mpert,mtheta,mthsurf,complex_flag,
     $     kernelsignin,wall_flag,farwal_flag,vgrri,vxzpts,ahg_file)
      kernelsignin=1.0
      CALL mscvac(vwv,mpert,mtheta,mthsurf,complex_flag,
     $     kernelsignin,wall_flag,farwal_flag,vgrre,vxzpts,ahg_file)
c-----------------------------------------------------------------------
c     construct surface inductance matrix for specified boundary.
c-----------------------------------------------------------------------
      ALLOCATE(grri_real(nths2),grri_imag(nths2),
     $     grre_real(nths2),grre_imag(nths2))

      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         dphi(itheta)=rzphi%f(3)
      ENDDO

      DO i=1,mpert
         vbwp_mn=0
         vbwp_mn(i)=1.0
         rbwp_mn=CONJG(vbwp_mn)
         grri_real=MATMUL(vgrri,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grri_imag=MATMUL(vgrri,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
         grre_real=MATMUL(vgrre,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grre_imag=MATMUL(vgrre,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))            
         
         DO itheta=0,mthsurf-1
            rtheta=mthsurf-itheta
            chi_fun(itheta+1)=(grri_real(rtheta)-
     $           ifac*grri_imag(rtheta))*EXP(-ifac*nn*dphi(itheta+1))
            che_fun(itheta+1)=(grre_real(rtheta)-
     $           ifac*grre_imag(rtheta))*EXP(-ifac*nn*dphi(itheta+1))
         ENDDO
         chi_fun(0)=chi_fun(mthsurf)
         che_fun(0)=che_fun(mthsurf)
         !jacobian with twopi angles.
         chi_fun=chi_fun/(twopi**2)
         che_fun=-che_fun/(twopi**2)         
         kax_fun=(chi_fun-che_fun)/mu0

         fflxmats(:,i)=vbwp_mn   
         CALL iscdftf(mfac,mpert,kax_fun,mthsurf,fkaxmats(:,i)) 
      ENDDO
      temp1=TRANSPOSE(fkaxmats)
      temp2=TRANSPOSE(fflxmats)
      CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
      CALL zgetrs('N',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
      temp1=TRANSPOSE(temp2)
      temp1=0.5*(temp1+CONJG(TRANSPOSE(temp1)))
      fsurf_indmats=temp1
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,temp1,mpert,fsurf_indev,work,
     $     lwork,rwork,info)

      DEALLOCATE(grri_real,grri_imag,grre_real,grre_imag,vgrri,vgrre,
     $     vxzpts)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpvacuum_flxsurf
c-----------------------------------------------------------------------
c     subprogram 3. gpvacuum_bnormal.
c     create bnormal input for vacuum code.  
c-----------------------------------------------------------------------
      SUBROUTINE gpvacuum_bnormal(psi,finmn,nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn

      INTEGER :: ipert,itheta,vn,vac_unit
      REAL(r8) :: qa,x1,x2,z1,z2
      CHARACTER(1), PARAMETER :: tab=CHAR(9)

      REAL(r8), DIMENSION(0:mthsurf) :: etas,dphi,delte

      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
c-----------------------------------------------------------------------
c     specify flux surface given by equilibrium file.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),0)
         rfac=SQRT(rzphi%f(1))
         etas(itheta)=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(etas(itheta))
         z(itheta)=zo+rfac*SIN(etas(itheta))
         dphi(itheta)=rzphi%f(3)
      ENDDO
      IF(ALLOCATED(gdr))THEN
         x1 = gdr(0,0)
         x2 = gdr(nr,0)
         z1 = gdz(0,0)
         z2 = gdz(0,nz)
      ELSE
         x1 = rmin
         x2 = rmax
         z1 = -zlim
         z2 = zlim      
      ENDIF

      delte=-dphi/sq%f(4)
c-----------------------------------------------------------------------
c     invert values for vn < 0.
c-----------------------------------------------------------------------
      vn=nn
      qa=sq%f(4)
      IF(nn <0) THEN
         qa=-qa
         delte=-delte
         vn=-vn
      ENDIF
c-----------------------------------------------------------------------
c     reverse poloidal and torodial coordinates and normalize.
c-----------------------------------------------------------------------
      rbwp_mn=CONJG(finmn)/(twopi**2)
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      vac_unit=7
      CALL ascii_open(vac_unit,'vacin5',"UNKNOWN")
      WRITE(vac_unit,'(a/)')"scalars"
      WRITE(vac_unit,'(i4,a)')nr+1,tab//tab//"Number of x grid"
      WRITE(vac_unit,'(i4,a)')nz+1,tab//tab//"Number of z grid"
      WRITE(vac_unit,'(f18.10,a)')x1,tab//tab//
     $     "Left x position in rectangular grid"
      WRITE(vac_unit,'(f18.10,a)')x2,tab//tab//
     $     "Right x position in rectangular grid"
      WRITE(vac_unit,'(f18.10,a)')z1,tab//tab//
     $     "Lower z position in rectangular grid"
      WRITE(vac_unit,'(f18.10,a)')z2,tab//tab//
     $     "Upper z position in rectangular grid"
      WRITE(vac_unit,'(i4,a)')mthsurf,tab//tab//"mthsurf"//tab//"mthin"
     $     //tab//"Number of poloidal nodes"
      WRITE(vac_unit,'(i4,a)')mlow,tab//tab//"mlow"//tab//"lmin"//tab
     $     //"Lowest poloidal harmonic"
      WRITE(vac_unit,'(i4,a,a)')mhigh,tab//tab//"mhigh"//tab//"lmax"
     $     //tab//"Highest poloidal harmonic"
      WRITE(vac_unit,'(i4,a)')vn,tab//tab//"nn"//tab//tab//"ntor"//tab
     $     //"Toroidal harmonic"
      WRITE(vac_unit,'(f13.10,a)')qa,tab//"qa"//tab//"qa1"//tab
     $     //"Safety factor at plasma edge"
c-----------------------------------------------------------------------
c     write arrays.
c-----------------------------------------------------------------------
      WRITE(vac_unit,'(a/)')"Poloidal Coordinate Theta:"
      WRITE(vac_unit,'(1p,4e18.10)')(1-theta(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(vac_unit,'(a/)')"Polar Angle Eta:"
      WRITE(vac_unit,'(1p,4e18.10)')(twopi-etas(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(vac_unit,'(a/)')"Radial Coordinate X:"
      WRITE(vac_unit,'(1p,4e18.10)')(r(itheta),itheta=mthsurf,0,-1)
      WRITE(vac_unit,'(a/)')"Axial Coordinate Z:"
      WRITE(vac_unit,'(1p,4e18.10)')(z(itheta),itheta=mthsurf,0,-1)
      WRITE(vac_unit,'(a/)')"Toroidal Angle Difference Delta:"
      WRITE(vac_unit,'(1p,4e18.10)')(delte(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(vac_unit,'(a/)')"Real component of normal b field:"
      WRITE(vac_unit,'(1p,4e18.10)')(REAL(rbwp_mn(ipert)),
     $     ipert=1,mpert)
      WRITE(vac_unit,'(a/)')"Imaginary component of normal b field:"
      WRITE(vac_unit,'(1p,4e18.10)')(AIMAG(rbwp_mn(ipert)),
     $     ipert=1,mpert)
      CALL ascii_close(vac_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpvacuum_bnormal
c-----------------------------------------------------------------------
c     subprogram 4. gpvacuum_ideal_mutuals.
c     calculate ideal mutual inductances.  
c-----------------------------------------------------------------------
      SUBROUTINE gpvacuum_ideal_mutuals
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8) :: majr,minr,psi,dist,err1,err2,rval
      INTEGER :: vmpert,i,info,lwork
      INTEGER, DIMENSION(:), POINTER :: ipiv
      REAL(r8), DIMENSION(:), POINTER :: d1,rwork
      COMPLEX(r8), DIMENSION(:), POINTER :: d2,work
      COMPLEX(r8), DIMENSION(:,:),POINTER :: lmat1,lmat2,lmat12,lmat21,
     $     ilmat1,ilmat2,mmat,immat,dmat,temp1,temp2,vl,vr

      CALL gpvacuum_arbsurf(psi=psilim)
      vmpert=SIZE(vsurf_indev)
      ALLOCATE(lmat1(vmpert,vmpert))
      lmat1=vsurf_indmats
      DEALLOCATE(vsurf_indmats,vsurf_indev)

      dist=prad*0.25
      CALL gpvacuum_arbsurf(psi=psilim,dist=dist)
      vmpert=SIZE(vsurf_indev)
      ALLOCATE(lmat2(vmpert,vmpert))
      lmat2=vsurf_indmats
      DEALLOCATE(vsurf_indmats,vsurf_indev)

      ALLOCATE(lmat12(vmpert,vmpert),lmat21(vmpert,vmpert),
     $     ilmat1(vmpert,vmpert),ilmat2(vmpert,vmpert),
     $     mmat(vmpert,vmpert),immat(vmpert,vmpert),
     $     temp1(vmpert,vmpert),temp2(vmpert,vmpert))
      lmat12=MATMUL(lmat1,lmat2)
      lmat21=MATMUL(lmat2,lmat1)
      err1=MAXVAL(ABS(lmat2-TRANSPOSE(CONJG(lmat2))))/MAXVAL(ABS(lmat2))
      err2=MAXVAL(ABS(lmat12-lmat21))/MAXVAL(ABS(lmat12))
c-----------------------------------------------------------------------
c     invert self inductance matrices.
c-----------------------------------------------------------------------
      temp1=lmat1
      CALL zpotrf('U',vmpert,temp1,vmpert,info)
      CALL zpotri('U',vmpert,temp1,vmpert,info)
      ilmat1=temp1
      temp1=lmat2
      CALL zpotrf('U',vmpert,temp1,vmpert,info)
      CALL zpotri('U',vmpert,temp1,vmpert,info)
      ilmat2=temp1
c-----------------------------------------------------------------------
c     calculate mutual inductance matrices.
c-----------------------------------------------------------------------
      mmat=MATMUL(lmat2,ilmat1)
      lwork=2*vmpert+1
      ALLOCATE(ipiv(vmpert))
      ALLOCATE(d2(vmpert),rwork(2*vmpert),work(lwork))
      ALLOCATE(dmat(vmpert,vmpert),vl(vmpert,vmpert),vr(vmpert,vmpert))
      temp1=mmat
      CALL zgeev('V','V',vmpert,temp1,vmpert,d2,
     $        vl,vmpert,vr,vmpert,work,lwork,rwork,info)
      dmat=0
      ! approximate for dmat^2=(a+bi)^2=d2 for a>0, a>>b
      DO i=1,vmpert
                  rval=SQRT(REAL(d2(i)))
         dmat(i,i)=rval+ifac*AIMAG(d2(i))/(2.0*rval)
      ENDDO
      DEALLOCATE(work)
      lwork=vmpert
      ALLOCATE(work(lwork))
      temp1=vr
      ipiv=0
      CALL zgetrf(vmpert,vmpert,temp1,vmpert,ipiv,info)
      CALL zgetri(vmpert,temp1,vmpert,ipiv,work,lwork,info)      
      mmat=MATMUL(vr,MATMUL(dmat,temp1))
      mmat=MATMUL(mmat,lmat1)
      ! Here we complete M=(L2#L1^(-1))^(1/2)#L_1 in a matrix form.

      DEALLOCATE(work)
      lwork=2*vmpert+1
      ALLOCATE(work(lwork))
      immat=MATMUL(lmat1,ilmat2)
      temp1=immat
      CALL zgeev('V','V',vmpert,temp1,vmpert,d2,
     $        vl,vmpert,vr,vmpert,work,lwork,rwork,info)
      dmat=0
      DO i=1,vmpert
         WRITE(*,*)d2(i)
         rval=SQRT(REAL(d2(i)))
         dmat(i,i)=rval+ifac*AIMAG(d2(i))/(2.0*rval)
      ENDDO
      DEALLOCATE(work)
      lwork=vmpert
      ALLOCATE(work(lwork))
      temp1=vr
      ipiv=0
      CALL zgetrf(vmpert,vmpert,temp1,vmpert,ipiv,info)
      CALL zgetri(vmpert,temp1,vmpert,ipiv,work,lwork,info)     
      immat=MATMUL(vr,MATMUL(dmat,temp1))
      immat=MATMUL(immat,lmat2)
      ! Here we complete M=(L1#L2^(-1))^(1/2)#L_2 in a matrix form.

      err1=MAXVAL(ABS(mmat-immat))/MAXVAL(ABS(mmat))
      WRITE(*,*)"mutual_err=",err1
          
      DEALLOCATE(lmat1,lmat2,lmat12,lmat21)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpvacuum_ideal_mutuals

      END MODULE gpvacuum_mod
      

