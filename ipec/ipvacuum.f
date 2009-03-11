c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     IPVACUUM: manipulate vacuum code separately.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipvacuum
c      1. ipvacuum_arbsurf
c      2. ipvacuum_flxsurf
c      3. ipvacuum_bnormal
c-----------------------------------------------------------------------
c     subprogram 0. ipvacuum_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipvacuum_mod
      USE ipeq_mod

      IMPLICIT NONE
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipvacuum_arbsurf.
c     compute arbitrary surface inductance.
c-----------------------------------------------------------------------
      SUBROUTINE ipvacuum_arbsurf(majr,minr)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: majr,minr

      INTEGER :: vmtheta,vmlow,vmhigh,vmpert,m,itheta,rtheta,i,vn,lwork
      REAL(r8) :: qa,kernelsignin
      CHARACTER(1), PARAMETER :: tab=CHAR(9)
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(:), POINTER :: ipiv,vmfac
      REAL(r8), DIMENSION(:), POINTER :: vtheta,vrfac,veta,vr,vz,
     $     jac,dphi,delte,delpsi,wgtfun,grri_real,grri_imag,
     $     grre_real,grre_imag,rwork
      COMPLEX(r8), DIMENSION(:), POINTER :: vbwp_mn,rbwp_mn,vbwp_fun,
     $     rbwp_fun,chi_fun,che_fun,flx_fun,kax_fun,work
      REAL(r8), DIMENSION(:,:), POINTER :: vgrri,vgrre
      COMPLEX(r8), DIMENSION(:,:), POINTER :: vflxmats,vkaxmats,
     $     temp1,temp2,vwv

c-----------------------------------------------------------------------
c     specify whatever boundary here with normal polar angles.
c-----------------------------------------------------------------------
      vmtheta=256
      vmlow=-10
      vmhigh=10
      vmpert=vmhigh-vmlow+1
      nn=1
      qlim=1.0

      ALLOCATE(vmfac(vmpert))
      ALLOCATE(vtheta(0:vmtheta),vrfac(0:vmtheta),veta(0:vmtheta),
     $     vr(0:vmtheta),vz(0:vmtheta),jac(0:vmtheta),
     $     dphi(0:vmtheta),delte(0:vmtheta),
     $     delpsi(0:vmtheta),wgtfun(0:vmtheta))

      vmfac=(/(m,m=vmlow,vmhigh)/)      
      vtheta=(/(itheta,itheta=0,vmtheta)/)/REAL(vmtheta,r8)
      vrfac=minr
      veta=twopi*vtheta
      vr=majr+vrfac*COS(veta)
      vz=0.0+vrfac*SIN(veta)
      jac=vrfac*vr
      dphi=0
      delpsi=1.0
      wgtfun=1.0/(jac*delpsi)
      delte=-dphi/qlim
c-----------------------------------------------------------------------
c     invert values for vn < 0.
c-----------------------------------------------------------------------
      vn=nn
      qa=qlim 
      IF(vn <0) THEN
         qa=-qa
         delte=-delte
         vn=-vn
      ENDIF  
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      CALL ascii_open(bin_unit,'ahg2msc.out',"UNKNOWN")
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
      WRITE(bin_unit,'(/a/)')"Toroidal Angle Difference Delte:"
      WRITE(bin_unit,'(1p,4e18.10)')(delte(itheta),
     $     itheta=vmtheta,0,-1)
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     get grri,grre matrices by calling mscvac function.
c-----------------------------------------------------------------------
      ALLOCATE(vwv(vmpert,vmpert))
      kernelsignin=-1.0
      CALL mscvac(vwv,vmpert,vmtheta,vmtheta,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(vgrri(nths2,nfm2))
      CALL grrget(nfm2,nths2,vgrri)
      kernelsignin=1.0
      CALL mscvac(vwv,vmpert,vmtheta,vmtheta,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(vgrre(nths2,nfm2))
      CALL grrget(nfm2,nths2,vgrre)
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
         vbwp_fun=vbwp_fun*jac
         DO itheta=0,vmtheta-1
            rtheta=vmtheta-itheta
            rbwp_fun(itheta+1)=CONJG(vbwp_fun(rtheta))
         ENDDO
         rbwp_fun(0)=rbwp_fun(vmtheta)
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
            WRITE(*,*) chi_fun(0)
            WRITE(*,*) che_fun(0)
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
      WRITE(*,*)vsurf_indev(1),vsurf_indev(2),
     $     vsurf_indev(2)/vsurf_indev(1)
      DEALLOCATE(vmfac,vtheta,vrfac,veta,vr,vz,jac,dphi,delte,
     $     delpsi,wgtfun)
      DEALLOCATE(vbwp_mn,rbwp_mn,vbwp_fun,rbwp_fun,chi_fun,che_fun,
     $     flx_fun,kax_fun,grri_real,grri_imag,grre_real,grre_imag,
     $     vflxmats,vkaxmats,ipiv,rwork,work,temp1,temp2,
     $     vgrri,vgrre,vwv)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipvacuum_arbsurf
c-----------------------------------------------------------------------
c     subprogram 2. ipvacuum_flxsurf.
c     compute flux surface inductances.
c     need to use ipeq_alloc and dealloc routines externally. 
c-----------------------------------------------------------------------
      SUBROUTINE ipvacuum_flxsurf(psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      INTEGER :: i,itheta,rtheta,vn,lwork
      REAL(r8) :: jac,qa,kernelsignin
      CHARACTER(1), PARAMETER :: tab=CHAR(9)
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work

      REAL(r8), DIMENSION(0:mtheta) :: vtheta,vrfac,veta,vr,vz,
     $     vdphi,delte
      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,wgtfun
      
      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,rbwp_fun,chi_fun,
     $     che_fun,flx_fun,kax_fun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: fflxmats,fkaxmats,
     $     temp1,temp2,vwv

      REAL(r8), DIMENSION(:), POINTER :: grri_real,grri_imag,
     $     grre_real,grre_imag
      REAL(r8), DIMENSION(:,:), POINTER :: vgrri,vgrre
c-----------------------------------------------------------------------
c     specify flux surface in hamada given by equilibrium file.
c-----------------------------------------------------------------------
      vtheta=rzphi%ys
      CALL spline_eval(sq,psi,0)
      DO itheta=0,mtheta
         CALL bicube_eval(rzphi,psi,vtheta(itheta),0)
         vrfac(itheta)=SQRT(rzphi%f(1))
         veta(itheta)=twopi*(vtheta(itheta)+rzphi%f(2))
         vr(itheta)=ro+vrfac(itheta)*COS(veta(itheta))
         vz(itheta)=zo+vrfac(itheta)*SIN(veta(itheta))
         vdphi(itheta)=rzphi%f(3)
      ENDDO
      delte=-vdphi/sq%f(4)
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
c     write scalars.
c-----------------------------------------------------------------------
      CALL ascii_open(bin_unit,'ahg2msc.out',"UNKNOWN")
      WRITE(bin_unit,'(i4,a)')mtheta,tab//tab//"mtheta"//tab//"mthin"
     $     //tab//"Number of poloidal nodes"
      WRITE(bin_unit,'(i4,a)')mlow,tab//tab//"mlow"//tab//"lmin"//tab
     $     //"Lowest poloidal harmonic"
      WRITE(bin_unit,'(i4,a,a)')mhigh,tab//tab//"mhigh"//tab//"lmax"
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
     $     itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Polar Angle Eta:"
      WRITE(bin_unit,'(1p,4e18.10)')(twopi-veta(itheta),
     $     itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Radial Coordinate X:"
      WRITE(bin_unit,'(1p,4e18.10)')(vr(itheta),itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Axial Coordinate Z:"
      WRITE(bin_unit,'(1p,4e18.10)')(vz(itheta),itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Toroidal Angle Difference Delte:"
      WRITE(bin_unit,'(1p,4e18.10)')(delte(itheta),
     $     itheta=mtheta,0,-1)
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     get grri,grre matrices from call mscvac functions.
c-----------------------------------------------------------------------
      kernelsignin=-1.0
      CALL mscvac(vwv,mpert,mtheta,mthvac,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(vgrri(nths2,nfm2))
      CALL grrget(nfm2,nths2,vgrri)
      kernelsignin=1.0
      CALL mscvac(vwv,mpert,mtheta,mthvac,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(vgrre(nths2,nfm2))
      CALL grrget(nfm2,nths2,vgrre)
c-----------------------------------------------------------------------
c     construct surface inductance matrix for specified boundary.
c     what is the last point? check!!
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
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         wgtfun(itheta)=1.0/(jac*delpsi(itheta))
      ENDDO

      DO i=1,mpert
         bwp_mn=0
         bwp_mn(i)=1.0
         CALL iscdftb(mfac,mpert,bwp_fun,mthvac,bwp_mn)
         DO itheta=0,mthvac-1
            rtheta=mthvac-itheta
            rbwp_fun(itheta+1)=CONJG(bwp_fun(rtheta))
         ENDDO
         rbwp_fun(0)=rbwp_fun(mthvac)
         CALL iscdftf(mfac,mpert,rbwp_fun,mthvac,rbwp_mn)
         grri_real=MATMUL(vgrri,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grri_imag=MATMUL(vgrri,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))
         grre_real=MATMUL(vgrre,(/REAL(rbwp_mn),-AIMAG(rbwp_mn)/))
         grre_imag=MATMUL(vgrre,(/AIMAG(rbwp_mn),REAL(rbwp_mn)/))            
         
         DO itheta=0,mthvac-1
            rtheta=mthvac-itheta
            chi_fun(itheta+1)=(grri_real(rtheta)-
     $           ifac*grri_imag(rtheta))*EXP(-ifac*nn*dphi(itheta+1))
            che_fun(itheta+1)=(grre_real(rtheta)-
     $           ifac*grre_imag(rtheta))*EXP(-ifac*nn*dphi(itheta+1))
         ENDDO
         chi_fun(0)=chi_fun(mthvac)
         che_fun(0)=che_fun(mthvac)
         chi_fun=chi_fun*jac/(twopi**2)
         che_fun=-che_fun*jac/(twopi**2)         
         flx_fun=bwp_fun*jac
         kax_fun=(chi_fun-che_fun)/mu0
         CALL iscdftf(mfac,mpert,flx_fun,mthsurf,fflxmats(:,i))   
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

      DEALLOCATE(grri_real,grri_imag,grre_real,grre_imag,vgrri,vgrre)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipvacuum_flxsurf
c-----------------------------------------------------------------------
c     subprogram 3. ipvacuum_bnormal.
c     create bnormal input for vacuum code.  
c-----------------------------------------------------------------------
      SUBROUTINE ipvacuum_bnormal(psi,bnomn,polo,toro,nr,nz,acpl,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: polo,toro,nr,nz,labl
      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnomn
      CHARACTER(6) :: acpl
      CHARACTER(1) :: slabl

      INTEGER :: ipert,itheta,rtheta,vn
      REAL(r8) :: qa
      CHARACTER(1), PARAMETER :: tab=CHAR(9)

      REAL(r8), DIMENSION(0:mthsurf) :: etas,dphi,delpsi,delte,jacs

      COMPLEX(r8), DIMENSION(mpert) :: rbwp_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,rbwp_fun
c-----------------------------------------------------------------------
c     specify flux surface in hamada given by equilibrium file.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         etas(itheta)=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(etas(itheta))
         z(itheta)=zo+rfac*SIN(etas(itheta))
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         dphi(itheta)=rzphi%f(3)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         jacs(itheta)=jac
      ENDDO
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
c     change and reverse poloidal and torodial coordinates in hamada.
c-----------------------------------------------------------------------
      IF ((polo /= 1).OR.(toro /=1)) THEN
         CALL ipeq_cotoha(psi,bnomn,mfac,mpert,polo,toro)
      ENDIF
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bnomn)
      bwp_fun=bwp_fun*delpsi*jacs/(twopi**2)
      DO itheta=0,mthvac-1
         rtheta=mthvac-itheta
         rbwp_fun(itheta+1)=CONJG(bwp_fun(rtheta))
      ENDDO
      rbwp_fun(0)=rbwp_fun(mthvac)
      CALL iscdftf(mfac,mpert,rbwp_fun,mthvac,rbwp_mn)
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl,FMT='(I1)')labl
      CALL ascii_open(bin_unit,'vacin_'//acpl//'_l'//slabl,"UNKNOWN")
      WRITE(bin_unit,'(a/)')"scalars"
      WRITE(bin_unit,'(i4,a)')nr+1,tab//tab//"Number of x grid"
      WRITE(bin_unit,'(i4,a)')nz+1,tab//tab//"Number of z grid"
      WRITE(bin_unit,'(f18.10,a)')psi_in%xs(0),tab//tab//
     $     "Left x position in rectangular grid"
      WRITE(bin_unit,'(f18.10,a)')psi_in%xs(mr),tab//tab//
     $     "Right x position in rectangular grid"
      WRITE(bin_unit,'(f18.10,a)')psi_in%ys(0),tab//tab//
     $     "Lower z position in rectangular grid"
      WRITE(bin_unit,'(f18.10,a)')psi_in%ys(mz),tab//tab//
     $     "Upper z position in rectangular grid"
      WRITE(bin_unit,'(i4,a)')mthsurf,tab//tab//"mthsurf"//tab//"mthin"
     $     //tab//"Number of poloidal nodes"
      WRITE(bin_unit,'(i4,a)')mlow,tab//tab//"mlow"//tab//"lmin"//tab
     $     //"Lowest poloidal harmonic"
      WRITE(bin_unit,'(i4,a,a)')mhigh,tab//tab//"mhigh"//tab//"lmax"
     $     //tab//"Highest poloidal harmonic"
      WRITE(bin_unit,'(i4,a)')vn,tab//tab//"nn"//tab//tab//"ntor"//tab
     $     //"Toroidal harmonic"
      WRITE(bin_unit,'(f13.10,a)')qa,tab//"qa"//tab//"qa1"//tab
     $     //"Safety factor at plasma edge"
c-----------------------------------------------------------------------
c     write arrays.
c-----------------------------------------------------------------------
      WRITE(bin_unit,'(a/)')"Poloidal Coordinate Theta:"
      WRITE(bin_unit,'(1p,4e18.10)')(1-theta(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(bin_unit,'(a/)')"Polar Angle Eta:"
      WRITE(bin_unit,'(1p,4e18.10)')(twopi-etas(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(bin_unit,'(a/)')"Radial Coordinate X:"
      WRITE(bin_unit,'(1p,4e18.10)')(r(itheta),itheta=mthsurf,0,-1)
      WRITE(bin_unit,'(a/)')"Axial Coordinate Z:"
      WRITE(bin_unit,'(1p,4e18.10)')(z(itheta),itheta=mthsurf,0,-1)
      WRITE(bin_unit,'(a/)')"Toroidal Angle Difference Delte:"
      WRITE(bin_unit,'(1p,4e18.10)')(delte(itheta),
     $     itheta=mthsurf,0,-1)
      WRITE(bin_unit,'(a/)')"Real component of normal b field:"
      WRITE(bin_unit,'(1p,4e18.10)')(REAL(rbwp_mn(ipert)),
     $     ipert=1,mpert)
      WRITE(bin_unit,'(a/)')"Imaginary component of normal b field:"
      WRITE(bin_unit,'(1p,4e18.10)')(AIMAG(rbwp_mn(ipert)),
     $     ipert=1,mpert)
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipvacuum_bnormal

      END MODULE ipvacuum_mod
      
