c-----------------------------------------------------------------------
c     file free.f.
c     stability of free-boundary modes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. free_mod.
c     1. free_run.
c     2. free_write_msc.
c     3. free_ahb_prep.
c     4. free_ahb_write.
c     5. free_get_wvac
c-----------------------------------------------------------------------
c     subprogram 0. free_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE free_mod
      USE ode_mod
      USE dcon_netcdf_mod
      IMPLICIT NONE

      INTEGER :: msol_ahb=-1
      INTEGER, PRIVATE :: mthsurf
      REAL(r8), PRIVATE :: qsurf,q1surf
      REAL(r8), DIMENSION(:), POINTER, PRIVATE :: theta,dphi,r,z
      REAL(r8), DIMENSION(:,:), POINTER, PRIVATE :: thetas
      REAL(r8), DIMENSION(:,:,:), POINTER, PRIVATE :: project

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: wvac
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: galwt

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. free_run.
c     computes plasma, vacuum, and total potential energies.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE free_run(plasma1,vacuum1,total1,nzero)

      REAL(r8), INTENT(OUT) :: plasma1,vacuum1,total1
      INTEGER, INTENT(IN) :: nzero

      LOGICAL, PARAMETER :: normalize=.TRUE.
      CHARACTER(1), DIMENSION(mpert,msol) :: star
      INTEGER :: ipert,jpert,isol,info,lwork
      INTEGER, DIMENSION(mpert) :: ipiv,m,eindex
      INTEGER, DIMENSION(1,msol) :: imax
      REAL(r8) :: v1
      REAL(r8), DIMENSION(mpert) :: singfac
      ! eigenvalues are complex for gpec
      COMPLEX(r8), DIMENSION(mpert) :: ep,ev,et,tt

      REAL(r8), DIMENSION(3*mpert-1) :: rwork
      REAL(r8), DIMENSION(2*mpert) :: rwork2
      COMPLEX(r8) :: phase,norm
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work2
      COMPLEX(r8), DIMENSION(mpert,mpert) :: wp,wv,wt,wt0,temp,wpt,wvt
      COMPLEX(r8), DIMENSION(mpert,mpert) :: nmat,smat
      COMPLEX(r8), DIMENSION(mpert,mpert) :: vl,vr
      CHARACTER(24), DIMENSION(mpert) :: message
      LOGICAL, PARAMETER :: complex_flag=.TRUE.,wall_flag=.FALSE.
      ! RDCON mutual inductance calculations not implemented
      LOGICAL :: wv_farwall_flag=.FALSE.
      LOGICAL :: farwal_flag=.FALSE.

      REAL(r8) :: kernelsignin
      INTEGER :: vac_unit
      REAL(r8), DIMENSION(:,:), POINTER :: grri,xzpts
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"Energies: plasma = ",es10.3,", vacuum = ",es10.3,
     $     ", real = ",es10.3,", imaginary = ",es10.3)
 20   FORMAT(/3x,"isol",3x,"plasma",5x,"vacuum",2x,"re total",2x,
     $     "im total"/)
 30   FORMAT(i6,1p,4e11.3,a)
 40   FORMAT(/3x,"isol",2x,"imax",3x,"plasma",5x,"vacuum",2x,"re total",
     $     2x,"im total"/)
 50   FORMAT(2i6,1p,4e11.3,a)
 60   FORMAT(/2x,"ipert",4x,"m",4x,"re wt",6x,"im wt",6x,"abs wt"/)
 70   FORMAT(2i6,1p,3e11.3,2x,a)
 80   FORMAT(/3x,"isol",3x,"plasma",5x,"vacuum"/)
 90   FORMAT(i6,1p,2e11.3)
c-----------------------------------------------------------------------
c     compute plasma response matrix.
c-----------------------------------------------------------------------
      IF(ode_flag)THEN
         temp=CONJG(TRANSPOSE(u(:,1:mpert,1)))
         wp=u(:,1:mpert,2)
         wp=CONJG(TRANSPOSE(wp))
         CALL zgetrf(mpert,mpert,temp,mpert,ipiv,info)  
         CALL zgetrs('N',mpert,mpert,temp,mpert,ipiv,wp,mpert,info)
         wp=(wp+CONJG(TRANSPOSE(wp)))/(2*psio**2)
      ELSE
         wp=0
      ENDIF
c-----------------------------------------------------------------------
c     calculate vacuum response matrix.
c-----------------------------------------------------------------------
      CALL free_write_msc
      IF(ahb_flag)CALL free_ahb_prep(wp,nmat,smat)
      CALL spline_eval(sq,psilim,0)
      v1=sq%f(3)
c-----------------------------------------------------------------------
c     compute vacuum response matrix.
c-----------------------------------------------------------------------
      vac_unit=4
      farwal_flag=.TRUE. ! self-inductance for plasma boundary.
      kernelsignin=-1.0
      ALLOCATE(grri(2*(mthvac+5),mpert*2),xzpts(mthvac+5,4))
      CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $     wall_flag,farwal_flag,grri,xzpts)
      CALL bin_open(vac_unit,"vacuum.bin","UNKNOWN","REWIND","none")
      WRITE(vac_unit)grri

      kernelsignin=1.0
      CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $     wall_flag,farwal_flag,grri,xzpts)
      WRITE(vac_unit)grri

      IF(wv_farwall_flag)THEN
         temp=wv
      ENDIF         

      farwal_flag=.FALSE. ! self-inductance with the wall.
      kernelsignin=-1.0
      CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $     wall_flag,farwal_flag,grri,xzpts)
      WRITE(vac_unit)grri

      kernelsignin=1.0
      CALL mscvac(wv,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $     wall_flag,farwal_flag,grri,xzpts)
      WRITE(vac_unit)grri
      WRITE(vac_unit)xzpts

      CALL bin_close(vac_unit)
      DEALLOCATE(grri,xzpts)
      singfac=mlow-nn*qlim+(/(ipert,ipert=0,mpert-1)/)
      DO ipert=1,mpert
         wv(ipert,:)=wv(ipert,:)*singfac
         wv(:,ipert)=wv(:,ipert)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     compute complex energy eigenvalues.
c-----------------------------------------------------------------------
      wt=wp+wv
      wt0=wt
      lwork=2*mpert+1
      CALL zgeev('V','V',mpert,wt,mpert,et,
     $        vl,mpert,vr,mpert,work2,lwork,rwork2,info)
      eindex(1:mpert)=(/(ipert,ipert=1,mpert)/)
      CALL bubble(REAL(et),eindex,1,mpert)
      tt=et
      DO ipert=1,mpert
         wt(:,ipert)=vr(:,eindex(mpert+1-ipert))
         et(ipert)=tt(eindex(mpert+1-ipert))
      ENDDO
c-----------------------------------------------------------------------
c     normalize eigenfunction and energy.
c-----------------------------------------------------------------------
      IF(normalize)THEN
         DO isol=1,mpert
            norm=0
            DO ipert=1,mpert
               DO jpert=1,mpert
                  norm=norm+jmat(jpert-ipert)
     $                 *wt(ipert,isol)*CONJG(wt(jpert,isol))
               ENDDO
            ENDDO
            norm=norm/v1
            wt(:,isol)=wt(:,isol)/SQRT(norm)
            et(isol)=et(isol)/norm
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     normalize phase and label largest component.
c-----------------------------------------------------------------------
      DO isol=1,mpert
         imax(:,isol)=MAXLOC(ABS(wt(:,isol)))
         phase=ABS(wt(imax(1,isol),isol))/wt(imax(1,isol),isol)
         wt(:,isol)=wt(:,isol)*phase
         star(:,isol)=' '
         star(imax(1,isol),isol)='*'
      ENDDO      
c-----------------------------------------------------------------------
c     compute plasma and vacuum contributions.
c-----------------------------------------------------------------------
      wpt=MATMUL(CONJG(TRANSPOSE(wt)),MATMUL(wp,wt))
      wvt=MATMUL(CONJG(TRANSPOSE(wt)),MATMUL(wv,wt))
      DO ipert=1,mpert
         ep(ipert)=wpt(ipert,ipert)
         ev(ipert)=wvt(ipert,ipert)
      ENDDO
      plasma1=REAL(ep(1))
      vacuum1=REAL(ev(1))
      total1=REAL(et(1))
c-----------------------------------------------------------------------
c     write data for ahb and deallocate.
c-----------------------------------------------------------------------
      IF(ahb_flag)THEN
         CALL free_ahb_write(nmat,smat,wt,et)
         DEALLOCATE(r,z,theta,dphi,thetas,project)
      ENDIF
c-----------------------------------------------------------------------
c     save eigenvalues and eigenvectors to file.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         WRITE(euler_bin_unit)3
         WRITE(euler_bin_unit)ep
         WRITE(euler_bin_unit)et
         WRITE(euler_bin_unit)wt
         WRITE(euler_bin_unit)wt0
         WRITE(euler_bin_unit)wv_farwall_flag
      ENDIF
c-----------------------------------------------------------------------
c     write to screen and copy to output.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,10) REAL(ep(1)),REAL(ev(1)),
     $     REAL(et(1)),AIMAG(et(1))
c-----------------------------------------------------------------------
c     write eigenvalues to file.
c-----------------------------------------------------------------------
      message=""
      WRITE(out_unit,'(/1x,a)')"Total Energy Eigenvalues:"
      WRITE(out_unit,20)
      WRITE(out_unit,30)(isol,REAL(ep(isol)),REAL(ev(isol)),
     $     REAL(et(isol)),AIMAG(et(isol)),
     $     TRIM(message(isol)),isol=1,mpert)
      WRITE(out_unit,20)
c-----------------------------------------------------------------------
c     write eigenvectors to file.
c-----------------------------------------------------------------------
      WRITE(out_unit,*)"Total Energy Eigenvectors:"
      m=mlow+(/(isol,isol=0,mpert-1)/)
      DO isol=1,mpert
         WRITE(out_unit,40)
         WRITE(out_unit,50)isol,imax(1,isol),REAL(ep(isol)),
     $         REAL(ev(isol)),REAL(et(isol)),
     $         AIMAG(et(isol)),TRIM(message(isol))
         WRITE(out_unit,60)
         WRITE(out_unit,70)(ipert,m(ipert),wt(ipert,isol),
     $        ABS(wt(ipert,isol)),star(ipert,isol),ipert=1,mpert)
         WRITE(out_unit,60)
      ENDDO
c-----------------------------------------------------------------------
c     write the plasma matrix.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/1x,a/)')"Plasma Energy Matrix:"
      DO isol=1,mpert
         WRITE(out_unit,'(1x,2(a,i3))')"isol = ",isol,", m = ",m(isol)
         WRITE(out_unit,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
         WRITE(out_unit,'(i3,1p,3e13.5)')
     $        (ipert,wp(ipert,isol),ABS(wp(ipert,isol)),ipert=1,mpert)
         WRITE(out_unit,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
      ENDDO
c-----------------------------------------------------------------------
c     compute separate plasma and vacuum eigenvalues.
c-----------------------------------------------------------------------
      ! CALL zheev('V','U',mpert,wp,mpert,ep,work,lwork,rwork,info)
      ! CALL zheev('V','U',mpert,wv,mpert,ev,work,lwork,rwork,info)
      ! WRITE(out_unit,*)"Separate Energy Eigenvalues:"
      ! WRITE(out_unit,80)
      ! WRITE(out_unit,90)(isol,ep(isol),ev(isol),isol=1,mpert)
      ! WRITE(out_unit,80)

      ! IF (.NOT.ALLOCATED(galwt)) ALLOCATE(galwt(mpert,mpert))
      ! galwt=wt*psio*twopi*1e-3
      







      lwork=2*mpert+1      
      CALL zgeev('V','V',mpert,wp,mpert,ep,
     $     vl,mpert,vr,mpert,work2,lwork,rwork2,info)
      eindex(1:mpert)=(/(ipert,ipert=1,mpert)/)
      CALL bubble(REAL(ep),eindex,1,mpert)
      tt=ep
      DO ipert=1,mpert
         wp(:,ipert)=vr(:,eindex(mpert+1-ipert))
         ep(ipert)=tt(eindex(mpert+1-ipert))
      ENDDO
      CALL zheev('V','U',mpert,wv,mpert,ev,work,lwork,rwork,info) 

c-----------------------------------------------------------------------
c     optionally write netcdf file.
c-----------------------------------------------------------------------
      IF(netcdf_out) CALL dcon_netcdf_out(wp,wv,wt,wt0,ep,ev,et)
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      ! CALL dcon_dealloc(0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE free_run
c-----------------------------------------------------------------------
c     subprogram 2. free_write_msc.
c     writes boundary data for Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE free_write_msc

      CHARACTER(1), PARAMETER :: tab=CHAR(9)
      INTEGER :: itheta,n
      REAL(r8) :: qa
      REAL(r8), DIMENSION(0:mtheta) :: angle,r,z,delta,rfac,theta
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      theta=rzphi%ys
      DO itheta=0,mtheta
         CALL bicube_eval(rzphi,psilim,theta(itheta),0)
         rfac(itheta)=SQRT(rzphi%f(1))
         angle(itheta)=twopi*(theta(itheta)+rzphi%f(2))
         delta(itheta)=-rzphi%f(3)/qlim
      ENDDO
      r=ro+rfac*COS(angle)
      z=zo+rfac*SIN(angle)
c-----------------------------------------------------------------------
c     invert values for nn < 0.
c-----------------------------------------------------------------------
      n=nn
      qa=qlim
      IF(nn < 0)THEN
         qa=-qa
         delta=-delta
         n=-n
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
      WRITE(bin_unit,'(i4,a)')n,tab//tab//"nn"//tab//"nadj"//tab
     $     //"Toroidal harmonic"
      WRITE(bin_unit,'(f13.10,a)')qa,tab//"qa"//tab//"qa1"//tab
     $     //"Safety factor at plasma edge"
c-----------------------------------------------------------------------
c     write arrays.
c-----------------------------------------------------------------------
      WRITE(bin_unit,'(/a/)')"Poloidal Coordinate Theta:"
      WRITE(bin_unit,'(1p,4e18.10)')(1-theta(itheta),
     $     itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Polar Angle Eta:"
      WRITE(bin_unit,'(1p,4e18.10)')(twopi-angle(itheta),
     $     itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Radial Coordinate X:"
      WRITE(bin_unit,'(1p,4e18.10)')(r(itheta),itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Axial Coordinate Z:"
      WRITE(bin_unit,'(1p,4e18.10)')(z(itheta),itheta=mtheta,0,-1)
      WRITE(bin_unit,'(/a/)')"Toroidal Angle Difference Delta:"
      WRITE(bin_unit,'(1p,4e18.10)')(delta(itheta),
     $     itheta=mtheta,0,-1)
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE free_write_msc
c-----------------------------------------------------------------------
c     subprogram 3. free_ahb_prep.
c     prepares boundary data for Allen H. Boozers's feedback problem.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE free_ahb_prep(wp,nmat,smat)

      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: wp
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: nmat,smat

      INTEGER :: itheta,iqty,ipert,info,neq
      REAL(r8) :: rfac,angle,bpfac,btfac,bfac,fac,jac,delpsi,psifac
      REAL(r8), DIMENSION(2,2) :: w
      COMPLEX(r8), DIMENSION(mpert,mpert,2) :: u,du
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      mthsurf=40*mthsurf0*MAX(ABS(mlow),ABS(mhigh))
      ALLOCATE(r(0:mthsurf),z(0:mthsurf),theta(0:mthsurf),
     $     dphi(0:mthsurf),thetas(0:mthsurf,4),project(3,3,0:mthsurf))
      CALL spline_alloc(spl,mthsurf,4)
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      spl%xs=theta
      psifac=psilim
      qsurf=qlim
      q1surf=q1lim
      CALL spline_eval(sq,psilim,0)
c-----------------------------------------------------------------------
c     compute geometric factors.
c-----------------------------------------------------------------------
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psilim,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         angle=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(angle)
         z(itheta)=zo+rfac*SIN(angle)
         dphi(itheta)=rzphi%f(3)
         jac=rzphi%f(4)
c-----------------------------------------------------------------------
c     compute covariant basis vectors.
c-----------------------------------------------------------------------
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
         w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
c-----------------------------------------------------------------------
c     compute project.
c-----------------------------------------------------------------------
         delpsi=SQRT(w(1,1)**2+w(1,2)**2)
         project(1,1,itheta)=1/(delpsi*jac)
         project(2,1,itheta)
     $        =-(w(1,1)*w(2,1)+w(1,2)*w(2,2))/(twopi*r(itheta)*delpsi)
         project(2,2,itheta)=delpsi/(twopi*r(itheta))
         project(3,1,itheta)=rzphi%fx(3)*r(itheta)/jac
         project(3,2,itheta)=rzphi%fy(3)*r(itheta)/jac
         project(3,3,itheta)=twopi*r(itheta)/jac
c-----------------------------------------------------------------------
c     compute alternative theta coordinates.
c-----------------------------------------------------------------------
         bpfac=psio*delpsi/r(itheta)
         btfac=sq%f(1)/(twopi*r(itheta))
         bfac=SQRT(bpfac*bpfac+btfac*btfac)
         fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
         spl%fs(itheta,1)=fac
         spl%fs(itheta,2)=fac/r(itheta)**2
         spl%fs(itheta,3)=fac*bpfac
         !modified for boozer coordinates.
         spl%fs(itheta,4)=fac*bfac**2
      ENDDO
c-----------------------------------------------------------------------
c     compute alternative poloidal coordinates.
c-----------------------------------------------------------------------
      CALL spline_fit(spl,"periodic")
      CALL spline_int(spl)
      DO iqty=1,4
         thetas(:,iqty)=spl%fsi(:,iqty)/spl%fsi(mthsurf,iqty)
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     compute displacements and their derivatives.
c-----------------------------------------------------------------------
      msol=mpert
      neq=4*mpert*msol
      psifac=psilim
      u(:,:,1)=0
      DO ipert=1,mpert
         u(ipert,ipert,1)=1
      ENDDO
      u(:,:,2)=wp*psio**2
      CALL sing_der(neq,psifac,u,du)
      nmat=du(:,:,1)
      smat=-(MATMUL(bmat,du(:,:,1))+MATMUL(cmat,u(:,:,1)))
      CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,smat,mpert,info)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE free_ahb_prep
c-----------------------------------------------------------------------
c     subprogram 4. free_ahb_write.
c     writes boundary data for Allen H. Boozers's feedback problem.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE free_ahb_write(nmat,smat,wt,et)

      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: nmat,smat,wt
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: et

      INTEGER :: itheta,ipert,isol,m
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: chi1
      COMPLEX(r8) :: expfac,expfac0,expfac1,phase
      COMPLEX(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), DIMENSION(mpert,mpert) :: xin,xis
      COMPLEX(r8), DIMENSION(mpert,mpert,3) :: jqvec
      COMPLEX(r8), DIMENSION(mpert,0:mthsurf,3) :: bvec0,bvec
c-----------------------------------------------------------------------
c     compute Fourier components of magnetic field.
c-----------------------------------------------------------------------
      xin=MATMUL(nmat,wt)
      xis=MATMUL(smat,wt)
      mvec=(/(m,m=mlow,mhigh)/)
      chi1=twopi*psio
      singfac=chi1*twopi*ifac*(mvec-nn*qsurf)
      DO isol=1,mpert
         jqvec(:,isol,1)=wt(:,isol)*singfac
         jqvec(:,isol,2)=-chi1*xin(:,isol)-twopi*ifac*nn*xis(:,isol)
         jqvec(:,isol,3)=-chi1*(qsurf*xin(:,isol)+q1surf*wt(:,isol))
     $        -twopi*ifac*mvec*xis(:,isol)
      ENDDO
c-----------------------------------------------------------------------
c     transform to configuration space.
c-----------------------------------------------------------------------
      expfac0=EXP(twopi*ifac*theta(1))
      expfac1=1
      bvec0=0
      DO itheta=0,mthsurf-1
         expfac=EXP(ifac*(twopi*mlow*theta(itheta)+nn*dphi(itheta)))
         DO ipert=1,mpert
            bvec0(:,itheta,:)=bvec0(:,itheta,:)+jqvec(ipert,:,:)*expfac
            expfac=expfac*expfac1
         ENDDO
         expfac1=expfac1*expfac0
      ENDDO
      bvec0(:,mthsurf,:)=bvec0(:,0,:)
c-----------------------------------------------------------------------
c     compute orthogonal components of perturbed magnetic field.
c-----------------------------------------------------------------------
      DO isol=1,msol
         bvec(isol,:,1)=project(1,1,:)*bvec0(isol,:,1)
         bvec(isol,:,2)
     $        =project(2,1,:)*bvec0(isol,:,1)
     $        +project(2,2,:)*bvec0(isol,:,2)
         bvec(isol,:,3)
     $        =project(3,1,:)*bvec0(isol,:,1)
     $        +project(3,2,:)*bvec0(isol,:,2)
     $        +project(3,3,:)*bvec0(isol,:,3)
         phase=bvec(isol,0,1)/ABS(bvec(isol,0,1))
         bvec(isol,:,:)=bvec(isol,:,:)/phase
      ENDDO
c-----------------------------------------------------------------------
c     open dcon_surf.out, write geometry and eigenvalues.
c-----------------------------------------------------------------------
      CALL ascii_open(bin_unit,"dcon_surf.out","UNKNOWN")
      WRITE(bin_unit,'(a/a/a/)')"Output from DCON, ",
     $     "total perturbed energy eigenvalues",
     $     "and normal magnetic field eigenvectors"
      WRITE(bin_unit,'(3(a,i4),a,1p,e16.8)')
     $     "mthsurf = ",mthsurf,", msol = ",msol,", nn = ",nn,
     $     ", qsurf = ",qsurf
      WRITE(bin_unit,'(/a//(1p,5e16.8))')"radial position r:",r
      WRITE(bin_unit,'(/a//(1p,5e16.8))')"axial position z:",z
      WRITE(bin_unit,'(/a//(1p,5e16.8))')
     $     "total energy eigenvalues et:",et*psio**2/(2*mu0)
c-----------------------------------------------------------------------
c     write perturbed normal magnetic field.
c-----------------------------------------------------------------------
      DO isol=1,mpert
         WRITE(bin_unit,'(/a,i3,a/)')"normal magnetic eigenvector"
     $        //", isol = ",isol,", cos factor:"
         WRITE(bin_unit,'(1p,5e16.8)')REAL(bvec(isol,:,1))
         WRITE(bin_unit,'(/a,i3,a/)')"normal magnetic eigenvector"
     $        //", isol = ",isol,", sin factor:"
         WRITE(bin_unit,'(1p,5e16.8)')AIMAG(bvec(isol,:,1))
      ENDDO
c-----------------------------------------------------------------------
c     write perturbed poloidal magnetic field.
c-----------------------------------------------------------------------
      DO isol=1,mpert
         WRITE(bin_unit,'(/a,i3,a/)')"poloidal magnetic eigenvector"
     $        //", isol = ",isol,", cos factor:"
         WRITE(bin_unit,'(1p,5e16.8)')REAL(bvec(isol,:,2))
         WRITE(bin_unit,'(/a,i3,a/)')"poloidal magnetic eigenvector,"
     $        //" isol = ",isol,", sin factor:"
         WRITE(bin_unit,'(1p,5e16.8)')AIMAG(bvec(isol,:,2))
      ENDDO
c-----------------------------------------------------------------------
c     write perturbed toroidal magnetic field.
c-----------------------------------------------------------------------
      DO isol=1,mpert
         WRITE(bin_unit,'(/a,i3,a/)')"toroidal magnetic eigenvector,"
     $        //" isol = ",isol,", cos factor:"
         WRITE(bin_unit,'(1p,5e16.8)')REAL(bvec(isol,:,3))
         WRITE(bin_unit,'(/a,i3,a/)')"toroidal magnetic eigenvector"
     $        //", isol = ",isol,", sin factor:"
         WRITE(bin_unit,'(1p,5e16.8)')AIMAG(bvec(isol,:,3))
      ENDDO
      CALL ascii_close(bin_unit)
c-----------------------------------------------------------------------
c     draw graphs.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"ahb.bin","UNKNOWN","REWIND","none")
      DO isol=1,mpert
         DO itheta=0,mthsurf
            WRITE(bin_unit)REAL(theta(itheta),4),
     $           REAL(REAL(bvec(isol,itheta,:)),4),
     $           REAL(AIMAG(bvec(isol,itheta,:)),4),
     $           REAL(thetas(itheta,:),4),
     $           REAL(thetas(itheta,:)-theta(itheta),4)
         ENDDO
         WRITE(bin_unit)
         IF(isol == msol_ahb)EXIT
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE free_ahb_write

c-----------------------------------------------------------------------
c     subprogram 5. free_get_wvac.
c     get vacuum response matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE free_get_wvac()
      LOGICAL, PARAMETER :: complex_flag=.FALSE.
      INTEGER :: ipert
      REAL(r8) :: kernelsignin
      REAL(r8), DIMENSION(mpert) :: singfac
      REAL(r8), DIMENSION(:,:), POINTER :: grri,xzpts
      LOGICAL :: wall_flag=.FALSE., farwal_flag=.FALSE.

c-----------------------------------------------------------------------
c     allocate vacuum matrix.
c-----------------------------------------------------------------------      
      IF (ALLOCATED(wvac)) RETURN
      ALLOCATE (wvac(mpert,mpert))
c-----------------------------------------------------------------------
c     write file for mscvac, prepare input for ahb, and deallocate.
c-----------------------------------------------------------------------
      CALL free_write_msc
c-----------------------------------------------------------------------
c     compute vacuum response matrix.
c-----------------------------------------------------------------------
      kernelsignin=1.0
      ALLOCATE(grri(2*(mthvac+5),mpert*2),xzpts(mthvac+5,4))
      CALL mscvac(wvac,mpert,mtheta,mthvac,complex_flag,kernelsignin,
     $     wall_flag,farwal_flag,grri,xzpts)
      DEALLOCATE(grri,xzpts)

      CALL system("rm -f ahg2msc.out")
      singfac=mlow-nn*qlim+(/(ipert,ipert=0,mpert-1)/)
      DO ipert=1,mpert
         wvac(ipert,:)=wvac(ipert,:)*singfac
         wvac(:,ipert)=wvac(:,ipert)*singfac
      ENDDO      
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------      
      RETURN
      END SUBROUTINE free_get_wvac
      END MODULE free_mod
