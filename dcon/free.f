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
c-----------------------------------------------------------------------
c     subprogram 0. free_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE free_mod
      USE ode_mod
      IMPLICIT NONE

      INTEGER :: msol_ahb=-1
      INTEGER, PRIVATE :: mthsurf
      REAL(r8), PRIVATE :: qsurf,q1surf
      REAL(r8), DIMENSION(:), POINTER, PRIVATE :: theta,dphi,r,z
      REAL(r8), DIMENSION(:,:), POINTER, PRIVATE :: thetas
      REAL(r8), DIMENSION(:,:,:), POINTER, PRIVATE :: project
c-----------------------------------------------------------------------
c     <MODIFIED>
c-----------------------------------------------------------------------
      INTEGER :: nfm2,nths2
      REAL(8), DIMENSION(:,:), pointer :: grri
c-----------------------------------------------------------------------

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. free_run.
c     computes plasma, vacuum, and total potential energies.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations. <MODIFIED>
c-----------------------------------------------------------------------
      SUBROUTINE free_run(plasma1,vacuum1,total1,nzero)

      REAL(r8), INTENT(OUT) :: plasma1,vacuum1,total1
      INTEGER, INTENT(IN) :: nzero

      LOGICAL, PARAMETER :: normalize=.TRUE.
      CHARACTER(1), DIMENSION(mpert,msol) :: star
      INTEGER :: ipert,jpert,isol,info,lwork
      INTEGER, DIMENSION(mpert) :: ipiv,m
      INTEGER, DIMENSION(1) :: imax
      REAL(r8) :: v1
      REAL(r8), DIMENSION(mpert) :: ep,ev,et,singfac
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8) :: phase,norm
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: wp,wv,wt,rt,temp,wpt,wvt
      COMPLEX(r8), DIMENSION(mpert,mpert) :: nmat,smat
      CHARACTER(24), DIMENSION(mpert) :: message
      LOGICAL, PARAMETER :: complex_flag=.TRUE.
      REAL(r8) :: kernelsignin
      INTEGER :: vac_unit
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"Energies: plasma = ",1p,e10.3,", vacuum = ",e10.3,
     $     ", total = ",e10.3)
 20   FORMAT(/3x,"isol",3x,"plasma",5x,"vacuum",5x,"total"/)
 30   FORMAT(i6,1p,3e11.3,a)
 40   FORMAT(/3x,"isol",2x,"imax",3x,"plasma",5x,"vacuum",5x,"total"/)
 50   FORMAT(2i6,1p,3e11.3,a)
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
c     write file for mscvac, prepare input for ahb, and deallocate.
c-----------------------------------------------------------------------
      CALL free_write_msc
      IF(ahb_flag)CALL free_ahb_prep(wp,nmat,smat)
      CALL spline_eval(sq,psilim,0)
      v1=sq%f(3)
      CALL dcon_dealloc
c-----------------------------------------------------------------------
c     compute vacuum response matrix.<MODIFIED>
c-----------------------------------------------------------------------
      vac_unit=4
      kernelsignin=-1.0
      CALL mscvac(wv,mpert,mtheta,mthvac,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(grri(nths2,nfm2))
      CALL grrget(nfm2,nths2,grri)
      CALL bin_open(vac_unit,"vacuum.bin","UNKNOWN","REWIND","none")
      WRITE(vac_unit)SIZE(grri,1),SIZE(grri,2)
      WRITE(vac_unit)grri
      DEALLOCATE(grri)

      kernelsignin=1.0
      CALL mscvac(wv,mpert,mtheta,mthvac,nfm2,nths2,complex_flag,
     $     kernelsignin)
      ALLOCATE(grri(nths2,nfm2))
      CALL grrget(nfm2,nths2,grri)
      WRITE(vac_unit)SIZE(grri,1),SIZE(grri,2)
      WRITE(vac_unit)grri
      DEALLOCATE(grri)
      CALL bin_close(vac_unit)

      singfac=mlow-nn*qlim+(/(ipert,ipert=0,mpert-1)/)
      DO ipert=1,mpert
         wv(ipert,:)=wv(ipert,:)*singfac
         wv(:,ipert)=wv(:,ipert)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     compute energy eigenvalues.<MODIFIED>
c-----------------------------------------------------------------------
      wt=wp+wv
      rt=wt
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,wt,mpert,et,work,lwork,rwork,info)
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
         imax=MAXLOC(ABS(wt(:,isol)))
         phase=ABS(wt(imax(1),isol))/wt(imax(1),isol)
         wt(:,isol)=wt(:,isol)*phase
         star(:,isol)=' '
         star(imax(1),isol)='*'
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
c-----------------------------------------------------------------------
c     write data for ahb and deallocate.
c-----------------------------------------------------------------------
      IF(ahb_flag)THEN
         CALL free_ahb_write(nmat,smat,wt,et)
         DEALLOCATE(r,z,theta,dphi,thetas,project)
      ENDIF
c-----------------------------------------------------------------------
c     save eigenvalues and eigenvectors to file. <MODIFIED>
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         WRITE(euler_bin_unit)3
         WRITE(euler_bin_unit)et
         WRITE(euler_bin_unit)wt
         WRITE(euler_bin_unit)ep
         WRITE(euler_bin_unit)rt
      ENDIF
c-----------------------------------------------------------------------
c     write to screen and copy to output.
c-----------------------------------------------------------------------
      WRITE(*,10)ep(1),ev(1),et(1)
      plasma1=ep(1)
      vacuum1=ev(1)
      total1=et(1)
c-----------------------------------------------------------------------
c     write eigenvalues to file.
c-----------------------------------------------------------------------
      message=""
c$$$      message(1:nzero)="  internal instability"
      WRITE(out_unit,'(/1x,a)')"Total Energy Eigenvalues:"
      WRITE(out_unit,20)
      WRITE(out_unit,30)(isol,ep(isol),ev(isol),et(isol),
     $     TRIM(message(isol)),isol=1,mpert)
      WRITE(out_unit,20)
c-----------------------------------------------------------------------
c     write eigenvectors to file.
c-----------------------------------------------------------------------
      WRITE(out_unit,*)"Total Energy Eigenvectors:"
      m=mlow+(/(isol,isol=0,mpert-1)/)
      DO isol=1,mpert
         WRITE(out_unit,40)
         WRITE(out_unit,50)isol,imax(1),ep(isol),ev(isol),et(isol),
     $        TRIM(message(isol))
         WRITE(out_unit,60)
         WRITE(out_unit,70)(ipert,m(ipert),wt(ipert,isol),
     $        ABS(wt(ipert,isol)),star(ipert,isol),ipert=1,mpert)
         WRITE(out_unit,60)
      ENDDO
c-----------------------------------------------------------------------
c     write the plasma matrix. (MODIFIED)
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/1x,a/)')"Plasma Energy Matrix:"
c-----Temporalily changed
c      DO isol=1,mpert
      DO isol=1,2
         WRITE(out_unit,'(1x,2(a,i3))')"isol = ",isol,", m = ",m(isol)
         WRITE(out_unit,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
         WRITE(out_unit,'(i3,1p,3e13.5)')
     $        (ipert,wp(ipert,isol),ABS(wp(ipert,isol)),ipert=1,mpert)
         WRITE(out_unit,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
      ENDDO
c-----------------------------------------------------------------------
c     compute and print separate plasma and vacuum eigenvalues. (MOD)
c-----------------------------------------------------------------------
      CALL zheev('V','U',mpert,wp,mpert,ep,work,lwork,rwork,info)
      CALL zheev('V','U',mpert,wv,mpert,ev,work,lwork,rwork,info)
c      WRITE(out_unit,*)"Separate Energy Eigenvalues:"
c      WRITE(out_unit,80)
c      WRITE(out_unit,90)(isol,ep(isol),ev(isol),isol=1,mpert)
c      WRITE(out_unit,80)
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
         spl%fs(itheta,4)=fac*bpfac**2
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
      REAL(r8), DIMENSION(mpert), INTENT(IN) :: et

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
      END MODULE free_mod
