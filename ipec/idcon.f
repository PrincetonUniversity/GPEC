c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     read and preprocess ideal DCON/VACUUM data
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. idcon_mod.
c     1. idcon_read.
c     2. idcon_transform
c     3. idcon_build
c     4. idcon_matrices
c     5. idcon_vacuum
c-----------------------------------------------------------------------
c     subprogram 0. idcon_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE idcon_mod
      USE ipglobal_mod
      USE ismath_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. idcon_read.
c     reads dcon output, psi_in.bin and euler.bin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_read(lpsixy)

      INTEGER, INTENT(IN) :: lpsixy
      CHARACTER(128) :: message
      INTEGER :: m,data_type,ifix,ios,msol,istep,ising,itheta
      REAL(r8) :: sfac0

      REAL(r4), DIMENSION(:,:), POINTER :: rgarr,zgarr,psigarr
c-----------------------------------------------------------------------
c     open euler.bin and read header.
c-----------------------------------------------------------------------
      WRITE(*,*)"Reading dcon eigenfuctions"
      CALL bin_open(in_unit,idconfile,"OLD","REWIND","none")
      READ(in_unit)mlow,mhigh,nn,mpsi,mtheta,ro,zo
      READ(in_unit)mband,mthsurf0,mthvac,psio,psilow,psilim,qlim,
     $     singfac_min
      READ(in_unit)power_b,power_r,power_bp
      IF ((power_b==0).AND.(power_bp==0).AND.(power_r==0)) THEN
         jac_type="hamada"
      ELSE IF ((power_b==0).AND.(power_bp==0).AND.(power_r==2)) THEN
         jac_type="pest"
      ELSE IF ((power_b==0).AND.(power_bp==1).AND.(power_r==0)) THEN
         jac_type="equal_arc"
      ELSE IF ((power_b==2).AND.(power_bp==0).AND.(power_r==0)) THEN  
         jac_type="boozer"
      ELSE 
         jac_type="other"
      ENDIF
      IF (jac_in=="") jac_in=jac_type
      IF (jac_out=="") jac_out=jac_type
      chi1=twopi*psio
      mpert=mhigh-mlow+1
      lmpert=lmhigh-lmlow+1
      mthsurf=mthvac
      ALLOCATE(r(0:mthsurf),z(0:mthsurf),theta(0:mthsurf))
      ALLOCATE(mfac(mpert),singfac(mpert))
      ALLOCATE(lmfac(lmpert))
      ALLOCATE(edge_mn(mpert),edge_fun(0:mthsurf))
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      mfac=(/(m,m=mlow,mhigh)/)
      lmfac=(/(m,m=lmlow,lmhigh)/)
      IF (nn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)')nn
         sn=TRIM(ADJUSTR(sn))
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')nn
      ENDIF
c-----------------------------------------------------------------------
c     only accept mband=0.
c-----------------------------------------------------------------------
      IF (mband /= (mpert-1)) THEN
         WRITE(message,'(a)')"IPEC needs full band matrix"
         CALL ipec_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     read equilibrium on flux coordinates.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq,mpsi,4)
      CALL bicube_alloc(rzphi,mpsi,mtheta,4)

      rzphi%periodic(2)=.TRUE.
      READ(in_unit)sq%xs,sq%fs,sq%fs1
      READ(in_unit)rzphi%xs,rzphi%ys,
     $     rzphi%fs,rzphi%fsx,rzphi%fsy,rzphi%fsxy,
     $     rzphi%x0,rzphi%y0,rzphi%xpower,rzphi%ypower
      mstep=-1
      mfix=0
      msing=0
c-----------------------------------------------------------------------
c     count solutions in euler.bin.
c-----------------------------------------------------------------------
      WRITE(*,*)"Counting and reading dcon solutions"
      DO
         READ(UNIT=in_unit,IOSTAT=ios)data_type
         IF(ios /= 0)EXIT
         SELECT CASE(data_type)
         CASE(1)
            mstep=mstep+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(2)
            mfix=mfix+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(3)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(4)
            msing=msing+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE DEFAULT
            WRITE(message,'(a,i1,a,i4)')"Cannot recognize data_type = ",
     $           data_type,", at istep = ",istep
            CALL ipec_stop(message)
         END SELECT
      ENDDO
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays and prepare to read data.
c-----------------------------------------------------------------------
      WRITE(*,*)"mlow = ",mlow,", mhigh = ",mhigh,", mpert = ",mpert
      WRITE(*,*)"mstep = ",mstep,", mfix = ",mfix,", msing = ",msing
      ALLOCATE(psifac(0:mstep),rhofac(0:mstep),qfac(0:mstep),
     $     soltype(0:mstep),singtype(msing))
      ALLOCATE(fixstep(0:mfix+1),fixtype(0:mfix),sing_flag(mfix))
      ALLOCATE(et(mpert),ep(mpert),ee(mpert),wt(mpert,mpert))
      fixstep(0)=0
      fixstep(mfix+1)=mstep
      CALL bin_open(in_unit,idconfile,"OLD","REWIND","none")
      READ(in_unit)mlow,mhigh,nn
      READ(in_unit)mband,mthsurf0,mthvac,psio,psilow,psilim,qlim,
     $     singfac_min
      READ(in_unit)power_b,power_r,power_bp
      istep=-1
      ifix=0
      ising=0
c-----------------------------------------------------------------------
c     read solution data in euler.bin.
c-----------------------------------------------------------------------
      DO
         READ(UNIT=in_unit,IOSTAT=ios)data_type
         IF(ios /= 0)EXIT
         SELECT CASE(data_type)
         CASE(1)
            istep=istep+1
            READ(UNIT=in_unit)psifac(istep),qfac(istep),
     $           soltype(istep)%msol
            ALLOCATE(soltype(istep)%u(mpert,soltype(istep)%msol,2))
            READ(UNIT=in_unit)soltype(istep)%u
         CASE(2)
            ifix=ifix+1
            fixstep(ifix)=istep
            READ(UNIT=in_unit)sing_flag(ifix),fixtype(ifix)%msol
            ALLOCATE(fixtype(ifix)%fixfac
     $           (fixtype(ifix)%msol,fixtype(ifix)%msol))
            ALLOCATE(fixtype(ifix)%index(fixtype(ifix)%msol))
            READ(UNIT=in_unit)fixtype(ifix)%fixfac,fixtype(ifix)%index
         CASE(3)
            READ(UNIT=in_unit)ep
            READ(UNIT=in_unit)et
            READ(UNIT=in_unit)wt
         CASE(4)
            ising=ising+1
            singtype(ising)%jfix=ifix
            singtype(ising)%jpert=INT(nn*qfac(istep)+.5)-mlow+1
            READ(UNIT=in_unit)singtype(ising)%psifac,
     $           singtype(ising)%q,singtype(ising)%q1
c-----------------------------------------------------------------------
c     information required for resistive MHD computations.
c-----------------------------------------------------------------------
            READ(UNIT=in_unit)msol
            singtype(ising)%msol_l=msol
            ALLOCATE(singtype(ising)%ca_l(mpert,msol,2))
            READ(UNIT=in_unit)singtype(ising)%ca_l
            READ(UNIT=in_unit)msol
            singtype(ising)%msol_r=msol
            ALLOCATE(singtype(ising)%ca_r(mpert,msol,2))
            READ(UNIT=in_unit)singtype(ising)%ca_r
            READ(UNIT=in_unit)
     $           singtype(ising)%restype%e,
     $           singtype(ising)%restype%f,
     $           singtype(ising)%restype%h,
     $           singtype(ising)%restype%m,
     $           singtype(ising)%restype%g,
     $           singtype(ising)%restype%k,
     $           singtype(ising)%restype%eta,
     $           singtype(ising)%restype%rho,
     $           singtype(ising)%restype%taua,
     $           singtype(ising)%restype%taur
         END SELECT
      ENDDO
      IF (psifac(mstep)<psilim-(1e-4)) THEN
         WRITE(message,'(a)')"Terminated by zero crossing"
         CALL ipec_stop(message)
      ENDIF
      rhofac=SQRT(psifac)
c-----------------------------------------------------------------------
c     normalize plasma/vacuum eigenenergy and eigenfunctions.
c-----------------------------------------------------------------------
      et=et/(mu0*2.0)*psio**2*(chi1*mili)**2
      ep=ep/(mu0*2.0)*psio**2*(chi1*mili)**2
      ee=et-ep
      wt=wt*chi1*mili
c-----------------------------------------------------------------------
c     modify Lundquist numbers.
c-----------------------------------------------------------------------
      sfac0=1
      DO ising=1,msing
         singtype(ising)%restype%taur=singtype(ising)%restype%taur*sfac0
      ENDDO
c-----------------------------------------------------------------------
c     close data file.
c-----------------------------------------------------------------------
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     read psi_in.bin.
c-----------------------------------------------------------------------
      IF (psixy == 1) THEN
         WRITE(*,*)"Reading axisymmetric equilibrium solutions"
         CALL bin_open(in_unit,ieqfile,"OLD","REWIND","none")
         READ(in_unit)
         READ(in_unit)mr,mz
         CALL bicube_alloc(psi_in,mr,mz,1)
         psi_in%name="equilibrium psi"
         psi_in%xtitle="r"
         psi_in%ytitle="z"
         psi_in%title=" normalized psi "
         ALLOCATE(rgarr(0:mr,0:mz),zgarr(0:mr,0:mz),psigarr(0:mr,0:mz))
         READ(in_unit)rgarr,zgarr
         psi_in%xs=rgarr(:,0)
         psi_in%ys=zgarr(0,:)
         READ(in_unit)psigarr
         psi_in%fs(:,:,1)=1-psigarr/psio
         CALL bin_close(in_unit)
         DEALLOCATE(rgarr,zgarr,psigarr)
         CALL bicube_fit(psi_in,"extrap","extrap")
         WRITE(*,*)"mr = ",mr,", mz = ",mz
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_read
c-----------------------------------------------------------------------
c     subprogram 2. idcon_transform.
c     build fixup matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_transform

      INTEGER :: ifix,isol,jsol,ksol
      LOGICAL, DIMENSION(mpert) :: mask
      COMPLEX(r8), DIMENSION(mpert,mpert) :: ident,temp
c-----------------------------------------------------------------------
c     create identity matrix.
c-----------------------------------------------------------------------
      ident=0
      DO isol=1,mpert
         ident(isol,isol)=1
      ENDDO
c-----------------------------------------------------------------------
c     compute gaussian reduction matrices.
c-----------------------------------------------------------------------
      DO ifix=1,mfix
         ALLOCATE(fixtype(ifix)%gauss(mpert,mpert))
         fixtype(ifix)%gauss=ident
         mask=.TRUE.
         DO isol=1,mpert
            ksol=fixtype(ifix)%index(isol)
            mask(ksol)=.FALSE.
            temp=ident
            DO jsol=1,mpert
               IF(mask(jsol))
     $              temp(ksol,jsol)=fixtype(ifix)%fixfac(ksol,jsol)
            ENDDO
            fixtype(ifix)%gauss=MATMUL(fixtype(ifix)%gauss,temp)
         ENDDO
         IF(sing_flag(ifix))
     $        fixtype(ifix)%gauss(:,fixtype(ifix)%index(1))=0
      ENDDO
c-----------------------------------------------------------------------
c     concatenate gaussian reduction matrices.
c-----------------------------------------------------------------------
      ALLOCATE(fixtype(mfix)%transform(mpert,mpert))
      fixtype(mfix)%transform=ident
      DO ifix=mfix-1,0,-1
         ALLOCATE(fixtype(ifix)%transform(mpert,mpert))
         fixtype(ifix)%transform
     $        =MATMUL(fixtype(ifix+1)%gauss,
     $        fixtype(ifix+1)%transform)
      ENDDO
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DO ifix=1,mfix
         DEALLOCATE(fixtype(ifix)%gauss)
      ENDDO
      CALL cspline_alloc(u1,mstep,mpert)
      CALL cspline_alloc(u2,mstep,mpert)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_transform
c-----------------------------------------------------------------------
c     subprogram 3. idcon_build.
c     builds ideal euler-Lagrange solutions.
c     __________________________________________________________________
c     egnum  : a label of eigenmode
c     xwpimn : arbitrary xipsi function on the control surface.
c            : only works when edge_flag is activated.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_build(egnum,xspimn)

      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspimn

      INTEGER :: istep,ifix,jfix,kfix,ieq,info
      INTEGER, DIMENSION(mpert) :: ipiv
      COMPLEX(r8), DIMENSION(mpert) :: uedge,temp1
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp2
c-----------------------------------------------------------------------
c     construct uedge.
c-----------------------------------------------------------------------
      IF(edge_flag)THEN
         uedge=xspimn
      ELSE
         uedge=wt(:,egnum)
      ENDIF
      temp2=soltype(mstep)%u(:,1:mpert,1)
      CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
      CALL zgetrs('N',mpert,1,temp2,mpert,ipiv,uedge,mpert,info)
c-----------------------------------------------------------------------
c     construct eigenfunctions.
c-----------------------------------------------------------------------
      jfix=0
      DO ifix=0,mfix
         temp1=MATMUL(fixtype(ifix)%transform,uedge)
         kfix=fixstep(ifix+1)
         DO istep=jfix,kfix
            u1%fs(istep,:)
     $           =MATMUL(soltype(istep)%u(:,1:mpert,1),temp1)
            u2%fs(istep,:)
     $           =MATMUL(soltype(istep)%u(:,1:mpert,2),temp1)
         ENDDO         
         jfix=kfix+1
      ENDDO
c-----------------------------------------------------------------------
c     fit the functions.
c-----------------------------------------------------------------------
      u1%xs=psifac
      u2%xs=psifac
      CALL cspline_fit(u1,"extrap")
      CALL cspline_fit(u2,"extrap")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_build
c-----------------------------------------------------------------------
c     subprogram 4. idcon_matrices.
c     reconstructs metric tensors from dcon.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_metric
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: ipsi,itheta
      REAL(r8) :: rr,w11,w12,delpsi

      WRITE(*,*)"Recontructing flux functions and metric tensors"
c-----------------------------------------------------------------------
c     set up fourier-spline type for metric tensors.
c-----------------------------------------------------------------------
      CALL fspline_alloc(metric,mpsi,mtheta,mband,8)
      metric%xs=rzphi%xs
      metric%ys=rzphi%ys*twopi
      metric%name="metric"
      metric%xtitle="psi"
      metric%ytitle="theta"
      metric%title=(/" g11  "," g22  "," g33  "," g23  "," g31  ",
     $     " g12  "," jmat "," jmat1"/)
c-----------------------------------------------------------------------
c     set up bicube type for equilibrium mod b.
c-----------------------------------------------------------------------
      CALL bicube_alloc(eqfun,mpsi,mtheta,1)
      eqfun%xs=rzphi%xs
      eqfun%ys=rzphi%ys
      eqfun%name="equilibrium funs"
      eqfun%xtitle="psi"
      eqfun%ytitle="theta"
      eqfun%title=" modb "      
c-----------------------------------------------------------------------
c     begin loop over nodes.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         CALL spline_eval(sq,sq%xs(ipsi),1)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(itheta/REAL(mtheta,r8)+rzphi%f(2))
            rr=ro+rfac*COS(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
            w11=(1+rzphi%fy(2))*twopi**2*rfac*rr/jac
            w12=-rzphi%fy(1)*pi*rr/(rfac*jac)
            delpsi=SQRT(w11**2+w12**2)
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(1,3)=rzphi%fx(3)*rr
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(2,3)=rzphi%fy(3)*rr
            v(3,3)=twopi*rr
c-----------------------------------------------------------------------
c     compute metric tensor components.
c-----------------------------------------------------------------------
            metric%fs(ipsi,itheta,1)=SUM(v(1,:)**2)/jac
            metric%fs(ipsi,itheta,2)=SUM(v(2,:)**2)/jac
            metric%fs(ipsi,itheta,3)=v(3,3)*v(3,3)/jac
            metric%fs(ipsi,itheta,4)=v(2,3)*v(3,3)/jac
            metric%fs(ipsi,itheta,5)=v(3,3)*v(1,3)/jac
            metric%fs(ipsi,itheta,6)=SUM(v(1,:)*v(2,:))/jac
            metric%fs(ipsi,itheta,7)=jac
            metric%fs(ipsi,itheta,8)=jac1
c-----------------------------------------------------------------------
c     compute equilibrium mod b.
c     sq%f(1) = twopi*f
c-----------------------------------------------------------------------
            eqfun%fs(ipsi,itheta,1)=SQRT(
     $           (sq%f(1)**2+chi1**2*delpsi**2)/(twopi*rr)**2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit spline types.
c-----------------------------------------------------------------------
      IF(fft_flag)THEN
         CALL fspline_fit_2(metric,"extrap",.TRUE.)
      ELSE
         CALL fspline_fit_1(metric,"extrap",.TRUE.)
      ENDIF
c-----------------------------------------------------------------------
c     fit bicube types.
c-----------------------------------------------------------------------
      CALL bicube_fit(eqfun,"extrap","periodic")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_metric
c-----------------------------------------------------------------------
c     subprogram 5. idcon_matrix
c     compute abcfgk matrix.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_matrix(psi)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi

      CHARACTER(128) :: message
      INTEGER :: ipert,jpert,m1,m2,m,dm,info,iqty
      REAL(r8) :: jtheta,nq,singfac1,singfac2,rm
      INTEGER, DIMENSION(mpert) :: ipiva
      COMPLEX(r8), DIMENSION(mpert*mpert) :: work

      COMPLEX(r8), DIMENSION(-mband:mband) ::
     $     g11,g22,g33,g23,g31,g12,jmat,jmat1,imat
      COMPLEX(r8), DIMENSION((mband+1)*(2*mpert-mband)/2) :: fvec,gvec
      COMPLEX(r8), DIMENSION((2*mband+1)*mpert) :: kvec

      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb
      COMPLEX(r8), DIMENSION(mpert,mpert) :: dmat,emat,fmat,gmat,hmat,
     $     kmat,temp1,temp2,temp3

      imat=0
      imat(0)=1
      CALL spline_eval(sq,psi,1)
      CALL cspline_eval(metric%cs,psi,0)
      p1=sq%f1(2)
      q=sq%f(4)
      q1=sq%f1(4)
      nq=nn*q
      jtheta=-sq%f1(1)
c-----------------------------------------------------------------------
c     compute lower half of matrices.
c-----------------------------------------------------------------------
      g11(0:-mband:-1)=metric%cs%f(1:mband+1)
      g22(0:-mband:-1)=metric%cs%f(mband+2:2*mband+2)
      g33(0:-mband:-1)=metric%cs%f(2*mband+3:3*mband+3)
      g23(0:-mband:-1)=metric%cs%f(3*mband+4:4*mband+4)
      g31(0:-mband:-1)=metric%cs%f(4*mband+5:5*mband+5)
      g12(0:-mband:-1)=metric%cs%f(5*mband+6:6*mband+6)
      jmat(0:-mband:-1)=metric%cs%f(6*mband+7:7*mband+7)
      jmat1(0:-mband:-1)=metric%cs%f(7*mband+8:8*mband+8)
c-----------------------------------------------------------------------
c     compute upper half of matrices.
c-----------------------------------------------------------------------
      g11(1:mband)=CONJG(g11(-1:-mband:-1))
      g22(1:mband)=CONJG(g22(-1:-mband:-1))
      g33(1:mband)=CONJG(g33(-1:-mband:-1))
      g23(1:mband)=CONJG(g23(-1:-mband:-1))
      g31(1:mband)=CONJG(g31(-1:-mband:-1))
      g12(1:mband)=CONJG(g12(-1:-mband:-1))
      jmat(1:mband)=CONJG(jmat(-1:-mband:-1))
      jmat1(1:mband)=CONJG(jmat1(-1:-mband:-1))
c-----------------------------------------------------------------------
c     begin loops over perturbed fourier components.
c-----------------------------------------------------------------------
      ipert=0
      DO m1=mlow,mhigh
         ipert=ipert+1
         singfac1=m1-nq
         DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
            m2=m1+dm
            singfac2=m2-nq
            jpert=ipert+dm
c-----------------------------------------------------------------------
c     construct primitive matrices by using metric tensors.
c-----------------------------------------------------------------------
            amat(ipert,jpert)=twopi**2*(nn*nn*g22(dm)
     $           +nn*(m1+m2)*g23(dm)+m1*m2*g33(dm))
            bmat(ipert,jpert)=-twopi*ifac*chi1
     $           *(nn*g22(dm)+(m1+nq)*g23(dm)+m1*q*g33(dm))
            cmat(ipert,jpert)=twopi*ifac*(
     $           twopi*ifac*chi1*singfac2*(nn*g12(dm)+m1*g31(dm))
     $           -q1*chi1*(nn*g23(dm)+m1*g33(dm)))
     $           -twopi*ifac*(jtheta*singfac1*imat(dm)
     $           +nn*p1/chi1*jmat(dm))
            dmat(ipert,jpert)=twopi*chi1*(g23(dm)+g33(dm)*m1/nn)
            emat(ipert,jpert)=-chi1/nn*(q1*chi1*g33(dm)
     $           -twopi*ifac*chi1*g31(dm)*singfac2
     $           +jtheta*imat(dm))
            hmat(ipert,jpert)=(q1*chi1)**2*g33(dm)
     $           +(twopi*chi1)**2*singfac1*singfac2*g11(dm)
     $           -twopi*ifac*chi1*dm*q1*chi1*g31(dm)
     $           +jtheta*q1*chi1*imat(dm)+p1*jmat1(dm)
            fmat(ipert,jpert)=(chi1/nn)**2*g33(dm)
            kmat(ipert,jpert)=twopi*ifac*chi1*(g23(dm)+g33(dm)*m1/nn)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     factor a.
c-----------------------------------------------------------------------
      temp1=amat
      CALL zhetrf('L',mpert,temp1,mpert,ipiva,work,mpert*mpert,info)
      IF(info /= 0)THEN
         WRITE(message,'(a,e12.3,a,i3,a)')
     $        "zhetrf: amat singular at psi = ",psi,
     $        ", ipert = ",info,", reduce delta_mband"
         CALL ipec_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     compute composite matrices fgk, different from dcon notes.
c-----------------------------------------------------------------------
      temp2=dmat
      temp3=cmat
      CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiva,temp2,mpert,info)
      CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiva,temp3,mpert,info)
      fmat=fmat-MATMUL(CONJG(TRANSPOSE(dmat)),temp2)
      kmat=emat-MATMUL(CONJG(TRANSPOSE(kmat)),temp3)
      gmat=hmat-MATMUL(CONJG(TRANSPOSE(cmat)),temp3)
c-----------------------------------------------------------------------
c     transfer f to banded matrix.
c-----------------------------------------------------------------------
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            fmatb(1+ipert-jpert,jpert)=fmat(ipert,jpert)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     factor f.
c-----------------------------------------------------------------------
      CALL zpbtrf('L',mpert,mband,fmatb,mband+1,info)
      IF(info /= 0) THEN
         WRITE(message,'(a,e12.3,a,i3,a)')
     $        "zpbtrf: fmat singular at psi = ",psi,
     $        ", ipert = ",info,", reduce delta_mband"
         CALL ipec_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     store hermitian matrices fg.
c-----------------------------------------------------------------------
      fvec=0
      gvec=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            fvec(iqty)=fmatb(1+ipert-jpert,jpert)
            gvec(iqty)=gmat(ipert,jpert)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     store non-hermitian matrix k.
c-----------------------------------------------------------------------
      kvec=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            kvec(iqty)=kmat(ipert,jpert)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy hermitian banded matrices fg.
c-----------------------------------------------------------------------
      fmats=0
      gmats=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            fmats(1+ipert-jpert,jpert)=fvec(iqty)
            gmats(1+ipert-jpert,jpert)=gvec(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy non-Hermitian banded matrix k.
c-----------------------------------------------------------------------
      kmats=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            kmats(1+mband+ipert-jpert,jpert)=kvec(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_matrix
c-----------------------------------------------------------------------
c     subprogram 6. idcon_vacuum.
c     read vacuum.bin from vacuum code.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_vacuum
c-----------------------------------------------------------------------
c     read vacuum data.
c-----------------------------------------------------------------------
      WRITE(*,*)"Reading vacuum energy matrices"
      CALL bin_open(bin_unit,ivacuumfile,"OLD","REWIND","none")
      READ(bin_unit)nths2,nfm2
      ALLOCATE(grri(nths2,nfm2))
      READ(bin_unit)grri
      READ(bin_unit)nths2,nfm2
      ALLOCATE(grre(nths2,nfm2))
      READ(bin_unit)grre
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_vacuum

      END MODULE idcon_mod
