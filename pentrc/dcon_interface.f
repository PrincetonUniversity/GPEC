!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

      MODULE dcon_interface
c-----------------------------------------------------------------------
c     !! SELECT ROUTINES FROM IDEAL PERTURBED EQUILIBRIUM CONTROL !!
c     Subroutines by Jong-Kyu Park for interfacing with DCON outputs.
c     For the majority, these are identical to counterparts in the
c     IPEC idcon module. There are a few additions from ipeq and ismath
c     modules.
c     Major changes are:
c     - psixy=1 IF statement that optionally read psi_in.bin removed (no ieqfile)
c     - No vacuum data read (no ivacuumfile)
c     - Major culling of global variables
c      -- Explicit imports
c      -- Local variables moved to their respective subroutines
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      USE spline_mod, only : spline_type,spline_eval,spline_alloc,
     $                       spline_dealloc,spline_fit,spline_int
      USE cspline_mod, only : cspline_type,cspline_eval,cspline_alloc,
     $                       cspline_dealloc,cspline_fit,cspline_int
      USE fspline_mod, only : fspline_type,fspline_eval,fspline_alloc,
     $                       fspline_dealloc,fspline_fit_1,fspline_fit_2
      USE bicube_mod, only : bicube_type,bicube_eval,bicube_alloc,
     $                       bicube_dealloc,bicube_fit
      USE params, only : r4,r8,pi,twopi,mu0
      USE utilities, only : get_free_file_unit,iscdftb,iscdftf
      
      IMPLICIT NONE

      LOGICAL :: edge_flag=.FALSE.,fft_flag=.FALSE.,
     $     kin_flag=.FALSE.,con_flag=.FALSE.,verbose=.TRUE.
      INTEGER :: mr,mz,mpsi,mstep,mpert,mband,mtheta,mthvac,mthsurf,
     $     mfix,mhigh,mlow,msing,nfm2,nths2,lmpert,lmlow,lmhigh,
     $     power_b,power_r,power_bp,
     $     nn,mmin,mmax

      REAL(r8) :: bo,ro,zo,psio,chi1,mthsurf0,psilow,psilim,qlim,
     $     singfac_min,   amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,betap1,betap2,betap3,betat,betan,bt0,
     $     q0,qmin,qmax,qa,crnt,q95,shotnum,shottime
      COMPLEX(r8), PARAMETER  :: ifac = (0,1)
      CHARACTER(16) :: jac_type,machine
      CHARACTER(128) :: idconfile

      LOGICAL, DIMENSION(:), POINTER :: sing_flag
      INTEGER, DIMENSION(:), POINTER :: fixstep,mfac,lmfac

      REAL(r8), DIMENSION(:), POINTER :: psifac,rhofac,qfac,singfac,
     $     r,z,theta,et,ep,ee,eft,efp

      COMPLEX(r8), DIMENSION(:), POINTER ::
     $     edge_mn,edge_fun
      COMPLEX(r8), DIMENSION(:,:), POINTER :: wt,wt0,wft,wtraw,
     $     amat,bmat,cmat,fmats,gmats,kmats

      TYPE(spline_type) :: sq,geom
      TYPE(bicube_type) :: psi_in,eqfun,rzphi
      TYPE(cspline_type) :: u1,u2,u3,u4
      TYPE(cspline_type) :: smats,tmats,xmats,ymats,zmats
      TYPE(fspline_type) :: metric

      TYPE :: resist_type
      REAL(r8) :: e,f,h,m,g,k,eta,rho,taua,taur,di,dr,sfac,deltac
      END TYPE resist_type

      TYPE :: solution_type
      INTEGER :: msol
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u
      END TYPE solution_type

      TYPE :: fixfac_type
      INTEGER :: msol
      INTEGER, DIMENSION(:), POINTER :: index
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fixfac,transform,gauss
      END TYPE fixfac_type

      TYPE :: sing_type
      INTEGER :: msol_l,msol_r,jfix,jpert
      REAL(r8) :: psifac,q,q1
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: ca_l,ca_r
      TYPE(resist_type) :: restype
      END TYPE sing_type

      TYPE(solution_type), DIMENSION(:), POINTER :: soltype
      TYPE(fixfac_type), DIMENSION(:), POINTER :: fixtype
      TYPE(sing_type), DIMENSION(:), POINTER :: singtype

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
      CHARACTER(2) :: sn
      INTEGER :: m,data_type,ifix,ios,msol,istep,ising,itheta,i,in_unit
      REAL(r8) :: sfac0

      REAL(r4), DIMENSION(:,:), POINTER :: rgarr,zgarr,psigarr
c-----------------------------------------------------------------------
c     open euler.bin and read header.
c-----------------------------------------------------------------------
      WRITE(*,*)"Reading dcon eigenfuctions"
      in_unit = get_free_file_unit(-1)
      OPEN(UNIT=in_unit,FILE=idconfile,STATUS="OLD",POSITION="REWIND",
     $        FORM="UNFORMATTED")
      READ(in_unit)mlow,mhigh,nn,mpsi,mtheta,ro,zo
      READ(in_unit)mband,mthsurf0,mthvac,psio,psilow,psilim,qlim,
     $     singfac_min
      READ(in_unit)power_b,power_r,power_bp
      READ(in_unit)kin_flag,con_flag
      READ(in_unit) amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,betap1,betap2,betap3,betat,betan,bt0,
     $     q0,qmin,qmax,qa,crnt,q95,shotnum,shottime
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
      chi1=twopi*psio
      mpert=mhigh-mlow+1
      lmlow=mmin
      lmhigh=mmax
      IF (mlow<mmin) lmlow=mlow
      IF (mhigh>mmax) lmhigh=mhigh
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
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')nn
      ENDIF
c-----------------------------------------------------------------------
c     only accept mband=0.
c-----------------------------------------------------------------------
      IF (mband /= (mpert-1)) THEN
         STOP 'PENTRC needs full band matrix'
      ENDIF
c-----------------------------------------------------------------------
c     read equilibrium on flux coordinates.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq,mpsi,4)
      CALL bicube_alloc(rzphi,mpsi,mtheta,4)

      rzphi%periodic(2)=.TRUE.
      READ(in_unit)sq%xs,sq%fs,sq%fs1,sq%xpower
      READ(in_unit)rzphi%xs,rzphi%ys,
     $     rzphi%fs,rzphi%fsx,rzphi%fsy,rzphi%fsxy,
     $     rzphi%x0,rzphi%y0,rzphi%xpower,rzphi%ypower
      mstep=-1
      mfix=0
      msing=0
c-----------------------------------------------------------------------
c     count solutions in euler.bin.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)"Counting and reading dcon solutions"
      DO
         READ(UNIT=in_unit,IOSTAT=ios)data_type
         IF(ios /= 0)EXIT
         SELECT CASE(data_type)
         CASE(1)
            mstep=mstep+1
            READ(UNIT=in_unit)
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
            READ(UNIT=in_unit)
         CASE(4)
            msing=msing+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(5)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE DEFAULT
            WRITE(message,'(a,i1,a,i4)')"Cannot recognize data_type = ",
     $           data_type,", at istep = ",istep
            PRINT *,message
            STOP
         END SELECT
      ENDDO
      CLOSE(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays and prepare to read data.
c-----------------------------------------------------------------------
      IF(verbose)THEN
        WRITE(*,*)"mlow = ",mlow,", mhigh = ",mhigh,", mpert = ",mpert
        WRITE(*,*)"mstep = ",mstep,", mfix = ",mfix,", msing = ",msing
      ENDIF
      ALLOCATE(psifac(0:mstep),rhofac(0:mstep),qfac(0:mstep),
     $     soltype(0:mstep),singtype(msing))
      ALLOCATE(fixstep(0:mfix+1),fixtype(0:mfix),sing_flag(mfix))
      ALLOCATE(et(mpert),ep(mpert),ee(mpert))
      ALLOCATE(wt(mpert,mpert),wt0(mpert,mpert))
      ALLOCATE(wft(mpert,mpert),eft(mpert),efp(mpert)) !LOGAN
      eft = -1 ! LOGAN - used to tell if it is read (actual vals >0)
      wft = 0
      fixstep(0)=0
      fixstep(mfix+1)=mstep
      in_unit = get_free_file_unit(-1)
      OPEN(UNIT=in_unit,FILE=idconfile,STATUS="OLD",POSITION="REWIND",
     $        FORM="UNFORMATTED")
      READ(in_unit)mlow,mhigh,nn
      READ(in_unit)mband,mthsurf0,mthvac,psio,psilow,psilim,qlim,
     $     singfac_min
      READ(in_unit)power_b,power_r,power_bp
      READ(in_unit)kin_flag,con_flag
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
            ALLOCATE(soltype(istep)%u(mpert,soltype(istep)%msol,4))
            READ(UNIT=in_unit)soltype(istep)%u(:,:,1:2)
            READ(UNIT=in_unit)soltype(istep)%u(:,:,3:4)
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
            READ(UNIT=in_unit)wt0
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
         CASE(5)
            ALLOCATE(wft(mpert,mpert),eft(mpert),efp(mpert))
            ALLOCATE(wtraw(mpert,mpert))
            READ(UNIT=in_unit)ep
            READ(UNIT=in_unit)et
            READ(UNIT=in_unit)wt
            READ(UNIT=in_unit)efp
            READ(UNIT=in_unit)eft
            READ(UNIT=in_unit)wft
            READ(UNIT=in_unit)wtraw
         END SELECT
      ENDDO
      IF (psifac(mstep)<psilim-(1e-4)) THEN
         WRITE(message,'(a)')"Terminated by zero crossing"
         STOP "Terminated by zero crossing"
      ENDIF
      rhofac=SQRT(psifac)
c-----------------------------------------------------------------------
c     normalize plasma/vacuum eigenvalues and eigenfunctions.
c-----------------------------------------------------------------------
      et=et/(mu0*2.0)*psio**2*(chi1*1e-3)**2
      ep=ep/(mu0*2.0)*psio**2*(chi1*1e-3)**2
      ee=et-ep
      wt=wt*(chi1*1e-3)
      wt0=wt0/(mu0*2.0)*psio**2
      IF(data_type==5)THEN
         wtraw=wtraw*(chi1*1e-3)
         eft=eft/(mu0*2.0)*psio**2*(chi1*1e-3)**2
         efp=efp/(mu0*2.0)*psio**2*(chi1*1e-3)**2
         wft=wft*(chi1*1e-3)
      ENDIF
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
      CLOSE(in_unit)
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
      CALL cspline_alloc(u3,mstep,mpert)
      CALL cspline_alloc(u4,mstep,mpert)
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
      SUBROUTINE idcon_build(egnum,xspmn)

      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ifix,jfix,kfix,ieq,info
      INTEGER, DIMENSION(mpert) :: ipiv
      COMPLEX(r8), DIMENSION(mpert) :: uedge,temp1
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp2
c-----------------------------------------------------------------------
c     construct uedge.
c-----------------------------------------------------------------------
      IF(edge_flag)THEN
         uedge=xspmn
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
            u1%fs(istep,:)=MATMUL(soltype(istep)%u(:,1:mpert,1),temp1)
            u2%fs(istep,:)=MATMUL(soltype(istep)%u(:,1:mpert,2),temp1)
            u3%fs(istep,:)=MATMUL(soltype(istep)%u(:,1:mpert,3),temp1)
            u4%fs(istep,:)=MATMUL(soltype(istep)%u(:,1:mpert,4),temp1)
         ENDDO         
         jfix=kfix+1
      ENDDO
c-----------------------------------------------------------------------
c     fit the functions.
c-----------------------------------------------------------------------
      u1%xs=psifac
      u2%xs=psifac
      u3%xs=psifac
      u4%xs=psifac
      CALL cspline_fit(u1,"extrap")
      CALL cspline_fit(u2,"extrap")
      CALL cspline_fit(u3,"extrap")
      CALL cspline_fit(u4,"extrap")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_build
c-----------------------------------------------------------------------
c     subprogram 4. idcon_metric.
c     reconstructs metric tensors from dcon.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_metric
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: ipsi,itheta,out_unit
      REAL(r8) :: w11,w12,delpsi,qintb,val1,val2,
     $      q,rfac,eta,jac,jac1
      
      REAL(r8), DIMENSION(0:mpsi) :: psitor
      REAL(r8), DIMENSION(0:mpsi,0:mtheta) :: rs,zs,jacb2
      REAL(r8), DIMENSION(3,3) :: w,v
      TYPE(spline_type) :: qs

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
      CALL bicube_alloc(eqfun,mpsi,mtheta,3)
      eqfun%xs=rzphi%xs
      eqfun%ys=rzphi%ys
      eqfun%name="equilibrium funs"
      eqfun%xtitle="psi"
      eqfun%ytitle="theta"
      eqfun%title=(/" modb  "," divx1 "," divx2 "/)      
c-----------------------------------------------------------------------
c     begin loop over nodes.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         q=sq%f(4)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(itheta/REAL(mtheta,r8)+rzphi%f(2))
            rs(ipsi,itheta)=ro+rfac*COS(eta)
            zs(ipsi,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
            w11=(1+rzphi%fy(2))*twopi**2*rfac*rs(ipsi,itheta)/jac
            w12=-rzphi%fy(1)*pi*rs(ipsi,itheta)/(rfac*jac)
            delpsi=SQRT(w11**2+w12**2)
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(1,3)=rzphi%fx(3)*rs(ipsi,itheta)
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(2,3)=rzphi%fy(3)*rs(ipsi,itheta)
            v(3,3)=twopi*rs(ipsi,itheta)
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
c-----------------------------------------------------------------------
            eqfun%fs(ipsi,itheta,1)=SQRT((chi1**2*delpsi**2+
     $           sq%f(1)**2)/(twopi*rs(ipsi,itheta))**2)
            eqfun%fs(ipsi,itheta,2)=(SUM(v(1,:)*v(2,:))+q*v(3,3)*v(1,3))
     $           /(jac*eqfun%fs(ipsi,itheta,1)**2)
            eqfun%fs(ipsi,itheta,3)=(v(2,3)*v(3,3)+q*v(3,3)*v(3,3))
     $           /(jac*eqfun%fs(ipsi,itheta,1)**2)
            jacb2(ipsi,itheta)=jac*eqfun%fs(ipsi,itheta,1)**2
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
c     write equilibrium quantities.
c-----------------------------------------------------------------------
      val1=0.0
      val2=1.0
      CALL spline_alloc(qs,mpsi+2,1)
      qs%xs(0)=val1
      qs%xs(1:mpsi+1)=sq%xs
      qs%xs(mpsi+2)=val2
      CALL spline_eval(sq,val1,0)
      qs%fs(0,1)=sq%f(4)
      qs%fs(1:mpsi+1,1)=sq%fs(:,4)
      CALL spline_eval(sq,val2,0)
      qs%fs(mpsi+2,1)=sq%f(4)
      CALL spline_fit(qs,"extrap") 
      CALL spline_int(qs)
      qintb=qs%fsi(mpsi+2,1)
      psitor(:)=qs%fsi(1:mpsi+1,1)/qintb
      CALL spline_dealloc(qs)

      out_unit = get_free_file_unit(-1)
      OPEN(UNIT=out_unit,FILE="idcon_equil.out",STATUS="UNKNOWN")
      WRITE(out_unit,*)"IDCON_EQUIL: "//
     $     "Various equilibrium quantities"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_type = ",jac_type
      WRITE(out_unit,'(2(1x,a8,1x,I6))')"mpsi =",mpsi,"mtheta =",mtheta
      WRITE(out_unit,'(2(1x,a14,1x,es16.8))')"psi_edge =",psio,
     $     "psitor_edge =",qintb*psio
      WRITE(out_unit,*)     
      WRITE(out_unit,*)" Flux functions:"
      WRITE(out_unit,*) 
      WRITE(out_unit,'(6(1x,a16))')"psi","psitor","p","q","g","I"
      DO ipsi=0,mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         CALL spline_alloc(qs,mtheta,1)
         qs%xs=rzphi%ys
         qs%fs(:,1)=jacb2(ipsi,:)
         CALL spline_fit(qs,"periodic") 
         CALL spline_int(qs)          
         WRITE(out_unit,'(6(1x,es16.8))')sq%xs(ipsi),psitor(ipsi),
     $        sq%f(2)/mu0,sq%f(4),sq%f(1)/(twopi*mu0),
     $        (qs%fsi(mtheta,1)/(twopi*chi1)-sq%f(4)*sq%f(1)/twopi)/mu0
         CALL spline_dealloc(qs)
      ENDDO
      WRITE(out_unit,*)     
      WRITE(out_unit,*)" 2D functions:"
      WRITE(out_unit,*)  
      WRITE(out_unit,'(8(1x,a16))')"psi","theta","r","z",
     $     "eta","dphi","jac","b0"
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),0)
            CALL bicube_eval(eqfun,eqfun%xs(ipsi),eqfun%ys(itheta),0)
            WRITE(out_unit,'(8(1x,es16.8))')rzphi%xs(ipsi),
     $           rzphi%ys(itheta)*twopi,rs(ipsi,itheta),zs(ipsi,itheta),
     $           rzphi%f(2)*twopi,rzphi%f(3)*twopi,
     $           rzphi%f(4)/(twopi*sq%f(4)*chi1),eqfun%f(1)
         ENDDO
      ENDDO
      CLOSE(out_unit)
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
      REAL(r8) :: q,q1,p,p1
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
         PRINT *,message
         STOP
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
         !CALL ipec_stop(message)
         PRINT *,message
         STOP         
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
c     subprogram 6. idcon_action_matrices.
c     Equilibrium matrices necessary to calc perturbed mod b for gpec.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_action_matrices()!(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: ipsi,istep,ipert,jpert,itheta,dm,m1,m2
      REAL(r8) :: psi,angle,rs,
     $     g12,g22,g13,g23,g33,singfac2,b2h,b2hp,b2ht,
     $     p1,q,rfac,eta,jac,jac1
      COMPLEX(r8), DIMENSION(-mband:mband) :: 
     $     sband,tband,xband,yband1,yband2,zband1,zband2,zband3
      COMPLEX(r8), DIMENSION(mpert,mpert) :: smat,tmat,xmat,ymat,zmat
      REAL(r8), DIMENSION(3,3) :: w,v

      TYPE(fspline_type) :: fmodb
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing perturbed b field for gpec"
c-----------------------------------------------------------------------
c     set up fourier-spline type.
c-----------------------------------------------------------------------
      CALL cspline_alloc(smats,mpsi,mpert**2)
      CALL cspline_alloc(tmats,mpsi,mpert**2)
      CALL cspline_alloc(xmats,mpsi,mpert**2)
      CALL cspline_alloc(ymats,mpsi,mpert**2)
      CALL cspline_alloc(zmats,mpsi,mpert**2)
      smats%xs=sq%xs
      tmats%xs=sq%xs
      xmats%xs=sq%xs
      ymats%xs=sq%xs
      zmats%xs=sq%xs

      CALL fspline_alloc(fmodb,mpsi,mtheta,mband,8)
      fmodb%xs=rzphi%xs
      fmodb%ys=rzphi%ys*twopi
      fmodb%name="fmodb"
      fmodb%xtitle=" psi  "
      fmodb%ytitle="theta "
      fmodb%title=(/" smat  "," tmat  "," xmat  ",
     $     " ymat1 "," ymat2 "," zmat1 ", " zmat2 "," zmat3 "/)
c-----------------------------------------------------------------------
c     computes fourier series of geometric tensors.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psi=sq%xs(ipsi)
         p1=sq%fs1(ipsi,2)
         q=sq%fs(ipsi,4)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            CALL bicube_eval(eqfun,psi,rzphi%ys(itheta),1)
            angle=rzphi%ys(itheta)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(angle+rzphi%f(2))
            rs=ro+rfac*cos(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
            b2h=eqfun%f(1)**2/2
            b2hp=eqfun%f(1)*eqfun%fx(1)
            b2ht=eqfun%f(1)*eqfun%fy(1)
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(1,3)=rzphi%fx(3)*rs
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(2,3)=rzphi%fy(3)*rs
            v(3,3)=twopi*rs

            g12=SUM(v(1,:)*v(2,:))
            g13=v(3,3)*v(1,3)
            g22=SUM(v(2,:)**2)
            g23=v(2,3)*v(3,3)
            g33=v(3,3)*v(3,3)
            
            fmodb%fs(ipsi,itheta,1)=jac*(p1+b2hp)
     $           -chi1**2*b2ht*(g12+q*g13)/(jac*b2h*2)
            fmodb%fs(ipsi,itheta,2)=
     $           chi1**2*b2ht*(g23+q*g33)/(jac*b2h*2)
            fmodb%fs(ipsi,itheta,3)=jac*b2h*2
            fmodb%fs(ipsi,itheta,4)=jac1*b2h*2-chi1**2*b2h*2*eqfun%fy(2)
            fmodb%fs(ipsi,itheta,5)=-twopi*chi1**2/jac*(g12+q*g13)
            fmodb%fs(ipsi,itheta,6)=chi1**2*b2h*2*eqfun%fy(3)
            fmodb%fs(ipsi,itheta,7)=twopi*chi1**2/jac*(g23+q*g33)
            fmodb%fs(ipsi,itheta,8)=twopi*chi1**2/jac*(g22+q*g23)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit Fourier-spline type.
c-----------------------------------------------------------------------
      IF(fft_flag)THEN
         CALL fspline_fit_2(fmodb,"extrap",.FALSE.)
      ELSE
         CALL fspline_fit_1(fmodb,"extrap",.FALSE.)
      ENDIF

      DO ipsi=0,mpsi
         
         q=sq%fs(ipsi,4)
         sband(0:-mband:-1)=fmodb%cs%fs(ipsi,1:mband+1)
         tband(0:-mband:-1)=fmodb%cs%fs(ipsi,mband+2:2*mband+2)
         xband(0:-mband:-1)=fmodb%cs%fs(ipsi,2*mband+3:3*mband+3)
         yband1(0:-mband:-1)=fmodb%cs%fs(ipsi,3*mband+4:4*mband+4)
         yband2(0:-mband:-1)=fmodb%cs%fs(ipsi,4*mband+5:5*mband+5)
         zband1(0:-mband:-1)=fmodb%cs%fs(ipsi,5*mband+6:6*mband+6)
         zband2(0:-mband:-1)=fmodb%cs%fs(ipsi,6*mband+7:7*mband+7)
         zband3(0:-mband:-1)=fmodb%cs%fs(ipsi,7*mband+8:8*mband+8)

         sband(1:mband)=CONJG(sband(-1:-mband:-1))
         tband(1:mband)=CONJG(tband(-1:-mband:-1))
         xband(1:mband)=CONJG(xband(-1:-mband:-1))
         yband1(1:mband)=CONJG(yband1(-1:-mband:-1))
         yband2(1:mband)=CONJG(yband2(-1:-mband:-1))
         zband1(1:mband)=CONJG(zband1(-1:-mband:-1))
         zband2(1:mband)=CONJG(zband2(-1:-mband:-1))
         zband3(1:mband)=CONJG(zband3(-1:-mband:-1))

         ipert=0
         DO m1=mlow,mhigh
            ipert=ipert+1
            DO dm=MAX(1-ipert,-mband),MIN(mpert-ipert,mband)
               m2=m1+dm
               singfac2=m2-nn*q
               jpert=ipert+dm
c-----------------------------------------------------------------------
c     construct primitive matrices.
c-----------------------------------------------------------------------
               smat(ipert,jpert)=sband(dm)
               tmat(ipert,jpert)=tband(dm)
               xmat(ipert,jpert)=xband(dm)
               ymat(ipert,jpert)=yband1(dm)+ifac*singfac2*yband2(dm)
               zmat(ipert,jpert)=zband1(dm)+
     $              ifac*(m2*zband2(dm)+nn*zband3(dm))
            ENDDO
         ENDDO

         smats%fs(ipsi,:)=RESHAPE(smat,(/mpert**2/))
         tmats%fs(ipsi,:)=RESHAPE(tmat,(/mpert**2/))
         xmats%fs(ipsi,:)=RESHAPE(xmat,(/mpert**2/))
         ymats%fs(ipsi,:)=RESHAPE(ymat,(/mpert**2/))
         zmats%fs(ipsi,:)=RESHAPE(zmat,(/mpert**2/))

      ENDDO

      ! global splines like sq,rzphi,eqfun etc.
      CALL cspline_fit(smats,"extrap")
      CALL cspline_fit(tmats,"extrap")
      CALL cspline_fit(xmats,"extrap")
      CALL cspline_fit(ymats,"extrap")
      CALL cspline_fit(zmats,"extrap")

      !CALL ipeq_dealloc
      CALL fspline_dealloc(fmodb)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_action_matrices
c-----------------------------------------------------------------------
c     compute root of discrete function by the secant method.
c     use when only function is monotonic.
c-----------------------------------------------------------------------
      FUNCTION issect(gridnum,xvec,yvec,yval)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: gridnum
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: xvec
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: yvec
      REAL(r8), INTENT(IN) :: yval

      REAL(r8) :: xval1,xval2,xval3,yval1,yval2,yval3,issect
      TYPE(spline_type) :: yfun

      CALL spline_alloc(yfun,gridnum,1)
      yfun%xs=xvec
      yfun%fs(:,1)=yvec
      CALL spline_fit(yfun,"extrap")

      xval1=xvec(0)
      yval1=yvec(0)
      xval2=xvec(gridnum)
      yval2=yvec(gridnum)

      xval3=yval
      CALL spline_eval(yfun,xval3,0)
      yval3=yfun%f(1)

      DO
         IF (ABS(yval3-yval) .LT. 1e-6) EXIT
         IF (yval3 .GT. yval) THEN
            xval2=xval3
            yval2=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ELSE
            xval1=xval3
            yval1=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ENDIF
         CALL spline_eval(yfun,xval3,0)
         yval3=yfun%f(1)
      ENDDO
      issect=xval3
      CALL spline_dealloc(yfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION issect
c-----------------------------------------------------------------------
c     subprogram 9. idcon_coords.
c     transform coordinates to dcon coordinates. 
c-----------------------------------------------------------------------
      SUBROUTINE idcon_coords(psi,ftnmn,amf,amp,ri,bpi,bi,rci,ti,ji)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: amp,ri,bpi,bi,rci,ti,ji
      REAL(r8), INTENT(IN) :: psi
      INTEGER, DIMENSION(amp), INTENT(IN) :: amf
      COMPLEX(r8), DIMENSION(amp), INTENT(INOUT) :: ftnmn

      INTEGER :: i,ising,itheta
      REAL(r8) :: thetai,jarea,
     $      rfac,eta,jac,bpfac,btfac,bfac,fac

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacfac
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun
      REAL(r8), DIMENSION(3,3) :: w

      TYPE(spline_type) :: spl     
      
      CALL spline_alloc(spl,mthsurf,2)
      spl%xs=theta

      CALL spline_eval(sq,psi,0)
      dphi=0
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         bpfac=psio*delpsi(itheta)/r(itheta)
         btfac=sq%f(1)/(twopi*r(itheta))
         bfac=SQRT(bpfac*bpfac+btfac*btfac)
         fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
         spl%fs(itheta,1)=fac/(r(itheta)**ri*rfac**rci)*
     $        bpfac**bpi*bfac**bi
         ! jacobian for coordinate angle at dcon angle
         spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**ri*rfac**rci/
     $        (bpfac**bpi*bfac**bi)  
         IF (ti .EQ. 0) THEN
            dphi(itheta)=rzphi%f(3)
         ENDIF
      ENDDO      

      CALL spline_fit(spl,"periodic")
      CALL spline_int(spl)

      ! coordinate angle at dcon angle
      thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
      IF (ji .EQ. 1) THEN
         DO itheta=0,mthsurf
            ! jacobian at coordinate angle
            thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
            CALL spline_eval(spl,thetai,0)
            jacfac(itheta)=spl%f(2)
         ENDDO
         jarea=0
         DO itheta=0,mthsurf-1
            jarea=jarea+jacfac(itheta)/mthsurf
         ENDDO
         CALL iscdftb(amf,amp,ftnfun,mthsurf,ftnmn)
         ftnfun=ftnfun/jacfac*jarea
         CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
      ENDIF         
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------      
      ! compute given function in dcon angle
      DO itheta=0,mthsurf
         ftnfun(itheta)=0
         DO i=1,amp
            ftnfun(itheta)=ftnfun(itheta)+
     $           ftnmn(i)*EXP(ifac*twopi*amf(i)*thetas(itheta))
         ENDDO
      ENDDO

      ! multiply toroidal factor for dcon angle
      IF (ti .EQ. 0) THEN
         ftnfun(:)=ftnfun(:)*EXP(-ifac*nn*dphi(:))
      ELSE
         ftnfun(:)=ftnfun(:)*
     $        EXP(-twopi*ifac*nn*sq%f(4)*(thetas(:)-theta(:)))
      ENDIF

      CALL spline_dealloc(spl)

      CALL iscdftf(amf,amp,ftnfun,mthsurf,ftnmn)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_coords

c-----------------------------------------------------------------------
c     subprogram 4. issurfint.
c     surface integration by simple method.
c-----------------------------------------------------------------------
      FUNCTION issurfint(func,fs,psi,wegt,ave)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: fs,wegt,ave
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), DIMENSION(0:fs), INTENT(IN) :: func

      LOGICAL  :: first = .TRUE. ! Only set at first call, implicitly saved
      INTEGER  :: itheta
      REAL(r8) :: issurfint
      REAL(r8) :: rfac,eta,jac,area
      REAL(r8), DIMENSION(1,2) :: w
      REAL(r8), DIMENSION(0:fs) :: z,thetas

      ! note we had to make arrays allocatable to be allowed to save
      INTEGER  :: fsave
      REAL(r8) :: psave
      REAL(r8), DIMENSION(:), ALLOCATABLE :: jacs,delpsi,r,a
      SAVE :: psave,fsave,jacs,delpsi,r,a

      issurfint=0
      area=0
      IF(first .OR. psi/=psave .OR. fs/=fsave)THEN
         psave = psi
         fsave = fs
         IF(.NOT.first) DEALLOCATE(jacs,delpsi,r,a)
         first = .FALSE.
         ALLOCATE(jacs(0:fs),delpsi(0:fs),r(0:fs),a(0:fs))
         thetas=(/(itheta,itheta=0,fs)/)/REAL(fs,r8)
         DO itheta=0,fs-1
            CALL bicube_eval(rzphi,psi,thetas(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(thetas(itheta)+rzphi%f(2))
            a(itheta)=rfac
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         ENDDO
      ENDIF

      IF (wegt==0) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==1) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           r(itheta)*jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==2) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/r(itheta)/fs
         ENDDO
      ELSE IF (wegt==3) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           a(itheta)*jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE
         STOP "ERROR: issurfint wegt must be in [0,1,2,3]"
      ENDIF

      IF (ave==1) THEN
         DO itheta=0,fs-1
            area=area+jacs(itheta)*delpsi(itheta)/fs
         ENDDO
         issurfint=issurfint/area
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfint




      !=======================================================================
      subroutine set_geom
      !-----------------------------------------------------------------------
      !*DESCRIPTION:
      !   Form a spline of additional geometric flux surface functions.
      !
      !*ARGUMENTS:
      !
      !-----------------------------------------------------------------------
        ! declare variables
        integer  :: ipsi
        real(r8) :: psifac
        real(r8), dimension(0:mthsurf) :: unitfun

        unitfun = 1
        call spline_alloc(geom,mpsi,3)
        geom%title=(/"area  ","<r>   ","<R>   "/)
        geom%xs = sq%xs
        do ipsi=0,mpsi
            psifac = geom%xs(ipsi)
            geom%fs(ipsi,1) = issurfint(unitfun,mthsurf,psifac,0,0)
            geom%fs(ipsi,2) = issurfint(unitfun,mthsurf,psifac,3,1)
            geom%fs(ipsi,3) = issurfint(unitfun,mthsurf,psifac,1,1)
        enddo
        call spline_fit(geom,"extrap")
        !call spline_int(geom) ! not necessary yet

      end subroutine set_geom


      !=======================================================================
      subroutine set_eq(set_eqfun,set_sq,set_rzphi,
     $              set_smats,set_tmats,set_xmats,set_ymats,set_zmats,
     $              set_chi1,set_ro,set_nn,set_jac_type,
     $              set_mlow,set_mhigh,set_mpert,set_mthsurf)
      !----------------------------------------------------------------------- 
      !*DESCRIPTION: 
      !   Set the dcon equilibrium global variables directly. For internal use
      !   within kinetic dcon.
      !
      !*ARGUMENTS:
      !   set_eqfun : bicube_type
      !   set_sq : spline_type
      !   set_rzphi : bicube_type
      !       Equilibrium matrix
      !   set_chi1 : real
      !       Derivative of poloidal flux on psifac
      !   set_ro : real
      !       Major radius
      !   set_nn : integer
      !       toroidal mode number
      !   set_jac_type : character
      !       DCON working jacobian
      !   set_mfac : integer array
      !       Poloidal modes
      !   set_psifac : real array
      !       psi grid
      !   set_mpert : integer
      !       Number of poloidal modes
      !   set_mstep : integer
      !       Number grid intervals in psi
      !   set_mthsurf : integer
      !       Number grid intervals in theta
      !
      !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        integer :: set_mlow,set_mhigh,set_mpert,set_nn,set_mthsurf
        real(r8) :: set_ro,set_chi1
        !real(r8), dimension(:) :: set_psifac
        character(*), intent(in) :: set_jac_type
        
        type(spline_type) :: set_sq
        type(bicube_type) :: set_eqfun,set_rzphi
        type(cspline_type) :: set_smats,set_tmats,set_xmats,
     $      set_ymats,set_zmats
        
        integer :: m
        
        ! directly transfer dcon global equilibrium variables to PENTRC
        eqfun   =set_eqfun
        sq      =set_sq
        rzphi   =set_rzphi
        
        ! needed to create w_i^T*w_j coefficient matices
        smats   =set_smats 
        tmats   =set_tmats
        xmats   =set_xmats
        ymats   =set_ymats
        zmats   =set_zmats
        !! only needed if using old ipec_o1 inputs
        !mstep   =set_mstep
        !allocate(psifac(0:mstep))
        !psifac(:) = set_psifac(:)
        
        chi1    =set_chi1
        ro      =set_ro
        nn      =set_nn
        jac_type=set_jac_type
        mpert   =set_mpert
        allocate(mfac(mpert))
        mfac    =(/(m,m=set_mlow,set_mhigh)/)
        mthsurf = set_mthsurf

        ! evaluate field on axis
        call spline_eval(sq,0.0_r8,0)
        bo = abs(sq%f(1))/(twopi*ro)
        
        ! set additional geometric spline
        call set_geom

      end subroutine set_eq
      
c-----------------------------------------------------------------------
c     subprogram idcon_harvest.
c     log inputs with harvest.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE idcon_harvest(hlog)
c-----------------------------------------------------------------------
c     no declarations. All module variables.
c-----------------------------------------------------------------------
      ! harvest variables
      include 'harvest_lib.inc77'
      integer :: ierr
      character(len=65507) :: hlog
      character, parameter :: nul = char(0)

      ! log inputs with harvest
      ! standard CODEDB records
      if(machine/='') ierr=set_harvest_payload_str(hlog,'MACHINE'//nul,
     $   trim(machine)//nul)
      if(shotnum>0)
     $   ierr=set_harvest_payload_int(hlog,'SHOT'//nul,int(shotnum))
      if(shottime>0)
     $   ierr=set_harvest_payload_int(hlog,'TIME'//nul,int(shottime))
      ! DCON specifc records
      ierr=set_harvest_payload_int(hlog,'mpsi'//nul,mpsi)
      ierr=set_harvest_payload_int(hlog,'mtheta'//nul,mtheta)
      ierr=set_harvest_payload_int(hlog,'mlow'//nul,mlow)
      ierr=set_harvest_payload_int(hlog,'mhigh'//nul,mhigh)
      ierr=set_harvest_payload_int(hlog,'mpert'//nul,mpert)
      ierr=set_harvest_payload_int(hlog,'mband'//nul,mband)
      ierr=set_harvest_payload_dbl(hlog,'psilow'//nul,psilow)
      ierr=set_harvest_payload_dbl(hlog,'amean'//nul,amean)
      ierr=set_harvest_payload_dbl(hlog,'rmean'//nul,rmean)
      ierr=set_harvest_payload_dbl(hlog,'aratio'//nul,aratio)
      ierr=set_harvest_payload_dbl(hlog,'kappa'//nul,kappa)
      ierr=set_harvest_payload_dbl(hlog,'delta1'//nul,delta1)
      ierr=set_harvest_payload_dbl(hlog,'delta2'//nul,delta2)
      ierr=set_harvest_payload_dbl(hlog,'li1'//nul,li1)
      ierr=set_harvest_payload_dbl(hlog,'li2'//nul,li2)
      ierr=set_harvest_payload_dbl(hlog,'li3'//nul,li3)
      ierr=set_harvest_payload_dbl(hlog,'ro'//nul,ro)
      ierr=set_harvest_payload_dbl(hlog,'zo'//nul,zo)
      ierr=set_harvest_payload_dbl(hlog,'psio'//nul,psio)
      ierr=set_harvest_payload_dbl(hlog,'betap1'//nul,betap1)
      ierr=set_harvest_payload_dbl(hlog,'betap2'//nul,betap2)
      ierr=set_harvest_payload_dbl(hlog,'betap3'//nul,betap3)
      ierr=set_harvest_payload_dbl(hlog,'betat'//nul,betat)
      ierr=set_harvest_payload_dbl(hlog,'betan'//nul,betan)
      ierr=set_harvest_payload_dbl(hlog,'bt0'//nul,bt0)
      ierr=set_harvest_payload_dbl(hlog,'q0'//nul,q0)
      ierr=set_harvest_payload_dbl(hlog,'qmin'//nul,qmin)
      ierr=set_harvest_payload_dbl(hlog,'qmax'//nul,qmax)
      ierr=set_harvest_payload_dbl(hlog,'qa'//nul,qa)
      ierr=set_harvest_payload_dbl(hlog,'crnt'//nul,crnt)
      ierr=set_harvest_payload_dbl(hlog,'q95'//nul,q95)
      ierr=set_harvest_payload_dbl(hlog,'betan'//nul,betan)
      ierr=set_harvest_payload_dbl_array(hlog,'et'//nul,et,mpert)
      ierr=set_harvest_payload_dbl_array(hlog,'ep'//nul,ep,mpert)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE idcon_harvest

      END MODULE dcon_interface
