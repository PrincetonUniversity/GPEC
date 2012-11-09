c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     calculate ntv torque
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipntv_mod
c      1. ipntv_perts
c      2. ipntv_ntv
c      3. ipntv_roots
c-----------------------------------------------------------------------
c     subprogram 0. ipntv_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipntv_mod
      USE ipresp_mod      

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipntv_perts.
c     compute perturbed mod b and divergance of xi perp for ntv calc.
c-----------------------------------------------------------------------
      SUBROUTINE ipntv_perts(egnum,xspmn,rout,bpout,bout,rcout,tout,
     $     npsi,ntheta)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) ::egnum,rout,bpout,bout,rcout,tout,npsi,ntheta
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,itheta,iindex
      REAL(r8) :: ileft,jac,jac1
      COMPLEX(r8) :: expdphi

      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,divxprp_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xsp1_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmz_fun,xvt_fun,xvz_fun,xmt_fun
      COMPLEX(r8), DIMENSION(mstep,mpert) :: eulbparmns,lagbparmns,
     $     divxprpmns,eulbparmns_out,lagbparmns_out,divxprpmns_out
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: eulbparfuns,lagbparfuns
     $     ,divxprpfuns,eulbparfuns_out,lagbparfuns_out,divxprpfuns_out

      REAL(r8), DIMENSION(:), POINTER :: psis,ches,chex,chetheta
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: chea,chepertmns,
     $     chepertfuns

      TYPE(cspline_type) :: chespl, cspl
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing delta-B and div normal-displacement for NTV"

      CALL bicube_alloc(pertfuns,mstep-1,mthsurf,4)
      pertfuns%title = (/"Rpmodb","Ipmodb","Rdivxn","Idivxn"/)
      pertfuns%xs(0:mstep-1) = psifac(1:mstep)
      pertfuns%ys(:) = theta(:)
      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta
      CALL idcon_build(egnum,xspmn)
      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a28)')
     $        "volume = ",iindex,"% perturbation computations"
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))         
         CALL ipeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b and div x variations in hamada.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         CALL iscdftb(mfac,mpert,xmz_fun,mthsurf,xmz_mn)
         CALL iscdftb(mfac,mpert,xmt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn) 

         singfac=mfac-nn*sq%f(4)
         xsp1_mn=xsp1_mn*(singfac**2/(singfac**2+reg_spot**2))
         CALL iscdftb(mfac,mpert,xsp1_fun,mthsurf,xsp1_mn)

         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),0)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),0)
            jac=rzphi%f(4)
            cspl%fs(itheta,1)=xmt_fun(itheta)-
     $           (chi1/eqfun%f(1))**2/jac*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))
            cspl%fs(itheta,2)=xmz_fun(itheta)-
     $           sq%f(4)*(chi1/eqfun%f(1))**2/jac*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))
         ENDDO
         CALL cspline_fit(cspl,"periodic")
 
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            CALL cspline_eval(cspl,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)
            eulbparfuns(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(jac*eqfun%f(1))
            lagbparfuns(istep,itheta)=
     $           eulbparfuns(istep,itheta)+
     $           eqfun%fx(1)*xsp_fun(itheta)+eqfun%fy(1)*
     $           (xms_fun(itheta)/(chi1*sq%f(4))+
     $           xmz_fun(itheta)/(jac*sq%f(4))-
     $           (chi1/(jac*eqfun%f(1)))**2*
     $           (sq%f(4)*xvz_fun(itheta)+xvt_fun(itheta)))
            divxprpfuns(istep,itheta)=xsp1_fun(itheta)+
     $           (jac1/jac)*xsp_fun(itheta)+
     $           cspl%f1(1)/jac-(twopi*ifac*nn)*cspl%f(2)/jac
            dphi(itheta)=rzphi%f(3)
         ENDDO
         CALL iscdftf(mfac,mpert,eulbparfuns(istep,:),
     $        mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbparfuns(istep,:),
     $        mthsurf,lagbpar_mn)
         CALL iscdftf(mfac,mpert,divxprpfuns(istep,:),
     $        mthsurf,divxprp_mn)
         eulbparmns(istep,:)=eulbpar_mn
         lagbparmns(istep,:)=lagbpar_mn
         divxprpmns(istep,:)=divxprp_mn
c-----------------------------------------------------------------------
c     decompose components on the given coordinates.
c-----------------------------------------------------------------------
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL ipeq_bcoords(psifac(istep),eulbpar_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),lagbpar_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),divxprp_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
         ENDIF
         eulbparmns_out(istep,:)=eulbpar_mn
         lagbparmns_out(istep,:)=lagbpar_mn
         divxprpmns_out(istep,:)=divxprp_mn
         eulbparfuns_out(istep,:)=eulbparfuns(istep,:)*EXP(ifac*nn*dphi)
         lagbparfuns_out(istep,:)=lagbparfuns(istep,:)*EXP(ifac*nn*dphi)
         divxprpfuns_out(istep,:)=divxprpfuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     globally accessible splines in hamada.
c-----------------------------------------------------------------------
      pertfuns%fs(0:mstep-1,:,1) = REAL(lagbparfuns(1:mstep,:))
      pertfuns%fs(0:mstep-1,:,2) = AIMAG(lagbparfuns(1:mstep,:))
      pertfuns%fs(0:mstep-1,:,3) = REAL(divxprpfuns(1:mstep,:))
      pertfuns%fs(0:mstep-1,:,4) = AIMAG(divxprpfuns(1:mstep,:))
      CALL bicube_fit(pertfuns,"extrap","periodic")


      CALL ascii_open(out_unit,"ipec_pmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_PMODB: "//
     $     "Components in perturbed mod b"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     

      WRITE(out_unit,'(1x,a16,1x,a4,6(1x,a16))')"psi","m",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $     "real(divxprp)","imag(divxprp)"
      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es16.8,1x,I4,6(1x,es16.8))')
     $           psifac(istep),mfac(ipert),
     $           REAL(eulbparmns_out(istep,ipert)),
     $           AIMAG(eulbparmns_out(istep,ipert)),
     $           REAL(lagbparmns(istep,ipert)),
     $           AIMAG(lagbparmns(istep,ipert)),
     $           REAL(divxprpmns(istep,ipert)),
     $           AIMAG(divxprpmns(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_pmodb_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_FUN: "//
     $        "Components in perturbed mod b in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(10(1x,a16))')"r","z","psi","theta",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $        "real(divxprp)","imag(divxprp)"
         DO istep=1,mstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(10(1x,es16.8))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              psifac(istep),theta(itheta),
     $              REAL(eulbparfuns_out(istep,itheta)),
     $              AIMAG(eulbparfuns_out(istep,itheta)),
     $              REAL(lagbparfuns_out(istep,itheta)),
     $              AIMAG(lagbparfuns_out(istep,itheta)),
     $              REAL(divxprpfuns_out(istep,itheta)),
     $              AIMAG(divxprpfuns_out(istep,itheta))
               IF(pertfuns%fs(istep-1,itheta,1) .NE. 
     $              REAL(lagbparfuns(istep,itheta)) .OR. 
     $              pertfuns%fs(istep-1,itheta,4) .NE. 
     $              AIMAG(divxprpfuns(istep,itheta)))
     $              CALL ipec_stop("Error in bicubic fit.")
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      
      
      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "pmodb.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
     $              REAL(REAL(eulbparmns_out(istep,ipert)),4),
     $              REAL(AIMAG(eulbparmns_out(istep,ipert)),4),
     $              REAL(REAL(lagbparmns_out(istep,ipert)),4),
     $              REAL(AIMAG(lagbparmns_out(istep,ipert)),4)     
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"pmodb_2d.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf-1
         WRITE(bin_2d_unit)REAL(rs(9:mstep,1:mthsurf),4),
     $        REAL(zs(9:mstep,1:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(eulbparfuns_out(
     $        9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(eulbparfuns_out(
     $        9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(lagbparfuns_out(
     $        9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(lagbparfuns_out(
     $        9:mstep,1:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF



c-----------------------------------------------------------------------
c     Calculate Chebyshev polynomials in hamada.
c-----------------------------------------------------------------------
      IF (chebyshev_flag) THEN
         ALLOCATE(chea(mpert,0:nche,2),chex(0:mstep))
         CALL cspline_alloc(chespl,mstep,2)
         chespl%xs=psifac
         chex=2.0*(psifac-0.5)
         DO ipert=1,mpert
            DO i=0,nche
               chespl%fs(:,1)=cos(REAL(i,r8)*acos(chex))/
     $              sqrt(1.0-chex**2.0)*lagbparmns(:,ipert)*4.0/pi
               chespl%fs(:,2)=cos(REAL(i,r8)*acos(chex))/
     $              sqrt(1.0-chex**2.0)*divxprpmns(:,ipert)*4.0/pi
               CALL cspline_fit(chespl,"extrap")
               CALL cspline_int(chespl)
               chea(ipert,i,:)=chespl%fsi(mstep,:)
            ENDDO
         ENDDO
         chea(:,0,:)=chea(:,0,:)/2.0
         CALL cspline_dealloc(chespl)

         ALLOCATE(psis(npsi),ches(npsi))
         ALLOCATE(chepertmns(npsi,mpert,2))
         chepertmns=0
         psis=(/(istep,istep=0,npsi-1)/)/REAL(npsi-1,r8)*
     $        (psilim-psilow)+psilow
         ches=2.0*(psis-0.5)
         
         DO ipert=1,mpert
            DO i=0,nche
               chepertmns(:,ipert,1)=chepertmns(:,ipert,1)+
     $              chea(ipert,i,1)*cos(REAL(i,r8)*acos(ches))
               chepertmns(:,ipert,2)=chepertmns(:,ipert,2)+
     $              chea(ipert,i,2)*cos(REAL(i,r8)*acos(ches))
            ENDDO
         ENDDO

         ALLOCATE(chetheta(0:ntheta))
         ALLOCATE(chepertfuns(npsi,0:ntheta,2))
         chepertfuns=0
         chetheta=(/(itheta,itheta=0,ntheta)/)/REAL(ntheta,r8)
         DO itheta=0,ntheta
            DO ipert=1,mpert
               chepertfuns(:,itheta,:)=chepertfuns(:,itheta,:)+
     $              chepertmns(:,ipert,:)*EXP(twopi*ifac*mfac(ipert)
     $              *chetheta(itheta))
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     globally accessible chebyshev splines in hamada.
c-----------------------------------------------------------------------
         CALL bicube_alloc(pertchebs,npsi-1,ntheta,4)
         pertchebs%title = (/"Rpmodb","Ipmodb","Rdivxn","Idivxn"/)
         pertchebs%xs(0:npsi-1) = psis(1:npsi)
         pertchebs%ys = chetheta
         pertchebs%fs(0:npsi-1,:,1) = REAL(chepertfuns(1:npsi,:,1))
         pertchebs%fs(0:npsi-1,:,2) = AIMAG(chepertfuns(1:npsi,:,1))
         pertchebs%fs(0:npsi-1,:,3) = REAL(chepertfuns(1:npsi,:,2))
         pertchebs%fs(0:npsi-1,:,4) = AIMAG(chepertfuns(1:npsi,:,2))
         CALL bicube_fit(pertchebs,"extrap","periodic")           

c-----------------------------------------------------------------------
c     write chebyshev splines in output coordinates.
c-----------------------------------------------------------------------
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            DO istep=1,npsi
               CALL ipeq_bcoords(psis(istep),chepertmns(istep,:,1),
     $              mfac,mpert,rout,bpout,bout,rcout,tout,0)
               CALL ipeq_bcoords(psis(istep),chepertmns(istep,:,2),
     $              mfac,mpert,rout,bpout,bout,rcout,tout,0)
            ENDDO
         ENDIF
         CALL ascii_open(out_unit,"ipec_pmodb_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_CHEBYSHEV: "//
     $        "Reconstructed perturbed mod b by chebyshev"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "cstep =",npsi,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(1x,a16,1x,a4,4(1x,a16))')"psi","m",
     $        "real(lagb)","imag(lagb)","real(divx)","imag(divx)"
         DO istep=1,npsi
            DO ipert=1,mpert
               WRITE(out_unit,'(1x,es16.8,1x,I4,4(1x,es16.8))')
     $              psis(istep),mfac(ipert),
     $              REAL(chepertmns(istep,ipert,1)),
     $              AIMAG(chepertmns(istep,ipert,1)),
     $              REAL(chepertmns(istep,ipert,2)),
     $              AIMAG(chepertmns(istep,ipert,2))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (fun_flag) THEN
            CALL ascii_open(out_unit,"ipec_pmodb_chebyshev_fun_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_PMODB_CHEBYSHEV_FUN: "//
     $           "Components in perturbed mod b in chebyshev functions"
            WRITE(out_unit,*)     
            WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
            WRITE(out_unit,'(1x,2(a12,I4))')
     $           "npsi =",npsi,"ntheta =",ntheta
            WRITE(out_unit,*)     
            
            WRITE(out_unit,'(8(1x,a16))')"r","z","psi","theta",
     $           "real(lagb)","imag(lagb)","real(divx)","imag(divx)"
            DO istep=1,npsi
               DO itheta=0,ntheta
                  CALL bicube_eval(rzphi,psis(istep),chetheta(itheta),0)
                  rfac=SQRT(rzphi%f(1))
                  eta=twopi*(theta(itheta)+rzphi%f(2))
                  expdphi = EXP(ifac*nn*rzphi%f(3))
                  WRITE(out_unit,'(8(1x,es16.8))')
     $                 ro+rfac*COS(eta),zo+rfac*SIN(eta),
     $                 psis(istep),chetheta(itheta),
     $                 REAL(chepertfuns(istep,itheta,1)*expdphi),
     $                 AIMAG(chepertfuns(istep,itheta,1)*expdphi),
     $                 REAL(chepertfuns(istep,itheta,2)*expdphi),
     $                 AIMAG(chepertfuns(istep,itheta,2)*expdphi)
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipntv_perts
c-----------------------------------------------------------------------
c     subprogram 2. ipntv_ntv.
c     compute ntv torque as function of psi.
c-----------------------------------------------------------------------
      SUBROUTINE ipntv_ntv(kfile,npsi,ntheta,nx,nlmda,nl,
     $     imass,icharge,zimp,zmass,passing_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: kfile
      INTEGER, INTENT(IN) :: npsi,ntheta,nx,nlmda,nl,
     $     imass,icharge,zimp,zmass
      LOGICAL, INTENT(IN) :: passing_flag

      INTEGER :: istep,itheta,iindex,ilmda,ll,zeroloc,i,nlines,ibmin
      INTEGER, DIMENSION(2,1) :: jroot
      REAL(r8) :: ileft,epsr,inu,enu,magb,jac,g21,g21y,g22x,epst,
     $     aconst,bconst,wbnce,wgrad,wcurv,dj2,dt
      REAL(r8), DIMENSION(4) :: pert,turn_int
      REAL(r8), DIMENSION(2,1) :: xroot
      REAL(r8), DIMENSION(0:npsi) :: psis
      REAL(r8), DIMENSION(0:ntheta) :: thetas,bncthetas,intthetas,htheta
      REAL(r8), DIMENSION(0:eqfun%my) :: btest     
      REAL(r8), DIMENSION(-1:nlmda+1) :: lmdas,lmda
      REAL(r8), DIMENSION(0:nx) :: x,resfac
      REAL(r8), DIMENSION(:), POINTER :: zeff,zpitch
      REAL(r8), DIMENSION(0:npsi,0:nlmda,0:nx/10) :: checks
      REAL(r8), DIMENSION(0:npsi,0:nlmda,5) :: bnc_checks
      COMPLEX(r8), DIMENSION(0:ntheta) :: pl,djdt
      COMPLEX(r8), DIMENSION(0:nx) :: ntvdlde

      LOGICAL :: trapped

      TYPE(spline_type) :: kindat,spl_vpar, spl_wbdt, spl_djdt,
     $     spl_intvdlde, spl_entvdlde, spl_intvdl, spl_entvdl,
     $     spl_idwdlde, spl_edwdlde, spl_idwdl, spl_edwdl

      WRITE(*,*)"NTV torque calculation v1.0"

c-----------------------------------------------------------------------
c     read kinetic inputs
c-----------------------------------------------------------------------
      WRITE(*,*)"Reading kinetic input file"

      nlines=0
      CALL ascii_open(in_unit,kfile,"OLD")
      READ(in_unit,*) 
      DO
         READ(in_unit,*,END=999)
         nlines=nlines+1
      ENDDO
 999  REWIND(in_unit)

      CALL spline_alloc(kindat,nlines-1,8)
      kindat%title = (/"idens ","edens ","itemp ","etemp ","wexb  ",
     $     "loglam","  nuii","  nuei"/)
      READ(in_unit,*)
      DO istep=0,nlines-1
         READ(in_unit,*) kindat%xs(istep),kindat%fs(istep,1:5)
      ENDDO
      CALL ascii_close(in_unit)
      WRITE(*,*)"Have kinetic profiles from psi ",kindat%xs(0)," to ",
     $     kindat%xs(kindat%mx)

      ALLOCATE(zeff(0:nlines-1),zpitch(0:nlines-1))
      zeff = zimp-(kindat%fs(:,1)/kindat%fs(:,2))*icharge
     $     *(zimp-icharge)
      zpitch = 1.0+(1.0+zmass)/(2.0*zmass)*zimp
     $     *(zeff-1.0)/(zimp-zeff)
c-----------------------------------------------------------------------
c     in SI units.
c-----------------------------------------------------------------------
      kindat%fs(:,3:4) = kindat%fs(:,3:4)*1.602E-19
      kindat%fs(:,6) = 17.3-0.5*LOG(kindat%fs(:,2)/1.0e20)
     $     +1.5*LOG(kindat%fs(:,4)/1.602e-16)
      kindat%fs(:,7) = (zpitch/3.5E17)*kindat%fs(:,1)*kindat%fs(:,6)
     $     /(SQRT(1.0*imass)*(kindat%fs(:,3)/1.602e-16)**1.5)
      kindat%fs(:,8) = (zpitch/3.5E17)*kindat%fs(:,2)*kindat%fs(:,6)
     $     /(SQRT(emass/pmass)*(kindat%fs(:,4)/1.602e-16)**1.5)
      CALL spline_fit(kindat,"extrap")
      

c-----------------------------------------------------------------------
c     basic grids and splines for NTV
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing NTV torque"

      thetas=(/(itheta,itheta=1,ntheta+1)/)/REAL(ntheta+2,r8)
      lmdas = (/(ilmda,ilmda=0,nlmda+2)/)/REAL(nlmda+2,r8)
      psis=(/(istep,istep=0,npsi)/)/REAL(npsi,r8)*
     $     (psilim-psilow)+psilow
      x=(/(20.0*ilmda,ilmda=0,nx)/)/REAL(nx,r8)

      CALL spline_alloc(intv,npsi,2*nl+1)
      intv%xs = psis
      CALL spline_alloc(entv,npsi,2*nl+1)
      entv%xs = psis
      CALL spline_alloc(idw,npsi,2*nl+1)
      idw%xs = psis
      CALL spline_alloc(edw,npsi,2*nl+1)
      edw%xs = psis

      CALL spline_alloc(spl_intvdl,nlmda,2*nl+1)
      CALL spline_alloc(spl_entvdl,nlmda,2*nl+1)
      CALL spline_alloc(spl_idwdl,nlmda,2*nl+1)
      CALL spline_alloc(spl_edwdl,nlmda,2*nl+1)

      CALL spline_alloc(spl_intvdlde,nx,2*nl+1)
      spl_intvdlde%xs = x
      CALL spline_alloc(spl_entvdlde,nx,2*nl+1)
      spl_entvdlde%xs = x
      CALL spline_alloc(spl_idwdlde,nx,2*nl+1)
      spl_idwdlde%xs = x
      CALL spline_alloc(spl_edwdlde,nx,2*nl+1)
      spl_edwdlde%xs = x
  
      CALL spline_alloc(spl_wbdt,ntheta,3)
      spl_wbdt%title = (/"wbnce ","wgrad ","wcurv "/)
      CALL spline_alloc(spl_djdt,ntheta,4)
      spl_djdt%title = (/"Redj  ","Imdj  ","Redjpl","Imdjpl"/)

      DO istep=0,npsi
         iindex = FLOOR(REAL(istep,8)/FLOOR(SIZE(psis)/10.0))*10
         ileft = REAL(istep,8)/FLOOR(SIZE(psis)/10.0)*10-iindex
         IF (ileft == 0 .AND. istep /= 0)
     $        WRITE(*,'(1x,a9,i3,a18)')
     $        "volume = ",iindex,"% ntv computations"
c-----------------------------------------------------------------------
c     functions on magnetic surfaces
c-----------------------------------------------------------------------
         CALL spline_eval(kindat,psis(istep),1)
         CALL bicube_eval(rzphi,psis(istep),0.0_r8,0)
         CALL spline_eval(sq,psis(istep),0)
         inu = kindat%f(7)*0.5*ro/SQRT(rzphi%f(1))
         enu = kindat%f(8)*0.5*ro/SQRT(rzphi%f(1))
c-----------------------------------------------------------------------
c     create lambda grid concentrated near trapped/passing boundry
c-----------------------------------------------------------------------
         DO itheta=0,eqfun%my
            CALL bicube_eval(eqfun,psis(istep),eqfun%ys(itheta),0)
            btest(itheta) = eqfun%f(1)
         ENDDO
         IF(passing_flag) THEN
            WHERE(lmdas < 1.0/3.0)
               lmda=3.0*lmdas/MAXVAL(btest,1)
            ELSEWHERE
               lmda = (1.5*lmdas-0.5)*(1.0/MINVAL(btest,1)
     $           -1.0/MAXVAL(btest,1))+1.0/MAXVAL(btest,1)
            END WHERE
         ELSE
            lmda=lmdas(:)*(1.0/MINVAL(btest,1)-1.0/MAXVAL(btest,1))
     $           +1.0/MAXVAL(btest,1)
         ENDIF
         spl_intvdl%xs(0:nlmda) = lmda(0:nlmda)
         spl_entvdl%xs(0:nlmda) = lmda(0:nlmda)
         spl_idwdl%xs(0:nlmda) = lmda(0:nlmda)
         spl_edwdl%xs(0:nlmda) = lmda(0:nlmda)
         DO ilmda=0,nlmda
c-----------------------------------------------------------------------
c     regrid theta to banana range if trapped
c-----------------------------------------------------------------------
            trapped = 1.0-lmda(ilmda)*MAXVAL(btest) .LT. 0.0
            IF (trapped) THEN
               jroot(:,1)=(/ 0,eqfun%my /)
               xroot(:,1)=(/ 0.0,1.0 /)
               CALL spline_alloc(spl_vpar,eqfun%my,1)
               spl_vpar%xs =eqfun%ys
               spl_vpar%fs(:,1) = 1.0-lmda(ilmda)*btest
               CALL spline_fit(spl_vpar,"periodic")
               CALL ipntv_roots(spl_vpar,jroot,xroot)
               IF(spl_vpar%fs(jroot(1,1),1) .LE. 0.0) 
     $              xroot(:,1) = (/ xroot(2,1)-1.0, xroot(1,1) /)
               intthetas = thetas*(xroot(2,1)-xroot(1,1))+xroot(1,1)
               CALL spline_dealloc(spl_vpar)
            ELSE
               intthetas = (thetas-thetas(0))/thetas(ntheta)               
            ENDIF
            bncthetas = intthetas
            WHERE(intthetas .LT. 0.0) bncthetas = intthetas+1.0

            zeroloc = MINLOC(ABS(intthetas),DIM=1)-1
            spl_wbdt%xs(0:ntheta)=intthetas(0:ntheta)
            spl_djdt%xs(0:ntheta)=intthetas(0:ntheta)
            DO itheta=0,ntheta
c-----------------------------------------------------------------------
c     compute bounce and precession frequencies, and variation in action
c-----------------------------------------------------------------------
               CALL bicube_eval(eqfun,psis(istep),bncthetas(itheta),1)
               magb = eqfun%f(1)
               CALL fspline_eval(metric,psis(istep),
     $              twopi*bncthetas(itheta),1)
               jac = metric%f(7)
               g22x = metric%fx(2)
               g21y = twopi*metric%fy(6)
               g21 = metric%f(6)
               IF(chebyshev_flag) THEN
                  CALL bicube_eval(pertchebs,
     $                 psis(istep),bncthetas(itheta),0)
                  pert = pertchebs%f(:)
               ELSE
                  CALL bicube_eval(pertfuns,psis(istep),
     $                 bncthetas(itheta),0)
                  pert = pertfuns%f(:)
               ENDIF
               spl_wbdt%fs(itheta,1) = jac*magb
     $              /(1.0-lmda(ilmda)*magb)**0.5 
               spl_wbdt%fs(itheta,2) = (jac*magb**2.0*eqfun%fx(1)
     $              -g21*chi1**2.0*eqfun%fy(1)) 
     $              / (chi1*magb*(1.0-lmda(ilmda)*magb)**0.5)
               spl_wbdt%fs(itheta,3) = (1.0-lmda(ilmda)*magb)**0.5
     $              *( (jac/chi1)*eqfun%fx(1) - (chi1/magb**2.0)*g21
     $              *eqfun%fy(1) + (chi1/magb)*(g21y-g22x) )
               spl_djdt%fs(itheta,1) = spl_wbdt%fs(itheta,1)
     $              *(1.0-1.5*lmda(ilmda)*magb)*pert(1)/magb
     $              +jac*magb*(1.0-lmda(ilmda)*magb)**0.5*pert(3)
               spl_djdt%fs(itheta,2) = spl_wbdt%fs(itheta,1)
     $              *(1.0-1.5*lmda(ilmda)*magb)*pert(2)/magb
     $              +jac*magb*(1.0-lmda(ilmda)*magb)**0.5*pert(4)
            ENDDO
c-----------------------------------------------------------------------
c     bounce integration
c-----------------------------------------------------------------------
            IF (trapped) THEN
               htheta = 0.5
               spl_wbdt%fs(:,:) = 2.0*spl_wbdt%fs(:,:)
               CALL spline_fit(spl_wbdt,"extrap")
               dt = intthetas(1)-intthetas(0)
               Call spline_eval(spl_wbdt,intthetas(0)-dt,0)
               turn_int(1:3) = 0.5*dt*(spl_wbdt%f(:)+spl_wbdt%fs(0,:)
     $              +spl_wbdt%fs(ntheta,:))
               CAll spline_eval(spl_wbdt,intthetas(ntheta)+dt,0)
               turn_int(1:3) = turn_int(1:3)+0.5*dt*(spl_wbdt%f(:))
            ELSE
               htheta = 1.0
               CALL spline_fit(spl_wbdt,"periodic")
               turn_int = (/ 0.0, 0.0, 0.0, 0.0 /)
            ENDIF
            CALL spline_int(spl_wbdt)
            wbnce =twopi*chi1/(spl_wbdt%fsi(ntheta,1)+turn_int(1))
            wgrad =wbnce/chi1*lmda(ilmda)*(spl_wbdt%fsi(ntheta,2)
     $           +turn_int(2))
            wcurv =twopi*wbnce/chi1*(spl_wbdt%fsi(ntheta,3)+turn_int(3))
            IF(wbnce .NE. wbnce) CALL ipec_stop("Bounce integration "//
     $           "failure. Try decreasing ntheta.")

            htheta = htheta*(lmda(ilmda)*(spl_wbdt%fsi(:,2)-
     $           spl_wbdt%fsi(zeroloc,2))+2.0*(spl_wbdt%fsi(:,3)-
     $           spl_wbdt%fsi(zeroloc,3)))/(lmda(ilmda)*
     $           spl_wbdt%fsi(ntheta,2)+2.0*spl_wbdt%fsi(ntheta,3))
            DO ll=-nl,nl
               pl = EXP(twopi*ifac*ll*htheta)
               IF(trapped) THEN
                  spl_djdt%fs(:,3) = spl_djdt%fs(:,1)*2.0*REAL(pl)/chi1
                  spl_djdt%fs(:,4) = spl_djdt%fs(:,2)*2.0*REAL(pl)/chi1
                  CALL spline_fit(spl_djdt,"extrap")
                  CALL spline_eval(spl_djdt,intthetas(0)-dt,0)
                  turn_int=0.5*dt*(spl_djdt%f(:)+spl_djdt%fs(0,:)
     $                 +spl_djdt%fs(ntheta,:))
                  CALL spline_eval(spl_djdt,intthetas(ntheta)+dt,0)
                  turn_int=turn_int+0.5*dt*spl_djdt%f(:)
               ELSE
                  djdt = CMPLX(spl_djdt%fs(:,1),spl_djdt%fs(:,2))
                  spl_djdt%fs(:,3) = REAL(djdt*pl)/chi1
                  spl_djdt%fs(:,4) = AIMAG(djdt*pl)/chi1
                  CALL spline_fit(spl_djdt,"periodic")
                  turn_int=(/ 0.0, 0.0, 0.0, 0.0 /)
               ENDIF
               CALL spline_int(spl_djdt)
               dj2 = (spl_djdt%fsi(ntheta,3)+turn_int(3))**2.0
     $              +(spl_djdt%fsi(ntheta,4)+turn_int(4))**2.0
               IF(ll .EQ. 0) bnc_checks(istep,ilmda,:) = 
     $              (/ lmda(ilmda),wbnce,wgrad,wcurv,dj2 /)
c-----------------------------------------------------------------------
c     ntv integrand - energy integration
c-----------------------------------------------------------------------
               resfac = ll*wbnce*(2.0*x*kindat%f(3)/(imass*pmass))**0.5
     $              -nn*x*kindat%f(3)*(wgrad+wcurv)/(-echarge*icharge)
     $              -nn*kindat%f(5)
               ntvdlde = -wbnce*dj2*x**2.5*EXP(-x)/
     $              (ifac*resfac-(1.0+(0.5*ll)**2.0)*inu*x**(-1.5))
               spl_intvdlde%fs(:,ll+nl+1) = REAL(ntvdlde)
               spl_idwdlde%fs(:,ll+nl+1) = REAL(ntvdlde/(2.0*ifac*nn))

               resfac = ll*wbnce*(2.0*x*kindat%f(4)/emass)**0.5-nn*x
     $              *kindat%f(4)*(wgrad+wcurv)/echarge-nn*kindat%f(5)
               ntvdlde = -wbnce*dj2*x**2.5*EXP(-x)/
     $              (ifac*resfac-(1.0+(0.5*ll)**2.0)*enu*x**(-1.5))
               spl_entvdlde%fs(:,ll+nl+1) = REAL(ntvdlde)
               spl_edwdlde%fs(:,ll+nl+1) = REAL(ntvdlde/(2.0*ifac*nn))
            ENDDO

            DO i=0,nx/10
               checks(istep,ilmda,i)=spl_intvdlde%fs(i*10,nl+1)
            ENDDO
            CALL spline_fit(spl_intvdlde,'extrap')
            CALL spline_int(spl_intvdlde)
            spl_intvdl%fs(ilmda,:) = spl_intvdlde%fsi(nx,:)
            CALL spline_fit(spl_entvdlde,'extrap')
            CALL spline_int(spl_entvdlde)
            spl_entvdl%fs(ilmda,:) = spl_entvdlde%fsi(nx,:)
            CALL spline_fit(spl_idwdlde,'extrap')
            CALL spline_int(spl_idwdlde)
            spl_idwdl%fs(ilmda,:) = spl_idwdlde%fsi(nx,:)
            CALL spline_fit(spl_edwdlde,'extrap')
            CALL spline_int(spl_edwdlde)
            spl_edwdl%fs(ilmda,:) = spl_edwdlde%fsi(nx,:)
         ENDDO

c-----------------------------------------------------------------------
c     ntv integrand - lambda integration
c-----------------------------------------------------------------------
         CALL spline_fit(spl_intvdl,'extrap')
         CALL spline_int(spl_intvdl)
         intv%fs(istep,:) = spl_intvdl%fsi(nlmda,:)*kindat%f1(1)*
     $        2.0*(nn*kindat%f(3))**2.0/(-echarge*icharge*pi**0.5)
         CALL spline_fit(spl_entvdl,'extrap')
         CALL spline_int(spl_entvdl)
         entv%fs(istep,:) = spl_entvdl%fsi(nlmda,:)*kindat%f1(2)*
     $        2.0*(nn*kindat%f(4))**2.0/(echarge*pi**0.5)
         CALL spline_fit(spl_idwdl,'extrap')
         CALL spline_int(spl_idwdl)
         idw%fs(istep,:) = spl_idwdl%fsi(nlmda,:)*kindat%f1(1)*
     $        2.0*(nn*kindat%f(3))**2.0/(-echarge*icharge*pi**0.5)
         CALL spline_fit(spl_edwdl,'extrap')
         CALL spline_int(spl_edwdl)
         edw%fs(istep,:) = spl_edwdl%fsi(nlmda,:)*kindat%f1(2)*
     $        2.0*(nn*kindat%f(4))**2.0/(echarge*pi**0.5)
      ENDDO
c-----------------------------------------------------------------------
c     ntv integrand - psi integration
c-----------------------------------------------------------------------
      CALL spline_fit(intv,"extrap")
      CALL spline_int(intv)
      CALL spline_fit(entv,"extrap")
      CALL spline_int(entv)
      CALL spline_fit(idw,"extrap")
      CALL spline_int(idw)
      CALL spline_fit(edw,"extrap")
      CALL spline_int(edw)
      
      CALL spline_dealloc(spl_wbdt)
      CALL spline_dealloc(spl_djdt)
      CALL spline_dealloc(spl_intvdlde)
      CALL spline_dealloc(spl_entvdlde)
      CALL spline_dealloc(spl_idwdlde)
      CALL spline_dealloc(spl_edwdlde)
      CALL spline_dealloc(spl_intvdl)
      CALL spline_dealloc(spl_entvdl)
      CALL spline_dealloc(spl_idwdl)
      CALL spline_dealloc(spl_edwdl)


      CALL ascii_open(out_unit,"ipec_ntv_checks_n"//
     $        TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"NTV TORQUE PROFILES:"
      WRITE(out_unit,*)"l = ",0,"npsi = ",npsi/5
      WRITE(OUT_unit,*)"nlmda = ",nlmda/5,"nx = ",nx/10
      WRITE(out_unit,*) ''
      WRITE(out_unit,'(1x,4(1x,a16))') "psi","lmda","x",
     $     "intv"
      DO istep=0,npsi/5
         DO ilmda=0,nlmda/5
            DO i=0,nx/10
               WRITE(out_unit,'(1x,4(1x,es16.8E3))') psis(istep*5),
     $              lmda(ilmda*5),x(i*10),
     $              checks(istep*5,ilmda*5,i)
            ENDDO
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)


      CALL ascii_open(out_unit,"ipec_ntv_n"//
     $        TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"NTV TORQUE PROFILES:"
      DO ll=-nl,nl
         WRITE(out_unit,*) ''
         WRITE(out_unit,*) 'l = ',ll,'itot =',intv%fsi(intv%mx,ll+nl+1),
     $        'etot =',entv%fsi(intv%mx,ll+nl+1)
         WRITE(out_unit,*) ''
         WRITE(out_unit,'(1x,6(1x,a16))') "psi","intv","entv",
     $        "idw","edw","dvdp"
          DO istep=0,npsi
             CALL spline_eval(sq,psis(istep),0)
             WRITE(out_unit,'(1x,6(1x,es16.8E3))') psis(istep),
     $            intv%fs(istep,ll+nl+1),entv%fs(istep,ll+nl+1),
     $            idw%fs(istep,ll+nl+1),edw%fs(istep,ll+nl+1),sq%f(3)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      CALL ascii_open(out_unit,"ipec_ntv_integrated_n"//
     $        TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"INTEGRATED NTV TORQUE PROFILES:"
      DO ll=-nl,nl
         WRITE(out_unit,*) ''
         WRITE(out_unit,*) 'l = ',ll,'itot =',intv%fsi(intv%mx,ll+nl+1),
     $        'etot =',entv%fsi(intv%mx,ll+nl+1)
         WRITE(out_unit,*) ''
         WRITE(out_unit,'(1x,6(1x,a16))') "psi","intv","entv",
     $        "idw","edw","dvdp"
          DO istep=0,npsi
             CALL spline_eval(sq,psis(istep),0)
             WRITE(out_unit,'(1x,6(1x,es16.8E3))') psis(istep),
     $            intv%fsi(istep-1,ll+nl+1),entv%fsi(istep-1,ll+nl+1),
     $            idw%fsi(istep-1,ll+nl+1),edw%fsi(istep-1,ll+nl+1),
     $            sq%f(3)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      CALL ascii_open(out_unit,"ipec_ntv_kinetics.out","UNKNOWN")
      WRITE(out_unit,*)"PROFILES USED IN NTV CALCULATION:"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a12,1x,I4)') "npsi =",npsi
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,13(1x,a16))') "psi", "epsr","ni","ne","ti",
     $     "te","wexb","loglam","nuii","nuei","q","dvdp","dndp"
      DO istep=0,npsi
         CALL spline_eval(sq,psis(istep),0)
         CALL spline_eval(kindat,psis(istep),1)
         CALL bicube_eval(rzphi,psis(istep),0.0_r8,0)
         WRITE(out_unit,'(1x,13(1x,es16.8E3))') psis(istep),
     $        SQRT(rzphi%f(1))/ro,kindat%f(:),sq%f(4),sq%f(3),
     $        kindat%f1(1)
      ENDDO
      CALL ascii_close(out_unit)

      CALL ascii_open(out_unit,"ipec_ntv_bounce_n"//
     $        TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"BOUNCE AVARAGED QUANTITIES:"
      WRITE(out_unit,*)"l = ",0,"npsi = ",npsi,"nlmda = ",nlmda
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,6(1x,a16))') "psi","lmda","wbnce",
     $     "wgrad","wcurv","dj2"
      DO istep=0,npsi
         DO ilmda=0,nlmda
               WRITE(out_unit,'(1x,6(1x,es16.8E3))') psis(istep),
     $              bnc_checks(istep,ilmda,:)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipntv_ntv
c-----------------------------------------------------------------------
c     subprogram 3. ipntv_roots.
c     finds roots of a spline.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ipntv_roots(spl,jroot,xroot)

      TYPE(spline_type) :: spl
      
      REAL(r8), DIMENSION(:,:) :: xroot
      INTEGER, DIMENSION(:,:) :: jroot
      CHARACTER(80) :: message
      INTEGER :: iroot,iqty,it,ix
      INTEGER, PARAMETER :: itmax=20
      REAL(r8), PARAMETER :: eps=1e-10
      REAL(r8) :: x,dx,lx,lf,f,df

c-----------------------------------------------------------------------
c     compute length.
c-----------------------------------------------------------------------
      lx=spl%xs(spl%mx)-spl%xs(0)
c-----------------------------------------------------------------------
c     start loops over iqty.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         iroot = 1
         lf = MAXVAL(spl%fs(:,iqty))-MINVAL(spl%fs(:,iqty))
         IF(MAXVAL(spl%fs(:,iqty))*MINVAL(spl%fs(:,iqty))>=0) CYCLE 
         DO ix=1,spl%mx
c-----------------------------------------------------------------------
c     find all zero passings, intialize at larger gradient.
c-----------------------------------------------------------------------
            IF (spl%fs(ix,iqty)*spl%fs(ix-1,iqty) .LE. 0.0) THEN
               x=spl%xs(ix-1)-spl%fs(ix-1,iqty)*(spl%xs(ix)
     $              -spl%xs(ix-1))/(spl%fs(ix,iqty)-spl%fs(ix-1,iqty))
               f=HUGE(f)
               dx=lx
               it=0
c-----------------------------------------------------------------------
c     locate roots by newton iteration.
c-----------------------------------------------------------------------
               DO
                  CALL spline_eval(spl,x,1)
                  df=spl%f(iqty)-f
                  IF(ABS(dx) < eps*lx 
c     $                 .OR. ABS(df) < eps*lf 
     $                 .OR. it >= itmax)EXIT
                  it=it+1
                  f=spl%f(iqty)
                  dx=-spl%f(iqty)/spl%f1(iqty)
                  x=x+dx
               ENDDO
c-----------------------------------------------------------------------
c     abort on failure.
c-----------------------------------------------------------------------
               IF(it >= itmax)THEN
                  x=spl%xs(ix)
                  IF(spl%fs(ix,iqty) .LT. 0.0) x = spl%xs(ix-1)
                  WRITE(*,*) "WARNING: roots convergence failure!"
                  WRITE(*,*) "-- failure accured at index ",ix
                  WRITE(*,*) "-- indexed x value ",spl%xs(ix)
                  WRITE(*,*) "-- estimated root ",x
               ENDIF
c-----------------------------------------------------------------------
c     store each root, and allocate it to the spline.
c-----------------------------------------------------------------------
               xroot(iroot,iqty) = x
               jroot(iroot,iqty) = ix
               iroot = iroot+1
            ENDIF
            IF (iroot .GT. SIZE(xroot,1)) EXIT
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipntv_roots
c-----------------------------------------------------------------------
c     subprogram 4. ipntv_prepert.
c     Read precalculated plasma perturbations necessary for ntv.
c-----------------------------------------------------------------------
      SUBROUTINE ipntv_readpertfuns(pertfile)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      integer :: istep,itheta,npsi,ntheta
      CHARACTER(128) :: pertfile
      CHARACTER(12) :: label1,label2,label3

      WRITE(*,*) "Reading previously calculated B field perturbations."
      CALL ascii_open(in_unit,pertfile,"OLD")
      READ(in_unit,*) npsi
      READ(in_unit,*) ntheta
      CALL bicube_alloc(pertchebs,npsi-1,ntheta,4)
      DO istep=1,npsi
         DO itheta=0,ntheta
            READ(in_unit,'(6(1x,es16.8))') pertchebs%xs(istep-1),
     $           pertchebs%ys(itheta),pertchebs%fs(istep-1,itheta,:)
         ENDDO
      ENDDO
      CALL ascii_close(in_unit)
      CALL bicube_fit(pertchebs,"extrap","periodic")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipntv_readpertfuns

      END MODULE ipntv_mod
