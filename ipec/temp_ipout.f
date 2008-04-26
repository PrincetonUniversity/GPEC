c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     write various output results of ipec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipout_mod
c      0. ipout_contra
c      0. ipout_ortho
c      0. ipout_bstrength
c      0. ipout-peq
c      0. ipout_contour
c      0. ipout_surfun
c      0. ipout_energy
c      0. ipout_model
c      0. ipout_surfmode
c      0. ipout_delta
c      0. ipout_jumpcoup
c      0. ipout_rawdata
c      0. ipout_3dsurface
c      0. ipeq_orth
c-----------------------------------------------------------------------
c     subprogram 0. ipout_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipout_mod
      USE idcon_mod
      USE ipeq_mod

      IMPLICIT NONE
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipout_contraxi.
c     write contravariant componets of six components in hamada system.
c     write by-products of contravariant components.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c     lobm   : lowest mode number to be observed
c     hobm   : highest mode number to be observed
c     rstep  : radial points to be observed
c     bin    : 0: off binary output for xdraw/(or unix graph)
c            : 1: on  binary output for xdraw/(or unix.graph)
c-----------------------------------------------------------------------
      SUBROUTINE ipout_contra(osol,edgemn,lobm,hobm,rstep,bin)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,lobm,hobm,rstep,bin
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn

      CHARACTER(128) :: message
      
      INTEGER :: i,ipert,istep
      REAL(r8), DIMENSION(:), POINTER :: psistep,rhostep,qstep

      COMPLEX(r8), DIMENSION(:,:), POINTER :: xipsis,xithetas,xizetas,
     $     bpsis,bthetas,bzetas,xiss,dxipsis,dbpsis,ndbpsis,absndbpsis
c-----------------------------------------------------------------------
c     check input parameters.
c-----------------------------------------------------------------------
      message="lobm or hobm should be readjusted"
      IF (((lobm.LT.mlow).OR.(hobm.GT.mhigh)).OR.(lobm.GT.hobm)) THEN
         CALL ipec_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     allocate big matrices.
c-----------------------------------------------------------------------
      ALLOCATE(psistep(0:rstep),rhostep(0:rstep),qstep(0:rstep))
      ALLOCATE(xipsis(0:rstep,mpert),xithetas(0:rstep,mpert),
     $     xizetas(0:rstep,mpert),bpsis(0:rstep,mpert),
     $     bthetas(0:rstep,mpert),bzetas(0:rstep,mpert),
     $     xiss(0:rstep,mpert),
     $     dxipsis(0:rstep,mpert),dbpsis(0:rstep,mpert),
     $     ndbpsis(0:rstep,mpert),absndbpsis(0:rstep,mpert))
c-----------------------------------------------------------------------
c     if rstep=mstep, use original integration points.
c-----------------------------------------------------------------------
      IF (rstep .EQ. mstep) THEN
         psistep=psifac
         rhostep=rho
         qstep=q
      ELSE
         psistep=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
         rhostep=SQRT(psistep)         
      ENDIF
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(osol,edgemn)
      WRITE(*,*)"calculating contravariant components in detail..."
      
      DO istep=0,rstep
         IF (rstep .NE. mstep) THEN
            CALL spline_eval(sq,psistep(istep),0)
            qstep(istep)=sq%f(4)
         ENDIF
         CALL ipeq_contra(psistep(istep))
         xipsis(istep,:)=xipsi
         xithetas(istep,:)=xitheta
         xizetas(istep,:)=xizeta
         bpsis(istep,:)=bpsi
         bthetas(istep,:)=btheta
         bzetas(istep,:)=bzeta
         xiss(istep,:)=xis
         dxipsis(istep,:)=dxipsi
         dbpsis(istep,:)=dbpsi
         ndbpsis(istep,:)=ndbpsi
         absndbpsis(istep,:)=ABS(ndbpsi)
      ENDDO
c-----------------------------------------------------------------------
c     write binary output for xdraw or unix graph.
c-----------------------------------------------------------------------
      IF (bin .EQ. 1) THEN
         CALL bin_open(bin_unit,"ipout_contra.bin",
     $        "UNKNOWN","REWIND","none")
         DO ipert=lobm-mlow+1,hobm-mlow+1
            DO istep=0,rstep
               WRITE(bin_unit)REAL(psistep(istep),4),
     $              REAL(rhostep(istep),4),REAL(qstep(istep),4),
     $              REAL(REAL(xipsis(istep,ipert)),4),
     $              REAL(AIMAG(xipsis(istep,ipert)),4),
     $              REAL(REAL(xithetas(istep,ipert)),4),
     $              REAL(AIMAG(xithetas(istep,ipert)),4),
     $              REAL(REAL(xizetas(istep,ipert)),4),
     $              REAL(AIMAG(xizetas(istep,ipert)),4),
     $              REAL(REAL(bpsis(istep,ipert)),4),
     $              REAL(AIMAG(bpsis(istep,ipert)),4),
     $              REAL(REAL(bthetas(istep,ipert)),4),
     $              REAL(AIMAG(bthetas(istep,ipert)),4),
     $              REAL(REAL(bzetas(istep,ipert)),4),
     $              REAL(AIMAG(bzetas(istep,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
         CALL bin_open(bin_unit,"ipout_contra_add.bin",
     $        "UNKNOWN","REWIND","none")
         DO ipert=lobm-mlow+1,hobm-mlow+1
            DO istep=0,rstep
               WRITE(bin_unit)REAL(psistep(istep),4),
     $              REAL(REAL(xiss(istep,ipert)),4),
     $              REAL(AIMAG(xiss(istep,ipert)),4),
     $              REAL(REAL(dxipsis(istep,ipert)),4),
     $              REAL(AIMAG(dxipsis(istep,ipert)),4),
     $              REAL(REAL(dbpsis(istep,ipert)),4),
     $              REAL(AIMAG(dbpsis(istep,ipert)),4),
     $              REAL(REAL(ndbpsis(istep,ipert)),4),
     $              REAL(AIMAG(ndbpsis(istep,ipert)),4),
     $              REAL(absndbpsis(istep,ipert),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write ascii output for IDL.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     flux functions.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_contra_fun.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTRA_FUN: FLUX FUNCTIONS"
      WRITE(out_unit,'(3(2x,a10,i3))')"mpert=",mpert,
     $     "mlow=",lobm,"mhigh=",hobm
      WRITE(out_unit,'(3(2x,a12))')"psi","rho","q"
      DO istep=0,rstep
         WRITE(out_unit,'(3(2x,es12.3))')
     $        REAL(psistep(istep),4),REAL(rhostep(istep),4),
     $        REAL(qstep(istep),4)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     plasma displacement vectors.
c-----------------------------------------------------------------------    
      CALL ascii_open(out_unit,"ipout_contra_xi.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTRA_XI: DISPLACEMENT VECTORS" 
      WRITE(out_unit,'(7(2x,a12))')"psi","re(tpsi)","im(tpsi)",
     $     "re(ttheta)","im(ttheta)","re(tzeta)","im(tzeta)"
      DO ipert=lobm-mlow+1,hobm-mlow+1
         WRITE(out_unit,'(2x,a3,i3)')"m=",ipert-1+mlow
         DO istep=0,rstep
            WRITE(out_unit,'(7(2x,es12.3))')
     $           REAL(psistep(istep),4),            
     $           REAL(REAL(xipsis(istep,ipert)),4),
     $           REAL(AIMAG(xipsis(istep,ipert)),4),
     $           REAL(REAL(xithetas(istep,ipert)),4),
     $           REAL(AIMAG(xithetas(istep,ipert)),4),
     $           REAL(REAL(xizetas(istep,ipert)),4),
     $           REAL(AIMAG(xizetas(istep,ipert)),4)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)      
c-----------------------------------------------------------------------
c     perturbed magnetic fields.
c----------------------------------------------------------------------- 
      CALL ascii_open(out_unit,"ipout_contra_b.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTRA_B: PERTURBED MAGNETIC FIELDS" 
      WRITE(out_unit,'(7(2x,a12))')"psi","re(bpsi)","im(bpsi)",
     $     "re(btheta)","im(btheta)","re(bzeta)","im(bzeta)" 
      DO ipert=lobm-mlow+1,hobm-mlow+1
         WRITE(out_unit,'(2x,a3,i3)')"m=",ipert-1+mlow
         DO istep=0,rstep
            WRITE(out_unit,'(7(2x,es12.3))')
     $           REAL(psistep(istep),4),    
     $           REAL(REAL(bpsis(istep,ipert)),4),
     $           REAL(AIMAG(bpsis(istep,ipert)),4),
     $           REAL(REAL(bthetas(istep,ipert)),4),
     $           REAL(AIMAG(bthetas(istep,ipert)),4),
     $           REAL(REAL(bzetas(istep,ipert)),4),
     $           REAL(AIMAG(bzetas(istep,ipert)),4)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)   
c-----------------------------------------------------------------------
c     additional quantities.
c----------------------------------------------------------------------- 
      CALL ascii_open(out_unit,"ipout_contra_add.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTRA_ADD: ADDITIONAL INFORMATIONS" 

      WRITE(out_unit,'(10(2x,a12))')"psi","re(xis)","im(xis)",
     $     "re(dxpsi)","im(dxpsi)","re(dbpsi)","im(dbpsi)",
     $     "re(delta)","im(delta)","abs(delta)"
      DO ipert=lobm-mlow+1,hobm-mlow+1
         WRITE(out_unit,'(2x,a3,i3)')"m=",ipert-1+mlow
         DO istep=0,rstep
            WRITE(out_unit,'(10(2x,es12.3))')
     $           REAL(psistep(istep),4),    
     $           REAL(REAL(xiss(istep,ipert)),4),
     $           REAL(AIMAG(xiss(istep,ipert)),4),
     $           REAL(REAL(dxipsis(istep,ipert)),4),
     $           REAL(AIMAG(dxipsis(istep,ipert)),4),
     $           REAL(REAL(dbpsis(istep,ipert)),4),
     $           REAL(AIMAG(dbpsis(istep,ipert)),4),
     $           REAL(REAL(ndbpsis(istep,ipert)),4),
     $           REAL(AIMAG(ndbpsis(istep,ipert)),4),
     $           REAL(absndbpsis(istep,ipert),4)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_contra
c-----------------------------------------------------------------------
c     subprogram 2. ipout_ortho.
c     write orthonormal components.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c     lobm   : lowest mode number to be observed
c     hobm   : highest mode number to be observed
c     rstep  : radial points to be observed
c     bin    : 0: off binary output for xdraw/(or unix graph)
c            : 1: on  binary output for xdraw/(or unix.graph)
c-----------------------------------------------------------------------
      SUBROUTINE ipout_ortho(osol,edgemn,lobm,hobm,rstep,bin)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,lobm,hobm,rstep,bin
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn

      CHARACTER(128) :: message
      
      INTEGER :: i,ipert,istep
      REAL(r8), DIMENSION(:), POINTER :: psistep

      COMPLEX(r8), DIMENSION(:,:), POINTER :: xinorms,bnorms
c-----------------------------------------------------------------------
c     check input parameters.
c-----------------------------------------------------------------------
      message="lobm or hobm should be readjusted"
      IF (((lobm.LT.mlow).OR.(hobm.GT.mhigh)).OR.(lobm.GT.hobm)) THEN
         CALL ipec_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     allocate big matrices.
c-----------------------------------------------------------------------
      ALLOCATE(psistep(0:rstep))
      ALLOCATE(xinorms(0:rstep,mpert),bnorms(0:rstep,mpert))
c-----------------------------------------------------------------------
c     if rstep=mstep, use original integration points.
c-----------------------------------------------------------------------
      IF (rstep .EQ. mstep) THEN
         psistep=psifac
      ELSE
         psistep=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(osol,edgemn)
      WRITE(*,*)"calculating orthonormal components in detail..."
      
      DO istep=0,rstep
         CALL ipeq_contra(psistep(istep))
         CALL ipeq_ortho(psistep(istep))
         xinorms(istep,:)=xinorm
         bnorms(istep,:)=bnorm
         CALL ipeq_tocoord(psistep(istep),xinorms(istep,:),1,0,1,0)
         CALL ipeq_tocoord(psistep(istep),bnorms(istep,:),1,0,1,0)
      ENDDO
c-----------------------------------------------------------------------
c     write binary output for xdraw or unix graph.
c-----------------------------------------------------------------------
      IF (bin .EQ. 1) THEN
         CALL bin_open(bin_unit,"ipout_ortho_norm.bin",
     $        "UNKNOWN","REWIND","none")
         DO ipert=lobm-mlow+1,hobm-mlow+1
            DO istep=0,rstep
               WRITE(bin_unit)REAL(psistep(istep),4),
     $              REAL(REAL(xinorms(istep,ipert)),4),
     $              REAL(AIMAG(xinorms(istep,ipert)),4),
     $              REAL(REAL(bnorms(istep,ipert)),4),
     $              REAL(AIMAG(bnorms(istep,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write ascii output for IDL.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_ortho_norm.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_ORTHO_NORM: NORMAL PERTURBED FIELDS" 
      WRITE(out_unit,'(5(2x,a12))')"psi","re(xinorm)","im(xinorm)",
     $     "re(bnorm)","im(bnorm)"
      DO ipert=lobm-mlow+1,hobm-mlow+1
         WRITE(out_unit,'(2x,a3,i3)')"m=",ipert-1+mlow
         DO istep=0,rstep
            WRITE(out_unit,'(5(2x,es12.3))')
     $           REAL(psistep(istep),4),    
     $           REAL(REAL(xinorms(istep,ipert)),4),
     $           REAL(AIMAG(xinorms(istep,ipert)),4),
     $           REAL(REAL(bnorms(istep,ipert)),4),
     $           REAL(AIMAG(bnorms(istep,ipert)),4)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)   
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_ortho
c-----------------------------------------------------------------------
c     subprogram 3. ipout_bstrength.
c     compute perturbed magnetic field strength.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c     rstep  : radial points to be observed
c     bin    : 0: off binary output for xdraw/(or unix graph)
c            : 1: on  binary output for xdraw/(or unix.graph)
c-----------------------------------------------------------------------
      SUBROUTINE ipout_bstrength(osol,edgemn,rstep,bin)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,rstep,bin
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn
      
      INTEGER :: i,ipert,istep,itheta
      REAL(r8) :: r,z,rfac,angle,jac
      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,wgtfun
      REAL(r8), DIMENSION(1,2) :: w

      REAL(r8), DIMENSION(:), POINTER :: psistep,r0step
      REAL(r8), DIMENSION(:), POINTER :: area,btot

      COMPLEX(r8), DIMENSION(mpert) :: fpsi,ftheta,fzeta,swgtmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bpsif,bthetaf,bzetaf,swgtfun
c-----------------------------------------------------------------------
c     allocate big matrices.
c-----------------------------------------------------------------------
      ALLOCATE(psistep(0:rstep),area(0:rstep),btot(0:rstep),
     $     r0step(0:rstep))
c-----------------------------------------------------------------------
c     if rstep=mstep, use original integration points.
c-----------------------------------------------------------------------
      IF (rstep .EQ. mstep) THEN
         psistep=psifac
      ELSE
         psistep=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(osol,edgemn)
      WRITE(*,*)"calculating bfield strength..."
      
      DO istep=0,rstep
         CALL ipeq_contra(psistep(istep))
         CALL iscdftb(mfac,mpert,bpsif,mthsurf,bpsi)
         CALL iscdftb(mfac,mpert,bthetaf,mthsurf,btheta)
         CALL iscdftb(mfac,mpert,bzetaf,mthsurf,bzeta)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psistep(istep),
     $           itheta/REAL(mthsurf,r8),1)
            rfac=SQRT(rzphi%f(1))
            angle=twopi*(itheta/REAL(mthsurf,r8)+rzphi%f(2))
            r=ro+rfac*COS(angle)
            z=zo+rfac*SIN(angle)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r/jac
            w(1,2)=-rzphi%fy(1)*pi*r/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            wgtfun(itheta)=1.0/(jac*delpsi(itheta))
            IF (itheta .EQ. 0) THEN 
               r0step(istep)=rfac
            ENDIF
         ENDDO
         bpsif=bpsif/wgtfun
         bthetaf=bthetaf/wgtfun
         bzetaf=bzetaf/wgtfun
         swgtfun=1.0/SQRT(wgtfun)
         CALL iscdftf(mfac,mpert,bpsif,mthsurf,fpsi)
         CALL iscdftf(mfac,mpert,bthetaf,mthsurf,ftheta)
         CALL iscdftf(mfac,mpert,bzetaf,mthsurf,fzeta)
         CALL iscdftf(mfac,mpert,swgtfun,mthsurf,swgtmn)
         CALL ipeq_cova(psistep(istep))
         btot(istep)=0
         area(istep)=0
         DO i=1,mpert
            btot(istep)=btot(istep)+
     $           REAL(CONJG(bcpsi(i))*fpsi(i)+
     $           CONJG(bctheta(i))*ftheta(i)+CONJG(bczeta(i))*fzeta(i))
            area(istep)=area(istep)+
     $           REAL(CONJG(swgtmn(i))*swgtmn(i))
         ENDDO
         btot(istep)=SQRT(ABS(btot(istep))/area(istep))
      ENDDO
c-----------------------------------------------------------------------
c     write binary output for xdraw or unix graph.
c-----------------------------------------------------------------------
      IF (bin .EQ. 1) THEN
         CALL bin_open(bin_unit,"ipout_bstrength.bin",
     $        "UNKNOWN","REWIND","none")
         DO istep=0,rstep
            WRITE(bin_unit)REAL(psistep(istep),4),REAL(btot(istep),4)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write ascii output for IDL.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_strength.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_BSTRENGTH: STRENGTH OF PERTURBED FIELDS" 
      WRITE(out_unit,'(3(2x,a12))')"psi","bstrength","rfac"
      DO istep=0,rstep
         WRITE(out_unit,'(3(2x,es12.3))')
     $        REAL(psistep(istep),4),
     $        REAL(btot(istep),4),REAL(r0step(istep),4)
      ENDDO
      CALL ascii_close(out_unit)   
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_bstrength
c-----------------------------------------------------------------------
c     subprogram 4. ipout_peq.
c     write data for perturbed equilibrium plot.
c     final solutions use polar-phi in 3d picture.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipout_peq(osol,edgemn,weight)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,weight
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn

      INTEGER :: i,j,istep,itheta,ipert,ising,rstep,stepnum,estepnum
      REAL(r8) :: rfac,angle,jac,delpsi,lpsi,rpsi,
     $     step,estep,range,erange,a1,a2,c1,delphi

      COMPLEX(r8) :: xipsif,xithetaf,xizetaf,
     $     t_xipsif,t_xithetaf,xitorof,
     $     bpsif,bthetaf,xinormf,bnormf

      REAL(r8), DIMENSION(:), POINTER :: psi
      REAL(r8), DIMENSION(:), POINTER :: psis

      REAL(r8), DIMENSION(0:mthsurf) :: theta,dphi
      REAL(r8), DIMENSION(1,2) :: nc,nr
      REAL(r8), DIMENSION(2,2) :: vc,vr

      REAL(r8), DIMENSION(:,:), POINTER :: r,z,cxvr,cxvz,sxvr,sxvz,
     $     cxnr,cxnz,sxnr,sxnz
      COMPLEX(r8), DIMENSION(:,:), POINTER :: xivr,xivz,
     $     t_xivr,t_xivz,xinr,xinz,bvr,bvz,bnr,bnz
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(osol,edgemn)
c-----------------------------------------------------------------------
c     assign proper radial points for contour plot.
c-----------------------------------------------------------------------
      i=1
      ALLOCATE(psis(100))
      psis=0
      psis(i)=0.002
      lpsi=psis(i)
      DO ising=1,msing
         rpsi=singtype(ising)%psifac-0.01
         range=rpsi-lpsi
         IF (range .GE. 0.1) THEN
            stepnum=INT(range/0.1)+1
            step=range/stepnum
            DO j=1,stepnum
               i=i+1
               psis(i)=psis(i-1)+step
            ENDDO
            lpsi=psis(i)
         ELSE IF ((range .GE. 0.04) .AND. (range .LT. 0.1)) THEN
            i=i+1
            psis(i)=rpsi
            lpsi=psis(i)
         ELSE
            EXIT
         ENDIF
      ENDDO
      i=i+1
      psis(i)=(psilim-0.001+lpsi)/2.0
      i=i+1
      psis(i)=psilim-0.001
      rstep=i
      ALLOCATE(psi(rstep))
      DO istep=1,rstep
         psi(istep)=psis(istep)
      ENDDO
      DEALLOCATE(psis)
      ALLOCATE(r(rstep,0:mthsurf),z(rstep,0:mthsurf),
     $     cxvr(rstep,0:mthsurf),cxvz(rstep,0:mthsurf),
     $     sxvr(rstep,0:mthsurf),sxvz(rstep,0:mthsurf),
     $     cxnr(rstep,0:mthsurf),cxnz(rstep,0:mthsurf),
     $     sxnr(rstep,0:mthsurf),sxnz(rstep,0:mthsurf),
     $     xivr(rstep,0:mthsurf),xivz(rstep,0:mthsurf),
     $     t_xivr(rstep,0:mthsurf),t_xivz(rstep,0:mthsurf),
     $     xinr(rstep,0:mthsurf),xinz(rstep,0:mthsurf),
     $     bvr(rstep,0:mthsurf),bvz(rstep,0:mthsurf),
     $     bnr(rstep,0:mthsurf),bnz(rstep,0:mthsurf))
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing perturbed equilibrium 2d plot..."
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      DO istep=1,rstep
         CALL ipeq_contra(psi(istep))  
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            angle=twopi*(theta(itheta)+rzphi%f(2))
            r(istep,itheta)=ro+rfac*COS(angle)
            z(istep,itheta)=zo+rfac*SIN(angle)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            vc(1,1)=1/(2*rfac)*rzphi%fx(1)
            vc(1,2)=twopi*rfac*rzphi%fx(2)
            vc(2,1)=1/(2*rfac)*rzphi%fy(1)
            vc(2,2)=twopi*rfac*rzphi%fy(2)
            vr(1,1)=cos(angle)*vc(1,1)-sin(angle)*vc(1,2)
            vr(1,2)=sin(angle)*vc(1,1)+cos(angle)*vc(1,2)
            vr(2,1)=cos(angle)*vc(2,1)-sin(angle)*vc(2,2)
            vr(2,2)=sin(angle)*vc(2,1)+cos(angle)*vc(2,2)
            
            nc(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(istep,itheta)/jac
            nc(1,2)=-rzphi%fy(1)*pi*r(istep,itheta)/(rfac*jac)
            delpsi=SQRT(nc(1,1)**2+nc(1,2)**2)
            nc=nc/delpsi
            nr(1,1)=cos(angle)*nc(1,1)-sin(angle)*nc(1,2)
            nr(1,2)=sin(angle)*nc(1,1)+cos(angle)*nc(1,2)

            xipsif=0
            xithetaf=0
            xizetaf=0
            bpsif=0
            bthetaf=0
            DO ipert=1,mpert
               xipsif=xipsif+xipsi(ipert)*EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
               xithetaf=xithetaf+xitheta(ipert)*
     $              EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
               xizetaf=xizetaf+xizeta(ipert)*
     $              EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
               bpsif=bpsif+bpsi(ipert)*EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
               bthetaf=bthetaf+btheta(ipert)*
     $              EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
            ENDDO
            xipsif=xipsif*weight
            xithetaf=xithetaf*weight
            xizetaf=xizetaf*weight
            xinormf=xipsif/delpsi
            xitorof=xizetaf*twopi*r(istep,itheta)
c-----------------------------------------------------------------------
c     delphi calculation in linear regime.
c-----------------------------------------------------------------------
            a1=REAL(xitorof)
            a2=AIMAG(xitorof)
            c1=(a2-r(istep,itheta))/a1
            IF (c1 .GT. 0) THEN
               delphi=-c1+SQRT(c1**2+2)
            ELSE
               delphi=-c1-SQRT(c1**2+2)
            ENDIF
c-----------------------------------------------------------------------
c     corrections by toroidal distortion.
c-----------------------------------------------------------------------
            t_xipsif=xipsif*EXP(ifac*nn*delphi)
            t_xithetaf=xithetaf*EXP(ifac*nn*delphi)

            bpsif=bpsif*weight
            bthetaf=bthetaf*weight
            bnormf=bpsif/delpsi

            xivr(istep,itheta)=xipsif*vr(1,1)+xithetaf*vr(2,1)
            xivz(istep,itheta)=xipsif*vr(1,2)+xithetaf*vr(2,2)

            t_xivr(istep,itheta)=t_xipsif*vr(1,1)+t_xithetaf*vr(2,1)
            t_xivz(istep,itheta)=t_xipsif*vr(1,2)+t_xithetaf*vr(2,2)
            t_xivr(istep,itheta)=t_xivr(istep,itheta)+
     $           r(istep,itheta)*(1.0/cos(delphi)-1)

            cxvr(istep,itheta)=r(istep,itheta)+REAL(xivr(istep,itheta))
            cxvz(istep,itheta)=z(istep,itheta)+REAL(xivz(istep,itheta))
            sxvr(istep,itheta)=r(istep,itheta)+AIMAG(xivr(istep,itheta))
            sxvz(istep,itheta)=z(istep,itheta)+AIMAG(xivz(istep,itheta))

            xinr(istep,itheta)=xinormf*nr(1,1)
            xinz(istep,itheta)=xinormf*nr(1,2)

            cxnr(istep,itheta)=r(istep,itheta)+REAL(xinr(istep,itheta))
            cxnz(istep,itheta)=z(istep,itheta)+REAL(xinz(istep,itheta))
            sxnr(istep,itheta)=r(istep,itheta)+AIMAG(xinr(istep,itheta))
            sxnz(istep,itheta)=z(istep,itheta)+AIMAG(xinz(istep,itheta))     

            bvr(istep,itheta)=bpsif*vr(1,1)+bthetaf*vr(2,1)
            bvz(istep,itheta)=bpsif*vr(1,2)+bthetaf*vr(2,2)

            bnr(istep,itheta)=bnormf*nr(1,1)
            bnz(istep,itheta)=bnormf*nr(1,2)          
          
         ENDDO
      ENDDO      
      CALL ascii_open(out_unit,"ipout_peq_xinorm_cont.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PEQ_XINORM_CONT: "//
     $     "contour plots by concerning normal displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","real(xinr)","real(xinz)",
     $     "imag(xinr)","imag(xinz)","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           REAL(xinr(istep,itheta)),REAL(xinz(istep,itheta)),
     $           AIMAG(xinr(istep,itheta)),AIMAG(xinz(istep,itheta)),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_peq_xivec_cont.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PEQ_XIVEC_CONT: "//
     $     "contour plots by concerning actual displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","real(xivr)","real(xivz)",
     $     "imag(xivr)","imag(xivz)","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           REAL(xivr(istep,itheta)),REAL(xivz(istep,itheta)),
     $           AIMAG(xivr(istep,itheta)),AIMAG(xivz(istep,itheta)),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_peq_txivec_cont.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PEQ_TXIVEC_CONT: "//
     $     "contour plots by concerning actual displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","real(txivr)","real(txivz)",
     $     "imag(txivr)","imag(txivz)","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           REAL(t_xivr(istep,itheta)),REAL(t_xivz(istep,itheta)),
     $           AIMAG(t_xivr(istep,itheta)),
     $           AIMAG(t_xivz(istep,itheta)),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_peq_bnorm_cont.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PEQ_BNORM_CONT: "//
     $     "contour plots by concerning normal displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","real(bnr)","real(bnz)",
     $     "imag(bnr)","imag(bnz)","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           REAL(bnr(istep,itheta)),REAL(bnz(istep,itheta)),
     $           AIMAG(bnr(istep,itheta)),AIMAG(bnz(istep,itheta)),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_peq_bvec_cont.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PEQ_BVEC_CONT: "//
     $     "contour plots by concerning actual displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","real(bvr)","real(bvz)",
     $     "imag(bvr)","imag(bvz)","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           REAL(bvr(istep,itheta)),REAL(bvz(istep,itheta)),
     $           AIMAG(bvr(istep,itheta)),AIMAG(bvz(istep,itheta)),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)

      DEALLOCATE(r,z,cxvr,cxvz,sxvr,sxvz,cxnr,cxnz,sxnr,sxnz,
     $     xivr,xivz,t_xivr,t_xivz,xinr,xinz,bvr,bvz,bnr,bnz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_peq
c-----------------------------------------------------------------------
c     subprogram 5. ipout_contour.
c     write data for 2d plots of real perturbed quantities.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipout_contour(osol,edgemn,rstep,weight)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,weight,rstep
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn

      INTEGER :: i,istep,itheta,ipert
      REAL(r8) :: rfac,angle,jac,delpsi
      COMPLEX(r8) :: xinormf,bnormf

      REAL(r8), DIMENSION(:), POINTER :: psi
      REAL(r8), DIMENSION(0:mthsurf) :: theta,dphi
      
      REAL(r8), DIMENSION(1,2) :: nc
      REAL(r8), DIMENSION(:,:), POINTER :: r,z,xiabsf,xiangf,
     $     babsf,bangf
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(osol,edgemn)
c-----------------------------------------------------------------------
c     allocate big matrices.
c-----------------------------------------------------------------------
      ALLOCATE(psi(0:rstep))
      ALLOCATE(r(0:rstep,0:mthsurf),z(0:rstep,0:mthsurf),
     $     xiabsf(0:rstep,0:mthsurf),xiangf(0:rstep,0:mthsurf),
     $     babsf(0:rstep,0:mthsurf),bangf(0:rstep,0:mthsurf))
c-----------------------------------------------------------------------
c     if rstep=mstep, use original integration points.
c-----------------------------------------------------------------------
      IF (rstep .EQ. mstep) THEN
         psi=psifac
      ELSE
         psi=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing perturbed equilibrium 2d contour.."
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      DO istep=0,rstep
         CALL ipeq_contra(psi(istep))
         CALL ipeq_ortho(psi(istep))
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            angle=twopi*(theta(itheta)+rzphi%f(2))
            r(istep,itheta)=ro+rfac*COS(angle)
            z(istep,itheta)=zo+rfac*SIN(angle)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            
            nc(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(istep,itheta)/jac
            nc(1,2)=-rzphi%fy(1)*pi*r(istep,itheta)/(rfac*jac)
            delpsi=SQRT(nc(1,1)**2+nc(1,2)**2)

            xinormf=0
            bnormf=0

            DO ipert=1,mpert
               xinormf=xinormf+xipsi(ipert)*
     $              EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
               bnormf=bnormf+bpsi(ipert)*
     $              EXP(twopi*ifac*(ipert+mlow-1)*
     $              theta(itheta)+nn*ifac*dphi(itheta))
            ENDDO
            xinormf=xinormf/delpsi*weight
            bnormf=bnormf/delpsi*weight

            xiabsf(istep,itheta)=
     $           SQRT(REAL(xinormf)**2+AIMAG(xinormf)**2)
            xiangf(istep,itheta)=ATAN(AIMAG(xinormf)/REAL(xinormf))
            babsf(istep,itheta)=
     $           SQRT(REAL(bnormf)**2+AIMAG(bnormf)**2)
            bangf(istep,itheta)=ATAN(AIMAG(bnormf)/REAL(bnormf))
         ENDDO
      ENDDO      
      
      CALL ascii_open(out_unit,"ipout_contour.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTOUR: "//
     $     "contour plots for only normal displacements"
      WRITE(out_unit,'(7(2x,a12))')"r","z","xiabsf","xiangf",
     $     "babsf","bangf","psi"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')
     $           r(istep,itheta),z(istep,itheta),
     $           xiabsf(istep,itheta),xiangf(istep,itheta),
     $           babsf(istep,itheta),bangf(istep,itheta),
     $           psi(istep)
         ENDDO
      ENDDO 
      CALL ascii_close(out_unit)

      DEALLOCATE(r,z,xiabsf,xiangf,babsf,bangf,psi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_contour
c-----------------------------------------------------------------------
c     subprogram 6. ipout_surfun.
c     write the surface functions for eigenmodes
c     __________________________________________________________________
c     lowmode  : lowest lable number to be observed
c     highmode : highest lable number to be observed
c     polo     : 0: polar angle
c                1: hamada
c                2: pest
c                3: equal arc length
c                4: boozer
c     toro     : 0: polar angle
c                1: hamada
c-----------------------------------------------------------------------
      SUBROUTINE ipout_surfun(lowmode,highmode,polo,toro)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: lowmode,highmode,polo,toro
      INTEGER :: i,j

      CHARACTER(1) :: spolo,storo

      REAL(r8), DIMENSION(0:mthsurf) :: theta
      COMPLEX(r8), DIMENSION(mpert) :: bnomn,chimn,chemn,chpmn,
     $     kapmn,kaxmn
      COMPLEX(r8), DIMENSION(lowmode:highmode,0:mthsurf) :: chi_funmats,
     $     che_funmats,chp_funmats,kap_funmats,kax_funmats,
     $     bno_funmats,flx_funmats

      theta=twopi*(/(i,i=0,mthsurf)/)/REAL(mthsurf,r8)

      DO i=lowmode,highmode
         bnomn=bnomats(:,i)
         chimn=chimats(:,i)
         chemn=chemats(:,i)
         chpmn=chpmats(:,i)
         kapmn=kapmats(:,i)
         kaxmn=kaxmats(:,i)
         CALL ipeq_tocoord(psilim,bnomn,polo,toro,1,0)
         CALL ipeq_tocoord(psilim,chimn,polo,toro,1,0)
         CALL ipeq_tocoord(psilim,chemn,polo,toro,1,0)
         CALL ipeq_tocoord(psilim,chpmn,polo,toro,1,0)
         CALL ipeq_tocoord(psilim,kapmn,polo,toro,1,0)
         CALL ipeq_tocoord(psilim,kaxmn,polo,toro,1,0)
         CALL iscdftb(mfac,mpert,bno_funmats(:,i),mthsurf,bnomn)
         CALL iscdftb(mfac,mpert,chi_funmats(:,i),mthsurf,chimn)
         CALL iscdftb(mfac,mpert,che_funmats(:,i),mthsurf,chemn)
         CALL iscdftb(mfac,mpert,chp_funmats(:,i),mthsurf,chpmn)
         CALL iscdftb(mfac,mpert,kap_funmats(:,i),mthsurf,kapmn)
         CALL iscdftb(mfac,mpert,kax_funmats(:,i),mthsurf,kaxmn)
      ENDDO
         
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro

      CALL ascii_open(out_unit,"ipout_surfun_p"//spolo//"_t"//
     $     storo//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_SURFUN: surface functions"
      WRITE(out_unit,*)"poloidal coordinate:",polo
      WRITE(out_unit,*)"toroidal coordinate:",toro
      WRITE(out_unit,'(7(2x,a12))')"theta","flux","internal","external",
     $     "plasma","plasmasurf","vacuumsurf" 
      DO i=lowmode,highmode
         WRITE(out_unit,'(2x,i2,a10)')i,"th mode"
         DO j=0,mthsurf
            WRITE(out_unit,'(7(2x,es12.3))')theta(j),
     $           REAL(bno_funmats(j,i)),REAL(chi_funmats(j,i)),
     $           REAL(che_funmats(j,i)),REAL(chp_funmats(j,i)),
     $           REAL(kap_funmats(j,i)),REAL(kax_funmats(j,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_surfun
c-----------------------------------------------------------------------
c     subprogram 7 ipout_energy.
c     write energy information.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_energy
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i

      CALL ascii_open(out_unit,"ipout_energy.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_ENERGY: comparison between energy from D
     $CON eigenmodes and IPEC surface eigenmodes"
      WRITE(out_unit,'(3(2x,a10,i3))')"mpert=",mpert,
     $     "mlow=",mlow,"mhigh=",mhigh
      WRITE(out_unit,*)"EIGENENERGIES"      
      WRITE(out_unit,'(6(2x,a12))')"surfep","ep","surfev","ev",
     $     "surfet","et"
c-----------------------------------------------------------------------
c     eigenenergies from (4+1)th lines.
c-----------------------------------------------------------------------
      DO i=1,mpert
         WRITE(out_unit,'(6(2x,es12.3))')surfep(i),ep(i),
     $        surfew(i),ew(i),surfet(i),et(i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_energy
c-----------------------------------------------------------------------
c     subprogram 8. ipout_resp.
c     write inductances and permeability matrix information.
c     use hamada-weighted coordinate system.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_resp
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j

      CALL ascii_open(out_unit,"ipout_resp.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_RESP: inductances, permeability matri
     $ces and their eigenvalues"
      WRITE(out_unit,'(3(2x,a10,i3))')"mpert=",mpert,
     $     "mlow=",mlow,"mhigh=",mhigh
      WRITE(out_unit,*)"EIGENVALUES"
      WRITE(out_unit,'(7(2x,a12))')"surf_ind","plas_ind","engy_ind",
     $     "permeab","eermeab","reluct","eeluct"
c-----------------------------------------------------------------------
c     eigenvalues from (4+1)th lines.
c-----------------------------------------------------------------------
      DO i=1,mpert
         WRITE(out_unit,'(7(2x,es12.3))')surf_indev(i),
     $        plas_indev(i),engy_indev(i),permeabev(i),eermeabev(i),
     $        reluctev(i),eeluctev(i)
      ENDDO
c      WRITE(out_unit,*)
c      WRITE(out_unit,*)"MATRIX COMPONENTS"
c      WRITE(out_unit,'(4(2x,a26))')"energy","plasma",
c     $     "surface","permeability"
c      WRITE(out_unit,'(8(2x,a12))')"re","im","re","im","re","im",
c     $     "re","im"
c-----------------------------------------------------------------------
c     matrix components from (8+mpert+1)th lines.
c-----------------------------------------------------------------------
c      DO i=1,mpert
c         WRITE(out_unit,'(2x,a12,i2)')"row number=",i
c         DO j=1,mpert
c            WRITE(out_unit,'(8(2x,es12.3))')
c     $           REAL(engy_indmats(j,i)),AIMAG(engy_indmats(j,i)),
c     $           REAL(plas_indmats(j,i)),AIMAG(plas_indmats(j,i)),
c     $           REAL(surf_indmats(j,i)),AIMAG(surf_indmats(j,i)),
c     $           REAL(permeabmats(j,i)),AIMAG(permeabmats(j,i))
c         ENDDO
c      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_resp
c-----------------------------------------------------------------------
c     subprogram 9. ipout_surfmode.
c     reponse to fourier modes for the control surface.
c     __________________________________________________________________
c     lowmode  : lowest number of m fourier mode applied
c     highmode : highest number of m fourier mode applied
c     polo     : 0: polar angle
c                1: hamada
c                2: pest
c                3: equal arc length
c                4: boozer
c     toro     : 0: polar angle
c                1: hamada
c-----------------------------------------------------------------------
      SUBROUTINE ipout_surfmode(lowmode,highmode,polo,toro)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: lowmode,highmode,polo,toro

      INTEGER :: i,j,mnum
      CHARACTER(1) :: spolo,storo

      REAL(r8), DIMENSION(0:mthsurf) :: theta
      COMPLEX(r8), DIMENSION(:,:), POINTER :: binmn,boutmn,finmn,
     $     binfun,boutfun

      theta=twopi*(/(i,i=0,mthsurf)/)/REAL(mthsurf,r8)
      mnum=highmode-lowmode+1
      ALLOCATE(binmn(mnum,mpert),boutmn(mnum,mpert),finmn(mnum,mpert),
     $     binfun(mnum,0:mthsurf),boutfun(mnum,0:mthsurf))
     
      DO i=lowmode,highmode
         j=i-lowmode+1
         binmn(j,:)=0
         binmn(j,i-mlow+1)=1
         finmn(j,:)=binmn(j,:)
         CALL iscdftb(mfac,mpert,binfun(j,:),mthsurf,finmn(j,:))         
         CALL ipeq_tobasis(psilim,finmn(j,:),polo,toro,1,1)
         boutmn(j,:)=MATMUL(permeabmats,finmn(j,:))
         CALL ipeq_tocoord(psilim,boutmn(j,:),polo,toro,1,1)
         CALL iscdftb(mfac,mpert,boutfun(j,:),mthsurf,boutmn(j,:))
      ENDDO
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      CALL ascii_open(out_unit,"ipout_surfmode_p"//spolo//"_t"
     $     //storo//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_SURFMODE: plasma response for the fourier 
     $modes on the control surface"
      WRITE(out_unit,*)"poloidal coordinate:",polo
      WRITE(out_unit,*)"toroidal coordinate:",toro
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(5(2x,a12))')"m",
     $     "rebin","imbin","rebout","imbout"
      DO j=1,mnum
         WRITE(out_unit,'(2x,a6,i3)')"mode=",j-1+lowmode
         DO i=1,mpert
            WRITE(out_unit,'(2x,I12,4(2x,es12.3))')mfac(i),
     $           REAL(binmn(j,i)),AIMAG(binmn(j,i)),
     $           REAL(boutmn(j,i)),AIMAG(boutmn(j,i))
         ENDDO
      ENDDO
      WRITE(out_unit,*)"FUNCTIONS"
      WRITE(out_unit,'(5(2x,a12))')"theta",
     $     "rebin","imbin","rebout","imbout"
      DO j=1,mnum
         WRITE(out_unit,'(2x,a6,i3)')"mode=",j-1+lowmode
         DO i=0,mthsurf
            WRITE(out_unit,'(5(2x,es12.3))')theta(i),
     $           REAL(binfun(j,i)),AIMAG(binfun(j,i)),
     $           REAL(boutfun(j,i)),AIMAG(boutfun(j,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_surfmode
c-----------------------------------------------------------------------
c     subprogram 10. ipout_delta.
c     measure asymtotic values of delta quantities
c     check every singular surfaces
c-----------------------------------------------------------------------
      SUBROUTINE ipout_delta(osol,edgemn,rsing,resol)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: osol,resol,rsing
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: edgemn

      INTEGER :: ising,i,itheta,resnum
      REAL(r8) :: respsi,resrpsi,reslpsi,spsi,astep,lpsi,rpsi,
     $     rfac,angle,r,z,jac,sqrpsi,correc
      COMPLEX(r8) :: lndbpsi,rndbpsi,lcorrec,rcorrec
     
      REAL(r8), DIMENSION(0:mthsurf) :: theta
      COMPLEX(r8), DIMENSION(mpert) :: cormn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bpsifun,corfun

      REAL(r8), DIMENSION(3,2) :: w
      REAL(r8), DIMENSION(rsing,resol) :: dist,absdelta,abstotdelta
      COMPLEX(r8), DIMENSION(rsing,resol) :: delta,deldelta,totdelta

      CALL idcon_build(osol,edgemn)

      WRITE(*,*)"checking asymtotic values of deltas..."
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      DO ising=1,rsing
         resnum=NINT(singtype(ising)%q)*nn-mlow+1
         respsi=singtype(ising)%psifac
        
         DO i=1,resol
            astep=0.1*10.0**(-1-(i-1)*5.0/(resol-1))
            lpsi=respsi-astep
            rpsi=respsi+astep
            dist(ising,i)=astep

            CALL ipeq_contra(lpsi)
            lndbpsi=ndbpsi(resnum)
c-----------------------------------------------------------------------
c     calculate the correction term in boozer's theory for lpsi.
c-----------------------------------------------------------------------
            CALL iscdftb(mfac,mpert,bpsifun,mthsurf,bpsi)    
            CALL spline_eval(sq,lpsi,0)
            DO itheta=0,mthsurf
               CALL bicube_eval(rzphi,lpsi,theta(itheta),1)
               rfac=SQRT(rzphi%f(1))
               angle=twopi*(theta(itheta)+rzphi%f(2))
               r=ro+rfac*COS(angle)
               z=zo+rfac*SIN(angle)
               jac=rzphi%f(4)
               w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r/jac
               w(1,2)=-rzphi%fy(1)*pi*r/(rfac*jac)
               w(2,1)=-rzphi%fx(2)*twopi**2*r*rfac/jac
               w(2,2)=rzphi%fx(1)*pi*r/(rfac*jac)
               w(3,1)=rzphi%fx(3)*2*rfac
               w(3,2)=rzphi%fy(3)/(twopi*rfac)
               sqrpsi=w(1,1)**2+w(1,2)**2
               correc=twopi*jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
               corfun(itheta)=bpsifun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,corfun,mthsurf,cormn)
            lcorrec=ifac*mfac(resnum)*cormn(resnum)/(chi1*sq%f(4))

            CALL ipeq_contra(rpsi)
            rndbpsi=ndbpsi(resnum)
c-----------------------------------------------------------------------
c     calculate the correction term in boozer's theory for rpsi.
c-----------------------------------------------------------------------
            CALL iscdftb(mfac,mpert,bpsifun,mthsurf,bpsi)    
            CALL spline_eval(sq,rpsi,0)
            DO itheta=0,mthsurf
               CALL bicube_eval(rzphi,rpsi,theta(itheta),1)
               rfac=SQRT(rzphi%f(1))
               angle=twopi*(theta(itheta)+rzphi%f(2))
               r=ro+rfac*COS(angle)
               z=zo+rfac*SIN(angle)
               jac=rzphi%f(4)
               w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r/jac
               w(1,2)=-rzphi%fy(1)*pi*r/(rfac*jac)
               w(2,1)=-rzphi%fx(2)*twopi**2*r*rfac/jac
               w(2,2)=rzphi%fx(1)*pi*r/(rfac*jac)
               w(3,1)=rzphi%fx(3)*2*rfac
               w(3,2)=rzphi%fy(3)/(twopi*rfac)
               sqrpsi=w(1,1)**2+w(1,2)**2
               correc=twopi*jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
               corfun(itheta)=bpsifun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,corfun,mthsurf,cormn)
            rcorrec=ifac*mfac(resnum)*cormn(resnum)/(chi1*sq%f(4))

            delta(ising,i)=rndbpsi-lndbpsi
            deldelta(ising,i)=rcorrec-lcorrec
            totdelta(ising,i)=delta(ising,i)+deldelta(ising,i)

            absdelta(ising,i)=ABS(delta(ising,i))
            abstotdelta(ising,i)=ABS(totdelta(ising,i))
         ENDDO
         WRITE(*,*)"finish a singular surface for the q=",
     $        singtype(ising)%q
      ENDDO
c-----------------------------------------------------------------------
c     write the results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_delta.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_DELTA: asymtotic deltas"
      WRITE(out_unit,'(5(2x,a12))')"distance","re(delta)",
     $     "im(delta)","abs(delta)","abs(tdelta)"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,a2,f6.3)')"q=",singtype(ising)%q
         DO i=1,resol
            WRITE(out_unit,'(5(2x,es12.3))')dist(ising,i),
     $           REAL(delta(ising,i)),AIMAG(delta(ising,i)),
     $           absdelta(ising,i),abstotdelta(ising,i)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_delta
c-----------------------------------------------------------------------
c     subprogram 11. ipout_jumpcoup.
c     obtain mode coupling between magnetic jumps on the singular
c     surfaces and fourier modes on the control surface
c     __________________________________________________________________
c     eps    : distance from the rational surface
c     rsing  : number of rational surfaces from the lowest
c     polo   : 0: polar angle
c              1: hamada
c              2: pest
c              3: equal arc length
c              4: boozer
c     toro   : 0: polar angle
c              1: hamada
c     resp   : 0: without plasma response
c              1: with plasma response
c-----------------------------------------------------------------------
      SUBROUTINE ipout_jumpcoup(eps,rsing,polo,toro,resp)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rsing,polo,toro,resp
      REAL(r8), INTENT(IN) :: eps

      INTEGER :: i,ising,resnum,lwork
      REAL(r8) :: jac,respsi,lpsi,rpsi
      COMPLEX(r8) :: lndbpsi,rndbpsi
      CHARACTER(1) :: spolo,storo,sresp

      REAL(r8), DIMENSION(rsing) :: s
      REAL(r8), DIMENSION(5*mpert) :: rwork
      REAL(r8), DIMENSION(0:mthsurf) :: theta

      COMPLEX(r8), DIMENSION(3*rsing+mpert) :: work
      COMPLEX(r8), DIMENSION(rsing,rsing) :: u
      COMPLEX(r8), DIMENSION(rsing,mpert) :: delta,a
      COMPLEX(r8), DIMENSION(mpert,mpert) :: vt,bvt,pvt
      COMPLEX(r8), DIMENSION(rsing,0:mthsurf) :: vtfun,pvtfun

      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn
      
      REAL(r8), DIMENSION(rsing,mpert) :: absdelta

      CALL spline_eval(sq,psilim,0)
      CALL bicube_eval(rzphi,psilim,thetam,0)
      singfac=mfac-nn*sq%f(4)
      jac=rzphi%f(4)
c-----------------------------------------------------------------------
c     evaluate the coupling matrix.
c-----------------------------------------------------------------------
      DO i=1,mpert
         binmn=0
         binmn(i)=1.0*gauss
         CALL ipeq_tobasis(psilim,binmn,polo,toro,1,1)
         IF (resp .EQ. 1) THEN
            boutmn=MATMUL(permeabmats,binmn)
         ELSE
            boutmn=binmn
         ENDIF
         CALL ipeq_tocoord(psilim,boutmn,1,1,1,1)
         CALL ipeq_bntox(psilim,boutmn)
         edge_flag=.TRUE.
         CALL idcon_build(0,boutmn)

         DO ising=1,rsing
            resnum=NINT(singtype(ising)%q)*nn-mlow+1
            respsi=singtype(ising)%psifac

            lpsi=respsi-eps
            CALL ipeq_contra(lpsi)
            lndbpsi=ndbpsi(resnum)
            rpsi=respsi+eps
            CALL ipeq_contra(rpsi)
            rndbpsi=ndbpsi(resnum)
            delta(ising,i)=rndbpsi-lndbpsi
            absdelta(ising,i)=ABS(delta(ising,i))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     svd analysis.
c-----------------------------------------------------------------------
      a=delta
      work=0
      s=0
      u=0
      vt=0
      bvt=0
      pvt=0
      lwork=5*rsing+mpert      
      CALL zgesvd('S','S',rsing,mpert,a,rsing,s,u,rsing,vt,mpert,
     $     work,lwork,rwork,info)
      vt=CONJG(vt)
      bvt=vt
      DO ising=1,rsing
         CALL iscdftb(mfac,mpert,vtfun(ising,:),mthsurf,vt(ising,:))
         CALL ipeq_tobasis(psilim,bvt(ising,:),polo,toro,1,1)
         pvt(ising,:)=MATMUL(permeabmats,bvt(ising,:))
         CALL ipeq_tocoord(psilim,pvt(ising,:),polo,toro,1,1)
         CALL iscdftb(mfac,mpert,pvtfun(ising,:),mthsurf,pvt(ising,:))
      ENDDO      
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=sresp, FMT='(I1)')resp

      CALL ascii_open(out_unit,"ipout_jumpcoup_p"//spolo//"_t"
     $     //storo//"_r"//sresp//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_JUMP COUPLINGS: "
      WRITE(out_unit,*)"epsilon: ",eps
      WRITE(out_unit,*)"poloidal coordinate:",polo
      WRITE(out_unit,*)"toroidal coordinate:",toro
      WRITE(out_unit,*)"plasma response:",resp

      DO ising=1,rsing
         WRITE(out_unit,'(2x,a2,f6.3)')"q=",singtype(ising)%q 
         WRITE(out_unit,'(2x,a3,3(2x,a12))')"m",
     $        "real","imaginary","absolute"
         DO i=1,mpert
            WRITE(out_unit,'(2x,I3,3(2x,es12.3))')mfac(i),
     $           REAL(delta(ising,i)),AIMAG(delta(ising,i)),
     $           absdelta(ising,i)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      theta=(/(i,i=0,mthsurf)/)/REAL(mthsurf,r8)
      CALL ascii_open(out_unit,"ipout_jumpsvd_p"//spolo//"_t"
     $     //storo//"_r"//sresp//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_JUMP SVD: "
      WRITE(out_unit,*)"rsing: ",rsing
      WRITE(out_unit,*)"mpert: ",mpert
      
      DO ising=1,rsing
         WRITE(out_unit,'(4x,I2,a20,es12.3)')ising,
     $        "th singular values:",s(ising)
         WRITE(out_unit,'(2x,a24)')"left singular vector"
         DO i=1,rsing
            WRITE(out_unit,'(2x,I12,2(2x,e12.3))')ising,
     $           REAL(u(i,ising)),AIMAG(u(i,ising))
         ENDDO
         WRITE(out_unit,'(2x,a24)')"right singular vector"
         DO i=1,mpert
            WRITE(out_unit,'(2x,I12,2(2x,e12.3))')mfac(i),
     $           REAL(vt(ising,i)),AIMAG(vt(ising,i))
         ENDDO
         WRITE(out_unit,'(2x,a24)')"right singular function"
         DO i=0,mthsurf
            WRITE(out_unit,'(3(2x,e12.3))')theta(i),
     $           REAL(vtfun(ising,i)),AIMAG(vtfun(ising,i))
         ENDDO
         WRITE(out_unit,'(2x,a24)')"right singular response"
         DO i=0,mthsurf
            WRITE(out_unit,'(3(2x,e12.3))')theta(i),
     $           REAL(pvtfun(ising,i)),AIMAG(pvtfun(ising,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_jumpcoup
c-----------------------------------------------------------------------
c     subprogram 12. ipout_rawdata
c     reponse to the rawdata
c     __________________________________________________________________
c     infile   : file name containing rawdata
c     polo     : 0: polar angle
c                1: hamada
c                2: pest
c                3: equal arc length
c                4: boozer
c     toro     : 0: polar angle
c                1: hamada
c     boutmn   : normal field in the used coordinate
c-----------------------------------------------------------------------
      SUBROUTINE ipout_rawdata(surfmode,infile,polo,toro,boutmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: polo,toro,surfmode
      CHARACTER(128), INTENT(IN) :: infile
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: boutmn

      INTEGER :: inputst,mstart,mtot,i1,i,j,mnum,nnum
      REAL(r8) :: r1,r2
      CHARACTER(1) :: spolo,storo

      REAL(r8), DIMENSION(0:mthsurf) :: theta
      COMPLEX(r8), DIMENSION(mpert) :: binmn,finmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun

      REAL(r8), DIMENSION(:,:), POINTER :: cosmn,sinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn

      theta=twopi*(/(i,i=0,mthsurf)/)/REAL(mthsurf,r8)  
c-----------------------------------------------------------------------
c     read data from file by jon.
c-----------------------------------------------------------------------   
c      CALL ascii_open(in_unit,infile,"OLD")
c      i=0
c      DO
c         i=i+1
c         READ(in_unit,'(2x,I3,1x,es15.8,1x,es15.8)',IOSTAT=inputst)
c     $        i1,r1,r2
c         IF (inputst<0) EXIT
c      ENDDO
c      mtot=i-1
c      ALLOCATE(rawmn(mtot,nn))
c      CALL ascii_close(in_unit)

c      CALL ascii_open(in_unit,infile,"OLD")
c      i=0
c      DO
c         i=i+1
c         READ(in_unit,'(2x,I3,1x,es15.8,1x,es15.8)',IOSTAT=inputst)
c     $        i1,r1,r2
c         IF (i.EQ.1) mstart=i1
c         rawmn(i,nn)=r1+ifac*r2
c         IF (inputst<0) EXIT         
c      ENDDO
c      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     read data from file for by michael.
c-----------------------------------------------------------------------
c      IF (surfmode .EQ. 100) THEN
c         mnum=64
c         nnum=64
c         ALLOCATE(cosmn(-mnum:mnum,0:nnum),sinmn(-mnum:mnum,0:nnum),
c     $        rawmn(-mnum:mnum,0:nnum))
c         CALL ascii_open(in_unit,infile,"old")
c 1111    FORMAT(1x,25f12.6)

c         DO i=-mnum,mnum
c            READ(in_unit,1111) (cosmn(-i,j),j=0,nnum)
c            READ(in_unit,1111) (sinmn(-i,j),j=0,nnum)
c         ENDDO
c         CALL ascii_close(in_unit)
c         rawmn=(cosmn+ifac*sinmn)*gauss
c         binmn=rawmn(mlow:mhigh,nn)
c      ELSE
c         binmn=0
c         binmn(surfmode-mlow+1)=1.0*gauss
c      ENDIF  
c-----------------------------------------------------------------------
c     read data from file for by myself.
c-----------------------------------------------------------------------
      IF (surfmode .EQ. 100) THEN
         mnum=30
         nnum=1
         ALLOCATE(cosmn(-mnum:mnum,nnum),sinmn(-mnum:mnum,nnum),
     $        rawmn(-mnum:mnum,nnum))
         CALL ascii_open(in_unit,infile,"old")
 1111    FORMAT(2x,I3,2(1x,e15.8))
         
         DO i=1,nnum
            DO j=-mnum,mnum
               READ(in_unit,1111) i1,cosmn(j,i),sinmn(j,i)
            ENDDO
         ENDDO
         sinmn=-sinmn
         CALL ascii_close(in_unit)
         rawmn=cosmn+ifac*sinmn
         binmn=rawmn(mlow:mhigh,1)
      ELSE
         binmn=0
         binmn(surfmode-mlow+1)=1.0*gauss
      ENDIF         
c-----------------------------------------------------------------------
c     get the plasma response on the control surface.
c-----------------------------------------------------------------------
      finmn=binmn
      CALL ipeq_tobasis(psilim,finmn,polo,toro,1,1)
      boutmn=MATMUL(permeabmats,finmn)
      CALL ipeq_tocoord(psilim,boutmn,polo,toro,1,1)
      CALL iscdftb(mfac,mpert,binfun,mthsurf,binmn)
      CALL iscdftb(mfac,mpert,boutfun,mthsurf,boutmn)  
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------   
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro  
      CALL ascii_open(out_unit,"ipout_rawdata_p"//spolo//"_t"//
     $     storo//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_RAWDATA: plasma response for an external 
     $perturbation on the control surface"
      WRITE(out_unit,*)"poloidal coordinate:",polo
      WRITE(out_unit,*)"toroidal coordinate:",toro
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(5(2x,a12))')"m",
     $     "rebin","imbin","rebout","imbout"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I12,4(2x,es12.3))')mfac(i),
     $        REAL(binmn(i)),AIMAG(binmn(i)),
     $        REAL(boutmn(i)),AIMAG(boutmn(i))
      ENDDO
      WRITE(out_unit,*)"FUNCTIONS"
      WRITE(out_unit,'(5(2x,a12))')"theta",
     $     "rebin","imbin","rebout","imbout"
      DO i=0,mthsurf
         WRITE(out_unit,'(5(2x,es12.3))')theta(i),
     $        REAL(binfun(i)),AIMAG(binfun(i)),
     $        REAL(boutfun(i)),AIMAG(boutfun(i))
      ENDDO
      CALL ascii_close(out_unit)      
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_rawdata
c-----------------------------------------------------------------------
c     subprogram 10. ipout_dataman.
c     manipulate error field data.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_dataman(outfile1,outfile2,outfile,factor)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: outfile1,outfile2,outfile
      REAL(r8), INTENT(IN) :: factor

      INTEGER :: i,j,mnum,nnum,ms

      REAL(r8), DIMENSION(:,:), POINTER :: cosmn,sinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn1,rawmn2,rawmn

      mnum=30
      nnum=1
      ALLOCATE(cosmn(-mnum:mnum,nnum),sinmn(-mnum:mnum,nnum),
     $     rawmn1(-mnum:mnum,nnum),rawmn2(-mnum:mnum,nnum),
     $     rawmn(-mnum:mnum,nnum))
 1002 FORMAT(2x,I3,2(1x,e15.8))

      CALL ascii_open(in_unit,outfile1,"old")
      DO i=1,nnum
         DO j=-mnum,mnum
            READ(in_unit,1002)ms,cosmn(j,i),sinmn(j,i)
         ENDDO
      ENDDO
      CALL ascii_close(in_unit)
      rawmn1=cosmn+ifac*sinmn

      CALL ascii_open(in_unit,outfile2,"old")
      DO i=1,nnum
         DO j=-mnum,mnum
            READ(in_unit,1002)ms,cosmn(j,i),sinmn(j,i)
         ENDDO
      ENDDO
      CALL ascii_close(in_unit)
      rawmn2=cosmn+ifac*sinmn

      rawmn=rawmn1+factor*(rawmn2-rawmn1)

      CALL ascii_open(out_unit,outfile,"UNKNOWN")
      DO i=1,nnum
         DO j=-mnum,mnum
            WRITE(out_unit,1002)j,REAL(rawmn(j,i)),AIMAG(rawmn(j,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_dataman
c-----------------------------------------------------------------------
c     subprogram 10. ipout_contra.
c     plot contravariant components of xi and b.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xipsi components on the control surface
c     rstep  : radial points to be observed
c-----------------------------------------------------------------------
      SUBROUTINE ipout_contra(egnum,xwpimn,rstep)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rstep
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn
      
      INTEGER :: istep,ipert

      REAL(r8), DIMENSION(:), POINTER :: psis,qs
      COMPLEX(r8), DIMENSION(:,:), POINTER :: xwpmns,bwpmns

      WRITE(*,*)"calculating contravariant components"

      ALLOCATE(psis(0:rstep),qs(0:rstep))
      ALLOCATE(xwpmns(mpert,0:rstep),bwpmns(mpert,0:rstep))
      IF (rstep .EQ. mstep) THEN
         psis=psifac
         qs=qfac
      ELSE
         psis=(/(istep,istep=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
         DO istep=0,rstep
            CALL spline_eval(sq,psis(istep),0)
            qs(istep)=sq%f(4)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     compute contravariant components.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      
      CALL ipeq_alloc
      DO istep=0,rstep
         CALL ipeq_contra(psis(istep))
         xwpmns(:,istep)=xwp_mn
         bwpmns(:,istep)=bwp_mn
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write data
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_contra.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CONTRA: "//
     $     "write contravariant components"
      WRITE(out_unit,'(6(2x,a12))')"psi","q","re(xwp_mn)","im(xwp_mn)",
     $     "re(bwp_mn)","im(bwp_mn)"
      DO ipert=1,mpert
         WRITE(out_unit, '(2x,a3,2x,I3)')"m=",mfac(ipert)
         DO istep=0,rstep
            WRITE(out_unit,'(6(2x,es12.3))')
     $           psis(istep),qs(istep),
     $           REAL(xwpmns(ipert,istep)),AIMAG(xwpmns(ipert,istep)),
     $           REAL(bwpmns(ipert,istep)),AIMAG(bwpmns(ipert,istep))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_contra
c-----------------------------------------------------------------------
c     subprogram 9. ipout_eulerbst.
c     compute strength of eulerian perturbed b field.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xipsi components on the control surface
c     rstep  : radial points to be observed
c-----------------------------------------------------------------------
      SUBROUTINE ipout_eulerbst(egnum,xwpimn,rstep)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rstep
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn
      
      INTEGER :: istep
      REAL(r8), DIMENSION(0:mthsurf) :: eulerbsqf,unitf
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,bwt_fun,bwz_fun,
     $     bvp_fun,bvt_fun,bvz_fun

      REAL(r8), DIMENSION(:), POINTER :: psis,rfacs,eulerbst,area

      WRITE(*,*)"calculating eulerian bfield strength"

      eulerbsqf=0
      unitf=1.0
      ALLOCATE(psis(0:rstep),rfacs(0:rstep),
     $     eulerbst(0:rstep),area(0:rstep))
      IF (rstep .EQ. mstep) THEN
         psis=psifac
      ELSE
         psis=(/(istep,istep=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      
      CALL ipeq_alloc
      DO istep=0,rstep
         CALL bicube_eval(rzphi,psis(istep),REAL(0,r8),0)
         rfacs(istep)=SQRT(rzphi%f(1))
         
         CALL ipeq_contra(psis(istep))
         CALL ipeq_cova(psis(istep))
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
         CALL iscdftb(mfac,mpert,bwz_fun,mthsurf,bwz_mn)
         CALL iscdftb(mfac,mpert,bvp_fun,mthsurf,bvp_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         eulerbsqf=REAL(CONJG(bwp_fun)*bvp_fun+
     $        CONJG(bwt_fun)*bvt_fun+CONJG(bwz_fun)*bvz_fun)
         eulerbst(istep)=
     $        issurfint(lmfac,lmpert,eulerbsqf,mthsurf,psis(istep))/2.0
         area(istep)=issurfint(lmfac,lmpert,unitf,mthsurf,psis(istep))
      ENDDO
      CALL ipeq_dealloc
      eulerbst=SQRT(ABS(eulerbst))/area
c-----------------------------------------------------------------------
c     write data
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_eulerbst.out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_EULERBST: "//
     $     "eulerian b field strength on flux surfaces"
      WRITE(out_unit,'(3(2x,a12))')"psi","rfac(0)","eulebst"
      DO istep=0,rstep
         WRITE(out_unit,'(3(2x,es12.3))')
     $        REAL(psis(istep),4),REAL(rfacs(istep),4),
     $        REAL(eulerbst(istep),4)
      ENDDO
      CALL ascii_close(out_unit)   
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_eulerbst

      END MODULE ipout_mod

c      CALL ipout_resp
c      CALL ipout_surfun(1,10,2,0)
c      CALL ipout_surfmode(0,26,2,0)
c-----------------------------------------------------------------------
c     simulations with various features.
c-----------------------------------------------------------------------
c      infile="nstx_ef.dat "
c      CALL ipout_rawdata(100,infile,1,1,edge_mn)
c      CALL ipeq_bntox(psilim,edge_mn)
c      edge_flag=.TRUE.
c      CALL ipout_peq(0,edge_mn,1)
c      CALL ipout_contour(0,edge_mn,100,1)
c      CALL ipout_delta(0,edge_mn,5,100)
c      eps=0.001
c      CALL ipout_jumpcoup(eps,5,1,1,1) 
c      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     analysis for DIII-D data.
c-----------------------------------------------------------------------
c      infile="surfmn.out.cs.erri"
c      CALL ipout_rawdata(100,infile,2,0,edge_mn)
c      CALL ipeq_tobasis(psilim,edge_mn,2,0,1,0)
c      CALL ipeq_bntox(psilim,edge_mn)
c      edge_flag=.TRUE.
c      CALL ipout_ortho(0,edge_mn,2,4,100,0)
c      CALL ipout_bstrength(0,edge_mn,100,0)
c      CALL ipout_peq(0,edge_mn,1)
c      CALL ipout_contour(0,edge_mn,100,1)
c      CALL ipout_delta(0,edge_mn,3,1000)
c      eps=0.001
c      CALL ipout_jumpcoup(eps,3,2,0,1)
c      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     analysis for NSTX data.
c-----------------------------------------------------------------------
c      ALLOCATE(edge_fun(0:mthsurf),edge_mn(mpert))
c      infile="nstx_ef.dat"
c      CALL ipout_rawdata(2,infile,1,1,edge_mn)
c      CALL ipeq_tobasis(psilim,edge_mn,1,1,1,0)
c      CALL ipeq_bntox(psilim,edge_mn)
c      edge_flag=.TRUE.
c      CALL ipout_bstrength(0,edge_mn,1000,0)
c      CALL ipout_peq(0,edge_mn,1)
c      CALL ipout_contour(0,edge_mn,100,1)
c      CALL ipout_delta(0,edge_mn,5,1000)
c      eps=0.001
c      CALL ipout_jumpcoup(eps,2,1,1,1)

c-----------------------------------------------------------------------
c     temporal routines for error field control.
c-----------------------------------------------------------------------
      sfactor=""
      ssfactor="p01"
      infile1="nstx_ef_uncorrected.dat"
      infile2="nstx_ef_corrected.dat"
c      DO in=1,21
c         factor=-1.0+0.1*REAL((in-1),r8)
c         WRITE(UNIT=sfactor,FMT='(I2)')NINT((ABS(10*factor)))
c         IF (factor >= 0) THEN
c            ssfactor="p"//sfactor
c         ELSE
c            ssfactor="m"//sfactor
c         ENDIF
c         
c         CALL ipout_dataman(infile1,infile2,infile,factor)
c         CALL ipout_errfld(infile,errtype,edge_mn,polo,toro,resp,
c     $        bnomn,xwpmn,ssfactor)
c         CALL ipout_singfld(0,xwpmn,dist,msing,ssfactor)
c      ENDDO


c-----------------------------------------------------------------------
c     create 3d surface data for controlling distributions.
c-----------------------------------------------------------------------
c      CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
c      WRITE(out_unit,'(1p,5e16.8)')rs
c      WRITE(out_unit,'(1p,5e16.8)')zs
c      WRITE(out_unit,'(1p,5e16.8)')REAL(hextfun(1,:))
c      WRITE(out_unit,'(1p,5e16.8)')AIMAG(hextfun(1,:))      
c      CALL ascii_close(out_unit)

c      filein="iptemp.txt"
c      file1="ipidl_3dsurf_cextfld.out"
c      ddist=0.1
c      nnn=nn
c      mmtheta=mthnumb
c      np=72
c      nt=72
c      CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1)   

c-----------------------------------------------------------------------
c     subprogram 2. ipout_singcoup.
c     compute coupling between singular surfaces and external fields.
c     __________________________________________________________________
c     dist     : measuring distance from rational surfaces.
c     lsing    : lowest singular surfaces to be analyzed.
c     rsing    : number of rational surfaces to be analyzed.
c     coordinates used in analysis.
c     polo     : 0: polar angle
c                1: hamada
c                2: pest
c                3: equal arc length
c                4: boozer
c     toro     : 0: polar angle
c                1: hamada
c     meas     : 1: deltas
c                2: singcurs
c                3: total resonant field
c                4: hwidth of magnetic island
c     mthnumb  : number of boundary points
c     __________________________________________________________________
c     bextmn   : most singular external distributions.
c     cextmn   : external distributions for an individual island.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singcoup(dist,lsing,nsing,rsing,polo,toro,meas,
     $     mthnumb,bextmn,cextmn,svdvec_flag,d3_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: lsing,nsing,rsing,polo,toro,meas,mthnumb
      REAL(r8), INTENT(IN) :: dist
      LOGICAL, INTENT(IN) :: svdvec_flag,d3_flag
      COMPLEX(r8), DIMENSION(rsing,mpert), INTENT(OUT) :: bextmn,cextmn

      INTEGER :: i,j,itheta,ising,ithnum,resnum,lwork,hfsurf,rank,d3sing
      INTEGER(r8) :: nnn,mmtheta,np,nt
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,htheta,rcond,
     $     rfacmax,rbextmax,rcextmax,ddist,shear,weighfac
      COMPLEX(r8) :: lnbwp1mn,rnbwp1mn,lcorrec,rcorrec,phasefac
      CHARACTER(1) :: spolo,storo,sresp,smeas,ssing,sising
      CHARACTER(2) :: sssing,ssising
      CHARACTER(72) :: filein,file1

      REAL(r8), DIMENSION(rsing) :: j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs
      REAL(r8), DIMENSION(0:mthnumb) :: rs,zs,mthang,delpsib,
     $     norvec,nozvec
      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,lcormn,rcormn,
     $     fkaxmn,singflx_mn,singbno_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,jcfun,lcorfun,rcorfun

      COMPLEX(r8), DIMENSION(rsing,mpert) :: deltas,delcurs,corcurs,
     $     singcurs,singbnoflds,islandhwids,measures

      REAL(r8), DIMENSION(rsing) :: s,cs
      COMPLEX(r8), DIMENSION(rsing,rsing) :: u
      COMPLEX(r8), DIMENSION(rsing,mpert) :: a,btotmn,rbextmn,rcextmn
      COMPLEX(r8), DIMENSION(mpert,rsing) :: b
      COMPLEX(r8), DIMENSION(mpert,mpert) :: vt
      COMPLEX(r8), DIMENSION(rsing,0:mthsurf) :: bextfun,btotfun,cextfun
      COMPLEX(r8), DIMENSION(rsing,0:mthnumb) :: bno_rvc,bno_zvc,
     $     rbextfun,rcextfun

      REAL(r8), DIMENSION(5*mpert) :: rwork
      REAL(r8), DIMENSION(5*rsing-4) :: rwork2
      COMPLEX(r8), DIMENSION(3*rsing+mpert) :: work
      COMPLEX(r8), DIMENSION(2*rsing+mpert+1) :: work2      
c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing singular coupling with external perturbation"

      DO ising=1,rsing
         resnum=NINT(singtype(lsing+ising-1)%q)*nn-mlow+1
         respsi=singtype(lsing+ising-1)%psifac
         CALL spline_eval(sq,respsi,0)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            rfacs(itheta)=rfac
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            sqreqb(itheta)=(sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $           /(twopi*rfac)**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
         ENDDO
         j_c(ising)=1.0/issurfint(lmfac,lmpert,jcfun,mthsurf,respsi)*
     $        (chi1*sq%f(4))**2/mu0
      ENDDO
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      deltas=0
      delcurs=0
      corcurs=0
      singcurs=0
      CALL ipeq_alloc
      DO i=1,mpert
         WRITE(*,'(a12,I3,a24)')"computing m=",mfac(i),
     $        "poloidal mode coupling"
         binmn=0
         binmn(i)=1.0*gauss
         CALL ipeq_cotoha(psilim,binmn,polo,toro)
         CALL ipeq_weight(psilim,binmn,1)
         boutmn=MATMUL(permeabmats(modelnum,:,:),binmn)
         CALL ipeq_weight(psilim,boutmn,0)
         CALL ipeq_bntoxp(psilim,boutmn,xwp_mn)
         edge_flag=.TRUE.
         CALL idcon_build(0,xwp_mn)
c-----------------------------------------------------------------------
c     evaluate delta and singular currents.
c-----------------------------------------------------------------------
         DO ising=1,rsing
            resnum=NINT(singtype(lsing+ising-1)%q)*nn-mlow+1
            respsi=singtype(lsing+ising-1)%psifac
            lpsi=respsi-dist
            CALL ipeq_contra(lpsi)
            lnbwp1mn=nbwp1_mn(resnum)
            CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
            CALL spline_eval(sq,lpsi,0)
            shear=mfac(resnum)*sq%f1(4)/sq%f(4)**2
            DO itheta=0,mthsurf
               CALL bicube_eval(rzphi,lpsi,theta(itheta),1)
               rfac=SQRT(rzphi%f(1))
               eta=twopi*(theta(itheta)+rzphi%f(2))
               r(itheta)=ro+rfac*COS(eta)
               z(itheta)=zo+rfac*SIN(eta)
               jac=rzphi%f(4)
               w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
               w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
               w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
               w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
               w(3,1)=rzphi%fx(3)*2*rfac
               w(3,2)=rzphi%fy(3)/(twopi*rfac)
               sqrpsi=w(1,1)**2+w(1,2)**2
               correc=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $              (sqrpsi*sq%f(4)*chi1)
               lcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
            
            rpsi=respsi+dist
            CALL ipeq_contra(rpsi)
            rnbwp1mn=nbwp1_mn(resnum)
            CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
            CALL spline_eval(sq,rpsi,0)
            DO itheta=0,mthsurf
               CALL bicube_eval(rzphi,rpsi,theta(itheta),1)
               rfac=SQRT(rzphi%f(1))
               eta=twopi*(theta(itheta)+rzphi%f(2))
               r(itheta)=ro+rfac*COS(eta)
               z(itheta)=zo+rfac*SIN(eta)
               jac=rzphi%f(4)
               w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
               w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
               w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
               w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
               w(3,1)=rzphi%fx(3)*2*rfac
               w(3,2)=rzphi%fy(3)/(twopi*rfac)
               sqrpsi=w(1,1)**2+w(1,2)**2
               correc=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $              (sqrpsi*sq%f(4)*chi1)
               rcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
            deltas(ising,i)=rnbwp1mn-lnbwp1mn
            delcurs(ising,i)=-ifac/mfac(resnum)*deltas(ising,i)
            corcurs(ising,i)=rcormn(resnum)-lcormn(resnum)
            singcurs(ising,i)=j_c(ising)*
     $         (delcurs(ising,i)-corcurs(ising,i))
            fkaxmn=0
            fkaxmn(resnum)=-singcurs(ising,i)*chi1/(twopi*ifac*nn)

            IF (meas == 3) THEN         
               ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
               CALL ipvacuum_flxsurf(respsi)
               singflx_mn=MATMUL(fsurf_indmats,fkaxmn)
               singbno_mn=singflx_mn
               CALL ipeq_weight(respsi,singbno_mn,0)
               DEALLOCATE(fsurf_indev,fsurf_indmats)
               singbnoflds(ising,i)=singbno_mn(resnum)
            ENDIF
            
            IF (meas == 4) THEN
               ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
               CALL ipvacuum_flxsurf(respsi)
               singflx_mn=MATMUL(fsurf_indmats,fkaxmn)
               DEALLOCATE(fsurf_indev,fsurf_indmats)
               islandhwids(ising,i)=
     $              SQRT(ABS(4*singflx_mn(resnum)/(shear*sq%f(4)*chi1)))
            ENDIF
         ENDDO
      ENDDO
      CALL ipeq_dealloc

      SELECT CASE(meas)
      CASE(1)
         measures=deltas
      CASE(2)
         measures=singcurs
      CASE(3)
         measures=singbnoflds
      CASE(4)
         measures=islandhwids
      END SELECT
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=smeas, FMT='(I1)')meas
      IF (nsing < 10) THEN
         WRITE(UNIT=ssing, FMT='(I1)')nsing
         sssing="0"//ssing
      ELSE
         WRITE(UNIT=sssing, FMT='(I2)')nsing
      ENDIF

      CALL ascii_open(out_unit,"ipout_coup_p"//spolo//"_t"
     $     //storo//"_m"//smeas//"_s"//sssing//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_COUP: "//
     $     "coupling between delta/singcurs/bnoflds and "//
     $     "external perturbation"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
      WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mpert: ",mpert
      WRITE(out_unit,'(2x,a4,3(2x,a12))')"mfac",
     $     "real","imaginary","absolute"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,a4,f12.3)')"q =",singtype(lsing+ising-1)%q 
         DO i=1,mpert
            WRITE(out_unit,'(2x,I4,3(2x,e12.3))')mfac(i),
     $           REAL(measures(ising,i)),AIMAG(measures(ising,i)),
     $           ABS(measures(ising,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     find the singular distributions and write the results.
c-----------------------------------------------------------------------
      a = measures
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=5*rsing+mpert      
      CALL zgesvd('S','S',rsing,mpert,a,rsing,s,u,rsing,vt,mpert,
     $     work,lwork,rwork,info)
      vt=CONJG(vt)
      DO ising=1,rsing
         bextmn(ising,:)=vt(ising,:) ! right singular vector
         CALL ipeq_cotoha(psilim,vt(ising,:),polo,toro)
         CALL ipeq_weight(psilim,vt(ising,:),1)
         btotmn(ising,:)=MATMUL(permeabmats(modelnum,:,:),vt(ising,:))
         CALL ipeq_weight(psilim,btotmn(ising,:),0)
         CALL ipeq_hatoco(psilim,btotmn(ising,:),polo,toro)
         CALL iscdftb(mfac,mpert,bextfun(ising,:),mthsurf,
     $        bextmn(ising,:))
         CALL iscdftb(mfac,mpert,btotfun(ising,:),mthsurf,
     $        btotmn(ising,:))
      ENDDO   
c-----------------------------------------------------------------------
c     compute the distribution to control indivisual magnetic islands.
c-----------------------------------------------------------------------
      lwork=2*rsing+mpert+1
      b=0
      DO ising=1,rsing
         b(ising,ising)=1.0
      ENDDO
      work2=0
      rwork2=0
      cs=0
      a=measures

      CALL zgelss(rsing,mpert,rsing,a,rsing,b,mpert,cs,rcond,rank,work2,
     $        lwork,rwork2,info) ! minimize square norms
      DO ising=1,rsing
         cextmn(ising,:)=b(:,ising)
         CALL iscdftb(mfac,mpert,cextfun(ising,:),mthsurf,
     $        cextmn(ising,:))
      ENDDO
c-----------------------------------------------------------------------
c     compute the plasma boundary information.
c-----------------------------------------------------------------------
      mthang=(/(ithnum,ithnum=0,mthnumb)/)/REAL(mthnumb,r8)
      DO ithnum=0,mthnumb
         CALL bicube_eval(rzphi,psilim,mthang(ithnum),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(mthang(ithnum)+rzphi%f(2))
         rs(ithnum)=ro+rfac*COS(eta)
         zs(ithnum)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(ithnum)/jac
         w(1,2)=-rzphi%fy(1)*pi*rs(ithnum)/(rfac*jac)
         delpsib(ithnum)=SQRT(w(1,1)**2+w(1,2)**2)
         norvec(ithnum)=(cos(eta)*w(1,1)-sin(eta)*w(1,2))/
     $        delpsib(ithnum)
         nozvec(ithnum)=(sin(eta)*w(1,1)+cos(eta)*w(1,2))/
     $        delpsib(ithnum)
      ENDDO
      rfacmax=MAXVAL(rfacs)
c-----------------------------------------------------------------------
c     plot the singular distributions on the plasma boundary.
c-----------------------------------------------------------------------
      DO ising=1,rsing
         rbextmn(ising,:)=bextmn(ising,:)
c-----------------------------------------------------------------------
c     do not change toroidal phi!
c-----------------------------------------------------------------------
         IF (toro == 0) THEN
            CALL ipeq_cotoha(psilim,rbextmn(ising,:),polo,1)
         ELSEIF ((polo == 1).and.(toro == 1)) THEN 
            CALL ipeq_hatoco(psilim,rbextmn(ising,:),1,0)
         ELSE
            CALL ipec_stop("sorry, the used coordinates have not"//
     $           " been implemented.")
         ENDIF
         CALL iscdftb(mfac,mpert,rbextfun(ising,:),mthnumb,
     $        rbextmn(ising,:))
c-----------------------------------------------------------------------
c     normalize with phase.
c-----------------------------------------------------------------------
         phasefac=EXP(ifac*ATAN(AIMAG(rbextfun(ising,0))/
     $        REAL(rbextfun(ising,0))))
         weighfac=rfacmax/MAXVAL(ABS(rbextfun(ising,:)))*0.2
         rbextfun(ising,:)=rbextfun(ising,:)*weighfac*phasefac
         bno_rvc(ising,:)=rbextfun(ising,:)*norvec
         bno_zvc(ising,:)=rbextfun(ising,:)*nozvec
      ENDDO

      CALL ascii_open(out_unit,"ipout_bextfld_m"//smeas//"_s"//sssing//
     $     "_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_BEXTFLD: "//
     $     "singular external field distribution on the boundary"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
      WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mthnumb: ",mthnumb
      WRITE(out_unit,'(6(2x,a12))')"r","z","real dr","real dz",
     $     "imag dr","imag dz"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,I4,a24,e12.3,a8,f12.3)')ising,
     $        "th singular values:",s(ising)," for q =",
     $        singtype(lsing+ising-1)%q        
         DO ithnum=0,mthnumb
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(ithnum),zs(ithnum),
     $           REAL(bno_rvc(ising,ithnum)),
     $           REAL(bno_zvc(ising,ithnum)),
     $           AIMAG(bno_rvc(ising,ithnum)),
     $           AIMAG(bno_zvc(ising,ithnum))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     plot the controlling distributions on the plasma boundary.
c-----------------------------------------------------------------------
      DO ising=1,rsing
         rcextmn(ising,:)=cextmn(ising,:)
c-----------------------------------------------------------------------
c     do not change toroidal phi!
c-----------------------------------------------------------------------
         IF (toro == 0) THEN
            CALL ipeq_cotoha(psilim,rcextmn(ising,:),polo,1)
         ELSEIF ((polo == 1).and.(toro == 1)) THEN 
            CALL ipeq_hatoco(psilim,rcextmn(ising,:),1,0)
         ELSE
            CALL ipec_stop("sorry, the used coordinates have not"//
     $           " been implemented.")
         ENDIF
         CALL iscdftb(mfac,mpert,rcextfun(ising,:),mthnumb,
     $        rcextmn(ising,:))
c-----------------------------------------------------------------------
c     normalize with phase.
c-----------------------------------------------------------------------
         phasefac=EXP(ifac*ATAN(AIMAG(rcextfun(ising,0))/
     $        REAL(rcextfun(ising,0))))
         weighfac=rfacmax/MAXVAL(ABS(rcextfun(ising,:)))*0.2
         rcextfun(ising,:)=rcextfun(ising,:)*weighfac*phasefac
         bno_rvc(ising,:)=rcextfun(ising,:)*norvec
         bno_zvc(ising,:)=rcextfun(ising,:)*nozvec
      ENDDO

      CALL ascii_open(out_unit,"ipout_cextfld_m"//smeas//"_s"//sssing//
     $     "_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CEXTFLD: "//
     $     "external field controlling each magnetic island"//
     $     " on the boundary"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
      WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mthnumb: ",mthnumb
      WRITE(out_unit,'(6(2x,a12))')"r","z","real dr","real dz",
     $     "imag dr","imag dz"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,I4,a24,e12.3,a8,f12.3)')ising,
     $        "th control values:",s(ising)," at q =",
     $        singtype(lsing+ising-1)%q     
         DO ithnum=0,mthnumb
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(ithnum),zs(ithnum),
     $           REAL(bno_rvc(ising,ithnum)),
     $           REAL(bno_zvc(ising,ithnum)),
     $           AIMAG(bno_rvc(ising,ithnum)),
     $           AIMAG(bno_zvc(ising,ithnum))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     write svd vector and functions.
c-----------------------------------------------------------------------
      IF (svdvec_flag) THEN
      
         CALL ascii_open(out_unit,"ipout_coupsvd_exvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_EXVEC: "//
     $        "right external singular vector"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert: ",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(bextmn(ising,i)),AIMAG(bextmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_exfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_EXFUN: "//
     $        "right external singular function"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf: ",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta","real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(bextfun(ising,j)),AIMAG(bextfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_ttvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_TTVEC: "//
     $        "right total singular vector"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert: ",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(btotmn(ising,i)),AIMAG(btotmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_ttfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_TTFUN: "//
     $        "right total singular function"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf: ",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(btotfun(ising,j)),AIMAG(btotfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupcvd_exvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPCVD_EXVEC: "//
     $        "external field vectors to control each magnetic island"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert: ",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a38)')ising,
     $           "th magnetic-island-control external vector"
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(cextmn(ising,i)),AIMAG(cextmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupcvd_exfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_s"//sssing//"_n"
     $        //sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPCVD_EXFUN: "//
     $        "external field functions to control each magnetic island"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance: ",dist
         WRITE(out_unit,'(2x,a12,2x,f12.3)')"lsing: ",singtype(lsing)%q 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing: ",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf: ",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a38)')ising,
     $           "th magnetic-island-control external function"
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(cextfun(ising,j)),AIMAG(cextfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write data for 3d surface plot.
c-----------------------------------------------------------------------
      IF (d3_flag) THEN
         d3sing=2
         DO ising=lsing,d3sing
            IF (ising < 10) THEN
               WRITE(UNIT=sising, FMT='(I1)')ising
               ssising="0"//sising
            ELSE
               WRITE(UNIT=ssising, FMT='(I2)')ising    
            ENDIF
            CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
            WRITE(out_unit,'(1p,5e16.8)')rs
            WRITE(out_unit,'(1p,5e16.8)')zs
            WRITE(out_unit,'(1p,5e16.8)')REAL(rbextfun(ising,:))
            WRITE(out_unit,'(1p,5e16.8)')AIMAG(rbextfun(ising,:))      
            CALL ascii_close(out_unit)
         
            filein="iptemp.txt"
            file1="ipidl_3dsurf_bextfld_m"//smeas//"_s"//sssing//
     $           "_i"//ssising//"_n"//sn//".out"
            ddist=0.3
            nnn=nn
            mmtheta=mthnumb
            np=72
            nt=72*nnn
            CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
         ENDDO
         DO ising=lsing,d3sing
            IF (ising < 10) THEN
               WRITE(UNIT=sising, FMT='(I1)')ising
               ssising="0"//sising
            ELSE
               WRITE(UNIT=ssising, FMT='(I2)')ising    
            ENDIF
            CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
            WRITE(out_unit,'(1p,5e16.8)')rs
            WRITE(out_unit,'(1p,5e16.8)')zs
            WRITE(out_unit,'(1p,5e16.8)')REAL(rcextfun(ising,:))
            WRITE(out_unit,'(1p,5e16.8)')AIMAG(rcextfun(ising,:))      
            CALL ascii_close(out_unit)
         
            filein="iptemp.txt"
            file1="ipidl_3dsurf_cextfld_m"//smeas//"_s"//sssing//
     $            "_i"//ssising//"_n"//sn//".out"
            ddist=0.3
            nnn=nn
            mmtheta=mthnumb
            np=72
            nt=72*nnn
            CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     normalize output.
c-----------------------------------------------------------------------
      bextmn=bextmn*gauss
      cextmn=cextmn*gauss
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singcoup
c-----------------------------------------------------------------------
c     subprogram 3. ipeq_orth.
c     compute orthonormal components for xi and b.
c     need to run ipeq_contra and ipeq_cova before running.
c-----------------------------------------------------------------------
      SUBROUTINE ipeq_orth(psi,tag_flag,vec_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: psi
      LOGICAL, INTENT(IN) :: tag_flag,vec_flag
   
      INTEGER :: itheta

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,tangth,normth
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,xwt_fun,
     $     bwt_fun,xvz_fun,bvz_fun,xno_fun,bno_fun,xta_fun,bta_fun,
     $     xpa_fun,bpa_fun
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
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
         IF (vec_flag) THEN
            no_rvec(itheta)=(cos(eta)*w(1,1)-sin(eta)*w(1,2))/
     $           delpsi(itheta)
            no_zvec(itheta)=(sin(eta)*w(1,1)+cos(eta)*w(1,2))/
     $           delpsi(itheta)
         ENDIF
         IF (tag_flag) THEN
            w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
            w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
            v(3,3)=twopi*r(itheta)    
            tangth(itheta)=jac*delpsi(itheta)**2
            normth(itheta)=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2))
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
      CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
      CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
      xno_fun=xwp_fun/delpsi
      bno_fun=bwp_fun/delpsi
      CALL iscdftf(mfac,mpert,xno_fun,mthsurf,xno_mn)
      CALL iscdftf(mfac,mpert,bno_fun,mthsurf,bno_mn)
      IF (tag_flag) THEN
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         xta_fun=(xwt_fun*tangth-xwp_fun*normth)/(delpsi*v(3,3))
         bta_fun=(bwt_fun*tangth-bwp_fun*normth)/(delpsi*v(3,3))
         xpa_fun=xvz_fun/v(3,3)
         bpa_fun=bvz_fun/v(3,3)
         CALL iscdftf(mfac,mpert,xta_fun,mthsurf,xta_mn)
         CALL iscdftf(mfac,mpert,bta_fun,mthsurf,bta_mn)
         CALL iscdftf(mfac,mpert,xpa_fun,mthsurf,xpa_mn)
         CALL iscdftf(mfac,mpert,bpa_fun,mthsurf,bpa_mn)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipeq_orth

c-----------------------------------------------------------------------
c     subprogram 2. ipout_singcoup.
c     compute coupling between singular surfaces and external fields.
c     __________________________________________________________________
c     dist     : measuring distance from rational surfaces.
c     lsing    : lowest singular surfaces to be analyzed.
c     rsing    : number of rational surfaces to be analyzed.
c     polo     : poloidal angle
c     toro     : toroidal angle
c     mthnumb  : number of boundary points
c     __________________________________________________________________
c     bextmn   : most singular external distributions.
c     cextmn   : external distributions for an individual island.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singcoup(dist,portion,polo,toro,meas,
     $     mthnumb,obexmn,ocexmn,svdvec_flag,d3_flag,osing,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: polo,toro,meas,mthnumb,osing,labl
      REAL(r8), INTENT(IN) :: dist,lportion,hportion
      LOGICAL, INTENT(IN) :: svdvec_flag,d3_flag
      COMPLEX(r8), DIMENSION(osing,mpert), INTENT(OUT) :: obexmn,ocexmn

      INTEGER :: i,j,itheta,ising,ithnum,resnum,lwork,hfsurf,rank,
     $     rsing,rpert,a0,sing1,sing2
      INTEGER(r8) :: nnn,mmtheta,np,nt
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,htheta,rcond,
     $     rfacmax,rbextmax,rbtotmax,rcextmax,ddist,shear,weighfac,a1,a2
      COMPLEX(r8) :: lnbwp1mn,rnbwp1mn,lcorrec,rcorrec,phasefac
      CHARACTER(1) :: spolo,storo,sresp,smeas,ssing,sising,slabl
      CHARACTER(2) :: sssing,ssising
      CHARACTER(72) :: filein,file1

      REAL(r8), DIMENSION(5*mpert) :: rwork
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun
      REAL(r8), DIMENSION(0:mthnumb) :: rs,zs,mthang,delpsib,
     $     norvec,nozvec
      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,lcormn,rcormn,
     $     fkaxmn,singflx_mn,singbno_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: vt

      REAL(r8), DIMENSION(:), POINTER :: j_c,s,cs,rwork2
      COMPLEX(r8), DIMENSION(:), POINTER :: work,work2
      COMPLEX(r8), DIMENSION(:,:), POINTER :: deltas,delcurs,corcurs,
     $     singcurs,singbnoflds,islandhwids,measures,bextmn,cextmn,u,
     $     a,btotmn,rbextmn,rbtotmn,rcextmn,b,bextfun,btotfun,cextfun,
     $     bnorvc,bnozvc,rbextfun,rbtotfun,rcextfun
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: fsurfindmats
c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing singular coupling with external perturbation"

      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=smeas, FMT='(I1)')meas
      WRITE(UNIT=slabl, FMT='(I1)')labl

      IF (meas < 5) THEN
         i=0
         DO ising=1,msing
            IF ((singtype(ising)%psifac > lportion) .AND.
     $           (singtype(ising)%psifac < hportion)) THEN 
               i=i+1
               IF (i==1) sing1=ising
            ENDIF
         ENDDO
         rsing=i
         sing2=sing1+rsing-1

         ALLOCATE(j_c(rsing),s(rsing),cs(rsing),rwork2(5*rsing-4),
     $        work(3*rsing+mpert),work2(2*rsing+mpert+1),
     $        deltas(rsing,mpert),delcurs(rsing,mpert),
     $        corcurs(rsing,mpert),singcurs(rsing,mpert),
     $        singbnoflds(rsing,mpert),islandhwids(rsing,mpert),
     $        measures(rsing,mpert),cextmn(rsing,mpert),
     $        bextmn(rsing,mpert),btotmn(rsing,mpert),
     $        u(rsing,rsing),b(mpert,rsing),a(rsing,mpert),
     $        rbextmn(rsing,mpert),rbtotmn(rsing,mpert),
     $        rcextmn(rsing,mpert),fsurfindmats(rsing,mpert,mpert),
     $        bextfun(rsing,0:mthsurf),btotfun(rsing,0:mthsurf),
     $        cextfun(rsing,0:mthsurf),
     $        bnorvc(rsing,0:mthnumb),bnozvc(rsing,0:mthnumb),
     $        rbextfun(rsing,0:mthnumb),rbtotfun(rsing,0:mthnumb),
     $        rcextfun(rsing,0:mthnumb))  
    
         WRITE(*,*)"computing surface currents and inductances"
         DO ising=sing1,sing2
            respsi=singtype(ising)%psifac
            CALL spline_eval(sq,respsi,0)
            DO itheta=0,mthsurf
               CALL bicube_eval(rzphi,respsi,theta(itheta),1)
               rfac=SQRT(rzphi%f(1))
               rfacs(itheta)=rfac
               eta=twopi*(theta(itheta)+rzphi%f(2))
               r(itheta)=ro+rfac*COS(eta)
               z(itheta)=zo+rfac*SIN(eta)
               jac=rzphi%f(4)
               w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
               w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
               delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
               sqreqb(itheta)=(sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $              /(twopi*r(itheta))**2
               jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
            ENDDO
            j_c(ising-sing1+1)=1.0/issurfint(jcfun,mthsurf,respsi,0,0)*
     $           (chi1*sq%f(4))**2/mu0
            IF ((meas == 3) .OR. (meas == 4)) THEN
               CALL ipeq_alloc
               ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
               CALL ipvacuum_flxsurf(respsi)
               fsurfindmats(ising-sing1+1,:,:)=fsurf_indmats
               DEALLOCATE(fsurf_indev,fsurf_indmats)
               CALL ipeq_dealloc
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
         deltas=0
         delcurs=0
         corcurs=0
         singcurs=0
         CALL ipeq_alloc
         DO i=1,mpert
            WRITE(*,'(a12,I3,a24)')"computing m=",mfac(i),
     $           "poloidal mode coupling"
            binmn=0
            binmn(i)=1.0
            CALL ipeq_cotoha(psilim,binmn,polo,toro)
            CALL ipeq_weight(psilim,binmn,1)
            boutmn=MATMUL(permeabmats(modelnum,:,:),binmn)
            CALL ipeq_weight(psilim,boutmn,0)
            CALL ipeq_bntoxp(psilim,boutmn,edge_mn)
            edge_flag=.TRUE.
            CALL idcon_build(0,edge_mn)
c-----------------------------------------------------------------------
c     evaluate delta/singular current/normal field/islands.
c-----------------------------------------------------------------------
            DO ising=sing1,sing2
               resnum=NINT(singtype(ising)%q*nn)-mlow+1
               respsi=singtype(ising)%psifac
               lpsi=respsi-dist/(nn*ABS(singtype(ising)%q1))
               CALL ipeq_contra(lpsi)
               lnbwp1mn=nbwp1_mn(resnum)
               CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
               CALL spline_eval(sq,lpsi,0)
               shear=mfac(resnum)*sq%f1(4)/sq%f(4)**2
               DO itheta=0,mthsurf
                  CALL bicube_eval(rzphi,lpsi,theta(itheta),1)
                  rfac=SQRT(rzphi%f(1))
                  eta=twopi*(theta(itheta)+rzphi%f(2))
                  r(itheta)=ro+rfac*COS(eta)
                  z(itheta)=zo+rfac*SIN(eta)
                  jac=rzphi%f(4)
                  w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
                  w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
                  w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
                  w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
                  w(3,1)=rzphi%fx(3)*2*rfac
                  w(3,2)=rzphi%fy(3)/(twopi*rfac)
                  sqrpsi=w(1,1)**2+w(1,2)**2
                  correc=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $                 (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $                 (sqrpsi*sq%f(4)*chi1)
                  lcorfun(itheta)=bwp_fun(itheta)*correc
               ENDDO
               CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
               
               rpsi=respsi+dist/(nn*ABS(singtype(ising)%q1)) 
               CALL ipeq_contra(rpsi)
               rnbwp1mn=nbwp1_mn(resnum)
               CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
               CALL spline_eval(sq,rpsi,0)
               DO itheta=0,mthsurf
                  CALL bicube_eval(rzphi,rpsi,theta(itheta),1)
                  rfac=SQRT(rzphi%f(1))
                  eta=twopi*(theta(itheta)+rzphi%f(2))
                  r(itheta)=ro+rfac*COS(eta)
                  z(itheta)=zo+rfac*SIN(eta)
                  jac=rzphi%f(4)
                  w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
                  w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
                  w(2,1)=-rzphi%fx(2)*twopi**2*r(itheta)*rfac/jac
                  w(2,2)=rzphi%fx(1)*pi*r(itheta)/(rfac*jac)
                  w(3,1)=rzphi%fx(3)*2*rfac
                  w(3,2)=rzphi%fy(3)/(twopi*rfac)
                  sqrpsi=w(1,1)**2+w(1,2)**2
                  correc=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $                 (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $                 (sqrpsi*sq%f(4)*chi1)
                  rcorfun(itheta)=bwp_fun(itheta)*correc
               ENDDO
               CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
               deltas(ising-sing1+1,i)=(rnbwp1mn-lnbwp1mn)/twopi
               delcurs(ising-sing1+1,i)=
     $              -ifac/mfac(resnum)*deltas(ising-sing1+1,i)
               corcurs(ising-sing1+1,i)=(rcormn(resnum)-lcormn(resnum))/
     $              twopi
               singcurs(ising-sing1+1,i)=j_c(ising-sing1+1)*
     $              (delcurs(ising-sing1+1,i)-corcurs(ising-sing1+1,i))
               fkaxmn=0
               fkaxmn(resnum)=-singcurs(ising-sing1+1,i)*
     $              chi1/(twopi*ifac*nn)

               IF (meas == 1) THEN
                  measures(ising-sing1+1,i)=deltas(ising-sing1+1,i)
               ELSE IF (meas == 2) THEN
                  measures(ising-sing1+1,i)=singcurs(ising-sing1+1,i)
               ELSE IF (meas == 3) THEN         
                  singflx_mn=
     $                 MATMUL(fsurfindmats(ising-sing1+1,:,:),fkaxmn)
                  singbno_mn=singflx_mn
                  CALL ipeq_weight(respsi,singbno_mn,0)
                  singbnoflds(ising-sing1+1,i)=singbno_mn(resnum)
                  measures(ising-sing1+1,i)=singbnoflds(ising-sing1+1,i)
               ELSE IF (meas == 4) THEN
                  singflx_mn=
     $                 MATMUL(fsurfindmats(ising-sing1+1,:,:),fkaxmn)
                  islandhwids(ising-sing1+1,i)= 
     $                 SQRT(ABS(4*singflx_mn(resnum)/
     $                 (shear*sq%f(4)*chi1)))
                  measures(ising-sing1+1,i)=4*singflx_mn(resnum)/
     $                 (shear*sq%f(4)*chi1)
               ENDIF
            ENDDO
         ENDDO
         CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
         CALL ascii_open(out_unit,"ipout_coup_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUP: "//
     $        "coupling between delta/singcurs/bnoflds and "//
     $        "external perturbation"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion 
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert:",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac","real","imaginary"
         DO ising=sing1,sing2
            WRITE(out_unit,'(2x,a4,f12.3)')"q =",singtype(ising)%q 
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(measures(ising-sing1+1,i)),
     $              AIMAG(measures(ising-sing1+1,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     write matrix.
c-----------------------------------------------------------------------
         CALL ascii_open(out_unit,"ipec_singcoup_matrix_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_n"//sn//".out","UNKNOWN")
         WRITE(out_unit,'(1x,4(a8,I4),2(a8,e16.8))')"rsing:",rsing,
     $        "mpert:",mpert,"mlow:",mlow,"mhigh:",mhigh,
     $        "psilim:",psilim,"qlim:",qlim
         WRITE(out_unit,'(1x,200(e16.8))')
     $        (singtype(ising)%psifac,singtype(ising)%q,
     $        ising=1,rsing)
         DO i=1,mpert
            WRITE(out_unit,'(1x,200(e16.8))')
     $           (REAL(measures(ising,i)),AIMAG(measures(ising,i)),
     $           ising=1,rsing)
         ENDDO
         CALL ascii_close(out_unit)
         DEALLOCATE(j_c,deltas,delcurs,corcurs,singcurs,singbnoflds,
     $        islandhwids)
      ELSE
         
         CALL ascii_open(in_unit,"measures_n"//sn//".in","old")
         READ(in_unit,'(I12)')rsing
         READ(in_unit,'(I12)')rpert
         
         ALLOCATE(s(rsing),cs(rsing),rwork2(5*rsing-4),
     $        work(3*rsing+mpert),work2(2*rsing+mpert+1),
     $        measures(rsing,mpert),cextmn(rsing,mpert),
     $        bextmn(rsing,mpert),btotmn(rsing,mpert),
     $        u(rsing,rsing),b(mpert,rsing),a(rsing,mpert),
     $        rbextmn(rsing,mpert),rbtotmn(rsing,mpert),
     $        rcextmn(rsing,mpert),
     $        bextfun(rsing,0:mthsurf),btotfun(rsing,0:mthsurf),
     $        cextfun(rsing,0:mthsurf),
     $        bnorvc(rsing,0:mthnumb),bnozvc(rsing,0:mthnumb),
     $        rbextfun(rsing,0:mthnumb),rbtotfun(rsing,0:mthnumb),
     $        rcextfun(rsing,0:mthnumb))  
         
         measures=0
         DO ising=1,rsing
            DO i=1,rpert
               READ(in_unit,'(I4,2e12.3)')a0,a1,a2
               IF ((a0 .GE. mlow) .AND. (a0 .LE. mhigh)) THEN
                  measures(ising,i)=a1+ifac*a2
               ENDIF
            ENDDO
         ENDDO
         CALL ascii_close(in_unit)
      ENDIF
c-----------------------------------------------------------------------
c     find the singular distributions and write the results.
c-----------------------------------------------------------------------
      a = measures
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=3*rsing+mpert      
      CALL zgesvd('S','S',rsing,mpert,a,rsing,s,u,rsing,vt,mpert,
     $     work,lwork,rwork,info)
      vt=CONJG(vt)
      DO ising=1,rsing
         bextmn(ising,:)=vt(ising,:) ! right singular vector
         CALL ipeq_cotoha(psilim,vt(ising,:),polo,toro)
         CALL ipeq_weight(psilim,vt(ising,:),1)
         btotmn(ising,:)=MATMUL(permeabmats(modelnum,:,:),vt(ising,:))
         CALL ipeq_weight(psilim,btotmn(ising,:),0)
         CALL ipeq_hatoco(psilim,btotmn(ising,:),polo,toro)
         CALL iscdftb(mfac,mpert,bextfun(ising,:),mthsurf,
     $        bextmn(ising,:))
         CALL iscdftb(mfac,mpert,btotfun(ising,:),mthsurf,
     $        btotmn(ising,:))
      ENDDO   
c-----------------------------------------------------------------------
c     compute the distribution to control indivisual magnetic islands.
c-----------------------------------------------------------------------
      lwork=2*rsing+mpert+1
      b=0
      DO ising=1,rsing
         b(ising,ising)=1.0*gauss
      ENDDO
      work2=0
      rwork2=0
      cs=0
      a=measures

      CALL zgelss(rsing,mpert,rsing,a,rsing,b,mpert,cs,rcond,rank,work2,
     $     lwork,rwork2,info)   ! minimize square norms
      DO ising=1,rsing
         cextmn(ising,:)=b(:,ising)
         CALL iscdftb(mfac,mpert,cextfun(ising,:),mthsurf,
     $        cextmn(ising,:))
      ENDDO
c-----------------------------------------------------------------------
c     compute the plasma boundary information.
c-----------------------------------------------------------------------
      mthang=(/(ithnum,ithnum=0,mthnumb)/)/REAL(mthnumb,r8)
      DO ithnum=0,mthnumb
         CALL bicube_eval(rzphi,psilim,mthang(ithnum),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(mthang(ithnum)+rzphi%f(2))
         rs(ithnum)=ro+rfac*COS(eta)
         zs(ithnum)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(ithnum)/jac
         w(1,2)=-rzphi%fy(1)*pi*rs(ithnum)/(rfac*jac)
         delpsib(ithnum)=SQRT(w(1,1)**2+w(1,2)**2)
         norvec(ithnum)=(cos(eta)*w(1,1)-sin(eta)*w(1,2))/
     $        delpsib(ithnum)
         nozvec(ithnum)=(sin(eta)*w(1,1)+cos(eta)*w(1,2))/
     $        delpsib(ithnum)
      ENDDO
      rfacmax=MAXVAL(rfacs)
c-----------------------------------------------------------------------
c     plot the singular distributions on the plasma boundary.
c-----------------------------------------------------------------------
      DO ising=1,rsing
         rbextmn(ising,:)=bextmn(ising,:)
c-----------------------------------------------------------------------
c     change coordinates to polar phi with hamada theta.
c-----------------------------------------------------------------------
         CALL ipeq_cotoha(psilim,rbextmn(ising,:),polo,toro)
         CALL ipeq_hatoco(psilim,rbextmn(ising,:),1,0)
         CALL iscdftb(mfac,mpert,rbextfun(ising,:),mthnumb,
     $        rbextmn(ising,:))
c-----------------------------------------------------------------------
c     normalize with phase.
c-----------------------------------------------------------------------
         phasefac=EXP(ifac*ATAN2(AIMAG(rbextfun(ising,0)),
     $        REAL(rbextfun(ising,0))))
         weighfac=rfacmax/MAXVAL(ABS(rbextfun(ising,:)))*0.2
         rbextfun(ising,:)=rbextfun(ising,:)*weighfac*phasefac
         bnorvc(ising,:)=rbextfun(ising,:)*norvec
         bnozvc(ising,:)=rbextfun(ising,:)*nozvec
      ENDDO
c-----------------------------------------------------------------------
c     write the result.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_bextfld_m"//smeas//
     $     "_l"//slabl//"_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_BEXTFLD: "//
     $     "singular external field distribution on the boundary"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mthnumb:",mthnumb
      WRITE(out_unit,'(6(2x,a12))')"r","z","real dr","real dz",
     $     "imag dr","imag dz"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $        "th singular values:",s(ising)        
         DO ithnum=0,mthnumb
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(ithnum),zs(ithnum),
     $           REAL(bnorvc(ising,ithnum)),
     $           REAL(bnozvc(ising,ithnum)),
     $           AIMAG(bnorvc(ising,ithnum)),
     $           AIMAG(bnozvc(ising,ithnum))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     plot the total distributions on the plasma boundary.
c-----------------------------------------------------------------------
      DO ising=1,rsing
         rbtotmn(ising,:)=btotmn(ising,:)
c-----------------------------------------------------------------------
c     change coordinates to polar phi with hamada theta.
c-----------------------------------------------------------------------
         CALL ipeq_cotoha(psilim,rbtotmn(ising,:),polo,toro)
         CALL ipeq_hatoco(psilim,rbtotmn(ising,:),1,0)
         CALL iscdftb(mfac,mpert,rbtotfun(ising,:),mthnumb,
     $        rbtotmn(ising,:))
c-----------------------------------------------------------------------
c     normalize with phase.
c-----------------------------------------------------------------------
         phasefac=EXP(ifac*ATAN2(AIMAG(rbtotfun(ising,0)),
     $        REAL(rbtotfun(ising,0))))
         weighfac=rfacmax/MAXVAL(ABS(rbtotfun(ising,:)))*0.2
         rbtotfun(ising,:)=rbtotfun(ising,:)*weighfac*phasefac
         bnorvc(ising,:)=rbtotfun(ising,:)*norvec
         bnozvc(ising,:)=rbtotfun(ising,:)*nozvec
      ENDDO
c-----------------------------------------------------------------------
c     write the result.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_btotfld_m"//smeas//
     $     "_l"//slabl//"_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_BTOTFLD: "//
     $     "singular total field distribution on the boundary"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mthnumb:",mthnumb
      WRITE(out_unit,'(6(2x,a12))')"r","z","real dr","real dz",
     $     "imag dr","imag dz"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $        "th singular values:",s(ising)        
         DO ithnum=0,mthnumb
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(ithnum),zs(ithnum),
     $           REAL(bnorvc(ising,ithnum)),
     $           REAL(bnozvc(ising,ithnum)),
     $           AIMAG(bnorvc(ising,ithnum)),
     $           AIMAG(bnozvc(ising,ithnum))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     plot the control distributions on the plasma boundary.
c-----------------------------------------------------------------------
      DO ising=1,rsing
         rcextmn(ising,:)=cextmn(ising,:)
c-----------------------------------------------------------------------
c     change coordinates to polar phi with hamada theta.
c-----------------------------------------------------------------------
         CALL ipeq_cotoha(psilim,rcextmn(ising,:),polo,toro)
         CALL ipeq_hatoco(psilim,rcextmn(ising,:),1,0)
         CALL iscdftb(mfac,mpert,rcextfun(ising,:),mthnumb,
     $        rcextmn(ising,:))
c-----------------------------------------------------------------------
c     normalize with phase.
c-----------------------------------------------------------------------
         phasefac=EXP(ifac*ATAN2(AIMAG(rcextfun(ising,0)),
     $        REAL(rcextfun(ising,0))))
         weighfac=rfacmax/MAXVAL(ABS(rcextfun(ising,:)))*0.2
         rcextfun(ising,:)=rcextfun(ising,:)*weighfac*phasefac
         bnorvc(ising,:)=rcextfun(ising,:)*norvec
         bnozvc(ising,:)=rcextfun(ising,:)*nozvec
      ENDDO
c-----------------------------------------------------------------------
c     write the result.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipout_cextfld_m"//smeas//
     $     "_l"//slabl//"_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_CEXTFLD: "//
     $     "external field controlling each magnetic island"//
     $     " on the boundary"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
      WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
      WRITE(out_unit,'(2x,a12,2x,I12)')"mthnumb:",mthnumb
      WRITE(out_unit,'(6(2x,a12))')"r","z","real dr","real dz",
     $     "imag dr","imag dz"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $        "th control values:",cs(ising)    
         DO ithnum=0,mthnumb
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(ithnum),zs(ithnum),
     $           REAL(bnorvc(ising,ithnum)),
     $           REAL(bnozvc(ising,ithnum)),
     $           AIMAG(bnorvc(ising,ithnum)),
     $           AIMAG(bnozvc(ising,ithnum))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     write svd vector and functions.
c-----------------------------------------------------------------------
      IF (svdvec_flag) THEN
      
         CALL ascii_open(out_unit,"ipout_coupsvd_exvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_EXVEC: "//
     $        "right external singular vector"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert:",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(bextmn(ising,i)),AIMAG(bextmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_exfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_EXFUN: "//
     $        "right external singular function"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf:",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta","real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(bextfun(ising,j)),AIMAG(bextfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_ttvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_TTVEC: "//
     $        "right total singular vector"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert:",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=1,mpert
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(btotmn(ising,i)),AIMAG(btotmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupsvd_ttfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPSVD_TTFUN: "//
     $        "right total singular function"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf:",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a24,e12.3)')ising,
     $           "th singular value:",s(ising)
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(btotfun(ising,j)),AIMAG(btotfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupcvd_exvec_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPCVD_EXVEC: "//
     $        "external field vectors to control each magnetic island"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mpert:",mpert
         WRITE(out_unit,'(2x,a4,2(2x,a12))')"mfac",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a38)')ising,
     $           "th magnetic-island-control external vector"
            DO i=1,mpert   
               WRITE(out_unit,'(2x,I4,2(2x,e12.3))')mfac(i),
     $              REAL(cextmn(ising,i)),AIMAG(cextmn(ising,i))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
         CALL ascii_open(out_unit,"ipout_coupcvd_exfun_p"//spolo//"_t"
     $        //storo//"_m"//smeas//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
         WRITE(out_unit,*)"IPOUT_COUPCVD_EXFUN: "//
     $        "external field functions to control each magnetic island"
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"lportion:",lportion
         WRITE(out_unit,'(2x,a12,2x,e12.3)')"hportion:",hportion 
         WRITE(out_unit,'(2x,a12,2x,I12)')"rsing:",rsing
         WRITE(out_unit,'(2x,a12,2x,I12)')"mthsurf:",mthsurf
         WRITE(out_unit,'(3(2x,a12))')"theta",
     $        "real","imaginary"
         DO ising=1,rsing
            WRITE(out_unit,'(2x,I4,a38)')ising,
     $           "th magnetic-island-control external function"
            DO i=0,mthsurf
               hfsurf=INT(mthsurf/2.0)
               IF (i <= hfsurf) THEN
                  j=i+hfsurf
                  htheta=twopi*theta(j)-twopi
               ELSE
                  j=i-hfsurf
                  htheta=twopi*theta(j)
               ENDIF
               WRITE(out_unit,'(3(2x,e12.3))')htheta,
     $              REAL(cextfun(ising,j)),AIMAG(cextfun(ising,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write data for 3d surface plot.
c-----------------------------------------------------------------------
      IF (d3_flag) THEN
         DO ising=1,osing
            IF (ising < 10) THEN
               WRITE(UNIT=sising, FMT='(I1)')ising
               ssising="0"//sising
            ELSE
               WRITE(UNIT=ssising, FMT='(I2)')ising    
            ENDIF
            CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
            WRITE(out_unit,'(1p,5e16.8)')rs
            WRITE(out_unit,'(1p,5e16.8)')zs
            WRITE(out_unit,'(1p,5e16.8)')REAL(rbextfun(ising,:))
            WRITE(out_unit,'(1p,5e16.8)')AIMAG(rbextfun(ising,:))      
            CALL ascii_close(out_unit)
         
            filein="iptemp.txt"
            file1="ipidl_3dsurf_bextfld_m"//smeas//
     $           "_i"//ssising//"_n"//sn//".out"
            ddist=0.3
            nnn=nn
            mmtheta=mthnumb
            np=72
            nt=72*nnn
            CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
         ENDDO
         DO ising=1,osing
            IF (ising < 10) THEN
               WRITE(UNIT=sising, FMT='(I1)')ising
               ssising="0"//sising
            ELSE
               WRITE(UNIT=ssising, FMT='(I2)')ising    
            ENDIF
            CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
            WRITE(out_unit,'(1p,5e16.8)')rs
            WRITE(out_unit,'(1p,5e16.8)')zs
            WRITE(out_unit,'(1p,5e16.8)')REAL(rbtotfun(ising,:))
            WRITE(out_unit,'(1p,5e16.8)')AIMAG(rbtotfun(ising,:))      
            CALL ascii_close(out_unit)
         
            filein="iptemp.txt"
            file1="ipidl_3dsurf_btotfld_m"//smeas//
     $           "_i"//ssising//"_l"//slabl//"_n"//sn//".out"
            ddist=0.3
            nnn=nn
            mmtheta=mthnumb
            np=72
            nt=72*nnn
            CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
         ENDDO
         DO ising=1,osing 
            IF (ising < 10) THEN
               WRITE(UNIT=sising, FMT='(I1)')ising
               ssising="0"//sising
            ELSE
               WRITE(UNIT=ssising, FMT='(I2)')ising    
            ENDIF
            CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
            WRITE(out_unit,'(1p,5e16.8)')rs
            WRITE(out_unit,'(1p,5e16.8)')zs
            WRITE(out_unit,'(1p,5e16.8)')REAL(rcextfun(ising,:))
            WRITE(out_unit,'(1p,5e16.8)')AIMAG(rcextfun(ising,:))      
            CALL ascii_close(out_unit)
         
            filein="iptemp.txt"
            file1="ipidl_3dsurf_cextfld_m"//smeas//
     $            "_i"//ssising//"_l"//slabl//"_n"//sn//".out"
            ddist=0.3
            nnn=nn
            mmtheta=mthnumb
            np=72
            nt=72*nnn
            CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     normalize output.
c-----------------------------------------------------------------------
      obexmn=bextmn(1:osing,:)*gauss
      ocexmn=cextmn(1:osing,:)*gauss

      DEALLOCATE(s,cs,rwork2,work,work2,
     $     measures,cextmn,bextmn,btotmn,u,b,a,fsurfindmats,
     $     rbextmn,rbtotmn,rcextmn,bextfun,btotfun,cextfun,
     $     bnorvc,bnozvc,rbextfun,rbtotfun,rcextfun)
c-----------------------------------------------------------------------
c     diagnose routines and solutions.
c-----------------------------------------------------------------------
c      svdvec_flag=.TRUE.
c      d3_flag=.FALSE.
c      rdist=0.01
c      bdist=0.01
c      bigstep=0.15
c      smlstep=0.01
c      wegtfac1=0.15
c      wegtfac2=10.0
c      smallwidth=singfac_min
c      rstep1=100
c      rstep2=1000
c      resol=100
c      meas=3
c      mthnumb=500
c      angnum=10
c      osing=1
c      ALLOCATE(obextmn(osing,mpert),ocextmn(osing,mpert))
c      CALL ipdiag_eigen
c      majr=1.0
c      minr=0.2
c      CALL ipdiag_arbsurf(majr,minr)
c      minr=0.1
c      CALL ipdiag_arbsurf(majr,minr)
c      minr=0.01
c      CALL ipdiag_arbsurf(majr,minr)
c      majr=5.0
c      minr=0.1
c      CALL ipdiag_arbsurf(majr,minr)
c      angnum=10
c      rstep1=100
c      CALL ipdiag_angles(angnum,rstep1,3)
c      CALL ipdiag_surfmode(-2,20,1,0)
c      CALL ipdiag_magpot
c      CALL ipdiag_energy
c      CALL ipdiag_respmat
c      CALL ipdiag_energy
c      CALL ipdiag_respmat
c      CALL ipout_errfld(infile,errtype,bermn,2,0,resp,bnomn,xwpmn,0)
c      CALL ipdiag_xicontra(0,xwpmn,0)
c      CALL ipdiag_singcurs(0,xwpmn,rsing,resol,smallwidth)
c      CALL ipdiag_rzpgrid(nr,nz)
c      CALL ipdiag_angles(angnum,rstep1)
c      CALL ipdiag_xbnovc(1,xwpmn)
c      CALL ipdiag_magpot
c      CALL ipdiag_energy
c      CALL ipdiag_respmat
c      edge_flag=.FALSE.
c      CALL ipdiag_xicontra(1,xwpmn,1)
c      CALL ipdiag_xicontra(2,xwpmn,2)
c      CALL ipdiag_xicontra(3,xwpmn,3)
c      CALL ipdiag_singcurs(1,xwpmn,msing,resol,smallwidth)      
c      edge_flag=.TRUE.
c      errtype="nstx2"
c      infile="nstx2_ef.dat"
c      CALL ipdiag_extfld(infile,errtype,bexmn,2,0,1,1,bermn,1)
c      CALL ipdiag_extfld(infile,errtype,bexmn,2,0,3,1,bermn,3)
c      CALL ipdiag_extfld(infile,errtype,bexmn,2,0,4,1,bermn,4) 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singcoup






