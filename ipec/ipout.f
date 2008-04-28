c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     write various output results of ipec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipout_mod
c      1. ipout_resp
c      2. ipout_singcoup
c      3. ipout_pmodcoup
c      4. ipout_ipeccoup
c      5. ipout_errfld
c      6. ipout_singfld
c      7. ipout_pmodb
c      8. ipout_xbnobo
c      9. ipout_xbnorm
c     10. ipout_xbrzphi
c-----------------------------------------------------------------------
c     subprogram 0. ipout_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipout_mod
      USE ipresp_mod
      USE ipvacuum_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipout_resp.
c     write basic information of response.
c     __________________________________________________________________
c     modes     : number of modes, less than mpert, in output.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_resp(modes)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: lwork,modes,i
      REAL(r8), DIMENSION(modes) :: sts
      REAL(r8), DIMENSION(mpert) :: s
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,u,vt

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(4*mpert) :: work

      IF (modes > mpert) modes=mpert
      DO i=1,modes 
         sts(i)=-et(i)/(surfei(i)+surfee(i))
      ENDDO
c-----------------------------------------------------------------------
c     svd analysis for permeability matrix.
c-----------------------------------------------------------------------
      a=permeabmats(modelnum,:,:)
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=4*mpert     
      CALL zgesvd('S','S',mpert,mpert,a,mpert,s,u,mpert,vt,mpert,
     $     work,lwork,rwork,info)

      CALL ascii_open(out_unit,"ipout_resp_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_RESP: response parameters"
      WRITE(out_unit,'(3(2x,a12))')"mode","s_energy","s_permeab"
      DO i=1,modes
         IF (sts(i) > 0) s(i)=-s(i)
         WRITE(out_unit,'(2x,I12,2(2x,e12.3))')i,sts(i),-1/s(i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_resp
c-----------------------------------------------------------------------
c     subprogram 2. ipout_singcoup.
c     compute coupling between singular surfaces and external fields.
c     __________________________________________________________________
c     dist     : measuring distance from rational surfaces
c     polo     : poloidal angle coordinate
c     toro     : toroidal angle coordinate
c     svd_flag : option for svd analysis of resonant normal field
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singcoup(dist,polo,toro)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: polo,toro
      REAL(r8), INTENT(IN) :: dist

      INTEGER :: i,j,itheta,ising,resnum,rsing,rpert,tmlow,tmhigh,tmpert
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,shear
      COMPLEX(r8) :: lnbwp1mn,rnbwp1mn,lcorrec,rcorrec
      CHARACTER(1) :: spolo,storo

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun,wcfun,
     $     dphi,thetas

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,lcormn,rcormn,
     $     fkaxmn,singflx_mn,singbno_mn,ftnmn
      COMPLEX(r8), DIMENSION(lmpert) :: lftnmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun,
     $     ftnfun

      REAL(r8), DIMENSION(msing) :: j_c,w_c

      COMPLEX(r8), DIMENSION(msing,mpert) :: deltas,delcurs,corcurs,
     $     singcurs,singpowers,singbnoflds,islandhwids
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: convmat
      COMPLEX(r8), DIMENSION(msing,mpert,mpert) :: fsurfindmats
      COMPLEX(r8), DIMENSION(:,:), POINTER :: t1mat,t2mat,t3mat,t4mat
      TYPE(spline_type) :: spl 

c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing singular coupling with external perturbation"

      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      
      WRITE(*,*)"computing surface currents and inductances"
      DO ising=1,msing
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
     $           /(twopi*r(itheta))**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
            wcfun(itheta)=SQRT(jac/delpsi(itheta))
         ENDDO
         j_c(ising)=1.0/issurfint(jcfun,mthsurf,respsi,0,0)*
     $        (chi1*sq%f(4))**2/mu0
         w_c(ising)=issurfave(wcfun,mthsurf,respsi)
         
         CALL ipeq_alloc
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
         fsurfindmats(ising,:,:)=fsurf_indmats
         DEALLOCATE(fsurf_indev,fsurf_indmats)
         CALL ipeq_dealloc
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
         binmn(i)=1.0
         CALL ipeq_weight(psilim,binmn,1)
         boutmn=MATMUL(permeabmats(modelnum,:,:),binmn)
         CALL ipeq_weight(psilim,boutmn,0)
         CALL ipeq_bntoxp(psilim,boutmn,edge_mn)
         edge_flag=.TRUE.
         CALL idcon_build(0,edge_mn)
c-----------------------------------------------------------------------
c     evaluate delta/singular current/normal field/islands.
c-----------------------------------------------------------------------
         DO ising=1,msing
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
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $              (sqrpsi*sq%f(4)*chi1)
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
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $              (sqrpsi*sq%f(4)*chi1)
               rcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
            deltas(ising,i)=(rnbwp1mn-lnbwp1mn)/twopi
            delcurs(ising,i)=-ifac/mfac(resnum)*deltas(ising,i)
            corcurs(ising,i)=(rcormn(resnum)-lcormn(resnum))/twopi
            singcurs(ising,i)=j_c(ising)*
     $           (delcurs(ising,i)-corcurs(ising,i))
            singpowers(ising,i)=w_c(ising)*singcurs(ising,i)
            fkaxmn=0
            fkaxmn(resnum)=-singcurs(ising,i)*chi1/(twopi*ifac*nn)

            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
            singbno_mn=singflx_mn
            CALL ipeq_weight(respsi,singbno_mn,0)
            ! change resonant field at desired coordinates
            IF ((polo /= 1) .OR. (toro /= 1)) THEN
               CALL ipeq_hatoco(respsi,singbno_mn,polo,toro)
            ENDIF
            singbnoflds(ising,i)=singbno_mn(resnum)
            islandhwids(ising,i)=4*singflx_mn(resnum)/
     $           (twopi*shear*sq%f(4)*chi1)
         ENDDO
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     convert coordinates for matrix on the plasma boundary.
c-----------------------------------------------------------------------
      IF ((polo /= 1) .OR. (toro /= 1)) THEN
         WRITE(*,*)"convert coordinates"
         CALL spline_alloc(spl,mthsurf,1)
         spl%xs=theta

         CALL spline_eval(sq,psilim,0)
         dphi=0
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psilim,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            bpfac=chi1*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            SELECT CASE(polo)
            CASE(0)
            spl%fs(itheta,1)=fac*bpfac/rfac
            CASE(1)
            spl%fs(itheta,1)=fac
            CASE(2)
            spl%fs(itheta,1)=fac/r(itheta)**2
            CASE(3)
            spl%fs(itheta,1)=fac*bpfac
            CASE(4)
            spl%fs(itheta,1)=fac*bfac**2
            END SELECT
            IF (toro .EQ. 0) THEN
               dphi(itheta)=rzphi%f(3)
            ENDIF
         ENDDO      

         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)
         ! coordinate angle at hamada
         thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
         CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
         DO i=1,lmpert
            lftnmn=0
            lftnmn(i)=1.0
         ! compute given function in hamada angle
            DO itheta=0,mthsurf
               ftnfun(itheta)=0
               ftnfun(itheta)=
     $              lftnmn(i)*EXP(ifac*twopi*lmfac(i)*thetas(itheta))               
            ENDDO
         ! multiply toroidal factor in hamada angle
            IF (toro .EQ. 0) THEN
               ftnfun(:)=ftnfun(:)*EXP(-ifac*nn*dphi(:))
            ELSE
               ftnfun(:)=ftnfun(:)*
     $              EXP(twopi*ifac*nn*sq%f(4)*(thetas(:)-theta(:)))
            ENDIF
            CALL iscdftf(mfac,mpert,ftnfun,mthsurf,ftnmn)
            convmat(:,i)=ftnmn
         ENDDO
         tmlow = lmlow
         tmhigh = lmhigh
         tmpert = lmpert
         ALLOCATE(t1mat(msing,tmpert),t2mat(msing,tmpert),
     $        t3mat(msing,tmpert),t4mat(msing,tmpert))
         t1mat = MATMUL(singbnoflds,convmat)
         t2mat = MATMUL(islandhwids,convmat)
         t3mat = MATMUL(singcurs,convmat)
         t4mat = MATMUL(singpowers,convmat)
      ELSE

         tmlow = mlow
         tmhigh = mhigh
         tmpert = mpert
         ALLOCATE(t1mat(msing,tmpert),t2mat(msing,tmpert),
     $        t3mat(msing,tmpert),t4mat(msing,tmpert))
         t1mat = singbnoflds
         t2mat = islandhwids
         t3mat = singcurs
         t4mat = singpowers
      ENDIF
c-----------------------------------------------------------------------
c     write matrix.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singcoup_matrix_p"//spolo//"_t"
     $     //storo//"_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,'(1x,4(a8,I4),2(a8,e16.8))')"msing:",msing,
     $     "mpert:",tmpert,"mlow:",tmlow,"mhigh:",tmhigh,
     $     "psilim:",psilim,"qlim:",qlim
      WRITE(out_unit,'(1x,400(e16.8))')
     $     (singtype(ising)%psifac,singtype(ising)%q,ising=1,msing)
      DO i=1,tmpert
         WRITE(out_unit,'(1x,400(e16.8))')
     $        (REAL(t1mat(ising,i)),AIMAG(t1mat(ising,i)),ising=1,msing)
      ENDDO
      DO i=1,tmpert
         WRITE(out_unit,'(1x,400(e16.8))')
     $        (REAL(t2mat(ising,i)),AIMAG(t2mat(ising,i)),ising=1,msing)
      ENDDO
      DO i=1,tmpert
         WRITE(out_unit,'(1x,400(e16.8))')
     $        (REAL(t3mat(ising,i)),AIMAG(t3mat(ising,i)),ising=1,msing)
      ENDDO
      DO i=1,tmpert
         WRITE(out_unit,'(1x,400(e16.8))')
     $        (REAL(t4mat(ising,i)),AIMAG(t4mat(ising,i)),ising=1,msing)
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(t1mat,t2mat,t3mat,t4mat)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singcoup
c-----------------------------------------------------------------------
c     subprogram 3. ipout_pmodcoup
c     compute coupling between pmodbst and external field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 4. ipout_ipeccoup
c     compute coupling between perturbed plasma and external field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 5. ipout_errfld
c     error field response on the boundary.
c     __________________________________________________________________
c     infile   : input file name containing rawdata
c     errtype  : input file format
c     binmn    : given error field spectrum
c     polo     : poloidal angle
c     toro     : toroidal angle
c     resp     : 0: actual perturbation
c                1: external perturbation
c     __________________________________________________________________
c     boutmn   : actual normal field in the used coordinate
c     xwpomn   : contravariant normal xi in hamada
c-----------------------------------------------------------------------
      SUBROUTINE ipout_errfld(rinfile,rerrtype,
     $     binmn,polo,toro,resp,boutmn,xwpomn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: polo,toro,resp,labl
      CHARACTER(128), INTENT(IN) :: rinfile,rerrtype
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: binmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: boutmn,xwpomn

      INTEGER :: i,j,mnum,nnum,ms,hfsurf
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: htheta,sengy,pengy,maxtor,
     $     binfunst,boutfunst,bplafunst
      COMPLEX(r8) :: py,sy,mt
      CHARACTER(1) :: spolo,storo,sresp,slabl
      CHARACTER(2) :: slabl2

      REAL(r8), DIMENSION(0:mthsurf) :: 
     $     sbinfun,sboutfun,sbplafun,sbtorfun

      COMPLEX(r8), DIMENSION(mpert) :: ftnmn,foutmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,bplafun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2

      REAL(r8), DIMENSION(:,:), POINTER :: cosmn,sinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn

      IF (edge_flag) THEN
c-----------------------------------------------------------------------
c     read data from file given by d3d, Mike
c-----------------------------------------------------------------------
         IF (rerrtype == "d3d") THEN
            mnum=64
            nnum=64
            ALLOCATE(cosmn(-mnum:mnum,0:nnum),sinmn(-mnum:mnum,0:nnum),
     $           rawmn(-mnum:mnum,0:nnum))
            CALL ascii_open(in_unit,rinfile,"old")
 1000       FORMAT(1x,25f12.6)

            DO i=-mnum,mnum
               READ(in_unit,1000) (cosmn(-i,j),j=0,nnum)
               READ(in_unit,1000) (sinmn(-i,j),j=0,nnum)
            ENDDO
            CALL ascii_close(in_unit)
            rawmn=(cosmn+ifac*sinmn)*gauss
            binmn=rawmn(mlow:mhigh,nn)
            DEALLOCATE(cosmn,sinmn,rawmn)
c-----------------------------------------------------------------------
c     read data from file given by d3d, Ilon.
c-----------------------------------------------------------------------
         ELSE IF (rerrtype == "d3d2") THEN
            mnum=32
            nnum=32
            ALLOCATE(cosmn(-mnum:mnum,0:nnum),sinmn(-mnum:mnum,0:nnum),
     $           rawmn(-mnum:mnum,0:nnum))
            CALL ascii_open(in_unit,rinfile,"old")
 1001       FORMAT(1x,33f12.6)

            DO i=-mnum,mnum
               READ(in_unit,1001) (cosmn(-i,j),j=0,nnum)
               READ(in_unit,1001) (sinmn(-i,j),j=0,nnum)
            ENDDO
            CALL ascii_close(in_unit)
            rawmn=(cosmn+ifac*sinmn)*gauss
            binmn=rawmn(mlow:mhigh,nn)
            DEALLOCATE(cosmn,sinmn,rawmn)
c-----------------------------------------------------------------------
c     read data from file given by nstx.
c-----------------------------------------------------------------------
         ELSE IF (rerrtype == "nstx") THEN
            mnum=60
            nnum=3
            ALLOCATE(cosmn(-mnum:mnum,nnum),sinmn(-mnum:mnum,nnum),
     $           rawmn(-mnum:mnum,nnum))
            CALL ascii_open(in_unit,rinfile,"old")
 1002       FORMAT(1x,I4,2(1x,e15.8))

            DO i=1,nnum
               DO j=-mnum,mnum
                  READ(in_unit,1002)ms,cosmn(j,i),sinmn(j,i)
               ENDDO
            ENDDO
            sinmn=sinmn
            CALL ascii_close(in_unit)
            rawmn=cosmn+ifac*sinmn
            binmn=rawmn(mlow:mhigh,nn)
            DEALLOCATE(cosmn,sinmn,rawmn)
c-----------------------------------------------------------------------
c     read data from file given by nstx.
c-----------------------------------------------------------------------
         ELSE IF (rerrtype == "nstx2") THEN
            mnum=60
            nnum=3
            ALLOCATE(cosmn(-mnum:mnum,nnum),sinmn(-mnum:mnum,nnum),
     $           rawmn(-mnum:mnum,nnum))
            CALL ascii_open(in_unit,rinfile,"old")
 1003       FORMAT(11(1x,e15.8))

            DO i=-mnum,mnum
               READ(in_unit,1003) (cosmn(i,j),j=1,nnum)
               READ(in_unit,1003) (sinmn(i,j),j=1,nnum)
            ENDDO
            CALL ascii_close(in_unit)
            rawmn=cosmn+ifac*sinmn
            binmn=rawmn(mlow:mhigh,nn)
            DEALLOCATE(cosmn,sinmn,rawmn)
c-----------------------------------------------------------------------
c     read data from file given by optima.
c-----------------------------------------------------------------------
         ELSE IF (rerrtype == "optima") THEN
            ALLOCATE(cosmn(mpert,1),sinmn(mpert,1),rawmn(mpert,1))
            CALL ascii_open(in_unit,rinfile,"old")
 1004       FORMAT(1x,I4,2(1x,e15.8))
            DO j=1,mpert
               READ(in_unit,1004)ms,cosmn(j,1),sinmn(j,1)
               IF (ms/=mfac(j)) WRITE(*,*)"different mfac!"
            ENDDO
            CALL ascii_close(in_unit)
            rawmn=cosmn+ifac*sinmn
            binmn=rawmn(:,1)
            DEALLOCATE(cosmn,sinmn,rawmn)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     get the plasma response on the control surface.
c-----------------------------------------------------------------------
      ftnmn=binmn
      CALL ipeq_cotoha(psilim,ftnmn,polo,toro)
      IF (resp == 1) THEN
         CALL ipeq_weight(psilim,ftnmn,1)
         boutmn=MATMUL(permeabmats(modelnum,:,:),ftnmn)
         foutmn=boutmn
         CALL ipeq_weight(psilim,boutmn,0)
         CALL ipeq_bntoxp(psilim,boutmn,xwpomn)
         CALL ipeq_hatoco(psilim,boutmn,polo,toro)
      ELSE
         boutmn=ftnmn
         foutmn=boutmn
         CALL ipeq_bntoxp(psilim,boutmn,xwpomn)
         CALL ipeq_hatoco(psilim,boutmn,polo,toro)
      ENDIF

      CALL iscdftb(mfac,mpert,binfun,mthsurf,binmn)
      CALL iscdftb(mfac,mpert,boutfun,mthsurf,boutmn)
c-----------------------------------------------------------------------
c     compute perturbed energy and maximum torque.
c-----------------------------------------------------------------------
      temp1=0
      work2=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2 = surf_indmats
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
      sy = 0.5*SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))
      sengy = REAL(sy)

      mt = nn*SUM(CONJG(ftnmn)*MATMUL(TRANSPOSE(temp1),foutmn-ftnmn))
      maxtor = ABS(mt)

      temp1=0
      work2=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2 = plas_indmats(modelnum,:,:)
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
      py = 0.5*SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))
      pengy = REAL(py)
c-----------------------------------------------------------------------
c     compute amplifications and maximal torque on the boundary surface.
c-----------------------------------------------------------------------
      bplafun=boutfun-binfun
    
      sbinfun=ABS(binfun)**2/2.0
      sboutfun=ABS(boutfun)**2/2.0
      sbplafun=ABS(bplafun)**2/2.0

      binfunst=SQRT(issurfint(sbinfun,mthsurf,psilim,0,1))
      boutfunst=SQRT(issurfint(sboutfun,mthsurf,psilim,0,1))
      bplafunst=SQRT(issurfint(sbplafun,mthsurf,psilim,0,1))
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=sresp, FMT='(I1)')resp
      IF (labl < 10) THEN 
         WRITE(UNIT=slabl, FMT='(I1)')labl      
         CALL ascii_open(out_unit,"ipout_errfld_p"//spolo//"_t"//
     $        storo//"_r"//sresp//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl
         CALL ascii_open(out_unit,"ipout_errfld_p"//spolo//"_t"//
     $        storo//"_r"//sresp//"_l"//slabl2//"_n"//sn//".out",
     $        "UNKNOWN")
      ENDIF
      WRITE(out_unit,*)"IPOUT_ERRFLD: "//
     $     "plasma response for an external perturbation on the "//
     $     "control surface"
      WRITE(out_unit,'(2x,a12,2x,I6)')"mpert:",mpert
      WRITE(out_unit,'(2x,a12,2x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"perturbed senergy:",sengy
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"perturbed penergy:",pengy
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"maximum torque:",maxtor
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"binfun strength:",binfunst
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"bplafun strength:",bplafunst
      WRITE(out_unit,'(2x,a24,2x,e12.3)')"boutfun strength:",boutfunst
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(2x,a6,4(2x,a12))')"m","rebin","imbin",
     $     "rebout","imbout"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I6,4(2x,e12.3))')mfac(i),
     $        REAL(binmn(i)),AIMAG(binmn(i)),
     $        REAL(boutmn(i)),AIMAG(boutmn(i))
      ENDDO
      WRITE(out_unit,*)"FUNCTIONS"
      WRITE(out_unit,'(5(2x,a12))')"theta",
     $     "rebin","imbin","rebout","imbout"
      DO i=0,mthsurf
         hfsurf=INT(mthsurf/2.0)
         IF (i <= hfsurf) THEN
            j=i+hfsurf
            htheta=twopi*theta(j)-twopi
         ELSE
            j=i-hfsurf
            htheta=twopi*theta(j)
         ENDIF
         WRITE(out_unit,'(5(2x,e12.3))')htheta,
     $        REAL(binfun(j)),AIMAG(binfun(j)),
     $        REAL(boutfun(j)),AIMAG(boutfun(j))
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_errfld
c-----------------------------------------------------------------------
c     subprogram 6. ipout_singfld.
c     compute current and field on rational surfaces.
c     __________________________________________________________________
c     dist   : distance from each rational surface
c     rsing  : number of rational surfaces from the lowest
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singfld(egnum,xwpimn,dist,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      REAL(r8), INTENT(IN) :: dist
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: i,itheta,ising
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,shear,hdist
      COMPLEX(r8) :: lnbwp1mn,rnbwp1mn,lcorrec,rcorrec
      CHARACTER(1) :: slabl
      CHARACTER(2) :: slabl2      

      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun,
     $     tempcos,tempsin,tempfun
      COMPLEX(r8), DIMENSION(mpert) :: lcormn,rcormn,fkaxmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun

      REAL(r8), DIMENSION(msing) :: singpower,singbnofunst,
     $     singmaxtor,island_hwidth,chirikov
      REAL(r8), DIMENSION(0:mthsurf,msing) :: ssingbnofun
      COMPLEX(r8), DIMENSION(msing) :: delta,delcur,corcur,singcur
      COMPLEX(r8), DIMENSION(mpert,msing) :: singbno_mn,singflx_mn
      COMPLEX(r8), DIMENSION(0:mthsurf,msing) :: singbnofun
      
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing delta,singular current and field"
      CALL ipeq_alloc
      CALL idcon_build(egnum,xwpimn)
c-----------------------------------------------------------------------
c     evaluate delta and singular currents.
c-----------------------------------------------------------------------
      DO ising=1,msing
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            sqreqb(itheta)=(sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $           /(twopi*r(itheta))**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
         ENDDO
         j_c(ising)=1.0/issurfint(jcfun,mthsurf,respsi,0,0)*
     $        (chi1*sq%f(4))**2/mu0
         shear=mfac(resnum(ising))*sq%f1(4)/sq%f(4)**2

         lpsi=respsi-dist/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_contra(lpsi)
         lnbwp1mn=nbwp1_mn(resnum(ising))
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL spline_eval(sq,lpsi,0)
         singfac=mfac-nn*sq%f(4)
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
     $           (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $           (sqrpsi*sq%f(4)*chi1)
            lcorfun(itheta)=bwp_fun(itheta)*correc
         ENDDO
         CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
         
         rpsi=respsi+dist/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_contra(rpsi)
         rnbwp1mn=nbwp1_mn(resnum(ising))
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
     $           (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $           (sqrpsi*sq%f(4)*chi1)
            rcorfun(itheta)=bwp_fun(itheta)*correc
         ENDDO
         CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
         delta(ising)=(rnbwp1mn-lnbwp1mn)/twopi
         delcur(ising)=-ifac/mfac(resnum(ising))*delta(ising)*j_c(ising)
         corcur(ising)=(rcormn(resnum(ising))-lcormn(resnum(ising)))*
     $        j_c(ising)/twopi
         singcur(ising)=delcur(ising)-corcur(ising)

         fkaxmn=0
         fkaxmn(resnum(ising))=-singcur(ising)*chi1/(twopi*ifac*nn)
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
         singflx_mn(:,ising)=MATMUL(fsurf_indmats,fkaxmn)
         singbno_mn(:,ising)=singflx_mn(:,ising)
         CALL ipeq_weight(respsi,singbno_mn(:,ising),0)
         DEALLOCATE(fsurf_indmats,fsurf_indev)
c-----------------------------------------------------------------------
c     compute ideal amplification factor.
c----------------------------------------------------------------------- 
         CALL iscdftb(mfac,mpert,singbnofun(:,ising),mthsurf,
     $        singbno_mn(:,ising))
         ssingbnofun(:,ising)=ABS(singbnofun(:,ising))**2/2.0
         singbnofunst(ising)=
     $        SQRT(issurfint(ssingbnofun(:,ising),mthsurf,respsi,0,1))
c-----------------------------------------------------------------------
c     compute half-width of magnetic island.
c-----------------------------------------------------------------------
         island_hwidth(ising)=
     $        SQRT(ABS(4*singflx_mn(resnum(ising),ising)/
     $        (twopi*shear*sq%f(4)*chi1)))
c-----------------------------------------------------------------------
c     compute chirikov paramter.
c-----------------------------------------------------------------------
         IF (ising==1) THEN 
            hdist=(singtype(ising+1)%psifac-respsi)/2.0
         ELSE IF (ising==msing) THEN
            hdist=(respsi-singtype(ising-1)%psifac)/2.0
         ELSE IF ((ising/=1).AND.(ising/=msing)) THEN
            hdist=MIN(singtype(ising+1)%psifac-respsi,
     $           respsi-singtype(ising-1)%psifac)/2.0
         ENDIF
         chirikov(ising)=island_hwidth(ising)/hdist
c-----------------------------------------------------------------------
c     compute normalized power by singular current.
c-----------------------------------------------------------------------
         tempfun=jac*delpsi*
     $        (REAL(singcur(ising))**2+AIMAG(singcur(ising))**2)
         singpower(ising)=issurfave(tempfun,mthsurf,respsi)
c-----------------------------------------------------------------------
c     compute maximum torque by singular current.
c-----------------------------------------------------------------------
         singmaxtor(ising)=nn*ABS(REAL(SUM(CONJG(fkaxmn)*
     $        singflx_mn(:,ising))))
c-----------------------------------------------------------------------
c     change resonant field to desired coordinates.
c-----------------------------------------------------------------------
         IF ((polo /= 1).OR.(toro /= 1)) THEN
            CALL ipeq_hatoco(respsi,singbno_mn(:,ising),polo,toro)
         ENDIF            
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipout_singfld_l"//slabl//
     $        "_n"//sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipout_singfld_l"//slabl2//
     $        "_n"//sn//".out","UNKNOWN")
      ENDIF
      WRITE(out_unit,*)"IPOUT_SINGFLD: "//
     $     "deltas, singular currents and normal fields"
      WRITE(out_unit,'(2x,a12,2x,e12.3)')"distance:",dist
      WRITE(out_unit,'(2x,a12,2x,I12)')"msing:",msing
      WRITE(out_unit,'(2x,a6,11(2x,a12))')"q","psi",
     $     "abs(delta)","abs(delcur)","abs(corcur)","abs(singcur)",
     $     "abs(sgpower)","abs(singbno)","singbnost","singmaxtor",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(2x,f6.3,11(2x,e12.3))')
     $        singtype(ising)%q,singtype(ising)%psifac,
     $        ABS(delta(ising)),ABS(delcur(ising)),ABS(corcur(ising)),
     $        ABS(singcur(ising)),singpower(ising),
     $        ABS(singbno_mn(resnum(ising),ising)),
     $        singbnofunst(ising),singmaxtor(ising),
     $        island_hwidth(ising),chirikov(ising)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singfld
c-----------------------------------------------------------------------
c     subprogram 7. ipout_pmodb.
c     compute perturbed mod b.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xipsi components on the control surface
c     polo   : poloidal angle
c     toro   : toroidal angle
c-----------------------------------------------------------------------
      SUBROUTINE ipout_pmodb(egnum,xwpimn,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert,ising,itheta,dsing,resmn
      REAL(r8) :: sdist,edist,limfac,maxdb
      CHARACTER(1) :: spolo,storo,slabl

      REAL(r8), DIMENSION(0:mthsurf) :: tempeqfun
      COMPLEX(r8), DIMENSION(mpert) :: xwptmn,xwttmn,bvttmn,bvztmn,
     $     eulbpar,lagbpar
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwptfun,xwttfun,
     $     bvttfun,bvztfun

      REAL(r8), DIMENSION(:), POINTER :: psis
      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs
      COMPLEX(r8), DIMENSION(:,:), POINTER :: eulbpar_mn,lagbpar_mn,
     $     eulbparfun,lagbparfun
c-----------------------------------------------------------------------
c     regularize things!!! important consideration!!!
c-----------------------------------------------------------------------
      maxdb=maxdbratio*eqfun%fs(0,0,1)
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing perturbed b field"
      ALLOCATE(psis(rstep),rs(rstep,0:mthsurf),zs(rstep,0:mthsurf),
     $     eulbpar_mn(rstep,mpert),lagbpar_mn(rstep,mpert),
     $     eulbparfun(rstep,0:mthsurf),lagbparfun(rstep,0:mthsurf))

      IF (rstep .EQ. mstep) THEN
         psis=psifac
      ELSE
         psis=(/(istep,istep=1,rstep)/)/REAL(rstep,r8)*(psilim-psilow)
      ENDIF

      CALL idcon_build(egnum,xwpimn)
      
      
      CALL ipeq_alloc
      DO istep=1,rstep
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(ffun,psis(istep),0)
         CALL spline_eval(sq,psis(istep),1)
         CALL ipeq_contra(psis(istep))
         CALL ipeq_cova(psis(istep))
         rs(istep,:)=r
         zs(istep,:)=z
         singfac=mfac-nn*sq%f(4)
         
         xwptmn=xwp_mn
         bvttmn=bvt_mn
         bvztmn=bvz_mn
         DO ipert=1,mpert
            limfac = 1.0
            IF (ABS(singfac(ipert))<bdist) THEN
               limfac = singfac(ipert)/SIGN(bdist,singfac(ipert))
            ENDIF
            xwttmn(ipert)=xwt_mn(ipert)*limfac
         ENDDO
c-----------------------------------------------------------------------
c     compute mod b variations.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwptfun,mthsurf,xwptmn)
         CALL iscdftb(mfac,mpert,xwttfun,mthsurf,xwttmn)
         CALL iscdftb(mfac,mpert,bvttfun,mthsurf,bvttmn)
         CALL iscdftb(mfac,mpert,bvztfun,mthsurf,bvztmn)
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psis(istep),theta(itheta),1)
            eulbparfun(istep,itheta)=
     $           chi1*(bvttfun(itheta)+sq%f(4)*bvztfun(itheta))
     $           /(ffun%f(1)*eqfun%f(1))
            lagbparfun(istep,itheta)=
     $           eulbparfun(istep,itheta)+
     $           xwptfun(itheta)*eqfun%fx(1)+xwttfun(itheta)*eqfun%fy(1)
            IF (abs(eulbparfun(istep,itheta)) > maxdb) 
     $           eulbparfun(istep,itheta) = eulbparfun(istep,itheta)*
     $           maxdb/abs(eulbparfun(istep,itheta))
            IF (abs(lagbparfun(istep,itheta)) > maxdb) 
     $           lagbparfun(istep,itheta) = lagbparfun(istep,itheta)*
     $           maxdb/abs(lagbparfun(istep,itheta))           
         ENDDO
c-----------------------------------------------------------------------
c     decompose components on the given coordinates.
c-----------------------------------------------------------------------
         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
     $        mthsurf,eulbpar_mn(istep,:))
         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
     $        mthsurf,lagbpar_mn(istep,:))
         eulbpar=eulbpar_mn(istep,:)
         lagbpar=lagbpar_mn(istep,:)
         CALL ipeq_hatoco(psis(istep),eulbpar_mn(istep,:),polo,toro)
         CALL ipeq_hatoco(psis(istep),lagbpar_mn(istep,:),polo,toro)
         CALL ipeq_hatoco(psis(istep),eulbpar,1,0)
         CALL ipeq_hatoco(psis(istep),lagbpar,1,0)
         CALL iscdftb(mfac,mpert,eulbparfun(istep,:),mthsurf,eulbpar)
         CALL iscdftb(mfac,mpert,lagbparfun(istep,:),mthsurf,lagbpar)         
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=slabl, FMT='(I1)')labl
      CALL ascii_open(out_unit,"ipout_pmodbmn_p"//spolo//"_t"
     $     //storo//"_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PMODBMN: "//
     $     "components in perturbed mod b on flux surfaces"
      WRITE(out_unit,'(2x,a8,2x,I4)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I4)')"mpert:",mpert
      WRITE(out_unit,'(6(2x,a12))')"psi","mfac",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)"
      DO istep=1,rstep
         DO ipert=1,mpert
            WRITE(out_unit,'(2x,e12.3,2x,I12,4(2x,e12.3))')
     $           psis(istep),mfac(ipert),
     $           REAL(eulbpar_mn(istep,ipert)),
     $           AIMAG(eulbpar_mn(istep,ipert)),
     $           REAL(lagbpar_mn(istep,ipert)),
     $           AIMAG(lagbpar_mn(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_pmodbrz_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_PMODBRZ: "//
     $     "perturbed mod b on a poloidal plane"
      WRITE(out_unit,'(2x,a8,2x,I4)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I4)')"mthsurf:",mthsurf
      WRITE(out_unit,'(6(2x,a12))')"r","z",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(eulbparfun(istep,itheta)),
     $           AIMAG(eulbparfun(istep,itheta)),
     $           REAL(lagbparfun(istep,itheta)),
     $           AIMAG(lagbparfun(istep,itheta))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(psis,rs,zs,eulbpar_mn,lagbpar_mn,
     $     eulbparfun,lagbparfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_pmodb
c-----------------------------------------------------------------------
c     subprogram 8. ipout_xbnorm.
c     compute normal perturbed quantities.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xwp_mn components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbnorm(egnum,xwpimn,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: i,istep,ipert,itheta,rstep
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xno_fun,bno_fun
      CHARACTER(1) :: spolo,storo,slabl

      REAL(r8), DIMENSION(:), POINTER :: psis
      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs
      COMPLEX(r8), DIMENSION(:,:), POINTER :: xnomns,bnomns,
     $     xnofuns,bnofuns

      rstep=100
      ALLOCATE(psis(rstep),rs(rstep,0:mthsurf),zs(rstep,0:mthsurf),
     $     xnomns(rstep,mpert),bnomns(rstep,mpert),
     $     xnofuns(rstep,0:mthsurf),bnofuns(rstep,0:mthsurf))
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
c-----------------------------------------------------------------------
c     if rstep=mstep, use original integration points.
c-----------------------------------------------------------------------
      IF (rstep == mstep) THEN
         psis=psifac
      ELSE
         psis=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing normal quantities"
      CALL ipeq_alloc
      DO istep=1,rstep
         CALL ipeq_contra(psis(istep))
         CALL ipeq_normal(psis(istep))
         rs(istep,:)=r
         zs(istep,:)=z
         xnomns(istep,:)=xno_mn
         bnomns(istep,:)=bno_mn
         CALL ipeq_hatoco(psis(istep),xnomns(istep,:),polo,toro)
         CALL ipeq_hatoco(psis(istep),bnomns(istep,:),polo,toro)
         CALL ipeq_hatoco(psis(istep),xno_mn,1,0)
         CALL ipeq_hatoco(psis(istep),bno_mn,1,0)
         CALL iscdftb(mfac,mpert,xnofuns(istep,:),mthsurf,xno_mn)
         CALL iscdftb(mfac,mpert,bnofuns(istep,:),mthsurf,bno_mn)
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=slabl, FMT='(I1)')labl   
      CALL ascii_open(out_unit,"ipout_xbnormmn_p"//spolo//"_t"//storo
     $     //"_l"//slabl//"_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_XBNORMMN: "//
     $     "components of normal perturbed quantities"
      WRITE(out_unit,'(2x,a8,2x,I4)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I4)')"mpert:",mpert
      WRITE(out_unit,'(6(2x,a12))')"psi","mfac",
     $     "real(xno)","imag(xno)","real(bno)","imag(bno)"
      DO istep=1,rstep
         DO ipert=1,mpert
            WRITE(out_unit,'(2x,e12.3,2x,I12,4(2x,e12.3))')
     $           psis(istep),mfac(ipert),
     $           REAL(xnomns(istep,ipert)),
     $           AIMAG(xnomns(istep,ipert)),
     $           REAL(bnomns(istep,ipert)),
     $           AIMAG(bnomns(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipout_xbnormrz_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_XBNORMRZ: "//
     $     "normal perturbed quantities on a poloidal plane"
      WRITE(out_unit,'(2x,a8,2x,I4)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I4)')"mthsurf:",mthsurf
      WRITE(out_unit,'(6(2x,a12))')"r","z","real(xno)","imag(xno)",
     $     "real(bno)","imag(bno)"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(6(2x,e12.3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(xnofuns(istep,itheta)),
     $           AIMAG(xnofuns(istep,itheta)),
     $           REAL(bnofuns(istep,itheta)),
     $           AIMAG(bnofuns(istep,itheta))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(psis,rs,zs,xnomns,bnomns,xnofuns,bnofuns)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbnorm
c-----------------------------------------------------------------------
c     subprogram 9. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xwp_mn components on the control surface
c     nr     : r grid number
c     nz     : z grid number
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphi(egnum,xwpimn,nr,nz,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: i,j,ipert
      REAL(r8) :: mid,bt0
      CHARACTER(1) :: slabl

      REAL(r8), DIMENSION(0:nr,0:nz) :: ebr,ebz,ebp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: xrr,xrz,xrp,brr,brz,brp
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      WRITE(*,*)"mapping rzphi components on rzphi grid"
      xrr=0
      xrz=0
      xrp=0
      brr=0
      brz=0
      brp=0
      CALL idcon_build(egnum,xwpimn)

      ! evaluate f value for vacuum
      mid = 0.0
      CALL spline_eval(sq,mid,0)
      bt0 = abs(sq%f(1))/(twopi*ro)
          
      CALL ipeq_alloc
      DO i=0,nr
         WRITE(*,*)"finish",i,"th radial mapping"
         DO j=0,nz
            CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),1)
            ebr(i,j) = -psi_in%fy(1)/gdr(i,j)*psio
            ebz(i,j) = psi_in%fx(1)/gdr(i,j)*psio

            IF (gdl(i,j) == 1) THEN
               CALL ipeq_contra(gdpsi(i,j))
               CALL ipeq_cova(gdpsi(i,j))
               CALL ipeq_rzphi(gdpsi(i,j))
               DO ipert=1,mpert
                  xrr(i,j)=xrr(i,j)+xrr_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                  xrz(i,j)=xrz(i,j)+xrz_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                  xrp(i,j)=xrp(i,j)+xrp_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                  brr(i,j)=brr(i,j)+brr_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                  brz(i,j)=brz(i,j)+brz_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                  brp(i,j)=brp(i,j)+brp_mn(ipert)*
     $                 EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
               ENDDO
               xrr(i,j)=
     $              REAL(xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(xrr(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))
               xrz(i,j)=
     $              REAL(xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(xrz(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))
               xrp(i,j)=
     $              REAL(xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(xrp(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))
               brr(i,j)=
     $              REAL(brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(brr(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))
               brz(i,j)=
     $              REAL(brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(brz(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))
               brp(i,j)=
     $              REAL(brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+ifac*
     $              REAL(brp(i,j)*EXP(twopi*ifac*nn*(0.25-gdphi(i,j))))

               CALL spline_eval(sq,gdpsi(i,j),0)
               ebp(i,j) = abs(sq%f(1))/(twopi*gdr(i,j))

            ELSE
               ebp(i,j) = bt0*ro/gdr(i,j)
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl, FMT='(I1)')labl

      CALL ascii_open(out_unit,"ipout_eqbrzphi_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_EQBRZPHI: eq b field in rzphi grid"
      WRITE(out_unit,'(1x,2(a4,I4))')"nr:",nr+1,"nz:",nz+1
      WRITE(out_unit,'(1x,a2,5(a16))')"l","r","z","re(eb_r)",
     $     "re(eb_z)","re(eb_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,5(e16.8))')
     $           gdl(i,j),gdr(i,j),gdz(i,j),
     $           ebr(i,j),ebz(i,j),ebp(i,j)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
   
      CALL ascii_open(out_unit,"ipout_xrzphi_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_XRZPHI: displacements in rzphi grid"
      WRITE(out_unit,'(1x,2(a4,I4))')"nr:",nr+1,"nz:",nz+1
      WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z","re(xi_r)","im(xi_r)",
     $     "re(xi_z)","im(xi_z)","re(xi_phi)","im(xi_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,8(e16.8))')
     $           gdl(i,j),gdr(i,j),gdz(i,j),
     $           REAL(xrr(i,j)),AIMAG(xrr(i,j)),
     $           REAL(xrz(i,j)),AIMAG(xrz(i,j)),
     $           REAL(xrp(i,j)),AIMAG(xrp(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      CALL ascii_open(out_unit,"ipout_brzphi_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPOUT_BRZPHI: perturbed field in rzphi grid"
      WRITE(out_unit,'(1x,2(a4,I4))')"nr:",nr+1,"nz:",nz+1
      WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z","re(b_r)","im(b_r)",
     $     "re(b_z)","im(b_z)","re(b_phi)","im(b_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,8(e16.8))')
     $           gdl(i,j),gdr(i,j),gdz(i,j),
     $           REAL(brr(i,j)),AIMAG(brr(i,j)),
     $           REAL(brz(i,j)),AIMAG(brz(i,j)),
     $           REAL(brp(i,j)),AIMAG(brp(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      RETURN
      END SUBROUTINE ipout_xbrzphi
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
   
      END MODULE ipout_mod
