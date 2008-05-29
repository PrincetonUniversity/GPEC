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
c-----------------------------------------------------------------------
      SUBROUTINE ipout_resp
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: lwork,i,j
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(mpert) :: sts,s
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,b,c,u,vt,temp,eigen,work2

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(4*mpert) :: work
c-----------------------------------------------------------------------
c     svd analysis for permeability matrix.
c-----------------------------------------------------------------------
      a=permeabmats(modelnum,:,:)
      b=a
      c=a
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=4*mpert
      CALL zgesvd('S','S',mpert,mpert,a,mpert,s,u,mpert,vt,mpert,
     $     work,lwork,rwork,info)
      temp=0
      DO i=1,mpert 
         sts(i)=-et(i)/(surfei(i)+surfee(i))
         IF (sts(i) > 0) s(i)=-s(i)
         s(i)=-1/s(i)
         CALL ipeq_xptobn(psilim,wt(:,i),eigen(:,i))
         temp(i,i)=1
      ENDDO
      lwork=2*mpert-1
      CALL zhetrf('L',mpert,b,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,b,mpert,ipiv,temp,mpert,info)
      eigen=MATMUL(temp,eigen)

      CALL ascii_open(out_unit,"ipec_resp_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_RESP: response parameters"
      WRITE(out_unit,'(1x,3(a8,I4),2(a8,e16.8))')
     $     "mpert:",mpert,"mlow:",mlow,"mhigh:",mhigh,
     $     "psilim:",psilim,"qlim:",qlim
      WRITE(out_unit,'(400(1x,e16.8))')(sts(i),i=1,mpert)
      DO j=1,mpert
         WRITE(out_unit,'(400(1x,e16.8))')
     $        (REAL(eigen(j,i)),AIMAG(eigen(j,i)),i=1,mpert)
      ENDDO
      WRITE(out_unit,'(400(1x,e16.8))')(s(i),i=1,mpert)    
      DO j=1,mpert
         WRITE(out_unit,'(400(1x,e16.8))')
     $        (REAL(c(j,i)),AIMAG(c(j,i)),i=1,mpert)
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
c     polo     : poloidal angle coordinate for output
c     toro     : toroidal angle coordinate for output
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
         CALL ipeq_weight(psilim,binmn,mfac,mpert,1)
         boutmn=MATMUL(permeabmats(modelnum,:,:),binmn)
         CALL ipeq_weight(psilim,boutmn,mfac,mpert,0)
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
            delcurs(ising,i)=j_c(ising)*
     $           ifac/mfac(resnum)*deltas(ising,i)
            corcurs(ising,i)=-j_c(ising)*(rcormn(resnum)-lcormn(resnum))
            singcurs(ising,i)=delcurs(ising,i)-corcurs(ising,i)
            singpowers(ising,i)=w_c(ising)*singcurs(ising,i)
            fkaxmn=0
            fkaxmn(resnum)=-singcurs(ising,i)*chi1/(twopi*ifac*nn)

            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
            singbno_mn=singflx_mn
            CALL ipeq_weight(respsi,singbno_mn,mfac,mpert,0)
            ! change resonant field at desired coordinates
            IF ((polo /= 1) .OR. (toro /= 1)) THEN
               CALL ipeq_hatoco(respsi,singbno_mn,mfac,mpert,polo,toro)
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
c     subprogram 3. ipout_errfld
c     error field response on the boundary.
c     __________________________________________________________________
c     infile     : input file name containing rawdata when edge_flag
c     formattype : input file format when edge_flag
c     right      : right-handedness
c     scale      : scale of input
c     binmn      : given error field spectrum without edge_flag
c     polo       : poloidal angle coordinates
c     toro       : toroidal angle coordinates
c     resp       : include plasma response
c     __________________________________________________________________
c     boutmn     : normal field in coordinates
c     xwpomn     : contravariant normal xi in hamada for idcon input
c     __________________________________________________________________
c     labl       : label for multiple run
c-----------------------------------------------------------------------
      SUBROUTINE ipout_errfld(rinfile,formattype,left,scale,
     $     binmn,polo,toro,resp,boutmn,xwpomn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: left,polo,toro,resp,labl
      REAL(r8), INTENT(IN) :: scale
      CHARACTER(128), INTENT(IN) :: rinfile,formattype
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: binmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: boutmn,xwpomn

      INTEGER :: i,j,i1,i2,i3,ms
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: vengy,sengy,pengy,binfunst,boutfunst
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(1) :: spolo,storo,sresp,slabl
      CHARACTER(2) :: slabl2
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: sbinfun,sboutfun

      COMPLEX(r8), DIMENSION(mpert) :: ftnmn,foutmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2

      REAL(r8), DIMENSION(:,:), POINTER :: cosmn,sinmn
      COMPLEX(r8), DIMENSION(:), POINTER :: hawmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn
c-----------------------------------------------------------------------
c     check formattype
c-----------------------------------------------------------------------
 1000 FORMAT(1x,25f12.6)         
 1001 FORMAT(1x,33f12.6)
 1010 FORMAT(11(1x,e15.8))
 1020 FORMAT(1x,I4,2(1x,e15.8))

      IF (edge_flag) THEN
         IF (left == 1) THEN
            i1 = errmmax
            i2 = errmmin
            i3 = -1
         ELSE
            i1 = errmmin
            i2 = errmmax
            i3 = 1
         ENDIF
c-----------------------------------------------------------------------
c     read data.
c-----------------------------------------------------------------------
         ALLOCATE(cosmn(errmmin:errmmax,errnmin:errnmax),
     $        sinmn(errmmin:errmmax,errnmin:errnmax),
     $        rawmn(errmmin:errmmax,errnmin:errnmax),
     $        hawmn(errmmax-errmmin+1))
         CALL ascii_open(in_unit,rinfile,"old")
            
         DO i=i1,i2,i3
            IF (formattype == '1x,25f12.6') THEN
               READ(in_unit,1000) (cosmn(i,j),j=errnmin,errnmax)
               READ(in_unit,1000) (sinmn(i,j),j=errnmin,errnmax)
            ELSE IF (formattype == '1x,33f12.6') THEN
               READ(in_unit,1001) (cosmn(i,j),j=errnmin,errnmax)
               READ(in_unit,1001) (sinmn(i,j),j=errnmin,errnmax)
            ELSE IF (formattype == '11(1x,e15.8)') THEN
               READ(in_unit,1010) (cosmn(i,j),j=errnmin,errnmax)
               READ(in_unit,1010) (sinmn(i,j),j=errnmin,errnmax)
            ELSE IF (formattype == '1x,I4,2(1x,e15.8)') THEN
               READ(in_unit,1020) ms,cosmn(i,1),sinmn(i,1)               
            ELSE
               WRITE(message,'(a)')"can't recognize input format"
               CALL ipec_stop(message)
            ENDIF            
         ENDDO
         CALL ascii_close(in_unit)
         rawmn=cosmn+ifac*sinmn
         hawmn=rawmn(:,nn)
         CALL ipeq_cotoha(psilim,hawmn,lmfac,lmpert,polo,toro)
         binmn=hawmn(mlow-lmlow+1:mhigh-lmlow+1)
         DEALLOCATE(cosmn,sinmn,rawmn,hawmn)
      ENDIF
c-----------------------------------------------------------------------
c     get the plasma response on the control surface.
c-----------------------------------------------------------------------
      IF (scale /= 0) binmn=binmn*scale
      ftnmn = binmn
      IF (resp == 1) THEN
         CALL ipeq_weight(psilim,ftnmn,mfac,mpert,1)
         boutmn=MATMUL(permeabmats(modelnum,:,:),ftnmn)
         foutmn=boutmn
         CALL ipeq_weight(psilim,boutmn,mfac,mpert,0)
         CALL ipeq_bntoxp(psilim,boutmn,xwpomn)
      ELSE
         boutmn=ftnmn
         foutmn=boutmn
         CALL ipeq_bntoxp(psilim,boutmn,xwpomn)
      ENDIF

      CALL iscdftb(mfac,mpert,binfun,mthsurf,binmn)
      CALL iscdftb(mfac,mpert,boutfun,mthsurf,boutmn)
c-----------------------------------------------------------------------
c     compute perturbed energy.
c-----------------------------------------------------------------------
      temp1=0
      work2=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2 = surf_indmats
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)

      vy = 0.5*SUM(CONJG(ftnmn)*MATMUL(temp1,ftnmn))
      sy = 0.5*SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))
      vengy = REAL(vy)
      sengy = REAL(sy)

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
c     compute amplifications on the boundary surface.
c-----------------------------------------------------------------------
      sbinfun=ABS(binfun)**2/2.0
      sboutfun=ABS(boutfun)**2/2.0
      binfunst=SQRT(issurfint(sbinfun,mthsurf,psilim,0,1))
      boutfunst=SQRT(issurfint(sboutfun,mthsurf,psilim,0,1))
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      WRITE(UNIT=sresp, FMT='(I1)')resp
      IF (labl < 10) THEN 
         WRITE(UNIT=slabl, FMT='(I1)')labl      
         CALL ascii_open(out_unit,"ipec_errfld_p"//spolo//"_t"//
     $        storo//"_r"//sresp//"_l"//slabl//"_n"//sn//".out",
     $        "UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl
         CALL ascii_open(out_unit,"ipec_errfld_p"//spolo//"_t"//
     $        storo//"_r"//sresp//"_l"//slabl2//"_n"//sn//".out",
     $        "UNKNOWN")
      ENDIF
      WRITE(out_unit,*)"IPEC_ERRFLD: "//
     $     "plasma response for an external perturbation on the "//
     $     "control surface"
      WRITE(out_unit,'(1x,a12,1x,I6)')"mpert:",mpert
      WRITE(out_unit,'(1x,a12,1x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,'(1x,a24,1x,e16.8)')"perturbed venergy:",vengy
      WRITE(out_unit,'(1x,a24,1x,e16.8)')"perturbed senergy:",sengy
      WRITE(out_unit,'(1x,a24,1x,e16.8)')"perturbed penergy:",pengy
      WRITE(out_unit,'(1x,a24,1x,e16.8)')"perturbed inputst:",binfunst
      WRITE(out_unit,'(1x,a24,1x,e16.8)')"perturbed totalst:",boutfunst
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(1x,a6,4(1x,a16))')"m","rebin","imbin",
     $     "rebout","imbout"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I6,4(1x,e16.8))')mfac(i),
     $        REAL(binmn(i)),AIMAG(binmn(i)),
     $        REAL(boutmn(i)),AIMAG(boutmn(i))
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_errfld
c-----------------------------------------------------------------------
c     subprogram 4. ipout_singfld.
c     compute current and field on rational surfaces.
c     __________________________________________________________________
c     egnum  : eigenmode number without edge_flag
c     xwpimn : edge deformation input when edge_flag
c     dist   : distance from each rational surface
c     polo   : poloidal angle coordinates for resonant field
c     toro   : toroidal angle coordinates for resonant field
c     __________________________________________________________________
c     labl   : label for multiple runs     
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
         delcur(ising)=j_c(ising)*ifac/mfac(resnum(ising))*delta(ising)
         corcur(ising)=-j_c(ising)*
     $        (rcormn(resnum(ising))-lcormn(resnum(ising)))
         singcur(ising)=delcur(ising)-corcur(ising)

         fkaxmn=0
         fkaxmn(resnum(ising))=-singcur(ising)*chi1/(twopi*ifac*nn)
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
         singflx_mn(:,ising)=MATMUL(fsurf_indmats,fkaxmn)
         singbno_mn(:,ising)=singflx_mn(:,ising)
         CALL ipeq_weight(respsi,singbno_mn(:,ising),mfac,mpert,0)
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
c     compute chirikov parameter.
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
            CALL ipeq_hatoco(respsi,singbno_mn(:,ising),mfac,mpert,
     $           polo,toro)
         ENDIF            
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipec_singfld_l"//slabl//
     $        "_n"//sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipec_singfld_l"//slabl2//
     $        "_n"//sn//".out","UNKNOWN")
      ENDIF
      WRITE(out_unit,*)"IPEC_SINGFLD: "//
     $     "deltas, singular currents and normal fields"
      WRITE(out_unit,'(1x,a12,1x,e16.8)')"distance:",dist
      WRITE(out_unit,'(1x,a12,1x,I16)')"msing:",msing
      WRITE(out_unit,'(1x,a6,11(1x,a16))')"q","psi",
     $     "abs(delta)","abs(delcur)","abs(corcur)","abs(singcur)",
     $     "abs(sgpower)","abs(singbno)","singbnost","singmaxtor",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(1x,f6.3,11(1x,e16.8))')
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
c     subprogram 5. ipout_pmodb.
c     compute perturbed mod b.
c     __________________________________________________________________
c     egnum  : eigenmode number without edge_flag
c     xwpimn : edge deformation input when edge_flag
c     polo   : poloidal angle coordinates for decomposition
c     toro   : toroidal angle coordinates for decomposition
c     __________________________________________________________________
c     labl   : label for multiple runs   
c-----------------------------------------------------------------------
      SUBROUTINE ipout_pmodb(egnum,xwpimn,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert,itheta
      CHARACTER(1) :: spolo,storo,slabl
      CHARACTER(2) :: slabl2

      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,xwt_fun,
     $     bvt_fun,bvz_fun

      REAL(r8), DIMENSION(:), POINTER :: psis
      COMPLEX(r8), DIMENSION(:,:), POINTER :: eulbpar_mn,lagbpar_mn,
     $     eulbparfun,lagbparfun
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"computing perturbed b field"
      ALLOCATE(psis(rstep),
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
c-----------------------------------------------------------------------
c     compute mod b variations.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psis(istep),theta(itheta),1)
            eulbparfun(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(ffun%f(1)*eqfun%f(1))
            lagbparfun(istep,itheta)=
     $           eulbparfun(istep,itheta)+
     $           xwp_fun(itheta)*eqfun%fx(1)+xwt_fun(itheta)*eqfun%fy(1)
         ENDDO
c-----------------------------------------------------------------------
c     decompose components on the given coordinates.
c-----------------------------------------------------------------------
         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
     $        mthsurf,eulbpar_mn(istep,:))
         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
     $        mthsurf,lagbpar_mn(istep,:))
         IF ((polo /= 1).OR.(toro /= 1)) THEN 
            CALL ipeq_hatoco(psis(istep),eulbpar_mn(istep,:),mfac,mpert,
     $           polo,toro)
            CALL ipeq_hatoco(psis(istep),lagbpar_mn(istep,:),mfac,mpert,
     $           polo,toro)
         ENDIF
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipec_pmodbmn_p"//spolo//"_t"
     $        //storo//"_l"//slabl//"_n"//sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipec_pmodbmn_p"//spolo//"_t"
     $        //storo//"_l"//slabl2//"_n"//sn//".out","UNKNOWN")
      ENDIF      

      WRITE(out_unit,*)"IPEC_PMODBMN: "//
     $     "components in perturbed mod b on flux surfaces"
      WRITE(out_unit,'(1x,a8,1x,I6)')"rstep:",rstep
      WRITE(out_unit,'(1x,a8,1x,I6)')"mpert:",mpert
      WRITE(out_unit,'(6(1x,a12))')"psi","mfac",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)"
      DO istep=1,rstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,e16.8,1x,I16,4(2x,e16.8))')
     $           psis(istep),mfac(ipert),
     $           REAL(eulbpar_mn(istep,ipert)),
     $           AIMAG(eulbpar_mn(istep,ipert)),
     $           REAL(lagbpar_mn(istep,ipert)),
     $           AIMAG(lagbpar_mn(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(psis,eulbpar_mn,lagbpar_mn,eulbparfun,lagbparfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_pmodb
c-----------------------------------------------------------------------
c     subprogram 6. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c     __________________________________________________________________
c     egnum  : eigenmode number without edge_flag
c     xwpimn : edge deformation when edge_flag
c     nr     : r grid number
c     nz     : z grid number
c     __________________________________________________________________
c     labl   : label for multiple runs   
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
      CHARACTER(2) :: slabl2

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
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipec_eqbrzphi_l"//slabl//"_n"
     $        //sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipec_eqbrzphi_l"//slabl2//"_n"
     $        //sn//".out","UNKNOWN")
      ENDIF    

      WRITE(out_unit,*)"IPEC_EQBRZPHI: eq b field in rzphi grid"
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

      IF (labl < 10) THEN
         CALL ascii_open(out_unit,"ipec_xrzphi_l"//slabl//"_n"
     $        //sn//".out","UNKNOWN")
      ELSE
         CALL ascii_open(out_unit,"ipec_xrzphi_l"//slabl2//"_n"
     $        //sn//".out","UNKNOWN")
      ENDIF       

      WRITE(out_unit,*)"IPEC_XRZPHI: displacements in rzphi grid"
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

      IF (labl < 10) THEN
         CALL ascii_open(out_unit,"ipec_brzphi_l"//slabl//"_n"
     $        //sn//".out","UNKNOWN")
      ELSE
         CALL ascii_open(out_unit,"ipec_brzphi_l"//slabl2//"_n"
     $        //sn//".out","UNKNOWN")
      ENDIF   

      WRITE(out_unit,*)"IPEC_BRZPHI: perturbed field in rzphi grid"
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
