c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     write various output results of ipec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipout_mod
c      1. ipout_response
c      2. ipout_singcoup
c      3. ipout_control
c      4. ipout_singfld
c      5. ipout_pmodb
c      6. ipout_xbnormal
c      7. ipout_xbrzphi
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
      USE ipdiag_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipout_response.
c     write basic information.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_response
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: lwork,i,j,k
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(mpert) :: sts,s,ev
      COMPLEX(r8), DIMENSION(mpert) :: cev
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,b,u,vt

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(3*mpert) :: work
c-----------------------------------------------------------------------
c     svd analysis for permeability matrix.
c-----------------------------------------------------------------------
      a=permeabmats(resp_index,:,:)
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=3*mpert
      CALL zgesvd('S','S',mpert,mpert,a,mpert,s,u,mpert,vt,mpert,
     $     work,lwork,rwork,info)
      DO i=1,mpert 
         sts(i)=-et(i)/(surfei(i)+surfee(i))
         IF (sts(i) > 0) s(i)=-s(i)
         s(i)=-1/s(i)
      ENDDO

      CALL ascii_open(out_unit,"ipec_response_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_RESPONSE: Response parameters"
      WRITE(out_unit,*)
      WRITE(out_unit,'(3(1x,a12,I4))')
     $     "mpert =",mpert,"mlow =",mlow,"mhigh =",mhigh
      WRITE(out_unit,'(2(1x,a12,es16.8))')
     $     "psilim =",psilim,"qlim =",qlim

      WRITE(out_unit,*)
      WRITE(out_unit,*)"Energy for dcon eigenmodes"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,9(1x,a12))')"mode","ev0","ev1","iv1",
     $     "ep0","ep1","ep2","ep3","ep4","et0"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,9(1x,es12.3))')i,ee(i),surfee(i),
     $        surfei(i),ep(i),surfep(1,i),surfep(2,i),surfep(3,i),
     $        surfep(4,i),et(i)
      ENDDO
      WRITE(out_unit,*)

      WRITE(out_unit,*)"Stability indices"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,2(1x,a12))')"mode","s","se"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,2(1x,es12.3))')i,s(i),sts(i)
      ENDDO
      WRITE(out_unit,*)

      b=surf_indevmats
      ev=surf_indev

      WRITE(out_unit,*)"Vacuum inductance eigenvectors"
      WRITE(out_unit,*)
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,'(1x,a12,es12.3)')"eigenvalue =",ev(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*)
      ENDDO

      b=plas_indevmats(resp_index,:,:)
      ev=plas_indev(resp_index,:)

      WRITE(out_unit,*)"Plasma inductance eigenvectors"
      WRITE(out_unit,*)
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,'(1x,a12,es12.3)')"eigenvalue =",ev(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO
         WRITE(out_unit,*) 
      ENDDO

      b=reluctevmats(resp_index,:,:)
      ev=reluctev(resp_index,:)

      WRITE(out_unit,*)"Reluctance eigenvectors"
      WRITE(out_unit,*) 
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,'(1x,a12,es12.3)')"eigenvalue =",ev(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      b=permeabevmats(resp_index,:,:)
      cev=permeabev(resp_index,:)

      WRITE(out_unit,*)"Permeability eigenvectors"
      WRITE(out_unit,*) 
      DO i=1,mpert
         k=permeabindex(resp_index,i)
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,'(2(1x,a12,es12.3))')
     $        "real(eigen) =",REAL(cev(k)),"imag(eigen) =",AIMAG(cev(k))
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,k)),AIMAG(b(j,k))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      b=surf_indmats

      WRITE(out_unit,*)"Vacuum inductance matrix"
      WRITE(out_unit,*) 
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      b=plas_indmats(resp_index,:,:)

      WRITE(out_unit,*)"Plasma inductance matrix"
      WRITE(out_unit,*) 
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      b=reluctmats(resp_index,:,:)

      WRITE(out_unit,*)"Reluctance matrix"
      WRITE(out_unit,*) 
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      b=a
      WRITE(out_unit,*)"Permeability matrix"
      WRITE(out_unit,*) 
      DO i=1,mpert
         WRITE(out_unit,'(1x,a6,I4)')"mode =",i
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a12))')"m","real","imag"
         DO j=1,mpert
            WRITE(out_unit,'(1x,I4,2(1x,es12.3))')
     $           mfac(j),REAL(b(j,i)),AIMAG(b(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_response
c-----------------------------------------------------------------------
c     subprogram 2. ipout_singcoup.
c     compute coupling between singular surfaces and external fields.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singcoup(spot,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      REAL(r8), INTENT(IN) :: spot

      REAL(r8), DIMENSION(msing) :: s,s1,s2,s3
      COMPLEX(r8), DIMENSION(msing,msing) :: u

      INTEGER :: i,j,itheta,ising,resnum,rsing,rpert,
     $     tmlow,tmhigh,tmpert,lwork,info
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,shear,jarea,thetai
      COMPLEX(r8) :: lbwp1mn,rbwp1mn,lcorrec,rcorrec

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun,wcfun,
     $     dphi,thetas,units,jacs,jacfac

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,lcormn,rcormn,
     $     fkaxmn,singflx_mn,ftnmn
      COMPLEX(r8), DIMENSION(lmpert) :: lftnmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun,
     $     ftnfun

      REAL(r8), DIMENSION(msing) :: area,j_c,w_c

      COMPLEX(r8), DIMENSION(msing,mpert) :: deltas,delcurs,corcurs,
     $     singcurs,singbnoflxs,islandhwids
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: convmat
      COMPLEX(r8), DIMENSION(msing,mpert,mpert) :: fsurfindmats
      COMPLEX(r8), DIMENSION(:), POINTER :: fldflxmn
      COMPLEX(r8), DIMENSION(:,:), POINTER ::  temp1,temp2,
     $     t1mat,t2mat,t3mat

      INTEGER, DIMENSION(:), POINTER :: ipiv,tmfac
      REAL(r8), DIMENSION(:), POINTER :: rwork
      COMPLEX(r8), DIMENSION(:), POINTER :: work
      COMPLEX(r8), DIMENSION(:,:), POINTER :: a,vt    

      TYPE(spline_type) :: spl 
c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      WRITE(*,*)
     $     "Computing coupling between total resonant fields and "//
     $     "external fields"
      
      WRITE(*,*)"Computing surface inductances at resonant surfaces"
      DO ising=1,msing
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,0)
         area(ising)=0
         j_c(ising)=0
         w_c(ising)=0
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
            wcfun(itheta)=SQRT(jac*sqreqb(itheta))
            area(ising)=area(ising)+jac*delpsi(itheta)/mthsurf
            j_c(ising)=j_c(ising)+jac*delpsi(itheta)
     $           *jcfun(itheta)/mthsurf
            w_c(ising)=w_c(ising)+0.5*wcfun(itheta)/mthsurf
         ENDDO
         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
         fsurfindmats(ising,:,:)=fsurf_indmats
         DEALLOCATE(fsurf_indev,fsurf_indmats)
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
         binmn=0
         binmn(i)=1.0
         CALL ipeq_weight(psilim,binmn,mfac,mpert,1)
         IF (fixed_boundary_flag) THEN
            boutmn=binmn
         ELSE
            boutmn=MATMUL(permeabmats(resp_index,:,:),binmn)
         ENDIF
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
            lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
            CALL ipeq_sol(lpsi)
            lbwp1mn=bwp1_mn(resnum)
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
               correc=(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
               lcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
            
            rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1)) 
            CALL ipeq_sol(rpsi)
            rbwp1mn=bwp1_mn(resnum)
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
               correc=(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
               rcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
            deltas(ising,i)=rbwp1mn-lbwp1mn
            delcurs(ising,i)=j_c(ising)*ifac/(twopi*mfac(resnum))*
     $           deltas(ising,i)
            corcurs(ising,i)=-j_c(ising)*
     $           (rcormn(resnum)-lcormn(resnum))
            delcurs(ising,i)=-delcurs(ising,i)/(twopi*ifac*nn)
            corcurs(ising,i)=-corcurs(ising,i)/(twopi*ifac*nn)
            singcurs(ising,i)=delcurs(ising,i)-corcurs(ising,i)

            fkaxmn=0
            fkaxmn(resnum)=singcurs(ising,i)

            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
            singbnoflxs(ising,i)=singflx_mn(resnum)/area(ising)
            islandhwids(ising,i)=4*singflx_mn(resnum)/
     $           (twopi*shear*sq%f(4)*chi1)
         ENDDO
         WRITE(*,'(1x,a16,i4,a22,es10.3,a10,es10.3)')
     $        "poloidal mode = ",mfac(i),", resonant coupling = ",
     $        SUM(ABS(singbnoflxs(:,i)))/msing,", error = ",
     $        SUM(ABS(corcurs(:,i)))/SUM(ABS(singcurs(:,i)))
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     convert coordinates for matrix on the plasma boundary.
c-----------------------------------------------------------------------
      IF ((jac_out /= jac_type).OR.(tout==0)) THEN
         WRITE(*,*)"Converting coordinates"
         CALL spline_eval(sq,psilim,0)
         dphi=0
         
         CALL spline_alloc(spl,mthsurf,2)
         spl%xs=theta
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psilim,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            bpfac=psio*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            spl%fs(itheta,1)=fac/(r(itheta)**rout*rfac**rcout)*
     $           bpfac**bpout*bfac**bout
            ! jacobian for coordinate angle at dcon angle
            spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**rout*rfac**rcout/
     $           (bpfac**bpout*bfac**bout)  
            IF (tout .EQ. 0) THEN
               dphi(itheta)=rzphi%f(3)
            ENDIF
         ENDDO      
         
         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)

         ! coordinate angle at dcon angle
         thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
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

         CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     convert coordinates. 
c-----------------------------------------------------------------------
         ALLOCATE(fldflxmn(lmpert),fldflxmat(lmpert,lmpert),
     $        t1mat(msing,lmpert),t2mat(msing,lmpert),
     $        t3mat(msing,lmpert),tmfac(lmpert))
         DO i=1,lmpert
            lftnmn=0
            lftnmn(i)=1.0
            ! compute given function in dcon angle
            DO itheta=0,mthsurf
               ftnfun(itheta)=0
               ftnfun(itheta)=
     $              lftnmn(i)*EXP(ifac*twopi*lmfac(i)*thetas(itheta))  
            ENDDO
            ! multiply toroidal factor for dcon angle
            IF (tout .EQ. 0) THEN
               ftnfun(:)=ftnfun(:)*EXP(-ifac*nn*dphi(:))
            ELSE
               ftnfun(:)=ftnfun(:)*
     $              EXP(-twopi*ifac*nn*sq%f(4)*(thetas(:)-theta(:)))
            ENDIF
            CALL iscdftf(mfac,mpert,ftnfun,mthsurf,ftnmn)
            convmat(:,i)=ftnmn
         
            CALL iscdftb(lmfac,lmpert,ftnfun,mthsurf,lftnmn)
            ftnfun(:)=ftnfun(:)*sqrt(jacfac(:))
            CALL iscdftf(lmfac,lmpert,ftnfun,mthsurf,fldflxmn)            
            fldflxmat(:,i)=fldflxmn/sqrt(jarea)
         ENDDO
         tmlow = lmlow
         tmhigh = lmhigh
         tmpert = lmpert
         tmfac= lmfac
         t1mat = MATMUL(singbnoflxs,convmat)
         t2mat = MATMUL(singcurs,convmat)*twopi*nn
         t3mat = MATMUL(islandhwids,convmat)
      ELSE
         ALLOCATE(fldflxmn(mpert),fldflxmat(mpert,mpert),
     $        t1mat(msing,mpert),t2mat(msing,mpert),
     $        t3mat(msing,mpert),tmfac(mpert))
         units = (/(1.0,itheta=0,mthsurf)/)
         jarea = issurfint(units,mthsurf,psilim,0,0)
         DO i=1,mpert
            ftnmn=0
            ftnmn(i)=1.0
            fldflxmn=ftnmn
            CALL ipeq_weight(psilim,fldflxmn,mfac,mpert,2)
            fldflxmat(:,i)=fldflxmn/sqrt(jarea)            
         ENDDO
         tmlow = mlow
         tmhigh = mhigh
         tmpert = mpert
         tmfac = mfac
         t1mat = singbnoflxs
         t2mat = singcurs*twopi*nn
         t3mat = islandhwids
      ENDIF
c-----------------------------------------------------------------------
c     svd analysis.
c-----------------------------------------------------------------------
      lwork=3*tmpert
      ALLOCATE(ipiv(tmpert),rwork(5*msing),work(lwork),
     $     a(msing,tmpert),vt(msing,tmpert),
     $     t1v(tmpert,msing),t2v(tmpert,msing),t3v(tmpert,msing),
     $     temp1(tmpert,tmpert),temp2(tmpert,tmpert))
      
      temp1=fldflxmat
      temp2=0
      DO i=1,tmpert
         temp2(i,i)=1
      ENDDO
      CALL zgetrf(tmpert,tmpert,temp1,tmpert,ipiv,info)
      CALL zgetrs('N',tmpert,tmpert,temp1,tmpert,ipiv,temp2,tmpert,info)

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t1mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s1=s
      t1v=MATMUL(fldflxmat,CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t2mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s2=s
      t2v=MATMUL(fldflxmat,CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t3mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s3=s
      t3v=MATMUL(fldflxmat,CONJG(TRANSPOSE(vt)))  
c-----------------------------------------------------------------------
c     write matrix.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singcoup_matrix_n"//
     $     sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGCOUP_MATRIX: Coupling matrices"//
     $     " between resonant field and external field"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout  
      WRITE(out_unit,'(4(1x,a12,I4))')
     $     "msing =",msing,"mpert =",tmpert,
     $     "mlow =",tmlow,"mhigh =",tmhigh
      WRITE(out_unit,'(2(1x,a12,es16.8))')
     $     "psilim =",psilim,"qlim =",qlim
      WRITE(out_unit,*)

      WRITE(out_unit,*)"Coupling matrix to resonant fields"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es16.8)')"q =",singtype(i)%q,
     $        "psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t1mat(i,j)),AIMAG(t1mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"Coupling matrix to singular currents"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es16.8)')"q =",singtype(i)%q,
     $        "psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t2mat(i,j)),AIMAG(t2mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"Coupling matrix to island half-widths"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es16.8)')"q =",singtype(i)%q,
     $        "psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t3mat(i,j)),AIMAG(t3mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO
      CALL ascii_close(out_unit)
      
      CALL ascii_open(out_unit,"ipec_singcoup_svd_n"//
     $     sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGCOUP_SVD: SVD analysis"//
     $     " for coupling matrices"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a12,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tmag_out  
      WRITE(out_unit,'(4(1x,a12,I4))')
     $     "msing =",msing,"mpert =",tmpert,
     $     "mlow =",tmlow,"mhigh =",tmhigh
      WRITE(out_unit,'(2(1x,a12,es16.8))')
     $     "psilim =",psilim,"qlim =",qlim
      WRITE(out_unit,*)

      WRITE(out_unit,*)"SVD for coupling matrix to resonant fields"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es16.8)')
     $        "mode =",i,"s =",s1(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t1v(j,i)),AIMAG(t1v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"SVD for coupling matrix to singular currents"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es16.8)')
     $        "mode =",i,"s =",s2(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t2v(j,i)),AIMAG(t2v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"SVD for coupling matrix to island half-widths"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es16.8)')
     $        "mode =",i,"s =",s3(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(1x,es16.8))')
     $           tmfac(j),REAL(t3v(j,i)),AIMAG(t3v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO
      CALL ascii_close(out_unit)

      DEALLOCATE(tmfac,t1mat,t2mat,t3mat,fldflxmn)
      DEALLOCATE(ipiv,rwork,work,a,vt,temp1,temp2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singcoup
c-----------------------------------------------------------------------
c     subprogram 3. ipout_control
c     calculate response from external field on the control surface.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_control(ifile,binmn,boutmn,xspomn,rin,bpin,bin,
     $     rcin,tin,jin,rout,bpout,bout,rcout,tout,svd_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: ifile
      INTEGER, INTENT(IN) :: rin,bpin,bin,rcin,tin,jin,
     $     rout,bpout,bout,rcout,tout
      LOGICAL, INTENT(IN) :: svd_flag
      
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: binmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: boutmn,xspomn

      INTEGER :: i,j,i1,i2,i3,ms,itheta,jout
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: vengy,sengy,pengy,area,thetai,scale
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacfac,
     $     sbinfun,sboutfun

      COMPLEX(r8), DIMENSION(mpert) :: finmn,foutmn,
     $     bninmn,bnoutmn,xinmn,xoutmn,xspimn
      COMPLEX(r8), DIMENSION(lmpert) :: cinmn,coutmn,hawmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,xinfun,xoutfun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2

      REAL(r8), DIMENSION(:,:), POINTER :: dcosmn,dsinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn

      TYPE(spline_type) :: spl       
c-----------------------------------------------------------------------
c     check data_type and read data.
c-----------------------------------------------------------------------
      scale=1.0
      IF (data_flag) THEN
         WRITE(*,*)"Reading external fields on the control surface"
         ALLOCATE(dcosmn(mmin:mmax,nmin:nmax),
     $        dsinmn(mmin:mmax,nmin:nmax),rawmn(mmin:mmax,nmin:nmax))
         IF ((data_type == "surfmn1").OR.(data_type == "surfmn2")) THEN
            i1 = mmax
            i2 = mmin
            i3 = -1
            scale = 1e-4
         ELSE
            i1 = mmin
            i2 = mmax
            i3 = 1
         ENDIF
         CALL ascii_open(in_unit,ifile,"old")
                        
         DO i=i1,i2,i3
            IF (data_type == "surfmn1") THEN
               READ(in_unit,1000) (dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,1000) (dsinmn(i,j),j=nmin,nmax)
            ELSE IF (data_type == "surfmn2") THEN
               READ(in_unit,1001) (dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,1001) (dsinmn(i,j),j=nmin,nmax)
            ELSE IF (data_type == "vac3d1") THEN
               READ(in_unit,1010) (dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,1010) (dsinmn(i,j),j=nmin,nmax)
            ELSE IF (data_type == "vac3d2") THEN
               READ(in_unit,1020) ms,dcosmn(i,nn),dsinmn(i,nn)               
            ELSE
               WRITE(message,'(a)')"Can't recognize data format"
               CALL ipec_stop(message)
            ENDIF            
         ENDDO

         CALL ascii_close(in_unit)
         rawmn=dcosmn+ifac*dsinmn
         hawmn=rawmn(:,nn)

         DEALLOCATE(dcosmn,dsinmn,rawmn)
         
      ELSE IF (harmonic_flag) THEN

         DO i=-hmnum,hmnum
            IF ((-mmin+i>=1).AND.(-mmin+i<=lmpert)) THEN
               hawmn(-mmin+i+1)=cosmn(i)+ifac*sinmn(i)
            ENDIF
         ENDDO

      ENDIF

 1000 FORMAT(1x,25f12.6)         
 1001 FORMAT(1x,33f12.6)
 1010 FORMAT(11(1x,e15.8))
 1020 FORMAT(1x,I4,2(1x,e15.8))
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      IF ((data_flag) .OR. (harmonic_flag)) THEN
         CALL ipeq_fcoords(psilim,hawmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)             
         binmn=0
         DO i=1,lmpert
            IF ((lmlow-mlow+i>=1).AND.(lmlow-mlow+i<=mpert)) THEN
               binmn(lmlow-mlow+i)=hawmn(i)
            ENDIF
         ENDDO                 
      ENDIF         
c-----------------------------------------------------------------------
c     get the plasma response at the control surface.
c-----------------------------------------------------------------------
      binmn=binmn*scale
      finmn=binmn
      CALL ipeq_weight(psilim,finmn,mfac,mpert,1)
      IF (fixed_boundary_flag) THEN
         boutmn=finmn
      ELSE
         boutmn=MATMUL(permeabmats(resp_index,:,:),finmn)
      ENDIF
      foutmn=boutmn
      CALL ipeq_weight(psilim,boutmn,mfac,mpert,0)
      bninmn=binmn
      bnoutmn=boutmn
      CALL ipeq_bntoxp(psilim,bninmn,xspimn)
      CALL ipeq_bntoxp(psilim,bnoutmn,xspomn)
      xinmn=xspimn
      xoutmn=xspomn
      CALL ipeq_weight(psilim,xinmn,mfac,mpert,4)
      CALL ipeq_weight(psilim,xoutmn,mfac,mpert,4)
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

      vy = 0.5*SUM(CONJG(finmn)*MATMUL(temp1,finmn))
      sy = 0.5*SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))
      vengy = REAL(vy)
      sengy = REAL(sy)

      temp1=0
      work2=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2 = plas_indmats(resp_index,:,:)
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
      py = 0.5*SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))
      pengy = REAL(py)
      WRITE(*,'(1x,a,es10.3)')"required energy to perturb vacuum = ",
     $     sengy
      WRITE(*,'(1x,a,es10.3)')"required energy to perturb plasma = ",
     $     pengy
      WRITE(*,'(1x,a,es10.3)')"amplification factor = ",
     $     sengy/pengy
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_control_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_CONTROL: "//
     $     "Plasma response for an external perturbation on the "//
     $     "control surface"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a12,I4)')"mpert =",mpert
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"vacuum energy =",vengy
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"surface energy =",sengy
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"plasma energy =",pengy
      WRITE(out_unit,*)

      WRITE(out_unit,*)"jac_type = "//jac_type
      WRITE(out_unit,*)

      WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)","imag(bin)",
     $     "real(bout)","imag(bout)"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,4(1x,es16.8))')mfac(i),
     $        REAL(binmn(i)),AIMAG(binmn(i)),
     $        REAL(boutmn(i)),AIMAG(boutmn(i))
      ENDDO
      WRITE(out_unit,*)

      IF ((jac_in /= jac_type).OR.(tin==0).OR.(jin==1)) THEN
         cinmn=0
         coutmn=0
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               cinmn(mlow-lmlow+i)=binmn(i)
               coutmn(mlow-lmlow+i)=boutmn(i)
            ENDIF
         ENDDO         
         CALL ipeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin) 
         CALL ipeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)
 
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_in = ",jac_in,
     $        "jsurf_in =",jin,"tmag_in =",tin
         WRITE(out_unit,*)             
         WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)","imag(bin)",
     $        "real(bout)","imag(bout)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,4(1x,es16.8))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i))
         ENDDO 
         WRITE(out_unit,*)       
      ENDIF

      IF ((jac_out /= jac_type).OR.(tout==0)) THEN
         cinmn=0
         coutmn=0
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               cinmn(mlow-lmlow+i)=binmn(i)
               coutmn(mlow-lmlow+i)=boutmn(i)
            ENDIF
         ENDDO  
         jout=0
         CALL ipeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout) 
         CALL ipeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout)
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_out = ",jac_out,
     $        "jsurf_out =",jout,"tmag_out =",tout  
         WRITE(out_unit,*)          
         WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)","imag(bin)",
     $        "real(bout)","imag(bout)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,4(1x,es16.8))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i))
         ENDDO 
         WRITE(out_unit,*)
         IF (svd_flag) THEN
            ALLOCATE(sbno_mn(lmpert),sbno_fun(0:mthsurf))
            sbno_mn=MATMUL(fldflxmat,cinmn)
            CALL iscdftb(lmfac,lmpert,sbno_fun,mthsurf,sbno_mn)
            sbno_mn=cinmn
         ENDIF
      ELSE
         IF (svd_flag) THEN
            ALLOCATE(sbno_mn(mpert),sbno_fun(0:mthsurf))
            sbno_mn=MATMUL(fldflxmat,binmn)
            CALL iscdftb(mfac,mpert,sbno_fun,mthsurf,sbno_mn)
            sbno_mn=binmn
         ENDIF
      ENDIF
      CALL ascii_close(out_unit)

      IF (fun_flag) THEN
         CALL ipeq_bcoords(psilim,xinmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,xoutmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,bninmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,bnoutmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)

         CALL iscdftb(mfac,mpert,xinfun,mthsurf,xinmn)    
         CALL iscdftb(mfac,mpert,xoutfun,mthsurf,xoutmn)     
         CALL iscdftb(mfac,mpert,binfun,mthsurf,bninmn)     
         CALL iscdftb(mfac,mpert,boutfun,mthsurf,bnoutmn)    
         
         CALL ascii_open(out_unit,"ipec_control_fun_n"//sn//
     $     ".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_CONTROL_FUN: "//
     $        "Plasma response for an external perturbation on the "//
     $        "control surface in functions"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a12,I4)')"mthsurf =",mthsurf
         WRITE(out_unit,'(1x,a16,1x,es16.8)')"vacuum energy =",vengy
         WRITE(out_unit,'(1x,a16,1x,es16.8)')"surface energy =",sengy
         WRITE(out_unit,'(1x,a16,1x,es16.8)')"plasma energy =",pengy
         WRITE(out_unit,*)
         
         WRITE(out_unit,*)"jac_type = "//jac_type
         WRITE(out_unit,*)
         
         WRITE(out_unit,'(10(1x,a16))')"r","z",
     $        "real(xin)","imag(xin)","real(xout)","imag(xout)",
     $        "real(bin)","imag(bin)","real(bout)","imag(bout)"
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psilim,theta(itheta),0)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)         
            WRITE(out_unit,'(10(1x,es16.8))')r(itheta),z(itheta),
     $           REAL(xinfun(itheta)),AIMAG(xinfun(itheta)),
     $           REAL(xoutfun(itheta)),AIMAG(xoutfun(itheta)),
     $           REAL(binfun(itheta)),AIMAG(binfun(itheta)),
     $           REAL(boutfun(itheta)),AIMAG(boutfun(itheta))
         ENDDO
         CALL ascii_close(out_unit)  
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_control
c-----------------------------------------------------------------------
c     subprogram 4. ipout_singfld.
c     compute current and field on rational surfaces.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singfld(egnum,xspimn,spot,
     $     rout,bpout,bout,rcout,tout,svd_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      REAL(r8), INTENT(IN) :: spot
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspimn
      LOGICAL, INTENT(IN) :: svd_flag

      INTEGER :: i,itheta,ising
      REAL(r8) :: respsi,lpsi,rpsi,sqrpsi,correc,shear,hdist,sbnosurf
      COMPLEX(r8) :: lbwp1mn,rbwp1mn,lcorrec,rcorrec

      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: area,j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: lcormn,rcormn,fkaxmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun

      REAL(r8), DIMENSION(msing) :: island_hwidth,chirikov,
     $     novf,novs,novi
      COMPLEX(r8), DIMENSION(msing) :: delta,delcur,corcur,singcur,
     $     ovf,ovs,ovi
      COMPLEX(r8), DIMENSION(mpert,msing) :: singflx_mn
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing total resonant fields"
      CALL ipeq_alloc
      CALL idcon_build(egnum,xspimn)
c-----------------------------------------------------------------------
c     evaluate delta and singular currents.
c     delta is delta*chi1*sq%f(4) and j_c is j_c/(chi1*sq%f(4))
c-----------------------------------------------------------------------
      DO ising=1,msing
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)
         area(ising)=0
         j_c(ising)=0
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
            area(ising)=area(ising)+jac*delpsi(itheta)/mthsurf
            j_c(ising)=j_c(ising)+jac*delpsi(itheta)
     $           *jcfun(itheta)/mthsurf
         ENDDO
         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0
         shear=mfac(resnum(ising))*sq%f1(4)/sq%f(4)**2

         lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(lpsi)
         lbwp1mn=bwp1_mn(resnum(ising))
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
            correc=(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $           (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
            lcorfun(itheta)=bwp_fun(itheta)*correc
         ENDDO
         CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
         
         rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(rpsi)
         rbwp1mn=bwp1_mn(resnum(ising))
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
            correc=(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $           (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/sqrpsi
            rcorfun(itheta)=bwp_fun(itheta)*correc
         ENDDO
         CALL iscdftf(mfac,mpert,rcorfun,mthsurf,rcormn)
         delta(ising)=rbwp1mn-lbwp1mn
         delcur(ising)=j_c(ising)*ifac/(twopi*mfac(resnum(ising)))*
     $        delta(ising)
         corcur(ising)=-j_c(ising)*
     $        (rcormn(resnum(ising))-lcormn(resnum(ising)))
         delcur(ising)=-delcur(ising)/(twopi*ifac*nn)
         corcur(ising)=-corcur(ising)/(twopi*ifac*nn)
         singcur(ising)=delcur(ising)-corcur(ising)

         fkaxmn=0
         fkaxmn(resnum(ising))=singcur(ising)
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
         singflx_mn(:,ising)=MATMUL(fsurf_indmats,fkaxmn)
         DEALLOCATE(fsurf_indmats,fsurf_indev)
c-----------------------------------------------------------------------
c     compute half-width of magnetic island.
c-----------------------------------------------------------------------
         island_hwidth(ising)=
     $        SQRT(ABS(4*singflx_mn(resnum(ising),ising)/
     $        (twopi*shear*sq%f(4)*chi1)))
c-----------------------------------------------------------------------
c     compute coordinate-independent resonant field.
c----------------------------------------------------------------------- 
         singflx_mn(:,ising)=singflx_mn(:,ising)/area(ising)
c-----------------------------------------------------------------------
c     compute pseudo-chirikov parameter.
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
         WRITE(*,'(1x,a6,es10.3,a6,f6.3,a26,es10.3)')
     $        "psi = ",singtype(ising)%psifac,
     $        ", q = ",singtype(ising)%q,
     $        ", total resonant fields = ",
     $        ABS(singflx_mn(resnum(ising),ising))
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singfld_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGFLD: "//
     $     "Resonant fields, singular currents, and islands"
      WRITE(out_unit,*)       
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout 
      WRITE(out_unit,'(1x,a12,1x,es16.8)')"sweet-spot =",spot
      WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a6,7(1x,a16))')"q","psi",
     $     "real(singflx)","imag(singflx)",
     $     "real(singcur)","imag(singcur)",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(1x,f6.3,7(1x,es16.8))')
     $        singtype(ising)%q,singtype(ising)%psifac,
     $        REAL(singflx_mn(resnum(ising),ising)),
     $        AIMAG(singflx_mn(resnum(ising),ising)),
     $        REAL(singcur(ising))*twopi*nn,
     $        AIMAG(singcur(ising))*twopi*nn,
     $        island_hwidth(ising),chirikov(ising)
      ENDDO
      WRITE(out_unit,*)
 
      IF (svd_flag) THEN
         sbnosurf=SQRT(ABS(DOT_PRODUCT(sbno_fun(1:mthsurf),
     $        sbno_fun(1:mthsurf)))/mthsurf)
         DO ising=1,msing
            ovf(ising)=DOT_PRODUCT(t1v(:,ising),sbno_mn(:))
            ovs(ising)=DOT_PRODUCT(t2v(:,ising),sbno_mn(:))
            ovi(ising)=DOT_PRODUCT(t3v(:,ising),sbno_mn(:))
         ENDDO
         DO ising=1,msing
            novf(ising)=ABS(ovf(ising))/sbnosurf*1e2
            novs(ising)=ABS(ovs(ising))/sbnosurf*1e2
            novi(ising)=ABS(ovi(ising))/sbnosurf*1e2
         ENDDO

         WRITE(out_unit,*)"Overlap fields, overlap singular "//
     $        "currents, and overlap islands"   
         WRITE(out_unit,*)        
         WRITE(out_unit,'(1x,a6,9(1x,a16))')
     $        "mode","real(ovf)","imag(ovf)","overlap(%)",
     $        "real(ovs)","imag(ovs)","overlap(%)",
     $        "real(ovi)","imag(ovi)","overlap(%)"
         DO ising=1,msing
            WRITE(out_unit,'(1x,I6,9(1x,es16.8))')ising,
     $           REAL(ovf(ising)),AIMAG(ovf(ising)),novf(ising),
     $           REAL(ovs(ising)),AIMAG(ovs(ising)),novs(ising),
     $           REAL(ovi(ising)),AIMAG(ovi(ising)),novi(ising)
         ENDDO
         WRITE(out_unit,*)
      ENDIF
        
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_singfld
c-----------------------------------------------------------------------
c     subprogram 5. ipout_pmodb.
c     compute perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_pmodb(egnum,xspimn,rout,bpout,bout,rcout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspimn

      INTEGER :: istep,ipert,itheta,iindex
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(rstep) :: psis
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      REAL(r8), DIMENSION(rstep,0:mthsurf) :: rs,zs
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun
      COMPLEX(r8), DIMENSION(rstep,mpert) :: eulbparmns,lagbparmns
      COMPLEX(r8), DIMENSION(rstep,0:mthsurf) :: eulbparfuns,lagbparfuns
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing total |b| fields"

      IF (rstep .EQ. mstep) THEN
         psis=psifac
      ELSE
         psis=(/(istep,istep=1,rstep)/)/REAL(rstep,r8)*(psilim-psilow)
      ENDIF
      CALL idcon_build(egnum,xspimn)

      CALL ipeq_alloc
      DO istep=1,rstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(rstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(rstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a18)')
     $        "volume = ",iindex,"% |b| computations"
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psis(istep),1)
         CALL ipeq_sol(psis(istep))
         CALL ipeq_contra(psis(istep))         
         CALL ipeq_cova(psis(istep))
c-----------------------------------------------------------------------
c     compute mod b variations in hamada.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psis(istep),theta(itheta),1)
            CALL bicube_eval(rzphi,psis(istep),theta(itheta),0)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            eulbparfuns(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(rzphi%f(4)*eqfun%f(1))
            lagbparfuns(istep,itheta)=
     $           eulbparfuns(istep,itheta)+
     $           xsp_fun(itheta)*eqfun%fx(1)+
     $           xms_fun(itheta)/(chi1*sq%f(4))*eqfun%fy(1)
            dphi(itheta)=rzphi%f(3)
         ENDDO
         CALL iscdftf(mfac,mpert,eulbparfuns(istep,:),
     $        mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbparfuns(istep,:),
     $        mthsurf,lagbpar_mn)
c-----------------------------------------------------------------------
c     decompose components on the given coordinates.
c-----------------------------------------------------------------------
         IF ((jac_out /= jac_type)) THEN
            CALL ipeq_bcoords(psis(istep),eulbpar_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,1,0)
            CALL ipeq_bcoords(psis(istep),lagbpar_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,1,0)
         ENDIF
         eulbparmns(istep,:)=eulbpar_mn
         lagbparmns(istep,:)=lagbpar_mn
         eulbparfuns(istep,:)=eulbparfuns(istep,:)*EXP(ifac*nn*dphi)
         lagbparfuns(istep,:)=lagbparfuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_pmodb_n"//sn//".out","UNKNOWN")
  
      WRITE(out_unit,*)"IPEC_PMODB: "//
     $     "Components in perturbed mod b"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "rstep =",rstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     

      WRITE(out_unit,'(1x,a16,1x,a4,4(1x,a16))')"psi","m",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)"
      DO istep=1,rstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es16.8,1x,I4,4(1x,es16.8))')
     $           psis(istep),mfac(ipert),
     $           REAL(eulbparmns(istep,ipert)),
     $           AIMAG(eulbparmns(istep,ipert)),
     $           REAL(lagbparmns(istep,ipert)),
     $           AIMAG(lagbparmns(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_pmodb_fun_n"//sn//
     $        ".out","UNKNOWN")
  
         WRITE(out_unit,*)"IPEC_PMODB_FUN: "//
     $        "Components in perturbed mod b in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "rstep =",rstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(6(1x,a16))')"r","z",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)"
         DO istep=1,rstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(6(1x,es16.8))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(eulbparfuns(istep,itheta)),
     $              AIMAG(eulbparfuns(istep,itheta)),
     $              REAL(lagbparfuns(istep,itheta)),
     $              AIMAG(lagbparfuns(istep,itheta))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
      
      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "pmodb.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO istep=1,rstep
               WRITE(bin_unit)REAL(psis(istep),4),
     $              REAL(REAL(eulbparmns(istep,ipert)),4),
     $              REAL(AIMAG(eulbparmns(istep,ipert)),4),
     $              REAL(REAL(lagbparmns(istep,ipert)),4),
     $              REAL(AIMAG(lagbparmns(istep,ipert)),4)     
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"pmodb_2d.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)rstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:rstep,0:mthsurf),4),
     $        REAL(zs(9:rstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(eulbparfuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(eulbparfuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(lagbparfuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(lagbparfuns(9:rstep,0:mthsurf)),4)

         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_pmodb
c-----------------------------------------------------------------------
c     subprogram 6. ipout_xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbnormal(egnum,xspimn,rout,bpout,bout,rcout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspimn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(rstep) :: psis,qs
      REAL(r8), DIMENSION(rstep,0:mthsurf) :: rs,zs,psiss,rvecs,zvecs
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs,dphi
      COMPLEX(r8), DIMENSION(rstep,mpert) :: xmns,ymns,
     $     xnomns,bnomns,bwpmns
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun
      COMPLEX(r8), DIMENSION(rstep,0:mthsurf) :: rss,zss,
     $     xnofuns,bnofuns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing x and b normal components"

      IF (rstep .EQ. mstep) THEN
         psis=psifac
         qs=qfac
      ELSE
         psis=(/(istep,istep=1,rstep)/)/REAL(rstep,r8)*(psilim-psilow)
      ENDIF
      CALL idcon_build(egnum,xspimn)

      CALL ipeq_alloc
      DO istep=1,rstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(rstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(rstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL ipeq_sol(psis(istep))
         CALL ipeq_contra(psis(istep))

         IF (rstep .NE. mstep) THEN
            CALL spline_eval(sq,psis(istep),0)
            qs(istep)=sq%f(4)
         ENDIF

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psis(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            dphi(itheta)=rzphi%f(3)
            rvecs(istep,itheta)=
     $           (cos(eta)*w(1,1)-sin(eta)*w(1,2))/delpsi(itheta)
            zvecs(istep,itheta)=
     $           (sin(eta)*w(1,1)+cos(eta)*w(1,2))/delpsi(itheta)
         ENDDO
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         xnofuns(istep,:)=xwp_fun/(jacs*delpsi)
         bnofuns(istep,:)=bwp_fun/(jacs*delpsi)
         CALL iscdftf(mfac,mpert,xnofuns(istep,:),mthsurf,xno_mn)
         CALL iscdftf(mfac,mpert,bnofuns(istep,:),mthsurf,bno_mn)
         IF ((jac_out /= jac_type)) THEN
            CALL ipeq_bcoords(psis(istep),xno_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,1,0)
            CALL ipeq_bcoords(psis(istep),bno_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,1,0)
            CALL ipeq_bcoords(psis(istep),bwp_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,1,0)            
         ENDIF
         xnomns(istep,:)=xno_mn
         bnomns(istep,:)=bno_mn
         bwpmns(istep,:)=bwp_mn
         xnofuns(istep,:)=xnofuns(istep,:)*EXP(ifac*nn*dphi)
         bnofuns(istep,:)=bnofuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      rss=CMPLX(rs,rs)+xnofuns*rvecs
      zss=CMPLX(zs,zs)+xnofuns*zvecs
      DO itheta=0,mthsurf
         psiss(:,itheta)=psis(:)
      ENDDO

      CALL ascii_open(out_unit,"ipec_xbnormal_n"//
     $     sn//".out","UNKNOWN")

      WRITE(out_unit,*)"IPOUT_XBNROMAL: "//
     $     "Normal components of displacement and field"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "rstep =",rstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $     "real(xno)","imag(xno)","real(bno)","imag(bno)",
     $     "real(bwp)","imag(bwp)"

      DO istep=1,rstep
         DO ipert=1,mpert
            WRITE(out_unit,'(2(1x,es16.8),1x,I4,6(1x,es16.8))')
     $           psis(istep),qs(istep),mfac(ipert),
     $           REAL(xnomns(istep,ipert)),AIMAG(xnomns(istep,ipert)),
     $           REAL(bnomns(istep,ipert)),AIMAG(bnomns(istep,ipert)),
     $           REAL(bwpmns(istep,ipert)),AIMAG(bwpmns(istep,ipert))        
         ENDDO
      ENDDO

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbnormal_fun_n"//
     $        sn//".out","UNKNOWN")
         
         WRITE(out_unit,*)"IPOUT_XBNROMAL_FUN: "//
     $        "Normal components of displacement and field in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "rstep =",rstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(8(1x,a16))')"r","z","rvec","zvec",
     $        "real(xno)","imag(xno)","real(bno)","imag(bno)"
         
         DO istep=1,rstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(8(1x,es16.8))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              rvecs(istep,itheta),zvecs(istep,itheta),
     $              REAL(xnofuns(istep,itheta)),
     $              AIMAG(xnofuns(istep,itheta)),
     $              REAL(bnofuns(istep,itheta)),
     $              AIMAG(bnofuns(istep,itheta))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (flux_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbnormal_flux_n"//
     $        sn//".out","UNKNOWN")
         
         WRITE(out_unit,*)"IPOUT_XBNROMAL_FLUX: "//
     $        "Perturbed flux surfaces"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "rstep =",rstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(3(1x,a16))')"real(r)","imag(r)",
     $        "real(z)","imag(z)","psi"
         
         DO istep=1,rstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(5(1x,es16.8))')
     $              REAL(rss(istep,itheta)),AIMAG(rss(istep,itheta)),
     $              REAL(zss(istep,itheta)),AIMAG(zss(istep,itheta)),
     $              psiss(istep,itheta)
            ENDDO
         ENDDO
         
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "xbnormal.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO istep=1,rstep
               WRITE(bin_unit)REAL(psis(istep),4),
     $              REAL(REAL(xnomns(istep,ipert)),4),
     $              REAL(AIMAG(xnomns(istep,ipert)),4),
     $              REAL(REAL(bnomns(istep,ipert)),4),
     $              REAL(AIMAG(bnomns(istep,ipert)),4),
     $              REAL(REAL(bwpmns(istep,ipert)),4),
     $              REAL(AIMAG(bwpmns(istep,ipert)),4)    
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF
  
      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"xbnormal_2d.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)rstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:rstep,0:mthsurf),4),
     $        REAL(zs(9:rstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(xnofuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(xnofuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(bnofuns(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(bnofuns(9:rstep,0:mthsurf)),4)
         CALL bin_close(bin_2d_unit)

         ximax=MAXVAL(ABS(xnofuns))
         CALL bicube_eval(rzphi,psilim,theta(0),0)
         rmax=SQRT(rzphi%f(1))
         xnofuns=xnofuns/ximax*rmax/6.0
         
         WRITE(*,'(1x,a,es10.3)')"maximum displacement = ",ximax
         WRITE(*,'(1x,a,es10.3)')"scale factor for 2d plots = ",
     $        rmax/(ximax*6.0)
      
         rss=CMPLX(rs,rs)+xnofuns*rvecs
         zss=CMPLX(zs,zs)+xnofuns*zvecs
         DO itheta=0,mthsurf
            psiss(:,itheta)=psis(:)
         ENDDO

         CALL bin_open(bin_2d_unit,
     $        "pflux_re_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)rstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(REAL(rss(9:rstep,0:mthsurf)),4),
     $        REAL(REAL(zss(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psiss(9:rstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         CALL bin_open(bin_2d_unit,
     $        "pflux_im_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)rstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(AIMAG(rss(9:rstep,0:mthsurf)),4),
     $        REAL(AIMAG(zss(9:rstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psiss(9:rstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         DO istep=1,mstep
            xmns(istep,:)=mfac
            ymns(istep,:)=psis(istep)
         ENDDO
         CALL bin_open(bin_2d_unit,"bnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)rstep-9,mpert-1
         WRITE(bin_2d_unit)REAL(xmns(9:rstep,:),4),
     $        REAL(ymns(9:rstep,:),4)
         WRITE(bin_2d_unit)REAL(ABS(bwpmns(9:rstep,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbnormal
c-----------------------------------------------------------------------
c     subprogram 7. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphi(egnum,xspimn,nr,nz,bnimn,bnomn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspimn
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnimn,bnomn

      INTEGER :: i,j,ipert,iindex
      REAL(r8) :: mid,bt0,ileft

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz,ebr,ebz,ebp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: xrr,xrz,xrp,brr,brz,brp,
     $     bpr,bpz,bpp,vbr,vbz,vbp,vpbr,vpbz,vpbp,vvbr,vvbz,vvbp
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      ! initialization
      ebr = 0
      ebz = 0
      ebp = 0
      xrr = 0
      xrz = 0
      xrp = 0
      brr = 0
      brz = 0
      brp = 0
      bpr = 0
      bpz = 0
      bpp = 0
      vbr = 0
      vbz = 0
      vbp = 0
      vpbr = 0
      vpbz = 0
      vpbp = 0
      vvbr = 0
      vvbz = 0
      vvbp = 0

      CALL idcon_build(egnum,xspimn)

      IF (eqbrzphi_flag) THEN
         WRITE(*,*)"Computing equilibrium magnetic fields"
         ! evaluate f value for vacuum
         mid = 0.0
         CALL spline_eval(sq,mid,0)
         bt0 = abs(sq%f(1))/(twopi*ro)
         DO i=0,nr
            DO j=0,nz
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),1)
               ebr(i,j) = -psi_in%fy(1)/gdr(i,j)*psio
               ebz(i,j) = psi_in%fx(1)/gdr(i,j)*psio
               IF (gdl(i,j) == 1) THEN  
                  CALL spline_eval(sq,gdpsi(i,j),0)
                  ebp(i,j) = abs(sq%f(1))/(twopi*gdr(i,j))
               ELSE
                  ebp(i,j) = bt0*ro/gdr(i,j)  
               ENDIF
            ENDDO   
         ENDDO
      ENDIF
      
      CALL ipeq_alloc

      WRITE(*,*)"Mapping fields to cylindrical coordinates"

      IF (brzphi_flag .OR. xrzphi_flag) THEN
         DO i=0,nr
         iindex = FLOOR(REAL(i,8)/FLOOR((nr-1)/10.0))*10
         ileft = REAL(i,8)/FLOOR((nr-1)/10.0)*10-iindex
         IF ((i /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a10)')
     $        "volume = ",iindex,"% mappings"
            DO j=0,nz
               IF (gdl(i,j) == 1) THEN
                  CALL ipeq_sol(gdpsi(i,j))
                  CALL ipeq_contra(gdpsi(i,j))
                  CALL ipeq_cova(gdpsi(i,j))
                  CALL ipeq_rzphi(gdpsi(i,j))
                  DO ipert=1,mpert
                     IF (brzphi_flag) THEN
                        brr(i,j)=brr(i,j)+brr_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                        brz(i,j)=brz(i,j)+brz_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                        brp(i,j)=brp(i,j)+brp_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                     ENDIF
                     IF (xrzphi_flag) THEN
                        xrr(i,j)=xrr(i,j)+xrr_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                        xrz(i,j)=xrz(i,j)+xrz_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                        xrp(i,j)=xrp(i,j)+xrp_mn(ipert)*
     $                       EXP(twopi*ifac*mfac(ipert)*gdthe(i,j))
                     ENDIF
                  ENDDO
                  IF (brzphi_flag) THEN
                     brr(i,j)=
     $                    REAL(brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(brr(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                     brz(i,j)=
     $                    REAL(brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(brz(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                     brp(i,j)=
     $                    REAL(brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(brp(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                  ENDIF
                  IF (xrzphi_flag) THEN
                     xrr(i,j)=
     $                    REAL(xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(xrr(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                     xrz(i,j)=
     $                    REAL(xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(xrz(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                     xrp(i,j)=
     $                    REAL(xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))+
     $                    ifac*REAL(xrp(i,j)*EXP(-twopi*ifac*
     $                    (0.25+nn*gdphi(i,j))))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL ipeq_dealloc

      IF (vbrzphi_flag) THEN
         WRITE(*,*)"Computing vacuum magnetic fields"
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
         IF (brzphi_flag) THEN
            DO i=0,nr
               DO j=0,nz
       
                  IF (gdl(i,j) /= 1) THEN
                     gdl(i,j)=vgdl(i,j)
                     brr(i,j)=vbr(i,j)
                     brz(i,j)=vbz(i,j)
                     brp(i,j)=vbp(i,j)
                  ENDIF

               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (divzero_flag) THEN
         CALL ipeq_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF
      IF (div_flag) THEN
         CALL ipdiag_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF

      IF (vpbrzphi_flag) THEN
         WRITE(*,*)"Computing vacuum fields with plasma response"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (brzphi_flag) THEN
            DO i=0,nr
               DO j=0,nz
                  
                  IF (gdl(i,j) /= 1) THEN
                     gdl(i,j)=vgdl(i,j)
                     bpr(i,j)=vpbr(i,j)
                     bpz(i,j)=vpbz(i,j)
                     bpp(i,j)=vpbp(i,j)
                  ELSE
                     bpr(i,j)=brr(i,j)
                     bpz(i,j)=brz(i,j)
                     bpp(i,j)=brp(i,j)
                  ENDIF
                  
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF (vvbrzphi_flag) THEN
         WRITE(*,*)"Computing vacuum fields without plasma response"
         CALL ipvacuum_bnormal(psilim,bnimn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vvbr,vvbz,vvbp)
      ENDIF
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF (eqbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_eqbrzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_EQBRZPHI: Eq. b field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,5(a16))')"l","r","z","eb_r",
     $        "eb_z","eb_phi"
      
         DO i=0,nr
            DO j=0,nz 
               WRITE(out_unit,'(1x,I2,5(es16.8))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              ebr(i,j),ebz(i,j),ebp(i,j)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (brzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_brzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BRZPHI: Perturbed field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(brr(i,j)),AIMAG(brr(i,j)),
     $              REAL(brz(i,j)),AIMAG(brz(i,j)),
     $              REAL(brp(i,j)),AIMAG(brp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

      ENDIF

      IF (brzphi_flag .AND. vpbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_pbrzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PBRZPHI: Perturbed field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(bpr(i,j)),AIMAG(bpr(i,j)),
     $              REAL(bpz(i,j)),AIMAG(bpz(i,j)),
     $              REAL(bpp(i,j)),AIMAG(bpp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

      ENDIF

      IF (xrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_xrzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_XRZPHI: Displacements in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(xi_r)","imag(xi_r)","real(xi_z)","imag(xi_z)",
     $        "real(xi_phi)","imag(xi_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(xrr(i,j)),AIMAG(xrr(i,j)),
     $              REAL(xrz(i,j)),AIMAG(xrz(i,j)),
     $              REAL(xrp(i,j)),AIMAG(xrp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

      ENDIF

      IF (vbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vbrzphi_n"//
     $        sn//".out","UNKNOWN")

         WRITE(out_unit,*)"IPEC_VBRZPHI: Vacuum field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(vb_r)","imag(vb_r)","real(vb_z)","imag(vb_z)",
     $        "real(vb_phi)","imag(vb_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $              REAL(vbr(i,j)),AIMAG(vbr(i,j)),
     $              REAL(vbz(i,j)),AIMAG(vbz(i,j)),
     $              REAL(vbp(i,j)),AIMAG(vbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
      ENDIF

      IF (vpbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vpbrzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VPBRZPHI: Vacuum field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(vpb_r)","imag(vpb_r)","real(vpb_z)","imag(vpb_z)",
     $        "real(vpb_phi)","imag(vpb_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $              REAL(vpbr(i,j)),AIMAG(vpbr(i,j)),
     $              REAL(vpbz(i,j)),AIMAG(vpbz(i,j)),
     $              REAL(vpbp(i,j)),AIMAG(vpbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
      ENDIF

      IF (vvbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vvbrzphi_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VVBRZPHI: Vacuum field in rzphi grid"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(vvb_r)","imag(vvb_r)","real(vvb_z)","imag(vvb_z)",
     $        "real(vvb_phi)","imag(vvb_phi)"
         
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es16.8))')
     $              vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $              REAL(vvbr(i,j)),AIMAG(vvbr(i,j)),
     $              REAL(vvbz(i,j)),AIMAG(vvbz(i,j)),
     $              REAL(vvbp(i,j)),AIMAG(vvbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
         
      ENDIF

      DEALLOCATE(gdr,gdz,gdl,gdpsi,gdthe,gdphi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbrzphi
   
      END MODULE ipout_mod
