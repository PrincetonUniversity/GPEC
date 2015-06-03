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
c      5. ipout_vsingfld
c      6. ipout_pmodb
c      7. ipout_xbnormal
c      8. ipout_vbnormal
c      9. ipout_xtangent
c     10. ipout_xbrzphi
c     11. ipout_vsbrzphi
c     12. ipout_xbrzphifun
c     13. ipout_arzphifun
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
      USE field_mod

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

      CALL ascii_open(out_unit,"ipec_response_n"//
     $	   TRIM(sn)//".out","UNKNOWN")
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

      b=permeabmats(resp_index,:,:)
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
      REAL(r8) :: respsi,lpsi,rpsi,jarea,thetai
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun,wcfun,
     $     dphi,thetas,units,jacs,jacfac

      COMPLEX(r8), DIMENSION(mpert) :: finmn,foutmn,
     $     fkaxmn,singflx_mn,ftnmn
      COMPLEX(r8), DIMENSION(lmpert) :: lftnmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      REAL(r8), DIMENSION(msing) :: area,j_c,w_c,shear

      COMPLEX(r8), DIMENSION(msing,mpert) :: deltas,delcurs,
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
         resnum=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)
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
         area(ising)=area(ising)-jac*delpsi(mthsurf)/mthsurf
         j_c(ising)=j_c(ising)-jac*delpsi(mthsurf)*
     $        jcfun(mthsurf)/mthsurf
         w_c(ising)=w_c(ising)-0.5*wcfun(mthsurf)/mthsurf

         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0  
         shear(ising)=mfac(resnum)*sq%f1(4)/sq%f(4)**2

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
      singcurs=0
      CALL ipeq_alloc
      DO i=1,mpert
         finmn=0
         finmn(i)=1.0
         CALL ipeq_weight(psilim,finmn,mfac,mpert,1)
         IF (fixed_boundary_flag) THEN
            foutmn=finmn
         ELSE
            foutmn=MATMUL(permeabmats(resp_index,:,:),finmn)
         ENDIF
         singfac=mfac-nn*qlim
         edge_mn=foutmn/(chi1*singfac*twopi*ifac)
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
            
            rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1)) 
            CALL ipeq_sol(rpsi)
            rbwp1mn=bwp1_mn(resnum)

            deltas(ising,i)=rbwp1mn-lbwp1mn
            delcurs(ising,i)=j_c(ising)*deltas(ising,i)*ifac/
     $           (twopi*mfac(resnum))
            singcurs(ising,i)=-delcurs(ising,i)/ifac
            fkaxmn=0
            fkaxmn(resnum)=singcurs(ising,i)/(twopi*nn)

            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
            singbnoflxs(ising,i)=singflx_mn(resnum)/area(ising)
            islandhwids(ising,i)=4*singflx_mn(resnum)/
     $           (twopi*shear(ising)*sq%f(4)*chi1)
         ENDDO
         WRITE(*,'(1x,a16,i4,a22,es10.3)')
     $        "poloidal mode = ",mfac(i),", resonant coupling = ",
     $        SUM(ABS(singbnoflxs(:,i)))/msing
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
         t2mat = MATMUL(singcurs,convmat)
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
      t1v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t2mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s2=s
      t2v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t3mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s3=s
      t3v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))  
c-----------------------------------------------------------------------
c     write matrix.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singcoup_matrix_n"//
     $     TRIM(sn)//".out","UNKNOWN")
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
     $     TRIM(sn)//".out","UNKNOWN")
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
      SUBROUTINE ipout_control(ifile,finmn,foutmn,xspmn,rin,bpin,bin,
     $     rcin,tin,jin,rout,bpout,bout,rcout,tout,svd_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: ifile
      INTEGER, INTENT(IN) :: rin,bpin,bin,rcin,tin,jin,
     $     rout,bpout,bout,rcout,tout
      LOGICAL, INTENT(IN) :: svd_flag
      
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: foutmn,xspmn

      INTEGER :: i,j,i1,i2,i3,ms,itheta,jout
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: vengy,sengy,pengy,area,thetai,scale
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacfac,
     $     sbinfun,sboutfun

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,xinmn,xoutmn,tempmn
      COMPLEX(r8), DIMENSION(lmpert) :: cinmn,coutmn,cawmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,xinfun,xoutfun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2

      REAL(r8), DIMENSION(:,:), POINTER :: dcosmn,dsinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn

      TYPE(spline_type) :: spl       
c-----------------------------------------------------------------------
c     check data_type and read data.
c-----------------------------------------------------------------------
      scale=1.0
      cawmn=0
      tempmn=0
      IF (data_flag) THEN
         WRITE(*,*)"Reading external fields on the control surface"
         ALLOCATE(dcosmn(mmin:mmax,nmin:nmax),
     $        dsinmn(mmin:mmax,nmin:nmax),rawmn(mmin:mmax,nmin:nmax))
         IF (data_type=="surfmn") THEN
            i1=mmax
            i2=mmin
            i3=-1
            scale=1e-4
         ELSE
            i1=mmin
            i2=mmax
            i3=1
         ENDIF
         CALL ascii_open(in_unit,ifile,"old")
                        
         DO i=i1,i2,i3
            IF (data_type=="surfmn") THEN
               READ(in_unit,'(1x,25f12.6)')(dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,'(1x,25f12.6)')(dsinmn(i,j),j=nmin,nmax)
            ELSE IF (data_type=="vac3d") THEN
               READ(in_unit,'(11(1x,e15.8))')(dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,'(11(1x,e15.8))')(dsinmn(i,j),j=nmin,nmax)
            ELSE
               WRITE(message,'(a)')"Can't recognize data format"
               CALL ipec_stop(message)
            ENDIF            
         ENDDO
         
         CALL ascii_close(in_unit)
         rawmn=dcosmn+ifac*dsinmn
         cawmn=rawmn(:,nn)
         DEALLOCATE(dcosmn,dsinmn,rawmn)
      ELSE IF (harmonic_flag) THEN
         DO i=-hmnum,hmnum
            IF ((-mmin+i>=1).AND.(-mmin+i<=lmpert)) THEN
               cawmn(-mmin+i+1)=cosmn(i)+ifac*sinmn(i)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      IF (data_flag .OR. harmonic_flag) THEN
         CALL ipeq_fcoords(psilim,cawmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)             
         binmn=0
         DO i=1,lmpert
            IF ((lmlow-mlow+i>=1).AND.(lmlow-mlow+i<=mpert)) THEN
               binmn(lmlow-mlow+i)=cawmn(i)
            ENDIF
         ENDDO   
c-----------------------------------------------------------------------
c     convert to field if displacement is given.
c-----------------------------------------------------------------------
         IF (displacement_flag) THEN
            CALL ipeq_weight(psilim,binmn,mfac,mpert,5)
            binmn=twopi*ifac*chi1*(mfac-nn*qlim)*binmn
            CALL ipeq_weight(psilim,binmn,mfac,mpert,0)
         ENDIF 
         binmn=binmn*scale
         tempmn=binmn
         CALL ipeq_weight(psilim,tempmn,mfac,mpert,1)              
      ENDIF 
c-----------------------------------------------------------------------
c     get plasma response on the control surface.
c-----------------------------------------------------------------------
      finmn=finmn+tempmn
      IF (fixed_boundary_flag) THEN
         foutmn=finmn
      ELSE
         foutmn=MATMUL(permeabmats(resp_index,:,:),finmn)
      ENDIF
      xspmn=foutmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      binmn=finmn
      boutmn=foutmn
      CALL ipeq_weight(psilim,binmn,mfac,mpert,0)
      CALL ipeq_weight(psilim,boutmn,mfac,mpert,0)
      xinmn=finmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      xoutmn=xspmn
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

      vy = SUM(CONJG(finmn)*MATMUL(temp1,finmn))/4.0
      sy = SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))/4.0
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
      py = SUM(CONJG(foutmn)*MATMUL(temp1,foutmn))/4.0
      pengy = REAL(py)
      WRITE(*,'(1x,a,es10.3)')"Required energy to perturb vacuum = ",
     $     sengy
      WRITE(*,'(1x,a,es10.3)')"Required energy to perturb plasma = ",
     $     pengy
      WRITE(*,'(1x,a,es10.3)')"Amplification factor = ",
     $     sengy/pengy
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_control_n"//
     $     TRIM(sn)//".out","UNKNOWN")
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

      CALL ipeq_bcoords(psilim,binmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL ipeq_bcoords(psilim,boutmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL ipeq_bcoords(psilim,xinmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL ipeq_bcoords(psilim,xoutmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      
      CALL iscdftb(mfac,mpert,binfun,mthsurf,binmn)     
      CALL iscdftb(mfac,mpert,boutfun,mthsurf,boutmn)    
      CALL iscdftb(mfac,mpert,xinfun,mthsurf,xinmn)    
      CALL iscdftb(mfac,mpert,xoutfun,mthsurf,xoutmn)     
      
      CALL ascii_open(out_unit,"ipec_control_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_CONTROL_FUN: "//
     $     "Plasma response for an external perturbation on the "//
     $     "control surface in functions"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a12,I4)')"mthsurf =",mthsurf
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"vacuum energy =",vengy
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"surface energy =",sengy
      WRITE(out_unit,'(1x,a16,1x,es16.8)')"plasma energy =",pengy
      WRITE(out_unit,*)
      
      WRITE(out_unit,*)"jac_type = "//jac_type
      WRITE(out_unit,*)
      
      WRITE(out_unit,'(10(1x,a16))')"r","z",
     $     "real(xin)","imag(xin)","real(xout)","imag(xout)",
     $     "real(bin)","imag(bin)","real(bout)","imag(bout)"
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psilim,theta(itheta),0)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)         
         WRITE(out_unit,'(10(1x,es16.8))')r(itheta),z(itheta),
     $        REAL(xinfun(itheta)),-helicity*AIMAG(xinfun(itheta)),
     $        REAL(xoutfun(itheta)),-helicity*AIMAG(xoutfun(itheta)),
     $        REAL(binfun(itheta)),-helicity*AIMAG(binfun(itheta)),
     $        REAL(boutfun(itheta)),-helicity*AIMAG(boutfun(itheta))
      ENDDO
      CALL ascii_close(out_unit) 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_control
c-----------------------------------------------------------------------
c     subprogram 4. ipout_singfld.
c     compute current and field on rational surfaces.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singfld(egnum,xspmn,spot,
     $     rout,bpout,bout,rcout,tout,svd_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      REAL(r8), INTENT(IN) :: spot
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      LOGICAL, INTENT(IN) :: svd_flag

      INTEGER :: i,itheta,ising
      REAL(r8) :: respsi,lpsi,rpsi,shear,hdist,sbnosurf
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: area,j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: fkaxmn

      REAL(r8), DIMENSION(msing) :: island_hwidth,chirikov,
     $     novf,novs,novi
      COMPLEX(r8), DIMENSION(msing) :: delta,delcur,singcur,
     $     ovf,ovs,ovi
      COMPLEX(r8), DIMENSION(mpert,msing) :: singflx_mn
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing total resonant fields"
      CALL ipeq_alloc
      CALL idcon_build(egnum,xspmn)
      IF (vsbrzphi_flag) ALLOCATE(singbno_mn(mpert,msing))
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
         area(ising)=area(ising)-jac*delpsi(mthsurf)/mthsurf
         j_c(ising)=j_c(ising)-jac*delpsi(mthsurf)*
     $        jcfun(mthsurf)/mthsurf

         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0
         shear=mfac(resnum(ising))*sq%f1(4)/sq%f(4)**2

         lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(lpsi)
         lbwp1mn=bwp1_mn(resnum(ising))
         
         rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(rpsi)
         rbwp1mn=bwp1_mn(resnum(ising))

         delta(ising)=rbwp1mn-lbwp1mn
         delcur(ising)=j_c(ising)*delta(ising)*ifac/
     $        (twopi*mfac(resnum(ising)))
         singcur(ising)=-delcur(ising)/ifac

         fkaxmn=0
         fkaxmn(resnum(ising))=singcur(ising)/(twopi*nn)
         
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
         IF (vsbrzphi_flag) THEN
            singbno_mn(:,ising)=-singflx_mn(:,ising)
            CALL ipeq_weight(respsi,singbno_mn(:,ising),mfac,mpert,0)
         ENDIF
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
         WRITE(*,'(1x,a6,es10.3,a6,f6.3,a25,es10.3)')
     $        "psi = ",singtype(ising)%psifac,
     $        ", q = ",singtype(ising)%q,
     $        ", total resonant field = ",
     $        ABS(singflx_mn(resnum(ising),ising))
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
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
     $        REAL(singcur(ising)),AIMAG(singcur(ising)),
     $        island_hwidth(ising),chirikov(ising)
      ENDDO
      WRITE(out_unit,*)
 
      IF (svd_flag) THEN
         sbnosurf=SQRT(ABS(DOT_PRODUCT(sbno_fun(1:mthsurf),
     $        sbno_fun(1:mthsurf)))/mthsurf/2.0)
         DO ising=1,msing
            ovf(ising)=DOT_PRODUCT(t1v(:,ising),sbno_mn(:))/SQRT(2.0)
            ovs(ising)=DOT_PRODUCT(t2v(:,ising),sbno_mn(:))/SQRT(2.0)
            ovi(ising)=DOT_PRODUCT(t3v(:,ising),sbno_mn(:))/SQRT(2.0)
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
c     subprogram 5. ipout_vsingfld.
c     compute resonant field by coils.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vsingfld(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout

      INTEGER :: ipert,ising,i
      REAL(r8) :: hdist,shear,area
      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: visland_hwidth,vchirikov
      COMPLEX(r8), DIMENSION(:), POINTER :: vcmn

      COMPLEX(r8), DIMENSION(msing) :: vflxmn
      REAL(r8), DIMENSION(0:mthsurf) :: unitfun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing resonant field from coils"
      ALLOCATE(vcmn(cmpert))
      vcmn=0
      DO ising=1,msing
         CALL field_bs_psi(singtype(ising)%psifac,vcmn,2)
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         DO i=1,cmpert
            IF (cmlow-mlow+i==resnum(ising)) THEN
               vflxmn(ising)=vcmn(i)
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     compute half-width of magnetic island.
c-----------------------------------------------------------------------
         shear=mfac(resnum(ising))*
     $        singtype(ising)%q1/singtype(ising)%q**2
         unitfun=1.0
         area=issurfint(unitfun,mthsurf,singtype(ising)%psifac,0,0)
         visland_hwidth(ising)=
     $        SQRT(ABS(4*vflxmn(ising)*area/
     $        (twopi*shear*singtype(ising)%q*chi1)))
         IF (ising==1) THEN 
            hdist=(singtype(ising+1)%psifac-singtype(ising)%psifac)/2.0
         ELSE IF (ising==msing) THEN
            hdist=(singtype(ising)%psifac-singtype(ising-1)%psifac)/2.0
         ELSE IF ((ising/=1).AND.(ising/=msing)) THEN
            hdist=MIN(singtype(ising+1)%psifac-singtype(ising)%psifac,
     $           singtype(ising)%psifac-singtype(ising-1)%psifac)/2.0
         ENDIF
         vchirikov(ising)=visland_hwidth(ising)/hdist
         WRITE(*,'(1x,a6,es10.3,a6,f6.3,a25,es10.3)')
     $        "psi = ",singtype(ising)%psifac,
     $        ", q = ",singtype(ising)%q,
     $        ", vacuum resonant field = ",
     $        ABS(vflxmn(ising))         
         
      ENDDO
      DEALLOCATE(vcmn)
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_vsingfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGFLD: "//
     $     "Resonant fields and islands from coils"
      WRITE(out_unit,*)       
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout 
      WRITE(out_unit,'(1x,a14)')"sweet-spot = 0"
      WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a6,5(1x,a16))')"q","psi",
     $     "real(singflx)","imag(singflx)",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(1x,f6.3,5(1x,es16.8))')
     $        singtype(ising)%q,singtype(ising)%psifac,
     $        REAL(vflxmn(ising)),AIMAG(vflxmn(ising)),
     $        visland_hwidth(ising),vchirikov(ising)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_vsingfld
c-----------------------------------------------------------------------
c     subprogram 6. ipout_pmodb.
c     compute perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_pmodb(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,itheta,iindex,cstep
      REAL(r8) :: ileft,jac

      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,equilbfun
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,
     $     divxprp_mn,curv_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmz_fun,xvt_fun,xvz_fun,xmt_fun,xsp1_fun
      COMPLEX(r8), DIMENSION(mstep,mpert) :: eulbparmns,lagbparmns,
     $     divxprpmns,curvmns,eulbparmout,lagbparmout,divxprpmout,
     $     curvmout
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) ::eulbparfun,lagbparfun,
     $     divxprpfun,curvfun,eulbparfout,lagbparfout,divxprpfout,
     $     curvfout

      REAL(r8), DIMENSION(:), POINTER :: psis,ches,chex
      COMPLEX(r8), DIMENSION(:,:), POINTER ::chea,chelagbmns,chelagbmout
      COMPLEX(r8), DIMENSION(mpert) :: chelagb_mn

      TYPE(cspline_type) :: cspl,chespl 

c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing |b| and del.x_prp"
      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta

      CALL idcon_build(egnum,xspmn)
      CALL ipeq_alloc

      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a32)')
     $        "volume = ",iindex,"% |b| and del.x_prp computations"
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))         
         CALL ipeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b variations in hamada.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         CALL iscdftb(mfac,mpert,xmt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xmz_fun,mthsurf,xmz_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)

         CALL spline_eval(sq,psifac(istep),1)

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

            eulbparfun(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(rzphi%f(4)*eqfun%f(1))
            lagbparfun(istep,itheta)=
     $           eulbparfun(istep,itheta)+
     $           xsp_fun(itheta)*eqfun%fx(1)+
     $           xms_fun(itheta)/(chi1*sq%f(4))*eqfun%fy(1)
            divxprpfun(istep,itheta)=xsp1_fun(itheta)+
     $           (jac1/jac)*xsp_fun(itheta)+
     $           cspl%f1(1)/jac-(twopi*ifac*nn)*cspl%f(2)/jac
            divxprpfun(istep,itheta)=-eqfun%f(1)*
     $           divxprpfun(istep,itheta)
            curvfun(istep,itheta)=-xsp_fun(itheta)/eqfun%f(1)*sq%f1(2)-
     $           (xsp_fun(itheta)*eqfun%fx(1)+eqfun%fy(1)*
     $           (xms_fun(itheta)/(chi1*sq%f(4))+
     $           xmz_fun(itheta)/(jac*sq%f(4))-
     $           (chi1/(jac*eqfun%f(1)))**2*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))))
            equilbfun(istep,itheta) = eqfun%f(1)
            dphi(itheta) = rzphi%f(3)
         ENDDO
         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
     $        mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
     $        mthsurf,lagbpar_mn)
         CALL iscdftf(mfac,mpert,divxprpfun(istep,:),
     $        mthsurf,divxprp_mn)
         CALL iscdftf(mfac,mpert,curvfun(istep,:),
     $        mthsurf,curv_mn)
c-----------------------------------------------------------------------
c     Save working coordinate verions for pent files.
c-----------------------------------------------------------------------
         lagbparmns(istep,:) = lagbpar_mn
         divxprpmns(istep,:) = divxprp_mn
         curvmns(istep,:) = curv_mn
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
            CALL ipeq_bcoords(psifac(istep),curv_mn,
     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
         ENDIF
         eulbparmout(istep,:)=eulbpar_mn
         lagbparmout(istep,:)=lagbpar_mn
         divxprpmout(istep,:)=divxprp_mn
         curvmout(istep,:)  =curv_mn
         eulbparfout(istep,:)=eulbparfun(istep,:)*EXP(ifac*nn*dphi)
         lagbparfout(istep,:)=lagbparfun(istep,:)*EXP(ifac*nn*dphi)
         divxprpfout(istep,:)=divxprpfun(istep,:)*EXP(ifac*nn*dphi)
         curvfout(istep,:)=curvfun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      CALL cspline_dealloc(cspl)

      CALL ascii_open(out_unit,"ipec_pmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_PMODB: "//
     $     "Components in perturbed mod b and del.x_prp"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     

      WRITE(out_unit,'(1x,a16,1x,a4,10(1x,a16))')"psi","m",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $     "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $     "imag(Bkxprp)"
      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es16.8,1x,I4,8(1x,es16.8))')
     $           psifac(istep),mfac(ipert),
     $           REAL(eulbparmout(istep,ipert)),
     $           AIMAG(eulbparmout(istep,ipert)),
     $           REAL(lagbparmout(istep,ipert)),
     $           AIMAG(lagbparmout(istep,ipert)),
     $           REAL(-divxprpmout(istep,ipert)),
     $           AIMAG(-divxprpmout(istep,ipert)),
     $           REAL(-curvmout(istep,ipert)),
     $           AIMAG(-curvmout(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      IF (chebyshev_flag) THEN
         WRITE(*,*)"Computing chebyshev for pmodb"
         ALLOCATE(chea(mpert,0:nche),chex(0:mstep))
         CALL cspline_alloc(chespl,mstep,1)
         chespl%xs=psifac
         chex=2.0*(psifac-0.5)
         DO ipert=1,mpert
            DO i=0,nche
               chespl%fs(:,1)=cos(REAL(i,r8)*acos(chex))/
     $              sqrt(1.0-chex**2.0)*lagbparmns(:,ipert)*4.0/pi
               CALL cspline_fit(chespl,"extrap")
               CALL cspline_int(chespl)
               chea(ipert,i)=chespl%fsi(mstep,1)
            ENDDO
         ENDDO
	 
	 chea(:,0)=chea(:,0)/2.0

         CALL ascii_open(out_unit,"ipec_pmodb_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_CHEBYSHEV: "//
     $        "Chebyshev coefficients for perturbed mod b"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nche =",nche,"mpert =",mpert
         WRITE(out_unit,*)     
         
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,a6,1x,I4)')"mfac =",mfac(ipert)
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a6,2(1x,a16))')"chea",
     $           "real(chea)","imag(chea)"
            DO i=0,nche
               WRITE(out_unit,'(1x,I6,2(1x,es16.8))')i,
     $              REAL(chea(ipert,i)),AIMAG(chea(ipert,i))
            ENDDO
            WRITE(out_unit,*)
         ENDDO
         CALL ascii_close(out_unit)

         cstep=200
         ALLOCATE(psis(cstep),ches(cstep))
         ALLOCATE(chelagbmns(cstep,mpert))
         chelagbmns=0
         psis=(/(istep,istep=0,cstep-1)/)/REAL(cstep-1,r8)*
     $        (psilim-psilow)+psilow
         ches=2.0*(psis-0.5)
         
         DO ipert=1,mpert
            DO i=0,nche
               chelagbmns(:,ipert)=chelagbmns(:,ipert)+
     $              chea(ipert,i)*cos(REAL(i,r8)*acos(ches))
            ENDDO
         ENDDO

         CALL ascii_open(out_unit,"ipec_chepmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_CHEPMODB: "//
     $        "Reconstructed perturbed mod b by chebyshev"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "cstep =",cstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(1x,a16,1x,a4,2(1x,a16))')"psi","m",
     $        "real(lagb)","imag(lagb)"
         DO istep=1,cstep
            DO ipert=1,mpert
               WRITE(out_unit,'(1x,es16.8,1x,I4,2(1x,es16.8))')
     $              psis(istep),mfac(ipert),
     $              REAL(chelagbmns(istep,ipert)),
     $              AIMAG(chelagbmns(istep,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         CALL cspline_dealloc(chespl)
      ENDIF

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_pmodb_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_FUN: "//
     $        "Components in perturbed mod b in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8,a13,es16.8)')
     $        "jac_out = ",jac_out,"R0 = ",ro
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(17(1x,a16))')"psi","theta","r","z",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $        "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $        "imag(Bkxprp)","equilb","dequilbdpsi","dequilbdtheta",
     $        "real(xms)","imag(xms)"
         DO istep=1,mstep
            CALL ipeq_sol(psifac(istep))
            CALL ipeq_contra(psifac(istep))         
            CALL ipeq_cova(psifac(istep))
            CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
            DO itheta=0,mthsurf
               CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
               WRITE(out_unit,'(17(1x,es16.8))')
     $              psifac(istep),theta(itheta),
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(eulbparfout(istep,itheta)),
     $              -helicity*AIMAG(eulbparfout(istep,itheta)),
     $              REAL(lagbparfout(istep,itheta)),
     $              -helicity*AIMAG(lagbparfout(istep,itheta)),
     $              REAL(-divxprpfout(istep,itheta)),
     $              -helicity*AIMAG(-divxprpfout(istep,itheta)),
     $              REAL(-curvfout(istep,itheta)),
     $              -helicity*AIMAG(-curvfout(istep,itheta)),
     $              equilbfun(istep,itheta),
     $              eqfun%fx(1),eqfun%fy(1),
     $              REAL(xms_fun(itheta)),
     $              -helicity*AIMAG(xms_fun(itheta))
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
     $              REAL(REAL(eulbparmout(istep,ipert)),4),
     $              REAL(AIMAG(eulbparmout(istep,ipert)),4),
     $              REAL(REAL(lagbparmout(istep,ipert)),4),
     $              REAL(AIMAG(lagbparmout(istep,ipert)),4)     
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
         WRITE(bin_2d_unit)REAL(REAL(eulbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(eulbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(lagbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(lagbparfout(9:mstep,1:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF

      IF(ntv_flag) THEN
         CALL bin_open(bin_2d_unit,"ipec_order1_fun_n"//
     $        TRIM(sn)//".bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)mstep,mthsurf
         WRITE(bin_2d_unit)lagbparfun/equilbfun
         WRITE(bin_2d_unit)-(lagbparfun-curvfun)/equilbfun
         CALL bin_close(bin_2d_unit)

         CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
         CALL bin_open(bin_2d_unit,"ipec_order1_n"//TRIM(sn)//".bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)mstep,mpert
         WRITE(bin_2d_unit)lagbparmns/eqfun%f(1)
         WRITE(bin_2d_unit)-(lagbparmns-curvmns)/eqfun%f(1)
         CALL bin_close(bin_2d_unit)
      ENDIF

      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_pmodb
c-----------------------------------------------------------------------
c     subprogram 7. ipout_xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbnormal(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax,area


      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs,dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun

      COMPLEX(r8), DIMENSION(mstep,mpert) :: xmns,ymns,
     $     xnomns,bnomns,bwpmns
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn
      COMPLEX(r8), DIMENSION(mstep,lmpert) :: pwpmns
      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,psis,rvecs,zvecs
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss,
     $     xnofuns,bnofuns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing x and b normal components"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))

         area=0
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
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
            area=area+jac*delpsi(itheta)/mthsurf
         ENDDO
         area=area-jac*delpsi(mthsurf)/mthsurf
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         xnofuns(istep,:)=xwp_fun/(jacs*delpsi)
         bnofuns(istep,:)=bwp_fun/(jacs*delpsi)
         CALL iscdftf(mfac,mpert,xnofuns(istep,:),mthsurf,xno_mn)
         CALL iscdftf(mfac,mpert,bnofuns(istep,:),mthsurf,bno_mn)
         bwp_mn=bwp_mn/area

         IF (bwp_pest_flag) THEN
            pwpmn=0
            DO i=1,mpert
               IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
                  pwpmn(mlow-lmlow+i)=bno_mn(i)
               ENDIF
            ENDDO
            CALL ipeq_bcoords(psifac(istep),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)
            pwpmns(istep,:)=pwpmn
         ENDIF            

         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            bwp_mn=bno_mn
            CALL ipeq_bcoords(psifac(istep),xno_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),bno_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),bwp_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,1)            
         ENDIF
         xnomns(istep,:)=xno_mn
         bnomns(istep,:)=bno_mn
         bwpmns(istep,:)=bwp_mn
         xnofuns(istep,:)=xnofuns(istep,:)*EXP(ifac*nn*dphi)
         bnofuns(istep,:)=bnofuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBNROMAL: "//
     $     "Normal components of displacement and field"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $     "real(xno)","imag(xno)","real(bno)","imag(bno)",
     $     "real(bwp)","imag(bwp)"

      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(2(1x,es16.8),1x,I4,6(1x,es16.8))')
     $           psifac(istep),qfac(istep),mfac(ipert),
     $           REAL(xnomns(istep,ipert)),AIMAG(xnomns(istep,ipert)),
     $           REAL(bnomns(istep,ipert)),AIMAG(bnomns(istep,ipert)),
     $           REAL(bwpmns(istep,ipert)),AIMAG(bwpmns(istep,ipert))        
         ENDDO
      ENDDO

      IF (bwp_pest_flag) THEN
         CALL ascii_open(out_unit,"ipec_bnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BNORMAL_PEST: "//
     $        "Normal components of field in pest coordinates"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         
         DO istep=1,mstep
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(1x,es16.8),1x,I4,6(1x,es16.8))')
     $              psifac(istep),qfac(istep),lmfac(ipert),
     $              REAL(pwpmns(istep,ipert)),AIMAG(pwpmns(istep,ipert))        
            ENDDO
         ENDDO
      ENDIF

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbnormal_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBNROMAL_FUN: "//
     $        "Normal components of displacement and field in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(10(1x,a16))')"psi","theta","r","z","rvec",
     $        "zvec","real(xno)","imag(xno)","real(bno)","imag(bno)"
         
         DO istep=1,mstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(10(1x,es16.8))')
     $              psifac(istep),theta(itheta),
     $              rs(istep,itheta),zs(istep,itheta),
     $              rvecs(istep,itheta),zvecs(istep,itheta),
     $              REAL(xnofuns(istep,itheta)),
     $              -helicity*AIMAG(xnofuns(istep,itheta)),
     $              REAL(bnofuns(istep,itheta)),
     $              -helicity*AIMAG(bnofuns(istep,itheta))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "xbnormal.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
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
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:mstep,0:mthsurf),4),
     $        REAL(zs(9:mstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(xnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(xnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(bnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(bnofuns(9:mstep,0:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
         
         ximax=MAXVAL(ABS(xnofuns))
         CALL bicube_eval(rzphi,psilim,theta(0),0)
         rmax=SQRT(rzphi%f(1))
         xnofuns=xnofuns/ximax*rmax/6.0

         WRITE(*,'(1x,a,es10.3)')"Maximum displacement = ",ximax
         WRITE(*,'(1x,a,es10.3)')"Scale factor for 2d plots = ",
     $        rmax/(ximax*6.0)

         rss=CMPLX(rs,rs)+xnofuns*rvecs
         zss=CMPLX(zs,zs)+xnofuns*zvecs
         DO itheta=0,mthsurf
            psis(:,itheta)=psifac(:)
         ENDDO

         CALL bin_open(bin_2d_unit,
     $        "pflux_re_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(REAL(rss(9:mstep,0:mthsurf)),4),
     $        REAL(REAL(zss(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psis(9:mstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         CALL bin_open(bin_2d_unit,
     $        "pflux_im_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(rss(9:mstep,0:mthsurf)),4),
     $        REAL(-helicity*AIMAG(zss(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psis(9:mstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         DO istep=1,mstep
            xmns(istep,:)=mfac
            ymns(istep,:)=psifac(istep)
         ENDDO
         CALL bin_open(bin_2d_unit,"bnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mpert-1
         WRITE(bin_2d_unit)REAL(xmns(9:mstep,:),4),
     $        REAL(ymns(9:mstep,:),4)
         WRITE(bin_2d_unit)REAL(ABS(bwpmns(9:mstep,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF

      IF (flux_flag) THEN
         IF (.NOT. bin_2d_flag) THEN
            ximax=MAXVAL(ABS(xnofuns))
            CALL bicube_eval(rzphi,psilim,theta(0),0)
            rmax=SQRT(rzphi%f(1))
            xnofuns=xnofuns/ximax*rmax/6.0
         ENDIF
            
         rss=xnofuns*rvecs
         zss=xnofuns*zvecs
         DO itheta=0,mthsurf
            psis(:,itheta)=psifac(:)
         ENDDO
         
         CALL ascii_open(out_unit,"ipec_xbnormal_flux_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBNROMAL_FLUX: "//
     $        "Perturbed flux surfaces"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,'(1x,a12,1x,es16.8)')"scale =",rmax/(ximax*6.0)
         WRITE(out_unit,*)  
         WRITE(out_unit,'(7(1x,a16))')"r","z","real(r)","imag(r)",
     $        "real(z)","imag(z)","psi"   
         
         DO istep=1,mstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(7(1x,es16.8))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(rss(istep,itheta)),AIMAG(rss(istep,itheta)),
     $              REAL(zss(istep,itheta)),AIMAG(zss(istep,itheta)),
     $              psis(istep,itheta)
            ENDDO
         ENDDO
         
         CALL ascii_close(out_unit)
      ENDIF         
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbnormal
c-----------------------------------------------------------------------
c     subprogram 8. ipout_vbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vbnormal(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      
      INTEGER :: ipsi,ipert,i
      REAL(r8), DIMENSION(0:cmpsi) :: psi
      COMPLEX(r8), DIMENSION(:), POINTER :: vcmn

      COMPLEX(r8), DIMENSION(mpert) :: vwpmn
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn
      REAL(r8), DIMENSION(cmpsi) :: qs
      REAL(r8), DIMENSION(cmpsi,mpert) :: xmns,ymns
      COMPLEX(r8), DIMENSION(cmpsi,mpert) :: vmn
      COMPLEX(r8), DIMENSION(cmpsi,lmpert) :: pmn
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing b normal components from coils"

      psi=(/(ipsi,ipsi=0,cmpsi)/)/REAL(cmpsi,r8)
      psi=psilow+(psilim-psilow)*SIN(psi*pi/2)**2
      ALLOCATE(vcmn(cmpert))
      vmn=0
      pmn=0

      DO ipsi=1,cmpsi
         CALL spline_eval(sq,psi(ipsi),0)
         qs(ipsi)=sq%f(4)

         IF (bwp_pest_flag) THEN
            CALL field_bs_psi(psi(ipsi),vcmn,0)
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  pmn(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
            pwpmn=0
            pwpmn=pmn(ipsi,:)
            CALL ipeq_bcoords(psi(ipsi),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)             
            pmn(ipsi,:)=pwpmn            
         ENDIF
         
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL field_bs_psi(psi(ipsi),vcmn,0)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  vmn(ipsi,cmlow-mlow+i)=vcmn(i)
               ENDIF
            ENDDO
            vwpmn=0
            vwpmn=vmn(ipsi,:)
            CALL ipeq_bcoords(psi(ipsi),vwpmn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,1)  
            vmn(ipsi,:)=vwpmn         
         ELSE
            CALL field_bs_psi(psi(ipsi),vcmn,2)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  vmn(ipsi,cmlow-mlow+i)=vcmn(i)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(vcmn)

      CALL ascii_open(out_unit,"ipec_vbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_VBNROMAL: "//
     $     "Normal components of coil field"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mpsi =",cmpsi,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,2(1x,a16))')"psi","q","m",
     $     "real(vb)","imag(vb)"

      DO ipsi=1,cmpsi
         DO ipert=1,mpert
            WRITE(out_unit,'(2(1x,es16.8),1x,I4,2(1x,es16.8))')
     $           psi(ipsi),qs(ipsi),mfac(ipert),
     $           REAL(vmn(ipsi,ipert)),AIMAG(vmn(ipsi,ipert)) 
         ENDDO
      ENDDO

      IF (bwp_pest_flag) THEN
         CALL ascii_open(out_unit,"ipec_vbnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VBNROMAL_PEST: "//
     $        "Normal components of coil field in pest coordinates"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mpsi =",cmpsi,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,2(1x,a16))')"psi","q","m",
     $        "real(vb)","imag(vb)"
         
         DO ipsi=1,cmpsi
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(1x,es16.8),1x,I4,2(1x,es16.8))')
     $              psi(ipsi),qs(ipsi),lmfac(ipert),
     $              REAL(pmn(ipsi,ipert)),AIMAG(pmn(ipsi,ipert)) 
            ENDDO
         ENDDO
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "vbnormal.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO ipsi=1,cmpsi
               WRITE(bin_unit)REAL(psi(ipsi),4),
     $              REAL(REAL(vmn(ipsi,ipert)),4),
     $              REAL(AIMAG(vmn(ipsi,ipert)),4)
               
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)

         DO ipsi=1,cmpsi
            xmns(ipsi,:)=mfac
            ymns(ipsi,:)=psi(ipsi)
         ENDDO
         CALL bin_open(bin_2d_unit,"vbnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)cmpsi-1,mpert-1
         WRITE(bin_2d_unit)REAL(xmns(1:cmpsi,:),4),
     $        REAL(ymns(1:cmpsi,:),4)
         WRITE(bin_2d_unit)REAL(ABS(vmn(1:cmpsi,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_vbnormal
c-----------------------------------------------------------------------
c     subprogram 9. ipout_xbtangent.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbtangent(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,psis,
     $     rvecs,zvecs,vecs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,bs,dphi
      COMPLEX(r8), DIMENSION(mstep,mpert) :: xmns,ymns,xtamns,btamns
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwt_fun,xvt_fun,xvz_fun,
     $     bwt_fun,bvt_fun,bvz_fun,xta_fun,bta_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss,
     $     xtafuns,btafuns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing x and b tangential components"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))

         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),0)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jacs(itheta)=rzphi%f(4)
            bs(itheta)=eqfun%f(1)
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            rvecs(istep,itheta)=v(2,1)*cos(eta)-v(2,2)*sin(eta)
            zvecs(istep,itheta)=v(2,1)*sin(eta)+v(2,2)*cos(eta)
         ENDDO
         vecs(istep,:)=sqrt(rvecs(istep,:)**2+zvecs(istep,:)**2)
         rvecs(istep,:)=rvecs(istep,:)/vecs(istep,:)
         zvecs(istep,:)=zvecs(istep,:)/vecs(istep,:)
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         xta_fun=xwt_fun/jacs-(chi1/(jacs*bs))**2*(xvt_fun+q*xvz_fun)
         bta_fun=bwt_fun/jacs-(chi1/(jacs*bs))**2*(bvt_fun+q*bvz_fun)
         xtafuns(istep,:)=xta_fun*vecs(istep,:)
         btafuns(istep,:)=bta_fun*vecs(istep,:)
         CALL iscdftf(mfac,mpert,xtafuns(istep,:),mthsurf,xta_mn)
         CALL iscdftf(mfac,mpert,btafuns(istep,:),mthsurf,bta_mn)

         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL ipeq_bcoords(psifac(istep),xta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),bta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
         ENDIF
         xtamns(istep,:)=xta_mn
         btamns(istep,:)=bta_mn
         xtafuns(istep,:)=xtafuns(istep,:)*EXP(ifac*nn*dphi)
         btafuns(istep,:)=btafuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbtangent_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBNROMAL: "//
     $     "Tangential components of displacement and field"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $     "real(xta)","imag(xta)","real(bta)","imag(bta)"

      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(2(1x,es16.8),1x,I4,4(1x,es16.8))')
     $           psifac(istep),qfac(istep),mfac(ipert),
     $           REAL(xtamns(istep,ipert)),AIMAG(xtamns(istep,ipert)),
     $           REAL(btamns(istep,ipert)),AIMAG(btamns(istep,ipert))
         ENDDO
      ENDDO

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbtangent_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBTANGENT_FUN: "//
     $        "Tangential components of displacement "//
     $        "and field in functions"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(8(1x,a16))')"r","z","rvec","zvec",
     $        "real(xta)","imag(xta)","real(bta)","imag(bta)"
         
         DO istep=1,mstep
            DO itheta=0,mthsurf
               WRITE(out_unit,'(8(1x,es16.8))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              rvecs(istep,itheta),zvecs(istep,itheta),
     $              REAL(xtafuns(istep,itheta)),
     $              AIMAG(xtafuns(istep,itheta)),
     $              REAL(btafuns(istep,itheta)),
     $              AIMAG(btafuns(istep,itheta))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "xbtangent.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
     $              REAL(REAL(xtamns(istep,ipert)),4),
     $              REAL(AIMAG(xtamns(istep,ipert)),4),
     $              REAL(REAL(btamns(istep,ipert)),4),
     $              REAL(AIMAG(btamns(istep,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"xbtangent_2d.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:mstep,0:mthsurf),4),
     $        REAL(zs(9:mstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(xtafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(xtafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(btafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(btafuns(9:mstep,0:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF         
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbtangent
c-----------------------------------------------------------------------
c     subprogram 10. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphi(egnum,xspmn,nr,nz,bnimn,bnomn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnimn,bnomn

      INTEGER :: i,j,k,l,ipert,iindex,np
      REAL(r8) :: mid,btlim,rlim,ileft,delr,delz,cha,chb,chc,chd

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz,ebr,ebz,ebp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: xrr,xrz,xrp,brr,brz,brp,
     $     bpr,bpz,bpp,vbr,vbz,vbp,vpbr,vpbz,vpbp,vvbr,vvbz,vvbp,
     $     btr,btz,btp,vcbr,vcbz,vcbp,xcr,xcz,xcp,bcr,bcz,bcp

      REAL(r8), DIMENSION(:), POINTER :: chex,chey
      COMPLEX(r8), DIMENSION(:,:), POINTER :: chear,cheaz,
     $     chxar,chxaz
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
      xcr = 0
      xcz = 0
      xcp = 0
      brr = 0
      brz = 0
      brp = 0
      btr = 0
      btz = 0
      btp = 0
      bcr = 0
      bcz = 0
      bcp = 0
      bpr = 0
      bpz = 0
      bpp = 0
      vbr = 0
      vbz = 0
      vbp = 0
      vcbr = 0
      vcbz = 0
      vcbp = 0
      vpbr = 0
      vpbz = 0
      vpbp = 0
      vvbr = 0
      vvbz = 0
      vvbp = 0

      CALL idcon_build(egnum,xspmn)

      IF (eqbrzphi_flag) THEN
         WRITE(*,*)"Computing equilibrium magnetic fields"
         ! evaluate f value for vacuum
         mid = 0.0
         CALL spline_eval(sq,psilim,0)
         CALL bicube_eval(rzphi,psilim,mid,0)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*rzphi%f(2)
         rlim=ro+rfac*COS(eta)
         btlim = abs(sq%f(1))/(twopi*rlim)
         DO i=0,nr
            DO j=0,nz
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),1)
               ebr(i,j) = -psi_in%fy(1)/gdr(i,j)*psio
               ebz(i,j) = psi_in%fx(1)/gdr(i,j)*psio
               IF (gdl(i,j) == 1) THEN  
                  CALL spline_eval(sq,gdpsi(i,j),0)
                  ebp(i,j) = abs(sq%f(1))/(twopi*gdr(i,j))
               ELSE
                  ebp(i,j) = btlim*rlim/gdr(i,j)  
               ENDIF
            ENDDO   
         ENDDO

      IF(ipd>0)THEN
         ebr=-ebr
         ebz=-ebz
      ENDIF
      IF(btd<0)ebp=-ebp
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
               IF (gdl(i,j)==1) THEN
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
     $                    REAL(brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(brr(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                     brz(i,j)=
     $                    REAL(brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(brz(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                     brp(i,j)=
     $                    REAL(brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(brp(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                  ENDIF
                  IF (xrzphi_flag) THEN
                     xrr(i,j)=
     $                    REAL(xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(xrr(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                     xrz(i,j)=
     $                    REAL(xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(xrz(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                     xrp(i,j)=
     $                    REAL(xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j)))-
     $                    helicity*ifac*AIMAG(xrp(i,j)*EXP(-twopi*ifac*
     $                    nn*gdphi(i,j)))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL ipeq_dealloc

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         WRITE(*,*)"Computing vacuum fields by surface currents"
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
         IF (helicity<0) THEN
            vbr=CONJG(vbr)
            vbz=CONJG(vbz)
            vbp=CONJG(vbp)
         ENDIF
         DO i=0,nr
            DO j=0,nz
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  brr(i,j)=vbr(i,j)
                  brz(i,j)=vbz(i,j)
                  brp(i,j)=vbp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF
      
      IF (divzero_flag) THEN
         CALL ipeq_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF
      IF (div_flag) THEN
         CALL ipdiag_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF

      IF (brzphi_flag) THEN
         WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF

         IF (coil_flag) THEN
            WRITE(*,*)"Computing vacuum fields by coils"
            np=nn*16
            CALL field_bs_rzphi(nr,nz,np,gdr,gdz,vcbr,vcbz,vcbp)
         ENDIF
         
         DO i=0,nr
            DO j=0,nz                  
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  bpr(i,j)=vpbr(i,j)
                  bpz(i,j)=vpbz(i,j)
                  bpp(i,j)=vpbp(i,j)
                  btr(i,j)=vpbr(i,j)+vcbr(i,j)
                  btz(i,j)=vpbz(i,j)+vcbz(i,j)
                  btp(i,j)=vpbp(i,j)+vcbp(i,j)
               ELSE
                  bpr(i,j)=brr(i,j)-vcbr(i,j)
                  bpz(i,j)=brz(i,j)-vcbz(i,j)
                  bpp(i,j)=brp(i,j)-vcbp(i,j)
                  btr(i,j)=brr(i,j)
                  btz(i,j)=brz(i,j)
                  btp(i,j)=brp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
      ENDIF

      IF (chebyshev_flag) THEN
         WRITE(*,*)"Computing chebyshev for xbrzphi"
         ALLOCATE(chex(0:nr),chey(0:nz),
     $        chear(0:nchr,0:nchz),cheaz(0:nchr,0:nchz),
     $        chxar(0:nchr,0:nchz),chxaz(0:nchr,0:nchz))
         delr = gdr(1,0)-gdr(0,0)
         delz = gdz(0,1)-gdz(0,0)
         cha = (gdr(nr,0)-gdr(0,0))/2.0+delr
         chb = (gdr(nr,0)+gdr(0,0))/2.0
         chc = (gdz(0,nz)-gdz(0,0))/2.0+delz
         chd = (gdz(0,nz)+gdz(0,0))/2.0
         chex = (gdr(:,0)-chb)/cha
         chey = (gdz(0,:)-chd)/chc
         chear=0
         cheaz=0
         chxar=0
         chxaz=0

         DO i=0,nchr
            DO j=0,nchz
               DO k=0,nr
                  DO l=0,nz
                     chear(i,j)=chear(i,j)+btr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     cheaz(i,j)=cheaz(i,j)+btz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                     chxar(i,j)=chxar(i,j)+xrr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     chxaz(i,j)=chxaz(i,j)+xrz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                  ENDDO
               ENDDO
               chear(i,j)=chear(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               cheaz(i,j)=cheaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxar(i,j)=chxar(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxaz(i,j)=chxaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
            ENDDO
         ENDDO
         chear(0,:)=chear(0,:)/2.0
         cheaz(0,:)=cheaz(0,:)/2.0        
         chear(:,0)=chear(:,0)/2.0
         cheaz(:,0)=cheaz(:,0)/2.0  
         chxar(0,:)=chxar(0,:)/2.0
         chxaz(0,:)=chxaz(0,:)/2.0        
         chxar(:,0)=chxar(:,0)/2.0
         chxaz(:,0)=chxaz(:,0)/2.0   

         CALL ascii_open(out_unit,"ipec_brzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for brzphi field"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,'(4(1x,a4,f12.3))')
     $        "a =",cha,"b =",chb,"c =",chc,"d =",chd
         WRITE(out_unit,*)   
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chear","cheaz",
     $        "real(chear)","imag(chear)","real(cheaz)","imag(cheaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(1x,es16.8))')i,j,
     $              REAL(chear(i,j)),AIMAG(chear(i,j)),
     $              REAL(cheaz(i,j)),AIMAG(cheaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

         CALL ascii_open(out_unit,"ipec_xrzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_XRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for displacement"
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chxar","chxaz",
     $        "real(chxar)","imag(chxar)","real(chxaz)","imag(chxaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(1x,es16.8))')i,j,
     $              REAL(chxar(i,j)),AIMAG(chxar(i,j)),
     $              REAL(chxaz(i,j)),AIMAG(chxaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

         WRITE(*,*)"Recontructing xbrzphi by chebyshev"

         DO i=0,nr
            DO j=0,nz
               DO k=0,nchr
                  DO l=0,nchz
                     bcr(i,j)=bcr(i,j)+chear(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcz(i,j)=bcz(i,j)+cheaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcp(i,j)=bcp(i,j)+chear(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+cheaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                     xcr(i,j)=xcr(i,j)+chxar(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcz(i,j)=xcz(i,j)+chxaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcp(i,j)=xcp(i,j)+chxar(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+chxaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                  ENDDO
               ENDDO
               bcp(i,j)=bcp(i,j)+bcr(i,j)
               xcp(i,j)=xcp(i,j)+xcr(i,j)
            ENDDO
         ENDDO
         bcp=bcp*(-ifac/nn)
         xcp=xcp*(-ifac/nn)

         DEALLOCATE(chex,chey,chear,cheaz,chxar,chxaz)

      ENDIF

      IF (pbrzphi_flag) THEN
         WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF
      ENDIF

      IF (vvbrzphi_flag) THEN
         WRITE(*,*)"Computing vacuum fields without plasma response"
         CALL ipvacuum_bnormal(psilim,bnimn,nr,nz)
         CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vvbr,vvbz,vvbp)
         IF (helicity<0) THEN
            vvbr=CONJG(vvbr)
            vvbz=CONJG(vvbz)
            vvbp=CONJG(vvbp)
         ENDIF
      ENDIF  
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF (eqbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_eqbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
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
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BRZPHI: Total perturbed field"
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
     $              REAL(btr(i,j)),AIMAG(btr(i,j)),
     $              REAL(btz(i,j)),AIMAG(btz(i,j)),
     $              REAL(btp(i,j)),AIMAG(btp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"ipec_chebrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CEHBRZPHI: Total perturbed field"
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es16.8))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(bcr(i,j)),AIMAG(bcr(i,j)),
     $                 REAL(bcz(i,j)),AIMAG(bcz(i,j)),
     $                 REAL(bcp(i,j)),AIMAG(bcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         CALL ascii_open(out_unit,"ipec_pbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PBRZPHI: Perturbed field by plasma"
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

         IF (coil_flag) THEN
            CALL ascii_open(out_unit,"ipec_cbrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CBRZPHI: External field by coils"
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es16.8))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(vcbr(i,j)),AIMAG(vcbr(i,j)),
     $                 REAL(vcbz(i,j)),AIMAG(vcbz(i,j)),
     $                 REAL(vcbp(i,j)),AIMAG(vcbp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VBRZPHI: Total perturbed field "//
     $        "by surface currents"
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

      IF (xrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_xrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_XRZPHI: Displacement"
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

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"ipec_chexrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CEHXRZPHI: Displacement"
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $           "real(x_r)","imag(x_r)","real(x_z)","imag(x_z)",
     $           "real(x_phi)","imag(x_phi)"
            
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es16.8))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(xcr(i,j)),AIMAG(xcr(i,j)),
     $                 REAL(xcz(i,j)),AIMAG(xcz(i,j)),
     $                 REAL(xcp(i,j)),AIMAG(xcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

      IF (pbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_pbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
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
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VVBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         
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

      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag) DEALLOCATE(gdr,gdz,gdl,gdpsi,gdthe,gdphi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbrzphi
c-----------------------------------------------------------------------
c     subprogram 10. ipout_vsbrzphi.
c     write brzphi components restored by removing shielding currents.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vsbrzphi(snum,nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: snum,nr,nz

      INTEGER :: i,j,ipert,iindex
      REAL(r8) :: mid,bt0,ileft

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: vbr,vbz,vbp
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      vbr = 0
      vbz = 0
      vbp = 0

      IF (snum<10) THEN
         WRITE(UNIT=ss,FMT='(I1)')snum
         ss=TRIM(ADJUSTL(ss))
      ELSE
         WRITE(UNIT=ss,FMT='(I2)')snum
      ENDIF

      WRITE(*,*)"Computing vacuum fields by "//
     $     TRIM(ss)//"th resonant field"
      CALL ipvacuum_bnormal(singtype(snum)%psifac,
     $     singbno_mn(:,snum),nr,nz)
      CALL mscfld(wv,mpert,mthvac,mthvac,nfm2,nths2,complex_flag,
     $     nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
      IF (helicity<0) THEN
         vbr=CONJG(vbr)
         vbz=CONJG(vbz)
         vbp=CONJG(vbp)
      ENDIF
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_vsbrzphi_n"//
     $     TRIM(sn)//"_s"//TRIM(ss)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_VSBRZPHI: Vacuum field in rzphi grid by "//
     $     TRIM(ss)//"th resonant field"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a2,8(a16))')"l","r","z",
     $     "real(vsb_r)","imag(vsb_r)","real(vsb_z)","imag(vsb_z)",
     $     "real(vsb_phi)","imag(vsb_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,8(es16.8))')
     $           vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $           REAL(vbr(i,j)),AIMAG(vbr(i,j)),
     $           REAL(vbz(i,j)),AIMAG(vbz(i,j)),
     $           REAL(vbp(i,j)),AIMAG(vbp(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_vsbrzphi
c-----------------------------------------------------------------------
c     subprogram 11. ipout_xbrzphifun
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,dphi,t11,t12,t21,t22,t33
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing x and b rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a29)')
     $        "volume = ",iindex,"% xi and b rzphi computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(3,3)=twopi*rs(istep,itheta)
            t11(itheta)=cos(eta)*v(1,1)-sin(eta)*v(1,2)
            t12(itheta)=cos(eta)*v(2,1)-sin(eta)*v(2,2)
            t21(itheta)=sin(eta)*v(1,1)+cos(eta)*v(1,2)
            t22(itheta)=sin(eta)*v(2,1)+cos(eta)*v(2,2)
            t33(itheta)=-1.0/v(3,3)
         ENDDO

         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         xrr_fun(istep,:)=(t11*xwp_fun+t12*xwt_fun)/jacs
         brr_fun(istep,:)=(t11*bwp_fun+t12*bwt_fun)/jacs
         xrz_fun(istep,:)=(t21*xwp_fun+t22*xwt_fun)/jacs
         brz_fun(istep,:)=(t21*bwp_fun+t22*bwt_fun)/jacs
         xrp_fun(istep,:)=t33*xvz_fun
         brp_fun(istep,:)=t33*bvz_fun

         xrr_fun(istep,:)=xrr_fun(istep,:)*EXP(ifac*nn*dphi)
         brr_fun(istep,:)=brr_fun(istep,:)*EXP(ifac*nn*dphi)
         xrz_fun(istep,:)=xrz_fun(istep,:)*EXP(ifac*nn*dphi)
         brz_fun(istep,:)=brz_fun(istep,:)*EXP(ifac*nn*dphi)
         xrp_fun(istep,:)=xrp_fun(istep,:)*EXP(ifac*nn*dphi)
         brp_fun(istep,:)=brp_fun(istep,:)*EXP(ifac*nn*dphi)         

      ENDDO

      IF(helicity>0)THEN
         xrr_fun=CONJG(xrr_fun)
         xrz_fun=CONJG(xrz_fun)
         xrp_fun=CONJG(xrp_fun)
         brr_fun=CONJG(brr_fun)
         brz_fun=CONJG(brz_fun)
         brp_fun=CONJG(brp_fun)
      ENDIF

      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbrzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBRZPHI_FUN: "//
     $     "Rzphi components of displacement and field in functions"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $     "real(xp)","imag(xp)","real(br)","imag(br)",
     $     "real(bz)","imag(bz)","real(bp)","imag(bp)"

      DO istep=1,mstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(1x,es16.8))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(xrr_fun(istep,itheta)),
     $           AIMAG(xrr_fun(istep,itheta)),
     $           REAL(xrz_fun(istep,itheta)),
     $           AIMAG(xrz_fun(istep,itheta)),
     $           REAL(xrp_fun(istep,itheta)),
     $           AIMAG(xrp_fun(istep,itheta)),
     $           REAL(brr_fun(istep,itheta)),
     $           AIMAG(brr_fun(istep,itheta)),
     $           REAL(brz_fun(istep,itheta)),
     $           AIMAG(brz_fun(istep,itheta)),
     $           REAL(brp_fun(istep,itheta)),
     $           AIMAG(brp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_xbrzphifun
c-----------------------------------------------------------------------
c     subprogram 12. ipout_arzphifun
c-----------------------------------------------------------------------
      SUBROUTINE ipout_arzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xss_fun,
     $     ear,eat,eap,arr,art,arp
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: ear_fun,eaz_fun,
     $     eap_fun,arr_fun,arz_fun,arp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing vector potential rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0))
     $        WRITE(*,'(1x,a9,i3,a37)')
     $        "volume = ",iindex,"% vector potential rzphi computations"
         CALL ipeq_sol(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xss_fun,mthsurf,xss_mn)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(istep,itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(2,1)=-rzphi%fx(2)*twopi**2*rs(istep,itheta)*rfac/jac
            w(2,2)=rzphi%fx(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(3,1)=rzphi%fx(3)*2*rfac
            w(3,2)=rzphi%fy(3)/(twopi*rfac)
            w(3,3)=1.0/(twopi*rs(istep,itheta))
            ear(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,1)-w(3,1))
            eat(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,2)-w(3,2))
            eap(itheta)=-chi1*psifac(istep)*w(3,3)
            ear_fun(istep,itheta)=ear(itheta)*cos(eta)-
     $           eat(itheta)*sin(eta)
            eaz_fun(istep,itheta)=ear(itheta)*sin(eta)+
     $           eat(itheta)*cos(eta)
            eap_fun(istep,itheta)=-eap(itheta)            
            arr(itheta)=xss_fun(itheta)*w(1,1)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,1)+xsp_fun(itheta)*w(3,1))
            art(itheta)=xss_fun(itheta)*w(1,2)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,2)+xsp_fun(itheta)*w(3,2))
            arp(itheta)=chi1*xsp_fun(itheta)*w(3,3)
            arr_fun(istep,itheta)=arr(itheta)*cos(eta)-
     $           art(itheta)*sin(eta)
            arz_fun(istep,itheta)=arr(itheta)*sin(eta)+
     $           art(itheta)*cos(eta)
            arp_fun(istep,itheta)=-arp(itheta)
         ENDDO
         ear_fun(istep,:)=ear_fun(istep,:)*EXP(ifac*nn*dphi)
         eaz_fun(istep,:)=eaz_fun(istep,:)*EXP(ifac*nn*dphi)
         eap_fun(istep,:)=eap_fun(istep,:)*EXP(ifac*nn*dphi)
         arr_fun(istep,:)=arr_fun(istep,:)*EXP(ifac*nn*dphi)
         arz_fun(istep,:)=arz_fun(istep,:)*EXP(ifac*nn*dphi)
         arp_fun(istep,:)=arp_fun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      IF(helicity>0)THEN
         ear_fun=CONJG(ear_fun)
         eaz_fun=CONJG(eaz_fun)
         eap_fun=CONJG(eap_fun)
         arr_fun=CONJG(arr_fun)
         arz_fun=CONJG(arz_fun)
         arp_fun=CONJG(arp_fun)
      ENDIF

      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_arzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_ARZPHI_FUN: "//
     $     "Rzphi components of vector potential in functions"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(ear)","imag(ear)","real(eaz)","imag(eaz)",
     $     "real(eap)","imag(eap)","real(arr)","imag(arr)",
     $     "real(arz)","imag(arz)","real(arp)","imag(arp)"

      DO istep=1,mstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(1x,es16.8))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(ear_fun(istep,itheta)),
     $           AIMAG(ear_fun(istep,itheta)),
     $           REAL(eaz_fun(istep,itheta)),
     $           AIMAG(eaz_fun(istep,itheta)),
     $           REAL(eap_fun(istep,itheta)),
     $           AIMAG(eap_fun(istep,itheta)),
     $           REAL(arr_fun(istep,itheta)),
     $           AIMAG(arr_fun(istep,itheta)),
     $           REAL(arz_fun(istep,itheta)),
     $           AIMAG(arz_fun(istep,itheta)),
     $           REAL(arp_fun(istep,itheta)),
     $           AIMAG(arp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_arzphifun
   
      END MODULE ipout_mod
