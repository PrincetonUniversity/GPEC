c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     diagnose various features.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipdiag_mod
c      1. ipdiag_eigen
c      2. ipdiag_magpot
c      3. ipdiag_arbsurf
c      4. ipdiag_angles
c      5. ipdiag_surfmode
c      6. ipdiag_singcurs
c      7. ipdiag_xbcontra
c      8. ipdiag_xbnobo
c      9. ipdiag_xbst
c     10. ipdiag_pmodbrz
c     11. ipdiag_rzphibx
c     12. ipdiag_rzpgrid
c     13. ipdiag_rzpdiv
c     14. ipdiag_radvar
c-----------------------------------------------------------------------
c     subprogram 0. ipdiag_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipdiag_mod
      USE ipresp_mod
      USE ipvacuum_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipdiag_eigen.
c     diagnose eigenvectors and eigenenergies.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_eigen
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j
      REAL(r8), DIMENSION(mpert) :: xrmn,ximn
      COMPLEX(r8), DIMENSION(mpert) :: xinmn
      REAL(r8), DIMENSION(mpert,mpert) :: potengy
c-----------------------------------------------------------------------
c     construct 2d eigenvector sets in fourier space.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_eigen.out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_EIGEN: "//
     $     "Diagnose DCON eigenvalues for eigenvectors"
      WRITE(out_unit,'(2(1x,a3),1x,a16)')"m","m","pot energy"

      DO i=1,mpert
         DO j=1,mpert
            xrmn=0
            ximn=0        
            xrmn(i)=1.0
            ximn(j)=1.0
            xinmn=xrmn+ifac*ximn

            potengy(i,j)=SUM(CONJG(xinmn)*MATMUL(wt,xinmn))
            WRITE(out_unit,'(2(1x,I3),1x,es16.8)')
     $           mfac(i),mfac(j),potengy(i,j)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_eigen
c-----------------------------------------------------------------------
c     subprogram 2. ipdiag_magpot.
c     diagnose magnetic potential errors.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_magpot
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j

      CALL ascii_open(out_unit,"ipdiag_magpot_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_MAGPOT: Magnetic potential errors"
      WRITE(out_unit,'(1x,a8,1x,I4)')"mpert=",mpert
      WRITE(out_unit,'(1x,a4,2(1x,a16))')"mode","chperr1","chperr2"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,2(1x,es16.8))')i,chperr(1,i),chperr(2,i)
      ENDDO
      WRITE(out_unit,'(1x,a4,4(1x,a16))')
     $     "mode","chpsqr1","chpsqr2","chpsqr3","chpsqr4"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,4(1x,es16.8))')
     $        i,chpsqr(1,i),chpsqr(2,i),chpsqr(3,i),chpsqr(4,i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_magpot
c-----------------------------------------------------------------------
c     subprogram 3. ipdiag_arbsurf.
c     diagnose surface inductance.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_arbsurf(majr,minr)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      REAL(r8), INTENT(IN) :: majr,minr
      INTEGER :: i,vmpert
      CHARACTER(4) :: smajr,sminr
c-----------------------------------------------------------------------
c     construct 2d eigenvector sets in fourier space.
c-----------------------------------------------------------------------
      CALL ipvacuum_arbsurf(majr,minr)
      vmpert=SIZE(vsurf_indev)
      
      WRITE(UNIT=smajr, FMT='(I4)')INT(100*majr)
      WRITE(UNIT=sminr, FMT='(I4)')INT(100*minr)
      CALL ascii_open(out_unit,"ipdiag_arbsurf_R"//smajr//"_r"//sminr//
     $     ".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_VACUUM: "//
     $     "Diagnose surface inductances for arbitrary shape"
      WRITE(out_unit,'(1x,a4,1x,a16)')"mode","surf_ind"
      DO i=1,vmpert
         WRITE(out_unit,'(1x,I4,1x,es16.8)')i,vsurf_indev(i)
      ENDDO      
      CALL ascii_close(out_unit)
      DEALLOCATE(vsurf_indmats,vsurf_indev)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_arbsurf
c-----------------------------------------------------------------------
c     subprogram 4. ipdiag_angles.
c     diagnose and visualize magnetic angles.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_angles
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: inum,istep,itheta,iqty,angnum,rstep

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,etas,dphi
      REAL(r8), DIMENSION(0:mthsurf,0:4) :: thetas,thetais
      
      REAL(r8), DIMENSION(:), POINTER :: psi,thetai,angles
      REAL(r8), DIMENSION(:,:), POINTER :: plerror,haerror
      REAL(r8), DIMENSION(:,:,:), POINTER :: rs,zs,omega
      
      TYPE(spline_type) :: spl,sdphi    
c-----------------------------------------------------------------------
c     write data for indicating used angles.
c-----------------------------------------------------------------------
      angnum=15
      rstep=400
      WRITE(*,*)"Diagnosing magnetic angles"
      ALLOCATE(psi(1:rstep))
      ALLOCATE(thetai(0:angnum-1),angles(0:angnum-1))
      ALLOCATE(plerror(1:rstep,0:mthsurf),haerror(1:rstep,0:mthsurf))
      ALLOCATE(rs(1:rstep,0:angnum-1,5),zs(1:rstep,0:angnum-1,5))
      ALLOCATE(omega(1:rstep,0:mthsurf,0:4))
      angles=(/(inum,inum=0,angnum-1)/)/REAL(angnum,r8)
      DO istep=1,rstep
         psi(istep)=(istep/REAL(rstep,r8))*psilim
         
         CALL spline_alloc(spl,mthsurf,5)
         CALL spline_alloc(sdphi,mthsurf,1)
         spl%xs=theta
         sdphi%xs=theta
         CALL spline_eval(sq,psi(istep),0)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psi(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            etas(itheta)=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(etas(itheta))
            z(itheta)=zo+rfac*SIN(etas(itheta))
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            dphi(itheta)=rzphi%f(3)
            bpfac=psio*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            spl%fs(itheta,1)=fac*twopi*bpfac/rfac
            spl%fs(itheta,2)=fac                        
            spl%fs(itheta,3)=fac/r(itheta)**2
            spl%fs(itheta,4)=fac*bpfac
            spl%fs(itheta,5)=fac*bfac**2
            sdphi%fs(itheta,1)=dphi(itheta)
         ENDDO
         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)
         CALL spline_fit(sdphi,"periodic")
         DO iqty=1,5
            thetas(:,iqty-1)=spl%fsi(:,iqty)/spl%fsi(mthsurf,iqty)
            DO itheta=0,mthsurf
               thetais(itheta,iqty-1)=issect(mthsurf,theta(:),
     $              thetas(:,iqty-1),theta(itheta))
               CALL spline_eval(sdphi,thetais(itheta,iqty-1),0)
               omega(istep,itheta,iqty-1)=-sdphi%f(1)+twopi*sq%f(4)*
     $              (theta(itheta)-thetais(itheta,iqty-1))
            ENDDO
         ENDDO
         CALL spline_dealloc(spl)
         CALL spline_dealloc(sdphi)
c-----------------------------------------------------------------------
c     store errors.
c-----------------------------------------------------------------------
         plerror(istep,:)=(etas(:)/twopi-thetas(:,0))
         haerror(istep,:)=thetas(:,1)-thetais(:,1)
c-----------------------------------------------------------------------
c     store visualizing data.
c-----------------------------------------------------------------------         
         DO iqty=1,5
            DO inum=0,angnum-1
               thetai(inum)=issect(mthsurf,theta(:),thetas(:,iqty-1),
     $              angles(inum))
            ENDDO
            DO inum=0,angnum-1
               CALL bicube_eval(rzphi,psi(istep),thetai(inum),0)
               rfac=SQRT(rzphi%f(1)) 
               eta=twopi*(thetai(inum)+rzphi%f(2))
               rs(istep,inum,iqty)=ro+rfac*COS(eta)
               zs(istep,inum,iqty)=zo+rfac*SIN(eta)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------      
      CALL ascii_open(out_unit,"ipdiag_angles.out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_ANGLES: "//
     $     "Diagnose and visualize magnetic angles"
      WRITE(out_unit,'(1x,a8,1x,I6)')"angnum=",angnum
      WRITE(out_unit,'(1x,a8,1x,I6)')"rstep=",rstep
      WRITE(out_unit,'(1x,a8,1x,I6)')"mthsurf=",mthsurf
      WRITE(out_unit,'(3(1x,a16))')"psi","plerror","haerror"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(3(1x,es16.8))')
     $           psi(istep),plerror(istep,itheta),haerror(istep,itheta)
         ENDDO
      ENDDO
      
      DO iqty=1,5
         WRITE(out_unit,'(1x,a8,1x,I2)')"iqty=",iqty
         WRITE(out_unit,'(3(1x,a16))')"r","z","angles"
         DO inum=0,angnum-1
            DO istep=1,rstep
               WRITE(out_unit,'(3(1x,es16.8))')
     $              rs(istep,inum,iqty),zs(istep,inum,iqty),angles(inum)
            ENDDO
         ENDDO
      ENDDO

      WRITE(out_unit,'(7(1x,a16))')"psi","theta","polar","hamada",
     $     "pest","equalarc","boozer"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(1x,es16.8))')
     $           psi(istep),theta(itheta),omega(istep,itheta,0),
     $           omega(istep,itheta,1),omega(istep,itheta,2),
     $           omega(istep,itheta,3),omega(istep,itheta,4)
         ENDDO      
      ENDDO
      CALL ascii_close(out_unit)      

      DEALLOCATE(psi,thetai,angles,plerror,haerror,rs,zs,omega)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_angles
c-----------------------------------------------------------------------
c     subprogram 5. ipdiag_surfmode.
c     response to fourier modes for the control surface.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_surfmode(lowmode,highmode,
     $     rin,bpin,bin,rcin,tin,jin)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: lowmode,highmode,rin,bpin,bin,rcin,tin,jin

      INTEGER :: i,j,mnum

      COMPLEX(r8), DIMENSION(:,:), POINTER :: binmn,boutmn,finmn,
     $     binfun,boutfun

      mnum=highmode-lowmode+1
      ALLOCATE(binmn(mnum,mpert),boutmn(mnum,mpert),finmn(mnum,mpert),
     $     binfun(mnum,0:mthsurf),boutfun(mnum,0:mthsurf))
     
      DO i=lowmode,highmode
         j=i-lowmode+1
         binmn(j,:)=0
         binmn(j,i-mlow+1)=1
         finmn(j,:)=binmn(j,:)
         CALL iscdftb(mfac,mpert,binfun(j,:),mthsurf,finmn(j,:))         
         CALL ipeq_fcoords(psilim,finmn(j,:),mfac,mpert,
     $        rin,bpin,bin,rcin,tin,jin)
         CALL ipeq_weight(psilim,finmn(j,:),mfac,mpert,1)  
         IF (fixed_boundary_flag) THEN
            boutmn(j,:)=finmn(j,:)
         ELSE
            boutmn(j,:)=MATMUL(permeabmats(resp_index,:,:),finmn(j,:))
         ENDIF       
         CALL ipeq_weight(psilim,finmn(j,:),mfac,mpert,0)  
         CALL ipeq_bcoords(psilim,boutmn(j,:),mfac,mpert,
     $        rin,bpin,bin,rcin,tin,jin)
         CALL iscdftb(mfac,mpert,boutfun(j,:),mthsurf,boutmn(j,:))
      ENDDO
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_surfmode.out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_SURFMODE: Plasma response for the "//
     $     "fourier modes on the control surface"
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(1x,a4,4(1x,a16))')"m",
     $     "rebin","imbin","rebout","imbout"
      DO j=1,mnum
         WRITE(out_unit,'(1x,a6,i3)')"mode=",j-1+lowmode
         DO i=1,mpert
            WRITE(out_unit,'(1x,I4,4(1x,es16.8))')mfac(i),
     $           REAL(binmn(j,i)),AIMAG(binmn(j,i)),
     $           REAL(boutmn(j,i)),AIMAG(boutmn(j,i))
         ENDDO
      ENDDO
      WRITE(out_unit,*)"FUNCTIONS"
      WRITE(out_unit,'(5(1x,a16))')"theta",
     $     "rebin","imbin","rebout","imbout"
      DO j=1,mnum
         WRITE(out_unit,'(1x,a6,i3)')"mode=",j-1+lowmode
         DO i=0,mthsurf
            WRITE(out_unit,'(5(1x,es16.8))')theta(i),
     $           REAL(binfun(j,i)),AIMAG(binfun(j,i)),
     $           REAL(boutfun(j,i)),AIMAG(boutfun(j,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(binmn,boutmn,finmn,binfun,boutfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_surfmode
c-----------------------------------------------------------------------
c     subprogram 6. ipdiag_singcurs.
c     diagnose asymtotic values of singular currents.
c     __________________________________________________________________
c     egnjm      : eigenmode number without edge_flag
c     xwpimn     : edge deformation input when edge_flag
c     rsing      : number of rational surfaces
c     resol      : resolution by number of grid points
c     smallwidth : closest point to approach
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_singcurs(egnum,xwpimn,rsing,resol,smallwidth)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,resol,rsing
      REAL(r8), INTENT(IN) :: smallwidth
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: ising,i,itheta,resnum
      REAL(r8) :: respsi,rrespsi,lrespsi,bigwidth,astep,lpsi,rpsi,
     $     sqrpsi,correc,lq1,rq1
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      INTEGER, DIMENSION(rsing) :: logsteps
      REAL(r8), DIMENSION(rsing) :: j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: lcormn,rcormn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun

      REAL(r8), DIMENSION(rsing,resol*10) :: spot
      COMPLEX(r8), DIMENSION(rsing,resol*10) :: deltas,delcurs,corcurs,
     $     singcurs

      spot=0
      deltas=0
      delcurs=0
      corcurs=0
      singcurs=0
      CALL ipeq_alloc
      CALL idcon_build(egnum,xwpimn)
c-----------------------------------------------------------------------
c     compute singular currents with logarithmic approach.
c-----------------------------------------------------------------------
      WRITE(*,*)"Diagnosing asymtotic values of singular currents"
      DO ising=1,rsing
         resnum=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac

         IF (ising == 1) THEN
            lrespsi=0.001
         ELSE
            lrespsi=singtype(ising-1)%psifac
         ENDIF
         IF (ising == msing) THEN
            rrespsi=psilim
         ELSE
            rrespsi=singtype(ising+1)%psifac
         ENDIF

         bigwidth=MIN(respsi-lrespsi,rrespsi-respsi)/2.0
         CALL spline_eval(sq,lrespsi,1)
         lq1=ABS(sq%f1(4))
         CALL spline_eval(sq,rrespsi,1)
         rq1=ABS(sq%f1(4))
         bigwidth=bigwidth*nn*MIN(lq1,rq1)
         logsteps(ising)=1+INT((resol-1)*log10(bigwidth/smallwidth))
c-----------------------------------------------------------------------
c     calculate surface averaged characteristic current.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,respsi,0)
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
     $              /(twopi*r(itheta))**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
         ENDDO
         j_c(ising)=1.0/issurfint(jcfun,mthsurf,respsi,0,0)*
     $        chi1**2*sq%f(4)/mu0
c-----------------------------------------------------------------------
c     main loop for surface current on the rational surfaces.
c-----------------------------------------------------------------------
         DO i=1,logsteps(ising)
            astep=smallwidth*10.0**((i-1.)/(resol-1.))
            lpsi=respsi-astep/(nn*ABS(singtype(ising)%q1))
            rpsi=respsi+astep/(nn*ABS(singtype(ising)%q1))
c-----------------------------------------------------------------------
c     calculate delta and correction term at lpsi.
c-----------------------------------------------------------------------
            CALL ipeq_sol(lpsi)
            lbwp1mn=bwp1_mn(resnum)
            CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
            CALL spline_eval(sq,lpsi,0)
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
c-----------------------------------------------------------------------
c     calculate delta and correction term at rpsi.
c-----------------------------------------------------------------------
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
            spot(ising,i)=astep/(nn*ABS(singtype(ising)%q1))

            deltas(ising,i)=rbwp1mn-lbwp1mn
            delcurs(ising,i)=j_c(ising)*ifac/(twopi*mfac(resnum))*
     $           deltas(ising,i)
            corcurs(ising,i)=-j_c(ising)*
     $           (rcormn(resnum)-lcormn(resnum))
            delcurs(ising,i)=-delcurs(ising,i)/(twopi*ifac*nn)
            corcurs(ising,i)=-corcurs(ising,i)/(twopi*ifac*nn)
            singcurs(ising,i)=delcurs(ising,i)-corcurs(ising,i)
         ENDDO

         WRITE(*,*)"Finished the analysis for the q =",singtype(ising)%q
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write the results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_deltas_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_DELTAS: "//
     $     "Asymtotic analysis for deltas at rational surfaces"
      WRITE(out_unit,'(1x,a8,1x,I6)')"rsing=",rsing
      WRITE(out_unit,'(1x,a8,1x,I6)')"resol=",resol
      WRITE(out_unit,'(3(1x,a16))')"distance","re(delta)","im(delta)"
      DO ising=1,rsing
         WRITE(out_unit,'(1x,a2,1x,f6.3)')"q=",singtype(ising)%q
         WRITE(out_unit,'(1x,a12,1x,I6)')"logsteps=",logsteps(ising)
         DO i=1,logsteps(ising)
            WRITE(out_unit,'(3(1x,es16.8))')spot(ising,i),
     $           REAL(deltas(ising,i)),AIMAG(deltas(ising,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipdiag_singcurs_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_SINGCURS: "//
     $     "asymtotic analysis for singular currents at "//
     $     "rational surfaces"
      WRITE(out_unit,'(1x,a8,1x,I6)')"rsing=",rsing
      WRITE(out_unit,'(1x,a8,1x,I6)')"resol=",resol
      WRITE(out_unit,'(7(1x,a16))')"distance","re(delcur)","im(delcur)",
     $     "re(corcur)","im(corcur)","re(singcur)","im(singcur)"
      DO ising=1,rsing
         WRITE(out_unit,'(1x,a2,1x,f6.3)')"q=",singtype(ising)%q
         WRITE(out_unit,'(1x,a12,1x,I6)')"logsteps=",logsteps(ising)
         DO i=1,logsteps(ising)
            WRITE(out_unit,'(7(1x,es16.8))')spot(ising,i),
     $           REAL(delcurs(ising,i)),AIMAG(delcurs(ising,i)),
     $           REAL(corcurs(ising,i)),AIMAG(corcurs(ising,i)),
     $           REAL(singcurs(ising,i)),AIMAG(singcurs(ising,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_singcurs
c-----------------------------------------------------------------------
c     subprogram 7. ipdiag_xbcontra.
c     diagnose various components of xi and b.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbcontra(egnum,xwpimn,rin,bpin,bin,rcin,tin)    
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rin,bpin,bin,rcin,tin
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert
      COMPLEX(r8), DIMENSION(mpert) :: xwd_mn,bwd_mn

      TYPE(cspline_type) :: u3         
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)

      WRITE(*,*)"Computing contravariant components in detail"
      CALL ascii_open(out_unit,"ipdiag_xbcontra_n"//
     $        TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_XBCONTRA: "//
     $     "Contravariant components of displacement and field"
      WRITE(out_unit,'(1x,a8,1x,I6)')"mstep=",mstep
      WRITE(out_unit,'(1x,a8,1x,I6)')"mlow=",mlow
      WRITE(out_unit,'(1x,a8,1x,I6)')"mhigh=",mhigh
      WRITE(out_unit,'(2(1x,a16))')"psifac","qfac"
      CALL ipeq_alloc
      DO istep=0,mstep
         WRITE(out_unit,'(2(1x,es16.8))')psifac(istep),qfac(istep)
      ENDDO
      WRITE(out_unit,'(22(1x,a16))')
     $     "xwp(real)","xwp(imag)","xwt(real)","xwt(imag)",
     $     "xwz(real)","xwz(imag)","bwp(real)","bwp(imag)",
     $     "bwt(real)","bwt(imag)","bwz(real)","bwz(imag)",
     $     "xvs(real)","xvs(imag)","xwd(real)","xwd(imag)",
     $     "xwp1(real)","xwp1(imag)","bwd(real)","bwd(imag)",
     $     "bwp1(real)","bwp1(imag)"

      CALL cspline_alloc(u3,mstep,mpert)
      u3%xs=psifac
      DO istep=0,mstep
         CALL ipeq_sol(psifac(istep))
         u3%fs(istep,:)=bwp_mn
      ENDDO
      CALL cspline_fit(u3,"extrap")

      DO istep=0,mstep
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL cspline_eval(u1,psifac(istep),1)
         CALL cspline_eval(u3,psifac(istep),1)
         xwd_mn(:)=u1%f1(:)
         bwd_mn(:)=u3%f1(:)
         IF (jac_type /= jac_out) THEN 
            CALL ipeq_bcoords(psifac(istep),xwp_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),xwt_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),xwz_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),bwp_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),bwt_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),bwz_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),xss_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),xwd_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),xsp1_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),bwd_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
            CALL ipeq_bcoords(psifac(istep),bwp1_mn,mfac,mpert,
     $           rin,bpin,bin,rcin,tin,0)
         ENDIF
         DO ipert=1,mpert
            WRITE(out_unit,'(22(1x,es16.8))')
     $           REAL(xwp_mn(ipert)),AIMAG(xwp_mn(ipert)),
     $           REAL(xwt_mn(ipert)),AIMAG(xwt_mn(ipert)),
     $           REAL(xwz_mn(ipert)),AIMAG(xwz_mn(ipert)),
     $           REAL(bwp_mn(ipert)),AIMAG(bwp_mn(ipert)),
     $           REAL(bwt_mn(ipert)),AIMAG(bwt_mn(ipert)),
     $           REAL(bwz_mn(ipert)),AIMAG(bwz_mn(ipert)),
     $           REAL(xss_mn(ipert)),AIMAG(xss_mn(ipert)),
     $           REAL(xwd_mn(ipert)),AIMAG(xwd_mn(ipert)),
     $           REAL(xsp1_mn(ipert)),AIMAG(xsp1_mn(ipert)),
     $           REAL(bwd_mn(ipert)),AIMAG(bwd_mn(ipert)),
     $           REAL(bwp1_mn(ipert)),AIMAG(bwp1_mn(ipert))
         ENDDO
      ENDDO

      CALL ipeq_dealloc
      CALL cspline_dealloc(u3)
      CALL ascii_close(out_unit)

      RETURN
      END SUBROUTINE ipdiag_xbcontra
c-----------------------------------------------------------------------
c     subprogram 8. ipdiag_xbnobo.
c     write normal perturbed quantities on the boundary.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbnobo(egnum,xwpimn,d3_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) ::xwpimn
      LOGICAL, INTENT(IN) :: d3_flag

      INTEGER :: i,ithnum,mthnum
      INTEGER(r8) :: nnn,mmtheta,np,nt
      REAL(r8) :: ddist
      CHARACTER(72) :: filein,file1

      REAL(r8), DIMENSION(:), POINTER :: mthang,rs,zs,delpsi,
     $     norvec,nozvec
      COMPLEX(r8), DIMENSION(:), POINTER :: xno_fun,bno_fun,
     $     xnorvc,xnozvc,bnorvc,bnozvc


      WRITE(*,*)"Computing normal quantities on the boundary"
      mthnum=1000
      ALLOCATE(rs(0:mthnum),zs(0:mthnum),
     $     mthang(0:mthnum),delpsi(0:mthnum),
     $     xnorvc(0:mthnum),xnozvc(0:mthnum),
     $     bnorvc(0:mthnum),bnozvc(0:mthnum),
     $     norvec(0:mthnum),nozvec(0:mthnum),
     $     xno_fun(0:mthnum),bno_fun(0:mthnum))
      mthang=(/(i,i=0,mthnum)/)/REAL(mthnum,r8)
      
      DO ithnum=0,mthnum
         CALL bicube_eval(rzphi,psilim,mthang(ithnum),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(mthang(ithnum)+rzphi%f(2))
         rs(ithnum)=ro+rfac*COS(eta)
         zs(ithnum)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(ithnum)/jac
         w(1,2)=-rzphi%fy(1)*pi*rs(ithnum)/(rfac*jac)
         delpsi(ithnum)=SQRT(w(1,1)**2+w(1,2)**2)
         norvec(ithnum)=(cos(eta)*w(1,1)-sin(eta)*w(1,2))/delpsi(ithnum)
         nozvec(ithnum)=(sin(eta)*w(1,1)+cos(eta)*w(1,2))/delpsi(ithnum)
      ENDDO
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      CALL ipeq_alloc
      CALL ipeq_sol(psilim)
      CALL ipeq_contra(psilim)
      CALL ipeq_normal(psilim)
      CALL ipeq_bcoords(psilim,xno_mn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL ipeq_bcoords(psilim,bno_mn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL iscdftb(mfac,mpert,xno_fun,mthnum,xno_mn)
      CALL iscdftb(mfac,mpert,bno_fun,mthnum,bno_mn)
      CALL ipeq_dealloc
      
      xnorvc=xno_fun*norvec
      xnozvc=xno_fun*nozvec
      bnorvc=bno_fun*norvec
      bnozvc=bno_fun*nozvec
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_xnobo_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_XNOBO: "//
     $     "Perturbed normal displacement vectors on the boundary"      
      WRITE(out_unit,'(6(1x,a16))')"r","z",
     $     "real(xnor)","real(xnoz)","imag(xnor)","imag(xnoz)"
      DO ithnum=0,mthnum
         WRITE(out_unit,'(6(1x,es16.8))')rs(ithnum),zs(ithnum),
     $        REAL(xnorvc(ithnum)),REAL(xnozvc(ithnum)),
     $        AIMAG(xnorvc(ithnum)),AIMAG(xnozvc(ithnum))
      ENDDO
      CALL ascii_open(out_unit,"ipdiag_bnobo_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_BNOBO: "//
     $     "perturbed normal field vectors on the boundary"      
      WRITE(out_unit,'(6(1x,a16))')"r","z",
     $     "real(bnor)","real(bnoz)","imag(bnor)","imag(bnoz)"
      DO ithnum=0,mthnum
         WRITE(out_unit,'(6(1x,es16.8))')rs(ithnum),zs(ithnum),
     $        REAL(bnorvc(ithnum)),REAL(bnozvc(ithnum)),
     $        AIMAG(bnorvc(ithnum)),AIMAG(bnozvc(ithnum))
      ENDDO
      DEALLOCATE(mthang,delpsi,
     $     norvec,nozvec,xnorvc,xnozvc,bnorvc,bnozvc)
c-----------------------------------------------------------------------
c     write data for 3d surface plot.
c-----------------------------------------------------------------------
      IF (d3_flag) THEN
         CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
         WRITE(out_unit,'(1p,5es16.8)')rs
         WRITE(out_unit,'(1p,5es16.8)')zs
         WRITE(out_unit,'(1p,5es16.8)')REAL(xno_fun)
         WRITE(out_unit,'(1p,5es16.8)')AIMAG(xno_fun)      
         CALL ascii_close(out_unit)
         
         filein="iptemp.txt"
         file1="ipidl_3dsurf_xnobo_n"//TRIM(sn)//".out"

         ddist=0.2
         mmtheta=mthnum
         nnn=nn
         np=144
         nt=144
         CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1)

         CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
         WRITE(out_unit,'(1p,5es16.8)')rs
         WRITE(out_unit,'(1p,5es16.8)')zs
         WRITE(out_unit,'(1p,5es16.8)')REAL(bno_fun)
         WRITE(out_unit,'(1p,5es16.8)')AIMAG(bno_fun)      
         CALL ascii_close(out_unit)
         
         filein="iptemp.txt"
         file1="ipidl_3dsurf_bnobo_n"//TRIM(sn)//".out"

         ddist=0.2
         mmtheta=mthnum
         nnn=nn
         np=144
         nt=144
         CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1) 
      ENDIF
      DEALLOCATE(rs,zs,xno_fun,bno_fun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_xbnobo
c-----------------------------------------------------------------------
c     subprogram 9. ipdiag_xbst.
c     diagnose strength of x and b.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbst(egnum,xwpimn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,itheta,ipert

      REAL(r8), DIMENSION(0:mthsurf) :: jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,bwt_fun,bwz_fun,
     $     bvp_fun,bvt_fun,bvz_fun,brr_fun,brz_fun,brp_fun,
     $     xwp_fun,xwt_fun,xwz_fun,
     $     xvp_fun,xvt_fun,xvz_fun,xrr_fun,xrz_fun,xrp_fun

      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: x1,x2,b1,b2
      COMPLEX(r8), DIMENSION(mstep,mpert) :: x1mns,x2mns,b1mns,b2mns
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing x and b field strength"

      CALL idcon_build(egnum,xwpimn)
      
      CALL ipeq_alloc
      DO istep=1,mstep
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces.
c-----------------------------------------------------------------------
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))
         CALL ipeq_rzphi(psifac(istep))

         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
         CALL iscdftb(mfac,mpert,bwz_fun,mthsurf,bwz_mn)
         CALL iscdftb(mfac,mpert,bvp_fun,mthsurf,bvp_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         CALL iscdftb(mfac,mpert,brr_fun,mthsurf,brr_mn)
         CALL iscdftb(mfac,mpert,brz_fun,mthsurf,brz_mn)
         CALL iscdftb(mfac,mpert,brp_fun,mthsurf,brp_mn)

         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
         CALL iscdftb(mfac,mpert,xwz_fun,mthsurf,xwz_mn)
         CALL iscdftb(mfac,mpert,xvp_fun,mthsurf,xvp_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,xrr_fun,mthsurf,xrr_mn)
         CALL iscdftb(mfac,mpert,xrz_fun,mthsurf,xrz_mn)
         CALL iscdftb(mfac,mpert,xrp_fun,mthsurf,xrp_mn)

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),0)
            jacs(itheta)=rzphi%f(4)
         ENDDO
         x1(istep,:)=(CONJG(xwp_fun)*xvp_fun+CONJG(xwt_fun)*xvt_fun+
     $        CONJG(xwz_fun)*xvz_fun)/jacs
         x2(istep,:)=CONJG(xrr_fun)*xrr_fun+CONJG(xrz_fun)*xrz_fun+
     $        CONJG(xrp_fun)*xrp_fun
         b1(istep,:)=(CONJG(bwp_fun)*bvp_fun+CONJG(bwt_fun)*bvt_fun+
     $        CONJG(bwz_fun)*bvz_fun)/jacs
         b2(istep,:)=CONJG(brr_fun)*brr_fun+CONJG(brz_fun)*brz_fun+
     $        CONJG(brp_fun)*brp_fun
         CALL iscdftf(mfac,mpert,x1(istep,:),mthsurf,x1mns(istep,:))
         CALL iscdftf(mfac,mpert,x2(istep,:),mthsurf,x2mns(istep,:))
         CALL iscdftf(mfac,mpert,b1(istep,:),mthsurf,b1mns(istep,:))
         CALL iscdftf(mfac,mpert,b2(istep,:),mthsurf,b2mns(istep,:))
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------      
      CALL ascii_open(out_unit,"ipdiag_xbst_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_XBST: "//
     $     "Perturbed x and b strength on flux surfaces"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a16,1x,a4,8(1x,a16))')"psi","m",
     $     "real(x1)","imag(x1)","real(x2)","imag(x2)",
     $     "real(b1)","imag(b2)","real(b2)","imag(b2)"
      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es16.8,1x,I4,8(1x,es16.8))')
     $           psifac(istep),mfac(ipert),
     $           REAL(x1(istep,ipert)),
     $           AIMAG(x1(istep,ipert)),
     $           REAL(x2(istep,ipert)),
     $           AIMAG(x2(istep,ipert)),
     $           REAL(b1(istep,ipert)),
     $           AIMAG(b1(istep,ipert)),
     $           REAL(b2(istep,ipert)),
     $           AIMAG(b2(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_xbst
c-----------------------------------------------------------------------
c     subprogram 10. ipdiag_pmodbrz.
c     plot perturbed mod b at rz coordinates.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_pmodbrz(egnum,xwpimn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert,itheta

      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,
     $     llagbpar_mn,divx_mn,curv_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmt_fun,xmz_fun,xvt_fun,xvz_fun,xsp1_fun

      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs,eqfunx,eqfuny,eqfuns
      COMPLEX(r8), DIMENSION(:,:), POINTER :: eulbparmns,lagbparmns,
     $     eulbparfun,lagbparfun,llagbparmns,llagbparfun,
     $     xspfun,xmsfun,bvtfun,bvzfun,xmtfun,xmzfun,xvtfun,xvzfun,
     $     divxfun,curvfun,divxmns,curvmns,lllagbparmns,cdeltamns

      TYPE(cspline_type) :: cspl       

c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing perturbed b field"
      ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf),
     $     eulbparmns(mstep,mpert),lagbparmns(mstep,mpert),
     $     eulbparfun(mstep,0:mthsurf),lagbparfun(mstep,0:mthsurf),
     $     llagbparmns(mstep,mpert),llagbparfun(mstep,0:mthsurf),
     $     divxmns(mstep,mpert),curvmns(mstep,mpert),
     $     lllagbparmns(mstep,mpert),cdeltamns(mstep,mpert))
      ALLOCATE(xspfun(mstep,0:mthsurf),xmsfun(mstep,0:mthsurf),
     $     bvtfun(mstep,0:mthsurf),bvzfun(mstep,0:mthsurf),
     $     xmzfun(mstep,0:mthsurf),xvtfun(mstep,0:mthsurf),
     $     xmtfun(mstep,0:mthsurf),xvzfun(mstep,0:mthsurf),
     $     eqfunx(mstep,0:mthsurf),eqfuny(mstep,0:mthsurf),
     $     eqfuns(mstep,0:mthsurf),
     $     divxfun(mstep,0:mthsurf),curvfun(mstep,0:mthsurf))

      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta

      CALL idcon_build(egnum,xwpimn)
      CALL ipeq_alloc
      DO istep=1,mstep
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b variations.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         CALL iscdftb(mfac,mpert,xmt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xmz_fun,mthsurf,xmz_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)

         xspfun(istep,:)=xsp_fun
         xmsfun(istep,:)=xms_fun
         bvtfun(istep,:)=bvt_fun
         bvzfun(istep,:)=bvz_fun
         xmtfun(istep,:)=xmt_fun
         xmzfun(istep,:)=xmz_fun
         xvtfun(istep,:)=xvt_fun
         xvzfun(istep,:)=xvz_fun

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

            eqfuns(istep,itheta)=eqfun%f(1)
            eqfunx(istep,itheta)=eqfun%fx(1)
            eqfuny(istep,itheta)=eqfun%fy(1)

            eulbparfun(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(rzphi%f(4)*eqfun%f(1))
            lagbparfun(istep,itheta)=
     $           eulbparfun(istep,itheta)+
     $           xsp_fun(itheta)*eqfun%fx(1)+
     $           xms_fun(itheta)/(chi1*sq%f(4))*eqfun%fy(1)
            llagbparfun(istep,itheta)=
     $           lagbparfun(istep,itheta)+eqfun%fy(1)*
     $           (xmz_fun(itheta)/(jac*sq%f(4))-
     $           (chi1/(jac*eqfun%f(1)))**2*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta)))
            divxfun(istep,itheta)=xsp1_fun(itheta)+
     $           (jac1/jac)*xsp_fun(itheta)+
     $           cspl%f1(1)/jac-(twopi*ifac*nn)*cspl%f(2)/jac
            divxfun(istep,itheta)=-eqfun%f(1)*divxfun(istep,itheta)
            curvfun(istep,itheta)=-xsp_fun(itheta)/eqfun%f(1)*sq%f1(2)-
     $           (xsp_fun(itheta)*eqfun%fx(1)+eqfun%fy(1)*
     $           (xms_fun(itheta)/(chi1*sq%f(4))+
     $           xmz_fun(itheta)/(jac*sq%f(4))-
     $           (chi1/(jac*eqfun%f(1)))**2*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))))
         ENDDO
         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
     $        mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
     $        mthsurf,lagbpar_mn)
         CALL iscdftf(mfac,mpert,llagbparfun(istep,:),
     $        mthsurf,llagbpar_mn)
         CALL iscdftf(mfac,mpert,divxfun(istep,:),
     $        mthsurf,divx_mn)
         CALL iscdftf(mfac,mpert,curvfun(istep,:),
     $        mthsurf,curv_mn)

         eulbparmns(istep,:)=eulbpar_mn
         lagbparmns(istep,:)=lagbpar_mn
         llagbparmns(istep,:)=llagbpar_mn
         divxmns(istep,:)=divx_mn
         curvmns(istep,:)=curv_mn
         lllagbparmns(istep,:)=divx_mn+curv_mn
         cdeltamns(istep,:)=divx_mn+2*curv_mn
      ENDDO

      CALL cspline_dealloc(cspl)
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_pmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IDIAG_PMODB: "//
     $     "Components in perturbed mod b"
      WRITE(out_unit,*)     
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     

      WRITE(out_unit,'(1x,a16,1x,a4,10(1x,a16))')"psi","m",
     $     "real(lagb)","imag(lagb)","real(llagb)","imag(llagb)",
     $     "real(lllagb)","imag(lllagb)","real(divx)","imag(divx)",
     $     "real(cdelta)","imag(cdelta)"
      DO istep=1,mstep
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es16.8,1x,I4,10(1x,es16.8))')
     $           psifac(istep),mfac(ipert),
     $           REAL(lagbparmns(istep,ipert)),
     $           AIMAG(lagbparmns(istep,ipert)),
     $           REAL(llagbparmns(istep,ipert)),
     $           AIMAG(llagbparmns(istep,ipert)),
     $           REAL(lllagbparmns(istep,ipert)),
     $           AIMAG(lllagbparmns(istep,ipert)),
     $           REAL(divxmns(istep,ipert)),
     $           AIMAG(divxmns(istep,ipert)),
     $           REAL(cdeltamns(istep,ipert)),
     $           AIMAG(cdeltamns(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      CALL ascii_open(out_unit,"ipdiag_pmodbrz_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_PMODBRZ: "//
     $     "Perturbed mod b on a poloidal plane"
      WRITE(out_unit,'(1x,a8,1x,I4)')"mstep=",mstep
      WRITE(out_unit,'(1x,a8,1x,I4)')"mthsurf=",mthsurf
      WRITE(out_unit,'(17(1x,a16))')"r","z",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $     "real(xwp)","imag(xwp)","real(xwt)","imag(xwt)",
     $     "real(bvt)","imag(bvt)","real(bvz)","imag(bvz)",
     $     "eqfunx","eqfuny","eqfuns"
      DO istep=1,mstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(17(1x,es16.8))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(eulbparfun(istep,itheta)),
     $           AIMAG(eulbparfun(istep,itheta)),
     $           REAL(llagbparfun(istep,itheta)),
     $           AIMAG(llagbparfun(istep,itheta)),
     $           REAL(xspfun(istep,itheta)),
     $           AIMAG(xspfun(istep,itheta)),
     $           REAL(xmsfun(istep,itheta)),
     $           AIMAG(xmsfun(istep,itheta)),
     $           REAL(bvtfun(istep,itheta)),
     $           AIMAG(bvtfun(istep,itheta)),
     $           REAL(bvzfun(istep,itheta)),
     $           AIMAG(bvzfun(istep,itheta)),
     $           eqfunx(istep,itheta),eqfuny(istep,itheta),
     $           eqfuns(istep,itheta)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      CALL ipeq_dealloc
      DEALLOCATE(rs,zs,eulbparmns,lagbparmns,llagbparmns,cdeltamns,
     $     eulbparfun,lagbparfun,llagbparfun,eqfunx,eqfuny,eqfuns,
     $     xspfun,xmsfun,bvtfun,bvzfun,xmzfun,xvtfun,xvzfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_pmodbrz
c-----------------------------------------------------------------------
c     subprogram 11. ipdiag_rzphibx.
c     write r,z,phi,b,x on hamada coordinates.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_rzphibx(egnum,xwpimn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: i,istep,itheta,ipert,rstep

      REAL(r8) :: dist,limfac
      REAL(r8), DIMENSION(0:mthsurf) :: t11,t12,t21,t22,t33,jacs
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun

      REAL(r8), DIMENSION(:), POINTER :: psis
      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs,phs
      COMPLEX(r8), DIMENSION(:,:), POINTER :: brs,bzs,bps,xrs,xzs,xps
c-----------------------------------------------------------------------
c     adjust parameters temporarily, dist is important.
c-----------------------------------------------------------------------
      rstep=200
      dist=0.1
c-----------------------------------------------------------------------
c     allocate memory.
c-----------------------------------------------------------------------
      ALLOCATE(psis(rstep),rs(rstep,0:mthsurf),zs(rstep,0:mthsurf),
     $     phs(rstep,0:mthsurf),
     $     brs(rstep,0:mthsurf),bzs(rstep,0:mthsurf),
     $     bps(rstep,0:mthsurf),
     $     xrs(rstep,0:mthsurf),xzs(rstep,0:mthsurf),
     $     xps(rstep,0:mthsurf))
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      IF (rstep == mstep) THEN
         psis=psifac
      ELSE
         psis=(/(i,i=1,rstep+1)/)/REAL(rstep+1,r8)*psilim
      ENDIF
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing rzphibx components"
      CALL ipeq_alloc
      DO istep=1,rstep
         CALL spline_eval(sq,psis(istep),0)
         CALL ipeq_sol(psis(istep))
         CALL ipeq_contra(psis(istep))
         CALL ipeq_cova(psis(istep))
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psis(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            phs(istep,itheta)=rzphi%f(3)
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(3,3)=twopi*rs(istep,itheta)
            t11(itheta)=cos(eta)*v(1,1)-sin(eta)*v(1,2)
            t12(itheta)=cos(eta)*v(2,1)-sin(eta)*v(2,2)
            t21(itheta)=sin(eta)*v(1,1)+cos(eta)*v(1,2)
            t22(itheta)=sin(eta)*v(2,1)+cos(eta)*v(2,2)
            t33(itheta)=1.0/v(3,3)
         ENDDO
c-----------------------------------------------------------------------
c     regulations.
c-----------------------------------------------------------------------
         singfac=mfac-nn*sq%f(4)
         DO ipert=1,mpert
            limfac = 1.0
            IF (ABS(singfac(ipert))<dist) THEN
               limfac = singfac(ipert)/SIGN(dist,singfac(ipert))
            ENDIF
            xwt_mn(ipert)=xwt_mn(ipert)*limfac
            xvz_mn(ipert)=xvz_mn(ipert)*limfac
         ENDDO
c-----------------------------------------------------------------------
c     three vector components.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         xrs(istep,:)=(t11*xwp_fun+t12*xwt_fun)/jacs
         brs(istep,:)=(t11*bwp_fun+t12*bwt_fun)/jacs
         xzs(istep,:)=(t21*xwp_fun+t22*xwt_fun)/jacs
         bzs(istep,:)=(t21*bwp_fun+t22*bwt_fun)/jacs
         xps(istep,:)=t33*xvz_fun
         bps(istep,:)=t33*bvz_fun
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_rzphibx_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RZPHIBX: "//
     $     "Give rzphibx information in hamada coordinates"
      WRITE(out_unit,'(1x,a8,1x,I6)')"rstep=",rstep
      WRITE(out_unit,'(1x,a8,1x,I6)')"mthsurf=",mthsurf
      WRITE(out_unit,'(17(1x,a16))')"psi","theta","r","z","phi",
     $     "real(br)","imag(br)","real(bz)","imag(bz)",
     $     "real(bp)","imag(bp)",
     $     "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $     "real(xp)","imag(xp)"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(17(1x,es16.8))')
     $           psis(istep),theta(itheta),
     $           rs(istep,itheta),zs(istep,itheta),phs(istep,itheta),
     $           REAL(brs(istep,itheta)),AIMAG(brs(istep,itheta)),
     $           REAL(bzs(istep,itheta)),AIMAG(bzs(istep,itheta)),
     $           REAL(bps(istep,itheta)),AIMAG(bps(istep,itheta)),
     $           REAL(xrs(istep,itheta)),AIMAG(xrs(istep,itheta)),
     $           REAL(xzs(istep,itheta)),AIMAG(xzs(istep,itheta)),
     $           REAL(xps(istep,itheta)),AIMAG(xps(istep,itheta))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      DEALLOCATE(psis,rs,zs,phs,brs,bzs,bps,xrs,xzs,xps)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_rzphibx
c-----------------------------------------------------------------------
c     subprogram 12. ipdiag_rzpgrid.
c     diagnose hamada coordinates inverted from rzphi.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_rzpgrid(nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz

      INTEGER :: i,j,itheta
      REAL(r8) :: xint,zint,ttheta,ptheta

      REAL(r8), DIMENSION(0:mthsurf) :: thetas,rbarr,etarr

      TYPE(spline_type) :: rbeta

      ALLOCATE(gdr(0:nr,0:nz),gdz(0:nr,0:nz),gdl(0:nr,0:nz),
     $     gdpsi(0:nr,0:nz),gdthe(0:nr,0:nz),gdphi(0:nr,0:nz))
c-----------------------------------------------------------------------
c     invert given rzphi to hamada coordinates.
c-----------------------------------------------------------------------
      WRITE(*,*)"Diagnosing mapping procedure for rzpgrid"
      gdr=0
      gdz=0
      gdl=0
      gdpsi=0
      gdthe=0
      gdphi=0

      xint=(psi_in%xs(mr)-psi_in%xs(0))/nr
      zint=(psi_in%ys(mz)-psi_in%ys(0))/nz

      CALL ascii_open(out_unit,"ipdiag_rzpgrid.out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RZPGRID: "//
     $     "Give hamada coordinates as functions of rzphi"
      WRITE(out_unit,'(6(1x,a12))')"r","z","limit","psi","theta","phi"

      ! information of plasma boundary
      CALL spline_alloc(rbeta,mthsurf,1)
      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psilim,theta(itheta),0)
         rbarr(itheta)=SQRT(rzphi%f(1))
         etarr(itheta)=theta(itheta)+rzphi%f(2)
      ENDDO
      rbeta%xs=etarr
      rbeta%fs(:,1)=rbarr
      CALL spline_fit(rbeta,"periodic")

      DO i=0,nr 
         DO j=0,nz
            gdr(i,j)=psi_in%xs(0)+i*xint
            gdz(i,j)=psi_in%ys(0)+j*zint

            ! compare grid with equilibrium input
            IF ((nr==mr) .AND. (nz==mz)) THEN
               gdpsi(i,j)=psi_in%fs(i,j,1)
            ELSE
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),0)
               gdpsi(i,j)=psi_in%f(1)
            ENDIF

            IF (gdpsi(i,j)<psilow) gdpsi(i,j)=psilow
            ttheta=ATAN2((gdz(i,j)-zo),(gdr(i,j)-ro))
            IF (ttheta >= 0) THEN 
               ptheta=ttheta/twopi
            ELSE
               ptheta=1+ttheta/twopi
            ENDIF
            
            IF (gdpsi(i,j)<psilim) THEN
               CALL spline_eval(rbeta,ptheta,0)
               IF (SQRT((gdr(i,j)-ro)**2+(gdz(i,j)-zo)**2)<rbeta%f(1))
     $              THEN
                  gdl(i,j)=1
                  DO itheta=0,mthsurf
                     CALL bicube_eval(rzphi,gdpsi(i,j),theta(itheta),0)
                     thetas(itheta)=(theta(itheta)+rzphi%f(2))
                  ENDDO
                  gdthe(i,j)=issect(mthsurf,theta(:),thetas(:),ptheta)
                  CALL bicube_eval(rzphi,gdpsi(i,j),gdthe(i,j),0)
                  gdphi(i,j)=-rzphi%f(3)/twopi
               ENDIF
            ENDIF
            
            WRITE(out_unit,'(6(1x,es12.3))')gdr(i,j),gdz(i,j),gdl(i,j),
     $           gdpsi(i,j),gdthe(i,j),gdphi(i,j)
            
         ENDDO
      ENDDO
         
      CALL ascii_close(out_unit)
      DEALLOCATE(gdr,gdz,gdl,gdpsi,gdthe,gdphi)
      CALL spline_dealloc(rbeta)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_rzpgrid
c-----------------------------------------------------------------------
c     subprogram 13. ipdiag_rzpdiv.
c     check divergence of rzphi functions.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_rzpdiv(nr,nz,lval,rval,zval,fr,fz,fp)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz
      INTEGER, DIMENSION(0:nr,0:nz), INTENT(IN) :: lval
      REAL(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: rval,zval
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: fr,fz,fp
 
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: div

      TYPE(bicube_type) :: rfr,ifr,rfz,ifz

      WRITE(*,*)"Diagnosing divergence for brzphi"

      CALL bicube_alloc(rfr,nr,nz,1)
      CALL bicube_alloc(ifr,nr,nz,1)
      CALL bicube_alloc(rfz,nr,nz,1)
      CALL bicube_alloc(ifz,nr,nz,1)

      rfr%xs=rval(:,0)
      ifr%xs=rval(:,0)
      rfz%xs=rval(:,0)
      ifz%xs=rval(:,0) 

      rfr%ys=zval(0,:)
      ifr%ys=zval(0,:)
      rfz%ys=zval(0,:)
      ifz%ys=zval(0,:) 

      rfr%fs(:,:,1)=rval*REAL(fr)
      ifr%fs(:,:,1)=rval*AIMAG(fr)
      rfz%fs(:,:,1)=REAL(fz)
      ifz%fs(:,:,1)=AIMAG(fz)

      CALL bicube_fit(rfr,"extrap","extrap")
      CALL bicube_fit(ifr,"extrap","extrap")
      CALL bicube_fit(rfz,"extrap","extrap")
      CALL bicube_fit(ifz,"extrap","extrap")      

      DO i=0,nr
         DO j=0,nz
            CALL bicube_eval(rfr,rval(i,j),zval(i,j),1)
            CALL bicube_eval(ifr,rval(i,j),zval(i,j),1)
            CALL bicube_eval(rfz,rval(i,j),zval(i,j),1)
            CALL bicube_eval(ifz,rval(i,j),zval(i,j),1)

            div(i,j)=rfr%fx(1)/rval(i,j)+rfz%fy(1)+
     $           nn*AIMAG(fp(i,j))/rval(i,j)+ifac*
     $           (ifr%fx(1)/rval(i,j)+ifz%fy(1)-
     $           nn*REAL(fp(i,j))/rval(i,j))

         ENDDO
      ENDDO
       
      CALL bicube_dealloc(rfr)
      CALL bicube_dealloc(ifr)
      CALL bicube_dealloc(rfz)
      CALL bicube_dealloc(ifz)

      CALL ascii_open(out_unit,"ipdiag_rzpdiv_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      
      WRITE(out_unit,*)"IPEC_RZPDIV: Divergence in rzphi grid"
      WRITE(out_unit,'(1x,a2,5(a16))')"l","r","z",
     $     "re(div)","im(div)","abs(div)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,5(es16.8))')
     $           lval(i,j),rval(i,j),zval(i,j),
     $           REAL(div(i,j)),AIMAG(div(i,j)),ABS(div(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)   
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_rzpdiv     
c-----------------------------------------------------------------------
c     subprogram 14. ipdiag_radvar.
c     generate various radial variables.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_radvar
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: istep
      REAL(r8) :: qintb
      REAL(r8), DIMENSION(0:mstep) :: psitor,rhotor
      TYPE(spline_type) :: qs

      CALL spline_alloc(qs,mstep,1)      
      qs%xs=psifac
      qs%fs(:,1)=qfac
      CALL spline_fit(qs,"extrap") 
      CALL spline_int(qs)
      qintb=qs%fsi(mstep,1)
      psitor(:)=qs%fsi(:,1)/qintb
      rhotor(:)=SQRT(psitor(:))
      CALL spline_dealloc(qs)
      CALL ascii_open(out_unit,"ipdiag_radvar_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RADVAR: "//
     $     "Various radial variables"
      WRITE(out_unit,'(1x,a8,1x,I6)')"mstep=",mstep
      WRITE(out_unit,'(1x,a8,1x,es16.8)')"qintb=",qintb
      WRITE(out_unit,'(4(1x,a16))')"psi","rho","psitor","rhotor"
      DO istep=1,mstep
         WRITE(out_unit,'(4(1x,es16.8))')psifac(istep),rhofac(istep),
     $        psitor(istep),rhotor(istep)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_radvar

      END MODULE ipdiag_mod
