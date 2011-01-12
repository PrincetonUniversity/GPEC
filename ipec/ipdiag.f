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
c      3. ipdiag_energy
c      4. ipdiag_respmat
c      5. ipdiag_arbsurf
c      6. ipdiag_angles
c      7. ipdiag_surfmode
c      8. ipdiag_extfld
c      9. ipdiag_singcurs
c     10. ipdiag_xbcontra
c     11. ipdiag_xbnormal
c     12. ipdiag_xbnovc
c     13. ipdiag_xbnobo
c     14. ipdiag_xbnorm
c     15. ipdiag_pmodbst
c     16. ipdiag_pmodbrz
c     17. ipdiag_rzphibx
c     18. ipdiag_rzpgrid
c     19. ipdiag_rzpdiv
c     20. ipdiag_radvar
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
     $     "diagnose DCON eigenvalues for eigenvectors"
      WRITE(out_unit,'(2(2x,a3),2x,a12)')"m","m","pot energy"

      DO i=1,mpert
         DO j=1,mpert
            xrmn=0
            ximn=0        
            xrmn(i)=1.0
            ximn(j)=1.0
            xinmn=xrmn+ifac*ximn

            potengy(i,j)=SUM(CONJG(xinmn)*MATMUL(rt,xinmn))
            WRITE(out_unit,'(2(2x,I3),2x,es12.3)')
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

      CALL ascii_open(out_unit,"ipdiag_magpot_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_MAGPOT: magnetic potential errors"
      WRITE(out_unit,'(2x,a8,2x,I4)')"mpert:",mpert
      WRITE(out_unit,'(2x,a4,2(2x,a12))')"mode","chperr1","chperr2"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,2(2x,e12.3))')i,chperr(1,i),chperr(2,i)
      ENDDO
      WRITE(out_unit,'(2x,a4,4(2x,a12))')
     $     "mode","chpsqr1","chpsqr2","chpsqr3","chpsqr4"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,4(2x,e12.3))')
     $        i,chpsqr(1,i),chpsqr(2,i),chpsqr(3,i),chpsqr(4,i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_magpot
c-----------------------------------------------------------------------
c     subprogram 3. ipdiag_energy.
c     write and diagnose energy information.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_energy
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i
c-----------------------------------------------------------------------
c     write eigenenergies.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipdiag_energy_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_ENERGY: "//
     $     "comparison between energy from DCON eigenmodes and "//
     $     "IPEC surface eigenmodes"
      WRITE(out_unit,'(2x,a8,2x,I4)')"mpert:",mpert
      WRITE(out_unit,'(2x,a4,7(2x,a12))')"mode","ee","surfee",
     $     "ep","surfep1","surfep2","surfep3","surfep4"

      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,7(2x,e12.3))')i,ee(i),surfee(i),
     $        ep(i),surfep(1,i),surfep(2,i),surfep(3,i),surfep(4,i)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_energy
c-----------------------------------------------------------------------
c     subprogram 4. ipdiag_respmat.
c     write information about plasma response matrices.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_respmat
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j

      CALL ipresp_reluct
      CALL ipresp_indrel

      CALL ascii_open(out_unit,"ipdiag_respmat_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RESPMAT:"//
     $     "eigenvalues of inductances, permeability and "//
     $     "reluctance matrices"
      WRITE(out_unit,'(2x,a8,2x,I4)')"mpert: ",mpert
      WRITE(out_unit,*)"INDUCTANCE EIGENVALES"
      WRITE(out_unit,'(2x,a4,6(2x,a12))')"mode","surf_ind","plas_ind0",
     $     "plas_ind1","plas_ind2","plas_ind3","plas_ind4"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,6(2x,e12.3))')i,surf_indev(i),
     $        plas_indev(0,i),plas_indev(1,i),plas_indev(2,i),
     $        plas_indev(3,i),plas_indev(4,i)
      ENDDO
      WRITE(out_unit,*)"PERMEABILITY EIGENVALUES"
      WRITE(out_unit,'(2x,a4,5(2x,a12))')"mode","permeab0",
     $     "permeab1","permeab2","permeab3","permeab4"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,5(2x,e12.3))')i,
     $        REAL(permeabev(0,permeabindex(0,i))),
     $        REAL(permeabev(1,permeabindex(1,i))),
     $        REAL(permeabev(2,permeabindex(2,i))),
     $        REAL(permeabev(3,permeabindex(3,i))),
     $        REAL(permeabev(4,permeabindex(4,i)))
      ENDDO

      WRITE(out_unit,*)"RELUCTANCE EIGENVALUES"     
      WRITE(out_unit,'(2x,a4,5(2x,a12))')"mode","reluct0",
     $     "reluct1","reluct2","reluct3","reluct4"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,5(2x,e12.3))')i,reluctev(0,i),
     $        reluctev(1,i),reluctev(2,i),reluctev(3,i),reluctev(4,i)
      ENDDO

      WRITE(out_unit,*)"INDRELUCTANCE EIGENVALUES"     
      WRITE(out_unit,'(2x,a4,5(2x,a12))')"mode","indrel0",
     $     "indrel1","indrel2","indrel3","indrel4"
      DO i=1,mpert
         WRITE(out_unit,'(2x,I4,5(2x,e12.3))')i,indrelev(0,i),
     $        indrelev(1,i),indrelev(2,i),indrelev(3,i),indrelev(4,i)
      ENDDO

      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_respmat      
c-----------------------------------------------------------------------
c     subprogram 5. ipdiag_arbsurf.
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
     $     "diagnose surface inductances for arbitrary shape"
      WRITE(out_unit,'(2x,a4,2x,a12)')"mode","surf_ind"
      DO i=1,vmpert
         WRITE(out_unit,'(2x,I4,2x,es12.3)')i,vsurf_indev(i)
      ENDDO      
      CALL ascii_close(out_unit)
      DEALLOCATE(vsurf_indmats,vsurf_indev)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_arbsurf
c-----------------------------------------------------------------------
c     subprogram 6. ipdiag_angles.
c     diagnose and visualize magnetic angles.
c     __________________________________________________________________
c     rstep   : number of radial points in a fixed angle
c     angnum  : number of angles for in a fixed radial point
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
     $     "diagnose and visualize magnetic angles"
      WRITE(out_unit,'(2x,a8,2x,I6)')"angnum:",angnum
      WRITE(out_unit,'(2x,a8,2x,I6)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,'(3(2x,a12))')"psi","plerror","haerror"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(3(2x,e12.3))')
     $           psi(istep),plerror(istep,itheta),haerror(istep,itheta)
         ENDDO
      ENDDO
      
      DO iqty=1,5
         WRITE(out_unit,'(2x,a8,2x,I2)')"iqty:",iqty
         WRITE(out_unit,'(3(2x,a12))')"r","z","angles"
         DO inum=0,angnum-1
            DO istep=1,rstep
               WRITE(out_unit,'(3(2x,e12.3))')
     $              rs(istep,inum,iqty),zs(istep,inum,iqty),angles(inum)
            ENDDO
         ENDDO
      ENDDO

      WRITE(out_unit,'(7(2x,a12))')"psi","theta","polar","hamada",
     $     "pest","equalarc","boozer"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(7(2x,e12.3))')
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
c     subprogram 7. ipdiag_surfmode.
c     reponse to fourier modes for the control surface.
c     __________________________________________________________________
c     lowmode  : lowest number of m fourier mode applied
c     highmode : highest number of m fourier mode applied
c     polo     : poloidal angle coordinates
c     toro     : toroidal angle coordinates
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_surfmode(lowmode,highmode,polo,toro)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: lowmode,highmode,polo,toro

      INTEGER :: i,j,mnum
      CHARACTER(1) :: spolo,storo

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
         CALL ipeq_cotoha(psilim,finmn(j,:),mfac,mpert,polo,toro)
         CALL ipeq_weight(psilim,finmn(j,:),mfac,mpert,1)         
         boutmn(j,:)=MATMUL(permeabmats(modelnum,:,:),finmn(j,:))
         CALL ipeq_weight(psilim,finmn(j,:),mfac,mpert,0)  
         CALL ipeq_hatoco(psilim,boutmn(j,:),mfac,mpert,polo,toro)
         CALL iscdftb(mfac,mpert,boutfun(j,:),mthsurf,boutmn(j,:))
      ENDDO
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      CALL ascii_open(out_unit,"ipdiag_surfmode_p"//spolo//"_t"
     $     //storo//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_SURFMODE: plasma response for the fourier 
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
      DEALLOCATE(binmn,boutmn,finmn,binfun,boutfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_surfmode
c-----------------------------------------------------------------------
c     subprogram 8. ipdiag_extfld
c     change coordinates on the boundary.
c     __________________________________________________________________
c     infile   : input file name containing rawdata
c     errtype  : input file format
c     binmn    : given error field spectrum.
c     poloin   : input poloidal angle
c     toroin   : input toroidal angle
c     poloout  : output poloidal angle
c     toroout  : output toroidal angle
c     __________________________________________________________________
c     boutmn   : new error field spectrum
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_extfld(rinfile,formattype,left,scale,
     $     binmn,poloin,toroin,poloout,toroout,boutmn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: left,poloin,toroin,poloout,toroout,labl
      REAL(r8), INTENT(IN) :: scale      
      CHARACTER(128), INTENT(IN) :: rinfile,formattype
      COMPLEX(r8), DIMENSION(lmpert), INTENT(INOUT) :: binmn
      COMPLEX(r8), DIMENSION(lmpert), INTENT(OUT) :: boutmn

      INTEGER :: i,j,i1,i2,i3,ms,hfsurf
      REAL(r8) :: htheta
      CHARACTER(1) :: spoloin,storoin,spoloout,storoout,slabl
      CHARACTER(128) :: message

      COMPLEX(r8), DIMENSION(lmpert) :: ftnmn,finmn,foutmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun

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
         binmn=rawmn(:,nn)
         DEALLOCATE(cosmn,sinmn,rawmn,hawmn)
      ENDIF
c-----------------------------------------------------------------------
c     get the plasma response on the control surface (CHECK!!!)
c-----------------------------------------------------------------------
      IF (scale /= 0) binmn=binmn*scale
      ftnmn=binmn
      finmn=binmn
      CALL ipeq_cotoha(psilim,ftnmn,lmfac,lmpert,poloin,toroin)
      CALL ipeq_cotoha(psilim,finmn,lmfac,lmpert,poloin,toroin)
      foutmn=ftnmn
      CALL ipeq_hatoco(psilim,ftnmn,lmfac,lmpert,poloout,toroout)

      CALL ipeq_weight(psilim,finmn,lmfac,lmpert,poloin)
      CALL ipeq_weight(psilim,foutmn,lmfac,lmpert,poloout)
      CALL ipeq_hatoco(psilim,finmn,lmfac,lmpert,poloin,toroin)
      CALL ipeq_hatoco(psilim,foutmn,lmfac,lmpert,poloout,toroout)
      boutmn=ftnmn

      CALL iscdftb(lmfac,lmpert,binfun,mthsurf,binmn)
      CALL iscdftb(lmfac,lmpert,boutfun,mthsurf,boutmn)
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=spoloin, FMT='(I1)')poloin
      WRITE(UNIT=storoin, FMT='(I1)')toroin
      WRITE(UNIT=spoloout, FMT='(I1)')poloout
      WRITE(UNIT=storoout, FMT='(I1)')toroout
      WRITE(UNIT=slabl, FMT='(I1)')labl

      CALL ascii_open(out_unit,"ipdiag_extfld_p"//spoloin//spoloout//
     $     "_t"//storoin//storoout//"_l"//slabl//
     $     "_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_EXTFLD: "//
     $     "external perturbations on the control surface"
      WRITE(out_unit,'(2x,a12,2x,I6)')"lmpert:",lmpert
      WRITE(out_unit,'(2x,a12,2x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(2x,a6,4(2x,a12))')"m","rebin","imbin",
     $     "rebout","imbout"
      DO i=1,lmpert
         WRITE(out_unit,'(2x,I6,4(2x,e12.3))')lmfac(i),
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

      CALL ascii_open(out_unit,"ipdiag_extflx_p"//spoloin//spoloout//
     $     "_t"//storoin//storoout//"_l"//slabl//
     $     "_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_EXTFLX: "//
     $     "external fluxes on the control surface"
      WRITE(out_unit,'(2x,a12,2x,I6)')"lmpert:",lmpert
      WRITE(out_unit,'(2x,a12,2x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,*)"MODES"
      WRITE(out_unit,'(2x,a6,4(2x,a12))')"m","refin","imfin",
     $     "refout","imfout"
      DO i=1,lmpert
         WRITE(out_unit,'(2x,I6,4(2x,e12.3))')lmfac(i),
     $        REAL(finmn(i)),AIMAG(finmn(i)),
     $        REAL(foutmn(i)),AIMAG(foutmn(i))
      ENDDO
      CALL ascii_close(out_unit)     
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_extfld
c-----------------------------------------------------------------------
c     subprogram 9. ipdiag_singcurs.
c     diagnose asymtotic values of singular currents.
c     __________________________________________________________________
c     egnjm      : eigenmode number without edge_flag
c     xwpimn     : edge deformation input when edge_flag
c     rsing      : number of rational surfaces
c     resol      : resolution by number of grid points
c     smallwidth : closest point to approach
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_singcurs(egnum,xwpimn,rsing,resol,smallwidth,
     $     labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,resol,rsing,labl
      REAL(r8), INTENT(IN) :: smallwidth
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn
      CHARACTER(1) :: slabl

      INTEGER :: ising,i,itheta,resnum
      REAL(r8) :: respsi,rrespsi,lrespsi,bigwidth,astep,lpsi,rpsi,
     $     sqrpsi,correc,lq1,rq1
      COMPLEX(r8) :: lnbwp1mn,rnbwp1mn

      INTEGER, DIMENSION(rsing) :: logsteps
      REAL(r8), DIMENSION(rsing) :: j_c
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: lcormn,rcormn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: bwp_fun,lcorfun,rcorfun

      REAL(r8), DIMENSION(rsing,resol*10) :: dist
      COMPLEX(r8), DIMENSION(rsing,resol*10) :: deltas,delcurs,corcurs,
     $     singcurs

      dist=0
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
     $        (chi1*sq%f(4))**2/mu0
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
            CALL ipeq_contra(lpsi)
            lnbwp1mn=nbwp1_mn(resnum)
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
               correc=jac*(w(1,1)*w(2,1)+w(1,2)*w(2,2)-
     $              (w(1,1)*w(3,1)+w(1,2)*w(3,2))/sq%f(4))/
     $              (sqrpsi*sq%f(4)*chi1)
               lcorfun(itheta)=bwp_fun(itheta)*correc
            ENDDO
            CALL iscdftf(mfac,mpert,lcorfun,mthsurf,lcormn)
c-----------------------------------------------------------------------
c     calculate delta and correction term at rpsi.
c-----------------------------------------------------------------------
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

            dist(ising,i)=astep/(nn*ABS(singtype(ising)%q1))
            deltas(ising,i)=(rnbwp1mn-lnbwp1mn)/twopi
            delcurs(ising,i)=j_c(ising)*
     $           ifac/mfac(resnum)*deltas(ising,i)
            corcurs(ising,i)=-j_c(ising)*
     $           (rcormn(resnum)-lcormn(resnum))
            singcurs(ising,i)=delcurs(ising,i)-corcurs(ising,i)
         ENDDO

         WRITE(*,*)"Finished the analysis for the q =",singtype(ising)%q
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write the results.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl, FMT='(I1)')labl    
      CALL ascii_open(out_unit,"ipdiag_deltas_l"//slabl//"_n"//
     $     sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_DELTAS: "//
     $     "asymtotic analysis for deltas at rational surfaces"
      WRITE(out_unit,'(2x,a8,2x,I6)')"rsing:",rsing
      WRITE(out_unit,'(2x,a8,2x,I6)')"resol:",resol
      WRITE(out_unit,'(3(2x,a16))')"distance","re(delta)","im(delta)"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,a2,2x,f6.3)')"q=",singtype(ising)%q
         WRITE(out_unit,'(2x,a12,2x,I6)')"logsteps:",logsteps(ising)
         DO i=1,logsteps(ising)
            WRITE(out_unit,'(3(2x,e16.9))')dist(ising,i),
     $           REAL(deltas(ising,i)),AIMAG(deltas(ising,i))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
      CALL ascii_open(out_unit,"ipdiag_singcurs_l"//slabl//"_n"//
     $     sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_SINGCURS: "//
     $     "asymtotic analysis for singular currents at "//
     $     "rational surfaces"
      WRITE(out_unit,'(2x,a8,2x,I6)')"rsing:",rsing
      WRITE(out_unit,'(2x,a8,2x,I6)')"resol:",resol
      WRITE(out_unit,'(7(2x,a16))')"distance","re(delcur)","im(delcur)",
     $     "re(corcur)","im(corcur)","re(singcur)","im(singcur)"
      DO ising=1,rsing
         WRITE(out_unit,'(2x,a2,2x,f6.3)')"q=",singtype(ising)%q
         WRITE(out_unit,'(2x,a12,2x,I6)')"logsteps:",logsteps(ising)
         DO i=1,logsteps(ising)
            WRITE(out_unit,'(7(2x,e16.9))')dist(ising,i),
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
c     subprogram 10. ipdiag_xbcontra.
c     write contravariant componets of xi.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbcontra(egnum,xwpimn,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert
      CHARACTER(1) :: spolo,storo,slabl
      CHARACTER(2) :: slabl2
      COMPLEX(r8), DIMENSION(mpert) :: xwd_mn
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      WRITE(*,*)"Computing contravariant components in detail"
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipdiag_xbcontra_p"//spolo//"_t"//
     $        storo//"_l"//slabl//"_n"//sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipdiag_xbcontra_p"//spolo//"_t"//
     $        storo//"_l"//slabl2//"_n"//sn//".out","UNKNOWN")
      ENDIF  
      WRITE(out_unit,*)"IPDIAG_XBCONTRA: "//
     $     "contravariant components of displacement and field"
      WRITE(out_unit,'(2x,a8,2x,I6)')"mstep:",mstep
      WRITE(out_unit,'(2x,a8,2x,I6)')"mlow:",mlow
      WRITE(out_unit,'(2x,a8,2x,I6)')"mhigh:",mhigh
      WRITE(out_unit,'(2(2x,a16))')"psifac","qfac"
      CALL ipeq_alloc
      DO istep=0,mstep
         WRITE(out_unit,'(2(2x,e16.9))')psifac(istep),qfac(istep)
      ENDDO
      WRITE(out_unit,'(18(2x,a16))')
     $     "xwp(real)","xwp(imag)","xwt(real)","xwt(imag)",
     $     "xwz(real)","xwz(imag)","bwp(real)","bwp(imag)",
     $     "bwt(real)","bwt(imag)","bwz(real)","bwz(imag)",
     $     "xvs(real)","xvs(imag)","xwd(real)","xwd(imag)",
     $     "xwp1(real)","xwp1(imag)"
      DO istep=0,mstep
         CALL ipeq_contra(psifac(istep))
         CALL cspline_eval(u1,psifac(istep),1)
         xwd_mn(:)=u1%f1(:)
         IF ((polo /= 1).OR.(toro /= 1)) THEN 
            CALL ipeq_hatoco(psifac(istep),xwp_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),xwt_mn,mfac,mpert,polo,toro) 
            CALL ipeq_hatoco(psifac(istep),xwz_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),bwp_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),bwt_mn,mfac,mpert,polo,toro) 
            CALL ipeq_hatoco(psifac(istep),bwz_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),xvs_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),xwd_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),xwp1_mn,mfac,mpert,polo,toro)
         ENDIF
         DO ipert=1,mpert
            WRITE(out_unit,'(18(2x,e16.9))')
     $           REAL(xwp_mn(ipert)),AIMAG(xwp_mn(ipert)),
     $           REAL(xwt_mn(ipert)),AIMAG(xwt_mn(ipert)),
     $           REAL(xwz_mn(ipert)),AIMAG(xwz_mn(ipert)),
     $           REAL(bwp_mn(ipert)),AIMAG(bwp_mn(ipert)),
     $           REAL(bwt_mn(ipert)),AIMAG(bwt_mn(ipert)),
     $           REAL(bwz_mn(ipert)),AIMAG(bwz_mn(ipert)),
     $           REAL(xvs_mn(ipert)),AIMAG(xvs_mn(ipert)),
     $           REAL(xwd_mn(ipert)),AIMAG(xwd_mn(ipert)),
     $           REAL(xwp1_mn(ipert)),AIMAG(xwp1_mn(ipert))
         ENDDO
      ENDDO
      CALL ipeq_dealloc
      CALL ascii_close(out_unit)

      RETURN
      END SUBROUTINE ipdiag_xbcontra
c-----------------------------------------------------------------------
c     subprogram 11. ipdiag_xbnormal.
c     write normal componets of xi and b.
c     __________________________________________________________________
c     osol   : label of an eigenmode
c     edgemn : xipsi components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbnormal(egnum,xwpimn,polo,toro,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,polo,toro,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert
      CHARACTER(1) :: spolo,storo,slabl
      CHARACTER(2) :: slabl2      
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
      WRITE(*,*)"Computing normal components in detail"
      WRITE(UNIT=spolo, FMT='(I1)')polo
      WRITE(UNIT=storo, FMT='(I1)')toro
      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl            
         CALL ascii_open(out_unit,"ipdiag_xbnormal_p"//spolo//"_t"//
     $        storo//"_l"//slabl//"_n"//sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl            
         CALL ascii_open(out_unit,"ipdiag_xbnormal_p"//spolo//"_t"//
     $        storo//"_l"//slabl2//"_n"//sn//".out","UNKNOWN")
      ENDIF
      WRITE(out_unit,*)"IPDIAG_XBNROMAL: "//
     $     "normal components of displacement and field"
      WRITE(out_unit,'(2x,a8,2x,I6)')"mstep:",mstep
      WRITE(out_unit,'(2x,a8,2x,I6)')"mlow:",mlow
      WRITE(out_unit,'(2x,a8,2x,I6)')"mhigh:",mhigh
      WRITE(out_unit,'(2(2x,a16))')"psifac","qfac"
      DO istep=0,mstep
         WRITE(out_unit,'(2(2x,e16.9))')psifac(istep),qfac(istep)
      ENDDO
      WRITE(out_unit,'(4(2x,a16))')
     $     "xno(real)","xno(imag)","bno(real)","bno(imag)"
      CALL ipeq_alloc
      DO istep=0,mstep
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_normal(psifac(istep))
         IF ((polo /= 1).OR.(toro /= 1)) THEN 
            CALL ipeq_hatoco(psifac(istep),xno_mn,mfac,mpert,polo,toro)
            CALL ipeq_hatoco(psifac(istep),bno_mn,mfac,mpert,polo,toro)
         ENDIF
         DO ipert=1,mpert
            WRITE(out_unit,'(4(2x,e16.9))')
     $           REAL(xno_mn(ipert)),AIMAG(xno_mn(ipert)),
     $           REAL(bno_mn(ipert)),AIMAG(bno_mn(ipert))        
         ENDDO
      ENDDO
      CALL ipeq_dealloc
      CALL ascii_close(out_unit)

      RETURN
      END SUBROUTINE ipdiag_xbnormal
c-----------------------------------------------------------------------
c     subprogram 12. ipdiag_xbnovc.
c     write data for perturbed flux surfaces.
c     __________________________________________________________________
c     egnum   : label of an eigenmode
c     xwpimn  : xwp_mn components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbnovc(egnum,xwpimn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      CHARACTER(1) :: slabl
      INTEGER :: i,j,istep,ithnum,ising,rstep,stepnum,mthnumb
      REAL(r8) :: lpsi,rpsi,step,range,idist,rdist,bigstep,smlstep,
     $     wegtfac

      INTEGER, DIMENSION(:), POINTER :: mthnum
      REAL(r8), DIMENSION(:), POINTER :: psis,ppsis,rs,zs,
     $     mthang,delpsi,no_rvec,no_zvec
      COMPLEX(r8), DIMENSION(:), POINTER :: xno_fun,bno_fun,
     $     xno_rvc,xno_zvc,bno_rvc,bno_zvc


      idist=0.002
      rdist=0.01
      bigstep=0.10
      smlstep=0.01
      mthnumb=500
      wegtfac=10.0      
      ALLOCATE(ppsis(100))
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xwpimn)
c-----------------------------------------------------------------------
c     assign proper radial points.
c-----------------------------------------------------------------------
      i=1
      ppsis=0
      ppsis(i)=idist
      lpsi=ppsis(i)
      DO ising=1,msing
         rpsi=singtype(ising)%psifac-rdist
         range=rpsi-lpsi
         IF (range .GE. bigstep) THEN
            stepnum=INT(range/bigstep)+1
            step=range/stepnum
            DO j=1,stepnum
               i=i+1
               ppsis(i)=ppsis(i-1)+step
            ENDDO
            lpsi=ppsis(i)
         ELSE IF ((range .GE. smlstep) .AND. (range .LT. bigstep)) THEN
            i=i+1
            ppsis(i)=rpsi
            lpsi=ppsis(i)
         ELSE
            EXIT
         ENDIF
      ENDDO
      i=i+1
      ppsis(i)=psilim
      rstep=i
      ALLOCATE(psis(rstep),mthnum(rstep))
      DO istep=1,rstep
         psis(istep)=ppsis(istep)
      ENDDO
      DEALLOCATE(ppsis)
c-----------------------------------------------------------------------
c     compute movement of flux surfaces.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl, FMT='(I1)')labl
      CALL ascii_open(20,"ipdiag_xbnovc_x_l"//slabl//"_n"//
     $     sn//".out","UNKNOWN")
      WRITE(20,*)"IPDIAG_XBNOVC_X: "//
     $     "perturbed normal displacement vectors"
      WRITE(20,'(7(2x,a12))')"r","z","psi",
     $     "real(xnor)","real(xnoz)","imag(xnor)","imag(xnoz)"

      CALL ascii_open(21,"ipdiag_xbnovc_b_l"//slabl//"_n"//
     $     sn//".out","UNKNOWN")
      WRITE(21,*)"IPDIAG_XBNOVC_B: "//
     $     "perturbed normal b field vectors"
      WRITE(21,'(7(2x,a12))')"r","z","psi",
     $     "real(bnor)","real(bnoz)","imag(bnor)","imag(bnoz)"

      DO istep=1,rstep
         CALL bicube_eval(rzphi,psis(istep),pi/twopi,0)
         rfac=SQRT(rzphi%f(1))
         mthnum(istep)=INT(mthnumb*rfac**2)
         ALLOCATE(rs(0:mthnum(istep)),zs(0:mthnum(istep)),
     $        mthang(0:mthnum(istep)),delpsi(0:mthnum(istep)),
     $        xno_fun(0:mthnum(istep)),bno_fun(0:mthnum(istep)),
     $        xno_rvc(0:mthnum(istep)),xno_zvc(0:mthnum(istep)),
     $        bno_rvc(0:mthnum(istep)),bno_zvc(0:mthnum(istep)),
     $        no_rvec(0:mthnum(istep)),no_zvec(0:mthnum(istep)))
         mthang=(/(i,i=0,mthnum(istep))/)/REAL(mthnum(istep),r8)
         DO ithnum=0,mthnum(istep)
            CALL bicube_eval(rzphi,psis(istep),mthang(ithnum),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(mthang(ithnum)+rzphi%f(2))
            rs(ithnum)=ro+rfac*COS(eta)
            zs(ithnum)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(ithnum)/jac
            w(1,2)=-rzphi%fy(1)*pi*rs(ithnum)/(rfac*jac)
            delpsi(ithnum)=SQRT(w(1,1)**2+w(1,2)**2)
            no_rvec(ithnum)=(cos(eta)*w(1,1)-sin(eta)*w(1,2))/
     $           delpsi(ithnum)
            no_zvec(ithnum)=(sin(eta)*w(1,1)+cos(eta)*w(1,2))/
     $           delpsi(ithnum)
         ENDDO
         CALL ipeq_alloc
         CALL ipeq_contra(psis(istep))
         CALL ipeq_normal(psis(istep))
         CALL ipeq_hatoco(psis(istep),xno_mn,mfac,mpert,1,0)
         CALL ipeq_hatoco(psis(istep),bno_mn,mfac,mpert,1,0)
         CALL iscdftb(mfac,mpert,xno_fun,mthnum(istep),xno_mn)
         CALL iscdftb(mfac,mpert,bno_fun,mthnum(istep),bno_mn)
         CALL ipeq_dealloc
         xno_rvc=xno_fun*no_rvec*wegtfac
         xno_zvc=xno_fun*no_zvec*wegtfac
         bno_rvc=bno_fun*no_rvec*wegtfac
         bno_zvc=bno_fun*no_zvec*wegtfac

         DO ithnum=0,mthnum(istep)
            WRITE(20,'(7(2x,es12.3))')
     $           rs(ithnum),zs(ithnum),psis(istep),
     $           REAL(xno_rvc(ithnum)),REAL(xno_zvc(ithnum)),
     $           AIMAG(xno_rvc(ithnum)),AIMAG(xno_zvc(ithnum))
            WRITE(21,'(7(2x,es12.3))')
     $           rs(ithnum),zs(ithnum),psis(istep),
     $           REAL(bno_rvc(ithnum)),REAL(bno_zvc(ithnum)),
     $           AIMAG(bno_rvc(ithnum)),AIMAG(bno_zvc(ithnum))
         ENDDO
         DEALLOCATE(rs,zs,mthang,delpsi,xno_fun,bno_fun,
     $        xno_rvc,xno_zvc,bno_rvc,bno_zvc,no_rvec,no_zvec)
      ENDDO
      CALL ascii_close(20)
      CALL ascii_close(21)
      DEALLOCATE(psis,mthnum)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_xbnovc
c-----------------------------------------------------------------------
c     subprogram 13. ipdiag_xbnobo.
c     write normal perturbed quantities on the boundary.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xwp_mn components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbnobo(egnum,xwpimn,d3_flag,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) ::xwpimn
      LOGICAL, INTENT(IN) :: d3_flag

      INTEGER :: i,ithnum,mthnum
      INTEGER(r8) :: nnn,mmtheta,np,nt
      REAL(r8) :: ddist
      CHARACTER(1) :: slabl
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
      CALL ipeq_contra(psilim)
      CALL ipeq_normal(psilim)
      CALL ipeq_hatoco(psilim,xno_mn,mfac,mpert,1,0)
      CALL ipeq_hatoco(psilim,bno_mn,mfac,mpert,1,0)
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
      WRITE(UNIT=slabl, FMT='(I1)')labl      
      CALL ascii_open(out_unit,"ipdiag_xnobo_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_XNOBO: "//
     $     "perturbed normal displacement vectors on the boundary"      
      WRITE(out_unit,'(6(2x,a16))')"r","z",
     $     "real(xnor)","real(xnoz)","imag(xnor)","imag(xnoz)"
      DO ithnum=0,mthnum
         WRITE(out_unit,'(6(2x,e16.9))')rs(ithnum),zs(ithnum),
     $        REAL(xnorvc(ithnum)),REAL(xnozvc(ithnum)),
     $        AIMAG(xnorvc(ithnum)),AIMAG(xnozvc(ithnum))
      ENDDO
      CALL ascii_open(out_unit,"ipdiag_bnobo_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_BNOBO: "//
     $     "perturbed normal field vectors on the boundary"      
      WRITE(out_unit,'(6(2x,a16))')"r","z",
     $     "real(bnor)","real(bnoz)","imag(bnor)","imag(bnoz)"
      DO ithnum=0,mthnum
         WRITE(out_unit,'(6(2x,e16.9))')rs(ithnum),zs(ithnum),
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
         WRITE(out_unit,'(1p,5e16.8)')rs
         WRITE(out_unit,'(1p,5e16.8)')zs
         WRITE(out_unit,'(1p,5e16.8)')REAL(xno_fun)
         WRITE(out_unit,'(1p,5e16.8)')AIMAG(xno_fun)      
         CALL ascii_close(out_unit)
         
         filein="iptemp.txt"
         file1="ipidl_3dsurf_xnobo_l"//slabl//"_n"//sn//".out"

         ddist=0.1
         mmtheta=mthnum
         nnn=nn
         np=144
         nt=144
         CALL ipidl_3dsurf(filein,nnn,mmtheta,np,nt,ddist,file1)

         CALL ascii_open(out_unit,"iptemp.txt","UNKNOWN")
         WRITE(out_unit,'(1p,5e16.8)')rs
         WRITE(out_unit,'(1p,5e16.8)')zs
         WRITE(out_unit,'(1p,5e16.8)')REAL(bno_fun)
         WRITE(out_unit,'(1p,5e16.8)')AIMAG(bno_fun)      
         CALL ascii_close(out_unit)
         
         filein="iptemp.txt"
         file1="ipidl_3dsurf_bnobo_l"//slabl//"_n"//sn//".out"

         ddist=0.1
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
c     subprogram 14. ipdiag_xbnorm.
c     compute normal perturbed quantities.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xwp_mn components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_xbnorm(egnum,xwpimn,polo,toro,labl)
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
      WRITE(*,*)"Computing normal quantities"
      CALL ipeq_alloc
      DO istep=1,rstep
         CALL ipeq_contra(psis(istep))
         CALL ipeq_normal(psis(istep))
         rs(istep,:)=r
         zs(istep,:)=z
         xnomns(istep,:)=xno_mn
         bnomns(istep,:)=bno_mn
         CALL ipeq_hatoco(psis(istep),xnomns(istep,:),mfac,mpert,
     $        polo,toro)
         CALL ipeq_hatoco(psis(istep),bnomns(istep,:),mfac,mpert,
     $        polo,toro)
         CALL ipeq_hatoco(psis(istep),xno_mn,mfac,mpert,1,0)
         CALL ipeq_hatoco(psis(istep),bno_mn,mfac,mpert,1,0)
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
      CALL ascii_open(out_unit,"ipdiag_xbnormmn_p"//spolo//"_t"//storo
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
      CALL ascii_open(out_unit,"ipdiag_xbnormrz_l"//slabl//"_n"
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
      END SUBROUTINE ipdiag_xbnorm
c-----------------------------------------------------------------------
c     subprogram 15. ipdiag_pmodbst.
c     compute strength of perturbed mod b.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xipsi components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_pmodbst(egnum,xwpimn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,itheta
      CHARACTER(1) :: slabl


      REAL(r8), DIMENSION(0:mthsurf) :: seulmodb,seulbpar,
     $     sxiwobbp,sxiwobbr,sxiwobbt,slagbpar
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: eulmodb,eulbpar,
     $     xiwobbp,xiwobbr,xiwobbt,lagbpar,xwp_fun,xwt_fun,
     $     bwp_fun,bwt_fun,bwz_fun,bvp_fun,bvt_fun,bvz_fun

      REAL(r8), DIMENSION(:), POINTER :: psis,
     $     eulmodbst,eulbparst,xiwobbpst,xiwobbrst,xiwobbtst,lagbparst,
     $     invlagbst,ntv_eulbparst,ntv_lagbparst

c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing perturbed b field strength"
      ALLOCATE(psis(rstep),eulmodbst(rstep),
     $     eulbparst(rstep),xiwobbpst(rstep),xiwobbrst(rstep),
     $     xiwobbtst(rstep),lagbparst(rstep),invlagbst(rstep),
     $     ntv_eulbparst(rstep),ntv_lagbparst(rstep))

      IF (rstep .EQ. mstep) THEN
         psis=psifac
      ELSE
         psis=(/(istep,istep=1,rstep)/)/REAL(rstep,r8)*
     $        (psilim-bdist/singtype(msing)%q)
      ENDIF

      CALL idcon_build(egnum,xwpimn)
      
      CALL ipeq_alloc
      DO istep=1,rstep
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces.
c-----------------------------------------------------------------------
         CALL spline_eval(ffun,psis(istep),0)
         CALL spline_eval(sq,psis(istep),1)
         CALL ipeq_contra(psis(istep))
         CALL ipeq_cova(psis(istep))

         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xwt_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bwt_mn)
         CALL iscdftb(mfac,mpert,bwz_fun,mthsurf,bwz_mn)
         CALL iscdftb(mfac,mpert,bvp_fun,mthsurf,bvp_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         eulmodb=SQRT(ABS(REAL(CONJG(bwp_fun)*bvp_fun+
     $        CONJG(bwt_fun)*bvt_fun+CONJG(bwz_fun)*bvz_fun)))
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psis(istep),theta(itheta),1)
            eulbpar(itheta)=(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(ffun%f(1)*eqfun%f(1))
            xiwobbp(itheta)=xwp_fun(itheta)*sq%f1(2)/eqfun%f(1)
            xiwobbr(itheta)=xwp_fun(itheta)*eqfun%fx(1)
            xiwobbt(itheta)=xwt_fun(itheta)*eqfun%fy(1)
         ENDDO
         lagbpar=eulbpar+xiwobbr+xiwobbt
c-----------------------------------------------------------------------
c     compute the ntv strength terms in hamada.
c-----------------------------------------------------------------------
         CALL iscdftf(mfac,mpert,eulbpar,mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbpar,mthsurf,lagbpar_mn)
         ntv_eulbparst(istep)=nn**2*SUM(ABS(eulbpar_mn)**2)
         ntv_lagbparst(istep)=nn**2*SUM(ABS(lagbpar_mn)**2)
c-----------------------------------------------------------------------
c     compute each strength.
c-----------------------------------------------------------------------
         seulmodb=ABS(eulmodb)**2
         seulbpar=ABS(eulbpar)**2
         sxiwobbp=ABS(xiwobbp)**2
         sxiwobbr=ABS(xiwobbr)**2
         sxiwobbt=ABS(xiwobbt)**2       
         slagbpar=ABS(lagbpar)**2
         eulmodbst(istep)=issurfint(seulmodb,mthsurf,psis(istep),0,1)
         eulbparst(istep)=issurfint(seulbpar,mthsurf,psis(istep),0,1)
         xiwobbpst(istep)=issurfint(sxiwobbp,mthsurf,psis(istep),0,1)  
         xiwobbrst(istep)=issurfint(sxiwobbr,mthsurf,psis(istep),0,1) 
         xiwobbtst(istep)=issurfint(sxiwobbt,mthsurf,psis(istep),0,1)
         lagbparst(istep)=issurfint(slagbpar,mthsurf,psis(istep),0,1)
         invlagbst(istep)=1.0/
     $        issurfint(1.0/slagbpar,mthsurf,psis(istep),0,1)
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------      
      WRITE(UNIT=slabl, FMT='(I1)')labl
      CALL ascii_open(out_unit,"ipdiag_pmodbst_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_PMODBST: "//
     $     "perturbed mod b strength on flux surfaces"
      WRITE(out_unit,'(10(2x,a12))')"psi","eulmodbst","eulbparst",
     $     "ntveulbst","xiwobbpst","xiwobbrst","xiwobbtst",
     $     "lagbparst","invlagbst","ntvlagbst"
      DO istep=1,rstep
         WRITE(out_unit,'(10(2x,e12.3))')psis(istep),
     $        eulmodbst(istep),eulbparst(istep),ntv_eulbparst(istep),
     $        xiwobbpst(istep),xiwobbrst(istep),xiwobbtst(istep),
     $        lagbparst(istep),invlagbst(istep),ntv_lagbparst(istep)
      ENDDO
      CALL ascii_close(out_unit)

      DEALLOCATE(psis,eulmodbst,eulbparst,xiwobbpst,xiwobbrst,
     $     xiwobbtst,lagbparst,invlagbst,ntv_eulbparst,ntv_lagbparst)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_pmodbst
c-----------------------------------------------------------------------
c     subprogram 16. ipdiag_pmodbrz.
c     plot perturbed mod b at rz coordinates.
c     __________________________________________________________________
c     egnum  : eigenmode number without edge_flag
c     xwpimn : edge deformation input when edge_flag
c     __________________________________________________________________
c     labl   : label for multiple runs   
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_pmodbrz(egnum,xwpimn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: istep,ipert,itheta
      CHARACTER(1) :: slabl

      COMPLEX(r8), DIMENSION(mpert) :: eulbpar,lagbpar
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,xwt_fun,
     $     bvt_fun,bvz_fun

      REAL(r8), DIMENSION(:), POINTER :: psis
      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs,eqfunx,eqfuny,eqfuns
      COMPLEX(r8), DIMENSION(:,:), POINTER :: eulbpar_mn,lagbpar_mn,
     $     eulbparfun,lagbparfun,xwpfun,xwtfun,bvtfun,bvzfun
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      WRITE(*,*)"Computing perturbed b field"
      ALLOCATE(psis(rstep),rs(rstep,0:mthsurf),zs(rstep,0:mthsurf),
     $     eulbpar_mn(rstep,mpert),lagbpar_mn(rstep,mpert),
     $     eulbparfun(rstep,0:mthsurf),lagbparfun(rstep,0:mthsurf),
     $     xwpfun(rstep,0:mthsurf),xwtfun(rstep,0:mthsurf),
     $     bvtfun(rstep,0:mthsurf),bvzfun(rstep,0:mthsurf),
     $     eqfunx(rstep,0:mthsurf),eqfuny(rstep,0:mthsurf),
     $     eqfuns(rstep,0:mthsurf))

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
         xwpfun(istep,:)=xwp_fun
         xwtfun(istep,:)=xwt_fun
         bvtfun(istep,:)=bvt_fun
         bvzfun(istep,:)=bvz_fun         
         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psis(istep),theta(itheta),1)
            eqfuns(istep,itheta)=eqfun%f(1)
            eqfunx(istep,itheta)=eqfun%fx(1)
            eqfuny(istep,itheta)=eqfun%fy(1)
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
c         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
c     $        mthsurf,eulbpar_mn(istep,:))
c         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
c     $        mthsurf,lagbpar_mn(istep,:))
c         eulbpar=eulbpar_mn(istep,:)
c         lagbpar=lagbpar_mn(istep,:)
c         CALL ipeq_hatoco(psis(istep),eulbpar,mfac,mpert,1,0)
c         CALL ipeq_hatoco(psis(istep),lagbpar,mfac,mpert,1,0)
c         CALL iscdftb(mfac,mpert,eulbparfun(istep,:),mthsurf,eulbpar)
c         CALL iscdftb(mfac,mpert,lagbparfun(istep,:),mthsurf,lagbpar)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psis(istep),theta(itheta),0)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta) 
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write data.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl, FMT='(I1)')labl 
      CALL ascii_open(out_unit,"ipdiag_pmodbrz_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_PMODBRZ: "//
     $     "perturbed mod b on a poloidal plane"
      WRITE(out_unit,'(2x,a8,2x,I4)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I4)')"mthsurf:",mthsurf
      WRITE(out_unit,'(17(2x,a12))')"r","z",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $     "real(xwp)","imag(xwp)","real(xwt)","imag(xwt)",
     $     "real(bvt)","imag(bvt)","real(bvz)","imag(bvz)",
     $     "eqfunx","eqfuny","eqfuns"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(17(2x,e12.3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(eulbparfun(istep,itheta)),
     $           AIMAG(eulbparfun(istep,itheta)),
     $           REAL(lagbparfun(istep,itheta)),
     $           AIMAG(lagbparfun(istep,itheta)),
     $           REAL(xwpfun(istep,itheta)),
     $           AIMAG(xwpfun(istep,itheta)),
     $           REAL(xwtfun(istep,itheta)),
     $           AIMAG(xwtfun(istep,itheta)),
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
      DEALLOCATE(psis,rs,zs,eulbpar_mn,lagbpar_mn,
     $     eulbparfun,lagbparfun,eqfunx,eqfuny,eqfuns,
     $     xwpfun,xwtfun,bvtfun,bvzfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_pmodbrz
c-----------------------------------------------------------------------
c     subprogram 17. ipdiag_rzphibx.
c     write r,z,phi,b,x on hamada coordinates.
c     __________________________________________________________________
c     egnum  : label of an eigenmode
c     xwpimn : xwp_mn components on the control surface
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_rzphibx(egnum,xwpimn,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,labl
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xwpimn

      INTEGER :: i,istep,itheta,ipert,rstep

      REAL(r8) :: bdist,limfac
      REAL(r8), DIMENSION(0:mthsurf) :: t11,t12,t21,t22,t33
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun
      CHARACTER(1) :: slabl

      REAL(r8), DIMENSION(:), POINTER :: psis
      REAL(r8), DIMENSION(:,:), POINTER :: rs,zs,phs
      COMPLEX(r8), DIMENSION(:,:), POINTER :: brs,bzs,bps,xrs,xzs,xps
c-----------------------------------------------------------------------
c     adjust parameters temporarily, bdist is important.
c-----------------------------------------------------------------------
      rstep=200
      bdist=0.1
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
            IF (ABS(singfac(ipert))<bdist) THEN
               limfac = singfac(ipert)/SIGN(bdist,singfac(ipert))
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

         xrs(istep,:)=t11*xwp_fun+t12*xwt_fun
         brs(istep,:)=t11*bwp_fun+t12*bwt_fun
         xzs(istep,:)=t21*xwp_fun+t22*xwt_fun
         bzs(istep,:)=t21*bwp_fun+t22*bwt_fun
         xps(istep,:)=t33*xvz_fun
         bps(istep,:)=t33*bvz_fun
      ENDDO
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      WRITE(UNIT=slabl, FMT='(I1)')labl   
      CALL ascii_open(out_unit,"ipdiag_rzphibx_l"//slabl//"_n"
     $     //sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RZPHIBX: "//
     $     "write down rzphibx informations in hamada coordinates"
      WRITE(out_unit,'(2x,a8,2x,I6)')"rstep:",rstep
      WRITE(out_unit,'(2x,a8,2x,I6)')"mthsurf:",mthsurf
      WRITE(out_unit,'(17(2x,a12))')"psi","theta","r","z","phi",
     $     "real(br)","imag(br)","real(bz)","imag(bz)",
     $     "real(bp)","imag(bp)",
     $     "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $     "real(xp)","imag(xp)"
      DO istep=1,rstep
         DO itheta=0,mthsurf
            WRITE(out_unit,'(17(2x,e12.3))')
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
c     subprogram 18. ipdiag_rzpgrid.
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
     $     "hamada coordinates as functions of rzphi"
      WRITE(out_unit,'(6(2x,a12))')"r","z","limit","psi","theta","phi"

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
            
            WRITE(out_unit,'(6(2x,e12.3))')gdr(i,j),gdz(i,j),gdl(i,j),
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
c     subprogram 19. ipdiag_rzpdiv.
c     check divergence of rzphi functions.
c-----------------------------------------------------------------------
      SUBROUTINE ipdiag_rzpdiv(nr,nz,lval,rval,zval,fr,fz,fp,labl)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nr,nz,labl
      INTEGER, DIMENSION(0:nr,0:nz), INTENT(IN) :: lval
      REAL(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: rval,zval
      COMPLEX(r8), DIMENSION(0:nr,0:nz), INTENT(IN) :: fr,fz,fp
 
      INTEGER :: i,j
      CHARACTER(1) :: slabl
      CHARACTER(2) :: slabl2
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

      IF (labl < 10) THEN
         WRITE(UNIT=slabl, FMT='(I1)')labl    
         CALL ascii_open(out_unit,"ipdiag_rzpdiv_l"//slabl//"_n"
     $        //sn//".out","UNKNOWN")
      ELSE
         WRITE(UNIT=slabl2, FMT='(I2)')labl    
         CALL ascii_open(out_unit,"ipdiag_rzpdiv_l"//slabl2//"_n"
     $        //sn//".out","UNKNOWN")
      ENDIF   
      
      WRITE(out_unit,*)"IPEC_RZPDIV: divergence in rzphi grid"
      WRITE(out_unit,'(1x,a2,5(a16))')"l","r","z",
     $     "re(div)","im(div)","abs(div)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,5(e16.8))')
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
c     subprogram 20. ipdiag_radvar.
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
      qs%xs = psifac
      qs%fs(:,1) = qfac
      CALL spline_fit(qs,"extrap") 
      CALL spline_int(qs)
      qintb = qs%fsi(mstep,1)
      psitor(:) = qs%fsi(:,1)/qintb
      rhotor(:) = SQRT(psitor(:))
      CALL spline_dealloc(qs)
      CALL ascii_open(out_unit,"ipdiag_radvar_n"//sn//".out","UNKNOWN")
      WRITE(out_unit,*)"IPDIAG_RADVAR: "//
     $     "various radial variables"
      WRITE(out_unit,'(2x,a8,2x,I6)')"mstep:",mstep
      WRITE(out_unit,'(2x,a8,2x,e16.9)')"qintb:",qintb
      WRITE(out_unit,'(4(2x,a16))')"psi","rho","psitor","rhotor"
      DO istep=0,mstep
         WRITE(out_unit,'(4(2x,e16.9))')psifac(istep),rhofac(istep),
     $        psitor(istep),rhotor(istep)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipdiag_radvar

      END MODULE ipdiag_mod
