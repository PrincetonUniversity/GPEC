c-----------------------------------------------------------------------
c     file equil_out.f.
c     writes diagnostic output for equilibria.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. equil_out_mod
c     1. equil_out_diagnose.
c     2. equil_out_write_2d.
c     3. equil_out_global.
c     4. equil_out_qfind.
c     5. equil_out_sep_find.
c     6. equil_out_gse.
c     7. equil_out_dump.
c-----------------------------------------------------------------------
c     subprogram 0. equil_out_mod
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE equil_out_mod
      USE global_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. equil_out_diagnose.
c     diagnoses equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_diagnose(sq_flag,unit)

      LOGICAL, INTENT(IN) :: sq_flag
      INTEGER, INTENT(IN) :: unit

      INTEGER :: ipsi
      TYPE(spline_type) :: sq_out
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(//2x,"mpsi",1x,"mtheta",3x,"psilow",4x,"psihigh",
     $     4x,"rs_right",3x,"rs_left",5x,"zs_top",5x,"zs_bot"
     $     //2i6,1p,6e11.3/)
 20   FORMAT(/4x,"amean",6x,"rmean",6x,"aratio",5x,"kappa",6x,"delta1",
     $     5x,"delta2",5x,"volume"//1p,7e11.3/)
 30   FORMAT(/5x,"li1",8x,"li2",8x,"li3",9x,"ro",9x,"zo",8x,"psio"
     $     //1p,6e11.3/)
 40   FORMAT(/4x,"betap1",5x,"betap2",5x,"betap3",5x,"betat",6x,
     $     "betan",5x,"ppeakfac",5x,"bt0"//1p,7e11.3/)
 50   FORMAT(/6x,"q0",8x,"qmin",8x,"q95",7x,"qmax",8x,"qa",8x,"crnt",
     $     7x,"I/aB"//1p,7e11.3/)
 60   FORMAT(/'ipsi',2x,'psifac',6x,'f',7x,'mu0 p',5x,'dvdpsi',6x,'q'/)
 70   FORMAT(i3,1p,5e10.3)
 80   FORMAT(/5x,"i",5x,"m",6x,"q",7x,"dq/dpsi",6x,"psi",8x,"rho"/)
 90   FORMAT(2i6,1p,4e11.3)
 100  FORMAT(/"ipsi",3x,"psifac",7x,"f",8x,"mu0 p",8x,"q"/)
 110  FORMAT(i3,1p,4e11.3)
c-----------------------------------------------------------------------
c     diagnose global data.
c-----------------------------------------------------------------------
      WRITE(unit,10)mpsi,mtheta,psilow,psihigh,rsep,zsep
      WRITE(unit,20)amean,rmean,aratio,kappa,delta1,delta2,volume
      WRITE(unit,30)li1,li2,li3,ro,zo,psio
      WRITE(unit,40)betap1,betap2,betap3,betat,betan,ppeakfac,bt0
      WRITE(unit,50)q0,qmin,q95,qmax,qa,crnt,crnt/(amean*bt0)
c-----------------------------------------------------------------------
c     output surface quantities.
c-----------------------------------------------------------------------
      IF(sq_flag)THEN
         CALL bin_open(bin_unit,"prof.bin","UNKNOWN","REWIND","none")
         WRITE(unit,60)
         DO ipsi=0,mpsi
            WRITE(unit,70)ipsi,sq%xs(ipsi),sq%fs(ipsi,1)/twopi,
     $           sq%fs(ipsi,2),sq%fs(ipsi,3),sq%fs(ipsi,4)
            WRITE(bin_unit)
     $           REAL(sq%xs(ipsi),4),
     $           REAL(SQRT(sq%xs(ipsi)),4),
     $           REAL(sq%fs(ipsi,1)/twopi,4),
     $           REAL(sq%fs(ipsi,2),4),
     $           REAL(sq%fs(ipsi,4),4)
         ENDDO
         WRITE(unit,60)
         WRITE(bin_unit)
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose sq_out.
c-----------------------------------------------------------------------
      IF(out_eq_1d .OR. bin_eq_1d)THEN
         CALL spline_alloc(sq_out,mpsi,4)
         sq_out%xs=sq%xs
         sq_out%fs(:,1)=sq%fs(:,1)/twopi
         sq_out%fs(:,2)=sq%fs(:,2)
         sq_out%fs(:,3)=sq%fs(:,4)
         sq_out%fs(:,4)=SQRT(sq_out%xs)
         CALL spline_fit(sq_out,"extrap")
         sq_out%name="  sq  "
         sq_out%title=(/"psifac","  f   ","mu0 p ","  q   "," rho  "/)
         IF(bin_eq_1d)CALL bin_open(bin_2d_unit,"sq_out.bin","UNKNOWN",
     $        "REWIND","none")
         IF(out_eq_1d)WRITE(out_2d_unit,'(a)')
     $        "output surface quantities:"
         CALL spline_write1(sq_out,out_eq_1d,bin_eq_1d,
     $        out_2d_unit,bin_2d_unit,interp)
         IF(bin_eq_1d)CALL bin_close(bin_2d_unit)
         CALL spline_dealloc(sq_out)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose Grad-Shafranov error.
c-----------------------------------------------------------------------
      IF(gse_flag)CALL equil_out_gse
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_diagnose
c-----------------------------------------------------------------------
c     subprogram 2. equil_out_write_2d.
c     produces ascii and binary output for R,Z(tau,a).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_write_2d

      INTEGER :: ia,itau,ja,jtau
      REAL(r8) :: a,da,dtau,tau,tau1,r,z,deta,eta,rfac,r2,dphi,rho,jac

      REAL(r8), DIMENSION(:), POINTER :: xs,ys
      REAL(r8), DIMENSION(:,:,:), POINTER :: fs
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"ia = ",i3,", a = ",1p,e11.3)
 20   FORMAT(1x,"ia = ",i3,", ja = ",i1,", a = ",1p,e11.3)
 30   FORMAT(/4x,"it",5x,"tau",9x,"r2",8x,"deta",7x,"eta",8x,"dphi",
     $     8x,"r",10x,"z",9x,"jac"/)
 40   FORMAT(i6,1p,8e11.3)
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(.NOT. (out_2d .OR. bin_2d))RETURN
      IF(out_2d)CALL ascii_open(out_2d_unit,"2d.out","UNKNOWN")
      IF(bin_2d)CALL bin_open(bin_2d_unit,"2d.bin","UNKNOWN",
     $     "REWIND","none")
      xs => rzphi%xs
      ys => rzphi%ys
      fs => rzphi%fs
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out_2d)WRITE(out_2d_unit,'(1x,a/)')"input data"
      DO ia=0,mpsi
         rho=SQRT(rzphi%xs(ia))
         IF(out_2d)THEN
            WRITE(out_2d_unit,10)ia,rzphi%xs(ia)
            WRITE(out_2d_unit,30)
         ENDIF
         DO itau=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ia),rzphi%ys(itau),0)
            tau=rzphi%ys(itau)
            r2=rzphi%f(1)
            deta=rzphi%f(2)
            dphi=rzphi%f(3)
            jac=rzphi%f(4)
            eta=tau+deta
            rfac=SQRT(r2)
            r=ro+rfac*COS(twopi*eta)
            z=zo+rfac*SIN(twopi*eta)
            IF(out_2d)WRITE(out_2d_unit,40)
     $           itau,tau,r2,deta,eta,dphi,r,z,jac
            IF(bin_2d)WRITE(bin_2d_unit)REAL(tau,4),REAL(r2,4),
     $           REAL(deta,4),REAL(eta,4),REAL(dphi,4),
     $           REAL(r,4),REAL(z,4),REAL(jac,4)
         ENDDO
         IF(out_2d)WRITE(out_2d_unit,30)
         IF(bin_2d)WRITE(bin_2d_unit)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops for interpolated data.
c-----------------------------------------------------------------------
      IF(interp)THEN
         IF(out_2d)WRITE(out_2d_unit,'(1x,a/)')"interpolated data"
         DO ia=0,mpsi-1
            da=(rzphi%xs(ia+1)-rzphi%xs(ia))/4
            DO ja=0,4
               a=rzphi%xs(ia)+da*ja
               rho=SQRT(a)
               IF(out_2d)THEN
                  WRITE(out_2d_unit,20)ia,ja,a
                  WRITE(out_2d_unit,30)
               ENDIF
               DO itau=0,mtheta-1
                  dtau=(rzphi%ys(itau+1)-rzphi%ys(itau))/4
                  DO jtau=0,4
                     tau1=dtau*jtau
                     tau=rzphi%ys(itau)+tau1
c-----------------------------------------------------------------------
c     compute and print coordinates.
c-----------------------------------------------------------------------
                     CALL bicube_eval(rzphi,a,tau,0)
                     r2=rzphi%f(1)
                     deta=rzphi%f(2)
                     dphi=rzphi%f(3)
                     jac=rzphi%f(4)
                     rfac=SQRT(r2)
                     eta=tau+deta
                     r=ro+rfac*COS(twopi*eta)
                     z=zo+rfac*SIN(twopi*eta)
                     IF(out_2d)WRITE(out_2d_unit,40)
     $                    itau,tau,r2,deta,eta,dphi,r,z,jac
                     IF(bin_2d)WRITE(bin_2d_unit)
     $                    REAL(tau,4),REAL(r2,4),REAL(deta,4),
     $                    REAL(eta,4),REAL(dphi,4),REAL(r,4),REAL(z,4),
     $                    REAL(jac,4)
c-----------------------------------------------------------------------
c     finish loops for interpolated data.
c-----------------------------------------------------------------------
                  ENDDO
               ENDDO
               IF(out_2d)WRITE(out_2d_unit,30)
               IF(bin_2d)WRITE(bin_2d_unit)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(out_2d)CALL ascii_close(out_2d_unit)
      IF(bin_2d)CALL bin_close(bin_2d_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_write_2d
c-----------------------------------------------------------------------
c     subprogram 3. equil_out_global.
c     computes global equilibrium parameters.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_global
      
      INTEGER :: itheta,ipsi,itmax=20
      REAL(r8) :: bp0,chi1,dpsi,dvsq,jacfac,r,rfac,eta,v21,v22,v33,jac
      TYPE(spline_type) :: gs,hs

      LOGICAL, PARAMETER :: diagnose=.FALSE.
c-----------------------------------------------------------------------
c     find inboard and outboard separatrix positions.
c-----------------------------------------------------------------------
      CALL equil_out_sep_find
c-----------------------------------------------------------------------
c     compute shape parameters.
c-----------------------------------------------------------------------
      rmean=(rsep(1)+rsep(2))/2
      amean=(rsep(1)-rsep(2))/2
      aratio=rmean/amean
      kappa=(zsep(1)-zsep(2))/(rsep(1)-rsep(2))
      delta1=(rmean-rext(1))/amean
      delta2=(rmean-rext(2))/amean
      dpsi=1-rzphi%xs(mpsi)
      bt0=(sq%fs(mpsi,1)+sq%fs1(mpsi,1)*dpsi)/(twopi*rmean)
c-----------------------------------------------------------------------
c     prepare cubics spline for integration.
c-----------------------------------------------------------------------
      CALL spline_alloc(gs,mtheta,2)
      CALL spline_alloc(hs,mpsi,3)
      gs%xs=rzphi%ys
      hs%xs=rzphi%xs
c-----------------------------------------------------------------------
c     compute poloidal integrands at plasma surface.
c-----------------------------------------------------------------------
      DO itheta=0,mtheta
         CALL bicube_eval(rzphi,rzphi%xs(mpsi),rzphi%ys(itheta),1)
         jac=rzphi%f(4)
         chi1=twopi*psio/jac
         jacfac=pi/jac
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(rzphi%ys(itheta)+rzphi%f(2))
         r=ro+rfac*COS(eta)
         v21=jacfac*rzphi%fy(1)/(twopi*rfac)
         v22=jacfac*(1+rzphi%fy(2))*(2*rfac)
         v33=jacfac*twopi*(r/pi)
         dvsq=(v21**2+v22**2)*(v33*jac**2)**2
         gs%fs(itheta,1)=SQRT(dvsq)/(twopi*r)
         gs%fs(itheta,2)=chi1*dvsq/(twopi*r)**2
      ENDDO
c-----------------------------------------------------------------------
c     compute plasma current and average surface poloidal field.
c-----------------------------------------------------------------------
      CALL spline_fit(gs,"periodic")
      CALL spline_int(gs)
      crnt=gs%fsi(mtheta,2)/(1e6*mu0)
      bp0=gs%fsi(mtheta,2)/gs%fsi(mtheta,1)
      CALL spline_dealloc(gs)
c-----------------------------------------------------------------------
c     integrate poloidal magnetic field over each flux surface.
c-----------------------------------------------------------------------
      CALL spline_alloc(gs,mtheta,1)
      gs%xs=rzphi%ys
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            jac=rzphi%f(4)
            jacfac=pi/jac
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(rzphi%ys(itheta)+rzphi%f(2))
            r=ro+rfac*COS(eta)
            v21=jacfac*rzphi%fy(1)/(twopi*rfac)
            v22=jacfac*(1+rzphi%fy(2))*(2*rfac)
            v33=jacfac*twopi*(r/pi)
            dvsq=(v21**2+v22**2)*(v33*jac**2)**2
            gs%fs(itheta,1)=dvsq/(r*r)/jac
         ENDDO
         CALL spline_fit(gs,"periodic")
         CALL spline_int(gs)
         hs%fs(ipsi,3)=gs%fsi(mtheta,1)*psio**2
      ENDDO
c-----------------------------------------------------------------------
c     compute flux functions and fit to cubic splines.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         hs%fs(ipsi,1)=sq%fs(ipsi,2)*sq%fs(ipsi,3)
         hs%fs(ipsi,2)=sq%fs(ipsi,3)
      ENDDO
      hs%xpower(1,3)=-1
      CALL spline_fit(hs,"extrap")
c-----------------------------------------------------------------------
c     diagnose flux functions.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(debug_unit,"hs.out","UNKNOWN")
         CALL bin_open(bin_unit,"hs.bin","UNKNOWN","REWIND","none")
         hs%title=(/" psi  "," hs1  "," hs2  "," hs3  "/)
         CALL spline_write1(hs,.TRUE.,.TRUE.,debug_unit,bin_unit,.TRUE.)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)
      ENDIF
c-----------------------------------------------------------------------
c     integrate flux functions.
c-----------------------------------------------------------------------
      CALL spline_int(hs)
      hs%fsi(mpsi,:)=hs%fsi(mpsi,:)
     $     +(hs%fs(mpsi,:)+hs%fs1(mpsi,:)*dpsi/2)*dpsi
      CALL spline_int(sq)
      volume=sq%fsi(mpsi,3)+(sq%fs(mpsi,3)+sq%fs1(mpsi,3)*dpsi/2)*dpsi
c-----------------------------------------------------------------------
c     compute remaining global quantities from integrals.
c-----------------------------------------------------------------------
      p0=sq%fs(0,2)-sq%fs1(0,2)*sq%xs(0)
      ppeakfac=p0*hs%fsi(mpsi,2)/hs%fsi(mpsi,1)
      betat=2*(hs%fsi(mpsi,1)/hs%fsi(mpsi,2))/bt0**2
      betan=100*amean*bt0*betat/crnt
      betap1=2*(hs%fsi(mpsi,1)/hs%fsi(mpsi,2))/bp0**2
      betap2=4*hs%fsi(mpsi,1)/((1e6*mu0*crnt)**2*ro)
      betap3=4*hs%fsi(mpsi,1)/((1e6*mu0*crnt)**2*rmean)
      li1=hs%fsi(mpsi,3)/hs%fsi(mpsi,2)/bp0**2
      li2=2*hs%fsi(mpsi,3)/((1e6*mu0*crnt)**2*ro)
      li3=2*hs%fsi(mpsi,3)/((1e6*mu0*crnt)**2*rmean)
c-----------------------------------------------------------------------
c     deallocate spline_types.
c-----------------------------------------------------------------------
      CALL spline_dealloc(gs)
      CALL spline_dealloc(hs)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_global
c-----------------------------------------------------------------------
c     subprogram 4. equil_out_qfind.
c     finds special values of q.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_qfind

      INTEGER :: ipsi,iex,i
      INTEGER, PARAMETER :: itmax=200
      REAL(r8) :: a,b,c,d,dx,x0,x,xmax
      REAL(r8), DIMENSION(0:100) :: psiexl,qexl
c-----------------------------------------------------------------------
c     store left end point.
c-----------------------------------------------------------------------
      iex=0
      psiexl(iex)=sq%xs(0)
      qexl(iex)=sq%fs(0,4)
c-----------------------------------------------------------------------
c     set up search for extrema of q.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi-1
         CALL spline_eval(sq,sq%xs(ipsi),3)
         xmax=sq%xs(ipsi+1)-sq%xs(ipsi)
         a=sq%f(4)
         b=sq%f1(4)
         c=sq%f2(4)
         d=sq%f3(4)
         x0=-c/d
         dx=x0*x0-2*b/d
c-----------------------------------------------------------------------
c     store extremum.
c-----------------------------------------------------------------------
         IF(dx >= 0)THEN
            dx=SQRT(dx)
            DO i=1,2
               x=x0-dx
               IF(x >= 0 .AND. x < xmax)THEN
                  iex=iex+1
                  psiexl(iex)=sq%xs(ipsi)+x
                  CALL spline_eval(sq,psiexl(iex),0)
                  qexl(iex)=sq%f(4)
               ENDIF
               dx=-dx
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c     complete search and store right end point.
c-----------------------------------------------------------------------
      ENDDO
      iex=iex+1
      psiexl(iex)=sq%xs(mpsi)
      qexl(iex)=sq%fs(mpsi,4)
      mex=iex
c-----------------------------------------------------------------------
c     copy to permanent storage
c-----------------------------------------------------------------------
      ALLOCATE(psiex(0:mex),qex(0:mex))
      psiex(0:mex)=psiexl(0:mex)
      qex(0:mex)=qexl(0:mex)
c-----------------------------------------------------------------------
c     find special q values.
c-----------------------------------------------------------------------
      q0=sq%fs(0,4)-sq%fs1(0,4)*sq%xs(0)
      qmax=sq%fs(mpsi,4)
      qmin=MIN(MINVAL(qex(0:mex)),q0)
      qmax=MAX(MAXVAL(qex(0:mex)),qmax)
      qa=sq%fs(mpsi,4)+sq%fs1(mpsi,4)*(1-sq%xs(mpsi))
      CALL spline_eval(sq,0.95_r8,0)
      q95=sq%f(4)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_qfind
c-----------------------------------------------------------------------
c     subprogram 5. equil_out_sep_find.
c     finds separatrix locations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_sep_find

      CHARACTER(64) :: message
      INTEGER :: iside,it,itmax=100
      INTEGER, DIMENSION(1) :: index
      REAL(r8) :: cosfac,dtheta,eta,eta0,eta1,eta2,eta_theta,
     $     phase,phase1,phase2,psifac,r2,r21,r22,rfac,rfac1,rfac2,
     $     sinfac,theta,z,z1,z2
      REAL(r8), DIMENSION(0:mtheta) :: vector
c-----------------------------------------------------------------------
c     start loop over sides.
c-----------------------------------------------------------------------
      DO it=0,rzphi%my
         CALL bicube_eval(rzphi,rzphi%xs(mpsi),rzphi%ys(it),0)
         vector(it)=rzphi%ys(it)+rzphi%f(2)
      ENDDO
      psifac=rzphi%xs(mpsi)
      eta0=0
      index=MINLOC(ABS(vector-eta0))
      theta=rzphi%ys(index(1))
      DO iside=1,2
c-----------------------------------------------------------------------
c     newton iteration.
c-----------------------------------------------------------------------
         it=0
         DO
            it=it+1
            CALL bicube_eval(rzphi,psifac,theta,1)
            eta=theta+rzphi%f(2)-eta0
            eta_theta=1+rzphi%fy(2)
            dtheta=-eta/eta_theta
            theta=theta+dtheta
            IF(ABS(eta) <= 1e-10)EXIT
            IF(it > itmax)THEN
               WRITE(message,'(a,i1)')
     $              "Can't find separatrix position, iside = ",iside
               CALL program_stop(message)
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     compute major radius of separatrix.
c-----------------------------------------------------------------------
         CALL bicube_eval(rzphi,psifac,theta,0)
         rsep(iside)=ro+SQRT(rzphi%f(1))*COS(twopi*(theta+rzphi%f(2)))
         eta0=.5
         index=MINLOC(ABS(vector-eta0))
         theta=rzphi%ys(index(1))
      ENDDO
c-----------------------------------------------------------------------
c     find top and bottom.
c-----------------------------------------------------------------------
      DO it=0,rzphi%my
         CALL bicube_eval(rzphi,rzphi%xs(mpsi),rzphi%ys(it),0)
         vector(it)=SQRT(rzphi%f(1))
     $        *SIN(twopi*(rzphi%ys(it)+rzphi%f(2)))
      ENDDO
      index=MAXLOC(vector)
      theta=rzphi%ys(index(1)-1)
      psifac=rzphi%xs(mpsi)
      DO iside=1,2
         DO
            CALL bicube_eval(rzphi,psifac,theta,2)
            r2=rzphi%f(1)
            r21=rzphi%fy(1)
            r22=rzphi%fyy(1)
            eta=rzphi%f(2)
            eta1=rzphi%fy(2)
            eta2=rzphi%fyy(2)
            rfac=SQRT(r2)
            rfac1=r21/(2*rfac)
            rfac2=(r22-r21*rfac1/rfac)/(2*rfac)
            phase=twopi*(theta+eta)
            phase1=twopi*(1+eta1)
            phase2=twopi*eta2
            cosfac=COS(phase)
            sinfac=SIN(phase)
            z=zo+rfac*sinfac
            z1=rfac*phase1*cosfac+rfac1*sinfac
            z2=(2*rfac1*phase1+rfac*phase2)*cosfac
     $           +(rfac2-rfac*phase1**2)*sinfac
            dtheta=-z1/z2
            theta=theta+dtheta
            IF(ABS(dtheta) < 1E-12*ABS(theta))EXIT
         ENDDO
         rext(iside)=ro+rfac*cosfac
         zsep(iside)=z
         index=MINLOC(vector)
         theta=rzphi%ys(index(1)-1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_sep_find
c-----------------------------------------------------------------------
c     subprogram 6. equil_out_gse.
c     diagnoses grad-shafranov solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_gse

      INTEGER :: ipsi,itheta,iqty
      REAL(r8), DIMENSION(0:mpsi,2) :: term
      REAL(r8), DIMENSION(0:mpsi) :: totali,errori,errlogi
      REAL(r8), DIMENSION(0:mtheta) :: rfac,angle
      REAL(r8), DIMENSION(0:mpsi,0:mtheta) :: r,z
      REAL(r8), DIMENSION(0:mpsi,0:mtheta) :: source,total,error,errlog
      TYPE(bicube_type) :: flux
      TYPE(spline_type) :: temp
c-----------------------------------------------------------------------
c     compute coordinates.
c-----------------------------------------------------------------------
      WRITE(*,*)"Diagnosing Grad-Shafranov solution"
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),0)
            rfac(itheta)=SQRT(rzphi%f(1))
            angle(itheta)=twopi*(rzphi%ys(itheta)+rzphi%f(2))
         ENDDO
         r(ipsi,:)=ro+rfac*COS(angle)
         z(ipsi,:)=zo+rfac*SIN(angle)
      ENDDO
c-----------------------------------------------------------------------
c     allocate and compute fluxes.
c-----------------------------------------------------------------------
      CALL bicube_alloc(flux,mpsi,mtheta,2)
      flux%xs=rzphi%xs
      flux%ys=rzphi%ys
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            flux%fs(ipsi,itheta,1)=rzphi%fy(1)**2/(twopi**2*rzphi%f(1))
     $           +(1+rzphi%fy(2))**2*4*rzphi%f(1)
            flux%fs(ipsi,itheta,2)=rzphi%fx(1)*rzphi%fy(1)
     $           /(twopi**2*rzphi%f(1))+rzphi%fx(2)
     $           *(1+rzphi%fy(2))*4*rzphi%f(1)
            DO iqty=1,2
               flux%fs(ipsi,itheta,iqty)=flux%fs(ipsi,itheta,iqty)*twopi
     $              *psio/rzphi%f(4)
            ENDDO
         ENDDO
      ENDDO
      CALL bicube_fit(flux,"extrap","periodic")
c-----------------------------------------------------------------------
c     compute source.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),0)
            source(ipsi,itheta)=rzphi%f(4)/(twopi*psio*pi**2)
     $           *(sq%fs(ipsi,1)*sq%fs1(ipsi,1)
     $           /(twopi*r(ipsi,itheta))**2+sq%fs1(ipsi,2))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute total and error.
c-----------------------------------------------------------------------
      total=flux%fsx(:,:,1)-flux%fsy(:,:,2)+source
      error=ABS(total)/MAX(MAXVAL(ABS(flux%fsx(:,:,1))),
     $     MAXVAL(ABS(flux%fsy(:,:,2))),MAXVAL(ABS(source)))
      WHERE(error > 0)
         errlog=LOG10(error)
      ELSEWHERE
         errlog=0
      ENDWHERE
c-----------------------------------------------------------------------
c     write contour plot.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"gsec.bin","UNKNOWN","REWIND","none")
      WRITE(bin_unit)1,0
      WRITE(bin_unit)mpsi,mtheta
      WRITE(bin_unit)REAL(r,4),REAL(z,4)
      WRITE(bin_unit)REAL(flux%fsx(:,:,1),4)
      WRITE(bin_unit)REAL(flux%fsy(:,:,2),4)
      WRITE(bin_unit)REAL(source,4)
      WRITE(bin_unit)REAL(total,4)
      WRITE(bin_unit)REAL(error,4)
      WRITE(bin_unit)REAL(errlog,4)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     write xy plot.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"gse.bin","UNKNOWN","REWIND","none")
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            WRITE(bin_unit)
     $           REAL(flux%ys(itheta),4),
     $           REAL(flux%xs(ipsi),4),
     $           REAL(flux%fs(ipsi,itheta,:),4),
     $           REAL(source(ipsi,itheta),4),
     $           REAL(total(ipsi,itheta),4),
     $           REAL(error(ipsi,itheta),4),
     $           REAL(errlog(ipsi,itheta),4)
         ENDDO
         WRITE(bin_unit)
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     compute integrated error criterion.
c-----------------------------------------------------------------------
      CALL spline_alloc(temp,mtheta,2)
      temp%xs=flux%ys
      DO ipsi=0,mpsi
         temp%fs(:,1)=flux%fsx(ipsi,:,1)
         temp%fs(:,2)=source(ipsi,:)
         CALL spline_fit(temp,"periodic")
         CALL spline_int(temp)
         term(ipsi,:)=temp%fsi(mtheta,:)
      ENDDO
      CALL spline_dealloc(temp)
c-----------------------------------------------------------------------
c     compute total and error.
c-----------------------------------------------------------------------
      totali=SUM(term,2)
      errori=ABS(totali)
      WHERE(errori > 0)
         errlogi=LOG10(errori)
      ELSEWHERE
         errlogi=0
      ENDWHERE
c-----------------------------------------------------------------------
c     write integrated error criterion.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"gsei.bin","UNKNOWN","REWIND","none")
      DO ipsi=0,mpsi
         WRITE(bin_unit)
     $        REAL(flux%xs(ipsi),4),
     $        REAL(term(ipsi,:),4),
     $        REAL(totali(ipsi),4),
     $        REAL(errori(ipsi),4),
     $        REAL(errlogi(ipsi),4)
      ENDDO
      WRITE(bin_unit)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     deallocate bicube_type.
c-----------------------------------------------------------------------
      CALL bicube_dealloc(flux)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_gse
c-----------------------------------------------------------------------
c     subprogram 7. equil_out_dump.
c     writes dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_out_dump

      CHARACTER(133) :: filename

      LOGICAL, PARAMETER :: diagnose=.TRUE.
      INTEGER :: ix,iy,iqty
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(2x,"nqty",4x,"mx",4x,"my",
     $     5x,"x01",8x,"x02",8x,"y01",8x,"y02"//3i6,1p,4e11.3/)
 20   FORMAT(/2x,"iqty",5x,"xpr1",7x,"xpr2",7x,"ypr1",7x,"ypr2"/)
 30   FORMAT(i6,1p,4e11.3)
 40   FORMAT(/4x,"ix",6x,"xs"/)
 50   FORMAT(i6,1p,e11.3)
 60   FORMAT(/4x,"iy",6x,"ys"/)
 70   FORMAT(/4x,"iy",4x,"ix",4(5x,"fs",i1,3x)/)
 80   FORMAT(2i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     open dump file.
c-----------------------------------------------------------------------
      filename=TRIM(eq_filename)//".dump"
      CALL bin_open(dump_unit,TRIM(filename),"UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     write scalars.
c-----------------------------------------------------------------------
      WRITE(dump_unit)mpsi,mtheta,jac_type,power_bp,power_b,power_r,
     $     grid_type,psilow,psihigh,ro,zo,psio,q0,qa
c-----------------------------------------------------------------------
c     write 1D splines.
c-----------------------------------------------------------------------
      WRITE(dump_unit)sq%mx,sq%nqty
      WRITE(dump_unit)sq%xs,sq%fs,
     $     sq%xpower,sq%x0,sq%title,sq%name,sq%periodic
c-----------------------------------------------------------------------
c     write 2D splines.
c-----------------------------------------------------------------------
      WRITE(dump_unit)rzphi%mx,rzphi%my,rzphi%nqty
      WRITE(dump_unit)rzphi%xs
      WRITE(dump_unit)rzphi%ys
      WRITE(dump_unit)rzphi%fs
      WRITE(dump_unit)rzphi%x0
      WRITE(dump_unit)rzphi%y0
      WRITE(dump_unit)rzphi%xpower
      WRITE(dump_unit)rzphi%ypower
      WRITE(dump_unit)rzphi%xtitle
      WRITE(dump_unit)rzphi%ytitle
      WRITE(dump_unit)rzphi%title
      WRITE(dump_unit)rzphi%name
      WRITE(dump_unit)rzphi%periodic
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL bin_close(dump_unit)
c-----------------------------------------------------------------------
c     diagnose rzphi.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(dump_unit,"rzphi0.out","UNKNOWN")
         WRITE(dump_unit,10)rzphi%nqty,rzphi%mx,rzphi%my,
     $        rzphi%x0,rzphi%y0
         WRITE(dump_unit,20)
         WRITE(dump_unit,30)(iqty,rzphi%xpower(:,iqty),
     $        rzphi%ypower(:,iqty),iqty=1,rzphi%nqty)
         WRITE(dump_unit,20)
         WRITE(dump_unit,40)
         WRITE(dump_unit,50)(ix,rzphi%xs(ix),ix=0,rzphi%mx)
         WRITE(dump_unit,40)
         WRITE(dump_unit,60)
         WRITE(dump_unit,50)(iy,rzphi%ys(iy),iy=0,rzphi%my)
         WRITE(dump_unit,60)
         WRITE(dump_unit,70)(iqty,iqty=1,rzphi%nqty)
         DO iy=0,rzphi%my
            DO ix=0,rzphi%mx
               WRITE(dump_unit,80)iy,ix,rzphi%fs(ix,iy,:)
            ENDDO
            WRITE(dump_unit,70)(iqty,iqty=1,rzphi%nqty)
         ENDDO
         CALL bin_close(dump_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_out_dump
      END MODULE equil_out_mod
