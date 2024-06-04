c-----------------------------------------------------------------------
c     file direct.f.
c     processes direct equilibria.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c----------------------------------------------------------------------
c     0. direct_mod.
c     1. direct_run.
c     2. direct_get_bfield.
c     3. direct_position.
c     4. direct_fl_int.
c     4.5 find_fl_surface.
c     5. direct_fl_der.
c     6. direct_refine.
c     7. direct_output.
c     8. direct_local_xpoint.
c     9. direct_saddle_angle.
c     10. direct_psisaddle.
c     11. direct_xpoint.
c     12. direct_initialise_xpoints.
c     17. direct_Blocal
c-----------------------------------------------------------------------
c     subprogram 0. direct_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE direct_mod
      USE global_mod
      USE utils_mod
      IMPLICIT NONE

      INTEGER, PRIVATE :: istep
      REAL(r8) :: rmin,rmax,zmin,zmax,rs1,rs2
      REAL(r8), DIMENSION(1:2) :: xpt_etas, rxs, zxs, xpt_b11s 
      REAL(r8), DIMENSION(1:2) :: xpt_gammas, xpt_varthetas
      REAL(r8), DIMENSION(2,2) :: xpt_brackets
      TYPE(bicube_type) :: psi_in
      LOGICAL :: direct_infinite_loop_flag
      INTEGER :: direct_infinite_loop_count = 2000
      INTEGER :: num_xpts

      TYPE :: direct_bfield_type
      REAL(r8) :: psi,psir,psiz,psirz,psirr,psizz,f,f1,p,p1
      REAL(r8) :: br,bz,brr,brz,bzr,bzz
      END TYPE direct_bfield_type

      REAL(r8) :: etol=1e-8
      INTEGER :: nstepd=2048

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. direct_run.
c     gets equilibrium data and massages it.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_run

      INTEGER :: ir,iz,itheta,ipsi,len_y_out,len_y_last
      INTEGER :: maxima_count,i
      REAL(r8) :: f0fac,f0,ffac,rfac,eta,r,jacfac,w11,w12,delpsi,q,flast
      REAL(r8), DIMENSION(0:nstepd,0:4) :: y_out, y_out_last
      REAL(r8), DIMENSION(2, mpsi+1) :: xdx
      REAL(r8), DIMENSION(3,3) :: v
      REAL(r8), DIMENSION(3,2) :: eta_brackets
      REAL(r8), DIMENSION(3) :: eta_maxes

      LOGICAL :: use_analytic,run_xpt

      REAL(r8) :: xm,dx,rholow,rhohigh,rx,zx
      TYPE(direct_bfield_type) :: bf
      TYPE(spline_type) :: ff

      use_analytic=.FALSE.
      run_xpt=.TRUE.
      xpt_etas=0.0
      eta_maxes=0.0
      eta_brackets=0.0 
      num_xpts=0
      maxima_count=0
c-----------------------------------------------------------------------
c     warning.
c-----------------------------------------------------------------------
      direct_flag=.TRUE.
      IF(psihigh >= 1-1e-6)WRITE(*,'(1x,a,es10.3,a)')
     $        "Warning: direct equilibrium with psihigh =",psihigh,
     $        " could hang on separatrix."
      direct_infinite_loop_flag = .FALSE.
c-----------------------------------------------------------------------
c     fit input to cubic splines and diagnose.
c-----------------------------------------------------------------------
      sq_in%fs(:,4)=SQRT(sq_in%xs)
      sq_in%name="  sq  "
      sq_in%title=(/"psifac","  f   ","mu0 p ","  q   "," rho  "/)
      CALL spline_fit(sq_in,"extrap")
      psi_in%xs=rmin+(/(ir,ir=0,psi_in%mx)/)*(rmax-rmin)/psi_in%mx
      psi_in%ys=zmin+(/(iz,iz=0,psi_in%my)/)*(zmax-zmin)/psi_in%my
      CALL bicube_fit(psi_in,"extrap","extrap")
      CALL direct_output
c-----------------------------------------------------------------------
c     prepare new spline type for surface quantities.
c-----------------------------------------------------------------------
      IF(grid_type == "original" .OR. grid_type == "orig")THEN
         IF(sq_in%xs(sq_in%mx) < 1-1e-6)THEN
            mpsi=sq_in%mx-1
         ELSE
            mpsi=sq_in%mx-2
         ENDIF
      ENDIF
      CALL spline_alloc(sq,mpsi,4)
      sq%name="  sq  "
      sq%title=(/"psifac","twopif","mu0 p ","dvdpsi","  q   "/)
c-----------------------------------------------------------------------
c     set up radial grid
c-----------------------------------------------------------------------
      SELECT CASE(grid_type)
      CASE("ldp")
         sq%xs=(/(ipsi,ipsi=0,mpsi)/)/REAL(mpsi,r8)
         sq%xs=psilow+(psihigh-psilow)*SIN(sq%xs*pi/2)**2
      CASE("pow1")
         xdx = powspace(psilow, psihigh, 1, mpsi+1, "upper")
         sq%xs=xdx(1,:)
      CASE("pow2")
         xdx = powspace(psilow, psihigh, 2, mpsi+1, "upper")
         sq%xs=xdx(1,:)
      CASE("rho")
         sq%xs=psihigh*(/(ipsi**2,ipsi=1,mpsi+1)/)/(mpsi+1)**2
      CASE("original","orig")
         sq%xs=sq_in%xs(1:mpsi)
      CASE("mypack")
         rholow=SQRT(psilow)
         rhohigh=SQRT(psihigh)
         xm=(rhohigh+rholow)/2
         dx=(rhohigh-rholow)/2
         sq%xs=xm+dx*mypack(mpsi/2,sp_pfac,"both")
         sq%xs=sq%xs**2
      CASE default
         CALL program_stop("Cannot recognize grid_type "//grid_type)
      END SELECT
c-----------------------------------------------------------------------
c     define positions and open output files.
c-----------------------------------------------------------------------
      CALL direct_position
      IF(out_fl)CALL ascii_open(out_2d_unit,"flint.out","UNKNOWN")
      IF(bin_fl)CALL bin_open(bin_2d_unit,"flint.bin","UNKNOWN",
     $     "REWIND","none")
c-----------------------------------------------------------------------
c     start loop over flux surfaces.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,'(a,1p,es10.3)')" etol = ",etol
      IF(verbose) WRITE(*,'(a,1p,i6)')" nstepd = ",nstepd
      DO ipsi=0,mpsi,+1
c-----------------------------------------------------------------------
c     logic whether to integrate around whole field line or use analytic
c     integral formulas near separatrix.
c-----------------------------------------------------------------------
         IF(use_analytic .AND. run_xpt)THEN
            CALL direct_mixed_spline_builder(sq%xs(ipsi))
         ELSE
            CALL direct_fl_int(sq%xs(ipsi),zero,twopi,y_out,bf,
     $                                                        len_y_out)
c-----------------------------------------------------------------------
c     checks whether q-integral is diverging.  
c-----------------------------------------------------------------------
            IF(sq%xs(ipsi)>0.999 .AND. run_xpt)THEN
               CALL direct_initialise_xpoints(y_out,len_y_out,.TRUE.,
     $             .TRUE.,bf,10*one,eta_maxes,eta_brackets,maxima_count)

               IF(maxima_count > 0 .AND. run_xpt)THEN
                  use_analytic = .TRUE.
                  num_xpts=maxima_count

                  DO i=1,maxima_count,+1
                     xpt_etas(i)=eta_maxes(i) !updated by direct_xpoint
                     xpt_brackets(i,1)=eta_brackets(i,1)
                     xpt_brackets(i,2)=eta_brackets(i,2)

                     CALL find_fl_surface(one,xpt_etas(i),rx,zx)
                     CALL direct_xpoint(rx,zx,i)
                  ENDDO
                  CYCLE
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c     fit data to cubic splines.
c-----------------------------------------------------------------------
            CALL spline_alloc(ff,istep,4)
            ff%xs(0:istep)=y_out(0:istep,4)/y_out(istep,4)
            ff%fs(0:istep,1)=y_out(0:istep,2)**2
            ff%fs(0:istep,2)=y_out(0:istep,0)/twopi-ff%xs(0:istep)
            ff%fs(0:istep,3)=bf%f*
     $        (y_out(0:istep,3)-ff%xs(0:istep)*y_out(istep,3))
            ff%fs(0:istep,4)=y_out(0:istep,1)/y_out(istep,1)-ff%xs
            CALL spline_fit(ff,"periodic")
c-----------------------------------------------------------------------
c     allocate space for rzphi and define grids.
c-----------------------------------------------------------------------
            IF(ipsi == 0)THEN
               IF(mtheta == 0)mtheta=istep
               CALL bicube_alloc(rzphi,mpsi,mtheta,4) !change mtheta
               CALL bicube_alloc(eqfun,mpsi,mtheta,3) ! new eq information
               rzphi%xs=sq%xs
               rzphi%ys=(/(itheta,itheta=0,mtheta)/)/REAL(mtheta,r8)
               rzphi%xtitle="psifac"
               rzphi%ytitle="theta "
               rzphi%title=(/"  r2  "," deta "," dphi ","  jac "/) ! deta = eta - theta, dphi = phi - zeta
               eqfun%title=(/"  b0  ","      ","      " /)
               eqfun%xs=sq%xs
               eqfun%ys=(/(itheta,itheta=0,mtheta)/)/REAL(mtheta,r8)
            ENDIF
c-----------------------------------------------------------------------
c     interpolate to uniform grid.
c-----------------------------------------------------------------------
            DO itheta=0,mtheta
               CALL spline_eval(ff,rzphi%ys(itheta),1)
               rzphi%fs(ipsi,itheta,1:3)=ff%f(1:3)
               rzphi%fs(ipsi,itheta,4)=(1+ff%f1(4))
     $           *y_out(istep,1)*twopi*psio
            ENDDO
c-----------------------------------------------------------------------
c     store surface quantities.
c-----------------------------------------------------------------------
            sq%fs(ipsi,1)=bf%f*twopi
            sq%fs(ipsi,2)=bf%p
            sq%fs(ipsi,3)=y_out(istep,1)*twopi*psio
            sq%fs(ipsi,4)=y_out(istep,3)*bf%f/twopi
            CALL spline_dealloc(ff)
c-----------------------------------------------------------------------
c     log maximum y_out, flast for x-point plotting.
c-----------------------------------------------------------------------
            IF (ipsi == mpsi) THEN
                  y_out_last = y_out
                  len_y_last=len_y_out
                  flast = bf%f*twopi
            ENDIF
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(out_fl)CALL ascii_close(out_2d_unit)
      IF(bin_fl)CALL bin_close(bin_2d_unit)
c-----------------------------------------------------------------------
c     fit surface quantities to cubic splines.
c-----------------------------------------------------------------------
      CALL spline_fit(sq,"extrap")
      sq%name="  sq  "
      sq%title=(/" psi  ","  f   ","  p   ","  q   "/)
      q0=sq%fs(0,4)-sq%fs1(0,4)*sq%xs(0)
      IF(newq0 == -1)newq0=-q0
c-----------------------------------------------------------------------
c     revise q profile.
c-----------------------------------------------------------------------
      IF(newq0 /= 0)THEN
         f0=sq%fs(0,1)-sq%fs1(0,1)*sq%xs(0)
         f0fac=f0**2*((newq0/q0)**2-1)
         q0=newq0
         DO ipsi=0,mpsi
            ffac=SQRT(1+f0fac/sq%fs(ipsi,1)**2)*SIGN(one,newq0)
            sq%fs(ipsi,1)=sq%fs(ipsi,1)*ffac
            sq%fs(ipsi,4)=sq%fs(ipsi,4)*ffac
            rzphi%fs(ipsi,:,3)=rzphi%fs(ipsi,:,3)*ffac
         ENDDO
         CALL spline_fit(sq,"extrap")
      ENDIF
      qa=sq%fs(mpsi,4)+sq%fs1(mpsi,4)*(1-sq%xs(mpsi))
c-----------------------------------------------------------------------
c     fit rzphi to bicubic splines.
c-----------------------------------------------------------------------
      IF(power_flag)rzphi%xpower(1,:)=(/1._r8,0._r8,.5_r8,0._r8/)
      CALL bicube_fit(rzphi,"extrap","periodic")
      CALL spline_dealloc(sq_in)
c-----------------------------------------------------------------------
c     evaluate eqfun.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         CALL spline_eval(sq,sq%xs(ipsi),0)
         q=sq%f(4)
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(itheta/REAL(mtheta,r8)+rzphi%f(2))
            r=ro+rfac*COS(eta)
            jacfac=rzphi%f(4)
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(1,3)=rzphi%fx(3)*r
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(2,3)=rzphi%fy(3)*r
            v(3,3)=twopi*r
            w11=(1+rzphi%fy(2))*twopi**2*rfac*r/jacfac
            w12=-rzphi%fy(1)*pi*r/(rfac*jacfac)

            delpsi=SQRT(w11**2+w12**2)
            eqfun%fs(ipsi,itheta,1)=SQRT(((twopi*psio*delpsi)**2+
     $           sq%f(1)**2)/(twopi*r)**2)
            eqfun%fs(ipsi,itheta,2)=(SUM(v(1,:)*v(2,:))+q*v(3,3)*v(1,3))
     $           /(jacfac*eqfun%fs(ipsi,itheta,1)**2)
            eqfun%fs(ipsi,itheta,3)=(v(2,3)*v(3,3)+q*v(3,3)*v(3,3))
     $           /(jacfac*eqfun%fs(ipsi,itheta,1)**2)
         ENDDO
      ENDDO
      CALL bicube_fit(eqfun,"extrap","periodic")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_run
c-----------------------------------------------------------------------
c     subprogram 2. direct_get_bfield.
c     evaluates bicubic splines for field components and derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_get_bfield(r,z,bf,mode)

      INTEGER, INTENT(IN) :: mode
      REAL(r8), INTENT(IN) :: r,z
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     compute spline interpolations.
c-----------------------------------------------------------------------
      CALL bicube_eval(psi_in,r,z,mode)
      bf%psi=psi_in%f(1)
      CALL spline_eval(sq_in,1-bf%psi/psio,1)
      bf%f=sq_in%f(1)
      bf%f1=sq_in%f1(1)
      bf%p=sq_in%f(2)
      bf%p1=sq_in%f1(2)
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     evaluate magnetic fields.
c-----------------------------------------------------------------------
      bf%psir=psi_in%fx(1)
      bf%psiz=psi_in%fy(1)
      bf%br=bf%psiz/r
      bf%bz=-bf%psir/r
      IF(mode == 1)RETURN
c-----------------------------------------------------------------------
c     evaluate derivatives of magnetic fields.
c-----------------------------------------------------------------------
      bf%psirr=psi_in%fxx(1)
      bf%psirz=psi_in%fxy(1)
      bf%psizz=psi_in%fyy(1)
      bf%brr=(bf%psirz-bf%br)/r
      bf%brz=bf%psizz/r
      bf%bzr=-(bf%psirr+bf%bz)/r
      bf%bzz=-bf%psirz/r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_get_bfield
c-----------------------------------------------------------------------
c     subprogram 3. direct_position.
c     finds radial positions of o-point and separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_position

      REAL(r8), PARAMETER :: eps=1e-12
      REAL(r8) :: ajac(2,2),det,dr,dz,fac,r,z
      TYPE(direct_bfield_type) :: bf
      INTEGER :: ir,ird
c-----------------------------------------------------------------------
c     scan to find zero crossing of bz on midplane.
c-----------------------------------------------------------------------
      IF(ro == 0)THEN
         r=(rmax+rmin)/2
         z=(zmax+zmin)/2
         dr=(rmax-rmin)/20
         ir = 0
         DO
            CALL direct_get_bfield(r,z,bf,1)
            IF(bf%bz >= 0)EXIT
            r=r+dr

            ir = ir+1
            IF (ir  > direct_infinite_loop_count) THEN
               direct_infinite_loop_flag = .TRUE.
               CALL program_stop("Took too many steps to get bz=0.")
            ENDIF
         ENDDO
      ELSE
         r=ro
         z=zo
      ENDIF
c-----------------------------------------------------------------------
c     use newton iteration to find o-point.
c-----------------------------------------------------------------------
      ir = 0
      DO
         CALL direct_get_bfield(r,z,bf,2)
         ajac(1,1)=bf%brr
         ajac(1,2)=bf%brz
         ajac(2,1)=bf%bzr
         ajac(2,2)=bf%bzz
         det=ajac(1,1)*ajac(2,2)-ajac(1,2)*ajac(2,1)
         dr=(ajac(1,2)*bf%bz-ajac(2,2)*bf%br)/det
         dz=(ajac(2,1)*bf%br-ajac(1,1)*bf%bz)/det
         r=r+dr
         z=z+dz
         IF(ABS(dr) <= eps*r .AND. ABS(dz) <= eps*r)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find o-point.")
         ENDIF
      ENDDO
      ro=r
      zo=z
      ro=r
      zo=z
c-----------------------------------------------------------------------
c     renormalize psi.
c-----------------------------------------------------------------------
      fac=psio/bf%psi
      psi_in%fs=psi_in%fs*fac
      psi_in%fsx=psi_in%fsx*fac
      psi_in%fsy=psi_in%fsy*fac
      psi_in%fsxy=psi_in%fsxy*fac
c-----------------------------------------------------------------------
c     use newton iteration to find inboard separatrix position.
c-----------------------------------------------------------------------
      ird=0
      r=((3-0.5*ird)*rmin+ro)/(1+3-0.5*ird)
      z=zo
      ir=0
      DO
         CALL direct_get_bfield(r,z,bf,1)
         dr=-bf%psi/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r) EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            ird=ird+1
            r=((3-0.5*ird)*rmin+ro)/(1+3-0.5*ird)
            ir=0
            IF (ird==6) THEN 
               direct_infinite_loop_flag = .TRUE.
               CALL program_stop("Took too many steps to find inb spx.")
            ENDIF
         ENDIF
      ENDDO
      rs1=r
c-----------------------------------------------------------------------
c     use newton iteration to find outboard separatrix position.
c-----------------------------------------------------------------------
      ird=0
      r=(ro+(3-0.5*ird)*rmax)/(1+3-0.5*ird)
      ir=0
      DO
         CALL direct_get_bfield(r,z,bf,1)
         dr=-bf%psi/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            ird=ird+1
            r=(ro+(3-0.5*ird)*rmax)/(1+3-0.5*ird)
            ir=0
            IF (ird==6) THEN 
               direct_infinite_loop_flag = .TRUE.
               CALL program_stop
     $              ("Took too many steps to find outb spx.")
            ENDIF
         ENDIF
      ENDDO
      rs2=r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_position
c-----------------------------------------------------------------------
c     subprogram 4. direct_fl_int.
c     integrates along field line.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_fl_int(psifac,eta1,eta2,y_out,bf,len_y_out)

      REAL(r8), INTENT(IN) :: psifac,eta1,eta2
      INTEGER, INTENT(OUT) :: len_y_out
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: y_out
      TYPE(direct_bfield_type), INTENT(OUT) :: bf

      CHARACTER(64) :: message,message2

      INTEGER, PARAMETER :: neq=4,liw=30,lrw=22+neq*16
      INTEGER :: iopt,istate,itask,itol,jac,mf
      INTEGER, DIMENSION(liw) :: iwork
      REAL(r8), PARAMETER :: eps=1e-12
      REAL(r8) :: atol,rtol,rfac,deta,r,z,eta,err,psi0
      REAL(r8), DIMENSION(neq) :: y
      REAL(r8), DIMENSION(lrw) :: rwork
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 20   FORMAT(/2x,"is",5x,"eta",8x,"deta",8x,"s",9x,"rfac",8x,"r",10x,
     $     "z",9x,"psi",8x,"err"/)
 30   FORMAT(i4,1p,8e11.3)
 40   FORMAT(a,i4,a,es10.3,a,i3)
 51   FORMAT(1x,"psifac =",es10.3)
 61   FORMAT(1x,"direct_int:",i6," steps taken of max",i6,".")
 11   FORMAT(1x,"Incomplete: eta=",es10.2," of [",es10.2,",",es10.2,"]")
c-----------------------------------------------------------------------
c     find flux surface.
c-----------------------------------------------------------------------
      CALL find_fl_surface(psifac,eta1,r,z)
      CALL direct_get_bfield(r,z,bf,1)
      psi0=bf%psi
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      istep=0
      eta=eta1
      deta=twopi/mtheta
      y=0
      y(2)=SQRT((r-ro)**2+(z-zo)**2)
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istate=1
      itask=5
      iopt=1
      mf=10
      itol=1
      rtol=etol
      atol=etol*y(2)
      iwork=0
      rwork=0
      rwork(1)=eta2
      rwork(11)=0
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      IF(out_fl)THEN
         WRITE(out_2d_unit,51)psifac
         WRITE(out_2d_unit,20)
      ENDIF
c-----------------------------------------------------------------------
c     store results for each step.
c-----------------------------------------------------------------------
      DO
         rfac=y(2)
         CALL direct_refine(rfac,eta,psi0)
         r=ro+rfac*COS(eta)
         z=zo+rfac*SIN(eta)
         CALL direct_get_bfield(r,z,bf,2)
         y_out(istep,:)=(/eta,y/)
         err=(bf%psi-psi0)/bf%psi
c-----------------------------------------------------------------------
c     compute and print output for each step.
c-----------------------------------------------------------------------
         IF(out_fl)WRITE(out_2d_unit,30)
     $        istep,eta,rwork(11),y(1:2),r,z,bf%psi,err
         IF(bin_fl)WRITE(bin_2d_unit)
     $        REAL(eta,4),REAL(rwork(11),4),REAL(y(1:4),4),REAL(r,4),
     $        REAL(z,4),REAL(psifac,4),REAL(err,4)
c-----------------------------------------------------------------------
c     advance differential equations.
c-----------------------------------------------------------------------
         IF(eta >= eta2 .OR. istep >= nstepd  .OR.  istate < 0
     $        .OR. ABS(err) >= 1)EXIT
         istep=istep+1
         CALL lsode(direct_fl_der,neq,y,eta,eta2,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
      IF(out_fl)WRITE(out_2d_unit,20)
      IF(bin_fl)WRITE(bin_2d_unit)
c-----------------------------------------------------------------------
c     abort if istep > nstepd.
c-----------------------------------------------------------------------
      IF(eta < eta2)THEN
         WRITE(message,61)istep,nstepd
         WRITE(message2,11)eta,eta1,eta2
         PRINT "(A)", message
         PRINT "(A)", "Increase nstepd or decrease etol."
         CALL program_stop(message2)
      ELSE
         WRITE(message,61)istep,nstepd
         IF(verbose)PRINT "(A)", message
      ENDIF
      len_y_out = istep
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_fl_int
c-----------------------------------------------------------------------
c     subprogram 4.5. find_fl_surface.
c     finds r,z given psi, eta.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE find_fl_surface(psifac,eta,r,z)

      REAL(r8), INTENT(IN) :: psifac,eta
      REAL(r8), INTENT(OUT) :: r,z

      TYPE(direct_bfield_type) :: bf
      REAL(r8), PARAMETER :: eps=1e-13
      INTEGER :: ir
      REAL(r8) :: cosfac,sinfac,radius,dradius,dfdradius,psi0
c-----------------------------------------------------------------------
c     find flux surface.
c-----------------------------------------------------------------------
      cosfac=COS(eta)
      sinfac=SIN(eta)

      psi0=psio*(1-psifac)
      radius=SQRT(psifac)*(rs2-ro)
      r=ro+cosfac*radius
      z=zo+sinfac*radius
      ir = 0
      DO
         CALL direct_get_bfield(r,z,bf,1)
         dfdradius = -bf%psir*cosfac-bf%psiz*sinfac
         dradius = -(psi0-bf%psi)/dfdradius

         radius=radius+dradius
         r=ro+cosfac*radius
         z=zo+sinfac*radius
         IF(ABS(dradius) <= eps*r)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find flux surf.")
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE find_fl_surface
c-----------------------------------------------------------------------
c     subprogram 5. direct_fl_der.
c     contains differential equations for field line averages.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_fl_der(neq,eta,y,dy)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: eta
      REAL(r8), INTENT(IN) :: y(neq)
      REAL(r8), INTENT(OUT) :: dy(neq)

      REAL(r8) :: cosfac,sinfac,bp,r,z,jacfac,bt,b
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     preliminary computations.
c     here magnetic coordinates are defined
c-----------------------------------------------------------------------
      cosfac=COS(eta)
      sinfac=SIN(eta)
      r=ro+y(2)*cosfac
      z=zo+y(2)*sinfac
      CALL direct_get_bfield(r,z,bf,1)
      bp=SQRT(bf%br**2+bf%bz**2)
      bt=bf%f/r
      b=SQRT(bp*bp+bt*bt)
      jacfac=bp**power_bp*b**power_b/r**power_r
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1)=y(2)/(bf%bz*cosfac-bf%br*sinfac)
      dy(2)=dy(1)*(bf%br*cosfac+bf%bz*sinfac)
      dy(3)=dy(1)/(r*r)
      dy(4)=dy(1)*jacfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_fl_der
c-----------------------------------------------------------------------
c     subprogram 6. direct_refine.
c     moves a point orthogonally to a specified flux surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_refine(rfac,eta,psi0)

      REAL(r8) :: rfac,eta,psi0

      REAL(r8) :: dpsi,cosfac,sinfac,drfac,r,z
      REAL(r8), PARAMETER :: eps=1e-12
      TYPE(direct_bfield_type) :: bf
      INTEGER :: ir
c-----------------------------------------------------------------------
c     initialize iteration.
c-----------------------------------------------------------------------
      cosfac=COS(eta)
      sinfac=SIN(eta)
      r=ro+rfac*cosfac
      z=zo+rfac*sinfac
      CALL direct_get_bfield(r,z,bf,1)
      dpsi=bf%psi-psi0
c-----------------------------------------------------------------------
c     refine rfac by newton iteration.
c-----------------------------------------------------------------------
      ir = 0
      DO
         drfac=-dpsi/(bf%psir*cosfac+bf%psiz*sinfac)
         rfac=rfac+drfac
         r=ro+rfac*cosfac
         z=zo+rfac*sinfac
         CALL direct_get_bfield(r,z,bf,1)
         dpsi=bf%psi-psi0
         IF(ABS(dpsi) <= eps*psi0 .OR. ABS(drfac) <= eps*rfac)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to refine rfac.")
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_refine
c-----------------------------------------------------------------------
c     subprogram 7. direct_output.
c     diagnoses input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_output

      INTEGER :: ix,iy
      REAL(r8), DIMENSION(:,:), POINTER :: x,y
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"ix",4x,"iy",6x,"r",10x,"z",9x,"psi"/)
 20   FORMAT(2i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(.NOT. (out_eq_1d .OR. bin_eq_1d .OR. out_eq_2d .OR. bin_eq_2d))
     $     RETURN
      IF(out_eq_1d .OR. out_eq_2d)
     $     CALL ascii_open(out_2d_unit,"input.out","UNKNOWN")
c-----------------------------------------------------------------------
c     diagnose 1d output.
c-----------------------------------------------------------------------
      IF(out_eq_1d)WRITE(out_2d_unit,'(a)')"input surface quantities:"
      IF(bin_eq_1d)CALL bin_open(bin_2d_unit,"sq_in.bin","UNKNOWN",
     $     "REWIND","none")
      CALL spline_write1(sq_in,out_eq_1d,bin_eq_1d,
     $     out_2d_unit,bin_2d_unit,interp)
      IF(bin_eq_1d)CALL bin_close(bin_2d_unit)
c-----------------------------------------------------------------------
c     ascii table of psi.
c-----------------------------------------------------------------------
      IF(out_eq_2d)THEN
         DO iy=0,psi_in%my
            WRITE(out_2d_unit,10)
            DO ix=0,psi_in%mx
               WRITE(out_2d_unit,20)
     $              ix,iy,psi_in%xs(ix),psi_in%ys(iy),psi_in%fs(ix,iy,1)
            ENDDO
         ENDDO
         WRITE(out_2d_unit,10)
      ENDIF
c-----------------------------------------------------------------------
c     draw contour plot of psi.
c-----------------------------------------------------------------------
      IF(bin_eq_2d)THEN
         CALL bin_open(bin_2d_unit,"psi_in.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)psi_in%mx,psi_in%my
         ALLOCATE(x(0:psi_in%mx,0:psi_in%my),y(0:psi_in%mx,0:psi_in%my))
         DO ix=0,psi_in%mx
            DO iy=0,psi_in%my
               x(ix,iy)=psi_in%xs(ix)
               y(ix,iy)=psi_in%ys(iy)
            ENDDO
         ENDDO
         WRITE(bin_2d_unit)REAL(x,4),REAL(y,4)
         WRITE(bin_2d_unit)REAL(psi_in%fs,4)
         DEALLOCATE(x,y)
         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(out_eq_1d .OR. out_eq_2d)CALL ascii_close(out_2d_unit)
      IF(input_only)CALL program_stop("Termination by direct_output.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_output
c-----------------------------------------------------------------------
c     subprogram 8. direct_local_xpoint.
c     finds location of nearby x-point where |Bp|=0 using Newton method.
c     can and will search outside separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_local_xpoint(r,z)

      REAL(r8), INTENT(INOUT) :: r,z
      REAL(r8) :: ajac(2,2),det,dr,dz
      REAL(r8), PARAMETER :: eps=1e-13
      TYPE(direct_bfield_type) :: bf
      INTEGER :: ir
c-----------------------------------------------------------------------
c     use newton iteration to find x-point.
c-----------------------------------------------------------------------
      ir = 0
      DO
         CALL direct_get_bfield(r,z,bf,2)
         ajac(1,1)=bf%brr
         ajac(1,2)=bf%brz
         ajac(2,1)=bf%bzr
         ajac(2,2)=bf%bzz
         det=ajac(1,1)*ajac(2,2)-ajac(1,2)*ajac(2,1)
         dr=(ajac(1,2)*bf%bz-ajac(2,2)*bf%br)/det
         dz=(ajac(2,1)*bf%br-ajac(1,1)*bf%bz)/det
         r=r+dr
         z=z+dz
         IF(ABS(dr) <= eps*r .AND. ABS(dz) <= eps*r)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find x-point.")
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_local_xpoint
c-----------------------------------------------------------------------
c     subprogram 9. direct_saddle_angle_DEPRECATED.
c     finds angle location of nearby saddle-node eigenvector using 
c     Newton iteration. Ends up infinite looping... I suspect the 
c     analytic form of the zero crosssing of Bout is unfriendly to 
c     the Newton method (I double checked formulas etc)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_saddle_angle_DEPRECATED(rx,zx,nu,Bnorm)

      REAL(r8), INTENT(IN) :: rx,zx,Bnorm
      REAL(r8), INTENT(INOUT) :: nu
      REAL(r8), PARAMETER :: r_eps1=1e-9,r_eps2=1e-10,nu_eps=1e-13
      INTEGER :: ir
      REAL(r8) :: Bout,dBout_dnu,cosfac,sinfac,dnu
      TYPE(direct_bfield_type) :: bf
      
c-----------------------------------------------------------------------
c     use newton iteration to find point where Bout = 0.
c-----------------------------------------------------------------------
c     vtheta = -sin(theta)*Rhat+cos(theta)*Zhat
c     B = bf%br*Rhat+bf%bz*Zhat
c     CHJECK LINEARITY OF LEG?
c     ONLY TWO X-pts, can check if need the good treatment for them both
      ir=0
      DO 
         cosfac=COS(nu)
         sinfac=SIN(nu)
         CALL direct_get_bfield(rx+r_eps1*rx*cosfac,
     $                          zx+r_eps1*rx*sinfac,bf,1)
         Bout = (-sinfac*bf%br+cosfac*bf%bz)
         dBout_dnu = ((-cosfac*bf%br-sinfac*bf%bz) +
     $      (-sinfac*bf%brr+cosfac*bf%bzr)*(-r_eps1*rx*sinfac) +
     $      (-sinfac*bf%brz+cosfac*bf%bzz)*(r_eps1*rx*cosfac))
         dnu = -Bout/dBout_dnu
         nu=nu+dnu
         IF(ABS(dnu) <= nu_eps)EXIT

         PRINT "(e16.3)", Bout/Bnorm

         ir = ir+1
         IF (ir  > 50) THEN !direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find x-pt angle.")
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     make sure nu is in [0,2pi). works even if nu is negative
c-----------------------------------------------------------------------
      nu = nu - twopi*floor(nu/twopi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_saddle_angle_DEPRECATED
c-----------------------------------------------------------------------
c     subprogram 9. direct_saddle_angle.
c     finds angle location of nearby saddle-node eigenvector using 
c     a binary search type algorithm. can find point where Bnu = 0 or 
c     Brho = 0.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_saddle_angle(rx,zx,rho,nustart_in,nu_var_in,nu,
     $                                                     Bcase,debug)

      REAL(r8), INTENT(IN) :: rx,zx,rho,nustart_in,nu_var_in
      CHARACTER, INTENT(IN) :: Bcase
      LOGICAL, INTENT(IN) :: debug
      REAL(r8), INTENT(OUT) :: nu

      INTEGER, PARAMETER :: ird=4
      REAL(r8), PARAMETER :: nu_eps=1e-13

      REAL(r8), DIMENSION(ird) :: nus,nu_tmp
      REAL(r8) :: bnorm,Bout,nustep,pos,nustart,nu_var
      INTEGER :: i,inuh


      nustart=nustart_in
      nu_var=nu_var_in
      pos=1.0
      inuh=0
c-----------------------------------------------------------------------
c     looping over narrower nu-intervals
c-----------------------------------------------------------------------
      CALL direct_Blocal(rx,zx,nustart_in,rho,Bcase,bnorm)
      IF(bnorm<zero)pos=-1.0
c-----------------------------------------------------------------------
c     looping to narrow nu-intervals
c-----------------------------------------------------------------------
      DO 
c-----------------------------------------------------------------------
c     generating vector of nus to search along for Bout sign change
c-----------------------------------------------------------------------
         nustep = nu_var/(ird-1)
         nu_tmp = (/(i-1,i=1,ird)/)
         nu_tmp = nu_tmp*nustep
         nus = nu_tmp+nustart
         IF(debug)PRINT "(A)", "inuh"
         IF(debug)PRINT "(i6)", inuh
c-----------------------------------------------------------------------
c     loop along nus from nustart to nustart+nu_var, seaching for Bout 
c     sign change
c-----------------------------------------------------------------------
         DO i=1,ird,+1
            CALL direct_Blocal(rx,zx,nus(i),rho,Bcase,Bout)

c-----------------------------------------------------------------------
c     debug print statement
c-----------------------------------------------------------------------
            IF(debug)THEN
               PRINT "(A)", "i"
               PRINT "(i6)", i
               PRINT "(A)", "nu"
               PRINT "(es16.10)", nus(i)
               PRINT "(A)", "Bout"
               PRINT "(f16.9)", Bout/bnorm
            ENDIF
c-----------------------------------------------------------------------
c     tightening search bracket, cyling original do loop
c-----------------------------------------------------------------------
            IF(Bout*pos<0.0)THEN
               nustart=nus(i-1)
               nu_var=nus(i)-nus(i-1)
               nu=(nus(i)+nus(i-1))/2
               CALL direct_Blocal(rx,zx,nu,rho,Bcase,Bout)

               IF(debug)THEN
                  PRINT "(A)", "Delta nu"
                  PRINT "(es16.10)", ABS(nu-nustart)
               ENDIF
               EXIT
            ENDIF
            IF(i==ird) CALL program_stop("couldn't find x-pt angle")
         ENDDO

         IF(ABS(nu-nustart)<nu_eps)EXIT

         inuh = inuh+1
         IF (inuh  > 80)THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find x-pt angle.")
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_saddle_angle
c-----------------------------------------------------------------------
c     subprogram 10. direct_psisaddle.
c     calculates the linear term of psi_in at the saddle point, as well
c     as gamma, and vartheta
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_psisaddle(rx,zx,oangle,
     $                              nu,b11,gamma,vartheta,lincheck)

      REAL(r8), INTENT(IN) :: rx,zx,oangle
      REAL(r8), DIMENSION(2), INTENT(IN) :: nu
      REAL(r8), INTENT(OUT) :: b11,lincheck,gamma,vartheta

      REAL(r8), PARAMETER :: nuh_eps=1e-13, nuh_eps2=1e-6
      INTEGER :: ir
      REAL(r8) :: nuh,nuperf,cosfac,sinfac,cosfact,sinfact
      REAL(r8) :: psix,psinuh,r,Rlocal,Zlocal,x,y,chi
      REAL(r8), DIMENSION(4) :: r_eps
      TYPE(direct_bfield_type) :: bf
      LOGICAL :: debug=.FALSE.

      r_eps(1) = 1e-7
      r_eps(2) = 1e-8
      r_eps(3) = 1e-9
      r_eps(4) = 1e-10
      ir=3
c-----------------------------------------------------------------------
c     finding correct initialisation point for nuh (read nu-half).
c     gamma should be less than pi for a real x-point
c-----------------------------------------------------------------------
      nuperf = (nu(1)+nu(2))/2
      gamma = nu(1)-nu(2)
      nuh = nuperf

      IF (gamma  > pi) THEN
            CALL program_stop("Angle between separatrix legs > pi.")
      ENDIF
c-----------------------------------------------------------------------
c     finding angle where Brho = 0.
c-----------------------------------------------------------------------
      CALL direct_saddle_angle(rx,zx,rx*r_eps(ir),nu(2),nu(1)-nu(2)
     $                                                 ,nuh,'r',.FALSE.)
c-----------------------------------------------------------------------
c     make sure nuh is in [0,2pi). works even if nuh is negative.
c-----------------------------------------------------------------------
      nuh = nuh - twopi*floor(nuh/twopi)
c-----------------------------------------------------------------------
c     in the linear approaximation of psi at the x-point, Brho=0 halfway
c     between the two separatrix legs. we check if we are within 
c     nuh_eps2 of this case
c-----------------------------------------------------------------------
      lincheck = abs(nuh-nuperf)
      IF(debug .OR. lincheck > nuh_eps2)THEN
         PRINT"(A)","nu value where Brho=0 deviates from linear case by"
         PRINT "(es16.10)", lincheck
      ENDIF
      IF(lincheck > nuh_eps2)THEN
         CALL program_stop("psi isn't well approximated at the x-point")
      ENDIF
c-----------------------------------------------------------------------
c     defining rotated saddle-point coordinate frame to extract linear
c     component.
c-----------------------------------------------------------------------
      IF (oangle>nu(2)) THEN
         vartheta = nu(1)-pi/2.0  
      ELSE 
         vartheta = nu(2)-pi/2.0  
      ENDIF
      vartheta = vartheta - twopi*floor(vartheta/twopi)

      r = r_eps(ir)*rx
      cosfac=COS(nuh)
      sinfac=SIN(nuh)
      cosfact=COS(vartheta)
      sinfact=SIN(vartheta)

      Rlocal = r*cosfac
      Zlocal = r*sinfac
      x = cosfact*Rlocal + sinfact*Zlocal
      y = -sinfact*Rlocal + cosfact*Zlocal
      chi = -COS(gamma)*x+SIN(gamma)*y
c-----------------------------------------------------------------------
c     extracting linear component b11, where psi = psi(rx,zx)+b11*x*chi.
c-----------------------------------------------------------------------
      CALL direct_get_bfield(rx,zx,bf,1)
      psix = bf%psi
      CALL direct_get_bfield(rx+r*cosfac,zx+r*sinfac,bf,1)
      psinuh = bf%psi

      b11 = (psinuh-psix)/(x*chi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_psisaddle
c-----------------------------------------------------------------------
c     subprogram 11. direct_xpoint.
c     finds location and angles of nearby x-point, checks if it's inside
C     separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_xpoint(rin,zin,x_i)

      REAL(r8), INTENT(IN) :: rin,zin
      INTEGER, INTENT(IN) :: x_i

      INTEGER, PARAMETER :: ird=4
      REAL(r8), PARAMETER :: psi_eps=1e-4, r_eps1=1e-9
      REAL(r8) :: r,z
      REAL(r8) :: b11,lincheck,gamma,vartheta
      REAL(r8), DIMENSION(1:2) :: nu
      REAL(r8) :: oangle,nu_var,Bnua,Bnub,Bnuc
      TYPE(direct_bfield_type) :: bf
      LOGICAL :: test_direct_local_xpoint,test_direct_saddle_angle

      r=rin
      z=zin
      nu=0.0
      test_direct_local_xpoint=.FALSE.
      test_direct_saddle_angle=.FALSE.
c-----------------------------------------------------------------------
c     testing direct_local_xpoint.
c-----------------------------------------------------------------------
      IF(test_direct_local_xpoint)THEN
         CALL direct_get_bfield(ro,zo,bf,1)
         PRINT "(A)", "psi at origin"
         PRINT "(f16.10)", bf%psi
         PRINT "(A)", "r,z b4 direct_local_xpoint"
         PRINT "(f16.5)", r
         PRINT "(f16.5)", z
         CALL direct_get_bfield(r,z,bf,1)
         PRINT "(A)", "|Bp| b4 direct_local_xpoint"
         PRINT "(f16.10)", SQRT(bf%br**2+bf%bz**2)
         PRINT "(A)", "psi b4 direct_local_xpoint"
         PRINT "(f16.10)", bf%psi
         CALL direct_local_xpoint(r,z)
         PRINT "(A)", "r,z aftr direct_local_xpoint"
         PRINT "(f16.5)", r
         PRINT "(f16.5)", z
         CALL direct_get_bfield(r,z,bf,1)
         PRINT "(A)", "|Bp| after direct_local_xpoint"
         PRINT "(f16.10)", SQRT(bf%br**2+bf%bz**2)
         PRINT "(A)", "psi after direct_local_xpoint"
         PRINT "(f16.10)", bf%psi
      ENDIF
c-----------------------------------------------------------------------
c     finds x-point and fills out global module variables
c-----------------------------------------------------------------------
      CALL direct_local_xpoint(r,z)
      CALL direct_get_bfield(r,z,bf,1)
      rxs(x_i) = r
      zxs(x_i) = z
c-----------------------------------------------------------------------
c     checks if outside separatrix. ended up not needing...
c-----------------------------------------------------------------------
      !IF (bf%psi < -psi_eps*psio) THEN
      !   near_sep = .FALSE. 
      !   RETURN
      !ENDIF
c-----------------------------------------------------------------------
c     updating xpt_etas with a more exact value than that calculated by
c     direct_initialise_xpoints.
c-----------------------------------------------------------------------
      xpt_etas(x_i)=ATAN2(zxs(x_i)-zo,rxs(x_i)-ro)
      xpt_etas(x_i)=xpt_etas(x_i) - twopi*floor(xpt_etas(x_i)/twopi)
c-----------------------------------------------------------------------
c     finding angle of o-point from perspective of x-point (oangle).
c-----------------------------------------------------------------------
      oangle = xpt_etas(x_i)+pi
      oangle = oangle - twopi*floor(oangle/twopi)
      nu_var = pi/2
c-----------------------------------------------------------------------
c     calculating separatrix leg angles (nu). nu(1) is always more than 
c     oangle, nu(2) is always less than oangle.
c-----------------------------------------------------------------------
      CALL direct_saddle_angle(rxs(x_i),zxs(x_i),r_eps1*rxs(x_i),oangle,
     $                                     nu_var,nu(1),'n',.FALSE.)
      CALL direct_saddle_angle(rxs(x_i),zxs(x_i),r_eps1*rxs(x_i),oangle,
     $                                     -nu_var,nu(2),'n',.FALSE.)
c-----------------------------------------------------------------------
c     testing direct_saddle_angle AGAGAG
c-----------------------------------------------------------------------
      IF(test_direct_saddle_angle)THEN
         CALL direct_Blocal(rxs(x_i),zxs(x_i),nu(1),r_eps1*rxs(x_i),
     $                                                     'n',Bnua)
         CALL direct_Blocal(rxs(x_i),zxs(x_i),oangle,r_eps1*rxs(x_i),
     $                                                     'n',Bnub)
         CALL direct_Blocal(rxs(x_i),zxs(x_i),nu(2),r_eps1*rxs(x_i),
     $                                                     'n',Bnuc)
         PRINT "(A)", "First x-point leg's angle nu:"
         PRINT "(f16.3)", nu(1)/pi
         PRINT "(A)", "First leg B_nu"
         PRINT "(e16.3)", Bnua
         PRINT "(A)", "nu angle pointing from x-pt to mag. axis:"
         PRINT "(f16.3)", oangle/pi
         PRINT "(A)", "B_nu at angle pointing from x-pt to mag. axis:"
         PRINT "(e16.3)", Bnub
         PRINT "(A)", "Second x-point leg's angle nu:"
         PRINT "(f16.3)", nu(2)/pi
         PRINT "(A)", "Second leg B_nu"
         PRINT "(e16.3)", Bnuc
         PRINT "(A)", "Angle between x-point legs (gamma)"
         PRINT "(f16.3)", (nu(1)-nu(2))/pi
      ENDIF
c-----------------------------------------------------------------------
c     calculating xpoint angles gamma, vartheta, and linear term b11.
c-----------------------------------------------------------------------
      CALL direct_psisaddle(r,z,oangle,nu,b11,gamma,vartheta,lincheck)
c-----------------------------------------------------------------------
c     filling out x-point module variables.
c-----------------------------------------------------------------------
      xpt_b11s(x_i) = b11
      xpt_gammas(x_i) = gamma
      xpt_varthetas(x_i) = vartheta
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_xpoint
c-----------------------------------------------------------------------
c     subprogram 12. direct_initialise_xpoints.
c     scans the field-line integral to find likely x-point locations,
c     denoted by diverging regions in the q-integral
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_initialise_xpoints(y_out,len_y_out,wrap,debug,
     $                    bf,dq_eps,eta_maxes,eta_brackets,maxima_count)

      INTEGER, INTENT(IN) :: len_y_out
      LOGICAL, INTENT(IN) :: wrap,debug
      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: y_out
      REAL(r8), INTENT(IN) :: dq_eps
      TYPE(direct_bfield_type), INTENT(IN) :: bf

      INTEGER, INTENT(OUT) :: maxima_count
      REAL(r8), DIMENSION(1:), INTENT(OUT) :: eta_maxes
      REAL(r8), DIMENSION(1:,1:), INTENT(OUT) :: eta_brackets

      INTEGER :: i
      LOGICAL :: prev_above_threshold,wrap_max,wrap_
      REAL(r8), DIMENSION(1:len_y_out) :: dqdeta
      REAL(r8), DIMENSION(100) :: max_dqdeta,max_locations

      max_locations=0.0
      eta_maxes=0.0
      maxima_count=0
      wrap_max=.FALSE.
      wrap_=wrap
      prev_above_threshold=.FALSE.
c-----------------------------------------------------------------------
c     sanity check that final eta is less than twopi. 
c     only called when wrap = .TRUE. and we are analysing the whole
c     [0,2pi] interval. lsode integrator makes sure final eta is twopi
c-----------------------------------------------------------------------
      IF ((y_out(len_y_out,0)-y_out(0,0)) > twopi .AND. wrap) THEN
         CALL program_stop("eta wrap error... debug direct.f")
      ENDIF
c-----------------------------------------------------------------------
c     turn y_out into dqdeta.
c-----------------------------------------------------------------------
      DO i=0,(len_y_out-1),+1
         dqdeta(i+1) = (bf%f)*(y_out(i+1,3)-y_out(i,3))
     $                 /(y_out(i+1,0)-y_out(i,0))!
      ENDDO
c-----------------------------------------------------------------------
c     identifies all regions of eta above threshold. calculates the 
c     maximum value of d_q_int/d_eta as well as eta for each region.
c-----------------------------------------------------------------------
      DO i=1,len_y_out-1,+1
         IF((dqdeta(i)>dq_eps).AND.(.NOT.prev_above_threshold)) THEN
            IF(i==1)THEN
               wrap_max=.TRUE.
            ENDIF
            maxima_count=maxima_count+1

            eta_brackets(maxima_count,1)=y_out(i-1,0)
            eta_brackets(maxima_count,2)=y_out(i,0)

            max_dqdeta(maxima_count)=dqdeta(i)
            max_locations(maxima_count) = 0.5*(y_out(i,0)+
     $                                         y_out(i-1,0))
            prev_above_threshold=.TRUE.
         ELSEIF ((dqdeta(i)>dq_eps) .AND. prev_above_threshold) THEN
            eta_brackets(maxima_count,2)=y_out(i,0)
            
            IF (dqdeta(i)>max_dqdeta(maxima_count)) THEN
               max_dqdeta(maxima_count)=dqdeta(i)
               max_locations(maxima_count) = 0.5*(y_out(i,0)+
     $                                            y_out(i-1,0))
            ENDIF
            prev_above_threshold=.TRUE.
         ELSE
            prev_above_threshold=.FALSE.
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     case where whole range is under threshold doesn't require wrapping
c-----------------------------------------------------------------------
      IF(MINVAL(dqdeta) > dq_eps)THEN
            wrap_=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     final datapoint is a special case due to the possibility of the 
c     region wrapping around to eta=0 region.
c-----------------------------------------------------------------------
      IF (dqdeta(len_y_out)>dq_eps .AND. wrap_) THEN
            IF (prev_above_threshold) THEN 
               IF (dqdeta(len_y_out)>max_dqdeta(maxima_count)) THEN
                  max_dqdeta(maxima_count)=dqdeta(len_y_out)
                  max_locations(maxima_count) = 
     $                  0.5*(y_out(len_y_out-1,0)+y_out(len_y_out,0))
               ENDIF
               eta_brackets(maxima_count,2)=y_out(0,0)

               IF (wrap_max) THEN
                  eta_brackets(1,1)=eta_brackets(maxima_count,1)

                  IF (max_dqdeta(maxima_count)>=max_dqdeta(1)) THEN
                     max_dqdeta(1)=max_dqdeta(maxima_count)
                     max_locations(1)=max_locations(maxima_count)
                  ENDIF

                  eta_brackets(maxima_count,1)=0.0
                  eta_brackets(maxima_count,2)=0.0

                  max_dqdeta(maxima_count)=0.0
                  max_locations(maxima_count)=0.0
                  maxima_count=maxima_count-1
               ENDIF

            ELSEIF (wrap_max) THEN !.AND. !prev_above_threshold)
               eta_brackets(1,1)=y_out(len_y_out-1,0)

               IF (dqdeta(len_y_out)>max_dqdeta(1)) THEN
                  max_dqdeta(1)=dqdeta(len_y_out)
                  max_locations(1)=0.5*(y_out(len_y_out-1,0)
     $                                  +y_out(len_y_out,0))        
               ENDIF
            ELSE
               maxima_count=maxima_count+1

               eta_brackets(maxima_count,1)=y_out(len_y_out-1,0)
               eta_brackets(maxima_count,2)=y_out(0,0)

               max_dqdeta(maxima_count)=dqdeta(len_y_out)
               max_locations(maxima_count) = 
     $                  0.5*(y_out(len_y_out-1,0)+y_out(len_y_out,0))
            ENDIF
      ELSEIF (dqdeta(len_y_out)>dq_eps) THEN
c-----------------------------------------------------------------------
c     non wrap case treated same as i=1:len_y_out-1
c-----------------------------------------------------------------------
         IF(prev_above_threshold) THEN
            eta_brackets(maxima_count,2)=y_out(len_y_out,0)
            
            IF (dqdeta(len_y_out)>max_dqdeta(maxima_count)) THEN
               max_dqdeta(maxima_count)=dqdeta(len_y_out)
               max_locations(maxima_count) = 0.5*(y_out(len_y_out,0)+
     $                                            y_out(len_y_out-1,0))
            ENDIF
            prev_above_threshold=.TRUE.
         ELSE
            maxima_count=maxima_count+1

            eta_brackets(maxima_count,1)=y_out(len_y_out-1,0)
            eta_brackets(maxima_count,2)=y_out(len_y_out,0)

            max_dqdeta(maxima_count)=dqdeta(len_y_out)
            max_locations(maxima_count) = 0.5*(y_out(len_y_out,0)+
     $                                         y_out(len_y_out-1,0))
            prev_above_threshold=.TRUE.
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     warning if theres more than 2 xpoints.
c-----------------------------------------------------------------------
      IF(maxima_count>2 .AND. wrap)THEN
         CALL program_stop("More than 2 x-points detected.")
      ENDIF
c-----------------------------------------------------------------------
c     filling out xpt_etas output.
c-----------------------------------------------------------------------
      IF (maxima_count==0) THEN
         xpt_etas(1)=0.0
      ELSE
         DO i=1,maxima_count,+1
            eta_maxes(i)=max_locations(i)
            IF(wrap)THEN
               xpt_etas(i)=max_locations(i)
               xpt_etas(i)=xpt_etas(i)-twopi*floor(xpt_etas(i)/twopi)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     debug print statements.
c-----------------------------------------------------------------------

      IF(debug)PRINT "(A)", "Num points::::::::::::::::::::::::::::::::"
      IF(debug)PRINT "(i6)", maxima_count
      IF(maxima_count>0 .AND. debug)THEN
         DO i=1,maxima_count
            PRINT "(A)", "Point "
            PRINT "(i6)", i
            PRINT "(f16.5)", (eta_brackets(i,1)
     $            )!-twopi*floor(eta_brackets(i,1)/twopi))
            PRINT "(f16.5)", (eta_maxes(i)
     $            )!-twopi*floor(eta_maxes(i)/twopi))
            PRINT "(f16.5)", (eta_brackets(i,2)
     $            )!-twopi*floor(eta_brackets(i,2)/twopi))
         ENDDO
      ENDIF
      IF (debug)THEN
         PRINT "(A)", "Total theta values over interval "
         PRINT "(f16.5)", y_out(0,0)
         PRINT "(f16.5)", y_out(len_y_out,0)
         PRINT "(A)", "Total q values over interval "
         PRINT "(es16.5)", y_out(0,3)
         PRINT "(es16.5)", y_out(len_y_out,3)
         PRINT "(es16.5)", (y_out(len_y_out,3)/(y_out(len_y_out,0)-
     $                                          y_out(0,0)))
         PRINT "(A)", "max and min dqdeta"
         PRINT "(es16.5)", MAXVAL(dqdeta)
         PRINT "(es16.5)", MINVAL(dqdeta)
      ENDIF
      IF(debug)PRINT "(A)", "::::::::::::::::::::::::::::::::::::::::::"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_initialise_xpoints
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 17. direct_Blocal.
c     calculates local B-field displaced from some r,z point in polar
c     coordinates. either Bnu or Brho depending on Bcase
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_Blocal(r,z,nu,rho,Bcase,Bout)

      REAL(r8), INTENT(IN) :: r,z,nu,rho
      CHARACTER, INTENT(IN) :: Bcase
      REAL(r8), INTENT(OUT) :: Bout
      REAL(r8) :: cosfac,sinfac
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     calculation of Bout.
c-----------------------------------------------------------------------
      cosfac=COS(nu)
      sinfac=SIN(nu)
      CALL direct_get_bfield(r+rho*cosfac,z+rho*sinfac,bf,1)
      SELECT CASE (Bcase)
         CASE ('n')
         Bout=-sinfac*bf%br+cosfac*bf%bz
         CASE ('r')
         Bout=cosfac*bf%br+sinfac*bf%bz
         CASE DEFAULT
         CALL program_stop("Invalid Bcase in direct_Blocal")  
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_Blocal
      END MODULE direct_mod
