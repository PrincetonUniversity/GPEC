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
c     5. direct_fl_der.
c     6. direct_refine.
c     7. direct_output.
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
      TYPE(bicube_type) :: psi_in
      LOGICAL :: direct_infinite_loop_flag
      INTEGER :: direct_infinite_loop_count = 2000

      TYPE :: direct_bfield_type
      REAL(r8) :: psi,psir,psiz,psirz,psirr,psizz,f,f1,p,p1
      REAL(r8) :: br,bz,brr,brz,bzr,bzz
      END TYPE direct_bfield_type

      REAL(r8) :: etol=1e-8

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. direct_run.
c     gets equilibrium data and massages it.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE direct_run

      INTEGER :: ir,iz,itheta,ipsi
      INTEGER, PARAMETER :: nstep=2048
      REAL(r8) :: f0fac,f0,ffac,rfac,eta,r,jacfac,w11,w12,delpsi,q
      REAL(r8), DIMENSION(0:nstep,0:4) :: y_out
      REAL(r8), DIMENSION(2, mpsi+1) :: xdx
      REAL(r8), DIMENSION(3,3) :: v

      REAL(r8) :: xm,dx,rholow,rhohigh
      TYPE(direct_bfield_type) :: bf
      TYPE(spline_type) :: ff
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
c     start loop over flux surfaces and integrate over field line.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,'(a,1p,es10.3)')" etol = ",etol
      DO ipsi=mpsi,0,-1
         CALL direct_fl_int(ipsi,y_out,bf)
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
         IF(ipsi == mpsi)THEN
            IF(mtheta == 0)mtheta=istep
            CALL bicube_alloc(rzphi,mpsi,mtheta,4)
            CALL bicube_alloc(eqfun,mpsi,mtheta,3) ! new eq information
            rzphi%xs=sq%xs
            rzphi%ys=(/(itheta,itheta=0,mtheta)/)/REAL(mtheta,r8)
            rzphi%xtitle="psifac"
            rzphi%ytitle="theta "
            rzphi%title=(/"  r2  "," deta "," dphi ","  jac "/)
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
      SUBROUTINE direct_fl_int(ipsi,y_out,bf)

      INTEGER, INTENT(IN) :: ipsi
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: y_out
      TYPE(direct_bfield_type), INTENT(OUT) :: bf

      CHARACTER(64) :: message

      INTEGER, PARAMETER :: neq=4,liw=30,lrw=22+neq*16
      INTEGER :: iopt,istate,itask,itol,jac,mf,ir
      INTEGER, DIMENSION(liw) :: iwork
      INTEGER, PARAMETER :: nstep=2048
      REAL(r8), PARAMETER :: eps=1e-12
      REAL(r8) :: atol,rtol,rfac,deta,r,z,eta,err,psi0,psifac,dr
      REAL(r8), DIMENSION(neq) :: y
      REAL(r8), DIMENSION(lrw) :: rwork
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"ipsi = ",i3,", psifac =",es10.3)
 20   FORMAT(/2x,"is",5x,"eta",8x,"deta",8x,"s",9x,"rfac",8x,"r",10x,
     $     "z",9x,"psi",8x,"err"/)
 30   FORMAT(i4,1p,8e11.3)
 40   FORMAT(a,i4,a,es10.3,a,i3)
c-----------------------------------------------------------------------
c     find flux surface.
c-----------------------------------------------------------------------
      psifac=sq%xs(ipsi)
      psi0=psio*(1-psifac)
      r=ro+SQRT(psifac)*(rs2-ro)
      z=zo
      ir = 0
      DO
         CALL direct_get_bfield(r,z,bf,1)
         dr=(psi0-bf%psi)/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r)EXIT

         ir = ir+1
         IF (ir  > direct_infinite_loop_count) THEN
            direct_infinite_loop_flag = .TRUE.
            CALL program_stop("Took too many steps to find flux surf.")
         ENDIF
      ENDDO
      psi0=bf%psi
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      istep=0
      eta=0
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
      rwork(1)=twopi
      rwork(11)=0
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      IF(out_fl)THEN
         WRITE(out_2d_unit,10)ipsi,psifac
         WRITE(out_2d_unit,20)
      ENDIF
c-----------------------------------------------------------------------
c     store results for each step.
c-----------------------------------------------------------------------
      DO
         rfac=y(2)
         CALL direct_refine_itp(rfac,eta,psi0)
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
     $        REAL(eta,4),REAL(y(1:2),4),REAL(r,4),REAL(z,4),REAL(err,4)
c-----------------------------------------------------------------------
c     advance differential equations.
c-----------------------------------------------------------------------
         IF(eta >= twopi .OR. istep >= nstep  .OR.  istate < 0
     $        .OR. ABS(err) >= 1)EXIT
         istep=istep+1
         CALL lsode(direct_fl_der,neq,y,eta,twopi,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
      IF(out_fl)WRITE(out_2d_unit,20)
      IF(bin_fl)WRITE(bin_2d_unit)
c-----------------------------------------------------------------------
c     abort if istep > nstep.
c-----------------------------------------------------------------------
      IF(eta < twopi)THEN
         WRITE(message,40)"direct_int: istep = nstep = ",nstep,
     $        " at eta = ",eta,", ipsi = ",ipsi
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE direct_fl_int
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
c     subprogram 6. psi_at_polar_loc
c     returns psi at a position in r,eta coordinate space
c-----------------------------------------------------------------------
      FUNCTION psi_at_polar_loc(rfac,eta)
        REAL(r8) :: rfac,eta
        REAL(r8) :: cosfac,sinfac,r,z
        REAL(r8) :: psi_at_polar_loc
        TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
        cosfac=COS(eta)
        sinfac=SIN(eta)
        r=ro+rfac*cosfac
        z=zo+rfac*sinfac
        CALL direct_get_bfield(r,z,bf,1)

        psi_at_polar_loc = bf%psi
      END FUNCTION psi_at_polar_loc

c-----------------------------------------------------------------------
c     subprogram 6. direct_refine_itp
c     moves a point orthogonally to a specified flux surface using ITP
c-----------------------------------------------------------------------
      SUBROUTINE direct_refine_itp(rfac,eta,psi0)
        REAL(8), INTENT(INOUT) :: rfac
        REAL(8), INTENT(IN) :: eta, psi0
        REAL(8), PARAMETER :: margin = 0.99

        REAL(8) :: a, b

        a = rfac*(1-margin)
        b = rfac*(1+margin)
        CALL ITP( a, b, eta, psi0, rfac )

      END SUBROUTINE direct_refine_itp

c-----------------------------------------------------------------------
c     subprogram 4. ITP
c     ITP (Interpolate-Truncate-Project) root-finding algorithm
c-----------------------------------------------------------------------
      SUBROUTINE ITP(a, b, eta, psi0, root)
        REAL(8), INTENT(IN) :: a, b, eta, psi0
        REAL(8), INTENT(OUT) :: root
        REAL(8), PARAMETER :: tol = 1e-12
    
        REAL(8) :: xL, xU, fL, fU, xM, fM, xITP
        INTEGER :: iter
        REAL(8) :: half, sigma, adjust, lambda
    
        xL = a
        xU = b
        fL = psi_at_polar_loc(xL,eta) - psi0
        fU = psi_at_polar_loc(xU,eta) - psi0
c        success = .FALSE.
    
        IF (fL * fU > 0.0) THEN
          PRINT *, 'ITP Error: f(a) and f(b) must have opposite signs'
          PRINT *, 'a = ', a, ' f(a) = ', fL
          PRINT *, 'b = ', b, ' f(b) = ', fU
          RETURN
        END IF
    
        DO iter = 1, direct_infinite_loop_count
          xM = (xL + xU) / 2.0  ! Midpoint
          fM = psi_at_polar_loc(xM,eta) - psi0
      
          half = (xU - xL) / 2.0
          sigma = SIGN(1.0d0, fL - fU)
          adjust = (xU - xL) / 2.0
          lambda = MAX(0.0, MIN(adjust, half - tol))
          xITP = xM - sigma * lambda
      
          IF (ABS(fM) < tol .OR. (xU - xL) < tol) THEN
            root = xM
c            success = .TRUE.
            RETURN
          END IF
      
          IF (fM * fL > 0.0) THEN
            xL = xM
            fL = fM
          ELSE
            xU = xM
            fU = fM
          END IF
        END DO
        CALL program_stop("Took too many steps to refine rfac via ITP.")
    
c        root = xM
      END SUBROUTINE ITP
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
      END MODULE direct_mod
