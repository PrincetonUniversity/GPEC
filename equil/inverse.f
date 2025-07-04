c-----------------------------------------------------------------------
c     file inverse.f
c     converts inverse equilibrium to straight-fieldline coordinates.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inverse_mod.
c     1. inverse_run.
c     2. inverse_output.
c     3. inverse_extrap.
c     4. inverse_chease4_run.
c-----------------------------------------------------------------------
c     subprogram 0. inverse_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inverse_mod
      USE global_mod
      USE utils_mod
      USE grid_mod
      IMPLICIT NONE

      TYPE(bicube_type) :: rz_in

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. inverse_run.
c     gets equilibrium data and massages it.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inverse_run

      LOGICAL, PARAMETER :: interp=.FALSE.,
     $     diagnose_rz_in=.FALSE.,diagnose_rzphi=.FALSE.
      INTEGER :: ipsi,itheta,mx,my
      INTEGER, PARAMETER :: me=3
      REAL(r8) :: psifac,theta,f0,f0fac,ffac,r,rfac,jacfac,w11,w12,bp,
     $     bt,b,xm,dx,rholow,rhohigh,eta,delpsi,q
      REAL(r8), DIMENSION(3,3) :: v
      REAL(r8), DIMENSION(2, mpsi+1) :: xdx
      REAL(r8),DIMENSION(:,:),ALLOCATABLE :: x,y,r2,deta
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     fit input to cubic splines and diagnose.
c-----------------------------------------------------------------------
      ALLOCATE(x(0:rz_in%mx,0:rz_in%my),y(0:rz_in%mx,0:rz_in%my))
      ALLOCATE(r2(0:rz_in%mx,0:rz_in%my),deta(0:rz_in%mx,0:rz_in%my))
      direct_flag=.FALSE.
      sq_in%fs(:,4)=SQRT(sq_in%xs)
      sq_in%name="  sq  "
      sq_in%title=(/"psifac","  f   ","mu0 p ","  q   "," rho  "/)
      CALL spline_fit(sq_in,"extrap")
      rz_in%xs=sq_in%xs
      rz_in%ys=(/(itheta,itheta=0,rz_in%my)/)/REAL(rz_in%my,r8)
      rz_in%name="  rz  "
      rz_in%xtitle="psifac"
      rz_in%ytitle="theta "
      rz_in%title=(/"  r   ","  z   "/)
      CALL inverse_output
c-----------------------------------------------------------------------
c     allocate and define local arrays.
c-----------------------------------------------------------------------
      mx=rz_in%mx
      my=rz_in%my
      x=rz_in%fs(:,:,1)-ro
      y=rz_in%fs(:,:,2)-zo
      r2=x*x+y*y
      ! atan2 cannot take both zeros. <MODIFIED>
      !deta=ATAN2(y,x)/twopi
      DO ipsi=0,rz_in%mx
         DO itheta=0,rz_in%my
            IF (r2(ipsi,itheta) == 0.0) THEN
               deta(ipsi,itheta)=0.0
            ELSE
               deta(ipsi,itheta)=
     $              ATAN2(y(ipsi,itheta),x(ipsi,itheta))/twopi
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     transform input coordinates from cartesian to polar.
c-----------------------------------------------------------------------
      DO ipsi=0,rz_in%mx
         DO itheta=1,rz_in%my
            IF(deta(ipsi,itheta)-deta(ipsi,itheta-1) > .5)
     $           deta(ipsi,itheta)=deta(ipsi,itheta)-1
            IF(deta(ipsi,itheta)-deta(ipsi,itheta-1) < -.5)
     $           deta(ipsi,itheta)=deta(ipsi,itheta)+1
         ENDDO
         WHERE(r2(ipsi,:) > 0)deta(ipsi,:)=deta(ipsi,:)-rz_in%ys
      ENDDO
      deta(0,:)=inverse_extrap(r2(1:me,:),deta(1:me,:),zero)
      rz_in%fs(:,:,1)=r2
      rz_in%fs(:,:,2)=deta
      CALL bicube_fit(rz_in,"extrap","periodic")
c-----------------------------------------------------------------------
c     diagnose rz_in.
c-----------------------------------------------------------------------
      IF(diagnose_rz_in)THEN
         WRITE(*,*)"Diagnose rz_in."

         CALL ascii_open(debug_unit,"xy.out","UNKNOWN")
         CALL bin_open(bin_unit,"xy.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_xy(rz_in,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL ascii_open(debug_unit,"yx.out","UNKNOWN")
         CALL bin_open(bin_unit,"yx.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_yx(rz_in,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL program_stop
     $        ("inverse_run: termination after diagnose rz_in.")
      ENDIF
c-----------------------------------------------------------------------
c     prepare new spline type for surface quantities.
c-----------------------------------------------------------------------
      IF(grid_type == "original" .OR. grid_type == "orig")mpsi=sq_in%mx
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
         sq%xs=sq_in%xs
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
c     prepare new bicube type for coordinates.
c-----------------------------------------------------------------------
      IF(mtheta == 0)mtheta=rz_in%my
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
c-----------------------------------------------------------------------
c     prepare local splines and start loop over new surfaces.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,mtheta,5)
      DO ipsi=0,mpsi
         psifac=rzphi%xs(ipsi)
         CALL spline_eval(sq_in,psifac,0)
         spl%xs=rzphi%ys
c-----------------------------------------------------------------------
c     interpolate to new surface.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
            CALL bicube_eval(rz_in,psifac,theta,1)
            CALL spline_eval(sq_in,psifac,0)
            IF(rz_in%f(1) < 0)CALL program_stop("Invalid extrapolation"
     $           //" near axis, rerun with larger value of psilow")
            rfac=SQRT(rz_in%f(1))
            r=ro+rfac*COS(twopi*(theta+rz_in%f(2)))
            jacfac=rz_in%fx(1)*(1+rz_in%fy(2))-rz_in%fy(1)*rz_in%fx(2)
            w11=(1+rz_in%fy(2))*twopi**2*rfac/jacfac
            w12=-rz_in%fy(1)*pi/(rfac*jacfac)
            bp=psio*SQRT(w11*w11+w12*w12)/r
            bt=sq_in%f(1)/r
            b=SQRT(bp*bp+bt*bt)
            spl%fs(itheta,1)=rz_in%f(1)
            spl%fs(itheta,2)=rz_in%f(2)
            spl%fs(itheta,3)=r*jacfac
            spl%fs(itheta,4)=spl%fs(itheta,3)/(r*r)
            spl%fs(itheta,5)=spl%fs(itheta,3)
     $           *bp**power_bp*b**power_b/r**power_r
         ENDDO
c-----------------------------------------------------------------------
c     fit to cubic splines and integrate.
c-----------------------------------------------------------------------
         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)
         spl%xs=spl%fsi(:,5)/spl%fsi(mtheta,5)
         spl%fs(:,2)=spl%fs(:,2)+rzphi%ys-spl%xs
         spl%fs(:,4)=(spl%fs(:,3)/spl%fsi(mtheta,3))
     $        /(spl%fs(:,5)/spl%fsi(mtheta,5))
     $        *spl%fsi(mtheta,3)*twopi*pi
         spl%fs(:,3)=sq_in%f(1)*pi/psio
     $        *(spl%fsi(:,4)-spl%fsi(mtheta,4)*spl%xs)
         CALL spline_fit(spl,"periodic")
c-----------------------------------------------------------------------
c     interpolate to final storage.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
            CALL spline_eval(spl,theta,0)
            rzphi%fs(ipsi,itheta,:)=spl%f(1:4)
         ENDDO
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
         sq%fs(ipsi,1)=sq_in%f(1)*twopi
         sq%fs(ipsi,2)=sq_in%f(2)
         sq%fs(ipsi,3)=spl%fsi(mtheta,3)*twopi*pi
         sq%fs(ipsi,4)=spl%fsi(mtheta,4)*sq%fs(ipsi,1)/(2*twopi*psio)
      ENDDO
      CALL spline_fit(sq,"extrap")
      q0=sq%fs(0,4)-sq%fs1(0,4)*sq%xs(0)
      IF(newq0 == -1)newq0=-q0
c-----------------------------------------------------------------------
c     revise q profile.
c-----------------------------------------------------------------------
      IF(newq0 /= 0)THEN
         f0=sq%fs(0,1)-sq%fs1(0,1)*sq%xs(0)
         f0fac=f0**2*((newq0/q0)**2-one)
         q0=newq0
         DO ipsi=0,mpsi
            ffac=SQRT(1+f0fac/sq%fs(ipsi,1)**2)*SIGN(one,newq0)
            sq%fs(ipsi,1)=sq%fs(ipsi,1)*ffac
            sq%fs(ipsi,4)=sq%fs(ipsi,4)*ffac
            rzphi%fs(ipsi,:,3)=rzphi%fs(ipsi,:,3)*ffac
         ENDDO
         CALL spline_fit(sq,"extrap")
      ENDIF
      qa=sq%fs(mpsi,4)+sq%fs1(mpsi,4)*(one-sq%xs(mpsi))
c-----------------------------------------------------------------------
c     fit rzphi to bicubic splines and deallocate input arrays.
c-----------------------------------------------------------------------
      IF(power_flag)rzphi%xpower(1,:)=(/1._r8,0._r8,.5_r8,0._r8/)
      CALL bicube_fit(rzphi,"extrap","periodic")
      CALL spline_dealloc(sq_in)
      CALL bicube_dealloc(rz_in)
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     diagnose rzphi.
c-----------------------------------------------------------------------
      IF(diagnose_rzphi)THEN
         WRITE(*,*)"Diagnose rzphi."

         CALL ascii_open(debug_unit,"xy.out","UNKNOWN")
         CALL bin_open(bin_unit,"xy.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_xy(rzphi,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL ascii_open(debug_unit,"yx.out","UNKNOWN")
         CALL bin_open(bin_unit,"yx.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_yx(rzphi,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL program_stop
     $        ("inverse_run: termination after diagnose rzphi.")
      ENDIF

      DEALLOCATE(x,y,r2,deta)
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
      END SUBROUTINE inverse_run
c-----------------------------------------------------------------------
c     subprogram 2. inverse_output.
c     diagnoses input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inverse_output
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(.NOT. (out_eq_1d .OR. bin_eq_1d .OR. out_eq_2d .OR. bin_eq_2d))
     $     RETURN
      IF(out_eq_1d .OR. out_eq_2d)
     $     CALL ascii_open(out_2d_unit,"input.out","UNKNOWN")
      IF(interp)CALL bicube_fit(rz_in,"extrap","periodic")
c-----------------------------------------------------------------------
c     diagnose 1d output.
c-----------------------------------------------------------------------
      IF(bin_eq_1d)CALL bin_open(bin_2d_unit,"sq_in.bin","UNKNOWN",
     $     "REWIND","none")
      IF(out_eq_1d)WRITE(out_2d_unit,'(a)')"input surface quantities:"
      CALL spline_write1(sq_in,out_eq_1d,bin_eq_1d,
     $     out_2d_unit,bin_2d_unit,interp)
      IF(bin_eq_1d)CALL bin_close(bin_2d_unit)
c-----------------------------------------------------------------------
c     diagnose 2d output.
c-----------------------------------------------------------------------
      IF(bin_eq_2d)CALL bin_open(bin_2d_unit,"rz_in.bin","UNKNOWN",
     $     "REWIND","none")
      IF(out_eq_2d)WRITE(out_2d_unit,'(1x,a/)')"input coordinates:"
      CALL bicube_write_yx(rz_in,out_eq_2d,bin_eq_2d,
     $     out_2d_unit,bin_2d_unit,interp)
      IF(bin_eq_2d)CALL bin_close(bin_2d_unit)
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(out_eq_1d .OR. out_eq_2d)CALL ascii_close(out_2d_unit)
      IF(input_only)CALL program_stop("Termination by inverse_output.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inverse_output
c-----------------------------------------------------------------------
c     subprogram 3. inverse_extrap.
c     extrapolates function to f(x).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION inverse_extrap(xx,ff,x) RESULT(f)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: xx,ff
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(SIZE(ff,2)) :: f
      
      INTEGER :: i,j,m
      REAL(r8), DIMENSION(SIZE(ff,2)) :: term
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      m=SIZE(xx,1)
      f=0
      DO i=1,m
         term=ff(i,:)
         DO j=1,m
            IF(j == i)CYCLE
            term=term*(x-xx(j,:))/(xx(i,:)-xx(j,:))
         ENDDO
         f=f+term
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION inverse_extrap
c-----------------------------------------------------------------------
c     subprogram 4. inverse_chease4_run.
c     gets chease equilibrium data and use in the computation directly.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inverse_chease4_run(jaceq)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: jaceq
      LOGICAL, PARAMETER :: interp=.FALSE.,
     $     diagnose_rz_in=.FALSE.,diagnose_rzphi=.FALSE.
      INTEGER :: ipsi,itheta,mx,my
      INTEGER, PARAMETER :: me=3
      REAL(r8) :: psifac,theta,f0,f0fac,ffac,r,rfac,jacfac,w11,w12,bp,
     $     bt,b,xm,dx,rholow,rhohigh,eta,delpsi,q
      REAL(r8), DIMENSION(3,3) :: v
      REAL(r8),DIMENSION(:,:),ALLOCATABLE :: x,y,r2,deta
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     fit input to cubic splines and diagnose.
c-----------------------------------------------------------------------
      ALLOCATE(x(0:rz_in%mx,0:rz_in%my),y(0:rz_in%mx,0:rz_in%my))
      ALLOCATE(r2(0:rz_in%mx,0:rz_in%my),deta(0:rz_in%mx,0:rz_in%my))
      direct_flag=.FALSE.
      sq_in%fs(:,4)=SQRT(sq_in%xs)
      sq_in%name="  sq  "
      sq_in%title=(/"psifac","  f   ","mu0 p ","  q   "," rho  "/)
      CALL spline_fit(sq_in,"extrap")
      rz_in%xs=sq_in%xs
      rz_in%ys=(/(itheta,itheta=0,rz_in%my)/)/REAL(rz_in%my,r8)
      rz_in%name="  rz  "
      rz_in%xtitle="psifac"
      rz_in%ytitle="theta "
      rz_in%title=(/"  r   ","  z   "/)
      CALL inverse_output
c-----------------------------------------------------------------------
c     allocate and define local arrays.
c-----------------------------------------------------------------------
      mx=rz_in%mx
      my=rz_in%my
      x=rz_in%fs(:,:,1)-ro
      y=rz_in%fs(:,:,2)-zo
      r2=x*x+y*y
      deta=ATAN2(y,x)/twopi
c-----------------------------------------------------------------------
c     transform input coordinates from cartesian to polar.
c-----------------------------------------------------------------------
      DO ipsi=0,rz_in%mx
         DO itheta=1,rz_in%my
            IF(deta(ipsi,itheta)-deta(ipsi,itheta-1) > .5)
     $           deta(ipsi,itheta)=deta(ipsi,itheta)-1
            IF(deta(ipsi,itheta)-deta(ipsi,itheta-1) < -.5)
     $           deta(ipsi,itheta)=deta(ipsi,itheta)+1
         ENDDO
         WHERE(r2(ipsi,:) > 0)deta(ipsi,:)=deta(ipsi,:)-rz_in%ys
      ENDDO
ccc      deta(0,:)=inverse_extrap(r2(1:me,:),deta(1:me,:),zero)
      rz_in%fs(:,:,1)=r2
      rz_in%fs(:,:,2)=deta
      CALL bicube_fit(rz_in,"extrap","periodic")
c-----------------------------------------------------------------------
c     prepare new spline type for surface quantities.
c-----------------------------------------------------------------------
      IF(grid_type == "original" .OR. grid_type == "orig")mpsi=sq_in%mx
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
      CASE("rho")
         sq%xs=psihigh*(/(ipsi**2,ipsi=1,mpsi+1)/)/(mpsi+1)**2
      CASE("original","orig")
         sq%xs=sq_in%xs
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
c     prepare new bicube type for coordinates.
c-----------------------------------------------------------------------
      IF(mtheta == 0) mtheta=rz_in%my
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
c-----------------------------------------------------------------------
c     prepare local splines and start loop over new surfaces.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,mtheta,6)
      DO ipsi=0,mpsi
         psifac=rzphi%xs(ipsi)
         CALL spline_eval(sq_in,psifac,0)
         spl%xs=rzphi%ys
c-----------------------------------------------------------------------
c     interpolate to new surface.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
            CALL bicube_eval(rz_in,psifac,theta,1)
            CALL spline_eval(sq_in,psifac,0)
            IF(rz_in%f(1) < 0)CALL program_stop("Invalid extrapolation"
     $           //" near axis, rerun with larger value of psilow")
            rfac=SQRT(rz_in%f(1))
            r=ro+rfac*COS(twopi*(theta+rz_in%f(2)))
            r=rz_in%f(5)
            jacfac=rz_in%fx(1)*(1+rz_in%fy(2))-rz_in%fy(1)*rz_in%fx(2)
ccc            jacfac=rz_in%f(4)
            w11=(1+rz_in%fy(2))*twopi**2*rfac/jacfac
            w12=-rz_in%fy(1)*pi/(rfac*jacfac)
            bp=psio*SQRT(w11*w11+w12*w12)/r
            bt=sq_in%f(1)/r
            b=SQRT(bp*bp+bt*bt)
            spl%fs(itheta,1)=rz_in%f(1)
            spl%fs(itheta,2)=rz_in%f(2)
            spl%fs(itheta,3)=r*jacfac
            spl%fs(itheta,4)=spl%fs(itheta,3)/(r*r)
            spl%fs(itheta,5)=spl%fs(itheta,3)
     $           *bp**power_bp*b**power_b/r**power_r
            spl%fs(itheta,6)=rz_in%f(3)
         ENDDO
c-----------------------------------------------------------------------
c     fit to cubic splines and integrate.
c-----------------------------------------------------------------------
         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)
         spl%xs=spl%fsi(:,5)/spl%fsi(mtheta,5)
         spl%fs(:,2)=spl%fs(:,2)
ccc      solve for jac
         spl%fs(:,4)=spl%fs(:,6)
ccc         spl%fs(:,4)=(spl%fs(:,3)/spl%fsi(mtheta,3))
ccc     $        /(spl%fs(:,5)/spl%fsi(mtheta,5))
ccc     $        *spl%fsi(mtheta,3)*twopi*pi
ccc      dphi=0 for PEST
         spl%fs(:,3)=0
ccc         spl%fs(:,3)=sq_in%f(1)*pi/psio
ccc     $        *(spl%fsi(:,4)-spl%fsi(mtheta,4)*spl%xs)
         CALL spline_fit(spl,"periodic")
c-----------------------------------------------------------------------
c     interpolate to final storage.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
            CALL spline_eval(spl,theta,0)
            rzphi%fs(ipsi,itheta,:)=spl%f(1:4)
         ENDDO
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
         sq%fs(ipsi,1)=sq_in%f(1)*twopi
         sq%fs(ipsi,2)=sq_in%f(2)
         sq%fs(ipsi,3)=spl%fsi(mtheta,3)*twopi*pi
         sq%fs(ipsi,4)=sq_in%f(3)
      ENDDO
      CALL spline_fit(sq,"extrap")
      q0=sq%fs(0,4)-sq%fs1(0,4)*sq%xs(0)
      IF(newq0 == -1)newq0=-q0
c-----------------------------------------------------------------------
c     revise q profile.
c-----------------------------------------------------------------------
      IF(newq0 /= 0)THEN
         f0=sq%fs(0,1)-sq%fs1(0,1)*sq%xs(0)
         f0fac=f0**2*((newq0/q0)**2-one)
         q0=newq0
         DO ipsi=0,mpsi
            ffac=SQRT(1+f0fac/sq%fs(ipsi,1)**2)*SIGN(one,newq0)
            sq%fs(ipsi,1)=sq%fs(ipsi,1)*ffac
            sq%fs(ipsi,4)=sq%fs(ipsi,4)*ffac
            rzphi%fs(ipsi,:,3)=rzphi%fs(ipsi,:,3)*ffac
         ENDDO
         CALL spline_fit(sq,"extrap")
      ENDIF
      qa=sq%fs(mpsi,4)+sq%fs1(mpsi,4)*(one-sq%xs(mpsi))
c-----------------------------------------------------------------------
c     fit rzphi to bicubic splines and deallocate input arrays.
c-----------------------------------------------------------------------
      IF(power_flag)rzphi%xpower(1,:)=(/1._r8,0._r8,.5_r8,0._r8/)
      CALL bicube_fit(rzphi,"extrap","periodic")
      CALL spline_dealloc(sq_in)
      CALL bicube_dealloc(rz_in)
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     diagnose rzphi.
c-----------------------------------------------------------------------
      IF(diagnose_rzphi)THEN
         WRITE(*,*)"Diagnose rzphi."

         CALL ascii_open(debug_unit,"xy.out","UNKNOWN")
         CALL bin_open(bin_unit,"xy.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_xy(rzphi,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL ascii_open(debug_unit,"yx.out","UNKNOWN")
         CALL bin_open(bin_unit,"yx.bin","UNKNOWN","REWIND","none")
         CALL bicube_write_yx(rzphi,.TRUE.,.TRUE.,
     $        debug_unit,bin_unit,interp)
         CALL bin_close(bin_unit)
         CALL ascii_close(debug_unit)

         CALL program_stop
     $        ("inverse_run: termination after diagnose rzphi.")
      ENDIF

      DEALLOCATE(x,y,r2,deta)
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
      END SUBROUTINE inverse_chease4_run
      END MODULE inverse_mod
