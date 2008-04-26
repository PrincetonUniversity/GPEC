c-----------------------------------------------------------------------
c     file field.f.
c     computes equilibrium fields.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. field_mod.
c     1. field_1d.
c     2. field_2d.
c     3. field_get_jac.
c     4. field_gc_in.
c     5. field_gc_out.
c-----------------------------------------------------------------------
c     subprogram 0. field_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE field_mod
      USE orbit_mod
      IMPLICIT NONE

      LOGICAL :: field_out=.FALSE.,field_bin=.FALSE.,
     $     interp_flag=.FALSE.,prefit=.TRUE.
      TYPE(bicube_type) :: metric,gc_infield,gc_outfield
      TYPE(spline_type) :: avec

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. field_1d.
c     computes vector potential components.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_1d
c-----------------------------------------------------------------------
c     prepare cubic splines.
c-----------------------------------------------------------------------
      CALL spline_alloc(avec,mpsi,2)
      avec%xs=rzphi%xs
      avec%name=" avec "
      avec%title=(/" psi  "," -Phi "," chi  ","1/V1^2"/)
c-----------------------------------------------------------------------
c     compute functions.
c-----------------------------------------------------------------------
      avec%fs(:,2)=-twopi*psio
      avec%fs(:,1)=-avec%fs(:,2)*sq%fs(:,4)
c-----------------------------------------------------------------------
c     fit and integrate cubic splnies.
c-----------------------------------------------------------------------
      CALL spline_fit(avec,"extrap")
      CALL spline_int(avec)
      avec%fs(:,1:2)=avec%fsi(:,1:2)
      CALL spline_fit(avec,"extrap")
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"avec.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"avec.bin","UNKNOWN",
     $        "REWIND","none")
         CALL spline_write1(avec,field_out,field_bin,
     $        field_out_unit,field_bin_unit,.TRUE.)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_1d
c-----------------------------------------------------------------------
c     subprogram 2. field_2d.
c     computes metric tensor components.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_2d

      INTEGER :: ipsi,itheta
      REAL(r8) :: cosfac,sinfac,rfac,rfacx,rfacy,eta,etax,etay,
     $     phix,phiy,r,theta
c-----------------------------------------------------------------------
c     prepare bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_alloc(metric,mpsi,mtheta,7)
      metric%xs=rzphi%xs
      metric%ys=rzphi%ys
      metric%name="metric"
      metric%xtitle=" psi  "
      metric%ytitle="theta "
      metric%title=(/" j11  "," j12  "," j21  "," j22  ",
     $     " j31  "," j32  "," j33  "/)
c-----------------------------------------------------------------------
c     begin loop over nodes.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
c-----------------------------------------------------------------------
c     compute cylindrical coordinates and their derivatives.
c-----------------------------------------------------------------------
            rfac=SQRT(rzphi%fs(ipsi,itheta,1))
            rfacx=rzphi%fsx(ipsi,itheta,1)/(2*rfac)
            rfacy=rzphi%fsy(ipsi,itheta,1)/(2*rfac)
            eta=twopi*(theta+rzphi%fs(ipsi,itheta,2))
            etax=twopi*rzphi%fsx(ipsi,itheta,2)
            etay=twopi*(1+rzphi%fsy(ipsi,itheta,2))
            phix=rzphi%fsx(ipsi,itheta,3)
            phiy=rzphi%fsy(ipsi,itheta,3)
            cosfac=COS(eta)
            sinfac=SIN(eta)
            r=ro+rfac*cosfac
c-----------------------------------------------------------------------
c     set up metric tensor.
c-----------------------------------------------------------------------
            metric%fs(ipsi,itheta,1)=rfacx*cosfac-rfac*etax*sinfac
            metric%fs(ipsi,itheta,2)=rfacy*cosfac-rfac*etay*sinfac
            metric%fs(ipsi,itheta,3)=rfacx*sinfac+rfac*etax*cosfac
            metric%fs(ipsi,itheta,4)=rfacy*sinfac+rfac*etay*cosfac
            metric%fs(ipsi,itheta,5)=phix*r
            metric%fs(ipsi,itheta,6)=phiy*r
            metric%fs(ipsi,itheta,7)=twopi*r
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit to bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_fit(metric,"extrap","periodic")
      IF(prefit)THEN
         CALL bicube_all_getco(metric)
         CALL bicube_all_getco(psi_in)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, xy.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metxy.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metxy.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_xy(metric,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, yx.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metyx.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metyx.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_yx(metric,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_2d
c-----------------------------------------------------------------------
c     subprogram 3. field_get_jac.
c     local evaluation of Jacobian matrix and its derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_get_jac(psi,theta,mode,jac,jacx,jacy)

      REAL(r8), INTENT(IN) :: psi,theta
      INTEGER, INTENT(IN) :: mode
      REAL(r8), DIMENSION(3,3), INTENT(OUT) :: jac,jacx,jacy
c-----------------------------------------------------------------------
c     evaluate splines.
c-----------------------------------------------------------------------
      CALL bicube_eval(metric,psi,theta,mode)
c-----------------------------------------------------------------------
c     set up jac matrix.
c-----------------------------------------------------------------------
      jac(1,1)=metric%f(1)
      jac(1,2)=metric%f(2)
      jac(1,3)=0
      jac(2,1)=metric%f(3)
      jac(2,2)=metric%f(4)
      jac(2,3)=0
      jac(3,1)=metric%f(5)
      jac(3,2)=metric%f(6)
      jac(3,3)=metric%f(7)
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     set up x derivatives of jac.
c-----------------------------------------------------------------------
      jacx(1,1)=metric%fx(1)
      jacx(1,2)=metric%fx(2)
      jacx(1,3)=0
      jacx(2,1)=metric%fx(3)
      jacx(2,2)=metric%fx(4)
      jacx(2,3)=0
      jacx(3,1)=metric%fx(5)
      jacx(3,2)=metric%fx(6)
      jacx(3,3)=metric%fx(7)
c-----------------------------------------------------------------------
c     set up y derivatives of jac.
c-----------------------------------------------------------------------
      jacy(1,1)=metric%fy(1)
      jacy(1,2)=metric%fy(2)
      jacy(1,3)=0
      jacy(2,1)=metric%fy(3)
      jacy(2,2)=metric%fy(4)
      jacy(2,3)=0
      jacy(3,1)=metric%fy(5)
      jacy(3,2)=metric%fy(6)
      jacy(3,3)=metric%fy(7)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_get_jac
c-----------------------------------------------------------------------
c     subprogram 4. field_gc_in.
c     computes gc_infield tensor components.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_gc_in

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ipsi,itheta,iqty
      REAL(r8) :: rfac,eta,r,chi1,jac,q,bmod
      REAL(r8), DIMENSION(3) :: bvec
      REAL(r8), DIMENSION(3,3) :: v
c-----------------------------------------------------------------------
c     prepare bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_alloc(gc_infield,mpsi,mtheta,5)
      gc_infield%xs=rzphi%xs
      gc_infield%ys=rzphi%ys
      gc_infield%name="gc_infield"
      gc_infield%xtitle=" psi  "
      gc_infield%ytitle="theta "
      gc_infield%title=(/"b1fac ","b2fac ","b3fac "," bmod "," jac  "/)
c-----------------------------------------------------------------------
c     begin loop over nodes.
c-----------------------------------------------------------------------
      chi1=twopi*psio
      v=0
      DO ipsi=0,mpsi
         q=sq%fs(ipsi,4)
         DO itheta=0,mtheta
c-----------------------------------------------------------------------
c     compute cylindrical coordinates.
c-----------------------------------------------------------------------
            rfac=SQRT(rzphi%fs(ipsi,itheta,1))
            eta=twopi*(rzphi%ys(itheta)+rzphi%fs(ipsi,itheta,2))
            r=ro+rfac*COS(eta)
            jac=rzphi%fs(ipsi,itheta,4)
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
            v(1,1)=rzphi%fsx(ipsi,itheta,1)/(2*rfac*jac)
            v(1,2)=rzphi%fsx(ipsi,itheta,2)*twopi*rfac/jac
            v(1,3)=rzphi%fsx(ipsi,itheta,3)*r/jac
            v(2,1)=rzphi%fsy(ipsi,itheta,1)/(2*rfac*jac)
            v(2,2)=(1+rzphi%fsy(ipsi,itheta,2))*twopi*rfac/jac
            v(2,3)=rzphi%fsy(ipsi,itheta,3)*r/jac
            v(3,3)=twopi*r/jac
c-----------------------------------------------------------------------
c     compute gc_infield.
c-----------------------------------------------------------------------
            bvec=chi1*(v(2,:)+q*v(3,:))
            bmod=SQRT(SUM(bvec**2))
            gc_infield%fs(ipsi,itheta,1)=SUM(bvec*v(1,:))*jac/bmod
            gc_infield%fs(ipsi,itheta,2)=SUM(bvec*v(2,:))*jac/bmod
            gc_infield%fs(ipsi,itheta,3)=SUM(bvec*v(3,:))*jac/bmod
            gc_infield%fs(ipsi,itheta,4)=bmod
            gc_infield%fs(ipsi,itheta,5)=jac
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit to bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_fit(gc_infield,"extrap","periodic")
      IF(prefit)CALL bicube_all_getco(gc_infield)
c-----------------------------------------------------------------------
c     diagnose arrays.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(field_out_unit,"gc_infield.out","UNKNOWN")
         DO iqty=1,4
            CALL bicube_write_arrays(gc_infield,.TRUE.,
     $           field_out_unit,iqty)
         ENDDO
         CALL ascii_close(field_out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, xy.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metxy.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metxy.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_xy(gc_infield,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(UNIT=field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, yx.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metyx.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metyx.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_yx(gc_infield,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_gc_in
c-----------------------------------------------------------------------
c     subprogram 5. field_gc_out.
c     computes gc_outfield tensor components.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_gc_out

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ir,iz,mr,mz,iqty
      REAL(r8) :: r,dr,br,bz,bt,bmod,f
c-----------------------------------------------------------------------
c     prepare bicubic splines.
c-----------------------------------------------------------------------
      mr=psi_in%mx
      mz=psi_in%my
      CALL bicube_alloc(gc_outfield,mr,mz,4)
      gc_outfield%xs=psi_in%xs
      gc_outfield%ys=psi_in%ys
      gc_outfield%name="gc_outfield"
      gc_outfield%xtitle="  r   "
      gc_outfield%ytitle="  z   "
      gc_outfield%title=(/"b1fac ","b2fac ","b3fac "," bmod "/)
c-----------------------------------------------------------------------
c     compute cubic spline values.
c-----------------------------------------------------------------------
      f=sq%fs(mpsi,1)/twopi
      dr=(rmax-rmin)/mr
      DO ir=0,mr
         r=rmin+dr*ir
         DO iz=0,mz
            br=psi_in%fsy(ir,iz,1)/r
            bz=-psi_in%fsx(ir,iz,1)/r
            bt=f/r
            bmod=SQRT(br*br+bz*bz+bt*bt)
            gc_outfield%fs(ir,iz,1)=br/bmod
            gc_outfield%fs(ir,iz,2)=bz/bmod
            gc_outfield%fs(ir,iz,3)=f/bmod
            gc_outfield%fs(ir,iz,4)=bmod
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit to bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_fit(gc_outfield,"extrap","extrap")
      IF(prefit)CALL bicube_all_getco(gc_outfield)
c-----------------------------------------------------------------------
c     diagnose arrays.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(field_out_unit,"gc_outfield.out","UNKNOWN")
         DO iqty=1,4
            CALL bicube_write_arrays(gc_outfield,.TRUE.,
     $           field_out_unit,iqty)
         ENDDO
         CALL ascii_close(field_out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, xy.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metxy.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metxy.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_xy(gc_outfield,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose, yx.
c-----------------------------------------------------------------------
      IF(field_out .OR. field_bin)THEN
         IF(field_out)CALL ascii_open(field_out_unit,"metyx.out",
     $        "UNKNOWN")
         IF(field_bin)CALL bin_open(field_bin_unit,"metyx.bin",
     $        "UNKNOWN","REWIND","none")
         CALL bicube_write_yx(gc_outfield,.FALSE.,.TRUE.,
     $        field_out_unit,field_bin_unit,interp_flag)
         IF(field_out)CALL ascii_close(field_out_unit)
         IF(field_bin)CALL bin_close(field_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_gc_out
      END MODULE field_mod
