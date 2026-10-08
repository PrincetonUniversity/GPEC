c-----------------------------------------------------------------------
c     file bal.f.
c     determines high-n ideal marginal MHD ballooning stability.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. bal_mod.
c     1. bal_scan.
c     2. bal_prep.
c     3. bal_der.
c     4. bal_int.
c     5. bal_pse.
c-----------------------------------------------------------------------
c     subprogram 0. bal_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE bal_mod
      USE dcon_mod
      IMPLICIT NONE

      INTEGER, PRIVATE :: ipsi
      REAL(r8), PRIVATE :: di,alpha,theta0=0
      REAL(r8), DIMENSION(2,2), PRIVATE :: v0
      TYPE(spline_type), PRIVATE :: bf,bg

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. bal_scan.
c     computes ballooning stability criterion.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bal_scan
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(out_bal1)CALL ascii_open(out_bal1_unit,"bal1.out","UNKNOWN")
      IF(out_bal2)CALL ascii_open(out_bal2_unit,"bal2.out","UNKNOWN")
      IF(bin_bal1)CALL bin_open(bin_bal1_unit,"bal1.bin","UNKNOWN",
     $     "REWIND","none")
      IF(bin_bal2)CALL bin_open(bin_bal2_unit,"bal2.bin","UNKNOWN",
     $     "REWIND","none")
c-----------------------------------------------------------------------
c     test ballooning stability.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         CALL bal_prep
         IF(di <= 0 .AND. di >= -1e4 .AND. sq%xs(ipsi) <= 1)CALL bal_int
         CALL spline_dealloc(bf)
         CALL spline_dealloc(bg)
      ENDDO
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(out_bal1)CALL ascii_close(out_bal1_unit)
      IF(out_bal2)CALL ascii_close(out_bal2_unit)
      IF(bin_bal1)CALL bin_close(bin_bal1_unit)
      IF(bin_bal2)CALL bin_close(bin_bal2_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bal_scan
c-----------------------------------------------------------------------
c     subprogram 2. bal_prep.
c     computes coefficients for the ideal marginal ballooning equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bal_prep

      LOGICAL, PARAMETER :: inverse_flag=.FALSE.
      INTEGER :: itheta,info
      INTEGER, DIMENSION(3) :: ipiv
      REAL(r8) :: chi1,det,eta,fexp,p1,psifac,q,q1,r,rfac,theta,
     $     twopif,v1
      REAL(r8), DIMENSION(2,2) :: c0,d0bar,d,v10
      REAL(r8), DIMENSION(2,2,0:mtheta) :: d1
      REAL(r8), DIMENSION(0:mtheta) :: b1,bsq,dbdb0,dbdb1,dbdb2,
     $     kappan,kappas,jac,jacfac
      REAL(r8), DIMENSION(3,3) :: v,w,temp
      TYPE(spline_type) :: spl0,spl1,spl2,spl3
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
      psifac=sq%xs(ipsi)
      twopif=sq%fs(ipsi,1)
      p1=sq%fs1(ipsi,2)
      v1=sq%fs(ipsi,3)
      q=sq%fs(ipsi,4)
      q1=sq%fs1(ipsi,4)
      chi1=twopi*psio
      v=0
      w=0
c-----------------------------------------------------------------------
c     compute coordinates and jacobian.
c-----------------------------------------------------------------------
      DO itheta=0,mtheta
         CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
         theta=rzphi%ys(itheta)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta+rzphi%f(2))
         r=ro+rfac*COS(eta)
         jac(itheta)=rzphi%f(4)
c-----------------------------------------------------------------------
c     compute contravariant basis vectors.
c-----------------------------------------------------------------------
         v(1,1)=rzphi%fx(1)/(2*rfac*jac(itheta))
         v(1,2)=rzphi%fx(2)*twopi*rfac/jac(itheta)
         v(1,3)=rzphi%fx(3)*r/jac(itheta)
         v(2,1)=rzphi%fy(1)/(2*rfac*jac(itheta))
         v(2,2)=(1+rzphi%fy(2))*twopi*rfac/jac(itheta)
         v(2,3)=rzphi%fy(3)*r/jac(itheta)
         v(3,3)=twopi*r/jac(itheta)
c-----------------------------------------------------------------------
c     compute covariant basis vectors, inverse method.
c-----------------------------------------------------------------------
         IF(inverse_flag)THEN
            temp=TRANSPOSE(v)*jac(itheta)
            w=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
            CALL dgetrf(3,3,temp,3,ipiv,info)
            CALL dgetrs('N',3,3,temp,3,ipiv,w,3,info)
c-----------------------------------------------------------------------
c     compute covariant basis vectors, direct method.
c-----------------------------------------------------------------------
         ELSE
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r
     $           /jac(itheta)
            w(1,2)=-rzphi%fy(1)*pi*r/(rfac*jac(itheta))
            w(2,1)=-rzphi%fx(2)*twopi**2*r*rfac/jac(itheta)
            w(2,2)=rzphi%fx(1)*pi*r/(rfac*jac(itheta))
            w(3,1)=(rzphi%fx(2)*rzphi%fy(3)
     $           -rzphi%fx(3)*(1+rzphi%fy(2)))
     $           *twopi*r*rfac/jac(itheta)
            w(3,2)=(rzphi%fx(3)*rzphi%fy(1)
     $           -rzphi%fx(1)*rzphi%fy(3))
     $           *r/(2*rfac*jac(itheta))
            w(3,3)=1/(twopi*r)
         ENDIF
c-----------------------------------------------------------------------
c     store physical quantities.
c-----------------------------------------------------------------------
         b1(itheta)=SUM((v(2,:)+q*v(3,:))*v(1,:))*chi1
         bsq(itheta)=SUM((v(2,:)+q*v(3,:))**2)*chi1**2
         dbdb0(itheta)=SUM(w(1,:)**2)*q1**2
         dbdb1(itheta)=-SUM((w(3,:)-q*w(2,:))*w(1,:))*2*q1
         dbdb2(itheta)=SUM((w(3,:)-q*w(2,:))**2)
      ENDDO
c-----------------------------------------------------------------------
c     compute curvature terms.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl0,mtheta,2)
      spl0%xs=rzphi%ys
      spl0%fs(:,1)=1/bsq
      spl0%fs(:,2)=jac*b1/bsq
      CALL spline_fit(spl0,"periodic")
      kappas=-spl0%fs1(:,1)*twopif/(2*jac)
      DO itheta=0,mtheta
         CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
         kappan(itheta)=((p1/bsq(itheta)-rzphi%fx(4)/jac(itheta))/chi1
     $        +twopif*q1/(bsq(itheta)*jac(itheta))+spl0%fs1(itheta,2)
     $        /jac(itheta))/2
      ENDDO
      CALL spline_dealloc(spl0)
c-----------------------------------------------------------------------
c     compute coefficients for ballooning equation.
c-----------------------------------------------------------------------
      CALL spline_alloc(bf,mtheta,5)
      bf%xs=rzphi%ys
      jacfac=jac/chi1
      bf%fs(:,1)=dbdb0/(bsq*jacfac)
      bf%fs(:,2)=dbdb1/(2*dbdb0)
      bf%fs(:,3)=(dbdb2-dbdb1**2/(4*dbdb0))/(bsq*jacfac)
      bf%fs(:,4)=2*kappan*p1/chi1*jacfac
      bf%fs(:,5)=-2*kappas*p1/chi1*q1*jacfac
      CALL spline_fit(bf,"periodic")
c-----------------------------------------------------------------------
c     diagnose bf.
c-----------------------------------------------------------------------
      bf%title=(/"theta ","dbdb0 ","dbdb1 ","dbdb2 ","kappas","kappan"/)
      CALL spline_write1(bf,out_bal1,bin_bal1,
     $     out_bal1_unit,bin_bal1_unit,.TRUE.)
c-----------------------------------------------------------------------
c     initialize bg.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl0,mtheta,1)
      spl0%xs=rzphi%ys
      spl0%fs(:,1)=-p1*q1*twopif/(bsq*chi1**2)
      CALL spline_fit(spl0,"periodic")
      CALL spline_int(spl0)
      CALL spline_alloc(bg,mtheta,5)
      bg%xs=rzphi%ys
      bg%fs(:,5)=spl0%fs(:,1)-spl0%fsi(mtheta,1)
      CALL spline_dealloc(spl0)
c-----------------------------------------------------------------------
c     compute c0.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl1,mtheta,4)
      spl1%xs=rzphi%ys
      DO itheta=0,mtheta
         c0(1,1)=.5
         c0(1,2)=1/bf%fs(itheta,1)
         c0(2,1)=-bf%fs(itheta,4)
         c0(2,2)=-.5
c-----------------------------------------------------------------------
c     compute d0.
c-----------------------------------------------------------------------
         spl1%fs(itheta,1)=c0(1,1)+c0(1,2)*bg%fs(itheta,5)
         spl1%fs(itheta,2)=c0(1,2)
         spl1%fs(itheta,3)=c0(2,1)
     $        +(c0(2,2)-c0(1,1)-c0(1,2)*bg%fs(itheta,5))*bg%fs(itheta,5)
         spl1%fs(itheta,4)=-spl1%fs(itheta,1)
c-----------------------------------------------------------------------
c     compute d1.
c-----------------------------------------------------------------------
         fexp=-c0(1,2)*bf%fs(itheta,2)*2
         d1(1,1,itheta)=fexp*bg%fs(itheta,5)
         d1(1,2,itheta)=fexp
         d1(2,1,itheta)=-fexp*bg%fs(itheta,5)**2
         d1(2,2,itheta)=-d1(1,1,itheta)
      ENDDO
c-----------------------------------------------------------------------
c     fit to cubic splines and write output.
c-----------------------------------------------------------------------
      CALL spline_fit(spl1,"periodic")
      CALL spline_int(spl1)
      d0bar(1,1)=spl1%fsi(mtheta,1)
      d0bar(1,2)=spl1%fsi(mtheta,2)
      d0bar(2,1)=spl1%fsi(mtheta,3)
      d0bar(2,2)=spl1%fsi(mtheta,4)
c-----------------------------------------------------------------------
c     compute mercier criterion and alpha.
c-----------------------------------------------------------------------
      di=d0bar(1,1)*d0bar(2,2)-d0bar(1,2)*d0bar(2,1)
      IF(di > 0)THEN
         CALL spline_dealloc(spl1)
         RETURN
      ENDIF
      alpha=SQRT(-di)
c-----------------------------------------------------------------------
c     compute zeroth-order eigenfunctions.
c-----------------------------------------------------------------------
      v0(1,1)=1
      v0(1,2)=1
      v0(2,1)=-(d0bar(1,1)-alpha)/d0bar(1,2)
      v0(2,2)=-(d0bar(1,1)+alpha)/d0bar(1,2)
c-----------------------------------------------------------------------
c     store derivatives for first-order terms.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl2,mtheta,4)
      spl2%xs=rzphi%ys
      DO itheta=0,mtheta
         spl2%fs(itheta,1)
     $        =(spl1%fs(itheta,1)-alpha)*v0(1,1)
     $        +spl1%fs(itheta,2)*v0(2,1)
         spl2%fs(itheta,2)
     $        =spl1%fs(itheta,3)*v0(1,1)
     $        +(spl1%fs(itheta,4)-alpha)*v0(2,1)
         spl2%fs(itheta,3)
     $        =(spl1%fs(itheta,1)+alpha)*v0(1,2)
     $        +spl1%fs(itheta,2)*v0(2,2)
         spl2%fs(itheta,4)
     $        =spl1%fs(itheta,3)*v0(1,2)
     $        +(spl1%fs(itheta,4)+alpha)*v0(2,2)
      ENDDO
c-----------------------------------------------------------------------
c     fit to cubic spines and integrate.
c-----------------------------------------------------------------------
      CALL spline_fit(spl2,"periodic")
      CALL spline_int(spl2)
      spl2%fs=spl2%fsi
c-----------------------------------------------------------------------
c     store derivatives for second-order terms.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl3,mtheta,4)
      spl3%xs=rzphi%ys
      DO itheta=0,mtheta
         spl3%fs(itheta,1)
     $        =(spl1%fs1(itheta,1)+1-alpha)*spl2%fs(itheta,1)
     $        +spl1%fs(itheta,2)*spl2%fs(itheta,2)
     $        +d1(1,1,itheta)*v0(1,1)+d1(1,2,itheta)*v0(2,1)
         spl3%fs(itheta,2)
     $        =spl1%fs(itheta,3)*spl2%fs(itheta,1)
     $        +(spl1%fs(itheta,4)+1-alpha)*spl2%fs(itheta,2)
     $        +d1(2,1,itheta)*v0(1,1)+d1(2,2,itheta)*v0(2,1)
         spl3%fs(itheta,3)
     $        =(spl1%fs(itheta,1)+1+alpha)*spl2%fs(itheta,3)
     $        +spl1%fs(itheta,2)*spl2%fs(itheta,4)
     $        +d1(1,1,itheta)*v0(1,2)+d1(1,2,itheta)*v0(2,2)
         spl3%fs(itheta,4)
     $        =spl1%fs(itheta,3)*spl2%fs(itheta,3)
     $        +(spl1%fs(itheta,4)+1+alpha)*spl2%fs(itheta,4)
     $        +d1(2,1,itheta)*v0(1,2)+d1(2,2,itheta)*v0(2,2)
      ENDDO
c-----------------------------------------------------------------------
c     fit to cubic spines and integrate.
c-----------------------------------------------------------------------
      CALL spline_fit(spl3,"periodic")
      CALL spline_int(spl3)
c-----------------------------------------------------------------------
c     compute first order constants, first eigenfunction.
c-----------------------------------------------------------------------
      d(2,1)=d0bar(2,1)
      d(1,2)=d0bar(1,2)
      d(1,1)=d0bar(1,1)+1-alpha
      d(2,2)=d0bar(2,2)+1-alpha
      det=d(1,1)*d(2,2)-d(1,2)*d(2,1)
      v10(1,1)=(d(1,2)*spl3%fsi(mtheta,2)-d(2,2)*spl3%fsi(mtheta,1))/det
      v10(2,1)=(d(2,1)*spl3%fsi(mtheta,1)-d(1,1)*spl3%fsi(mtheta,2))/det
c-----------------------------------------------------------------------
c     compute first order constants, second eigenfunction.
c-----------------------------------------------------------------------
      d(1,1)=d0bar(1,1)+1+alpha
      d(2,2)=d0bar(2,2)+1+alpha
      det=d(1,1)*d(2,2)-d(1,2)*d(2,1)
      v10(1,2)=(d(1,2)*spl3%fsi(mtheta,4)-d(2,2)*spl3%fsi(mtheta,3))/det
      v10(2,2)=(d(2,1)*spl3%fsi(mtheta,3)-d(1,1)*spl3%fsi(mtheta,4))/det
c-----------------------------------------------------------------------
c     add constants to v1 and fit to cubic splines.
c-----------------------------------------------------------------------
      bg%fs(:,1)=spl2%fs(:,1)+v10(1,1)
      bg%fs(:,2)=spl2%fs(:,2)+v10(2,1)
      bg%fs(:,3)=spl2%fs(:,3)+v10(1,2)
      bg%fs(:,4)=spl2%fs(:,4)+v10(2,2)
      CALL spline_fit(bg,"periodic")
c-----------------------------------------------------------------------
c     diagnose bg.
c-----------------------------------------------------------------------
c     bg%title=(/" bg1  "," bg2  "," bg2  "," bg2  "," bg2  "/)
c     CALL spline_write1(bg,out_bal1,bin_bal1,
c     $     out_bal1_unit,bin_bal1_unit,.FALSE.)
c-----------------------------------------------------------------------
c     deallocate local splines.
c-----------------------------------------------------------------------
      CALL spline_dealloc(spl1)
      CALL spline_dealloc(spl2)
      CALL spline_dealloc(spl3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bal_prep
c-----------------------------------------------------------------------
c     subprogram 3. bal_der
c     differential equations for the marginal ballooning.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bal_der(neq,theta,y,dy)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: theta
      REAL(r8), DIMENSION(neq), INTENT(IN) :: y
      REAL(r8), DIMENSION(neq), INTENT(OUT) :: dy

      REAL(r8) :: dbdb,kappaw,thfac
c-----------------------------------------------------------------------
c     compute coefficients.
c-----------------------------------------------------------------------
      CALL spline_eval(bf,theta,0)
      thfac=theta-theta0
      dbdb=bf%f(1)*(thfac+bf%f(2))**2+bf%f(3)
      kappaw=bf%f(4)+thfac*bf%f(5)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1)=y(2)/dbdb
      dy(2)=-y(1)*kappaw
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bal_der
c-----------------------------------------------------------------------
c     subprogram 4. bal_int
c     integrates the ideal marginal ballooning equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bal_int

      INTEGER, PARAMETER :: neq=2,liw=20,lrw=22+neq*16
      INTEGER, DIMENSION(liw) :: iwork=0
      INTEGER :: iopt,iside,istate,itask,itol,jac,mf
      INTEGER :: istep,nstep=HUGE(0),icross
      REAL(r8), PARAMETER :: tol=1e-5,thmin=1,hu_min=1e-10
      REAL(r8) :: thmax,atol,ca1,ca2,det,rtol,theta,thout
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8), DIMENSION(neq) :: y,dy
      REAL(r8), DIMENSION(2,2) :: u
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"ist",2x,"nqu",1x,"icr",4x,"theta",6x,"dtheta",4x,
     $     "asinh y1",2x,"asinh ca1",2x,"asinh ca2"/)
 20   FORMAT(i5,2i4,1x,1p,5e11.3)
c-----------------------------------------------------------------------
c     initialize independent variables.
c-----------------------------------------------------------------------
      thmax=SQRT(MAXVAL(ABS(bf%fs(:,3)/bf%fs(:,1))))*10*thmax0
      thmax=MIN(thmax,100._r8)
      theta=-thmax
      thout=thmax
      istep=0
c-----------------------------------------------------------------------
c     compute equilibrium quantities and asymptotic solutions.
c-----------------------------------------------------------------------
      CALL bal_pse(theta,u)
      y(1)=u(1,2)*sinh(1.)
      y(2)=u(2,2)*sinh(1.)
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      icross=0
      iside=1
      istate=1
      itask=5
      mf=10
      rwork(1)=thout
      iopt=1
      itol=1
      rtol=tol
      atol=tol*tol
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      IF(out_bal2)WRITE(out_bal2_unit,10)
c-----------------------------------------------------------------------
c     check for crossings.
c-----------------------------------------------------------------------
      DO
         IF(iside*y(1) < 0)THEN
            icross=icross+1
            iside=-iside
         ENDIF
c-----------------------------------------------------------------------
c     compute and print asymptotic coefficient.
c-----------------------------------------------------------------------
         IF(out_bal2 .OR. bin_bal2)THEN
            IF(ABS(theta) > thmin)THEN
               CALL bal_pse(theta,u)
               det=u(1,1)*u(2,2)-u(1,2)*u(2,1)
               ca1=(y(1)*u(2,2)-y(2)*u(1,2))/det
               ca2=(y(2)*u(1,1)-y(1)*u(2,1))/det
               CALL bal_der(neq,theta,y,dy)
               IF(out_bal2)
     $              WRITE(out_bal2_unit,20)istep,iwork(14),icross,
     $              theta,rwork(11),asinh(y(1)),asinh(ca1),asinh(ca2)
               IF(bin_bal2)
     $              WRITE(bin_bal2_unit)REAL(theta,4),
     $              REAL(asinh(y(1)),4),REAL(asinh(y(2)),4),
     $              REAL(y(2)/dy(1),4),REAL(dy(2)/y(1),4),
     $              REAL(asinh(ca1),4),REAL(asinh(ca2),4),REAL(u,4)
            ELSE
               ca1=0
               ca2=0
               det=0
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     advance.
c-----------------------------------------------------------------------
         IF(istep >= nstep .OR. theta >= thmax .OR. istate < 0)EXIT
         istep=istep+1
         CALL lsode(bal_der,neq,y,theta,thout,itol,rtol,atol,itask,
     $        istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         IF(rwork(11) < hu_min)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     write trailer.
c-----------------------------------------------------------------------
      IF(out_bal2)WRITE(out_bal2_unit,10)
      IF(bin_bal2)WRITE(bin_bal2_unit)
c-----------------------------------------------------------------------
c     compute and store final asymptotic coefficients.
c-----------------------------------------------------------------------
      IF(theta >= thmax)THEN
         CALL bal_pse(theta,u)
         det=u(1,1)*u(2,2)-u(1,2)*u(2,1)
         ca1=(y(1)*u(2,2)-y(2)*u(1,2))/det
         ca2=(y(2)*u(1,1)-y(1)*u(2,1))/det
         locstab%fs(ipsi,4:5)=(/ca1,ca2/)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bal_int
c-----------------------------------------------------------------------
c     subprogram 5. bal_pse.
c     computes alpha series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bal_pse(theta,u)

      REAL(r8), INTENT(IN) :: theta
      REAL(r8), INTENT(OUT), DIMENSION(2,2) :: u

      REAL(r8) :: dtheta,thfac
      REAL(r8), DIMENSION(2,2) :: v
c-----------------------------------------------------------------------
c     compute alpha series.
c-----------------------------------------------------------------------
      CALL spline_eval(bg,theta,0)
      dtheta=theta-theta0
      v(1,1)=v0(1,1)+bg%f(1)/dtheta
      v(2,1)=v0(2,1)+bg%f(2)/dtheta
      v(1,2)=v0(1,2)+bg%f(3)/dtheta
      v(2,2)=v0(2,2)+bg%f(4)/dtheta
c-----------------------------------------------------------------------
c     compute transformed solutions.
c-----------------------------------------------------------------------
      u(1,1)=v(1,1)
      u(1,2)=v(1,2)
      u(2,1)=bg%f(5)*v(1,1)+v(2,1)
      u(2,2)=bg%f(5)*v(1,2)+v(2,2)
c-----------------------------------------------------------------------
c     multiply by alphas.
c-----------------------------------------------------------------------
      thfac=ABS(dtheta)**(alpha+.5)/dtheta
      u(1,1)=u(1,1)*thfac
      u(2,1)=u(2,1)*thfac*dtheta
      u(1,2)=u(1,2)/(thfac*dtheta)
      u(2,2)=u(2,2)/thfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bal_pse
      END MODULE bal_mod
