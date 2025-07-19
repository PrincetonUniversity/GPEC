c-----------------------------------------------------------------------
c     file mercier.f.
c     computes mercier criterion and related quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. rdcon_mercier_mod.
c     1. mercier_scan.
c-----------------------------------------------------------------------
c     subprogram 0. rdcon_mercier_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rdcon_mercier_mod
      USE rdcon_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. mercier_scan.
c     evaluates mercier criterion and related quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE mercier_scan

      INTEGER :: ipsi,itheta
      REAL(r8) :: bsq,chi1,di,dpsisq,eta,h,jac,p1,psifac,q,q1,r,
     $     rfac,term,theta,twopif,v1,v2,v21,v22,v23,v33,
     $     bt,f1,Dnc,Dnc_prefac,Wc_prefac,eps_loc,Jpara,ftr,Mloc,
     $     Hbs_prefac,Jboot_dot_B,mufrac,taua_prefac,taur_prefac
      REAL(r8), DIMENSION(:), POINTER :: avg
      TYPE(spline_type), TARGET :: ff
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      IF(compute_MRE_terms)THEN
         CALL spline_alloc(ff,mtheta,21)
      ELSE
         CALL spline_alloc(ff,mtheta,5)
      ENDIF
      ff%xs=rzphi%ys
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         twopif=sq%fs(ipsi,1)
         f1=sq%fs1(ipsi,1)/twopi
         p1=sq%fs1(ipsi,2)
         v1=sq%fs(ipsi,3)
         v2=sq%fs1(ipsi,3)
         q=sq%fs(ipsi,4)
         q1=sq%fs1(ipsi,4)
         chi1=twopi*psio
c-----------------------------------------------------------------------
c     evaluate coordinates and jacobian.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            theta=rzphi%ys(itheta)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta+rzphi%f(2))
            r=ro+rfac*COS(eta)
            jac=rzphi%f(4)
            bt=twopif/(twopi*r) !This is toroidal B field
c-----------------------------------------------------------------------
c     evaluate other local quantities.
c-----------------------------------------------------------------------
            v21=rzphi%fy(1)/(2*rfac*jac)
            v22=(1+rzphi%fy(2))*twopi*rfac/jac
            v23=rzphi%fy(3)*r/jac
            v33=twopi*r/jac
            bsq=chi1**2*(v21**2+v22**2+(v23+q*v33)**2) !|B|^2
            dpsisq=(twopi*r)**2*(v21**2+v22**2)        !|nabla psi|^2
c-----------------------------------------------------------------------
c     evaluate integrands.
c-----------------------------------------------------------------------
            ff%fs(itheta,1)=bsq/dpsisq
            ff%fs(itheta,2)=1/dpsisq
            ff%fs(itheta,3)=1/bsq
            ff%fs(itheta,4)=1/(bsq*dpsisq)
            ff%fs(itheta,5)=bsq
            IF(compute_MRE_terms)THEN
               ff%fs(itheta,6)=dpsisq/bsq
               ff%fs(itheta,7)=dpsisq
               ff%fs(itheta,8)=SQRT(bsq)        ! |B|
               ff%fs(itheta,9)=bt               ! |B_toroidal|
               ff%fs(itheta,10)=SQRT(bsq-bt**2) ! |B_poloidal|
               ff%fs(itheta,11)=rfac            ! minor radius
               ff%fs(itheta,12)=r               ! major radius
               ff%fs(itheta,13)=1.d0/r          ! 1/major radius
               ff%fs(itheta,14)=dpsisq/(r**2)
               ff%fs(itheta,15)=1.d0/(r**2)
               ff%fs(itheta,16)=dpsisq/(r**2*SQRT(bsq))
               ff%fs(itheta,17)=1.d0/(r**2*SQRT(bsq))
               ff%fs(itheta,18)=1.d0/SQRT(bsq)
               ff%fs(itheta,19)=SQRT(dpsisq)/r
               ff%fs(itheta,20)=r**2*v1/jac !overbar{R^2}     (Hegna 1999)
               ff%fs(itheta,21)=r**2        !avg{R^2} ~ [m^2] (Hegna 1999)
            ENDIF
            ff%fs(itheta,:)=ff%fs(itheta,:)*jac/v1
         ENDDO
c-----------------------------------------------------------------------
c     integrate quantities with respect to theta.
c-----------------------------------------------------------------------
         CALL spline_fit(ff,"periodic")
         CALL spline_int(ff)
         avg => ff%fsi(mtheta,:)
c-----------------------------------------------------------------------
c     evaluate mercier criterion and related quantities.
c-----------------------------------------------------------------------
         term=twopif*p1*v1/(q1*chi1**3)*avg(2)
         di=-.25+term*(1-term)+p1*(v1/(q1*chi1**2))**2*avg(1)
     $        *(p1*(avg(3)+(twopif/chi1)**2*avg(4))-v2/v1)
         h=twopif*p1*v1/(q1*chi1**3)*(avg(2)-avg(1)/avg(5))
         locstab%fs(ipsi,1)=di*locstab%xs(ipsi)
         locstab%fs(ipsi,2)=(di+(h-0.5)**2)*locstab%xs(ipsi)
         locstab%fs(ipsi,3)=h
c-----------------------------------------------------------------------
c     MRE term calculations
c-----------------------------------------------------------------------
         IF(compute_MRE_terms)THEN
c-----------------------------------------------------------------------
c     computes mass factor M from Glasser 2016 eq. A8, as in resist.f
c-----------------------------------------------------------------------
            M=avg(1)*(avg(6)+(twopif/chi1)**2*(avg(3)-1/avg(5)))
c-----------------------------------------------------------------------
c     computes geometric prefactors of Alfven and resistive time scales
c     from Glasser 2016 eqs. A12, A13
c-----------------------------------------------------------------------
            taua_prefac=SQRT(M*mu0)/ABS(twopi*q1*chi1/v1) 
            !to get taua, multiply by local sqrt(rho) and divide by 
            !toroidal mode number nn (see resist.f)
            taur_prefac=avg(1)/avg(5)*mu0
            !to get taur, divide by local resistivity (see resist.f)
c-----------------------------------------------------------------------
c     simple estimates of trapped fraction from Sauter 2002:
c     <https://infoscience.epfl.ch/server/api/core/bitstreams/c42baba0-9
c     909-4f21-978a-0d0c2646c3ad/content>
c-----------------------------------------------------------------------
            eps_loc=avg(11)/avg(12)
            ftr=1.d0-(1-eps_loc)**2/
     $       (SQRT(1-eps_loc**2)*(1.d0+1.46d0*SQRT(eps_loc)))
            ftr=MIN(ftr,1.0d0)
         ENDIF
120      FORMAT(19(E30.15,1X))
      ENDDO
      CALL spline_dealloc(ff)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mercier_scan
      END MODULE rdcon_mercier_mod
