c-----------------------------------------------------------------------
c     file mercier.f.
c     computes mercier criterion and related quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. mercier_mod.
c     1. mercier_scan.
c-----------------------------------------------------------------------
c     subprogram 0. mercier_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE mercier_mod
      USE dcon_mod
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
     $     rfac,term,theta,twopif,v1,v2,v21,v22,v23,v33
      REAL(r8), DIMENSION(:), POINTER :: avg
      TYPE(spline_type), TARGET :: ff
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      CALL spline_alloc(ff,mtheta,5)
      ff%xs=rzphi%ys
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         twopif=sq%fs(ipsi,1)
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
c-----------------------------------------------------------------------
c     evaluate other local quantities.
c-----------------------------------------------------------------------
            v21=rzphi%fy(1)/(2*rfac*jac)
            v22=(1+rzphi%fy(2))*twopi*rfac/jac
            v23=rzphi%fy(3)*r/jac
            v33=twopi*r/jac
            bsq=chi1**2*(v21**2+v22**2+(v23+q*v33)**2)
            dpsisq=(twopi*r)**2*(v21**2+v22**2)
c-----------------------------------------------------------------------
c     evaluate integrands.
c-----------------------------------------------------------------------
            ff%fs(itheta,1)=bsq/dpsisq
            ff%fs(itheta,2)=1/dpsisq
            ff%fs(itheta,3)=1/bsq
            ff%fs(itheta,4)=1/(bsq*dpsisq)
            ff%fs(itheta,5)=bsq
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
      ENDDO
      CALL spline_dealloc(ff)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mercier_scan
      END MODULE mercier_mod
