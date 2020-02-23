c-----------------------------------------------------------------------
c     file resist.f.
c     computes resistive surface quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. resist_mod.
c     1. resist_eval.
c-----------------------------------------------------------------------
c     subprogram 0. resist_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE resist_mod
      USE dcon_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. resist_eval.
c     computes resistive surface quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_eval(sing)

      TYPE(sing_type), INTENT(INOUT) :: sing

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: itheta
      REAL(r8) :: bsq,chi1,dpsisq,eta,jac,lambda,p,p1,psifac,q,q1,r,
     $     rfac,taue,theta,twopif,v1,v2,v21,v22,v23,v33
      REAL(r8), PARAMETER :: ne=1e14,te=3e3,mi=2*mp,gamma=5._r8/3._r8
      REAL(r8), DIMENSION(:), POINTER :: avg
      TYPE(spline_type), TARGET :: ff
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(a,1p,e10.3)
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      CALL spline_alloc(ff,mtheta,6)
      ff%xs=rzphi%ys
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
      psifac=sing%psifac
      CALL spline_eval(sq,psifac,1)
      twopif=sq%f(1)
      p=sq%f(2)
      p1=sq%f1(2)
      v1=sq%f(3)
      v2=sq%f1(3)
      q=sq%f(4)
      q1=sq%f1(4)
      chi1=twopi*psio
c-----------------------------------------------------------------------
c     evaluate coordinates and jacobian.
c-----------------------------------------------------------------------
      DO itheta=0,mtheta
         theta=rzphi%ys(itheta)
         CALL bicube_eval(rzphi,psifac,theta,1)
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
         ff%fs(itheta,6)=dpsisq/bsq
         ff%fs(itheta,:)=ff%fs(itheta,:)*jac/v1
      ENDDO
c-----------------------------------------------------------------------
c     integrate quantities with respect to theta.
c-----------------------------------------------------------------------
      CALL spline_fit(ff,"periodic")
      CALL spline_int(ff)
      avg => ff%fsi(mtheta,:)
c-----------------------------------------------------------------------
c     compute resistive curvature terms E, F, and H.
c-----------------------------------------------------------------------
      sing%restype%e=p1*v1/(q1*chi1**2)**2*avg(1)
     $     *(twopif*q1*chi1/avg(5)-v2)
      sing%restype%f=(p1*v1/(q1*chi1**2))**2*(avg(1)*avg(3)
     $     +(twopif/chi1)**2*(avg(1)*avg(4)-avg(2)**2))
      sing%restype%h=twopif*p1*v1/(q1*chi1**3)*(avg(2)-avg(1)/avg(5))
c-----------------------------------------------------------------------
c     computes mass factor and related terms, M, G, and K.
c-----------------------------------------------------------------------
      sing%restype%m=avg(1)*(avg(6)+(twopif/chi1)**2*(avg(3)-1/avg(5)))
      sing%restype%g=avg(5)/(sing%restype%m*gamma*p)
      sing%restype%k=(q1*chi1**2/(p1*v1))**2
     $     *avg(5)/(sing%restype%m*avg(1))
c-----------------------------------------------------------------------
c     compute resistivity and density.
c-----------------------------------------------------------------------
      lambda=24-.5*LOG(ne)+LOG(te)
      taue=3.44e5*te**1.5/(ne*lambda)
      sing%restype%eta=me/(ne*1e6*e**2*taue*1.96)
      sing%restype%rho=ne*mi*1e6
c-----------------------------------------------------------------------
c     compute Alfven and resistive time scales.
c-----------------------------------------------------------------------
      sing%restype%taua=SQRT(sing%restype%rho*sing%restype%m*mu0)
     $     /ABS(twopi*nn*sing%q1*chi1/v1)
      sing%restype%taur=avg(1)/avg(5)*mu0/sing%restype%eta
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(*,10)"q = ",sing%q
         WRITE(*,10)"q1 = ",sing%q1
         WRITE(*,10)"e = ",sing%restype%e
         WRITE(*,10)"f = ",sing%restype%f
         WRITE(*,10)"h = ",sing%restype%h
         WRITE(*,10)"m = ",sing%restype%m
         WRITE(*,10)"g = ",sing%restype%g
         WRITE(*,10)"k = ",sing%restype%k
         WRITE(*,10)"rho = ",sing%restype%rho
         WRITE(*,10)"taue = ",taue
         WRITE(*,10)"taua = ",sing%restype%taua
         WRITE(*,10)"taur = ",sing%restype%taur
         WRITE(*,10)"sfac = ",sing%restype%taur/sing%restype%taua
         WRITE(*,10)"eta = ",sing%restype%eta
         CALL program_stop("Termination by resist_eval.")
      ENDIF
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      CALL spline_dealloc(ff)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_eval
      END MODULE resist_mod
