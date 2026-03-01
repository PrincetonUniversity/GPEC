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
      REAL(r8) :: bsq,chi1,di,dpsisq,eta,h,jac,p1,psifac,q,q1,r,z,
     $     rfac,term,theta,twopif,v1,v2,v21,v22,v23,v33,
     $     bt,f1,Dnc,Dnc_prefac,Wc_prefac,eps_loc,Jpara,ftr,M,Jtor,
     $     Hbs_prefac,avg_Jboot_dot_B,mufrac,taua_prefac,taur_prefac,
     $     rmean_loc,amean_loc,deltatop_loc,deltabot_loc,delta_loc,
     $     eps_eff
      REAL(r8), DIMENSION(:), POINTER :: avg
      REAL(r8), DIMENSION(2) :: rmax_loc, rmin_loc, zmax_loc, zmin_loc
      TYPE(spline_type), TARGET :: ff
c-----------------------------------------------------------------------
c     prepare spline types.
c-----------------------------------------------------------------------
      IF(MRE_flag)THEN
         CALL spline_alloc(ff,mtheta,22)
      ELSE
         CALL spline_alloc(ff,mtheta,5)
      ENDIF
      ff%xs=rzphi%ys
c-----------------------------------------------------------------------
c     compute surface quantities.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         twopif=sq%fs(ipsi,1)    ! f = toroidal B field * major radius
         f1=sq%fs1(ipsi,1)/twopi ! df/dpsi, 
         p1=sq%fs1(ipsi,2)
         v1=sq%fs(ipsi,3)     ! d(volume inside flux surface)/d(psi)
         v2=sq%fs1(ipsi,3)
         q=sq%fs(ipsi,4)
         q1=sq%fs1(ipsi,4)
         chi1=twopi*psio      ! d(poloidal flux)/dpsi, chi=poloidal flux
c-----------------------------------------------------------------------
c     initialise extremum surface locations.
c-----------------------------------------------------------------------
         rmax_loc=0.0
         rmin_loc=HUGE(1.0_8)
         zmax_loc=-HUGE(1.0_8)
         zmin_loc=HUGE(1.0_8)
c-----------------------------------------------------------------------
c     evaluate coordinates and jacobian.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            theta=rzphi%ys(itheta)       ! magnetic poloidal angle
            rfac=SQRT(rzphi%f(1))        ! minor radius
            eta=twopi*(theta+rzphi%f(2)) ! machine poloidal angle
            r=ro+rfac*COS(eta)           ! major radius R
            z=zo+rfac*SIN(eta)           ! vertical coordinate Z
            jac=rzphi%f(4)               ! jacobian of mag. coordinates
            bt=twopif/(twopi*r)          ! toroidal B field
c-----------------------------------------------------------------------
c     update extremum surface locations.
c-----------------------------------------------------------------------
            IF(r>rmax_loc(1))THEN
               rmax_loc(1)=r
               rmax_loc(2)=z
            ENDIF
            IF(r<rmin_loc(1))THEN
               rmin_loc(1)=r
               rmin_loc(2)=z
            ENDIF
            IF(z>zmax_loc(2))THEN
               zmax_loc(1)=r
               zmax_loc(2)=z
            ENDIF
            IF(z<zmin_loc(2))THEN
               zmin_loc(1)=r
               zmin_loc(2)=z
            ENDIF
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
            IF(MRE_flag)THEN
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
               ff%fs(itheta,22)=SQRT(dpsisq)
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
         IF(MRE_flag)THEN
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
c     compute surface geometric quantities
c-----------------------------------------------------------------------
            rmean_loc=0.5d0*(rmax_loc(1)+rmin_loc(1))
            amean_loc=0.5d0*(rmax_loc(1)-rmin_loc(1))
            eps_loc=amean_loc/rmean_loc
            deltatop_loc=(rmean_loc-zmax_loc(1))/amean_loc
            deltabot_loc=(rmean_loc-zmin_loc(1))/amean_loc
            delta_loc=0.5d0*(deltatop_loc+deltabot_loc)
            eps_eff=MAX(0.d0, 
     $             0.67d0*(1.d0-1.4d0*delta_loc*ABS(delta_loc))*eps_loc)
c-----------------------------------------------------------------------
c     simple estimates of trapped fraction from Sauter 2002:
c     <https://infoscience.epfl.ch/server/api/core/bitstreams/c42baba0-9
c     909-4f21-978a-0d0c2646c3ad/content>
c-----------------------------------------------------------------------
            ftr=1.d0-SQRT((1.d0-eps_loc)/(1.d0+eps_loc))*(1.d0-eps_eff)/
     $                                        (1.d0+2.d0*SQRT(eps_eff))   
            ftr=MIN(ftr,1.0d0)
c-----------------------------------------------------------------------
c     simple estimates of bootstrap drive-related terms from 
c     Callen, 2010 UW-CPTC 09-6R, and Hegna PoP 1999
c-----------------------------------------------------------------------
            ! mu_e_on_nu_e=(ftr/(1.d0-ftr))*(1.d0+0.533d0/Zeff)
            ! ^ Taking banana limit of eq. B17 (& B14) in Callen 
            mufrac=ftr*(1.d0+0.533d0/Zeff)/
     $                            ((1.d0-ftr)+ftr*(1.d0+0.533d0/Zeff))
            avg_Jboot_dot_B=-mufrac*(twopif/chi1)*p1 !mu0 included in p1
            ! ^Broken, not printing the above out in netcdf
c-----------------------------------------------------------------------
c     evaluate geometric prefactors of MRE stability terms from 
c     Hegna 1999 https://doi.org/10.1063/1.873661
c-----------------------------------------------------------------------
            Dnc_prefac=-q*(p1/(q1*avg(5)))*  !unitless
     $      avg(20)*                         ! \overbar{R^2} ~ [m^2]
     $      avg(1)/(psio**2)                 ! [1/m^2]  
            Dnc=Dnc_prefac*mufrac
c-----------------------------------------------------------------------
c     evaluate geometric prefactor of Wc from 
c     Schlutt and Hegna 2012 https://doi.org/10.1063/1.4747500
c-----------------------------------------------------------------------
            Wc_prefac=(v1*avg(5)/(q*psio))*! overbar{J B^2}           [in psi_tor]
     $      v1*avg(7)*(q*psio)*            ! overbar{J |nabla psi|^2} [in psi_tor]
     $      (q**6/((q1/psio)**2))*         ! {d iota / d psi_tor}^(-2)[in psi_tor]
     $      (1/(q*psio))**4                ! *Conversion from psi_tor^4 to psi_pol_norm^4
            !Divide by mode num. m^2 to finish, units are (psi_pol_norm)^4
            !Note, psi coordinate in original paper is toroidal flux 
            !per radian `psi_tor'. *Using small island approximation.
c-----------------------------------------------------------------------
c     parallel current density from Freidberg Ideal MHD eqs. 6.15, 6.16.
c     note psi in eqs. 6.15, 6.16 is poloidal flux/(2pi), same as psi_in
c-----------------------------------------------------------------------
            Jpara=psio*f1*avg(16)+p1*avg(18)*twopif/(twopi*psio) + 
     $      (twopif/twopi)**2*f1*avg(17)/psio !(mu0 included in p1))
            Jtor=p1*avg(12)/psio + (twopif/twopi)*f1*avg(13)/psio 
c-----------------------------------------------------------------------
c     compute Hbs_prefac from Shi et al. 2024, using identity from 
c     Glasser et al. 1975. Shi doi -> https://doi.org/10.1063/5.0183474 
c-----------------------------------------------------------------------
            Hbs_prefac=(avg(1)/avg(5))*
     $  (-v1/(twopi**2*psio**2*q1))! mult by avg_Jboot_dot_B to get Hbs
c-----------------------------------------------------------------------
c     save terms.
c-----------------------------------------------------------------------
            mreterms%fs(ipsi,1)=Hbs_prefac
            mreterms%fs(ipsi,2)=taua_prefac
            mreterms%fs(ipsi,3)=taur_prefac
            mreterms%fs(ipsi,4)=ftr
            mreterms%fs(ipsi,5)=mufrac
            mreterms%fs(ipsi,6)=avg(22) ! <|nabla psi|>
            mreterms%fs(ipsi,7)=Dnc
            mreterms%fs(ipsi,8)=Wc_prefac
            mreterms%fs(ipsi,9)=Jpara
            mreterms%fs(ipsi,10)=avg(8) !Avg B field
            mreterms%fs(ipsi,11)=avg(9) !Avg toroidal B field
            mreterms%fs(ipsi,12)=avg(10) !Avg poloidal B field
            mreterms%fs(ipsi,13)=avg(11) !Avg minor radius
            mreterms%fs(ipsi,14)=avg(12) !Avg major radius
            mreterms%fs(ipsi,15)=avg(13) !Avg 1/major radius
            mreterms%fs(ipsi,16)=avg(20) !overbar{R^2}     (Hegna 1999)
            mreterms%fs(ipsi,17)=avg(21) !avg{R^2} ~ [m^2] (Hegna 1999)
            ! will only print out the following if geom_flag is true:
            mreterms%fs(ipsi,18)=avg(1) ! <bsq/dpsisq>
            mreterms%fs(ipsi,19)=avg(2) ! <1/dpsisq>
            mreterms%fs(ipsi,20)=avg(3) ! <1/bsq>
            mreterms%fs(ipsi,21)=avg(4) ! <1/(bsq*dpsisq)>
            mreterms%fs(ipsi,22)=avg(5) ! <bsq>
            mreterms%fs(ipsi,23)=avg(6) ! <dpsisq/bsq>
            mreterms%fs(ipsi,24)=avg(7) ! <dpsisq>
            mreterms%fs(ipsi,25)=avg(14)! <dpsisq/(R^2)>
            mreterms%fs(ipsi,26)=avg(15)! <1/(R^2)> 
            mreterms%fs(ipsi,27)=avg(16)! <dpsisq/(R^2 B))>
            mreterms%fs(ipsi,28)=avg(17)! <1/(R^2 B)>
            mreterms%fs(ipsi,29)=avg(18)! <1/B>
            mreterms%fs(ipsi,30)=avg(19)! <sqrt(dpsisq)/R>
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
