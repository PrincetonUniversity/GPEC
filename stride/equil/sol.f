c-----------------------------------------------------------------------
c     file sol.f.
c     computes Soloviev's analytical equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. sol_mod.
c     1. sol_run.
c     2. sol_galkin.
c-----------------------------------------------------------------------
c     subprogram 0. sol_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE sol_mod
      USE direct_mod
      USE inverse_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. sol_run.
c     computes Soloviev's analytical equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sol_run
      
      INTEGER :: ir,iz,ia
      REAL(r8) :: psifac,efac,f0,pfac
      REAL(r8), DIMENSION(:,:), POINTER :: rg,zg
      REAL(r8), DIMENSION(:), POINTER :: r,z

      INTEGER :: mr=65      ! number of radial grid zones
      INTEGER :: mz=65      ! number of axial grid zones
      INTEGER :: ma=64      ! number of flux grid zones
      REAL(r8) :: e=1           ! elongation
      REAL(r8) :: a=1           ! minor radius
      REAL(r8) :: r0=3          ! major radius
      REAL(r8) :: q0=1.26       ! safety factor at the o-point

      NAMELIST/sol_input/mr,mz,ma,e,a,r0,q0
c-----------------------------------------------------------------------
c     read input data.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,eq_filename,"OLD")
      READ(in_unit,NML=sol_input)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(r(0:mr),z(0:mz),rg(0:mr,0:mz),zg(0:mr,0:mz))
      CALL bicube_alloc(psi_in,mr,mz,1)
      CALL spline_alloc(sq_in,ma,4)
c-----------------------------------------------------------------------
c     compute scalar data.
c-----------------------------------------------------------------------
      ro=0
      zo=0
      f0=r0
      psio=e*f0*a*a/(2*q0*r0)
      psifac=psio/(a*r0)**2
      efac=1/(e*e)
      pfac=2*psio**2*(e*e+1)/(a*r0*e)**2
      rmin=r0-1.5*a
      rmax=r0+1.5*a
      zmax=1.5*e*a
      zmin=-zmax
c-----------------------------------------------------------------------
c     compute 1D data.
c-----------------------------------------------------------------------
      sq_in%xs=(/(ia,ia=1,ma+1)/)/REAL(ma+1,r8)
      sq_in%xs=sq_in%xs**2
      sq_in%fs(:,1)=f0
      sq_in%fs(:,2)=pfac*(1-sq_in%xs)
      sq_in%fs(:,3)=0
c-----------------------------------------------------------------------
c     compute 2D data.
c-----------------------------------------------------------------------
      r=rmin+(/(ir,ir=0,mr)/)*(rmax-rmin)/mr
      z=zmin+(/(iz,iz=0,mz)/)*(zmax-zmin)/mz
      DO iz=0,mz
         DO ir=0,mr
            rg(ir,iz)=r(ir)
            zg(ir,iz)=z(iz)
            psi_in%fs(ir,iz,1)=psio
     $           -psifac*(efac*(r(ir)*z(iz))**2+(r(ir)**2-r0**2)**2/4)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL direct_run
      DEALLOCATE(r,z,rg,zg)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sol_run
c-----------------------------------------------------------------------
c     subprogram 2. sol_galkin.
c     computes Sergei Galkin's inverse Soloviev equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sol_galkin

      INTEGER :: ma=128,mt=128,ia,it
      REAL(r8), DIMENSION(:), POINTER :: f,p,q
      REAL(r8), DIMENSION(:,:), POINTER :: rsq,r,z

      REAL(r8) :: del0,del1,sigma,psimax,q0,c0,alpha,f0
      REAL(r8), DIMENSION(:), POINTER :: omega,theta,costh,sinth,a,
     $     lambda

      NAMELIST/galkin_input/ma,mt,del0,del1,sigma,psimax,q0
c-----------------------------------------------------------------------
c     read control data.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,TRIM(eq_filename),"OLD")
      READ(in_unit,NML=galkin_input)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(f(0:ma),p(0:ma),q(0:ma),a(0:ma),omega(0:ma),lambda(0:ma),
     $     theta(0:mt),costh(0:mt),sinth(0:mt),rsq(0:mt,0:ma),
     $     r(0:mt,0:ma),z(0:mt,0:ma))
c-----------------------------------------------------------------------
c     scalars.
c-----------------------------------------------------------------------
      alpha=1/SQRT(1-2*sigma**2/(del1**2+del0**2))
      c0=4*psimax/(del1**2-del0**2)**2
      f0=SQRT(8*(del1**2+del0**2-2*sigma**2))
     $     *(del1**2+del0**2)*c0*q0*alpha
c-----------------------------------------------------------------------
c     radial arrays.
c-----------------------------------------------------------------------
      a=((/(ia,ia=0,ma)/)/REAL(ma,4))**2
      omega=(del1**2-del0**2)/(del1**2+del0**2)*SQRT(a)
      lambda=.5*(sigma/q0)**2*(del1**2-del0**2)**2
     $     /((del1**2+sigma**2)**2*(del1**2+del0**2-2*sigma**2))
      p=2*(1+alpha**2)*(del1**2-del0**2)**2*c0**2*(1-a)
      f=f0*SQRT(1+lambda*a)
c-----------------------------------------------------------------------
c     angular arrays.
c-----------------------------------------------------------------------
      theta=twopi*(/(it,it=0,mt)/)/mt
      costh=COS(theta)
      sinth=SIN(theta)
c-----------------------------------------------------------------------
c     position arrays.
c-----------------------------------------------------------------------
      DO it=0,mt
         rsq(:,it)=(del1**2+del0**2)/2*(1+omega*costh(it))
         z(:,it)=(del1**2+del0**2)/(4*alpha)*omega*sinth(it)
     $        /SQRT(rsq(:,it)-sigma**2)
      ENDDO
      r=SQRT(rsq)
c-----------------------------------------------------------------------
c     fill 1D arrays.
c-----------------------------------------------------------------------
      psio=psimax
      CALL spline_alloc(sq_in,ma,4)
      sq_in%xs=a
      sq_in%fs(:,1)=f
      sq_in%fs(:,2)=p
      sq_in%fs(:,3)=0
c-----------------------------------------------------------------------
c     fill 2D arrays.
c-----------------------------------------------------------------------
      ro=SQRT((del1**2+del0**2)/2)
      zo=0
      CALL bicube_alloc(rz_in,ma,mt,2)
      rz_in%fs(:,:,1)=r
      rz_in%fs(:,:,2)=z
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(f,p,q,a,omega,lambda,theta,costh,sinth,rsq,r,z)
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sol_galkin
      END MODULE sol_mod
