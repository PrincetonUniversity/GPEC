c-----------------------------------------------------------------------
c     orbit_mod.f.
c     module declarations for orbit parameters.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE orbit_mod
      USE equil_mod
      USE equil_out_mod
      IMPLICIT NONE

      INTEGER :: field_out_unit=25
      INTEGER :: field_bin_unit=26
      INTEGER :: ode_out_unit=27
      INTEGER :: ode_bin_unit=28

      LOGICAL :: orbit_out=.FALSE.,orbit_bin=.FALSE.,
     $     break_flag=.TRUE.,out_flag=.FALSE.,
     $     restart_flag=.FALSE.,ode_bin=.TRUE.,
     $     gc_flag=.TRUE.
      CHARACTER(9) :: particle="ion",orbit_type="inner"
      INTEGER :: stride=1,ai=2,zi=1,nstep=HUGE(0)
      REAL(r8) :: alpha0,energy0,energy,error,f0,taumax,m,omega0,
     $     phase0,phi0,psi0,q,r0,report=100,rtol,tau_report,theta0,
     $     tmax,z0,zeta0,errmax=1e-2,tol=1e-8,mu_0,mu

      INTEGER :: jac,mf,itol,itask,istate,iopt
      INTEGER, PARAMETER :: liw=20
      INTEGER, DIMENSION(liw) :: iwork
      REAL(r8), DIMENSION(:), POINTER :: rwork

      END MODULE orbit_mod
