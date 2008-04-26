c-----------------------------------------------------------------------
c     file global.f.
c     global declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. global_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE global_mod
      USE bicube_mod
      IMPLICIT NONE

      CHARACTER(16) :: eq_type="fluxgrid"
      CHARACTER(128) :: eq_filename="fluxgrid.dat"

      LOGICAL :: direct_flag=.FALSE.
      LOGICAL :: out_eq_1d=.FALSE.,bin_eq_1d=.FALSE.
      LOGICAL :: out_eq_2d=.FALSE.,bin_eq_2d=.FALSE.
      LOGICAL :: out_2d=.FALSE.,bin_2d=.FALSE.,interp=.FALSE.
      LOGICAL :: input_only=.FALSE.
      LOGICAL :: out_fl=.FALSE.
      LOGICAL :: bin_fl=.FALSE.
      LOGICAL :: gse_flag=.FALSE.
      LOGICAL :: power_flag=.TRUE.
      LOGICAL :: dump_flag=.FALSE.

      REAL(r8) :: psilow=1e-4
      REAL(r8) :: psihigh=1-1e-6
      REAL(r8) :: newq0=0

      INTEGER :: mex
      REAL(r8), DIMENSION(:), POINTER :: qex,psiex

      CHARACTER(16) :: jac_type="hamada"
      INTEGER :: power_bp=0,power_b=0,power_r=0,jac_method=1

      CHARACTER(16) :: grid_type="ldp"
      INTEGER :: mpsi=128,mtheta=128
      REAL(r8) :: ro,zo,psio,q0,qa,qmin,qmax,amean,rmean,aratio,kappa,
     $     delta1,delta2,bt0,crnt,betat,betan,betap1,betap2,betap3,
     $     li1,li2,li3,volume,p0,ppeakfac,q95
      REAL(r8), DIMENSION(2) :: rext,rsep,zsep
      TYPE(spline_type) :: sq,sq_in
      TYPE(bicube_type) :: rzphi

      END MODULE global_mod
