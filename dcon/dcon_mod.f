c-----------------------------------------------------------------------
c     file dcon_mod.f.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. dcon_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE dcon_mod
      USE spline_mod
      USE global_mod
      USE equil_mod
      USE equil_out_mod
      IMPLICIT NONE

      INTEGER :: out_bal1_unit=21
      INTEGER :: out_bal2_unit=22
      INTEGER :: bin_bal1_unit=23
      INTEGER :: bin_bal2_unit=24
      INTEGER :: fourfit_out_unit=21
      INTEGER :: fourfit_bin_unit=22
      INTEGER :: evals_out_unit=23
      INTEGER :: evals_bin_unit=24
      INTEGER :: crit_out_unit=25
      INTEGER :: crit_bin_unit=26
      INTEGER :: euler_bin_unit=27
      INTEGER :: init_out_unit=28
      INTEGER :: dcon_unit=1
      INTEGER :: unorm_unit=29
      INTEGER :: ca_unit=30
      INTEGER :: err_unit=31

      LOGICAL :: out_bal1=.FALSE.
      LOGICAL :: out_bal2=.FALSE.
      LOGICAL :: bin_bal1=.FALSE.
      LOGICAL :: bin_bal2=.FALSE.
      LOGICAL :: out_fmat=.FALSE.
      LOGICAL :: out_gmat=.FALSE.
      LOGICAL :: out_kmat=.FALSE.
      LOGICAL :: out_metric=.FALSE.
      LOGICAL :: bin_fmat=.FALSE.
      LOGICAL :: bin_gmat=.FALSE.
      LOGICAL :: bin_kmat=.FALSE.
      LOGICAL :: bin_metric=.FALSE.
      LOGICAL :: feval_flag=.FALSE.
      LOGICAL :: fft_flag=.FALSE.
      LOGICAL :: bin_euler=.FALSE.
      LOGICAL :: out_evals=.FALSE.
      LOGICAL :: bin_evals=.FALSE.
      LOGICAL :: out_sol=.FALSE.
      LOGICAL :: bin_sol=.FALSE.

      LOGICAL :: bal_flag=.FALSE.
      LOGICAL :: mat_flag=.FALSE.
      LOGICAL :: ode_flag=.FALSE.
      LOGICAL :: vac_flag=.FALSE.

      LOGICAL :: crit_break=.TRUE.
      LOGICAL :: node_flag=.FALSE.
      LOGICAL :: res_flag=.FALSE.
      LOGICAL :: ahb_flag=.FALSE.

      INTEGER, PARAMETER :: sol_base=50
      INTEGER :: mlow,mhigh,mpert,mband,nn,nstep=HUGE(0),bin_sol_min,
     $     bin_sol_max,euler_stride=1,mthvac=480,ksing=-1,delta_mlow=0,
     $     delta_mhigh=0,delta_mband=0,out_sol_min,out_sol_max
      REAL(r8) :: thmax0=1,ucrit=1e4,tol_r=1e-5,tol_nr=1e-5,
     $     crossover=1e-2,sing_start=0,mthsurf0=1

      TYPE(spline_type) :: locstab

      TYPE :: resist_type
      REAL(r8) :: e,f,h,m,g,k,eta,rho,taua,taur
      END TYPE resist_type

      TYPE ::  sing_type
      INTEGER :: m
      INTEGER, DIMENSION(1) :: r1
      INTEGER, DIMENSION(2) :: r2
      INTEGER, DIMENSION(:), POINTER :: n1,n2
      REAL(r8) :: psifac,rho,q,q1,di
      COMPLEX(r8) :: alpha
      COMPLEX(r8), DIMENSION(:), POINTER :: power
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat,mmat
      TYPE(resist_type) :: restype
      END TYPE sing_type

      INTEGER :: msing
      TYPE(sing_type), DIMENSION(:), POINTER :: sing

      LOGICAL :: sas_flag=.FALSE.,lim_flag
      EQUIVALENCE (sas_flag,lim_flag)
      REAL(r8) :: psilim,qlim,q1lim,dmlim=.5_r8

      END MODULE dcon_mod
