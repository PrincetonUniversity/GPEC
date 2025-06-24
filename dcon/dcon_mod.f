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
      USE pentrc_interface, ONLY: version

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
      LOGICAL :: out_fund=.FALSE.
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
      LOGICAL :: bin_vac=.FALSE.
      LOGICAL :: out_evals=.FALSE.
      LOGICAL :: bin_evals=.FALSE.
      LOGICAL :: out_sol=.FALSE.
      LOGICAL :: bin_sol=.FALSE.
      LOGICAL :: netcdf_out=.TRUE.

      LOGICAL :: bal_flag=.FALSE.
      LOGICAL :: mat_flag=.FALSE.
      LOGICAL :: ode_flag=.FALSE.
      LOGICAL :: vac_flag=.FALSE.
      LOGICAL :: mer_flag=.FALSE.

      LOGICAL :: crit_break=.TRUE.
      LOGICAL :: node_flag=.FALSE.
      LOGICAL :: res_flag=.FALSE.
      LOGICAL :: ahb_flag=.FALSE.
      LOGICAL :: out_ahg2msc=.TRUE.
      LOGICAL :: vac_memory=.FALSE.

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
      INTEGER :: order
      INTEGER, DIMENSION(1) :: r1
      INTEGER, DIMENSION(2) :: r2
      INTEGER, DIMENSION(:), POINTER :: n1,n2
      REAL(r8) :: psifac,rho,q,q1,di
      COMPLEX(r8) :: alpha
      COMPLEX(r8), DIMENSION(:), POINTER :: power
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat,mmat
      TYPE(resist_type) :: restype
      END TYPE sing_type

      INTEGER :: msing,kmsing
      TYPE(sing_type), DIMENSION(:), POINTER :: sing
      TYPE(sing_type), DIMENSION(:), POINTER :: kinsing

      LOGICAL :: sas_flag=.FALSE.,lim_flag
      EQUIVALENCE (sas_flag,lim_flag)
      REAL(r8) :: psilim,qlim,q1lim,dmlim=.5_r8,qhigh=1e3,qlow=0

      INTEGER :: kingridtype=0
      REAL(r8) :: kinfac1=1.0,kinfac2=1.0,ktc=0.1,ktw=50
      LOGICAL :: kin_flag = .FALSE.
      LOGICAL :: con_flag = .FALSE.
      LOGICAL :: ktanh_flag = .FALSE.
      LOGICAL :: passing_flag = .FALSE.
      LOGICAL :: trapped_flag = .TRUE.
      LOGICAL :: electron_flag = .FALSE.
      LOGICAL :: ion_flag = .TRUE.
      LOGICAL :: fkg_kmats_flag = .FALSE.
      LOGICAL :: keq_out = .FALSE.
      LOGICAL :: theta_out = .FALSE.
      LOGICAL :: xlmda_out = .FALSE.

      ! for recording dw(q) near the boundary
      REAL(r8) :: psiedge = 1.0
      REAL(r8), DIMENSION(:), ALLOCATABLE :: q_edge, psi_edge
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: dw_edge
      INTEGER :: nperq_edge=20, size_edge=0, pre_edge=1, i_edge=1

      END MODULE dcon_mod
