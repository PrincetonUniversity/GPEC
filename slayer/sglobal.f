c=======================================================================
c     MODULE sglobal_mod
c
c     Global shared state for the SLAYER layer-physics package.
c     Contains:
c       - physical and mathematical constants
c       - module-level scalar variables written by params() and
c         consumed by the dispersion solvers (delta, slayer, gslayer)
c       - AMR scanner storage (hash-based v1, cell-based v2)
c       - derived types: slayer_inputs_type, slayer_outputs_type,
c         deltas_outputs_type, result_type, amr_cell_type
c
c     TODO: `pr` and `pe` (magnetic Prandtl number and its
c       electron analogue) are declared but never explicitly set
c       within this module or in params().  They must be initialised
c       before params() is called (tau_v = tau_r / pr).  Consider
c       adding a dedicated init routine or making them INTENT(IN)
c       arguments of params().
c
c=======================================================================
      MODULE sglobal_mod

      USE local_mod, ONLY: r8     ! double-precision kind parameter

      IMPLICIT NONE

c-----------------------------------------------------------------------
c     I/O unit numbers and trace control.
c-----------------------------------------------------------------------
      INTEGER :: in_unit           ! standard input  unit
      INTEGER :: out_unit          ! primary output  unit
      INTEGER :: out2_unit         ! secondary output unit
      INTEGER :: out3_unit         ! tertiary output  unit
      INTEGER :: bin_unit          ! binary output    unit
      INTEGER :: bin_2d_unit       ! 2-D binary output unit
      INTEGER :: input_unit        ! namelist / input  unit
      INTEGER :: n_trace           ! trace / debug verbosity level

c-----------------------------------------------------------------------
c     Mode numbers (integer and real representations).
c-----------------------------------------------------------------------
      INTEGER  :: mm               ! poloidal mode number (integer)
      INTEGER  :: nn               ! toroidal mode number (integer)
      REAL(r8) :: mr               ! poloidal mode number (real copy)
      REAL(r8) :: nr               ! toroidal mode number (real copy)
      CHARACTER(2) :: sn_str       ! toroidal n as string
      CHARACTER(2) :: sm_str       ! poloidal m as string

c-----------------------------------------------------------------------
c     Layer-physics scalars (set by params(), read by solvers).
c     These are the normalised parameters that enter the SLAYER
c     dispersion relation.
c-----------------------------------------------------------------------
c --- temperature / collisionality
      REAL(r8) :: tau              ! T_i / T_e
      REAL(r8) :: eta              ! Spitzer resistivity [Ohm*m]
      REAL(r8) :: visc             ! anomalous viscosity [m^2/s]
      REAL(r8) :: lnLamb           ! Coulomb logarithm (updated at runtime)
c --- length scales
      REAL(r8) :: rho_s            ! ion Larmor radius at T_e [m]
      REAL(r8) :: d_i              ! ion skin depth [m]
      REAL(r8) :: d_beta           ! beta-weighted ion scale d_beta = c_beta * d_i
c --- timescales
      REAL(r8) :: tau_r            ! resistive diffusion time [s]
      REAL(r8) :: tauk             ! Q-conversion factor (= Qconv)
c --- Lundquist and Prandtl numbers
      REAL(r8) :: lu               ! Lundquist number S = tau_r / tau_h
      REAL(r8) :: pr               ! magnetic Prandtl number (TODO: see header)
      REAL(r8) :: pe               ! electron Prandtl number (TODO: see header)
      REAL(r8) :: P_perp           ! perpendicular magnetic Prandtl number
      REAL(r8) :: P_tor            ! toroidal magnetic Prandtl number
c --- normalised layer parameters
      REAL(r8) :: ds               ! normalised ion Larmor radius
      REAL(r8) :: c_beta           ! compressional beta parameter
      REAL(r8) :: D_norm           ! normalised beta-weighted ion scale
      REAL(r8) :: delta_n          ! Delta normalisation factor
      REAL(r8) :: Qconv            ! frequency normalisation (Cole)
c --- diamagnetic and rotation frequencies
      REAL(r8) :: omega_e          ! electron diamagnetic frequency [rad/s]
      REAL(r8) :: omega_i          ! ion diamagnetic frequency [rad/s]
      COMPLEX(r8) :: Q             ! normalised complex rotation frequency
      REAL(r8) :: Q_e              ! normalised electron diamagnetic Q
      REAL(r8) :: Q_i              ! normalised ion diamagnetic Q
c --- stability / Delta_crit
      REAL(r8) :: dc_tmp           ! computed Delta_crit
      REAL(r8) :: delta_eff        ! effective Deltaprime shift
      CHARACTER(20) :: dc_type     ! dc formula selector ('lar','rfitzp','toroidal')
c --- solver workspace / results
      COMPLEX(r8) :: g_tmp         ! temporary complex growth rate
      REAL(r8) :: gamma_fac        ! growth-rate conversion factor
c --- miscellaneous
      REAL(r8) :: iota_e           ! electron iota
      REAL(r8) :: layfac           ! layer singularity guard factor

c-----------------------------------------------------------------------
c     Physical and mathematical constants.
c-----------------------------------------------------------------------
      REAL(r8), PARAMETER :: pi   = 3.1415926535897932385d0
      REAL(r8), PARAMETER :: mu0  = 4.0d-7 * pi       ! vacuum permeability [H/m]
      REAL(r8), PARAMETER :: m_e  = 9.1094d-31         ! electron mass [kg]
      REAL(r8), PARAMETER :: m_p  = 1.6726d-27         ! proton mass   [kg]
      REAL(r8), PARAMETER :: chag = 1.6021917d-19      ! elementary charge [C]
      REAL(r8), PARAMETER :: kval = 1.3807d-23         ! Boltzmann constant [J/K]
      REAL(r8), PARAMETER :: eps0 = 8.8542d-12         ! vacuum permittivity [F/m]
      COMPLEX(r8), PARAMETER :: ifac = (0.0d0, 1.0d0) ! imaginary unit

c-----------------------------------------------------------------------
c     AMR scanner storage -- hash-based deduplication (v1).
c-----------------------------------------------------------------------
      INTEGER, PARAMETER  :: MAX_PTS    = 500000  ! max unique eval points
      INTEGER, PARAMETER  :: HASH_SZ    = 500009  ! hash table size (prime)
      REAL(r8), PARAMETER :: HASH_SCALE = 1.0d5   ! Re/Im quantisation scale
      INTEGER, ALLOCATABLE :: hash_head(:)         ! bucket heads  (HASH_SZ)
      INTEGER, ALLOCATABLE :: hash_next(:)         ! chain pointers (MAX_PTS)

c-----------------------------------------------------------------------
c     AMR scanner storage -- cell-based refinement (v2).
c-----------------------------------------------------------------------
      INTEGER, PARAMETER :: MAX_CELLS = 500000    ! max AMR cells

      TYPE :: amr_cell_type
          COMPLEX(r8) :: Q(4)          ! corner Q-values (BL, BR, TL, TR)
          COMPLEX(r8) :: D(4)          ! corner dispersion values
          LOGICAL      :: needs_refine ! flagged for subdivision
      END TYPE amr_cell_type

      TYPE(amr_cell_type), ALLOCATABLE :: amr_cells(:)
      INTEGER :: n_amr_cells           ! current number of active cells

c-----------------------------------------------------------------------
c     AMR output arrays (shared by v1 and v2).
c-----------------------------------------------------------------------
      COMPLEX(r8), ALLOCATABLE :: Q_store(:)  ! unique Q-points
      COMPLEX(r8), ALLOCATABLE :: D_store(:)  ! corresponding D-values
      INTEGER :: n_pts                         ! number of stored points

c-----------------------------------------------------------------------
c     Derived types: solver I/O and scan results.
c-----------------------------------------------------------------------

c     result_type -- torque-scan output bucket
      TYPE result_type
          REAL(r8), ALLOCATABLE :: inQs(:)       ! Re(Q) scan values
          REAL(r8), ALLOCATABLE :: iinQs(:)      ! Im(Q) scan values
          REAL(r8), ALLOCATABLE :: Re_deltas(:)  ! Re(Delta) results
          REAL(r8), ALLOCATABLE :: Im_deltas(:)  ! Im(Delta) results
          INTEGER :: count                        ! number of entries
      END TYPE result_type

c     slayer_inputs_type -- per-surface input arrays for SLAYER
      TYPE slayer_inputs_type
          INTEGER, ALLOCATABLE  :: qval_arr(:)        ! safety-factor integers
          REAL(r8), ALLOCATABLE :: chi_p_arr(:)       ! chi_perp  [m^2/s]
          REAL(r8), ALLOCATABLE :: chi_t_arr(:)       ! chi_tor   [m^2/s]
          REAL(r8), ALLOCATABLE :: kappa_arr(:)       ! kappa (thermal cond.)
          REAL(r8), ALLOCATABLE :: psi_n_arr(:)       ! normalised psi
          REAL(r8), ALLOCATABLE :: lu_arr(:)          ! Lundquist number
          REAL(r8), ALLOCATABLE :: Qconv_arr(:)       ! Q-conversion factor
          REAL(r8), ALLOCATABLE :: Q_e_arr(:)         ! normalised Q_e
          REAL(r8), ALLOCATABLE :: Q_i_arr(:)         ! normalised Q_i
          REAL(r8), ALLOCATABLE :: c_beta_arr(:)      ! compressional beta
          REAL(r8), ALLOCATABLE :: d_beta_arr(:)      ! beta-related width
          REAL(r8), ALLOCATABLE :: D_norm_arr(:)      ! normalised D
          REAL(r8), ALLOCATABLE :: tau_arr(:)         ! T_i / T_e
          REAL(r8), ALLOCATABLE :: P_perp_arr(:)      ! perp Prandtl number
          REAL(r8), ALLOCATABLE :: P_tor_arr(:)       ! toroidal Prandtl
          REAL(r8), ALLOCATABLE :: omegas_arr(:)      ! rotation [rad/s]
          REAL(r8), ALLOCATABLE :: omegas_e_arr(:)    ! omega_e  [rad/s]
          REAL(r8), ALLOCATABLE :: omegas_i_arr(:)    ! omega_i  [rad/s]
          REAL(r8), ALLOCATABLE :: gammafac_arr(:)    ! gamma conversion
          REAL(r8), ALLOCATABLE :: Re_dp_arr(:)       ! Re(Deltaprime)
          REAL(r8), ALLOCATABLE :: Im_dp_arr(:)       ! Im(Deltaprime)
          REAL(r8), ALLOCATABLE :: d_crit_arr(:)      ! Delta_crit
          COMPLEX(r8), ALLOCATABLE :: dp_matrix(:,:)  ! full Deltaprime matrix
      END TYPE slayer_inputs_type

c     slayer_outputs_type -- per-surface solver results
      TYPE slayer_outputs_type
          COMPLEX(r8), ALLOCATABLE :: dels_db_arr(:)   ! Delta from d_beta
          COMPLEX(r8), ALLOCATABLE :: gamma_sol_arr(:) ! solved growth rate
          COMPLEX(r8), ALLOCATABLE :: gamma_est_arr(:) ! estimated growth rate
          REAL(r8), ALLOCATABLE :: br_th_arr(:)        ! Br threshold
      END TYPE slayer_outputs_type

c     deltas_outputs_type -- scan output (Q vs Delta)
      TYPE deltas_outputs_type
          REAL(r8), ALLOCATABLE :: inQs(:)        ! Re(Q) values
          REAL(r8), ALLOCATABLE :: iinQs(:)       ! Im(Q) values
          REAL(r8), ALLOCATABLE :: real_deltas(:)  ! Re(Delta)
          REAL(r8), ALLOCATABLE :: imag_deltas(:)  ! Im(Delta)
      END TYPE deltas_outputs_type

      END MODULE sglobal_mod
