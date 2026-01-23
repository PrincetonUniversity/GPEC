      MODULE sglobal_mod
      USE local_mod, ONLY: r8

      IMPLICIT NONE
      INTEGER :: mm,nn
      INTEGER :: in_unit,out_unit,out2_unit,out3_unit,
     $     bin_unit,bin_2d_unit,input_unit,n_trace
c      INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(13,307)

      REAL(r8) :: mr,nr
      REAL(r8) :: Q_e,Q_i,pr,pe,c_beta,ds,tau,d_i,
     $            d_beta,D_norm,P_perp,P_tor,gamma_fac
      REAL(r8) :: eta,visc,rho_s,lu,omega_e,omega_i,iota_e,
     $            delta_n,layfac,Qconv,lnLamb,deltaprim,dc_tmp,
     $            d_crit,tau_r,tauk,g_r,g_i,delta_eff
      COMPLEX(r8) :: Q,g_tmp,delta_det
      CHARACTER(20) :: dc_type
     
      REAL(r8), PARAMETER :: pi=3.1415926535897932385, mu0=4e-7*pi,
     $     m_e=9.1094e-31,m_p=1.6726e-27,chag=1.6021917e-19,
     $     kval=1.3807e-23,eps0 = 8.8542e-12

      INTEGER, PARAMETER :: MAX_PTS = 500000   ! Max unique points allowed
      INTEGER, PARAMETER :: HASH_SZ = 500009   ! Prime number for hash table
      REAL(r8), PARAMETER :: HASH_SCALE = 1.0d5                  ! Scaling factor for integer hashing
      INTEGER, ALLOCATABLE :: hash_head(:)     ! Hash bucket heads
      INTEGER, ALLOCATABLE :: hash_next(:)     ! Linked list next pointers
      INTEGER :: n_pts
      COMPLEX(r8), ALLOCATABLE :: Q_store(:)    ! Stores Q coordinates
      COMPLEX(r8), ALLOCATABLE :: D_store(:)    ! Stores Result Delta

      TYPE result_type
          REAL(r8), ALLOCATABLE :: inQs(:), iinQs(:),
     $     Re_deltas(:), Im_deltas(:)
          INTEGER :: count
      END TYPE result_type

      TYPE slayer_inputs_type
          INTEGER, ALLOCATABLE :: qval_arr(:)
          REAL(r8), ALLOCATABLE :: chi_p_arr(:),chi_t_arr(:),
     $      kappa_arr(:),psi_n_arr(:),
     $      lu_arr(:),Qconv_arr(:),Q_e_arr(:),Q_i_arr(:),c_beta_arr(:),
     $      d_beta_arr(:),D_norm_arr(:),tau_arr(:),P_perp_arr(:),
     $      P_tor_arr(:),omegas_arr(:),omegas_e_arr(:),omegas_i_arr(:),
     $      gammafac_arr(:),Re_dp_arr(:),Im_dp_arr(:),d_crit_arr(:)
          COMPLEX(r8), ALLOCATABLE :: dp_matrix(:,:)
      END TYPE slayer_inputs_type

      TYPE slayer_outputs_type
          COMPLEX(r8), ALLOCATABLE :: dels_db_arr(:),gamma_sol_arr(:),
     $      gamma_est_arr(:)
      END TYPE slayer_outputs_type

      TYPE deltas_outputs_type
          REAL(r8), ALLOCATABLE :: inQs(:)
          REAL(r8), ALLOCATABLE :: iinQs(:)
          REAL(r8), ALLOCATABLE :: real_deltas(:)
          REAL(r8), ALLOCATABLE :: imag_deltas(:)
      END TYPE deltas_outputs_type

      ! lnLamb will be updated.

      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      CHARACTER(2) :: sn,sm

      END MODULE sglobal_mod
