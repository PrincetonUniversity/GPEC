      MODULE sglobal_mod
      USE local_mod, ONLY: r8

      IMPLICIT NONE
      INTEGER :: mm,nn
      INTEGER :: in_unit,out_unit,out2_unit,out3_unit,
     $     bin_unit,bin_2d_unit,input_unit,n_trace
c      INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(13,307)

      REAL(r8) :: mr,nr
      REAL(r8) :: Q_e,Q_i,pr,pe,c_beta,ds,tau,d_i,
     $            d_beta,D_norm,P_perp,gamma_fac
      REAL(r8) :: eta,visc,rho_s,lu,omega_e,omega_i,
     $            delta_n,layfac,Qconv,lnLamb,deltaprim,dc_tmp,
     $            d_crit,tau_r,tauk,g_r,g_i,delta_eff
      REAL(r8), DIMENSION(:), ALLOCATABLE :: re_trace,im_trace
      COMPLEX(r8) :: Q,g_tmp
      CHARACTER(20) :: dc_type
     
      REAL(r8), PARAMETER :: pi=3.1415926535897932385, mu0=4e-7*pi,
     $     m_e=9.1094e-31,m_p=1.6726e-27,chag=1.6021917e-19,
     $     kval=1.3807e-23,eps0 = 8.8542e-12

      TYPE result_type
          REAL(r8), ALLOCATABLE :: inQs(:), iinQs(:),
     $     Re_deltas(:), Im_deltas(:)
          INTEGER :: count
      END TYPE result_type

      ! lnLamb will be updated.

      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      CHARACTER(2) :: sn
      
      END MODULE sglobal_mod
