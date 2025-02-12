c-----------------------------------------------------------------------
c     file slayer_netcdf.f
c     writes slayer information to a netcdf file
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. slayer_netcdf_mod
c     1. check
c     2. stride_netcdf_out
c-----------------------------------------------------------------------
c     subprogram 0. slayer_netcdf_mod
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE slayer_netcdf_mod
      USE sglobal_mod
      USE netcdf
      IMPLICIT NONE
      CONTAINS
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c      subprogram 2. stride_netcdf_out.
c      Replicate stride.out information in netcdf format.
c -----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 1. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE check(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
        PRINT *, TRIM(nf90_strerror(stat))
         !STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c -----------------------------------------------------------------------
c      declarations.
c -----------------------------------------------------------------------
      SUBROUTINE slayer_netcdf_out(msing,lar_gamma_eq_flag,
     $  lar_gamma_flag,stabscan_eq_flag,stabscan_flag,br_th_flag,
     $  qval_arr,omegas_arr,inQ_arr,inQ_e_arr,inQ_i_arr,psi_n_rational,
     $  inpr_arr,br_th,Re_deltaprime_arr,Im_deltaprime_arr,dels_db,
     $  d_beta,D_beta_norm,lar_gamma,inQs,iinQs,results)

        ! ds = D_beta_norm for lar growth rate routines

c        OPTIONAL
c        br_th,Re_deltaprime_arr,Im_deltaprime_arr,dels_db,d_b,ds,lar_gamma,
c        inQs,iinQs,results,

      INTEGER, INTENT(IN) :: msing
      LOGICAL, INTENT(IN) :: lar_gamma_eq_flag,lar_gamma_flag,
     $   stabscan_eq_flag,stabscan_flag,br_th_flag

      INTEGER, ALLOCATABLE, DIMENSION(:) :: qval_arr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: omegas_arr,
     $ inQ_arr,inQ_e_arr,inQ_i_arr,psi_n_rational,inpr_arr

      REAL(r8), DIMENSION(:), ALLOCATABLE :: Re_deltaprime_arr,
     $         Im_deltaprime_arr,inQs,iinQs

      REAL(r8) :: br_th, d_beta, D_beta_norm
      COMPLEX(r8) :: dels_db, lar_gamma

      TYPE(result_type), INTENT(IN) :: results(8)

      INTEGER :: i,ncid,r_id,ReQ_dim,ImQ_dim,qsing_dim,
     $    i_dim, m_dim, mo_dim, p_dim, i_id, m_id, mo_id, p_id,
     $    ReQ_id,ImQ_id,gamma_id,omegas_id,Q_id,Q_e_id,Q_i_id,
     $    r_dim,pr_id, qr_id,shear_id,slice_id,inQs_id,S_id,
     $    gamma_err_id,gamma_loc_id,roots_dim,Re_dp_id,Im_dp_id,
     $    rdpp_id,idpp_id,inpr_id,br_th_id,dels_db_id,d_b_id,
     $    inds_id,lar_gamma_id,qsing_id

      INTEGER :: run, run_dimid, point_dimid, varids(4),
     $     max_points

      CHARACTER(64) :: ncfile
      LOGICAL, PARAMETER :: debug_flag = .FALSE.
      CHARACTER(len=*), PARAMETER :: version ='v1.0.0-99-gc873bd6'
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
c      i=0; ncid=0;r_id=0;ReQ_dim=0;ImQ_dim=0;qsing_dim=0;qsing_id=0;
c     $    i_dim=0;m_dim=0;mo_dim=0;p_dim=0;i_id=0;m_id=0;mo_id=0;p_id=0;
c     $    ReQ_id=0;ImQ_id=0;gamma_id=0;omegas_id=0;Q_id=0;Q_e_id=0;
c     $    Q_i_id=0;r_dim=0;pr_id=0;qr_id=0;shear_id=0;slice_id=0;
c     $    inQs_id=0;S_id=0;gamma_err_id=0;gamma_loc_id=0;roots_dim=0;
c     $    Re_dp_id=0;Im_dp_id=0;rdpp_id=0;idpp_id=0;inpr_id=0;
c     $    br_th_id=0;dels_db_id=0;d_b_id=0;inds_id=0;lar_gamma_id=0
      IF(debug_flag) PRINT *,"Called slayer_netcdf_out"
      IF (nn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)')nn
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')nn
      ENDIF
      ncfile = "slayer_output_n"//TRIM(sn)//".nc"
      IF(debug_flag) PRINT *, ncfile
c -----------------------------------------------------------------------
c      open files
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Creating netcdf files"
      CALL check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

      max_points = maxval([(results(run)%count, run=1,msing)])
c -----------------------------------------------------------------------
c      define global file attributes
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining netcdf globals"
      CALL check( nf90_put_att(ncid,nf90_global,"title",
     $     "SLAYER outputs"))
      !CALL check( nf90_put_att(ncid,nf90_global,"shot", INT(shotnum)) )
      !CALL check( nf90_put_att(ncid,nf90_global,"time",INT(shottime)) )
      !CALL check( nf90_put_att(ncid,nf90_global,"n", nn))
      CALL check( nf90_put_att(ncid,nf90_global,"version", version))
      ! define global attributes
      ! define dimensions
      IF(debug_flag) PRINT *," - Defining dimensions in netcdf"

      WRITE(*,*)"netcdf msing=",msing
      !WRITE(*,*)"netcdf qval=",qval
      WRITE(*,*)"netcdf qval_arr=",qval_arr

      IF(msing>0)THEN
         CALL check( nf90_def_dim(ncid,"qsing",msing,qsing_dim) ) !r_dim = q_rational
         CALL check( nf90_def_dim(ncid, "i", 2, i_dim) )
         CALL check( nf90_def_var(ncid,"qsing",nf90_int,qsing_dim,
     $    qsing_id))
         CALL check( nf90_def_var(ncid,"omegas",nf90_double,
     $    qsing_dim,omegas_id))
         CALL check( nf90_def_var(ncid,"Q",nf90_double,
     $    qsing_dim,Q_id))
         CALL check( nf90_def_var(ncid,"Q_e",nf90_double,
     $    qsing_dim,Q_e_id))
         CALL check( nf90_def_var(ncid,"Q_i",nf90_double,
     $    qsing_dim,Q_i_id))
         CALL check( nf90_def_var(ncid,"S",nf90_double,
     $    qsing_dim,S_id))
         CALL check( nf90_def_var(ncid,"psi_n_rational",nf90_double,
     $                            qsing_dim,pr_id) )
         CALL check( nf90_def_var(ncid,"P",nf90_double,
     $                            qsing_dim,inpr_id) )
         CALL check( nf90_def_var(ncid,"q_rational",nf90_double,
     $                            qsing_dim,qr_id) )
         CALL check( nf90_def_dim(ncid, "points", max_points,
     $    point_dimid) )

         IF ((stabscan_eq_flag) .OR. (stabscan_flag)) THEN
            CALL check( nf90_def_var(ncid,"growthrates",nf90_double,
     $      qsing_dim,gamma_id))
            CALL check(nf90_def_var(ncid,"growthrate_locs",nf90_double,
     $      qsing_dim,gamma_loc_id))
            CALL check( nf90_def_var(ncid, "Re_Qs", nf90_double,
     $      [point_dimid, qsing_dim], varids(1)) )
            CALL check( nf90_def_var(ncid, "Im_Qs", nf90_double,
     $       [point_dimid, qsing_dim], varids(2)) )
            CALL check( nf90_def_var(ncid, "Re_deltas", nf90_double,
     $      [point_dimid, qsing_dim], varids(3)) )
            CALL check( nf90_def_var(ncid, "Im_deltas", nf90_double,
     $      [point_dimid, qsing_dim], varids(4)) )
         END IF

         IF (br_th_flag) THEN
            CALL check( nf90_def_var(ncid,"br_th",nf90_double,
     $                            qsing_dim,br_th_id) )
         END IF
      END IF

      IF ((stabscan_eq_flag) .OR. (stabscan_flag)) THEN
        CALL check( nf90_def_var(ncid,"Re_deltaprime",nf90_double,
     $      qsing_dim,rdpp_id) )
        CALL check( nf90_def_var(ncid,"Im_deltaprime",nf90_double,
     $      qsing_dim,idpp_id) )
      END IF
      IF ((lar_gamma_flag) .OR. (lar_gamma_eq_flag)) THEN
        CALL check( nf90_def_var(ncid,"delta_s_d_b",nf90_double,
     $      (/qsing_dim,i_dim/),dels_db_id) )
        CALL check( nf90_def_var(ncid,"d_beta",nf90_double,
     $      qsing_dim,d_b_id) )
        CALL check( nf90_def_var(ncid,"D_beta_norm",nf90_double,
     $      qsing_dim,inds_id) )
        CALL check( nf90_def_var(ncid,"growthrate",nf90_double,
     $      (/qsing_dim,i_dim/),lar_gamma_id) )
      END IF
      ! end definitions
      CALL check( nf90_enddef(ncid) )
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
      CALL check( nf90_put_var(ncid,qsing_id, qval_arr))
      CALL check( nf90_put_var(ncid,omegas_id, omegas_arr))
      CALL check( nf90_put_var(ncid,Q_id, inQ_arr))
      CALL check( nf90_put_var(ncid,Q_e_id, inQ_e_arr))
      CALL check( nf90_put_var(ncid,Q_i_id, inQ_i_arr))
      CALL check( nf90_put_var(ncid,S_id, (/lu/)))
      CALL check( nf90_put_var(ncid,pr_id, psi_n_rational))
      CALL check( nf90_put_var(ncid,inpr_id, inpr_arr))
      !CALL check( nf90_put_var(ncid,qr_id, qval_arr))

      IF ((stabscan_eq_flag) .OR. (stabscan_flag)) THEN
        CALL check( nf90_put_var(ncid,rdpp_id, Re_deltaprime_arr))
        CALL check( nf90_put_var(ncid,idpp_id, Im_deltaprime_arr))
        DO run = 1, msing
            CALL check( nf90_put_var(ncid,varids(1),results(run)%inQs,
     $       start=[1, run], count=[results(run)%count, 1]) )
            CALL check( nf90_put_var(ncid,varids(2),results(run)%iinQs,
     $       start=[1, run], count=[results(run)%count, 1]) )
            CALL check( nf90_put_var(ncid, varids(3),
     $       results(run)%Re_deltas, start=[1, run],
     $       count=[results(run)%count, 1]) )
            CALL check( nf90_put_var(ncid, varids(4),
     $       results(run)%Im_deltas, start=[1, run],
     $       count=[results(run)%count, 1]) )
        END DO
      END IF

      IF (br_th_flag) THEN
        CALL check( nf90_put_var(ncid,br_th_id, (/ br_th /)))
      END IF

      IF ((lar_gamma_flag) .OR. (lar_gamma_eq_flag)) THEN
        CALL check( nf90_put_var(ncid,dels_db_id, 
     $  RESHAPE((/REAL(dels_db),AIMAG(dels_db)/),(/qsing_dim,2/))))
        CALL check( nf90_put_var(ncid,d_b_id, (/ d_beta /)))
        CALL check( nf90_put_var(ncid,inds_id, (/ D_beta_norm /)))
        CALL check( nf90_put_var(ncid,lar_gamma_id, 
     $  RESHAPE((/REAL(lar_gamma),AIMAG(lar_gamma)/),(/qsing_dim,2/))))
      END IF

c -----------------------------------------------------------------------
c      close file
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Closing netcdf file"
      CALL check( nf90_close(ncid) )
c -----------------------------------------------------------------------
c      terminate.
c -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slayer_netcdf_out
c -----------------------------------------------------------------------
c      declarations.
c -----------------------------------------------------------------------
      SUBROUTINE slayer_netcdf_inputs(msing,qval_arr,ne_arr,te_arr,
     $  ni_arr,ti_arr,zeff_arr,shear,bt_arr,rs_arr,R0_arr,
     $  resm,nns_arr,inc_beta_arr,inds_arr,intau_arr,inpr_arr,inpe_arr,
     $  inQ_arr,omegas_arr,omegas_e_arr,omegas_i_arr,
     $  Re_deltaprime_arr,Im_deltaprime_arr)

      INTEGER, INTENT(IN) :: msing
      REAL(r8), DIMENSION(:), INTENT(IN) :: ne_arr,te_arr,ni_arr,
     $   ti_arr,zeff_arr,shear,bt_arr,rs_arr,
     $   R0_arr,inc_beta_arr,inds_arr,intau_arr,inpr_arr,
     $   inpe_arr,omegas_arr,omegas_e_arr,omegas_i_arr,
     $   Re_deltaprime_arr,Im_deltaprime_arr,inQ_arr
      INTEGER, DIMENSION(:), INTENT(IN) :: qval_arr,resm,nns_arr

      INTEGER :: i, ncid,r_id,qsing_dim,qsing_id,msing_id,
     $    i_dim,ne_id,te_id,ni_id,ti_id,zeff_id,shear_id,bt_id,rs_id,
     $    R0_id,resm_id,nns_id,inQ_id,inQ_e_id,inc_beta_id,
     $    inds_id,qval_id,inQ_i_id,qr_id,
     $    intau_id,inpr_id,inpe_id,omegas_id,Re_delta_id,Im_delta_id,
     $    omegas_e_id,omegas_i_id

      CHARACTER(64) :: ncfile
      LOGICAL, PARAMETER :: debug_flag = .FALSE.
      CHARACTER(len=*), PARAMETER :: version ='v1.0.0-99-gc873bd6'
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
      ne_id=0
      te_id=0
      ni_id=0
      ti_id=0
      zeff_id=0
      shear_id=0
      bt_id=0
      rs_id=0
      R0_id=0
      resm_id=0
      nns_id=0
      inQ_id=0
      inQ_e_id=0
      inc_beta_id=0
      inds_id=0
      qval_id=0
      inQ_i_id=0
      qr_id=0
      intau_id=0
      inpr_id=0
      inpe_id=0
      omegas_id=0
      Re_delta_id=0
      Im_delta_id=0
      omegas_e_id=0
      omegas_i_id=0

      IF(debug_flag) PRINT *,"Called slayer_netcdf_inputs"
      IF (nn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)')nn
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')nn
      ENDIF
      ncfile = "slayer_inputs_n"//TRIM(sn)//".nc"
      IF(debug_flag) PRINT *, ncfile
c -----------------------------------------------------------------------
c      open files
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Creating netcdf files"
      CALL check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
c -----------------------------------------------------------------------
c      define global file attributes
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining netcdf globals"
      CALL check( nf90_put_att(ncid,nf90_global,"title",
     $     "SLAYER outputs"))
      CALL check( nf90_put_att(ncid,nf90_global,"version", version))

      IF(debug_flag) PRINT *," - Defining dimensions in netcdf"

      IF(msing>0)THEN
         !CALL check( nf90_def_dim(ncid,"qsing",msing,qsing_dim) ) !r_dim = q_rational
         CALL check( nf90_def_var(ncid,"qsing",nf90_int,qsing_dim,
     $    qsing_id))

         CALL check( nf90_def_var(ncid,"ne",nf90_double,
     $    qsing_dim,ne_id))
         CALL check( nf90_def_var(ncid,"te",nf90_double,
     $    qsing_dim,te_id))
         CALL check( nf90_def_var(ncid,"ni",nf90_double,
     $    qsing_dim,ni_id))
         CALL check( nf90_def_var(ncid,"ti",nf90_double,
     $    qsing_dim,ti_id))
         CALL check( nf90_def_var(ncid,"zeff",nf90_double,
     $    qsing_dim,zeff_id))
         CALL check( nf90_def_var(ncid,"shear",nf90_double,
     $    qsing_dim,shear_id))
         CALL check( nf90_def_var(ncid,"bt",nf90_double,
     $    qsing_dim,bt_id))
         CALL check( nf90_def_var(ncid,"rs",nf90_double,
     $    qsing_dim,rs_id))
         CALL check( nf90_def_var(ncid,"R0",nf90_double,
     $    qsing_dim,R0_id))
         !CALL check( nf90_def_var(ncid,"mu_i",nf90_double,
      !$    qsing_dim,mu_i_id))
         CALL check( nf90_def_var(ncid,"resm",nf90_int,
     $    qsing_dim,resm_id))
         CALL check( nf90_def_var(ncid,"nns_arr",nf90_int,
     $    qsing_dim,nns_id))
         !CALL check( nf90_def_var(ncid,"qval",nf90_int,
      !$    qsing_dim,qval_id))
         CALL check( nf90_def_var(ncid,"Q",nf90_double,
     $    qsing_dim,inQ_id))
       !  CALL check( nf90_def_var(ncid,"Q_e",nf90_double,
      !$    qsing_dim,inQ_e_id))
       !  CALL check( nf90_def_var(ncid,"Q_i",nf90_double,
       !$    qsing_dim,inQ_i_id))
         CALL check( nf90_def_var(ncid,"c_beta",nf90_double,
     $    qsing_dim,inc_beta_id))
         CALL check( nf90_def_var(ncid,"ds",nf90_double,
     $    qsing_dim,inds_id))
         CALL check( nf90_def_var(ncid,"tau",nf90_double,
     $    qsing_dim,intau_id))
         CALL check( nf90_def_var(ncid,"pr",nf90_double,
     $    qsing_dim,inpr_id))
         CALL check( nf90_def_var(ncid,"pe",nf90_double,
     $    qsing_dim,inpe_id))
         CALL check( nf90_def_var(ncid,"omegas",nf90_double,
     $    qsing_dim,omegas_id))
         CALL check( nf90_def_var(ncid,"omegas_e",nf90_double,
     $    qsing_dim,omegas_e_id))
         CALL check( nf90_def_var(ncid,"omegas_i",nf90_double,
     $    qsing_dim,omegas_i_id))
         CALL check( nf90_def_var(ncid,"Re_deltaprime",nf90_double,
     $    qsing_dim,Re_delta_id))
         CALL check( nf90_def_var(ncid,"Im_deltaprime",nf90_double,
     $    qsing_dim,Im_delta_id))
      ENDIF
      ! define variables
      IF(debug_flag) PRINT *," - Defining variables in netcdf"
      ! end definitions
      CALL check( nf90_enddef(ncid) )
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
 !     IF(debug_flag) PRINT *," - Putting profile variables in netcdf"
      CALL check( nf90_put_var(ncid,qsing_id, qval_arr))
      CALL check( nf90_put_var(ncid,ne_id, ne_arr))
      CALL check( nf90_put_var(ncid,ni_id, ni_arr))
      CALL check( nf90_put_var(ncid,te_id, te_arr))
      CALL check( nf90_put_var(ncid,ti_id, ti_arr))
      CALL check( nf90_put_var(ncid,zeff_id, zeff_arr))
      CALL check( nf90_put_var(ncid,shear_id, shear))
      CALL check( nf90_put_var(ncid,bt_id, bt_arr))
      CALL check( nf90_put_var(ncid,rs_id, rs_arr))
      CALL check( nf90_put_var(ncid,R0_id, R0_arr))
      !CALL check( nf90_put_var(ncid,mu_i_id, mu_i_arr))
      CALL check( nf90_put_var(ncid,resm_id, resm))
      CALL check( nf90_put_var(ncid,nns_id, nns_arr))
      !CALL check( nf90_put_var(ncid,qval_id, qval_arr))
      CALL check( nf90_put_var(ncid,inQ_id, inQ_arr))
      !CALL check( nf90_put_var(ncid,inQ_e_id, inQ_e_arr))
      !CALL check( nf90_put_var(ncid,inQ_i_id, inQ_i_arr))
      CALL check( nf90_put_var(ncid,inc_beta_id, inc_beta_arr))
      CALL check( nf90_put_var(ncid,inds_id, inds_arr))
      CALL check( nf90_put_var(ncid,intau_id, intau_arr))
      CALL check( nf90_put_var(ncid,inpr_id, inpr_arr))
      CALL check( nf90_put_var(ncid,inpe_id, inpe_arr))
      CALL check( nf90_put_var(ncid,omegas_id, omegas_arr))
      CALL check( nf90_put_var(ncid,omegas_e_id, omegas_e_arr))
      CALL check( nf90_put_var(ncid,omegas_i_id, omegas_i_arr))
      CALL check( nf90_put_var(ncid,Re_delta_id,Re_deltaprime_arr))
      CALL check( nf90_put_var(ncid,Im_delta_id,Im_deltaprime_arr))

c -----------------------------------------------------------------------
c      close file
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Closing netcdf file"
      CALL check( nf90_close(ncid) )
c -----------------------------------------------------------------------
c      terminate.
c -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slayer_netcdf_inputs
      END MODULE slayer_netcdf_mod