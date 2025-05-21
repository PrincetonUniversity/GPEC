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
c-----------------------------------------------------------------------
c     subprogram 1. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE sl_check(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sl_check
c -----------------------------------------------------------------------
c      subprogram 2. slayer_netcdf_out.
c      Replicate stride.out information in netcdf format.
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c      declarations.
c -----------------------------------------------------------------------
      SUBROUTINE slayer_netcdf_out(msing,est_gamma_flag,
     $         qval_arr,omegas_arr,Q_arr,Q_e_arr,Q_i_arr,d_beta_arr,
     $         c_beta_arr,D_norm_arr,P_perp_arr,lu_arr,psi_n_rational,
     $         Re_deltaprime_arr,Im_deltaprime_arr,delta_crit_arr,
     $         dels_db_arr,gamma_sol_arr,gamma_est_arr,Qconv_arr,
     $         re_trace,im_trace)

      INTEGER, INTENT(IN) :: msing
      LOGICAL, INTENT(IN) :: est_gamma_flag
      INTEGER, INTENT(IN), DIMENSION(:) :: qval_arr
      REAL(r8), INTENT(IN), DIMENSION(:) :: omegas_arr,
     $      Q_arr,Q_e_arr,Q_i_arr,d_beta_arr,c_beta_arr,D_norm_arr,
     $      P_perp_arr,lu_arr,psi_n_rational,Re_deltaprime_arr,
     $      Im_deltaprime_arr,delta_crit_arr,re_trace,im_trace,
     $      Qconv_arr
      COMPLEX(r8),INTENT(IN),DIMENSION(:) :: dels_db_arr,
     $                                    gamma_sol_arr,gamma_est_arr

      INTEGER :: i,ncid,r_id,qsing_dim,i_dim,r_dim,qr_id,omegas_id,
     $    Q_id,Q_e_id,Q_i_id,d_b_id,c_b_id,Dnorm_id,inpr_id,S_id,
     $    pr_id,dpp_id,dc_id,dels_db_id,gs_id,ge_id,tr_id,tr_dim,
     $    qsing_id,qc_id

      INTEGER :: run, run_dimid, point_dimid, varids(4)

      CHARACTER(64) :: ncfile
      LOGICAL, PARAMETER :: debug_flag = .FALSE.
      CHARACTER(len=*), PARAMETER :: version ='v1.0.0-99-gc873bd6'
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
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
      CALL sl_check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

c -----------------------------------------------------------------------
c      define global file attributes
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining netcdf globals"
      CALL sl_check( nf90_put_att(ncid,nf90_global,"title",
     $     "SLAYER outputs"))
      !CALL sl_check( nf90_put_att(ncid,nf90_global,"shot", INT(shotnum)) )
      !CALL sl_check( nf90_put_att(ncid,nf90_global,"time",INT(shottime)) )
      !CALL sl_check( nf90_put_att(ncid,nf90_global,"n", nn))
      CALL sl_check( nf90_put_att(ncid,nf90_global,"version", version))
      ! define global attributes
      ! define dimensions
      IF(debug_flag) PRINT *," - Defining dimensions in netcdf"

      !WRITE(*,*)"netcdf qval=",qval
      WRITE(*,*)">>> Writing results to NetCDF output file"

      IF(msing>0)THEN
         CALL sl_check( nf90_def_dim(ncid,"qsing",msing,qsing_dim) ) !r_dim = q_rational
         CALL sl_check( nf90_def_dim(ncid, "i", 2, i_dim) )
         CALL sl_check( nf90_def_var(ncid,"qsing",nf90_int,qsing_dim,
     $    qsing_id))
         CALL sl_check( nf90_def_dim(ncid, "step", SIZE(re_trace), 
     $                  tr_dim) )
         CALL sl_check( nf90_def_var(ncid,"omegas",nf90_double,
     $    qsing_dim,omegas_id))
         CALL sl_check( nf90_def_var(ncid,"tau_k",nf90_double,
     $    qsing_dim,qc_id))
         CALL sl_check( nf90_def_var(ncid,"Q",nf90_double,
     $    qsing_dim,Q_id))
         CALL sl_check( nf90_def_var(ncid,"Q_e",nf90_double,
     $    qsing_dim,Q_e_id))
         CALL sl_check( nf90_def_var(ncid,"Q_i",nf90_double,
     $    qsing_dim,Q_i_id))
         CALL sl_check( nf90_def_var(ncid,"S",nf90_double,
     $    qsing_dim,S_id))
         CALL sl_check( nf90_def_var(ncid,"psi_n_rational",
     $                            nf90_double,qsing_dim,pr_id) )
         CALL sl_check( nf90_def_var(ncid,"P",nf90_double,
     $                            qsing_dim,inpr_id) )
         CALL sl_check( nf90_def_var(ncid,"q_rational",nf90_double,
     $                            qsing_dim,qr_id) )
      END IF

      CALL sl_check( nf90_def_var(ncid,"D",nf90_double,
     $      qsing_dim,Dnorm_id) )
      CALL sl_check( nf90_def_var(ncid,"Delta_prime_rational",
     $      nf90_double,(/qsing_dim,i_dim/),dpp_id) )
      CALL sl_check( nf90_def_var(ncid,"Delta_crit_rational",
     $      nf90_double,qsing_dim,dc_id) )

      IF (est_gamma_flag) THEN
        CALL sl_check( nf90_def_var(ncid,"delta_s_d_b",nf90_double,
     $      (/qsing_dim,i_dim/),dels_db_id) )
        CALL sl_check( nf90_def_var(ncid,"d_beta",nf90_double,
     $      qsing_dim,d_b_id) )
        CALL sl_check( nf90_def_var(ncid,"est. growth rate",
     $      nf90_double,(/qsing_dim,i_dim/),ge_id) )
      END IF

      CALL sl_check( nf90_def_var(ncid,"growth rate",
     $      nf90_double,(/qsing_dim,i_dim/),gs_id) )
      CALL sl_check( nf90_def_var(ncid,"growth rate trace",
     $      nf90_double,(/tr_dim,i_dim/),tr_id) )
      ! end definitions
      CALL sl_check( nf90_enddef(ncid) )
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
      CALL sl_check( nf90_put_var(ncid,qsing_id, qval_arr))
      CALL sl_check( nf90_put_var(ncid,qr_id, qval_arr))
      CALL sl_check( nf90_put_var(ncid,pr_id, psi_n_rational))
      CALL sl_check( nf90_put_var(ncid,omegas_id, omegas_arr))
      CALL sl_check( nf90_put_var(ncid,S_id, lu_arr))
      CALL sl_check( nf90_put_var(ncid,qc_id, Qconv_arr))
      CALL sl_check( nf90_put_var(ncid,Q_id, Q_arr))
      CALL sl_check( nf90_put_var(ncid,Q_e_id, Q_e_arr))
      CALL sl_check( nf90_put_var(ncid,Q_i_id, Q_i_arr))
      CALL sl_check( nf90_put_var(ncid,inpr_id, P_perp_arr))
      CALL sl_check( nf90_put_var(ncid,Dnorm_id, D_norm_arr))

      CALL sl_check( nf90_put_var(ncid,dpp_id, 
     $      RESHAPE((/Re_deltaprime_arr,Im_deltaprime_arr/),
     $      (/msing,2/))))
      CALL sl_check( nf90_put_var(ncid,dc_id, delta_crit_arr))

      IF (est_gamma_flag) THEN
        CALL sl_check( nf90_put_var(ncid,dels_db_id, 
     $      RESHAPE((/REAL(dels_db_arr),AIMAG(dels_db_arr)/),
     $      (/msing,2/))))
        CALL sl_check( nf90_put_var(ncid,d_b_id,d_beta_arr))
        CALL sl_check( nf90_put_var(ncid,ge_id, 
     $      RESHAPE((/REAL(gamma_est_arr),AIMAG(gamma_est_arr)/),
     $      (/msing,2/))))
      END IF

      CALL sl_check( nf90_put_var(ncid,gs_id, 
     $      RESHAPE((/REAL(gamma_sol_arr),AIMAG(gamma_sol_arr)/),
     $      (/msing,2/))))

      CALL sl_check( nf90_put_var(ncid,tr_id, 
     $      RESHAPE((/re_trace,im_trace/),
     $      (/SIZE(re_trace),2/))))

c -----------------------------------------------------------------------
c      close file
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Closing netcdf file"
      CALL sl_check( nf90_close(ncid) )
c -----------------------------------------------------------------------
c      terminate.
c -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slayer_netcdf_out
      END MODULE slayer_netcdf_mod