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
      !USE dcon_mod
      USE sglobal_mod
      !USE layerinputs_mod
      !USE gslayer_mod
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
      SUBROUTINE slayer_netcdf_out(msing,ReQ_n,ImQ_n,qval_arr,
     $ inQs_log,iinQs,growthrates,omegas_arr,inQ_arr,psi_n_rational,
     $ all_Re_deltas,all_slices)
      INTEGER, INTENT(IN) :: msing,ReQ_n,ImQ_n
      REAL(r8), DIMENSION(:), INTENT(IN) :: qval_arr,
     $ inQs_log,iinQs,growthrates,omegas_arr,inQ_arr,
     $ psi_n_rational!,shear
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: all_slices
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: all_Re_deltas
      INTEGER :: i, ncid,r_id,ReQ_dim,ImQ_dim,qsing_dim,qsing_id,
     $    i_dim, m_dim, mo_dim, p_dim, i_id, m_id, mo_id, p_id,
     $    ReQ_id,ImQ_id,gamma_id,omegas_id,Q_id,
     $    r_dim,pr_id, qr_id, dp_id,shear_id,slice_id
      !COMPLEX(r8), DIMENSION(mpert) :: ep,ev,et
      !CHARACTER(2) :: sn
      CHARACTER(64) :: ncfile
 !     INTEGER :: ising,jsing
      !COMPLEX(r8), DIMENSION(msing,msing) :: ap,bp,gammap,deltap
      LOGICAL, PARAMETER :: debug_flag = .FALSE.
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
      CALL check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
c -----------------------------------------------------------------------
c      define global file attributes
c -----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining netcdf globals"
      CALL check( nf90_put_att(ncid,nf90_global,"title",
     $     "SLAYER outputs"))
      ! define global attributes
      !CALL check( nf90_put_att(ncid,nf90_global,'ro', ro))
      !CALL check( nf90_put_att(ncid,nf90_global,'psio', psio))
      !CALL check( nf90_put_att(ncid,nf90_global,'bt0', bt0))
      !CALL check( nf90_put_att(ncid,nf90_global,"shot", INT(shotnum)) )
      !CALL check( nf90_put_att(ncid,nf90_global,"time",INT(shottime)) )
      !CALL check( nf90_put_att(ncid,nf90_global,"n", nn))
      !CALL check( nf90_put_att(ncid,nf90_global,"version", version))
      ! define dimensions
      IF(debug_flag) PRINT *," - Defining dimensions in netcdf"
      !CALL check( nf90_def_dim(ncid, "i", 2, i_dim) )
      !CALL check( nf90_def_var(ncid, "i", nf90_int, i_dim, i_id) )
      !CALL check( nf90_def_dim(ncid, "m", mpert,  m_dim) )
      !CALL check( nf90_def_var(ncid, "m", nf90_int, m_dim, m_id) )
      !CALL check( nf90_def_dim(ncid, "mode", mpert, mo_dim) )
      !CALL check( nf90_def_var(ncid, "mode", nf90_int, mo_dim, mo_id))
      !CALL check( nf90_def_dim(ncid, "psi_n", sq%mx+1, p_dim) )
      !CALL check( nf90_def_var(ncid, "psi_n", nf90_double, p_dim, p_id))
      IF(msing>0)THEN
         CALL check( nf90_def_dim(ncid,"qsing",msing,qsing_dim) ) !r_dim = q_rational
         CALL check( nf90_def_var(ncid,"qsing",nf90_int,qsing_dim,
     $    qsing_id))
         CALL check( nf90_def_dim(ncid,"ReQ_arr",ReQ_n,
     $    ReQ_dim) ) !r_dim = q_rational
         CALL check( nf90_def_var(ncid,"ReQ_arr",nf90_double,ReQ_dim,
     $    ReQ_id))
         CALL check( nf90_def_dim(ncid,"ImQ_arr",ImQ_n,ImQ_dim) ) !r_dim = q_rational
         CALL check( nf90_def_var(ncid,"ImQ_arr",nf90_double,ImQ_dim,
     $    ImQ_id))
         CALL check( nf90_def_var(ncid,"growthrates",nf90_double,
     $    qsing_dim,gamma_id))
         CALL check( nf90_def_var(ncid,"omegas",nf90_double,
     $    qsing_dim,omegas_id))
         CALL check( nf90_def_var(ncid,"Q",nf90_double,
     $    qsing_dim,Q_id))
         CALL check( nf90_def_var(ncid,"psi_n_rational",nf90_double,
     $                            qsing_dim,pr_id) )
         CALL check( nf90_def_var(ncid,"q_rational",nf90_double,
     $                            qsing_dim,qr_id) )
      !   CALL check( nf90_def_var(ncid,"shear",nf90_double,q_rational,
      !$                            shear_id) )
      ENDIF
      ! define variables
      IF(debug_flag) PRINT *," - Defining variables in netcdf"
      !CALL check( nf90_def_var(ncid, "f", nf90_double, p_dim, f_id) )
      !CALL check( nf90_def_var(ncid, "mu0p", nf90_double, p_dim, mu_id))
      !CALL check( nf90_def_var(ncid, "dvdpsi", nf90_double,p_dim,dv_id))
      !CALL check( nf90_def_var(ncid, "q", nf90_double, p_dim, q_id) )
      !CALL check( nf90_def_var(ncid, "di", nf90_double, p_dim, di_id) )
      !CALL check( nf90_def_var(ncid, "dr", nf90_double, p_dim, dr_id) )
      !CALL check( nf90_def_var(ncid, "ca1", nf90_double, p_dim, ca_id))
      !CALL check( nf90_def_var(ncid, "W_p_eigenvector", nf90_double,
      !$    (/m_dim, mo_dim, i_dim/), wp_id) )
      !CALL check( nf90_def_var(ncid, "W_p_eigenvalue", nf90_double,
      !$     (/mo_dim, i_dim/), wpv_id) )
      IF(msing>0)THEN
         CALL check( nf90_def_var(ncid, "Re_Delta", nf90_double,
     $       (/ReQ_dim, ImQ_dim, qsing_dim/), dp_id) )
         CALL check( nf90_def_var(ncid, "slices", nf90_double,
     $       (/ImQ_dim, qsing_dim/), slice_id) )
      ENDIF
      ! end definitions
      CALL check( nf90_enddef(ncid) )
c -----------------------------------------------------------------------
c      set variables
c -----------------------------------------------------------------------
 !     IF(debug_flag) PRINT *," - Putting profile variables in netcdf"
      CALL check( nf90_put_var(ncid,qsing_id, qval_arr))
      CALL check( nf90_put_var(ncid,ReQ_id, inQs_log))
      CALL check( nf90_put_var(ncid,ImQ_id, iinQs))
      CALL check( nf90_put_var(ncid,gamma_id, growthrates))
      CALL check( nf90_put_var(ncid,omegas_id, omegas_arr))
      CALL check( nf90_put_var(ncid,Q_id, inQ_arr))
      CALL check( nf90_put_var(ncid,pr_id, psi_n_rational))
      CALL check( nf90_put_var(ncid,qr_id, qval_arr))
      !CALL check( nf90_put_var(ncid,shear_id, shear))
 !     IF(debug_flag) PRINT *," - Putting matrix variables in netcdf"
      CALL check( nf90_put_var(ncid,dp_id,all_Re_deltas))!,
      !$(/ReQ_n,ImQ_n,msing/)))
      CALL check( nf90_put_var(ncid,slice_id,all_slices))!,
      !$(/ImQ_n,msing/)))
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
      END MODULE slayer_netcdf_mod