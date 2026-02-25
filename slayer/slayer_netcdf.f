c=======================================================================
c     file slayer_netcdf.f
c     Writes SLAYER solver results to a NetCDF output file.
c
c     The per-surface scalar inputs (Lundquist number, Q-normalisation,
c     Prandtl numbers, …) and the solver outputs (growth rates, Delta
c     values) are stored in a single NetCDF-3/64-bit-offset file named
c     slayer_output_n<n>.nc.
c
c     Ragged AMR scan data (variable number of evaluation points per
c     surface) are zero-padded into rectangular arrays before writing.
c
c     (Resolved: FLAG 1 -- early RETURN for msing == 0.
c      FLAG 2 -- buffers now use fill_val; _FillValue attributes added.
c      FLAG 3 -- Q_id removed (unused).
c      FLAG 4 -- c_b_id removed (unused).
c      FLAG 5 -- stale unused-variable comment removed.
c      FLAG 6 -- version from INCLUDE "version.inc".
c      FLAG 7 -- local sn_local replaces global sn_str.)
c=======================================================================
c-----------------------------------------------------------------------
c     code organisation.
c-----------------------------------------------------------------------
c     0. slayer_netcdf_mod   -- module declarations
c     1. sl_check            -- NetCDF status checker
c     2. slayer_netcdf_out   -- main output routine
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     subprogram 0. slayer_netcdf_mod.
c     Module wrapper — imports sglobal_mod (shared types and globals)
c     and the NetCDF Fortran-90 API.
c-----------------------------------------------------------------------
      MODULE slayer_netcdf_mod

      USE sglobal_mod
      USE netcdf
      USE ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan

      IMPLICIT NONE

      CONTAINS
c
c-----------------------------------------------------------------------
c     subprogram 1. sl_check.
c     Assert that a NetCDF operation succeeded; abort with a message
c     if it did not.
c-----------------------------------------------------------------------
      SUBROUTINE sl_check(stat)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: stat    ! return code from any nf90_* call
c-----------------------------------------------------------------------
c     check status and abort on error.
c-----------------------------------------------------------------------
      IF (stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF

      RETURN
      END SUBROUTINE sl_check
c
c-----------------------------------------------------------------------
c     subprogram 2. slayer_netcdf_out.
c     Write per-surface SLAYER inputs and solver outputs to the
c     NetCDF file  slayer_output_n<n>.nc .
c
c     Arguments:
c       msing           -- number of rational surfaces
c       m_AMR           -- number of AMR-scanned surfaces
c       est_gamma_flag  -- .TRUE. to include estimated growth rates
c       sl_in           -- slayer_inputs_type  (per-surface inputs)
c       sl_out          -- slayer_outputs_type (solver results)
c       all_deltas_out  -- array(m_AMR) of deltas_outputs_type (AMR
c                          scan results, potentially ragged)
c-----------------------------------------------------------------------
      SUBROUTINE slayer_netcdf_out(msing, m_AMR, est_gamma_flag,
     $                             sl_in, sl_out, all_deltas_out)
c-----------------------------------------------------------------------
c     declarations -- subroutine arguments.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: msing          ! number of rational surfaces
      INTEGER, INTENT(IN) :: m_AMR          ! number of AMR surfaces
      LOGICAL, INTENT(IN) :: est_gamma_flag ! include estimated gammas?

      TYPE(slayer_inputs_type),  INTENT(IN) :: sl_in
      TYPE(slayer_outputs_type), INTENT(IN) :: sl_out
      TYPE(deltas_outputs_type), INTENT(IN) :: all_deltas_out(m_AMR)

c-----------------------------------------------------------------------
c     declarations -- NetCDF file and dimension IDs.
c-----------------------------------------------------------------------
      INTEGER :: ncid              ! NetCDF file ID
      INTEGER :: qsing_dim        ! dim: rational surfaces  (msing)
      INTEGER :: nAMR_dim         ! dim: AMR surfaces       (m_AMR)
      INTEGER :: i_dim            ! dim: Re/Im component    (2)
      INTEGER :: dim_pts_id       ! dim: max AMR eval points

c-----------------------------------------------------------------------
c     declarations -- NetCDF variable IDs.
c     Each *_id holds the handle returned by nf90_def_var and is
c     later passed to the matching nf90_put_var call.
c-----------------------------------------------------------------------
      INTEGER :: qsing_id         ! "r"                  — surface index
      INTEGER :: qr_id            ! "q_rational"         — safety factor
      INTEGER :: omegas_id        ! "omegas"             — rotation freq
      INTEGER :: qc_id            ! "tau_k"              — Q-conversion
      INTEGER :: Q_e_id           ! "Q_e"                — norm Q_e
      INTEGER :: Q_i_id           ! "Q_i"                — norm Q_i
      INTEGER :: S_id             ! "S"                  — Lundquist
      INTEGER :: pr_id            ! "psi_n_rational"     — norm psi
      INTEGER :: p_perp_id        ! "P_perp"
      INTEGER :: p_tor_id         ! "P_tor"
      INTEGER :: Dnorm_id         ! "D"                  — normalised D
      INTEGER :: dpp_id           ! "Delta_prime_rational" (complex)
      INTEGER :: dc_id            ! "Delta_crit_rational"
      INTEGER :: dels_db_id       ! "delta_s_d_b"        (complex)
      INTEGER :: d_b_id           ! "d_beta"
      INTEGER :: gs_id            ! "growth rate"        (complex)
      INTEGER :: ge_id            ! "est. growth rate"   (complex)
      INTEGER :: br_th_id         ! "br_th"              — Br threshold

c     AMR variable IDs
      INTEGER :: var_q_id         ! "Q_AMR"    — scan Q-points
      INTEGER :: var_d_id         ! "Deltas_AMR" — scan Delta values
      INTEGER :: var_npts_id      ! "n_amr_pts"  — points per surface

c-----------------------------------------------------------------------
c     declarations -- AMR rectangular-buffer workspace.
c     The ragged per-surface scan data are padded into fixed-size
c     rectangular arrays before writing to NetCDF.
c-----------------------------------------------------------------------
      INTEGER :: max_pts_all             ! max points across surfaces
      INTEGER :: s                       ! surface loop index
      INTEGER :: n_curr                  ! points on current surface
      REAL(r8), ALLOCATABLE :: buffer_q(:,:,:)  ! (pts, surf, Re/Im)
      REAL(r8), ALLOCATABLE :: buffer_d(:,:,:)  ! (pts, surf, Re/Im)
      INTEGER,  ALLOCATABLE :: n_pts_arr(:)     ! points per surface
      REAL(r8) :: fill_val               ! NaN padding for ragged arrays

c-----------------------------------------------------------------------
c     declarations -- miscellaneous locals.
c-----------------------------------------------------------------------
      CHARACTER(64) :: ncfile                      ! output file name
      CHARACTER(2)  :: sn_local                     ! local n-string for filename
      LOGICAL, PARAMETER :: debug_flag = .FALSE.   ! verbose trace
      INCLUDE "version.inc"

c-----------------------------------------------------------------------
c     build the output filename from the toroidal mode number.
c-----------------------------------------------------------------------
      IF (debug_flag) PRINT *, "Called slayer_netcdf_out"

      IF (nn < 10) THEN
         WRITE(UNIT=sn_local, FMT='(I1)') nn
         sn_local = ADJUSTL(sn_local)
      ELSE
         WRITE(UNIT=sn_local, FMT='(I2)') nn
      ENDIF
      ncfile = "slayer_output_n"//TRIM(sn_local)//".nc"
      IF (debug_flag) PRINT *, ncfile

c-----------------------------------------------------------------------
c     create the NetCDF file (clobber any existing file).
c-----------------------------------------------------------------------
      IF (debug_flag) PRINT *, " - Creating netcdf file"
      CALL sl_check( nf90_create(ncfile,
     $     cmode=OR(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid=ncid) )

c-----------------------------------------------------------------------
c     reform ragged AMR Delta outputs into rectangular arrays.
c
c     Each surface may have a different number of AMR scan points.
c     We find the maximum, allocate rectangular buffers of that size,
c     and copy in the per-surface data.  Unused trailing slots are
c     filled with fill_val (-9.99E33).
c-----------------------------------------------------------------------

c     step 1: find the maximum AMR grid size across all surfaces.
      max_pts_all = 0
      IF (ALLOCATED(all_deltas_out(1)%inQs)) THEN
         DO s = 1, m_AMR
            max_pts_all = MAX(max_pts_all,
     $                        SIZE(all_deltas_out(s)%inQs))
         END DO
      END IF

c     step 2: allocate rectangular buffers (pts × surfaces × Re/Im).
      ALLOCATE(buffer_q(max_pts_all, m_AMR, 2))
      ALLOCATE(buffer_d(max_pts_all, m_AMR, 2))
      ALLOCATE(n_pts_arr(m_AMR))

      fill_val  = ieee_value(1.0d0, ieee_quiet_nan)
      buffer_q  = fill_val        ! NaN padding for ragged arrays
      buffer_d  = fill_val
      n_pts_arr = 0

c     step 3: flatten the ragged data into the buffers.
      IF (ALLOCATED(all_deltas_out(1)%inQs)) THEN
         DO s = 1, m_AMR
            n_curr       = SIZE(all_deltas_out(s)%inQs)
            n_pts_arr(s) = n_curr

c           Q-coordinate (Re and Im parts)
            buffer_q(1:n_curr, s, 1) =
     $           all_deltas_out(s)%inQs(1:n_curr)
            buffer_q(1:n_curr, s, 2) =
     $           all_deltas_out(s)%iinQs(1:n_curr)

c           Delta result (Re and Im parts)
            buffer_d(1:n_curr, s, 1) =
     $           all_deltas_out(s)%real_deltas(1:n_curr)
            buffer_d(1:n_curr, s, 2) =
     $           all_deltas_out(s)%imag_deltas(1:n_curr)
         END DO
      END IF

c-----------------------------------------------------------------------
c     define global file attributes.
c-----------------------------------------------------------------------
      IF (debug_flag) PRINT *, " - Defining netcdf globals"
      CALL sl_check( nf90_put_att(ncid, nf90_global,
     $                            "title", "SLAYER outputs") )
      CALL sl_check( nf90_put_att(ncid, nf90_global,
     $                            "version", version) )

c-----------------------------------------------------------------------
c     define dimensions and per-surface NetCDF variables.
c-----------------------------------------------------------------------
      IF (debug_flag) PRINT *, " - Defining dimensions in netcdf"
      WRITE(*,*) ">>> Writing results to NetCDF output file"

      IF (msing == 0) THEN
         WRITE(*,*) "WARNING: msing == 0, skipping NetCDF output"
         DEALLOCATE(buffer_q, buffer_d, n_pts_arr)
         CALL sl_check( nf90_close(ncid) )
         RETURN
      END IF

c     -- core dimensions --
      CALL sl_check( nf90_def_dim(ncid, "r",     msing, qsing_dim) )
      CALL sl_check( nf90_def_dim(ncid, "r_AMR", m_AMR, nAMR_dim)  )
      CALL sl_check( nf90_def_dim(ncid, "i",     2,     i_dim)     )

c     -- scalar per-surface variables --
      CALL sl_check( nf90_def_var(ncid, "r",         nf90_int,
     $     qsing_dim, qsing_id)  )
      CALL sl_check( nf90_def_var(ncid, "q_rational",
     $     nf90_double, qsing_dim, qr_id)     )
      CALL sl_check( nf90_def_var(ncid, "omegas",    nf90_double,
     $     qsing_dim, omegas_id) )
      CALL sl_check( nf90_def_var(ncid, "tau_k",     nf90_double,
     $     qsing_dim, qc_id)     )
      CALL sl_check( nf90_def_var(ncid, "Q_e",       nf90_double,
     $     qsing_dim, Q_e_id)    )
      CALL sl_check( nf90_def_var(ncid, "Q_i",       nf90_double,
     $     qsing_dim, Q_i_id)    )
      CALL sl_check( nf90_def_var(ncid, "S",         nf90_double,
     $     qsing_dim, S_id)      )
      CALL sl_check( nf90_def_var(ncid, "psi_n_rational",
     $     nf90_double, qsing_dim, pr_id)     )
      CALL sl_check( nf90_def_var(ncid, "P_perp",    nf90_double,
     $     qsing_dim, p_perp_id) )
      CALL sl_check( nf90_def_var(ncid, "P_tor",     nf90_double,
     $     qsing_dim, p_tor_id)  )

c-----------------------------------------------------------------------
c     define additional variables (D, Delta', growth rates).
c-----------------------------------------------------------------------
      CALL sl_check( nf90_def_var(ncid, "D", nf90_double,
     $     qsing_dim, Dnorm_id) )
      CALL sl_check( nf90_def_var(ncid, "Delta_prime_rational",
     $     nf90_double, (/qsing_dim, i_dim/), dpp_id) )
      CALL sl_check( nf90_def_var(ncid, "Delta_crit_rational",
     $     nf90_double, qsing_dim, dc_id) )

      IF (est_gamma_flag) THEN
         CALL sl_check( nf90_def_var(ncid, "delta_s_d_b",
     $        nf90_double, (/qsing_dim, i_dim/), dels_db_id) )
         CALL sl_check( nf90_def_var(ncid, "d_beta", nf90_double,
     $        qsing_dim, d_b_id) )
         CALL sl_check( nf90_def_var(ncid, "est. growth rate",
     $        nf90_double, (/qsing_dim, i_dim/), ge_id) )
      END IF

      CALL sl_check( nf90_def_var(ncid, "growth rate",
     $     nf90_double, (/qsing_dim, i_dim/), gs_id) )

      IF (ALLOCATED(sl_out%br_th_arr)) THEN
         CALL sl_check( nf90_def_var(ncid, "br_th",
     $        nf90_double, qsing_dim, br_th_id) )
      END IF

c-----------------------------------------------------------------------
c     define AMR scan dimensions and variables.
c-----------------------------------------------------------------------
      CALL sl_check( nf90_def_dim(ncid, "amr_pts",
     $     max_pts_all, dim_pts_id) )
      CALL sl_check( nf90_def_var(ncid, "n_amr_pts", NF90_INT,
     $     (/nAMR_dim/), var_npts_id) )
      CALL sl_check( nf90_def_var(ncid, "Q_AMR", NF90_DOUBLE,
     $     (/dim_pts_id, nAMR_dim, i_dim/), var_q_id) )
      CALL sl_check( nf90_def_var(ncid, "Deltas_AMR", NF90_DOUBLE,
     $     (/dim_pts_id, nAMR_dim, i_dim/), var_d_id) )

c     set _FillValue attribute on AMR arrays for proper padding.
      CALL sl_check( nf90_put_att(ncid, var_q_id,
     $     "_FillValue", fill_val) )
      CALL sl_check( nf90_put_att(ncid, var_d_id,
     $     "_FillValue", fill_val) )

c-----------------------------------------------------------------------
c     end NetCDF define mode.
c-----------------------------------------------------------------------
      CALL sl_check( nf90_enddef(ncid) )

c-----------------------------------------------------------------------
c     write per-surface scalar variables.
c-----------------------------------------------------------------------
      CALL sl_check( nf90_put_var(ncid, qsing_id,
     $     sl_in%qval_arr)  )
      CALL sl_check( nf90_put_var(ncid, qr_id,
     $     sl_in%qval_arr)  )
      CALL sl_check( nf90_put_var(ncid, pr_id,
     $     sl_in%psi_n_arr) )
      CALL sl_check( nf90_put_var(ncid, omegas_id,
     $     sl_in%omegas_arr))
      CALL sl_check( nf90_put_var(ncid, S_id,
     $     sl_in%lu_arr)    )
      CALL sl_check( nf90_put_var(ncid, qc_id,
     $     sl_in%Qconv_arr) )
      CALL sl_check( nf90_put_var(ncid, Q_e_id,
     $     sl_in%Q_e_arr)   )
      CALL sl_check( nf90_put_var(ncid, Q_i_id,
     $     sl_in%Q_i_arr)   )
      CALL sl_check( nf90_put_var(ncid, p_perp_id,
     $     sl_in%P_perp_arr))
      CALL sl_check( nf90_put_var(ncid, p_tor_id,
     $     sl_in%P_tor_arr) )
      CALL sl_check( nf90_put_var(ncid, Dnorm_id,
     $     sl_in%D_norm_arr))

c-----------------------------------------------------------------------
c     write complex Delta' as a RESHAPE'd (msing, 2) real array.
c-----------------------------------------------------------------------
      CALL sl_check( nf90_put_var(ncid, dpp_id,
     $     RESHAPE( (/sl_in%Re_dp_arr, sl_in%Im_dp_arr/),
     $              (/msing, 2/) )) )
      CALL sl_check( nf90_put_var(ncid, dc_id,
     $     sl_in%d_crit_arr) )

c-----------------------------------------------------------------------
c     write estimated growth-rate outputs (only when requested).
c-----------------------------------------------------------------------
      IF (est_gamma_flag) THEN
         CALL sl_check( nf90_put_var(ncid, dels_db_id,
     $        RESHAPE( (/REAL(sl_out%dels_db_arr),
     $                    AIMAG(sl_out%dels_db_arr)/),
     $                 (/msing, 2/) )) )

         CALL sl_check( nf90_put_var(ncid, d_b_id,
     $        sl_in%d_beta_arr) )

         CALL sl_check( nf90_put_var(ncid, ge_id,
     $        RESHAPE( (/REAL(sl_out%gamma_est_arr),
     $                    AIMAG(sl_out%gamma_est_arr)/),
     $                 (/msing, 2/) )) )
      END IF

c-----------------------------------------------------------------------
c     write solved growth rate (always present).
c-----------------------------------------------------------------------
      CALL sl_check( nf90_put_var(ncid, gs_id,
     $     RESHAPE( (/REAL(sl_out%gamma_sol_arr),
     $                 AIMAG(sl_out%gamma_sol_arr)/),
     $              (/msing, 2/) )) )

      IF (ALLOCATED(sl_out%br_th_arr)) THEN
         CALL sl_check( nf90_put_var(ncid, br_th_id,
     $        sl_out%br_th_arr) )
      END IF

c-----------------------------------------------------------------------
c     write AMR scan arrays.
c-----------------------------------------------------------------------
      CALL sl_check( nf90_put_var(ncid, var_npts_id, n_pts_arr) )
      CALL sl_check( nf90_put_var(ncid, var_q_id,    buffer_q)  )
      CALL sl_check( nf90_put_var(ncid, var_d_id,    buffer_d)  )

c-----------------------------------------------------------------------
c     deallocate AMR rectangular buffers.
c-----------------------------------------------------------------------
      DEALLOCATE(buffer_q, buffer_d, n_pts_arr)

c-----------------------------------------------------------------------
c     close the NetCDF file.
c-----------------------------------------------------------------------
      IF (debug_flag) PRINT *, " - Closing netcdf file"
      CALL sl_check( nf90_close(ncid) )

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slayer_netcdf_out
      END MODULE slayer_netcdf_mod