c-----------------------------------------------------------------------
c     file rdcon_netcdf.f
c     writes dcon.out information to a netcdf file
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. rdcon_netcdf_mod
c     1. check
c     2. rdcon_netcdf_out
c-----------------------------------------------------------------------
c     subprogram 0. rdcon_netcdf_mod
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rdcon_netcdf_mod
      USE rdcon_mod
      ! USE sing_mod, ONLY: sing_detf
      USE netcdf
      USE local_mod, ONLY: r4
      IMPLICIT NONE

      CONTAINS
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
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c-----------------------------------------------------------------------
c     subprogram 2. rdcon_netcdf_out.
c     Replicate dcon.out information in netcdf format.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE rdcon_netcdf_out(wp,wv,wt,wt0,ep,ev,et)

      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: ep,ev,et
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: wp,wv,wt,wt0

      INTEGER :: i, ncid,
     $    i_dim, m_dim, mo_dim, p_dim, i_id, m_id, mo_id, p_id,
     $    f_id, q_id, dv_id, mu_id, di_id, dr_id, ca_id,
     $    wp_id, wpv_id, wv_id, wvv_id, wt_id, wtv_id, wt0_id,
     $    l_dim, l_id, coil_dim, coil_id, dpc_id, dc_id,
     $    lp_dim, lp_id, r_dim, r_id, rp_dim, rp_id, pr_id, qr_id,
     $    dp_id, ap_id, bp_id, gp_id, dpp_id, lrc_dim, lrc_id

      REAL(r4) :: cpusec, wallsec
      CHARACTER(2) :: sn
      CHARACTER(64) :: ncfile

      INTEGER :: ising,jsing
      COMPLEX(r8), DIMENSION(msing,msing) :: ap,bp,gammap,deltap

      LOGICAL, PARAMETER :: debug_flag = .FALSE.
c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *,"Called rdcon_netcdf_out"
      IF (nn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)')nn
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')nn
      ENDIF
      ncfile = "rdcon_output_n"//TRIM(sn)//".nc"
      CALL timer(1,0,cpusec,wallsec)
c-----------------------------------------------------------------------
c     open files
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Creating netcdf files"
      CALL check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
c-----------------------------------------------------------------------
c     define global file attributes
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining modal netcdf globals"
      CALL check( nf90_put_att(ncid,nf90_global,"title",
     $     "RDCON outputs"))
      ! define global attributes
      CALL check( nf90_put_att(ncid,nf90_global,"jacobian",
     $     trim(jac_type)))
      CALL check( nf90_put_att(ncid,nf90_global,"power_bp", power_bp))
      CALL check( nf90_put_att(ncid,nf90_global,"power_b", power_b))
      CALL check( nf90_put_att(ncid,nf90_global,"power_r", power_r))
      CALL check( nf90_put_att(ncid,nf90_global,'mpsi', mpsi))
      CALL check( nf90_put_att(ncid,nf90_global,'mtheta', mtheta))
      CALL check( nf90_put_att(ncid,nf90_global,'mlow', mlow))
      CALL check( nf90_put_att(ncid,nf90_global,'mhigh', mhigh))
      CALL check( nf90_put_att(ncid,nf90_global,'mpert', mpert))
      CALL check( nf90_put_att(ncid,nf90_global,'mband', mband))
      CALL check( nf90_put_att(ncid,nf90_global,'psilow', psilow))
      CALL check( nf90_put_att(ncid,nf90_global,'amean', amean))
      CALL check( nf90_put_att(ncid,nf90_global,'rmean', rmean))
      CALL check( nf90_put_att(ncid,nf90_global,'aratio', aratio))
      CALL check( nf90_put_att(ncid,nf90_global,'kappa', kappa))
      CALL check( nf90_put_att(ncid,nf90_global,'delta1', delta1))
      CALL check( nf90_put_att(ncid,nf90_global,'delta2', delta2))
      CALL check( nf90_put_att(ncid,nf90_global,'li1', li1))
      CALL check( nf90_put_att(ncid,nf90_global,'li2', li2))
      CALL check( nf90_put_att(ncid,nf90_global,'li3', li3))
      CALL check( nf90_put_att(ncid,nf90_global,'ro', ro))
      CALL check( nf90_put_att(ncid,nf90_global,'zo', zo))
      CALL check( nf90_put_att(ncid,nf90_global,'psio', psio))
      CALL check( nf90_put_att(ncid,nf90_global,'betap1', betap1))
      CALL check( nf90_put_att(ncid,nf90_global,'betap2', betap2))
      CALL check( nf90_put_att(ncid,nf90_global,'betap3', betap3))
      CALL check( nf90_put_att(ncid,nf90_global,'betat', betat))
      CALL check( nf90_put_att(ncid,nf90_global,'betan', betan))
      CALL check( nf90_put_att(ncid,nf90_global,'bt0', bt0))
      CALL check( nf90_put_att(ncid,nf90_global,'q0', q0))
      CALL check( nf90_put_att(ncid,nf90_global,'qmin', qmin))
      CALL check( nf90_put_att(ncid,nf90_global,'qmax', qmax))
      CALL check( nf90_put_att(ncid,nf90_global,'qa', qa))
      CALL check( nf90_put_att(ncid,nf90_global,'crnt', crnt))
      CALL check( nf90_put_att(ncid,nf90_global,'q95', q95))
      CALL check( nf90_put_att(ncid,nf90_global,'betan', betan))
      CALL check( nf90_put_att(ncid,nf90_global,"qlim", qlim))
      CALL check( nf90_put_att(ncid,nf90_global,"psilim", psilim))
      CALL check( nf90_put_att(ncid,nf90_global,"shot", INT(shotnum)) )
      CALL check( nf90_put_att(ncid,nf90_global,"time",INT(shottime)) )
      CALL check( nf90_put_att(ncid,nf90_global,"n", nn))
      CALL check( nf90_put_att(ncid,nf90_global,"version", version))
      CALL check( nf90_put_att(ncid,nf90_global,"cpu_time",cpusec) )
      CALL check( nf90_put_att(ncid,nf90_global,"wall_time",wallsec))
      ! define dimensions
      CALL check( nf90_def_dim(ncid, "i", 2, i_dim) )
      CALL check( nf90_def_var(ncid, "i", nf90_int, i_dim, i_id) )
      CALL check( nf90_def_dim(ncid, "m", mpert,  m_dim) )
      CALL check( nf90_def_var(ncid, "m", nf90_int, m_dim, m_id) )
      CALL check( nf90_def_dim(ncid, "mode", mpert, mo_dim) )
      CALL check( nf90_def_var(ncid, "mode", nf90_int, mo_dim, mo_id))
      CALL check( nf90_def_dim(ncid, "psi_n", sq%mx+1, p_dim) )
      CALL check( nf90_def_var(ncid, "psi_n", nf90_double, p_dim, p_id))
      CALL check( nf90_put_att(ncid,p_id, "long_name",
     $       "Normalized Poloidal Flux") )

      IF(msing>0)THEN
         CALL check( nf90_def_dim(ncid,"lr_index",2*msing,l_dim) )
         CALL check( nf90_def_var(ncid,"lr_index",nf90_int,l_dim,l_id))
         CALL check( nf90_put_att(ncid,l_dim,"long_name",
     $       "Left and Right of Rational Surfaces Index") )
         CALL check( nf90_def_dim(ncid,"lr_prime",2*msing,lp_dim) )
         CALL check( nf90_def_var(ncid,"lr_prime",nf90_int,lp_dim,
     $                           lp_id) )
         CALL check( nf90_put_att(ncid,lp_id,"long_name",
     $       "Left and Right of Rational Surfaces Prime Index") )
         CALL check( nf90_def_dim(ncid,"r",msing,r_dim) )
         CALL check( nf90_def_var(ncid,"r",nf90_double,r_dim,r_id) )
         CALL check( nf90_put_att(ncid,r_id,"long_name",
     $       "Rational Surface Index") )
         CALL check( nf90_def_dim(ncid,"r_prime",msing,rp_dim) )
         CALL check( nf90_def_var(ncid,"r_prime",nf90_double,rp_dim,
     $                           rp_id) )
         CALL check( nf90_put_att(ncid,rp_id,"long_name",
     $       "Rational Surface Prime Index") )
         CALL check( nf90_def_var(ncid,"psi_n_rational",nf90_double,
     $                 r_dim, pr_id) )
         CALL check( nf90_put_att(ncid,pr_id,"long_name",
     $       "Normalized Poloidal Flux at Rational Surfaces") )    
         CALL check( nf90_def_var(ncid,"q_rational",nf90_double,
     $                 r_dim, qr_id) )
         CALL check( nf90_put_att(ncid,qr_id,"long_name",
     $       "Safety Factor at Rational Surfaces") )
         IF(coil%rpec_flag)THEN
            CALL check( nf90_def_dim(ncid,"lrc_index",
     $                              2*msing+coil%mcoil,lrc_dim) )
            CALL check( nf90_def_var(ncid,"lrc_index",nf90_int,lrc_dim,
     $                              lrc_id) )
            CALL check( nf90_put_att(ncid,lrc_id,"long_name",
     $          "Left and Right of Rational Surfaces and Coil Index") )
            CALL check( nf90_def_dim(ncid,"coil_index",coil%mcoil,
     $                              coil_dim) )
            CALL check( nf90_def_var(ncid,"coil_index",nf90_int,
     $                              coil_dim,coil_id) )
            CALL check( nf90_put_att(ncid,coil_id,"long_name",
     $          "Coil Index") )
         ENDIF

      ENDIF
      CALL check( nf90_def_var(ncid, "f", nf90_double, p_dim, f_id) )
      CALL check( nf90_def_var(ncid, "mu0p", nf90_double, p_dim, mu_id))
      CALL check( nf90_def_var(ncid, "dvdpsi", nf90_double,p_dim,dv_id))
      CALL check( nf90_def_var(ncid, "q", nf90_double, p_dim, q_id) )
         CALL check( nf90_put_att(ncid,q_id,"long_name",
     $       "Safety Factor") )
      CALL check( nf90_def_var(ncid, "di", nf90_double, p_dim, di_id) )
      CALL check( nf90_def_var(ncid, "dr", nf90_double, p_dim, dr_id) )
      CALL check( nf90_def_var(ncid, "ca1", nf90_double, p_dim, ca_id))
      CALL check( nf90_def_var(ncid, "W_p_eigenvector", nf90_double,
     $    (/m_dim, mo_dim, i_dim/), wp_id) )
      CALL check( nf90_def_var(ncid, "W_p_eigenvalue", nf90_double,
     $    (/mo_dim, i_dim/), wpv_id) )
      CALL check( nf90_put_att(ncid,wp_id,"long_name",
     $    "Plasma Energy Eigenmodes") )
      CALL check( nf90_put_att(ncid,wpv_id,"long_name",
     $    "Plasma Energy Eigenvalues") )
      CALL check( nf90_def_var(ncid, "W_v_eigenvector", nf90_double,
     $    (/m_dim, mo_dim, i_dim/), wv_id) )
      CALL check( nf90_def_var(ncid, "W_v_eigenvalue", nf90_double,
     $    (/mo_dim, i_dim/), wvv_id) )
      CALL check( nf90_put_att(ncid,wv_id,"long_name",
     $    "Vacuum Energy Eigenmodes") )
      CALL check( nf90_put_att(ncid,wvv_id,"long_name",
     $    "Vacuum Energy Eigenvalues") )
      CALL check( nf90_def_var(ncid, "W_t_eigenvector", nf90_double,
     $    (/m_dim, mo_dim, i_dim/), wt_id) )
      CALL check( nf90_def_var(ncid, "W_t_eigenvalue", nf90_double,
     $    (/mo_dim, i_dim/), wtv_id) )
      CALL check( nf90_put_att(ncid,wt_id,"long_name",
     $    "Total Energy Eigenmodes") )
      CALL check( nf90_put_att(ncid,wtv_id,"long_name",
     $    "Total Energy Eigenvalues") )
      CALL check( nf90_def_var(ncid, "W_t", nf90_double,
     $    (/m_dim, mo_dim, i_dim/), wt0_id) )
      CALL check( nf90_put_att(ncid,wt0_id,"long_name",
     $    "Total Energy Matrix") )
      IF(msing>0 .AND. ALLOCATED(delta))THEN
         CALL check( nf90_def_var(ncid, "Delta", nf90_double,
     $       (/l_dim, lp_dim, i_dim/), dp_id) )
         CALL check( nf90_put_att(ncid,dp_id,"long_name",
     $       "Galerkin Solution Delta Matrix") )
         CALL check( nf90_def_var(ncid, "A_prime", nf90_double,
     $       (/r_dim, rp_dim, i_dim/), ap_id) )
         CALL check( nf90_def_var(ncid, "B_prime", nf90_double,
     $       (/r_dim, rp_dim, i_dim/), bp_id) )
         CALL check( nf90_def_var(ncid, "Gamma_prime", nf90_double,
     $       (/r_dim, rp_dim, i_dim/), gp_id) )
         CALL check( nf90_def_var(ncid, "Delta_prime", nf90_double,
     $       (/r_dim, rp_dim, i_dim/), dpp_id) )
         CALL check( nf90_put_att(ncid,dpp_id,"long_name",
     $       "PEST3 Delta Prime Matrix" ))
         CALL check( nf90_def_var(ncid, "Delta_gw", nf90_double,
     $       (/lrc_dim, lp_dim, i_dim/), dc_id) )
         CALL check( nf90_put_att(ncid,dc_id,"long_name",
     $ "Galerkin Solution Delta Matrix with Coil Matching Terms"))
         IF(coil%rpec_flag)THEN
            CALL check( nf90_def_var(ncid, "Delta_coil", nf90_double,
     $         (/coil_dim, lp_dim, i_dim/), dpc_id) )
            CALL check( nf90_put_att(ncid,dpc_id,"long_name",
     $      "Galerkin Solution Delta Matrix Coil Matching Terms Only"))
         ENDIF
      ENDIF

      ! end definitions
      CALL check( nf90_enddef(ncid) )
c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Putting dimension variables in netcdf"
      CALL check( nf90_put_var(ncid,i_id, (/0,1/)) )
      CALL check( nf90_put_var(ncid,m_id, (/(i+mlow, i=0,mpert-1)/)) )
      CALL check( nf90_put_var(ncid,mo_id, (/(i, i=1,mpert)/)) )
      CALL check( nf90_put_var(ncid,p_id, sq%xs(:)))
      IF(msing>0)THEN
         CALL check( nf90_put_var(ncid,l_id, (/(i,i=1,2*msing)/)) )
         CALL check( nf90_put_var(ncid,lp_id, (/(i,i=1,2*msing)/)) )
         CALL check( nf90_put_var(ncid, r_id, (/(sing(i)%m,
     $                                           i=1,msing)/)) )
         CALL check( nf90_put_var(ncid,rp_id, (/(sing(i)%m,
     $                                           i=1,msing)/)) )
         CALL check( nf90_put_var(ncid,pr_id, (/(sing(i)%psifac,
     $                                           i=1,msing)/)) )
         CALL check( nf90_put_var(ncid,qr_id, (/(sing(i)%q,
     $                                           i=1,msing)/)) )
         CALL check( nf90_put_var(ncid,lrc_id,
     $       (/(i,i=1,2*msing+coil%mcoil)/)))
         CALL check( nf90_put_var(ncid,coil_id,
     $       (/(i,i=1,coil%mcoil)/)) )
      ENDIF

      IF(debug_flag) PRINT *," - Putting profile variables in netcdf"
      CALL check( nf90_put_var(ncid,f_id, sq%fs(:,1)/twopi))
      CALL check( nf90_put_var(ncid,mu_id, sq%fs(:,2)))
      CALL check( nf90_put_var(ncid,dv_id, sq%fs(:,3)))
      CALL check( nf90_put_var(ncid,q_id, sq%fs(:,4)))
      CALL check( nf90_put_var(ncid,di_id, locstab%fs(:,1)/sq%xs(:)))
      CALL check( nf90_put_var(ncid,dr_id, locstab%fs(:,2)/sq%xs(:)))
      CALL check( nf90_put_var(ncid,ca_id, locstab%fs(:,4)))

      IF(debug_flag) PRINT *," - Putting matrix variables in netcdf"
      CALL check( nf90_put_var(ncid,wp_id,RESHAPE((/REAL(wp),
     $             AIMAG(wp)/),(/mpert,mpert,2/))) )
      CALL check( nf90_put_var(ncid,wpv_id,RESHAPE((/REAL(ep),
     $             AIMAG(ep)/),(/mpert,2/))) )
      CALL check( nf90_put_var(ncid,wv_id,RESHAPE((/REAL(wv),
     $             AIMAG(wv)/),(/mpert,mpert,2/))) )
      CALL check( nf90_put_var(ncid,wvv_id,RESHAPE((/REAL(ev),
     $             AIMAG(ev)/),(/mpert,2/))) )
      CALL check( nf90_put_var(ncid,wt_id,RESHAPE((/REAL(wt),
     $             AIMAG(wt)/),(/mpert,mpert,2/))) )
      CALL check( nf90_put_var(ncid,wtv_id,RESHAPE((/REAL(et),
     $             AIMAG(et)/),(/mpert,2/))) )
      CALL check( nf90_put_var(ncid,wt0_id,RESHAPE((/REAL(wt0),
     $             AIMAG(wt0)/),(/mpert,mpert,2/))) )


      IF(msing>0 .AND. ALLOCATED(delta))THEN
         ! construct PEST3 matching data
         ap=0.0_r8
         bp=0.0_r8
         gammap=0.0_r8
         deltap=0.0_r8
         DO ising=1,msing
            DO jsing=1,msing
               ap(ising,jsing) = delta(2*ising,2*jsing)
     $            +delta(2*ising,2*jsing-1)
     $            +delta(2*ising-1,2*jsing)
     $            +delta(2*ising-1,2*jsing-1)
               bp(ising,jsing) = delta(2*ising,2*jsing)
     $            -delta(2*ising,2*jsing-1)
     $            +delta(2*ising-1,2*jsing)
     $            -delta(2*ising-1,2*jsing-1)
               gammap(ising,jsing) = delta(2*ising,2*jsing)
     $            +delta(2*ising,2*jsing-1)
     $            -delta(2*ising-1,2*jsing)
     $            -delta(2*ising-1,2*jsing-1)
               deltap(ising,jsing) = delta(2*ising,2*jsing)
     $            -delta(2*ising,2*jsing-1)
     $            -delta(2*ising-1,2*jsing)
     $            +delta(2*ising-1,2*jsing-1)
            ENDDO
         ENDDO
         CALL check( nf90_put_var(ncid,dp_id,
     $    RESHAPE((/REAL(delta(1:2*msing,:)),IMAG(delta(1:2*msing,:))/),
     $    (/2*msing,2*msing,2/))) )
         CALL check( nf90_put_var(ncid,ap_id,
     $    RESHAPE((/REAL(ap),IMAG(ap)/),(/msing,msing,2/)) ) )
         CALL check( nf90_put_var(ncid,bp_id,
     $    RESHAPE((/REAL(bp),IMAG(bp)/),(/msing,msing,2/)) ) )
         CALL check( nf90_put_var(ncid,gp_id,
     $    RESHAPE((/REAL(gammap),IMAG(gammap)/),(/msing,msing,2/)) ) )
         CALL check( nf90_put_var(ncid,dpp_id,
     $    RESHAPE((/REAL(deltap),IMAG(deltap)/),(/msing,msing,2/)) ) )

         IF(coil%rpec_flag)THEN
            CALL check( nf90_put_var(ncid,dc_id,
     $       RESHAPE((/REAL(delta),IMAG(delta)/),(/2*msing+coil%mcoil,
     $                                               2*msing,2/)) ) )
            CALL check( nf90_put_var(ncid,dpc_id,
     $       RESHAPE((/REAL(delta(2*msing+1:2*msing+coil%mcoil,:)),
     $        IMAG(delta(2*msing+1:2*msing+coil%mcoil,:))/),
     $       (/coil%mcoil,2*msing,2/)) ) )


         ENDIF
      ENDIF   
c-----------------------------------------------------------------------
c     close file
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Closing netcdf file"
      CALL check( nf90_close(ncid) )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rdcon_netcdf_out

      END MODULE rdcon_netcdf_mod