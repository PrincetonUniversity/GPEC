c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     GPEC: main program
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     gpec_main.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM gpec_main
      USE gpdiag_mod
      USE gpout_mod
      USE rdcon_mod

      IMPLICIT NONE

      INTEGER :: i,j,in,resol,
     $     mode,m3mode,lowmode,highmode,filter_modes,
     $     sing_npsi
      INTEGER :: pmode,p1mode,rmode,dmode,d1mode,fmode,smode,tmp_outs(5)
      INTEGER, DIMENSION(:), POINTER :: ipiv
      REAL(r8) :: sing_spot,majr,minr,smallwidth,fp,normpsi
      CHARACTER(8) :: filter_types
      CHARACTER(128) :: infile
      LOGICAL :: singcoup_flag,singfld_flag,vsingfld_flag,pmodb_flag,
     $     xbcontra_flag,xbnormal_flag,vbnormal_flag,xbnobo_flag,
     $     d3_flag,xbst_flag,rzphibx_flag,dw_flag,
     $     radvar_flag,eigen_flag,magpot_flag,xbtangent_flag,
     $     arbsurf_flag,angles_flag,surfmode_flag,rzpgrid_flag,
     $     singcurs_flag,m3d_flag,cas3d_flag,test_flag,nrzeq_flag,
     $     arzphifun_flag,xbrzphifun_flag,pmodbmn_flag,xclebsch_flag,
     $     filter_flag,gal_flag
      LOGICAL, DIMENSION(100) :: ss_flag
      COMPLEX(r8), DIMENSION(:), POINTER :: finmn,foutmn,xspmn,
     $     fxmn,fxfun,coilmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: invmats,temp1

      NAMELIST/gpec_input/dcon_dir,ieqfile,idconfile,ivacuumfile,
     $     power_flag,fft_flag,mthsurf0,fixed_boundary_flag,
     $     data_flag,data_type,nmin,nmax,mmin,mmax,jsurf_in,mthsurf,
     $     jac_in,power_bin,power_rin,power_bpin,power_rcin,tmag_in,
     $     infile,harmonic_flag,mode_flag,sinmn,cosmn,
     $     displacement_flag,mode,coil_flag,
     $     ip_direction,bt_direction,rdconfile,
     $     pmode,p1mode,dmode,d1mode,fmode,rmode,smode,
     $     filter_types,filter_modes,gal_flag
      NAMELIST/gpec_control/resp_index,resp_induct_flag,
     $     sing_spot,sing_npsi,reg_flag,reg_spot,
     $     chebyshev_flag,nche,nchr,nchz,use_classic_splines
      NAMELIST/gpec_output/resp_flag,singcoup_flag,nrzeq_flag,nr,nz,
     $     singfld_flag,pmodb_flag,xbnormal_flag,rstep,jsurf_out,
     $     jac_out,power_bout,power_rout,power_bpout,power_rcout,
     $     tmag_out,mlim_out,eqbrzphi_flag,brzphi_flag,xrzphi_flag,
     $     vbrzphi_flag,vvbrzphi_flag,divzero_flag,dw_flag,opsi1,opsi2,
     $     bin_flag,bin_2d_flag,fun_flag,flux_flag,bwp_pest_flag,
     $     vsbrzphi_flag,ss_flag,arzphifun_flag,xbrzphifun_flag,
     $     vsingfld_flag,vbnormal_flag,eigm_flag,xbtangent_flag,
     $     xclebsch_flag,pbrzphi_flag,verbose,max_linesout,filter_flag,
     $     netcdf_flag,ascii_flag
      NAMELIST/gpec_diagnose/singcurs_flag,xbcontra_flag,
     $     xbnobo_flag,d3_flag,div_flag,xbst_flag,
     $     pmodbmn_flag,rzphibx_flag,radvar_flag,eigen_flag,magpot_flag,
     $     arbsurf_flag,majr,minr,angles_flag,surfmode_flag,
     $     lowmode,highmode,rzpgrid_flag,m3d_flag,m3mode,
     $     cas3d_flag,test_flag,resol,smallwidth,debug_flag,timeit,
     $     malias
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      sinmn=0
      cosmn=0
      jsurf_in=0
      tmag_in=1
      jac_in=""
      dcon_dir=""
      ieqfile="psi_in.bin"
      idconfile="euler.bin"
      ivacuumfile="vacuum.bin"
      rdconfile="globalsol.bin"
      power_flag=.TRUE.
      fft_flag=.FALSE.
      fixed_boundary_flag=.FALSE.
      data_flag=.FALSE.
      harmonic_flag=.FALSE.
      mode_flag=.FALSE.
      displacement_flag=.FALSE.
      gal_flag=.FALSE.
      mthsurf0=1
      nmin=1
      nmax=1
      mmin=-128
      mmax=128
      mthsurf=0
      pmode =0
      p1mode =0
      rmode =0
      smode =0
      fmode =0
      dmode =0
      d1mode=0
      filter_modes = 0
      filter_types = '   '

      resp_index=0
      resp_induct_flag=.TRUE.
      sing_spot=5e-4
      sing_npsi=1e2
      reg_flag=.TRUE.
      reg_spot=5e-2
      chebyshev_flag=.TRUE.
      nche=20
      nchr=20
      nchz=20

      jsurf_out=0
      tmag_out=1
      mlim_out=64
      jac_out=""
      resp_flag=.TRUE.
      singcoup_flag=.FALSE.
      singfld_flag=.TRUE.
      vsingfld_flag=.FALSE.
      pmodb_flag=.FALSE.
      xbnormal_flag=.TRUE.
      vbnormal_flag=.FALSE.
      xbtangent_flag=.FALSE.
      xclebsch_flag=.FALSE.
      filter_flag  = .TRUE.
      rstep=0
      nrzeq_flag=.FALSE.
      nr=64
      nz=64
      eqbrzphi_flag=.FALSE.
      brzphi_flag=.FALSE.
      xrzphi_flag=.FALSE.
      vbrzphi_flag=.FALSE.
      pbrzphi_flag=.FALSE.      
      vvbrzphi_flag=.FALSE.
      bin_flag=.TRUE.
      bin_2d_flag=.TRUE.
      netcdf_flag=.TRUE.
      ascii_flag=.TRUE.
      fun_flag=.FALSE.
      flux_flag=.FALSE.
      max_linesout=0
      vsbrzphi_flag=.FALSE.
      DO i=1,100 
         ss_flag(i)=.FALSE.
      ENDDO
      arzphifun_flag=.FALSE.
      xbrzphifun_flag=.FALSE.
      bwp_pest_flag=.TRUE.
      dw_flag = .FALSE.

      singcurs_flag=.FALSE.
      xbcontra_flag=.FALSE.
      xbnobo_flag=.FALSE.
      d3_flag=.FALSE.
      div_flag=.FALSE.
      xbst_flag=.FALSE.
      pmodbmn_flag=.FALSE.
      rzphibx_flag=.FALSE.
      radvar_flag=.TRUE.
      eigen_flag=.FALSE.
      magpot_flag=.FALSE.
      arbsurf_flag=.FALSE.
      angles_flag=.FALSE.
      surfmode_flag=.FALSE.
      rzpgrid_flag=.FALSE.
      m3d_flag=.FALSE.
      cas3d_flag=.FALSE.
      test_flag=.FALSE.
      eigm_flag=.FALSE.

      majr=10.0
      minr=1.0
      mode=1
      lowmode=-10
      highmode=10
      m3mode=2
      resol=INT(1e4)
      smallwidth=1e-6
      
      timeit = .FALSE.
      verbose = .TRUE.
      debug_flag = .FALSE.
      malias=0

      psixy = 0
      opsi1=0.0
      opsi2=1.0
c-----------------------------------------------------------------------
c     read gpec.in.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)""
      IF(verbose) WRITE(*,*)"GPEC START => "//TRIM(version)
      IF(verbose) WRITE(*,*)"__________________________________________"
      CALL ascii_open(in_unit,"gpec.in","OLD")
      READ(in_unit,NML=gpec_input)
      READ(in_unit,NML=gpec_control)
      READ(in_unit,NML=gpec_output)
      READ(in_unit,NML=gpec_diagnose)
      CALL ascii_close(in_unit)
      galsol%gal_flag=gal_flag
      IF(timeit) CALL gpec_timer(0)
c-----------------------------------------------------------------------
c     Deprecated variable errors
c-----------------------------------------------------------------------
      IF((pmode/=0).or.(p1mode/=0).or.(dmode/=0).or.(d1mode/=0).or.
     $   (fmode/=0).or.(rmode/=0).or.(smode/=0))THEN
         PRINT *,"WARNING: p/d/f/r/smode syntax is a deprecated!"
         PRINT *,"  Use filter_types to filter external spectrum."
         CALL gpec_stop("Deprecated input.")
      ENDIF
      IF(malias/=0) THEN
       PRINT *,"WARNING: malias may not be supported in future versions"
      ENDIF
c-----------------------------------------------------------------------
c     define relative file paths.
c-----------------------------------------------------------------------
      IF(dcon_dir/="")THEN
         CALL setahgdir(dcon_dir)
         idconfile = TRIM(dcon_dir)//"/"//TRIM(idconfile)
         ieqfile = TRIM(dcon_dir)//"/"//TRIM(ieqfile)
         ivacuumfile = TRIM(dcon_dir)//"/"//TRIM(ivacuumfile)
         rdconfile = TRIM(dcon_dir)//"/"//TRIM(rdconfile)
      ENDIF
c-----------------------------------------------------------------------
c     define coordinates.
c-----------------------------------------------------------------------
      SELECT CASE(jac_in)
      CASE("hamada")
         power_bin=0
         power_bpin=0
         power_rin=0
         power_rcin=0
      CASE("pest")
         power_bin=0
         power_bpin=0
         power_rin=2
         power_rcin=0
      CASE("equal_arc")
         power_bin=0
         power_bpin=1
         power_rin=0
         power_rcin=0
      CASE("boozer")
         power_bin=2
         power_bpin=0
         power_rin=0
         power_rcin=0
      CASE("park")
         power_bin=1
         power_bpin=0
         power_rin=0
         power_rcin=0
      CASE("polar")
         power_bin=0
         power_bpin=1
         power_rin=0
         power_rcin=1         
      CASE("other")
      CASE DEFAULT
      END SELECT

      SELECT CASE(jac_out)
      CASE("hamada")
         power_bout=0
         power_bpout=0
         power_rout=0
         power_rcout=0
      CASE("pest")
         power_bout=0
         power_bpout=0
         power_rout=2
         power_rcout=0
      CASE("equal_arc")
         power_bout=0
         power_bpout=1
         power_rout=0
         power_rcout=0
      CASE("boozer")
         power_bout=2
         power_bpout=0
         power_rout=0
         power_rcout=0
      CASE("park")
         power_bout=1
         power_bpout=0
         power_rout=0
         power_rcout=0
      CASE("polar")
         power_bout=0
         power_bpout=1
         power_rout=0
         power_rcout=1         
      CASE("other")
      CASE DEFAULT
      END SELECT
c-----------------------------------------------------------------------
c     set parameters from inputs.
c-----------------------------------------------------------------------
      IF (brzphi_flag .OR. vbrzphi_flag .OR. eqbrzphi_flag) psixy=1
      IF (xrzphi_flag) psixy=1
      IF (rstep==0) rstep=mstep
      IF (nchr==20) nchr=nche
      IF (nchz==20) nchz=nche
      IF (max_linesout==0) max_linesout=-1
c-----------------------------------------------------------------------
c     check time.
c-----------------------------------------------------------------------
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
c-----------------------------------------------------------------------
c     prepare DCON solutions.
c-----------------------------------------------------------------------
      CALL idcon_read(psixy)
      CALL idcon_transform
c-----------------------------------------------------------------------
c     replace radial grids and u1-u4 with RDON solutions.
c-----------------------------------------------------------------------
      IF(gal_flag) CALL rdcon_read_solution
c-----------------------------------------------------------------------
c     reconstruct metric tensors.
c-----------------------------------------------------------------------
      CALL idcon_metric
c-----------------------------------------------------------------------
c     read vacuum data.
c-----------------------------------------------------------------------
      CALL idcon_vacuum
      IF(timeit) CALL gpec_timer(2)
c-----------------------------------------------------------------------
c     set parameters from dcon.
c-----------------------------------------------------------------------
      ALLOCATE(xspmn(mpert),finmn(mpert),foutmn(mpert),
     $     fxmn(mpert),fxfun(mthsurf))
      xspmn = 0
      finmn = 0
      foutmn = 0
      fxmn = 0
      fxfun = 0
c-----------------------------------------------------------------------
c     read coil data.
c-----------------------------------------------------------------------
      cmlow=mlow
      cmhigh=mhigh
      cmpert=mpert
      finmn=0
      ipd=1.0
      btd=1.0
      IF(ip_direction=="negative")ipd=-1.0
      IF(bt_direction=="negative")btd=-1.0
      helicity=ipd*btd
      IF (coil_flag) THEN
         IF(verbose) WRITE(*,*)
     $     "Calculating field on the boundary from coils"
         CALL coil_read(idconfile)
         ALLOCATE(coilmn(cmpert))
         ALLOCATE(coil_indmat(mpert,coil_num))
         DO j=1,coil_num
            coilmn=0
            CALL field_bs_psi(psilim,coilmn,1,op_start=j,op_stop=j)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  coil_indmat(cmlow-mlow+i,j)=coilmn(i)
                  finmn(cmlow-mlow+i)=finmn(cmlow-mlow+i)+coilmn(i)
               ENDIF
            ENDDO
         ENDDO
         DEALLOCATE(coilmn)
         IF(timeit) CALL gpec_timer(2)
      ENDIF
c-----------------------------------------------------------------------
c     log inputs with harvest
c-----------------------------------------------------------------------
      ierr=init_harvest('CODEDB_GPEC'//NUL,hlog,len(hlog))
      ierr=set_harvest_verbose(0)
      ! standard CODEDB records
      ierr=set_harvest_payload_str(hlog,'CODE'//nul,'GPEC'//nul)
      IF (machine=='') then
         machine = "UNKNOWN"
      ELSEIF (machine=='d3d') then
         machine = "DIII-D"
      ENDIF
      machine = to_upper(machine)
      ierr=set_harvest_payload_str(hlog,'MACHINE'//nul,
     $                             trim(machine)//nul)
      ierr=set_harvest_payload_str(hlog,'VERSION'//nul,version//nul)
      if(shotnum>0)
     $   ierr=set_harvest_payload_int(hlog,'SHOT'//nul,INT(shotnum))
      if(shottime>0)
     $   ierr=set_harvest_payload_int(hlog,'TIME'//nul,INT(shottime))
      ! DCON equilibrium descriptors
      ierr=set_harvest_payload_int(hlog,'mpsi'//nul,mpsi)
      ierr=set_harvest_payload_int(hlog,'mtheta'//nul,mtheta)
      ierr=set_harvest_payload_int(hlog,'mlow'//nul,mlow)
      ierr=set_harvest_payload_int(hlog,'mhigh'//nul,mhigh)
      ierr=set_harvest_payload_int(hlog,'mpert'//nul,mpert)
      ierr=set_harvest_payload_int(hlog,'mband'//nul,mband)
      ierr=set_harvest_payload_dbl(hlog,'psilow'//nul,psilow)
      ierr=set_harvest_payload_dbl(hlog,'psilim'//nul,psilim)
      ierr=set_harvest_payload_dbl(hlog,'amean'//nul,amean)
      ierr=set_harvest_payload_dbl(hlog,'rmean'//nul,rmean)
      ierr=set_harvest_payload_dbl(hlog,'aratio'//nul,aratio)
      ierr=set_harvest_payload_dbl(hlog,'kappa'//nul,kappa)
      ierr=set_harvest_payload_dbl(hlog,'delta1'//nul,delta1)
      ierr=set_harvest_payload_dbl(hlog,'delta2'//nul,delta2)
      ierr=set_harvest_payload_dbl(hlog,'li1'//nul,li1)
      ierr=set_harvest_payload_dbl(hlog,'li2'//nul,li2)
      ierr=set_harvest_payload_dbl(hlog,'li3'//nul,li3)
      ierr=set_harvest_payload_dbl(hlog,'ro'//nul,ro)
      ierr=set_harvest_payload_dbl(hlog,'zo'//nul,zo)
      ierr=set_harvest_payload_dbl(hlog,'psio'//nul,psio)
      ierr=set_harvest_payload_dbl(hlog,'betap1'//nul,betap1)
      ierr=set_harvest_payload_dbl(hlog,'betap2'//nul,betap2)
      ierr=set_harvest_payload_dbl(hlog,'betap3'//nul,betap3)
      ierr=set_harvest_payload_dbl(hlog,'betat'//nul,betat)
      ierr=set_harvest_payload_dbl(hlog,'betan'//nul,betan)
      ierr=set_harvest_payload_dbl(hlog,'bt0'//nul,bt0)
      ierr=set_harvest_payload_dbl(hlog,'q0'//nul,q0)
      ierr=set_harvest_payload_dbl(hlog,'qmin'//nul,qmin)
      ierr=set_harvest_payload_dbl(hlog,'qmax'//nul,qmax)
      ierr=set_harvest_payload_dbl(hlog,'qa'//nul,qa)
      ierr=set_harvest_payload_dbl(hlog,'crnt'//nul,crnt)
      ierr=set_harvest_payload_dbl(hlog,'q95'//nul,q95)
      ierr=set_harvest_payload_dbl(hlog,'betan'//nul,betan)
      ierr=set_harvest_payload_dbl_array(hlog,'et'//nul,et,mpert)
      ierr=set_harvest_payload_dbl_array(hlog,'ep'//nul,ep,mpert)
      ! gpec_input
      ierr=set_harvest_payload_bol(hlog,'fixed_boundary_flag'//nul,
     $                             fixed_boundary_flag)
      ierr=set_harvest_payload_bol(hlog,'mode_flag'//nul,mode_flag)
      ierr=set_harvest_payload_int(hlog,'mode'//nul,mode)
      ierr=set_harvest_payload_str(hlog,'filter_types'//nul,
     $                             filter_types//nul)
      ierr=set_harvest_payload_int(hlog,'filter_modes'//nul,
     $                             filter_modes)
      ! gpec_control
      ierr=set_harvest_payload_int(hlog,'resp_index'//nul,resp_index)
      ierr=set_harvest_payload_dbl(hlog,'sing_spot'//nul,sing_spot)
      ierr=set_harvest_payload_dbl(hlog,'sing_npsi'//nul,sing_npsi)
      ierr=set_harvest_payload_bol(hlog,'reg_flag'//nul,reg_flag)
      ierr=set_harvest_payload_dbl(hlog,'reg_spot'//nul,reg_spot)
      ! gpec_output
      ierr=set_harvest_payload_str(hlog,'jac_out'//nul,jac_out//nul)
      ierr=set_harvest_payload_int(hlog,'jsurf_out'//nul,jsurf_out)
      ierr=set_harvest_payload_int(hlog,'tmag_out'//nul,tmag_out)
c-----------------------------------------------------------------------
c     compute plasma response.
c-----------------------------------------------------------------------
      CALL gpresp_eigen
      IF(timeit) CALL gpec_timer(2)
      CALL gpresp_pinduct
      CALL gpresp_sinduct
      CALL gpresp_permeab
      CALL gpresp_reluct
      IF(timeit) CALL gpec_timer(2)
c-----------------------------------------------------------------------
c     Set parameters for outputs.
c-----------------------------------------------------------------------
      IF (nrzeq_flag) THEN
         nr=mr
         nz=mz
      ENDIF
      CALL gpeq_rzpgrid(nr,nz,psixy)
c-----------------------------------------------------------------------
c     full analysis.
c-----------------------------------------------------------------------
      IF (netcdf_flag) CALL gpout_init_netcdf
      IF (resp_flag) THEN
         CALL gpout_response(power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out,jsurf_out)
      ENDIF
      DO i=1,LEN_TRIM(filter_types)
         IF(filter_types(i:i)=='s') singcoup_flag=.TRUE.
      ENDDO
      IF (singcoup_flag) THEN
         IF (msing==0) THEN
            PRINT *,"WARNING: no rationals for singcoup_flag"
            singcoup_flag = .FALSE.
         ELSE         
            CALL gpout_singcoup(sing_spot,sing_npsi,power_rout,
     $           power_bpout,power_bout,power_rcout,tmag_out)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and external field.
c-----------------------------------------------------------------------
      IF (data_flag.OR.harmonic_flag.OR.coil_flag)THEN
         edge_flag=.TRUE.
         CALL gpout_control(infile,finmn,foutmn,xspmn,power_rin,
     $        power_bpin,power_bin,power_rcin,tmag_in,jsurf_in,
     $        power_rout,power_bpout,power_bout,power_rcout,tmag_out,
     $        jsurf_out,filter_types,filter_modes,filter_flag)
      ELSE IF (mode_flag) THEN
         edge_flag=.FALSE.
      ENDIF

      IF (singfld_flag) THEN
         IF (con_flag) THEN
            PRINT *,"WARNING: singfld_flag not supported with con_flag"
            singfld_flag = .FALSE.
            vsingfld_flag = .FALSE.
         ELSEIF (msing==0) THEN
            PRINT *,"WARNING: no rationals for singfld_flag"
            singfld_flag = .FALSE.
            vsingfld_flag = .FALSE.
         ELSE
            CALL gpout_singfld(mode,xspmn,sing_spot,sing_npsi)
         ENDIF
      ENDIF
      IF (coil_flag .AND. vsingfld_flag) THEN
         CALL gpout_vsingfld()
      ENDIF
      ! here we see the subroutine is simply called in series with other
      ! similar subroutines by the driving program here
      ! this and gpec_pmodb, for example, are completely independent
      ! (just different ways of breaking up the components of the fields/displacements)
      IF(netcdf_flag.and.(xclebsch_flag.or.dw_flag.or.pmodb_flag
     $   .or.xbnormal_flag.or.xbtangent_flag.or.vbnormal_flag))THEN
         CALL gpout_qrv
      ENDIF
      IF (xclebsch_flag) THEN
         CALL gpout_xclebsch(mode,xspmn)
      ENDIF
      IF (kin_flag .AND. dw_flag) THEN
         CALL gpout_dw(mode,xspmn)
         CALL gpout_dw_matrix(coil_flag)
      ENDIF
      IF (pmodb_flag) THEN
         CALL gpout_pmodb(mode,xspmn)
      ENDIF
      IF (xbnormal_flag) THEN
         CALL gpout_xbnormal(mode,xspmn,sing_spot,sing_npsi)
      ENDIF
      IF (xbtangent_flag) THEN
         CALL gpout_xbtangent(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (coil_flag .AND. vbnormal_flag) THEN
         CALL gpout_vbnormal(power_rout,power_bpout,power_bout,
     $        power_rcout,tmag_out)
      ENDIF
      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag .OR. vvbrzphi_flag .OR. pbrzphi_flag) THEN
         IF (.NOT.mode_flag) THEN
            CALL gpout_xbrzphi(mode,xspmn,nr,nz,finmn,foutmn)
         ELSE
            ALLOCATE(ipiv(mpert),
     $           invmats(mpert,mpert),temp1(mpert,mpert))
            DO i=1,mpert
               invmats(i,i)=1.0
            ENDDO
            temp1=TRANSPOSE(permeabmats(resp_index,:,:))
            CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
            CALL zgetrs('N',mpert,mpert,temp1,mpert,
     $           ipiv,invmats,mpert,info)
            invmats=TRANSPOSE(invmats)
            CALL idcon_build(mode,xspmn)
            CALL gpeq_alloc
            CALL gpeq_sol(psilim)
            CALL gpeq_contra(psilim)
            CALL gpeq_dealloc
            CALL gpeq_weight(psilim,foutmn,mfac,mpert,1)
            finmn = MATMUL(invmats,foutmn)
            CALL gpeq_weight(psilim,foutmn,mfac,mpert,0)
            CALL gpeq_weight(psilim,finmn,mfac,mpert,0)
            CALL gpout_xbrzphi(mode,xspmn,nr,nz,finmn,foutmn)
            DEALLOCATE(ipiv,invmats,temp1)
         ENDIF
      ENDIF
      IF (singfld_flag .AND. vsbrzphi_flag) THEN
         DO i=1,msing
            IF (ss_flag(i)) CALL gpout_vsbrzphi(i,nr,nz)
         ENDDO
      ENDIF

      IF (xbrzphifun_flag) THEN
         CALL gpout_xbrzphifun(mode,xspmn)
      ENDIF
      IF (arzphifun_flag) THEN
         CALL gpout_arzphifun(mode,xspmn)
      ENDIF

c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF (singcurs_flag) THEN
         CALL gpdiag_singcurs(mode,xspmn,msing,resol,smallwidth)
      ENDIF
      IF (xbcontra_flag) THEN
         CALL gpdiag_xbcontra(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (xbnobo_flag) THEN
         CALL gpdiag_xbnobo(mode,xspmn,d3_flag)
      ENDIF
      IF (xbst_flag) THEN
         CALL gpdiag_xbst(mode,xspmn)
      ENDIF
      IF (pmodbmn_flag) THEN
         CALL gpdiag_pmodb(mode,xspmn)
         CALL gpdiag_pmodbmn(mode,xspmn)
      ENDIF            
      IF (rzphibx_flag) THEN
         CALL gpdiag_rzphibx(mode,xspmn)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose without a given error field.
c-----------------------------------------------------------------------
      IF (radvar_flag) THEN
         CALL gpdiag_radvar
      ENDIF
      IF (eigen_flag) THEN
         CALL gpdiag_eigen
      ENDIF
      IF (magpot_flag) THEN
         CALL gpdiag_magpot
      ENDIF
      IF (arbsurf_flag) THEN
         CALL gpdiag_arbsurf(majr,minr)
      ENDIF         
      IF (angles_flag) THEN
         CALL gpdiag_angles
      ENDIF

      IF (surfmode_flag) THEN
         CALL gpdiag_surfmode(lowmode,highmode,power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out,jsurf_out)
      ENDIF
      IF (rzpgrid_flag) THEN
         CALL gpdiag_rzpgrid(nr,nz)
      ENDIF
      
      IF (m3d_flag) THEN
         normpsi=1.0
         fp=1e-3         
         fxmn=0
         fxmn(m3mode-mlow+1)=fp*normpsi
         CALL gpeq_fcoords(psilim,fxmn,mfac,mpert,0,1,0,1,0,0)
         fxmn=-twopi*ifac*chi1*(mfac-nn*qlim)*fxmn
         CALL gpeq_weight(psilim,fxmn,mfac,mpert,0)
         CALL gpout_control(infile,fxmn,foutmn,xspmn,
     $        0,0,0,0,1,0,0,0,0,0,1,0,'   ',0,.FALSE.)
         edge_flag=.TRUE.
         CALL gpout_singfld(mode,xspmn,sing_spot,sing_npsi)
      ENDIF

      IF (cas3d_flag) THEN
         fp = -1e-2
         
         fxmn=0
         fxmn(m3mode-mlow+1)=fp
         ! temporary override of output options
         tmp_outs = (/power_rout,power_bpout,power_bout,power_rcout,
     $                tmag_out/)
         power_rout = 0
         power_bpout = 0
         power_bout = 2
         power_rcout = 0
         tmag_out = 1
         CALL gpeq_fcoords(psilim,fxmn,mfac,mpert,0,0,2,0,1,0)
         CALL gpout_control(infile,finmn,foutmn,xspmn,power_rin,
     $        power_bpin,power_bin,power_rcin,tmag_in,jsurf_in,
     $        power_rout,power_bpout,power_bout,power_rcout,
     $        tmag_out,jsurf_out,'   ',0,.FALSE.)
         edge_flag=.TRUE.
         CALL gpout_singfld(mode,xspmn,sing_spot,sing_npsi)
         CALL gpdiag_xbcontra(mode,xspmn,0,0,2,0,1)
         CALL gpout_xbnormal(mode,xspmn,sing_spot,sing_npsi)
         CALL gpdiag_xbnobo(mode,xspmn,d3_flag)
         CALL gpdiag_radvar
         ! reset output options
         power_rout = tmp_outs(1)
         power_bpout = tmp_outs(2)
         power_bout = tmp_outs(3)
         power_rcout = tmp_outs(4)
         tmag_out = tmp_outs(5)
      ENDIF
c-----------------------------------------------------------------------
c     various simple test.
c-----------------------------------------------------------------------
      IF (test_flag) THEN
         fxmn=0
         fxmn(3-mlow+1)=1.0
         fxmn(7-mlow+1)=-3.0*ifac
         CALL ascii_open(out_unit,"gptest_coordtrans_n3.out","UNKNOWN")
         WRITE(out_unit,*)
     $        "GPTEST_COORDTRANS: hamada to pest and to hamada back"
         WRITE(out_unit,'(6(1x,a16))')"normpsi","mfac","real1",
     $        "real2","imag1","imag2"
         DO i=1,99
            WRITE(*,*)i
            foutmn=fxmn 
            normpsi = REAL(i)/100.0
            CALL gpeq_bcoords(normpsi,foutmn,mfac,mpert,2,0,0,0,0,0)
            CALL gpeq_fcoords(normpsi,foutmn,mfac,mpert,2,0,0,0,0,0)
            DO in=1,mpert
               WRITE(out_unit,'(6(1x,es16.8))')
     $              REAL(normpsi),REAL(mfac(in)),
     $              REAL(fxmn(in)),REAL(foutmn(in)),
     $              AIMAG(fxmn(in)),AIMAG(foutmn(in))
               WRITE(*,*)REAL(fxmn(in)),REAL(foutmn(in)),
     $              AIMAG(fxmn(in)),AIMAG(foutmn(in))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

c         fxfun=cos(twopi*3*theta)
c         CALL iscdftf(mfac,mpert,fxfun,mthsurf,fxmn)
c         CALL iscdftb(mfac,mpert,fxfun,mthsurf,fxmn)
c         CALL iscdftf(mfac,mpert,fxfun,mthsurf,fxmn)
c
         ! test orthoganality of permeabev eigenvectors
c         CALL gpdiag_permeabev_orthogonality
c         ! Test coordinate independence of power eigenvectors
c         CALL gpdiag_reluctpowout(power_rout,power_bpout,power_bout,
c     $        power_rcout)
      ENDIF
c-----------------------------------------------------------------------
c     send harvest record.
c-----------------------------------------------------------------------
      ierr=harvest_send(hlog)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL gpec_dealloc
      CALL gpec_stop("Normal termination.")
      END PROGRAM gpec_main

