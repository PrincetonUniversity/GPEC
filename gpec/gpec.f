c-----------------------------------------------------------------------
c     GENERALIZED PERTURBED EQUILIBRIUM CODE
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

      USE diagnostics
      USE gpec_global
      USE local_mod, ONLY: ascii_open, ascii_close, bin_open, bin_close
      USE response
      USE peq
      USE dcon_interface, ONLY: dcon_read, dcon_build, dcon_metric,
     $    dcon_transform, dcon_vacuum
      USE rdcon_interface, ONLY: rdcon_read
      USE coil_mod, ONLY : ipd,btd,helicity,cmlow,cmhigh,cmpert,
     $    machine,ip_direction,bt_direction,
     $    coil_read
      USE field_mod, ONLY : field_bs_psi
      USE output_control, only : response_matrices, singcoup, control,
     $    singfld, vsingfld, init_netcdf
      USE output_cylindrical, only : xbrzphi, vsbrzphi, xbrzphifun,
     $    arzphifun
      USE output_profile, only : dw_profile, dw_matrix, pmodb,
     $    xbnormal, vbnormal, xbtangent, xclebsch

      IMPLICIT NONE

      INTEGER :: i,in,osing,resol,angnum,
     $     mthnumb,meas,mode,m3mode,lowmode,highmode,filter_modes
      INTEGER :: pmode,p1mode,rmode,dmode,d1mode,fmode,smode
      INTEGER, DIMENSION(:), POINTER :: ipiv
      REAL(r8) :: majr,minr,rdist,smallwidth,factor,fp,normpsi
      CHARACTER(8) :: filter_types
      CHARACTER(128) :: infile
      LOGICAL :: singcoup_flag,singfld_flag,vsingfld_flag,pmodb_flag,
     $     xbcontra_flag,xbnormal_flag,vbnormal_flag,xbnobo_flag,
     $     d3_flag,xbst_flag,pmodbrz_flag,rzphibx_flag,dw_flag,
     $     radvar_flag,eigen_flag,magpot_flag,xbtangent_flag,
     $     arbsurf_flag,angles_flag,surfmode_flag,rzpgrid_flag,
     $     singcurs_flag,m3d_flag,cas3d_flag,test_flag,nrzeq_flag,
     $     arzphifun_flag,xbrzphifun_flag,pmodbmn_flag,xclebsch_flag,
     $     filter_flag
      LOGICAL, DIMENSION(100) :: ss_flag
      COMPLEX(r8), DIMENSION(:), POINTER :: finmn,foutmn,xspmn,
     $     fxmn,fxfun,coilmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: invmats,temp1

      ! harvest variables
      INCLUDE 'harvest_lib.inc77'
      INTEGER  :: ierr
      CHARACTER(LEN=50000) :: hnml
      CHARACTER(LEN=65507) :: hlog
      CHARACTER, PARAMETER :: nul = char(0)

      ! input file namelists
      NAMELIST/gpec_input/dcon_dir,ieqfile,idconfile,ivacuumfile,
     $     power_flag,fft_flag,mthsurf0,fixed_boundary_flag,
     $     data_flag,data_type,nmin,nmax,mmin,mmax,jsurf_in,mthsurf,
     $     jac_in,power_bin,power_rin,power_bpin,power_rcin,tmag_in,
     $     infile,harmonic_flag,mode_flag,sinmn,cosmn,
     $     displacement_flag,mode,coil_flag,
     $     ip_direction,bt_direction,rdconfile,
     $     pmode,p1mode,dmode,d1mode,fmode,rmode,smode,
     $     filter_types,filter_modes
      NAMELIST/gpec_control/resp_index,sing_spot,reg_flag,reg_spot,
     $     chebyshev_flag,nche,nchr,nchz,resp_induct_flag
      NAMELIST/gpec_output/resp_flag,singcoup_flag,nrzeq_flag,nr,nz,
     $     singfld_flag,pmodb_flag,xbnormal_flag,rstep,jsurf_out,
     $     jac_out,power_bout,power_rout,power_bpout,power_rcout,
     $     tmag_out,eqbrzphi_flag,brzphi_flag,xrzphi_flag,
     $     vbrzphi_flag,vvbrzphi_flag,divzero_flag,dw_flag,
     $     bin_flag,bin_2d_flag,fun_flag,flux_flag,bwp_pest_flag,
     $     vsbrzphi_flag,ss_flag,arzphifun_flag,xbrzphifun_flag,
     $     vsingfld_flag,vbnormal_flag,eigm_flag,xbtangent_flag,
     $     xclebsch_flag,pbrzphi_flag,verbose,max_linesout,filter_flag
      NAMELIST/gpec_diagnose/singcurs_flag,xbcontra_flag,
     $     xbnobo_flag,d3_flag,div_flag,xbst_flag,pmodbrz_flag,
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
      rdconfile="gal_solution.bin"
      power_flag=.TRUE.
      fft_flag=.FALSE.
      fixed_boundary_flag=.FALSE.
      data_flag=.FALSE.
      harmonic_flag=.FALSE.
      mode_flag=.FALSE.
      displacement_flag=.FALSE.
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
      reg_flag=.TRUE.
      reg_spot=5e-2
      chebyshev_flag=.TRUE.
      nche=20
      nchr=20
      nchz=20

      jsurf_out=0
      tmag_out=1
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
      fun_flag=.FALSE.
      flux_flag=.FALSE.
      max_linesout=0
      vsbrzphi_flag=.FALSE.
      DO i=1,100 
         ss_flag(i)=.FALSE.
      ENDDO
      arzphifun_flag=.FALSE.
      xbrzphifun_flag=.FALSE.
      bwp_pest_flag=.FALSE.
      dw_flag = .FALSE.

      singcurs_flag=.FALSE.
      xbcontra_flag=.FALSE.
      xbnobo_flag=.FALSE.
      d3_flag=.FALSE.
      div_flag=.FALSE.
      xbst_flag=.FALSE.
      pmodbrz_flag=.FALSE.
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
      bwp_pest_flag=.FALSE.

      majr=10.0
      minr=1.0
      mode=1
      lowmode=-10
      highmode=10
      m3mode=2
      resol=1e4
      smallwidth=1e-6
      
      timeit = .FALSE.
      verbose = .TRUE.
      debug_flag = .FALSE.
      malias=0
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
c     prepare ideal solutions.
c-----------------------------------------------------------------------
      CALL dcon_read(psixy)
      CALL dcon_transform
c-----------------------------------------------------------------------
c     reconstruct metric tensors.
c-----------------------------------------------------------------------
      CALL dcon_metric
c-----------------------------------------------------------------------
c     read vacuum data.
c-----------------------------------------------------------------------
      CALL dcon_vacuum
      IF(timeit) CALL gpec_timer(2)
c-----------------------------------------------------------------------
c     set parameters from dcon.
c-----------------------------------------------------------------------
      ALLOCATE(xspmn(mpert),finmn(mpert),foutmn(mpert),
     $     fxmn(mpert),fxfun(mthsurf))
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
         coilmn=0
         CALL field_bs_psi(psilim,coilmn,1)
         DO i=1,cmpert
            IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
               finmn(cmlow-mlow+i)=coilmn(i)
            ENDIF
         ENDDO
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
      !machine = to_upper(machine)
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
      ! gpec inputs
      write(hnml,nml=gpec_input)
      ierr=set_harvest_payload_nam(hlog,'GPEC_INPUT'//nul,
     $                             trim(hnml)//nul)
      write(hnml,nml=gpec_control)
      ierr=set_harvest_payload_nam(hlog,'GPEC_CONTROL'//nul,
     $                             trim(hnml)//nul)
      write(hnml,nml=gpec_output)
      ierr=set_harvest_payload_nam(hlog,'GPEC_OUTPUT'//nul,
     $                             trim(hnml)//nul)
c-----------------------------------------------------------------------
c     compute plasma response.
c-----------------------------------------------------------------------
      CALL response_eigen
      IF(timeit) CALL gpec_timer(2)
      CALL response_pinduct
      CALL response_sinduct
      CALL response_permeab
      CALL response_reluct
      IF(timeit) CALL gpec_timer(2)
c-----------------------------------------------------------------------
c     run and test rdcon.
c-----------------------------------------------------------------------
c      CALL rdcon_read
c-----------------------------------------------------------------------
c     Set parameters for outputs.
c-----------------------------------------------------------------------
      IF (nrzeq_flag) THEN
         nr=mr
         nz=mz
      ENDIF
      CALL peq_rzpgrid(nr,nz,psixy)
c-----------------------------------------------------------------------
c     full analysis.
c-----------------------------------------------------------------------
      CALL init_netcdf
      IF (resp_flag) THEN
         CALL response_matrices(power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out,jsurf_out)
      ENDIF
      DO i=1,LEN_TRIM(filter_types)
         IF(filter_types(i:i)=='s') singcoup_flag=.TRUE.
      ENDDO
      IF (singcoup_flag) THEN
         CALL singcoup(sing_spot,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and external field.
c-----------------------------------------------------------------------
      IF (data_flag.OR.harmonic_flag.OR.coil_flag)THEN
         edge_flag=.TRUE.
         CALL control(infile,finmn,foutmn,xspmn,power_rin,
     $        power_bpin,power_bin,power_rcin,tmag_in,jsurf_in,
     $        power_rout,power_bpout,power_bout,power_rcout,tmag_out,
     $        filter_types,filter_modes,filter_flag)
      ELSE IF (mode_flag) THEN
         edge_flag=.FALSE.
      ENDIF

      IF (singfld_flag) THEN
         CALL singfld(mode,xspmn,sing_spot,power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out,singcoup_flag)
      ENDIF
      IF (coil_flag .AND. vsingfld_flag) THEN
         CALL vsingfld(power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out)
      ENDIF
      IF (xclebsch_flag) THEN
         CALL xclebsch(mode,xspmn)
      ENDIF
      IF (kin_flag .AND. dw_flag) THEN
         CALL dw_profile(mode,xspmn)
         CALL dw_matrix
      ENDIF
      IF (pmodb_flag) THEN
         CALL pmodb(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (xbnormal_flag) THEN
         CALL xbnormal(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (xbtangent_flag) THEN
         CALL xbtangent(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (coil_flag .AND. vbnormal_flag) THEN
         CALL vbnormal(power_rout,power_bpout,power_bout,
     $        power_rcout,tmag_out)
      ENDIF
      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag .OR. vvbrzphi_flag .OR. pbrzphi_flag) THEN
         IF (.NOT.mode_flag) THEN
            CALL xbrzphi(mode,xspmn,nr,nz,finmn,foutmn)
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
            CALL dcon_build(mode,xspmn)
            CALL peq_alloc
            CALL peq_sol(psilim)
            CALL peq_contra(psilim)
            CALL peq_dealloc
            CALL peq_weight(psilim,foutmn,mfac,mpert,1)
            finmn = MATMUL(invmats,foutmn)
            CALL peq_weight(psilim,foutmn,mfac,mpert,0)
            CALL peq_weight(psilim,finmn,mfac,mpert,0)
            CALL xbrzphi(mode,xspmn,nr,nz,finmn,foutmn)
            DEALLOCATE(ipiv,invmats,temp1)
         ENDIF
      ENDIF
      IF (singfld_flag .AND. vsbrzphi_flag) THEN
         DO i=1,msing
            IF (ss_flag(i)) CALL vsbrzphi(i,nr,nz)
         ENDDO
      ENDIF

      IF (xbrzphifun_flag) THEN
         CALL xbrzphifun(mode,xspmn)
      ENDIF
      IF (arzphifun_flag) THEN
         CALL arzphifun(mode,xspmn)
      ENDIF

c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF (singcurs_flag) THEN
         CALL diagnose_singcurs(mode,xspmn,msing,resol,smallwidth)
      ENDIF
      IF (xbcontra_flag) THEN
         CALL diagnose_xbcontra(mode,xspmn,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out)
      ENDIF
      IF (xbnobo_flag) THEN
         CALL diagnose_xbnobo(mode,xspmn,d3_flag)
      ENDIF
      IF (xbst_flag) THEN
         CALL diagnose_xbst(mode,xspmn)
      ENDIF
      IF (pmodbrz_flag) THEN
         CALL diagnose_pmodbrz(mode,xspmn)
         CALL diagnose_pmodbmn(mode,xspmn)
      ENDIF            
      IF (rzphibx_flag) THEN
         CALL diagnose_rzphibx(mode,xspmn)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose without a given error field.
c-----------------------------------------------------------------------
      IF (radvar_flag) THEN
         CALL diagnose_radvar
      ENDIF
      IF (eigen_flag) THEN
         CALL diagnose_eigen
      ENDIF
      IF (magpot_flag) THEN
         CALL diagnose_magpot
      ENDIF
      IF (arbsurf_flag) THEN
         CALL diagnose_arbsurf(majr,minr)
      ENDIF         
      IF (angles_flag) THEN
         CALL diagnose_angles
      ENDIF

      IF (surfmode_flag) THEN
         CALL diagnose_surfmode(lowmode,highmode,power_rout,power_bpout,
     $        power_bout,power_rcout,tmag_out,jsurf_out)
      ENDIF
      IF (rzpgrid_flag) THEN
         CALL diagnose_rzpgrid(nr,nz)
      ENDIF
      
      IF (m3d_flag) THEN
         normpsi=1.0
         fp=1e-3         
         fxmn=0
         fxmn(m3mode-mlow+1)=fp*normpsi
         CALL peq_fcoords(psilim,fxmn,mfac,mpert,0,1,0,1,0,0)
         fxmn=-twopi*ifac*chi1*(mfac-nn*qlim)*fxmn
         CALL peq_weight(psilim,fxmn,mfac,mpert,0)
         CALL control(infile,fxmn,foutmn,xspmn,
     $        0,0,0,0,1,0,0,0,0,0,1,'   ',0,.FALSE.)
         edge_flag=.TRUE.
         CALL singfld(mode,xspmn,sing_spot,power_rout,
     $        power_bpout,power_bout,power_rcout,tmag_out,.FALSE.)
      ENDIF

      IF (cas3d_flag) THEN
         fp = -1e-2
         
         fxmn=0
         fxmn(m3mode-mlow+1)=fp
         CALL peq_fcoords(psilim,fxmn,mfac,mpert,0,0,2,0,1,0)
         CALL control(infile,finmn,foutmn,xspmn,power_rin,
     $        power_bpin,power_bin,power_rcin,tmag_in,jsurf_in,
     $        power_rout,power_bpout,power_bout,power_rcout,
     $        tmag_out,'   ',0,.FALSE.)
         edge_flag=.TRUE.
         CALL singfld(mode,xspmn,sing_spot,0,0,0,0,1,.FALSE.)
         CALL diagnose_xbcontra(mode,xspmn,0,0,0,0,1)
         CALL diagnose_xbcontra(mode,xspmn,0,0,2,0,1)
         CALL xbnormal(mode,xspmn,0,0,0,0,1)
         CALL xbnormal(mode,xspmn,0,0,2,0,1)
         CALL diagnose_xbnobo(mode,xspmn,d3_flag)
         CALL diagnose_radvar
      ENDIF
c-----------------------------------------------------------------------
c     various simple test.
c-----------------------------------------------------------------------
      IF (test_flag) THEN
         fxmn=0
         fxmn(7-mlow+1)=1.0
         CALL ascii_open(out_unit,"iptest_coordtrans_n3.out","UNKNOWN")
         WRITE(out_unit,*)
     $        "IPTEST_COORDTRANS: cotoha to hatoco, co(0,0)"
         WRITE(out_unit,'(5(1x,a16))')"normpsi","mfac","fxmn",
     $        "real","imag"
         DO i=1,100
            foutmn=fxmn 
            normpsi = REAL(i)/100.0
            CALL peq_bcoords(normpsi,foutmn,mfac,mpert,0,0,0,0,1,0)
            DO in=1,mpert
               WRITE(out_unit,'(5(1x,es16.8))')
     $              REAL(normpsi),REAL(mfac(in)),
     $              REAL(fxmn(in)),REAL(foutmn(in)),AIMAG(foutmn(in))
            ENDDO            
         ENDDO
         CALL ascii_close(out_unit)
        
         fxfun=cos(twopi*3*theta)
         CALL iscdftf(mfac,mpert,fxfun,mthsurf,fxmn)
         CALL iscdftb(mfac,mpert,fxfun,mthsurf,fxmn)
         CALL iscdftf(mfac,mpert,fxfun,mthsurf,fxmn)

         ! test orthoganality of permeabev eigenvectors
         CALL diagnose_permeabev_orthogonality
         ! Test coordinate independence of power eigenvectors
         CALL diagnose_reluctpowout(power_rout,power_bpout,power_bout,
     $        power_rcout)
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
