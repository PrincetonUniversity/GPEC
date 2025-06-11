c-----------------------------------------------------------------------
c     file dcon.f.
c     performs ideal MHD stability analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. dcon.
c     2. dcon_dealloc.
c-----------------------------------------------------------------------
c     subprogram 1. dcon.
c     performs ideal & kinetic MHD stability analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM dcon
      USE equil_mod
      USE equil_out_mod
      USE bal_mod
      USE mercier_mod
      USE ode_mod
      USE free_mod
      USE resist_mod
      USE pentrc_interface,         ! rename overlapping names
     $    pentrc_verbose=>verbose,  ! should get a more fundamental fix
     $    pentrc_mpert=>mpert,
     $    pentrc_nn=>nn,
     $    pentrc_r8=>r8,
     $    pentrc_timer=>timer
      IMPLICIT NONE

      LOGICAL :: use_classic_splines_for_dcon=.FALSE., cyl_flag=.FALSE.,
     $ reform_eq_with_psilim=.FALSE.
      INTEGER :: mmin,ipsi,m,peak_index
      COMPLEX(r8) :: plasma1,vacuum1,total1

      INTEGER, DIMENSION(:), ALLOCATABLE :: mtmp
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psitmp
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xtmp

      ! harvest variables
      INCLUDE 'harvest_lib.inc77'
      INTEGER  :: ierr
      INTEGER  :: hint = 0
      REAL(r8) :: hdbl = 0
      CHARACTER(LEN=16)    :: hkey
      CHARACTER(LEN=50000) :: hnml
      CHARACTER(LEN=65507) :: hlog
      CHARACTER, PARAMETER :: nul = char(0)

      NAMELIST/dcon_control/bal_flag,mat_flag,ode_flag,vac_flag,
     $     mer_flag,res_flag,fft_flag,node_flag,mthvac,sing_start,nn,
     $     delta_mlow,delta_mhigh,delta_mband,thmax0,nstep,ksing,
     $     tol_nr,tol_r,crossover,ucrit,singfac_min,singfac_max,
     $     cyl_flag,dmlim,lim_flag,sas_flag,sing_order,sort_type,
     $     termbycross_flag,qhigh,kin_flag,con_flag,kinfac1,kinfac2,
     $     kingridtype,ktanh_flag,passing_flag,trapped_flag,
     $     ion_flag,electron_flag,ktc,ktw,qlow,
     $     use_classic_splines,reform_eq_with_psilim,
     $     psiedge,nperq_edge,wv_farwall_flag,
     $     dcon_kin_threads,parallel_threads
      NAMELIST/dcon_output/interp,crit_break,out_bal1,
     $     bin_bal1,out_bal2,bin_bal2,out_metric,bin_metric,out_fmat,
     $     bin_fmat,out_gmat,bin_gmat,out_kmat,bin_kmat,out_sol,
     $     out_sol_min,out_sol_max,bin_sol,bin_sol_min,bin_sol_max,
     $     out_fl,bin_fl,out_evals,bin_evals,bin_euler,euler_stride,
     $     bin_vac,ahb_flag,mthsurf0,msol_ahb,netcdf_out,out_fund,
     $     out_ahg2msc
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/'ipsi',2x,'psifac',6x,'f',7x,'mu0 p',5x,'dvdpsi',
     $     6x,'q',10x,'di',9x,'dr',8x,'ca1'/)
 20   FORMAT(i4,1p,8e11.3)
 30   FORMAT(/3x,"mlow",1x,"mhigh",1x,"mpert",1x,"mband",3x,"nn",2x,
     $     "lim_fl",2x,"dmlim",7x,"qlim",6x,"psilim"//5i6,l6,1p,3e11.3/)
c-----------------------------------------------------------------------
c     read input data.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)""
      IF(verbose) WRITE(*,*)"DCON START => "//TRIM(version)
      IF(verbose) WRITE(*,*)"__________________________________________"
      CALL timer(0,out_unit)
      CALL ascii_open(in_unit,"dcon.in","OLD")
      READ(UNIT=in_unit,NML=dcon_control)
      READ(UNIT=in_unit,NML=dcon_output)
      CALL ascii_close(in_unit)

      IF (out_ahg2msc) THEN
         WRITE(*,*) "WARNING: ahg2msc.out is deprecated and will be " //
     $          "removed in a future version. Set out_ahg2msc = .FALSE."
         WRITE(*,*) "         to disable this warning."
         vac_memory=.FALSE.
      ELSE
         vac_memory=.TRUE.
      ENDIF

c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(parallel_threads > 0)THEN
        WRITE(*,*) "!! WARNING: parallel_threads has been deprecated. "
     $    //"Use dcon_kin_threads"
      ENDIF
#ifdef _OPENMP
c The code is compiled with the OpenMP option so the preprocessor
c will pick include the following line
      IF(dcon_kin_threads > 0)THEN
        CALL OMP_SET_NUM_THREADS(dcon_kin_threads)
      ENDIF
#else
c No OpenMP so the following line will appear in source code
      dcon_kin_threads = 1
#endif
      delta_mhigh=delta_mhigh*2
      use_classic_splines_for_dcon = use_classic_splines
c-----------------------------------------------------------------------
c     open output files, read, process, and diagnose equilibrium.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"dcon.out","UNKNOWN")
      CALL equil_read(out_unit)
      CALL equil_out_global
      CALL equil_out_qfind
c-----------------------------------------------------------------------
c     optionally reform the eq splines to concentrate at true truncation
c-----------------------------------------------------------------------
      CALL sing_lim  ! determine if qhigh is truncating before psihigh
      IF(psilim /= psihigh .AND. reform_eq_with_psilim)THEN
         CALL equil_read(out_unit, psilim)
         CALL equil_out_global
         CALL equil_out_qfind
      ENDIF
c-----------------------------------------------------------------------
c     record the equilibrium properties
c-----------------------------------------------------------------------
      CALL equil_out_diagnose(.FALSE.,out_unit)
      CALL equil_out_write_2d
      IF(dump_flag .AND. eq_type /= "dump")CALL equil_out_dump
      IF(direct_flag)CALL bicube_dealloc(psi_in)
c-----------------------------------------------------------------------
c     enable different settings for dcon and equilibrium splines
c-----------------------------------------------------------------------
      use_classic_splines = use_classic_splines_for_dcon
c-----------------------------------------------------------------------
c     prepare local stability criteria.
c-----------------------------------------------------------------------
      CALL spline_alloc(locstab,mpsi,5)
      locstab%xs=sq%xs
      locstab%fs=0
      locstab%name="locstb"
      locstab%title=(/"  di  ","  dr  ","  h   "," ca1  "," ca2  "/)
      IF (mer_flag) THEN
         WRITE(*,*)"Evaluating Mercier criterion"
         CALL mercier_scan
      ENDIF
      IF(bal_flag)THEN
         IF(verbose) WRITE(*,*)"Evaluating ballooning criterion"
         CALL bal_scan
      ENDIF
      CALL spline_fit(locstab,"extrap")
c-----------------------------------------------------------------------
c     output surface quantities.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,"dcon.bin","UNKNOWN","REWIND","none")
      WRITE(out_unit,10)
      DO ipsi=0,mpsi
         WRITE(out_unit,20)ipsi,sq%xs(ipsi),sq%fs(ipsi,1)/twopi,
     $        sq%fs(ipsi,2),sq%fs(ipsi,3),sq%fs(ipsi,4),
     $        locstab%fs(ipsi,1)/sq%xs(ipsi),
     $        locstab%fs(ipsi,2)/sq%xs(ipsi),locstab%fs(ipsi,4)
         WRITE(bin_unit)
     $        REAL(sq%xs(ipsi),4),
     $        REAL(SQRT(sq%xs(ipsi)),4),
     $        REAL(sq%fs(ipsi,1)/twopi,4),
     $        REAL(sq%fs(ipsi,2),4),
     $        REAL(sq%fs(ipsi,4),4),
     $        REAL(asinh(locstab%fs(ipsi,1)/sq%xs(ipsi)),4),
     $        REAL(asinh(locstab%fs(ipsi,2)/sq%xs(ipsi)),4),
     $        REAL(asinh(locstab%fs(ipsi,3)),4),
     $        REAL(asinh(locstab%fs(ipsi,4)),4)
      ENDDO
      WRITE(out_unit,10)
      WRITE(bin_unit)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     define poloidal mode numbers.
c-----------------------------------------------------------------------
      CALL sing_find
      !CALL sing_lim  ! now done in above test for re-forming equil
      IF(cyl_flag)THEN
         mlow=delta_mlow
         mhigh=delta_mhigh
      ELSEIF(sing_start == 0)THEN
         mlow=MIN(nn*qmin,zero)-4-delta_mlow
         mhigh=nn*qmax+delta_mhigh
      ELSE
         mmin=HUGE(mmin)
         DO ising=INT(sing_start),msing
            mmin=MIN(mmin,sing(ising)%m)
         ENDDO
         mlow=mmin-delta_mlow
         mhigh=nn*qmax+delta_mhigh
      ENDIF
      mpert=mhigh-mlow+1
      mband=mpert-1-delta_mband
      mband=MIN(MAX(mband,0),mpert-1)
c-----------------------------------------------------------------------
c     fit equilibrium quantities to Fourier-spline functions.
c-----------------------------------------------------------------------
      IF(mat_flag .OR. ode_flag)THEN
         IF(verbose) WRITE(*,'(3x,1p,4(a,es10.3))')"q0 = ",q0,
     $        ", qmin = ",qmin,", qmax = ",qmax,", q95 = ",q95
         IF(verbose) WRITE(*,'(3x,a,l1,1p,3(a,es10.3))')
     $        "sas_flag = ",sas_flag,", dmlim = ",dmlim,
     $        ", qlim = ",qlim,", psilim = ",psilim
         IF(verbose) WRITE(*,'(3x,1p,3(a,es10.3))')"betat = ",betat,
     $        ", betan = ",betan,", betap1 = ",betap1
         IF(verbose) WRITE(*,'(3x,5(a,i4))')"nn = ",nn,", mlow = ",mlow,
     $        ", mhigh = ",mhigh,", mpert = ",mpert,", mband = ",mband
         IF(verbose) WRITE(*,'(1x,a)')
     $        "Fourier analysis of metric tensor components"
         CALL fourfit_make_metric
         IF(verbose) WRITE(*,*)"Computing F, G, and K Matrices"
         CALL fourfit_make_matrix(out_fund)
         WRITE(out_unit,30)mlow,mhigh,mpert,mband,nn,sas_flag,dmlim,
     $        qlim,psilim
         IF(kin_flag)THEN
            CALL fourfit_action_matrix
            IF(verbose) WRITE(*,*) "Initializing PENTRC"
            ! call the automatic reading and distributing of inputs
            CALL initialize_pentrc(op_kin=.FALSE.,op_deq=.FALSE.,
     $          op_peq=.FALSE.)
            ! manually set the pentrc equilibrium description
            CALL set_eq(eqfun,sq,rzphi,smats,tmats,xmats,ymats,zmats,
     $           twopi*psio,ro,nn,jac_type,mlow,mhigh,mpert,mthvac)
            ! manually set the kinetic profiles
            CALL read_kin(kinetic_file,zi,zimp,mi,mimp,nfac,
     $          tfac,wefac,wpfac,indebug)
            ! manually set the perturbed equilibrium displacements
            ! use false flat xi and xi' for equal weighting
            ALLOCATE(psitmp(sq%mx+1),mtmp(mpert),xtmp(sq%mx+1,mpert))
            psitmp(:) = sq%xs(0:)
            mtmp = (/(m,m=mlow,mhigh)/)
            xtmp = 1e-4
            CALL set_peq(psitmp,mtmp,xtmp,xtmp,xtmp,.false.,tdebug)
            DEALLOCATE(xtmp,mtmp,psitmp)
            IF(verbose) WRITE(*,*)"Computing Kinetic Matrices"
            CALL fourfit_kinetic_matrix(kingridtype,out_fund)
         ENDIF
         CALL sing_scan
         DO ising=1,msing
            CALL resist_eval(sing(ising))
         ENDDO
         IF (kin_flag)THEN
            CALL ksing_find
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     integrate main ODE's.
c-----------------------------------------------------------------------
      ALLOCATE(ud(mpert,mpert,2))
      IF(ode_flag)THEN
         IF(verbose) WRITE(*,*)"Starting integration of ODE's"
         CALL ode_run
         IF(size_edge > 0)THEN
            peak_index = MAXLOC(REAL(dw_edge(pre_edge:i_edge)),dim=1) +
     $         (pre_edge - 1)
            qhigh = q_edge(peak_index)
            sas_flag = .FALSE.
            psiedge = psihigh
            CALL sing_lim
            WRITE(*,'(1x,a34,f10.6,a10,f10.6)') "Re-Integrating to "//
     $            "peak dW @ qlim =",qlim,", psilim =",psilim
            ! full re-run because outputs were written to disk each step
            ! making it hard to backtrack
            IF(bin_euler)CALL bin_close(euler_bin_unit)
            CALL ode_run
         ENDIF
      ENDIF
      DEALLOCATE(ud)
c-----------------------------------------------------------------------
c     compute free boundary energies.
c-----------------------------------------------------------------------
      IF(vac_flag .AND. .NOT.
     $     (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
         IF(verbose) WRITE(*,*)"Computing free boundary energies"
         ALLOCATE(ud(mpert,mpert,2))
         CALL free_run(plasma1,vacuum1,total1,nzero,netcdf_out)
         DEALLOCATE(ud)
      ELSE
         plasma1=0
         vacuum1=0
         total1=0
         CALL dcon_dealloc
      ENDIF
      IF(mat_flag .OR. ode_flag)DEALLOCATE(asmat,bsmat,csmat,ipiva,jmat)
      IF(bin_euler)CALL bin_close(euler_bin_unit)
c-----------------------------------------------------------------------
c     the bottom line.
c-----------------------------------------------------------------------
      IF(nzero /= 0)THEN
         IF(verbose) WRITE(*,'(1x,a,i2,".")')
     $        "Fixed-boundary mode unstable for nn = ",nn
      ENDIF
      IF(vac_flag .AND. .NOT.
     $        (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
         IF(REAL(total1) < 0)THEN
            IF(verbose) WRITE(*,'(1x,a,i2,".")')
     $           "Free-boundary mode unstable for nn = ",nn
         ELSE
            IF(verbose) WRITE(*,'(1x,a,i2,".")')
     $           "All free-boundary modes stable for nn = ",nn
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     save output in sum1.bin.
c-----------------------------------------------------------------------
      CALL bin_open(sum_unit,"sum1.bin","UNKNOWN","REWIND","none")
      WRITE(sum_unit)mpsi,mtheta,mlow,mhigh,mpert,mband,
     $     REAL(psilow,4),REAL(psihigh,4),REAL(amean,4),REAL(rmean,4),
     $     REAL(aratio,4),REAL(kappa,4),REAL(delta1,4),REAL(delta2,4),
     $     REAL(li1,4),REAL(li2,4),REAL(li3,4),REAL(ro,4),REAL(zo,4),
     $     REAL(psio,4),REAL(betap1,4),REAL(betap2,4),REAL(betap3,4),
     $     REAL(betat,4),REAL(betan,4),REAL(bt0,4),REAL(q0,4),
     $     REAL(qmin,4),REAL(qmax,4),REAL(qa,4),REAL(crnt,4),
     $     REAL(plasma1,4),REAL(vacuum1,4),REAL(total1,4),REAL(q95,4),
     $     REAL(bwall,4),REAL(betaj,4)
      CALL bin_close(sum_unit)
c-----------------------------------------------------------------------
c     log inputs/outputs with harvest
c-----------------------------------------------------------------------
      ierr=init_harvest('CODEDB_DCON'//nul,hlog,LEN(hlog))
      ierr=set_harvest_verbose(0)
      ! standard CODEDB records
      ierr=set_harvest_payload_str(hlog,'CODE'//nul,'DCON'//nul)
      ierr=set_harvest_payload_str(hlog,'VERSION'//nul,version//nul)
      if(shotnum>0)
     $   ierr=set_harvest_payload_int(hlog,'SHOT'//nul,INT(shotnum))
      if(shottime>0)
     $   ierr=set_harvest_payload_int(hlog,'TIME'//nul,INT(shottime))
      ! DCON specifc records
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
      ierr=set_harvest_payload_dbl(hlog,'plasma1'//nul,plasma1)
      ierr=set_harvest_payload_dbl(hlog,'vacuum1'//nul,vacuum1)
      ierr=set_harvest_payload_dbl(hlog,'total1'//nul,total1)
      ierr=harvest_send(hlog)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END PROGRAM dcon
c-----------------------------------------------------------------------
c     subprogram 2. dcon_dealloc.
c     deallocates internal memory.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dcon_dealloc
      USE ode_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     deallocate internal memory.
c-----------------------------------------------------------------------
      CALL spline_dealloc(sq)
      CALL spline_dealloc(locstab)
      CALL bicube_dealloc(rzphi)
      IF(mat_flag .OR. ode_flag)THEN
         CALL cspline_dealloc(amats)
         CALL cspline_dealloc(bmats)
         CALL cspline_dealloc(cmats)
         CALL cspline_dealloc(dmats)
         CALL cspline_dealloc(emats)
         CALL cspline_dealloc(hmats)
         CALL cspline_dealloc(dbats)
         CALL cspline_dealloc(ebats)
         CALL cspline_dealloc(fbats)
         CALL cspline_dealloc(fmats)
         CALL cspline_dealloc(gmats)
         CALL cspline_dealloc(kmats)
      
         DO ising=1,msing
            DEALLOCATE(sing(ising)%vmat)
            DEALLOCATE(sing(ising)%mmat)
            DEALLOCATE(sing(ising)%power)
            DEALLOCATE(sing(ising)%n1)
            DEALLOCATE(sing(ising)%n2)
         ENDDO
         DEALLOCATE(sing)
      ENDIF
      IF(ode_flag)DEALLOCATE(u,du,u_save)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dcon_dealloc
