c-----------------------------------------------------------------------
c     file stride.f.
c     performs ideal MHD stability analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. stride.
c     2. stride_dealloc.
c-----------------------------------------------------------------------
c     subprogram 1. stride.
c     performs ideal MHD stability analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM stride
      USE equil_mod
      USE equil_out_mod
      USE stride_bal_mod
      USE stride_mercier_mod
      USE stride_ode_mod
      USE stride_free_mod
      USE omp_lib
      IMPLICIT NONE

      LOGICAL :: cyl_flag=.FALSE., use_notaknot_splines=.TRUE.,
     $           use_classic_splines_for_stride,
     $           reform_eq_with_psilim=.FALSE.
      INTEGER :: mmin
      REAL(r8) :: plasma1,vacuum1,total1
      REAL(r8) :: psilow_tmp, psilim_tmp

      CHARACTER(len=64) :: arg_str1, arg_str2
      INTEGER :: nArg
      INTEGER :: s0Time, sTime, fTime, cr

      INTEGER :: ithread, mainlevelThreads

      NAMELIST/stride_control/bal_flag,mat_flag,ode_flag,vac_flag,
     $     mer_flag,fft_flag,node_flag,mthvac,sing_start,nn,
     $     delta_mlow,delta_mhigh,delta_mband,thmax0,nstep,ksing,
     $     tol_nr,tol_r,crossover,ucrit,singfac_min,singfac_max,
     $     cyl_flag,dmlim,lim_flag,sas_flag,sing_order,
     $     use_classic_splines,use_notaknot_splines,qlow,qhigh,
     $     reform_eq_with_psilim
      NAMELIST/stride_output/interp,crit_break,out_bal1,
     $     bin_bal1,out_bal2,bin_bal2,out_metric,bin_metric,out_fmat,
     $     bin_fmat,out_gmat,bin_gmat,out_kmat,bin_kmat,out_sol,
     $     out_sol_min,out_sol_max,bin_sol,bin_sol_min,bin_sol_max,
     $     out_fl,bin_fl,out_evals,bin_evals,bin_euler,euler_stride,
     $     ahb_flag,mthsurf0,msol_ahb,netcdf_out,out_ahg2msc
      NAMELIST/stride_params/grid_packing,asymp_at_sing,
     $     integrate_riccati,calc_delta_prime,calc_dp_with_vac,
     $     axis_mid_pt_skew,
     $     big_soln_err_tol, kill_big_soln_for_ideal_dW,
     $     ric_dt,ric_tol,riccati_bounce,verbose_riccati_output,
     $     riccati_match_hamiltonian_evals,verbose_performance_output,
     $     fourfit_metric_parallel,vac_parallel,nIntervalsTot,nThreads
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/'ipsi',2x,'psifac',6x,'f',7x,'mu0 p',5x,'dvdpsi',
     $     6x,'q',10x,'di',9x,'dr',8x,'ca1'/)
 20   FORMAT(i3,1p,5e10.3,3e11.3)
 30   FORMAT(/3x,"mlow",1x,"mhigh",1x,"mpert",1x,"mband",3x,"nn",2x,
     $     "lim_fl",2x,"dmlim",7x,"qlim",6x,"psilim"//5i6,l6,1p,3e11.3/)

      CALL SYSTEM_CLOCK(COUNT_RATE=cr)
      CALL SYSTEM_CLOCK(COUNT=s0Time)
      sTime = s0Time
c-----------------------------------------------------------------------
c     read input data.
c-----------------------------------------------------------------------
      PRINT *, ""
      PRINT *, "STRIDE START => "//TRIM(version)
      PRINT *, "__________________________________________"
      CALL timer(0,out_unit)
      CALL ascii_open(in_unit,"stride.in","OLD")
      READ(UNIT=in_unit,NML=stride_control)
      REWIND(UNIT=in_unit)
      READ(UNIT=in_unit,NML=stride_output)
      REWIND(UNIT=in_unit)
      READ(UNIT=in_unit,NML=stride_params)
      CALL ascii_close(in_unit)

      IF (out_ahg2msc) THEN
         WRITE(*,*) "WARNING: ahg2msc.out is deprecated and will be " //
     $          "removed in a future version. Set out_ahg2msc = .FALSE."
         WRITE(*,*) "         to disable this warning."
         stride_vac_memory=.FALSE.
      ELSE
         stride_vac_memory=.TRUE.
      ENDIF
         
c-----------------------------------------------------------------------
c     check datatypes of command line inputs.
c-----------------------------------------------------------------------
      nArg = iargc()
      IF (nArg /=0 .AND. nArg /= 2) THEN
         CALL program_stop("Must input two integer arguments"//
     $        "--nIntervalsTot nThreads")
      ELSEIF (nArg == 2) THEN
         CALL getarg(1,arg_str1)
         IF (VERIFY(arg_str1,"0123456789") /= LEN_TRIM(arg_str1)+1) THEN
            CALL program_stop("Argument must be an integer.")
         ENDIF
         CALL getarg(2,arg_str2)
         IF (VERIFY(arg_str2,"0123456789") /= LEN_TRIM(arg_str2)+1) THEN
            CALL program_stop("Argument must be an integer.")
         ENDIF
         READ (arg_str1,'(I10)') nIntervalsTot
         READ (arg_str2,'(I10)') nThreads
      ENDIF
c-----------------------------------------------------------------------
c     set spline methodology for program.
c-----------------------------------------------------------------------
      IF (use_notaknot_splines) THEN
         spline_str="not-a-knot"
      ELSE
         ! ok because fortran standard fills the shorter string with spaces when checking equality
         spline_str="extrap"
      ENDIF
c-----------------------------------------------------------------------
c     check consistency of input settings.
c-----------------------------------------------------------------------
      IF (calc_delta_prime) THEN
         IF (.NOT. asymp_at_sing) THEN
            CALL program_stop("Cannot calculate delta prime w/o " //
     $           "asymptotic expansions at singular surfaces.")
         ELSEIF (integrate_riccati) THEN
            CALL program_stop("Cannot calculate delta prime with " //
     $           "riccati integration.  Turn one of them off.")
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     open output files, read, process, and diagnose equilibrium.
c-----------------------------------------------------------------------
      use_classic_splines_for_stride = use_classic_splines ! allow equil to use or not
      CALL ascii_open(out_unit,"stride.out","UNKNOWN")
      CALL equil_read(out_unit)
      IF(dump_flag .AND. eq_type /= "dump")CALL equil_out_dump
      CALL equil_out_global
      CALL equil_out_qfind
c-----------------------------------------------------------------------
c     optionally reform the eq splines to concentrate at true truncation
c-----------------------------------------------------------------------
      CALL sing_lim  ! determine if qhigh is truncating before psihigh
      CALL sing_min  ! dettermine if qlow excludes more of the core
      ! Unlike DCON, we force a resplining.
      ! The equil_out_qfind propogates to sing_find and then to the
      ! parallelized intervals in a complicated web. Thus, it is not
      ! sufficient to simply set the axisPsi to "start" the ODE somewhere new
      IF(.NOT. reform_eq_with_psilim)THEN
         PRINT *, "** STRIDE requires reformation of equil splines "//
     $            "on q-based sub-interval."
         PRINT *, "  > Forcing reform_eq_with_psilim=t"
      ENDIF
      IF(psilim /= psihigh .OR. psilow /= sq%xs(0))THEN
         psilow_tmp = psilow  ! if we feed psilow directly, it gets overwritten by namelist read
         psilim_tmp = psilim
         CALL equil_read(out_unit, psilim_tmp, psilow_tmp)
         CALL equil_out_global
         CALL equil_out_qfind
      ENDIF
c-----------------------------------------------------------------------
c     record the equilibrium properties
c-----------------------------------------------------------------------
      CALL equil_out_diagnose(.FALSE.,out_unit)
      CALL equil_out_write_2d
      IF(direct_flag)CALL bicube_dealloc(psi_in)
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** equil-input time=",
     $                REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
      use_classic_splines = use_classic_splines_for_stride  ! optionaly different from equil
c-----------------------------------------------------------------------
c     prepare local stability criteria.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
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
         WRITE(*,*)"Evaluating ballooning criterion"
         CALL bal_scan
      ENDIF
      CALL spline_fit(locstab,spline_str)
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** locstab time=",REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     define poloidal mode numbers.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
      CALL sing_find
      delta_mhigh=delta_mhigh*2  ! added for consistency with dcon
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
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** sing time=",REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     fit equilibrium quantities to Fourier-spline functions.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
      IF(mat_flag .OR. ode_flag)THEN
         WRITE(*,'(3x,1p,4(a,e10.3))')"q0 = ",q0,", qmin = ",qmin,
     $        ", qmax = ",qmax,", qa = ",qa
         WRITE(*,'(3x,1p,3(a,e10.3))')"betat = ",betat,
     $        ", betan = ",betan,", betaj = ",betaj
         WRITE(*,'(3x,5(a,i3))')"nn = ",nn,", mlow = ",mlow,
     $        ", mhigh = ",mhigh,", mpert = ",mpert,", mband = ",mband
         WRITE(*,'(1x,a)')"Fourier analysis of metric tensor components"
         CALL fourfit_make_metric
         WRITE(*,*)"Computing F, G, and K Matrices"
         CALL SYSTEM_CLOCK(COUNT=fTime)
         IF (verbose_performance_output) THEN
            print *,"*** fourfit-mid time=",
     $           REAL(fTime-sTime,r8)/REAL(cr,r8)
         ENDIF
         CALL fourfit_make_matrix
         WRITE(out_unit,30)mlow,mhigh,mpert,mband,nn,sas_flag,dmlim,
     $        qlim,psilim
      ENDIF
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** fourfit-tot time=",
     $                REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     prepare code to split: [vac][all else]
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)

#ifdef _OPENMP

#else
      nThreads = 0
      mainlevelThreads = 1
#endif

      !Note: nThreads+1 code is faster than nThreads.  It allows other
      !      modules to spawn nThreads, despite vac_parallel using 1!
      !      Even if nThread+1 > nProc, experiments suggest it's faster.
      CALL sing_parallel_alloc(nThreads+1)

#ifdef _OPENMP
      CALL OMP_SET_NESTED(.TRUE.)
      IF (vac_parallel .AND. nThreads > 1) THEN
         mainlevelThreads = 2
      ELSE
         mainlevelThreads = 1
      ENDIF
      CALL OMP_SET_NUM_THREADS(mainlevelThreads)
#endif

      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** thread-alloc time=",
     $             REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     split code.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithread)
!$OMP DO SCHEDULE(STATIC,1)
      DO ithread = 0,1
         IF (ithread == 0) THEN
c-----------------------------------------------------------------------
c     generate asymptotic expansions at singular surfaces.
c-----------------------------------------------------------------------
            IF (asymp_at_sing) THEN
               CALL sing_asymp
            ENDIF
c-----------------------------------------------------------------------
c     integrate main ODE's.
c-----------------------------------------------------------------------
            IF(ode_flag)THEN
               WRITE(*,*)"Starting ODE integration..."
               WRITE(*, '(1x,a5,2(a11))') "ising",'left','right'
               WRITE(*, '(1x,a27)')'---------------------------'
               CALL ode_run
               IF (integrate_riccati) THEN
                  CALL free_calc_wp
               ELSE
                  CALL ode_propagate_LR
                  CALL ode_calc_modes_for_wp
                  CALL free_calc_wp
                  IF (calc_delta_prime .AND..NOT. calc_dp_with_vac) THEN
                     CALL ode_set_delta_prime_edge
                     CALL ode_calc_delta_prime
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ithread == 1) THEN
c-----------------------------------------------------------------------
c     compute wv vacuum response matrix.
c-----------------------------------------------------------------------
            IF(vac_flag .AND. .NOT.
     $           (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
               CALL free_calc_wv
            ENDIF
         ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print*,"*** tot-parallel time=",
     $                   REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
      IF (ode_flag .AND. calc_delta_prime .AND. calc_dp_with_vac) THEN
c-----------------------------------------------------------------------
c     get Delta' matrix of L/R small soln coeffs driven by L/R big solns
c-----------------------------------------------------------------------
         CALL ode_set_delta_prime_edge(wv)
         CALL ode_calc_delta_prime
      ENDIF
c-----------------------------------------------------------------------
c     calculate free energies.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
      IF(vac_flag .AND. .NOT.
     $     (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
         WRITE(*,*)"Computing free boundary energies"
         CALL free_run(plasma1,vacuum1,total1,netcdf_out)
      ELSE
         plasma1=0
         vacuum1=0
         total1=0
         CALL stride_dealloc
      ENDIF
      IF(mat_flag .OR. ode_flag)DEALLOCATE(amat,bmat,cmat,ipiva,jmat)
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** free-run time=",REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     the bottom line.
c-----------------------------------------------------------------------
      IF(nzero /= 0)THEN
         WRITE(*,'(1x,a,i2,".")')
     $        "Fixed-boundary mode unstable for nn = ",nn
      ENDIF
      IF(vac_flag .AND. .NOT.
     $        (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
         IF(total1 < 0)THEN
            WRITE(*,'(1x,a,i2,".")')
     $           "Free-boundary mode unstable for nn = ",nn
         ELSE
            WRITE(*,'(1x,a,i2,".")')
     $           "All free-boundary modes stable for nn = ",nn
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** complete-run time=",
     $        REAL(fTime-s0Time,r8)/REAL(cr,r8)
      ENDIF
      CALL program_stop("Normal termination.")
      END PROGRAM stride
c-----------------------------------------------------------------------
c     subprogram 2. stride_dealloc.
c     deallocates internal memory.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE stride_dealloc
      USE stride_ode_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     deallocate internal memory.
c-----------------------------------------------------------------------
      CALL spline_dealloc(sq)
      CALL spline_dealloc(locstab)
      CALL bicube_dealloc(rzphi)
      IF(mat_flag .OR. ode_flag)THEN
         CALL cspline_dealloc(fmats)
         CALL cspline_dealloc(gmats)
         CALL cspline_dealloc(kmats)
         IF (asymp_at_sing) THEN
            DO ising=1,msing
               DEALLOCATE(sing(ising)%vmatl)
               DEALLOCATE(sing(ising)%mmatl)
               DEALLOCATE(sing(ising)%vmatr)
               DEALLOCATE(sing(ising)%mmatr)
               DEALLOCATE(sing(ising)%n1)
               DEALLOCATE(sing(ising)%n2)
               DEALLOCATE(sing(ising)%power)
            ENDDO
         ENDIF
         DEALLOCATE(sing)
         IF (.NOT. integrate_riccati) THEN
            DEALLOCATE(uAxis,uAxisD)
            IF (calc_delta_prime) THEN
               DEALLOCATE(uEdge)
            ENDIF
         ENDIF
      ENDIF
      IF(ode_flag)THEN
         DEALLOCATE(u,uFM_all,uFM_sing_inv,uFM_sing_init,
     $        psiInters,psiDirs,psiPoints,scalc)
         IF (calc_delta_prime) THEN
            DEALLOCATE(delta_prime_mat)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stride_dealloc
