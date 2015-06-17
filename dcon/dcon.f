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
c     performs ideal MHD stability analysis.
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
      IMPLICIT NONE

      LOGICAL :: cyl_flag=.FALSE.
      INTEGER :: mmin,ipsi
      REAL(r8) :: plasma1,vacuum1,total1

      NAMELIST/dcon_control/bal_flag,mat_flag,ode_flag,vac_flag,
     $     res_flag,fft_flag,node_flag,mthvac,sing_start,nn,
     $     delta_mlow,delta_mhigh,delta_mband,thmax0,nstep,ksing,
     $     tol_nr,tol_r,crossover,ucrit,singfac_min,singfac_max,
     $     cyl_flag,dmlim,lim_flag,sas_flag,sing_order,sort_type,
     $     termbycross_flag
      NAMELIST/dcon_output/interp,crit_break,out_bal1,
     $     bin_bal1,out_bal2,bin_bal2,out_metric,bin_metric,out_fmat,
     $     bin_fmat,out_gmat,bin_gmat,out_kmat,bin_kmat,out_sol,
     $     out_sol_min,out_sol_max,bin_sol,bin_sol_min,bin_sol_max,
     $     out_fl,bin_fl,out_evals,bin_evals,bin_euler,euler_stride,
     $     ahb_flag,mthsurf0,msol_ahb
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/'ipsi',2x,'psifac',6x,'f',7x,'mu0 p',5x,'dvdpsi',
     $     6x,'q',10x,'di',9x,'dr',8x,'ca1'/)
 20   FORMAT(i3,1p,5e10.3,3e11.3)
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
      delta_mhigh=delta_mhigh*2
c-----------------------------------------------------------------------
c     open output files, read, process, and diagnose equilibrium.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"dcon.out","UNKNOWN")
      CALL equil_read(out_unit)
      IF(dump_flag .AND. eq_type /= "dump")CALL equil_out_dump
      CALL equil_out_global
      CALL equil_out_qfind
      CALL equil_out_diagnose(.FALSE.,out_unit)
      CALL equil_out_write_2d
      IF(direct_flag)CALL bicube_dealloc(psi_in)
c-----------------------------------------------------------------------
c     prepare local stability criteria.
c-----------------------------------------------------------------------
      CALL spline_alloc(locstab,mpsi,5)
      locstab%xs=sq%xs
      locstab%fs=0
      locstab%name="locstb"
      locstab%title=(/"  di  ","  dr  ","  h   "," ca1  "," ca2  "/)
      IF(verbose) WRITE(*,*)"Evaluating Mercier criterion"
      CALL mercier_scan
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
      CALL sing_lim
      IF(cyl_flag)THEN
         mlow=delta_mlow
         mhigh=delta_mhigh
      ELSEIF(sing_start == 0)THEN
         mlow=MIN(nn*qmin,zero)-4-delta_mlow
         mhigh=nn*qmax+delta_mhigh
      ELSE
         mmin=HUGE(mmin)
         DO ising=sing_start,msing
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
         IF(verbose) WRITE(*,'(1x,a)')
     $        "Fourier analysis of metric tensor components"
         IF(verbose) WRITE(*,'(1x,1p,4(a,e10.3))')"q0 = ",q0,
     $        ", qmin = ",qmin,", qmax = ",qmax,", q95 = ",q95
         IF(verbose) WRITE(*,'(1x,a,l1,1p,3(a,e10.3))')
     $        "sas_flag = ",sas_flag,", dmlim = ",dmlim,
     $        ", qlim = ",qlim,", psilim = ",psilim
         IF(verbose) WRITE(*,'(1x,1p,3(a,e10.3))')"betat = ",betat,
     $        ", betan = ",betan,", betap1 = ",betap1
         IF(verbose) WRITE(*,'(1x,5(a,i3))')"nn = ",nn,", mlow = ",mlow,
     $        ", mhigh = ",mhigh,", mpert = ",mpert,", mband = ",mband
         CALL fourfit_make_metric
         IF(verbose) WRITE(*,*)"Computing F, G, and K Matrices"
         CALL fourfit_make_matrix
         WRITE(out_unit,30)mlow,mhigh,mpert,mband,nn,sas_flag,dmlim,
     $        qlim,psilim
         CALL sing_scan
         DO ising=1,msing
            CALL resist_eval(sing(ising))
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     integrate main ODE's.
c-----------------------------------------------------------------------
      IF(ode_flag)THEN
         IF(verbose) WRITE(*,*)"Starting integration of ODE's"
         CALL ode_run
      ENDIF
c-----------------------------------------------------------------------
c     compute free boundary energies.
c-----------------------------------------------------------------------
      IF(vac_flag .AND. .NOT.
     $     (ksing > 0 .AND. ksing <= msing+1 .AND. bin_sol))THEN
         IF(verbose) WRITE(*,*)"Computing free boundary energies"
         CALL free_run(plasma1,vacuum1,total1,nzero)
      ELSE
         plasma1=0
         vacuum1=0
         total1=0
         CALL dcon_dealloc
      ENDIF
      IF(mat_flag .OR. ode_flag)DEALLOCATE(amat,bmat,cmat,ipiva,jmat)
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
         IF(total1 < 0)THEN
            IF(verbose) WRITE(*,'(1x,a,i2,".")')
     $           "Free-boundary mode unstable for nn = ",nn
         ELSE
            IF(verbose) WRITE(*,'(1x,a,i2,".")')
     $           "All free-boundary modes stable for nn = ",nn
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     save output in sum1.dat.
c-----------------------------------------------------------------------
      CALL bin_open(sum_unit,"sum1.dat","UNKNOWN","REWIND","none")
      WRITE(sum_unit)mpsi,mtheta,mlow,mhigh,mpert,mband,
     $     REAL(psilow,4),REAL(psihigh,4),REAL(amean,4),REAL(rmean,4),
     $     REAL(aratio,4),REAL(kappa,4),REAL(delta1,4),REAL(delta2,4),
     $     REAL(li1,4),REAL(li2,4),REAL(li3,4),REAL(ro,4),REAL(zo,4),
     $     REAL(psio,4),REAL(betap1,4),REAL(betap2,4),REAL(betap3,4),
     $     REAL(betat,4),REAL(betan,4),REAL(bt0,4),REAL(q0,4),
     $     REAL(qmin,4),REAL(qmax,4),REAL(qa,4),REAL(crnt,4),
     $     REAL(plasma1,4),REAL(vacuum1,4),REAL(total1),REAL(q95)
      CALL bin_close(sum_unit)
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
         CALL cspline_dealloc(fmats)
         CALL cspline_dealloc(gmats)
         CALL cspline_dealloc(kmats)
         DO ising=1,msing
            DEALLOCATE(sing(ising)%vmat)
            DEALLOCATE(sing(ising)%mmat)
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
