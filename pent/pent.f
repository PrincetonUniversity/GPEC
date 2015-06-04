c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     main program 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     pent_main.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM pent_main
      USE profile_mod
      USE diagnostic_mod
      USE idcon_mod

      IMPLICIT NONE

      LOGICAL :: gen_flag,rlar_flag,clar_flag,cgl_flag,o1native,
     $     fmnl_flag,diag_flag
      INTEGER :: maskpsi,status,i
      REAL(r8) :: ptfac,ximag,xmax,wefac,wdfac,wpfac,diag_psi,gafac
      REAL(r8), DIMENSION(10) :: psiout
      CHARACTER(128) :: kinetic_file,o1fun_file,o1mn_file,buffer
      
      NAMELIST/pent_input/kinetic_file,o1fun_file,o1mn_file,o1native,
     $     idconfile,imass,icharge,zimp,zmass,nl,
     $     collision
      NAMELIST/pent_control/maskpsi,wefac,wdfac,wpfac,
     $     neorot_flag,divxprp_flag,mdc2_flag,kolim_flag,xlsode_flag,
     $     lmdalsode_flag,lsode_rtol,lsode_atol,ximag,xmax,nx,nt,nlmda,
     $     ptfac,treg_flag,gafac
      NAMELIST/pent_output/kinetic_flag,pitch_flag,bounce_flag,
     $     energy_flag,passing_flag,rlar_flag,gen_flag,clar_flag,
     $     cgl_flag,offset_flag,psiout
      NAMELIST/pent_admin/fmnl_flag,diag_flag
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      jsurf_in=0
      tmag_in=1
      jac_in=""
      ieqfile="psi_in.bin"
      idconfile="euler.bin"
      ivacuumfile="vacuum.bin"
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


      psixy = 0

      kinetic_file="kin.dat"
      o1fun_file ="ipec_order1_fun_n1.bin"
      o1mn_file  ="ipec_order1_n1.bin"
      o1native = .TRUE.
      imass=2
      icharge=1
      zimp=6
      zmass=12
      nl=0
      maskpsi=1
      collision = "harmonic"

      mdc2_flag = .FALSE.
      treg_flag = .FALSE.
      kolim_flag = .FALSE.
      neorot_flag = .FALSE.
      divxprp_flag = .TRUE.
      xlsode_flag = .FALSE.
      lmdalsode_flag = .FALSE.
      lsode_rtol = 1e-9
      lsode_atol = 1e-12
      wefac = 1
      wdfac = 1
      wpfac = 1
      gafac = 0
      ptfac = 1e-3
      ximag = 0
      xmax = 36
      nx = 200
      nt = 101
      nlmda = 101

      kinetic_flag=.FALSE.
      bounce_flag=.FALSE.
      pitch_flag=.FALSE.
      energy_flag=.FALSE.
      passing_flag=.FALSE.
      rlar_flag = .FALSE.
      clar_flag = .FALSE.
      cgl_flag = .FALSE.
      gen_flag = .TRUE.
      psiout(:) = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
      
      fmnl_flag = .FALSE.
      diag_flag = .FALSE.

c-----------------------------------------------------------------------
c     read pent.in.
c-----------------------------------------------------------------------
      PRINT *,"PENT START=> v1.00"
c-----------------------------------------------------------------------
c     read pent.in.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"pent.in","OLD")
      READ(in_unit,NML=pent_input)
      READ(in_unit,NML=pent_control)
      READ(in_unit,NML=pent_output)
      READ(in_unit,NML=pent_admin)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     check time.
c-----------------------------------------------------------------------
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
c-----------------------------------------------------------------------
c     prepare ideal solutions.
c-----------------------------------------------------------------------
      CALL idcon_read(psixy)
      CALL idcon_transform
c-----------------------------------------------------------------------
c     reconstruct metric tensors.
c-----------------------------------------------------------------------
      CALL idcon_metric
c-----------------------------------------------------------------------
c     read vacuum data.
c-----------------------------------------------------------------------
c      CALL idcon_vacuum


c-----------------------------------------------------------------------
c     Some intial consistency checks.
c-----------------------------------------------------------------------
      IF(collision=="zero" .AND. ximag==0) CALL pentio_stop("Must use "
     $     //"complex energy contour (xiamg/=0) when collisionless.")
c-----------------------------------------------------------------------
c     Set some global standards.
c-----------------------------------------------------------------------
      ALLOCATE(xnorm(0:nx),tnorm(0:nt),lnorm(0:nlmda))
      xnorm =(/(i*1.0/nx,i=0,nx)/)
      tnorm =(/(i*1.0/nx,i=0,nt)/)
      lnorm =(/(i*1.0/nx,i=0,nlmda)/)
c-----------------------------------------------------------------------
c     Initial io.
c-----------------------------------------------------------------------
      CALL pentio_cleardir
      CALL pentio_headers(nl,ximag,gen_flag,rlar_flag,clar_flag)
c-----------------------------------------------------------------------
c     Computer Non-ambipolar transport toroidal torque.
c-----------------------------------------------------------------------
      IF(.NOT. diag_flag)THEN
         IF(gen_flag)THEN
            ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
            CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
            CALL pentio_read_profile(kinetic_file,wefac,wpfac)
            CALL profile_gen(bpar,xprp,nl,maskpsi,ptfac,ximag,xmax,
     $           wefac,wdfac,wpfac,gafac,collision,psiout,.FALSE.)
            CALL pentg_dealloc
         ENDIF
         IF(clar_flag)THEN
            ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
            CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
            CALL pentio_read_profile(kinetic_file,wefac,wpfac)
            CALL pentio_read_special_functions
            CALL profile_gen(bpar,xprp,nl,maskpsi,ptfac,ximag,xmax,
     $           wefac,wdfac,wpfac,gafac,collision,psiout,.TRUE.)
            CALL pentg_dealloc
         ENDIF
         IF(rlar_flag)THEN         
            ALLOCATE(bpar(mstep,mpert),xprp(mstep,mpert))
            CALL pentio_read_order1_mns(bpar,xprp,o1mn_file,o1native)
            CALL pentio_read_profile(kinetic_file,wefac,wpfac)
            CALL pentio_read_special_functions
            CALL profile_rlar(bpar,xprp,nl,maskpsi,ximag,xmax,wefac,
     $           wdfac,wpfac,gafac,collision,psiout)
            CALL pentg_dealloc
         ENDIF
         IF(cgl_flag)THEN         
            ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
            CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
            CALL pentio_read_profile(kinetic_file,wefac,wpfac)
            CALL profile_cgl(bpar,xprp,maskpsi)
            CALL pentg_dealloc
         ENDIF
      ENDIF

c-----------------------------------------------------------------------
c     Administrative setup/diagnostics/debugging.
c-----------------------------------------------------------------------
      IF(fmnl_flag) CALL pentg_fmnl
      IF(diag_flag)THEN
         CALL diag_complexpow
         CALL diag_elliptics
c         i = diag_outer(1)
         CALL diag_grid
         i = 1
         call diag_sub1(i)
c         CALL diag_lsode_example
         CALL diag_lsode_custom
         ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
         CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
         diag_psi = 0.95
         CALL diag_singlesurf(bpar,xprp,nl,diag_psi,ptfac,ximag,xmax,
     $           wefac,wdfac,collision,psiout)
         DEALLOCATE(bpar,xprp)
      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL pentio_stop("Normal termination.")
      
      END PROGRAM pent_main
