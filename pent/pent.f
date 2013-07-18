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

      LOGICAL :: lar_flag,gen_flag,o1native,fmnl_flag,diag_flag
      INTEGER :: nl,maskpsi,status,i
      REAL(r8) :: ptfac,ximag,xmax,wefac,wdfac
      REAL(r8), DIMENSION(10) :: psiout
      CHARACTER(16)  :: collision 
      CHARACTER(128) :: kinetic_file,o1fun_file,o1mn_file,buffer
      
      NAMELIST/pent_input/kinetic_file,o1fun_file,o1mn_file,o1native,
     $     imass,icharge,zimp,zmass,nl,
     $     collision
      NAMELIST/pent_control/maskpsi,wefac,wdfac,
     $     neorot_flag,divxprp_flag,kolim_flag,xlsode_flag,
     $     ximag,xmax,nx,nt,nlmda,ptfac
      NAMELIST/pent_output/kinetic_flag,pitch_flag,bounce_flag,
     $     energy_flag,passing_flag,lar_flag,gen_flag,offset_flag,psiout
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

      kolim_flag = .FALSE.
      neorot_flag = .FALSE.
      divxprp_flag = .TRUE.
      xlsode_flag = .FALSE.
      wefac = 1
      wdfac = 1
      ptfac = 1e-3
      ximag = 1e-6
      xmax = 36
      nx = 200
      nt = 101
      nlmda = 101

      kinetic_flag=.FALSE.
      bounce_flag=.FALSE.
      pitch_flag=.FALSE.
      energy_flag=.FALSE.
      passing_flag=.FALSE.
      lar_flag = .FALSE.
      gen_flag = .TRUE.
      psiout = 0
      psiout(:) = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
      
      fmnl_flag = .FALSE.
      diag_flag = .FALSE.

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
      CALL idcon_vacuum

c-----------------------------------------------------------------------
c     Set some global standards.
c-----------------------------------------------------------------------
      ALLOCATE(xnorm(0:nx),tnorm(0:nt),lnorm(0:nlmda))
      xnorm =(/(i*1.0/nx,i=0,nx)/)
      tnorm =(/(i*1.0/nx,i=0,nt)/)
      lnorm =(/(i*1.0/nx,i=0,nlmda)/)
c-----------------------------------------------------------------------
c     Computer Non-ambipolar transport toroidal torque.
c-----------------------------------------------------------------------
      CALL pentio_read_profile(kinetic_file)
      IF(.NOT. diag_flag)THEN
         IF(gen_flag)THEN
            ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
            CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
            CALL profile_gen(bpar,xprp,nl,maskpsi,ptfac,ximag,xmax,
     $           wefac,wdfac,collision,psiout)
            DEALLOCATE(bpar,xprp)
         ENDIF
         IF(lar_flag)THEN         
            ALLOCATE(bpar(mstep,mpert),xprp(mstep,mpert))
            CALL pentio_read_order1_mns(bpar,xprp,o1mn_file,o1native)
            CALL pentio_read_special_functions
            CALL profile_lar(bpar,xprp,nl,maskpsi,ximag,xmax,wefac,
     $           wdfac,collision,psiout)
            DEALLOCATE(bpar,xprp)
         ENDIF
         CALL pentg_dealloc(gen_flag,lar_flag)
      ENDIF

c-----------------------------------------------------------------------
c     Administrative setup/diagnostics/debugging.
c-----------------------------------------------------------------------
      IF(fmnl_flag) CALL pentg_fmnl
      IF(diag_flag)THEN
c         i = diag_outer(1)
         CALL diag_unpack
         CALL diag_grid
c         CALL diag_lsode_example
c         CALL diag_lsode_custom
         ALLOCATE(bpar(mstep,0:mthsurf),xprp(mstep,0:mthsurf))
         CALL pentio_read_order1_funs(bpar,xprp,o1fun_file,o1native)
         CALL diag_lsode_gen(bpar,xprp,nl,maskpsi,ptfac,ximag,xmax,
     $           wefac,wdfac,collision,psiout)
         DEALLOCATE(bpar,xprp)
      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL pentio_stop("Normal termination.")
      
      END PROGRAM pent_main
