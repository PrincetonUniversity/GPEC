c-----------------------------------------------------------------------
c     file equil_.f.
c     reads ascii input for equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. equil_mod.
c     1. equil_read.
c-----------------------------------------------------------------------
c     subprogram 0. equil_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE equil_mod
      USE read_eq_mod
      USE lar_mod
      USE gsec_mod
      USE sol_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. equil_read.
c     reads input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE equil_read(unit, op_psihigh)

      LOGICAL :: file_stat
      INTEGER, INTENT(IN) :: unit
      REAL(r8), OPTIONAL, INTENT(IN) :: op_psihigh

      NAMELIST/equil_control/eq_filename,eq_type,grid_type,mpsi,mtheta,
     $     newq0,psihigh,psilow,input_only,jac_type,power_bp,power_r,
     $     power_b,jac_method,convert_type,power_flag,ns1,
     $     sp_pfac,sp_nx,sp_dx1,sp_dx2,use_galgrid,etol,
     $     use_classic_splines
      NAMELIST/equil_output/bin_2d,bin_eq_1d,bin_eq_2d,out_2d,out_eq_1d,
     $     out_eq_2d,bin_fl,out_fl,interp,gse_flag,dump_flag,verbose
c-----------------------------------------------------------------------
c     read input data.
c-----------------------------------------------------------------------
      INQUIRE(FILE="equil.in",EXIST=file_stat)
      IF(.NOT.file_stat)CALL program_stop
     $     ("Can't open input file equil.in")
      CALL ascii_open(in_unit,"equil.in","OLD")
      IF (eq_type.NE."repack_eq") THEN
         READ(UNIT=in_unit,NML=equil_control)
         READ(UNIT=in_unit,NML=equil_output)
      ENDIF
      CALL ascii_close(in_unit)
      IF(PRESENT(op_psihigh))THEN
         psihigh = op_psihigh
         IF(verbose) WRITE(*,*) "Reforming equilibrium with new psihigh"
      ENDIF
      IF(verbose) WRITE(*,'(1x,a,es10.3)')"psihigh =",psihigh
      psihigh=MIN(psihigh,1._r8)
c-----------------------------------------------------------------------
c     define Jacobian.
c-----------------------------------------------------------------------
      SELECT CASE(jac_type)
      CASE("hamada")
         power_b=0
         power_bp=0
         power_r=0
      CASE("pest")
         power_b=0
         power_bp=0
         power_r=2
      CASE("equal_arc")
         power_b=0
         power_bp=1
         power_r=0
      CASE("boozer")
         power_b=2
         power_bp=0
         power_r=0
      CASE("park")
         power_b=1
         power_bp=0
         power_r=0
      CASE("other")
      CASE DEFAULT
         CALL program_stop
     $        ("Cannot recognize jac_type = "//TRIM(jac_type))
      END SELECT
c-----------------------------------------------------------------------
c     write equilibrium and jacobian data.
c-----------------------------------------------------------------------
      IF(.NOT. PRESENT(op_psihigh))THEN  ! only need this once
        IF(verbose) WRITE(*,'(1x,4a/1x,2a,3(a,i1))')
     $       "Equilibrium: ",TRIM(eq_filename),", Type: ",TRIM(eq_type),
     $       "Jac_type = ",TRIM(jac_type),", power_bp = ",power_bp,
     $       ", power_b = ",power_b,", power_r = ",power_r
        WRITE(unit,'(1x,4a//1x,2a,3(a,i1))')
     $       "Equilibrium: ",TRIM(eq_filename),", Type: ",TRIM(eq_type),
     $       "Jac_type = ",TRIM(jac_type),", power_bp = ",power_bp,
     $       ", power_b = ",power_b,", power_r = ",power_r
      ENDIF
c-----------------------------------------------------------------------
c     read equilibrium data and diagnose.
c-----------------------------------------------------------------------
      SELECT CASE(TRIM(eq_type))
      CASE("fluxgrid")
         CALL read_eq_fluxgrid
      CASE("miller")
         CALL read_eq_miller
      CASE("miller4")
         CALL read_eq_miller4
      CASE("galkin")
         CALL read_eq_galkin
      CASE("chease")
         CALL read_eq_chease
      CASE("chease2")
         CALL read_eq_chease2
      CASE("chease3")
         CALL read_eq_chease3
      CASE("chease4")
         CALL read_eq_chease4
      CASE("chum")
         CALL read_eq_chum
      CASE("efit")
         CALL read_eq_efit
      CASE("rsteq")
         CALL read_eq_rsteq
      CASE("ldp_d")
         CALL read_eq_ldp_d
      CASE("ldp_i")
         CALL read_eq_ldp_i
      CASE("jsolver")
         CALL read_eq_jsolver
      CASE("lez")
         CALL read_eq_lez
      CASE("sontag")
         CALL read_eq_sontag
      CASE("tokamac")
         CALL read_eq_tokamac
      CASE("lar")
         CALL lar_run
      CASE("gsec")
         CALL gsec_run
      CASE("soloviev")
         CALL sol_run
      CASE("sol_galkin")
         CALL sol_galkin
      CASE("transp")
         CALL read_eq_transp
      CASE("popov1")
         CALL read_eq_popov1
      CASE("popov2")
         CALL read_eq_popov2
      CASE("rtaylor")
         CALL read_eq_rtaylor
      CASE("wdn")
         CALL read_eq_wdn
      CASE("dump")
         CALL read_eq_dump
      CASE("vmec")
         CALL read_eq_vmec
      CASE("t7_old")
         CALL read_eq_t7_old
      CASE("t7")
         CALL read_eq_t7
      CASE("marklin_direct")
         CALL read_eq_marklin_direct
      CASE("marklin_inverse")
         CALL read_eq_marklin_inverse
      CASE("hansen_inverse")
         CALL read_eq_hansen_inverse
      CASE("pfrc")
         CALL read_eq_pfrc
      CASE DEFAULT
         CALL program_stop("Cannot recognize eq_type "//TRIM(eq_type))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_read
      END MODULE equil_mod
