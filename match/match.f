c-----------------------------------------------------------------------
c     file match.f.
c     main program for matched asymptotic modes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. match.
c-----------------------------------------------------------------------
c     subprogram 1. match.
c     main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM match_main
      USE ideal_mod
      USE resist_mod
      IMPLICIT NONE

      LOGICAL :: scan_flag=.FALSE.,root_flag=.FALSE.,
     $     contour_flag=.FALSE.
      CHARACTER(128) :: filename="euler.bin"
      INTEGER :: ising,ns
      REAL(r8) :: big,small
      COMPLEX(r8) :: root

      NAMELIST/match_input/filename,ideal_flag,res_flag,ff0,ff1,
     $     big,small,ns,root,scan_flag,root_flag,contour_flag,full,
     $     z_level_type,phi_level,transform_flag,poly_form_type,
     $     poly_test_type,scan0,scan1,nscan,delta0,delta1,condense,
     $     ripple_flag,mripple
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(4x,"is",4x,"re f1",6x,"im f1",6x,"re f2",6x,"im f2",
     $     6x,"re f3",6x,"im f3",3x,"re deltap1",1x,"im deltap1",
     $     1x,"re deltap2",1x,"im deltap2"
     $     //(i6,1p,10e11.3/))
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"match.in","OLD")
      READ(in_unit,NML=match_input)
      CALL ascii_close(in_unit)
      phi_level=phi_level*pi/180
c-----------------------------------------------------------------------
c     read data file.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"match.out","UNKNOWN")
      CALL timer(0,out_unit)
      CALL ideal_read(filename)
c-----------------------------------------------------------------------
c     reconstruct most unstable ideal solution.
c-----------------------------------------------------------------------
      IF(ideal_flag)THEN
         CALL ideal_transform
         CALL ideal_build
         CALL ideal_write
         CALL ideal_chord
         IF(contour_flag)CALL ideal_contour
      ENDIF
c-----------------------------------------------------------------------
c     study resistive modes.
c-----------------------------------------------------------------------
      IF(res_flag)THEN
         CALL matrix_run
         DO jsing=1,msing
            CALL resist_deltac_find
         ENDDO
         CALL resist_write
         IF(scan_flag)CALL resist_nyquist(big,small,ns)
         IF(root_flag)CALL resist_root_find(root)
      ENDIF
c-----------------------------------------------------------------------
c     write sum2.dat.
c-----------------------------------------------------------------------
      IF(res_flag)THEN
         CALL bin_open(sum_unit,"sum2.dat","UNKNOWN","REWIND","none")
         WRITE(sum_unit)msing,mterm
         WRITE(sum_unit)
     $        (REAL(singtype(ising)%psifac,4),ising=1,msing),
     $        (REAL(singtype(ising)%q,4),ising=1,msing),
     $        (REAL(singtype(ising)%restype%di,4),ising=1,msing),
     $        REAL(deltap1,4),REAL(deltap2,4),REAL(ff,4)
         CALL bin_close(sum_unit)
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      CALL match_dealloc
      CALL program_stop("Normal termination.")
      END PROGRAM match_main
