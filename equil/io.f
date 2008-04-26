c-----------------------------------------------------------------------
c     file io.f.
c     input and output unit declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
      MODULE io_mod
      IMPLICIT NONE

      INTEGER :: in_unit=1
      INTEGER :: out_unit=2
      INTEGER :: bin_unit=3
      INTEGER :: equil_unit=4
      INTEGER :: term_unit=6

      INTEGER :: dump_unit=9
      INTEGER :: sum_unit=10

      INTEGER :: out_2d_unit=11
      INTEGER :: bin_2d_unit=12

      INTEGER :: lar_out_unit=13
      INTEGER :: lar_bin_unit=14

      INTEGER :: debug_unit=99

      END MODULE io_mod
