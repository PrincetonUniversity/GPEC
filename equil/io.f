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
      INTEGER :: sing_unit=7

      INTEGER :: dump_unit=9
      INTEGER :: sum_unit=10

      INTEGER :: out_2d_unit=11
      INTEGER :: bin_2d_unit=12

      INTEGER :: lar_out_unit=13
      INTEGER :: lar_bin_unit=14

      INTEGER :: out1_unit=21
      INTEGER :: out2_unit=22
      INTEGER :: out3_unit=23
      INTEGER :: gal_out_unit=24

      INTEGER :: bin1_unit=31
      INTEGER :: bin2_unit=32
      INTEGER :: bin3_unit=33
      INTEGER :: gal_bin_unit=34

      ! rmatch specific units
      INTEGER :: match_unit=41
      INTEGER :: vsol_out_unit=42
      INTEGER :: vsol_bin_unit=43
      INTEGER :: jmat_unit=52
      INTEGER :: coefs_unit=54
      INTEGER :: split_unit=53
      INTEGER :: deltac_out_unit=57
      INTEGER :: deltac_bin_unit=58
      INTEGER :: inpso_out_unit=59
      INTEGER :: inpso_bin_unit=60
      INTEGER :: array_unit=61
      INTEGER :: grid_out_unit=62

      INTEGER :: debug_unit=99
      END MODULE io_mod
