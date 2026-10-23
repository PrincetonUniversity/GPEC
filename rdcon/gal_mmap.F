c-----------------------------------------------------------------------
c     file gal_mmap.F.
c     lets glibc return the freed Galerkin element matrices to the
c     system.
c-----------------------------------------------------------------------
c     gal_assemble_mat frees each element matrix once it is in the
c     band. glibc raises its mmap threshold as large blocks are freed,
c     so later blocks of that size come from the heap, and their memory
c     is kept after they are freed. fixing the threshold at 128 kB,
c     below the size of an element matrix, keeps each element matrix in
c     its own mapping, which is returned to the system when it is
c     freed. this is the same as setting MALLOC_MMAP_THRESHOLD_=131072.
c     the build defines HAVE_MALLOPT on Linux. macOS has no mallopt and
c     already returns such blocks, so there this does nothing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_mmap_threshold
#ifdef HAVE_MALLOPT
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int
      IMPLICIT NONE

      INTERFACE
         FUNCTION mallopt(param,value) BIND(C,name="mallopt")
         IMPORT :: c_int
         INTEGER(c_int), VALUE :: param,value
         INTEGER(c_int) :: mallopt
         END FUNCTION mallopt
      END INTERFACE

      INTEGER(c_int), PARAMETER :: m_mmap_threshold=-3
      INTEGER(c_int) :: info
c-----------------------------------------------------------------------
c     set the threshold.
c-----------------------------------------------------------------------
      info=mallopt(m_mmap_threshold,131072_c_int)
#endif
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_mmap_threshold
