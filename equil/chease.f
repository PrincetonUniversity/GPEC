c-----------------------------------------------------------------------
c     subprogram 4. read_eq_chease.
c     reads data from chease.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_eq_chease
      USE inverse_mod

      INTEGER :: ntnova,npsi1,isym,ma,mtau
      REAL(r8), DIMENSION(5) :: axx
      REAL(r8), DIMENSION(:), POINTER :: zcpr,zcppr,zq,zdq,ztmf,
     $     ztp,zfb,zfbp,zpsi,zpsim
      REAL(r8), DIMENSION(:,:), POINTER :: zrcp,zzcp,zjacm,zjac
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,TRIM(eq_filename),"OLD","REWIND")
      READ(in_unit)ntnova,npsi1,isym
      READ(in_unit)axx
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(zcpr(npsi),zcppr(npsi1),zq(npsi1),zdq(npsi1),ztmf(npsi1),
     $     ztp(npsi1),zfb(npsi1),zfbp(npsi1),zpsi(npsi1),zpsim(npsi1))
      ALLOCATE(zrcp(ntnova+3,npsi1),zzcp(ntnova+3,npsi1),
     $     zjacm(ntnova+3,npsi1),zjac(ntnova+3,npsi1))
c-----------------------------------------------------------------------
c     read arrays and close file.
c-----------------------------------------------------------------------
      READ(in_unit)zcpr
      READ(in_unit)zcppr
      READ(in_unit)zq
      READ(in_unit)zdq
      READ(in_unit)ztmf
      READ(in_unit)ztp
      READ(in_unit)zfb
      READ(in_unit)zfbp
      READ(in_unit)zpsi
      READ(in_unit)zpsim
      READ(in_unit)zrcp
      READ(in_unit)zzcp
      READ(in_unit)zjacm
      READ(in_unit)zjac
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     allocate native arrays.
c-----------------------------------------------------------------------
      mtau=ntnova+2
      ma=npsi1-1
      CALL spline_alloc(sq_in,ma,4)
      CALL bicube_alloc(rz_in,ma,mtau,2)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      sq_in%xs=zpsi(2:)
      sq_in%fs(:,1)=ztmf(2:)
      sq_in%fs(:,2)=zcpr(2:)
      sq_in%fs(:,3)=zq(2:)
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(zcpr,zcppr,zq,zdq,ztmf,ztp,zfb,zfbp,zpsi,zpsim)
      DEALLOCATE(zrcp,zzcp,zjacm,zjac)
      CALL program_stop("Termination by read_eq_chease.")
      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_eq_chease
