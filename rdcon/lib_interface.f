c-----------------------------------------------------------------------
c     file lib_interface.f.
c     interfact between dcon and lib
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. lib_interface_mod.
c     1. lib_interface_init.
c     2. lib_interface_set_equil.
c-----------------------------------------------------------------------
c     subprogram 0. lib_interface_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE lib_interface_mod
      USE utils_mod
      USE global_mod
      USE local_mod
      USE inverse_mod
      IMPLICIT NONE

      TYPE transpeq         
      INTEGER :: nsm1, ntm1     !R_minor, poloidal angle dimensions respectively
      INTEGER :: nr, nz         !Cylindrical dims for inverse equilibrium
      REAL(r8), DIMENSION(:),   POINTER :: psis, qs, fs, ps  !Flux, safety factor, R*B_T, pressure (0:nsm1)
      REAL(r8), DIMENSION(:,:), POINTER :: rg, zg            !Cylindrical (R,Z) coords (0:nsm1,0:ntm1)
      REAL(r8), DIMENSION(:,:), POINTER :: psirz             !Inverse equilibrium psi(R,z) (1:nr,1:nz)
      END TYPE transpeq

      TYPE :: lib_interface_data
      REAL(r8) :: time
      INTEGER :: ntheta,nsurf
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psis,ps,qs,fs
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rg,zg
      END TYPE lib_interface_data

      TYPE(transpeq) :: teq
      TYPE(lib_interface_data) :: lib_input

      PUBLIC :: teq

      LOGICAL :: run_lib_interface
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. lib_interface_init.
c     init interface for external call.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lib_interface_init()
         run_lib_interface=.TRUE.
      END SUBROUTINE lib_interface_init
c-----------------------------------------------------------------------
c     subprogram 2. lib_interface_input.
c     set equilibrium input to prepare DCON run.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lib_interface_input()
         CALL lib_interface_input_memory()
C         CALL lib_interface_input_reorder()
C         CALL lib_interface_input_1()
C         CALL lib_interface_input_chease()
      END SUBROUTINE lib_interface_input
c-----------------------------------------------------------------------
c     subprogram 3. lib_interface_set_equil.
c     set equilibrium spline to prepare DCON run.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lib_interface_set_equil()
         CALL lib_interface_set_equil_transp()
C         CALL lib_interface_set_equil_chease()
      END SUBROUTINE lib_interface_set_equil


      SUBROUTINE lib_interface_input_memory()
      INTEGER :: isurf,itheta,ntheta,nsurf,nr,nz,idx,idx1,dir
      REAL(r8) :: time,maxr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq,jaceq,g12l,rg,zg
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------

      lib_input%time=0
      lib_input%ntheta=teq%nsm1
      lib_input%nsurf=teq%ntm1
      ntheta=lib_input%ntheta
      nsurf=lib_input%nsurf
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lib_input%psis(nsurf),lib_input%ps(nsurf),
     $         lib_input%qs(nsurf),lib_input%fs(nsurf))
      ALLOCATE(lib_input%rg(nsurf,ntheta),
     $         lib_input%zg(nsurf,ntheta))
      ALLOCATE(rg(nsurf,ntheta),zg(nsurf,ntheta))
c-----------------------------------------------------------------------
c     assign profile
c-----------------------------------------------------------------------
      lib_input%psis=teq%psis
      lib_input%ps=teq%ps
      lib_input%qs=teq%qs
      lib_input%fs=teq%fs
      rg=teq%rg
      zg=teq%zg

      maxr=rg(nsurf,1)
      idx=1
      DO itheta=2,ntheta
         IF (rg(nsurf,itheta)>maxr) THEN 
            maxr=rg(nsurf,itheta)
            idx=itheta
         ENDIF
      ENDDO
      idx1=MOD(idx,ntheta)+1
      IF (zg(nsurf,idx1)>zg(nsurf,idx)) THEN
         dir=1
      ELSE
         dir=-1
      ENDIF
      DO itheta=0,ntheta-1
         idx1=mod(idx-1+itheta*dir+ntheta,ntheta)+1
         lib_input%rg(:,itheta+1)=rg(:,idx1)
         lib_input%zg(:,itheta+1)=zg(:,idx1)
      ENDDO
      DEALLOCATE(rg,zg)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_input_memory

      SUBROUTINE lib_interface_input_reorder()
      INTEGER :: isurf,itheta,ntheta,nsurf,nr,nz,idx,idx1,dir
      REAL(r8) :: time,maxr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq,jaceq,g12l,rg,zg
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
10    FORMAT(1p,3e24.16)
      CALL ascii_open(in_unit,"dconeq_tftr.dat","OLD")
c      OPEN(in_unit,"dconeq.dat","OLD")
      READ(in_unit,'(1e14.6)') time
      READ(in_unit,'(2I7)') ntheta,nsurf
      READ(in_unit,'(2I7)') nr,nz
      lib_input%time=time
      lib_input%ntheta=ntheta
      lib_input%nsurf=nsurf
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lib_input%psis(nsurf),lib_input%ps(nsurf),
     $         lib_input%qs(nsurf),lib_input%fs(nsurf))
      ALLOCATE(lib_input%rg(nsurf,ntheta),
     $         lib_input%zg(nsurf,ntheta))
      ALLOCATE(rg(nsurf,ntheta),zg(nsurf,ntheta))
c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      READ(in_unit,10) (lib_input%psis(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%ps(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%qs(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%fs(isurf),isurf=1,nsurf)
      READ(in_unit,10) 
     $     (rg(1:nsurf,itheta),itheta=1,ntheta)
      READ(in_unit,10) 
     $     (zg(1:nsurf,itheta),itheta=1,ntheta)

      CALL ascii_close(in_unit)
      maxr=rg(nsurf,1)
      idx=1
      DO itheta=2,ntheta
         IF (rg(nsurf,itheta)>maxr) THEN 
            maxr=rg(nsurf,itheta)
            idx=itheta
         ENDIF
      ENDDO
      idx1=MOD(idx,ntheta)+1
      IF (zg(nsurf,idx1)>zg(nsurf,idx)) THEN
         dir=1
      ELSE
         dir=-1
      ENDIF
      DO itheta=0,ntheta-1
         idx1=mod(idx-1+itheta*dir+ntheta,ntheta)+1
         lib_input%rg(:,itheta+1)=rg(:,idx1)
         lib_input%zg(:,itheta+1)=zg(:,idx1)
      ENDDO
      DEALLOCATE(rg,zg)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_input_reorder

      SUBROUTINE lib_interface_input_reorder_good()
      INTEGER :: isurf,itheta,ntheta,nsurf,nr,nz,idx,idx1,dir
      REAL(r8) :: time,maxr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq,jaceq,g12l,rg,zg
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
10    FORMAT(1p,4e20.12)
      CALL ascii_open(in_unit,"dconeq_inv.dat","OLD")
c      OPEN(in_unit,"dconeq.dat","OLD")
      READ(in_unit,'(1e14.6)') time
      READ(in_unit,'(I7)') ntheta
      READ(in_unit,'(I7)') nsurf
      READ(in_unit,'(2I7)') nr,nz
      lib_input%time=time
      lib_input%ntheta=ntheta
      lib_input%nsurf=nsurf
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lib_input%psis(nsurf),lib_input%ps(nsurf),
     $         lib_input%qs(nsurf),lib_input%fs(nsurf))
      ALLOCATE(lib_input%rg(nsurf,ntheta),
     $         lib_input%zg(nsurf,ntheta))
      ALLOCATE(rg(nsurf,ntheta),zg(nsurf,ntheta))
c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      READ(in_unit,10) (lib_input%psis(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%ps(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%qs(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%fs(isurf),isurf=1,nsurf)
      READ(in_unit,10) 
     $     (rg(1:nsurf,itheta),itheta=1,ntheta)
      READ(in_unit,10) 
     $     (zg(1:nsurf,itheta),itheta=1,ntheta)

      CALL ascii_close(in_unit)
      maxr=rg(nsurf,1)
      idx=1
      DO itheta=2,ntheta
         IF (rg(nsurf,itheta)>maxr) THEN 
            maxr=rg(nsurf,itheta)
            idx=itheta
         ENDIF
      ENDDO
      idx1=MOD(idx,ntheta)+1
      IF (zg(nsurf,idx1)>zg(nsurf,idx)) THEN
         dir=1
      ELSE
         dir=-1
      ENDIF
      DO itheta=0,ntheta-1
         idx1=mod(idx-1+itheta*dir+ntheta,ntheta)+1
         lib_input%rg(:,itheta+1)=rg(:,idx1)
         lib_input%zg(:,itheta+1)=zg(:,idx1)
      ENDDO
      DEALLOCATE(rg,zg)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_input_reorder_good




c-----------------------------------------------------------------------
c     subprogram 1. lib_interface_input.
c     read lib input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lib_interface_input_1()
      INTEGER :: isurf,itheta,ntheta,nsurf,nr,nz
      REAL(r8) :: time
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq,jaceq,g12l
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
10    FORMAT(1p,4e20.12)
      CALL ascii_open(in_unit,"dconeq_inv.dat","OLD")
c      OPEN(in_unit,"dconeq.dat","OLD")
      READ(in_unit,'(1e14.6)') time
      READ(in_unit,'(I7)') ntheta
      READ(in_unit,'(I7)') nsurf
      READ(in_unit,'(2I7)') nr,nz
      lib_input%time=time
      lib_input%ntheta=ntheta
      lib_input%nsurf=nsurf
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lib_input%psis(nsurf),lib_input%ps(nsurf),
     $         lib_input%qs(nsurf),lib_input%fs(nsurf))
      ALLOCATE(lib_input%rg(nsurf,ntheta),
     $         lib_input%zg(nsurf,ntheta))
c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      READ(in_unit,10) (lib_input%psis(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%ps(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%qs(isurf),isurf=1,nsurf)
      READ(in_unit,10) (lib_input%fs(isurf),isurf=1,nsurf)
      READ(in_unit,10) 
     $     (lib_input%rg(1:nsurf,itheta),itheta=1,ntheta)
      READ(in_unit,10) 
     $     (lib_input%zg(1:nsurf,itheta),itheta=1,ntheta)


      CALL ascii_close(in_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_input_1
c-----------------------------------------------------------------------
c     subprogram 2. lib_interface_set_equil.
c     
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lib_interface_set_equil_transp()
      INTEGER :: isurf,itheta
c-----------------------------------------------------------------------
c     copy variables.
c-----------------------------------------------------------------------
c      ro=SUM(lib_input%rg(1,1))/(lib_input%ntheta+1)
c      zo=SUM(lib_input%zg(1,1))/(lib_input%ntheta+1)
      ro=lib_input%rg(1,1)
      zo=lib_input%zg(1,1)
      psio=abs(lib_input%psis(lib_input%nsurf))
      lib_input%psis=abs(lib_input%psis/psio)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,lib_input%nsurf-1,4)
      sq_in%xs=lib_input%psis
      sq_in%fs(:,1)=lib_input%fs
      sq_in%fs(:,2)=mu0*(lib_input%ps-lib_input%ps(lib_input%nsurf))
      sq_in%fs(:,3)=lib_input%qs
c      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      CALL bicube_alloc(rz_in,lib_input%nsurf-1,lib_input%ntheta-1,2)
      DO isurf=1,lib_input%nsurf
         DO itheta=1,lib_input%ntheta
            rz_in%fs(isurf-1,itheta-1,1)=lib_input%rg(isurf,itheta)
            rz_in%fs(isurf-1,itheta-1,2)=lib_input%zg(isurf,itheta)
         ENDDO
      ENDDO
20    FORMAT(81e14.6)
      CALL ascii_open(in_unit,"rg.dat","REPLACE")
      DO isurf=1,lib_input%nsurf
         WRITE(in_unit,20) lib_input%rg(isurf,:)
      ENDDO
      CALL ascii_close(in_unit)
      CALL ascii_open(in_unit,"zg.dat","REPLACE")
      DO isurf=1,lib_input%nsurf
         WRITE(in_unit,20) lib_input%zg(isurf,:)
      ENDDO
      CALL ascii_close(in_unit)

c      CALL bicube_fit(rz_in,"extrap","periodic")
c      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_set_equil_transp




      SUBROUTINE lib_interface_input_chease

      INTEGER :: npsi,ipsi,nchi
      REAL(r8) :: b0exp,r0exp,rmag,zmag,tpsio
      REAL(r8), DIMENSION(:), ALLOCATABLE :: norpsi,psi,t,mu0p,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: req,zeq,jaceq,g12l
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="OUTRDCON",STATUS="OLD",
     $     FORM="UNFORMATTED")
      READ(in_unit)npsi,nchi,r0exp,b0exp,rmag,zmag,tpsio

      lib_input%time=0
      lib_input%ntheta=nchi
      lib_input%nsurf=npsi
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(norpsi(npsi),psi(npsi),t(npsi),mu0p(npsi),q(npsi))
      ALLOCATE(req(npsi,nchi),zeq(npsi,nchi),
     $        jaceq(npsi,nchi),g12l(npsi,nchi))


c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lib_input%psis(npsi),lib_input%ps(npsi),
     $         lib_input%qs(npsi),lib_input%fs(npsi))
      ALLOCATE(lib_input%rg(npsi,nchi),
     $         lib_input%zg(npsi,nchi))


c-----------------------------------------------------------------------
c     read local arrays and close file.
c-----------------------------------------------------------------------
      DO ipsi=1,npsi
         READ(in_unit)lib_input%psis(ipsi),psi(ipsi),
     $                lib_input%fs(ipsi),lib_input%ps(ipsi),
     $                lib_input%qs(ipsi)
         READ(in_unit)lib_input%rg(ipsi,:)
         READ(in_unit)lib_input%zg(ipsi,:)
         READ(in_unit)jaceq(ipsi,:)
         READ(in_unit)g12l(ipsi,:)
      ENDDO
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     copy variables.
c-----------------------------------------------------------------------
      ro=rmag
      zo=zmag
      psio=abs(tpsio)

c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(norpsi,psi,t,mu0p,q)
      DEALLOCATE(req,zeq,jaceq,g12l)


c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_input_chease

      SUBROUTINE lib_interface_set_equil_chease()
      INTEGER :: isurf,itheta
c-----------------------------------------------------------------------
c     copy variables.
c-----------------------------------------------------------------------
c      ro=SUM(lib_input%rg(1,1))/(lib_input%ntheta+1)
c      zo=SUM(lib_input%zg(1,1))/(lib_input%ntheta+1)
c      ro=lib_input%rg(1,1)
c      zo=lib_input%zg(1,1)
c      psio=abs(lib_input%psis(lib_input%nsurf))
c      lib_input%psis=abs(lib_input%psis/psio)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      CALL spline_alloc(sq_in,lib_input%nsurf-1,4)
      sq_in%xs=lib_input%psis
      sq_in%fs(:,1)=lib_input%fs
      sq_in%fs(:,2)=lib_input%ps
      sq_in%fs(:,3)=lib_input%qs
c      CALL spline_fit(sq_in,"extrap")
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      CALL bicube_alloc(rz_in,lib_input%nsurf-1,lib_input%ntheta-1,2)
      DO isurf=1,lib_input%nsurf
         DO itheta=1,lib_input%ntheta
            rz_in%fs(isurf-1,itheta-1,1)=lib_input%rg(isurf,itheta)
            rz_in%fs(isurf-1,itheta-1,2)=lib_input%zg(isurf,itheta)
         ENDDO
      ENDDO
20    FORMAT(481e14.6)
      CALL ascii_open(in_unit,"rg.dat","REPLACE")
      DO isurf=1,lib_input%nsurf
         WRITE(in_unit,20) lib_input%rg(isurf,:)
      ENDDO
      CALL ascii_close(in_unit)
      CALL ascii_open(in_unit,"zg.dat","REPLACE")
      DO isurf=1,lib_input%nsurf
         WRITE(in_unit,20) lib_input%zg(isurf,:)
      ENDDO
      CALL ascii_close(in_unit)

c      CALL bicube_fit(rz_in,"extrap","periodic")
c      CALL inverse_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lib_interface_set_equil_chease


      END MODULE lib_interface_mod
