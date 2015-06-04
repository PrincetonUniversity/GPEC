c-----------------------------------------------------------------------
c     NON AMBIPOLAR TRANSPORT and TOROIDAL TORQUE
c     global variable initialization
c-----------------------------------------------------------------------
      MODULE pentglobal_mod
      USE spline_mod
      USE bicube_mod
      USE cspline_mod
      USE ipglobal_mod
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE

!     This overrides io_mod values from equil/dcon/ipec
c      INTEGER, PARAMETER ::
c     $     in_unit=1,
c     $     out_unit=2,
c     $     bin_unit=3,
c     $     bin_2d_unit=4,
c     $     bounce_unit=5,
c     $     pitch_unit =6,
c     $     energy_unit=7,

c     !ALREADY IN equil's local_mod used by spline modules     
c      INTEGER, PARAMETER ::
c     $     r4=SELECTED_REAL_KIND(6,37), ! already in spline_mod->local_mod
c     $     r8=SELECTED_REAL_KIND(13,307)
c      REAL(r8), PARAMETER :: 
c     $     mp=1.672614e-27,
c     $     me=9.1091e-31,
c     $     e=1.6021917e-19,
c     $     ev=e
c      REAL(r8), PARAMETER :: 
c     $     pi=3.1415926535897932385_r8,
c     $     twopi=2*pi,pisq=pi*pi,
c     $     mu0=4e-7_r8*pi,
c     $     rtod=180/pi,dtor=pi/180,
c     $     alog10=2.302585093_r8
c      COMPLEX(r8), PARAMETER :: ifac=(0,1)    
  

      LOGICAL :: kinetic_flag,bounce_flag,energy_flag,pitch_flag,
     $     passing_flag,kolim_flag,neorot_flag,divxprp_flag,
     $     offset_flag,xlsode_flag,lmdalsode_flag,mdc2_flag,
     $     treg_flag

      INTEGER :: imass,icharge,zimp,zmass,nl

      REAL(r8) :: lsode_atol,lsode_rtol

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: bpar,xprp

      INTEGER :: nx, nt, nlmda
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xnorm,tnorm,lnorm

      CHARACTER(2) :: sl
      CHARACTER(32)  :: collision 
!      CHARACTER(128) :: outfile,infile,kfile
!      CHARACTER(512) :: outtext,outlbls

      TYPE(spline_type)  :: kin,ellipk,ellipe
      TYPE(bicube_type)  :: fmnl
      TYPE(cspline_type) :: ntv,ntvp

      CONTAINS

c-----------------------------------------------------------------------
c     function 1. specialf
c     special function from Park PRL 2009.
c-----------------------------------------------------------------------
      FUNCTION specialf(k,mnql)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(r8) :: specialf
      REAL(r8), INTENT(IN) :: k,mnql
      INTEGER :: ntheta,i
      REAL(r8) :: thetab
      TYPE(spline_type) :: fspl         
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      thetab = 2*ASIN(k)
      ntheta = 100+20*INT(ABS(mnql*2*thetab/twopi))
      CALL spline_alloc(fspl,ntheta,1)
      fspl%xs(:) = (/(i,i=1,ntheta+1)/)/(ntheta+2.0)
      fspl%xs(:) = COS(pi*fspl%xs(:)+pi)*thetab
      fspl%fs(:,1) = COS(mnql*fspl%xs)/
     $     SQRT(k*k-SIN(fspl%xs/2)*SIN(fspl%xs/2))
      CALL spline_fit(fspl,'extrap')
      CALL spline_int(fspl)
      specialf = fspl%fsi(fspl%mx,1)
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION specialf

c-----------------------------------------------------------------------
c     subroutine 1. pentg_fmnl.
c     precalculated matrix of special function in (kappa,m-nq-l) space.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentg_fmnl
      
      INTEGER :: k,t
      INTEGER, PARAMETER :: nk = 500
      INTEGER, PARAMETER :: nt = 500
      REAL(r8), DIMENSION(0:nk) :: kappa
      REAL(r8), DIMENSION(-nt:nt) :: mnql
      REAL(r8), DIMENSION(0:nk,-nt:nt) :: fkmnql

c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      WRITE(*,*) "Calculating F^1/2_mnl: kappa [0,1], m-nq-l [-100,100]"
      DO t = -nt,nt
         mnql(t)  = t/(nt/100.0)
      ENDDO
      DO k = 0,nk
         kappa(k) = (k+1.0)/(nk+2.0)
         DO t = -nt,nt
            fkmnql(k,t) = specialf(kappa(k),mnql(t))
         ENDDO
      ENDDO

      CALL bin_open(bin_2d_unit,"fkmnql.bin","UNKNOWN","REWIND","none")
      WRITE(bin_2d_unit)nk,2*nt
      WRITE(bin_2d_unit)kappa(:)
      WRITE(bin_2d_unit)mnql(:)
      WRITE(bin_2d_unit)fkmnql(:,:)
      CALL bin_close(bin_2d_unit)

      CALL ascii_open(out_unit,"fkmnql.out","UNKNOWN")
      WRITE(out_unit,*)"LARGE ASPECT RATIO F^1/2_mnl: "
      WRITE(out_unit,*)
      WRITE(out_unit,*)"nk = ",nk,"nt = ",2*nt
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,3(1x,a16))') "kappa","m-nq-l","f_mnl"
      DO k = 0,nk
         DO t = -nt,nt
            WRITE(out_unit,'(1x,3(1x,es16.8E3))') kappa(k),mnql(t),
     $           fkmnql(k,t)
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
     
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentg_fmnl

c-----------------------------------------------------------------------
c     subroutine 2. pentg_dealloc.
c     deallocate storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentg_dealloc
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(ntv%xs)) CALL cspline_dealloc(ntv)
      IF(ASSOCIATED(ntvp%xs)) CALL cspline_dealloc(ntvp)

      IF(ASSOCIATED(kin%xs)) CALL spline_dealloc(kin)
      IF(ASSOCIATED(ellipk%xs)) CALL spline_dealloc(ellipk)
      IF(ASSOCIATED(ellipe%xs)) CALL spline_dealloc(ellipe)
      IF(ASSOCIATED(fmnl%xs)) CALL bicube_dealloc(fmnl)

      IF(ALLOCATED(bpar)) DEALLOCATE(bpar)
      IF(ALLOCATED(xprp)) DEALLOCATE(xprp)
            
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentg_dealloc

      END MODULE pentglobal_mod
