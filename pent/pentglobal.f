c-----------------------------------------------------------------------
c     NON AMBIPOLAR TRANSPORT and TOROIDAL TORQUE
c     global variable initialization
c-----------------------------------------------------------------------
      MODULE pentglobal_mod
      USE spline_mod
      USE bicube_mod
      USE cspline_mod
      USE local_mod
      USE ipglobal_mod
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      LOGICAL :: kinetic_flag,bounce_flag,energy_flag,pitch_flag,
     $     passing_flag,kolim_flag,neorot_flag,divxprp_flag,
     $     offset_flag,xlsode_flag

      INTEGER :: imass,icharge,zimp,zmass

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: bpar,xprp
      REAL(r8), PARAMETER :: emass = 9.109e-31, echarge = -1.602e-19,
     $     pmass = 1.673e-27

c$$$      INTEGER, PARAMETER :: nx = 200, nt = 101, nlmda = 101
c$$$      REAL(r8), DIMENSION(0:nx) :: xnorm =(/(i*1.0/nx,i=0,nx)/)
c$$$      REAL(r8), DIMENSION(0:nt) :: tnorm =(/(i*1.0/nt,i=0,nt)/)
c$$$      REAL(r8), DIMENSION(0:nlmda) :: lnorm =(/(i*1.0/nlmda,i=0,nlmda)/)
      INTEGER :: nx, nt, nlmda
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xnorm,tnorm,lnorm

      CHARACTER(128) :: outfile,infile,kfile
      CHARACTER(512) :: outtext,outlbls

      TYPE(spline_type)  :: kin,ellipk
      TYPE(cspline_type) :: ntv,ntvp,lar,larp
      TYPE(bicube_type)  :: fmnl

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
      SUBROUTINE pentg_dealloc(gen_flag,lar_flag)
      LOGICAL, INTENT(IN) :: gen_flag,lar_flag
      INTEGER :: istep,ifix,ising
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      IF(gen_flag)THEN
         CALL cspline_dealloc(ntv)
         CALL cspline_dealloc(ntvp)
      ENDIF
      IF(lar_flag)THEN
         CALL cspline_dealloc(lar)
         CALL cspline_dealloc(larp)
         CALL spline_dealloc(ellipk)
         CALL bicube_dealloc(fmnl)
      ENDIF

      IF(ALLOCATED(bpar)) DEALLOCATE(bpar)
      IF(ALLOCATED(xprp)) DEALLOCATE(xprp)
      CALL spline_dealloc(kin)
            
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentg_dealloc

      END MODULE pentglobal_mod
