c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     global definition
c-----------------------------------------------------------------------
      MODULE gpglobal_mod
      USE bicube_mod
      USE fspline_mod
      USE omp_lib
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      LOGICAL :: power_flag,fft_flag,edge_flag,pbrzphi_flag,
     $     eqbrzphi_flag,brzphi_flag,xrzphi_flag,divzero_flag,
     $     vbrzphi_flag,vvbrzphi_flag,div_flag,surface_flag,
     $     data_flag,harmonic_flag,mode_flag,resp_flag,
     $     bin_flag,bin_2d_flag,fixed_boundary_flag,reg_flag,
     $     fun_flag,flux_flag,vsbrzphi_flag,displacement_flag,
     $     chebyshev_flag,coil_flag,eigm_flag,bwp_pest_flag,verbose,
     $     debug_flag,timeit,kin_flag,con_flag,resp_induct_flag,
     $     netcdf_flag,ascii_flag,wv_farwall_flag,mutual_test_flag

      INTEGER :: mr,mz,mpsi,mstep,mpert,mband,mtheta,mthvac,mthsurf,
     $     mfix,mhigh,mlow,msing,nfm2,nths,nths2,lmpert,lmlow,lmhigh,
     $     power_b,power_r,power_bp,jsurf_in,jsurf_out,mlim_out,
     $     power_bin,power_rin,power_bpin,power_rcin,tmag_in,
     $     power_bout,power_rout,power_bpout,power_rcout,tmag_out,
     $     nn,info,resp_index,rstep,resp,psixy,nmin,nmax,mmin,mmax,
     $     nche,nchr,nchz,rsing,rnqty,rnx,max_linesout,
     $     nr,nz,malias,ol,ou,osing

      REAL(r8) :: ro,zo,psio,chi1,mthsurf0,psilow,psilim,qlim,
     $     qmin,qmax,seconds,rfac,eta,singfac_min,rmin,rmax,zlim,
     $     jac,jac1,q,q1,p,p1,bpfac,btfac,bfac,fac,reg_spot,
     $     amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,betap1,betap2,betap3,betat,betan,bt0,
     $     q0,qa,crnt,q95,shotnum,shottime,opsi1,opsi2

      CHARACTER(2) :: sn,ss
      CHARACTER(10) :: date,time,zone
      CHARACTER(16) :: jac_type,jac_in,jac_out,data_type
      CHARACTER(128) :: ieqfile,idconfile,rdconfile,dcon_dir

      INTEGER, PARAMETER :: hmnum=128
      REAL(r8), PARAMETER :: gauss=0.0001
      COMPLEX(r8), PARAMETER :: ione=1

      REAL(r8), DIMENSION(-hmnum:hmnum) :: sinmn
      REAL(r8), DIMENSION(-hmnum:hmnum) :: cosmn

      INTEGER, DIMENSION(8) :: values

      LOGICAL, DIMENSION(:), ALLOCATABLE :: sing_flag
      INTEGER, DIMENSION(:), ALLOCATABLE :: fixstep,mfac,lmfac
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: permeabindex,gdl

      REAL(r8), DIMENSION(:), ALLOCATABLE :: psifac,rhofac,qfac,singfac,
     $     r,z,theta,ee,surfee,surfei,rpsifac,
     $     surf_indev,vsurf_indev,fsurf_indev,surf_indinvev
     $     ,eft,efp
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xzpts,
     $     chperr,chpsqr,grri,grre,griw,grrw,gdr,gdz,gdpsi,gdthe,gdphi
      REAL(r8), DIMENSION(3,3) :: w,v


      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: et,ep,
     $     xsp_mn,xsp1_mn,xss_mn,xms_mn,bwp1_mn,xmp1_mn,
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,xmt_mn,bmt_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,xmz_mn,bmz_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn,
     $     chi_mn,che_mn,kax_mn,sbno_mn,sbno_fun,edge_mn,edge_fun,
     $     chy_mn,chx_mn,chw_mn,kaw_mn,mutual_indev,mutual_indinvev
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: wt,wt0,chp_mn,kap_mn,
     $     permeabev,chimats,chemats,flxmats,kaxmats,singbno_mn,
     $     plas_indev,plas_indinvev,reluctev,indrelev,permeabsv,
     $     surf_indmats,surf_indevmats,vsurf_indmats,fsurf_indmats,
     $     surf_indinvmats,surf_indinvevmats,surfet,surfep,
     $     amat,bmat,cmat,fmats,gmats,kmats,coil_indmat,coilmn,
     $     wft,wtraw,chymats,chxmats,chwmats,kawmats,
     $     mutual_indmats,mutual_indevmats,
     $     mutual_indinvmats,mutual_indinvevmats
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: chpmats,kapmats,
     $     plas_indmats,plas_indinvmats,permeabmats,diff_indmats,
     $     plas_indevmats,plas_indinvevmats,indrelmats,indrelevmats,
     $     reluctmats,reluctevmats,permeabevmats,permeabsvmats,
     $     permeabinvmats

      TYPE(spline_type) :: sq
      TYPE(bicube_type) :: psi_in,eqfun,rzphi
      TYPE(cspline_type) :: u1,u2,u3,u4,u5
      TYPE(cspline_type) :: smats,tmats,xmats,ymats,zmats
      TYPE(fspline_type) :: metric

      TYPE :: resist_type
      REAL(r8) :: e,f,h,m,g,k,eta,rho,taua,taur,di,dr,sfac,deltac
      END TYPE resist_type

      TYPE :: solution_type
      INTEGER :: msol
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: u
      END TYPE solution_type

      TYPE :: galsol_type
      LOGICAL :: gal_flag
      INTEGER :: mpert, tot_grids,mlow,mhigh,mtot
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: u
      END TYPE galsol_type

      TYPE :: fixfac_type
      INTEGER :: msol
      INTEGER, DIMENSION(:), ALLOCATABLE :: index
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: fixfac,transform,gauss
      END TYPE fixfac_type

      TYPE :: sing_type
      INTEGER :: msol_l,msol_r,jfix,jpert
      REAL(r8) :: psifac,q,q1
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: ca_l,ca_r
      TYPE(resist_type) :: restype
      END TYPE sing_type

      INTEGER :: gal_pert,gal_unit=105
      LOGICAL :: gal_check
      TYPE(galsol_type) :: galsol
      TYPE(solution_type), DIMENSION(:), ALLOCATABLE :: soltype,rsoltype
      TYPE(fixfac_type), DIMENSION(:), ALLOCATABLE :: fixtype
      TYPE(sing_type), DIMENSION(:), ALLOCATABLE :: singtype

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gpec_dealloc.
c     deallocate storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_dealloc

      INTEGER :: istep,ifix
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(sing_flag,fixstep,mfac,psifac,rhofac,qfac,singfac,
     $     et,ep,ee,surfet,surfep,surfee,surfei,lmfac,r,z,theta,
     $     surf_indev,plas_indev,permeabev,edge_mn,edge_fun,chperr)
      DEALLOCATE(wt,chimats,chemats,chpmats,kapmats,kaxmats,flxmats,
     $     surf_indmats,surf_indevmats,plas_indmats,permeabmats,
     $     plas_indevmats,permeabevmats,chxmats,chwmats,
     $     mutual_indmats,mutual_indevmats,mutual_indinvmats,
     $     mutual_indinvevmats)

      CALL spline_dealloc(sq)
      CALL bicube_dealloc(eqfun)
      CALL bicube_dealloc(rzphi)
      CALL cspline_dealloc(u1)
      CALL cspline_dealloc(u2)
      CALL fspline_dealloc(metric)

      IF(ALLOCATED(soltype))THEN
         DO istep=0,mstep
            DEALLOCATE(soltype(istep)%u)
         ENDDO
         DEALLOCATE(soltype)
      ENDIF
      IF(ALLOCATED(fixtype))THEN
         DO ifix=0,mfix
            IF(ALLOCATED(fixtype(ifix)%fixfac))
     $         DEALLOCATE(fixtype(ifix)%fixfac)
            IF(ALLOCATED(fixtype(ifix)%index))
     $         DEALLOCATE(fixtype(ifix)%index)
            IF(ALLOCATED(fixtype(ifix)%transform))
     $         DEALLOCATE(fixtype(ifix)%transform)
            IF(ALLOCATED(fixtype(ifix)%gauss))
     $         DEALLOCATE(fixtype(ifix)%gauss)
         ENDDO
         DEALLOCATE(fixtype)
      ENDIF

      IF (psixy == 1) CALL bicube_dealloc(psi_in)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpec_dealloc
c-----------------------------------------------------------------------
c     subprogram 2. gpec_stop.
c     program termination.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_stop(message)

      CHARACTER(*), INTENT(IN) :: message
      INTEGER :: hrs,mins,secs,cpuhrs,cpumins,cpusecs
      REAL(4) :: cpu_time

      cpu_time = 0
      CALL gpec_timer(-2,cpu_time)
      cpusecs = INT(cpu_time)
      cpuhrs = cpusecs/(60*60)
      cpumins = (cpusecs-cpuhrs*60*60)/60
      cpusecs = cpusecs-cpuhrs*60*60-cpumins*60
      
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
     $     -seconds
      secs = int(seconds)
      hrs = secs/(60*60)
      mins = (secs-hrs*60*60)/60
      secs = secs-hrs*60*60-mins*60
      IF (verbose) THEN
         ! write total time since start

#ifdef _OPENMP
         IF (OMP_GET_MAX_THREADS() > 1) THEN
            IF(cpuhrs>0)THEN
               PRINT *,"Total CPU time for GPEC = ",cpuhrs," hours, ",
     $           cpumins," minutes, ",cpusecs," seconds"
            ELSEIF(cpumins>0)then
               PRINT *,"Total CPU time for GPEC = ",
     $           cpumins," minutes, ",cpusecs," seconds"
            ELSE
               PRINT *,"Total CPU time for GPEC = ",cpusecs," seconds"
            ENDIF
         ENDIF
#endif
         
         IF(hrs>0)THEN
             PRINT *,"Total wall time for GPEC = ",hrs," hours, ",
     $           mins," minutes, ",secs," seconds"
         ELSEIF(mins>0)then
             PRINT *,"Total wall time for GPEC = ",
     $           mins," minutes, ",secs," seconds"
         ELSE
             PRINT *,"Total wall time for GPEC = ",secs," seconds"
         ENDIF
         WRITE(*,'(1x,2a)')'GPEC STOP => ',TRIM(message)
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE gpec_stop
c-----------------------------------------------------------------------
c     subprogram 3. gpec_timer.
!-----------------------------------------------------------------------
c     *DESCRIPTION:
c        Handles machine-dependent timing statistics.
c
c     *ARGUMENTS:
c        mode : integer, in
c            mode = 0 starts timer
c            mode = 1 writes total time from last start
c            mode = 2 writes split & total time from start
c            mode = -2 resets the split without writing
c
c     *OPTIONAL ARGUMENTS:
c        opunit : integer, in
c            Output written to this unit (default to terminal)
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gpec_timer(mode,timeout,opunit)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mode
      REAL(4), INTENT(OUT), OPTIONAL :: timeout
      INTEGER, INTENT(IN), OPTIONAL :: opunit

      REAL(4) :: time
      REAL(4), SAVE :: start,split
      INTEGER :: hrs,mins,secs,msec

      ! get time
      CALL CPU_TIME(time)

      IF(mode==0)THEN ! start timer
         start=time
         split=0
      ELSEIF(mode==-2)THEN ! Reset split
         split=time-start
         IF(PRESENT(timeout))THEN
            timeout=split
         ENDIF
      ELSEIF(mode>0)THEN ! outputs
         IF(mode==2)THEN ! write split (time since last call)
            split=time-start-split
            secs = INT(split)
            hrs = secs/(60*60)
            mins = (secs-hrs*60*60)/60
            secs = secs-hrs*60*60-mins*60
            msec = int((split-secs)*1e3)
            IF(PRESENT(opunit))THEN
               write(opunit,"(1x,a,1p,es10.3,a)")"  split cpu time = ",
     $            split," seconds"
            ELSE
               IF(hrs>0)THEN
                  PRINT *,"  split cpu time = ",hrs," hours, ",
     $               mins," minutes, ",secs," seconds"
               ELSEIF(mins>0)THEN
                  PRINT *,"  split cpu time = ",
     $               mins," minutes, ",secs," seconds"
               ELSEIF(secs>0)THEN
                  PRINT *,"  split cpu time = ",secs," seconds"
               ELSE
                  PRINT *,"  split cpu time = ",msec," miliseconds"
               ENDIF
            ENDIF
         ENDIF
         ! write total time since start
         split=time-start
         secs = INT(split)
         hrs = secs/(60*60)
         mins = (secs-hrs*60*60)/60
         secs = secs-hrs*60*60-mins*60
         msec = int((split-real(secs,4))*1e3)
         IF(PRESENT(opunit))THEN
            write(opunit,"(1x,a,1p,es10.3,a)")"  total cpu time = ",
     $         split," seconds"
         ELSE
            IF(hrs>0)THEN
               PRINT *,"  total cpu time = ",hrs," hours, ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(mins>0)THEN
               PRINT *,"  total cpu time = ",
     $            mins," minutes, ",secs," seconds"
            ELSEIF(secs>0)THEN
               PRINT *,"  total cpu time = ",secs," seconds"
            ELSE
               PRINT *,"  total cpu time = ",msec," miliseconds"
            ENDIF
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpec_timer

      END MODULE gpglobal_mod
