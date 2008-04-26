c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     IPGLOBAL: global definition
c-----------------------------------------------------------------------
      MODULE ipglobal_mod
      USE bicube_mod
      USE fspline_mod
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      LOGICAL :: power_flag,fft_flag,edge_flag
      INTEGER :: mr,mz,mpsi,mstep,mpert,mband,mtheta,mthvac,mthsurf,
     $     mfix,mhigh,mlow,msing,nfm2,nths2,lmpert,lmlow,lmhigh,
     $     power_b,power_r,power_bp,nn,info,modelnum,
     $     rstep,resp

      REAL(r8) :: ro,zo,psio,chi1,mthsurf0,psilow,psilim,qlim,
     $     qmin,qmax,seconds,rfac,eta,singfac_min,maxdbratio,
     $     jac,jac1,q,q1,p,p1,bpfac,btfac,bfac,fac,dist,bdist

      CHARACTER(1) :: sn
      CHARACTER(10) :: date,time,zone
      CHARACTER(128) :: ieqfile,idconfile,ivacuumfile,errtype
      CHARACTER(128), DIMENSION(10) :: infiles

      REAL(r8), PARAMETER :: mili=0.001,gauss=0.0001
      COMPLEX(r8), PARAMETER :: ione=1

      INTEGER, DIMENSION(8) :: values

      LOGICAL, DIMENSION(:), POINTER :: sing_flag
      INTEGER, DIMENSION(:), POINTER :: fixstep,mfac,lmfac
      INTEGER, DIMENSION(:,:), POINTER :: permeabindex,gdl

      REAL(r8), DIMENSION(:), POINTER :: psifac,rhofac,qfac,singfac,
     $     r,z,theta,et,ep,ee,surfee,surfei,
     $     surf_indev,vsurf_indev,fsurf_indev
      REAL(r8), DIMENSION(:,:), POINTER :: surfet,surfep,surfes,
     $     chperr,plas_indev,reluctev,grri,grre,
     $     gdr,gdz,gdpsi,gdthe,gdphi
      REAL(r8), DIMENSION(3,3) :: w,v


      COMPLEX(r8), DIMENSION(:), POINTER ::
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,
     $     xvs_mn,xwp1_mn,bwp1_mn,nbwp1_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn,
     $     chi_mn,che_mn,flx_mn,kax_mn,
     $     edge_mn,edge_fun
      COMPLEX(r8), DIMENSION(:,:), POINTER :: wt,rt,
     $     chp_mn,kap_mn,pjp_mn,permeabev,
     $     chimats,chemats,flxmats,kaxmats,
     $     surf_indmats,surf_indevmats,vsurf_indmats,fsurf_indmats,
     $     amat,bmat,cmat,fmats,gmats,kmats
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: chpmats,kapmats,
     $     plas_indmats,permeabmats,diff_indmats,reluctmats,
     $     plas_indevmats,permeabevmats,reluctevmats,pjpmats

      TYPE(spline_type) :: sq,ffun
      TYPE(bicube_type) :: psi_in,eqfun,rzphi
      TYPE(cspline_type) :: u1,u2
      TYPE(fspline_type) :: metric

      TYPE :: resist_type
      REAL(r8) :: e,f,h,m,g,k,eta,rho,taua,taur,di,dr,sfac,deltac
      END TYPE resist_type

      TYPE :: solution_type
      INTEGER :: msol
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u
      END TYPE solution_type

      TYPE :: fixfac_type
      INTEGER :: msol
      INTEGER, DIMENSION(:), POINTER :: index
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fixfac,transform,gauss
      END TYPE fixfac_type

      TYPE :: sing_type
      INTEGER :: msol_l,msol_r,jfix,jpert
      REAL(r8) :: psifac,q,q1

      TYPE(resist_type) :: restype
      END TYPE sing_type

      TYPE(solution_type), DIMENSION(:), POINTER :: soltype
      TYPE(fixfac_type), DIMENSION(:), POINTER :: fixtype
      TYPE(sing_type), DIMENSION(:), POINTER :: singtype

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipec_dealloc.
c     deallocate storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ipec_dealloc

      INTEGER :: istep,ifix,ising
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(sing_flag,fixstep,mfac,psifac,rhofac,qfac,singfac,
     $     et,ep,ee,surfet,surfep,surfee,surfei,surfes,lmfac,r,z,theta,
     $     surf_indev,plas_indev,permeabev,reluctev,edge_mn,edge_fun,
     $     chperr)
      DEALLOCATE(wt,rt,
     $     chimats,chemats,chpmats,kapmats,kaxmats,flxmats,pjpmats,
     $     surf_indmats,surf_indevmats,
     $     plas_indmats,diff_indmats,permeabmats,reluctmats,
     $     plas_indevmats,permeabevmats,reluctevmats)

      CALL spline_dealloc(sq)
      CALL spline_dealloc(ffun)
      CALL bicube_dealloc(psi_in)
      CALL bicube_dealloc(eqfun)
      CALL bicube_dealloc(rzphi)
      CALL cspline_dealloc(u1)
      CALL cspline_dealloc(u2)
      CALL fspline_dealloc(metric)

      DO istep=0,mstep
         DEALLOCATE(soltype(istep)%u)
      ENDDO
      DO ifix=1,mfix
         DEALLOCATE(fixtype(ifix)%fixfac,fixtype(ifix)%index)
      ENDDO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipec_dealloc
c-----------------------------------------------------------------------
c     subprogram 2. ipec_stop.
c     program termination.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ipec_stop(message)

      CHARACTER(*), INTENT(IN) :: message
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
     $     -seconds
      WRITE(*,*)"Total cpu time for ipec=",seconds,"seconds"
      WRITE(*,'(1x,2a)')'IPEC_STOP=>',TRIM(message)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE ipec_stop

      END MODULE ipglobal_mod
