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
      LOGICAL :: power_flag,fft_flag,edge_flag,pbrzphi_flag,
     $     eqbrzphi_flag,brzphi_flag,xrzphi_flag,divzero_flag,
     $     vbrzphi_flag,vvbrzphi_flag,pbrzphi_flag,div_flag,
     $     data_flag,harmonic_flag,mode_flag,svd_flag,resp_flag,
     $     bin_flag,bin_2d_flag,fixed_boundary_flag,reg_flag,
     $     fun_flag,flux_flag,vsbrzphi_flag,displacement_flag,
     $     chebyshev_flag,coil_flag,eigm_flag,bwp_pest_flag,verbose
      INTEGER :: mr,mz,mpsi,mstep,mpert,mband,mtheta,mthvac,mthsurf,
     $     mfix,mhigh,mlow,msing,nfm2,nths2,lmpert,lmlow,lmhigh,
     $     power_b,power_r,power_bp,jsurf_in,jsurf_out,
     $     power_bin,power_rin,power_bpin,power_rcin,tmag_in,
     $     power_bout,power_rout,power_bpout,power_rcout,tmag_out,
     $     nn,info,resp_index,rstep,resp,psixy,nmin,nmax,mmin,mmax,
     $     nche,nchr,nchz,rsing,rnqty,rnx,
     $     pmode,p1mode,rmode,dmode,d1mode,fmode,smode ! LOGAN

      REAL(r8) :: ro,zo,psio,chi1,mthsurf0,psilow,psilim,qlim,
     $     qmin,qmax,seconds,rfac,eta,singfac_min,rmin,rmax,zlim,
     $     jac,jac1,q,q1,p,p1,bpfac,btfac,bfac,fac,sing_spot,reg_spot

      CHARACTER(2) :: sn,ss
      CHARACTER(10) :: date,time,zone
      CHARACTER(16) :: jac_type,jac_in,jac_out,data_type
      CHARACTER(128) :: ieqfile,idconfile,ivacuumfile,rdconfile

      INTEGER, PARAMETER :: hmnum=128
      REAL(r8), PARAMETER :: gauss=0.0001
      COMPLEX(r8), PARAMETER :: ione=1

      REAL(r8), DIMENSION(-hmnum:hmnum) :: sinmn
      REAL(r8), DIMENSION(-hmnum:hmnum) :: cosmn
      REAL(r8), DIMENSION(0:2*hmnum) :: svdfac = 0
      

      INTEGER, DIMENSION(8) :: values

      LOGICAL, DIMENSION(:), POINTER :: sing_flag
      INTEGER, DIMENSION(:), POINTER :: fixstep,mfac,lmfac
      INTEGER, DIMENSION(:,:), POINTER :: permeabindex,gdl

      REAL(r8), DIMENSION(:), POINTER :: psifac,rhofac,qfac,singfac,
     $     r,z,theta,et,ep,ee,surfee,surfei,rpsifac,
     $     surf_indev,vsurf_indev,fsurf_indev
     $     ,eft,efp,perms ! LOGAN
      REAL(r8), DIMENSION(:,:), POINTER :: surfet,surfep,
     $     chperr,chpsqr,plas_indev,reluctev,indrelev,grri,grre,
     $     gdr,gdz,gdpsi,gdthe,gdphi
      REAL(r8), DIMENSION(3,3) :: w,v


      COMPLEX(r8), DIMENSION(:), POINTER ::
     $     xsp_mn,xsp1_mn,xss_mn,xms_mn,bwp1_mn,xmp1_mn,
     $     xwp_mn,xwt_mn,xwz_mn,bwp_mn,bwt_mn,bwz_mn,xmt_mn,bmt_mn,
     $     xvp_mn,xvt_mn,xvz_mn,bvp_mn,bvt_mn,bvz_mn,xmz_mn,bmz_mn,
     $     xno_mn,xta_mn,xpa_mn,bno_mn,bta_mn,bpa_mn,
     $     xrr_mn,xrz_mn,xrp_mn,brr_mn,brz_mn,brp_mn,
     $     chi_mn,che_mn,kax_mn,sbno_mn,sbno_fun,
     $     edge_mn,edge_fun
      COMPLEX(r8), DIMENSION(:,:), POINTER :: wt,chp_mn,kap_mn,
     $     permeabev,chimats,chemats,flxmats,kaxmats,singbno_mn,
     $     surf_indmats,surf_indevmats,vsurf_indmats,fsurf_indmats,
     $     amat,bmat,cmat,fmats,gmats,kmats,t1v,t2v,t3v,w1v,fldflxmat
     $     ,wft,permv ! LOGAN
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: chpmats,kapmats,
     $     plas_indmats,permeabmats,diff_indmats,reluctmats,
     $     plas_indevmats,permeabevmats,reluctevmats,
     $     indrelmats,indrelevmats

      TYPE(spline_type) :: sq
      TYPE(bicube_type) :: psi_in,eqfun,rzphi
      TYPE(cspline_type) :: u1,u2
      TYPE(cspline_type) :: smats,tmats,xmats,ymats,zmats
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
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: ca_l,ca_r
      TYPE(resist_type) :: restype
      END TYPE sing_type

      TYPE(solution_type), DIMENSION(:), POINTER :: soltype,rsoltype
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
     $     et,ep,ee,surfet,surfep,surfee,surfei,lmfac,r,z,theta,
     $     surf_indev,plas_indev,permeabev,edge_mn,edge_fun,chperr)
      DEALLOCATE(wt,chimats,chemats,chpmats,kapmats,kaxmats,flxmats,
     $     surf_indmats,surf_indevmats,plas_indmats,permeabmats,
     $     plas_indevmats,permeabevmats)

      CALL spline_dealloc(sq)
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

      IF (psixy == 1) CALL bicube_dealloc(psi_in)
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
      IF(verbose) WRITE(*,*)"Total cpu time for ipec=",seconds,"seconds"
      IF(verbose) WRITE(*,'(1x,2a)')'IPEC_STOP=>',TRIM(message)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE ipec_stop

      END MODULE ipglobal_mod
