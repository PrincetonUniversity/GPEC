c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     write various output results of gpec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. gpout_mod
c      1. gpout_response
c      2. gpout_singcoup
c      3. gpout_control
c      4. gpout_singfld
c      5. gpout_vsingfld
c      6. gpout_dw
c      7. gpout_dw_matrix
c      8. gpout_pmodb
c      9. gpout_xbnormal
c     10. gpout_vbnormal
c     11. gpout_xtangent
c     12. gpout_xbrzphi
c     13. gpout_vsbrzphi
c     14. gpout_xbrzphifun
c     15. gpout_arzphifun
c     16. gpout_clebsch
c     17. gpout_control_filter
c     18. gpout_qrv
c     19. check
c     20. gpout_init_netcdf
c     21. gpout_close_netcdf
c-----------------------------------------------------------------------
c     subprogram 0. gpout_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gpout_mod
      USE gpresp_mod
      USE gpvacuum_mod
      USE gpdiag_mod
      USE field_mod
      USE netcdf

      IMPLICIT NONE

      ! harvest variables
      INCLUDE 'harvest_lib.inc77'
      INTEGER  :: ierr
      INTEGER  :: hint = 0
      REAL(r8) :: hdbl = 0
      CHARACTER(LEN=16)    :: hkey
      CHARACTER(LEN=50000) :: hnml
      CHARACTER(LEN=65507) :: hlog
      CHARACTER, PARAMETER :: nul = char(0)

      ! netcdf ids
      INTEGER :: mncid,fncid,cncid
      CHARACTER(64) :: mncfile,fncfile,cncfile

      ! module wide output variables
      INCLUDE "version.inc"
      INTEGER, PARAMETER :: nsingcoup=5             ! number of resonant coupling models
      LOGICAL :: singcoup_set = .FALSE.             ! whether singcoup subroutine has been run
      REAL(r8) :: jarea                             ! control surface area
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE ::
     $   singcoup,                                  ! Resonant coupling to external field (b_x) in working coordinates
     $   singcoup_out_vecs,                         ! Right singular vectors of power normalized resonant coupling matrices in output coordinates
     $   localcoup_out_vecs                         ! Right singular vectors of local power normalized resonant coupling matrices in output coordinates
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::
     $   fldflxmat                                  ! convert power normalized field to area normalized flux (i.e. increase by sqrt(A_m)/sqrt(A) weighting)

      CONTAINS
      
      !-----------------------------------------------------------------
      function str(k,fmt)
      !-----------------------------------------------------------------
      !*DESCRIPTION:
      !   Convert an integer to string.
      !*ARGUMENTS:
      !   k : integer. Maximum unit.
      !   fmt : str. Format specification.
      !*RETURNS:
      !     str. length 20.
      !-----------------------------------------------------------------
          character(len=20) str
          integer, intent(in) :: k
          character(*),intent(in) :: fmt
          
          write (str, trim(fmt)) k
          str = adjustl(str)
      end function str
      !-----------------------------------------------------------------
      function to_upper(strIn) result(strOut)
      !-----------------------------------------------------------------
      !*DESCRIPTION:
      !   Capitalize a string.
      !   Adapted from http://www.star.le.ac.uk/~cgp/fortran.html
      !   Original author: Clive Page
      !*ARGUMENTS:
      !   strIn : str. Original string.
      !*RETURNS:
      !   strOut : str. New capitalized string.
      !-----------------------------------------------------------------
         implicit none
         character(len=*), intent(in) :: strIn
         character(len=len(strIn)) :: strOut
         integer :: i,j

         do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j>= iachar("a") .and. j<=iachar("z") ) then
                   strOut(i:i) = achar(iachar(strIn(i:i))-32)
              else
                   strOut(i:i) = strIn(i:i)
              end if
         end do
      end function to_upper

c-----------------------------------------------------------------------
c     subprogram 1. gpout_response.
c     write basic information.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_response(rout,bpout,bout,rcout,tout,jout)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout,jout
      INTEGER :: lwork,i,j
      REAL(r8), DIMENSION(mpert) :: sts,s
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,b,u,vt

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(3*mpert) :: work

      INTEGER :: idid,mdid,edid,l_id,r_id,la_id,p_id,w_id,wi_id
      COMPLEX(r8), DIMENSION(lmpert) :: vL,vL1,vLi,vP,vP1,vR,vW,templ
      COMPLEX(r8), DIMENSION(mpert,mpert) :: matmm
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat
      COMPLEX(r8), DIMENSION(0:mthsurf,mpert) :: wtfun,ilfun,rfun,pfun

      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*) "Writing response matrices"
c-----------------------------------------------------------------------
c     svd analysis for permeability matrix.
c-----------------------------------------------------------------------
      a=permeabmats(resp_index,:,:)
      work=0
      rwork=0
      s=0
      u=0
      vt=0
      lwork=3*mpert
      CALL zgesvd('S','S',mpert,mpert,a,mpert,s,u,mpert,vt,mpert,
     $     work,lwork,rwork,info)
      DO i=1,mpert
         sts(i)=-REAL(et(i),r8)/(surfei(i)+surfee(i))
         IF (sts(i) > 0) s(i)=-s(i)
         s(i)=-1/s(i)
      ENDDO
c-----------------------------------------------------------------------
c     log eigenvalues with harvest
c-----------------------------------------------------------------------
      ierr=set_harvest_payload_dbl_array(hlog,"s_P"//nul,
     $     REAL(permeabev(resp_index,:)),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_L"//nul,
     $     surf_indev(:),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_Lambda"//nul,
     $     plas_indev(resp_index,:),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_rho"//nul,
     $     reluctev(resp_index,:),mpert)
c-----------------------------------------------------------------------
c     fundamental matrices in netcdf
c-----------------------------------------------------------------------
      IF(netcdf_flag)THEN
         CALL check( nf90_open(mncfile,nf90_write,mncid) )
         CALL check( nf90_inq_dimid(mncid,"i",idid) )
         CALL check( nf90_inq_dimid(mncid,"m",mdid) )
         CALL check( nf90_inq_dimid(mncid,"mode",edid) )
         CALL check( nf90_redef(mncid))
         CALL check( nf90_def_var(mncid,"L",nf90_double,
     $                  (/mdid,edid,idid/),l_id) )
         CALL check( nf90_put_att(mncid,l_id,"long_name",
     $       "Surface inductance") )
         CALL check( nf90_def_var(mncid,"Lambda",nf90_double,
     $                  (/mdid,edid,idid/),la_id) )
         CALL check( nf90_put_att(mncid,la_id,"long_name",
     $       "Plasma inductance") )
         CALL check( nf90_def_var(mncid,"P",nf90_double,
     $                  (/mdid,edid,idid/),p_id) )
         CALL check( nf90_put_att(mncid,p_id,"long_name",
     $       "Permeability") )
         CALL check( nf90_def_var(mncid,"rho",nf90_double,
     $                  (/mdid,edid,idid/),r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $       "Reluctance") )
         CALL check( nf90_def_var(mncid,"M_w",nf90_double,
     $                  (/mdid,edid,idid/),w_id) )
         CALL check( nf90_put_att(mncid,w_id,"long_name",
     $       "Mutual Inductance with Wall") )
         CALL check( nf90_def_var(mncid,"M_w_inv",nf90_double,
     $                  (/mdid,edid,idid/),wi_id) )
         CALL check( nf90_put_att(mncid,wi_id,"long_name",
     $       "Inverse of Mutual Inductance with Wall") )
         CALL check( nf90_enddef(mncid) )
         matmm = surf_indmats
         CALL check( nf90_put_var(mncid,l_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         matmm = plas_indmats(resp_index,:,:)
         CALL check( nf90_put_var(mncid,la_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         matmm = permeabmats(resp_index,:,:)
         CALL check( nf90_put_var(mncid,p_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         matmm = reluctmats(resp_index,:,:)
         CALL check( nf90_put_var(mncid,r_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         matmm = mutual_indmats
         CALL check( nf90_put_var(mncid,w_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         matmm = mutual_indinvmats
         CALL check( nf90_put_var(mncid,wi_id,RESHAPE((/REAL(matmm),
     $               AIMAG(matmm)/),(/mpert,mpert,2/))) )
         CALL check( nf90_close(mncid) )
      ENDIF
c-----------------------------------------------------------------------
c     All possible descriptions in large ascii files
c-----------------------------------------------------------------------
      IF(.NOT. ascii_flag)THEN
         IF(timeit) CALL gpec_timer(2)
         RETURN
      ENDIF

      CALL ascii_open(out_unit,"gpec_response_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_RESPONSE: Response parameters"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(3(1x,a12,I4))')
     $     "mpert =",mpert,"mlow =",mlow,"mhigh =",mhigh
      WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $     "psilim =",psilim,"qlim =",qlim

      WRITE(out_unit,*)
      WRITE(out_unit,*)"Energy for dcon eigenmodes"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,9(1x,a12))')"mode","ev0","ev1","iv1",
     $     "ep0","ep1","ep2","ep3","ep4","et0"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,9(1x,es12.3))')i,ee(i),surfee(i),
     $        surfei(i),REAL(ep(i)),REAL(surfep(1,i)),REAL(surfep(2,i)),
     $        REAL(surfep(3,i)),REAL(surfep(4,i)),REAL(et(i))
      ENDDO
      WRITE(out_unit,*)

      WRITE(out_unit,*)
      WRITE(out_unit,*)"Torque for dcon eigenmodes"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,5(1x,a12))')"mode",
     $     "ep0","ep1","ep2","ep3","ep4"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,5(1x,es12.3))')i,-2.0*nn*AIMAG(ep(i)),
     $        -2.0*nn*AIMAG(surfep(1,i)),-2.0*nn*AIMAG(surfep(2,i)),
     $        -2.0*nn*AIMAG(surfep(3,i)),-2.0*nn*AIMAG(surfep(4,i))
      ENDDO
      WRITE(out_unit,*)

      WRITE(out_unit,*)"Stability indices"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,2(1x,a12))')"mode","s","se"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,2(1x,es12.3))')i,s(i),sts(i)
      ENDDO
      WRITE(out_unit,*)

      WRITE(out_unit,*)"Eigenvalues (e) and Singular Values (s)"
      WRITE(out_unit,*)" jac_type = ",jac_type
      WRITE(out_unit,*)"  L = Vacuum Inductance"
      WRITE(out_unit,*)"  Lambda = Inductance"
      WRITE(out_unit,*)"  P = Permeability   *Complex (not Hermitian)"
      WRITE(out_unit,*)"  P1= Response (P-I) *Singular values"
      WRITE(out_unit,*)"  rho = Reluctance (power norm)"
      WRITE(out_unit,*)"  W = Displacement energy"
      WRITE(out_unit,*)"  F = Flux energy"
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a4,8(1x,a16))')"mode",
     $  "e_L","e_Lambda","real(e_P)","imag(e_P)","s_P","e_rho",
     $  "e_W","e_iLmda"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,8(1x,es17.8))') i,
     $     surf_indev(i),REAL(plas_indev(resp_index,i)),
     $     REAL(permeabev(resp_index,i)),AIMAG(permeabev(resp_index,i)),
     $     REAL(permeabsv(resp_index,i)),REAL(reluctev(resp_index,i)),
     $     REAL(et(i)),REAL(plas_indinvev(resp_index,i))
      ENDDO
      WRITE(out_unit,*)

      a = wt
      DO i=1,mpert
         a(:,i) = a(:,i)*(chi1*twopi*ifac*(mfac-nn*qlim))
      ENDDO
      WRITE(out_unit,*)"Eigenvectors"
      WRITE(out_unit,*)" jac_type = ",jac_type
      WRITE(out_unit,*)"  L = Vacuum Inductance"
      WRITE(out_unit,*)"  Lambda = Inductance (iLmda = Lambda^-1)"
      WRITE(out_unit,*)"  P = Permeability"
      WRITE(out_unit,*)"  rho = Reluctance (power norm)"
      WRITE(out_unit,*)"  W = Displacement energy"
      WRITE(out_unit,*)"  F = Flux energy"
      WRITE(out_unit,*)
      WRITE(out_unit,'(2(1x,a4),16(1x,a16))')"mode","m",
     $   "real(K_x^L)","imag(K_x^L)",
     $   "real(K_x^Lambda)","imag(K_x^Lambda)",
     $   "real(Phi_x^P)","imag(Phi_x^P)",
     $   "real(Phi_x^PS)","imag(Phi_x^PS)",
     $   "real(Phi_x^rho)","imag(Phi_x^rho)",
     $   "real(X^W)","imag(X^W)","real(Phi^W)","imag(Phi^W)",
     $   "real(Phi^iLmda)","imag(Phi^iLmda)"
      DO i=1,mpert
         DO j=1,mpert
            WRITE(out_unit,'(2(1x,I4),16(es17.8e3))')i,mfac(j),
     $           surf_indevmats(j,i),
     $           plas_indevmats(resp_index,j,i),
     $           permeabevmats(resp_index,j,i),
     $           permeabsvmats(resp_index,j,i),
     $           reluctevmats(resp_index,j,i),
     $           wt(j,i),a(j,i),
     $           plas_indinvevmats(resp_index,j,i)
         ENDDO
      ENDDO
      WRITE(out_unit,*)

      IF ((jac_out /= jac_type).OR.(tout==0)) THEN
         WRITE(out_unit,*)"Eigenvectors"
         WRITE(out_unit,*)" jac_out  = ",jac_out
         WRITE(out_unit,*)" tmag_out = ",tout
         WRITE(out_unit,*)" jsurf_out = ",jout
         WRITE(out_unit,*)"  L = Vacuum Inductance"
         WRITE(out_unit,*)"  Lambda = Inductance"
         WRITE(out_unit,*)"  P = Permeability"
         WRITE(out_unit,*)"  rho = Reluctance (power norm)"
         WRITE(out_unit,*)"  W = Displacement energy"
         WRITE(out_unit,*)"  F = Flux energy"
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a4),14(1x,a16))')"mode","m",
     $     "real(V_L)","imag(V_L)","real(V_Lambda)","imag(V_Lambda)",
     $     "real(V_P)","imag(V_P)","real(V_PS)","imag(V_PS)",
     $     "real(V_rho)","imag(V_rho)",
     $     "real(V_W)","imag(V_W)",
     $     "real(Phi^iLmda)","imag(Phi^iLmda)"
         DO i=1,mpert
            DO j=1,mpert
                IF ((mlow-lmlow+j>=1).AND.(mlow-lmlow+j<=lmpert)) THEN
                   vL(mlow-lmlow+j) =surf_indevmats(j,i)
                   vL1(mlow-lmlow+j)=plas_indevmats(resp_index,j,i)
                   vLi(mlow-lmlow+j)=plas_indinvevmats(resp_index,j,i)
                   vP(mlow-lmlow+j) =permeabevmats(resp_index,j,i)
                   vP1(mlow-lmlow+j)=permeabsvmats(resp_index,j,i)
                   vR(mlow-lmlow+j) =reluctevmats(resp_index,j,i)
                   vW(mlow-lmlow+j) =wt(j,i)
                ENDIF
             ENDDO
             CALL gpeq_bcoords(psilim,vL,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vL1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vLi,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vP,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vP1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vR,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL gpeq_bcoords(psilim,vW,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
            DO j=1,lmpert
               WRITE(out_unit,'(2(1x,I4),14(es17.8e3))')i,lmfac(j),
     $              vL(j),vL1(j),vP(j),vP1(j),vR(j),vW(j),vLi(j)
            ENDDO
         ENDDO
         WRITE(out_unit,*)
      ENDIF

      ! start with DCON's total displacement vectors
      a = plas_indinvmats(resp_index,:,:)
      b = 0
      DO i=1,mpert
         b(i,i) = (chi1*twopi*ifac*(mfac(i)-nn*qlim))
      ENDDO
      ! convert to total flux
      a = MATMUL(CONJG(b),MATMUL(a,b))
      WRITE(out_unit,*)"Raw Matrices"
      WRITE(out_unit,*)" jac_type = ",jac_type
      WRITE(out_unit,*)" L = Vacuum Inductance"
      WRITE(out_unit,*)" Lambda = Inductance"
      WRITE(out_unit,*)" P = Permeability"
      WRITE(out_unit,*)" rho = Reluctance"
      WRITE(out_unit,*)
      WRITE(out_unit,'(2(1x,a4),12(1x,a12))')"m_i","m_j",
     $  "real(L)","imag(L)","real(Lambda)","imag(Lambda)",
     $  "real(P)","imag(P)","real(rho)","imag(rho)",
     $  "real(iLmda)","imag(iLmda)","real(W)","imag(W)"
      DO i=1,mpert
         DO j=1,mpert
            WRITE(out_unit,'(2(1x,I4),12(1x,es12.3))')mfac(i),mfac(j),
     $           surf_indmats(j,i),
     $           plas_indmats(resp_index,j,i),
     $           permeabmats(resp_index,j,i),
     $           reluctmats(resp_index,j,i),
     $           plas_indinvmats(resp_index,j,i),
     $           a(j,i)
         ENDDO
      ENDDO
      WRITE(out_unit,*)

      CALL ascii_close(out_unit)

      ! LOGAN - Eigenvectors on control surface in functions
      IF(fun_flag) THEN
        DO i=1,mpert
          CALL iscdftb(mfac,mpert,wtfun(:,i),mthsurf,wt(:,i))
          CALL iscdftb(mfac,mpert,ilfun(:,i),mthsurf,
     $      plas_indinvmats(resp_index,:,i))
          CALL iscdftb(mfac,mpert,rfun(:,i),mthsurf,
     $      reluctevmats(resp_index,:,i))
          CALL iscdftb(mfac,mpert,pfun(:,i),mthsurf,
     $      permeabevmats(resp_index,:,i))
        ENDDO
        DO i=0,mthsurf
          CALL bicube_eval(rzphi,psilim,theta(i),0)
          rfac=SQRT(rzphi%f(1))
          eta=twopi*(theta(i)+rzphi%f(2))
          r(i)=ro+rfac*COS(eta)
          z(i)=zo+rfac*SIN(eta)
        ENDDO
        IF(ascii_flag)THEN
           CALL ascii_open(out_unit,"gpec_response_fun_n"//
     $       TRIM(sn)//".out","UNKNOWN")
           WRITE(out_unit,*)"GPEC_RESPONSE_FUN: "//
     $       "Eigenmodes on control surface in functions"
           WRITE(out_unit,*)version
           WRITE(out_unit,*)" P = Permeability"
           WRITE(out_unit,*)" rho = Reluctance (power norm)"
           WRITE(out_unit,*)" W = Displacement energy"
           WRITE(out_unit,*)" F = Flux energy"
           WRITE(out_unit,*)
           WRITE(out_unit,'(1x,a12,I4)')"mthsurf =",mthsurf
           WRITE(out_unit,'(1x,a12,I4)')"mpert =",mpert
           WRITE(out_unit,*)
           WRITE(out_unit,*)"jac_type = "//jac_type
           WRITE(out_unit,*)
           DO j=1,mpert
             WRITE(out_unit,'(1x,2(a12,I4))')"isol =",j," n =",nn
             WRITE(out_unit,*)
             WRITE(out_unit,'(a16,9(1x,a16))') "r","z",
     $       "real(P)","imag(P)","real(rho)","imag(rho)",
     $       "real(W)","imag(W)","real(iLmda)","imag(iLmda)"
             DO i=0,mthsurf
               WRITE(out_unit,'(10(es17.8e3))') r(i),z(i),
     $           pfun(i,j),rfun(i,j),wtfun(i,j),ilfun(i,j)
             ENDDO
             WRITE(out_unit,*)
           ENDDO
           CALL ascii_close(out_unit)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_response
c-----------------------------------------------------------------------
c     subprogram 2. gpout_singcoup.
c     compute coupling between singular surfaces and external fields.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_singcoup(spot,nspot,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout,nspot
      REAL(r8), INTENT(IN) :: spot

      INTEGER :: i,j,itheta,ising,resnum,rsing,rpert,
     $     tmlow,tmhigh,tmpert,lwork,info
      REAL(r8) :: respsi,lpsi,rpsi,jarea,thetai
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun,wcfun,
     $     dphi,thetas,ones,jacs,jacfac

      COMPLEX(r8), DIMENSION(mpert) :: finmn,foutmn,
     $     fkaxmn,singflx_mn,ftnmn
      COMPLEX(r8), DIMENSION(lmpert) :: lftnmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      REAL(r8), DIMENSION(msing) :: area,j_c,w_c,shear

      COMPLEX(r8), DIMENSION(msing,mpert) :: deltas,delcurs,
     $     singcurs,islandhwids
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: convmat
      COMPLEX(r8), DIMENSION(msing,mpert,mpert) :: fsurfindmats
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: fldflxmn,bmn
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::  temp1, flxtofld,
     $     singbnoflxs, singbwp

      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE ::
     $     singcoup_out, localcoup_out,                       ! Coupling to external field in output coordinates
     $     singcoup_out_bvecs,localcoup_out_bvecs             ! Right singular vectors of power normalized coupling in output coordinates weighted by A_m^1/2/A^1/2
      REAL(r8), DIMENSION(:,:), ALLOCATABLE ::
     $     singcoup_out_vals, localcoup_out_vals              ! Singular values of energy normalized coupling in output coordinates

      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv,tmfac
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork,s
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: work
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: u,a,vt

      CHARACTER(72), DIMENSION(nsingcoup) :: titles, names
      CHARACTER(4), DIMENSION(nsingcoup) :: units
      CHARACTER(1), DIMENSION(nsingcoup) :: tags
      INTEGER :: idid, mdid, odid, cdid
      INTEGER, DIMENSION(nsingcoup) :: c_id, r_id, s_id,
     $     cl_id, rl_id, sl_id

      TYPE(spline_type) :: spl
      TYPE(cspline_type) :: fsp_sol
      COMPLEX(r8), DIMENSION(mpert) :: interpbwn
c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)
     $     "Computing coupling between total resonant fields and "//
     $     "external fields"
      IF(verbose) WRITE(*,*)"Computing surface inductances at "//
     $     "resonant surfaces"
      DO ising=1,msing
         resnum=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)
         area(ising)=0
         j_c(ising)=0
         w_c(ising)=0
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            rfacs(itheta)=rfac
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            sqreqb(itheta)=(sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $           /(twopi*r(itheta))**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
            wcfun(itheta)=SQRT(jac*sqreqb(itheta))
            area(ising)=area(ising)+jac*delpsi(itheta)/mthsurf
            j_c(ising)=j_c(ising)+jac*delpsi(itheta)
     $           *jcfun(itheta)/mthsurf
            w_c(ising)=w_c(ising)+0.5*wcfun(itheta)/mthsurf
         ENDDO
         area(ising)=area(ising)-jac*delpsi(mthsurf)/mthsurf
         j_c(ising)=j_c(ising)-jac*delpsi(mthsurf)*
     $        jcfun(mthsurf)/mthsurf
         w_c(ising)=w_c(ising)-0.5*wcfun(mthsurf)/mthsurf

         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0  
         shear(ising)=mfac(resnum)*sq%f1(4)/sq%f(4)**2

         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL gpvacuum_flxsurf(respsi)
         fsurfindmats(ising,:,:)=fsurf_indmats
         DEALLOCATE(fsurf_indev,fsurf_indmats)
      ENDDO
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      deltas=0
      delcurs=0
      singcurs=0
      CALL gpeq_alloc
      ALLOCATE(singbnoflxs(msing,mpert),singbwp(msing,mpert))
      DO i=1,mpert
         finmn=0
         finmn(i)=1.0                                 ! unit field
         CALL gpeq_weight(psilim,finmn,mfac,mpert,1)  ! converted to flux to get the solution
         IF (fixed_boundary_flag) THEN
            foutmn=finmn
         ELSE
            foutmn=MATMUL(permeabmats(resp_index,:,:),finmn)
         ENDIF
         singfac=mfac-nn*qlim
         edge_mn=foutmn/(chi1*singfac*twopi*ifac)
         edge_flag=.TRUE.
         CALL idcon_build(0,edge_mn)
c-----------------------------------------------------------------------
c     construct fsp_sol.
c-----------------------------------------------------------------------
         CALL gpeq_interp_singsurf(fsp_sol,spot,nspot)         
c-----------------------------------------------------------------------
c     evaluate delta/singular current/normal field/islands.
c-----------------------------------------------------------------------
         DO ising=1,msing
            resnum=NINT(singtype(ising)%q*nn)-mlow+1
            respsi=singtype(ising)%psifac
            lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
            CALL gpeq_sol(lpsi)
            lbwp1mn=bwp1_mn(resnum)
            
            rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1)) 
            CALL gpeq_sol(rpsi)
            rbwp1mn=bwp1_mn(resnum)

            deltas(ising,i)=(rbwp1mn - lbwp1mn) / (twopi * chi1)  ! like Park Phys. Plasmas 2007, but normalized by B.grad(theta)
            delcurs(ising,i)= (rbwp1mn - lbwp1mn) *
     $           j_c(ising) * ifac / (twopi * mfac(resnum))
            singcurs(ising,i)=-delcurs(ising,i)/ifac
            fkaxmn=0
            fkaxmn(resnum)=singcurs(ising,i)/(twopi*nn)
            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
c-----------------------------------------------------------------------
c     evaluation based on the interpolation.
c-----------------------------------------------------------------------
            CALL gpeq_interp_sol(fsp_sol,respsi,interpbwn)
            singbwp(ising,i)=interpbwn(resnum) / area(ising) ! flux normalized by area for units Tesla

            singbnoflxs(ising,i)=singflx_mn(resnum)/area(ising)  ! flux normalized by area for units Tesla
            islandhwids(ising,i)=4*singflx_mn(resnum)/
     $           (twopi*shear(ising)*singtype(ising)%q*chi1)     ! square of the width to be linear
         ENDDO
c-----------------------------------------------------------------------
c     deallocate fsp_sol.
c-----------------------------------------------------------------------
         CALL cspline_dealloc(fsp_sol)

         IF(verbose) WRITE(*,'(1x,a16,i4,a22,es10.3)')
     $        "poloidal mode =",mfac(i),", resonant coupling =",
     $        SUM(ABS(singbnoflxs(:,i)))/msing
      ENDDO
      CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     agregate coupling matrices
c-----------------------------------------------------------------------
      ALLOCATE(singcoup(nsingcoup,msing,mpert))  ! coupling field to external field b_x
      singcoup(1,:,:) = singbnoflxs        ! effective resonant area-normalized flux (Phi_r/A) shielded by current
      singcoup(2,:,:) = singcurs*twopi*nn  ! resonant current (amps)
      singcoup(3,:,:) = islandhwids        ! square of the penetrated island half-width in normalized poloidal flux
      singcoup(4,:,:) = singbwp            ! interpolated (penetrated) resonant field
      singcoup(5,:,:) = deltas             ! jump in db/dpsi modified from [Park, Phys. Plasmas 2007], using B.grad(theta) norm instread of B.grad(phi)
      singcoup_set  = .TRUE.
c-----------------------------------------------------------------------
c     convert coordinates for matrix on the plasma boundary.
c-----------------------------------------------------------------------
      IF ((jac_out /= jac_type).OR.(tout==0)) THEN
         IF(verbose) WRITE(*,*)"Converting coordinates"
         CALL spline_eval(sq,psilim,0)
         dphi=0
         
         CALL spline_alloc(spl,mthsurf,2)
         spl%xs=theta
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psilim,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            bpfac=psio*delpsi(itheta)/r(itheta)
            btfac=sq%f(1)/(twopi*r(itheta))
            bfac=SQRT(bpfac*bpfac+btfac*btfac)
            fac=r(itheta)**power_r/(bpfac**power_bp*bfac**power_b)
            spl%fs(itheta,1)=fac/(r(itheta)**rout*rfac**rcout)*
     $           bpfac**bpout*bfac**bout
            ! jacobian for coordinate angle at dcon angle
            spl%fs(itheta,2)=delpsi(itheta)*r(itheta)**rout*rfac**rcout/
     $           (bpfac**bpout*bfac**bout)  
            IF (tout .EQ. 0) THEN
               dphi(itheta)=rzphi%f(3)
            ENDIF
         ENDDO      
         
         CALL spline_fit(spl,"periodic")
         CALL spline_int(spl)

         ! coordinate angle at dcon angle
         thetas(:)=spl%fsi(:,1)/spl%fsi(mthsurf,1)
         DO itheta=0,mthsurf
            ! jacobian at coordinate angle
            thetai=issect(mthsurf,theta(:),thetas(:),theta(itheta))
            CALL spline_eval(spl,thetai,0)
            jacfac(itheta)=spl%f(2)
         ENDDO
         jarea=0
         DO itheta=0,mthsurf-1
            jarea=jarea+jacfac(itheta)/mthsurf
         ENDDO       

         CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     convert coordinates. 
c-----------------------------------------------------------------------
         ALLOCATE(fldflxmn(lmpert), fldflxmat(lmpert,lmpert),
     $        singcoup_out(nsingcoup,msing,lmpert), tmfac(lmpert))
         DO i=1,lmpert
            lftnmn=0
            lftnmn(i)=1.0
            ! compute given function in dcon angle
            DO itheta=0,mthsurf
               ftnfun(itheta)=0
               ftnfun(itheta)=
     $              lftnmn(i)*EXP(ifac*twopi*lmfac(i)*thetas(itheta))  
            ENDDO
            ! multiply toroidal factor for dcon angle
            IF (tout .EQ. 0) THEN
               ftnfun(:)=ftnfun(:)*EXP(-ifac*nn*dphi(:))
            ELSE
               ftnfun(:)=ftnfun(:)*
     $              EXP(-twopi*ifac*nn*sq%f(4)*(thetas(:)-theta(:)))
            ENDIF
            CALL iscdftf(mfac,mpert,ftnfun,mthsurf,ftnmn)
            convmat(:,i)=ftnmn
         
            CALL iscdftb(lmfac,lmpert,ftnfun,mthsurf,lftnmn)
            ftnfun(:)=ftnfun(:)*sqrt(jacfac(:))
            CALL iscdftf(lmfac,lmpert,ftnfun,mthsurf,fldflxmn)            
            fldflxmat(:,i)=fldflxmn/sqrt(jarea)
         ENDDO
         tmlow = lmlow
         tmhigh = lmhigh
         tmpert = lmpert
         tmfac= lmfac
         DO i=1,nsingcoup
            singcoup_out(i,:,:) = MATMUL(singcoup(i,:,:),convmat)
         ENDDO
      ELSE
         ALLOCATE(fldflxmn(mpert), fldflxmat(mpert,mpert),
     $        singcoup_out(nsingcoup,msing,mpert), tmfac(mpert))
         ones = (/(1.0,itheta=0,mthsurf)/)
         jarea = issurfint(ones,mthsurf,psilim,0,0)
         DO i=1,mpert
            ftnmn=0
            ftnmn(i)=1.0
            fldflxmn=ftnmn
            CALL gpeq_weight(psilim,fldflxmn,mfac,mpert,2)
            fldflxmat(:,i)=fldflxmn/sqrt(jarea)            
         ENDDO
         tmlow = mlow
         tmhigh = mhigh
         tmpert = mpert
         tmfac = mfac
         singcoup_out = singcoup
      ENDIF

c-----------------------------------------------------------------------
c     svd analysis in output coordinates.
c-----------------------------------------------------------------------
      lwork=3*tmpert
      ALLOCATE(s(msing),u(msing,msing),a(msing,tmpert),vt(msing,tmpert),
     $   work(lwork),rwork(5*msing),ipiv(tmpert),
     $   temp1(tmpert,tmpert), flxtofld(tmpert,tmpert),
     $   singcoup_out_vals(nsingcoup,msing),
     $   singcoup_out_vecs(nsingcoup,tmpert,msing),
     $   singcoup_out_bvecs(nsingcoup,tmpert,msing))

      ! get inverse of fldflxmat for converting to flux to power normalized field ((bA^1/2)_m / A^1/2)
      CALL iszhinv(fldflxmat,tmpert,flxtofld)

      ! calculate the output coodinat matrix SVD 
      ! re-normalize such that overlap dot products operate on unweighted fields
      DO i=1,nsingcoup
         work=0
         rwork=0
         s=0
         u=0
         vt=0
         a=MATMUL(singcoup_out(i, :, :), flxtofld)             ! coupling to power normalized external field
         CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $        work,lwork,rwork,info)                           ! note zgesvd vt is actually V^dagger
         singcoup_out_vals(i,:) = s                            ! singular values of power normalized coupling
         singcoup_out_vecs(i,:,:) = CONJG(TRANSPOSE(vt))       ! right singular vectors of power normalized coupling
         ! overlap = V.Phi_xe = V.W.b, and we record V.W so we can dot with b
         singcoup_out_bvecs(i,:,:) = MATMUL(CONJG(fldflxmat),  ! note fldflxmat is hermitian (W*=W^T)
     $                                  CONJG(TRANSPOSE(vt)))  ! vectors converted so overlap is with b
      ENDDO

      DEALLOCATE(s,u,a,vt,work,rwork,ipiv)

c-----------------------------------------------------------------------
c     svd analysis of local coupling in output coordinates.
c-----------------------------------------------------------------------
      IF (osing<msing) THEN
         lwork=3*tmpert
         ALLOCATE(s(osing),u(osing,osing),a(osing,tmpert),
     $      vt(osing,tmpert),work(lwork),rwork(5*osing),ipiv(tmpert),
     $      localcoup_out_vals(nsingcoup, osing),
     $      localcoup_out_vecs(nsingcoup,tmpert,osing),
     $      localcoup_out_bvecs(nsingcoup,tmpert,osing))

         DO i=1,nsingcoup
            work=0
            rwork=0
            s=0
            u=0
            vt=0
            a=MATMUL(singcoup_out(i,ol:ou,:),flxtofld)              ! coupling to power normalized external field
            CALL zgesvd('S','S',osing,tmpert,a,osing,s,u,osing,vt,osing,
     $           work,lwork,rwork,info)                             ! note zgesvd vt is actually V^dagger
            localcoup_out_vals(i,:) = s                             ! singular values of power normalized coupling
            localcoup_out_vecs(i,:,:) = CONJG(TRANSPOSE(vt))        ! right singular vectors of power normalized coupling
            ! overlap = V.Phi_xe = V.W.b, and we record V.W so we can dot with b
            localcoup_out_bvecs(i,:,:) = MATMUL(CONJG(fldflxmat),   ! note fldflxmat is hermitian (W*=W^T)
     $                                       CONJG(TRANSPOSE(vt)))  ! vectors converted so overlap is with b
         ENDDO
         DEALLOCATE(s,u,a,vt,work,rwork,ipiv)
      ENDIF

c-----------------------------------------------------------------------
c     write the outputs
c-----------------------------------------------------------------------
         titles(1) = "coupling matrix to effective resonant fields"
         titles(2) = "coupling matrix to singular currents"
         titles(3) = "coupling matrix to square of island half-widths"
         titles(4) = "coupling matrix to to penetrated resonant fields"
         titles(5) = "coupling matrix to unitless Delta"
         names(1) = "Effective energy normalized resonant flux"
         names(2) = "Resonant current"
         names(3) = "Square of island half-width in psi_n"
         names(4) = "Energy normalized penetrated flux"
         names(5) = "Delta"
         tags = (/'f', 's', 'i', 'p', 'd'/)
         units = (/'none', 'A/T ', '1/T ', 'none', '1/T '/)

c-----------------------------------------------------------------------
c     write the output coordinate coupling matrices.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_singcoup_matrix_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_SINGCOUP_MATRIX: Coupling matrices"//
     $        " between resonant field and external field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $        "jac_out = ",jac_out,"tmag_out =",tout
         WRITE(out_unit,'(4(1x,a12,I4))')
     $        "msing =",msing,"mpert =",tmpert,
     $        "mlow =",tmlow,"mhigh =",tmhigh
         WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $        "psilim =",psilim,"qlim =",qlim
         WRITE(out_unit,*)
         
         DO i=1,nsingcoup
            WRITE(out_unit,*) "The "//trim(titles(i))
            WRITE(out_unit,*)
            DO ising=1,msing
               WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es17.8e3)')
     $              "q =",singtype(ising)%q,
     $              "psi =",singtype(ising)%psifac
               WRITE(out_unit,*)
               WRITE(out_unit,'(1x,a4,2(1x,a16))')"m",
     $         "real(C_"//tags(i)//")","imag(C_"//tags(i)//")"
               DO j=1,tmpert
                  WRITE(out_unit,'(1x,I4,2(es17.8e3))') tmfac(j),
     $               REAL(singcoup_out(i,ising,j)),
     $               AIMAG(singcoup_out(i,ising,j))
               ENDDO
               WRITE(out_unit,*)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

c-----------------------------------------------------------------------
c     write the output coordinate coupling matrix SVDs.
c-----------------------------------------------------------------------
         CALL ascii_open(out_unit,"gpec_singcoup_svd_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_SINGCOUP_SVD: SVD analysis"//
     $        " for coupling matrices"
         WRITE(out_unit,*)"Note the right singular vectors have an "//
     $        "additional half-area weighting,"
         WRITE(out_unit,*)"such that their inner product with b_x"//
     $        " is the overlap field."
         WRITE(out_unit,*)"Dividing by Phi_xe**2 then gives the overlap"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a12,a8,1x,a12,I2)')
     $        "jac_out = ",jac_out,"tmag_out =",tmag_out
         WRITE(out_unit,'(4(1x,a12,I4))')
     $        "msing =",msing,"mpert =",tmpert,
     $        "mlow =",tmlow,"mhigh =",tmhigh
         WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $        "psilim =",psilim,"qlim =",qlim
         WRITE(out_unit,*)

         DO i=1,nsingcoup
            WRITE(out_unit,*) "Right singular vectors of the "
     $         //trim(titles(i))
            WRITE(out_unit,*)
            DO ising=1,msing
            WRITE(out_unit,'(1x,a6,I4,1x,a6,es17.8e3)')
     $           "mode =",ising,"s =",singcoup_out_vals(i,ising)
               WRITE(out_unit,*)
               WRITE(out_unit,'(1x,a4,2(1x,a16))') "m",
     $            "real(V_"//tags(i)//")", "imag(V_"//tags(i)//")"
               DO j=1,tmpert
                  WRITE(out_unit,'(1x,I4,2(es17.8e3))') tmfac(j),
     $               REAL(singcoup_out_bvecs(i,j,ising)),
     $               AIMAG(singcoup_out_bvecs(i,j,ising))
               ENDDO
               WRITE(out_unit,*)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

c-----------------------------------------------------------------------
c     write the output coordinate local coupling matrix SVDs.
c-----------------------------------------------------------------------
         IF (osing<msing) THEN
            CALL ascii_open(out_unit,"gpec_singcoup_svd_local_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_SINGCOUP_LOCAL_SVD: SVD"//
     $           " analysis for local coupling matrices"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a12,a8,1x,a12,I2)')
     $           "jac_out = ",jac_out,"tmag_out =",tmag_out
            WRITE(out_unit,'(4(1x,a12,I4))')
     $           "msing =",osing,"mpert =",tmpert,
     $           "mlow =",tmlow,"mhigh =",tmhigh
            WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $           "psilim =",psilim,"qlim =",qlim
            WRITE(out_unit,*)

            DO i=1,nsingcoup
               WRITE(out_unit,*) "Right singular vectors of the local "
     $            //trim(titles(i))
               WRITE(out_unit,*)
               DO ising=1,osing
               WRITE(out_unit,'(1x,a6,I4,1x,a6,es17.8e3)')
     $              "mode =",ising,"s =",localcoup_out_vals(i,ising)
                  WRITE(out_unit,*)
                  WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
                  DO j=1,tmpert
                     WRITE(out_unit,'(1x,I4,2(es17.8e3))') tmfac(j),
     $                  REAL(localcoup_out_bvecs(i,j,ising)),
     $                  AIMAG(localcoup_out_bvecs(i,j,ising))
                  ENDDO
                  WRITE(out_unit,*)
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF
      ENDIF

c-----------------------------------------------------------------------
c     Write the same output coordinate data to netcdf
c-----------------------------------------------------------------------
      IF(netcdf_flag)THEN
         ! append netcdf file
         CALL check( nf90_open(mncfile,nf90_write,mncid) )

         ! references to pre-defined dimensions
         CALL check( nf90_inq_dimid(mncid,"i",idid) )
         CALL check( nf90_inq_dimid(mncid,"m_out",mdid) )
         CALL check( nf90_inq_dimid(mncid,"mode_C",cdid) )
         CALL check( nf90_inq_dimid(mncid,"mode_C_local",odid) )

         ! start definitions
         CALL check( nf90_redef(mncid))
         DO i=1,nsingcoup
            CALL check( nf90_def_var(mncid,"C_"//tags(i)//"_x_out",
     $           nf90_double, (/mdid,cdid,idid/), c_id(i)) )
            CALL check( nf90_put_att(mncid,c_id(i),"long_name",
     $           trim(names(i))//"coupling to external flux") )
            CALL check( nf90_put_att(mncid,c_id(i),"units",units(i)) )
            CALL check( nf90_def_var(mncid,
     $           "C_"//tags(i)//"_xb_out_singvec",
     $           nf90_double,(/mdid,cdid,idid/),r_id(i)) )
            CALL check( nf90_put_att(mncid,r_id(i),"long_name",
     $           trim(names(i))//" coupling right-singular vectors") )
            CALL check( nf90_def_var(mncid,
     $           "C_"//tags(i)//"_xe_out_singval",
     $           nf90_double,(/cdid/),s_id(i)) )
            CALL check( nf90_put_att(mncid,s_id(i),"long_name",
     $           trim(names(i))//" coupling singular values") )
            CALL check( nf90_put_att(mncid,s_id(i),"units",units(i)) )
            IF(osing<msing)THEN
               CALL check( nf90_def_var(mncid,"C_l"//tags(i)//"_x_out",
     $              nf90_double, (/mdid,odid,idid/), cl_id(i)) )
               CALL check( nf90_put_att(mncid,cl_id(i),"long_name",
     $           "Local "//trim(names(i))//"coupling to external flux"))
               CALL check(nf90_put_att(mncid,cl_id(i),"units",units(i)))
               CALL check( nf90_def_var(mncid,
     $              "C_l"//tags(i)//"_xb_out_singvec",
     $              nf90_double,(/mdid,odid,idid/),rl_id(i)) )
               CALL check( nf90_put_att(mncid,rl_id(i),"long_name",
     $          "Local "//trim(names(i))//" coupling right singular"))
               CALL check( nf90_def_var(mncid,
     $              "C_l"//tags(i)//"_xe_out_singval",
     $              nf90_double,(/odid/),sl_id(i)) )
               CALL check( nf90_put_att(mncid,sl_id(i),"long_name",
     $          "Local "//trim(names(i))//" coupling singular values"))
               CALL check(nf90_put_att(mncid,sl_id(i),"units",units(i)))
            ENDIF 
         ENDDO

         ! end definitions
         CALL check( nf90_enddef(mncid) )

         ! write variables
         DO i=1,nsingcoup
            CALL check( nf90_put_var(mncid,c_id(i),RESHAPE(
     $           (/REAL(singcoup_out(i,:,:)),
     $             AIMAG(singcoup_out(i,:,:))/),(/tmpert,msing,2/))))
            CALL check( nf90_put_var(mncid,r_id(i),RESHAPE(
     $           (/REAL(singcoup_out_bvecs(i,:,:)),
     $             AIMAG(singcoup_out_bvecs(i,:,:))/),
     $           (/tmpert,msing,2/))))
            CALL check( nf90_put_var(mncid,s_id(i),
     $           singcoup_out_vals(i,:)) )
            IF(osing<msing)THEN
               CALL check( nf90_put_var(mncid,cl_id(i),RESHAPE(
     $              (/REAL(localcoup_out(i,:,:)),
     $              AIMAG(localcoup_out(i,:,:))/),(/tmpert,osing,2/))))
               CALL check( nf90_put_var(mncid,rl_id(i),RESHAPE(
     $              (/REAL(localcoup_out_bvecs(i,:,:)),
     $              AIMAG(localcoup_out_bvecs(i,:,:))/),
     $              (/tmpert,osing,2/))))
               CALL check( nf90_put_var(mncid,sl_id(i),
     $              localcoup_out_vals(i,:)))
            ENDIF
         ENDDO

         ! close the file
         CALL check( nf90_close(mncid) )
      ENDIF

c-----------------------------------------------------------------------
c     Deallocate local variables.
c-----------------------------------------------------------------------
      DEALLOCATE(singcoup_out, singcoup_out_bvecs, singcoup_out_vals,
     $     tmfac, fldflxmn, temp1, flxtofld)
      IF (osing<msing) THEN
         DEALLOCATE(localcoup_out_bvecs, localcoup_out_vals)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_singcoup
c-----------------------------------------------------------------------
c     subprogram 3. gpout_control
c     calculate response from external field on the control surface.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_control(ifile,finmn,foutmn,xspmn,
     $     rin,bpin,bin,rcin,tin,jin,rout,bpout,bout,rcout,tout,jout,
     $     filter_types,filter_modes,filter_out)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: ifile
      INTEGER, INTENT(IN) :: rin,bpin,bin,rcin,tin,jin,
     $     rout,bpout,bout,rcout,tout,jout,filter_modes
      LOGICAL, INTENT(IN) :: filter_out
      CHARACTER(len=*), INTENT(IN) :: filter_types
      
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: foutmn,xspmn

      INTEGER :: i,j,i1,i2,i3,itheta,mpert_in
      INTEGER, DIMENSION(:), ALLOCATABLE :: mfac_in
      REAL(r8) :: vengy,sengy,pengy,scale,norm,jarea
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,jacs,rvecs,
     $     zvecs,units

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,xinmn,xoutmn,tempmn
      COMPLEX(r8), DIMENSION(lmpert) :: cinmn,coutmn,templ
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,xinfun,xoutfun
      COMPLEX(r8), DIMENSION(mpert,mpert) :: sqrtamat,j2mat,tempmm
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: tempml
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: cawmn

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: dcosmn,dsinmn
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rawmn
      
      INTEGER :: i_id,m_id,modid,t_id,r_id,z_id,rn_id,zn_id,p_id,
     $    x_id,xx_id,xm_id,xxm_id,bm_id,bxm_id,b_id,bx_id,jo_id,j2_id,
     $    mpdid,mp_id
c-----------------------------------------------------------------------
c     basic definitions
c-----------------------------------------------------------------------
      mpert_in = ABS(mmax)+ABS(mmin)+1
      ALLOCATE(mfac_in(mpert_in),cawmn(mpert_in))
      mfac_in = (/(i,i=mmin,mmax)/)
c-----------------------------------------------------------------------
c     check data_type and read data.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      scale=1.0
      cawmn=0
      tempmn=0
      IF (data_flag) THEN
         IF(verbose) WRITE(*,*)"Reading external fields on the "//
     $      "control surface"
         ALLOCATE(dcosmn(mmin:mmax,nmin:nmax),
     $        dsinmn(mmin:mmax,nmin:nmax),rawmn(mmin:mmax,nmin:nmax))
         IF (data_type=="surfmn") THEN
            i1=mmax
            i2=mmin
            i3=-1
            scale=1e-4
         ELSE
            i1=mmin
            i2=mmax
            i3=1
         ENDIF
         CALL ascii_open(in_unit,ifile,"old")
                        
         DO i=i1,i2,i3
            IF (data_type=="surfmn") THEN
               READ(in_unit,'(1x,25f12.6)')(dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,'(1x,25f12.6)')(dsinmn(i,j),j=nmin,nmax)
            ELSE IF (data_type=="vac3d") THEN
               READ(in_unit,'(11(1x,e15.8))')(dcosmn(i,j),j=nmin,nmax)
               READ(in_unit,'(11(1x,e15.8))')(dsinmn(i,j),j=nmin,nmax)
            ELSE
               WRITE(message,'(a)')"Can't recognize data format"
               CALL gpec_stop(message)
            ENDIF            
         ENDDO
         
         CALL ascii_close(in_unit)
         rawmn=dcosmn+ifac*dsinmn
         cawmn=rawmn(:,nn)
         DEALLOCATE(dcosmn,dsinmn,rawmn)
      ELSE IF (harmonic_flag) THEN
         DO i=-hmnum,hmnum
            IF ((-mmin+i+1>=1).AND.(-mmin+i+1<=mpert_in)) THEN
               cawmn(-mmin+i+1)=cosmn(i)+ifac*sinmn(i)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      IF (data_flag .OR. harmonic_flag) THEN
         CALL gpeq_fcoords(psilim,cawmn,mfac_in,mpert_in,
     $        rin,bpin,bin,rcin,tin,jin)
         binmn=0
         DO i=1,mpert_in
            IF ((mmin-mlow+i>=1).AND.(mmin-mlow+i<=mpert)) THEN
               binmn(mmin-mlow+i)=cawmn(i)
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     convert to field if displacement is given.
c-----------------------------------------------------------------------
         IF (displacement_flag) THEN
            CALL gpeq_weight(psilim,binmn,mfac,mpert,5)
            binmn=twopi*ifac*chi1*(mfac-nn*qlim)*binmn
            CALL gpeq_weight(psilim,binmn,mfac,mpert,0)
         ENDIF 
         binmn=binmn*scale
         tempmn=binmn
         CALL gpeq_weight(psilim,tempmn,mfac,mpert,1)
      ENDIF
      finmn=finmn+tempmn
c-----------------------------------------------------------------------
c     filter external flux
c-----------------------------------------------------------------------
      CALL gpout_control_filter(finmn,foutmn,filter_types,filter_modes,
     $           rout,bpout,bout,rcout,tout,filter_out)
c-----------------------------------------------------------------------
c     get plasma response on the control surface.
c-----------------------------------------------------------------------
      IF (fixed_boundary_flag) THEN
         foutmn=finmn
      ELSE
         foutmn=MATMUL(permeabmats(resp_index,:,:),finmn)
      ENDIF
      xspmn=foutmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      binmn=finmn
      boutmn=foutmn
      CALL gpeq_weight(psilim,binmn,mfac,mpert,0)
      CALL gpeq_weight(psilim,boutmn,mfac,mpert,0)
      xinmn=finmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      xoutmn=xspmn
      CALL gpeq_weight(psilim,xinmn,mfac,mpert,4)
      CALL gpeq_weight(psilim,xoutmn,mfac,mpert,4)
c-----------------------------------------------------------------------
c     compute perturbed energy.
c-----------------------------------------------------------------------
      vy = SUM(CONJG(finmn)*MATMUL(surf_indinvmats,finmn))/4.0
      sy = SUM(CONJG(foutmn)*MATMUL(surf_indinvmats,foutmn))/4.0
      vengy = REAL(vy)
      sengy = REAL(sy)

      py = SUM(CONJG(foutmn)*MATMUL(plas_indinvmats(resp_index,:,:),
     $     foutmn))/4.0
      pengy = REAL(py)
      IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $  "Required energy to perturb vacuum = ",sengy
      IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $  "Required energy to perturb plasma = ",REAL(py)
      IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $  "Required torque to perturb plasma = ",-2.0*nn*AIMAG(py)
      IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $  "Amplification factor = ",sengy/pengy

      ! log results with harvest
      norm=SQRT(ABS(DOT_PRODUCT(finmn,finmn)))
      ierr=set_harvest_payload_dbl(hlog,"MODPhi_x"//nul,norm)
      norm=SQRT(ABS(DOT_PRODUCT(foutmn,foutmn)))
      ierr=set_harvest_payload_dbl(hlog,"MODPhi"//nul,norm)
      ierr=set_harvest_payload_dbl(hlog,"energy_vacuum"//nul,vengy)
      ierr=set_harvest_payload_dbl(hlog,"energy_surface"//nul,sengy)
      ierr=set_harvest_payload_dbl(hlog,"energy_plasma"//nul,pengy)
      ierr=set_harvest_payload_dbl(hlog,
     $     "toroidal_torque"//nul,-2.0*nn*AIMAG(py))
      ierr=set_harvest_payload_dbl(hlog,
     $     "amplification"//nul,sengy/pengy)
c-----------------------------------------------------------------------
c     calculate the scalar surface area for normalizations.
c-----------------------------------------------------------------------
      units = 1
      jarea = issurfint(units,mthsurf,psilim,0,0)
c-----------------------------------------------------------------------
c     calculate dimensionless half area weighting matrix W
c-----------------------------------------------------------------------
!     W_m,m' = int{sqrt(J|delpsi|)exp[-i*(m-m')t]dt}/
!              int{sqrt(J|delpsi|)dt}
      sqrtamat = 0
      DO i=1,mpert
         sqrtamat(i,i) = 1.0
         CALL gpeq_weight(psilim,sqrtamat(:,i),mfac,mpert,2) ! A^1/2
      ENDDO
      j2mat = sqrtamat/sqrt(jarea)
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_control_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_CONTROL: "//
     $        "Plasma response for an external perturbation on the "//
     $        "control surface"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a12,I4)')"mpert =",mpert
         WRITE(out_unit,'(1x,a16,es17.8e3)')"vacuum energy =",vengy
         WRITE(out_unit,'(1x,a16,es17.8e3)')"surface energy =",sengy
         WRITE(out_unit,'(1x,a16,es17.8e3)')"plasma energy =",pengy
         WRITE(out_unit,'(1x,a18,es17.8e3)')"toroidal torque =",
     $        -2*nn*AIMAG(py)
         WRITE(out_unit,*)
         WRITE(out_unit,*)"jac_type = "//jac_type
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,8(1x,a16))') "m",
     $        "real(bin)","imag(bin)","real(bout)","imag(bout)",
     $        "real(Phi^x)","imag(Phi^x)","real(Phi)","imag(Phi)"
         DO i=1,mpert
            WRITE(out_unit,'(1x,I4,8(es17.8e3))')mfac(i),
     $           REAL(binmn(i)),AIMAG(binmn(i)),
     $           REAL(boutmn(i)),AIMAG(boutmn(i)),
     $           REAL(finmn(i)),AIMAG(finmn(i)),
     $           REAL(foutmn(i)),AIMAG(foutmn(i))
         ENDDO
         WRITE(out_unit,*)

         IF ((jac_in /= jac_type).OR.(tin==0).OR.(jin/=0)) THEN
            cinmn=0
            coutmn=0
            DO i=1,mpert
               IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
                  cinmn(mlow-lmlow+i)=binmn(i)
                  coutmn(mlow-lmlow+i)=boutmn(i)
               ENDIF
            ENDDO
            CALL gpeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $           rin,bpin,bin,rcin,tin,jin)
            CALL gpeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $           rin,bpin,bin,rcin,tin,jin)
            WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')
     $           "jac_in = ",jac_in,"jsurf_in =",jin,"tmag_in =",tin
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a4,4(1x,a16))')"m",
     $           "real(bin)","imag(bin)","real(bout)","imag(bout)"
            DO i=1,lmpert
               WRITE(out_unit,'(1x,I4,4(es17.8e3))')lmfac(i),
     $              REAL(cinmn(i)),AIMAG(cinmn(i)),
     $              REAL(coutmn(i)),AIMAG(coutmn(i))
            ENDDO
            WRITE(out_unit,*)
         ENDIF
      ENDIF

      IF ((jac_out /= jac_type).OR.(tout==0).OR.(jout/=0)) THEN
         cinmn=0
         coutmn=0
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               cinmn(mlow-lmlow+i)=binmn(i)
               coutmn(mlow-lmlow+i)=boutmn(i)
            ENDIF
         ENDDO
         CALL gpeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout) 
         CALL gpeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout)
         IF(ascii_flag)THEN
            WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_out = ",
     $           jac_out,"jsurf_out =",jout,"tmag_out =",tout
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)",
     $           "imag(bin)","real(bout)","imag(bout)"
            DO i=1,lmpert
               WRITE(out_unit,'(1x,I4,6(es17.8e3))')lmfac(i),
     $              REAL(cinmn(i)),AIMAG(cinmn(i)),
     $              REAL(coutmn(i)),AIMAG(coutmn(i))
            ENDDO
            WRITE(out_unit,*)
         ENDIF
         IF (singcoup_set) THEN
            ALLOCATE(sbno_mn(lmpert),sbno_fun(0:mthsurf))
            sbno_mn=MATMUL(fldflxmat,cinmn)
            CALL iscdftb(lmfac,lmpert,sbno_fun,mthsurf,sbno_mn)
            sbno_mn=cinmn
         ENDIF
      ELSE
         IF (singcoup_set) THEN
            ALLOCATE(sbno_mn(mpert),sbno_fun(0:mthsurf))
            sbno_mn=MATMUL(fldflxmat,binmn)
            CALL iscdftb(mfac,mpert,sbno_fun,mthsurf,sbno_mn)
            sbno_mn=binmn
         ENDIF
      ENDIF

      IF(ascii_flag) CALL ascii_close(out_unit)

      ! Convert to output coords
      DO i=1,mpert
         IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
            templ = 0
            templ(mlow-lmlow+i) = 1.0
            IF((jac_out /= jac_type).OR.(tout==0)) CALL gpeq_bcoords(
     $        psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,jout)
            coordmat(:,i) = templ
         ENDIF
      ENDDO

      IF(netcdf_flag)THEN
         CALL check( nf90_open(mncfile,nf90_write,mncid) )
         CALL check( nf90_inq_dimid(mncid,"i",i_id) )
         CALL check( nf90_inq_dimid(mncid,"m",m_id) )
         CALL check( nf90_inq_dimid(mncid,"m_out",modid) )
         CALL check( nf90_redef(mncid))
         CALL check( nf90_def_dim(mncid, "m_prime", mpert, mpdid) )
         CALL check( nf90_def_var(mncid, "m_prime", nf90_int, mpdid,
     $      mp_id))
         CALL check( nf90_put_att(mncid,nf90_global,
     $               "energy_vacuum",vengy) )
         CALL check( nf90_put_att(mncid,nf90_global,
     $               "energy_surface",sengy) )
         CALL check( nf90_put_att(mncid,nf90_global,
     $               "energy_plasma",pengy) )
         CALL check( nf90_put_att(mncid,nf90_global,
     $               "toroidal_torque",-2.0*nn*AIMAG(py)) )
         CALL check( nf90_put_att(mncid,nf90_global,
     $               "area",jarea) )
         CALL check( nf90_def_var(mncid,"b_n",nf90_double,
     $               (/m_id,i_id/),bm_id) )
         CALL check( nf90_put_att(mncid,bm_id,"units","Tesla") )
         CALL check( nf90_put_att(mncid,bm_id,"long_name",
     $               "Normal field") )
         CALL check( nf90_def_var(mncid,"b_n_x",nf90_double,
     $               (/m_id,i_id/),bxm_id) )
         CALL check( nf90_put_att(mncid,bxm_id,"units","Tesla") )
         CALL check( nf90_put_att(mncid,bxm_id,"long_name",
     $               "Externally applied normal field") )
         CALL check( nf90_def_var(mncid,"xi_n",nf90_double,
     $               (/m_id,i_id/),xm_id) )
         CALL check( nf90_put_att(mncid,xm_id,"units","m") )
         CALL check( nf90_put_att(mncid,xm_id,"long_name",
     $               "Normal displacement") )
         CALL check( nf90_def_var(mncid,"xi_n_x",nf90_double,
     $               (/m_id,i_id/),xxm_id) )
         CALL check( nf90_put_att(mncid,xxm_id,"units","m") )
         CALL check( nf90_put_att(mncid,xxm_id,"long_name",
     $               "Externally applied normal displacement") )
         CALL check( nf90_def_var(mncid,"J_out",nf90_double,
     $               (/m_id,modid,i_id/),jo_id) )
         CALL check( nf90_put_att(mncid,jo_id,"long_name",
     $      "Transform to jac_out, tmag_out, jsurf_out") )
         CALL check( nf90_put_att(mncid,jo_id,"jacobian",jac_type) )
         CALL check( nf90_put_att(mncid,jo_id,"jac_out",jac_out) )
         CALL check( nf90_put_att(mncid,jo_id,"tmag_out",tmag_out) )
         CALL check( nf90_put_att(mncid,jo_id,"jsurf_out",jsurf_out) )
         CALL check( nf90_def_var(mncid,"J_surf_2",nf90_double,
     $               (/m_id,mpdid,i_id/),j2_id) )
         CALL check( nf90_put_att(mncid,j2_id,"long_name",
     $      "Transform to jsurf_out 2") )
         CALL check( nf90_enddef(mncid) )
         CALL check( nf90_put_var(mncid,bxm_id,RESHAPE((/REAL(binmn),
     $             AIMAG(binmn)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,bm_id,RESHAPE((/REAL(boutmn),
     $             AIMAG(boutmn)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,xxm_id,RESHAPE((/REAL(xinmn),
     $             AIMAG(xinmn)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,xm_id,RESHAPE((/REAL(xoutmn),
     $             AIMAG(xoutmn)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,mp_id,mfac) )
         tempml = TRANSPOSE(coordmat)
         CALL check( nf90_put_var(mncid,jo_id,RESHAPE((/REAL(tempml),
     $             AIMAG(tempml)/),(/mpert,lmpert,2/))) )
         tempmm = TRANSPOSE(j2mat)
         CALL check( nf90_put_var(mncid,j2_id,RESHAPE((/REAL(tempmm),
     $             AIMAG(tempmm)/),(/mpert,mpert,2/))) )
         CALL check( nf90_close(mncid) )
      ENDIF

      CALL gpeq_bcoords(psilim,binmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL gpeq_bcoords(psilim,boutmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL gpeq_bcoords(psilim,xinmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)
      CALL gpeq_bcoords(psilim,xoutmn,mfac,mpert,
     $     power_r,power_bp,power_b,0,0,0)

      CALL iscdftb(mfac,mpert,binfun,mthsurf,binmn)
      CALL iscdftb(mfac,mpert,boutfun,mthsurf,boutmn)
      CALL iscdftb(mfac,mpert,xinfun,mthsurf,xinmn)
      CALL iscdftb(mfac,mpert,xoutfun,mthsurf,xoutmn)

      DO itheta=0,mthsurf
         CALL bicube_eval(rzphi,psilim,theta(itheta),0)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         dphi(itheta)=rzphi%f(3)
         jacs(itheta)=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jacs(itheta)
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jacs(itheta))
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
         rvecs(itheta)=
     $        (cos(eta)*w(1,1)-sin(eta)*w(1,2))/delpsi(itheta)
         zvecs(itheta)=
     $        (sin(eta)*w(1,1)+cos(eta)*w(1,2))/delpsi(itheta)
      ENDDO

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_control_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_CONTROL_FUN: "//
     $        "Plasma response for an external perturbation on the "//
     $        "control surface in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,I4)')"mthsurf =",mthsurf
         WRITE(out_unit,'(1x,a16,es17.8e3)')"vacuum energy =",vengy
         WRITE(out_unit,'(1x,a16,es17.8e3)')"surface energy =",sengy
         WRITE(out_unit,'(1x,a16,es17.8e3)')"plasma energy =",pengy
         WRITE(out_unit,'(1x,a16,es17.8e3)')"toroidal torque =",
     $        -2.0*nn*AIMAG(py)
         WRITE(out_unit,*)
         WRITE(out_unit,*)"jac_type = "//jac_type
         WRITE(out_unit,*)
         WRITE(out_unit,'(12(1x,a16))')
     $        "r","z","theta","dphi",
     $        "real(xin)","imag(xin)","real(xout)","imag(xout)",
     $        "real(bin)","imag(bin)","real(bout)","imag(bout)"
         DO itheta=0,mthsurf
            WRITE(out_unit,'(12(es17.8e3))')r(itheta),z(itheta),
     $           theta(itheta),dphi(itheta),
     $           REAL(xinfun(itheta)),-helicity*AIMAG(xinfun(itheta)),
     $           REAL(xoutfun(itheta)),-helicity*AIMAG(xoutfun(itheta)),
     $           REAL(binfun(itheta)),-helicity*AIMAG(binfun(itheta)),
     $           REAL(boutfun(itheta)),-helicity*AIMAG(boutfun(itheta))
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(mncfile,nf90_write,mncid) )
         CALL check( nf90_inq_dimid(mncid,"i",i_id) )
         CALL check( nf90_inq_dimid(mncid,"theta",t_id) )
         CALL check( nf90_redef(mncid))
         CALL check( nf90_def_var(mncid, "xi_n_fun", nf90_double,
     $                    (/t_id,i_id/),x_id) )
         CALL check( nf90_put_att(mncid,x_id,"long_name",
     $               "Normal displacement") )
         CALL check( nf90_put_att(mncid,x_id,"units","m") )
         CALL check( nf90_def_var(mncid, "xi_n_x_fun", nf90_double,
     $                    (/t_id,i_id/),xx_id) )
         CALL check( nf90_put_att(mncid,xx_id,"long_name",
     $               "Externally applied displacement") )
         CALL check( nf90_put_att(mncid,xx_id,"units","m") )
         CALL check( nf90_def_var(mncid, "b_n_fun", nf90_double,
     $                    (/t_id,i_id/),b_id) )
         CALL check( nf90_put_att(mncid,b_id,"long_name",
     $               "Normal field") )
         CALL check( nf90_put_att(mncid,b_id,"units","Tesla") )
         CALL check( nf90_def_var(mncid, "b_n_x_fun", nf90_double,
     $                    (/t_id,i_id/),bx_id) )
         CALL check( nf90_put_att(mncid,bx_id,"long_name",
     $               "Externally applied field") )
         CALL check( nf90_put_att(mncid,bx_id,"units","Tesla") )
         CALL check( nf90_def_var(mncid, "R", nf90_double,t_id,r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $               "Major radius") )
         CALL check( nf90_put_att(mncid,r_id,"units","m") )
         CALL check( nf90_def_var(mncid, "z", nf90_double,t_id,z_id) )
         CALL check( nf90_put_att(mncid,z_id,"long_name",
     $               "Vertical position") )
         CALL check( nf90_put_att(mncid,z_id,"units","m") )
         CALL check( nf90_def_var(mncid, "R_n", nf90_double,t_id,rn_id))
         CALL check( nf90_put_att(mncid,rn_id,"long_name",
     $               "Major radius component of normal unit vector") )
         CALL check( nf90_put_att(mncid,rn_id,"units","m") )
         CALL check( nf90_def_var(mncid, "z_n", nf90_double,t_id,zn_id))
         CALL check( nf90_put_att(mncid,zn_id,"long_name",
     $               "Vertical component of normal unit vector") )
         CALL check( nf90_put_att(mncid,zn_id,"units","m") )
         CALL check( nf90_def_var(mncid, "delta_phi", nf90_double,
     $               t_id,p_id))
         CALL check( nf90_put_att(mncid,p_id,"long_name",
     $               "Toroidal angle - magnetic angle") )
         CALL check( nf90_enddef(mncid) )
         CALL check( nf90_put_var(mncid,xx_id,RESHAPE((/REAL(xinfun),
     $             -helicity*AIMAG(xinfun)/),(/mthsurf+1,2/))) )
         CALL check( nf90_put_var(mncid,x_id,RESHAPE((/REAL(xoutfun),
     $             -helicity*AIMAG(xoutfun)/),(/mthsurf+1,2/))) )
         CALL check( nf90_put_var(mncid,bx_id,RESHAPE((/REAL(binfun),
     $             -helicity*AIMAG(binfun)/),(/mthsurf+1,2/))) )
         CALL check( nf90_put_var(mncid,b_id,RESHAPE((/REAL(boutfun),
     $             -helicity*AIMAG(boutfun)/),(/mthsurf+1,2/))) )
         CALL check( nf90_put_var(mncid,r_id,r) )
         CALL check( nf90_put_var(mncid,z_id,z) )
         CALL check( nf90_put_var(mncid,rn_id,rvecs) )
         CALL check( nf90_put_var(mncid,zn_id,zvecs) )
         CALL check( nf90_put_var(mncid,p_id,dphi) )
         CALL check( nf90_close(mncid) )
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_control
c-----------------------------------------------------------------------
c     subprogram 4. gpout_singfld.
c     compute current and field on rational surfaces.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_singfld(egnum,xspmn,spot,nspot)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nspot
      REAL(r8), INTENT(IN) :: spot
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i_id,q_id,m_id,p_id,c_id,w_id,k_id,n_id,d_id,
     $           pp_id,cp_id,wp_id,np_id,dp_id

      INTEGER :: itheta,ising,icoup
      REAL(r8) :: respsi,lpsi,rpsi,shear,hdist,sbnosurf
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: area,j_c,aq,asingflx
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: fkaxmn

      REAL(r8), DIMENSION(msing) :: island_hwidth,chirikov
      REAL(r8), DIMENSION(nsingcoup,msing) :: op
      COMPLEX(r8), DIMENSION(msing) :: delta,delcur,singcur,
     $     singflx,singbwp
      COMPLEX(r8), DIMENSION(nsingcoup, msing) :: olap
      COMPLEX(r8), DIMENSION(mpert,msing) :: singflx_mn

      TYPE(spline_type) :: spl
      TYPE(cspline_type) :: fsp_sol
      COMPLEX(r8), DIMENSION(mpert) :: interpbwn
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing total resonant fields"
      CALL gpeq_alloc
      CALL idcon_build(egnum,xspmn)

      CALL gpeq_interp_singsurf(fsp_sol,spot,nspot)

      IF (vsbrzphi_flag) ALLOCATE(singbno_mn(mpert,msing))
c-----------------------------------------------------------------------
c     evaluate delta and singular currents.
c-----------------------------------------------------------------------
!     j_c is j_c/(chi1*sq%f(4))
      DO ising=1,msing
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)
         area(ising)=0
         j_c(ising)=0
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(eta)
            z(itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            sqreqb(itheta)=(sq%f(1)**2+chi1**2*delpsi(itheta)**2)
     $           /(twopi*r(itheta))**2
            jcfun(itheta)=sqreqb(itheta)/(delpsi(itheta)**3)
            area(ising)=area(ising)+jac*delpsi(itheta)/mthsurf
            j_c(ising)=j_c(ising)+jac*delpsi(itheta)
     $           *jcfun(itheta)/mthsurf
         ENDDO
         area(ising)=area(ising)-jac*delpsi(mthsurf)/mthsurf
         j_c(ising)=j_c(ising)-jac*delpsi(mthsurf)*
     $        jcfun(mthsurf)/mthsurf

         j_c(ising)=1.0/j_c(ising)*chi1**2*sq%f(4)/mu0
         shear=mfac(resnum(ising))*sq%f1(4)/sq%f(4)**2

         CALL gpeq_interp_sol(fsp_sol,respsi,interpbwn)
         singbwp(ising)=interpbwn(resnum(ising))

         lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
         CALL gpeq_sol(lpsi)
         lbwp1mn=bwp1_mn(resnum(ising))
         
         rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1))
         CALL gpeq_sol(rpsi)
         rbwp1mn=bwp1_mn(resnum(ising))

         delta(ising) = (rbwp1mn - lbwp1mn) / (twopi * chi1)
         delcur(ising)= (rbwp1mn - lbwp1mn) * j_c(ising) * ifac /
     $        (twopi*mfac(resnum(ising)))
         singcur(ising)=-delcur(ising)/ifac

         fkaxmn=0
         fkaxmn(resnum(ising))=singcur(ising)/(twopi*nn)
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL gpvacuum_flxsurf(respsi)
         singflx_mn(:,ising)=MATMUL(fsurf_indmats,fkaxmn)
         DEALLOCATE(fsurf_indmats,fsurf_indev)
c-----------------------------------------------------------------------
c     compute half-width of magnetic island.
c-----------------------------------------------------------------------
         island_hwidth(ising)=
     $        SQRT(ABS(4*singflx_mn(resnum(ising),ising)/
     $        (twopi*shear*sq%f(4)*chi1)))
c-----------------------------------------------------------------------
c     compute coordinate-independent resonant field.
c----------------------------------------------------------------------- 
         IF (vsbrzphi_flag) THEN
            singbno_mn(:,ising)=-singflx_mn(:,ising)
!            CALL gpeq_weight(respsi,singbno_mn(:,ising),mfac,mpert,0)
         ENDIF
         singflx_mn(:,ising)=singflx_mn(:,ising)/area(ising)
c-----------------------------------------------------------------------
c     compute pseudo-chirikov parameter.
c-----------------------------------------------------------------------
         IF (ising==1) THEN 
            hdist=(singtype(ising+1)%psifac-respsi)/2.0
         ELSE IF (ising==msing) THEN
            hdist=(respsi-singtype(ising-1)%psifac)/2.0
         ELSE IF ((ising/=1).AND.(ising/=msing)) THEN
            hdist=MIN(singtype(ising+1)%psifac-respsi,
     $           respsi-singtype(ising-1)%psifac)/2.0
         ENDIF
         chirikov(ising)=island_hwidth(ising)/hdist
         IF(verbose) WRITE(*,'(1x,a6,es10.3,a6,f6.3,a25,es10.3)')
     $        "psi = ",singtype(ising)%psifac,
     $        ", q = ",singtype(ising)%q,
     $        ", total resonant field = ",
     $        ABS(singflx_mn(resnum(ising),ising))
      ENDDO

      CALL cspline_dealloc(fsp_sol)
      CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_singfld_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_SINGFLD: "//
     $        "Resonant fields, singular currents, and islands"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $        "jac_out = ",jac_out,"tmag_out =",tmag_out
         WRITE(out_unit,'(1x,a12,es17.8e3)')"sweet-spot =",spot
         WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a6,13(1x,a16))')"q","psi",
     $        "real(singflx)","imag(singflx)",
     $        "real(singcur)","imag(singcur)",
     $        "real(singbwp)","imag(singbwp)",
     $        "real(Delta)","imag(Delta)",
     $        "islandhwidth","chirikov"
         DO ising=1,msing
            WRITE(out_unit,'(1x,f6.3,13(es17.8e3))')
     $           singtype(ising)%q,singtype(ising)%psifac,
     $           REAL(singflx_mn(resnum(ising),ising)),
     $           AIMAG(singflx_mn(resnum(ising),ising)),
     $           REAL(singcur(ising)),AIMAG(singcur(ising)),
     $           REAL(singbwp(ising)),AIMAG(singbwp(ising)),
     $           REAL(delta(ising)),AIMAG(delta(ising)),
     $           island_hwidth(ising),chirikov(ising)
         ENDDO
         WRITE(out_unit,*)
      ENDIF

      IF(netcdf_flag .AND. msing>0)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n_rational",q_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "Phi_res", nf90_double,
     $      (/q_id,i_id/), p_id) )
         CALL check( nf90_put_att(fncid, p_id, "units", "T") )
         CALL check( nf90_put_att(fncid, p_id, "long_name",
     $     "Pitch resonant flux normalized by the surface area") )
         CALL check( nf90_def_var(fncid, "Delta", nf90_double,
     $      (/q_id,i_id/), d_id) )
         CALL check( nf90_put_att(fncid, d_id, "long_name",
     $     "Unitless Resonance Parameter $\partial_\psi \frac{"//
     $     "\delta B \cdot \nabla \psi}{B \cdot \nabla \theta}$") )
         CALL check( nf90_def_var(fncid, "I_res", nf90_double,
     $      (/q_id,i_id/), c_id) )
         CALL check( nf90_put_att(fncid, c_id, "units", "A") )
         CALL check( nf90_put_att(fncid, c_id, "long_name",
     $     "Pitch resonant current") )
         CALL check( nf90_def_var(fncid, "w_isl", nf90_double,
     $      (/q_id/), w_id) )
         CALL check( nf90_put_att(fncid, w_id, "units", "psi_n") )
         CALL check( nf90_put_att(fncid, w_id, "long_name",
     $     "Full width of saturated island") )
         CALL check( nf90_def_var(fncid, "K_isl", nf90_double,
     $      (/q_id/), k_id) )
         CALL check( nf90_put_att(fncid, k_id, "long_name",
     $     "Chirikov parameter of fully saturated islands") )
         CALL check( nf90_enddef(fncid) )
         singflx = (/(singflx_mn(resnum(ising),ising), ising=1,msing)/)
         CALL check( nf90_put_var(fncid, p_id,
     $      RESHAPE((/REAL(singflx), AIMAG(singflx)/), (/msing,2/))) )
         CALL check( nf90_put_var(fncid, d_id,
     $      RESHAPE((/REAL(delta), AIMAG(delta)/), (/msing,2/))) )
         CALL check( nf90_put_var(fncid, c_id,
     $      RESHAPE((/REAL(singcur), AIMAG(singcur)/), (/msing,2/))) )
         CALL check( nf90_put_var(fncid, w_id, 2*island_hwidth) )
         CALL check( nf90_put_var(fncid, k_id, chirikov) )
         CALL check( nf90_close(fncid) )
      ENDIF

      ! log singular response with harvest
      DO ising=1,msing
         aq(ising) = singtype(ising)%q
         asingflx(ising) = ABS(singflx_mn(resnum(ising),ising))
      ENDDO
      ierr=set_harvest_payload_dbl_array(hlog,"q"//nul,aq,msing)
      ierr=set_harvest_payload_dbl_array(hlog,"singcur"//nul,
     $     ABS(singcur),msing)
      ierr=set_harvest_payload_dbl_array(hlog,"singflx"//nul,
     $     asingflx,msing)

      IF (singcoup_set .AND. ALLOCATED(sbno_fun)) THEN
         sbnosurf=SQRT(ABS(DOT_PRODUCT(sbno_fun(1:mthsurf),
     $        sbno_fun(1:mthsurf)))/mthsurf)
         sbno_mn = MATMUL(fldflxmat,sbno_mn)
         DO icoup=1,nsingcoup
            DO ising=1,msing
               olap(icoup, ising) = DOT_PRODUCT(
     $            singcoup_out_vecs(icoup,:,ising),
     $            sbno_mn(:))
            ENDDO
         ENDDO
         op = ABS(olap) / sbnosurf * 1e2

         IF(ascii_flag)THEN
            WRITE(out_unit,*)"Overlap fields, overlap singular "//
     $           "currents, and overlap islands"
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a6,15(1x,a16))') "mode",
     $           "real(ovf)","imag(ovf)","overlap(%)",  ! Resonant flux RSV overlap
     $           "real(ovs)","imag(ovs)","overlap(%)",  ! Resonant current RSV overlap
     $           "real(ovi)","imag(ovi)","overlap(%)",  ! Island width RSV overlap
     $           "real(ovp)","imag(ovp)","overlap(%)",  ! Penetrated flux RSV overlap
     $           "real(ovd)","imag(ovd)","overlap(%)"   ! Unitless Delta RSV overlap
            DO ising=1,msing
               WRITE(out_unit,'(1x,I6,15(es17.8e3))') ising,
     $             REAL(olap(1,ising)),AIMAG(olap(1,ising)),op(1,ising),
     $             REAL(olap(2,ising)),AIMAG(olap(2,ising)),op(2,ising),
     $             REAL(olap(3,ising)),AIMAG(olap(3,ising)),op(3,ising),
     $             REAL(olap(4,ising)),AIMAG(olap(4,ising)),op(4,ising),
     $             REAL(olap(5,ising)),AIMAG(olap(5,ising)),op(5,ising)
            ENDDO
            WRITE(out_unit,*)
         ENDIF

         IF(netcdf_flag .AND. msing>0)THEN
            CALL check( nf90_open(mncfile,nf90_write,mncid) )
            CALL check( nf90_inq_dimid(mncid,"i",i_id) )
            CALL check( nf90_inq_dimid(mncid,"mode_C",m_id) )
            CALL check( nf90_redef(mncid))
            CALL check( nf90_def_var(mncid, "Phi_overlap", nf90_double,
     $         (/m_id,i_id/), p_id) )
            CALL check( nf90_put_att(mncid, p_id, "units", "unitless") )
            CALL check( nf90_put_att(mncid, p_id, "long_name",
     $        "Pitch resonant flux overlap") )
            CALL check( nf90_def_var(mncid, "Phi_overlap_norm",
     $         nf90_double,(/m_id/), pp_id) )
            CALL check( nf90_put_att(mncid, pp_id, "units", "unitless"))
            CALL check( nf90_put_att(mncid, pp_id, "long_name",
     $        "Pitch resonant flux overlap percentage") )
            CALL check( nf90_def_var(mncid, "I_overlap", nf90_double,
     $         (/m_id,i_id/), c_id) )
            CALL check( nf90_put_att(mncid, c_id, "units", "unitless") )
            CALL check( nf90_put_att(mncid, c_id, "long_name",
     $        "Pitch resonant current overlap") )
            CALL check( nf90_def_var(mncid, "I_overlap_norm",
     $         nf90_double,(/m_id,i_id/), cp_id) )
            CALL check( nf90_put_att(mncid, cp_id, "units", "unitless"))
            CALL check( nf90_put_att(mncid, cp_id, "long_name",
     $        "Pitch resonant current overlap percentage") )
            CALL check( nf90_def_var(mncid, "w_overlap", nf90_double,
     $         (/m_id,i_id/), w_id) )
            CALL check( nf90_put_att(mncid, w_id, "units", "untiless") )
            CALL check( nf90_put_att(mncid, w_id, "long_name",
     $        "Full width of saturated island overlap") )
            CALL check( nf90_def_var(mncid, "w_overlap_norm",
     $         nf90_double,(/m_id/), wp_id) )
            CALL check( nf90_put_att(mncid, wp_id, "units", "untiless"))
            CALL check( nf90_put_att(mncid, wp_id, "long_name",
     $        "Full width of saturated island overlap percentage") )
            CALL check( nf90_def_var(mncid, "Pen_overlap", nf90_double,
     $         (/m_id, i_id/), n_id) )
            CALL check( nf90_put_att(mncid, n_id, "units", "untiless") )
            CALL check( nf90_put_att(mncid, n_id, "long_name",
     $        "Penetrated flux overlap") )
            CALL check( nf90_def_var(mncid, "Pen_overlap_norm",
     $         nf90_double,(/m_id/), np_id) )
            CALL check( nf90_put_att(mncid, np_id, "units", "untiless"))
            CALL check( nf90_put_att(mncid, np_id, "long_name",
     $        "Penetrated flux overlap percentage") )
            CALL check( nf90_def_var(mncid, "Delta_overlap",nf90_double,
     $         (/m_id, i_id/), d_id) )
            CALL check( nf90_put_att(mncid, d_id, "units", "untiless") )
            CALL check( nf90_put_att(mncid, d_id, "long_name",
     $        "Extrenal Delta prime overlap") )
            CALL check( nf90_def_var(mncid, "Delta_overlap_norm",
     $         nf90_double,(/m_id/), dp_id) )
            CALL check( nf90_put_att(mncid, dp_id, "units", "untiless"))
            CALL check( nf90_put_att(mncid, dp_id, "long_name",
     $        "Extrenal Delta prime overlap percentage") )
            CALL check( nf90_enddef(mncid) )
            CALL check( nf90_put_var(mncid, p_id, RESHAPE((/
     $         REAL(olap(1,:)),AIMAG(olap(1,:))/), (/msing,2/))) )
            CALL check( nf90_put_var(mncid, c_id, RESHAPE((/
     $         REAL(olap(2,:)),AIMAG(olap(2,:))/), (/msing,2/))) )
            CALL check( nf90_put_var(mncid, w_id, RESHAPE((/
     $         REAL(olap(3,:)),AIMAG(olap(3,:))/), (/msing,2/))) )
            CALL check( nf90_put_var(mncid, n_id, RESHAPE((/
     $         REAL(olap(4,:)),AIMAG(olap(4,:))/), (/msing,2/))) )
            CALL check( nf90_put_var(mncid, d_id, RESHAPE((/
     $         REAL(olap(5,:)),AIMAG(olap(5,:))/), (/msing,2/))) )
            CALL check( nf90_put_var(mncid, pp_id, op(1,:)) )
            CALL check( nf90_put_var(mncid, cp_id, op(2,:)) )
            CALL check( nf90_put_var(mncid, wp_id, op(3,:)) )
            CALL check( nf90_put_var(mncid, np_id, op(4,:)) )
            CALL check( nf90_put_var(mncid, dp_id, op(5,:)) )
            CALL check( nf90_close(mncid) )
         ENDIF

         ! information already in netcdf from control_filter subroutine
         ! log svd overlap with harvest
         ierr=set_harvest_payload_dbl_array(hlog,'overlap_percent'//nul,
     $        op(1,:),msing)

         IF(osing<msing)THEN
            DO icoup=1,nsingcoup
               DO ising=1,osing
                  olap(icoup, ising) = DOT_PRODUCT(
     $               localcoup_out_vecs(icoup,:,ising),
     $               sbno_mn(:))
               ENDDO
            ENDDO
            op = ABS(olap) / sbnosurf * 1e2

            IF(ascii_flag)THEN
               WRITE(out_unit,*)"Local overlap fields, overlap "//
     $              "singular currents, and overlap islands"
               WRITE(out_unit,*)
               WRITE(out_unit,'(1x,a6,15(1x,a16))') "mode",
     $           "real(ovf)","imag(ovf)","overlap(%)",  ! Resonant flux RSV overlap
     $           "real(ovs)","imag(ovs)","overlap(%)",  ! Resonant current RSV overlap
     $           "real(ovi)","imag(ovi)","overlap(%)",  ! Island width RSV overlap
     $           "real(ovp)","imag(ovp)","overlap(%)",  ! Penetrated flux RSV overlap
     $           "real(ovd)","imag(ovd)","overlap(%)"   ! Unitless Delta RSV overlap
               DO ising=1,osing
                  WRITE(out_unit,'(1x,I6,15(es17.8e3))') ising,
     $             REAL(olap(1,ising)),AIMAG(olap(1,ising)),op(1,ising),
     $             REAL(olap(2,ising)),AIMAG(olap(2,ising)),op(2,ising),
     $             REAL(olap(3,ising)),AIMAG(olap(3,ising)),op(3,ising),
     $             REAL(olap(4,ising)),AIMAG(olap(4,ising)),op(4,ising),
     $             REAL(olap(5,ising)),AIMAG(olap(5,ising)),op(5,ising)
               ENDDO
               WRITE(out_unit,*)
            ENDIF

            IF(netcdf_flag .AND. osing>0)THEN
               CALL check( nf90_open(mncfile,nf90_write,mncid) )
               CALL check( nf90_inq_dimid(mncid,"i",i_id) )
               CALL check( nf90_inq_dimid(mncid,"mode_C_local",m_id) )
               CALL check( nf90_redef(mncid))
               CALL check( nf90_def_var(mncid, "Phi_local_overlap",
     $            nf90_double,(/m_id,i_id/), p_id) )
               CALL check( nf90_put_att(mncid, p_id,"units","unitless"))
               CALL check( nf90_put_att(mncid, p_id, "long_name",
     $           "Local pitch resonant flux overlap") )
               CALL check( nf90_def_var(mncid, "Phi_local_overlap_norm",
     $            nf90_double,(/m_id/), pp_id) )
               CALL check( nf90_put_att(mncid,pp_id,"units","unitless"))
               CALL check( nf90_put_att(mncid, pp_id, "long_name",
     $           "Local pitch resonant flux overlap percentage") )
               CALL check( nf90_def_var(mncid, "I_local_overlap",
     $            nf90_double,(/m_id,i_id/), c_id) )
               CALL check( nf90_put_att(mncid, c_id,"units","unitless"))
               CALL check( nf90_put_att(mncid, c_id, "long_name",
     $           "Local pitch resonant current overlap") )
               CALL check( nf90_def_var(mncid, "I_local_overlap_norm",
     $            nf90_double,(/m_id/), cp_id) )
               CALL check( nf90_put_att(mncid,cp_id,"units","unitless"))
               CALL check( nf90_put_att(mncid, cp_id, "long_name",
     $           "Local pitch resonant current overlap percentage") )
               CALL check( nf90_def_var(mncid, "w_local_overlap",
     $            nf90_double,(/m_id,i_id/), w_id) )
               CALL check( nf90_put_att(mncid, w_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, w_id, "long_name",
     $           "Local full width of saturated island overlap") )
               CALL check( nf90_def_var(mncid, "w_local_overlap_norm",
     $            nf90_double,(/m_id/), wp_id) )
               CALL check( nf90_put_att(mncid,wp_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, wp_id, "long_name",
     $      "Local full width of saturated island overlap percentage") )
               CALL check( nf90_def_var(mncid, "Pen_local_overlap",
     $            nf90_double,(/m_id, i_id/), n_id) )
               CALL check( nf90_put_att(mncid, n_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, n_id, "long_name",
     $           "Local penetrated flux overlap") )
               CALL check( nf90_def_var(mncid, "Pen_local_overlap_norm",
     $            nf90_double,(/m_id/), np_id) )
               CALL check( nf90_put_att(mncid,np_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, np_id, "long_name",
     $           "Local penetrated flux overlap percentage") )
               CALL check( nf90_def_var(mncid, "Delta_local_overlap",
     $            nf90_double,(/m_id, i_id/), d_id) )
               CALL check( nf90_put_att(mncid, d_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, d_id, "long_name",
     $           "Local extrenal Delta prime overlap") )
              CALL check( nf90_def_var(mncid,"Delta_local_overlap_norm",
     $            nf90_double,(/m_id/), dp_id) )
               CALL check( nf90_put_att(mncid,dp_id,"units","untiless"))
               CALL check( nf90_put_att(mncid, dp_id, "long_name",
     $           "Local extrenal Delta prime overlap percentage") )
               CALL check( nf90_enddef(mncid) )
               CALL check( nf90_put_var(mncid, p_id, RESHAPE((/
     $            REAL(olap(1,:osing)),AIMAG(olap(1,:osing))/),
     $            (/osing,2/))) )
               CALL check( nf90_put_var(mncid, c_id, RESHAPE((/
     $            REAL(olap(2,:osing)),AIMAG(olap(2,:osing))/),
     $            (/osing,2/))) )
               CALL check( nf90_put_var(mncid, w_id, RESHAPE((/
     $            REAL(olap(3,:osing)),AIMAG(olap(3,:osing))/),
     $            (/osing,2/))) )
               CALL check( nf90_put_var(mncid, n_id, RESHAPE((/
     $            REAL(olap(4,:osing)),AIMAG(olap(4,:osing))/),
     $            (/osing,2/))) )
               CALL check( nf90_put_var(mncid, d_id, RESHAPE((/
     $            REAL(olap(5,:osing)),AIMAG(olap(5,:osing))/),
     $            (/osing,2/))) )
               CALL check( nf90_put_var(mncid, pp_id, op(1,1:osing)) )
               CALL check( nf90_put_var(mncid, cp_id, op(2,1:osing)) )
               CALL check( nf90_put_var(mncid, wp_id, op(3,1:osing)) )
               CALL check( nf90_put_var(mncid, np_id, op(4,1:osing)) )
               CALL check( nf90_put_var(mncid, dp_id, op(5,1:osing)) )
               CALL check( nf90_close(mncid) )
            ENDIF

         ENDIF
      ENDIF
        
      IF(ascii_flag) CALL ascii_close(out_unit)
      IF(singcoup_set) THEN
         DEALLOCATE(singcoup_out_vecs)
         IF(osing<msing) THEN
            DEALLOCATE(localcoup_out_vecs)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_singfld
c-----------------------------------------------------------------------
c     subprogram 5. gpout_vsingfld.
c     compute resonant field by coils.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_vsingfld()
c-----------------------------------------------------------------------
      INTEGER :: ising,i
      REAL(r8) :: hdist,shear,area
      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: visland_hwidth,vchirikov
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: vcmn

      COMPLEX(r8), DIMENSION(msing) :: vflxmn
      REAL(r8), DIMENSION(0:mthsurf) :: unitfun

      INTEGER :: i_id,q_id,p_id,w_id,k_id
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing resonant field from coils"
      ALLOCATE(vcmn(cmpert))
      vcmn=0
      DO ising=1,msing
         CALL field_bs_psi(singtype(ising)%psifac,vcmn,2)
         resnum(ising)=NINT(singtype(ising)%q*nn)-mlow+1
         DO i=1,cmpert
            IF (cmlow-mlow+i==resnum(ising)) THEN
               vflxmn(ising)=vcmn(i)
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     compute half-width of magnetic island.
c-----------------------------------------------------------------------
         shear=mfac(resnum(ising))*
     $        singtype(ising)%q1/singtype(ising)%q**2
         unitfun=1.0
         area=issurfint(unitfun,mthsurf,singtype(ising)%psifac,0,0)
         visland_hwidth(ising)=
     $        SQRT(ABS(4*vflxmn(ising)*area/
     $        (twopi*shear*singtype(ising)%q*chi1)))
         IF (ising==1) THEN 
            hdist=(singtype(ising+1)%psifac-singtype(ising)%psifac)/2.0
         ELSE IF (ising==msing) THEN
            hdist=(singtype(ising)%psifac-singtype(ising-1)%psifac)/2.0
         ELSE IF ((ising/=1).AND.(ising/=msing)) THEN
            hdist=MIN(singtype(ising+1)%psifac-singtype(ising)%psifac,
     $           singtype(ising)%psifac-singtype(ising-1)%psifac)/2.0
         ENDIF
         vchirikov(ising)=visland_hwidth(ising)/hdist
         IF(verbose) WRITE(*,'(1x,a6,es10.3,a6,f6.3,a25,es10.3)')
     $        "psi = ",singtype(ising)%psifac,
     $        ", q = ",singtype(ising)%q,
     $        ", vacuum resonant field = ",
     $        ABS(vflxmn(ising))         
         
      ENDDO
      DEALLOCATE(vcmn)
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_vsingfld_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_SINGFLD: "//
     $        "Resonant fields and islands from coils"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $        "jac_out = ",jac_out,"tmag_out =",tmag_out
         WRITE(out_unit,'(1x,a14)')"sweet-spot = 0"
         WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a6,5(1x,a16))')"q","psi",
     $        "real(singflx)","imag(singflx)",
     $        "islandhwidth","chirikov"
         DO ising=1,msing
            WRITE(out_unit,'(1x,f6.3,5(es17.8e3))')
     $           singtype(ising)%q,singtype(ising)%psifac,
     $           REAL(vflxmn(ising)),AIMAG(vflxmn(ising)),
     $           visland_hwidth(ising),vchirikov(ising)
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag .AND. msing>0)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n_rational",q_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "Phi_res_v", nf90_double,
     $      (/q_id,i_id/), p_id) )
         CALL check( nf90_put_att(fncid, p_id, "units", "Wb") )
         CALL check( nf90_put_att(fncid, p_id, "long_name",
     $     "Pitch resonant vacuum flux") )
         CALL check( nf90_def_var(fncid, "w_isl_v", nf90_double,
     $      (/q_id/), w_id) )
         CALL check( nf90_put_att(fncid, w_id, "units", "psi_n") )
         CALL check( nf90_put_att(fncid, w_id, "long_name",
     $     "Full width of vacuum island") )
         CALL check( nf90_def_var(fncid, "K_isl_v", nf90_double,
     $      (/q_id/), k_id) )
         CALL check( nf90_put_att(fncid, k_id, "long_name",
     $     "Chirikov parameter of vacuum islands") )
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid, p_id,
     $      RESHAPE((/REAL(vflxmn), AIMAG(vflxmn)/), (/msing,2/))) )
         CALL check( nf90_put_var(fncid, w_id, 2*visland_hwidth) )
         CALL check( nf90_put_var(fncid, k_id, vchirikov) )
         CALL check( nf90_close(fncid) )
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_vsingfld
c-----------------------------------------------------------------------
c     subprogram 6. gpout_dw.
c     restore energy and torque profiles from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_dw(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep
      INTEGER :: i_id, p_id, t_id, d_id
      TYPE(cspline_type) :: dwk
c-----------------------------------------------------------------------
c     build and spline energy and torque profiles.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xspmn)
      CALL cspline_alloc(dwk,mstep,1)
      dwk%xs=psifac
      ! division by 2 corrects quadratic use of 2A_+n instead of proper A_+n + A_-n
      dwk%fs(:,1) = SUM(CONJG(u1%fs) * u2%fs, DIM=2) / 2
     $   / (2 * mu0)        ! convert to Jules
     $   * 2 * nn * ifac    ! convert from energy to torque
      CALL cspline_fit(dwk,"extrap")

      WRITE(*,*)"Restoring energy and torque profiles"

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_dw_profile_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_dw: "//
     $        "Energy and torque profiles by self-consistent solutions."
         WRITE(out_unit,*)
         WRITE(out_unit,'(6(1x,a16))')"psi","T_phi","2ndeltaW",
     $        "int(T_phi)","int(2ndeltaW)","dv/dpsi_n"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            CALL cspline_eval(dwk,psifac(istep),1)
            CALL spline_eval(sq,psifac(istep),0)
            WRITE(out_unit,'(6(1x,es16.8))')
     $           psifac(istep),REAL(dwk%f1(1)),AIMAG(dwk%f1(1)),
     $           REAL(dwk%fs(istep,1)),AIMAG(dwk%fs(istep,1)),sq%f(3)
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "T", nf90_double,
     $               (/p_id,i_id/), t_id) )
         CALL check( nf90_put_att(fncid, t_id, "long_name",
     $               "Integrated toroidal torque") )
         CALL check( nf90_put_att(fncid, t_id, "units", "Nm") )
         CALL check( nf90_def_var(fncid, "dTdpsi", nf90_double,
     $               (/p_id,i_id/), d_id) )
         CALL check( nf90_put_att(fncid, d_id, "long_name",
     $               "Toroidal torque profile") )
         CALL check( nf90_put_att(fncid, d_id, "units", "Nm") )
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid, t_id, RESHAPE(
     $      (/REAL(dwk%fs(:,1)), AIMAG(dwk%fs(:,1))/), (/mstep,2/))) )
         CALL check( nf90_put_var(fncid, d_id,RESHAPE(
     $      (/REAL(dwk%fs1(:,1)), AIMAG(dwk%fs1(:,1))/), (/mstep,2/))) )
         CALL check( nf90_close(fncid) )
      ENDIF


      CALL cspline_dealloc(dwk)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpout_dw
c-----------------------------------------------------------------------
c     subprogram 7. gpout_dw_matrix.
c     restore energy and torque response matrix from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_dw_matrix(coil_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      LOGICAL, INTENT(IN) :: coil_flag

      INTEGER :: ipert,istep,i,j,iindex
      REAL(r8) :: jarea,ileft
      COMPLEX(r8) :: t1,t2

      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8), DIMENSION(0:mthsurf) :: units
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,sqrtamat,ptof
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: bsurfmat,dwks,dwk,
     $     gind,gindp,gres,gresp

      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: coilmn
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: tmat,mmat,mdagger
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: gcoil

      INTEGER :: p_id,j_id,m_id,c_id,k_id,i_id,ci_id,cp_id,m1_id,m2_id,
     $    te_id,tc_id,sl_id,cn_id

c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(bsurfmat(mstep,mpert,mpert),dwks(mstep,mpert,mpert),
     $   dwk(mstep,mpert,mpert),gind(mstep,mpert,mpert),
     $   gindp(mstep,mpert,mpert),gres(mstep,mpert,mpert),
     $   gresp(mstep,mpert,mpert))
c-----------------------------------------------------------------------
c     calculate the scalar surface area for normalizations.
c-----------------------------------------------------------------------
      units = 1
      jarea = issurfint(units,mthsurf,psilim,0,0)
c-----------------------------------------------------------------------
c     calculate dimensionless half area weighting matrix W.
c-----------------------------------------------------------------------
!     W_m,m' = int{sqrt(J|delpsi|)exp[-i*(m-m')t]dt}/int{sqrt(J|delpsi|)dt}
      sqrtamat = 0
      DO i=1,mpert
         sqrtamat(i,i) = 1.0
         CALL gpeq_weight(psilim,sqrtamat(:,i),mfac,mpert,2) ! A^1/2
      ENDDO
      ptof = sqrtamat * sqrt(jarea)
c-----------------------------------------------------------------------
c     call solutions.
c-----------------------------------------------------------------------
      WRITE(*,*) "Call solutions for general response matrix functions"
      edge_flag=.TRUE.
      DO ipert=1,mpert
         edge_mn=0
         edge_mn(ipert)=1.0
         CALL idcon_build(0,edge_mn)
         DO istep=1,mstep
            bsurfmat(istep,:,ipert)=u1%fs(istep,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     construct dws response matrix (normalized by xi).
c-----------------------------------------------------------------------
      WRITE(*,*)"Build general response matrix functions"
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a30)')
     $        "volume = ",iindex,"% response matrix calculations"
         temp1=CONJG(TRANSPOSE(soltype(istep)%u(:,:,1)))
         temp2=CONJG(TRANSPOSE(soltype(istep)%u(:,:,2)))
         CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         dwks(istep,:,:)=CONJG(TRANSPOSE(temp2))/(2.0*mu0)*2*nn*ifac
c-----------------------------------------------------------------------
c     construct dwk response matrix (normalized by xi boundary).
c-----------------------------------------------------------------------
         temp1=MATMUL(dwks(istep,:,:),bsurfmat(istep,:,:))
         temp2=CONJG(TRANSPOSE(bsurfmat(istep,:,:)))
         dwk(istep,:,:)=MATMUL(temp2,temp1)
c-----------------------------------------------------------------------
c     make general (2x) inductance matrix (normalized by phi boundary).
c-----------------------------------------------------------------------
         DO i=1,mpert
            DO j=1,mpert
               t1=-1/(chi1*(mfac(i)-nn*qlim)*twopi*ifac)
               t2=1/(chi1*(mfac(j)-nn*qlim)*twopi*ifac)
               gind(istep,i,j)=t1*dwk(istep,i,j)*t2
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     make general response matrix (normalized by phix boundary).
c-----------------------------------------------------------------------
         temp1=MATMUL(gind(istep,:,:),permeabmats(resp_index,:,:))
         temp2=CONJG(TRANSPOSE(permeabmats(resp_index,:,:)))
         gres(istep,:,:)=MATMUL(temp2,temp1)
c-----------------------------------------------------------------------
c     make coordinate-independent matrix (normalized by phi power).
c-----------------------------------------------------------------------
         temp1=MATMUL(gind(istep,:,:),ptof)
         temp2=CONJG(TRANSPOSE(ptof))
         gindp(istep,:,:)=MATMUL(temp2,temp1)
c-----------------------------------------------------------------------
c     make coordinate-independent matrix (normalized by phix power).
c-----------------------------------------------------------------------
         temp1=MATMUL(gres(istep,:,:),ptof)
         temp2=CONJG(TRANSPOSE(ptof))
         gresp(istep,:,:)=MATMUL(temp2,temp1)
      ENDDO
c-----------------------------------------------------------------------
c     write general response matrix functions.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_dw_matrix_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_dw_matrix: "//
     $        "Self-consistent response matrix functions."
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a16,2(1x,a4),10(1x,a16))')"psi","i","j",
     $        "real(T_x)","imag(T_x)","real(T_f)","imag(T_f)",
     $        "real(T_ef)","imag(T_ef)","real(T_fp)","imag(T_fp)",
     $        "real(T_efp)","imag(T_efp)"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            DO i=1,mpert
               DO j=1,mpert
                  WRITE(out_unit,'(1x,es16.8,2(1x,I4),10(1x,es16.8))')
     $                 psifac(istep),i,j,
     $                 REAL(dwk(istep,i,j)),AIMAG(dwk(istep,i,j)),
     $                 REAL(gind(istep,i,j)),AIMAG(gind(istep,i,j)),
     $                 REAL(gres(istep,i,j)),AIMAG(gres(istep,i,j)),
     $                 REAL(gindp(istep,i,j)),AIMAG(gindp(istep,i,j)),
     $                 REAL(gresp(istep,i,j)),AIMAG(gresp(istep,i,j))
               ENDDO
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_dim(fncid, "m", mpert, m_id) )
         CALL check( nf90_def_var(fncid, "m", nf90_int,m_id,m1_id))
         CALL check( nf90_def_dim(fncid, "m_prime", mpert, j_id) )
         CALL check( nf90_def_var(fncid, "m_prime",nf90_int,j_id,m2_id))
         CALL check( nf90_def_var(fncid, "T_xe", nf90_double,
     $               (/p_id,m_id,j_id,i_id/), te_id) )
         CALL check( nf90_put_att(fncid, te_id, "long_name",
     $               "Energy-norm external flux torque matrix") )
         CALL check( nf90_put_att(fncid, te_id, "units", "Nm/T^2") )
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid, m1_id, mfac) )
         CALL check( nf90_put_var(fncid, m2_id, mfac) )
         CALL check( nf90_put_var(fncid, te_id, RESHAPE((/REAL(gresp),
     $               AIMAG(gresp)/), (/mstep,mpert,mpert,2/))) )
         CALL check( nf90_close(fncid) )
      ENDIF

c-----------------------------------------------------------------------
c     construct coil-torque response matrix.
c----------------------------------------------------------------------
      IF (coil_flag) THEN
         ! form mutual inductance between coils and plasma surface
         ALLOCATE(mmat(mpert,coil_num),mdagger(coil_num,mpert))
         ALLOCATE(coilmn(cmpert))
         DO j=1,coil_num
            CALL field_bs_psi(psilim,coilmn,1,op_start=j,op_stop=j)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  mmat(cmlow-mlow+i,j)=coilmn(i)
               ENDIF
            ENDDO
         ENDDO
         mdagger = CONJG(TRANSPOSE(mmat))
         DEALLOCATE(coilmn)

         WRITE(*,*)"Build coil response matrix functions"
         ALLOCATE(gcoil(mstep,coil_num,coil_num),tmat(mpert,mpert))
         DO istep=1,mstep
            tmat = gres(istep,:,:)
            gcoil(istep,:,:)=MATMUL(mdagger,MATMUL(tmat,mmat))
         ENDDO
c-----------------------------------------------------------------------
c     write coil response matrix functions.
c-----------------------------------------------------------------------
         IF(ascii_flag)THEN
            CALL ascii_open(out_unit,"gpec_dw_coil_matrix_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_dw_coil_matrix: "//
     $           "Self-consistent coil response matrix functions."
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a16,2(1x,a4),2(1x,a16))')"psi","m",
     $           "m_prime","real(T_coil)","imag(T_coil)"
            DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
               DO i=1,coil_num
                  DO j=1,coil_num
                     WRITE(out_unit,'(es17.8,2(1x,I4),2(es17.8))')
     $                    psifac(istep),i,j,
     $                    REAL(gcoil(istep,i,j)),AIMAG(gcoil(istep,i,j))
                  ENDDO
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         IF(netcdf_flag)THEN
            CALL check( nf90_open(fncfile,nf90_write,fncid) )
            CALL check( nf90_inq_dimid(fncid,"i",i_id) )
            CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
            CALL check( nf90_redef(fncid))
            CALL check( nf90_def_dim(fncid,"coil_index",coil_num,c_id))
            CALL check( nf90_def_var(fncid, "coil_index", nf90_int,
     $                  (/c_id/), ci_id) )
            CALL check( nf90_def_dim(fncid, "coil_index_prime",
     $                  coil_num, k_id) )
            CALL check( nf90_def_var(fncid, "coil_index_prime",
     $                  nf90_int, (/k_id/), cp_id) )
            CALL check( nf90_def_dim(fncid,"coil_strlen",24,sl_id))
            CALL check( nf90_def_var(fncid, "coil_name", nf90_char,
     $                  (/sl_id, c_id/), cn_id) )
            CALL check( nf90_def_var(fncid, "T_coil", nf90_double,
     $                  (/p_id, c_id, k_id, i_id/), tc_id) )
            CALL check( nf90_put_att(fncid,tc_id,"long_name",
     $                  "Coil-space torque matrix") )
            CALL check( nf90_put_att(fncid, tc_id, "units", "~Nm/A^2") )
            CALL check( nf90_enddef(fncid) )
            CALL check( nf90_put_var(fncid, ci_id,(/(i,i=1,coil_num)/)))
            CALL check( nf90_put_var(fncid, cp_id,(/(i,i=1,coil_num)/)))
            CALL check( nf90_put_var(fncid,cn_id,coil_name(1:coil_num)))
            CALL check( nf90_put_var(fncid, tc_id,RESHAPE((/REAL(gcoil),
     $                  AIMAG(gcoil)/), (/mstep,coil_num,coil_num,2/))))
            CALL check( nf90_close(fncid) )
         ENDIF
         DEALLOCATE(gcoil,tmat,mmat,mdagger)
      ENDIF
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(bsurfmat,dwks,dwk,gind,gindp,gres,gresp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpout_dw_matrix
c-----------------------------------------------------------------------
c     subprogram 8. gpout_pmodb.
c     compute perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_pmodb(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,itheta,iindex,cstep,tout
      REAL(r8) :: ileft,jac,psi

      INTEGER :: p_id,t_id,i_id,m_id,r_id,z_id,b_id,bme_id,be_id,
     $   bml_id,bl_id,xm_id,x_id,km_id,k_id,rzstat

      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs,zs,equilbfun
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,
     $     divxprp_mn,curv_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmz_fun,xvt_fun,xvz_fun,xmt_fun,xsp1_fun
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: eulbparmns,
     $     lagbparmns,divxprpmns,curvmns,
     $     xmp1mns,xspmns,xmsmns,xmtmns,xmzmns
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: eulbparmout,
     $     lagbparmout,divxprpmout,curvmout
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: eulbparfun,lagbparfun,
     $     divxprpfun,curvfun,eulbparfout,lagbparfout,divxprpfout,
     $     curvfout

      REAL(r8), DIMENSION(:), ALLOCATABLE :: psis,ches,chex
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::chea,chelagbmns

      TYPE(cspline_type) :: cspl,chespl
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf),
     $   equilbfun(mstep,0:mthsurf))
      ALLOCATE(eulbparmns(mstep,mpert),lagbparmns(mstep,mpert),
     $   divxprpmns(mstep,mpert),curvmns(mstep,mpert),
     $   xmp1mns(mstep,mpert),xspmns(mstep,mpert),xmsmns(mstep,mpert),
     $   xmtmns(mstep,mpert),xmzmns(mstep,mpert))
      ALLOCATE(eulbparmout(mstep,lmpert),lagbparmout(mstep,lmpert),
     $   divxprpmout(mstep,lmpert),curvmout(mstep,lmpert))
      ALLOCATE(eulbparfun(mstep,0:mthsurf),lagbparfun(mstep,0:mthsurf),
     $   divxprpfun(mstep,0:mthsurf),curvfun(mstep,0:mthsurf),
     $   eulbparfout(mstep,0:mthsurf),lagbparfout(mstep,0:mthsurf),
     $   divxprpfout(mstep,0:mthsurf),curvfout(mstep,0:mthsurf))
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing |b| and del.x_prp"
      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta

      CALL idcon_build(egnum,xspmn)
      CALL gpeq_alloc
      tout = tmag_out

      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a32)')
     $        "volume = ",iindex,"% |b| and del.x_prp computations"
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL gpeq_sol(psifac(istep))
         CALL gpeq_contra(psifac(istep))
         CALL gpeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b variations in hamada.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         CALL iscdftb(mfac,mpert,xmt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xmz_fun,mthsurf,xmz_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)

         CALL spline_eval(sq,psifac(istep),1)

         singfac=mfac-nn*sq%f(4)
         xsp1_mn=xsp1_mn*(singfac**2/(singfac**2+reg_spot**2))
         CALL iscdftb(mfac,mpert,xsp1_fun,mthsurf,xsp1_mn)

         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),0)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),0)
            jac=rzphi%f(4)

            cspl%fs(itheta,1)=xmt_fun(itheta)-
     $           (chi1/eqfun%f(1))**2/jac*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))
            cspl%fs(itheta,2)=xmz_fun(itheta)-
     $           sq%f(4)*(chi1/eqfun%f(1))**2/jac*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))
         ENDDO

         CALL cspline_fit(cspl,"periodic")

         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            CALL cspline_eval(cspl,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jac1=rzphi%fx(4)

            eulbparfun(istep,itheta)=
     $           chi1*(bvt_fun(itheta)+sq%f(4)*bvz_fun(itheta))
     $           /(rzphi%f(4)*eqfun%f(1))
            lagbparfun(istep,itheta)=
     $           eulbparfun(istep,itheta)+
     $           xsp_fun(itheta)*eqfun%fx(1)+
     $           xms_fun(itheta)/(chi1*sq%f(4))*eqfun%fy(1)
            divxprpfun(istep,itheta)=xsp1_fun(itheta)+
     $           (jac1/jac)*xsp_fun(itheta)+
     $           cspl%f1(1)/jac-(twopi*ifac*nn)*cspl%f(2)/jac
            divxprpfun(istep,itheta)=-eqfun%f(1)*
     $           divxprpfun(istep,itheta)
            curvfun(istep,itheta)=-xsp_fun(itheta)/eqfun%f(1)*sq%f1(2)-
     $           (xsp_fun(itheta)*eqfun%fx(1)+eqfun%fy(1)*
     $           (xms_fun(itheta)/(chi1*sq%f(4))+
     $           xmz_fun(itheta)/(jac*sq%f(4))-
     $           (chi1/(jac*eqfun%f(1)))**2*
     $           (xvt_fun(itheta)+sq%f(4)*xvz_fun(itheta))))
            equilbfun(istep,itheta) = eqfun%f(1)
            dphi(itheta) = rzphi%f(3)
         ENDDO
         CALL iscdftf(mfac,mpert,eulbparfun(istep,:),
     $        mthsurf,eulbpar_mn)
         CALL iscdftf(mfac,mpert,lagbparfun(istep,:),
     $        mthsurf,lagbpar_mn)
         CALL iscdftf(mfac,mpert,divxprpfun(istep,:),
     $        mthsurf,divxprp_mn)
         CALL iscdftf(mfac,mpert,curvfun(istep,:),
     $        mthsurf,curv_mn)
c-----------------------------------------------------------------------
c     Save working coordinate verions for pent files.
c-----------------------------------------------------------------------
         lagbparmns(istep,:) = lagbpar_mn
         divxprpmns(istep,:) = divxprp_mn
         curvmns(istep,:) = curv_mn
         xmp1mns(istep,:) = xmp1_mn
         xspmns(istep,:) = xsp_mn
         xmsmns(istep,:) = xms_mn
         xmtmns(istep,:) = xmt_mn
         xmzmns(istep,:) = xmz_mn   
c-----------------------------------------------------------------------
c     decompose components on the given coordinates.
c-----------------------------------------------------------------------
         i = istep
         psi = psifac(istep)
         CALL gpeq_bcoordsout(eulbparmout(i,:),eulbpar_mn,psi,tout,0)
         CALL gpeq_bcoordsout(lagbparmout(i,:),lagbpar_mn,psi,tout,0)
         CALL gpeq_bcoordsout(divxprpmout(i,:),divxprp_mn,psi,tout,0)
         CALL gpeq_bcoordsout(curvmout(i,:)   ,curv_mn,psi,tout,0)
         eulbparfout(istep,:)=eulbparfun(istep,:)*EXP(ifac*nn*dphi)
         lagbparfout(istep,:)=lagbparfun(istep,:)*EXP(ifac*nn*dphi)
         divxprpfout(istep,:)=divxprpfun(istep,:)*EXP(ifac*nn*dphi)
         curvfout(istep,:)=curvfun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      CALL cspline_dealloc(cspl)

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"m_out",m_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_inq_dimid(fncid,"theta_dcon",t_id) )
         CALL check( nf90_redef(fncid))
         rzstat = nf90_inq_varid(fncid, "R", r_id) ! check if R,z already stored
         IF(rzstat/=nf90_noerr)THEN
            CALL check( nf90_def_var(fncid, "R", nf90_double,
     $                  (/p_id,t_id/),r_id) )
            CALL check( nf90_put_att(fncid,r_id,"long_name",
     $                  "Major radius") )
            CALL check( nf90_put_att(fncid,r_id,"units","m") )
            CALL check( nf90_def_var(fncid, "z", nf90_double,
     $                  (/p_id,t_id/),z_id) )
            CALL check( nf90_put_att(fncid,z_id,"long_name",
     $                  "Vertical position") )
            CALL check( nf90_put_att(fncid,z_id,"units","m") )
         ENDIF
         CALL check( nf90_def_var(fncid, "b_eul_fun", nf90_double,
     $               (/p_id,t_id,i_id/),be_id) )
         CALL check( nf90_put_att(fncid,be_id,"long_name",
     $               "Eulerian perturbed field") )
         CALL check( nf90_put_att(fncid,be_id,"units","Tesla") )
         CALL check( nf90_def_var(fncid, "b_lag_fun", nf90_double,
     $               (/p_id,t_id,i_id/),bl_id) )
         CALL check( nf90_put_att(fncid,bl_id,"long_name",
     $               "Lagrangian perturbed field") )
         CALL check( nf90_put_att(fncid,bl_id,"units","Tesla") )
         CALL check( nf90_def_var(fncid, "b_eul", nf90_double,
     $               (/p_id,m_id,i_id/),bme_id) )
         CALL check( nf90_put_att(fncid,bme_id,"long_name",
     $               "Eulerian perturbed field") )
         CALL check( nf90_put_att(fncid,bme_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,bme_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "b_lag", nf90_double,
     $               (/p_id,m_id,i_id/),bml_id) )
         CALL check( nf90_put_att(fncid,bml_id,"long_name",
     $               "Lagrangian perturbed field") )
         CALL check( nf90_put_att(fncid,bml_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,bml_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "Bdivxi_perp_fun", nf90_double,
     $               (/p_id,t_id,i_id/),x_id) )
         CALL check( nf90_put_att(fncid,x_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_def_var(fncid, "Bdivxi_perp", nf90_double,
     $               (/p_id,m_id,i_id/),xm_id) )
         CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_put_att(fncid,xm_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "Bkappaxi_perp_fun",
     $               nf90_double, (/p_id,t_id,i_id/),k_id) )
         CALL check( nf90_put_att(fncid,k_id,"long_name",
     $               "Curvature component of the normal displacement") )
         CALL check( nf90_def_var(fncid, "Bkappaxi_perp", nf90_double,
     $               (/p_id,m_id,i_id/),km_id) )
         CALL check( nf90_put_att(fncid,km_id,"long_name",
     $               "Curvature component of the normal displacement") )
         CALL check( nf90_put_att(fncid,km_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "B", nf90_double,
     $               (/p_id,t_id/),b_id) )
         CALL check( nf90_put_att(fncid,b_id,"long_name",
     $               "Equilibrium field strength") )
         CALL check( nf90_enddef(fncid) )
         IF(rzstat/=nf90_noerr)THEN
            CALL check( nf90_put_var(fncid,r_id,rs) )
            CALL check( nf90_put_var(fncid,z_id,zs) )
         ENDIF
         CALL check( nf90_put_var(fncid,bme_id,RESHAPE((/REAL(
     $      eulbparmout),AIMAG(eulbparmout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,bml_id,RESHAPE((/REAL(
     $      lagbparmout),AIMAG(lagbparmout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,xm_id,RESHAPE((/REAL(
     $      -divxprpmout),AIMAG(-divxprpmout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,km_id,RESHAPE((/REAL(
     $      -curvmout),AIMAG(-curvmout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,be_id,RESHAPE(
     $      (/REAL(eulbparfout),-helicity*AIMAG(eulbparfout)/),
     $      (/mstep,mthsurf,2/))) )
         CALL check( nf90_put_var(fncid,bl_id,RESHAPE(
     $      (/REAL(lagbparfout),-helicity*AIMAG(lagbparfout)/),
     $      (/mstep,mthsurf,2/))) )
         CALL check( nf90_put_var(fncid,x_id,RESHAPE(
     $      (/REAL(-divxprpfout),-helicity*AIMAG(-divxprpfout)/),
     $      (/mstep,mthsurf,2/))))
         CALL check( nf90_put_var(fncid,k_id,RESHAPE(
     $      (/REAL(-curvfout),-helicity*AIMAG(-curvfout)/),
     $      (/mstep,mthsurf,2/))) )
         CALL check( nf90_put_var(fncid,b_id,equilbfun) )
         CALL check( nf90_close(fncid) )
      ENDIF

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_pmodb_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PMODB: "//
     $        "Components in perturbed mod b and del.x_prp"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a16,1x,a4,8(1x,a16))')"psi","m",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $        "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $        "imag(Bkxprp)"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            DO ipert=1,lmpert
               WRITE(out_unit,'(es17.8e3,1x,I4,8(es17.8e3))')
     $              psifac(istep),lmfac(ipert),
     $              REAL(eulbparmout(istep,ipert)),
     $              AIMAG(eulbparmout(istep,ipert)),
     $              REAL(lagbparmout(istep,ipert)),
     $              AIMAG(lagbparmout(istep,ipert)),
     $              REAL(-divxprpmout(istep,ipert)),
     $              AIMAG(-divxprpmout(istep,ipert)),
     $              REAL(-curvmout(istep,ipert)),
     $              AIMAG(-curvmout(istep,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (chebyshev_flag) THEN
         IF(verbose) WRITE(*,*)"Computing chebyshev for pmodb"
         ALLOCATE(chea(mpert,0:nche),chex(0:mstep))
         CALL cspline_alloc(chespl,mstep,1)
         chespl%xs=psifac
         chex=2.0*(psifac-0.5)
         DO ipert=1,mpert
            DO i=0,nche
               chespl%fs(:,1)=cos(REAL(i,r8)*acos(chex))/
     $              sqrt(1.0-chex**2.0)*lagbparmns(:,ipert)*4.0/pi
               CALL cspline_fit(chespl,"extrap")
               CALL cspline_int(chespl)
               chea(ipert,i)=chespl%fsi(mstep,1)
            ENDDO
         ENDDO

         chea(:,0)=chea(:,0)/2.0

         IF(ascii_flag)THEN
            CALL ascii_open(out_unit,"gpec_pmodb_chebyshev_n"//
     $        TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_PMODB_CHEBYSHEV: "//
     $           "Chebyshev coefficients for perturbed mod b"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(2(1x,a12,I4))')
     $           "nche =",nche,"mpert =",mpert
            WRITE(out_unit,*)
            DO ipert=1,mpert
               WRITE(out_unit,'(1x,a6,1x,I4)')"mfac =",mfac(ipert)
               WRITE(out_unit,*)
               WRITE(out_unit,'(1x,a6,2(1x,a16))')"chea",
     $              "real(chea)","imag(chea)"
               DO i=0,nche
                  WRITE(out_unit,'(1x,I6,2(es17.8e3))')i,
     $                 REAL(chea(ipert,i)),AIMAG(chea(ipert,i))
               ENDDO
               WRITE(out_unit,*)
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         cstep=200
         ALLOCATE(psis(cstep),ches(cstep))
         ALLOCATE(chelagbmns(cstep,mpert))
         chelagbmns=0
         psis=(/(istep,istep=0,cstep-1)/)/REAL(cstep-1,r8)*
     $        (psilim-psilow)+psilow
         ches=2.0*(psis-0.5)
         
         DO ipert=1,mpert
            DO i=0,nche
               chelagbmns(:,ipert)=chelagbmns(:,ipert)+
     $              chea(ipert,i)*cos(REAL(i,r8)*acos(ches))
            ENDDO
         ENDDO

         IF(ascii_flag)THEN
            CALL ascii_open(out_unit,"gpec_chepmodb_n"//
     $        TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CHEPMODB: "//
     $           "Reconstructed perturbed mod b by chebyshev"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
            WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $           "cstep =",cstep,"mpert =",mpert,"mthsurf =",mthsurf
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a16,1x,a4,2(1x,a16))')"psi","m",
     $           "real(lagb)","imag(lagb)"
            DO istep=1,cstep
               DO ipert=1,mpert
                  WRITE(out_unit,'(es17.8e3,1x,I4,2(es17.8e3))')
     $                 psis(istep),mfac(ipert),
     $                 REAL(chelagbmns(istep,ipert)),
     $                 AIMAG(chelagbmns(istep,ipert))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         CALL cspline_dealloc(chespl)
      ENDIF

      IF (fun_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_pmodb_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PMODB_FUN: "//
     $        "Components in perturbed mod b in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8,a13,es17.8e3)')
     $        "jac_out = ",jac_out,"R0 = ",ro
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(17(1x,a16))')"psi","theta","r","z",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $        "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $        "imag(Bkxprp)","equilb","dequilbdpsi","dequilbdtheta",
     $        "real(xms)","imag(xms)"
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            CALL gpeq_sol(psifac(istep))
            CALL gpeq_contra(psifac(istep))
            CALL gpeq_cova(psifac(istep))
            CALL iscdftb(mfac,mpert,xms_fun,mthsurf,xms_mn)
            DO itheta=0,mthsurf
               CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
               WRITE(out_unit,'(17(es17.8e3))')
     $              psifac(istep),theta(itheta),
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(eulbparfout(istep,itheta)),
     $              -helicity*AIMAG(eulbparfout(istep,itheta)),
     $              REAL(lagbparfout(istep,itheta)),
     $              -helicity*AIMAG(lagbparfout(istep,itheta)),
     $              REAL(-divxprpfout(istep,itheta)),
     $              -helicity*AIMAG(-divxprpfout(istep,itheta)),
     $              REAL(-curvfout(istep,itheta)),
     $              -helicity*AIMAG(-curvfout(istep,itheta)),
     $              equilbfun(istep,itheta),
     $              eqfun%fx(1),eqfun%fy(1),
     $              REAL(xms_fun(itheta)),
     $              -helicity*AIMAG(xms_fun(itheta))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
      
      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "pmodb.bin","UNKNOWN","REWIND","none")
         DO ipert=1,lmpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
     $              REAL(REAL(eulbparmout(istep,ipert)),4),
     $              REAL(AIMAG(eulbparmout(istep,ipert)),4),
     $              REAL(REAL(lagbparmout(istep,ipert)),4),
     $              REAL(AIMAG(lagbparmout(istep,ipert)),4)     
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"pmodb_2d.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf-1
         WRITE(bin_2d_unit)REAL(rs(9:mstep,1:mthsurf),4),
     $        REAL(zs(9:mstep,1:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(eulbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(eulbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(lagbparfout(9:mstep,1:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(lagbparfout(9:mstep,1:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF

      CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(rs,zs,equilbfun)
      DEALLOCATE(eulbparmns,lagbparmns,divxprpmns,curvmns,
     $   xmp1mns,xspmns,xmsmns,xmtmns,xmzmns)
      DEALLOCATE(eulbparmout,lagbparmout,divxprpmout,curvmout)
      DEALLOCATE(eulbparfun,lagbparfun,divxprpfun,curvfun,
     $   eulbparfout,lagbparfout,divxprpfout,curvfout)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_pmodb
c-----------------------------------------------------------------------
c     subprogram 9. gpout_xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_xbnormal(egnum,xspmn,spot,nspot)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nspot
      REAL(r8), INTENT(IN) :: spot
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: p_id,t_id,i_id,m_id,r_id,z_id,bm_id,b_id,
     $   wm_id,xm_id,x_id,rv_id,zv_id,rzstat

      INTEGER :: i,istep,ipert,iindex,itheta,tout
      REAL(r8) :: ileft,ximax,rmax,area

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs,dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xmns,ymns,
     $     xnomns,bnomns,bwpmns

      INTEGER :: mlow_pest,mhigh_pest,mpert_pest
      INTEGER, DIMENSION(:), ALLOCATABLE :: mfac_pest
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: pwpmns

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs,zs,psis,rvecs,zvecs
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rss,zss,
     $     xnofuns,bnofuns,intbwpmns

      COMPLEX(r8), DIMENSION(mpert) :: interpbwn
      TYPE(cspline_type) :: fsp_sol
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(rvecs(mstep,0:mthsurf),zvecs(mstep,0:mthsurf),
     $   rs(mstep,0:mthsurf),zs(mstep,0:mthsurf),psis(mstep,0:mthsurf),
     $   rss(mstep,0:mthsurf),zss(mstep,0:mthsurf),
     $   xnofuns(mstep,0:mthsurf),bnofuns(mstep,0:mthsurf))
      ALLOCATE(xmns(mstep,lmpert),ymns(mstep,lmpert),
     $   xnomns(mstep,lmpert),bnomns(mstep,lmpert),bwpmns(mstep,lmpert))
      IF (msing>0) ALLOCATE(intbwpmns(mstep,lmpert))
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b normal components"

      CALL idcon_build(egnum,xspmn)
      tout = tmag_out

      ! set up pest grid
      IF(TRIM(jac_out)=="pest")THEN ! will be mlow,mhigh if jac_type pest
        mlow_pest = lmlow
        mhigh_pest = lmhigh
      ELSE
        mlow_pest = MIN(lmlow,-64)
        mhigh_pest = MAX(lmhigh,64)
      ENDIF
      mpert_pest = ABS(mhigh_pest)+ABS(mlow_pest)+1
      ALLOCATE(mfac_pest(mpert_pest),pwpmns(mstep,mpert_pest))
      mfac_pest = (/(i,i=mlow_pest,mhigh_pest)/)
      pwpmns = 0

      CALL gpeq_alloc
      IF (msing>0) THEN
         CALL gpeq_interp_singsurf(fsp_sol,spot,nspot)
      ENDIF

      ! these surfaces are phsyically independent
      ! just repeating some generic geometric and coordinate conversion stuff in a loop
      ! filling in the xnomns, bnomns, xnofuns, bnofuns matrices
      ! The parallelization catch might be that we use these common fourier (iscdftb, iscdftf)
      ! and coordinate converting (gpeq_bcoordsout) subroutines? Will they confuse themselves in parallel?
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL gpeq_sol(psifac(istep))
         CALL gpeq_contra(psifac(istep))

         area=0
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
            dphi(itheta)=rzphi%f(3)
            rvecs(istep,itheta)=
     $           (cos(eta)*w(1,1)-sin(eta)*w(1,2))/delpsi(itheta)
            zvecs(istep,itheta)=
     $           (sin(eta)*w(1,1)+cos(eta)*w(1,2))/delpsi(itheta)
            area=area+jac*delpsi(itheta)/mthsurf
         ENDDO
         area=area-jac*delpsi(mthsurf)/mthsurf
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         xnofuns(istep,:)=xwp_fun/(jacs*delpsi)
         bnofuns(istep,:)=bwp_fun/(jacs*delpsi)
         CALL iscdftf(mfac,mpert,xnofuns(istep,:),mthsurf,xno_mn)
         CALL iscdftf(mfac,mpert,bnofuns(istep,:),mthsurf,bno_mn)
         IF (bwp_pest_flag) THEN
            ! distribute working decomposition on pest mrange
            DO i=1,mpert
               IF ((mlow-mlow_pest+i>=1).AND.
     $              (mlow-mlow_pest+i<=mpert_pest))
     $              pwpmns(istep,mlow-mlow_pest+i)=bno_mn(i)
            ENDDO
            ! convert to pest with magnetic angle
            CALL gpeq_bcoords(psifac(istep),pwpmns(istep,:),
     $           mfac_pest,mpert_pest,2,0,0,0,0,1)
         ENDIF            

         CALL gpeq_bcoordsout(xnomns(istep,:),xno_mn,psifac(istep),ji=0)
         CALL gpeq_bcoordsout(bnomns(istep,:),bno_mn,psifac(istep),ji=0)

         IF (msing>0) THEN
            CALL gpeq_interp_sol(fsp_sol,psifac(istep),interpbwn)
         ENDIF
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL gpeq_bcoordsout(bwpmns(istep,:),bno_mn,
     $                           psifac(istep),ji=1)
            IF (msing>0) THEN
               CALL gpeq_bcoordsout(intbwpmns(istep,:),interpbwn,
     $              psifac(istep),ji=0)
            ENDIF
         ELSE ! no need to re-weight bno_mn with expensive invfft and fft
            bwp_mn=bwp_mn/area
            bwpmns(istep,:)=0
            bwpmns(istep,mlow-lmlow+1:mlow-lmlow+mpert)=bwp_mn

            IF (msing>0) THEN
               interpbwn=interpbwn/area
               intbwpmns(istep,:)=0
               intbwpmns(istep,mlow-lmlow+1:mlow-lmlow+mpert)=interpbwn
            ENDIF
         ENDIF
         xnofuns(istep,:)=xnofuns(istep,:)*EXP(ifac*nn*dphi)
         bnofuns(istep,:)=bnofuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL gpeq_dealloc
      IF (msing>0) THEN
         CALL cspline_dealloc(fsp_sol)
      ENDIF

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_xbnormal_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBNORMAL: "//
     $        "Normal components of displacement and field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(xno)","imag(xno)","real(bno)","imag(bno)",
     $        "real(bwp)","imag(bwp)"
         DO istep=1,mstep,MAX(1,(mstep*lmpert-1)/max_linesout+1)
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $            psifac(istep),qfac(istep),lmfac(ipert),
     $            REAL(xnomns(istep,ipert)),AIMAG(xnomns(istep,ipert)),
     $            REAL(bnomns(istep,ipert)),AIMAG(bnomns(istep,ipert)),
     $            REAL(bwpmns(istep,ipert)),AIMAG(bwpmns(istep,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bwp_pest_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_bnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BNORMAL_PEST: "//
     $        "Normal components of field in pest coordinates"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",mpert_pest,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         
         DO istep=1,mstep,MAX(1,(mstep*lmpert-1)/max_linesout+1)
            DO ipert=1,mpert_pest
               WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $              psifac(istep),qfac(istep),mfac_pest(ipert),
     $              REAL(pwpmns(istep,ipert)),AIMAG(pwpmns(istep,ipert))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (fun_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_xbnormal_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"GPEC_XBNORMAL_FUN: "//
     $        "Normal components of displacement and field in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(10(1x,a16))')"psi","theta","r","z","rvec",
     $        "zvec","real(xno)","imag(xno)","real(bno)","imag(bno)"
         
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            DO itheta=0,mthsurf
               WRITE(out_unit,'(10(es17.8e3))')
     $              psifac(istep),theta(itheta),
     $              rs(istep,itheta),zs(istep,itheta),
     $              rvecs(istep,itheta),zvecs(istep,itheta),
     $              REAL(xnofuns(istep,itheta)),
     $              -helicity*AIMAG(xnofuns(istep,itheta)),
     $              REAL(bnofuns(istep,itheta)),
     $              -helicity*AIMAG(bnofuns(istep,itheta))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "xbnormal.bin","UNKNOWN","REWIND","none")
         IF (msing>0) THEN
            DO ipert=1,lmpert
               DO istep=1,mstep
                  WRITE(bin_unit)REAL(psifac(istep),4),
     $                 REAL(REAL(xnomns(istep,ipert)),4),
     $                 REAL(AIMAG(xnomns(istep,ipert)),4),
     $                 REAL(REAL(bnomns(istep,ipert)),4),
     $                 REAL(AIMAG(bnomns(istep,ipert)),4),
     $                 REAL(REAL(bwpmns(istep,ipert)),4),
     $                 REAL(AIMAG(bwpmns(istep,ipert)),4),
     $                 REAL(REAL(intbwpmns(istep,ipert)),4),
     $                 REAL(AIMAG(intbwpmns(istep,ipert)),4)
               ENDDO
               WRITE(bin_unit)
            ENDDO
         ELSE
            DO ipert=1,lmpert
               DO istep=1,mstep
                  WRITE(bin_unit)REAL(psifac(istep),4),
     $                 REAL(REAL(xnomns(istep,ipert)),4),
     $                 REAL(AIMAG(xnomns(istep,ipert)),4),
     $                 REAL(REAL(bnomns(istep,ipert)),4),
     $                 REAL(AIMAG(bnomns(istep,ipert)),4),
     $                 REAL(REAL(bwpmns(istep,ipert)),4),
     $                 REAL(AIMAG(bwpmns(istep,ipert)),4)
               ENDDO
               WRITE(bin_unit)
            ENDDO
         ENDIF
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"xbnormal_2d.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:mstep,0:mthsurf),4),
     $        REAL(zs(9:mstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(xnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(xnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(bnofuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(bnofuns(9:mstep,0:mthsurf)),4)

         CALL bin_close(bin_2d_unit)
         
         ximax=MAXVAL(ABS(xnofuns))
         CALL bicube_eval(rzphi,psilim,theta(0),0)
         rmax=SQRT(rzphi%f(1))
         xnofuns=xnofuns/ximax*rmax/6.0

         IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $     "Maximum displacement = ",ximax
         IF(verbose) WRITE(*,'(1x,a,es10.3)')
     $     "Scale factor for 2d plots = ",rmax/(ximax*6.0)

         rss=CMPLX(rs,rs)+xnofuns*rvecs
         zss=CMPLX(zs,zs)+xnofuns*zvecs
         DO itheta=0,mthsurf
            psis(:,itheta)=psifac(:)
         ENDDO

         CALL bin_open(bin_2d_unit,
     $        "pflux_re_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(REAL(rss(9:mstep,0:mthsurf)),4),
     $        REAL(REAL(zss(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psis(9:mstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         CALL bin_open(bin_2d_unit,
     $        "pflux_im_2d.bin","UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)
     $        REAL(-helicity*AIMAG(rss(9:mstep,0:mthsurf)),4),
     $        REAL(-helicity*AIMAG(zss(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(psis(9:mstep,0:mthsurf),4)
         CALL bin_close(bin_2d_unit)

         DO istep=1,mstep
            xmns(istep,:)=lmfac
            ymns(istep,:)=psifac(istep)
         ENDDO
         CALL bin_open(bin_2d_unit,"bnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,lmpert-1
         WRITE(bin_2d_unit)REAL(xmns(9:mstep,:),4),
     $        REAL(ymns(9:mstep,:),4)
         WRITE(bin_2d_unit)REAL(ABS(bwpmns(9:mstep,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"m_out",m_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_inq_dimid(fncid,"theta_dcon",t_id) )
         CALL check( nf90_redef(fncid))
         rzstat = nf90_inq_varid(fncid, "R", r_id) ! check if R,z already stored
         IF(rzstat/=nf90_noerr)THEN
            CALL check( nf90_def_var(fncid, "R", nf90_double,
     $                  (/p_id,t_id/),r_id) )
            CALL check( nf90_put_att(fncid,r_id,"long_name",
     $                  "Major radius") )
            CALL check( nf90_put_att(fncid,r_id,"units","m") )
            CALL check( nf90_def_var(fncid, "z", nf90_double,
     $                  (/p_id,t_id/),z_id) )
            CALL check( nf90_put_att(fncid,z_id,"long_name",
     $                  "Vertical position") )
            CALL check( nf90_put_att(fncid,z_id,"units","m") )
         ENDIF
         CALL check( nf90_def_var(fncid, "b_n_fun", nf90_double,
     $               (/p_id,t_id,i_id/),b_id) )
         CALL check( nf90_put_att(fncid,b_id,"long_name",
     $               "Perturbed field normal to the flux surface") )
         CALL check( nf90_put_att(fncid,b_id,"units","Tesla") )
         CALL check( nf90_def_var(fncid, "b_n", nf90_double,
     $               (/p_id,m_id,i_id/),bm_id) )
         CALL check( nf90_put_att(fncid,bm_id,"long_name",
     $               "Perturbed field normal to the flux surface") )
         CALL check( nf90_put_att(fncid,bm_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,bm_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "Jbgradpsi", nf90_double,
     $               (/p_id,m_id,i_id/),wm_id) )
         CALL check( nf90_put_att(fncid,wm_id,"long_name",
     $      "Jaconbian weighted contravariant psi component of "//
     $      "the perturbed field") )
         CALL check( nf90_put_att(fncid,wm_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,wm_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "xi_n_fun", nf90_double,
     $               (/p_id,t_id,i_id/),x_id) )
         CALL check( nf90_put_att(fncid,x_id,"long_name",
     $               "Displacement normal to the flux surface") )
         CALL check( nf90_put_att(fncid,x_id,"units","m") )
         CALL check( nf90_def_var(fncid, "xi_n", nf90_double,
     $               (/p_id,m_id,i_id/),xm_id) )
         CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $               "Displacement normal to the flux surface") )
         CALL check( nf90_put_att(fncid,xm_id,"units","m") )
         CALL check( nf90_put_att(fncid,xm_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "R_n", nf90_double,
     $               (/p_id,t_id/),rv_id) )
         CALL check( nf90_put_att(fncid,rv_id,"long_name",
     $               "Radial unit normal: R_3D = R_sym+x_n*R_n") )
         CALL check( nf90_def_var(fncid, "z_n", nf90_double,
     $               (/p_id,t_id/),zv_id) )
         CALL check( nf90_put_att(fncid,zv_id,"long_name",
     $               "Vertical unit normal: z_3D = z_sym+x_n*z_n") )
         CALL check( nf90_enddef(fncid) )
         IF(rzstat/=nf90_noerr)THEN
            CALL check( nf90_put_var(fncid,r_id,rs) )
            CALL check( nf90_put_var(fncid,z_id,zs) )
         ENDIF
         CALL check( nf90_put_var(fncid,xm_id,RESHAPE((/REAL(xnomns),
     $                AIMAG(xnomns)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,bm_id,RESHAPE((/REAL(bnomns),
     $                AIMAG(bnomns)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,wm_id,RESHAPE((/REAL(bwpmns),
     $                AIMAG(bwpmns)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,x_id,RESHAPE((/REAL(xnofuns),
     $               -helicity*AIMAG(xnofuns)/),(/mstep,mthsurf+1,2/))))
         CALL check( nf90_put_var(fncid,b_id,RESHAPE((/REAL(bnofuns),
     $               -helicity*AIMAG(bnofuns)/),(/mstep,mthsurf+1,2/))))
         CALL check( nf90_put_var(fncid,rv_id,rvecs) )
         CALL check( nf90_put_var(fncid,zv_id,zvecs) )
         CALL check( nf90_close(fncid) )
      ENDIF

      IF (flux_flag) THEN
         IF (.NOT. bin_2d_flag) THEN
            ximax=MAXVAL(ABS(xnofuns))
            CALL bicube_eval(rzphi,psilim,theta(0),0)
            rmax=SQRT(rzphi%f(1))
            xnofuns=xnofuns/ximax*rmax/6.0
         ENDIF
            
         rss=xnofuns*rvecs
         zss=xnofuns*zvecs
         DO itheta=0,mthsurf
            psis(:,itheta)=psifac(1:mstep)
         ENDDO

         IF(ascii_flag)THEN
            CALL ascii_open(out_unit,"gpec_xbnormal_flux_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_XBNROMAL_FLUX: "//
     $           "Perturbed flux surfaces"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
            WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $           "mstep =",mstep,"mthsurf =",mthsurf
            WRITE(out_unit,'(1x,a12,es17.8e3)')"scale =",rmax/(ximax*6)
            WRITE(out_unit,*)
            WRITE(out_unit,'(7(1x,a16))')"r","z","real(r)","imag(r)",
     $           "real(z)","imag(z)","psi"
            DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
               DO itheta=0,mthsurf
                  WRITE(out_unit,'(7(es17.8e3))')
     $                 rs(istep,itheta),zs(istep,itheta),
     $                 REAL(rss(istep,itheta)),AIMAG(rss(istep,itheta)),
     $                 REAL(zss(istep,itheta)),AIMAG(zss(istep,itheta)),
     $                 psis(istep,itheta)
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF
      ENDIF         
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(rvecs,zvecs,rs,zs,psis,rss,zss,xnofuns,bnofuns)
      DEALLOCATE(xmns,ymns,xnomns,bnomns,bwpmns)
      IF (msing>0) DEALLOCATE(intbwpmns)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_xbnormal
c-----------------------------------------------------------------------
c     subprogram 10. gpout_vbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_vbnormal(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      
      INTEGER :: ipsi,ipert,i,iindex
      REAL(r8) :: ileft
      REAL(r8), DIMENSION(0:cmpsi) :: psi
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: vcmn

      INTEGER :: p_id,m_id,t_id,i_id,mp_id,
     $    mpv_id,qs_id,vn_id,vw_id,pw_id

      REAL(r8), DIMENSION(cmpsi) :: qs
      REAL(r8), DIMENSION(cmpsi,lmpert) :: xmns,ymns
      COMPLEX(r8), DIMENSION(cmpsi,lmpert) :: vnomns,vwpmns
      COMPLEX(r8), DIMENSION(mstep,lmpert) :: vnomns_mstep,vwpmns_mstep

      INTEGER :: mlow_pest,mhigh_pest,mpert_pest
      INTEGER, DIMENSION(:), ALLOCATABLE :: mfac_pest
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: pwpmns, pwpmns_mstep

      logical :: debug_flag = .true.
      TYPE(cspline_type) :: vnos, vwps, pwps
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vb normal components from coils"

      IF(TRIM(jac_out)=="pest")THEN ! will be mlow,mhigh if jac_type pest
        mlow_pest = lmlow
        mhigh_pest = lmhigh
      ELSE
        mlow_pest = MIN(lmlow,-64)
        mhigh_pest = MAX(lmhigh,64)
      ENDIF
      mpert_pest = ABS(mhigh_pest)+ABS(mlow_pest)+1
      ALLOCATE( mfac_pest(mpert_pest), pwpmns(cmpsi,mpert_pest),
     $    pwpmns_mstep(mstep,mpert_pest) )
      mfac_pest = (/(i,i=mlow_pest,mhigh_pest)/)
      pwpmns=0

      psi=(/(ipsi,ipsi=0,cmpsi)/)/REAL(cmpsi,r8)
      psi=psilow+(psilim-psilow)*SIN(psi*pi/2)**2
      ALLOCATE(vcmn(cmpert))
      vcmn=0
      vnomns=0
      vwpmns=0

      DO ipsi=1,cmpsi
         iindex = FLOOR(REAL(ipsi,8)/FLOOR(cmpsi/10.0))*10
         ileft = REAL(ipsi,8)/FLOOR(cmpsi/10.0)*10-iindex
         IF ((ipsi-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a24)')
     $        "volume = ",iindex,"% vacuum b computations"
         CALL spline_eval(sq,psi(ipsi),0)
         qs(ipsi)=sq%f(4)
         CALL field_bs_psi(psi(ipsi),vcmn,0)

         IF (bwp_pest_flag) THEN
            DO i=1,cmpert
               IF ((cmlow-mlow_pest+i>=1).AND.
     $              (cmlow-mlow_pest+i<=mpert_pest)) THEN
                  pwpmns(ipsi,cmlow-mlow_pest+i)=vcmn(i)
               ENDIF
            ENDDO
            CALL gpeq_bcoords(psi(ipsi),pwpmns(ipsi,:),mfac_pest,
     $           mpert_pest,2,0,0,0,0,1)
         ENDIF
         
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  vwpmns(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
            vnomns(ipsi,:)=vwpmns(ipsi,:)
            CALL gpeq_bcoords(psi(ipsi),vnomns(ipsi,:),lmfac,lmpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL gpeq_bcoords(psi(ipsi),vwpmns(ipsi,:),lmfac,lmpert,
     $           rout,bpout,bout,rcout,tout,1)
         ELSE
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  vnomns(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
            CALL field_bs_psi(psi(ipsi),vcmn,2)
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  vwpmns(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(vcmn)

      ! put onto full radial grid
      CALL cspline_alloc(vnos,cmpsi-1,lmpert)
      vnos%xs(0:) = psi(1:)
      vnos%fs(0:,:) = vnomns(:,:)
      call cspline_fit(vnos,"extrap")
      CALL cspline_alloc(vwps,cmpsi-1,lmpert)
      vwps%xs(0:) = psi(1:)
      vwps%fs(0:,:) = vwpmns(:,:)
      call cspline_fit(vwps,"extrap")
      CALL cspline_alloc(pwps,cmpsi-1,mpert_pest)
      pwps%xs(0:) = psi(1:)
      pwps%fs(0:,:) = pwpmns(:,:)
      call cspline_fit(pwps,"extrap")
      DO ipsi=1,mstep
         CALL cspline_eval(vnos,psifac(ipsi),0)
         CALL cspline_eval(vwps,psifac(ipsi),0)
         CALL cspline_eval(pwps,psifac(ipsi),0)
         vnomns_mstep(ipsi,:) = vnos%f(:)
         vwpmns_mstep(ipsi,:) = vwps%f(:)
         pwpmns_mstep(ipsi,:) = pwps%f(:)
      ENDDO

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"m_out",m_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_inq_dimid(fncid,"theta_dcon",t_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "b_n_x", nf90_double,
     $      (/p_id,m_id,i_id/),vn_id) )
         CALL check( nf90_put_att(fncid,vn_id,"long_name",
     $      "Externally applied normal field") )
         CALL check( nf90_put_att(fncid,vn_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,vn_id,"jacobian",jac_out) )
         CALL check( nf90_def_var(fncid, "Jbgradpsi_x", nf90_double,
     $      (/p_id,m_id,i_id/),vw_id) )
         CALL check( nf90_put_att(fncid,vw_id,"long_name",
     $      "Jacobian weighted contravariant psi component "//
     $      "of the applied field") )
         CALL check( nf90_put_att(fncid,vw_id,"units","Tesla") )
         CALL check( nf90_put_att(fncid,vw_id,"jacobian",jac_out) )
         IF(TRIM(jac_out)/="pest" .and. bwp_pest_flag)THEN
            CALL check( nf90_def_dim(fncid,"m_pest",mpert_pest,mp_id) )
            CALL check( nf90_def_var(fncid, "m_pest", nf90_int, mp_id,
     $         mpv_id) )
            CALL check( nf90_def_var(fncid, "Jbgradpsi_x_pest",
     $         nf90_double, (/p_id, mp_id, i_id/), pw_id) )
            CALL check( nf90_put_att(fncid,pw_id,"long_name",
     $         "Jacobian weighted contravariant psi component "//
     $         "of the applied field") )
            CALL check( nf90_put_att(fncid,pw_id,"units","Tesla") )
            CALL check( nf90_put_att(fncid,pw_id,"jacobian","pest") )
         ENDIF
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid,vn_id,RESHAPE((/REAL(
     $      vnomns_mstep), AIMAG(vnomns_mstep)/), (/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,vw_id,RESHAPE((/REAL(
     $      vwpmns_mstep), AIMAG(vwpmns_mstep)/), (/mstep,lmpert,2/))) )
         IF(TRIM(jac_out)/="pest")THEN
            CALL check( nf90_put_var(fncid,mpv_id,mfac_pest) )
            CALL check( nf90_put_var(fncid,pw_id,RESHAPE(
     $         (/REAL(pwpmns_mstep),AIMAG(pwpmns_mstep)/),
     $         (/mstep,mpert_pest,2/))) )
         ENDIF
         CALL check( nf90_close(fncid) )
      ENDIF

      ! write to ascii table
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_vbnormal_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBNROMAL: "//
     $        "Normal components of coil field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mpsi =",cmpsi,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $        "real(bno)","imag(bno)","real(bwp)","imag(bwp)"
         DO ipsi=1,cmpsi
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,4(es17.8e3))')
     $              psi(ipsi),qs(ipsi),lmfac(ipert),
     $              REAL(vnomns(ipsi,ipert)),AIMAG(vnomns(ipsi,ipert)),
     $              REAL(vwpmns(ipsi,ipert)),AIMAG(vwpmns(ipsi,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bwp_pest_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_vbnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBNROMAL_PEST: "//
     $        "Normal components of coil field in pest coordinates"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mpsi =",cmpsi,"mpert =",mpert_pest,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,2(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         DO ipsi=1,cmpsi
            DO ipert=1,mpert_pest
               WRITE(out_unit,'(2(es17.8e3),1x,I4,2(es17.8e3))')
     $              psi(ipsi),qs(ipsi),mfac_pest(ipert),
     $              REAL(pwpmns(ipsi,ipert)),AIMAG(pwpmns(ipsi,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "vbnormal.bin","UNKNOWN","REWIND","none")
         DO ipert=1,lmpert
            DO ipsi=1,cmpsi
               WRITE(bin_unit)REAL(psi(ipsi),4),
     $              REAL(REAL(vnomns(ipsi,ipert)),4),
     $              REAL(AIMAG(vnomns(ipsi,ipert)),4),
     $              REAL(REAL(vwpmns(ipsi,ipert)),4),
     $              REAL(AIMAG(vwpmns(ipsi,ipert)),4)
               
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)

         DO ipsi=1,cmpsi
            xmns(ipsi,:)=lmfac
            ymns(ipsi,:)=psi(ipsi)
         ENDDO
         CALL bin_open(bin_2d_unit,"vbnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)cmpsi-1,lmpert-1
         WRITE(bin_2d_unit)REAL(xmns(1:cmpsi,:),4),
     $        REAL(ymns(1:cmpsi,:),4)
         WRITE(bin_2d_unit)REAL(ABS(vwpmns(1:cmpsi,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      RETURN
      END SUBROUTINE gpout_vbnormal
c-----------------------------------------------------------------------
c     subprogram 11. gpout_xbtangent.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_xbtangent(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs,zs,psis,
     $     rvecs,zvecs,vecs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,bs,dphi
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xtamns,btamns
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwt_fun,xvt_fun,xvz_fun,
     $     bwt_fun,bvt_fun,bvz_fun,xta_fun,bta_fun
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rss,zss,
     $     xtafuns,btafuns
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf),
     $   psis(mstep,0:mthsurf),vecs(mstep,0:mthsurf),
     $   rvecs(mstep,0:mthsurf),zvecs(mstep,0:mthsurf))
      ALLOCATE(xtamns(mstep,lmpert),btamns(mstep,lmpert))
      ALLOCATE(rss(mstep,0:mthsurf),zss(mstep,0:mthsurf),
     $   xtafuns(mstep,0:mthsurf),btafuns(mstep,0:mthsurf))
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b tangential components"

      CALL idcon_build(egnum,xspmn)

      CALL gpeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL gpeq_sol(psifac(istep))
         CALL gpeq_contra(psifac(istep))
         CALL gpeq_cova(psifac(istep))

         DO itheta=0,mthsurf
            CALL bicube_eval(eqfun,psifac(istep),theta(itheta),0)
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jacs(itheta)=rzphi%f(4)
            bs(itheta)=eqfun%f(1)
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            rvecs(istep,itheta)=v(2,1)*cos(eta)-v(2,2)*sin(eta)
            zvecs(istep,itheta)=v(2,1)*sin(eta)+v(2,2)*cos(eta)
         ENDDO
         vecs(istep,:)=sqrt(rvecs(istep,:)**2+zvecs(istep,:)**2)
         rvecs(istep,:)=rvecs(istep,:)/vecs(istep,:)
         zvecs(istep,:)=zvecs(istep,:)/vecs(istep,:)
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,xvt_fun,mthsurf,xvt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
         CALL iscdftb(mfac,mpert,bvt_fun,mthsurf,bvt_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)
         xta_fun=xwt_fun/jacs-(chi1/(jacs*bs))**2*(xvt_fun+q*xvz_fun)
         bta_fun=bwt_fun/jacs-(chi1/(jacs*bs))**2*(bvt_fun+q*bvz_fun)
         xtafuns(istep,:)=xta_fun*vecs(istep,:)
         btafuns(istep,:)=bta_fun*vecs(istep,:)
         CALL iscdftf(mfac,mpert,xtafuns(istep,:),mthsurf,xta_mn)
         CALL iscdftf(mfac,mpert,btafuns(istep,:),mthsurf,bta_mn)
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL gpeq_bcoords(psifac(istep),xta_mn,mfac,lmpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL gpeq_bcoords(psifac(istep),bta_mn,mfac,lmpert,
     $           rout,bpout,bout,rcout,tout,0)
         ENDIF
         xtamns(istep,:)=xta_mn
         btamns(istep,:)=bta_mn
         xtafuns(istep,:)=xtafuns(istep,:)*EXP(ifac*nn*dphi)
         btafuns(istep,:)=btafuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL gpeq_dealloc

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_xbtangent_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBNROMAL: "//
     $        "Tangential components of displacement and field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $        "real(xta)","imag(xta)","real(bta)","imag(bta)"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,4(es17.8e3))')
     $            psifac(istep),qfac(istep),lmfac(ipert),
     $            REAL(xtamns(istep,ipert)),AIMAG(xtamns(istep,ipert)),
     $            REAL(btamns(istep,ipert)),AIMAG(btamns(istep,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (fun_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_xbtangent_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"GPEC_XBTANGENT_FUN: "//
     $        "Tangential components of displacement "//
     $        "and field in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(8(1x,a16))')"r","z","rvec","zvec",
     $        "real(xta)","imag(xta)","real(bta)","imag(bta)"
         
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            DO itheta=0,mthsurf
               WRITE(out_unit,'(8(es17.8e3))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              rvecs(istep,itheta),zvecs(istep,itheta),
     $              REAL(xtafuns(istep,itheta)),
     $              AIMAG(xtafuns(istep,itheta)),
     $              REAL(btafuns(istep,itheta)),
     $              AIMAG(btafuns(istep,itheta))        
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "xbtangent.bin","UNKNOWN","REWIND","none")
         DO ipert=1,lmpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
     $              REAL(REAL(xtamns(istep,ipert)),4),
     $              REAL(AIMAG(xtamns(istep,ipert)),4),
     $              REAL(REAL(btamns(istep,ipert)),4),
     $              REAL(AIMAG(btamns(istep,ipert)),4)
            ENDDO
            WRITE(bin_unit)
         ENDDO
         CALL bin_close(bin_unit)
      ENDIF

      IF (bin_2d_flag) THEN
         CALL bin_open(bin_2d_unit,"xbtangent_2d.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)mstep-9,mthsurf
         WRITE(bin_2d_unit)REAL(rs(9:mstep,0:mthsurf),4),
     $        REAL(zs(9:mstep,0:mthsurf),4)
         WRITE(bin_2d_unit)REAL(REAL(xtafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(xtafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(REAL(btafuns(9:mstep,0:mthsurf)),4)
         WRITE(bin_2d_unit)REAL(AIMAG(btafuns(9:mstep,0:mthsurf)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF         
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(rs,zs,psis,vecs,rvecs,zvecs)
      DEALLOCATE(xtamns,btamns)
      DEALLOCATE(rss,zss,xtafuns,btafuns)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_xbtangent
c-----------------------------------------------------------------------
c     subprogram 12. gpout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 10. gpout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_xbrzphi(egnum,xspmn,nr,nz,bnimn,bnomn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnimn,bnomn

      INTEGER :: i,j,k,l,iindex,np
      REAL(r8) :: mid,btlim,rlim,ileft,delr,delz,cha,chb,chc,chd,
     $   rij,zij,t11,t12,t21,t22,t33
      COMPLEX(r8) :: xwp,bwp,xwt,bwt,xvz,bvz

      INTEGER :: r_id,z_id,i_id,xr_id,xz_id,xp_id,br_id,bz_id,bp_id,
     $   bre_id,bze_id,bpe_id,brp_id,bzp_id,bpp_id

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz,ebr,ebz,ebp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: xrr,xrz,xrp,brr,brz,brp,
     $     bpr,bpz,bpp,vbr,vbz,vbp,vpbr,vpbz,vpbp,vvbr,vvbz,vvbp,
     $     btr,btz,btp,vcbr,vcbz,vcbp,xcr,xcz,xcp,bcr,bcz,bcp

      REAL(r8), DIMENSION(:), ALLOCATABLE :: chex,chey
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: chear,cheaz,
     $     chxar,chxaz

c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      ! initialization
      ebr = 0
      ebz = 0
      ebp = 0
      xrr = 0
      xrz = 0
      xrp = 0
      xcr = 0
      xcz = 0
      xcp = 0
      brr = 0
      brz = 0
      brp = 0
      btr = 0
      btz = 0
      btp = 0
      bcr = 0
      bcz = 0
      bcp = 0
      bpr = 0
      bpz = 0
      bpp = 0
      vbr = 0
      vbz = 0
      vbp = 0
      vcbr = 0
      vcbz = 0
      vcbp = 0
      vpbr = 0
      vpbz = 0
      vpbp = 0
      vvbr = 0
      vvbz = 0
      vvbp = 0

      CALL idcon_build(egnum,xspmn)

      IF (mode_flag) THEN
         CALL gpeq_alloc
         CALL gpeq_sol(psilim)
         bnomn=bwp_mn
         CALL gpeq_dealloc
      ENDIF

      IF (eqbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing equilibrium magnetic fields"
         IF(timeit) CALL gpec_timer(-2)
         ! evaluate f value for vacuum
         mid = 0.0
         CALL spline_eval(sq,psilim,0)
         CALL bicube_eval(rzphi,psilim,mid,0)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*rzphi%f(2)
         rlim=ro+rfac*COS(eta)
         btlim = abs(sq%f(1))/(twopi*rlim)
         DO i=0,nr
            DO j=0,nz
               CALL bicube_eval(psi_in,gdr(i,j),gdz(i,j),1)
               ebr(i,j) = -psi_in%fy(1)/gdr(i,j)*psio
               ebz(i,j) = psi_in%fx(1)/gdr(i,j)*psio
               IF (gdl(i,j) == 1) THEN  
                  CALL spline_eval(sq,gdpsi(i,j),0)
                  ebp(i,j) = abs(sq%f(1))/(twopi*gdr(i,j))
               ELSE
                  ebp(i,j) = btlim*rlim/gdr(i,j)  
               ENDIF
            ENDDO   
         ENDDO

         IF(ipd>0)THEN
            ebr=-ebr
            ebz=-ebz
         ENDIF
         IF(btd<0)ebp=-ebp
         IF(timeit) CALL gpec_timer(2)
      ENDIF

      CALL gpeq_alloc

      IF(verbose) WRITE(*,*)"Mapping fields to cylindrical coordinates"
      IF(timeit) CALL gpec_timer(-2)

      IF (brzphi_flag .OR. xrzphi_flag) THEN
         DO i=0,nr
         iindex = FLOOR(REAL(i,8)/FLOOR((nr-1)/10.0))*10
         ileft = REAL(i,8)/FLOOR((nr-1)/10.0)*10-iindex
         IF ((i /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a10)')
     $        "volume = ",iindex,"% mappings"
            DO j=0,nz
               IF (gdl(i,j)==1) THEN
                  CALL gpeq_sol(gdpsi(i,j))
                  CALL gpeq_contra(gdpsi(i,j))
                  CALL gpeq_cova(gdpsi(i,j))
                  ! compute matric tensor components.
                  CALL bicube_eval(rzphi,gdpsi(i,j),gdthe(i,j),1)
                  rfac=SQRT(rzphi%f(1))
                  eta=twopi*(gdthe(i,j)+rzphi%f(2))
                  rij=ro+rfac*COS(eta)
                  zij=zo+rfac*SIN(eta)
                  jac=rzphi%f(4)
                  v(1,1)=rzphi%fx(1)/(2*rfac)
                  v(1,2)=rzphi%fx(2)*twopi*rfac
                  v(2,1)=rzphi%fy(1)/(2*rfac)
                  v(2,2)=(1+rzphi%fy(2))*twopi*rfac
                  v(3,3)=twopi*rij
                  t11=cos(eta)*v(1,1)-sin(eta)*v(1,2)
                  t12=cos(eta)*v(2,1)-sin(eta)*v(2,2)
                  t21=sin(eta)*v(1,1)+cos(eta)*v(1,2)
                  t22=sin(eta)*v(2,1)+cos(eta)*v(2,2)
                  t33=-1.0/v(3,3)
                  ! three vector components.
                  xwp = SUM(xwp_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bwp = SUM(bwp_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xwt = SUM(xwt_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bwt = SUM(bwt_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xvz = SUM(xvz_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  bvz = SUM(bvz_mn*EXP(ifac*mfac*twopi*gdthe(i,j)))
                  xrr(i,j) = (t11*xwp+t12*xwt)/jac
                  brr(i,j) = (t11*bwp+t12*bwt)/jac
                  xrz(i,j) = (t21*xwp+t22*xwt)/jac
                  brz(i,j) = (t21*bwp+t22*bwt)/jac
                  xrp(i,j) = t33*xvz
                  brp(i,j) = t33*bvz
                  ! machine toroidal angle ! it matters.
                  xrr(i,j) = xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brr(i,j) = brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrz(i,j) = xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brz(i,j) = brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrp(i,j) = xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brp(i,j) = brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  ! Correction for helicity
                  IF(helicity>0)THEN
                     xrr(i,j) = CONJG(xrr(i,j))
                     brr(i,j) = CONJG(brr(i,j))
                     xrz(i,j) = CONJG(xrz(i,j))
                     brz(i,j) = CONJG(brz(i,j))
                     xrp(i,j) = CONJG(xrp(i,j))
                     brp(i,j) = CONJG(brp(i,j))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      IF(timeit) CALL gpec_timer(2)

      CALL gpeq_dealloc

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)
     $      "Computing vacuum fields by surface currents"
         CALL gpvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
         IF (helicity<0) THEN
            vbr=CONJG(vbr)
            vbz=CONJG(vbz)
            vbp=CONJG(vbp)
         ENDIF
         DO i=0,nr
            DO j=0,nz
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  brr(i,j)=vbr(i,j)
                  brz(i,j)=vbz(i,j)
                  brp(i,j)=vbp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
         IF(timeit) CALL gpec_timer(2)
      ENDIF
      
      IF (brzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL gpvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF

         IF (coil_flag) THEN
            IF(verbose) WRITE(*,*)"Computing vacuum fields by coils"
            np=nn*48 ! make it consistent with cmzeta later.
            CALL field_bs_rzphi(nr,nz,np,gdr,gdz,vcbr,vcbz,vcbp)
            IF (divzero_flag) THEN
               CALL gpeq_rzpdiv(nr,nz,gdr,gdz,vcbr,vcbz,vcbp)
            ENDIF
            IF (div_flag) THEN
               CALL gpdiag_rzpdiv(nr,nz,gdl,gdr,gdz,
     $              vcbr,vcbz,vcbp,"c")
            ENDIF
         ENDIF
         
         DO i=0,nr
            DO j=0,nz                  
               IF (gdl(i,j)/=1) THEN
                  gdl(i,j)=vgdl(i,j)
                  bpr(i,j)=vpbr(i,j)
                  bpz(i,j)=vpbz(i,j)
                  bpp(i,j)=vpbp(i,j)
                  btr(i,j)=vpbr(i,j)+vcbr(i,j)
                  btz(i,j)=vpbz(i,j)+vcbz(i,j)
                  btp(i,j)=vpbp(i,j)+vcbp(i,j)
               ELSE
                  bpr(i,j)=brr(i,j)-vcbr(i,j)
                  bpz(i,j)=brz(i,j)-vcbz(i,j)
                  bpp(i,j)=brp(i,j)-vcbp(i,j)
                  btr(i,j)=brr(i,j)
                  btz(i,j)=brz(i,j)
                  btp(i,j)=brp(i,j)
               ENDIF
               
            ENDDO
         ENDDO
         IF (divzero_flag) THEN
            CALL gpeq_rzpdiv(nr,nz,gdr,gdz,btr,btz,btp)
         ENDIF
         IF (div_flag) THEN
            CALL gpdiag_rzpdiv(nr,nz,gdl,gdr,gdz,btr,btz,btp,"b")
         ENDIF
      ENDIF

      IF (chebyshev_flag) THEN
         IF(verbose) WRITE(*,*)"Computing chebyshev for xbrzphi"
         ALLOCATE(chex(0:nr),chey(0:nz),
     $        chear(0:nchr,0:nchz),cheaz(0:nchr,0:nchz),
     $        chxar(0:nchr,0:nchz),chxaz(0:nchr,0:nchz))
         delr = gdr(1,0)-gdr(0,0)
         delz = gdz(0,1)-gdz(0,0)
         cha = (gdr(nr,0)-gdr(0,0))/2.0+delr
         chb = (gdr(nr,0)+gdr(0,0))/2.0
         chc = (gdz(0,nz)-gdz(0,0))/2.0+delz
         chd = (gdz(0,nz)+gdz(0,0))/2.0
         chex = (gdr(:,0)-chb)/cha
         chey = (gdz(0,:)-chd)/chc
         chear=0
         cheaz=0
         chxar=0
         chxaz=0

         DO i=0,nchr
            DO j=0,nchz
               DO k=0,nr
                  DO l=0,nz
                     chear(i,j)=chear(i,j)+btr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     cheaz(i,j)=cheaz(i,j)+btz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                     chxar(i,j)=chxar(i,j)+xrr(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))
                     chxaz(i,j)=chxaz(i,j)+xrz(k,l)*
     $                    cos(REAL(i,r8)*acos(chex(k)))*
     $                    cos(REAL(j,r8)*acos(chey(l)))/
     $                    sqrt((1.0-chex(k)**2.0)*(1.0-chey(l)**2.0))                    
                  ENDDO
               ENDDO
               chear(i,j)=chear(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               cheaz(i,j)=cheaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxar(i,j)=chxar(i,j)*4*delr*delz/(cha*chc*pi**2.0)
               chxaz(i,j)=chxaz(i,j)*4*delr*delz/(cha*chc*pi**2.0)
            ENDDO
         ENDDO
         chear(0,:)=chear(0,:)/2.0
         cheaz(0,:)=cheaz(0,:)/2.0        
         chear(:,0)=chear(:,0)/2.0
         cheaz(:,0)=cheaz(:,0)/2.0  
         chxar(0,:)=chxar(0,:)/2.0
         chxaz(0,:)=chxaz(0,:)/2.0        
         chxar(:,0)=chxar(:,0)/2.0
         chxaz(:,0)=chxaz(:,0)/2.0   

         IF(ascii_flag)THEN
            CALL ascii_open(out_unit,"gpec_brzphi_chebyshev_n"//
     $        TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_BRZPHI_CHEBYSHEV: "//
     $           "Chebyshev coefficients for brzphi field"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(2(1x,a12,I4))')
     $           "nchr =",nchr,"nchz =",nchz
            WRITE(out_unit,'(4(1x,a4,f12.3))')
     $           "a =",cha,"b =",chb,"c =",chc,"d =",chd
            WRITE(out_unit,*)
            WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chear","cheaz",
     $           "real(chear)","imag(chear)","real(cheaz)","imag(cheaz)"
            DO i=0,nchr
               DO j=0,nchz
                  WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $                 REAL(chear(i,j)),AIMAG(chear(i,j)),
     $                 REAL(cheaz(i,j)),AIMAG(cheaz(i,j))
               ENDDO
             ENDDO
            CALL ascii_close(out_unit)

            CALL ascii_open(out_unit,"gpec_xrzphi_chebyshev_n"//
     $        TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_XRZPHI_CHEBYSHEV: "//
     $           "Chebyshev coefficients for displacement"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(2(1x,a12,I4))')
     $           "nchr =",nchr,"nchz =",nchz
            WRITE(out_unit,*)
            WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chxar","chxaz",
     $           "real(chxar)","imag(chxar)","real(chxaz)","imag(chxaz)"
            DO i=0,nchr
               DO j=0,nchz
                  WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $                 REAL(chxar(i,j)),AIMAG(chxar(i,j)),
     $                 REAL(chxaz(i,j)),AIMAG(chxaz(i,j))
               ENDDO
             ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         IF(verbose) WRITE(*,*)"Recontructing xbrzphi by chebyshev"

         DO i=0,nr
            DO j=0,nz
               DO k=0,nchr
                  DO l=0,nchz
                     bcr(i,j)=bcr(i,j)+chear(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcz(i,j)=bcz(i,j)+cheaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     bcp(i,j)=bcp(i,j)+chear(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+cheaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                     xcr(i,j)=xcr(i,j)+chxar(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcz(i,j)=xcz(i,j)+chxaz(k,l)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))
                     xcp(i,j)=xcp(i,j)+chxar(k,l)*
     $                    gdr(i,j)/cha/sqrt(1.0-chex(i)**2.0)*
     $                    sin(REAL(k,r8)*acos(chex(i)))*
     $                    cos(REAL(l,r8)*acos(chey(j)))+chxaz(k,l)*
     $                    gdr(i,j)/chc/sqrt(1.0-chey(i)**2.0)*
     $                    cos(REAL(k,r8)*acos(chex(i)))*
     $                    sin(REAL(l,r8)*acos(chey(j)))
                  ENDDO
               ENDDO
               bcp(i,j)=bcp(i,j)+bcr(i,j)
               xcp(i,j)=xcp(i,j)+xcr(i,j)
            ENDDO
         ENDDO
         bcp=bcp*(-ifac/nn)
         xcp=xcp*(-ifac/nn)

         DEALLOCATE(chex,chey,chear,cheaz,chxar,chxaz)

      ENDIF

      IF (pbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL gpvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF
      ENDIF

      IF (vvbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)
     $      "Computing vacuum fields without plasma response"
         CALL gpvacuum_bnormal(psilim,bnimn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vvbr,vvbz,vvbp)
         IF (helicity<0) THEN
            vvbr=CONJG(vvbr)
            vvbz=CONJG(vvbz)
            vvbp=CONJG(vvbp)
         ENDIF
      ENDIF  
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF(netcdf_flag)THEN
         CALL check( nf90_open(cncfile,nf90_write,cncid) )
         CALL check( nf90_inq_dimid(cncid,"i",i_id) )
         CALL check( nf90_inq_dimid(cncid,"R",r_id) )
         CALL check( nf90_inq_dimid(cncid,"z",z_id) )
         CALL check( nf90_redef(cncid))
         CALL check( nf90_def_var(cncid, "b_r_equil", nf90_double,
     $               (/r_id,z_id/),bre_id) )
         CALL check( nf90_put_att(cncid,bre_id,"long_name",
     $               "Radial equilibrium field") )
         CALL check( nf90_put_att(cncid,bre_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_z_equil", nf90_double,
     $               (/r_id,z_id/),bze_id) )
         CALL check( nf90_put_att(cncid,bze_id,"long_name",
     $               "Vertical equilibrium field") )
         CALL check( nf90_put_att(cncid,bze_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_t_equil", nf90_double,
     $               (/r_id,z_id/),bpe_id) )
         CALL check( nf90_put_att(cncid,bpe_id,"long_name",
     $               "Toroidal equilibrium field") )
         CALL check( nf90_put_att(cncid,bpe_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_r_plasma", nf90_double,
     $               (/r_id,z_id,i_id/),brp_id) )
         CALL check( nf90_put_att(cncid,brp_id,"long_name",
     $               "Radial plasma field") )
         CALL check( nf90_put_att(cncid,brp_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_z_plasma", nf90_double,
     $               (/r_id,z_id,i_id/),bzp_id) )
         CALL check( nf90_put_att(cncid,bzp_id,"long_name",
     $               "Vertical plasma field") )
         CALL check( nf90_put_att(cncid,bzp_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_t_plasma", nf90_double,
     $               (/r_id,z_id,i_id/),bpp_id) )
         CALL check( nf90_put_att(cncid,bpp_id,"long_name",
     $               "Toroidal plasma field") )
         CALL check( nf90_put_att(cncid,bpp_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_r", nf90_double,
     $               (/r_id,z_id,i_id/),br_id) )
         CALL check( nf90_put_att(cncid,br_id,"long_name",
     $               "Radial nonaxisymmetric field") )
         CALL check( nf90_put_att(cncid,br_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_z", nf90_double,
     $               (/r_id,z_id,i_id/),bz_id) )
         CALL check( nf90_put_att(cncid,bz_id,"long_name",
     $               "Vertical nonaxisymmetric field") )
         CALL check( nf90_put_att(cncid,bz_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "b_t", nf90_double,
     $               (/r_id,z_id,i_id/),bp_id) )
         CALL check( nf90_put_att(cncid,bp_id,"long_name",
     $               "Toroidal nonaxisymmetric field") )
         CALL check( nf90_put_att(cncid,bp_id,"units","Tesla") )
         CALL check( nf90_def_var(cncid, "xi_r", nf90_double,
     $               (/r_id,z_id,i_id/),xr_id) )
         CALL check( nf90_put_att(cncid,xr_id,"long_name",
     $               "Radial nonaxisymmetric displacement") )
         CALL check( nf90_put_att(cncid,xr_id,"units","m") )
         CALL check( nf90_def_var(cncid, "xi_z", nf90_double,
     $               (/r_id,z_id,i_id/),xz_id) )
         CALL check( nf90_put_att(cncid,xz_id,"long_name",
     $               "Vertical nonaxisymmetric displacement") )
         CALL check( nf90_put_att(cncid,bz_id,"units","m") )
         CALL check( nf90_def_var(cncid, "xi_t", nf90_double,
     $               (/r_id,z_id,i_id/),xp_id) )
         CALL check( nf90_put_att(cncid,xp_id,"long_name",
     $               "Toroidal nonaxisymmetric displacement") )
         CALL check( nf90_put_att(cncid,bp_id,"units","m") )
         CALL check( nf90_enddef(cncid) )
         CALL check( nf90_put_var(cncid,bre_id,ebr) )
         CALL check( nf90_put_var(cncid,bze_id,ebz) )
         CALL check( nf90_put_var(cncid,bpe_id,ebp) )
         CALL check( nf90_put_var(cncid,br_id,RESHAPE((/REAL(btr),
     $                AIMAG(btr)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bz_id,RESHAPE((/REAL(btz),
     $                AIMAG(btz)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bp_id,RESHAPE((/REAL(btp),
     $                AIMAG(btp)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,brp_id,RESHAPE((/REAL(bpr),
     $                AIMAG(bpr)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bzp_id,RESHAPE((/REAL(bpz),
     $                AIMAG(bpz)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bpp_id,RESHAPE((/REAL(bpp),
     $                AIMAG(bpp)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,xr_id,RESHAPE((/REAL(xrr),
     $                AIMAG(xrr)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,xz_id,RESHAPE((/REAL(xrz),
     $                AIMAG(xrz)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,xp_id,RESHAPE((/REAL(xrp),
     $                AIMAG(xrp)/),(/nr+1,nz+1,2/))) )
         CALL check( nf90_close(cncid) )
      ENDIF

      IF (eqbrzphi_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_eqbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_EQBRZPHI: Eq. b field in rzphi grid"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,5(a17))')"l","r","z","eb_r",
     $        "eb_z","eb_phi"
      
         DO i=0,nr
            DO j=0,nz 
               WRITE(out_unit,'(1x,I2,5(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              ebr(i,j),ebz(i,j),ebp(i,j)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (brzphi_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_pbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PBRZPHI: Perturbed field by plasma"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(bpr(i,j)),AIMAG(bpr(i,j)),
     $              REAL(bpz(i,j)),AIMAG(bpz(i,j)),
     $              REAL(bpp(i,j)),AIMAG(bpp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         CALL ascii_open(out_unit,"gpec_brzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BRZPHI: Total perturbed field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(btr(i,j)),AIMAG(btr(i,j)),
     $              REAL(btz(i,j)),AIMAG(btz(i,j)),
     $              REAL(btp(i,j)),AIMAG(btp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"gpec_chebrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CEHBRZPHI: Total perturbed field"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(bcr(i,j)),AIMAG(bcr(i,j)),
     $                 REAL(bcz(i,j)),AIMAG(bcz(i,j)),
     $                 REAL(bcp(i,j)),AIMAG(bcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         IF (coil_flag) THEN
            CALL ascii_open(out_unit,"gpec_cbrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CBRZPHI: External field by coils"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(vcbr(i,j)),AIMAG(vcbr(i,j)),
     $                 REAL(vcbz(i,j)),AIMAG(vcbz(i,j)),
     $                 REAL(vcbp(i,j)),AIMAG(vcbp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

         IF (vbrzphi_flag) THEN
            CALL ascii_open(out_unit,"gpec_vbrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_VBRZPHI: Total perturbed field "//
     $           "by surface currents"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $           "real(b_phi)","imag(b_phi)"
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(brr(i,j)),AIMAG(brr(i,j)),
     $                 REAL(brz(i,j)),AIMAG(brz(i,j)),
     $                 REAL(brp(i,j)),AIMAG(brp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF

      ENDIF

      IF (xrzphi_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_xrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XRZPHI: Displacement"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(xi_r)","imag(xi_r)","real(xi_z)","imag(xi_z)",
     $        "real(xi_phi)","imag(xi_phi)"
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              gdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(xrr(i,j)),AIMAG(xrr(i,j)),
     $              REAL(xrz(i,j)),AIMAG(xrz(i,j)),
     $              REAL(xrp(i,j)),AIMAG(xrp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         IF (chebyshev_flag) THEN
            CALL ascii_open(out_unit,"gpec_chexrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_CEHXRZPHI: Displacement"
            WRITE(out_unit,*)version
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $           "real(x_r)","imag(x_r)","real(x_z)","imag(x_z)",
     $           "real(x_phi)","imag(x_phi)"
            DO i=0,nr
               DO j=0,nz
                  WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $                 gdl(i,j),gdr(i,j),gdz(i,j),
     $                 REAL(xcr(i,j)),AIMAG(xcr(i,j)),
     $                 REAL(xcz(i,j)),AIMAG(xcz(i,j)),
     $                 REAL(xcp(i,j)),AIMAG(xcp(i,j))
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF
      ENDIF

      IF (pbrzphi_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_vpbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VPBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              vgdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(vpbr(i,j)),AIMAG(vpbr(i,j)),
     $              REAL(vpbz(i,j)),AIMAG(vpbz(i,j)),
     $              REAL(vpbp(i,j)),AIMAG(vpbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (vvbrzphi_flag .AND. ascii_flag) THEN
         CALL ascii_open(out_unit,"gpec_vvbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VVBRZPHI: Vacuum field by "//
     $        "surface currents"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $        "real(b_r)","imag(b_r)","real(b_z)","imag(b_z)",
     $        "real(b_phi)","imag(b_phi)"
         DO i=0,nr
            DO j=0,nz
               WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $              vgdl(i,j),gdr(i,j),gdz(i,j),
     $              REAL(vvbr(i,j)),AIMAG(vvbr(i,j)),
     $              REAL(vvbz(i,j)),AIMAG(vvbz(i,j)),
     $              REAL(vvbp(i,j)),AIMAG(vvbp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag) DEALLOCATE(gdr,gdz,gdl,gdpsi,gdthe,gdphi)
      CALL gpeq_rzpgrid(nr,nz,psixy) ! reset the grid
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_xbrzphi
c-----------------------------------------------------------------------
c     subprogram 13. gpout_vsbrzphi.
c     write brzphi components restored by removing shielding currents.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_vsbrzphi(ss_flag,nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      LOGICAL, DIMENSION(100), INTENT(IN) :: ss_flag
      INTEGER, INTENT(IN) :: nr,nz

      INTEGER :: i,j,k,snum,iss,ns
      INTEGER :: i_id, r_id, z_id, br_id, bz_id, bp_id, s_id, sdid

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(:), ALLOCATABLE :: snums
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: vgdl
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: vgdr,vgdz
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: vbr,vbz,vbp

c-----------------------------------------------------------------------
c     Find out how many surfaces we are calculating
c-----------------------------------------------------------------------
      ns = 0
      DO i=1,msing
         IF (ss_flag(i)) ns = ns + 1
      ENDDO
      IF (ns == 0) RETURN ! only bother is any ss_flags are true
      ALLOCATE(snums(ns), vgdl(ns,0:nr,0:nz), vgdr(ns,0:nr,0:nz),
     $         vgdz(ns,0:nr,0:nz), vbr(ns,0:nr,0:nz), vbz(ns,0:nr,0:nz),
     $         vbp(ns,0:nr,0:nz))
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      vbr = 0
      vbz = 0
      vbp = 0

      iss = 0
      DO snum=1,msing
         IF (ss_flag(snum)) THEN
            iss = iss + 1
            snums(iss) = iss
            IF(verbose) WRITE(*,'(1x,a55,I2)')"Computing vacuum "//
     $           "fields from resonant surface current #",snum
            CALL gpvacuum_bnormal(singtype(snum)%psifac,
     $           singbno_mn(:,snum),nr,nz)
            CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,
     $           complex_flag,nr,nz,vgdl(iss,:,:),vgdr(iss,:,:),
     $           vgdz(iss,:,:),vbr(iss,:,:),vbz(iss,:,:),vbp(iss,:,:))
         ENDIF
      ENDDO
      IF (helicity<0) THEN
         vbr=CONJG(vbr)
         vbz=CONJG(vbz)
         vbp=CONJG(vbp)
      ENDIF
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_vsbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VSBRZPHI: Vacuum field in rzphi grid "//
     $        "by "//TRIM(ss)//"th resonant field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a2,2(1x,a16),a2,6(1x,a16))')"s","r","z",
     $        "l","real(vsb_r)","imag(vsb_r)","real(vsb_z)",
     $        "imag(vsb_z)","real(vsb_phi)","imag(vsb_phi)"
         DO i=1,ns
            DO j=0,nr
               DO k=0,nz
                  WRITE(out_unit,'(1x,I2,2(es17.8e3),I2,6(es17.8e3))')
     $                 snums(i),vgdr(i,j,k),vgdz(i,j,k),vgdl(i,j,k),
     $                 REAL(vbr(i,j,k)),AIMAG(vbr(i,j,k)),
     $                 REAL(vbz(i,j,k)),AIMAG(vbz(i,j,k)),
     $                 REAL(vbp(i,j,k)),AIMAG(vbp(i,j,k))
               ENDDO
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(cncfile,nf90_write,cncid) )
         CALL check( nf90_inq_dimid(cncid,"i",i_id) )
         CALL check( nf90_inq_dimid(cncid,"R",r_id) )
         CALL check( nf90_inq_dimid(cncid,"z",z_id) )

         CALL check( nf90_redef(cncid))
         CALL check( nf90_def_dim(cncid,"s",ns,sdid) )
         CALL check( nf90_def_var(cncid,"s",nf90_int,sdid,s_id) )
         CALL check( nf90_put_att(cncid, s_id ,"long_name",
     $    "Singular surface number from core") )
         CALL check( nf90_def_var(cncid, "b_r_rational"//trim(ss),
     $       nf90_double, (/sdid,r_id,z_id,i_id/),br_id) )
         CALL check( nf90_put_att(cncid,br_id,"units","T") )
         CALL check( nf90_put_att(cncid,br_id,"long_name",
     $      "Radial field from rational surface "//trim(ss)) )
         CALL check( nf90_def_var(cncid, "b_z_rational"//trim(ss),
     $       nf90_double, (/sdid,r_id,z_id,i_id/),bz_id) )
         CALL check( nf90_put_att(cncid,bz_id,"units","T") )
         CALL check( nf90_put_att(cncid,bz_id,"long_name",
     $      "Vertical field from rational surface "//trim(ss)) )
         CALL check( nf90_def_var(cncid, "b_t_rational"//trim(ss),
     $       nf90_double, (/sdid,r_id,z_id,i_id/),bp_id) )
         CALL check( nf90_put_att(cncid,bp_id,"units","T") )
         CALL check( nf90_put_att(cncid,bp_id,"long_name",
     $      "Toroidal field from rational surface "//trim(ss)) )
         CALL check( nf90_enddef(cncid) )
         CALL check( nf90_put_var(cncid,s_id,snums) )
         CALL check( nf90_put_var(cncid,br_id,RESHAPE((/REAL(vbr),
     $                AIMAG(vbr)/),(/ns,nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bz_id,RESHAPE((/REAL(vbz),
     $                AIMAG(vbz)/),(/ns,nr+1,nz+1,2/))) )
         CALL check( nf90_put_var(cncid,bp_id,RESHAPE((/REAL(vbp),
     $                AIMAG(vbp)/),(/ns,nr+1,nz+1,2/))) )
         CALL check( nf90_close(cncid) )
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      DEALLOCATE(snums, vgdl, vgdr, vgdz, vbr, vbz, vbp)
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_vsbrzphi
c-----------------------------------------------------------------------
c     subprogram 14. gpout_xbrzphifun
c-----------------------------------------------------------------------
      SUBROUTINE gpout_xbrzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i_id, p_id, t_id, xr_id,xz_id,xp_id, br_id,bz_id,bp_id
      INTEGER :: istep,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,dphi,t11,t12,t21,t22,t33
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf))
      ALLOCATE(xrr_fun(mstep,0:mthsurf),brr_fun(mstep,0:mthsurf),
     $   xrz_fun(mstep,0:mthsurf),brz_fun(mstep,0:mthsurf),
     $   xrp_fun(mstep,0:mthsurf),brp_fun(mstep,0:mthsurf))
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      IF(verbose) WRITE(*,*)"Computing x and b rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL gpeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a29)')
     $        "volume = ",iindex,"% xi and b rzphi computations"
         CALL gpeq_sol(psifac(istep))
         CALL gpeq_contra(psifac(istep))
         CALL gpeq_cova(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            jacs(itheta)=jac
            v(1,1)=rzphi%fx(1)/(2*rfac)
            v(1,2)=rzphi%fx(2)*twopi*rfac
            v(2,1)=rzphi%fy(1)/(2*rfac)
            v(2,2)=(1+rzphi%fy(2))*twopi*rfac
            v(3,3)=twopi*rs(istep,itheta)
            t11(itheta)=cos(eta)*v(1,1)-sin(eta)*v(1,2)
            t12(itheta)=cos(eta)*v(2,1)-sin(eta)*v(2,2)
            t21(itheta)=sin(eta)*v(1,1)+cos(eta)*v(1,2)
            t22(itheta)=sin(eta)*v(2,1)+cos(eta)*v(2,2)
            t33(itheta)=-1.0/v(3,3)
         ENDDO

         CALL iscdftb(mfac,mpert,xwp_fun,mthsurf,xwp_mn)
         CALL iscdftb(mfac,mpert,bwp_fun,mthsurf,bwp_mn)
         CALL iscdftb(mfac,mpert,xwt_fun,mthsurf,xmt_mn)
         CALL iscdftb(mfac,mpert,bwt_fun,mthsurf,bmt_mn)
         CALL iscdftb(mfac,mpert,xvz_fun,mthsurf,xvz_mn)
         CALL iscdftb(mfac,mpert,bvz_fun,mthsurf,bvz_mn)

         xrr_fun(istep,:)=(t11*xwp_fun+t12*xwt_fun)/jacs
         brr_fun(istep,:)=(t11*bwp_fun+t12*bwt_fun)/jacs
         xrz_fun(istep,:)=(t21*xwp_fun+t22*xwt_fun)/jacs
         brz_fun(istep,:)=(t21*bwp_fun+t22*bwt_fun)/jacs
         xrp_fun(istep,:)=t33*xvz_fun
         brp_fun(istep,:)=t33*bvz_fun

         xrr_fun(istep,:)=xrr_fun(istep,:)*EXP(ifac*nn*dphi)
         brr_fun(istep,:)=brr_fun(istep,:)*EXP(ifac*nn*dphi)
         xrz_fun(istep,:)=xrz_fun(istep,:)*EXP(ifac*nn*dphi)
         brz_fun(istep,:)=brz_fun(istep,:)*EXP(ifac*nn*dphi)
         xrp_fun(istep,:)=xrp_fun(istep,:)*EXP(ifac*nn*dphi)
         brp_fun(istep,:)=brp_fun(istep,:)*EXP(ifac*nn*dphi)         

      ENDDO

      IF(helicity>0)THEN
         xrr_fun=CONJG(xrr_fun)
         xrz_fun=CONJG(xrz_fun)
         xrp_fun=CONJG(xrp_fun)
         brr_fun=CONJG(brr_fun)
         brz_fun=CONJG(brz_fun)
         brp_fun=CONJG(brp_fun)
      ENDIF

      CALL gpeq_dealloc

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_xbrzphi_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBRZPHI_FUN: "//
     $        "Rzphi components of displacement and field in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(14(1x,a16))')"r","z",
     $        "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $        "real(xp)","imag(xp)","real(br)","imag(br)",
     $        "real(bz)","imag(bz)","real(bp)","imag(bp)"
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            DO itheta=0,mthsurf
               WRITE(out_unit,'(14(es17.8e3))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(xrr_fun(istep,itheta)),
     $              AIMAG(xrr_fun(istep,itheta)),
     $              REAL(xrz_fun(istep,itheta)),
     $              AIMAG(xrz_fun(istep,itheta)),
     $              REAL(xrp_fun(istep,itheta)),
     $              AIMAG(xrp_fun(istep,itheta)),
     $              REAL(brr_fun(istep,itheta)),
     $              AIMAG(brr_fun(istep,itheta)),
     $              REAL(brz_fun(istep,itheta)),
     $              AIMAG(brz_fun(istep,itheta)),
     $              REAL(brp_fun(istep,itheta)),
     $              AIMAG(brp_fun(istep,itheta))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         IF(debug_flag) PRINT *,"Opening "//TRIM(fncfile)
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         IF(debug_flag) PRINT *,"  Inquiring about dimensions"
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"theta_dcon",t_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         IF(debug_flag) PRINT *,"  Defining variables"
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "xi_r", nf90_double,
     $               (/p_id,t_id,i_id/),xr_id) )
         CALL check( nf90_put_att(fncid,xr_id,"units","m") )
         CALL check( nf90_put_att(fncid,xr_id,"long_name",
     $               "Major radial displacement") )
         CALL check( nf90_def_var(fncid, "xi_z", nf90_double,
     $               (/p_id,t_id,i_id/),xz_id) )
         CALL check( nf90_put_att(fncid,xz_id,"units","m") )
         CALL check( nf90_put_att(fncid,xz_id,"long_name",
     $               "Vertical displacement") )
         CALL check( nf90_def_var(fncid, "xi_phi", nf90_double,
     $               (/p_id,t_id,i_id/),xp_id) )
         CALL check( nf90_put_att(fncid,xp_id,"units","m") )
         CALL check( nf90_put_att(fncid,xp_id,"long_name",
     $               "Toroidal displacement") )
         CALL check( nf90_def_var(fncid, "b_r", nf90_double,
     $               (/p_id,t_id,i_id/),br_id) )
         CALL check( nf90_put_att(fncid,br_id,"units","T") )
         CALL check( nf90_put_att(fncid,br_id,"long_name",
     $               "Major radial field") )
         CALL check( nf90_def_var(fncid, "b_z", nf90_double,
     $               (/p_id,t_id,i_id/),bz_id) )
         CALL check( nf90_put_att(fncid,bz_id,"units","T") )
         CALL check( nf90_put_att(fncid,bz_id,"long_name",
     $               "Vertical field") )
         CALL check( nf90_def_var(fncid, "b_phi", nf90_double,
     $               (/p_id,t_id,i_id/),bp_id) )
         CALL check( nf90_put_att(fncid,bp_id,"units","T") )
         CALL check( nf90_put_att(fncid,bp_id,"long_name",
     $               "Toroidal field") )
         CALL check( nf90_enddef(fncid) )
         IF(debug_flag) PRINT *,"  Writing variables"
         CALL check( nf90_put_var(fncid,xr_id,RESHAPE((/REAL(xrr_fun),
     $                AIMAG(xrr_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,xz_id,RESHAPE((/REAL(xrz_fun),
     $                AIMAG(xrz_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,xp_id,RESHAPE((/REAL(xrp_fun),
     $                AIMAG(xrp_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,br_id,RESHAPE((/REAL(brr_fun),
     $                AIMAG(brr_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,bz_id,RESHAPE((/REAL(brz_fun),
     $                AIMAG(brz_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,bp_id,RESHAPE((/REAL(brp_fun),
     $                AIMAG(brp_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_close(fncid) )
         IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)
      ENDIF
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(rs,zs)
      DEALLOCATE(xrr_fun,brr_fun,xrz_fun,brz_fun,xrp_fun,brp_fun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_xbrzphifun
c-----------------------------------------------------------------------
c     subprogram 15. gpout_arzphifun
c-----------------------------------------------------------------------
      SUBROUTINE gpout_arzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i_id,p_id,t_id, er_id,ez_id,ep_id, ar_id,az_id,ap_id

      INTEGER :: istep,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xss_fun,
     $     ear,eat,eap,arr,art,arp
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: ear_fun,eaz_fun,
     $     eap_fun,arr_fun,arz_fun,arp_fun
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf))
      ALLOCATE(ear_fun(mstep,0:mthsurf),eaz_fun(mstep,0:mthsurf),
     $   eap_fun(mstep,0:mthsurf),arr_fun(mstep,0:mthsurf),
     $   arz_fun(mstep,0:mthsurf),arp_fun(mstep,0:mthsurf))
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vector potential rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL gpeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a37)')
     $        "volume = ",iindex,"% vector potential rzphi computations"
         CALL gpeq_sol(psifac(istep))
c-----------------------------------------------------------------------
c     normal and two tangent components to flux surface.
c-----------------------------------------------------------------------
         CALL iscdftb(mfac,mpert,xsp_fun,mthsurf,xsp_mn)
         CALL iscdftb(mfac,mpert,xss_fun,mthsurf,xss_mn)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psifac(istep),theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta(itheta)+rzphi%f(2))
            rs(istep,itheta)=ro+rfac*COS(eta)
            zs(istep,itheta)=zo+rfac*SIN(eta)
            dphi(itheta)=rzphi%f(3)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*rs(istep,itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(2,1)=-rzphi%fx(2)*twopi**2*rs(istep,itheta)*rfac/jac
            w(2,2)=rzphi%fx(1)*pi*rs(istep,itheta)/(rfac*jac)
            w(3,1)=rzphi%fx(3)*2*rfac
            w(3,2)=rzphi%fy(3)/(twopi*rfac)
            w(3,3)=1.0/(twopi*rs(istep,itheta))
            ear(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,1)-w(3,1))
            eat(itheta)=chi1*psifac(istep)*(qfac(istep)*w(2,2)-w(3,2))
            eap(itheta)=-chi1*psifac(istep)*w(3,3)
            ear_fun(istep,itheta)=ear(itheta)*cos(eta)-
     $           eat(itheta)*sin(eta)
            eaz_fun(istep,itheta)=ear(itheta)*sin(eta)+
     $           eat(itheta)*cos(eta)
            eap_fun(istep,itheta)=-eap(itheta)            
            arr(itheta)=xss_fun(itheta)*w(1,1)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,1)+xsp_fun(itheta)*w(3,1))
            art(itheta)=xss_fun(itheta)*w(1,2)-chi1*(qfac(istep)*
     $           xsp_fun(itheta)*w(2,2)+xsp_fun(itheta)*w(3,2))
            arp(itheta)=chi1*xsp_fun(itheta)*w(3,3)
            arr_fun(istep,itheta)=arr(itheta)*cos(eta)-
     $           art(itheta)*sin(eta)
            arz_fun(istep,itheta)=arr(itheta)*sin(eta)+
     $           art(itheta)*cos(eta)
            arp_fun(istep,itheta)=-arp(itheta)
         ENDDO
         ear_fun(istep,:)=ear_fun(istep,:)*EXP(ifac*nn*dphi)
         eaz_fun(istep,:)=eaz_fun(istep,:)*EXP(ifac*nn*dphi)
         eap_fun(istep,:)=eap_fun(istep,:)*EXP(ifac*nn*dphi)
         arr_fun(istep,:)=arr_fun(istep,:)*EXP(ifac*nn*dphi)
         arz_fun(istep,:)=arz_fun(istep,:)*EXP(ifac*nn*dphi)
         arp_fun(istep,:)=arp_fun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      IF(helicity>0)THEN
         ear_fun=CONJG(ear_fun)
         eaz_fun=CONJG(eaz_fun)
         eap_fun=CONJG(eap_fun)
         arr_fun=CONJG(arr_fun)
         arz_fun=CONJG(arz_fun)
         arp_fun=CONJG(arp_fun)
      ENDIF

      CALL gpeq_dealloc

      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_arzphi_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_ARZPHI_FUN: "//
     $        "Rzphi components of vector potential in functions"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,*)
         WRITE(out_unit,'(14(1x,a16))')"r","z",
     $        "real(ear)","imag(ear)","real(eaz)","imag(eaz)",
     $        "real(eap)","imag(eap)","real(arr)","imag(arr)",
     $        "real(arz)","imag(arz)","real(arp)","imag(arp)"
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            DO itheta=0,mthsurf
               WRITE(out_unit,'(14(es17.8e3))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(ear_fun(istep,itheta)),
     $              AIMAG(ear_fun(istep,itheta)),
     $              REAL(eaz_fun(istep,itheta)),
     $              AIMAG(eaz_fun(istep,itheta)),
     $              REAL(eap_fun(istep,itheta)),
     $              AIMAG(eap_fun(istep,itheta)),
     $              REAL(arr_fun(istep,itheta)),
     $              AIMAG(arr_fun(istep,itheta)),
     $              REAL(arz_fun(istep,itheta)),
     $              AIMAG(arz_fun(istep,itheta)),
     $              REAL(arp_fun(istep,itheta)),
     $              AIMAG(arp_fun(istep,itheta))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"theta_dcon",t_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "e_r", nf90_double,
     $               (/p_id,t_id,i_id/),er_id) )
         CALL check( nf90_put_att(fncid,er_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,er_id,"long_name",
     $               "Major radial vector potential") )
         CALL check( nf90_def_var(fncid, "e_z", nf90_double,
     $               (/p_id,t_id,i_id/),ez_id) )
         CALL check( nf90_put_att(fncid,ez_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,ez_id,"long_name",
     $               "Vertical vector potential") )
         CALL check( nf90_def_var(fncid, "e_phi", nf90_double,
     $               (/p_id,t_id,i_id/),ep_id) )
         CALL check( nf90_put_att(fncid,ep_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,ep_id,"long_name",
     $               "Toroidal vector potential") )
         CALL check( nf90_def_var(fncid, "a_r", nf90_double,
     $               (/p_id,t_id,i_id/),ar_id) )
         CALL check( nf90_put_att(fncid,ar_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,ar_id,"long_name",
     $               "Major radial vector potential") )
         CALL check( nf90_def_var(fncid, "a_z", nf90_double,
     $               (/p_id,t_id,i_id/),az_id) )
         CALL check( nf90_put_att(fncid,az_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,az_id,"long_name",
     $               "Vertical vector potential") )
         CALL check( nf90_def_var(fncid, "a_phi", nf90_double,
     $               (/p_id,t_id,i_id/),ap_id) )
         CALL check( nf90_put_att(fncid,ap_id,"units","Vs/m") )
         CALL check( nf90_put_att(fncid,ap_id,"long_name",
     $               "Toroidal vector potential") )
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid,er_id,RESHAPE((/REAL(ear_fun),
     $                AIMAG(ear_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,ez_id,RESHAPE((/REAL(eaz_fun),
     $                AIMAG(eaz_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,ep_id,RESHAPE((/REAL(eap_fun),
     $                AIMAG(eap_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,ar_id,RESHAPE((/REAL(arr_fun),
     $                AIMAG(arr_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,az_id,RESHAPE((/REAL(arz_fun),
     $                AIMAG(arz_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(fncid,ap_id,RESHAPE((/REAL(arp_fun),
     $                AIMAG(arp_fun)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_close(fncid) )
      ENDIF
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(rs,zs)
      DEALLOCATE(ear_fun,eaz_fun,eap_fun,arr_fun,arz_fun,arp_fun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_arzphifun
c-----------------------------------------------------------------------
c     subprogram 16. gpout_xclebsch.
c     Write Clebsch coordinate displacements.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_xclebsch(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,itheta,iindex,ids(3)
      INTEGER :: i_id,m_id,p_id,dp_id,xp_id,xa_id
      REAL(r8) :: ileft, psi, rfac, eta, rs, zs

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xmp1out,xspout,xmsout
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xmp1funs,xspfuns,
     $    xmsfuns
c-----------------------------------------------------------------------
c     allocation puts memory in heap, avoiding stack overfill
c-----------------------------------------------------------------------
      ALLOCATE(xmp1out(mstep,lmpert),xspout(mstep,lmpert),
     $   xmsout(mstep,lmpert))
      ALLOCATE(xmp1funs(mstep,0:mthsurf),xspfuns(mstep,0:mthsurf),
     $   xmsfuns(mstep,0:mthsurf))
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing Clebsch displacements"

      CALL idcon_build(egnum,xspmn)
      CALL gpeq_alloc

      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% Clebsch decomposition"
         ! compute contravarient displacement on surface with regulation
         psi = psifac(istep)
         CALL spline_eval(sq,psi,1)
         CALL gpeq_sol(psi)

         ! convert to jac_out
         CALL gpeq_bcoordsout(xmp1out(istep,:),xmp1_mn,psi)
         CALL gpeq_bcoordsout(xspout(istep,:),xsp_mn,psi)
         CALL gpeq_bcoordsout(xmsout(istep,:),xms_mn,psi)

         ! convert to real space
         IF (fun_flag) THEN
            CALL iscdftb(mfac,mpert,xmp1funs(istep,:),mthsurf,xmp1_mn)
            CALL iscdftb(mfac,mpert,xspfuns(istep,:),mthsurf,xsp_mn)
            CALL iscdftb(mfac,mpert,xmsfuns(istep,:),mthsurf,xms_mn)
         ENDIF

      ENDDO

      ! remove DCON normalization
      xmsout = xmsout/chi1
      IF(fun_flag) xmsfuns = xmsfuns/chi1


      IF(ascii_flag)THEN
         CALL ascii_open(out_unit,"gpec_xclebsch_n"//
     $      TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XCLEBSCH: "//
     $      "Clebsch Components of the displacement."
         WRITE(out_unit,*)version
         WRITE(out_unit,'(1/,1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1/,1x,a12,1x,I6,1x,2(a12,I4),1/)')
     $      "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,'(1x,a23,1x,a4,6(1x,a16))')"psi","m",
     $      "real(derxi^psi)","imag(derxi^psi)",
     $      "real(xi^psi)","imag(xi^psi)",
     $      "real(xi^alpha)","imag(xi^alpha)"
         DO istep=1,mstep,MAX(1,(mstep*(lmpert)-1)/max_linesout+1)
            DO ipert=1,lmpert
               WRITE(out_unit,'(1x,es23.15,1x,I4,6(es17.8e3))')
     $            psifac(istep),lmfac(ipert),xmp1out(istep,ipert),
     $            xspout(istep,ipert),xmsout(istep,ipert)
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      IF(fun_flag .AND. ascii_flag)THEN
          CALL ascii_open(out_unit,"gpec_xclebsch_fun_n"//
     $       TRIM(sn)//".out","UNKNOWN")
          WRITE(out_unit,*)"GPEC_XCLEBSCH: "//
     $       "Clebsch Components of the displacement."
          WRITE(out_unit,*)version
          WRITE(out_unit,'(1/,1x,a13,a8)')"jac_out = ",jac_out
          WRITE(out_unit,'(1/,1x,a12,1x,I6,1x,2(a12,I4),1/)')
     $       "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
          WRITE(out_unit,'(1x,a23,9(1x,a16))')"psi","theta","r","z",
     $       "real(derxi^psi)","imag(derxi^psi)",
     $       "real(xi^psi)","imag(xi^psi)",
     $       "real(xi^alpha)","imag(xi^alpha)"
          DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
             DO itheta=0,mthsurf
                CALL bicube_eval(rzphi,psifac(istep),theta(itheta),0)
                rfac = SQRT(rzphi%f(1))
                eta = twopi * (theta(itheta) + rzphi%f(2))
                rs = ro + rfac * COS(eta)
                zs = zo + rfac * SIN(eta)
                WRITE(out_unit,'(1x,es23.15,9(es17.8e3))')psifac(istep),
     $             theta(itheta),rs,zs,xmp1funs(istep,itheta),
     $             xspfuns(istep,itheta),xmsfuns(istep,itheta)
             ENDDO
          ENDDO
          CALL ascii_close(out_unit)
      ENDIF

      IF(netcdf_flag)THEN
         CALL check( nf90_open(fncfile,nf90_write,fncid) )
         CALL check( nf90_inq_dimid(fncid,"i",i_id) )
         CALL check( nf90_inq_dimid(fncid,"m_out",m_id) )
         CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
         CALL check( nf90_redef(fncid))
         CALL check( nf90_def_var(fncid, "xigradpsi_dpsi", nf90_double,
     $               (/p_id,m_id,i_id/),dp_id) )
         CALL check( nf90_put_att(fncid,dp_id,"long_name",
     $      "Psi derivative of contravariant psi displacement") )
         CALL check( nf90_def_var(fncid, "xigradpsi", nf90_double,
     $               (/p_id,m_id,i_id/),xp_id) )
         CALL check( nf90_put_att(fncid,xp_id,"long_name",
     $      "Contravariant psi displacement") )
         CALL check( nf90_def_var(fncid, "xigradalpha", nf90_double,
     $               (/p_id,m_id,i_id/),xa_id) )
         CALL check( nf90_put_att(fncid,xa_id,"long_name",
     $      "Contravariant Clebsch angle displacement") )
         ids = (/xp_id,dp_id,xa_id/)
         DO i=1,3
            CALL check( nf90_put_att(fncid,ids(i),"units","m") )
            CALL check( nf90_put_att(fncid,ids(i),"jacobian",jac_out) )
         ENDDO
         CALL check( nf90_enddef(fncid) )
         CALL check( nf90_put_var(fncid,dp_id,RESHAPE((/REAL(xmp1out),
     $                AIMAG(xmp1out)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,xp_id,RESHAPE((/REAL(xspout),
     $                AIMAG(xspout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(fncid,xa_id,RESHAPE((/REAL(xmsout),
     $                AIMAG(xmsout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_close(fncid) )
      ENDIF

      CALL gpeq_dealloc
c-----------------------------------------------------------------------
c     deallocation cleans memory in heap
c-----------------------------------------------------------------------
      DEALLOCATE(xmp1out,xspout,xmsout)
      DEALLOCATE(xmp1funs,xspfuns,xmsfuns)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_xclebsch
c-----------------------------------------------------------------------
c     subprogram 17. gpout_control_filter.
c     Filter control surface flux vector in flux bases with energy norms
c-----------------------------------------------------------------------
      SUBROUTINE gpout_control_filter(finmn,foutmn,ftypes,fmodes,
     $           rout,bpout,bout,rcout,tout,op_write)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: fmodes,rout,bpout,bout,rcout,tout
      CHARACTER(len=*), INTENT(IN) :: ftypes
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn,foutmn
      LOGICAL, INTENT(IN), OPTIONAL :: op_write

      ! eigendecompositions variables
      INTEGER :: info,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      ! SVD variables
      REAL(r8), DIMENSION(5*mpert) :: rworksvd
      REAL(r8), DIMENSION(5*msing) :: sworksvd
      COMPLEX(r8), DIMENSION(3*mpert) :: worksvd

      LOGICAL :: output
      INTEGER :: i,j,k,maxmode
      INTEGER :: idid,mdid,xdid,wdid,rdid,pdid,sdid,tdid,cdid,
     $   mx_id,mw_id,mr_id,mp_id,mc_id,
     $   we_id,re_id,pe_id,se_id,fc_id,fcf_id,sl_id,cn_id,
     $   w_id,r_id,p_id,s_id, sc_id,wr_id,wp_id,rp_id,ws_id,rs_id,ps_id,
     $   ft_id,fx_id,wx_id,rx_id,px_id,sx_id,wa_id,rl_id,
     $   x_id,xe_id,xt_id,wf_id,rf_id,sf_id,ex_id,et_id,
     $   wev_id,wes_id,wep_id,rev_id,res_id,rep_id,sev_id,ses_id,sep_id,
     $   etf_id,ftf_id,exf_id,fxf_id,rm_id,wm_id,pm_id
      REAL(r8) :: norm
      REAL(r8), DIMENSION(0:mthsurf) :: units
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: tempfun
      COMPLEX(r8), DIMENSION(0:mthsurf, coil_num) :: tempfunc
      COMPLEX(r8), DIMENSION(msing) :: temps
      COMPLEX(r8), DIMENSION(mpert) :: temp,tempm,eigmn,filmn
      COMPLEX(r8), DIMENSION(mpert,mpert) :: mat,matmm,tempmm,singmat,
     $    sqrtamat,ftop,ptof
      COMPLEX(r8), DIMENSION(mpert,msing) :: matms
      COMPLEX(r8), DIMENSION(msing,msing) :: matss
      COMPLEX(r8), DIMENSION(msing,mpert) :: matsm,singcoupmat
      CHARACTER(64) :: message

      INTEGER,  DIMENSION(mpert) :: aindx,indx
      REAL(r8), DIMENSION(mpert) :: xvals,wvals,rvals,pvals,avals,
     $    rlvals,sengys,vengys,pengys
      REAL(r8), DIMENSION(msing) :: svals
      COMPLEX(r8), DIMENSION(mpert,mpert)::xvecs,wvecs,rvecs,pvecs,
     $    avecs,wmat,rmat,pmat
      COMPLEX(r8), DIMENSION(mpert,msing) :: svecs
      COMPLEX(r8), DIMENSION(0:mthsurf,mpert)::wfuns,rfuns
      COMPLEX(r8), DIMENSION(0:mthsurf,msing)::sfuns

      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing Energy-Normalized flux bases"
c-----------------------------------------------------------------------
c     basic definitions
c-----------------------------------------------------------------------
      IF(.NOT. PRESENT(op_write)) THEN
         output = .FALSE.
      ELSE
        output = op_write
      ENDIF
c-----------------------------------------------------------------------
c     calculate singular-factor matrix
c-----------------------------------------------------------------------
      singmat = 0
      DO i=1,mpert
         singmat(i,i) = chi1*twopi*ifac*(mfac(i)-nn*qlim)
      ENDDO
c-----------------------------------------------------------------------
c     calculate the scalar surface area for normalizations
c-----------------------------------------------------------------------
      units = 1
      jarea = issurfint(units,mthsurf,psilim,0,0)
c-----------------------------------------------------------------------
c     calculate area dimensioned half-area weighting matrix W as,
c     W_m,m' = int{sqrt(J|delpsi|)exp[-i*(m-m')t]dt} * int{sqrt(J|delpsi|)dt}
c-----------------------------------------------------------------------
      sqrtamat = 0
      DO i=1,mpert
         sqrtamat(i,i) = 1.0
         CALL gpeq_weight(psilim,sqrtamat(:,i),mfac,mpert,2) ! A^1/2
      ENDDO
      ptof = sqrtamat * sqrt(jarea) ! transform power-norm field to flux
      CALL iszhinv(ptof,mpert,ftop)
c-----------------------------------------------------------------------
c     compute DCON energy per displacement eigenvalues and eigenvectors.
c-----------------------------------------------------------------------
      ! start with the inverse of the inductance matrix (dW = 0.5 Phi Lambda^-1 Phi)
      ! remove border of modes/solutions (diagnostic only)
      i = malias+1
      j = mpert-malias
      xvecs = 0
      xvecs(i:j,i:j) = 0.5*plas_indinvmats(resp_index,i:j,i:j)
      ! convert to displacement
      xvecs = MATMUL(MATMUL(singmat,xvecs),CONJG(singmat))
      ! get eigenvalues and eigenvectors
      work = 0
      rwork = 0
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,xvecs,mpert,xvals,work,lwork,rwork,info)
c-----------------------------------------------------------------------
c     compute energy eigenvalues and eigenvectors.
c      - We want to use Phi'=bsqrtA/|sqrtA| basis so a eigenvector norm means
c        int{b^2da}/int{da} = 1
c     - we start with E = xi* Wt xi where * is the conjugate transpose
c     - convert to total flux so E = Phi* F Phi
c     - convert to external flux using permeability E = Phi* P* FP Phi
c     - then use our weighting matrix to get E = Phi'* W*P*FPW Phi'
c-----------------------------------------------------------------------
      ! start with the inverse of the inductance matrix (dW = 0.5 Phi Lambda^-1 Phi)
      ! remove border of modes/solutions (diagnostic only)
      i = malias+1
      j = mpert-malias
      wvecs = 0
      wvecs(i:j,i:j) = 0.5*plas_indinvmats(resp_index,i:j,i:j)
      ! convert to external flux
      mat = permeabmats(resp_index,:,:)
      wvecs=MATMUL(MATMUL(CONJG(TRANSPOSE(mat)),wvecs),mat)
      ! convert to bsqrtA/|sqrtA|
      wmat = MATMUL(MATMUL(ptof,wvecs),ptof)*2*mu0
      wvecs = wmat
      ! get eigenvalues and eigenvectors
      work = 0
      rwork = 0
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,wvecs,mpert,wvals,work,lwork,rwork,info)
      ! put in descending order like zgesvd
      wvals(:)   = wvals(mpert:1:-1)
      wvecs(:,:) = wvecs(:,mpert:1:-1)
c-----------------------------------------------------------------------
c     re-order energy eigenmodes by amplification dW_vac/dW.
c-----------------------------------------------------------------------
      mat=MATMUL(MATMUL(ptof,surf_indinvmats),ptof)
      DO i=1,mpert
         tempm = wvecs(:,i)
         avals(i) = 0.5*REAL(DOT_PRODUCT(tempm,MATMUL(mat,tempm)))
         avals(i) = avals(i)/wvals(i)
      ENDDO
      aindx = (/(i,i=1,mpert)/)
      CALL isbubble(avals,aindx,1,mpert)
      DO i=1,mpert
         avecs(:,i) = wvecs(:,aindx(i))
      ENDDO
c-----------------------------------------------------------------------
c     Calculate reluctance eigenvectors and eigenvalues
c      - We want to use Phi'=BsqrtA/|sqrtA| basis so a eigenvector norm means
c        int{b^2da}/int{da} = 1
c        - This is a physically meaningful quantity (energy), and thus
c          independent of coordinate system
c      - Normalizing the input, I = RPhi becomes I = RWPhi'
c      - Normalizing the output gives W^dagger I = W^daggerRW Phi'
c      - Since W is Hermitian (W^dagger=W) the new operator is too
c      - Eigenvalues correspond to int{I^2da} (power) per int{b^2da} (energy)
c-----------------------------------------------------------------------
      ! start with GPEC flux matrix
      ! remove border of modes/solutions (diagnostic only)
      i = malias+1
      j = mpert-malias
      rvecs = 0
      rvecs(i:j,i:j) = reluctmats(resp_index,i:j,i:j)
      ! convert to bsqrtA/|sqrtA|
      rmat = MATMUL(MATMUL(ptof,rvecs),ptof)
      rvecs = rmat
      work = 0
      rwork = 0
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,rvecs,mpert,rvals,work,lwork,rwork,info)
      ! put in descending order like zgesvd
      rvals(:)   = rvals(mpert:1:-1)
      rvecs(:,:) = rvecs(:,mpert:1:-1)
c-----------------------------------------------------------------------
c     Calculate permeability right-singular vectors and singular values
c      - We want to use Phi'=bsqrtA/|sqrtA| basis so a eigenvector norm means
c        int{b^2da}/int{da} = 1
c        - This is a physically meaningful quantity (energy)
c      - Phi = P Phix -> WPhi' = PW Phix' -> Phi' = W*PW Phix'
c      - Eigenvalues correspond to int{Phi^2da} (energy) per int{Phix'^2da} (energy)
c-----------------------------------------------------------------------
      ! remove border of modes/solutions (diagnostic only)
      i = malias+1
      j = mpert-malias
      mat = 0
      mat(i:j,i:j) = permeabmats(resp_index,i:j,i:j)
      ! convert to bsqrtA/|sqrtA|
      pmat = MATMUL(MATMUL(ptof,mat),ptof)
      mat = TRANSPOSE(pmat)
      worksvd=0
      rworksvd=0
      lwork=3*mpert
      CALL zgesvd('S','S',mpert,mpert,mat,mpert,pvals,matmm,mpert,
     $     pvecs,mpert,worksvd,lwork,rworksvd,info)
      pvecs=CONJG(TRANSPOSE(pvecs))
c-----------------------------------------------------------------------
c     Convert singcoup vectors to working coordinates.
c-----------------------------------------------------------------------
      IF (singcoup_set) THEN
         matsm = 0
         ! field to sqrt flux transform matrix in DCON coordinate
         matmm=0
         DO i=1,mpert
            temp=0
            temp(i)=1.0
            CALL gpeq_weight(psilim,temp,mfac,mpert,2)
            matmm(:,i)=temp/sqrt(jarea)
         ENDDO
         ! inverse
         tempmm=0
         DO i=1,mpert
            tempmm(i,i)=1
         ENDDO
         CALL zgetrf(mpert,mpert,matmm,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,matmm,mpert,ipiv,tempmm,mpert,info)
         ! singular coupling in DCON coordinate
         matsm=MATMUL(singcoup(1,:,:),tempmm)
         singcoupmat = matsm
         lwork = 3*mpert
         worksvd  = 0
         sworksvd = 0
         CALL zgesvd('S','O',msing,mpert,matsm,msing,svals, !'O' writes VT to A
     $        matss,msing,matsm,msing,worksvd,lwork,sworksvd,info)
         svecs=CONJG(TRANSPOSE(matsm))
      ENDIF
c-----------------------------------------------------------------------
c     Filter to keep desired physics modes
c-----------------------------------------------------------------------
      foutmn = MATMUL(permeabmats(resp_index,:,:),finmn) ! total flux
      IF(fmodes/=0)THEN
         DO k=1,LEN_TRIM(ftypes)
            filmn = 0
            temp = MATMUL(ftop, finmn) ! flux (Tm^2) to power-normalized field (T)
            SELECT CASE(ftypes(k:k))
            CASE('x')
               temp = foutmn/(chi1*twopi*ifac*(mfac-nn*qlim)) ! total displacement
               mat = xvecs
               message = "DCON eigenmodes"
               maxmode = mpert
            CASE('w')
               mat = wvecs
               message = "energy eigenmodes"
               maxmode = mpert
            CASE('a')
               mat = avecs
               message = "amplification ordered energy eigenmodes"
               maxmode = mpert
            CASE('r')
               mat = rvecs
               message = "reluctance eigenmodes"
               maxmode = mpert
           CASE('p')
               mat = pvecs
               message = "permeability RS vectors"
               maxmode = mpert
            CASE('s')
               mat = svecs
               message = "singular coupling RS vectors"
               maxmode = msing
            CASE DEFAULT
               CYCLE
            END SELECT
            IF (fmodes>0) THEN
               PRINT *,'Isolating perturbation || to largest ',fmodes,
     $            TRIM(message)
               DO i=1,MIN(maxmode,ABS(fmodes))
                  eigmn = mat(:,i)
                  norm=SQRT(ABS(DOT_PRODUCT(eigmn,eigmn)))
                  filmn=filmn+eigmn*DOT_PRODUCT(eigmn,temp)/(norm*norm)
               ENDDO
            ELSE !(fmodes<0) THEN
               PRINT *,'Isolating perturbation || to smallest ',fmodes,
     $            TRIM(message)
               DO j=1,MIN(maxmode,ABS(fmodes))
                  i = maxmode+1-j
                  eigmn = mat(:,i)
                  norm=SQRT(ABS(DOT_PRODUCT(eigmn,eigmn)))
                  filmn=filmn+eigmn*DOT_PRODUCT(eigmn,temp)/(norm*norm)
               ENDDO
            ENDIF
            IF(ftypes(k:k)=='x')THEN
               temp = filmn
               filmn = filmn*(chi1*twopi*ifac*(mfac-nn*qlim)) ! total flux
               filmn = MATMUL(permeabinvmats(resp_index,:,:),filmn) ! external flux
            ELSE
               filmn = MATMUL(ptof, filmn) ! power-normalized field back to true flux
            ENDIF
            finmn = filmn
         ENDDO
      ENDIF
      foutmn = MATMUL(permeabmats(resp_index,:,:),finmn) ! total flux
c-----------------------------------------------------------------------
c     Write outputs
c-----------------------------------------------------------------------
      IF(output .AND. netcdf_flag)THEN
         ! append netcdf file
         CALL check( nf90_open(mncfile,nf90_write,mncid) )

         ! references to pre-defined dimensions
         CALL check( nf90_inq_dimid(mncid,"i",idid) )
         CALL check( nf90_inq_dimid(mncid,"m",mdid) )
         CALL check( nf90_inq_dimid(mncid,"theta",tdid) )
         IF (singcoup_set) THEN
            CALL check( nf90_inq_dimid(mncid,"mode_C",sdid) )
         ENDIF

         ! Start definitions
         CALL check( nf90_redef(mncid))

         IF(debug_flag) PRINT *,"  Defining vecs"
         CALL check( nf90_def_dim(mncid,"mode_X",mpert,   xdid) )
         CALL check( nf90_def_var(mncid,"mode_X",nf90_int,xdid,mx_id))
         CALL check( nf90_put_att(mncid, mx_id ,"long_name",
     $    "Total displacement energy eigenmode index") )
         CALL check( nf90_def_var(mncid,"X_eigenvector",
     $               nf90_double,(/mdid,xdid,idid/),x_id) )
         CALL check( nf90_put_att(mncid,x_id,"long_name",
     $    "Total displacement energy eigendecomposition") )
         CALL check( nf90_def_var(mncid,"X_eigenvalue",
     $               nf90_double,(/xdid/),xe_id) )
         CALL check( nf90_put_att(mncid,xe_id,"units","J/m^2") )
         CALL check( nf90_put_att(mncid,xe_id,"long_name",
     $    "Total displacement energy eigenvalues") )

         CALL check( nf90_def_dim(mncid,"mode_W",mpert,   wdid) )
         CALL check( nf90_def_var(mncid,"mode_W",nf90_int,wdid,mw_id))
         CALL check( nf90_put_att(mncid, mw_id ,"long_name",
     $    "Energy-norm external flux energy eigenmode index") )
         CALL check( nf90_def_var(mncid,"W_xe",
     $               nf90_double,(/mdid,wdid,idid/),wm_id) )
         CALL check( nf90_put_att(mncid,wm_id,"long_name",
     $    "Energy-norm external flux energy matrix") )
         CALL check( nf90_def_var(mncid,"W_xe_eigenvector",
     $               nf90_double,(/mdid,wdid,idid/),w_id) )
         CALL check( nf90_put_att(mncid,w_id,"long_name",
     $    "Energy-norm external flux energy eigendecomposition") )
         CALL check( nf90_def_var(mncid,"W_xe_eigenvalue",
     $               nf90_double,(/wdid/),we_id) )
         CALL check( nf90_put_att(mncid,we_id,"units","J/T^2") )
         CALL check( nf90_put_att(mncid,we_id,"long_name",
     $    "Energy-norm external flux energy eigenvalues") )
         CALL check( nf90_def_var(mncid,"W_xe_amp",nf90_double,
     $                         (/wdid/),wa_id) )
         CALL check( nf90_put_att(mncid,wa_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode amplifications") )
         CALL check( nf90_def_var(mncid,"W_xe_energyv",nf90_double,
     $                         (/wdid/),wev_id) )
         CALL check( nf90_put_att(mncid,wev_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode vacuum energy") )
         CALL check( nf90_def_var(mncid,"W_xe_energys",nf90_double,
     $                         (/wdid/),wes_id) )
         CALL check( nf90_put_att(mncid,wes_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode surface energy") )
         CALL check( nf90_def_var(mncid,"W_xe_energyp",nf90_double,
     $                         (/wdid/),wep_id) )
         CALL check( nf90_put_att(mncid,wep_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode total energy") )

         CALL check( nf90_def_dim(mncid,"mode_R",mpert,   rdid) )
         CALL check( nf90_def_var(mncid,"mode_R",nf90_int,rdid,mr_id))
         CALL check( nf90_put_att(mncid, mr_id ,"long_name",
     $    "Energy-norm external flux reluctance eigenmode index"))
         CALL check( nf90_def_var(mncid,"R_xe",nf90_double,
     $               (/mdid,rdid,idid/),rm_id) )
         CALL check( nf90_put_att(mncid,rm_id,"long_name",
     $    "Energy-norm external flux reluctance matrix") )
         CALL check( nf90_def_var(mncid,"R_xe_eigenvector",nf90_double,
     $               (/mdid,rdid,idid/),r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $    "Energy-norm external flux reluctance eigendecomposition") )
         CALL check( nf90_def_var(mncid,"R_xe_eigenvalue",nf90_double,
     $                         (/rdid/),re_id) )
         CALL check( nf90_put_att(mncid,re_id,"long_name",
     $    "Energy-norm external flux reluctance eigenvalues") )
         CALL check( nf90_put_att(mncid,re_id,"units","A/T") )
         CALL check( nf90_def_var(mncid,"R_xe_RL",nf90_double,
     $                         (/rdid/),rl_id) )
         CALL check( nf90_put_att(mncid,rl_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode RL-normalized") )
         CALL check( nf90_def_var(mncid,"R_xe_energyv",nf90_double,
     $                         (/rdid/),rev_id) )
         CALL check( nf90_put_att(mncid,rev_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode vacuum energy") )
         CALL check( nf90_def_var(mncid,"R_xe_energys",nf90_double,
     $                         (/rdid/),res_id) )
         CALL check( nf90_put_att(mncid,res_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode surface energy") )
         CALL check( nf90_def_var(mncid,"R_xe_energyp",nf90_double,
     $                         (/rdid/),rep_id) )
         CALL check( nf90_put_att(mncid,rep_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode total energy") )

         CALL check( nf90_def_dim(mncid,"mode_P",mpert,   pdid) )
         CALL check( nf90_def_var(mncid,"mode_P",nf90_int,pdid,mp_id))
         CALL check( nf90_put_att(mncid, mp_id ,"long_name",
     $    "Energy-norm external flux permeability eigenmode index") )
         CALL check( nf90_def_var(mncid,"P_xe",nf90_double,
     $               (/mdid,pdid,idid/),pm_id) )
         CALL check( nf90_put_att(mncid,pm_id,"long_name",
     $    "Energy-norm external flux permeability matrix") )
         CALL check( nf90_def_var(mncid,"P_xe_eigenvector",nf90_double,
     $               (/mdid,pdid,idid/),p_id) )
         CALL check( nf90_put_att(mncid,p_id,"long_name",
     $    "Energy-norm external flux permeability eigendecomposition") )
         CALL check( nf90_def_var(mncid,"P_xe_eigenvalue",nf90_double,
     $               (/pdid/),pe_id) )
         CALL check( nf90_put_att(mncid,pe_id,"long_name",
     $    "Energy-norm external flux permeability eigenvalues") )

         CALL check( nf90_def_var(mncid,"O_WR",nf90_double,
     $               (/wdid,rdid/),wr_id) )
         CALL check( nf90_put_att(mncid,wr_id,"long_name",
     $    "Overlap of energy and reluctance eigendecompositions") )
         CALL check( nf90_def_var(mncid,"O_WP",nf90_double,
     $               (/wdid,pdid/),wp_id) )
         CALL check( nf90_put_att(mncid,wp_id,"long_name",
     $    "Overlap of energy and permeability eigendecompositions") )
         CALL check( nf90_def_var(mncid,"O_RP",nf90_double,
     $               (/rdid,pdid/),rp_id) )
         CALL check( nf90_put_att(mncid,rp_id,"long_name",
     $    "Overlap of reluctance and permeability eigendecompositions"))

         CALL check( nf90_def_var(mncid,"O_Xxi_n",nf90_double,
     $               (/xdid,idid/),xt_id) )
         CALL check( nf90_put_att(mncid,xt_id,"long_name","Total "//
     $    "displacement decomposed in energy eigenmodes") )
         CALL check( nf90_def_var(mncid,"O_WPhi_xe",nf90_double,
     $               (/xdid,idid/),wx_id) )
         CALL check( nf90_put_att(mncid,wx_id,"long_name","Energy "//
     $    "normalized external flux decomposed in energy eigenmodes") )
         CALL check( nf90_def_var(mncid,"O_RPhi_xe",nf90_double,
     $               (/rdid,idid/),rx_id) )
         CALL check( nf90_put_att(mncid,rx_id,"long_name","Energy no"//
     $    "rmalized external flux decomposed in reluctance eigenmodes"))
         CALL check( nf90_def_var(mncid,"O_PPhi_xe",nf90_double,
     $               (/pdid,idid/),px_id) )
         CALL check( nf90_put_att(mncid,px_id,"long_name","Energy no"//
     $  "rmalized external flux decomposed in permeability eigenmodes"))

         IF(singcoup_set)THEN
            CALL check( nf90_def_var(mncid,"C_xe",nf90_double,
     $                  (/mdid,sdid,idid/),sc_id) )
            CALL check( nf90_put_att(mncid,sc_id,"long_name",
     $       "Energy normalized external flux singular-coupling") )
            CALL check( nf90_def_var(mncid,"C_xe_eigenvector",
     $                  nf90_double,(/mdid,sdid,idid/),s_id) )
            CALL check( nf90_put_att(mncid,s_id,"long_name",
     $       "Energy normalized external flux singular-coupling "//
     $       "SVD right-singular vectors") )
            CALL check( nf90_def_var(mncid,"C_xe_eigenvalue",
     $                            nf90_double,(/sdid/),se_id) )
            CALL check( nf90_put_att(mncid,se_id,"long_name",
     $       "Energy normalized external flux singular-coupling "//
     $       "SVD singular values") )
            CALL check( nf90_put_att(mncid,se_id,"units","unitless") )
            CALL check( nf90_def_var(mncid,"C_xe_energyv",nf90_double,
     $                            (/sdid/),sev_id) )
            CALL check( nf90_put_att(mncid,sev_id,"long_name",
     $       "Singular-coupling eigenmode vacuum energy") )
            CALL check( nf90_def_var(mncid,"C_xe_energys",nf90_double,
     $                            (/sdid/),ses_id) )
            CALL check( nf90_put_att(mncid,ses_id,"long_name",
     $       "Singular-coupling eigenmode surface energy") )
            CALL check( nf90_def_var(mncid,"C_xe_energyp",nf90_double,
     $                            (/sdid/),sep_id) )
            CALL check( nf90_put_att(mncid,sep_id,"long_name",
     $       "Singular-coupling eigenmode total energy") )
            CALL check( nf90_def_var(mncid,"O_WC",nf90_double,
     $                       (/wdid,sdid/),ws_id) )
            CALL check( nf90_put_att(mncid,ws_id,"long_name",
     $       "Overlap of energy and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_RC",nf90_double,
     $                       (/rdid,sdid/),rs_id) )
            CALL check( nf90_put_att(mncid,rs_id,"long_name",
     $       "Overlap of reluctance and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_PC",nf90_double,
     $                       (/pdid,sdid/),ps_id) )
            CALL check( nf90_put_att(mncid,ps_id,"long_name",
     $       "Overlap of permeability and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_CPhi_xe",nf90_double,
     $                       (/sdid,idid/),sx_id) )
            CALL check( nf90_put_att(mncid,sx_id,"long_name",
     $       "Energy normalized external flux decomposed in singular "//
     $       "coupling modes"))

         ENDIF

         CALL check( nf90_def_var(mncid,"Phi_xe",nf90_double,
     $               (/mdid,idid/),ex_id) )
         CALL check( nf90_put_att(mncid,ex_id,"units","T") )
         CALL check( nf90_put_att(mncid,ex_id,"long_name",
     $    "Energy-norm external field") )
         CALL check( nf90_def_var(mncid,"Phi_x",nf90_double,
     $               (/mdid,idid/),fx_id) )
         CALL check( nf90_put_att(mncid,fx_id,"units","Wb") )
         CALL check( nf90_put_att(mncid,fx_id,"long_name",
     $    "External flux") )
         CALL check( nf90_def_var(mncid,"Phi_e",nf90_double,
     $               (/mdid,idid/),et_id) )
         CALL check( nf90_put_att(mncid,et_id,"units","T") )
         CALL check( nf90_put_att(mncid,et_id,"long_name",
     $    "Energy-norm total field") )
         CALL check( nf90_def_var(mncid,"Phi",nf90_double,
     $               (/mdid,idid/),ft_id) )
         CALL check( nf90_put_att(mncid,ft_id,"units","Wb") )
         CALL check( nf90_put_att(mncid,ft_id,"long_name",
     $    "Total flux") )
         IF(coil_flag) THEN
           CALL check( nf90_def_dim(mncid,"coil_index",coil_num,cdid) )
           CALL check( nf90_def_var(mncid,"coil_index",
     $                 nf90_int,cdid,mc_id))
           CALL check( nf90_def_dim(mncid,"coil_strlen",24,sl_id) )
           CALL check( nf90_def_var(mncid,"coil_name",nf90_char,
     $                 (/sl_id,cdid/),cn_id) )
           CALL check( nf90_def_var(mncid,"Phi_coil",nf90_double,
     $                 (/mdid,cdid,idid/),fc_id) )
           CALL check( nf90_put_att(mncid,fc_id,"units","Wb") )
           CALL check( nf90_put_att(mncid,fc_id,"long_name",
     $      "Coil flux") )
         ENDIF

         IF(fun_flag)THEN
            CALL check( nf90_def_var(mncid,"Phi_xe_fun",nf90_double,
     $                  (/tdid,idid/),exf_id) )
            CALL check( nf90_put_att(mncid,exf_id,"units","T") )
            CALL check( nf90_put_att(mncid,exf_id,"long_name",
     $       "Energy-norm external flux") )
            CALL check( nf90_def_var(mncid,"Phi_x_fun",nf90_double,
     $                  (/tdid,idid/),fxf_id) )
            CALL check( nf90_put_att(mncid,fxf_id,"units","Wb") )
            CALL check( nf90_put_att(mncid,fxf_id,"long_name",
     $       "External flux") )
            CALL check( nf90_def_var(mncid,"Phi_e_fun",nf90_double,
     $                  (/tdid,idid/),etf_id) )
            CALL check( nf90_put_att(mncid,etf_id,"units","T") )
            CALL check( nf90_put_att(mncid,etf_id,"long_name",
     $       "Energy-norm total flux") )
            CALL check( nf90_def_var(mncid,"Phi_fun",nf90_double,
     $                  (/tdid,idid/),ftf_id) )
            CALL check( nf90_put_att(mncid,ftf_id,"units","Wb") )
            CALL check( nf90_put_att(mncid,ftf_id,"long_name",
     $       "Total flux") )
            IF(coil_flag) THEN
              CALL check( nf90_def_var(mncid,"Phi_coil_fun",nf90_double,
     $                    (/tdid, cdid, idid/),fcf_id) )
              CALL check( nf90_put_att(mncid,fcf_id,"units","Wb") )
              CALL check( nf90_put_att(mncid,fcf_id,"long_name",
     $         "Coil flux") )
            ENDIF
            CALL check( nf90_def_var(mncid,"W_xe_eigenvector_fun",
     $               nf90_double,(/tdid,wdid,idid/),wf_id) )
            CALL check( nf90_put_att(mncid,wf_id,"long_name",
     $         "Energy-norm external flux energy eigenmodes") )
            CALL check( nf90_def_var(mncid,"R_xe_eigenvector_fun",
     $               nf90_double,(/tdid,rdid,idid/),rf_id) )
            CALL check( nf90_put_att(mncid,rf_id,"long_name",
     $         "Energy-norm external flux reluctance eigenmodes") )
            IF(singcoup_set)THEN
               CALL check( nf90_def_var(mncid,"C_xe_eigenvector_fun",
     $                  nf90_double,(/tdid,sdid,idid/),sf_id) )
               CALL check( nf90_put_att(mncid,sf_id,"long_name",
     $            "Energy-norm external flux resonant-coupling modes") )
            ENDIF
         ENDIF
         ! End definitions
         CALL check( nf90_enddef(mncid) )

         IF(debug_flag) PRINT *,"  Putting variables"
         ! dimensions
         indx = (/(i,i=1,mpert)/)
         CALL check( nf90_put_var(mncid,mx_id,indx) )
         CALL check( nf90_put_var(mncid,mw_id,indx) )
         CALL check( nf90_put_var(mncid,mr_id,indx) )
         CALL check( nf90_put_var(mncid,mp_id,indx) )

         ! energy normalized matrices
         CALL check( nf90_put_var(mncid,wm_id,RESHAPE((/REAL(wmat),
     $               AIMAG(wmat)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,rm_id,RESHAPE((/REAL(rmat),
     $               AIMAG(rmat)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,pm_id,RESHAPE((/REAL(pmat),
     $               AIMAG(pmat)/),(/mpert,mpert,2/))) )

         ! Basis vectors and values
         CALL check( nf90_put_var(mncid,x_id,RESHAPE((/REAL(xvecs),
     $               AIMAG(xvecs)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,xe_id,xvals) )
         CALL check( nf90_put_var(mncid,w_id,RESHAPE((/REAL(wvecs),
     $               AIMAG(wvecs)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,we_id,wvals) )
         CALL check( nf90_put_var(mncid,wa_id,avals) )
         CALL check( nf90_put_var(mncid,r_id,RESHAPE((/REAL(rvecs),
     $               AIMAG(rvecs)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,re_id,rvals) )
         CALL check( nf90_put_var(mncid,p_id,RESHAPE((/REAL(pvecs),
     $               AIMAG(pvecs)/),(/mpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,pe_id,pvals) )
         IF(singcoup_set) THEN
            matms = TRANSPOSE(singcoupmat)
            CALL check( nf90_put_var(mncid,sc_id,RESHAPE(
     $       (/REAL(matms),AIMAG(matms)/),(/mpert,msing,2/))))
            CALL check( nf90_put_var(mncid,s_id,RESHAPE(
     $       (/REAL(svecs),AIMAG(svecs)/),(/mpert,msing,2/))))
            CALL check( nf90_put_var(mncid,se_id,svals) )
         ENDIF

         ! Energies for postprocessing re-normalization
         mat=MATMUL(MATMUL(ptof,plas_indinvmats(resp_index,:,:)),ptof)
         matmm=MATMUL(MATMUL(ptof,surf_indinvmats),ptof)
         tempmm=plas_indmats(resp_index,:,:)
         tempmm=MATMUL(CONJG(TRANSPOSE(tempmm)),surf_indinvmats)
         tempmm=MATMUL(tempmm,plas_indinvmats(resp_index,:,:))
         tempmm=MATMUL(MATMUL(ptof,tempmm),ptof) ! P-dagger.L^-1.P
         DO i=1,mpert
            tempm = rvecs(:,i)
            sengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(tempmm,tempm)))
            vengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(matmm,tempm)))
            pengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(mat,tempm)))
         ENDDO
         CALL check( nf90_put_var(mncid,res_id,sengys) )
         CALL check( nf90_put_var(mncid,rev_id,vengys) )
         CALL check( nf90_put_var(mncid,rep_id,pengys) )
         DO i=1,mpert
            tempm = wvecs(:,i)
            sengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(tempmm,tempm)))
            vengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(matmm,tempm)))
            pengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(mat,tempm)))
         ENDDO
         CALL check( nf90_put_var(mncid,wes_id,sengys) )
         CALL check( nf90_put_var(mncid,wev_id,vengys) )
         CALL check( nf90_put_var(mncid,wep_id,pengys) )
         IF(singcoup_set) THEN
            DO i=1,msing
               tempm = svecs(:,i)
               sengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(tempmm,tempm)))
               vengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(matmm,tempm)))
               pengys(i) =REAL(DOT_PRODUCT(tempm,MATMUL(mat,tempm)))
            ENDDO
            CALL check( nf90_put_var(mncid,ses_id,sengys(1:msing)) )
            CALL check( nf90_put_var(mncid,sev_id,vengys(1:msing)) )
            CALL check( nf90_put_var(mncid,sep_id,pengys(1:msing)) )
         ENDIF

         ! L.rho = -(1+s)/s (Boozer, single mode model)
         ! We calculate the "effective energy" Phi'.Phi_eff' for each mode where
         ! Phi_eff' = (W^-1.L.W^-1.W.R.W) .W^-1.Phi
         ! Note, the () term is unitless as it should be and we use R'.Phi' = e*Phi'
         matmm=MATMUL(MATMUL(ftop,surf_indmats),ftop)
         DO i=1,mpert
            tempm = rvecs(:,i)
            rlvals(i) = rvals(i) *
     $          REAL(DOT_PRODUCT(tempm,MATMUL(matmm,tempm)))
         ENDDO
         CALL check( nf90_put_var(mncid,rl_id,rlvals) )

         ! Basis intersections
         mat = MATMUL(CONJG(TRANSPOSE(wvecs)),rvecs)
         CALL check( nf90_put_var(mncid,wr_id,ABS(mat)) )
         mat = MATMUL(CONJG(TRANSPOSE(wvecs)),pvecs)
         CALL check( nf90_put_var(mncid,wp_id,ABS(mat)) )
         mat = MATMUL(CONJG(TRANSPOSE(rvecs)),pvecs)
         CALL check( nf90_put_var(mncid,rp_id,ABS(mat)) )
         IF(singcoup_set) THEN
            matms = MATMUL(CONJG(TRANSPOSE(wvecs)),svecs)
            CALL check( nf90_put_var(mncid,ws_id,ABS(matms)) )
            matms = MATMUL(CONJG(TRANSPOSE(rvecs)),svecs)
            CALL check( nf90_put_var(mncid,rs_id,ABS(matms)) )
            matms = MATMUL(CONJG(TRANSPOSE(pvecs)),svecs)
            CALL check( nf90_put_var(mncid,ps_id,ABS(matms)) )
         ENDIF

         ! Energy normalized flux
         temp = MATMUL(ftop,foutmn)
         CALL check( nf90_put_var(mncid,et_id,RESHAPE((/REAL(temp),
     $               AIMAG(temp)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,ft_id,RESHAPE((/REAL(foutmn),
     $               AIMAG(foutmn)/),(/mpert,2/))) )
         IF(fun_flag)THEN
           DO i=0,mthsurf
             CALL bicube_eval(rzphi,psilim,REAL(i,r8)/mthsurf,0)
             dphi(i) = rzphi%f(3)
           ENDDO
           CALL iscdftb(mfac,mpert,tempfun,mthsurf,temp)
           tempfun = tempfun * EXP(ifac * nn* dphi)
           CALL check(nf90_put_var(mncid,etf_id,RESHAPE((/REAL(tempfun),
     $         -helicity*AIMAG(tempfun)/),(/mthsurf+1,2/))) )
           CALL iscdftb(mfac,mpert,tempfun,mthsurf,foutmn)
           tempfun = tempfun * EXP(ifac * nn* dphi)
           CALL check(nf90_put_var(mncid,ftf_id,RESHAPE((/REAL(tempfun),
     $         -helicity*AIMAG(tempfun)/),(/mthsurf+1,2/))) )
         ENDIF
         temp = MATMUL(ftop,finmn)
         CALL check( nf90_put_var(mncid,ex_id,RESHAPE((/REAL(temp),
     $               AIMAG(temp)/),(/mpert,2/))) )
         CALL check( nf90_put_var(mncid,fx_id,RESHAPE((/REAL(finmn),
     $               AIMAG(finmn)/),(/mpert,2/))) )
         IF(fun_flag)THEN
           CALL iscdftb(mfac,mpert,tempfun,mthsurf,foutmn)
           tempfun = tempfun * EXP(ifac * nn* dphi)
           CALL check(nf90_put_var(mncid,exf_id,RESHAPE((/REAL(tempfun),
     $         -helicity*AIMAG(tempfun)/),(/mthsurf+1,2/))) )
           CALL iscdftb(mfac,mpert,tempfun,mthsurf,finmn)
           tempfun = tempfun * EXP(ifac * nn* dphi)
           CALL check(nf90_put_var(mncid,fxf_id,RESHAPE((/REAL(tempfun),
     $         -helicity*AIMAG(tempfun)/),(/mthsurf+1,2/))) )
         ENDIF
         IF(coil_flag) THEN
           CALL check( nf90_put_var(mncid,mc_id,(/(i,i=1,coil_num)/)) )
           CALL check( nf90_put_var(mncid,cn_id,coil_name(1:coil_num)))
           CALL check( nf90_put_var(mncid,fc_id,
     $                 RESHAPE((/REAL(coil_indmat),AIMAG(coil_indmat)/),
     $                         (/mpert,coil_num,2/))) )
           IF(fun_flag)THEN
             DO i=1,coil_num
               CALL iscdftb(mfac,mpert,tempfun,mthsurf,coil_indmat(:,i))
               tempfunc(:,i) = tempfun * EXP(ifac * nn* dphi)
             ENDDO
             CALL check( nf90_put_var(mncid, fcf_id,
     $            RESHAPE((/REAL(tempfunc), -helicity*AIMAG(tempfunc)/),
     $            (/mthsurf+1, coil_num, 2/))) )
           ENDIF
         ENDIF


         ! Decomposition of the applied flux
         temp = MATMUL(ftop, finmn) ! flux to power-normalized field (T)
         tempm = MATMUL(CONJG(TRANSPOSE(wvecs)),temp)
         CALL check( nf90_put_var(mncid,wx_id,RESHAPE((/REAL(tempm),
     $               AIMAG(tempm)/),(/mpert,2/))) )
         tempm = MATMUL(CONJG(TRANSPOSE(rvecs)),temp)
         CALL check( nf90_put_var(mncid,rx_id,RESHAPE((/REAL(tempm),
     $               AIMAG(tempm)/),(/mpert,2/))) )
         tempm = MATMUL(CONJG(TRANSPOSE(pvecs)),temp)
         CALL check( nf90_put_var(mncid,px_id,RESHAPE((/REAL(tempm),
     $               AIMAG(tempm)/),(/mpert,2/))) )
         IF(singcoup_set) THEN
            temps =MATMUL(CONJG(TRANSPOSE(svecs)),temp)
            CALL check( nf90_put_var(mncid,sx_id,RESHAPE((/REAL(temps),
     $               AIMAG(temps)/),(/msing,2/))) )
         ENDIF

         ! Decomposition of displacement
         temp = foutmn/(chi1*twopi*ifac*(mfac-nn*qlim)) ! total displacement
         tempm = MATMUL(CONJG(TRANSPOSE(xvecs)),temp)
         CALL check( nf90_put_var(mncid,xt_id,RESHAPE((/REAL(tempm),
     $               AIMAG(tempm)/),(/mpert,2/))) )

         ! Eigenmodes in real space (R,z) written in gpout_control
         IF(fun_flag)THEN
            DO i=1,mpert
               CALL iscdftb(mfac,mpert,wfuns(:,i),mthsurf,wvecs(:,i))
               wfuns(:,i) = wfuns(:,i) * EXP(ifac * nn * dphi)
               CALL iscdftb(mfac,mpert,rfuns(:,i),mthsurf,rvecs(:,i))
               rfuns(:,i) = rfuns(:,i) * EXP(ifac * nn * dphi)
               IF(singcoup_set .AND. (i<=msing))THEN
                 CALL iscdftb(mfac,mpert,sfuns(:,i),mthsurf,svecs(:,i))
                 sfuns(:,i) = sfuns(:,i) * EXP(ifac * nn * dphi)
               ENDIF
            ENDDO
            CALL check( nf90_put_var(mncid,wf_id,RESHAPE((/REAL(wfuns),
     $             -helicity*AIMAG(wfuns)/),(/mthsurf+1,mpert,2/))) )
            CALL check( nf90_put_var(mncid,rf_id,RESHAPE((/REAL(rfuns),
     $             -helicity*AIMAG(rfuns)/),(/mthsurf+1,mpert,2/))) )
            IF(singcoup_set)
     $       CALL check( nf90_put_var(mncid,sf_id,RESHAPE((/REAL(sfuns),
     $             -helicity*AIMAG(sfuns)/),(/mthsurf+1,msing,2/))) )
         ENDIF

         ! close the file
         CALL check( nf90_close(mncid) )
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE gpout_control_filter
c-----------------------------------------------------------------------
c     subprogram 18. gpout_qrv.
c     Add some basic alternative x coordinates into the profile outputs.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_qrv
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,p_id, qs_id,rn_id,dv_id

      REAL(r8), DIMENSION(mstep) :: qs_mstep, rn_mstep, dv_mstep
      REAL(r8), DIMENSION(0:mthsurf) :: unitfun
c-----------------------------------------------------------------------
c     calculation
c-----------------------------------------------------------------------
      unitfun = 1
      DO i=1,mstep
         CALL spline_eval(sq, psifac(i),0)
         qs_mstep(i) = sq%f(4)
         dv_mstep(i) = sq%f(3)
         rn_mstep(i) = issurfint(unitfun,mthsurf,psifac(i),3,1)
      ENDDO
      rn_mstep = rn_mstep / issurfint(unitfun,mthsurf,REAL(1.0,r8),3,1)
c-----------------------------------------------------------------------
c     ouput
c-----------------------------------------------------------------------
      CALL check( nf90_open(fncfile,nf90_write,fncid) )
      CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
      CALL check( nf90_redef(fncid))
      CALL check( nf90_def_var(fncid,"q",nf90_double,p_id,qs_id) )
      CALL check( nf90_put_att(fncid,qs_id,"long_name",
     $   "Safety factor") )
      CALL check( nf90_def_var(fncid,"rho",nf90_double,p_id,rn_id) )
      CALL check( nf90_put_att(fncid,rn_id,"long_name",
     $   "Normalized flux surface average minor radius") )
      CALL check( nf90_def_var(fncid,"dvdpsi_n",nf90_double,p_id,dv_id))
      CALL check( nf90_put_att(fncid,dv_id,"long_name",
     $   "Differential volume per normalized poloidal flux") )
      CALL check( nf90_put_att(fncid,dv_id,"units","m^3") )
      CALL check( nf90_enddef(fncid) )

      CALL check( nf90_put_var(fncid,dv_id,dv_mstep) )
      CALL check( nf90_put_var(fncid,qs_id,qs_mstep) )
      CALL check( nf90_put_var(fncid,rn_id,rn_mstep) )
      CALL check( nf90_close(fncid) )

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpout_qrv
c-----------------------------------------------------------------------
c     subprogram 19. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE check(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c-----------------------------------------------------------------------
c     subprogram 20. gpout_init_netcdf.
c     Initialize the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_init_netcdf
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER:: i,midid,mmdid,mcdid,mldid,modid,medid,mtdid,mjdid,
     $   mivid,mmvid,mcvid,mlvid,movid,mevid,mtvid,mrvid,mjvid,mqvid,
     $   fidid,fmdid,fpdid,ftdid,
     $   fivid,fmvid,fpvid,ftvid,frvid,frdid,fqvid,fq1vid,
     $   cidid,crdid,czdid,
     $   clvid,civid,crvid,czvid,
     $   id,fileids(3),ising
      INTEGER, DIMENSION(mpert) :: mmodes
c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *,"Initializing NETCDF files"
      mmodes = (/(i,i=1,mpert)/)
      mncfile = "gpec_control_output_n"//TRIM(sn)//".nc"
      fncfile = "gpec_profile_output_n"//TRIM(sn)//".nc"
      cncfile = "gpec_cylindrical_output_n"//TRIM(sn)//".nc"
c-----------------------------------------------------------------------
c     open files
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Creating netcdf files"
      CALL check( nf90_create(mncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=mncid) )
      CALL check( nf90_create(fncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=fncid) )
      CALL check( nf90_create(cncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=cncid) )
c-----------------------------------------------------------------------
c     define global file attributes
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Defining modal netcdf globals"
      CALL check( nf90_put_att(mncid,nf90_global,"title",
     $     "GPEC outputs in Fourier or alternate modal bases"))
      CALL check( nf90_def_dim(mncid,"i",2,       midid) )
      CALL check( nf90_def_var(mncid,"i",nf90_int,midid,mivid) )
      CALL check( nf90_def_dim(mncid,"m",mpert,  mmdid) )
      CALL check( nf90_def_var(mncid,"m",nf90_int,mmdid,mmvid) )
      CALL check( nf90_def_dim(mncid,"m_out",lmpert, modid) )
      CALL check( nf90_def_var(mncid,"m_out",nf90_int,modid,movid) )
      CALL check( nf90_def_dim(mncid,"mode",mpert,   medid) )
      CALL check( nf90_def_var(mncid,"mode",nf90_int,medid,mevid))
      CALL check( nf90_def_dim(mncid,"theta",mthsurf+1,  mtdid) )
      CALL check( nf90_def_var(mncid,"theta",nf90_double,mtdid,mtvid) )
      IF(msing>0)THEN
         CALL check( nf90_def_dim(mncid,"psi_n_rational", msing, mjdid))
         CALL check( nf90_def_var(mncid,"psi_n_rational", nf90_double,
     $      mjdid, mjvid) )
         CALL check( nf90_def_var(mncid,"q_rational",nf90_double, mjdid,
     $      mqvid) )
         CALL check( nf90_def_var(mncid,"m_rational", nf90_int, mjdid,
     $      mrvid) )
         CALL check( nf90_def_dim(mncid,"mode_C", msing, mcdid) )
         CALL check( nf90_def_var(mncid,"mode_C", nf90_double,
     $      mcdid, mcvid) )
         CALL check( nf90_put_att(mncid, mcvid ,"long_name",
     $    "Singular coupling mode index") )
         CALL check( nf90_def_dim(mncid,"mode_C_local", osing, mldid) )
         CALL check( nf90_def_var(mncid,"mode_C_local", nf90_double,
     $      mldid,mlvid) )
         CALL check( nf90_put_att(mncid, mlvid ,"long_name",
     $    "Local singular coupling mode index") )
      ENDIF
      CALL check( nf90_put_att(mncid,nf90_global,"jacobian",jac_type))
      CALL check( nf90_put_att(mncid,nf90_global,"helicity",helicity))
      CALL check( nf90_put_att(mncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(mncid,nf90_global,"psi_n_lim",psilim))
      CALL check( nf90_put_att(mncid,nf90_global,"chi1",chi1))

      IF(debug_flag) PRINT *," - Defining flux netcdf globals"
      CALL check( nf90_put_att(fncid,nf90_global,"title",
     $     "GPEC outputs in magnetic coordinate systems"))
      CALL check( nf90_def_dim(fncid,"i",2, fidid) )
      CALL check( nf90_def_var(fncid,"i",nf90_int,fidid,fivid) )
      CALL check( nf90_def_dim(fncid,"m_out",lmpert, fmdid) )
      CALL check( nf90_def_var(fncid,"m_out",nf90_int,fmdid,fmvid) )
      CALL check( nf90_def_dim(fncid,"psi_n",mstep, fpdid) )
      CALL check( nf90_def_var(fncid,"psi_n",nf90_double,fpdid,fpvid) )
      CALL check( nf90_def_dim(fncid,"theta_dcon",mthsurf+1, ftdid) )
      CALL check( nf90_def_var(fncid,"theta_dcon",nf90_double,
     $     ftdid,ftvid) )
      IF(msing>0)THEN
         CALL check( nf90_def_dim(fncid,"psi_n_rational",msing, frdid) )
         CALL check( nf90_def_var(fncid,"psi_n_rational",nf90_double,
     $        frdid,frvid))
         CALL check(nf90_def_var(fncid,"q_rational", nf90_double, frdid,
     $        fqvid))
         CALL check(nf90_def_var(fncid,"dqdpsi_n_rational",nf90_double,
     $        frdid,fq1vid))
      ENDIF
      CALL check( nf90_put_att(fncid,nf90_global,"helicity",helicity))

      IF(debug_flag) PRINT *," - Defining cylindrical netcdf globals"
      CALL check( nf90_put_att(cncid,nf90_global,"title",
     $     "GPEC outputs in (R,z) coordinates"))
      CALL check( nf90_def_dim(cncid,"i",2,       cidid) )
      CALL check( nf90_def_var(cncid,"i",nf90_int,cidid,civid) )
      CALL check( nf90_def_dim(cncid,"R",nr+1,crdid) )
      CALL check( nf90_def_var(cncid,"R",nf90_double,crdid,crvid) )
      CALL check( nf90_put_att(cncid,crvid,"units","m") )
      CALL check( nf90_def_dim(cncid,"z",nz+1,czdid) )
      CALL check( nf90_def_var(cncid,"z",nf90_double,czdid,czvid) )
      CALL check( nf90_put_att(cncid,czvid,"units","m") )
      CALL check( nf90_def_var(cncid,"l",nf90_double,
     $                         (/crdid,czdid/),clvid) )

      ! Add common global attributes
      fileids = (/mncid,fncid,cncid/)
      DO i=1,3
         id = fileids(i)
         CALL check( nf90_put_att(id,nf90_global,"shot",INT(shotnum)) )
         CALL check( nf90_put_att(id,nf90_global,"time",INT(shottime)) )
         CALL check( nf90_put_att(id,nf90_global,"machine",machine) )
         CALL check( nf90_put_att(id,nf90_global,"n",nn) )
         CALL check( nf90_put_att(id,nf90_global,"version",version))
      ENDDO

      CALL check( nf90_enddef(mncid) )
      CALL check( nf90_enddef(fncid) )
      CALL check( nf90_enddef(cncid) )
c-----------------------------------------------------------------------
c     set dimensional variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *," - Putting coordinates in control netcdfs"
      CALL check( nf90_put_var(mncid,mivid,(/0,1/)) )
      CALL check( nf90_put_var(mncid,mmvid,mfac) )
      CALL check( nf90_put_var(mncid,movid,lmfac) )
      CALL check( nf90_put_var(mncid,mevid,mmodes) )
      CALL check( nf90_put_var(mncid,mtvid,theta) )
      IF(msing>0)THEN
         CALL check( nf90_put_var(mncid,mjvid,
     $      (/(singtype(i)%psifac,i=1,msing)/)) )
         CALL check( nf90_put_var(mncid,mrvid,
     $      (/(INT(singtype(i)%q*nn),i=1,msing)/)) )
         CALL check( nf90_put_var(mncid,mqvid,
     $      (/(singtype(i)%q,i=1,msing)/)) )
         CALL check( nf90_put_var(mncid, mcvid, (/(i,i=1,msing)/)) )
         CALL check( nf90_put_var(mncid, mlvid, (/(i,i=1,osing)/)) )
      ENDIF

      IF(debug_flag) PRINT *," - Putting coordinates in flux netcdfs"
      CALL check( nf90_put_var(fncid,fivid,(/0,1/)) )
      CALL check( nf90_put_var(fncid,fmvid,lmfac) )
      CALL check( nf90_put_var(fncid,fpvid,psifac(1:mstep)) )
      CALL check( nf90_put_var(fncid,ftvid,theta) )
      IF(msing>0)THEN
         CALL check( nf90_put_var(fncid,frvid,
     $      (/(singtype(i)%psifac,i=1,msing)/)) )
         CALL check( nf90_put_var(fncid,fqvid,
     $      (/(singtype(i)%q,i=1,msing)/)) )
         CALL check( nf90_put_var(fncid,fq1vid,
     $      (/(singtype(i)%q1,i=1,msing)/)) )
      ENDIF

      IF(debug_flag) PRINT *," - Putting coordinates in cyl. netcdfs"
      CALL check( nf90_put_var(cncid,civid,(/0,1/)) )
      CALL check( nf90_put_var(cncid,crvid,gdr(:,0)) )
      CALL check( nf90_put_var(cncid,czvid,gdz(0,:)) )
      CALL check( nf90_put_var(cncid,clvid,gdl(:,:)) )

      IF(debug_flag) PRINT *," - Closing netcdf files"
      CALL check( nf90_close(mncid) )
      CALL check( nf90_close(fncid) )
      CALL check( nf90_close(cncid) )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpout_init_netcdf
c-----------------------------------------------------------------------
c     subprogram 21. gpout_close_netcdf.
c     Close the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE gpout_close_netcdf
c-----------------------------------------------------------------------
c     close files
c-----------------------------------------------------------------------
      CALL check( nf90_close(mncid) )
      CALL check( nf90_close(fncid) )
      CALL check( nf90_close(cncid) )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpout_close_netcdf

      END MODULE gpout_mod
