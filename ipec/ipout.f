c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     write various output results of ipec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipout_mod
c      1. ipout_response
c      2. ipout_singcoup
c      3. ipout_control
c      4. ipout_singfld
c      5. ipout_vsingfld
c      6. ipout_pmodb
c      7. ipout_xbnormal
c      8. ipout_vbnormal
c      9. ipout_xtangent
c     10. ipout_xbrzphi
c     11. ipout_vsbrzphi
c     12. ipout_xbrzphifun
c     13. ipout_arzphifun
c-----------------------------------------------------------------------
c     subprogram 0. ipout_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipout_mod
      USE ipresp_mod
      USE ipvacuum_mod
      USE ipdiag_mod
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
      LOGICAL :: singcoup_set = .FALSE.
      REAL(r8) :: jarea
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: singcoup1mat,
     $   w1v,w2v,w3v,t1v,t2v,t3v,fldflxmat,singbnoflxs

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
c     subprogram 1. ipout_response.
c     write basic information.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_response(rout,bpout,bout,rcout,tout,jout)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout,jout
      INTEGER :: lwork,i,j,k
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(mpert) :: sts,s,ev
      COMPLEX(r8), DIMENSION(mpert) :: cev
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,b,u,vt

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(3*mpert) :: work

      INTEGER :: idid,mdid,edid,l_id,r_id,la_id,p_id
      COMPLEX(r8), DIMENSION(lmpert) :: vL,vL1,vLi,vP,vP1,vR,vW,templ
      COMPLEX(r8), DIMENSION(mpert,mpert) :: matmm
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat
      COMPLEX(r8), DIMENSION(0:mthsurf,mpert) :: wtfun,ilfun,rfun,pfun
      CHARACTER(2048) :: header

      IF(timeit) CALL ipec_timer(-2)
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
         sts(i)=-et(i)/(surfei(i)+surfee(i))
         IF (sts(i) > 0) s(i)=-s(i)
         s(i)=-1/s(i)
      ENDDO

      CALL ascii_open(out_unit,"ipec_response_n"//
     $	   TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_RESPONSE: Response parameters"
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
      ! torque = -2*nn*imag(ep)
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
             CALL ipeq_bcoords(psilim,vL,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vL1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vLi,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vP,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vP1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vR,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL ipeq_bcoords(psilim,vW,lmfac,lmpert,
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
        CALL ascii_open(out_unit,"ipec_response_fun_n"//
     $    TRIM(sn)//".out","UNKNOWN")
        WRITE(out_unit,*)"IPEC_RESPONSE_FUN: "//
     $    "Eigenmodes on control surface in functions"
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
        WRITE(header,'(a16,9(1x,a16))')"r","z",
     $    "real(P)","imag(P)","real(rho)","imag(rho)",
     $    "real(W)","imag(W)","real(iLmda)","imag(iLmda)"
        DO j=1,mpert
          WRITE(out_unit,'(1x,2(a12,I4))')"isol =",j," n =",nn
          WRITE(out_unit,*)
          WRITE(out_unit,*) header
          DO i=0,mthsurf
            WRITE(out_unit,'(10(es17.8e3))') r(i),z(i),
     $        pfun(i,j),rfun(i,j),wtfun(i,j),ilfun(i,j)
          ENDDO
          WRITE(out_unit,*)
        ENDDO
        CALL ascii_close(out_unit)
      ENDIF


      ! Write to netCDF file
      IF(debug_flag) PRINT *,"Opening "//TRIM(mncfile)
      CALL check( nf90_open(mncfile,nf90_write,mncid) )
      CALL check( nf90_inq_dimid(mncid,"i",idid) )
      CALL check( nf90_inq_dimid(mncid,"m",mdid) )
      CALL check( nf90_inq_dimid(mncid,"mode",edid) )

      ! Start definitions
      CALL check( nf90_redef(mncid))
      CALL check( nf90_def_var(mncid,"L",nf90_double,
     $               (/mdid,edid,idid/),l_id) )
      CALL check( nf90_put_att(mncid,l_id,"long_name",
     $    "Surface Inductance") )      
      CALL check( nf90_def_var(mncid,"Lambda",nf90_double,
     $               (/mdid,edid,idid/),la_id) )
      CALL check( nf90_put_att(mncid,la_id,"long_name",
     $    "Plasma inductance") )
      CALL check( nf90_def_var(mncid,"P",nf90_double,
     $               (/mdid,edid,idid/),p_id) )
      CALL check( nf90_put_att(mncid,p_id,"long_name",
     $    "Permeability") )
      CALL check( nf90_def_var(mncid,"rho",nf90_double,
     $               (/mdid,edid,idid/),r_id) )
      CALL check( nf90_put_att(mncid,r_id,"long_name",
     $    "Reluctance") )
      ! End definitions
      CALL check( nf90_enddef(mncid) )
      ! Convert to output coords
      DO i=1,mpert
        templ = 0
        templ(mlow-lmlow+i) = 1.0
        IF((jac_out /= jac_type).OR.(tout==0)) CALL ipeq_bcoords(
     $     psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,0)
        coordmat(:,i) = templ
      ENDDO
      !matlm = MATMUL(coordmat,surf_indmats)
      matmm = surf_indmats
      CALL check( nf90_put_var(mncid,l_id,RESHAPE((/REAL(matmm),
     $            AIMAG(matmm)/),(/mpert,mpert,2/))) )
      !matlm = MATMUL(coordmat,plas_indmats(resp_index,:,:))
      matmm = plas_indmats(resp_index,:,:)
      CALL check( nf90_put_var(mncid,la_id,RESHAPE((/REAL(matmm),
     $            AIMAG(matmm)/),(/mpert,mpert,2/))) )
      !matlm = MATMUL(coordmat,permeabmats(resp_index,:,:))
      matmm = permeabmats(resp_index,:,:)
      CALL check( nf90_put_var(mncid,p_id,RESHAPE((/REAL(matmm),
     $            AIMAG(matmm)/),(/mpert,mpert,2/))) )
      !matlm = MATMUL(coordmat,reluctmats(resp_index,:,:))
      matmm = reluctmats(resp_index,:,:)
      CALL check( nf90_put_var(mncid,r_id,RESHAPE((/REAL(matmm),
     $            AIMAG(matmm)/),(/mpert,mpert,2/))) )
      CALL check( nf90_close(mncid) )

      ! log eigenvalues with harvest
      ierr=set_harvest_payload_dbl_array(hlog,"s_P"//nul,
     $     REAL(permeabev(resp_index,:)),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_L"//nul,
     $     surf_indev(:),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_Lambda"//nul,
     $     plas_indev(resp_index,:),mpert)
      ierr=set_harvest_payload_dbl_array(hlog,"s_rho"//nul,
     $     reluctev(resp_index,:),mpert)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_response
c-----------------------------------------------------------------------
c     subprogram 2. ipout_singcoup.
c     compute coupling between singular surfaces and external fields.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singcoup(spot,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      REAL(r8), INTENT(IN) :: spot

      REAL(r8), DIMENSION(msing) :: s,s1,s2,s3
      COMPLEX(r8), DIMENSION(msing,msing) :: u

      INTEGER :: i,j,itheta,ising,resnum,rsing,rpert,
     $     tmlow,tmhigh,tmpert,lwork,info
      REAL(r8) :: respsi,lpsi,rpsi,jarea,thetai
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,rfacs,jcfun,wcfun,
     $     dphi,thetas,units,jacs,jacfac

      COMPLEX(r8), DIMENSION(mpert) :: finmn,foutmn,
     $     fkaxmn,singflx_mn,ftnmn
      COMPLEX(r8), DIMENSION(lmpert) :: lftnmn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: ftnfun

      REAL(r8), DIMENSION(msing) :: area,j_c,w_c,shear

      COMPLEX(r8), DIMENSION(msing,mpert) :: deltas,delcurs,
     $     singcurs,islandhwids
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: convmat
      COMPLEX(r8), DIMENSION(msing,mpert,mpert) :: fsurfindmats
      COMPLEX(r8), DIMENSION(:), POINTER :: fldflxmn,bmn
      COMPLEX(r8), DIMENSION(:,:), POINTER ::  temp1,temp2,
     $     t1mat,t2mat,t3mat

      INTEGER, DIMENSION(:), POINTER :: ipiv,tmfac
      REAL(r8), DIMENSION(:), POINTER :: rwork
      COMPLEX(r8), DIMENSION(:), POINTER :: work
      COMPLEX(r8), DIMENSION(:,:), POINTER :: a,vt    

      TYPE(spline_type) :: spl 
c-----------------------------------------------------------------------
c     compute characteristic currents with normalization.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
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
         CALL ipvacuum_flxsurf(respsi)
         fsurfindmats(ising,:,:)=fsurf_indmats
         DEALLOCATE(fsurf_indev,fsurf_indmats)
      ENDDO
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      deltas=0
      delcurs=0
      singcurs=0
      CALL ipeq_alloc
      ALLOCATE(singbnoflxs(msing,mpert))
      DO i=1,mpert
         finmn=0
         finmn(i)=1.0
         CALL ipeq_weight(psilim,finmn,mfac,mpert,1)
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
c     evaluate delta/singular current/normal field/islands.
c-----------------------------------------------------------------------
         DO ising=1,msing
            resnum=NINT(singtype(ising)%q*nn)-mlow+1
            respsi=singtype(ising)%psifac
            lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
            CALL ipeq_sol(lpsi)
            lbwp1mn=bwp1_mn(resnum)
            
            rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1)) 
            CALL ipeq_sol(rpsi)
            rbwp1mn=bwp1_mn(resnum)

            deltas(ising,i)=rbwp1mn-lbwp1mn
            delcurs(ising,i)=j_c(ising)*deltas(ising,i)*ifac/
     $           (twopi*mfac(resnum))
            singcurs(ising,i)=-delcurs(ising,i)/ifac
            fkaxmn=0
            fkaxmn(resnum)=singcurs(ising,i)/(twopi*nn)

            singflx_mn=MATMUL(fsurfindmats(ising,:,:),fkaxmn)
            singbnoflxs(ising,i)=singflx_mn(resnum)/area(ising)
            islandhwids(ising,i)=4*singflx_mn(resnum)/
     $           (twopi*shear(ising)*sq%f(4)*chi1)
         ENDDO
         IF(verbose) WRITE(*,'(1x,a16,i4,a22,es10.3)')
     $        "poloidal mode = ",mfac(i),", resonant coupling = ",
     $        SUM(ABS(singbnoflxs(:,i)))/msing
      ENDDO
      CALL ipeq_dealloc
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
         ALLOCATE(fldflxmn(lmpert),fldflxmat(lmpert,lmpert),
     $        t1mat(msing,lmpert),t2mat(msing,lmpert),
     $        t3mat(msing,lmpert),tmfac(lmpert))
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
         t1mat = MATMUL(singbnoflxs,convmat)
         t2mat = MATMUL(singcurs,convmat)
         t3mat = MATMUL(islandhwids,convmat)
      ELSE
         ALLOCATE(fldflxmn(mpert),fldflxmat(mpert,mpert),
     $        t1mat(msing,mpert),t2mat(msing,mpert),
     $        t3mat(msing,mpert),tmfac(mpert))
         units = (/(1.0,itheta=0,mthsurf)/)
         jarea = issurfint(units,mthsurf,psilim,0,0)
         DO i=1,mpert
            ftnmn=0
            ftnmn(i)=1.0
            fldflxmn=ftnmn
            CALL ipeq_weight(psilim,fldflxmn,mfac,mpert,2)
            fldflxmat(:,i)=fldflxmn/sqrt(jarea)            
         ENDDO
         tmlow = mlow
         tmhigh = mhigh
         tmpert = mpert
         tmfac = mfac
         t1mat = singbnoflxs
         t2mat = singcurs*twopi*nn
         t3mat = islandhwids
      ENDIF
c-----------------------------------------------------------------------
c     svd analysis.
c-----------------------------------------------------------------------
      lwork=3*tmpert
      ALLOCATE(ipiv(tmpert),rwork(5*msing),work(lwork),
     $     a(msing,tmpert),vt(msing,tmpert),
     $     w1v(tmpert,msing),w2v(tmpert,msing),w3v(tmpert,msing),
     $     t1v(tmpert,msing),t2v(tmpert,msing),t3v(tmpert,msing),
     $     temp1(tmpert,tmpert),temp2(tmpert,tmpert),bmn(tmpert),
     $     singcoup1mat(msing,tmpert))
      
      temp1=fldflxmat
      temp2=0
      DO i=1,tmpert
         temp2(i,i)=1
      ENDDO

      CALL zgetrf(tmpert,tmpert,temp1,tmpert,ipiv,info)
      CALL zgetrs('N',tmpert,tmpert,temp1,tmpert,ipiv,temp2,tmpert,info)

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t1mat,temp2)
      singcoup1mat = a
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s1=s
      w1v=CONJG(TRANSPOSE(vt))
      t1v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t2mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s2=s
      w2v=CONJG(TRANSPOSE(vt))
      t2v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))

      work=0
      rwork=0
      s=0
      u=0
      vt=0
      a=MATMUL(t3mat,temp2)
      CALL zgesvd('S','S',msing,tmpert,a,msing,s,u,msing,vt,msing,
     $     work,lwork,rwork,info)    
      s3=s
      w3v=CONJG(TRANSPOSE(vt))
      t3v=MATMUL(CONJG(fldflxmat),CONJG(TRANSPOSE(vt)))
c-----------------------------------------------------------------------
c     save module wide variables
c-----------------------------------------------------------------------
      singcoup_set  = .TRUE.
c-----------------------------------------------------------------------
c     write matrix.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singcoup_matrix_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGCOUP_MATRIX: Coupling matrices"//
     $     " between resonant field and external field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout  
      WRITE(out_unit,'(4(1x,a12,I4))')
     $     "msing =",msing,"mpert =",tmpert,
     $     "mlow =",tmlow,"mhigh =",tmhigh
      WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $     "psilim =",psilim,"qlim =",qlim
      WRITE(out_unit,*)

      WRITE(out_unit,*)"Coupling matrix to resonant fields"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es17.8e3)')
     $        "q =",singtype(i)%q,"psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(es17.8e3))')
     $           tmfac(j),REAL(t1mat(i,j)),AIMAG(t1mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"Coupling matrix to singular currents"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es17.8e3)')
     $        "q =",singtype(i)%q,"psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(es17.8e3))')
     $           tmfac(j),REAL(t2mat(i,j)),AIMAG(t2mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"Coupling matrix to island half-widths"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a4,f6.3,1x,a6,es17.8e3)')
     $        "q =",singtype(i)%q,"psi =",singtype(i)%psifac
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real","imag"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(es17.8e3))')
     $           tmfac(j),REAL(t3mat(i,j)),AIMAG(t3mat(i,j))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO
      CALL ascii_close(out_unit)
      
      CALL ascii_open(out_unit,"ipec_singcoup_svd_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGCOUP_SVD: SVD analysis"//
     $     " for coupling matrices"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a12,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tmag_out  
      WRITE(out_unit,'(4(1x,a12,I4))')
     $     "msing =",msing,"mpert =",tmpert,
     $     "mlow =",tmlow,"mhigh =",tmhigh
      WRITE(out_unit,'(2(1x,a12,es17.8e3))')
     $     "psilim =",psilim,"qlim =",qlim
      WRITE(out_unit,*)

      WRITE(out_unit,*)"SVD for coupling matrix to resonant fields"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es17.8e3)')
     $        "mode =",i,"s =",s1(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(Phi)","imag(Phi)"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,4(es17.8e3))')
     $           tmfac(j),REAL(t1v(j,i)),AIMAG(t1v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"SVD for coupling matrix to singular currents"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es17.8e3)')
     $        "mode =",i,"s =",s2(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real(Phi)","imag(Phi)"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(es17.8e3))')
     $           tmfac(j),REAL(t2v(j,i)),AIMAG(t2v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO

      WRITE(out_unit,*)"SVD for coupling matrix to island half-widths"
      WRITE(out_unit,*) 
      DO i=1,msing
         WRITE(out_unit,'(1x,a6,I4,1x,a6,es17.8e3)')
     $        "mode =",i,"s =",s3(i)
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a4,2(1x,a16))')"m","real(Phi)","imag(Phi)"
         DO j=1,tmpert
            WRITE(out_unit,'(1x,I4,2(es17.8e3))')
     $           tmfac(j),REAL(t3v(j,i)),AIMAG(t3v(j,i))
         ENDDO 
         WRITE(out_unit,*) 
      ENDDO
      CALL ascii_close(out_unit)

      DEALLOCATE(tmfac,t1mat,t2mat,t3mat,fldflxmn)
      DEALLOCATE(ipiv,rwork,work,a,vt,temp1,temp2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_singcoup
c-----------------------------------------------------------------------
c     subprogram 3. ipout_control
c     calculate response from external field on the control surface.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_control(ifile,finmn,foutmn,xspmn,
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

      INTEGER :: i,j,i1,i2,i3,ms,itheta
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: vengy,sengy,pengy,area,thetai,scale,norm
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacs,rvecs,
     $     zvecs,sbinfun,sboutfun,units

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,xinmn,xoutmn,tempmn,
     $      abinmn
      COMPLEX(r8), DIMENSION(lmpert) :: cinmn,coutmn,cawmn,templ
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,xinfun,xoutfun
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat
      COMPLEX(r8), DIMENSION(mpert,lmpert) :: tempml

      REAL(r8), DIMENSION(:,:), POINTER :: dcosmn,dsinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn
      
      INTEGER :: i_id,m_id,modid,mo_id,t_id,r_id,z_id,rn_id,zn_id,p_id,
     $    x_id,xx_id,xm_id,xxm_id,bm_id,bxm_id,b_id,bx_id,jo_id
c-----------------------------------------------------------------------
c     compute control surface area.
c-----------------------------------------------------------------------
      units = (/(1.0,itheta=0,mthsurf)/)
      jarea = issurfint(units,mthsurf,psilim,0,0)
c-----------------------------------------------------------------------
c     check data_type and read data.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
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
               CALL ipec_stop(message)
            ENDIF            
         ENDDO
         
         CALL ascii_close(in_unit)
         rawmn=dcosmn+ifac*dsinmn
         cawmn=rawmn(:,nn)
         DEALLOCATE(dcosmn,dsinmn,rawmn)
      ELSE IF (harmonic_flag) THEN
         DO i=-hmnum,hmnum
            IF ((-mmin+i>=1).AND.(-mmin+i<=lmpert)) THEN
               cawmn(-mmin+i+1)=cosmn(i)+ifac*sinmn(i)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      IF (data_flag .OR. harmonic_flag) THEN
         CALL ipeq_fcoords(psilim,cawmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)             
         binmn=0
         DO i=1,lmpert
            IF ((lmlow-mlow+i>=1).AND.(lmlow-mlow+i<=mpert)) THEN
               binmn(lmlow-mlow+i)=cawmn(i)
            ENDIF
         ENDDO   
c-----------------------------------------------------------------------
c     convert to field if displacement is given.
c-----------------------------------------------------------------------
         IF (displacement_flag) THEN
            CALL ipeq_weight(psilim,binmn,mfac,mpert,5)
            binmn=twopi*ifac*chi1*(mfac-nn*qlim)*binmn
            CALL ipeq_weight(psilim,binmn,mfac,mpert,0)
         ENDIF 
         binmn=binmn*scale
         tempmn=binmn
         CALL ipeq_weight(psilim,tempmn,mfac,mpert,1)              
      ENDIF
      finmn=finmn+tempmn
c-----------------------------------------------------------------------
c     filter external flux
c-----------------------------------------------------------------------
      CALL ipout_control_filter(finmn,foutmn,filter_types,filter_modes,
     $           rout,bpout,bout,rcout,tout,jout,filter_out)
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
      CALL ipeq_weight(psilim,binmn,mfac,mpert,0)
      CALL ipeq_weight(psilim,boutmn,mfac,mpert,0)
      xinmn=finmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      xoutmn=xspmn
      CALL ipeq_weight(psilim,xinmn,mfac,mpert,4)
      CALL ipeq_weight(psilim,xoutmn,mfac,mpert,4)
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
      ierr=set_harvest_payload_dbl(hlog,"amplification"//nul,
     $                             sengy/pengy)
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_control_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_CONTROL: "//
     $     "Plasma response for an external perturbation on the "//
     $     "control surface"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a12,I4)')"mpert =",mpert
      WRITE(out_unit,'(1x,a16,es17.8e3)')"vacuum energy =",vengy
      WRITE(out_unit,'(1x,a16,es17.8e3)')"surface energy =",sengy
      WRITE(out_unit,'(1x,a16,es17.8e3)')"plasma energy =",pengy
      WRITE(out_unit,'(1x,a18,es17.8e3)')"toroidal torque =",
     $     -2*nn*AIMAG(py)
      WRITE(out_unit,*)

      WRITE(out_unit,*)"jac_type = "//jac_type
      WRITE(out_unit,*)

      WRITE(out_unit,'(1x,a4,8(1x,a16))') "m",
     $     "real(bin)","imag(bin)","real(bout)","imag(bout)",
     $     "real(Phi^x)","imag(Phi^x)","real(Phi)","imag(Phi)"
      DO i=1,mpert
         WRITE(out_unit,'(1x,I4,8(es17.8e3))')mfac(i),
     $        REAL(binmn(i)),AIMAG(binmn(i)),
     $        REAL(boutmn(i)),AIMAG(boutmn(i)),
     $        REAL(finmn(i)),AIMAG(finmn(i)),
     $        REAL(foutmn(i)),AIMAG(foutmn(i))
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
         CALL ipeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin) 
         CALL ipeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_in = ",jac_in,
     $        "jsurf_in =",jin,"tmag_in =",tin
         WRITE(out_unit,*)             
         WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)","imag(bin)",
     $     "real(bout)","imag(bout)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,4(es17.8e3))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i))
         ENDDO 
         WRITE(out_unit,*)       
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
         CALL ipeq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout) 
         CALL ipeq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout)
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_out = ",jac_out,
     $        "jsurf_out =",jout,"tmag_out =",tout  
         WRITE(out_unit,*)          
         WRITE(out_unit,'(1x,a4,4(1x,a16))')"m","real(bin)","imag(bin)",
     $        "real(bout)","imag(bout)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,6(es17.8e3))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i))
         ENDDO 
         WRITE(out_unit,*)
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
      CALL ascii_close(out_unit)


      ! netcdf output
      ! Convert to output coords
      DO i=1,mpert
        templ = 0
        templ(mlow-lmlow+i) = 1.0
        IF((jac_out /= jac_type).OR.(tout==0)) CALL ipeq_bcoords(
     $     psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,jout)
        coordmat(:,i) = templ
      ENDDO

      IF(debug_flag) PRINT *,"Opening "//TRIM(mncfile)
      CALL check( nf90_open(mncfile,nf90_write,mncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(mncid,"i",i_id) )
      CALL check( nf90_inq_dimid(mncid,"m",m_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(mncid))
      CALL check( nf90_def_dim(mncid,"m_out",lmpert,modid) )
      CALL check( nf90_def_var(mncid,"m_out",nf90_int,modid,mo_id))
      CALL check( nf90_put_att(mncid, mo_id ,"long_name",
     $    "Poloidal mode number of jacobian conversion") )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_vacuum",vengy) )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_surface",sengy) )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_plasma",pengy) )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "area",jarea) )
      CALL check( nf90_def_var(mncid,"b_n",nf90_double,
     $            (/m_id,i_id/),bm_id) )
      CALL check( nf90_put_att(mncid,bm_id,"units","Tesla") )
      CALL check( nf90_put_att(mncid,bm_id,"long_name",
     $            "Normal Field") )
      CALL check( nf90_def_var(mncid,"b_n_x",nf90_double,
     $            (/m_id,i_id/),bxm_id) )
      CALL check( nf90_put_att(mncid,bxm_id,"units","Tesla") )
      CALL check( nf90_put_att(mncid,bxm_id,"long_name",
     $            "Externally Applied Normal Field") )
      CALL check( nf90_def_var(mncid,"xi_n",nf90_double,
     $            (/m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(mncid,xm_id,"units","m") )
      CALL check( nf90_put_att(mncid,xm_id,"long_name",
     $            "Normal Displacement") )
      CALL check( nf90_def_var(mncid,"xi_n_x",nf90_double,
     $            (/m_id,i_id/),xxm_id) )
      CALL check( nf90_put_att(mncid,xxm_id,"units","m") )
      CALL check( nf90_put_att(mncid,xxm_id,"long_name",
     $            "Externally Applied Normal Displacement") )
      CALL check( nf90_def_var(mncid,"J_out",nf90_double,
     $            (/m_id,modid,i_id/),jo_id) )
      CALL check( nf90_put_att(mncid,jo_id,"long_name",
     $            "Transform to jac_out, tmag_out, jsurf_out (beta)") )
      CALL check( nf90_enddef(mncid) )
      templ = MATMUL(coordmat,binmn)
      CALL check( nf90_put_var(mncid,bxm_id,RESHAPE((/REAL(binmn),
     $          AIMAG(binmn)/),(/mpert,2/))) )
      templ = MATMUL(coordmat,boutmn)
      CALL check( nf90_put_var(mncid,bm_id,RESHAPE((/REAL(boutmn),
     $          AIMAG(boutmn)/),(/mpert,2/))) )
      templ = MATMUL(coordmat,xinmn)
      CALL check( nf90_put_var(mncid,xm_id,RESHAPE((/REAL(xinmn),
     $          AIMAG(xinmn)/),(/mpert,2/))) )
      templ = MATMUL(coordmat,xoutmn)
      CALL check( nf90_put_var(mncid,bm_id,RESHAPE((/REAL(xoutmn),
     $          AIMAG(xoutmn)/),(/mpert,2/))) )
      CALL check( nf90_put_var(mncid,mo_id,lmfac) )
      tempml = TRANSPOSE(coordmat)
      CALL check( nf90_put_var(mncid,jo_id,RESHAPE((/REAL(tempml),
     $          AIMAG(tempml)/),(/mpert,lmpert,2/))) )
      CALL check( nf90_close(mncid) )


      IF (fun_flag) THEN
         CALL ipeq_bcoords(psilim,binmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,boutmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,xinmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL ipeq_bcoords(psilim,xoutmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)

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
     $           (cos(eta)*w(1,1)-sin(eta)*w(1,2))/delpsi(itheta)
            zvecs(itheta)=
     $           (sin(eta)*w(1,1)+cos(eta)*w(1,2))/delpsi(itheta)
         ENDDO

         CALL ascii_open(out_unit,"ipec_control_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_CONTROL_FUN: "//
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
     $        REAL(xinfun(itheta)),-helicity*AIMAG(xinfun(itheta)),
     $        REAL(xoutfun(itheta)),-helicity*AIMAG(xoutfun(itheta)),
     $        REAL(binfun(itheta)),-helicity*AIMAG(binfun(itheta)),
     $        REAL(boutfun(itheta)),-helicity*AIMAG(boutfun(itheta))
         ENDDO
         CALL ascii_close(out_unit)
      
         IF(debug_flag) PRINT *,"Opening "//TRIM(mncfile)
         CALL check( nf90_open(mncfile,nf90_write,mncid) )
         IF(debug_flag) PRINT *,"  Inquiring about dimensions"
         CALL check( nf90_inq_dimid(mncid,"i",i_id) )
         CALL check( nf90_inq_dimid(mncid,"theta",t_id) )
         IF(debug_flag) PRINT *,"  Defining variables"
         CALL check( nf90_redef(mncid))
         CALL check( nf90_def_var(mncid, "xi_n_fun", nf90_double,
     $                    (/t_id,i_id/),x_id) )
         CALL check( nf90_put_att(mncid,x_id,"long_name",
     $               "Displacement") )
         CALL check( nf90_put_att(mncid,x_id,"units","m") )
         CALL check( nf90_def_var(mncid, "xi_n_x_fun", nf90_double,
     $                    (/t_id,i_id/),xx_id) )
         CALL check( nf90_put_att(mncid,xx_id,"long_name",
     $               "Externally Applied Displacement") )
         CALL check( nf90_put_att(mncid,xx_id,"units","m") )
         CALL check( nf90_def_var(mncid, "b_n_fun", nf90_double,
     $                    (/t_id,i_id/),b_id) )
         CALL check( nf90_put_att(mncid,b_id,"long_name",
     $               "Field") )
         CALL check( nf90_put_att(mncid,x_id,"units","Tesla") )
         CALL check( nf90_def_var(mncid, "b_n_x_fun", nf90_double,
     $                    (/t_id,i_id/),bx_id) )
         CALL check( nf90_put_att(mncid,bx_id,"long_name",
     $               "Externally Applied Field") )
         CALL check( nf90_put_att(mncid,bx_id,"units","Tesla") )
         CALL check( nf90_def_var(mncid, "R", nf90_double,t_id,r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $               "Major Radius") )
         CALL check( nf90_put_att(mncid,r_id,"units","m") )
         CALL check( nf90_def_var(mncid, "z", nf90_double,t_id,z_id) )
         CALL check( nf90_put_att(mncid,z_id,"long_name",
     $               "Vertical Position") )
         CALL check( nf90_put_att(mncid,z_id,"units","m") )
         CALL check( nf90_def_var(mncid, "R_n", nf90_double,t_id,rn_id))
         CALL check( nf90_put_att(mncid,rn_id,"long_name",
     $               "Major radius component of normal unit vector") )
         CALL check( nf90_put_att(mncid,rn_id,"units","m") )
         CALL check( nf90_def_var(mncid, "z_n", nf90_double,t_id,zn_id))
         CALL check( nf90_put_att(mncid,zn_id,"long_name",
     $               "Vertical component of normal unit vector") )
         CALL check( nf90_put_att(mncid,zn_id,"units","m") )
         CALL check( nf90_def_var(mncid, "dphi", nf90_double,t_id,p_id))
         CALL check( nf90_put_att(mncid,p_id,"long_name",
     $               "Toroidal - Magnetic Angle") )
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
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_control
c-----------------------------------------------------------------------
c     subprogram 4. ipout_singfld.
c     compute current and field on rational surfaces.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_singfld(egnum,xspmn,spot,
     $     rout,bpout,bout,rcout,tout,svd_flag)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      REAL(r8), INTENT(IN) :: spot
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      LOGICAL, INTENT(IN) :: svd_flag

      INTEGER :: i,itheta,ising
      REAL(r8) :: respsi,lpsi,rpsi,shear,hdist,sbnosurf
      COMPLEX(r8) :: lbwp1mn,rbwp1mn

      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: area,j_c,aq,asingflx
      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,sqreqb,jcfun
      COMPLEX(r8), DIMENSION(mpert) :: fkaxmn

      REAL(r8), DIMENSION(msing) :: island_hwidth,chirikov,
     $     novf,novs,novi
      COMPLEX(r8), DIMENSION(msing) :: delta,delcur,singcur,
     $     ovf,ovs,ovi
      COMPLEX(r8), DIMENSION(mpert,msing) :: singflx_mn
c-----------------------------------------------------------------------
c     solve equation from the given poloidal perturbation.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing total resonant fields"
      CALL ipeq_alloc
      CALL idcon_build(egnum,xspmn)
      IF (vsbrzphi_flag) ALLOCATE(singbno_mn(mpert,msing))
c-----------------------------------------------------------------------
c     evaluate delta and singular currents.
c     delta is delta*chi1*sq%f(4) and j_c is j_c/(chi1*sq%f(4))
c-----------------------------------------------------------------------
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

         lpsi=respsi-spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(lpsi)
         lbwp1mn=bwp1_mn(resnum(ising))
         
         rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1))
         CALL ipeq_sol(rpsi)
         rbwp1mn=bwp1_mn(resnum(ising))

         delta(ising)=rbwp1mn-lbwp1mn
         delcur(ising)=j_c(ising)*delta(ising)*ifac/
     $        (twopi*mfac(resnum(ising)))
         singcur(ising)=-delcur(ising)/ifac

         fkaxmn=0
         fkaxmn(resnum(ising))=singcur(ising)/(twopi*nn)
         
         ALLOCATE(fsurf_indev(mpert),fsurf_indmats(mpert,mpert))         
         CALL ipvacuum_flxsurf(respsi)
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
            CALL ipeq_weight(respsi,singbno_mn(:,ising),mfac,mpert,0)
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
      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_singfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGFLD: "//
     $     "Resonant fields, singular currents, and islands"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout 
      WRITE(out_unit,'(1x,a12,es17.8e3)')"sweet-spot =",spot
      WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a6,7(1x,a16))')"q","psi",
     $     "real(singflx)","imag(singflx)",
     $     "real(singcur)","imag(singcur)",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(1x,f6.3,7(es17.8e3))')
     $        singtype(ising)%q,singtype(ising)%psifac,
     $        REAL(singflx_mn(resnum(ising),ising)),
     $        AIMAG(singflx_mn(resnum(ising),ising)),
     $        REAL(singcur(ising)),AIMAG(singcur(ising)),
     $        island_hwidth(ising),chirikov(ising)
      ENDDO
      WRITE(out_unit,*)

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

      IF (svd_flag) THEN
         sbnosurf=SQRT(ABS(DOT_PRODUCT(sbno_fun(1:mthsurf),
     $        sbno_fun(1:mthsurf)))/mthsurf/2.0)
         sbno_mn = MATMUL(fldflxmat,sbno_mn)
         DO ising=1,msing
            ovf(ising)=DOT_PRODUCT(w1v(:,ising),sbno_mn(:))/SQRT(2.0)
            ovs(ising)=DOT_PRODUCT(w2v(:,ising),sbno_mn(:))/SQRT(2.0)
            ovi(ising)=DOT_PRODUCT(w3v(:,ising),sbno_mn(:))/SQRT(2.0)
         ENDDO
         DO ising=1,msing
            novf(ising)=ABS(ovf(ising))/sbnosurf*1e2
            novs(ising)=ABS(ovs(ising))/sbnosurf*1e2
            novi(ising)=ABS(ovi(ising))/sbnosurf*1e2
         ENDDO

         WRITE(out_unit,*)"Overlap fields, overlap singular "//
     $        "currents, and overlap islands"   
         WRITE(out_unit,*)        
         WRITE(out_unit,'(1x,a6,9(1x,a16))')
     $        "mode","real(ovf)","imag(ovf)","overlap(%)",
     $        "real(ovs)","imag(ovs)","overlap(%)",
     $        "real(ovi)","imag(ovi)","overlap(%)"
         DO ising=1,msing
            WRITE(out_unit,'(1x,I6,9(es17.8e3))')ising,
     $           REAL(ovf(ising)),AIMAG(ovf(ising)),novf(ising),
     $           REAL(ovs(ising)),AIMAG(ovs(ising)),novs(ising),
     $           REAL(ovi(ising)),AIMAG(ovi(ising)),novi(ising)
         ENDDO
         WRITE(out_unit,*)

         ! log svd response with harvest
         ierr=set_harvest_payload_dbl_array(hlog,'overlap_percent'//nul,
     $        novf,msing)
      ENDIF
        
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_singfld
c-----------------------------------------------------------------------
c     subprogram 5. ipout_vsingfld.
c     compute resonant field by coils.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vsingfld(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout

      INTEGER :: ipert,ising,i
      REAL(r8) :: hdist,shear,area
      INTEGER, DIMENSION(msing) :: resnum
      REAL(r8), DIMENSION(msing) :: visland_hwidth,vchirikov
      COMPLEX(r8), DIMENSION(:), POINTER :: vcmn

      COMPLEX(r8), DIMENSION(msing) :: vflxmn
      REAL(r8), DIMENSION(0:mthsurf) :: unitfun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
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
      CALL ascii_open(out_unit,"ipec_vsingfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_SINGFLD: "//
     $     "Resonant fields and islands from coils"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8,1x,a12,I2)')
     $     "jac_out = ",jac_out,"tmag_out =",tout 
      WRITE(out_unit,'(1x,a14)')"sweet-spot = 0"
      WRITE(out_unit,'(1x,a12,1x,I4)')"msing =",msing
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a6,5(1x,a16))')"q","psi",
     $     "real(singflx)","imag(singflx)",
     $     "islandhwidth","chirikov"
      DO ising=1,msing
         WRITE(out_unit,'(1x,f6.3,5(es17.8e3))')
     $        singtype(ising)%q,singtype(ising)%psifac,
     $        REAL(vflxmn(ising)),AIMAG(vflxmn(ising)),
     $        visland_hwidth(ising),vchirikov(ising)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_vsingfld
c-----------------------------------------------------------------------
c     subprogram 6. ipout_dw.
c     restore energy and torque profiles from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_dw(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep
      TYPE(cspline_type) :: dwk
c-----------------------------------------------------------------------
c     build and spline energy and torque profiles.
c-----------------------------------------------------------------------
      CALL idcon_build(egnum,xspmn)
      CALL cspline_alloc(dwk,mstep,1)
      dwk%xs=psifac
      DO istep=0,mstep
         dwk%fs(istep,1)=2*nn*ifac*
     $        SUM(CONJG(u1%fs(istep,:))*u2%fs(istep,:))/(2.0*mu0)
      ENDDO
      CALL cspline_fit(dwk,"extrap")

      WRITE(*,*)"Restoring energy and torque profiles"

      CALL ascii_open(out_unit,"ipec_dw_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_dw: "//
     $     "Energy and torque profiles by self-consistent solutions."
      WRITE(out_unit,*)

      WRITE(out_unit,'(6(1x,a16))')"psi","T_phi","2ndeltaW",
     $     "int(T_phi)","int(2ndeltaW)","dv/dpsi_n"
      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         CALL cspline_eval(dwk,psifac(istep),1)
         CALL spline_eval(sq,psifac(istep),0)
         WRITE(out_unit,'(6(1x,es16.8))')
     $        psifac(istep),REAL(dwk%f1(1)),AIMAG(dwk%f1(1)),
     $        REAL(dwk%fs(istep,1)),AIMAG(dwk%fs(istep,1)),sq%f(3)
      ENDDO
      CALL ascii_close(out_unit)
      CALL cspline_dealloc(dwk)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_dw
c-----------------------------------------------------------------------
c     subprogram 7. ipout_dw_matrix.
c     restore energy and torque response matrix from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_dw_matrix
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: ipert,istep,itheta,i,j,mval,iindex,lwork
      REAL(r8) :: jarea,ileft
      COMPLEX(r8) :: t1,t2
      LOGICAL :: diagnose1=.FALSE.,diagnose2=.TRUE.
      LOGICAL :: kstar_coil_optimization=.FALSE.

      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8), DIMENSION(mpert) :: teigen
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      REAL(r8), DIMENSION(2*mpert) :: rwork2
      REAL(r8), DIMENSION(0:mthsurf) :: units
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work2
      COMPLEX(r8), DIMENSION(mpert*mpert) :: work3
      COMPLEX(r8), DIMENSION(mpert) :: fldflxmn,tvec1,tvec2,reigen
      COMPLEX(r8), DIMENSION(mstep) :: tprof
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,fldflxmat,vr,vl
      COMPLEX(r8), DIMENSION(mstep,mpert,mpert) :: bsurfmat,dwks,dwk,
     $     gind,gindp,gres,gresp
      COMPLEX(r8), DIMENSION(5,mpert) :: tvecs
      COMPLEX(r8), DIMENSION(9,mpert,mpert) :: tgresp

      REAL(r8), DIMENSION(10,48) :: icoilcur
      CHARACTER(24), DIMENSION(10) :: icoilname
      COMPLEX(r8), DIMENSION(:), POINTER :: coilmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: coilmns,coil1,
     $     coil2,coil3,coil4
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: gcoil

      INTEGER, DIMENSION(3) :: cipiv
      REAL(r8), DIMENSION(3) :: ceigen
      REAL(r8), DIMENSION(3*3-2) :: crwork
      REAL(r8), DIMENSION(2*3) :: crwork2
      COMPLEX(r8), DIMENSION(2*3-1) :: cwork
      COMPLEX(r8), DIMENSION(2*3+1) :: cwork2
      COMPLEX(r8), DIMENSION(3*3) :: cwork3
      COMPLEX(r8), DIMENSION(3) :: creigen
      COMPLEX(r8), DIMENSION(3,3) :: cvr,cvl

      TYPE(cspline_type) :: optorq
c-----------------------------------------------------------------------
c     call solutions.
c-----------------------------------------------------------------------
      units = (/(1.0,itheta=0,mthsurf)/)
      jarea = issurfint(units,mthsurf,psilim,0,0)
      DO ipert=1,mpert
         fldflxmn=0
         fldflxmn(ipert)=1.0
         CALL ipeq_weight(psilim,fldflxmn,mfac,mpert,2)
         fldflxmat(:,ipert)=fldflxmn*sqrt(jarea)
      ENDDO
      WRITE(*,*)
     $     "Call solutions for general response matrix functions"
      edge_flag=.TRUE.
      DO ipert=1,mpert
         edge_mn=0
         edge_mn(ipert)=1.0
         CALL idcon_build(0,edge_mn)
         DO istep=0,mstep
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
         temp1=MATMUL(gind(istep,:,:),fldflxmat)
         temp2=CONJG(TRANSPOSE(fldflxmat))
         gindp(istep,:,:)=MATMUL(temp2,temp1)
c-----------------------------------------------------------------------
c     make coordinate-independent matrix (normalized by phix power).
c-----------------------------------------------------------------------
         temp1=MATMUL(gres(istep,:,:),fldflxmat)
         temp2=CONJG(TRANSPOSE(fldflxmat))
         gresp(istep,:,:)=MATMUL(temp2,temp1)
      ENDDO
c-----------------------------------------------------------------------
c     test by field applications.
c-----------------------------------------------------------------------
      IF (diagnose1) THEN
         mval=1
         tvec1=0
         tvec1(mval-mlow+1)=1e-3

         CALL cspline_alloc(optorq,mstep,1)
         optorq%xs=psifac

         DO istep=1,mstep
            temp1=MATMUL(gresp(istep,:,:),fldflxmat)
            temp2=MATMUL(CONJG(TRANSPOSE(fldflxmat)),temp1)
            temp1=(temp2+CONJG(TRANSPOSE(temp2)))/2.0
            tvec2=MATMUL(temp1,tvec1)
            tprof(istep)=DOT_PRODUCT(tvec1,tvec2)/jarea**2
            optorq%fs(istep,1)=tprof(istep)
         ENDDO

         optorq%fs(0,1)=tprof(1)-(tprof(2)-tprof(1))/
     $        (psifac(2)-psifac(1))*(psifac(1)-psifac(0))
         CALL cspline_fit(optorq,"extrap")
         CALL ascii_open(out_unit,"ipec_dw_diag_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_dw_diag: "//
     $        "Energy and torque profiles by self-consistent solutions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(6(1x,a16))')"psi","T_phi","2ndeltaW",
     $        "int(T_phi)","int(2ndeltaW)","dv/dpsi_n"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            CALL cspline_eval(optorq,psifac(istep),1)
            CALL spline_eval(sq,psifac(istep),0)
            WRITE(out_unit,'(6(1x,es16.8))')
     $           psifac(istep),REAL(optorq%f1(1)),AIMAG(optorq%f1(1)),
     $           REAL(optorq%f(1)),AIMAG(optorq%f(1)),sq%f(3)
         ENDDO
         CALL ascii_close(out_unit)
         CALL cspline_dealloc(optorq)
      ENDIF
c-----------------------------------------------------------------------
c     write general response matrix functions.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_dw_matrix_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_dw_matrix: "//
     $     "Self-consistent response matrix functions."
      WRITE(out_unit,*)

      WRITE(out_unit,'(1x,a16,2(1x,a4),10(1x,a16))')"psi","i","j",
     $     "Re(T_x)","Im(T_x)","Re(T_f)","Im(T_f)",
     $     "Re(T_ef)","Im(T_ef)","Re(T_fp)","Im(T_fp)",
     $     "Re(T_efp)","Im(T_efp)"
      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         DO i=1,mpert
            DO j=1,mpert
               WRITE(out_unit,'(1x,es16.8,2(1x,I4),10(1x,es16.8))')
     $              psifac(istep),i,j,
     $              REAL(dwk(istep,i,j)),AIMAG(dwk(istep,i,j)),
     $              REAL(gind(istep,i,j)),AIMAG(gind(istep,i,j)),
     $              REAL(gres(istep,i,j)),AIMAG(gres(istep,i,j)),
     $              REAL(gindp(istep,i,j)),AIMAG(gindp(istep,i,j)),
     $              REAL(gresp(istep,i,j)),AIMAG(gresp(istep,i,j))
            ENDDO
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     calculate maximum torque eigenvalues.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_dw_eigen_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_dw_eigen: "//
     $     "Eigenvalue profiles by self-consistent solutions."
      WRITE(out_unit,*)

      WRITE(out_unit,'(3(1x,a16))')"psi","MAX(teigen)","MIN(teigen)"
      lwork=2*mpert-1
      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         temp1=(gresp(istep,:,:)+CONJG(TRANSPOSE(gresp(istep,:,:))))/2.0
         CALL zheev('V','U',mpert,temp1,mpert,teigen,
     $        work,lwork,rwork,info)
         WRITE(out_unit,'(3(1x,es16.8))')
     $        psifac(istep),MAXVAL(teigen),MINVAL(teigen)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     diagnose eigenvector.
c-----------------------------------------------------------------------
      IF (diagnose1) THEN
         tvec1=temp1(:,1)
         DO istep=1,mstep
            temp2=gresp(istep,:,:)
            temp1=(temp2+CONJG(TRANSPOSE(temp2)))/2.0
            tvec2=MATMUL(temp1,tvec1)
            tprof(istep)=DOT_PRODUCT(tvec1,tvec2)
         ENDDO
         CALL ascii_open(out_unit,"ipec_dw_minimum_torque_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_dw_diag: "//
     $        "Energy and torque profiles by self-consistent solutions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(3(1x,a16))')"psi","int(T_phi)","int(2ndeltaW)"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            WRITE(out_unit,'(3(1x,es16.8))')psifac(istep),
     $           REAL(tprof(istep)),AIMAG(tprof(istep))
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     calculate maximum torque-ratio TT^{-1}_b eigenvalues.
c-----------------------------------------------------------------------
      temp1=(gresp(mstep,:,:)+CONJG(TRANSPOSE(gresp(mstep,:,:))))/2.0
      CALL zhetrf('L',mpert,temp1,mpert,ipiv,work3,mpert*mpert,info)

      CALL ascii_open(out_unit,"ipec_dw_ratio_eigen_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_dw_ratio_eigen: "//
     $     "Eigenvalue profiles by self-consistent solutions."
      WRITE(out_unit,*)

      WRITE(out_unit,'(3(1x,a16))')"psi","MAX(teigen)","MIN(teigen)"
      lwork=2*mpert+1
      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         temp2=(gresp(istep,:,:)+CONJG(TRANSPOSE(gresp(istep,:,:))))/2.0
         DO i=1,9
            IF (ABS(psifac(istep)-0.1*i)<1e-3) THEN
               tgresp(i,:,:)=temp2
            ENDIF
         ENDDO
         CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         CALL zgeev('V','V',mpert,temp2,mpert,reigen,
     $        vl,mpert,vr,mpert,work2,lwork,rwork2,info)
         WRITE(out_unit,'(3(1x,es16.8))')
     $        psifac(istep),MAXVAL(ABS(reigen)),MINVAL(ABS(reigen))
         IF (ABS(psifac(istep)-0.5)<1e-3) THEN
            tvecs(1,:)=vr(:,1)
         ENDIF
      ENDDO
      DO i=1,4
         temp2=tgresp(i+5,:,:)-tgresp(i,:,:)
         CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         CALL zgeev('V','V',mpert,temp2,mpert,reigen,
     $        vl,mpert,vr,mpert,work2,lwork,rwork2,info)
         tvecs(i+1,:)=vr(:,1)
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     optimized torque profile (0-0.5,0.1-0.6,0.2-0.7..0.4-0.9)
c-----------------------------------------------------------------------
      IF (diagnose2) THEN
         CALL cspline_alloc(optorq,mstep,5)
         optorq%xs=psifac
         DO i=1,5
            tvec1=tvecs(i,:)
            DO istep=1,mstep
               temp2=gresp(istep,:,:)
               temp1=(temp2+CONJG(TRANSPOSE(temp2)))/2.0
               tvec2=MATMUL(temp1,tvec1)
               tprof(istep)=DOT_PRODUCT(tvec1,tvec2)
               optorq%fs(istep,i)=tprof(istep)
            ENDDO
            optorq%fs(0,i)=tprof(1)-(tprof(2)-tprof(1))/
     $           (psifac(2)-psifac(1))*(psifac(1)-psifac(0))
         ENDDO
         CALL cspline_fit(optorq,"extrap")
         CALL ascii_open(out_unit,"ipec_dw_optimized_torque_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_dw_diag: "//
     $        "Energy and torque profiles by self-consistent solutions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(7(1x,a16))')"psi","T_phi05",
     $        "T_phi16","T_phi27","T_phi38","T_phi49","dv/dpsi_n"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            CALL cspline_eval(optorq,psifac(istep),1)
            CALL spline_eval(sq,psifac(istep),0)
            WRITE(out_unit,'(7(1x,es16.8))')
     $           psifac(istep),REAL(optorq%f1(1)),REAL(optorq%f1(2)),
     $           REAL(optorq%f1(3)),REAL(optorq%f1(4)),
     $           REAL(optorq%f1(5)),sq%f(3)
         ENDDO
         CALL ascii_close(out_unit)
         CALL cspline_dealloc(optorq)
      ENDIF
c-----------------------------------------------------------------------
c     construct coil-torque response matrix.
c----------------------------------------------------------------------
      IF (kstar_coil_optimization) THEN
         ALLOCATE(coilmns(mpert,3),coil1(mpert,3))
         ALLOCATE(coil2(3,mpert),coil3(3,3),coil4(3,3))
         icoilname='fecu'
         icoilcur(1,1)=1e3
         icoilcur(1,2)=0
         icoilcur(1,3)=-1e3
         icoilcur(1,4)=0
         CALL coil_read(idconfile,1,icoilname,icoilcur)
         ALLOCATE(coilmn(cmpert))
         CALL field_bs_psi(psilim,coilmn,1)
         DO i=1,cmpert
            IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
               coilmns(cmlow-mlow+i,1)=coilmn(i)
            ENDIF
         ENDDO
         DEALLOCATE(coilmn)
         CALL coil_dealloc
         icoilname='fecm'
         icoilcur(1,1)=1e3
         icoilcur(1,2)=0
         icoilcur(1,3)=-1e3
         icoilcur(1,4)=0
         CALL coil_read(idconfile,1,icoilname,icoilcur)
         ALLOCATE(coilmn(cmpert))
         CALL field_bs_psi(psilim,coilmn,1)
         DO i=1,cmpert
            IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
               coilmns(cmlow-mlow+i,2)=coilmn(i)
            ENDIF
         ENDDO
         DEALLOCATE(coilmn)
         CALL coil_dealloc
         icoilname='fecl'
         icoilcur(1,1)=1e3
         icoilcur(1,2)=0
         icoilcur(1,3)=-1e3
         icoilcur(1,4)=0
         CALL coil_read(idconfile,1,icoilname,icoilcur)
         ALLOCATE(coilmn(cmpert))
         CALL field_bs_psi(psilim,coilmn,1)
         DO i=1,cmpert
            IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
               coilmns(cmlow-mlow+i,3)=coilmn(i)
            ENDIF
         ENDDO
         DEALLOCATE(coilmn)
         CALL coil_dealloc

         ALLOCATE(gcoil(mstep,3,3))

         WRITE(*,*)"Build coil response matrix functions"
         DO istep=1,mstep
            iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
            ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
            IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $           WRITE(*,'(1x,a9,i3,a26)')
     $           "volume = ",iindex,"% coil matrix calculations"

            coil1=MATMUL(gres(istep,:,:),coilmns(:,:))
            coil2=CONJG(TRANSPOSE(coilmns(:,:)))
            gcoil(istep,:,:)=MATMUL(coil2,coil1)
         ENDDO
c-----------------------------------------------------------------------
c     write coil response matrix functions.
c-----------------------------------------------------------------------
         CALL ascii_open(out_unit,"ipec_coil_matrix_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_coil_matrix: "//
     $        "Self-consistent coil response matrix functions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(1x,a16,2(1x,a4),2(1x,a16))')"psi","i","j",
     $        "Re(T_x)","Im(T_x)"
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            DO i=1,3
               DO j=1,3
                  WRITE(out_unit,'(1x,es16.8,2(1x,I4),2(1x,es16.8))')
     $                 psifac(istep),i,j,
     $                 REAL(gcoil(istep,i,j)),AIMAG(gcoil(istep,i,j))
               ENDDO
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     calculate maximum coil eigenvalues.
c-----------------------------------------------------------------------
         CALL ascii_open(out_unit,"ipec_coil_eigen_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_coil_eigen: "//
     $        "Coil eigenvalue profiles by self-consistent solutions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(15(1x,a16))')"psi",
     $        "MAX(ceigen)","MIN(ceigen)",
     $        "Real(MIN1)","Imag(MIN1)","Real(MIN2)","Imag(MIN2)",
     $        "Real(MIN3)","Imag(MIN3)","Real(MAX1)","Imag(MAX1)",
     $        "Real(MAX2)","Imag(MAX2)","Real(MAX3)","Imag(MAX3)"
         lwork=2*3-1
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            coil3=(gcoil(istep,:,:)+
     $           CONJG(TRANSPOSE(gcoil(istep,:,:))))/2.0
            CALL zheev('V','U',3,coil3,3,ceigen,cwork,lwork,crwork,info)
            WRITE(out_unit,'(15(1x,es16.8))')
     $           psifac(istep),MAXVAL(ceigen),MINVAL(ceigen),
     $           REAL(coil3(1,1)),AIMAG(coil3(1,1)),
     $           REAL(coil3(2,1)),AIMAG(coil3(2,1)),
     $           REAL(coil3(3,1)),AIMAG(coil3(3,1)),
     $           REAL(coil3(1,3)),AIMAG(coil3(1,3)),
     $           REAL(coil3(2,3)),AIMAG(coil3(2,3)),
     $           REAL(coil3(3,3)),AIMAG(coil3(3,3))
         ENDDO
         CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     calculate maximum torque-ratio TT^{-1}_b coil eigenvalues.
c-----------------------------------------------------------------------
         coil3=(gcoil(mstep,:,:)+CONJG(TRANSPOSE(gcoil(mstep,:,:))))/2.0
         CALL zhetrf('L',3,coil3,3,cipiv,cwork3,3*3,info)

         CALL ascii_open(out_unit,"ipec_coil_ratio_eigen_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_coil_ratio_eigen: "//
     $        "Eigenvalue profiles by self-consistent solutions."
         WRITE(out_unit,*)

         WRITE(out_unit,'(15(1x,a16))')"psi",
     $        "MAX(ceigen)","MIN(ceigen)",
     $        "Real(MIN1)","Imag(MIN1)","Real(MIN2)","Imag(MIN2)",
     $        "Real(MIN3)","Imag(MIN3)","Real(MAX1)","Imag(MAX1)",
     $        "Real(MAX2)","Imag(MAX2)","Real(MAX3)","Imag(MAX3)"

         lwork=2*3+1
         DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
            coil4=(gcoil(istep,:,:)+
     $           CONJG(TRANSPOSE(gcoil(istep,:,:))))/2.0
            CALL zhetrs('L',3,3,coil3,3,cipiv,coil4,3,info)
            CALL zgeev('V','V',3,coil4,3,creigen,
     $           cvl,3,cvr,3,cwork2,lwork,crwork2,info)

            WRITE(out_unit,'(15(1x,es16.8))')psifac(istep),
     $           MAXVAL(ABS(creigen)),MINVAL(ABS(creigen)),
     $           REAL(cvr(1,1)),AIMAG(cvr(1,1)),
     $           REAL(cvr(2,1)),AIMAG(cvr(2,1)),
     $           REAL(cvr(3,1)),AIMAG(cvr(3,1)),
     $           REAL(cvr(1,3)),AIMAG(cvr(1,3)),
     $           REAL(cvr(2,3)),AIMAG(cvr(2,3)),
     $           REAL(cvr(3,3)),AIMAG(cvr(3,3))
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipout_dw_matrix
c-----------------------------------------------------------------------
c     subprogram 6. ipout_pmodb.
c     compute perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_pmodb(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,itheta,iindex,cstep
      REAL(r8) :: ileft,jac,psi

      INTEGER :: p_id,t_id,i_id,m_id,r_id,z_id,b_id,bme_id,be_id,
     $   bml_id,bl_id,xm_id,x_id,km_id,k_id,rzstat

      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,equilbfun
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,
     $     divxprp_mn,curv_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmz_fun,xvt_fun,xvz_fun,xmt_fun,xsp1_fun
      COMPLEX(r8), DIMENSION(mstep,mpert) :: eulbparmns,lagbparmns,
     $     divxprpmns,curvmns,xmp1mns,xspmns,xmsmns,xmtmns,xmzmns
      COMPLEX(r8), DIMENSION(mstep,lmpert) :: eulbparmout,lagbparmout,
     $     divxprpmout,curvmout
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) ::eulbparfun,lagbparfun,
     $     divxprpfun,curvfun,eulbparfout,lagbparfout,divxprpfout,
     $     curvfout

      REAL(r8), DIMENSION(:), POINTER :: psis,ches,chex
      COMPLEX(r8), DIMENSION(:,:), POINTER ::chea,chelagbmns,chelagbmout
      COMPLEX(r8), DIMENSION(mpert) :: chelagb_mn

      TYPE(cspline_type) :: cspl,chespl
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing |b| and del.x_prp"
      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta

      CALL idcon_build(egnum,xspmn)
      CALL ipeq_alloc

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
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))         
         CALL ipeq_cova(psifac(istep))
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
         CALL ipeq_bcoordsout(eulbparmout(i,:),eulbpar_mn,psi,tout,0)
         CALL ipeq_bcoordsout(lagbparmout(i,:),lagbpar_mn,psi,tout,0)
         CALL ipeq_bcoordsout(divxprpmout(i,:),divxprp_mn,psi,tout,0)
         CALL ipeq_bcoordsout(curvmout(i,:)   ,curv_mn,psi,tout,0)
         eulbparfout(istep,:)=eulbparfun(istep,:)*EXP(ifac*nn*dphi)
         lagbparfout(istep,:)=lagbparfun(istep,:)*EXP(ifac*nn*dphi)
         divxprpfout(istep,:)=divxprpfun(istep,:)*EXP(ifac*nn*dphi)
         curvfout(istep,:)=curvfun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      CALL cspline_dealloc(cspl)

      ! append to netcdf file
      IF(debug_flag) PRINT *,"Opening "//TRIM(fncfile)
      CALL check( nf90_open(fncfile,nf90_write,fncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(fncid,"i",i_id) )
      CALL check( nf90_inq_dimid(fncid,"m",m_id) )
      CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
      CALL check( nf90_inq_dimid(fncid,"theta",t_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(fncid))
      rzstat = nf90_inq_varid(fncid, "R", r_id) ! check if R,z already stored
      IF(rzstat/=nf90_noerr)THEN
         CALL check( nf90_def_var(fncid, "R", nf90_double,
     $               (/p_id,t_id/),r_id) )
         CALL check( nf90_put_att(fncid,r_id,"long_name",
     $               "Major Radius") )
         CALL check( nf90_put_att(fncid,r_id,"units","m") )
         CALL check( nf90_def_var(fncid, "z", nf90_double,
     $               (/p_id,t_id/),z_id) )
         CALL check( nf90_put_att(fncid,z_id,"long_name",
     $               "Vertical Position") )
         CALL check( nf90_put_att(fncid,z_id,"units","m") )
      ENDIF
      CALL check( nf90_def_var(fncid, "b_eul_fun", nf90_double,
     $            (/p_id,t_id,i_id/),be_id) )
      CALL check( nf90_put_att(fncid,be_id,"long_name",
     $            "Eulerian Perturbed Field") )
      CALL check( nf90_put_att(fncid,be_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_lag_fun", nf90_double,
     $            (/p_id,t_id,i_id/),bl_id) )
      CALL check( nf90_put_att(fncid,bl_id,"long_name",
     $            "Lagrangian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bl_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_eul", nf90_double,
     $            (/p_id,m_id,i_id/),bme_id) )
      CALL check( nf90_put_att(fncid,bme_id,"long_name",
     $            "Eulerian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bme_id,"units","Tesla") )
      CALL check( nf90_put_att(fncid,bme_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "b_lag", nf90_double,
     $            (/p_id,m_id,i_id/),bml_id) )
      CALL check( nf90_put_att(fncid,bml_id,"long_name",
     $            "Lagrangian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bml_id,"units","Tesla") )
      CALL check( nf90_put_att(fncid,bml_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "Bdivxi_perp_fun", nf90_double,
     $            (/p_id,t_id,i_id/),x_id) )
      CALL check( nf90_put_att(fncid,x_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_def_var(fncid, "Bdivxi_perp", nf90_double,
     $            (/p_id,m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_put_att(fncid,xm_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "Bkappaxi_perp_fun", nf90_double,
     $            (/p_id,t_id,i_id/),k_id) )
      CALL check( nf90_put_att(fncid,k_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_def_var(fncid, "Bkappaxi_perp", nf90_double,
     $            (/p_id,m_id,i_id/),km_id) )
      CALL check( nf90_put_att(fncid,km_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_put_att(fncid,km_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "B", nf90_double,
     $            (/p_id,t_id/),b_id) )
      CALL check( nf90_put_att(fncid,b_id,"long_name",
     $            "Equilibrium Field Strength") )
      CALL check( nf90_enddef(fncid) )
      IF(debug_flag) PRINT *,"  Writting variables"
      IF(rzstat/=nf90_noerr)THEN
         CALL check( nf90_put_var(fncid,r_id,rs) )
         CALL check( nf90_put_var(fncid,z_id,zs) )
      ENDIF
      CALL check( nf90_put_var(fncid,bme_id,RESHAPE((/REAL(eulbparmout),
     $             AIMAG(eulbparmout)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,bml_id,RESHAPE((/REAL(lagbparmout),
     $             AIMAG(lagbparmout)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,xm_id,RESHAPE((/REAL(-divxprpmout),
     $             AIMAG(-divxprpmout)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,km_id,RESHAPE((/REAL(-curvmout),
     $             AIMAG(-curvmout)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,be_id,RESHAPE((/REAL(eulbparfout),
     $            -helicity*AIMAG(eulbparfout)/),(/mstep,mthsurf,2/))) )
      CALL check( nf90_put_var(fncid,bl_id,RESHAPE((/REAL(lagbparfout),
     $            -helicity*AIMAG(lagbparfout)/),(/mstep,mthsurf,2/))) )
      CALL check( nf90_put_var(fncid,x_id,RESHAPE((/REAL(-divxprpfout),
     $            -helicity*AIMAG(-divxprpfout)/),(/mstep,mthsurf,2/))))
      CALL check( nf90_put_var(fncid,k_id,RESHAPE((/REAL(-curvfout),
     $            -helicity*AIMAG(-curvfout)/),(/mstep,mthsurf,2/))) )
      CALL check( nf90_put_var(fncid,b_id,equilbfun) )

      CALL check( nf90_close(fncid) )
      IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)


      CALL ascii_open(out_unit,"ipec_pmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_PMODB: "//
     $     "Components in perturbed mod b and del.x_prp"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     

      WRITE(out_unit,'(1x,a16,1x,a4,10(1x,a16))')"psi","m",
     $     "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $     "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $     "imag(Bkxprp)"
      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         DO ipert=1,lmpert
            WRITE(out_unit,'(es17.8e3,1x,I4,8(es17.8e3))')
     $           psifac(istep),lmfac(ipert),
     $           REAL(eulbparmout(istep,ipert)),
     $           AIMAG(eulbparmout(istep,ipert)),
     $           REAL(lagbparmout(istep,ipert)),
     $           AIMAG(lagbparmout(istep,ipert)),
     $           REAL(-divxprpmout(istep,ipert)),
     $           AIMAG(-divxprpmout(istep,ipert)),
     $           REAL(-curvmout(istep,ipert)),
     $           AIMAG(-curvmout(istep,ipert))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

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

         CALL ascii_open(out_unit,"ipec_pmodb_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_CHEBYSHEV: "//
     $        "Chebyshev coefficients for perturbed mod b"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nche =",nche,"mpert =",mpert
         WRITE(out_unit,*)     
         
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,a6,1x,I4)')"mfac =",mfac(ipert)
            WRITE(out_unit,*)
            WRITE(out_unit,'(1x,a6,2(1x,a16))')"chea",
     $           "real(chea)","imag(chea)"
            DO i=0,nche
               WRITE(out_unit,'(1x,I6,2(es17.8e3))')i,
     $              REAL(chea(ipert,i)),AIMAG(chea(ipert,i))
            ENDDO
            WRITE(out_unit,*)
         ENDDO
         CALL ascii_close(out_unit)

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

         CALL ascii_open(out_unit,"ipec_chepmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_CHEPMODB: "//
     $        "Reconstructed perturbed mod b by chebyshev"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "cstep =",cstep,"mpert =",mpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         
         WRITE(out_unit,'(1x,a16,1x,a4,2(1x,a16))')"psi","m",
     $        "real(lagb)","imag(lagb)"
         DO istep=1,cstep
            DO ipert=1,mpert
               WRITE(out_unit,'(es17.8e3,1x,I4,2(es17.8e3))')
     $              psis(istep),mfac(ipert),
     $              REAL(chelagbmns(istep,ipert)),
     $              AIMAG(chelagbmns(istep,ipert))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

         CALL cspline_dealloc(chespl)
      ENDIF

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_pmodb_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PMODB_FUN: "//
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
            CALL ipeq_sol(psifac(istep))
            CALL ipeq_contra(psifac(istep))         
            CALL ipeq_cova(psifac(istep))
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

      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_pmodb
c-----------------------------------------------------------------------
c     subprogram 7. ipout_xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbnormal(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: p_id,t_id,i_id,m_id,r_id,z_id,bm_id,b_id,xm_id,x_id,
     $   rv_id,zv_id,rzstat

      INTEGER :: i,istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax,area

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs,dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun

      COMPLEX(r8), DIMENSION(mstep,lmpert) :: xmns,ymns,
     $     xnomns,bnomns,bwpmns
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn
      COMPLEX(r8), DIMENSION(mstep,lmpert) :: pwpmns

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,psis,rvecs,zvecs
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss,
     $     xnofuns,bnofuns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b normal components"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))

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
            pwpmn=0
            DO i=1,mpert
               IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
                  pwpmn(mlow-lmlow+i)=bno_mn(i)
               ENDIF
            ENDDO
            CALL ipeq_bcoords(psifac(istep),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)
            pwpmns(istep,:)=pwpmn
         ENDIF            

         CALL ipeq_bcoordsout(xnomns(istep,:),xno_mn,psifac(istep),ji=0)
         CALL ipeq_bcoordsout(bnomns(istep,:),bno_mn,psifac(istep),ji=0)
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            CALL ipeq_bcoordsout(bwpmns(istep,:),bno_mn,psifac(istep),
     $                           ji=1)
         ELSE ! no need to re-weight bno_mn with expensive invfft and fft
            bwp_mn = bwp_mn/area
         ENDIF
         xnofuns(istep,:)=xnofuns(istep,:)*EXP(ifac*nn*dphi)
         bnofuns(istep,:)=bnofuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBNORMAL: "//
     $     "Normal components of displacement and field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $     "real(xno)","imag(xno)","real(bno)","imag(bno)",
     $     "real(bwp)","imag(bwp)"

      DO istep=1,mstep,MAX(1,(mstep*lmpert-1)/max_linesout+1)
         DO ipert=1,lmpert
            WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $           psifac(istep),qfac(istep),lmfac(ipert),
     $           REAL(xnomns(istep,ipert)),AIMAG(xnomns(istep,ipert)),
     $           REAL(bnomns(istep,ipert)),AIMAG(bnomns(istep,ipert)),
     $           REAL(bwpmns(istep,ipert)),AIMAG(bwpmns(istep,ipert))        
         ENDDO
      ENDDO

      IF (bwp_pest_flag) THEN
         CALL ascii_open(out_unit,"ipec_bnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BNORMAL_PEST: "//
     $        "Normal components of field in pest coordinates"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mstep =",mstep,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         
         DO istep=1,mstep,MAX(1,(mstep*lmpert-1)/max_linesout+1)
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $              psifac(istep),qfac(istep),lmfac(ipert),
     $              REAL(pwpmns(istep,ipert)),AIMAG(pwpmns(istep,ipert))        
            ENDDO
         ENDDO
      ENDIF

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbnormal_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBNORMAL_FUN: "//
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
         DO ipert=1,lmpert
            DO istep=1,mstep
               WRITE(bin_unit)REAL(psifac(istep),4),
     $              REAL(REAL(xnomns(istep,ipert)),4),
     $              REAL(AIMAG(xnomns(istep,ipert)),4),
     $              REAL(REAL(bnomns(istep,ipert)),4),
     $              REAL(AIMAG(bnomns(istep,ipert)),4),
     $              REAL(REAL(bwpmns(istep,ipert)),4),
     $              REAL(AIMAG(bwpmns(istep,ipert)),4)    
            ENDDO
            WRITE(bin_unit)
         ENDDO
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

      ! append to netcdf file
      IF(debug_flag) PRINT *,"Opening "//TRIM(fncfile)
      CALL check( nf90_open(fncfile,nf90_write,fncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(fncid,"i",i_id) )
      CALL check( nf90_inq_dimid(fncid,"m",m_id) )
      CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
      CALL check( nf90_inq_dimid(fncid,"theta",t_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(fncid))
      rzstat = nf90_inq_varid(fncid, "R", r_id) ! check if R,z already stored
      IF(rzstat/=nf90_noerr)THEN
         CALL check( nf90_def_var(fncid, "R", nf90_double,
     $               (/p_id,t_id/),r_id) )
         CALL check( nf90_put_att(fncid,r_id,"long_name",
     $               "Major Radius") )
         CALL check( nf90_put_att(fncid,r_id,"units","m") )
         CALL check( nf90_def_var(fncid, "z", nf90_double,
     $               (/p_id,t_id/),z_id) )
         CALL check( nf90_put_att(fncid,z_id,"long_name",
     $               "Vertical Position") )
         CALL check( nf90_put_att(fncid,z_id,"units","m") )
      ENDIF
      CALL check( nf90_def_var(fncid, "b_n_fun", nf90_double,
     $            (/p_id,t_id,i_id/),b_id) )
      CALL check( nf90_put_att(fncid,b_id,"long_name",
     $            "Perturbed Field Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,b_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_n", nf90_double,
     $            (/p_id,m_id,i_id/),bm_id) )
      CALL check( nf90_put_att(fncid,bm_id,"long_name",
     $            "Perturbed Field Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,bm_id,"units","Tesla") )
      CALL check( nf90_put_att(fncid,bm_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "xi_n_fun", nf90_double,
     $            (/p_id,t_id,i_id/),x_id) )
      CALL check( nf90_put_att(fncid,x_id,"long_name",
     $            "Displacement Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,x_id,"units","m") )
      CALL check( nf90_def_var(fncid, "xi_n", nf90_double,
     $            (/p_id,m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $            "Displacement Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,xm_id,"units","m") )
      CALL check( nf90_put_att(fncid,xm_id,"jacobian",jac_out) )
      CALL check( nf90_def_var(fncid, "R_n", nf90_double,
     $            (/p_id,t_id/),rv_id) )
      CALL check( nf90_put_att(fncid,rv_id,"long_name",
     $            "Radial unit normal: R_3D = R_sym+x_n*R_n") )
      CALL check( nf90_def_var(fncid, "z_n", nf90_double,
     $            (/p_id,t_id/),zv_id) )
      CALL check( nf90_put_att(fncid,zv_id,"long_name",
     $            "Vertical unit normal: z_3D = z_sym+x_n*z_n") )
      CALL check( nf90_enddef(fncid) )
      IF(debug_flag) PRINT *,"  Writting variables"
      IF(rzstat/=nf90_noerr)THEN
         CALL check( nf90_put_var(fncid,r_id,rs) )
         CALL check( nf90_put_var(fncid,z_id,zs) )
      ENDIF
      CALL check( nf90_put_var(fncid,xm_id,RESHAPE((/REAL(xnomns),
     $             AIMAG(xnomns)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,bm_id,RESHAPE((/REAL(bnomns),
     $             AIMAG(bnomns)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,x_id,RESHAPE((/REAL(xnofuns),
     $             -helicity*AIMAG(xnofuns)/),(/mstep,mthsurf+1,2/))) )
      CALL check( nf90_put_var(fncid,b_id,RESHAPE((/REAL(bnofuns),
     $             -helicity*AIMAG(bnofuns)/),(/mstep,mthsurf+1,2/))) )
      CALL check( nf90_put_var(fncid,rv_id,rvecs) )
      CALL check( nf90_put_var(fncid,zv_id,zvecs) )
      CALL check( nf90_close(fncid) )
      IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)

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
            psis(:,itheta)=psifac(:)
         ENDDO
         
         CALL ascii_open(out_unit,"ipec_xbnormal_flux_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBNROMAL_FLUX: "//
     $        "Perturbed flux surfaces"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
         WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $        "mstep =",mstep,"mthsurf =",mthsurf
         WRITE(out_unit,'(1x,a12,es17.8e3)')"scale =",rmax/(ximax*6.0)
         WRITE(out_unit,*)  
         WRITE(out_unit,'(7(1x,a16))')"r","z","real(r)","imag(r)",
     $        "real(z)","imag(z)","psi"   
         
         DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
            DO itheta=0,mthsurf
               WRITE(out_unit,'(7(es17.8e3))')
     $              rs(istep,itheta),zs(istep,itheta),
     $              REAL(rss(istep,itheta)),AIMAG(rss(istep,itheta)),
     $              REAL(zss(istep,itheta)),AIMAG(zss(istep,itheta)),
     $              psis(istep,itheta)
            ENDDO
         ENDDO
         
         CALL ascii_close(out_unit)
      ENDIF         
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_xbnormal
c-----------------------------------------------------------------------
c     subprogram 8. ipout_vbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vbnormal(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      
      INTEGER :: ipsi,ipert,i
      REAL(r8), DIMENSION(0:cmpsi) :: psi
      COMPLEX(r8), DIMENSION(:), POINTER :: vcmn

      INTEGER :: p_id,t_id,i_id,m_id,bm_id,q_id

      COMPLEX(r8), DIMENSION(mpert) :: vwpmn
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn
      REAL(r8), DIMENSION(cmpsi) :: qs
      REAL(r8), DIMENSION(cmpsi,mpert) :: xmns,ymns
      COMPLEX(r8), DIMENSION(cmpsi,mpert) :: vnomns,vwpmns
      COMPLEX(r8), DIMENSION(cmpsi,lmpert) :: pwpmns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing b normal components from coils"

      psi=(/(ipsi,ipsi=0,cmpsi)/)/REAL(cmpsi,r8)
      psi=psilow+(psilim-psilow)*SIN(psi*pi/2)**2
      ALLOCATE(vcmn(cmpert))
      vnomns = 0
      vwpmns = 0
      pwpmns = 0

      DO ipsi=1,cmpsi
         CALL spline_eval(sq,psi(ipsi),0)
         qs(ipsi)=sq%f(4)
         CALL field_bs_psi(psi(ipsi),vcmn,0)

         IF (bwp_pest_flag) THEN
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  pwpmns(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
            pwpmn=0
            pwpmn=pwpmns(ipsi,:)
            CALL ipeq_bcoords(psi(ipsi),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)             
            pwpmns(ipsi,:)=pwpmn
         ENDIF
         
         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  vwpmns(ipsi,cmlow-mlow+i)=vcmn(i)
               ENDIF
            ENDDO
            vwpmn=0
            vwpmn=vwpmns(ipsi,:)
            CALL ipeq_bcoords(psi(ipsi),vwpmn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            vnomns(ipsi,:)=vwpmn
            vwpmn=0
            vwpmn=vwpmns(ipsi,:)
            CALL ipeq_bcoords(psi(ipsi),vwpmn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,1)  
            vwpmns(ipsi,:)=vwpmn
         ELSE
            CALL field_bs_psi(psi(ipsi),vcmn,0)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  vnomns(ipsi,cmlow-mlow+i)=vcmn(i)
               ENDIF
            ENDDO

            CALL field_bs_psi(psi(ipsi),vcmn,2)
            DO i=1,cmpert
               IF ((cmlow-mlow+i>=1).AND.(cmlow-mlow+i<=mpert)) THEN
                  vwpmns(ipsi,cmlow-mlow+i)=vcmn(i)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(vcmn)
      ! append to netcdf file once this is (mstep,mpert)
c      IF(debug_flag) PRINT *,"Opening "//TRIM(fncfile)
c      CALL check( nf90_open(fncfile,nf90_write,fncid) )
c      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
c      CALL check( nf90_inq_dimid(fncid,"i",i_id) )
c      CALL check( nf90_inq_dimid(fncid,"m",m_id) )
c      CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
c      CALL check( nf90_inq_dimid(fncid,"theta",t_id) )
c      IF(debug_flag) PRINT *,"  Defining variables"
c      CALL check( nf90_redef(fncid))
c      CALL check( nf90_def_var(fncid, "q", nf90_double,(/p_id/),q_id) )
c      CALL check( nf90_put_att(fncid,q_id,"long_name","Safety Factor") )
c      CALL check( nf90_def_var(fncid, "b_n_x", nf90_double,
c     $            (/p_id,m_id,i_id/),bm_id) )
c      CALL check( nf90_put_att(fncid,bm_id,"long_name",
c     $            "Vacuum field normal to the flux surface") )
c      CALL check( nf90_put_att(fncid,bm_id,"units","Tesla") )
c      CALL check( nf90_put_att(fncid,bm_id,"jacobian",jac_out) )
c      CALL check( nf90_enddef(fncid) )
c      IF(debug_flag) PRINT *,"  Writting variables"
c      CALL check( nf90_put_var(fncid,bm_id,qs) )
c      CALL check( nf90_put_var(fncid,bm_id,RESHAPE((/REAL(vnomns),
c     $             AIMAG(vnomns)/),(/mstep,mpert,2/))) )
c      CALL check( nf90_close(fncid) )
c      IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)

      CALL ascii_open(out_unit,"ipec_vbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_VBNROMAL: "//
     $     "Normal components of coil field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mpsi =",cmpsi,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $     "real(bno)","imag(bno)","real(bwp)","imag(bwp)"

      DO ipsi=1,cmpsi
         DO ipert=1,mpert
            WRITE(out_unit,'(2(1x,es16.8),1x,I4,4(1x,es16.8))')
     $           psi(ipsi),qs(ipsi),mfac(ipert),
     $           REAL(vnomns(ipsi,ipert)),AIMAG(vnomns(ipsi,ipert)),
     $           REAL(vwpmns(ipsi,ipert)),AIMAG(vwpmns(ipsi,ipert))
         ENDDO
      ENDDO

      IF (bwp_pest_flag) THEN
         CALL ascii_open(out_unit,"ipec_vbnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VBNROMAL_PEST: "//
     $        "Normal components of coil field in pest coordinates"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(1x,a13,a8)')"jac_out = ","pest"
         WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $        "mpsi =",cmpsi,"mpert =",lmpert,"mthsurf =",mthsurf
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a16),1x,a4,2(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         
         DO ipsi=1,cmpsi
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,2(es17.8e3))')
     $              psi(ipsi),qs(ipsi),lmfac(ipert),
     $              REAL(pwpmns(ipsi,ipert)),AIMAG(pwpmns(ipsi,ipert))
            ENDDO
         ENDDO
      ENDIF

      IF (bin_flag) THEN
         CALL bin_open(bin_unit,
     $        "vbnormal.bin","UNKNOWN","REWIND","none")
         DO ipert=1,mpert
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
            xmns(ipsi,:)=mfac
            ymns(ipsi,:)=psi(ipsi)
         ENDDO
         CALL bin_open(bin_2d_unit,"vbnormal_spectrum.bin",
     $        "UNKNOWN","REWIND","none")
         WRITE(bin_2d_unit)1,0
         WRITE(bin_2d_unit)cmpsi-1,mpert-1
         WRITE(bin_2d_unit)REAL(xmns(1:cmpsi,:),4),
     $        REAL(ymns(1:cmpsi,:),4)
         WRITE(bin_2d_unit)REAL(ABS(vwpmns(1:cmpsi,:)),4)
         CALL bin_close(bin_2d_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      RETURN
      END SUBROUTINE ipout_vbnormal
c-----------------------------------------------------------------------
c     subprogram 9. ipout_xbtangent.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbtangent(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs,psis,
     $     rvecs,zvecs,vecs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,bs,dphi
      COMPLEX(r8), DIMENSION(mstep,mpert) :: xmns,ymns,xtamns,btamns
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwt_fun,xvt_fun,xvz_fun,
     $     bwt_fun,bvt_fun,bvz_fun,xta_fun,bta_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss,
     $     xtafuns,btafuns
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b tangential components"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))

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
            CALL ipeq_bcoords(psifac(istep),xta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL ipeq_bcoords(psifac(istep),bta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
         ENDIF
         xtamns(istep,:)=xta_mn
         btamns(istep,:)=bta_mn
         xtafuns(istep,:)=xtafuns(istep,:)*EXP(ifac*nn*dphi)
         btafuns(istep,:)=btafuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbtangent_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBNROMAL: "//
     $     "Tangential components of displacement and field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,a12,1x,I6,1x,2(a12,I4))')
     $     "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $     "real(xta)","imag(xta)","real(bta)","imag(bta)"

      DO istep=1,mstep,MAX(1,(mstep*mpert-1)/max_linesout+1)
         DO ipert=1,mpert
            WRITE(out_unit,'(2(es17.8e3),1x,I4,4(es17.8e3))')
     $           psifac(istep),qfac(istep),mfac(ipert),
     $           REAL(xtamns(istep,ipert)),AIMAG(xtamns(istep,ipert)),
     $           REAL(btamns(istep,ipert)),AIMAG(btamns(istep,ipert))
         ENDDO
      ENDDO

      IF (fun_flag) THEN
         CALL ascii_open(out_unit,"ipec_xbtangent_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"IPEC_XBTANGENT_FUN: "//
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
         DO ipert=1,mpert
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
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_xbtangent
c-----------------------------------------------------------------------
c     subprogram 10. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 10. ipout_xbrzphi.
c     write perturbed rzphi components on rzphi grid.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphi(egnum,xspmn,nr,nz,bnimn,bnomn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,nr,nz
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: bnimn,bnomn

      INTEGER :: i,j,k,l,ipert,iindex,np
      REAL(r8) :: mid,btlim,rlim,ileft,delr,delz,cha,chb,chc,chd,
     $   rij,zij,tij,t11,t12,t21,t22,t33
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

      REAL(r8), DIMENSION(:), POINTER :: chex,chey
      COMPLEX(r8), DIMENSION(:,:), POINTER :: chear,cheaz,
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

      IF (eqbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing equilibrium magnetic fields"
         IF(timeit) CALL ipec_timer(-2)
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
         IF(timeit) CALL ipec_timer(2)
      ENDIF

      CALL ipeq_alloc

      IF(verbose) WRITE(*,*)"Mapping fields to cylindrical coordinates"
      IF(timeit) CALL ipec_timer(-2)

      IF (brzphi_flag .OR. xrzphi_flag) THEN
         DO i=0,nr
         iindex = FLOOR(REAL(i,8)/FLOOR((nr-1)/10.0))*10
         ileft = REAL(i,8)/FLOOR((nr-1)/10.0)*10-iindex
         IF ((i /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a10)')
     $        "volume = ",iindex,"% mappings"
            DO j=0,nz
               IF (gdl(i,j)==1) THEN
                  CALL ipeq_sol(gdpsi(i,j))
                  CALL ipeq_contra(gdpsi(i,j))
                  CALL ipeq_cova(gdpsi(i,j))
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
                  ! Correction for helicity
                  IF(helicity>0)THEN
                     xrr(i,j) = CONJG(xrr(i,j))
                     brr(i,j) = CONJG(brr(i,j))
                     xrz(i,j) = CONJG(xrz(i,j))
                     brz(i,j) = CONJG(brz(i,j))
                     xrp(i,j) = CONJG(xrp(i,j))
                     brp(i,j) = CONJG(brp(i,j))
                  ENDIF
                  ! machine toroidal angle
                  xrr(i,j) = xrr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brr(i,j) = brr(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrz(i,j) = xrz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brz(i,j) = brz(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  xrp(i,j) = xrp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
                  brp(i,j) = brp(i,j)*EXP(-twopi*ifac*nn*gdphi(i,j))
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      IF(timeit) CALL ipec_timer(2)

      CALL ipeq_dealloc

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         IF(verbose) WRITE(*,*)
     $      "Computing vacuum fields by surface currents"
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
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
         IF(timeit) CALL ipec_timer(2)
      ENDIF
      
      IF (divzero_flag) THEN
         CALL ipeq_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF
      IF (div_flag) THEN
         CALL ipdiag_rzpdiv(nr,nz,gdl,gdr,gdz,brr,brz,brp)
      ENDIF

      IF (brzphi_flag) THEN
         IF(verbose) WRITE(*,*)"Computing total perturbed fields"
         bnomn=bnomn-bnimn
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
         CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $        nr,nz,vgdl,vgdr,vgdz,vpbr,vpbz,vpbp)
         IF (helicity<0) THEN
            vpbr=CONJG(vpbr)
            vpbz=CONJG(vpbz)
            vpbp=CONJG(vpbp)
         ENDIF

         IF (coil_flag) THEN
            IF(verbose) WRITE(*,*)"Computing vacuum fields by coils"
            np=nn*16
            CALL field_bs_rzphi(nr,nz,np,gdr,gdz,vcbr,vcbz,vcbp)
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

         CALL ascii_open(out_unit,"ipec_brzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for brzphi field"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,'(4(1x,a4,f12.3))')
     $        "a =",cha,"b =",chb,"c =",chc,"d =",chd
         WRITE(out_unit,*)   
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chear","cheaz",
     $        "real(chear)","imag(chear)","real(cheaz)","imag(cheaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $              REAL(chear(i,j)),AIMAG(chear(i,j)),
     $              REAL(cheaz(i,j)),AIMAG(cheaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

         CALL ascii_open(out_unit,"ipec_xrzphi_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_XRZPHI_CHEBYSHEV: "//
     $        "Chebyshev coefficients for displacement"
         WRITE(out_unit,*)version
         WRITE(out_unit,*)
         WRITE(out_unit,'(2(1x,a12,I4))')
     $        "nchr =",nchr,"nchz =",nchz
         WRITE(out_unit,*)     
         WRITE(out_unit,'(2(1x,a6),4(1x,a16))')"chxar","chxaz",
     $        "real(chxar)","imag(chxar)","real(chxaz)","imag(chxaz)"         
         DO i=0,nchr
            DO j=0,nchz
               WRITE(out_unit,'(2(1x,I6),4(es17.8e3))')i,j,
     $              REAL(chxar(i,j)),AIMAG(chxar(i,j)),
     $              REAL(chxaz(i,j)),AIMAG(chxaz(i,j))
            ENDDO
          ENDDO
         CALL ascii_close(out_unit)

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
         CALL ipvacuum_bnormal(psilim,bnomn,nr,nz)
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
         CALL ipvacuum_bnormal(psilim,bnimn,nr,nz)
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
      ! append to netcdf file once this is (mstep,mpert)
      IF(debug_flag) PRINT *,"Opening "//TRIM(cncfile)
      CALL check( nf90_open(cncfile,nf90_write,cncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(cncid,"i",i_id) )
      CALL check( nf90_inq_dimid(cncid,"R",r_id) )
      CALL check( nf90_inq_dimid(cncid,"z",z_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(cncid))
      CALL check( nf90_def_var(cncid, "b_r_equil", nf90_double,
     $            (/r_id,z_id/),bre_id) )
      CALL check( nf90_put_att(cncid,bre_id,"long_name",
     $            "Radial Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bre_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z_equil", nf90_double,
     $            (/r_id,z_id/),bze_id) )
      CALL check( nf90_put_att(cncid,bze_id,"long_name",
     $            "Vertical Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bze_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t_equil", nf90_double,
     $            (/r_id,z_id/),bpe_id) )
      CALL check( nf90_put_att(cncid,bpe_id,"long_name",
     $            "Toroidal Equilibrium Field") )
      CALL check( nf90_put_att(cncid,bpe_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_r_plasma", nf90_double,
     $            (/r_id,z_id,i_id/),brp_id) )
      CALL check( nf90_put_att(cncid,brp_id,"long_name",
     $            "Radial Plasma Field") )
      CALL check( nf90_put_att(cncid,brp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z_plasma", nf90_double,
     $            (/r_id,z_id,i_id/),bzp_id) )
      CALL check( nf90_put_att(cncid,bzp_id,"long_name",
     $            "Vertical Plasma Field") )
      CALL check( nf90_put_att(cncid,bzp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t_plasma", nf90_double,
     $            (/r_id,z_id,i_id/),bpp_id) )
      CALL check( nf90_put_att(cncid,bpp_id,"long_name",
     $            "Toroidal Plasma Field") )
      CALL check( nf90_put_att(cncid,bpp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_r", nf90_double,
     $            (/r_id,z_id,i_id/),br_id) )
      CALL check( nf90_put_att(cncid,br_id,"long_name",
     $            "Radial Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,br_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_z", nf90_double,
     $            (/r_id,z_id,i_id/),bz_id) )
      CALL check( nf90_put_att(cncid,bz_id,"long_name",
     $            "Vertical Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,bz_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "b_t", nf90_double,
     $            (/r_id,z_id,i_id/),bp_id) )
      CALL check( nf90_put_att(cncid,bp_id,"long_name",
     $            "Toroidal Nonaxisymmetric Field") )
      CALL check( nf90_put_att(cncid,bp_id,"units","Tesla") )
      CALL check( nf90_def_var(cncid, "xi_r", nf90_double,
     $            (/r_id,z_id,i_id/),xr_id) )
      CALL check( nf90_put_att(cncid,xr_id,"long_name",
     $            "Radial Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,xr_id,"units","m") )
      CALL check( nf90_def_var(cncid, "xi_z", nf90_double,
     $            (/r_id,z_id,i_id/),xz_id) )
      CALL check( nf90_put_att(cncid,xz_id,"long_name",
     $            "Vertical Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,bz_id,"units","m") )
      CALL check( nf90_def_var(cncid, "xi_t", nf90_double,
     $            (/r_id,z_id,i_id/),xp_id) )
      CALL check( nf90_put_att(cncid,xp_id,"long_name",
     $            "Toroidal Nonaxisymmetric Displacement") )
      CALL check( nf90_put_att(cncid,bp_id,"units","m") )
      CALL check( nf90_enddef(cncid) )
      IF(debug_flag) PRINT *,"  Writting variables"
      CALL check( nf90_put_var(cncid,bre_id,ebr) )
      CALL check( nf90_put_var(cncid,bze_id,ebz) )
      CALL check( nf90_put_var(cncid,bpe_id,ebp) )
      CALL check( nf90_put_var(cncid,br_id,RESHAPE((/REAL(btr),
     $             AIMAG(btr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bz_id,RESHAPE((/REAL(btz),
     $             AIMAG(btz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bp_id,RESHAPE((/REAL(btp),
     $             AIMAG(btp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,brp_id,RESHAPE((/REAL(bpr),
     $             AIMAG(bpr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bzp_id,RESHAPE((/REAL(bpz),
     $             AIMAG(bpz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,bpp_id,RESHAPE((/REAL(bpp),
     $             AIMAG(bpp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xr_id,RESHAPE((/REAL(xrr),
     $             AIMAG(xrr)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xz_id,RESHAPE((/REAL(xrz),
     $             AIMAG(xrz)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_put_var(cncid,xp_id,RESHAPE((/REAL(xrp),
     $             AIMAG(xrp)/),(/nr+1,nz+1,2/))) )
      CALL check( nf90_close(cncid) )
      IF(debug_flag) PRINT *,"Closed "//TRIM(cncfile)


      IF (eqbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_eqbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_EQBRZPHI: Eq. b field in rzphi grid"
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

      IF (brzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_brzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_BRZPHI: Total perturbed field"
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
            CALL ascii_open(out_unit,"ipec_chebrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CEHBRZPHI: Total perturbed field"
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

         CALL ascii_open(out_unit,"ipec_pbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_PBRZPHI: Perturbed field by plasma"
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

         IF (coil_flag) THEN
            CALL ascii_open(out_unit,"ipec_cbrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CBRZPHI: External field by coils"
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

      ENDIF

      IF (brzphi_flag .AND. vbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VBRZPHI: Total perturbed field "//
     $        "by surface currents"
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
     $              REAL(brr(i,j)),AIMAG(brr(i,j)),
     $              REAL(brz(i,j)),AIMAG(brz(i,j)),
     $              REAL(brp(i,j)),AIMAG(brp(i,j))
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)

      ENDIF

      IF (xrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_xrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_XRZPHI: Displacement"
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
            CALL ascii_open(out_unit,"ipec_chexrzphi_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"IPEC_CEHXRZPHI: Displacement"
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

      IF (pbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vpbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VPBRZPHI: Vacuum field by "//
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

      IF (vvbrzphi_flag) THEN
         CALL ascii_open(out_unit,"ipec_vvbrzphi_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"IPEC_VVBRZPHI: Vacuum field by "//
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
      CALL ipeq_rzpgrid(nr,nz,psixy) ! reset the grid
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_xbrzphi
c-----------------------------------------------------------------------
c     subprogram 10. ipout_vsbrzphi.
c     write brzphi components restored by removing shielding currents.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_vsbrzphi(snum,nr,nz)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: snum,nr,nz

      INTEGER :: i,j,ipert,iindex
      REAL(r8) :: mid,bt0,ileft

      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      LOGICAL, PARAMETER :: complex_flag=.TRUE.      

      INTEGER, DIMENSION(0:nr,0:nz) :: vgdl
      REAL(r8), DIMENSION(0:nr,0:nz) :: vgdr,vgdz
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: vbr,vbz,vbp
c-----------------------------------------------------------------------
c     build solutions.
c-----------------------------------------------------------------------
      vbr = 0
      vbz = 0
      vbp = 0

      IF (snum<10) THEN
         WRITE(UNIT=ss,FMT='(I1)')snum
         ss=TRIM(ADJUSTL(ss))
      ELSE
         WRITE(UNIT=ss,FMT='(I2)')snum
      ENDIF

      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vacuum fields by "//
     $     TRIM(ss)//"th resonant field"
      CALL ipvacuum_bnormal(singtype(snum)%psifac,
     $     singbno_mn(:,snum),nr,nz)
      CALL mscfld(wv,mpert,mthsurf,mthsurf,nfm2,nths2,complex_flag,
     $     nr,nz,vgdl,vgdr,vgdz,vbr,vbz,vbp)
      IF (helicity<0) THEN
         vbr=CONJG(vbr)
         vbz=CONJG(vbz)
         vbp=CONJG(vbp)
      ENDIF
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"ipec_vsbrzphi_n"//
     $     TRIM(sn)//"_s"//TRIM(ss)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_VSBRZPHI: Vacuum field in rzphi grid by "//
     $     TRIM(ss)//"th resonant field"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,2(a6,I6))')"nr =",nr+1,"nz =",nz+1
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a2,8(1x,a16))')"l","r","z",
     $     "real(vsb_r)","imag(vsb_r)","real(vsb_z)","imag(vsb_z)",
     $     "real(vsb_phi)","imag(vsb_phi)"
      
      DO i=0,nr
         DO j=0,nz
            WRITE(out_unit,'(1x,I2,8(es17.8e3))')
     $           vgdl(i,j),vgdr(i,j),vgdz(i,j),
     $           REAL(vbr(i,j)),AIMAG(vbr(i,j)),
     $           REAL(vbz(i,j)),AIMAG(vbz(i,j)),
     $           REAL(vbp(i,j)),AIMAG(vbp(i,j))
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_vsbrzphi
c-----------------------------------------------------------------------
c     subprogram 11. ipout_xbrzphifun
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xbrzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,dphi,t11,t12,t21,t22,t33
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun,
     $     xwt_fun,bwt_fun,xvz_fun,bvz_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: xrr_fun,brr_fun,
     $     xrz_fun,brz_fun,xrp_fun,brp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      IF(verbose) WRITE(*,*)"Computing x and b rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a29)')
     $        "volume = ",iindex,"% xi and b rzphi computations"
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
         CALL ipeq_cova(psifac(istep))
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

      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_xbrzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_XBRZPHI_FUN: "//
     $     "Rzphi components of displacement and field in functions"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(xr)","imag(xr)","real(xz)","imag(xz)",
     $     "real(xp)","imag(xp)","real(br)","imag(br)",
     $     "real(bz)","imag(bz)","real(bp)","imag(bp)"

      DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(es17.8e3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(xrr_fun(istep,itheta)),
     $           AIMAG(xrr_fun(istep,itheta)),
     $           REAL(xrz_fun(istep,itheta)),
     $           AIMAG(xrz_fun(istep,itheta)),
     $           REAL(xrp_fun(istep,itheta)),
     $           AIMAG(xrp_fun(istep,itheta)),
     $           REAL(brr_fun(istep,itheta)),
     $           AIMAG(brr_fun(istep,itheta)),
     $           REAL(brz_fun(istep,itheta)),
     $           AIMAG(brz_fun(istep,itheta)),
     $           REAL(brp_fun(istep,itheta)),
     $           AIMAG(brp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_xbrzphifun
c-----------------------------------------------------------------------
c     subprogram 12. ipout_arzphifun
c-----------------------------------------------------------------------
      SUBROUTINE ipout_arzphifun(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: rs,zs
      REAL(r8), DIMENSION(0:mthsurf) :: dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xss_fun,
     $     ear,eat,eap,arr,art,arp
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: ear_fun,eaz_fun,
     $     eap_fun,arr_fun,arz_fun,arp_fun
c-----------------------------------------------------------------------
c     compute solutions and contravariant/additional components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing vector potential rzphi functions"

      CALL idcon_build(egnum,xspmn)

      CALL ipeq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a37)')
     $        "volume = ",iindex,"% vector potential rzphi computations"
         CALL ipeq_sol(psifac(istep))
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

      CALL ipeq_dealloc

      CALL ascii_open(out_unit,"ipec_arzphi_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_ARZPHI_FUN: "//
     $     "Rzphi components of vector potential in functions"
      WRITE(out_unit,*)version
      WRITE(out_unit,*)
      WRITE(out_unit,'(1x,a13,a8)')"jac_out = ",jac_out
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1x,a12,1x,I6,1x,a12,I4)')
     $     "mstep =",mstep,"mthsurf =",mthsurf
      WRITE(out_unit,*)     
      WRITE(out_unit,'(14(1x,a16))')"r","z",
     $     "real(ear)","imag(ear)","real(eaz)","imag(eaz)",
     $     "real(eap)","imag(eap)","real(arr)","imag(arr)",
     $     "real(arz)","imag(arz)","real(arp)","imag(arp)"

      DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
         DO itheta=0,mthsurf
            WRITE(out_unit,'(14(es17.8e3))')
     $           rs(istep,itheta),zs(istep,itheta),
     $           REAL(ear_fun(istep,itheta)),
     $           AIMAG(ear_fun(istep,itheta)),
     $           REAL(eaz_fun(istep,itheta)),
     $           AIMAG(eaz_fun(istep,itheta)),
     $           REAL(eap_fun(istep,itheta)),
     $           AIMAG(eap_fun(istep,itheta)),
     $           REAL(arr_fun(istep,itheta)),
     $           AIMAG(arr_fun(istep,itheta)),
     $           REAL(arz_fun(istep,itheta)),
     $           AIMAG(arz_fun(istep,itheta)),
     $           REAL(arp_fun(istep,itheta)),
     $           AIMAG(arp_fun(istep,itheta))
         ENDDO
      ENDDO
      
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_arzphifun
c-----------------------------------------------------------------------
c     subprogram 13. ipout_xclebsch.
c     Write Clebsch coordinate displacements.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_xclebsch(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,j,istep,ipert,iindex,ids(3)
      INTEGER :: i_id,m_id,p_id,dp_id,xp_id,xa_id
      REAL(r8) :: ileft

      COMPLEX(r8), DIMENSION(mstep,mpert) :: xmp1mns,xspmns,xmsmns
      COMPLEX(r8), DIMENSION(mstep,lmpert) :: xmp1out,xspout,xmsout
c-----------------------------------------------------------------------
c     compute necessary components.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing Clebsch displacements"

      CALL idcon_build(egnum,xspmn)
      CALL ipeq_alloc

      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% Clebsch decomposition"
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psifac(istep),1)
         CALL ipeq_sol(psifac(istep))
         CALL ipeq_contra(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b variations in hamada.
c-----------------------------------------------------------------------
         xmp1mns(istep,:) = xmp1_mn
         xspmns(istep,:) = xsp_mn
         xmsmns(istep,:) = xms_mn
      ENDDO
         
      CALL ascii_open(out_unit,"ipec_xclebsch_n"//
     $   TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"IPEC_PMODB: "//
     $   "Clebsch Components of the displacement."
      WRITE(out_unit,*)version
      WRITE(out_unit,'(1/,1x,a13,a8)')"jac_type = ",jac_type
      WRITE(out_unit,'(1/,1x,a12,1x,I6,1x,2(a12,I4),1/)')
     $   "mstep =",mstep,"mpert =",mpert,"mthsurf =",mthsurf
      WRITE(out_unit,'(1x,a23,1x,a4,6(1x,a16))')"psi","m",
     $   "real(derxi^psi)","imag(derxi^psi)",
     $   "real(xi^psi)","imag(xi^psi)",
     $   "real(xi^alpha)","imag(xi^alpha)"
      DO istep=1,mstep,MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
         DO ipert=1,mpert
            WRITE(out_unit,'(1x,es23.15,1x,I4,6(es17.8e3))')
     $         psifac(istep),mfac(ipert),xmp1mns(istep,ipert),
     $         xspmns(istep,ipert),xmsmns(istep,ipert)/chi1
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)

      ! convert to jac_out
      xmp1out = 0
      xspout = 0
      xmsout = 0
      DO istep=1,mstep
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               xmp1out(istep,mlow-lmlow+i) = xmp1mns(istep,i)
               xspout(istep,mlow-lmlow+i) = xspmns(istep,i)
               xmsout(istep,mlow-lmlow+i) = xmsmns(istep,i)
            ENDIF
         ENDDO
      ENDDO

      ! append to netcdf file
      IF(debug_flag) PRINT *,"Opening "//TRIM(fncfile)
      CALL check( nf90_open(fncfile,nf90_write,fncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(fncid,"i",i_id) )
      CALL check( nf90_inq_dimid(fncid,"m",m_id) )
      CALL check( nf90_inq_dimid(fncid,"psi_n",p_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(fncid))
      CALL check( nf90_def_var(fncid, "xi_psi1", nf90_double,
     $            (/p_id,m_id,i_id/),dp_id) )
      CALL check( nf90_put_att(fncid,dp_id,"long_name",
     $            "Psi derivative of contravariant psi displacement") )
      CALL check( nf90_def_var(fncid, "xi_psi", nf90_double,
     $            (/p_id,m_id,i_id/),xp_id) )
      CALL check( nf90_put_att(fncid,xp_id,"long_name",
     $            "contravariant psi displacement") )
      CALL check( nf90_def_var(fncid, "xi_alpha", nf90_double,
     $            (/p_id,m_id,i_id/),xa_id) )
      CALL check( nf90_put_att(fncid,xa_id,"long_name",
     $            "contravariant Clebsch angle displacement") )
      ids = (/xp_id,dp_id,xa_id/)
      DO i=1,3
         CALL check( nf90_put_att(fncid,ids(i),"units","m") )
         CALL check( nf90_put_att(fncid,ids(i),"jacobian",jac_type) )
      ENDDO
      CALL check( nf90_enddef(fncid) )
      IF(debug_flag) PRINT *,"  Writting variables"
      CALL check( nf90_put_var(fncid,dp_id,RESHAPE((/REAL(xmp1out),
     $             AIMAG(xmp1out)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,xp_id,RESHAPE((/REAL(xspout),
     $             AIMAG(xspout)/),(/mstep,lmpert,2/))) )
      CALL check( nf90_put_var(fncid,xa_id,RESHAPE((/REAL(xmsout),
     $             AIMAG(xmsout)/),(/mstep,lmpert,2/))) )

      CALL check( nf90_close(fncid) )
      IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)


      CALL ipeq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_xclebsch
c-----------------------------------------------------------------------
c     subprogram 14. ipout_control_filter.
c     Filter control surface flux vector in flux bases with energy norms
c-----------------------------------------------------------------------
      SUBROUTINE ipout_control_filter(finmn,foutmn,ftypes,fmodes,
     $           rout,bpout,bout,rcout,tout,jout,op_write)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: fmodes,rout,bpout,bout,rcout,tout,jout
      CHARACTER(len=*), INTENT(IN) :: ftypes
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn,foutmn
      LOGICAL, INTENT(IN), OPTIONAL :: op_write

      ! eigendecompositions variables
      INTEGER :: info,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: work2
      ! SVD variables
      REAL(r8), DIMENSION(5*mpert) :: rworksvd
      REAL(r8), DIMENSION(5*msing) :: sworksvd
      COMPLEX(r8), DIMENSION(3*mpert) :: worksvd

      LOGICAL :: output
      INTEGER :: i,j,k,ipert,maxmode
      INTEGER :: idid,mdid,xdid,wdid,rdid,pdid,sdid,tdid,
     $   mx_id,mw_id,mr_id,mp_id,ms_id,
     $   we_id,re_id,pe_id,se_id,
     $   w_id,r_id,p_id,s_id, sc_id,wr_id,wp_id,rp_id,ws_id,rs_id,ps_id,
     $   ft_id,fx_id,wx_id,rx_id,px_id,sx_id,wa_id,ra_id,rl_id,
     $   x_id,xe_id,xt_id,wf_id,rf_id,sf_id,sm_id,ex_id,et_id,
     $   wev_id,wes_id,wep_id,rev_id,res_id,rep_id,sev_id,ses_id,sep_id
      REAL(r8) :: norm
      REAL(r8), DIMENSION(mpert) :: singfac
      REAL(r8), DIMENSION(mpert,2) :: tempmi
      REAL(r8), DIMENSION(msing,2) :: tempsi
      COMPLEX(r8), DIMENSION(msing) :: temps
      COMPLEX(r8), DIMENSION(mpert) :: temp,tempm,eigmn,filmn
      COMPLEX(r8), DIMENSION(lmpert) :: templ
      COMPLEX(r8), DIMENSION(mpert,mpert) :: mat,matmm,tempmm,singmat,
     $    sqrta,sqrtainv
      COMPLEX(r8), DIMENSION(mpert,msing) :: matms
      COMPLEX(r8), DIMENSION(msing,msing) :: matss
      COMPLEX(r8), DIMENSION(msing,mpert) :: matsm,singcoupmat
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat
      CHARACTER(64) :: message

      INTEGER,  DIMENSION(mpert) :: aindx,indx,raindx,rlindx
      REAL(r8), DIMENSION(mpert) :: xvals,wvals,rvals,pvals,avals,
     $    rlvals,sengys,vengys,pengys
      REAL(r8), DIMENSION(msing) :: svals
      COMPLEX(r8), DIMENSION(mpert,mpert)::xvecs,wvecs,rvecs,pvecs,avecs
      COMPLEX(r8), DIMENSION(mpert,msing) :: svecs
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: xveco,wveco,rveco,pveco
      COMPLEX(r8), DIMENSION(lmpert,msing) :: sveco
      COMPLEX(r8), DIMENSION(0:mthsurf,mpert)::wfuns,rfuns
      COMPLEX(r8), DIMENSION(0:mthsurf,msing)::sfuns

      IF(timeit) CALL ipec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing Energy-Normalized flux bases"
c-----------------------------------------------------------------------
c     basic definitions
c-----------------------------------------------------------------------
      singfac=mfac-nn*qlim
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
c     calculate sqrt(A) weighting matrix.
c      - Define sqrt(A) weighting matrix as
c        W_m,m' = int{sqrt(J|delpsi|)exp[-i*(m-m')t]dt}/int{sqrt(J|delpsi|)dt}
c-----------------------------------------------------------------------
      DO i=1,mpert
         temp = 0
         temp(i) = 1.0
         CALL ipeq_weight(psilim,temp,mfac,mpert,2) ! A^1/2
         sqrta(:,i) = temp
      ENDDO
      ! inverse of the 1/2 area weighting
      lwork=2*mpert-1
      work=0
      work2=0
      rwork=0
      mat=0
      DO i=1,mpert
         mat(i,i)=1
      ENDDO
      matmm=sqrta
      CALL zhetrf('L',mpert,matmm,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,matmm,mpert,ipiv,mat,mpert,info)
      sqrtainv = mat
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
c      - We want to use Phi'=Bsqrt(A) basis so a eigenvector norm means
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
      ! convert to bsqrt(A)
      wvecs=MATMUL(MATMUL(sqrta,wvecs),sqrta)*2*mu0
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
      mat=MATMUL(MATMUL(sqrta,surf_indinvmats),sqrta)
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
c      - We want to use Phi'=Bsqrt(A) basis so a eigenvector norm means
c        int{b^2da}/int{da} = 1
c        - This is a physiCALLy meaningful quantity (energy), and thus
c          independent of coordinate system
c      - Normalizing the input, I = RPhi becomes I = RWPhi'
c      - Normalizing the output gives W^dagger I = W^daggerRW Phi'
c      - Since W is Hermitian (W^dagger=W) the new operator is too
c      - Eigenvalues correspond to int{I^2da} (power) per int{b^2da} (energy)
c      - NOTE: No need to include 1/jarea=1/int{da} (gets normalized)
c-----------------------------------------------------------------------
      ! start with IPEC flux matrix
      ! remove border of modes/solutions (diagnostic only)
      i = malias+1
      j = mpert-malias
      rvecs = 0
      rvecs(i:j,i:j) = reluctmats(resp_index,i:j,i:j)
      ! convert to bsqrt(A)
      rvecs=MATMUL(MATMUL(sqrta,rvecs),sqrta)
      work = 0
      rwork = 0
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,rvecs,mpert,rvals,work,lwork,rwork,info)
      ! put in descending order like zgesvd
      rvals(:)   = rvals(mpert:1:-1)
      rvecs(:,:) = rvecs(:,mpert:1:-1)
c-----------------------------------------------------------------------
c     Calculate permeability right-singular vectors and singular values
c      - We want to use Phi'=Bsqrt(A) basis so a eigenvector norm means
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
      ! convert to bsqrt(A)
      mat = MATMUL(MATMUL(sqrta,mat),sqrta)
      mat = TRANSPOSE(mat)
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
            CALL ipeq_weight(psilim,temp,mfac,mpert,2)
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
         matsm=MATMUL(singbnoflxs,tempmm)
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
            temp = finmn
            filmn = 0
            CALL ipeq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
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
               CALL ipeq_weight(psilim,filmn,mfac,mpert,2) ! sqrt(A)b to flux
            ENDIF
            finmn = filmn
         ENDDO
      ENDIF
      foutmn = MATMUL(permeabmats(resp_index,:,:),finmn) ! total flux
c-----------------------------------------------------------------------
c     Write outputs
c-----------------------------------------------------------------------
      IF(output)THEN
         ! Convert to output coords
         DO i=1,mpert
           templ = 0
           templ(mlow-lmlow+i) = 1.0
           IF((jac_out /= jac_type).OR.(tout==0)) CALL ipeq_bcoords(
     $        psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,0)
           coordmat(:,i) = templ
         ENDDO
         wveco = MATMUL(coordmat,wvecs)
         rveco = MATMUL(coordmat,rvecs)
         pveco = MATMUL(coordmat,pvecs)
         IF(singcoup_set) sveco = MATMUL(coordmat,svecs)
         !xveco = MATMUL(coordmat,MATMUL(singmat,xvecs))
         !DO i=1,mpert
         !  xveco(:,i) = xveco(:,i)/(chi1*twopi*ifac*(lmfac-nn*qlim)) ! not right if tout=0?
         !ENDDO
         xveco = xvecs

         ! Write to netCDF file
         IF(debug_flag) PRINT *,"Opening "//TRIM(mncfile)
         CALL check( nf90_open(mncfile,nf90_write,mncid) )
         IF(debug_flag) PRINT *,"  Inquiring about dimensions"
         CALL check( nf90_inq_dimid(mncid,"i",idid) )
         CALL check( nf90_inq_dimid(mncid,"m",mdid) )
         CALL check( nf90_inq_dimid(mncid,"theta",tdid) )

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
         CALL check( nf90_put_att(mncid,xe_id,"units","J/(Wb/m)^2") )
         CALL check( nf90_put_att(mncid,xe_id,"long_name",
     $    "Total displacement energy eigenvalues") )

         CALL check( nf90_def_dim(mncid,"mode_W",mpert,   wdid) )
         CALL check( nf90_def_var(mncid,"mode_W",nf90_int,wdid,mw_id))
         CALL check( nf90_put_att(mncid, mw_id ,"long_name",
     $    "Energy-norm external flux energy eigenmode index") )
         CALL check( nf90_def_var(mncid,"W_x_eigenvector",
     $               nf90_double,(/mdid,wdid,idid/),w_id) )
         CALL check( nf90_put_att(mncid,w_id,"long_name",
     $    "Energy-norm external flux energy eigendecomposition") )
         CALL check( nf90_def_var(mncid,"W_x_eigenvalue",
     $               nf90_double,(/wdid/),we_id) )
         CALL check( nf90_put_att(mncid,we_id,"units","J/(Wb/m)^2") )
         CALL check( nf90_put_att(mncid,we_id,"long_name",
     $    "Energy-norm external flux energy eigenvalues") )
         CALL check( nf90_def_var(mncid,"W_x_amp",nf90_double,
     $                         (/wdid/),wa_id) )
         CALL check( nf90_put_att(mncid,wa_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode amplifications") )
         CALL check( nf90_def_var(mncid,"W_x_energyv",nf90_double,
     $                         (/wdid/),wev_id) )
         CALL check( nf90_put_att(mncid,wev_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode vacuum energy") )
         CALL check( nf90_def_var(mncid,"W_x_energys",nf90_double,
     $                         (/wdid/),wes_id) )
         CALL check( nf90_put_att(mncid,wes_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode surface energy") )
         CALL check( nf90_def_var(mncid,"W_x_energyp",nf90_double,
     $                         (/wdid/),wep_id) )
         CALL check( nf90_put_att(mncid,wep_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode total energy") )

         CALL check( nf90_def_dim(mncid,"mode_R",mpert,   rdid) )
         CALL check( nf90_def_var(mncid,"mode_R",nf90_int,rdid,mr_id))
         CALL check( nf90_put_att(mncid, mr_id ,"long_name",
     $    "Energy-norm external flux reluctance eigenmode index"))
         CALL check( nf90_def_var(mncid,"R_xe_eigenvector",nf90_double,
     $               (/mdid,rdid,idid/),r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $    "Energy-norm external flux reluctance eigendecomposition") )
         CALL check( nf90_def_var(mncid,"R_xe_eigenvalue",nf90_double,
     $                         (/rdid/),re_id) )
         CALL check( nf90_put_att(mncid,re_id,"long_name",
     $    "Energy-norm external flux reluctance eigenvalues") )
         CALL check( nf90_put_att(mncid,re_id,"units","A/(Wb/m)") )
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

         CALL check( nf90_def_var(mncid,"O_XPhi",nf90_double,
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
            CALL check( nf90_def_dim(mncid,"mode_CX",msing,   sdid) )
            CALL check( nf90_def_var(mncid,"mode_CX",nf90_int,sdid,
     $                  ms_id) )
            CALL check( nf90_put_att(mncid, ms_id ,"long_name",
     $       "Singular coupling mode index") )
            CALL check( nf90_def_var(mncid,"C",nf90_double,
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

         CALL check( nf90_def_var(mncid,"Phi_xe_m",nf90_double,
     $               (/mdid,idid/),ex_id) )
         CALL check( nf90_put_att(mncid,ex_id,"units","Wb/m") )
         CALL check( nf90_put_att(mncid,ex_id,"long_name",
     $    "Energy-norm external flux") )
         CALL check( nf90_def_var(mncid,"Phi_x_m",nf90_double,
     $               (/mdid,idid/),fx_id) )
         CALL check( nf90_put_att(mncid,fx_id,"units","Wb") )
         CALL check( nf90_put_att(mncid,fx_id,"long_name",
     $    "External flux") )
         CALL check( nf90_def_var(mncid,"Phi_e_m",nf90_double,
     $               (/mdid,idid/),et_id) )
         CALL check( nf90_put_att(mncid,et_id,"units","Wb/m") )
         CALL check( nf90_put_att(mncid,et_id,"long_name",
     $    "Energy-norm total flux") )
         CALL check( nf90_def_var(mncid,"Phi_m",nf90_double,
     $               (/mdid,idid/),ft_id) )
         CALL check( nf90_put_att(mncid,ft_id,"units","Wb") )
         CALL check( nf90_put_att(mncid,ft_id,"long_name",
     $    "Total flux") )

         IF(fun_flag)THEN
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
            CALL check( nf90_put_var(mncid,ms_id,indx(:msing)) )
            CALL check( nf90_put_var(mncid,sc_id,RESHAPE(
     $       (/REAL(matms),AIMAG(matms)/),(/mpert,msing,2/))))
            CALL check( nf90_put_var(mncid,s_id,RESHAPE(
     $       (/REAL(svecs),AIMAG(svecs)/),(/mpert,msing,2/))))
            CALL check( nf90_put_var(mncid,se_id,svals) )
         ENDIF

         ! Energies for postprocessing re-normalization
         mat=MATMUL(MATMUL(sqrta,plas_indinvmats(resp_index,:,:)),sqrta)
         matmm=MATMUL(MATMUL(sqrta,surf_indinvmats),sqrta)
         tempmm=plas_indmats(resp_index,:,:)
         tempmm=MATMUL(CONJG(TRANSPOSE(tempmm)),surf_indinvmats)
         tempmm=MATMUL(tempmm,plas_indinvmats(resp_index,:,:))
         tempmm=MATMUL(MATMUL(sqrta,tempmm),sqrta) ! P-dagger.L^-1.P
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

         ! L.rho = -(1+s)/s ??
         matmm=MATMUL(MATMUL(sqrtainv,surf_indmats),sqrtainv)
         DO i=1,mpert
            tempm = rvecs(:,i)
            rlvals(i) = rvals(i)*DOT_PRODUCT(tempm,MATMUL(matmm,tempm))
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
         temp = foutmn
         CALL ipeq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
         templ = MATMUL(coordmat,temp) ! output coords
         !CALL ipeq_bcoordsout(templ,temp,psilim,tout,0)
         CALL check( nf90_put_var(mncid,et_id,RESHAPE((/REAL(temp),
     $               AIMAG(temp)/),(/mpert,2/))) )
         templ = MATMUL(coordmat,foutmn) ! output coords
         !CALL ipeq_bcoordsout(templ,foutmn,psilim,tout,0)
         CALL check( nf90_put_var(mncid,ft_id,RESHAPE((/REAL(foutmn),
     $               AIMAG(foutmn)/),(/mpert,2/))) )
         temp = finmn
         CALL ipeq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
         templ = MATMUL(coordmat,temp) ! output coords
         !CALL ipeq_bcoordsout(templ,temp,psilim,tout,0)
         CALL check( nf90_put_var(mncid,ex_id,RESHAPE((/REAL(temp),
     $               AIMAG(temp)/),(/mpert,2/))) )
         templ = MATMUL(coordmat,finmn) ! output coords
         !CALL ipeq_bcoordsout(templ,finmn,psilim,tout,0)
         CALL check( nf90_put_var(mncid,fx_id,RESHAPE((/REAL(finmn),
     $               AIMAG(finmn)/),(/mpert,2/))) )

         ! Decomposition of the applied flux
         temp = finmn
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
     $               AIMAG(tempm)/),(/msing,2/))) )
         ENDIF

         ! Decomposition of displacement
         temp = foutmn/(chi1*twopi*ifac*(mfac-nn*qlim)) ! total displacement
         tempm = MATMUL(CONJG(TRANSPOSE(xvecs)),temp)
         CALL check( nf90_put_var(mncid,xt_id,RESHAPE((/REAL(tempm),
     $               AIMAG(tempm)/),(/mpert,2/))) )

         ! Eigenmodes in real space (R,z) written in ipout_control
         IF(fun_flag)THEN
            DO i=1,mpert
                CALL iscdftb(mfac,mpert,wfuns(:,i),mthsurf,wvecs(:,i))
                CALL iscdftb(mfac,mpert,rfuns(:,i),mthsurf,rvecs(:,i))
                IF(singcoup_set .AND. i<=msing)
     $            CALL iscdftb(mfac,mpert,sfuns(:,i),mthsurf,svecs(:,i))
            ENDDO
            CALL check( nf90_put_var(mncid,wf_id,RESHAPE((/REAL(wfuns),
     $             -helicity*AIMAG(wfuns)/),(/mthsurf+1,mpert,2/))) )
            CALL check( nf90_put_var(mncid,rf_id,RESHAPE((/REAL(rfuns),
     $             -helicity*AIMAG(rfuns)/),(/mthsurf+1,mpert,2/))) )
            IF(singcoup_set)
     $       CALL check( nf90_put_var(mncid,sf_id,RESHAPE((/REAL(sfuns),
     $             -helicity*AIMAG(sfuns)/),(/mthsurf+1,msing,2/))) )
         ENDIF

         CALL check( nf90_close(mncid) )
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL ipec_timer(2)
      RETURN
      END SUBROUTINE ipout_control_filter
c-----------------------------------------------------------------------
c     subprogram 15. check.
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
c     subprogram 16. ipout_init_netcdf.
c     Initialize the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_init_netcdf
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER:: i,midid,mmdid,medid,mvdid,msdid,mtdid,
     $   fidid,fmdid,fpdid,ftdid,cidid,crdid,czdid,clvid,
     $   mivid,mmvid,mevid,mvvid,msvid,mtvid,fivid,fmvid,fpvid,ftvid,
     $   civid,crvid,czvid,id,fileids(3)
      INTEGER, DIMENSION(mpert) :: mmodes
c-----------------------------------------------------------------------
c     set variables
c-----------------------------------------------------------------------
      IF(debug_flag) PRINT *,"Initializing NETCDF files"
      mmodes = (/(i,i=1,mpert)/)
      mncfile = "ipec_control_output_n"//TRIM(sn)//".nc"
      fncfile = "ipec_profile_output_n"//TRIM(sn)//".nc"
      cncfile = "ipec_cylindrical_output_n"//TRIM(sn)//".nc"
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
     $     "IPEC outputs in Fourier or alternate modal bases"))
      CALL check( nf90_def_dim(mncid,"i",2,       midid) )
      CALL check( nf90_def_var(mncid,"i",nf90_int,midid,mivid) )
      CALL check( nf90_def_dim(mncid,"m",mpert,  mmdid) )
      CALL check( nf90_def_var(mncid,"m",nf90_int,mmdid,mmvid) )
      CALL check( nf90_def_dim(mncid,"mode",mpert,   medid) )
      CALL check( nf90_def_var(mncid,"mode",nf90_int,medid,mevid))
      CALL check( nf90_def_dim(mncid,"theta",mthsurf+1,  mtdid) )
      CALL check( nf90_def_var(mncid,"theta",nf90_double,mtdid,mtvid) )
      CALL check( nf90_put_att(mncid,nf90_global,"shot",INT(shotnum)) )
      CALL check( nf90_put_att(mncid,nf90_global,"time",INT(shottime)) )
      CALL check( nf90_put_att(mncid,nf90_global,"machine",machine) )
      CALL check( nf90_put_att(mncid,nf90_global,"n",nn) )
      CALL check( nf90_put_att(mncid,nf90_global,"jacobian",jac_type))
      CALL check( nf90_put_att(mncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(mncid,nf90_global,"psi_n_lim",psilim))
      CALL check( nf90_put_att(mncid,nf90_global,"version",version))
      
      IF(debug_flag) PRINT *," - Defining flux netcdf globals"
      CALL check( nf90_put_att(fncid,nf90_global,"title",
     $     "IPEC outputs in magnetic coordinate systems"))
      CALL check( nf90_def_dim(fncid,"i",2,       fidid) )
      CALL check( nf90_def_var(fncid,"i",nf90_int,fidid,fivid) )
      CALL check( nf90_def_dim(fncid,"m",lmpert,  fmdid) )
      CALL check( nf90_def_var(fncid,"m",NF90_INT,fmdid,fmvid) )
      CALL check( nf90_def_dim(fncid,"psi_n",mstep,    fpdid) )
      CALL check( nf90_def_var(fncid,"psi_n",nf90_double,fpdid,fpvid) )
      CALL check( nf90_def_dim(fncid,"theta",mthsurf+1,  ftdid) )
      CALL check( nf90_def_var(fncid,"theta",nf90_double,ftdid,ftvid) )
      CALL check( nf90_put_att(fncid,nf90_global,"shot",INT(shotnum)) )
      CALL check( nf90_put_att(fncid,nf90_global,"time",INT(shottime)) )
      CALL check( nf90_put_att(fncid,nf90_global,"machine",machine) )
      CALL check( nf90_put_att(fncid,nf90_global,"n",nn) )
      CALL check( nf90_put_att(fncid,nf90_global,"version",version))

      IF(debug_flag) PRINT *," - Defining cylindrical netcdf globals"
      CALL check( nf90_put_att(cncid,nf90_global,"title",
     $     "IPEC outputs in (R,z) coordinates"))
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
      CALL check( nf90_put_var(mncid,mevid,mmodes) )
      CALL check( nf90_put_var(mncid,mtvid,theta) )

      IF(debug_flag) PRINT *," - Putting coordinates in flux netcdfs"
      CALL check( nf90_put_var(fncid,fivid,(/0,1/)) )
      CALL check( nf90_put_var(fncid,fmvid,lmfac) )
      CALL check( nf90_put_var(fncid,fpvid,psifac(1:)) )
      CALL check( nf90_put_var(fncid,ftvid,theta) )

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
      END SUBROUTINE ipout_init_netcdf
c-----------------------------------------------------------------------
c     subprogram 17. ipout_close_netcdf.
c     Close the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE ipout_close_netcdf
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
      END SUBROUTINE ipout_close_netcdf

      END MODULE ipout_mod
