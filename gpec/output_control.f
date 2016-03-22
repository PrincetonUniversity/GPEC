c-----------------------------------------------------------------------
c     GENERALIZED PERTURBED EQUILIBRIUM CODE
c     write various output results of gpec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. output_control
c      1. response
c      2. singcoup
c      3. control
c      4. singfld
c      5. vsingfld
c-----------------------------------------------------------------------
c     subprogram 0. output_control.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE output_control
      USE gpec_response
      USE ipvacuum_mod
      USE gpec_diagnostic
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
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: singcoup1mat,
     $   w1v,w2v,w3v,t1v,t2v,t3v,fldflxmat

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
c     subprogram 1. response.
c     write basic information.
c-----------------------------------------------------------------------
      SUBROUTINE response(rout,bpout,bout,rcout,tout,jout)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      INTEGER :: lwork,i,j,k,rout,bpout,bout,rcout,tout,jout
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(mpert) :: sts,s,ev
      COMPLEX(r8), DIMENSION(mpert) :: cev
      COMPLEX(r8), DIMENSION(mpert,mpert) :: a,b,u,vt

      REAL(r8), DIMENSION(5*mpert) :: rwork
      COMPLEX(r8), DIMENSION(3*mpert) :: work

      INTEGER :: idid,mdid,edid,l_id,r_id,la_id,p_id
      COMPLEX(r8), DIMENSION(lmpert) :: vL,vL1,vLi,vP,vP1,vR,vW,templ
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: matlm,coordmat
      COMPLEX(r8), DIMENSION(0:mthsurf,mpert) :: wtfun,ilfun,rfun,pfun
      CHARACTER(2048) :: header

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
         sts(i)=-et(i)/(surfei(i)+surfee(i))
         IF (sts(i) > 0) s(i)=-s(i)
         s(i)=-1/s(i)
      ENDDO

      CALL ascii_open(out_unit,"gpec_response_n"//
     $	   TRIM(sn)//".out","UNKNOWN")
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
             CALL peq_bcoords(psilim,vL,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vL1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vLi,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vP,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vP1,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vR,lmfac,lmpert,
     $            rout,bpout,bout,rcout,tout,jout)
             CALL peq_bcoords(psilim,vW,lmfac,lmpert,
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
        CALL ascii_open(out_unit,"gpec_response_fun_n"//
     $    TRIM(sn)//".out","UNKNOWN")
        WRITE(out_unit,*)"GPEC_RESPONSE_FUN: "//
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
      CALL check( nf90_put_att(mncid,r_id,"long_name",
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
        IF((jac_out /= jac_type).OR.(tout==0)) CALL peq_bcoords(
     $     psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,0)
        coordmat(:,i) = templ
      ENDDO
      matlm = MATMUL(coordmat,surf_indmats)
      CALL check( nf90_put_var(mncid,l_id,RESHAPE((/REAL(matlm),
     $            AIMAG(matlm)/),(/lmpert,mpert,2/))) )
      matlm = MATMUL(coordmat,plas_indmats(resp_index,:,:))
      CALL check( nf90_put_var(mncid,la_id,RESHAPE((/REAL(matlm),
     $            AIMAG(matlm)/),(/lmpert,mpert,2/))) )
      matlm = MATMUL(coordmat,permeabmats(resp_index,:,:))
      CALL check( nf90_put_var(mncid,p_id,RESHAPE((/REAL(matlm),
     $            AIMAG(matlm)/),(/lmpert,mpert,2/))) )
      matlm = MATMUL(coordmat,reluctmats(resp_index,:,:))
      CALL check( nf90_put_var(mncid,l_id,RESHAPE((/REAL(matlm),
     $            AIMAG(matlm)/),(/lmpert,mpert,2/))) )
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE response
c-----------------------------------------------------------------------
c     subprogram 2. singcoup.
c     compute coupling between singular surfaces and external fields.
c-----------------------------------------------------------------------
      SUBROUTINE singcoup(spot,rout,bpout,bout,rcout,tout)
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
     $     singcurs,singbnoflxs,islandhwids
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
      CALL peq_alloc
      DO i=1,mpert
         finmn=0
         finmn(i)=1.0
         CALL peq_weight(psilim,finmn,mfac,mpert,1)
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
            CALL peq_sol(lpsi)
            lbwp1mn=bwp1_mn(resnum)
            
            rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1)) 
            CALL peq_sol(rpsi)
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
      CALL peq_dealloc
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
            CALL peq_weight(psilim,fldflxmn,mfac,mpert,2)
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
      CALL ascii_open(out_unit,"gpec_singcoup_matrix_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_SINGCOUP_MATRIX: Coupling matrices"//
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
      
      CALL ascii_open(out_unit,"gpec_singcoup_svd_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_SINGCOUP_SVD: SVD analysis"//
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE singcoup
c-----------------------------------------------------------------------
c     subprogram 3. control
c     calculate response from external field on the control surface.
c-----------------------------------------------------------------------
      SUBROUTINE control(ifile,finmn,foutmn,xspmn,
     $     rin,bpin,bin,rcin,tin,jin,rout,bpout,bout,rcout,tout,
     $     filter_types,filter_modes,filter_out)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(128), INTENT(IN) :: ifile
      INTEGER, INTENT(IN) :: rin,bpin,bin,rcin,tin,jin,
     $     rout,bpout,bout,rcout,tout,filter_modes
      LOGICAL, INTENT(IN) :: filter_out
      CHARACTER(len=*), INTENT(IN) :: filter_types
      
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: finmn
      COMPLEX(r8), DIMENSION(mpert), INTENT(OUT) :: foutmn,xspmn

      INTEGER :: i,j,i1,i2,i3,ms,itheta,jout
      INTEGER, DIMENSION(mpert) :: ipiv
      REAL(r8) :: vengy,sengy,pengy,area,thetai,scale,norm
      COMPLEX(r8) :: vy,sy,py
      CHARACTER(128) :: message

      REAL(r8), DIMENSION(0:mthsurf) :: dphi,delpsi,thetas,jacs,rvecs,
     $     zvecs,sbinfun,sboutfun

      COMPLEX(r8), DIMENSION(mpert) :: binmn,boutmn,xinmn,xoutmn,tempmn,
     $      abinmn
      COMPLEX(r8), DIMENSION(lmpert) :: cinmn,coutmn,cawmn,acinmn,templ
      COMPLEX(r8), DIMENSION(0:mthsurf) :: binfun,boutfun,xinfun,xoutfun
      COMPLEX(r8), DIMENSION(lmpert,mpert) :: coordmat

      REAL(r8), DIMENSION(:,:), POINTER :: dcosmn,dsinmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: rawmn
      
      INTEGER :: i_id,m_id,t_id,r_id,z_id,rn_id,zn_id,p_id,
     $    x_id,xx_id,xm_id,xxm_id,bm_id,bxm_id,b_id,bx_id
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
            IF ((-mmin+i>=1).AND.(-mmin+i<=lmpert)) THEN
               cawmn(-mmin+i+1)=cosmn(i)+ifac*sinmn(i)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     convert coordinates.
c-----------------------------------------------------------------------
      IF (data_flag .OR. harmonic_flag) THEN
         CALL peq_fcoords(psilim,cawmn,lmfac,lmpert,
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
            CALL peq_weight(psilim,binmn,mfac,mpert,5)
            binmn=twopi*ifac*chi1*(mfac-nn*qlim)*binmn
            CALL peq_weight(psilim,binmn,mfac,mpert,0)
         ENDIF 
         binmn=binmn*scale
         tempmn=binmn
         CALL peq_weight(psilim,tempmn,mfac,mpert,1)
      ENDIF
      finmn=finmn+tempmn
c-----------------------------------------------------------------------
c     filter external flux
c-----------------------------------------------------------------------
      CALL control_filter(finmn,foutmn,filter_types,filter_modes,
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
      CALL peq_weight(psilim,binmn,mfac,mpert,0)
      CALL peq_weight(psilim,boutmn,mfac,mpert,0)
      xinmn=finmn/(chi1*twopi*ifac*(mfac-nn*qlim))
      xoutmn=xspmn
      CALL peq_weight(psilim,xinmn,mfac,mpert,4)
      CALL peq_weight(psilim,xoutmn,mfac,mpert,4)
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
      CALL ascii_open(out_unit,"gpec_control_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_CONTROL: "//
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

      IF ((jac_in /= jac_type).OR.(tin==0).OR.(jin==1)) THEN
         cinmn=0
         coutmn=0
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               cinmn(mlow-lmlow+i)=binmn(i)
               coutmn(mlow-lmlow+i)=boutmn(i)
            ENDIF
         ENDDO         
         acinmn = cinmn
         CALL peq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin) 
         CALL peq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,jin)
         CALL peq_bcoords(psilim,acinmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,1)      
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_in = ",jac_in,
     $        "jsurf_in =",jin,"tmag_in =",tin
         WRITE(out_unit,*)             
         WRITE(out_unit,'(1x,a4,6(1x,a16))')"m","real(bin)","imag(bin)",
     $     "real(bout)","imag(bout)","real(Phi^x)","imag(Phi^x)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,6(es17.8e3))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i)),
     $           REAL(acinmn(i)),AIMAG(acinmn(i))
         ENDDO 
         WRITE(out_unit,*)       
      ENDIF

      IF ((jac_out /= jac_type).OR.(tout==0)) THEN
         cinmn=0
         coutmn=0
         DO i=1,mpert
            IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
               cinmn(mlow-lmlow+i)=binmn(i)
               coutmn(mlow-lmlow+i)=boutmn(i)
            ENDIF
         ENDDO  
         jout=0
         acinmn = cinmn
         CALL peq_bcoords(psilim,cinmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout) 
         CALL peq_bcoords(psilim,coutmn,lmfac,lmpert,
     $        rout,bpout,bout,rcout,tout,jout)
         CALL peq_bcoords(psilim,acinmn,lmfac,lmpert,
     $        rin,bpin,bin,rcin,tin,1)
         WRITE(out_unit,'(1x,a13,a8,1x,2(a12,I2))')"jac_out = ",jac_out,
     $        "jsurf_out =",jout,"tmag_out =",tout  
         WRITE(out_unit,*)          
         WRITE(out_unit,'(1x,a4,6(1x,a16))')"m","real(bin)","imag(bin)",
     $        "real(bout)","imag(bout)","real(Phi^x)","imag(Phi^x)"
         DO i=1,lmpert
            WRITE(out_unit,'(1x,I4,6(es17.8e3))')lmfac(i),
     $           REAL(cinmn(i)),AIMAG(cinmn(i)),
     $           REAL(coutmn(i)),AIMAG(coutmn(i)),
     $           REAL(acinmn(i)),AIMAG(acinmn(i))
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
        IF((jac_out /= jac_type).OR.(tout==0)) CALL peq_bcoords(
     $     psilim,templ,lmfac,lmpert,rout,bpout,bout,rcout,tout,0)
        coordmat(:,i) = templ
      ENDDO

      IF(debug_flag) PRINT *,"Opening "//TRIM(mncfile)
      CALL check( nf90_open(mncfile,nf90_write,mncid) )
      IF(debug_flag) PRINT *,"  Inquiring about dimensions"
      CALL check( nf90_inq_dimid(mncid,"i",i_id) )
      CALL check( nf90_inq_dimid(mncid,"m",m_id) )
      IF(debug_flag) PRINT *,"  Defining variables"
      CALL check( nf90_redef(mncid))
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_vacuum",vengy) )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_surface",sengy) )
      CALL check( nf90_put_att(mncid,nf90_global,
     $            "energy_plasma",pengy) )
      CALL check( nf90_def_var(mncid,"b_nm",nf90_double,
     $            (/m_id,i_id/),bm_id) )
      CALL check( nf90_put_att(mncid,bm_id,"units","Tesla") )
      CALL check( nf90_put_att(mncid,bm_id,"long_name",
     $            "Normal Field") )
      CALL check( nf90_def_var(mncid,"b_xnm",nf90_double,
     $            (/m_id,i_id/),bxm_id) )
      CALL check( nf90_put_att(mncid,bxm_id,"units","Tesla") )
      CALL check( nf90_put_att(mncid,bxm_id,"long_name",
     $            "Externally Applied Normal Field") )
      CALL check( nf90_def_var(mncid,"xi_nm",nf90_double,
     $            (/m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(mncid,xm_id,"units","m") )
      CALL check( nf90_put_att(mncid,xm_id,"long_name",
     $            "Normal Displacement") )
      CALL check( nf90_def_var(mncid,"xi_xnm",nf90_double,
     $            (/m_id,i_id/),xxm_id) )
      CALL check( nf90_put_att(mncid,xxm_id,"units","m") )
      CALL check( nf90_put_att(mncid,xxm_id,"long_name",
     $            "Externally Applied Normal Displacement") )
      CALL check( nf90_enddef(mncid) )
      templ = MATMUL(coordmat,binmn)
      CALL check( nf90_put_var(mncid,bxm_id,RESHAPE((/REAL(templ),
     $          AIMAG(templ)/),(/lmpert,2/))) )
      templ = MATMUL(coordmat,boutmn)
      CALL check( nf90_put_var(mncid,bm_id,RESHAPE((/REAL(templ),
     $          AIMAG(templ)/),(/lmpert,2/))) )
      templ = MATMUL(coordmat,xinmn)
      CALL check( nf90_put_var(mncid,xm_id,RESHAPE((/REAL(templ),
     $          AIMAG(templ)/),(/lmpert,2/))) )
      templ = MATMUL(coordmat,xoutmn)
      CALL check( nf90_put_var(mncid,xxm_id,RESHAPE((/REAL(templ),
     $          AIMAG(templ)/),(/lmpert,2/))) )
      CALL check( nf90_close(mncid) )


      IF (fun_flag) THEN
         CALL peq_bcoords(psilim,binmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL peq_bcoords(psilim,boutmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL peq_bcoords(psilim,xinmn,mfac,mpert,
     $        power_r,power_bp,power_b,0,0,0)
         CALL peq_bcoords(psilim,xoutmn,mfac,mpert,
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

         CALL ascii_open(out_unit,"gpec_control_fun_n"//
     $     TRIM(sn)//".out","UNKNOWN")
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
         CALL check( nf90_def_var(mncid, "xi_n", nf90_double,
     $                    (/t_id,i_id/),x_id) )
         CALL check( nf90_put_att(mncid,x_id,"long_name",
     $               "Displacement") )
         CALL check( nf90_put_att(mncid,x_id,"units","m") )
         CALL check( nf90_def_var(mncid, "xi_xn", nf90_double,
     $                    (/t_id,i_id/),xx_id) )
         CALL check( nf90_put_att(mncid,xx_id,"long_name",
     $               "Externally Applied Displacement") )
         CALL check( nf90_put_att(mncid,xx_id,"units","m") )
         CALL check( nf90_def_var(mncid, "b_n", nf90_double,
     $                    (/t_id,i_id/),b_id) )
         CALL check( nf90_put_att(mncid,b_id,"long_name",
     $               "Field") )
         CALL check( nf90_put_att(mncid,x_id,"units","Tesla") )
         CALL check( nf90_def_var(mncid, "b_xn", nf90_double,
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE control
c-----------------------------------------------------------------------
c     subprogram 4. singfld.
c     compute current and field on rational surfaces.
c-----------------------------------------------------------------------
      SUBROUTINE singfld(egnum,xspmn,spot,
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
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing total resonant fields"
      CALL peq_alloc
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
         CALL peq_sol(lpsi)
         lbwp1mn=bwp1_mn(resnum(ising))
         
         rpsi=respsi+spot/(nn*ABS(singtype(ising)%q1))
         CALL peq_sol(rpsi)
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
            CALL peq_weight(respsi,singbno_mn(:,ising),mfac,mpert,0)
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
      CALL peq_dealloc
c-----------------------------------------------------------------------
c     write results.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"gpec_singfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_SINGFLD: "//
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE singfld
c-----------------------------------------------------------------------
c     subprogram 5. vsingfld.
c     compute resonant field by coils.
c-----------------------------------------------------------------------
      SUBROUTINE vsingfld(rout,bpout,bout,rcout,tout)
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
      CALL ascii_open(out_unit,"gpec_vsingfld_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_SINGFLD: "//
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE vsingfld
c-----------------------------------------------------------------------
c     subprogram 14. control_filter.
c     Filter control surface flux vector in flux bases with energy norms
c-----------------------------------------------------------------------
      SUBROUTINE control_filter(finmn,foutmn,ftypes,fmodes,
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
     $   mx_id,mw_id,mr_id,mp_id,
     $   we_id,re_id,pe_id,se_id,
     $   w_id,r_id,p_id,s_id, wr_id,wp_id,rp_id,ws_id,rs_id,ps_id,
     $   ft_id,fx_id,wx_id,rx_id,px_id,sx_id,wa_id,ra_id,rl_id,
     $   x_id,xe_id,xt_id,wf_id,rf_id,sf_id,sm_id,
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
      COMPLEX(r8), DIMENSION(msing,mpert) :: matsm
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

      IF(timeit) CALL gpec_timer(-2)
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
         CALL peq_weight(psilim,temp,mfac,mpert,2) ! A^1/2
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
      ! start with GPEC flux matrix
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
         DO j=1,msing
            DO i=1,mpert
              IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
                 matsm(j,i) = singcoup1mat(j,mlow-lmlow+i)
              ENDIF
            ENDDO
            CALL peq_fcoords(psilim,matsm(j,:),mfac,mpert,
     $           rout,bpout,bout,rcout,tout,2) ! removes sqrt(A) wieght
            CALL peq_weight(psilim,matsm(j,:),mfac,mpert,2) ! field to sqrt(A)b
         ENDDO
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
            CALL peq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
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
               CALL peq_weight(psilim,filmn,mfac,mpert,2) ! sqrt(A)b to flux
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
           IF((jac_out /= jac_type).OR.(tout==0)) CALL peq_bcoords(
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
         CALL check( nf90_inq_dimid(mncid,"mode_SC",sdid) )
         CALL check( nf90_inq_dimid(mncid,"theta",tdid) )

         ! Start definitions
         CALL check( nf90_redef(mncid))

         IF(debug_flag) PRINT *,"  Defining vecs"
         CALL check( nf90_def_dim(mncid,"mode_XT",mpert,   xdid) )
         CALL check( nf90_def_var(mncid,"mode_XT",nf90_int,xdid,mx_id))
         CALL check( nf90_put_att(mncid, mx_id ,"long_name",
     $    "Total displacement energy eigenmode index") )
         CALL check( nf90_def_var(mncid,"X_EDT",nf90_double,
     $               (/mdid,xdid,idid/),x_id) )
         CALL check( nf90_put_att(mncid,x_id,"long_name",
     $    "Total displacement energy eigendecomposition") )
         CALL check( nf90_def_var(mncid,"X_EVT",nf90_double,
     $               (/xdid/),xe_id) )
         CALL check( nf90_put_att(mncid,xe_id,"units","J/(Wb/m)^2") )
         CALL check( nf90_put_att(mncid,xe_id,"long_name",
     $    "Total displacement energy eigenvalues") )

         CALL check( nf90_def_dim(mncid,"mode_WX",mpert,   wdid) )
         CALL check( nf90_def_var(mncid,"mode_WX",nf90_int,wdid,mw_id))
         CALL check( nf90_put_att(mncid, mw_id ,"long_name",
     $    "Energy-norm external flux energy eigenmode index") )
         CALL check( nf90_def_var(mncid,"W_EDX",nf90_double,
     $               (/mdid,wdid,idid/),w_id) )
         CALL check( nf90_put_att(mncid,w_id,"long_name",
     $    "Energy-norm external flux energy eigendecomposition") )
         CALL check( nf90_def_var(mncid,"W_EVX",nf90_double,
     $               (/wdid/),we_id) )
         CALL check( nf90_put_att(mncid,we_id,"units","J/(Wb/m)^2") )
         CALL check( nf90_put_att(mncid,we_id,"long_name",
     $    "Energy-norm external flux energy eigenvalues") )
         CALL check( nf90_def_var(mncid,"W_EVX_A",nf90_double,
     $                         (/wdid/),wa_id) )
         CALL check( nf90_put_att(mncid,wa_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode amplifications") )
         CALL check( nf90_def_var(mncid,"W_EVX_energyv",nf90_double,
     $                         (/wdid/),wev_id) )
         CALL check( nf90_put_att(mncid,wev_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode vacuum energy") )
         CALL check( nf90_def_var(mncid,"W_EVX_energys",nf90_double,
     $                         (/wdid/),wes_id) )
         CALL check( nf90_put_att(mncid,wes_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode surface energy") )
         CALL check( nf90_def_var(mncid,"W_EVX_energyp",nf90_double,
     $                         (/wdid/),wep_id) )
         CALL check( nf90_put_att(mncid,wep_id,"long_name",
     $    "Energy-norm ex. flux energy eigenmode total energy") )

         CALL check( nf90_def_dim(mncid,"mode_RX",mpert,   rdid) )
         CALL check( nf90_def_var(mncid,"mode_RX",nf90_int,rdid,mr_id))
         CALL check( nf90_put_att(mncid, mr_id ,"long_name",
     $    "Energy-norm external flux reluctance eigenmode index"))
         CALL check( nf90_def_var(mncid,"R_EDX",nf90_double,
     $               (/mdid,rdid,idid/),r_id) )
         CALL check( nf90_put_att(mncid,r_id,"long_name",
     $    "Energy-norm external flux reluctance eigendecomposition") )
         CALL check( nf90_def_var(mncid,"R_EVX",nf90_double,
     $                         (/rdid/),re_id) )
         CALL check( nf90_put_att(mncid,re_id,"long_name",
     $    "Energy-norm external flux reluctance eigenvalues") )
         CALL check( nf90_put_att(mncid,re_id,"units","A/(Wb/m)") )
         CALL check( nf90_def_var(mncid,"R_EVX_RL",nf90_double,
     $                         (/rdid/),rl_id) )
         CALL check( nf90_put_att(mncid,rl_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode RL-normalized") )
         CALL check( nf90_def_var(mncid,"R_EVX_energyv",nf90_double,
     $                         (/rdid/),rev_id) )
         CALL check( nf90_put_att(mncid,rev_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode vacuum energy") )
         CALL check( nf90_def_var(mncid,"R_EVX_energys",nf90_double,
     $                         (/rdid/),res_id) )
         CALL check( nf90_put_att(mncid,res_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode surface energy") )
         CALL check( nf90_def_var(mncid,"R_EVX_energyp",nf90_double,
     $                         (/rdid/),rep_id) )
         CALL check( nf90_put_att(mncid,rep_id,"long_name",
     $    "Energy-norm ex. flux reluctance eigenmode total energy") )

         CALL check( nf90_def_dim(mncid,"mode_PX",mpert,   pdid) )
         CALL check( nf90_def_var(mncid,"mode_PX",nf90_int,pdid,mp_id))
         CALL check( nf90_put_att(mncid, mp_id ,"long_name",
     $    "Energy-norm external flux permeability eigenmode index") )
         CALL check( nf90_def_var(mncid,"P_EDX",nf90_double,
     $               (/mdid,pdid,idid/),p_id) )
         CALL check( nf90_put_att(mncid,p_id,"long_name",
     $    "Energy-norm external flux permeability eigendecomposition") )
         CALL check( nf90_def_var(mncid,"P_EVX",nf90_double,
     $               (/pdid/),pe_id) )
         CALL check( nf90_put_att(mncid,pe_id,"long_name",
     $    "Energy-norm external flux permeability eigenvalues") )

         CALL check( nf90_def_var(mncid,"O_WRX",nf90_double,
     $               (/wdid,rdid/),wr_id) )
         CALL check( nf90_put_att(mncid,wr_id,"long_name",
     $    "Overlap of energy and reluctance eigendecompositions") )
         CALL check( nf90_def_var(mncid,"O_WPX",nf90_double,
     $               (/wdid,pdid/),wp_id) )
         CALL check( nf90_put_att(mncid,wp_id,"long_name",
     $    "Overlap of energy and permeability eigendecompositions") )
         CALL check( nf90_def_var(mncid,"O_RPX",nf90_double,
     $               (/rdid,pdid/),rp_id) )
         CALL check( nf90_put_att(mncid,rp_id,"long_name",
     $    "Overlap of reluctance and permeability eigendecompositions"))

         CALL check( nf90_def_var(mncid,"O_XT",nf90_double,
     $               (/xdid,idid/),xt_id) )
         CALL check( nf90_put_att(mncid,xt_id,"long_name","Total "//
     $    "displacement decomposed in energy eigenmodes") )
         CALL check( nf90_def_var(mncid,"O_WX",nf90_double,
     $               (/xdid,idid/),wx_id) )
         CALL check( nf90_put_att(mncid,wx_id,"long_name","Energy "//
     $    "normalized external flux decomposed in energy eigenmodes") )
         CALL check( nf90_def_var(mncid,"O_RX",nf90_double,
     $               (/rdid,idid/),rx_id) )
         CALL check( nf90_put_att(mncid,rx_id,"long_name","Energy no"//
     $    "rmalized external flux decomposed in reluctance eigenmodes"))
         CALL check( nf90_def_var(mncid,"O_PX",nf90_double,
     $               (/pdid,idid/),px_id) )
         CALL check( nf90_put_att(mncid,px_id,"long_name","Energy no"//
     $  "rmalized external flux decomposed in permeability eigenmodes"))

         IF(singcoup_set)THEN
            CALL check( nf90_def_var(mncid,"C",nf90_double,
     $                  (/mdid,sdid,idid/),sm_id) )
            CALL check( nf90_put_att(mncid,sm_id,"long_name",
     $       "Energy normalized external flux singular-coupling "//
     $       "matrix") )
            CALL check( nf90_def_var(mncid,"C_EDX",nf90_double,
     $                  (/mdid,sdid,idid/),s_id) )
            CALL check( nf90_put_att(mncid,s_id,"long_name",
     $       "Energy normalized external flux singular-coupling "//
     $       "SVD right-singular vectors") )
            CALL check( nf90_def_var(mncid,"C_EVX",nf90_double,
     $                            (/sdid/),se_id) )
            CALL check( nf90_put_att(mncid,se_id,"long_name",
     $       "Energy normalized external flux singular-coupling "//
     $       "SVD singular values") )
            CALL check( nf90_put_att(mncid,se_id,"units","unitless") )
            CALL check( nf90_def_var(mncid,"C_EVX_energyv",nf90_double,
     $                            (/sdid/),sev_id) )
            CALL check( nf90_put_att(mncid,sev_id,"long_name",
     $       "Singular-coupling eigenmode vacuum energy") )
            CALL check( nf90_def_var(mncid,"C_EVX_energys",nf90_double,
     $                            (/sdid/),ses_id) )
            CALL check( nf90_put_att(mncid,ses_id,"long_name",
     $       "Singular-coupling eigenmode surface energy") )
            CALL check( nf90_def_var(mncid,"C_EVX_energyp",nf90_double,
     $                            (/sdid/),sep_id) )
            CALL check( nf90_put_att(mncid,sep_id,"long_name",
     $       "Singular-coupling eigenmode total energy") )
            CALL check( nf90_def_var(mncid,"O_WCX",nf90_double,
     $                       (/wdid,sdid/),ws_id) )
            CALL check( nf90_put_att(mncid,ws_id,"long_name",
     $       "Overlap of energy and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_RCX",nf90_double,
     $                       (/rdid,sdid/),rs_id) )
            CALL check( nf90_put_att(mncid,rs_id,"long_name",
     $       "Overlap of reluctance and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_PCX",nf90_double,
     $                       (/pdid,sdid/),ps_id) )
            CALL check( nf90_put_att(mncid,ps_id,"long_name",
     $       "Overlap of permeability and singular-coupling modes") )
            CALL check( nf90_def_var(mncid,"O_CX",nf90_double,
     $                       (/sdid,idid/),sx_id) )
            CALL check( nf90_put_att(mncid,sx_id,"long_name",
     $       "Energy normalized external flux decomposed in singular "//
     $       "coupling modes"))

         ENDIF

         CALL check( nf90_def_var(mncid,"Phi_EX",nf90_double,
     $               (/mdid,idid/),fx_id) )
         CALL check( nf90_put_att(mncid,fx_id,"units","Wb/m") )
         CALL check( nf90_put_att(mncid,fx_id,"long_name",
     $    "Energy-norm external flux") )
         CALL check( nf90_def_var(mncid,"Phi_ET",nf90_double,
     $               (/mdid,idid/),ft_id) )
         CALL check( nf90_put_att(mncid,ft_id,"units","Wb/m") )
         CALL check( nf90_put_att(mncid,ft_id,"long_name",
     $    "Energy-norm total flux") )

         IF(fun_flag)THEN
            CALL check( nf90_def_var(mncid,"W_EDX_FUN",nf90_double,
     $               (/tdid,wdid,idid/),wf_id) )
            CALL check( nf90_put_att(mncid,wf_id,"long_name",
     $         "Energy-norm external flux energy eigenmodes") )
            CALL check( nf90_def_var(mncid,"R_EDX_FUN",nf90_double,
     $               (/tdid,rdid,idid/),rf_id) )
            CALL check( nf90_put_att(mncid,rf_id,"long_name",
     $         "Energy-norm external flux reluctance eigenmodes") )
            IF(singcoup_set)THEN
               CALL check( nf90_def_var(mncid,"C_EDX_FUN",nf90_double,
     $                  (/tdid,sdid,idid/),sf_id) )
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
         CALL check( nf90_put_var(mncid,x_id,RESHAPE((/REAL(xveco),
     $               AIMAG(xveco)/),(/lmpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,xe_id,xvals) )
         CALL check( nf90_put_var(mncid,w_id,RESHAPE((/REAL(wveco),
     $               AIMAG(wveco)/),(/lmpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,we_id,wvals) )
         CALL check( nf90_put_var(mncid,wa_id,avals) )
         CALL check( nf90_put_var(mncid,r_id,RESHAPE((/REAL(rveco),
     $               AIMAG(rveco)/),(/lmpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,re_id,rvals) )
         CALL check( nf90_put_var(mncid,p_id,RESHAPE((/REAL(pveco),
     $               AIMAG(pveco)/),(/lmpert,mpert,2/))) )
         CALL check( nf90_put_var(mncid,pe_id,pvals) )
         IF(singcoup_set) THEN
            CALL check( nf90_put_var(mncid,sm_id,RESHAPE( (/REAL(
     $       singcoup1mat),AIMAG(singcoup1mat)/),(/lmpert,msing,2/))))
            CALL check( nf90_put_var(mncid,s_id,RESHAPE(
     $       (/REAL(sveco),AIMAG(sveco)/),(/lmpert,msing,2/))))
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
         matms = MATMUL(CONJG(TRANSPOSE(rvecs)),pvecs)
         CALL check( nf90_put_var(mncid,rp_id,ABS(matms)) )
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
         CALL peq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
         templ = MATMUL(coordmat,temp) ! output coords
         CALL check( nf90_put_var(mncid,ft_id,RESHAPE((/REAL(templ),
     $               AIMAG(templ)/),(/lmpert,2/))) )
         temp = finmn
         CALL peq_weight(psilim,temp,mfac,mpert,6) ! flux to sqrt(A)b
         templ = MATMUL(coordmat,temp) ! output coords
         CALL check( nf90_put_var(mncid,fx_id,RESHAPE((/REAL(templ),
     $               AIMAG(templ)/),(/lmpert,2/))) )

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

         ! Eigenmodes in real space (R,z) written in control
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE control_filter
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
c     subprogram 16. init_netcdf.
c     Initialize the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE init_netcdf
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
      CALL check( nf90_def_dim(mncid,"m",lmpert,  mmdid) )
      CALL check( nf90_def_var(mncid,"m",nf90_int,mmdid,mmvid) )
      CALL check( nf90_def_dim(mncid,"mode",mpert,   medid) )
      CALL check( nf90_def_var(mncid,"mode",nf90_int,medid,mevid))
      CALL check( nf90_def_dim(mncid,"mode_SC",msing,   msdid) )
      CALL check( nf90_def_var(mncid,"mode_SC",nf90_int,msdid,msvid))
      CALL check( nf90_def_dim(mncid,"theta",mthsurf+1,  mtdid) )
      CALL check( nf90_def_var(mncid,"theta",nf90_double,mtdid,mtvid) )
      CALL check( nf90_put_att(mncid,nf90_global,"Jacobian",jac_out))
      CALL check( nf90_put_att(mncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(mncid,nf90_global,"psi_n_lim",psilim))

      IF(debug_flag) PRINT *," - Defining flux netcdf globals"
      CALL check( nf90_put_att(fncid,nf90_global,"title",
     $     "GPEC outputs in magnetic coordinate systems"))
      CALL check( nf90_def_dim(fncid,"i",2,       fidid) )
      CALL check( nf90_def_var(fncid,"i",nf90_int,fidid,fivid) )
      CALL check( nf90_def_dim(fncid,"m",lmpert,  fmdid) )
      CALL check( nf90_def_var(fncid,"m",NF90_INT,fmdid,fmvid) )
      CALL check( nf90_def_dim(fncid,"psi_n",mstep,    fpdid) )
      CALL check( nf90_def_var(fncid,"psi_n",nf90_double,fpdid,fpvid) )
      CALL check( nf90_def_dim(fncid,"theta",mthsurf+1,  ftdid) )
      CALL check( nf90_def_var(fncid,"theta",nf90_double,ftdid,ftvid) )
      CALL check( nf90_put_att(fncid,nf90_global,"Jacobian",jac_type))
      CALL check( nf90_put_att(fncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(fncid,nf90_global,"psi_n_lim",psilim))

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
      CALL check( nf90_put_var(mncid,mmvid,lmfac) )
      CALL check( nf90_put_var(mncid,mevid,mmodes) )
      CALL check( nf90_put_var(mncid,msvid,mmodes(:msing)) )
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
      END SUBROUTINE init_netcdf
c-----------------------------------------------------------------------
c     subprogram 17. close_netcdf.
c     Close the netcdf files used for module outputs.
c-----------------------------------------------------------------------
      SUBROUTINE close_netcdf
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
      END SUBROUTINE close_netcdf

      END MODULE output_control
