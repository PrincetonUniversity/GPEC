c-----------------------------------------------------------------------
c     GENERALIZED PERTURBED EQUILIBRIUM CODE
c     write various output results of gpec routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. output_profile
c      1. dw_profile
c      2. dw_matrix
c      3. pmodb
c      4. xbnormal
c      5. vbnormal
c      6. xtangent
c-----------------------------------------------------------------------
c     subprogram 0. output_profile.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE output_profile
      USE gpec_global
      USE local_mod, ONLY: ascii_open, ascii_close, bin_open, bin_close
      USE spline_mod, ONLY: spline_type, spline_alloc, spline_eval,
     $    spline_fit, spline_int
      USE cspline_mod, ONLY: cspline_alloc, cspline_eval, cspline_fit,
     $    cspline_int
      USE bicube_mod, ONLY: bicube_eval
      USE gpec_math, ONLY : iscdftf,iscdftb,issurfint
      USE peq, ONLY: peq_sol, peq_contra, peq_cova, peq_alloc,
     $    peq_dealloc, peq_bcoords, peq_fcoords, peq_bcoordsout,
     $    peq_weight
      USE dcon_interface, ONLY: dcon_build
      USE coil_mod, ONLY : helicity, cmlow, cmhigh, cmpert, cmpsi,
     $    machine, coil_read, coil_dealloc
      USE field_mod, ONLY: field_bs_psi
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

c-----------------------------------------------------------------------
c     subprogram 6. w_profile.
c     restore energy and torque profiles from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE dw_profile(egnum,xspmn)
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
      CALL dcon_build(egnum,xspmn)
      CALL cspline_alloc(dwk,mstep,1)
      dwk%xs=psifac
      DO istep=0,mstep
         dwk%fs(istep,1)=2*nn*ifac*
     $        SUM(CONJG(u1%fs(istep,:))*u2%fs(istep,:))/(2.0*mu0)
      ENDDO
      CALL cspline_fit(dwk,"extrap")

      WRITE(*,*)"Restoring energy and torque profiles"

      CALL ascii_open(out_unit,"gpec_dw_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_dw: "//
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
      END SUBROUTINE dw_profile
c-----------------------------------------------------------------------
c     subprogram 7. dw_matrix.
c     restore energy and torque response matrix from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE dw_matrix
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: ipert,istep,itheta,i,j,mval,iindex,lwork
      REAL(r8) :: jarea,ileft
      COMPLEX(r8) :: t1,t2
      LOGICAL :: diagnose1=.FALSE.,diagnose2=.TRUE.
      LOGICAL :: kstar_coil_optimization=.TRUE.

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
         CALL peq_weight(psilim,fldflxmn,mfac,mpert,2)
         fldflxmat(:,ipert)=fldflxmn*sqrt(jarea)
      ENDDO
      WRITE(*,*)
     $     "Call solutions for general response matrix functions"
      edge_flag=.TRUE.
      DO ipert=1,mpert
         edge_mn=0
         edge_mn(ipert)=1.0
         CALL dcon_build(0,edge_mn)
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
         CALL ascii_open(out_unit,"gpec_dw_diag_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_dw_diag: "//
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
      CALL ascii_open(out_unit,"gpec_dw_matrix_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_dw_matrix: "//
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
      CALL ascii_open(out_unit,"gpec_dw_eigen_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_dw_eigen: "//
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
         CALL ascii_open(out_unit,"gpec_dw_minimum_torque_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_dw_diag: "//
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

      CALL ascii_open(out_unit,"gpec_dw_ratio_eigen_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_dw_ratio_eigen: "//
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
         CALL ascii_open(out_unit,"gpec_dw_optimized_torque_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_dw_diag: "//
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
         CALL ascii_open(out_unit,"gpec_coil_matrix_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_coil_matrix: "//
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
         CALL ascii_open(out_unit,"gpec_coil_eigen_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_coil_eigen: "//
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

         CALL ascii_open(out_unit,"gpec_coil_ratio_eigen_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_coil_ratio_eigen: "//
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
      END SUBROUTINE dw_matrix
c-----------------------------------------------------------------------
c     subprogram 6. pmodb.
c     compute perturbed mod b.
c-----------------------------------------------------------------------
      SUBROUTINE pmodb(egnum,xspmn,rout,bpout,bout,rcout,tout)
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
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing |b| and del.x_prp"
      CALL cspline_alloc(cspl,mthsurf,2)
      cspl%xs=theta

      CALL dcon_build(egnum,xspmn)
      CALL peq_alloc

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
         CALL peq_sol(psifac(istep))
         CALL peq_contra(psifac(istep))
         CALL peq_cova(psifac(istep))
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
         CALL peq_bcoordsout(eulbparmout(i,:),eulbpar_mn,psi,tout,0)
         CALL peq_bcoordsout(lagbparmout(i,:),eulbpar_mn,psi,tout,0)
         CALL peq_bcoordsout(divxprpmout(i,:),eulbpar_mn,psi,tout,0)
         CALL peq_bcoordsout(curvmout(i,:)   ,eulbpar_mn,psi,tout,0)
c         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
c            CALL peq_bcoords(psifac(istep),eulbpar_mn,
c     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
c            CALL peq_bcoords(psifac(istep),lagbpar_mn,
c     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
c            CALL peq_bcoords(psifac(istep),divxprp_mn,
c     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
c            CALL peq_bcoords(psifac(istep),curv_mn,
c     $           mfac,mpert,rout,bpout,bout,rcout,tout,0)
c         ENDIF
c         eulbparmout(istep,:)=eulbpar_mn
c         lagbparmout(istep,:)=lagbpar_mn
c         divxprpmout(istep,:)=divxprp_mn
c         curvmout(istep,:)  =curv_mn
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
      CALL check( nf90_def_var(fncid, "b_eul", nf90_double,
     $            (/p_id,t_id,i_id/),be_id) )
      CALL check( nf90_put_att(fncid,be_id,"long_name",
     $            "Eulerian Perturbed Field") )
      CALL check( nf90_put_att(fncid,be_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_lag", nf90_double,
     $            (/p_id,t_id,i_id/),bl_id) )
      CALL check( nf90_put_att(fncid,bl_id,"long_name",
     $            "Lagrangian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bl_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_m-eul", nf90_double,
     $            (/p_id,m_id,i_id/),bme_id) )
      CALL check( nf90_put_att(fncid,bme_id,"long_name",
     $            "Eulerian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bme_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_m-lag", nf90_double,
     $            (/p_id,m_id,i_id/),bml_id) )
      CALL check( nf90_put_att(fncid,bml_id,"long_name",
     $            "Lagrangian Perturbed Field") )
      CALL check( nf90_put_att(fncid,bml_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "Bdivxi_perp", nf90_double,
     $            (/p_id,t_id,i_id/),x_id) )
      CALL check( nf90_put_att(fncid,x_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_def_var(fncid, "Bdivxi_m-perp", nf90_double,
     $            (/p_id,m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_def_var(fncid, "Bkappaxi_perp", nf90_double,
     $            (/p_id,t_id,i_id/),k_id) )
      CALL check( nf90_put_att(fncid,k_id,"long_name",
     $            "Divergence of the normal displacement") )
      CALL check( nf90_def_var(fncid, "Bkappaxi_m-perp", nf90_double,
     $            (/p_id,m_id,i_id/),km_id) )
      CALL check( nf90_put_att(fncid,km_id,"long_name",
     $            "Divergence of the normal displacement") )
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


      CALL ascii_open(out_unit,"gpec_pmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_PMODB: "//
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

         CALL ascii_open(out_unit,"gpec_pmodb_chebyshev_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PMODB_CHEBYSHEV: "//
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

         CALL ascii_open(out_unit,"gpec_chepmodb_n"//
     $     TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_CHEPMODB: "//
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
            CALL peq_sol(psifac(istep))
            CALL peq_contra(psifac(istep))
            CALL peq_cova(psifac(istep))
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

      CALL peq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE pmodb
c-----------------------------------------------------------------------
c     subprogram 7. xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE xbnormal(egnum,xspmn,rout,bpout,bout,rcout,tout)
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
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b normal components"

      CALL dcon_build(egnum,xspmn)

      CALL peq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL peq_sol(psifac(istep))
         CALL peq_contra(psifac(istep))

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
         bwp_mn=bwp_mn/area
         IF (bwp_pest_flag) THEN
            pwpmn=0
            DO i=1,mpert
               IF ((mlow-lmlow+i>=1).AND.(mlow-lmlow+i<=lmpert)) THEN
                  pwpmn(mlow-lmlow+i)=bno_mn(i)
               ENDIF
            ENDDO
            CALL peq_bcoords(psifac(istep),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)
            pwpmns(istep,:)=pwpmn
         ENDIF            

         CALL peq_bcoordsout(xnomns(istep,:),xno_mn,psifac(istep))
         CALL peq_bcoordsout(bnomns(istep,:),xno_mn,psifac(istep))
         CALL peq_bcoordsout(bwpmns(istep,:),xno_mn,psifac(istep),ji=1)
c         IF ((jac_out /= jac_type).OR.(tout==0)) THEN
c            bwp_mn=bno_mn
c            CALL peq_bcoords(psifac(istep),xno_mn,mfac,mpert,
c     $           rout,bpout,bout,rcout,tout,0)
c            CALL peq_bcoords(psifac(istep),bno_mn,mfac,mpert,
c     $           rout,bpout,bout,rcout,tout,0)
c            CALL peq_bcoords(psifac(istep),bwp_mn,mfac,mpert,
c     $           rout,bpout,bout,rcout,tout,1)
c         ENDIF
c         xnomns(istep,:)=xno_mn
c         bnomns(istep,:)=bno_mn
c         bwpmns(istep,:)=bwp_mn
         xnofuns(istep,:)=xnofuns(istep,:)*EXP(ifac*nn*dphi)
         bnofuns(istep,:)=bnofuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL peq_dealloc

      CALL ascii_open(out_unit,"gpec_xbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_XBNORMAL: "//
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
         CALL ascii_open(out_unit,"gpec_bnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BNORMAL_PEST: "//
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
      CALL check( nf90_def_var(fncid, "b_n", nf90_double,
     $            (/p_id,t_id,i_id/),b_id) )
      CALL check( nf90_put_att(fncid,b_id,"long_name",
     $            "Perturbed Field Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,b_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "b_nm", nf90_double,
     $            (/p_id,m_id,i_id/),bm_id) )
      CALL check( nf90_put_att(fncid,bm_id,"long_name",
     $            "Perturbed Field Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,bm_id,"units","Tesla") )
      CALL check( nf90_def_var(fncid, "xi_n", nf90_double,
     $            (/p_id,t_id,i_id/),x_id) )
      CALL check( nf90_put_att(fncid,x_id,"long_name",
     $            "Displacement Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,x_id,"units","m") )
      CALL check( nf90_def_var(fncid, "xi_nm", nf90_double,
     $            (/p_id,m_id,i_id/),xm_id) )
      CALL check( nf90_put_att(fncid,xm_id,"long_name",
     $            "Displacement Normal to the Flux Surface") )
      CALL check( nf90_put_att(fncid,xm_id,"units","m") )
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
         
         CALL ascii_open(out_unit,"gpec_xbnormal_flux_n"//
     $        TRIM(sn)//".out","UNKNOWN")         
         WRITE(out_unit,*)"GPEC_XBNROMAL_FLUX: "//
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbnormal
c-----------------------------------------------------------------------
c     subprogram 8. vbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE vbnormal(rout,bpout,bout,rcout,tout)
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
      IF(timeit) CALL gpec_timer(-2)
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
            CALL peq_bcoords(psi(ipsi),pwpmn,lmfac,lmpert,
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
            CALL peq_bcoords(psi(ipsi),vwpmn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            vnomns(ipsi,:)=vwpmn
            vwpmn=0
            vwpmn=vwpmns(ipsi,:)
            CALL peq_bcoords(psi(ipsi),vwpmn,mfac,mpert,
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
c      CALL check( nf90_def_var(fncid, "b_m-vperp", nf90_double,
c     $            (/p_id,m_id,i_id/),bm_id) )
c      CALL check( nf90_put_att(fncid,bm_id,"long_name",
c     $            "Vacuum Field Normal to the Flux Surface") )
c      CALL check( nf90_put_att(fncid,bm_id,"units","Tesla") )
c      CALL check( nf90_enddef(fncid) )
c      IF(debug_flag) PRINT *,"  Writting variables"
c      CALL check( nf90_put_var(fncid,bm_id,qs) )
c      CALL check( nf90_put_var(fncid,bm_id,RESHAPE((/REAL(vnomns),
c     $             AIMAG(vnomns)/),(/mstep,mpert,2/))) )
c      CALL check( nf90_close(fncid) )
c      IF(debug_flag) PRINT *,"Closed "//TRIM(fncfile)

      CALL ascii_open(out_unit,"gpec_vbnormal_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_VBNROMAL: "//
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
         CALL ascii_open(out_unit,"gpec_vbnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBNROMAL_PEST: "//
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
      IF(timeit) CALL gpec_timer(-2)
      RETURN
      END SUBROUTINE vbnormal
c-----------------------------------------------------------------------
c     subprogram 9. xbtangent.
c-----------------------------------------------------------------------
      SUBROUTINE xbtangent(egnum,xspmn,rout,bpout,bout,rcout,tout)
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
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing x and b tangential components"

      CALL dcon_build(egnum,xspmn)

      CALL peq_alloc
      DO istep=1,mstep
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% xi and b computations"
         CALL peq_sol(psifac(istep))
         CALL peq_contra(psifac(istep))
         CALL peq_cova(psifac(istep))

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
            CALL peq_bcoords(psifac(istep),xta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
            CALL peq_bcoords(psifac(istep),bta_mn,mfac,mpert,
     $           rout,bpout,bout,rcout,tout,0)
         ENDIF
         xtamns(istep,:)=xta_mn
         btamns(istep,:)=bta_mn
         xtafuns(istep,:)=xtafuns(istep,:)*EXP(ifac*nn*dphi)
         btafuns(istep,:)=btafuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL peq_dealloc

      CALL ascii_open(out_unit,"gpec_xbtangent_n"//
     $     TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_XBNROMAL: "//
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
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbtangent
c-----------------------------------------------------------------------
c     subprogram 13. xclebsch.
c     Write Clebsch coordinate displacements.
c-----------------------------------------------------------------------
      SUBROUTINE xclebsch(egnum,xspmn)
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
      IF(timeit) CALL gpec_timer(-2)
      IF(verbose) WRITE(*,*)"Computing Clebsch displacements"

      CALL dcon_build(egnum,xspmn)
      CALL peq_alloc

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
         CALL peq_sol(psifac(istep))
         CALL peq_contra(psifac(istep))
c-----------------------------------------------------------------------
c     compute mod b variations in hamada.
c-----------------------------------------------------------------------
         xmp1mns(istep,:) = xmp1_mn
         xspmns(istep,:) = xsp_mn
         xmsmns(istep,:) = xms_mn
      ENDDO
         
      CALL ascii_open(out_unit,"gpec_xclebsch_n"//
     $   TRIM(sn)//".out","UNKNOWN")
      WRITE(out_unit,*)"GPEC_PMODB: "//
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
      CALL check( nf90_def_var(fncid, "derxi_m_contrapsi", nf90_double,
     $            (/p_id,m_id,i_id/),dp_id) )
      CALL check( nf90_put_att(fncid,dp_id,"long_name",
     $            "Psi derivative of contravarient psi displacement") )
      CALL check( nf90_def_var(fncid, "xi_m_contrapsi", nf90_double,
     $            (/p_id,m_id,i_id/),xp_id) )
      CALL check( nf90_put_att(fncid,xp_id,"long_name",
     $            "Contravarient psi displacement") )
      CALL check( nf90_def_var(fncid, "xi_m_contraalpha", nf90_double,
     $            (/p_id,m_id,i_id/),xa_id) )
      CALL check( nf90_put_att(fncid,xa_id,"long_name",
     $            "Contravarient clebsch angle displacement") )
      ids = (/xp_id,dp_id,xa_id/)
      DO i=1,3
         CALL check( nf90_put_att(fncid,ids(i),"units","m") )
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


      CALL peq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xclebsch
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

      END MODULE output_profile
