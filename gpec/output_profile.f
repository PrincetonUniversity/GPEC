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
      !USE utilities, ONLY: progressbar
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

      ! public variables
      REAL(r8), DIMENSION(:), ALLOCATABLE :: qs
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rs, zs, rvecs, zvecs

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::
     $   xtamns, xtafuns, btamns, btafuns

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: equilbfun
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::
     $   eulbparmout, lagbparmout, divxprpmout, curvmout,
     $   eulbparfout, lagbparfout, divxprpfout, curvfout

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: xnofuns, bnofuns,
     $   xnomns,bnomns, bwpmns, pwpmns

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: vnomns,vwpmns,vpwpmns

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE ::
     $   xmp1out, xspout, xmsout,
     $   xmp1funs, xspfuns, xmsfuns

      TYPE(cspline_type) :: dwk

      ! private things that might have namespace clashes if used elsewhere
      PRIVATE :: ascii_header

      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. dw_profile.
c     restore energy and torque profiles from solutions.
c-----------------------------------------------------------------------
      SUBROUTINE dw_profile(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep

c-----------------------------------------------------------------------
c     build and spline energy and torque profiles.
c-----------------------------------------------------------------------
      WRITE(*,*)"Restoring energy and torque profiles"

      CALL dcon_build(egnum,xspmn)
      CALL cspline_alloc(dwk,mstep,1)
      dwk%xs=psifac
      DO istep=0,mstep
         dwk%fs(istep,1)=2*nn*ifac*
     $        SUM(CONJG(u1%fs(istep,:))*u2%fs(istep,:))/(2.0*mu0)
      ENDDO
      CALL cspline_fit(dwk,"extrap")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dw_profile
c-----------------------------------------------------------------------
c     subprogram 2. dw_matrix.
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
c     subprogram 3. pmodb.
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
      COMPLEX(r8), DIMENSION(mpert) :: eulbpar_mn,lagbpar_mn,
     $     divxprp_mn,curv_mn
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xsp_fun,xms_fun,
     $     bvt_fun,bvz_fun,xmz_fun,xvt_fun,xvz_fun,xmt_fun,xsp1_fun
      COMPLEX(r8), DIMENSION(mstep,mpert) :: eulbparmns,lagbparmns,
     $     divxprpmns,curvmns,xmp1mns,xspmns,xmsmns,xmtmns,xmzmns
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) ::eulbparfun,lagbparfun,
     $     divxprpfun,curvfun

      REAL(r8), DIMENSION(:), POINTER :: psis,ches,chex
      COMPLEX(r8), DIMENSION(:,:), POINTER ::chea,chelagbmns,chelagbmout
      COMPLEX(r8), DIMENSION(mpert) :: chelagb_mn

      TYPE(cspline_type) :: cspl,chespl
c-----------------------------------------------------------------------
c     allocation of module wide variables.
c-----------------------------------------------------------------------
      ALLOCATE( eulbparmout(mstep,lmpert), lagbparmout(mstep,lmpert),
     $   divxprpmout(mstep,lmpert), curvmout(mstep,lmpert) )
      ALLOCATE( eulbparfout(mstep,0:mthsurf), curvfout(mstep,0:mthsurf),
     $   lagbparfout(mstep,0:mthsurf), divxprpfout(mstep,0:mthsurf),
     $   equilbfun(mstep,0:mthsurf),
     $   rs(mstep,0:mthsurf), zs(mstep,0:mthsurf))
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
         CALL peq_bcoordsout(lagbparmout(i,:),lagbpar_mn,psi,tout,0)
         CALL peq_bcoordsout(divxprpmout(i,:),divxprp_mn,psi,tout,0)
         CALL peq_bcoordsout(curvmout(i,:),curv_mn,psi,tout,0)
         eulbparfout(istep,:)=eulbparfun(istep,:)*EXP(ifac*nn*dphi)
         lagbparfout(istep,:)=lagbparfun(istep,:)*EXP(ifac*nn*dphi)
         divxprpfout(istep,:)=divxprpfun(istep,:)*EXP(ifac*nn*dphi)
         curvfout(istep,:)=curvfun(istep,:)*EXP(ifac*nn*dphi)
      ENDDO

      CALL cspline_dealloc(cspl)
      CALL peq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE pmodb
c-----------------------------------------------------------------------
c     subprogram 4. xbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE xbnormal(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: i,istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax,area

      REAL(r8), DIMENSION(0:mthsurf) :: delpsi,jacs,dphi
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwp_fun,bwp_fun

      COMPLEX(r8), DIMENSION(mstep,lmpert) :: xmns,ymns
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: psis
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss
c-----------------------------------------------------------------------
c     allocation of module wide variables.
c-----------------------------------------------------------------------
      ALLOCATE( xnomns(mstep,lmpert), bnomns(mstep,lmpert),
     $   bwpmns(mstep,lmpert), pwpmns(mstep,lmpert) )
      ALLOCATE( xnofuns(mstep,lmpert), bnofuns(mstep,lmpert) )
      IF ( .NOT. ALLOCATED(rs) )
     $   ALLOCATE(rs(mstep,0:mthsurf), zs(mstep,0:mthsurf))
      IF ( .NOT. ALLOCATED(rvecs) )
     $   ALLOCATE(rvecs(mstep,0:mthsurf), zvecs(mstep,0:mthsurf))
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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbnormal
c-----------------------------------------------------------------------
c     subprogram 5. vbnormal.
c-----------------------------------------------------------------------
      SUBROUTINE vbnormal(rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: rout,bpout,bout,rcout,tout
      
      INTEGER :: ipsi,ipert,i
      REAL(r8), DIMENSION(0:cmpsi) :: psi
      COMPLEX(r8), DIMENSION(:), POINTER :: vcmn

      COMPLEX(r8), DIMENSION(mpert) :: vwpmn
      COMPLEX(r8), DIMENSION(lmpert) :: pwpmn
      REAL(r8), DIMENSION(cmpsi) :: qs
      REAL(r8), DIMENSION(cmpsi,mpert) :: xmns,ymns
c-----------------------------------------------------------------------
c     allocation of module wide variables.
c-----------------------------------------------------------------------
      ALLOCATE( vnomns(cmpsi,mpert), vwpmns(cmpsi,mpert),
     $   vpwpmns(cmpsi,lmpert) )
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
      vpwpmns = 0

      DO ipsi=1,cmpsi
         CALL spline_eval(sq,psi(ipsi),0)
         qs(ipsi)=sq%f(4)
         CALL field_bs_psi(psi(ipsi),vcmn,0)

         IF (bwp_pest_flag) THEN
            DO i=1,cmpert
               IF ((cmlow-lmlow+i>=1).AND.(cmlow-lmlow+i<=lmpert)) THEN
                  vpwpmns(ipsi,cmlow-lmlow+i)=vcmn(i)
               ENDIF
            ENDDO
            pwpmn=0
            pwpmn=vpwpmns(ipsi,:)
            CALL peq_bcoords(psi(ipsi),pwpmn,lmfac,lmpert,
     $           2,0,0,0,0,1)             
            vpwpmns(ipsi,:)=pwpmn
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
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(-2)
      RETURN
      END SUBROUTINE vbnormal
c-----------------------------------------------------------------------
c     subprogram 6. xbtangent.
c-----------------------------------------------------------------------
      SUBROUTINE xbtangent(egnum,xspmn,rout,bpout,bout,rcout,tout)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum,rout,bpout,bout,rcout,tout
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep,ipert,iindex,itheta
      REAL(r8) :: ileft,ximax,rmax

      REAL(r8), DIMENSION(mstep,0:mthsurf) :: vecs
      REAL(r8), DIMENSION(0:mthsurf) :: jacs,bs,dphi
      COMPLEX(r8), DIMENSION(mstep,mpert) :: xmns,ymns
      COMPLEX(r8), DIMENSION(0:mthsurf) :: xwt_fun,xvt_fun,xvz_fun,
     $     bwt_fun,bvt_fun,bvz_fun,xta_fun,bta_fun
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf) :: rss,zss
c-----------------------------------------------------------------------
c     allocation of module wide variables.
c-----------------------------------------------------------------------
      ALLOCATE(xtamns(mstep,lmpert),btamns(mstep,lmpert))
      ALLOCATE(xtafuns(mstep,0:mthsurf),btafuns(mstep,0:mthsurf))
      ALLOCATE(rvecs(mstep,0:mthsurf),zvecs(mstep,0:mthsurf))
      IF ( .NOT. ALLOCATED(rs) )
     $   ALLOCATE(rs(mstep,0:mthsurf),zs(mstep,0:mthsurf))
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

         ! convert to output coordinates
         CALL peq_bcoordsout(xtamns(istep,:),xta_mn,psifac(istep))
         CALL peq_bcoordsout(btamns(istep,:),xta_mn,psifac(istep))
         xtafuns(istep,:)=xtafuns(istep,:)*EXP(ifac*nn*dphi)
         btafuns(istep,:)=btafuns(istep,:)*EXP(ifac*nn*dphi)
      ENDDO
      CALL peq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xbtangent
c-----------------------------------------------------------------------
c     subprogram 7. xclebsch.
c     Write Clebsch coordinate displacements.
c-----------------------------------------------------------------------
      SUBROUTINE xclebsch(egnum,xspmn)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: egnum
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xspmn

      INTEGER :: istep, iindex
      REAL(r8) :: psi, ileft

c-----------------------------------------------------------------------
c     start calculations if not already already calculated
c-----------------------------------------------------------------------
      IF( ALLOCATED(xmp1out) ) RETURN
      IF(verbose) WRITE(*,*)"Computing Clebsch displacements"
      IF(timeit) CALL gpec_timer(-2)
c-----------------------------------------------------------------------
c     allocation.
c-----------------------------------------------------------------------
      ALLOCATE( xmp1out(mstep,lmpert), xspout(mstep,lmpert),
     $    xmsout(mstep,lmpert) )
      ALLOCATE( xmp1funs(mstep,0:mthsurf), xspfuns(mstep,0:mthsurf),
     $    xmsfuns(mstep,0:mthsurf) )
c-----------------------------------------------------------------------
c     compute functions on magnetic surfaces with regulation.
c-----------------------------------------------------------------------
      CALL dcon_build(egnum,xspmn)
      CALL peq_alloc

      DO istep=1,mstep
         !IF(verbose) call progressbar(istep,1,mstep)
         iindex = FLOOR(REAL(istep,8)/FLOOR(mstep/10.0))*10
         ileft = REAL(istep,8)/FLOOR(mstep/10.0)*10-iindex
         IF ((istep-1 /= 0) .AND. (ileft == 0) .AND. verbose)
     $        WRITE(*,'(1x,a9,i3,a23)')
     $        "volume = ",iindex,"% Clebsch decomposition"
         psi = psifac(istep)
         CALL spline_eval(sq,psi,1)
         CALL peq_sol(psi)
         CALL peq_contra(psi) ! not necessary ??
         ! convert to real space
         !IF (fun_flag) THEN
         !   CALL iscdftb(mfac,mpert,xmp1funs(istep,:),mthsurf,xmp1_mn)
         !   CALL iscdftb(mfac,mpert,xspfuns(istep,:),mthsurf,xsp_mn)
         !   CALL iscdftb(mfac,mpert,xmsfuns(istep,:),mthsurf,xms_mn)
         !ENDIF
         ! convert to jac_out
         CALL peq_bcoordsout(xmp1out(istep,:),xmp1_mn,psi)
         CALL peq_bcoordsout(xspout(istep,:),xsp_mn,psi)
         CALL peq_bcoordsout(xmsout(istep,:),xms_mn,psi)
      ENDDO

      CALL peq_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE xclebsch


c-----------------------------------------------------------------------
      SUBROUTINE ascii_header( op_jac )
c-----------------------------------------------------------------------
c     *DESCRIPTION:
c        Write common profile ascii file headers.
c        Having this in one place makes it easier to modify consistently.
c     *ARGUMENTS:
c        op_jac: character(16) optional. Jacobian of file variables
c                (default is jac_out)
c
c-----------------------------------------------------------------------
      CHARACTER(16), INTENT(IN), OPTIONAL :: op_jac
      CHARACTER(16) :: jac

      IF(PRESENT(OP_JAC))THEN
         jac = op_jac
      ELSE
         jac = jac_out
      ENDIF

      WRITE(out_unit,*) version
      WRITE(out_unit,'(1/,1x,a13,a8)') "jacobian = ", TRIM(jac)
      WRITE(out_unit,'(1x,1(a6,I6))')"n  =",nn
      WRITE(out_unit,'(1/,1x,a12,1x,I6,1x,2(a12,I4),1/)')
     $   "mstep =", mstep, "mpert =", mpert, "mthsurf =", mthsurf

      RETURN
      END SUBROUTINE ascii_header


c-----------------------------------------------------------------------
      SUBROUTINE output_profile_ascii( )
c-----------------------------------------------------------------------
c     *DESCRIPTION:
c        Write profile ascii files.
c     *ARGUMENTS:
c        None.
c
c-----------------------------------------------------------------------
      INTEGER :: istep, ipert, itheta, ipsi,  stepsize, funsize
      REAL(r8), DIMENSION(0:cmpsi) :: psi
      LOGICAL :: debug = .TRUE.

      ! start timer
      IF(verbose) WRITE(*,*)"Writing profile ascii output"
      IF(timeit) CALL gpec_timer(-2)

      ! common variables
      stepsize = MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
      funsize = MAX(1,(mstep*(mthsurf+1)-1)/max_linesout+1)
      psi = (/(ipsi,ipsi=0,cmpsi)/)/REAL(cmpsi,r8)

      ! write dw_profile
      IF ( dwk%mx>0 ) THEN
         IF(debug) WRITE(*,*) "Writing dw"
         CALL ascii_open(out_unit,"gpec_dw_n"//TRIM(sn)//".out",
     $      "UNKNOWN")
         WRITE(out_unit,*)"GPEC_DW_PROFILE: "//
     $        "Energy and torque profiles by self-consistent solutions."
         CALL ascii_header
         WRITE(out_unit,'(6(1x,a16))')"psi","T_phi","2ndeltaW",
     $        "int(T_phi)","int(2ndeltaW)","dv/dpsi_n"
         DO istep=1,mstep,stepsize
            WRITE(out_unit,'(6(1x,es16.8))') psifac(istep),
     $         dwk%fs1(istep,1),dwk%fs(istep,1),sq%fs(istep,3)
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      ! write pmodb
      IF ( ALLOCATED(eulbparmout) ) THEN
         IF(debug) WRITE(*,*) "  > Writing pmodb"
         CALL ascii_open(out_unit,"gpec_pmodb_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_PMODB: "//
     $        "Components in perturbed mod b and del.x_prp"
         CALL ascii_header
         WRITE(out_unit,'(1x,a16,1x,a4,10(1x,a16))')"psi","m",
     $        "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $        "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $        "imag(Bkxprp)"
         DO istep=1,mstep,stepsize
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
         IF (fun_flag) THEN
            IF(debug) WRITE(*,*) "  > Writing pmodb functions"
            CALL ascii_open(out_unit,"gpec_pmodb_fun_n"//
     $           TRIM(sn)//".out","UNKNOWN")
            WRITE(out_unit,*)"GPEC_PMODB_FUN: "//
     $           "Components in perturbed mod b in functions"
            CALL ascii_header
            WRITE(out_unit,'(15(1x,a16))')"psi","theta","r","z",
     $           "real(eulb)","imag(eulb)","real(lagb)","imag(lagb)",
     $           "real(Bdivxprp)","imag(Bdivxprp)","real(Bkxprp)",
     $           "imag(Bkxprp)","equilb","dequilbdpsi","dequilbdtheta"
            DO istep=1,mstep,funsize
               DO itheta=0,mthsurf
                  CALL bicube_eval(eqfun,psifac(istep),theta(itheta),1)
                  WRITE(out_unit,'(15(es17.8e3))')
     $                 psifac(istep),theta(itheta),
     $                 rs(istep,itheta),zs(istep,itheta),
     $                 REAL(eulbparfout(istep,itheta)),
     $                 -helicity*AIMAG(eulbparfout(istep,itheta)),
     $                 REAL(lagbparfout(istep,itheta)),
     $                 -helicity*AIMAG(lagbparfout(istep,itheta)),
     $                 REAL(-divxprpfout(istep,itheta)),
     $                 -helicity*AIMAG(-divxprpfout(istep,itheta)),
     $                 REAL(-curvfout(istep,itheta)),
     $                 -helicity*AIMAG(-curvfout(istep,itheta)),
     $                 equilbfun(istep,itheta),
     $                 eqfun%fx(1),eqfun%fy(1)
               ENDDO
            ENDDO
            CALL ascii_close(out_unit)
         ENDIF
      ENDIF

      ! write xbnormal
      IF ( ALLOCATED(bnomns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing xbnormal"
         CALL ascii_open(out_unit,"gpec_xbnormal_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBNORMAL: "//
     $        "Normal components of displacement and field"
         CALL ascii_header
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(xno)","imag(xno)","real(bno)","imag(bno)",
     $        "real(bwp)","imag(bwp)"
         DO istep=1,mstep,stepsize
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $            psifac(istep),qfac(istep),lmfac(ipert),
     $            REAL(xnomns(istep,ipert)),AIMAG(xnomns(istep,ipert)),
     $            REAL(bnomns(istep,ipert)),AIMAG(bnomns(istep,ipert)),
     $            REAL(bwpmns(istep,ipert)),AIMAG(bwpmns(istep,ipert))
            ENDDO
         ENDDO
      ENDIF
      IF (bwp_pest_flag .AND. ALLOCATED(pwpmns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing pest bnormal"
         CALL ascii_open(out_unit,"gpec_bnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_BNORMAL_PEST: "//
     $        "Normal components of field in pest coordinates"
         CALL ascii_header
         WRITE(out_unit,'(2(1x,a16),1x,a4,6(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         DO istep=1,mstep,stepsize
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,6(es17.8e3))')
     $              psifac(istep),qfac(istep),lmfac(ipert),
     $              REAL(pwpmns(istep,ipert)),AIMAG(pwpmns(istep,ipert))
            ENDDO
         ENDDO
      ENDIF
      IF (fun_flag .AND. ALLOCATED(bnofuns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing xbnormal functions"
         CALL ascii_open(out_unit,"gpec_xbnormal_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBNORMAL_FUN: "//
     $        "Normal components of displacement and field in functions"
         CALL ascii_header
         WRITE(out_unit,'(10(1x,a16))')"psi","theta","r","z","rvec",
     $        "zvec","real(xno)","imag(xno)","real(bno)","imag(bno)"
         DO istep=1,mstep,funsize
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

      ! write xbtangent
      IF ( ALLOCATED(xtamns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing xbtangent"
         CALL ascii_open(out_unit,"gpec_xbtangent_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBNROMAL: "//
     $        "Tangential components of displacement and field"
         CALL ascii_header
         WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $        "real(xta)","imag(xta)","real(bta)","imag(bta)"
         DO istep=1,mstep,stepsize
            DO ipert=1,mpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,4(es17.8e3))')
     $            psifac(istep),qfac(istep),mfac(ipert),
     $            REAL(xtamns(istep,ipert)),AIMAG(xtamns(istep,ipert)),
     $            REAL(btamns(istep,ipert)),AIMAG(btamns(istep,ipert))
            ENDDO
         ENDDO
      ENDIF
      IF ( fun_flag .AND. ALLOCATED(xtafuns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing xbtangent functions"
         CALL ascii_open(out_unit,"gpec_xbtangent_fun_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XBTANGENT_FUN: "//
     $        "Tangential components of displacement "//
     $        "and field in functions"
         CALL ascii_header
         WRITE(out_unit,'(8(1x,a16))')"r","z","rvec","zvec",
     $        "real(xta)","imag(xta)","real(bta)","imag(bta)"
         DO istep=1,mstep,funsize
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

      ! write vbnormal
      IF ( ALLOCATED(vnomns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing vbnormal"
         CALL ascii_open(out_unit,"gpec_vbnormal_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBNROMAL: "//
     $        "Normal components of coil field"
         CALL ascii_header
         WRITE(out_unit,'(2(1x,a16),1x,a4,4(1x,a16))')"psi","q","m",
     $        "real(bno)","imag(bno)","real(bwp)","imag(bwp)"
         DO ipsi=1,cmpsi
            CALL spline_eval(sq,psi(ipsi),0)
            DO ipert=1,mpert
               WRITE(out_unit,'(2(1x,es16.8),1x,I4,4(1x,es16.8))')
     $              psi(ipsi),sq%f(4),mfac(ipert),
     $              REAL(vnomns(ipsi,ipert)),AIMAG(vnomns(ipsi,ipert)),
     $              REAL(vwpmns(ipsi,ipert)),AIMAG(vwpmns(ipsi,ipert))
            ENDDO
         ENDDO
      ENDIF
      IF ( ALLOCATED(vpwpmns) ) THEN
         IF(debug) WRITE(*,*) "  > Writing vbnormal pest"
         CALL ascii_open(out_unit,"gpec_vbnormal_pest_n"//
     $        TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_VBNROMAL_PEST: "//
     $        "Normal components of coil field in pest coordinates"
         CALL ascii_header
         WRITE(out_unit,'(2(1x,a16),1x,a4,2(1x,a16))')"psi","q","m",
     $        "real(bwp)","imag(bwp)"
         DO ipsi=1,cmpsi
            CALL spline_eval(sq,psi(ipsi),0)
            DO ipert=1,lmpert
               WRITE(out_unit,'(2(es17.8e3),1x,I4,2(es17.8e3))')
     $              psi(ipsi),sq%f(4),lmfac(ipert),
     $              REAL(vpwpmns(ipsi,ipert)),AIMAG(vpwpmns(ipsi,ipert))
            ENDDO
         ENDDO
      ENDIF

      ! write xclebsch
      IF ( ALLOCATED(xmp1out) ) THEN
         IF(debug) WRITE(*,*) "  > Writing xclebsch"
         CALL ascii_open(out_unit,"gpec_xclebsch_n"//
     $      TRIM(sn)//".out","UNKNOWN")
         WRITE(out_unit,*)"GPEC_XCLEBSCH: "//
     $      "Clebsch Components of the displacement."
         CALL ascii_header
         WRITE(out_unit,'(1x,a23,1x,a4,6(1x,a16))')"psi","m",
     $      "real(derxi^psi)","imag(derxi^psi)",
     $      "real(xi^psi)","imag(xi^psi)",
     $      "real(xi^alpha)","imag(xi^alpha)"
         DO istep=1,mstep,stepsize
            DO ipert=1,lmpert
               WRITE(out_unit,'(1x,es23.15,1x,I4,6(es17.8e3))')
     $            psifac(istep),mfac(ipert),xmp1out(istep,ipert),
     $            xspout(istep,ipert),xmsout(istep,ipert)/chi1
            ENDDO
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF

      ! end timer
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE output_profile_ascii


c-----------------------------------------------------------------------
      SUBROUTINE output_profile_netcdf( )
c-----------------------------------------------------------------------
c     *DESCRIPTION:
c        Write netcdf file with all calculated profiles.
c     *ARGUMENTS:
c        None.
c
c-----------------------------------------------------------------------
      INTEGER :: ncid
      CHARACTER(64) :: ncfile
      INTEGER :: idid, mdid, pdid, tdid, i_id, m_id, p_id, t_id

      INTEGER :: x1_id, xp_id, xa_id
      INTEGER :: r_id, z_id, rv_id, zv_id
      INTEGER :: be_id, bme_id, bl_id, bml_id, eb_id,
     $   dx_id, dxm_id, k_id, km_id
      INTEGER :: b_id, bm_id, x_id, xm_id

      ! the usual messages
      IF(verbose) WRITE(*,*)"Writing profile netcdf output"
      IF(timeit) CALL gpec_timer(-2)

      ! initialize file
      IF(debug_flag) PRINT *,"Initializing profile NETCDF file"
      ncfile = "gpec_profile_output_n"//TRIM(sn)//".nc"
      CALL check( nf90_create(ncfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )
     
      ! define global file attributes
      IF(debug_flag) PRINT *," - Defining profile netcdf globals"
      CALL check( nf90_put_att(ncid,nf90_global,"title",
     $     "GPEC profile outputs in magnetic coordinate system"))
      CALL check( nf90_put_att(ncid,nf90_global,"Jacobian",jac_type))
      CALL check( nf90_put_att(ncid,nf90_global,"q_lim",qlim))
      CALL check( nf90_put_att(ncid,nf90_global,"psi_n_lim",psilim))
      CALL check( nf90_put_att(ncid,nf90_global,"shot",INT(shotnum)) )
      CALL check( nf90_put_att(ncid,nf90_global,"time",INT(shottime)) )
      CALL check( nf90_put_att(ncid,nf90_global,"machine",machine) )
      CALL check( nf90_put_att(ncid,nf90_global,"n",nn) )
      CALL check( nf90_put_att(ncid,nf90_global,"version",version))
      
      ! define dimensions
      CALL check( nf90_def_dim(ncid,"i",2,       idid) )
      CALL check( nf90_def_var(ncid,"i",nf90_int,idid,i_id) )
      CALL check( nf90_def_dim(ncid,"m",lmpert,  mdid) )
      CALL check( nf90_def_var(ncid,"m",NF90_INT,mdid,m_id) )
      CALL check( nf90_def_dim(ncid,"psi_n",mstep,      pdid) )
      CALL check( nf90_def_var(ncid,"psi_n",nf90_double,pdid,p_id) )
      CALL check( nf90_def_dim(ncid,"theta",mthsurf+1,  tdid) )
      CALL check( nf90_def_var(ncid,"theta",nf90_double,tdid,t_id) )

      ! define variables
      IF(debug_flag) PRINT *,"  Defining variables"
      IF ( ALLOCATED(xmp1out) ) THEN
         CALL check( nf90_def_var(ncid, "derxi_m_contrapsi",
     $        nf90_double, (/p_id,m_id,i_id/), x1_id) )
         CALL check( nf90_put_att(ncid, x1_id, "units", "m") )
         CALL check( nf90_put_att(ncid, x1_id, "long_name",
     $        "Psi derivative of contravarient psi displacement") )
         CALL check( nf90_def_var(ncid, "xi_m_contrapsi", nf90_double,
     $        (/p_id,m_id,i_id/), xp_id) )
         CALL check( nf90_put_att(ncid, xp_id, "units", "m") )
         CALL check( nf90_put_att(ncid, xp_id, "long_name",
     $        "Contravarient psi displacement") )
         CALL check( nf90_def_var(ncid, "xi_m_contraalpha", nf90_double,
     $        (/p_id,m_id,i_id/), xa_id) )
         CALL check( nf90_put_att(ncid, xa_id, "units", "m") )
         CALL check( nf90_put_att(ncid, xa_id, "long_name",
     $        "Contravarient clebsch angle displacement") )
      ENDIF

      IF( ALLOCATED(rs) )THEN
         CALL check( nf90_def_var(ncid, "R", nf90_double,
     $               (/p_id,t_id/),r_id) )
         CALL check( nf90_put_att(ncid,r_id,"long_name",
     $               "Major Radius") )
         CALL check( nf90_put_att(ncid,r_id,"units","m") )
         CALL check( nf90_def_var(ncid, "z", nf90_double,
     $               (/p_id,t_id/),z_id) )
         CALL check( nf90_put_att(ncid,z_id,"long_name",
     $               "Vertical Position") )
         CALL check( nf90_put_att(ncid,z_id,"units","m") )
      ENDIF

      IF ( ALLOCATED(rvecs) ) THEN
         CALL check( nf90_def_var(ncid, "R_n", nf90_double,
     $               (/p_id,t_id/),rv_id) )
         CALL check( nf90_put_att(ncid,rv_id,"long_name",
     $               "Radial unit normal: R_3D = R_sym+x_n*R_n") )
         CALL check( nf90_def_var(ncid, "z_n", nf90_double,
     $               (/p_id,t_id/),zv_id) )
         CALL check( nf90_put_att(ncid,zv_id,"long_name",
     $               "Vertical unit normal: z_3D = z_sym+x_n*z_n") )
      ENDIF

      IF( ALLOCATED(eulbparmout) ) THEN
         CALL check( nf90_def_var(ncid, "b_eul", nf90_double,
     $               (/p_id,t_id,i_id/),be_id) )
         CALL check( nf90_put_att(ncid,be_id,"long_name",
     $               "Eulerian Perturbed Field") )
         CALL check( nf90_put_att(ncid,be_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "b_lag", nf90_double,
     $               (/p_id,t_id,i_id/),bl_id) )
         CALL check( nf90_put_att(ncid,bl_id,"long_name",
     $               "Lagrangian Perturbed Field") )
         CALL check( nf90_put_att(ncid,bl_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "b_m-eul", nf90_double,
     $               (/p_id,m_id,i_id/),bme_id) )
         CALL check( nf90_put_att(ncid,bme_id,"long_name",
     $               "Eulerian Perturbed Field") )
         CALL check( nf90_put_att(ncid,bme_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "b_m-lag", nf90_double,
     $               (/p_id,m_id,i_id/),bml_id) )
         CALL check( nf90_put_att(ncid,bml_id,"long_name",
     $               "Lagrangian Perturbed Field") )
         CALL check( nf90_put_att(ncid,bml_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "Bdivxi_perp", nf90_double,
     $               (/p_id,t_id,i_id/),dx_id) )
         CALL check( nf90_put_att(ncid,dx_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_def_var(ncid, "Bdivxi_m-perp", nf90_double,
     $               (/p_id,m_id,i_id/),dxm_id) )
         CALL check( nf90_put_att(ncid,dxm_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_def_var(ncid, "Bkappaxi_perp", nf90_double,
     $               (/p_id,t_id,i_id/),k_id) )
         CALL check( nf90_put_att(ncid,k_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_def_var(ncid, "Bkappaxi_m-perp", nf90_double,
     $               (/p_id,m_id,i_id/),km_id) )
         CALL check( nf90_put_att(ncid,km_id,"long_name",
     $               "Divergence of the normal displacement") )
         CALL check( nf90_def_var(ncid, "B", nf90_double,
     $               (/p_id,t_id/),eb_id) )
         CALL check( nf90_put_att(ncid,eb_id,"long_name",
     $               "Equilibrium Field Strength") )
      ENDIF

      IF ( ALLOCATED(bnomns) ) THEN
         CALL check( nf90_def_var(ncid, "b_n", nf90_double,
     $               (/p_id,t_id,i_id/),b_id) )
         CALL check( nf90_put_att(ncid,b_id,"long_name",
     $               "Perturbed Field Normal to the Flux Surface") )
         CALL check( nf90_put_att(ncid,b_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "b_nm", nf90_double,
     $               (/p_id,m_id,i_id/),bm_id) )
         CALL check( nf90_put_att(ncid,bm_id,"long_name",
     $               "Perturbed Field Normal to the Flux Surface") )
         CALL check( nf90_put_att(ncid,bm_id,"units","Tesla") )
         CALL check( nf90_def_var(ncid, "xi_n", nf90_double,
     $               (/p_id,t_id,i_id/),x_id) )
         CALL check( nf90_put_att(ncid,x_id,"long_name",
     $               "Displacement Normal to the Flux Surface") )
         CALL check( nf90_put_att(ncid,x_id,"units","m") )
         CALL check( nf90_def_var(ncid, "xi_nm", nf90_double,
     $               (/p_id,m_id,i_id/),xm_id) )
         CALL check( nf90_put_att(ncid,xm_id,"long_name",
     $               "Displacement Normal to the Flux Surface") )
         CALL check( nf90_put_att(ncid,xm_id,"units","m") )
      ENDIF

      ! end of definitions
      CALL check( nf90_enddef(ncid) )
      
      ! put coordinate dimensions in file
      IF(debug_flag) PRINT *," - Putting coordinates in flux netcdfs"
      CALL check( nf90_put_var(ncid,i_id,(/0,1/)) )
      CALL check( nf90_put_var(ncid,m_id,lmfac) )
      CALL check( nf90_put_var(ncid,p_id,psifac(1:)) )
      CALL check( nf90_put_var(ncid,t_id,theta) )

      ! put variables in file
      IF ( ALLOCATED(xmp1out) ) THEN
         CALL check( nf90_put_var(ncid, x1_id, RESHAPE((/REAL(xmp1out),
     $        AIMAG(xmp1out)/), (/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(ncid, xp_id, RESHAPE((/REAL(xspout),
     $        AIMAG(xspout)/), (/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(ncid, xa_id, RESHAPE((/REAL(xmsout),
     $        AIMAG(xmsout)/), (/mstep,lmpert,2/))) )
      ENDIF
      IF ( ALLOCATED(rs) ) THEN
         CALL check( nf90_put_var(ncid,r_id,rs) )
         CALL check( nf90_put_var(ncid,z_id,zs) )
      ENDIF
      IF ( ALLOCATED(eulbparmout) ) THEN
         CALL check( nf90_put_var(ncid,bme_id,RESHAPE((/
     $      REAL(eulbparmout),AIMAG(eulbparmout)/),(/mstep,lmpert,2/))))
         CALL check( nf90_put_var(ncid,bml_id,RESHAPE((/
     $      REAL(lagbparmout),AIMAG(lagbparmout)/),(/mstep,lmpert,2/))))
         CALL check( nf90_put_var(ncid,dxm_id,RESHAPE((/REAL(
     $      -divxprpmout),AIMAG(-divxprpmout)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(ncid,km_id,RESHAPE((/
     $      REAL(-curvmout),AIMAG(-curvmout)/),(/mstep,lmpert,2/))) )
         CALL check(nf90_put_var(ncid,be_id,RESHAPE((/REAL(eulbparfout),
     $      -helicity*AIMAG(eulbparfout)/),(/mstep,mthsurf,2/))) )
         CALL check(nf90_put_var(ncid,bl_id,RESHAPE((/REAL(lagbparfout),
     $      -helicity*AIMAG(lagbparfout)/),(/mstep,mthsurf,2/))) )
         CALL check(nf90_put_var(ncid,dx_id,RESHAPE((/REAL(-divxprpfout)
     $      ,-helicity*AIMAG(-divxprpfout)/),(/mstep,mthsurf,2/))))
         CALL check( nf90_put_var(ncid,k_id,RESHAPE((/REAL(-curvfout),
     $      -helicity*AIMAG(-curvfout)/),(/mstep,mthsurf,2/))) )
         CALL check( nf90_put_var(ncid,eb_id,equilbfun) )
      ENDIF
      IF ( ALLOCATED(rvecs) ) THEN
         CALL check( nf90_put_var(ncid,rv_id,rvecs) )
         CALL check( nf90_put_var(ncid,zv_id,zvecs) )
      ENDIF
      IF ( ALLOCATED(bnomns) ) THEN
         CALL check( nf90_put_var(ncid,xm_id,RESHAPE((/REAL(xnomns),
     $        AIMAG(xnomns)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(ncid,bm_id,RESHAPE((/REAL(bnomns),
     $        AIMAG(bnomns)/),(/mstep,lmpert,2/))) )
         CALL check( nf90_put_var(ncid,x_id,RESHAPE((/REAL(xnofuns),
     $        -helicity*AIMAG(xnofuns)/),(/mstep,mthsurf+1,2/))) )
         CALL check( nf90_put_var(ncid,b_id,RESHAPE((/REAL(bnofuns),
     $        -helicity*AIMAG(bnofuns)/),(/mstep,mthsurf+1,2/))) )
      ENDIF

      ! close file
      IF(debug_flag) PRINT *," - Closing "//TRIM(ncfile)
      CALL check( nf90_close(ncid) )

      ! end timer
      IF(timeit) CALL gpec_timer(2)
      RETURN
      END SUBROUTINE output_profile_netcdf

c-----------------------------------------------------------------------
      SUBROUTINE check(stat)
c-----------------------------------------------------------------------
c     *DESCRIPTION:
c        Check status of netcdf file and stop program if there is an 
c        error.
c     *ARGUMENTS:
c        stat: integer. Status returned by a netcdf module function.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
      
      ! stop if it is an error.
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
      
      RETURN
      END SUBROUTINE check


      END MODULE output_profile
