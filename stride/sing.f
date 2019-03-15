c-----------------------------------------------------------------------
c     file sing.f.
c     computations relating to singular surfaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. sing_mod.
c     1. sing_parallel_alloc.
c     2. sing_asymp.
c     3. sing_find.
c     4. sing_lim.
c     5. sing_vmat.
c     6. sing_mmat.
c     7. sing_solve.
c     8. sing_matmul.
c     9. sing_get_ua.
c     10. sing_der.
c     11. sing_derFM.
c     12. sing_derL.
c-----------------------------------------------------------------------
c     subprogram 0. sing_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sing_mod
      USE fourfit_mod
      IMPLICIT NONE

      INTEGER :: sing_order=2

      TYPE sing_calculator
         REAL(r8), DIMENSION(2) :: singEdgesLR
         INTEGER :: singIntervalL,singIntervalR,singInterval0
         INTEGER :: singMidPt, sing_col
         COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: uShootL, uShootR
      END TYPE sing_calculator
c-----------------------------------------------------------------------
c     declarations for parallelization (of sing.f variables).
c-----------------------------------------------------------------------
      INTEGER :: locstab_s_ix
      REAL(r8), DIMENSION(:), ALLOCATABLE :: locstab_s_f
      REAL(r8), DIMENSION(:), ALLOCATABLE :: sq_s_f1,sq_s_f2,sq_s_f3
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: fmats_s_f1,fmats_s_f2,
     $     fmats_s_f3
c-----------------------------------------------------------------------
c     declarations for parallelization (of global variables)
c-----------------------------------------------------------------------
      INTEGER :: sq_s_ix,fmats_s_ix,gmats_s_ix,kmats_s_ix
      REAL(r8), DIMENSION(:), ALLOCATABLE :: sq_s_f
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE ::
     $     fmats_s_f,gmats_s_f,kmats_s_f
!$OMP THREADPRIVATE(
!$OMP& sq_s_ix,fmats_s_ix,gmats_s_ix,kmats_s_ix,
!$OMP& sq_s_f,fmats_s_f,gmats_s_f,kmats_s_f,
!$OMP& locstab_s_ix,locstab_s_f,sq_s_f1,sq_s_f2,sq_s_f3,
!$OMP& fmats_s_f1,fmats_s_f2,fmats_s_f3)
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. sing_parallel_alloc.
c     allocate global threadprivate variables in parallel region
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_parallel_alloc(allocThreads)
      INTEGER, INTENT(IN) :: allocThreads
c-----------------------------------------------------------------------
c     allocate.
c-----------------------------------------------------------------------
      IF (verbose_performance_output) THEN
         print *,"sing_parallel_alloc: #Threads=",allocThreads
      ENDIF
      CALL OMP_SET_NUM_THREADS(allocThreads)
!$OMP PARALLEL DEFAULT(NONE) SHARED(locstab,fmats,sq)
!$OMP& COPYIN(sq_s_f,fmats_s_f,
!$OMP& locstab_s_f,sq_s_f1,sq_s_f2,sq_s_f3,
!$OMP& fmats_s_f1,fmats_s_f2,fmats_s_f3)
      ALLOCATE(locstab_s_f(SIZE(locstab%f)))
      ALLOCATE(sq_s_f(SIZE(sq%f)),sq_s_f1(SIZE(sq%f1)),
     $     sq_s_f2(SIZE(sq%f2)),sq_s_f3(SIZE(sq%f3)))
      ALLOCATE(fmats_s_f(SIZE(fmats%f)),fmats_s_f1(SIZE(fmats%f1)),
     $     fmats_s_f2(SIZE(fmats%f2)),fmats_s_f3(SIZE(fmats%f3)))
!$OMP END PARALLEL
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_parallel_alloc
c-----------------------------------------------------------------------
c     subprogram 2. sing_asymp.
c     generates asymptotic data at each singular surface
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_asymp
      INTEGER :: ising, sTime, fTime, cr
      REAL(r8) :: di0

 10   FORMAT(/2x,"i",5x,"psi",8x,"rho",9x,"q",10x,"q1",8x,"di0",9x,"di",
     $     8x,"err"/)

      CALL SYSTEM_CLOCK(COUNT=sTime)
      CALL SYSTEM_CLOCK(COUNT_RATE=cr)
c-----------------------------------------------------------------------
c     set up parallel loop.  NOTE: ALLOCATEs occur w/in parallel region!
c-----------------------------------------------------------------------
      locstab_s_ix = locstab%ix
      sq_s_ix = sq%ix
      fmats_s_ix = fmats%ix

      !Note: nThreads faster than nThreads-1, despite creating with
      !      vac_parallel more threads than processors.
      IF (verbose_performance_output) THEN
         print *,"sing_asymp: #Threads=",nThreads
      ENDIF
      CALL OMP_SET_NUM_THREADS(nThreads)
!$OMP PARALLEL DEFAULT(NONE)
!$OMP& SHARED(msing,sing,sing_order,mlow,mhigh,locstab,
!$OMP& nn,fmats,gmats,kmats,sTime,cr,mpert,mband,
!$OMP& verbose_performance_output)
!$OMP& PRIVATE(ising,fTime)
!$OMP& COPYIN(sq,sq_s_f,fmats_s_f,locstab_s_f,sq_s_f1,sq_s_f2,sq_s_f3,
!$OMP& fmats_s_f1,fmats_s_f2,fmats_s_f3,locstab_s_ix,sq_s_ix,fmats_s_ix)
!$OMP DO SCHEDULE(GUIDED)
c-----------------------------------------------------------------------
c     loop over singular surfaces.
c-----------------------------------------------------------------------
      DO ising = 1,msing
         CALL sing_vmat(ising)
         CALL SYSTEM_CLOCK(COUNT=fTime)
         IF (verbose_performance_output) THEN
            print *,"*** asymp time #",ising,"=",
     $           REAL(fTime-sTime,8)/REAL(cr,8)
         ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
c-----------------------------------------------------------------------
c     write to table of singular surface quantities.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/1x,a)')"Singular Surfaces:"
      WRITE(out_unit,10)
      DO ising=1,msing
         CALL spline_eval(locstab,sing(ising)%psifac,0)
         di0=locstab%f(1)/sing(ising)%psifac
         WRITE(out_unit,'(i3,1p,7e11.3)')ising,sing(ising)%psifac,
     $         sing(ising)%rho,sing(ising)%q,sing(ising)%q1,
     $         di0,sing(ising)%di,sing(ising)%di/di0-1
      ENDDO
      WRITE(out_unit,10)
c-----------------------------------------------------------------------
c     timer.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print *,"*** asymp time=",REAL(fTime-sTime,8)/REAL(cr,8)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_asymp
c-----------------------------------------------------------------------
c     subprogram 3. sing_find.
c     finds positions of singular values of q.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_find

      INTEGER, PARAMETER :: itmax=200,nsing=1000
      INTEGER :: iex,ising,m,dm,it
      INTEGER, DIMENSION(nsing) :: m_sing
      REAL(r8) :: dq,psifac,psifac0,psifac1,singfac
      REAL(r8), DIMENSION(nsing) :: psising,qsing,q1sing
c-----------------------------------------------------------------------
c     start loop over extrema to find singular surfaces.
c-----------------------------------------------------------------------
      ising=0
      DO iex=1,mex
         dq=qex(iex)-qex(iex-1)
         m=nn*qex(iex-1)
         IF(dq > 0)m=m+1
         dm=SIGN(one,dq*nn)
c-----------------------------------------------------------------------
c     find singular surfaces by binary search.
c-----------------------------------------------------------------------
         DO
            IF((m-nn*qex(iex-1))*(m-nn*qex(iex)) > 0)EXIT
            it=0
            psifac0=psiex(iex-1)
            psifac1=psiex(iex)
            DO
               it=it+1
               psifac=(psifac0+psifac1)/2
               CALL spline_eval(sq,psifac,0)
               singfac=(m-nn*sq%f(4))*dm
               IF(singfac > 0)THEN
                  psifac0=psifac
                  psifac=(psifac+psifac1)/2
               ELSE
                  psifac1=psifac
                  psifac=(psifac+psifac0)/2
               ENDIF
               IF(ABS(singfac) <= 1e-12)EXIT
               IF(it > itmax)
     $              CALL program_stop("sing_find can't find root")
            ENDDO
c-----------------------------------------------------------------------
c     store singular surfaces.
c-----------------------------------------------------------------------
            ising=ising+1
            CALL spline_eval(sq,psifac,1)
            m_sing(ising)=m
            qsing(ising)=REAL(m,r8)/nn
            q1sing(ising)=sq%f1(4)
            psising(ising)=psifac
            m=m+dm
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     transfer to permanent storage.
c-----------------------------------------------------------------------
      msing=ising
      ALLOCATE(sing(msing))
      DO ising=1,msing
         sing(ising)%m=m_sing(ising)
         sing(ising)%psifac=psising(ising)
         sing(ising)%rho=SQRT(psising(ising))
         sing(ising)%q=qsing(ising)
         sing(ising)%q1=q1sing(ising)
c-----------------------------------------------------------------------
c     print singular surfaces details.
c-----------------------------------------------------------------------
         print '("Sing# = ",i12)',ising
         print '("--------------------")'
         print '("m = ",i12)',sing(ising)%m
         print '("psifac = ",f15.9)',sing(ising)%psifac
         print '("q = ",f15.9)',sing(ising)%q
         print '("--------------------")'
      ENDDO
      DEALLOCATE(psiex,qex)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_find
c-----------------------------------------------------------------------
c     subprogram 4. sing_lim.
c     computes limiter values.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_lim

      INTEGER :: it,itmax=50
      INTEGER, DIMENSION(1) :: jpsi
      REAL(r8) :: dpsi,q,q1,eps=1e-10
c-----------------------------------------------------------------------
c     compute psilim and qlim.
c-----------------------------------------------------------------------
      qlim=MIN(qmax,qhigh)
      q1lim=sq%fs1(mpsi,4)
      psilim=psihigh
      IF(sas_flag)THEN
c-----------------------------------------------------------------------
c     normalze dmlim to interval [0,1).
c-----------------------------------------------------------------------
         DO
            IF(dmlim < 1)EXIT
            dmlim=dmlim-1
         ENDDO
         DO
            IF(dmlim >= 0)EXIT
            dmlim=dmlim+1
         ENDDO
c-----------------------------------------------------------------------
c     compute qlim.
c-----------------------------------------------------------------------
         qlim=(INT(nn*qlim)+dmlim)/nn
         DO
            IF(qlim <= qmax)EXIT
            qlim=qlim-1._r8/nn
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     use newton iteration to find psilim.
c-----------------------------------------------------------------------
      IF(qlim/=qmax)THEN
         jpsi=MINLOC(ABS(sq%fs(:,4)-qlim))
         IF (jpsi(1)>= mpsi) jpsi(1)=mpsi-1
         psilim=sq%xs(jpsi(1))
         it=0
         DO
            it=it+1
            CALL spline_eval(sq,psilim,1)
            q=sq%f(4)
            q1=sq%f1(4)
            dpsi=(qlim-q)/q1
            psilim=psilim+dpsi
            IF(ABS(dpsi) < eps*ABS(psilim) .OR. it > itmax)EXIT
         ENDDO
         q1lim=q1
c-----------------------------------------------------------------------
c     abort if not found.
c-----------------------------------------------------------------------
         IF(it > itmax)THEN
            CALL program_stop("Can't find psilim.")
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_lim
c-----------------------------------------------------------------------
c     subprogram 5. sing_vmat.
c     computes asymptotic behavior at the singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_vmat(ising)

      INTEGER, INTENT(IN) :: ising

      CHARACTER(5), DIMENSION(2) :: side=(/"left ","right"/)
      INTEGER :: ipert0,ipert,k,iside
      REAL(r8) :: psifac,di,di0,q,q1,rho,sig
      REAL(r8), PARAMETER :: half=.5_r8
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat,mmat
      COMPLEX(r8), DIMENSION(2,2) :: m0mat
      INTEGER, DIMENSION(:), POINTER :: r1,r2,n1,n2
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointer to singular surface.
c-----------------------------------------------------------------------
      IF(ising < 1 .OR. ising > msing)RETURN
      singp => sing(ising)
c-----------------------------------------------------------------------
c     zero di if out of range.
c-----------------------------------------------------------------------
      ipert0=NINT(nn*singp%q)-mlow+1
      q=singp%q
      IF(ipert0 <= 0 .OR. mlow > nn*q .OR. mhigh < nn*q)THEN
         singp%di=0
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     allocate and compute ranges.
c-----------------------------------------------------------------------
      ALLOCATE(sing(ising)%n1(mpert-1),sing(ising)%n2(2*mpert-2))
      singp%r1=(/ipert0/)
      singp%r2=(/ipert0,ipert0+mpert/)
      singp%n1=(/(ipert,ipert=1,ipert0-1),(ipert,ipert=ipert0+1,mpert)/)
      singp%n2=(/singp%n1,singp%n1+mpert/)
      r1 => singp%r1
      r2 => singp%r2
      n1 => singp%n1
      n2 => singp%n2
c-----------------------------------------------------------------------
c     interpolate local values.
c-----------------------------------------------------------------------
      psifac=singp%psifac
      CALL spline_eval_external(locstab,psifac,locstab_s_ix,
     $     locstab_s_f)
      di0=locstab_s_f(1)/psifac
      q1=singp%q1
      rho=singp%rho
c-----------------------------------------------------------------------
c     compute Mercier criterion.
c-----------------------------------------------------------------------
      singp%order=0
      CALL sing_mmat(ising,r1,r2,n1,n2)
      m0mat=TRANSPOSE(singp%mmatr(r1(1),r2,:,0))
      di=m0mat(1,1)*m0mat(2,2)-m0mat(2,1)*m0mat(1,2)
      WRITE (*,*) "di=",di,"ising=",ising
      singp%di=di
      singp%alpha=SQRT(-CMPLX(singp%di))
c-----------------------------------------------------------------------
c     reset order and compute higher-order terms.
c-----------------------------------------------------------------------
      singp%order=sing_order
      DEALLOCATE(singp%mmatl,singp%mmatr)
      CALL sing_mmat(ising,r1,r2,n1,n2)
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      ALLOCATE(singp%power(2*mpert))
      singp%power=0
      singp%power(ipert0)=-singp%alpha
      singp%power(ipert0+mpert)=singp%alpha
c-----------------------------------------------------------------------
c     allocate vmatl and vmatr at singular surface.
c-----------------------------------------------------------------------
      ALLOCATE(singp%vmatl(mpert,2*mpert,2,0:2*singp%order))
      ALLOCATE(singp%vmatr(mpert,2*mpert,2,0:2*singp%order))
c-----------------------------------------------------------------------
c     set pointers and compute m0mat.
c-----------------------------------------------------------------------
      DO iside=1,2
         SELECT CASE(side(iside))
         CASE ("left")
            vmat => singp%vmatl
            mmat => singp%mmatl
            m0mat=TRANSPOSE(singp%mmatl(r1(1),r2,:,0))
            sig=-1
         CASE ("right")
            vmat => singp%vmatr
            mmat => singp%mmatr
            m0mat=TRANSPOSE(singp%mmatr(r1(1),r2,:,0))
            sig=1
         END SELECT
c-----------------------------------------------------------------------
c     zeroth-order non-resonant solutions.
c-----------------------------------------------------------------------
         vmat=0
         DO ipert=1,mpert
            vmat(ipert,ipert,1,0)=1
            vmat(ipert,ipert+mpert,2,0)=1
         ENDDO
c-----------------------------------------------------------------------
c     zeroth-order resonant solutions.
c-----------------------------------------------------------------------
         vmat(ipert0,ipert0,1,0)=1
         vmat(ipert0,ipert0+mpert,1,0)=1
         vmat(ipert0,ipert0,2,0)
     $        =-(m0mat(1,1)+sig*singp%alpha)/m0mat(1,2)
         vmat(ipert0,ipert0+mpert,2,0)
     $        =-(m0mat(1,1)-sig*singp%alpha)/m0mat(1,2)
c-----------------------------------------------------------------------
c     compute higher-order solutions.
c-----------------------------------------------------------------------
         DO k=1,2*singp%order
            CALL sing_solve(k,sig,singp%power,mmat,vmat,m0mat,r1,n1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_vmat
c-----------------------------------------------------------------------
c     subprogram 6. sing_mmat.
c     computes series expansion of coefficient matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_mmat(ising,r1,r2,n1,n2)

      INTEGER, INTENT(IN) :: ising
      INTEGER, DIMENSION(:), POINTER, INTENT(IN) :: r1,r2,n1,n2

      CHARACTER(5), DIMENSION(2) :: side=(/"left ","right"/)
      INTEGER :: ipert0,ipert,jpert,kpert,isol,iqty,info,msol,m,i,j,n,
     $     fac0,fac1,iside
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: psifac,sig
      REAL(r8), DIMENSION(0:3) :: q
      REAL(r8), DIMENSION(mpert,0:3) :: singfac
      COMPLEX(r8), PARAMETER :: one=1,half=.5_r8
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: f0
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: v
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: mmat

      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: f,ff,g,k
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: x

      TYPE(sing_type), POINTER :: singp
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: gpot_s_f,gpot_s_f1,
     $     gpot_s_f2,gpot_s_f3,kpot_s_f,kpot_s_f1,kpot_s_f2,kpot_s_f3
      INTEGER :: gpot_s_ix,kpot_s_ix
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"io",4x,"is",4x,"ip",2(4x,"re x",i1,6x,"im x",i1,2x)/)
 20   FORMAT(3i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     evaluate cubic splines.
c-----------------------------------------------------------------------
      ALLOCATE(gpot_s_f(SIZE(gmats%f)),gpot_s_f1(SIZE(gmats%f1)),
     $     gpot_s_f2(SIZE(gmats%f2)),gpot_s_f3(SIZE(gmats%f3)),
     $     kpot_s_f(SIZE(kmats%f)),kpot_s_f1(SIZE(kmats%f1)),
     $     kpot_s_f2(SIZE(kmats%f2)),kpot_s_f3(SIZE(kmats%f3)))
      gpot_s_ix = gmats%ix
      kpot_s_ix = kmats%ix

      msol=2*mpert
      singp => sing(ising)
      psifac=singp%psifac

      CALL spline_eval_external(sq,psifac,sq_s_ix,
     $     sq_s_f,sq_s_f1,sq_s_f2,sq_s_f3)
      CALL cspline_eval_external(fmats,psifac,fmats_s_ix,
     $     fmats_s_f,fmats_s_f1,fmats_s_f2,fmats_s_f3)
      CALL cspline_eval_external(gmats,psifac,gpot_s_ix,
     $     gpot_s_f,gpot_s_f1,gpot_s_f2,gpot_s_f3)
      CALL cspline_eval_external(kmats,psifac,kpot_s_ix,
     $     kpot_s_f,kpot_s_f1,kpot_s_f2,kpot_s_f3)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(f(mband+1,mpert,0:singp%order))
      ALLOCATE(ff(mband+1,mpert,0:singp%order))
      ALLOCATE(g(mband+1,mpert,0:singp%order))
      ALLOCATE(k(2*mband+1,mpert,0:singp%order))
      ALLOCATE(x(mpert,2*mpert,2,0:singp%order))
      ALLOCATE(singp%mmatl(mpert,2*mpert,2,0:2*singp%order+2))
      ALLOCATE(singp%mmatr(mpert,2*mpert,2,0:2*singp%order+2))
c-----------------------------------------------------------------------
c     computation of M matrix on both side of resonant surface.
c-----------------------------------------------------------------------
      DO iside=1,2
         SELECT CASE(side(iside))
         CASE("left")
            sig=-1.0
            mmat => singp%mmatl
         CASE("right")
            sig=1.0
            mmat => singp%mmatr
         END SELECT
c-----------------------------------------------------------------------
c     evaluate safety factor and its derivatives.
c-----------------------------------------------------------------------
         q(0)=sq_s_f(4)
         q(1)=sig*sq_s_f1(4)
         q(2)=sq_s_f2(4)
         q(3)=sig*sq_s_f3(4)
c-----------------------------------------------------------------------
c     evaluate singfac and its derivatives.
c-----------------------------------------------------------------------
         ipert0=singp%m-mlow+1
         mvec=(/(m,m=mlow,mhigh)/)
         singfac=0
         singfac(:,0)=mvec-nn*q(0)
         singfac(:,1)=-nn*q(1)
         singfac(:,2)=-nn*q(2)
         singfac(:,3)=-nn*q(3)
         singfac(ipert0,0)=-nn*q(1)
         singfac(ipert0,1)=-nn*q(2)/2
         singfac(ipert0,2)=-nn*q(3)/3
         singfac(ipert0,3)=0
c-----------------------------------------------------------------------
c     compute factored Hermitian banded matrix F and its derivatives.
c-----------------------------------------------------------------------
         f=0
         iqty=0
         DO jpert=1,mpert
            DO ipert=jpert,MIN(mpert,jpert+mband)
               iqty=iqty+1
               f(1+ipert-jpert,jpert,0)
     $              =singfac(ipert,0)*fmats_s_f(iqty)
               IF(singp%order < 1)CYCLE
               f(1+ipert-jpert,jpert,1)
     $              =singfac(ipert,0)*fmats_s_f1(iqty)*sig
     $              +singfac(ipert,1)*fmats_s_f(iqty)
               IF(singp%order < 2)CYCLE
               f(1+ipert-jpert,jpert,2)
     $              =singfac(ipert,0)*fmats_s_f2(iqty)
     $              +2*singfac(ipert,1)*fmats_s_f1(iqty)*sig
     $              +singfac(ipert,2)*fmats_s_f(iqty)
               IF(singp%order < 3)CYCLE
               f(1+ipert-jpert,jpert,3)
     $              =singfac(ipert,0)*fmats_s_f3(iqty)*sig
     $              +3*singfac(ipert,1)*fmats_s_f2(iqty)
     $              +3*singfac(ipert,2)*fmats_s_f1(iqty)*sig
     $              +singfac(ipert,3)*fmats_s_f(iqty)
               IF(singp%order < 4)CYCLE
               f(1+ipert-jpert,jpert,4)
     $              =4*singfac(ipert,1)*fmats_s_f3(iqty)*sig
     $              +6*singfac(ipert,2)*fmats_s_f2(iqty)
     $              +4*singfac(ipert,3)*fmats_s_f1(iqty)*sig
               IF(singp%order < 5)CYCLE
               f(1+ipert-jpert,jpert,5)
     $              =10*singfac(ipert,2)*fmats_s_f3(iqty)*sig
     $              +10*singfac(ipert,3)*fmats_s_f2(iqty)
               IF(singp%order < 6)CYCLE
               f(1+ipert-jpert,jpert,6)
     $              =20*singfac(ipert,3)*fmats_s_f3(iqty)*sig
            ENDDO
         ENDDO
         f0=f(:,:,0)
c-----------------------------------------------------------------------
c     compute product of factored Hermitian banded matrix F.
c-----------------------------------------------------------------------
         ff=0
         fac0=1
         DO n=0,singp%order
            fac1=1
            DO j=0,n
               DO jpert=1,mpert
                  DO ipert=jpert,MIN(mpert,jpert+mband)
                     DO kpert=MAX(1,ipert-mband),jpert
                        ff(1+ipert-jpert,jpert,n)
     $                       =ff(1+ipert-jpert,jpert,n)
     $                       +fac1*f(1+ipert-kpert,kpert,j)
     $                       *CONJG(f(1+jpert-kpert,kpert,n-j))
                     ENDDO
                  ENDDO
               ENDDO
               fac1=fac1*(n-j)/(j+1)
            ENDDO
            ff(:,:,n)=ff(:,:,n)/fac0
            fac0=fac0*(n+1)
         ENDDO
c-----------------------------------------------------------------------
c     compute non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
         k=0
         iqty=0
         DO jpert=1,mpert
            DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
               iqty=iqty+1
               k(1+mband+ipert-jpert,jpert,0)
     $              =singfac(ipert,0)*kpot_s_f(iqty)
               IF(singp%order < 1)CYCLE
               k(1+mband+ipert-jpert,jpert,1)
     $              =singfac(ipert,0)*kpot_s_f1(iqty)*sig
     $              +singfac(ipert,1)*kpot_s_f(iqty)
               IF(singp%order < 2)CYCLE
               k(1+mband+ipert-jpert,jpert,2)
     $              =singfac(ipert,0)*kpot_s_f2(iqty)/2
     $              +singfac(ipert,1)*kpot_s_f1(iqty)*sig
     $              +singfac(ipert,2)*kpot_s_f(iqty)/2
               IF(singp%order < 3)CYCLE
               k(1+mband+ipert-jpert,jpert,3)
     $              =singfac(ipert,0)*kpot_s_f3(iqty)*sig/6
     $              +singfac(ipert,1)*kpot_s_f2(iqty)/2
     $              +singfac(ipert,2)*kpot_s_f1(iqty)*sig/2
     $              +singfac(ipert,3)*kpot_s_f(iqty)/6
               IF(singp%order < 4)CYCLE
               k(1+mband+ipert-jpert,jpert,4)
     $              =singfac(ipert,1)*kpot_s_f3(iqty)*sig/6
     $              +singfac(ipert,2)*kpot_s_f2(iqty)/4
     $              +singfac(ipert,3)*kpot_s_f1(iqty)*sig/6
               IF(singp%order < 5)CYCLE
               k(1+mband+ipert-jpert,jpert,5)
     $              =singfac(ipert,2)*kpot_s_f3(iqty)*sig/12
     $              +singfac(ipert,3)*kpot_s_f2(iqty)/12
               IF(singp%order < 6)CYCLE
               k(1+mband+ipert-jpert,jpert,6)
     $              =singfac(ipert,3)*kpot_s_f3(iqty)*sig/36
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     compute Hermitian banded matrix G.
c-----------------------------------------------------------------------
         g=0
         iqty=0
         DO jpert=1,mpert
            DO ipert=jpert,MIN(mpert,jpert+mband)
               iqty=iqty+1
               g(1+ipert-jpert,jpert,0)=gpot_s_f(iqty)
               IF(singp%order < 1)CYCLE
               g(1+ipert-jpert,jpert,1)=gpot_s_f1(iqty)*sig
               IF(singp%order < 2)CYCLE
               g(1+ipert-jpert,jpert,2)=gpot_s_f2(iqty)/2
               IF(singp%order < 3)CYCLE
               g(1+ipert-jpert,jpert,3)=gpot_s_f3(iqty)*sig/6
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     compute identity.
c-----------------------------------------------------------------------
         v=0
         DO ipert=1,mpert
            v(ipert,ipert,1)=1
            v(ipert,ipert+mpert,2)=1
         ENDDO
c-----------------------------------------------------------------------
c     compute zeroth-order x1.
c-----------------------------------------------------------------------
         x=0
         DO isol=1,msol
            x(:,isol,1,0)=v(:,isol,2)
            CALL zgbmv('N',mpert,mpert,mband,mband,-one,k(:,:,0),
     $           2*mband+1,v(:,isol,1),1,one,x(:,isol,1,0),1)
         ENDDO
         CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,0),
     $        mpert,info)
c-----------------------------------------------------------------------
c     compute higher-order x1.
c-----------------------------------------------------------------------
         DO i=1,singp%order
            DO isol=1,msol
               DO j=1,i
                  CALL zhbmv('L',mpert,mband,-one,ff(:,:,j),
     $                 mband+1,x(:,isol,1,i-j),1,one,x(:,isol,1,i),1)
               ENDDO
               CALL zgbmv('N',mpert,mpert,mband,mband,-one,k(:,:,i),
     $              2*mband+1,v(:,isol,1),1,one,x(:,isol,1,i),1)
            ENDDO
            CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,i),
     $           mpert,info)
         ENDDO
c-----------------------------------------------------------------------
c     compute x2.
c-----------------------------------------------------------------------
         DO i=0,singp%order
            DO isol=1,msol
               DO j=0,i
                  CALL zgbmv('C',mpert,mpert,mband,mband,one,k(:,:,j),
     $                 2*mband+1,x(:,isol,1,i-j),1,one,x(:,isol,2,i),1)
               ENDDO
               CALL zhbmv('L',mpert,mband,one,g(:,:,i),
     $              mband+1,v(:,isol,1),1,one,x(:,isol,2,i),1)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     principal terms of mmat.
c-----------------------------------------------------------------------
         mmat=0
         j=0
         DO i=0,singp%order
            mmat(r1,r2,:,j)=x(r1,r2,:,i)
            mmat(r1,n2,:,j+1)=x(r1,n2,:,i)
            mmat(n1,r2,:,j+1)=x(n1,r2,:,i)
            mmat(n1,n2,:,j+2)=x(n1,n2,:,i)
            j=j+2
         ENDDO
c-----------------------------------------------------------------------
c     shearing terms.
c-----------------------------------------------------------------------
         mmat(r1,r2(1),1,0)=mmat(r1,r2(1),1,0)+half*sig
         mmat(r1,r2(2),2,0)=mmat(r1,r2(2),2,0)-half*sig
      ENDDO
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(f,ff,g,k,x)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_mmat
c-----------------------------------------------------------------------
c     subprogram 7. sing_solve.
c     solves iteratively for the next order in the power series vmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_solve(k,sig,power,mmat,vmat,m0mat,r1,n1)

      INTEGER, INTENT(IN) :: k
      REAL(r8), INTENT(IN) :: sig
      INTEGER, DIMENSION(:), POINTER, INTENT(IN) :: r1,n1
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: power
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(IN) :: mmat
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(INOUT) :: vmat
      COMPLEX(r8), DIMENSION(2,2), INTENT(IN) :: m0mat

      INTEGER :: l,isol
      REAL(r8), PARAMETER :: two=2
      COMPLEX(r8) :: det
      COMPLEX(r8), DIMENSION(2,2) :: a
      COMPLEX(r8), DIMENSION(2) :: x
c-----------------------------------------------------------------------
c     compute rhs.
c-----------------------------------------------------------------------
      DO l=1,k
         vmat(:,:,:,k)=vmat(:,:,:,k)
     $        +sing_matmul(mmat(:,:,:,l),vmat(:,:,:,k-l))
      ENDDO
c-----------------------------------------------------------------------
c     evaluate solutions.
c-----------------------------------------------------------------------
      DO isol=1,2*mpert
         a=m0mat
         a(1,1)=a(1,1)-sig*(k/two+power(isol))
         a(2,2)=a(2,2)-sig*(k/two+power(isol))
         det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
         x=-vmat(r1(1),isol,:,k)
         vmat(r1(1),isol,1,k)=(a(2,2)*x(1)-a(1,2)*x(2))/det
         vmat(r1(1),isol,2,k)=(a(1,1)*x(2)-a(2,1)*x(1))/det
         vmat(n1,isol,:,k)=vmat(n1,isol,:,k)*sig/(power(isol)+k/two)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_solve
c-----------------------------------------------------------------------
c     subprogram 8. sing_matmul.
c     multiplies matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sing_matmul(a,b) RESULT(c)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: a,b
      COMPLEX(r8), DIMENSION(SIZE(a,1),SIZE(b,2),2) :: c

      CHARACTER(64) :: message
      INTEGER :: i,j,m,n
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      m=SIZE(b,1)
      n=SIZE(b,2)
      IF(SIZE(a,2) /= 2*m)THEN
         WRITE(message,'(2(a,i3))')
     $        "Sing_matmul: SIZE(a,2) = ",SIZE(a,2),
     $        " /= 2*SIZE(b,1) = ",2*SIZE(b,1)
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     main computations.
c-----------------------------------------------------------------------
      c=0
      DO i=1,n
         DO j=1,2
            c(:,i,j)=c(:,i,j)
     $           +MATMUL(a(:,1:m,j),b(:,i,1))
     $           +MATMUL(a(:,1+m:2*m,j),b(:,i,2))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION sing_matmul
c-----------------------------------------------------------------------
c     subprogram 9. sing_get_ua.
c     computes asymptotic series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_get_ua(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ua

      INTEGER :: iorder
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: dpsi,sqrtfac
      COMPLEX(r8) :: pfac
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointers and compute distance from singular surface.
c-----------------------------------------------------------------------
      singp => sing(ising)
      IF(psifac < singp%psifac)THEN
         vmat => singp%vmatl
         dpsi=singp%psifac-psifac
      ELSE
         vmat => singp%vmatr
         dpsi=psifac-singp%psifac
      ENDIF
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute powers of distance from singular surface.
c-----------------------------------------------------------------------
      sqrtfac=SQRT(dpsi)
      pfac=dpsi**singp%alpha
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      ua=vmat(:,:,:,2*singp%order)
      DO iorder=2*singp%order-1,0,-1
         ua=ua*sqrtfac+vmat(:,:,:,iorder)
      ENDDO
c-----------------------------------------------------------------------
c     restore powers (unshear v->u).
c-----------------------------------------------------------------------
      ua(r1,:,1)=ua(r1,:,1)/sqrtfac    !<---mult row by shear z^(-1/2)
      ua(r1,:,2)=ua(r1,:,2)*sqrtfac    !<---mult row by shear z^(+1/2)
      ua(:,r2(1),:)=ua(:,r2(1),:)/pfac !<---mult column by z^alpha-
      ua(:,r2(2),:)=ua(:,r2(2),:)*pfac !<---mult column by z^alpha+
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_get_ua
c-----------------------------------------------------------------------
c     subprogram 10. sing_der.
c     evaluates the differential equations of DCON.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_der(neq,psifac,u,du)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,mpert,2), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(mpert,mpert,2), INTENT(OUT) :: du

      INTEGER :: ipert,jpert,isol,iqty,info
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb,gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psifac,0)
      CALL cspline_eval(fmats,psifac,0)
      CALL cspline_eval(gmats,psifac,0)
      CALL cspline_eval(kmats,psifac,0)
c-----------------------------------------------------------------------
c     define local scalars.
c-----------------------------------------------------------------------
      q=sq%f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
      singfac=1/singfac
c-----------------------------------------------------------------------
c     copy Hermitian banded matrices F and G.
c-----------------------------------------------------------------------
      fmatb=0
      gmatb=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            fmatb(1+ipert-jpert,jpert)=fmats%f(iqty)
            gmatb(1+ipert-jpert,jpert)=gmats%f(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
      kmatb=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            kmatb(1+mband+ipert-jpert,jpert)=kmats%f(iqty)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute du1.
c-----------------------------------------------------------------------
      du=0
      DO isol=1,mpert
         du(:,isol,1)=u(:,isol,2)*singfac
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,u(:,isol,1),1,one,du(:,isol,1),1)
      ENDDO
      CALL zpbtrs('L',mpert,mband,mpert,fmatb,mband+1,du,mpert,info)
c-----------------------------------------------------------------------
c     compute du2.
c-----------------------------------------------------------------------
      DO isol=1,mpert
         CALL zhbmv('L',mpert,mband,one,gmatb,
     $        mband+1,u(:,isol,1),1,one,du(:,isol,2),1 )
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,du(:,isol,1),1,one,du(:,isol,2),1)
         du(:,isol,1)=du(:,isol,1)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_der
c-----------------------------------------------------------------------
c     subprogram 11. sing_derFM.
c     evaluates the differential equations of DCON
c     using Fundametal Matrix of Solutions
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_derFM(neqFM,startPsi,uFM,duFM,rpar,ipar)

      INTEGER, INTENT(IN) :: neqFM
      REAL(r8), INTENT(IN) :: startPsi
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert), INTENT(IN) :: uFM
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert), INTENT(OUT) :: duFM
      COMPLEX(r8), INTENT(IN) :: rpar, ipar

      INTEGER :: ipert,jpert,isol,iqty,info
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1
      COMPLEX(r8), DIMENSION(mpert,mpert) :: fmatb,gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
      COMPLEX(r8), DIMENSION(mpert,mpert) :: duTemp
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      CALL spline_eval_external(sq,startPsi,sq_s_ix,sq_s_f)
      CALL cspline_eval_external(fmats,startPsi,
     $     fmats_s_ix,fmats_s_f)
      CALL cspline_eval_external(gmats,startPsi,
     $     gmats_s_ix,gmats_s_f)
      CALL cspline_eval_external(kmats,startPsi,
     $     kmats_s_ix,kmats_s_f)
c-----------------------------------------------------------------------
c     define local scalars.
c-----------------------------------------------------------------------
      q=sq_s_f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
      singfac=1/singfac
c-----------------------------------------------------------------------
c     copy Hermitian banded matrices F and G.
c-----------------------------------------------------------------------
      fmatb=0
      gmatb=0
      DO jpert=1,mpert
         DO ipert=jpert,mpert
            iqty=1+(jpert-1)*(mpert+1)-jpert*(jpert-1)/2+(ipert-jpert)
            fmatb(1+ipert-jpert,jpert)=fmats_s_f(iqty)
            gmatb(1+ipert-jpert,jpert)=gmats_s_f(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
      kmatb=0
      DO jpert=1,mpert
         DO ipert=1,mpert
            iqty=1+(jpert-1)*mpert+(ipert-1)
            kmatb(1+mband+ipert-jpert,jpert)=kmats_s_f(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = [-K*u11+Q^{-1}*u21 | -K*u12+Q^{-1}*u22]
c-----------------------------------------------------------------------
      duFM=0
      DO isol=1,mpert
         duFM(1:mpert,isol)=uFM(mpert+1:2*mpert,isol)*singfac
         duFM(1:mpert,isol+mpert)=uFM(mpert+1:2*mpert,isol+mpert)*
     $        singfac
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,uFM(1:mpert,isol),1,one,duFM(1:mpert,isol),1)
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,uFM(1:mpert,isol+mpert),1,one,
     $        duFM(1:mpert,isol+mpert),1)
      ENDDO
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = [F^{-1}*du11 | F^{-1}*du12]
c-----------------------------------------------------------------------
      duTemp = duFM(1:mpert,1:mpert)
      CALL zpbtrs('L',mpert,mband,mpert,fmatb,mband+1,duTemp,mpert,info)
      duFM(1:mpert,1:mpert) = duTemp
      duTemp = duFM(1:mpert,mpert+1:2*mpert)
      CALL zpbtrs('L',mpert,mband,mpert,fmatb,mband+1,duTemp,mpert,info)
      duFM(1:mpert,mpert+1:2*mpert) = duTemp
c-----------------------------------------------------------------------
c     Compute: Final loop.
c-----------------------------------------------------------------------
      DO isol=1,mpert
c-----------------------------------------------------------------------
c     Compute: du21 = du21(==0) + G*u11
c-----------------------------------------------------------------------
         CALL zhbmv('L',mpert,mband,one,gmatb,
     $        mband+1,uFM(1:mpert,isol),1,one,
     $        duFM(mpert+1:2*mpert,isol),1)
c-----------------------------------------------------------------------
c     Compute: du22 = du22(==0) + G*u12
c-----------------------------------------------------------------------
         CALL zhbmv('L',mpert,mband,one,gmatb,
     $        mband+1,uFM(1:mpert,isol+mpert),1,one,
     $        duFM(mpert+1:2*mpert,isol+mpert),1)
c-----------------------------------------------------------------------
c     Compute: du21 = du21 + Kt*du11
c-----------------------------------------------------------------------
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,duFM(1:mpert,isol),1,one,
     $        duFM(mpert+1:2*mpert,isol),1)
c-----------------------------------------------------------------------
c     Compute: du22 = du22 + Kt*du12
c-----------------------------------------------------------------------
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,duFM(1:mpert,isol+mpert),1,one,
     $        duFM(mpert+1:2*mpert,isol+mpert),1)
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = Q^{-1} * [du11 | du12]
c-----------------------------------------------------------------------
         duFM(1:mpert,isol)=duFM(1:mpert,isol)*singfac
         duFM(1:mpert,isol+mpert)=duFM(1:mpert,isol+mpert)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_derFM
c-----------------------------------------------------------------------
c     subprogram 12. sing_derL.
c     evaluates the differential equations of DCON: U'=LU
c     using FM of Sol'ns. Gives output L, (not L*x as in x'=Lx).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing_derL(startPsi,L)

      REAL(r8), INTENT(IN) :: startPsi
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert), INTENT(OUT) :: L

      INTEGER :: ipert,jpert,isol,iqty,info
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1
      COMPLEX(r8), DIMENSION(mpert,mpert) :: fmatb,gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
      COMPLEX(r8), DIMENSION(mpert,mpert) :: dTemp
      COMPLEX(r8), DIMENSION(mpert,mpert) :: idMat
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      CALL spline_eval_external(sq,startPsi,sq_s_ix,sq_s_f)
      CALL cspline_eval_external(fmats,startPsi,
     $     fmats_s_ix,fmats_s_f)
      CALL cspline_eval_external(gmats,startPsi,
     $     gmats_s_ix,gmats_s_f)
      CALL cspline_eval_external(kmats,startPsi,
     $     kmats_s_ix,kmats_s_f)
c-----------------------------------------------------------------------
c     define local scalars.
c-----------------------------------------------------------------------
      q=sq_s_f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
      singfac=1/singfac
c-----------------------------------------------------------------------
c     copy Hermitian banded matrices F and G.
c-----------------------------------------------------------------------
      fmatb=0
      gmatb=0
      DO jpert=1,mpert
         DO ipert=jpert,mpert
            iqty=1+(jpert-1)*(mpert+1)-jpert*(jpert-1)/2+(ipert-jpert)
            fmatb(1+ipert-jpert,jpert)=fmats_s_f(iqty)
            gmatb(1+ipert-jpert,jpert)=gmats_s_f(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     copy non-Hermitian banded matrix K.
c-----------------------------------------------------------------------
      kmatb=0
      DO jpert=1,mpert
         DO ipert=1,mpert
            iqty=1+(jpert-1)*mpert+(ipert-1)
            kmatb(1+mband+ipert-jpert,jpert)=kmats_s_f(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = [-bar(K) | Q^{-1}]
c-----------------------------------------------------------------------
      L=0
      idMat=0
      DO isol=1,mpert
         idMat(isol,isol) = 1.0_r8
         L(isol,mpert+isol) = singfac(isol)
      ENDDO
      DO isol=1,mpert
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,idMat(:,isol),1,zero,L(1:mpert,isol),1)
      ENDDO
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = [bar(F)^{-1}*du11 | bar(F)^{-1}*du12]
c-----------------------------------------------------------------------
      dTemp = L(1:mpert,1:mpert)
      CALL zpbtrs('L',mpert,mband,mpert,fmatb,mband+1,dTemp,mpert,info)
      L(1:mpert,1:mpert) = dTemp
      dTemp = L(1:mpert,mpert+1:2*mpert)
      CALL zpbtrs('L',mpert,mband,mpert,fmatb,mband+1,dTemp,mpert,info)
      L(1:mpert,mpert+1:2*mpert) = dTemp
c-----------------------------------------------------------------------
c     Compute: Final loop.
c-----------------------------------------------------------------------
      DO isol=1,mpert
c-----------------------------------------------------------------------
c     Compute: du21 = du21(==0) + G
c-----------------------------------------------------------------------
         CALL zhbmv('L',mpert,mband,one,gmatb,mband+1,idMat(:,isol),1,
     $        zero,L(mpert+1:2*mpert,isol),1)
c-----------------------------------------------------------------------
c     Compute: du21 = du21 + Kt*du11
c-----------------------------------------------------------------------
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,L(1:mpert,isol),1,one,
     $        L(mpert+1:2*mpert,isol),1)
c-----------------------------------------------------------------------
c     Compute: du22 = du22 + Kt*du12
c-----------------------------------------------------------------------
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,L(1:mpert,isol+mpert),1,one,
     $        L(mpert+1:2*mpert,isol+mpert),1)
c-----------------------------------------------------------------------
c     Compute: [du11 | du12] = Q^{-1} * [du11 | du12]
c-----------------------------------------------------------------------
         L(1:mpert,isol)=L(1:mpert,isol)*singfac
         L(1:mpert,isol+mpert)=L(1:mpert,isol+mpert)*singfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing_derL

      END MODULE sing_mod
