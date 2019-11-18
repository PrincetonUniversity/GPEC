c-----------------------------------------------------------------------
c     file sing1.f.
c     computations relating to singular surfaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. sing1_mod.
c     1. sing1_mmat.
c     2. sing1_matmul.
c     3. sing1_solve.
c     4. sing1_vmat.
c     5. sing1_get_ua.
c     6. sing1_get_dua.
c     7. sing1_der.
c     8. sing1_matvec.
c     9. sing1_log.
c     10. sing1_write_mat.
c     11. sing1_vmat_diagnose.
c     12. sing1_ua_diagnose.
c     13. sing1_vmat_d.
c     14. sing1_solve_d.
c     15. sing1_get_ua_d.
c     16. sing1_get_dua_d.
c     17. sing1_get_ua_d_cut.
c     18. sing1_get_ua_cut.
c     19. sing1_delta.
c     20. sing1_kxscan.
c     21. sing1_xmin.
c-----------------------------------------------------------------------
c     subprogram 0. sing1_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sing1_mod
      USE fourfit_mod
      IMPLICIT NONE

      INTEGER :: msol,sing_order=2,sing_order_save=2
      INTEGER, DIMENSION(:), POINTER :: r1,r2,n1,n2
      COMPLEX(r8), DIMENSION(2,2) :: m0mat

      LOGICAL :: sing1_flag=.FALSE.,sing_order_ceiling =.TRUE.
      REAL(r8) :: degen_tol=5e-5

      TYPE :: ua_diagnose_type
      LOGICAL :: flag,neg,phase
      INTEGER :: ising,ipert,isol,iqty
      REAL(r8) :: zlogmin,zlogmax,eps
      END TYPE ua_diagnose_type

      TYPE(ua_diagnose_type) :: uad
      NAMELIST/ua_diagnose_list/uad

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. sing1_mmat.
c     computes series expansion of coefficient matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_mmat(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER :: ipert0,ipert,jpert,kpert,isol,iqty,info,msol,m,i,j,n,
     $     fac0,fac1
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: psifac
      REAL(r8), DIMENSION(0:3) :: q
      REAL(r8), DIMENSION(mpert,0:3) :: singfac
      COMPLEX(r8), PARAMETER :: one=1,half=.5_r8
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: f0
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: v
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: mmat

      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: f,ff,g,k
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: x

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      TYPE(sing_type), POINTER :: singp
      TYPE(cspline_type), POINTER :: kpot,gpot
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"io",4x,"is",4x,"ip",2(4x,"re x",i1,6x,"im x",i1,2x)/)
 20   FORMAT(3i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     evaluate cubic splines.
c-----------------------------------------------------------------------
      IF(1) THEN
         kpot=>kmats
         gpot=>gmats
      ELSE
         kpot=>kmatsp
         gpot=>gmatsp
      ENDIF
      msol=2*mpert
      singp => sing(ising)
      psifac=singp%psifac
      CALL spline_eval(sq,psifac,3)
      CALL cspline_eval(fmats,psifac,3)
      CALL cspline_eval(gpot,psifac,3)
      CALL cspline_eval(kpot,psifac,3)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(f(mband+1,mpert,0:singp%order))
      ALLOCATE(ff(mband+1,mpert,0:singp%order))
      ALLOCATE(g(mband+1,mpert,0:singp%order))
      ALLOCATE(k(2*mband+1,mpert,0:singp%order))
      ALLOCATE(x(mpert,2*mpert,2,0:singp%order))
      ALLOCATE(singp%mmat(mpert,2*mpert,2,0:2*singp%order+2))
c-----------------------------------------------------------------------
c     evaluate safety factor and its derivatives.
c-----------------------------------------------------------------------
      mmat => singp%mmat
      q(0)=sq%f(4)
      q(1)=sq%f1(4)
      q(2)=sq%f2(4)
      q(3)=sq%f3(4)
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
     $           =singfac(ipert,0)*fmats%f(iqty)
            IF(singp%order < 1)CYCLE
            f(1+ipert-jpert,jpert,1)
     $           =singfac(ipert,0)*fmats%f1(iqty)
     $           +singfac(ipert,1)*fmats%f(iqty)
            IF(singp%order < 2)CYCLE
            f(1+ipert-jpert,jpert,2)
     $           =singfac(ipert,0)*fmats%f2(iqty)
     $           +2*singfac(ipert,1)*fmats%f1(iqty)
     $           +singfac(ipert,2)*fmats%f(iqty)
            IF(singp%order < 3)CYCLE
            f(1+ipert-jpert,jpert,3)
     $           =singfac(ipert,0)*fmats%f3(iqty)
     $           +3*singfac(ipert,1)*fmats%f2(iqty)
     $           +3*singfac(ipert,2)*fmats%f1(iqty)
     $           +singfac(ipert,3)*fmats%f(iqty)
            IF(singp%order < 4)CYCLE
            f(1+ipert-jpert,jpert,4)
     $           =4*singfac(ipert,1)*fmats%f3(iqty)
     $           +6*singfac(ipert,2)*fmats%f2(iqty)
     $           +4*singfac(ipert,3)*fmats%f1(iqty)
            IF(singp%order < 5)CYCLE
            f(1+ipert-jpert,jpert,5)
     $           =10*singfac(ipert,2)*fmats%f3(iqty)
     $           +10*singfac(ipert,3)*fmats%f2(iqty)
            IF(singp%order < 6)CYCLE
            f(1+ipert-jpert,jpert,5)
     $           =20*singfac(ipert,3)*fmats%f3(iqty)
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
     $                    =ff(1+ipert-jpert,jpert,n)
     $                    +fac1*f(1+ipert-kpert,kpert,j)
     $                    *CONJG(f(1+jpert-kpert,kpert,n-j))
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
     $           =singfac(ipert,0)*kpot%f(iqty)
            IF(singp%order < 1)CYCLE
            k(1+mband+ipert-jpert,jpert,1)
     $           =singfac(ipert,0)*kpot%f1(iqty)
     $           +singfac(ipert,1)*kpot%f(iqty)
            IF(singp%order < 2)CYCLE
            k(1+mband+ipert-jpert,jpert,2)
     $           =singfac(ipert,0)*kpot%f2(iqty)/2
     $           +singfac(ipert,1)*kpot%f1(iqty)
     $           +singfac(ipert,2)*kpot%f(iqty)/2
            IF(singp%order < 3)CYCLE
            k(1+mband+ipert-jpert,jpert,3)
     $           =singfac(ipert,0)*kpot%f3(iqty)/6
     $           +singfac(ipert,1)*kpot%f2(iqty)/2
     $           +singfac(ipert,2)*kpot%f1(iqty)/2
     $           +singfac(ipert,3)*kpot%f(iqty)/6
            IF(singp%order < 4)CYCLE
            k(1+mband+ipert-jpert,jpert,4)
     $           =singfac(ipert,1)*kpot%f3(iqty)/6
     $           +singfac(ipert,2)*kpot%f2(iqty)/4
     $           +singfac(ipert,3)*kpot%f1(iqty)/6
            IF(singp%order < 5)CYCLE
            k(1+mband+ipert-jpert,jpert,5)
     $           =singfac(ipert,2)*kpot%f3(iqty)/12
     $           +singfac(ipert,3)*kpot%f2(iqty)/12
            IF(singp%order < 6)CYCLE
            k(1+mband+ipert-jpert,jpert,6)
     $           =singfac(ipert,3)*kpot%f3(iqty)/36
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
            g(1+ipert-jpert,jpert,0)=gpot%f(iqty)
            IF(singp%order < 1)CYCLE
            g(1+ipert-jpert,jpert,1)=gpot%f1(iqty)
            IF(singp%order < 2)CYCLE
            g(1+ipert-jpert,jpert,2)=gpot%f2(iqty)/2
            IF(singp%order < 3)CYCLE
            g(1+ipert-jpert,jpert,3)=gpot%f3(iqty)/6
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
     $        2*mband+1,v(:,isol,1),1,one,x(:,isol,1,0),1)
      ENDDO
      CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,0),
     $     mpert,info)
c-----------------------------------------------------------------------
c     compute higher-order x1.
c-----------------------------------------------------------------------
      DO i=1,singp%order
         DO isol=1,msol
            DO j=1,i
               CALL zhbmv('L',mpert,mband,-one,ff(:,:,j),
     $              mband+1,x(:,isol,1,i-j),1,one,x(:,isol,1,i),1)
            ENDDO
            CALL zgbmv('N',mpert,mpert,mband,mband,-one,k(:,:,i),
     $           2*mband+1,v(:,isol,1),1,one,x(:,isol,1,i),1)
         ENDDO
         CALL zpbtrs('L',mpert,mband,msol,f0,mband+1,x(:,:,1,i),
     $        mpert,info)
      ENDDO
c-----------------------------------------------------------------------
c     compute x2.
c-----------------------------------------------------------------------
      DO i=0,singp%order
         DO isol=1,msol
            DO j=0,i
               CALL zgbmv('C',mpert,mpert,mband,mband,one,k(:,:,j),
     $              2*mband+1,x(:,isol,1,i-j),1,one,x(:,isol,2,i),1)
            ENDDO
            CALL zhbmv('L',mpert,mband,one,g(:,:,i),
     $           mband+1,v(:,isol,1),1,one,x(:,isol,2,i),1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(debug_unit,"xmat.out","UNKNOWN")
         DO j=0,singp%order
            DO isol=1,2*mpert
               WRITE(debug_unit,10)(i,i,i=1,2)
               DO ipert=1,mpert
                  WRITE(debug_unit,20)j,isol,ipert,
     $                 x(ipert,isol,:,j)
               ENDDO
            ENDDO
         ENDDO
         WRITE(debug_unit,10)(i,i,i=1,2)
         CALL ascii_close(debug_unit)
      ENDIF
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
      mmat(r1,r2(1),1,0)=mmat(r1,r2(1),1,0)+half
      mmat(r1,r2(2),2,0)=mmat(r1,r2(2),2,0)-half
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(f,ff,g,k,x)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_mmat
c-----------------------------------------------------------------------
c     subprogram 2. sing1_matmul.
c     multiplies matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sing1_matmul(mat,vec) RESULT(matvec)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: mat,vec
      COMPLEX(r8), DIMENSION(SIZE(mat,1),SIZE(vec,2),2) :: matvec

      CHARACTER(64) :: message
      INTEGER :: isol,j,mpert,msol
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      mpert=SIZE(vec,1)
      msol=SIZE(vec,2)
c-----------------------------------------------------------------------
c     mismatch abort.
c-----------------------------------------------------------------------
      IF(SIZE(mat,2) /= 2*mpert)THEN
         WRITE(message,'(2(a,i3))')
     $        "sing1_matmul: SIZE(mat,2) = ",SIZE(mat,2),
     $        " /= 2*SIZE(vec,1) = ",2*SIZE(vec,1)
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     main computations.
c-----------------------------------------------------------------------
         matvec=0
         DO isol=1,msol
            DO j=1,2
               matvec(:,isol,j)=matvec(:,isol,j)
     $              +MATMUL(mat(:,1:mpert,j),vec(:,isol,1))
     $              +MATMUL(mat(:,mpert+1:2*mpert,j),vec(:,isol,2))
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
         RETURN
         END FUNCTION sing1_matmul
c-----------------------------------------------------------------------
c     subprogram 3. sing1_solve.
c     solves iteratively for the next order in the power series vmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_solve(k,power,mmat,vmat)

      INTEGER, INTENT(IN) :: k
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: power
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(IN) :: mmat
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,0:k), INTENT(INOUT) :: vmat

      INTEGER :: l,isol,info
      INTEGER, DIMENSION(2) :: ipiv
      REAL(r8), PARAMETER :: two=2
      COMPLEX(r8), DIMENSION(2,2) :: a
      COMPLEX(r8), DIMENSION(2) :: x
c-----------------------------------------------------------------------
c     compute rhs.
c-----------------------------------------------------------------------
      DO l=1,k
         vmat(:,:,:,k)=vmat(:,:,:,k)
     $        +sing1_matmul(mmat(:,:,:,l),vmat(:,:,:,k-l))
      ENDDO
c-----------------------------------------------------------------------
c     solve.
c-----------------------------------------------------------------------
      DO isol=1,2*mpert
         a=m0mat-(k/two+power(isol))*RESHAPE((/1,0,0,1/),(/2,2/))
         x=-vmat(r1(1),isol,:,k)
         CALL zgetrf(2,2,a,2,ipiv,info)
         CALL zgetrs('N',2,1,a,2,ipiv,x,2,info)
         vmat(r1(1),isol,:,k)=x
         vmat(n1,isol,:,k)=vmat(n1,isol,:,k)/(power(isol)+k/two)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_solve
c-----------------------------------------------------------------------
c     subprogram 4. sing1_vmat.
c     computes asymptotic behavior at the singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_vmat(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER :: ipert0,ipert,k
      REAL(r8) :: psifac,di,di0,q,q1,rho
      REAL(r8), PARAMETER :: half=.5_r8
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat,mmat
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: amat,bmat
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
      IF(ipert0 <= 0 .OR. mlow > nn*q .OR. mhigh < nn*q)THEN
         singp%di=0
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     allocate and compute ranges.
c-----------------------------------------------------------------------
      ALLOCATE(sing(ising)%n1(mpert-1),sing(ising)%n2(2*mpert-2))
      q=singp%q
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
      CALL spline_eval(locstab,psifac,0)
      di0=locstab%f(1)/psifac
      q1=singp%q1
      rho=singp%rho
c-----------------------------------------------------------------------
c     compute Mercier index.
c-----------------------------------------------------------------------
      singp%order=0
      CALL sing1_mmat(ising)
      m0mat=TRANSPOSE(singp%mmat(r1(1),r2,:,0))
      di=m0mat(1,1)*m0mat(2,2)-m0mat(2,1)*m0mat(1,2)
      singp%di=di
      singp%alpha=SQRT(-CMPLX(singp%di))
c      singp%degen=.TRUE.
      singp%degen=ABS(di+0.25_r8) < degen_tol
c$$$      WRITE(*,*)'SINGP%DEGEN=',singp%degen,'di+0.25=',di+0.25_r8
c$$$      WRITE(*,*) 'degen_tol=',degen_tol
c-----------------------------------------------------------------------
c     reset order and compute higher-order terms.
c-----------------------------------------------------------------------
      singp%order=sing_order
      IF(sing_order_ceiling)
     $     singp%order=singp%order+CEILING(2*REAL(singp%alpha))
      DEALLOCATE(singp%mmat)
      CALL sing1_mmat(ising)
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      ALLOCATE(singp%power(2*mpert))
      singp%power=0
      singp%power(ipert0)=-singp%alpha
      singp%power(ipert0+mpert)=singp%alpha
c-----------------------------------------------------------------------
c     write to table of singular surface quantities.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(i3,1p,7e11.3)')ising,psifac,rho,q,q1,di0,
     $     singp%di,singp%di/di0-1
c-----------------------------------------------------------------------
c     allocate vmat and set pointers.
c-----------------------------------------------------------------------
      ALLOCATE(singp%vmat(mpert,2*mpert,2,0:2*singp%order))
      vmat => singp%vmat
      mmat => singp%mmat
      m0mat=TRANSPOSE(singp%mmat(r1(1),r2,:,0))
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
     $     =-(m0mat(1,1)+singp%alpha)/m0mat(1,2)
      vmat(ipert0,ipert0+mpert,2,0)
     $     =-(m0mat(1,1)-singp%alpha)/m0mat(1,2)
c-----------------------------------------------------------------------
c     compute higher-order solutions.
c-----------------------------------------------------------------------
      DO k=1,2*singp%order
         CALL sing1_solve(k,singp%power,mmat,vmat)
      ENDDO
c-----------------------------------------------------------------------
c     compute degenerate solutions.
c-----------------------------------------------------------------------
      IF(singp%degen)THEN
         ALLOCATE(amat(mpert,2*mpert,2,0:2*singp%order))
         ALLOCATE(bmat(mpert,2*mpert,2,0:2*singp%order))
         CALL sing1_vmat_d(ising,amat,bmat)
      ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(uad%flag .AND. ising == uad%ising)THEN
         IF(singp%degen)THEN
            CALL sing1_vmat_diagnose(ising,mmat,vmat,amat,bmat)
            DEALLOCATE(amat,bmat)
         ELSE
            CALL sing1_vmat_diagnose(ising,mmat,vmat)
         ENDIF
         CALL sing1_ua_diagnose(ising)
         CALL program_stop
     $        ("sing1_vmat: abort after ua_diagnose")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_vmat
c-----------------------------------------------------------------------
c     subprogram 5. sing1_get_ua.
c     computes asymptotic series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_ua(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ua

      INTEGER :: iorder,isol
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: pfac,z,sqrtfac,phase
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      singp => sing(ising)
      vmat => singp%vmat
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      ua=vmat(:,:,:,2*singp%order)
      DO iorder=2*singp%order-1,0,-1
         ua=ua*sqrtfac+vmat(:,:,:,iorder)
      ENDDO
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      pfac=z**singp%alpha
      ua(r1,:,1)=ua(r1,:,1)/sqrtfac
      ua(r1,:,2)=ua(r1,:,2)*sqrtfac
      ua(:,r2(1),:)=ua(:,r2(1),:)/pfac
      ua(:,r2(2),:)=ua(:,r2(2),:)*pfac
c-----------------------------------------------------------------------
c     compute degenerate solutions.
c-----------------------------------------------------------------------
      IF(singp%degen)CALL sing1_get_ua_d(ising,psifac,ua)
c-----------------------------------------------------------------------
c     adjust phase.
c-----------------------------------------------------------------------
      IF(psifac < singp%psifac .AND. uad%phase)THEN
         phase=EXP(ifac*pi*(.5+singp%alpha))
         DO isol=1,2
            ua(:,r2(isol),:)=ua(:,r2(isol),:)*phase
            phase=EXP(ifac*pi*(.5-singp%alpha))
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_ua
c-----------------------------------------------------------------------
c     subprogram 6. sing1_get_dua.
c     computes asymptotic series derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_dua(ising,psifac,dua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: dua

      INTEGER :: iorder,ipert0,isol
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: pfac,z,sqrtfac,phase
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: power
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      singp => sing(ising)
      vmat => singp%vmat
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      ipert0=NINT(nn*singp%q)-mlow+1
      power=2*singp%order
      power(r1,:,1)=power(r1,:,1)-1
      power(r1,:,2)=power(r1,:,2)+1
      power(:,r2(1),:)=power(:,r2(1),:)-2*singp%alpha
      power(:,r2(2),:)=power(:,r2(2),:)+2*singp%alpha
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      dua=vmat(:,:,:,2*singp%order)*power
      DO iorder=2*singp%order-1,0,-1
         power=power-1
         dua=dua*sqrtfac+vmat(:,:,:,iorder)*power
      ENDDO
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      pfac=z**singp%alpha
      dua(r1,:,1)=dua(r1,:,1)/sqrtfac
      dua(r1,:,2)=dua(r1,:,2)*sqrtfac
      dua(:,r2(1),:)=dua(:,r2(1),:)/pfac
      dua(:,r2(2),:)=dua(:,r2(2),:)*pfac
      dua=dua/(2*z)
c-----------------------------------------------------------------------
c     compute degenerate solutions.
c-----------------------------------------------------------------------
      IF(singp%degen)CALL sing1_get_dua_d(ising,psifac,dua)
c-----------------------------------------------------------------------
c     adjust phase.
c-----------------------------------------------------------------------
      IF(psifac < singp%psifac .AND. uad%phase)THEN
         phase=EXP(ifac*pi*(.5+singp%alpha))
         DO isol=1,2
            dua(:,r2(isol),:)=dua(:,r2(isol),:)*phase
            phase=EXP(ifac*pi*(.5-singp%alpha))
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_dua
c-----------------------------------------------------------------------
c     subprogram 7. sing1_der.
c     evaluates the differential equations of DCON.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_der(neq,psifac,u,du)

      INTEGER, INTENT(INOUT) :: neq
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(OUT) :: du

      INTEGER :: ipert,jpert,isol,iqty,info
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: fmatb,gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      neq=neq
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
      DO isol=1,msol
         du(:,isol,1)=u(:,isol,2)*singfac
         CALL zgbmv('N',mpert,mpert,mband,mband,-one,kmatb,
     $        2*mband+1,u(:,isol,1),1,one,du(:,isol,1),1)
      ENDDO
      CALL zpbtrs('L',mpert,mband,msol,fmatb,mband+1,du,mpert,info)
c-----------------------------------------------------------------------
c     compute du2.
c-----------------------------------------------------------------------
      DO isol=1,msol
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
      END SUBROUTINE sing1_der
c-----------------------------------------------------------------------
c     subprogram 8. sing1_matvec.
c     computes -(F u' + K u)' + (K^\dagger u' + G u)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_matvec(psifac,u,du,term1,term2,term3,matvec)

      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(IN) :: u,du
      COMPLEX(r8), DIMENSION(mpert,msol), INTENT(OUT) ::
     $     term1,term2,term3,matvec

      INTEGER :: ipert,jpert,isol,iqty
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), PARAMETER :: one=1,zero=0
      COMPLEX(r8), DIMENSION(mband+1,mpert) :: gmatb
      COMPLEX(r8), DIMENSION(2*mband+1,mpert) :: kmatb
c-----------------------------------------------------------------------
c     cubic spline evaluation.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psifac,0)
      CALL cspline_eval(gmats,psifac,0)
      CALL cspline_eval(kmats,psifac,0)
c-----------------------------------------------------------------------
c     define local scalars.
c-----------------------------------------------------------------------
      q=sq%f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
c-----------------------------------------------------------------------
c     copy Hermitian banded matrix G.
c-----------------------------------------------------------------------
      gmatb=0
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
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
c     compute matvec.
c-----------------------------------------------------------------------
      DO isol=1,msol
         matvec(:,isol)=du(:,isol,1)*singfac
         CALL zgbmv('C',mpert,mpert,mband,mband,one,kmatb,
     $        2*mband+1,matvec(:,isol),1,zero,term1(:,isol),1)
         CALL zhbmv('L',mpert,mband,one,gmatb,
     $        mband+1,u(:,isol,1),1,zero,term2(:,isol),1 )
      ENDDO
      term3=-du(:,:,2)
      matvec=term1+term2+term3
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_matvec
c-----------------------------------------------------------------------
c     subprogram 9. sing1_log.
c     returns bounded log of specific component.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sing1_log(u) RESULT(ulog)

      COMPLEX(r8), INTENT(IN) :: u
      COMPLEX(r4) :: ulog

      REAL(R8), PARAMETER :: eps=1e-1
      REAL(r4), PARAMETER :: minlog=-12*alog10
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(u == 0)THEN
         ulog=minlog
      ELSE
         ulog=LOG(u)
         IF(REAL(ulog) < minlog)ulog=CMPLX(minlog,IMAG(ulog))
         IF(IMAG(ulog) < -eps)ulog=CMPLX(REAL(ulog),IMAG(ulog)+twopi)
      ENDIF
      ulog=CMPLX(REAL(ulog)/alog10,IMAG(ulog))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION sing1_log
c-----------------------------------------------------------------------
c     subprogram 10. sing1_write_mat.
c     writes one matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_write_mat(r1,mat_in,name,unit)

      INTEGER :: r1
      COMPLEX(r8), DIMENSION(:,:,:,0:), INTENT(IN) :: mat_in
      CHARACTER(*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: unit

      LOGICAL :: degen
      INTEGER :: isol,msol,iorder,norder,ipert,m,iqty,i
      INTEGER, DIMENSION(2) :: r2
      REAL(r8), PARAMETER :: tol=1e-15
      COMPLEX(r8) :: detfac
      COMPLEX(r8), DIMENSION(2,2) :: mat2
      COMPLEX(r8), DIMENSION(SIZE(mat_in,1),SIZE(mat_in,2),
     $     SIZE(mat_in,3),0:SIZE(mat_in,4)-1) :: mat
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"i/j",2(2x,"re ",a4,"2",i1,2x,"im ",a4,"2",i1)/)
 20   FORMAT(i6,1p,4e11.3)
 30   FORMAT(/2x,"is",2x,"ip",3x,"m",2(3x,"re ",a4,i1,3x,"im ",a4,i1)/)
 40   FORMAT(3i4,1p,4e11.3)
c-----------------------------------------------------------------------
c     set sizes.
c-----------------------------------------------------------------------
      norder=SIZE(mat,4)-1
      msol=SIZE(mat,2)
c-----------------------------------------------------------------------
c     copy matrix.
c-----------------------------------------------------------------------
      mat=mat_in
      WHERE(ABS(IMAG(mat_in)) < tol*ABS(REAL(mat_in)))
         mat=REAL(mat_in)
      ELSEWHERE
         mat=mat_in
      END WHERE
c-----------------------------------------------------------------------
c     compute mmat2 and detfac.
c-----------------------------------------------------------------------
      r2=(/r1,r1+mpert/)
      mat2=TRANSPOSE(mat(r1,r2,:,0))
      detfac=mat2(1,1)*mat2(2,2)-mat2(1,2)*mat2(2,1)+.25_r8
      degen=ABS(detfac) < degen_tol
c-----------------------------------------------------------------------
c     diagnose mat2 and detfac.
c-----------------------------------------------------------------------
      WRITE(debug_unit,'(a)')TRIM(name)//"2:"
      WRITE(debug_unit,10)(TRIM(name),iqty,TRIM(name),iqty,iqty=1,2)
      WRITE(debug_unit,20)(i,mat2(i,:),i=1,2)
      WRITE(debug_unit,10)(TRIM(name),iqty,TRIM(name),iqty,iqty=1,2)
      WRITE(debug_unit,'(1p,2(a,e10.3),a,l1/)')
     $     " ABS(detfac) = ",ABS(detfac),", degen_tol = ",degen_tol,
     $     ", degen = ",degen
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      WRITE(unit,'(a/)')TRIM(name)//":"
      DO iorder=0,norder
         WRITE(debug_unit,'(a,i2,a,1p,e10.3)')" iorder = ",iorder,
     $        ", MAXVAL("//TRIM(name)//") = ",
     $        MAXVAL(ABS(mat(:,:,:,iorder)))
         WRITE(debug_unit,30)(TRIM(name),iqty,TRIM(name),iqty,iqty=1,2)
         DO isol=1,msol
            m=mlow
            DO ipert=1,mpert
               WRITE(debug_unit,40)isol,ipert,m,mat(ipert,isol,:,iorder)
               m=m+1
            ENDDO
            WRITE(debug_unit,30)
     $           (TRIM(name),iqty,TRIM(name),iqty,iqty=1,2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_write_mat
c-----------------------------------------------------------------------
c     subprogram 11. sing1_vmat_diagnose.
c     diagnoses asymptotic behavior at the singular surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_vmat_diagnose(ising,mmat,vmat,amat,bmat)

      INTEGER, INTENT(IN) :: ising
      COMPLEX(r8), DIMENSION(:,:,:,0:) :: mmat,vmat
      COMPLEX(r8), DIMENSION(:,:,:,0:), OPTIONAL :: amat,bmat

      INTEGER :: r1
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointer and open file.
c-----------------------------------------------------------------------
      singp => sing(ising)
      r1=singp%r1(1)
      CALL ascii_open(debug_unit,"vmat.out","UNKNOWN")
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      WRITE(debug_unit,'(a,l1)')" sing1_flag = ",sing1_flag
      WRITE(debug_unit,'(4(a,i2))')
     $     " mpert = ",mpert,", msol = ",2*mpert,
     $     ", sing_order = ",sing_order,", singp%order = ",singp%order
      WRITE(debug_unit,'(1p,2(a,e10.3),a,2i3)')" q = ",singp%q,
     $     ", di = ",singp%di,", r2 = ",singp%r2
      WRITE(debug_unit,'(1p,3(a,e10.3)/)')
     $     " alpha = ",REAL(singp%alpha),
     $     ", alpha+0.5 = ",REAL(singp%alpha)+.5,
     $     ", alpha-0.5 = ",REAL(singp%alpha)-.5
c-----------------------------------------------------------------------
c     write matrices.
c-----------------------------------------------------------------------
      CALL sing1_write_mat(r1,mmat,"mmat",debug_unit)
      IF(PRESENT(amat))CALL sing1_write_mat(r1,amat,"amat",debug_unit)
      IF(PRESENT(bmat))CALL sing1_write_mat(r1,bmat,"bmat",debug_unit)
      CALL sing1_write_mat(r1,vmat,"vmat",debug_unit)
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_vmat_diagnose
c-----------------------------------------------------------------------
c     subprogram 12. sing1_ua_diagnose.
c     diagnoses asymptotic solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_ua_diagnose(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER, PARAMETER :: nz=10
      INTEGER :: isol,iqty,iz,mz,msol_old,r1
      REAL(r8) :: psi0,psi,psim,psip,z,zlog,dzlog,singfac,singlog
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,-1:1) :: u
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2,3) :: du
      COMPLEX(r8), DIMENSION(mpert,2*mpert) :: term1,term2,term3,matvec
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     open output files and set pointer.
c-----------------------------------------------------------------------
      CALL bin_open(bin1_unit,"ua.bin","UNKNOWN","REWIND","none")
      singp => sing(ising)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      psi0=sing(ising)%psifac
      mz=nz*(uad%zlogmax-uad%zlogmin)
      dzlog=(uad%zlogmax-uad%zlogmin)/mz
      msol_old=msol
      msol=2*mpert
      r1=sing(ising)%r1(1)
      WRITE(*,'(2(a,i2))')
     $     " sing_order = ",sing_order,", singp%order = ",singp%order
c-----------------------------------------------------------------------
c     start loops over isol, iqty, and iz.
c-----------------------------------------------------------------------
      zlog=uad%zlogmin
      DO iz=0,mz
c-----------------------------------------------------------------------
c     compute positions.
c-----------------------------------------------------------------------
         z=10**zlog
         IF(uad%neg)z=-z
         psi=psi0+z
         psim=psi0+z*(1-uad%eps)
         psip=psi0+z*(1+uad%eps)
c-----------------------------------------------------------------------
c     compute singfac and singlog.
c-----------------------------------------------------------------------
         CALL spline_eval(sq,psi,0)
         singfac=(sing(ising)%m-nn*sq%f(4))
         singlog=LOG10(ABS(singfac))
c-----------------------------------------------------------------------
c     compute power series solutions.
c-----------------------------------------------------------------------
         CALL sing1_get_ua(ising,psi,u(:,:,:,0))
         CALL sing1_get_ua(ising,psim,u(:,:,:,-1))
         CALL sing1_get_ua(ising,psip,u(:,:,:,1))
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         du(:,:,:,1)=(u(:,:,:,1)-u(:,:,:,-1))/(2*z*uad%eps)
         CALL sing1_get_dua(ising,psi,du(:,:,:,2))
c-----------------------------------------------------------------------
c     compute -(F u' + K u)' + (K^\dagger u' + G u)
c-----------------------------------------------------------------------
         CALL sing1_matvec(psi,u(:,:,:,0),du(:,:,:,2),
     $        term1,term2,term3,matvec)
c-----------------------------------------------------------------------
c     write graphical output.
c-----------------------------------------------------------------------
         WRITE(bin1_unit)REAL(zlog,4),REAL(singlog,4),
     $        ((sing1_log(u(r1,r2(isol),iqty,0)),
     $        sing1_log(du(r1,r2(isol),iqty,1)),
     $        sing1_log(du(r1,r2(isol),iqty,2)),iqty=1,2),
     $        sing1_log(term1(r1,r2(isol))),
     $        sing1_log(term2(r1,r2(isol))),
     $        sing1_log(term3(r1,r2(isol))),
     $        sing1_log(matvec(r1,r2(isol))),isol=1,2)
c-----------------------------------------------------------------------
c     finish loops over isol, iqty, and iz.
c-----------------------------------------------------------------------
         zlog=zlog+dzlog
      ENDDO
      WRITE(bin1_unit)
c-----------------------------------------------------------------------
c     close files and restore msol.
c-----------------------------------------------------------------------
      CALL bin_close(bin1_unit)
      msol=msol_old
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_ua_diagnose
c-----------------------------------------------------------------------
c     subprogram 13. sing1_vmat_d.
c     computes degenerate power series coefficients
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_vmat_d(ising,amat,bmat)

      INTEGER, INTENT(IN) :: ising
      COMPLEX(r8), DIMENSION(:,:,:,0:), INTENT(OUT) :: amat,bmat
      
      INTEGER :: k,k1,info,order,r1,ipert,jpert
      INTEGER, DIMENSION(2) :: ipiv
      INTEGER, DIMENSION(2*mpert) :: p
      INTEGER, DIMENSION(:), POINTER :: r2
      COMPLEX(r8), PARAMETER :: half=.5_r8,zero=0,one=1
      COMPLEX(r8), DIMENSION(2,2) :: mmat2,smat,sinv,mat2
      COMPLEX(r8), DIMENSION(mpert,2,2,0:2*sing(ising)%order) :: vmat
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert,0:2*sing(ising)%order) ::
     $     aflat,bflat,mflat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     select resonant vectors.
c-----------------------------------------------------------------------
      singp => sing(ising)
      r1=singp%r1(1)
      r2 => singp%r2
      order=singp%order
c-----------------------------------------------------------------------
c     flatten block-structured matrix mmat.
c-----------------------------------------------------------------------
      mflat(1:mpert,:,:)=singp%mmat(:,:,1,:)
      mflat(mpert+1:2*mpert,:,:)=singp%mmat(:,:,2,:)
c-----------------------------------------------------------------------
c     compute mmat2, smat, and sinv.
c-----------------------------------------------------------------------
      mmat2=mflat(r2,r2,0)
      smat=RESHAPE((/mmat2(1,2),-half-mmat2(1,1),
     $     mmat2(1,2),half-mmat2(1,1)/),(/2,2/))/mmat2(1,2)
      mat2=smat
      sinv=ident(2)
      CALL zgetrf(2,2,mat2,2,ipiv,info)
      CALL zgetrs('N',2,2,mat2,2,ipiv,sinv,2,info)
c-----------------------------------------------------------------------
c     compute aflat.
c-----------------------------------------------------------------------
      aflat=mflat
      DO k=0,2*order
         aflat(:,r2,k)=MATMUL(aflat(:,r2,k),smat)
         aflat(r2,:,k)=MATMUL(sinv,aflat(r2,:,k))
      ENDDO
      aflat(:,:,0)=0
      aflat(r2(1),r2(1),0)=-half
      aflat(r2(2),r2(2),0)=half
      singp%beta=aflat(r2(2),r2(1),2)
c-----------------------------------------------------------------------
c     compute bflat.
c-----------------------------------------------------------------------
      p=0
      p(r2)=(/-1,1/)
      bflat=0
      DO ipert=1,2*mpert
         DO jpert=1,2*mpert
            DO k=0,2*order
               k1=k+p(ipert)-p(jpert)
               IF(k1 < 0 .OR. k1 > 2*order)CYCLE
               bflat(ipert,jpert,k)=aflat(ipert,jpert,k1)
            ENDDO
         ENDDO
         bflat(ipert,ipert,0)=bflat(ipert,ipert,0)-half*p(ipert)
      ENDDO
c-----------------------------------------------------------------------
c     copy to block structured matrices.
c-----------------------------------------------------------------------
      amat(:,:,1,:)=aflat(1:mpert,:,:)
      amat(:,:,2,:)=aflat(mpert+1:2*mpert,:,:)
      bmat(:,:,1,:)=bflat(1:mpert,:,:)
      bmat(:,:,2,:)=bflat(mpert+1:2*mpert,:,:)
c-----------------------------------------------------------------------
c     compute zeroth-order resonant solutions.
c-----------------------------------------------------------------------
      vmat=0
      vmat(r1,1,1,0)=1
      vmat(r1,2,2,0)=1
c-----------------------------------------------------------------------
c     compute higher-order resonant solutions.
c-----------------------------------------------------------------------
      DO k=1,2*order
         CALL sing1_solve_d(ising,k,bmat,vmat)
      ENDDO
c-----------------------------------------------------------------------
c     copy resonant solutions to output.
c-----------------------------------------------------------------------
      singp%vmat(:,r2,:,:)=vmat
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_vmat_d
c-----------------------------------------------------------------------
c     subprogram 14. sing1_solve_d.
c     solves iteratively for the next order in the power series vmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_solve_d(ising,k,matrix,vector)

      INTEGER, INTENT(IN) :: ising,k
      COMPLEX(r8), DIMENSION(:,:,:,0:), INTENT(IN) :: matrix
      COMPLEX(r8), DIMENSION(:,:,:,0:), INTENT(INOUT) :: vector

      INTEGER :: l,info
      INTEGER, DIMENSION(2*mpert) :: ipiv
      REAL(r8), PARAMETER :: half=.5
      COMPLEX(r8), DIMENSION(2*mpert,2) :: vec
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: cmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     compute rhs.
c-----------------------------------------------------------------------
      DO l=1,k
         vector(:,:,:,k)=vector(:,:,:,k)
     $        -sing1_matmul(matrix(:,:,:,l),vector(:,:,:,k-l))
      ENDDO
c-----------------------------------------------------------------------
c     compute and factor matrix.
c-----------------------------------------------------------------------
      singp => sing(ising)
      cmat=-ident(2*mpert)*k*half
      cmat(singp%r2(2),singp%r2(1))=singp%beta
      CALL zgetrf(2*mpert,2*mpert,cmat,2*mpert,ipiv,info)
c-----------------------------------------------------------------------
c     solve.
c-----------------------------------------------------------------------
      vec(1:mpert,:)=vector(:,:,1,k)
      vec(mpert+1:2*mpert,:)=vector(:,:,2,k)
      CALL zgetrs('N',2*mpert,1,cmat,2*mpert,ipiv,vec(:,2),2*mpert,info)
      vec(:,1)=vec(:,1)+singp%beta*vec(:,2)
      CALL zgetrs('N',2*mpert,1,cmat,2*mpert,ipiv,vec(:,1),2*mpert,info)
c-----------------------------------------------------------------------
c     transfer solutions to output.
c-----------------------------------------------------------------------
      vector(:,:,1,k)=vec(1:mpert,:)
      vector(:,:,2,k)=vec(mpert+1:2*mpert,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_solve_d
c-----------------------------------------------------------------------
c     subprogram 15. sing1_get_ua_d.
c     computes degenerate power series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_ua_d(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: ua

      INTEGER :: iorder,r1
      INTEGER, DIMENSION(:), POINTER :: r2
      COMPLEX(r8), PARAMETER :: half=.5_r8,zero=0,one=1
      COMPLEX(r8) :: z,sqrtfac
      COMPLEX(r8), DIMENSION(2,2) :: b0exp,smat,mmat2,tmat,shear
      COMPLEX(r8), DIMENSION(2*mpert,2) :: ur
      COMPLEX(r8), DIMENSION(2*mpert,2,0:2*sing(ising)%order) :: vflat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     select resonant vectors.
c-----------------------------------------------------------------------
      singp => sing(ising)
      r1=singp%r1(1)
      r2 => singp%r2
c-----------------------------------------------------------------------
c     flatten block-structured matrix mmat.
c-----------------------------------------------------------------------
      vflat(1:mpert,:,:)=singp%vmat(:,r2,1,:)
      vflat(mpert+1:2*mpert,:,:)=singp%mmat(:,r2,2,:)
c-----------------------------------------------------------------------
c     compute power series by horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      ur=vflat(:,:,2*singp%order)
      DO iorder=2*singp%order-1,0,-1
         ur=ur*sqrtfac+vflat(:,:,iorder)
      ENDDO
c-----------------------------------------------------------------------
c     compute transformation matrices.
c-----------------------------------------------------------------------
c      b0exp=RESHAPE((/one,singp%beta*LOG(z),zero,one/),(/2,2/))
      b0exp=RESHAPE((/one,singp%beta*LOG(ABS(z)),zero,one/),(/2,2/))
      tmat=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
      mmat2=TRANSPOSE(singp%mmat(r1,r2,:,0))
      smat=RESHAPE((/mmat2(1,2),-half-mmat2(1,1),
     $     mmat2(1,2),half-mmat2(1,1)/),(/2,2/))/mmat2(1,2)
      shear=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
c-----------------------------------------------------------------------
c     apply transformation matrices.
c-----------------------------------------------------------------------
      ur=MATMUL(ur,b0exp)
      ur(r2,:)=MATMUL(tmat,ur(r2,:))
      ur(r2,:)=MATMUL(smat,ur(r2,:))
      ur(r2,:)=MATMUL(shear,ur(r2,:))
c-----------------------------------------------------------------------
c     transfer to return value.
c-----------------------------------------------------------------------
      ua(:,r2,1)=ur(1:mpert,:)
      ua(:,r2,2)=ur(mpert+1:2*mpert,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_ua_d
c-----------------------------------------------------------------------
c     subprogram 16. sing1_get_dua_d.
c     computes degenerate power series resonant solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_dua_d(ising,psifac,dua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: dua

      INTEGER :: r1,iorder
      INTEGER, DIMENSION(:), POINTER :: r2
      COMPLEX(r8), PARAMETER :: half=.5_r8,zero=0,one=1
      COMPLEX(r8) :: z,sqrtfac
      COMPLEX(r8), DIMENSION(2,2) :: b0mat,b0exp,smat,mmat2,
     $     tmat,zdb0exp,zdtmat,shear,zdshear
      COMPLEX(r8), DIMENSION(mpert,2,2,0:2*sing(ising)%order) :: vmat
      COMPLEX(r8), DIMENSION(2*mpert,2) :: xa
      COMPLEX(r8), DIMENSION(2*mpert,2,4) :: zdxai
      COMPLEX(r8), DIMENSION(2*mpert,2,0:2*sing(ising)%order) :: vflat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     select resonant vectors.
c-----------------------------------------------------------------------
      singp => sing(ising)
      r1=singp%r1(1)
      r2 => singp%r2
      vmat=singp%vmat(:,r2,:,:)
c-----------------------------------------------------------------------
c     flatten block-structured matrix mmat.
c-----------------------------------------------------------------------
      vflat(1:mpert,:,:)=singp%vmat(:,r2,1,:)
      vflat(mpert+1:2*mpert,:,:)=singp%mmat(:,r2,2,:)
c-----------------------------------------------------------------------
c     compute power series by horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      xa=vflat(:,:,2*singp%order)
      zdxai=0
      zdxai(:,:,1)=vflat(:,:,2*singp%order)*singp%order
      DO iorder=2*singp%order-1,0,-1
         xa=xa*sqrtfac+vflat(:,:,iorder)
         zdxai(:,:,1)=zdxai(:,:,1)*sqrtfac+vflat(:,:,iorder)*iorder/2
      ENDDO
c-----------------------------------------------------------------------
c     compute transformation matrices.
c-----------------------------------------------------------------------
      b0mat=RESHAPE((/zero,singp%beta,zero,zero/),(/2,2/))
c      b0exp=RESHAPE((/one,singp%beta*LOG(z),zero,one/),(/2,2/))
      b0exp=RESHAPE((/one,singp%beta*LOG(ABS(z)),zero,one/),(/2,2/))
      tmat=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
      mmat2=TRANSPOSE(singp%mmat(r1,r2,:,0))
      smat=RESHAPE((/mmat2(1,2),-half-mmat2(1,1),
     $     mmat2(1,2),half-mmat2(1,1)/),(/2,2/))/mmat2(1,2)
      shear=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
c-----------------------------------------------------------------------
c     compute derivatives of transposed transformation matrices.
c-----------------------------------------------------------------------
      zdb0exp=MATMUL(b0exp,b0mat)
      zdtmat=RESHAPE((/-half/sqrtfac,zero,zero,half*sqrtfac/),(/2,2/))
      zdshear=RESHAPE((/-half/sqrtfac,zero,zero,half*sqrtfac/),(/2,2/))
c-----------------------------------------------------------------------
c     apply b0exp.
c-----------------------------------------------------------------------
      xa=MATMUL(xa,b0exp)
      zdxai(:,:,1)=MATMUL(zdxai(:,:,1),b0exp)
      zdxai(:,:,2)=MATMUL(xa(:,:),b0mat)
c-----------------------------------------------------------------------
c     compute zdxa1.
c-----------------------------------------------------------------------
      zdxai(r2,:,1)=MATMUL(tmat,zdxai(r2,:,1))
      zdxai(r2,:,1)=MATMUL(smat,zdxai(r2,:,1))
      zdxai(r2,:,1)=MATMUL(shear,zdxai(r2,:,1))
c-----------------------------------------------------------------------
c     compute zdxa2.
c-----------------------------------------------------------------------
      zdxai(r2,:,2)=MATMUL(tmat,zdxai(r2,:,2))
      zdxai(r2,:,2)=MATMUL(smat,zdxai(r2,:,2))
      zdxai(r2,:,2)=MATMUL(shear,zdxai(r2,:,2))
c-----------------------------------------------------------------------
c     compute zdxa3.
c-----------------------------------------------------------------------
      zdxai(r2,:,3)=MATMUL(zdtmat,xa(r2,:))
      zdxai(r2,:,3)=MATMUL(smat,zdxai(r2,:,3))
      zdxai(r2,:,3)=MATMUL(shear,zdxai(r2,:,3))
c-----------------------------------------------------------------------
c     compute zdxa4.
c-----------------------------------------------------------------------
      zdxai(r2,:,4)=MATMUL(tmat,xa(r2,:))
      zdxai(r2,:,4)=MATMUL(smat,zdxai(r2,:,4))
      zdxai(r2,:,4)=MATMUL(zdshear,zdxai(r2,:,4))
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      dua(:,r2,1)=SUM(zdxai(1:mpert,:,:),3)/z
      dua(:,r2,2)=SUM(zdxai(mpert+1:2*mpert,:,:),3)/z
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_dua_d
c-----------------------------------------------------------------------
c     subprogram 17. sing1_get_ua_d_cut.
c     computes degenerate power series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_ua_d_cut(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: ua

      INTEGER :: r1
      INTEGER, DIMENSION(:), POINTER :: r2
      COMPLEX(r8), PARAMETER :: half=.5_r8,zero=0,one=1
      COMPLEX(r8) :: z,sqrtfac
      COMPLEX(r8), DIMENSION(2,2) :: b0exp,smat,mmat2,tmat,shear
      COMPLEX(r8), DIMENSION(2*mpert,2) :: ur
      COMPLEX(r8), DIMENSION(2*mpert,2,0:2*sing(ising)%order) :: vflat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     select resonant vectors.
c-----------------------------------------------------------------------
      singp => sing(ising)
      r1=singp%r1(1)
      r2 => singp%r2
c-----------------------------------------------------------------------
c     flatten block-structured matrix mmat.
c-----------------------------------------------------------------------
      vflat(1:mpert,:,:)=singp%vmat(:,r2,1,:)
      vflat(mpert+1:2*mpert,:,:)=singp%mmat(:,r2,2,:)
c-----------------------------------------------------------------------
c     compute power series by horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      ur=vflat(:,:,0)
c      DO iorder=2*singp%order-1,1,-1
c         ur=ur*sqrtfac+vflat(:,:,iorder)
c      ENDDO
c-----------------------------------------------------------------------
c     compute transformation matrices.
c-----------------------------------------------------------------------
c      b0exp=RESHAPE((/one,singp%beta*LOG(z),zero,one/),(/2,2/))
      b0exp=RESHAPE((/one,singp%beta*LOG(ABS(z)),zero,one/),(/2,2/))
      tmat=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
      mmat2=TRANSPOSE(singp%mmat(r1,r2,:,0))
      smat=RESHAPE((/mmat2(1,2),-half-mmat2(1,1),
     $     mmat2(1,2),half-mmat2(1,1)/),(/2,2/))/mmat2(1,2)
      shear=RESHAPE((/one/sqrtfac,zero,zero,sqrtfac/),(/2,2/))
c-----------------------------------------------------------------------
c     apply transformation matrices.
c-----------------------------------------------------------------------
      ur=MATMUL(ur,b0exp)
      ur(r2,:)=MATMUL(tmat,ur(r2,:))
      ur(r2,:)=MATMUL(smat,ur(r2,:))
      ur(r2,:)=MATMUL(shear,ur(r2,:))
c-----------------------------------------------------------------------
c     transfer to return value.
c-----------------------------------------------------------------------
      ua(:,r2,1)=ur(1:mpert,:)
      ua(:,r2,2)=ur(mpert+1:2*mpert,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_ua_d_cut
c-----------------------------------------------------------------------
c     subprogram 18. sing1_get_ua_cut.
c     computes asymptotic series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_get_ua_cut(ising,psifac,ua)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: ua

      INTEGER :: isol
      INTEGER, DIMENSION(:), POINTER :: r1,r2
      COMPLEX(r8) :: pfac,z,sqrtfac,phase
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: vmat
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      singp => sing(ising)
      vmat => singp%vmat
      r1 => singp%r1
      r2 => singp%r2
c-----------------------------------------------------------------------
c     compute power series by Horner's method.
c-----------------------------------------------------------------------
      z=psifac-singp%psifac
      sqrtfac=SQRT(z)
      ua=vmat(:,:,:,0)
c      DO iorder=2*singp%order-1,1,-1
c         ua=ua*sqrtfac+vmat(:,:,:,iorder)
c      ENDDO
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      pfac=z**singp%alpha
      ua(r1,:,1)=ua(r1,:,1)/sqrtfac
      ua(r1,:,2)=ua(r1,:,2)*sqrtfac
      ua(:,r2(1),:)=ua(:,r2(1),:)/pfac
      ua(:,r2(2),:)=ua(:,r2(2),:)*pfac
c-----------------------------------------------------------------------
c     compute degenerate solutions.
c-----------------------------------------------------------------------
      IF(singp%degen)CALL sing1_get_ua_d_cut(ising,psifac,ua)
c-----------------------------------------------------------------------
c     adjust phase.
c-----------------------------------------------------------------------
      IF(psifac < singp%psifac .AND. uad%phase)THEN
         phase=EXP(ifac*pi*(.5+singp%alpha))
         DO isol=1,2
            ua(:,r2(isol),:)=ua(:,r2(isol),:)*phase
            phase=EXP(ifac*pi*(.5-singp%alpha))
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_get_ua_cut
c-----------------------------------------------------------------------
c     subprogram 19. sing1_delta.
c     computes delta.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_delta(ising,psifac,delta)
      
      INTEGER, INTENT(IN) :: ising
      REAL(r8), INTENT(IN) :: psifac
      REAL(r8), DIMENSION(2), INTENT(OUT) :: delta

      INTEGER :: j,ipert0,neq=1
      INTEGER :: msol_save
      INTEGER :: range(2)
      REAL(r8), DIMENSION(0:2) :: norm
      COMPLEX(r8), DIMENSION(mpert,2,2,0:2) ::matvec
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua,dua
      COMPLEX(r8), DIMENSION(mpert,2,2) :: du
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     choose singular solutions.
c-----------------------------------------------------------------------
      singp => sing(ising)
      ipert0=singp%m-mlow+1
      range=(/ipert0,ipert0+mpert/)
      msol_save=msol
      msol=2
c-----------------------------------------------------------------------
c     compute ua, dua, and du.
c-----------------------------------------------------------------------
      CALL sing1_get_ua(ising,psifac,ua)
      CALL sing1_get_dua(ising,psifac,dua)
      CALL sing1_der(neq,psifac,ua(:,range,:),du)
c-----------------------------------------------------------------------
c     compute matvec.
c-----------------------------------------------------------------------
      matvec(:,:,:,1)=dua(:,range,:)
      matvec(:,:,:,2)=-du
      matvec(:,:,:,0)=SUM(matvec(:,:,:,1:2),4)
c-----------------------------------------------------------------------
c     compute delta.
c-----------------------------------------------------------------------
      DO j=1,2
         norm=MAXVAL(MAXVAL(ABS(matvec(:,j,:,:)),1),1)
         delta(j)=norm(0)/MAXVAL(norm(1:2))
      ENDDO
      msol=msol_save
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_delta 
c-----------------------------------------------------------------------
c     subprogram 20. sing1_kxscan.
c     tests convergence of differential equation.
c----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_kxscan

      INTEGER :: ua_bin_unit=51
      INTEGER :: klow,khigh,kmax,dkmax,ixlog,mxlog,nxlog,ik,nk,
     $     ising,order_save
      REAL(r8) :: xlogmin,xlogmax,xmax,xlog,dxlog,x,dxfac,psifac
      REAL(r8), DIMENSION(2) :: delta
      TYPE(sing_type), POINTER :: singp

      NAMELIST/kxscan_list/ising,klow,khigh,dkmax,xlogmin,xlogmax,mxlog
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="dcon.in",STATUS="OLD")
      READ(in_unit,NML=kxscan_list)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     set pointer and save order.
c-----------------------------------------------------------------------
      singp => sing(ising)
      order_save=singp%order
c-----------------------------------------------------------------------
c     prepare for scans.
c-----------------------------------------------------------------------
      dxlog=-one/mxlog
      nxlog=(xlogmax-xlogmin)*mxlog
      dxfac=10**dxlog
      xmax=10**xlogmax
      nk=(khigh-klow)/dkmax
      kmax=klow
c-----------------------------------------------------------------------
c     scans over kmax and x.
c-----------------------------------------------------------------------
      OPEN(UNIT=ua_bin_unit,FILE="kxscan.bin",STATUS="REPLACE",
     $     FORM="UNFORMATTED")
      DO ik=0,nk
         sing_order=kmax
         CALL sing1_vmat(ising)
         x=xmax
         xlog=xlogmax
         DO ixlog=0,nxlog
            psifac=singp%psifac+x
            CALL sing1_delta(ising,psifac,delta)
            WRITE(ua_bin_unit)REAL(xlog,4),REAL(LOG10(delta),4)
            xlog=xlog+dxlog
            x=x*dxfac
         ENDDO
         WRITE(ua_bin_unit)
         DEALLOCATE(singp%vmat,singp%mmat)
         kmax=kmax+dkmax
      ENDDO
      CLOSE(UNIT=ua_bin_unit)
c-----------------------------------------------------------------------
c     restore order.
c-----------------------------------------------------------------------
      singp%order=order_save
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_kxscan
c-----------------------------------------------------------------------
c     subprogram 21. sing1_xmin.
c     computes xmin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sing1_xmin(ising,eps,xmin)

      INTEGER, INTENT(IN) :: ising
      REAL(r8), DIMENSION(0:2), INTENT(IN) :: eps
      REAL(r8), DIMENSION(0:2), INTENT(OUT) :: xmin

      LOGICAL, DIMENSION(0:1) :: set
      CHARACTER(64) :: message
      REAL(r8), PARAMETER :: xlogmax=-1,dxlog=-.01
      INTEGER :: ixlog,nxlog=5000,i
      REAL(r8) :: xlog,x,dxfac,psifac
      REAL(r8), DIMENSION(2) :: delta
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     prepare for scan.
c-----------------------------------------------------------------------
      set=.TRUE.
      dxfac=10**dxlog
      xlog=xlogmax
      x=10**xlog
      ixlog=0
      singp => sing(ising)
c-----------------------------------------------------------------------
c     start loops over x.
c-----------------------------------------------------------------------
      DO
         psifac=singp%psifac+x
         CALL sing1_delta(ising,psifac,delta)
         DO i=0,1
            IF(MAXVAL(delta) < eps(i) .AND. set(i))THEN
               xmin(i)=x
               set(i)=.FALSE.
            ENDIF
         ENDDO
         IF(MAXVAL(delta) < eps(2))THEN
            xmin(2)=x
            EXIT
         ENDIF
c-----------------------------------------------------------------------
c     abort if ixlog = nxlog.
c-----------------------------------------------------------------------
         IF(ixlog == nxlog)THEN
            WRITE(message,'(a,g0,a,es10.3)')
     $           "sing1_xmin: abort with ixlog = nxlog = ",nxlog,
     $           ", xlog = ",xlog
            CALL program_stop(message)
         ENDIF
c-----------------------------------------------------------------------
c     finish loop over x.
c-----------------------------------------------------------------------
         ixlog=ixlog+1
         xlog=xlog+dxlog
         x=x*dxfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sing1_xmin
      END MODULE sing1_mod
