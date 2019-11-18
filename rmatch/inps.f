c-----------------------------------------------------------------------
c     file inps.f.
c     solves for inner region power series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inps_mod.
c     1. inps_init.
c     2. inps_dealloc.
c     3. inps_cmatwrite.
c     4. inps_tjmat.
c     5. inps_lyap_solve.
c     6. inps_split.
c     7. inps_coefs.
c     8. inps_ua.
c     9. inps_horner.
c     10. inps_xmax.
c     11. inps_delta.
c-----------------------------------------------------------------------
c     subprogram 0. inps_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inps_mod
      USE local_mod
      IMPLICIT NONE

      TYPE :: resist_type_inps
      REAL(r8) :: e,f,h,di,g,k
      REAL(r8), DIMENSION(2) :: r
      COMPLEX(r8) :: q
      END TYPE resist_type_inps

      LOGICAL :: extract=.FALSE.
      INTEGER, PRIVATE :: kmax=8
      INTEGER, DIMENSION(2), PARAMETER, PRIVATE :: r1=(/1,2/)
      INTEGER, DIMENSION(4), PARAMETER, PRIVATE :: r2=(/3,4,5,6/)
      COMPLEX(r8), PRIVATE :: lambda
      REAL(r8), DIMENSION(2,2), PRIVATE :: rmat
      COMPLEX(r8), DIMENSION(6,2), PRIVATE :: pp,dpp
      COMPLEX(r8), DIMENSION(6,6), PRIVATE :: tmat,tinv
      COMPLEX(r8), DIMENSION(6,6,0:2), PRIVATE :: amat,jmat
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: 
     $     k2mat,k6mat,bmat,pmat,qmat,cmat,dmat,ymat,emat,zmat
      TYPE(resist_type_inps) :: rt

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. inps_init.
c     initializes amat, tmat, tinv, and jmat matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_init(rt_in,diagnose)

      LOGICAL, INTENT(IN) :: diagnose
      TYPE(resist_type_inps), INTENT(IN) :: rt_in

      REAL(r8) :: sqrtfac
c-----------------------------------------------------------------------
c     copy input values.
c-----------------------------------------------------------------------
      rt%e=rt_in%e 
      rt%f=rt_in%f
      rt%h=rt_in%h
      rt%g=rt_in%g
      rt%k=rt_in%k 
      rt%q=rt_in%q 
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      rt%di=rt%e+rt%f+rt%h-.25
      sqrtfac=SQRT(-rt%di)
      rt%r=(/sqrtfac,-sqrtfac/)+1.5
c-----------------------------------------------------------------------
c     create arrays.
c-----------------------------------------------------------------------
      CALL inps_tjmat(diagnose)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_init
c-----------------------------------------------------------------------
c     subprogram 2. inps_dealloc.
c     deallocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_dealloc
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(bmat,pmat,k6mat,qmat,cmat,ymat,k2mat,dmat,emat,zmat)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. inps_cmatwrite.
c     writes complex matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_cmatwrite(mat,name,unit)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: mat
      CHARACTER(*), INTENT(IN) :: name 
      INTEGER, INTENT(IN) :: unit

      INTEGER :: m,n,i,j
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",<n>(5x,"re ",i1,7x,"im ",i1,2x)/)
 20   FORMAT(i6,<2*n>es11.3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      m=SIZE(mat,1)
      n=SIZE(mat,2)
      WRITE(unit,'(1x,a,":")')TRIM(name)
      WRITE(unit,10)(j,j,j=1,n)
      WRITE(unit,20)(i,mat(i,:),i=1,m)
      WRITE(UNIT,'()')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_cmatwrite
c-----------------------------------------------------------------------
c     subprogram 4. inps_tjmat.
c     computes coefficient matrices; inps.pdf pages 1 and 2.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_tjmat(diagnose)

      LOGICAL, INTENT(IN) :: diagnose

      INTEGER :: i
      COMPLEX(r8) :: e,f,h,g,k,q,q2
      COMPLEX(r8), DIMENSION(2) :: r
      COMPLEX(r8), DIMENSION(2,2) :: rmat
      COMPLEX(r8), PARAMETER :: zero=0,one=1,two=2,half=one/2
c-----------------------------------------------------------------------
c     compute scalars.
c-----------------------------------------------------------------------
      q=rt%q
      e=rt%e
      f=rt%f
      h=rt%h
      g=rt%g
      k=rt%k
      r=r
      lambda=one/SQRT(rt%q)
      q2=rt%q**2
      rmat=RESHAPE((/r(1),zero,zero,r(2)/),(/2,2/))
c-----------------------------------------------------------------------
c     compute tmat.
c-----------------------------------------------------------------------
      tmat=RESHAPE(
     $     (/one,zero,zero,zero,zero,zero,
     $     zero,zero,zero,one,zero,zero,
     $     h*q,zero,-one/lambda,-h/lambda,zero,one,
     $     q2/lambda,-one/lambda,zero,-q2,one,zero,
     $     h*q,zero,one/lambda,h/lambda,zero,one,
     $     -q2/lambda,one/lambda,zero,-q2,one,zero/),(/6,6/))
c-----------------------------------------------------------------------
c     compute tinv.
c-----------------------------------------------------------------------
      tinv=RESHAPE(
     $     (/one,zero,zero,zero,zero,zero,
     $     q2,zero,zero,-half*lambda,zero,half*lambda,
     $     zero,-h,-half*lambda,zero,half*lambda,zero,
     $     zero,one,zero,zero,zero,zero,
     $     zero,q2,zero,half,zero,half,
     $     -h*q,zero,half,zero,half,zero/),(/6,6/))
c-----------------------------------------------------------------------
c     compute amat.
c-----------------------------------------------------------------------
      amat=0
      amat(4:5,2,0)=(/-q,one/q/)
      amat(6,3,0)=one/q
      amat(1,4,0)=one
      amat(2,5,0)=one
      amat(3:4,6,0)=(/one,h/)
      amat(4:6,1,1)=(/q,-one/q,-one/q/)
      amat(6,2,1)=-(g-k*e)*q
      amat(5:6,3,1)=(/-(e+f)/q2,(g+k*f)*q/)
      amat(4:6,4,1)=(/one,-h/q2,h*k*q/)
      amat(5,5,1)=-1
      amat(6,6,1)=-1
      amat(4:6,1,2)=(/-two,h/q2,-h*k*q/)
c-----------------------------------------------------------------------
c     compute jmat.
c-----------------------------------------------------------------------
      DO i=0,2
         jmat(:,:,i)=MATMUL(tinv,MATMUL(amat(:,:,i),tmat))
      ENDDO
      IF(.NOT. diagnose)RETURN
c-----------------------------------------------------------------------
c     write matrices.
c-----------------------------------------------------------------------
      OPEN(UNIT=jmat_unit,FILE="jmat.out",STATUS="REPLACE")
      WRITE(jmat_unit,
     $     '(4(a,es10.3)/2(a,es10.3),a,2es11.3/a,es10.3,a,2es10.3/)')
     $     " e = ",rt%e,", f = ",rt%f,", h = ",rt%h,", di = ",rt%di,
     $     " g = ",rt%g,", k = ",rt%k,", r = ",rt%r,
     $     " r1 - r2 = ",rt%r(1)-rt%r(2),", q = ",rt%q
      CALL inps_cmatwrite(tmat,"tmat",jmat_unit)
      CALL inps_cmatwrite(tinv,"tinv",jmat_unit)
      CALL inps_cmatwrite(MATMUL(tmat,tinv),"tmat.tinv",jmat_unit)
      CALL inps_cmatwrite(amat(:,:,0),"a0mat",jmat_unit)
      CALL inps_cmatwrite(amat(:,:,1),"a1mat",jmat_unit)
      CALL inps_cmatwrite(amat(:,:,2),"a2mat",jmat_unit)
      CALL inps_cmatwrite(jmat(:,:,0),"j0mat",jmat_unit)
      CALL inps_cmatwrite(jmat(:,:,1),"j1mat",jmat_unit)
      CALL inps_cmatwrite(jmat(:,:,2),"j2mat",jmat_unit)
      CLOSE(UNIT=jmat_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_tjmat
c-----------------------------------------------------------------------
c     subprogram 5. inps_lyap_solve.
c     solves Lyapunov equation; inps.pdf Eq. (22).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_lyap_solve(kmat,bmat,pmat)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: kmat
      COMPLEX(r8), DIMENSION(:,:), INTENT(OUT) :: pmat,bmat
c-----------------------------------------------------------------------
c     solve for bmat.
c-----------------------------------------------------------------------
      bmat=0
      bmat(r1,r1)=kmat(r1,r1)
      bmat(r2,r2)=kmat(r2,r2)
c-----------------------------------------------------------------------
c     solve for pmat(r1,r2), new way.
c-----------------------------------------------------------------------
      pmat=0
      pmat(2,3:4)=-kmat(2,3:4)/lambda
      pmat(2,5:6)=kmat(2,5:6)/lambda
      pmat(1,3:4)=-(kmat(1,3:4)+pmat(2,3:4))/lambda
      pmat(1,5:6)=(kmat(1,5:6)+pmat(2,5:6))/lambda
c-----------------------------------------------------------------------
c     solve for pmat(r2,r1), new way.
c-----------------------------------------------------------------------
      pmat(3:4,1)=kmat(3:4,1)/lambda
      pmat(5:6,1)=-kmat(5:6,1)/lambda
      pmat(3:4,2)=(kmat(3:4,2)-pmat(3:4,1))/lambda
      pmat(5:6,2)=-(kmat(5:6,2)-pmat(5:6,1))/lambda
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_lyap_solve
c-----------------------------------------------------------------------
c     subprogram 6. inps_split.
c     splitting transformation, inps.pdf page 3.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_split(kmax_in,diagnose)

      INTEGER, INTENT(IN) :: kmax_in
      LOGICAL, INTENT(IN) :: diagnose

      INTEGER :: k,l
      CHARACTER(4) :: order
      REAL(r8), DIMENSION(kmax+2) :: error,norm
      COMPLEX(r8), DIMENSION(6,6,0:kmax+2) :: errmat
c-----------------------------------------------------------------------
c     compute coefficient matrices.
c-----------------------------------------------------------------------
      CALL inps_tjmat(diagnose)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      kmax=kmax_in
      ALLOCATE(bmat(6,6,0:kmax+2),pmat(6,6,0:kmax+2),
     $     k6mat(6,6,0:kmax+2),qmat(2,2,0:kmax+2),cmat(2,2,0:kmax+2),
     $     ymat(2,2,0:kmax),k2mat(2,2,kmax+2),dmat(2,2,0:kmax),
     $     emat(2,2,kmax),zmat(2,2,0:kmax))
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      bmat=0
      pmat=0
      k6mat=0
      bmat(:,:,0)=jmat(:,:,0)
      pmat(:,:,0)=ident(6)
      k6mat(:,:,0)=jmat(:,:,0)
c-----------------------------------------------------------------------
c     compute higher order matrices.
c-----------------------------------------------------------------------
      DO k=1,kmax+2
         IF(k <= 2)k6mat(:,:,k)=jmat(:,:,k)
         IF(k > 1)k6mat(:,:,k)=k6mat(:,:,k)+2*(k-1)*pmat(:,:,k-1)
         DO l=1,k-1
            IF(k-l <= 2)k6mat(:,:,k)=k6mat(:,:,k)
     $           +MATMUL(jmat(:,:,k-l),pmat(:,:,l))
            k6mat(:,:,k)=k6mat(:,:,k)-MATMUL(pmat(:,:,l),bmat(:,:,k-l))
         ENDDO
         CALL inps_lyap_solve(k6mat(:,:,k),bmat(:,:,k),pmat(:,:,k))
      ENDDO
c-----------------------------------------------------------------------
c     compute coefs.
c-----------------------------------------------------------------------
      CALL inps_coefs(diagnose)
c-----------------------------------------------------------------------
c     compute errors.
c-----------------------------------------------------------------------
      IF(.NOT. diagnose)RETURN
      norm=0
      errmat=0
      DO k=1,kmax+2
         norm(k)=MAXVAL(ABS(k6mat(:,:,k)))
         errmat(:,:,k)=MATMUL(jmat(:,:,0),pmat(:,:,k))
     $        -MATMUL(pmat(:,:,k),jmat(:,:,0))+k6mat(:,:,k)-bmat(:,:,k)
         errmat(:,:,k)=errmat(:,:,k)/norm(k)
         error(k)=MAXVAL(ABS(errmat(:,:,k)))
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      OPEN(UNIT=split_unit,FILE="split.out",STATUS="REPLACE")
      WRITE(split_unit,'(a,i2/)')" kmax = ",kmax
      DO k=1,kmax+2
         WRITE(split_unit,'(a,g0/)')" k = ",k
         WRITE(order,'(a,g0,a)')"(",k,")"
         IF(k <= 2)
     $        CALL inps_cmatwrite(jmat(:,:,k),"jmat"//order,split_unit)
         CALL inps_cmatwrite(k6mat(:,:,k),"k6mat"//order,split_unit)
         CALL inps_cmatwrite(bmat(:,:,k),"bmat"//order,split_unit) 
         CALL inps_cmatwrite(pmat(:,:,k),"pmat"//order,split_unit)
         CALL inps_cmatwrite(errmat(:,:,k),"errmat"//order,split_unit) 
         WRITE(order,'(a,i2.2,a)')"(",k,")"
         WRITE(split_unit,'(a,es10.3/)')
     $        " error"//order//" = ",error(k)
      ENDDO
      CLOSE(UNIT=split_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_split
c-----------------------------------------------------------------------
c     subprogram 7. inps_coefs.
c     computes qmat, cmat, dmat, emat, ymat, inps.pdf pages 6-8.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_coefs(diagnose)

      LOGICAL, INTENT(IN) :: diagnose

      CHARACTER(4) :: order
      INTEGER :: i,j,k,l
      REAL(r8) :: yerror
      REAL(r8), DIMENSION(kmax+2) :: qcerr
      COMPLEX(r8), PARAMETER :: zero=0,one=1,two=2,half=one/2,three=3
      COMPLEX(r8), DIMENSION(2) :: r
      COMPLEX(r8), DIMENSION(2,2) :: y0mat,y0inv,idmat
      COMPLEX(r8), DIMENSION(2,2) :: yerrmat
      COMPLEX(r8), DIMENSION(2,2,kmax+2) :: qcemat
c-----------------------------------------------------------------------
c     zeroth order terms.
c-----------------------------------------------------------------------
      qmat=0
      cmat=0
      qmat(:,:,0)=ident(2)
      cmat(:,:,0)=bmat(r1,r1,0)
c-----------------------------------------------------------------------
c     higher order terms.
c-----------------------------------------------------------------------
      DO k=1,kmax+2
         k2mat(:,:,k)=bmat(r1,r1,k)+2*(k-1)*qmat(:,:,k-1)
         DO l=1,k-1
            k2mat(:,:,k)=k2mat(:,:,k)
     $           +MATMUL(bmat(r1,r1,k-l),qmat(:,:,l))
     $           -MATMUL(qmat(:,:,l),cmat(:,:,k-l))
         ENDDO
         qmat(:,:,k)=RESHAPE((/zero,-k2mat(1,1,k),zero,-k2mat(1,2,k)/),
     $        (/2,2/))
         cmat(:,:,k)=RESHAPE(
     $        (/zero,k2mat(2,1,k),zero,k2mat(1,1,k)+k2mat(2,2,k)/),
     $        (/2,2/))
      ENDDO
c-----------------------------------------------------------------------
c     compute dmat.
c-----------------------------------------------------------------------
      dmat=0
      dmat(:,:,0)=RESHAPE((/zero,cmat(2,1,2),one,three/),(/2,2/))
      DO k=1,kmax
         dmat(:,:,k)=RESHAPE((/zero,cmat(2,1,k+2),zero,cmat(2,2,k+1)/),
     $        (/2,2/))
      ENDDO
c-----------------------------------------------------------------------
c     lowest-order solution and inverse.
c-----------------------------------------------------------------------
      r=rt%r
      y0mat=RESHAPE((/one,r(1),one,r(2)/),(/2,2/))
      y0inv=RESHAPE((/-r(2),r(1),one,-one/),(/2,2/))
     $     /(r(1)-r(2))
c-----------------------------------------------------------------------
c     compute emat, zmat, and ymat.
c-----------------------------------------------------------------------
      zmat=0
      ymat=0
      zmat(:,:,0)=ident(2)
      ymat(:,:,0)=y0mat
      DO k=1,kmax
         emat(:,:,k)=MATMUL(y0inv,MATMUL(dmat(:,:,k),y0mat))
         DO l=1,k
            zmat(:,:,k)=zmat(:,:,k)+MATMUL(emat(:,:,l),zmat(:,:,k-l))
         ENDDO
         DO i=1,2
            DO j=1,2
               zmat(i,j,k)=zmat(i,j,k)/(rt%r(j)-rt%r(i)-2*k)
            ENDDO
         ENDDO
         ymat(:,:,k)=MATMUL(y0mat,zmat(:,:,k))
      ENDDO
      IF(.NOT. diagnose)RETURN
c-----------------------------------------------------------------------
c     compute qcemat and qcerr.
c-----------------------------------------------------------------------
      DO k=1,kmax+2
         qcemat(:,:,k)=MATMUL(bmat(r1,r1,0),qmat(:,:,k))
     $        -MATMUL(qmat(:,:,k),bmat(r1,r1,0))
     $        +k2mat(:,:,k)-cmat(:,:,k)
         qcerr=MAXVAL(ABS(qcemat))
      ENDDO
c-----------------------------------------------------------------------
c     diagnose qcemat and qcerr.
c-----------------------------------------------------------------------
      OPEN(UNIT=coefs_unit,FILE="qc.out",STATUS="REPLACE")
      WRITE(coefs_unit,
     $     '(4(a,es10.3)/2(a,es10.3),a,2es11.3/a,es10.3,a,2es10.3/)')
     $     " e = ",rt%e,", f = ",rt%f,", h = ",rt%h,", di = ",rt%di,
     $     " g = ",rt%g,", k = ",rt%k,", r = ",rt%r,
     $     " r1 - r2 = ",rt%r(1)-rt%r(2),", q = ",rt%q
      WRITE(coefs_unit,'(a,g0/)')" kmax = ",kmax
      DO k=1,kmax+2
         WRITE(coefs_unit,'(a,g0/)')" k = ",k
         WRITE(order,'(a,g0,a)')"(",k,")"
         CALL inps_cmatwrite(k2mat(:,:,k),"k2mat"//order,coefs_unit) 
         CALL inps_cmatwrite(qmat(:,:,k),"qmat"//order,coefs_unit) 
         CALL inps_cmatwrite(cmat(:,:,k),"cmat"//order,coefs_unit) 
         CALL inps_cmatwrite(qcemat(:,:,k),"qcemat"//order,coefs_unit) 
         WRITE(order,'(a,i2.2,a)')"(",k,")"
         WRITE(coefs_unit,'(a,es10.3/)')" qcerr"//order//" = ",qcerr(k)
      ENDDO
      CLOSE(UNIT=coefs_unit)
c-----------------------------------------------------------------------
c     diagnose ymat and yerrmat.
c-----------------------------------------------------------------------
      OPEN(UNIT=coefs_unit,FILE="yz.out",STATUS="REPLACE")
      WRITE(coefs_unit,'(a,g0/)')" kmax = ",kmax
      idmat=ident(2)
      DO k=0,kmax
         yerrmat=MATMUL(ymat(:,:,k),(rmat-2*k*idmat))
         DO l=0,k
            yerrmat=yerrmat-MATMUL(dmat(:,:,l),ymat(:,:,k-l))
         ENDDO
         yerror=MAXVAL(ABS(yerrmat))
         WRITE(coefs_unit,'(a,g0/)')" k = ",k
         WRITE(order,'(a,g0,a)')"(",k,")"
         CALL inps_cmatwrite(zmat(:,:,k),"zmat"//order,coefs_unit) 
         CALL inps_cmatwrite(ymat(:,:,k),"ymat"//order,coefs_unit) 
         CALL inps_cmatwrite(yerrmat,"yerrmat"//order,coefs_unit) 
         WRITE(order,'(a,g0,a)')"(",k,")"
         WRITE(coefs_unit,'(a,es10.3/)')" yerror"//order//" = ",yerror
      ENDDO
      CLOSE(UNIT=coefs_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_coefs
c-----------------------------------------------------------------------
c     subprogram 8. inps_ua.
c     computes power series solutions, inps.pdf Eq. (62).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_ua(x,ua,dua,tflag)

      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(6,2), INTENT(OUT) :: ua
      COMPLEX(r8), DIMENSION(6,2), INTENT(OUT), OPTIONAL :: dua
      LOGICAL, INTENT(IN), OPTIONAL :: tflag

      REAL(r8) :: xfac,power
      REAL(r8), DIMENSION(4) :: rvec
      COMPLEX(r8), DIMENSION(4) :: yy,dyy
      COMPLEX(r8), DIMENSION(12) :: zz,dzz
      COMPLEX(r8), DIMENSION(4,2) :: p21,dp21
      COMPLEX(r8), DIMENSION(2,2) :: smat,dsmat,y,dy,q,dq,qsy,dqsy
      COMPLEX(r8), DIMENSION(4,0:kmax) :: cc
      COMPLEX(r8), DIMENSION(12,0:kmax) :: dd
c-----------------------------------------------------------------------
c     compute y and its derivatives with rvec.
c-----------------------------------------------------------------------
      xfac=1/x**2
      rvec=-(/rt%r(1),rt%r(1),rt%r(2),rt%r(2)/)/2
      cc=RESHAPE(ymat,(/4,kmax/))
      IF(PRESENT(dua))THEN
         CALL inps_horner(xfac,cc,yy,dyy,rvec)
         dy=-RESHAPE(dyy,(/2,2/))*2*xfac/x
      ELSE
         CALL inps_horner(xfac,cc,yy,rvec=rvec)
      ENDIF
      y=RESHAPE(yy,(/2,2/))
c-----------------------------------------------------------------------
c     compute q and p21 and their derivatives without rvec.
c-----------------------------------------------------------------------
      dd(1:4,:)=RESHAPE(qmat(:,:,0:kmax),(/4,kmax+1/))
      dd(5:12,:)=RESHAPE(pmat(r2,r1,0:kmax),(/8,kmax+1/))
      IF(PRESENT(dua))THEN
         CALL inps_horner(xfac,dd,zz,dzz)
         dzz=-dzz*2*xfac/x
         dq=RESHAPE(dzz(1:4),(/2,2/))
         dp21=RESHAPE(dzz(5:12),(/4,2/))
      ELSE
         CALL inps_horner(xfac,dd,zz)
      ENDIF
      q=RESHAPE(zz(1:4),(/2,2/))
      p21=RESHAPE(zz(5:12),(/4,2/))
c-----------------------------------------------------------------------
c     construct splitting matrix and its derivative.
c-----------------------------------------------------------------------
      pp=0
      dpp=0
      pp(r1,r1)=ident(2)
      pp(r2,r1)=p21
      dpp(r2,r1)=dp21
c-----------------------------------------------------------------------
c     back substitution for ua.
c-----------------------------------------------------------------------
      smat=RESHAPE((/one,zero,zero,xfac/),(/2,2/))
      qsy=MATMUL(q,MATMUL(smat,y))
      ua=MATMUL(pp,qsy)
      IF(PRESENT(tflag) .AND. tflag .OR. .NOT. PRESENT(tflag))
     $     ua=MATMUL(tmat,ua)
c-----------------------------------------------------------------------
c     back substitution for dua.
c-----------------------------------------------------------------------
      IF(PRESENT(dua))THEN
         dsmat=RESHAPE((/zero,zero,zero,-2*xfac/x/),(/2,2/))
         dqsy=MATMUL(q,MATMUL(smat,dy))
     $        +MATMUL(q,MATMUL(dsmat,y))
     $        +MATMUL(dq,MATMUL(smat,y))
         dua=MATMUL(pp,dqsy)+MATMUL(dpp,qsy)
         IF(PRESENT(tflag) .AND. tflag .OR. .NOT. PRESENT(tflag))
     $        dua=MATMUL(tmat,dua)
      ENDIF
c-----------------------------------------------------------------------
c     extract mercier power.
c-----------------------------------------------------------------------
      IF(extract)THEN
         power=rt%r(1)-1.5
         xfac=x**power
         ua(:,1)=ua(:,1)/xfac
         ua(:,2)=ua(:,2)*xfac
         IF(PRESENT(dua))THEN
            dua(:,1)=dua(:,1)/xfac
            dua(:,2)=dua(:,2)*xfac
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_ua
c-----------------------------------------------------------------------
c     subprogram 9. inps_horner.
c     horner's method for evaluation of polynomials and derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_horner(x,c,y,dy,rvec)

      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(:,0:), INTENT(IN) :: c
      COMPLEX(r8), DIMENSION(:), INTENT(OUT) :: y
      COMPLEX(r8), DIMENSION(:), INTENT(OUT), OPTIONAL :: dy
      REAL(r8), DIMENSION(:), INTENT(IN), OPTIONAL :: rvec

      INTEGER :: k,n
      REAL(r8), DIMENSION(SIZE(y)) :: xrvec
c-----------------------------------------------------------------------
c     compute y.
c-----------------------------------------------------------------------
      n=SIZE(c,2)-1
      y=c(:,n)
      DO k=n-1,0,-1
         y=y*x+c(:,k)
      ENDDO
      IF(PRESENT(rvec))THEN
         xrvec=x**rvec
         y=y*xrvec
      ENDIF
c-----------------------------------------------------------------------
c     compute dy.
c-----------------------------------------------------------------------
      IF(PRESENT(dy))THEN
         IF(PRESENT(rvec))THEN
            dy=c(:,n)*(rvec+n)
            DO k=n-1,0,-1
               dy=dy*x+c(:,k)*(rvec+k)
            ENDDO
            dy=dy*xrvec/x
         ELSE
            dy=c(:,n)*n
            DO k=n-1,0,-1
               dy=dy*x+c(:,k)*k
            ENDDO
            dy=dy/x
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_horner
c-----------------------------------------------------------------------
c     subprogram 10. inps_xmax.
c     computes xmax.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_xmax(eps1,eps2,x1,x2)

      REAL(r8), INTENT(IN) :: eps1,eps2
      REAL(r8), INTENT(OUT) :: x1,x2

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      LOGICAL :: set1
      CHARACTER(64) :: message
      REAL(r8), PARAMETER :: xlogmin=-1,dxlog=.01
      INTEGER :: ixlog,nxlog=1000
      REAL(r8) :: xlog,x,dxfac
      REAL(r8), DIMENSION(2) :: delta
c-----------------------------------------------------------------------
c     start loops over x.
c-----------------------------------------------------------------------
      set1=.TRUE.
      dxfac=10**dxlog
      xlog=xlogmin
      x=10**xlog
      ixlog=0
      DO
c-----------------------------------------------------------------------
c     compute delta and diagnose.
c-----------------------------------------------------------------------
         CALL inps_delta(x,delta)
         IF(diagnose)THEN
            WRITE(*,'(3(a,es10.3))')
     $           " e = ",rt%e,", f = ",rt%f,", h = ",rt%h
            WRITE(*,'(4(a,es10.3))')
     $           " g = ",rt%g,", k = ",rt%k,", q = ",REAL(rt%q)
            WRITE(*,'(a,es10.3,a,2es10.3)')" x = ",x,", delta = ",delta
            CALL program_stop("inps_xmax: abort after diagnose.")
         ENDIF
c-----------------------------------------------------------------------
c     compute x1 and x2.
c-----------------------------------------------------------------------
         IF(MAXVAL(delta) < eps1 .AND. set1)THEN
            x1=x
            set1=.FALSE.
         ENDIF
         IF(MAXVAL(delta) < eps2)THEN
            x2=x
            EXIT
         ENDIF
c-----------------------------------------------------------------------
c     abort if ixlog = nxlog.
c-----------------------------------------------------------------------
         IF(ixlog == nxlog)THEN
            WRITE(message,'(a,g0,a,es10.3)')
     $           "inps_xmax: abort with ixlog = nxlog = ",nxlog,
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
      END SUBROUTINE inps_xmax
c-----------------------------------------------------------------------
c     subprogram 11. inps_delta.
c     computes delta.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inps_delta(x,delta)
      
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(2), INTENT(OUT) :: delta

      INTEGER :: j
      REAL(r8) :: xfac
      REAL(r8), DIMENSION(2) :: xrfac
      REAL(r8), DIMENSION(2,0:2) :: norm
      COMPLEX(r8), DIMENSION(6,2,0:2) ::matvec
      COMPLEX(r8), DIMENSION(6,2) :: ua,dua
      COMPLEX(r8), DIMENSION(6,6) :: matrix
c-----------------------------------------------------------------------
c     compute ua and dua and remove mercier powers.
c-----------------------------------------------------------------------
      CALL inps_ua(x,ua,dua,.FALSE.)
      xrfac=x**rt%r
      DO j=1,2
         ua(:,j)=ua(:,j)/xrfac(j)
         dua(:,j)=dua(:,j)/xrfac(j)
      ENDDO
c-----------------------------------------------------------------------
c     compute matrix, matvec, and errmat.
c-----------------------------------------------------------------------
      xfac=1/x**2
      matrix=jmat(:,:,0)
      IF(kmax > 0)matrix=matrix+xfac*jmat(:,:,1)
      IF(kmax > 1)matrix=matrix+xfac*xfac*jmat(:,:,2)
c-----------------------------------------------------------------------
c     compute and apply norms.
c-----------------------------------------------------------------------
      matvec(:,:,1)=dua
      matvec(:,:,2)=-x*MATMUL(matrix,ua)
      matvec(:,:,0)=SUM(matvec(:,:,1:2),3)
c-----------------------------------------------------------------------
c     compute and apply norms.
c-----------------------------------------------------------------------
      DO j=1,2
         norm(j,:)=MAXVAL(ABS(matvec(:,j,:)),1)
         delta(j)=norm(j,0)/MAXVAL(norm(j,1:2))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inps_delta
      END MODULE inps_mod
