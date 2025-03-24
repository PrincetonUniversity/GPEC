c-----------------------------------------------------------------------
c     file spline.f
c     fits functions to cubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. spline_mod.
c     1. spline_alloc.
c     2. spline_dealloc.
c     3. spline_fit.
c     4. spline_fit_ahg.
c     5. spline_fit_classic.
c     6. spline_fit_ha.
c     7. spline_fac.
c     8. spline_eval.
c     9. spline_eval_external.
c     10. spline_all_eval.
c     11. spline_write1.
c     12. spline_write2.
c     13. spline_int.
c     14. spline_triluf.
c     15. spline_trilus.
c     16. spline_sherman.
c     17. spline_morrison.
c     18. spline_copy.
c     19. spline_fill_matrix.
c     20. spline_get_yp.
c     21. spline_thomas.
c     22. spline_roots.
c     23. spline_refine_root.
c-----------------------------------------------------------------------
c     subprogram 0. spline_type definition.
c     defines spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE spline_mod
      USE local_mod
      USE utils_mod, ONLY: bubble
      IMPLICIT NONE

      LOGICAL :: use_classic_splines = .FALSE.
      TYPE :: spline_type
      INTEGER :: mx,nqty,ix
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xs,f,f1,f2,f3
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: fs,fs1,fsi,xpower
      REAL(r8), DIMENSION(2) :: x0
      CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
      CHARACTER(6) :: name
      LOGICAL :: periodic=.FALSE., allocated=.FALSE.
      END TYPE spline_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. spline_alloc.
c     allocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_alloc(spl,mx,nqty)

      INTEGER, INTENT(IN) :: mx,nqty
      TYPE(spline_type), INTENT(INOUT) :: spl
      ! IF(spl%allocated) CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     safety check.
c-----------------------------------------------------------------------
      IF(spl%allocated)
     $   CALL program_stop("spline_dealloc: spline already allocated")
      
c-----------------------------------------------------------------------
c     set scalars.
c-----------------------------------------------------------------------
      spl%mx=mx
      spl%nqty=nqty
      spl%ix=0
      spl%periodic=.FALSE.
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(spl%xs(0:mx))
      ALLOCATE(spl%f(nqty))
      ALLOCATE(spl%f1(nqty))
      ALLOCATE(spl%f2(nqty))
      ALLOCATE(spl%f3(nqty))
      ALLOCATE(spl%title(0:nqty))
      ALLOCATE(spl%fs(0:mx,nqty))
      ALLOCATE(spl%fs1(0:mx,nqty))
      ALLOCATE(spl%xpower(2,nqty))
      spl%xpower=0
      spl%x0=0
      spl%allocated=.TRUE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_alloc
c-----------------------------------------------------------------------
c     subprogram 2. spline_dealloc.
c     deallocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_dealloc(spl)

      TYPE(spline_type), INTENT(INOUT) :: spl
c-----------------------------------------------------------------------
c     safety check.
c-----------------------------------------------------------------------
      IF(.NOT.spl%allocated)
     $   CALL program_stop("spline_dealloc: spline not allocated")
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(spl%xs)
      DEALLOCATE(spl%f)
      DEALLOCATE(spl%f1)
      DEALLOCATE(spl%f2)
      DEALLOCATE(spl%f3)
      DEALLOCATE(spl%title)
      DEALLOCATE(spl%fs)
      DEALLOCATE(spl%fs1)
      DEALLOCATE(spl%xpower)
      IF(ALLOCATED(spl%fsi))DEALLOCATE(spl%fsi)
      spl%allocated=.FALSE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. spline_fit.
c     switch between cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fit(spl,endmode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode
c-----------------------------------------------------------------------
c     switch between two spline_fit.
c-----------------------------------------------------------------------
      IF (use_classic_splines .AND.
     $    (endmode.EQ."extrap".OR.endmode.EQ."natural"))THEN
         CALL spline_fit_classic(spl,endmode)
      ELSE
         CALL spline_fit_ahg(spl,endmode)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fit
c-----------------------------------------------------------------------
c     subprogram 4. spline_fit_ahg.
c     fits real functions to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fit_ahg(spl,endmode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode

      INTEGER :: iqty,iside
      REAL(r8), DIMENSION(-1:1,0:spl%mx) :: a
      REAL(r8), DIMENSION(spl%mx) :: b
      REAL(r8), DIMENSION(4) :: cl,cr
      REAL(r8), DIMENSION(0:spl%mx) :: xfac
c-----------------------------------------------------------------------
c     extract powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(spl%xs-spl%x0(iside))**spl%xpower(iside,iqty)
               spl%fs(:,iqty)=spl%fs(:,iqty)*xfac
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set up grid matrix.
c-----------------------------------------------------------------------
      CALL spline_fac(spl,a,b,cl,cr,endmode)
c-----------------------------------------------------------------------
c     compute first derivatives, interior.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         spl%fs1(1:spl%mx-1,iqty)=
     $        3*((spl%fs(2:spl%mx,iqty)-spl%fs(1:spl%mx-1,iqty))
     $        *b(2:spl%mx)
     $        +(spl%fs(1:spl%mx-1,iqty)-spl%fs(0:spl%mx-2,iqty))
     $        *b(1:spl%mx-1))
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         DO iqty=1,spl%nqty
            spl%fs1(0,iqty)=SUM(cl(1:4)*spl%fs(0:3,iqty))
            spl%fs1(spl%mx,iqty)=SUM(cr(1:4)
     $           *spl%fs(spl%mx:spl%mx-3:-1,iqty))
            spl%fs1(1,iqty)=spl%fs1(1,iqty)-spl%fs1(0,iqty)
     $           /(spl%xs(1)-spl%xs(0))
            spl%fs1(spl%mx-1,iqty)=
     $           spl%fs1(spl%mx-1,iqty)-spl%fs1(spl%mx,iqty)
     $           /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
         ENDDO
         CALL spline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         spl%fs1(1,:)=spl%fs1(1,:)-(2*spl%fs(1,:)
     $        -spl%fs(0,:)-spl%fs(2,:))*2*b(1)
         spl%fs1(spl%mx-1,:)=spl%fs1(spl%mx-1,:)
     $        +(2*spl%fs(spl%mx-1,:)-spl%fs(spl%mx,:)
     $        -spl%fs(spl%mx-2,:))*2*b(spl%mx)
         CALL spline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
         spl%fs1(0,:)=(2*(2*spl%fs(1,:)-spl%fs(0,:)-spl%fs(2,:))
     $        +(spl%fs1(1,:)+spl%fs1(2,:))*(spl%xs(2)-spl%xs(1))
     $        -spl%fs1(1,:)*(spl%xs(1)-spl%xs(0)))/(spl%xs(1)-spl%xs(0))
         spl%fs1(spl%mx,:)=
     $        (2*(spl%fs(spl%mx-2,:)+spl%fs(spl%mx,:)
     $        -2*spl%fs(spl%mx-1,:))
     $        +(spl%fs1(spl%mx-1,:)+spl%fs1(spl%mx-2,:))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        -spl%fs1(spl%mx-1,:)
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-1)))
     $        /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boudary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         spl%periodic=.TRUE.
         spl%fs1(0,:)=3*((spl%fs(1,:)-spl%fs(0,:))*b(1)
     $        +(spl%fs(0,:)-spl%fs(spl%mx-1,:))*b(spl%mx))
         CALL spline_morrison(a(:,0:spl%mx-1),spl%fs1(0:spl%mx-1,:))
         spl%fs1(spl%mx,:)=spl%fs1(0,:)
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fit_ahg
c-----------------------------------------------------------------------
c     subprogram 5. spline_fit_classic.
c     classical way to solve spline coefficients
c     does not support periodic and not-a-knot boundary conditions
c     g=a+b(x-xi)+c(x-xi)^2+d(x-xi)^3
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fit_classic(spl,endmode)
      TYPE(spline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode
      REAL(r8), DIMENSION(:),ALLOCATABLE :: d,l,u,h
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: r
      REAL(r8), DIMENSION(0:spl%mx) :: xfac

      INTEGER :: iside,iqty,i
      REAL(r8),DIMENSION(spl%nqty) :: bs,cs,ds
c-----------------------------------------------------------------------
c     extract powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(spl%xs-spl%x0(iside))**spl%xpower(iside,iqty)
               spl%fs(:,iqty)=spl%fs(:,iqty)*xfac
            ENDIF
         ENDDO
      ENDDO
      ALLOCATE (d(0:spl%mx),l(spl%mx),u(spl%mx),r(0:spl%mx,spl%nqty))
      ALLOCATE (h(0:spl%mx-1))
c-----------------------------------------------------------------------
c     compute tridiagnol matrix for natural B.C.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         h(i)=spl%xs(i+1)-spl%xs(i)
      ENDDO

      d(0)=1
      DO i=1,spl%mx-1
         d(i)=2*(h(i-1)+h(i))
      ENDDO
      d(spl%mx)=1

      DO i=1,spl%mx-1
         l(i)=h(i-1)
      ENDDO
      l(spl%mx)=0

      u(1)=0
      DO i=2,spl%mx
         u(i)=h(i-1)
      ENDDO

      r(0,:)=0
      DO i=1,spl%mx-1
         r(i,:)=( (spl%fs(i+1,:)-spl%fs(i,:))/h(i)
     $           -(spl%fs(i,:)-spl%fs(i-1,:))/h(i-1) )*6
      ENDDO
      r(spl%mx,:)=0

      IF (endmode=="extrap") THEN
         CALL spline_get_yp(spl%xs(0:3),spl%fs(0:3,:),
     $                      spl%xs(0),r(0,:),spl%nqty)
         CALL spline_get_yp(spl%xs(spl%mx-3:spl%mx),
     $        spl%fs(spl%mx-3:spl%mx,:),spl%xs(spl%mx),
     $        r(spl%mx,:),spl%nqty)
         d(0)=2*h(0)
         d(spl%mx)=2*h(spl%mx-1)
         u(1)=h(0)
         l(spl%mx)=h(spl%mx-1)
         r(0,:)=( (spl%fs(1,:)-spl%fs(0,:))/h(0) - r(0,:) )*6
         r(spl%mx,:)=( r(spl%mx,:)
     $  -(spl%fs(spl%mx,:)-spl%fs(spl%mx-1,:))/h(spl%mx-1) )*6

      ENDIF
c-----------------------------------------------------------------------
c     solve and contrruct spline.
c-----------------------------------------------------------------------

      CALL spline_thomas(l,d,u,r,spl%mx+1,spl%nqty)

      DO i=0, spl%mx-1
         bs=(spl%fs(i+1,:)-spl%fs(i,:))/h(i)
     $    - 0.5*h(i)*r(i,:)
     $    - h(i)*(r(i+1,:)-r(i,:))/6
         spl%fs1(i,:)=bs
      ENDDO
      ds=(r(spl%mx,:)-r(spl%mx-1,:))/(h(spl%mx-1)*6)
      cs=r(spl%mx-1,:)*0.5
      i=spl%mx-1
      bs=(spl%fs(i+1,:)-spl%fs(i,:))/h(i)
     $    - 0.5*h(i)*r(i,:)
     $    - h(i)*(r(i+1,:)-r(i,:))/6
      i=spl%mx
      spl%fs1(i,:)=bs+h(i-1)*(cs*2+h(i-1)*ds*3)
      DEALLOCATE (d,l,u,r,h)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fit_classic
c-----------------------------------------------------------------------
c     subprogram 6. spline_fit_ha.
c     fits real functions to highly accurate cubic splines
c     using generalized matrix solves: sacrificing speed for accuracy
c     and generalizability
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fit_ha(spl,endmode)
      TYPE(spline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode

      INTEGER ::icount,icoef,imx,iqty,istart,jstart,info,iside
      INTEGER :: ndim,nqty,kl,ku,ldab,nvar
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      INTEGER,  DIMENSION(:,:), ALLOCATABLE :: imap
      REAL(r8) :: x0,x1,x2,dx
      REAL(r8), DIMENSION(0:spl%mx) :: xfac
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rhs,locrhs,fs,locrhs0,
     $                                         locrhs1,tmpmat
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: mat,locmat,locmat0,locmat1
      REAL(r8), DIMENSION(:,:,:),ALLOCATABLE :: coef
c-----------------------------------------------------------------------
c     extract powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(spl%xs-spl%x0(iside))**spl%xpower(iside,iqty)
               spl%fs(:,iqty)=spl%fs(:,iqty)*xfac
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     contruct grid and related matrix block.
c     five unknow variables at each interval a,b,c,e,f
c     a(x-xi)^3+b(x-xi)^2+c(x-xi)+d, where d=yi
c     e and f are used for passing periodic conditions
c-----------------------------------------------------------------------
      nqty=spl%nqty
      nvar=5
      ALLOCATE (coef(nvar,0:spl%mx-1,nqty),imap(nvar,0:spl%mx-1))
      ALLOCATE (fs(0:spl%mx,nqty))
      icount=1
      DO imx=0, spl%mx-1
         DO icoef=1,nvar
            imap(icoef,imx)=icount
            icount=icount+1
         ENDDO
      ENDDO
      ndim=icount-1
      kl=nvar*3
      ku=kl
      ldab=2*kl+ku+1
      ALLOCATE(mat(ldab,ndim),rhs(ndim,nqty),ipiv(ndim))
      mat=0
      rhs=0
      fs=spl%fs(:,:)
c-----------------------------------------------------------------------
c     contruct local matrix in each interval.
c-----------------------------------------------------------------------
      ALLOCATE(locmat(nvar,nvar*3),locmat0(nvar,nvar*3),
     &         locmat1(nvar,nvar*2),tmpmat(nvar,nvar*2))
      ALLOCATE(locrhs(nvar,nqty),locrhs0(nvar,nqty),locrhs1(nvar,nqty))
      DO imx=1,spl%mx-2
         x0=spl%xs(imx-1)
         x1=spl%xs(imx)
         x2=spl%xs(imx+1)
         locmat=0
         locrhs=0
c-----------------------------------------------------------------------
c     S'0(x1)=S'1(x1)
c-----------------------------------------------------------------------
         dx=x1-x0
         locmat(1,1)=3*dx*dx
         locmat(1,2)=2*dx
         locmat(1,3)=1
         locmat(1,nvar+3)=-1
c-----------------------------------------------------------------------
c     S''1(x2)=S''2(x2)
c-----------------------------------------------------------------------
         dx=x2-x1
         locmat(2,nvar+1)=6*dx
         locmat(2,nvar+2)=2
         locmat(2,2*nvar+2)=-2
c-----------------------------------------------------------------------
c     S1(x2)=S2(x2)
c-----------------------------------------------------------------------
         locmat(3,nvar+1)=dx*dx*dx
         locmat(3,nvar+2)=dx*dx
         locmat(3,nvar+3)=dx
         locrhs(3,:)=fs(imx+1,:)-fs(imx,:)
c-----------------------------------------------------------------------
c     e1=e2
c-----------------------------------------------------------------------
         locmat(4,nvar+4)=1
         locmat(4,2*nvar+4)=-1
c-----------------------------------------------------------------------
c     f0=f1
c-----------------------------------------------------------------------
         locmat(5,5)=1
         locmat(5,nvar+5)=-1
c-----------------------------------------------------------------------
c     fill global matrix
c-----------------------------------------------------------------------
         istart=imap(1,imx)
         jstart=imap(1,imx-1)
         CALL spline_fill_matrix(mat,locmat,istart,jstart,kl,ku)
         rhs(istart:istart+nvar-1,:)=locrhs
      ENDDO
c-----------------------------------------------------------------------
c     boundary condition
c-----------------------------------------------------------------------
      locmat0=0
      locrhs0=0
      dx=spl%xs(1)-spl%xs(0)

      locmat0(2,nvar+1)=6*dx
      locmat0(2,nvar+2)=2
      locmat0(2,2*nvar+2)=-2

      locmat0(3,nvar+1)=dx*dx*dx
      locmat0(3,nvar+2)=dx*dx
      locmat0(3,nvar+3)=dx
      locrhs0(3,:)=fs(1,:)-fs(0,:)

      locmat0(4,nvar+4)=1
      locmat0(4,2*nvar+4)=-1

      locmat1=0
      locrhs1=0
      dx=spl%xs(spl%mx-1)-spl%xs(spl%mx-2)

      locmat1(1,1)=3*dx*dx
      locmat1(1,2)=2*dx
      locmat1(1,3)=1
      locmat1(1,nvar+3)=-1

      dx=spl%xs(spl%mx)-spl%xs(spl%mx-1)
      locmat1(2,nvar+1)=dx*dx*dx
      locmat1(2,nvar+2)=dx*dx
      locmat1(2,nvar+3)=dx
      locrhs1(2,:)=fs(spl%mx,:)-fs(spl%mx-1,:)

      locmat1(5,5)=1
      locmat1(5,nvar+5)=-1

      SELECT CASE(endmode)
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         locmat0(1,nvar+1)=6
         locmat0(1,2*nvar+1)=-6
         locmat0(5,nvar+5)=1
         locrhs0(5,:)=1

         locmat1(3,1)=6
         locmat1(3,nvar+1)=-6

         locmat1(4,nvar+4)=1
         locrhs1(4,:)=1
c-----------------------------------------------------------------------
c     extrap boudary conditions, use first and last four points to
c     calculate y'(0) and y'(1).
c-----------------------------------------------------------------------
      CASE("extrap")
         locmat0(1,nvar+3)=1
         CALL spline_get_yp(spl%xs(0:3),spl%fs(0:3,:),
     $                      spl%xs(0),locrhs0(1,:),nqty)

         locmat0(5,nvar+5)=1
         locrhs0(5,:)=1

         dx=spl%xs(spl%mx)-spl%xs(spl%mx-1)
         locmat1(3,nvar+1)=3*dx*dx
         locmat1(3,nvar+2)=2*dx
         locmat1(3,nvar+3)=1
         CALL spline_get_yp(spl%xs(spl%mx-3:spl%mx),
     $   spl%fs(spl%mx-3:spl%mx,:),spl%xs(spl%mx),locrhs1(3,:),nqty)

         locmat1(4,nvar+4)=1
         locrhs1(4,:)=1
c-----------------------------------------------------------------------
c     natural boudary conditions.
c-----------------------------------------------------------------------
      CASE("natural")
         locmat0(1,nvar+2)=2

         locmat0(5,nvar+5)=1
         locrhs0(5,:)=1

         dx=spl%xs(spl%mx)-spl%xs(spl%mx-1)
         locmat1(3,nvar+1)=6*dx
         locmat1(3,nvar+2)=2

         locmat1(4,nvar+4)=1
         locrhs1(4,:)=1
c-----------------------------------------------------------------------
c     periodic boudary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
c-----------------------------------------------------------------------
c     s'0(x0)=s'm-1(xm).
c-----------------------------------------------------------------------
         DO iqty=1,nqty
            IF (ABS(spl%fs(0,iqty)-spl%fs(spl%mx,iqty)) > 1E-15) THEN
               WRITE(*,*)
     $             "Warning: first and last points are different.
     $              IQTY= ",IQTY,",  averaged value is used."//
     $              TRIM(endmode)
              spl%fs(0,iqty)=(spl%fs(0,iqty)+spl%fs(spl%mx,iqty))*0.5
              spl%fs(spl%mx,iqty)=spl%fs(0,iqty)
            ENDIF
         ENDDO
         locmat0(1,nvar+3)=1
         locmat0(1,nvar+4)=-1

         dx=spl%xs(spl%mx)-spl%xs(spl%mx-1)
         locmat1(4,nvar+1)=3*dx*dx
         locmat1(4,nvar+2)=2*dx
         locmat1(4,nvar+3)=1
         locmat1(4,nvar+4)=-1
c-----------------------------------------------------------------------
c     s''0(x0)=s''m-1(xm).
c-----------------------------------------------------------------------
         locmat1(3,nvar+1)=6*dx
         locmat1(3,nvar+1)=2
         locmat1(3,nvar+5)=-1

         locmat0(5,nvar+2)=2
         locmat0(5,nvar+5)=-1
         spl%periodic=.TRUE.
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop
     $       ("Cannot recognize endmode = "//TRIM(endmode))
      END SELECT
c-----------------------------------------------------------------------
c     fill global matrix at x0
c-----------------------------------------------------------------------
      istart=imap(1,0)
      jstart=imap(1,0)
      tmpmat=locmat0(:,nvar+1:3*nvar)
      CALL spline_fill_matrix(mat,tmpmat,istart,jstart,kl,ku)
      rhs(istart:istart+nvar-1,:)=locrhs0
c-----------------------------------------------------------------------
c     fill global matrix at xm
c-----------------------------------------------------------------------
      istart=imap(1,spl%mx-1)
      jstart=imap(1,spl%mx-2)
      tmpmat=locmat1(:,1:2*nvar)
      CALL spline_fill_matrix(mat,tmpmat,istart,jstart,kl,ku)
      rhs(istart:istart+nvar-1,:)=locrhs1
c-----------------------------------------------------------------------
c     solve global matrix
c-----------------------------------------------------------------------
      CALL dgbtrf(ndim,ndim,kl,ku,mat,ldab,ipiv,info)
      IF (info .NE. 0) THEN
         CALL program_stop
     $           ("Error: LU factorization of spline matrix info ne 0")
      ENDIF
      CALL dgbtrs("N",ndim,kl,ku,nqty,mat,ldab,ipiv,rhs,ndim,info )
      IF (info .NE. 0) THEN
         CALL program_stop
     $        ("Error: solve spline matrix info ne 0")
      ENDIF
c-----------------------------------------------------------------------
c     get spl%fs1.
c-----------------------------------------------------------------------
      DO imx=0, spl%mx-1
         DO icoef=1,nvar
            coef(icoef,imx,:)=rhs(imap(icoef,imx),:)
         ENDDO
         spl%fs1(imx,:)=coef(3,imx,:)
      ENDDO
      dx=spl%xs(imx)-spl%xs(imx-1)
      spl%fs1(imx,:)=coef(3,imx-1,:)
     $                   +dx*(coef(2,imx-1,:)*2+dx*coef(1,imx-1,:)*3)

      DEALLOCATE (coef,imap)
      DEALLOCATE (fs)
      DEALLOCATE(mat,rhs,ipiv)
      DEALLOCATE(locmat,locmat0,locmat1,tmpmat)
      DEALLOCATE(locrhs,locrhs0,locrhs1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fit_ha
c-----------------------------------------------------------------------
c     subprogram 7. spline_fac.
c     sets up matrix for cubic spline fitting.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fac(spl,a,b,cl,cr,endmode)

      TYPE(spline_type), INTENT(IN) :: spl
      REAL(r8), DIMENSION(-1:1,0:spl%mx), INTENT(OUT) :: a
      REAL(r8), DIMENSION(spl%mx), INTENT(OUT) :: b
      REAL(r8), DIMENSION(4), INTENT(OUT) :: cl,cr
      CHARACTER(*), INTENT(IN) :: endmode

      INTEGER :: j
c-----------------------------------------------------------------------
c     compute interior matrix.
c-----------------------------------------------------------------------
      b=1/(spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1))
      DO j=1,spl%mx-1
         a(-1,j)=b(j)
         a(0,j)=2*(b(j)+b(j+1))
         a(1,j)=b(j+1)
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         b=b*b
         cl(1)=(spl%xs(0)*(3*spl%xs(0)
     $        -2*(spl%xs(1)+spl%xs(2)+spl%xs(3)))
     $        +spl%xs(1)*spl%xs(2)+spl%xs(1)*spl%xs(3)
     $        +spl%xs(2)*spl%xs(3))
     $        /((spl%xs(0)-spl%xs(1))*(spl%xs(0)-spl%xs(2))
     $        *(spl%xs(0)-spl%xs(3)))
         cl(2)=((spl%xs(2)-spl%xs(0))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(1)-spl%xs(0))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(1)-spl%xs(3)))
         cl(3)=((spl%xs(0)-spl%xs(1))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(0)-spl%xs(2))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(3)-spl%xs(2)))
         cl(4)=((spl%xs(1)-spl%xs(0))*(spl%xs(2)-spl%xs(0)))
     $        /((spl%xs(3)-spl%xs(0))*(spl%xs(3)-spl%xs(1))
     $        *(spl%xs(3)-spl%xs(2)))
         cr(1)=(spl%xs(spl%mx)*(3*spl%xs(spl%mx)
     $        -2*(spl%xs(spl%mx-1)+spl%xs(spl%mx-2)
     $        +spl%xs(spl%mx-3)))+spl%xs(spl%mx-1)
     $        *spl%xs(spl%mx-2)+spl%xs(spl%mx-1)*spl%xs(spl%mx
     $        -3)+spl%xs(spl%mx-2)*spl%xs(spl%mx-3))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-2))*(spl%xs(spl%mx
     $        )-spl%xs(spl%mx-3)))
         cr(2)=((spl%xs(spl%mx-2)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-3)))
         cr(3)=((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         cr(4)=((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-2)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-3)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         CALL spline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         b=b*b
         a(0,1)=a(0,1)+(spl%xs(2)+spl%xs(0)-2*spl%xs(1))*b(1)
         a(1,1)=a(1,1)+(spl%xs(2)-spl%xs(1))*b(1)
         a(0,spl%mx-1)=a(0,spl%mx-1)
     $        +(2*spl%xs(spl%mx-1)-spl%xs(spl%mx-2)
     $        -spl%xs(spl%mx))*b(spl%mx)
         a(-1,spl%mx-1)=a(-1,spl%mx-1)
     $        +(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))*b(spl%mx)
         CALL spline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boundary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         a(0,0:spl%mx:spl%mx)=2*(b(spl%mx)+b(1))
         a(1,0)=b(1)
         a(-1,0)=b(spl%mx)
         b=b*b
         CALL spline_sherman(a(:,0:spl%mx-1))
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fac
c-----------------------------------------------------------------------
c     subprogram 8. spline_eval.
c     evaluates real cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_eval(spl,x,mode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: mode

      INTEGER :: iqty,iside
      REAL(r8) :: xx,d,z,z1,xfac,dx
      REAL(r8) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      xx=x
      spl%ix=MAX(spl%ix,0)
      spl%ix=MIN(spl%ix,spl%mx-1)
c-----------------------------------------------------------------------
c     normalize interval for periodic splines.
c-----------------------------------------------------------------------
      IF(spl%periodic)THEN
         DO
            IF(xx < spl%xs(spl%mx))EXIT
            xx=xx-spl%xs(spl%mx)
         ENDDO
         DO
            IF(xx >= spl%xs(0))EXIT
            xx=xx+spl%xs(spl%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find cubic spline interval.
c-----------------------------------------------------------------------
      DO
         IF(xx >= spl%xs(spl%ix).OR.spl%ix <= 0)EXIT
         spl%ix=spl%ix-1
      ENDDO
      DO
         IF(xx < spl%xs(spl%ix+1).OR.spl%ix >= spl%mx-1)EXIT
         spl%ix=spl%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      d=spl%xs(spl%ix+1)-spl%xs(spl%ix)
      z=(xx-spl%xs(spl%ix))/d
      z1=1-z
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      spl%f=spl%fs(spl%ix,1:spl%nqty)*z1*z1*(3-2*z1)
     $     +spl%fs(spl%ix+1,1:spl%nqty)*z*z*(3-2*z)
     $     +d*z*z1*(spl%fs1(spl%ix,1:spl%nqty)*z1
     $     -spl%fs1(spl%ix+1,1:spl%nqty)*z)
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         spl%f1=6*(spl%fs(spl%ix+1,1:spl%nqty)
     $        -spl%fs(spl%ix,1:spl%nqty))*z*z1/d
     $        +spl%fs1(spl%ix,1:spl%nqty)*z1*(3*z1-2)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)*z*(3*z-2)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         spl%f2=(6*(spl%fs(spl%ix+1,1:spl%nqty)
     $        -spl%fs(spl%ix,1:spl%nqty))*(z1-z)/d
     $        -spl%fs1(spl%ix,1:spl%nqty)*(6*z1-2)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)*(6*z-2))/d
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         spl%f3=(12*(spl%fs(spl%ix,1:spl%nqty)
     $        -spl%fs(spl%ix+1,1:spl%nqty))/d
     $        +6*(spl%fs1(spl%ix,1:spl%nqty)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)))/(d*d)
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=spl%f(iqty)*xfac
            IF(mode > 0)g1=(spl%f1(iqty)+spl%f(iqty)
     $           *spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(spl%f2(iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*spl%f1(iqty)+(spl%xpower(iside,iqty)-1)
     $           *spl%f(iqty)/dx))*xfac
            IF(mode > 2)g3=(spl%f3(iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*spl%f2(iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*spl%f1(iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *spl%f(iqty))))*xfac
            spl%f(iqty)=g
            IF(mode > 0)spl%f1(iqty)=g1
            IF(mode > 1)spl%f2(iqty)=g2
            IF(mode > 2)spl%f3(iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_eval
c-----------------------------------------------------------------------
c     subprogram 9. spline_eval_external
c     evaluates real cubic spline with external arrays (parallel).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_eval_external(spl,x,s_ix,s_f,s_f1,s_f2,s_f3)

      TYPE(spline_type), INTENT(IN) :: spl
      REAL(r8), INTENT(IN) :: x

      INTEGER :: iqty,iside
      REAL(r8) :: xx,d,z,z1,xfac,dx
      REAL(r8) :: g,g1,g2,g3

      INTEGER, INTENT(INOUT) :: s_ix
      REAL(r8), DIMENSION(:), INTENT(OUT) :: s_f
      REAL(r8), DIMENSION(:), OPTIONAL, INTENT(OUT) :: s_f1,s_f2,s_f3
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      xx=x
      s_ix=MAX(s_ix,0)
      s_ix=MIN(s_ix,spl%mx-1)
c-----------------------------------------------------------------------
c     normalize interval for periodic splines.
c-----------------------------------------------------------------------
      IF(spl%periodic)THEN
         DO
            IF(xx < spl%xs(spl%mx))EXIT
            xx=xx-spl%xs(spl%mx)
         ENDDO
         DO
            IF(xx >= spl%xs(0))EXIT
            xx=xx+spl%xs(spl%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find cubic spline interval.
c-----------------------------------------------------------------------
      DO
         IF(xx >= spl%xs(s_ix).OR.s_ix <= 0)EXIT
         s_ix=s_ix-1
      ENDDO
      DO
         IF(xx < spl%xs(s_ix+1).OR.s_ix >= spl%mx-1)EXIT
         s_ix=s_ix+1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      d=spl%xs(s_ix+1)-spl%xs(s_ix)
      z=(xx-spl%xs(s_ix))/d
      z1=1-z
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      s_f=spl%fs(s_ix,1:spl%nqty)*z1*z1*(3-2*z1)
     $     +spl%fs(s_ix+1,1:spl%nqty)*z*z*(3-2*z)
     $     +d*z*z1*(spl%fs1(s_ix,1:spl%nqty)*z1
     $     -spl%fs1(s_ix+1,1:spl%nqty)*z)
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(PRESENT(s_f1))THEN
         s_f1=6*(spl%fs(s_ix+1,1:spl%nqty)
     $        -spl%fs(s_ix,1:spl%nqty))*z*z1/d
     $        +spl%fs1(s_ix,1:spl%nqty)*z1*(3*z1-2)
     $        +spl%fs1(s_ix+1,1:spl%nqty)*z*(3*z-2)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(PRESENT(s_f2))THEN
         s_f2=(6*(spl%fs(s_ix+1,1:spl%nqty)
     $        -spl%fs(s_ix,1:spl%nqty))*(z1-z)/d
     $        -spl%fs1(s_ix,1:spl%nqty)*(6*z1-2)
     $        +spl%fs1(s_ix+1,1:spl%nqty)*(6*z-2))/d
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(PRESENT(s_f3))THEN
         s_f3=(12*(spl%fs(s_ix,1:spl%nqty)
     $        -spl%fs(s_ix+1,1:spl%nqty))/d
     $        +6*(spl%fs1(s_ix,1:spl%nqty)
     $        +spl%fs1(s_ix+1,1:spl%nqty)))/(d*d)
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=s_f(iqty)*xfac
            IF(PRESENT(s_f1))g1=(s_f1(iqty)+s_f(iqty)
     $           *spl%xpower(iside,iqty)/dx)*xfac
            IF(PRESENT(s_f2))g2=(s_f2(iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*s_f1(iqty)+(spl%xpower(iside,iqty)-1)
     $           *s_f(iqty)/dx))*xfac
            IF(PRESENT(s_f3))g3=(s_f3(iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*s_f2(iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*s_f1(iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *s_f(iqty))))*xfac
            s_f(iqty)=g
            IF(PRESENT(s_f1))s_f1(iqty)=g1
            IF(PRESENT(s_f2))s_f2(iqty)=g2
            IF(PRESENT(s_f3))s_f3(iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_eval_external
c-----------------------------------------------------------------------
c     subprogram 10. spline_all_eval.
c     evaluates cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_all_eval(spl,z,f,f1,f2,f3,mode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: z
      REAL(r8), DIMENSION(spl%mx,spl%nqty), INTENT(OUT) :: f,f1,f2,f3
      INTEGER, INTENT(IN) :: mode

      INTEGER :: iqty,nqty,n,iside
      REAL(r8) :: z1
      REAL(r8), DIMENSION(spl%mx) :: d,xfac,dx
      REAL(r8), DIMENSION(spl%mx) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      n=spl%mx
      nqty=spl%nqty
      z1=1-z
      d=spl%xs(1:n)-spl%xs(0:n-1)
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         f(:,iqty)=spl%fs(0:n-1,iqty)*z1*z1*(3-2*z1)
     $        +spl%fs(1:n,iqty)*z*z*(3-2*z)
     $        +d*z*z1*(spl%fs1(0:n-1,iqty)*z1-spl%fs1(1:n,iqty)*z)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         DO iqty=1,nqty
            f1(:,iqty)=6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*z*z1/d
     $           +spl%fs1(0:n-1,iqty)*z1*(3*z1-2)
     $           +spl%fs1(1:n,iqty)*z*(3*z-2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         DO iqty=1,nqty
            f2(:,iqty)=(6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*(z1-z)/d
     $           -spl%fs1(0:n-1,iqty)*(6*z1-2)
     $           +spl%fs1(1:n,iqty)*(6*z-2))/d
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         DO iqty=1,nqty
            f3(:,iqty)=(12*(spl%fs(0:n-1,iqty)-spl%fs(1:n,iqty))/d
     $           +6*(spl%fs1(0:n-1,iqty)+spl%fs1(1:n,iqty)))/(d*d)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=(spl%xs(0:spl%mx-1)+z*d(1:spl%mx))-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=f(:,iqty)*xfac
            IF(mode > 0)g1=(f1(:,iqty)
     $           +f(:,iqty)*spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(f2(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*f1(:,iqty)+(spl%xpower(iside,iqty)-1)
     $           *f(:,iqty)/dx))*xfac
     $
            IF(mode > 2)g3=(f3(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*f2(:,iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*f1(:,iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *f(:,iqty))))*xfac
            f(:,iqty)=g
            IF(mode > 0)f1(:,iqty)=g1
            IF(mode > 1)f2(:,iqty)=g2
            IF(mode > 2)f3(:,iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_all_eval
c-----------------------------------------------------------------------
c     subprogram 11. spline_write1.
c     produces ascii and binary output for real cubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_write1(spl,out,bin,iua,iub,interp)

      TYPE(spline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      CHARACTER(30) :: format1,format2
      INTEGER :: i,j
      REAL(r8) :: x,dx
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"ix",',i2.2,'(4x,a6,1x)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
 30   FORMAT('(/3x,"ix",2x,"j",',i2.2,'(4x,a6,1x)/)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN
c-----------------------------------------------------------------------
c     print node values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(format1,10)spl%nqty+1
         WRITE(format2,20)spl%nqty+1
         WRITE(iua,'(/1x,a)')'node values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      DO i=0,spl%mx
         CALL spline_eval(spl,spl%xs(i),0)
         IF(out)WRITE(iua,format2)i,spl%xs(i),spl%f
         IF(bin)WRITE(iub)REAL(spl%xs(i),4),REAL(spl%f,4)
      ENDDO
      IF(out)WRITE(iua,format1)spl%title(0:spl%nqty)
      IF(bin)WRITE(iub)
      IF(.NOT. interp)RETURN
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,'(/1x,a)')'interpolated values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         dx=(spl%xs(i+1)-spl%xs(i))/4
         DO j=0,4
            x=spl%xs(i)+j*dx
            CALL spline_eval(spl,x,0)
            IF(out)WRITE(iua,format2)i,x,spl%f
            IF(bin)WRITE(iub)REAL(x,4),REAL(spl%f,4)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
      x=spl%xs(spl%mx)
      CALL spline_eval(spl,x,0)
      IF(out)THEN
         WRITE(iua,format2)i,x,spl%f
         WRITE(iua,format1)spl%title
      ENDIF
      IF(bin)THEN
         WRITE(iub)REAL(x,4),REAL(spl%f,4)
         WRITE(iub)
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_write1
c-----------------------------------------------------------------------
c     subprogram 12. spline_write2.
c     produces ascii and binary output for real cubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_write2(spl,out,bin,iua,iub,interp)

      TYPE(spline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      CHARACTER(30) :: format1,format2
      INTEGER :: i,j,iz
      REAL(r8) :: x,dx,z
      REAL(r8), DIMENSION(spl%mx,spl%nqty) :: f
      REAL(r8), DIMENSION(0:4*spl%mx,spl%nqty) :: g
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"i",',i2.2,'(4x,a6,1x)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
 30   FORMAT('(/4x,"i",2x,"j",',i2.2,'(4x,a6,1x)/)')
c-----------------------------------------------------------------------
c     compute values.
c-----------------------------------------------------------------------
      z=0
      DO iz=0,3
         CALL spline_all_eval(spl,z,f,f,f,f,0)
         g(iz:4*spl%mx-1:4,:)=f
         z=z+.25_r8
      ENDDO
      CALL spline_eval(spl,spl%xs(spl%mx),0)
      g(4*spl%mx,:)=spl%f
c-----------------------------------------------------------------------
c     print node values.
c-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN
      IF(out)THEN
         WRITE(format1,10)spl%nqty+1
         WRITE(format2,20)spl%nqty+1
         WRITE(iua,'(/1x,a)')'node values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      DO i=0,spl%mx
         IF(out)WRITE(iua,format2)i,spl%xs(i),g(4*i,:)
         IF(bin)WRITE(iub)REAL(spl%xs(i),4),REAL(g(4*i,:),4)
      ENDDO
      IF(out)WRITE(iua,format1)spl%title(0:spl%nqty)
      IF(bin)WRITE(iub)
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
      IF(.NOT. interp)RETURN
      IF(out)THEN
         WRITE(iua,'(/1x,a)')'interpolated values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         dx=(spl%xs(i+1)-spl%xs(i))/4
         DO j=0,4
            x=spl%xs(i)+j*dx
            IF(out)WRITE(iua,format2)i,x,g(4*i+j,:)
            IF(bin)WRITE(iub)REAL(x,4),REAL(g(4*i+j,:),4)
         ENDDO
         IF(out)WRITE(iua,'(1x)')
      ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,format2)i,x,g(spl%mx*4,:)
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      IF(bin)THEN
         WRITE(iub)REAL(x,4),REAL(g(spl%mx*4,:),4)
         WRITE(iub)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_write2
c-----------------------------------------------------------------------
c     subprogram 13. spline_int.
c     integrates real cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_int(spl)

      TYPE(spline_type), INTENT(INOUT) :: spl

      INTEGER :: ix,iqty,ig
      REAL(r8), DIMENSION(spl%mx) :: dx
      REAL(r8), DIMENSION(spl%mx,spl%nqty) :: term,f

      INTEGER, PARAMETER :: mg=4
      REAL(r8), DIMENSION(mg) :: xg=(1+(/-0.861136311594053_r8,
     $     -0.339981043584856_r8,0.339981043584856_r8,
     $     0.861136311594053_r8/))/2
      REAL(r8), DIMENSION(mg) :: wg=(/0.347854845137454_r8,
     $     0.652145154862546_r8,0.652145154862546_r8,
     $     0.347854845137454_r8/)/2
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      IF(.NOT.ALLOCATED(spl%fsi))ALLOCATE(spl%fsi(0:spl%mx,spl%nqty))
      dx=spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1)
      term=0
c-----------------------------------------------------------------------
c     compute integrals over intervals.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         IF(spl%xpower(1,iqty) == 0 .AND. spl%xpower(2,iqty) == 0)THEN
            term(:,iqty)=dx/12
     $           *(6*(spl%fs(0:spl%mx-1,iqty)+spl%fs(1:spl%mx,iqty))
     $           +dx*(spl%fs1(0:spl%mx-1,iqty)-spl%fs1(1:spl%mx,iqty)))
         ELSE
            DO ig=1,mg
               CALL spline_all_eval(spl,xg(ig),f,f,f,f,0)
               term(:,iqty)=term(:,iqty)+dx*wg(ig)*f(:,iqty)
            ENDDO
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     accumulate over intervals.
c-----------------------------------------------------------------------
      spl%fsi(0,:)=0
      DO ix=1,spl%mx
         spl%fsi(ix,:)=spl%fsi(ix-1,:)+term(ix,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_int
c-----------------------------------------------------------------------
c     subprogram 14. spline_triluf.
c     performs tridiagonal LU factorization.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_triluf(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: i,j,k,jmin,jmax,n
c-----------------------------------------------------------------------
c     begin loop over rows and define limits.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         jmin=MAX(1-i,-1)
         jmax=MIN(n-i,1)
c-----------------------------------------------------------------------
c     compute lower elements.
c-----------------------------------------------------------------------
         DO j=jmin,-1
            DO k=MAX(jmin,j-1),j-1
               a(j,i)=a(j,i)-a(k,i)*a(j-k,i+k)
            ENDDO
            a(j,i)=a(j,i)*a(0,i+j)
         ENDDO
c-----------------------------------------------------------------------
c     compute diagonal element
c-----------------------------------------------------------------------
         DO k=MAX(jmin,-1),-1
            a(0,i)=a(0,i)-a(k,i)*a(-k,i+k)
         ENDDO
         a(0,i)=1/a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_triluf
c-----------------------------------------------------------------------
c     subprogram 15. spline_trilus.
c     performs tridiagonal LU solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_trilus(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: i,j,n
c-----------------------------------------------------------------------
c     down sweep.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         DO j=MAX(1-i,-1),-1
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     up sweep.
c-----------------------------------------------------------------------
      DO i=n,1,-1
         DO j=1,MIN(n-i,1)
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
         x(i,:)=x(i,:)*a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_trilus
c-----------------------------------------------------------------------
c     subprogram 16. spline_sherman.
c     uses Sherman-Morrison formula to factor periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_sherman(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: j,n
      REAL(r8), DIMENSION(SIZE(a,2),1) :: u
c-----------------------------------------------------------------------
c     prepare matrices.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      a(0,1)=a(0,1)-a(-1,1)
      a(0,n)=a(0,n)-a(-1,1)
      u=RESHAPE((/one,(zero,j=2,n-1),one/),SHAPE(u))
      CALL spline_triluf(a)
      CALL spline_trilus(a,u)
      a(-1,1)=a(-1,1)/(1+a(-1,1)*(u(1,1)+u(n,1)))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_sherman
c-----------------------------------------------------------------------
c     subprogram 17. spline_morrison.
c     uses Sherman-Morrison formula to solve periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_morrison(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: n
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: y
c-----------------------------------------------------------------------
c     solve for x.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      y=x
      CALL spline_trilus(a,y)
      x(1,:)=x(1,:)-a(-1,1)*(y(1,:)+y(n,:))
      x(n,:)=x(n,:)-a(-1,1)*(y(1,:)+y(n,:))
      CALL spline_trilus(a,x)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_morrison
c-----------------------------------------------------------------------
c     subprogram 18. spline_copy.
c     copies one spline_type to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_copy(spl1,spl2)

      TYPE(spline_type), INTENT(IN) :: spl1
      TYPE(spline_type), INTENT(INOUT) :: spl2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ALLOCATED(spl2%xs))CALL spline_dealloc(spl2)
      CALL spline_alloc(spl2,spl1%mx,spl1%nqty)
      spl2%xs=spl1%xs
      spl2%fs=spl1%fs
      spl2%fs1=spl1%fs1
      spl2%name=spl1%name
      spl2%title=spl1%title
      spl2%periodic=spl1%periodic
      spl2%xpower=spl1%xpower
      spl2%x0=spl1%x0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_copy
c-----------------------------------------------------------------------
c     subprogram 19. spline_fill_matrix.
c     fill local matrix into global matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fill_matrix(mat,locmat,istart,jstart,kl,ku)
      INTEGER :: istart,jstart,itot,jtot,i,j,kl,ku,m,n,offset
      REAL(r8), DIMENSION(:,:), ALLOCATABLE,INTENT(INOUT) :: mat
      REAL(r8), DIMENSION(:,:), ALLOCATABLE,INTENT(IN) :: locmat
c-----------------------------------------------------------------------
c     fill matrix.
c-----------------------------------------------------------------------
      itot=SIZE(locmat,1)
      jtot=SIZE(locmat,2)
      offset=kl+ku+1
      DO m=1, itot
         i=m+istart-1
         DO n=1, jtot
            j=n+jstart-1
            mat(offset+i-j,j)=mat(offset+i-j,j)+locmat(m,n)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fill_matrix
c-----------------------------------------------------------------------
c     subprogram 20. spline_get_yp.
c     get yi' with four points for spline boundary condtion .
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_get_yp(x,y,xi,yip,nqty)
      INTEGER, INTENT(IN) :: nqty
      INTEGER :: n,nrhs,lda,info,ldb,i
      INTEGER, DIMENSION(4) :: ipiv
      REAL(r8) :: dx
      REAL(r8), INTENT(IN) :: xi
      REAL(r8), DIMENSION(nqty),INTENT(OUT) :: yip
      REAL(r8), DIMENSION(4), INTENT(IN) :: x
      REAL(r8), DIMENSION(4,nqty), INTENT(IN) :: y
      REAL(r8), DIMENSION(4,nqty) :: b
      REAL(r8), DIMENSION(4,4) :: a
      n=4
      nrhs=nqty
      lda=N
      ldb=N
      a=0
      b=0

      a(1,4)=1
      b(1,:)=y(1,:)
      DO i=2,n
         dx=x(i)-x(1)
         a(i,1)=dx*dx*dx
         a(i,2)=dx*dx
         a(i,3)=dx
         a(i,4)=1
         b(i,:)=y(i,:)
      ENDDO
      CALL dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      dx=xi-x(1)
      yip=(3*b(1,:)*dx+2*b(2,:))*dx+b(3,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_get_yp
c-----------------------------------------------------------------------
c     subprogram 21. spline_thomas.
c     thomas method to solve tri-diagnol matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_thomas(l,d,u,b,n,m)

      INTEGER, INTENT(IN):: n,m
      REAL(r8), DIMENSION(n), INTENT(INOUT):: d
      REAL(r8), DIMENSION(n-1), INTENT(INOUT):: l,u
      REAL(r8), DIMENSION(n,m), INTENT(INOUT):: b
      INTEGER:: i
c-----------------------------------------------------------------------
c     calculate tri-diagno matrix
c     l=[A(1,2),A(2,3),...,A(n-1,n)];
c     d=[A(1,1),A(2,2),...,A(n,n)];
c     u=[A(2,1),A(3,2),...,A(n,n-1)];
c     b is n row m column matrix
c-----------------------------------------------------------------------
      DO i = 2, n
         l(i-1) = l(i-1)/d(i-1)
         d(i) = d(i) - u(i-1) * l(i-1)
         b(i,:) = b(i,:) - b(i-1,:) * l(i-1)
      ENDDO

      b(n,:) = b(n,:) / d(n);
      DO i = n-1, 1, -1
         b(i,:) = (b(i,:) - u(i) * b(i+1,:)) / d(i);
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_thomas
c-----------------------------------------------------------------------
c     subprogram 22. spline_roots.
c     Calculate all the roots of a cubic spline.
c     Necessary because simple iterative techniques may miss cases with
c     multiple roots between knots.
c     Analytics from https://www.e-education.psu.edu/png520/m11_p6.html
c
c     Note: Doesn't handle powers (what are those?)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_roots(spl, iqty, nroots, roots,op_extrap,op_eps)

      TYPE(spline_type), INTENT(INOUT) :: spl
      INTEGER, INTENT(IN):: iqty
      INTEGER, INTENT(OUT) :: nroots
      REAL(r8), DIMENSION(*), INTENT(OUT):: roots ! should be DIMENSION(1:spl%mx*3)
      LOGICAL, INTENT(IN), OPTIONAL :: op_extrap
      REAL(r8), INTENT(IN), OPTIONAL :: op_eps

      LOGICAL, PARAMETER :: debug = .FALSE.

      LOGICAL :: extrap
      INTEGER:: ix, jroot, lx, nzvalid
      REAL(r8)::  f_0, f1_0, f_1, f1_1,
     $  c3, c2, c1, c0, a, b, c, q, r, m, s, t, theta,
     $  z1, z2, z3, x1, x2, x3, last1, last2, last3,
     $  dx, eps, tol
      INTEGER, DIMENSION(:), ALLOCATABLE :: index
      REAL(r8), DIMENSION(:), ALLOCATABLE :: tmproots

      ! allow optional accuracy manipulation repeated solutions
      ! within eps*stepsize are rejected
      IF(PRESENT(op_eps))THEN
         eps = op_eps
      ELSE
         eps = 1e-6
      ENDIF

      ! allow optional extrapolation beyond ends of the spline
      IF(PRESENT(op_extrap))THEN
         extrap = op_extrap
      ELSE
         extrap = .FALSE.
      ENDIF

      roots(1:spl%mx*3) = spl%xs(0) - HUGE(last1)
      last1 = spl%xs(0) - HUGE(last1)
      last2 = spl%xs(0) - HUGE(last2)
      last3 = spl%xs(0) - HUGE(last3)
      lx = spl%mx-1
      jroot = 0
      ! find the analytic roots of a cubic polynomial between each knot
      DO ix=0, lx
         nzvalid = 0
         ! step size and associated tolerance for repeated roots
         dx=spl%xs(ix+1)-spl%xs(ix)
         tol = eps * dx
         ! get the value and derivatives for this step span
         f_0 = spl%fs(ix,iqty)
         f1_0 = spl%fs1(ix,iqty)
         f_1 = spl%fs(ix+1,iqty)
         f1_1 = spl%fs1(ix+1,iqty)
         ! starting from the equation for f(x) in spline_eval,
         ! we reorganize into a classic c3 z^3 + c2 z^2 + c1 z + c0
         c3 = -2*f_0 - 2*f_1 +   dx*f1_0 + dx*f1_1
         c2 =  5*f_0 + 3*f_1 - 2*dx*f1_0 - dx*f1_1
         c1 = -4*f_0 + dx*f1_0
         c0 = f_0
         c3 =  2*f_0 - 2*f_1 +   dx*f1_0 - dx*f1_1
         c2 = -3*f_0 + 3*f_1 - 2*dx*f1_0 + dx*f1_1
         c1 = dx*f1_0
         c0 = f_0
         c3 =  2*f_0 - 2*f_1 +   dx*f1_0 + dx*f1_1
         c2 = -3*f_0 + 3*f_1 - 2*dx*f1_0 - dx*f1_1
         c1 = dx*f1_0
         c0 = f_0
         ! we have to watch out for secret reductions!
         IF(c3==0)THEN
            IF(debug) PRINT *, "  >> Spline quadratic @",spl%xs(ix)
            a = c2
            b = c1
            c = c0
            IF(b**2 - 4.0*a*c >= 0)THEN
               z1 = (-b + sqrt(b**2 - 4.0*a*c)) / (2.0 * a)
               z2 = (-b - sqrt(b**2 - 4.0*a*c)) / (2.0 * a)
               nzvalid = 2
            ELSEIF(a==0)THEN
               IF(debug) PRINT *, "  >> Spline linear @",spl%xs(ix)
               z1 = -c / b
               z2 = -huge(z2)
               nzvalid = 1
            ELSEIF(a==0 .and. b==0 .and. c==0)THEN
               IF(debug) PRINT *, "  >> Spline all 0 @",spl%xs(ix)
               z1 = 0
               nzvalid = 1
            ELSE
               z1 = -huge(z1)
               z2 = -huge(z2)
               nzvalid = 0
            ENDIF
           z3 = -huge(z3)
         ELSE
            ! truely a cubic case
            ! since we want f = 0 we can normalize out the c3
            ! so 0 = z^3 + a z^2 + b z + c
            a = c2 / c3
            b = c1 / c3
            c = c0 / c3
            ! now we just get the analytic solutions
            q = (a * a - 3.0 * b) / 9.0
            r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0
            m = r**2 - q**3  ! discrimenent
            IF(m > 0)THEN ! one real root
               s = -sign(1.0_r8, r) * (abs(r) + sqrt(m))**(1.0/3.0)
               t = q / s
               z1 = s + t - (a / 3.0)
               z2 = -huge(z2)
               z3 = -huge(z3)
               nzvalid = 1
            ELSE  ! three real roots (watch out for repeates)
               theta = acos(r / sqrt(q**3))
               z1 = -2*sqrt(q)*cos(theta/3.0) - a/3.0
               z2 = -2*sqrt(q)*cos((theta+twopi)/3.0) - a/3.0
               z3 = -2*sqrt(q)*cos((theta-twopi)/3.0) - a/3.0
               nzvalid = 3
               IF(abs(z2-z3)<tol) nzvalid = 2  ! double root
            ENDIF
         ENDIF
         ! note spline_eval uses step-size normalized coordinate z,
         ! but we want real x values
         x1 = spl%xs(ix) + z1 * dx
         x2 = spl%xs(ix) + z2 * dx
         x3 = spl%xs(ix) + z3 * dx
         ! check if they are in the knot interval (exclude repeats if perioic)
         IF(debug .AND. (f_0*f_1 <= 0))THEN
            PRINT *," > f zero crossing between",spl%xs(ix),spl%xs(ix+1)
            PRINT *,"  >> f_0, f_1 =",f_0,f_1
            PRINT *,"  >> nzvalid =",nzvalid
            PRINT *,"  >> z1, z2, z3 = ",z1,z2,z3
            PRINT *,"  >> c3,c2,c1,c0 =",c3,c2,c1,c0
         ENDIF
         IF(nzvalid>0 .AND. (z1>-eps .OR. (ix==0 .AND. extrap)))THEN
            IF(z1<=1+eps .OR. (ix==lx .AND. extrap)) THEN
               IF(debug) PRINT '(a12,es17.8e3,a4,es17.8e3,a15,'//
     $            'I3,a1,I3,a11,es13.4e3,a1,es13.4)',
     $            "  > Found z1",z1,", x1",x1," between knots ",ix,",",
     $            ix+1," where x =",spl%xs(ix),",",spl%xs(ix+1)
               ! refine solution numerically (analytics vs reality of interp)
               CALL spline_refine_root(spl,iqty,x1)
               ! avoid repeated left/right solutions at a f=0 knot
               IF(abs(last1-x1)>tol .AND. abs(last2-x1)>tol .AND.
     $            abs(last3-x1)>tol) THEN
                  jroot = jroot + 1
                  roots(jroot) = x1
               ENDIF
            ENDIF
         ENDIF
         IF(nzvalid>1 .AND. (z2>-eps .OR. (ix==0 .AND. extrap)))THEN
            IF(z2<=1+eps .OR. (ix==lx .AND. extrap)) THEN
               IF(debug) PRINT '(a12,es17.8e3,a4,es17.8e3,a15,'//
     $            'I3,a1,I3,a11,es13.4e3,a1,es13.4)',
     $            "  > Found z2",z2,", x2",x2," between knots ",ix,
     $            ",",ix+1," where x =",spl%xs(ix),",",spl%xs(ix+1)
               ! refine solution numerically (analytics vs reality of interp)
               CALL spline_refine_root(spl,iqty,x2)
               ! avoid double roots or repeats right at a knot location
               IF(abs(last1-x2)>tol .AND. abs(last2-x2)>tol .AND.
     $            abs(last3-x2)>tol .AND. abs(x1-x2)>tol) THEN
                  jroot = jroot + 1
                  roots(jroot) = x2
               ENDIF
            ENDIF
         ENDIF
         IF(nzvalid>2 .AND. (z3>-eps .OR. (ix==0 .AND. extrap)))THEN
            IF(z3<=1+eps .OR. (ix==lx .AND. extrap)) THEN
               IF(debug) PRINT '(a12,es17.8e3,a4,es17.8e3,a15,'//
     $            'I3,a1,I3,a11,es13.4e3,a1,es13.4)',
     $            "  > Found z3",z3,", x3",x3," between knots ",ix,
     $            ",",ix+1," where x =",spl%xs(ix),",",spl%xs(ix+1)
               ! refine solution numerically (analytics vs reality of interp)
               CALL spline_refine_root(spl,iqty,x3)
               ! avoid double roots or repeats right at a knot location
               IF(abs(last1-x3)>tol .AND. abs(last2-x3)>tol .AND.
     $            abs(last3-x3)>tol .AND. abs(x1-x3)>tol .AND.
     $            abs(x2-x3)>tol) THEN
                  jroot = jroot + 1
                  roots(jroot) = x3
               ENDIF
            ENDIF
         ENDIF
         last1 = roots(max(1,jroot))
         last2 = roots(max(1,jroot-1))
         last3 = roots(max(1,jroot-2))
      ENDDO
      nroots = jroot
      ! sort the roots lowest to highest
      ALLOCATE(index(nroots), tmproots(nroots))
      index=(/(ix,ix=1,nroots)/)
      tmproots(1:nroots) = roots(1:nroots)
      CALL bubble(tmproots,index,1,nroots)
      DO ix=1,nroots
         tmproots(ix) = roots(index(nroots + 1 - ix))
      ENDDO
      roots(1:nroots) = tmproots(1:nroots)
      IF(nroots>0 .AND. debug)
     $   PRINT *,"  > Sorted roots are",roots(1:nroots)
      DEALLOCATE(index, tmproots)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_roots

c-----------------------------------------------------------------------
c     subprogram 23. spline_refine_root.
c     Calculate the root of a spline using newton iteration from an
c     initial guess.
c     Doesn't handle powers (what are those?)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_refine_root(spl, iqty, x, op_tol)

      TYPE(spline_type), INTENT(INOUT) :: spl
      INTEGER, INTENT(IN):: iqty
      REAL(r8), INTENT(INOUT) :: x
      REAL(r8), INTENT(IN), OPTIONAL :: op_tol
      ! declare variables
      INTEGER :: iroot,it,nroots
      INTEGER, PARAMETER :: itmax=500
      REAL(r8) :: tol=1e-12
      REAL(r8) :: dx,lx,lf,f,df

      ! set optional tolerance
      tol = 1e-12
      IF(PRESENT(op_tol)) tol = op_tol

      ! if its exact, we don't need to do anything
      CALL spline_eval(spl,x,1)
      IF(spl%f(iqty) == 0) RETURN
      ! otherwise, we'll iterate
      lx=spl%xs(spl%ix+1)-spl%xs(spl%ix)  ! note ix set in above eval
      lf = maxval(spl%fs(:,iqty))-minval(spl%fs(:,iqty))
      dx=lx
      f=huge(f)
      it=0
      DO
         CALL spline_eval(spl,x,1)
         df=spl%f(iqty)-f
         !IF(abs(dx) < tol*lx .OR. abs(df) < tol*lf .OR. it >= itmax)EXIT
         IF(abs(dx) <= abs(tol*lx) .OR. it >= itmax)EXIT
         it=it+1
         f=spl%f(iqty)
         dx=-spl%f(iqty)/spl%f1(iqty)
         x=x+dx
      ENDDO
      ! abort on failure.
      IF(it >= itmax)THEN
         WRITE(*,*)
         WRITE(*,*) "!! WARNING: root refining convergence failure"
         WRITE(*,*) " -> estimated root ",x,
     $              " has error/tol ",abs(dx)/(tol*lx),abs(df)/(tol*lf)
         WRITE(*,*)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_refine_root

      END MODULE spline_mod
