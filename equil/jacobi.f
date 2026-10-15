c-----------------------------------------------------------------------
c     file jacobi.f.
c     computes jacobi polynomials and related functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. jacobi_mod.
c     1. jacobi_alloc.
c     2. jacobi_dealloc.
c     3. jacobi_pzero.
c     4. jacobi_qzero.
c     5. jacobi_rzero.
c     6. jacobi_eval.
c     7. jacobi_basis.
c     8. jacobi_basis_modal.
c     9. jacobi_basis_nodal.
c     10. jacobi_interp.
c     11. jacobi_product.
c-----------------------------------------------------------------------
c     subprogram 0. jacobi_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE jacobi_mod
      USE local_mod
      IMPLICIT NONE

      TYPE :: jacobi_type
      CHARACTER(4) :: quadr
      LOGICAL :: nodal,uniform,massoc
      INTEGER :: np,order
      REAL(r8) :: alpha,beta
      REAL(r8), DIMENSION(:), POINTER :: p,q,q2,q3,pb,qb,qb2,qb3,
     $     pzero,qzero,rzero,qweight,pweight,rweight,weight,node
      REAL(r8), DIMENSION(:,:), POINTER :: mass,stiff
      END TYPE jacobi_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. jacobi_alloc.
c     allocates jacobi_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_alloc(jac,np,nodal,uniform,quadr)

      TYPE(jacobi_type) :: jac
      INTEGER, INTENT(IN) :: np
      LOGICAL, INTENT(IN) :: nodal,uniform
      CHARACTER(*), INTENT(IN) :: quadr
c-----------------------------------------------------------------------
c     allocate jacobi_type.
c-----------------------------------------------------------------------
      jac%np=np
      jac%nodal=nodal
      jac%uniform=uniform
      jac%quadr=quadr
      jac%massoc=.FALSE.
      IF(nodal)THEN
         jac%alpha=0
      ELSE
         jac%alpha=1
      ENDIF
      jac%beta=jac%alpha
      jac%order=2
      ALLOCATE(jac%p(0:np+1))
      ALLOCATE(jac%q(0:np+1))
      ALLOCATE(jac%q2(0:np+1))
      ALLOCATE(jac%q3(0:np+1))
      ALLOCATE(jac%pb(0:np))
      ALLOCATE(jac%qb(0:np))
      ALLOCATE(jac%qb2(0:np))
      ALLOCATE(jac%qb3(0:np))
      ALLOCATE(jac%pzero(0:np))
      ALLOCATE(jac%qzero(0:np))
      ALLOCATE(jac%rzero(0:np))
      ALLOCATE(jac%qweight(0:np))
      ALLOCATE(jac%pweight(0:np))
      ALLOCATE(jac%rweight(0:np))
      ALLOCATE(jac%node(0:np))
      ALLOCATE(jac%weight(0:np))
c-----------------------------------------------------------------------
c     evaluate weights and nodes.
c-----------------------------------------------------------------------
      CALL jacobi_pzero(jac)
      CALL jacobi_qzero(jac)
      CALL jacobi_rzero(jac)
c-----------------------------------------------------------------------
c     choose weights and nodes.
c-----------------------------------------------------------------------
      SELECT CASE(jac%quadr)
      CASE("gll")
         jac%node=jac%qzero
         jac%weight=jac%qweight
      CASE("gl0")
         jac%node=jac%pzero
         jac%weight=jac%pweight
      CASE("none")
         jac%node=0
         jac%weight=0
      CASE DEFAULT
         CALL program_stop("Jacobi_alloc: cannot recognize quad_type = "
     $        //TRIM(jac%quadr)//".")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_alloc
c-----------------------------------------------------------------------
c     subprogram 2. jacobi_dealloc.
c     deallocates jacobi_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_dealloc(jac)

      TYPE(jacobi_type) :: jac
c-----------------------------------------------------------------------
c     deallocate jacobi_type.
c-----------------------------------------------------------------------
      DEALLOCATE(jac%p,jac%q,jac%q2,jac%q3,jac%pb,jac%qb,jac%qb2,
     $     jac%qb3,jac%pzero,jac%qzero,jac%rzero,
     $     jac%qweight,jac%pweight,jac%rweight,jac%node,jac%weight)
      IF(jac%massoc)DEALLOCATE(jac%mass,jac%stiff)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. jacobi_pzero.
c     finds zeros of jacobi polynomials.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_pzero(jac)

      TYPE(jacobi_type) :: jac

      INTEGER :: j,nh,it,itmax=20,n
      REAL(r8) :: dth,pder,poly,recsum,x,dx,eps=1e-14
      REAL(r8), DIMENSION(jac%np/2+1) :: theta
c-----------------------------------------------------------------------
c     set up recursion relation for initial guess for the roots.
c-----------------------------------------------------------------------
      nh=jac%np/2
      dth=pi/(2*jac%np+3)
      theta=(/(2*n-1,n=1,nh+1)/)*dth
c-----------------------------------------------------------------------
c     compute first half of roots by polynomial deflation.
c-----------------------------------------------------------------------
      DO j=0,nh
         it=0
         x=COS(theta(j+1))
         DO
            it=it+1
            CALL jacobi_eval(x,jac)
            poly=jac%p(jac%np+1)
            pder=jac%q(jac%np+1)
            IF(j==0)THEN
               recsum=0
            ELSE
               recsum=SUM(1/(x-jac%pzero(0:j-1)))
            ENDIF
            dx=-poly/(pder-recsum*poly)
            x=x+dx
            IF(ABS(dx) < eps .OR. it >= itmax)EXIT
         ENDDO
         jac%pzero(j)=x
         jac%pweight(j)=2/((1-x**2)*pder**2)
      ENDDO
c-----------------------------------------------------------------------
c     use symmetry for second half of roots and reverse order.
c-----------------------------------------------------------------------
      DO j=0,nh
         jac%pzero(jac%np-j)=-jac%pzero(j)
         jac%pweight(jac%np-j)=jac%pweight(j)
      ENDDO
      jac%pzero=jac%pzero(jac%np:0:-1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_pzero
c-----------------------------------------------------------------------
c     subprogram 4. jacobi_qzero.
c     finds zeros of jacobi polynomial derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_qzero(jac)

      TYPE(jacobi_type) :: jac

      INTEGER :: j,nh,it,itmax=20,n
      REAL(r8) :: a,b,c,dth,pder,poly,recsum,x,dx,eps=1e-14
      REAL(r8), DIMENSION((jac%np-1)/2) :: theta
c-----------------------------------------------------------------------
c     uniform spacing.
c-----------------------------------------------------------------------
      IF(jac%uniform)THEN
         jac%qzero=(/(j,j=0,jac%np)/)*2/REAL(jac%np,r8)-1
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     compute coefficients of polynomial whose roots are desired.
c-----------------------------------------------------------------------
      c=2*jac%np+jac%alpha+jac%beta
      a=jac%np*(jac%alpha-jac%beta)/c
      b=2*(jac%np+jac%alpha)*(jac%np+jac%beta)/c
c-----------------------------------------------------------------------
c     set up recursion relation for initial guess for the roots.
c-----------------------------------------------------------------------
      nh=(jac%np+1)/2
      dth=pi/(2*jac%np+1)
      theta=(/(1+2*n,n=1,nh-1)/)*dth
c-----------------------------------------------------------------------
c     compute first half of roots by polynomial deflation.
c-----------------------------------------------------------------------
      DO j=1,nh-1
         x=COS(theta(j))
         it=0
         DO
            it=it+1
            CALL jacobi_eval(x,jac)
            poly=((a-jac%np*x)*jac%p(jac%np)
     $           +b*jac%p(jac%np-1))/(1-x**2)
            pder=((a-jac%np*x)*jac%q(jac%np)
     $           +b*jac%q(jac%np-1)
     $           -jac%np*jac%p(jac%np)+2*x*poly)/(1-x**2)
            IF(j==1)THEN
               recsum=0
            ELSE
               recsum=SUM(1/(x-jac%qzero(1:j-1)))
            ENDIF
            dx=-poly/(pder-recsum*poly)
            x=x+dx
            IF(ABS(dx) < eps .OR. it >= itmax)EXIT
         ENDDO
         jac%qzero(j)=x
         CALL jacobi_eval(x,jac)
         jac%qweight(j)=2/(jac%np*(jac%np+1)*jac%p(jac%np)**2)
      ENDDO
c-----------------------------------------------------------------------
c     use symmetry for second half of roots.
c-----------------------------------------------------------------------
      DO j=1,nh-1
         jac%qzero(jac%np-j)=-jac%qzero(j)
         jac%qweight(jac%np-j)=jac%qweight(j)
      ENDDO
      IF(mod(jac%np,2) == 0)THEN
         jac%qzero(nh)=0
         CALL jacobi_eval(0._r8,jac)
         jac%qweight(nh)=2/(jac%np*(jac%np+1)*jac%p(jac%np)**2)
      ENDIF
c-----------------------------------------------------------------------
c     end points values and reverse order.
c-----------------------------------------------------------------------
      jac%qzero(0)=1
      jac%qzero(jac%np)=-1
      IF(jac%quadr=="gll")THEN
         jac%qweight(0)=2.0_r8/(jac%np*(jac%np+1))
         jac%qweight(jac%np)=2.0_r8/(jac%np*(jac%np+1))
      ENDIF
      jac%qzero=jac%qzero(jac%np:0:-1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_qzero
c-----------------------------------------------------------------------
c     subprogram 5. jacobi_rzero.
c     finds zeros of (P_(n-1)+P_(n))/(x+1).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_rzero(jac)

      TYPE(jacobi_type) :: jac

      INTEGER :: j,it,itmax=20,n
      REAL(r8) :: dth,pder,poly,recsum,x,dx,eps=1e-14
      REAL(r8), DIMENSION(jac%np) :: theta
c-----------------------------------------------------------------------
c     set up recursion relation for initial guess for the roots.
c-----------------------------------------------------------------------
      dth=pi/(2*jac%np+3)
      theta=(/(1+2*n,n=1,jac%np)/)*dth
c-----------------------------------------------------------------------
c     compute first half of roots by polynomial deflation.
c-----------------------------------------------------------------------
      DO j=1,jac%np
         it=0
         x=COS(theta(j))
         DO
            it=it+1
            CALL jacobi_eval(x,jac)
            poly=(jac%p(jac%np)+jac%p(jac%np+1))/(x+1)
            pder=(jac%q(jac%np)+jac%q(jac%np+1))/(x+1)
     $           -(jac%p(jac%np)+jac%p(jac%np+1))/(x+1)**2
            IF(j==0)THEN
               recsum=0
            ELSE
               recsum=SUM(1/(x-jac%rzero(1:j-1)))
            ENDIF
            dx=-poly/(pder-recsum*poly)
            x=x+dx
            IF(ABS(dx) < eps .OR. it >= itmax)EXIT
         ENDDO
         jac%rzero(j)=x
         CALL jacobi_eval(x,jac)
         jac%rweight(j)=(1-x)/((jac%np+1)*jac%p(jac%np))**2
      ENDDO
c-----------------------------------------------------------------------
c     end point value.
c-----------------------------------------------------------------------
      jac%rzero(0)=-1
      jac%rweight(0)=2.0_r8/(jac%np+1)**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_rzero
c-----------------------------------------------------------------------
c     subprogram 6. jacobi_eval.
c     computes jacobi polynomials and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_eval(x,jac)

      REAL(r8), INTENT(IN) :: x
      TYPE(jacobi_type) :: jac

      INTEGER :: n
      REAL(r8) :: ab1,ab2,a1,a2,a3,a4,b3
c-----------------------------------------------------------------------
c     compute jacobi polynomials and their derivatives.
c-----------------------------------------------------------------------
      jac%p(0)=1
      jac%q(0)=0
      jac%q2(0)=0
      jac%q3(0)=0
      IF(jac%np > 0)THEN
         jac%p(1)=((jac%alpha-jac%beta)+(2+jac%alpha+jac%beta)*x)/2
         jac%q(1)=1+(jac%alpha+jac%beta)/2
         jac%q2(1)=0
         jac%q3(1)=0
         ab1=jac%alpha+jac%beta
         ab2=jac%alpha**2-jac%beta**2
         DO n=1,jac%np
            a1=2*(n+1)*(n+ab1+1)*(2*n+ab1)
            a2=(2*n+ab1+1)*ab2
            b3=2*n+ab1
            a3=b3*(b3+1)*(b3+2)
            a4=2*(n+jac%alpha)*(n+jac%beta)*(2*n+ab1+2)
            jac%p(n+1)=((a2+a3*x)*jac%p(n)-a4*jac%p(n-1))/a1
            jac%q(n+1)=((a2+a3*x)*jac%q(n)
     $           -a4*jac%q(n-1)+a3*jac%p(n))/a1
            IF(jac%order < 2)CYCLE
            jac%q2(n+1)=((a2+a3*x)*jac%q2(n)
     $           -a4*jac%q2(n-1)+2*a3*jac%q(n))/a1
            IF(jac%order < 3)CYCLE
            jac%q3(n+1)=((a2+a3*x)*jac%q3(n)
     $           -a4*jac%q3(n-1)+3*a3*jac%q2(n))/a1
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_eval
c-----------------------------------------------------------------------
c     subprogram 7. jacobi_basis.
c     computes jacobi basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_basis(x,jac)

      REAL(r8), INTENT(IN) :: x
      TYPE(jacobi_type) :: jac
c-----------------------------------------------------------------------
c     compute nodal basis.
c-----------------------------------------------------------------------
      IF(jac%nodal)THEN
         CALL jacobi_basis_nodal(x,jac)
      ELSE
         CALL jacobi_basis_modal(x,jac)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_basis
c-----------------------------------------------------------------------
c     subprogram 8. jacobi_basis_modal.
c     computes jacobi modal basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_basis_modal(x,jac)

      REAL(r8), INTENT(IN) :: x
      TYPE(jacobi_type) :: jac
c-----------------------------------------------------------------------
c     compute modal basis.
c-----------------------------------------------------------------------
      CALL jacobi_eval(x,jac)
      jac%pb(0)=(1-x)/2
      jac%pb(1:jac%np-1)=(1-x**2)*jac%p(0:jac%np-2)
      jac%pb(jac%np)=(1+x)/2
      jac%qb(0)=-.5
      jac%qb(1:jac%np-1)=(1-x**2)*jac%q(0:jac%np-2)
     $     -2*x*jac%p(0:jac%np-2)
      jac%qb(jac%np)=.5
      IF(jac%order >= 2)THEN
         jac%qb2(0)=0
         jac%qb2(1:jac%np-1)=(1-x**2)*jac%q2(0:jac%np-2)
     $        -4*x*jac%q(0:jac%np-2)-2*jac%p(0:jac%np-2)
         jac%qb2(jac%np)=0
      ENDIF
      IF(jac%order >= 3)THEN
         jac%qb3(0)=0
         jac%qb3(1:jac%np-1)=(1-x**2)*jac%q3(0:jac%np-2)
     $        -6*x*jac%q2(0:jac%np-2)-6*jac%q(0:jac%np-2)
         jac%qb3(jac%np)=0
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_basis_modal
c-----------------------------------------------------------------------
c     subprogram 9. jacobi_basis_nodal.
c     computes jacobi nodal basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_basis_nodal(x,jac)

      REAL(r8), INTENT(IN) :: x
      TYPE(jacobi_type) :: jac

      LOGICAL, DIMENSION(0:jac%np) :: mask
      INTEGER :: l,m,n
      REAL(r8) :: fac
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      mask=.TRUE.
      jac%qb=0
      jac%qb2=0
c-----------------------------------------------------------------------
c     start loop over nodes and compute basis function.
c-----------------------------------------------------------------------
      DO n=0,jac%np
         mask(n)=.FALSE.
         fac=1/PRODUCT(jac%qzero(n)-jac%qzero,1,mask)
         jac%pb(n)=PRODUCT(x-jac%qzero,1,mask)*fac
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         DO m=0,jac%np
            IF(m == n)CYCLE
            mask(m)=.FALSE.
            jac%qb(n)=jac%qb(n)+PRODUCT(x-jac%qzero,1,mask)
            DO l=0,jac%np
               IF(l == m .OR. l == n)CYCLE
               mask(l)=.FALSE.
               jac%qb2(n)=jac%qb2(n)+PRODUCT(x-jac%qzero,1,mask)
               mask(l)=.TRUE.
            ENDDO
            mask(m)=.TRUE.
         ENDDO
         jac%qb(n)=jac%qb(n)*fac
         jac%qb2(n)=jac%qb2(n)*fac
c-----------------------------------------------------------------------
c     finish loop over basis functions.
c-----------------------------------------------------------------------
         mask(n)=.TRUE.
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_basis_nodal
c-----------------------------------------------------------------------
c     subprogram 10. jacobi_interp.
c     find values of jacobi polynomials at interpolatory points.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_interp(nx,x,ix,jac,bqx)

      INTEGER, INTENT(IN) :: nx
      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: ix
      TYPE(jacobi_type) :: jac
      REAL(r8), DIMENSION(0:,:), INTENT(OUT) :: bqx

      INTEGER :: j
      REAL(r8) :: loc
      REAL(r8), DIMENSION(0:nx) :: xx
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      xx = (/(j,j=0,nx)/)
c-----------------------------------------------------------------------
c     locate x interval, compute values of basis functions.
c-----------------------------------------------------------------------
      ix=nx-1
      DO
         IF(x >= xx(ix) .OR. ix == 0)EXIT
         ix=ix-1
      ENDDO
      DO
         IF(x < xx(ix+1) .OR. ix == nx-1)EXIT
         ix=ix+1
      ENDDO

      loc=two*(x-xx(ix))-one
      CALL jacobi_basis(loc,jac)

      bqx(:,1)=jac%pb
      IF(SIZE(bqx,2) > 1)bqx(:,2)=jac%qb*two
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_interp
c-----------------------------------------------------------------------
c     subprogram 11. jacobi_product.
c     computes mass and stiffness matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE jacobi_product(basis,quad)

      TYPE(jacobi_type) :: basis,quad

      INTEGER :: m,n,ix,np
      REAL(r8), PARAMETER :: tol=1e-15
      REAL(r8) :: x,w
c-----------------------------------------------------------------------
c     compute scalar products, using Gaussian quadrature.
c-----------------------------------------------------------------------
      np=basis%np
      IF(.NOT. basis%massoc)THEN
         ALLOCATE(basis%mass(0:np,0:np),basis%stiff(0:np,0:np))
         basis%massoc=.TRUE.
      ENDIF
      basis%mass=0
      basis%stiff=0
      DO ix=0,quad%np
         w=quad%weight(ix)
         x=quad%node(ix)
         CALL jacobi_eval(x,basis)
         CALL jacobi_basis(x,basis)
         DO m=0,np
            DO n=m,np
               basis%mass(m,n)=basis%mass(m,n)
     $              +w*basis%pb(m)*basis%pb(n)
               basis%stiff(m,n)=basis%stiff(m,n)
     $              +w*basis%qb(m)*basis%qb(n)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     symmetrize.
c-----------------------------------------------------------------------
      DO m=0,np
         DO n=0,m-1
            basis%mass(m,n)=basis%mass(n,m)
            basis%stiff(m,n)=basis%stiff(n,m)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     truncate small terms.
c-----------------------------------------------------------------------
      WHERE(ABS(basis%mass) < tol)
         basis%mass = 0
      ENDWHERE
      WHERE(ABS(basis%stiff) < tol)
         basis%stiff = 0
      ENDWHERE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jacobi_product
      END MODULE jacobi_mod
