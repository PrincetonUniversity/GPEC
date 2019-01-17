c-----------------------------------------------------------------------
c     program deltac.f.
c     Galerkin method for the inner layer model in configuration space.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. delmatch_mod.
c     1. deltac_run.
c     2. deltac_init_asymp.
c     3. deltac_get_ua.
c     4. deltac_get_dua.
c     5. deltac_solve.
c     6. deltac_alloc.
c     7. deltac_dealloc.
c     8. deltac_make_grid.
c     9. deltac_pack.
c     10. deltac_hermite.
c     11. deltac_make_map_hermite.
c     12. deltac_make_arrays.
c     13. deltac_gauss_quad.
c     14. deltac_get_uv.
c     15. deltac_assemble_mat.
c     16. deltac_assemble_rhs.
c     17. deltac_set_boundary.
c     18. deltac_extension.
c     19. deltac_lsode_int.
c     20. deltac_lsode_der.
c     21. deltac_get_d2ua.
c     22. deltac_ua_diagnose.
c     23. deltac_get_solution.
c     24. deltac_output_solution.
c     25. deltac_read_parameters.
c     26. deltac_estimate_zi.
c-----------------------------------------------------------------------
c     subprogram 0. delmatch_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE deltac_mod
      USE gamma_mod
      USE jacobi_mod
      USE deltar_mod, ONLY : resist_type
      IMPLICIT NONE

      TYPE :: inner_type
      INTEGER :: ising
      REAL(r8) :: e,f,g,h,k,m,taua,taur,v1
      REAL(r8) :: di,dr,p1,sfac
      COMPLEX(r8) :: eig,q,x0
      COMPLEX(r8),DIMENSION(2) :: deltar
      END TYPE inner_type
      
      TYPE :: asymp_type
      REAL(r8), DIMENSION(2) :: p
      COMPLEX(r8), DIMENSION(2):: s
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: v
      END TYPE asymp_type
      
      TYPE :: hermite_type
      REAL(r8), DIMENSION(0:3) :: pb,qb
      END TYPE hermite_type

      TYPE :: cell_type
      CHARACTER(6) :: etype
      INTEGER :: np
      INTEGER, DIMENSION(:), ALLOCATABLE :: emap
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: map
      REAL(r8) :: x_lsode
      REAL(r8), DIMENSION(2) :: x 
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rhs
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat
      END TYPE cell_type

      TYPE :: interval_type
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,dx
      TYPE(cell_type), DIMENSION(:), POINTER :: cell
      END TYPE interval_type
      
      TYPE :: solution_type
      REAL(r8),DIMENSION(:), ALLOCATABLE :: xvar
      COMPLEX(r8), DIMENSION(:,:,:),ALLOCATABLE :: sol
      END TYPE solution_type
      
      TYPE :: gal_type
      INTEGER :: nx,nq,ndim,kl,ku,ldab,msol
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      REAL(r8) :: pfac,dx1,dx2
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rhs,sol
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: mat
      TYPE(interval_type), POINTER :: intvl
      TYPE(jacobi_type) :: quad
      TYPE(hermite_type) :: hermite
      END TYPE gal_type

      LOGICAL, PRIVATE :: diagnose_res=.FALSE., method=.true.
      lOGICAL :: deltac_bin_sol=.false.,deltac_out_sol=.false.
      LOGICAL :: dx1dx2_flag=.false.,rescale=.true.
      LOGICAL :: restore_uh=.false.,output_sol=.false.
      LOGICAL :: restore_us=.false.,restore_ul=.false.
      LOGICAL :: noexp=.false.
      LOGICAL :: diagnose_params=.FALSE.
      CHARACTER(10) :: gal_method="normal",side="right"
      CHARACTER(256) :: deltabin_filename,galsol_filename
      CHARACTER(256) :: galsol_filename_cut
      INTEGER :: interp_np=10,xmax_method
      INTEGER :: basis_type=0
      INTEGER, PRIVATE :: np=3,mpert=3,msol=2
      INTEGER :: order_pow=10,order_exp=3
      INTEGER :: nx=128,nq=4,fulldomain=0
      INTEGER :: cutoff=5,outt=3,nx_ua=100
      INTEGER, DIMENSION(4), PRIVATE:: tid=(/3,5,6,4/)
      INTEGER :: deltac_bin_unit=120
      INTEGER :: deltac_out_unit=121
      REAL(r8) :: dx1=1e-2,dx2=1e-2
      REAL(r8) :: xmin=0,deltac_tol=1e-5,pfac=1
      REAL(r8) :: xmax=1,xfac=1
      REAL(r8) :: x0_ua=0.01,x1_ua=1
      REAL(r8) :: zi_deltac
      COMPLEX(r8) :: q_deltac
      TYPE(inner_type), PRIVATE :: in
      TYPE(asymp_type), PRIVATE :: asp
      TYPE(cell_type), POINTER :: cell  
      TYPE(solution_type), POINTER :: sol

      CONTAINS    
c-----------------------------------------------------------------------
c     subprogram 1. deltac_run.
c     sets up and solve inner layer.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_run(restype,eig,delta,fulldelta)

      TYPE(resist_type), INTENT(IN) :: restype
      COMPLEX(r8), INTENT(IN) :: eig
      COMPLEX(r8), DIMENSION(2), INTENT(OUT) :: delta
      COMPLEX(r8), DIMENSION(2,2), INTENT(OUT) :: fulldelta
      CHARACTER(80) :: message
      REAL(r8) :: x0, q0,tmp1,tmp2,tmp3,tmp4
      COMPLEX(r8) :: tmp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(6x,"E",10x,"F",10x,"H",10x,"G",10x,"K",10x,"M"//1p,6e11.3/)
 20   FORMAT(6x,"DI",9x,"DR",9x,"p1",8x,"taua",7x,"taur",7x,"sfac"//
     $     1p,6e11.3/)
 30   FORMAT(6x,"X0",9x,"Q0",9X,"V1"//1p,3e11.3)
c-----------------------------------------------------------------------
c     set the domain to be solved
c-----------------------------------------------------------------------
      SELECT CASE (fulldomain)
      CASE(0)
         xmin=0
      CASE(1,2)
         xmin=-xmax
      END SELECT
c-----------------------------------------------------------------------
c     copy input values.
c-----------------------------------------------------------------------
      in%e=restype%e
      in%f=restype%f
      in%g=restype%g
      in%h=restype%h
      in%k=restype%k
      in%m=restype%m
      in%taua=restype%taua
      in%taur=restype%taur
      in%v1=restype%v1
      in%ising=restype%ising
      in%eig=eig
      
      in%dr=in%e+in%f+in%h*in%h
      in%di=in%e+in%f+in%h-0.25
      in%p1=SQRT(-in%di)
c-----------------------------------------------------------------------
c     define scale factors.
c-----------------------------------------------------------------------
      in%sfac=in%taur/in%taua
      x0=in%sfac**(-1._r8/3._r8)
      q0=x0/in%taua
      in%q=in%eig/q0
      in%x0=x0
c-----------------------------------------------------------------------
c     write output file with inner region parameters.
c-----------------------------------------------------------------------
      IF(diagnose_params)THEN
         OPEN(UNIT=debug_unit,FILE="params.out",STATUS="REPLACE")
         WRITE(debug_unit,'(2a/)')
     $        " deltabin_filename = ",TRIM(deltabin_filename)
         WRITE(debug_unit,10)in%e,in%f,in%h,in%g,in%k,in%m
         WRITE(debug_unit,20)in%di,in%dr,in%p1,in%taua,in%taur,in%sfac
         WRITE(debug_unit,30)x0,q0,in%v1
         CLOSE(UNIT=debug_unit)
      endif
c-----------------------------------------------------------------------
c     setup asymptotic solutions at large x.
c-----------------------------------------------------------------------
      CALL deltac_init_asymp
      SELECT CASE(xmax_method)
      CASE(1)
         tmp1=SQRT(ABS(in%q)**.5*(ABS(asp%s(1))+1))*xfac
         tmp2=SQRT(ABS(in%q)**.5*(ABS(asp%s(2))+1))*xfac
         IF (tmp2 > tmp1) tmp1=tmp2
         xmax=tmp2
      CASE(2)
         tmp1=MAXVAL(ABS(asp%v(:,5,1)))
         tmp2=MAXVAL(ABS(asp%v(:,5,2)))
         tmp3=tmp1/tmp2
         tmp1=MAXVAL(ABS(asp%v(:,6,1)))
         tmp2=MAXVAL(ABS(asp%v(:,6,2)))
         tmp4=tmp1/tmp2
         IF (tmp4 > tmp3) tmp3=tmp4
         xmax=tmp3*xfac
      CASE(3)
         xmax=ABS(in%q)**.25*xfac
      CASE(4)
         xmax=xfac*ABS(in%q)**.25
     $        *MAXVAL(1+SQRT(ABS(asp%s(:))))
      CASE DEFAULT
         WRITE(message,'(a,i2)')
     $        "deltac_run: invalid value xmax_method = ",xmax_method
         CALL program_stop(message)
      END SELECT
c-----------------------------------------------------------------------
c     set the domain to be solved
c-----------------------------------------------------------------------
      SELECT CASE (fulldomain)
      CASE(0)
         xmin=0
      CASE(1,2)
         xmin=-xmax
      CASE DEFAULT
         WRITE(message,'(a,i2)') 
     $        "deltac_run: invalide value fulldomain = ",fulldomain
      END SELECT
      IF(diagnose_res)CALL deltac_ua_diagnose
c-----------------------------------------------------------------------
c     estimate zi.
c-----------------------------------------------------------------------
      CALL deltac_estimate_zi(zi_deltac)
      q_deltac=in%q
c-----------------------------------------------------------------------
c     run galerkin method to solve the inner layer.
c-----------------------------------------------------------------------
      CALL deltac_solve(delta,fulldelta)
      delta=delta*in%sfac**(2.0*in%p1/3.0)*in%v1**(2.0*in%p1)
      fulldelta=fulldelta*in%sfac**(2.0*in%p1/3.0)*in%v1**(2.0*in%p1)
      tmp=delta(1)
      delta(1)=delta(2)
      delta(2)=tmp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_run
c-----------------------------------------------------------------------
c     subprogram 2. deltac_init_asymp.
c     initialize the asymptotic solutions at large x.
c     coefficients of two power series and two small exponential  
c     solutions are solved.
c     GJT Appendix A.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_init_asymp
      
      INTEGER :: i,j,l,t,order
      INTEGER :: info
      INTEGER, DIMENSION(3) :: ipiv
      
      REAL(r8) :: e,f,k,g,h
      COMPLEX(r8),DIMENSION(2) :: p,s
      COMPLEX(r8) :: sq,tmp,q
      
      COMPLEX(r8), DIMENSION(4) :: term
      COMPLEX(r8), DIMENSION(3) :: matb
      COMPLEX(r8), DIMENSION(3,3) :: mata,mat
c-----------------------------------------------------------------------
c     initialize the parameters (f,h,k,g,h,q) and the coefficients v
c-----------------------------------------------------------------------
      order=order_pow
      IF (order < order_exp) order=order_exp
      ALLOCATE(asp%v(3,6,order+1))
      e=in%e
      f=in%f
      k=in%k
      g=in%g
      h=in%h
      q=in%q
      asp%v=0.0
c-----------------------------------------------------------------------
c     p power of large and small solutions (T_3,4).
c-----------------------------------------------------------------------
      p(1)=-0.5+in%p1
      p(2)=-0.5-in%p1
      asp%p=p
c-----------------------------------------------------------------------
c     solve large and small power-like solutions.
c-----------------------------------------------------------------------
      DO i=1,2
         t=2+i
c-----------------------------------------------------------------------
c     zeroth order term.
c-----------------------------------------------------------------------
         j=0
         l=j+1
         asp%v(:,t,l)=1
c-----------------------------------------------------------------------
c     solve high order term.
c-----------------------------------------------------------------------
         mata(1,1)=q
         mata(1,2)=-q
         mata(1,3)=0.0
         mata(2,2)=0.0
         mata(3,1)=1
         mata(3,2)=0.0
         mata(3,3)=-1
         DO j=1,order_pow
            l=j+1
            mata(2,1)=(p(i)-2*j+h)*(p(i)-2*j+1)
            mata(2,3)=-(h*(p(i)-2*j)-e-f)

            matb(1)=(p(i)-2*j+2)*(p(i)-2*j+3)*asp%v(1,t,l-1)
     $           -h*(p(i)-2*j+2)*asp%v(3,t,l-1)
            matb(2)=-q*q*(p(i)-2*j+2)*(p(i)-2*j+1)*asp%v(2,t,l-1)
            matb(3)=h*k*q*q*(p(i)-2*j+3)*asp%v(1,t,l-1)
     $           -(g-e*k)*q*q*asp%v(2,t,l-1)
     $           +(g+f*k)*q*q*asp%v(3,t,l-1)
            IF (j > 1) matb(3)=matb(3)
     $           -q*(p(i)-2*j+4)*(p(i)-2*j+3)*asp%v(3,t,l-2)
            mat=mata
            asp%v(:,t,l)=matb
            CALL zgesv(3,1,mat,3,ipiv,asp%v(:,t,l),3,info)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     s power of T_5,6, exponential power series.
c-----------------------------------------------------------------------
      sq=sqrt(q)
      term(1)=-sq*q*(1+g+k*(f+h*h))*0.25
      term(2)=(g+k*f-1)**2.0
      tmp=k*h*h
      term(3)=tmp*(tmp+(g+k*f+1)*2.0)
      term(4)=((g-k*e-1.0)*in%dr+h*h)*4.0
      tmp=sqrt(q**3.0*(term(2)+term(3))+term(4))*0.25
      s(1)=term(1)+tmp-0.5
      s(2)=term(1)-tmp-0.5
      asp%s=s
c-----------------------------------------------------------------------
c     solve coefficients of small exponential solutions T5 and T6.
c-----------------------------------------------------------------------
      DO i=1,2
         mata(1,1)=1.0/q
         mata(1,2)=q
         mata(1,3)=h/sq
         mata(2,1)=q-h/sq
         mata(2,2)=-q*sq*(2*s(i)+1)
         mata(2,3)=e+f
         mata(3,1)=q*sq*k*h+1
         mata(3,2)=q*q*(g-k*e)
         mata(3,3)=-q*q*(g+k*f)-sq*(2.0*s(i)+1.0)
c-----------------------------------------------------------------------
c     solve zeroth order term.
c-----------------------------------------------------------------------
         t=4+i
         j=0
         l=j+1

         asp%v(1,t,l)=1
         
         matb(1)=-1.0/q
         matb(2)=-(q-h/sq)
         
         asp%v(2:3,t,l)=matb(1:2)
         mat=mata
         CALL zgesv(2,1,mat(1:2,2:3),2,ipiv(1:2),
     $        asp%v(2:3,t,l),2,info)
c-----------------------------------------------------------------------
c     solve high order.
c-----------------------------------------------------------------------
         DO j=1,order_exp
            l=j+1
            
            mata(2,2)=-q*sq*(2*s(i)-4*j+1)
            mata(3,3)=-q*q*(g+k*f)-sq*(2.0*s(i)-4*j+1.0)
            
            matb(1)=((2.0*s(i)-4.0*j+3)/sq+q)*asp%v(1,t,l-1)
            matb(1)=matb(1)+h*(s(i)-2*j+2)*asp%v(3,t,l-1)
            IF (j > 1) matb(1)=-(s(i)+3-2*j)*(s(i)+2-2*j)
     $           *asp%v(1,t,l-2)
            matb(2)=-h*(s(i)+1.0-2.0*j)*asp%v(1,t,l-1)
            matb(2)=matb(2)-q*q*(s(i)-2*j+2)*(s(i)-2*j+1)
     $           *asp%v(2,t,l-1)
            matb(3)=q*q*k*h*(s(i)+1.0-2.0*j)
     $           *asp%v(1,t,l-1)
     $           -q*(s(i)-2*j+2)*(s(i)-2*j+1)
     $           *asp%v(3,t,l-1)
            mat=mata
            asp%v(:,t,l)=matb
            CALL zgesv(3,1,mat,3,ipiv,asp%v(:,t,l),3,info)
         ENDDO
      ENDDO     
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE deltac_init_asymp
c-----------------------------------------------------------------------
c     subprogram 3. deltac_get_ua.
c     get the value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_get_ua(xin,ua)

      REAL(r8),INTENT(IN) :: xin
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: ua
      
      INTEGER :: i,j,t,l
      REAL(r8) :: x,x2
      COMPLEX(r8) :: xp,xs,xj
c-----------------------------------------------------------------------
c     change sign of x.
c-----------------------------------------------------------------------
      IF (xin<0) THEN
         x=-xin
      ELSE
         x=xin
      ENDIF
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      ua=0.0
      x2=x*x
      DO i=1,2
         xp=x**asp%p(i)
         IF (rescale) THEN
            xs=EXP((-x2+xmax*xmax)/(SQRT(in%q)*2.0) )
     $           *((x/xmax)**asp%s(i))
         ELSE
            xs=EXP(-x2/(SQRT(in%q)*2.0) )*(x**asp%s(i))
         ENDIF
         DO j=0, order_pow
            l=j+1
            xj=x**(-2.0*j)
            t=2+i
            ua(:,t)=ua(:,t)+asp%v(:,t,l)*xj         
         ENDDO

         DO j=0, order_exp
            l=j+1
            xj=x**(-2.0*j)
            t=4+i
            ua(:,t)=ua(:,t)+asp%v(:,t,l)*xj
         ENDDO

         t=2+i
         ua(:,t)=xp*ua(:,t)
         ua(1,t)=x*ua(1,t)
         t=4+i
         ua(:,t)=xs*ua(:,t)
         ua(1,t)=ua(1,t)/x
      ENDDO
c-----------------------------------------------------------------------
c     extend to negtive x with even parity.
c-----------------------------------------------------------------------
      IF (xin<0 .AND. fulldomain .EQ. 1) THEN
         ua(1,:)=-ua(1,:)
c         ua(2,:)=-ua(2,:)
c         ua(3,:)=-ua(3,:)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_get_ua
c-----------------------------------------------------------------------
c     subprogram 4. deltac_get_dua.
c     get the derivative value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_get_dua(xin,dua)

      REAL(r8),INTENT(IN) :: xin
      REAL(r8) :: x
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: dua
      
      INTEGER i,j,t,l
      REAL(r8),DIMENSION(2) :: p
      REAL(r8) :: x2
      COMPLEX(r8) :: xp,xs,xj,q
      COMPLEX(r8),DIMENSION(2) :: s
c-----------------------------------------------------------------------
c     change sign of x.
c-----------------------------------------------------------------------
      IF (xin<0) THEN
         x=-xin
      ELSE
         x=xin
      ENDIF
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      q=in%q
      dua=0.0
      p=asp%p
      s=asp%s
      x2=x*x
      DO i=1,2
         xp=x**p(i)
         IF (rescale) THEN
            xs=EXP((-x2+xmax*xmax)/(SQRT(in%q)*2.0) )
     $           *((x/xmax)**asp%s(i))
         ELSE
            xs=EXP(-x2/(SQRT(in%q)*2.0) )*(x**asp%s(i))
         ENDIF
         DO j=0, order_pow
            l=j+1
            t=2+i
            xj=x**(-1.0-2.0*j)
            dua(1,t)=dua(1,t)+asp%v(1,t,l)*(p(i)+1-2*j)*xj
            dua(2:3,t)=dua(2:3,t)+asp%v(2:3,t,l)*(p(i)-2*j)*xj
         ENDDO
         DO j=0, order_exp            
            l=j+1
            t=4+i
            xj=x**(-2.0*j)
            dua(1,t)=dua(1,t)
     $           +asp%v(1,t,l)*((s(i)-1-2*j)/x2-1.0/SQRT(q))*xj
            dua(2:3,t)=dua(2:3,t)
     $           +asp%v(2:3,t,l)*((s(i)-2*j)/x-1.0/SQRT(q)*x)*xj
         ENDDO
         t=2+i
         dua(:,t)=xp*dua(:,t)
         dua(1,t)=x*dua(1,t)
         t=4+i
         dua(:,t)=xs*dua(:,t)
      ENDDO
c-----------------------------------------------------------------------
c     extend to negtive x with even parity.
c-----------------------------------------------------------------------
      IF (xin<0 .AND. fulldomain .EQ. 1) THEN
c         dua(1,:)=-dua(1,:)
         dua(2,:)=-dua(2,:)
         dua(3,:)=-dua(3,:)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_get_dua
c-----------------------------------------------------------------------
c     subprogram 5. deltac_solve.
c     solves the inner layer model with galerkin method.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_solve(delta,fulldelta)

      COMPLEX(r8), DIMENSION(2), INTENT(OUT) :: delta
      COMPLEX(r8), DIMENSION(2,2), INTENT(OUT) :: fulldelta
      INTEGER :: info,isol,imap,ix
      TYPE(gal_type) :: gal
10    FORMAT("D11=",1p,2e11.3,"D12=",1p,2e11.3,
     $       "D21=",1p,2e11.3,"D22=",1p,2e11.3/)
c-----------------------------------------------------------------------
c     allocate and compute arrays.
c-----------------------------------------------------------------------
      IF (basis_type == 0) np=3
      IF (np <= 3)np=3
      msol=2
      CALL deltac_alloc(gal,nx,nq,dx1,dx2,pfac)
      CALL deltac_make_arrays(gal)
c-----------------------------------------------------------------------
c     solve the galerkin matrix.
c-----------672------------------------------------------------------------
      DO isol=1,msol
         gal%sol(:,isol)=gal%rhs(:,isol)
         CALL zgbtrf(gal%ndim,gal%ndim,gal%kl,gal%ku,gal%mat(:,:,isol),
     $        gal%ldab,gal%ipiv,info)
         CALL zgbtrs("N",gal%ndim,gal%kl,gal%ku,1,gal%mat(:,:,isol),
     $        gal%ldab,gal%ipiv,gal%sol(:,isol),gal%ndim,info)
      ENDDO
c-----------------------------------------------------------------------
c     compute and write delta.
c-----------------------------------------------------------------------
      SELECT CASE (gal_method)
      CASE ("normal")
         SELECT CASE (fulldomain)
         CASE (0)
            DO isol=1,gal%msol
               imap=gal%intvl%cell(nx)%map(1,4)
               delta(isol)=gal%sol(imap,isol)
            ENDDO
c            WRITE(*,*)"delta1=",delta(1),"delta2=",delta(2)
c            CALL program_stop("delta+- stop")
         CASE (1,2)
            DO isol=1,gal%msol
               imap=gal%intvl%cell(-nx)%map(1,-1)
               fulldelta(isol,1)=gal%sol(imap,isol)
               imap=gal%intvl%cell(nx)%map(1,4)
               fulldelta(isol,2)=gal%sol(imap,isol)
            ENDDO
c            WRITE (*,10) fulldelta(1,1),fulldelta(1,2),
c     $                   fulldelta(2,1),fulldelta(2,2)
c            CALL program_stop
c     $        ("Finish full domain comptation with normal method.")
         END SELECT
      CASE ("resonant")
         DO isol=1,gal%msol
            DO ix=1,gal%nx
               cell => gal%intvl%cell(ix)
               IF(cell%etype == "res")THEN
                  imap=cell%emap(1)
                  delta(isol)=gal%sol(imap,isol)
               ENDIF
            ENDDO
         ENDDO
      END SELECT
      
      IF(output_sol)CALL deltac_output_solution(gal)
c-----------------------------------------------------------------------
c     deallocate arrays and finish.
c-----------------------------------------------------------------------
      CALL deltac_dealloc(gal)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_solve
c-----------------------------------------------------------------------
c     subprogram 6. deltac_alloc.
c     allocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_alloc(gal,nx,nq,dx1,dx2,pfac)
      
      TYPE(gal_type), INTENT(OUT) :: gal
      INTEGER, INTENT(IN) :: nx,nq
      REAL(r8), INTENT(IN) :: dx1,dx2,pfac

      INTEGER :: ix,ixmin
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     set integers.
c-----------------------------------------------------------------------
	SELECT CASE (fulldomain)
      CASE(0)
         gal%nx=nx
      CASE(1,2)
         gal%nx=nx
      END SELECT
      gal%nq=nq
      gal%dx1=dx1
      gal%dx2=dx2
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(gal%intvl)
      CALL jacobi_alloc(gal%quad,nq,.TRUE.,.FALSE.,"gll")
      intvl => gal%intvl
      SELECT CASE (fulldomain)
      CASE(0)
         ALLOCATE(intvl%x(0:nx),intvl%dx(0:nx),intvl%cell(nx))
         ixmin=1
      CASE(1,2)
         ALLOCATE(intvl%x(-nx:nx),intvl%dx(-nx:nx),intvl%cell(-nx:nx))
         ixmin=-nx
      CASE DEFAULT
         WRITE(*,'(a,i2)') 
     $        "deltac_run: invalide value fulldomain = ",fulldomain
         STOP
      END SELECT
         
      CALL deltac_make_grid(nx,dx1,dx2,pfac,intvl)
c-----------------------------------------------------------------------
c     allocate map and local matrix for each element.
c-----------------------------------------------------------------------
      DO ix=ixmin,nx
         IF (ix == 0) CYCLE
         cell => intvl%cell(ix)
         IF (basis_type == 0) THEN
            SELECT CASE(cell%etype)
c-----------------------------------------------------------------------
c     allocate normal element.
c-----------------------------------------------------------------------
            CASE ("none")
               ALLOCATE(cell%map(mpert,0:np))
               ALLOCATE(cell%mat(mpert,mpert,0:np,0:np))
               cell%np=np
               cell%x_lsode=0.0
c-----------------------------------------------------------------------
c     allocate resonant element.
c-----------------------------------------------------------------------
            CASE ("res")
               IF (method) THEN
c-----------------------------------------------------------------------
c     no hermite basis in resonant element.
c-----------------------------------------------------------------------
                  ALLOCATE(cell%map(3,0:0),cell%emap(3))
                  ALLOCATE(cell%mat(3,3,0:0,0:0))
                  ALLOCATE(cell%rhs(3,0:0))
                  cell%np=-1
               ELSE
c-----------------------------------------------------------------------
c     include hermite basis in resonant element.
c-----------------------------------------------------------------------
                  ALLOCATE(cell%map(mpert,0:2),cell%emap(3))
                  ALLOCATE(cell%mat(mpert,mpert,0:2,0:2))
                  ALLOCATE(cell%rhs(mpert,0:2))
                  cell%np=1
               ENDIF
               cell%rhs=0.0
               cell%emap=0.0
               cell%x_lsode=cell%x(2)
               IF (ix<0) cell%x_lsode=cell%x(1)
c-----------------------------------------------------------------------
c     allocate extension element connect to resonant element.
c-----------------------------------------------------------------------
            CASE ("ext")
               IF (basis_type == 0) THEN
                  IF (method) THEN
c-----------------------------------------------------------------------
c     no hermite basis in resonant element.
c-----------------------------------------------------------------------
                  ALLOCATE(cell%map(mpert,0:2),cell%emap(3))
                  ALLOCATE(cell%mat(mpert,mpert,0:2,0:2))
                  ALLOCATE(cell%rhs(mpert,0:2))
                  cell%np=1
                  ELSE
c-----------------------------------------------------------------------
c     include hermite basis in resonant element.
c-----------------------------------------------------------------------
                  ALLOCATE(cell%map(mpert,0:np+1),cell%emap(3))
                  ALLOCATE(cell%mat(mpert,mpert,0:np+1,0:np+1))
                  ALLOCATE(cell%rhs(mpert,0:np+1))
                  cell%np=np
                  ENDIF
               ELSE 
                  ALLOCATE(cell%map(mpert,0:np),cell%emap(3))
                  ALLOCATE(cell%mat(mpert,mpert,0:np,0:np))
                  ALLOCATE(cell%rhs(mpert,0:np))

                  cell%np=np-1               
               ENDIF
               cell%emap=0.0
               cell%rhs=0.0               
               cell%x_lsode=0.0
c-----------------------------------------------------------------------
c     allocate extension element only has driving term.
c-----------------------------------------------------------------------
            CASE ("ext1","ext2")
               SELECT CASE (gal_method)
               CASE("resonant")
                  ALLOCATE(cell%map(mpert,0:np))
                  ALLOCATE(cell%mat(mpert,mpert,0:np,0:np))
                  ALLOCATE(cell%rhs(mpert,0:np))
               CASE("normal")
                  IF (ix == nx) THEN
                     ALLOCATE(cell%map(mpert,0:np+1))
                     ALLOCATE(cell%mat(mpert,mpert,0:np+1,0:np+1))
                     ALLOCATE(cell%rhs(mpert,0:np+1))
                  ELSEIF (ix==-nx) THEN
                     ALLOCATE(cell%map(mpert,-1:np))
                     ALLOCATE(cell%mat(mpert,mpert,-1:np,-1:np))
                     ALLOCATE(cell%rhs(mpert,-1:np))

                  ELSE
                     ALLOCATE(cell%map(mpert,0:np))
                     ALLOCATE(cell%mat(mpert,mpert,0:np,0:np))
                     ALLOCATE(cell%rhs(mpert,0:np))
                  ENDIF               
               END SELECT
               cell%np=np
               cell%rhs=0.0
               cell%x_lsode=0.0
            END SELECT
         ELSE
         
         ENDIF
         cell%map=0.0
         cell%mat=0.0         
      ENDDO

c      IF(diagnose_grid)CALL deltac_diagnose_grid(gal)
c-----------------------------------------------------------------------
c     create and diagnose local-to-global mapping.
c-----------------------------------------------------------------------
      IF (basis_type==0) THEN
         CALL deltac_make_map_hermite(gal)
      ELSE
      
      ENDIF
c-----------------------------------------------------------------------
c     allocate global arrays for LU solver.
c     note: kl considers both hermite and spectrum element's cases
c-----------------------------------------------------------------------
      SELECT CASE (gal_method)
      CASE ("normal")
         gal%kl=mpert*(np+2)-1
      CASE ("resonant")
         IF (noexp) THEN
            gal%kl=mpert*(np+1)+1-1
         ELSE
            gal%kl=mpert*(np+2)-1
         ENDIF
      END SELECT
      gal%ku=gal%kl
      gal%ldab=2*gal%kl+gal%ku+1    
      gal%msol=msol
      ALLOCATE(gal%rhs(gal%ndim,gal%msol),gal%sol(gal%ndim,gal%msol))
      ALLOCATE(gal%mat(gal%ldab,gal%ndim,2),gal%ipiv(gal%ndim))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_alloc      
c-----------------------------------------------------------------------
c     subprogram 7. deltac_dealloc.
c     deallocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_dealloc(gal)
      
      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ix,ixmin
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(asp%v)
      intvl => gal%intvl
      ixmin=1
      IF (fulldomain>0) ixmin=-gal%nx
      DO ix=ixmin,gal%nx
         IF (ix==0) CYCLE
         cell => intvl%cell(ix)
         DEALLOCATE(cell%map,cell%mat)
         IF(cell%etype /= "none") DEALLOCATE(cell%rhs)
         IF(cell%etype == "res"  .OR.  cell%etype == "ext")
     $      DEALLOCATE(cell%emap)            
      ENDDO
      DEALLOCATE(intvl%x,intvl%dx,intvl%cell)
      DEALLOCATE(gal%rhs,gal%sol,gal%mat,gal%ipiv,gal%intvl)
      CALL jacobi_dealloc(gal%quad)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_dealloc      
c-----------------------------------------------------------------------
c     subprogram 8. deltac_make_grid.
c     sets up grid in the interval.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_make_grid(nx,dx1,dx2,pfac,intvl)

      INTEGER :: nx
      REAL(r8), INTENT(IN) :: dx1,dx2,pfac
      TYPE(interval_type), INTENT(INOUT) :: intvl

      INTEGER :: ixmin,ixmax,ix,mx
      REAL(r8) :: x0,x1,xm,dx
c-----------------------------------------------------------------------
c     set default extra.
c----------------------------------------------------------------------- 
      DO ix=1,nx
         intvl%cell(ix)%etype="none"
      ENDDO
c-----------------------------------------------------------------------
c     set lower bound.
c-----------------------------------------------------------------------
      ixmin=0
      intvl%x(ixmin)=0
c-----------------------------------------------------------------------
c     set upper bound.
c-----------------------------------------------------------------------
      x1=xmax
      intvl%x(nx)=x1
      IF (dx1dx2_flag) THEN
         intvl%x(nx-1)=x1-dx1
         intvl%x(nx-2)=x1-(dx1+dx2)
         ixmax=nx-2
      ELSE
         ixmax=nx
      ENDIF
      
      SELECT CASE (gal_method)
      CASE("normal")
         intvl%cell(nx)%etype="ext1"
         DO ix=nx-1,nx-cutoff,-1
            intvl%cell(ix)%etype="ext1"
         ENDDO
      CASE("resonant")
         intvl%cell(nx-1)%etype="ext"
         intvl%cell(nx)%etype="res"      
         DO ix=nx-2,nx-cutoff,-1
            intvl%cell(ix)%etype="ext1"
         ENDDO
      CASE DEFAULT
         CALL program_stop
     $        ("invalid galerkin method.")
      END SELECT
      intvl%cell(ix)%etype="ext2"
c-----------------------------------------------------------------------
c     compute interior packed grid.
c-----------------------------------------------------------------------
      x0=intvl%x(ixmin)
      x1=intvl%x(ixmax)
      xm=(intvl%x(ixmax)+intvl%x(ixmin))/2
      dx=(intvl%x(ixmax)-intvl%x(ixmin))/2
      mx=(ixmax-ixmin)/2
c     check pack output      
      IF(pfac < 1) side="left"
      intvl%x(ixmin:ixmax)=xm+dx*deltac_pack(mx,pfac,side)
      intvl%x(ixmin)=x0
      intvl%x(ixmax)=x1
      intvl%dx(0)=0
      intvl%dx(1:nx)=intvl%x(1:nx)-intvl%x(0:nx-1)
      DO ix=1,nx
         intvl%cell(ix)%x=(/intvl%x(ix-1),intvl%x(ix)/)
      ENDDO
      IF (fulldomain>0 .AND. gal_method=="normal") THEN
         DO ix=-nx,-1
            intvl%cell(ix)%etype=intvl%cell(-ix)%etype
            intvl%x(ix)=intvl%x(-ix)
            intvl%dx(ix)=intvl%dx(-ix)
            intvl%cell(ix)%x(1)=-intvl%cell(-ix)%x(2)
            intvl%cell(ix)%x(2)=-intvl%cell(-ix)%x(1)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_make_grid      
c-----------------------------------------------------------------------
c     subprogram 9. deltac_pack.
c     computes packed grid on (0,1).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION deltac_pack(nx,pfac,side) RESULT(x)

      INTEGER, INTENT(IN) :: nx
      REAL(r8), INTENT(IN) :: pfac
      CHARACTER(*), INTENT(IN) :: side
      REAL(r8), DIMENSION(-nx:nx) :: x

      INTEGER :: ix
      REAL(r8) :: lambda,a
      REAL(r8), DIMENSION(-nx:nx) :: xi,num,den
c-----------------------------------------------------------------------
c     fill logical grid.
c-----------------------------------------------------------------------
      SELECT CASE(side)
      CASE("left")
         xi=(/(ix,ix=-2*nx,0)/)/REAL(2*nx,r8)
      CASE("right")
         xi=(/(ix,ix=0,2*nx)/)/REAL(2*nx,r8)
      CASE("both")
         xi=(/(ix,ix=-nx,nx)/)/REAL(nx,r8)
      CASE DEFAULT
         CALL program_stop
     $        ("deltac_pack: cannot recognize side = "//TRIM(side))
      END SELECT
c-----------------------------------------------------------------------
c     compute grid.
c-----------------------------------------------------------------------
      IF(pfac > 1)THEN
         lambda=SQRT(1-1/pfac)
         num=LOG((1+lambda*xi)/(1-lambda*xi))
         den=LOG((1+lambda)/(1-lambda))
         x=num/den
      ELSEIF(pfac == 1)THEN
         x=xi
      ELSEIF(pfac > 0)THEN
         lambda=SQRT(1-pfac)
         a=LOG((1+lambda)/(1-lambda))
         num=EXP(a*xi)-1
         den=EXP(a*xi)+1
         x=num/den/lambda
      ELSE
         x=xi**(-pfac)
      ENDIF
c-----------------------------------------------------------------------
c     shift output to (-1,1)
c-----------------------------------------------------------------------
      SELECT CASE(side)
      CASE("left")
         x=2*x+1
      CASE("right")
         x=2*x-1
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION deltac_pack      
c-----------------------------------------------------------------------
c     subprogram 10. deltac_hermite.
c     computes hermite cubic basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_hermite(x,x0,x1,hermite)

      REAL(r8),INTENT(IN) :: x,x0,x1
      TYPE(hermite_type),INTENT(INOUT) :: hermite 
      
      REAL(r8) :: dx,t0,t1,t02,t12
c-----------------------------------------------------------------------
c     compute variables.
c-----------------------------------------------------------------------
      dx=x1-x0
      t0=(x-x0)/dx
      t1=1-t0
      t02=t0*t0
      t12=t1*t1
c-----------------------------------------------------------------------
c     compute hermite basis.
c-----------------------------------------------------------------------
      hermite%pb(0)=t12*(1+2*t0)
      hermite%pb(1)=t12*t0*dx
      hermite%pb(2)=t02*(1+2*t1)
      hermite%pb(3)=-t02*t1*dx
c-----------------------------------------------------------------------
c     compute derivative of hermite basis.
c-----------------------------------------------------------------------
      hermite%qb(0)=-6*t0*t1/dx
      hermite%qb(1)=t0*(3*t0-4)+1
      hermite%qb(2)=6*t1*t0/dx
      hermite%qb(3)=t0*(3*t0-2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_hermite      
c-----------------------------------------------------------------------
c     subprogram 11. deltac_make_map_hermite.
c     creates local-to-global mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_make_map_hermite(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ix,ip,imap,ipert,nx,ixmin
      INTEGER,DIMENSION(3) :: emap
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over grid cells.
c-----------------------------------------------------------------------
      intvl => gal%intvl
      SELECT CASE(gal_method)
      CASE ("normal")
         nx=gal%nx
      CASE ("resonant")
         nx=gal%nx-1
      END SELECT
      imap=1
      ixmin=1
      IF (fulldomain>0) THEN
         IF (gal_method=="normal") THEN
            ixmin=-nx
            cell => intvl%cell(ixmin)
            emap=(/imap,imap,imap/)+(/0,1,2/)
            cell%map(:,-1)=emap
            imap=imap+3
            DO ip=0,cell%np
               DO ipert=1,mpert
                  cell%map(ipert,ip)=imap
                  imap=imap+1
               ENDDO   
            ENDDO
            ixmin=-nx+1
         ENDIF
      ENDIF
      DO ix=ixmin,nx
         IF (ix==0) CYCLE
         cell => intvl%cell(ix)
c-----------------------------------------------------------------------
c     nonresonant basis functions.
c-----------------------------------------------------------------------
         IF(imap > 1)imap=imap-2*mpert
         DO ip=0,cell%np
            DO ipert=1,mpert
               cell%map(ipert,ip)=imap
                  imap=imap+1
            ENDDO
         ENDDO
      ENDDO
      SELECT CASE(gal_method)
      CASE("normal")
         cell => intvl%cell(gal%nx)
         emap=(/imap,imap,imap/)+(/0,1,2/)
         cell%map(:,cell%np+1)=emap
         gal%ndim=imap+2
      CASE("resonant")
         IF (method) THEN
c-----------------------------------------------------------------------
c     resonant basis functions.
c-----------------------------------------------------------------------
         emap=(/imap,imap,imap/)+(/0,1,2/)
         intvl%cell(gal%nx-1)%map(:,2)=emap
         intvl%cell(gal%nx-1)%emap=emap
         intvl%cell(gal%nx)%emap=emap
         intvl%cell(gal%nx)%map(:,0)=emap
         IF (noexp) THEN
            gal%ndim=imap
         ELSE
            gal%ndim=imap+2
         ENDIF
         ELSE
c-----------------------------------------------------------------------
c     resonant basis functions.
c-----------------------------------------------------------------------
         emap=(/imap,imap,imap/)+(/0,1,2/)
         intvl%cell(gal%nx-1)%map(:,4)=emap
         intvl%cell(gal%nx-1)%emap=emap
         intvl%cell(gal%nx)%emap=emap
         intvl%cell(gal%nx)%map(:,:)=intvl%cell(gal%nx-1)%map(:,2:4)
         gal%ndim=imap+2
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_make_map_hermite
            
c-----------------------------------------------------------------------
c     subprogram 12. deltac_make_arrays.
c     computes matrix and rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_make_arrays(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ix,ip,jp,ep,ixmin
      REAL(r8) :: x1
      COMPLEX(r8), DIMENSION(3,4) :: u_res
      COMPLEX(r8), DIMENSION(mpert,0:1,1:4) :: u_h1
      COMPLEX(r8), DIMENSION(mpert,1:3,0:1) :: u_h2  
      COMPLEX(r8), DIMENSION(1,1) :: term
      COMPLEX(r8), DIMENSION(1,3) :: ua1
      COMPLEX(r8), DIMENSION(3,3) :: imat,umat,vmat
      COMPLEX(r8), DIMENSION(3,1) :: dua1
      COMPLEX(r8), DIMENSION(3,6) :: ua,dua      
      TYPE(interval_type), POINTER :: intvl
c-----------------------------------------------------------------------
c     start loops over grid cells.
c-----------------------------------------------------------------------
      intvl => gal%intvl
      ixmin=1
      IF (fulldomain>0.AND.gal_method=="normal") THEN
         ixmin=-gal%nx
      ENDIF
      DO ix=ixmin,gal%nx
         IF (ix == 0) CYCLE
         cell => intvl%cell(ix)
c-----------------------------------------------------------------------
c     nonresonant arrays.
c-----------------------------------------------------------------------
         CALL deltac_gauss_quad(gal%quad)
c-----------------------------------------------------------------------
c     resonant arrrays.
c-----------------------------------------------------------------------
         SELECT CASE(cell%etype)
         CASE("res")
            CALL deltac_lsode_int(u_res,u_h1,u_h2)
            ep=cell%np+1
            IF (method) THEN
               cell%mat(:,:,ep,ep)=cell%mat(:,:,ep,ep)+u_res(:,1:3)
               cell%rhs(:,ep)=cell%rhs(:,ep)-u_res(:,4)
            ELSE
               DO ip=0,1
                  DO jp=1,3
                     cell%mat(:,jp,ip,ep)=cell%mat(:,jp,ip,ep)
     $                    +u_h1(:,ip,jp)
                  ENDDO
               ENDDO
               DO ip=1,3
                  DO jp=0,1
                     cell%mat(ip,:,ep,jp)=cell%mat(ip,:,ep,jp)
     $                    +u_h2(:,ip,jp)
                  ENDDO
               ENDDO
               
               cell%mat(:,:,ep,ep)=cell%mat(:,:,ep,ep)+u_res(:,1:3)
               cell%rhs(:,0:1)=cell%rhs(:,0:1)-u_h1(:,:,4)
               cell%rhs(:,ep)=cell%rhs(:,ep)-u_res(:,4)
            ENDIF
         CASE("ext","ext1","ext2")
            CALL deltac_extension(gal%quad)
         END SELECT
c-----------------------------------------------------------------------
c     finish loops over grid cells.
c-----------------------------------------------------------------------
      ENDDO
c-----------------------------------------------------------------------
c     implement boundary condition at xmax.
c-----------------------------------------------------------------------
      cell => intvl%cell(gal%nx)
      SELECT CASE(gal_method)
      CASE ("normal")
         x1=cell%x(2)
         ep=cell%np+1
         CALL deltac_get_ua(x1,ua)
         CALL deltac_get_dua(x1,dua)
         cell%mat(:,:,2,:)=0
         cell%mat(1,1,2,2)=1
         cell%mat(2,2,2,2)=1
         cell%mat(3,3,2,2)=1
         cell%mat(:,1,2,ep)=-ua(:,tid(1))         
         cell%mat(:,2,2,ep)=-ua(:,tid(2))         
         cell%mat(:,3,2,ep)=-ua(:,tid(3))         
         
         cell%mat(1,1,ep,3)=1
         cell%mat(2,2,ep,3)=1
         cell%mat(3,3,ep,3)=1
         cell%mat(:,1,ep,ep)=-dua(:,tid(1))         
         cell%mat(:,2,ep,ep)=-dua(:,tid(2))         
         cell%mat(:,3,ep,ep)=-dua(:,tid(3))         

         cell%rhs(:,2)=0
         cell%rhs(:,ep)=0

         IF (fulldomain>0) THEN
            cell => intvl%cell(-gal%nx)
            x1=cell%x(1)
C           CHECK EP
            ep=-1
            CALL deltac_get_ua(x1,ua)
            CALL deltac_get_dua(x1,dua)
            cell%mat(:,:,0,:)=0
            cell%mat(1,1,0,0)=1
            cell%mat(2,2,0,0)=1
            cell%mat(3,3,0,0)=1
            cell%mat(:,1,0,ep)=-ua(:,tid(1))         
            cell%mat(:,2,0,ep)=-ua(:,tid(2))         
            cell%mat(:,3,0,ep)=-ua(:,tid(3))         
         
            cell%mat(1,1,ep,1)=1
            cell%mat(2,2,ep,1)=1
            cell%mat(3,3,ep,1)=1
            cell%mat(:,1,ep,ep)=-dua(:,tid(1))         
            cell%mat(:,2,ep,ep)=-dua(:,tid(2))         
            cell%mat(:,3,ep,ep)=-dua(:,tid(3))         

            cell%rhs(:,0)=0
            cell%rhs(:,ep)=0

         ENDIF
c-----------------------------------------------------------------------
c     restore the surface term at xmax.
c-----------------------------------------------------------------------
      CASE ("resonant")
         x1=cell%x_lsode
         ep=cell%np+1
         CALL deltac_get_ua(x1,ua)
         CALL deltac_get_dua(x1,dua)
         CALL deltac_get_uv(x1,imat,umat,vmat)
         DO ip=1,3
            ua1(1,:)=ua(:,tid(ip))
            DO jp=1,3
               dua1(:,1)=dua(:,tid(jp))
               term=MATMUL(ua1,MATMUL(imat,dua1))
               cell%mat(ip,jp,ep,ep)=cell%mat(ip,jp,ep,ep)-term(1,1)
            ENDDO
            dua1(:,1)=dua(:,tid(4))
            term=MATMUL(ua1,MATMUL(imat,dua1))
            cell%rhs(ip,ep)=cell%rhs(ip,ep)+term(1,1)
         ENDDO      
      END SELECT
c-----------------------------------------------------------------------
c     assemble matrix and rhs.
c-----------------------------------------------------------------------
      CALL deltac_assemble_mat(gal)
      CALL deltac_assemble_rhs(gal)
c-----------------------------------------------------------------------
c     set boundary condition for odd and even parity.
c-----------------------------------------------------------------------
      IF (fulldomain.EQ.0) THEN
         CALL deltac_set_boundary(gal)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_make_arrays
c-----------------------------------------------------------------------
c     subprogram 13. deltac_gauss_quad.
c     Gauss quadrature evaluation of nonresonant matrix elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_gauss_quad(quad)

      TYPE(jacobi_type) :: quad

      INTEGER :: ip,jp,iq
      REAL(r8) :: x0,dx,x,w
      REAL(r8), DIMENSION(0:np) :: pb,qb
      COMPLEX(r8), DIMENSION(mpert,mpert) :: imat,umat,vmat
      TYPE(hermite_type) :: hermite
c-----------------------------------------------------------------------
c     start loop over quadrature points.
c-----------------------------------------------------------------------
      cell%mat=0.0
      DO iq=0,quad%np
         x0=(cell%x(2)+cell%x(1))/2
         dx=(cell%x(2)-cell%x(1))/2
         x=x0+dx*quad%node(iq)
         w=dx*quad%weight(iq)
         CALL deltac_get_uv(x,imat,umat,vmat)
         IF (basis_type == 0) THEN
            CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)
            pb=hermite%pb
            qb=hermite%qb
         ELSE
         
         ENDIF
c-----------------------------------------------------------------------
c     compute gaussian quadratures.
c-----------------------------------------------------------------------
         DO ip=0,cell%np
            DO jp=0,cell%np
               cell%mat(:,:,ip,jp)
     $              =cell%mat(:,:,ip,jp)
     $              +w*(imat*qb(ip)*qb(jp)
     $              +vmat*pb(ip)*qb(jp)
     $              +umat*pb(ip)*pb(jp))     
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     finish loop over quadrature points.
c-----------------------------------------------------------------------
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_gauss_quad   
c-----------------------------------------------------------------------
c     subprogram 14. deltac_get_uv.
c     computes U, V matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_get_uv(x,imat,umat,vmat)

      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: imat,umat,vmat

      COMPLEX(r8) :: zero=0
c-----------------------------------------------------------------------
c     compute each coffeicent matrix.
c-----------------------------------------------------------------------
      imat=0.0
      imat(1,1)=1.0
      imat(2,2)=in%q*in%q
      imat(3,3)=in%q
      
      umat=RESHAPE((/ in%q,   -x/in%q,             -x/in%q,
     $   -in%q*x, x*x/in%q,   -(in%g-in%k*in%e)*in%q,
     $   zero,-(in%e+in%f)/(in%q*in%q), x*x/in%q+(in%g+in%k*in%f)*in%q/)
     $  ,SHAPE(umat))

          
      vmat=RESHAPE((/  zero,   -in%h/(in%q*in%q),  in%h*in%k*in%q,   
     $     zero,            zero,   zero,
     $     dcmplx(in%h),    zero,   zero/)
     $  ,SHAPE(vmat))
     
      umat(2,:)=umat(2,:)*in%q*in%q
      umat(3,:)=umat(3,:)*in%q    
      vmat(2,:)=vmat(2,:)*in%q*in%q
      vmat(3,:)=vmat(3,:)*in%q    
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_get_uv
c-----------------------------------------------------------------------
c     subprogram 15. deltac_assemble_mat.
c     assembles global matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_assemble_mat(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ix,ip,ipert,i,jp,jpert,j,offset,npp
      INTEGER :: ix0,np0
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals, grid cells, and rows.
c-----------------------------------------------------------------------
      offset=gal%kl+gal%ku+1
      gal%mat=0
      intvl => gal%intvl
      ix0=1
      IF (fulldomain>0) ix0=-gal%nx
      DO ix=ix0,gal%nx
         IF (ix==0) CYCLE
         cell => intvl%cell(ix)
         npp=cell%np
         np0=0
         SELECT CASE (gal_method)
         CASE("resonant")
            IF (cell%etype == "ext"  .OR.  cell%etype == "res") THEN
               npp=cell%np+1 
            ENDIF
         CASE("normal")
            IF (cell%etype == "ext1") THEN
               IF (ix == gal%nx)   npp=cell%np+1
               IF (ix == -gal%nx)  np0=-1
            ENDIF
         END SELECT
         DO ip=np0,npp
            DO ipert=1,mpert
               i=cell%map(ipert,ip)
               DO jp=np0,npp
                  DO jpert=1,mpert
                     j=cell%map(jpert,jp)
                     IF (i>gal%ndim.OR.j>gal%ndim)THEN
                        CYCLE
                     ENDIF
                     gal%mat(offset+i-j,j,1)=gal%mat(offset+i-j,j,1)
     $                    +cell%mat(ipert,jpert,ip,jp)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      gal%mat(:,:,2)=gal%mat(:,:,1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_assemble_mat
c-----------------------------------------------------------------------
c     subprogram 16. deltac_assemble_rhs.
c     assembles global rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_assemble_rhs(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ix,imap,ipert,ip,npp
      INTEGER :: ix0,np0
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      gal%rhs=0
      intvl => gal%intvl
      ix0=1
      IF (fulldomain>0) ix0=-gal%nx
      DO ix=ix0,gal%nx
         IF (ix==0) CYCLE
         cell => intvl%cell(ix)
         IF (cell%etype == "none") CYCLE
         npp=cell%np
         np0=0
         SELECT CASE(gal_method)
         CASE("resonant")
            IF (cell%etype == "ext"  .OR.  cell%etype == "res") THEN
               npp=cell%np+1
            ENDIF
         CASE("normal")
            IF (cell%etype == "ext1") THEN
               IF (ix == gal%nx)  npp=cell%np+1
               IF (ix == -gal%nx) np0=-1
            ENDIF
         END SELECT
         DO ip=np0,npp
            DO ipert=1,mpert
               imap=cell%map(ipert,ip)
               IF (imap>gal%ndim) THEN
                  WRITE (*,*) "ERROR: imap error."
                  CYCLE
               ENDIF
               SELECT CASE (fulldomain)
               CASE (0)
                  gal%rhs(imap,1)=gal%rhs(imap,1)+cell%rhs(ipert,ip)
                  gal%rhs(imap,2)=gal%rhs(imap,1)
               CASE (1,2)
                  IF (ix<0) THEN
                     gal%rhs(imap,2)=gal%rhs(imap,2)+cell%rhs(ipert,ip)
                  ELSEIF (ix>0) THEN
                     gal%rhs(imap,1)=gal%rhs(imap,1)+cell%rhs(ipert,ip)
                  ENDIF
               END SELECT
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_assemble_rhs
c-----------------------------------------------------------------------
c     subprogram 17. deltac_set_boundary.
c     sets boundary conditions at x=0.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_set_boundary(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      COMPLEX(r8), DIMENSION(mpert,mpert,0:np,0:np) :: mat
      INTEGER :: isol,i,j,ip,jp,ipert,jpert,offset
      COMPLEX(r8), DIMENSION(mpert,mpert) :: umat,vmat,imat
c-----------------------------------------------------------------------
c     restore surface term at x=0.
c-----------------------------------------------------------------------
      cell => gal%intvl%cell(1)
      IF (basis_type == 0) THEN
         CALL deltac_get_uv(cell%x(1),imat,umat,vmat)
         cell%mat(:,:,0,1)=cell%mat(:,:,0,1)+imat
      ELSE
      ENDIF
c-----------------------------------------------------------------------
c     set boundary condition at u(0).
c-----------------------------------------------------------------------
      DO isol=1,2
         mat=cell%mat
         IF (basis_type == 0) THEN
            IF (isol == 1) THEN
c-----------------------------------------------------------------------
c     set boundary condition for odd modes.
c-----------------------------------------------------------------------
               mat(:,:,0,:)=0.0
               mat(1,1,0,1)=1.0
               mat(2,2,0,0)=1.0
               mat(3,3,0,0)=1.0               
            ELSE
c-----------------------------------------------------------------------
c     set boundary condition for even modes.
c-----------------------------------------------------------------------
               mat(:,:,0,:)=0.0
               mat(1,1,0,0)=1.0
               mat(2,2,0,1)=1.0
               mat(3,3,0,1)=1.0
            ENDIF
c-----------------------------------------------------------------------
c     update boundary condition in global matrix.
c-----------------------------------------------------------------------
            offset=gal%kl+gal%ku+1
            DO ip=0,1
               DO ipert=1,mpert
                  i=cell%map(ipert,ip)
                  DO jp=0,cell%np
                     DO jpert=1,mpert
                        j=cell%map(jpert,jp)
                        gal%mat(offset+i-j,j,isol)=
     $                       mat(ipert,jpert,ip,jp)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO      
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN        
      END SUBROUTINE deltac_set_boundary    
c-----------------------------------------------------------------------
c     subprogram 18. deltac_extension.
c     computes extension terms of matrix and rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_extension(quad)

      TYPE(jacobi_type) :: quad

      INTEGER :: ep
      INTEGER :: jp,ip,iq
      REAL(r8) :: x,x0,dx,w
      REAL(r8), DIMENSION(0:np) :: pb,qb
      COMPLEX(r8), DIMENSION(1,1) :: term
      COMPLEX(r8), DIMENSION(mpert,1) :: term1,ua1,dua1
      COMPLEX(r8), DIMENSION(1,mpert) :: term2,ua2,dua2
      COMPLEX(r8), DIMENSION(3,6) :: ua,dua,uax,duax
      COMPLEX(r8), DIMENSION(mpert,mpert) :: umat,vmat,imat
      
      TYPE(hermite_type) :: hermite
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      IF (cell%etype == "ext") THEN
c-----------------------------------------------------------------------
c     get asymptotic solutions and derivatives at cell boundary.
c-----------------------------------------------------------------------
         x=cell%x(2)
         CALL deltac_get_ua(x,ua)
         CALL deltac_get_dua(x,dua)      
c-----------------------------------------------------------------------
c     set extension element.
c-----------------------------------------------------------------------
         ep=cell%np+1
         DO iq=0,quad%np
            x0=(cell%x(2)+cell%x(1))/2
            dx=(cell%x(2)-cell%x(1))/2
            x=x0+dx*quad%node(iq)
            w=dx*quad%weight(iq)
            CALL deltac_get_ua(x,uax)
            CALL deltac_get_dua(x,duax)               
            CALL deltac_get_uv(x,imat,umat,vmat)
            CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)
            pb=hermite%pb
            qb=hermite%qb
c-----------------------------------------------------------------------
c     fill matrix for (hi,Luaj) and (ua*^Ti,Lhj).
c-----------------------------------------------------------------------
            DO ip=0,cell%np
               DO jp=1,3
                  ua1(:,1)=ua(:,tid(jp))*pb(2)+dua(:,tid(jp))*pb(3)
                  dua1(:,1)=ua(:,tid(jp))*qb(2)+dua(:,tid(jp))*qb(3)
                  term1=qb(ip)*MATMUL(imat,dua1)
     $                 +pb(ip)*MATMUL(vmat,dua1)
     $                 +pb(ip)*MATMUL(umat,ua1)
                  cell%mat(:,jp,ip,ep)=cell%mat(:,jp,ip,ep)+w*term1(:,1)

c                  ua2(1,:)=CONJG(ua(:,tid(jp)))*pb(2)
c     $               +CONJG(dua(:,tid(jp)))*pb(3)
c                  dua2(1,:)=CONJG(ua(:,tid(jp)))*qb(2)
c     $                +CONJG(dua(:,tid(jp)))*qb(3)

                  ua2(1,:)=ua(:,tid(jp))*pb(2)
     $                 +dua(:,tid(jp))*pb(3)
                  dua2(1,:)=ua(:,tid(jp))*qb(2)
     $                 +dua(:,tid(jp))*qb(3)
                  term2=MATMUL(dua2,imat)*qb(ip)
     $                 +MATMUL(ua2,vmat)*qb(ip)
     $                 +MATMUL(ua2,umat)*pb(ip)
                  cell%mat(jp,:,ep,ip)=cell%mat(jp,:,ep,ip)+w*term2(1,:)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     fill matrix for (ua*Ti,Luaj)
c-----------------------------------------------------------------------
            DO ip=1,3
               DO jp=1,3
c                  ua2(1,:)=CONJG(ua(:,tid(ip)))*pb(2)
c     $                    +CONJG(dua(:,tid(ip)))*pb(3)
c                  dua2(1,:)=CONJG(ua(:,tid(ip)))*qb(2)
c     $                     +CONJG(dua(:,tid(ip)))*qb(3)
                  ua2(1,:)=ua(:,tid(ip))*pb(2)+dua(:,tid(ip))*pb(3)
                  dua2(1,:)=ua(:,tid(ip))*qb(2)+dua(:,tid(ip))*qb(3)
                  ua1(:,1)=ua(:,tid(jp))*pb(2)+dua(:,tid(jp))*pb(3)
                  dua1(:,1)=ua(:,tid(jp))*qb(2)+dua(:,tid(jp))*qb(3)
                  
                  term=MATMUL(dua2,MATMUL(imat,dua1))
     $                +MATMUL(ua2,MATMUL(vmat,dua1))
     $                +MATMUL(ua2,MATMUL(umat,ua1))
                  cell%mat(ip,jp,ep,ep)=cell%mat(ip,jp,ep,ep)
     $                 +w*term(1,1)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     fill rhs for ext
c-----------------------------------------------------------------------
            ua1(:,1)=uax(:,tid(4))
            dua1(:,1)=duax(:,tid(4)) 
            DO ip=0,cell%np   
               term1=qb(ip)*MATMUL(imat,dua1)+pb(ip)*MATMUL(vmat,dua1)
     $              +pb(ip)*MATMUL(umat,ua1)
               cell%rhs(:,ip)=cell%rhs(:,ip)-term1(:,1)*w
            ENDDO
            DO ip=1,3
c               ua2(1,:)=CONJG(ua(:,tid(ip)))*pb(2)
c     $                 +CONJG(dua(:,tid(ip)))*pb(3)
c               dua2(1,:)=CONJG(ua(:,tid(ip)))*qb(2)
c     $                  +CONJG(dua(:,tid(ip)))*qb(3)
               ua2(1,:)=ua(:,tid(ip))*pb(2)+dua(:,tid(ip))*pb(3)
               dua2(1,:)=ua(:,tid(ip))*qb(2)+dua(:,tid(ip))*qb(3)

               term=MATMUL(dua2,MATMUL(imat,dua1))
     $             +MATMUL(ua2,MATMUL(vmat,dua1))
     $             +MATMUL(ua2,MATMUL(umat,ua1))   
               cell%rhs(ip,ep)=cell%rhs(ip,ep)-term(1,1)*w
            ENDDO
         ENDDO                 
c-----------------------------------------------------------------------
c     surface term.
c-----------------------------------------------------------------------
c     note: here assume the surface term in the resonant element
c     consider the resonant element is integrated by parts.
      ELSEIF (cell%etype == "ext1"  .OR.  cell%etype == "ext2") THEN
c-----------------------------------------------------------------------
c     compute rhs for ext1 and ext2.
c-----------------------------------------------------------------------
         DO iq=0,quad%np
            x0=(cell%x(2)+cell%x(1))/2
            dx=(cell%x(2)-cell%x(1))/2
            x=x0+dx*quad%node(iq)
            w=dx*quad%weight(iq)
            CALL deltac_get_uv(x,imat,umat,vmat)
            CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)
            pb=hermite%pb
            qb=hermite%qb
            SELECT CASE (cell%etype)
            CASE ("ext1")
               CALL deltac_get_ua(x,ua)
               CALL deltac_get_dua(x,dua)    
               ua1(:,1)=ua(:,tid(4))
               dua1(:,1)=dua(:,tid(4))
            CASE ("ext2")
               IF (x>=0) THEN
                  CALL deltac_get_ua(cell%x(2),ua)
                  CALL deltac_get_dua(cell%x(2),dua)  
                  ua1(:,1)=ua(:,tid(4))*pb(2)+dua(:,tid(4))*pb(3)
                  dua1(:,1)=ua(:,tid(4))*qb(2)+dua(:,tid(4))*qb(3)
               ELSE
                  CALL deltac_get_ua(cell%x(1),ua)
                  CALL deltac_get_dua(cell%x(1),dua)  
                  ua1(:,1)=ua(:,tid(4))*pb(0)+dua(:,tid(4))*pb(1)
                  dua1(:,1)=ua(:,tid(4))*qb(0)+dua(:,tid(4))*qb(1)
               ENDIF
            END SELECT
            DO ip=0,cell%np   
               term1=qb(ip)*MATMUL(imat,dua1)+pb(ip)*MATMUL(vmat,dua1)
     $              +pb(ip)*MATMUL(umat,ua1)
               cell%rhs(:,ip)=cell%rhs(:,ip)-term1(:,1)*w
            ENDDO
         ENDDO
      ELSE
         CALL program_stop("unknown element type")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_extension    
c-----------------------------------------------------------------------
c     subprogram 19. deltac_lsode_int.
c     computes resonant quadratures with lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_lsode_int(u_res,u_h1,u_h2)

      COMPLEX(r8), DIMENSION(3,4), INTENT(INOUT) :: u_res
      COMPLEX(r8), DIMENSION(mpert,0:1,1:4),INTENT(INOUT) :: u_h1
      COMPLEX(r8), DIMENSION(mpert,1:3,0:1),INTENT(INOUT) :: u_h2  
      INTEGER :: ip
      INTEGER :: nstep
      INTEGER :: neq,itol,itask,istate,iopt,lrw,liw,jac,mf,istep
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8) :: rtol,tmp
      REAL(r8) :: x0,x1,x,dx,t,dt
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: u,atol
c-----------------------------------------------------------------------
c     set integrator parameters.
c-----------------------------------------------------------------------
      IF (method) THEN
         neq=2*(3*4)
      ELSE
         neq=2*(3*4+mpert*2*4+mpert*3*2)
      ENDIF
      itol=2
      mf=10
      liw=20
      lrw=20+16*neq     
      istate=1
      itask=5
      iopt=1
      istep=0
      nstep=5000
      rtol=deltac_tol
c-----------------------------------------------------------------------
c     initialize independent variable.
c-----------------------------------------------------------------------
      x0=cell%x(1)
      x1=cell%x_lsode
      IF (x0 > x1) CALL program_stop
     $     ("gal_lsode_int: left resonant element too small.")
      dx=x1-x0
      x=x0
c-----------------------------------------------------------------------
c     allocate and initialize work arrays.
c-----------------------------------------------------------------------
      ALLOCATE(u(neq/2),iwork(liw),rwork(lrw),atol(neq/2))
c-----------------------------------------------------------------------
c    set atol for each component.
c-----------------------------------------------------------------------
      CALL deltac_lsode_der(neq,x0,u,u)
      u=u*deltac_tol
      DO ip=1,neq/2
c         tmp=MAX( DABS(REAL(u(ip),8)),DABS(IMAG(u(ip),8)) )
         tmp=CDABS(u(ip))
         IF (tmp <deltac_tol) tmp=deltac_tol
         atol(ip)=CMPLX(tmp,tmp,8)
      ENDDO
      u=0
      iwork=0
      rwork=0
      rwork(1)=x1
c-----------------------------------------------------------------------
c     advance and diagnose differential equations.
c-----------------------------------------------------------------------
      DO
         t=(x-x0)/dx
         dt=rwork(11)/dx
         IF(ABS(t-1) < deltac_tol  .OR.  istep >= nstep)EXIT
         istep=istep+1
         CALL lsode(deltac_lsode_der,neq,u,x,x1,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
      IF (istep >= nstep) THEN 
         WRITE (*,*)"Warning: LSODE exceeds nstep."
         READ (*,*)
      ENDIF
c-----------------------------------------------------------------------
c     output and deallocate.
c-----------------------------------------------------------------------
      IF (method) THEN
         u_res=RESHAPE(u,SHAPE(u_res))
      ELSE
         u_res=RESHAPE(u(1:12),SHAPE(u_res))
         u_h1=RESHAPE(u(13:36),SHAPE(u_h1))
         u_h2=RESHAPE(u(37:neq/2),SHAPE(u_h2))
      ENDIF
      DEALLOCATE(u,iwork,rwork,atol)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_lsode_int
c-----------------------------------------------------------------------
c     subprogram 20. deltac_lsode_der.
c     derivatives for resonant quadratures with lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_lsode_der(neq,x,u,du)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq/2), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(neq/2), INTENT(OUT) :: du

      INTEGER :: ip,jp
      COMPLEX(r8), DIMENSION(3,4) :: du_res
      COMPLEX(r8), DIMENSION(1,1) :: term
      COMPLEX(r8), DIMENSION(1,3) :: ua1,dua1
      COMPLEX(r8), DIMENSION(3,1) :: ua2,dua2
      COMPLEX(r8), DIMENSION(3,6) :: ua,dua
      COMPLEX(r8), DIMENSION(mpert,0:1,1:4) :: du_h1
      COMPLEX(r8), DIMENSION(mpert,1:3,0:1) :: du_h2
      
      COMPLEX(r8), DIMENSION(3,3) :: umat,vmat,imat
      TYPE(hermite_type) :: hermite      
c-----------------------------------------------------------------------
c     compute basis functions
c-----------------------------------------------------------------------
      CALL deltac_get_ua(x,ua)
      CALL deltac_get_dua(x,dua)
      CALL deltac_get_uv(x,imat,umat,vmat)
      CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)

c-----------------------------------------------------------------------
c     nonresonant terms.
c-----------------------------------------------------------------------
      du_h1=0
      du_h2=0
      du_res=0
      IF(.NOT.method) THEN
         DO ip=0,1
            DO jp=1,4
               ua2(:,1)=ua(:,tid(jp))
               dua2(:,1)=dua(:,tid(jp))
               du_h1(:,ip:ip,jp)=hermite%qb(ip)*MATMUL(imat,dua2)
     $              +hermite%pb(ip)*MATMUL(vmat,dua2)
     $              +hermite%pb(ip)*MATMUL(umat,ua2)
            ENDDO
         ENDDO
         DO ip=1,3
            DO jp=0,1
c              ua1(1,:)=CONJG(ua(:,tid(ip)))
c              dua1(1,:)=CONJG(dua(:,tid(ip)))    
               ua1(1,:)=ua(:,tid(ip))
               dua1(1,:)=dua(:,tid(ip))         
               du_h2(:,ip,jp:jp)=MATMUL(dua1,imat)*hermite%qb(jp)
     $              +MATMUL(ua1,vmat)*hermite%qb(jp)
     $              +MATMUL(ua1,umat)*hermite%pb(jp)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     resonant terms.
c-----------------------------------------------------------------------
      DO ip=1,3
c     ua1(1,:)=CONJG(ua(:,tid(ip)))
c     dua1(1,:)=CONJG(dua(:,tid(ip)))
         ua1(1,:)=ua(:,tid(ip))
         dua1(1,:)=dua(:,tid(ip))
         DO jp=1,4
            ua2(:,1)=ua(:,tid(jp))
            dua2(:,1)=dua(:,tid(jp))
            term=MATMUL(dua1,MATMUL(imat,dua2))
     $           +MATMUL(ua1,MATMUL(vmat,dua2))
     $           +MATMUL(ua1,MATMUL(umat,ua2))
            du_res(ip,jp)=term(1,1)
         ENDDO
      ENDDO
      IF (gal_method=="resonant".AND.noexp) THEN
         DO ip=1,3
            DO jp=1,4
               IF (ip==1.AND.jp==1) CYCLE
               IF (ip==1.AND.jp==4) CYCLE
               du_res(ip,jp)=0
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     output.
c-----------------------------------------------------------------------
      IF (method) THEN
         du=RESHAPE((du_res),SHAPE(du))
      ELSE
         du=RESHAPE((/du_res,du_h1,du_h2/),SHAPE(du))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_lsode_der
c-----------------------------------------------------------------------
c     subprogram 21. deltac_get_d2ua.
c     get 2nd derivative value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_get_d2ua(x,d2ua)

      REAL(r8),INTENT(IN) :: x
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: d2ua
      
      INTEGER i,j,t,l
      REAL(r8),DIMENSION(2) :: p
      REAL(r8) :: x2,x3
      COMPLEX(r8) :: xp,xs,xj,q
      COMPLEX(r8),DIMENSION(2) :: s

      q=in%q
      d2ua=0.0
      p=asp%p
      s=asp%s
      x2=x*x
      x3=x2*x
      DO i=1,2
         xp=x**p(i)
         IF (rescale) THEN
            xs=EXP((-x2+xmax*xmax)/(SQRT(in%q)*2.0) )
     $           *((x/xmax)**asp%s(i))
         ELSE
            xs=EXP(-x2/(SQRT(in%q)*2.0) )*(x**asp%s(i))
         ENDIF
         DO j=0, order_pow
            l=j+1
            xj=x**(-2.0-2.0*j)            
            t=2+i
            d2ua(1,t)=d2ua(1,t)+asp%v(1,t,l)*(p(i)+1-2*j)*(p(i)-2*j)*xj
            d2ua(2:3,t)=d2ua(2:3,t)
     $           +asp%v(2:3,t,l)*(p(i)-2*j)*(p(i)-2*j-1)*xj
         ENDDO
         DO j=0, order_exp
            l=j+1   
            t=4+i
            xj=x**(-2.0*j)
            d2ua(1,t)=d2ua(1,t)
     $           +asp%v(1,t,l)*(
     $           (s(i)-1-2*j)*(s(i)-2-2*j)/x3
     $           -1.0/SQRT(q)*(s(i)-2*j)/x
     $           -1.0/SQRT(q)*(s(i)-1-2*j)/x
     $           +1.0/q*x
     $           )*xj
            d2ua(2:3,t)=d2ua(2:3,t)
     $           +asp%v(2:3,t,l)*(
     $           (s(i)-2*j)*(s(i)-1-2*j)/x2
     $           -1.0/SQRT(q)*(s(i)-2*j+1)
     $           -1.0/SQRT(q)*(-2*j+s(i))
     $           +1.0/q*x*x
     $           )*xj
         ENDDO
         t=2+i
         d2ua(:,t)=xp*d2ua(:,t)
         d2ua(1,t)=x*d2ua(1,t)
         t=4+i
         d2ua(:,t)=xs*d2ua(:,t)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_get_d2ua      
c-----------------------------------------------------------------------
c     subprogram 22. deltac_ua_diagnose.
c     diagnoses asymptotic solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_ua_diagnose

      INTEGER :: ix,t
      REAL(r8) ::dx,x,logx,logx0,logx1
      COMPLEX(r8), DIMENSION(mpert,6) :: ua,dua,d2ua
      COMPLEX(r8), DIMENSION(mpert,mpert) :: imat,vmat,umat
      COMPLEX(r8), DIMENSION(3,6) :: res
c-----------------------------------------------------------------------
c     open output files and set pointer.
c-----------------------------------------------------------------------
      OPEN(UNIT=deltac_out_unit,FILE="ua.out",STATUS="REPLACE")
      OPEN(UNIT=deltac_bin_unit,FILE="ua.bin",STATUS="REPLACE",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      logx0=log10(x0_ua)
      logx1=log10(x1_ua)
      dx=(logx1-logx0)/nx_ua
c-----------------------------------------------------------------------
c     start loops over isol, iqty, and iz.
c-----------------------------------------------------------------------
      DO ix=0,nx_ua
c-----------------------------------------------------------------------
c     compute positions.
c-----------------------------------------------------------------------
         x=10**(logx0+dx*ix)
c-----------------------------------------------------------------------
c     compute power series solutions.
c-----------------------------------------------------------------------
         CALL deltac_get_ua(x,ua)
         CALL deltac_get_dua(x,dua)
         CALL deltac_get_d2ua(x,d2ua)
         CALL deltac_get_uv(x,imat,umat,vmat)
c         d2ua(2,:)=d2ua(2,:)*in%q*in%q
c         d2ua(3,:)=d2ua(3,:)*in%q
c-----------------------------------------------------------------------
c     compute u''-vu'-v
c-----------------------------------------------------------------------
         DO t=3,6
            res(:,t)=MATMUL(imat,d2ua(:,t))
     $           -MATMUL(vmat,dua(:,t))-MATMUL(umat,ua(:,t))
         ENDDO
c-----------------------------------------------------------------------
c     write graphical output.
c-----------------------------------------------------------------------
         logx=log10(x)
         WRITE(deltac_bin_unit)REAL(x,4),REAL(logx,4),
     $        floored_log(res(1,outt)),
     $        floored_log(res(2,outt)),
     $        floored_log(res(3,outt))
         WRITE(deltac_out_unit,'(5(e30.10))')x,logx,
     $        floored_log(res(1,outt)),
     $        floored_log(res(2,outt)),
     $        floored_log(res(3,outt))
      ENDDO
      WRITE(deltac_bin_unit)
c-----------------------------------------------------------------------
c     close files and restore msol.
c-----------------------------------------------------------------------
      CLOSE(deltac_out_unit)
      CLOSE(deltac_bin_unit)
      CALL program_stop("ua diagnostic.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_ua_diagnose
c-----------------------------------------------------------------------
c     subprogram 23. deltac_get_solution.
c     get the solutions of inner layer.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_get_solution(x,gal,icell,isol,sol)

      INTEGER, INTENT(IN) :: isol
      INTEGER, INTENT(INOUT) :: icell
      REAL(r8), INTENT(IN) :: x
      TYPE(gal_type), INTENT(IN) :: gal
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: sol
      
      INTEGER :: i,ndelta
      REAL(r8) :: xext
      COMPLEX(r8), DIMENSION(mpert) :: delta
      REAL(r8), DIMENSION(2) :: epb
      REAL(r8), DIMENSION(0:np) :: pb
      COMPLEX(r8),DIMENSION(mpert,0:np) :: u 
      COMPLEX(r8), DIMENSION(mpert,6) :: ua,uaext,duaext
      TYPE(cell_type), POINTER :: cell
      TYPE(hermite_type) :: hermite
c-----------------------------------------------------------------------
c     find the cell and interval contain x.
c-----------------------------------------------------------------------
      IF (x < xmin .OR. x > xmax) THEN
         CALL program_stop("x is out of range.")
      ENDIF
      DO
         cell => gal%intvl%cell(icell)
         IF (x >= cell%x(1) .AND. x <= cell%x(2)) THEN
            EXIT
         ENDIF
         icell=icell+1
         IF (icell > gal%nx) THEN
            icell=1
         ENDIF
      ENDDO
      sol=0.0
c-----------------------------------------------------------------------
c     construct non-resonant solution (normal solution).
c-----------------------------------------------------------------------
      IF (restore_uh) THEN
         CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)
         pb=hermite%pb
         DO i=0,cell%np
            u(:,i)=gal%sol(cell%map(:,i),isol)
            sol=sol+u(:,i)*pb(i)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     construct small and large solutions.
c-----------------------------------------------------------------------
      IF (.NOT.(restore_us .OR. restore_ul)) RETURN
      IF (cell%etype == "none") RETURN
      CALL deltac_hermite(x,cell%x(1),cell%x(2),hermite)
      epb(1)=hermite%pb(2)
      epb(2)=hermite%pb(3)
      xext=cell%x(2)

      IF (cell%etype == "res" .OR. cell%etype == "ext"
     $     .OR. cell%etype == "ext1")
     $     CALL deltac_get_ua(x,ua)
      IF (cell%etype == "ext" .OR. cell%etype == "ext2") THEN
         CALL deltac_get_ua(xext,uaext)
         CALL deltac_get_dua(xext,duaext)
      ENDIF
      IF (restore_us .AND.
     $     (cell%etype == "ext" .OR. cell%etype == "res"))THEN
         IF (gal_method=="resonant".AND.noexp) THEN
            delta(1)=gal%sol(cell%emap(1),isol)
            ndelta=1            
         ELSE
            delta=gal%sol(cell%emap,isol)
            ndelta=3
         ENDIF
         SELECT CASE(cell%etype)
         CASE("res")
            DO i=1,ndelta
               sol=sol+delta(i)*ua(:,tid(i))
            ENDDO
         CASE("ext")
            DO i=1,ndelta
               sol=sol+delta(i)*
     $             (epb(1)*uaext(:,tid(i))
     $            +epb(2)*duaext(:,tid(i)))
            ENDDO 
         END SELECT
      ENDIF
      IF (restore_ul) THEN   
         i=4
         SELECT CASE(cell%etype)
      CASE ("res","ext","ext1")
         sol=sol+ua(:,tid(i))
      CASE("ext2")
         sol=sol
     $        +(epb(1)*uaext(:,tid(i))
     $        +epb(2)*duaext(:,tid(i)))
      END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_get_solution
c-----------------------------------------------------------------------
c     subprogram 24. deltac_output_solution.
c     output the solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_output_solution (gal)

      TYPE(gal_type), INTENT(IN) :: gal
      
      INTEGER :: icell, ip,tot_grids,isol,ipert,m
      CHARACTER(100) :: comp_tittle,tmp
      CHARACTER(100),DIMENSION(3) :: filename
      REAL(r8),PARAMETER :: eps=1e-4
      REAL(r8), DIMENSION(interp_np) :: x
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     generate the interpolated grid.
c-----------------------------------------------------------------------
      tot_grids=interp_np*gal%nx
      x=(1.0/interp_np)*(/(ip,ip=0,interp_np-1)/)
      ip=0
      DO icell=1,gal%nx
         cell => gal%intvl%cell(icell)
         sol%xvar(ip:ip+interp_np-1)=cell%x(1)+(cell%x(2)-cell%x(1))*x
         ip=ip+interp_np
      ENDDO
      sol%xvar(ip)=xmax
      sol%sol=0
c-----------------------------------------------------------------------
c     get solution.
c-----------------------------------------------------------------------
      DO isol=1,msol
         icell=1
         DO ip=0,tot_grids
            CALL deltac_get_solution(sol%xvar(ip),gal,icell,isol,
     $           sol%sol(:,ip,isol))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     normalize solution to match out region.
c-----------------------------------------------------------------------
      sol%xvar=sol%xvar*in%x0/in%v1
      sol%sol=sol%sol*(in%v1/in%x0)**(0.5+in%p1)
c-----------------------------------------------------------------------
c     output solution for ascii.
c-----------------------------------------------------------------------
      IF (deltac_out_sol) THEN
         DO isol=1,msol
            WRITE (filename(2),"(I4)")in%ising
            filename(2)=ADJUSTL(filename(2))
            WRITE (filename(3),"(I4)")isol
            filename(3)=ADJUSTL(filename(3))
            WRITE(filename(1),*) 'deltac_sur_'
     $           //TRIM(filename(2))//'_sol_'//TRIM(filename(3))//'.out'
            OPEN(UNIT=deltac_out_unit,FILE=TRIM(filename(1)),
     $           STATUS="REPLACE")
            WRITE (deltac_out_unit,10) 'xvar'
            DO m=1,3
               WRITE (tmp,"(I4)") m
               tmp=ADJUSTL(tmp)
               WRITE (comp_tittle,*) 'REAL(F',TRIM(tmp),')'
               WRITE (deltac_out_unit,10) TRIM(comp_tittle)
               WRITE (comp_tittle,*) 'IMAG(F',TRIM(tmp),')'
               WRITE (deltac_out_unit,10) TRIM(comp_tittle)
            ENDDO
 10         FORMAT (1P,A15,$)
            WRITE (deltac_out_unit,*)
            DO ip=0,tot_grids
               WRITE (deltac_out_unit,20) sol%xvar(ip)
 20            FORMAT (1P,E15.5,$)
               DO ipert=1,mpert
                  WRITE (deltac_out_unit,20)
     $                 REAL(sol%sol(ipert,ip,isol)),
     $                 IMAG(sol%sol(ipert,ip,isol))
               ENDDO
               WRITE (deltac_out_unit,*)
            ENDDO
            CLOSE(deltac_out_unit)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     xdraw output.
c-----------------------------------------------------------------------
      IF (deltac_bin_sol) THEN
         DO isol=1,msol
            WRITE (filename(2),"(I4)")in%ising
            filename(2)=ADJUSTL(filename(2))
            WRITE (filename(3),"(I4)")isol
            filename(3)=ADJUSTL(filename(3))
            WRITE(filename(1),*) 'deltac_sur_'
     $           //TRIM(filename(2))//'_sol_'//TRIM(filename(3))//'.bin'
            OPEN(UNIT=deltac_bin_unit,FILE=TRIM(filename(1)),
     $           STATUS="REPLACE",FORM="UNFORMATTED")
            DO ip=0,tot_grids                
               WRITE (deltac_bin_unit) REAL(sol%xvar(ip),4),
     $              REAL(sol%sol(:,ip,isol),4)
            ENDDO
            WRITE(deltac_bin_unit)
            CLOSE(deltac_bin_unit)
         ENDDO
         
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltac_output_solution  
c-----------------------------------------------------------------------
c     subprogram 25. deltac_read_parameters.
c     read the parameters for deltac run.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltac_read_parameters(filename)

      CHARACTER(*), INTENT(IN) :: filename

      NAMELIST/deltac_list/
     $     deltac_bin_sol,deltac_out_sol,dx1dx2_flag,restore_uh,
     $     restore_us,restore_ul,interp_np,nx,nq,order_pow,order_exp,
     $     cutoff,deltac_tol,pfac,xmax,diagnose_res,outt,fulldomain,
     $     nx_ua,x0_ua,x1_ua,tid,dx1dx2_flag,dx1,dx2,gal_method,
     $     side,xfac,rescale,xmax_method,diagnose_params,noexp
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE=TRIM(filename),STATUS="OLD")
      READ(in_unit,NML=deltac_list)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE deltac_read_parameters  
c-----------------------------------------------------------------------
c     subprogram 26. deltac_estimate_zi.
c     estimate zi value for zeroth order term dominant 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE deltac_estimate_zi(zi)
      REAL(r8), INTENT (OUT) :: zi
      REAL(r8) :: tmp1,tmp2,tmp3,tmp4
c-----------------------------------------------------------------------
c     computation.
c----------------------------------------------------------------------- 
      tmp1=ABS(asp%v(2,3,1))
      tmp2=ABS(asp%v(2,3,2))
      tmp3=tmp2/tmp1
      tmp1=ABS(asp%v(2,4,1))
      tmp2=ABS(asp%v(2,4,2))
      tmp4=tmp2/tmp1
      IF (tmp4 > tmp3) tmp3=tmp4
      zi=SQRT(tmp3)*in%x0/in%v1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------      
      END SUBROUTINE deltac_estimate_zi     
      END MODULE deltac_mod
