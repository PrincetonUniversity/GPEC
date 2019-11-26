c-----------------------------------------------------------------------
c     program inpso.f.
c     gjt power series solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inpso_mod.
c     1. inpso_init.
c     2. inpso_get_ua.
c     3. inpso_get_dua.
c     4. inpso_get_d2ua.
c     5. inpso_get_uv.
c     6. inpso_ua_diagnose.
c     7. inpso_xmax.
c     8. inpso_dealloc.
c     9. inpso_delta.
c-----------------------------------------------------------------------
c     subprogram 0. inpso_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inpso_mod
      USE jacobi_mod
      USE inps_mod
      USE deltar_mod, ONLY : resist_type

      IMPLICIT NONE

      
c      TYPE :: resist_type
c      INTEGER :: ising
c      REAL(r8) :: e,f,g,h,k,m,taua,taur,v1
c      END TYPE resist_type
      
      TYPE :: asymp_type
      REAL(r8), DIMENSION(2) :: p
      COMPLEX(r8), DIMENSION(2):: s
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: v
      END TYPE asymp_type

      TYPE :: inner_type
      INTEGER :: ising
      REAL(r8) :: e,f,g,h,k,m,taua,taur,v1
      REAL(r8) :: di,dr,p1,sfac
      COMPLEX(r8) :: eig,q,x0
      COMPLEX(r8),DIMENSION(2) :: deltar
      END TYPE inner_type

      LOGICAL :: rescale=.TRUE.,dx1dx2_flag=.false.
      LOGICAL :: grid_diagnose=.FALSE.
      CHARACTER(4) :: inps_type="gjt"
      CHARACTER(10) :: gal_method="normal"
      INTEGER, DIMENSION(2) :: rp=(/3,4/)
      INTEGER :: order_pow=10,order_exp=3,fulldomain=0,
     $     mpert=3,kmax=8,nx_ua=100,outt=3
      REAL(r8) :: inps_xfac
      REAL(r8), DIMENSION(0:2) :: inps_eps=(/1e-2,5e-7,1e-7/),inps_xmaxx
      REAL(r8) :: xmax=1,xfac=1
      REAL(r8) :: x0_ua=0.01,x1_ua=1
      REAL(r8) :: dx1,dx2
      TYPE(inner_type) :: in
      TYPE(asymp_type) :: asp
      TYPE(resist_type_inps) :: rt_in

      CONTAINS    
c-----------------------------------------------------------------------
c     subprogram 1. inpso_init.
c     initialize the asymptotic solutions at large x.
c     coefficients of two power series and two small exponential  
c     solutions are solved.
c     GJT Appendix A.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_init
      
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
c     initialize inps.
c-----------------------------------------------------------------------
      SELECT CASE(inps_type)
      CASE("inps")
         rt_in%e=in%e
         rt_in%f=in%f
         rt_in%h=in%h
         rt_in%g=in%g
         rt_in%k=in%k
         rt_in%q=in%q
         CALL inps_init(rt_in,.FALSE.)
         CALL inps_split(kmax,.FALSE.)
c-----------------------------------------------------------------------
c     initialize the parameters (f,h,k,g,h,q) and the coefficients v
c-----------------------------------------------------------------------
      CASE("gjt")
         e=in%e
         f=in%f
         k=in%k
         g=in%g
         h=in%h
         q=in%q
         in%dr=in%e+in%f+in%h*in%h
         in%di=in%e+in%f+in%h-0.25
         in%p1=SQRT(-in%di)
c-----------------------------------------------------------------------
c     p power of large and small solutions (T_3,4).
c-----------------------------------------------------------------------
         order=order_pow
         IF (order < order_exp) order=order_exp
         ALLOCATE(asp%v(3,6,order+1))
         asp%v=0
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
     $              -h*(p(i)-2*j+2)*asp%v(3,t,l-1)
               matb(2)=-q*q*(p(i)-2*j+2)*(p(i)-2*j+1)*asp%v(2,t,l-1)
               matb(3)=h*k*q*q*(p(i)-2*j+3)*asp%v(1,t,l-1)
     $              -(g-e*k)*q*q*asp%v(2,t,l-1)
     $              +(g+f*k)*q*q*asp%v(3,t,l-1)
               IF (j > 1) matb(3)=matb(3)
     $              -q*(p(i)-2*j+4)*(p(i)-2*j+3)*asp%v(3,t,l-2)
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
     $           asp%v(2:3,t,l),2,info)
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
     $              *asp%v(1,t,l-2)
               matb(2)=-h*(s(i)+1.0-2.0*j)*asp%v(1,t,l-1)
               matb(2)=matb(2)-q*q*(s(i)-2*j+2)*(s(i)-2*j+1)
     $              *asp%v(2,t,l-1)
               matb(3)=q*q*k*h*(s(i)+1.0-2.0*j)
     $              *asp%v(1,t,l-1)
     $              -q*(s(i)-2*j+2)*(s(i)-2*j+1)
     $              *asp%v(3,t,l-1)
               mat=mata
               asp%v(:,t,l)=matb
               CALL zgesv(3,1,mat,3,ipiv,asp%v(:,t,l),3,info)
            ENDDO
         ENDDO     
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("inpso_init: cannot recognize inps_type = "
     $        //TRIM(inps_type))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
         RETURN      
         END SUBROUTINE inpso_init
c-----------------------------------------------------------------------
c     subprogram 2. inpso_get_ua.
c     get the value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_get_ua(xin,ua)

      REAL(r8),INTENT(IN) :: xin
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: ua
      COMPLEX(r8), DIMENSION(6,2) :: inps_ua1
      INTEGER :: i,j,t,l
      REAL(r8) :: x,x2
      COMPLEX(r8) :: xp,xs,xj
c-----------------------------------------------------------------------
c     change sign of x.
c-----------------------------------------------------------------------
      IF (xin < 0) THEN
         x=-xin
      ELSE
         x=xin
      ENDIF
      ua=0
c-----------------------------------------------------------------------
c     new method.
c-----------------------------------------------------------------------
      SELECT CASE(inps_type)
      CASE("inps")
         CALL inps_ua(x,inps_ua1)
         ua(1,rp)=inps_ua1(1,:)/x
         ua(2:3,rp)=inps_ua1(2:3,:)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      CASE("gjt")
         x2=x*x
         DO i=1,2
            xp=x**asp%p(i)
            IF (rescale) THEN
               xs=EXP((-x2+xmax*xmax)/(SQRT(in%q)*2.0) )
     $              *((x/xmax)**asp%s(i))
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
      END SELECT
c-----------------------------------------------------------------------
c     extend to negtive x with even parity.
c-----------------------------------------------------------------------
      IF (xin < 0 .AND. fulldomain == 0)ua(1,:)=-ua(1,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_get_ua
c-----------------------------------------------------------------------
c     subprogram 3. inpso_get_dua.
c     get the derivative value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_get_dua(xin,dua)

      REAL(r8),INTENT(IN) :: xin
      REAL(r8) :: x
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: dua
      COMPLEX(r8), DIMENSION(6,2) :: inps_ua1
      INTEGER i,j,t,l
      REAL(r8),DIMENSION(2) :: p
      REAL(r8) :: x2
      COMPLEX(r8) :: xp,xs,xj,q
      COMPLEX(r8),DIMENSION(2) :: s
c-----------------------------------------------------------------------
c     change sign of x.
c-----------------------------------------------------------------------
      IF (xin < 0) THEN
         x=-xin
      ELSE
         x=xin
      ENDIF
      dua=0
c-----------------------------------------------------------------------
c     new method.
c-----------------------------------------------------------------------
      SELECT CASE(inps_type)
      CASE("inps")
         CALL inps_ua(x,inps_ua1)
         dua(1,rp)=inps_ua1(4,:)-inps_ua1(1,:)/x**2
         dua(2:3,rp)=inps_ua1(5:6,:)*x
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      CASE("gjt")
         q=in%q
         p=asp%p
         s=asp%s
         x2=x*x
         DO i=1,2
            xp=x**p(i)
            IF (rescale) THEN
               xs=EXP((-x2+xmax*xmax)/(SQRT(in%q)*2.0) )
     $              *((x/xmax)**asp%s(i))
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
     $              +asp%v(1,t,l)*((s(i)-1-2*j)/x2-1.0/SQRT(q))*xj
               dua(2:3,t)=dua(2:3,t)
     $              +asp%v(2:3,t,l)*((s(i)-2*j)/x-1.0/SQRT(q)*x)*xj
            ENDDO
            t=2+i
            dua(:,t)=xp*dua(:,t)
            dua(1,t)=x*dua(1,t)
            t=4+i
            dua(:,t)=xs*dua(:,t)
         ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     extend to negtive x with even parity.
c-----------------------------------------------------------------------
      IF(xin < 0 .AND. fulldomain == 1)dua(2:3,:)=-dua(2:3,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_get_dua
c-----------------------------------------------------------------------
c     subprogram 4. inpso_get_d2ua.
c     get 2nd derivative value of the asymptotic solutions at large x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_get_d2ua(x,d2ua)

      REAL(r8),INTENT(IN) :: x
      COMPLEX(r8),DIMENSION(3,6),INTENT(OUT) :: d2ua
      
      INTEGER i,j,t,l
      REAL(r8),DIMENSION(2) :: p
      REAL(r8) :: x2,x3
      COMPLEX(r8), DIMENSION(2) :: psi,dpsi
      COMPLEX(r8), DIMENSION(6,2) :: inps_ua1,inps_dua1
      COMPLEX(r8) :: xp,xs,xj,q
      COMPLEX(r8),DIMENSION(2) :: s
c-----------------------------------------------------------------------
c     new method.
c-----------------------------------------------------------------------
      IF(inps_type == "inps")THEN
         CALL inps_ua(x,inps_ua1,inps_dua1)
         psi=inps_ua1(1,:)/x
         dpsi=inps_ua1(4,:)-inps_ua1(1,:)/x**2
         d2ua(1,rp)=inps_dua1(4,:)-dpsi/x+psi/x**2
         d2ua(2:3,rp)=inps_dua1(5:6,:)*x+inps_ua1(5:6,:)
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      q=in%q
      d2ua=0.0
      p=asp%p
      s=asp%s
      x2=x*x
      x3=x2*x
c-----------------------------------------------------------------------
c     compute power series.
c-----------------------------------------------------------------------
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
            d2ua(1,t)=d2ua(1,t)+asp%v(1,t,l)*(p(i)+1-2*j)*(p(i)-2*j)
     $           *xj
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
      END SUBROUTINE inpso_get_d2ua      
c-----------------------------------------------------------------------
c     subprogram 5. inpso_get_uv.
c     computes U, V matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_get_uv(x,imat,umat,vmat)

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
     $     -in%q*x, x*x/in%q,   -(in%g-in%k*in%e)*in%q,
     $     zero,-(in%e+in%f)/(in%q*in%q), x*x/in%q+(in%g+in%k*in%f)*in%q
     $     /),SHAPE(umat))

      
      vmat=RESHAPE((/  zero,   -in%h/(in%q*in%q),  in%h*in%k*in%q,   
     $     zero,            zero,   zero,
     $     dcmplx(in%h),    zero,   zero/)
     $     ,SHAPE(vmat))
      
      umat(2,:)=umat(2,:)*in%q*in%q
      umat(3,:)=umat(3,:)*in%q    
      vmat(2,:)=vmat(2,:)*in%q*in%q
      vmat(3,:)=vmat(3,:)*in%q    
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_get_uv
c-----------------------------------------------------------------------
c     subprogram 6. inpso_ua_diagnose.
c     diagnoses asymptotic solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_ua_diagnose

      INTEGER :: ix,t
      REAL(r8) ::dx,x,logx,logx0,logx1
      COMPLEX(r8), DIMENSION(mpert,6) :: ua,dua,d2ua
      COMPLEX(r8), DIMENSION(mpert,mpert) :: imat,vmat,umat
      COMPLEX(r8), DIMENSION(3,6) :: res
c-----------------------------------------------------------------------
c     open output files and set pointer.
c-----------------------------------------------------------------------
      OPEN(UNIT=inpso_out_unit,FILE="ua.out",STATUS="REPLACE")
      OPEN(UNIT=inpso_bin_unit,FILE="ua.bin",STATUS="REPLACE",
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
         CALL inpso_get_ua(x,ua)
         CALL inpso_get_dua(x,dua)
         CALL inpso_get_d2ua(x,d2ua)
         CALL inpso_get_uv(x,imat,umat,vmat)
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
         WRITE(inpso_bin_unit)REAL(x,4),REAL(logx,4),
     $        mylog(res(1,outt)),
     $        mylog(res(2,outt)),
     $        mylog(res(3,outt))
         WRITE(inpso_out_unit,'(5(e30.10))')x,logx,
     $        mylog(res(1,outt)),
     $        mylog(res(2,outt)),
     $        mylog(res(3,outt))
      ENDDO
      WRITE(inpso_bin_unit)
c-----------------------------------------------------------------------
c     close files and restore msol.
c-----------------------------------------------------------------------
      CLOSE(inpso_out_unit)
      CLOSE(inpso_bin_unit)
      CALL program_stop("ua diagnostic.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_ua_diagnose
c-----------------------------------------------------------------------
c     subprogram 7. inpso_xmax.
c     computes xmax.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_xmax(xmax)

      REAL(r8), INTENT(OUT) :: xmax
c-----------------------------------------------------------------------
c     compute xmax for inps.
c-----------------------------------------------------------------------
      SELECT CASE(inps_type)
      CASE("inps")
         dx1dx2_flag=.TRUE.
         gal_method="resonant"
         CALL inps_xmax(inps_eps,inps_xmaxx)
         inps_xmaxx(1)=inps_xmaxx(1)*inps_xfac
         inps_xmaxx(2)=inps_xmaxx(2)*inps_xfac
         xmax=inps_xmaxx(2)
         dx1=inps_xmaxx(2)-inps_xmaxx(1)
         dx2=dx1/2
         IF(grid_diagnose)WRITE(*,'(a,3es10.3)')
     $        " inps_xmaxx = ",inps_xmaxx
c-----------------------------------------------------------------------
c     compute xmax for inps.
c-----------------------------------------------------------------------
      CASE("gjt")
         xmax=MAXVAL(SQRT(ABS(in%q)**.5*(ABS(asp%s)+1))*xfac)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_xmax
c-----------------------------------------------------------------------
c     subprogram 8. inpso_dealloc.
c     deallocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_dealloc
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      SELECT CASE(inps_type)
      CASE("inps")
         CALL inps_dealloc
      CASE("gjt")
         DEALLOCATE(asp%v)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_dealloc
c-----------------------------------------------------------------------
c     subprogram 9. inpso_delta.
c     computes error in differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE inpso_delta(x,delta)
      
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(2), INTENT(OUT) :: delta

      INTEGER :: j
      REAL(r8), DIMENSION(2,0:3) :: norm
      COMPLEX(r8), DIMENSION(3,6) :: ua6,dua6,d2ua6
      COMPLEX(r8), DIMENSION(3,2) :: ua,dua,d2ua
      COMPLEX(r8), DIMENSION(3,3) :: imat,umat,vmat
      COMPLEX(r8), DIMENSION(3,2,0:3) :: matvec
c-----------------------------------------------------------------------
c     call subroutines.
c-----------------------------------------------------------------------
      CALL inpso_get_ua(x,ua6)
      CALL inpso_get_dua(x,dua6)
      CALL inpso_get_d2ua(x,d2ua6)
      CALL inpso_get_uv(x,imat,umat,vmat)
c-----------------------------------------------------------------------
c     extract power-like solutions.
c-----------------------------------------------------------------------
      ua=ua6(:,rp)
      dua=dua6(:,rp)
      d2ua=d2ua6(:,rp)
c-----------------------------------------------------------------------
c     compute matvec.
c-----------------------------------------------------------------------
      matvec(:,:,1)=MATMUL(imat,d2ua)
      matvec(:,:,2)=-MATMUL(vmat,dua)
      matvec(:,:,3)=-MATMUL(umat,ua)
      matvec(:,:,0)=SUM(matvec(:,:,1:3),3)
c-----------------------------------------------------------------------
c     compute norm and delta.
c-----------------------------------------------------------------------
      DO j=1,2
         norm(j,:)=MAXVAL(ABS(matvec(:,j,:)),1)
         delta(j)=norm(j,0)/MAXVAL(norm(j,1:3))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE inpso_delta
      END MODULE inpso_mod
