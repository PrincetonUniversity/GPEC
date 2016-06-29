c-----------------------------------------------------------------------
c     file gal.f.
c     Dewar's Galerkin method for singular modes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. gal_mod.
c     1. gal_alloc.
c     2. gal_dealloc.
c     3. gal_hermite.
c     4. gal_pack.
c     5. gal_make_grid.
c     6. gal_diagnose_grid.
c     7. gal_make_map.
c     8. gal_diagnose_map.
c     9. gal_assemble_rhs.
c     10. gal_assemble_mat.
c     11. gal_lsode_int.
c     12. gal_lsode_der.
c     13. gal_lsode_diagnose.
c     14. gal_get_fkg.
c     15. gal_gauss_quad.
c     16. gal_extension.
c     17. gal_make_arrays.
c     18. gal_solve.
c     19. gal_write_delta.
c     20. gal_matmul.
c     21. gal_make_delta.
c-----------------------------------------------------------------------
c     subprogram 0. gal_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gal_mod
      USE free_mod
      USE jacobi_mod
      USE sing_mod
      IMPLICIT NONE

      TYPE :: hermite2_type
      REAL(r8), DIMENSION(0:3) :: pb,qb
      END TYPE hermite2_type

      TYPE :: cell_type
      CHARACTER(6) :: extra,etype
      INTEGER :: emap
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: map
      REAL(r8), DIMENSION(2) :: x
      COMPLEX(r8) :: erhs,ediag
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: rhs,emat
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat
      END TYPE cell_type

      TYPE :: interval_type
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,dx
      TYPE(cell_type), DIMENSION(:), POINTER :: cell
      END TYPE interval_type

      TYPE :: gal_type
      INTEGER :: nx,nq,ndim,kl,ku,ldab,msol
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      REAL(r8) :: pfac,dx0,dx1,dx2
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: mat,rhs,sol
      TYPE(interval_type), DIMENSION(:), POINTER :: intvl
      TYPE(jacobi_type) :: quad
      TYPE(hermite2_type) :: hermite
      END TYPE gal_type

      LOGICAL, PRIVATE :: diagnose_map=.FALSE.,diagnose_grid=.FALSE.,
     $     diagnose_lsode=.FALSE.,diagnose_integrand=.FALSE.
      CHARACTER(1), DIMENSION(2) :: side=(/"l","r"/)
      CHARACTER(128), PRIVATE :: format1,format2
      INTEGER :: ndiagnose=12
      INTEGER, PRIVATE :: jsing,np=3
      TYPE(cell_type), POINTER, PRIVATE :: cell

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gal_alloc.
c     allocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_alloc(gal,nx,nq,dx0,dx1,dx2,pfac)
      
      TYPE(gal_type), INTENT(OUT) :: gal
      INTEGER, INTENT(IN) :: nx,nq
      REAL(r8), INTENT(IN) :: dx0,dx1,dx2,pfac

      INTEGER :: ix,ising
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     set integers.
c-----------------------------------------------------------------------
      gal%nx=nx
      gal%nq=nq
      gal%dx0=dx0
      gal%dx1=dx1
      gal%dx2=dx2
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(gal%intvl(0:msing))
      CALL jacobi_alloc(gal%quad,nq,.TRUE.,.FALSE.,"gll")
      DO ising=0,msing
         intvl => gal%intvl(ising)
         ALLOCATE(intvl%x(0:nx),intvl%dx(0:nx),intvl%cell(nx))
         DO ix=1,nx
            cell => intvl%cell(ix)
            ALLOCATE(cell%map(mpert,0:np),
     $           cell%mat(mpert,mpert,0:np,0:np))
            IF(cell%extra /= "none")
     $           ALLOCATE(cell%emat(mpert,0:np),cell%rhs(mpert,0:np))
         ENDDO
         CALL gal_make_grid(ising,nx,dx0,dx1,dx2,pfac,intvl)
      ENDDO
      IF(diagnose_grid)CALL gal_diagnose_grid(gal)
c-----------------------------------------------------------------------
c     create and diagnose local-to-global mapping.
c-----------------------------------------------------------------------
      CALL gal_make_map(gal)
      IF(diagnose_map)CALL gal_diagnose_map(gal)
c-----------------------------------------------------------------------
c     allocate global arrays.
c-----------------------------------------------------------------------
      gal%kl=mpert*(np+1)+1
      gal%ku=gal%kl
      gal%ldab=2*gal%kl+gal%ku+1
      gal%msol=2*msing
      ALLOCATE(gal%rhs(gal%ndim,gal%msol),gal%sol(gal%ndim,gal%msol),
     $     gal%mat(gal%ldab,gal%ndim),gal%ipiv(gal%ndim))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_alloc
c-----------------------------------------------------------------------
c     subprogram 2. gal_dealloc.
c     deallocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_dealloc(gal)
      
      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ising,ix
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
            DEALLOCATE(cell%map,cell%mat)
            IF(cell%extra /= "none")DEALLOCATE(cell%emat,cell%rhs)
         ENDDO
         DEALLOCATE(intvl%x,intvl%dx,intvl%cell)
      ENDDO
      DEALLOCATE(gal%rhs,gal%sol,gal%mat,gal%ipiv,gal%intvl)
      CALL jacobi_dealloc(gal%quad)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. gal_hermite.
c     computes hermite cubic basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_hermite(x,x0,x1,hermite)

      REAL(r8),INTENT(IN) :: x,x0,x1
      TYPE(hermite2_type),INTENT(INOUT) :: hermite 
      
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
      hermite%pb(1)=t12*t0
      hermite%pb(2)=t02*(1+2*t1)
      hermite%pb(3)=-t02*t1
c-----------------------------------------------------------------------
c     compute derivative of hermite basis.
c-----------------------------------------------------------------------
      hermite%qb(0)=-6*t0*t1/dx
      hermite%qb(1)=-t1*(3*t0-1)/dx
      hermite%qb(2)=6*t1*t0/dx
      hermite%qb(3)=t0*(3*t1-1)/dx
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_hermite
c-----------------------------------------------------------------------
c     subprogram 4. gal_pack.
c     computes packed grid on (0,1).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION gal_pack(nx,pfac,side) RESULT(x)

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
     $        ("gal_pack: cannot recognize side = "//TRIM(side))
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
      END FUNCTION gal_pack
c-----------------------------------------------------------------------
c     subprogram 5. gal_make_grid.
c     sets up grid in one interval.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_make_grid(ising,nx,dx0,dx1,dx2,pfac,intvl)

      INTEGER :: ising,nx
      REAL(r8), INTENT(IN) :: dx0,dx1,dx2,pfac
      TYPE(interval_type), INTENT(INOUT) :: intvl

      INTEGER :: ixmin,ixmax,ix,mx
      REAL(r8) :: x0,x1,xm,dx,nq1
c-----------------------------------------------------------------------
c     set default extra.
c-----------------------------------------------------------------------
      DO ix=1,nx
         intvl%cell(ix)%extra="none"
         intvl%cell(ix)%etype="none"
      ENDDO
c-----------------------------------------------------------------------
c     set lower bound.
c-----------------------------------------------------------------------
      IF(ising == 0)THEN
         ixmin=0
         intvl%x(ixmin)=psilow
      ELSE
         x0=sing(ising)%psifac
         nq1=ABS(nn*sing(ising)%q1)
         intvl%x(0)=x0+dx0/nq1
         ixmin=0
c$$$         intvl%x(1)=x0+(dx0+dx1)/nq1
c$$$         intvl%x(2)=x0+(dx0+dx1+dx2)/nq1
c$$$         ixmin=2
         intvl%cell(1)%extra="right"
         intvl%cell(2)%extra="right"
         intvl%cell(1)%etype="res"
         intvl%cell(2)%etype="ext"
      ENDIF
c-----------------------------------------------------------------------
c     set upper bound.
c-----------------------------------------------------------------------
      IF(ising == msing)THEN
         ixmax=nx
         intvl%x(ixmax)=psihigh
      ELSE
         x1=sing(ising+1)%psifac
         nq1=ABS(sing(ising+1)%q1)
         intvl%x(nx)=x1-dx0/nq1
         ixmax=nx
c$$$         intvl%x(nx-1)=x1-(dx0+dx1)/nq1
c$$$         intvl%x(nx-2)=x1-(dx0+dx1+dx2)/nq1
c$$$         ixmax=nx-2
         intvl%cell(nx-1)%extra="left"
         intvl%cell(nx)%extra="left"
         intvl%cell(nx-1)%etype="ext"
         intvl%cell(nx)%etype="res"
      ENDIF
c-----------------------------------------------------------------------
c     compute interior packed grid.
c-----------------------------------------------------------------------
      xm=(intvl%x(ixmax)+intvl%x(ixmin))/2
      dx=(intvl%x(ixmax)-intvl%x(ixmin))/2
      mx=(ixmax-ixmin)/2
c$$$      intvl%x(ixmin:ixmax)=xm+dx*pack(mx,pfac,"both")
      intvl%dx(0)=0
      intvl%dx(1:nx)=intvl%x(1:nx)-intvl%x(0:nx-1)
      DO ix=1,nx
         intvl%cell(ix)%x=(/intvl%x(ix-1),intvl%x(ix)/)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_make_grid
c-----------------------------------------------------------------------
c     subprogram 6. gal_diagnose_grid.
c     diagnoses grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_diagnose_grid(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ising,ix,ipert
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: q,singfac,psi
      TYPE(interval_type), POINTER :: intvl
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"ix",8x,"x",12x,"dx",10x,"singfac"/)
 20   FORMAT(i6,1p,3e14.6)
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      mvec=mlow+(/(ipert,ipert=0,mpert-1)/)
      OPEN(UNIT=gal_out_unit,FILE="grid.out",STATUS="UNKNOWN")
      OPEN(UNIT=gal_bin_unit,FILE="grid.bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
      DO ising=0,msing
         intvl => gal%intvl(ising)
         WRITE(gal_out_unit,'(a,i2)')"ising = ",ising
         WRITE(gal_out_unit,10)
         DO ix=0,gal%nx
            psi=intvl%x(ix)
            CALL spline_eval(sq,psi,0)
            q=sq%f(4)
            singfac=MINVAL(ABS(mvec-nn*q))
            WRITE(gal_out_unit,20)ix,intvl%x(ix),intvl%dx(ix),singfac
            WRITE(gal_bin_unit)REAL(intvl%x(ix),4),REAL(one,4),
     $           REAL(intvl%dx(ix),4)
         ENDDO
         WRITE(gal_out_unit,10)
         WRITE(gal_bin_unit)
      ENDDO
      CLOSE(UNIT=gal_out_unit)
      CLOSE(UNIT=gal_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_diagnose_grid
c-----------------------------------------------------------------------
c     subprogram 7. gal_make_map.
c     creates local-to-global mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_make_map(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      LOGICAL :: skip=.FALSE.
      INTEGER :: ising,ix,ip,imap,ipert,emap
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      imap=1
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
c-----------------------------------------------------------------------
c     nonresonant basis functions.
c-----------------------------------------------------------------------
            IF(imap > 1)imap=imap-2*mpert
            DO ip=0,np
               DO ipert=1,mpert
                  cell%map(ipert,ip)=imap
                  imap=imap+1
               ENDDO
               IF(ip == 1 .AND. skip)THEN
                  imap=imap+1
                  skip=.FALSE.
               ENDIF
            ENDDO
c-----------------------------------------------------------------------
c     resonant basis functions, left of singular point.
c-----------------------------------------------------------------------
            SELECT CASE(cell%extra)
            CASE("left")
               SELECT CASE(cell%etype)
               CASE("ext")
                  emap=imap
                  cell%emap=emap
                  skip=.TRUE.
               CASE("res")
                  cell%emap=emap
               END SELECT
c-----------------------------------------------------------------------
c     resonant basis functions, right of singular point.
c-----------------------------------------------------------------------
            CASE("right")
               SELECT CASE(cell%etype)
               CASE("res")
                  emap=imap
                  cell%emap=emap
                  skip=.TRUE.
               CASE("ext")
                  cell%emap=emap
               END SELECT
            END SELECT
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
         ENDDO
      ENDDO
      gal%ndim=imap-1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_make_map
c-----------------------------------------------------------------------
c     subprogram 8. gal_diagnose_map.
c     diagnoses local-to-global mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_diagnose_map(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ising,ix,ip,ipert
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"ip",1x,"ipert",5x,"map"/)
 20   FORMAT(2i6,i8)
 30   FORMAT(/2a6,i8,2x,"*")
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      OPEN(UNIT=debug_unit,FILE="map.out",STATUS="UNKNOWN")
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
            WRITE(debug_unit,'(2(a,i2))')"ising = ",ising,", ix = ",ix
c-----------------------------------------------------------------------
c     write map.
c-----------------------------------------------------------------------
            DO ip=0,np
               WRITE(debug_unit,10)
               DO ipert=1,mpert
                  WRITE(debug_unit,20)ip,ipert,cell%map(ipert,ip)
               ENDDO
            ENDDO
            IF(cell%extra /= "none")WRITE(debug_unit,30)
     $           ADJUSTR(cell%extra),ADJUSTR(cell%etype),cell%emap
            WRITE(debug_unit,10)
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_diagnose_map
c-----------------------------------------------------------------------
c     subprogram 9. gal_assemble_rhs.
c     assembles global rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_assemble_rhs(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      LOGICAL :: first
      INTEGER :: ising,ix,imap,isol
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,*)"start rhs_assemble"
      gal%rhs=0
      isol=0
      first=.TRUE.
      IF(diagnose)WRITE(*,'(a,2i5)')" SHAPE(gal%rhs) = ",SHAPE(gal%rhs)
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
            IF(cell%extra == "none")CYCLE
            IF(first)THEN
               isol=isol+1
               first=.FALSE.
            ELSE
               first=.TRUE.
            ENDIF
            imap=cell%emap
            IF(diagnose)WRITE(*,'(3(a,i2),a,i5)')" ising = ",ising,
     $           ", ix = ",ix,", isol = ",isol,", imap = ",imap
            gal%rhs(imap,isol)=gal%rhs(imap,isol)+cell%erhs
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,*)"finish rhs_ assemble"
      RETURN
      END SUBROUTINE gal_assemble_rhs
c-----------------------------------------------------------------------
c     subprogram 10. gal_assemble_mat.
c     assembles global mat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_assemble_mat(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ising,ix,ip,ipert,i,jp,jpert,j,offset
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals, grid cells, and rows.
c-----------------------------------------------------------------------
      offset=gal%kl+gal%ku+1
      gal%mat=0
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
            DO ip=0,np
               DO ipert=1,mpert
                  i=cell%map(ipert,ip)
c-----------------------------------------------------------------------
c     nonresonant basis functions.
c-----------------------------------------------------------------------
                  DO jp=0,np
                     DO jpert=1,mpert
                        j=cell%map(jpert,jp)
                        gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                       +cell%mat(ipert,jpert,ip,jp)
                     ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c     resonant basis functions.
c-----------------------------------------------------------------------
                  IF(cell%extra /= "none")THEN
                     j=cell%emap
                     gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                    +cell%emat(ipert,ip)
                     gal%mat(offset+j-i,i)=gal%mat(offset+j-i,i)
     $                    +CONJG(cell%emat(ipert,ip))
                     gal%mat(offset,j)=gal%mat(offset,j)+cell%ediag
                  ENDIF
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_assemble_mat
c-----------------------------------------------------------------------
c     subprogram 11. gal_lsode_int.
c     computes resonant quadratures with lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_lsode_int(u_res,u_hermite)

      COMPLEX(r8), DIMENSION(2), INTENT(INOUT) :: u_res
      COMPLEX(r8), DIMENSION(0:np,mpert,2), INTENT(INOUT) :: u_hermite

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(16) :: name
      INTEGER, PARAMETER :: nstep=500
      INTEGER :: neq,itol,itask,istate,iopt,lrw,liw,jac,mf,istep
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), PARAMETER :: tol=1e-6,rtol=tol,atol=tol*tol
      REAL(r8) :: x,x0,x1,dx,t,dt
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: u
c-----------------------------------------------------------------------
c     set integrator parameters.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,'(a,i2,2a)')" start gal_lsode_int with "
     $     //"jsing = ",jsing,", side = ",cell%extra
      neq=4*(mpert*(np+1)+1)
      liw=20
      lrw=22*neq+16
      itol=1
      itask=5
      istate=1
      iopt=1
      mf=10
      istep=0
c-----------------------------------------------------------------------
c     initialize independent variable.
c-----------------------------------------------------------------------
      IF(cell%extra == "left")THEN
         x0=cell%x(1)
         x1=cell%x(2)
      ELSE
         x0=cell%x(2)
         x1=cell%x(1)
      ENDIF
      dx=x1-x0
      x=x0
c-----------------------------------------------------------------------
c     diagnose integrand.
c-----------------------------------------------------------------------
      IF(diagnose_integrand .AND. jsing == uad%ising .AND.
     $     ((cell%extra == "left" .AND. uad%neg) .OR.
     $     (cell%extra == "right" .AND. .NOT. uad%neg)))THEN
         uad%ipert=sing(jsing)%r1(1)
         uad%iqty=1
         uad%ising=jsing
         uad%zlogmax=LOG10(ABS(x0-sing(jsing)%psifac))
         uad%zlogmin=LOG10(ABS(x1-sing(jsing)%psifac))
         CALL sing_ua_diagnose(jsing)
         CALL program_stop
     $        ("gal_lsode_int: abort after integrand_diagnose.")
      ENDIF
c-----------------------------------------------------------------------
c     allocate and initialize work arrays.
c-----------------------------------------------------------------------
      ALLOCATE(u(neq/2),iwork(liw),rwork(lrw))
      u=0
      iwork=0
      rwork=0
      rwork(1)=x1
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(diagnose_lsode)THEN
         WRITE(name,'(a2,i2.2,a1)')"ls",jsing,cell%extra
         OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//".out",
     $        STATUS="UNKNOWN")
         OPEN(UNIT=gal_bin_unit,FILE=TRIM(name)//".bin",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
      ENDIF
c-----------------------------------------------------------------------
c     advance and diagnose differential equations.
c-----------------------------------------------------------------------
      DO
         t=(x-x0)/dx
         dt=rwork(11)/dx
         IF(diagnose_lsode)
     $        CALL gal_lsode_diagnose(neq,istep,x0,x1,x,t,dt,u)
         IF(ABS(t-1) < tol .OR. istep >= nstep)EXIT
         istep=istep+1
         CALL lsode(gal_lsode_der,neq,u,x,x1,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
      IF(cell%extra == "right")u=-u
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(diagnose_lsode)THEN
         CLOSE(UNIT=gal_out_unit)
         CLOSE(UNIT=gal_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     compute output and deallocate.
c-----------------------------------------------------------------------
      u_res=u(1:2)
      u_hermite=RESHAPE(u(3:neq),SHAPE(u_hermite))
      DEALLOCATE(u,iwork,rwork)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_lsode_int
c-----------------------------------------------------------------------
c     subprogram 12. gal_lsode_der.
c     derivatives for resonant quadratures with lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_lsode_der(neq,x,u,du)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq/2), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(neq/2), INTENT(OUT) :: du

      INTEGER :: ip,ipert0,i
      INTEGER, DIMENSION(2) :: isol
      COMPLEX(r8), DIMENSION(2) :: du_res
      COMPLEX(r8), DIMENSION(mpert,0:np,2) :: du_hermite
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua,dua
      COMPLEX(r8), DIMENSION(mpert,msol) :: term1,term2,term3,matvec
      TYPE(hermite2_type) :: hermite
c-----------------------------------------------------------------------
c     set isol.
c-----------------------------------------------------------------------
      ipert0=NINT(nn*sing(jsing)%q)-mlow+1
      isol=(/ipert0,ipert0+mpert/)
c-----------------------------------------------------------------------
c     compute basis functions and matvec.
c-----------------------------------------------------------------------
      CALL sing_get_ua(jsing,x,ua)
      CALL sing_get_dua(jsing,x,dua)
      CALL sing_matvec(x,ua,dua,term1,term2,term3,matvec)
      CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
c-----------------------------------------------------------------------
c     resonant terms.
c-----------------------------------------------------------------------
      DO i=1,2
         du_res(i)=SUM(CONJG(ua(:,isol(2),1))*matvec(:,isol(i)))
      ENDDO
c-----------------------------------------------------------------------
c     nonresonant terms.
c-----------------------------------------------------------------------
      DO ip=0,np
         du_hermite(:,ip,1:2)=hermite%pb(ip)*matvec(:,isol)
      ENDDO
c-----------------------------------------------------------------------
c     output.
c-----------------------------------------------------------------------
      du=RESHAPE((/du_res,du_hermite/),SHAPE(du))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_lsode_der
c-----------------------------------------------------------------------
c     subprogram 13. gal_lsode_diagnose.
c     derivatives for quadratures with lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_lsode_diagnose(neq,istep,x0,x1,x,t,dt,u)

      INTEGER, INTENT(IN) :: neq,istep
      REAL(r8), INTENT(IN) :: x0,x1,x,t,dt
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: u

      INTEGER :: i
      REAL(r8) :: q,singfac,singlog
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/1x,"istep",8x,"x",11x,"t",10x,"dt",6x,"singfac",1x,',
     $     i2.2,'(4x,"re u",i2.2,5x,"im u",i2.2,1x)/)')
 20   FORMAT('(i6,1p,e14.6,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      IF(istep == 0)THEN
         WRITE(format1,10)ndiagnose
         WRITE(format2,20)2*ndiagnose+4
         WRITE(gal_out_unit,'(a,i4,1p,2(a,e13.6))')
     $        "neq = ",neq,", x0 = ",x0,", x1 = ",x1
      ENDIF
c-----------------------------------------------------------------------
c     compute singfac.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,x,0)
      q=sq%f(4)
      singfac=sing(jsing)%m-nn*q
      singlog=LOG10(ABS(singfac))
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      IF(istep == 0 .OR. (MOD(istep,10) == 1 .AND. istep > 2))
     $     WRITE(gal_out_unit,format1)(i,i,i=1,ndiagnose)
      WRITE(gal_out_unit,format2)istep,x,t,dt,singfac,u(1:ndiagnose)
      WRITE(gal_bin_unit)REAL(t,4),REAL(dt,4),
     $     (REAL(u(i),4),REAL(IMAG(u(i)),4),i=1,ndiagnose),
     $     REAL(singfac,4),REAL(singlog,4)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_lsode_diagnose
c-----------------------------------------------------------------------
c     subprogram 14. gal_get_fkg.
c     computes nonresonant quadratures with Gaussian quadrature.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_get_fkg(psifac,f,k,g)

      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: f,k,g

      INTEGER :: ipert,jpert,iqty
      REAL(r8) :: q
      REAL(r8), DIMENSION(mpert) :: singfac
c-----------------------------------------------------------------------
c     compute q and singfac.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psifac,0)
      q=sq%f(4)
      singfac=mlow-nn*q+(/(ipert,ipert=0,mpert-1)/)
c-----------------------------------------------------------------------
c     compute f.
c-----------------------------------------------------------------------
      CALL cspline_eval(fmats_gal,psifac,0)
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            f(ipert,jpert)=fmats_gal%f(iqty)
     $          *singfac(ipert)*singfac(jpert)
            IF(ipert /= jpert)THEN
               f(jpert,ipert)=CONJG(f(ipert,jpert))
            ENDIF
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute k.
c-----------------------------------------------------------------------
      CALL cspline_eval(kmats,psifac,0)
      iqty=1
      DO jpert=1,mpert
         DO ipert=MAX(1,jpert-mband),MIN(mpert,jpert+mband)
            k(ipert,jpert)=kmats%f(iqty)*singfac(ipert)
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute g.
c-----------------------------------------------------------------------
      CALL cspline_eval(gmats,psifac,0)
      iqty=1
      DO jpert=1,mpert
         DO ipert=jpert,MIN(mpert,jpert+mband)
            g(ipert,jpert)=gmats%f(iqty)
            IF(ipert /= jpert)THEN
               g(jpert,ipert)=CONJG(g(ipert,jpert))
            ENDIF
            iqty=iqty+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_get_fkg
c-----------------------------------------------------------------------
c     subprogram 15. gal_gauss_quad.
c     computes nonresonant matrix elements in one cell.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_gauss_quad(ising,quad)

      INTEGER, INTENT(IN) :: ising
      TYPE(jacobi_type) :: quad

      INTEGER :: ip,jp,iq,ipert,jpert
      REAL(r8) :: x0,dx,x,w
      REAL(r8), DIMENSION(0:np) :: pb,qb
      COMPLEX(r8), DIMENSION(mpert,mpert) :: f,k,g
      TYPE(hermite2_type) :: hermite
c-----------------------------------------------------------------------
c     start loop over quadrature points.
c-----------------------------------------------------------------------
      cell%mat=0
      DO iq=0,quad%np
         x0=(cell%x(2)+cell%x(1))/2
         dx=(cell%x(2)-cell%x(1))/2
         x=x0+dx*quad%node(iq)
         w=dx*quad%weight(iq)
         CALL gal_get_fkg(x,f,k,g)
         CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
         pb=hermite%pb
         qb=hermite%qb
c-----------------------------------------------------------------------
c     compute gaussian quadratures.
c-----------------------------------------------------------------------
         DO ipert=1,mpert
            DO ip=0,np
               DO jpert=1,mpert
                  DO jp=0,np
                     cell%mat(ipert,jpert,ip,jp)
     $                    =cell%mat(ipert,jpert,ip,jp)
     $                    +f(ipert,jpert)*qb(ip)*qb(jp)
     $                    +k(ipert,jpert)*qb(ip)*pb(jp)
     $                    +CONJG(k(jpert,ipert))*pb(ip)*qb(jp)
     $                    +g(ipert,jpert)*pb(ip)*pb(jp)
                  ENDDO
               ENDDO
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
      END SUBROUTINE gal_gauss_quad
c-----------------------------------------------------------------------
c     subprogram 16. gal_extension.
c     computes extension terms of matrix and rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_extension(ising)

      INTEGER, INTENT(IN) :: ising

      INTEGER :: ipert,isol,jp,jsing
      REAL(r8) :: x
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua,dua
c-----------------------------------------------------------------------
c     set values at cell boundary.
c-----------------------------------------------------------------------
      IF(cell%extra == "left")THEN
         jp=2
         jsing=ising+1
         x=cell%x(2)
      ELSEIF(cell%extra == "right")THEN
         jp=0
         jsing=ising
         x=cell%x(1)
      ENDIF
c-----------------------------------------------------------------------
c     get asymptotic solutions and derivatives.
c-----------------------------------------------------------------------
      CALL sing_get_ua(jsing,x,ua)
      CALL sing_get_dua(jsing,x,dua)
c-----------------------------------------------------------------------
c     compute emat and erhs.
c-----------------------------------------------------------------------
      isol=NINT(nn*sing(jsing)%q)-mlow+1
      DO ipert=1,mpert
         cell%rhs=cell%rhs
     $        -cell%mat(:,ipert,:,jp)*ua(ipert,isol,1)
     $        -cell%mat(:,ipert,:,jp+1)*dua(ipert,isol,1)
         cell%emat=cell%emat
     $        +cell%mat(:,ipert,:,jp)*ua(ipert,isol+mpert,1)
     $        +cell%mat(:,ipert,:,jp+1)*dua(ipert,isol+mpert,1)
      ENDDO
c-----------------------------------------------------------------------
c     compute ediag.
c-----------------------------------------------------------------------
      cell%ediag=0
      DO ipert=1,mpert
         cell%ediag=cell%ediag
     $        +CONJG(ua(ipert,isol+mpert,1))*cell%emat(ipert,jp)
     $        +CONJG(dua(ipert,isol+mpert,1))*cell%emat(ipert,jp+1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_extension
c-----------------------------------------------------------------------
c     subprogram 17. gal_make_arrays.
c     computes matrix and rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_make_arrays(gal)

      TYPE(gal_type), INTENT(INOUT) :: gal

      INTEGER :: ising,ix
      COMPLEX(r8), DIMENSION(2) :: u_res
      COMPLEX(r8), DIMENSION(mpert,0:np,2) :: u_hermite
      TYPE(interval_type), POINTER :: intvl
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      WRITE(*,*)"galerkin make arrays"
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
c-----------------------------------------------------------------------
c     nonresonant arrays.
c-----------------------------------------------------------------------
            CALL gal_gauss_quad(ising,gal%quad)
c-----------------------------------------------------------------------
c     resonant arrrays.
c-----------------------------------------------------------------------
            SELECT CASE(cell%etype)
            CASE("res")
               SELECT CASE(cell%extra)
               CASE("left")
                  jsing=ising+1
               CASE("right")
                  jsing=ising
               END SELECT
               CALL gal_lsode_int(u_res,u_hermite)
               cell%erhs=-u_res(1)
               cell%ediag=u_res(2)
               cell%emat=u_hermite(:,:,2)
            CASE("ext")
               CALL gal_extension(ising)
            END SELECT
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     assemble matrix and rhs.
c-----------------------------------------------------------------------
      CALL gal_assemble_mat(gal)
      CALL gal_assemble_rhs(gal)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_make_arrays
c-----------------------------------------------------------------------
c     subprogram 18. gal_solve.
c     solves for Delta' matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_solve

      INTEGER :: nx,nq,info,isol,jsol,imap
      REAL(r8) :: dx0,dx1,dx2,pfac
      COMPLEX(r8), DIMENSION(2*msing,2*msing) :: delta
      TYPE(gal_type) :: gal

      NAMELIST/gal_input/nx,nq,dx0,dx1,dx2,pfac,diagnose_map,
     $     diagnose_grid,diagnose_lsode,ndiagnose,diagnose_integrand
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      WRITE(*,*)"start galerkin"
      OPEN(UNIT=in_unit,FILE="rdcon.in",STATUS="OLD")
      READ(UNIT=in_unit,NML=gal_input)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     allocate and compute arrays.
c-----------------------------------------------------------------------
      msol=2*mpert
      CALL gal_alloc(gal,nx,nq,dx0,dx1,dx2,pfac)
      CALL gal_make_arrays(gal)
c-----------------------------------------------------------------------
c     solve.
c-----------------------------------------------------------------------
      WRITE(*,*)"galerkin matrix factorization"
      CALL zgbtrf(gal%ndim,gal%ndim,gal%kl,gal%ku,gal%mat,gal%ldab,
     $     gal%ipiv,info)
      WRITE(*,*)"galerkin matrix solve"
      CALL zgbtrs("N",gal%ndim,gal%kl,gal%ku,gal%msol,gal%mat,gal%ldab,
     $     gal%ipiv,gal%sol,gal%ndim,info )
c-----------------------------------------------------------------------
c     compute and write delta from small resonant coefficients.
c-----------------------------------------------------------------------
      isol=0
      DO jsol=1,gal%msol
         IF(cell%etype == "res")THEN
            isol=isol+1
            imap=cell%emap
            delta(isol,jsol)=gal%rhs(imap,isol)
         ENDIF
      ENDDO
      CALL gal_write_delta(delta,"delta1")
c-----------------------------------------------------------------------
c     deallocate arrays and finish.
c-----------------------------------------------------------------------
      CALL gal_dealloc(gal)
      WRITE(*,*)"finish galerkin"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_solve
c-----------------------------------------------------------------------
c     subprogram 19. gal_write_delta.
c     writes delta matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_write_delta(delta,name)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: delta
      CHARACTER(*), INTENT(IN) :: name

      CHARACTER(128) :: format1,format2
      INTEGER :: iside,ising,i
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/1x,"i\j",',i2.2,'(4x,"re ",i1,a1,6x,"im ",i1,a1,2x)/)')
 20   FORMAT('(i3,a1,1p,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      WRITE(format1,10)2*msing
      WRITE(format2,20)4*msing
c-----------------------------------------------------------------------
c     write delta.
c-----------------------------------------------------------------------
      OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//".out",STATUS="UNKNOWN")
      WRITE(gal_out_unit,'(a)')" Delta matrix:"
      WRITE(gal_out_unit,format1)((ising,side(iside),ising,side(iside),
     $     iside=1,2),ising=1,msing)
      i=0
      DO ising=1,msing
         DO iside=1,2
            i=i+1
            WRITE(gal_out_unit,format2)ising,side(iside),delta(i,:)
         ENDDO
      ENDDO
      WRITE(gal_out_unit,format1)((ising,side(iside),ising,side(iside),
     $     iside=1,2),ising=1,msing)
      CLOSE(UNIT=gal_out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE gal_write_delta
c-----------------------------------------------------------------------
c     subprogram 20. gal_matmul.
c     computes product of banded matrix matrix times vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION gal_matmul(gal,u) RESULT(v)

      TYPE(gal_type), INTENT(IN) :: gal
      COMPLEX(r8), DIMENSION(gal%ndim), INTENT(IN) :: u
      COMPLEX(r8), DIMENSION(gal%ndim) :: v

      COMPLEX(r8), PARAMETER :: one=1,zero=0
c-----------------------------------------------------------------------
c     compute matrix-vector product.
c-----------------------------------------------------------------------
      CALL zgbmv('N',gal%ndim,gal%ndim,gal%kl,gal%kl,one,gal%mat,
     $     gal%ldab,u,1,zero,v,1)
c--------------------------------g---------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END FUNCTION gal_matmul
c-----------------------------------------------------------------------
c     subprogram 21. gal_make_delta.
c     computes symmetrized matrix elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_make_delta(gal,delta)

      TYPE(gal_type), INTENT(IN) :: gal
      COMPLEX(r8), DIMENSION(2*msing,2*msing), INTENT(OUT) :: delta

      INTEGER :: isol,jsol
c-----------------------------------------------------------------------
c     compute symmetrized matrix elements.
c-----------------------------------------------------------------------
      DO isol=1,gal%msol
         DO jsol=1,gal%msol
            delta(isol,jsol)
     $           =SUM(gal%sol(:,jsol)*gal_matmul(gal,gal%sol(:,isol)))
     $           -SUM(gal%sol(:,jsol)*gal%rhs(:,isol))
     $           -CONJG(SUM(gal%sol(:,jsol)*gal%rhs(:,isol)))
     $           +SUM(gal%rhs(:,jsol)*gal%rhs(:,jsol))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE gal_make_delta
      END MODULE gal_mod
