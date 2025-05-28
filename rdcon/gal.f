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
c     20. gal_set_boundary.
c     21. gal_diagnose_mat.
c     22. gal_spline_pack.
c     23. gal_read_input.
c     24. gal_write_pest3_data.
c     25. gal_get_solution.
c     26. gal_output_solution.
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
      REAL(r8) :: x_lsode
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
      INTEGER :: nx,nq,ndim,kl,ku,ldab,nsol
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      REAL(r8) :: pfac,dx0,dx1,dx2
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: mat,rhs,sol
      TYPE(interval_type), DIMENSION(:), POINTER :: intvl
      TYPE(jacobi_type) :: quad
      TYPE(hermite2_type) :: hermite
      END TYPE gal_type

      TYPE :: coil_type
      LOGICAL :: rpec_flag=.FALSE.
      INTEGER :: mcoil,m1,m2
      END TYPE coil_type
      
      LOGICAL, PRIVATE :: diagnose_map=.FALSE.,diagnose_grid=.FALSE.,
     $     diagnose_lsode=.FALSE.,diagnose_integrand=.FALSE.,
     $     diagnose_mat=.FALSE.,dx1dx2_flag=.TRUE.,
     $     restore_uh=.TRUE.,restore_us=.TRUE.,restore_ul=.TRUE.,
     $     bin_delmatch=.FALSE.,b_flag=.FALSE.,out_galsol=.FALSE.,
     $     bin_galsol=.FALSE.,bin_coilsol=.FALSE.
      CHARACTER(1), DIMENSION(2) :: side=(/"l","r"/)
      CHARACTER(128), PRIVATE :: format1,format2
      CHARACTER(20) :: solver="cholesky"
      INTEGER :: nx,nq,cutoff=1 
      INTEGER :: ndiagnose=12,interp_np=3
      INTEGER, PRIVATE :: jsing,np=3
      REAL(r8) :: dx0,dx1,dx2,pfac
      REAL(r8) :: gal_tol=1e-10
      LOGICAL :: gal_xmin_flag = .FALSE.
      REAL(r8), DIMENSION(0:2) :: gal_eps_xmin = (/1e-2,1e-6,5e-7/)
      TYPE(cell_type), POINTER, PRIVATE :: cell
      TYPE(coil_type), SAVE, PRIVATE :: coil
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
            cell%emap=0
            cell%map=0.0
            cell%mat=0.0
            cell%x_lsode=0.0
            IF(cell%extra /= "none") THEN
               ALLOCATE(cell%emat(mpert,0:np),cell%rhs(mpert,0:np))
               cell%emat=0.0
               cell%rhs=0.0
            ENDIF
         ENDDO
         CALL gal_make_grid(ising,nx,dx0,dx1,dx2,pfac,intvl)
      ENDDO
c-----------------------------------------------------------------------
c     resistive perturbed equilibrium solution.
c-----------------------------------------------------------------------
      IF (coil%rpec_flag ) THEN
         coil%mcoil=mpert
         coil%m1=2*msing+1
         coil%m2=2*msing+coil%mcoil
      ENDIF
      IF(diagnose_grid)CALL gal_diagnose_grid(gal)
c-----------------------------------------------------------------------
c     create and diagnose local-to-global mapping.
c-----------------------------------------------------------------------
      CALL gal_make_map(gal)
      IF(diagnose_map)CALL gal_diagnose_map(gal)
c-----------------------------------------------------------------------
c     allocate global arrays.
c-----------------------------------------------------------------------
      gal%kl=mpert*(np+1)
      gal%ku=gal%kl
      IF (solver == "LU") THEN
         gal%ldab=2*gal%kl+gal%ku+1
      ELSEIF (solver == "cholesky") THEN
         gal%ldab=gal%kl+1
      ENDIF
      gal%nsol=2*msing
      IF (coil%rpec_flag) gal%nsol=gal%nsol+coil%mcoil
      ALLOCATE(gal%rhs(gal%ndim,gal%nsol),gal%sol(gal%ndim,gal%nsol),
     $     gal%mat(gal%ldab,gal%ndim),gal%ipiv(gal%ndim))
      gal%rhs=0
      gal%sol=0
      gal%mat=0
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
            IF(ALLOCATED(cell%emat))DEALLOCATE(cell%emat)
            IF(ALLOCATED(cell%rhs))DEALLOCATE(cell%rhs)
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
      REAL(r8), DIMENSION(0:2) :: myxmin
      REAL(r8) :: x0,x1,xm,dx,nq1
c-----------------------------------------------------------------------
c     set default extra.
c-----------------------------------------------------------------------
      DO ix=1,nx
         intvl%cell(ix)%extra="none"
         intvl%cell(ix)%etype="none"
      ENDDO
      IF (ising>0 .AND. gal_xmin_flag .AND. sing1_flag) THEN
         CALL sing1_xmin(ising,gal_eps_xmin,myxmin)
         WRITE(*,*) "ising=",ising
         WRITE(*,'(a,3es10.3)') "xmin=",myxmin
         WRITE(*,*) "=================="
      ENDIF
c-----------------------------------------------------------------------
c     set lower bound.
c-----------------------------------------------------------------------
      IF(ising == 0)THEN
         ixmin=0
         intvl%x(ixmin)=psilow
      ELSE
         x0=sing(ising)%psifac
         nq1=ABS(nn*sing(ising)%q1)
         intvl%x(0)=x0
         IF (gal_xmin_flag .AND. sing1_flag) THEN
            CALL sing1_xmin(ising,gal_eps_xmin,myxmin)
            intvl%cell(1)%x_lsode=x0+myxmin(2)/nq1
            intvl%x(1)=x0+myxmin(1)/nq1
            intvl%x(2)=x0+2*myxmin(1)/nq1
            ixmin=2
         ELSE
            intvl%cell(1)%x_lsode=x0+dx0/nq1
            IF (dx1dx2_flag) THEN
               intvl%x(1)=x0+(dx1)/nq1
               intvl%x(2)=x0+(dx1+dx2)/nq1
               ixmin=2
            ELSE
               ixmin=0
            ENDIF
         ENDIF
         intvl%cell(1)%extra="right"
         intvl%cell(2)%extra="right"
         intvl%cell(1)%etype="res"
         intvl%cell(2)%etype="ext"
         DO ix=3,cutoff+1
            intvl%cell(ix)%etype="ext1"
            intvl%cell(ix)%extra="right"
         ENDDO
         IF (ix > nx) THEN
            CALL program_stop("Too many elements include big solution.")
         ELSEIF (intvl%cell(ix)%etype .NE. "none") THEN
            CALL program_stop("Too many elements include big solution.")
         ENDIF
         intvl%cell(ix)%etype="ext2"
         intvl%cell(ix)%extra="right"
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
         intvl%x(nx)=x1
         IF (gal_xmin_flag .AND. sing1_flag) THEN
            CALL sing1_xmin(ising+1,gal_eps_xmin,myxmin)
            intvl%cell(nx)%x_lsode=x1-myxmin(2)/nq1
            intvl%x(nx-1)=x1-myxmin(1)/nq1
            intvl%x(nx-2)=x1-2.0_r8*myxmin(1)/nq1
            ixmax=nx-2
         ELSE
            intvl%cell(nx)%x_lsode=x1-dx0/nq1
            IF (dx1dx2_flag) THEN
               intvl%x(nx-1)=x1-(dx1)/nq1
               intvl%x(nx-2)=x1-(dx1+dx2)/nq1
               ixmax=nx-2
            ELSE
               ixmax=nx
            ENDIF
         ENDIF
         intvl%cell(nx-1)%extra="left"
         intvl%cell(nx)%extra="left"
         intvl%cell(nx-1)%etype="ext"
         intvl%cell(nx)%etype="res"
         DO ix=nx-2,nx-cutoff,-1
            intvl%cell(ix)%etype="ext1"
            intvl%cell(ix)%extra="left"
         ENDDO
         IF (intvl%cell(ix)%etype.NE."none".OR.ix < 1) THEN
            CALL program_stop("Too many elements include big solution.")
         ENDIF
         intvl%cell(ix)%etype="ext2"
         intvl%cell(ix)%extra="left"
      ENDIF
c-----------------------------------------------------------------------
c     compute interior packed grid.
c-----------------------------------------------------------------------
      x0=intvl%x(ixmin)
      x1=intvl%x(ixmax)
      xm=(intvl%x(ixmax)+intvl%x(ixmin))/2.0_r8
      dx=(intvl%x(ixmax)-intvl%x(ixmin))/2.0_r8
      mx=(ixmax-ixmin)/2.0_r8
      intvl%x(ixmin:ixmax)=xm+dx*gal_pack(mx,pfac,"both")
      intvl%x(ixmin)=x0
      intvl%x(ixmax)=x1
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
      INTEGER :: ising,ix,imap,isol,ipert,ip,ipert0,harmo
      TYPE(interval_type), POINTER :: intvl
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,*)"start rhs_assemble"
      isol=0
      IF(diagnose)WRITE(*,'(a,2i5)')" SHAPE(gal%rhs) = ",SHAPE(gal%rhs)
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
            IF(cell%etype == "none") CYCLE
            IF(cell%etype == "ext2".AND.cell%extra == "left") THEN
               isol=isol+1
            ELSEIF(cell%etype == "res".AND.cell%extra == "right") THEN
               isol=isol+1
            ENDIF
            DO ip=0,np
               DO ipert=1,mpert
                  imap=cell%map(ipert,ip)
                  gal%rhs(imap,isol)=gal%rhs(imap,isol)
     $                 +cell%rhs(ipert,ip)
               ENDDO
            ENDDO
            IF (cell%etype == 'ext'.OR.cell%etype == 'res') THEN
               imap=cell%emap
               gal%rhs(imap,isol)=gal%rhs(imap,isol)
     $              +cell%erhs
            ENDIF
            
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,*)"finish rhs_assemble"
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
      IF (solver == "LU") THEN
         offset=gal%kl+gal%ku+1
      ELSEIF (solver == "cholesky") THEN
         offset=1
      ENDIF
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
                        IF (solver == "LU") THEN
                           gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                          +cell%mat(ipert,jpert,ip,jp)
                        ELSEIF (solver == "cholesky") THEN
                           IF (j>i) CYCLE
                           gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                          +cell%mat(ipert,jpert,ip,jp)
                        ENDIF                        
                     ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c     resonant basis functions.
c-----------------------------------------------------------------------
                  IF(cell%extra /= "none")THEN
                     IF (cell%etype == "ext1" .OR. cell%etype == "ext2")
     $                    CYCLE
                     j=cell%emap
                     IF (solver == "LU") THEN
                        gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                       +cell%emat(ipert,ip)
                        gal%mat(offset+j-i,i)=gal%mat(offset+j-i,i)
     $                       +CONJG(cell%emat(ipert,ip))
                     ELSEIF (solver == "cholesky") THEN
                        IF (j<i) THEN
                           gal%mat(offset+i-j,j)=gal%mat(offset+i-j,j)
     $                          +cell%emat(ipert,ip)
                        ELSEIF (j>i)THEN
                           gal%mat(offset+j-i,i)=gal%mat(offset+j-i,i)
     $                          +CONJG(cell%emat(ipert,ip))
                        ELSE
                           STOP "WRONG CASE"
                        ENDIF
                     ENDIF
                  ENDIF
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
               ENDDO
            ENDDO
            IF (cell%extra /= "none") THEN
               IF (cell%etype == "ext1".OR.cell%etype == "ext2") CYCLE
               j=cell%emap
               gal%mat(offset,j)=gal%mat(offset,j)+cell%ediag
            ENDIF
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
      INTEGER :: nstep
      INTEGER :: neq,itol,itask,istate,iopt,lrw,liw,jac,mf,istep
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8) :: rtol,atol
      REAL(r8) :: x,x0,x1,dx,t,dt
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: u
c-----------------------------------------------------------------------
c     set integrator parameters.
c-----------------------------------------------------------------------
      IF(diagnose)WRITE(*,'(a,i2,2a)')" start gal_lsode_int with "
     $     //"jsing = ",jsing,", side = ",cell%extra
      neq=4*(mpert*(np+1)+1)
      itol=1
      mf=10
      liw=20
      lrw=20+16*neq     
      istate=1
      itask=5
      iopt=1
      
      istep=0
      nstep=5000
      rtol=gal_tol
      atol=gal_tol
c-----------------------------------------------------------------------
c     initialize independent variable.
c-----------------------------------------------------------------------
      IF(cell%extra == "left")THEN
         x0=cell%x(1)
         x1=cell%x_lsode
         IF (x0 > x1) CALL program_stop
     $        ("gal_lsode_int: left resonant element too small.")
      ELSE
         x0=cell%x(2)
         x1=cell%x_lsode
         IF (x0 < x1) CALL program_stop
     $        ("gal_lsode_int: right resonant element too small.")
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
         IF(ABS(t-1) < gal_tol .OR. istep >= nstep)EXIT
         istep=istep+1
         CALL lsode(gal_lsode_der,neq,u,x,x1,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
      IF (istep >= nstep) THEN 
         WRITE (*,*)"Warning: LSODE exceeds nstep."
         WRITE (*,*)"  ABS(t-1)=",ABS(t-1)," > gal_tol=",gal_tol
      ENDIF
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
      u_hermite=RESHAPE(u(3:neq/2),SHAPE(u_hermite))
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
      COMPLEX(r8), DIMENSION(mpert,2*mpert) :: term1,term2,term3,matvec
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
      SUBROUTINE gal_gauss_quad(ising,last,quad)

      INTEGER, INTENT(IN) :: ising,last
      TYPE(jacobi_type) :: quad

      INTEGER :: ip,jp,iq,ipert,jpert
      REAL(r8) :: x0,dx,x,w,tmp
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
c     swap values.
c-----------------------------------------------------------------------
         IF (ising == msing .AND. last == 1) THEN
            tmp=pb(2)
            pb(2)=pb(3)
            pb(3)=tmp
            tmp=qb(2)
            qb(2)=qb(3)
            qb(3)=tmp
         ENDIF
c-----------------------------------------------------------------------
c     compute gaussian quadratures.
c-----------------------------------------------------------------------
         DO ipert=1,mpert
            DO ip=0,np
               DO jpert=1,mpert
                  DO jp=0,np
                     cell%mat(ipert,jpert,ip,jp)
     $                    =cell%mat(ipert,jpert,ip,jp)
     $                    +w*(f(ipert,jpert)*qb(ip)*qb(jp)
     $                    +k(ipert,jpert)*qb(ip)*pb(jp)
     $                    +CONJG(k(jpert,ipert))*pb(ip)*qb(jp)
     $                    +g(ipert,jpert)*pb(ip)*pb(jp))
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
      SUBROUTINE gal_extension(ising,quad)

      INTEGER, INTENT(IN) :: ising
      TYPE(jacobi_type) :: quad

      INTEGER :: ipert,isol,jp,jsing,ip,iq
      REAL(r8) :: x,ws,w1,w2,x0,dx,w
      COMPLEX(r8), DIMENSION(mpert) :: u,du,u1,du1,u2,du2,ub,dub
      COMPLEX(r8), DIMENSION(mpert) :: term1,term2
      COMPLEX(r8), DIMENSION(mpert,mpert) :: f,k,g
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua,dua
      TYPE(hermite2_type) :: hermite
c-----------------------------------------------------------------------
c     set values at cell boundary.
c-----------------------------------------------------------------------
      IF(cell%extra == "left")THEN
         jp=2
         jsing=ising+1
         x=cell%x(2)
         ws=-1.0
      ELSEIF(cell%extra == "right")THEN
         jp=0
         jsing=ising
         x=cell%x(1)
         ws=1.0
      ENDIF
c-----------------------------------------------------------------------
c     set values of big solution at cell boundary.
c-----------------------------------------------------------------------
      IF (cell%etype == "ext".OR.cell%etype == "ext1") THEN
         w1=1.0
         w2=1.0
      ELSEIF (cell%etype == "ext2") THEN
         IF (cell%extra == "left") THEN
            w1=0.0
            w2=1.0
         ELSEIF (cell%extra == "right") THEN
            w1=1.0
            w2=0.0
         ENDIF         
      ENDIF
      isol=NINT(nn*sing(jsing)%q)-mlow+1
c-----------------------------------------------------------------------
c     compute emat and ediag.
c-----------------------------------------------------------------------
      IF (cell%etype == "ext") THEN
c-----------------------------------------------------------------------
c     get asymptotic solutions and derivatives.
c-----------------------------------------------------------------------
         CALL sing_get_ua(jsing,x,ua)
         CALL sing_get_dua(jsing,x,dua)
         u=ua(:,isol+mpert,1)
         du=dua(:,isol+mpert,1)  
         cell%emat=0.0    
         cell%ediag=0.0            
         DO ipert=1,mpert
            cell%emat=cell%emat
     $           +cell%mat(:,ipert,:,jp)*u(ipert)
     $           +cell%mat(:,ipert,:,jp+1)*du(ipert)
         ENDDO
         DO ipert=1,mpert
            cell%ediag=cell%ediag
     $           +CONJG(u(ipert))*cell%emat(ipert,jp)
     $           +CONJG(du(ipert))*cell%emat(ipert,jp+1)
         ENDDO
c-----------------------------------------------------------------------
c     surface term.
c-----------------------------------------------------------------------            
         CALL gal_hermite(x,cell%x(1),cell%x(2),hermite) 
         CALL gal_get_fkg(x,f,k,g)
         DO ip=0,np
            cell%emat(:,ip)=cell%emat(:,ip)
     $     +(MATMUL(f,du)+MATMUL(k,u))*hermite%pb(ip)*ws
         ENDDO
         cell%ediag=cell%ediag
     $   +SUM(CONJG(u)*(MATMUL(f,du)+MATMUL(k,u)))*ws
      ENDIF
c-----------------------------------------------------------------------
c     compute rhs.
c-----------------------------------------------------------------------
      CALL sing_get_ua(jsing,cell%x(1),ua)
      CALL sing_get_dua(jsing,cell%x(1),dua)
      u1=ua(:,isol,1)
      du1=dua(:,isol,1)         
      CALL sing_get_ua(jsing,cell%x(2),ua)
      CALL sing_get_dua(jsing,cell%x(2),dua)
      u2=ua(:,isol,1)
      du2=dua(:,isol,1)
      cell%rhs=0.0   
      IF (cell%etype == "ext" .OR. cell%etype == "ext1") THEN
         cell%erhs=0.0
         DO iq=0,quad%np
            x0=(cell%x(2)+cell%x(1))/2
            dx=(cell%x(2)-cell%x(1))/2
            x=x0+dx*quad%node(iq)
            w=dx*quad%weight(iq)
            CALL gal_get_fkg(x,f,k,g)
            CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
            CALL sing_get_ua(jsing,x,ua)
            CALL sing_get_dua(jsing,x,dua)
            ub=ua(:,isol,1)
            dub=dua(:,isol,1)
            term1=MATMUL(f,dub)+MATMUL(k,ub)
            term2=MATMUL(CONJG(TRANSPOSE(k)),dub)+MATMUL(g,ub)
            DO ip=0,np
               cell%rhs(:,ip)=cell%rhs(:,ip)
     $        -(term1*hermite%qb(ip)+term2*hermite%pb(ip))*w
            ENDDO
         ENDDO
         IF (cell%etype == "ext") THEN
            DO ipert=1,mpert
               cell%erhs=cell%erhs
     $              +(CONJG(u(ipert))*cell%rhs(ipert,jp)
     $              +CONJG(du(ipert))*cell%rhs(ipert,jp+1))
            ENDDO
         ENDIF
      ELSEIF (cell%etype == "ext2") THEN
         DO ipert=1,mpert
            cell%rhs=cell%rhs
     $           -w1
     $           *(cell%mat(:,ipert,:,0)*u1(ipert)
     $           +cell%mat(:,ipert,:,1)*du1(ipert))
            cell%rhs=cell%rhs
     $           -w2
     $           *(cell%mat(:,ipert,:,2)*u2(ipert)
     $           +cell%mat(:,ipert,:,3)*du2(ipert))
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     surface term.
c-----------------------------------------------------------------------            
      CALL gal_hermite(cell%x(1),cell%x(1),cell%x(2),hermite) 
      CALL gal_get_fkg(cell%x(1),f,k,g)
      CALL sing_get_ua(jsing,cell%x(1),ua)
      DO ip=0,np
         cell%rhs(:,ip)=cell%rhs(:,ip)
     $  -(MATMUL(f,du1)+MATMUL(k,u1))*hermite%pb(ip)*w1
      ENDDO
      IF (cell%etype == "ext".AND.cell%extra == "right") THEN
         cell%erhs=cell%erhs
     $            -SUM(CONJG(ua(:,isol+mpert,1))
     $            *(MATMUL(f,du1)+MATMUL(k,u1)))
      ENDIF
      CALL gal_hermite(cell%x(2),cell%x(1),cell%x(2),hermite) 
      CALL gal_get_fkg(cell%x(2),f,k,g)
      CALL sing_get_ua(jsing,cell%x(2),ua)      
      DO ip=0,np
         cell%rhs(:,ip)=cell%rhs(:,ip)
     $  +(MATMUL(f,du2)+MATMUL(k,u2))*hermite%pb(ip)*w2
      ENDDO
      IF (cell%etype == "ext".AND.cell%extra == "left") THEN
         cell%erhs=cell%erhs
     $            +SUM(CONJG(ua(:,isol+mpert,1))
     $            *(MATMUL(f,du2)+MATMUL(k,u2)))
      ENDIF

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

      INTEGER :: ising,ix,last
      COMPLEX(r8), DIMENSION(2) :: u_res
      COMPLEX(r8), DIMENSION(mpert,0:np,2) :: u_hermite
      TYPE(interval_type), POINTER :: intvl
c-----------------------------------------------------------------------
c     start loops over intervals and grid cells.
c-----------------------------------------------------------------------
      WRITE(*,*)"Making galerkin arrays"
      DO ising=0,msing
         intvl => gal%intvl(ising)
         DO ix=1,gal%nx
            cell => intvl%cell(ix)
c-----------------------------------------------------------------------
c     nonresonant arrays.
c-----------------------------------------------------------------------
            last=0
            IF (ix == gal%nx) last=1
            CALL gal_gauss_quad(ising,last,gal%quad)
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
               cell%rhs=-u_hermite(:,:,1)
               cell%emat=u_hermite(:,:,2)
            CASE("ext","ext1","ext2")
               CALL gal_extension(ising,gal%quad)
            END SELECT
c-----------------------------------------------------------------------
c     finish loops over intervals and grid cells.
c-----------------------------------------------------------------------
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set boundary condition.
c-----------------------------------------------------------------------
      CALL gal_set_boundary(gal)
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

      INTEGER :: info,isol,jsol,imap,ising,ix,msol_old
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: delta
      TYPE(gal_type) :: gal
c-----------------------------------------------------------------------
c     allocate and compute arrays.
c-----------------------------------------------------------------------
      WRITE(*,*)"Starting Galerkin method"
      msol_old=msol
      msol=2*mpert
      CALL gal_alloc(gal,nx,nq,dx0,dx1,dx2,pfac)
      CALL gal_make_arrays(gal)
c-----------------------------------------------------------------------
c     solve global matrix.
c-----------------------------------------------------------------------
      gal%sol=gal%rhs
      IF (solver == "LU") THEN
         WRITE(*,*)"Performing Galerkin matrix LU factorization"
         CALL zgbtrf(gal%ndim,gal%ndim,gal%kl,gal%ku,gal%mat,gal%ldab,
     $        gal%ipiv,info)
         WRITE(*,*)"Calculating Galerkin matrix solution"
         CALL zgbtrs("N",gal%ndim,gal%kl,gal%ku,gal%nsol,gal%mat,
     $        gal%ldab,gal%ipiv,gal%sol,gal%ndim,info )
      ELSEIF (solver == "cholesky") THEN
         WRITE(*,*)"Performing Galerkin matrix Cholesky factorization"
         CALL zpbtrf('L',gal%ndim,gal%kl,gal%mat,gal%ldab,info)
         WRITE(*,*)"Calculating Galerkin matrix solution"
         CALL zpbtrs('L',gal%ndim,gal%kl,gal%nsol,gal%mat,gal%ldab,
     $        gal%sol,gal%ndim,info)     
      ENDIF
c-----------------------------------------------------------------------
c     compute and write delta from small resonant coefficients.
c-----------------------------------------------------------------------
      ALLOCATE (delta(gal%nsol,2*msing))
      DO isol=1,gal%nsol
         jsol=0
         DO ising=0,msing
            DO ix=1,gal%nx
               cell=>gal%intvl(ising)%cell(ix)
               IF(cell%etype == "res")THEN
                  jsol=jsol+1
                  imap=cell%emap
                  delta(isol,jsol)=gal%sol(imap,isol)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CALL gal_write_delta(delta,"delta_gw")
      CALL gal_write_pest3_data(delta,"pest3_data")
      CALL gal_output_solution(delta,gal)
c-----------------------------------------------------------------------
c     deallocate arrays and finish.
c-----------------------------------------------------------------------
      CALL gal_dealloc(gal)
      DEALLOCATE (delta)
      msol=msol_old
      WRITE(*,*)"Finished Galerkin method"
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

      CHARACTER(128) :: format1,format2,format3
      INTEGER :: iside,ising,i,nsol
      TYPE(sing_type), POINTER :: singp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/1x,"i\j",',i2.2,'(4x,"re ",i1,a1,6x,"im ",i1,a1,2x)/)')
 20   FORMAT('(i3,a1,1p,',i2.2,'e11.3)')
 30   FORMAT('(i3,"c",1p,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      WRITE(format1,10)2*msing
      WRITE(format2,20)4*msing
      WRITE(format3,30)4*msing
c-----------------------------------------------------------------------
c     write delta to ascii file.
c-----------------------------------------------------------------------
      OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//".out",STATUS="UNKNOWN")
      WRITE(gal_out_unit,'(a)')" Delta matrix:"
      IF( msing>0 ) THEN
         WRITE(gal_out_unit,format1)((ising,side(iside),ising,
     $     side(iside),iside=1,2),ising=1,msing)
      ENDIF
      i=0
      DO ising=1,msing
         DO iside=1,2
            i=i+1
            WRITE(gal_out_unit,format2)ising,side(iside),delta(i,:)
         ENDDO
      ENDDO
      IF (coil%rpec_flag) THEN
         DO ising=coil%m1,coil%m2
            i=i+1
            WRITE(gal_out_unit,format3)ising-2*msing,delta(i,:)
         ENDDO
      ENDIF
      IF( msing>0 ) THEN
         WRITE(gal_out_unit,format1)((ising,side(iside),ising,
     $     side(iside),iside=1,2),ising=1,msing)
      ENDIF
      CLOSE(UNIT=gal_out_unit)
c-----------------------------------------------------------------------
c     write binary data for matching.
c-----------------------------------------------------------------------
      OPEN(UNIT=gal_bin_unit,FILE=TRIM(name)//".dat",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
      nsol=SIZE(delta,1)
      WRITE(gal_bin_unit)msing,nsol,coil%rpec_flag
      WRITE(gal_bin_unit)delta
      DO ising=1,msing
         singp => sing(ising)
         WRITE(gal_bin_unit)
     $        singp%restype%e,singp%restype%f,
     $        singp%restype%h,singp%restype%m,
     $        singp%restype%g,singp%restype%k,
     $        singp%restype%eta,singp%restype%rho,
     $        singp%restype%taua,singp%restype%taur,
     $        singp%restype%v1
      ENDDO      
      CLOSE(UNIT=gal_bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE gal_write_delta
c-----------------------------------------------------------------------
c     subprogram 20. gal_set_boundary.
c     set boundary condition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE gal_set_boundary(gal)
      TYPE(gal_type), INTENT(INOUT) :: gal
      COMPLEX(r8), DIMENSION(mpert,mpert) :: wv
      INTEGER :: idx,imap,ipert,egnum
c-----------------------------------------------------------------------
c     set boundary u(0).
c-----------------------------------------------------------------------      
      cell=>gal%intvl(0)%cell(1)
      cell%mat(:,:,0,:)=0.0
      IF (solver == "cholesky") THEN
         cell%mat(:,:,:,0)=0.0
      ENDIF
      DO idx=1,mpert
         cell%mat(idx,idx,0,0)=1.0
      ENDDO    
c-----------------------------------------------------------------------
c     set boundary u(1).
c-----------------------------------------------------------------------      
      cell=>gal%intvl(msing)%cell(gal%nx)
      IF (coil%rpec_flag) THEN
         cell%mat(:,:,3,:)=0.0
         IF (solver == "cholesky") THEN
            cell%mat(:,:,:,3)=0.0
         ENDIF
         DO idx=1,mpert
            cell%mat(idx,idx,3,3)=1.0
         ENDDO
         ipert=1
         DO idx=coil%m1,coil%m2
            imap=cell%map(ipert,3)
            gal%rhs(imap,idx)=1
            ipert=ipert+1
         ENDDO
c         CALL bin_open(gal_bin_unit,'galwt.dat',"OLD",
c     $                 "REWIND","none")
c         READ(gal_bin_unit) galwt

c         DO ipert=1,mpert
c            DO egnum=1,mpert
c               imap=cell%map(ipert,3)
c               gal%rhs(imap,coil%m1+egnum-1)=galwt(ipert,egnum)
c            ENDDO
c         ENDDO
c         WRITE(*,*) "GALWT=",galwt

c          CALL bin_open(gal_bin_unit,'galwt.dat',"REPLACE",
c     $                 "REWIND","none")
c         WRITE(gal_bin_unit) galwt
c         CALL bin_close(gal_bin_unit)
      ELSEIF (vac_flag) THEN
         CALL free_get_wvac
         wv=wvac*psio**2
         cell%mat(:,:,3,3)=cell%mat(:,:,3,3)+wv
      ELSE
         cell%mat(:,:,3,:)=0.0
         IF (solver == "cholesky") THEN
            cell%mat(:,:,:,3)=0.0
         ENDIF
         DO idx=1,mpert
            cell%mat(idx,idx,3,3)=1.0
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN        
      END SUBROUTINE gal_set_boundary
c-----------------------------------------------------------------------
c     subprogram 21. gal_diagnose_mat.
c     diagnose the global matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------            
      SUBROUTINE gal_diagnose_mat(gal,m1,m2,name)
      INTEGER, INTENT(IN) :: m1,m2
      CHARACTER(*), INTENT(IN) :: name
      TYPE(gal_type), INTENT(INOUT) :: gal
      
      INTEGER :: idx,jdx,i,j,offset
      COMPLEX(r8) :: val
      IF (.NOT.diagnose_mat) RETURN
c-----------------------------------------------------------------------
c     output gal%mat and gal%rhs.
c-----------------------------------------------------------------------       
      IF (solver == "cholesky") RETURN
      OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//".out",STATUS="UNKNOWN")
      offset=gal%kl+gal%ku+1
      WRITE(gal_out_unit,'(a)')"global_matrix:"
      i=1
      DO idx=m1,m2
         j=1
         DO jdx=m1,m2
            val=gal%mat(offset+idx-jdx,jdx)
            WRITE (gal_out_unit,10) i,j,REAL(val),IMAG(val),idx,jdx
            j=j+1
         ENDDO
         i=i+1
      ENDDO
10    FORMAT('M ',I5,I5,2(1P,E11.3),I10,I10)
      WRITE(gal_out_unit,*)
      WRITE(gal_out_unit,'(a)')"rhs:"
      i=1
      DO idx=m1,m2
         val=SUM(gal%rhs(idx,:))
         WRITE (gal_out_unit,20) i,REAL(val),IMAG(val),idx
         i=i+1
      ENDDO
20    FORMAT('R ',I5,2(1P,E11.3),I10)
      CLOSE (UNIT=gal_out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN              
      END SUBROUTINE gal_diagnose_mat
c-----------------------------------------------------------------------
c     subprogram 22. gal_spline_pack.
c     sets up grid for spline.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE gal_spline_pack (nx,dx1,dx2,pfac,sx)
      INTEGER,INTENT(IN) :: nx
      REAL(r8),INTENT(IN) :: dx1,dx2,pfac
      REAL(r8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: sx

      INTEGER :: ising,ix
      TYPE(interval_type),DIMENSION(:),ALLOCATABLE :: intvl
c-----------------------------------------------------------------------
c     distribute nodes between two resonant surface.
c-----------------------------------------------------------------------           
      ALLOCATE(intvl(0:msing))
      DO ising=0,msing
         ALLOCATE(intvl(ising)%x(0:nx),
     $            intvl(ising)%dx(0:nx),intvl(ising)%cell(nx))
         CALL gal_make_grid(ising,nx,0.0_8,dx1,dx2,pfac,intvl(ising))
      ENDDO 
      ix=0
      sx=0.0
      DO ising=0,msing
         IF (ising == 0) THEN
            sx(ix:ix+nx-1)=intvl(ising)%x(0:nx-1)
            ix=ix+nx
         ELSEIF (ising == msing) THEN
            sx(ix:ix+nx-1)=intvl(ising)%x(1:nx)
            ix=ix+nx
         ELSE
            sx(ix:ix+nx-2)=intvl(ising)%x(1:nx-1)
            ix=ix+nx-1
         ENDIF
      ENDDO
      DO ising=0,msing
         DEALLOCATE(intvl(ising)%x,
     $              intvl(ising)%dx,intvl(ising)%cell)
      ENDDO
      DEALLOCATE(intvl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE gal_spline_pack
c-----------------------------------------------------------------------
c     subprogram 23. gal_read_input.
c     read inputs of gal_input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE gal_read_input
      NAMELIST/gal_input/nx,nq,dx0,dx1,dx2,pfac,diagnose_map,solver,
     $     diagnose_grid,diagnose_lsode,ndiagnose,diagnose_integrand,
     $     diagnose_mat,gal_tol,dx1dx2_flag,cutoff,prefac,dpsi_intvl,
     $     dpsi1_intvl,gal_xmin_flag,gal_eps_xmin
      NAMELIST /gal_output/interp_np,restore_uh,restore_us,
     $     restore_ul,bin_delmatch,out_galsol,bin_galsol,b_flag,coil,
     $     bin_coilsol
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      READ(UNIT=in_unit,NML=gal_input)
      REWIND(UNIT=in_unit)
      READ(UNIT=in_unit,NML=gal_output)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_read_input   
c-----------------------------------------------------------------------
c     subprogram 24. gal_write_pest3_data.
c     writes delta matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_write_pest3_data(delta,name)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: delta
      CHARACTER(*), INTENT(IN) :: name
      
      INTEGER :: ising,jsing
      COMPLEX(r8), DIMENSION(msing,msing) :: ap,bp,gammap,deltap
c-----------------------------------------------------------------------
c     construct PEST3 matching data.
c-----------------------------------------------------------------------
      ap=0.0
      bp=0.0
      gammap=0.0
      deltap=0.0
      DO ising=1,msing
         DO jsing=1,msing
            ap(ising,jsing)=delta(2*ising,2*jsing)
     $           +delta(2*ising,2*jsing-1)
     $           +delta(2*ising-1,2*jsing)
     $           +delta(2*ising-1,2*jsing-1)
            bp(ising,jsing)=delta(2*ising,2*jsing)
     $           -delta(2*ising,2*jsing-1)
     $           +delta(2*ising-1,2*jsing)
     $           -delta(2*ising-1,2*jsing-1)
            gammap(ising,jsing)=delta(2*ising,2*jsing)
     $           +delta(2*ising,2*jsing-1)
     $           -delta(2*ising-1,2*jsing)
     $           -delta(2*ising-1,2*jsing-1)
            deltap(ising,jsing)=delta(2*ising,2*jsing)
     $           -delta(2*ising,2*jsing-1)
     $           -delta(2*ising-1,2*jsing)
     $           +delta(2*ising-1,2*jsing-1)
         ENDDO
      ENDDO
c      ap=ap*0.5
c      bp=bp*0.5
c      gammap=gammap*0.5
c      deltap=deltap*0.5
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(2x,"i",3x,"j",2x,"Aij",8x,"Bij",8x,"Gammaij",4x,"Deltaij")
 20   FORMAT(i3,i4,1p,e11.3,1p,e11.3,1p,e11.3,1p,e11.3)
c-----------------------------------------------------------------------
c     write PEST3 matching data to ascii file.
c-----------------------------------------------------------------------
      OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//"_re.out",
     $     STATUS="UNKNOWN")
      WRITE(gal_out_unit,'(a)')" PEST3 matching data real part:"
      WRITE(gal_out_unit,10)
      DO ising=1,msing
         DO jsing=1,msing
            WRITE(gal_out_unit,20)ising,jsing,
     $           REAL(ap(ising,jsing)),REAL(bp(ising,jsing)),
     $           REAL(gammap(ising,jsing)),REAL(deltap(ising,jsing))
         ENDDO
      ENDDO
      CLOSE(UNIT=gal_out_unit)      

      OPEN(UNIT=gal_out_unit,FILE=TRIM(name)//"_im.out",
     $     STATUS="UNKNOWN")
      WRITE(gal_out_unit,'(a)')" PEST3 matching data imag part:"
      WRITE(gal_out_unit,10)
      DO ising=1,msing
         DO jsing=1,msing

            WRITE(gal_out_unit,20)ising,jsing,
     $           IMAG(ap(ising,jsing)),IMAG(bp(ising,jsing)),
     $           IMAG(gammap(ising,jsing)),IMAG(deltap(ising,jsing))
         ENDDO
      ENDDO
      CLOSE(UNIT=gal_out_unit)  
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE gal_write_pest3_data
c-----------------------------------------------------------------------
c     subprogram 25. gal_get_solution.
c     writes delta matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_get_solution(x,gal,icell,iintvl,isol,sol,
     $                            deltaij,cut_flag)
      LOGICAL, INTENT(IN) :: cut_flag
      INTEGER, INTENT(IN) :: isol
      INTEGER, INTENT(INOUT) :: icell,iintvl
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: deltaij
      COMPLEX(r8), DIMENSION(mpert), INTENT(INOUT) :: sol
      TYPE(gal_type), INTENT(IN) :: gal
      
      INTEGER :: i,ising,ipert0,jsol
      REAL(r8) :: tmp,xext
      COMPLEX(r8) :: delta
      INTEGER, DIMENSION(mpert,0:np) :: umap
      REAL(r8), DIMENSION(2) :: epb
      REAL(r8), DIMENSION(0:np) :: pb
      COMPLEX(r8),DIMENSION(mpert,0:np) :: u 
      COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua,uaext,duaext
      TYPE(cell_type), POINTER :: cell
      TYPE(hermite2_type) :: hermite
c-----------------------------------------------------------------------
c     find the cell and interval containing x.
c-----------------------------------------------------------------------      
      IF (x.LT.psilow.OR.x.GT.psihigh) THEN
         CALL program_stop("x is out of range.")
      ENDIF
      DO
         cell=>gal%intvl(iintvl)%cell(icell)
         IF (x.GE.cell%x(1).AND.x.LE.cell%x(2)) THEN
            EXIT
         ENDIF
         icell=icell+1
         IF (icell.GT.gal%nx) THEN
            icell=1
            iintvl=iintvl+1
         ENDIF
         IF (iintvl.GT.msing) THEN
            iintvl=0
         ENDIF
      ENDDO
      sol=0.0
c-----------------------------------------------------------------------
c     construct non-resonant solution (normal solution).
c-----------------------------------------------------------------------
      IF (restore_uh) THEN
         umap=cell%map
         CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
         pb=hermite%pb
         IF (iintvl.EQ.msing.AND.icell.EQ.gal%nx) THEN
            tmp=pb(2)
            pb(2)=pb(3)
            pb(3)=tmp   
         ENDIF
         DO i=0,np
            u(:,i)=gal%sol(umap(:,i),isol)
            sol=sol+u(:,i)*pb(i)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     construct small and large solutions.
c-----------------------------------------------------------------------
      IF (.NOT.(restore_us.OR.restore_ul)) RETURN
      IF (cell%etype == "none") RETURN
      SELECT CASE(cell%extra)
      CASE("left")
         ising=iintvl+1
         CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
         epb(1)=hermite%pb(2)
         epb(2)=hermite%pb(3)
         xext=cell%x(2)
      CASE("right")
         ising=iintvl
         CALL gal_hermite(x,cell%x(1),cell%x(2),hermite)
         epb(1)=hermite%pb(0)
         epb(2)=hermite%pb(1)
         xext=cell%x(1)
      END SELECT      
      IF (cell%etype=="res".OR.cell%etype=="ext".OR.cell%etype=="ext1")
     $   CALL sing_get_ua(ising,x,ua)
      IF (cell%etype=="ext".OR.cell%etype=="ext2") THEN
         CALL sing_get_ua(ising,xext,uaext)
         CALL sing_get_dua(ising,xext,duaext)
      ENDIF
      ipert0=NINT(nn*sing(ising)%q)-mlow+1
      IF (restore_us.AND.(cell%etype=="ext".OR.cell%etype=="res")) THEN
         delta=gal%sol(cell%emap,isol)
         SELECT CASE(cell%etype)
         CASE("res")
            sol=sol+delta*ua(:,ipert0+mpert,1)
         CASE("ext")
            sol=sol+delta*
     $     (epb(1)*uaext(:,ipert0+mpert,1)
     $     +epb(2)*duaext(:,ipert0+mpert,1))
         END SELECT
      ENDIF
      IF (restore_ul) THEN
         IF ( (isol==ising*2-1.AND.cell%extra=="left")
     $        .OR.
     $        (isol==ising*2.AND.cell%extra=="right")
     $      ) THEN     
         SELECT CASE(cell%etype)
            CASE ("res","ext","ext1")
               sol=sol+ua(:,ipert0,1)
            CASE("ext2")
               sol=sol
     $            +(epb(1)*uaext(:,ipert0,1)
     $            +epb(2)*duaext(:,ipert0,1))
         END SELECT
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     prepare for uniform solution matching.
c-----------------------------------------------------------------------      
      IF (cut_flag) THEN
         SELECT CASE (cell%extra)
         CASE("left")
            jsol=2*ising-1
         CASE("right")
            jsol=2*ising
         END SELECT
         delta=deltaij(isol,jsol)
         CALL sing_get_ua_cut(ising,x,ua)
         IF (restore_us) sol=sol-delta*ua(:,ipert0+mpert,1)
         IF (restore_ul) THEN
            IF ( (isol==ising*2-1.AND.cell%extra=="left")
     $           .OR.
     $        (isol==ising*2.AND.cell%extra=="right")
     $      ) THEN
               sol=sol-ua(:,ipert0,1)
            ENDIF
         ENDIF 
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------      
      RETURN
      END SUBROUTINE gal_get_solution
c-----------------------------------------------------------------------
c     subprogram 26. gal_output_solution.
c     output the solution
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gal_output_solution (delta,gal)
      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: delta
      TYPE(gal_type), INTENT(IN) :: gal
      
      INTEGER :: icell, iintvl,ip,tot_grids,isol,ising,j,ipert,m
      LOGICAL, DIMENSION(:), ALLOCATABLE :: issing
      CHARACTER(100) :: comp_tittle,tmp
      CHARACTER(100),DIMENSION(2) :: filename
      REAL(r4) :: realsinh,imagsinh
      REAL(r8), PARAMETER :: eps=1e-4
      REAL(r8), DIMENSION(interp_np) :: xvar
      REAL(r8), DIMENSION(mpert) :: singfac      
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psi,q
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xext
      COMPLEX(r8),DIMENSION(:,:,:), ALLOCATABLE :: sol,sol_cut
      TYPE(cell_type), POINTER :: cell
c-----------------------------------------------------------------------
c     check flag
c-----------------------------------------------------------------------
      IF (.NOT.(bin_delmatch.OR.out_galsol.OR.bin_galsol)) RETURN
c-----------------------------------------------------------------------
c     generate the interpolated grid.
c-----------------------------------------------------------------------
      tot_grids=interp_np*(msing+1)*gal%nx
      ALLOCATE (issing(0:tot_grids),psi(0:tot_grids),q(0:tot_grids),
     $          sol(mpert,0:tot_grids,gal%nsol),
     $          sol_cut(mpert,0:tot_grids,gal%nsol))
      ALLOCATE (xext(msing,2))
      xvar=(1.0/interp_np)*(/(ip,ip=0,interp_np-1)/)
      ip=0
      issing=.FALSE.
      DO iintvl=0,msing
         DO icell=1,gal%nx
            cell=>gal%intvl(iintvl)%cell(icell)
            psi(ip:ip+interp_np-1)=cell%x(1)+(cell%x(2)-cell%x(1))*xvar
            IF (cell%extra=="right".AND.cell%etype=="res") 
     $         issing(ip)=.TRUE.
            ip=ip+interp_np
         ENDDO
      ENDDO
      psi(ip)=psihigh
      sol=0
      sol_cut=0
      DO ip=0,tot_grids
         CALL spline_eval(sq,psi(ip),0)
         q(ip)=sq%f(4)
      ENDDO
c-----------------------------------------------------------------------
c     get solution.
c-----------------------------------------------------------------------       
      DO isol=1,gal%nsol
         iintvl=0
         icell=1
         DO ip=0,tot_grids
            IF (.NOT.issing(ip)) THEN
               CALL gal_get_solution(psi(ip),gal,icell,iintvl,isol,
     $                               sol(:,ip,isol),delta,.FALSE.)
               CALL gal_get_solution(psi(ip),gal,icell,iintvl,isol,
     $                               sol_cut(:,ip,isol),delta,.TRUE.)

            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write outer region solution for match.
c-----------------------------------------------------------------------
      IF (bin_delmatch) THEN
         CALL bin_open(gal_bin_unit,'gal_solution.dat',"UNKNOWN",
     $                 "REWIND","none")
         WRITE (gal_bin_unit) gal%nsol,msing,mpert,mhigh,mlow,nn,
     $                        psio,tot_grids
         DO isol=1,msing
            WRITE(gal_bin_unit) sing(isol)%psifac,sing(isol)%q
         ENDDO 
         WRITE (gal_bin_unit) psi,issing,q
         DO isol=1,gal%nsol
            WRITE (gal_bin_unit) sol(:,:,isol)
         ENDDO
         CALL bin_close(gal_bin_unit)
c-----------------------------------------------------------------------
c     find large solution truncate region.
c-----------------------------------------------------------------------
          DO iintvl=0,msing
            DO icell=1,gal%nx
               cell=>gal%intvl(iintvl)%cell(icell)
               IF (cell%etype=="ext2") THEN
                  IF (cell%extra=="left") THEN
                     xext(iintvl+1,1) = cell%x(2)
                  ENDIF
                  IF (cell%extra=="right") THEN
                     xext(iintvl,2) = cell%x(1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         CALL bin_open(gal_bin_unit,'gal_solution_cut.dat',"UNKNOWN",
     $                 "REWIND","none")
         WRITE (gal_bin_unit) xext
         DO isol=1,gal%nsol
            WRITE (gal_bin_unit) sol_cut(:,:,isol)
         ENDDO
         CALL bin_close(gal_bin_unit)                
      ENDIF
c-----------------------------------------------------------------------
c     convert to perturbed b in radial direction.
c-----------------------------------------------------------------------
      IF (b_flag) THEN
         DO ip=0,tot_grids
            DO isol=1,gal%nsol
               singfac=mlow-nn*q(ip)+(/(ipert,ipert=0,mpert-1)/)
               sol(:,ip,isol)=twopi*ifac*psio*singfac*sol(:,ip,isol)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     output solution for ascii.
c-----------------------------------------------------------------------
      IF (out_galsol) THEN
         DO ising=1,msing
            WRITE (filename(2),"(I4)")ising
            filename(2)=ADJUSTL(filename(2))
            WRITE(filename(1),*) 'galsol_left_'
     $                          //TRIM(filename(2))//'.out'
            WRITE(filename(2),*) 'galsol_right_'
     $                          //TRIM(filename(2))//'.out'
            DO j=1,2
               isol=2*(ising-1)+j
               CALL ascii_open(gal_out_unit,TRIM(ADJUSTL(filename(j))),
     $              "REPLACE")
               WRITE (gal_out_unit,10) 'psifac'
               DO m=mlow,mhigh
                  WRITE (tmp,"(I4)") m
                  tmp=ADJUSTL(tmp)
                  WRITE (comp_tittle,*) 'REAL(',TRIM(tmp),')'
                  WRITE (gal_out_unit,10) TRIM(comp_tittle)
                  WRITE (comp_tittle,*) 'IMAG(',TRIM(tmp),')'
                  WRITE (gal_out_unit,10) TRIM(comp_tittle)
               ENDDO
!10             FORMAT (1P,A11,$)
10            FORMAT (1P,A15,$)
               WRITE (gal_out_unit,*)
               DO ip=0,tot_grids
                  IF (issing(ip)) CYCLE
                  WRITE (gal_out_unit,20) psi(ip)
!20                FORMAT (1P,E11.3,$)
!20               FORMAT (1P,E15.5,$)
20               FORMAT (1P,E20.10,$)
                  DO ipert=1,mpert
                     WRITE (gal_out_unit,20)
     $                     REAL(sol(ipert,ip,isol)),
     $                     IMAG(sol(ipert,ip,isol))
                  ENDDO
                  WRITE (gal_out_unit,*)
               ENDDO
               CALL ascii_close(gal_out_unit)
            ENDDO
         ENDDO
       ENDIF
c-----------------------------------------------------------------------
c     xdraw output.
c-----------------------------------------------------------------------
      IF (bin_galsol) THEN
         DO ising=1,msing
            WRITE (filename(2),"(I4)")ising
            filename(2)=ADJUSTL(filename(2))
            WRITE(filename(1),*) 'galsol_left_'
     $                           //TRIM(filename(2))//'.bin'
            WRITE(filename(2),*) 'galsol_right_'
     $                           //TRIM(filename(2))//'.bin'
            DO j=1,2
               isol=2*(ising-1)+j
               CALL bin_open(gal_bin_unit,TRIM(filename(j)),"UNKNOWN",
     $                        "REWIND","none")
               DO ipert=1,mpert
                  DO ip=0,tot_grids                
                     IF (issing(ip)) THEN
                        WRITE(gal_bin_unit)
                        CYCLE
                     ENDIF
                     realsinh=ASINH(REAL(sol(ipert,ip,isol))/eps)
                     imagsinh=ASINH(IMAG(sol(ipert,ip,isol))/eps)
                     WRITE (gal_bin_unit) REAL(psi(ip),4),
     $                     REAL(sol(ipert,ip,isol),4),
     $                     REAL(IMAG(sol(ipert,ip,isol)),4),
     $                     realsinh,imagsinh,
     $                     sing_log(sol(ipert,ip,isol))
                  ENDDO
                  WRITE(gal_bin_unit)
               ENDDO
               CALL bin_close(gal_bin_unit)
            ENDDO
         ENDDO
      ENDIF
      IF(bin_coilsol)THEN
         IF (coil%rpec_flag) THEN
            DO isol=coil%m1,coil%m2
               WRITE (filename(2),"(I4)")isol-2*msing
               filename(2)=ADJUSTL(filename(2))
               WRITE(filename(1),*) 'coilsol_'
     $                              //TRIM(filename(2))//'.bin'
               CALL bin_open(gal_bin_unit,TRIM(filename(1)),"UNKNOWN",
     $                       "REWIND","none")
               DO ipert=1,mpert
                  DO ip=0,tot_grids                
                     IF (issing(ip)) THEN
                        WRITE(gal_bin_unit)
                        CYCLE
                     ENDIF
                     realsinh=ASINH(REAL(sol(ipert,ip,isol))/eps)
                     imagsinh=ASINH(IMAG(sol(ipert,ip,isol))/eps)
                     WRITE (gal_bin_unit) REAL(psi(ip),4),
     $                      REAL(sol(ipert,ip,isol),4),
     $                      REAL(IMAG(sol(ipert,ip,isol)),4),
     $                      realsinh,imagsinh,
     $                      sing_log(sol(ipert,ip,isol))
                  ENDDO
                  WRITE(gal_bin_unit)
               ENDDO
               CALL bin_close(gal_bin_unit)
            ENDDO
         ENDIF
      ENDIF
      DEALLOCATE (issing,psi,q,sol,sol_cut)
      DEALLOCATE (xext)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gal_output_solution 
      END MODULE gal_mod
