c-----------------------------------------------------------------------
c     file riccati.f.
c     implement Riccati integration of the Euler Lagrange equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     -. riccati_mod.
c     0. riccati_controller.
c     1. riccati_ode.
c     2. riccati_der.
c     3. riccati_inv_der.
c     4. riccati_inv_hamil_complex.
c     5. riccati_form_dx.
c     6. riccati_nojac.
c     7. riccati_inv.
c     8. riccati_eig.
c     9. riccati_mylog10.
c     10. riccati_get_p_from_u.
c-----------------------------------------------------------------------
c     subprogram 0. riccati_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE riccati_mod
      USE sing_mod
      USE zvode1_mod
      IMPLICIT NONE

      INTEGER :: peig_out_unit=53
      INTEGER :: peig_bin_unit=54
      INTEGER :: peig_size_unit=55

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 0. riccati_controller.
c     control the integration of the riccati ode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_controller(t0,te,X0,ric_dt,ric_tol,
     $     init_mode,first_call,last_call)

      REAL(r8), INTENT (INOUT) :: t0
      REAL(r8), INTENT(IN) :: te,ric_dt,ric_tol
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(INOUT) :: X0
      CHARACTER(7), INTENT(IN) :: init_mode
      LOGICAL, INTENT(IN) :: first_call, last_call

      CHARACTER(80) :: format1
      CHARACTER(7) :: ode_ps_mode
      LOGICAL :: bounce_it
      INTEGER :: ric_dir,i
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/1x,"istep",5x,"psi",3x,',i2.2,'(5x,"pe",i2.2,2x)/)')
c-----------------------------------------------------------------------
c     prep verbose output.
c-----------------------------------------------------------------------
      IF ( verbose_riccati_output ) THEN
         WRITE(format1,10)mpert
         IF (first_call) THEN
            OPEN(UNIT=peig_out_unit,FILE="peig.out",STATUS="REPLACE")
            OPEN(UNIT=peig_size_unit,FILE="Psize.out",STATUS="REPLACE")
            WRITE(peig_out_unit,format1)(i,i=1,mpert)
            OPEN(UNIT=peig_bin_unit,FILE="peig.bin",STATUS="REPLACE",
     $           FORM="UNFORMATTED")
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     initialize integration.
c-----------------------------------------------------------------------
      IF ( t0 > te ) THEN
         ric_dir = -1.0_r8
      ELSE
         ric_dir = 1.0_r8
      ENDIF
c-----------------------------------------------------------------------
c     call riccati integrator.
c-----------------------------------------------------------------------
      ode_ps_mode = init_mode
      DO
         IF (t0*ric_dir .GE. te*ric_dir) EXIT
         CALL riccati_ode(t0,te,X0,ric_dt,ric_tol,ode_ps_mode,bounce_it)
         IF (bounce_it) THEN
            SELECT CASE(ode_ps_mode)
            CASE("int_ric")
               ode_ps_mode = "inv_ric"
            CASE("inv_ric")
               ode_ps_mode = "int_ric"
            CASE DEFAULT
               CALL program_stop("Unknown riccati integration mode.")
            END SELECT
            CALL riccati_inv(X0)
         ENDIF
      ENDDO
      IF ( ode_ps_mode .NE. init_mode ) THEN
         CALL riccati_inv(X0)
      ENDIF
c-----------------------------------------------------------------------
c     symmetrize result for accuracy.
c-----------------------------------------------------------------------
      X0 = (X0 + CONJG(TRANSPOSE(X0))) / 2.0_r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      IF (verbose_riccati_output .AND. last_call) THEN
         WRITE(peig_out_unit,format1)(i,i=1,mpert)
         CLOSE(UNIT=peig_out_unit)
         CLOSE(UNIT=peig_bin_unit)
         CLOSE(UNIT=peig_size_unit)
      ENDIF

      RETURN
      END SUBROUTINE riccati_controller
c-----------------------------------------------------------------------
c     subprogram 1. riccati_ode.
c     integrate the Riccati equations with ZVODE.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_ode(t0,te,X0,ric_dt,ric_tol,ps_mode,
     $     bounce_flag)

      REAL(r8), INTENT (INOUT) :: t0
      REAL(r8), INTENT(IN) :: te,ric_dt,ric_tol
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(INOUT) :: X0
      CHARACTER(7), INTENT(IN) :: ps_mode !inv_ric, int_ric
      LOGICAL, INTENT(OUT) :: bounce_flag

      INTEGER :: neq,liw,lrw,lzw,info
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, ipiv
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: zwork, work
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: atol
      INTEGER :: iopt,itask,itol,mf,istate,ipert
      REAL(r8) :: rtol(1), atol0, maxatol
      COMPLEX(r8), DIMENSION(1) :: rpar
      INTEGER, DIMENSION(1) :: ipar

      REAL(r8) :: xcrit = 1.0E6
      REAL(r8) :: ric_dir
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      bounce_flag = .FALSE.
      neq = mpert*mpert
      liw = 30
      lrw = 20+neq
      lzw = 15*neq
      ALLOCATE(rwork(lrw), zwork(lzw), iwork(liw), atol(mpert,mpert),
     $     ipiv(mpert), work(mpert))
c-----------------------------------------------------------------------
c     prep the integration constants.
c-----------------------------------------------------------------------
      iopt=1
      itask=5
      itol=2
      mf=10
      istate=1
      iwork=0
      rwork=0.0_r8
      rwork(1)=te
      IF ( t0 > te ) THEN
         ric_dir = -1.0_r8
      ELSE
         ric_dir = 1.0_r8
      ENDIF
      rwork(5)=ric_dt * ric_dir
      rwork(11)=rwork(5)
      rtol(1) = ric_tol
c-----------------------------------------------------------------------
c     main riccati integration loop.
c-----------------------------------------------------------------------
      DO
c-----------------------------------------------------------------------
c     check integration exit condition.
c-----------------------------------------------------------------------
         IF (t0*ric_dir .GE. te*ric_dir) EXIT
c-----------------------------------------------------------------------
c     take integration step.
c-----------------------------------------------------------------------
         DO ipert=1,mpert
            atol0=MAXVAL(ABS(X0(:,ipert)))*ric_tol
            atol(:,ipert)=atol0
         ENDDO
         maxatol=MAX(MAXVAL(ABS(atol))*100,ric_tol)
         WHERE(atol == 0)
            atol=maxatol
         ENDWHERE

         SELECT CASE(ps_mode)
         CASE("int_ric")
            CALL ZVODE1(riccati_der,neq,X0,t0,te,
     $           itol,rtol,atol,itask,istate,iopt,zwork,lzw,
     $           rwork,lrw,iwork,liw,riccati_nojac,mf,rpar,ipar)
         CASE("inv_ric")
            CALL ZVODE1(riccati_inv_der,neq,X0,t0,te,
     $           itol,rtol,atol,itask,istate,iopt,zwork,lzw,
     $           rwork,lrw,iwork,liw,riccati_nojac,mf,rpar,ipar)
         CASE DEFAULT
            CALL program_stop("Unknown riccati integration mode.")
         END SELECT
c-----------------------------------------------------------------------
c     toggle integration mode to control growth of matrix.
c-----------------------------------------------------------------------
         IF (riccati_bounce) THEN
            IF (SUM(ABS(X0)) > xcrit) THEN
               bounce_flag = .TRUE.
               RETURN
            ENDIF
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_ode
c-----------------------------------------------------------------------
c     subprogram 2. riccati_der.
c     gets the derivative for the riccati equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_der(neq,startPsi,X,dX,rpar,ipar)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: startPsi
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: X
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: dX
      COMPLEX(r8), INTENT(IN) :: rpar, ipar

      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: ac

      COMPLEX(r8), DIMENSION(mpert,mpert) :: pmat,pmat1
      INTEGER :: istep
c-----------------------------------------------------------------------
c     calculate derivative.
c-----------------------------------------------------------------------
      CALL sing_derL(startPsi,ac)
      CALL riccati_form_dx(mpert,ac,X,dX)
c-----------------------------------------------------------------------
c     write verbose riccati output.
c-----------------------------------------------------------------------
         IF ( verbose_riccati_output ) THEN
            IF (SUM(ABS(X)) > 0) THEN
               pmat = X
               pmat1 = dX

               istep = 0
               CALL riccati_eig(istep,startPsi,pmat,peig_out_unit,
     $              peig_bin_unit,pmat)
               CALL riccati_modesize(startPsi,pmat)
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_der
c-----------------------------------------------------------------------
c     subprogram 3. riccati_inv_der.
c     gets the derivative for the inverse riccati equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_inv_der(neq,startPsi,X,dX,rpar,ipar)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: startPsi
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: X
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: dX
      COMPLEX(r8), INTENT(IN) :: rpar, ipar

      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: ac

      COMPLEX(r8), DIMENSION(mpert,mpert) :: pmat,pmat1,pmat1t
      INTEGER :: istep
c-----------------------------------------------------------------------
c     calculate derivative.
c-----------------------------------------------------------------------
      CALL sing_derL(startPsi,ac)
      CALL riccati_inv_hamil_complex(mpert,ac)
      CALL riccati_form_dx(mpert,ac,X,dX)
c-----------------------------------------------------------------------
c     write verbose riccati output.
c-----------------------------------------------------------------------
         IF ( verbose_riccati_output ) THEN
            IF (SUM(ABS(X)) > 0) THEN
               pmat = X
               pmat1 = dX
               CALL riccati_inv(pmat)
               ! (A^(-1))' = -A^(-1)*A'*A^(-1)
               pmat1t = 0
               CALL ZGEMM('N','N',mpert,mpert,mpert,one_c,pmat,mpert,
     $              dX,mpert,zero_c,pmat1t,mpert)
               CALL ZGEMM('N','N',mpert,mpert,mpert,one_c,pmat1t,mpert,
     $              pmat,mpert,zero_c,pmat1,mpert)
               pmat1 = -pmat1

               istep = 0
               CALL riccati_eig(istep,startPsi,pmat,peig_out_unit,
     $              peig_bin_unit,pmat)
               CALL riccati_modesize(startPsi,pmat)
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_inv_der
c-----------------------------------------------------------------------
c     subprogram 4. riccati_inv_hamil_complex.
c     converts complex hamiltonian to its "riccati inverse."
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_inv_hamil_complex(m1,a)

      INTEGER, INTENT(IN) :: m1
      COMPLEX(r8), DIMENSION(1:2*m1,1:2*m1), INTENT(INOUT) :: a

      COMPLEX(r8), DIMENSION(1:m1,1:m1) :: a_ij
      INTEGER :: m2
c-----------------------------------------------------------------------
c     initialize constants.
c-----------------------------------------------------------------------
      m2 = 2*m1
c-----------------------------------------------------------------------
c     compute the hamiltonian for the inverse riccati equation.
c-----------------------------------------------------------------------
      !Swap A11 and A22
      a_ij = a(1:m1,1:m1)
      a(1:m1,1:m1) = a(m1+1:m2,m1+1:m2)
      a(m1+1:m2,m1+1:m2) = a_ij
      !Swap A12 and A21
      a_ij = a(1:m1,m1+1:m2)
      a(1:m1,m1+1:m2) = a(m1+1:m2,1:m1)
      a(m1+1:m2,1:m1) = a_ij
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_inv_hamil_complex
c-----------------------------------------------------------------------
c     subprogram 5. riccati_form_dx.
c     form the RHS of the Riccati equation from the hamiltonian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_form_dx(m1,ac,X,dX)

      INTEGER, INTENT(IN) :: m1
      COMPLEX(r8), DIMENSION(2*m1,2*m1), INTENT(IN) :: ac
      COMPLEX(r8), DIMENSION(m1,m1), INTENT(IN) :: X
      COMPLEX(r8), DIMENSION(m1,m1), INTENT(OUT) :: dX

      COMPLEX(r8), DIMENSION(m1,m1) :: Xtemp, acTemp
      INTEGER :: m2
c-----------------------------------------------------------------------
c     initialize constants.
c-----------------------------------------------------------------------
      m2 = 2*m1
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
      !dX = A21
      dX = ac(m1+1:m2,1:m1)
      !dX = A21 + A22*X
      acTemp = ac(m1+1:m2,m1+1:m2)
      CALL zgemm('N','N',m1,m1,m1,one_c,acTemp,m1,X,m1,one_c,dX,m1)
      !dX = A21 + A22*X - X*A11
      acTemp = ac(1:m1,1:m1)
      CALL zgemm('N','N',m1,m1,m1,negone_c,X,m1,acTemp,m1,one_c,dX,m1)
      !dX = A21 + A22*X - X*A11 - X*A12*X
      Xtemp = 0.0_r8
      acTemp = ac(1:m1,m1+1:m2)
      CALL zgemm('N','N',m1,m1,m1,one_c,acTemp,m1,X,m1,one_c,Xtemp,m1)
      CALL zgemm('N','N',m1,m1,m1,negone_c,X,m1,Xtemp,m1,one_c,dX,m1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_form_dx
c-----------------------------------------------------------------------
c     subprogram 6. riccati_nojac.
c     a dummy Jacobian call for the LSODE/ZVODE package
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_nojac(neq_jac, t, y, ml, mu, pd, nrpd)
         INTEGER  neq_jac, ml, mu, nrpd
         REAL(r8)  t, y, pd(nrpd,2)
         pd(:,:) = 0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_nojac
c-----------------------------------------------------------------------
c     subprogram 7. riccati_inv.
c     invert riccati matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_inv(X)
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(INOUT) :: X

      COMPLEX(r8), DIMENSION(mpert,mpert) :: ident,X_init,Y
      COMPLEX(r8), DIMENSION(mpert*mpert) :: zwork
      INTEGER, DIMENSION(mpert) :: ipiv
      INTEGER :: lzw,info,i
c-----------------------------------------------------------------------
c     preliminaries.
c-----------------------------------------------------------------------
      ident = 0.0_r8
      DO i = 1,mpert
         ident(i,i) = 1.0_r8
      ENDDO
      lzw = mpert*mpert
      zwork = 0.0_r8
      X_init = X
      Y = ident
c-----------------------------------------------------------------------
c     invert matrix.
c-----------------------------------------------------------------------
      CALL ZHETRF('U',mpert,X,mpert,ipiv,zwork,lzw,info)
      CALL ZHETRS('U',mpert,mpert,X,mpert,ipiv,Y,mpert,info)
      X = Y
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_inv
c-----------------------------------------------------------------------
c     subprogram 8. riccati_eig.
c     get and print eigenvalues of riccati matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_eig(nstep,t,X,out_unit,bin_unit,dX)
      INTEGER, INTENT(IN) :: nstep, out_unit, bin_unit
      REAL(r8), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: X
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: dX

      CHARACTER(80) :: format2
      REAL(r8), DIMENSION(mpert) :: eigval
      REAL(r8), DIMENSION(3*mpert-1) :: rwork_diagnose
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work_diagnose
      INTEGER, DIMENSION(1) :: jmin
      INTEGER :: lwork,info,i
      REAL(r8) :: crit1, crit2
      REAL(r8) :: eloglim = 1e-5
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 20   FORMAT('(i6,1p,',i8.8,'e15.7)')
      WRITE(format2,20)mpert+1
c-----------------------------------------------------------------------
c     preliminaries.
c-----------------------------------------------------------------------
      lwork=2*mpert-1
      CALL zheev('N','U',mpert,X,mpert,eigval,work_diagnose,
     $     lwork,rwork_diagnose,info)
      jmin=MINLOC(ABS(eigval))
      crit1=eigval(jmin(1)) * sq%f(3)**2
      crit2=LOG10(ABS(PRODUCT(eigval)))
      WRITE(out_unit,format2)nstep,t,eigval
      DO i=1,mpert
         WRITE(bin_unit)REAL(t,4),REAL(eigval(i),4),
     $        riccati_mylog10(eigval(i),eloglim)
      ENDDO
      WRITE(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_eig
c-----------------------------------------------------------------------
c     subprogram 9. riccati_modesize.
c     just temporary for presentation purposes
c-----------------------------------------------------------------------
      SUBROUTINE riccati_modesize(t,X)
      REAL(r8), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN) :: X
      REAL(r8), DIMENSION(mpert) :: col_mod
      CHARACTER(80) :: format3
 70   FORMAT('(',i8.8,'e15.7)')
      WRITE(format3,70)mpert+1
      col_mod = SUM(ABS(X),DIM=1)
      WRITE(peig_size_unit,format3)t,col_mod
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_modesize
c-----------------------------------------------------------------------
c     subprogram 9. riccati_mylog10.
c     computes and returns integer riccati_mylog10.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION riccati_mylog10(x,xmin) RESULT(y)

      REAL(r8), INTENT(IN) :: x,xmin
      REAL(r4) :: y
c-----------------------------------------------------------------------
c     compute riccati_mylog10.
c-----------------------------------------------------------------------
      IF(ABS(x) < xmin)THEN
         y=LOG10(ABS(xmin))
      ELSE
         y=LOG10(ABS(x))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION riccati_mylog10
c-----------------------------------------------------------------------
c     subprogram 10. riccati_get_p_from_u.
c     computes the riccati matrix from (half) the u matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE riccati_get_p_from_u(uhalf,P)

      COMPLEX(r8), DIMENSION(2*mpert,mpert), INTENT(IN) :: uhalf
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(OUT) :: P
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp
      INTEGER, DIMENSION(mpert) :: ipiv
      INTEGER :: info
c-----------------------------------------------------------------------
c     compute riccati matrix P.
c-----------------------------------------------------------------------
      temp = CONJG(TRANSPOSE(uhalf(1:mpert,:)))
      P = CONJG(TRANSPOSE(uhalf(mpert+1:2*mpert,:)))
      !P = temp^{-1}*P
      CALL zgetrf(mpert,mpert,temp,mpert,ipiv,info)
      CALL zgetrs('N',mpert,mpert,temp,mpert,ipiv,P,mpert,info)
      P = (P + CONJG(TRANSPOSE(P)))/2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE riccati_get_p_from_u

      END MODULE riccati_mod
