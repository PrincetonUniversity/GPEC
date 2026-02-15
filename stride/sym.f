c-----------------------------------------------------------------------
c     file stride_netcdf.f
c     writes stride.out information to a netcdf file
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. stride_netcdf_mod
c     1. check
c     2. stride_netcdf_out
c-----------------------------------------------------------------------
c     subprogram 0. stride_netcdf_mod
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sym_mod

      USE stride_dcon_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. symmetrize.
c     Generate and symmetrize Delta matrices.
c-----------------------------------------------------------------------
      SUBROUTINE symmetrize(dp)

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: dp

      INTEGER :: m,ising,jsing,itheta,mthsurf,ipsi,jpsi,npsi,idx,
     $           si,sj
      REAL(r8) :: bsq,chi1,dpsisq,eta,jac,p1,psifac,q,q1,r,
     $     rfac,theta,twopif,v1,v2,v21,v22,v23,v33
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: resnum,resm_sing,min_diff,dtheta
      INTEGER, DIMENSION(msing) :: resm      
      LOGICAL :: is_hermitian
      REAL(r8) :: tolerance,dx,angle
      !REAL(r8), DIMENSION(:), ALLOCATABLE :: r,z,theta
      REAL(r8), DIMENSION(sq%mx+1) :: I_psin,ln_q
      REAL(r8),DIMENSION(msing) :: L,f_L,f_S,nu_L,nu_S,DI,J,rho,shear

      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime,B_prime,
     $ Gamma_prime,Delta_prime
      COMPLEX(r8), DIMENSION(msing) :: A_prime_tmp,B_prime_tmp,
     $ Gamma_prime_tmp,Delta_prime_tmp
      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime_sym,B_prime_sym,
     $ Gamma_prime_sym,Delta_prime_sym,tmp_arr

      REAL(r8) :: delpsi
      REAL(r8) :: respsi,psi_a
      REAL(r8), DIMENSION(2,2) :: w
      TYPE(spline_type) :: spl,psi_t,I_spl,shr_spl,J_spl

      REAL(r8), DIMENSION(:), POINTER :: avg
      TYPE(spline_type), TARGET :: fspl

      !WRITE(*,*)"mthsurf=",mthsurf

      !mthsurf=40*mthsurf0*MAX(ABS(mlow),ABS(mhigh))
      !ALLOCATE(r(0:mthsurf),z(0:mthsurf),theta(0:mthsurf))

      npsi = sq%mx

      !chi1=twopi*psio
      psi_a=twopi*psio
      WRITE(*,*)"psi_a=",psi_a

      WRITE(*,*)"sq%mx=",sq%mx
      WRITE(*,*)"npsi=",npsi

      ! construct PEST3 matching data (keep synced with RDCON!)
      A_prime=0.0
      B_prime=0.0
      Gamma_prime=0.0
      Delta_prime=0.0
      DO ising=1,msing
         DO jsing=1,msing
            A_prime(ising,jsing)=dp(2*ising,2*jsing)
     $              +dp(2*ising,2*jsing-1)
     $              +dp(2*ising-1,2*jsing)
     $              +dp(2*ising-1,2*jsing-1)
            B_prime(ising,jsing)=dp(2*ising,2*jsing)
     $              -dp(2*ising,2*jsing-1)
     $              +dp(2*ising-1,2*jsing)
     $              -dp(2*ising-1,2*jsing-1)
            Gamma_prime(ising,jsing)=dp(2*ising,2*jsing)
     $              +dp(2*ising,2*jsing-1)
     $              -dp(2*ising-1,2*jsing)
     $              -dp(2*ising-1,2*jsing-1)
            Delta_prime(ising,jsing)=dp(2*ising,2*jsing)
     $              -dp(2*ising,2*jsing-1)
     $              -dp(2*ising-1,2*jsing)
     $              +dp(2*ising-1,2*jsing-1)
         ENDDO
      ENDDO

      mvec=(/(m,m=mlow,mhigh)/)
      DO ising=1,msing
         respsi=sing(ising)%psifac
         resnum=NINT(sing(ising)%q*nn)-mlow+1
         resm_sing=mvec(resnum)
         resm(ising)=resm_sing
      ENDDO

      WRITE(*,*)"resnum=",resnum
      WRITE(*,*)"resm=",resm

      ! Prepare toroidal flux spline
      CALL spline_alloc(psi_t,SIZE(sq%fsi(:, 4)),1)
      psi_t%xs=sq%xs(:)
      psi_t%fs(:,1)=sq%fsi(:,4)*twopi*psio ! Un-normalize toroidal flux
      CALL spline_fit(psi_t,"extrap")

      ! Prepare shear spline
      CALL spline_alloc(shr_spl,SIZE(sq%xs(:)),1)
      shr_spl%xs=sq%xs(:)
      ln_q=LOG(sq%fs(:,4)) 
      shr_spl%fs(:,1)=ln_q! log(q)
      CALL spline_fit(shr_spl,"extrap")

      ! Calculate I(psi_N) integral
      I_psin(1) = 0.0  ! integral at x=0 is 0
      DO ipsi = 2, npsi
         I_psin(ipsi) = 0.0
         DO jpsi = 1, ipsi-1
            dx = sq%xs(jpsi+1) - sq%xs(jpsi)
            CALL spline_eval(psi_t,sq%xs(jpsi),1)
            I_psin(ipsi) = I_psin(ipsi) + 
     $       (2.0*sq%fs(jpsi,4)/(psi_t%f(1)))*dx
         ENDDO
      ENDDO

      ! Prepare I(psiN) spline
      CALL spline_alloc(I_spl,SIZE(sq%xs(:)),1)
      I_spl%xs=sq%xs(:)
      I_spl%fs(:,1)=I_psin
      CALL spline_fit(I_spl,"extrap")

      ! Compute J integral
      CALL spline_alloc(fspl,mtheta,2)
      fspl%xs=rzphi%ys
      CALL spline_alloc(J_spl,mpsi,1)
      J_spl%xs=sq%xs(:)
      WRITE(*,*)'(mtheta)=',mtheta

      DO ipsi=0,mpsi
         psifac=sq%xs(ipsi)
         twopif=sq%fs(ipsi,1)
         p1=sq%fs1(ipsi,2)
         v1=sq%fs(ipsi,3)
         v2=sq%fs1(ipsi,3)
         q=sq%fs(ipsi,4)
         q1=sq%fs1(ipsi,4)
         chi1=twopi*psio
c-----------------------------------------------------------------------
c     evaluate coordinates and jacobian.
c-----------------------------------------------------------------------
         DO itheta=0,mtheta
            CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            theta=rzphi%ys(itheta)
            rfac=SQRT(rzphi%f(1))
            eta=twopi*(theta+rzphi%f(2))
            r=ro+rfac*COS(eta)
            jac=rzphi%f(4)
c-----------------------------------------------------------------------
c     evaluate other local quantities.
c-----------------------------------------------------------------------
            v21=rzphi%fy(1)/(2*rfac*jac)
            v22=(1+rzphi%fy(2))*twopi*rfac/jac
            v23=rzphi%fy(3)*r/jac
            v33=twopi*r/jac
            bsq=chi1**2*(v21**2+v22**2+(v23+q*v33)**2)
            dpsisq=(twopi*r)**2*(v21**2+v22**2)
c-----------------------------------------------------------------------
c     evaluate integrands.
c-----------------------------------------------------------------------
            fspl%fs(itheta,1)=1/dpsisq!*(chi1**2.0)) dpsisq is grad(psi_N)**2 ??
            fspl%fs(itheta,2)=bsq
            fspl%fs(itheta,:)=fspl%fs(itheta,:)*(jac/v1)
         ENDDO
c-----------------------------------------------------------------------
c     integrate quantities with respect to theta.
c-----------------------------------------------------------------------
         CALL spline_fit(fspl,"periodic")
         CALL spline_int(fspl)
         avg => fspl%fsi(mtheta,:)
         J_spl%fs(ipsi,1)=avg(1)

      ENDDO
      CALL spline_dealloc(fspl)
      CALL spline_fit(J_spl,"extrap")

      ! Loop across rational surfaces to evaluate remaining quantities
      DO ising=1,msing
         respsi=sing(ising)%psifac
         WRITE(*,*)"respsi=",respsi

         ! Evaluate splines on rational surface
         CALL spline_eval(sq,respsi,1)
         CALL spline_eval(locstab,respsi,1)
         CALL spline_eval(psi_t,respsi,1)
         CALL spline_eval(I_spl,respsi,1)
         CALL spline_eval(shr_spl,respsi,1)
         CALL spline_eval(J_spl,respsi,1)

         J(ising)=J_spl%f(1)
         rho(ising) = (J(ising) * (psi_t%f(1)))/sq%f(4) !DB

         WRITE(*,*)"GPEC shear=",sing(ising)%q1
         WRITE(*,*)"d(ln(q))/dpsi(ising)=",shr_spl%f1(1)

         ! Combine integrated quantities into L(psi_N)
         L(ising)=I_spl%f(1)*(J(ising)*(((resm(ising)*
     $    psi_t%f(1))/sq%f(4))**2.0) + (nn * psi_a)**2.0 ) !DB

         DI(ising) = locstab%f(1)/respsi !DB

         nu_L(ising) = 0.5 - SQRT(-DI(ising)) !DB
         nu_S(ising) = 0.5 + SQRT(-DI(ising)) !DB

         f_L(ising) = (rho(ising) ** nu_L(ising)) * (((nu_S(ising)- 
     $       nu_L(ising)) /L(ising))**0.5)*shr_spl%f1(1)*resm(ising) !DB
         f_S(ising) = (rho(ising) ** nu_S(ising)) * (((nu_S(ising)- 
     $       nu_L(ising)) /L(ising))**0.5)*shr_spl%f1(1)*resm(ising) !DB

         ! First component of vector*matrix*vector multiply
      !   A_prime_tmp = MATMUL(A_prime,f_L) !DB
      !   B_prime_tmp = MATMUL(B_prime,f_L) !DB
      !   Gamma_prime_tmp = MATMUL(Gamma_prime,f_L) !DB

c         Delta_prime_tmp = MATMUL(Delta_prime,f_L) !DB

         ! Second component of vector*matrix*vector multiply
      !   A_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
      !$    RESHAPE([A_prime_tmp], [1,msing])) !DB
      !   B_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
      !$    RESHAPE([B_prime_tmp], [1,msing])) !DB
      !   Gamma_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
      !$!    RESHAPE([Gamma_prime_tmp], [1,msing])) !DB

c         Delta_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
c     $    RESHAPE([Delta_prime_tmp], [1,msing])) !DB
      ENDDO

      ! New cosine method from Richard's paper
      DO si=1,msing
         DO sj=1,msing
            Delta_prime_sym(si,sj) = COS((si + sj)*pi)*
     $                     (f_S(si)/f_L(sj))*Delta_prime(si,sj)
         ENDDO
      ENDDO

      tmp_arr=(Delta_prime - transpose(conjg(Delta_prime)))
      WRITE(*,*)"Original delta diff="
      DO si = 1, msing
         DO sj = 1, msing
         WRITE(*,'(F8.3)',advance='no') REAL(tmp_arr(si,sj))
         IF (sj < msing) WRITE(*,'(A)', advance='no') '  '
         ENDDO
         WRITE(*,*)  ! New line after each row
      ENDDO

      WRITE(*,*)"Delta_prime="
      DO si = 1, msing
         DO sj = 1, msing
            WRITE(*,'(F8.3)', advance='no') REAL(Delta_prime(si,sj))
            IF (sj < msing) WRITE(*,'(A)', advance='no') '  '
         ENDDO
         WRITE(*,*)  ! New line after each row
      ENDDO
      WRITE(*,*)"Delta_prime_sym="
      DO si = 1, msing
         DO sj = 1, msing
         WRITE(*,'(F8.3)',advance='no') REAL(Delta_prime_sym(si,sj))
         IF (sj < msing) WRITE(*,'(A)', advance='no') '  '
         ENDDO
         WRITE(*,*)  ! New line after each row
      ENDDO

c      WRITE(*, '(3(F5.3, 1X))') ((abs(Delta_prime_sym - 
c     $      transpose(conjg(Delta_prime_sym)))(i, j), j = 1, 
c     $                            msing), i = 1, msing)

      !WRITE(*,*)"B_prime=",B_prime
      !WRITE(*,*)"Gamma_prime_sym=",Gamma_prime_sym
      !WRITE(*,*)"B_Gamma_diff=",(B_prime_sym - 
      !$      transpose(conjg(Gamma_prime_sym)))
      tolerance = 1.0!e-01

c      is_hermitian = all(abs(A_prime_sym - 
c     $      transpose(conjg(A_prime_sym))) < tolerance)
c      IF (is_hermitian) then
c         WRITE(*,*), "A_prime is Hermitian"
c      ELSE
c         WRITE(*,*), "A_prime is not Hermitian"
c      END IF
c      is_hermitian = all(abs(Gamma_prime_sym - 
c     $      transpose(conjg(B_prime_sym))) < tolerance)
c      IF (is_hermitian) then
c         WRITE(*,*), "B_prime is Hermitian"
c      ELSE
c         WRITE(*,*), "B_prime is not Hermitian"
c      END IF
c      is_hermitian = all(abs(B_prime_sym - 
c     $      transpose(conjg(Gamma_prime_sym))) < tolerance)
c      IF (is_hermitian) then
c         WRITE(*,*), "Gamma_prime_sym is Hermitian"
c      ELSE
c         WRITE(*,*), "Gamma_prime_sym is not Hermitian"
c      END IF

      tmp_arr=(Delta_prime_sym - transpose(conjg(Delta_prime_sym)))
      WRITE(*,*)"New delta diff="
      DO si = 1, msing
         DO sj = 1, msing
         WRITE(*,'(F8.3)',advance='no') REAL(tmp_arr(si,sj))
         IF (sj < msing) WRITE(*,'(A)', advance='no') '  '
         ENDDO
         WRITE(*,*)  ! New line after each row
      ENDDO

      is_hermitian = all(abs(Delta_prime_sym - 
     $      transpose(conjg(Delta_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "Delta_prime_sym is Hermitian"
      ELSE
         WRITE(*,*), "Delta_prime_sym is not Hermitian"
      END IF

      RETURN
      END SUBROUTINE symmetrize
      END MODULE sym_mod