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

      USE dcon_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. symmetrize.
c     Generate and symmetrize Delta matrices.
c-----------------------------------------------------------------------
      SUBROUTINE symmetrize(dp)

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: dp

      INTEGER :: m,ising,jsing,itheta,mthsurf,npsi,ipsi      
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: resnum,resm_sing, I
      INTEGER, DIMENSION(msing) :: resm      
      LOGICAL :: is_hermitian
      REAL(r8) :: tolerance
      REAL(r8), DIMENSION(0:mtheta) :: r,z,theta

      REAL(r8),DIMENSION(msing) :: L,f_L,f_S,nu_L,nu_S,DI,J,rho,shear

      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime,B_prime,
     $ Gamma_prime,Delta_prime
      COMPLEX(r8), DIMENSION(msing) :: A_prime_tmp,B_prime_tmp,
     $ Gamma_prime_tmp,Delta_prime_tmp
      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime_sym,B_prime_sym,
     $ Gamma_prime_sym,Delta_prime_sym

      REAL(r8) :: delpsi,jac,psifac,q,q1,rfac
      REAL(r8) :: respsi,psi_a
      REAL(r8), DIMENSION(2,2) :: w
      TYPE(spline_type) :: spl

      npsi = sq%mx+1

      psi_a=0.1 !!!!!!!!!!!!!!

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

      ! STILL NEED psi_a
      I = 0
      DO ipsi=1,npsi
         !CALL spline_eval(sq,psi(ipsi),0)
         I = I + 2.0*(sq%f(4) / (sq%fs(ipsi,1)/twopi)) !DB  this is F!
      ENDDO

      DO ising=1,msing
         !resnum=NINT(sing(ising)%q*nn)-mlow+1
         respsi=sing(ising)%psifac
         CALL spline_eval(sq,respsi,1)

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi=SQRT(w(1,1)**2+w(1,2)**2)

            J(ising)=J(ising)+(1.0/delpsi**2)*jac!*(1.0/twopi) !DB !!!
         ENDDO

         rho(ising) = (J(ising) * (sq%fs(respsi,1)/twopi))/sq%f(4) !DB
         shear(ising)=sing(ising)%q1

         L(ising)=I*(J(ising)*((resm(ising)*(sq%fs(respsi,1)/twopi) 
     $    )/sq%f(4))**2+(nn * psi_a)**2 ) !DB

         ! CALL mercier_mod
         DI(ising) = locstab%fs(respsi,1)/sq%xs(respsi) !DB

         nu_L(ising) = 0.5 - SQRT(-DI(ising)) !DB
         nu_S(ising) = 0.5 + SQRT(-DI(ising)) !DB

         f_L(ising) = (rho(ising) ** nu_L(ising)) * (((nu_S(ising) - 
     $       nu_L(ising)) /L(ising))**0.5)*shear(ising)*resm(ising) !DB
         f_S(ising) = (rho(ising) ** nu_S(ising)) * (((nu_S(ising) - 
     $       nu_L(ising)) /L(ising))**0.5)*shear(ising)*resm(ising) !DB

         A_prime_tmp = MATMUL(A_prime,f_L) !DB
         B_prime_tmp = MATMUL(B_prime,f_L) !DB
         Gamma_prime_tmp = MATMUL(Gamma_prime,f_L) !DB
         Delta_prime_tmp = MATMUL(Delta_prime,f_L) !DB

         A_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([A_prime_tmp], [1,msing])) !DB
         B_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([B_prime_tmp], [1,msing])) !DB
         Gamma_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([Gamma_prime_tmp], [1,msing])) !DB
         Delta_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([Delta_prime_tmp], [1,msing])) !DB
      ENDDO

      tolerance = 1.0e-2

      is_hermitian = all(abs(A_prime_sym - 
     $      transpose(conjg(A_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "A_prime is Hermitian"
      ELSE
         WRITE(*,*), "A_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(B_prime_sym - 
     $      transpose(conjg(B_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "B_prime is Hermitian"
      ELSE
         WRITE(*,*), "B_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(Gamma_prime_sym - 
     $      transpose(conjg(Gamma_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "Gamma_prime_sym is Hermitian"
      ELSE
         WRITE(*,*), "Gamma_prime_sym is not Hermitian"
      END IF
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