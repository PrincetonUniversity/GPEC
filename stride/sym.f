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
      SUBROUTINE symmetrize(msing,dp)

      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: dp

      INTEGER, INTENT(IN) :: msing, cmpsi
      COMPLEX(r8), INTENT(IN) :: dp
      sq
      

      INTEGER :: ising,jsing
      LOGICAL :: is_hermitian
      REAL(r8) :: tolerance
      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime,B_prime,
     $ Gamma_prime,Delta_prime
      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime_tmp,B_prime_tmp,
     $ Gamma_prime_tmp,Delta_prime_tmp
      COMPLEX(r8), DIMENSION(msing,msing) :: A_prime_sym,B_prime_sym,
     $ Gamma_prime_sym,Delta_prime_sym

      ! Construct PEST3 matrices
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

      ! STILL NEED G(ipsi), psi_a
      I = 0
      DO ipsi=1,cmpsi
         CALL spline_eval(sq,psi(ipsi),0)
         I = I + 2.0*(sq%f(4) / G(ipsi)) !DB !!!

      ENDDO

      DO ising=1,msing
         resnum=NINT(singtype(ising)%q*nn)-mlow+1
         respsi=singtype(ising)%psifac
         CALL spline_eval(sq,respsi,1)

         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,respsi,theta(itheta),1)
            jac=rzphi%f(4)
            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
            delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)

            J(ising)=J(ising)+(1.0/delpsi(itheta)**2)*jac!*(1.0/twopi) !DB !!!
         ENDDO

         rho(ising) = (J(ising) * G(respsi))/sq%f(4) !DB
         shear(ising)=mfac(resnum)*sq%f1(4)/sq%f(4)**2

         L(ising) = I * (J(ising) * (( resm * G(ising) )/sq%f(4))**2+
     $    (nn * psi_a)**2 ) !DB

         ! CALL mercier_mod
         DI(ising) = locstab%fs(respsi,1)/sq%xs(respsi) !DB

         nu_L(ising) = 0.5 - SQRT(-DI(ising)) !DB
         nu_S(ising) = 0.5 + SQRT(-DI(ising)) !DB

         f_L(ising) = (rho(ising) ** nu_L(ising)) * (((nu_S(ising) - 
     $       nu_L(ising)) / L(ising))**0.5)*shear(ising)*resm !DB
         f_S(ising) = (rho(ising) ** nu_S(ising)) * (((nu_S(ising) - 
     $       nu_L(ising)) / L(ising))**0.5)*shear(ising)*resm !DB

         A_prime_tmp = MATMUL(A_prime,f_L) !DB
         B_prime_tmp = MATMUL(B_prime,f_L) !DB
         Gamma_prime_tmp = MATMUL(Gamma_prime,f_L) !DB
         Delta_prime_tmp = MATMUL(Delta_prime,f_L) !DB

         A_prime_sym = MATMUL(RESHAPE([f_S], [1,n]),A_prime_tmp) !DB
         B_prime_sym = MATMUL(RESHAPE([f_S], [1,n]),B_prime_tmp) !DB
         Gamma_prime_sym = MATMUL(RESHAPE([f_S], 
     $       [1,n]),Gamma_prime_tmp) !DB
         Delta_prime_sym = MATMUL(RESHAPE([f_S], 
     $       [1,n]),Delta_prime_tmp) !DB
      ENDDO

      tolerance = 1.0e-2

      is_hermitian = all(abs(A_prime_sym - 
     $      transpose(conjg(A_prime_sym))) < tolerance)
      IF (is_hermitian) then
         PRINT(*,*), "A_prime is Hermitian"
      ELSE
         PRINT(*,*), "A_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(B_prime_sym - 
     $      transpose(conjg(B_prime_sym))) < tolerance)
      IF (is_hermitian) then
         PRINT(*,*), "B_prime is Hermitian"
      ELSE
         PRINT(*,*), "B_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(Gamma_prime_sym - 
     $      transpose(conjg(Gamma_prime_sym))) < tolerance)
      IF (is_hermitian) then
         PRINT(*,*), "Gamma_prime_sym is Hermitian"
      ELSE
         PRINT(*,*), "Gamma_prime_sym is not Hermitian"
      END IF
      is_hermitian = all(abs(Delta_prime_sym - 
     $      transpose(conjg(Delta_prime_sym))) < tolerance)
      IF (is_hermitian) then
         PRINT(*,*), "Delta_prime_sym is Hermitian"
      ELSE
         PRINT(*,*), "Delta_prime_sym is not Hermitian"
      END IF