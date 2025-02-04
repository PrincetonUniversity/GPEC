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

      INTEGER :: m,ising,jsing,itheta,mthsurf,ipsi,jpsi,npsi,idx
      INTEGER, DIMENSION(mpert) :: mvec
      REAL(r8) :: resnum,resm_sing,min_diff,dtheta
      INTEGER, DIMENSION(msing) :: resm      
      LOGICAL :: is_hermitian
      REAL(r8) :: tolerance,dx,angle,chi1
      REAL(r8), DIMENSION(:), ALLOCATABLE :: r,z,theta
      REAL(r8), DIMENSION(sq%mx+1) :: I_psin
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

      WRITE(*,*)"mthsurf=",mthsurf

      mthsurf=40*mthsurf0*MAX(ABS(mlow),ABS(mhigh))
      ALLOCATE(r(0:mthsurf),z(0:mthsurf),theta(0:mthsurf))

      npsi = sq%mx

      !chi1=twopi*psio
      psi_a=psio
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

      ! Calculate I(psi_N) integral
      I_psin(1) = 0.0  ! integral at x=0 is 0
      DO ipsi = 2, npsi
         I_psin(ipsi) = 0.0
         DO jpsi = 1, ipsi-1
            dx = sq%xs(jpsi+1) - sq%xs(jpsi)
            I_psin(ipsi) = I_psin(ipsi) + 
     $       (2.0*sq%fs(jpsi,4)/(sq%fs(jpsi,1)/twopi))*dx ! check q
         END DO
      END DO

      ! Prepare theta quantities for J(psi_N) integral
      CALL spline_alloc(spl,mthsurf,4)
      theta=(/(itheta,itheta=0,mthsurf)/)/REAL(mthsurf,r8)
      spl%xs=theta
      psifac=psilim

      ! Loop across rational surfaces to evaluate remaining quantities
      DO ising=1,msing
         !resnum=NINT(sing(ising)%q*nn)-mlow+1
         respsi=sing(ising)%psifac
         WRITE(*,*)"respsi=",respsi

         ! Find correct I(psiN) index
         min_diff = abs(sq%xs(1) - respsi)
         idx = 1
    
         do ipsi = 2, npsi
            if (abs(sq%xs(ipsi) - respsi) < min_diff) then
               min_diff = abs(sq%xs(ipsi) - respsi)
               idx = ipsi
            end if
         end do

         WRITE(*,*)"idx=",idx 

         ! Evaluate splines on rational surface
         CALL spline_eval(sq,respsi,1)
         CALL spline_eval(locstab,respsi,1)

         WRITE(*,*)"SIZE(theta)=",SIZE(theta)
         WRITE(*,*)"mthsurf=",mthsurf

         ! Calculate J(psi_N) contour integral
         dtheta = twopi / real(mthsurf)
         DO itheta=0,mthsurf
            CALL bicube_eval(rzphi,psilim,theta(itheta),1)
            rfac=SQRT(rzphi%f(1))
            angle=twopi*(theta(itheta)+rzphi%f(2))
            r(itheta)=ro+rfac*COS(angle)
            z(itheta)=zo+rfac*SIN(angle)
            jac=rzphi%f(4)

            w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
            w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)

            delpsi=SQRT(w(1,1)**2+w(1,2)**2)

            J(ising)=J(ising)+
     $       (1.0/(delpsi**2))*(dtheta/twopi)!/jac !DB !!
         ENDDO
         
         WRITE(*,*)"J=",J

         rho(ising) = (J(ising) * (sq%f(1)/twopi))/sq%f(4) !DB
         shear(ising) = sing(ising)%q1
         WRITE(*,*)"rho=",rho
         WRITE(*,*)"shear(ising)=",shear(ising)

         WRITE(*,*)"sq%f(4)=",sq%f(4)

         ! Combine integrated quantities into L(psi_N)
         L(ising)=I_psin(idx)*(J(ising)*((resm(ising)*
     $    (sq%f(1)/twopi))/sq%f(4))**2+(nn * psi_a)**2 ) !DB

         WRITE(*,*)"L=",L
         ! CALL mercier_mod
         DI(ising) = locstab%f(1)/respsi !DB
         WRITE(*,*)"DI=",DI

         nu_L(ising) = 0.5 - SQRT(-DI(ising)) !DB
         nu_S(ising) = 0.5 + SQRT(-DI(ising)) !DB

         f_L(ising) = (rho(ising) ** nu_L(ising)) * (((nu_S(ising)- 
     $       nu_L(ising)) /L(ising))**0.5)*shear(ising)*resm(ising) !DB
         f_S(ising) = (rho(ising) ** nu_S(ising)) * (((nu_S(ising)- 
     $       nu_L(ising)) /L(ising))**0.5)*shear(ising)*resm(ising) !DB

         ! First component of vector*matrix*vector multiply
         A_prime_tmp = MATMUL(A_prime,f_L) !DB
         B_prime_tmp = MATMUL(B_prime,f_L) !DB
         Gamma_prime_tmp = MATMUL(Gamma_prime,f_L) !DB
         Delta_prime_tmp = MATMUL(Delta_prime,f_L) !DB

         ! Second component of vector*matrix*vector multiply
         A_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([A_prime_tmp], [1,msing])) !DB
         B_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([B_prime_tmp], [1,msing])) !DB
         Gamma_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([Gamma_prime_tmp], [1,msing])) !DB
         Delta_prime_sym = MATMUL(RESHAPE([f_S], [1,msing]),
     $    RESHAPE([Delta_prime_tmp], [1,msing])) !DB
      ENDDO

      WRITE(*,*)"Delta_prime=",Delta_prime
      WRITE(*,*)"Delta_prime_sym=",Delta_prime_sym
      WRITE(*,*)"Delta_diff=",abs(Delta_prime_sym - 
     $      transpose(conjg(Delta_prime_sym)))
      !WRITE(*,*)"B_prime=",B_prime
      !WRITE(*,*)"Gamma_prime_sym=",Gamma_prime_sym
      WRITE(*,*)"B_Gamma_diff=",(B_prime_sym - 
     $      transpose(conjg(Gamma_prime_sym)))
      tolerance = 1.0e-01

      is_hermitian = all(abs(A_prime_sym - 
     $      transpose(conjg(A_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "A_prime is Hermitian"
      ELSE
         WRITE(*,*), "A_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(Gamma_prime_sym - 
     $      transpose(conjg(B_prime_sym))) < tolerance)
      IF (is_hermitian) then
         WRITE(*,*), "B_prime is Hermitian"
      ELSE
         WRITE(*,*), "B_prime is not Hermitian"
      END IF
      is_hermitian = all(abs(B_prime_sym - 
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