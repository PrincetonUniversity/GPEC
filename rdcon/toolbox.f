c-----------------------------------------------------------------------
c     file toolbox.f.
c     tools for different purpose of reuse code
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. toolbox_mod.
c     1. toolbox_hermite.
c     2. toolbox_preflatcof.
c-----------------------------------------------------------------------
c     subprogram 0. toolbox_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE toolbox_mod
      USE rdcon_mod

      IMPLICIT NONE

      TYPE :: toolbox_hermite_type
      REAL(r8), DIMENSION(0:3) :: pb,qb
      END TYPE toolbox_hermite_type

      TYPE :: toolbox_interval_type
      INTEGER :: int_type
      REAL(r8), DIMENSION(2) :: x
      END TYPE toolbox_interval_type

      TYPE(toolbox_hermite_type) :: tb_hermite

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. toolbox_hermite.
c     computes hermite cubic basis functions and their derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE toolbox_hermite(x,x0,x1,hermite)

      REAL(r8),INTENT(IN) :: x,x0,x1
      TYPE(toolbox_hermite_type),INTENT(INOUT) :: hermite 
      
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
      END SUBROUTINE toolbox_hermite
c-----------------------------------------------------------------------
c     subprogram 2. toolbox_preflatcof.
c     coefficients for flattening pressure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE toolbox_preflatcof(psi,dpsi_fac,dpsi1_fac,fac,cof,cof1)
      REAL(r8),INTENT(IN) :: psi,dpsi_fac,dpsi1_fac
      REAL(r8),INTENT(OUT) :: cof, cof1
      REAL(r8),INTENT(IN) :: fac
      INTEGER ising
      REAL(r8), DIMENSION(2) :: dx
      REAL(r8), DIMENSION(3) :: x
      REAL(r8), DIMENSION(0:msing+1) :: singpsi
      TYPE(toolbox_interval_type),DIMENSION(4*msing+1) :: tb_intvl
c-----------------------------------------------------------------------
c     generate the flatten interval.
c-----------------------------------------------------------------------
      singpsi(0)=psilow
      DO ising=1,msing
         singpsi(ising)=sing(ising)%psifac
      ENDDO
      singpsi(ising)=psihigh
      tb_intvl(1)%x(1)=psilow
      tb_intvl(1)%int_type=1
      DO ising=1,msing
         x=singpsi(ising-1:ising+1)
         dx(1)=(x(2)-x(1))*dpsi_fac
         dx(2)=(x(3)-x(2))*dpsi_fac
         IF (dx(2)<dx(1)) dx(1)=dx(2)
         tb_intvl(ising*4-3)%x(2)=x(2)-dx(1)

         tb_intvl(ising*4-2)%x(1)=x(2)-dx(1)
         tb_intvl(ising*4-2)%x(2)=x(2)-dx(1)*dpsi1_fac
         tb_intvl(ising*4-2)%int_type=2

         tb_intvl(ising*4-1)%x(1)=x(2)-dx(1)*dpsi1_fac
         tb_intvl(ising*4-1)%x(2)=x(2)+dx(1)*dpsi1_fac
         tb_intvl(ising*4-1)%int_type=4


         tb_intvl(ising*4)%x(1)=x(2)+dx(1)*dpsi1_fac
         tb_intvl(ising*4)%x(2)=x(2)+dx(1)
         tb_intvl(ising*4)%int_type=3

         tb_intvl(ising*4+1)%x(1)=x(2)+dx(1)
         tb_intvl(ising*4+1)%int_type=1
      ENDDO
      tb_intvl(ising*4-3)%x(2)=psihigh

      DO ising=1,4*msing+1
         IF (psi>=tb_intvl(ising)%x(1).AND.psi<=tb_intvl(ising)%x(2))
     $   EXIT
      ENDDO
      IF (ising>4*msing+1) THEN
         CALL program_stop("toolbox_preflatcof: can't find inteval.")
      ENDIF
      x(1)=tb_intvl(ising)%x(1)
      x(2)=tb_intvl(ising)%x(2)
      CALL toolbox_hermite(psi,x(1),x(2),tb_hermite)
      SELECT CASE(tb_intvl(ising)%int_type)
      CASE(1)
         cof=1
         cof1=0
      CASE(2)
         cof=tb_hermite%pb(0)+fac*tb_hermite%pb(2)
         cof1=tb_hermite%qb(0)+fac*tb_hermite%qb(2)
      CASE(3)
         cof=fac*tb_hermite%pb(0)+tb_hermite%pb(2)
         cof1=fac*tb_hermite%qb(0)+tb_hermite%qb(2)
      CASE(4)
         cof=fac
         cof1=0
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE toolbox_preflatcof
      END MODULE toolbox_mod
