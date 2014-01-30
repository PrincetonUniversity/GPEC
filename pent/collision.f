c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     Collision operators
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. collision_mod
c     A. nu_x
c     B. nu_all
c-----------------------------------------------------------------------
c     subprogram 0. collision_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE collision_mod
      USE pentglobal_mod

      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     Function A. nu_x
c     Collision operator for single energy
c-----------------------------------------------------------------------
      FUNCTION nu_x(nutype,s,sigma,l,x)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8) :: nu_x
      INTEGER, INTENT(IN) :: s,sigma
      REAL(r8), INTENT(IN) :: l
      COMPLEX(r8), INTENT(IN) :: x
      CHARACTER(32), INTENT(IN) :: nutype       
      INTEGER :: trap
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      trap = 1-sigma
      SELECT CASE (nutype)
         CASE ("zero")
            nu_x = 0.0
         CASE ("small")
            nu_x = 1E-5*kin%f(5)
         CASE ("krook")
            nu_x = kin%f(s+6)
     $           /(sigma+trap*2.0*kin%f(9))
         CASE DEFAULT ! "harmonic"
            IF(x==0)THEN        ! (0+0i)^-3 = (Nan,Nan) -> causes error
               nu_x=0.0**(-1.5) ! 0^-3 = Inf -> smooth arithmatic
            ELSE
               nu_x = kin%f(s+6)*(1+0.25*l*l)*x**(-1.5)
     $              /(sigma+trap*2.0*kin%f(9))
            ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION nu_x

c-----------------------------------------------------------------------
c     Function B. nu_all
c     Collision operator for x array of dimension 0:nx
c-----------------------------------------------------------------------
      FUNCTION nu_all(nutype,s,sigma,l,x)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      COMPLEX(r8), DIMENSION(0:nx) :: nu_all
      INTEGER, INTENT(IN) :: s,sigma
      REAL(r8), INTENT(IN) :: l
      COMPLEX(r8), DIMENSION(0:nx), INTENT(IN) :: x
      CHARACTER(32), INTENT(IN) :: nutype       
      INTEGER :: trap
c-----------------------------------------------------------------------
c     Calculations.
c-----------------------------------------------------------------------
      trap = 1-sigma
      SELECT CASE (nutype)
         CASE ("zero")
            nu_all = 0.0
         CASE ("small")
            nu_all = 1E-5*kin%f(5)
         CASE ("krook")
            nu_all = kin%f(s+6)
     $           /(sigma+trap*2.0*kin%f(9))
         CASE ("harmonic")
            nu_all = kin%f(s+6)*(1+0.25*l*l)*x(:)**(-1.5)
     $           /(sigma+trap*2.0*kin%f(9))
         CASE DEFAULT           ! harmonic
            nu_all = kin%f(s+6)*(1+0.25*l*l)*x(:)**(-1.5)
     $           /(sigma+trap*2.0*kin%f(9))
      END SELECT
c-----------------------------------------------------------------------
c     Terminate Function.
c-----------------------------------------------------------------------
      END FUNCTION nu_all


c-----------------------------------------------------------------------
c     terminate module.
c-----------------------------------------------------------------------
      END MODULE collision_mod
