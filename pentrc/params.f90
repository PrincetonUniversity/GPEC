!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module params
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Collection of basic contants.
    !
    !*PUBLIC MEMBER FUNCTIONS: 
    !
    !*PUBLIC DATA MEMBERS:
    !   All variables are public
    !
    !*REVISION HISTORY:
    !     2014.03.06 -Logan- initial writting. 
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
        
    implicit none
    
    ! fortran
    integer, parameter :: &
        r4=selected_real_kind(6,37), & ! originally in spline_mod->local_mod
        r8=selected_real_kind(13,307)
    real(r8) :: &
        inf = huge(0.0_r8)
        
    ! physics
    real(r8), parameter :: &
        mp=1.672614e-27,   &
        me=9.1091e-31,     &
        e=1.6021917e-19,   &
        ev=e
    
    ! math
    real(r8), parameter ::  &
        pi= 3.1415926535897932385_r8, &
        twopi = 2*pi, &
        mu0=4e-7_r8*pi,     &
        rtod=180/pi,        &
        dtor=pi/180
    complex(r8), parameter ::   &
        xj = (0,1)
    
end module params