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
    !     2014.03.06 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
        
    implicit none
    
    include "version.inc"

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

    ! code
    integer, parameter :: npsi_out = 30, nlambda_out = 5, nell_out = 3, &
        nmethods = 18, ngrids = 3
    character(4), dimension(nmethods), parameter :: methods = (/ &
        "fgar","tgar","pgar","rlar","clar","fcgl",&
        "fwmm","twmm","pwmm","ftmm","ttmm","ptmm",&
        "fkmm","tkmm","pkmm","frmm","trmm","prmm" /)
    character(5), dimension(ngrids), parameter :: grids = (/'lsode','equil','input'/)
    character(64), dimension(nmethods), parameter :: docs = (/&
        "Full general-aspect-ratio calculation                       ",&
        "Trapped particle general-aspect-ratio calculation           ",&
        "Passing particle general-aspect-ratio calculation           ",&
        "Trapped particle large-aspect-ratio calculation             ",&
        "Trapped particle cylindrical large-aspect-ratio calculation ",&
        "Fluid Chew-Goldberger-Low calculation                       ",&
        "Full    energy calculation using MXM euler lagrange matrix  ",&
        "Trapped energy calculation using MXM euler lagrange matrix  ",&
        "Passing energy calculation using MXM euler lagrange matrix  ",&
        "Full    torque calculation using MXM euler lagrange matrix  ",&
        "Trapped torque calculation using MXM euler lagrange matrix  ",&
        "Passing torque calculation using MXM euler lagrange matrix  ",&
        "Full    MXM euler lagrange energy matrix norm calculation   ",&
        "Trapped MXM euler lagrange energy matrix norm calculation   ",&
        "Passing MXM euler lagrange energy matrix norm calculation   ",&
        "Full    MXM euler lagrange torque matrix norm calculation   ",&
        "Trapped MXM euler lagrange torque matrix norm calculation   ",&
        "Passing MXM euler lagrange torque matrix norm calculation   " &
        /)

end module params