!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module diagnostics
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !     A subroutine that calls many argument free test functions
    !     to check consistency of PENTRC components.
    !
    !*PUBLIC MEMBER FUNCTIONS:
    !
    !*PUBLIC DATA MEMBERS:
    ! 
    !*REVISION HISTORY:
    !     2014.03.05 -Logan- initial writing.
    !
    !-----------------------------------------------------------------------
    ! AUTHOR: Logan
    ! EMAIL: nlogan@pppl.gov
    !-----------------------------------------------------------------------
    
    use params, only: r8, pi,xj
    use special, only: wofz,ellipe,ellipk
    use spline_mod, only: spline_type,spline_alloc,spline_fit
    use cspline_mod, only: cspline_type,cspline_alloc,cspline_fit
    use energy_integration, only : xintgrl_lsode,output_energy_netcdf,xf0type,xnutype
    use pitch_integration, only: lambdaintgrl_lsode, output_pitch_netcdf
    
    implicit none
    private
    
    public diagnose_all,xtest
    
    contains
    
    
    subroutine diagnose_all()
        !----------------------------------------------------------------------- 
        !*DESCRIPTION: 
        !     Run all the diagnostics.
        !
        !*ARGUMENTS:
        !
        !-----------------------------------------------------------------------
        
        implicit none
        
        call xtest
        call pitchtest
    
    end subroutine diagnose_all
    



    subroutine xtest()
        !----------------------------------------------------------------------- 
        !*DESCRIPTION: 
        !     Run energy integration.
        !
        !*ARGUMENTS:
        !
        !-----------------------------------------------------------------------
        implicit none
        
        ! declare internal variables      
        real(r8) :: wn=0.0, wt=0.0, we=1e4, wd=1e2, wb=1e4, psi=0, lambda=0, &
            nuk=1e2, lnq = 0
        integer :: n = 1, ell = 0
        real(r8) :: zr,zi,wr,wi
        complex(r8) :: xint,analytic,omegan,omegastara,omegastarb,zfun,z
        logical :: wflag
        
        print *,"------------------------------------------------------"
        print *,"Energy integration diagnostics"
        print *,"------------------------------------------------------"
        
        ! We have a analytic solution for the CGL limit: 15 sqrt(pi)/8
        xf0type = "cgl"
        xint = xintgrl_lsode(wn,wt,we,wd,wb,nuk,ell,lnq,n,psi,lambda,"fcgl",op_record=.true.)
        call output_energy_netcdf(n,op_label="diagnostic_cgl")
        analytic = -xj*15.0*sqrt(pi)/8.0
        print *,"CGL limit: "
        print *," -> numerical  ",xint
        print *," -> analytical ",analytic
        print *,''
        
        ! We have a analytic solution for the Krook limit 
        ! ** note I changed the sign of nuk compared to MISK mdc2 document **
        xf0type = "maxwellian"
        xnutype = "krook"
        xint = xintgrl_lsode(wn,wt,we,wd,wb,nuk,ell,lnq,n,psi,lambda,"fcgl",op_record=.true.)
        call output_energy_netcdf(n,op_label="diagnostic_krook")
        omegan = (n*we+xj*nuk)/(n*wd)
        omegastara = (n*wn-1.5*n*wt-xj*nuk)/(n*wd)
        omegastarb = wt/wd
        ! ** need plasma displersion function (see mdc2 benchmark with MISK) **
        ! from http://w3.pppl.gov/~hammett/comp/src/wofz_readme.html we have
        z = xj*sqrt(omegan)
        zr = real(z)
        zi = aimag(z)
        !zfun = xj*sqrt(pi)*exp(-z**2)*erfc(-xj*z)
        call wofz(zr,zi,wr,wi,wflag)
        zfun = wr+xj*wi
        analytic = -xj*( 15*sqrt(pi)*omegastarb/8 &
            +2*sqrt(pi)*(omegan+omegastara-omegan*omegastarb)&
            *(3.0/8-0.25*omegan+0.5*omegan**2 +xj*0.5*omegan**2.5*zfun) )
        print *,"Krook limit:"
        if(.not.wflag)then ! dispersion function worked
            print *," -> numerical  ",xint
            print *," -> analytical ",analytic
            print *,''
        else
            print *,' -> ERROR: dispersion function overflow error'
            print *,''
        endif
        
        ! An easier analytic solution for the Krook limit is the wd=0 case
        ! ** note I changed the sign of nuk compared to MISK mdc2 document **
        xf0type = "maxwellian"
        xnutype = 'krook'
        wd = 0
        xint = xintgrl_lsode(wn,wt,we,wd,wb,nuk,ell,lnq,n,psi,lambda,"fcgl",op_record=.true.)
        call output_energy_netcdf(n,op_label="diagnostic_wd0")
        analytic = -xj*(15*sqrt(pi)/8) * (n*(wn+2*wt + we)) / (n*we+xj*nuk)
        print *,"Krook, omegaD=0 limit:"
        print *," -> numerical  ",xint
        print *," -> analytical ",analytic
        
        print *,"------------------------------------------------------"       
    end subroutine xtest




    
    
    subroutine pitchtest()
        !----------------------------------------------------------------------- 
        !*DESCRIPTION: 
        !   Run pitch integration.
        !   -> The most simplified "analytic" limit would be the form in
        !       [Park, Boozer, Menard, PRL 2009], which has wb,wd constant in
        !       Lambda... bu tthis still requires numerical integration.
        !
        !*ARGUMENTS:
        !
        !-----------------------------------------------------------------------
        implicit none
        
        ! declare internal variables      
        integer, parameter :: nlmda = 1000
        real(r8) :: wn=0.0, wt=0.0, we=1e4, wd=1e2,wb=1e4,nuk=1e2, &
            bobmax,epsr,lmax,q,psi,rex,imx
        integer :: i,n=1,l=0
        real(r8), dimension(0:nlmda) :: kappa2,lmda
        complex(r8), dimension(:), allocatable :: lint
        type(cspline_type) :: bdf_spl
        type(spline_type) :: turns
        
        print *,"------------------------------------------------------"
        print *,"Lambda integration diagnostics"
        print *,"------------------------------------------------------"
        
        ! use a cylindrical plasma B = B0(1-epsr*cos(theta))
        psi = 1.0
        epsr = 0.33
        bobmax = 1/(1+epsr)
        lmax = 1/(1-epsr)
        q = 1.0
        wd = q**3*wd
        rex = 1.0
        imx = 1.0
        !lmda = (/ bobmax+i*(lmax-bobmax)/nlmda,i=0,nlmda /) ! trapped only
        do i=0,nlmda
            lmda(i) = bobmax+i*(lmax-bobmax)/nlmda
        enddo
        lmda(0) = (lmda(1)+lmda(0))/2.0 ! prevent elliptic integral errors
        
        ! use a general distribution and simple collision op.
        xf0type = "maxwellian"
        xnutype = "krook"
        
        ! Easy case is wb,wd constant in Lambda
        ! Give fl some Lambda dependence so its not just a test of many xintgrls
        call cspline_alloc(bdf_spl,nlmda,3)
        bdf_spl%xs(:) = lmda(:)
        bdf_spl%fs(:,1) = wb
        bdf_spl%fs(:,2) = wd
        bdf_spl%fs(:,3) = lmda(:)**2
        call cspline_fit(bdf_spl,'extrap')
        ! placeholder for turn point info
        call spline_alloc(turns,nlmda,6)
        turns%xs(:) = lmda(:)
        turns%fs = 0
        call spline_fit(turns,'extrap')
        lint = lambdaintgrl_lsode(wn,wt,we,nuk,bobmax,epsr,q,bdf_spl,l,n, &
            rex,imx,psi,turns,'rlar',op_record=.true.)
        print *,"Reduced LAR approximations: "
        print *," -> numerical  ",lint
        
        ! Circular large aspect ratio freqs
        kappa2 = (1-lmda*(1-epsr))/(2*epsr*lmda)
        do i=0,nlmda
            !print *,i,lmda(i),kappa2(1)
            bdf_spl%fs(i,1) = wb*2.0/ellipk(kappa2(i))
            bdf_spl%fs(i,2) = (q**3*lmda(i)/epsr)*(ellipe(kappa2(i))/ellipk(kappa2(i))-0.5)
        enddo
        call cspline_fit(bdf_spl,'extrap')
        lint = lambdaintgrl_lsode(wn,wt,we,nuk,bobmax,epsr,q,bdf_spl,l,n, &
            rex,imx,psi,turns,'clar',op_record=.true.)
        call output_energy_netcdf(n,op_label="diagnostic_pitch")
        call output_pitch_netcdf(n,op_label="diagnostic_pitch")
        print *,"Circular LAR approximations: "
        print *," -> numerical  ",lint
        
        print *,"------------------------------------------------------"       
    end subroutine pitchtest

    
end module diagnostics
