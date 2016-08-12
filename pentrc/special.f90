!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module special
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Collection of special functions.
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

    use utilities, only: get_free_file_unit
    use params, only: r8, pi, twopi

    use spline_mod, only: spline_type,spline_fit,spline_int, spline_alloc ! in equil library
    
    implicit none
    
    contains
    
    !=======================================================================
    function ellipk(m)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Complete elliptical integral of the first kind using polynomial
    !   approximation taken from
    !   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680022044.pdf
    !
    !   Range is 0<=m<1, error is expected to be < 2E-8.
    !
    !*ARGUMENTS:
    !   k : real
    !       normalized pitch angle variable (0 to 1)
    !
    !*RETURNS:
    !   First elliptical integral of the first kind.
    !
    !-----------------------------------------------------------------------  
        implicit none
        real(r8) :: ellipk
        real(r8), intent(in) :: m
        real(r8) :: m1
        
        real(r8), dimension(5) :: &
            a = (/ 1.38629436112, 0.09666344259, 0.03590092383, 0.03742563713, 0.01451196212 /),&
            b = (/ 0.50000000000, 0.12498593597, 0.06880248576, 0.03328355346, 0.00441787012 /)
            
        ! polynomial approximation
        m1 = 1.0-m
        ellipk = (a(1)+a(2)*m1+a(3)*m1**2+a(4)*m1**3+a(5)*m1**4) &
            +(b(1)+b(2)*m1+b(3)*m1**2+b(4)*m1**3+b(5)*m1**4)*log(1/m1)
    end function ellipk
    
    !=======================================================================
    function ellipe(m)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Complete elliptical integral of the second kind using polynomial
    !   approximation taken from
    !   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680022044.pdf
    !
    !   Range is 0<=m<1, error is expected to be < 2E-8.
    !
    !*ARGUMENTS:
    !   k : real
    !       normalized pitch angle variable (0 to 1)
    !
    !*RETURNS:
    !   Second elliptical integral of the first kind.
    !
    !-----------------------------------------------------------------------  
        implicit none
        real(r8) :: ellipe
        real(r8), intent(in) :: m
        real(r8) :: m1
        
        real(r8), dimension(4) :: &
            a = (/ 0.44325141463, 0.06260601220, 0.04757383546, 0.01736506451 /),&
            b = (/ 0.24998368310, 0.09200180037, 0.04069697526, 0.00526449639 /)
            
        ! polynomial approximation
        m1 = 1.0-m
        ellipe = (1+a(1)*m1+a(2)*m1**2+a(3)*m1**3+a(4)*m1**4) &
            +(b(1)*m1+b(2)*m1**2+b(3)*m1**3+b(4)*m1**4)*log(1/m1)
    end function ellipe
    
    !=======================================================================
    subroutine set_ellip
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Pre-calculate elliptic integrals of the first and second kind 0->1.
    !
    !*ARGUMENTS:
    !
    !----------------------------------------------------------------------- 
        integer :: k,fout_unit
        integer, parameter :: nk = 1000
        real(r8), dimension(1:nk) :: kappa2,ek,ee
    
        ! calculations
        print *, "Calculating Elliptic integrals: kappa2 [0,1]"
        do k = 1,nk
            kappa2(k) = k/(nk+1.0)
            ek(k) = ellipk(kappa2(k))
            ee(k) = ellipe(kappa2(k))
        enddo
        print *,"Writting results to unformatted ellip.bin"
        fout_unit = get_free_file_unit(-1)
        open(unit=fout_unit,file="ellip.bin",status="unknown",&
            position="rewind",form="unformatted")
        write(fout_unit) nk
        write(fout_unit) kappa2(:)
        write(fout_unit) ek(:)
        write(fout_unit) ee(:)
        close(fout_unit)
    
        print *,"Writting results to ascii ellip.out"
        open(unit=fout_unit,file="ellip.out",status="unknown")
        write(fout_unit,'(1x,3(1x,a16))') "m","ellipk","ellipe"
        do k = 1,nk
            write(fout_unit,'(1x,3(1x,es16.8e3))') kappa2(k),ek(k),ee(k)
        enddo
        close(fout_unit)
        
        return
    end subroutine set_ellip
    
    !=======================================================================
    function fynml(k,mnql)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Special function from [Park, Boozer, Menard, Phys. Rev. Let., 2009]
    !   Eq. (13).
    !
    !*ARGUMENTS:
    !   k : real
    !       normalized pitch angle variable
    !   mnql : real
    !       m-nq-l where m is the poloidal mode, n is the toroidal mode,
    !       q is the safety factor and l is the bounce harmonic.
    !*RETURNS:
    !   Special function F^y_mnl
    !
    !-----------------------------------------------------------------------  
        implicit none
        real(r8) :: fynml
        real(r8), intent(in) :: k,mnql
        integer :: ntheta,i
        real(r8) :: thetab
        type(spline_type) :: fspl         
    
        ! find bounce point
        thetab = 2*asin(k)
        ! set spline length based on bounce point
        ntheta = 100+20*int(abs(mnql*2*thetab/twopi))
        call spline_alloc(fspl,ntheta,1)
        ! fill spline
        fspl%xs(:) = (/(i,i=1,ntheta+1)/)/(ntheta+2.0)
        fspl%xs(:) = cos(pi*fspl%xs(:)+pi)*thetab
        fspl%fs(:,1) = cos(mnql*fspl%xs)/sqrt(k*k-sin(fspl%xs/2)*sin(fspl%xs/2))
        ! integrate
        call spline_fit(fspl,'extrap')
        call spline_int(fspl)
        fynml = fspl%fsi(fspl%mx,1)
    end function fynml

    
    !=======================================================================
    subroutine set_fymnl
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Pre-calculate (kappa,m-nq-l) space matrix of the special function 
    !   from [Park, Boozer, Menard, Phys. Rev. Let., 2009] Eq. (13)
    !
    !*ARGUMENTS:
    !
    !----------------------------------------------------------------------- 
        integer :: k,t,fout_unit
        integer, parameter :: nk = 500
        integer, parameter :: nt = 500
        real(r8), dimension(0:nk) :: kappa
        real(r8), dimension(-nt:nt) :: mnql
        real(r8), dimension(0:nk,-nt:nt) :: fkmnql
    
        ! calculations
        print *, "Calculating F^1/2_mnl: kappa [0,1], m-nq-l [-100,100]"
        !call timer(1)
        do t = -nt,nt
           mnql(t)  = t/(nt/100.0)
        enddo
        do k = 0,nk
           kappa(k) = (k+1.0)/(nk+2.0)
           do t = -nt,nt
              fkmnql(k,t) = fynml(kappa(k),mnql(t))
           enddo
        enddo
        !call timer(1)
        print *,"Writting bin map to fkmnql.bin"
        fout_unit = get_free_file_unit(-1)
        open(unit=fout_unit,file="fkmnql.bin",status="unknown",&
            position="rewind",form="unformatted")
        write(fout_unit) nk,2*nt
        write(fout_unit) kappa(:)
        write(fout_unit) mnql(:)
        write(fout_unit) fkmnql(:,:)
        close(fout_unit)
    
        print *,"Writting ascii map to fkmnql.out"
        open(unit=fout_unit,file="fkmnql.out",status="unknown")
        write(fout_unit,*)"LARGE ASPECT RATIO F^1/2_mnl: "
        write(fout_unit,*)
        write(fout_unit,*)"nk = ",nk,"nt = ",2*nt
        write(fout_unit,*)
        write(fout_unit,'(1x,3(1x,a16))') "kappa","m-nq-l","f_mnl"
        do k = 0,nk
            do t = -nt,nt
                write(fout_unit,'(1x,3(1x,es16.8e3))') kappa(k),mnql(t),fkmnql(k,t)
           enddo
        enddo
        close(fout_unit)
        
        return
    end subroutine set_fymnl
    
    
    !=======================================================================
    subroutine wofz (xi, yi, u, v, flag)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
    !   THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !   VOL. 16, NO. 1, PP. 47.
    !
    !   ** downloaded from http://www.netlib.org/toms/680 **
    !
    ! GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
    ! THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
    ! WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
    ! MEANS SQRT(-1).
    ! THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
    ! IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
    ! DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
    ! OF THE FUNCTION.
    ! ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
    !
    ! THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
    !    RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
    !               RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
    !               IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
    !               FLOATING-POINT ARITHMETIC
    !    RMAXEXP  = LN(RMAX) - LN(2)
    !    RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
    !               GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
    ! THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
    ! BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
    !
    ! PARAMETER LIST
    !    XI     = REAL      PART OF Z
    !    YI     = IMAGINARY PART OF Z
    !    U      = REAL      PART OF W(Z)
    !    V      = IMAGINARY PART OF W(Z)
    !    FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
    !             OCCUR OR NOT; TYPE LOGICAL;
    !             THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
    !             MEANING :
    !             FLAG=.FALSE. : NO ERROR CONDITION
    !             FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
    !                            BECOMES INACTIVE
    ! XI, YI      ARE THE INPUT-PARAMETERS
    ! U, V, FLAG  ARE THE OUTPUT-PARAMETERS
    !
    ! FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
    !
    ! THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
    ! PUT TO 0 UPON UNDERFLOW;
    !
    ! REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
    ! THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
    !----------------------------------------------------------------------- 
        
        implicit double precision (a-h, o-z)
        implicit integer (i-n)

        logical a, b, flag
        real(r8), parameter :: factor   = 1.12837916709551257388d0,&
                   rmaxreal = 0.5d+154,&
                   rmaxexp  = 708.503061461606d0,&
                   rmaxgoni = 3.53711887601422d+15
        
        flag = .false.
        c = 0
        
        xabs = dabs(xi)
        yabs = dabs(yi)
        x    = xabs/6.3
        y    = yabs/4.4
        
        !    the following if-statement protects
        !    qrho = (x**2 + y**2) against overflow
        
        if ((xabs.gt.rmaxreal).or.(yabs.gt.rmaxreal)) then
          flag = .true.
          return
        endif
        
        qrho = x**2 + y**2
        
        xabsq = xabs**2
        xquad = xabsq - yabs**2
        yquad = 2*xabs*yabs
        
        a     = qrho.lt.0.085264d0
        
        if (a) then
        
        ! if (qrho.lt.0.085264d0) then the faddeeva-function is evaluated
        ! using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
        ! n is the minimum number of terms needed to obtain the required
        ! accuracy
        
          qrho  = (1-0.85*y)*dsqrt(qrho)
          n     = idnint(6 + 72*qrho)
          j     = 2*n+1
          xsum  = 1.0/j
          ysum  = 0.0d0
          do i=n, 1, -1
            j    = j - 2
            xaux = (xsum*xquad - ysum*yquad)/i
            ysum = (xsum*yquad + ysum*xquad)/i
            xsum = xaux + 1.0/j
          enddo
          u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0
          v1   =  factor*(xsum*xabs - ysum*yabs)
          daux =  dexp(-xquad)
          u2   =  daux*dcos(yquad)
          v2   = -daux*dsin(yquad)
        
          u    = u1*u2 - v1*v2
          v    = u1*v2 + v1*u2
        
        else
        
        ! if (qrho.gt.1.o) then w(z) is evaluated using the laplace
        ! continued fraction
        ! nu is the minimum number of terms needed to obtain the required
        ! accuracy
        
        ! if ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
        ! by a truncated taylor expansion, where the laplace continued fraction
        ! is used to calculate the derivatives of w(z)
        ! kapn is the minimum number of terms in the taylor expansion needed
        ! to obtain the required accuracy
        ! nu is the minimum number of terms of the continued fraction needed
        ! to calculate the derivatives with the required accuracy
        
          if (qrho.gt.1.0) then
            h    = 0.0d0
            kapn = 0
            qrho = dsqrt(qrho)
            nu   = idint(3 + (1442/(26*qrho+77)))
          else
            qrho = (1-y)*dsqrt(1-qrho)
            h    = 1.88*qrho
            h2   = 2*h
            kapn = idnint(7  + 34*qrho)
            nu   = idnint(16 + 26*qrho)
          endif
        
          b = (h.gt.0.0)
        
          if (b) qlambda = h2**kapn
        
          rx = 0.0
          ry = 0.0
          sx = 0.0
          sy = 0.0
        
          do n=nu, 0, -1
            np1 = n + 1
            tx  = yabs + h + np1*rx
            ty  = xabs - np1*ry
                !  = 0.5/(tx**2 + ty**2)
            rx  = c*tx
            ry  = c*ty
            if ((b).and.(n.le.kapn)) then
              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda/h2
            endif
          enddo
          if (h.eq.0.0) then
            u = factor*rx
            v = factor*ry
          else
            u = factor*sx
            v = factor*sy
          end if
        
          if (yabs.eq.0.0) u = dexp(-xabs**2)
        
        end if
        
        ! evaluation of w(z) in the other quadrants
        
        if (yi.lt.0.0) then
        
          if (a) then
            u2    = 2*u2
            v2    = 2*v2
          else
            xquad =  -xquad
        
        !        the following if-statement protects 2*exp(-z**2)
        !        against overflow
        
            if ((yquad.gt.rmaxgoni).or.&
                (xquad.gt.rmaxexp)) then
                flag = .true.
                return
            endif
        
            w1 =  2*dexp(xquad)
            u2  =  w1*dcos(yquad)
            v2  = -w1*dsin(yquad)
          end if
        
          u = u2 - u
          v = v2 - v
          if (xi.gt.0.0) v = -v
        else
          if (xi.lt.0.0) v = -v
        end if
        
        return
        
        flag = .true.
        return

    end subroutine wofz

      
end module special