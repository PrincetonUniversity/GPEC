!@PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT CODE

module utilities
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Collection of utility procedures for general fortran needs.
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
    use params, only: r8,twopi
    use netcdf
    
    implicit none
    
    contains
    
    !=======================================================================
    function get_free_file_unit (lu_max)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Find a unit number that is not in use.
    !
    !*ARGUMENTS:
    !   lu_max : integer.
    !       Maximum unit.
    !
    !*RETURNS:
    !     integer.
    !        Free unit.
    !-----------------------------------------------------------------------
        implicit none
        integer get_free_file_unit
        integer lu_max,  lu, m, iostat
        logical opened
        
        m = lu_max
        if(m < 1) m = 97
        
        do lu = m,1,-1
           inquire (unit=lu, opened=opened, iostat=iostat)
           if (iostat.ne.0) cycle
           if (.not.opened) exit
        end do
        get_free_file_unit = lu
        
        return
    end function get_free_file_unit
    
    !=======================================================================
    function to_upper(strIn) result(strOut)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Capitalize a string.
    !   Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    !   Original author: Clive Page
    !
    !*ARGUMENTS:
    !   lu_max : integer.
    !       Maximum unit.
    !
    !*RETURNS:
    !     integer.
    !        Free unit.
    !-----------------------------------------------------------------------
    
         implicit none
    
         character(len=*), intent(in) :: strIn
         character(len=len(strIn)) :: strOut
         integer :: i,j
    
         do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j>= iachar("a") .and. j<=iachar("z") ) then
                   strOut(i:i) = achar(iachar(strIn(i:i))-32)
              else
                   strOut(i:i) = strIn(i:i)
              end if
         end do
    
    end function to_upper

    !=======================================================================
    function btoi(bool) result(ibool)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Convert logical variable to integer.
    !
    !*ARGUMENTS:
    !   bool : logical.
    !       True of False.
    !
    !*RETURNS:
    !     integer.
    !        0 or 1.
    !-----------------------------------------------------------------------

         implicit none

         logical, intent(in) :: bool
         integer :: ibool

         ibool = 0
         if(bool) ibool=1

    end function btoi

    !=======================================================================
    function itob(ibool) result(bool)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !   Convert integer variable to logical.
    !
    !*ARGUMENTS:
    !   ibool : integer.
    !       Any integer variable
    !
    !*RETURNS:
    !     logical.
    !        True if ibool>0.
    !-----------------------------------------------------------------------

         implicit none

         integer, intent(in) :: ibool
         logical :: bool

         bool = .false.
         if(ibool>0) bool=.true.

    end function itob

    !=======================================================================
    subroutine progressbar (j,jstart,jstop,op_step,op_percent)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Print progress bar for do loop every 10% of iterations.
    !
    !*ARGUMENTS:
    !   j : integer, in
    !       Cuurent iteration variable.
    !   jstart : integer, in
    !       Start of do loop
    !   jstop : integer, in
    !       End of do loop
    !*OPTIONAL ARGUMENTS:
    !   op_tep : integer, optional in
    !       Do loop increment (default 1)
    !   op_percent : integer, optional in
    !       Prints progress every percent of the iterations (default 10).
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: j,jstart,jstop
        integer, optional :: op_step,op_percent
    
        integer :: iteration,iterations,k,done,jstep,percent
        character(64) :: bar
        
        ! set up
        jstep = 1
        percent = 10
        if(present(op_percent)) percent = op_percent
        if(present(op_step)) jstep = op_step
        iteration = (j-jstart)/jstep + 1
        iterations= (jstop-jstart)/jstep + 1
        bar = "  |----------| ???% iterations complete"
        
        if(iteration==1 .or. mod(iteration,iterations*percent/100)==0)then
            done = iteration/(iterations/10)
            write(unit=bar(16:18),fmt="(i3)") 10*done
            do k=1, done
               bar(3+k:3+k)="o"
            enddo  
            print *,trim(bar)
        endif
        
        return
    end subroutine progressbar
    
    !=======================================================================
    subroutine timer (mode,unit)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Handles machine-dependent timing statistics.
    !
    !*ARGUMENTS:
    !   mode : integer, in
    !       mode = 0 starts timer
    !       mode < 0 writes time from last start and restarts timer
    !       mode > 0 writes time from last start
    !
    !*OPTIONAL ARGUMENTS:
    !   unit : integer, in
    !       Output written to this unit (default to terminal)
    !
    !-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: mode
        integer, intent(in), optional :: unit
    
        character(10) :: date,time,zone
        integer, dimension(8) :: values
        real(4), save :: seconds
        integer :: hrs,mins,secs
        
        ! write split
        if(mode/=0)then
            call date_and_time(date,time,zone,values)
            seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3-seconds
            secs = int(seconds)
            hrs = secs/(60*60)
            mins = (secs-hrs*60*60)/60
            secs = secs-hrs*60*60-mins*60
            if(present(unit))then
                write(unit,"(1x,a,1p,e10.3,a)")"total cpu time = ",seconds," seconds"
            else
                if(hrs>0)then
                    print *,"total cpu time = ",hrs," hours, ", &
                        mins," minutes, ",secs," seconds"
                elseif(mins>0)then
                    print *,"total cpu time = ", &
                        mins," minutes, ",secs," seconds"
                else
                    print *,"total cpu time = ",secs," seconds"
                endif
            endif
        endif
        ! start timmer
        if(mode<=0)then
            call date_and_time(date,time,zone,values)
            seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
        endif
        
        return
    end subroutine timer
    
    
    !=======================================================================
    subroutine splitstring(str,words,sep,maxsplit,debug)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Split string into substrings using sep as delimiter.
    !
    !*ARGUMENTS:
    !    str : character(*) (in)
    !       Original string.
    !   words : character 1D array, allocatable (inout)
    !       List of substrings, length 32.
    !       -> If allocated, then filled with words from index 1
    !       -> Allocated to be length of words.
    !   sep : character(*)
    !       Deliminator seperating words in original string
    !
    !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        character(*), intent(in) :: str,sep
        integer, optional :: maxsplit
        logical, optional :: debug
        character(32), dimension(:), allocatable :: words
        integer :: pos1 = 1, pos2, pmax, wmax, wmin,n = 0, m=0
        
        ! setup
        if(.not. present(debug)) debug = .false.
        pmax = len(str)
        if(.not. present(maxsplit)) then
            wmax = pmax        
        else
            wmax = maxsplit 
        endif
        pos1 = 1
        n = 0
        m = 0
        
        ! first loop determines number of words
        do
            pos2 = 1
            do while (pos2==1)! .and. pos1<pmax) ! ignores repeated delimiters
                pos2 = index(str(pos1:), sep)
                if(pos2==1) pos1=pos1+1
            enddo
            !print *,pos1,pos2,str(pos1:min(pmax,pos1+pos2))
            if (pos2 == 0) then
                exit
            end if
            pos1 = pos2+pos1
            n = n + 1
        end do
        
        ! form words array
        if(allocated(words))then
            wmin = lbound(words,1)
            wmax = min(wmax,size(words,1))
            wmax=min(n,wmax)
        else
            wmin = 1
            wmax=min(n,wmax)
            allocate(words(wmin:wmax))
        endif
        if(debug) print *,'number of words is ',wmax,' of possible ',n
        words(:) = ' '
        
        ! second loop fills words
        pos1 = 1
        do m=wmin,wmin+wmax-1
            pos2=1
            do while (pos2==1)! .and. pos1<pmax) ! ignores repeated delimiters
                pos2 = index(str(pos1:), sep)
                if(pos2==1) pos1=pos1+1
            enddo
            if (pos2 == 0 .or. m==wmax) then
                words(m) = str(pos1:)
            else
                words(m) = str(pos1:pos1+pos2)
            end if
            pos1 = pos2+pos1
            if(debug) print *,m,words(m)
        end do
        
    end subroutine splitstring
    
    !=======================================================================
    subroutine replacestring(str,subout,subin,debug)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Replace substring with an alternate substring onf equal length
    !
    !*ARGUMENTS:
    !   str : character(*) (inout)
    !       Original string.
    !   subout : character(*) (in)
    !       Substring to be replaced.
    !   subin : character(*) (in)
    !       Substring replacing subout.
    !
    !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        character(*), intent(inout) :: str
        character(*), intent(in) :: subout,subin
        logical, optional :: debug
        integer :: i,ilast,lstr,lout
        
        ! setup
        if(.not. present(debug)) debug = .false.
        lstr = len(str)
        lout = len(subout)-1
        ilast = -lout
        
        if(debug) print *,'  Replacing ',subout,' with ',subin
        if(debug) print *,'  substring length ',lout+1,' in string length ',lstr
        
        ! search and replace
        do i = 1,lstr-lout
            if((str(i:i+lout)==subout) .and. (i-ilast>lout))then
                if(debug) print *,'  found substring ',subout,' at index ',i
                str(i:i+lout) = subin
                ilast = i
            endif
        end do
        
    end subroutine replacestring
    
    !=======================================================================
    subroutine readtable(file,table,titles,verbose,debug)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !   Read ascii file containing a single table.
    !       - Assumes data starts on the first line starting with a number
    !       - Assumes column titles are on the line directly above the data
    !       - Ignores all other lines
    !
    !*ARGUMENTS:
    !    file : character(256) (in)
    !       File path.
    !   table : real, 2D allocatable (out)
    !       The data
    !   titles : character, 1D allocatable (out)
    !       The column headers
    !*OPTIONAL ARGUMENTS:
    !   debug : logical
    !       Print intermediate status messages to terminal.
    !
    !-----------------------------------------------------------------------
    
        implicit none
        ! declare arguments
        character(512), intent(in) :: file
        real(r8), dimension(:,:), allocatable :: table
        character(32), dimension(:), allocatable :: titles,dummies
        logical, optional :: verbose,debug
        ! declare local variables
        integer :: i,j,k,l,in_unit,iostat,startline,endline,ncol
        character(1024) :: line
        character(10) :: numbers = '1234567890'
        character(2) :: signs = '+-'
        logical :: fexists,isnum
        
        ! setup
        if(.not. present(debug)) debug = .false.
        !inquire(file=file,exist=fexists)
        !if(.not. fexists) stop "File does not exist"
        if(verbose) print *, "Reading table from file: "
        if(verbose) print *, '  '//trim(file)
        if(allocated(table)) deallocate(table)
        if(allocated(titles)) deallocate(titles)
        in_unit = get_free_file_unit(-1)
        open(unit=in_unit,file=file,status="unknown")
        if(debug) print *,"-> File opened"
        
        ! find data lines
        startline = 0
        endline = 0
        l = 1
        do while (startline==0 .or. endline==0)
            read(in_unit,'(a)',iostat=k) line
            line = adjustl(line)
            j = 1
            do i=1,2    ! pass over plus/minus signs
                if(line(j:j)==signs(i:i)) j=j+1
            enddo
            isnum = .false.
            do i=1,10   ! see if line starts with number
                if(line(j:j)==numbers(i:i)) isnum = .true.
            enddo
            if(isnum .and. startline==0) startline = l
            if(.not. isnum .and. startline/=0) endline = l-1
            !if(endline/=0) exit
            if(k<0) then
                endline = l-1
                if(startline==0) stop "ERROR: read_table - No data found"
            endif
            !if(debug) print *,"--> l,iostat,start,end = ",l,k,startline,endline
            l=l+1
        enddo
        rewind(in_unit)
        if(debug) print *,"-> Table start,end lines = ",startline,endline
        
        ! determine size of data
        do l=1,startline-1
            read(in_unit,*)
        enddo
        read(in_unit,'(a)') line
        call replacestring(line,char(9),' ',debug=debug) !tabs
        call splitstring(line,dummies,' ',debug=debug)
        ncol = size(dummies,1)
        rewind(in_unit)
        allocate(table(1:1+endline-startline,1:ncol))
        allocate(titles(1:ncol))
        if(debug) print *,' -> Found ',1+endline-startline,' by ',ncol,' table'
        
        ! read titles
        if(startline>1)then
            do l=1,startline-2
                read(in_unit,*)
            enddo
            line(:) = ' '
            read(in_unit,'(a)') line
            call splitstring(line,titles,' ',debug=debug)
        endif
        !if(debug) print *,"-> Table titles are ",titles
        
        ! read table
        do l=1,1+endline-startline
            read(in_unit,*) table(l,:)
        enddo
        close(in_unit)
        if(debug) print *,"Done"
    
    end subroutine readtable

    
    !=======================================================================
    function nunique(array,op_sorted)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Find number of unique elements in an array.
    !
    !*ARGUMENTS:
    !    array : real
    !       Original array.
    !*OPTIONAL ARGUMENTS:
    !   op_sorted : bool
    !       Common values or all unique values are assumed to be grouped.
    !
    !-----------------------------------------------------------------------
    
        IMPLICIT NONE
        !declare function
        integer :: nunique
        ! declare arguments
        real(r8), dimension(:) ::array
        logical, optional :: op_sorted
        ! declare local variables
        logical,dimension(:), allocatable::mask
        real(r8),dimension(:), allocatable::tmp
        integer::i,j,num
        logical :: sorted
        sorted = .FALSE.
        if(present(op_sorted)) sorted=op_sorted
        
        num=size(array)
        if(sorted)then
            nunique=0
            if(array(num)==array(num-1))then !common values are grouped
                do i=num,2,-1
                    if(array(i)/=array(i-1)) nunique = nunique+1
                enddo
            else ! unique values assumed grouped
                do i=2,num
                    if(array(i)==array(1))then
                        nunique=i-1
                        exit
                    endif
                enddo
            endif
        else ! general approach (can be slow for huge arrays)
            ! form mask
            allocate(mask(num)); mask = .true.
            ! check uniqueness of each element
            do i=num,2,-1
               mask(i)=.not.(any(array(:i-1)==array(i)))
            end do
            ! make an index vector
            nunique = count(mask,1)
            deallocate(mask)
        endif
    end function nunique
    
    
    !=======================================================================
    function median(array)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Compute the median of a 1D array, by sorting a copy of the array.
    !
    !*ARGUMENTS:
    !    array : real
    !       Original array.
    !
    !-----------------------------------------------------------------------
        implicit  none
        real(r8), dimension(1:), intent(in) :: array
        real(r8) :: median
        real(r8), dimension(:), allocatable :: temp
        integer :: i,n
        
  
        ! make a copy of array that is sorted min to max
        n = size(array)
        allocate(temp(n))
        temp = sorted(array)
        if (mod(n,2) == 0) then           ! compute the median
            median = (temp(n/2) + temp(n/2+1)) / 2.0
        else
            median = temp(n/2+1)
        end if
        deallocate(temp)
    end function  median
    
    !=======================================================================
    function sorted(array)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Return a sorted copy of the array.
    !
    !*ARGUMENTS:
    !    array : real
    !       Original array.
    !   n : integer
    !       Length of array.
    !   
    !-----------------------------------------------------------------------
        implicit  none
        real(r8), dimension(:), intent(in) :: array
        real(r8), dimension(:), allocatable :: sorted
        integer :: i,minIndex,lindx,uindx
        real(r8) :: temp
  
        ! make a copy of array
        lindx = lbound(array,1)
        uindx = ubound(array,1)
        allocate(sorted(lindx:uindx))
        sorted(:) = array(:)
        do i = lindx,uindx
            minindex = minloc(sorted(i:), 1) + i - 1
            if (sorted(i) > sorted(minindex)) then
                temp = sorted(i)
                sorted(i) = sorted(minindex)
                sorted(minindex) = temp
            end if
         end do
    end function  sorted
    
    !=======================================================================
    function pnorm(array,p)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Return the p-norm of a 1D array.
    !
    !*ARGUMENTS:
    !    array : real
    !       Original array.
    !   p : real
    !       Must be greater than 1.
    !       - 1 is taxicab norm (max absolute column sum)
    !       - 2 is euclidean norm
    !       - inf is maximum norm (max absolute row sum)
    !-----------------------------------------------------------------------
        implicit  none
        real(r8), dimension(:), intent(in) :: array
        real(r8), intent(in) :: p
        integer :: i,lindx,uindx
        real(r8) :: pnorm
        
        if(p<1) stop "ERROR: p-norm with p<1 is not a norm" 
  
        ! make a copy of array
        lindx = lbound(array,1)
        uindx = ubound(array,1)
        pnorm = 0
        do i = lindx,uindx
            pnorm = pnorm + abs(array(i))**p
        end do
        pnorm = pnorm**(1/p)
    end function  pnorm
    
    !=======================================================================
    function cpnorm(array,p)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Return the p-norm of a 1D array.
    !
    !*ARGUMENTS:
    !    array : complex
    !       Original array.
    !   p : real
    !       Must be greater than 1.
    !       - 1 is taxicab norm (max absolute column sum)
    !       - 2 is euclidean norm
    !       - inf is maximum norm (max absolute row sum)
    !-----------------------------------------------------------------------
        implicit  none
        complex(r8), dimension(:), intent(in) :: array
        real(r8), intent(in) :: p
        integer :: i,lindx,uindx
        real(r8) :: cpnorm
        
        if(p<1) stop "ERROR: p-norm with p<1 is not a norm"  
  
        ! make a copy of array
        lindx = lbound(array,1)
        uindx = ubound(array,1)
        cpnorm = 0
        do i = lindx,uindx
            cpnorm = cpnorm + abs(array(i))**p
        end do
        cpnorm = cpnorm**(1/p)
    end function  cpnorm
    
    !=======================================================================
    function ri(c)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !  Return the real and imaginary components of a complex value as a tuple.
    !
    !*ARGUMENTS:
    !    c : complex
    !       Complex value.
    !
    !-----------------------------------------------------------------------
        implicit  none
        complex(r8), intent(in) :: c

        real(r8), dimension(2) :: ri

        ri = (/ REAL(c), AIMAG(c) /)
    end function  ri

    !=======================================================================
    subroutine append_1d(list, element)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !  Append a 1D allocatable array.
    !
    !*ARGUMENTS:
    !    list : real, allocatable
    !       1D array to be appended.
    !    element: real
    !       Number appended to end of list.
    !
    !-----------------------------------------------------------------------
        implicit none
        ! declarations
        integer :: i, isize
        real(r8), intent(in) :: element
        real(r8), dimension(:), allocatable, intent(inout) :: list
        real(r8), dimension(:), allocatable :: tmp

        if(allocated(list)) then
            isize = size(list)
            allocate(tmp(isize+1))
            tmp(1:isize) = list(:)
            tmp(isize+1) = element
            deallocate(list)
            !call move_alloc(tmp, list) ! fortran 2003 standard
            allocate(list(isize+1))
            list(:) = tmp(:)
            deallocate(tmp)
        else
            allocate(list(1))
            list(1) = element
        end if

      end subroutine append_1d

    !=======================================================================
    subroutine append_2d(list, elements)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !  Append a 2D allocatable array.
    !
    !*ARGUMENTS:
    !    list : real, allocatable
    !       2D array to be appended.
    !    elements: real
    !       1D array appended to 1st dimension of list.
    !
    !-----------------------------------------------------------------------
        implicit none
        ! declarations

        integer :: isize,jsize
        real(r8), intent(in) :: elements(:)
        real(r8), dimension(:,:), allocatable, intent(inout) :: list
        real(r8), dimension(:,:), allocatable :: tmp

        if(allocated(list)) then
            isize = size(list,dim=1)
            jsize = size(list,dim=2)
            if (isize /= size(elements,dim=1)) stop "Cannot append arrays without matched dimensions"
            allocate(tmp(isize,jsize+1))
            tmp(:,1:jsize) = list(:,:)
            tmp(:,jsize+1) = elements(1:isize)
            deallocate(list)
            !call move_alloc(tmp, list) ! fortran 2003 standard
            allocate(list(isize,jsize+1))
            list(:,:) = tmp(:,:)
            deallocate(tmp)
        else
            isize = size(elements,dim=1)
            jsize = 1
            allocate(list(isize,jsize))
            list(:,1) = elements(1:isize)
        end if

      end subroutine append_2d

    !=======================================================================
    subroutine iscdftf(m,ms,func,fs,funcm)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Compute 1D truncated-discrete Fourier forward transform.
    !
    !*ARGUMENTS:
    !    m : integer
    !       Fourier modes
    !    ms: integer
    !       Number of modes
    !   func : real array (in)
    !       Function on a regular grid 0-1
    !   fs : integer
    !       Number of function points
    !   funcm : real array (out)
    !       Fourier components
    !
    !-----------------------------------------------------------------------
        implicit none
        ! declarations
        integer :: i,j,ms,fs
        integer, dimension(ms), intent(in) :: m
        complex(r8), parameter :: xj = (0,1)
        complex(r8), dimension(0:fs), intent(in) :: func
        complex(r8), dimension(ms), intent(out) :: funcm
    
        funcm=0
        do i=1,ms
           do j=0,fs-1
              funcm(i)=funcm(i)+(1.0/fs)*func(j)*exp(-(twopi*xj*m(i)*j)/fs)
           enddo
        enddo
        return
    end subroutine iscdftf

    !=======================================================================
    subroutine iscdftb(m,ms,func,fs,funcm)
    !----------------------------------------------------------------------- 
    !*DESCRIPTION: 
    !  Compute 1D truncated-discrete Fourier backward transform.
    !
    !*ARGUMENTS:
    !    m : integer
    !       Fourier modes
    !    ms: integer
    !       Number of modes
    !   func : real array (in)
    !       Function on a regular grid 0-1
    !   fs : integer
    !       Number of function points
    !   funcm : real array (out)
    !       Fourier components
    !
    !-----------------------------------------------------------------------
        implicit none
        integer :: i,j,ms,fs
        integer, dimension(ms), intent(in) :: m
        complex(r8), parameter :: xj = (0,1)
        complex(r8), dimension(ms), intent(in) :: funcm
        complex(r8), dimension(0:fs), intent(out) :: func
    
        func=0
        do i=0,fs-1
           do j=1,ms
              func(i)=func(i)+funcm(j)*exp((twopi*xj*m(j)*i)/fs)
           enddo
        enddo
        func(fs)=func(0)
        return
    end subroutine iscdftb

    !=======================================================================
    subroutine check(stat)
    !-----------------------------------------------------------------------
    !*DESCRIPTION:
    !  Check if a netcdf call was successful. If not, raise an error.
    !
    !*ARGUMENTS:
    !    stat : integer
    !       Status returned by a netcdf library function
    !
    !-----------------------------------------------------------------------
        integer, intent (in) :: stat
        !stop if it is an error.
        if(stat /= nf90_noerr) then
          print *, trim(nf90_strerror(stat))
          stop "ERROR: failed to write/read netcdf file"
        endif
        return
    end subroutine check

end module utilities
