c-----------------------------------------------------------------------
c     file vacuum_ma.f.
c     Main routines for Morrell Chance's vacuum eigenvalue code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. mscvac
c     2. defglo
c     3. ent33
c     4. funint
c-----------------------------------------------------------------------
c     subprogram 1. mscvac.
c     Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mscvac(wv,mpert,mtheta,mthvac,complex_flag)
      USE vglobal_mod
      implicit real*8 (a-h,o-z)

      integer mpert,mtheta,mthvac
      complex*16 wv(mpert,mpert)
      logical, intent(in) :: complex_flag

      complex(8), parameter :: ifac=(0,1)
      dimension xi(nfm), xii(nfm), lfm(nfm), xilnq(nfm), xiilnq(nfm)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   format( //,4x, "omsq=",e12.5,/,
     $     4x, "i",3x,"l",3x,"xi-edge",5x,"xi*(l-nq)",
     $     /, 101(1x,2i4,1p2e12.4,/),/ )
 20   format( i5 )
 30   format( e12.5,/,(8e12.5) )
 40   format( //,4x, "omsq=",e12.5,/,
     $     4x, "i",3x,"l",3x,"xi-edge",5x,"xi*(l-nq)",2x,
     $     "xii-edge",5x,"xii*(l-nq)",
     $     /, 101(1x,2i4,1p4e12.4,/),/ )
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ntsin0=mtheta+1
      nths0=mthvac
      nfm=mpert
      mtot=mpert
      call global_alloc(nths0,nfm,mtot,ntsin0)
c-----------------------------------------------------------------------
c     initialization.
c-----------------------------------------------------------------------
      call defglo
      open (iotty,file='mscvac.out',status='unknown')
      open (outpest,file='pestotv',status='unknown',form='formatted')
      open (inmode,file='vac.in',status='old', form='formatted' )
      open (outmod,file='modovmc',status='unknown', form='formatted' )
      call msctimer ( outmod, "top of main" )
      call ent33
      If ( lspark .ne. 0 ) call testvec
      if ( ieig .eq. 0 ) goto 99
      jmax1 = lmax(1) - lmin(1) + 1
      do j1 = 1, jmax1
         l1 = j1 + lmin(1) - 1
         lfm(j1) = l1
         xi(j1)  = 0.0
         xii(j1) = 0.0
      enddo
c-----------------------------------------------------------------------
c     main conditional.
c-----------------------------------------------------------------------
      if(ieig .eq. 1)then
         do  l = 1, jmax1
            xilnq(l) = ( lfm(l)-n*qa1 ) * xi(l)
         enddo
         write(outmod,10)omsq,( l,lfm(l),xi(l),xilnq(l),l=1,jmax1)
         goto 99
      elseif(ieig .eq. 4)then
         open ( 60, file='outidst', status='old', form='formatted' )
         mflag=1
         do while(mflag .ne. 0)
            read ( 60,20 ) mflag
c            if ( mflag .eq. 0 )exit
         enddo
         read ( 60,30 ) omsq, ( xi(l),l = 1,jmax1 )
      elseif(ieig .eq. 5)then
         do j1 = 1, jmax1
            xi(j1) = xiin(j1)
         enddo
      elseif(ieig .eq. 8)then
         do j1 = 1, jmax1
            xii(j1) = xirc(j1) - xiis(j1)
            xi(j1)  = xiic(j1) + xirs(j1)
         enddo
      endif
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      do l = 1, jmax1
         xilnq(l)  = ( lfm(l)-n*qa1 ) * xi(l)
         xiilnq(l) = ( lfm(l)-n*qa1 ) * xii(l)
      enddo
      write ( outmod,40 ) omsq, ( l,lfm(l),xi(l),xilnq(l),
     $     xii(l),xiilnq(l), l = 1,jmax1 )
c-----------------------------------------------------------------------
c     copy vacuum response matrix to output.
c-----------------------------------------------------------------------
 99   continue
      IF(complex_flag)THEN
         wv=vacmat+ifac*vacmtiu
      ELSE
         wv=vacmat
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      call msctimer ( outmod, "end of main" )
      close(iotty)
      call global_dealloc
      call cleanup
      return
      end
c-----------------------------------------------------------------------
c     subprogram 2. defglo.
c     defines logical unit numbers and control switches.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine defglo
      USE vglobal_mod
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c     define constants.
c-----------------------------------------------------------------------
      zero   = 0.0e0
      pt1    = 1.e-1
      half   = 0.5e0
      one    = 1.0e0
      two    = 2.0e0
      three  = 3.0e0
      four   = 4.0e0
      five   = 5.0e0
      seven  = 7.0e0
      epsq = 1.0e-5
      pye    = 3.1415926535897931_8
      twopi  = two * pye
      twopi2 = twopi * twopi
      alx     = 1
      alz     = 1
      n       = 1
      m       = 2
      mp     = 3
      minc   = 10
      mth    = 128
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthin  = 128
      mthin1 = mthin + 1
      mthin2 = mthin + 2
      mdiv   = 2
      idgt = 0
      nosurf = ( mp - 1 ) * mdiv + 1
      n1surf = 4
      npsurf = 6
      dpsi   = one / mdiv / m
      bit    = 1.0e-8
      amu0   = four * pye * 1.e-7
      dth    = twopi / mth
      dthinc = dth  / minc
      gamma  = five / three
      p0     = zero
      upsiln = one
      r      = 20.0e0
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 *r2
      rgato = r
      fa1 = 1.0
      ga1 = 1.0
      qa1 = 1.0
      isymz  = 2
      xiin(1) = 1.0
      xma = 1.0
      zma = 0.0
      xzero = 1.0
      ntloop = 8
      deloop = .001
c-----------------------------------------------------------------------
c     define logical unit numbers.
c-----------------------------------------------------------------------
      idsk = 1
      intty = 5
      iotty = 86
      inpest = 2
      outpest = 3
      iomode = 19
      inmode = 22
      outmod = 23
      iodsk = 16
      outmap1 = 18
      iovac = 36
      do i = 1, 3
         nout(i) = 0
         nout0(1) = 0
      enddo
      nout(1) = 6
      nout(2) = outmod
c      mp0 = 'mapdsk'
c      mp1 = 'mpout1'
c-----------------------------------------------------------------------
c     define logicals.
c-----------------------------------------------------------------------
      lzio = 1
      lsymz  = .false.
      check1 = .false.
      check2 = .false.
      lanal   = .false.
      lkdis = .true.
      lpest1 = .false.
      lnova = .false.
      lspark= 0
      ladj = 0
      ldcon = 0
      wall   = .false.
      lkplt = 0
      do ich = 1, 60
         seps(ich:ich) = '.'
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. ent33.
c     entry point for segment 33.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine ent33
      USE vglobal_mod
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      call cardmo
      call inglo
      call dskmd1
      call funint
      if(.not. wall) call vaccal
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. funint.
c     function initialization.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine funint
      USE vglobal_mod
      implicit real*8 (a-h,o-z)

      dimension zork1(nths), zork2(nths), dlenth(nths)
      dimension the(nths)
c-----------------------------------------------------------------------
c     functions not fixed by the specific equilibrium.
c-----------------------------------------------------------------------
      upsil2 = upsiln**2
      if ( lpest1 ) qa1 = upsiln*ga1/(twopi*fa1)
      if ( ipshp .eq. 1 ) qa1 = qain
      f02 = fa1**2
c-----------------------------------------------------------------------
c     calculate arc length on surface.
c-----------------------------------------------------------------------
      mth1 = mth + 1
      mth2 = mth + 2
      mth3 = mth + 3
      mth4 = mth + 4
      mth5 = mth + 5
      do i = 1, mth1
         the(i) = (i-1) * dth
      enddo
c-----------------------------------------------------------------------
c     get derivative on plasma points.
c-----------------------------------------------------------------------
      call difspl ( mth, the, xpla, xplap )
      call difspl ( mth, the, zpla, zplap )
      do i = 1, mth1
         dlenth(i) = sqrt ( xplap(i)**2 + zplap(i)**2 )
      enddo
      do i = 1, mth1
         zork1(i+2) = dlenth(i)
      enddo
      zork1(1) = dlenth(mth-1)
      zork1(2) = dlenth(mth) 
      zork1(mth3) = dlenth(1)
      zork1(mth4) = dlenth(2)
      zork1(mth5) = dlenth(3)
      call indef4 ( zork1, zork2, dth, 3,mth3, alen, 0 )
      do i = 1, mth1
         slngth(i) = zork2(i+2)
      enddo
      write( iotty,'("Circumference of Plasma = ",1pe11.3)')
     $     slngth(mth1)
      write(outmod,'("Circumference of Plasma = ",1pe11.3)')
     $     slngth(mth1)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
