c-----------------------------------------------------------------------
c     file vacuum_io.f.
c     input and output.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. inglo
c     2. cardmo
c     3. dskmd1
c     4. readahg
c     5. readvacin
c     6. readvacin5
c     7. adjustm
c     8. setahgdir
c-----------------------------------------------------------------------
c     subprogram 1. inglo.
c     read data from inadjv.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine inglo
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nmap1 = 0
      nmpdsk = 0
      ss = 0.
      if ( ladj .eq. 1 ) then
         open ( 50, file='inadjv', status='old', form='formatted' )
         return
      endif
      if ( ldcon .eq. 1 ) return
      if ( lgpec .eq. 1 ) return
      if ( lrgato .eq. 1 ) return
      if ( lzio .eq. 1 ) then
         call shellb
         call zop ( outmap1, mp1, nmap1, ndsk, ss, 100 )
         call zop ( iomode, mp0, nmpdsk, ndsk, ss, 100 )
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
 100  call errmes ( outmod, 'inglo' )
      end
c-----------------------------------------------------------------------
c     subprogram 2. cardmo.
c     read data from modivmc.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine cardmo
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      
      logical, save :: warned = .false.

      character(8) under
      data under / "--------" /
      namelist / modes  / mfel,m,mth,n,mdiv,lsymz,lfunin,xiin,
     .     leqarcw, lpest1, lnova, ladj, ldcon, lgato, lrgato, lspark,
     $     ismth, lzio, mp0,mp1
      namelist / debugs / checkd, checke, check1, check2, checks,
     $     wall, lkplt, verbose_timer_output
      namelist / vacdat / ishape,aw,bw,cw,dw,tw,nsing,epsq,noutv,delg,
     .     idgt, idot, delfac, idsk, cn0
      namelist / diagns / lkdis, ieig, iloop,
     $     nloop,nloopr,
     .     lpsub, nphil, nphse, mx, mz, nph, xofsl,
     $     aloop, bloop, dloop, rloop, ntloop, deloop,
     $     nxlpin,nzlpin,epslp,xlpmin,xlpmax,zlpmin,zlpmax,linterior
      namelist / shape  / ipshp, xpl, apl,bpl, dpl,  a, b, r,
     $     abulg, bbulg, tbulg, qain
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 8001 format ( a20 )
 8002 format ( a60 )
 600  format ( 2i5, 1pe12.5 )
 9000 format ( 1x, 20a4, / 1x, 2a10 / )
 9001 format ( 1x, " form of data input" )
 9002 format ( 1x, " card data input" )
 9003 format(1x," lmax,lmin=",2i4," m,mdiv=",2i4,3x," n=",e12.4,/)
 9100 format ( 20a4 )
c-----------------------------------------------------------------------
c     read and write input data.
c-----------------------------------------------------------------------
      rewind inmode
      read ( inmode, 9100 )   (ntitle(i),i=1,20)
      write ( outmod, 9000 )   ntitle, ( under,i=1,2 )
      write ( outmod, 9001 )
      write ( outmod,9002 )
      rewind inmode ! make robust to namelists with no title line
      read(inmode,modes)
      read(inmode,debugs)
      read(inmode,vacdat)
      read(inmode,shape)
      read(inmode,diagns)
      write(outmod,modes)
      write(outmod,debugs)
      write(outmod,vacdat)
      write(outmod,shape)
      write(outmod,diagns)
c-----------------------------------------------------------------------
c     override namelist r in favor of r=0 for DCON
c-----------------------------------------------------------------------
      if(r/=0.0)then
         if(.not. warned) write(*,'(1x,a,es9.2,a,i2)')
     $        "  > Vacuum code overriding r",r," from vac.in, to be",0
         r=0.0
         warned = .true.
      endif 
c-----------------------------------------------------------------------
c     subsidiary computations.
c-----------------------------------------------------------------------
      write ( outpest, 9003 )lmax(1),lmin(1),m,mdiv,n
      mp     = m + 1
      nosurf = ( mp - 1 ) * mdiv + 1
      mth    = nths0
      mth1   = mth + 1
      mth2   = mth1 + 1
      no2pi  = n / twopi
      no2pi2 = no2pi * no2pi
      dth    = twopi / mth
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 * r2
      lfour = 1
      lfele = 0
      if ( lgato .eq. 1 ) then
         lfour = 0
         lfele = 1
      endif
      if ( (a .ge. 10.0) .or. (lspark .ne. 0) ) then
         farwal = .true.
c      else
c         farwal = .false.
      endif
c-----------------------------------------------------------------------
c     special stuff for lspark .ne. 0.
c-----------------------------------------------------------------------
      if ( lspark .ne. 0 ) then
         open (iodsk,file='vdata',status='unknown',form='formatted' )
     $
         write ( iodsk,8001) mp1
         write ( iodsk,8002)
         write ( iodsk,8002)
         write ( iodsk, 600 ) lmin(1), lmax(1), n
         lnsav = lmin(1)
         lxsav = lmax(1)
         lmin(1) = - lmax(1)
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. dskmd1.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine dskmd1
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension vecin(ntsin), xigr_(ntsin), xigi_(ntsin)
      dimension zerov(nths), thgr(nths)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 701  format (/, 'mfel, rgato,ndum2, ngato, ga1, fa1, qa1 = ',/,
     $     i5,1pe13.5,2i5,1p3e13.5, / )
 601  format ( 4i5, e13.5 )
 602  format ( /,'mthin, lmin,lmax, nadj, ga1, fa1, qa1 = ',/,
     $     4i5, 1p3e13.5, / )
 605  format ( 10e13.5 )
 702  format (/, 'mthin, lmin,lmax, ndcon, ga1, fa1, qa1 = ',/,
     $     4i5, 1p3e13.5, / )
 5702 format (/, 'mthin, lmin,lmax, n, qa1 = ',/,
     $        3i5, 2e13.5, / )
 8011 format ( 5i4, 1p5e14.6 )
 8021 format ( 1p10e14.6 )
 9600 format ( /, 1x, " mp1, qa1, fa1, ga1 = ", a,1p3e12.5,/ )
 8031 format ( 1p3e14.6 )
 9000 format (1x,20a4,/,
     $     1x, " equilibrium from disk, calculated on date",a)
 9100 format(20a4,a10)
 9200 format(10i5)
 9300 format(4e20.13)
 111  format ( /,"<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>",/,
     $     1x, "ipshp = ", i3, " qa1 = ", e13.5,/,
     $     "<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>",/ )
 9500 format ( /, 1x, "Mapping parameters, nsf0, ntsin0 = ", 2i5,/,
     $     1x, "Working parameter, nths0  = ", i5,/,
     $     1x, "nfm, mtot = ", 2i5,/,
     $     1x, "r, upsiln, mthin, nosurf = ", 1p2e12.5, 2i5,/ )
c-----------------------------------------------------------------------
c     defaults.
c-----------------------------------------------------------------------
      lcdf = 0
      do i = 1, mth2
         delta(i) = 0.0
         xjacob(i) = 0.0
      enddo
      r = sqrt(r2)
c-----------------------------------------------------------------------
c     gato inputs.
c-----------------------------------------------------------------------
      if ( lrgato .eq. 1 ) then
         lzio = 0
         dx0 = 0.5
         call readvacin ( mfel,rgato,ndum2,ngato,qa1,xinf,zinf,
     $        delta, vecin, xigr_,xigi_, mth,mth1,mth2, ndfel,dx0,
     $        ieig, outmod, iotty )
         call wrtout ( mfel, xigr_, "xigr", 1, mfel )
         l11 = lmin(1)
         l22 = lmax(1)
         call fanal ( xigr_, mfel, xirc, xirs, l11,l22, pye,-0.5_8 )
         call fanal ( xigi_, mfel, xiic, xiis, l11,l22, pye,-0.5_8 )
         llnn = l22 - l11 + 1
         call vecwrt ( llnn, xirc, "xirc(l)", 1, llnn, outmod,0 )
         call vecwrt ( llnn, xirs, "xirs(l)", 1, llnn, outmod,0 )
         call vecwrt ( llnn, xiic, "xiic(l)", 1, llnn, outmod,0 )
         call vecwrt ( llnn, xiis, "xiis(l)", 1, llnn, outmod,0 )
         dth = twopi / mth
         n = ngato
         do ii = 1, mth+1
            delta(ii) = 0.0
         enddo
         write ( outmod, '(/,"******* DELTA set to 0.0 ******" )' )
         write ( iotty, 701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
         write ( outmod,701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
      endif
c-----------------------------------------------------------------------
c     more gato computations.
c-----------------------------------------------------------------------
      if ( lgato .ne. 0 ) then
         if ( lgato .eq. 1 ) then
            lfele = 1
            lmin(1) = 1
            lmax(1) = mfel
         endif
         if (  lrgato .eq. 0 )
     $        call adjustm ( mth, mfel, mth1,mth2, ndfel, iotty,outmod )
         dth = twopi/mth
         write ( outmod, '(/,5x, "mth, mfel, ndfel, dth = ",
     $        3i5, 1pe13.5 )' ) mth, mfel, ndfel, dth
         write ( iotty,  '(/,5x, "mth, mfel, ndfel, dth = ",
     $        3i5, 1pe13.5 )' ) mth, mfel, ndfel, dth
      endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
      if ( ladj .eq. 1 ) then
         read ( 50, 601 ) mthin1,lmin(1),lmax(1),nadj, qa1
         mthin = mthin1 - 1
         mthin2 = mthin1 + 1
         n = nadj
         write ( iotty, 602 ) mthin,lmin(1),lmax(1),nadj, ga1,fa1,qa1
         write ( outmod,602 ) mthin,lmin(1),lmax(1),nadj, ga1,fa1,qa1
         read ( 50,605 ) ( vecin(i), i=1,mthin1 )
         call trans ( vecin,mthin, xinf,mth )
         read ( 50,605 ) ( vecin(i), i = 1,mthin1 )
         call trans ( vecin,mthin, zinf,mth )
         read ( 50,605 ) ( vecin(i), i = 1,mthin1 )
         call trans ( vecin,mthin, delta,mth )
         close (50)
         go to 1111
      endif
c-----------------------------------------------------------------------
c     dcon inputs.
c-----------------------------------------------------------------------
      if ( ldcon .eq. 1 ) then
         lzio = 1
         call readahg (ahgfile, ahgdir, mthin,lmin(1),lmax(1),ndcon,qa1,
     $      xinf, zinf, delta, vecin, mth )
         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ndcon
         write ( iotty, 702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
         write ( outmod,702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
         go to 1111
      endif
c-----------------------------------------------------------------------
c     gpec inputs.
c-----------------------------------------------------------------------
      if ( lgpec .eq. 1 ) then
         lzio = 1
         ieig = 0
c
         dx0 = 0.0
         dx1 = 0.0

         call readvacin5 ( nxlpin,nzlpin, xlpmin,xlpmax, zlpmin,zlpmax,
     $        mthin, lmin(1),lmax(1),ntor, qa1,
     $        xinf,zinf,delta, vecin, bnlr,bnli,
     $        mth,dx0,dx1, ieig)

         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ntor
         write ( iotty, 5702 ) mthin,lmin(1),lmax(1),n, qa1
         write ( outmod,5702 ) mthin,lmin(1),lmax(1),n, qa1
         go to 1111
      endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
      if ( lzio .eq. 1 ) then
         r=0
         lgivup=1
         nadres = 1
         call zrd(iomode,real(ntitle(1),r8),43,nadres,lgivup,999)
         write ( outmod, 9000 )   ntitle,dat
         lj = 0
         zma = 0.0
         write ( iodsk, 8011 ) nosurf,mthin, lj,mj,nj, xzero, r,
     $        upsiln, xma, zma
         r2     = r * r
         r4     = r2 * r2
         r6     = r4 * r2
         mthin1 = mthin + 1
         mthin2 = mthin + 2
      endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
      if ( lzio .eq. 1 ) then
         nadres = 50
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
         call trans ( vecin,mthin, xinf,mth )
         nadres = nadres + mthin2
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
         call trans ( vecin,mthin, zinf,mth )
         length = ntsin * nsf
         ladres = 50 + 2*mthin2 + nosurf
         nadres = ladres + (nosurf-1)*ntsin
         lgivup = 1
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
         call trans ( vecin,mthin, grpssq,mth )
         if ( .not. lpest1 ) then
            nadres = nadres + 8*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
            call trans ( vecin,mthin, xjacob,mth )
            length = ntsin * nsf
            ladres = 50 + 2*mthin2 + nosurf
            nadres = ladres + (nosurf-1)*ntsin + 10*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
            call trans ( vecin,mthin, delta,mth )
         endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
         if ( lpest1 ) then
            do i = 1, mth2
               xjacob(i) = upsiln * xinf(i)**2 / ( twopi*r )
            enddo
         endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
         nadres = 50
         nadres = nadres + nosurf*3 - 1
         call zrd(outmap1,qa1,1,nadres,lgivup,999)
         nadres = nadres + nosurf*2
         call zrd(outmap1,ga1,1,nadres,lgivup,999)
         nadres = nadres + nosurf*2
         call zrd(outmap1,fa1,1,nadres,lgivup,999)
         write ( outmod,9600 ) mp1, qa1, fa1, ga1
         write ( iotty, 9600 ) mp1, qa1, fa1, ga1
         write ( iodsk, 8031 ) qa1, ga1, fa1
c-----------------------------------------------------------------------
c     nova inputs.
c-----------------------------------------------------------------------
         if ( lnova ) then
            length = ntsin * nsf
            ladres = 50 + 2*mthin2 + nosurf
            nadres = ladres + (nosurf-1)*ntsin
            lgivup = 1
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, grpssq,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xsq,mth )
            nadres=nadres+5*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, gpsdth,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xsqdth,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xjacob,mth )
            mthd2p1=mth/2+1
            do i=1,mthd2p1
               xsdtxs(i)=xsqdth(i)/xsq(i)
               gpdtgp(i)=gpsdth(i)/grpssq(i)
               xjdtxj(i)=0.5*(mj*xsdtxs(i)-nj*gpdtgp(i))
               xjacob(mth2-i)=xjacob(i)
               xjdtxj(mth2-i)=-xjdtxj(i)
               delta(mth2-i)=-delta(i)
            enddo
            xjdtxj(1)=0.
            xjdtxj(mthd2p1)=0.
            delta(1)=0.
            delta(mthd2p1)=0.
         endif
         call zcl ( outmap1, 999 )
         call zcl ( iomode, 999 )
      endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
 1111 continue
      if ( ipshp .eq. 1 ) then
         qa1  = qain
         rgato = xpl
         write ( outmod, 111 ) ipshp, qa1
         write ( iotty,  111 ) ipshp, qa1
         do i = 1, mth2
            theta = (i-1) * dth
            xinf(i) = xpl + apl * cos(theta+dpl*sin(theta))
            zinf(i) =     - bpl* apl * sin(theta)
            delta(i) = 0.0
            xjacob(i) = 0.0
         enddo
      endif
      write ( outmod,9500 ) nsf0, ntsin0, nths0, nfm, mtot,
     $     r, upsiln, mthin, nosurf
      write ( iotty,9500 ) nsf0, ntsin0, nths0, nfm, mtot,
     $     r, upsiln, mthin, nosurf
      do i = 1, mth1
         zerov(i) = 0.0
         thgr(i) = (i-1)*dth
      enddo

      call arrays
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
 999  call errmes(outpest,'dskmd1')
      end
c-----------------------------------------------------------------------
c     subprogram 4. readahg.
c     read data from dcon.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine readahg (ahgfile, ahgdir, mthin,lmin,lmax,ndcon,qa1,
     $     xinf, zinf, delta, vecin, mth )
      use local_mod, only: r8
      use vglobal_mod, only: dcon_set, mth_dcon, lmin_dcon,lmax_dcon,
     $     nn_dcon, qa1_dcon, x_dcon, z_dcon, delta_dcon
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character(128) ahgdir
      character(128) ahgfile
      integer mthin,lmin,lmax,ndcon,ith
      dimension xinf(*), zinf(*), delta(*), vecin(*)
c-----------------------------------------------------------------------
c     don't bother with file io is passed in memory.
c-----------------------------------------------------------------------
      if(dcon_set)then
         mthin = mth_dcon
         lmin = lmin_dcon
         lmax = lmax_dcon
         ndcon = nn_dcon
         qa1 = qa1_dcon
         mthin1 = mthin + 1
         ! note mthin is 1 less than len(x_dcon). bug? or vac doesn't close theta?
         call trans ( x_dcon,mthin, xinf,mth )
         call trans ( z_dcon,mthin, zinf,mth )
         call trans ( delta_dcon,mthin, delta,mth )
         vecin(1:mthin) = delta_dcon(1:mthin)
      else
c-----------------------------------------------------------------------
c     read data.
c-----------------------------------------------------------------------
         print *, "Using ahgfile: ", ahgfile
         open(unit=3,file=trim(ahgdir)//'/'//trim(ahgfile))
         read(3,*)mthin
         read(3,*)lmin
         read(3,*)lmax
         read(3,*)ndcon
         read(3,*)qa1
         mthin1 = mthin + 1
         read(3,'(//)')
         read(3,*)(vecin(ith),ith=1,mthin1)
         read(3,'(//)')
         read(3,*)(vecin(ith),ith=1,mthin1)
         read(3,'(//)')
         read(3,*)(vecin(ith),ith=1,mthin1)
         call trans ( vecin,mthin, xinf,mth )
         read(3,'(//)')
         read(3,*)(vecin(ith),ith=1,mthin1)
         call trans ( vecin,mthin, zinf,mth )
         read(3,'(//)')
         read(3,*)(vecin(ith),ith=1,mthin1)
         call trans ( vecin,mthin, delta,mth )
         close(unit=3)
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 5. readvacin.
c     read gato input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine readvacin ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $     delta, vecin,xigr,xigi, mth,mth1,mth2, ndfel, dx0,
     $     ireig, nout1, nout2 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer mthin,ndum2,ngato,ith
      dimension xinf(*), zinf(*), delta(*), vecin(*), xigr(*),xigi(*)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      open(unit=3,file='vacin')
      write ( nout1, '(/"reading VACIN",/)' )
      write ( nout2, '(/"reading VACIN",/)' )
      read(3,'(//)')
      read(3,*)mthin
      read(3,*)rgato
      read(3,*)ndum2
      read(3,*)ngato
      read(3,*)qa1
      mthin1 = mthin + 1
      call adjustm ( mth, mthin, mth1,mth2,  ndfel, nout1,nout2 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xinf,mth, dx0 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zinf,mth, dx0 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, delta,mth, dx0 )
      if ( ireig .eq. 8 ) then
         read(3,'(//)')
         read(3,*)(xigr(ith),ith=1,mthin1)
         read(3,'(//)')
         read(3,*)(xigi(ith),ith=1,mthin1)
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      close(unit=3)
 2    return
      end
c-----------------------------------------------------------------------
c     subprogram 6. readvacin5.
c     read vacin5 input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine readvacin5 ( nlx,nlz, xll,xlr, zlb,zlt,
     $     mthin, lmin,lmax,ntor, qa1,
     $     xinf,zinf, delta, vecin, bnlr,bnli,
     $     mth, dx0,dx1, ireig5)
      use local_mod, only: r8
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      INTEGER mthin, lmin, lmax, ntor, ith,jl, jmax1
      DIMENSION xinf(*),zinf(*),delta(*), vecin(*),bnlr(*), bnli(*)
c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN ( UNIT=7, FILE='vacin5' )
      READ(7,'(/)')
      READ(7,*) nlx
      READ(7,*) nlz
      READ(7,*) xll
      READ(7,*) xlr
      READ(7,*) zlb
      READ(7,*) zlt
      READ(7,*) mthin
      READ(7,*) lmin
      READ(7,*) lmax
      READ(7,*) ntor
      READ(7,*) qa1

      mthin1 = mthin + 1
      jmax1 = lmax - lmin + 1
c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------

      READ(7,'(/)')
      READ(7,*)(vecin(ith),ith=1,mthin1)
      READ(7,'(/)')
      READ(7,*)(vecin(ith),ith=1,mthin1)
      READ(7,'(/)')
      READ(7,*)(vecin(ith),ith=1,mthin1)
      CALL transdxx ( vecin,mthin, xinf,mth, dx0,dx1)
      READ(7,'(/)')
      READ(7,*)(vecin(ith),ith=1,mthin1)
      CALL transdxx ( vecin,mthin, zinf,mth, dx0,dx1)

      !Xwall
      ! READ(7,*)(vecin(ith),ith=1,mthin1)
      ! CALL transdxx ( vecin,mthin, xwin,mth, dx0,dx1 )

      !Zwall
      ! READ(7,*)(vecin(ith),ith=1,mthin1)
      ! CALL transdxx ( vecin,mthin, zwin,mth, dx0,dx1 )

      READ(7,'(/)')
      READ(7,*)(vecin(ith),ith=1,mthin1)
      CALL transdxx ( vecin,mthin, delta,mth, dx0,dx1)
      READ(7,'(/)')
      READ(7,*)(bnlr(jl),jl=1,jmax1)
      READ(7,'(/)')
      READ(7,*)(bnli(jl),jl=1,jmax1)
      CLOSE ( UNIT=7 )
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 7. adjustm.
c     read netcdf data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine adjustm ( mth, mfel, mth1,mth2, ndfel, nout1,nout2 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if ( mth .eq. mfel ) then
         ndfel = 1
         return
      endif
      mth00 = mth
      ndfel = mth / mfel
      if ( ndfel .lt. 1 ) ndfel = 1
      ndfel = ( (ndfel+1)/2 ) * 2
      mth = ndfel * mfel
      if ( mth00 .ne. mth ) then
         mth1 = mth + 1
         mth2 = mth1 + 1
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
c     subprogram 7. setahgdir.
c     set path to dcon outputs in case they aren't in the current dir.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine setahgdir ( directory )
      use vglobal_mod
      character(len=*), intent(in) :: directory
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      ahgdir = directory
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
c     subprogram 7. set_dcon_params.
c     In memory io for DCON
c     Speedup over original ascii io (readahg) for multiple calls
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine set_dcon_params ( mthin,lmin,lmax,nnin,qa1in,xin,
     $     zin, deltain )
      use local_mod, only: r8
      use vglobal_mod, only: dcon_set, mth_dcon, lmin_dcon,lmax_dcon,
     $     nn_dcon, qa1_dcon, x_dcon, z_dcon, delta_dcon
      integer :: mthin,lmin,lmax,nnin,mthin1,mthin5,i
      real(r8) :: qa1in
      real(r8), dimension(*) :: xin, zin, deltain
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mthin1 = mthin + 1
      mthin5 = mthin1 + 4
      if(dcon_set) call unset_dcon_params
      allocate(x_dcon(mthin5), z_dcon(mthin5), delta_dcon(mthin5))
      mth_dcon = mthin
      lmin_dcon = lmin
      lmax_dcon = lmax
      nn_dcon = nnin
      qa1_dcon = qa1in
      x_dcon(1:mthin1) = xin(1:mthin1)
      z_dcon(1:mthin1) = zin(1:mthin1)
      delta_dcon(1:mthin1) = deltain(1:mthin1)

      do i=1,5
         x_dcon(mthin+i) = x_dcon(i)
         z_dcon(mthin+i) = z_dcon(i)
         delta_dcon(mthin+i) = delta_dcon(i)
      enddo

      dcon_set = .TRUE.
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
c     subprogram 8. unset_dcon_params.
c     Tell vacuum to ignore in-memory dcon params and read ascii values
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine unset_dcon_params
      use vglobal_mod, only: dcon_set, x_dcon, z_dcon, delta_dcon
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      deallocate(x_dcon, z_dcon, delta_dcon)
      dcon_set = .FALSE.
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
