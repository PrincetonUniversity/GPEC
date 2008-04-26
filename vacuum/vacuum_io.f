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
c     6. mdskrd0
c     7. mdskrd1
c     8. adjustm
c-----------------------------------------------------------------------
c     subprogram 1. inglo.
c     read data from inadjv.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine inglo
      USE vglobal_mod
      implicit real*8 (a-h,o-z)
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
      implicit real*8 (a-h,o-z)

      real*8 under
      data under / "--------" /
      namelist / modes  / mfel,m,mth,n,mdiv,lsymz,lfunin,xiin,
     .     leqarcw, lpest1, lnova, ladj, ldcon, lgato, lrgato, lspark, 
     $     ismth, lzio, mp0,mp1
      namelist / debugs / checkd, checke, check1, check2, checks,
     $     wall, lkplt
      namelist / vacdat / ishape,aw,bw,cw,dw,tw,nsing,epsq,noutv,delg,
     .     idgt, idot, delfac, idsk, cn0
      namelist / diagns / lkdis, ieig, iloop, xloop, zloop,
     $     nloop,nloopr, 
     .     lpsub, nphil, nphse, mx, mz, nph, xofsl,
     $     aloop, bloop, dloop, rloop, ntloop, deloop
      namelist / shape  / ipshp, xpl, apl,bpl, dpl,  a, b, r,
     $     abulg, bbulg, tbulg, qain
c$$$      namelist / sprk / nminus, nplus, mphi, lwrt11,civ,
c$$$     $     sp2sgn1, sp2sgn2, sp2sgn3, sp2sgn4, sp2sgn5,
c$$$     $     sp3sgn1, sp3sgn2, sp3sgn3, sp3sgn4, sp3sgn5,
c$$$     $     lff, ff, fv
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
      rsave  = r
      write ( outmod,9002 )
      read(inmode,modes)
      read(inmode,debugs)
      read(inmode,vacdat)
      read(inmode,shape)
      read(inmode,diagns)
c      read(inmode,sprk)
      write(outmod,modes)
      write(outmod,debugs)
      write(outmod,vacdat)
      write(outmod,shape)
      write(outmod,diagns)
c      write(outmod,sprk)
c-----------------------------------------------------------------------
c     subsidiary computations.
c-----------------------------------------------------------------------
      r      = rsave
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
      implicit real*8 (a-h,o-z)

      dimension vecin(ntsin), xigr(ntsin), xigi(ntsin)
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
c     zero arrays.
c-----------------------------------------------------------------------
      lcdf = 0
      do i = 1, mth2
         delta(i) = 0.0
         xjacob(i) = 0.0
      enddo
c-----------------------------------------------------------------------
c     gato inputs.
c-----------------------------------------------------------------------
      if ( lrgato .eq. 1 ) then
         lzio = 0
         dx0 = 0.5
         call readvacin ( mfel,rgato,ndum2,ngato,qa1,xinf,zinf,
     $        delta, vecin, xigr,xigi, mth,mth1,mth2, ndfel,dx0,
     $        ieig, outmod, iotty )
         call wrtout ( mfel, xigr, "xigr", 1, mfel )
         l11 = lmin(1)
         l22 = lmax(1)
         call fanal ( xigr, mfel, xirc, xirs, l11,l22, pye,-0.5_8 )
         call fanal ( xigi, mfel, xiic, xiis, l11,l22, pye,-0.5_8 )
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
         call readahg ( mthin,lmin(1),lmax(1),ndcon,qa1,xinf,zinf,
     $        delta, vecin, mth )
         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ndcon
         write ( iotty, 702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
         write ( outmod,702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
         go to 1111
      endif
c-----------------------------------------------------------------------
c     more computations.
c-----------------------------------------------------------------------
      if ( lzio .eq. 1 ) then
         r=0
         lgivup=1
         nadres = 1
         call zrd(iomode,ntitle(1),43,nadres,lgivup,999)
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
      else
c         if ( lcdf .eq. 1 ) call mdskrd0
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
      else
c         if ( lcdf .eq. 1 ) call mdskrd1
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
      subroutine readahg ( mthin,lmin,lmax,ndcon,qa1,xinf,zinf,
     $     delta, vecin, mth )
      implicit real*8 (a-h,o-z)
      integer mthin,lmin,lmax,ndcon,ith
      dimension xinf(*), zinf(*), delta(*), vecin(*)
c-----------------------------------------------------------------------
c     read data.
c-----------------------------------------------------------------------
      open(unit=3,file='ahg2msc.out')
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
      implicit real*8 (a-h,o-z)
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
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 6. mdskrd0.
c$$$c     read netcdf data.
c$$$c-----------------------------------------------------------------------
c$$$c-----------------------------------------------------------------------
c$$$c     declarations.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine mdskrd0
c$$$      USE vglobal_mod
c$$$      implicit real*8 (a-h,o-z)
c$$$
c$$$      include 'netcdf.inc'
c$$$      character*72 datype(6)
c$$$c-----------------------------------------------------------------------
c$$$c     format statements.
c$$$c-----------------------------------------------------------------------
c$$$ 1001 format( /,3x, "Id #", 2x, "Variable", 3x, "V-Type", 1x,
c$$$     $     "Attributes" 1x, "Dimensions" )
c$$$ 1002 format( 2x, i5, 2x, a10, i5,4x, i5,3x, 5i5 ) 
c$$$ 9000 format(/, 1x, a,/, 1x, " jobid = ", a )
c$$$ 2002 format(/,'           VACUUM on ', A10, ' at ', a10,/,
c$$$     $     '          MAPPING on ', a10, ' at ', a10,/,
c$$$     $     ' From EQUILIBRIUM on ', a10, ' at ', a10,
c$$$     $     ' to dsk: ', a ,/,
c$$$     $     ' with NOSURF = ',i3,' MTHIN = ',i4,/,
c$$$     $     '      MJACX, NJACG, LJACB =',3i2,/,
c$$$     $     ' Jacobian =  X ** ',i1,
c$$$     $     '/ ( Grad-PSI ** ',i1,' * B ** ',i1,')',/ )
c$$$ 3066 format(/," betap,betat,beta=",f6.2,f8.2,"%   ",f8.2,"% ")
c$$$ 3067 format(" betap,betat,beta(old def.)=",f6.2,f8.2,"%   ",f8.2,"% ")
c$$$ 3068 format(" beta star =",f8.2,"%")
c$$$ 3777 format(" tor field=  ",f6.2,"(T)     IP=",f7.3,"(MA) ",/,
c$$$     1     " I(MA)/A(m)B(T)= ",f6.2,"   troyon factor",f6.2,/,
c$$$     2     " q(axis)  =  ",f6.2,"        q(edge)",f6.2)
c$$$ 3778 format(" qstar    =  ",f6.2,"  qstar/q(1) = ",f6.2)
c$$$ 3779 format(" bt2dv= ",e12.4," pdv= ",e12.4," dv= ",e12.4)
c$$$ 3780 format(" li(GA)= ",e12.4," dlp= ",e12.4)
c$$$ 9600 format( /, 1x, " qa1, fa1, ga1 = ", 1p3e12.5,/ )
c$$$ 9500 format( /, 1x, "r, upsiln, xma  = ", 1p3e12.5,/ )
c$$$c-----------------------------------------------------------------------
c$$$c     define names for data types.
c$$$c-----------------------------------------------------------------------
c$$$      datype(1) = "byte:      Eight-bit data. For saving space."
c$$$      datype(2) = "character: Synonymous with byte. ASCII characters."
c$$$      datype(3) = "short:     16-bit integers."
c$$$      datype(4) = "long:      32-bit integers."
c$$$      datype(5) = "float:     32-bit IEEE floating-point."
c$$$      datype(6) = "double:    64-bit IEEE floating-point."
c$$$c-----------------------------------------------------------------------
c$$$c     computations.
c$$$c-----------------------------------------------------------------------
c$$$      call writg1 ( "Opening a NetCDF File (5.3).", seps, 3, nout ) 
c$$$      cdfid = ncopn ( cdfin, ncnowrit, rcode )
c$$$      call nwopn ( cdfid, cdfin, "ncnowrit", rcode, nout )
c$$$      call writg1 ( "Inquiring about a NetCDF File (5.7).",
c$$$     $     seps, 3, nout )
c$$$      call ncinq ( cdfid, ndims, nvars, ngatts, recid, rcode )
c$$$      call nwinq ( cdfid, ndims, nvars, ngatts, recid, rcode, nout )
c$$$      call writg1 ( "Dimensions' Information (6.3).", seps, 1, nout )
c$$$      do idims = 1, ndims
c$$$         call ncdinq ( cdfid, idims, vname, dimsiz(idims), rcode )
c$$$         call nwdinq ( cdfid, idims, vname, dimsiz(idims), rcode,
c$$$     $        nout )
c$$$      enddo
c$$$      call writg1 ( "The Record Dimension info., if any.. ",
c$$$     $     seps, 3, nout )
c$$$      if ( recid .ne. -1 ) then
c$$$         call ncdinq ( cdfid, recid, recnam, nrecs, rcode )
c$$$         call nwdinq ( cdfid, recid, recnam, nrecs, rcode, nout )
c$$$      endif
c$$$      if ( check1 ) then
c$$$         call writg1 ( "The Global Atrtributes (8.3): ", seps, 3, nout )
c$$$         do igatts = 1, ngatts
c$$$            call ncanam ( cdfid, ncglobal, igatts, attnam, rcode )
c$$$            call ncainq ( cdfid, ncglobal, attnam, attype, attlen,
c$$$     $           rcode )
c$$$            if ( attype .eq. 2 ) then
c$$$               call ncagtc ( cdfid, ncglobal, attnam, astrng, maxa1,
c$$$     $              rcode )
c$$$               call nwagtc ( cdfid, ncglobal, attnam, astrng, attlen,
c$$$     $              maxa1, rcode, nout )
c$$$            else
c$$$               call ncagt ( cdfid, ncglobal, attnam, attval, rcode )
c$$$               call nwagt ( cdfid, ncglobal, attnam, attval, attlen,
c$$$     $              rcode, nout )
c$$$            endif
c$$$            call writg1 ( ".............................", seps, 0,
c$$$     $           nout )
c$$$         enddo
c$$$         call writg1 ( "Variable's Information (7.3): ", seps, 1, nout )
c$$$         call writg1 ( "        Data Types:", seps, 1, nout )
c$$$         write ( 6, '(/,3x,"No.",15x,"Type",/, (i4, 2x, a72) )' )
c$$$     $        ( idat, datype(idat), idat = 1, 6 )
c$$$         write ( outmod, '(/,3x,"No.",15x,"Type",/, (i4, 2x, a72) )' )
c$$$     $        ( idat, datype(idat), idat = 1, 6 )
c$$$      endif
c$$$      do ivars  = 1, nvars
c$$$         call ncvinq ( cdfid, ivars, vname, vtype, vn(ivars),
c$$$     $        vdims, vnatt, rcode )
c$$$         vrname(ivars) = vname
c$$$         vrtype(ivars) = vtype
c$$$         vrnat(ivars) = vnatt
c$$$         if ( vn(ivars) .ne. 0 ) then
c$$$            do ivn = 1, vn(ivars)
c$$$               vvdims(ivars,ivn) = vdims(ivn)
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$      if ( check1 ) then
c$$$         write ( 6, 1001 )
c$$$         write ( outmod, 1001 )
c$$$         do ivars = 1, nvars
c$$$            write ( 6, 1002 ) ivars, vrname(ivars), vrtype(ivars),
c$$$     $           vrnat(ivars),
c$$$     $           ( dimsiz(vvdims(ivars,ivn)), ivn = 1, vn(ivars) )
c$$$            write ( outmod, 1002 ) ivars, vrname(ivars), vrtype(ivars),
c$$$     $           vrnat(ivars),
c$$$     $           ( dimsiz(vvdims(ivars,ivn)), ivn = 1, vn(ivars) )
c$$$         enddo
c$$$      endif
c$$$      call writg1 ( "Reading the Variables.", seps, 3, nout )
c$$$      vid = ncvid ( cdfid, 'ctitle', rcode )
c$$$      call nwvid ( cdfid, 'ctitle', vid, rcode, nout )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = maxc1
c$$$      call ncvgtc ( cdfid, vid, start,count, ctitle, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'date0', rcode )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = maxct
c$$$      call ncvgtc ( cdfid, vid, start,count, date0, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'time0', rcode )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = maxct
c$$$      call ncvgtc ( cdfid, vid, start,count, time0, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'datem', rcode )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = maxct
c$$$      call ncvgtc ( cdfid, vid, start,count, datem, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'timem', rcode )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = maxct
c$$$      call ncvgtc ( cdfid, vid, start,count, timem, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'dskout', rcode )
c$$$      start(1) = 1
c$$$      count(1) = dimsiz( vvdims(vid,1) )
c$$$      nd = 1
c$$$      length = nccl3
c$$$      call ncvgtc ( cdfid, vid, start,count, dskout, length,
c$$$     $     rcode )
c$$$      vid = ncvid ( cdfid, 'nx', rcode )
c$$$      vindx(1) = 1
c$$$      nd = vn(vid)
c$$$      call ncvgt1 ( cdfid, vid, vindx, nx, rcode )
c$$$      vid = ncvid ( cdfid, 'nz', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, nz, rcode )
c$$$      vid = ncvid ( cdfid, 'nosurf', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, nosurf, rcode )
c$$$      vid = ncvid ( cdfid, 'mth', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, mthin, rcode )
c$$$      vid = ncvid ( cdfid, 'lpless', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, lpless, rcode )
c$$$      vid = ncvid ( cdfid, 'alx', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, alx, rcode )
c$$$      vid = ncvid ( cdfid, 'alz', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, alz, rcode )
c$$$      vid = ncvid ( cdfid, 'xzpst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, xzero, rcode )
c$$$      vid = ncvid ( cdfid, 'xmag', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, xma, rcode )
c$$$      vid = ncvid ( cdfid, 'rpst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, r, rcode )
c$$$      vid = ncvid ( cdfid, 'p0pst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, p0, rcode )
c$$$      vid = ncvid ( cdfid, 'gppst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, gp0, rcode )
c$$$      vid = ncvid ( cdfid, 'pminpst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, psimin, rcode )
c$$$      vid = ncvid ( cdfid, 'plimpst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, psilim, rcode )
c$$$      vid = ncvid ( cdfid, 'plspst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, psipls, rcode )
c$$$      vid = ncvid ( cdfid, 'betapst', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betag, rcode )
c$$$      vid = ncvid ( cdfid, 'betag', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betap, rcode )
c$$$      vid = ncvid ( cdfid, 'rjacps', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, rjacps, rcode )
c$$$      if ( lpest1 ) then
c$$$         vid = ncvid ( cdfid, 'upsiln', rcode )
c$$$         call ncvgt1 ( cdfid, vid, vindx, upsiln, rcode )
c$$$      endif
c$$$      vid = ncvid ( cdfid, 'mjacx', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, mj, rcode )
c$$$      vid = ncvid ( cdfid, 'njacg', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, nj, rcode )
c$$$      vid = ncvid ( cdfid, 'ljacb', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, lj, rcode )
c$$$      vid = ncvid ( cdfid, 'nzd1', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, nzd1map, rcode )
c$$$      vid = ncvid ( cdfid, 'nzd2', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, nzd2map, rcode )
c$$$      vid = ncvid ( cdfid, 'dth', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, dtmap, rcode )
c$$$      vid = ncvid ( cdfid, 'dr', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, drmap, rcode )
c$$$      vid = ncvid ( cdfid, 'pi', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, pimap, rcode )
c$$$      vid = ncvid ( cdfid, 'xzero', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, xzeromap, rcode )
c$$$      vid = ncvid ( cdfid, 'zmag', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, zma, rcode )
c$$$      vid = ncvid ( cdfid, 'upsiln', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, upsiln, rcode )
c$$$      vid = ncvid ( cdfid, 'betat', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betat, rcode )
c$$$      vid = ncvid ( cdfid, 'betap', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betap, rcode )
c$$$      vid = ncvid ( cdfid, 'ctroy', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, ctroy, rcode )
c$$$      vid = ncvid ( cdfid, 'betatot', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betatot, rcode )
c$$$      vid = ncvid ( cdfid, 'betapo', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betapo, rcode )
c$$$      vid = ncvid ( cdfid, 'betato', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betato, rcode )
c$$$      vid = ncvid ( cdfid, 'betats', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, betats, rcode )
c$$$      vid = ncvid ( cdfid, 'btor', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, btor, rcode )
c$$$      vid = ncvid ( cdfid, 'aima', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, aima, rcode )
c$$$      vid = ncvid ( cdfid, 'aioab', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, aioab, rcode )
c$$$      vid = ncvid ( cdfid, 'bt2dv', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, bt2dv, rcode )
c$$$      vid = ncvid ( cdfid, 'pdvamu0', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, pdvamu0, rcode )
c$$$      vid = ncvid ( cdfid, 'dv', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, dv, rcode )
c$$$      vid = ncvid ( cdfid, 'xliga', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, xliga, rcode )
c$$$      vid = ncvid ( cdfid, 'dlp', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, dlp, rcode )
c$$$      vid = ncvid ( cdfid, 'qstar', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, qstar, rcode )
c$$$      vid = ncvid ( cdfid, 'qratio', rcode )
c$$$      call ncvgt1 ( cdfid, vid, vindx, qratio, rcode )
c$$$      vid = ncvid ( cdfid, 'q', rcode )
c$$$      call ncvgt1 ( cdfid, vid, 1, q00, rcode )
c$$$      vid = ncvid ( cdfid, 'p', rcode )
c$$$      call ncvgt1 ( cdfid, vid, 1, p01, rcode )
c$$$      vid = ncvid ( cdfid, 'fb', rcode )
c$$$      call ncvgt1 ( cdfid, vid, 1, f0, rcode )
c$$$      vid = ncvid ( cdfid, 'fb', rcode )
c$$$      call ncvgt1 ( cdfid, vid, nosurf, fa1, rcode )
c$$$      vid = ncvid ( cdfid, 'g', rcode )
c$$$      call ncvgt1 ( cdfid, vid, 1, g0, rcode )
c$$$      vid = ncvid ( cdfid, 'g', rcode )
c$$$      call ncvgt1 ( cdfid, vid, nosurf, ga1, rcode )
c$$$      vid = ncvid ( cdfid, 'q', rcode )
c$$$      call ncvgt1 ( cdfid, vid, nosurf, qa1, rcode )
c$$$      q0 = q00
c$$$      q1 = qa1
c$$$      ljacb = lj
c$$$      mjacx = mj
c$$$      njacg = nj
c$$$      write ( outmod, 9000 ) ctitle, jobid
c$$$      write ( iotty,  9000 ) ctitle, jobid
c$$$      write(outmod,2002) datev,timev, datem,timem, date0,time0,
c$$$     $     dskout, nosurf, mthin,
c$$$     $     mjacx, njacg, ljacb,
c$$$     $     mjacx, njacg, ljacb
c$$$      write( 6,2002) datev,timev, datem,timem, date0,time0,
c$$$     $     dskout, nosurf, mthin,
c$$$     $     mjacx, njacg, ljacb,
c$$$     $     mjacx, njacg, ljacb
c$$$      write(nout(2),3066) betap,betat,betatot
c$$$      write(6,3066) betap,betat,betatot
c$$$      write(outmod,3067) betapo,betato,betatot
c$$$      write(6,3067) betapo,betato,betatot
c$$$      write(outmod,3068) betats
c$$$      write(6,3068) betats
c$$$      write(outmod,3777) btor,aima,aioab,ctroy,q00,q1
c$$$      write(6,3777) btor,aima,aioab,ctroy,q00,q1
c$$$      write(outmod,3778) qstar,qratio
c$$$      write(6,3778) qstar,qratio
c$$$      write(outmod,3779) bt2dv,pdvamu0,dv
c$$$      write(6,3779) bt2dv,pdvamu0,dv
c$$$      write(outmod,3780) xliga,dlp
c$$$      write(6,3780) xliga,dlp
c$$$      write ( outmod,9600 ) qa1, fa1, ga1
c$$$      write ( iotty, 9600 ) qa1, fa1, ga1
c$$$      write ( outmod,9500 ) r, upsiln, xma
c$$$      write ( iotty, 9500 ) r, upsiln, xma
c$$$      mthin1 = mthin  + 1
c$$$      mthin2 = mthin1 + 1
c$$$      r2     = r * r
c$$$      r4     = r2 * r2
c$$$      r6     = r4 * r2
c$$$c-----------------------------------------------------------------------
c$$$c     termination.
c$$$c-----------------------------------------------------------------------
c$$$      return
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 7. mdskrd1.
c$$$c     read netcdf data.
c$$$c-----------------------------------------------------------------------
c$$$c-----------------------------------------------------------------------
c$$$c     declarations.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine mdskrd1
c$$$      USE vglobal_mod
c$$$      implicit real*8 (a-h,o-z)
c$$$
c$$$      include 'netcdf.inc'
c$$$      dimension vecin(ntsin)
c$$$c-----------------------------------------------------------------------
c$$$c     computations.
c$$$c-----------------------------------------------------------------------
c$$$      iref = 0
c$$$      mthin1 = mthin + 1
c$$$      mthin2 = mthin1 + 1
c$$$      vid = ncvid ( cdfid, 'x', rcode )
c$$$      nd = 2
c$$$      start(1) = 1
c$$$      start(2) = nosurf
c$$$      count(1) = dimsiz ( vvdims(vid,1) )
c$$$      count(2) = 1
c$$$      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $     rcode, nout0 )
c$$$      call trans ( vecin,mthin, xinf,mth )
c$$$      vid = ncvid ( cdfid, 'z', rcode )
c$$$      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $     rcode, nout0 )
c$$$      call trans ( vecin,mthin, zinf,mth )
c$$$      vid = ncvid ( cdfid, 'grpssq', rcode )
c$$$      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $     rcode, nout0 )
c$$$      call trans ( vecin,mthin, grpssq,mth )
c$$$      if ( .not. lpest1 ) then
c$$$         vid = ncvid ( cdfid, 'xjacob', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, xjacob,mth )
c$$$         vid = ncvid ( cdfid, 'delta', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, delta,mth )
c$$$      endif
c$$$      if ( lpest1 ) then
c$$$         do i = 1, mth2
c$$$            xjacob(i) = upsiln * xinf(i)**2 / ( twopi*r )
c$$$         enddo
c$$$      endif
c$$$      if ( lnova ) then
c$$$         vid = ncvid ( cdfid, 'xsq', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, xsq,mth )
c$$$         vid = ncvid ( cdfid, 'gpsdth', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, gpsdth,mth )
c$$$         vid = ncvid ( cdfid, 'xsqdth', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, xsqdth,mth )
c$$$         vid = ncvid ( cdfid, 'xjacob', rcode )
c$$$         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c$$$         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c$$$     $        rcode, nout0 )
c$$$         call trans ( vecin,mthin, xjacob,mth )
c$$$         mthd2p1=mth/2+1
c$$$         do i = 1, mthd2p1
c$$$            zbsq = ( grpssq(i) + r2*ga1 ) / xsq(i)
c$$$            zdbsqb = ( gpsdth(i) - zbsq*xsqdth(i) ) / (xsq(i)*zbsq )
c$$$            xsdtxs(i) = xsqdth(i)/xsq(i)
c$$$            gpdtgp(i) = gpsdth(i)/grpssq(i)
c$$$            xjdtxj(i) = 0.5 * ( mj*xsdtxs(i)-nj*gpdtgp(i) - lj*zdbsqb )
c$$$            xjacob(mth2-i)=xjacob(i)
c$$$            xjdtxj(mth2-i)=-xjdtxj(i)
c$$$            delta(mth2-i)=-delta(i)
c$$$         enddo
c$$$         xjdtxj(1)=0.
c$$$         xjdtxj(mthd2p1)=0.
c$$$         delta(1)=0.
c$$$         delta(mthd2p1)=0.
c$$$      endif
c$$$c-----------------------------------------------------------------------
c$$$c     termination.
c$$$c-----------------------------------------------------------------------
c$$$      return
c$$$      end
c-----------------------------------------------------------------------
c     subprogram 8. adjustm.
c     read netcdf data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine adjustm ( mth, mfel, mth1,mth2, ndfel, nout1,nout2 )
      implicit real*8 (a-h,o-z)
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
         write ( nout1,
     $        '(/,1x,"****** WARNING: mth00 .ne. mth ******",/ )' )
         write ( nout2,
     $        '(/,1x,"****** WARNING: mth00 .ne. mth ******",/ )' )
         mth1 = mth + 1
         mth2 = mth1 + 1
      endif
      write ( nout1, '(/,5x, "mth00, mth, mfel, ndfel = ", 4i5 )' )
     $     mth00, mth, mfel, ndfel
      write ( nout2, '(/,5x, "mth00, mth, mfel, ndfel = ", 4i5 )' )
     $     mth00, mth, mfel, ndfel
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
