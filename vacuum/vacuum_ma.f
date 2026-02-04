c-----------------------------------------------------------------------
c     file vacuum_ma.f.
c     Main routines for Morrell Chance's vacuum eigenvalue code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. mscvac
c     2. mscfld
c     3. defglo
c     4. ent33
c     5. funint
c     6. make_bltobp
c     7. diaplt
c     8. pickup
c     9. loop
c    10. chi

      module vacuum_mod
      ! implicit none
      contains


c-----------------------------------------------------------------------
c     subprogram 1. mscvac.
c     calculate vacuum response matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mscvac(wv,mpert,mtheta,mthvac,complex_flag,
     $     kernelsignin,wall_flag,farwal_flag,grrio,xzptso,op_ahgfile)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      REAL(r8) :: kernelsignin
      integer mpert,mtheta,mthvac
      complex(r8) wv(mpert,mpert)
      logical, intent(in) :: complex_flag,wall_flag,farwal_flag
      REAL(r8) :: grrio(2*(mthvac+5),mpert*2),xzptso(mthvac+5,4)

      complex(r8), parameter :: ifac=(0,1)
      dimension xi(nfm), xii(nfm), xilnq(nfm), xiilnq(nfm)
      character(128), intent(in), optional :: op_ahgfile

      if (present(op_ahgfile)) then
         ahgfile = trim(op_ahgfile)
      else
         ahgfile = 'ahg2msc.out'
      endif


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
      kernelsign=kernelsignin
      ntsin0=mtheta+1
      nths0=mthvac
      nfm=mpert
      mtot=mpert
      call global_alloc(nths0,nfm,mtot,ntsin0)
      farwal=.false.
      IF (farwal_flag) farwal=.true.
c-----------------------------------------------------------------------
c     initialization.
c-----------------------------------------------------------------------
      call defglo(mthvac)
      ldcon = 1
      lgpec = 0
      open (iotty,file='mscvac.out',status='unknown')
      open (outpest,file='pestotv',status='unknown',form='formatted')
      IF (wall_flag) THEN
         open (inmode,file='vac_wall.in',status='old', form='formatted')
      ELSE
         open (inmode,file='vac.in',status='old', form='formatted' )
      ENDIF
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
      close(outpest)
      close(inmode)
      close(outmod)
      grrio(:,:)=grri(:,:)
      xzptso(:,1)=xinf(:)
      xzptso(:,2)=zinf(:)
      xzptso(:,3)=xwal(:)
      xzptso(:,4)=zwal(:)
      call global_dealloc
      call cleanup
      return
      end
c-----------------------------------------------------------------------
c     subprogram 2. mscfld.
c     calculate vacuum field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mscfld(wv,mpert,mtheta,mthvac,complex_flag,
     $     lx,lz,vgdl,vgdx,vgdz,vbx,vbz,vbp,op_ahgfile)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      integer mpert,mtheta,mthvac
      complex(r8) wv(mpert,mpert)
      logical, intent(in) :: complex_flag

      integer, dimension(0:lx,0:lz) :: vgdl
      REAL(r8), dimension(0:lx,0:lz) :: vgdx,vgdz
      complex(r8), dimension(0:lx,0:lz) :: vbx,vbz,vbp

      complex(r8), parameter :: ifac=(0,1)
      dimension xi(nfm), xii(nfm), xilnq(nfm), xiilnq(nfm)
      character(128), intent(in), optional :: op_ahgfile

      if (present(op_ahgfile)) then
         ahgfile = trim(op_ahgfile)
      else
         ahgfile = 'ahg2msc.out'
      endif

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
      kernelsign=1
      ntsin0=mtheta+1
      nths0=mthvac
      nfm=mpert
      mtot=mpert
      call global_alloc(nths0,nfm,mtot,ntsin0)
      ndimlp=(lx+1)*(lz+1)
      ALLOCATE(xobp(ndimlp),zobp(ndimlp),xloop(ndimlp),zloop(ndimlp))
c-----------------------------------------------------------------------
c     initialization.
c-----------------------------------------------------------------------
      call defglo(mthvac)
      ldcon = 0
      lgpec = 1
      ieig = 0
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
      IF(complex_flag)THEN
         wv=vacmat+ifac*vacmtiu
      ELSE
         wv=vacmat
      ENDIF
c-----------------------------------------------------------------------
c     compute b field.
c-----------------------------------------------------------------------
 99   continue
      call make_bltobp
      call diaplt
      call pickup(bnlr,bnli,lx,lz,vgdl,vgdx,vgdz,vbx,vbz,vbp)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      call msctimer ( outmod, "end of main" )
      close(iotty)
      close(outpest)
      close(inmode)
      close(outmod)
      call global_dealloc
      call cleanup
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. defglo.
c     defines logical unit numbers and control switches.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine defglo(mthvac)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer, intent(in) :: mthvac
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
      mth    = mthvac
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthin  = mthvac
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
c     subprogram 4. ent33.
c     entry point for segment 33.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine ent33
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 5. funint.
c     function initialization.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine funint
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

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
         gpsjp(i) = xpla(i) * dlenth(i)
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
c-----------------------------------------------------------------------
c     subprogram 6. make_bltobp
c     calculate normal field in real space.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine make_bltobp
      USE vglobal_mod
      IMPLICIT REAL (a-h,o-z)

      ! DIMENSION chlagdy(nths,nfm)
      DIMENSION thmgr(nths), z1tmp(nths), z2tmp(nths)

      psipr = 1.0
      lrnge = nfm
      DO  l = 1, lrnge
         psilnq = psipr * ( lfm(l)-n*qa1 )
         xilr(l) =  bnli(l) / psilnq
         xili(l) = -bnlr(l) / psilnq
      END DO

      DO i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         DO l = 1, lrnge
            zgr = zgr + xilr(l)*cslth(i,l) - xili(l)*snlth(i,l)
            zgi = zgi + xilr(l)*snlth(i,l) + xili(l)*cslth(i,l)
         END DO
         xigr(i) = zgr
         xigi(i) = zgi
      END DO

      DO i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         DO l = 1, lrnge
            zgr = zgr + bnlr(l)*cslth(i,l) - bnli(l)*snlth(i,l)
            zgi = zgi + bnlr(l)*snlth(i,l) + bnli(l)*cslth(i,l)
         END DO
         bnkr(i) = zgr
         bnki(i) = zgi
      END DO

      l11 = lmin(1)
      l22 = lmax(1)
      CALL fanal ( xigr, mth, xirc, xirs, l11,l22, pye,0.0 )
      CALL fanal ( xigi, mth, xiic, xiis, l11,l22, pye,0.0 )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 7. diaplt
c     calculate tangential b field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine diaplt
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

   !    DIMENSION z1tmp(nths), z2tmp(nths), zorkr(nths),zorki(nths),
   !   $     zorkpr(nths), zorkpi(nths), zork3(nths), !chlagdy(nths,nfm),
   !   $     thph(nths), cppgr(nths),cppgi(nths),
   !   $     cplgr(nths), cplgi(nths), cplgtr(nths), cplgti(nths),
   !   $     chwr1(nths),chwi1(nths),
   !   $     dxdt(nths), dzdt(nths), zkt(nths,2), zkp(nths,2)
      DIMENSION zork3(nths), zorkpr(nths), zorkpi(nths)

      lrnge = nfm
      jmax1 = lrnge
      laxis = lrnge
      IF ( laxis == 1 ) laxis =2

      l11 = lmin(1)

      do  i = 1, mth1
         bnpptr(i) = 0.0
         bnppti(i) = 0.0
         do j1 = 1, lrnge
            bnpptr(i) = bnpptr(i) + ( bnlr(j1)*cslth(i,j1)
     $           - bnli(j1)*snlth(i,j1) ) / gpsjp(i)
            bnppti(i) = bnppti(i) + ( bnlr(j1)*snlth(i,j1)
     $           + bnli(j1)*cslth(i,j1) )/ gpsjp(i)
         end do
      end do

      do i = 1, mth
         chipr(i) = 0.0
         chipi(i) = 0.0
         do  j1 = 1, jmax1
            zrrc = cplar(i,j1)
            zrrs = cplai(i,j1)
            chipr(i) = chipr(i) + zrrc*bnlr(j1) - zrrs*bnli(j1)
            chipi(i) = chipi(i) + zrrs*bnlr(j1) + zrrc*bnli(j1)
         end do
      end do
      chipr(mth1) = chipr(1)
      chipi(mth1) = chipi(1)

      DO i = 1, mth1
         CALL lagp ( slngth,chipr, mth1,5, slngth(i), f,zorkpr(i), 1,1 )
         CALL lagp ( slngth,chipi, mth1,5, slngth(i), f,zorkpi(i), 1,1 )
         zork3(i) = sqrt ( zorkpr(i)**2 + zorkpi(i)**2 )
      END DO
      bthpr = zorkpr
      bthpi = zorkpi
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 8. pickup
c     calculate vacuum b field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine pickup(blr,bli,lx,lz,vgdl,vgdx,vgdz,vbx,vbz,vbp)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      INTEGER, DIMENSION(0:lx,0:lz) :: vgdl
      REAL(r8), DIMENSION(0:lx,0:lz) :: vgdx,vgdz
      COMPLEX(r8), DIMENSION(0:lx,0:lz) :: vbx,vbz,vbp

      DIMENSION blr(*), bli(*)
      ! CHARACTER(130), DIMENSION(10) :: string
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      DIMENSION xloops(ndimlp),zloops(ndimlp),
     $     chir(5,ndimlp),chii(5,ndimlp),
     $     zchipr(ndimlp),zchipi(ndimlp),chirr(ndimlp),
     $     cwrkr(5,ndimlp),cwrki(5,ndimlp),
     $     igdl(ndimlp),rgdl(ndimlp),rwall(ndimlp),
     $     bxr(ndimlp),bxi(ndimlp),bzr(ndimlp),bzi(ndimlp),
     $     btr(ndimlp),bti(ndimlp), bpr(ndimlp),bpi(ndimlp),
     $     bphir(ndimlp),bphii(ndimlp)

      ndlp = mth / ntloop

      jmax1 = lmax(1) - lmin(1) + 1
      mth2 = 2*mth
      IF ( farwal ) mth2 = mth

      call bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp

      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )

      do i = 1, ndimlp
         bxr(i) = 0.0
         bxi(i) = 0.0
         bzr(i) = 0.0
         bzi(i) = 0.0
         bphir = 0.0
         bphii = 0.0
         chirr = 0.0
         btr(i) = 0.0
         bti(i) = 0.0
         bpr(i) = 0.0
         bpi(i) = 0.0
      end do

      if ( .not. farwal ) then
        call bounds(xwal,zwal,1,mw,xmnw,xmxw,zmnw,zmxw)
        xmin = min(xmnp,xmnw)
        xmax = max(xmxp,xmxw)
        zmin = min(zmnp,zmnw)
        zmax = max(zmxp,zmxw)
      endif

      CALL loops

      nobs = nloop + 3*nloopr

      call bounds ( xloop,zloop,1, nobs, xmn,xmx, zmn,zmx )
      xmin = min( xmin,xmn )
      xmax = max( xmax,xmx )
      zmin = min( zmin,zmn )
      zmax = max( zmax,zmx )

      dtpw = dth
      ns = mth
      delx = plrad * deloop
      delz = plrad * deloop

      igdl = 8
      isgchi = -1

      DO i = 1, nobs

         xloops(i) = xloop(i)
         zloops(i) = zloop(i)

         fintjj = 0.0
         DO jj  = 1, mth
            dxjj = xloop(i) - xinf(jj)
            dzjj = zloop(i) - zinf(jj)
            rhoj2 = dxjj**2 + dzjj**2

            IF ( rhoj2 < 1.0e-16 ) GO TO 238
            fintjj = fintjj +
     $           ( zplap(jj)*dxjj - xplap(jj)*dzjj ) / rhoj2
         END DO
         fintjj = fintjj / mth

         IF ( fintjj > 0.1 ) THEN ! interior
            igdl(i) = 1
            rgdl(i) = -1.0

            IF ( linterior == 0 ) GO TO 239
         END IF
         IF ( fintjj < 0.1 ) THEN ! exterior
            igdl(i) = 0
            rgdl(i) = 1.0

            IF ( linterior == 1 ) GO TO 239
         END IF
         rwall(i) = -1.0

 238     CONTINUE

         DO j = 1, mth
            IF ( (xinf(j)-xloop(i))**2 + (zinf(j)-zloop(i))**2
     $           < (epslp*plrad)**2 )  THEN ! points on surface for these.
               igdl(i) = -1

               xloops(i) = xinf(j)
               zloops(i) = zinf(j)

               isgchix = int(rgdl(i))*isgchi
               k = max(1,j-1)
               alph = atan2m (xinf(j+1)-xinf(k),zinf(k)-zinf(j+1))
               cosalph = COS(alph)
               sinalph = SIN(alph)
               bxr(i) = bnkr(j)*cosalph/gpsjp(j) +
     $              isgchix*bthpr(j)*sinalph
               bxi(i) = bnki(j)*cosalph/gpsjp(j) +
     $              isgchix*bthpi(j)*sinalph
               bzr(i) = bnkr(j)*sinalph/gpsjp(j) -
     $              isgchix*bthpr(j)*cosalph
               bzi(i) = bnki(j)*sinalph/gpsjp(j) -
     $              isgchix*bthpi(j)*cosalph
               zchipr(i) = isgchix*chipr(j)
               zchipi(i) = isgchix*chipi(j)
               bphir(i) =   n * zchipi(i) / xloops(i)
               bphii(i) = - n * zchipr(i) / xloops(i)
            END IF
         END DO

 239     CONTINUE

      ENDDO

      chir = 0.0
      chii = 0.0

      DO NSEW = 1, 5

         DO i = 1, nobs

            chir(nsew,i) = 0.0
            chii(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0

            go to ( 51, 52, 53, 54, 55), nsew

 51         continue

            xobp(i) = xloop(i)
            zobp(i) = zloop(i) + delz
            go to 90

 52         continue

            xobp(i) = xloop(i)
            zobp(i) = zloop(i) - delz
            go to 90

 53         continue

            xobp(i) = xloop(i) + delx
            zobp(i) = zloop(i)
            go to 90

 54         continue

            xobp(i) = xloop(i) - delx
            zobp(i) = zloop(i)
            GO TO 90

 55         CONTINUE

            xobp(i) = xloop(i)
            zobp(i) = zloop(i)

 90         CONTINUE

         END DO

         do l1 = 1, jmax1
            do i = 1, mth
               chiwc(i,l1) = grri(i,l1)
               chiws(i,l1) = grri(i,jmax1+l1)
            end do
         end do

         isg = -1

         call chi ( xpla,zpla,xplap,zplap,isg, chiwc,chiws, ns,1,
     $        cwrkr,cwrki,nsew, blr,bli,rgdl )

         do  i = 1, nobs
            chir(nsew,i) = cwrkr(nsew,i)
            chii(nsew,i) = cwrki(nsew,i)
         end do

         IF ( .not. farwal ) THEN
           do l1 = 1, jmax1
             do i = 1, mw
               chiwc(i,l1) = grri(mth+i,l1)
               chiws(i,l1) = grri(mth+i,jmax1+l1)
             end do
           end do

           do i = 1, nobs
             cwrkr(nsew,i) = 0.0
             cwrki(nsew,i) = 0.0
           end do
           
           isg = 1
           
           call chi ( xwal,zwal,xwalp,zwalp,isg,chiwc,chiws, ns,0,
     $          cwrkr,cwrki,nsew, blr,bli,rwall )

           do  i = 1, nobs
             chir(nsew,i) = chir(nsew,i) + cwrkr(nsew,i)
             chii(nsew,i) = chii(nsew,i) + cwrki(nsew,i)
           end do
         ENDIF !.not. farwal

      END DO

      DO i = 1, nobs
        IF (igdl(i) .NE. -1) THEN
          bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
          bxi(i) = ( chii(3,i) - chii(4,i) ) / (2.0*delx)
          bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
          bzi(i) = ( chii(1,i) - chii(2,i) ) / (2.0*delz)
          zchipr(i) = chir(5,i)
          zchipi(i) = chii(5,i)
          bphir(i) =   n * chii(5,i) / xloop(i)
          bphii(i) = - n * chir(5,i) / xloop(i)
        ENDIF
      END DO

      DO i = 1, nxlpin
         DO j = 1, nzlpin
            indxl = (i-1)*nzlpin + j
            vgdl(i-1,j-1) = igdl(indxl)
            vgdx(i-1,j-1) = xloops(indxl)
            vgdz(i-1,j-1) = zloops(indxl)
            vbx(i-1,j-1) = bxr(indxl)+ifac*bxi(indxl)
            vbz(i-1,j-1) = bzr(indxl)+ifac*bzi(indxl)
            vbp(i-1,j-1) = bphir(indxl)+ifac*bphii(indxl)
         END DO
      END DO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 9. loops.
c     grid loops.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine loops
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      REAL, DIMENSION(:,:), ALLOCATABLE :: xloopin, zloopin
      REAL, DIMENSION(:), ALLOCATABLE :: sourcemat


      dxlin = 1.0 / (nxlpin-1)
      dzlin = 1.0 / (nzlpin-1)
      nxzlin = nxlpin * nzlpin
      nloop = nxzlin

      ALLOCATE ( xloopin(nxlpin,nzlpin), zloopin(nxlpin,nzlpin) )
      ALLOCATE ( sourcemat(nxzlin) )

      sourcemat = (/ (((i-1)*dxlin, i=1,nxlpin),j=1,nzlpin) /)
      sourcemat = sourcemat * (xlpmax-xlpmin) + xlpmin

      xloopin(1:nxlpin,1:nzlpin) =
     $     RESHAPE ( sourcemat(1:nxzlin), SHAPE = (/ nxlpin,nzlpin /) )

      sourcemat = (/ (((j-1)*dzlin, i=1,nxlpin),j=1,nzlpin) /)
      sourcemat = sourcemat * (zlpmax-zlpmin) + zlpmin

      zloopin(1:nxlpin,1:nzlpin) =
     $     RESHAPE ( sourcemat(1:nxzlin), SHAPE = (/ nxlpin,nzlpin /) )

      xloop(1:nxzlin) = (/ ( (xloopin(i,j), j=1,nzlpin), i=1,nxlpin ) /)
      zloop(1:nxzlin) = (/ ( (zloopin(i,j), j=1,nzlpin), i=1,nxlpin ) /)

      DEALLOCATE ( xloopin, zloopin, sourcemat )

      GO TO 150

 150  CONTINUE

      IF ( nloopr == 0 ) go to 200

      IF ( a < 10.0 ) THEN
         drl = plrad * a / (nloopr)
      ELSE
         drl = 2.0 * plrad / nloopr
      ENDIF

      DO i = 1, nloopr
         xloop(nloop+i) = xmaj + plrad + i * drl - drl/2.0
         zloop(nloop+i) = zpla(ixx)
      END DO
c
      DO i = 1, nloopr
         il = nloop + nloopr + i
         xloop(il) = xmaj - plrad - i * drl + drl/2.0
         IF ( xloop(il) <= .01*xmaj ) xloop(il) = .01*xmaj
         zloop(il) = zpla(ixn)
      END DO
c
      DO i = 1, nloopr
         il = nloop + 2*nloopr + i
         xloop(il) = xpla(izx)
         zloop(il) = zmxp + i * drl - drl/2.0
      END DO

 200  CONTINUE
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 10. chi.
c     magnetic scalar potential.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine chi(xsce,zsce,xscp,zscp,isg,creal,cimag,ns,ip,
     $     chir,chii,nsew,blr,bli,rgdl)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      DIMENSION blr(*),bli(*),xsce(*),zsce(*),xscp(*),zscp(*)
      DIMENSION creal(nths,nfm), cimag(nths,nfm)
      DIMENSION chir(5,ndimlp), chii(5,ndimlp), rgdl(ndimlp)
      REAL nq

      factpi = twopi
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1

      nq = n * q
      dtpw = twopi / ns

      ns1 = ns + 1
      nobs = nloop + 3*nloopr

      do io = 1, nobs

         xs = xobp(io)
         zs = zobp(io)

         do is = 1, ns

            xt = xsce(is)
            zt = zsce(is)
            xtp = xscp(is)
            ztp = zscp(is)

            call green
            bval = bval / factpi

            do l1 = 1, jmax1
               zbr = blr(l1)
               zbi = bli(l1)
               chir(nsew,io) = chir(nsew,io) +
     $              aval * ( creal(is,l1)*zbr - cimag(is,l1)*zbi )
               chii(nsew,io) = chii(nsew,io) +
     $              aval * ( cimag(is,l1)*zbr + creal(is,l1)*zbi )

               if ( ip .eq. 0 ) go to 60

               chir(nsew,io) = chir(nsew,io) + rgdl(io) * bval
     $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
               chii(nsew,io) = chii(nsew,io) + rgdl(io) * bval
     $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )

 60            continue
            end do
         end do

         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
         chii(nsew,io) = 0.5 * isg*dtpw * chii(nsew,io)

      end do
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end

      end module vacuum_mod