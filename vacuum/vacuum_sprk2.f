c-----------------------------------------------------------------------
c     file vacuum_sprk2.f
c     seems to contain miscellaneous utilities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. testvec
c     2. vecpot
c     3. grchi
c     4. diff5
c     5. difspl
c     6. matwrtn
c     7. matwrt9
c     8. mtrans
c     9. mtransr
c-----------------------------------------------------------------------
c     subprogram 1. testvec.
c     tests some sort of vector for some sort of property.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine testvec 
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension xob(nths),zob(nths), the(nths)
      dimension axrl(nths,nfm), axil(nths,nfm),
     $     azrl(nths,nfm), azil(nths,nfm), bdotnr(nfm), bdotni(nfm)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 8    format ( "n = ",f8.2, " phi = ",f8.3, " mdeld = ", i3,
     $     " jdeld = ", i3, / )
 11   format ( i4,("   axrl ",1p10e11.4,/) )
 12   format ( i4,("   axil ",1p10e11.4,/) )
 15   format ( i4,("   azrl ",1p10e11.4,/) )
 16   format ( i4,("   azil ",1p10e11.4,/) )
 13   format ( i4,(" bdotnr ",1p10e11.4,/) )
 33   format ( i4,(" bdotni ",1p10e11.4,/) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
      mdel = 16
      jdel = 8
      do i = 1, mth1
         the(i) = (i-1) * dth
      enddo
      call bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj  = 0.5 * ( xmxp + xmnp )
      call bounds(xwal,zwal,1,mw,xmnw,xmxw,zmnw,zmxw)
      xmin = min(xmnp,xmnw)
      xmax = max(xmxp,xmxw)
      zmin = min(zmnp,zmnw)
      zmax = max(zmxp,zmxw)
      dx = xmax - xmin
      dz = 2.0 * zmax
      dxz = max(dx,dz)
      xmaxw = 0.5 * ( xmax + xmin ) + 0.50*dxz
      xminw = 0.5 * ( xmax + xmin ) - 0.50*dxz
      zmaxw = 0.50 * dxz
      zminw = - zmaxw
      xming = xminw - 0.05*(xmaxw-xminw)
      xmaxg = xmaxw + 0.05*(xmaxw-xminw)
      zming = zminw - 0.05*(zmaxw-zminw)
      zmaxg = zmaxw + 0.05*(zmaxw-zminw)
      phi = 30.0 * pye / 180.0
      nobs = 0
      mthdel = mth/mdel
      jmxdel = jmax1/jdel
      mdeld = max0(1,mthdel)
      jdeld = max0(1,jmxdel)
      write ( outmod,'(a)' ) char(12)
      write ( outmod, 8 ) n, phi, mdeld, jdeld
      write ( iotty, 8 ) n, phi, mdeld, jdeld
      do iobs = 1, mth, mdeld
         nobs = nobs + 1
         xob(nobs) = xwal(iobs)
         zob(nobs) = zwal(iobs)
      enddo
      if (checks )then
         isg = -1
         call vecpot ( xwal,zwal,mth, xpla,zpla,xplap,zplap,
     $        mth, isg, phi,axrl,axil, azrl,azil )
         do jj = 1, jmax1
            axrl(mth1,jj) = axrl(1,jj)
            axil(mth1,jj) = axil(1,jj)
            azrl(mth1,jj) = azrl(1,jj)
            azil(mth1,jj) = azil(1,jj)
         enddo
         do iobs = 1, mth, mdeld
            write ( outmod, 11 ) iobs, (axrl(iobs,ll), ll = 1,jmax1
     $           ,jdeld)
            write ( outmod, 12 ) iobs, (axil(iobs,ll), ll = 1,jmax1
     $           ,jdeld)
            write ( outmod, 15 ) iobs, (azrl(iobs,ll), ll = 1,jmax1
     $           ,jdeld)
            write ( outmod, 16 ) iobs, (azil(iobs,ll), ll = 1,jmax1
     $           ,jdeld)
            do ll = 1, jmax1
               bdotnr(ll) = - n * ( axil(iobs,ll)*xwalp(iobs)
     $              + azil(iobs,ll)*zwalp(iobs) )
               bdotni(ll) =   n * ( axrl(iobs,ll)*xwalp(iobs)
     $              + azrl(iobs,ll)*zwalp(iobs) )
            enddo
            write ( outmod, 13 ) iobs, (bdotnr(ll), ll = 1,jmax1,jdeld)
            write ( outmod, 33 ) iobs, (bdotni(ll), ll = 1,jmax1,jdeld)
         enddo
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 2. vecpot.
c     compute vector potential.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vecpot ( xobs,zobs,nobs, xsce,zsce,xscp,zscp,nsce,
     $     isg, phi, axrl,axil, azrl,azil )
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension xobs(*),zobs(*),xsce(*),zsce(*),xscp(*),zscp(*),
     $     axrl(nths,nfm),axil(nths,nfm),
     $     azrl(nths,nfm),azil(nths,nfm),
     $     dchxr(nths,nfm),dchxi(nths,nfm),
     $     dchzr(nths,nfm),dchzi(nths,nfm)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 11   format ( /,1x, "!!! N= 0 in subroutine VECPOT !!! " )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      sinnph = sin(n*phi)
      cosnph = cos(n*phi)
      call grchi(xobs,zobs,nobs, xsce,zsce,xscp,zscp,nsce,
     $     isg,cplar,cplai,cslth,snlth, 1,1,
     $     cpwr,cpwi, dchxr,dchxi,dchzr,dchzi)
      if ( abs(n) .lt. 1.e-10 ) then
         write ( outmod, 11 ) 
         write ( iotty,  11 ) 
         return
      endif
      do iobs = 1, nobs
         do ll = 1, jmax1
            axrl(iobs,ll) = - xobs(iobs) * ( dchzr(iobs,ll)*sinnph
     $           - dchzi(iobs,ll)*cosnph ) / n
            axil(iobs,ll) = - xobs(iobs) * ( dchzr(iobs,ll)*cosnph
     $           + dchzi(iobs,ll)*sinnph ) / n
            azrl(iobs,ll) =   xobs(iobs) * ( dchxr(iobs,ll)*sinnph
     $           - dchxi(iobs,ll)*cosnph ) / n
            azil(iobs,ll) =   xobs(iobs) * ( dchxr(iobs,ll)*cosnph
     $           + dchxi(iobs,ll)*sinnph ) / n
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. grchi.
c     2nd application of Green's identity.  Returns chi, Bx, Bz.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine grchi(xobs,zobs,nobs, xsce,zsce,xscp,zscp,ns,
     $     isg,creal,cimag,cdriv,sdriv, ip,ipm, cpwr1,cpwi1,
     $     dchxr,dchxi,dchzr,dchzi)
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension xobs(*),zobs(*), xsce(*),zsce(*),xscp(*),zscp(*)
      dimension creal(nths,nfm), cimag(nths,nfm)
      dimension cdriv(nths,nfm), sdriv(nths,nfm)
      dimension chrl(5,nfm), chil(5,nfm), wrkr(5),wrki(5)
      dimension cpwr1(nths,nfm),cpwi1(nths,nfm),
     $     dchxr(nths,nfm),dchxi(nths,nfm),
     $     dchzr(nths,nfm),dchzi(nths,nfm)
      real nq
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      factpi = twopi
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1
      nq = n * q
      dtpw = twopi / ns
      do io = 1, nobs
         do ixz = 1, 2
            yes1 = 2 - ixz
            yes2 = ixz - 1
            do io5 = 1, 5
               do ll = 1, jmax1
                  chrl(io5,ll) = 0.0
                  chil(io5,ll) = 0.0
               enddo
            enddo
            do io5 = 1, 5
               xs = yes1*(xobs(io) + (io5-3)*delx) + yes2*xobs(io)
               zs = yes2*(zobs(io) + (io5-3)*delz) + yes1*zobs(io)
     $              
               if ( (io5 .ne. 3) .or. (ixz .ne. 2) ) then
                  do is = 1, ns
                     xt = xsce(is)
                     zt = zsce(is)
                     xtp = xscp(is)
                     ztp = zscp(is)
                     call green
                     bval = bval / factpi
                     do l1 = 1, jmax1
                        ll = lmin(1) - 1 + l1
                        if ( ip .ne. 3 ) then
                           chrl(io5,l1) = chrl(io5,l1)  + creal(is,l1)
     $                          *aval
                           chil(io5,l1) = chil(io5,l1)  + cimag(is,l1)
     $                          *aval
                        endif
                        if ( ip .ne. 0 ) then
                           chrl(io5,l1)  = chrl(io5,l1)  +
     $                          ipm * bval * cdriv(is,l1)
                           chil(io5,l1)  = chil(io5,l1)  +
     $                          ipm * bval * sdriv(is,l1)
                        endif
                     enddo
                  enddo
                  do l1 = 1, jmax1
                     chrl(io5,l1)  = 0.5 * isg*dtpw * chrl(io5,l1) 
                     chil(io5,l1)  = 0.5 * isg*dtpw * chil(io5,l1) 
                  enddo
               endif
            enddo
            do l1 = 1, jmax1
               do id = 1, 5
                  wrkr(id) = chrl(id,l1)
                  wrki(id) = chil(id,l1)
               enddo
               if ( yes1 .eq. 1 ) then
                  cpwr1(io,l1) = chrl(3,l1)
                  cpwi1(io,l1) = chil(3,l1)
                  call diff5 ( wrkr, delx, dchxr(io,l1) )
                  call diff5 ( wrki, delx, dchxi(io,l1) )
               endif
               if ( yes2 .eq. 1 ) then
                  call diff5 ( wrkr, delz, dchzr(io,l1) )
                  call diff5 ( wrki, delz, dchzi(io,l1) )
               endif
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. diff5.
c     5-point difference stencil.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine diff5 ( fin, h, df5 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension fin(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      fm2 = fin(1)
      fm1 = fin(2)
      fp1 = fin(4)
      fp2 = fin(5)
      df5 = 2.0 * ( fp1-fm1 -(fp2-fm2)/8.0 ) / (3.0*h)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 5. difspl.
c     unknown purpose.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine difspl ( nsp,the,xin,xout )
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension  the(*), xin(*), xout(*)
      dimension xpp(nths)
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      ns1 = nsp + 1
      iop(1) = 4
      iop(2) = 4
      call spl1d1(ns1,the,xin,xpp,iop,1,ww1,ww2,ww3)
      do i = 1, ns1
         theta = (i-1) * dth
         call spl1d2 ( ns1,the, xin,xpp, 1, theta, tab )
         xout(i) = tab(2)
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 6. matwrtn.
c     write a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matwrtn ( a, maxj1,maxj2,l1,l2,jmax1,jmax2,jn1,jn2,
     $     label, nout1,nout2 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
 20   format ( i4,10(1x,1p10e11.4,/) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nwrt2 = 9
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, jn1, jn2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, jn1, jn2
      jd1 = jmax1 / jn1
      jd2 = jmax2 / jn2
      if ( jd1 .lt. 1 ) jd1 = 1
      if ( jd2 .lt. 1 ) jd2 = 1
c      njd1 = min ( jmax1,(nwrt1-1)*jd1 + 1 )
      njd2 = min ( jmax2,(nwrt2-1)*jd2 + 1 )
      if ( jn2 .le. 101 ) then
         j2m = 0
         do j2 = 1, jmax2,jd2
            j2m = j2m + 1
            jw2(j2m) = j2 - 1 + l2
         enddo
         j2m = min(nwrt2,j2m)
         if ( nout1. ne. 0 )
     $        write ( nout1, '(10i11)' ) ( jw2(i),i = 1, j2m )
         if ( nout2. ne. 0 )
     $        write ( nout2, '(10i11)' ) ( jw2(i),i = 1, j2m )
      endif
      do j1 = 1, jmax1,jd1
         jw = j1 -1 + l1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 ) jw, ( a(j1,j2),j2 = 1, njd2,jd2 )
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 ) jw, ( a(j1,j2),j2 = 1, njd2,jd2 )
      enddo
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 7. matwrt9.
c     write a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matwrt9 ( a, maxj1,maxj2,l1,l2, label, nout1,nout2 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
 20   format ( i4,10(1x,1p10e11.4,/) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, l1, l2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, l1, l2
      ln4 = max ( l1, -4 )
      lx4 = min ( l2, +4 )
      jnx = lx4 - ln4 + 1
      do j2 = 1, jnx
         jw2(j2) = j2 - 1 + ln4
      enddo
      if ( nout1. ne. 0 )
     $     write ( nout1, '(10i11)' ) ( jw2(i),i = 1, jnx )
      if ( nout2. ne. 0 )
     $     write ( nout2, '(10i11)' ) ( jw2(i),i = 1, jnx )
      js = max ( -4-l1, 0 ) + 1
      lw = 0
      do j1 = js, js+jnx-1
         lw = lw + 1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 ) jw2(lw),(a(j1,j2),j2 = js, js+jnx-1)
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 ) jw2(Lw),(a(j1,j2),j2 = js, js+jnx-1)
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 8. mtrans.
c     transpose a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mtrans ( a, nd, n )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension a(nd,nd)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do i = 1, n
         do j = 1, i
            z = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = z
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 9. mtransr.
c     transpose a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mtransr ( a, nd1,nd2, n1,n2, b )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension a(nd1,nd2), b(nd2,nd1)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do i = 1, n1
         do j = 1, n2
            b(j,i) = a(i,j)
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
