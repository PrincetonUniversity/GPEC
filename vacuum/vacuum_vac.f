c-----------------------------------------------------------------------
c     file vacuum_vac.f.
c     main computation of vacuum energy in Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. vaccal
c     2. kernel
c     3. mateig
c     4. mateig2
c     5. arrays
c     6. wwall
c     7. d3dwall
c     8. d3dvesl
c     9. eqarcw
c     10. gatonorm
c     11. adjustb
c     12. fouran
c     13. foranv
c     14. foura2
c     15. fanal
c     16. fanal1
c     17. felang
c     18. felanv
c     19. fotofi
c     20. orchek
c     21. tmat
c     22. wtopest
c     23. savewall
c-----------------------------------------------------------------------
c     subprogram 1. vaccal.
c     solution of the vacuum integral equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vaccal
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      real nq
      integer tmth

      REAL(r8), DIMENSION(2) :: summ
      REAL(r8), DIMENSION(nfm) :: work0
      REAL(r8), DIMENSION(nfmsq) :: vacpstr,vacpsti,work,work1
      REAL(r8), DIMENSION(nfm,nfm) :: vacmti,wrkvr,wrkvi
      REAL(r8), DIMENSION(mtot,mtot) :: ajll,rmatr,rmati
      REAL(r8), DIMENSION(nths2,nfm2) :: dummy

      REAL(r8), DIMENSION(:,:), POINTER :: grdgre,arr,aii,ari,air
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE kernel(xobs,zobs,xsce,zsce,grdgre,gren,
     $        j1,j2,isgn,iopw,iops,ischk)
         USE vglobal_mod
         REAL(r8), DIMENSION(:), INTENT(IN) :: xobs,zobs,xsce,zsce
         REAL(r8), DIMENSION(:,:), INTENT(OUT) :: grdgre
         REAL(r8), DIMENSION(:,:), INTENT(OUT), TARGET :: gren
         INTEGER, INTENT(IN) :: j1,j2,isgn,iopw,iops,ischk
         INTEGER :: sTime, fTime, cr
         END SUBROUTINE kernel
      END INTERFACE
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 9002 format (  1x, i4, 3x, 1p5e14.5 )
 8180 format ( /, 5x, "constant,cn0, added to grdgre = ", f7.3 )
 8050 format (/, 1x, " ier in f04aae = ", i5 )
 8312 format ( " Writing Chi-i-l " )
 8311 format ( 2i5 )
 8313 format ( 1p10e14.6 )
 215  format ( /, 3x,"l",5x,"sumpc",5x,"sumps",5x,
     $     "sumwc",5x,"sumws" )
 217  format ( 1x, i4, 1p4e12.4 )
 8020 format ( 1x,/, "jtop, jbot, mw, icount = ", 4i5, / )
 500  format (//,1x,'n,q, nj,mj,lj = ',1p2e12.4,3i5,/ )
 501  format ( 1x, " delta =",/, (1x,10e11.4) )
 554  format ( 4i5, e13.5 )
 555  format ( 10e13.5 )
c-----------------------------------------------------------------------
c     zero local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(grdgre(nths2,nths2))
      grdgre=0
      grri=0
      ajll=0
      rmatr=0
      rmati=0
      dummy=0
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      write ( outmod, '(/,10x, "Matrix Storage: K(obs_ji,sou_ji):",//,
     $     10x, "j = observer points.  :i = source points.",/,
     $     10x, "  ie. K operates on chi from the left."//,
     $     10x, "Observer Source  Block",/,
     $     10x, " plasma  plasma   1  1",/,
     $     10x, " plasma  wall     1  2",/,
     $     10x, " wall    plasma   2  1",/,
     $     10x, " wall    wall     2  2",/ )' )
      xwal(1) = 0.
      zwal(1) = 0.
      ier = 0
      factpi = twopi
      jmax = 2*lmax(1) + 1
      jmax1 = lmax(1) - lmin(1) + 1
      lmax1 = lmax(1) + 1
      ln = lmin(1)
      lx = lmax(1)
      jdel = 8
      q = qa1
      nq = n*q
      tmth = 2*mth
      mthsq = tmth * tmth
      lmth = tmth * 2*jmax1
      j1v = nfm
      j2v = nfm
      if ( check1 )
     $     call msctimer ( outmod, "before kernels" )
      j1 = 1
      j2 = 1
      ksgn = 2*j2 - 3
      call kernel(xpla,zpla,xpla,zpla,grdgre,grwp,j1,j2,ksgn,1,1,0)
      if ( checkd .and. (lfele .eq. 0) ) then
         call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
     $        "grwp at 1,1", outmod, iotty )         
         call matwrtn ( grdgre,nths2,nths2,1,1,mth,mth,mth,mth,
     $        "grdgre at 1,1", outmod, iotty ) 
      endif
      do i = 1, mth2
         grwp(i,mth1) = grwp(i,1)
         grwp(i,mth2) = grwp(i,2)
      enddo
      if ( lfele .eq. 1 ) then
         call felang ( grwp, grri, cnqd, 0,0 )
         call felang ( grwp, grri, snqd, 0,jmax1 )
      endif
      if ( lfour .eq. 1 ) then
         call fouran ( grwp, grri, cslth, 0,0 )
         call fouran ( grwp, grri, snlth, 0,jmax1 )
      endif
      if ( checkd .and. (lfele .eq. 0) ) then
         call foura2 ( grdgre,0,0, grwp, 0 )
         call fanal1 ( grwp,nths,nths,0, wrkvr,wrkvi,nfm,nfm )
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Kpp(l,l)", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Kpp(l,l)", outmod,0 )
         call fanal1 ( grri,nths2,nfm2,0, wrkvr,wrkvi,nfm,nfm )
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real GPP_ll", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag GPP_ll", outmod,0 )
      endif
      iwp = 1
      if ( .not. farwal ) then
         call wwall ( mth1, xwal, zwal )
         j1 = 1
         j2 = 2
         ksgn = 2*j2 - 3
         call kernel(xpla,zpla,xwal,zwal,grdgre,grwp,j1,j2,ksgn,0,0,1)
         j1 = 2
         j2 = 2
         ksgn = 2*j2 - 3
         call kernel(xwal,zwal,xwal,zwal,grdgre,grwp,j1,j2,ksgn,0,0,1)
         if ( checkd .and. (lfele .eq. 0) ) then
            call foura2 ( grdgre,mth,mth, grwp, 0 )
            call fanal1 ( grwp,nths,nths,0, wrkvr,wrkvi,nfm,nfm )
            call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Real Kww(l,l)", outmod,0 )
            call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Imag Kww(l,l)", outmod,0 )
         endif
         j1 = 2
         j2 = 1
         ksgn = 2*j2 - 3
         call kernel(xwal,zwal,xpla,zpla,grdgre,grwp,j1,j2,ksgn,1,0,1)
         do i = 1, mth2
            grwp(i,mth1) = grwp(i,1)
            grwp(i,mth2) = grwp(i,2)
         enddo
         if ( lfele .eq. 1 ) then
            call felang ( grwp, grri, cnqd,  mth,0 )
            call felang ( grwp, grri, snqd,  mth,jmax1 )
         endif
         if ( lfour .eq. 1 ) then
            call fouran ( grwp, grri, cslth, mth,0 )
            call fouran ( grwp, grri, snlth, mth,jmax1 )
         endif
         if ( checkd .and. (lfele .eq. 0) ) then
            call fanal1 ( grri,nths2,nfm2,mth, wrkvr,wrkvi,nfm,nfm )
            call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Real GWP_ll", outmod,0 )
            call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Imag GWP_ll", outmod,0 )
         endif
         iwp = 2
      endif
      if ( check1 )
     $     call msctimer ( outmod, "aft kern and fourier" )
      mth12 = iwp*mth
      if ( n .le. 0.1 ) then
         do j = 1, mth12
            do ii = 1, 2
               summ(ii) = zero
               i1 = (ii-1) * mth + 1
               i2 = ii * mth
               do i = i1, i2
                  summ(ii) = summ(ii) + grdgre ( j, i )
               enddo
            enddo
            sum11 = summ(1)
            sum12 = summ(2)
            if ( check2 )
     $           write ( outmod, 9002 ) j, sum11, sum12
         enddo
      endif
      if ( check2 ) then
      endif
      if ( checkd .and. (lfele .eq. 0) ) then
         call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
     $        "grwp at end", outmod, iotty )         
         call matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
     $        "grdgre at end", outmod, iotty )   
      endif
      write ( iotty,  '(/,
     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
      write ( outmod, '(/,
     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
      do j1 = 1,2
         do j2 = 1,2
            sumg = 0.0
            do i = 1, mth
               sumg = sumg + grdgre( (j1-1)*mth+i,(j2-1)*mth+1 )
            enddo
            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
         enddo
      enddo
      write ( iotty,  '(/,
     $     "Sum of first ROW in each block of GRDGRE:",/)')
      write ( outmod, '(/,
     $     "Sum of first ROW in each block of GRDGRE:",/)')
      do j1 = 1,2
         do j2 = 1,2
            sumg = 0.0
            do i = 1, mth
               sumg = sumg + grdgre( (j1-1)*mth+1,(j2-1)*mth+i )
            enddo
            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
         enddo
      enddo
      if ( check1 )
     $     call msctimer ( outmod, "before leqt1f" )
      lmax2 = 2*jmax1
      if ( (abs(n) .le. 1.e-5) .and. .not. farwal .and. (ishape .le. 10)
     $     )then
         write (outmod,8180) cn0
         do i = 1, mth12
            do j = 1, mth12
               grdgre(i,j) = grdgre(i,j) + cn0
            enddo
         enddo
      endif
      ier = 0

      IF (kernelsign < 0) THEN
         grdgre=kernelsign*grdgre
         DO i=1,nths2
            grdgre(i,i)=grdgre(i,i)+2.0
         ENDDO
      ENDIF

      call gelimb (grdgre,nths2,dummy,nths2,mth12,lmax2,
     $        grri,nths2,work1,ier )
      deallocate(grdgre)

      allocate(arr(nfm,nfm),ari(nfm,nfm),air(nfm,nfm),aii(nfm,nfm))
      if ( check1 )
     $     call msctimer ( outmod, "after leqt1f" )
      write ( outmod,8050 ) ier
      write ( iotty, 8050 ) ier
      if ( lspark .ne. 0 ) then
         write ( iodsk, 8312 ) 
         write ( iodsk, 8311 ) jmax1, mth12
         do jwdsk = 1, lmax2
            write ( iodsk, 8313 ) ( grri(iwdsk,jwdsk),iwdsk=1,mth12 )
         enddo
      endif
      if ( check2 ) write (outmod,215)
      do l1 = 1, jmax1
         sumpc = 0.0
         sumps = 0.0
         sumwc = 0.0
         sumws = 0.0
         do i = 1, mth
            sumpc = sumpc + grri(i,l1)
            sumps = sumps + grri(i,jmax1+l1)
            sumwc = sumwc + grri(mth+i,l1)
            sumws = sumws + grri(mth+i,jmax1+l1)
         enddo
         ll = l1 - 1 + lmin(1)
         if ( check2 ) write ( outmod,217 ) ll, sumpc,sumps,
     $        sumwc,sumws
      enddo
      call wwall ( mth1, xwal, zwal )
      mw = mth - (jtop-jbot-1)
      do l1 = 1, jmax1
         iw = mth - jtop + 2
         icnt = 0
         do i = 1, mth
            il11 = (l1-1)*mth12 + mth + i
            il12 = jmax1*mth12 + il11
            if (  (i.le.jbot) .or. (i.ge.jtop) ) then
               if ( i .eq. jtop ) iw = 1
               chiwc(iw,l1) = grri(mth+i,l1)
               chiws(iw,l1) = grri(mth+i,jmax1+l1)
               xpass(iw) = xwal(i)
               zpass(iw) = zwal(i)
               iw = iw + 1
               icnt = icnt + 1
            endif
         enddo
         if ( jtop .eq. jbot ) then
            chiwc(mw,l1) = chiwc(1,l1)
            chiws(mw,l1) = chiws(1,l1)
            xpass(mw) = xpass(1)
            zpass(mw) = zpass(1)
         endif
      enddo
      if ( jtop .eq. jbot ) then
         if ( checkd .and. (lfele .eq. 0) ) then
            call fanal1 ( grri,nths2,nfm2,0, wrkvr,wrkvi,nfm,nfm )
            call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Chiplar", outmod,0 )
            call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Chiplai", outmod,0 )
            call fanal1 ( grri,nths2,nfm2,mth, wrkvr,wrkvi,nfm,nfm )
            call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Chiwal", outmod,0 )
            call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "Chiwali", outmod,0 )
         endif
      endif
      if ( check2 )
     $     write ( outmod,8020 ) jtop, jbot, mw, icnt
      if ( lfele .ne.  0 ) then
         write ( outmod, '(/,20x,"Finite Elements",/)' )
         call felanv ( grri,arr, cnqd, 0,0 )
         call felanv ( grri,aii, snqd, 0,jmax1 )
         call felanv ( grri,ari, snqd, 0,0 )
         call felanv ( grri,air, cnqd, 0,jmax1 )
         do j1 = 1, jmax1
            do j2 = 1, jmax1
               vacmat(j1,j2) = ( arr(j1,j2) + aii(j1,j2) )
               vacmti(j1,j2) = ( air(j1,j2) - ari(j1,j2) )
            enddo
         enddo
      else
         call foranv ( grri,arr, cslth, 0,0 )
         call foranv ( grri,aii, snlth, 0,jmax1 )
         call foranv ( grri,ari, snlth, 0,0 )
         call foranv ( grri,air, cslth, 0,jmax1 )
         if ( lnova ) then
            do l1 = 1, jmax1
               ll1 = l1 - 1 + lmin(1)
               al1nq = ll1 - nq
               do l2 = 1, jmax1
                  ll = l2 - 1 + lmin(1)
                  al2nq = ll - nq
                  do i = 1, mth
                     theta = (i-1)*dth
                     ar = ( xjdtxj(i)*grri(i,l1)-
     $                    al2nq*grri(i,jmax1+l1) ) * factpi
                     ai = ( xjdtxj(i)*grri(i,jmax1+l1)+
     $                    al2nq*grri(i,l1) ) * factpi
                     ar1 = (ai*cslth(i,l2)-ar*snlth(i,l2))
     $                    /xjacob(i)
                     ai1 = -(ar*cslth(i,l2)+ai*snlth(i,l2))
     $                    /xjacob(i)
                     rmatr(l2,l1) = rmatr(l2,l1)
     $                    + ar1
                     rmati(l2,l1) = rmati(l2,l1)
     $                    + ai1
                     el1l2t = ( ll1-ll ) * theta
                     ajll(l2,l1) = ajll(l2,l1)
     $                    + cos(el1l2t)
     $                    / (twopi*xjacob(i))
                  enddo
                  ajll(l2,l1) = dth * ajll(l2,l1)
                  rmatr(l2,l1) = rmatr(l2,l1)
     $                 * dth*al1nq/twopi2
                  rmati(l2,l1) = rmati(l2,l1)
     $                 * dth*al1nq/twopi2
               enddo
            enddo
         endif
         do j1 = 1, jmax1
            do j2 = 1, jmax1
               vacmat(j1,j2) = arr(j1,j2) + aii(j1,j2)
               vacmti(j1,j2) = air(j1,j2) - ari(j1,j2) 
            enddo
         enddo
      endif
      if ( checkd ) then
         call matwrtn ( arr, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "arr", outmod,0 )
         call matwrtn ( aii, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "aii", outmod,0 )
         call matwrtn ( ari, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "ari", outmod,0 )
         call matwrtn ( air, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "air", outmod,0 )
      endif
      write ( iotty,  '(/,1x," n, mth, lmin, lmax = ",f4.1, 3i5)')
     $     n, mth, lmin(1), lmax(1)
      write ( outmod, '(/,1x," n, mth, lmin, lmax = ",f4.1, 3i5)')
     $     n, mth, lmin(1), lmax(1)
      if ( lgato .eq. 1 ) then
         call gatonorm ( vacmat, gatovac, nfm, rgato,mfel,mth,
     $        qa1,twopi )
         call matwrtn ( gatovac,nfm,nfm,ln,ln,jmax1,jmax1,8,8,
     $        "GATOVAC", outmod,iotty )      
      endif
      if ( check1 )
     $     call msctimer ( outmod, "end of vacmat" )
      if ( lgato .eq. 2 )  then 
         call orchek ( air, ari, rmatr, rmati, work, work1 )
         do j1 = 1, jmax1
            do j2 = 1, jmax1
               rmatr(j1,j2) = (arr(j1,j2) + aii(j1,j2))
            enddo
         enddo
         iopc = 1
         iops = 2
         call fotofi ( rmatr, rmati, sinlt, air, ari, iops )
         call matwrtn ( rmati,nfm,nfm,1,1,mfel,mfel,8,8,"VFinels",
     $        outmod,iotty )
         call fotofi ( rmatr, rmatr, coslt, air, ari, iopc )
         call matwrtn ( rmatr,nfm,nfm,1,1,mfel,mfel,8,8,"VFinelc",
     $        outmod,iotty )
         do ma = 1, mfel
            do mb = 1, mfel
               rmatr(ma,mb) = rmatr(ma,mb) + rmati(ma,mb)
            enddo
         enddo
         call matwrtn ( rmatr,nfm,nfm,1,1,mfel,mfel,8,8,"VFinel",
     $        outmod,iotty )
         call vacasym ( rmatr, nfm,mfel,"VFinel", outmod,iotty )
         call gatonorm ( rmatr, gatovac, nfm, rgato,mfel,mth,
     $        qa1,twopi )
         call matwrtn ( gatovac,nfm,nfm,1,1,mfel,mfel,8,8,
     $        "GATOVAC", outmod,iotty )      
      endif
      deallocate(arr,ari,air,aii)
      write ( outmod, 500 ) n,q, nj,mj,lj
      if ( (.not. lpest1) .and. check2 )
     $     write ( outmod, 501 )  ( delta(i),i=1,mth1 )
      if ( checke ) then
         if ( lgato .eq. 2 ) then
            write ( outmod, '(/,20x, "lgato = 2:")' )
            call mateig2 ( rmatr, nfm,mfel, wrkvr,work0,work,work1,
     $           ln,lx,ff, 0, jobid, outmod )
         endif
         if ( lgato .eq. 1 ) then
            write ( outmod, '(/,20x, "lgato = 1:")' )
            call mateig2 ( vacmat, nfm,mfel, wrkvr,work0,work,work1,
     $           ln,lx,ff, 0, jobid, outmod )
         endif
      endif
      if ( lgato .eq. 0 ) then
         do l1 = 1, jmax1
            do l2 = 1, jmax1
               ajll(l1,l2) = ( rmati(l1,l2)-rmatr(l1,l2) )
            enddo
         enddo
         do j1 = 1, jmax1
            l1 = j1 + lmin(1) - 1
            alnq1 = l1 - nq
            do j2 = 1, jmax1
               l2 = j2 + lmin(1) - 1
               alnq2 = l2 - nq

               vacmatu(j1,j2) = vacmat(j1,j2)
               vacmtiu(j1,j2) = vacmti(j1,j2)
            enddo
         enddo
         if ( lnova ) then
            call matwrtn(rmati,mtot,mtot,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "rmati", outmod,0 )
            call matmul1 ( ajll,vacmat, mtot,nfm,jmax1, rmati,mtot )
            do l1 = 1, jmax1
               do l2 = 1, jmax1
                  rmati(l1,l2) = rmati(l1,l2) / twopi2
               enddo
            enddo
            call matwrtn ( ajll,mtot,mtot,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "inverse j matrix", outmod,0 )
            call matwrtn ( rmati,mtot,mtot,ln,ln,jmax1,jmax1,jdel,jdel,
     $           "p-xi matrix ", outmod,0 )
         endif
         write ( outmod, '(a)' ) char(12)
         call matwrtn ( vacmat,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Vacmat with l-nq", outmod,iotty )
         call matwrt9 ( vacmat, nfm,nfm,ln,lx, "Vacmat with l-nq",
     $        outmod,iotty )
         if ( checkd )
     $        call matwrtn(vacmti,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "vacmti", outmod,0 )
         call vacasym ( vacmat, nfm,jmax1,"Vacmat", outmod,iotty )
         call masym0 ( vacmat, nfm,jmax1,ln,ln,0,0,
     $        "Vacmat", outmod,iotty )
         if ( checke ) then
            write ( outmod, '(/,20x, "Fourier With l-nq :")' )
            call mateig2 ( vacmat, nfm,jmax1, wrkvr,work0,work,work1,
     $           ln,lx, ff, 2, jobid, outmod )
         endif
      endif
      IF ( lsymz ) THEN

         do 601 l1 = 1, jmax1
            do 600 l2 = l1, jmax1
               vacmat(l1,l2) = 0.5 * ( vacmat(l1,l2)+vacmat(l2,l1) )
               vacmti(l1,l2) = 0.5 * ( vacmti(l1,l2)-vacmti(l2,l1) )
               rmatr(l1,l2) = 0.5 * ( rmatr(l1,l2)+rmatr(l2,l1) )
 600        continue
 601     continue
         do 621 l1 = 1, jmax1
            do 620 l2 = l1, jmax1
               vacmat(l2,l1) = vacmat(l1,l2)
               IF ( l1 /= l2 ) vacmti(l2,l1) = - vacmti(l1,l2)
               IF ( L1 == L2 ) vacmti(l2,l1) = 0.0
               rmatr(l2,l1) = rmatr(l1,l2)
 620        continue
 621     continue
      ENDIF
      if ( ladj .eq. 1 ) then
         open ( 60, file='vacadj', status='unknown', form='formatted'
     .        ,recl=160)
         nadj = n + 0.001
         write ( 60, 554 ) mthin1,lmin(1),lmax(1),nadj,qa1
         do j1 = 1, jmax1
            write ( 60, 555 )  ( vacmat(j1,j2), j2=1,jmax1 )
         enddo
         close ( 60 )
      elseif ( (ldcon .eq. 1).or.(lgpec .eq. 1)) then
         open ( 62, file='vacdcon', status='unknown', form='formatted'
     .        ,recl=160)
         ndcon = n + 0.001
         write ( 62,'("mthin1, lmin, lmax, n, q = ")' )
         write ( 62, '(4i5, e15.7)' ) mthin1,lmin(1),lmax(1),ndcon,qa1
         write ( 62, '("Real Vacmat:")' )
         do j1 = 1, jmax1
            write ( 62,'(10e15.7)' )  ( vacmat(j1,j2), j2=1,jmax1 )
         enddo
         write ( 62, '("Imag. Vacmat:")' )
         do j1 = 1, jmax1
            write ( 62,'(10e15.7)' )  ( vacmti(j1,j2), j2=1,jmax1 )
         enddo
         close ( 62 )
      elseif ( lnova ) then
         call matwrtn(rmatr,mtot,mtot,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "rmatr", outmod,0 )
         call matwrtn ( ajll,mtot,mtot,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "pll-mc minus pll-fc ", outmod,0 )
         mtots = mtot**2
         len=mtots+5
         ndsk=1
         call zop (iovac,"vacout",len,ndsk,iiff,999)
         lgivup=1
         call zwr(iovac,rmatr,mtots,1,lgivup,999)
         call zcl ( iovac, 999 ) 
      else
         j12 = 1
         do j2 = 1, jmax1
            do j1 = 1, jmax1
               vacpstr(j12) = vacmat(j1,j2)
               vacpsti(j12) = vacmti(j1,j2)
               j12 = j12 + 1
            enddo
         enddo
         call wtopest ( vacpstr, vacpsti, xwal, zwal )
      endif
      call arrays
      if ( check2 ) then
         write ( outmod, '("I, xp, zp, xw, zw, xpp, zpp, xwp, zwp =")' )
         write ( iotty,  '("I, xp, zp, xw, zw, xpp, zpp, xwp, zwp =")' )
         do i = 1, mth1, 8
            write ( outmod, '(i3, 1x, 8f8.2)' )
     $           i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
            write ( iotty,  '(i3, 1x, 8f8.2)' )
     $           i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
         enddo
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
 999  call errmes ( outpest,'vacuum' )
      END SUBROUTINE vaccal
c-----------------------------------------------------------------------
c     subprogram 2. kernel.
c     computes kernels of integral equation for Laplace's
c     equation for a torus.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE kernel(xobs,zobs,xsce,zsce,grdgre,gren,
     $     j1,j2,isgn,iopw,iops,ischk)
      USE vglobal_mod
      USE local_mod, ONLY: r8
      IMPLICIT NONE

      REAL(r8), DIMENSION(:), INTENT(IN) :: xobs,zobs,xsce,zsce
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: grdgre
      REAL(r8), DIMENSION(:,:), INTENT(OUT), TARGET :: gren
      INTEGER, INTENT(IN) :: j1,j2,isgn,iopw,iops,ischk

      INTEGER  :: i,ic,iend,ig,ilr,isph,istart,j,j1j2,jres,js1,js2,
     $     js3,js4,js5,mthm,mths
      INTEGER, DIMENSION(2) :: iop
      INTEGER :: sTime, fTime, cr

      REAL(r8) :: a0,ak0i,alg,alg0,alg1,alg2,algdth,am,amm,ap,app,
     $     aval1,resdg,residu,resk0,slog0,slog1m,slog1p,tdth,
     $     thes,theta,third,wgbg,wsimpa,wsimpa1,wsimpa2,wsimpa4,wsimpb,
     $     wsimpb1,wsimpb2,wsimpb4,xl,xu
      REAL(r8), DIMENSION(3) :: tab
      REAL(r8), DIMENSION(nths) :: the,xpp,zpp,work,ww1,ww2,ww3,xpr,zpr
c-----------------------------------------------------------------------
c     declarations of Gaussian quadrature quantities.
c-----------------------------------------------------------------------
      REAL(r8) :: tgaus0,agaus,bgaus,pgaus,pgaus2
      REAL(r8), DIMENSION(8) :: tgaus
      REAL(r8), DIMENSION(8), PARAMETER :: wgaus=(/
     $     0.101228536290376_8,0.222381034453374_8,
     $     0.313706645877887_8,0.362683783378362_8,
     $     0.362683783378362_8,0.313706645877887_8,
     $     0.222381034453374_8,0.101228536290376_8/)
      REAL(r8), DIMENSION(8), PARAMETER :: xgaus=(/
     $     -0.960289856497536_8,-0.796666477413627_8,
     $     -0.525532409916329_8,-0.183434642495650_8,
     $     0.183434642495650_8,0.525532409916329_8,
     $     0.796666477413627_8,0.960289856497536_8/)
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      DO i=1,nths
         xpp(i)=0
         zpp(i)=0
         ww1(i)=0
         ww2(i)=0
         ww3(i)=0
      ENDDO
      tab(1:3)=0
      ak0i=0
      jres=1
      isph = 0
      mthm=mth-1
      gren(1:mth1,1:mth1)=0
      the(1:mth1)=(/(i,i=0,mth)/)*dth
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(lfele /= 0)THEN
         wsimpb1=dth/two
         wsimpb2=dth
         wsimpb4=dth 
      ELSE
         wsimpb1=dth/three
         wsimpb2=two*dth/three
         wsimpb4=four*dth/three
      ENDIF

      wsimpa1=dth/three
      wsimpa2=two*dth/three
      wsimpa4=four*dth/three
      third=one/three
      algdth=LOG(dth)
      slog1m=third*dth*(algdth-third)
      slog0=four*third*dth*(algdth-four*third)
      slog1p=slog1m
      tdth=two*dth
      alg=LOG(tdth)
      alg0=16.0*dth*(alg-68.0/15.0)/15.0
      alg1=128.0*dth*(alg-8.0/15.0)/45.0
      alg2=4.0*dth*(7.0*alg-11.0/15.0)/45.0

      IF(ischk == 0)THEN
         jbot=mth/2+1
         jtop=mth/2+1
         CALL wwall(mth1,ww1,ww2)
         DO i=1,mth1
            IF(i == mth1)CYCLE
            IF(ww1(i)*ww1(i+1) > zero)CYCLE
            IF(ww1(i) > zero)jbot=i
            IF(ww1(i) < zero)jtop=i+1
            isph=1
         ENDDO
      ENDIF

      iop(1)=4
      iop(2)=4
      CALL spl1d1(mth1,the,xsce,xpp,iop,1,ww1,ww2,ww3)
      CALL spl1d1(mth1,the,zsce,zpp,iop,1,ww1,ww2,ww3)

      DO i=1,mth1
         theta=(i-1)*dth
         CALL spl1d2(mth1,the,xsce,xpp,1,theta,tab)
         xpr(i)=tab(2)
         CALL spl1d2(mth1,the,zsce,zpp,1,theta,tab)
         zpr(i)=tab(2)
      ENDDO
c-----------------------------------------------------------------------
c     begin loop over intervals.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT_RATE=cr)
      CALL SYSTEM_CLOCK(COUNT=sTime)
      DO j=1,mth
         xs=xobs(j)
         zs=zobs(j)
         thes=the(j)
         work(1:mth1)=0

         IF(xs < zero)GO TO 175

         aval1=zero
         iend=2

         IF(isph == 1 .AND. j2 == 2)THEN
            IF(jbot-j == 1)iend=3
            IF(jbot-j == 0)iend=4
            IF(j-jtop == 0)iend=0
            IF(j-jtop == 1)iend=1
         ENDIF

         istart=4-iend
         mths=mth-(istart+iend-1)

         DO i=1,mths
            ic=j+i+istart-1
            IF(ic >= mth1)ic=ic-mth
            theta=(ic-1)*dth
            xt=xsce(ic)
            zt=zsce(ic)
            IF(xt < zero)CYCLE
            xtp=xpr(ic)
            ztp=zpr(ic)
            IF(ic == j)CYCLE
            CALL green
            wsimpb=wsimpb2
            IF((i/2)*2 == i)wsimpb=wsimpb4
            IF(i == 1 .OR. i == mths)wsimpb=wsimpb1
            wsimpa=wsimpa2
            IF((i/2)*2 == i)wsimpa=wsimpa4
            IF(i == 1 .OR. i == mths)wsimpa=wsimpa1
            work(ic)=work(ic)+isgn*aval*wsimpa
            gren(j,ic)=gren(j,ic)+bval*wsimpb
            aval1=aval1+aval0*wsimpa
         ENDDO

         j1j2=j1+j2

         IF(j1j2 /= 2 .AND. isph == 1 .AND. j > jbot .AND. j < jtop)
     $        go to 175

         thes=the(j)
         js1=MOD(j-iend+mth-1,mth)+1
         js2=MOD(j-iend+mth,mth)+1
         js3=MOD(j-iend+mth+1,mth)+1
         js4=MOD(j-iend+mth+2,mth)+1
         js5=MOD(j-iend+mth+3,mth)+1
c-----------------------------------------------------------------------
c     perform gaussian quadratures.
c-----------------------------------------------------------------------
         DO ilr=1,2
            xl=thes+(2*ilr-iend-2)*dth
            xu=xl+tdth
            agaus=half*(xu+xl)
            bgaus=half*(xu-xl)
            tgaus=agaus+xgaus*bgaus
            DO ig=1,8
               tgaus0=tgaus(ig)
               IF(tgaus0 < zero)tgaus0=twopi+tgaus0
               IF(tgaus0 >= twopi)tgaus0=tgaus0-twopi
               CALL spl1d2(mth1,the,xsce,xpp,1,tgaus0,tab)
               xt=tab(1)
               xtp=tab(2)
               CALL spl1d2(mth1,the,zsce,zpp,1,tgaus0,tab)
               zt=tab(1)
               ztp=tab(2)
               CALL green
               bval=bval+iops*LOG((thes-tgaus(ig))**2)/xs
               pgaus=(tgaus(ig)-thes-(2-iend)*dth)/dth
               pgaus2=pgaus*pgaus
               wgbg=wgaus(ig)*bgaus
               amm=(pgaus2-one)*pgaus*(pgaus-two)/24.0
               amm=amm*wgbg
               work(js1)=work(js1)+isgn*aval*amm
               am=-(pgaus-one)*pgaus*(pgaus2-four)/6.0
               am=am*wgbg
               work(js2)=work(js2)+isgn*aval*am
               a0=(pgaus2-one)*(pgaus2-four)/four
               work(js3)=work(js3)+isgn*a0*aval*wgbg
               ap=-(pgaus+one)*pgaus*(pgaus2-four)/6.0
               ap=ap*wgbg
               work(js4)=work(js4)+isgn*aval*ap
               app=(pgaus2-one)*pgaus*(pgaus+two)/24.0
               app=app*wgbg
               work(js5)=work(js5)+isgn*aval*app
               work(j)=work(j)-isgn*aval0*wgbg
               IF(j == jres)ak0i=ak0i-isgn*aval0*wgbg
               IF(iopw == 0)CYCLE
               gren(j,js1)=gren(j,js1)+bval*amm
               gren(j,js2)=gren(j,js2)+bval*am
               gren(j,js3)=gren(j,js3)+bval*a0*wgbg
               gren(j,js4)=gren(j,js4)+bval*ap
               gren(j,js5)=gren(j,js5)+bval*app
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     remaining computations.
c-----------------------------------------------------------------------
         residu=0.0
         IF(j1 == j2)residu=two

         IF(ishape < 10 .OR. ishape == 41 .OR. ishape == 42)THEN
            resdg=(2-j1)*(2-j2)+(j1-1)*(j2-1)
            resk0=(2-j1)*(2-j2)+(j1-3)*(j2-1)
            residu=resdg+resk0
         ENDIF

         work(j)=work(j)-isgn*aval1+residu
         IF(j == jres)ak0i=ak0i-isgn*aval1

         IF(iops == 1 .AND. iopw /= 0)THEN
            gren(j,js1)=gren(j,js1)-alg2/xs
            gren(j,js2)=gren(j,js2)-alg1/xs
            gren(j,js3)=gren(j,js3)-alg0/xs
            gren(j,js4)=gren(j,js4)-alg1/xs
            gren(j,js5)=gren(j,js5)-alg2/xs
         ENDIF

 175     CONTINUE

         IF((xs < zero).AND.(j2 == 2))work(j)=1.0

         DO ic=1,mth
            grdgre((j1-1)*mth+j,(j2-1)*mth+ic)=work(ic)
            gren(j,ic)=gren(j,ic)/twopi
         ENDDO
      ENDDO
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_timer_output) THEN
         print *,"vac loop time=",REAL(fTime-sTime,r8)/REAL(cr,r8)
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE kernel
c-----------------------------------------------------------------------
c     subprogram 3. mateig.
c     computes matrix eigenvalues.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mateig ( zvec, work, work1 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension zvec(nfm,nfm), work(*), work1(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 9003 format(1x,/,1x, "eigenvector no.",i4,
     $     "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      do i = 1, nfm*nfm
         work1(i) = zero
      enddo
      kk = 1
      do i = 1, jmax1
         do j = 1, i
            work1(kk) = zvec(i,j)
            kk = kk + 1
         enddo
      enddo
      call eigen ( work1,work, jmax1, 0 )
      kk = 1
      do i = 1, jmax1
         ii = ( i*  (i+1) ) / 2
         j = (i-1)*jmax1 + 1
         jj = j + jmax1 - 1
         write(outmod,9003) i, work1(ii), ( work(j1), j1=j,jj )
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine mateig
c-----------------------------------------------------------------------
c     subprogram 4. mateig2.
c     computes matrix eigenvalues.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mateig2 ( zvec, nd,msiz, zwk,work0,work,work1, l1,l2,
     $     ff1,lcone, jobid1, nout1 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      character(*) jobid1
      dimension zvec(nd,nd), zwk(nd,nd), work(*), work1(*), work0(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 9003 format(1x,/,1x, "eigenvector no.",i4,
     $     "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mfsq = msiz*msiz
      do i = 1, mfsq
         work1(i) = 0
      enddo
      kk = 1
      do i = 1, msiz
         do j = 1, i
            work1(kk) = zvec(i,j)
            kk = kk + 1
         enddo
      enddo
      call eigen ( work1,work, msiz, 0 )
      kk = 1
      do i = 1, msiz
         ii = ( i*  (i+1) ) / 2
         j = (i-1)*msiz + 1
         jj = j + msiz - 1
         work0(i) = work1(ii)
         write(nout1,9003) i, work1(ii), ( work(j1), j1=j,jj )
      enddo
      do i = 1, msiz
         ii = ( i*(i+1) ) / 2
         j = (i-1)*msiz + 1
         jj = j + msiz - 1
         kk = 1
         do j1 = j, jj
            zwk(kk,i) = work(j1)
            kk = kk + 1
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine mateig2
c-----------------------------------------------------------------------
c     subprogram 5. arrays.
c     computes arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine arrays
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      real nq
      dimension the(nths)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 3    format ( 1x, "  delx, dely = ", 1p2e13.5 )
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      mth12 = 2*mth
      jmax1 = lmax(1) -lmin(1) + 1
      q = qa1
      nq = n*q
      call bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )
      delx = plrad * delfac
      delz = plrad * delfac
      write ( iotty, 3 ) delx, delz
      write ( outmod, 3 ) delx, delz
      do i = 1, mth1
         the(i) = (i-1) * dth
      enddo
      call wwall ( mth1, xwal, zwal )
      call difspl ( mth, the, xwal, xwalp )
      call difspl ( mth, the, zwal, zwalp )
      xwalp(mth+1) = xwalp(1)
      zwalp(mth+1) = zwalp(1)
c-----------------------------------------------------------------------
c     big do loop.
c-----------------------------------------------------------------------
      do l1 = 1, jmax1
         do i = 1, mth
            cplar(i,l1) = grri(i,l1)
            cplai(i,l1) = grri(i,jmax1+l1)
            if ( .not. farwal ) then
               cwallr(i,l1) = grri(mth+i,l1)
               cwalli(i,l1) = grri(mth+i,jmax1+l1)
            endif
         enddo
         cplar(mth1,l1) = cplar(1,l1)
         cplai(mth1,l1) = cplai(1,l1)
         if ( .not. farwal ) then
            cwallr(mth1,l1) = cwallr(1,l1)
            cwalli(mth1,l1) = cwalli(1,l1)
         endif
      enddo
c-----------------------------------------------------------------------
c     another big do loop.
c-----------------------------------------------------------------------
      do is = 1, mth1
         theta = (is-1) * dth
         znqd = nq*delta(is)
         cnqd(is) = cos(znqd)
         snqd(is) = sin(znqd)
         do l1 =  1, jmax1
            ll = lmin(1) - 1 + l1
            elth = ll * theta
            elthnq = ll * theta + znqd
            sinlt(is,l1) = sin(elth)
            coslt(is,l1) = cos(elth)
            snlth(is,l1) = sin(elthnq)
            cslth(is,l1) = cos(elthnq)
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine arrays
c-----------------------------------------------------------------------
c     subprogram 6. wwall.
c     conducting wall specification.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine wwall(nqnqnq,xwal1,zwal1)
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      INTEGER:: npots0,npots
      logical lfix, insect
      REAL(r8) :: csmin
      REAL(r8), DIMENSION(:),POINTER :: thetatmp,xwaltmp,xpptmp,
     $            ww1tmp,ww2tmp,ww3tmp,tabtmp,zwaltmp
      INTEGER, DIMENSION(:),POINTER :: rioptmp
      dimension xwal1(*), zwal1(*)
      dimension iop(2),xpp(nths),zpp(nths),ww1(nths),ww2(nths),
     $    ww3(nths),thet(nths),tabx(3),tabz(3)
c      character*(80) stringv(5)
      data iplt/0/
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 1000 format(" error bw is less than 1 =",e12.5," wall will intersect
     $  the plasma"/)
 1100 format(" **error** b is less than zero=",e12.5)
 1430 format(" newton scheme for finding the nearest wall point",/,
     $  "  did not converge ",/,
     $  " i=",i4,"x/zplasma=",2e15.8,"x/zwal1l=",2e15.8,
     $ " theta=",1e10.3," f=",1e15.8,"fp=",1e15.8)
 1460 format(i3," xinf,zinf=",2e10.3," xwal1,zwal1=",2e10.3)
 1450 format(" there are at least ",i3," wall points in the plasma")
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      xpp(1) = 1.
      zpp(1) = 1.
      ww3(1) = 1.
      tabx(1) = 1.
      tabz(1) = 1.
      awsave = aw
      bwsave = bw
      insect = .false.
      isph = 0
c-----------------------------------------------------------------------
c     ishape< 0 Spherical topology
c     ishape<10 Closed toroidal topology
c     ishape<20 Solid conductors not linking plasma
c     ishape<30 Toroidal conductor with a toroidally symmetric gap.
c               Geometry correlated to plasma position and geometry.
c     ishape<40 Toroidal conductor with a toroidally symmetric gap.
c               Geometry independent of plasma.
c-----------------------------------------------------------------------
      if ( a < -100.0 ) then
        isph = 1
        ishape = -10
        call bounds(xinf,zinf,1,mth,xmnp,xmxp,zmnp,zmxp)
        xmin = xmnp
        xmax = xmxp
        zmin = zmnp
        zmax = zmxp
        plrad = 0.5 * ( xmxp - xmnp )
        xmaj = 0.5 * ( xmxp + xmnp )
        zmid = 0.5 * ( zmnp + zmxp )
        hrad = xmax + aw*(xmax-xmaj)
        vrad = zmax + bw*(zmax-zmid)
        do i = 1, mth1
          xi = xinf(i) - xmaj
          zeta = zinf(i) - zmid
          bbb = (xi*vrad)**2 + (zeta*hrad)**2
          ccc = -xmaj*vrad*xi + hrad*sqrt( bbb - (zeta*xmaj)**2 )
          xwal1(i) = xmaj + xi*vrad*ccc/bbb
          zwal1(i) = zmid + zeta*vrad*ccc/bbb
        enddo
        go to 145 !go directly to cleanup
      endif
      if(a .gt. -10.)lfix=.true.
      if( farwal ) then
        return !early exit because there is no wall
      endif
      xshift = a

      mthalf = mth2 / 2
      zmin = 3.4e38
      zmax = -3.4e38
      xmin = 3.4e38
      xmax = -3.4e38
      do i = 1, mth
        if(xmax .lt. xinf(i))xmax = xinf(i)
        if(xmin .gt. xinf(i))xmin = xinf(i)
        if(zmax .lt. zinf(i)) zmax = zinf(i)
        if(zmin .gt. zinf(i))zmin = zinf(i)
      enddo
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      zmid = 0.5 * ( zmax + zmin )
      zrad = 0.5 * ( zmax - zmin )
      scale =  ( zmax - zmin )
      if(((xmax-xmin)/2.0).gt.scale) scale = (xmax-xmin)/2.0
      scale = 1.0
      aw = aw * scale
      bw = bw * scale
      delta1 = dw * (xinf(1) - xma)
c-----------------------------------------------------------------------
c     ishape=2 Elliptical shell confocal to the plasma's radius and
c              height. The radius of the shell is a.
c-----------------------------------------------------------------------
      if ( ishape==2 ) then
        zh = sqrt ( abs(zrad**2 - plrad**2) ) ! h metric
        zah = a / zh
        zph = plrad / zh
        zmup = 0.5*dlog ((zrad+plrad)/(zrad-plrad)) ! mu-plas
        zmuw = dlog( zah + sqrt(zah**2 + 1) ) ! mu-wall
        zxmup = exp(zmup)
        zxmuw = exp(zmuw)
        zbwal = zh * cosh ( zmuw ) ! Major radius of wall
        bw = zbwal / a          ! Elongation of wall
        do i = 1, mth2
          the = (i-1) * dth
          xwal1(i) =   xmaj + a*cos( the )
          zwal1(i) = - bw * a*sin( the )
        enddo
        write ( outmod, '(/,"Confocal Ellipse:"/,
     $     "mup, expmup = ", 1p2e13.5,/,
     $     "muw, expmuw = ", 1p2e13.5,/)') zmup, zxmup, zmuw, zxmuw
        write ( iotty,  '(/,"Confocal Ellipse:"/,
     $     "mup, expmup = ", 1p2e13.5,/,
     $     "muw, expmuw = ", 1p2e13.5,/)') zmup, zxmup, zmuw, zxmuw
      endif
c-----------------------------------------------------------------------
c     ishape=3
c-----------------------------------------------------------------------
      if ( ishape==3 ) then
        do i = 1, mth2
          rr = (xinf(i)-xma)**2 + (zinf(i)-zma)**2
          ro = sqrt(rr)
          the = atan((zinf(i)-zma)/(xinf(i)-xma))
          thex = abs(the)
          lsgn = 1
          if(xma .gt. xinf(i)) the = the + pye
          if(i <= mthalf) then
            if(xma .gt. xinf(i)) thex = pye - thex
            thet(i) =abs(thex)
          endif
          if(.not. lfix) then
            ro = ro + delta1
            xwal1(i) = xma + lsgn * ro * cos ( the )
            zwal1(i) = zma + lsgn * ro * sin ( the )
          else
            xshift = ( xmax + xmin ) / 2.0
            xshift = a
            the = (i-1) * dth
            xwal1(i) = xshift + aw * cos(the + dw*sin(the) )
            zwal1(i) = zma - bw * sin(the)
          endif
          if(i .gt. mthalf) cycle
          if(zwal1(i) .lt. zmin) cycle
          j = i
          insect = .false.
          jsmall = j
          jlarge = j
          ics = 1
          if(xma .ge. xinf(i))ics = -1
          do
            if(zinf(j) .ge. zwal1(i)) then
              jsmall = j
            else
              exit
            endif
            if(j.ge. mthalf) cycle
            if(j.lt.1) cycle
            j = j + ics
          enddo
          jlarge = j
          if(abs(xinf(jsmall)-xma).ge.abs(xwal1(i)-xma))insect=.true.
          if(abs(xinf(jlarge)-xma).ge.abs(xwal1(i)-xma))insect=.true.
          if(.not. insect) cycle
          inside = inside + 1
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=4 Modified dee-shaped wall independent of plasma gerometry
c              with triangularity, dw, squareness and 2nd harmonic of zwall
c              in aw and tw. Centered at cw, radius a and elongation bw.
c              dee-shaped wall scaled to the radius and geometric center of the
c-----------------------------------------------------------------------
      if ( ishape==4 ) then
        wcentr = cw
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          sn2th = sin(2.0*the)
          xwal1(i) =   cw + a*cos( the + dw*sin(the) )
          zwal1(i) = - bw * a*sin( the + tw*sn2th ) - aw*sn2th
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=5 Dee-shaped wall scaled to the radius and geometric center of the
c              plasma. Offset of cw. Other variables as ishape=4.
c-----------------------------------------------------------------------
      if ( ishape==5 ) then
        wcentr = xmaj + cw*plrad
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          sn2th = sin(2.0*the)
          xwal1(i) = xmaj + cw*plrad +
     $         plrad*(1.0+a-cw)*cos(the+dw*sin(the))
          zwal1(i) = - bw*plrad*(1.0+a-cw) * sin( the + tw*sn2th ) -
     $         aw*plrad*sn2th
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=6 Conforming shell at distance a*prad. This is best option for a
c              close fitting shell since the plasma and shell nodes are aligned.
c-----------------------------------------------------------------------
      if ( ishape==6 ) then
        wcentr = xmaj
        ! force a 10cm center stack or one tenth of the plasma R if it is a super small tokamak
        csmin = min(0.1, 1e-1 * minval(xinf(2:mth1)))
        do i = 2, mth1
          alph = atan2 ( xinf(i+1)-xinf(i-1), zinf(i-1)-zinf(i+1) )
          xwal1(i) = xinf(i) + a*plrad * cos(alph)
          zwal1(i) = zinf(i) + a*plrad * sin(alph)
          ! if the wall crosses the R=0 axis, force a thin center stack
          if ( xwal1(i)<=csmin ) then
            xwal1(i) = csmin
          endif
        enddo
        xwal1(1) = xwal1(mth1)
        zwal1(1) = zwal1(mth1)
        xwal1(mth2) = xwal1(2)
        zwal1(mth2) = zwal1(2)
      endif
c-----------------------------------------------------------------------
c     ishape=7 Enclosing bean-shaped wall. Consult paper.
c              [Physics of Plasmas, 4, June 1997, 2161].
c-----------------------------------------------------------------------
      if ( ishape==7 ) then
        cwr = cw * pye / 180.0
        do i = 1, mth2
          the0 = (i-1)*dth
          the = the0
          rho = aw * ( 1.0 + bw*cos(the) )
          the2 = cwr * sin(the)
          xofsw = xmax + a*plrad - aw*(1.0+bw)
          xwal1(i) = xofsw + rho*cos(the2)
          zwal1(i) = - b*rho*sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=8 Wall of DIII-D.
c-----------------------------------------------------------------------
      if ( ishape==8 ) then
        call d3dwall ( xwal1, zwal1, mth, outmod, iotty )
      endif
c-----------------------------------------------------------------------
c     ishape=11 Dee-shaped conductor.
c-----------------------------------------------------------------------
      if ( ishape==11 ) then
        do i = 1, mth2
          the = (i-1) * dth
          plrad = 0.5 * ( xmax - xmin )
          xwal1(i) = xmax + plrad*( a + aw - aw*cos(the + dw*sin(the)) )
          zwal1(i) = - plrad * bw * sin(the)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=12 Solid bean-shaped conductor on right.
c-----------------------------------------------------------------------
      if ( ishape==12 ) then
        plrad = 0.5 * ( xmax-xmin )
        xmaj = 0.5 * ( xmax + xmin )
        a0 = plrad * ( 1.0 + aw - cw + a )
        brad = b * pye / 180.0
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          rho = a0 - aw*plrad*cos(the)
          the2 = brad * sin(the)
          xwal1(i) = xmaj + cw*plrad + rho * cos(the2)
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=13 Solid bean-shaped conductor on left.
c-----------------------------------------------------------------------
      if ( ishape==13 ) then
        plrad = 0.5 * ( xmax-xmin )
        xmaj = 0.5 * ( xmax + xmin )
        a0 = plrad * ( 1.0 + aw - cw + a )
        brad = b * pye / 180.0
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          rho = a0 + aw*plrad*cos(the)
          the2 = brad * sin(the)
          xwal1(i) = xmaj + cw*plrad - rho * cos(the2)
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=21 Shell scaled to plasma. Gap on the inner side.
c-----------------------------------------------------------------------
      if ( ishape==21 ) then
        plrad = 0.5 * ( xmax-xmin )
        xmaj = 0.5 * ( xmax + xmin )
        a0 = plrad * ( 1.0 + aw - cw + a )
        a0b = (a0 + plrad*aw)*bw
        brad0 = b * pye / 180.0
        brad = brad0
        blgrad0 = bbulg * pye / 180.0
        wcentr = xmaj + cw*plrad
        call adjustb ( blgrad0, blgrado, a, bw, cw, dw, xmaj, plrad,
     $       ishape )
        dthb = ( 2.0*aw*plrad / a0b )
     $       * ( 1.0 - sin(blgrado) ) / cos(blgrado)
        blgrad0 = blgrad0 - dthb
        call adjustb ( blgrad0, blgradi, a, bw, cw, dw, xmaj, plrad,
     $       ishape )
        do i = 1, mth2
          the0 = (i-1) * dth
          if ( the0 .gt. 0.5*pye .and. the0 .lt. 1.5*pye ) then
            thbulg = blgrado
          else
            thbulg = blgradi
          endif
          cost2b = cos(2*thbulg)
          the = the0
          cost = cos(the)
          ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
          rho = a0 - aw*plrad*ferm
          the2 = brad * sin(the)
          cost2 = cos(2.0*the2)
          fermb = 1.0 / ( exp( (cost2b - cost2)/tbulg ) + 1.0 )
          bulge = abulg*plrad*fermb
          xwal1(i) = xmaj + cw*plrad + rho * cos(the2+dw*sin(the2))
     $         + bulge
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=24 Shell scaled to plasma. Gap on the outer side.
c-----------------------------------------------------------------------
      if ( ishape==24 ) then
        plrad = 0.5 * ( xmax-xmin )
        xmaj = 0.5 * ( xmax + xmin )
        a0 = plrad * ( 1.0 + aw - cw + a )
        brad = b * pye / 180.0
        wcentr = xmaj + cw*plrad
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          cost = cos(the)
          ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
          rho = a0 + aw*plrad*ferm
          the2 = brad * sin(the)
          xwal1(i) = xmaj + cw*plrad - rho * cos(the2-dw*sin(the2))
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=31 Shell independent of plasma. Gap on the inner side.
c-----------------------------------------------------------------------
      if ( ishape==31 ) then
        a0 = a + aw
        a0b = (a0 + aw)*bw
        brad0 = b * pye / 180.0
        brad = brad0
        blgrad0 = bbulg * pye / 180.0
        wcentr = cw
        call adjustb ( blgrad0, blgrado, a, bw, cw, dw, xmaj, plrad,
     $       ishape )
        dthb = ( 2.0*aw / a0b )
     $       * ( 1.0 - sin(blgrado) ) / cos(blgrado)
        blgrad0 = blgrad0 - dthb
        call adjustb ( blgrad0, blgradi, a, bw, cw, dw, xmaj, plrad,
     $       ishape )
        do i = 1, mth2
          the0 = (i-1) * dth
          if ( the0 .gt. 0.5*pye .and. the0 .lt. 1.5*pye ) then
            thbulg = blgrado
          else
            thbulg = blgradi
          endif
          cost2b = cos(2.0*thbulg)
          the = the0
          cost = cos(the)
          ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
          rho = a0 - aw*ferm
          the2 = brad * sin(the)
          cost2 = cos(2.0*the2)
          fermb = 1.0 / ( exp( (cost2b - cost2)/tbulg ) + 1.0 )
          bulge = abulg*fermb
          xwal1(i) = cw + rho * cos(the2+dw*sin(the2)) + bulge
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=34 Shell independent of plasma. Gap on the outer side.
c-----------------------------------------------------------------------
      if ( ishape==34 ) then
        a0 = a + aw
        brad = b * pye / 180.0
        wcentr = cw
        do i = 1, mth2
          the0 = (i-1) * dth
          the = the0
          cost = cos(the)
          ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
          rho = a0 + aw*ferm
          the2 = brad * sin(the)
          xwal1(i) = cw - rho * cos(the2-dw*sin(the2))
          zwal1(i) = - bw * rho * sin(the2)
        enddo
      endif
c-----------------------------------------------------------------------
c     ishape=41 Arbitrary wall generated by spline data in wall_geo.in.
c
c               npots0 is number of points
c               wcentr is radial offset
c               thetatmp, xwaltmp, zwaltmp are columns containing angle
c               and position of each point for generating splines
c-----------------------------------------------------------------------
      if ( ishape==41 ) then
        open(unit=41,file="wall_geo.in",status='old',
     $             form='FORMATTED')
        read(41,*) npots0
        read(41,*) wcentr
        read(41,*)
c       follow nth=nth0+5  nth0=mth  mth2=mth+2 mth1=mth+1=npots0
        npots=npots0+5-1
        allocate (thetatmp(npots),xwaltmp(npots),zwaltmp(npots),
     $       rioptmp(2),xpptmp(npots),ww1tmp(npots),ww2tmp(npots),
     $       ww3tmp(npots),tabtmp(3))
        do i = 1, npots0
          read (41,*) thetatmp(i),xwaltmp(i),zwaltmp(i)
          xwaltmp(i)=xwaltmp(i)-wcentr
        enddo
        close (41)
        rioptmp(1)=4
        rioptmp(2)=4
        call spl1d1(npots0,thetatmp,xwaltmp,xpptmp,rioptmp,1,
     $              ww1tmp,ww2tmp,ww3tmp)
        do i = 1, mth1
          the0 = (i-1)*dth
          call spl1d2(npots0,thetatmp,xwaltmp,xpptmp,1,the0,tabtmp)
          xwal1(i)=tabtmp(1)*a+wcentr
        enddo
        rioptmp(1)=4
        rioptmp(2)=4
        call spl1d1(npots0,thetatmp,zwaltmp,xpptmp,rioptmp,1,
     $              ww1tmp,ww2tmp,ww3tmp)
        do i = 1, mth1
          the0 = (i-1)*dth
          call spl1d2(npots0,thetatmp,zwaltmp,xpptmp,1,the0,tabtmp)
          zwal1(i)=tabtmp(1)*a
        enddo
        xwal1(1)=xwal1(mth1)
        zwal1(1)=zwal1(mth1)
        xwal1(mth2)=xwal1(2)
        zwal1(mth2)=zwal1(2)
        deallocate(thetatmp,xwaltmp,zwaltmp,xpptmp,ww1tmp,ww2tmp,ww3tmp,
     $             rioptmp,tabtmp)
      endif
c-----------------------------------------------------------------------
c     ishape=42 Arbitrary wall generated by position data in wall_geo.in.
c               Similar to ishape=42, except points are used directly without
c               splining.  Therefore, thetatmp is read but ignored.
c
c               npots0 is number of points and must be equal to mth+2.
c               wcentr is radial offset
c               thetatmp, xwaltmp, zwaltmp are columns containing angle
c               and position of each point in wall
c-----------------------------------------------------------------------
      if ( ishape==42 ) then
        open(unit=41,file="wall_geo.in",status='old',
     $       form='FORMATTED')
        read(41,*) npots0
        read(41,*) wcentr
        read(41,*)
        if( npots0 .ne. mth+2 ) then
          write(*,'(A,A)') 'ERROR: Number of points in wall_geo.in',
     $         ' must be equal to mth+2.'
          stop
        endif
        allocate (thetatmp(npots0))
        do i = 1, npots0
          read (41,*) thetatmp(i),xwal1(i),zwal1(i)
        enddo
        close (41)
        deallocate(thetatmp)
        if( xwal1(mth1).ne.xwal1(1) .or. zwal1(mth1).ne.zwal1(1) ) then
          write(*,'(A,A)') 'ERROR: First point in wall_geo.in',
     $         ' must be equal to 2nd last point.'
          stop
        endif
        if( xwal1(mth2).ne.xwal1(2) .or. zwal1(mth2).ne.zwal1(2) ) then
          write(*,'(A,A)') 'ERROR: Last point in wall_geo.in',
     $         ' must be equal to 2nd point.'
          stop
        endif
      endif

      xmx = xma + xshift

c-----------------------------------------------------------------------
c     This section is never called
c-----------------------------------------------------------------------
c      if( .false. ) then
c        call spl1d1(mth1,thet,xwal1,xpp,iop,1,ww1,ww2,ww3)
c        call spl1d1(mth1,thet,zwal1,zpp,iop,1,ww1,ww2,ww3)
c        do i=2,mthalf-1
c          xs = xinf(i)
c          zs = zinf(i)
c          xt = xwal1(i)
c          zt = zwal1(i)
c          tt = thet(i)
c          do k=1,20
c            call spl1d2(mth1,thet,xwal1,xpp,1,tt,tabx)
c            call spl1d2(mth1,thet,zwal1,zpp,1,tt,tabz)
c            xt = tabx(1)
c            zt = tabz(1)
c            xmx1 = xs - xt
c            zmz1 = zs - zt
c            f = xmx1 * tabx(2) + zmz1 * tabz(2)
c            fp = xmx1*tabx(3) +zmz1*tabz(3) -(tabx(2))**2 -(tabz(2))**2
c            delt = f / fp
c            if(abs(delt).lt.1.0e-4 .and. abs(f) .lt. 1.0e-4) exit
c            tt = tt - delt
c            if( k==20 ) then
c              write(iotty,1430)i,xs,zs,xt,zt,tt,f,fp
c              write(outpest,1430)i,xs,zs,xt,zt,tt,f,fp
c              call errmes(outpest,'vacdat')
c            endif
c          enddo
c          ww1(i) = xt
c          ww2(i) = zt
c        enddo
c        do i = 2, mthalf - 1
c          iq1 = mth1 - i + 1
c          xwal1(i) = ww1(i)
c          xwal1(iq1) = ww1(i)
c          zwal1(i) = ww2(i)
c          zwal1(iq1) = - ww2(i)
c        enddo
c      endif

c-----------------------------------------------------------------------
c     Do cleanup
c-----------------------------------------------------------------------
  145 continue
      if ( leqarcw .eq. 1 ) then
         call  eqarcw ( xwal1,zwal1, xpp,zpp, ww1,ww2,ww3, mth1 )
         do i = 1, mth1
            xwal1(i) = xpp(i)
            zwal1(i) = zpp(i)
         enddo
      endif
      call savewall(wcentr,xwal1,zwal1) !Write wall shape to disk
      if ( iplt .le. 0 ) then
        xmx = xmaj
        zma = 0.0
        iplt = 1
      endif
      if(insect) then
        write(iotty,1450)inside
        write(outpest,1450)inside
        write(outmod,1450)inside
        call errmes(outpest,'vacdat')
      endif
      aw = awsave !restore value of aw that was saved at beginning of function
      bw = bwsave !restore value of bw that was saved at beginning of function
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine wwall
c-----------------------------------------------------------------------
c     subprogram 7. d3dwall.
c     defines the shape of the DIII-D wall.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine d3dwall ( xwall, ywall, mthh, iomod, iotty1 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      integer, parameter :: ncdf=26
      dimension xwall(*), ywall(*), rwi(ncdf), zwi(ncdf)
      integer :: nwcoef=ncdf
      data rwall0/ 1.6400000_8/, zwall0/ 0.0000000_8/,
     $     awall0/ 0.8839410_8/, ewall0/ 1.4037020_8/
      data (rwi(k),k=1,ncdf)
     $     / 0.1000000d+01, 0.5526794d-01,-0.1738114d+00, 0.1850757d-01
     $     , 0.3714965d-01,-0.2882647d-01,-0.2357329d-02, 0.9548103d-02
     $     ,-0.1214923d-01,-0.1853416d-02, 0.6837493d-02,-0.1711245d-02
     $     , 0.2270762d-02, 0.3689963d-02,-0.3959393d-02,-0.1098017d-02
     $     , 0.3745465d-02,-0.2157904d-03,-0.3977743d-03,-0.2725623d-03
     $     ,-0.1005857d-02,-0.4579016d-05, 0.2396789d-02,-0.7057043d-03
     $     , 0.1158347d-02, 0.3552319d-03/
      data (zwi(k),k=1,ncdf)
     $     / 0.1000000d+01,-0.3236632d-01,-0.1629422d+00, 0.6013983d-01
     $     , 0.1167756d-01,-0.2579542d-01, 0.1626464d-01,-0.2085857d-02
     $     ,-0.9098639d-02, 0.1022163d-01,-0.4388253d-02,-0.9367258d-02
     $     , 0.8308497d-02, 0.4765150d-02,-0.4611675d-02,-0.1121423d-02
     $     ,-0.2501100d-03, 0.4282634d-03, 0.2669702d-02,-0.1073800d-02
     $     ,-0.2191338d-02, 0.1328267d-02, 0.5050959d-03,-0.5758863d-03
     $     , 0.9348883d-03, 0.7094351d-03/
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      rwll      = rwall0
      zwll      = zwall0
      awll      = awall0*rext
      ewll      = ewall0
      zstart    = zlim
      zlim = 0.0
      zstart = 0.0
      nwalp = mthh
      rext = 1.0
      awll = awall0*rext
      call d3dvesl(rwll,zwll,awll,ewll,rwi,zwi,nwcoef,zstart
     $     ,xwall,ywall,nwalp,ier)
      write ( iomod, '("ier in d3dwall = ", i3)' ) ier
      write ( iotty1,  '("ier in d3dwall = ", i3)' ) ier
      xwall(mthh+1) = xwall(1)
      ywall(mthh+1) = ywall(1)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine d3dwall
c-----------------------------------------------------------------------
c     subprogram 8. d3dvesl.
c     defines the shape of the DIII-D vacuum vessel.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine d3dvesl(r0,z0,a0,e0,ar,az,nval,zst,r,z,npts,ier)
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension ar(nval),az(nval)
      dimension r (npts),z (npts)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      pii    =  3.1415926535897932385_8
      ier    = 0
      if ( abs(z0-zst) .le. 1.0e-6 ) then
      else
         if(z0 .lt. zst) isgn  = +1
         if(z0 .gt. zst) isgn  = -1
         zpcmp    = (zst-z0)/(a0*e0)
         arci     = 0.0
         arcf     = 0.5*pii
         dfi      = (arcf-arci)/npts
         arca     = arci
         zza      = 0.0
         do 20 j  = 2,npts
            arcb     = arca + isgn*dfi
            zzb      = 0.0
            do  5 k  = 1,nval
               ackb     = k*arcb
               zzb      = zzb + az(k)*sin(ackb)
 5          continue
            if((zza-zpcmp)*(zzb-zpcmp) .le. 0.0) then
               arcc     = arcb + isgn*dfi
               zzc      = 0.0
               do 10 k  = 1,nval
                  ackc     = k*arcc
                  zzc      = zzc + az(k)*sin(ackc)
 10            continue
               go to 25
            else
               arca     = arcb
               zza      = zzb
            endif
 20      continue
         ier    = 1
         return
 25      continue
         dzp    = zzc - zzb
         dzm    = zzb - zza
         dzt    = zzc - zza
         dcf1   = dfi*(dzm/dzp + dzp/dzm)/dzt
         dcf2   = dfi*(1.0/dzp - 1.0/dzm)/dzt
         zdf    = zpcmp - zzb
         arcst  = arcb + dcf1*zdf + dcf2*zdf*zdf
      endif
      arc0     =  arcst
      arc1     =  arcst + 2.0*pii
      darc     = (arc1-arc0)/npts
      do 100 j = 1,npts
         arc      = arc0 + (j-1.)*darc
         sumr     = 0.0
         sumz     = 0.0
         do 50 k  = 1,nval
            arck     = k*arc
            sumr     = sumr + ar(k)*cos(arck)
            sumz     = sumz + az(k)*sin(arck)
 50      continue
         rpval    = r0  +    a0*sumr
         zpval    = z0  - e0*a0*sumz
         r(j)     = rpval
         z(j)     = zpval
 100  continue
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine d3dvesl
c-----------------------------------------------------------------------
c     subprogram 9. eqarcw
c     computes equal-arc-length coordinates for wall.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine eqarcw ( xin, zin, xout, zout, ell, thgr, thlag, mw1 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension xin(*), xout(*), zin(*), zout(*),
     $     ell(*), thgr(*), thlag(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do  iw = 1, mw1
         thlag(iw) = (1.0/(mw1-1))*(iw-1)
      enddo
      ell(1) = 1.e-8
      do  iw = 1, mw1
         if ( iw .ne. 1 ) then
            thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
            call lag ( thlag,xin,mw1,3,thet,f,df,1 )
            xtzt = df
            call lag ( thlag,zin,mw1,3,thet,f,df,1 )
            xtzt = sqrt ( xtzt**2 + df**2 )
            ell(iw) = ell(iw-1) + xtzt/(mw1-1)
         endif
      enddo
      do i = 1, mw1
         elgr = ( ell(mw1)/(mw1-1) ) * (i-1)
         call lag ( ell,thlag,mw1,3,elgr,f,df,0 )
         thgr(i) = f
      enddo
      do  i = 1, mw1
         ttt = thgr(i)
         call lag ( thlag,xin,mw1,3,ttt,f,df,0 )
         xout(i) = f
         call lag ( thlag,zin,mw1,3,ttt,f,df,0 )
         zout(i) = f
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine eqarcw
c-----------------------------------------------------------------------
c     subprogram 10. gatonorm.
c     normalize vacuum matrix to GATO.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine gatonorm ( vacin, gatovac_, nd, rgato_,mfel_,mth_,
     $     qa1_,twopi_ )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension  vacin(nd,nd), gatovac_(nd,nd)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do m1 = 1, mfel_
         do m2 = 1, mfel_
            gatovac_(m1,m2) = rgato_*vacin(m1,m2) / (twopi_*mfel_)
         enddo
      enddo
      open ( unit=36,file='vacgato' )
      write ( 36, '("VACGATO: rgato_, mfel_,mth_, qa1_ = ",/,
     $     1pe13.5, 2i4, 1pe13.5 / )' ) rgato_, mfel_,mth_, qa1_
      do m1 = 1, mfel_
         write ( 36,'(/, i5 )') m1
         write ( 36,'(10e13.5)') (gatovac_(m1,m2), m2 = 1, mfel_)
      enddo
      close ( unit=36 )
      return
      end subroutine gatonorm
c-----------------------------------------------------------------------
c     subprogram 11. adjustb.
c     Adjusts for the skewing of the subtending bulge angle due to
c     elongation and triangularity. Only good for small angles.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine adjustb(betin,betout,a_,bw_,cw_,dw_,xmaj_,plrad_,
     $   ishape_)
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if ( ishape_ .eq. 31 ) then
         r0 = cw_
         r  = a_
      endif
      if ( ishape_ .eq. 21 ) then
         r0 = xmaj_ + cw_*plrad_
         r  = plrad_ * ( 1.0 + a_ - cw_ )
      endif
      bet2 = betin
      betout = abs ( atan ( tan(bet2) / bw_ ) )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine adjustb
c-----------------------------------------------------------------------
c     subprogram 12. fouran.
c     Fourier analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine fouran ( gij, gil, cs, m00,l00 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      REAL(r8), DIMENSION(nths,nths) :: gij
      REAL(r8), DIMENSION(nths2,nfm2) :: gil
      REAL(r8), DIMENSION(nths,nfm) :: cs
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      do l1 = 1, jmax1
         do i = 1, mth
            gil(m00+i,l00+l1) = 0.0
         enddo
      enddo
      do l1 = 1, jmax1
         ll = l1 - 1 + lmin(1)
         do j = 1, mth
            do i = 1, mth
               gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $              cs(j,l1)*gij(i,j)
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine fouran
c-----------------------------------------------------------------------
c     subprogram 13. fouranv.
c     inverse Fourier analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine foranv ( gil, gll, cs, m00,l00 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension gil(nths2,nfm2), gll(nfm,nfm), cs(nths,nfm)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            gll(l1,l2) = 0.0
         enddo
      enddo
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            do i = 1, mth
               gll(l2,l1) = gll(l2,l1) +
     $              dth*cs(i,l2)*gil(m00+i,l00+l1)*twopi
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine foranv
c-----------------------------------------------------------------------
c     subprogram 14. foura2.
c     more Fourier analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine foura2 ( gij,m01,m02, gil, m00 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      real nq
      dimension gij(nths2,nths2), gil(nths,nths)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
      do l1 = 1, jmax1
         do i = 1, mth
            gil(m00+i,l1) = 0.0
            gil(m00+i,jmax1+l1) = 0.0
         enddo
      enddo
      do l1 = 1, jmax1
         ll = l1 - 1 + lmin(1)
         do j = 1, mth
            theta = (j-1) * dth
            elth = ll * theta
            elthnq = elth + nq*delta(j)
            sinlth = sin(elthnq)
            coslth = cos(elthnq)
            do i = 1, mth
               gil(m00+i,l1) = gil(m00+i,l1) + coslth*gij(i+m01,j+m02)
               gil(m00+i,jmax1+l1) = gil(m00+i,jmax1+l1)
     $              + sinlth*gij(i+m01,j+m02)
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      return
      end subroutine foura2
c-----------------------------------------------------------------------
c     subprogram 15. fanal.
c     yet more Fourier analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine fanal ( fth, nt, flc,fls, l1,l2, pi,ddt )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension fth(*), flc(*), fls(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      dt0 = 1.0 / nt
      dt = 2.0 * pi * dt0
      fth(nt+1) = fth(1)
      nl = l2 - l1 + 1
      do l = 1, nl
         ll = l1 - 1 + l
         flc(l) = 0.0
         fls(l) = 0.0
         do i = 1, nt
            th = dt*(i-1 + ddt)
            flc(l) = flc(l) + fth(i)*cos(ll*th)
            fls(l) = fls(l) + fth(i)*sin(ll*th)
         enddo
         flc(l) = dt0*flc(l)
         fls(l) = dt0*fls(l)
      enddo
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      return
      end subroutine fanal
c-----------------------------------------------------------------------
c     subprogram 16. fanal1.
c     yet more Fourier analysis.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine fanal1 ( gi,ndi1,ndi2,mi1, gor,goi,ndo1,ndo2 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      real nq
      dimension gi(ndi1,ndi2), gor(ndo1,ndo2), goi(ndo1,ndo2)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
      do  l1 = 1, jmax1
         do l2 = 1, jmax1
            gor(l2,l1) = 0.0
            goi(l2,l1) = 0.0
         enddo
      enddo
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            ll2 = l2 - 1 + lmin(1)
            do i = 1, mth
               elth = ll2*(i-1)*dth 
               gor(l1,l2) = gor(l1,l2)
     $              + cos(elth) * gi(mi1+i,l1)
     $              + sin(elth) * gi(mi1+i,jmax1+l1)
               goi(l1,l2) = goi(l1,l2)
     $              + cos(elth) * gi(mi1+i,jmax1+l1)
     $              - sin(elth) * gi(mi1+i,l1)
            enddo
            gor(l1,l2) = gor(l1,l2) * dth / twopi
            goi(l1,l2) = goi(l1,l2) * dth / twopi
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine fanal1
c-----------------------------------------------------------------------
c     subprogram 17. felang.
c     normalize finite elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine felang ( gij, gil, cs, m00,l00 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      data izcal / 0 /
      real nq
      dimension gij(nths,nths), gil(nths2,nfm2), cs(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      izcal = izcal + 1
      nzwrt = 8
      nzwdel = mth / nzwrt
      znorm = sqrt( float(mfel) )
      zwt1 = 0.5 * znorm
      zwt2 = 1.0 * znorm
      zws1 = 1.0 * znorm /3.0
      zws2 = 2.0 * zws1
      zws4 = 2.0 * zws2
      nzdel = ndfel
      nzdel1 = nzdel + 1
      q = qa1
      nq = n*q
      do l1 = 1, mfel
         do i = 1, mth
            gil(m00+i,l00+l1) = 0.0
         enddo
      enddo
      do i = 1, mth
         do l1 = 1, mfel
            izwrt = 0
            if ( (i .eq. (i/nzwdel * nzwdel + 1))
     $           .and. ( l1 .eq. (l1/8 * 8 +1 ))
     $           .and. izcal .eq. 1 ) then
               izwrt = 1
            endif
            mzl = (l1-1) * nzdel + 1
            mzr = mzl + nzdel
            do j = 1, nzdel1
               jth0 = mzl + j - 2
               jth = mod(jth0,mth)
               if ( jth .lt. 0 ) jth = mth + jth
               jth1 = jth + 1
               if ( nzdel .eq. 1 ) then
                  gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $                 zwt1 * gij(i,jth1) * cs(jth1)
               else
                  zwt = zws2
                  if ( (j/2)*2 .eq. j ) zwt = zws4
                  if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                  gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $                 zwt * gij(i,jth1) * cs(jth1)
               endif
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine felang
c-----------------------------------------------------------------------
c     subprogram 18. felanv.
c     normalize finite elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine felanv ( gil, gll, cs, m00,l00 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      data izcal / 0 /
      real nq
      dimension gil(nths2,nfm2), gll(nfm,nfm), cs(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      izcal = izcal + 1
      znorm = sqrt( float(mfel) )
      zwt1 = 0.5 * dth * znorm
      zwt2 = dth * znorm
      zws1 = dth * znorm /3.0
      zws2 = 2.0 * zws1 
      zws4 = 2.0 * zws2 
      nzdel = ndfel
      nzdel1 = nzdel + 1
      q = qa1
      nq = n*q
      do l1 = 1, mfel
         do l2 = 1, mfel
            gll(l1,l2) = 0.0
         enddo
      enddo
      do l2 = 1, mfel
         do l1 = 1, mfel
            izwrt = 0
            if ((l2 .eq. (l2/8 * 8 + 1))
     $           .and. (l1 .eq. (l1/8 * 8 +1))
     $           .and. (izcal .eq. 1)) then
               izwrt = 1
            endif
            mzl = (l1-1) * nzdel + 1
            mzr = mzl + nzdel
            do j = 1, nzdel1
               jth0 = mzl + j - 2
               jth = mod(jth0,mth)
               if ( jth .lt. 0 ) jth = mth + jth
               jth1 = jth + 1
               if ( nzdel .eq. 1 ) then
                  gll(l1,l2) = gll(l1,l2) +
     $                 zwt1 * gil(m00+jth1,l00+l2) * cs(jth1) * twopi
               else
                  zwt = zws2
                  if ( (j/2)*2 .eq. j ) zwt = zws4
                  if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                  if ( jth1 .ne. mth ) then
                     gll(l1,l2) = gll(l1,l2) +
     $                    zwt * gil(m00+jth1,l00+l2) * cs(jth1) * twopi
                  else 
                     gll(l1,l2) = gll(l1,l2)
     $                    + zwt * gil(m00+1,l00+l2) * cs(jth1) * twopi
                  endif
               endif
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine felanv
c-----------------------------------------------------------------------
c     subprogram 19. fotofi.
c     Transforms the Vacuum matrix form Fourier-l space to finite
c     elements-i space.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine fotofi ( vin,vout, scnlth, wrk1, wrk2, iopsc )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension vin(nfm,nfm), vout(nfm,nfm), scnlth(nths,nfm),
     $     wrk1(nfm,nfm), wrk2(nfm,nfm)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      call tmat ( scnlth, wrk1, iopsc )
      call matmul3 ( wrk1, vin, nfm,nfm, mfel,jmax1,jmax1, wrk2,nfm ) 
      call mtrans ( wrk1, nfm, nfm )
      call matmul3 ( wrk2, wrk1, nfm,nfm, mfel,jmax1,mfel, vout,nfm )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine fotofi
c-----------------------------------------------------------------------
c     subprogram 20. orchek.
c     checks for orthogonality of transformation matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine orchek ( wrkr, wrki, wrkrt, wrkit, wrko1,wrko2 )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension  wrkr(nfm,nfm), wrki(nfm,nfm),
     $     wrkrt(nfm,nfm), wrkit(nfm,nfm),
     $     wrko1(nfm,nfm), wrko2(nfm,nfm)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      iopc = 1
      iops = 2
      call tmat ( coslt, wrkrt, iopc )
      call tmat ( sinlt, wrkit, iops )

      do i = 1, mfel
         do m = 1, jmax1
            wrkr(i,m) = wrkrt(i,m)
            wrki(i,m) = wrkit(i,m)
         end do
      end do
      call mtrans ( wrkrt, nfm, nfm )
      call mtrans ( wrkit, nfm, nfm )
      call matmul3 ( wrkr,wrkrt, nfm,nfm, mfel,jmax1,mfel, wrko1,nfm )
      call matmul3 ( wrki,wrkit, nfm,nfm, mfel,jmax1,mfel, wrko2,nfm )
      do i = 1, mfel
         do j = 1, mfel
            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
         end do
      end do
      call matwrtn ( wrko1,nfm,nfm,1,1,mfel,mfel,8,8,
     $     "Real TMAT*trans-TMAT: kk", outmod,iotty )
      call matmul3 ( wrkrt,wrkr, nfm,nfm, jmax1,mfel,jmax1, wrko1,nfm )
      call matmul3 ( wrkit,wrki, nfm,nfm, jmax1,mfel,jmax1, wrko2,nfm )
      do i = 1, jmax1
         do j = 1, jmax1
            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
         end do
      end do
      call matwrtn ( wrko1,nfm,nfm,1,1,jmax1,jmax1,8,8,
     $     "Real TMAT*trans-TMAT: ll", outmod,iotty )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine orchek
c-----------------------------------------------------------------------
c     subprogram 21. tmat.
c     calculates the matrix for transforming from Fourier to finite
c     elements. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine tmat ( sil, tll, iop )
      USE vglobal_mod
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension sil(nths,nfm), tll(nfm,nfm)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      znorm = sqrt( float(mfel) )
      jmax1 = lmax(1) - lmin(1) + 1
      pi = pye
      do i = 1, nfm
         do j = 1, nfm
            tll(i,j) = 0.0 
         end do
      end do
      if ( iop .eq. 0 ) then
         zwt1 = 0.5 * dth * znorm
         zwt2 = dth * znorm
         zws1 = dth * znorm /3.0
         zws2 = 2.0 * zws1
         zws4 = 2.0 * zws2
         nzdel = ndfel
         nzdel1 = nzdel + 1
         do l1 = 1, mfel
            mzl = (l1-1) * nzdel + 1
            mzr = mzl + nzdel
            do j = 1, nzdel1
               jth0 = mzl + j - 2
               jth = mod(jth0,mth)
               if ( jth .lt. 0 ) jth = mth + jth
               jth1 = jth + 1
               do l2 = 1, jmax1
                  if ( nzdel .eq. 1 ) then
                     tll(l1,l2) = tll(l1,l2) + zwt1 * sil(jth1,l2)
                  else
                     zwt = zws2
                     if ( (j/2)*2 .eq. j ) zwt = zws4
                     if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                     tll(l1,l2) = tll(l1,l2) + zwt * sil(jth1,l2)
                  endif
               enddo
            enddo
         enddo
         do mm = 1, mfel
            do ll = 1, jmax1
               tll(mm,ll) = tll(mm,ll) / twopi
            end do
         end do
      endif
      do ll = 1, jmax1
         la = lmin(1) - 1 + ll
         if ( la .ne. 0 )
     $        zslpn = znorm * sin(la*pi/mfel) / ( la * pi )
         do mm = 1, mfel
            if ( (iop .eq. 1) .and. (la .ne. 0) )
     $           tll(mm,ll) = cos(la*pi*(2*mm-1)/mfel) * zslpn
            if ( (iop .eq. 1) .and. (la .eq. 0) )
     $           tll(mm,ll) = 1.0 / znorm
            if ( (iop .eq. 2) .and. (la .ne. 0 ) )
     $           tll(mm,ll) = sin(la*pi*(2*mm-1)/mfel) * zslpn
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine tmat
c-----------------------------------------------------------------------
c     subprogram 22. wtopest.
c     conversion to PEST?
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine wtopest ( vacpstr, vacpsti, xwal_, zwal_ )
            USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      dimension vacpstr(*), vacpsti(*), xwal_(*), zwal_(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jmax1 = lmax(1) - lmin(1) + 1
      len = 2*nfm*nfm + 2*mth2 + 10
      ndsk = 1
      write (  iotty, '("Writing to unit iovac: iovac, jmax1, len = ",
     $     /, 3i6 )' ) iovac, jmax1, len
      write ( outmod, '("writing to unit iovac: iovac, jmax1, len = ",
     $     /, 3i6 )' ) iovac, jmax1, len
      call zop(iovac,"vacmat",len,ndsk,iiff,999)
      lmn = lmin(1)
      lmx = lmax(1)
      lgivup = 1
      nadres = 1
      call zwr(iovac,real(lmn,r8),1,nadres,lgivup,999)
      nadres = nadres + 1
      call zwr(iovac,real(lmx,r8),1,nadres,lgivup,999)
      nadres = nadres + 1
      call trans ( xwal_,mth, xjdtxj,mthin )
      call zwr(iovac,xjdtxj(1),mthin2,nadres,lgivup,999)
      nadres = nadres + mthin2
      call trans ( zwal_,mth, xjdtxj,mthin )
      call zwr(iovac,xjdtxj(1),mthin2,nadres,lgivup,999)
      nadres = nadres + mthin2
      length = jmax1**2
      call zwr(iovac,vacpstr(1),length,nadres,lgivup,999)
      nadres = nadres + jmax1**2
      call zwr(iovac,vacpsti(1),length,nadres,lgivup,999)
      write (  iotty, '("Wrote to iovac, lmn, lmx, nadres = ",
     $     /, 3i6 )' ) lmn, lmx, nadres
      write ( outmod, '("Wrote to iovac, lmn, lmx, nadres = ",
     $     /, 3i6 )' ) lmn, lmx, nadres
      call zcl( iovac, 999 )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
 999  call errmes ( outpest, 'wtopest' )
      end subroutine wtopest
c-----------------------------------------------------------------------
c     subprogram 23. savewall
c     After processing the wall, write the coordinates to be used
c     to wall_geo.out.
c-----------------------------------------------------------------------
      subroutine savewall(wcentr,xwal1,zwal1)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      USE vglobal_mod
      !implicit none
      IMPLICIT REAL(r8) (a-h,o-z)

      real(r8) :: wcentr
      dimension xwal1(*), zwal1(*)
      INTEGER :: firstpos, i
      real(r8) :: thwal0,thwal1
      logical :: retro
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      retro=.false.

      open(unit=41,file="wall_geo.out",status='unknown',
     $             form='FORMATTED')
      write(41,'(1x,a7,i6)') 'mth2 = ', mth2
      write(41,'(1x,a7,g14.6)') 'ro = ',wcentr
      write(41,'(A)') '  theta(rad)    R(m)          Z(m)'

c     find first point above angle zero (slightly clockwise from R-direction)
      thwal0=atan2(-zwal1(mth1),xwal1(mth1)-wcentr)
      do i=1,mth1
         thwal1 = atan2(-zwal1(i),xwal1(i)-wcentr)
         if(thwal1>0.0 .and. thwal0<0.0) then
            firstpos=i
            exit
         endif
         thwal0=thwal1
      enddo

c     first point
      i=firstpos-1
      thwal1 = atan2(-zwal1(i),xwal1(i)-wcentr)
      write(41,'(3e14.6)') thwal1,xwal1(i),zwal1(i)
c     points in clockwise direction at end of array
      do i=firstpos,mth
         thwal0 = thwal1
         thwal1 = atan2(-zwal1(i),xwal1(i)-wcentr)
         if(thwal1<0.0) thwal1=thwal1+twopi
         if( thwal1<thwal0 ) retro=.true.
         write(41,'(3e14.6)') thwal1,xwal1(i),zwal1(i)
      enddo
c     points in clockwise direction at beginning of array
      do i=1,firstpos-1
         thwal0 = thwal1
         thwal1 = atan2(-zwal1(i),xwal1(i)-wcentr)
         if(thwal1<0.0) thwal1=thwal1+twopi
         if( thwal1<thwal0 ) retro=.true.
         write(41,'(3e14.6)') thwal1,xwal1(i),zwal1(i)
      enddo
c     last point
      i=firstpos
      thwal0 = thwal1
      thwal1 = atan2(-zwal1(i),xwal1(i)-wcentr)
      thwal1=thwal1+twopi
      if( thwal1<thwal0 ) retro=.true.
      write(41,'(3e14.6)') thwal1,xwal1(i),zwal1(i)
      close (41)

      if( retro ) then
         write(*,*) "!! WARNING: Retrograde angles present in wall"
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end subroutine
