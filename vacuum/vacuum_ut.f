c-----------------------------------------------------------------------
c     file vacuum_ut.f.
c     utility routines used by Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. errmes
c      2. wrtout
c      3. vecwrt
c      4. writg1
c      5. errnc
c      6. nwopn
c      7. nwinq
c      8. nwdid
c      9. nwdinq
c     10. nwvid
c     11. nwvinq
c     12. nwvg1i
c     13. nwvgt1
c     14. nwvg1c
c     15. nwvgti
c     16. nwvgt
c     17. nwvgtc
c     18. nwainq
c     19. nwagtc
c     20. nwagt
c     21. nwanam
c     22. writn1
c     23. writnn
c     24. writa1
c     25. writv1
c     26. writis
c     27. writvs
c     28. vfrmt1
c     29. vfrmt2
c     30. vfrmt3
c     31. bounds
c     32. boundsi
c     33. macopy
c     34. matwrt
c     35. vacasym
c     36. masym0
c     37. vacasymi
c     38. msctimer
c     39. shftpi
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 1. errmes.
c     prints error message and terminates program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine errmes ( nout,mesage )
      use local_mod, only: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer, parameter :: iotty=86
      character*(*) mesage
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
  100 format(" error in subroutine",1x, a10/)
c-----------------------------------------------------------------------
c     write error message.
c-----------------------------------------------------------------------
      write(iotty,100)mesage
      write(nout,100)mesage
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      stop
      end
c-----------------------------------------------------------------------
c     subprogram 2. wrtout.
c     prints message.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine wrtout ( nsrf, vec, wvec, m1,m2 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension vec(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
  100 format ( //,1x, 'nsrf = ',i3, 2x, a8, / )
  200 format ( 1x, 1p10e13.5 )
c-----------------------------------------------------------------------
c     write message.
c-----------------------------------------------------------------------
      write ( 16, 100 ) nsrf, wvec
      write ( 16, 200 ) ( vec(i), i = m1,m2 )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. vecwrt.
c     writes long message.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vecwrt ( mm, vec, wvec, m1,m2, nout1,nout2 )
      use local_mod, only: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension vec(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
  100 format ( //,1x, 'mm = ',i3, 2x, a8, / )
  200 format ( 1x, 1p10e13.5 )
c-----------------------------------------------------------------------
c     write statements.
c-----------------------------------------------------------------------
      if ( nout1 .ne. 0 ) then
         write ( nout1, 100 ) mm, wvec
         write ( nout1, 200 ) ( vec(i), i = m1,m2 )
      endif
      if ( nout2 .ne. 0 ) then
         write ( nout2, 100 ) mm, wvec
         write ( nout2, 200 ) ( vec(i), i = m1,m2 )
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. writg1.
c     prints character string.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writg1 ( string, seps, iseps, nout )
      use local_mod, only: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension nout(*)
      character*(*) string
      character*(*) seps
c-----------------------------------------------------------------------
c     write message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 )then
            if ( iseps .eq. 0 )
     $           write ( nout(io), '(/,a)' ) string
            if ( iseps .eq. 1 )
     $           write ( nout(io), '(/, a,/, 5x, a)' ) seps, string
            if ( iseps .eq. 2 )
     $           write ( nout(io), '(/, 5x,a, /,a )' ) string, seps
            if ( iseps .eq. 3 )
     $           write ( nout(io), '( /, a ,/, 5x,a,/,a )' )
     $           seps, string, seps
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 5. errnc.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine errnc ( mesage, nrcode, nout )
      use local_mod, only: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      integer nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( "nrcode in ", a16, " = ", i5 )
c-----------------------------------------------------------------------
c     print messasge.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, nrcode
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 6. nwopn.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwopn ( cdfid, file, nrw, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, nout(*)
      character*(*) file, nrw
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"cdfid = ", i3, /,
     $     "ncopn called for file, ",a, " ---nrw = ", a )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncopn", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, file, nrw
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 7. newinq
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwinq ( cdfid, ndims, nvars, ngatts, recid, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, recid, nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /, "ncinq: cdfid = ", i3, /,
     $     "number of dimensions, ndims is ...", i5, /,
     $     "number of variables, nvars is ...", i5, /,
     $     "number of global attributes, ngatts is ...", i5, /,
     $     "ID of the unlimited dimension, recid is ...", i5 )
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncinq", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, ndims, nvars,
     $           ngatts, recid
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 8. nwdid.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwdid ( cdfid, dimnam, ndid, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, nout(*)
      character*(*) dimnam
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncdid: cdfid = ", i3,
     $     " dimension name is ...", a, /, "dimension id is ...",i5 )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncdid", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, dimnam, ndid
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 9. nwdinq.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwdinq ( cdfid, dimid, dimnam, dimsiz, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, dimid, dimsiz, nout(*)
      character*(*) dimnam
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncdinq: cdfid = ", i3,
     $     ".  Dimension name, id, and size = ",/ ,a, 2i5 )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncdinq", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, dimnam, dimid, dimsiz
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 10. nwvid.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvid ( cdfid, varnam, vid, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, nout(*)
      character*(*) varnam
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncvid: cdfid = ", i3, " varnam= ",a ,/,
     $     "variable id is ...", i5 )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncvid", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, varnam, vid
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 11. nwvinq.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvinq ( cdfid, vid, varnam, vartyp, nvdims,
     $     vdims, nvatts, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, vartyp, vdims(*), nout(*)
      character*(*) varnam
      character*8 fmt
      character*72 fmt1
      character*72 fmt2
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      call vfrmt2 ( nvdims, "i5", fmt )
      fmt1 = '/,"ncvinq: cdfid = ", i3, ".  Vid= ",i5,/,'
      fmt2 = '( "Vector of Dimension ID''s is ...",'//fmt//')'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"variable type and name is ...", i5,'
      fmtn(3) = '20x, ".....  ", a16, /,'
      fmtn(4) = '"no. of dimensions and attributes are ... ", 2i5 )'
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncvinq", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid, vartyp, varnam,
     $        nvdims, nvatts
            if ( nvdims .gt. 0 ) then
               write ( nout(io), fmt2 ) (vdims(j),j=1,nvdims)
            endif
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 12. nwvg1i.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvg1i ( cdfid, vid, vindx,nd, val1, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, vindx(*), val1, nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncvg1i", rcode, nout )
      if (  nd .eq. 0 ) nd = 1
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgli: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the integer value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the integer is ...",i5 )'
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(i),i=1,nd), val1
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 13. nwvgt1.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvgt1 ( cdfid, vid, vindx,nd, val1, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, vindx(*), nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncvgt1", rcode, nout )
      if ( nd .eq. 0 ) nd = 1
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgt1: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the data value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the data is ...",1pe12.4 )'
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(i),i=1,nd), val1
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 14. nwvg1c.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvg1c ( cdfid, vid, vindx,nd, chval, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, vindx(*), nout(*)
      character chval
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvg1c: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the character value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the character is ... ",a )'
      if ( rcode .gt. 0 ) call errnc ( "ncvg1c", rcode, nout )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(j),j=1,nd), chval
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 15. nwgti.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvgti ( cdfid, vid, start, count,nd, val0, vals,
     $     rcode, nd1,nd2, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, start(*), count(*),
     $     val0(*), vals(nd1,nd2), nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncvgt1: cdfid = ", i3, " vid= ",i5,/,
     $     "start, first corner of hyperslab is ...",5i5,/,
     $     "count, edge lengths of the hyperslab is ...", 5i5 )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do i2 = 1, count(2)
         do i1 = 1, count(1)
            vals(i1,i2) = val0( i1+(i2-1)*count(1) ) 
         enddo
      enddo      
      if ( rcode .gt. 0 ) call errnc ( "ncvgt", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid,
     $           (start(i),i=1,5), (count(j),j=1,5)
         endif
      enddo
      call writis ( "Hyperslab of integers is..", vals, nd1,nd2,
     $     count, nd, nout )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 16. nwvgt.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvgt ( cdfid, vid, start,count,nd, val0, vals,
     $     nd1,nd2, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, start(*), count(*), nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
      dimension vals(nd1,nd2), val0(*)
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( nd .eq. 1 ) count(2) = 1
      do i2 = 1, count(2)
         do i1 = 1, count(1)
            vals(i1,i2) = val0( i1+(i2-1)*count(1) )
         enddo
      enddo
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgt: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"start, first corner of hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"count, edge lengths of the hyperslab is ...",'
     $     //fmt//')'
      if ( rcode .gt. 0 ) call errnc ( "ncvgt", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (start(i),i=1,nd), (count(j),j=1,nd)
         endif
      enddo
      call writvs ( "Hyperslab of Values is..", vals, nd1,nd2,
     $     count, nd, nout )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 17. nwvgtc.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwvgtc ( cdfid, vid, start,count,nd, string, lenstr,
     $     rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, start(*), count(*), nout(*)
      character*(*) string
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgtc: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"start, first corner of hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"count, edge lengths of the hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(4) = '"string, character string is ...",/, a )'
      if ( rcode .gt. 0 ) call errnc ( "ncvgtc", rcode, nout )
c-----------------------------------------------------------------------
c     write message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (start(i),i=1,nd),(count(j),j=1,nd), string
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 18. nwainq.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwainq ( cdfid, vid, attnam, attype, attlen, rcode,
     $     nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode, vid, attype, attlen, nout(*)
      character*(*) attnam
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncainq: cdfid = ", i3, ". Vid = ",i5,/,
     $     "attribute's name, type, and length: ",/,"... ", a,
     $     "... ", 2i5 )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "ncainq", rcode, nout ) 
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid, attnam,
     $           attype, attlen
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 19. nwagtc.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwagtc ( cdfid, vid, attnam, astrng, attlen, maxa1,
     $     rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode,vid, attlen, nout(*)
      character*(*) attnam, astrng
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 50   format ( /, "*** ", a, "'s length too long: ", i4, " >",i4 )
 100  format ( /,"ncagtc: cdfid =  ", i3,
     $     " attribute's Vid, name and value  = ", i5,/,
     $     "... ", a," ...",/,
     $     ">> ", a, " <<",/,
     $     "attribute length is ... ", i5 )
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      call vfrmt3 ( attlen, "a", fmt )
      fmt1 = '/,"ncagtc: cdfid = ", i3,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"attribute''s name, Vid, length, and value: ",/,'
      fmtn(3) = '"... ", a," ...", 2i5,/,'
      fmtn(4) = '">> ",'//fmt//' " <<" )'
      if ( rcode .gt. 0 ) call errnc ( "ncatgc", rcode, nout )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            if ( attlen .gt. maxa1 ) then
               write ( nout(io),50 ) attnam, attlen, maxa1
               return
            else
               write ( nout(io), fmtn ) cdfid, attnam, vid, attlen,
     $              astrng
            endif
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 20. nwagt.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwagt ( cdfid, vid, attnam, attval, attlen,
     $     rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode,vid, attlen, nout(*)
      character*(*) attnam
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncagt: cdfid =  ", i3,
     $     " attribute's id, name and value  = ", i5,/,
     $     "... ", a," ...",/,
     $     ">> ", 1pe12.4, " <<",/,
     $     "attribute length is ... ", i5 )
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      call vfrmt2 ( attlen, "f13.5", fmt )
      fmt1 = '/,"ncagt: cdfid = ", i3,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '" attribute''s name, Vid, length, and value: ",/,'
      fmtn(3) = '"... ", a," ...", 2i5,/,'
      fmtn(4) = '">> ",'//fmt//' " <<" )'
      if ( rcode .gt. 0 ) call errnc ( "ncatgc", rcode, nout )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
               write ( nout(io), fmtn ) cdfid, attnam, vid, attlen,
     $           attval
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 21. nwanam.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine nwanam ( cdfid, vid, n, attnam, rcode, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer cdfid, rcode,vid, nout(*)
      character*(*) attnam
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( /,"ncanam: cdfid = ", i3,
     $     ".  Attribute's variable ID is ... ", i5,/,
     $     "attribute's no. and name is ... ",
     $     i4," -- ", a )
c-----------------------------------------------------------------------
c     print messasge.
c-----------------------------------------------------------------------
      if ( rcode .gt. 0 ) call errnc ( "nwanam", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid, n, attnam
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 22. writn1.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writn1 ( mesage, n, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      integer nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( a, i5 ) 
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, n
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 23. writnn.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writnn ( mesage, n, nd, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      dimension n(*), nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( a, 5i5 ) 
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, (n(i), i=1,nd)
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 24. writa1.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writa1 ( mesage, string1, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      character*(*) string1
      integer nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( a, a ) 
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, string1
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 25. writv1.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writv1 ( mesage, a, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      integer nout(*)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 100  format ( a, 1pe12.4 ) 
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, a
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 26. writis.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writis ( mesage, a, nd1, nd2, nw, nd, nout )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      integer a(nd1,nd2), nout(*), nw(*)
c-----------------------------------------------------------------------
c     format statments.
c-----------------------------------------------------------------------
 10   format (/, a )
 50   format ( 10(1x,10i5,/) )
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 10 ) mesage
            do j1 = 1, nw(1)
               write ( nout(io), 50 ) ( a(j1,j2),j2 = 1, nw(2) )
            enddo
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 27. writvs.
c     prints message and numerical code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine writvs ( mesage, a, nd1, nd2, nw, nd, nout )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      dimension a(nd1,nd2), nw(*)
      integer nout(*)
c-----------------------------------------------------------------------
c     format statments.
c-----------------------------------------------------------------------
 10   format (/, a )
 50   format ( 10(1x,5e13.5,/) )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nw1 = nw(1)
      nw2 = nw(2)
      if ( nd .eq. 1 ) then
         nw1 = nw(2)
         nw2 = nw(1)
      endif
c-----------------------------------------------------------------------
c     print message.
c-----------------------------------------------------------------------
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), '(/,a)' ) mesage
            do j1 = 1, nw1
               if ( nd .ne. 1 ) then
                  write ( nout(io), 50 ) ( a(j1,j2),j2 = 1, nw2 )
               else
                  write ( nout(io), 50 ) ( a(j11,1), j11 = 1, nw2 )
               endif
            enddo
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 28. vfrmt1.
c     creates a format statement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vfrmt1 ( n, field1, fmt )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) fmt, field1
      character*8 nn
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   format ( i4 )
c-----------------------------------------------------------------------
c     create format statement.
c-----------------------------------------------------------------------
      write ( nn, 10 ) n
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
      fmt = '('//nn(1:index(nn, ' ')-1)//field1//')'
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 29. vfrmt2.
c     creates a format statement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vfrmt2 ( n, field1, fmt )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) fmt, field1
      character*8 nn
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   format ( i4 )
c-----------------------------------------------------------------------
c     create format statement.
c-----------------------------------------------------------------------
      write ( nn, 10 ) n
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
      fmt = nn(1:index(nn, ' ')-1)//field1
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 30. vfrmt3.
c     creates a format statement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vfrmt3 ( n, field1, fmt )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) fmt, field1
      character*8 nn
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   format ( i4 )
c-----------------------------------------------------------------------
c     create format statement.
c-----------------------------------------------------------------------
      write ( nn, 10 ) n
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
      fmt = field1//nn(1:index(nn, ' ')-1)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 31. bounds.
c     defines bounds.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine bounds ( x1,z1, n1,n2, xmin,xmax, zmin,zmax )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension x1(*), z1(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      xmax = -1.0e30
      zmax = -1.0e30
      xmin =  1.0e30
      zmin =  1.0e30
      do i = n1, n2
         if ( xmin .gt. x1(i) ) xmin = x1(i)
         if ( zmin .gt. z1(i) ) zmin = z1(i)
         if ( xmax .lt. x1(i) ) xmax = x1(i)
         if ( zmax .lt. z1(i) ) zmax = z1(i)
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 32. boundsi.
c     defines bounds.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine boundsi ( x1,z1, n1,n2, xmin,xmax, zmin,zmax,
     $     ixn,ixx,izn,izx )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension x1(*), z1(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      xmax = -1.0e30
      zmax = -1.0e30
      xmin =  1.0e30
      zmin =  1.0e30
      do i = n1, n2
         if ( xmin .gt. x1(i) ) then 
            xmin = x1(i)
            ixn = i
         endif
         if ( zmin .gt. z1(i) ) then
            zmin = z1(i)
            izn = i
         endif
         if ( xmax .lt. x1(i) ) then
            xmax = x1(i)
            ixx = i
         endif
         if ( zmax .lt. z1(i) ) then
            zmax = z1(i)
            izx = i
         endif
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 33. macopy.
c     copies a real array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine macopy ( ain,ndin, aout,ndout )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension ain(ndin,ndin), aout(ndout,ndout)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do i = 1, ndin
         do j = 1, ndin
            aout(i,j) = ain(i,j)
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 34. matwrt.
c     writes a real array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matwrt ( a, maxj1, maxj2, jmax1, jmax2, label )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) label
      integer out
      dimension a ( maxj1, maxj2 )
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   format ( //, 5x, "matrix elements of  ", a )
 20   format ( i4,10(1x,1p10e11.4,/) )
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
      out = 23
      write ( out, 10 ) label
      do j1 = 1, jmax1
         write ( out, 20 ) j1, ( a(j1,j2),j2 = 1, jmax2 )
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 35. vacasym.
c     writes a real array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vacasym ( vacin, nd,nnj,label, nout1,nout2 )
      USE local_mod, ONLY: r8

      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) label
      dimension vacin(nd,nd)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 111  format (/, 1x,"       Asymmetries in ", a )
 112  format(  /,1x," sum of diagonals, sd =            ",1pe12.4,
     $     /,1x," sum of off-diagonals, sod =       ",1pe12.4,
     $     /,1x," difference in off diagonals, dod =",1pe12.4,
     $     //,1x," dod/sod =                         ",1pe12.4,
     $     /,1x," dod/(sd+sod) =                    ",1pe12.4,/ )
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      vacrd = 0.0
      vacrso = 0.0
      vacrsd = 0.0
      write ( nout1, 111 ) label
      write ( nout2, 111 ) label
      do l1 = 1, nnj
         vacrsd = vacrsd + vacin(l1,l1)
         do l2 = 1, l1
            vacrd = vacrd + abs ( vacin(l1,l2) - vacin(l2,l1) )
            if ( l1 .ne. l2 )
     .           vacrso = vacrso + abs ( vacin(l1,l2)+vacin(l2,l1) )
         enddo
      enddo
      asymv1 = vacrd / vacrso
      asymv2 = vacrd / ( vacrso + vacrsd )
      write ( nout1,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
      write ( nout2,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 36. masym0.
c     writes a real array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine masym0 ( vacin, nd,nnj,n1,n2,n10,n20,
     $     label, nout1,nout2 )
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) label
      dimension vacin(nd,nd)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 111  format (/, 1x,"       Asymmetries in ", a,
     $     " -- Ignores la,lb = ", 2i3 )
 112  format(  /,1x," sum of diagonals, sd =            ",1pe12.4,
     $     /,1x," sum of off-diagonals, sod =       ",1pe12.4,
     $     /,1x," difference in off diagonals, dod =",1pe12.4,
     $     //,1x," dod/sod =                         ",1pe12.4,
     $     /,1x," dod/(sd+sod) =                    ",1pe12.4,/ )
c-----------------------------------------------------------------------
c     compuations.
c-----------------------------------------------------------------------
      vacrd = 0.0
      vacrso = 0.0
      vacrsd = 0.0
      write ( nout1, 111 ) label, n10, n20
      write ( nout2, 111 ) label, n10, n20
      do l1 = 1, nnj
         lla = l1 - 1 + n1
         if ( lla .ne. n10 ) then
            vacrsd = vacrsd + vacin(l1,l1)
            do l2 = 1, l1
               llb =  l2 - 1 + n2
               if ( llb .ne. n20 ) then
                  vacrd = vacrd + abs ( vacin(l1,l2) - vacin(l2,l1) )
                  if ( l1 .ne. l2 )
     $                 vacrso=vacrso+abs(vacin(l1,l2)+vacin(l2,l1))
               endif
            enddo
         endif
      enddo
      asymv1 = vacrd / vacrso
      asymv2 = vacrd / ( vacrso + vacrsd )
      write ( nout1,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
      write ( nout2,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 37. vacasymi.
c     writes a real array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine vacasymi ( vacin, nd,nnj,label, nout1,nout2 )
      USE local_mod, ONLY: r8
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      character*(*) label
      dimension vacin(nd,nd), wrkrd(nfm),wrkrso(nfm),wrkrsd(nfm),
     $     zsymv1(nfm),zsymv2(nfm)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 111  format (/, 1x,"  Asymmetries in nested sub-matrices of ", a )
 221  format(/,1x," sum of diagonals, sd =            ",/,(1p8e12.4))
 222  format(1x," sum of off-diagonals, sod =       ",/,(1p8e12.4))
 223  format(1x," difference in off diagonals, dod =",/,(1p8e12.4))
 224  format(1x," dod/sod =                         ",/,(1p8e12.4))
 225  format(1x," dod/(sd+sod) =                    ",/,(1p8e12.4))
c-----------------------------------------------------------------------
c     zero arrays.
c-----------------------------------------------------------------------
      do j = 1, nnj
         wrkrd(j) = 0.0
         wrkrso(j) = 0.0
         wrkrsd(j) = 0.0
      enddo
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      write ( nout1, 111 ) label
      write ( nout2, 111 ) label
      ncnt = 0
      do l1 = 1, nnj/2 + 1
         ja = l1
         jb = nnj - l1 + 1
         jab = jb - ja + 1
         if ( jab .ge. 2 )then
            do jin1 = 1, jab
               ja1 = ja + jin1 - 1
               wrkrsd(l1) = wrkrsd(l1) + vacin(ja1,ja1)
               do jin2 = 1, jin1
                  ja2 = ja + jin2 - 1
                  wrkrd(l1) = wrkrd(l1) +
     $                 abs ( vacin(ja1,ja2) - vacin(ja2,ja1) )
                  if ( ja1 .ne. ja2 )
     $                 wrkrso(l1) = wrkrso(l1) +
     $                 abs ( vacin(ja1,ja2)+vacin(ja2,ja1) )
               enddo
            enddo
            ncnt = ncnt + 1
         endif
      enddo
      do j = 1, ncnt
         zsymv1(j) = wrkrd(j) / wrkrso(j)
         zsymv2(j) = wrkrd(j) / ( wrkrso(j) + wrkrsd(j) )
      enddo
c-----------------------------------------------------------------------
c     write statements.
c-----------------------------------------------------------------------
      write ( nout1,'("ncnt = ",i3)') ncnt
      write ( nout2,'("ncnt = ",i3)') ncnt
      write ( nout1,221 ) ( wrkrsd(i), i = 1, ncnt )
      write ( nout1,222 ) ( wrkrso(i), i = 1, ncnt )
      write ( nout1,223 ) ( wrkrd(i), i = 1, ncnt )
      write ( nout1,224 ) ( zsymv1(i), i = 1, ncnt )
      write ( nout1,225 ) ( zsymv2(i), i = 1, ncnt )
      write ( nout2,221 ) ( wrkrsd(i), i = 1, ncnt )
      write ( nout2,222 ) ( wrkrso(i), i = 1, ncnt )
      write ( nout2,223 ) ( wrkrd(i), i = 1, ncnt )
      write ( nout2,224 ) ( zsymv1(i), i = 1, ncnt )
      write ( nout2,225 ) ( zsymv2(i), i = 1, ncnt )
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 38. msctimer.
c     computes elapsed cpu time; dummy routine.
c-----------------------------------------------------------------------
      subroutine msctimer ( nout,mesage )
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      character*(*) mesage
      return
      end
c-----------------------------------------------------------------------
c     subprogram 39. shftpi.
c     shifts angular variable by pi?
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine shftpi ( xin,fin, xout,fout, n )
      USE local_mod, ONLY: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension xin(*), fin(*), xout(*), fout(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nm = n-1
      n2 = nm/2
      do i = 1, n
         im = i-1 + n2
         ig = mod(im,nm) + 1
         xout(ig) = xin(i)
         if ( i .gt. n2 ) xout(ig) = xin(i) - xin(n)
         fout(ig) = fin(i)
      enddo
      fout(n) = fout(1)
      xout(n) = xin(n2+1)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      return
      end
