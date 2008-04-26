c-----------------------------------------------------------------------
c     file vacuum_penn.f.
c     local routines for penning.lanl.gov HP 735.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. date_time
c      2. clock
c      3. cleanup
c      4. shellb
c      5. gelima
c      6. gelimb
c      7. zop
c      8. zcl
c      9. zwr
c     10. zrd
c     11. skipeof
c     12. timedate
c     13. userinfo
c     14. close
c-----------------------------------------------------------------------
c     subprogram 1. date_time.
c-----------------------------------------------------------------------
      subroutine date_time(date_array,datex,timex)
      implicit real*8 (a-h,o-z)
      integer date_array(*)
      character*(*) datex,timex
      return
      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 2 clock.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine clock(ntim)
c$$$      implicit real*8 (a-h,o-z)
c$$$      character*(10) ntim
c$$$      return
c$$$      end
c-----------------------------------------------------------------------
c     subprogram 3. cleanup.
c-----------------------------------------------------------------------
      subroutine cleanup
      implicit real*8 (a-h,o-z)
      call system('rm -f modovmc')
      call system('rm -f pestotv')
      call system('rm -f vacdcon')
c      call system('rm -f ahg2msc.out')
      call system('rm -f mscvac.out')
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. shellb
c-----------------------------------------------------------------------
      subroutine shellb
      implicit real*8 (a-h,o-z)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 5. gelima.
c-----------------------------------------------------------------------
      subroutine gelima(copmat,nfm,uvpr,nfm1,jmax1,jmax2,uvp0,nfm2,
     $     wrki,waa,nfm3,wbb,nfm4,ifail)
      implicit real*8 (a-h,o-z)
      integer ipiv(jmax1),info
      call dgetrf(jmax1,jmax1,copmat,nfm,ipiv,info)
      call dgetrs('N',jmax1,jmax12evp0,copmat,nfm,ipiv,uvpr,nfm1,info)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 6. gelimb.
c-----------------------------------------------------------------------
      subroutine gelimb(copmat,nfm,uvpwr,nfm1,jmax1,jmax2,uvpw0,
     $     nfm2,wrki,ifail)
      implicit real*8 (a-h,o-z)
      integer ipiv(jmax1),info
      call dgetrf(jmax1,jmax1,copmat,nfm,ipiv,info)
      call dgetrs('N',jmax1,jmax2,copmat,nfm,ipiv,uvpw0,nfm2,info)
      return
      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 7. zop.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zop(ioc,name,nsize,idisk,icode,ilab)
c$$$      implicit real*8 (a-h,o-z)
c$$$      character*(8) name
c$$$      open(unit=ioc,file=name,form='unformatted',status='unknown')
c$$$      return
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 8. zcl.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zcl(ioc,ierr)
c$$$      implicit real*8 (a-h,o-z)
c$$$      close(ioc)
c$$$      return
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 9. zwr.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zwr(ioc,a,nwords,nadres,lgivup,irr)
c$$$      implicit real*8 (a-h,o-z)
c$$$      dimension a(1)
c$$$      integer*4 iocs,noffsets,fseek,ftell
c$$$      nbytes = 8
c$$$      iocs = ioc
c$$$      ncurr = ftell(iocs)
c$$$      ncurr = 0
c$$$      noffsets = nadres * nbytes - ncurr
c$$$      ierr = fseek(iocs,noffsets,0)
c$$$      write(ioc)(a(i),i=1,nwords)
c$$$      if(ierr .eq. 0) return
c$$$      write(6,100)ierr      
c$$$ 100  format(' Error in ZWR,error code = ',i4,
c$$$     $     ' Check man 3f perror ')
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 10. zrd.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zrd(ioc,a,nwords,nadres,lgivup,irr)
c$$$      implicit real*8 (a-h,o-z)
c$$$      dimension a(1)
c$$$      integer*4 iocs,noffsets,fseek,ftell
c$$$      nbytes = 8
c$$$      iocs = ioc
c$$$      ncurr = ftell(iocs)
c$$$      ncurr = 0
c$$$      noffsets = nadres * nbytes - ncurr
c$$$      ierr = fseek(iocs,noffsets,0)
c$$$      read(ioc)(a(i),i=1,nwords)
c$$$      if(ierr .eq. 0) return
c$$$      write(6,100)ierr
c$$$ 100  format(' Error in ZWR,error code = ',i4,
c$$$     $     ' Check man 3f perror ')
c$$$      stop
c$$$      end
c-----------------------------------------------------------------------
c     subprogram 11. skipeof.
c-----------------------------------------------------------------------
	subroutine skipeof(iva,iva1)
      implicit real*8 (a-h,o-z)
	return
	end
c-----------------------------------------------------------------------
c     subprogram 12. timedate.
c-----------------------------------------------------------------------
      subroutine timedate(ntim,ndat,mach,nsfx)
      implicit real*8 (a-h,o-z)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 13. userinfo.
c-----------------------------------------------------------------------
      subroutine userinfo(nuser,nacct,ndrop,nsfx)
      implicit real*8 (a-h,o-z)
      return
      end
c$$$c-----------------------------------------------------------------------
c$$$c     subprogram 14. close.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine close(iun)
c$$$      implicit real*8 (a-h,o-z)
c$$$      close(iun)
c$$$      return
c$$$      end
