c-----------------------------------------------------------------------
c     file vacuum_vms.f.
c     local routines for OpenVMS Alpha platforms
c     based on vacuum_penn.f 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. date_time
c      2. clock
c      3. shella
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
      call date_and_time(datex,timex)
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
c     subprogram 3. shella.
c-----------------------------------------------------------------------
      subroutine shella
      implicit real*8 (a-h,o-z)
      integer*4 lib$delete_file,status
      status=lib$delete_file('modovmc')
      status=lib$delete_file('pestotv')
      status=lib$delete_file('vdata')
      status=lib$delete_file('vac.cgm')
      status=lib$delete_file('vacadj')
      status=lib$delete_file('vacdcon')
      status=lib$delete_file('vacgato')
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
c-----------------------------------------------------------------------
c     subprogram 7. zop.
c-----------------------------------------------------------------------
      subroutine zop(ioc,name,nsize,idisk,icode,ilab)
c He apparently wants to control the writes directly. 
c VMS Fortran lets you do this, but it's a real kludge. 
c Since his reads and writes all pretend to be real*8
c I open the file for direct access with 8byte records
c (recl is in units of 4bytes)
      implicit real*8 (a-h,o-z)
      character*(8) name
      open(unit=ioc,file=name,
     .     form='unformatted',recl=2,
     .     status='unknown',
     .     access='direct')
      return
      end
c-----------------------------------------------------------------------
c     subprogram 8. zcl.
c-----------------------------------------------------------------------
      subroutine zcl(ioc,ierr)
      implicit real*8 (a-h,o-z)
      close(ioc)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 9. zwr.
c-----------------------------------------------------------------------
      subroutine zwr(ioc,a,nwords,nadres,lgivup,irr)
      implicit real*8 (a-h,o-z)
      dimension a(1)
      integer*4 iocs,noffsets,fseek,ftell
      nbytes = 1
      iocs = ioc
c      ncurr = ftell(iocs)
c Note the next line was present in the original, so the call to ftell did nothing?
      ncurr = 0
      noffsets = nadres * nbytes - ncurr
c      ierr = fseek(iocs,noffsets,0)
c Now do the writes one 8 byte record at a time.
c Probably adds a lot of overhead to the file in RMS. 
c I hope he doesn't need compatibility
      do i=1,nwords
         write(iocs,rec=noffsets,iostat=ierr) a(i)
         noffsets=noffsets+1
      enddo
      if(ierr .eq. 0) return
      write(6,100)ierr      
 100  format(' Error in ZWR,error code = ',i4,
     $     ' Check man 3f perror ')
      end
c-----------------------------------------------------------------------
c     subprogram 10. zrd.
c-----------------------------------------------------------------------
      subroutine zrd(ioc,a,nwords,nadres,lgivup,irr)
      implicit real*8 (a-h,o-z)
      dimension a(1)
      integer*4 iocs,noffsets,fseek,ftell
      nbytes = 1
      iocs = ioc
c      ncurr = ftell(iocs)
      ncurr = 0
      noffsets = nadres * nbytes - ncurr
c      ierr = fseek(iocs,noffsets,0)
c read the same way it's written, one 8 byte record at
c at time. 
      do i=1,nwords
         read(ioc,rec=noffsets,iostat=ierr) a(i)
         noffsets=noffsets+1
      enddo
      if(ierr .eq. 0) return
      write(6,100)ierr
 100  format(' Error in ZRD,error code = ',i4,
     $     ' Check man 3f perror ')
      return
      stop
      end
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
