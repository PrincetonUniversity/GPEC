c-----------------------------------------------------------------------
c     file vacuum_penn.f.
c     local routines for penning.lanl.gov HP 735.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. date_time
c      2. cleanup
c      3. shellb
c      4. gelima
c      5. gelimb
c      6. skipeof
c      7. timedate
c      8. userinfo
c-----------------------------------------------------------------------
c     subprogram 1. date_time.
c-----------------------------------------------------------------------
      subroutine date_time(date_array,datex,timex)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer date_array(*)
      character*(*) datex,timex
      return
      end
c-----------------------------------------------------------------------
c     subprogram 2. cleanup.
c-----------------------------------------------------------------------
      subroutine cleanup
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      call system('rm -f modovmc')
      call system('rm -f pestotv')
      call system('rm -f vacdcon')
      call system('rm -f mscvac.out')
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. shellb. 
c-----------------------------------------------------------------------
      subroutine shellb
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. gelima. presently unused.
c-----------------------------------------------------------------------
!       subroutine gelima(copmat,nfm,uvpr,nfm1,jmax1,jmax2,uvp0,nfm2,
!      $     wrki,waa,nfm3,wbb,nfm4,ifail)
!       use local_mod, only: r8
!       implicit real(r8) (a-h,o-z)
!       implicit integer (i-n)
!       integer ipiv(jmax1),info
!       call dgetrf(jmax1,jmax1,copmat,nfm,ipiv,info)
!       call dgetrs('N',jmax1,jmax12evp0,copmat,nfm,ipiv,uvpr,nfm1,info)
!       return
!       end
c-----------------------------------------------------------------------
c     subprogram 5. gelimb. Solves AX=B for A=copmat, B=uvpw0
c-----------------------------------------------------------------------
      subroutine gelimb(copmat,nfm,uvpwr,nfm1,jmax1,jmax2,uvpw0,
     $     nfm2,wrki,ifail)
      use local_mod, only: r8
      real(r8), dimension(jmax1,jmax1) :: copmat
      real(r8), dimension(nfm2,jmax2) :: uvpwr,uvpw0
      real(r8), dimension(*) :: wrki
      integer :: nfm, nfm1, nfm2, jmax1, jmax2, ifail, info
      integer, dimension(jmax1) :: ipiv
      call dgetrf(jmax1,jmax1,copmat,nfm,ipiv,info)
      call dgetrs('N',jmax1,jmax2,copmat,nfm,ipiv,uvpw0,nfm2,info)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 6. skipeof.
c-----------------------------------------------------------------------
	subroutine skipeof(iva,iva1)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
	return
	end
c-----------------------------------------------------------------------
c     subprogram 7. timedate.
c-----------------------------------------------------------------------
      subroutine timedate(ntim,ndat,mach,nsfx)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      return
      end
c-----------------------------------------------------------------------
c     subprogram 8. userinfo.
c-----------------------------------------------------------------------
      subroutine userinfo(nuser,nacct,ndrop,nsfx)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      return
      end
