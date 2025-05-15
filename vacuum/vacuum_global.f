c-----------------------------------------------------------------------
c     file vacuum_alloc.f.
c     dynamically allocates arrays.
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. vglobal_mod
c     1. global_alloc.
c     2. global_dealloc.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. vglobal_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE vglobal_mod
      USE local_mod, ONLY: r8
      IMPLICIT NONE

      LOGICAL :: lsymz,check1,check2,lanal,lkdis,lpest1,lpless,wall,
     $     lnova,checks,lfunin,checke,checkd,symvac,verbose_timer_output

      INTEGER, PARAMETER :: nc31=31,maxc1=100,maxa1=200,max2=3,
     $     nccl3=72,numvar=100,ndima=2,ndim0=5,neqv1=1

      CHARACTER(8) :: jobid
      CHARACTER(60) :: seps
      CHARACTER(128) :: ahgdir = "."
      CHARACTER(128) :: ahgfile = "ahg2msc.out"
      CHARACTER(nccl3) :: dskout,cdfin,cdfout
      CHARACTER(maxc1) :: ctitle
      CHARACTER(maxa1) :: astrng
      CHARACTER(nc31) :: vname, vrname(numvar), recnam, attnam

      INTEGER :: nsf0=0,nths0,nfm,mtot,ntsin0,ntsin,nsf,nfe,nths,nths2,
     $     nfm2,nwrkl,nthsq,nthsq4,nfnt4,nfm21,outpest,outmod,outmap1,
     $     idsk,iloop,inmode,inpest,intty,iodsk,iomode,iotty,
     $     iovac,ipshp,ladj,ldcon,leqarcw,lfele,lfour,lgato,lj,
     $     lkplt,lnsav,lrgato,lspark,lxsav,lzio,m,mdiv,mfel,mj,mp,mp0,
     $     mp1,mth,mth1,mth2,mthin,mthin1,mthin2,ndfel,neigvc,neigvl,
     $     nj,nosurf,nunst,nx,nz,nzd1,nzd2,idgt,idot,ieig,ieps,ishape,
     $     ismth,lff,lmax1,lpsub,lwrt11,mphi,mx,mz,nloop,nloopr,nminus,
     $     noutv,nph,nphil,nphse,nplus,nsing,ntloop,cdfid,nd1,nd2,nd12,
     $     neqv2,neqv3,neqv4,neqv5,nfmsq,ndimlp,
     $     nxlpin,nzlpin,lgpec,linterior,ixx,ixn,izx

      INTEGER, DIMENSION(3) :: nout,nout0
      INTEGER, DIMENSION(5) :: nw
      INTEGER, DIMENSION(20) :: ntitle
      INTEGER, DIMENSION(:), POINTER :: lmax,lmin,lfm

      REAL(r8) :: alx,alz,betag,betai,betap,betav,bit,ctroy,dat,date0,
     $     datem,datev,datime,dth,dthin,eqdr,eqdt,eqpi,fa1,five,
     $     four,ga1,gp0,half,one,p0,psilim,psimin,psipls,pye,qa1,
     $     r2,r4,r6,rgato,seven,three,time0,timem,timev,two,twopi,
     $     twopi2,upsiln,x000,xma,xzero,z47,z48,z49,z54,z55,z56,z57,
     $     zero,zma,n,no2pi,no2pi2,aval0
      REAL(r8) :: a,abulg,aloop,apl,aval,aw,b,bbulg,bloop,bpl,bval,bw,
     $     civ,cn0,cw,delfac,delg,deloop,delx,delz,dloop,dpl,dw,epsq,
     $     ff,qain,rloop,sp2sgn1,sp2sgn2,sp2sgn3,sp2sgn4,sp2sgn5,
     $     sp3sgn1,sp3sgn2,sp3sgn3,sp3sgn4,sp3sgn5,tbulg,tw,xofsl,xpl,
     $     xs,xt,xtp,zs,zt,ztp,
     $     epslp,xlpmin,xlpmax,zlpmin,zlpmax,plrad,xmaj,zmxp

      REAL(r8) :: kernelsign
      REAL(r8), DIMENSION(numvar) :: vn,vrnat,vrtype
      REAL(r8), DIMENSION(numvar,ndima) :: vvdims
      REAL(r8), DIMENSION(ndima) :: count,start,vdims,vindx
      REAL(r8), DIMENSION(ndim0) :: dimsiz

      REAL(r8), DIMENSION(9) :: xiin=(/0,0,0,0,0,0,0,1,0/)
      REAL(r8), DIMENSION(100) :: eigval
      REAL(r8), DIMENSION(:), POINTER :: xirc,xirs,xiic,xiis,
     $     grpssq,xsq,gpsdth,xsqdth,xjacob,delta,xjdtxj,xsdtxs,
     $     gpdtgp,slngth,xinf,zinf,xplap,zplap,fv,val0,
     $     xobp,zobp,xloop,zloop
      REAL(r8), DIMENSION(:,:), POINTER :: vacmat,vacmatu,vacmtiu,vals

      REAL(r8), DIMENSION(:), POINTER :: xwal,zwal,xwalp,zwalp,cnqd,snqd
      REAL(r8), DIMENSION(:,:), POINTER :: coslt,cplai,cplar,cpwi,cpwr,
     $     cslth,cwali,cwalli,cwallr,cwalr,cwpti,cwptr,cwwti,cwwtr,
     $     sinlt,snlth,wrktr1,wrktr2,gatovac

      LOGICAL :: farwal = .false.
      INTEGER :: mw,jtop,jbot
      REAL(r8), DIMENSION(:), POINTER :: xpass,zpass,xpla,zpla
      REAL(r8), DIMENSION(:,:), POINTER :: grwp,grri,chiwc,chiws

      REAL(r8), DIMENSION(:), POINTER :: bxpwtr,bxpwti,bzpwtr,bzpwti,
     $     bnpwtr,bnpwti,bnlr,bnli,xilr,xili,xigr,xigi,bnkr,bnki,
     $     bthpr,bthpi,gpsjp,bnpptr,bnppti,chipr,chipi
      REAL(r8), DIMENSION(:,:), POINTER :: bxpwr,bxpwi,bzpwr,bzpwi,
     $     bnpwr,bnpwi

      ! dcon interface
      LOGICAL :: dcon_set = .FALSE.
      INTEGER :: mth_dcon,lmin_dcon,lmax_dcon,nn_dcon
      REAL(r8) :: qa1_dcon
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x_dcon, z_dcon, delta_dcon

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. global_alloc.
c     dynamically allocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE global_alloc(nths0,nfm,mtot,ntsin0)

      INTEGER, INTENT(IN) :: nths0,nfm,mtot,ntsin0
c-----------------------------------------------------------------------
c     define derived sizes.
c-----------------------------------------------------------------------
      ntsin=ntsin0+5
      nsf=nsf0+1
      nfe=1+nsf/2
      nths=nths0+5 !!
      nths2=2*nths
      nfm2=2*nfm
      nwrkl=nths2**2+3*nths2
      nthsq=nths*nths
      nthsq4=4*nthsq
      nfnt4=4*nths*nfm
      nfm21=2*nfm+1
      nd1=ntsin
      nd2=nsf
      nd12=nd1*nd2
      neqv2=(nfm-1)/2+1
      neqv3=nfm+1
      neqv4=3*(nfm+1)/2
      neqv5=2*nfm+1
      nfmsq=nfm*nfm
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(lmax(nfe),lmin(nfe))
      ALLOCATE(xirc(nfm),xirs(nfm),xiic(nfm),xiis(nfm),fv(nfm))
      ALLOCATE(grpssq(nths),xsq(nths),gpsdth(nths),xsqdth(nths),
     $     xjacob(nths),delta(nths),xjdtxj(nths),xsdtxs(nths),
     $     gpdtgp(nths),slngth(nths),xinf(nths),zinf(nths),xplap(nths),
     $     zplap(nths))
      ALLOCATE(vacmat(nfm,nfm),vacmatu(nfm,nfm),vacmtiu(nfm,nfm))
      ALLOCATE(val0(nd12),vals(nd1,nd2))
      ALLOCATE(xwal(nths),zwal(nths),xwalp(nths),zwalp(nths),cnqd(nths)
     $     ,snqd(nths))
      ALLOCATE(cslth(nths,nfm),snlth(nths,nfm),coslt(nths,nfm),
     $     sinlt(nths,nfm),cplar(nths,nfm),cplai(nths,nfm),
     $     cwallr(nths,nfm),cwalli(nths,nfm),cwalr(nths,nfm),
     $     cwali(nths,nfm),cpwr(nths,nfm),cpwi(nths,nfm),
     $     cwwtr(nths,nfm),cwwti(nths,nfm),cwptr(nths,nfm),
     $     cwpti(nths,nfm),wrktr1(nfm,nths),wrktr2(nfm,nths))
      ALLOCATE(xpass(nths),zpass(nths),chiwc(nths,nfm),chiws(nths,nfm))
      ALLOCATE(grwp(nths,nths),grri(nths2,nfm2))
      ALLOCATE(gatovac(nfm,nfm))
      ALLOCATE(bxpwr(nths,nfm),bxpwi(nths,nfm),bzpwr(nths,nfm),
     $     bzpwi(nths,nfm),bnpwr(nths,nfm),bnpwi(nths,nfm))
      ALLOCATE(bxpwtr(nths),bxpwti(nths),bzpwtr(nths),bzpwti(nths),
     $     bnpwtr(nths),bnpwti(nths))
      ALLOCATE(xilr(nths),xili(nths),xigr(nths),xigi(nths),
     $     bnkr(nths),bnki(nths),bthpr(nths),bthpi(nths),
     $     gpsjp(nths),bnpptr(nths),bnppti(nths),
     $     chipr(nths),chipi(nths))
      ALLOCATE(lfm(nfm),bnlr(nfm),bnli(nfm))
c-----------------------------------------------------------------------
c     zero arrays.
c-----------------------------------------------------------------------
      lmax=0
      lmin=0
      xirc=0
      xirs=0
      xiic=0
      xiis=0
      fv=0
      grpssq=0
      xsq=0
      gpsdth=0
      xsqdth=0
      xjacob=0
      delta=0
      xjdtxj=0
      xsdtxs=0
      gpdtgp=0
      slngth=0
      xinf=0
      zinf=0
      xplap=0
      zplap=0
      vacmat=0
      vacmatu=0
      vacmtiu=0
      val0=0
      vals=0
      xwal=0
      zwal=0
      xwalp=0
      zwalp=0
      cnqd=0
      snqd=0
      cslth=0
      snlth=0
      coslt=0
      sinlt=0
      cplar=0
      cplai=0
      cwallr=0
      cwalli=0
      cwalr=0
      cwali=0
      cpwr=0
      cpwi=0
      cwwtr=0
      cwwti=0
      cwptr=0
      cwpti=0
      wrktr1=0
      wrktr2=0
      xpass=0
      zpass=0
      chiwc=0
      chiws=0
      grwp=0
      grri=0
      gatovac=0
      bxpwr=0
      bxpwi=0
      bzpwr=0
      bzpwi=0
      bnpwr=0
      bnpwi=0
      bxpwtr=0
      bxpwti=0
      bzpwtr=0
      bzpwti=0
      bnpwtr=0
      bnpwti=0
      bnlr=0
      bnli=0
c-----------------------------------------------------------------------
c     set pointers.
c-----------------------------------------------------------------------
      xpla => xinf
      zpla => zinf
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE global_alloc
c-----------------------------------------------------------------------
c     subprogram 2. global_dealloc.
c     deallocates arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE global_dealloc
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(lmax,lmin)
      DEALLOCATE(xirc,xirs,xiic,xiis,fv)
      DEALLOCATE(grpssq,xsq,gpsdth,xsqdth,xjacob,delta,xjdtxj,xsdtxs,
     $     gpdtgp,slngth,xplap,zplap,xinf,zinf)
      DEALLOCATE(vacmat,vacmatu,vacmtiu)
      DEALLOCATE(val0,vals)
      DEALLOCATE(xwal,zwal,xwalp,zwalp,cnqd,snqd)
      DEALLOCATE(cslth,snlth,coslt,sinlt,cplar,cplai,cwallr,cwalli,
     $     cwalr,cwali,cpwr,cpwi,cwwtr,cwwti,cwptr,cwpti,wrktr1,wrktr2)
      DEALLOCATE(xpass,zpass,chiwc,chiws)
c-----------------------------------------------------------------------
c     modification.
c-----------------------------------------------------------------------
      DEALLOCATE(grwp,grri)
c-----------------------------------------------------------------------
      DEALLOCATE(gatovac)
      DEALLOCATE(bxpwr,bxpwi,bzpwr,bzpwi,bnpwr,bnpwi)
      DEALLOCATE(bxpwtr,bxpwti,bzpwtr,bzpwti,bnpwtr,bnpwti)
      DEALLOCATE(xilr,xili,xigr,xigi,bnkr,bnki,bthpr,bthpi,
     $     gpsjp,bnpptr,bnppti,chipr,chipi)     
      DEALLOCATE(lfm,bnlr,bnli)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE global_dealloc
      END MODULE vglobal_mod
