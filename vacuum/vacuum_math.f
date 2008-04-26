c-----------------------------------------------------------------------
c     file vacuum_math.f.
c     mathematical subroutines used by Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. spl1d1
c      2. spl1d2
c      3. labrt
c      4. search
c      5. searchx
c      6. green
c      7. aleg
c      8. trans
c      9. transdx
c     10. smooth0
c     11. smooth
c     12. lagp
c     13. shft
c     14. lagpe4
c     15. lag
c     16. eigen
c     17. mult
c     18. matmul1
c     19. matmul3
c     20. indef4
c-----------------------------------------------------------------------
c     subprogram 1. spl1d1.
c     spline fitting routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine spl1d1(n,x,f,w,iop,ij,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension iop(*),x(*),f(*),w(*),a(*),b(*),c(*)
      real*8 comm(6)
      data  comm          /8hspl1d1 n,8h less th,8han 4. re,8hsults in,
     1 8hcorrect.,8h        /
      data zz,oz,tz,sz/0.0e0,1.0e0,3.0e0,6.0e0/
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      k=n-1
      a(2)=-(x(2)-x(1))/sz
      b(2)=(x(3)-x(1))/tz
      w(ij+1)=(f(2*ij+1)-f(ij+1))/(x(3)-x(2))-(f(ij+1)-f(1))
     1/(x(2)-x(1))
      if (n-3)3,4,3
    3 do 10 i=3,k
      m=(i-1)*ij+1
      j1=m+ij
      j2=m-ij
      con=(x(i+1)-x(i-1))/tz
      don=(x(i)-x(i-1))/sz
      b(i)=con-(don**2)/b(i-1)
      e=(f(j1)-f(m))/(x(i+1)-x(i))-(f(m)-f(j2))/
     1(x(i)-x(i-1))
      w(m)=e-(don*w(j2))/b(i-1)
   10 a(i)=-(don*a(i-1))/b(i-1)
    4 k1=(n-2)*ij+1
      c(n-1)=-((x(n)-x(n-1))/sz)/b(n-1)
      w(k1)=w(k1)/b(n-1)
      a(n-1)=a(n-1)/b(n-1)
      k2=k-1
      if (n-3)7,8,7
    7 do 20 i=2,k2
      j=n-i
      con=(x(j+1)-x(j))/sz
      a(j)=(a(j)-con*a(j+1))/b(j)
      c(j)=-(con*c(j+1))/b(j)
      k3=(j-1)*ij+1
      m=k3+ij
   20 w(k3)=(w(k3)-con*w(m))/b(j)
    8 k4=(n-1)*ij+1
      if (iop(1)-5) 201,200,201
  201 c1=w(1)
      if (iop(2)-5) 203,202,203
  203 c2=w(k4)
      go to 205
  200 if (n-4)300,302,302
  302 a1=x(1)-x(2)
      a2=x(1)-x(3)
      a3=x(1)-x(4)
      a4=x(2)-x(3)
      a5=x(2)-x(4)
      a6=x(3)-x(4)
      w(1)=f(1)*(oz/a1+oz/a2+oz/a3)-a2*a3*f(ij+1)/(a1*a4*a5)+
     1 a1*a3*f(2*ij+1)/(a2*a4*a6 )-a1*a2*f(3*ij+1)/(a3*a5*a6)
      go to 201
  202 if (n-4)300,303,303
  303 b1=x(n)-x(n-3)
      b2=x(n)-x(n-2)
      b3=x(n)-x(n-1)
      b4=x(n-1)-x(n-3)
      b5=x(n-1)-x(n-2)
      b6=x(n-2)-x(n-3)
      l1=k4-ij
      l2=l1-ij
      l3=l2-ij
      w(k4)=-b2*b3*f(l3)/(b6*b4*b1)+b1*b3*f(l2)/(b6*b5*b2)
     1 -b1*b2*f(l1)/(b4*b5*b3)+f(k4)*(oz/b1+oz/b2+oz/b3)
      go to 203
 205          i    =    1
 2051 continue
      m=(i-1)*ij+1
      go to 60
   70 if (i-1)80,50,80
   80 w(1)=w(1)-bob*w(m)
      w(k4)=w(k4)-bill*w(m)
      a(1)=a(1)-bob*a(i)
      a(n)=a(n)-bill*a(i)
      c(1)=c(1)-bob*c(i)
      c(n)=c(n)-bill*c(i)
   50 continue
      i=i+1
      if ( i .le. k )   go to 2051
      go to 100
   60 mk=iop(1)
      go to (62,64,66,68,66),mk
   62 if (i-1)71,63,71
   63 a(1)=-oz
      c(1)=zz
      go to 500
   71 bob=zz
      go to 500
   64 if (i-1)73,76,73
   76 a(1)=-oz
      c(1)=zz
      w(1)=zz
      go to 500
   73 if (i-2)81,81,82
   81 bob=-c1
      go to 500
   82 bob=zz
      go to 500
   66 if (i-1)83,84,83
   84 a(1)=-(x(2)-x(1))/tz
      c(1)=zz
      w(1)=-c1+(f(ij+1)-f(1))/(x(2)-x(1))
      go to 500
   83 if (i-2)85,85,86
   85 bob=(x(2)-x(1))/sz
      go to 500
   86 bob=zz
      go to 500
   68 if (i-1)87,88,87
   88 a(1)=-oz
      c(1)=oz
      w(1)=zz
      go to 500
   87 bob=zz
  500 ml=iop(2)
      go to (120,130,140,150,140),ml
  120 if (i-1)121,122,121
  122 a(n)=zz
      c(n)=-oz
      go to 70
  121 bill=zz
      go to 70
  130 if (i-1)131,132,131
  132 a(n)=zz
      c(n)=-oz
      w(k4)=zz
      go to 70
  131 if (i-k)134,133,134
  133 bill=-c2
      go to 70
  134 bill=zz
      go to 70
  140 if (i-1)141,142,141
  142 a(n)=zz
      c(n)=(x(n-1)-x(n))/tz
      w(k4)=c2-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  141 if (i-k)143,144,143
  144 bill=(x(n)-x(n-1))/sz
      go to 70
  143 bill=zz
      go to 70
  150 if (i-1)151,152,151
  152 a(n)=zz
      c(n)=(x(n-1)+x(1)-x(n)-x(2))/tz
      w(k4)=(f(ij+1)-f(1))/(x(2)-x(1))-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  151 if (i-2)153,154,153
  154 bill=(x(2)-x(1))/sz
      go to 70
  153 if (i-k)155,156,155
  156 bill=(x(n)-x(n-1))/sz
      go to 70
  155 bill=zz
      go to 70
  100 con=a(1)*c(n)-c(1)*a(n)
      d1=-w(1)
      d2=-w(k4)
      w(1)=(d1*c(n)-c(1)*d2)/con
      w(k4)=(a(1)*d2-d1*a(n))/con
      do 110 i=2,k
      m=(i-1)*ij+1
  110 w(m)=w(m)+a(i)*w(1)+c(i)*w(k4)
      go to 305
  300 call labrt(1,comm,1)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
  305 return
      end
c-----------------------------------------------------------------------
c     subprogram 2. spl1d2.
c     spline evaluation routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine spl1d2(n,x,f,w,ij,y,tab)
      implicit real*8 (a-h,o-z)
      data wz,sz/2.0e0,6.0e0/
      dimension x(3),f(3),w(3),tab(3)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mflag = 0
      if(y-x(1))10,10,20
   10 i=1
      go to 30
   20 if(y-x(n))15,40,40
   40 i=n-1
      go to 30
   15 call search(y,x,n,i,mflag)
   30 mi=(i-1)*ij+1
      k1=mi+ij
      flk=x(i+1)-x(i)
      a=(w(mi)*(x(i+1)-y)**3+w(k1)*(y-x(i))**3)/(sz*flk)
      b=(f(k1)/flk-w(k1)*flk/sz)*(y-x(i))
      c=(f(mi)/flk-flk*w(mi)/sz)*(x(i+1)-y)
      tab(1)=a+b+c
      a=(w(k1)*(y-x(i))**2-w(mi)*(x(i+1)-y)**2)/(wz*flk)
      b=(f(k1)-f(mi))/flk
      c=flk*(w(mi)-w(k1))/sz
      tab(2)=a+b+c
      tab(3)=(w(mi)*(x(i+1)-y)+w(k1)*(y-x(i)))/flk
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 3. labrt.
c     incomprehensible spaghetti code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine labrt(isw,lhol,inx)
      implicit real*8 (a-h,o-z)

      real*8 lhol(8)
      logical ps, ts
      data np/10/,ps/.true./,ts/.false./
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 27   format(1h0,9x,8a10,3x,z4)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if((isw.eq.0).or.(isw.gt.5))return
      go to ( 1,2,3,4,5 ), isw
    1 if ( ps .and. (np .gt. 0) )   write ( 3, 27 )   lhol, inx
      np=np-1
      if ( ts )   stop
      return
    2 ps=.false.
      return
    3 ps=.true.
      np=inx
      return
    4 ts=.true.
      return
    5 ts=.false.
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 4. search.
c     finds cubic spline interval.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine search(xbar,x,n,i,mflag)
      implicit real*8 (a-h,o-z)
      real*8 com1(5)
      dimension x(1)
      data com1/ 'search  ','xbar is ','outside ','range of',' table'/
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mflag=0
      i=n
      if(xbar.eq.x(n))return
      i=1
      if(n.le.1) return
      ixbar=sign(1.0d0,xbar)
      ix1=sign(1.0d0,x(1))
      ixn=sign(1.0d0,x(n))
      do 5 k=1,n
      j=i+i
      if(j.ge.n) go to 6
  5   i=j
 6    k=i
      mflag = 1
      do 115  l=2,n
      a=x(l-1)
      b=x(l)
      if(sign(1.0d0,a)-sign(1.0d0,b)) 7,113,8
 113   if(a-b)7,115,8
 115    continue
  7    j=1
      if(ixbar.lt.ix1.or.(ixbar.eq.ix1.and.xbar.lt.x(1)).or.ixbar
     1.gt.ixn.or.(ixbar.eq.ixn.and.xbar.gt.x(n))) go to 16
      go to 10
  8    j=2
      if(ixbar.lt.ixn.or.(ixbar.eq.ixn.and.xbar.lt.x(n)).or.ixbar
     1.gt.ix1.or.(ixbar.eq.ix1.and.xbar.gt.x(1))) go to 16
   10 k=k/2
       a=x(i)
      go to (11,20),j
  11   if(ixbar-sign(1.0d0,a)) 111,1111,2111
 1111 if(xbar-a)111,14,2111
 2111 b=x(i+1)
      if(ixbar-sign(1.0d0,b)) 2112,2113,12
 2113   if(xbar.ge.b) go to 12
 2112   return
 111  i = i-k
      go to 13
   12 i = i+k
   13 if (i.lt.n) go to 10
      k=k/2
      go to 111
   14 mflag=0
      return
   16 call labrt(1,com1,1)
      mflag=2
      return
  20   if(ixbar-sign(1.0d0,a) ) 2120,2121,111
 2121  if(xbar-a) 2120,14,111
 2120 b=x(i+1)
       if(ixbar-sign(1.0d0,b)) 12,2122,2112
 2122  if(xbar-b) 12,12,2112
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
 8888    return
      end
c-----------------------------------------------------------------------
c     subprogram 5. searchx.
c     another routine for finding cubic spline interval; never used.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine searchx(xbar,x,n,i,mflag)
      implicit real*8 (a-h,o-z)
      real*8 com1(5)
      integer xbar,x(1)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if(.false.) go to 8888
      data com1/ "search  ","xbar is ","outside ","range of"," table"/
      mflag=0
      i=n
      if(xbar.eq.x(n))return
      i=1
      if(n.le.1) return
      ixbar=isign(1,xbar)
      ix1=isign(1,x(1))
      ixn=isign(1,x(n))
      do 5 k=1,n
      j=i+i
      if(j.ge.n) go to 6
  5   i=j
 6    k=i
      mflag = 1
      do 115  l=2,n
      ia=x(l-1)
      ib=x(l)
      if(isign(1,ia)-isign(1,ib)) 7,113,8
 113   if(ia-ib)7,115,8
 115    continue
  7    j=1
      if(ixbar.lt.ix1.or.(ixbar.eq.ix1.and.xbar.lt.x(1)).or.ixbar
     1.gt.ixn.or.(ixbar.eq.ixn.and.xbar.gt.x(n))) go to 16
      go to 10
  8    j=2
      if(ixbar.lt.ixn.or.(ixbar.eq.ixn.and.xbar.lt.x(n)).or.ixbar
     1.gt.ix1.or.(ixbar.eq.ix1.and.xbar.gt.x(1))) go to 16
   10 k=k/2
       ia=x(i)
      go to (11,20),j
  11   if(ixbar-isign(1,ia)) 111,1111,2111
 1111 if(xbar-ia)111,14,2111
 2111 ib=x(i+1)
      if(ixbar-isign(1,ib)) 2112,2113,12
 2113   if(xbar.ge.ib) go to 12
 2112   return
 111  i = i-k
      go to 13
   12 i = i+k
   13 if (i.lt.n) go to 10
      k=k/2
      go to 111
   14 mflag=0
      return
   16 call labrt(1,com1,1)
      mflag=2
      return
  20   if(ixbar-isign(1,ia) ) 2120,2121,111
 2121  if(xbar-ia) 2120,14,111
 2120 ib=x(i+1)
       if(ixbar-isign(1,ib)) 12,2122,2112
 2122  if(xbar-ib) 12,12,2112
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
 8888    return
      end
c-----------------------------------------------------------------------
c     subprogram 6. green.
c     computes Green's function from Legendre functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine green
      USE vglobal_mod
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      pm = 0.
      pn = 0.
      pp = 0.
      aleg0 = 0.
      aleg1 = 0.
      sqpi  = sqrt(pye)
      pii  = two / pye
      gam  = sqpi
      nloc  = n + 0.1e0
      xs2  = xs * xs
      xt2  = xt * xt
      xp2  = xs2 + xt2
      xm2  = xt2 - xs2
      zm2  = ( zt - zs )**2
      r14  = xm2*xm2 + zm2*zm2 + two*xp2*zm2
      r1sq = sqrt( r14 )
      r1 = sqrt( r1sq )
      s  = (xp2 + zm2 )/r1sq
      call aleg ( s,nloc, pm,pn,pp, aleg0,aleg1 )
      kloc=0
      ak=zero
      if ( nloc .eq. 0 )  go to 10
    5 kloc = kloc+1
      ak = float(kloc)
      ak02 = half-ak
      gam = gam / ak02
      if ( kloc .ne. nloc )  go to 5
   10 gg  = -two * sqpi * gam / r1
      bval  = -gg*pn
      aval1 = ( n*(xs2+xt2+zm2)*(xt2-xs2-zm2)+xt2*(xm2+zm2))*pn
      aval2 = two*xt*xs*(xm2-zm2)*pp
      aval3 = ztp*(aval1+aval2) / xt
      aval4 = ( two*n+one)*(xp2+zm2)*pn+four*xt*xs*pp
      aval5 = xtp*(zt-zs)*aval4
      aval6 =(aval3-aval5) / ( xt*r14 )
      aval = - xt2*aval6 * gg / twopi
      aval0 = ztp*(two*xs*(zm2-xm2)*aleg1 - xt*(xm2+zm2)*aleg0)
      aval0 = aval0 + xtp*(zt-zs)*(four*xt*xs*aleg1+(xp2+zm2)*aleg0)
      aval0 = -aval0*xt / (r14*r1)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 7. aleg.
c     computes Legendre functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE aleg(x,nloc,pm,pn,pp, aleg0,aleg1)
      IMPLICIT NONE

      REAL(8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: nloc
      REAL(8), INTENT(OUT) :: pm,pn,pp,aleg0,aleg1

      INTEGER :: kloc
      REAL(8) :: ak02,elipe,elipk,gam,s,v,w,x1,x2,x3,x4,xxq,y,ysq

      REAL(8), PARAMETER :: pi=3.1415926535897931_8,pii=2/pi,
     $     sqpi=1.7724538509055159_8
      REAL(8), PARAMETER :: 
     $     ak0=1.38629436112_8,
     $     ak1=0.09666344259_8,
     $     ak2=0.03590092383_8,
     $     ak3=0.03742563713_8,
     $     ak4=0.01451196212_8,
     $     bk0=0.5_8,
     $     bk1=0.12498593597_8,
     $     bk2=0.06880248576_8,
     $     bk3=0.03328355346_8,
     $     bk4=0.00441787012_8,
     $     ae1=0.44325141463_8,
     $     ae2=0.0626060122_8,
     $     ae3=0.04757383546_8,
     $     ae4=0.01736506451_8,
     $     be1=0.2499836831_8,
     $     be2=0.09200180037_8,
     $     be3=0.04069697526_8,
     $     be4=0.00526449639_8
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      gam=sqpi
      xxq=x*x
      s=(x-1)/(x+1)
      ysq=xxq-1
      y=SQRT(ysq)
      w=x+y
      v=2*y/w
      x1=2/(x+1)
      x2=x1*x1
      x3=x2*x1
      x4=x3*x1
      elipk=ak0+ak1*x1+ak2*x2+ak3*x3+ak4*x4
     $     -(bk0+bk1*x1+bk2*x2+bk3*x3+bk4*x4)*dlog(x1)
      pn=pii*SQRT(2/(x+1))*elipk
      aleg0=pn
      x1=1/w**2
      x2=x1*x1
      x3=x2*x1
      x4=x3*x1
      elipe=1
      IF(ABS(x1)  >  1e-6)
     $     elipe=1+ae1*x1+ae2*x2+ae3*x3+ae4*x4
     $     -(be1*x1+be2*x2+be3*x3+be4*x4)*dlog(x1)
      pp=(pii*SQRT(w)*elipe-x*pn)/(2*y)
      aleg1=pp
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO kloc=1,nloc
         ak02=.5_8-kloc
         pm=pn
         pn=pp
         pp=-2*kloc*x*pn/y-ak02*ak02*pm
         gam=gam/ak02
      ENDDO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 8. trans.
c     transforms (translates?) something.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine trans ( vecin,mthin, vecout,mth )
      implicit real*8 (a-h,o-z)
      dimension vecin(*), vecout(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
      if ( mth .eq. mthin )then
         do i = 1, mth + 2
            vecout(i) = vecin(i)
         enddo
         return
      endif
      do i = 1, mth
         ai = i-1
         x = ai / mth
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
      enddo
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 9. transdx.
c     transforms (translates?) something.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine transdx ( vecin,mthin, vecout,mth, dx0 )
      implicit real*8 (a-h,o-z)
      dimension vecin(*), vecout(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
      if ( (mth .eq. mthin) .and. (abs(dx0) .le. 1.e-6) ) then
         do i = 1, mth + 2
            vecout(i) = vecin(i)
         enddo
         return
      endif
      do i = 1, mth
         ai = i-1
         x = ai / mth + dx0 / mthin
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
      enddo
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 10. smooth0.
c     smooth an array with a moving average.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine smooth0 ( g, n )
      implicit real*8 (a-h,o-z)
      dimension g(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nm = n - 1
      do j = 1, nm
         g(j) = g(j) + g(j+1)
      enddo
      do j = 1, nm
         g(n+1-j) = g(n+1-j) + g(n-j)
      enddo
      g(1) = 2.0 * g(1)
      g(n) = 2.0 * g(n)
      do j = 1, n
         g(j) = g(j) / 4.0
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 11. smooth.
c     smooth an array with a moving average.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine smooth ( g, n )
      implicit real*8 (a-h,o-z)
      dimension g(*), h(301)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      nm = n - 1
      do j = 1, nm
         h(j) = g(j)
         g(j) = g(j) + g(j+1)
      enddo
      h(n) = g(n)
      do j = 1, nm
         h(n+1-j) = h(n+1-j) + h(n-j)
      enddo
      h(1) = h(n)
      g(n) = g(1)
      do j = 1, n
         g(j) = ( g(j)+h(j) ) / 4.0
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 12. lagp.
c     some sort of interpolation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine lagp ( ax, af, m, nl, x, f, df, iop, iper )
      implicit real*8 (a-h,o-z)
      dimension ax(*), af(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      zero = 0.0
      one = 1.0
      dax = ax(m) - ax(1)
      in = 1
      do 20 i = 1, m
      if ( ax(i) .gt. x ) go to 25
   20 continue
   25 continue
      in = i - 1
      if ( (in .eq. 0) .or. (in .eq. m) ) go to 30
      if ( ax(in+1)-x .lt. x-ax(in) ) in = in + 1
   30 continue
      inmm = (nl-0.1) / 2.0
      inpp = (nl+0.1) / 2.0
      nll = in - inmm
      nlr = in + inpp
      if ( ( (nl/2)*2 .eq. nl) .and. ( ax(in) .gt. x ) ) then
        nll = nll - 1
        nlr = nlr - 1
      endif
      if ( (nll.ge.1) .and. (nlr.le.m) ) go to 35
      if ( iper .eq. 1 ) go to 35
      if ( nlr .le. m ) go to 33
      nlr = m
      nll = nlr - nl + 1
      go to 35
   33 continue
      nll = 1
      nlr = nl
   35 continue
      if ( iop .eq. 1 ) go to 120
      f = zero
      do 100 i0 = nll, nlr
        call shft ( i0,i, axi, ax, m, dax )
      alag = one
      do 50 j0 = nll, nlr
         call shft ( j0,j, axj, ax, m, dax )
      if ( i0 .eq. j0 ) go to 50
      alag = alag * ( x-axj ) / ( axi - axj )
   50 continue
      f = f + alag * af(i)
  100 continue
      if ( iop .eq. 0 ) return
  120 continue
      df = zero
      do 400 i0 = nll, nlr
         call shft ( i0, i, axi, ax, m, dax )
      slag = zero
      do 300 id0 = nll, nlr
         call shft ( id0, id, axid, ax, m, dax )
      if ( id0 .eq. i0 ) go to 300
      alag = one
      do 200 j0 = nll, nlr
         call shft ( j0, j, axj, ax, m, dax )
      if ( j0 .eq. i0 ) go to 200
      if ( j0 .eq. id0 ) go to 160
      alag = alag*( x-axj ) / ( axi-axj )
      go to 200
  160 continue
      alag = alag / ( axi-axid )
  200 continue
      slag = slag + alag
  300 continue
      df = df + slag * af(i)
  400 continue
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 13. shft.
c     some sort of shift operation, used by lagp.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine shft ( i0,i, axi, ax, m, dax )
      implicit real*8 (a-h,o-z)
      dimension ax(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if ( (i0 .ge. 1) .and. (i0 .le. m) ) then
         i = i0
         axi = ax(i)
      else
         if ( i0 .gt. m ) then
            i = mod(i0,m) + 1
            axi = ax(i) + dax
         endif
         if ( i0 .lt. 1 ) then
            i = m + i0 - 1
            axi = ax(i) - dax
         endif
      endif
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end    
c-----------------------------------------------------------------------
c     subprogram 14. lagpe4.
c     routine used by trans and transdx.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine lagpe4 ( f0,m, x,f, df, iop )
      implicit real*8 (a-h,o-z)
      dimension f0(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      h = 1.0 / m
      f0(m+1) = f0(1)
      f0(m+2) = f0(2)
      f0(m+3) = f0(3)
      m0 = x / h
      x0 = m0 * h
      p = (x-x0) / h
      pm1 = p - 1.0
      pm2 = p - 2.0
      pp1 = p + 1.0
      pp2 = p + 2.0
      am1 = - p*pm1*pm2 / 6.0
      a00 =  pm1*pp1*pm2 / 2.0
      ap1 = - p*pp1*pm2 / 2.0
      ap2 =  p*pp1*pm1 / 6.0
      if ( m0 .eq. 0 ) fm1 = f0(m)
      if ( m0 .ge. 1 ) fm1 = f0(m0)
      f00 = f0(m0+1)
      fp1 = f0(m0+2)
      fp2 = f0(m0+3)
      f = am1*fm1 + a00*f00 + ap1*fp1 + ap2*fp2
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 15. lag.
c     some sort of routine used by main vacuum computation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lag(ax,af,m,nl,x,f,df,iop)
      IMPLICIT NONE

      REAL(8), DIMENSION(*), INTENT(IN) :: ax,af
      INTEGER, INTENT(IN) :: m,nl,iop
      REAL(8), INTENT(OUT) :: f,df

      REAL(8) :: alag,slag,x
      INTEGER :: i,id,jn,jnmm,jnpp,j,nll,nlr
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      jn=1

      DO i=1,m
         IF(ax(i)  >=  x)EXIT
      ENDDO

      jn=MAX(i-1,1)
      IF(jn /= m .AND. ax(jn+1)-x  <  x-ax(jn))jn=jn+1

      jnmm=(nl-.1)/2
      jnpp=(nl+.1)/2
      nll=jn-jnmm
      nlr=jn+jnpp

      IF(((nl/2)*2  ==  nl)  .AND.  (ax(jn)  >  x))THEN
         nll=nll-1
         nlr=nlr+1
      ENDIF

      IF(nlr > m)THEN
         nlr=m
         nll=nlr-nl+1
      ELSE IF(nll < 1)THEN
         nll=1
         nlr=nl
      ENDIF

      IF(iop  /=  1)THEN
         f=0
         DO i=nll,nlr
            alag=1
            DO j=nll,nlr
               IF(i  ==  j)CYCLE
               alag=alag*(x-ax(j))/(ax(i)-ax(j))
            ENDDO
            f=f+alag*af(i)
         ENDDO
         IF(iop  ==  0)RETURN
      ENDIF

      df=0
      DO i=nll,nlr
         slag=0
         DO id=nll,nlr
            IF(id  ==  i)CYCLE
            alag=1
            DO j=nll,nlr
               IF(j == i)CYCLE
               IF(j  /=  id)THEN
                  alag=alag*(x-ax(j))/(ax(i)-ax(j))
               ELSE
                  alag=alag/(ax(i)-ax(id))
               ENDIF
            ENDDO
            slag=slag+alag
         ENDDO
         df=df+slag*af(i)
      ENDDO
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 16. eigen.
c     computes eigenvalues and eigenvectors of a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine eigen(a,r,n,mv)
      implicit real*8 (a-h,o-z)
      dimension a(*),r(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
    5 range=1.0e-12
      if(mv-1) 10,25,10
   10 iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0
      if(i-j) 20,15,20
   15 r(ij)=1.0
   20 continue
   25 anorm=0.0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
   30 ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
   35 continue
      if(anorm) 165,165,40
   40 anorm=1.414*sqrt(anorm)
      anrmx=anorm*range/float(n)
      ind=0
      thr=anorm
   45 thr=thr/float(n)
   50 l=1
   55 m=l+1
   60 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
   62 if( abs(a(lm))-thr) 130,65,65
   65 ind=1
      ll=l+lq
      mm=m+mq
      x=0.5*(a(ll)-a(mm))
   68 y=-a(lm)/ sqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
   70 y=-y
   75 continue
      yp = 1.0 - y*y
      if ( abs(yp) .lt. 1.0e-10 ) yp = 0.0
      sinx=y/ sqrt(2.0*(1.0+sqrt(yp+0.0)))
      sinx2=sinx*sinx
   78 cosx= sqrt(1.0-sinx2)
      cosx2=cosx*cosx
      sincs =sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
   80 if(i-m) 85,115,90
   85 im=i+mq
      go to 95
   90 im=m+iq
   95 if(i-l) 100,105,105
  100 il=i+lq
      go to 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  115 if(mv-1) 120,125,120
  120 ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
      r(ilr)=x
  125 continue
      x=2.0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
  130 if(m-n) 135,140,135
  135 m=m+1
      go to 60
  140 if(l-(n-1)) 145,150,145
  145 l=l+1
      go to 55
  150 if(ind-1) 160,155,160
  155 ind=0
      go to 50
  160 if(thr-anrmx) 165,165,45
  165 iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
  170 x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
  175 do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
  180 r(imr)=x
  185 continue
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 17. mult.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mult ( a, b, c, ndim, ndim2, m, l )
      implicit real*8 (a-h,o-z)
      dimension a(ndim,ndim), b(ndim,ndim2), c(ndim,ndim2)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do j1 = 1, m
         do j2 = 1, l
            sum = 0.0
            do i = 1, m
               sum = sum + a(j1,i) * b(i,j2)
            enddo
            c(j1,j2) = sum
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 18. matmul1.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matmul1 ( a, b, nda,ndb, n, c,ndc )
      implicit real*8 (a-h,o-z)
      dimension a(nda,nda), b(ndb,ndb), c(ndc,ndc)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      do i = 1, n
         do j = 1, n
            sum = 0.0
            do k = 1, n
               sum = sum + a(i,k)*b(k,j)
            enddo
            c(i,j) = sum
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 19. matmul3.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matmul3 ( a, b, nda,ndb, la,lb,lc, c,ndc )
      implicit real*8 (a-h,o-z)
      dimension a(nda,nda), b(ndb,ndb), c(ndc,ndc)
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      do i = 1, la
         do j = 1, lc
            sum = 0.0
            do k = 1, lb
               sum = sum + a(i,k)*b(k,j)
            enddo
            c(i,j) = sum
         enddo
      enddo
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 20. indef4.
c     computes indefinite and definite integral.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine indef4 ( f, fin, dx, n1, n2, defint, iend )
      implicit real*8 (a-h,o-z)
      dimension f(*), fin(*)
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      dendf = 0.0
      fin(n1) = 0.0
      ind = 1
      if ( iend .eq. 1 ) then
         fin(n1+1) = dx * ( 9.0*f(n1) + 19.0*f(n1+1) - 5.0*f(n1+2)
     $        + f(n1+3) ) / 24.0
         dendf = dx * ( 9.0*f(n2) + 19.0*f(n2-1) - 5.0*f(n2-2)
     $        + f(n2-3) ) / 24.0
         ind = 2
      endif
      do i = n1+ind, n2-ind+1
         ip1 = i+1
         im1 = i-1
         im2 = i-2
         fin(i) = fin(im1) + dx * ( -f(im2) + 13.0*f(im1) + 13.0*f(i)
     $        -f(ip1) ) / 24.0
      enddo
      fin(n2) = fin(n2-ind+1) + dendf
      defint = fin(n2)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
