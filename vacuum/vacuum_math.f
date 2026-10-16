c-----------------------------------------------------------------------
c     file vacuum_math.f.
c     mathematical subroutines used by Morrell Chance's vacuum code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      1. spl1d1
c      2. spl1d2
c      4. search
c      5. searchx
c      6. green
c      7. aleg_old
c      8. aleg
c      9. ek3
c     10. trans
c     11. transdx
c     12. transdxx
c     13. smooth0
c     14. smooth
c     15. lagp
c     16. shft
c     17. lagpe4
c     18. lag
c     19. eigen
c     20. mult
c     21. matmul1
c     22. matmul3
c     23. indef4
c     24. atan2m
c-----------------------------------------------------------------------
c     subprogram 1. spl1d1.
c     spline fitting routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine spl1d1(n,x,f,w,iop,ij,a,b,c)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension iop(*),x(*),f(*),w(*),a(*),b(*),c(*)
      data zz,oz,tz,sz/0.0e0,1.0e0,3.0e0,6.0e0/
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      k=n-1
      a(2)=-(x(2)-x(1))/sz
      b(2)=(x(3)-x(1))/tz
      w(ij+1)=(f(2*ij+1)-f(ij+1))/(x(3)-x(2))-(f(ij+1)-f(1))
     1     /(x(2)-x(1))
      SELECT CASE (n - 3)
      CASE (: -1)  ! Less than 0
          GOTO 3
      CASE (0)     ! Equal to 0
          GOTO 4
      CASE (1 :)   ! Greater than 0
          GOTO 3
      END SELECT
    3 do 10 i=3,k
      m=(i-1)*ij+1
      j1=m+ij
      j2=m-ij
      con=(x(i+1)-x(i-1))/tz
      don=(x(i)-x(i-1))/sz
      b(i)=con-(don**2)/b(i-1)
      e=(f(j1)-f(m))/(x(i+1)-x(i))-(f(m)-f(j2))/
     1     (x(i)-x(i-1))
      w(m)=e-(don*w(j2))/b(i-1)
      a(i)=-(don*a(i-1))/b(i-1)
   10 continue
    4 k1=(n-2)*ij+1
      c(n-1)=-((x(n)-x(n-1))/sz)/b(n-1)
      w(k1)=w(k1)/b(n-1)
      a(n-1)=a(n-1)/b(n-1)
      k2=k-1
      SELECT CASE (n - 3)
      CASE (: -1)  ! Less than 0
          GOTO 7
      CASE (0)     ! Equal to 0
          GOTO 8
      CASE (1 :)   ! Greater than 0
          GOTO 7
      END SELECT
    7 do 20 i=2,k2
      j=n-i
      con=(x(j+1)-x(j))/sz
      a(j)=(a(j)-con*a(j+1))/b(j)
      c(j)=-(con*c(j+1))/b(j)
      k3=(j-1)*ij+1
      m=k3+ij
      w(k3)=(w(k3)-con*w(m))/b(j)
   20 continue
    8 k4=(n-1)*ij+1
      SELECT CASE (iop(1)-5)
      CASE (: -1)  ! Less than 0
         GOTO 201
      CASE (0)     ! Equal to 0
         GOTO 200
      CASE (1 :)   ! Greater than 0
         GOTO 201
      END SELECT
  201 c1=w(1)
      SELECT CASE (iop(2)-5)
      CASE (: -1)  ! Less than 0
         GOTO 203
      CASE (0)     ! Equal to 0
         GOTO 202
      CASE (1 :)   ! Greater than 0
         GOTO 203
      END SELECT
  203 c2=w(k4)
      go to 205
  200 SELECT CASE (n - 4)
      CASE (: -1)  ! Less than 0
         GOTO 300
      CASE (0)     ! Equal to 0
         GOTO 302
      CASE (1 :)   ! Greater than 0
         GOTO 302
      END SELECT
  302 a1=x(1)-x(2)
      a2=x(1)-x(3)
      a3=x(1)-x(4)
      a4=x(2)-x(3)
      a5=x(2)-x(4)
      a6=x(3)-x(4)
      w(1)=f(1)*(oz/a1+oz/a2+oz/a3)-a2*a3*f(ij+1)/(a1*a4*a5)+
     1     a1*a3*f(2*ij+1)/(a2*a4*a6 )-a1*a2*f(3*ij+1)/(a3*a5*a6)
      go to 201
  202 SELECT CASE (n-4)
      CASE (: -1)  ! Less than 0
         GOTO 300
      CASE (0)     ! Equal to 0
         GOTO 303
      CASE (1 :)   ! Greater than 0
         GOTO 303
      END SELECT
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
   70 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 80
      CASE (0)     ! Equal to 0
         GOTO 50
      CASE (1 :)   ! Greater than 0
         GOTO 80
      END SELECT
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
   62 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 71
      CASE (0)     ! Equal to 0
         GOTO 63
      CASE (1 :)   ! Greater than 0
         GOTO 71
      END SELECT
   63 a(1)=-oz
      c(1)=zz
      go to 500
   71 bob=zz
      go to 500
   64 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 73
      CASE (0)     ! Equal to 0
         GOTO 76
      CASE (1 :)   ! Greater than 0
         GOTO 73
      END SELECT
   76 a(1)=-oz
      c(1)=zz
      w(1)=zz
      go to 500
   73 SELECT CASE (i-2)
      CASE (: -1)  ! Less than 0
         GOTO 81
      CASE (0)     ! Equal to 0
         GOTO 81
      CASE (1 :)   ! Greater than 0
         GOTO 82
      END SELECT
   81 bob=-c1
      go to 500
   82 bob=zz
      go to 500
   66 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 83
      CASE (0)     ! Equal to 0
         GOTO 84
      CASE (1 :)   ! Greater than 0
         GOTO 83
      END SELECT
   84 a(1)=-(x(2)-x(1))/tz
      c(1)=zz
      w(1)=-c1+(f(ij+1)-f(1))/(x(2)-x(1))
      go to 500
   83 SELECT CASE (i-2)
      CASE (: -1)  ! Less than 0
         GOTO 85
      CASE (0)     ! Equal to 0
         GOTO 85
      CASE (1 :)   ! Greater than 0
         GOTO 86
      END SELECT
   85 bob=(x(2)-x(1))/sz
      go to 500
   86 bob=zz
      go to 500
   68 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 87
      CASE (0)     ! Equal to 0
         GOTO 88
      CASE (1 :)   ! Greater than 0
         GOTO 87
      END SELECT
   88 a(1)=-oz
      c(1)=oz
      w(1)=zz
      go to 500
   87 bob=zz
  500 ml=iop(2)
      go to (120,130,140,150,140),ml
  120 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 121
      CASE (0)     ! Equal to 0
         GOTO 122
      CASE (1 :)   ! Greater than 0
         GOTO 121
      END SELECT
  122 a(n)=zz
      c(n)=-oz
      go to 70
  121 bill=zz
      go to 70
  130 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 131
      CASE (0)     ! Equal to 0
         GOTO 132
      CASE (1 :)   ! Greater than 0
         GOTO 131
      END SELECT
  132 a(n)=zz
      c(n)=-oz
      w(k4)=zz
      go to 70
  131 SELECT CASE (i-k)
      CASE (: -1)  ! Less than 0
         GOTO 134
      CASE (0)     ! Equal to 0
         GOTO 133
      CASE (1 :)   ! Greater than 0
         GOTO 134
      END SELECT
  133 bill=-c2
      go to 70
  134 bill=zz
      go to 70
  140 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 141
      CASE (0)     ! Equal to 0
         GOTO 142
      CASE (1 :)   ! Greater than 0
         GOTO 141
      END SELECT
  142 a(n)=zz
      c(n)=(x(n-1)-x(n))/tz
      w(k4)=c2-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  141 SELECT CASE (i-k)
      CASE (: -1)  ! Less than 0
         GOTO 143
      CASE (0)     ! Equal to 0
         GOTO 144
      CASE (1 :)   ! Greater than 0
         GOTO 143
      END SELECT
  144 bill=(x(n)-x(n-1))/sz
      go to 70
  143 bill=zz
      go to 70
  150 SELECT CASE (i-1)
      CASE (: -1)  ! Less than 0
         GOTO 151
      CASE (0)     ! Equal to 0
         GOTO 152
      CASE (1 :)   ! Greater than 0
         GOTO 151
      END SELECT
  152 a(n)=zz
      c(n)=(x(n-1)+x(1)-x(n)-x(2))/tz
      w(k4)=(f(ij+1)-f(1))/(x(2)-x(1))-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  151 SELECT CASE (i-2)
      CASE (: -1)  ! Less than 0
         GOTO 153
      CASE (0)     ! Equal to 0
         GOTO 154
      CASE (1 :)   ! Greater than 0
         GOTO 153
      END SELECT
  154 bill=(x(2)-x(1))/sz
      go to 70
  153 SELECT CASE (i-k)
      CASE (: -1)  ! Less than 0
         GOTO 155
      CASE (0)     ! Equal to 0
         GOTO 156
      CASE (1 :)   ! Greater than 0
         GOTO 155
      END SELECT
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
      w(m)=w(m)+a(i)*w(1)+c(i)*w(k4)
  110 continue
      go to 305
  300 write(3,*) ' spl1d1: Results incorrect because n<4.'
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
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      data wz,sz/2.0e0,6.0e0/
      dimension x(*),f(*),w(*)
      dimension tab(3)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mflag = 0
      IF (y-x(1)<0) THEN
         GOTO 10
      ELSE IF (y-x(1)==0) THEN
         GOTO 10
      ELSE
         GOTO 20
      END IF
   10 i=1
      go to 30
   20 IF (y-x(n)<0) THEN
         GOTO 15
      ELSE IF (y-x(n)==0) THEN
         GOTO 40
      ELSE
         GOTO 40
      END IF
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
c     subprogram 4. search.
c     finds cubic spline interval.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine search(xbar,x,n,i,mflag)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension x(*)
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
      i=j
  5   continue
 6    k=i
      mflag = 1
      do 115  l=2,n
      a=x(l-1)
      b=x(l)
      IF (sign(1.0d0,a)-sign(1.0d0,b)<0) THEN
         GOTO 7
      ELSE IF (sign(1.0d0,a)-sign(1.0d0,b)==0) THEN
         GOTO 113
      ELSE
         GOTO 8
      END IF
113   IF (a-b<0) THEN
         GOTO 7
      ELSE IF (a-b==0) THEN
         GOTO 115
      ELSE
         GOTO 8
      END IF
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
  11   IF (ixbar-sign(1.0d0,a)<0) THEN
          GOTO 111
       ELSE IF (ixbar-sign(1.0d0,a)==0) THEN
          GOTO 1111
       ELSE
          GOTO 2111
       END IF
 1111 IF (xbar-a<0) THEN
         GOTO 111
      ELSE IF (xbar-a==0) THEN
         GOTO 14
      ELSE
         GOTO 2111
      END IF
 2111 b=x(i+1)
      IF (ixbar-sign(1.0d0,b)<0) THEN
         GOTO 2112
      ELSE IF (ixbar-sign(1.0d0,b)==0) THEN
         GOTO 2113
      ELSE
         GOTO 12
      END IF
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
   16 write(3,*) ' search: xbar is outside range of table.'
      mflag=2
      return
  20   IF (ixbar-sign(1.0d0,a)<0) THEN
          GOTO 2120
       ELSE IF (ixbar-sign(1.0d0,a)==0) THEN
          GOTO 2121
       ELSE
          GOTO 111
       END IF
 2121  IF (xbar-a<0) THEN
          GOTO 2120
       ELSE IF (xbar-a==0) THEN
          GOTO 14
       ELSE
          GOTO 111
       END IF
 2120 b=x(i+1)
       IF (ixbar-sign(1.0d0,b)<0) THEN
          GOTO 12
       ELSE IF (ixbar-sign(1.0d0,b)==0) THEN
          GOTO 2122
       ELSE
          GOTO 2112
       END IF
 2122  IF (xbar-b<0) THEN
          GOTO 12
       ELSE IF (xbar-b==0) THEN
          GOTO 12
       ELSE
          GOTO 2112
       END IF
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
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      integer xbar,x(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      if(.false.) go to 8888
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
      i=j
  5   continue
 6    k=i
      mflag = 1
      do 115  l=2,n
      ia=x(l-1)
      ib=x(l)
      SELECT CASE (isign(1,ia)-isign(1,ib))
      CASE (: -1)  ! Less than 0
         GOTO 7
      CASE (0)     ! Equal to 0
         GOTO 113
      CASE (1 :)   ! Greater than 0
         GOTO 8
      END SELECT
 113   SELECT CASE (ia-ib)
       CASE (: -1)  ! Less than 0
          GOTO 7
       CASE (0)     ! Equal to 0
          GOTO 115
       CASE (1 :)   ! Greater than 0
          GOTO 8
       END SELECT
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
  11   SELECT CASE (ixbar-isign(1,ia))
       CASE (: -1)  ! Less than 0
          GOTO 111
       CASE (0)     ! Equal to 0
          GOTO 1111
       CASE (1 :)   ! Greater than 0
          GOTO 2111
       END SELECT
 1111 SELECT CASE (xbar-ia)
      CASE (: -1)  ! Less than 0
         GOTO 111
      CASE (0)     ! Equal to 0
         GOTO 14
      CASE (1 :)   ! Greater than 0
         GOTO 2111
      END SELECT
 2111 ib=x(i+1)
      SELECT CASE (ixbar-isign(1,ib))
      CASE (: -1)  ! Less than 0
         GOTO 2112
      CASE (0)     ! Equal to 0
         GOTO 2113
      CASE (1 :)   ! Greater than 0
         GOTO 12
      END SELECT
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
   16 write(3,*) ' searchx: xbar is outside range of table.'
      mflag=2
      return
  20   SELECT CASE (ixbar-isign(1,ia))
       CASE (: -1)  ! Less than 0
          GOTO 2120
       CASE (0)     ! Equal to 0
          GOTO 2121
       CASE (1 :)   ! Greater than 0
          GOTO 111
       END SELECT
 2121  SELECT CASE (xbar-ia)
       CASE (: -1)  ! Less than 0
          GOTO 2120
       CASE (0)     ! Equal to 0
          GOTO 14
       CASE (1 :)   ! Greater than 0
          GOTO 111
       END SELECT
 2120 ib=x(i+1)
       SELECT CASE (ixbar-isign(1,ib))
       CASE (: -1)  ! Less than 0
          GOTO 12
       CASE (0)     ! Equal to 0
          GOTO 2122
       CASE (1 :)   ! Greater than 0
          GOTO 2112
       END SELECT
 2122  SELECT CASE (xbar-ib)
       CASE (: -1)  ! Less than 0
          GOTO 12
       CASE (0)     ! Equal to 0
          GOTO 12
       CASE (1 :)   ! Greater than 0
          GOTO 2112
       END SELECT
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
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
      if (use_legacy_greens_function) then
         call aleg_old ( s,nloc, pm,pn,pp, aleg0,aleg1 )
      else
         call aleg ( s,nloc, pm,pn,pp, aleg0,aleg1 )
      endif
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
c     subprogram 7. aleg_old.
c     computes Legendre functions, using [Chance Phys. Plasmas 1997].
c
c     subroutine to calculate half integral legendre functions.
c     uses upwards recurrence relations starting from elliptic
c     integrals evaluated using Bulirsch's algorithm
c     these expressions are very bad for large values of nloc.
c     zkisq is ths the 1 - k**2 in Elliptic integeral parlance.
c.    This is now replaced by the new aleg subroutine below.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE aleg_old(x,nloc,pm,pn,pp, aleg0,aleg1)
      USE local_mod, only: r8
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: nloc
      REAL(r8), INTENT(OUT) :: pm,pn,pp,aleg0,aleg1

      INTEGER :: kloc
      REAL(r8) :: ak02,elipe,elipk,gam,s,v,w,x1,x2,x3,x4,xxq,y,ysq

      REAL(r8), PARAMETER :: pi=3.1415926535897931_8,pii=2/pi,
     $     sqpi=1.7724538509055159_8
      REAL(r8), PARAMETER ::
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
c     subprogram 8. aleg.
c     Computes half-integral Legendre functions, using methods outlined
c     in Chance J. Comp. Phys 221 (2007) 330-348.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE aleg(x,nloc,pm,pn,pp, aleg0,aleg1 )
c     subroutine to calculate half integral legendre functions.
c     uses upwards recurrence relations starting from elliptic
c     integrals evaluated using Bulirsch's algorithm
c     these expressions are very bad for large values of nloc.
c     zkisq is ths the 1 - k**2 in Elliptic integeral parlance.

c     This modified from the old aleg subroutine to use the
c     Bulirsch algorithms for the Elliptic functions.
c     The new integral representation of the Legendre function is used
c     here for n*rhohat >= 0.1

c     Reference: JCP 221 (2007) 330-348
      USE local_mod, only: r8
      IMPLICIT NONE
      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: nloc
      REAL(r8), INTENT(OUT) :: pm,pn,pp,aleg0,aleg1
      REAL(r8), PARAMETER :: pye=3.141592653589793_r8, pii=2.0_r8/pye,
     $            sqpi=SQRT(pye), sqtwo=SQRT(2.0_r8), half=0.5_r8

c...  Sum of ak_i = pi/2. Sum of ae_i = pi/2 - 1.0

!::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     This stuff for Gaussian Itegration:

      REAL(r8), DIMENSION(32):: tg32, wg32, xg32
      REAL(r8), DIMENSION(5):: xu, xl

      REAL(r8) :: gam, xxq, ysq, y, w, rhohatsq, rhohat, zk1i, zk1,
     $ zk1sq, zk1sqrt, zk1sqrti, errbu, elipk, elipe, convbu,
     $ pnp, ak, ak02, gint, gintp, agaus, bgaus, ginti, gintip, tg0,
     $ tg02, tg1, tg1p, sinhtg1, sinhtg1p, sinhtg12, sinhtg12p, dnom,
     $ dnomp, anumr, pcoef, twopi, gamn, gamp
      INTEGER :: kcbu, kloc, ngauss, nng, ing, i, ig, ierbu

!.... Weights and abscissae for 32 points gaussian quadrature.

      wg32(1)  =  0.007018610009470096600_r8
      wg32(2)  =  0.016274394730905670605_r8
      wg32(3)  =  0.025392065309262059456_r8
      wg32(4)  =  0.034273862913021433103_r8
      wg32(5)  =  0.042835898022226680657_r8
      wg32(6)  =  0.050998059262376176196_r8
      wg32(7)  =  0.058684093478535547145_r8
      wg32(8)  =  0.065822222776361846838_r8
      wg32(9)  =  0.072345794108848506225_r8
      wg32(10) =  0.078193895787070306472_r8
      wg32(11) =  0.083311924226946755222_r8
      wg32(12) =  0.087652093004403811143_r8
      wg32(13) =  0.091173878695763884713_r8
      wg32(14) =  0.093844399080804565639_r8
      wg32(15) =  0.095638720079274859419_r8
      wg32(16) =  0.096540088514727800567_r8

      DO i = 1, 16
         wg32(16+i) = wg32(17-i)
      END DO

      xg32(1:16) = (/ -0.997263861849481563545_r8,
     $ 0.985611511545268335400_r8,
     $ 0.964762255587506430774_r8,
     $ 0.934906075937739689171_r8,
     $ 0.896321155766052123965_r8,
     $ 0.849367613732569970134_r8,
     $ 0.794483795967942406963_r8,
     $ 0.732182118740289680387_r8,
     $ 0.663044266930215200975_r8,
     $ 0.587715757240762329041_r8,
     $ 0.506899908932229390024_r8,
     $ 0.421351276130635345364_r8,
     $ 0.331868602282127649780_r8,
     $ 0.239287362252137074545_r8,
     $ 0.144471961582796493485_r8,
     $ 0.048307665687738316235_r8 /)

!     xg32(17:32) = - (/ xg32(16:1) /)

      DO i = 1, 16
         xg32(16+i) = - xg32(17-i)
      END DO

c-----------------------------------------------------------------------
c     Computations.
c-----------------------------------------------------------------------

      gam = sqpi
      xxq = x*x
      ysq = xxq - 1.0_r8
      y = SQRT( ysq )
      w = x+y

      rhohatsq = 1.0_r8 / ( 2.0_r8 * y*w )
      rhohat = SQRT (rhohatsq)

      zk1i = w
      zk1 = 1.0_r8/w        ! This is k1 = SQRT(1-k**2) = SQRT(m_1)
      zk1sq = zk1**2        ! This is m_1
      zk1sqrt = SQRT(zk1)   ! This is m_1^(1/4)
      zk1sqrti = SQRT(zk1i) ! This is m_1^(-1/4)

      errbu = 1.0e-8_r8
      CALL ek3 ( zk1sq, ierbu, errbu, 10, elipk, elipe, convbu, kcbu )

      pn = pii * zk1sqrt * elipk
      pnp = pii * zk1sqrti * elipe

      aleg0 = pn

      pp = (  pnp - x*pn ) / (2.0_r8*y)

      aleg1 = pp

c... Use Gaussian Integration if ...
      IF ( nloc*rhohat >= 0.1 ) GO TO 100

      kloc=0
      ak = 0.0_r8

      IF ( nloc == 0 )  GO TO 10

    5 kloc=kloc+1
      ak = FLOAT(kloc)
      ak02 = 0.5_r8 - ak
      pm = pn
      pn = pp
      pp = -2.0*ak*x*pn/y - ak02*ak02*pm
      gam = gam / ak02
      IF ( kloc /= nloc )  GO TO 5

 10   CONTINUE

      GO TO 500

 100  CONTINUE

c...  use Gauss integration of the new integral representation
c     if n*rhohat >= 0.1
c...  The integration is done in nng segments [xl(ing),xu(ing)].
c     Each stored in gint.

      ngauss = 32
      nng = 1
      xl(1) = 0.0
      xu(1) = 5.0

      gint = 0.0
      gintp = 0.0

      DO 165 ing = 1, nng

!.....xl, xu are the lower and upper limits of the gaussian integration
!     The integration is done in nng sections
!     This will calculate P(n) and P(n+1) together.
!        variables for P(n+1) will usually have p appended.

         agaus = half*( xu(ing)+xl(ing) )
         bgaus = half*( xu(ing)-xl(ing) )

         tg32(1:32) = agaus + xg32(1:32) * bgaus

         ginti = 0.0
         gintip = 0.0

         DO ig = 1, ngauss
            tg0 = tg32(ig)
            tg02 = tg0**2
            tg1  = tg02 / (2.0*nloc)
            tg1p = tg02 / (2.0*nloc+2.0)
            sinhtg1  = SINH(tg1)
            sinhtg1p = SINH(tg1p)
            sinhtg12  = sinhtg1  * sinhtg1
            sinhtg12p = sinhtg1p * sinhtg1p
            dnom  = x * sinhtg12  + sinhtg1 *SQRT(1.0 + sinhtg12)
            dnomp = x * sinhtg12p + sinhtg1p*SQRT(1.0 + sinhtg12p)
            dnom  = SQRT(dnom)
            dnomp = SQRT(dnomp)
            anumr = tg0 * EXP(-tg02)
            ginti  = ginti  + wg32(ig)*anumr / dnom
            gintip = gintip + wg32(ig)*anumr / dnomp
         END DO                 ! 32 point Gaussian

         ginti  = bgaus * ginti
         gintip = bgaus * gintip
         gint  = gint  + ginti
         gintp = gintp + gintip

 165  CONTINUE                  !  Gaussian integration segments

c... Now calculate the coeficients for the Legendre functions.

      pcoef = SQRT ( (x-1.0)/(x+1.0) )
      twopi = 2.0_r8 * pye

c.. gamn is  Gamma[1/2-n]
c   gamp is  Gamma[1/2-(n+1)]

      gamn = sqpi
      gamp = - 2.0_r8 * sqpi

      IF ( nloc /= 0 ) THEN

         gamn = sqpi /
     $        PRODUCT( (/ ( -(i-1)-0.5, i = 1, nloc ) /) )
         gamp = - gamn / (nloc+0.5)

      END IF

      gint  = sqtwo * pcoef**nloc * gint / (nloc*sqpi*gamn)
      gintp = sqtwo * pcoef**(nloc+1) * gintp / ((nloc+1.0)*sqpi*gamp)
      pn = gint  ! P(n)
      pp = gintp  ! P(n+1)

 500  CONTINUE

c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 9. ek3.
c     Computes complete elliptic integrals of first and second kind.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ek3(eta,ier,error,maxit,cel1,cel2,convg, kounter)

!  Compute the complete elliptic integral of first and second kind
!      cel(kc,p,a,b).
!  Bulirsch's method. Numerical Recipes, modified by Turnbull to
!    calculate both K and E simultaneously.

!  Returns cel1 = K, cel2 = E.
!  Precision is error**2,

!  eta, the complementary parameter, (1 - k^2), is the square of
!          the argument kc
!  p   is 1
!  a   is 1
!  b   is 1 for the first kind and b is eta( = kc**2) for the second kind

      USE local_mod, only: r8
      IMPLICIT NONE

      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8 ,
     $                       pi2 = pi/2.0_r8

      REAL(r8), INTENT(IN) :: eta, error
      INTEGER, INTENT(IN) :: maxit
      REAL(r8), INTENT(OUT) :: cel1, cel2, convg
      INTEGER, INTENT(OUT) :: ier, kounter

      REAL(r8) :: pp, aa, bb1, bb2, qcval, aval0, bval1, bval2, pval0,
     $  eval, emval, pval, aval1, aval2, fval, tval, gval, qval1, qval2,
     $  hval1, hval2, rval, sval, snorm, cnvlog

      INTEGER :: logcnv

      pp     = 1.0_r8
      aa     = 1.0_r8
      bb1    = 1.0_r8
      bb2    = ABS(eta)

      ier    = 0
      IF(eta .LE. 0.0_r8  .OR.  eta > 1.0_r8) THEN
         IF(eta < 0.0_r8) ier   = 1
         IF(eta == 0.0_r8) ier   = 2
         IF(eta > 1.0_r8) ier   = 3
         cel1  = 0.0_r8
         cel2  = 0.0_r8
         RETURN
      END IF

c-----------------------------------------------------------------------
c     Computations.
c-----------------------------------------------------------------------

      qcval  = SQRT(ABS(eta))
      aval0  = aa
      bval1  = bb1
      bval2  = bb2
      pval0  = pp

      eval   = qcval
      emval  = 1.0_r8


      IF(pval0 > 0.0_r8) THEN
         pval  = SQRT(pval0)
         aval1 = aval0
         aval2 = aval0
         bval1 = bval1/pval
         bval2 = bval2/pval

      else
         fval  = qcval*qcval
         tval  = 1.0_r8  - fval
         gval  = 1.0_r8  - pval0
         fval  = fval - pval0
         qval1 = tval*(bval1 - aval0*pval0)
         qval2 = tval*(bval2 - aval0*pval0)

         pval  = SQRT(fval/gval)
         aval1 = (aval0 - bval1) / gval
         aval2 = (aval0 - bval2) / gval
         bval1 =  aval1*pval - qval1/(gval*gval*pval)
         bval2 =  aval2*pval - qval2/(gval*gval*pval)
      END IF


      kounter = 0
100   CONTINUE
      kounter = kounter + 1

      hval1  = aval1
      hval2  = aval2
      aval1  = aval1 + bval1/pval
      aval2  = aval2 + bval2/pval
      rval   = eval/pval
      bval1  = bval1 + hval1*rval
      bval1  = bval1 + bval1
      bval2  = bval2 + hval2*rval
      bval2  = bval2 + bval2
      pval   = rval + pval

      sval   = emval
      emval  = qcval + emval

      IF (ABS(sval-qcval) > sval*error) THEN
         qcval  = SQRT(eval)
         qcval  = qcval + qcval
         eval   = qcval*emval
         GO TO 100
      END IF

      IF(sval /= 0.0_r8) snorm = sval*sval
      IF(sval == 0.0_r8) snorm = 1.0_r8
      convg  = (sval-qcval)*(sval-qcval) / snorm
      IF ( convg <= 1.0e-100_r8 ) convg = 1.0e-100_r8
      cnvlog = LOG10(ABS(convg))
      logcnv = INT(cnvlog)

      IF (kounter > maxit) THEN
         IF(logcnv < 0) ier   = logcnv
         IF(logcnv >= 0) ier   = -1
         cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
         cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))
         RETURN
      END IF


      cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
      cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))

c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 10. trans.
c     Interpolates input variables of length mthin to length mth.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine trans ( vecin,mthin, vecout,mth )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 11. transdx.
c     transforms (translates?) something.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine transdx ( vecin,mthin, vecout,mth, dx0 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 12. transdxx.
c     transforms (translates?) something.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine transdxx ( vecin,mthin, vecout,mth, dx0,dx1 )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension vecin(*), vecout(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)

      if ( (mth .ne. mthin) .or. (abs(dx0) .gt. 1.e-6) ) go to 20

      do 10 i = 1, mth + 2
         vecout(i) = vecin(i)
 10   continue

      return

 20   continue

      do 50 i = 1, mth
         ai = i-1 + dx1
         x = ai / mth + dx0 / mthin
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
 50   continue

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
c     subprogram 13. smooth0.
c     smooth an array with a moving average.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine smooth0 ( g, n )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 14. smooth.
c     smooth an array with a moving average.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine smooth ( g, n )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 15. lagp.
c     some sort of interpolation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine lagp ( ax, af, m, nl, x, f, df, iop, iper )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 16. shft.
c     some sort of shift operation, used by lagp.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine shft ( i0,i, axi, ax, m, dax )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 17. lagpe4.
c     routine used by trans and transdx.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine lagpe4 ( f0,m, x,f, df, iop )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 18. lag.
c     some sort of routine used by main vacuum computation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE lag(ax,af,m,nl,x,f,df,iop)
      USE local_mod, ONLY: r8
      IMPLICIT NONE

      REAL(r8), DIMENSION(*), INTENT(IN) :: ax,af
      INTEGER, INTENT(IN) :: m,nl,iop
      REAL(r8), INTENT(OUT) :: f,df

      REAL(r8) :: alag,slag,x
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
c     subprogram 19. eigen.
c     computes eigenvalues and eigenvectors of a matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine eigen(a,r,n,mv)
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
      dimension a(*),r(*)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
    5 range=1.0e-12
      if(mv-1 == 0) goto 25
   10 iq=-n
      do 21 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0
      if(i-j .ne. 0) goto 20
   15 r(ij)=1.0
   20 continue
   21 continue
   25 anorm=0.0
      do 36 i=1,n
      do 35 j=i,n
      if(i-j .NE. 0) then
         ia=i+(j*j-j)/2
      endif
      anorm=anorm+a(ia)*a(ia)
   35 continue
   36 continue
      if(anorm .LE. 0) then
         goto 165
      else
         goto 40
      endif
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
   62 if( abs(a(lm))-thr .GE. 0) then
         goto 65
      else
         goto 130
      endif
   65 ind=1
      ll=l+lq
      mm=m+mq
      x=0.5*(a(ll)-a(mm))
   68 y=-a(lm)/ sqrt(a(lm)*a(lm)+x*x)
      if(x .ge. 0) go to 75
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
      if(i-l == 0) goto 115
   80 SELECT CASE (i-m)
      CASE (: -1)  ! Less than 0
         GOTO 85
      CASE (0)     ! Equal to 0
         GOTO 115
      CASE (1 :)   ! Greater than 0
         GOTO 90
      END SELECT
   85 im=i+mq
      go to 95
   90 im=m+iq
   95 if(i-l .ge. 0) goto 105
  100 il=i+lq
      go to 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  115 if(mv-1 == 0) goto 125
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
  130 if(m-n == 0) goto 140
  135 m=m+1
      go to 60
  140 if(l-(n-1) == 0) goto 150
  145 l=l+1
      go to 55
  150 if(ind-1 .ne. 0) goto 160
  155 ind=0
      go to 50
  160 if(thr-anrmx .gt. 0) go to 45
  165 iq=-n
      do 186 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm) .ge. 0) go to 185
  170 x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1 == 0) goto 185
  175 do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
      r(imr)=x
  180 continue
  185 continue
  186 continue
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c     subprogram 20. mult.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine mult ( a, b, c, ndim, ndim2, m, l )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 21. matmul1.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matmul1 ( a, b, nda,ndb, n, c,ndc )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 22. matmul3.
c     matrix times matrix multiplication.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine matmul3 ( a, b, nda,ndb, la,lb,lc, c,ndc )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c     subprogram 23. indef4.
c     computes indefinite and definite integral.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      subroutine indef4 ( f, fin, dx, n1, n2, defint, iend )
      use local_mod, only: r8
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)
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
c-----------------------------------------------------------------------
c     subprogram 24. atan2m.
c     computes indefinite and definite integral.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      function atan2m ( z,x )
      USE vglobal_mod
      implicit real(r8) (a-h,o-z)
      implicit integer (i-n)

      epszer = 1e-30

      if ( abs(x) .le. epszer ) then
         atan2m = sign ( pye/2.0, z )
         return
      end if

      atan2m = atan2 (z,x)
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      return
      end


