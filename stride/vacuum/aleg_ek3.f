c     
c.....................................................
      SUBROUTINE aleg(x,nloc,pm,pn,pp, aleg0,aleg1 )
c.....................................................
c     
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

      PARAMETER ( pye=3.141592653589793, pii=2.0/pye, sqpi=SQRT(pye),
     $            sqtwo=SQRT(2.0), half=0.5 )

c...  Sum of ak_i = pi/2. Sum of ae_i = pi/2 - 1.0

!::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     This stuff for Gaussian Itegration:

c$$$      REAL, DIMENSION(8):: tgaus, wgaus 
      REAL, DIMENSION(32):: tg32, wg32, xg32
      REAL, DIMENSION(5):: xu, xl
c$$$      REAL, DIMENSION(10):: cfac, wksp

!.... Weights and abscissae for 32 points gaussian quadrature.

      wg32(1)  =  0.007018610009470096600 
      wg32(2)  =  0.016274394730905670605
      wg32(3)  =  0.025392065309262059456
      wg32(4)  =  0.034273862913021433103
      wg32(5)  =  0.042835898022226680657
      wg32(6)  =  0.050998059262376176196
      wg32(7)  =  0.058684093478535547145
      wg32(8)  =  0.065822222776361846838
      wg32(9)  =  0.072345794108848506225
      wg32(10) =  0.078193895787070306472
      wg32(11) =  0.083311924226946755222
      wg32(12) =  0.087652093004403811143
      wg32(13) =  0.091173878695763884713
      wg32(14) =  0.093844399080804565639
      wg32(15) =  0.095638720079274859419
      wg32(16) =  0.096540088514727800567
      
      DO i = 1, 16
         wg32(16+i) = wg32(17-i)
      END DO

      xg32(1:16) = (/ -0.997263861849481563545,
     $ 0.985611511545268335400, 
     $ 0.964762255587506430774,
     $ 0.934906075937739689171, 
     $ 0.896321155766052123965, 
     $ 0.849367613732569970134, 
     $ 0.794483795967942406963, 
     $ 0.732182118740289680387, 
     $ 0.663044266930215200975, 
     $ 0.587715757240762329041, 
     $ 0.506899908932229390024, 
     $ 0.421351276130635345364, 
     $ 0.331868602282127649780, 
     $ 0.239287362252137074545, 
     $ 0.144471961582796493485, 
     $ 0.048307665687738316235 /)

!     xg32(17:32) = - (/ xg32(16:1) /)

      DO i = 1, 16
         xg32(16+i) = - xg32(17-i)
      END DO

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c$$$      pii = 2.0 / pye
c$$$      sqpi = SQRT ( pye )
c$$$      sqtwo = SQRT(2.0)
c$$$c      init = init + 1
c$$$      half = 0.5

      gam = sqpi
      xxq = x*x
      ysq = xxq - 1.0
      y = SQRT( ysq )
      w = x+y

      rhohatsq = 1.0 / ( 2.0 * y*w )
      rhohat = SQRT (rhohatsq)

      zk1i = w
      zk1 = 1.0/w
      zk1sq = zk1**2
      zk1sqrt = SQRT(zk1)
      zk1sqrti = SQRT(zk1i)

      errbu = 1.0e-8
      CALL ek3 ( zk1sq, ierbu, errbu, 10, elipk, elipe, convbu, kcbu )

      pn = pii * zk1sqrt * elipk
      pnp = pii * zk1sqrti * elipe

      aleg0 = pn

      pp = (  pnp - x*pn ) / (2.0*y)

      aleg1 = pp

c... Use Gaussian Integration if ...

      IF ( nloc*rhohat >= 0.1 ) GO TO 100

c      GO TO 100
c 85    CONTINUE

      kloc=0
      ak = 0.0

      IF ( nloc == 0 )  GO TO 10

    5 kloc=kloc+1
      ak = FLOAT(kloc)
      ak02 = 0.5 - ak
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
      twopi = 2.0 * pye

c.. gamn is  Gamma[1/2-n]
c   gamp is  Gamma[1/2-(n+1)]

      gamn = sqpi
      gamp = - 2.0 * sqpi

      IF ( nloc /= 0 ) THEN

         gamn = sqpi /
     $        PRODUCT( (/ ( -(i-1)-0.5, i = 1, nloc ) /) )
         gamp = - gamn / (nloc+0.5)

c$$$         nwrt = 1
c$$$         IF ( nwrt == 1 ) THEN
c$$$         WRITE (6, '("nloc, gamn, gamp = ", i3, 2es12.4)' )
c$$$     $        nloc, gamn, gamp
c$$$         nwrt = nwrt + 1
c$$$         END IF
         
      END IF

      gint  = sqtwo * pcoef**nloc * gint / (nloc*sqpi*gamn)
      gintp = sqtwo * pcoef**(nloc+1) * gintp / ((nloc+1.0)*sqpi*gamp)
      pn = gint  ! P(n)
      pp = gintp  ! P(n+1)

 500  CONTINUE
c$$$
c$$$         nwrt = 1
c$$$         IF ( nwrt == 1 ) THEN
c$$$            WRITE (23, '("nloc, x, rhohat, pn, pp = ", i3, 4es12.4)' )
c$$$     $        nloc, x, rhohat, pn, pp
c$$$         nwrt = nwrt + 1
c$$$         END IF

      RETURN
      END
c     

!...................................................
       SUBROUTINE ek3(eta,ier,error,maxit,cel1,cel2,convg, kounter)
!..................................................

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

       PARAMETER (pi=3.1415926535897932385 , pi2 = pi/2.0)


       pp     = 1.0
       aa     = 1.0
       bb1    = 1.0
       bb2    = ABS(eta)

       ier    = 0
       IF(eta .LE. 0.0  .OR.  eta > 1.0) THEN
          IF(eta < 0.0) ier   = 1
          IF(eta == 0.0) ier   = 2
          IF(eta > 1.0) ier   = 3
          cel1  = 0.0
          cel2  = 0.0
          RETURN
       END IF


       qcval  = SQRT(ABS(eta))
       aval0  = aa
       bval1  = bb1
       bval2  = bb2
       pval0  = pp

       eval   = qcval
       emval  = 1.0


       IF(pval0 > 0.0) THEN
          pval  = SQRT(pval0)
          aval1 = aval0
          aval2 = aval0
          bval1 = bval1/pval
          bval2 = bval2/pval

       else
          fval  = qcval*qcval
          tval  = 1.0  - fval
          gval  = 1.0  - pval0
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

       IF(sval /= 0.0) snorm = sval*sval
       IF(sval == 0.0) snorm = 1.0
       convg  = (sval-qcval)*(sval-qcval) / snorm
       cnvlog = ALOG10(ABS(convg))
       logcnv = IFIX(cnvlog)

       IF (kounter > maxit) THEN
          IF(logcnv < 0) ier   = logcnv
          IF(logcnv >= 0) ier   = -1
          cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
          cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))
          RETURN
       END IF


       cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
       cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))

       RETURN
       END
