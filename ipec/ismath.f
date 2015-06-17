c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     ISMATH: simple math routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ismath_mod
c     1. iscdftf
c     2. iscdftb
c     3. issect
c     4. issurfint
c     5. ipidl_3dsurf
c-----------------------------------------------------------------------
c     subprogram 0. ismath_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ismath_mod
      USE local_mod
      USE bicube_mod
      USE fspline_mod
      USE ipglobal_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iscdftf.
c     compute 1D truncated-discrete Fourier forward transform.
c-----------------------------------------------------------------------
      SUBROUTINE iscdftf(m,ms,func,fs,funcm)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,ms,fs
      INTEGER, DIMENSION(ms), INTENT(IN) :: m
      COMPLEX(r8), DIMENSION(0:fs), INTENT(IN) :: func
      COMPLEX(r8), DIMENSION(ms), INTENT(OUT) :: funcm

      funcm=0
      DO i=1,ms
         DO j=0,fs-1
            funcm(i)=funcm(i)+(1.0/fs)*
     $           func(j)*EXP(-(twopi*ifac*m(i)*j)/fs)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     complete Fourier transfrom, but slow.
c-----------------------------------------------------------------------
c      TYPE(cspline_type) :: cspl
c      DO i=1,ms
c         CALL cspline_alloc(cspl,fs,1)
c         cspl%xs=theta
c         cspl%fs(:,1)=func(:)*EXP(-(twopi*ifac*m(i)*theta))
c         CALL cspline_fit(cspl,"periodic")
c         CALL cspline_int(cspl)
c         funcm(i) = cspl%fsi(fs,1)
c         CALL cspline_dealloc(cspl)
c      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iscdftf
c-----------------------------------------------------------------------
c     subprogram 2. iscdfb.
c     compute 1D truncated-discrete Fourier backward transform.
c-----------------------------------------------------------------------
      SUBROUTINE iscdftb(m,ms,func,fs,funcm)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,ms,fs
      INTEGER, DIMENSION(ms), INTENT(IN) :: m
      COMPLEX(r8), DIMENSION(ms), INTENT(IN) :: funcm
      COMPLEX(r8), DIMENSION(0:fs), INTENT(OUT) :: func

      func=0

      DO i=0,fs-1
         DO j=1,ms
            func(i)=func(i)+funcm(j)*
     $           EXP((twopi*ifac*m(j)*i)/fs)
         ENDDO
      ENDDO

      func(fs)=func(0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iscdftb
c-----------------------------------------------------------------------
c     subprogram 3. issect.
c     compute root of discrete function by the secant method.
c     use when only function is monotonic.
c-----------------------------------------------------------------------
      FUNCTION issect(gridnum,xvec,yvec,yval)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: gridnum
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: xvec
      REAL(r8), DIMENSION(0:gridnum), INTENT(IN) :: yvec
      REAL(r8), INTENT(IN) :: yval

      REAL(r8) :: xval1,xval2,xval3,yval1,yval2,yval3,issect
      TYPE(spline_type) :: yfun

      CALL spline_alloc(yfun,gridnum,1)
      yfun%xs=xvec
      yfun%fs(:,1)=yvec
      CALL spline_fit(yfun,"extrap")

      xval1=xvec(0)
      yval1=yvec(0)
      xval2=xvec(gridnum)
      yval2=yvec(gridnum)

      xval3=yval
      CALL spline_eval(yfun,xval3,0)
      yval3=yfun%f(1)

      DO
         IF (ABS(yval3-yval) .LT. 1e-6) EXIT
         IF (yval3 .GT. yval) THEN
            xval2=xval3
            yval2=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ELSE
            xval1=xval3
            yval1=yval3
            xval3=xval1+(yval-yval1)*(xval2-xval1)/(yval2-yval1)
         ENDIF
         CALL spline_eval(yfun,xval3,0)
         yval3=yfun%f(1)
      ENDDO
      issect=xval3
      CALL spline_dealloc(yfun)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION issect
c-----------------------------------------------------------------------
c     subprogram 4. issurfint.
c     surface integration by simple method.
c-----------------------------------------------------------------------
      FUNCTION issurfint(func,fs,psi,wegt,ave)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: fs,wegt,ave
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), DIMENSION(0:fs), INTENT(IN) :: func

      INTEGER :: itheta
      REAL(r8) :: issurfint,area
      REAL(r8), DIMENSION(0:fs) :: thetas,delpsi,jacs

      issurfint=0
      area=0
      thetas=(/(itheta,itheta=0,fs)/)/REAL(fs,r8)

      DO itheta=0,fs-1
         CALL bicube_eval(rzphi,psi,thetas(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(thetas(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         jacs(itheta)=jac
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
      ENDDO

      IF (wegt==0) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==1) THEN
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           r(itheta)*jacs(itheta)*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE 
         DO itheta=0,fs-1
            issurfint=issurfint+
     $           jacs(itheta)*delpsi(itheta)*func(itheta)/r(itheta)/fs
         ENDDO
      ENDIF

      IF (ave==1) THEN
         DO itheta=0,fs-1
            area=area+jacs(itheta)*delpsi(itheta)/fs
         ENDDO
         issurfint=issurfint/area
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfint
c-----------------------------------------------------------------------
c     subprogram 5. issurfave.
c     surface average by simple method.
c-----------------------------------------------------------------------
      FUNCTION issurfave(func,fs,psi)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: fs
      REAL(r8), INTENT(IN) :: psi
      REAL(r8), DIMENSION(0:fs), INTENT(IN) :: func

      INTEGER :: itheta
      REAL(r8) :: issurfave

      issurfave=0

      DO itheta=0,fs-1
         issurfave=issurfave+0.5*func(itheta)/fs
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfave
c-----------------------------------------------------------------------
c     subprogram 5. isbubble.
c     performs a bubble sort in decreasing order of value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE isbubble(key,index,mmin,mmax)

      REAL(r8), DIMENSION(:), INTENT(IN) :: key
      INTEGER, DIMENSION(:), INTENT(INOUT) :: index
      INTEGER :: mmin,mmax

      LOGICAL :: switch
      INTEGER :: i,temp
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      switch= .TRUE.
      DO while(switch)
         switch= .FALSE.
         DO i=mmin,mmax-1
            IF(key(index(i)) < key(index(i+1)))THEN
               temp=index(i)
               index(i)=index(i+1)
               index(i+1)=temp
               switch= .TRUE.
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE isbubble
c-----------------------------------------------------------------------
c     subprogram 6. ipidl_3dsurf.
c     convert 2d data into 3d data.
c-----------------------------------------------------------------------
C     SIMPLE CODE TO READ/ECHO GLASSER DCON DATA
C     code to read/echo & gen input to cugd.f
C     MODIFIED TO READ ANY EIGENMODE
C     MODIFIED  19 DEC 2000    JMB
C     MODIFIED MARCH 2001  ( ADDED OPTION #6 )
C     MODIFIED JULY 2001  ( MADE WORK FOR N>1 )
C     MODIFIED DEC 2001  ( NOW READS IN POL & TOR EIGENVECTORS )
C     MODIFIED SEPT 2003  ( NOW ALLOWS ROTATION OF DCON DATA )
C                             WHEN GEN INPUT TO CUGD3
C                             ALSO IMPLICIT NONE
C     MODIFIED OCT 2004  PRODUCE DATA FOR IDL PLOTS
C     MODIFIED DEC 2004, PRODUCE PLOTS OF BN AT SEL.PHI FOR POLAR LIKE PLTS
C     LAST ALTERED MARCH 2007, PRODUCE DATA FOR IDL PLOTS (UNWRAPPED SURFACE)
c-----------------------------------------------------------------------
      SUBROUTINE ipidl_3dsurf(FILEIN,NNN,MMTHETA,NP,NT,DIST,FILE1)

      IMPLICIT NONE
C
      INTEGER(KIND=8), PARAMETER :: MS=80,MT=2100,MP=100000
C
      CHARACTER*8  :: DUMC
      CHARACTER*72 :: FILE1, FILE2, FILEIN
      CHARACTER*72 :: DATALINE
C     
      INTEGER(KIND=8) :: IFIELD, I, I0, I1, I2, IS, J1, J2, K,
     >     IGOTO, NP, NT, IMAT, N3, N8, NF, NE, MMTHETA, MSOL, NNN,
     >     MTHETAP, IMODE, NUMEL, I4, IC, ID
C     
      INTEGER(KIND=8) :: IELE(4,MP)
C     
      REAL(KIND=8) :: PIE, DTR, QSURF, AN, DUM, D2, PHIR,
     >     PHIDEG, CC, SS, ANS, BIG, DPHI, DARC, ARC, THK, PPHI,
     >     ANSR, ANSZ, CCC, SSS, ANSRD, ANSZD, AMAG,
     >     ANSBNC, ANSBNS, DIST, DR, DZ, DS, DRN, DZN, DXN,
     >     EXPHID, EXPHIR, BNMAX, PHASED, PHASER, RAVG
C     
      REAL(KIND=8) :: RR(MT), ZZ(MT), ARCLEN(MT), BNCOS(MT),BNSIN(MT),
     >     ET(MS), R2(MT), Z2(MT), BNC2(MT), BNS2(MT),
     >     XYZ(3,MP), XYZP(3,MP), XYZN(3,MP), BN(MP),
     >     ARCPOS(MT), RS(MT), ZS(MT), BPCOS(MT), BPSIN(MT),
     >     BTCOS(MT), BTSIN(MT), BPC2(MT), BPS2(MT), BTC2(MT),
     >     BTS2(MT)
C     
      PIE = 2.*ASIN(1.)
      DTR = 180./PIE
      IFIELD = 0                !   DEFAULT SETTING ( NOT FOUND )

      OPEN ( 10,FILE=FILEIN,STATUS='OLD')

      AN = REAL ( NNN )
      MTHETAP = MMTHETA + 0

      READ ( 10,1002) ( RR(I),I=MMTHETA,1,-1),DUM
 1002 FORMAT(5E16.8)

      READ ( 10,1002) ( ZZ(I),I=MMTHETA,1,-1),DUM

      READ ( 10,1002) ( BNCOS(I),I=MMTHETA,1,-1),DUM

      READ ( 10,1002) ( BNSIN(I),I=MMTHETA,1,-1),DUM

      IMODE=1

C     
 11   CLOSE ( 10 )

      DO  I=1,MTHETAP
         IF ( I .EQ. 1 ) THEN
            ARCLEN(1) = 0.
         ELSE
            D2 = (RR(I)-RR(I-1))*(RR(I)-RR(I-1))+
     $           (ZZ(I)-ZZ(I-1))*(ZZ(I)-ZZ(I-1))
            ARCLEN(I) = ARCLEN(I-1) + SQRT(D2)
         ENDIF
      ENDDO

 5001 FORMAT(A32)
 6011 FORMAT(I6,F16.8)
 6012 FORMAT(1X,'ENTER NAME OF FILE WITH ARC-LENGTH DATA')
 6013 FORMAT(1X,'PROBLEM  ARCPOS(NP) .GE. ARCLEN(MTHETAP)',
     >     /1X,'ARCPOS(NP),ARCLEN(MTHETAP) =',2E14.6,
     >     /1X,'CODE STOPS  ASAP')
 6014 FORMAT(1X,'PROBLEM  NP = ',I6,5X,'STOP ASAP')
 6024 FORMAT(2F12.6 )
 1105 FORMAT(I6,'     numnp')
 1110 FORMAT(5I6,F12.6)

      PHASED=0.0
      PHASER = PHASED / DTR
C     
      IF ( NT*NP .GT. MP ) THEN
         WRITE ( 6,*) ' SIZE PROBLEMS WITH NT,NP & MP  ',NT,NP,MP
         STOP
      ENDIF

C
C     DO SPLINE FITS TO MAKE INTERPOLATION AUTOMATIC
C     
      BIG = 1.E+32              !  SETS KEY IN SPLINE ( FROM N.R. )
C     
      CALL SPLINE ( ARCLEN, RR,MTHETAP,BIG,BIG, R2)
C     
      CALL SPLINE ( ARCLEN, ZZ,MTHETAP,BIG,BIG, Z2)
C     
      CALL SPLINE ( ARCLEN, BNCOS,MTHETAP,BIG,BIG, BNC2)
C     
      CALL SPLINE ( ARCLEN, BNSIN,MTHETAP,BIG,BIG, BNS2)

      DPHI = 2.*PIE / REAL(NT)

C     FIND MAX B-NORMAL

      BNMAX = 0.
      DO  I = 1, MTHETAP
         IF ( ABS(BNCOS(I)) .GT. BNMAX ) BNMAX = ABS(BNCOS(I))
         IF ( ABS(BNSIN(I)) .GT. BNMAX ) BNMAX = ABS(BNSIN(I))
      ENDDO

      DARC = ARCLEN(MTHETAP) / REAL(NP)
      ARC = -DARC
      DO  I=1,NP
         ARC = ARC + DARC
         ARCPOS(I) = ARC
      ENDDO

      DO  I = 1,NP
         ARC = ARCPOS(I)
         CALL SPLINT ( ARCLEN,RR,R2,MTHETAP, ARC,ANSR)
         CALL SPLINT ( ARCLEN,ZZ,Z2,MTHETAP, ARC,ANSZ)
         CALL SPLINTD ( ARCLEN,RR,R2,MTHETAP, ARC,DR )
         CALL SPLINTD ( ARCLEN,ZZ,Z2,MTHETAP, ARC,DZ )
C     
C     CALC NORMAL
C     
         DS = SQRT ( DR*DR + DZ*DZ )
         DRN = DR/DS
         DZN = DZ/DS
C     
         DXN = -DZN
         DZN = DRN              !   CROSS PRODUCT  PHI IS Y-DIR  R&Z (X) Y
C     
         RS(I) = ANSR + DIST*DXN
         ZS(I) = ANSZ + DIST*DZN !  THIS IS A SURF. CONFORMAL TO THE PLASMA
C     
      ENDDO

      NNN = 0
      PPHI = -DPHI
      DO  IS = 1,NT
         PPHI = PPHI + DPHI
         CC = COS( PPHI)
         SS = SIN( PPHI)
         CCC = COS( AN*PPHI+PHASER)
         SSS = SIN( AN*PPHI+PHASER)
         DO  I = 1,NP
            ARC = ARCPOS(I)
            NNN = NNN + 1
            CALL SPLINT ( ARCLEN, RR,R2,MTHETAP, ARC, ANSR)
            CALL SPLINT ( ARCLEN, ZZ,Z2,MTHETAP, ARC, ANSZ)

            CALL SPLINT ( ARCLEN, BNCOS, BNC2, MTHETAP, ARC, ANSBNC )
            CALL SPLINT ( ARCLEN, BNSIN, BNS2, MTHETAP, ARC, ANSBNS )
            BN(NNN) = ANSBNC*CCC + ANSBNS*SSS

C     MOVE SURFACE IF DIST .NE. 0.

            XYZ(1,NNN) = CC* ( ANSR + (RS(I)-ANSR)*BN(NNN)/BNMAX)
            XYZ(2,NNN) = SS* ( ANSR + (RS(I)-ANSR)*BN(NNN)/BNMAX)
            XYZ(3,NNN) =     ( ANSZ + (ZS(I)-ANSZ)*BN(NNN)/BNMAX)
         ENDDO
      ENDDO

      NE = 0
      DO  IS = 1,NT
         DO  I = 1,NP-1
            NE = NE + 1
            IELE(1,NE) = (IS-1)*NP +I
            IELE(2,NE) = IELE(1,NE) + NP
            IELE(3,NE) = IELE(2,NE) + 1
            IELE(4,NE) = IELE(1,NE) + 1
            IF ( IS .EQ. NT ) THEN
               IELE(2,NE) = I
               IELE(3,NE) = I+1
            ENDIF
         ENDDO
         IF ( IS .NE. NT ) THEN
            NE = NE + 1
            IELE(1,NE) = (IS-1)*NP + NP
            IELE(2,NE) = IELE(1,NE) + NP
            IELE(3,NE) =  IS   *NP  + 1
            IELE(4,NE) = (IS-1)*NP  + 1
         ELSEIF ( IS .EQ. NT ) THEN
            NE = NE + 1
            IELE(1,NE) = NP * NT
            IELE(2,NE) = NP
            IELE(3,NE) = 1
            IELE(4,NE) = NP*NT - NP + 1
         ENDIF
      ENDDO

      OPEN ( 11,FILE=FILE1,STATUS='UNKNOWN')
      I0 = 0
      I1 = 1                    !    CONVIENENT CONSTANTS LATER ON
      I4 = 4
C     
C     
      WRITE ( 11,1105)  NNN
      DO  I=1,NNN
         WRITE ( 11,2106) XYZ(1,I),XYZ(2,I),XYZ(3,I),BN(I),I
      ENDDO
 2106 FORMAT(3F10.5,E16.8,I6)
C     
      NUMEL = NE
      N3 = 0
      N8 = 0
      WRITE ( 11,2107) NUMEL
 2107 FORMAT(6I6,'    numel nmat n2 n3 n4 n8 ')
C     
      DO  I=1,NE
         WRITE ( 11,1110) I,IELE(1,I),IELE(2,I),IELE(3,I),IELE(4,I)
      ENDDO
C     
      CLOSE ( 11 )

      RETURN
      END SUBROUTINE ipidl_3dsurf

      SUBROUTINE SPLINE ( x, y, n, yp1, ypn, y2)
C     
      IMPLICIT NONE
C     
      INTEGER(KIND=8), PARAMETER :: NMAX=1800
C     
      INTEGER(KIND=8) :: N, I, K
C     
      REAL(KIND=8)    :: YP1, YPN, X(N), Y(N), Y2(N), P, QN, SIG,
     >     UN, U(NMAX)
C     
C     CHECK ADEQUACY ON NMAX
C     
      IF ( NMAX .LT. N ) THEN
 6001    FORMAT(1X,'IN SPLINE  NMAX < N  STOP ASAP',
     >        /1X,'NMAX N = ',2I5 )
         STOP
      ENDIF
C     
      IF (yp1  .GT.  .99e30) THEN
         y2(1) = 0.
         u(1) = 0.
      else
         y2(1) = -0.5
         u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF
      do 11 i = 2,n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.
         y2(i) = (sig-1.)/p
         u(i) = (6.*((y(i+1)-y(i))/(x(i+
     *        1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *        u(i-1))/p
 11   continue
      IF (ypn .GT. .99e30) THEN
         qn = 0.
         un = 0.
      else
         qn = 0.5
         un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k = n-1,1,-1
         y2(k) = y2(k)*y2(k+1)+u(k)
 12   continue
      return
      END SUBROUTINE SPLINE 

      SUBROUTINE SPLINT ( xa, ya, y2a, n, x, y)
C     
      IMPLICIT NONE
C     
      INTEGER(KIND=8) :: N, K, KHI, KLO
C     
      REAL(KIND=8)    :: X, Y, XA(N), Y2A(N), YA(N), A, B, H
C     
      klo = 1
      khi = n
 1    IF (khi-klo .GT. 1) THEN
         k = (khi+klo)/2
         IF (xa(k) .GT. x) THEN
            khi = k
         else
            klo = k
         ENDIF
         GO TO 1
      ENDIF
      h = xa(khi)-xa(klo)
      IF ( H .EQ. 0. ) THEN
         WRITE ( 6,*) 'bad xa input in splint   STOP ASAP'
         STOP
      ENDIF
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *     2)/6.
      return
      END SUBROUTINE SPLINT

      SUBROUTINE  SPLINTD ( xa, ya, y2a, n, x, YD)
C     
      IMPLICIT NONE
C     
      INTEGER(KIND=8) :: N, K, KHI, KLO
C     
      REAL(KIND=8)    :: X, Y, XA(N), Y2A(N), YA(N), A, B, H, YD
C     
      klo = 1
      khi = n
 1    IF (khi-klo .GT. 1) THEN
         k = (khi+klo)/2
         IF (xa(k) .GT. x) THEN
            khi = k
         else
            klo = k
         ENDIF
         GO TO 1
      ENDIF
      h = xa(khi)-xa(klo)
      IF ( H .EQ. 0. ) THEN
         WRITE ( 6,*) 'bad xa input in splintd   STOP ASAP'
         STOP
      ENDIF
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h

      YD  =  YA(KHI)/H - YA(KLO)/H +
     >     H*((3.*B*B-1.)*Y2A(KHI) - (3.*A*A-1.)*Y2A(KLO))/6.
C     
      return
      END SUBROUTINE SPLINTD
c-----------------------------------------------------------------------

      END MODULE ismath_mod
