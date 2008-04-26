c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     ISMATH: basic math routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ismath_mod
c     1. iscdftf
c     2. iscdftb
c     3. issect
c     4. issurfint
c     5. is3dsurf
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
            funcm(i) = funcm(i)+(1.0/fs)*
     $           func(j)*EXP(-(twopi*ifac*m(i)*j)/fs)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     complete Fourier transfrom does not work!
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
            func(i) = func(i)+funcm(j)*
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
c     surface integration by simple method in hamada.
c     __________________________________________________________________
c     wegt     : 0: no weight
c                1: weighted by r, that is, torque
c     ave      : 0: no average
c                1: average
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
      REAL(r8), DIMENSION(0:fs) :: delpsi

      issurfint = 0
      area = 0

      DO itheta=0,fs-1
         CALL bicube_eval(rzphi,psi,theta(itheta),1)
         rfac=SQRT(rzphi%f(1))
         eta=twopi*(theta(itheta)+rzphi%f(2))
         r(itheta)=ro+rfac*COS(eta)
         z(itheta)=zo+rfac*SIN(eta)
         jac=rzphi%f(4)
         w(1,1)=(1+rzphi%fy(2))*twopi**2*rfac*r(itheta)/jac
         w(1,2)=-rzphi%fy(1)*pi*r(itheta)/(rfac*jac)
         delpsi(itheta)=SQRT(w(1,1)**2+w(1,2)**2)
      ENDDO

      IF (wegt==0) THEN
         DO itheta=0,fs-1
            issurfint = issurfint + 
     $           jac*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ELSE IF (wegt==1) THEN
         DO itheta=0,fs-1
            issurfint = issurfint + 
     $           r(itheta)*jac*delpsi(itheta)*func(itheta)/fs
         ENDDO
      ENDIF

      IF (ave==1) THEN
         DO itheta=0,fs-1
            area = area + jac*delpsi(itheta)/fs
         ENDDO
         issurfint = issurfint/area
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION issurfint
c-----------------------------------------------------------------------
c     subprogram 5. issurfave.
c     surface average by simple method in hamada.
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

      issurfave = 0

      DO itheta=0,fs-1
         issurfave = issurfave + 0.5*func(itheta)/fs
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
C     
C      WRITE ( 6,6000)
C 6000 FORMAT(1X,'ENTER NAME OF INPUT FILE')
C      READ ( 5,5000) FILEIN
C 5000 FORMAT(A32)

      OPEN ( 10,FILE=FILEIN,STATUS='OLD')
C     
C      READ ( 10,1999) DATALINE
C 1999 FORMAT(////A72)
C      WRITE ( 6,1999) DATALINE
C      IF ( DATALINE(1:6) .EQ. 'mtheta' ) THEN
C         READ ( DATALINE,1001) MMTHETA,MSOL,NNN,QSURF
C 1001    FORMAT(9X,I4,9x,I4,7x,I4,10X,E16.8)
C      ELSEIF ( DATALINE(1:7) .EQ. 'mthsurf' ) THEN
C         READ ( DATALINE,1011) MMTHETA,MSOL,NNN,QSURF
C 1011    FORMAT(10X,I4,9x,I4,7x,I4,10X,E16.8)
C      ELSE
C         WRITE ( 6,6010) DATALINE
C 6010    FORMAT(1X,A72,/1X,'CONFUSED    STOP ASAP')
C         STOP
C      ENDIF

C      WRITE ( 6,1001) MMTHETA,MSOL,NNN,QSURF

C      IF ( MMTHETA .GT. MT .OR. MSOL .GT. MS ) THEN
C         WRITE ( 6,6009) MMTHETA,MT, MSOL,MS
C 6009    FORMAT(1X,'**** SIZE PROBLEM  CHECK PARAMETERS',
C     >        //1X,'MMTHETA,MT = ',2I5,
C     >        /1X,'MSOL,MS   = ',2I5,
C     >        //1X,'CODE STOPS ASAP')
C         STOP
C      ENDIF

      AN = REAL ( NNN )
      MTHETAP = MMTHETA + 0
C     
C     TRY GLASSER LOGIC
C     
C      READ ( 10,1000) DUMC
C      READ ( 10,1000) DUMC
C      READ ( 10,1000) DUMC
C 1000 FORMAT(A8)

      READ ( 10,1002) ( RR(I),I=MMTHETA,1,-1),DUM
 1002 FORMAT(5E16.8)

C      write (6,*) ' RR(I) read in '
C     
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM

      READ ( 10,1002) ( ZZ(I),I=MMTHETA,1,-1),DUM

C      write ( 6,*) ' ZZ(I) read in '
C     
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1002) ( ET(I),I=1,MSOL) ! was MSOL
C      write ( 6,*) ' ET(I) read in '
C      write ( 6,1002) ( ET(I),I=1,MSOL)
C     
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM

      READ ( 10,1002) ( BNCOS(I),I=MMTHETA,1,-1),DUM
C     
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
      READ ( 10,1002) ( BNSIN(I),I=MMTHETA,1,-1),DUM
C     
C      WRITE ( 6,6015)  1,MSOL
C 6015 FORMAT(1X,'MODES FROM',I3,' TO ',I5,' ENTER A VALID #')
C      READ ( 5,*) IMODE

      IMODE=1
C      IF ( IMODE .LT. 1 .OR. IMODE .GT. MSOL ) THEN
C         WRITE ( 6,6016) IMODE
C 6016    FORMAT(1X,' MODE # = ',I5,'   IS INVALID   CODE STOPS ASAP')
C         STOP
C      ENDIF
C     
C      IF ( IMODE .NE. 1 ) THEN  !  IMODE OK  NEED TO READ IN
C         DO  K = 2,IMODE
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1002) ( BNCOS(I),I=MMTHETA,1,-1),DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1002) ( BNSIN(I),I=MMTHETA,1,-1),DUM
C         ENDDO
C      ENDIF
C     NOW LOOK FOR START OF POLOIDAL & TOROIDAL FIELD INFO

C      DO  I = 1,100000
C         READ ( 10,1000,END=999 ) DUMC
C         IF ( DUMC .EQ. 'poloidal' ) GO TO 10
C      ENDDO
C     FALL THROUGH A PROBLEM
C      WRITE ( 6,6030)
C 6030 FORMAT(1X,'PROBLEM  FAILED TO FIND POL FIELD DATA')
C      GO TO 11
C 999  WRITE ( 6,6031)
C 6031 FORMAT(1X,'PROBLEM  FOUND END OF FILE  NO FIELD DATA')
C      GO TO 11
C     
C     FIRST PIECE OF POL FIELD DATA FOUND
C     
C 10   IFIELD = 1                !   THIS MEANS DATA HAS BEEN FOUND
C      READ ( 10,1000) DUM
C      READ ( 10,1002) ( BPCOS(I),I=MMTHETA,1,-1),DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1000) DUM
C      READ ( 10,1002) ( BPSIN(I),I=MMTHETA,1,-1),DUM
C      IF ( IMODE .NE. 1 ) THEN
C         DO  K = 2,IMODE
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1002) ( BPCOS(I),I=MMTHETA,1,-1),DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1000) DUM
C            READ ( 10,1002) ( BPSIN(I),I=MMTHETA,1,-1),DUM
C         ENDDO
C      ENDIF
C     
C      DO I = 1,100000
C         READ  ( 10,1000) DUMC
C         IF ( DUMC .EQ. 'toroidal' ) GO TO 12
C      ENDDO
C      WRITE ( 6,6032)
C 6032 FORMAT(1X,'PROBLEM  FOUND POL DATA FAILED TO FIND TOR DATA')
C      STOP
C 12   READ ( 10,1000) DUMC
C      READ ( 10,1002) ( BTCOS(I),I=MMTHETA,1,-1),DUM
C      READ ( 10,1000) DUMC
C      READ ( 10,1000) DUMC
C      READ ( 10,1000) DUMC
C      READ ( 10,1002) ( BTSIN(I),I=MMTHETA,1,-1),DUM
C     
C      IF ( IMODE .NE. 1 ) THEN
C         DO  K = 2,IMODE
C            READ ( 10,1000) DUMC
C            READ ( 10,1000) DUMC
C            READ ( 10,1000) DUMC
C            READ ( 10,1002) ( BTCOS(I),I=MMTHETA,1,-1),DUM
C            READ ( 10,1000) DUMC
C            READ ( 10,1000) DUMC
C            READ ( 10,1000) DUMC
C            READ ( 10,1002) ( BTSIN(I),I=MMTHETA,1,-1),DUM
C         ENDDO
C      ENDIF
C     
 11   CLOSE ( 10 )
C     
C      WRITE ( 6,6001) IMODE,ET(IMODE)
C 6001 FORMAT(/1X,'ENERGY IN MODE# ',I5,' = ',E16.8)
C     
      DO  I=1,MTHETAP
         IF ( I .EQ. 1 ) THEN
            ARCLEN(1) = 0.
         ELSE
            D2 = (RR(I)-RR(I-1))*(RR(I)-RR(I-1))+
     $           (ZZ(I)-ZZ(I-1))*(ZZ(I)-ZZ(I-1))
            ARCLEN(I) = ARCLEN(I-1) + SQRT(D2)
         ENDIF
      ENDDO

C      WRITE ( 6,6033) IFIELD
C 6033 FORMAT(1X,'IFIELD = ',I5)
C     
C 1    WRITE ( 6,6002)
C 6002 FORMAT(1X,'STOP ?  ENTER 1/0')
C      READ ( 5,*) IGOTO
C      IF ( IGOTO .NE. 1 ) THEN
C         WRITE ( 6,6003)
C 6003    FORMAT(1X,'MAKE A CHOICE ( ENTER 1 NUMBER )',
C     >        /5X,'ENTER 1  - LIST ALL R & Z',
C     >        /5X,'ENTER 2  - LIST ALL BNCOS & BNSIN',
C     >        /5X,'ENTER 3  - LIST B-NORMAL AT GIVEN PHI',
C     >        /5X,'ENTER 4  - LIST ALL BPCOS & BPSIN',
C     >        /5X,'ENTER 5  - LIST B-POLOIDAL AT GIVEN PHI',
C     >        /5X,'ENTER 6  - LIST ALL BTCOS & BTSIN',
C     >        /5X,'ENTER 7  - LIST B-TOROIDAL AT GIVEN PHI',
C     >        /5X,'ENTER 8  - LIST ARC-LENGTH',
C     >        /5X,'ENTER 9  - GENERATE-INPUT-TO-CUGD & STOP',
C     >        /5X,'             NEED DELTA PHI  ',
C     >        /5X,'ENTER 10 - GENERATE-CONFORMAL WALL & STOP',
C     >        /5X,'ENTER 12 - Bn AS NEW R & Z AT A PHI',
C     >        /5X,'ENTER 21 - GENERATE-FILES FOR 3-D IDL CONTOUR PLOTS',
C     >        /5X,'ENTER 22 - GENERATE-FILES FOR UNWRAPPED',
C     >        ' IDL CONTOUR PLOTS')
C         READ ( 5,*) IGOTO
C      ENDIF

 5001 FORMAT(A32)
 6011 FORMAT(I5,F16.8)
 6012 FORMAT(1X,'ENTER NAME OF FILE WITH ARC-LENGTH DATA')
 6013 FORMAT(1X,'PROBLEM  ARCPOS(NP) .GE. ARCLEN(MTHETAP)',
     >     /1X,'ARCPOS(NP),ARCLEN(MTHETAP) =',2E14.6,
     >     /1X,'CODE STOPS  ASAP')
 6014 FORMAT(1X,'PROBLEM  NP = ',I5,5X,'STOP ASAP')
 6024 FORMAT(2F12.6 )
 1105 FORMAT(I5,'     numnp')
 1110 FORMAT(6I5,F12.6)

C
C      WRITE ( 6,6028)
C 6028 FORMAT(1X,'CODE TO GEN IDL INPUT FOR COUNTOUR PLOTS',
C     >     /1X,'ENTER NEW FILENAME FOR THIS DATA')
C      READ ( 5,5001) FILE1
      
C      WRITE ( 6,6029)
C 6029 FORMAT(1X,'ENTER NPOL NTOR  MAX SURFACE DISPLACEMENT',
C     >     ' & ROTATION[DEG]')
C      READ ( 5,*) NP,NT,DIST,PHASED
C      NP=72
C      NT=72
C      DIST=0.0

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
C     
C     CALL SPLINE ( ARCLEN, BPCOS,MTHETAP,BIG,BIG, BPC2)
C     
C     CALL SPLINE ( ARCLEN, BPSIN,MTHETAP,BIG,BIG, BPS2)
C     
C     CALL SPLINE ( ARCLEN, BTCOS,MTHETAP,BIG,BIG, BTC2)
C     
C     CALL SPLINE ( ARCLEN, BTSIN,MTHETAP,BIG,BIG, BTS2)
C     
C      WRITE ( 6,*) ' SPLINE FITS COMPLETE '
C     
      DPHI = 2.*PIE / REAL(NT)

C     FIND MAX B-NORMAL

      BNMAX = 0.
      DO  I = 1, MTHETAP
         IF ( ABS(BNCOS(I)) .GT. BNMAX ) BNMAX = ABS(BNCOS(I))
         IF ( ABS(BNSIN(I)) .GT. BNMAX ) BNMAX = ABS(BNSIN(I))
      ENDDO

C      WRITE ( 6,*) ' BNMAX = ', BNMAX
C     
C      IF ( NP .GT. 0 ) THEN

      DARC = ARCLEN(MTHETAP) / REAL(NP)
      ARC = -DARC
      DO  I=1,NP
         ARC = ARC + DARC
         ARCPOS(I) = ARC
C         WRITE ( 6,6011) I,ARCPOS(I)
      ENDDO

C      ELSEIF ( NP .LT. 0 )  THEN
C         NP = -NP               ! CHANGE NP TO A POS. NUMBER
C         WRITE ( 6,6012)
C         READ ( 5,5001) FILE2
C         OPEN ( 12,FILE=FILE2,STATUS='OLD')
C         DO  I=1,NP
C            READ ( 12,*) ARCPOS(I)
C            WRITE ( 6,6011) I,ARCPOS(I)
C         ENDDO
C         CLOSE ( 12 )
C     CHECK LAST ARCPOS() < ARCLEN(MTHETAP)
C         IF ( ARCPOS(NP) .GE. ARCLEN(MTHETAP) ) THEN
C            WRITE ( 6,6013) ARCPOS(NP),ARCLEN(MTHETAP)
C            STOP
C         ENDIF
C      ELSE                      !  IF
C         WRITE ( 6,6014) NP
C         STOP
C      ENDIF
C     
      DO  I = 1,NP
         ARC = ARCPOS(I)
         CALL SPLINT ( ARCLEN,RR,R2,MTHETAP, ARC,ANSR)
         CALL SPLINT ( ARCLEN,ZZ,Z2,MTHETAP, ARC,ANSZ)
C         WRITE ( 6,6024 ) ANSR,ANSZ
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

C     WRITE ( 6,6025)
C     DO  I = 1,NP
C     WRITE ( 6,6024) RS(I),ZS(I)
C     ENDDO
C     
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

C     XYZ(1,NNN) = ANSR * CC
C     XYZ(2,NNN) = ANSR * SS
C     XYZ(3,NNN) = ANSZ

            CALL SPLINT ( ARCLEN, BNCOS, BNC2, MTHETAP, ARC, ANSBNC )
            CALL SPLINT ( ARCLEN, BNSIN, BNS2, MTHETAP, ARC, ANSBNS )
            BN(NNN) = ANSBNC*CCC + ANSBNS*SSS

C     MOVE SURFACE IF DIST .NE. 0.

            XYZ(1,NNN) = CC* ( ANSR + (RS(I)-ANSR)*BN(NNN)/BNMAX)
            XYZ(2,NNN) = SS* ( ANSR + (RS(I)-ANSR)*BN(NNN)/BNMAX)
            XYZ(3,NNN) =     ( ANSZ + (ZS(I)-ANSZ)*BN(NNN)/BNMAX)
         ENDDO
      ENDDO

C      WRITE ( 6,*) ' NNN = ',NNN
C     
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

C     NOW DO TOROIDAL STRIP AT + PIE
C     DO  IS = 1,NT-1
C     NE = NE + 1
C     IELE(1,NE) = (IS-1)*NP + NP
C     IELE(2,NE) = IELE(1,NE) + NP
C     IELE(3,NE) =  IS   *NP  + 1
C     IELE(4,NE) = (IS-1)*NP  + 1
C     ENDDO     !  BELOW IS THE LAST ELEMENT
C     NE = NE + 1
C     IELE(1,NE) = NP * NT
C     IELE(2,NE) = NP
C     IELE(3,NE) = 1
C     IELE(4,NE) = NP*NT - NP + 1
C     
C      WRITE ( 6,*) ' NE = ',NE
C     
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
 2106 FORMAT(3F10.5,E16.8,I5)
C     
      NUMEL = NE
      N3 = 0
      N8 = 0
      WRITE ( 11,2107) NUMEL
 2107 FORMAT(6I5,'    numel nmat n2 n3 n4 n8 ')
C     
      DO  I=1,NE
         WRITE ( 11,1110) I,IELE(1,I),IELE(2,I),IELE(3,I),IELE(4,I)
      ENDDO
C     
      CLOSE ( 11 )
C     
C      STOP
C     
C     
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
C         WRITE ( 6,6001) NMAX,N
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
C     y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
C     *2)/6.
C     
      YD  =  YA(KHI)/H - YA(KLO)/H +
     >     H*((3.*B*B-1.)*Y2A(KHI) - (3.*A*A-1.)*Y2A(KLO))/6. !  JMB99
C     
      return
      END SUBROUTINE SPLINTD
c-----------------------------------------------------------------------

      END MODULE ismath_mod
