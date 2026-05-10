      MODULE delta_mod

      USE sglobal_mod

      IMPLICIT NONE

      LOGICAL :: riccati_out,parflow_flag,PeOhmOnly_flag,verbose
      LOGICAL :: verbose_delta, IonScreening_flag, Pe_flag

      abstract interface
         function riccati_functions (inQ,inQ_e,inQ_i,inpr,inc_beta,inds,
     $     intau,inpe,inKp,iinQ,inx,iny) result(riccati_outcome)
            REAL(8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta
            REAL(8),INTENT(IN) :: inds,intau,inKp
            REAL(8),INTENT(IN),OPTIONAL :: iinQ,inx
            COMPLEX(8), INTENT(IN), OPTIONAL :: iny
            COMPLEX(8) :: riccati_outcome
         end function riccati_functions
      end interface

      procedure(riccati_functions), pointer :: riccati => null()

      CONTAINS

      subroutine select_riccati(IonScreening_flag, Pe_flag)
         LOGICAL :: IonScreening_flag, Pe_flag
         
         IF (IonScreening_flag) THEN
            IF (Pe_flag) THEN
               riccati => riccati_full
            ELSE
               riccati => riccati4
            ENDIF
         ELSE
            riccati => riccati3
         ENDIF
      end subroutine select_riccati

c-----------------------------------------------------------------------
c     calculate delta based on riccati w_der formulation.
c-----------------------------------------------------------------------
      FUNCTION riccati3(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,
     $     inKp,iinQ,inx,iny) result(riccati_outcome)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
      REAL(r8),INTENT(IN) :: intau,inKp
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ,inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati_outcome

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,iopt,mf

      REAL(r8) :: xintv,x,xout,rtol,jac,xmin
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,atol,rwork
      
      Q=inQ
      IF(present(iinQ)) Q=inQ+ifac*iinQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      pe=inpe
      c_beta=inc_beta
      ds=inds
      tau=intau
   
      IF ((layfac>0).AND.(ABS(Q-Q_e)<layfac)) THEN
         Q=Q_e+layfac*EXP(ifac*ATAN2(AIMAG(Q-Q_e),REAL(Q-Q_e)))
      ENDIF

      neq = 2
      itol = 2
      rtol = 1e-7            !1e-7*pr**0.4 ! !1e-7 at front 1e-6 !e-4
      ALLOCATE(atol(neq),y(1),dy(1))
      atol(:) = 1e-7*pr**0.4 ! 1e-8 !e-4
      itask = 2
      istate = 1
      iopt = 0
      mf = 10
      liw = 20    
      lrw = 22+16*neq
      ALLOCATE(iwork(liw),rwork(lrw))
      
!     MXSTEP? 
      iopt = 1
      iwork=0
      iwork(6)=10000 !5000 ! maximum step size, e.g. 50000
      rwork=0
!      x=10.0*(1.0+log10(Q/pr))
      x=20.0
      xmin=1e-3
      IF(present(inx)) x=inx
      xout=xmin
      y(1)=-c_beta/sqrt((1+tau))/ds*x**2.0 ! it was (1+tau*ds). To be updated.
      IF(present(iny)) y(1)=iny
!      y(1)=0.5-ifac*10.0
!      WRITE(*,*)y(1)
      

      IF (riccati_out) THEN
         istep = 1
         itask = 2
         OPEN(UNIT=bin_unit,FILE='slayer_riccati_profile_n'//
     $      TRIM(sn)//'.bin',STATUS='UNKNOWN',
     $      POSITION='REWIND',FORM='UNFORMATTED')
         
         OPEN(UNIT=out2_unit,FILE='slayer_riccati_profile_n'//
     $      TRIM(sn)//'.out',STATUS='UNKNOWN')
         WRITE(out2_unit,'(1x,3(a17))') "x","RE(y)","IM(y)"
         DO WHILE (x>xout)
            istep=istep+1
            CALL lsode(w_der,neq,y,x,xout,itol,rtol,atol,
     $           itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
            WRITE(bin_unit)REAL(x,4),REAL(REAL(y),4),REAL(AIMAG(y),4) 
            WRITE(out2_unit,'(1x,3(es17.8e3))') x,REAL(y),AIMAG(y)
         ENDDO        
         CLOSE(bin_unit)
         CLOSE(out2_unit)
      ELSE
         istep = 1
         itask = 1
         CALL lsode(w_der,neq,y,x,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)

      ENDIF

      ! w=0 when Q=Q_e. Why?
      
      CALL w_der(neq,x,y,dy)
      riccati_outcome=pi/dy(1)
      DEALLOCATE(atol,y,dy,iwork,rwork)      

      END FUNCTION riccati3
c-----------------------------------------------------------------------
c     calculate delta with the ion parallel flow (four-field model)
c     without electron viscosity and thermal conductivity.
c     Subroutines used in this function are w_derl, w_derr,
c     dw_der_wl, dw_der_wr, Transform_R, Update_Delta.
c     Reference : Y. Lee, J.-K. Park, & Y.-S. Na, (2024), Effect of
c         parallel flow on resonant layer responses in high beta plasmas.
c         Nucl. Fusion, 64(10), 106058.
c-----------------------------------------------------------------------
      FUNCTION riccati4(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,
     $     inKp,iinQ,inx,iny) result(riccati_outcome)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
      REAL(r8),INTENT(IN) :: intau,inKp
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ,inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati_outcome, Delta_old, Delta_new

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,lzw,iopt,mf
      INTEGER :: ml,mu,nrpd,ipar,ind,nR
      INTEGER :: nerr, jsv, meth, miter

      REAL(r8) :: xintv,x,xout,x_match,jac,xmin
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy,zwork,y1,dy1
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: pd, Rmatrix, dRmatrix

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,rwork,rtol,atol,RT
      
      Q=inQ
      IF(present(iinQ)) Q=inQ+ifac*iinQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      pe=inpe
      c_beta=inc_beta
      ds=inds
      tau=intau
      Kp=inKp
   
      IF ((layfac>0).AND.(ABS(Q-Q_e)<layfac)) THEN
         Q=Q_e+layfac*EXP(ifac*ATAN2(AIMAG(Q-Q_e),REAL(Q-Q_e)))
      ENDIF

      neq = 25
      itol = 4
      ALLOCATE(y(neq),dy(neq),pd(neq,neq),rtol(neq),atol(neq))
      rtol = 1e-12
      atol = 1e-12
      itask = 2
      istate = 1
      istep = 1
      itask = 2
      mf = 21
      liw = 60   
      lrw = 60
      lzw = 4400
      ALLOCATE(iwork(liw),rwork(lrw),zwork(lzw))
      
      iopt=0
      iwork=0
      iwork(1)=1
      iwork(2)=1
      iwork(5)=1e-9 !Initial H0.
      iwork(6)=100000 !5000 ! maximum step size, e.g. 50000
      iwork(17) = 8*neq + 2*neq**2
      iwork(18) = 29
      iwork(19) = 30+neq 

      rwork=0
      x=0.0
      xout=1e-9
      rwork(1)=1e-5

      y = 0.0
      dy = 0.0

      !--------------------------------------------------
      ! Forward sweeping
      !--------------------------------------------------
      ! leftside subregion
      DO WHILE (x<xout)
         istep=istep+1
         CALL ZVODE(w_derl,neq,y,x,xout,itol,rtol,atol,itask,
     $   istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_wl,mf,
     $   ipar)
      ENDDO        

      ! Transform R matrix from left to right region
      CALL Transform_R(y,x,neq)

      istate=1
      rtol=1e-9
      atol=1e-9
      Delta_old=0.0
      Delta_new=1.0

      ! rightside subregion
      DO WHILE (abs(Delta_new-Delta_old)/abs(Delta_new)>0.01)
         IF (riccati_out) THEN
            xout=xout+4.0
         ELSE
            xout=xout+1.0
         ENDIF
         Delta_old=Delta_new
         !WRITE(*,*) xout, Delta_new
         DO WHILE (x<xout)
            istep=istep+1
            ! Do Forward sweep by solving dR=A21+A22R-RA11-RA12R
            CALL ZVODE(w_derr,neq,y,x,xout,itol,rtol,atol,itask,
     $      istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_wr,mf,
     $      ipar)
            !WRITE(*,*) x, Delta_new
         ENDDO        
         CALL Update_Delta(Delta_new,y,x,neq)
         IF (verbose_delta) WRITE(*,*) xout, Delta_new
         !WRITE(*,*) xout, Delta_new
      END DO
      riccati_outcome=Delta_new

      ! Solution reconstruction routine if riccati_out is true
      IF (riccati_out) THEN
         ! Fix Riccati matrix size
         x_match=xout+4.0
         nR = int(x_match)*1000
         ALLOCATE(RT(nR),Rmatrix(nR,neq),y1(5),dy1(5))
         ind=1
         Rt=0
         Rmatrix=0
         RT(ind)=0
         Rmatrix(ind,:)=y(:)
   
         ! Redo Forward sweep to fill Riccati matrix
         istate=1
         x=0.0
         xout=1e-9
         y=0.0
         dy=0.0
         rtol = 1e-12
         atol = 1e-12
         rwork=0
         rwork(1)=1e-5

         ! leftside subregion
         DO WHILE (x<xout)
            istep=istep+1
            CALL ZVODE(w_derl,neq,y,x,xout,itol,rtol,atol,itask,
     $      istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_wl,mf,
     $      ipar)
         ENDDO        

         ! Transform R matrix from left to right region
         CALL Transform_R(y,x,neq)

         istate=1
         rtol=1e-9
         atol=1e-9

         DO WHILE (xout<x_match) 
           xout=xout+1e-3
           DO WHILE (x<xout)
              istep=istep+1
              CALL ZVODE(w_derr,neq,y,x,xout,itol,rtol,atol,itask,
     $        istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_wr,
     $        mf,ipar)
           ENDDO        
           IF (RT(ind) .lt. x) then
              ind=ind+1
              RT(ind) = x
              Rmatrix(ind,:) = y(:)
           endif
         ENDDO

         ! Allocate splR
         ! RT is empty from ind+1 to nR so using full RT raises an error
         ! since RT(i+1)-RT(i) = 0
         CALL cspline_alloc(splR,ind,25)
         splR%xs=RT(1:ind)
         splR%fs=Rmatrix(1:ind,:)
   
         ! Find spline fit
         CALL cspline_fit(splR,"extrap")
   
         !--------------------------------------------------
         ! Backward sweeping
         !--------------------------------------------------
         neq = 5
         itol = 4
         ! If the reconstruction takes too long time, relax tolerance
         ! If the solutions show oscillation, tighten tolerance
         rtol = 1e-8
         atol = 1e-8
!         DO i = 1,14
!           rtol(i) = 1.d-5
!           atol(i) = 1.d-5
!         END DO
         itask = 2
         iopt = 0
         lzw = 4400
         lrw = 60
         liw = 60
         iwork(1) = 1
         iwork(2) = 1
         iwork(6) = 10000
         nerr = 0
         
         IWORK(18) = 29
         JSV=1
         meth=2
         miter=1
         MF = JSV*(10*METH + MITER)
         
         IWORK(17) = 8*neq + 2*neq**2
         IWORK(19) = 30 + neq
         ISTATE = 1
         xout = 0.0
         RWORK(1) = xout*3
         
         ! Initialize y1 at X=XM
         CALL Init_yr1 (neq, x, y, y1)
         OPEN(UNIT=out4_unit,FILE='riccati_config_profile.out'
     $         ,STATUS='UNKNOWN')
         DO WHILE (x>xout)
            ! Do backward sweep by solving dy1 = (A11+A12R)y1
            CALL ZVODE(y_derr, NEQ, y1, x, xout, ITOL, RTOL, ATOL, 
     $         ITASK, ISTATE, IOPT, ZWORK, LZW, RWORK, LRW, IWORK, LIW,
     $         dy_der_yr, MF, IPAR)
            WRITE(out4_unit,'(es17.8e3)') x
         
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DREAL(y1(ind))
            END DO
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DIMAG(y1(ind))
            END DO
         
            ! Recover y2 by using y2 = Ry1
            CALL get_yr2 (NEQ, x, y1, dy1, IPAR)
         
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DREAL(dy1(ind))
            END DO
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DIMAG(dy1(ind))
            END DO
         END DO
         CLOSE(out4_unit)
      ENDIF
      END FUNCTION riccati4

      FUNCTION riccati_full(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,
     $     inKp,inpe,iinQ,inx,iny) result(riccati_outcome)
      USE global_mod

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
      REAL(r8),INTENT(IN) :: intau,inKp
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ,inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati_outcome, Delta_old, Delta_new

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,lzw,iopt,mf
      INTEGER :: ml,mu,nrpd,ipar,ind,nR
      INTEGER :: nerr, jsv, meth, miter

      REAL(r8) :: xintv,x,xout,x_match,jac,xmin
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy,zwork,y1,dy1
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: pd, Rmatrix, dRmatrix

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,rwork,rtol,atol,RT
      
      Q=inQ
      IF(present(iinQ)) Q=inQ+ifac*iinQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      pe=inpe
      c_beta=inc_beta
      ds=inds
      tau=intau
      Kp=inKp
   
      IF ((layfac>0).AND.(ABS(Q-Q_e)<layfac)) THEN
         Q=Q_e+layfac*EXP(ifac*ATAN2(AIMAG(Q-Q_e),REAL(Q-Q_e)))
      ENDIF

      neq = 49
      itol = 4
      ALLOCATE(y(neq),dy(neq),pd(neq,neq),rtol(neq),atol(neq)
     $        ,y1(7),dy1(7))
      rtol = 1e-15
      atol = 1e-15
      itask = 2
      istate = 1
      mf = 21
      liw = 100 !length of array IWORK
      lrw = 100 !length of RWORK
      lzw = 6000 !length of ZWORK
      ALLOCATE(iwork(liw),rwork(lrw),zwork(lzw))
      
      iopt=0
      iwork=0
      iwork(1)=1
      iwork(2)=1
      iwork(5)=1e-9 !H0 initial step size
      iwork(6)=100000 !5000 ! maximum step size, e.g. 50000
      iwork(17) = 8*neq + 2*neq**2
      iwork(18) = 29
      iwork(19) = 30+neq 

      rwork=0
      x=0.0
      xout=1e-4
      rwork(1)=1e-5

      y = 0.0
      dy = 0.0

      Delta_old=0.0
      Delta_new=1.0

      !--------------------------------------------------
      ! Forward sweeping
      !--------------------------------------------------
      DO WHILE (abs(Delta_new-Delta_old)/abs(Delta_new)>0.01)
         IF (riccati_out) THEN
            xout=xout+0.1
         ELSE
            xout=xout+10.0
         ENDIF
         Delta_old=Delta_new
         DO WHILE (x<xout)
            istep=istep+1
            ! Do Forward sweep by solving dR=A21+A22R-RA11-RA12R
            CALL ZVODE(w_der_full,neq,y,x,xout,itol,rtol,atol,itask,
     $         istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_w_full,
     $         mf,ipar)
         ENDDO        
         CALL Update_Delta_full(Delta_new,y,x,neq)
         write(*,*) xout, Delta_new
         IF (verbose_delta) WRITE(*,*) xout, Delta_new
      END DO
      riccati_outcome=Delta_new

      ! Solution reconstruction routine if riccati_out is true
      IF (riccati_out) THEN
         ! Fix Riccati matrix size
         x_match=xout+100.15
         nR = int(x_match)*10000
         ALLOCATE(RT(nR),Rmatrix(nR,neq))
         ind=1
         Rt=0
         Rmatrix=0
         RT(ind)=0
         Rmatrix(ind,:)=y(:)
   
         ! Redo Forward sweep to fill Riccati matrix
         istate=1
         x=0.0
         xout=0.0
         y=0.0
         dy=0.0
         DO WHILE (xout<x_match) 
           xout=xout+1e-5
           DO WHILE (x<xout)
              istep=istep+1
              CALL ZVODE(w_der_full,neq,y,x,xout,itol,rtol,atol,itask,
     $        istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dw_der_w_full,
     $        mf,ipar)
           ENDDO        
           IF (RT(ind) .lt. x) then
              ind=ind+1
              RT(ind) = x
              Rmatrix(ind,:) = y(:)
           endif
         ENDDO

         ! Allocate splR
         ! RT is empty from ind+1 to nR so using full RT raises an error
         ! since RT(i+1)-RT(i) = 0
         CALL cspline_alloc(splR,ind,49)
         splR%xs=RT(1:ind)
         splR%fs=Rmatrix(1:ind,:)
   
         ! Find spline fit
         CALL cspline_fit(splR,"extrap")
   
         !--------------------------------------------------
         ! Backward sweeping
         !--------------------------------------------------
         neq = 7
         itol = 4
         ! If the reconstruction takes too long time, relax tolerance
         ! If the solutions show oscillation, tighten tolerance
         rtol = 1e-13
         atol = 1e-13
!         DO i = 1,14
!           rtol(i) = 1.d-5
!           atol(i) = 1.d-5
!         END DO
         itask = 2
         iopt = 0
         lzw = 4400
         lrw = 60
         liw = 60
         iwork(1) = 1
         iwork(2) = 1
         iwork(6) = 10000
         nerr = 0
         
         IWORK(18) = 29
         JSV=1
         meth=2
         miter=1
         MF = JSV*(10*METH + MITER)
         
         IWORK(17) = 8*neq + 2*neq**2
         IWORK(19) = 30 + neq
         
         ISTATE = 1
         xout = 0.0
         RWORK(1) = xout*3
         
         ! Initialize y1 at X=XM
         CALL Init_y1_full (neq, x, y, y1)
         Delta_new=2.0/(y1(1)-x)
         riccati_outcome=Delta_new
         write(*,*) x, Delta_new
         RWORK(6) = 1e-5
         OPEN(UNIT=out4_unit,FILE='riccati_config_profile.out'
     $         ,STATUS='UNKNOWN')
         DO WHILE (x>xout)
            ! Do backward sweep by solving dy1 = (A11+A12R)y1
            CALL ZVODE(y_der_full, NEQ, y1, x, xout, ITOL, RTOL, ATOL, 
     $         ITASK, ISTATE, IOPT, ZWORK, LZW, RWORK, LRW, IWORK, LIW,
     $         dy_der_y_full, MF, IPAR)
            WRITE(out4_unit,'(es17.8e3)') x
         
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DREAL(y1(ind))
            END DO
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DIMAG(y1(ind))
            END DO
         
            ! Recover y2 by using y2 = Ry1
            CALL get_y2_full (NEQ, x, y1, dy1, IPAR)
         
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DREAL(dy1(ind))
            END DO
            DO ind = 1,7
               WRITE(out4_unit,'(es17.8e3)') DIMAG(dy1(ind))
            END DO
         END DO
         CLOSE(out4_unit)
      ENDIF
      END FUNCTION riccati_full

c-----------------------------------------------------------------------
c     riccati integration.
c-----------------------------------------------------------------------
      SUBROUTINE w_der(neq,x,y,dy)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8) :: G
      COMPLEX(r8) :: C1
      COMPLEX(r8) :: C1p
      COMPLEX(r8) :: C2
      COMPLEX(r8) :: C3
      COMPLEX(r8) :: C3p
      COMPLEX(r8) :: C2p
      COMPLEX(r8) :: A1
      COMPLEX(r8) :: A2
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      IF (parflow_flag) THEN
         C1=((1 + tau)*x**2*pe*
     $       (-(((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $              (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                   ifac*(Q - Q_e) + 
     $                   x**2*
     $                    (c_beta**2 + ifac*ds**2*(Q - Q_i)))/
     $                 (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $              (4*ds**2*pr*(1 + tau)*x**3 + 
     $                2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i))))
     $             /(ds**2*pr*(1 + tau)*x**4 + 
     $               ifac*(Q - Q_e) + 
     $               x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $           **2) + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $            (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                   2*x*(c_beta**2 + 
     $                      ifac*ds**2*(Q - Q_i)))/
     $              (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)) + 
     $              ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                 (ds**2*pr*(1 + tau)*x**4 + 
     $                   ifac*(Q - Q_e) + 
     $                   x**2*
     $                    (c_beta**2 + ifac*ds**2*(Q - Q_i))))
     $            /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))/
     $          (ds**2*pr*(1 + tau)*x**4 + ifac*(Q - Q_e) + 
     $            x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $         ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $            (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                 ifac*(Q - Q_e) + 
     $                 x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))
     $               )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $          (ds**2*pr*(1 + tau)*x**4 + ifac*(Q - Q_e) + 
     $            x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $     (ifac*Q + pr*x**2 + x**2*pe - 
     $       (ds**2*(1 + tau)*x**6*pe**2)/
     $        (c_beta**2*
     $          (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $            c_beta**2 + ifac*(Q - Q_e)))
     $        + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $        (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $           c_beta**2 + ifac*(Q - Q_e)))
	    
         C1p=((1 + tau)*x**2*pe*
     $        ((2*(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $             (4*ds**2*pr*(1 + tau)*x**3 + 
     $                2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               **2)/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $              ifac*(Q - Q_e) + 
     $              x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))**
     $            3 - ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $             (12*ds**2*pr*(1 + tau)*x**2 + 
     $               2*(c_beta**2 + ifac*ds**2*(Q - Q_i))))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $              ifac*(Q - Q_e) + 
     $              x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))**
     $          2 - (2*(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))+ 
     $               ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $            )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2)*
     $             (4*ds**2*pr*(1 + tau)*x**3 + 
     $               2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i))))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $              ifac*(Q - Q_e) + 
     $              x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))**
     $            2 - (2*(2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $             (4*ds**2*pr*(1 + tau)*x**3 + 
     $               2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i))))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $              ifac*(Q - Q_e) + 
     $              x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))**
     $            2 + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (-((12*ds**2*pr*(1 + tau)*x**2 + 
     $                    2*(c_beta**2 + ifac*ds**2*(Q - Q_i))
     $                 )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))
     $                + (2*(2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2
     $                - (2*(2*c_beta**2*x + 4*ds**2*pr*tau*x**3)**2*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $             )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**3
     $                + ((2*c_beta**2 + 12*ds**2*pr*tau*x**2)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $           )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))
     $            /(ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $          (2*(2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $             (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))+ 
     $               ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $            )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))
     $            /(ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $          ((2*c_beta**2 + 12*ds**2*pr*tau*x**2)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $        (ds**2*(1 + tau)*x**6*pe**2)/
     $         (c_beta**2*
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $          + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $         (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $            c_beta**2 + ifac*(Q - Q_e)))
     $      - ((1 + tau)*x**2*pe*
     $        (2*pr*x + 2*x*pe + 
     $          (ds**2*(1 + tau)*x**6*pe**2*
     $             (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                c_beta**2))/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                 c_beta**2 + 
     $                ifac*(Q - Q_e))**2) - 
     $          (6*ds**2*(1 + tau)*x**5*pe**2)/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                c_beta**2 + ifac*(Q - Q_e)
     $               )) - (ifac*(1 + tau)*x**2*pe*
     $             (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                c_beta**2)*Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $             **2 + (2*ifac*(1 + tau)*x*pe*
     $             Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $         *(-(((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $               (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $              /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $               (4*ds**2*pr*(1 + tau)*x**3 + 
     $                 2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               )/
     $             (ds**2*pr*(1 + tau)*x**4 + 
     $                ifac*(Q - Q_e) + 
     $                x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $          **2) + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $              /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)) + 
     $               ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $           )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))
     $            /(ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $          ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $         (ds**2*(1 + tau)*x**6*pe**2)/
     $          (c_beta**2*
     $            (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $            ) + (ifac*(1 + tau)*x**2*pe*
     $            Q_e)/
     $          (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $             c_beta**2 + ifac*(Q - Q_e)))
     $        **2 + (2*(1 + tau)*x*pe*
     $        (-(((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $               (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $                /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $               (4*ds**2*pr*(1 + tau)*x**3 + 
     $                 2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $               )/
     $             (ds**2*pr*(1 + tau)*x**4 + 
     $                ifac*(Q - Q_e) + 
     $                x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $         **2) + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $                /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))+ 
     $               ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $            )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))
     $            /(ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $          ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $        (ds**2*(1 + tau)*x**6*pe**2)/
     $         (c_beta**2*
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $          + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $         (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $            c_beta**2 + ifac*(Q - Q_e)))
	    
         C2=((1 + tau)*x**2*pe*
     $       ((ds**2*x**4*pe)/
     $          (c_beta**2*
     $            (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $            ) - (ifac*Q_e)/
     $          (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $             c_beta**2 + ifac*(Q - Q_e))
     $          + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $            (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                 ifac*(Q - Q_e) + 
     $                 x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))
     $               )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $          (ds**2*pr*(1 + tau)*x**4 + ifac*(Q - Q_e) + 
     $            x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $     (ifac*Q + pr*x**2 + x**2*pe - 
     $       (ds**2*(1 + tau)*x**6*pe**2)/
     $        (c_beta**2*
     $          (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $             c_beta**2 + ifac*(Q - Q_e)))
     $        + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $        (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $           c_beta**2 + ifac*(Q - Q_e)))
	    
         C2p=((1 + tau)*x**2*pe*
     $        (-((ds**2*x**4*pe*
     $               (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                  c_beta**2))/
     $             (c_beta**2*
     $               (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                   c_beta**2 + 
     $                  ifac*(Q - Q_e))**2)) + 
     $          (4*ds**2*x**3*pe)/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                c_beta**2 + ifac*(Q - Q_e)
     $               )) + (ifac*
     $             (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                c_beta**2)*Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $           **2 - ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4))*
     $             (4*ds**2*pr*(1 + tau)*x**3 + 
     $               2*x*(c_beta**2 + ifac*ds**2*(Q - Q_i))))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $              ifac*(Q - Q_e) + 
     $              x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))**
     $            2 + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (-((4*ds**2*pr*(1 + tau)*x**3 + 
     $                    2*x*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $              /(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)) + 
     $               ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $                  (ds**2*pr*(1 + tau)*x**4 + 
     $                    ifac*(Q - Q_e) + 
     $                    x**2*
     $                     (c_beta**2 + ifac*ds**2*(Q - Q_i)))
     $            )/(ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)**2))
     $            /(ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))) + 
     $          ((2*c_beta**2*x + 4*ds**2*pr*tau*x**3)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $        (ds**2*(1 + tau)*x**6*pe**2)/
     $         (c_beta**2*
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $          + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $         (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $            c_beta**2 + ifac*(Q - Q_e)))
     $      - ((1 + tau)*x**2*pe*
     $        (2*pr*x + 2*x*pe + 
     $          (ds**2*(1 + tau)*x**6*pe**2*
     $             (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                c_beta**2))/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                 c_beta**2 + 
     $                ifac*(Q - Q_e))**2) - 
     $          (6*ds**2*(1 + tau)*x**5*pe**2)/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                c_beta**2 + ifac*(Q - Q_e)
     $               )) - (ifac*(1 + tau)*x**2*pe*
     $             (2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $                c_beta**2)*Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $             **2 + (2*ifac*(1 + tau)*x*pe*
     $             Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $         *((ds**2*x**4*pe)/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                c_beta**2 + ifac*(Q - Q_e)
     $               )) - (ifac*Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e))
     $           + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $         (ds**2*(1 + tau)*x**6*pe**2)/
     $          (c_beta**2*
     $            (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $               c_beta**2 + ifac*(Q - Q_e))
     $            ) + (ifac*(1 + tau)*x**2*pe*
     $            Q_e)/
     $          (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $             c_beta**2 + ifac*(Q - Q_e)))
     $        **2 + (2*(1 + tau)*x*pe*
     $        ((ds**2*x**4*pe)/
     $           (c_beta**2*
     $             (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $                c_beta**2 + ifac*(Q - Q_e)
     $               )) - (ifac*Q_e)/
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e))
     $           + ((ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)*
     $             (1 - (ds**2*pr*(1 + tau)*x**4 + 
     $                  ifac*(Q - Q_e) + 
     $                  x**2*(c_beta**2 + 
     $                     ifac*ds**2*(Q - Q_i)))/
     $                (ifac*Q + c_beta**2*x**2 + ds**2*pr*tau*x**4)))/
     $           (ds**2*pr*(1 + tau)*x**4 + 
     $             ifac*(Q - Q_e) + 
     $             x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i)))))/
     $      (ifac*Q + pr*x**2 + x**2*pe - 
     $        (ds**2*(1 + tau)*x**6*pe**2)/
     $         (c_beta**2*
     $           (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $              c_beta**2 + ifac*(Q - Q_e)))
     $          + (ifac*(1 + tau)*x**2*pe*Q_e)/
     $         (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $            c_beta**2 + ifac*(Q - Q_e)))
      ELSE
         C1=0
         C1p=0
         C2=0
         C2p=0
      ENDIF		 
	
      IF (PeOhmOnly_flag) THEN
         G=((c_beta**2*pr*x**4 - Q*(Q - Q_i) + 
     $        ifac*(c_beta**2 + pr)*x**2*(Q - Q_i))/
     $      (ds**2*pr*(1 + tau)*x**4 + ifac*(Q - Q_e) + 
     $        x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))))*(x**2.0)
      ELSE
         G=(x**2*pe + 
     $     (c_beta**2*pr*x**4 - Q*(Q - Q_i) + 
     $        ifac*(c_beta**2 + pr)*x**2*(Q - Q_i))/
     $      (ds**2*pr*(1 + tau)*x**4 + ifac*(Q - Q_e) + 
     $        x**2*(c_beta**2 + ifac*ds**2*(Q - Q_i))))*(x**2.0)
      ENDIF
	 
      C3=x**2/(x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $     c_beta**2 + ifac*(Q - Q_e))
	 
      C3p=-((x**2*(2*x + (4*ds**2*(1 + tau)*x**3*pe)/
     $          c_beta**2))/
     $     (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $         c_beta**2 + (0,1)*(Q - Q_e))**2
     $     ) + (2*x)/
     $   (x**2 + (ds**2*(1 + tau)*x**4*pe)/
     $      c_beta**2 + (0,1)*(Q - Q_e))

      A1=(C1 + (C3p/C3)*(C2 + 1) + C2p)/(C2 + 1)
	  
      A2=(C1p + C1*(C3p/C3) - G/C3)/(C2 + 1)
	  
      dy(1)=(-A1 + 1/x)*y(1) - y(1)*y(1)/x - A2*x
	  
      RETURN
      END SUBROUTINE w_der

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in left region
c     dy (output) : Riccati differential matrix in left region
c     See Eqns. 42 and A.1-A.24
c-----------------------------------------------------------------------
      SUBROUTINE w_derl(neq,x,y,dy,ipar)

      INTEGER, INTENT(IN) :: neq,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), DIMENSION(5,5) :: ytemp
      COMPLEX(r8), DIMENSION(5,5) :: dytemp
      COMPLEX(r8), DIMENSION(10,10) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,1)=1.0
      A(2,8)=1.0
      A(3,1)=(Q-Q_e-tau*Q_e)*x/pr-tau*ds**4/c_beta**4
     $      *(Q-Q_e)*x**3-ifac*tau*ds**2/c_beta**2*(Q-Q_e)**2
     $      *x-ifac*tau*Q*ds**2/c_beta**4*(Q-Q_e)*x
      A(3,2)=4.0*tau*ds**2/c_beta**2*x
      A(3,4)=-4.0*tau*ds**2/c_beta**2*x
      A(3,5)=-tau*Q/c_beta**2*x-tau*Q/pr*x+ifac*tau*ds**2
     $      /c_beta**2*x**3
      A(3,6)=-2.0*tau*ds**2/c_beta**2*(Q-Q_e)
      A(3,7)=-tau*Q*Q_e/c_beta**4+2.0*tau*ds**2/c_beta**2
     $      -x**2/pr+tau*ds**4/c_beta**4*x**4+ifac*tau*ds**2
     $      /c_beta**2*(Q-Q_e)*x**2+ifac*tau*Q*ds**2/c_beta**4
     $      *x**2+ifac*tau*Q_e*ds**2/c_beta**4*x**2
      A(3,8)=tau*ds**2/c_beta**2*x**2+ifac*(Q-Q_i)/pr+ifac
     $      *tau*Q_e/c_beta**2
      A(3,9)=tau*Q**2/c_beta**4-2.0*tau*ds**2/c_beta**2+(1.0
     $      +tau)/pr*x**2-tau*ds**4/c_beta**4*x**4-ifac*tau
     $      *ds**2/c_beta**2*(Q-Q_e)*x**2-ifac*2.0*tau*Q*ds**2
     $      /c_beta**4*x**2
      A(3,10)=ifac*2.0*tau
      A(4,1)=ds**2/c_beta**2*(Q-Q_e)*x
      A(4,5)=-ifac*x
      A(4,7)=-ds**2/c_beta**2*x**2-ifac*Q_e/c_beta**2
      A(4,9)=ds**2/c_beta**2*x**2+ifac*Q/c_beta**2
      A(5,10)=1.0
      A(6,1)=ifac*(Q-Q_e)
      A(6,7)=-ifac*x
      A(6,9)=ifac*x
      A(7,2)=1.0
      A(8,3)=1.0
      A(9,4)=1.0
      A(10,1)=ifac*Q_e/pr
      A(10,5)=ifac*Q/pr
      A(10,9)=-ifac*x/pr
      
      DO i=1,5
        DO j=1,5
          ytemp(i,j)=y(i+5*(j-1)) 
        END DO
      END DO
      
      dytemp=A(6:,:5)+matmul(A(6:,6:),ytemp)-matmul(ytemp,A(:5,:5))
     $      -matmul(matmul(ytemp,A(:5,6:)),ytemp)

      DO i=1,5
        DO j=1,5
          dy(i+5*(j-1))=dytemp(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE w_derl

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in right region
c     dy (output) : Riccati differential matrix in right region
c     See Eqns. 42 and A.55-A.78
c-----------------------------------------------------------------------
      SUBROUTINE w_derr(neq,x,y,dy,ipar)

      INTEGER, INTENT(IN) :: neq,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), DIMENSION(5,5) :: ytemp
      COMPLEX(r8), DIMENSION(5,5) :: dytemp
      COMPLEX(r8), DIMENSION(10,10) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,6)=1.0
      A(2,7)=1.0
      A(3,2)=(Q_e*tau+Q)/c_beta**2*ds**2*x
      A(3,5)=-ifac*(Q_e*tau+Q)*x
      A(3,8)=ifac*(Q_e*tau+Q)/c_beta**2
      A(3,9)=-Q_e
      A(4,2)=x/pr-ifac*(Q-Q_i)/pr*ds**2/c_beta**2*tau*x
      A(4,5)=-1.0/pr*(Q-Q_i)*tau*x
      A(4,8)=1.0/pr*(Q-Q_i)*tau/c_beta**2
      A(4,9)=ifac*(Q-Q_i)/pr
      A(5,10)=1.0
      A(6,2)=ifac
      A(7,1)=2.0*(Q-Q_e)/x**2
      A(7,2)=-2.0/x**2+(tau+1.0)*ds**2/c_beta**2*x**2+ifac*(Q-Q_e)
      A(7,5)=-ifac*(tau+1.0)*x**2
      A(7,6)=-2.0*(Q-Q_e)/x
      A(7,7)=2.0/x
      A(7,8)=ifac*(tau+1.0)*x/c_beta**2
      A(7,9)=-x
      A(8,3)=1.0
      A(9,4)=1.0
      A(10,2)=ifac/pr/(Q-Q_e)*Q_e
      A(10,5)=ifac*Q/pr
      A(10,8)=-ifac/pr/(Q-Q_e)*x
      
      DO i=1,5
        DO j=1,5
          ytemp(i,j)=y(i+5*(j-1)) 
        END DO
      END DO
     
      dytemp=A(6:,:5)+matmul(A(6:,6:),ytemp)-matmul(ytemp,A(:5,:5))
     $      -matmul(matmul(ytemp,A(:5,6:)),ytemp)

      DO i=1,5
        DO j=1,5
          dy(i+5*(j-1))=dytemp(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE w_derr

      subroutine y_derr (neq,x,y,dy,ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ipar
      INTEGER :: i,j
      REAL(r8), INTENT(IN) :: x
      REAL(r8) :: Rtemp, dRtemp
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), DIMENSION(neq**2) :: R1D
      COMPLEX(r8), DIMENSION(neq,neq) :: R, U11
      COMPLEX(r8), DIMENSION(neq,neq*2) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,6) = 1.0
      A(2,7) = 1.0
      A(3,2) = (Q_e*tau+Q)/c_beta**2*ds**2*x
      A(3,5) = -ifac*(Q_e*tau+Q)*x
      A(3,8) = ifac*(Q_e*tau+Q)/c_beta**2
      A(3,9) = -Q_e
      A(4,2) = x/pr-ifac*(Q-Q_i)/pr*ds**2/c_beta**2*tau*x
      A(4,5) = -1.0/pr*(Q-Q_i)*tau*x
      A(4,8) = 1.0/pr*(Q-Q_i)*tau/c_beta**2
      A(4,9) = ifac*(Q-Q_i)/pr
      A(5,10) = 1.0
      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,5
        do j=1,5
          R(i,j) = R1D(i+5*(j-1))
        end do
      end do
      
      U11 = A(:,1:5) + matmul(A(:,6:10), R)
      dy = matmul(U11, y)
      return 
      END

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on the full model.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix
c     dy (output) : Riccati differential matrix
c     Documentation is not prepared.
c-----------------------------------------------------------------------
      SUBROUTINE w_der_full(neq,x,y,dy,ipar)
      INTEGER, INTENT(IN) :: neq,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), DIMENSION(7,7) :: ytemp
      COMPLEX(r8), DIMENSION(7,7) :: dytemp
      COMPLEX(r8), DIMENSION(14,14) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,8)  = 1.0
      A(2,9)  = c_beta**2/ds**2/Pe
      A(3,12) = tau/(1.0+tau)
      A(3,13) = 1.0/(1.0+tau)
      A(4,12) = -1.0/(1.0+tau)
      A(4,13) = 1.0/(1.0+tau)
      A(5,2)  = -ifac*x
      A(5,7)  = -ifac*x*c_beta**2/ds**2
      A(5,10) = -ifac*Q_e/ds**2
      A(5,11) = ifac*Q/ds**2
      A(5,12) = 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,12) = 1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(5,13) = -1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,13) = -1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(6,7)  = ifac*c_beta**2/ds**2*x
      A(6,10) = ifac*Q_e/ds**2
      A(6,11) = -ifac*Q/ds**2
      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
     $        - 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
!     $          - 1.0/(1.0+tau)*c_beta**2/ds**2
      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
     $        + 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
!     $        + 1.0/(1.0+tau)*c_beta**2/ds**2
      A(7,14) = 1.0/Pr
      A(8,2)  = 1.0
      A(9,1)  = -ifac/(1.0+tau)*((Q-Q_e)+(Q-Q_i)*Pe/Pr)
      A(9,2)  = (1.0+Pe/Pr)/(1.0+tau)
      A(9,7)  = -ifac*Q*Pe/Pr
      A(9,10) = ifac*(1.0+Pe/Pr)/(1.0+tau)*x
      A(9,11) = ifac*(tau*Pe/Pr-1.0)/(1.0+tau)*x
      A(10,3) = 1.0
      A(11,4) = 1.0
      A(12,5) = 1.0/Pe
      A(13,6) = 1.0/Pr
      A(14,1) = ifac*(Q-Q_i)/(1.0+tau)
      A(14,2) = -1.0/(1.0+tau)
      A(14,7) = ifac*Q
      A(14,10)= -ifac*x/(1.0+tau)
      A(14,11)= -ifac*tau*x/(1.0+tau)
      
      DO i=1,7
        DO j=1,7
          ytemp(i,j)=y(i+7*(j-1)) 
        END DO
      END DO
      
      dytemp=A(8:,:7)+matmul(A(8:,8:),ytemp)-matmul(ytemp,A(:7,:7))
     $      -matmul(matmul(ytemp,A(:7,8:)),ytemp)

      DO i=1,7
        DO j=1,7
          dy(i+7*(j-1))=dytemp(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE w_der_full

      subroutine y_der_full (neq,x,y,dy,ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ipar
      INTEGER :: i,j
      REAL(r8), INTENT(IN) :: x
      REAL(r8) :: Rtemp, dRtemp
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), DIMENSION(neq**2) :: R1D
      COMPLEX(r8), DIMENSION(neq,neq) :: R, U11
      COMPLEX(r8), DIMENSION(neq,neq*2) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,8)  = 1.0
      A(2,9)  = c_beta**2/ds**2/Pe
      A(3,12) = tau/(1.0+tau)
      A(3,13) = 1.0/(1.0+tau)
      A(4,12) = -1.0/(1.0+tau)
      A(4,13) = 1.0/(1.0+tau)
      A(5,2)  = -ifac*x
      A(5,7)  = -ifac*x*c_beta**2/ds**2
      A(5,10) = -ifac*Q_e/ds**2
      A(5,11) = ifac*Q/ds**2
      A(5,12) = 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,12) = 1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(5,13) = -1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,13) = -1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(6,7)  = ifac*c_beta**2/ds**2*x
      A(6,10) = ifac*Q_e/ds**2
      A(6,11) = -ifac*Q/ds**2
      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
     $        - 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
!     $        - 1.0/(1.0+tau)*c_beta**2/ds**2
      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
     $        + 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
!     $        + 1.0/(1.0+tau)*c_beta**2/ds**2
      A(7,14) = 1.0/Pr
      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,7
        do j=1,7
          R(i,j) = R1D(i+7*(j-1))
        end do
      end do
      
      U11 = A(:,1:7) + matmul(A(:,8:14), R)
      dy = matmul(U11, y)
      return 
      END

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in left region
c     pd (output) : analytic jacobian tensor d(dR/dX)/dR (Eqn. 44)
c-----------------------------------------------------------------------
      SUBROUTINE dw_der_wl(neq,x,y,ml,mu,pd,nrpd,ipar)

      INTEGER, INTENT(IN) :: neq,ml,mu,nrpd,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j,k,l
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq,neq), INTENT(OUT) :: pd
      COMPLEX(r8), DIMENSION(5,5) :: ytemp
      COMPLEX(r8), DIMENSION(5,5,5,5) :: pdtemp
      COMPLEX(r8), DIMENSION(10,10) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A = 0.0
      A(1,1)=1.0
      A(2,8)=1.0
      A(3,1)=(Q-Q_e-tau*Q_e)*x/pr-tau*ds**4/c_beta**4
     $      *(Q-Q_e)*x**3-ifac*tau*ds**2/c_beta**2*(Q-Q_e)**2
     $      *x-ifac*tau*Q*ds**2/c_beta**4*(Q-Q_e)*x
      A(3,2)=4.0*tau*ds**2/c_beta**2*x
      A(3,4)=-4.0*tau*ds**2/c_beta**2*x
      A(3,5)=-tau*Q/c_beta**2*x-tau*Q/pr*x+ifac*tau*ds**2
     $      /c_beta**2*x**3
      A(3,6)=-2.0*tau*ds**2/c_beta**2*(Q-Q_e)
      A(3,7)=-tau*Q*Q_e/c_beta**4+2.0*tau*ds**2/c_beta**2
     $      -x**2/pr+tau*ds**4/c_beta**4*x**4+ifac*tau*ds**2
     $      /c_beta**2*(Q-Q_e)*x**2+ifac*tau*Q*ds**2/c_beta**4
     $      *x**2+ifac*tau*Q_e*ds**2/c_beta**4*x**2
      A(3,8)=tau*ds**2/c_beta**2*x**2+ifac*(Q-Q_i)/pr+ifac
     $      *tau*Q_e/c_beta**2
      A(3,9)=tau*Q**2/c_beta**4-2.0*tau*ds**2/c_beta**2+(1.0
     $      +tau)/pr*x**2-tau*ds**4/c_beta**4*x**4-ifac*tau
     $      *ds**2/c_beta**2*(Q-Q_e)*x**2-ifac*2.0*tau*Q*ds**2
     $      /c_beta**4*x**2
      A(3,10)=ifac*2.0*tau
      A(4,1)=ds**2/c_beta**2*(Q-Q_e)*x
      A(4,5)=-ifac*x
      A(4,7)=-ds**2/c_beta**2*x**2-ifac*Q_e/c_beta**2
      A(4,9)=ds**2/c_beta**2*x**2+ifac*Q/c_beta**2
      A(5,10)=1.0
      A(6,1)=ifac*(Q-Q_e)
      A(6,7)=-ifac*x
      A(6,9)=ifac*x
      A(7,2)=1.0
      A(8,3)=1.0
      A(9,4)=1.0
      A(10,1)=ifac*Q_e/pr
      A(10,5)=ifac*Q/pr
      A(10,9)=-ifac*x/pr
      
      DO i=1,5
        DO j=1,5
          ytemp(i,j)=y(i+5*(j-1)) 
        END DO
      END DO

      DO i=1,5
        DO j=1,5
          DO k=1,5
            DO l=1,5
              pdtemp(i,j,k,l)=0.0
              if (j==l) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(5+i,5+k)
     $        +ytemp(i,1)*A(1,5+k)+ytemp(i,2)*A(2,5+k)+ytemp(i,3)
     $        *A(3,5+k)+ytemp(4,5+k)*A(4,5+k)+ytemp(i,5)*A(5,5+k)
              if (i==k) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(l,j)
     $        +A(l,6)*ytemp(1,j)+A(l,7)*ytemp(2,j)+A(l,8)
     $        *ytemp(3,j)+A(l,9)*ytemp(4,j)+A(l,10)*ytemp(5,j)
              pd(i+5*(j-1),k+5*(l-1))=pdtemp(i,j,k,l)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE dw_der_wl

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in right region
c     pd (output) : analytic jacobian tensor d(dR/dX)/dR (Eqn. 44)
c-----------------------------------------------------------------------
      SUBROUTINE dw_der_wr(neq,x,y,ml,mu,pd,nrpd,ipar)

      INTEGER, INTENT(IN) :: neq,ml,mu,nrpd,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j,k,l
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq,neq), INTENT(OUT) :: pd
      COMPLEX(r8), DIMENSION(5,5) :: ytemp
      COMPLEX(r8), DIMENSION(5,5,5,5) :: pdtemp
      COMPLEX(r8), DIMENSION(10,10) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,6)=1.0
      A(2,7)=1.0
      A(3,2)=(Q_e*tau+Q)/c_beta**2*ds**2*x
      A(3,5)=-ifac*(Q_e*tau+Q)*x
      A(3,8)=ifac*(Q_e*tau+Q)/c_beta**2
      A(3,9)=-Q_e
      A(4,2)=x/pr-ifac*(Q-Q_i)/pr*ds**2/c_beta**2*tau*x
      A(4,5)=-1.0/pr*(Q-Q_i)*tau*x
      A(4,8)=1.0/pr*(Q-Q_i)*tau/c_beta**2
      A(4,9)=ifac*(Q-Q_i)/pr
      A(5,10)=1.0
      A(6,2)=ifac
      A(7,1)=2.0*(Q-Q_e)/x**2
      A(7,2)=-2.0/x**2+(tau+1.0)*ds**2/c_beta**2*x**2+ifac*(Q-Q_e)
      A(7,5)=-ifac*(tau+1.0)*x**2
      A(7,6)=-2.0*(Q-Q_e)/x
      A(7,7)=2.0/x
      A(7,8)=ifac*(tau+1.0)*x/c_beta**2
      A(7,9)=-x
      A(8,3)=1.0
      A(9,4)=1.0
      A(10,2)=ifac/pr/(Q-Q_e)*Q_e
      A(10,5)=ifac*Q/pr
      A(10,8)=-ifac/pr/(Q-Q_e)*x
     
      DO i=1,5
        DO j=1,5
          ytemp(i,j)=y(i+5*(j-1)) 
        END DO
      END DO

      DO i=1,5
        DO j=1,5
          DO k=1,5
            DO l=1,5
              pdtemp(i,j,k,l)=0.0
              if (j==l) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(5+i,5+k)
     $        +ytemp(i,1)*A(1,5+k)+ytemp(i,2)*A(2,5+k)+ytemp(i,3)
     $        *A(3,5+k)+ytemp(4,5+k)*A(4,5+k)+ytemp(i,5)*A(5,5+k)
              if (i==k) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(l,j)
     $        +A(l,6)*ytemp(1,j)+A(l,7)*ytemp(2,j)+A(l,8)
     $        *ytemp(3,j)+A(l,9)*ytemp(4,j)+A(l,10)*ytemp(5,j)
              pd(i+5*(j-1),k+5*(l-1))=pdtemp(i,j,k,l)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE dw_der_wr

      SUBROUTINE dy_der_yr(neq,x,y,ml,mu,pd,nrpd,ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ml,mu,nrpd,ipar
      REAL(r8), INTENT(IN) :: x
      REAL(r8) :: Rtemp, dRtemp
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq,neq), INTENT(OUT) :: pd
      COMPLEX(r8), DIMENSION(5,5) :: R, U11
      COMPLEX(r8), DIMENSION(5,10) :: A
      COMPLEX(r8), DIMENSION(25) :: R1D
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,6) = 1.0
      A(2,7) = 1.0
      A(3,2) = (Q_e*tau+Q)/c_beta**2*ds**2*x
      A(3,5) = -ifac*(Q_e*tau+Q)*x
      A(3,8) = ifac*(Q_e*tau+Q)/c_beta**2
      A(3,9) = -Q_e
      A(4,2) = x/Pr-ifac*(Q-Q_i)/Pr*ds**2/c_beta**2*tau*x
      A(4,5) = -1.0/Pr*(Q-Q_i)*tau*x
      A(4,8) = 1.0/Pr*(Q-Q_i)*tau/c_beta**2
      A(4,9) = ifac*(Q-Q_i)/Pr
      A(5,10) = 1.0
      
      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,5
        do j=1,5
          R(i,j) = R1D(i+5*(j-1))
        end do
      end do
      
      U11 = A(:,1:5) + matmul(A(:,6:10), R)
      pd = transpose(U11)
      END

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on full model.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix
c     pd (output) : analytic jacobian tensor d(dR/dX)/dR
c-----------------------------------------------------------------------
      SUBROUTINE dw_der_w_full(neq,x,y,ml,mu,pd,nrpd,ipar)

      INTEGER, INTENT(IN) :: neq,ml,mu,nrpd,ipar
      REAL(r8), INTENT(IN) :: x
      INTEGER :: i,j,k,l
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq,neq), INTENT(OUT) :: pd
      COMPLEX(r8), DIMENSION(7,7) :: ytemp
      COMPLEX(r8), DIMENSION(7,7,7,7) :: pdtemp
      COMPLEX(r8), DIMENSION(14,14) :: A
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      parflow_flag=.TRUE.
      IF (parflow_flag) THEN
        A=0.0
        A(1,8)  = 1.0
        A(2,9)  = c_beta**2/ds**2/Pe
        A(3,12) = tau/(1.0+tau)
        A(3,13) = 1.0/(1.0+tau)
        A(4,12) = -1.0/(1.0+tau)
        A(4,13) = 1.0/(1.0+tau)
        A(5,2)  = -ifac*x
        A(5,7)  = -ifac*x*c_beta**2/ds**2
        A(5,10) = -ifac*Q_e/ds**2
        A(5,11) = ifac*Q/ds**2
        A(5,12) = 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
        !A(5,12) = 1.0/(1.0+tau)*(c_beta**2)/ds**2
        A(5,13) = -1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
        !A(5,13) = -1.0/(1.0+tau)*(c_beta**2)/ds**2
        A(6,7)  = ifac*c_beta**2/ds**2*x
        A(6,10) = ifac*Q_e/ds**2
        A(6,11) = -ifac*Q/ds**2
        A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
     $          - 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!        A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
!     $          - 1.0/(1.0+tau)*c_beta**2/ds**2
        A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
     $          + 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!        A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
!     $          + 1.0/(1.0+tau)*c_beta**2/ds**2
        A(7,14) = 1.0/Pr
        A(8,2)  = 1.0
        A(9,1)  = -ifac/(1.0+tau)*((Q-Q_e)+(Q-Q_i)*Pe/Pr)
        A(9,2)  = (1.0+Pe/Pr)/(1.0+tau)
        A(9,7)  = -ifac*Q*Pe/Pr
        A(9,10) = ifac*(1.0+Pe/Pr)/(1.0+tau)*x
        A(9,11) = ifac*(tau*Pe/Pr-1.0)/(1.0+tau)*x
        A(10,3) = 1.0
        A(11,4) = 1.0
        A(12,5) = 1.0/Pe
        A(13,6) = 1.0/Pr
        A(14,1) = ifac*(Q-Q_i)/(1.0+tau)
        A(14,2) = -1.0/(1.0+tau)
        A(14,7) = ifac*Q
        A(14,10)= -ifac*x/(1.0+tau)
        A(14,11)= -ifac*tau*x/(1.0+tau)
      ENDIF 
      
      DO i=1,7
        DO j=1,7
          ytemp(i,j)=y(i+7*(j-1)) 
        END DO
      END DO

      DO i=1,7
        DO j=1,7
          DO k=1,7
            DO l=1,7
              pdtemp(i,j,k,l)=0.0
              if (j==l) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(7+i,7+k)
     $        +ytemp(i,1)*A(1,7+k)+ytemp(i,2)*A(2,7+k)+ytemp(i,3)
     $        *A(3,7+k)+ytemp(4,7+k)*A(4,7+k)+ytemp(i,5)*A(5,7+k)
     $        +ytemp(i,6)*A(6,7+k)+ytemp(i,7)*A(7,7+k)
              if (i==k) pdtemp(i,j,k,l)=pdtemp(i,j,k,l)+A(l,j)
     $        +A(l,8)*ytemp(1,j)+A(l,9)*ytemp(2,j)+A(l,10)
     $        *ytemp(3,j)+A(l,11)*ytemp(4,j)+A(l,12)*ytemp(5,j)
     $        +ytemp(l,13)*ytemp(6,j)+ytemp(l,14)*ytemp(7,j)
              pd(i+7*(j-1),k+7*(l-1))=pdtemp(i,j,k,l)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE dw_der_w_full

      SUBROUTINE dy_der_y_full(neq,x,y,ml,mu,pd,nrpd,ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ml,mu,nrpd,ipar
      REAL(r8), INTENT(IN) :: x
      REAL(r8) :: Rtemp, dRtemp
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq,neq), INTENT(OUT) :: pd
      COMPLEX(r8), DIMENSION(7,7) :: R, U11
      COMPLEX(r8), DIMENSION(7,14) :: A
      COMPLEX(r8), DIMENSION(49) :: R1D
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      A=0.0
      A(1,8)  = 1.0
      A(2,9)  = c_beta**2/ds**2/Pe
      A(3,12) = tau/(1.0+tau)
      A(3,13) = 1.0/(1.0+tau)
      A(4,12) = -1.0/(1.0+tau)
      A(4,13) = 1.0/(1.0+tau)
      A(5,2)  = -ifac*x
      A(5,7)  = -ifac*x*c_beta**2/ds**2
      A(5,10) = -ifac*Q_e/ds**2
      A(5,11) = ifac*Q/ds**2
      A(5,12) = 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,12) = 1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(5,13) = -1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
      !A(5,13) = -1.0/(1.0+tau)*(c_beta**2)/ds**2
      A(6,7)  = ifac*c_beta**2/ds**2*x
      A(6,10) = ifac*Q_e/ds**2
      A(6,11) = -ifac*Q/ds**2
      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
     $          - 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,12) = ifac*tau/(1.0+tau)*(Q-Q_i) 
!     $        - 1.0/(1.0+tau)*c_beta**2/ds**2
      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
     $          + 1.0/(1.0+tau)*(c_beta**2+(1-c_beta**2)*Kp)/ds**2
!      A(6,13) = ifac*1.0/(1.0+tau)*(Q-Q_i) 
!     $        + 1.0/(1.0+tau)*c_beta**2/ds**2
      A(7,14) = 1.0/Pr
      
      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,7
        do j=1,7
          R(i,j) = R1D(i+7*(j-1))
        end do
      end do
      
      U11 = A(:,1:7) + matmul(A(:,8:14), R)
      pd = transpose(U11)
      END

      subroutine Init_yr1 (neq, x, y, y1)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: y1
      INTEGER :: i,j,info
      COMPLEX(r8), DIMENSION(5,5) :: R, B21, B22, B21B22R
      COMPLEX(r8), DIMENSION(5) :: beta
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      INTEGER, DIMENSION(5):: ipiv

      DO i=1,5
        beta(i)=0.0
        DO j=1,5
          B21(i,j)=0.0
          B22(i,j)=0.0
          R(i,j)=y(i+5*(j-1))
          B21B22R(i,j)=0.0
        END DO
      END DO

      beta(1) = 1.0
      B22(1,1) = 1.0
      B21(2,2) = 1.0
      B22(3,2) = 1.0
      B22(4,3) = 1.0
      B21(5,5) = 1.0

      B21B22R=B21+matmul(B22,R)
      CALL zgetrf(5,5,B21B22R,5,ipiv,info)
      CALL zgetrs('N',5,1,B21B22R,5,ipiv,beta,5,info)
      y1=beta
      END

      subroutine Init_y1_full (neq, x, y, y1)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: y1
      INTEGER :: i,j,info
      REAL(r8) :: lambda
      COMPLEX(r8), DIMENSION(4,7) :: Bcon1, Bcon2, Bcon, Bconn
      COMPLEX(r8), DIMENSION(14,7) :: Bsql1, Bsql2, Bsql, Bsqln
      REAL(r8), DIMENSION(7) :: scaler
      COMPLEX(r8), DIMENSION(4) :: betacon, betaconn
      COMPLEX(r8), DIMENSION(14) :: betasql, betasqln, weights
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      COMPLEX(r8), DIMENSION(7,7) :: R

      DO i=1,7
        DO j=1,7
          R(j,i)=y(j+7*(i-1))
        END DO
        DO j=1,4
          Bcon1(j,i)=0.0
          Bcon2(j,i)=0.0
        END DO
        DO j=1,14
          Bsql1(j,i)=0.0
          Bsql2(j,i)=0.0
        END DO
      END DO
      betacon=0.0
      betasql=0.0
      weights=1.0

      ! Asymptotic constraint
!      betacon(1)=1.0
!      betacon(4)=Q-Q_e
!      Bcon2(1,1)=1.0
!      Bcon1(2,1)=(Q-Q_e)
!      Bcon2(2,3)=-x
!      Bcon2(2,4)=x
!      Bcon2(3,3)=-Q_e
!      Bcon2(3,4)=Q
!      Bcon2(4,3)=1.0
!      Bcon2(4,4)=-1.0
!      Bcon1(4,3)=x
!      Bcon1(4,4)=-x
!      !Bcon2(5,5)=1.0
!      !Bcon1(6,7)=1.0
!      !Bcon2(7,7)=1.0

      ! Equation-based constraint
      betacon(1)=1.0
      Bcon2(1,1)=1.0
      Bcon1(2,1)=ifac*(Q-Q_e)
      Bcon1(2,2)=-1.0
      Bcon2(2,3)=-ifac*x
      Bcon2(2,4)=ifac*x
      !Bcon1(2,1)=(Q-Q_e)+Pe/Pr*(Q-Q_i)
      !Bcon1(2,2)=ifac*(1.0+Pe/Pr)
      !Bcon2(2,3)=-(1.0+Pe/Pr)*x
      !Bcon2(2,4)=(1.0-Pe/Pr*tau)*x
      !Bcon1(2,7)=(1.0+tau)*Q*Pe/Pr
      Bcon2(3,3)=Q_e
      Bcon2(3,4)=-Q
      Bcon1(3,2)=x*ds**2
      Bcon1(3,7)=x*c_beta**2
      Bcon2(3,5)=-(c_beta**2+(1.0-c_beta**2)*Kp)/(1.0+tau)
      Bcon2(3,6)=(c_beta**2+(1.0-c_beta**2)*Kp)/(1.0+tau)
      Bcon1(4,2)=x
      Bcon2(4,5)=-(Q-Q_i)*tau*(c_beta**2+(1.0-c_beta**2)*Kp)
      Bcon2(4,6)=-(Q-Q_i)*(c_beta**2+(1.0-c_beta**2)*Kp)
!      Bcon1(3,2)=-ifac*x
!      Bcon2(3,3)=-ifac*Q_e/ds**2
!      Bcon2(3,4)=ifac*Q/ds**2
!      Bcon1(3,7)=-ifac*c_beta**2/ds**2*x
!      Bcon2(3,5)=(c_beta**2+(1.0-c_beta**2)*Kp)/(1.0+tau)/ds**2
!      Bcon2(3,6)=-(c_beta**2+(1.0-c_beta**2)*Kp)/(1.0+tau)/ds**2
!      Bcon2(4,3)=ifac*Q_e/ds**2
!      Bcon2(4,4)=-ifac*Q/ds**2
!      Bcon1(4,7)=ifac*c_beta**2/ds**2*x
!      Bcon2(4,5)=-(c_beta**2+(1.0-c_beta**2)*Kp)/(1.0+tau)/ds**2 
!     $          +ifac*tau/(1.0+tau)*(Q-Q_i)
!      Bcon2(4,6)=(c_beta**2+(1-c_beta**2)*Kp)/(1.0+tau)/ds**2
!     $          +ifac/(1.0+tau)*(Q-Q_i)
      !Bcon1(5,1)=Q_e
      !Bcon2(5,4)=-x
      !Bcon1(5,7)=Q
      !Bcon1(5,1)=ifac*(Q-Q_i)
      !Bcon1(5,2)=-1.0
      !Bcon2(5,3)=-ifac*x
      !Bcon2(5,4)=-ifac*tau*x
      !Bcon1(5,7)=ifac*Q*(1+tau)
      Bcon=Bcon1+matmul(Bcon2,R)

      lambda = 10.0
!      betasql(10)=Q-Q_e
      betasql(12)=Q
      betasql(13)=Q_e
      Bsql1(1,2)=1.0
      Bsql2(2,2)=1.0/(Pe*ds**2/c_beta**2)
      Bsql1(3,3)=-Q_e
      Bsql1(3,4)=Q
      Bsql2(4,5)=1.0
      Bsql1(5,5)=1.0/Pe
      Bsql2(6,6)=1.0
      Bsql1(7,6)=1.0/Pr
      Bsql1(8,7)=1000.0
      Bsql2(9,7)=1000.0/Pr
      Bsql1(10,1)=Q
      Bsql2(10,3)=-x
      Bsql1(11,1)=Q_e
      Bsql2(11,4)=-x
      Bsql1(12,3)=x
      Bsql2(12,3)=1.0
      Bsql1(13,4)=x
      Bsql2(13,4)=1.0
      Bsql1(14,1)=Q_e
      Bsql2(14,4)=-x
      Bsql1(14,7)=Q
      !Bsql2(14,3)=-Q_e
      !Bsql2(14,4)=Q
      !Bsql1(14,7)=-3.0*Pr*1e3
      !Bsql2(14,7)=x*1e3
      Bsql1(8:,:)=Bsql1(8:,:)*lambda
      Bsql2(8:,:)=Bsql2(8:,:)*lambda
      betasql(8:)=betasql(8:)*lambda
!      Bsql2(10,3)=lambda*1.0
!      Bsql2(10,4)=-lambda*1.0
!      Bsql1(10,3)=lambda*x
!      Bsql1(10,4)=-lambda*x
!      Bsql2(11,3)=-lambda*Q_e
!      Bsql2(11,4)=lambda*Q
!      Bsql1(12,1)=lambda*(Q-Q_e)
!      Bsql2(12,3)=-lambda*x
!      Bsql2(12,4)=lambda*x
      Bsql=Bsql1+matmul(Bsql2,R)
      !weights(1)=1/Pe*100000
      !weights(7)=1/Pr
      !weights(9)=1/Pr

      ! row-wise normalization
      do i=1,7
         !scaler(i) = sqrt(sum(abs(Bcon(:,i))**2))
         scaler(i) = sqrt(sum(abs(Bsql(:,i))**2))
         if (scaler(i) > 1.0e-14) then
            Bconn(:,i) = Bcon(:,i) / scaler(i)
            Bsqln(:,i) = Bsql(:,i) / scaler(i)
            !betaconn(i) = betacon(i) / scaler(i)
            !betasqln(i) = betasql(i) / scaler(i)
         else
            scaler(i) = 1.0
            Bconn(:,i) = Bcon(:,i)
            Bsqln(:,i) = Bsql(:,i)
            !betaconn(i) = betacon(i)
            !betasqln(i) = betasql(i)
         end if
      end do

!      ! coloumn-wise normalization
!      do i=1,7
!         !scaler(i) = sqrt(sum(abs(Bcon(:,i))**2))
!         scaler(i) = sqrt(sum(abs(Bsql(i,:))**2))
!         if (scaler(i) > 1.0e-14) then
!            Bconn(i,:) = Bcon(i,:) / scaler(i)
!            Bsqln(i,:) = Bsql(i,:) / scaler(i)
!            betaconn(i) = betacon(i) / scaler(i)
!            betasqln(i) = betasql(i) / scaler(i)
!         else
!            scaler(i) = 1.0
!            Bconn(i,:) = Bcon(i,:)
!            Bsqln(i,:) = Bsql(i,:)
!            betaconn(i) = betacon(i)
!            betasqln(i) = betasql(i)
!         end if
!      end do

      write(*,*) 'Bconn norm is ', norm2(abs(Bconn))
      write(*,*) 'Bsqln norm is ', norm2(abs(Bsqln))
!      do i=1,9
!         Bsqln(i,:) = Bsqln(i,:) * weights(i)
!         betasql(i) = betasql(i) * weights(i)
!      end do

      lambda = 100.0*x
!      CALL constrained_ls_socp_admm(Bsqln,Bconn,betasql,betacon,y1,
!     $                        1.0e-8_r8,
!     $                        1.0_r8,
!     $                        100000,
!     $                        1.0e-6_r8)
      CALL constrained_ls_kkt(Bsqln,Bconn,betasql,betacon,y1)
!      CALL soft_constrained_ls_svd(Bsqln,Bconn,betasql,betacon,y1, 
!     $                             lambda)
      y1 = y1 / scaler

      END

      subroutine get_yr2 (neq, x, y, dy, ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ipar
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      REAL(r8) :: Rtemp, dRtemp
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(5,5) :: R
      COMPLEX(r8), DIMENSION(25) :: R1D

      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,5
        do j=1,5
          R(i,j) = R1D(i+5*(j-1))
        end do
      end do
      
      dy = matmul(R, y)
      return 
      END

      subroutine get_y2_full (neq, x, y, dy, ipar)
      USE sglobal_mod

      INTEGER, INTENT(IN) :: neq,ipar
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      REAL(r8) :: Rtemp, dRtemp
      INTEGER :: i,j
      COMPLEX(r8), DIMENSION(7,7) :: R
      COMPLEX(r8), DIMENSION(49) :: R1D

      CALL cspline_eval(splR,x,0)
      R1D=splR%f
      
      do i=1,7
        do j=1,7
          R(i,j) = R1D(i+7*(j-1))
        end do
      end do
      
      dy = matmul(R, y)
      return 
      END

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in left region
c     y (output)  : riccati matrix in right region
c     Connect left and right regions by transforming R
c     See Eqns. 59, A.25-A.54.
c-----------------------------------------------------------------------
      SUBROUTINE Transform_R(y,x,neq)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(INOUT) :: y
      INTEGER :: i,j,info
      COMPLEX(r8), DIMENSION(5,5) :: R, T22RT12, RT11T21, invT22RT12
      COMPLEX(r8), DIMENSION(10,10) :: Tx
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      INTEGER, DIMENSION(5):: ipiv

      Tx=0.0
      Tx(1,1)=1.0
      Tx(2,1)=-Q/x**2
      Tx(2,2)=Q/(Q-Q_e)/x**2
      Tx(2,3)=1.0/(Q-Q_e)
      Tx(2,6)=Q/x
      Tx(2,7)=-Q/(Q-Q_e)/x
      Tx(3,2)=-tau*ds**2/c_beta**2
      Tx(3,3)=-ifac*tau/c_beta**2
      Tx(3,4)=1.0
      Tx(3,5)=ifac*tau
      Tx(3,7)=-tau*ds**2/c_beta**2*x
      Tx(3,10)=ifac*tau*x
      Tx(4,1)=-Q_e/x**2
      Tx(4,2)=Q_e/(Q-Q_e)/x**2
      Tx(4,3)=1.0/(Q-Q_e)
      Tx(4,6)=Q_e/x
      Tx(4,7)=-Q_e/(Q-Q_e)/x
      Tx(5,5)=1.0
      Tx(6,6)=1.0
      Tx(7,1)=Q/x
      Tx(7,2)=-Q/(Q-Q_e)/x
      Tx(7,8)=1.0/(Q-Q_e)
      Tx(8,2)=-tau*ds**2/c_beta**2*x
      Tx(8,8)=-ifac*tau/c_beta**2
      Tx(8,9)=1.0
      Tx(8,10)=ifac*tau*x
      Tx(9,1)=Q_e/x
      Tx(9,2)=-Q_e/(Q-Q_e)/x
      Tx(9,8)=1.0/(Q-Q_e)
      Tx(10,10)=1.0

      DO i=1,5
        DO j=1,5
          R(i,j)=y(i+5*(j-1))
          invT22RT12(i,j)=0.0
        END DO
      END DO

      T22RT12=Tx(6:10,6:10)-matmul(R,Tx(1:5,6:10))
      RT11T21=matmul(R,Tx(1:5,1:5))-Tx(6:10,1:5)
      DO i=1,5
        invT22RT12(i,i)=1.0
      ENDDO
      CALL zgetrf(5,5,T22RT12,5,ipiv,info)
      CALL zgetrs('N',5,5,T22RT12,5,ipiv,invT22RT12,5,info)
      
      R=matmul(invT22RT12,RT11T21)

      DO i=1,5
        DO j=1,5
          y(i+5*(j-1))=R(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in right region
c     Delta_new (output)  : inner layer Delta
c     See Eqns. 46, 47, 65, 66
c-----------------------------------------------------------------------
      SUBROUTINE Update_Delta(Delta_new,y,x,neq)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), INTENT(OUT) :: Delta_new
      INTEGER :: i,j,info
      COMPLEX(r8), DIMENSION(5,5) :: R, B21, B22, B21B22R
      COMPLEX(r8), DIMENSION(5) :: beta, y1
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      INTEGER, DIMENSION(5):: ipiv

      DO i=1,5
        beta(i)=0.0
        DO j=1,5
          B21(i,j)=0.0
          B22(i,j)=0.0
          R(i,j)=y(i+5*(j-1))
        END DO
      END DO

      beta(1)=1.0
      B22(1,1)=1.0
      B21(2,2)=1.0
      B22(3,2)=1.0
      B22(4,3)=1.0
      B21(5,5)=1.0

      B21B22R=B21+matmul(B22,R)
      CALL zgetrf(5,5,B21B22R,5,ipiv,info)
      CALL zgetrs('N',5,1,B21B22R,5,ipiv,beta,5,info)
      y1=beta
      Delta_new=2.0/(y1(1)-x)
      RETURN
      END SUBROUTINE

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on the full model.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in right region
c     Delta_new (output)  : inner layer Delta
c-----------------------------------------------------------------------
      SUBROUTINE Update_Delta_full(Delta_new,y,x,neq)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), INTENT(OUT) :: Delta_new
      COMPLEX(r8), DIMENSION(neq) :: y1

      CALL Init_y1_full (neq, x, y, y1)
      Delta_new=2.0/(y1(1)-x)
      RETURN
      END SUBROUTINE

      subroutine constrained_ls_svd(B1, B2, beta1, beta2, y1)
      implicit none
    
      !----------------------------------
      ! input / output
      !----------------------------------
      complex(r8), intent(in)  :: B1(:,:), B2(:,:), beta1(:), beta2(:)
      complex(r8), intent(out) :: y1(:)
    
      !----------------------------------
      ! local variables
      !----------------------------------
      integer :: m1, m2, n
      integer :: lda, ldu, ldvt, info, lwork
      integer :: i, r
    
      complex(r8), allocatable :: U(:,:), VT(:,:), work(:)
      complex(r8), allocatable :: Vn(:,:), yp(:), z(:), Bred_copy(:,:)
      complex(r8), allocatable :: rhs(:), Bred(:,:), B2_copy(:,:)
    
      real(r8), allocatable :: S(:), rwork(:)
    
      !----------------------------------
      ! dimensions
      !----------------------------------
      m1 = size(B1,1)
      m2 = size(B2,1)
      n  = size(B2,2)
    
      allocate(B2_copy(m2,n))
      B2_copy = B2

      !lda  = max(m2,n)
      lda  = m2
      ldu  = m2
      ldvt = n
    
      !----------------------------------
      ! SVD of B2
      !----------------------------------
      allocate(U(ldu,m2), VT(ldvt,n))
      allocate(S(min(m2,n)))
      allocate(rwork(5*min(m2,n)))
    
      lwork = -1
      allocate(work(1))
      call zgesvd('S','S', m2, n, B2_copy, lda, S, U, ldu, VT, ldvt,
     $              work, lwork, rwork, info)
    
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))
      B2_copy = B2
    
      call zgesvd('S','S', m2, n, B2_copy, lda, S, U, ldu, VT, ldvt,
     $              work, lwork, rwork, info)
    
      if (info /= 0) stop 'zgesvd failed on B2'
    
      !----------------------------------
      ! rank determination
      !----------------------------------
      r = 0
      do i = 1, size(S)
         if (S(i) > 1.0e-14) r = r + 1
      end do
    
      !----------------------------------
      ! null space Vn
      !----------------------------------
      allocate(Vn(n, n-r))
      !Vn = transpose(conjg(VT(r+1:n, :)))
      do i = 1, n-r
         Vn(:,i) = conjg(VT(r+i,:))
      end do
    
      !----------------------------------
      ! particular solution yp = B2^+ b2
      !----------------------------------
      allocate(yp(n))
      yp = (0.0, 0.0)
    
      do i = 1, r
         yp = yp + ( dot_product(conjg(U(:,i)), beta2) / S(i) ) *
     $               conjg(VT(i,:))
      end do
    
      !----------------------------------
      ! reduced least squares
      !----------------------------------
      allocate(rhs(m1))
      rhs = beta1 - matmul(B1, yp)
    
      allocate(Bred(m1, n-r))
      Bred = matmul(B1, Vn)
    
      !----------------------------------
      ! SVD of reduced system
      !----------------------------------
      allocate(Bred_copy(m1, n-r))
      Bred_copy = Bred
      deallocate(U, VT, S, rwork)
      allocate(U(m1,m1), VT(n-r,n-r))
      allocate(S(min(m1,n-r)))
      allocate(rwork(5*min(m1,n-r)))
    
      lwork = -1
      deallocate(work)
      allocate(work(1))
      call zgesvd('S','S', m1, n-r, Bred_copy, m1, S, U, m1, VT, n-r,
     $              work, lwork, rwork, info)
    
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))
      Bred_copy = Bred
    
      call zgesvd('S','S', m1, n-r, Bred_copy, m1, S, U, m1, VT, n-r,
     $              work, lwork, rwork, info)
    
      if (info /= 0) stop 'zgesvd failed on reduced system'
    
      !----------------------------------
      ! solve for z
      !----------------------------------
      allocate(z(n-r))
      z = (0.0, 0.0)
    
      do i = 1, size(S)
         if (S(i) > 1.0e-14) then
            z = z + ( dot_product(conjg(U(:,i)), rhs) / S(i) ) *
     $               conjg(VT(i,:))
         end if
      end do
    
      !----------------------------------
      ! final solution
      !----------------------------------
      y1 = yp + matmul(Vn, z)
    
      end subroutine constrained_ls_svd

      subroutine soft_constrained_ls_svd(B1,B2,beta1,beta2,y1,lambda)
      implicit none
    
      !----------------------------------
      ! input / output
      !----------------------------------
      complex(r8), intent(in)  :: B1(:,:), B2(:,:), beta1(:), beta2(:)
      real(r8), intent(in)     :: lambda  ! penalty parameter
      complex(r8), intent(out) :: y1(:)
    
      !----------------------------------
      ! local variables
      !----------------------------------
      integer :: m1, m2, n
      integer :: lda, ldu, ldvt, info, lwork
      integer :: i
    
      complex(r8), allocatable :: A_aug(:,:), b_aug(:)
      complex(r8), allocatable :: U(:,:), VT(:,:), work(:)
      complex(r8), allocatable :: A_aug_copy(:,:)
    
      real(r8), allocatable :: S(:), rwork(:)
    
      !----------------------------------
      ! dimensions
      !----------------------------------
      m1 = size(B1,1)
      m2 = size(B2,1)
      n  = size(B2,2)
    
      !----------------------------------
      ! Augmented system construction
      ! Minimize: ||B1*y - beta1||^2 + lambda * ||B2*y - beta2||^2
      ! 
      ! This is equivalent to:
      ! [ B1           ]     [ beta1 ]
      ! [ sqrt(lambda)*B2 ] * y = [ sqrt(lambda)*beta2 ]
      !----------------------------------
      allocate(A_aug(m1+m2, n))
      allocate(b_aug(m1+m2))
    
      ! Upper block: B1
      A_aug(1:m1, :) = B1
      b_aug(1:m1) = beta1
    
      ! Lower block: sqrt(lambda)*B2
      A_aug(m1+1:m1+m2, :) = sqrt(lambda) * B2
      b_aug(m1+1:m1+m2) = sqrt(lambda) * beta2
    
      write(*,*) 'Soft constraint: lambda =', lambda
      write(*,*) 'Augmented system size:', m1+m2, 'x', n
    
      !----------------------------------
      ! SVD of augmented system
      !----------------------------------
      allocate(A_aug_copy(m1+m2, n))
      A_aug_copy = A_aug
    
      lda  = m1 + m2
      ldu  = m1 + m2
      ldvt = n
    
      allocate(U(ldu, m1+m2), VT(ldvt, n))
      allocate(S(min(m1+m2, n)))
      allocate(rwork(5*min(m1+m2, n)))
    
      ! Query optimal work size
      lwork = -1
      allocate(work(1))
      call zgesvd('S','S', m1+m2, n, A_aug_copy, lda, S, U, ldu, VT, 
     $              ldvt, work, lwork, rwork, info)
    
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))
    
      ! Actual SVD computation
      call zgesvd('S','S', m1+m2, n, A_aug, lda, S, U, ldu, VT, 
     $              ldvt, work, lwork, rwork, info)
    
      if (info /= 0) stop 'zgesvd failed on augmented system'
    
      write(*,*) 'SVD completed successfully'
      write(*,*) 'Singular values:', S
    
      !----------------------------------
      ! Solve using pseudoinverse
      ! y = A^+ * b_aug
      !----------------------------------
      y1 = (0.0, 0.0)
    
      do i = 1, size(S)
         if (S(i) > 1.0e-12) then
            y1 = y1 + ( dot_product(conjg(U(:,i)), b_aug) / S(i) ) *
     $                   conjg(VT(i,:))
         end if
      end do
    
      !----------------------------------
      ! Print residuals for diagnostics
      !----------------------------------
      write(*,*) 'Residual ||B1*y - beta1||:', 
     $           sqrt(sum(abs(matmul(B1,y1) - beta1)**2))
      write(*,*) 'Constraint ||B2*y - beta2||:', 
     $           sqrt(sum(abs(matmul(B2,y1) - beta2)**2))
    
      end subroutine soft_constrained_ls_svd

      subroutine constrained_ls_kkt(B1, B2, beta1, beta2, y1)
      implicit none
    
      !----------------------------------
      ! input / output
      !----------------------------------
      complex(r8), intent(in)  :: B1(:,:), B2(:,:), beta1(:), beta2(:)
      complex(r8), intent(out) :: y1(:)
    
      !----------------------------------
      ! local variables
      !----------------------------------
      integer :: m1, m2, n
      integer :: i, j, kkt_size, info
      integer, allocatable :: ipiv(:)
    
      complex(r8), allocatable :: KKT_matrix(:,:), KKT_rhs(:)
      complex(r8), allocatable :: KKT_sol(:), lambda(:)
      complex(r8), allocatable :: B1H(:,:), B2H(:,:)
      complex(r8), allocatable :: B1HB1(:,:), B1H_beta1(:)
    
      !----------------------------------
      ! dimensions
      !----------------------------------
      m1 = size(B1,1)  ! number of least squares equations
      m2 = size(B2,1)  ! number of constraints
      n  = size(B2,2)  ! number of unknowns
    
      !----------------------------------
      ! Setup KKT system:
      ! [ B1^H * B1    B2^H ] [ y1     ]   [ B1^H * beta1 ]
      ! [ B2           0    ] [ lambda ] = [ beta2        ]
      !----------------------------------
      kkt_size = n + m2
      allocate(KKT_matrix(kkt_size, kkt_size))
      allocate(KKT_rhs(kkt_size))
      allocate(KKT_sol(kkt_size))
      allocate(ipiv(kkt_size))
      
      ! Compute B1^H (Hermitian transpose)
      allocate(B1H(n, m1))
      B1H = transpose(conjg(B1))
      
      ! Compute B2^H (Hermitian transpose)
      allocate(B2H(n, m2))
      B2H = transpose(conjg(B2))
      
      ! Compute B1^H * B1
      allocate(B1HB1(n, n))
      B1HB1 = matmul(B1H, B1)
      
      ! Compute B1^H * beta1
      allocate(B1H_beta1(n))
      B1H_beta1 = matmul(B1H, beta1)
      
      !----------------------------------
      ! Assemble KKT matrix
      !----------------------------------
      KKT_matrix = (0.0d0, 0.0d0)
      
      ! Top-left block: B1^H * B1
      KKT_matrix(1:n, 1:n) = B1HB1
      
      ! Top-right block: B2^H
      KKT_matrix(1:n, n+1:kkt_size) = B2H
      
      ! Bottom-left block: B2
      KKT_matrix(n+1:kkt_size, 1:n) = B2
      
      ! Bottom-right block: 0 (already initialized)
      
      !----------------------------------
      ! Assemble KKT right-hand side
      !----------------------------------
      KKT_rhs = (0.0d0, 0.0d0)
      KKT_rhs(1:n) = B1H_beta1
      KKT_rhs(n+1:kkt_size) = beta2
      
      !----------------------------------
      ! Solve KKT system using ZGESV (LU decomposition)
      !----------------------------------
      KKT_sol = KKT_rhs
      
      call zgesv(kkt_size, 1, KKT_matrix, kkt_size, ipiv, 
     $           KKT_sol, kkt_size, info)
      
      if (info /= 0) then
         write(*,*) 'Error: zgesv failed with info = ', info
         stop 'KKT system solve failed'
      end if
      
      !----------------------------------
      ! Extract solution
      !----------------------------------
      y1 = KKT_sol(1:n)
      
      ! Optional: extract Lagrange multipliers if needed
      allocate(lambda(m2))
      lambda = KKT_sol(n+1:kkt_size)
      
      ! Clean up
      deallocate(KKT_matrix, KKT_rhs, KKT_sol, ipiv)
      deallocate(B1H, B2H, B1HB1, B1H_beta1, lambda)
      
      end subroutine constrained_ls_kkt

      subroutine constrained_ls_socp_admm(B1, B2, beta1, beta2, y1,
     $                                     mu_reg, rho, max_iter, tol)
      implicit none
      
      !----------------------------------
      ! ADMM for SOCP:
      ! minimize    ||z||₂ + (μ/2)||y1||²₂
      ! subject to  B1*y1 - beta1 = z
      !             B2*y1 = beta2
      !----------------------------------
      
      complex(r8), intent(in)  :: B1(:,:), B2(:,:), beta1(:), beta2(:)
      real(r8), intent(in), optional :: mu_reg, rho, tol
      integer, intent(in), optional :: max_iter
      complex(r8), intent(out) :: y1(:)
      
      integer :: m1, m2, n, iter, max_it, info, i
      real(r8) :: mu, rho_val, tol_val, primal_res, dual_res
      real(r8) :: z_norm, z_norm_prev
      
      complex(r8), allocatable :: z(:), u(:), y1_prev(:)
      complex(r8), allocatable :: KKT_matrix(:,:), KKT_rhs(:)
      complex(r8), allocatable :: KKT_sol(:)
      integer, allocatable :: ipiv(:)
      
      complex(r8), allocatable :: B1H(:,:), B2H(:,:)
      complex(r8), allocatable :: Gram(:,:), temp(:)
      complex(r8) :: mu_cmplx, rho_cmplx
      
      !----------------------------------
      ! Parameters
      !----------------------------------
      mu = 0.0_r8
      if (present(mu_reg)) mu = mu_reg
      
      rho_val = 1.0_r8
      if (present(rho)) rho_val = rho
      
      tol_val = 1.0e-6_r8
      if (present(tol)) tol_val = tol
      
      max_it = 100
      if (present(max_iter)) max_it = max_iter
      
      mu_cmplx = cmplx(mu, 0.0_r8, kind=r8)
      rho_cmplx = cmplx(rho_val, 0.0_r8, kind=r8)
      
      !----------------------------------
      ! Dimensions
      !----------------------------------
      m1 = size(B1,1)
      m2 = size(B2,1)
      n  = size(B2,2)
      
      !----------------------------------
      ! Initialize
      !----------------------------------
      allocate(z(m1), u(m1), y1_prev(n))
      z = (0.0d0, 0.0d0)
      u = (0.0d0, 0.0d0)
      y1 = (0.0d0, 0.0d0)
      
      !----------------------------------
      ! Pre-compute matrices for y1-update
      ! [ B1^H*B1 + ρI + μI    B2^H ] [ y1     ]   [ B1^H*(z-u+beta1) ]
      ! [ B2                   0    ] [ lambda ] = [ beta2            ]
      !----------------------------------
      
      allocate(B1H(n, m1))
      allocate(B2H(n, m2))
      B1H = transpose(conjg(B1))
      B2H = transpose(conjg(B2))
      
      allocate(Gram(n, n))
      Gram = matmul(B1H, B1)
      do i = 1, n
         Gram(i,i) = Gram(i,i) + (rho_cmplx + mu_cmplx)
      end do
      
      allocate(KKT_matrix(n+m2, n+m2))
      allocate(KKT_rhs(n+m2))
      allocate(KKT_sol(n+m2))
      allocate(ipiv(n+m2))
      allocate(temp(n))
      
      KKT_matrix = (0.0d0, 0.0d0)
      KKT_matrix(1:n, 1:n) = Gram
      KKT_matrix(1:n, n+1:n+m2) = B2H
      KKT_matrix(n+1:n+m2, 1:n) = B2
      
      !----------------------------------
      ! ADMM iterations
      !----------------------------------
      write(*,*) 'Starting ADMM iterations...'
      
      do iter = 1, max_it
         y1_prev = y1
         z_norm_prev = sqrt(sum(abs(z)**2))
         
         ! 1. y1-update (solve KKT system)
         temp = matmul(B1H, z - u + beta1)
         KKT_rhs(1:n) = temp
         KKT_rhs(n+1:n+m2) = beta2
         KKT_sol = KKT_rhs
         
         call zgesv(n+m2, 1, KKT_matrix, n+m2, ipiv,
     $              KKT_sol, n+m2, info)
         
         if (info /= 0) then
            write(*,*) 'ADMM y1-update failed at iter', iter
            exit
         end if
         
         y1 = KKT_sol(1:n)
         
         ! Reset KKT matrix for next iteration
         KKT_matrix(1:n, 1:n) = Gram
         KKT_matrix(1:n, n+1:n+m2) = B2H
         KKT_matrix(n+1:n+m2, 1:n) = B2
         
         ! 2. z-update (soft thresholding / projection onto ℓ2 ball)
         temp = matmul(B1, y1) - beta1 + u
         z_norm = sqrt(sum(abs(temp)**2))
         
         if (z_norm > 1.0_r8/rho_val) then
            z = (1.0_r8 - 1.0_r8/(rho_val*z_norm)) * temp
         else
            z = (0.0d0, 0.0d0)
         end if
         
         ! 3. u-update (dual variable)
         u = u + matmul(B1, y1) - beta1 - z
         
         ! Check convergence
         primal_res = sqrt(sum(abs(matmul(B1,y1) - beta1 - z)**2))
         dual_res = rho_val * sqrt(sum(abs(z - temp)**2))
         
         if (mod(iter, 10) == 0 .or. iter == 1) then
            write(*,'(A,I4,A,ES10.3,A,ES10.3)') 
     $         ' Iter ', iter, ': primal_res=', primal_res,
     $         ', dual_res=', dual_res
         end if
         
         if (primal_res < tol_val .and. dual_res < tol_val) then
            write(*,*) 'ADMM converged at iteration', iter
            exit
         end if
      end do
      
      if (iter >= max_it) then
         write(*,*) 'Warning: ADMM reached max iterations'
      end if
      
      deallocate(z, u, y1_prev, B1H, B2H, Gram)
      deallocate(KKT_matrix, KKT_rhs, KKT_sol, ipiv, temp)
      
      end subroutine constrained_ls_socp_admm

c-----------------------------------------------------------------------
c     riccati integration with simplified version for test.
c-----------------------------------------------------------------------
      SUBROUTINE w_der_temp(neq,x,y,dy)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      dy(1)=(2.0*x/(ifac*(Q-Q_e)+x**2.0)-1.0/x)*y(1)-y(1)*y(1)/x
     $     +x*(ifac*(Q-Q_e)+x**2.0)
     $     *(-Q*(Q-Q_i)+ifac*(Q-Q_i)*(pr+c_beta**2.0)*x**2.0
     $     +pr*c_beta**2.0*x**4.0)/(ifac*(Q-Q_e)
     $     +(c_beta**2.0+ifac*(Q-Q_i)*ds**2.0)*x**2.0
     $     +(1+tau)*pr*ds**2.0*x**4.0)

c      dy(1)=-y(1)/x-y(1)*y(1)/x+x*(ifac*(Q-Q_e))
c     $     *(-Q*(Q-Q_i)+ifac*(Q-Q_i)*(pr+c_beta**2.0)*x**2.0
c     $     +pr*c_beta**2.0*x**4.0)/(ifac*(Q-Q_e)
c     $     +(c_beta**2.0+ifac*(Q-Q_i)*ds**2.0)*x**2.0
c     $     +(1+tau)*pr*ds**2.0*x**4.0)

      RETURN
      END SUBROUTINE w_der_temp
c-----------------------------------------------------------------------
c     calculate delta based on direct phi_der (obsolete).
c-----------------------------------------------------------------------
      FUNCTION directred(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
      REAL(r8),INTENT(IN) :: intau
      COMPLEX(r8) :: directred

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,iopt,mf

      REAL(r8) :: xintv,x,xout,rtol,jac,xmax
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,atol,rwork
      
      Q=inQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      pe=inpe
      c_beta=inc_beta
      ds=inds
      tau=intau

      neq = 4
      itol = 2
      rtol = 1e-3
      ALLOCATE(atol(neq),y(2),dy(2))
      atol(:) = 1e-4
      itask = 1
      istate = 1
      iopt = 0
      mf = 10
      liw = 20    
      lrw = 22+16*neq
      ALLOCATE(iwork(liw),rwork(lrw))
      
      xintv = 0.1
      x=0.1
      istep=1
      xout=x+xintv
      xmax=1e3

      y(1)=1e-10
      y(2)=y(1)/(c_beta/(sqrt(1+tau)*ds)*(1+ifac*(Q-Q_e)*x**2.0)*x)

      DO
         istep=istep+1
         CALL lsode(phi_der,neq,y,x,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         xout=xout+xintv
         IF (xmax-xout<xintv) EXIT
      ENDDO
      CALL phi_der(neq,x,y,dy)
      directred=pi/(y(1)-x*dy(1))*dy(1)
      WRITE(*,*)directred
      DEALLOCATE(atol,y,dy,iwork,rwork)

      END FUNCTION directred
c-----------------------------------------------------------------------
c     direct integration (obsolete).
c-----------------------------------------------------------------------
      SUBROUTINE phi_der(neq,x,y,dy)
      
      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: x
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dy
      COMPLEX(r8), PARAMETER :: ifac=(0,1)      

      dy(1)=(1+ifac*(Q-Q_e)*x**2.0)/x**2.0*y(2)
      dy(2)=(-Q*(Q-Q_i)*x**4.0+ifac*(Q-Q_i)*(pr+c_beta**2.0)*x**2.0
     $     +pr*c_beta**2.0)/(ifac*(Q-Q_e)*x**8.0
     $     +(c_beta**2.0+ifac*(Q-Q_i)*ds**2.0)*x**6.0
     $     +(1+tau)*pr*ds**2.0*x**4.0)*y(1)
      RETURN
      END SUBROUTINE phi_der  
      
      END MODULE delta_mod
