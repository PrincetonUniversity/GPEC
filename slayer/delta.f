      MODULE delta_mod

      USE sglobal_mod

      IMPLICIT NONE

      LOGICAL :: riccati_out,parflow_flag,PeOhmOnly_flag

      CONTAINS
c-----------------------------------------------------------------------
c     calculate delta based on riccati w_der formulation.
c-----------------------------------------------------------------------
      FUNCTION riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,
     $     iinQ,inx,iny)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
	  REAL(r8),INTENT(IN) :: intau
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ,inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati

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
      riccati=pi/dy(1)
      DEALLOCATE(atol,y,dy,iwork,rwork)      

      END FUNCTION riccati
c-----------------------------------------------------------------------
c     calculate delta with the ion parallel flow (four-field model)
c     without electron viscosity and thermal conductivity.
c     Subroutines used in this function are w_derl, w_derr,
c     dy_der_yl, dy_der_yr, Transform_R, Update_Delta.
c     Reference : Y. Lee, J.-K. Park, & Y.-S. Na, (2024), Effect of
c         parallel flow on resonant layer responses in high beta plasmas.
c         Nucl. Fusion, 64(10), 106058.
c-----------------------------------------------------------------------
      FUNCTION riccati4(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,
     $     iinQ,inx,iny)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inpe,inc_beta,inds
      REAL(r8),INTENT(IN) :: intau
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ,inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati4, Delta_old, Delta_new

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,lzw,iopt,mf
      INTEGER :: ml,mu,nrpd,ipar

      REAL(r8) :: xintv,x,xout,jac,xmin
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy,zwork
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: pd

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,rwork,rtol,atol
      
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

      neq = 25
      itol = 4
      ALLOCATE(y(neq),dy(neq),pd(neq,neq),rtol(neq),atol(neq))
      rtol = 1e-9
      atol = 1e-9
      itask = 2
      istate = 1
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
      xout=1e-4
      rwork(1)=1e-5

      y = 0.0
      dy = 0.0

      IF (riccati_out) THEN
         istep = 1
         itask = 2
         ! leftside subregion
         DO WHILE (x<xout)
            istep=istep+1
            CALL ZVODE(w_derl,neq,y,x,xout,itol,rtol,atol,itask,
     $      istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dy_der_yl,mf,
     $      ipar)
         ENDDO        

         ! Transform R matrix from left to right region
         CALL Transform_R(y,x,neq)

         istate=1
         rtol=1e-9
         atol=1e-9
         Delta_old=0.0
         Delta_new=1.0

         ! rightside subregion
         OPEN(UNIT=bin_unit,FILE='riccati4_profile.bin',STATUS='UNKNOWN'
     $        ,POSITION='REWIND',FORM='UNFORMATTED')   
         
         OPEN(UNIT=out2_unit,FILE='riccati4_profile.out',
     $        STATUS='UNKNOWN')
         WRITE(out2_unit,'(1x,3(a17))'),"x","RE(delta(x))",
     $   "IM(delta(y))"      
         DO WHILE (abs(Delta_new-Delta_old)/abs(Delta_new)>0.01)
           xout=xout+10.0
           Delta_old=Delta_new
           DO WHILE (x<xout)
              istep=istep+1
              CALL ZVODE(w_derr,neq,y,x,xout,itol,rtol,atol,itask,
     $        istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dy_der_yr,mf,
     $        ipar)
              CALL Update_Delta(Delta_new,y,x,neq)
              WRITE(bin_unit)REAL(x,4),REAL(REAL(Delta_new),4),
     $        REAL(AIMAG(Delta_new),4) 
              WRITE(out2_unit,'(1x,3(es17.8e3))')x,REAL(Delta_new)
     $        ,AIMAG(Delta_new)
           ENDDO        
         END DO
         CLOSE(bin_unit)
         CLOSE(out2_unit)
      ELSE
         istep = 1
         itask = 2
         ! leftside subregion
         DO WHILE (x<xout)
            istep=istep+1
            CALL ZVODE(w_derl,neq,y,x,xout,itol,rtol,atol,itask,
     $      istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dy_der_yl,mf,
     $      ipar)
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
           xout=xout+10.0
           Delta_old=Delta_new
           DO WHILE (x<xout)
              istep=istep+1
              CALL ZVODE(w_derr,neq,y,x,xout,itol,rtol,atol,itask,
     $        istate,iopt,zwork,lzw,rwork,lrw,iwork,liw,dy_der_yr,mf,
     $        ipar)
           ENDDO        
           CALL Update_Delta(Delta_new,y,x,neq)
           WRITE(*,*) xout
         END DO
      ENDIF

      DEALLOCATE(y,dy,pd,rtol,atol,iwork,rwork,zwork)
      riccati4=Delta_new
      END FUNCTION riccati4

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

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in left region
c     pd (output) : analytic jacobian tensor d(dR/dX)/dR (Eqn. 44)
c-----------------------------------------------------------------------
      SUBROUTINE dy_der_yl(neq,x,y,ml,mu,pd,nrpd,ipar)

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
      END SUBROUTINE dy_der_yl

c-----------------------------------------------------------------------
c     Subroutine for riccati integration based on four-field models.
c     x (input)   : stretched variable
c     y (input)   : riccati matrix in right region
c     pd (output) : analytic jacobian tensor d(dR/dX)/dR (Eqn. 44)
c-----------------------------------------------------------------------
      SUBROUTINE dy_der_yr(neq,x,y,ml,mu,pd,nrpd,ipar)

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
      END SUBROUTINE dy_der_yr

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
