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
         WRITE(out2_unit,'(1x,3(a17))'),"x","RE(y)","IM(y)"
         DO WHILE (x>xout)
            istep=istep+1
            CALL lsode(w_der,neq,y,x,xout,itol,rtol,atol,
     $           itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
            WRITE(bin_unit)REAL(x,4),REAL(REAL(y),4),REAL(AIMAG(y),4)
            WRITE(out2_unit,'(1x,3(es17.8e3))')x,REAL(y),AIMAG(y)
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
c     calculate delta based on riccati w_der formulation.
c-----------------------------------------------------------------------
      FUNCTION riccati_del_s(inQ,inQ_e,inQ_i,inpr,inc_beta,ind_beta,
     $     intau,inx,iny)

      REAL(r8),INTENT(IN) :: inQ,inQ_e,inQ_i,inpr,inc_beta,ind_beta
	  REAL(r8),INTENT(IN) :: intau
      REAL(r8),INTENT(IN),OPTIONAL :: inx
      COMPLEX(r8), INTENT(IN), OPTIONAL :: iny
      COMPLEX(r8) :: riccati_del_s

      INTEGER :: istep,neq,itol,itask,istate,liw,lrw,iopt,mf
      INTEGER :: ml = 0, mu = 0, nrpd = 1

      REAL(r8) :: xintv,x,xout,rtol,jac,xmin,my_q,P_hat,alpha
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: W,dW_dq,y,dy

      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xfac,atol,rwork

      Q=inQ
      !IF(present(iinQ)) Q=inQ+ifac*iinQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      !pe=inpe
      c_beta=inc_beta
      !ds=inds
      tau=intau

      IF ((layfac>0).AND.(ABS(Q-Q_e)<layfac)) THEN
         Q=Q_e+layfac*EXP(ifac*ATAN2(AIMAG(Q-Q_e),REAL(Q-Q_e)))
      ENDIF

      neq = 2
      itol = 2
      rtol = 1e-08    !changed to 1e-08       !1e-7*pr**0.4 ! !1e-7 at front 1e-6 !e-4
      ALLOCATE(atol(neq),W(1),dW_dq(1))
      atol(:) = 1e-08!*pr**0.4 ! changed to 1e-08 1e-8 !e-4
      itask = 2
      istate = 1
      iopt = 0
      mf = 21!21 IS STIFF WITH USER-SPECIFIED JACOBIAN, 10 iS NON STIFF
      liw = 20*2
      lrw = 22+9*neq+neq**2 !just (22+16*neq) for mf=10
      ALLOCATE(iwork(liw+neq),rwork(lrw)) ! just iwork(liw) for mf=10

!     MXSTEP?
      iopt = 1
      iwork=0
      iwork(6)=50000 !5000 ! maximum # of steps per call, e.g. 50000
      rwork=0
!      x=10.0*(1.0+log10(Q/pr))

      !!!!!!!!
      my_q=10.0 ! "starting backwards integration at large q"
      !!!!!!!!

      xmin=1e-5
      IF(present(inx)) x=inx
      xout=xmin

      !y(1)=-c_beta/sqrt((1+tau))/ds*x**2.0 ! it was (1+tau*ds). To be updated.

      P_hat = inpr / ind_beta**6 ! P_perp, 0.377 for Pperp_hat benchmark

      WRITE(*,*)"riccati_del_s inpr = ",inpr
      WRITE(*,*)"riccati_del_s inQ = ",inQ

      alpha = (P_hat/(1+1/tau))**0.5 ! this is actually tau', we need tau
      W(1) = -alpha*my_q**2 - 0.5

      IF(present(iny)) W(1)=iny
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
         WRITE(out2_unit,'(1x,3(a17))'),"x","RE(y)","IM(y)"
         DO WHILE (my_q>xout)
            istep=istep+1
            CALL lsode(w_der_del_s,neq,W,my_q,xout,itol,rtol,atol,
     $           itask,istate,iopt,rwork,lrw,iwork,liw,my_jac,mf)
            WRITE(bin_unit)REAL(my_q,4),REAL(REAL(W),4),REAL(AIMAG(W),4)
            WRITE(out2_unit,'(1x,3(es17.8e3))')my_q,REAL(W),AIMAG(W)
         ENDDO
         CLOSE(bin_unit)
         CLOSE(out2_unit)
      ELSE
         istep = 1
         itask = 1
         CALL lsode(w_der_del_s,neq,W,my_q,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,my_jac,mf)

      ENDIF

      ! w=0 when Q=Q_e. Why?

      CALL w_der_del_s(neq,my_q,W,dW_dq)

      riccati_del_s=-( pi/((1+1/tau)**0.5) )*dW_dq(1)
      DEALLOCATE(atol,W,dW_dq,iwork,rwork)

      END FUNCTION riccati_del_s
c-----------------------------------------------------------------------
c     jacobian for riccati_del_s()
c------------------------------------------- ----------------------------
      SUBROUTINE my_jac(neq, my_q, W, ml, mu, pd, nrpd)
            INTEGER, INTENT(IN) :: neq, ml, mu, nrpd
            REAL(r8), INTENT(IN) :: my_q
            COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
            COMPLEX(r8), DIMENSION(nrpd,neq), INTENT(INOUT) :: pd
            pd(1,1) = 1.0/my_q - 2.0d0*W(1)/my_q
      END SUBROUTINE my_jac
c-----------------------------------------------------------------------
c     riccati integration.
c-----------------------------------------------------------------------
      SUBROUTINE w_der_del_s(neq,my_q,W,dW_dq)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: my_q
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dW_dq
      REAL(r8) :: Q_hat, P_tor_hat, P_perp_hat
      COMPLEX(r8) :: E,F

      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      !Q_hat = Q / ds**4
      Q_hat = (Q_e*(1+tau)/tau) / D_beta_norm**4 ! Q_star = Q_e * (1+tau), 2.4e-02 for benchmark
      P_perp_hat = pr / D_beta_norm**6 ! 0.377 for benchmark
      P_tor_hat = pr / D_beta_norm**6 ! 1.15 for benchmark

      E = (-(Q_hat**2)/(1+1/tau)) - ifac*Q_hat*(P_perp_hat+
     $  P_tor_hat)*(my_q**2) + P_perp_hat*P_tor_hat*(my_q**4) ! P_tor = P_perp
      F = P_perp_hat - ifac*Q_hat + (1+1/tau)*P_tor_hat*my_q**2

      !dy(1)=(-A1 + 1/x)*y(1) - y(1)*y(1)/x - A2*x
      dW_dq(1)=W(1)/my_q - (W(1)**2)/my_q + (my_q*E)/F !p*D = my_q
      RETURN
      END SUBROUTINE w_der_del_s
c
c
c
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