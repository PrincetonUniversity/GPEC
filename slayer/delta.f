      MODULE delta_mod
c-----------------------------------------------------------------------
c     delta_mod: Riccati-based tearing-mode layer Delta solvers.
c
c     Contains three Riccati formulations:
c       riccati       - standard formulation (non-stiff, lsode mf=10)
c       riccati_del_s - del_s formulation (stiff, lsode mf=21)
c       riccati_f     - Fitzpatrick P_perp/P_tor formulation (mf=21)
c
c     Each solver integrates a Riccati ODE for W(x) from a large-|x|
c     asymptotic boundary condition inward to x ~ 0, then extracts
c     Delta = pi / W'(0).
c
c     Associated ODE subroutines (w_der, w_der_del_s, w_der_f) and
c     Jacobian subroutines (jac_del_s, jac_f) follow below.
c
c     Module-level flags:
c       riccati_out    - write W(x) profile to binary + text file
c       parflow_flag   - include parallel electron flow terms in w_der
c       PeOhmOnly_flag - retain only Pe-Ohm coupling in w_der
c-----------------------------------------------------------------------

      USE sglobal_mod

      IMPLICIT NONE

c --- module-level control flags
      LOGICAL :: riccati_out    ! write W(x) profile to binary + text
      LOGICAL :: parflow_flag   ! enable parallel-flow terms in w_der
      LOGICAL :: PeOhmOnly_flag ! Pe-Ohm-only coupling in w_der

      CONTAINS
c-----------------------------------------------------------------------
c     riccati: compute Delta via Riccati integration of w_der.
c
c     Integrates W(x) from x_start inward to x_min using lsode
c     (non-stiff, mf=10).  Optionally applies the layfac singularity
c     guard when Q is near Q_e.  Returns Delta = pi / W'(x_min).
c
c     If riccati_out = .TRUE., writes the W(x) profile to
c     slayer_riccati_profile_n<sn>.{bin,out}.
c
c     BUG FLAG 1: xintv and xfac are declared but never used.
c       Remove them.
c-----------------------------------------------------------------------
      FUNCTION riccati(inQ,inQ_e,inQ_i,inpr,inc_beta,inds,intau,inpe,
     $     iinQ,inx,iny)

c --- input arguments: physical parameters for this surface
      REAL(r8),INTENT(IN) :: inQ       ! real part of Q
      REAL(r8),INTENT(IN) :: inQ_e     ! electron diamagnetic freq
      REAL(r8),INTENT(IN) :: inQ_i     ! ion diamagnetic freq
      REAL(r8),INTENT(IN) :: inpr      ! Prandtl number
      REAL(r8),INTENT(IN) :: inpe      ! electron Prandtl number
      REAL(r8),INTENT(IN) :: inc_beta  ! c_beta parameter
      REAL(r8),INTENT(IN) :: inds      ! magnetic shear ds
      REAL(r8),INTENT(IN) :: intau     ! ion-to-electron temp ratio
c --- optional arguments
      REAL(r8),INTENT(IN),OPTIONAL :: iinQ ! imaginary part of Q
      REAL(r8),INTENT(IN),OPTIONAL :: inx  ! override starting x
      COMPLEX(r8),INTENT(IN),OPTIONAL :: iny ! override starting W
c --- function result
      COMPLEX(r8) :: riccati

c --- lsode solver control
      INTEGER :: istep           ! integration step counter
      INTEGER :: neq             ! number of equations (=2)
      INTEGER :: itol,itask      ! lsode tolerance/task flags
      INTEGER :: istate,iopt,mf  ! lsode state/option/method flags
      INTEGER :: liw,lrw         ! lsode work array sizes
      REAL(r8) :: x,xout         ! current and target x
      REAL(r8) :: xmin           ! inner integration bound
      REAL(r8) :: rtol           ! relative tolerance
      REAL(r8) :: jac            ! dummy Jacobian (unused for mf=10)
c --- work arrays
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: y,dy      ! W and dW
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork          ! lsode int work
      REAL(r8), DIMENSION(:), ALLOCATABLE :: atol,rwork    ! lsode real work

c --- copy input arguments to module-level globals for w_der
      Q=inQ
      IF(present(iinQ)) Q=inQ+ifac*iinQ
      Q_e=inQ_e
      Q_i=inQ_i
      pr=inpr
      pe=inpe
      c_beta=inc_beta
      ds=inds
      tau=intau

c --- singularity guard: displace Q away from Q_e when too close
      IF ((layfac>0).AND.(ABS(Q-Q_e)<layfac)) THEN
         Q=Q_e+layfac*EXP(ifac*ATAN2(AIMAG(Q-Q_e),REAL(Q-Q_e)))
      ENDIF

c --- configure lsode: non-stiff Adams method (mf=10)
      neq = 2
      itol = 2
      rtol = 1e-7
      ALLOCATE(atol(neq),y(1),dy(1))
      atol(:) = 1e-7*pr**0.4
      itask = 2
      istate = 1
      iopt = 1              ! enable optional inputs (iwork(6))
      mf = 10               ! non-stiff, no Jacobian needed
      liw = 20
      lrw = 22+16*neq
      ALLOCATE(iwork(liw),rwork(lrw))

c --- set maximum internal steps
      iwork=0
      iwork(6)=10000         ! MXSTEP: max internal steps per call
      rwork=0

c --- boundary condition: asymptotic W at large x
c     BUG FLAG 2: inline comment said "To be updated" -- verify formula.
      x=20.0
      xmin=1e-3
      IF(present(inx)) x=inx
      xout=xmin
      y(1)=-c_beta/sqrt((1+tau))/ds*x**2.0
      IF(present(iny)) y(1)=iny


c --- integrate W(x) from x_start inward to x_min via lsode
      IF (riccati_out) THEN
c        profile output: step-by-step integration with file writes
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
c        single-shot integration to x_min
         istep = 1
         itask = 1
         CALL lsode(w_der,neq,y,x,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDIF

c --- extract Delta from final W derivative at x_min
c     NOTE: W -> 0 when Q -> Q_e (see layfac guard above).
      CALL w_der(neq,x,y,dy)
      riccati=pi/dy(1)
      DEALLOCATE(atol,y,dy,iwork,rwork)

      END FUNCTION riccati
c-----------------------------------------------------------------------
c     riccati_del_s: compute Delta via the del_s Riccati formulation.
c
c     Uses a stiff solver (lsode mf=21) with user-supplied Jacobian
c     (jac_del_s).  Integrates W(q) from large q inward to q_min.
c     Returns Delta = -(pi / sqrt(1+1/tau)) * W'(q_min).
c
c     BUG FLAG 3: arguments inQ, inc_beta, ind_beta, intau are never
c       used to set their module-level counterparts (Q, c_beta, d_beta,
c       tau).  Either add assignments (e.g. tau=intau) or remove the
c       unused arguments if the caller sets them beforehand.
c     BUG FLAG 4: variables y, dy, xfac, xintv, ml, mu, nrpd are
c       declared but never used -- remove them.
c-----------------------------------------------------------------------
      FUNCTION riccati_del_s(inQ_e,inQ_i,inpr,inx,iny)

c --- input arguments
      REAL(r8),INTENT(IN) :: inQ_e     ! electron diamagnetic freq
      REAL(r8),INTENT(IN) :: inQ_i     ! ion diamagnetic freq
      REAL(r8),INTENT(IN) :: inpr      ! mapped to P_perp (see below)
c --- optional arguments
c     BUG FLAG 5: inx is declared OPTIONAL but my_q=inx is accessed
c       unconditionally.  If inx is ever absent, this will crash.
c       Either make inx required or add IF(present(inx)) guard.
      REAL(r8),INTENT(IN) :: inx  ! starting q for integration
      COMPLEX(r8),INTENT(IN),OPTIONAL :: iny ! override starting W
c --- function result
      COMPLEX(r8) :: riccati_del_s

c --- lsode solver control
      INTEGER :: istep           ! integration step counter
      INTEGER :: neq             ! number of equations (=2)
      INTEGER :: itol,itask      ! lsode tolerance/task flags
      INTEGER :: istate,iopt,mf  ! lsode state/option/method flags
      INTEGER :: liw,lrw         ! lsode work array sizes
      REAL(r8) :: x              ! secondary x variable (set from inx)
      REAL(r8) :: xout           ! target integration endpoint
      REAL(r8) :: xmin           ! inner integration bound
      REAL(r8) :: rtol           ! relative tolerance
      REAL(r8) :: jac            ! dummy (overridden by jac_del_s)
      REAL(r8) :: my_q           ! integration variable (large -> small)
      REAL(r8) :: P_hat          ! normalized P_perp
      REAL(r8) :: alpha          ! boundary condition coefficient
c --- work arrays
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: W,dW_dq   ! W and dW/dq
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork          ! lsode int work
      REAL(r8), DIMENSION(:), ALLOCATABLE :: atol,rwork    ! lsode real work

c --- configure lsode: stiff BDF method with user Jacobian (mf=21)
      neq = 2
      itol = 2
      rtol = 1e-10
      ALLOCATE(atol(neq),W(1),dW_dq(1))
      atol(:) = 1e-10
      itask = 2
      istate = 1
      iopt = 1              ! enable optional inputs (iwork(6))
      mf = 21               ! stiff, user-supplied Jacobian (jac_del_s)
      liw = 20*2
      lrw = 22+9*neq+neq**2 ! stiff work array size
      ALLOCATE(iwork(liw+neq),rwork(lrw))

c --- set maximum internal steps
      iwork=0
      iwork(6)=50000         ! MXSTEP: max internal steps per call
      rwork=0

c --- set starting integration point
      my_q=inx               ! start backwards integration at large q
      xmin=1e-5
      IF(present(inx)) x=inx
      xout=xmin

c --- copy input arguments to module-level globals for w_der_del_s
      Q_e = inQ_e
      Q_i = inQ_i
      P_perp = inpr

c --- BUG FLAG 6: P_hat was computed BEFORE P_perp=inpr, so it used
c       the stale module-level P_perp.  Now moved after assignment.
      P_hat = P_perp / D_norm**6.0

c --- asymptotic boundary condition at large q
      alpha = (P_hat/(1+1/tau))**0.5
      W(1) = -alpha*my_q**2 - 0.5
      IF(present(iny)) W(1)=iny

c --- integrate W(q) from q_start inward to q_min via lsode
      IF (riccati_out) THEN
c        profile output: step-by-step integration with file writes
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
     $           itask,istate,iopt,rwork,lrw,iwork,liw,jac_del_s,mf)
            WRITE(bin_unit)REAL(my_q,4),REAL(REAL(W),4),REAL(AIMAG(W),4)
            WRITE(out2_unit,'(1x,3(es17.8e3))')my_q,REAL(W),AIMAG(W)
         ENDDO
         CLOSE(bin_unit)
         CLOSE(out2_unit)
      ELSE
c        single-shot integration to q_min
         istep = 1
         itask = 1
         CALL lsode(w_der_del_s,neq,W,my_q,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac_del_s,mf)
      ENDIF

c --- extract Delta from final W derivative at q_min
      CALL w_der_del_s(neq,my_q,W,dW_dq)
      riccati_del_s=-( pi/((1+1/tau)**0.5) )*dW_dq(1)
      DEALLOCATE(atol,W,dW_dq,iwork,rwork)

      END FUNCTION riccati_del_s
c-----------------------------------------------------------------------
c     jacobian for riccati_del_s(): pd = dF/dW for stiff lsode.
c-----------------------------------------------------------------------
      SUBROUTINE jac_del_s(neq, my_q, W, ml, mu, pd, nrpd)
            INTEGER, INTENT(IN) :: neq, ml, mu, nrpd
            REAL(r8), INTENT(IN) :: my_q
            COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
            COMPLEX(r8), DIMENSION(nrpd,neq), INTENT(INOUT) :: pd
            pd(1,1) = 1.0/my_q - 2.0d0*W(1)/my_q
      END SUBROUTINE jac_del_s
c-----------------------------------------------------------------------
c     w_der_del_s: ODE right-hand side dW/dq for riccati_del_s.
c     Implements the del_s dispersion relation using normalised
c     quantities Q_hat, P_perp_hat, P_tor_hat.
c-----------------------------------------------------------------------
      SUBROUTINE w_der_del_s(neq,my_q,W,dW_dq)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: my_q
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dW_dq
      REAL(r8) :: Q_hat, P_tor_hat, P_perp_hat
      COMPLEX(r8) :: E,F

c --- normalise physical quantities
      Q_hat = (Q_e*(1+tau)/tau) / D_norm**4.0
      P_perp_hat = P_perp / D_norm**6.0
c     BUG FLAG 7: P_tor_hat is assigned from P_perp, not P_tor.
c       Benchmark values differ (P_perp_hat=0.377, P_tor_hat=1.15),
c       suggesting this should be: P_tor_hat = P_tor / D_norm**6.0
      P_tor_hat = P_perp / D_norm**6.0
c --- build the E and F dispersion coefficients
      E = (-(Q_hat**2)/(1+1/tau)) - ifac*Q_hat*(P_perp_hat+
     $  P_tor_hat)*(my_q**2) + P_perp_hat*P_tor_hat*(my_q**4)
      F = P_perp_hat - ifac*Q_hat + (1+1/tau)*P_tor_hat*my_q**2

c --- Riccati ODE: dW/dq = W/q - W^2/q + q*E/F
      dW_dq(1)=W(1)/my_q - (W(1)**2)/my_q + (my_q*E)/F
      RETURN
      END SUBROUTINE w_der_del_s
c-----------------------------------------------------------------------
c     riccati_f: compute Delta via Fitzpatrick P_perp / P_tor
c     Riccati formulation.
c
c     Uses a stiff solver (lsode mf=21) with user-supplied Jacobian
c     (jac_f).  Boundary conditions are set analytically from the
c     large-p asymptotic behaviour; the branch depends on whether
c     D_norm^2 exceeds iota_e * P_perp / P_tor^(2/3).
c
c     BUG FLAG 8: argument tmp_g (COMPLEX) is never used in the
c       function body.  All references use the module-level variable
c       g_tmp instead.  The caller (slayer.f:785) passes g_tmp as
c       tmp_g, so in practice the values match -- but tmp_g is
c       redundant.  Likely needs: g_tmp = tmp_g  at the top, or
c       remove the argument and rely on the module variable.
c     BUG FLAG 9: variables xintv, xfac, y, dy, ck_1, ck_2, ml, mu,
c       nrpd, alpha are declared but never used -- remove them.
c       Optional argument inx is also never referenced in the body.
c-----------------------------------------------------------------------
      FUNCTION riccati_f(tmp_g,inx)

c --- input arguments
      COMPLEX(r8),INTENT(IN) :: tmp_g   ! growth rate (UNUSED -- BUG FLAG 8)
      REAL(r8),INTENT(IN),OPTIONAL :: inx ! (UNUSED -- BUG FLAG 9)
c --- function result
      COMPLEX(r8) :: riccati_f

c --- lsode solver control
      INTEGER :: istep           ! integration step counter
      INTEGER :: neq             ! number of equations (=2)
      INTEGER :: itol,itask      ! lsode tolerance/task flags
      INTEGER :: istate,iopt,mf  ! lsode state/option/method flags
      INTEGER :: liw,lrw         ! lsode work array sizes
      REAL(r8) :: xout           ! target integration endpoint
      REAL(r8) :: xmin           ! inner integration bound
      REAL(r8) :: rtol           ! relative tolerance
      REAL(r8) :: jac            ! dummy (overridden by jac_f)
      REAL(r8) :: my_p           ! integration variable p (large -> small)
      REAL(r8) :: bk             ! asymptotic coefficient b_k
c --- boundary-condition intermediates
      COMPLEX(r8) :: ak          ! asymptotic coefficient a_k
      COMPLEX(r8) :: ck          ! asymptotic coefficient c_k
      COMPLEX(r8) :: xk          ! asymptotic coefficient x_k
      COMPLEX(r8) :: W_bound     ! boundary value for W(p_start)
c --- work arrays
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: W,dWdp    ! W and dW/dp
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork          ! lsode int work
      REAL(r8), DIMENSION(:), ALLOCATABLE :: atol,rwork    ! lsode real work

c --- configure lsode: stiff BDF method with user Jacobian (mf=21)
      neq = 2
      itol = 2
      rtol = 1e-10
      ALLOCATE(atol(neq),W(1),dWdp(1))
      atol(:) = 1e-10
      itask = 2
      istate = 1
      iopt = 1              ! enable optional inputs (iwork(6))
      mf = 21               ! stiff, user-supplied Jacobian (jac_f)
      liw = 20*2
      lrw = 22+9*neq+neq**2 ! stiff work array size
      ALLOCATE(iwork(liw+neq),rwork(lrw))

c --- set maximum internal steps
      iwork=0
      iwork(6)=50000         ! MXSTEP: max internal steps per call
      rwork=0

      xmin=1e-6
      xout=xmin

c --- compute starting p and W boundary condition
c     Branch on asymptotic regime: D_norm^2 vs iota_e*P_perp/P_tor^(2/3)
c --- branch 1: D_norm^2 > iota_e * P_perp / P_tor^(2/3)
c     large-D_norm regime: p scales with (P_tor*D_norm^2/(iota_e*...))
      IF ((D_norm**2.0) > ((iota_e*P_perp)/(P_tor**(2.0/3.0)))) THEN
          my_p = ( (P_tor*D_norm**2)/(iota_e*P_tor*P_perp) )**0.25
          IF (my_p < 6.0) THEN
            my_p = 6.0
          END IF

          ak = -(g_tmp + ifac*Q_e)
          bk = (iota_e*P_perp*P_tor)/(P_tor*(D_norm**2.0))

          ck = bk*(1+(g_tmp+ifac*Q_i)*((P_tor+P_perp)/(P_tor*P_perp))-
     $       (P_perp+(g_tmp + 
     $       ifac*Q_i)*(D_norm**2.0) )*(iota_e/(P_tor*(D_norm**2.0))))

          xk = (ck - SQRT(bk)*(1 - SQRT(bk)*ak)) / (2.0*SQRT(bk))

          W_bound = xk - SQRT(bk)*my_p
      ELSE
c --- branch 2: D_norm^2 <= iota_e * P_perp / P_tor^(2/3)
c     small-D_norm regime: p scales with 1/P_tor^(1/6)
          my_p = 1.0/(P_tor**(1.0/6.0))
          IF (my_p < 6.0) THEN
            my_p = 6.0
          END IF

          ak = -(g_tmp + ifac*Q_e)
          bk = P_tor
          ck = -ifac*(Q_e - Q_i)*(P_tor/P_perp) + (g_tmp + ifac*Q_i)
          xk = (ak*bk - ck)/(2.0*SQRT(bk))

          W_bound = -1 + xk*my_p - SQRT(bk)*(my_p**3.0)
      END IF

      W(1) = W_bound

c --- integrate W(p) from p_start inward to p_min via lsode
      IF (riccati_out) THEN
c        profile output: step-by-step integration with file writes
         istep = 1
         itask = 2
         OPEN(UNIT=bin_unit,FILE='slayer_riccati_profile_n'//
     $      TRIM(sn)//'.bin',STATUS='UNKNOWN',
     $      POSITION='REWIND',FORM='UNFORMATTED')

         OPEN(UNIT=out2_unit,FILE='slayer_riccati_profile_n'//
     $      TRIM(sn)//'.out',STATUS='UNKNOWN')
         WRITE(out2_unit,'(1x,3(a17))'),"x","RE(y)","IM(y)"
         DO WHILE (my_p>xout)
            istep=istep+1
            CALL lsode(w_der_f,neq,W,my_p,xout,itol,rtol,atol,
     $           itask,istate,iopt,rwork,lrw,iwork,liw,jac_f,mf)
            WRITE(bin_unit)REAL(my_p,4),REAL(REAL(W),4),
     $                                    REAL(AIMAG(W),4)
            WRITE(out2_unit,'(1x,3(es17.8e3))')my_p,REAL(W),AIMAG(W)
         ENDDO
         CLOSE(bin_unit)
         CLOSE(out2_unit)
      ELSE
c        single-shot integration to p_min
         istep = 1
         itask = 1
         CALL lsode(w_der_f,neq,W,my_p,xout,itol,rtol,atol,
     $        itask,istate,iopt,rwork,lrw,iwork,liw,jac_f,mf)
      ENDIF

c --- extract Delta from final W derivative at p_min
      CALL w_der_f(neq,my_p,W,dWdp)
      riccati_f = pi / dWdp(1)
      DEALLOCATE(atol,W,dWdp,iwork,rwork)

      END FUNCTION riccati_f
c-----------------------------------------------------------------------
c     jacobian for riccati_f(): pd = dF/dW for stiff lsode.
c-----------------------------------------------------------------------
      SUBROUTINE jac_f(neq, my_p, W, ml, mu, pd, nrpd)
            INTEGER, INTENT(IN) :: neq, ml, mu, nrpd
            REAL(r8), INTENT(IN) :: my_p
            COMPLEX(r8) :: fA_p
            COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
            COMPLEX(r8), DIMENSION(nrpd,neq), INTENT(INOUT) :: pd

            fA_p = (g_tmp + ifac*Q_e - (my_p**2)) / (g_tmp + 
     $          ifac*Q_e + (my_p**2.0))

            pd(1,1) = (-fA_p/my_p) - (2.0*W(1))/my_p
      END SUBROUTINE jac_f
c-----------------------------------------------------------------------
c     w_der_f: ODE right-hand side dW/dp for riccati_f.
c     Implements the Fitzpatrick P_perp / P_tor dispersion relation.
c     Coefficients fA, fB, fC are evaluated at the current p.
c-----------------------------------------------------------------------
      SUBROUTINE w_der_f(neq,my_p,W,dWdp)

      INTEGER, INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: my_p
      COMPLEX(r8), DIMENSION(neq), INTENT(IN) :: W
      COMPLEX(r8), DIMENSION(neq), INTENT(OUT) :: dWdp
      COMPLEX(r8) :: fA, fA_prime, fB, fC
      
      ! Evaluate coefficients at the current p
      fA = (my_p**2)/(g_tmp + ifac*Q_e + (my_p**2.0))
      fA_prime = (g_tmp + ifac*Q_e - (my_p**2)) / (g_tmp + 
     $          ifac*Q_e + (my_p**2.0))
      fB = g_tmp*(g_tmp + ifac*Q_i) + (g_tmp + 
     $    ifac*Q_i)*(P_perp+P_tor)*(my_p**2.0) + 
     $    (P_perp*P_tor)*(my_p**4.0)
      fC = g_tmp + ifac*Q_e + ( P_perp + (g_tmp + 
     $    ifac*Q_i)*(D_norm**2.0))*(my_p**2.0) + 
     $    (1.0/iota_e)*P_tor*(D_norm**2.0)*(my_p**4.0)

      dWdp(1) = -(fA_prime/my_p)*W(1) - (W(1)**2.0)/my_p + 
     $          (fB/(fA*fC))*(my_p**3.0)

      RETURN
      END SUBROUTINE w_der_f
c-----------------------------------------------------------------------
c     W derivative for riccati()
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
      !COMPLEX(r8), PARAMETER :: ifac=(0,1)

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
      !COMPLEX(r8), PARAMETER :: ifac=(0,1)

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
      !COMPLEX(r8), PARAMETER :: ifac=(0,1)

      dy(1)=(1+ifac*(Q-Q_e)*x**2.0)/x**2.0*y(2)
      dy(2)=(-Q*(Q-Q_i)*x**4.0+ifac*(Q-Q_i)*(pr+c_beta**2.0)*x**2.0
     $     +pr*c_beta**2.0)/(ifac*(Q-Q_e)*x**8.0
     $     +(c_beta**2.0+ifac*(Q-Q_i)*ds**2.0)*x**6.0
     $     +(1+tau)*pr*ds**2.0*x**4.0)*y(1)
      RETURN
      END SUBROUTINE phi_der

      END MODULE delta_mod