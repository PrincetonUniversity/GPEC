c-----------------------------------------------------------------------
c     file ode.f.
c     sets up and integrates Euler-Lagrange differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ode_mod.
c     1. ode_run.
c     2. ode_set_intervals.
c     3. ode_set_interval_details.
c     4. ode_propagate_LR.
c     5. ode_calc_modes_for_wp.
c     6. ode_set_delta_prime_edge.
c     7. ode_calc_delta_prime.
c     8. ode_fixup.
c     9. ode_power_spacing.
c     10. ode_pseudo_inv.
c     11. ode_nojac.
c     12. ode_itime.
c     13. ode_etime.
c-----------------------------------------------------------------------
c     subprogram 0. ode_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ode_mod
      USE debug_mod
      USE sing_mod
      USE zvode1_mod
      USE sparse_mod
      USE riccati_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER, PRIVATE :: eps=1e-10

      REAL(r8) :: singfac_min=1e-5
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: u
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: uFM_all
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: uFM_sing_inv
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: uFM_sing_init
      INTEGER :: nIntervals
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: psiInters
      INTEGER, DIMENSION(:), ALLOCATABLE :: psiDirs
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: psiPoints

      !Grid-packing variables.
      CHARACTER(20) :: grid_packing
      INTEGER :: nIntervalsTot=33
      !Variables used in dcon.f.
      INTEGER :: ising,nzero
      REAL(r8) :: singfac_max=1e-4
      !Delta_prime variable
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: delta_prime_mat

      LOGICAL :: asymp_at_sing
      LOGICAL :: integrate_riccati
      LOGICAL :: calc_delta_prime
      LOGICAL :: solve_delta_prime_with_sparse_mat
      LOGICAL :: kill_big_soln_for_ideal_dW
      LOGICAL :: calc_dp_with_vac

      TYPE(sing_calculator), DIMENSION(:), TARGET, ALLOCATABLE :: scalc
      REAL(r8) :: axisPsi, outerPsi
      INTEGER :: edgeMidPt
      REAL(r8) :: axis_mid_pt_skew,big_soln_err_tol
      REAL(r8) :: ric_dt,ric_tol

      COMPLEX(r8),DIMENSION(:,:),ALLOCATABLE :: uAxis,uAxisD
      COMPLEX(r8),DIMENSION(:,:),ALLOCATABLE :: uEdge

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ode_run.
c     integrate main ODE's using Fundamental Matrix of Solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_run
         INTEGER :: istate,iopt,itask,itol,mf
         INTEGER :: liw, lrw, lzw
         REAL(r8) :: startPsi, endPsi, psi0, psi1
         REAL(r8) :: t0, t1
         REAL(r8) :: projectionFac
         !
         INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
         REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: zwork
         REAL(r8), DIMENSION(:,:), ALLOCATABLE :: atol

         INTEGER :: iInterval
         INTEGER :: sTime0, sTime, fTime, cr

         REAL(r8) :: rtol(1), atol0, maxatol
         INTEGER :: i, ipert, iS, info
         INTEGER :: neq !this (and ode_nojac) are zvode1 inputs

         COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: uFM,duFM
         COMPLEX(r8), DIMENSION(mpert,mpert) :: X0, eV
         REAL(r8), DIMENSION(mpert) :: eD
         COMPLEX(r8), DIMENSION(1) :: rpar,iparC
         INTEGER, DIMENSION(1) :: ipar, zIdx
         COMPLEX(r8), DIMENSION(2*mpert-1) :: ework
         REAL(r8), DIMENSION(3*mpert-2) :: erwork
         INTEGER :: elwork, einfo, ej

         ! variables used in initial qlow finder
         INTEGER :: it,itmax=50
         INTEGER, DIMENSION(1) :: jpsi
         REAL(r8) :: dpsi,q,q1,eps=1e-10

         !Variables related to asymptotic expansions at sing surfs
         INTEGER :: ipert0
         COMPLEX(r8), DIMENSION(mpert,2*mpert,2) :: ua
         !To invert non-Id fundamental matrices
         INTEGER, DIMENSION(2*mpert) :: ipiv
         COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: uFMInv
         COMPLEX(r8), DIMENSION(2*mpert) :: uwork
         COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: identityMat

!NOTE #1: INTEGER, DIMENSION(:), POINTER, PRIVATE :: iwork
!         REAL(r8), DIMENSION(:), POINTER, PRIVATE :: rwork
!         COMPLEX(r8), DIMENSION(:), POINTER, PRIVATE :: zwork
!         ...are ALLOCATABLE, and create difficulty at the module level.
!NOTE #2: A common block name that appears in a copyin clause must be
!         declared to be a Fortran common block in the same scoping
!         unit in which the copyin clause appears.
!NOTE #3: If a THREADPRIVATE directive specifying a common block name
!         appears in one program unit, then such a directive must also
!         appear in every other program unit that contains a COMMON
!         statement specifying the same name.  It must appear after
!         the last such COMMON statement in the PROGRAM unit.
         DOUBLE PRECISION DBL_ODE11, DBL_ODE21
         INTEGER INT_ODE11, INT_ODE21
         COMMON /ZVOD011/ DBL_ODE11(50), INT_ODE11(33)
         COMMON /ZVOD021/ DBL_ODE21(9), INT_ODE21(8)
!$OMP THREADPRIVATE(/ZVOD011/,/ZVOD021/)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10      FORMAT('Int# ',I4,': [',E12.5,',',E12.5,'] ',F6.3,'s, ',
     $        I4,'steps')
c-----------------------------------------------------------------------
c     start timer and initialize.
c-----------------------------------------------------------------------
         CALL SYSTEM_CLOCK(COUNT=sTime0)
         CALL SYSTEM_CLOCK(COUNT_RATE=cr)
         ALLOCATE(scalc(msing))
c-----------------------------------------------------------------------
c     set global integration interval parameters.
c-----------------------------------------------------------------------
         axisPsi = sq%xs(0)
         ! use newton iteration to find starting psi if qlow it is above q0
         IF(qlow > sq%fs(0, 4))THEN
            jpsi=MINLOC(ABS(sq%fs(:,4)-qlow))
            IF (jpsi(1)>= mpsi) jpsi(1)=mpsi-1
            axisPsi=sq%xs(jpsi(1))
            it=0
            DO
               it=it+1
               CALL spline_eval(sq,axisPsi,1)
               q=sq%f(4)
               q1=sq%f1(4)
               dpsi=(qlow-q)/q1
               axisPsi=axisPsi+dpsi
               IF(ABS(dpsi) < eps*ABS(axisPsi) .OR. it > itmax)EXIT
            ENDDO
         ENDIF
         outerPsi = psilim*(1-eps)
         DO iS = 1,msing
            scalc(iS)%singEdgesLR(1) = sing(iS)%psifac - singfac_min/
     $           ABS(nn*sing(iS)%q1)
            scalc(iS)%singEdgesLR(2) = sing(iS)%psifac + singfac_min/
     $           ABS(nn*sing(iS)%q1)
            print *,iS,": ",scalc(iS)%singEdgesLR(1)," ",
     $           scalc(iS)%singEdgesLR(2)

            !This finds the index of the singular column
            scalc(iS)%sing_col = NINT(nn*sing(iS)%q)-mlow+1
         ENDDO
c-----------------------------------------------------------------------
c     divy up subintervals. record psi and surface type of each.
c-----------------------------------------------------------------------
         CALL ode_set_intervals() !Output: psiInters, psiPoints
c-----------------------------------------------------------------------
c     set integration directions and intvl indices of singular surfaces.
c-----------------------------------------------------------------------
         CALL ode_set_interval_details() !Output: psiDirs, scalc%...
c-----------------------------------------------------------------------
c     initialize arrays for integration.
c-----------------------------------------------------------------------
         ALLOCATE(uFM_all(2*mpert,2*mpert,nIntervals),
     $        uFM_sing_inv(2*mpert,2*mpert,msing,2),
     $        uFM_sing_init(2*mpert,2*mpert,msing,2))
         ALLOCATE(u(mpert,mpert,2))

         uFM_all = 0.0_r8
         uFM_sing_inv = 0.0_r8
         uFM_sing_init = 0.0_r8
         u = 0.0_r8

         neq = (2*mpert)*(2*mpert)
         liw = 30
         lrw = 20+neq
         lzw = 15*neq
         ALLOCATE(rwork(lrw), zwork(lzw), iwork(liw),
     $        atol(2*mpert,2*mpert))

         identityMat = 0.0_r8
         DO i = 1,2*mpert
            identityMat(i,i) = 1.0_r8
         ENDDO
c-----------------------------------------------------------------------
c     proceed to integration.
c-----------------------------------------------------------------------
         IF (integrate_riccati) THEN
            print *,"Using Riccati for non-parallel integration..."
            CALL OMP_SET_NUM_THREADS(1) !Riccati must not be run in ||.
         ELSE
c-----------------------------------------------------------------------
c     enter parallel region.
c-----------------------------------------------------------------------
            !Note: nThreads faster than nThreads-1, despite creating
            !      with vac_parallel more threads than processors.
            IF (verbose_performance_output) THEN
               print *,"ode_run: #Threads=",nThreads
            ENDIF
            CALL OMP_SET_NUM_THREADS(nThreads)
         ENDIF
!$OMP PARALLEL DEFAULT(NONE)
!.......................................................................
!$OMP& SHARED(nIntervals,psiPoints,psiInters,psiDirs,sing,axisPsi,
!$OMP& tol_nr,neq,lzw,lrw,liw,uFM_all,sTime0,cr,nn,
!$OMP& asymp_at_sing,calc_delta_prime,integrate_riccati,
!$OMP& riccati_match_hamiltonian_evals,projectionFac,
!$OMP& uFM_sing_inv,uFM_sing_init,sing_order,outerPsi,psio,scalc,mlow,
!$OMP& fmats,gmats,kmats,mpert,mband,identityMat,grid_packing,sq,
!$OMP& u,X0,ric_dt,ric_tol,psi0,psi1,
!$OMP& zIdx,eV,eD,ework,elwork,erwork,einfo,ej)
!.......................................................................
!$OMP& FIRSTPRIVATE(uFM,ipert,startPsi,endPsi,iopt,itask,itol,mf,
!$OMP& istate,rtol,atol0,atol,maxatol,rpar,ipar,iparC,iwork,rwork,zwork)
!.......................................................................
!$OMP& PRIVATE(iInterval,i,t0,t1,sTime,fTime,ipert0,ua,ising,
!$OMP& ipiv,uFMInv,info,uwork,duFM)
!......FROM zvode1.f....................................................
!$OMP& COPYIN(/ZVOD011/,/ZVOD021/,DBL_ZVODE11,DBL_ZVODE21,
!$OMP& INT_ZVODE11,INT_ZVODE21)
c-----------------------------------------------------------------------
c     allocate arrays in parallel for use w/o allocatable subcomponents.
c-----------------------------------------------------------------------
         sq_s_ix = sq%ix
         IF (.NOT. ALLOCATED(sq_s_f)) THEN
            ALLOCATE(sq_s_f(SIZE(sq%f)))
         ENDIF
         sq_s_f = sq%f !Allocated in sing.f

         fmats_s_ix = fmats%ix
         IF (.NOT. ALLOCATED(fmats_s_f)) THEN
            ALLOCATE(fmats_s_f(SIZE(fmats%f)))
         ENDIF
         fmats_s_f=fmats%f !Allocated in sing.f

         gmats_s_ix = gmats%ix
         ALLOCATE(gmats_s_f(SIZE(gmats%f)))
         gmats_s_f=gmats%f

         kmats_s_ix = kmats%ix
         ALLOCATE(kmats_s_f(SIZE(kmats%f)))
         kmats_s_f=kmats%f
c-----------------------------------------------------------------------
c     OMP DO LOOP: execute parallel loop.
c-----------------------------------------------------------------------
!$OMP DO SCHEDULE(GUIDED)
         DO iInterval = 1,nIntervals
c-----------------------------------------------------------------------
c     set interval start and end "times" (radial coordinate psi).
c-----------------------------------------------------------------------
            IF (psiDirs(iInterval) == 1) THEN
               t0 = psiPoints(iInterval, 1)
               t1 = psiPoints(iInterval, 2)
            ELSEIF (psiDirs(iInterval) == -1) THEN
               t0 = psiPoints(iInterval, 2)
               t1 = psiPoints(iInterval, 1)
            ELSE
               CALL program_stop("Unknown integration direction.")
            ENDIF
            IF (t0 > t1 .AND. psiDirs(iInterval)/=-1) THEN
               CALL program_stop("Integration direction not reversed.")
            ELSEIF (t0 < t1 .AND. psiDirs(iInterval)/=1) THEN
               CALL program_stop("Integration direction not forward.")
            ENDIF
c-----------------------------------------------------------------------
c     prep interval details.
c-----------------------------------------------------------------------
            startPsi = t0
            endPsi = t1
            ising = psiInters(iInterval,3)
            IF (.NOT. ising == 0) THEN
               ipert0 = NINT(nn*sing(ising)%q)-mlow+1
            ENDIF
c-----------------------------------------------------------------------
c     approach #1: riccati. (serial riccati integration)
c-----------------------------------------------------------------------
            IF (integrate_riccati) THEN
               IF (iInterval == 1) THEN
                  X0 = 0.0_r8
                  CALL riccati_controller(startPsi,endPsi,X0,ric_dt,
     $                 ric_tol,"inv_ric",.TRUE.,.FALSE.)
                  CALL riccati_inv(X0)
               ELSE
                  IF (psiInters(iInterval,1) == 2) THEN
                     !This is a singular layer crossing.
                     !Perform eigendecomposition of X0.
                     eV = (X0 + CONJG(TRANSPOSE(X0))) / 2.0_r8
                     elwork = 2*mpert-1
                     CALL zheev('V','U',mpert,eV,mpert,eD,ework,elwork,
     $                    erwork,einfo)
                     !Locate 0-evec, substract 2x its projection.
                     zIdx = MINLOC(ABS(eD))
                     IF (riccati_match_hamiltonian_evals) THEN
                        projectionFac = 1.0_r8
                     ELSE
                        projectionFac = 2.0_r8
                     ENDIF
                     DO ej = 1,mpert
                        !Recall: DOT_PRODUCT takes Hermitian dagger
                        !of first argument for complex arguments.
                        X0(:,ej) = X0(:,ej) - projectionFac *
     $                       DOT_PRODUCT(eV(:,zIdx(1)),
     $                       X0(:,ej))*eV(:,zIdx(1))
                     ENDDO
                     X0 = (X0 + CONJG(TRANSPOSE(X0))) / 2.0_r8
                  ELSE
                     psi0 = MIN(startPsi,endPsi)
                     psi1 = MAX(startPsi,endPsi)
                     IF (iInterval == nIntervals) THEN
                        CALL riccati_controller(psi0,psi1,X0,ric_dt,
     $                       ric_tol,"int_ric",.FALSE.,.TRUE.)
                     ELSE
                        CALL riccati_controller(psi0,psi1,X0,ric_dt,
     $                       ric_tol,"int_ric",.FALSE.,.FALSE.)
                     ENDIF
                  ENDIF
               ENDIF
               !Overwrite u on each step, though only last matters.
               u(:,:,1) = 0.0_r8
               DO i = 1,mpert
                  u(i,i,1) = 1.0_r8
               ENDDO
               u(:,:,2) = X0
            ELSE
c-----------------------------------------------------------------------
c     approach #2: hamiltonian. (parallel hamiltonian integration)
c-----------------------------------------------------------------------
               !Initialize the state transition matrix as the Identity,
               !then modify it where necessary.
               uFM = identityMat
               IF (psiInters(iInterval,1) == 2) THEN
c-----------------------------------------------------------------------
c     the inner layer itself (crossing singularity): initialize.
c-----------------------------------------------------------------------
                  !For the inner layer itself, just omit integration
                  !of the big and small solutions. (Must be set to Id
                  !in those columns at the end.)  The projection of
                  !each mode on the resonant modes is zero'd.
                  uFM(ipert0,:) = 0
                  uFM(:,ipert0) = 0
                  uFM(ipert0+mpert,:) = 0
                  uFM(:,ipert0+mpert) = 0
c-----------------------------------------------------------------------
c     the inner layer itself (crossing singularity): integrate.
c-----------------------------------------------------------------------
                  CALL sing_derFM(neq,startPsi,uFM,duFM,rpar(1),
     $                 iparC(1))
                  uFM = uFM + duFM * (endPsi-startPsi)

                  !For the resonant modes, set to Id in the end,
                  !and remove all projections on these modes.
                  uFM(:,ipert0) = 0
                  uFM(:,ipert0+mpert) = 0
                  uFM(ipert0,:) = 0
                  uFM(ipert0+mpert,:) = 0
                  uFM(ipert0,ipert0) = 1
                  uFM(ipert0+mpert,ipert0+mpert) = 1
               ELSE
c-----------------------------------------------------------------------
c     now for other intervals, set integration initial conditions.
c-----------------------------------------------------------------------
                  IF(asymp_at_sing .AND. psiInters(iInterval,3)/=0)THEN
                     !Perform asymptotic expansion at singular surface
                     !Get q (big) and (p) small asymp expansions at t0
                     !for resonant AND nonresonant modes.
                     CALL sing_get_ua(ising,t0,ua)
                     DO i = 1,mpert
                        uFM(i,:) = ua(i,:,1)
                        uFM(i+mpert,:) = ua(i,:,2)
                     ENDDO

                     !Invert the init. fund. matrix (/=Id), to save it.
                     uFMInv = uFM
                     CALL ZGETRF(2*mpert,2*mpert,uFMInv,2*mpert,ipiv,
     $                    info)
                     CALL ZGETRI(2*mpert,uFMInv,2*mpert,ipiv,uwork,
     $                    2*mpert,info)
c-----------------------------------------------------------------------
c     save the interval's initial fundamental matrix and its inverse.
c-----------------------------------------------------------------------
                     IF (psiInters(iInterval,2) == 1) THEN
                     !Right edge of the the interval is a sing surface
                     !so it modifies the LEFT side of psi_s
                        uFM_sing_inv(:,:,ising,1) = uFMInv
                        uFM_sing_init(:,:,ising,1) = uFM
                     ELSEIF (psiInters(iInterval,1) == 1) THEN
                     !Left edge of the the interval is a sing surface
                     !so it modifies the RIGHT side of psi_s
                        uFM_sing_inv(:,:,ising,2) = uFMInv
                        uFM_sing_init(:,:,ising,2) = uFM
                     ELSE
                        CALL program_stop("Unexpected psiInters!")
                     ENDIF
                  ENDIF
c-----------------------------------------------------------------------
c     perform the integration: parallel hamiltonian integration.
c-----------------------------------------------------------------------
                  iopt=1
                  itask=5
                  itol=2
                  mf=10
                  istate=1
                  iwork=0
                  rwork=0
                  rwork(1)=endPsi
                  rwork(5)=axisPsi*1e-3*psiDirs(iInterval)
                  rwork(11)=rwork(5)
                  CALL SYSTEM_CLOCK(COUNT=sTime)
                  DO WHILE(psiDirs(iInterval)*startPsi
     $                 < psiDirs(iInterval)*endPsi)
                     rtol(1) = tol_nr

                     DO ipert=1,2*mpert
                        atol0=MAXVAL(ABS(uFM(:,ipert)))*tol_nr
                        atol(:,ipert)=atol0 !CMPLX(atol0,atol0)
                     ENDDO

                     maxatol=MAXVAL(ABS(atol))*100
                     WHERE(atol == 0)
                        atol=maxatol !CMPLX(maxatol,maxatol)
                     ENDWHERE

                     CALL ZVODE1(sing_derFM,neq,uFM,startPsi,endPsi,
     $                    itol,rtol,atol,itask,istate,iopt,zwork,lzw,
     $                    rwork,lrw,iwork,liw,ode_nojac,mf,rpar,ipar)
                  ENDDO
                  IF ( grid_packing == "naive" ) THEN
                     CALL SYSTEM_CLOCK(COUNT=fTime)
                     WRITE(*,'(i3,a,f8.6,a,f8.6,a,f8.6)')
     $                    iInterval,",",t0,",",t1,",",
     $                    REAL(fTime-sTime,8)/REAL(cr,8)
                  ENDIF
               ENDIF
c-----------------------------------------------------------------------
c     store integration output for hamiltonian ODE.
c-----------------------------------------------------------------------
               uFM_all(:,:,iInterval) = uFM
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     terminate parallel integration.
c-----------------------------------------------------------------------
!$OMP END DO
!$OMP END PARALLEL
         CALL SYSTEM_CLOCK(COUNT=fTime)
         IF (verbose_performance_output) THEN
            print *,"*** ode-parallel-integration time=",
     $           REAL(fTime-sTime0,8)/REAL(cr,8)
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_run
c-----------------------------------------------------------------------
c     subprogram 2. ode_set_intervals.
c     divy ode subintervals, set bounds and bndry class (e.g. singular)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_set_intervals
         INTEGER :: i, j, jsing, iS
         INTEGER :: newItvls, cumulItvls, innerItvls
         REAL(r8) :: s1, s2, sLast, alpha
         REAL(r8), DIMENSION(:), ALLOCATABLE :: aT
         LOGICAL, DIMENSION(:), ALLOCATABLE :: aMask
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
         !psiPoints = start and end points of intervals
         ALLOCATE(psiPoints(nIntervalsTot,2))
         !psiInters = indication of whether edge L/R/LR singular
         !e.g. [0,1,2] = intvl LEFT of singsurf #2
         !e.g. [1,1,3] = intvl sandwiched between singsurf #3 and #4
         ALLOCATE(psiInters(nIntervalsTot,3))
         !the direction of integration for each intvl
         ALLOCATE(psiDirs(nIntervalsTot))

         !aT tracks the runtime expected of each interval
         ALLOCATE(aT(nIntervalsTot))
         !aMask sets which intervals can be divided
         ALLOCATE(aMask(nIntervalsTot))
c-----------------------------------------------------------------------
c     first set minimal intervals defined by axis and singular surfaces.
c-----------------------------------------------------------------------
         psiPoints(1,1) = axisPsi
         psiPoints(1,2) = scalc(1)%singEdgesLR(1)
         psiInters(1,:) = (/ 0, 1, 1 /)
         j = 1
         DO jsing = 1,msing
            j=j+1
            psiPoints(j,1) = scalc(jsing)%singEdgesLR(1)
            psiPoints(j,2) = scalc(jsing)%singEdgesLR(2)
            psiInters(j,:) = (/ 2, 2, jsing /)
            !2,2 above indicates the interval should not be split
            !because it's the "inner layer" at psi_s itself!

            j=j+1
            psiPoints(j,1) = scalc(jsing)%singEdgesLR(2)
            IF (jsing==msing) THEN
               psiPoints(j,2) = outerPsi
               psiInters(j,:) = (/ 1, 0, jsing /)
            ELSE
               psiPoints(j,2) = scalc(jsing+1)%singEdgesLR(1)
               psiInters(j,:) = (/ 1, 1, jsing /)
            ENDIF
         ENDDO
         nIntervals = j
c-----------------------------------------------------------------------
c     set aMask to control division of intervals.
c-----------------------------------------------------------------------
         aMask(:) = .FALSE.
         aMask(1:nIntervals) = .TRUE.
         SELECT CASE(grid_packing)
c-----------------------------------------------------------------------
c     uniformly divy up intervals.
c-----------------------------------------------------------------------
         CASE("naive")
            cumulItvls = 0
            sLast = psiPoints(nIntervals,2)
            DO iS=1,nIntervals
               s1 = psiPoints(1+cumulItvls,1)
               s2 = psiPoints(1+cumulItvls,2)
               IF (psiInters(1+cumulItvls,1) == 2) THEN
                  newItvls = 1 !Inner layer case
               ELSE
                  innerItvls = msing-NINT((iS-1)/2.0) !Inner layers
                  IF (iS .LE. nIntervals) THEN
                     newItvls = MAX(1,NINT((nIntervalsTot-innerItvls-
     $                    cumulItvls) * (s2-s1)/(sLast-s1)))
                  ELSE
                     newItvls = nIntervalsTot-cumulItvls
                  ENDIF
               ENDIF
               IF (newItvls > 1 .AND. aMask(iS)) THEN
                  !Copy data forward in array to avoid overwriting
                  DO i = cumulItvls+1+nIntervals-iS,cumulItvls+2,-1
                     psiPoints(i+newItvls-1,:) = psiPoints(i,:)
                     psiInters(i+newItvls-1,:) = psiInters(i,:)
                  ENDDO
                  !Now expand the current interval into these slots.
                  IF (iS == 1) THEN
                     psiInters(newItvls,:) = psiInters(1,:)
                     psiInters(2:newItvls-1,:) = 0
                     psiInters(1,:) = 0
                  ELSEIF (iS == nIntervals) THEN
                     psiInters(cumulItvls+2:nIntervalsTot,:) = 0
                  ELSE
                     psiInters(cumulItvls+newItvls,1) = 0
                     psiInters(cumulItvls+newItvls,2) = 1
                     psiInters(cumulItvls+newItvls,3) =
     $                    psiInters(cumulItvls+1,3) + 1
                     psiInters(cumulItvls+1,2) = 0 !(,1) & (,3) are fine
                     psiInters(cumulItvls+2:cumulItvls+newItvls-1,:) = 0
                  ENDIF
                  psiPoints(cumulItvls+newItvls,2) = s2
                  DO i=1,newItvls-1
                     psiPoints(i+cumulItvls,2) = s1 + i*(s2-s1)/newItvls
                     psiPoints(i+1+cumulItvls,1) =
     $                    psiPoints(i+cumulItvls,2)

                  ENDDO
               ENDIF
               cumulItvls = cumulItvls + newItvls
            ENDDO
            nIntervals = nIntervalsTot
            !DO i = 1,nIntervalsTot
            !   WRITE(*,'(3i2)')(psiInters(i,j),j=1,3)
            !ENDDO
            !DO i = 1,nIntervalsTot
            !   WRITE(*,'(2f9.6)')(psiPoints(i,j),j=1,2)
            !ENDDO
c-----------------------------------------------------------------------
c     use axis and singular surfaces to divy up intervals.
c-----------------------------------------------------------------------
         CASE("singularities")
c-----------------------------------------------------------------------
c     calc. itvl i approx. compute time: sum_k(psi_i-psi_k) k<-{a,s,e}
c-----------------------------------------------------------------------
            aT(:) = ode_etime()        !The time to first "enter" ZVODE.
            DO i = 1,nIntervals
               !Inner layers require minimal time, and needn't be split,
               !so we zero out their time requirement.
               IF (psiInters(i,1) == 2) THEN
                  aT(i) = 0
               ELSE
                  aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                 0.0_r8,"a",psiPoints(i,2))
                  DO iS = 1,msing
                     aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                    sing(iS)%psifac,"s",psiPoints(i,2))
                  ENDDO
                  aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                 1.0_r8,"e",psiPoints(i,2))
               ENDIF
            ENDDO
c-----------------------------------------------------------------------
c     subdivide max-compute-time interval until nIntervals=nIntervalsTot
c-----------------------------------------------------------------------
            DO WHILE (nIntervals < nIntervalsTot)
               j = MAXLOC(aT(1:nIntervals), DIM=1, MASK=aMask)
               psiPoints(j+2:nIntervals+1,:)=psiPoints(j+1:nIntervals,:)
               psiInters(j+2:nIntervals+1,:)=psiInters(j+1:nIntervals,:)
               aT(j+2:nIntervals+1)=aT(j+1:nIntervals)
               aMask(nIntervals+1)=.TRUE.
c-----------------------------------------------------------------------
c     recompute interval times.
c-----------------------------------------------------------------------
               s1 = ode_itime(psiPoints(j,1),0.0_r8,"a")
               s2 = ode_itime(psiPoints(j,2),0.0_r8,"a")
               DO iS = 1,msing
                  s1 = s1+ode_itime(psiPoints(j,1),sing(iS)%psifac,"s")
                  s2 = s2+ode_itime(psiPoints(j,2),sing(iS)%psifac,"s")
               ENDDO
               s1 = s1+ode_itime(psiPoints(j,1),1.0_r8,"e")
               s2 = s2+ode_itime(psiPoints(j,2),1.0_r8,"e")
               !This finds the points along the parallelogram that
               !evenly splits the area underneath it.
               alpha = (2*s1-SQRT(2*s1**2.0_r8+2*s2**2_r8))/(2*(s1-s2))

               psiPoints(j+1,2) = psiPoints(j,2)
               psiPoints(j+1,1) = alpha*psiPoints(j,2)+
     $              (1-alpha)*psiPoints(j,1)
               psiPoints(j,2) = psiPoints(j+1,1)

               aT(j:j+1) = 0.093
               DO i = j,j+1
                  aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                 0.0_r8,"a",psiPoints(i,2))
                  DO iS = 1,msing
                     aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                    sing(iS)%psifac,"s",psiPoints(i,2))
                  ENDDO
                  aT(i) = aT(i) + ode_itime(psiPoints(i,1),
     $                 1.0_r8,"e",psiPoints(i,2))
               ENDDO
c-----------------------------------------------------------------------
c     reassign interval types.
c-----------------------------------------------------------------------
               IF (psiInters(j,1) == 1) THEN
                  IF (psiInters(j,2) == 1) THEN
                     psiInters(j,2)   = 0
                     psiInters(j+1,1) = 0
                     psiInters(j+1,2) = 1
                     psiInters(j+1,3) = psiInters(j,3)+1
                  ELSE
                     psiInters(j+1,:) = 0
                  ENDIF
               ELSEIF(psiInters(j,2) == 1) THEN
                  psiInters(j+1,:) = psiInters(j,:)
                  psiInters(j,:) = 0
               ELSE
                  psiInters(j+1,:) = 0
               ENDIF

               nIntervals = nIntervals+1
            ENDDO
c-----------------------------------------------------------------------
c     error handling for unrecognized grid packing preference.
c-----------------------------------------------------------------------
         CASE default
            CALL program_stop("Cannot recognize grid_packing "//
     $           grid_type)
         END SELECT
c-----------------------------------------------------------------------
c     error-check for sufficient # of intervals.
c-----------------------------------------------------------------------
         DO i=1,nIntervals
            IF (psiInters(i,1) == 1 .AND. psiInters(i,2) == 1) THEN
               CALL program_stop("Not enough intervals were "//
     $              "requested to have each interval abut only one "//
     $              "singular surface. Increase nIntervals argument.")
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_set_intervals
c-----------------------------------------------------------------------
c     subprogram 3. ode_set_interval_details.
c     set integration directions, and interval #s of singular surfaces.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_set_interval_details
         INTEGER :: j, k, jsing
c-----------------------------------------------------------------------
c     Get singular intvl indices, "midway pts.", and integration dirs.
c-----------------------------------------------------------------------
         DO j=1,nIntervals
            psiDirs(j) = 1  !This is the integration direction
            jsing = psiInters(j,3)  !Either 0 or the # of the sing surf
            IF (psiInters(j,2) == 1) THEN
               !Right edge of intvl is sing => intvl # to LEFT of psi_s
               scalc(jsing)%singIntervalL = j
               !Now get midpoint of the singsurf idx and the PRIOR idx
               IF (jsing == 1) THEN
                  !Skew intvls b/w axis and psi_1 => fixup for [0, >>0]
                  scalc(1)%singMidPt = FLOOR(REAL(1.0_r8 +
     $                 (axis_mid_pt_skew-1.0_r8)*scalc(1)%singIntervalL)
     $                 / axis_mid_pt_skew)
               ELSEIF (jsing > 1) THEN
                  scalc(jsing)%singMidPt =
     $                 FLOOR(REAL(scalc(jsing-1)%singIntervalR +
     $                 scalc(jsing)%singIntervalL)   / 2  )
               ENDIF
               DO k = j,scalc(jsing)%singMidPt+1,-1
                  !Reverse integration dir on intvls before sing surf
                  psiDirs(k) = -1
               ENDDO
            ELSEIF (psiInters(j,1) == 1) THEN
               !Left edge of intvl is sing => intvl # to RIGHT of psi_s
               scalc(jsing)%singIntervalR = j
            ELSEIF (psiInters(j,1) == 2) THEN
               scalc(jsing)%singInterval0 = j
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
         RETURN
         END SUBROUTINE ode_set_interval_details
c-----------------------------------------------------------------------
c     subprogram 4. ode_propagate_LR.
c     propagate interval solutions across full domain, with fixups
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_propagate_LR
         COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: u_i, identityMat, tM
         INTEGER :: i,jsing,m2
         TYPE(sing_calculator), POINTER :: s, s1
         INTEGER :: sTime, fTime, cr
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
         CALL SYSTEM_CLOCK(COUNT=sTime)
         CALL SYSTEM_CLOCK(COUNT_RATE=cr)

         m2 = 2*mpert

         ALLOCATE(uAxis(m2,m2),uAxisD(m2,m2))

         identityMat = 0
         DO i = 1,m2
            identityMat(i,i) = 1
         ENDDO
         uAxis = 0
         DO i = 1,mpert
            uAxis(i+mpert,i+mpert) = 1 !q=0 at axis
         ENDDO
c-----------------------------------------------------------------------
c     solve the BVP to make Xi=0 at psi=0,1.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     shoot propagator from axis
c-----------------------------------------------------------------------
         DO i = 1,scalc(1)%singMidPt
            u_i = uFM_all(:, :, i)
            tM = uAxis
            CALL ZGEMM('N','N',m2,m2,m2,one_c,u_i,m2,tM,m2,zero_c,uAxis,
     $           m2)
            !Can afford to do fixup for the axis because it's just the
            !entire 40-D subspace that matters--free to take lin. combo.
            CALL ode_fixup(uAxis,uAxisD)
         ENDDO
c-----------------------------------------------------------------------
c     shoot propagators from L and R from each sing surf.
c-----------------------------------------------------------------------
         DO jsing = 1,msing
            s => scalc(jsing)
            ALLOCATE(s%uShootL(m2,m2),s%uShootR(m2,m2))
            s%uShootL = identityMat
            DO i = s%singIntervalL,s%singMidPt+1,-1
               u_i = uFM_all(:, :, i)
               tM = s%uShootL
               CALL ZGEMM('N','N',m2,m2,m2,one_c,u_i,m2,tM,m2,zero_c,
     $              s%uShootL,m2)
            ENDDO
            IF (jsing < msing) THEN
               s%uShootR = identityMat
               s1 => scalc(jsing+1)
               DO i = s%singIntervalR,s1%singMidPt
                  u_i = uFM_all(:, :, i)
                  tM = s%uShootR
                  CALL ZGEMM('N','N',m2,m2,m2,one_c,u_i,m2,tM,m2,zero_c,
     $                 s%uShootR,m2)
               ENDDO
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     shoot propagator from last singular surface to the edge.
c-----------------------------------------------------------------------
         s => scalc(msing)
         s%uShootR = identityMat
         DO i = s%singIntervalR,nIntervals
            u_i = uFM_all(:, :, i)
            tM = s%uShootR
            CALL ZGEMM('N','N',m2,m2,m2,one_c,u_i,m2,tM,m2,zero_c,
     $           s%uShootR,m2)
         ENDDO
c-----------------------------------------------------------------------
c     print timer.
c-----------------------------------------------------------------------
         CALL SYSTEM_CLOCK(COUNT=fTime)
         IF (verbose_performance_output) THEN
            print*,"*** ode-propagateLR time=",
     $           REAL(fTime-sTime,8) /REAL(cr,8)
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_propagate_LR
c-----------------------------------------------------------------------
c     subprogram 5. ode_calc_modes_for_wp.
c     link all propagations to solve for the modes at plasma edge.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_calc_modes_for_wp

      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: u_tot, A, tM
      INTEGER :: m2,i,info,jsing
      INTEGER, DIMENSION(2*mpert) :: ipiv
      TYPE(sing_calculator), POINTER :: s
      COMPLEX(r8), DIMENSION(2*mpert) :: ubig
      COMPLEX(r8), DIMENSION(2*mpert,2*mpert) :: umat
      COMPLEX(r8) :: uB, zdotc

      INTEGER :: sTime,fTime,cr
c-----------------------------------------------------------------------
c     initialize STATE TRANSITION PROPAGATION across full domain.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
      CALL SYSTEM_CLOCK(COUNT_RATE=cr)

      !Initialize u_tot with uAxis, which maps: axis => midpt1
      u_tot = uAxis
      m2 = 2*mpert
      !Loop over each singular surface.
      DO jsing = 1,msing
c-----------------------------------------------------------------------
c     propagate over interval left of the singular surface with fixup.
c-----------------------------------------------------------------------
         s => scalc(jsing)
         A = s%uShootL

         !This maps u_tot --> A^(-1) * u_tot:
         CALL ZGESV(m2,m2,A,m2,ipiv,u_tot,m2,info)
         CALL ode_fixup(u_tot)

         IF (asymp_at_sing) THEN
            !Invert non-Id initial FM at singular surface (left)
            !(Here, because we just INVERTED the FM, the init FM
            !is simply multiplied by it...)
            tM = u_tot
            umat = uFM_sing_init(:,:,jsing,1)
            !This maps u_tot --> umat*tM
            CALL ZGEMM('N','N',m2,m2,m2,one_c,umat,m2,tM,m2,zero_c,
     $           u_tot,m2)
            CALL ode_fixup(u_tot)
         ENDIF
         IF (kill_big_soln_for_ideal_dW) THEN
c-----------------------------------------------------------------------
c     Remove big solution on left of singular surface.
c-----------------------------------------------------------------------
            !If explicitly instructed to remove the big solution.
            IF (asymp_at_sing) THEN
               ubig = uFM_sing_init(:,s%sing_col,jsing,1)
               uB = zdotc(2*mpert,ubig,1,ubig,1)
               DO i=1,2*mpert
                  u_tot(:,i) = u_tot(:,i) - ubig *
     $                 zdotc(2*mpert,ubig,1,u_tot(:,i),1) / uB
               ENDDO
            ELSE
               u_tot(:,s%sing_col) = 0
               u_tot(s%sing_col,:) = 0
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     propagate over INNER layer at the singular surface with fixup.
c-----------------------------------------------------------------------
         !Multiply by the inner layer's fundamental matrix
         A = uFM_all(:, :, s%singInterval0)
         tM = u_tot
         !This maps u_tot --> A*tM
         CALL ZGEMM('N','N',m2,m2,m2,one_c,A,m2,tM,m2,zero_c,u_tot,m2)
         CALL ode_fixup(u_tot)
         IF (kill_big_soln_for_ideal_dW) THEN
c-----------------------------------------------------------------------
c     Remove big solution on right of singular surface.
c-----------------------------------------------------------------------
            !If explicitly instructed to remove the big solution.
            IF (asymp_at_sing) THEN
               ubig = uFM_sing_init(:,s%sing_col,jsing,2)
               uB = zdotc(2*mpert,ubig,1,ubig,1)
               DO i=1,2*mpert
                  u_tot(:,i) = u_tot(:,i) - ubig *
     $                 zdotc(2*mpert,ubig,1,u_tot(:,i),1) / uB
               ENDDO
            ELSE
               u_tot(:,s%sing_col) = 0
               u_tot(s%sing_col,:) = 0
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     propagate over interval right of the singular surface with fixup.
c-----------------------------------------------------------------------
         IF (asymp_at_sing) THEN
            !Invert non-Id initial FM at singular surface (right)
            umat = uFM_sing_inv(:,:,jsing,2)
            CALL ZGEMM('N','N',m2,m2,m2,one_c,s%uShootR,m2,umat,m2,
     $           zero_c,A,m2)
         ELSE
            A = s%uShootR
         ENDIF

         !Multiply by the right-side (of psi_s) fundamental matrix.
         tM = u_tot
         CALL ZGEMM('N','N',m2,m2,m2,one_c,A,m2,tM,m2,zero_c,u_tot,m2)
         CALL ode_fixup(u_tot)
      ENDDO
c-----------------------------------------------------------------------
c     save the final fundamental matrix (for mpert modes) to array u.
c-----------------------------------------------------------------------
      u(:,:,1) = u_tot(1:mpert,mpert+1:m2)
      u(:,:,2) = u_tot(mpert+1:m2,mpert+1:m2)
c-----------------------------------------------------------------------
c     print timer.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print*,"*** calc-modes-for-wp time=",
     $        REAL(fTime-sTime,8) / REAL(cr,8)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_calc_modes_for_wp
c-----------------------------------------------------------------------
c     subprogram 6. ode_set_delta_prime_edge
c     set the edge boundary condition for delta prime solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_set_delta_prime_edge(wv)
      COMPLEX(r8), DIMENSION(mpert,mpert), INTENT(IN), OPTIONAL :: wv
      INTEGER :: m2,i
      INTEGER :: sTime, fTime, cr
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=sTime)
      CALL SYSTEM_CLOCK(COUNT_RATE=cr)
      m2 = 2*mpert
      ALLOCATE(uEdge(m2,m2))
c-----------------------------------------------------------------------
c     set the edge condition.
c-----------------------------------------------------------------------
      uEdge = 0.0_r8
      IF (calc_dp_with_vac) THEN
         IF (PRESENT(wv)) THEN
            ![Q; P] = [Id; -W_V] in right half of matrix
            DO i = 1,mpert
               uEdge(i,i+mpert) = 1
            ENDDO
            uEdge(mpert+1:m2,mpert+1:m2) = -wv * psio**2.0_r8
         ELSE
            CALL program_stop("Attempting to calculate delta prime "//
     $           "with vaccuum edge conditions, but no wv available!")
         ENDIF
      ELSE
         DO i = 1,mpert
            uEdge(i+mpert,i+mpert) = 1 !q=0 at edge
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     print timer.
c-----------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=fTime)
      IF (verbose_performance_output) THEN
         print*,"*** dprime-set-edge-time=",
     $        REAL(fTime-sTime,8) / REAL(cr,8)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_set_delta_prime_edge
c-----------------------------------------------------------------------
c     subprogram 7. ode_calc_delta_prime.
c     solve BVP to calculate delta prime matrix [2*sing x 2*sing].
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_calc_delta_prime
         REAL(r8), DIMENSION(:,:), ALLOCATABLE :: atol

         INTEGER :: nMat, nSp, i, j, k, jsing, ksing, s2, m2
         INTEGER :: iChg, sparse_mode, info, ipert0, dRow
         REAL(r8) :: rcond, mnorm, ZLANGE
         INTEGER, DIMENSION(:), ALLOCATABLE :: isp, jsp
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: asp, b, x
         COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: M, MInv, TempMat
         INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv, ipivTemp
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: uwork, cwork,
     $        cworkTemp
         REAL(r8), DIMENSION(:), ALLOCATABLE :: rwork, rworkTemp
         TYPE(sparse_array) :: A
         TYPE(sing_calculator), POINTER :: s

         CHARACTER(LEN=3),DIMENSION(:,:), ALLOCATABLE :: imag_unit
         INTEGER :: sTime, fTime, cr
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
         CALL SYSTEM_CLOCK(COUNT=sTime)
         CALL SYSTEM_CLOCK(COUNT_RATE=cr)

         m2 = 2*mpert
         s2 = 2*msing
         ALLOCATE(delta_prime_mat(s2,s2),imag_unit(s2,s2))
         ALLOCATE(atol(2*mpert,2*mpert))

         !matrix size
         nMat = (2+2*s2)*mpert  !mpert at 0 and 1, 4x mpert per sing
         !# of sparse elements
         nSp = 2 * (m2*m2+m2*mpert) + !axis and edge
     $        s2 * (m2-2) +     !each side of each sing
     $        (msing-1) * m2*m2*2 + !intervals sandwiched b/w sings
     $        s2  !driven big solution coefficients (2sides * msing)

         ALLOCATE(b(nMat),x(nMat),isp(nSp),jsp(nSp),asp(nSp))
         ALLOCATE(TempMat(m2,m2),ipivTemp(m2),rworkTemp(2*m2),
     $        cworkTemp(2*m2))
         k = 0
c-----------------------------------------------------------------------
c     form BVP matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     form axis fundamental matrix.
c-----------------------------------------------------------------------
         DO i = 1,m2
            DO j = 1,mpert
               k=k+1
               isp(k) = i
               jsp(k) = j
               asp(k) = -uAxisD(i,j+mpert) !Only take q=0 modes (RHS)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     construct matrix by looping over each singular surface.
c-----------------------------------------------------------------------
         DO jsing = 1,msing
            s => scalc(jsing)
            ipert0 = s%sing_col
c-----------------------------------------------------------------------
c     left shoot.
c-----------------------------------------------------------------------
            DO i = 1,m2
               DO j = 1,m2
                  k=k+1
                  isp(k) = i + (2*m2-2) * (jsing-1)
                  jsp(k) = j + mpert + (2*m2) * (jsing-1)
                  asp(k) = s%uShootL(i,j)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     non-singular modes' continuity at surface.
c-----------------------------------------------------------------------
            iChg = 0
            DO i = 1,m2
               IF (i/=ipert0 .AND. i/=ipert0+mpert) THEN
                  k=k+1
                  isp(k) = i + iChg + m2 + (2*m2-2) * (jsing-1)
                  jsp(k) = i + mpert + (2*m2) * (jsing-1)
                  asp(k) = 1.0

                  k=k+1
                  isp(k) = isp(k-1)
                  jsp(k) = jsp(k-1) + m2
                  asp(k) = -1.0
               ELSE
                  iChg = iChg - 1
               ENDIF
            ENDDO
c-----------------------------------------------------------------------
c     right shoot.
c-----------------------------------------------------------------------
            DO i = 1,m2
               DO j = 1,m2
                  k=k+1
                  isp(k) = i + (2*m2-2) * jsing
                  jsp(k) = j + mpert + m2 + (2*m2)*(jsing-1)
                  asp(k) = -s%uShootR(i,j)
               ENDDO
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     form edge fundamental matrix.
c-----------------------------------------------------------------------
         DO i = 1,m2
            DO j = 1,mpert
               k=k+1
               isp(k) = i + (2*m2-2) * msing
               jsp(k) = j + mpert + (2*m2) * msing
               asp(k) = uEdge(i,j+mpert)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     set big solution coefficients for each singular surface.
c-----------------------------------------------------------------------
         DO jsing = 1,msing
            s => scalc(jsing)
            ipert0 = s%sing_col
            !left side
            k=k+1
            isp(k) = isp(k-1)+1
            jsp(k) = ipert0 + mpert + (2*m2)*(jsing-1)
            asp(k) = 1.0
            !right side
            k=k+1
            isp(k) = isp(k-1)+1
            jsp(k) = ipert0 + mpert + m2 + (2*m2)*(jsing-1)
            asp(k) = 1.0
         ENDDO
c-----------------------------------------------------------------------
c     error check.
c-----------------------------------------------------------------------
         IF (k .NE. nSp) THEN
            CALL program_stop("Incorrect # of sparse elements filled.")
         ENDIF
c-----------------------------------------------------------------------
c     check condition number of BVP matrix.
c-----------------------------------------------------------------------
         ALLOCATE(M(nMat,nMat),MInv(nMat,nMat))
         ALLOCATE(ipiv(nMat),uwork(nMat),cwork(2*nMat),rwork(2*nMat))
         M=0.0
         DO k = 1,nSp
            M(isp(k),jsp(k)) = asp(k)
         ENDDO
         MInv = M
         CALL ZGETRF(nMat,nMat,MInv,nMat,ipiv,info)
         mnorm = ZLANGE('O',nMat,nMat,MInv,nMat,rwork)
         CALL ZGECON('O',nMat,MInv,nMat,mnorm,rcond,cwork,rwork,info)
         print *,"Delta' sparse matrix condition number = ",rcond
c-----------------------------------------------------------------------
c     solve BVP matrix with sparse matrix.
c-----------------------------------------------------------------------
         IF (solve_delta_prime_with_sparse_mat) THEN
            CALL sparse_form_array(nSp, nMat, asp, isp, jsp, A)
            !Loop over big solution driving terms.
            sparse_mode = 1
            DO jsing = 1,msing
               DO i = 1,2
                  k = -mpert !Thus easy to add m2 each time (past uAxis)
                  b = 0.0
                  IF (i == 1) THEN
                     !Left side driving term
                     b(nMat-s2+2*jsing-1) = 1.0
                     dRow = 2*jsing-1
                  ELSE
                     !Right side driving term
                     b(nMat-s2+2*jsing) = 1.0
                     dRow = 2*jsing
                  ENDIF
                  CALL sparse_solve_A_x_equals_b(A, b, x, sparse_mode)
                  sparse_mode = 2
                  DO ksing = 1,msing
                     s => scalc(ksing)
                     ipert0 = s%sing_col

                     !Find the small solution coeffs (p mode, so +mpert)
                     k = k+m2
                     delta_prime_mat(dRow,2*ksing-1) = x(k+ipert0+mpert)
                     k = k+m2
                     delta_prime_mat(dRow,2*ksing) = x(k+ipert0+mpert)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
c-----------------------------------------------------------------------
c     solve BVP matrix without sparse matrix.
c-----------------------------------------------------------------------
            !Loop over big solution driving terms.
            DO jsing = 1,msing
               DO i = 1,2
                  b = 0.0
                  IF (i == 1) THEN
                     !Left side driving term
                     b(nMat-s2+2*jsing-1) = 1.0
                     dRow = 2*jsing-1
                  ELSE
                     !Right side driving term
                     b(nMat-s2+2*jsing) = 1.0
                     dRow = 2*jsing
                  ENDIF
                  x = b
                  CALL ZGETRS('N',nMat,1,MInv,nMat,ipiv,x,nMat,info)
                  IF (info /= 0) THEN
                     print *,"info was not zero! info=",info
                  ENDIF

                  k = -mpert
                  DO ksing = 1,msing
                     s => scalc(ksing)
                     ipert0 = s%sing_col

                     !Find the small solution coeffs (p mode, so +mpert)
                     k = k+m2
                     delta_prime_mat(dRow,2*ksing-1) = x(k+ipert0+mpert)
                     k = k+m2
                     delta_prime_mat(dRow,2*ksing) = x(k+ipert0+mpert)
                  ENDDO
c-----------------------------------------------------------------------
c     error-check non-sparse matrix solution.
c-----------------------------------------------------------------------
                  !Error-check big solution elements by backward solve.
                  b = ABS(MATMUL(M,x))
                  IF (ABS(1-SUM(ABS(b))) > big_soln_err_tol) THEN
                     print *,"Big sol'n error, surface #",jsing,"side",i
                     print *,"=",ABS(1-SUM(ABS(b)))
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         CALL ascii_open(delta_prime_out_unit,"delta_prime.out",
     $        "UNKNOWN")
         WRITE(delta_prime_out_unit,'(a)') "Delta_prime_Matrix:"
         imag_unit = "+i*"
         WHERE(AIMAG(delta_prime_mat)<0.)imag_unit = '-i*'
         DO i=1,s2
            write(delta_prime_out_unit,'(30(g12.5,a,g12.5,2x,a))')
     $           ( REAL(delta_prime_mat(i,j)),imag_unit(i,j),
     $           ABS(AIMAG(delta_prime_mat(i,j))),"| ", j=1,s2 )
         ENDDO
         CALL ascii_close(delta_prime_out_unit)
c-----------------------------------------------------------------------
c     print timer.
c-----------------------------------------------------------------------
         CALL SYSTEM_CLOCK(COUNT=fTime)
         IF (verbose_performance_output) THEN
            print*,"*** calc-dprime-time=",
     $           REAL(fTime-sTime,8) / REAL(cr,8)
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
         RETURN
      END SUBROUTINE ode_calc_delta_prime
c-----------------------------------------------------------------------
c     subprogram 8. ode_fixup.
c     performs Gaussian reduction of solution matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_fixup(uT, uTNorm)
         COMPLEX(r8), DIMENSION(2*mpert,2*mpert), INTENT(INOUT) :: uT
         COMPLEX(r8), DIMENSION(2*mpert,2*mpert),
     $        INTENT(INOUT), OPTIONAL :: uTNorm
         COMPLEX(r8), DIMENSION(2*mpert,mpert) :: uT2
         REAL(r8), DIMENSION(mpert) :: unorm
         LOGICAL, DIMENSION(2,mpert) :: mask
         INTEGER, DIMENSION(mpert) :: index
         INTEGER :: ipert, isol, kpert, ksol, jsol, jmax(1), i
         REAL(r8) :: colNorm
c-----------------------------------------------------------------------
c     calculate and sort the Euclidean norms of the propagator's columns
c-----------------------------------------------------------------------
         !u = [q; p] (at the axis)
         !Look at the q(axis) = 0 modes only = RHS of uT
         uT2 = uT(:,mpert+1:2*mpert)
         !Get the norm of the "q" halves of the propagator's columns
         unorm = SQRT( SUM(ABS(uT2(1:mpert,:))**2.0_r8, DIM=1) )
         !Sort the columns by their norms => store in "index"
         index(1:mpert)=(/(ipert,ipert=1,mpert)/)
         CALL bubble(unorm,index,1,mpert)
c-----------------------------------------------------------------------
c     triangularize primary solutions--(take linear combos of columns)
c-----------------------------------------------------------------------
         mask=.TRUE.
         DO isol=1,mpert
            ksol=index(isol) !ksol = largest remaining mode (column)
            mask(2,ksol)=.FALSE.
            jmax=MAXLOC(ABS(uT2(1:mpert,ksol)),mask(1,1:mpert))
            kpert=jmax(1) !kpert = largest unmasked row element in ksol
            mask(1,kpert)=.FALSE.
            DO jsol=1,mpert
               IF(mask(2,jsol))THEN
                  IF (uT2(kpert,ksol) == 0) THEN
                     IF (uT2(kpert,jsol)/=0) THEN
                        !There is an all zero column.
                        CALL program_stop("Unable to Gauss-reduce!")
                     ENDIF
                  ELSE
                     !The actual linear combination step.
                     !Remove from all other (as yet unmaksed) modes
                     !the projection of the largest mode on them.
                     uT2(:,jsol) = uT2(:,jsol) - uT2(:,ksol) *
     $                    uT2(kpert,jsol) / uT2(kpert,ksol)
                     uT2(kpert,jsol) = 0 !Just 0's an element that is
                                         !already meant to be exactly 0.
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         uT(:, 1:mpert) = 0.0
         uT(:, mpert+1:2*mpert) = uT2

         IF (PRESENT(uTNorm)) THEN
            uTNorm = uT
            DO i = mpert+1,2*mpert
               colNorm = SQRT(SUM(ABS(uT(:,i))**2.0_r8, DIM=1))
               uTNorm(:,i) = uT(:,i) / colNorm
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_fixup
c-----------------------------------------------------------------------
c     subprogram 9. ode_power_spacing.
c     generate interval mesh points according to a power law
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_power_spacing(pi,pf,power,n,pout)
         REAL(r8), INTENT(IN) :: pi, pf
         INTEGER, INTENT(IN) :: power, n
         REAL(r8), DIMENSION(n), INTENT(OUT) :: pout

         REAL(r8), DIMENSION(n) :: x, y
         INTEGER :: i
c-----------------------------------------------------------------------
c     calculate.
c-----------------------------------------------------------------------
         x = REAL((/(i,i=0,n-1)/),8) / REAL(n-1,8)
         IF (pi < pf) THEN
            y = x**power
         ELSE
            y = (1-x)**power
         ENDIF
         pout = pi + (pf-pi)*y
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_power_spacing
c-----------------------------------------------------------------------
c     subprogram 10. ode_pseudo_inv.
c     calculate the pseudoinverse of a singular matrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_pseudo_inv(uIn,uOut)
         COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: uIn
         COMPLEX(r8), DIMENSION(:,:), INTENT(INOUT) :: uOut
         COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: U,Vt
         COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: WORK,RWORK
         REAL(r8), DIMENSION(:), ALLOCATABLE :: S
         INTEGER :: M,N,K,K5,L,LWORK,INFO,iS

         M = SIZE(uIn,1)
         N = SIZE(uIn,2)
         K = MIN(M,N)
         L = MAX(M,N)
         LWORK = MAX(1,2*K+L)
         K5 = 5*K
         ALLOCATE(S(K),U(M,K),Vt(K,N),WORK(LWORK),RWORK(K5))
c-----------------------------------------------------------------------
c     perform SVD, uIn = U*S*Vt.
c-----------------------------------------------------------------------
         CALL ZGESVD('A','A',M,N,uIn,M,S,U,M,Vt,K,WORK,LWORK,RWORK,INFO)
c-----------------------------------------------------------------------
c     zscal scales the vector U(1,iS) by a constant.
c-----------------------------------------------------------------------
         DO iS = 1,K
            IF(iS==K)THEN
               CALL ZSCAL(M, DCMPLX(0.0), U(1,iS), 1)
            ELSE
               CALL ZSCAL(M, DCMPLX(1/ S(iS)), U(1,iS), 1)
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     return uOut = 1.0 * ((V)t)t * (U)t + 0*uOut.
c-----------------------------------------------------------------------
         CALL ZGEMM('C','C',N,M,K,one_c,Vt,K,U,M,zero_c,uOut,N)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_pseudo_inv
c-----------------------------------------------------------------------
c     subprogram 11. ode_nojac.
c     a dummy Jacobian call for the LSODE/ZVODE package
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_nojac(neq_jac, t, y, ml, mu, pd, nrpd)
         INTEGER  neq_jac, ml, mu, nrpd
         REAL*8  t, y, pd(nrpd,2)
         pd(:,:) = 0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_nojac
c-----------------------------------------------------------------------
c     subprogram 12. ode_itime.
c     compute time for interval: [pt1, pt2] OR instantaneous time [~pt2]
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION ode_itime(pt1, refpt, mode, pt2)
         REAL(r8), INTENT(IN) :: pt1, refpt
         CHARACTER, INTENT(IN) :: mode
         REAL(r8), INTENT(IN), OPTIONAL :: pt2
         REAL(r8) :: ode_itime
         REAL(r8) :: a, b

         SELECT CASE(mode)
         CASE("a")
            !axis
            a = REAL(39695,8)
            b = REAL(212830,8)
         CASE("s")
            !singular surface
            a = REAL(17147,8)
            b = REAL(470710,8)
         CASE("e")
            !edge
            a = REAL(1646,8)
            b = REAL(4683,8)
         CASE DEFAULT
            CALL program_stop("Unknown itime mode")
         END SELECT

         IF (PRESENT(pt2)) THEN
            !Integration of the itime form factor on interval.
            ode_itime = (a/b) * ABS(LOG(1+b*ABS(pt2-refpt))-
     $           LOG(1+b*ABS(pt1-refpt)))
         ELSE
            ode_itime = a/(1+b*ABS(pt1-refpt))
         ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION ode_itime
c-----------------------------------------------------------------------
c     subprogram 13. ode_etime.
c     returns the unavoidable time spent entering a ZVODE1 call.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION ode_etime()
         REAL(r8) :: ode_etime
         ode_etime = 0.015523
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION ode_etime
      END MODULE ode_mod
