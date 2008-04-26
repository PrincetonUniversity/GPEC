c-----------------------------------------------------------------------
c     file resist.f
c     studies resistive modes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. resist_mod.
c     1. resist_deltac_eval.
c     2. resist_deltac_find.
c     3. resist_write.
c     4. resist_nyquist.
c     5. resist_scan.
c     6. resist_root_eval.
c     7. resist_root_find.
c-----------------------------------------------------------------------
c     subprogram 0. resist_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE resist_mod
      USE matrix_mod
      IMPLICIT NONE

      INTEGER :: jsing

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. resist_deltac_eval.
c     evaluates imaginary part of odd inner delta.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION resist_deltac_eval(s_im) RESULT(delta_im)

      REAL(r8), INTENT(IN) :: s_im
      REAL(r8) :: delta_im

      COMPLEX(r8) :: s
      COMPLEX(r8), DIMENSION(2) :: deltar
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      s=ifac*s_im
      CALL deltar_run(singtype(jsing)%restype,s,deltar)
      delta_im=AIMAG(deltar(2))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION resist_deltac_eval
c-----------------------------------------------------------------------
c     subprogram 2. resist_deltac_find.
c     finds critical delta.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_deltac_find

      INTEGER :: it
      REAL(r8) :: root,err
      REAL(r8), PARAMETER :: power=1._r8/3._r8
      COMPLEX(r8) :: s
      COMPLEX(r8), DIMENSION(2) :: deltar
c-----------------------------------------------------------------------
c     compute derived parameters.
c-----------------------------------------------------------------------
      singtype(jsing)%restype%dr
     $     =singtype(jsing)%restype%e
     $     +singtype(jsing)%restype%f
     $     +singtype(jsing)%restype%h**2
      singtype(jsing)%restype%di
     $     =singtype(jsing)%restype%dr
     $     -(singtype(jsing)%restype%h-.5)**2
      singtype(jsing)%restype%sfac
     $     =singtype(jsing)%restype%taur
     $     /singtype(jsing)%restype%taua
      RETURN
c-----------------------------------------------------------------------
c     compute derived parameters.
c-----------------------------------------------------------------------
      IF(singtype(jsing)%restype%dr < 0)THEN
         root=((singtype(jsing)%restype%dr
     $        /singtype(jsing)%restype%taua)**2
     $        /singtype(jsing)%restype%taur)**power
         CALL droot_newton(resist_deltac_eval,root,err,it)
         s=ifac*root
         CALL deltar_run(singtype(jsing)%restype,s,deltar)
         singtype(jsing)%restype%deltac=REAL(deltar(2))
      ELSE
         singtype(jsing)%restype%deltac=-1
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_deltac_find
c-----------------------------------------------------------------------
c     subprogram 3. resist_write.
c     writes singular surface parameters.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_write

      CHARACTER(132) :: form0
      CHARACTER(32) :: form1
      INTEGER :: ising,nsing
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(3x,"variable",2x,',i2.2,'(3x,"sing ",i2.2,1x)/)')
 20   FORMAT('(a11,2x,1p,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     format statements and header.
c-----------------------------------------------------------------------
      nsing=MIN(msing,11)
      WRITE(form0,10)nsing
      WRITE(form1,20)nsing
      WRITE(out_unit,form0)(ising,ising=1,nsing)
c-----------------------------------------------------------------------
c     location.
c-----------------------------------------------------------------------
      WRITE(out_unit,form1)"q",
     $     (singtype(ising)%q,ising=1,nsing)
      WRITE(out_unit,form1)"q1",
     $     (singtype(ising)%q1,ising=1,nsing)
      WRITE(out_unit,form1)"psi",
     $     (singtype(ising)%psifac,ising=1,nsing)
      WRITE(out_unit,form1)"rho",
     $     (SQRT(singtype(ising)%psifac),ising=1,nsing)
c$$$c-----------------------------------------------------------------------
c$$$c     pressure-curvature parameters E, F, H, DR, DI, power.
c$$$c-----------------------------------------------------------------------
c$$$      WRITE(out_unit,'()')
c$$$      WRITE(out_unit,form1)"E",
c$$$     $     (singtype(ising)%restype%e,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"F",
c$$$     $     (singtype(ising)%restype%f,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"H",
c$$$     $     (singtype(ising)%restype%h,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"DR",
c$$$     $     (singtype(ising)%restype%dr,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"DI",
c$$$     $     (singtype(ising)%restype%di,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"power",
c$$$     $     (SQRT(-singtype(ising)%restype%di),ising=1,nsing)
c$$$c-----------------------------------------------------------------------
c$$$c     other GGJ parameters M, G, K.
c$$$c-----------------------------------------------------------------------
c$$$      WRITE(out_unit,'()')
c$$$      WRITE(out_unit,form1)"M",
c$$$     $     (singtype(ising)%restype%m,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"G",
c$$$     $     (singtype(ising)%restype%g,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"K",
c$$$     $     (singtype(ising)%restype%k,ising=1,nsing)
c$$$c-----------------------------------------------------------------------
c$$$c     scale factors eta, rho, taua,taur, S.
c$$$c-----------------------------------------------------------------------
c$$$      WRITE(out_unit,'()')
c$$$      WRITE(out_unit,form1)"eta",
c$$$     $     (singtype(ising)%restype%eta,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"rho",
c$$$     $     (singtype(ising)%restype%rho,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"taua",
c$$$     $     (singtype(ising)%restype%taua,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"taur",
c$$$     $     (singtype(ising)%restype%taur,ising=1,nsing)
c$$$      WRITE(out_unit,form1)"S",
c$$$     $     (singtype(ising)%restype%sfac,ising=1,nsing)
c$$$c-----------------------------------------------------------------------
c$$$c     scale deltap1.
c$$$c-----------------------------------------------------------------------
c$$$      DO ising=1,nsing
c$$$         deltap1(ising)=deltap1(ising)*singtype(ising)%psifac
c$$$     $        **(2*SQRT(-singtype(ising)%restype%di))
c$$$      ENDDO
c-----------------------------------------------------------------------
c     deltap.
c-----------------------------------------------------------------------
      WRITE(out_unit,'()')
      WRITE(out_unit,form1)"re deltap1",
     $     (REAL(deltap1(ising)),ising=1,nsing)
      WRITE(out_unit,form1)"im deltap1",
     $     (AIMAG(deltap1(ising)),ising=1,nsing)
      WRITE(*,form1)" re deltap1",
     $     (REAL(deltap1(ising)),ising=1,nsing)
      WRITE(*,form1)" im deltap1",
     $     (AIMAG(deltap1(ising)),ising=1,nsing)
c$$$      WRITE(out_unit,form1)"re deltap2",
c$$$     $     (REAL(deltap2(ising)),ising=1,nsing)
c$$$      WRITE(out_unit,form1)"im deltap2",
c$$$     $     (AIMAG(deltap2(ising)),ising=1,nsing)
c$$$      WRITE(*,form1)" re deltap2",
c$$$     $     (REAL(deltap2(ising)),ising=1,nsing)
c$$$      WRITE(*,form1)" im deltap2",
c$$$     $     (AIMAG(deltap2(ising)),ising=1,nsing)
c$$$      WRITE(out_unit,form1)"deltac",
c$$$     $     (singtype(ising)%restype%deltac,ising=1,nsing)
c-----------------------------------------------------------------------
c     trailer.
c-----------------------------------------------------------------------
      WRITE(out_unit,'()')
      WRITE(out_unit,form0)(ising,ising=1,nsing)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_write
c-----------------------------------------------------------------------
c     subprogram 4. resist_nyquist.
c     draws nyquist plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_nyquist(big,small,ns)

      REAL(r8), INTENT(IN) :: big,small
      INTEGER, INTENT(IN) :: ns

      COMPLEX(r8) :: s,ds
c-----------------------------------------------------------------------
c     first branch, positive imaginary axis.
c-----------------------------------------------------------------------
      WRITE(*,*)"Constructing Nyquist Plot"
      CALL bin_open(bin_unit,"scan.bin","UNKNOWN","REWIND","none")
      s=big*EXP(ifac*pi/2)
      ds=(small/big)**(1._r8/ns)
      CALL resist_scan(s,ds,ns)
c-----------------------------------------------------------------------
c     second branch, inner semicircle.
c-----------------------------------------------------------------------
      s=small*EXP(ifac*pi/2)
      ds=EXP(-ifac*pi/ns)
      CALL resist_scan(s,ds,ns)
c-----------------------------------------------------------------------
c     third branch, negative imaginary axis.
c-----------------------------------------------------------------------
      s=small*EXP(-ifac*pi/2)
      ds=(big/small)**(1._r8/ns)
      CALL resist_scan(s,ds,ns)
c-----------------------------------------------------------------------
c     fourth branch, outer semicircle.
c-----------------------------------------------------------------------
      s=big*EXP(-ifac*pi/2)
      ds=EXP(ifac*pi/ns)
      CALL resist_scan(s,ds,ns)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_nyquist
c-----------------------------------------------------------------------
c     subprogram 5. resist_scan
c     draws one branch of nyquist plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_scan(s,ds,ns)

      COMPLEX(r8), INTENT(INOUT) :: s
      COMPLEX(r8), INTENT(IN) :: ds
      INTEGER, INTENT(IN) :: ns

      INTEGER :: is,ising
      COMPLEX(r8) :: det
      COMPLEX(r8), DIMENSION(msing,2) :: deltar
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      DO is=0,ns
         DO ising=1,msing
            CALL deltar_run(singtype(ising)%restype,s,deltar(ising,:))
         ENDDO
         CALL matrix_poly_eval(deltar,det)
         WRITE(bin_unit)
     $        REAL(REAL(s),4),
     $        REAL(AIMAG(s),4),
     $        REAL(REAL(det),4),
     $        REAL(AIMAG(det),4),
     $        REAL(asinh(REAL(det)),4),
     $        REAL(asinh(AIMAG(det)),4)
         s=s*ds
      ENDDO
      WRITE(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_scan
c-----------------------------------------------------------------------
c     subprogram 6. resist_root_eval.
c     evaluates matching determinant for root finding.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION resist_root_eval(s) RESULT(det)

      COMPLEX(r8), INTENT(IN) :: s
      COMPLEX(r8) :: det

      INTEGER :: ising
      COMPLEX(r8), DIMENSION(msing,2) :: deltar
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      DO ising=1,msing
         CALL deltar_run(singtype(ising)%restype,s,deltar(ising,:))
      ENDDO
      CALL matrix_poly_eval(deltar,det)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION resist_root_eval
c-----------------------------------------------------------------------
c     subprogram 7. resist_root_find.
c     finds root of matching determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE resist_root_find(root)

      COMPLEX(r8), INTENT(INOUT) :: root

      INTEGER :: it
      REAL(r8) :: err
c-----------------------------------------------------------------------
c     computation.
c-----------------------------------------------------------------------
      WRITE(*,*)"Finding root"
      CALL zroot_newton(resist_root_eval,root,err,it)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE resist_root_find
      END MODULE resist_mod
