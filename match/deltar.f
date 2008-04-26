c-----------------------------------------------------------------------
c     program deltar.
c     integrates the 4th-order fourier-transformed equations for the 
c     singular layer from large values to small values of the 
c     independent variable t, fits the solutions to asymptotic 
c     solutions to evaluate asymptotic coefficients, and uses these 
c     coefficients to construct the asymptotic ratios dlrp and dlrm.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. deltar_mod.
c     1. deltar_run.
c     2. deltar_der.
c     3. deltar_upsfit.
c     4. deltar_vpsfit.
c     5. deltar_origin.
c     6. deltar_infinity.
c     7. deltar_ratio.
c-----------------------------------------------------------------------
c     subprogram 0. deltar_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE deltar_mod
      USE gamma_mod
      IMPLICIT NONE

      LOGICAL, PARAMETER, PRIVATE :: out=.FALSE.,bin=.FALSE.,
     $     out2=.FALSE.
      INTEGER, PARAMETER, PRIVATE :: neqc=4,nsol=2,nmax=8,nstep=500,
     $     nps=8,nxps=1,ntmax=1
      REAL(r8), PARAMETER, PRIVATE :: fmin=1,fmax=1,rtol=1e-3,atol=1e-6

      REAL(r8), PRIVATE :: ee,ff,hh,mm,gg,kk,di,dr,tmin,tmax,pplus,p1
      REAL(r8), DIMENSION(neqc), PRIVATE :: al0
      COMPLEX(r8), PRIVATE :: q
      COMPLEX(r8), DIMENSION(6), PRIVATE :: dtr
      COMPLEX(r8), DIMENSION(neqc), PRIVATE :: pexp,bl1
      COMPLEX(r8), DIMENSION(neqc,neqc), PRIVATE :: v,d1,d1inv,u,up,
     $     al1,d0,d0inv
      COMPLEX(r8), DIMENSION(neqc,neqc,nmax), PRIVATE :: ups,vps
      
      TYPE :: resist_type
      REAL(r8) :: e,f,h,m,g,k,eta,rho,taua,taur,di,dr,sfac,deltac
      END TYPE resist_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. deltar_run.
c     sets up and integrates odes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_run(restype,s,deltar)
      
      TYPE(resist_type), INTENT(IN) :: restype
      COMPLEX(r8), INTENT(IN) :: s
      COMPLEX(r8), DIMENSION(2), INTENT(OUT) :: deltar

      INTEGER, PARAMETER :: neq=2*neqc*nsol,liw=20,lrw=22+2*neq**2*nsol
      INTEGER :: i,iflag,iopt,istate,itask,itol,j,jac,mf,istep
      INTEGER, DIMENSION(liw) :: iwork=0
      REAL(r8) :: t,tout,sfac,taua,taur,x0,q0
      REAL(r8), DIMENSION(lrw) :: rwork=0
      COMPLEX(r8), DIMENSION(4,2) :: c0
      COMPLEX(r8), DIMENSION(neqc,nsol) :: c1,y
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/1x,"ist",6x,"hu",9x,"t",4x,4(5x,"y(",i1,")",2x),
     $     4(4x,"c0(",i1,")",2x)/)
 20   FORMAT(i4,1p,10e11.3,i4)
c-----------------------------------------------------------------------
c     copy input parameters.
c-----------------------------------------------------------------------
      ee=restype%e
      ff=restype%f
      hh=restype%h
      mm=restype%m
      gg=restype%g
      kk=restype%k
      taua=restype%taua
      taur=restype%taur
c-----------------------------------------------------------------------
c     define curvatures and related parameters.
c-----------------------------------------------------------------------
      dr=ee+ff+hh*hh
      di=dr-(hh-.5)**2
      p1=SQRT(-di)
c-----------------------------------------------------------------------
c     define scale factors.
c-----------------------------------------------------------------------
      sfac=taur/taua
      x0=sfac**(-1._r8/3._r8)
      q0=x0/taua
      q=s/q0
c-----------------------------------------------------------------------
c     set up coefficient arrays.
c-----------------------------------------------------------------------
      CALL deltar_origin
      CALL deltar_infinity
c-----------------------------------------------------------------------
c     initialize at t=tmax.
c-----------------------------------------------------------------------
      CALL deltar_vpsfit(tmax,y,c1)
      t=tmax
      tout=tmin
      y=MATMUL(d1,v(:,1:nsol))
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istep=0
      istate=1
      iflag=1
      itask=5
      iopt=1
      itol=1
      mf=10
      rwork(1)=tmin
      rwork(11)=0
c-----------------------------------------------------------------------
c     print heading.
c-----------------------------------------------------------------------
      IF(out)WRITE(out_unit,10)(i,i=1,neqc),(i,i=1,neqc)
c-----------------------------------------------------------------------
c     compute and print output for each step.
c-----------------------------------------------------------------------
      DO
         IF(out .OR. bin .OR. t <= tmin)
     $        CALL deltar_upsfit(t,y,c0)
         IF(out .OR. bin)CALL deltar_vpsfit(t,y,c1)
         IF(out)WRITE(out_unit,20)istep,rwork(11),t,
     $        (REAL(y(j,1),4),j=1,neqc),(REAL(c0(j,1),4),j=1,neqc)
         IF(bin)WRITE(bin_unit)REAL(t,4),(REAL(y(j,1),4),j=1,neqc),
     $        (REAL(c0(j,1),4),j=1,neqc),(REAL(c1(j,1),4),j=1,neqc)
c-----------------------------------------------------------------------
c     advance ordinary differential equations.
c-----------------------------------------------------------------------
         IF(t <= tmin .OR. istep > nstep .OR. istate < 0)EXIT
         istep=istep+1
         CALL lsode(deltar_der,neq,y,t,tout,itol,rtol,atol,itask,istate,
     $        iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     compute asymptotic matching ratios and determinant factors.
c-----------------------------------------------------------------------
      IF(bin)WRITE(bin_unit)
      IF(out)WRITE(out_unit,10)(i,i=1,neqc),(i,i=1,neqc)
      CALL deltar_ratio(c0,deltar)
      deltar=deltar*sfac**(2*p1/3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_run
c-----------------------------------------------------------------------
c     subprogram 2. deltar_der.
c     contains the differential equations for the resistive layer.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_der(n,t,y,dy)
      
      INTEGER, INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(OUT) :: dy
      
      INTEGER :: i
c-----------------------------------------------------------------------
c     evaluate derivatives, new way.
c-----------------------------------------------------------------------
      dy=MATMUL(al1,y)*t
      DO i=1,neqc
         dy(i,:)=dy(i,:)+al0(i)*y(i,:)/t
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_der
c-----------------------------------------------------------------------
c     subprogram 3. deltar_upsfit.
c     evaluates convergent power series expansion about origin
c     and asympotic coefficients for this expansion.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_upsfit(t,y,c0)
      
      REAL(r8), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(OUT) :: c0
      
      INTEGER :: i,j,k,info
      INTEGER, DIMENSION(neqc) :: ipiv
      REAL(r8) :: t2
      COMPLEX(r8), DIMENSION(neqc,neqc) :: u1
c-----------------------------------------------------------------------
c     evaluate series solutions.
c-----------------------------------------------------------------------
      t2=t*t
      DO j=1,neqc
         DO i=1,neqc
            u(i,j)=ups(i,j,nps)
            DO k=nps-1,1,-1
               u(i,j)=u(i,j)*t2+ups(i,j,k)
            ENDDO
         ENDDO
         u(:,j)=u(:,j)*t**al0(j)
      ENDDO
c-----------------------------------------------------------------------
c     compute c0=uinv*y.
c-----------------------------------------------------------------------
      u1=u
      c0=y
      CALL zgetrf(neqc,neqc,u1,neqc,ipiv,info)
      CALL zgetrs('N',neqc,nsol,u1,neqc,ipiv,c0,neqc,info)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_upsfit
c-----------------------------------------------------------------------
c     subprogram 4. deltar_vpsfit.
c     evaluates asymptotic power series expansion about infinity
c     and asympotic coefficients for this expansion.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_vpsfit(t,y,c1)
      
      REAL(r8), INTENT(IN) :: t
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(IN) :: y
      COMPLEX(r8), DIMENSION(neqc,nsol), INTENT(OUT) :: c1
      
      INTEGER :: i,j,k,info
      INTEGER, DIMENSION(neqc) :: ipiv
      REAL(r8) :: err,t2,tfac1,trmrat,trmrat1,trmrat2
      COMPLEX(r8) :: tfac
      COMPLEX(r8), DIMENSION(neqc,neqc) :: term,term1,v1
c-----------------------------------------------------------------------
c     evaluate asymptotic series.
c-----------------------------------------------------------------------
      t2=1/(t*t)
      v=0
      DO j=1,neqc
         tfac=EXP(.5*bl1(j)*t*t+pexp(j)*LOG(t))
         DO i=1,neqc
            tfac1=tfac
            trmrat=0
            trmrat2=0
            k=1
            DO while(k <= nxps .AND. (trmrat <= 1 .OR. trmrat2 <= 1))
               term1(i,j)=term(i,j)
               term(i,j)=tfac1*vps(i,j,k)
               trmrat1=trmrat
               IF(term1(i,j) /= 0)trmrat=ABS(term(i,j)/term1(i,j))
               IF(trmrat1 /= 0)trmrat2=trmrat/trmrat1
               IF(trmrat <= 1 .OR. trmrat2 <= 1)THEN
                  v(i,j)=v(i,j)+term(i,j)
                  IF(v(i,j) /= 0)err=ABS(term(i,j)/v(i,j))
               ENDIF
               tfac1=tfac1*t2
               k=k+1
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute c1=vinv*y.
c-----------------------------------------------------------------------
      v1=v
      c1=MATMUL(d1inv,y)
      CALL zgetrf(neqc,neqc,v1,neqc,ipiv,info)
      CALL zgetrs('N',neqc,nsol,v1,neqc,ipiv,c1,neqc,info)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_vpsfit
c-----------------------------------------------------------------------
c     subprogram 5. deltar_origin.
c     evaluates arrays relating to the origin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_origin
      
      LOGICAL :: out
      INTEGER :: i,j,k,l
      REAL(r8) :: aminus,aplus,pminus,unmax
      COMPLEX(r8) :: a1,a2,kk1,q3
      COMPLEX(r8), DIMENSION(4,4) :: temp
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/46x,"al0"//24x,1p,4e11.3///46x,"al1"//
     $     4(4(1p,2e11.3,2x)/)/)
 20   FORMAT(47x,"d0"//4(4(1p,2e11.3,2x)/)/)
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      pplus=-.5+p1
      pminus=-.5-p1
      q3=q**3
      kk1=1+kk*q3
      a1=(gg-kk*ee)*q*q/kk1
      a2=SQRT(.5/(p1*q))
      aplus=hh-.5+p1
      aminus=hh-.5-p1
      out=out .OR. out2
c-----------------------------------------------------------------------
c     compute diagonalizing matrix d0.
c-----------------------------------------------------------------------
      d0(1,1)=q*a2
      d0(1,2)=1
      d0(1,3)=-q*a2
      d0(1,4)=0
      d0(2,1)=0
      d0(2,2)=1
      d0(2,3)=0
      d0(2,4)=0
      d0(3,1)=aplus*a2
      d0(3,2)=0
      d0(3,3)=-aminus*a2
      d0(3,4)=0
      d0(4,1)=-aplus*a2
      d0(4,2)=-a1
      d0(4,3)=aminus*a2
      d0(4,4)=1
c-----------------------------------------------------------------------
c     compute diagonalizing matrix d0inv.
c-----------------------------------------------------------------------
      d0inv(1,1)= d0(3,3)
      d0inv(2,2)= d0(4,4)
      d0inv(1,2)= d0(4,3)
      d0inv(2,1)= d0(3,4)
      d0inv(3,3)= d0(1,1)
      d0inv(4,4)= d0(2,2)
      d0inv(3,4)= d0(2,1)
      d0inv(4,3)= d0(1,2)
      d0inv(1,3)=-d0(1,3)
      d0inv(2,4)=-d0(2,4)
      d0inv(1,4)=-d0(2,3)
      d0inv(2,3)=-d0(1,4)
      d0inv(3,1)=-d0(3,1)
      d0inv(4,2)=-d0(4,2)
      d0inv(3,2)=-d0(4,1)
      d0inv(4,1)=-d0(3,2)
c-----------------------------------------------------------------------
c     compute zeroth order coefficient matrix elements al0(i).
c-----------------------------------------------------------------------
      al0(1)=p1
      al0(2)=.5
      al0(3)=-p1
      al0(4)=-.5
c-----------------------------------------------------------------------
c     zero first order coefficient matrix elements.
c-----------------------------------------------------------------------
      DO j=1,neqc
         DO i=1,neqc
            al1(i,j)=0
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute first order coefficient matrix elements.
c-----------------------------------------------------------------------
      al1(1,3)=1
      al1(2,4)=-kk1
      al1(3,1)=q
      al1(4,2)=-q/kk1
c-----------------------------------------------------------------------
c     transform first order coefficient matrix elements.
c-----------------------------------------------------------------------
      DO j=1,neqc
         DO i=1,neqc
            temp(i,j)=0
            DO k=1,neqc
               temp(i,j)=temp(i,j)+al1(i,k)*d0(k,j)
            ENDDO
         ENDDO
      ENDDO
      DO j=1,neqc
         DO i=1,neqc
            al1(i,j)=0
            DO k=1,neqc
               al1(i,j)=al1(i,j)+d0inv(i,k)*temp(k,j)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute lowest-order power series coefficient matrix elements ups.
c-----------------------------------------------------------------------
      DO i=1,neqc
         DO k=1,nps
            DO j=1,neqc
               ups(i,j,k)=0
            ENDDO
         ENDDO
         ups(i,i,1)=1
      ENDDO
c-----------------------------------------------------------------------
c     compute higher-power series coefficient matrix elements ups.
c-----------------------------------------------------------------------
      DO k=2,nps
         DO j=1,neqc
            DO i=1,neqc
               DO l=1,neqc
                  ups(i,j,k)=ups(i,j,k)+al1(i,l)*ups(l,j,k-1)
               ENDDO
               ups(i,j,k)=ups(i,j,k)/(al0(j)-al0(i)+2*(k-1))
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute tmin.
c-----------------------------------------------------------------------
      unmax=0
      DO j=1,neqc
         DO i=1,neqc
            unmax=MAX(unmax,ABS(ups(i,j,nps)))
         ENDDO
      ENDDO
      tmin=fmin*.5*(rtol/unmax)**(.5/nps)
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      IF(out2)THEN
         WRITE(out_unit,10)al0,((al1(i,j),j=1,neqc),i=1,neqc)
         WRITE(out_unit,20)((d0(i,j),j=1,neqc),i=1,neqc)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_origin
c-----------------------------------------------------------------------
c     subprogram 6. deltar_infinity.
c     evaluates arrays relating to the point at infinity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_infinity
      
      INTEGER :: i,j,j1,j2,j3,j4,jlam,k,l
      COMPLEX(r8) :: a1,bb,cc,dd,ddsq,kk1,lamda,lamdaq,lamfac,q3,sigma,
     $     stfac,tau
      COMPLEX(r8), DIMENSION(2,2) :: bmmat,bpmat,m
      COMPLEX(r8), DIMENSION(4,4) :: bl0
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(46x,"bl1"//4(1p,2e11.3,2x)///46x,"bl0"//
     $     4(4(1p,2e11.3,2x)/)/)
 20   FORMAT(47x,"d1"//4(4(1p,2e11.3,2x)/)/)
 30   FORMAT(47x,"m"//2(24x,2(1p,2e11.3,2x)/)/)
 40   FORMAT(45x,"bmmat"//2(24x,2(1p,2e11.3,2x)/)/)
 50   FORMAT(45x,"bpmat"//2(24x,2(1p,2e11.3,2x)/)/)
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      q3=q**3
      kk1=1+kk*q3
      a1=(gg-kk*ee)*q*q/kk1
      lamda=SQRT(q)
      lamdaq=lamda*q
      bb=-.25*lamdaq*(1+gg+kk*(dr-ee))
      cc=.25*kk1*(a1*q*(1-dr/q3)+dr-hh*hh)
      ddsq=bb*bb-cc
      dd=SQRT(ddsq)
c-----------------------------------------------------------------------
c     compute diagonalizing matrix m.
c-----------------------------------------------------------------------
      m(1,1)=-lamdaq+dr/lamdaq
      m(1,2)=hh-dr/lamdaq
      m(2,1)=kk1*(hh+dr/lamdaq)
      m(2,2)=-kk1*(a1/lamda+dr/lamdaq)
c-----------------------------------------------------------------------
c     compute powers.
c-----------------------------------------------------------------------
      pexp(1)= bb-dd
      pexp(2)= bb+dd
      pexp(3)=-bb+dd
      pexp(4)=-bb-dd
c-----------------------------------------------------------------------
c     compute diagonalizing matrices bmmat and d1.
c-----------------------------------------------------------------------
      DO j=1,2
         bmmat(1,j)=m(1,2)
         bmmat(2,j)=2.*pexp(j)-m(1,1)
         DO i=1,neqc
            d1(i,j)=d0inv(i,1)*bmmat(1,j)+d0inv(i,2)*bmmat(2,j)
     $           +lamda*(-d0inv(i,3)*bmmat(1,j)+d0inv(i,4)*bmmat(2,j)
     $           /kk1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute diagonalizing matrix bpmat and d1.
c-----------------------------------------------------------------------
      sigma=bmmat(1,2)/bmmat(2,2)
      tau=bmmat(2,1)/bmmat(1,1)
      stfac=.5/(1-sigma*tau)
      bpmat(1,1)=stfac/bmmat(1,1)
      bpmat(2,2)=stfac/bmmat(2,2)
      bpmat(1,2)=-tau*bpmat(2,2)
      bpmat(2,1)=-sigma*bpmat(1,1)
      DO j=1,2
         DO i=1,neqc
            d1(i,j+2)=(d0inv(i,1)*bpmat(1,j)-d0inv(i,2)*bpmat(2,j)*kk1)
     $           /lamda+d0inv(i,3)*bpmat(1,j)+d0inv(i,4)*bpmat(2,j)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute diagonalizing matrix d1inv.
c-----------------------------------------------------------------------
      d1inv(1,1)= d1(3,3)
      d1inv(2,2)= d1(4,4)
      d1inv(1,2)= d1(4,3)
      d1inv(2,1)= d1(3,4)
      d1inv(3,3)= d1(1,1)
      d1inv(4,4)= d1(2,2)
      d1inv(3,4)= d1(2,1)
      d1inv(4,3)= d1(1,2)
      d1inv(1,3)=-d1(1,3)
      d1inv(2,4)=-d1(2,4)
      d1inv(1,4)=-d1(2,3)
      d1inv(2,3)=-d1(1,4)
      d1inv(3,1)=-d1(3,1)
      d1inv(4,2)=-d1(4,2)
      d1inv(3,2)=-d1(4,1)
      d1inv(4,1)=-d1(3,2)
c-----------------------------------------------------------------------
c     compute transformed coefficient matrices bl0 and bl1.
c-----------------------------------------------------------------------
      DO i=1,neqc
         DO j=1,neqc
            bl0(i,j)=0
         ENDDO
      ENDDO
      DO k=1,neqc
         DO j=1,neqc
            DO i=1,neqc
               bl0(i,j)=bl0(i,j)+d1inv(i,k)*al0(k)*d1(k,j)
            ENDDO
         ENDDO
      ENDDO
      bl1(1)=-lamda
      bl1(2)=-lamda
      bl1(3)=lamda
      bl1(4)=lamda
c-----------------------------------------------------------------------
c     compute lowest-order power series coefficient matrix elements.
c-----------------------------------------------------------------------
      DO j=1,neqc
         DO i=1,neqc
            DO k=1,nxps
               vps(i,j,k)=0
            ENDDO
         ENDDO
         vps(j,j,1)=1
      ENDDO
c-----------------------------------------------------------------------
c     compute higher-order power series coefficient matrix elements.
c-----------------------------------------------------------------------
      j1=1
      j2=2
      j3=3
      j4=4
      lamfac=-.5/lamda
      DO jlam=1,2
         DO j=j1,j2
            DO k=2,nxps
               DO i=j3,j4
                  DO l=1,neqc
                     vps(i,j,k)=vps(i,j,k)+bl0(i,l)*vps(l,j,k-1)
                  ENDDO
                  vps(i,j,k)=vps(i,j,k)+(2.*k-4.-pexp(j))*vps(i,j,k-1)
                  vps(i,j,k)=lamfac*vps(i,j,k)
               ENDDO
               DO i=j1,j2
                  DO l=j3,j4
                     vps(i,j,k)=vps(i,j,k)+bl0(i,l)*vps(l,j,k)
                  ENDDO
                  vps(i,j,k)=vps(i,j,k)/(pexp(j)-pexp(i)-2.*k+2.)
               ENDDO
            ENDDO
         ENDDO
         j1=3
         j2=4
         j3=1
         j4=2
         lamfac=-lamfac
      ENDDO
c-----------------------------------------------------------------------
c     compute tmax.
c-----------------------------------------------------------------------
      IF(ntmax == 1)THEN
         tmax=MIN(SQRT((MAX(.5d0,p1)**2-LOG(rtol))/ABS(lamda)),6/ABS(q))
      ELSE IF(ntmax == 2)THEN
         tmax=SQRT(vps(2,2,2))
      ELSE
         tmax=SQRT((MAX(.5d0,p1)**2-LOG(rtol))/ABS(lamda))
      ENDIF
      tmax=tmax*fmax
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      IF(out2)THEN
         WRITE(out_unit,10)bl1,((bl0(i,j),j=1,neqc),i=1,neqc)
         WRITE(out_unit,20)((d1(i,j),j=1,neqc),i=1,neqc)
         WRITE(out_unit,30)((m(i,j),i=1,2),j=1,2)
         WRITE(out_unit,40)((bmmat(i,j),i=1,2),j=1,2)
         WRITE(out_unit,50)((bpmat(i,j),i=1,2),j=1,2)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_infinity
c-----------------------------------------------------------------------
c     subprogram 7. deltar_ratio.
c     evaluates asymptotic ratios.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE deltar_ratio(c0,deltar)

      COMPLEX(r8), DIMENSION(4,2), INTENT(IN) :: c0
      COMPLEX(r8), DIMENSION(2), INTENT(OUT) :: deltar

      REAL(r8) :: phase,gamfac
      COMPLEX(r8) :: dfac
c-----------------------------------------------------------------------
c     compute asymptotic matching ratios deltar.
c-----------------------------------------------------------------------
      dfac=-d0(1,3)/d0(1,1)
      phase=pplus*pi/2
      gamfac=-pi/(SIN(2.*phase)*gamma(1+pplus)**2)
      deltar(1)=dfac*gamfac/tan(phase)*(c0(4,1)*c0(3,2)-c0(4,2)*c0(3,1))
     $     /(c0(4,1)*c0(1,2)-c0(4,2)*c0(1,1))
      deltar(2)=dfac*gamfac*tan(phase)*(c0(2,1)*c0(3,2)-c0(2,2)*c0(3,1))
     $     /(c0(2,1)*c0(1,2)-c0(2,2)*c0(1,1))
c-----------------------------------------------------------------------
c     compute determinant ratios.
c-----------------------------------------------------------------------
      dtr(1)=(c0(1,1)*c0(2,2)-c0(1,2)*c0(2,1))
     $     /(c0(1,1)*c0(2,2)+c0(1,2)*c0(2,1))
      dtr(2)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))
     $     /(c0(1,1)*c0(3,2)+c0(1,2)*c0(3,1))
      dtr(3)=(c0(1,1)*c0(4,2)-c0(1,2)*c0(4,1))
     $     /(c0(1,1)*c0(4,2)+c0(1,2)*c0(4,1))
      dtr(4)=(c0(2,1)*c0(3,2)-c0(2,2)*c0(3,1))
     $     /(c0(2,1)*c0(3,2)+c0(2,2)*c0(3,1))
      dtr(5)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))
     $     /(c0(2,1)*c0(4,2)+c0(2,2)*c0(4,1))
      dtr(6)=(c0(3,1)*c0(4,2)-c0(3,2)*c0(4,1))
     $     /(c0(3,1)*c0(4,2)+c0(3,2)*c0(4,1))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deltar_ratio
      END MODULE deltar_mod
