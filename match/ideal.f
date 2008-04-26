c-----------------------------------------------------------------------
c     file ideal.f
c     reads data from dcon, computes and prints ideal eigenfunction.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ideal_mod.
c     1. ideal_read.
c     2. ideal_transform.
c     3. ideal_build.
c     4. ideal_write.
c     5. ideal_contour.
c     6. ideal_chord.
c     7. ideal_eval.
c-----------------------------------------------------------------------
c     subprogram 0. ideal_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ideal_mod
      USE match_mod
      USE bicube_mod
      IMPLICIT NONE

      LOGICAL :: ripple_flag=.FALSE.
      INTEGER :: mpsi,mtheta,mripple
      REAL(r8) :: ro,zo
      TYPE(spline_type) :: sq
      TYPE(bicube_type) :: rzphi

      CHARACTER(8) :: z_level_type="axis"
      REAL(r8) :: phi_level=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ideal_read.
c     reads Euler-Lagrange solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_read(filename)

      CHARACTER(*), INTENT(IN) :: filename

      CHARACTER(128) :: message
      INTEGER :: data_type,ifix,ios,msol,istep,ising

      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u
c-----------------------------------------------------------------------
c     open data file and read header.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,filename,"OLD","REWIND","none")
      READ(in_unit)mlow,mhigh,nn,mpsi,mtheta,ro,zo
      CALL spline_alloc(sq,mpsi,4)
      CALL bicube_alloc(rzphi,mpsi,mtheta,4)
      rzphi%periodic(2)=.TRUE.
      READ(in_unit)sq%xs,sq%fs,sq%fs1
      READ(in_unit)rzphi%xs,rzphi%ys,
     $        rzphi%fs,rzphi%fsx,rzphi%fsy,rzphi%fsxy,
     $        rzphi%x0,rzphi%y0,rzphi%xpower,rzphi%ypower
      mstep=-1
      mfix=0
      msing=0
      WRITE(*,*)"Count solutions."
c-----------------------------------------------------------------------
c     count.
c-----------------------------------------------------------------------
      DO
         READ(UNIT=in_unit,IOSTAT=ios)data_type
         IF(ios /= 0)EXIT
         SELECT CASE(data_type)
         CASE(1)
            mstep=mstep+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(2)
            mfix=mfix+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(3)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE(4)
            msing=msing+1
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
            READ(UNIT=in_unit)
         CASE DEFAULT
            WRITE(message,'(a,i1,a,i4)')"Cannot regognize data_type = ",
     $           data_type,", at istep = ",istep
            CALL program_stop(message)
         END SELECT
      ENDDO
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     allocate arrays, prepare to read data.
c-----------------------------------------------------------------------
      mpert=mhigh-mlow+1
      WRITE(*,*)"mlow = ",mlow,", mhigh = ",mhigh," mpert = ",mpert
      WRITE(*,*)"mstep = ",mstep,", mfix = ",mfix,", msing = ",msing
      ALLOCATE(psifac(0:mstep),rho(0:mstep),q(0:mstep),
     $     soltype(0:mstep),v(mpert,2,0:mstep),singtype(msing))
      ALLOCATE(fixstep(0:mfix+1),fixtype(0:mfix),sing_flag(mfix))
      ALLOCATE(et(mpert),wt(mpert,mpert))
      fixstep(0)=0
      fixstep(mfix+1)=mstep
      CALL bin_open(in_unit,filename,"OLD","REWIND","none")
      READ(in_unit)mlow,mhigh,nn
      istep=-1
      ifix=0
      ising=0
      WRITE(*,*)"Read solutions."
c-----------------------------------------------------------------------
c     read data.
c-----------------------------------------------------------------------
      DO
         READ(UNIT=in_unit,IOSTAT=ios)data_type
         IF(ios /= 0)EXIT
         SELECT CASE(data_type)
         CASE(1)
            istep=istep+1
            READ(UNIT=in_unit)psifac(istep),q(istep),soltype(istep)%msol
            ALLOCATE(soltype(istep)%u(mpert,soltype(istep)%msol,2))
            u => soltype(istep)%u
            READ(UNIT=in_unit)soltype(istep)%u
         CASE(2)
            ifix=ifix+1
            fixstep(ifix)=istep
            READ(UNIT=in_unit)sing_flag(ifix),fixtype(ifix)%msol
            ALLOCATE(fixtype(ifix)%fixfac
     $           (fixtype(ifix)%msol,fixtype(ifix)%msol))
            ALLOCATE(fixtype(ifix)%index(fixtype(ifix)%msol))
            READ(UNIT=in_unit)fixtype(ifix)%fixfac,fixtype(ifix)%index
         CASE(3)
            READ(UNIT=in_unit)et
            READ(UNIT=in_unit)wt
         CASE(4)
            ising=ising+1
            singtype(ising)%jfix=ifix
            singtype(ising)%jpert=INT(nn*q(istep)+.5)-mlow+1
            READ(UNIT=in_unit)singtype(ising)%psifac,
     $           singtype(ising)%q,singtype(ising)%q1
            READ(UNIT=in_unit)msol
            singtype(ising)%msol_l=msol
            ALLOCATE(singtype(ising)%ca_l(mpert,msol,2))
            READ(UNIT=in_unit)singtype(ising)%ca_l
            READ(UNIT=in_unit)msol
            singtype(ising)%msol_r=msol
            ALLOCATE(singtype(ising)%ca_r(mpert,msol,2))
            READ(UNIT=in_unit)singtype(ising)%ca_r
            READ(UNIT=in_unit)
     $           singtype(ising)%restype%e,
     $           singtype(ising)%restype%f,
     $           singtype(ising)%restype%h,
     $           singtype(ising)%restype%m,
     $           singtype(ising)%restype%g,
     $           singtype(ising)%restype%k,
     $           singtype(ising)%restype%eta,
     $           singtype(ising)%restype%rho,
     $           singtype(ising)%restype%taua,
     $           singtype(ising)%restype%taur
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     modify Lundquist numbers.
c-----------------------------------------------------------------------
      DO ising=1,msing
         singtype(ising)%restype%taur=singtype(ising)%restype%taur*sfac0
      ENDDO
c-----------------------------------------------------------------------
c     close data file.
c-----------------------------------------------------------------------
      CALL bin_close(in_unit)
      rho=SQRT(psifac)
c-----------------------------------------------------------------------
c     allocate matching matrix.
c-----------------------------------------------------------------------
      IF(singtype(msing)%msol_r > mpert .AND. res_flag)THEN
         mmatch=mpert
         DO ising=1,msing
            mmatch=mmatch+singtype(ising)%msol_r
         ENDDO
         ALLOCATE(match(mmatch,mmatch),deltap1(msing),deltap2(msing))
         WRITE(*,'(1x,a,i4)')"mmatch = ",mmatch
      ELSE
         res_flag=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_read
c-----------------------------------------------------------------------
c     subprogram 2. ideal_transform.
c     build fixup matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_transform

      INTEGER :: ifix,isol,jsol,ksol
      LOGICAL, DIMENSION(mpert) :: mask
      COMPLEX(r8), DIMENSION(mpert,mpert) :: ident,temp
c-----------------------------------------------------------------------
c     create identity matrix.
c-----------------------------------------------------------------------
      WRITE(*,*)"Compute transformation matrices."
      ident=0
      DO isol=1,mpert
         ident(isol,isol)=1
      ENDDO
c-----------------------------------------------------------------------
c     compute gaussian reduction matrices.
c-----------------------------------------------------------------------
      DO ifix=1,mfix
         ALLOCATE(fixtype(ifix)%gauss(mpert,mpert))
         fixtype(ifix)%gauss=ident
         mask=.TRUE.
         DO isol=1,mpert
            ksol=fixtype(ifix)%index(isol)
            mask(ksol)=.FALSE.
            temp=ident
            DO jsol=1,mpert
               IF(mask(jsol))
     $              temp(ksol,jsol)=fixtype(ifix)%fixfac(ksol,jsol)
            ENDDO
            fixtype(ifix)%gauss=MATMUL(fixtype(ifix)%gauss,temp)
         ENDDO
         IF(sing_flag(ifix))
     $        fixtype(ifix)%gauss(:,fixtype(ifix)%index(1))=0
      ENDDO
c-----------------------------------------------------------------------
c     concatenate gaussian reduction matrices.
c-----------------------------------------------------------------------
      ALLOCATE(fixtype(mfix)%transform(mpert,mpert))
      fixtype(mfix)%transform=ident
      DO ifix=mfix-1,0,-1
         ALLOCATE(fixtype(ifix)%transform(mpert,mpert))
         fixtype(ifix)%transform
     $        =MATMUL(fixtype(ifix+1)%gauss,
     $        fixtype(ifix+1)%transform)
      ENDDO
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DO ifix=1,mfix
         DEALLOCATE(fixtype(ifix)%gauss)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_transform
c-----------------------------------------------------------------------
c     subprogram 3. ideal_build.
c     builds ideal Euler-Lagrange solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_build

      INTEGER :: istep,ifix,jfix,kfix,ieq,info
      INTEGER, DIMENSION(mpert) :: ipiv
      COMPLEX(r8), DIMENSION(mpert) :: uedge,temp1
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp2
c-----------------------------------------------------------------------
c     construct uedge.
c-----------------------------------------------------------------------
      WRITE(*,*)"Construct ideal eigenfunctions."
      IF(ripple_flag)THEN
         uedge=0
         uedge(mripple-mlow+1)=1
      ELSE
         uedge=wt(:,1)
      ENDIF
      temp2=soltype(mstep)%u(:,1:mpert,1)
      CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
      CALL zgetrs('N',mpert,1,temp2,mpert,ipiv,uedge,mpert,info)
c-----------------------------------------------------------------------
c     construct eigenfunctions.
c-----------------------------------------------------------------------
      jfix=0
      DO ifix=0,mfix
         temp1=MATMUL(fixtype(ifix)%transform,uedge)
         kfix=fixstep(ifix+1)
         DO ieq=1,2
            DO istep=jfix,kfix
               v(:,ieq,istep)
     $              =MATMUL(soltype(istep)%u(:,1:mpert,ieq),temp1)
            ENDDO
         ENDDO
         jfix=kfix+1
      ENDDO
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DO ifix=0,mfix
         DEALLOCATE(fixtype(ifix)%transform)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_build
c-----------------------------------------------------------------------
c     subprogram 4. ideal_write.
c     writes Euler-Lagrange solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_write

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ipert,istep,m
      REAL(r8), DIMENSION(mpert) :: singfac
      COMPLEX(r8), DIMENSION(mpert,0:mstep) :: b
c-----------------------------------------------------------------------
c     compute normal perturbed magnetic field.
c-----------------------------------------------------------------------
      DO istep=0,mstep
         CALL spline_eval(sq,psifac(istep),0)
         singfac=(/(m,m=mlow,mhigh)/)-nn*sq%f(4)
         b(:,istep)=ifac*singfac*v(:,1,istep)
      ENDDO
c-----------------------------------------------------------------------
c     write binary output for graphs.
c-----------------------------------------------------------------------
      WRITE(*,*)"Write binary output for graphs."
      CALL bin_open(bin_unit,"solutions.bin","UNKNOWN","REWIND","none")
      DO ipert=1,mpert
         DO istep=0,mstep
            WRITE(bin_unit)REAL(psifac(istep),4),REAL(rho(istep),4),
     $           REAL(q(istep),4),
     $           REAL(REAL(v(ipert,1,istep)),4),
     $           REAL(AIMAG(v(ipert,1,istep)),4),
     $           REAL(REAL(b(ipert,istep)),4),
     $           REAL(AIMAG(b(ipert,istep)),4)
         ENDDO
         WRITE(bin_unit)
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_write
c-----------------------------------------------------------------------
c     subprogram 5. ideal_contour.
c     writes contour plots.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_contour

      INTEGER :: istep,itheta,nqty=4
      REAL(r8), DIMENSION(0:mstep,0:mtheta,2) :: rz
      COMPLEX(r8), DIMENSION(0:mstep,0:mtheta) :: xi_psi,b_psi,
     $     xi_norm,b_norm

      LOGICAL, PARAMETER :: diagnose=.FALSE.
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"ith",2x,"istep",3x,"theta",7x,"psi",9x,"r",
     $     10x,"z",8x,"xi_psi",5x,"b_psi",5x,"xi_norm",5x,"b_norm"/)
 20   FORMAT(2i6,1p,8e11.3)
 30   FORMAT(/4x,"name",5x,"min",8x,"max"/)
 40   FORMAT(a8,1p,2e11.3)
c-----------------------------------------------------------------------
c     compute contour values.
c-----------------------------------------------------------------------
      WRITE(*,*)"Compute values for contour plots"
      DO istep=0,mstep
         DO itheta=0,mtheta
            CALL ideal_eval(psifac(istep),rzphi%ys(itheta),
     $           v(:,1,istep),rz(istep,itheta,1),rz(istep,itheta,2),
     $           xi_psi(istep,itheta),b_psi(istep,itheta),
     $           xi_norm(istep,itheta),b_norm(istep,itheta))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write contour values.
c-----------------------------------------------------------------------
      WRITE(*,*)"Write contour plots"
      CALL bin_open(bin_unit,"contour.bin","UNKNOWN","REWIND","none")
      WRITE(bin_unit)1,0,nqty
      WRITE(bin_unit)mstep,mtheta
      WRITE(bin_unit)REAL(rz,4)
      WRITE(bin_unit)REAL(xi_psi,4)
      WRITE(bin_unit)REAL(b_psi,4)
      WRITE(bin_unit)REAL(xi_norm,4)
      WRITE(bin_unit)REAL(b_norm,4)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     diagnose contour values.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(*,*)"Diagnose values for contour plots"
         CALL ascii_open(98,"debug.out","UNKNOWN")
         CALL bin_open(99,"debug.bin","UNKNOWN","REWIND","none")
         WRITE(98,10)
         DO itheta=0,mtheta
            DO istep=0,mstep
               WRITE(98,20)itheta,istep,
     $              rzphi%ys(itheta),psifac(istep),
     $              rz(istep,itheta,:),
     $              REAL(xi_psi(istep,itheta)),
     $              REAL(b_psi(istep,itheta)),
     $              REAL(xi_norm(istep,itheta)),
     $              REAL(b_norm(istep,itheta))
               WRITE(99)REAL(psifac(istep),4),
     $              REAL(rzphi%ys(itheta),4),
     $              REAL(rz(istep,itheta,:),4),
     $              REAL(xi_psi(istep,itheta),4),
     $              REAL(b_psi(istep,itheta),4),
     $              REAL(xi_norm(istep,itheta),4),
     $              REAL(b_norm(istep,itheta),4)
            ENDDO
            WRITE(98,10)
            WRITE(99)
         ENDDO
         CALL bin_close(99)
c-----------------------------------------------------------------------
c     diagnose extrema.
c-----------------------------------------------------------------------
         WRITE(98,30)
         WRITE(98,40)"xi_psi",
     $        MINVAL(REAL(xi_psi)),
     $        MAXVAL(REAL(xi_psi))
         WRITE(98,40)"b_psi",
     $        MINVAL(REAL(b_psi)),
     $        MAXVAL(REAL(b_psi))
         WRITE(98,40)"xi_norm",
     $        MINVAL(REAL(xi_norm)),
     $        MAXVAL(REAL(xi_norm))
         WRITE(98,40)"b_norm",
     $        MINVAL(REAL(b_norm)),
     $        MAXVAL(REAL(b_norm))
         WRITE(98,30)
         CALL ascii_close(98)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_contour
c-----------------------------------------------------------------------
c     subprogram 6. ideal_chord.
c     diagnoses solutions along horizontal chord.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_chord

      INTEGER :: istep,jstep,itheta,iside,instep,it,itmax=100
      INTEGER, DIMENSION(mstep,2) :: jtheta
      REAL(r8) :: psi,theta,theta0,theta1,zlevel,z,zz0,zz1
      REAL(r8), DIMENSION(0:mtheta) :: zz
      REAL(r8), DIMENSION(mstep,2) :: r,th0,z0,q
      COMPLEX(r8), DIMENSION(mstep,2) :: xi_psi,b_psi,xi_norm,b_norm 
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"js",5x,"psi",8x,"th0",9x,"r",10x,"q",8x,
     $     "xi_psi",5x,"b_psi",5x,"xi_norm",5x,"b_norm"/)
 20   FORMAT(i6,1p,8e11.3)
c-----------------------------------------------------------------------
c     choose z position of horizontal chord.
c-----------------------------------------------------------------------
      WRITE(*,*)"Write chord diagnostic."
      SELECT CASE(z_level_type)
      CASE("axis")
         zlevel=zo
      CASE("midplane")
         zlevel=0
      CASE DEFAULT
         CALL program_stop("Cannot recognize z_level_type = "
     $        //TRIM(z_level_type))
      END SELECT
c-----------------------------------------------------------------------
c     find z positions on flux surface.
c-----------------------------------------------------------------------
      DO istep=mstep,1,-1
         psi=psifac(istep)
         DO itheta=0,mtheta
            theta=rzphi%ys(itheta)
            CALL bicube_eval(rzphi,psi,theta,0)
            zz(itheta)=zo-zlevel
     $           +SQRT(rzphi%f(1))*SIN(twopi*(theta+rzphi%f(2)))
         ENDDO         
c-----------------------------------------------------------------------
c     find theta indices intersecting chord.
c-----------------------------------------------------------------------
         iside=0
         DO itheta=0,mtheta-1
            IF(zz(itheta)*zz(itheta+1) <= 0)THEN
               iside=iside+1
               jtheta(istep,iside)=itheta
            ENDIF
         ENDDO
         IF(iside < 2)EXIT
c-----------------------------------------------------------------------
c     refine theta position by binary search.
c-----------------------------------------------------------------------
         DO iside=1,2
            theta0=rzphi%ys(jtheta(istep,iside))
            theta1=rzphi%ys(jtheta(istep,iside)+1)
            zz0=zz(jtheta(istep,iside))
            zz1=zz(jtheta(istep,iside)+1)
            it=0
            DO
               it=it+1
               theta=(theta0+theta1)/2
               CALL bicube_eval(rzphi,psi,theta,0)
               z=zo-zlevel
     $              +SQRT(rzphi%f(1))*SIN(twopi*(theta+rzphi%f(2)))
               IF(z*zz0 > 0)THEN
                  theta0=theta
                  theta=(theta+theta1)/2
               ELSE
                  theta1=theta
                  theta=(theta+theta0)/2
               ENDIF
               IF(ABS(z) <= 1e-8)EXIT
               IF(it > itmax)
     $              CALL program_stop("Ideal_chord can't find root")
            ENDDO
            th0(istep,iside)=theta
            z0(istep,iside)=z
c-----------------------------------------------------------------------
c     evaluate fourier series.
c-----------------------------------------------------------------------
            CALL ideal_eval(psi,theta,v(:,1,istep),r(istep,iside),z,
     $           xi_psi(istep,iside),b_psi(istep,iside),
     $           xi_norm(istep,iside),b_norm(istep,iside))
            q(istep,iside)=sq%f(4)
         ENDDO
      ENDDO
      instep=istep+1
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      CALL ascii_open(debug_unit,"chord.out","UNKNOWN")
      CALL bin_open(bin_unit,"chord.bin","UNKNOWN","REWIND","none")
      WRITE(debug_unit,'(1x,3a,1p,e10.3)')
     $     "zlevel_type = ",TRIM(z_level_type),", zlevel = ",zlevel
      WRITE(debug_unit,10)
      jstep=0
      DO istep=mstep,instep,-1
         WRITE(bin_unit)REAL(r(istep,1),4),REAL(q(istep,1),4),
     $        REAL(xi_psi(istep,1),4),REAL(b_psi(istep,1),4),
     $        REAL(xi_norm(istep,1),4),REAL(b_norm(istep,1),4)
         WRITE(debug_unit,20)jstep,psifac(istep),th0(istep,1),
     $        r(istep,1),q(istep,1),
     $        REAL(xi_psi(istep,1)),REAL(b_psi(istep,1)),
     $        REAL(xi_norm(istep,1)),REAL(b_norm(istep,1))
         jstep=jstep+1
      ENDDO
      WRITE(debug_unit,10)
      DO istep=instep,mstep
         WRITE(bin_unit)REAL(r(istep,2),4),REAL(q(istep,2),4),
     $        REAL(xi_psi(istep,2),4),REAL(b_psi(istep,2),4),
     $        REAL(xi_norm(istep,2),4),REAL(b_norm(istep,2),4)
         WRITE(debug_unit,20)jstep,psifac(istep),th0(istep,2),
     $        r(istep,2),q(istep,2),
     $        REAL(xi_psi(istep,2)),REAL(b_psi(istep,2)),
     $        REAL(xi_norm(istep,2)),REAL(b_norm(istep,2))
         jstep=jstep+1
      ENDDO
      WRITE(debug_unit,10)
      CALL bin_close(bin_unit) 
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_chord
c-----------------------------------------------------------------------
c     subprogram 7. ideal_eval.
c     diagnoses solutions along horizontal eval.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ideal_eval(psi,theta,xi_vec,r,z,
     $     xi_psi,b_psi,xi_norm,b_norm)

      REAL(r8), INTENT(IN) :: psi,theta
      COMPLEX(r8), DIMENSION(mpert), INTENT(IN) :: xi_vec
      REAL(r8), INTENT(OUT) :: r,z
      COMPLEX(r8), INTENT(OUT) :: xi_psi,b_psi,xi_norm,b_norm

      INTEGER :: ipert
      REAL(r8) :: r2,deta,dphi,rfac,eta,jac,v21,v22,dpsisq,norm,singfac
      COMPLEX(r8) :: expfac,expfac1
c-----------------------------------------------------------------------
c     evaluate radial position.
c-----------------------------------------------------------------------
      CALL bicube_eval(rzphi,psi,theta,1)
      r2=rzphi%f(1)
      deta=rzphi%f(2)
      dphi=rzphi%f(3)
      rfac=SQRT(r2)
      eta=theta+deta
      r=ro+rfac*COS(twopi*eta)
      z=zo+rfac*SIN(twopi*eta)
c-----------------------------------------------------------------------
c     evaluate normalizing factors.
c-----------------------------------------------------------------------
      jac=rzphi%f(4)
      v21=rzphi%fy(1)/(2*rfac*jac)
      v22=(1+rzphi%fy(2))*twopi*rfac/jac
      dpsisq=(twopi*r)**2*(v21**2+v22**2)
c-----------------------------------------------------------------------
c     evaluate fourier series.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      expfac=EXP(ifac*(mlow*twopi*theta-nn*(phi_level-dphi)))
      expfac1=EXP(ifac*twopi*theta)
      singfac=mlow-nn*sq%f(4)
      xi_psi=0
      b_psi=0
      DO ipert=1,mpert
         xi_psi=xi_psi+xi_vec(ipert)*expfac
         b_psi=b_psi+xi_vec(ipert)*expfac*singfac*ifac
         singfac=singfac+1
         expfac=expfac*expfac1
      ENDDO
c-----------------------------------------------------------------------
c     compute xi_norm and b_norm.
c-----------------------------------------------------------------------
      norm=1/SQRT(dpsisq)
      xi_norm=xi_psi*norm
      b_psi=b_psi/jac
      b_norm=b_psi*norm
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ideal_eval
      END MODULE ideal_mod
