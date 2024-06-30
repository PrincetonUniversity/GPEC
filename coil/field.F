c-----------------------------------------------------------------------
c     file field.f.
c     calculate 3d field from coils
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. field_mod.
c     1. field_bs_psi.
c     2. field_bs_rz.
c-----------------------------------------------------------------------
c     subprogram 0. field_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE field_mod
      USE coil_mod
      USE omp_lib

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. field_bs_psi.
c
c     *DESCRIPTION:
c        Calculate 3d field spectrum on a psi using biot-savart law.
c
c     *ARGUMENTS:
c         psi : real
c            Flux surface on which biat-savart field is decomposed
c            magnetic fourier components
c         bmn : real array
c            Fourier spectrum of normal field on psi. Coordinates are jac_type?
c         wegt : integer
c            Surface area weighting flag (like jsurf).
c         op_start : integer (optional)
c            Field from coils op_start to op_stop is calculated. Default is 1 to coil_num.
c         op_stop : integer (optional)
c            Field from coils op_start to op_stop is calculated. Default is 1 to coil_num.
c         op_verbose : logical (optional)
c            Print information about parallelization threads.
c
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_bs_psi(psi,bmn,wegt,op_start,op_stop, op_verbose)

      REAL(r8), INTENT(IN) :: psi
      INTEGER, INTENT(IN) :: wegt
      INTEGER, INTENT(IN), OPTIONAL :: op_start,op_stop
      LOGICAL, INTENT(IN), OPTIONAL :: op_verbose
      COMPLEX(r8), DIMENSION(cmpert) :: bmn

      LOGICAL :: verbose = .FALSE.
      INTEGER :: num_threads = 1


      INTEGER :: i,j,k,iseg,itheta,izeta,ipert,nseg,
     $    istart,istop,spline_itheta,spline_ipsi

      REAL(r8) :: rfac,eta,phi,jac,delpsi,rr,zz,rx,ry,rz,dl,dbx,dby,dbz
      REAL(r8) :: cosang,sinang,w11,w12,area
      REAL(r8), DIMENSION(0:cmtheta) :: ctheta
      REAL(r8), DIMENSION(0:cmzeta) :: czeta
      REAL(r8), DIMENSION(cmtheta,cmzeta) :: xobs,yobs,zobs,
     $     bx,by,bz,br,bp,bn,tbx,tby,tbz,tbr,tbp,tbn
      REAL(r8), DIMENSION(:,:,:), POINTER :: xa,ya,za,
     $     dx,dy,dz

      REAL(r8), DIMENSION(4) :: crzphi_f, crzphi_fx, crzphi_fy

      REAL(r8), DIMENSION(cmtheta) :: aspl_f_arr

      TYPE(spline_type) :: aspl
      TYPE(cspline_type) :: bntspl,bnzspl


      ! Optional verbosity tells user how many threads we are using
      IF(PRESENT(op_verbose)) THEN
        verbose = op_verbose
      ELSE
        verbose = .FALSE.
      ENDIF

      istart = 1
      istop = coil_num
      IF(PRESENT(op_start)) istart = op_start
      IF(PRESENT(op_stop)) istop = op_stop

      tbx=0
      tby=0
      tbz=0
      tbr=0
      tbp=0
      tbn=0

      ctheta=(/(itheta,itheta=0,cmtheta)/)/REAL(cmtheta,r8)
      czeta=(/(izeta,izeta=0,cmzeta)/)/REAL(cmzeta,r8)

      CALL spline_alloc(aspl,cmtheta,1)
      aspl%xs=ctheta
      area=1.0
      
c-----------------------------------------------------------------------
c     for each coil:
c-----------------------------------------------------------------------

      DO i=istart,istop
         nseg=coil(i)%nsec-1

         ALLOCATE(xa(coil(i)%ncoil,coil(i)%s,nseg),
     $        ya(coil(i)%ncoil,coil(i)%s,nseg),
     $        za(coil(i)%ncoil,coil(i)%s,nseg),
     $        dx(coil(i)%ncoil,coil(i)%s,nseg),
     $        dy(coil(i)%ncoil,coil(i)%s,nseg),
     $        dz(coil(i)%ncoil,coil(i)%s,nseg))

      IF(verbose) WRITE(*,'(1x,a,i6,a,i3,a,i3,a)'),"  > Mapping "//
     $      TRIM(coil(i)%coil_name)//" coil's ",
     $      coil(i)%ncoil*coil(i)%s*nseg,
     $      " filaments to a ",cmtheta,"x",cmzeta," plasma grid"

!$OMP PARALLEL DO PRIVATE(j,k,iseg),COLLAPSE(3)
         DO j=1,coil(i)%ncoil
            DO k=1,coil(i)%s
               DO iseg=1,nseg
                  xa(j,k,iseg)=(coil(i)%x(j,k,iseg+1)+
     $                 coil(i)%x(j,k,iseg))/2.0
                  ya(j,k,iseg)=(coil(i)%y(j,k,iseg+1)+
     $                 coil(i)%y(j,k,iseg))/2.0      
                  za(j,k,iseg)=(coil(i)%z(j,k,iseg+1)+
     $                 coil(i)%z(j,k,iseg))/2.0
                  
                  dx(j,k,iseg)=coil(i)%x(j,k,iseg+1)-coil(i)%x(j,k,iseg)
                  dy(j,k,iseg)=coil(i)%y(j,k,iseg+1)-coil(i)%y(j,k,iseg)
                  dz(j,k,iseg)=coil(i)%z(j,k,iseg+1)-coil(i)%z(j,k,iseg)

#ifdef _OPENMP
                  IF(j==1 .AND. k==1 .AND. iseg==1)THEN
                    num_threads = OMP_GET_NUM_THREADS()
                    IF(verbose) WRITE(*,'(1x,a,i4,a)'),"  >> "//
     $                "Running with ",num_threads," OMP threads"
                  ENDIF
#endif
               ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL PRIVATE(izeta,phi,dbx,dby,dbz,j,k,iseg,rx,ry,rz,dl)
!$OMP& PRIVATE(cosang,sinang)
!$OMP& FIRSTPRIVATE(w11,w12)
!$OMP& FIRSTPRIVATE(crzphi_f,crzphi_fx,crzphi_fy)
!$OMP& FIRSTPRIVATE(spline_itheta,spline_ipsi)
!$OMP& PRIVATE(itheta,rfac,eta,rr,zz,jac,delpsi)

!$OMP DO
         ! Loop over every theta, zeta point on the flux surface
         DO itheta=1,cmtheta

            CALL bicube_eval_external(crzphi,psi,ctheta(itheta),1,
     $                                spline_ipsi,spline_itheta,
     $                                crzphi_f,crzphi_fx,crzphi_fy)

            rfac=SQRT(crzphi_f(1))
            eta=twopi*(ctheta(itheta)+crzphi_f(2))
            rr=cro+rfac*COS(eta)
            zz=czo+rfac*SIN(eta)
            jac=crzphi_f(4)
            w11=(1+crzphi_fy(2))*twopi**2*rfac*rr/jac
            w12=-crzphi_fy(1)*pi*rr/(rfac*jac)
            delpsi=SQRT(w11**2+w12**2)
            IF(wegt/=0)THEN
               ! aspl%fs(itheta,1)=jac*delpsi
               aspl_f_arr(itheta)=jac*delpsi
            ELSE
               ! aspl%fs(itheta,1)=1.0
               aspl_f_arr(itheta)=1.0
            ENDIF
            w11=w11/delpsi
            w12=w12/delpsi

            DO izeta=1,cmzeta
               phi=-helicity*(twopi*czeta(izeta)+crzphi_f(3))
               xobs(itheta,izeta)=rr*COS(phi)
               yobs(itheta,izeta)=rr*SIN(phi)
               zobs(itheta,izeta)=zz
               dbx=0
               dby=0
               dbz=0
               
               DO j=1,coil(i)%ncoil
                  DO k=1,coil(i)%s
                     DO iseg=1,nseg
                        rx=xobs(itheta,izeta)-xa(j,k,iseg)
                        ry=yobs(itheta,izeta)-ya(j,k,iseg)
                        rz=zobs(itheta,izeta)-za(j,k,iseg)
                        
                        dl=SQRT(rx**2+ry**2+rz**2)
                        dbx=dbx+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                       dl**3*(dy(j,k,iseg)*rz-dz(j,k,iseg)*ry)
                        dby=dby+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                       dl**3*(dz(j,k,iseg)*rx-dx(j,k,iseg)*rz)
                        dbz=dbz+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                       dl**3*(dx(j,k,iseg)*ry-dy(j,k,iseg)*rx)               
                        
                     ENDDO
                  ENDDO
               ENDDO
               
               bx(itheta,izeta)=dbx
               by(itheta,izeta)=dby
               bz(itheta,izeta)=dbz
               cosang=xobs(itheta,izeta)/rr
               sinang=yobs(itheta,izeta)/rr
               br(itheta,izeta)=bx(itheta,izeta)*cosang+
     $              by(itheta,izeta)*sinang
               bp(itheta,izeta)=by(itheta,izeta)*cosang-
     $              bx(itheta,izeta)*sinang
               ! bn(itheta,izeta)=aspl%fs(itheta,1)*
               bn(itheta,izeta)=aspl_f_arr(itheta)*
     $              (br(itheta,izeta)*(w11*COS(eta)-w12*SIN(eta))+
     $              bz(itheta,izeta)*(w11*SIN(eta)+w12*COS(eta)))
            ENDDO


         ! OMP call


         ENDDO

!$OMP END DO
!$OMP END PARALLEL

         DEALLOCATE(xa,ya,za,dx,dy,dz)

         tbx=tbx+bx
         tby=tby+by
         tbz=tbz+bz
         tbr=tbr+br
         tbp=tbp+bp
         tbn=tbn+bn
      ENDDO

      aspl%fs(:,1)=aspl_f_arr
         
      CALL cspline_alloc(bnzspl,cmzeta,1)
      bnzspl%xs=czeta
      CALL cspline_alloc(bntspl,cmtheta,1)
      bntspl%xs=ctheta
      aspl%fs(0,1)=aspl%fs(cmtheta,1)

      IF(wegt==2)THEN
         CALL spline_fit(aspl,"periodic")
         CALL spline_int(aspl)
         area=aspl%fsi(cmtheta,1)
      ENDIF
            
      DO ipert=1,cmpert
         DO izeta=1,cmzeta
            DO itheta=1,cmtheta
               bntspl%fs(itheta,1)=tbn(itheta,izeta)*
     $              EXP(-twopi*ifac*cmfac(ipert)*ctheta(itheta))
            ENDDO
            bntspl%fs(0,1)=bntspl%fs(cmtheta,1)
            CALL cspline_fit(bntspl,"periodic")
            CALL cspline_int(bntspl)
            bnzspl%fs(izeta,1)=bntspl%fsi(cmtheta,1)*
     $           EXP(twopi*ifac*cnn*czeta(izeta))
         ENDDO
         bnzspl%fs(0,1)=bnzspl%fs(cmzeta,1)
         CALL cspline_fit(bnzspl,"periodic")
         CALL cspline_int(bnzspl)
         bmn(ipert)=2.0*bnzspl%fsi(cmzeta,1)/area

      ENDDO

      CALL spline_dealloc(aspl)
      CALL cspline_dealloc(bnzspl)
      CALL cspline_dealloc(bntspl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_bs_psi
c-----------------------------------------------------------------------
c     subprogram 2. field_bs_rzphi.
c     calculate rzphi 3d field given a grid using biot-savart law.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE field_bs_rzphi(nr,nz,np,gdr,gdz,cbr,cbz,cbp)

      INTEGER, INTENT(IN) :: nr,nz,np
      REAL(r8), DIMENSION(0:nr,0:nz) :: gdr,gdz
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: cbr,cbz,cbp

      INTEGER :: i,j,k,iseg,ir,iz,ip,nseg

      REAL(r8) :: phi,rx,ry,rz,dl,cosang,sinang,dlmin,dbx,dby,dbz
      REAL(r8), DIMENSION(0:nr,0:nz,0:np) :: xobs,yobs,zobs,
     $     bx,by,bz,br,bp
      COMPLEX(r8), DIMENSION(0:nr,0:nz) :: tbr,tbz,tbp
      REAL(r8), DIMENSION(:,:,:), POINTER :: xa,ya,za,
     $     dx,dy,dz

      TYPE(cspline_type) :: bspl

      dlmin=2e-3*cro
      
      DO i=1,coil_num
         nseg=coil(i)%nsec-1

         ALLOCATE(xa(coil(i)%ncoil,coil(i)%s,nseg),
     $        ya(coil(i)%ncoil,coil(i)%s,nseg),
     $        za(coil(i)%ncoil,coil(i)%s,nseg),
     $        dx(coil(i)%ncoil,coil(i)%s,nseg),
     $        dy(coil(i)%ncoil,coil(i)%s,nseg),
     $        dz(coil(i)%ncoil,coil(i)%s,nseg))
         
         DO j=1,coil(i)%ncoil
            DO k=1,coil(i)%s
               DO iseg=1,nseg
                  xa(j,k,iseg)=(coil(i)%x(j,k,iseg+1)+
     $                 coil(i)%x(j,k,iseg))/2.0
                  ya(j,k,iseg)=(coil(i)%y(j,k,iseg+1)+
     $                 coil(i)%y(j,k,iseg))/2.0      
                  za(j,k,iseg)=(coil(i)%z(j,k,iseg+1)+
     $                 coil(i)%z(j,k,iseg))/2.0
                  
                  dx(j,k,iseg)=coil(i)%x(j,k,iseg+1)-coil(i)%x(j,k,iseg)
                  dy(j,k,iseg)=coil(i)%y(j,k,iseg+1)-coil(i)%y(j,k,iseg)
                  dz(j,k,iseg)=coil(i)%z(j,k,iseg+1)-coil(i)%z(j,k,iseg)
               ENDDO
            ENDDO
         ENDDO
         
         CALL cspline_alloc(bspl,np,3)

         DO ir=0,nr
            DO iz=0,nz
               DO ip=0,np
                  phi=ip/REAL(np,8)
                  bspl%xs(ip)=phi
                  xobs(ir,iz,ip)=gdr(ir,iz)*COS(twopi*phi)
                  yobs(ir,iz,ip)=gdr(ir,iz)*SIN(twopi*phi)
                  zobs(ir,iz,ip)=gdz(ir,iz)
                  
                  dbx=0
                  dby=0
                  dbz=0

!$OMP PARALLEL PRIVATE(j,k,iseg,rx,ry,rz,dl) SHARED(i,ir,iz,ip)
!$OMP DO REDUCTION(+:dbx,dby,dbz) COLLAPSE(3)
                  DO j=1,coil(i)%ncoil
                     DO k=1,coil(i)%s
                        DO iseg=1,nseg
                           rx=xobs(ir,iz,ip)-xa(j,k,iseg)
                           ry=yobs(ir,iz,ip)-ya(j,k,iseg)
                           rz=zobs(ir,iz,ip)-za(j,k,iseg)
                           
                           dl=MAXVAL((/SQRT(rx**2+ry**2+rz**2),dlmin/))
                           
                           dbx=dbx+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                          dl**3*(dy(j,k,iseg)*rz-dz(j,k,iseg)*ry)
                           dby=dby+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                          dl**3*(dz(j,k,iseg)*rx-dx(j,k,iseg)*rz)
                           dbz=dbz+coil(i)%cur(j)*coil(i)%nw*1e-7/
     $                          dl**3*(dx(j,k,iseg)*ry-dy(j,k,iseg)*rx)            
                           
                        ENDDO
                     ENDDO
                  ENDDO
!$OMP END DO
!$OMP END PARALLEL
               
                  bx(ir,iz,ip)=dbx
                  by(ir,iz,ip)=dby
                  bz(ir,iz,ip)=dbz
                  cosang=xobs(ir,iz,ip)/gdr(ir,iz)
                  sinang=yobs(ir,iz,ip)/gdr(ir,iz)
                  br(ir,iz,ip)=bx(ir,iz,ip)*cosang+by(ir,iz,ip)*sinang
                  bp(ir,iz,ip)=by(ir,iz,ip)*cosang-bx(ir,iz,ip)*sinang

                  bspl%fs(ip,1)=br(ir,iz,ip)*EXP(twopi*ifac*cnn*phi)
                  bspl%fs(ip,2)=bz(ir,iz,ip)*EXP(twopi*ifac*cnn*phi)
                  bspl%fs(ip,3)=bp(ir,iz,ip)*EXP(twopi*ifac*cnn*phi)
                 
               ENDDO
               CALL cspline_fit(bspl,"periodic")
               CALL cspline_int(bspl)

               tbr(ir,iz)=2.0*bspl%fsi(np,1)
               tbz(ir,iz)=2.0*bspl%fsi(np,2)
               tbp(ir,iz)=2.0*bspl%fsi(np,3)

            ENDDO
         ENDDO

         CALL cspline_dealloc(bspl)

         DEALLOCATE(xa,ya,za,dx,dy,dz)

         cbr=cbr+tbr
         cbz=cbz+tbz
         cbp=cbp+tbp
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_bs_rzphi

      END MODULE field_mod
     







                     



