c-----------------------------------------------------------------------
c     GENERAL PERTURBED EQUILIBRIUM CONTROL
c     computation of plasma response
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. gpresp_mod
c      1. gpresp_eigen
c      2. gpresp_pinduct
c      3. gpresp_sinduct
c      4. gpresp_permeab
c      5. gpresp_reluct
c      6. gpresp_indrel
c-----------------------------------------------------------------------
c     subprogram 0. gpresp_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE gpresp_mod
      USE gpeq_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. gpresp_eigen.
c     construct flux and current matrices from eigenmodes.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_eigen
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j
      REAL(r8) :: chptsq,chpdsq,error

      COMPLEX(r8), DIMENSION(2,mpert) :: chpwmn,chpdif,chpwif
      CHARACTER(32) :: emessage
      IF(debug_flag) PRINT *, "Entering gpresp_eigen"
c-----------------------------------------------------------------------
c     build ideal solutions.
c-----------------------------------------------------------------------
      ALLOCATE(surfet(4,mpert),surfep(4,mpert),
     $     surfee(mpert),surfei(mpert),chperr(2,mpert),chpsqr(4,mpert))
      ALLOCATE(chimats(mpert,mpert),chemats(mpert,mpert),
     $     kaxmats(mpert,mpert),flxmats(mpert,mpert),
     $     chpmats(4,mpert,mpert),kapmats(4,mpert,mpert))
      IF(verbose) WRITE(*,*)"Building free boundary solutions"
      DO i=1,mpert
         edge_mn=0
         edge_flag=.FALSE.
         CALL idcon_build(i,edge_mn)
c-----------------------------------------------------------------------
c     compute the perturbed quantities and contruct hermitian matrices.
c-----------------------------------------------------------------------
         ALLOCATE(chi_mn(mpert),che_mn(mpert),chp_mn(4,mpert),
     $        kap_mn(4,mpert),kax_mn(mpert))
         CALL gpeq_alloc
         surface_flag=.FALSE.
         CALL gpeq_sol(psilim)
         surface_flag=.FALSE.
         CALL gpeq_contra(psilim)
         CALL gpeq_cova(psilim)
         CALL gpeq_surface(psilim)
c-----------------------------------------------------------------------
c     compute each fourier component.
c-----------------------------------------------------------------------
         flxmats(:,i)=bwp_mn
         chimats(:,i)=chi_mn
         chemats(:,i)=che_mn
         kaxmats(:,i)=kax_mn
         chpmats(:,:,i)=chp_mn
         kapmats(:,:,i)=kap_mn
c-----------------------------------------------------------------------
c     estimate errors in magnetic scalar potentials.
c-----------------------------------------------------------------------
         chpwmn(1,:)=chp_mn(1,:)
         chpdif(1,:)=chp_mn(1,:)-chp_mn(2,:)
         chpwmn(2,:)=chp_mn(3,:)
         chpdif(2,:)=chp_mn(3,:)-chp_mn(4,:)
         chpwif=chpdif
         DO j=1,2
            CALL gpeq_weight(psilim,chpwmn(j,:),mfac,mpert,1)
            CALL gpeq_weight(psilim,chpwif(j,:),mfac,mpert,1)
            chptsq=SUM(REAL(CONJG(chpwmn(j,:))*chp_mn(2*j-1,:),r8))
            chpdsq=SUM(REAL(CONJG(chpwif(j,:))*chpdif(j,:),r8))
            chperr(j,i)=SQRT(chpdsq/chptsq)
         ENDDO
         DO j=1,4
            chpsqr(j,i)=SUM(REAL(CONJG(chp_mn(j,:))*chp_mn(j,:),r8))
         ENDDO
c-----------------------------------------------------------------------
c     compute surface energy.
c-----------------------------------------------------------------------
         surfei(i)=0.5/mu0*SUM(REAL(chimats(:,i)*CONJG(flxmats(:,i))))
         surfee(i)=-0.5/mu0*SUM(REAL(chemats(:,i)*CONJG(flxmats(:,i))))

         IF(kin_flag)THEN
            DO j=1,4
               surfep(j,i)=0.5/mu0*SUM(chpmats(j,:,i)*
     $              CONJG(flxmats(:,i)))
               surfet(j,i)=surfep(j,i)+surfee(i)
            ENDDO
            error = ABS(1-surfet(2,i)/surfet(1,i))
            IF(error>1e-6)THEN
               emessage="   **WARNING: Large Error**"
            ELSE
               emessage=""
            ENDIF
            ! surface current needs dP_perp contributions
c            IF(verbose)WRITE(*,'(1x,a12,i3,2(a7,es10.3),a10,es10.3)')
c     $           "eigenmode = ",i,", dW = ",REAL(surfet(1,i)),
c     $           ", T = ",-2*nn*AIMAG(surfet(1,i)),", error = ",
c     $           ABS(1-surfet(2,i)/surfet(1,i))
            IF(verbose)
     $           WRITE(*,'(1x,a12,i3,2(a7,es10.3))')
     $           "eigenmode = ",i,", dW = ",REAL(et(i)),
     $           ", T = ",-2*nn*AIMAG(et(i))
         ELSE
            DO j=1,4
               surfep(j,i)=0.5/mu0*SUM(REAL(chpmats(j,:,i)*
     $              CONJG(flxmats(:,i))))
               surfet(j,i)=surfep(j,i)+surfee(i)
            ENDDO
            error = ABS(1-surfet(2,i)/surfet(1,i))
            IF(error>1e-6)THEN
               emessage="   **WARNING: Large Error**"
            ELSE
               emessage=""
            ENDIF
            IF(verbose.AND.((error>1e-6)))
     $           WRITE(*,'(1x,a12,i3,a7,es10.3,a10,es10.3)')
     $           "eigenmode = ",i,", dW = ",REAL(surfet(1,i)),
     $           ", error = ",ABS(1-surfet(2,i)/surfet(1,i))
         ENDIF

         CALL gpeq_dealloc
         DEALLOCATE(chi_mn,che_mn,chp_mn,kap_mn,kax_mn)
      ENDDO
      DEALLOCATE(grri,grre)
      IF(debug_flag) PRINT *, "->Leaving gpresp_eigen"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpresp_eigen
c-----------------------------------------------------------------------
c     subprogram 2. gpresp_pinduct.
c     construct plasma inductance matrices.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_pinduct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,k,lwork
      COMPLEX(r8) :: t1,t2
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(2*mpert) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,vr,vl
c-----------------------------------------------------------------------
c     calculate plasma inductance matrix by surface consideration.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)"Calculating inductance and permeability"
      ALLOCATE(plas_indmats(0:4,mpert,mpert),
     $   plas_indinvmats(0:4,mpert,mpert),
     $   plas_indev(0:4,mpert),plas_indevmats(0:4,mpert,mpert),
     $   plas_indinvev(0:4,mpert),plas_indinvevmats(0:4,mpert,mpert))
      plas_indmats = 0
      plas_indev = 0
      plas_indevmats = 0
      plas_indinvmats = 0
      plas_indinvev = 0
      plas_indinvevmats = 0

      lwork=2*mpert+1
      DO j=1,4
         ! plasma inductance
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         temp1=TRANSPOSE(kapmats(j,:,:))
         temp2=TRANSPOSE(flxmats)
         CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         temp1=TRANSPOSE(temp2)
         plas_indmats(j,:,:)=temp1
         work=0
         rwork=0
         CALL zgeev('V','V',mpert,temp1,mpert,plas_indev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         plas_indevmats(j,:,:)=vr
      ENDDO
c-----------------------------------------------------------------------
c     calculate energy inductance matrix by energy consideration.
c-----------------------------------------------------------------------
      IF(resp_induct_flag .AND. .NOT. galsol%gal_flag)THEN
         temp1=0
         DO i=1,mpert
            DO j=1,mpert
               t1=-1/(chi1*(mfac(i)-nn*qlim)*twopi*ifac)
               t2=1/(chi1*(mfac(j)-nn*qlim)*twopi*ifac)
               temp2(i,j)=2.0*t1*wt0(i,j)*t2
               IF(i==j)temp1(i,j)=1.0
            ENDDO
         ENDDO
         CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         plas_indmats(0,:,:)=temp1
      ELSEIF(.NOT. galsol%gal_flag)THEN
         temp1=0
         temp2=flxmats
         DO i=1,mpert
            DO j=1,mpert
               DO k=1,mpert
                  temp1(i,j)=temp1(i,j)+temp2(i,k)*CONJG(temp2(j,k))
     $                 /(et(k)*2.0)
               ENDDO
            ENDDO
         ENDDO
         plas_indmats(0,:,:)=temp1
      ENDIF
      CALL zgeev('V','V',mpert,temp1,mpert,plas_indev(0,:),
     $     vl,mpert,vr,mpert,work,lwork,rwork,info)
      plas_indevmats(0,:,:)=vr
c-----------------------------------------------------------------------
c     calculate inverse of the plasma inductance ~ Energy
c-----------------------------------------------------------------------
      plas_indinvevmats = 0
      DO j=0,4
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         temp1=0
         DO i=1,mpert
            temp1(i,i)=1
         ENDDO
         temp2=plas_indmats(j,:,:)
         CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         plas_indinvmats(j,:,:)=temp1
         work=0
         rwork=0
         CALL zgeev('V','V',mpert,temp1,mpert,plas_indinvev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         plas_indinvevmats(j,:,:)=vr
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpresp_pinduct
c-----------------------------------------------------------------------
c     subprogram 3. gpresp_sinduct.
c     construct surface inductance matrix.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_sinduct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2
c-----------------------------------------------------------------------
c     calculate surface inductance matrix by vacuum consideration.
c-----------------------------------------------------------------------
      ALLOCATE(surf_indev(mpert),surf_indmats(mpert,mpert),
     $     surf_indevmats(mpert,mpert),surf_indinvev(mpert),
     $     surf_indinvmats(mpert,mpert),surf_indinvevmats(mpert,mpert))
      temp1=TRANSPOSE(kaxmats)
      temp2=TRANSPOSE(flxmats)
      CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
      CALL zgetrs('N',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
      temp1=TRANSPOSE(temp2)
      temp1=0.5*(temp1+CONJG(TRANSPOSE(temp1)))
      surf_indmats=temp1
      lwork=2*mpert-1
      CALL zheev('V','U',mpert,temp1,mpert,surf_indev,work,
     $     lwork,rwork,info)
      surf_indevmats=temp1
c-----------------------------------------------------------------------
c     calculate inverse of surface inductance
c-----------------------------------------------------------------------
      work = 0
      rwork = 0
      work2=0
      temp1=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2 = surf_indmats
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
      surf_indinvmats = temp1
      CALL zheev('V','U',mpert,temp1,mpert,surf_indinvev,work,
     $     lwork,rwork,info)
      surf_indinvevmats=temp1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpresp_sinduct
c-----------------------------------------------------------------------
c     subprogram 4. gpresp_permeab.
c     construct permeability matrix.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_permeab
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,k,ii,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(2*mpert) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2,vr,vl
      COMPLEX(r8) :: ev

      REAL(r8), DIMENSION(5*mpert) :: rworks
      COMPLEX(r8), DIMENSION(3*mpert) :: works
c-----------------------------------------------------------------------
c     calculate permeability matrix.
c-----------------------------------------------------------------------
      ALLOCATE(permeabev(0:4,mpert),permeabmats(0:4,mpert,mpert),
     $     permeabevmats(0:4,mpert,mpert),permeabindex(0:4,mpert),
     $     permeabinvmats(0:4,mpert,mpert))
      lwork=2*mpert+1
      DO j=0,4
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         work=0
         rwork=0
         temp1=TRANSPOSE(surf_indmats)
         temp2=TRANSPOSE(plas_indmats(j,:,:))
         CALL zhetrf('L',mpert,temp1,mpert,ipiv,work2,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         temp1=TRANSPOSE(temp2)
         permeabmats(j,:,:)=temp1
         CALL zgeev('V','V',mpert,temp1,mpert,permeabev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         permeabevmats(j,:,:)=vr

         ! calculate inverse.
         temp1=0
         DO i=1,mpert
            temp1(i,i)=1
         ENDDO
         temp2=permeabmats(j,:,:)
         CALL zgetrf(mpert,mpert,temp2,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         permeabinvmats(j,:,:) = temp1
c-----------------------------------------------------------------------
c     sort permeability eigenvalues by absolute value
c      - sorting taken from http://stackoverflow.com/questions/8834585/sorting-eigensystem-obtained-from-zgeev
c      - which is from the end of zsteqr.f
c-----------------------------------------------------------------------
         temp1 = permeabmats(j,:,:)
         DO i = 1, mpert-1
            k = i
            ev = permeabev(j,i)
            DO ii = i+1, mpert
               IF( ABS(permeabev(j,ii)) > ABS(ev) ) THEN
                  K = ii
                  ev = permeabev(j,ii)
               ENDIF
            ENDDO
            IF( K.NE.I ) THEN
               permeabev(j,k) = permeabev(j,i)
               permeabev(j,i) = ev
               CALL zswap( mpert, permeabevmats(j,:,i), 1,
     $                            permeabevmats(j,:,k), 1)
            END IF
         ENDDO
      ENDDO
      IF(verbose) WRITE(*,'(1x,a,es10.3)')"Single mode permeability = ",
     $     MAXVAL(ABS(permeabev(resp_index,:)))
c-----------------------------------------------------------------------
c     SVD permeability matrix for orthonormal basis.
c-----------------------------------------------------------------------
      ALLOCATE(permeabsv(0:4,mpert),permeabsvmats(0:4,mpert,mpert))
      DO j=0,4
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         lwork=3*mpert
         works=0
         rworks=0
         temp1=permeabmats(resp_index,:,:)
         CALL zgesvd('S','S',mpert,mpert,temp1,mpert,permeabsv(j,:),
     $        vl,mpert,vr,mpert,works,lwork,rworks,info)
         permeabsvmats(j,:,:)=CONJG(TRANSPOSE(vr))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpresp_permeab
c-----------------------------------------------------------------------
c     subprogram 5. gpresp_reluct.
c     construct reluctance matrix.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_reluct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(2*mpert) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2,vl,vr
c-----------------------------------------------------------------------
c     calculate reluctance matrix.
c-----------------------------------------------------------------------
      ALLOCATE(reluctev(0:4,mpert),diff_indmats(0:4,mpert,mpert),
     $     reluctmats(0:4,mpert,mpert),reluctevmats(0:4,mpert,mpert))
      DO j=0,4
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         work=0
         work2=0
         rwork=0
         temp1=0
         DO i=1,mpert
            temp1(i,i)=1
         ENDDO
         temp2=surf_indmats
         CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         diff_indmats(j,:,:)=CONJG(TRANSPOSE(plas_indmats(j,:,:)))
     $        -surf_indmats
         reluctmats(j,:,:)=MATMUL(temp1,
     $        MATMUL(diff_indmats(j,:,:),temp1))
         temp1=reluctmats(j,:,:)

         lwork=2*mpert+1
         CALL zgeev('V','V',mpert,temp1,mpert,reluctev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         reluctevmats(j,:,:)=vr
      ENDDO
      RETURN
      END SUBROUTINE gpresp_reluct
c-----------------------------------------------------------------------
c     subprogram 6. gpresp_indrel.
c     construct combined matrix with surface inductance and reluctance.
c-----------------------------------------------------------------------
      SUBROUTINE gpresp_indrel
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(2*mpert) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2,vl,vr
c-----------------------------------------------------------------------
c     calculate reluctance matrix.
c-----------------------------------------------------------------------
      ALLOCATE(indrelev(0:4,mpert),indrelmats(0:4,mpert,mpert),
     $     indrelevmats(0:4,mpert,mpert))
      DO j=0,4
         IF(galsol%gal_flag .AND. j/=1) CYCLE
         work=0
         work2=0
         rwork=0
         temp1=0
         DO i=1,mpert
            temp1(i,i)=1
         ENDDO
         temp2=surf_indmats
         CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         indrelmats(j,:,:)=MATMUL(temp1,
     $        MATMUL(plas_indmats(j,:,:),temp1))
         temp1=indrelmats(j,:,:)
         lwork=2*mpert+1
         CALL zgeev('V','V',mpert,temp1,mpert,indrelev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         indrelevmats(j,:,:)=vr
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gpresp_indrel

      END MODULE gpresp_mod
