c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     IPRESP: computation of plasma response
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. ipresp_mod
c      1. ipresp_eigen
c      2. ipresp_pinduct
c      3. ipresp_sinduct
c      4. ipresp_permeab
c      5. ipresp_reluct
c      6. ipresp_indrel
c-----------------------------------------------------------------------
c     subprogram 0. ipresp_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ipresp_mod
      USE ipeq_mod

      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ipresp_eigen.
c     construct flux and current matrices from eigenmodes.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_eigen
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j
      REAL(r8) :: chptsq,chpdsq,error

      COMPLEX(r8), DIMENSION(2,mpert) :: chpwmn,chpdif,chpwif
      CHARACTER(32) :: emessage
      IF(debug_flag) PRINT *, "Entering ipresp_eigen"      
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
         CALL ipeq_alloc
         CALL ipeq_sol(psilim)
         CALL ipeq_contra(psilim)
         CALL ipeq_cova(psilim)
         CALL ipeq_surface(psilim)
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
            CALL ipeq_weight(psilim,chpwmn(j,:),mfac,mpert,1)
            CALL ipeq_weight(psilim,chpwif(j,:),mfac,mpert,1)
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
         surfei(i)=0.5/mu0*SUM(REAL(CONJG(chimats(:,i))*
     $        flxmats(:,i),r8))
         surfee(i)=-0.5/mu0*SUM(REAL(CONJG(chemats(:,i))*
     $        flxmats(:,i),r8))
         DO j=1,4
            surfep(j,i)=0.5/mu0*SUM(REAL(CONJG(chpmats(j,:,i))*
     $           flxmats(:,i),r8))
            surfet(j,i)=surfep(j,i)+surfee(i)
         ENDDO
         error = ABS(1-surfet(2,i)/surfet(1,i))
         IF(error>1e-6)THEN
            emessage="   **WARNING: Large Error**"
         ELSE
            emessage=""
         ENDIF
         IF(verbose.AND.((i==1).OR.(i==mpert).OR.(error>1e-6)))
     $      WRITE(*,'(1x,a12,i3,a7,es10.3,a10,es10.3,a32)')
     $      "eigenmode = ",i,", dw = ",surfet(1,i),", error = ",error,
     $      emessage
         CALL ipeq_dealloc
         DEALLOCATE(chi_mn,che_mn,chp_mn,kap_mn,kax_mn)
      ENDDO
      DEALLOCATE(grri,grre)
      IF(debug_flag) PRINT *, "->Leaving ipresp_eigen"      
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipresp_eigen
c-----------------------------------------------------------------------
c     subprogram 2. ipresp_pinduct.
c     construct plasma inductance matrices.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_pinduct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,k,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2
c-----------------------------------------------------------------------
c     calculate plasma inductance matrix by surface consideration.
c-----------------------------------------------------------------------
      IF(verbose) WRITE(*,*)"Calculating inductrances and permeability"
      ALLOCATE(plas_indmats(0:4,mpert,mpert),
     $   pinv_indmats(0:4,mpert,mpert),
     $   plas_indev(0:4,mpert),plas_indevmats(0:4,mpert,mpert),
     $   pinv_indev(0:4,mpert),pinv_indevmats(0:4,mpert,mpert))
      DO j=1,4
         work=0
         work2=0
         rwork=0
         temp1=TRANSPOSE(kapmats(j,:,:))
         temp2=TRANSPOSE(flxmats)
         CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
         CALL zgetrs('N',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         temp1=TRANSPOSE(temp2)
         temp1=0.5*(temp1+CONJG(TRANSPOSE(temp1)))
         plas_indmats(j,:,:)=temp1
         lwork=2*mpert-1
         CALL zheev('V','U',mpert,temp1,mpert,plas_indev(j,:),work,
     $        lwork,rwork,info)
         !temp2=0
         !DO i=1,mpert
         !   temp2(i,i)=temp1(i,i)
         !ENDDO
         !temp1=temp1+CONJG(TRANSPOSE(temp1))-temp2 !!whyyyyy??
         plas_indevmats(j,:,:)=temp1
         ! Inverse of the plasma inductance ~ Energy
         lwork=2*mpert-1
         work=0
         work2=0
         rwork=0
         temp1=0
         DO i=1,mpert
            temp1(i,i)=1
         ENDDO
         temp2=plas_indmats(j,:,:)
         CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
         pinv_indmats(j,:,:) = temp1
         CALL zheev('V','U',mpert,temp1,mpert,pinv_indev(j,:),work,
     $        lwork,rwork,info)
         pinv_indevmats(j,:,:)=temp1
      ENDDO
c-----------------------------------------------------------------------
c     calculate energy inductance matrix by energy consideration.
c-----------------------------------------------------------------------
      work=0
      rwork=0
      temp1=0
      temp2=flxmats
      DO i=1,mpert
         DO j=1,mpert
            DO k=1,mpert
               temp1(i,j)=temp1(i,j)+temp2(i,k)*CONJG(temp2(j,k))
     $              /(et(k)*2.0)
            ENDDO
         ENDDO
      ENDDO      
      lwork=2*mpert-1
      plas_indmats(0,:,:)=temp1
      CALL zheev('V','U',mpert,temp1,mpert,plas_indev(0,:),work,
     $     lwork,rwork,info)
      plas_indevmats(0,:,:)=temp1
      ! Inverse of the plasma inductance ~ Energy
      lwork=2*mpert-1
      work=0
      work2=0
      rwork=0
      temp1=0
      DO i=1,mpert
         temp1(i,i)=1
      ENDDO
      temp2=2*plas_indmats(0,:,:)
      CALL zhetrf('L',mpert,temp2,mpert,ipiv,work2,mpert*mpert,info)
      CALL zhetrs('L',mpert,mpert,temp2,mpert,ipiv,temp1,mpert,info)
      pinv_indmats(0,:,:) = temp1
      CALL zheev('V','U',mpert,temp1,mpert,pinv_indev(0,:),work,
     $     lwork,rwork,info)
      pinv_indevmats(0,:,:)=temp1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipresp_pinduct
c-----------------------------------------------------------------------
c     subprogram 3. ipresp_sinduct.
c     construct surface inductance matrix.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_sinduct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2
c-----------------------------------------------------------------------
c     calculate surface inductance matrix by vacuum consideration.
c-----------------------------------------------------------------------
      ALLOCATE(surf_indev(mpert),surf_indmats(mpert,mpert),
     $     surf_indevmats(mpert,mpert))
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
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipresp_sinduct
c-----------------------------------------------------------------------
c     subprogram 4. ipresp_permeab.
c     construct permeability matrix.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_permeab
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,k,ii,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(2*mpert) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert+1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2,vr,vl
      COMPLEX(r8) :: ev
      COMPLEX(r8), DIMENSION(mpert) :: vec
      
      REAL(r8), DIMENSION(5*mpert) :: rworks
      COMPLEX(r8), DIMENSION(3*mpert) :: works
      COMPLEX(r8), DIMENSION(mpert,mpert) :: ctp,a,permu
c-----------------------------------------------------------------------
c     calculate permeability matrix.
c-----------------------------------------------------------------------
      ALLOCATE(permeabev(0:4,mpert),permeabmats(0:4,mpert,mpert),
     $     permeabevmats(0:4,mpert,mpert),permeabindex(0:4,mpert))
      DO j=0,4
         work=0
         work2=0
         rwork=0
         temp1=TRANSPOSE(surf_indmats)
         temp2=TRANSPOSE(plas_indmats(j,:,:))
         CALL zhetrf('L',mpert,temp1,mpert,ipiv,work2,mpert*mpert,info)
         CALL zhetrs('L',mpert,mpert,temp1,mpert,ipiv,temp2,mpert,info)
         temp1=TRANSPOSE(temp2)
         permeabmats(j,:,:)=temp1
c-----------------------------------------------------------------------
c     calculate permeability eigenvalues and vectors, then sort them.
c      - sorting taken from http://stackoverflow.com/questions/8834585/sorting-eigensystem-obtained-from-zgeev
c      - which is from the end of zsteqr.f
c-----------------------------------------------------------------------
         lwork=2*mpert+1
         CALL zgeev('V','V',mpert,temp1,mpert,permeabev(j,:),
     $        vl,mpert,vr,mpert,work,lwork,rwork,info)
         permeabevmats(j,:,:)=vr
         ! sort by absolute value of eigenvector
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
c-----------------------------------------------------------------------
c     LOGAN - SVD permeability matrix for orthonormal basis.
c-----------------------------------------------------------------------
      ALLOCATE(permeabsv(0:4,mpert),permeabsvmats(0:4,mpert,mpert))
      DO j=0,4
         lwork=3*mpert
         works=0
         rworks=0
         temp1=permeabmats(resp_index,:,:)
         CALL zgesvd('S','S',mpert,mpert,temp1,mpert,permeabsv(j,:),
     $        vl,mpert,vr,mpert,works,lwork,rworks,info)
         permeabsvmats(j,:,:)=CONJG(TRANSPOSE(vr))
      ENDDO
c-----------------------------------------------------------------------
c     sort by amplitudes.
c-----------------------------------------------------------------------
      DO j=0,4
         DO i=1,mpert
            permeabindex(j,i)=i
         ENDDO
         CALL isbubble(REAL(permeabev(j,:)),
     $        permeabindex(j,:),1,mpert)
      ENDDO
      IF(verbose) WRITE(*,'(1x,a,es10.3)')"Single mode permeability = ",
     $     MAXVAL(ABS(permeabev(resp_index,:)))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipresp_permeab
c-----------------------------------------------------------------------
c     subprogram 5. ipresp_reluct.
c     construct reluctance matrix.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_reluct
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert) :: temp
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2,sqrta
c-----------------------------------------------------------------------
c     calculate reluctance matrix.
c-----------------------------------------------------------------------
      ALLOCATE(reluctev(0:4,mpert),diff_indmats(0:4,mpert,mpert),
     $     reluctmats(0:4,mpert,mpert),reluctevmats(0:4,mpert,mpert))
      DO j=0,4
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
         diff_indmats(j,:,:)=plas_indmats(j,:,:)-surf_indmats
         reluctmats(j,:,:)=MATMUL(temp1,
     $        MATMUL(diff_indmats(j,:,:),temp1))
         temp1=reluctmats(j,:,:)
         lwork=2*mpert-1
         CALL zheev('V','U',mpert,temp1,mpert,reluctev(j,:),work,
     $        lwork,rwork,info)
         reluctevmats(j,:,:)=temp1
      ENDDO
      RETURN
      END SUBROUTINE ipresp_reluct
c-----------------------------------------------------------------------
c     subprogram 6. ipresp_indrel.
c     construct combined matrix with surface inductance and reluctance.
c-----------------------------------------------------------------------
      SUBROUTINE ipresp_indrel
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: i,j,lwork
      INTEGER, DIMENSION(mpert):: ipiv
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp1,temp2,work2
c-----------------------------------------------------------------------
c     calculate reluctance matrix.
c-----------------------------------------------------------------------
      ALLOCATE(indrelev(0:4,mpert),indrelmats(0:4,mpert,mpert),
     $     indrelevmats(0:4,mpert,mpert))
      DO j=0,4
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
         lwork=2*mpert-1
         CALL zheev('V','U',mpert,temp1,mpert,indrelev(j,:),work,
     $        lwork,rwork,info)
         indrelevmats(j,:,:)=temp1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ipresp_indrel
      
      END MODULE ipresp_mod
