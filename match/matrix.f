c-----------------------------------------------------------------------
c     file matrix.f.
c     computes matching matrix and its determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. matrix_mod.
c     1. matrix_run.
c     2. matrix_fixup.
c     3. matrix_make.
c     4. matrix_layer.
c     5. matrix_diag.
c     6. matrix_poly_form1.
c     7. matrix_poly_form2.
c     8. matrix_poly_eval.
c     9. matrix_poly_test.
c     10. matrix_write1.
c     11. matrix_write2.
c     12. matrix_test_fixup.
c     13. matrix_scan.
c     14. matrix_condense.
c-----------------------------------------------------------------------
c     subprogram 0. matrix_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE matrix_mod
      USE match_mod
      USE debug_mod
      USE qsort_mod
      IMPLICIT NONE

      LOGICAL :: full=.TRUE.,transform_flag=.TRUE.,condense=.FALSE.
      INTEGER :: poly_form_type=1,poly_test_type=1,nscan,m1
      REAL(r8) :: ff0=1,ff1=1,scan0,scan1,delta0,delta1
      COMPLEX(r8) :: det0
      COMPLEX(r8), DIMENSION(:), POINTER :: diag0,diag
      COMPLEX(r8), DIMENSION(:,:), POINTER :: m11

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. matrix_run.
c     build fixup matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_run

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ising,ifix
c-----------------------------------------------------------------------
c     fixup solutions.
c-----------------------------------------------------------------------
      IF(transform_flag)THEN
         IF(diagnose)CALL matrix_test_fixup
         WRITE(*,*)"Transform asymptotic coefficients."
         DO ising=1,msing
            DO ifix=singtype(ising)%jfix+1,mfix
               CALL matrix_fixup(singtype(ising)%ca_l,ifix)
               CALL matrix_fixup(singtype(ising)%ca_r,ifix)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      IF(condense)THEN
         m1=2*msing
      ELSE
         m1=mmatch
      ENDIF
      ALLOCATE(diag0(m1),diag(m1),m11(m1,m1))
c-----------------------------------------------------------------------
c     form and test matching determinant.
c-----------------------------------------------------------------------
      WRITE(*,*)"Form polynomial."
      WRITE(*,'(2(a,i1),a,l1)')" poly_form_type = ",poly_form_type,
     $     ", poly_test_type = ",poly_test_type,", condense = ",condense
      IF(poly_form_type == 1)THEN
         CALL matrix_poly_form1
         CALL matrix_poly_test
      ELSEIF(poly_form_type == 2)THEN
         CALL matrix_poly_form2
         CALL matrix_poly_test
      ELSE
         CALL matrix_poly_scan
      ENDIF
      DEALLOCATE(diag0,diag)
      IF(condense)DEALLOCATE(m11)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_run
c-----------------------------------------------------------------------
c     subprogram 2. matrix_fixup.
c     performs Gaussian reduction of solution matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_fixup(ca,ifix)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: ca
      INTEGER, INTENT(IN) :: ifix

      INTEGER :: isol,jsol,ksol
      INTEGER, DIMENSION(:), POINTER :: index
      COMPLEX(r8), DIMENSION(:,:), POINTER :: fixfac
c-----------------------------------------------------------------------
c     triangularize primary solutions.
c-----------------------------------------------------------------------
      index => fixtype(ifix)%index
      fixfac => fixtype(ifix)%fixfac
      DO isol=1,mpert
         ksol=fixtype(ifix)%index(isol)
         DO jsol=1,SIZE(ca,2)
            IF(jsol == ksol)CYCLE
            ca(:,jsol,:)=ca(:,jsol,:)+ca(:,ksol,:)*fixfac(ksol,jsol)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_fixup
c-----------------------------------------------------------------------
c     subprogram 3. matrix_make.
c     make matching matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_make(deltar)

      COMPLEX(r8), DIMENSION(msing,2), INTENT(IN) :: deltar

      INTEGER :: ising,i,j,msol_l,msol_r
c-----------------------------------------------------------------------
c     fill matrix.
c-----------------------------------------------------------------------
      match=0
      i=0
      j=0
      DO ising=1,msing
         msol_l=singtype(ising)%msol_l
         msol_r=singtype(ising)%msol_r
         CALL matrix_layer(singtype(ising),deltar(ising,:),
     $        match(i+1:i+msol_r,j+1:j+msol_l+msol_r))
         i=i+msol_r
         j=j+msol_l
      ENDDO
      match(i+1:mmatch,j+1:mmatch)=soltype(mstep)%u(:,:,1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_make
c-----------------------------------------------------------------------
c     subprogram 4. matrix_layer.
c     sets up matching conditions across one singular layer.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_layer(st,deltar,match)

      TYPE(sing_type), INTENT(IN) :: st
      COMPLEX(r8), DIMENSION(2), INTENT(IN) :: deltar
      COMPLEX(r8), DIMENSION(:,:), INTENT(OUT) :: match

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(64) :: format1,format2
      INTEGER :: ipert,jpert,kpert,msol_l,msol_r,isol,index,m,n,i
      COMPLEX(r8), DIMENSION(2*(mpert-1),st%msol_l) :: cl
      COMPLEX(r8), DIMENSION(2*(mpert-1),st%msol_r) :: cr
      COMPLEX(r8), DIMENSION(st%msol_r,st%msol_l+st%msol_r) :: temp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"i",',i2,'(6x,i2,3x)/)')
 20   FORMAT('(i4,1p,',i2,'e11.3)')
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      msol_l=st%msol_l
      msol_r=st%msol_r
      jpert=st%jpert
      index=fixtype(st%jfix)%index(1)
      match=0
c-----------------------------------------------------------------------
c     resonant components.
c-----------------------------------------------------------------------
      match(1,1:msol_l)
     $     =st%ca_l(jpert,:,1)*deltar(2)-st%ca_l(jpert,:,2)
      match(1,msol_l+1:msol_l+msol_r)
     $     =-st%ca_r(jpert,:,1)*deltar(2)+st%ca_r(jpert,:,2)
      match(2,1:msol_l)
     $     =st%ca_l(jpert,:,1)*deltar(1)-st%ca_l(jpert,:,2)
      match(2,msol_l+1:msol_l+msol_r)
     $     =st%ca_r(jpert,:,1)*deltar(1)-st%ca_r(jpert,:,2)
      isol=3
c-----------------------------------------------------------------------
c     nonresonant components, msol_r < 2*mpert.
c-----------------------------------------------------------------------
      IF(msol_r < 2*mpert)THEN
         kpert=0
         DO ipert=1,mpert
            IF(ipert == jpert)CYCLE
            cl(kpert+1,:)=st%ca_l(ipert,:,1)
            cr(kpert+1,:)=st%ca_r(ipert,:,1)
            cl(kpert+2,:)=st%ca_l(ipert,:,2)
            cr(kpert+2,:)=st%ca_r(ipert,:,2)
            kpert=kpert+2
         ENDDO
         temp(:,1:msol_l)
     $        =MATMUL(CONJG(TRANSPOSE(cr)),cl)
         temp(:,msol_l+1:msol_l+msol_r)
     $        =-MATMUL(CONJG(TRANSPOSE(cr)),cr)
         match(3:msol_r,:)=temp(2:msol_r-1,:)
c-----------------------------------------------------------------------
c     nonresonant components, msol_r = 2*mpert.
c-----------------------------------------------------------------------
      ELSE
         DO ipert=1,mpert
            IF(ipert == jpert)CYCLE
            match(isol,1:msol_l)=st%ca_l(ipert,:,1)
            match(isol,msol_l+1:msol_l+msol_r)=-st%ca_r(ipert,:,1)
            isol=isol+1
            match(isol,1:msol_l)=st%ca_l(ipert,:,2)
            match(isol,msol_l+1:msol_l+msol_r)=-st%ca_r(ipert,:,2)
            isol=isol+1
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     diagnose asymptotic coefficients.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL ascii_open(debug_unit,"layer.out","UNKNOWN")
         CALL debug1(st%ca_l,"ca_l",debug_unit)
         CALL debug1(st%ca_r,"ca_r",debug_unit)
         m=SIZE(match,1)
         n=SIZE(match,2)
         WRITE(format1,10)n
         WRITE(format2,20)n
         WRITE(debug_unit,format1)(i,i=1,n)
         DO i=1,m
            WRITE(debug_unit,format2)i,REAL(match(i,:))
         ENDDO
         WRITE(debug_unit,format1)(i,i=1,n)
         CALL ascii_close(debug_unit)
         CALL program_stop("Termination by matrix_layer")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_layer
c-----------------------------------------------------------------------
c     subprogram 5. matrix_diag.
c     evaluates diagonal elements of matching matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_diag(diag)

      COMPLEX(r8), DIMENSION(:), INTENT(OUT) :: diag

      INTEGER :: i,info
      INTEGER, DIMENSION(m1) :: ipiv
c-----------------------------------------------------------------------
c     compute diagonal elements by LU factorization.
c-----------------------------------------------------------------------
      IF(condense)THEN
         CALL matrix_condense(m11)
      ELSE
         m11=match
      ENDIF
c-----------------------------------------------------------------------
c     compute diagonal elements by LU factorization.
c-----------------------------------------------------------------------
      CALL zgetrf(m1,m1,m11,m1,ipiv,info)
      diag=(/(m11(i,i),i=1,m1)/)
      DO i=1,m1
         IF(ipiv(i) /= i)diag(i)=-diag(i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_diag
c-----------------------------------------------------------------------
c     subprogram 6. matrix_poly_form1.
c     determines functional form of matching determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_poly_form1

      CHARACTER(80) :: format1,format2
      INTEGER :: iterm,ibit,mfac,ising
      COMPLEX(r8) :: deta,det
      COMPLEX(r8), DIMENSION(msing,2) :: deltar

      REAL(r8), DIMENSION(:,:), POINTER :: data
      INTEGER, DIMENSION(:), POINTER :: index
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"it"',i2,'x,"it",',i2,'x,"re ff",6x,"im ff"/)')
 20   FORMAT('(i6,2x,b',i2,'.',i2,',1p,2e11.3)')
c-----------------------------------------------------------------------
c     allocate and zero form factors.
c-----------------------------------------------------------------------
      mbit=2*msing
      mterm=ISHFT(1,mbit)-1
      ALLOCATE(ff(0:mterm),data(1,0:mterm),index(0:mterm))
      ff=0
c-----------------------------------------------------------------------
c     base case.
c-----------------------------------------------------------------------
      deltar=0
      CALL matrix_make(deltar)
      CALL matrix_diag(diag0)
      ff(0)=1
      det0=1
c-----------------------------------------------------------------------
c     compute form factors.
c-----------------------------------------------------------------------
      DO iterm=1,mterm
         deltar=0
         mfac=0
         DO ibit=0,mbit-1
            IF(BTEST(iterm,ibit))THEN
               deltar(ibit/2+1,MOD(ibit,2)+1)=ff0
               mfac=mfac+1
            ENDIF
         ENDDO
         CALL matrix_poly_eval(deltar,deta)
         CALL matrix_make(deltar)
         CALL matrix_diag(diag)
         det=PRODUCT(diag/diag0)-deta
         ff(iterm)=det/ff0**mfac
      ENDDO
c-----------------------------------------------------------------------
c     sort ff.
c-----------------------------------------------------------------------
      data(1,:)=ABS(ff)
      CALL qsort(data,index)
      index=index(mterm:0:-1)-1
c-----------------------------------------------------------------------
c     diagnose form factors.
c-----------------------------------------------------------------------
      WRITE(format1,10)mbit/2+1,mbit/2+3
      WRITE(format2,20)mbit,mbit
      WRITE(out_unit,'(4x,a,i1,1p,2(a,e9.3))')
     $     "poly_form_type = ",poly_form_type,
     $     ", ff0 = ",ff0,", ff1 = ",ff1
      WRITE(out_unit,format1)
      WRITE(out_unit,format2)(index(iterm),index(iterm),
     $     ff(index(iterm)),iterm=0,mterm)
      WRITE(out_unit,format1)
      DEALLOCATE(data,index)
c-----------------------------------------------------------------------
c     compute deltap.
c-----------------------------------------------------------------------
      DO ising=1,msing
         iterm=IBCLR(mterm,2*ising-1)
         deltap1(ising)=-2*ff(iterm)/ff(mterm)
         deltap2(ising)=-2*ff(iterm+1)/ff(mterm)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_poly_form1
c-----------------------------------------------------------------------
c     subprogram 7. matrix_poly_form2.
c     determines functional form of matching determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_poly_form2

      LOGICAL, DIMENSION(:,:), POINTER :: mask
      CHARACTER(80) :: format1,format2
      INTEGER :: iterm,jterm,ibit,ising,info
      INTEGER, DIMENSION(:), POINTER :: ipiv
      COMPLEX(r8), DIMENSION(msing,2) :: deltar
      COMPLEX(r8), DIMENSION(:,:), POINTER :: a

      REAL(r8), DIMENSION(:,:), POINTER :: data
      INTEGER, DIMENSION(:), POINTER :: index
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"it"',i2,'x,"it",',i2,'x,"re ff",6x,"im ff"/)')
 20   FORMAT('(i6,2x,b',i2,'.',i2,',1p,2e11.3)')
c-----------------------------------------------------------------------
c     base cases.
c-----------------------------------------------------------------------
      deltar=ff0
      CALL matrix_make(deltar)
      CALL matrix_diag(diag0)
      deltar=0
      CALL matrix_make(deltar)
      CALL matrix_diag(diag)
      det0=PRODUCT(diag/diag0)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      mbit=2*msing
      mterm=ISHFT(1,mbit)-1
      ALLOCATE(ff(0:mterm),ipiv(mterm),a(mterm,mterm),mask(mterm,mterm),
     $     data(1,0:mterm),index(0:mterm))
      ff(0)=1
c-----------------------------------------------------------------------
c     compute rhs and set mask.
c-----------------------------------------------------------------------
      mask=.TRUE.
      DO iterm=1,mterm
         DO ibit=0,mbit-1
            IF(BTEST(iterm,ibit))THEN
               deltar(ibit/2+1,MOD(ibit,2)+1)=ff0
            ELSE
               deltar(ibit/2+1,MOD(ibit,2)+1)=0
               DO jterm=1,mterm
                  IF(BTEST(jterm,ibit))mask(iterm,jterm)=.FALSE.
               ENDDO
            ENDIF
         ENDDO
         CALL matrix_make(deltar)
         CALL matrix_diag(diag)
         ff(iterm)=PRODUCT(diag/diag0)/det0-1
      ENDDO
c-----------------------------------------------------------------------
c     compute matrix.
c-----------------------------------------------------------------------
      WHERE(mask)
         a=1
      ELSEWHERE
         a=0
      ENDWHERE
c-----------------------------------------------------------------------
c     compute form factors.
c-----------------------------------------------------------------------
      CALL zgetrf(mterm,mterm,a,mterm,ipiv,info)
      CALL zgetrs("N",mterm,1,a,mterm,ipiv,ff(1:mterm),mterm,info)
      DEALLOCATE(a,ipiv)
c-----------------------------------------------------------------------
c     rescale form factors.
c-----------------------------------------------------------------------
      DO iterm=1,mterm
         DO ibit=0,mbit-1
            IF(BTEST(iterm,ibit))ff(iterm)=ff(iterm)/ff0
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     sort ff.
c-----------------------------------------------------------------------
      data(1,:)=ABS(ff)
      CALL qsort(data,index)
      index=index(mterm:0:-1)-1
c-----------------------------------------------------------------------
c     diagnose form factors.
c-----------------------------------------------------------------------
      WRITE(format1,10)mbit/2+1,mbit/2+3
      WRITE(format2,20)mbit,mbit
      WRITE(out_unit,'(4x,a,i1,1p,2(a,e9.3))')
     $     "poly_form_type = ",poly_form_type,
     $     ", ff0 = ",ff0,", ff1 = ",ff1
      WRITE(out_unit,format1)
      WRITE(out_unit,format2)(index(iterm),index(iterm),
     $     ff(index(iterm)),iterm=0,mterm)
      WRITE(out_unit,format1)
      DEALLOCATE(data,index)
c-----------------------------------------------------------------------
c     compute deltap.
c-----------------------------------------------------------------------
      DO ising=1,msing
         iterm=IBCLR(mterm,2*ising-1)
         deltap1(ising)=-2*ff(iterm)/ff(mterm)
         deltap2(ising)=-2*ff(iterm+1)/ff(mterm)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_poly_form2
c-----------------------------------------------------------------------
c     subprogram 8. matrix_poly_eval.
c     compute analytical form of determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_poly_eval(deltar,det)

      COMPLEX(r8), DIMENSION(msing,2), INTENT(IN) :: deltar
      COMPLEX(r8), INTENT(OUT) :: det

      INTEGER :: iterm,ibit
      COMPLEX(r8) :: term
c-----------------------------------------------------------------------
c     compute det.
c-----------------------------------------------------------------------
      det=ff(0)
      DO iterm=1,mterm
         term=ff(iterm)
         DO ibit=0,mbit-1
            IF(BTEST(iterm,ibit))THEN
               term=term*deltar(ibit/2+1,MOD(ibit,2)+1)
            ENDIF
         ENDDO
         det=det+term
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_poly_eval
c-----------------------------------------------------------------------
c     subprogram 9. matrix_poly_test
c     tests polynomial form of determinant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_poly_test

      INTEGER :: ibit,iterm
      REAL(r8) :: errmax,errmin
      REAL(r8), DIMENSION(0:mterm) :: error
      REAL(r8), DIMENSION(msing,2) :: parts
      COMPLEX(r8), DIMENSION(0:mterm) :: det,deta
      COMPLEX(r8), DIMENSION(msing,2) :: deltar

      REAL(r8), DIMENSION(1,0:mterm) :: data
      INTEGER, DIMENSION(0:mterm) :: index
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",5x,"j",8x,"det",12x,"deta",9x,"error"/)
 20   FORMAT(2i6,1p,2e16.8,e11.3)
c-----------------------------------------------------------------------
c     compute determinant both ways.
c-----------------------------------------------------------------------
      DO iterm=0,mterm
         IF(poly_test_type == 1)THEN
            CALL RANDOM_NUMBER(parts)
            deltar=(2*ff1-1)*parts
         ELSE
            DO ibit=0,mbit-1
               IF(BTEST(iterm,ibit))deltar(ibit/2+1,MOD(ibit,2)+1)=ff1
            ENDDO
         ENDIF
         CALL matrix_make(deltar)
         CALL matrix_diag(diag)
         det(iterm)=PRODUCT(diag/diag0)/det0
         CALL matrix_poly_eval(deltar,deta(iterm))
      ENDDO
c-----------------------------------------------------------------------
c     compute errors.
c-----------------------------------------------------------------------
      error=ABS(deta/det-1)
      errmax=MAXVAL(error)
      errmin=MINVAL(error)
c-----------------------------------------------------------------------
c     sort errors.
c-----------------------------------------------------------------------
      data(1,:)=error
      CALL qsort(data,index)
      index=index(mterm:0:-1)-1
c-----------------------------------------------------------------------
c     diagnose determinants and errors.
c-----------------------------------------------------------------------
      WRITE(out_unit,10)
      WRITE(out_unit,20)(iterm,index(iterm),REAL(det(index(iterm))),
     $     REAL(deta(index(iterm))),error(index(iterm)),iterm=0,mterm)
      WRITE(out_unit,10)
      WRITE(*,'(1x,1p,2(a,e9.3))')
     $     "errmax = ",errmax,", errmin = ",errmin
      WRITE(out_unit,'(1x,1p,2(a,e9.3)/)')
     $     "errmax = ",errmax,", errmin = ",errmin
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_poly_test
c-----------------------------------------------------------------------
c     subprogram 10. matrix_write1.
c     writes matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_write1(u)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: u

      CHARACTER(128) :: format1,format2
      INTEGER :: n,ipert,isol,ieq
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/1x,"ipert",2x,"ieq",',i2,'(4x,"re ",i2,2x)/)')
 20   FORMAT('(2i5,1p,',i3,'e11.3)')
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      n=SIZE(u,2)
c      n=MIN(n,9)
      WRITE(format1,10)n
      WRITE(format2,20)2*n
c-----------------------------------------------------------------------
c     output.
c-----------------------------------------------------------------------
      WRITE(out_unit,format1)(isol,isol=1,n)
      DO ieq=1,SIZE(u,3)
         DO ipert=1,SIZE(u,1)
            WRITE(out_unit,format2)ipert,ieq,REAL(u(ipert,1:n,ieq))
         ENDDO
      ENDDO
      WRITE(out_unit,format1)(isol,isol=1,n)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_write1
c-----------------------------------------------------------------------
c     subprogram 11. matrix_write2.
c     writes matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_write2(u)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: u

      CHARACTER(128) :: format1,format2
      INTEGER :: n,i,j
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/2x,"i",',i3,'(4x,"re ",i3,1x)/)')
 20   FORMAT('(i3,1p,',i3,'e11.3)')
c-----------------------------------------------------------------------
c     create formats.
c-----------------------------------------------------------------------
      n=SIZE(u,2)
c      n=MIN(n,9)
      WRITE(format1,10)n
      WRITE(format2,20)n
c-----------------------------------------------------------------------
c     output.
c-----------------------------------------------------------------------
      WRITE(out_unit,format1)(j,j=1,n)
      DO i=1,SIZE(u,1)
         WRITE(out_unit,format2)i,REAL(u(i,1:n))
      ENDDO
      WRITE(out_unit,format1)(j,j=1,n)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_write2
c-----------------------------------------------------------------------
c     subprogram 12. matrix_test_fixup.
c     tests fixup routine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_test_fixup

      INTEGER :: ifix,isol,msol,ipert
      REAL(r8), DIMENSION(:,:,:), POINTER :: re,im
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"isol",2x,"index"/)
 20   FORMAT(2i6)
 30   FORMAT(2(a,i3))
 40   FORMAT(/1x,"i",2(3x,"re u(",i1,")",4x,"im u(",i1,")",1x)/)
 50   FORMAT(i2,1p,4e11.3)
c-----------------------------------------------------------------------
c     create initial u.
c-----------------------------------------------------------------------
      ifix=mfix
      msol=fixtype(ifix)%msol
      ALLOCATE(u(mpert,msol,2),re(mpert,msol,2),im(mpert,msol,2))
      CALL put_seed(0)
      CALL RANDOM_NUMBER(re)
      CALL RANDOM_NUMBER(im)
      u=CMPLX(re,im)
      DEALLOCATE(re,im)
c-----------------------------------------------------------------------
c     diagnose index.
c-----------------------------------------------------------------------
      CALL ascii_open(debug_unit,"fixup_match.out","UNKNOWN")
      WRITE(debug_unit,10)
      DO isol=1,msol
         WRITE(debug_unit,20)isol,fixtype(ifix)%index(isol)
      ENDDO
      WRITE(debug_unit,10)
c-----------------------------------------------------------------------
c     diagnose initial u.
c-----------------------------------------------------------------------
      WRITE(debug_unit,'(/a/)')"initial values:"
      DO isol=1,msol
         WRITE(debug_unit,30)"isol = ",isol,", m = ",mlow+isol-1
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
         WRITE(debug_unit,50)(ipert,u(ipert,isol,:),ipert=1,mpert)
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
      ENDDO
c-----------------------------------------------------------------------
c     diagnose final u.
c-----------------------------------------------------------------------
      CALL matrix_fixup(u,ifix)
      WRITE(debug_unit,'(/a/)')"final values:"
      DO isol=1,msol
         WRITE(debug_unit,30)"isol = ",isol,", m = ",mlow+isol-1
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
         WRITE(debug_unit,50)(ipert,u(ipert,isol,:),ipert=1,mpert)
         WRITE(debug_unit,40)(ipert,ipert,ipert=1,2)
      ENDDO
      DEALLOCATE(u)
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_test_fixup
c-----------------------------------------------------------------------
c     subprogram 13. matrix_poly_scan.
c     graphs behavior of determinant as a function of deltas.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_poly_scan

      INTEGER :: ising,iscan
      REAL(r8) :: scan,dscan
      COMPLEX(r8) :: det0,det1,a,b,d
      COMPLEX(r8), DIMENSION(msing,2) :: deltar
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"i",5x,"scan",7x,
     $     "re a",7x,"im a",7x,"re b",7x,"im b",7x,"re d",7x,"im d"/)
 20   FORMAT(i3,1p,7e11.3)
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      ALLOCATE(ff(0:mterm))
      deltap1=0
      deltap2=0
      dscan=(scan1/scan0)**(1._r8/nscan)
c-----------------------------------------------------------------------
c     open files.
c-----------------------------------------------------------------------
      CALL ascii_open(debug_unit,"scan.out","UNKNOWN")
      CALL bin_open(bin_unit,"scan.bin","UNKNOWN","REWIND","none")
c-----------------------------------------------------------------------
c     perform scan.
c-----------------------------------------------------------------------
      DO ising=1,msing
         scan=scan0
         WRITE(debug_unit,'(2x,a,i2)')"ising = ",ising
         WRITE(debug_unit,10)
         DO iscan=0,nscan

            deltar=scan
            CALL matrix_make(deltar)
            CALL matrix_diag(diag0)

            deltar(ising,2)=delta0
            CALL matrix_make(deltar)
            CALL matrix_diag(diag)
            det0=PRODUCT(diag/diag0)

            deltar(ising,2)=delta1
            CALL matrix_make(deltar)
            CALL matrix_diag(diag)
            det1=PRODUCT(diag/diag0)

            a=scan*(det0*delta1-det1*delta0)/(delta1-delta0)
            b=scan*(det1-det0)/(delta1-delta0)
            d=-2*a/b

            WRITE(debug_unit,20)iscan,scan,a,b,d
            WRITE(bin_unit)REAL(LOG10(scan),4),
     $           REAL(REAL(a),4),REAL(REAL(b),4),REAL(REAL(d),4),
     $           REAL(AIMAG(a),4),REAL(AIMAG(b),4),REAL(AIMAG(d),4)

            scan=scan*dscan
         ENDDO
         deltap1(ising)=d
         WRITE(bin_unit)
         WRITE(debug_unit,10)
      ENDDO
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
      CALL bin_close(bin_unit)
      CALL ascii_close(debug_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_poly_scan
c-----------------------------------------------------------------------
c     subprogram 14. matrix_condense.
c     computes condensed matching matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_condense(m11)

      COMPLEX(r8), DIMENSION(2*msing,2*msing), INTENT(OUT) :: m11

      INTEGER :: info,ipmax,ipmin,ising,i,j,k,m2,mmin,msol,ir,jr,kr
      INTEGER, DIMENSION(mmatch) :: ipiv
      INTEGER, DIMENSION(2*msing) :: r1
      INTEGER, DIMENSION(mmatch-2*msing) :: r2
      COMPLEX(r8), DIMENSION(2*msing,mmatch-2*msing) :: m12
      COMPLEX(r8), DIMENSION(mmatch-2*msing,2*msing) :: m21
      COMPLEX(r8), DIMENSION(mmatch-2*msing,mmatch-2*msing) :: m22
c-----------------------------------------------------------------------
c     determine size of secondary matrix.
c-----------------------------------------------------------------------
      m2=mmatch-m1
      IF(m2 == 0)THEN
         m11=match
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     compute indices.
c-----------------------------------------------------------------------
      mmin=mhigh
      DO ising=1,msing
         mmin=MIN(mmin,NINT(nn*singtype(ising)%q))
      ENDDO
      ipmin=2*(mmin-mlow)+1
      ipmax=2*msing+ipmin-1
      r1=(/(ir,ir=ipmin,ipmax)/)
      r2=(/(jr,jr=1,ipmin-1),(kr,kr=ipmax+1,mmatch)/)
c-----------------------------------------------------------------------
c     distribute matching matrix.
c-----------------------------------------------------------------------
      i=0
      j=0
      k=0
      DO ising=1,msing
         msol=singtype(ising)%msol_r
         m11(j+1:j+2,:)=match(i+1:i+2,r1)
         m12(j+1:j+2,:)=match(i+1:i+2,r2)
         m21(k+1:k+msol-2,:)=match(i+3:i+msol,r1)
         m22(k+1:k+msol-2,:)=match(i+3:i+msol,r2)
         i=i+msol
         j=j+2
         k=k+msol-2
      ENDDO
      m21(k+1:m2,:)=match(i+1:mmatch,r1)
      m22(k+1:m2,:)=match(i+1:mmatch,r2)
c-----------------------------------------------------------------------
c     compute condensed matrix.
c-----------------------------------------------------------------------
      CALL zgetrf(m2,m2,m22,m2,ipiv,info)
      CALL zgetrs("N",m2,m1,m22,m2,ipiv,m21,m2,info)
      m11=m11-MATMUL(m12,m21)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_condense
      END MODULE matrix_mod
