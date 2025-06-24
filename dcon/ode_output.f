c-----------------------------------------------------------------------
c     file ode_output.f.
c     diagnoses main differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. ode_output_mod.
c     1. ode_output_open
c     2. ode_output_step
c     3. ode_output_close
c     4. ode_output_get_evals.
c     5. ode_output_monitor.
c     6. ode_output_get_crit.
c     7. ode_output_sol.
c-----------------------------------------------------------------------
c     subprogram 0. ode_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE ode_output_mod
      USE sing_mod
      USE dcon_mod, ONLY : shotnum,shottime
      IMPLICIT NONE

      CHARACTER(6), DIMENSION(:), POINTER :: name
      INTEGER, DIMENSION(:), POINTER :: sol_out_unit,sol_bin_unit
      
      INTEGER :: neq,ising,istep,m1,nzero
      REAL(r8) :: psifac,psizero,singfac,q,psi_save,psimax,singfac_old,
     $     psifac_old
      REAL(r8) :: singfac_min=1e-5,singfac_max=1e-4
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: u,u_old,u_save,du,
     $     ca,ca_old

      LOGICAL :: termbycross_flag=.FALSE.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. ode_output_open
c     integrates the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_open

      INTEGER :: isol,iunit
      REAL(r8), DIMENSION(:), POINTER :: xs,ys
      REAL(r8), DIMENSION(:,:,:), POINTER :: fs,fsx,fsy,fsxy
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"ising",3x,"psi",9x,"q",10x,"di",6x,"re alpha",
     $     3x,"im alpha"//i6,1p,5e11.3/)
 20   FORMAT(/3x,"is",4x,"psifac",6x,"dpsi",8x,"q",7x,"singfac",5x,
     $     "eval1"/)
 30   FORMAT(3x,"mlow",4x,"m",2x,"mhigh",2x,"isol",2x,"msol"//5i6/)
 40   FORMAT(1x,2x,a,es10.3,a,f6.3)
c-----------------------------------------------------------------------
c     allocate space for asymptotic coefficients.
c-----------------------------------------------------------------------
      nzero=0
      ALLOCATE(ca(mpert,msol,2),ca_old(mpert,msol,2))
c-----------------------------------------------------------------------
c     open file for solutions.
c-----------------------------------------------------------------------
      IF(bin_euler)THEN
         xs => rzphi%xs
         ys => rzphi%ys
         fs => rzphi%fs
         fsx => rzphi%fsx
         fsy => rzphi%fsy
         fsxy => rzphi%fsxy
         CALL bin_open(euler_bin_unit,"euler.bin","UNKNOWN","REWIND",
     $        "none")
         WRITE(euler_bin_unit)mlow,mhigh,nn,mpsi,mtheta,ro,zo
c-----------------------------------------------------------------------
c     <MODIFIED>
c-----------------------------------------------------------------------
         WRITE(euler_bin_unit)mband,mthsurf0,mthvac,psio,
     $        psilow,psilim,qlim,singfac_min
         WRITE(euler_bin_unit)power_b,power_r,power_bp
         WRITE(euler_bin_unit)kin_flag,con_flag
         ! from sum1 equilibrium descriptions
         WRITE(euler_bin_unit) amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,betap1,betap2,betap3,betat,betan,bt0,
     $     q0,qmin,qmax,qa,crnt,q95,shotnum,shottime
c-----------------------------------------------------------------------
         WRITE(euler_bin_unit)sq%xs,sq%fs,sq%fs1,sq%xpower
         WRITE(euler_bin_unit)rzphi%xs,rzphi%ys,
     $        rzphi%fs,rzphi%fsx,rzphi%fsy,rzphi%fsxy,
     $        rzphi%x0,rzphi%y0,rzphi%xpower,rzphi%ypower
      ENDIF
c-----------------------------------------------------------------------
c     open output files for crit and write header.
c-----------------------------------------------------------------------
      CALL ascii_open(crit_out_unit,"crit.out","UNKNOWN")
      CALL bin_open(crit_bin_unit,"crit.bin","UNKNOWN","REWIND","none")
      IF(out_evals)CALL ascii_open(evals_out_unit,"evals.out","UNKNOWN")
      IF(bin_evals)CALL bin_open(evals_bin_unit,"evals.bin","UNKNOWN",
     $     "REWIND","none")
      IF(ising > 0 .AND. ising <= msing)
     $     WRITE(crit_out_unit,10)ising,sing(ising)%psifac,
     $     sing(ising)%q,sing(ising)%di,sing(ising)%alpha
      WRITE(crit_out_unit,20)
c-----------------------------------------------------------------------
c     define file names.
c-----------------------------------------------------------------------
      IF(bin_sol)THEN
         ALLOCATE(name(2*mpert))
         DO isol=1,2*mpert
            WRITE(name(isol),'(a,i2.2)')"sol",isol
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     open binary output files for solutions.
c-----------------------------------------------------------------------
      IF(bin_sol)THEN
         bin_sol_min=1
         bin_sol_min=MAX(bin_sol_min,1)
         IF(res_flag)THEN
            bin_sol_max=MIN(bin_sol_max,mpert+2*msing,2*mpert)
         ELSE
            bin_sol_max=MIN(bin_sol_max,mpert)
         ENDIF
         IF(res_flag)THEN
            bin_sol_max=MIN(2*mpert,9)
         ELSE
            bin_sol_max=MIN(mpert,9)
         ENDIF
         bin_sol_max=MIN(bin_sol_max,98+bin_sol_min-sol_base)
         ALLOCATE(sol_bin_unit(bin_sol_min:bin_sol_max))
         sol_bin_unit=sol_base-bin_sol_min
     $        +(/(isol,isol=bin_sol_min,bin_sol_max)/)
         DO isol=bin_sol_min,bin_sol_max
            iunit=sol_bin_unit(isol)
            CALL bin_open(iunit,TRIM(name(isol))//".bin",
     $           "UNKNOWN","REWIND","none")
         ENDDO
         CALL bin_open(unorm_unit,"unorm.bin","UNKNOWN","REWIND","none")
      ENDIF
c-----------------------------------------------------------------------
c     print initial point of integration.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psifac,0)
      q=sq%f(4)
      IF(verbose) WRITE(*,40)"psi =",psifac,", q = ",q
c-----------------------------------------------------------------------
c     open file for error output.
c-----------------------------------------------------------------------
      IF(ksing > 0 .AND. ksing <= msing)THEN
         CALL bin_open(err_unit,"error.bin","UNKNOWN","REWIND","none")
         ALLOCATE(u_old(mpert,msol,2))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_open
c-----------------------------------------------------------------------
c     subprogram 2. ode_output_step
c     integrates the differential equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_step(unorm, op_force)

      REAL(r8), DIMENSION(:), INTENT(IN) :: unorm
      LOGICAL, OPTIONAL :: op_force
      LOGICAL :: force
      ! set optional parameters
      force = .FALSE.
      IF(PRESENT(op_force)) force = op_force
c-----------------------------------------------------------------------
c     compute and print critical data for each time step.
c-----------------------------------------------------------------------
      CALL ode_output_monitor
      IF(out_evals .OR. bin_evals)CALL ode_output_get_evals
c-----------------------------------------------------------------------
c     write solutions.
c-----------------------------------------------------------------------
      IF(bin_euler .AND. (mod(istep,euler_stride) == 0 .OR. force))THEN
         CALL sing_der(neq,psifac,u,du)
         WRITE(euler_bin_unit)1
         WRITE(euler_bin_unit)psifac,q,msol
         WRITE(euler_bin_unit)u
         WRITE(euler_bin_unit)ud
      ENDIF
c-----------------------------------------------------------------------
c     output solutions components for each time step.
c-----------------------------------------------------------------------
      IF(bin_sol)CALL ode_output_sol(psifac,u,unorm)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_step
c-----------------------------------------------------------------------
c     subprogram 3. ode_output_close.
c     produces ascii output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_close

      INTEGER :: isol
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
 10   FORMAT(1x,2x,a,es10.3,a,f6.3)
 20   FORMAT(/3x,"is",4x,"psifac",6x,"dpsi",8x,"q",7x,"singfac",5x,
     $     "eval1"/)
c-----------------------------------------------------------------------
c     close crit files.
c-----------------------------------------------------------------------
      WRITE(term_unit,10)"psi =",psifac,", q = ",q
      WRITE(crit_out_unit,20)
      CALL ascii_close(crit_out_unit)
      WRITE(crit_bin_unit)
      CALL bin_close(crit_bin_unit)
      DEALLOCATE(ca,ca_old)
c-----------------------------------------------------------------------
c     close other files.
c-----------------------------------------------------------------------
      IF(out_evals)CALL ascii_close(evals_out_unit)
      IF(bin_evals)CALL bin_close(evals_bin_unit)
      IF(ksing > 0 .AND. ksing <= msing)THEN
         CALL bin_close(err_unit)
         DEALLOCATE(u_old)
      ENDIF
c-----------------------------------------------------------------------
c     close binary output files.
c-----------------------------------------------------------------------
      IF(bin_sol)DEALLOCATE(name)
      IF(bin_sol)THEN
         DO isol=bin_sol_min,bin_sol_max
            CALL bin_close(sol_bin_unit(isol))
         ENDDO
         CALL bin_close(unorm_unit)
         DEALLOCATE(sol_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_close
c-----------------------------------------------------------------------
c     subprogram 4. ode_output_get_evals.
c     computes critical criterminant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_get_evals

      INTEGER :: info,ipert,lwork
      INTEGER, DIMENSION(mpert) :: ipiv,index,indexi

      REAL(r8) :: logpsi1,logpsi2
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      REAL(r8), DIMENSION(mpert) :: key,evals,evalsi

      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp,wp
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/a,i4,1p,2(a,es10.3))
 20   FORMAT(/5x,"i",5x,"eval",7x,"evali",7x,"err"/)
 30   FORMAT(i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     compute plasma response matrix.
c-----------------------------------------------------------------------
      temp=CONJG(TRANSPOSE(u(:,1:mpert,1)))
      wp=u(:,1:mpert,2)
      wp=CONJG(TRANSPOSE(wp))
      CALL zgetrf(mpert,mpert,temp,mpert,ipiv,info)  
      CALL zgetrs('N',mpert,mpert,temp,mpert,ipiv,wp,mpert,info)
      wp=(wp+CONJG(TRANSPOSE(wp)))/2
c-----------------------------------------------------------------------
c     compute and sort eigenvalues.
c-----------------------------------------------------------------------
      lwork=2*mpert-1  
      CALL zheev('N','U',mpert,wp,mpert,evals,work,lwork,rwork,info)
      index=(/(ipert,ipert=1,mpert)/)
      key=-ABS(1/evals)
      CALL bubble(key,index,1,mpert)
c-----------------------------------------------------------------------
c     compute and sort inverse eigenvalues.
c-----------------------------------------------------------------------
      lwork=2*mpert-1
      CALL zheev('N','U',mpert,wp,mpert,evalsi,work,lwork,rwork,info)
      indexi=(/(ipert,ipert=1,mpert)/)
      key=-ABS(evalsi)
      CALL bubble(key,indexi,1,mpert)
c-----------------------------------------------------------------------
c     write ascii eigenvalues.
c-----------------------------------------------------------------------
      IF(out_evals .AND. istep > 0)THEN
         WRITE(evals_out_unit,10)"istep = ",istep,
     $        ", psifac = ",psifac,", q = ",q
         WRITE(evals_out_unit,20)
         WRITE(evals_out_unit,30)(ipert,
     $        evalsi(indexi(ipert)),evals(index(ipert)),
     $        evalsi(indexi(ipert))*evals(index(ipert))-1,
     $        ipert=1,mpert)
         WRITE(evals_out_unit,20)
      ENDIF
c-----------------------------------------------------------------------
c     write binary eigenvalues.
c-----------------------------------------------------------------------
      IF(bin_evals .AND. psifac > .1)THEN
         CALL spline_eval(sq,psifac,0)
         q=sq%f(4)
         singfac=ABS(m1-nn*q)
         logpsi1=LOG10(psifac)
         logpsi2=LOG10(singfac)
         DO ipert=2,mpert
            WRITE(evals_bin_unit)REAL(ipert,4),REAL(psifac,4),
     $           REAL(logpsi1,4),REAL(logpsi2,4),
     $           REAL(q,4),REAL(1/evals(index(ipert)),4)
         ENDDO
         WRITE(evals_bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_get_evals
c-----------------------------------------------------------------------
c     subprogram 5. ode_output_monitor.
c     writes crit and monitors zero crossings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_monitor

      REAL(r8) :: logpsi1,logpsi2,crit,fac,
     $     psi_med,crit_med,q_med,singfac_med,logpsi1_med,logpsi2_med
      REAL(r8), SAVE :: crit_save=0,psi_save=0,dpsi
      COMPLEX(r8), DIMENSION(mpert,msol,2) :: u_med
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(i5,1p,6e11.3)
 20   FORMAT(1x,a,1p,e10.3,a,e10.3)
c-----------------------------------------------------------------------
c     compute new crit.
c-----------------------------------------------------------------------
      dpsi=psifac-psi_save
      CALL ode_output_get_crit(psifac,u,q,singfac,logpsi1,logpsi2,crit)
c-----------------------------------------------------------------------
c     check for zero crossing.
c-----------------------------------------------------------------------
      IF(crit*crit_save < 0)THEN
         fac=crit/(crit-crit_save)
         psi_med=psifac-fac*(psifac-psi_save)
         dpsi=psi_med-psi_save
         u_med=u-fac*(u-u_save)
         CALL ode_output_get_crit(psi_med,u_med,q_med,singfac_med,
     $        logpsi1_med,logpsi2_med,crit_med)
         IF((crit_med-crit)*(crit_med-crit_save) < 0 .AND.
     $        ABS(crit_med) < .5*MIN(ABS(crit),ABS(crit_save)))THEN
            WRITE(term_unit,20)
     $           "Zero crossing at psi =",psi_med,", q =",q_med
            WRITE(out_unit,20)
     $           "Zero crossing at psi =",psi_med,", q =",q_med
            WRITE(crit_out_unit,20)
     $           "Zero crossing at psi =",psi_med,", q =",q_med
            WRITE(crit_out_unit,10)istep,psi_med,dpsi,q_med,
     $           singfac_med,crit_med
            WRITE(crit_bin_unit)REAL(psi_med,4),REAL(logpsi1_med,4),
     $           REAL(logpsi2_med,4),REAL(q_med,4),REAL(crit_med,4)
            nzero=nzero+1
         ENDIF
         IF (termbycross_flag) THEN
            CALL program_stop("Terminated by zero crossing.")
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     write new crit.
c-----------------------------------------------------------------------
      WRITE(crit_out_unit,10)istep,psifac,dpsi,q,singfac,crit
      WRITE(crit_bin_unit)REAL(psifac,4),REAL(logpsi1,4),
     $     REAL(logpsi2,4),REAL(q,4),REAL(crit,4)
c-----------------------------------------------------------------------
c     update saved values.
c-----------------------------------------------------------------------
      psi_save=psifac
      crit_save=crit
      u_save=u
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_monitor
c-----------------------------------------------------------------------
c     subprogram 6. ode_output_get_crit.
c     computes critical criterminant.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_get_crit(psi,u,q,singfac,
     $     logpsi1,logpsi2,crit)

      REAL(r8), INTENT(IN) :: psi
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: u
      REAL(r8), INTENT(OUT) :: q,singfac,logpsi1,logpsi2,crit

      INTEGER :: info,ipert,lwork
      INTEGER, DIMENSION(mpert) :: ipiv,indexi
      REAL(r8), DIMENSION(3*mpert-2) :: rwork
      REAL(r8), DIMENSION(mpert) :: key,evalsi
      COMPLEX(r8), DIMENSION(mpert,mpert) :: temp,wp
      COMPLEX(r8), DIMENSION(2*mpert-1) :: work
      COMPLEX(r8), DIMENSION(mpert,mpert,2) :: uu
c-----------------------------------------------------------------------
c     compute dependent variables.
c-----------------------------------------------------------------------
      uu=u(:,1:mpert,:)
c-----------------------------------------------------------------------
c     compute inverse plasma response matrix.
c-----------------------------------------------------------------------
      wp=CONJG(TRANSPOSE(uu(:,:,1)))
      temp=uu(:,:,2)
      temp=CONJG(TRANSPOSE(temp))
      CALL zgetrf(mpert,mpert,temp,mpert,ipiv,info)  
      CALL zgetrs('N',mpert,mpert,temp,mpert,ipiv,wp,mpert,info)
      wp=(wp+CONJG(TRANSPOSE(wp)))/2
c-----------------------------------------------------------------------
c     compute and sort inverse eigenvalues.
c-----------------------------------------------------------------------
      lwork=2*mpert-1  
      CALL zheev('N','U',mpert,wp,mpert,evalsi,work,lwork,rwork,info)
      indexi=(/(ipert,ipert=1,mpert)/)
      key=-ABS(evalsi)
      CALL bubble(key,indexi,1,mpert)
c-----------------------------------------------------------------------
c     compute critical data for each time step.
c-----------------------------------------------------------------------
      CALL spline_eval(sq,psi,0)
      q=sq%f(4)
      singfac=ABS(m1-nn*q)
      logpsi1=LOG10(psi)
      logpsi2=LOG10(singfac)
      crit=evalsi(indexi(1))*sq%f(3)**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_get_crit
c-----------------------------------------------------------------------
c     subprogram 7. ode_output_sol.
c     produces binary output for each step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ode_output_sol(psifac,u,unorm)

      REAL(r8), INTENT(IN) :: psifac
      COMPLEX(r8), DIMENSION(mpert,msol,2), INTENT(IN) :: u
      REAL(r8), DIMENSION(:), INTENT(IN) :: unorm
      
      INTEGER :: isol,ipert,iunit,jsing
      REAL(r8) :: psil,psir,psiminl,psiminr,a,b,x,singlog
      REAL(r8), DIMENSION(mpert) :: logunorm
      COMPLEX(r8), DIMENSION(mpert,msol,2) :: uu,ca
      COMPLEX(r8), DIMENSION(mpert,2) :: logu,logca
c-----------------------------------------------------------------------
c     determine singular interval.
c-----------------------------------------------------------------------
      jsing=1
      DO
         IF(jsing > msing)EXIT
         IF(sing(jsing)%psifac > psifac)EXIT
         jsing=jsing+1
      ENDDO
c-----------------------------------------------------------------------
c     compute left edge of interval.
c-----------------------------------------------------------------------
      IF(jsing == 1)THEN
         psil=0
         psiminl=psilow
      ELSE
         psil=sing(jsing-1)%psifac
         psiminl=singfac_min/ABS(nn*sing(jsing-1)%q1)
      ENDIF
c-----------------------------------------------------------------------
c     compute right edge of interval.
c-----------------------------------------------------------------------
      IF(jsing > msing)THEN
         psiminr=psiminl
         psir=psilim+psiminr
      ELSE
         psir=sing(jsing)%psifac
         psiminr=singfac_min/ABS(nn*sing(jsing)%q1)
      ENDIF
c-----------------------------------------------------------------------
c     compute abscissa.
c-----------------------------------------------------------------------
      b=2/(LOG10((psir-psil-psiminl)/psiminl)
     $     +LOG10((psir-psil-psiminr)/psiminr))
      a=b*LOG10((psir-psil-psiminl)/psiminl)+2*jsing-2
      x=a+b*LOG10((psifac-psil)/(psir-psifac))
      singlog=LOG10(ABS(singfac))
c-----------------------------------------------------------------------
c     compute solutions and asymptotic coefficients.
c-----------------------------------------------------------------------
      uu=u
      IF(psil < psir .AND. jsing > 1 .AND. psifac-psil < psir-psifac
     $     .OR. jsing > msing)THEN
         CALL sing_get_ca(jsing-1,psifac,uu,ca)
      ELSE
         CALL sing_get_ca(jsing,psifac,uu,ca)
      ENDIF
c-----------------------------------------------------------------------
c     compute logu and logca
c-----------------------------------------------------------------------
      DO isol=MAX(1,bin_sol_min),MIN(bin_sol_max,msol)
         iunit=sol_bin_unit(isol)
         WHERE(uu(:,isol,:) /= 0)
            logu=LOG(uu(:,isol,:))
         ELSEWHERE
            logu=0
         END WHERE
         WHERE(ca(:,isol,:) /= 0)
            logca=LOG(ca(:,isol,:))
         ELSEWHERE
            logca=0
         END WHERE
c-----------------------------------------------------------------------
c     write binary output.
c-----------------------------------------------------------------------
         DO ipert=1,mpert
            WRITE(iunit)REAL(psifac,4),REAL(LOG10(psifac),4),REAL(x,4),
     $           REAL(DREAL(logu(ipert,:))/alog10,4),
     $           REAL(DREAL(logca(ipert,:))/alog10,4),
     $           REAL(LOG10(ABS((psifac-sing(1)%psifac)
     $           /(nn*sing(1)%q1)))),REAL(singlog,4)
         ENDDO
         WRITE(iunit)
      ENDDO
c-----------------------------------------------------------------------
c     diagnose unorm.
c-----------------------------------------------------------------------
      WHERE(unorm(1:mpert) > 1)
         logunorm=LOG10(unorm(1:mpert))
      ELSEWHERE
         logunorm=0
      ENDWHERE
      DO isol=1,mpert
         WRITE(unorm_unit)REAL(isol,4),REAL(psifac,4),
     $        REAL(LOG10(psifac),4),REAL(x,4),
     $        REAL(logunorm(isol),4)
      ENDDO
      WRITE(unorm_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ode_output_sol
      END MODULE ode_output_mod
