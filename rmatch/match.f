c-----------------------------------------------------------------------
c     file match.f.
c     match delta value between outer and inner regions;
c     find the mode root of growth rate.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. match_mod.
c     1. match_run.
c     2. match_init.
c     3. match_newton.
c     4. match_delta.
c     5. match_solution.
c     6. match_eta_scan.
c     7. match_qscan.
c     8. match_delta_jardin.
c     9. match_nyquist.
c     10. match_branch.
c     11. match_solve.
c     12. match_matrix_diagnose.
c     13. match_rpec.
c     14. match_alloc_sol.
c     15. match_dealloc_sol.
c     16. match_output_solution.
c     17. match_main.
c-----------------------------------------------------------------------
c     subprogram 0. match_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE match_mod

      USE fspline_mod
      USE deltar_mod
      USE innerc_module
      USE deltac_mod
      USE msing_mod
      IMPLICIT NONE
      
      TYPE :: branch_type
      CHARACTER(128) :: name
      COMPLEX(r8) :: s0,s1,ds,cycles
      END TYPE branch_type
      
      TYPE :: nyquist_type
      LOGICAL :: flag,out,bin
      INTEGER :: ns,narc
      REAL(r8) :: big,small
      END TYPE nyquist_type
     
      TYPE :: match_sol_type
      LOGICAL :: flag=.FALSE.
      LOGICAL :: uniform=.FALSE.
      LOGICAL :: auto_connect = .FALSE.
      LOGICAL :: b_flag = .FALSE.
      REAL(r8) :: qpert=2
      REAL(r8) :: connect_threshold=1e-3
      END TYPE match_sol_type 
      
      TYPE :: coil_type
         LOGICAL :: rpec_flag=.FALSE.
         LOGICAL :: ideal_flag=.FALSE.
         INTEGER :: mcoil,m1,m2      
      END TYPE coil_type
      
      TYPE :: outsol_type
         LOGICAL, DIMENSION(:),ALLOCATABLE :: issing
         INTEGER :: nsol,msing,mpert,mhigh,mlow,nn,tot_grids
         REAL(r8) :: psio
         REAL(r8), DIMENSION(:),ALLOCATABLE :: psi,q,qpsifac,qsing
         REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xext
         COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: sols,sols_cut
      END TYPE outsol_type
      
      TYPE :: insol_type
         INTEGER :: tot_g
         TYPE(solution_type), DIMENSION(:), POINTER :: sols
      END TYPE insol_type
      
      LOGICAL :: scan_flag=.FALSE.,sol_flag=.FALSE.,qscan_flag=.FALSE.,
     $     matrix_diagnose=.FALSE.
      LOGICAL :: qscan_out=.TRUE.,deltar_flag=.FALSE.,deflate=.FALSE.,
     $           deltac_flag=.FALSE.,deltaj_flag=.FALSE.,
     $           match_flag=.FALSE.
      CHARACTER(10) :: model="deltac"
      INTEGER :: msing,totmsing,nstep=32,scan_nstep,qscan_ising=1
      INTEGER :: nroot=1,iroot,totnsol,ising_output=1,itermax=500
      REAL(r8) :: eta(20),dlim=1000
      REAL(r8) :: scan_x0,scan_x1
      REAL(r8), DIMENSION(:), ALLOCATABLE :: taur_save
      REAL(r8), DIMENSION(:), ALLOCATABLE :: zo_out,zi_in
      COMPLEX(r8) :: initguess
      COMPLEX(r8), DIMENSION(:),ALLOCATABLE :: cofout,cofin,q_in
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: delta,deltar
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: deltaf
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: oldroots
            
      TYPE(resist_type),DIMENSION(:),POINTER :: restype
      TYPE(nyquist_type) :: nyquist
      TYPE(match_sol_type) :: match_sol
      TYPE(coil_type) :: coil
      TYPE(outsol_type),PRIVATE :: outs
      TYPE(insol_type),PRIVATE :: ins
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. match_run.
c     run delta matching between inner and outer regions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_run
      
      INTEGER :: ising,iter
      REAL(r8) :: rho,eta1,err
      COMPLEX(r8) :: eigval
      
      NAMELIST/match_input/ deltabin_filename,galsol_filename,
     $                         galsol_filename_cut,
     $                         initguess,msing,eta,sol_flag,
     $                         nstep,rtol,atol,fmin,fmax,lam,
     $                         scan_flag,scan_x0,scan_x1,scan_nstep,
     $                         model,qscan_ising,qscan_flag,qscan_out,
     $                         deltar_flag,deltac_flag,deltaj_flag,
     $                         deflate,nroot,match_flag,ising_output,
     $                         match_sol,matrix_diagnose,fulldomain,
     $                         coil,itermax 
      NAMELIST/nyquist_input/nyquist
10    FORMAT(1x,"Eigenvalue=",1p,2e11.3)
20    FORMAT(1x,"ising=",I2,1x,"q_in=",1p,2e11.3)
30    FORMAT(1x,"ising=",I2,1x,"zi=  ",1p,e11.3," zi*SQRT(10)=",e11.3)
40    FORMAT(1x,"ising=",I2,1x,"zo=  ",1p,e11.3," zo/10      =",e11.3)

c-----------------------------------------------------------------------
c     initialize and set parameters.
c-----------------------------------------------------------------------
      xvar=-1.0
      CALL ascii_open(out_unit,"match.out","REPLACE")
      CALL timer(0,out_unit)
      CALL match_init
      OPEN(UNIT=in_unit,FILE="match.in",STATUS="OLD")
      READ(in_unit,NML=match_input)
      REWIND(in_unit)
      READ(in_unit,NML=nyquist_input)
      REWIND(in_unit)
      CLOSE(UNIT=in_unit)      
      OPEN(UNIT=bin_unit,FILE=deltabin_filename,STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
      READ(bin_unit)totmsing,totnsol,coil%rpec_flag
      CALL deltac_read_parameters("match.in")
      ALLOCATE (delta(totnsol,2*totmsing))
      ALLOCATE (deltar(totmsing,2),restype(totmsing))
      ALLOCATE (deltaf(totmsing,2,2))
      ALLOCATE (taur_save(totmsing))
      ALLOCATE (zi_in(msing),zo_out(msing),q_in(msing))
      IF (totmsing.LT.msing) THEN
         WRITE(*,*)"msing is larger than totmsing."
         msing=totmsing
      ENDIF
      ALLOCATE (cofout(2*msing),cofin(2*msing))
      READ(bin_unit)delta
      DO ising=1,totmsing
         READ(bin_unit)
     $        restype(ising)%e,restype(ising)%f,
     $        restype(ising)%h,restype(ising)%m,
     $        restype(ising)%g,restype(ising)%k,
     $        eta1,rho,
     $        restype(ising)%taua,restype(ising)%taur,
     $        restype(ising)%v1
         taur_save(ising)=restype(ising)%taur*eta1
         restype(ising)%taur=taur_save(ising)/eta(ising)
         restype(ising)%ising=ising

         WRITE (*,*) 'ising=',ising
         WRITE (*,*) 'e=',restype(ising)%e,'f=',restype(ising)%f
         WRITE (*,*) 'h=',restype(ising)%h,'m=',restype(ising)%m
         WRITE (*,*) 'g=',restype(ising)%g,'k=',restype(ising)%k
         WRITE (*,*) 'taua=',restype(ising)%taua
         WRITE (*,*) 'taur=',restype(ising)%taur
         WRITE (*,*) 'v1=',restype(ising)%v1
      ENDDO
      CLOSE(UNIT=bin_unit)
      IF (fulldomain>0) THEN
         model="deltac"
      ENDIF
c-----------------------------------------------------------------------
c     msing of estimating zo.
c-----------------------------------------------------------------------
      CALL msing_init
      DO ising=1,msing
         CALL msing_estimate_zo(zo_out(ising),ising)
      ENDDO
c-----------------------------------------------------------------------
c     resistive perturbed equilibrium reconstruction.
c-----------------------------------------------------------------------
      IF (coil%rpec_flag) THEN
         coil%mcoil=totnsol-2*msing
         coil%m1=2*msing+1
         coil%m2=totnsol
         CALL match_rpec
         CALL program_stop("RPEC termination.")
      ENDIF      
c-----------------------------------------------------------------------
c     scan eigen value (Q) for different inner models.
c-----------------------------------------------------------------------
      IF(qscan_flag) CALL match_qscan        
c-----------------------------------------------------------------------
c     nyquist plot.
c-----------------------------------------------------------------------
      IF(nyquist%flag)CALL match_nyquist
c-----------------------------------------------------------------------
c     scan over resistivity.
c-----------------------------------------------------------------------
      IF (scan_flag) CALL match_eta_scan
c-----------------------------------------------------------------------
c     construct the global solution in both outer and inner regions.
c-----------------------------------------------------------------------
      IF (match_flag) THEN
         nroot=1
         eigval=initguess
         CALL match_newton(match_delta,eigval,err,iter)
         WRITE (*,10) eigval
         WRITE(out_unit,10) eigval
         DO ising=1,msing
            WRITE(*,20) ising,q_in(ising)
            WRITE(out_unit,20) ising,q_in(ising)
         ENDDO
         DO ising=1,msing
            WRITE(*,30) ising,zi_in(ising),zi_in(ising)*SQRT(10.0)
            WRITE(out_unit,30)ising,zi_in(ising),zi_in(ising)*SQRT(10.0)
         ENDDO         
c         CALL match_solution(eigval)
c         DO ising=1,msing
c            WRITE(*,40) ising,zo_out(ising),zo_out(ising)/10
c            WRITE(out_unit,40) ising,zo_out(ising),zo_out(ising)/10
c         ENDDO
         CALL ascii_close(match_unit)
         CALL program_stop("Normal termination for solution match.")
      ENDIF
c-----------------------------------------------------------------------
c     solve and output delta.
c-----------------------------------------------------------------------
      IF (sol_flag) THEN
         IF (.NOT.deflate) nroot=1
         CALL match_solve(eigval)
         CALL ascii_close(match_unit)
      ENDIF
c-----------------------------------------------------------------------
      DEALLOCATE (delta,deltar,deltaf,restype,cofin,cofout,taur_save)
      DEALLOCATE (zo_out,zi_in,q_in)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END SUBROUTINE match_run      
c-----------------------------------------------------------------------
c     subprogram 2. match_init.
c     initialize the parameters.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------         
      SUBROUTINE match_init
      initguess=0.0
      msing=-1
      eta=0.0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE match_init
c-----------------------------------------------------------------------
c     subprogram 3. match_newton.
c     uses Newton's method to find root of scalar complex f(z).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_newton(ff,z,err,it)

      COMPLEX(r8) :: ff
      COMPLEX(r8), INTENT(INOUT) :: z
      REAL(r8), INTENT(OUT) :: err
      INTEGER, INTENT(OUT) :: it

      INTEGER :: ising,info,nmat
      INTEGER, DIMENSION(4*msing-1) :: ipiv
      REAL(r8) :: dzfac=0.05, tol=1e-10, itmax=1000,relaxfac=0.1
      COMPLEX(r8) :: z_old,f_old,f,dz
      COMPLEX(r8), DIMENSION(4*msing) :: cof
      COMPLEX(r8), DIMENSION(4*msing,4*msing) :: mat
      COMPLEX(r8), DIMENSION(4*msing-1,4*msing-1) :: cmat
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
10    FORMAT("it",9x,"err",7x,"re z",7x,"im z",7x,"re f",7x,"im f",$)
20    FORMAT(i2,1x,5(1p,e11.3),$)
30    FORMAT(1x,"cof_out(",i2,") CL=",1p,2e11.3,"  CR=",1p,2e11.3)
40    FORMAT(1x,"cof_in(",i2,")  d+=",1p,2e11.3,"  d-=",1p,2e11.3)
      itmax=itermax
c-----------------------------------------------------------------------
c     find initial guess.
c-----------------------------------------------------------------------
      f=ff(z,mat)
      dz=z*dzfac
      it=0
c-----------------------------------------------------------------------
c     iterate.
c-----------------------------------------------------------------------
      DO
         it=it+1
         err=ABS(dz/z)
         IF (it==1) THEN
            WRITE (out_unit,10)
            WRITE(out_unit,*)
            WRITE(out_unit,*)            
         ENDIF
         WRITE(out_unit,20)it,err,REAL(z),AIMAG(z),REAL(f),AIMAG(f)
         WRITE(out_unit,*)

         IF(err < tol) EXIT
         IF(it > itmax) THEN
            it=-1
            WRITE(*,*) "Solution is not well converged."
            WRITE(out_unit,*) "Solution is not well converged."
            EXIT
         ENDIF
         z_old=z
         z=z+dz*relaxfac
         f_old=f
         f=ff(z,mat)
         dz=-f*(z-z_old)/(f-f_old)
      ENDDO
c-----------------------------------------------------------------------
c     compute the coefficients of outter and inner region solutions.
c-----------------------------------------------------------------------
      IF (match_flag) THEN
         cmat=mat(2:4*msing,2:4*msing)
         cof=0
         cof(1)=1
         nmat=4*msing
         cof(2:4*msing)=-mat(2:4*msing,1)
         CALL zgetrf(nmat-1,nmat-1,cmat,nmat-1,ipiv,info)
         CALL zgetrs('N',nmat-1,1,cmat,nmat-1,ipiv,
     $        cof(2:nmat),nmat-1,info)          
         IF(matrix_diagnose)CALL match_matrix_diagnose(mat,cof)
         DO ising=1,msing
            WRITE (*,30) ising,cof(2*ising-1),cof(2*ising)
            WRITE (*,40) ising,cof(2*msing+2*ising-1),
     $                         cof(2*msing+2*ising)
         ENDDO
         cofout=cof(1:2*msing)
         cofin=cof(2*msing+1:4*msing)
         
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_newton
c-----------------------------------------------------------------------
c     subprogram 4. match_delta.
c     get determinant of matching matrix (dispersion relation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION match_delta(guess,mat) RESULT(det)

      COMPLEX(r8), INTENT(IN) :: guess
      COMPLEX(r8), DIMENSION(4*msing,4*msing),INTENT(INOUT) :: mat
      COMPLEX(r8):: det

      INTEGER :: m,info,i,d,ising,idx1,idx2,idx3,idx4
      COMPLEX(r8) :: delta1,delta2,drl,drr,dll,dlr
      COMPLEX(r8), DIMENSION(4,2) :: sol
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: mattmp
      
      ALLOCATE (mattmp(4*msing,4*msing))
      ALLOCATE (ipiv(4*msing))
c-----------------------------------------------------------------------
c     preliminary computation.
c-----------------------------------------------------------------------
      mat=0
      mat(2*msing+1:4*msing,1:2*msing)
     $     =TRANSPOSE(delta(1:2*msing,1:2*msing))
      deltar=0
c-----------------------------------------------------------------------
c     start loop over singular surfaces.
c-----------------------------------------------------------------------
      DO ising=1,msing
         idx1=ising*2-1
         idx2=ising*2
         idx3=idx1+2*msing
         idx4=idx2+2*msing
c-----------------------------------------------------------------------
c     compute inner region matching data.
c-----------------------------------------------------------------------
         SELECT CASE(model)
         CASE ("deltaj")
            CALL match_delta_jardin(restype(ising),guess,
     $           deltar(ising,:),sol)     
         CASE ("deltar")
            CALL deltar_run(restype(ising),guess,deltar(ising,:),sol)
         CASE ("deltac")
            CALL deltac_run(restype(ising),guess,deltar(ising,:),
     $                      deltaf(ising,:,:))
            zi_in(ising)=zi_deltac
            q_in(ising)=q_deltac
            sol=0
         END SELECT
c-----------------------------------------------------------------------
c     construct the matching matrix.
c-----------------------------------------------------------------------
         SELECT CASE(fulldomain)
         CASE (0)
            delta1=deltar(ising,1)
            delta2=deltar(ising,2)
            mat(idx1,idx1)=1
            mat(idx2,idx2)=1
            mat(idx1,idx3)=-1
            mat(idx1,idx4)=1
            mat(idx2,idx3)=-1
            mat(idx2,idx4)=-1
            mat(idx3,idx3)=-delta1
            mat(idx3,idx4)=delta2
            mat(idx4,idx3)=-delta1
            mat(idx4,idx4)=-delta2
         CASE(1,2)
            drl=deltaf(ising,1,1)
            drr=deltaf(ising,1,2)
            dll=deltaf(ising,2,1)
            dlr=deltaf(ising,2,2)
            mat(idx1,idx1)=1
            mat(idx2,idx2)=1
            mat(idx1,idx3)=-1
            mat(idx2,idx4)=-1
            mat(idx3,idx3)=-dll
            mat(idx3,idx4)=-drl
            mat(idx4,idx3)=-dlr
            mat(idx4,idx4)=-drr
         END SELECT
c-----------------------------------------------------------------------
c     finish loop over singular surfaces.
c-----------------------------------------------------------------------
      ENDDO
c-----------------------------------------------------------------------
c     compute determinant of matching matrix.
c-----------------------------------------------------------------------
      m=4*msing
      info=-1
      mattmp=mat
      ipiv=0
      CALL zgetrf(m,m,mattmp,m,ipiv,info )
      d=1
      det=1.0
      DO i=1,m
         IF (ipiv(i).ne.i) d=-d
         det=det*mattmp(i,i)
      ENDDO
      det=det*d
c-----------------------------------------------------------------------
c     deflation: remove old roots.
c-----------------------------------------------------------------------
      IF(deflate .AND. iroot > 1)
     $     det=det/PRODUCT(guess-oldroots(1:iroot-1))
      DEALLOCATE (mattmp)
      DEALLOCATE (ipiv)
    
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION match_delta      
c-----------------------------------------------------------------------
c     subprogram 5. match_solution.
c     construct global solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_solution(eig)
      COMPLEX(r8), INTENT(IN) :: eig

      INTEGER :: galnsol,msing,tot_grids,isol,mhigh,mlow,ip,eff_grids
      INTEGER :: ipert,mpert,nn,ising,tot_g,jsing
      CHARACTER(100) :: comp_tittle,tmp
      REAL(r8),PARAMETER :: eps=1e-4
      REAL(r8) :: psio,inpsi

      LOGICAL, DIMENSION(:), ALLOCATABLE :: issing
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psi,q,singfac,psi_cut
      REAL(r8), DIMENSION(:), ALLOCATABLE :: qpsifac,qsing
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: inpsifac,xext
      COMPLEX(r8) :: insol
      COMPLEX(r8), DIMENSION(2) :: deltai
      COMPLEX(r8),DIMENSION(:,:), ALLOCATABLE :: outtotsol,intotsol
      COMPLEX(r8),DIMENSION(:,:), ALLOCATABLE :: outtotsol_cut,tmp_cut
      COMPLEX(r8),DIMENSION(:,:,:), ALLOCATABLE  :: outsol,outsol_cut
      TYPE(solution_type), DIMENSION(:),TARGET,ALLOCATABLE :: sols
      TYPE(cspline_type) :: outcut_sp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
10    FORMAT (1P,A15,$)
20    FORMAT (1P,E15.5,$)
30    FORMAT(/4x,"is",5x,"xvar",6x,"re x1",6x,"im x1",
     $     6x,"re x2",6x,"im x2"/)
40    FORMAT(i6,1p,5e11.3)
50    FORMAT(1x,"delta+=",1p,e11.3,1x"+",e11.3,1x,"i")
60    FORMAT(1x,"delta-=",1p,e11.3,1x"+",e11.3,1x,"i")
c-----------------------------------------------------------------------
c     read outer region solutions.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,galsol_filename,"OLD","REWIND","none")
      READ(bin_unit) galnsol,msing,mpert,mhigh,mlow,nn,psio,tot_grids
      ALLOCATE(psi(0:tot_grids),issing(0:tot_grids),q(0:tot_grids))
      ALLOCATE(qpsifac(msing),qsing(msing))
      ALLOCATE(singfac(mpert),outsol(mpert,0:tot_grids,galnsol),
     $         outtotsol(mpert,0:tot_grids))
      ALLOCATE(outsol_cut(mpert,0:tot_grids,galnsol),
     $         outtotsol_cut(mpert,0:tot_grids),
     $         tmp_cut(mpert,0:tot_grids),psi_cut(0:tot_grids))
      ALLOCATE(xext(msing,2))
      DO isol=1,msing
         READ(bin_unit) qpsifac(isol),qsing(isol)
      ENDDO
      READ(bin_unit) psi,issing,q
      DO isol=1,galnsol
         READ (bin_unit) outsol(:,:,isol)
      ENDDO
      CALL bin_close(bin_unit)

      CALL bin_open(bin_unit,galsol_filename_cut,"OLD",
     $                 "REWIND","none")
      READ (bin_unit) xext
      DO isol=1,galnsol
         READ (bin_unit) outsol_cut(:,:,isol)
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     construct linear combination of outer region solutions with cofout.
c-----------------------------------------------------------------------
      outtotsol=0
      outtotsol_cut=0
      DO isol=1,2*msing
         outtotsol=outtotsol+cofout(isol)*outsol(:,:,isol)
         outtotsol_cut=outtotsol_cut+cofout(isol)*outsol_cut(:,:,isol)
      ENDDO
c-----------------------------------------------------------------------
c     construct inner region solutions for each singular surface.
c-----------------------------------------------------------------------
      ALLOCATE(sols(msing))
      output_sol=.true.
      tot_g=interp_np*nx
      ALLOCATE(intotsol(-tot_g:tot_g,msing),
     $     inpsifac(-tot_g:tot_g,msing))
      DO ising=1,msing
         sol => sols(ising)
         ALLOCATE (sol%xvar(0:tot_g),sol%sol(3,0:tot_g,2))
         CALL deltac_run(restype(ising),eig,deltai,deltaf(ising,:,:))
c-----------------------------------------------------------------------
c    combine inner layer solutions
c-----------------------------------------------------------------------           
         DO ip=-tot_g,-1
            inpsifac(ip,ising)=qpsifac(ising)-sols(ising)%xvar(-ip)
            intotsol(ip,ising)=-sol%sol(2,-ip,1)*cofin(2*ising)
     $                         +sol%sol(2,-ip,2)*cofin(2*ising-1)
         ENDDO
         DO ip=0,tot_g
            inpsifac(ip,ising)=qpsifac(ising)+sols(ising)%xvar(ip)
            intotsol(ip,ising)=sol%sol(2,ip,1)*cofin(2*ising)
     $                        +sol%sol(2,ip,2)*cofin(2*ising-1)
         ENDDO
      ENDDO
!c-----------------------------------------------------------------------
!c     convert xi to b field.
!c-----------------------------------------------------------------------      
!      DO ip=0,tot_grids
!            singfac=mlow-nn*q(ip)+(/(ipert,ipert=0,mpert-1)/)
!            outtotsol(:,ip)=twopi*ifac*psio*singfac*outtotsol(:,ip)
!      ENDDO
c-----------------------------------------------------------------------
c     write full outer region solutions, ascii.
c-----------------------------------------------------------------------
      CALL ascii_open(match_unit,"outsol.out","REPLACE")
      WRITE (match_unit,10) 'psifac'
      DO ipert=mlow,mhigh
         WRITE (tmp,"(I)") ipert
         tmp=ADJUSTL(tmp)
         WRITE (comp_tittle,*) 'REAL(',TRIM(tmp),')'
         WRITE (match_unit,10) TRIM(comp_tittle)
         WRITE (comp_tittle,*) 'IMAG(',TRIM(tmp),')'
         WRITE (match_unit,10) TRIM(comp_tittle)
      ENDDO
      WRITE (match_unit,*)
      DO ip=0,tot_grids
         IF (issing(ip)) CYCLE
         WRITE (match_unit,20) psi(ip)
         DO ipert=1,mpert
            WRITE (match_unit,20)
     $            REAL(outtotsol(ipert,ip)),
     $            IMAG(outtotsol(ipert,ip))
         ENDDO
         WRITE (match_unit,*)
      ENDDO
      CALL ascii_close(match_unit)

c-----------------------------------------------------------------------
c     write full outer region solutions, binary.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,'outsol_tot.bin',"UNKNOWN","REWIND","none")
      DO ipert=1,mpert
         DO ip=0,tot_grids                
            IF (issing(ip)) THEN
               WRITE(bin_unit)
               CYCLE
            ENDIF
            WRITE (bin_unit) REAL(psi(ip),4),
     $                       REAL(outtotsol(ipert,ip),4),
     $                       REAL(IMAG(outtotsol(ipert,ip)),4),
     $                       mylog(outtotsol(ipert,ip))
         ENDDO
         WRITE(bin_unit)
      ENDDO  
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     write inner region solutions, binary.
c-----------------------------------------------------------------------      
      tmp_cut=outtotsol_cut
      outtotsol_cut=0
      eff_grids=-1
      DO ip=0,tot_grids
         IF (issing(ip)) THEN
            CYCLE
         ENDIF
         eff_grids=eff_grids+1
         outtotsol_cut(:,eff_grids)=tmp_cut(:,ip)
         psi_cut(eff_grids)=psi(ip)
      ENDDO
      CALL cspline_alloc(outcut_sp,eff_grids,mpert)
      outcut_sp%xs(0:eff_grids)=psi_cut(0:eff_grids)
      DO ip=0,eff_grids
         outcut_sp%fs(ip,:)=outtotsol_cut(:,ip)
      ENDDO
      CALL cspline_fit(outcut_sp,"extrap")
      CALL bin_open(bin_unit,'insol.bin',"UNKNOWN",
     $                 "REWIND","none")
	  CALL ascii_open(match_unit,"xi_in.out","REPLACE")
      DO ising=1,msing
         DO ip=-tot_g,tot_g
            IF (match_sol%uniform) THEN
               inpsi=inpsifac(ip,ising)
               DO jsing=1,msing
                  IF (xext(ising,1)<inpsi .AND. inpsi<xext(ising,2))
     $            THEN 
                     ipert=NINT(nn*qsing(ising))-mlow+1
                     CALL cspline_eval(outcut_sp,inpsi,0)
                     insol=outcut_sp%f(ipert)+intotsol(ip,ising)
                     WRITE (bin_unit) REAL(inpsifac(ip,ising),4),
     $                                REAL(insol,4),
     $                                REAL(IMAG(insol),4),
     $                                mylog(insol)
	                 WRITE (match_unit,20) inpsifac(ip,ising),
     $                                REAL(insol),
     $                                IMAG(insol)
	                 WRITE (match_unit,*)
                  ENDIF
               ENDDO
            ELSE
               insol=intotsol(ip,ising)
               WRITE (bin_unit) REAL(inpsifac(ip,ising),4),
     $                          REAL(insol,4),
     $                          REAL(IMAG(insol),4),
     $                          mylog(insol)
	           WRITE (match_unit,20) inpsifac(ip,ising),
     $                          REAL(insol),
     $                          IMAG(insol)
	           WRITE (match_unit,*)
            ENDIF
         ENDDO
         WRITE(bin_unit)
      ENDDO 
      CALL bin_close(bin_unit)
	  CALL ascii_close(match_unit)
c-----------------------------------------------------------------------
c     write resonant outer region solutions, binary.
c-----------------------------------------------------------------------
      IF (match_sol%flag) THEN
         CALL bin_open(bin_unit,'outsol_qpert.bin',"UNKNOWN",
     $                 "REWIND","none")
         
         DO ip=0,tot_grids                
            IF (issing(ip)) THEN
               WRITE(bin_unit)
               CYCLE
            ENDIF
            ipert=NINT(nn*match_sol%qpert)-mlow+1
            WRITE (bin_unit) REAL(psi(ip),4),
     $                       REAL(outtotsol(ipert,ip),4),
     $                       REAL(IMAG(outtotsol(ipert,ip)),4),
     $                       mylog(outtotsol(ipert,ip))
         ENDDO
         WRITE(bin_unit)
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write deltas to delta.out, ascii.
c-----------------------------------------------------------------------
      output_sol=.false.
      CALL ascii_open(match_unit,"delta.out","REPLACE")
      WRITE(match_unit,50)deltai(1)
      WRITE(match_unit,60)deltai(2)
      CALL ascii_close(match_unit)
c-----------------------------------------------------------------------
c     diagnostic of u_out.
c-----------------------------------------------------------------------
      CALL msing_run(cofout,delta)
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(psi,issing,q,qpsifac,qsing)
      DEALLOCATE(singfac,outsol,outtotsol)
      DEALLOCATE(outsol_cut,outtotsol_cut,tmp_cut,psi_cut)
      DEALLOCATE(intotsol,inpsifac)
      DO ising=1,msing
         DEALLOCATE(sols(ising)%xvar,sols(ising)%sol)
      ENDDO
      DEALLOCATE (sols)
      DEALLOCATE (xext)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_solution      
c-----------------------------------------------------------------------
c     subprogram 6. match_eta_scan.
c     scan parameter.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------            
      SUBROUTINE match_eta_scan
      INTEGER :: istep,ising,iter
      REAL(r8) :: step,eta_scan,log_scan_x0,err
      COMPLEX(r8) :: eigval
c-----------------------------------------------------------------------
c     format output.
c-----------------------------------------------------------------------                  
10    FORMAT(/12x,"eta",10x,"re_gr",10x,"im_gr",2x,"iter",3x,"ising",
     $       13x,"zi",13x,"zo",4x,"zi*SQRT(10)",10x,"zo/10",
     $       7x,"re(q_in)",7x,"im(q_in)"/)
20    FORMAT(1p,3e15.5,i6,i8,8e15.5)
c-----------------------------------------------------------------------
c     scan constant eta parameter.
c-----------------------------------------------------------------------                  
      log_scan_x0=log10(scan_x0)
      step=(log10(scan_x1)-log_scan_x0)/scan_nstep
      eigval=initguess
      CALL bin_open(bin_unit,"scanres.bin","UNKNOWN","REWIND","none")
      CALL ascii_open(debug_unit,"scanres.out","UNKNOWN")
      WRITE (debug_unit,10)
      DO istep=0,scan_nstep
         eta_scan=10**(log_scan_x0+istep*step)
         DO ising=1,totmsing
            restype(ising)%taur=taur_save(ising)
     $                         /(eta_scan*eta(ising)/eta(1))
         ENDDO
         CALL match_newton(match_delta,eigval,err,iter)
         ising=ising_output
         WRITE(bin_unit)REAL(eta_scan,4),REAL(log10(eta_scan),4),
     $      REAL(eigval,4),REAL(AIMAG(eigval),4),
     $      mylog(eigval),
     $      REAL(zi_in(ising),4),REAL(zo_out(ising),4),
     $      REAL(zi_in(ising)*SQRT(10.0),4),REAL(zo_out(ising)/10,4),
     $      REAL(q_in(ising),4),REAL(AIMAG(q_in(ising)),4),
     $      mylog(q_in(ising)),REAL(iter,4)
         WRITE(debug_unit,20)eta_scan,
     $    REAL(eigval),IMAG(eigval),iter,ising,zi_in(ising),
     $    zo_out(ising),zi_in(ising)*SQRT(10.0),zo_out(ising)/10,
     $    REAL(q_in(ising)),IMAG(q_in(ising))
      ENDDO
      WRITE(bin_unit)
      CALL ascii_close(debug_unit)
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination for eta scan.")
      END SUBROUTINE match_eta_scan
c-----------------------------------------------------------------------
c     subprogram 7. match_qscan.
c     scan q for different inner layer model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------            
      SUBROUTINE match_qscan
      INTEGER :: istep
      REAL(r8) :: log_scan_x0,step,qlog
      COMPLEX(r8) :: q_scan
      COMPLEX(r8), DIMENSION(2) :: deltac,deltar,deltaj
      COMPLEX(r8), DIMENSION(2,2) :: df
      COMPLEX(r8), DIMENSION(4,2) :: sol
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(deltar_flag)THEN
         OPEN(UNIT=bin1_unit,FILE="deltar.bin",STATUS="REPLACE",
     $        FORM="UNFORMATTED")
         CALL ascii_open(out1_unit,"deltar.out","REPLACE")
      ENDIF
      IF(deltac_flag)THEN
         OPEN(UNIT=bin2_unit,FILE="deltac.bin",STATUS="REPLACE",
     $        FORM="UNFORMATTED")
         CALL ascii_open(out2_unit,"deltac.out","REPLACE")
      ENDIF
      IF(deltaj_flag)THEN
         OPEN(UNIT=bin3_unit,FILE="deltaj.bin",STATUS="REPLACE",
     $        FORM="UNFORMATTED")
         CALL ascii_open(out3_unit,"deltaj.out","REPLACE")
      ENDIF
10    FORMAT(1x,3(1p,e11.3))      
c-----------------------------------------------------------------------
c     start loops over Q.
c-----------------------------------------------------------------------
      log_scan_x0=log10(scan_x0)
      step=(log10(scan_x1)-log_scan_x0)/scan_nstep
      DO istep=0,scan_nstep
         qlog=log_scan_x0+istep*step
         q_scan=10**(qlog)
c-----------------------------------------------------------------------
c     run deltar code and record output.
c-----------------------------------------------------------------------
         IF(deltar_flag)THEN
            CALL deltar_run(restype(qscan_ising),q_scan,deltar,sol)
            IF(qscan_out .AND. ABS(deltar(2)) < dlim .AND.
     $      ISNAN(ABS(deltar(1)))==.FALSE. .AND. 
     $      ISNAN(ABS(deltar(2)))==.FALSE.) THEN
                 WRITE(out1_unit,10)REAL(qlog),
     $           mylog(deltar(1)),REAL(deltar(2))
                 WRITE(bin1_unit)REAL(qlog,4),
     $           mylog(deltar(1)),REAL(deltar(2),4)
            ENDIF 
         ENDIF
c-----------------------------------------------------------------------
c     run deltac code and record output.
c-----------------------------------------------------------------------
         IF(deltac_flag)THEN
           CALL deltac_run(restype(qscan_ising),q_scan,deltac,df)
            IF(qscan_out .AND. ABS(deltac(2)) < dlim .AND.
     $      ISNAN(ABS(deltac(1)))==.FALSE. .AND.
     $      ISNAN(ABS(deltac(2)))==.FALSE.) THEN
                WRITE(out2_unit,10)REAL(qlog),
     $           mylog(deltac(1)),REAL(deltac(2))
                WRITE(bin2_unit)REAL(qlog,4),
     $          mylog(deltac(1)),REAL(deltac(2),4)
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     run inner code and record output.
c-----------------------------------------------------------------------
         IF(deltaj_flag)THEN
            CALL match_delta_jardin(restype(qscan_ising),q_scan,
     $                               deltaj,sol)
            IF(qscan_out .AND. ABS(deltaj(2)) < dlim .AND.
     $      ISNAN(ABS(deltaj(1)))==.FALSE. .AND.
     $      ISNAN(ABS(deltaj(2)))==.FALSE.) THEN
               WRITE(out3_unit,10)REAL(qlog),
     $           mylog(deltaj(1)),REAL(deltaj(2))
               WRITE(bin3_unit)REAL(qlog,4),
     $         mylog(deltaj(1)),REAL(deltaj(2),4)
            ENDIF
         ENDIF
      ENDDO
         IF(deltar_flag)WRITE(bin1_unit)
         IF(deltac_flag)WRITE(bin2_unit)
         IF(deltaj_flag)WRITE(bin3_unit)
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(deltar_flag)THEN
         CLOSE(UNIT=out1_unit)
         CLOSE(UNIT=bin1_unit)
      ENDIF
      IF(deltac_flag)THEN
         CLOSE(UNIT=out2_unit)
         CLOSE(UNIT=bin2_unit)
      ENDIF
      IF(deltaj_flag)THEN
         CLOSE(UNIT=out3_unit)
         CLOSE(UNIT=bin3_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination for q scan.")
      END SUBROUTINE match_qscan
      
c-----------------------------------------------------------------------
c     subprogram 8. match_delta_jardin.
c     finite differential method of GGJ.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------            
      SUBROUTINE match_delta_jardin(restype,s,deltar,sol)
      
      TYPE(resist_type), INTENT(IN) :: restype
      COMPLEX(r8), INTENT(IN) :: s
      COMPLEX(r8), DIMENSION(2), INTENT(OUT) :: deltar
      COMPLEX(r8), DIMENSION(4,2), INTENT(INOUT) :: sol
      
      INTEGER :: ifail
      REAL(r8) :: sfac,x0,q0,v1,taur,taua,dr,di,p1,ee,ff,hh
      REAL(8) :: epsd,rmatch
      COMPLEX(8) :: q,deltae,deltao

      ifail=0
      rmatch=3
      epsd=1e-4
      sol=0

      ee=restype%e   
      ff=restype%f
      hh=restype%h
      v1=restype%v1
      taur=restype%taur
      taua=restype%taua
      sfac=taur/taua
      x0=sfac**(-1._r8/3._r8)
      q0=x0/taua
      q=s/q0
      
      dr=ee+ff+hh*hh
      di=dr-(hh-.5)**2
      p1=SQRT(-di)
      
      CALL innerc(q,restype%e,restype%f,restype%g,
     $           restype%h,restype%k,deltae,deltao,ifail,epsd,rmatch)
      deltar(1)=deltao
      deltar(2)=deltae
      deltar=deltar*sfac**(2*p1/3)*v1**(2*p1)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------      
      END SUBROUTINE match_delta_jardin
c-----------------------------------------------------------------------
c     subprogram 9. match_nyquist.
c     draws nyquist plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_nyquist

      INTEGER :: iarc,roots,narc
      REAL(r8) :: sarc,dsarc,error
      COMPLEX(r8) :: cycles
      TYPE(branch_type), DIMENSION(:), POINTER ::
     $     pbranch,nbranch,cbranch
      TYPE(branch_type), POINTER :: bp
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"iarc",1x,"roots",4x,"error",7x,"|s|"/)
 20   FORMAT(2i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      WRITE(*,*)"Construct Nyquist Plot"
      narc=nyquist%narc
      ALLOCATE(pbranch(narc),nbranch(narc),cbranch(0:narc))
c-----------------------------------------------------------------------
c     set up positive radial branches.
c-----------------------------------------------------------------------
      sarc=nyquist%big
      dsarc=(nyquist%small/nyquist%big)**(one/nyquist%narc)
      DO iarc=1,narc
         bp => pbranch(iarc)
         bp%s0=sarc*EXP(ifac*pi/2)
         bp%s1=sarc*EXP(ifac*pi/2)*dsarc
         bp%ds=(bp%s1/bp%s0)**(one/nyquist%ns)
         WRITE(bp%name,'(a,i2,1p,2(a,2e11.3))')
     $        "pbranc ",iarc,", s0 = ",bp%s0,", s1 = ",bp%s1
         sarc=sarc*dsarc
      ENDDO
c-----------------------------------------------------------------------
c     set up negative radial branches.
c-----------------------------------------------------------------------
      sarc=nyquist%small
      dsarc=(nyquist%big/nyquist%small)**(one/nyquist%narc)
      DO iarc=1,narc
         bp => nbranch(iarc)
         bp%s0=sarc*EXP(-ifac*pi/2)
         bp%s1=sarc*EXP(-ifac*pi/2)*dsarc
         bp%ds=(bp%s1/bp%s0)**(one/nyquist%ns)
         WRITE(bp%name,'(a,i2,1p,2(a,2e11.3))')
     $        "nbranch ",iarc,", s0 = ",bp%s0,", s1 = ",bp%s1
         sarc=sarc*dsarc
      ENDDO
c-----------------------------------------------------------------------
c     set up semicircular branches.
c-----------------------------------------------------------------------
      sarc=nyquist%small
      DO iarc=0,narc
         bp => cbranch(iarc)
         bp%s0=sarc*EXP(-ifac*pi/2)
         bp%s1=sarc*EXP(ifac*pi/2)
         bp%ds=(bp%s1/bp%s0)**(one/nyquist%ns)
         WRITE(bp%name,'(a,i2,1p,2(a,2e11.3))')
     $        "cbranch ",iarc,", s0 = ",bp%s0,", s1 = ",bp%s1
         sarc=sarc*dsarc
      ENDDO
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(nyquist%out)OPEN(UNIT=debug_unit,FILE="nyquist.out",
     $     STATUS="REPLACE")
      IF(nyquist%bin)OPEN(UNIT=bin_unit,FILE="nyquist.bin",
     $     STATUS="REPLACE",FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     scan over branches.
c-----------------------------------------------------------------------
      DO iarc=1,narc
         CALL match_branch(pbranch(iarc))
      ENDDO
      IF(nyquist%bin)WRITE(bin_unit)
      DO iarc=1,narc
         CALL match_branch(nbranch(iarc))
      ENDDO
      IF(nyquist%bin)WRITE(bin_unit)
      DO iarc=0,narc
         CALL match_branch(cbranch(iarc))
         IF(nyquist%bin)WRITE(bin_unit)
      ENDDO
c-----------------------------------------------------------------------
c     close output files.
c-----------------------------------------------------------------------
      IF(nyquist%out)CLOSE(debug_unit)
      IF(nyquist%bin)CLOSE(bin_unit)
c-----------------------------------------------------------------------
c     compute and print number of roots.
c-----------------------------------------------------------------------
      WRITE(out_unit,10)
      WRITE(out_unit,20)2,0,zero,ABS(cbranch(0)%s0)
      DO iarc=1,narc
         cycles=(SUM(nbranch(1:iarc)%cycles)
     $        +SUM(pbranch(narc-iarc+1:narc)%cycles)
     $        +cbranch(iarc)%cycles-cbranch(0)%cycles)/(twopi*ifac)
         roots=NINT(REAL(cycles))
         error=ABS(cycles-roots)
         WRITE(out_unit,20)iarc+2,roots,error,ABS(cbranch(iarc)%s0)
      ENDDO
      WRITE(out_unit,10)
      DEALLOCATE(pbranch,nbranch,cbranch)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination for nyquist plot.")
      END SUBROUTINE match_nyquist
c-----------------------------------------------------------------------
c     subprogram 10. match_branch.
c     draws one branch of nyquist plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_branch(branch)

      TYPE(branch_type), INTENT(INOUT) :: branch

      INTEGER :: is
      COMPLEX(r8) :: det,s
      COMPLEX(r8), DIMENSION(msing) :: qfac
      COMPLEX(r8), DIMENSION(2*msing,2*msing) :: mat

      COMPLEX(r8), SAVE :: det_old
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"is",5x,"re s",7x,"im s",5x,"re qfac",4x,"im qfac",
     $     5x,"re det",5x,"im det"/)
 20   FORMAT(i6,1p,6e11.3)
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      IF(nyquist%out)THEN
         WRITE(debug_unit,'(1x,a,":")')TRIM(branch%name)
         WRITE(debug_unit,10)
      ENDIF
c-----------------------------------------------------------------------
c     compute det and cycles.
c-----------------------------------------------------------------------
      s=branch%s0
      DO is=0,nyquist%ns
         det=match_delta(s,mat)
         IF(is == 0)THEN
            branch%cycles=0
            det_old=det
         ELSE
            branch%cycles=branch%cycles
     $           +2*(det-det_old)/(det+det_old)
            det_old=det
         ENDIF
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
         IF(nyquist%bin)WRITE(bin_unit)
     $        REAL(REAL(s),4),REAL(IMAG(s),4),
     $        REAL(qfac(1),4),REAL(IMAG(qfac(1)),4),
     $        REAL(REAL(det),4),REAL(IMAG(det),4)
         IF(nyquist%out)WRITE(debug_unit,20)is,s,qfac(1),det
         s=s*branch%ds
      ENDDO
c-----------------------------------------------------------------------
c     write trailer.
c-----------------------------------------------------------------------
      IF(nyquist%out)WRITE(debug_unit,10)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_branch      
c-----------------------------------------------------------------------
c     subprogram 11. match_solve.
c     solves for one root.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_solve(eigval)

      COMPLEX(r8), INTENT(OUT) :: eigval
            
      INTEGER :: iter
      REAL(r8), PARAMETER :: eps=1e-4
      REAL(r8) :: err
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      eigval=initguess
      iroot=0
      ALLOCATE(oldroots(nroot-1))
      oldroots=0
c-----------------------------------------------------------------------
c     solve for root and write output.
c-----------------------------------------------------------------------
      DO
         iroot=iroot+1
         CALL match_newton(match_delta,eigval,err,iter)
         WRITE(*,'(a,1p,e10.3,a,i4,a,2e11.3)')
     $        " eta = ",eta(1),", it = ",iter,", eigval = ",eigval
         IF(iroot == nroot)EXIT
         oldroots(iroot)=eigval
         eigval=eigval*(1+eps)
      ENDDO
      DEALLOCATE(oldroots)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_solve     
c-----------------------------------------------------------------------
c     subprogram 12. match_matrix_diagnose
c     solves for one root.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_matrix_diagnose(m,z)

      COMPLEX(r8), DIMENSION(:,:), INTENT(IN) :: m
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: z

      CHARACTER(256) format1,format2
      INTEGER :: n,i
      COMPLEX(r8), DIMENSION(SIZE(z)) :: matvec
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"i",',i2.2,'(4x,"re m",i2.2,5x,"im m",i2.2,1x),5x,'
     $     '"re z",7x,"im z",4x,"re matvec",2x,"im matvec"/)')
 20   FORMAT('(i4,1p,',i2.2,'e11.3)')
c-----------------------------------------------------------------------
c     create formats and compute d.
c-----------------------------------------------------------------------
      n=SIZE(z)
      WRITE(format1,10)n
      WRITE(format2,20)2*(n+2)
      matvec=MATMUL(m,z)
c-----------------------------------------------------------------------
c     write data to file.
c-----------------------------------------------------------------------
      OPEN(UNIT=debug_unit,FILE="matrix.out",STATUS="REPLACE")
      WRITE(debug_unit,format1)(i,i,i=1,n)
      WRITE(debug_unit,format2)(i,m(i,:),z(i),matvec(i),i=1,n)
      WRITE(debug_unit,format1)(i,i,i=1,n)
      CLOSE(UNIT=debug_unit)
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      CALL program_stop("match_matrix_diagnose: abort after diagnose.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_matrix_diagnose
c-----------------------------------------------------------------------
c     subprogram 13. match_rpec.
c     construct resistive perturbed equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_rpec
      
      CHARACTER(100),DIMENSION(2) :: filename
      INTEGER :: nmat,ising,idx1,idx2,idx3,idx4,info,ip,ipert
      INTEGER :: jsol,isol,countsing
      INTEGER, DIMENSION(4*msing-1) :: ipiv
      COMPLEX(r8) :: delta1,delta2
      COMPLEX(r8), DIMENSION(4*msing,coil%mcoil) :: cof,rmat
      COMPLEX(r8), DIMENSION(4*msing,4*msing) :: mat,cmat
      COMPLEX(r8), DIMENSION(2*msing,coil%mcoil) :: cout,cin
      COMPLEX(r8), DIMENSION(:,:,:),ALLOCATABLE :: globalsol
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
10    FORMAT("it",9x,"err",7x,"re z",7x,"im z",7x,"re f",7x,"im f",$)
c-----------------------------------------------------------------------
c     compute matching matrix for rpec.
c-----------------------------------------------------------------------
      mat=0
      rmat=0
      mat(2*msing+1:4*msing,1:2*msing)
     $     =TRANSPOSE(delta(1:2*msing,1:2*msing))
      rmat(2*msing+1:4*msing,1:coil%mcoil)
     $     =-TRANSPOSE(delta(coil%m1:coil%m2,1:2*msing))
      deltar=0
c-----------------------------------------------------------------------
c     start loop over singular surfaces.
c-----------------------------------------------------------------------
      DO ising=1,msing
         idx1=ising*2-1
         idx2=ising*2
         idx3=idx1+2*msing
         idx4=idx2+2*msing
c-----------------------------------------------------------------------
c     compute inner region matching data.
c-----------------------------------------------------------------------
         CALL deltac_run(restype(ising),initguess,deltar(ising,:),
     $                   deltaf(ising,:,:))
         delta1=deltar(ising,1)
         delta2=deltar(ising,2)
c-----------------------------------------------------------------------
c     construct the matching matrix.
c-----------------------------------------------------------------------
         mat(idx1,idx1)=1
         mat(idx2,idx2)=1
         mat(idx1,idx3)=-1
         mat(idx1,idx4)=1
         mat(idx2,idx3)=-1
         mat(idx2,idx4)=-1
         mat(idx3,idx3)=-delta1
         mat(idx3,idx4)=delta2
         mat(idx4,idx3)=-delta1
         mat(idx4,idx4)=-delta2
c-----------------------------------------------------------------------
c     finish loop over singular surfaces.
c-----------------------------------------------------------------------
      ENDDO
c-----------------------------------------------------------------------
c     compute the coefficients of outter and inner region solutions.
c-----------------------------------------------------------------------
      cmat=mat
      cof=rmat
      nmat=4*msing
      CALL zgetrf(nmat,nmat,cmat,nmat,ipiv,info)
      CALL zgetrs('N',nmat,coil%mcoil,cmat,nmat,ipiv,cof,nmat,info)
      cout=cof(1:2*msing,:)
      cin=cof(2*msing+1:4*msing,:)
c-----------------------------------------------------------------------
c     output inner and outer regions' solutions.
c-----------------------------------------------------------------------
      CALL match_alloc_sol(initguess)
      ALLOCATE (globalsol(outs%mpert,0:outs%tot_grids,coil%mcoil))
      jsol=0
      DO isol=coil%m1,coil%m2
         jsol=jsol+1
         WRITE (filename(2),"(I)") jsol  
         filename(2)=ADJUSTL(filename(2))
         WRITE(filename(1),*) 'rpec_sol_'//TRIM(filename(2))
         CALL match_output_solution(cout(:,jsol),cin(:,jsol),isol,
     $                              globalsol(:,:,jsol),filename(1))
      ENDDO
c-----------------------------------------------------------------------
c     output for coupling to the coil (final rpec run).
c-----------------------------------------------------------------------
      countsing=0
      DO ip=0,outs%tot_grids
         IF (outs%issing(ip)) countsing=countsing+1
      ENDDO
      CALL bin_open(bin_unit,"globalsol.bin","REPLACE","REWIND","none")
      WRITE (bin_unit) outs%mpert,outs%tot_grids-countsing,coil%mcoil,
     $                 outs%mlow,outs%mhigh 
      WRITE (bin_unit) outs%psi
      DO isol=1,coil%mcoil
         DO ip=0,outs%tot_grids
            IF (outs%issing(ip)) THEN
               CYCLE
            ENDIF 
            DO ipert=1,outs%mpert
               WRITE (bin_unit) globalsol(ipert,ip,isol)
            ENDDO
         ENDDO
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------     
      CALL match_dealloc_sol
      DEALLOCATE (globalsol)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_rpec
c-----------------------------------------------------------------------
c     subprogram 14. match_alloc_sol.
c     allocate and read inner and outer region solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE match_alloc_sol(eig)
      COMPLEX(r8), INTENT(IN) :: eig
      
      INTEGER isol,ising
      COMPLEX(r8), DIMENSION(2) :: deltai
      COMPLEX(r8), DIMENSION(2,2) :: df
c-----------------------------------------------------------------------
c     read outer region solutions.
c-----------------------------------------------------------------------
      CALL bin_open(bin_unit,galsol_filename,"OLD","REWIND","none")
      READ(bin_unit) outs%nsol,outs%msing,outs%mpert,outs%mhigh,
     $               outs%mlow,outs%nn,outs%psio,outs%tot_grids
      ALLOCATE(outs%psi(0:outs%tot_grids),
     $         outs%issing(0:outs%tot_grids),
     $         outs%q(0:outs%tot_grids))
      ALLOCATE(outs%qpsifac(outs%msing),outs%qsing(outs%msing))
      ALLOCATE(outs%sols(outs%mpert,0:outs%tot_grids,outs%nsol))
      ALLOCATE(outs%sols_cut(outs%mpert,0:outs%tot_grids,outs%nsol))
      ALLOCATE(outs%xext(outs%msing,2))
      DO isol=1,outs%msing
         READ(bin_unit) outs%qpsifac(isol),outs%qsing(isol)
      ENDDO
      READ(bin_unit) outs%psi,outs%issing,outs%q
      DO isol=1,outs%nsol
         READ (bin_unit) outs%sols(:,:,isol)
      ENDDO
      CALL bin_close(bin_unit)

      CALL bin_open(bin_unit,galsol_filename_cut,"OLD",
     $                 "REWIND","none")
      READ (bin_unit) outs%xext
      DO isol=1,outs%nsol
         READ (bin_unit) outs%sols_cut(:,:,isol)
      ENDDO
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     get inner region solutions
c-----------------------------------------------------------------------      
      ALLOCATE(ins%sols(msing))
      ins%tot_g=interp_np*nx
      output_sol=.TRUE.
      DO ising=1,msing
         sol => ins%sols(ising)
         ALLOCATE (sol%xvar(0:ins%tot_g),sol%sol(3,0:ins%tot_g,2))
         CALL deltac_run(restype(ising),eig,deltai,df)
      ENDDO
      output_sol=.FALSE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_alloc_sol
c-----------------------------------------------------------------------
c     subprogram 15. match_dealloc_sol.
c     deallocate inner and outer region solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------            
      SUBROUTINE match_dealloc_sol
      INTEGER ising
c-----------------------------------------------------------------------
c     deallocate outsol.
c-----------------------------------------------------------------------
      DEALLOCATE(outs%psi,outs%issing,outs%q)
      DEALLOCATE(outs%qpsifac,outs%qsing)
      DEALLOCATE(outs%sols)
      DEALLOCATE(outs%sols_cut)
      DEALLOCATE(outs%xext)
c-----------------------------------------------------------------------
c     deallocate insol.
c-----------------------------------------------------------------------      
      DO ising=1,msing
         DEALLOCATE(ins%sols(ising)%xvar,ins%sols(ising)%sol)
      ENDDO
      DEALLOCATE (ins%sols)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE match_dealloc_sol
c-----------------------------------------------------------------------
c     subprogram 16. match_output_solution.
c     output inner and outer region solutions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE match_output_solution(cout,cin,csol,globalsol,filename)
      CHARACTER(*),INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: csol
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: cout,cin
      COMPLEX(r8), DIMENSION(:,0:), INTENT(OUT) :: globalsol
      
      CHARACTER(100) :: filename1
      CHARACTER(100) :: comp_tittle,tmp
      INTEGER :: isol,ising,ip,ipert,eff_grids,jsing,comp,m
      REAL(r8) :: inpsi,chi1,sfac,q1,sig
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psi_cut,singfac
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: inpsifac
      COMPLEX(r8) :: insol,x0
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: outtotsol
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: outtotsol_cut,tmp_cut
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: intotsol
      TYPE(cspline_type) :: outcut_sp,q_sp
c-----------------------------------------------------------------------
c     allocate.
c-----------------------------------------------------------------------      
      ALLOCATE(outtotsol(outs%mpert,0:outs%tot_grids))
      ALLOCATE(singfac(outs%mpert))
      ALLOCATE(outtotsol_cut(outs%mpert,0:outs%tot_grids),
     $   tmp_cut(outs%mpert,0:outs%tot_grids),psi_cut(0:outs%tot_grids))
      ALLOCATE(intotsol(-ins%tot_g:ins%tot_g,msing),
     $   inpsifac(-ins%tot_g:ins%tot_g,msing))
c-----------------------------------------------------------------------
c     construct linear combination of outer region solutions with cofout.
c-----------------------------------------------------------------------
      outtotsol=0
      outtotsol_cut=0
      DO isol=1,2*msing
         outtotsol=outtotsol+cout(isol)*outs%sols(:,:,isol)
         outtotsol_cut=outtotsol_cut+cout(isol)*outs%sols_cut(:,:,isol)
      ENDDO
      IF (coil%ideal_flag) THEN
         outtotsol=outs%sols(:,:,csol)
      ELSE
         outtotsol=outtotsol+outs%sols(:,:,csol)
      ENDIF
      outtotsol_cut=outtotsol_cut+outs%sols_cut(:,:,csol)
c-----------------------------------------------------------------------
c     construct inner region solutions for each singular surface.
c-----------------------------------------------------------------------
      comp=2
      sig=1
      IF (match_sol%b_flag) THEN
         sig=-1
         comp=1
      ENDIF
      DO ising=1,msing
         sol => ins%sols(ising)          
         DO ip=-ins%tot_g,-1
            inpsifac(ip,ising)=outs%qpsifac(ising)-sol%xvar(-ip)
            intotsol(ip,ising)=-sol%sol(comp,-ip,1)*cin(2*ising)*sig
     $                         +sol%sol(comp,-ip,2)*cin(2*ising-1)*sig
         ENDDO
         DO ip=0,ins%tot_g
            inpsifac(ip,ising)=outs%qpsifac(ising)+sol%xvar(ip)
            intotsol(ip,ising)=sol%sol(comp,ip,1)*cin(2*ising)
     $                        +sol%sol(comp,ip,2)*cin(2*ising-1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     convert xi to b field.
c-----------------------------------------------------------------------      
      IF (match_sol%b_flag) THEN
      CALL cspline_alloc(q_sp,outs%tot_grids,1)
      q_sp%xs=outs%psi
      q_sp%fs(:,1)=outs%q
      CALL cspline_fit(q_sp,"extrap")
      chi1=twopi*outs%psio
      DO ip=0,outs%tot_grids
            singfac=outs%mlow-outs%nn*outs%q(ip)
     $             +(/(ipert,ipert=0,outs%mpert-1)/)
            outtotsol(:,ip)=chi1*ifac*singfac*outtotsol(:,ip)
            outtotsol_cut(:,ip)=chi1*ifac*singfac*outtotsol_cut(:,ip)
      ENDDO
      DO ising=1,msing
         CALL cspline_eval(q_sp,outs%qpsifac(ising),1)
         q1=q_sp%f1(1)
         sfac=restype(ising)%taur/restype(ising)%taua
         x0=sfac**(-1._r8/3._r8)
         intotsol(:,ising)=intotsol(:,ising)*chi1*ifac*outs%nn*q1*x0
      ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     final construction of perturbed equilibrium.
c-----------------------------------------------------------------------      
      IF (coil%ideal_flag) THEN
      ELSE
         IF (match_sol%auto_connect) THEN
            CALL match_auto_connect (csol,cout,inpsifac,
     $                               intotsol,outtotsol)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     write full outer region solutions, binary.
c-----------------------------------------------------------------------
      WRITE(filename1,*) TRIM(filename)//'_out.bin'
      CALL bin_open(bin_unit,filename1,"REPLACE","REWIND","none")
      DO ipert=1,outs%mpert
         DO ip=0,outs%tot_grids                
            IF (outs%issing(ip)) THEN
c               WRITE(bin_unit)
               CYCLE
            ENDIF
            WRITE (bin_unit) REAL(outs%psi(ip),4),
     $                       REAL(outtotsol(ipert,ip),4),
     $                       REAL(IMAG(outtotsol(ipert,ip)),4),
     $                       mylog(outtotsol(ipert,ip))
         ENDDO
         WRITE(bin_unit)
      ENDDO  
      CALL bin_close(bin_unit)
c-----------------------------------------------------------------------
c     write full outer region solutions, ascii.
c-----------------------------------------------------------------------
      WRITE(filename1,*) TRIM(filename)//'_out.out'
      CALL ascii_open(match_unit,TRIM(filename1),"REPLACE")
      WRITE (match_unit,10) 'psifac'
               DO m=outs%mlow,outs%mhigh
                  WRITE (tmp,"(I)") m
                  tmp=ADJUSTL(tmp)
                  WRITE (comp_tittle,*) 'REAL(',TRIM(tmp),')'
                  WRITE (match_unit,10) TRIM(comp_tittle)
                  WRITE (comp_tittle,*) 'IMAG(',TRIM(tmp),')'
                  WRITE (match_unit,10) TRIM(comp_tittle)
               ENDDO
10            FORMAT (1P,A15,$)
               WRITE (match_unit,*)
               DO ip=0,outs%tot_grids
                  IF (outs%issing(ip)) CYCLE
                  WRITE (match_unit,20) outs%psi(ip)
20               FORMAT (1P,E20.10,$)
                  DO ipert=1,outs%mpert
                     WRITE (match_unit,20)
     $                     REAL(outtotsol(ipert,ip)),
     $                     IMAG(outtotsol(ipert,ip))
                  ENDDO
                  WRITE (match_unit,*)
               ENDDO
               CALL ascii_close(match_unit)

c-----------------------------------------------------------------------
c     write inner region solutions, binary.
c-----------------------------------------------------------------------      
      tmp_cut=outtotsol_cut
      outtotsol_cut=0
      eff_grids=-1
      DO ip=0,outs%tot_grids
         IF (outs%issing(ip)) THEN
            CYCLE
         ENDIF
         eff_grids=eff_grids+1
         outtotsol_cut(:,eff_grids)=tmp_cut(:,ip)
         psi_cut(eff_grids)=outs%psi(ip)
      ENDDO
      CALL cspline_alloc(outcut_sp,eff_grids,outs%mpert)
      outcut_sp%xs(0:eff_grids)=psi_cut(0:eff_grids)
      DO ip=0,eff_grids
         outcut_sp%fs(ip,:)=outtotsol_cut(:,ip)
      ENDDO
      CALL cspline_fit(outcut_sp,"extrap")
      WRITE(filename1,*) TRIM(filename)//'_in.bin'
      CALL bin_open(bin_unit,filename1,"REPLACE","REWIND","none")
      DO ising=1,msing
         DO ip=-ins%tot_g,ins%tot_g
            IF (match_sol%uniform) THEN
               inpsi=inpsifac(ip,ising)
               DO jsing=1,outs%msing
                  IF (outs%xext(ising,1)<inpsi .AND. 
     $                inpsi<outs%xext(ising,2)) THEN 
                     ipert=NINT(outs%nn*outs%qsing(ising))-outs%mlow+1
                     CALL cspline_eval(outcut_sp,inpsi,0)
                     insol=outcut_sp%f(ipert)+intotsol(ip,ising)
                     WRITE (bin_unit) REAL(inpsifac(ip,ising),4),
     $                                REAL(insol,4),
     $                                REAL(IMAG(insol),4),
     $                                mylog(insol)
                  ENDIF
               ENDDO
            ELSE
               insol=intotsol(ip,ising)
               WRITE (bin_unit) REAL(inpsifac(ip,ising),4),
     $                          REAL(insol,4),
     $                          REAL(IMAG(insol),4),
     $                          mylog(insol)

            ENDIF
         ENDDO
         WRITE(bin_unit)
      ENDDO 
      CALL bin_close(bin_unit)  
      CALL cspline_dealloc(outcut_sp)
c-----------------------------------------------------------------------
c     write resonant outer region solutions, binary.
c-----------------------------------------------------------------------
      IF (match_sol%flag) THEN
         WRITE(filename1,*) TRIM(filename)//'_out_qpert.bin'
         CALL bin_open(bin_unit,filename1,"REPLACE","REWIND","none")
         DO ip=0,outs%tot_grids                
            IF (outs%issing(ip)) THEN
c               WRITE(bin_unit)
               CYCLE
            ENDIF
            ipert=NINT(outs%nn*match_sol%qpert)-outs%mlow+1
            WRITE (bin_unit) REAL(outs%psi(ip),4),
     $                       REAL(outtotsol(ipert,ip),4),
     $                       REAL(IMAG(outtotsol(ipert,ip)),4),
     $                       mylog(outtotsol(ipert,ip))
         ENDDO
         WRITE(bin_unit)
         CALL bin_close(bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     save to global solution for coulping to coils.
c-----------------------------------------------------------------------      
      globalsol=outtotsol
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(outtotsol)
      DEALLOCATE(outtotsol_cut,tmp_cut,psi_cut)
      DEALLOCATE(intotsol,inpsifac,singfac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE match_output_solution
c-----------------------------------------------------------------------
c     subprogram 17. match_auto_connect.
c     auto connect the outer and inner region.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      SUBROUTINE match_auto_connect(csol,cout,inpsifac,
     $                              intotsol,outtotsol)
      INTEGER, INTENT(IN) :: csol
      REAL(r8), DIMENSION(-ins%tot_g:ins%tot_g,msing),
     $   INTENT(IN) :: inpsifac
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: cout
      COMPLEX(r8), DIMENSION(-ins%tot_g:ins%tot_g,msing),
     $   INTENT(IN):: intotsol
      COMPLEX(r8), DIMENSION(outs%mpert,0:outs%tot_grids),
     $   INTENT(INOUT):: outtotsol
      
      INTEGER :: ipsi,ngrids,ipert,idx,ising,jsing
      INTEGER, DIMENSION(0:msing+1) :: idxsing
      INTEGER, DIMENSION(msing,2) :: idxconnect
      REAL(r8) :: psifac,singfac,di,dpsi
      REAL(r8), DIMENSION(2) :: psibou,alpha
      REAL(r8), DIMENSION(0:msing+1) :: psising
      COMPLEX(r8) :: outsols,insols,diffsols,ul0,us0,u0
      TYPE(cspline_type), DIMENSION(msing) :: insp
c-----------------------------------------------------------------------
c     find auto connect interval.
c-----------------------------------------------------------------------      
      psising(0)=outs%psi(0)
      psising(1:msing)=outs%qpsifac      
      psising(msing+1)=outs%psi(outs%tot_grids)
      idxsing(0)=0
      ising=1
      DO ipsi=0,outs%tot_grids    
         IF (outs%issing(ipsi)) THEN
            idxsing(ising)=ipsi
            ising=ising+1
         ENDIF
      ENDDO
      idxsing(ising)=outs%tot_grids
c-----------------------------------------------------------------------
c     scan psi and find the point to connect the outer and inner region.
c-----------------------------------------------------------------------
      ngrids=2*ins%tot_g
      DO ising=1,msing
         CALL cspline_alloc(insp(ising),ngrids,1)
         insp(ising)%xs(0:ngrids)=inpsifac(-ins%tot_g:ins%tot_g,ising)
         insp(ising)%fs(:,1)=intotsol(:,ising)
         CALL cspline_fit(insp(ising),"extrap")
         psibou(1)=psising(ising-1)
         psibou(2)=psising(ising+1)
         IF (psibou(1) < inpsifac(-ins%tot_g,ising)) THEN
            psibou(1)=inpsifac(-ins%tot_g,ising)
         ENDIF
         IF (psibou(2) > inpsifac(ins%tot_g,ising)) THEN
            psibou(2)=inpsifac(ins%tot_g,ising)
         ENDIF
c-----------------------------------------------------------------------
c     compare the outer/inner region and generate global solutions.
c-----------------------------------------------------------------------
         singfac=outs%q(idxsing(ising))
         ipert=NINT(outs%nn*singfac)-outs%mlow+1
         di=restype(ising)%e+restype(ising)%f
     $           +restype(ising)%h-0.25
         alpha(1)=-0.5+SQRT(-di)
         alpha(2)=-0.5-SQRT(-di)   
         DO ipsi=idxsing(ising)-1,idxsing(ising-1),-1
            psifac=outs%psi(ipsi)
            dpsi=psifac-outs%psi(idxsing(ising))
            us0=ABS(dpsi)**alpha(1)
            ul0=ABS(dpsi)**alpha(2)
            idx=2*ising
            IF (dpsi<0) idx=2*ising-1
            u0=cout(idx)*ul0
            DO jsing=1,2*msing
               u0=u0+cout(jsing)*delta(jsing,idx)*us0
            ENDDO
            u0=u0+delta(csol,idx)*us0
            CALL cspline_eval(insp(ising),psifac,0)
            outsols=outtotsol(ipert,ipsi)
            insols=insp(ising)%f(1)+outsols-u0
            diffsols=(outsols-insols)/insols
            IF (psifac<psibou(1)) THEN
               WRITE(*,*) "ising=",ising,"csol=",csol
               WRITE(*,*) "threshold or inpsifac maybe too small."
               EXIT
c               CALL program_stop("psifac<psibou(1).")
            ENDIF
            IF (ABS(diffsols)<match_sol%connect_threshold) THEN
               idxconnect(ising,1)=ipsi
               EXIT
            ENDIF
            outtotsol(ipert,ipsi)=insols
         ENDDO
         DO ipsi=idxsing(ising)+1,idxsing(ising+1)
            psifac=outs%psi(ipsi)
            dpsi=psifac-outs%psi(idxsing(ising))
            us0=ABS(dpsi)**alpha(1)
            ul0=ABS(dpsi)**alpha(2)
            idx=2*ising
            IF (dpsi<0) idx=2*ising-1
            u0=cout(idx)*ul0
            DO jsing=1,2*msing
               u0=u0+cout(jsing)*delta(jsing,idx)*us0
            ENDDO
            u0=u0+delta(csol,idx)*us0
            CALL cspline_eval(insp(ising),psifac,0)
            outsols=outtotsol(ipert,ipsi)
            insols=insp(ising)%f(1)+outsols-u0
            diffsols=(outsols-insols)/insols            
            IF (psifac>psibou(2)) THEN
               WRITE(*,*) "ising=",ising,"csol=",csol
               WRITE(*,*) "threshold or inpsifac  maybe too small."
               EXIT
c               CALL program_stop("psifac>psibou(2).")
            ENDIF
            IF (ABS(diffsols)<match_sol%connect_threshold) THEN
               idxconnect(ising,2)=ipsi
               EXIT
            ENDIF
            outtotsol(ipert,ipsi)=insols
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------      
      DO ising=1,msing
         CALL cspline_dealloc(insp(ising))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN      
      END SUBROUTINE match_auto_connect
      
      END MODULE match_mod
c-----------------------------------------------------------------------
c     subprogram 18. match_main.
c     trivial main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------      
      PROGRAM match_main
      USE match_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     do it.
c-----------------------------------------------------------------------  
      CALL match_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------               
      END PROGRAM match_main
