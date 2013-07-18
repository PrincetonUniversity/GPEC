c-----------------------------------------------------------------------
c     PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT
c     input/output opperations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. pentio_mod
c     1. pentio_read_order1_funs
c     2. pentio_read_order1_mns
c     3. pentio_read_profile
c     4. pentio_read_fperz
c     5. pentio_specialfunctions
c     6. pentio_write_bxpflag
c     7. pentio_write_profile
c     8. pentio_xheader
c     9. pentio_lmdaheader
c    10. pentio_stop
c-----------------------------------------------------------------------
c     subprogram 0. pentio_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE pentio_mod
      USE pentglobal_mod
      USE bicube_mod

      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     Subroutine 1. pentio_read_order1_funs
c     read pmodb and divxprp functions (psi,theta) from input file.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_read_order1_funs(lagbpar,divxprp,infile,native)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      LOGICAL, INTENT(IN) :: native
      CHARACTER(*), INTENT(IN) :: infile
      COMPLEX(r8), DIMENSION(mstep,0:mthsurf), INTENT(INOUT) :: 
     $     lagbpar,divxprp

      INTEGER :: mt,ms,istep,itheta

      TYPE(bicube_type) o1spl
c-----------------------------------------------------------------------
c     read file
c-----------------------------------------------------------------------
      WRITE(*,*) "Reading pmodb and divxprp input file: "//TRIM(infile)

      IF(native)THEN
         CALL bin_open(in_unit,infile,"OLD","REWIND","none")
         READ(in_unit)ms,mt
         IF((ms .NE. mstep) .OR. (mt .NE. mthsurf)) CALL
     $        pentio_stop("Input perturbations don't match equilibrium")
         READ(in_unit)lagbpar
         READ(in_unit)divxprp
         CALL bin_close(in_unit)
      ELSE
         CALL ascii_open(in_unit,infile,"OLD")
         READ(in_unit,*) ms,mt
         CALL bicube_alloc(o1spl,ms-1,mt-1,2)
         DO istep = 0,ms-1
            DO itheta = 0,mt-1
               READ(in_unit,*) o1spl%xs(istep),o1spl%ys(itheta),
     $              o1spl%fs(istep,itheta,:)
            ENDDO
         ENDDO
         CALL ascii_close(in_unit)
         CALL bicube_fit(o1spl,"extrap","periodic")
         
         WRITE(*,*) "Extrapolating to DCON/IPEC/PENT (psi,theta)"
         DO istep=1,mstep
            DO itheta=0,mthsurf
               CALL bicube_eval(o1spl,psifac(istep),theta(itheta),0)
               lagbpar(istep,itheta) = o1spl%f(1)
               divxprp(istep,itheta) = o1spl%f(2)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_read_order1_funs

c-----------------------------------------------------------------------
c     Subroutine 2. pentio_read_order1_mns
c     read pmodb and divxprp fourier components (psi,m) from input file.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_read_order1_mns(lagbpar,divxprp,infile,native)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      LOGICAL, INTENT(IN) :: native
      CHARACTER(*), INTENT(IN) :: infile
      COMPLEX(r8), DIMENSION(mstep,mpert), INTENT(INOUT) :: 
     $     lagbpar,divxprp

      INTEGER :: ms,mp,istep,ipert,jpert
      INTEGER, DIMENSION(:), ALLOCATABLE :: mtemp
      TYPE(spline_type) :: o1bspl,o1xspl
c-----------------------------------------------------------------------
c     read file
c-----------------------------------------------------------------------
      WRITE(*,*) "Reading pmodb and divxprp input file: "//TRIM(infile)

      IF(native)THEN
         CALL bin_open(bin_unit,infile,"OLD","REWIND","none")
         READ(bin_unit) ms,mp
         IF((ms .NE. mstep) .OR. (mp .NE. mpert)) CALL pentio_stop(
     $        "Input perturbations don't match equilibrium")
         READ(bin_unit)lagbpar
         READ(bin_unit)divxprp
         CALL bin_close(bin_unit)
      ELSE
         CALL ascii_open(in_unit,infile,"OLD")
         READ(in_unit,*) ms,mp
         CALL spline_alloc(o1bspl,ms-1,mp)
         CALL spline_alloc(o1xspl,ms-1,mp)
         ALLOCATE(mtemp(mp))
         DO istep = 0,ms-1
            DO ipert = 1,mp
               READ(in_unit,*) o1bspl%xs(istep),mtemp(ipert),
     $              o1bspl%fs(istep,ipert),o1xspl%fs(istep,ipert)
            ENDDO
         ENDDO
         CALL ascii_close(in_unit)
         CALL spline_fit(o1bspl,"extrap")
         CALL spline_fit(o1xspl,"extrap")
         
         WRITE(*,*) "Extrapolating to DCON/IPEC/PENT (psi,m)"
         DO istep=1,mstep
            CALL spline_eval(o1bspl,psifac(istep),0)
            CALL spline_eval(o1xspl,psifac(istep),0)
            DO ipert=1,mpert
               lagbpar(istep,ipert) = 0
               divxprp(istep,ipert) = 0
               DO jpert = 1,mp
                  IF(mtemp(jpert)==mfac(ipert))THEN
                     lagbpar(istep,ipert) = o1bspl%f(jpert)
                     divxprp(istep,ipert) = o1xspl%f(jpert)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_read_order1_mns


      


c-----------------------------------------------------------------------
c     subprogram 3. pentio_read_profile.
c     read kinetic input file.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_read_profile(kfile)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(*), INTENT(IN) :: kfile

      INTEGER :: nlines,istep
      REAL(r8) :: psi,epsa
      REAL(r8), DIMENSION(:), ALLOCATABLE :: zeff,zpitch

c-----------------------------------------------------------------------
c     read kinetic inputs
c-----------------------------------------------------------------------
      WRITE(*,*) "Reading kinetic input file: "
      WRITE(*,*) TRIM(kfile)

      nlines=0
      CALL ascii_open(in_unit,kfile,"OLD")
      READ(in_unit,*) 
      DO
         READ(in_unit,*,END=999)
         nlines=nlines+1
      ENDDO
 999  REWIND(in_unit)
      
      CALL spline_alloc(kin,nlines-1,9)
      kin%title = (/"idens ","edens ","itemp ","etemp ","wexb  ",
     $     "loglam","  nuii","  nuei","  epsr"/)
      READ(in_unit,*)
      DO istep=0,nlines-1
         READ(in_unit,*) kin%xs(istep),kin%fs(istep,1:5)
      ENDDO
      CALL ascii_close(in_unit)
     
      IF(kin%xs(0)==0.0) THEN
         kin%xs(0) = kin%xs(1)/10.0
         WRITE(*,'(a40,F6.4)') "Warning: Forced kinetic input psi 0 to "
     $        ,kin%xs(0)
      ENDIF
      WRITE(*,'(a26,F6.4,a4,F6.4)') "Kinetic profiles from psi ",
     $     kin%xs(0)," to ",kin%xs(kin%mx)

      ALLOCATE(zeff(0:nlines-1),zpitch(0:nlines-1))
      zeff = zimp-(kin%fs(:,1)/kin%fs(:,2))*icharge
     $     *(zimp-icharge)
      zpitch = 1.0+(1.0+zmass)/(2.0*zmass)*zimp
     $     *(zeff-1.0)/(zimp-zeff)
      epsa = SQRT(SUM(rzphi%fs(rzphi%mx,:,1))/rzphi%my)/ro

c-----------------------------------------------------------------------
c     in SI units.
c-----------------------------------------------------------------------
      kin%fs(:,3:4) = kin%fs(:,3:4)*1.602E-19
      kin%fs(:,6) = 17.3-0.5*LOG(kin%fs(:,2)/1.0e20)
     $     +1.5*LOG(kin%fs(:,4)/1.602e-16)
      kin%fs(:,7) = (zpitch/3.5E17)*kin%fs(:,1)*kin%fs(:,6)
     $     /(SQRT(1.0*imass)*(kin%fs(:,3)/1.602e-16)**1.5)
      kin%fs(:,8) = (zpitch/3.5E17)*kin%fs(:,2)*kin%fs(:,6)
     $     /(SQRT(emass/pmass)*(kin%fs(:,4)/1.602e-16)**1.5)
      kin%fs(:,9) = epsa*SQRT(kin%xs(:))

      CALL spline_fit(kin,"extrap")
      
c-----------------------------------------------------------------------
c     Outputs if desired.
c-----------------------------------------------------------------------
      IF(kinetic_flag) THEN
         CALL ascii_open(out_unit,"pent_kinetics.out","UNKNOWN")
         WRITE(out_unit,*)"PROFILES USED IN PENT CALCULATION:"
         WRITE(out_unit,'(1/,a8,I4,2(a10,es16.8E3))') " npsi = ",1000,
     $        " charge = ",-echarge*icharge,"   mass = ",imass*pmass
         WRITE(out_unit,'(1/,1x,15(1x,a16))') "psi","eps_r","n_i","n_e",
     $        "T_i","T_e","omega_EXB","logLambda","nu_ii","nu_ei","q",
     $        "dconPmu_0","dvdpsi","omega_*n","omega_*t"
         DO istep=1,1000
            psi = istep/1000.0
            CALL spline_eval(sq,psi,0)
            CALL spline_eval(kin,psi,1)
            CALL bicube_eval(rzphi,psi,0.0_r8,0)
            WRITE(out_unit,'(1x,15(1x,es16.8E3))') psi,
     $           kin%f(9),kin%f(1:8),sq%f(4),sq%f(2),sq%f(3),-twopi*
     $           kin%f(3)*kin%f1(1)/(-echarge*icharge*chi1*kin%f(1)),
     $           -twopi*kin%f1(3)/(-echarge*icharge*chi1)
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_read_profile

c-----------------------------------------------------------------------
c     subprogram 4. pentio_feprz.
c     read special function input data.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_read_feprz(f,file)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      CHARACTER(*), INTENT(IN) :: file
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: f
      REAL(r8), DIMENSION(:), ALLOCATABLE :: r,z,e,mu,pitch
      INTEGER :: ne,np,nr

      CALL ascii_open(in_unit,file,"OLD")
      READ(in_unit,*) ne
      READ(in_unit,*) np
      READ(in_unit,*) nr

      CALL ascii_close(in_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_read_feprz

c-----------------------------------------------------------------------
c     subprogram 5. pentio_specialfunctions.
c     read special function input data.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_read_special_functions
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER :: nek,nfk,nft

c-----------------------------------------------------------------------
c     Complete elliptic Integral of the 1st kind from 0 to 1.
c-----------------------------------------------------------------------
      WRITE(*,*) "Reading elliptic integral from file: ellipk01.dat"
      CALL ascii_open(in_unit,"ellipk01.dat","OLD")
      READ(in_unit,*) nek
      CALL spline_alloc(ellipk,nek-1,1)
      READ(in_unit,*) ellipk%xs(:)
      READ(in_unit,*) ellipk%fs(:,1)
      CALL ascii_close(in_unit)
      CALL spline_fit(ellipk,'extrap')

c-----------------------------------------------------------------------
c     Integral Bounce Average function.
c-----------------------------------------------------------------------
      WRITE(*,*) "Reading F^-1/2_mnl from file: fkmnql.dat"
      CALL bin_open(in_unit,"fkmnl.dat","OLD","REWIND","none")
      READ(in_unit) nfk,nft
      CALL bicube_alloc(fmnl,nfk,nft,1)
      READ(in_unit)fmnl%xs(:)
      READ(in_unit)fmnl%ys(:)
      READ(in_unit)fmnl%fs(:,:,1)
      CALL bin_close(in_unit)
      CALL bicube_fit(fmnl,'extrap','extrap')

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_read_special_functions

c-----------------------------------------------------------------------
c     subprogram 6. pentio_write_bxpflag.
c     Write an output file for optional energy sampling.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_write_bxpflag(name,title,lbls,dims,dat,sumend)
      
      CHARACTER(*), INTENT(IN):: name,title,lbls
      INTEGER, OPTIONAL, INTENT(IN) :: sumend
      INTEGER, DIMENSION(:), INTENT(IN) :: dims
      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(IN) :: dat
      
      INTEGER :: l,i,j,k,nl,sumstart,sumstop
      REAL(r8) :: bo
      CHARACTER(20):: frmt
      WRITE(frmt,'(a4,I2,a14)') '(1x,',dims(5),'(1x,es16.8E3))'
      nl = dims(SIZE(dims)-1)
      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
c-----------------------------------------------------------------------
c     Writing.
c-----------------------------------------------------------------------
      WRITE(*,*) "Writing file "//TRIM(name)//".out"
      CALL ascii_open(out_unit,TRIM(name)//".out","UNKNOWN")
      WRITE(out_unit,*)"NON-AMBIPOLAR TRANSPORT and TOROIDAL TORQUE: "//
     $     TRIM(title)
      WRITE(out_unit,*)
      WRITE(out_unit,*)"R0 = ",ro," B0 = ",bo," chi1 = ",chi1," n = ",nn
      !Write sum over l
      IF(dims(4)>1)THEN
         sumstart = 4
         sumstop = dims(5)
         DO i=1,3
            IF(dims(i)==1) sumstart=sumstart-1
         ENDDO
         IF(PRESENT(sumend)) sumstop = sumend
         WRITE(out_unit,'(1/,a8,1/)') " l = sum"
         WRITE(out_unit,*) TRIM(lbls)
         DO i = 1,dims(3)
            DO j = 1,dims(2)
               DO k = 1,dims(1)
                  IF(dat(k,j,i,1,1)==0) EXIT
                  IF(sumstop<dims(5)) THEN
                     WRITE(out_unit,frmt) 
     $                    dat(k,j,i,1,1:sumstart-1),
     $                    SUM(dat(k,j,i,:,sumstart:sumstop),DIM=1),
     $                    dat(k,j,i,1,sumstop+1:dims(5))
                  ELSE
                     WRITE(out_unit,frmt) 
     $                    dat(k,j,i,1,1:sumstart-1),
     $                    SUM(dat(k,j,i,:,sumstart:sumstop),DIM=1) 
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !Write individual l contributions
      DO l = 1,dims(4)
         WRITE(out_unit,'(1/,a5,i4,1/)') " l = ",l-dims(4)/2-1
         WRITE(out_unit,*) TRIM(lbls)
         DO i = 1,dims(3)
            DO j = 1,dims(2)
               DO k = 1,dims(1)
                  IF(dat(k,j,i,l,1)==0) EXIT
                  WRITE(out_unit,frmt) dat(k,j,i,l,:)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_write_bxpflag

      
c-----------------------------------------------------------------------
c     subprogram 7. pentio_write_profile.
c     Write an output file for non ambipolar transport toroidal torque.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_write_profile(name,title,spl,wspl)
      
      CHARACTER(*), INTENT(IN) :: name,title
      TYPE(cspline_type) :: spl
      TYPE(cspline_type), OPTIONAL, INTENT(IN) :: wspl

      INTEGER :: ll,istep
      REAL(r8) :: bo

      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
c-----------------------------------------------------------------------
c     Writing.
c-----------------------------------------------------------------------
      WRITE(*,*) "Writing file "//TRIM(name)//".out"
      CALL ascii_open(out_unit,TRIM(name)//".out","UNKNOWN")
      WRITE(out_unit,*)"NON-AMBIPOLAR TRANSPORT and TOROIDAL TORQUE: "//
     $     TRIM(title)
      WRITE(out_unit,*)
      WRITE(out_unit,*)"R0 = ",ro," B0 = ",bo," chi1 = ",chi1," n = ",nn
      WRITE(out_unit,'(2/,a9,2(1x,a20,1x,es16.8E3),1/)') "l = sum  ",
     $     " total(T_phi) =",REAL(SUM(spl%fsi(spl%mx,:))),
     $     " total(2ndeltaW) =",AIMAG(SUM(spl%fsi(spl%mx,:)))
      WRITE(out_unit,'(1x,6(1x,a16))') "psi","T_phi","2ndeltaW",
     $        "int(T_phi)","int(2ndeltaW)","dvdpsi"
      DO istep=0,spl%mx
         CALL spline_eval(sq,spl%xs(istep),0)
         WRITE(out_unit,'(1x,6(1x,es16.8E3))') spl%xs(istep),
     $        SUM(spl%fs(istep,:)),SUM(spl%fsi(istep,:)),
     $        sq%f(3)
      ENDDO
      DO ll=1,spl%nqty
         WRITE(out_unit,'(1/,a7,i4,2(1x,a20,1x,es16.8E3),1/)') ' l = ',
     $        ll-spl%nqty/2-1,' total(T_phi) = ',REAL(spl%fsi(spl%mx,ll)
     $        ),' total(2ndeltaW) = ',AIMAG(spl%fsi(spl%mx,ll))
         IF(PRESENT(wspl))THEN
            WRITE(out_unit,'(1x,8(1x,a16))') "psi","T_phi","2ndeltaW",
     $           "int(T_phi)","int(2ndeltaW)","omega_rot",
     $           "omega_off_T","omega_off_W"
            DO istep=0,spl%mx
               WRITE(out_unit,'(1x,8(1x,es16.8E3))') spl%xs(istep),
     $              spl%fs(istep,ll),spl%fsi(istep,ll),
     $              REAL(wspl%fs(istep,1)),wspl%fs(istep,1+ll)
            ENDDO
         ELSE
            WRITE(out_unit,'(1x,5(1x,a16))') "psi","T_phi","2ndeltaW",
     $           "int(T_phi)","int(2ndeltaW)"
            DO istep=0,spl%mx
               WRITE(out_unit,'(1x,5(1x,es16.8E3))') spl%xs(istep),
     $              spl%fs(istep,ll),spl%fsi(istep,ll)
            ENDDO
         ENDIF
      ENDDO
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pentio_write_profile

c-----------------------------------------------------------------------
c     subprogram 8. pentio_xheader.
c     open/clear energy file and write header info.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_xheader()
      REAL(r8) :: bo
      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
c-----------------------------------------------------------------------
c     Get on-axis field for reference.
c-----------------------------------------------------------------------
      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
c-----------------------------------------------------------------------
c     Writing.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"pent_energy_n"//TRIM(sn)//".out",
     $     "UNKNOWN")
      WRITE(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT:"
      WRITE(out_unit,*) "Energy integrand"//newl//
     $     "Normalizations are as follows:"//newl//
     $     "Lambda: m*v_perp^2/(2B)"//newl//"x: E/T"//newl//
     $     "T_phi: T_phi/(-2n^2TN*chi'/sqrt(pi))"//newl//
     $     "2ndeltaW: 2*n*deltaW /(-2n^2TN*chi'/sqrt(pi))"
      WRITE(out_unit,*)
      WRITE(out_unit,*)"R0 = ",ro," B0 = ",bo," chi1 = ",chi1," n = ",nn
      WRITE(out_unit,'(7(1x,a16))') "psi","Lambda","x","T_phi",
     $        "2ndeltaW","int(T_phi)","int(2ndeltaW)"
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE pentio_xheader
      
c-----------------------------------------------------------------------
c     subprogram 8. pentio_lmdaheader.
c     open/clear pitch file and write header info.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_lmdaheader()
      REAL(r8) :: bo
      CHARACTER(2), PARAMETER :: newl = CHAR(13)//CHAR(10)
c-----------------------------------------------------------------------
c     Get on-axis field for reference.
c-----------------------------------------------------------------------
      CALL bicube_eval(eqfun,0.0_r8,0.0_r8,0)
      bo = eqfun%f(1)
c-----------------------------------------------------------------------
c     Writing.
c-----------------------------------------------------------------------
      CALL ascii_open(out_unit,"pent_pitch_n"//TRIM(sn)//".out",
     $     "UNKNOWN")
      WRITE(out_unit,*) "PERTURBED EQUILIBRIUM NONAMBIPOLAR TRANSPORT:"
      WRITE(out_unit,*) "Pitch angle integrand and functions"
     $        //newl//"Normalizations are as follows:"//newl//
     $     "omega_D: omega_D/(E/Ze)"//newl//
     $     "T_phi,2ndeltaW: T_phi,2ndeltaW/(-2n^2TN*chi'/sqrt(pi))"
      WRITE(out_unit,*)
      WRITE(out_unit,*)"R0 = ",ro," B0 = ",bo," chi1 = ",chi1," n = ",nn
      WRITE(outlbls,'(9(1x,a16))') "psi","Lambda","T_phi","2ndeltaW",
     $     "int(T_phi)","int(2ndeltaW)","omega_b","omega_D","|deltaJ|^2"
      CALL ascii_close(out_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE pentio_lmdaheader

c-----------------------------------------------------------------------
c     subprogram 10. pentio_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pentio_stop(message)

      CHARACTER(*), INTENT(IN) :: message
      CHARACTER(10) :: date,time,zone
      REAL(r8) :: seconds     
c-----------------------------------------------------------------------
c     write  message and stop program.
c-----------------------------------------------------------------------
      CALL program_stop("PENT: "//message)
      END SUBROUTINE pentio_stop
      


      END MODULE pentio_mod
