c-----------------------------------------------------------------------
c     file coil.f.
c     set up coil information.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. coil_mod.
c     1. coil_read.
c-----------------------------------------------------------------------
c     subprogram 0. coil_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE coil_mod
      USE bicube_mod
      USE fspline_mod

      IMPLICIT NONE

      INTEGER :: coil_unit=31,cmlow,cmhigh
      INTEGER :: cnn,cmpsi,cmtheta,cmzeta,cmpert,coil_num
      CHARACTER(24) :: machine,ceq_type,ip_direction,bt_direction
      CHARACTER(24), DIMENSION(10) :: coil_name
      REAL(r8), DIMENSION(10,48) :: coil_cur
      REAL(r8) :: cro,czo,cpsio,cpsilow,cpsilim,cqlim,
     $     ipd,btd,helicity

      LOGICAL :: gpec_interface
      CHARACTER(256) :: data_dir = 'default'
      CHARACTER(512) :: cfile
      INTEGER, DIMENSION(:), POINTER :: cmfac

      TYPE :: coil_type
      CHARACTER(24) :: coil_name
      INTEGER :: ncoil,s,nsec
      REAL(r8) :: nw
      REAL(r8), DIMENSION(:), POINTER :: cur 
      REAL(r8), DIMENSION(:,:,:), POINTER :: x,y,z
      END TYPE coil_type      

      TYPE(coil_type), DIMENSION(:), POINTER :: coil
      TYPE(spline_type) :: csq
      TYPE(bicube_type) :: crzphi

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. coil_read.
c     read coils.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE coil_read(cdconfile,icoil_num,icoil_name,icoil_cur)
      
      CHARACTER(128), INTENT(IN) :: cdconfile
      INTEGER, OPTIONAL :: icoil_num
      REAL(r8), DIMENSION(10,48), OPTIONAL :: icoil_cur
      CHARACTER(24), DIMENSION(10), OPTIONAL :: icoil_name

      NAMELIST/coil_control/ceq_type,cmpsi,cmtheta,cmzeta,cmlow,cmhigh,
     $     data_dir,machine,ip_direction,bt_direction,
     $     coil_num,coil_name,coil_cur
      NAMELIST/coil_output/gpec_interface

      INTEGER :: ci,cj,ck,cl,cm,ci1,ci2,ci3,ci4
      REAL(r8) :: cr1,cr2
c-----------------------------------------------------------------------
c     initialize and read input data.
c-----------------------------------------------------------------------
      cmpsi=64
      cmtheta=480
      cmzeta=32
      cmlow=-64
      cmhigh=64
      coil_cur=0
      CALL ascii_open(in_unit,"coil.in","OLD")
      READ(UNIT=in_unit,NML=coil_control)
      READ(UNIT=in_unit,NML=coil_output)
      CALL ascii_close(in_unit)
      IF (TRIM(data_dir)=='' .OR. TRIM(data_dir)=='default') THEN
         CALL getenv('GPECHOME',data_dir)
         IF(LEN(TRIM(data_dir))==0) stop
     $  "ERROR: Default coil dir requires GPECHOME environment variable"
         data_dir = TRIM(data_dir)//'/coil'
      ENDIF
      IF (present(icoil_num)) coil_num=icoil_num
      IF (present(icoil_name)) coil_name=icoil_name
      IF (present(icoil_cur)) coil_cur=icoil_cur
c-----------------------------------------------------------------------
c     read coils for each machine.
c-----------------------------------------------------------------------
      IF(coil_num>0)ALLOCATE(coil(coil_num))

      DO ci=1,coil_num
         cfile=TRIM(data_dir)//"/"//TRIM(machine)//"_"//
     $         TRIM(coil_name(ci))//".dat"
         CALL ascii_open(coil_unit,cfile,"old")
         coil(ci)%coil_name=TRIM(coil_name(ci))
         READ(coil_unit,'(3(1x,I4),1x,f7.2)')
     $        coil(ci)%ncoil,coil(ci)%s,coil(ci)%nsec,coil(ci)%nw
         ALLOCATE(coil(ci)%cur(coil(ci)%ncoil))
         ALLOCATE(coil(ci)%x(coil(ci)%ncoil,coil(ci)%s,coil(ci)%nsec),
     $        coil(ci)%y(coil(ci)%ncoil,coil(ci)%s,coil(ci)%nsec),
     $        coil(ci)%z(coil(ci)%ncoil,coil(ci)%s,coil(ci)%nsec))
         DO cj=1,coil(ci)%ncoil
            coil(ci)%cur(cj)=coil_cur(ci,cj)
            DO ck=1,coil(ci)%s
               DO cl=1,coil(ci)%nsec
                  READ(coil_unit,'(3(1x,e12.4))')coil(ci)%x(cj,ck,cl),
     $                 coil(ci)%y(cj,ck,cl),coil(ci)%z(cj,ck,cl)
               ENDDO
            ENDDO
         ENDDO
         CALL ascii_close(coil_unit)
      ENDDO
c-----------------------------------------------------------------------
c     read equilibrium information.
c-----------------------------------------------------------------------
      CALL bin_open(in_unit,cdconfile,"OLD","REWIND","none")
      READ(in_unit)ci1,ci2,cnn,ci3,ci4,cro,czo
      READ(in_unit)ci1,cr1,ci2,cpsio,cpsilow,cpsilim,cqlim,cr2
      READ(in_unit)
      READ(in_unit)
      READ(in_unit)

      CALL spline_alloc(csq,ci3,4)
      CALL bicube_alloc(crzphi,ci3,ci4,4)
      crzphi%periodic(2)=.TRUE.
      READ(in_unit)csq%xs,csq%fs,csq%fs1,csq%xpower
      READ(in_unit)crzphi%xs,crzphi%ys,
     $     crzphi%fs,crzphi%fsx,crzphi%fsy,crzphi%fsxy,
     $     crzphi%x0,crzphi%y0,crzphi%xpower,crzphi%ypower
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      cmpsi=MAXVAL((/64,cmpsi/))
      cmtheta=MAXVAL((/400,cmtheta/))
      cmzeta=cnn*MAXVAL((/40,cmzeta/))
      cmpert=cmhigh-cmlow+1
      ALLOCATE(cmfac(cmpert))
      cmfac=(/(cm,cm=cmlow,cmhigh)/)
      ipd=1.0
      btd=1.0
      IF(ip_direction=="negative")ipd=-1.0
      IF(bt_direction=="negative")btd=-1.0
      helicity=ipd*btd
      gpec_interface=.TRUE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE coil_read
c-----------------------------------------------------------------------
c     subprogram 2. coil_dealloc.
c     read coils.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE coil_dealloc

      INTEGER :: i

      DO i=1,coil_num
         DEALLOCATE(coil(i)%cur,coil(i)%x,coil(i)%y,coil(i)%z)
      ENDDO
      DEALLOCATE(coil)
      DEALLOCATE(cmfac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE coil_dealloc

      END MODULE coil_mod
