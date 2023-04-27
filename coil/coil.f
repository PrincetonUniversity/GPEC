c-----------------------------------------------------------------------
c     file coil.f.
c     set up coil information.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. coil_mod.
c     1. coil_read.
c     2. erchk.
c     3. write_coil_geometry.
c     4. coil_dealloc.
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
      USE netcdf

      IMPLICIT NONE

      INTEGER :: coil_unit=31,cmlow,cmhigh
      INTEGER :: cnn,cmpsi,cmtheta,cmzeta,cmpert,coil_num
      INTEGER :: cnpert=-1
      CHARACTER(24) :: machine,ceq_type,ip_direction,bt_direction
      CHARACTER(24), DIMENSION(50) :: coil_name
      REAL(r8), DIMENSION(50,48) :: coil_cur
      REAL(r8), DIMENSION(50,48) :: coil_shiftx,coil_shifty,coil_shiftz
      REAL(r8), DIMENSION(50,48) :: coil_tiltx,coil_tilty,coil_tiltz
      REAL(r8), DIMENSION(50,48) :: coil_xnom,coil_ynom,coil_znom
      REAL(r8) :: cro,czo,cpsio,cpsilow,cpsilim,cqlim,
     $     ipd,btd,helicity

      LOGICAL :: gpec_interface, tilt_meters=.FALSE., coil_out=.FALSE.
      CHARACTER(256) :: data_dir = 'default'
      CHARACTER(512) :: cfile
      INTEGER, DIMENSION(:), POINTER :: cmfac

      TYPE :: coil_type
      CHARACTER(24) :: coil_name
      INTEGER :: ncoil,s,nsec
      REAL(r8) :: nw
      REAL(r8), DIMENSION(:,:), POINTER  :: x0, y0, z0, r_nom
      REAL(r8), DIMENSION(:,:), POINTER  :: tiltx, tilty, tiltz
      REAL(r8), DIMENSION(:), POINTER :: shiftx, shifty, shiftz
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
      REAL(r8), DIMENSION(50,48), OPTIONAL :: icoil_cur
      CHARACTER(24), DIMENSION(50), OPTIONAL :: icoil_name

      NAMELIST/coil_control/ceq_type,cmpsi,cmtheta,cmzeta,cmlow,cmhigh,
     $     data_dir,machine,ip_direction,bt_direction,
     $     coil_num,coil_name,coil_cur,
     $     coil_shiftx, coil_shifty, coil_shiftz,
     $     coil_tiltx, coil_tilty, coil_tiltz,
     $     coil_xnom, coil_ynom, coil_znom, cnpert
      NAMELIST/coil_output/gpec_interface, coil_out

      INTEGER :: ci,cj,ck,cl,cm,ci1,ci2,ci3,ci4,nsec
      REAL(r8) :: cr1,cr2
      REAL(r8) :: x0, y0, z0, x, y, z, r, r_nom, length
      REAL(r8) :: phi, angle, tiltx, tilty, tiltz
      REAL(r8) :: dx, dy, dz, dr
      REAL(r8) :: dtor = pi/180
      REAL(r8), DIMENSION(:), ALLOCATABLE :: dl, xm, ym, zm, rm
c-----------------------------------------------------------------------
c     initialize and read input data.
c-----------------------------------------------------------------------
      cmpsi=64
      cmtheta=480
      cmzeta=32
      cmlow=-64
      cmhigh=64
      coil_cur=0
      coil_shiftx=0
      coil_shifty=0
      coil_shiftz=0
      coil_tiltx=0
      coil_tilty=0
      coil_tiltz=0
      coil_xnom=1e9  ! large number defaults it to the center of mass
      coil_ynom=1e9
      coil_znom=1e9
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

      ! always output geometry if modifying it
      IF (SUM(SUM(ABS(coil_shiftx) + ABS(coil_shifty) + ABS(coil_shiftz)
     $    + ABS(coil_tiltx) + ABS(coil_tilty) + ABS(coil_tiltz), DIM=2))
     $    > 0) THEN
         coil_out = .TRUE.
      ENDIF
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
      cmtheta=MAXVAL((/480,cmtheta/))
      cmzeta=cnn*MAXVAL((/32,cmzeta/))
      cmpert=cmhigh-cmlow+1
      ALLOCATE(cmfac(cmpert))
      cmfac=(/(cm,cm=cmlow,cmhigh)/)
      ipd=1.0
      btd=1.0
      IF(ip_direction=="negative")ipd=-1.0
      IF(bt_direction=="negative")btd=-1.0
      helicity=ipd*btd
      gpec_interface=.TRUE.
      IF(cnpert < 0) cnpert = cnn ! default geometry manipulations correspond to n from DCON
c-----------------------------------------------------------------------
c     read coils for each machine.
c-----------------------------------------------------------------------
      IF(coil_num>0)ALLOCATE(coil(coil_num))

      DO ci=1, coil_num
         cfile=TRIM(data_dir)//"/"//TRIM(machine)//"_"//
     $         TRIM(coil_name(ci))//".dat"
         CALL ascii_open(coil_unit, cfile, "old")
         coil(ci)%coil_name=TRIM(coil_name(ci))
         READ(coil_unit,  *)
     $        coil(ci)%ncoil, coil(ci)%s, nsec, coil(ci)%nw
         coil(ci)%nsec = nsec
         ALLOCATE(coil(ci)%cur(coil(ci)%ncoil))
         ALLOCATE(coil(ci)%x(coil(ci)%ncoil, coil(ci)%s, nsec),
     $        coil(ci)%y(coil(ci)%ncoil, coil(ci)%s, nsec),
     $        coil(ci)%z(coil(ci)%ncoil, coil(ci)%s, nsec))
         ALLOCATE(coil(ci)%x0(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%y0(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%z0(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%r_nom(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%tiltx(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%tilty(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%tiltz(coil(ci)%ncoil, coil(ci)%s),
     $            coil(ci)%shiftx(coil(ci)%ncoil),
     $            coil(ci)%shifty(coil(ci)%ncoil),
     $            coil(ci)%shiftz(coil(ci)%ncoil) )
         ALLOCATE(dl(nsec), xm(nsec), ym(nsec), zm(nsec), rm(nsec))
         DO cj=1, coil(ci)%ncoil
            coil(ci)%cur(cj)=coil_cur(ci, cj)

            ! read original geometry
            DO ck=1, coil(ci)%s
               DO cl=1, nsec
                  READ(coil_unit, *) coil(ci)%x(cj, ck, cl),
     $                 coil(ci)%y(cj, ck, cl), coil(ci)%z(cj, ck, cl)
               ENDDO
            ENDDO

            ! inform user of any coil modifications
            IF(ABS(coil_shiftx(ci, cj)) + ABS(coil_shifty(ci, cj))
     $         +ABS(coil_shiftz(ci, cj)) > 0) WRITE(*,'(1x,a,i2)'),
     $         " > Shifting "//TRIM(coil_name(ci))//" coil ",cj
            IF(ABS(coil_tiltx(ci, cj)) + ABS(coil_tilty(ci, cj))
     $         +ABS(coil_tiltz(ci, cj)) > 0) THEN
                 IF(cnpert /= 0)THEN
                    WRITE(*,'(1x,a,i2)'),
     $         " > Tilting "//TRIM(coil_name(ci))//" coil",cj
                 ELSE
                    coil_tiltx = 0
                    coil_tilty = 0
                    coil_tiltz = 0
                    WRITE(*,'(1x,a)'), " >>> Ignoring tilt requests, "//
     $         "which have no n=0 component"
                 ENDIF
            ENDIF

            ! center of mass is default coil center if user enters super large numbers
            x0 = coil_xnom(ci, cj)
            y0 = coil_ynom(ci, cj)
            z0 = coil_znom(ci, cj)
            DO ck=1, coil(ci)%s
               dl = sqrt( (coil(ci)%x(cj, ck, 2:nsec)
     $                     - coil(ci)%x(cj, ck, 1:nsec-1)) ** 2
     $                   + (coil(ci)%y(cj, ck, 2:nsec)
     $                      - coil(ci)%y(cj, ck, 1:nsec-1)) ** 2
     $                   + (coil(ci)%z(cj, ck, 2:nsec)
     $                      - coil(ci)%z(cj, ck, 1:nsec-1)) ** 2
     $                  )
               length = sum(dl)
               xm = (coil(ci)%x(cj, ck, 2:nsec)
     $               + coil(ci)%x(cj, ck, 1:nsec-1))/2
               ym = (coil(ci)%y(cj, ck, 2:nsec)
     $               + coil(ci)%y(cj, ck, 1:nsec-1))/2
               zm = (coil(ci)%z(cj, ck, 2:nsec)
     $               + coil(ci)%z(cj, ck, 1:nsec-1))/2
               rm = SQRT(xm ** 2 + ym ** 2)
               r_nom = sum(rm * dl) / length  ! nominal radius of coil (makes most sense for PF coils)
               IF(ABS(x0) > 1e3) x0 = sum(xm * dl) / length
               IF(ABS(y0) > 1e3) y0 = sum(ym * dl) / length
               IF(ABS(z0) > 1e3) z0 = sum(zm * dl) / length
               ! apply tilts around coil center (user supplies this in degrees or meters)
               IF(tilt_meters)THEN
                  tiltx = ASIN(coil_tiltx(ci, cj) / r_nom)
                  tilty = ASIN(coil_tilty(ci, cj) / r_nom)
                  tiltz = ASIN(coil_tiltz(ci, cj) / r_nom)
               ELSE
                  tiltx = coil_tiltx(ci, cj) * dtor
                  tilty = coil_tilty(ci, cj) * dtor
                  tiltz = coil_tiltz(ci, cj) * dtor
               ENDIF
               IF(ABS(tiltx) + ABS(tilty) + ABS(tiltz)
     $            + ABS(coil_shiftx(ci, cj)) + ABS(coil_shifty(ci, cj))
     $            + ABS(coil_shiftz(ci, cj)) > 0)THEN ! only change the data if really necessary
                  DO cl=1, nsec
                     x = coil(ci)%x(cj, ck, cl) - x0
                     y = coil(ci)%y(cj, ck, cl) - y0
                     z = coil(ci)%z(cj, ck, cl) - z0
                     dx = x0
                     dy = y0
                     dz = z0
                     phi = ATAN2(y, x)
                     ! ztilt (clocking toroidally around the axis)
                     IF(tiltz /= 0)THEN
                        r = sqrt( y ** 2 + x ** 2)
                        angle = ATAN2(y, x) + tiltz*cos((cnpert-1)*phi)
                        dx = dx + (r * cos(angle) - x)
                        dy = dy + (r * sin(angle) - y)
                     ENDIF
                     ! ytilt (rotating around the y axis)
                     IF(tilty /= 0)THEN
                        r = sqrt( x ** 2 + z ** 2)
                        angle = ATAN2(x, z) + tilty*cos((cnpert-1)*phi)
                        dz = dz + (r * cos(angle) - z)
                        dx = dx + (r * sin(angle) - x)
                     ENDIF
                     ! xtilt (rotating around the x axis)
                     IF(tiltx /= 0)THEN
                        r = sqrt( z ** 2 + y ** 2)
                        angle = ATAN2(z, y) + tiltx*cos((cnpert-1)*phi)
                        dy = dy + (r * cos(angle) - y)
                        dz = dz + (r * sin(angle) - z)
                     ENDIF
                     coil(ci)%x(cj, ck, cl) = x + dx
                     coil(ci)%y(cj, ck, cl) = y + dy
                     coil(ci)%z(cj, ck, cl) = z + dz

                     ! apply shifts (elongations, etc.)
                     x = coil(ci)%x(cj, ck, cl) - x0
                     y = coil(ci)%y(cj, ck, cl) - y0
                     r = sqrt( y ** 2 + x ** 2)
                     dr = coil_shiftx(ci,cj) * cos(cnpert * phi)
     $                   +coil_shifty(ci,cj) * sin(cnpert * phi)
                     coil(ci)%x(cj, ck, cl) = x0 + (r + dr) * cos(phi)
                     coil(ci)%y(cj, ck, cl) = y0 + (r + dr) * sin(phi)
                     coil(ci)%z(cj, ck, cl) = coil(ci)%z(cj, ck, cl)
     $                                       +coil_shiftz(ci,cj)
                  ENDDO
               ENDIF
               ! record the subsystem's nominal center and modifications
               coil(ci)%x0(cj, ck) = x0
               coil(ci)%y0(cj, ck) = y0
               coil(ci)%z0(cj, ck) = z0
               coil(ci)%tiltx(cj, ck) = tiltx
               coil(ci)%tilty(cj, ck) = tilty
               coil(ci)%tiltz(cj, ck) = tiltz
               coil(ci)%r_nom(cj, ck) = r_nom
            ENDDO
            ! save shift info
            coil(ci)%shiftx(cj) = coil_shiftx(ci, cj)
            coil(ci)%shifty(cj) = coil_shifty(ci, cj)
            coil(ci)%shiftz(cj) = coil_shiftz(ci, cj)

         ENDDO
         CALL ascii_close(coil_unit)
         DEALLOCATE(dl, xm, ym, zm, rm)
      ENDDO
c-----------------------------------------------------------------------
c     Write final coil geometry information.
c-----------------------------------------------------------------------
      IF(coil_out) CALL write_coil_geometry(cnpert)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE coil_read
c-----------------------------------------------------------------------
c     subprogram 2. check.
c     Check status of netcdf file.
c-----------------------------------------------------------------------
      SUBROUTINE erchk(stat)
c-----------------------------------------------------------------------
c     declaration.
c-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: stat
c-----------------------------------------------------------------------
c     stop if it is an error.
c-----------------------------------------------------------------------
      IF(stat /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(stat))
         STOP "ERROR: failed to write/read netcdf file"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE erchk
c-----------------------------------------------------------------------
c     subprogram 3. coil_output.
c     Output the coil geometry information
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_coil_geometry(cnpert)

      INTEGER, INTENT(IN) :: cnpert

      LOGICAL :: debug_flag=.FALSE.
      INTEGER :: i, ci
      INTEGER :: ncid
      CHARACTER(512) :: cfile
      CHARACTER(2) :: sn

      INTEGER :: i_dim, i_id, n_dim, n_id

      INTEGER, DIMENSION(:), ALLOCATABLE :: c_dim,c_id,
     $         s_dim,s_id,
     $         p_dim,p_id, 
     $         x_id,y_id,z_id, 
     $         x0_id,y0_id,z0_id, 
     $         r_id,cur_id,w_id,
     $         xs_id,ys_id,zs_id, 
     $         xt_id,yt_id,zt_id

      ! label the toroidal mode number for consistent file naming convention
      IF (cnn<10) THEN
         WRITE(UNIT=sn,FMT='(I1)')cnn
         sn=ADJUSTL(sn)
      ELSE
         WRITE(UNIT=sn,FMT='(I2)')cnn
      ENDIF

      ! open netcdf file
      cfile = "gpec_coil_output_n"//TRIM(sn)//".nc"
      CALL erchk( nf90_create(cfile,
     $     cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

      ! define global file attributes
      IF(debug_flag) PRINT *," - Defining modal netcdf globals"
      CALL erchk( nf90_put_att(ncid,nf90_global,"title",
     $     "GPEC coil geometries"))

      ! define dimensions
      ALLOCATE(c_dim(coil_num),c_id(coil_num),
     $         s_dim(coil_num),s_id(coil_num),
     $         p_dim(coil_num),p_id(coil_num), 
     $         x_id(coil_num),y_id(coil_num),z_id(coil_num),
     $         x0_id(coil_num),y0_id(coil_num),z0_id(coil_num), 
     $         r_id(coil_num),cur_id(coil_num),w_id(coil_num),
     $         xs_id(coil_num),ys_id(coil_num),zs_id(coil_num), 
     $         xt_id(coil_num),yt_id(coil_num),zt_id(coil_num)
     $        )
      IF(debug_flag) PRINT *," - Defining dimensions in netcdf"
      CALL erchk( nf90_def_dim(ncid,"coil_index",coil_num,i_dim))
      CALL erchk( nf90_def_var(ncid, "coil_index", nf90_int,
     $            (/i_dim/), i_id) )
      CALL erchk( nf90_def_dim(ncid,"coil_strlen",24,n_dim))
      CALL erchk( nf90_def_var(ncid, "coil_name", nf90_char,
     $            (/n_dim, i_dim/), n_id) )
      DO ci=1, coil_num
         IF(debug_flag) PRINT *, ci, TRIM(coil_name(ci))
         CALL erchk( nf90_def_dim(ncid, TRIM(coil_name(ci))//
     $              "_coil", coil(ci)%ncoil, c_dim(ci)) )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//
     $              "_coil", nf90_int, c_dim(ci), c_id(ci)) )
         CALL erchk( nf90_def_dim(ncid, TRIM(coil_name(ci))//
     $              "_subsystem", coil(ci)%s, s_dim(ci)) )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//
     $              "_subsystem", nf90_int, s_dim(ci), s_id(ci)) )
         CALL erchk( nf90_def_dim(ncid, TRIM(coil_name(ci))//
     $              "_point", coil(ci)%nsec, p_dim(ci)) )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//
     $              "_point", nf90_int, p_dim(ci), p_id(ci)) )

         ! define variables
         IF(debug_flag) PRINT *," - Defining variables in netcdf"
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_x",
     $      nf90_double, (/c_dim(ci), s_dim(ci), p_dim(ci)/), x_id(ci)))
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_y",
     $      nf90_double, (/c_dim(ci), s_dim(ci), p_dim(ci)/), y_id(ci)))
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_z",
     $      nf90_double, (/c_dim(ci), s_dim(ci), p_dim(ci)/), z_id(ci)))
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_x0",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), x0_id(ci)) )
         CALL erchk( nf90_put_att(ncid,x0_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,x0_id(ci),"long_name",
     $               "Nominal x-center of coil") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_y0",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), y0_id(ci)) )
         CALL erchk( nf90_put_att(ncid,y0_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,y0_id(ci),"long_name",
     $               "Nominal y-center of coil") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_z0",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), z0_id(ci)) )
         CALL erchk( nf90_put_att(ncid,z0_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,z0_id(ci),"long_name",
     $               "Nominal z-center of coil") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_R_nom",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), r_id(ci)) )
         CALL erchk( nf90_put_att(ncid,r_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,r_id(ci),"long_name",
     $               "Nominal radius of coil") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_tiltx",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), xt_id(ci)) )
         CALL erchk( nf90_put_att(ncid,xt_id(ci),"units","degrees") )
         CALL erchk( nf90_put_att(ncid,xt_id(ci),"long_name",
     $               "Tilt about the x-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_tilty",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), yt_id(ci)) )
         CALL erchk( nf90_put_att(ncid,yt_id(ci),"units","degrees") )
         CALL erchk( nf90_put_att(ncid,yt_id(ci),"long_name",
     $               "Tilt about the y-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_tiltz",
     $      nf90_double, (/c_dim(ci), s_dim(ci)/), zt_id(ci)) )
         CALL erchk( nf90_put_att(ncid,zt_id(ci),"units","degrees") )
         CALL erchk( nf90_put_att(ncid,zt_id(ci),"long_name",
     $               "Tilt about the z-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_shiftx",
     $      nf90_double, c_dim(ci), xs_id(ci)) )
         CALL erchk( nf90_put_att(ncid,xs_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,xs_id(ci),"long_name",
     $               "Shift along x-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_shifty",
     $      nf90_double, c_dim(ci), ys_id(ci)) )
         CALL erchk( nf90_put_att(ncid,ys_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,ys_id(ci),"long_name",
     $               "Shift along y-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_shiftz",
     $      nf90_double, c_dim(ci), zs_id(ci)) )
         CALL erchk( nf90_put_att(ncid,zs_id(ci),"units","m") )
         CALL erchk( nf90_put_att(ncid,zs_id(ci),"long_name",
     $               "Shift along z-axis") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_current",
     $      nf90_double, c_dim(ci), cur_id(ci)) )
         CALL erchk( nf90_put_att(ncid,cur_id(ci),"units","A") )
         CALL erchk( nf90_def_var(ncid, TRIM(coil_name(ci))//"_nw",
     $      nf90_double, varid=w_id(ci)) )
         CALL erchk( nf90_put_att(ncid,w_id(ci),"long_name",
     $               "Number of windings") )
      ENDDO

      ! add attributes
      CALL erchk( nf90_put_att(ncid,nf90_global,'cnpert', cnpert))

      ! end definitions
      CALL erchk( nf90_enddef(ncid) )

      ! set variables
      CALL erchk( nf90_put_var(ncid, i_id, (/(i,i=1, coil_num)/)) )
      CALL erchk( nf90_put_var(ncid, n_id, coil_name(1:coil_num)) )
      DO ci=1, coil_num
         IF(debug_flag) PRINT *," - Putting dimension variables"//
     $      " in netcdf for "//TRIM(coil_name(ci))
         CALL erchk( nf90_put_var(ncid,c_id(ci),
     $      (/(i, i=1,coil(ci)%ncoil)/)) )
         CALL erchk( nf90_put_var(ncid,s_id(ci),
     $      (/(i, i=1,coil(ci)%s)/)) )
         CALL erchk( nf90_put_var(ncid,p_id(ci),
     $      (/(i, i=1,coil(ci)%nsec)/)) )

         CALL erchk( nf90_put_var(ncid,x_id(ci), coil(ci)%x) )
         CALL erchk( nf90_put_var(ncid,y_id(ci), coil(ci)%y) )
         CALL erchk( nf90_put_var(ncid,z_id(ci), coil(ci)%z) )

         CALL erchk( nf90_put_var(ncid,x0_id(ci), coil(ci)%x0) )
         CALL erchk( nf90_put_var(ncid,y0_id(ci), coil(ci)%y0) )
         CALL erchk( nf90_put_var(ncid,z0_id(ci), coil(ci)%z0) )

         CALL erchk( nf90_put_var(ncid,xs_id(ci), coil(ci)%shiftx) )
         CALL erchk( nf90_put_var(ncid,ys_id(ci), coil(ci)%shifty) )
         CALL erchk( nf90_put_var(ncid,zs_id(ci), coil(ci)%shiftz) )

         CALL erchk( nf90_put_var(ncid,xt_id(ci), coil(ci)%tiltx) )
         CALL erchk( nf90_put_var(ncid,yt_id(ci), coil(ci)%tilty) )
         CALL erchk( nf90_put_var(ncid,zt_id(ci), coil(ci)%tiltz) )

         CALL erchk( nf90_put_var(ncid,r_id(ci), coil(ci)%r_nom) )

         CALL erchk( nf90_put_var(ncid,cur_id(ci), coil(ci)%cur) )
         CALL erchk( nf90_put_var(ncid,w_id(ci), coil(ci)%nw) )
      ENDDO

      ! close file
      IF(debug_flag) PRINT *," - Closing netcdf file"
      CALL erchk( nf90_close(ncid) )
      DEALLOCATE( c_dim, c_id, s_dim,s_id, p_dim, p_id, 
     $         x_id, y_id, z_id,
     $         x0_id, y0_id, z0_id,  
     $         r_id, cur_id, w_id, 
     $         xs_id, ys_id, zs_id,  
     $         xt_id, yt_id, zt_id )

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_coil_geometry
c-----------------------------------------------------------------------
c     subprogram 4. coil_dealloc.
c     read coils.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE coil_dealloc

      INTEGER :: ci

      DO ci=1,coil_num
         DEALLOCATE(coil(ci)%x,coil(ci)%y,coil(ci)%z,
     $            coil(ci)%cur, coil(ci)%r_nom,
     $            coil(ci)%x0, coil(ci)%y0, coil(ci)%z0,
     $            coil(ci)%tiltx, coil(ci)%tilty, coil(ci)%tiltz,
     $            coil(ci)%shiftx, coil(ci)%shifty, coil(ci)%shiftz )
      ENDDO
      DEALLOCATE(coil)
      DEALLOCATE(cmfac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE coil_dealloc

      END MODULE coil_mod
