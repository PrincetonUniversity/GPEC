c-----------------------------------------------------------------------
c     file multi.f.
c     controls multiple runs of dcon.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. multi_mod.
c     1. multi_expand
c     2. multi_case.
c     3. multi_dir.
c     4. multi_bubble.
c     5. multi_insert.
c     6. multi_main.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. multi_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE multi_mod
      USE local_mod
      IMPLICIT NONE

      TYPE :: loop_type
      CHARACTER(16) :: varname,format,type
      INTEGER :: nvalue
      REAL(r8) :: value0,dvalue
      CHARACTER(128), DIMENSION(:), POINTER :: valname,dirname
      REAL(r8), DIMENSION(:), POINTER :: value
      END TYPE loop_type

      LOGICAL :: match_flag=.FALSE.
      CHARACTER(128) :: indir

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. multi_expand.
c     fills in various arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE multi_expand(loop)

      TYPE(loop_type), INTENT(INOUT) :: loop

      INTEGER :: iv,nv

      CHARACTER(128), DIMENSION(:), POINTER :: valname
      REAL(r8), DIMENSION(:), POINTER :: value
c-----------------------------------------------------------------------
c     start case select.
c-----------------------------------------------------------------------
      nv=loop%nvalue
      SELECT CASE(loop%type)
c-----------------------------------------------------------------------
c     additive and multiplicative cases.
c-----------------------------------------------------------------------
      CASE("add","mult")
         ALLOCATE(loop%value(0:nv),loop%valname(0:nv),
     $        loop%dirname(0:nv))
         ALLOCATE(value(0:nv),valname(0:nv))
         loop%value(0)=loop%value0
         WRITE(loop%valname(0),'('//loop%format//')')ABS(loop%value(0))
         DO iv=1,nv
            IF(loop%type == "add")THEN
               loop%value(iv)=loop%value(iv-1)+loop%dvalue
            ELSE
               loop%value(iv)=loop%value(iv-1)*loop%dvalue
            ENDIF
            WRITE(loop%valname(iv),'('//loop%format//')')
     $           ABS(loop%value(iv))
         ENDDO
         loop%valname=ADJUSTL(loop%valname)
         loop%dirname=loop%valname
         WHERE(loop%value < 0)
            loop%valname="-"//loop%valname
            loop%dirname="_"//loop%dirname
         END WHERE
c-----------------------------------------------------------------------
c     directory case.
c-----------------------------------------------------------------------
      CASE("dir")
         CALL multi_dir(indir,loop%format,loop%dirname,loop%valname)
         loop%nvalue=SIZE(loop%dirname)-1
c-----------------------------------------------------------------------
c     finish case select.
c-----------------------------------------------------------------------
      CASE("none")
      CASE DEFAULT
         CALL program_stop("Cannot recognize type "//TRIM(loop%type))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE multi_expand
c-----------------------------------------------------------------------
c     subprogram 2. multi_case.
c     runs and stores one case.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE multi_case(dirname)

      CHARACTER(*), INTENT(IN) :: dirname

      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     run dcon and store output.
c-----------------------------------------------------------------------
      CALL system("rm -f log.out")
      CALL system("dcon > log.out")
      IF(match_flag)THEN
         CALL system("match >> log.out")
         CALL system("rm -f euler.bin")
      ENDIF
      INQUIRE(FILE=dirname,EXIST=file_stat)
      IF(file_stat)CALL system("rm -fr "//TRIM(dirname))
      CALL system("mkdir -p "//TRIM(dirname))
      CALL system("cp equil.in "//TRIM(dirname))
      CALL system("cp dcon.in "//TRIM(dirname))
      CALL system("cp vac.in "//TRIM(dirname))
      IF(match_flag)CALL system("cp match.in "//TRIM(dirname))
      CALL system("rm -f ahg2msc.out")
      CALL system("rm -f mscvac.out")
      CALL system("mv *.out "//TRIM(dirname))
      CALL system("mv *.bin "//TRIM(dirname))
      CALL system("mv *.dat "//TRIM(dirname))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE multi_case
c-----------------------------------------------------------------------
c     subprogram 3. multi_dir.
c     gets and sorts directory listing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE multi_dir(rootname,format,dirname,valname)

      CHARACTER(*) :: rootname,format
      CHARACTER(*), DIMENSION(:), POINTER :: dirname,valname

      INTEGER :: ios,id,nd,nbad
      REAL(r8), DIMENSION(:), POINTER :: key
      INTEGER, DIMENSION(:), POINTER :: index
      CHARACTER(16) :: keyname
      CHARACTER(16), DIMENSION(:), POINTER :: tempname
c-----------------------------------------------------------------------
c     read directory.
c-----------------------------------------------------------------------
      CALL system("rm -f ls.out")
      CALL system("ls "//TRIM(rootname)//" | wc -l > ls.out")
      CALL system("ls "//TRIM(rootname)//" >> ls.out")
      CALL ascii_open(dir_unit,"ls.out","OLD")
      READ(dir_unit,*)nd
      ALLOCATE(tempname(nd),key(nd),index(nd))
      nbad=0
      DO id=1,nd
         READ(dir_unit,*)tempname(id)
         keyname=tempname(id)
         IF(keyname(1:1) == "_")keyname="-"//keyname(2:)
         READ(keyname,'('//format//')',IOSTAT=ios)key(id)
         IF(ios /= 0)THEN
            key(id)=HUGE(key(id))
            nbad=nbad+1
         ENDIF
         index(id)=id
      ENDDO
      CALL ascii_close(dir_unit)
      CALL system("rm ls.out")
c-----------------------------------------------------------------------
c     sort directory.
c-----------------------------------------------------------------------
      CALL multi_bubble(key,index)
      nd=nd-nbad-1
      ALLOCATE(dirname(0:nd),valname(0:nd))
      DO id=0,nd
         dirname(id)=TRIM(tempname(index(id+1)))
         valname(id)='"'//TRIM(indir)//"\/"
     $        //TRIM(tempname(index(id+1)))//'"'
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE multi_dir
c-----------------------------------------------------------------------
c     subprogram 4. multi_bubble.
c     performs a bubble sort in increasing order of value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE multi_bubble(key,index)

      REAL(r8), DIMENSION(:), INTENT(IN) :: key
      INTEGER, DIMENSION(:), INTENT(INOUT) :: index

      LOGICAL :: switch
      INTEGER :: i,temp,m
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      switch= .TRUE.
      m=SIZE(key)
      DO while(switch)
         switch= .FALSE.
         DO i=1,m-1
            IF(key(index(i)) > key(index(i+1)))THEN
               temp=index(i)
               index(i)=index(i+1)
               index(i+1)=temp
               switch= .TRUE.
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE multi_bubble
c-----------------------------------------------------------------------
c     subprogram 5. multi_insert.
c     inserts data into dcon namelist input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE multi_insert(loop,iv)

      TYPE(loop_type), INTENT(IN) :: loop
      INTEGER, INTENT(IN) :: iv

      CHARACTER(32) :: string
      CHARACTER(128) :: command
c-----------------------------------------------------------------------
c     modify loop value in input file.
c-----------------------------------------------------------------------
      string='=[a-zA-Z0-9+-._"/ ]*/'
      command="sed 's/"//TRIM(loop%varname)//TRIM(string)
     $     //TRIM(loop%varname)//"="//TRIM(loop%valname(iv))//"/' "
      CALL system("rm -f temp")
      CALL system(TRIM(command)//" equil.in > temp")
      CALL system("mv temp equil.in")
      CALL system(TRIM(command)//" dcon.in > temp")
      CALL system("mv temp dcon.in")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE multi_insert
      END MODULE multi_mod
c-----------------------------------------------------------------------
c     subprogram 6. multi_main.
c     main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM multi_main
      USE multi_mod
      IMPLICIT NONE

      CHARACTER(16) :: varname,format,type
      INTEGER :: nvalue
      REAL(r8) :: value0,dvalue
      
      INTEGER :: index_in,index_out
      CHARACTER(128) :: dirname,outdir
      TYPE(loop_type) :: inner,outer

      NAMELIST/dir_input/indir,outdir,match_flag
      NAMELIST/inner_input/varname,format,type,nvalue,value0,dvalue
      NAMELIST/outer_input/varname,format,type,nvalue,value0,dvalue
c-----------------------------------------------------------------------
c     read control file and fill loops.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"multi.in","OLD")
      READ(in_unit,NML=dir_input)
      READ(in_unit,NML=outer_input)
      outer%varname=varname
      outer%format=format
      outer%type=type
      outer%nvalue=nvalue
      outer%value0=value0
      outer%dvalue=dvalue
      CALL multi_expand(outer)
      READ(in_unit,NML=inner_input)
      inner%type=type
      IF(type /= "none")THEN
         inner%varname=varname
         inner%format=format
         inner%nvalue=nvalue
         inner%value0=value0
         inner%dvalue=dvalue
      CALL multi_expand(inner)
      ENDIF
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     create output directory, store multi.in.
c-----------------------------------------------------------------------
      CALL system("mkdir -p "//TRIM(outdir))
      CALL system("cp multi.in "//TRIM(outdir))
c-----------------------------------------------------------------------
c     loop over values.
c-----------------------------------------------------------------------
      DO index_out=0,outer%nvalue
         CALL multi_insert(outer,index_out)
         IF(inner%type == "none")THEN
            dirname=TRIM(outdir)//"/"//TRIM(outer%dirname(index_out))
            CALL multi_case(dirname)
         ELSE
            DO index_in=0,inner%nvalue
               CALL multi_insert(inner,index_in)
               dirname=TRIM(outdir)
     $              //"/"//TRIM(outer%dirname(index_out))
     $              //"/"//TRIM(inner%dirname(index_in))
               CALL multi_case(dirname)
            ENDDO
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END PROGRAM multi_main
