c-----------------------------------------------------------------------
c     file sum.f.
c     extracts data from multiple runs of dcon.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. sum_mod.
c     1. sum_case1.
c     2. sum_case2.
c     3. sum_dir.
c     4. sum_bubble.
c     5. sum_deltap.
c     6. sum_fixfiles.
c     7. sum_asinh.
c     8. sum_log.
c     9. sum_main.
c-----------------------------------------------------------------------
c     subprogram 0. sum_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE sum_mod
      USE local_mod
      IMPLICIT NONE

      TYPE :: loop_type
      CHARACTER(16) :: varname,format,type
      INTEGER :: nvalue
      REAL(r4) :: value0,dvalue
      CHARACTER(256), DIMENSION(:), POINTER :: dirname
      REAL(r4), DIMENSION(:), POINTER :: value
      END TYPE loop_type

      TYPE :: sing_type
      REAL(r4) :: psifac,q,di,deltap1,deltap2
      END TYPE sing_type

      TYPE :: node_type
      INTEGER :: msing
      REAL(r4) :: inval,outval
      TYPE(sing_type), DIMENSION(:), POINTER :: sing
      END TYPE node_type

      CHARACTER(128) :: indir
      TYPE(loop_type) :: inner,outer

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. sum_case1.
c     runs and stores one case.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sum_case1(dirname,outval,inval)

      CHARACTER(*), INTENT(IN) :: dirname
      REAL(r4), INTENT(IN) :: outval,inval

      CHARACTER(128) :: filename
      INTEGER :: mpsi,mtheta,mlow,mhigh,mpert,mband
      REAL(r4) :: psilow,psihigh,amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,ro,zo,psio,betap1,betap2,betap3,betat,betan,
     $     bt0,q0,qmin,qmax,qa,crnt,plasma1,vacuum1,total1
c-----------------------------------------------------------------------
c     read sum1.dat.
c-----------------------------------------------------------------------
      filename=TRIM(dirname)//"/sum1.dat"
      CALL bin_open(in_unit,TRIM(filename),"OLD","REWIND")
      READ(in_unit)mpsi,mtheta,mlow,mhigh,mpert,mband,
     $     psilow,psihigh,amean,rmean,aratio,kappa,delta1,delta2,
     $     li1,li2,li3,ro,zo,psio,betap1,betap2,betap3,betat,betan,
     $     bt0,q0,qmin,qmax,qa,crnt,plasma1,vacuum1,total1
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     write graphics file.
c-----------------------------------------------------------------------
      WRITE(sum1_unit)
     $     sum_log(outval,outer%type),sum_log(inval,inner%type),
     $     psilow,psihigh,amean,rmean,aratio,kappa,delta1,delta2,li1,
     $     li2,li3,ro,zo,psio,betap1,betap2,betap3,betat,betan,bt0,q0,
     $     qmin,qmax,qa,crnt,plasma1,vacuum1,total1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sum_case1
c-----------------------------------------------------------------------
c     subprogram 2. sum_case2.
c     runs and stores one case.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sum_case2(dirname,outval,inval,node)

      CHARACTER(*), INTENT(IN) :: dirname
      REAL(r4), INTENT(IN) :: outval,inval
      TYPE(node_type), INTENT(OUT) :: node

      CHARACTER(128) :: filename
      INTEGER :: ising
c-----------------------------------------------------------------------
c     read sum2.dat.
c-----------------------------------------------------------------------
      node%outval=outval
      node%inval=inval
      filename=TRIM(dirname)//"/sum2.dat"
      CALL bin_open(in_unit,TRIM(filename),"OLD","REWIND")
      READ(in_unit)node%msing
      ALLOCATE(node%sing(node%msing))
      READ(in_unit)
     $     (node%sing(ising)%psifac,ising=1,node%msing),
     $     (node%sing(ising)%q,ising=1,node%msing),
     $     (node%sing(ising)%di,ising=1,node%msing),
     $     (node%sing(ising)%deltap1,ising=1,node%msing),
     $     (node%sing(ising)%deltap2,ising=1,node%msing)
      CALL bin_close(in_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sum_case2
c-----------------------------------------------------------------------
c     subprogram 3. sum_dir.
c     gets and sorts directory listing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sum_dir(rootname,format,dirname,value)

      CHARACTER(*) :: rootname,format
      CHARACTER(*), DIMENSION(:), POINTER :: dirname
      REAL(r4), DIMENSION(:), POINTER :: value

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: ios,id,nd,nbad
      REAL(r4), DIMENSION(:), POINTER :: key
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
      CALL sum_bubble(key,index)
      nd=nd-nbad-1
      ALLOCATE(dirname(0:nd),value(0:nd))
      DO id=0,nd
         dirname(id)=TRIM(tempname(index(id+1)))
         value(id)=key(index(id+1))
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(*,'(i3,2x,a,2x,'//TRIM(format)//')')
     $        (id,TRIM(dirname(id)),value(id),id=0,nd)
         CALL program_stop("Termination by sum_dir")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      DEALLOCATE(tempname,key,index)
      RETURN
      END SUBROUTINE sum_dir
c-----------------------------------------------------------------------
c     subprogram 4. sum_bubble.
c     performs a bubble sort in order of increasing value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sum_bubble(key,index)

      REAL(r4), DIMENSION(:), INTENT(IN) :: key
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
      END SUBROUTINE sum_bubble
c-----------------------------------------------------------------------
c     subprogram 5. sum_deltap.
c     writes binary files for drawing deltap.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE sum_deltap(node)

      TYPE(node_type), DIMENSION(0:,0:), INTENT(IN), TARGET :: node

      INTEGER :: i_in,i_out,m_in,m_out,mq,iq,ising
      INTEGER, PARAMETER :: nq=100
      REAL(r4), DIMENSION(nq) :: q
      TYPE(node_type), POINTER :: np
      TYPE(sing_type), POINTER :: sp
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      m_out=SIZE(node,1)-1
      m_in=SIZE(node,2)-1
      q=HUGE(q)
      mq=0
c-----------------------------------------------------------------------
c     count distinct singular values.
c-----------------------------------------------------------------------
      DO i_out=0,m_out
         DO i_in=0,m_in
            DO ising=1,node(i_out,i_in)%msing
               iq=0
               DO
                  iq=iq+1
                  IF(iq > mq .OR.
     $                 node(i_out,i_in)%sing(ising)%q == q(iq))EXIT
               ENDDO
               IF(iq > mq)THEN
                  mq=iq
                  q(iq)=node(i_out,i_in)%sing(ising)%q
               ENDIF
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write binary data.
c-----------------------------------------------------------------------
      CALL bin_open(sum2_unit,"sum2.bin","UNKNOWN","REWIND")
      DO iq=1,mq
         DO i_out=0,m_out
            DO i_in=0,m_in
               np => node(i_out,i_in)
               DO ising=1,np%msing
                  sp => np%sing(ising)
                  IF(sp%q /= q(iq))CYCLE
                  WRITE(sum2_unit)sum_log(np%outval,outer%type),
     $                 sum_log(np%inval,inner%type),sp%psifac,sp%di,
     $                 sp%deltap1,sum_asinh(sp%deltap1)
               ENDDO
            ENDDO
         ENDDO
         WRITE(sum2_unit)
      ENDDO
      CALL bin_close(sum2_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sum_deltap
c-----------------------------------------------------------------------
c     subprogram 6. sum_fixfiles.
c     modifies draw.in files.
c-----------------------------------------------------------------------
      SUBROUTINE sum_fixfiles
c-----------------------------------------------------------------------
c     modify variable names.
c-----------------------------------------------------------------------
      IF(inner%type == "mult")inner%varname="log10 "//inner%varname
      IF(outer%type == "mult")outer%varname="log10 "//outer%varname
c-----------------------------------------------------------------------
c     modify drawsum1.in.
c-----------------------------------------------------------------------
      IF(outer%type /= "none")THEN
         CALL system("sed '9,10s/0\	[ a-zA-Z0-9]*/"
     $        //"0\	"//TRIM(outer%varname)//"/'"
     $        //" drawsum1.in > temp")
         CALL system("mv temp drawsum1.in")
      ENDIF
      IF(inner%type /= "none")THEN
         CALL system("sed '9,10s/1\	[ a-zA-Z0-9]*/"
     $        //"1\	"//TRIM(inner%varname)//"/'"
     $        //" drawsum1.in > temp")
         CALL system("mv temp drawsum1.in")
      ENDIF
c-----------------------------------------------------------------------
c     modify drawsum2.in.
c-----------------------------------------------------------------------
      IF(outer%type /= "none")THEN
         CALL system("sed '9,10s/0\	[ a-zA-Z0-9]*/"
     $        //"0\	"//TRIM(outer%varname)//"/'"
     $        //" drawsum2.in > temp")
         CALL system("mv temp drawsum2.in")
      ENDIF
      IF(inner%type /= "none")THEN
         CALL system("sed '9,10s/1\	[ a-zA-Z0-9]*/"
     $        //"1\	"//TRIM(inner%varname)//"/'"
     $        //" drawsum2.in > temp")
         CALL system("mv temp drawsum2.in")
      ENDIF
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE sum_fixfiles
c-----------------------------------------------------------------------
c     subprogram 7. sum_asinh.
c     computes inverse hyperbolic sin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION sum_asinh(x) RESULT(y)

      REAL(r4), INTENT(IN) :: x
      REAL(r4) :: y

      REAL(r4) :: x2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      x2=x*x
      IF(x < -1e3)THEN
         y=-LOG(-2.*x)-.25/x2
      ELSE IF(ABS(x) < 1e-2)THEN
         y=x*(1-x2/6*(1-x2*9/20))
      ELSE IF(x > 1e3)THEN
         y=LOG(2.*x)+.25/x2
      ELSE
         y=LOG(x+SQRT(1.+x2))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION sum_asinh
c-----------------------------------------------------------------------
c     subprogram 8. sum_log.
c     takes log of argument if appropriate.
c-----------------------------------------------------------------------
      FUNCTION sum_log(x,type) RESULT(y)

      REAL(r4), INTENT(IN) :: x
      CHARACTER(*), INTENT(IN) :: type
      REAL(r4) :: y
c-----------------------------------------------------------------------
c     compute outval1 and inval1.
c-----------------------------------------------------------------------
      IF(type == "mult")THEN
         y=LOG10(ABS(x))
      ELSE
         y=x
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION sum_log
      END MODULE sum_mod
c-----------------------------------------------------------------------
c     subprogram 9. sum_main.
c     main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM sum_main
      USE sum_mod
      IMPLICIT NONE

      CHARACTER(16) :: varname,format,type
      INTEGER :: nvalue
      REAL(r4) :: value0,dvalue
      
      LOGICAL :: sum2_flag=.FALSE.
      INTEGER :: index_in,index_out
      CHARACTER(128) :: dirname,filename,rootname
      TYPE(node_type), DIMENSION(:,:), POINTER :: node
      
      NAMELIST/sum_input/rootname,sum2_flag
      NAMELIST/inner_input/varname,format,type,nvalue,value0,dvalue
      NAMELIST/outer_input/varname,format,type,nvalue,value0,dvalue
c-----------------------------------------------------------------------
c     read sum.in.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"sm.in","OLD")
      READ(in_unit,NML=sum_input)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     read multi.in.
c-----------------------------------------------------------------------
      filename=TRIM(rootname)//"/multi.in"
      CALL ascii_open(in_unit,TRIM(filename),"OLD")
      READ(in_unit,NML=outer_input)
      outer%varname=varname
      outer%format=format
      outer%type=type
      READ(in_unit,NML=inner_input)
      inner%type=type
      IF(type /= "none")THEN
         inner%varname=varname
         inner%format=format
      ENDIF
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      CALL bin_open(sum1_unit,"sum1.bin","UNKNOWN","REWIND")
c-----------------------------------------------------------------------
c     start loop over outer values.
c-----------------------------------------------------------------------
      dirname=rootname
      CALL sum_dir(dirname,outer%format,outer%dirname,outer%value)
      outer%nvalue=SIZE(outer%value)-1
      IF(inner%type == "none")THEN
         ALLOCATE(node(0:outer%nvalue,0:0))
      ELSE
         dirname=TRIM(rootname)//"/"//TRIM(outer%dirname(index_out))
         CALL sum_dir(dirname,inner%format,inner%dirname,inner%value)
         inner%nvalue=SIZE(inner%value)-1
         ALLOCATE(node(0:outer%nvalue,0:inner%nvalue))
      ENDIF
      DO index_out=0,outer%nvalue
         dirname=TRIM(rootname)//"/"//TRIM(outer%dirname(index_out))
c-----------------------------------------------------------------------
c     get outer values if no inner loop.
c-----------------------------------------------------------------------
         IF(inner%type == "none")THEN
            CALL sum_case1(dirname,outer%value(index_out),0._r4)
            IF(sum2_flag)CALL sum_case2(dirname,outer%value(index_out),
     $           0._r4,node(index_out,0))
c-----------------------------------------------------------------------
c     loop over inner values.
c-----------------------------------------------------------------------
         ELSE
            CALL sum_dir(dirname,inner%format,inner%dirname,inner%value)
            inner%nvalue=SIZE(inner%value)-1
            DO index_in=0,inner%nvalue
               dirname=TRIM(rootname)
     $              //"/"//TRIM(outer%dirname(index_out))
     $              //"/"//TRIM(inner%dirname(index_in))
               CALL sum_case1(dirname,outer%value(index_out),
     $              inner%value(index_in))
               IF(sum2_flag)CALL sum_case2(dirname,
     $              outer%value(index_out),inner%value(index_in),
     $              node(index_out,index_in))
            ENDDO
            WRITE(sum1_unit)
            WRITE(sum2_unit)
         ENDIF
      ENDDO
      IF(sum2_flag)CALL sum_deltap(node)
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
      WRITE(sum1_unit)
      CALL bin_close(sum1_unit)
      CALL sum_fixfiles
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      IF(TRIM(inner%type) /= "none")
     $     DEALLOCATE(inner%dirname,inner%value)
      IF(TRIM(outer%type) /= "none")
     $     DEALLOCATE(outer%dirname,outer%value)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END PROGRAM sum_main
