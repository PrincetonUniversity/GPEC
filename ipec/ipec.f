c-----------------------------------------------------------------------
c     IDEAL PERTURBED EQUILIBRIUM CONTROL
c     IPEC: main program
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     ipec_main.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM ipec_main
      USE ipdiag_mod
      USE ipout_mod

      IMPLICIT NONE

      INTEGER :: in,osing,rstep1,rstep2,resol,angnum,labl,infnum,
     $     mthnumb,meas,modem,label,nr,nz,modemin,modemax,
     $     errnmin,errnmax,errmmax,poloin,toroin,poloout,toroout
      REAL(r8) :: majr,minr
      REAL(r8) :: rdist,bigstep,smlstep,smallwidth,factor,
     $     wegtfac1,wegtfac2,fp,normpsi,lportion,hportion
      CHARACTER(128) :: infile
      LOGICAL :: erdata_flag,mode_flag,response_flag,singcoup_flag,
     $     singfld_flag,pmodb_flag,xbnorm_flag,xbrzphi_flag,
     $     nrzeq_flag,extp_flag,extt_flag,
     $     xbcontra_flag,energy_flag,respmat_flag,angles_flag,
     $     m3d_flag,svdvec_flag,d3_flag
      COMPLEX(r8), DIMENSION(:), POINTER :: bexmn,bermn,
     $     brrmn,bnomn,bpamn,fxmn,xwpmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: obextmn,ocextmn

      NAMELIST/ipec_input/ieqfile,idconfile,ivacuumfile,
     $     power_flag,fft_flag,mthsurf0,
     $     erdata_flag,errtype,errnmin,errnmax,errmmax,
     $     poloin,toroin,infnum,infiles,mode_flag,modemin,modemax
      NAMELIST/ipec_control/response_flag,
     $     dist,bdist,maxdbratio,rstep,modelnum
      NAMELIST/ipec_output/singcoup_flag,
     $     singfld_flag,pmodb_flag,xbnorm_flag,
     $     poloout,toroout,labl,
     $     xbrzphi_flag,nrzeq_flag,nr,nz,extp_flag,extt_flag
      NAMELIST/ipec_diagnose/xbcontra_flag,energy_flag,respmat_flag,
     $     angles_flag,m3d_flag
c-----------------------------------------------------------------------
c     read ipec.in.
c-----------------------------------------------------------------------
      CALL ascii_open(in_unit,"ipec.in","OLD")
      READ(in_unit,NML=ipec_input)
      READ(in_unit,NML=ipec_control)  
      READ(in_unit,NML=ipec_output)
      READ(in_unit,NML=ipec_diagnose)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     assign temporal values.
c-----------------------------------------------------------------------
      lmlow=-errmmax
      lmhigh=errmmax
      IF (response_flag) THEN
         resp=1
      ENDIF
c-----------------------------------------------------------------------
c     check time.
c-----------------------------------------------------------------------
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
c-----------------------------------------------------------------------
c     prepare for ideal solutions.
c-----------------------------------------------------------------------
      CALL idcon_read
      CALL idcon_transform
c-----------------------------------------------------------------------
c     reconstruct metric tensors.
c-----------------------------------------------------------------------
      CALL idcon_metric
c-----------------------------------------------------------------------
c     read vacuum data.
c-----------------------------------------------------------------------
      CALL idcon_vacuum
c-----------------------------------------------------------------------
c     compute plasma response.
c-----------------------------------------------------------------------
      CALL ipresp_eigen
      CALL ipresp_pinduct
      CALL ipresp_sinduct
      CALL ipresp_permeab
      CALL ipresp_reluct
c-----------------------------------------------------------------------
c     define variables and assign temporal values.
c-----------------------------------------------------------------------
      ALLOCATE(bexmn(lmpert),bermn(lmpert),
     $     brrmn(mpert),bnomn(mpert),bpamn(mpert),
     $     fxmn(mpert),xwpmn(mpert))
      IF (nrzeq_flag) THEN
         nr=mr
         nz=mz
      ENDIF
c-----------------------------------------------------------------------
c     full analysis.
c-----------------------------------------------------------------------
      IF (response_flag) THEN
         CALL ipout_resp(10)
      ENDIF
      IF (singcoup_flag) THEN
         CALL ipout_singcoup(dist,poloout,toroout)
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and external field.
c-----------------------------------------------------------------------
      IF (erdata_flag) THEN
         DO in=1,infnum
            label=labl+in-1
            infile=infiles(in)
            IF ((poloin /= toroin).OR.(poloout /= toroout)) THEN
               edge_flag=.TRUE.
               CALL ipdiag_extfld(infile,errtype,bexmn,
     $              poloin,toroin,poloout,toroout,bermn,label)
               brrmn=bermn(mlow-lmlow+1:mhigh-lmlow+1)
               edge_flag=.FALSE.
               CALL ipout_errfld(infile,errtype,brrmn,
     $              poloout,toroout,resp,bnomn,xwpmn,label)
               edge_flag=.TRUE.
            ELSE 
               edge_flag=.TRUE.
               CALL ipout_errfld(infile,errtype,brrmn,
     $              poloout,toroout,resp,bnomn,xwpmn,label)
            ENDIF
            IF (singfld_flag) THEN
               CALL ipout_singfld(0,xwpmn,dist,msing,label)
            ENDIF
            IF (pmodb_flag) THEN
               CALL ipout_pmodb(0,xwpmn,poloout,toroout,label)
            ENDIF
c-----------------------------------------------------------------------
c     additional information.
c-----------------------------------------------------------------------
            IF (xbnorm_flag) THEN
               CALL ipout_xbnorm(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbrzphi_flag) THEN
               CALL ipeq_rzpgrid(nr,nz)
               CALL ipout_xbrzphi(0,xwpmn,nr,nz,label)
            ENDIF
            IF (extt_flag) THEN
               CALL ipvacuum_bnormal(psilim,bnomn,
     $              poloout,toroout,nr,nz,'actual',label)
            ENDIF
            IF (extp_flag) THEN
               bpamn=bnomn-brrmn
               CALL ipvacuum_bnormal(psilim,bpamn,
     $              poloout,toroout,nr,nz,'plasma',label)
            ENDIF
            IF (xbcontra_flag) THEN
               CALL ipdiag_xbcontra(0,xwpmn,label)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and mode field.
c-----------------------------------------------------------------------
      IF (mode_flag) THEN
         infnum=modemax-modemin+1
         DO in=1,infnum
            label=labl+in-1
            edge_flag=.FALSE.
            brrmn=0
            brrmn(modemin+in-mlow)=1e-4
            CALL ipout_errfld(infile,errtype,brrmn,
     $           poloout,toroout,resp,bnomn,xwpmn,label)
            edge_flag=.TRUE.
            IF (singfld_flag) THEN
               CALL ipout_singfld(0,xwpmn,dist,msing,label)
            ENDIF
            IF (pmodb_flag) THEN
               CALL ipout_pmodb(0,xwpmn,poloout,toroout,label)
            ENDIF
c-----------------------------------------------------------------------
c     additional information.
c-----------------------------------------------------------------------
            IF (xbnorm_flag) THEN
               CALL ipout_xbnorm(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbrzphi_flag) THEN
               CALL ipeq_rzpgrid(nr,nz)
               CALL ipout_xbrzphi(0,xwpmn,nr,nz,label)
            ENDIF
            IF (extt_flag) THEN
               CALL ipvacuum_bnormal(psilim,bnomn,
     $              poloout,toroout,nr,nz,'actual',label)
            ENDIF
            IF (extp_flag) THEN
               bpamn=bnomn-brrmn
               CALL ipvacuum_bnormal(psilim,bpamn,
     $              poloout,toroout,nr,nz,'plasma',label)
            ENDIF
            IF (xbcontra_flag) THEN
               CALL ipdiag_xbcontra(0,xwpmn,label)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF (energy_flag) THEN
         CALL ipdiag_energy
      ENDIF
      IF (respmat_flag) THEN
         CALL ipdiag_respmat
      ENDIF
      IF (angles_flag) THEN
         CALL ipdiag_angles
      ENDIF
      IF (m3d_flag) THEN
         normpsi=1.0/0.5
         fp=1e-3
         
         DO in=-10,10
            label=in+10
            fxmn=0
            fxmn(in-mlow+1)=fp*normpsi
            CALL ipeq_cotoha(psilim,fxmn,0,0)
            fxmn=-twopi*ifac*chi1*(mfac-nn*qlim)*fxmn
            CALL ipeq_weight(psilim,fxmn,0)
            edge_flag=.FALSE.
            CALL ipout_errfld(infile,errtype,fxmn,
     $           poloout,toroout,resp,bnomn,xwpmn,label)
            edge_flag=.TRUE.
            CALL ipout_singfld(0,xwpmn,dist,msing,label)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL ipec_dealloc
      CALL ipec_stop("Normal termination.")
      END PROGRAM ipec_main
