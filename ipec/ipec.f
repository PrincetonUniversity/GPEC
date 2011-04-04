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

      INTEGER :: i,in,osing,resol,angnum,labl,infnum,left,wegt,
     $     mthnumb,meas,modem,label,nr,nz,modemin,modemax,
     $     poloin,toroin,poloout,toroout,lowmode,highmode,
     $     m3low,m3high
      INTEGER, DIMENSION(:), POINTER :: ipiv
      REAL(r8) :: majr,minr,scale,rdist,smallwidth,factor,fp,normpsi
      CHARACTER(128) :: infile,formattype
      LOGICAL :: erdata_flag,harmonic_flag,mode_flag,response_flag,
     $     singcoup_flag,singfld_flag,pmodb_flag,nrzeq_flag,
     $     xbcontra_flag,xbnormal_flag,xbnovc_flag,xbnobo_flag,
     $     d3_flag,xbnorm_flag,pmodbst_flag,pmodbrz_flag,rzphibx_flag,
     $     radvar_flag,eigen_flag,magpot_flag,energy_flag,respmat_flag,
     $     arbsurf_flag,angles_flag,surfmode_flag,rzpgrid_flag,
     $     singcurs_flag,m3d_flag,cas3d_flag,test_flag
      COMPLEX(r8), DIMENSION(:), POINTER :: bexmn,bermn,
     $     brrmn,bnomn,fxmn,xwpmn
      COMPLEX(r8), DIMENSION(:,:), POINTER :: invmats,temp1

      NAMELIST/ipec_input/ieqfile,idconfile,ivacuumfile,
     $     power_flag,fft_flag,mthsurf0,left,scale,wegt,
     $     erdata_flag,formattype,errnmin,errnmax,errmmin,errmmax,
     $     poloin,toroin,infnum,infiles,
     $     harmonic_flag,mode_flag,modemin,modemax,eqoff_flag
      NAMELIST/ipec_control/response_flag,dist,bdist,modelnum
      NAMELIST/ipec_output/singcoup_flag,nrzeq_flag,nr,nz,labl,
     $     singfld_flag,pmodb_flag,rstep,poloout,toroout,
     $     eqbrzphi_flag,brzphi_flag,xrzphi_flag,
     $     vbrzphi_flag,vpbrzphi_flag,vvbrzphi_flag,divzero_flag
      NAMELIST/ipec_diagnose/singcurs_flag,xbcontra_flag,xbnormal_flag,
     $     xbnovc_flag,xbnobo_flag,d3_flag,xbnorm_flag,div_flag,
     $     pmodbst_flag,pmodbrz_flag,rzphibx_flag,radvar_flag,
     $     eigen_flag,magpot_flag,energy_flag,respmat_flag,
     $     arbsurf_flag,majr,minr,angles_flag,surfmode_flag,
     $     lowmode,highmode,rzpgrid_flag,m3d_flag,m3low,m3high,
     $     cas3d_flag,test_flag,resol,smallwidth
c-----------------------------------------------------------------------
c     read ipec.in.
c-----------------------------------------------------------------------
      WRITE(*,*)"Starting ipec calculations - v1.0"
      CALL ascii_open(in_unit,"ipec.in","OLD")
      READ(in_unit,NML=ipec_input)
      READ(in_unit,NML=ipec_control)  
      READ(in_unit,NML=ipec_output)
      READ(in_unit,NML=ipec_diagnose)
      CALL ascii_close(in_unit)
c-----------------------------------------------------------------------
c     assign temporal values.
c-----------------------------------------------------------------------
      lmlow=errmmin
      lmhigh=errmmax
      IF (response_flag) resp=1
      IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $     vbrzphi_flag .OR. vpbrzphi_flag .OR. vvbrzphi_flag) psixy=1
      IF (eqoff_flag) psixy=0
c-----------------------------------------------------------------------
c     check time.
c-----------------------------------------------------------------------
      CALL DATE_AND_TIME(date,time,zone,values)
      seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
c-----------------------------------------------------------------------
c     prepare for ideal solutions.
c-----------------------------------------------------------------------
      CALL idcon_read(psixy)
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
c-----------------------------------------------------------------------
c     define variables and assign temporal values.
c-----------------------------------------------------------------------
      ALLOCATE(bexmn(lmpert),bermn(lmpert),
     $     brrmn(mpert),bnomn(mpert),fxmn(mpert),xwpmn(mpert))
c-----------------------------------------------------------------------
c     full analysis.
c-----------------------------------------------------------------------
      IF (response_flag) THEN
         CALL ipout_resp
      ENDIF
      IF (singcoup_flag) THEN
         CALL ipout_singcoup(dist,poloout,toroout)
      ENDIF
      IF (rstep .EQ. 0) rstep=mstep
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and external field.
c-----------------------------------------------------------------------
      IF (erdata_flag) THEN
         DO in=1,infnum
            label=labl+in-1
            infile=infiles(in)
            edge_flag=.TRUE.
            CALL ipout_errfld(infile,formattype,left,scale,brrmn,
     $           poloin,toroin,wegt,resp,bnomn,xwpmn,label)
            IF (singfld_flag) THEN
               CALL ipout_singfld(0,xwpmn,dist,poloout,toroout,label)
            ENDIF
            IF (pmodb_flag) THEN
               CALL ipout_pmodb(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $           vbrzphi_flag .OR. vpbrzphi_flag .OR. vvbrzphi_flag)
     $           THEN
               IF (nrzeq_flag) THEN
                  nr=mr
                  nz=mz
               ENDIF
               IF (.NOT. eqoff_flag) CALL ipeq_rzpgrid(nr,nz)
               CALL ipout_xbrzphi(0,xwpmn,nr,nz,brrmn,bnomn,label)
            ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
            IF (singcurs_flag) THEN
               CALL ipdiag_singcurs(0,xwpmn,msing,resol,
     $              smallwidth,label)
            ENDIF
            IF (xbcontra_flag) THEN
               CALL ipdiag_xbcontra(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnormal_flag) THEN
               CALL ipdiag_xbnormal(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnovc_flag) THEN
               CALL ipdiag_xbnovc(0,xwpmn,label)
            ENDIF
            IF (xbnobo_flag) THEN
               CALL ipdiag_xbnobo(0,xwpmn,d3_flag,label)
            ENDIF
            IF (xbnorm_flag) THEN
               CALL ipdiag_xbnorm(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (pmodbst_flag) THEN
               CALL ipdiag_pmodbst(0,xwpmn,label)
            ENDIF
            IF (pmodbrz_flag) THEN
               CALL ipdiag_pmodbrz(0,xwpmn,label)
            ENDIF            
            IF (rzphibx_flag) THEN
               CALL ipdiag_rzphibx(0,xwpmn,label)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with a given equilibrium and mode field.
c-----------------------------------------------------------------------
      IF (harmonic_flag) THEN
         infnum=modemax-modemin+1
         DO in=1,infnum
            label=labl+in-1
            edge_flag=.FALSE.
            brrmn=0
            brrmn(modemin+in-mlow)=1e-4
            CALL ipout_errfld(infile,errtype,left,scale,brrmn,
     $           poloin,toroin,wegt,resp,bnomn,xwpmn,label)
            edge_flag=.TRUE.
            IF (singfld_flag) THEN
               CALL ipout_singfld(0,xwpmn,dist,poloout,toroout,label)
            ENDIF
            IF (pmodb_flag) THEN
               CALL ipout_pmodb(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $           vbrzphi_flag .OR. vpbrzphi_flag .OR. vvbrzphi_flag)
     $           THEN
               IF (nrzeq_flag) THEN
                  nr=mr
                  nz=mz
               ENDIF
               IF (.NOT. eqoff_flag) CALL ipeq_rzpgrid(nr,nz)
               CALL ipout_xbrzphi(0,xwpmn,nr,nz,brrmn,bnomn,label)
            ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
            IF (singcurs_flag) THEN
               CALL ipdiag_singcurs(0,xwpmn,msing,resol,
     $              smallwidth,label)
            ENDIF
            IF (xbcontra_flag) THEN
               CALL ipdiag_xbcontra(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnormal_flag) THEN
               CALL ipdiag_xbnormal(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnovc_flag) THEN
               CALL ipdiag_xbnovc(0,xwpmn,label)
            ENDIF
            IF (xbnobo_flag) THEN
               CALL ipdiag_xbnobo(0,xwpmn,d3_flag,label)
            ENDIF
            IF (xbnorm_flag) THEN
               CALL ipdiag_xbnorm(0,xwpmn,poloout,toroout,label)
            ENDIF
            IF (pmodbst_flag) THEN
               CALL ipdiag_pmodbst(0,xwpmn,label)
            ENDIF
            IF (pmodbrz_flag) THEN
               CALL ipdiag_pmodbrz(0,xwpmn,label)
            ENDIF
            IF (rzphibx_flag) THEN
               CALL ipdiag_rzphibx(0,xwpmn,label)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     perturbed equilibria with an eigenmode.
c-----------------------------------------------------------------------
      IF (mode_flag) THEN
         infnum=modemax-modemin+1
         DO in=1,infnum
            label=labl+in-1
            edge_flag=.FALSE.
            IF (singfld_flag) THEN
               CALL ipout_singfld(in,xwpmn,dist,poloout,toroout,label)
            ENDIF
            IF (pmodb_flag) THEN
               CALL ipout_pmodb(in,xwpmn,poloout,toroout,label)
            ENDIF
            IF (eqbrzphi_flag .OR. brzphi_flag .OR. xrzphi_flag .OR. 
     $           vbrzphi_flag .OR. vpbrzphi_flag .OR. vvbrzphi_flag) 
     $           THEN
               IF (nrzeq_flag) THEN
                  nr=mr
                  nz=mz
               ENDIF
               IF (.NOT. eqoff_flag) CALL ipeq_rzpgrid(nr,nz)
               ALLOCATE(ipiv(mpert),
     $              invmats(mpert,mpert),temp1(mpert,mpert))
               DO i=1,mpert
                  invmats(i,i)=1.0
               ENDDO
               temp1=TRANSPOSE(permeabmats(modelnum,:,:))
               CALL zgetrf(mpert,mpert,temp1,mpert,ipiv,info)
               CALL zgetrs('N',mpert,mpert,temp1,mpert,
     $              ipiv,invmats,mpert,info)
               invmats=TRANSPOSE(invmats)
               CALL idcon_build(in,xwpmn)
               CALL ipeq_alloc
               CALL ipeq_contra(psilim)
               CALL ipeq_xptobn(psilim,xwpmn,bnomn)
               CALL ipeq_dealloc
               CALL ipeq_weight(psilim,bnomn,mfac,mpert,1)
               brrmn = MATMUL(invmats,bnomn)
               CALL ipeq_weight(psilim,bnomn,mfac,mpert,0)
               CALL ipeq_weight(psilim,brrmn,mfac,mpert,0)
               CALL ipout_xbrzphi(in,xwpmn,nr,nz,brrmn,bnomn,label)
               DEALLOCATE(ipiv,invmats,temp1)
            ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
            IF (singcurs_flag) THEN
               CALL ipdiag_singcurs(in,xwpmn,msing,resol,
     $              smallwidth,label)
            ENDIF
            IF (xbcontra_flag) THEN
               CALL ipdiag_xbcontra(in,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnormal_flag) THEN
               CALL ipdiag_xbnormal(in,xwpmn,poloout,toroout,label)
            ENDIF
            IF (xbnovc_flag) THEN
               CALL ipdiag_xbnovc(in,xwpmn,label)
            ENDIF
            IF (xbnobo_flag) THEN
               CALL ipdiag_xbnobo(in,xwpmn,d3_flag,label)
            ENDIF
            IF (xbnorm_flag) THEN
               CALL ipdiag_xbnorm(in,xwpmn,poloout,toroout,label)
            ENDIF
            IF (pmodbst_flag) THEN
               CALL ipdiag_pmodbst(in,xwpmn,label)
            ENDIF
            IF (pmodbrz_flag) THEN
               CALL ipdiag_pmodbrz(in,xwpmn,label)
            ENDIF
            IF (rzphibx_flag) THEN
               CALL ipdiag_rzphibx(in,xwpmn,label)
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     diagnose without a given error field.
c-----------------------------------------------------------------------
      IF (radvar_flag) THEN
         CALL ipdiag_radvar
      ENDIF
      IF (eigen_flag) THEN
         CALL ipdiag_eigen
      ENDIF
      IF (magpot_flag) THEN
         CALL ipdiag_magpot
      ENDIF
      IF (energy_flag) THEN
         CALL ipdiag_energy
      ENDIF
      IF (respmat_flag) THEN
         CALL ipdiag_respmat
      ENDIF
      IF (arbsurf_flag) THEN
         CALL ipdiag_arbsurf(majr,minr)
      ENDIF         
      IF (angles_flag) THEN
         CALL ipdiag_angles
      ENDIF
      IF (surfmode_flag) THEN
         CALL ipdiag_surfmode(lowmode,highmode,poloout,toroout)
      ENDIF
      IF (rzpgrid_flag) THEN
         IF (.NOT. eqoff_flag) CALL ipdiag_rzpgrid(nr,nz)
      ENDIF
      IF (m3d_flag) THEN
         normpsi=1.0
         fp=1e-3
         
         DO in=m3low,m3high
            label=in-m3low
            fxmn=0
            fxmn(in-mlow+1)=fp*normpsi
            CALL ipeq_cotoha(psilim,fxmn,mfac,mpert,0,0)
            fxmn=-twopi*ifac*chi1*(mfac-nn*qlim)*fxmn
            CALL ipeq_weight(psilim,fxmn,mfac,mpert,0)
            edge_flag=.FALSE.
            CALL ipout_errfld(infile,errtype,left,scale,fxmn,
     $           poloout,toroout,wegt,resp,bnomn,xwpmn,label)
            edge_flag=.TRUE.
            CALL ipout_singfld(0,xwpmn,dist,poloout,toroout,label)
         ENDDO
      ENDIF
      IF (cas3d_flag) THEN
c         fp=1e-4
c         DO in=m3low,m3high
c            label=in-m3low
c            fxmn=0
c            fxmn(in-mlow+1)=fp
c            edge_flag=.FALSE.
c            CALL ipout_errfld(infile,errtype,left,scale,fxmn,
c     $           poloout,toroout,resp,bnomn,xwpmn,label)
c            edge_flag=.TRUE.
c            CALL ipout_singfld(0,xwpmn,dist,1,1,label)
c            CALL ipdiag_xbcontra(0,xwpmn,1,1,label)
c            CALL ipdiag_xbcontra(0,xwpmn,4,1,label)
c            CALL ipdiag_xbnormal(0,xwpmn,1,1,label)
c            CALL ipdiag_xbnormal(0,xwpmn,4,1,label)
c            CALL ipdiag_xbnobo(0,xwpmn,d3_flag,label)
c            CALL ipdiag_radvar
c         ENDDO
         fp = -1e-2
         DO in=m3low,m3high
            label=in-m3low
            fxmn=0
            fxmn(in-mlow+1)=fp
            CALL ipeq_cotoha(psilim,fxmn,mfac,mpert,4,1)
            CALL ipeq_xptobn(psilim,fxmn,brrmn)
            edge_flag=.FALSE.
            CALL ipout_errfld(infile,errtype,left,scale,brrmn,
     $           poloout,toroout,wegt,resp,bnomn,xwpmn,label)
            edge_flag=.TRUE.
            CALL ipout_singfld(0,xwpmn,dist,1,1,label)
            CALL ipdiag_xbcontra(0,xwpmn,1,1,label)
            CALL ipdiag_xbcontra(0,xwpmn,4,1,label)
            CALL ipdiag_xbnormal(0,xwpmn,1,1,label)
            CALL ipdiag_xbnormal(0,xwpmn,4,1,label)
            CALL ipdiag_xbnobo(0,xwpmn,d3_flag,label)
            CALL ipdiag_radvar
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     only for test.
c-----------------------------------------------------------------------
      IF (test_flag) THEN
         fxmn=0
         fxmn(10-mlow+1)=1e-4
         bnomn=fxmn
         CALL ipeq_cotoha(psilim,bnomn,mfac,mpert,0,0)
         CALL ipeq_hatoco(psilim,bnomn,mfac,mpert,0,0)
         CALL ascii_open(out_unit,"iptest_coordtrans_n"//
     $        sn//".out","UNKNOWN")
         WRITE(out_unit,*)
     $        "IPTEST_COORDTRANS: cotoha to hatoco, co(0,0)"
         WRITE(out_unit,'(2(1x,a16))')"binmn","boutmn"
         DO in=1,mpert
            WRITE(out_unit,'(2(1x,e16.8))')
     $           REAL(fxmn(in)),REAL(bnomn(in))
         ENDDO
         CALL ascii_close(out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL ipec_dealloc
      CALL ipec_stop("Normal termination.")
      END PROGRAM ipec_main
