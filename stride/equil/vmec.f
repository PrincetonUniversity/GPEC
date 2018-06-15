      dcon_file = "dcon_" // TRIM(input_extension) // ".txt"
      OPEN (unit=51,FILE=dcon_file,FORM='FORMATTED',iostat=istat)
      IF (istat .ne. 0) STOP 'Error writing dcon output file'
 
      IF (mnmax .ne. mpol) STOP 'THIS IS NOT AXISYMMETRIC!'
 
      WRITE (51, *) ns                             !Number of flux surfaces
      WRITE (51, *) mpol                           !Number of poloidal modes, m=[0:mpol-1]
      WRITE (51, *) rmncc(1:ns,0,0:mpol1)          !r = sum [rmnc * cos (mu)], full mesh
      WRITE (51, *) zmnsc(1:ns,0,0:mpol1)          !z = sum [zmns * sin (mu)], full mesh
      WRITE (51, *) lmnsc(1:ns,0,0:mpol1)          !lam = sum[lmns * sin(mu)], half mesh
!     NOTE: u + lam give a straight magnetic field line
      IF (lasym) THEN
         WRITE (51, *) rmnsc(1:ns,0,0:mpol1)       !r = r+sum [rmns * sin (mu)], full mesh
         WRITE (51, *) zmncc(1:ns,0,0:mpol1)       !z = z+sum [zmnc * cos (mu)], full mesh
         WRITE (51, *) lmncc(1:ns,0,0:mpol1)       !lam = lam+sum[lmnc * cos(mu)], half mesh
      END IF
      WRITE (51, *) chi(1:ns)                      !pol flux, full mesh
      WRITE (51, *) fpsi(2:ns)                     !R*BT, half mesh
      WRITE (51, *) presf(1:ns)                    !pressure, full mesh (MKS)
      WRITE (51, *) 1/iotaf(1:ns)                  !q, full mesh
 
      CLOSE (unit=51)
