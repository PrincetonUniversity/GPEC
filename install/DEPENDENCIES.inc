#------------------------------------------------------------------------------
# DEPENDENCIES.inc
#
# Purpose:
#   Specifies dependencies between object files and modules.
#
# What it does:
#   - Ensures files and modules are built in the correct order
#   - Lists which object files depend on others (e.g. ../equil/grid.o depends on ../equil/local.o)
#
# How to use:
#   - If you add a new object file or module that depends on others, add its dependencies here.
#   - Keeping this file up-to-date helps avoid build errors from missing dependencies.
#------------------------------------------------------------------------------

# Equil
../equil/local.o: ../equil/io.o
../equil/utils.o: ../equil/local.o
../equil/jacobi.o: ../equil/local.o
../equil/grid.o: ../equil/local.o ../equil/spline.o
../equil/spline.o: ../equil/local.o ../equil/utils.o
../equil/bicube.o: ../equil/spline.o
../equil/cspline.o: ../equil/local.o ../equil/spline.o ../equil/utils.o
../equil/fft.o: ../equil/local.o
../equil/fspline.o: ../equil/local.o ../equil/spline.o ../equil/cspline.o ../equil/utils.o ../equil/fft.o
../equil/global.o: ../equil/bicube.o ../equil/cspline.o
../equil/inverse.o: ../equil/global.o ../equil/utils.o ../equil/grid.o
../equil/direct.o: ../equil/global.o
../equil/chease.o: ../equil/inverse.o
../equil/read_eq.o: ../equil/inverse.o ../equil/direct.o
../equil/lar.o: ../equil/inverse.o
../equil/gsec.o: ../equil/inverse.o
../equil/sol.o: ../equil/direct.o
../equil/equil.o: ../equil/read_eq.o ../equil/lar.o ../equil/gsec.o ../equil/sol.o
../equil/equil_out.o: ../equil/global.o

# LSODE
../lsode/daxpy.o: ../equil/local.o
../lsode/dcopy.o: ../equil/local.o
../lsode/dewset.o: ../equil/local.o
../lsode/dgbsl.o: ../equil/local.o
../lsode/dgesl.o: ../equil/local.o
../lsode/dprepj.o: ../equil/local.o
../lsode/dsolsy.o: ../equil/local.o
../lsode/dstode.o: ../equil/local.o
../lsode/dvnorm.o: ../equil/local.o
../lsode/iumach.o: ../equil/local.o
../lsode/lsode.o: ../equil/local.o
../lsode/xsetf.o: ../equil/local.o
../lsode/dcfode.o: ../equil/local.o
../lsode/ddot.o: ../equil/local.o
../lsode/dgbfa.o: ../equil/local.o
../lsode/dgefa.o: ../equil/local.o
../lsode/dintdy.o: ../equil/local.o
../lsode/dscal.o: ../equil/local.o
../lsode/dsrcom.o: ../equil/local.o
../lsode/dumach.o: ../equil/local.o
../lsode/idamax.o: ../equil/local.o
../lsode/ixsav.o: ../equil/local.o
../lsode/xerrwd.o: ../equil/local.o
../lsode/xsetun.o: ../equil/local.o

# ZLANGE
../zlange/zlange.o: ../equil/local.o

# ZVODE
../zvode/zvode.o: ../equil/local.o

# Orbit
../orbit/orbit_mod.o: ../equil/equil.o ../equil/equil_out.o
../orbit/restart.o: ../orbit/orbit_mod.o
../orbit/field.o: ../orbit/orbit_mod.o
../orbit/inner.o: ../orbit/field.o ../orbit/restart.o
../orbit/outer.o: ../orbit/field.o ../orbit/restart.o
../orbit/inner_gc.o: ../orbit/field.o ../orbit/restart.o
../orbit/outer_gc.o: ../orbit/field.o ../orbit/restart.o
../orbit/orbit.o: ../orbit/inner.o ../orbit/outer.o ../orbit/inner_gc.o ../orbit/outer_gc.o

# Vacuum
../vacuum/vacuum_global.o: ../equil/local.o
../vacuum/vacuum_diag.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_io.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_ma.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_math.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_sprk2.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_ut.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_vac.o: ../vacuum/vacuum_global.o
../vacuum/debug.o: ../vacuum/vacuum_global.o
../vacuum/vacuum_penn.o: ../equil/local.o

# PENTRC
../pentrc/lsode1.o: ../equil/local.o ../lsode/lsode.o
../pentrc/lsode2.o: ../equil/local.o ../lsode/lsode.o
../pentrc/params.o: ../equil/local.o ../pentrc/version.inc
../pentrc/utilities.o: ../equil/local.o ../pentrc/params.o
../pentrc/special.o: ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o
../pentrc/dcon_interface.o: ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../equil/spline.o ../equil/fspline.o ../equil/cspline.o ../equil/bicube.o
../pentrc/inputs.o:  ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../pentrc/special.o ../pentrc/dcon_interface.o
../pentrc/energy.o:  ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../pentrc/lsode2.o ../pentrc/dcon_interface.o
../pentrc/pitch.o:   ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../pentrc/lsode1.o ../pentrc/energy.o ../pentrc/special.o
../pentrc/torque.o:  ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../pentrc/lsode1.o ../pentrc/special.o ../pentrc/dcon_interface.o ../pentrc/inputs.o ../pentrc/energy.o ../pentrc/pitch.o
../pentrc/diagnostics.o: ../equil/local.o ../pentrc/params.o ../pentrc/special.o ../pentrc/energy.o ../pentrc/pitch.o
../pentrc/pentrc_interface.o: ../equil/local.o ../pentrc/params.o ../pentrc/utilities.o ../pentrc/special.o ../pentrc/inputs.o ../pentrc/energy.o ../pentrc/pitch.o ../pentrc/torque.o ../pentrc/diagnostics.o
../pentrc/pentrc.o:  ../equil/local.o ../pentrc/params.o  ../pentrc/pentrc_interface.o

# DCON
../dcon/dcon_mod.o: ../equil/local.o ../equil/spline.o ../equil/global.o \
	../equil/equil.o ../equil/equil_out.o ../pentrc/dcon_interface.o ../pentrc/pentrc_interface.o
../dcon/mercier.o: ../equil/local.o ../dcon/dcon_mod.o
../dcon/bal.o: ../equil/local.o ../dcon/dcon_mod.o
../dcon/fourfit.o: ../equil/local.o ../equil/fspline.o ../equil/global.o \
	../pentrc/pentrc_interface.o ../pentrc/utilities.o ../pentrc/torque.o ../pentrc/inputs.o \
	../pentrc/energy.o ../pentrc/pitch.o ../pentrc/dcon_interface.o \
	../dcon/dcon_mod.o
../dcon/sing.o: ../equil/local.o ../dcon/fourfit.o
../dcon/debug.o: ../equil/local.o
../dcon/ode_output.o: ../equil/local.o ../dcon/dcon_mod.o ../dcon/sing.o ../dcon/debug.o
../dcon/ode.o: ../equil/local.o ../dcon/ode_output.o ../dcon/free.o
../dcon/dcon_netcdf.o: ../equil/local.o ../dcon/dcon_mod.o ../dcon/sing.o
../dcon/free.o: ../equil/local.o ../equil/global.o ../dcon/dcon_netcdf.o ../dcon/ode_output.o
../dcon/resist.o: ../equil/local.o ../dcon/dcon_mod.o
../dcon/dcon.o: ../equil/local.o ../equil/equil.o ../equil/equil_out.o \
	../pentrc/inputs.o ../pentrc/pentrc_interface.o \
	../dcon/bal.o ../dcon/mercier.o ../dcon/ode.o ../dcon/free.o ../dcon/resist.o

# RDCON
../rdcon/dcon_mod.o: ../equil/spline.o ../equil/global.o \
	../equil/equil.o ../equil/equil_out.o ../rdcon/version.inc
../rdcon/toolbox.o: ../rdcon/dcon_mod.o
../rdcon/mercier.o: ../rdcon/dcon_mod.o
../rdcon/bal.o: ../rdcon/dcon_mod.o
../rdcon/fourfit.o: ../equil/fspline.o ../rdcon/dcon_mod.o ../rdcon/toolbox.o
../rdcon/sing1.o: ../rdcon/fourfit.o
../rdcon/sing.o: ../rdcon/sing1.o
../rdcon/ode_output.o: ../rdcon/sing.o ../rdcon/debug.o
../rdcon/ode.o: ../rdcon/ode_output.o
../rdcon/rdcon_netcdf.o: ../rdcon/dcon_mod.o
../rdcon/free.o: ../rdcon/ode.o ../rdcon/rdcon_netcdf.o
../rdcon/resist.o: ../rdcon/dcon_mod.o ../rdcon/toolbox.o
../rdcon/gal.o: ../equil/jacobi.o ../rdcon/sing.o ../rdcon/free.o
../rdcon/lib_interface.o: ../equil/local.o ../equil/global.o ../equil/utils.o ../equil/inverse.o
../rdcon/dcon.o: ../rdcon/bal.o ../rdcon/mercier.o ../rdcon/ode.o ../rdcon/resist.o ../rdcon/gal.o ../rdcon/lib_interface.o
../rdcon/debug.o: ../equil/local.o

# STRIDE
../stride/dcon_mod.o: ../equil/spline.o ../equil/global.o \
	../equil/equil.o ../equil/equil_out.o ../stride/version.inc
../stride/mercier.o: ../stride/dcon_mod.o
../stride/bal.o: ../stride/dcon_mod.o
../stride/fourfit.o: ../equil/fspline.o ../stride/dcon_mod.o
../stride/sing.o: ../stride/fourfit.o
../stride/zvode1.o: ../zvode/zvode.o ../zlange/zlange.o
../stride/riccati.o: ../stride/sing.o
../stride/ode.o: ../stride/sing.o ../stride/debug.o ../stride/zvode1.o ../stride/riccati.o
../stride/stride_netcdf.o: ../stride/dcon_mod.o
../stride/free.o: ../stride/ode.o ../stride/stride_netcdf.o
../stride/stride.o: ../stride/bal.o ../stride/free.o ../stride/mercier.o
../stride/debug.o: ../equil/local.o

# Match
../match/debug.o: ../equil/local.o
../match/gamma.o: ../equil/local.o
../match/lar.o: ../equil/inverse.o
../match/qsort.o: ../equil/local.o
../match/deltar.o: ../match/gamma.o
../match/root.o: ../equil/local.o
../match/zroot.o: ../equil/local.o
../match/match_mod.o: ../match/deltar.o ../match/root.o ../equil/utils.o
../match/ideal.o: ../match/match_mod.o ../equil/cspline.o
../match/matrix.o: ../match/match_mod.o ../match/debug.o ../match/qsort.o ../equil/local.o
../match/resist.o: ../match/matrix.o
../match/match.o: ../match/ideal.o ../match/resist.o

# Rmatch
../rmatch/msing.o: ../equil/local.o
../rmatch/gamma.o: ../equil/local.o
../rmatch/deltar.o: ../rmatch/gamma.o ../equil/spline.o
../rmatch/inps.o:  ../equil/local.o ../equil/utils.o
../rmatch/inpso.o: ../rmatch/inps.o ../rmatch/deltar.o
../rmatch/deltac.o: ../equil/jacobi.o ../rmatch/inpso.o
../rmatch/inner.o: ../rmatch/deltac.o ../equil/local.o
../rmatch/match.o: ../rmatch/deltar.o ../rmatch/deltac.o ../rmatch/inner.o ../rmatch/msing.o

# Sum
../sum/local.o: ../sum/io.o
../sum/utils.o: ../sum/local.o
../sum/sum.o: ../equil/utils.o ../equil/local.o

# Multi
../multi/local.o: ../multi/io.o
../multi/multi.o: ../equil/local.o

# SLAYER
../slayer/sglobal.o: ../equil/local.o
../slayer/params.o : ../slayer/sglobal.o
../slayer/delta.o : ../slayer/sglobal.o
../slayer/slayer_netcdf.o : ../slayer/sglobal.o
../slayer/layerinputs.o : ../slayer/sglobal.o ../slayer/params.o ../slayer/slayer_netcdf.o
../slayer/gslayer.o : ../slayer/sglobal.o ../slayer/params.o ../slayer/delta.o ../slayer/layerinputs.o ../slayer/slayer_netcdf.o
../slayer/slayer.o : ../slayer/sglobal.o ../slayer/params.o ../slayer/delta.o

# COIL
../coil/coil.o: ../equil/local.o ../equil/bicube.o ../coil/version.inc
../coil/field.o: ../coil/coil.o

# GPEC
../gpec/gpglobal.o: ../equil/fspline.o ../equil/bicube.o
../gpec/ismath.o: ../gpec/gpglobal.o  ../equil/local.o
../gpec/idcon.o: ../gpec/gpglobal.o ../gpec/ismath.o
../gpec/rdcon.o: ../gpec/gpglobal.o ../gpec/ismath.o ../gpec/idcon.o
../gpec/gpeq.o: ../gpec/idcon.o ../gpec/ismath.o
../gpec/gpresp.o: ../gpec/gpeq.o
../gpec/gpvacuum.o : ../gpec/gpeq.o
../gpec/gpdiag.o: ../gpec/gpresp.o ../gpec/gpvacuum.o
../gpec/gpout.o: ../gpec/gpdiag.o ../gpec/gpresp.o ../gpec/gpvacuum.o ../coil/field.o ../coil/coil.o \
	../equil/equil.o ../equil/equil_out.o \
	../pentrc/pentrc_interface.o ../pentrc/inputs.o \
	../slayer/gslayer.o \
	../gpec/version.inc
../gpec/gpec.o: ../gpec/gpdiag.o ../gpec/gpout.o ../gpec/rdcon.o

# XDraw
# no dependencies

# Full project dependencies
$(ALL_OBJS): | v neededdeps