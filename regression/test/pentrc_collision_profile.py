"""Native read_kin and torque query must preserve positive physical collision rates.

Compile the production profile reader, its query routines and the exact query
calls used by torque.F90, with the native spline/table modules.  A positive
linear-temperature profile has an independent analytic Coulomb-rate oracle.

Use the native make check-pentrc-collisions target to select platform linker
flags. Direct invocation accepts GPEC_TEST_LDFLAGS (use -Wl,-dead_strip on Darwin).
"""
import math
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get('GPEC_ROOT') or Path(__file__).resolve().parents[2])


def routine(source, name):
    m = re.search(r'^    subroutine ' + name + r'\(.*?^    end subroutine ' + name + r'\b',
                  source, re.M | re.S | re.I)
    return m.group() + '\n' if m else ''


class CollisionProfileTests(unittest.TestCase):
    def test_native_offgrid_rates_match_positive_coulomb_oracle(self):
        source = (ROOT/'pentrc/inputs.f90').read_text()
        torque = (ROOT/'pentrc/torque.F90').read_text()
        call_external = re.search(r'call (?:kin_eval_external\(|spline_eval_external\(kin,)[^\n]+', torque).group()
        call_diagnostic = re.search(r'call (?:kin_eval\(x,1\)|spline_eval\(kin,x,1\))', torque).group()
        decl = source[source.index('    integer :: kin_zi'):source.index('    type(cspline_type)')]
        parts = ''.join(routine(source, name) for name in
                        ('read_kin', 'kin_eval', 'kin_eval_external', 'update_collision_values',
                         'kinetic_collisions', 'validate_kinetic_values', 'supplied_collision_values'))
        inputs = '''module inputs
use, intrinsic :: iso_fortran_env, only: int64
use params, only: r8,mp,me,e,twopi
use utilities, only: readtable,get_free_file_unit
use spline_mod, only: spline_type,spline_alloc,spline_fit, &
    spline_eval,spline_eval_external,spline_write1
implicit none
logical :: verbose=.false.
real(r8) :: chi1=1.0_r8
''' + decl + '\ncontains\n' + parts + \
            (ROOT/'pentrc/supplied_rates.inc').read_text() + '\nend module inputs\n'
        driver = '''program collision_profile
use params, only: r8
use inputs
use spline_mod, only: spline_eval,spline_eval_external,spline_dealloc
implicit none
integer :: i,unit,ix,icase
real(r8) :: psi,x,nfac,tfac,kin_f(9),kin_f1(9)
character(512) :: file,mode
call get_command_argument(1,mode)
file='profile.kin'
open(newunit=unit,file=trim(file),status='replace')
write(unit,'(a)') 'psi ni ne ti te omegae'
do i=0,199
  psi=0.99485_r8*i/199
  write(unit,'(6es25.16e3)') psi,1.e19_r8,1.e19_r8, &
       100._r8+72000._r8*(1-psi),1.2_r8*(100._r8+72000._r8*(1-psi)),-3000._r8
end do
close(unit)
do icase=1,3
  nfac=1._r8
  tfac=1._r8
  if(icase==2) then
    nfac=0.5_r8
    tfac=2._r8
  end if
  if(icase==3) file='CI_PROFILE_FILE'
  call read_kin(file,1,6,2,12,nfac,tfac,1._r8,1._r8,.true.)
  if(icase==3.and.trim(mode)=='invalid-query') call kin_eval(1._r8,1)
  ix=0
  do i=0,104
    psi=i/100._r8
    if(i==101) psi=0.984411_r8
    if(i==102) psi=0.985957_r8
    if(i==103) psi=0.9756_r8
    if(i==104) psi=0.989_r8
    if(icase==3) then
      psi=0.997995992_r8*i/100
      if(i==101) psi=.25_r8
      if(i==102) psi=.5_r8
      if(i==103) psi=.75_r8
      if(i==104) psi=.99_r8
    end if
    QUERY_EXTERNAL
    call spline_eval(kin,psi,1)
    if(abs(kin%f(9)-1._r8)>1.e-12_r8) error stop 'raw slayer Zeff changed'
    x=psi
    QUERY_DIAGNOSTIC
    write(*,'(2i4,1x,28es25.16e3)') icase,i,psi,kin_f,kin%f,kin_f1(7:8),kin_f1(1:5)
  end do
  call spline_dealloc(kin)
  call spline_dealloc(kin_collision)
end do
end program collision_profile
'''.replace('QUERY_EXTERNAL',call_external).replace('QUERY_DIAGNOSTIC',call_diagnostic).replace(
    'CI_PROFILE_FILE',str(ROOT/'docs/examples/regression_solovev_ideal_example/kinetic.txt'))
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', ROOT)) as tmp:
            work=Path(tmp)
            (work/'inputs.f90').write_text(inputs)
            (work/'driver.f90').write_text(driver)
            (work/'version.inc').write_text("character(128), parameter :: version='collision-profile-test'\n")
            files=['equil/io.f','equil/local.f','equil/utils.f','equil/spline.f',
                   'pentrc/params.f90','pentrc/utilities.f90']
            nf=os.environ.get('NF_CONFIG','nf-config')
            netcdf_flags=shlex.split(os.environ['NETCDF_FFLAGS']) if 'NETCDF_FFLAGS' in os.environ else \
                shlex.split(subprocess.check_output([nf,'--fflags'],text=True))
            netcdf_libs=shlex.split(os.environ['NETCDF_FLIBS']) if 'NETCDF_FLIBS' in os.environ else \
                shlex.split(subprocess.check_output([nf,'--flibs'],text=True))
            math_libs=shlex.split(os.environ.get('GPEC_MATH_LIBS','-llapack -lblas'))
            linker_flags=shlex.split(os.environ.get('GPEC_TEST_LDFLAGS','-Wl,--gc-sections'))
            command=[os.environ.get('FC','gfortran'),'-O0','-fcheck=all','-std=legacy',
                     '-ffree-line-length-none','-ffixed-line-length-none','-fallow-argument-mismatch',
                     '-ffunction-sections','-fdata-sections',*netcdf_flags,'-I'+str(work),
                     *[str(ROOT/f) for f in files],str(work/'inputs.f90'),str(work/'driver.f90'),
                     *linker_flags,*netcdf_libs,*math_libs,'-o',str(work/'check')]
            compiled=subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=90)
            self.assertEqual(compiled.returncode,0,compiled.stdout+compiled.stderr)
            run=subprocess.run([str(work/'check')],cwd=work,text=True,capture_output=True,timeout=30)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            invalid=subprocess.run([str(work/'check'),'invalid-query'],cwd=work,
                                   text=True,capture_output=True,timeout=30)
            self.assertNotEqual(invalid.returncode,0,'invalid consumed extrapolation accepted')
            self.assertIn('PENTRC: non-positive collision model coefficient',
                          invalid.stdout+invalid.stderr)
            # Debug-exported values use the same physical query, too.
            debug_rates=[]
            for line in (work/'pentrc_kinetics.out').read_text().splitlines():
                tokens=line.split()
                if len(tokens)!=11:
                    continue
                try:
                    fields=list(map(float,tokens))
                except ValueError:
                    continue
                debug_rates.extend(fields[8:10])
            self.assertTrue(debug_rates)
            self.assertGreater(min(debug_rates),0.,'debug export retains a negative rate')
            minnu=math.inf
            maxerr=0.
            violations=[]
            row_count=0
            maxprimitiveerr=0.
            min_ci_loglambda=math.inf
            # Native primitive values/gradients from unmodified develop
            # 3c0b8f3b, GNU Fortran binary64, using the exact CI input file.
            # This parent-code receipt verifies the existing two-stage
            # primitive interpolation, including its cold-edge error.
            primitive_receipt={
                0: (1.51795e21,1.4097600405650037e-17,3958.999995,
                    -1.0625634632084995e21,-4.2292892968672215e-18,-3959.0000556058844),
                100: (4.575143868829413e20,9.373007088720216e-20,7.933867835582145,
                      -1.0625654363734113e21,-4.65544415687748e-17,-3959.000011093797),
                101: (1.2523087477594865e21,1.2816000366244083e-17,2969.2499995564554,
                      -1.0625652711557528e21,-6.213818814970546e-18,-3959.001265133593),
                102: (9.866675000002848e20,1.0844307997650295e-17,1979.5000065100269,
                      -1.0625650002070866e21,-1.0010130403478095e-17,-3958.9997417756244),
                103: (7.210262496942402e20,7.419789684897973e-18,989.7500003888853,
                      -1.0625649979732148e21,-1.874473135647245e-17,-3959.0005334947086),
                104: (4.660106501937482e20,4.592052280328285e-19,39.590000328186484,
                      -1.0625651093240409e21,-4.487429506331625e-17,-3959.000030368476),
            }
            for line in run.stdout.splitlines():
                values=line.split()
                if len(values)!=28:
                    continue
                case,index=int(values[0]),int(values[1])
                row_count+=1
                psi=float(values[2]); f=list(map(float,values[3:12])); diagnostic=list(map(float,values[12:21]))
                nfac,tfac=(0.5,2.) if case==2 else (1.,1.)
                ni=ne=1.e19*nfac
                ti=(100.+72000.*(1.-psi))*tfac
                te=1.2*ti
                if case==3:
                    # The CI profile is n=n0*(1-.7*psi), n*T proportional to
                    # 1-psi. Its finite supplied endpoint has logLambda>5;
                    # extrapolating it to psi=1 violates weak coupling.
                    ni=ne=1.51795e21*(1-.7*psi)
                    # Rate coefficients belong to the consumed primitives,
                    # rather than the analytic profile before interpolation.
                    ti=f[2]/1.602e-19
                    te=f[3]/1.602e-19
                    if index in primitive_receipt:
                        pf=list(map(float,values[23:28]))
                        actual=(f[0],f[2],f[4],pf[0],pf[2],pf[4])
                        for a,b in zip(actual,primitive_receipt[index]):
                            maxprimitiveerr=max(maxprimitiveerr,abs(a/b-1.))
                            self.assertLess(abs(a/b-1.),2.e-10)
                        self.assertEqual(f[0],f[1])
                        self.assertEqual(f[2],f[3])
                # Preserve nfac/tfac semantics: collisions belong to the input
                # profiles; these optional multipliers alter drive profiles.
                base_ni,base_ne=ni/nfac,ne/nfac
                base_ti,base_te=ti/tfac,te/tfac
                lnl=17.3-0.5*math.log(base_ne/1.e20)+1.5*math.log(base_te/1.e3)
                if case==3:
                    min_ci_loglambda=min(min_ci_loglambda,f[5])
                    # The original reader's unsuffixed constants are binary32
                    # literals promoted to binary64; preserve that convention.
                    native=lambda v: struct.unpack('f',struct.pack('f',v))[0]
                    native_lnl=native(17.3)-.5*math.log(f[1]/native(1.e20)) \
                        +1.5*math.log(f[3]/native(1.602e-16))
                    self.assertLess(abs(f[5]/native_lnl-1.),2.e-13)
                # Pure deuterium: Zeff=zpitch=1; temperatures here are eV.
                nui=base_ni*lnl/(3.5e17*math.sqrt(2.)*(base_ti/1.e3)**1.5)
                nue=base_ne*lnl/(3.5e17*math.sqrt(9.1091e-31/1.672614e-27)*(base_te/1.e3)**1.5)
                self.assertEqual(f,diagnostic,'integrand and exported profile queries disagree')
                df=list(map(float,values[21:23]))
                self.assertLess(abs(f[0]/ni-1.),1.e-8 if case==3 else 1.e-9)
                if case!=3:
                    self.assertLess(abs(f[2]/(ti*1.602e-19)-1.),1.e-6)
                minnu=min(minnu,f[6],f[7])
                error=max(abs(f[6]/nui-1.),abs(f[7]/nue-1.))
                maxerr=max(maxerr,error)
                derivative=(1.5/lnl-1.5)*(-72000.*tfac)/ti
                if case==3:
                    pf=list(map(float,values[23:28]))
                    nd=pf[0]/f[0]
                    td=pf[2]/f[2]
                    derivative=nd+(-.5*nd+1.5*td)/lnl-1.5*td
                derivative_error=max(abs(df[0]/(nui*derivative)-1.),abs(df[1]/(nue*derivative)-1.))
                if derivative_error>2.e-6:
                    violations.append((case,index,psi,'derivative',derivative_error))
                if min(f[6:8])<=0 or error>1.e-6:
                    violations.append((case,index,psi,f[6],nui,error))
            self.assertTrue(math.isfinite(minnu),'no native rows read')
            self.assertEqual(row_count,315,'all three native profile cases must complete')
            self.assertGreater(min_ci_loglambda,5.)
            print('native collision profile: minnu=',minnu,'maxerr=',maxerr,'violations=',violations[:8])
            print('CI primitive preservation: maxerr=',maxprimitiveerr,
                  'minimum logLambda=',min_ci_loglambda)
            self.assertFalse(bool(violations),'native collision-rate positivity/Coulomb oracle failed: '+str(violations[:8]))


if __name__=='__main__': unittest.main()
