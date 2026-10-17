"""Native supplied-rate D/T projection: analytic profile and invalid-input oracles."""
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from pentrc_collision_profile import ROOT, routine


class SuppliedRateTests(unittest.TestCase):
    def test_native_supplied_rates_and_fast_validation(self):
        source = (ROOT/'pentrc/inputs.f90').read_text()
        declaration = source[source.index('    integer :: kin_zi'):source.index('    type(cspline_type)')]
        parts = ''.join(routine(source, name) for name in (
            'read_kin', 'kin_eval', 'kin_eval_external', 'update_collision_values',
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
''' + declaration + '\ncontains\n' + parts
        generated = ROOT/'pentrc/supplied_rates.inc'
        if generated.exists():
            inputs += generated.read_text()
        inputs += '\nend module inputs\n'
        driver = '''program supplied_profile
use params, only: r8
use inputs
use spline_mod, only: spline_dealloc
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
implicit none
integer :: i,unit,ix,icase,mi,mimp
real(r8) :: psi,nD,nT,td,tt,omega,rows(9),f(9),f1(9)
character(512) :: mode,file
call get_command_argument(1,mode)
file='profile.kin'
do icase=1,2
  open(newunit=unit,file=trim(file),status='replace')
  write(unit,'(a)') 'psi ni ne ti te OmegaE logLambda nu_i nu_e'
  do i=0,24
    psi=i/24._r8
    nD=2.e19_r8*(1+.1_r8*psi)
    nT=1.e19_r8*(1+.2_r8*psi)
    td=1200._r8*(1+.3_r8*psi)
    tt=2100._r8*(1+.4_r8*psi)
    omega=0._r8
    if(icase==2) omega=-3500._r8
    rows=[psi,nD,nD+nT,td,900._r8*(1+.25_r8*psi),omega, &
          17.5_r8,700._r8*exp(2._r8*psi),1500._r8*exp(-1.5_r8*psi)]
    if(icase==2) rows(2:4:2)=[nT,tt]
    if(i==10) then
      select case(trim(mode))
      case('nan-rate'); rows(8)=ieee_value(0._r8,ieee_quiet_nan)
      case('inf-rate'); rows(9)=ieee_value(0._r8,ieee_positive_inf)
      case('nan-log'); rows(7)=ieee_value(0._r8,ieee_quiet_nan)
      case('inf-density'); rows(2)=ieee_value(0._r8,ieee_positive_inf)
      case('negative-rate'); rows(8)=-1._r8
      case('zero-rate'); rows(9)=0._r8
      end select
    end if
    if(trim(mode)=='columns') then
      write(unit,'(8es25.16e3)') rows(:8)
    else
      write(unit,'(9es25.16e3)') rows
    end if
  end do
  close(unit)
  mi=2; mimp=3
  if(icase==2) then
    mi=3; mimp=2
  end if
  call read_kin(file,1,1,mi,mimp,1._r8,1._r8,1._r8,1._r8,.true.,'rates-v1')
  ix=0
  do i=0,104
    psi=i/104._r8
    call kin_eval_external(psi,ix,f,f1)
    call kin_eval(psi,1)
    if(any(abs(f-kin%f)>1.e-12_r8*max(abs(f),abs(kin%f),tiny(1._r8)))) &
        error stop 'query values disagree'
    if(any(abs(f1-kin%f1)>1.e-12_r8*max(abs(f1),abs(kin%f1),tiny(1._r8)))) &
        error stop 'query derivatives disagree'
    write(*,'(i4,1x,13es25.16e3)') icase,psi,f,f1(6:8)
  end do
  call spline_dealloc(kin)
  call spline_dealloc(kin_collision)
  call spline_dealloc(kin_supplied)
end do
end program supplied_profile
'''
        nf=os.environ.get('NF_CONFIG','nf-config')
        ff=shlex.split(os.environ.get('NETCDF_FFLAGS') or subprocess.check_output([nf,'--fflags'],text=True))
        libs=shlex.split(os.environ.get('NETCDF_FLIBS') or subprocess.check_output([nf,'--flibs'],text=True))
        linker=shlex.split(os.environ.get('GPEC_TEST_LDFLAGS','-Wl,--gc-sections'))
        mathlibs=shlex.split(os.environ.get('GPEC_MATH_LIBS','-llapack -lblas'))
        files=['equil/io.f','equil/local.f','equil/utils.f','equil/spline.f',
               'pentrc/params.f90','pentrc/utilities.f90']
        with tempfile.TemporaryDirectory(prefix='.pentrc-rates-',dir=ROOT) as tmp:
            work=Path(tmp)
            (work/'inputs.f90').write_text(inputs)
            (work/'driver.f90').write_text(driver)
            (work/'version.inc').write_text("character(128), parameter :: version='supplied-rate-test'\n")
            for optimization in (['-O0','-fcheck=all'],['-O3','-ffast-math']):
                command=[os.environ.get('FC','gfortran'),*optimization,'-std=legacy',
                         '-ffree-line-length-none','-ffixed-line-length-none',
                         '-fallow-argument-mismatch','-ffunction-sections','-fdata-sections',
                         *ff,'-I'+str(work),*[str(ROOT/f) for f in files],
                         str(work/'inputs.f90'),str(work/'driver.f90'),*linker,*libs,*mathlibs,
                         '-o',str(work/'check')]
                result=subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=90)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                run=subprocess.run([str(work/'check')],cwd=work,text=True,capture_output=True,timeout=20)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                rows=0; maxerr=0.
                for line in run.stdout.splitlines():
                    tokens=line.split()
                    if len(tokens)!=14:
                        continue
                    case=int(tokens[0]); psi=float(tokens[1]); f=list(map(float,tokens[2:11]))
                    df=list(map(float,tokens[11:14]))
                    ni=(2.e19*(1+.1*psi) if case==1 else 1.e19*(1+.2*psi))
                    ti=(1200*(1+.3*psi) if case==1 else 2100*(1+.4*psi))
                    nu_i=700*math.exp(2*psi); nu_e=1500*math.exp(-1.5*psi)
                    expected=[ni,3.e19+4.e18*psi,ti*1.602e-19,900*(1+.25*psi)*1.602e-19,
                              0. if case==1 else -3500.,17.5,nu_i,nu_e,1.]
                    if expected[4]==0.:
                        self.assertEqual(f[4],0.,'zero physical rotation changed')
                    else:
                        self.assertLess(abs(f[4]/expected[4]-1.),1.e-14)
                    error=max(abs(a/b-1.) for i,(a,b) in enumerate(zip(f,expected)) if i!=4)
                    maxerr=max(maxerr,error,abs(df[1]/(2*nu_i)-1.),abs(df[2]/(-1.5*nu_e)-1.))
                    self.assertLess(abs(df[0]),1.e-10)
                    rows+=1
                self.assertEqual(rows,210,run.stdout[-1600:])
                self.assertLess(maxerr,2.e-6)
                print('supplied D/T/e profile',optimization,'rows=',rows,'maxerr=',maxerr)
                for mode in ('nan-rate','inf-rate','nan-log','inf-density','negative-rate','zero-rate','columns'):
                    bad=subprocess.run([str(work/'check'),mode],cwd=work,text=True,capture_output=True,timeout=20)
                    self.assertNotEqual(bad.returncode,0,optimization[0]+': accepted '+mode)
                    self.assertIn('PENTRC:',bad.stdout+bad.stderr)


if __name__=='__main__': unittest.main()
