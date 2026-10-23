"""Execute the production field-line failure formatter without an MHD solve."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


class DirectErrorMessage(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is required')
    def test_step_limit_preserves_complete_diagnostic(self):
        root=Path(__file__).resolve().parents[2]
        source=(root/'equil/direct.f').read_text()
        routine=source.split('      SUBROUTINE direct_fl_int(',1)[1].split(
            '      END SUBROUTINE direct_fl_int',1)[0]
        declaration=re.search(r'(?m)^      CHARACTER[^\n]*:: message',routine)[0]
        format_line=re.search(r'(?m)^ 40   FORMAT[^\n]*',routine)[0]
        write=re.search(r'(?m)^         WRITE\(message,40\)[^\n]*\n     \$[^\n]*',routine)[0]
        # Compile the actual declaration/FORMAT/WRITE from direct_fl_int. The
        # small driver reaches the error path deterministically with native I/O.
        program='\n'.join(['      PROGRAM diagnostic','      IMPLICIT NONE',declaration,
            '      INTEGER :: enstep=10000, ipsi=512','      REAL(8) :: eta=1.0D0',
            format_line,write,"      WRITE(*,'(a)') TRIM(message)",'      END',''])
        scratch=Path(os.environ.get('TMPDIR',root/'build-test'));scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as directory:
            path=Path(directory);(path/'diagnostic.f').write_text(program)
            subprocess.run(['gfortran','-ffixed-line-length-none',str(path/'diagnostic.f'),
                            '-o',str(path/'diagnostic')],check=True,capture_output=True)
            result=subprocess.run([str(path/'diagnostic')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('direct_int: istep = enstep = 10000',result.stdout)
            self.assertIn(' at eta =  1.000E+00, ipsi = 512',result.stdout)


if __name__=='__main__':unittest.main()
