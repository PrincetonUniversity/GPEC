"""Require the native guard's specific diagnostic, rather than any failure."""
import subprocess
import sys

executable, mode = sys.argv[1:]
assert mode in ("chebyshev", "divzero")
result = subprocess.run([executable, mode], capture_output=True, text=True, timeout=10)
output = result.stdout + result.stderr
expected = f"disable {mode}_flag for raw masked exports"
if result.returncode == 0 or expected not in output:
    print(output)
    raise SystemExit(f"Expected native {mode} rejection was not observed")
print(f"PASS native {mode} rejection with specific diagnostic")
