"""Prove and emit the real leaves; Fourier indexing stays in pmodb_gauge_mod."""
from pathlib import Path
import fortsym
import fortsym.sympy as fs
s, Tr, Ti = fs.symbols('s Tr Ti')
weight = 1/(2*fs.pi*s)
wr, wi = -weight*Ti, weight*Tr
assert (2*fs.pi*s*weight-1).is_zero is True
assert (2*fs.pi*fs.I*s*(wr+fs.I*wi)+Tr+fs.I*Ti).is_zero is True
a, b, bt, dr, di, vr, vi, ur, ui = fs.symbols('a b bt dr di vr vi ur ui')
source = {'tr': vr+a*dr/b, 'ti': vi+a*di/b}
target = {'out_r': dr+ur*bt, 'out_i': di+ui*bt}
assert (source['tr']+fs.I*source['ti']-(vr+fs.I*vi)-a*(dr+fs.I*di)/b).is_zero is True
assert (target['out_r']+fs.I*target['out_i']-(dr+fs.I*di)-(ur+fs.I*ui)*bt).is_zero is True
kernels = [
    fortsym.emit_kernel('pmodb_gauge_inverse', ['s','Tr','Ti'], {'wr':wr,'wi':wi}, generator='gen/pmodb_gauge/generate.py'),
    fortsym.emit_kernel('pmodb_gauge_source', ['a','b','dr','di','vr','vi'], source, generator='gen/pmodb_gauge/generate.py'),
    fortsym.emit_kernel('pmodb_gauge_apply', ['dr','di','ur','ui','bt'], target, generator='gen/pmodb_gauge/generate.py')]
root = Path(__file__).resolve().parents[2]
(root/'gpec/pmodb_gauge_kernel.f90').write_text('module pmodb_gauge_kernel\ncontains\n'+'\n'.join(kernels)+'\nend module pmodb_gauge_kernel\n')
print('PROVED and emitted: complex source, inverse field-line coefficient, scalar correction')
