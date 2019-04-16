from pyeff_optimizer_pbc import * 
import time

t0 = time.clock()

p_cfg = 'li_solid_222.cfg'
emax=1e-10
fmax=0.002
steps = 10000
fix_nuc = 'True'
method = 'BFGS' 
scale = 1.0025 

calc = pyeff_pbc_optimize(p_cfg,emax,fmax,steps,fix_nuc,method=method,scale=scale)
calc.run()
calc.write_xyz()
t = time.clock() - t0

print('Timing: %0.5f s' % t)

