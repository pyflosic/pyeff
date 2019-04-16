from pyeff_calculator_pbc import *
import time
t0 = time.clock()
# single point 
print 'Single Point Calculation'
p_cfg = 'pyeff.cfg'
calc = pyeff_pbc(p_cfg=p_cfg,scale=1.0)
calc.initialize()
calc.show_all()
print(calc.get_energy())
t = time.clock() - t0
print('Timing: %0.5f s' % t)

