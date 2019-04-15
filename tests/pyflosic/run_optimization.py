from pyflosic_io import *
from pyeff_optimizer import pyeffBFGS 

# structure 
#xyz2cfg(p_xyz='./structs/He10.xyz',er=0.7) 
#xyz2cfg(p_xyz='./structs/B.xyz',er=0.7) 
p_cfg = 'pyflosic.cfg'

# variables 
emax=0.00001
fmax=0.001
steps = 100000
fix_nuc = None

# optimization 
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()
calc.write_frmorb()
