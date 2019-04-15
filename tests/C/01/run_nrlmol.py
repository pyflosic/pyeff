from pyeff_optimizer import * 
from nrlmol_io import * 

# structure 
nrlmol2cfg(p_symbol='./structs/SYMBOL',p_frmorb='./structs/FRMORB',er=1.3)
p_cfg = './nrlmol.cfg'

# variables 
emax=0.0000001
fmax=0.0001
steps = 1000000
fix_nuc = 'True'

# optimization 
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()
calc.write_frmorb()
