from pyeff_optimizer import * 

# structure 
p_cfg = './structs/Li.cfg'

# variables 
emax=0.0000000000000000001
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
