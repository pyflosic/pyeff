from pyeff_optimizer import * 

# structure 
p_cfg = './structs/BH3.cfg'

# variables 
emax=0.0000001
fmax=0.005
steps = 10000
fix_nuc = None 

# optimization 
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()
