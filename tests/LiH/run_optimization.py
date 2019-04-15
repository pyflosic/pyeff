from pyeff_optimizer import * 

# structure 
p_cfg = './structs/LiH.cfg'

# variables 
emax=0.00000001
fmax=0.0001
steps = 100000
fix_nuc = None

# optimization 
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()
