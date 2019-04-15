from pyeff_optimizer import * 

# structure 
p_cfg = './structs/CH4.cfg'
# variables 
emax=0.000000001
fmax=0.002
steps = 10000
fix_nuc = None
# optimization 
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()
