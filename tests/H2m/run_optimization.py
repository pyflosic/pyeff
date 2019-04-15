from pyeff_optimizer import * 

# structure 
p_cfg = './structs/H2m.cfg'

# variables 
emax=0.0000001
fmax=0.005
steps = 10000
fix_nuc = 'True'

# optimization 
print 'Optimization'
calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)
calc.run()
# show final structure 
calc.view()
# write final structure as xyz file 
calc.write_xyz()

# single point 
print 'Single Point Calculation'
calc = pyeff(p_cfg=p_cfg)
calc.initialize()
print calc.show_all()
