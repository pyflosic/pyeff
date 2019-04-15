from pyeff_energy_force import * 
from pyeff_optimizer import * 
from pyeff_io import * 

from matplotlib.pyplot import * 
import numpy as np 
from ase.io import write 

from nrlmol_io import *

def set_x(types,chemical_symbols,spins,x,dist_nuc,dist_value,er):
        # moves nuclei and electrons positions according the dist value 
	# Note: the guess needs to be symmetric 
	# types			...		nuclei/electrons
	# chemical_symbols	...		chemical symbols 
	# spins			...		eff spins 
	# x 			...		1d positions arrary (e.g. x0 = calc.positions1d())
	# dist_nuc		...		array of nuclei which should be moved 
	# dist_value		...		distance value 
	# er			...		electron radius 
	idx = 0
        for p in range(len(types)):
                if types[p] == 'nuclei':
                        if int(chemical_symbols[p]) == int(dist_nuc[0]) or int(chemical_symbols[p]) == int(dist_nuc[1]):
                                x[idx+0] = x[idx+0] + np.sign(x[idx+0])*dist_value
                        print 'nuclei '+str(chemical_symbols[p])+' '+str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])
                        idx = idx+3
                if types[p] == 'electron':
                        x[idx+0] = x[idx+0] + np.sign(x[idx+0])*dist_value
			if er != None:
                                x[idx+3] = er
                        print 'electron '+str(chemical_symbols[p])+' '+str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])+' '+str(np.exp(x[idx+3]))
                        idx = idx+4
        return x


# initial structure 
p_cfg = './structs/LiH.cfg'

# read cfg 
[calc,Z,R,r,s,delta] = read_cfg(p_cfg)
# intital values 
types = calc.types()
chemical_symbols = calc.chemical_symbols()
spins = calc.spins()
x0 = calc.positions1d()

# optimization values 
# for a smooth curve we need tight parameters 
emax=0.00000000000001
fmax=0.001
steps = 10000
# fix the nuclei 
fix_nuc = True

# binding energy 
# relative displacement 
dist = [0.1,0.11,0.125,0.15,0.16,0.17,0.175,0.2,0.225,0.25,0.275,0.3]
# final energies 
E = []
# bonding length 
b_len = []
# collection of the optimized structures at b_len 
traj = []

# calculation 
for d in dist:
 	# new cfg 
	p_new = 'new.cfg'
	# new x 
	x = set_x(types,chemical_symbols,spins,x=x0,dist_nuc=['3','1'],dist_value=d,er=3)
	# write new x 
	write_cfg(p_new,types,chemical_symbols,spins,x)       
	# unit conversion 
        bohr = 0.5291772105638411
        print x[0]
	print x[3]
	print abs(x[0]-x[3])
	b_len.append(abs(x[0]-x[3])*bohr)
        write_cfg(p_new,types,chemical_symbols,spins,x)
        # calculate the new configuration 
        opt = pyeffBFGS(p_new,emax,fmax,steps,fix_nuc)
        opt.run()
        # save values 
        # unit conversion         
        atoms = opt.traj[-1]
        pos = atoms.get_positions()
        #pos = pos*bohr
        atoms.set_positions(pos)
        traj.append(atoms)
        opt.calc.show_all()
        E.append(opt.energies[-1])


# save binding energy curve geometries as trajectory 
write('trajectory.xyz',traj,'extxyz',write_results=False)

# plot binding energy curve 
fig1 = figure(1)
ax = subplot(111)
plot(b_len,E,'-o')
xlabel(r'$d_{Li-H}$ [$\AA$]')
ylabel(r'Energy [Hartree]')
show()

# save the calculated energies 
np.savetxt('BE.dat',np.transpose([b_len,E]))
