from numpy import * 
import numpy as np
import scipy.optimize as opt
from pyeff_calculator_pbc import *
import sys

Ha = 27.211386024367243
Bohr = 0.529177210563841

class pyeff_pbc_optimize:
	"""Optimization for PyEFF PBC"""
	"""for pyeff based on scipy and ase optimiters"""
	"""author: S. Schwalbe"""
	def __init__(self,p_cfg,emax,fmax,steps,fix_nuc=None,method='BFGS',scale=1.0):
		#self.x = x 		# x		: a one dimensional vector containing all position variables 
		#self.f = f 		# f		: pyeff energy 
		#self.x0 = x0 		# x0	 	: pyeff starting positions 
		#self.fprime = fprime 	# fprime 	: pyeff forces 
		self.fmax = fmax 	# fmax		: max force criteria 
		self.steps = steps	# steps		: maximal number of iteration steps 
		self.p_cfg = p_cfg 	# p_cfg 	: path to cfg file 
		self.H0 = 70		# H0		: has something to do with how hessian is treated in scipy 
		# emax 					: energy tolerance criteria 
		self.etol = emax 
		# optimizer 
		self.method = method 
		# nucleus fixed or not 
		self.fix_nuc = None
                if fix_nuc != None:
                        self.fix_nuc = fix_nuc
		self.scale = scale 
		# pyeff calculation object 
                calc = pyeff_pbc(p_cfg=self.p_cfg,fix_nuc=self.fix_nuc,scale=self.scale)
                calc.initialize()
		# Simulation box dimensions 
		self.Lx = calc.Lx
		self.Ly = calc.Ly
		self.Lz = calc.Lz
		# starting varibales 
                self.chemical_symbols = calc.chemical_symbols
               	self.types = calc.types 
		self.positions = np.array(calc.positions)*Bohr
                self.spins = calc.spins
                self.sizes = calc.sizes
		# actual calculation object 
		self.calc = calc 
		# print the intial values 
		calc.show_all()	
		# collection of all calculation atoms objects 
		# collect initial atoms object 
                atoms = calc.get_atoms()
                # unit conversion: bohr2angstroem 
                pos = atoms.get_positions()
                pos = pos*Bohr
                atoms.set_positions(pos)
		self.traj = [atoms] 
		# collection of all energies 
		self.energies = [calc.get_energy()*Ha]  
		# logfiles 
		logfile = open('opt.log','a')
		self.logfile = logfile
	def write_input(self,x):
		def wrap(xi,Li):
                     if xi < 0: 
                         xi = xi + Li
                     if xi > Li:
                         xi = xi - Li
                     return xi	
                x = np.array(x)/Bohr
		x = x.tolist()
		# write a new cfg file based on the x vector 
		chemical_symbols = self.chemical_symbols 
		spins = self.spins 
                o = open('pyeff.cfg','w')
                o.write('@params\n')
                o.write('x_bound = 0.000000 %0.6f\n' % (self.Lx))
                o.write('y_bound = 0.000000 %0.6f\n' % (self.Ly))
                o.write('z_bound = 0.000000 %0.6f\n' % (self.Lz))
		o.write('calc = minimize\n')
                o.write('@nuclei\n')
                idx = 0
		for p in range(len(self.types)):
                        if self.types[p] == 'nuclei':
				dx = x[idx+0]
				dy = x[idx+1]
				dz = x[idx+2] 
				#dx = wrap(dx,self.Lx)
				#dy = wrap(dy,self.Ly)
				#dz = wrap(dz,self.Lz)
                                o.write('%0.5f %0.5f %0.5f %s\n' %(dx,dy,dz,chemical_symbols[p]))
				idx = idx + 3 
                o.write('@electrons\n')
                for p in range(len(self.types)):
                        if self.types[p] == 'electron':
				dx = x[idx+0]
				dy = x[idx+1]
				dz = x[idx+2]
				#dx = wrap(dx,self.Lx)
				#dy = wrap(dz,self.Ly)
				#dz = wrap(dz,self.Lz)
                                o.write('%0.5f %0.5f %0.5f %i %0.5f\n' %(dx,dy,dz,spins[p],np.exp(x[idx+3])))
				idx = idx +4 
                o.close()
		self.write_xyz()


	def x0(self):
		# starting guess for all position like variables in a 1d vector 
		x0 = self.calc.get_positions1d()
		x0 = np.array(x0)*Bohr
		x0 = x0.tolist()
		return x0

	def f(self,x):
		# objective function, pyeff total energy 
		# may need to rewrite the io part 
		self.write_input(x)
		calc = pyeff_pbc(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
                calc.initialize()
                # Check boundary
                # calc.check_boundary()
		f = calc.get_energy()*Ha
		self.calc = calc
		self.logfile.write('Energy = %0.9f\n' % f)
		self.logfile.flush()
		return f/self.H0 

	def fprime(self,x):
		# 1st derivative of objective function, note: note forces, sign! 
		# may need to rewrite the io part 
		#self.write_input(x)
		#calc = pyeff_pbc(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
		#calc.initialize()
		#fprime =  np.array(calc.get_force1d(),dtype=np.float64)
		fprime =  np.array(self.calc.get_force1d(),dtype=np.float64)
		fprime = fprime*Ha/Bohr
		#print('Force = %0.9f eV/Ang' % f)
		fmax = np.sqrt((fprime**2).max())
		self.logfile.write('Fmax = %0.9f\n' % fmax)
                self.logfile.flush()
		return fprime/self.H0 

	def callback(self, x):
		# this function is like a output function called after each scipy iteration step 
		# may need to rewrite the io part 
		#self.write_input(x)
                #calc = pyeff_pbc(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
                #calc.initialize()
		calc = self.calc 
		# Simulation box dimensions 
                self.Lx = calc.Lx
                self.Ly = calc.Ly
                self.Lz = calc.Lz
		calc.show_all()	
		# update the calculation object 
		self.calc = calc
		# collect atoms object 
		atoms = calc.get_atoms() 
		# unit conversion: bohr2angstroem 
		bohr = 0.5291772105638411
		pos = atoms.get_positions()
		pos = pos*bohr 
		atoms.set_positions(pos)
		# append all structures to a xyz trajectory 
		self.traj.append(atoms) 
		# append all total energy per iteration 
		etot = calc.get_energy()*Ha
		# declare a additional stopping criteria 
		# based on the energy difference 
		ediff = abs(self.energies[-1] - etot)
		print 'Ediff = ',ediff
		if ediff <= self.etol: 
			print 'Converged: achieved Etol'
			self.write_xyz()
			sys.exit() 
		else:	self.energies.append(calc.get_energy()*Ha)

	def view(self):
		self.calc.view()

	def write_xyz(self):
		from ase.io import write 
		write('pyeff_movie.xyz',self.traj,'extxyz',write_results=False)

	def write_frmorb(self):
		self.calc.write_frmorb()

	def run(self):
		# the actual optimization function 
		if self.method == 'BFGS':
			output = opt.fmin_bfgs(self.f,
				self.x0(),
				fprime=self.fprime,
				# args=(),
				gtol=self.fmax * 0.1,  # Should never be reached
				norm=np.inf,
				#epsilon=1.4901161193847656e-08,
				maxiter=self.steps,
				full_output=1,
				disp=0,
				# retall=0,
				callback=self.callback
				)
		if self.method == 'CG':
			output = opt.fmin_cg(self.f,
                                self.x0(),
                                fprime=self.fprime,
                                # args=(),
                                gtol=self.fmax * 0.1,  # Should never be reached
                                norm=np.inf,
                                #epsilon=1.4901161193847656e-08,
                                maxiter=self.steps,
                                full_output=1,
                                disp=0,
                                # retall=0,
                                callback=self.callback
                                )

		#print(output)
		warnflag = output[-1]
		if warnflag == 2:
			print('Desired error not necessarily achieved (due to precision loss)')

if __name__ == "__main__":

	import time 
	t0 = time.clock()
	#p_cfg = 'Li_222_eff_opt.cfg'
	p_cfg = 'li_solid_222.cfg'
	#p_cfg = 'pyeff.cfg'
	emax=0.00001
	fmax=0.002
	steps = 10000
	fix_nuc = 'True'
	calc = pyeff_pbc_optimize(p_cfg,emax,fmax,steps,fix_nuc,method='BFGS')	
	calc.run()
	#calc.view()
	calc.write_xyz()
	t = time.clock() - t0
	print('Timing: %0.5f s' % t)

