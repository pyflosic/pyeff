#   Copyright 2019 PyEFF developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
from numpy import * 
import numpy as np
import scipy.optimize as opt
from pyeff_calculator import *
import sys

class pyeffBFGS:
	"""Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)"""
	"""for pyeff based on scipy and ase optimiters"""
	"""author: S. Schwalbe"""
	def __init__(self,p_cfg,emax,fmax,steps,fix_nuc=None):
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
		# nucleus fixed or not 
		self.fix_nuc = None
                if fix_nuc != None:
                        self.fix_nuc = fix_nuc
		# pyeff calculation object 
                calc = pyeff(p_cfg=self.p_cfg,fix_nuc=self.fix_nuc)
                calc.initialize()
		# starting varibales 
                self.chemical_symbols = calc.chemical_symbols
               	self.types = calc.types 
		self.positions = np.array(calc.positions)
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
                bohr = 0.5291772105638411
                pos = atoms.get_positions()
                pos = pos*bohr
                atoms.set_positions(pos)
		self.traj = [atoms] 
		# collection of all energies 
		self.energies = [calc.get_energy()] 

	def write_input(self,x):
		# write a new cfg file based on the x vector 
		chemical_symbols = self.chemical_symbols 
		spins = self.spins 
                o = open('pyeff.cfg','w')
                o.write('@params\n')
                o.write('calc = minimize\n')
                o.write('@nuclei\n')
                idx = 0
		for p in range(len(self.types)):
                        if self.types[p] == 'nuclei':
                                o.write(str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])+' '+str(chemical_symbols[p])+'\n')
				idx = idx + 3 
                o.write('@electrons\n')
                for p in range(len(self.types)):
                        if self.types[p] == 'electron':
                                o.write(str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])+' '+str(spins[p])+' '+str(np.exp(x[idx+3]))+'\n')
				idx = idx +4 
                o.close()

	def x0(self):
		# starting guess for all position like variables in a 1d vector 
		x0 = self.calc.get_positions1d()
		return x0

	def f(self,x):
		# objective function, pyeff total energy 
		# may need to rewrite the io part 
		self.write_input(x)
		calc = pyeff(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
		calc.initialize()
		f = calc.get_energy()
		return f/self.H0 

	def fprime(self,x):
		# 1st derivative of objective function, note: note forces, sign! 
		# may need to rewrite the io part 
		self.write_input(x)
		calc = pyeff(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
		calc.initialize()
		fprime =  np.array(calc.get_force1d(),dtype=np.float64)
		return fprime/self.H0 

	def callback(self, x):
		# this function is like a output function called after each scipy iteration step 
		# may need to rewrite the io part 
		self.write_input(x)
                calc = pyeff(p_cfg='pyeff.cfg',fix_nuc=self.fix_nuc)
                calc.initialize()
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
		etot = calc.get_energy()
		# declare a additional stopping criteria 
		# based on the energy difference 
		ediff = abs(self.energies[-1] - etot)
		print 'Ediff = ',ediff
		if ediff <= self.etol: 
			print 'Converged: achieved Etol'
			self.write_xyz()
			sys.exit() 
		else:	self.energies.append(calc.get_energy())

	def view(self):
		self.calc.view()

	def write_xyz(self):
		from ase.io import write 
		write('pyeff_movie.xyz',self.traj,'extxyz',write_results=False)

	def write_frmorb(self):
		self.calc.write_frmorb()

	def run(self):
		# the actual optimization function 
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
		warnflag = output[-1]
		if warnflag == 2:
			print('Desired error not necessarily achieved (due to precision loss)')

if __name__ == "__main__":

	p_cfg = '../tests/CH4/structs/CH4.cfg'
	emax=0.001
	fmax=0.003
	steps = 1000000 
	fix_nuc = 'True'
	calc = pyeffBFGS(p_cfg,emax,fmax,steps,fix_nuc)	
	calc.run()
	calc.view()
	calc.write_xyz()
