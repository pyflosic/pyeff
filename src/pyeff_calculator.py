from __future__ import division

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, \
    ReadError

from pyeff_energy_force import * 
from ase.atom import chemical_symbols
from ase.atoms import Atoms 

#from ase.units import Ha,Bohr

# global variables 
# numerical zero 
# note: do not make it to small !
zero = 0.0000001 # to avoid divide by zero error 

class pyeff(FileIOCalculator):
	implemented_properties = ['energy', 'forces']
    	default_parameters = {'p_cfg': None,'fix_nuc': None}

        def __init__(self, restart=None, ignore_bad_restart_file=False,
                label=os.curdir, atoms=None, **kwargs):
                """Constructor pyeff calculator.
                   by Sebastian Schwalbe                
                """

                FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                          label, atoms, **kwargs)
                valid_args = ('p_cfg','fix_nuc')
                # set any additional keyword arguments
                #print self.parameters.items()
		for arg, val in self.parameters.items():
                        if arg in valid_args:
                                setattr(self, arg, val)
                        else:
                                raise RuntimeError('unknown keyword arg "%s" : not in %s'% (arg, valid_args))
	
	def initialize(self, atoms=None,properties=['energy'],system_changes=['positions']):
       		Calculator.calculate(self, atoms, properties, system_changes)
	
		# pyeff_energy_force calculation 
		[pyeff_calc,Z,R,r,s,delta] = read_cfg(self.p_cfg)
		types = pyeff_calc.types()
		self.types = types
		chemical_symbols = pyeff_calc.chemical_symbols()
		self.chemical_symbols = chemical_symbols
		positions = pyeff_calc.positions() 
		self.positions =  positions 
		sizes = pyeff_calc.sizes()
		self.sizes = sizes 
		spins = pyeff_calc.spins()
		self.spins = spins 
		Eke = E_ke(pyeff_calc,s)
		Enucnuc = E_nuc_nuc(pyeff_calc,Z,R)
		Enucelec = E_nuc_elec(pyeff_calc,Z,R,r,s)
		Eelecelec = E_elec_elec(pyeff_calc,r,s)
		EPauli = E_Pauli(pyeff_calc,delta,s,r)
		self.results['energy'] = pyeff_calc.total_energy()
		self.results['forces'] = np.array(pyeff_calc.forces())
		# atoms object 
		bohr = 0.5291772105638411
		# unit conversion in pyeff positons have bohr: bohr2angstroem
		atoms = Atoms(''.join(self.get_chemical_symbols()),positions=np.array(pyeff_calc.positions()))# *bohr)
		self.atoms = atoms 
		fr = pyeff_calc.rforces()
                self.fr = fr
		# 1d vectors for pyeff optimization 
		positions1d = pyeff_calc.positions1d()
                self.positions1d = positions1d
		self.positions = positions1d
		forces1d = pyeff_calc.forces1d(fix_nuc=self.fix_nuc)
		self.forces1d = forces1d
		self.forces = forces1d 
		# general access of the pyeff object 
		self.pyeff_calc = pyeff_calc


	def calculate(self, atoms, properties=['energy'],system_changes=['positions']):
		Calculator.calculate(self, atoms, properties, system_changes)

		atoms = self.get_atoms()
		atoms = self.atoms 
		self.write_input(atoms)
		# pyeff_energy_force calculation 
		[pyeff_calc,Z,R,r,s,delta] = read_cfg(p_cfg='pyeff.cfg')
                types = pyeff_calc.types()
                self.types = types
                chemical_symbols = pyeff_calc.chemical_symbols()
                self.chemical_symbols = chemical_symbols
                positions = pyeff_calc.positions()
                self.positions =  positions
                sizes = pyeff_calc.sizes()
                self.sizes = sizes
                spins = pyeff_calc.spins()
                self.spins = spins
                Eke = E_ke(pyeff_calc,s)
                Enucnuc = E_nuc_nuc(pyeff_calc,Z,R)
                Enucelec = E_nuc_elec(pyeff_calc,Z,R,r,s)
                Eelecelec = E_elec_elec(pyeff_calc,r,s)
                EPauli = E_Pauli(pyeff_calc,delta,s,r)
		self.results['energy'] = pyeff_calc.total_energy()
                self.results['forces'] = np.array(pyeff_calc.forces())
                # radial forces 
		fr = pyeff_calc.rforces()
                self.fr = fr
		# 1d vectors for pyeff optimization  
                positions1d = pyeff_calc.positions1d()
                self.positions = positions1d
                forces1d = pyeff_calc.forces1d(fix_nuc=self.fix_nuc)
                self.forces = forces1d 
		# general access of the pyeff object  
		self.pyeff_calc = pyeff_calc 

	def write_input(self,atoms):
		# writes a cfg input file 
		# input variables 
		positions = atoms.positions 
		self.pos_transform(atoms) 
		sizes = self.sizes
		# write new file  
		o = open('pyeff.cfg','w')
		o.write('@params\n')
		o.write('calc = minimize\n')
		o.write('@nuclei\n')
		for p in range(len(self.types)):
			if self.types[p] == 'nuclei':
				o.write(str(positions[p][0])+' '+str(positions[p][1])+' '+str(positions[p][2])+' '+str(self.chemical_symbols[p])+'\n')
		o.write('@electrons\n')
		for p in range(len(self.types)):
                        if self.types[p] == 'electron':
                                o.write(str(positions[p][0])+' '+str(positions[p][1])+' '+str(positions[p][2])+' '+str(self.spins[p])+' '+str(sizes[p])+'\n')
		o.close()

        def get_energy(self):
		# return the pyeff total energy 
		atoms = self.atoms
		self.check_state(atoms) 
		self.calculate(atoms)
                return self.results['energy']

	def get_force(self):
		# return the pyeff forces for all particles 
		# in trible pairs 
		atoms = self.atoms 
		self.calculate(atoms)
		return self.results['forces']
	
	def get_force1d(self):
		# return the pyeff forces as 1d vector 
		return self.forces1d
	
	def get_chemical_symbols(self):
		# return chemical symbols
		# 1st spin channel = X symbol  
		# 2nd spin channel = He symbol 
		symbols = []
		for s in range(len(self.chemical_symbols)):
			if self.types[s] == 'electron' and int(self.chemical_symbols[s]) == int(1): 
				symbols.append('X')
			elif self.types[s] == 'electron' and int(self.chemical_symbols[s]) == int(-1):
				symbols.append('He')
			else:
				symbols.append(chemical_symbols[int(self.chemical_symbols[s])])
		return symbols
	
	def get_positions(self):
		# returns ase like positions object 
		return self.positions

	def get_positions1d(self):
		# returns the position as 1d vector 
                return self.positions1d 

	def check_state(self, atoms):	
		# from the standard ase.calculators 
                system_changes = FileIOCalculator.check_state(self, atoms)
                return system_changes

        def set(self, **kwargs):
		# from the standard ase.calculators 
                changed_parameters = FileIOCalculator.set(self, **kwargs)
                if changed_parameters:
                        self.reset()
 
	def get_atoms(self):
		# returns a ase.atoms object 
		# needed for the pyeff optimization 
        	if self.atoms is None:
            		raise ValueError('Calculator has no atoms')
        	atoms = self.atoms.copy()
        	atoms.calc = self
        	return atoms

	def pos_transform(self,atoms):
		# transform 1d properties to the correspondin ase format 
		sizes = [] 
		positions = []
		fr = [] 
		forces = [] 
		idx = 0
		for p in range(len(self.types)):
			if self.types[p] == 'nuclei':
				positions.append([self.positions[idx+0],self.positions[idx+1],self.positions[idx+2]])
				forces.append([self.forces[idx+0],self.forces[idx+1],self.forces[idx+2]])
				sizes.append(0)
				fr.append(0) 
				idx = idx + 3 
			if self.types[p] == 'electron':
				positions.append([self.positions[idx+0],self.positions[idx+1],self.positions[idx+2]])
				forces.append([self.forces[idx+0],self.forces[idx+1],self.forces[idx+2]])
				sizes.append(np.exp(self.positions[idx+3]))
				fr.append(self.forces[idx+3]/(self.positions[idx+3]+zero))
				idx = idx + 4 
		self.sizes = sizes 
		self.positions = np.array(positions)
		self.forces = forces 
		self.fr = fr 
	
	def show_all(self):
		# print formated output (forces and energy) 
		self.pyeff_calc.show_all()

	def view(self):
		# use ase.gui for visualization 
		# the import is only needed here 
		from ase.visualize import view
		bohr = 0.5291772105638411
                # unit conversion from bohr to angstroem 
		pos = np.array(self.pyeff_calc.positions())*bohr 
		sym = self.get_chemical_symbols()
		atoms = Atoms(sym,pos)
		view(atoms)
	
	def write_xyz(self):
		# writes a xyz file 
		# the import is only needed here
		from ase.io import write 
		bohr = 0.5291772105638411
		# unit conversion from bohr to angstroem 
		pos = np.array(self.pyeff_calc.positions())*bohr 
                sym = self.get_chemical_symbols()
                atoms = Atoms(sym,pos)
		write('pyeff.xyz',atoms,'xyz')

        def write_frmorb(self):
                # writes a cfg input file 
                # input variables 
		positions = np.array(self.pyeff_calc.positions())
		spin1 = self.spins.count(1)
		spin2 = self.spins.count(-1)
                # write frmorb 
                o = open('FRMORB.pyeff','w')
                o.write(str(spin1)+' '+str(spin2)+'\n')
		for p in range(len(self.types)):
                        if self.types[p] == 'electron':
                                o.write(str(positions[p][0])+' '+str(positions[p][1])+' '+str(positions[p][2])+'\n')
                o.close()

	def set_positions(self,positions):
		# test for BE 
		self.positions = positions


 	def set_positions1d(self,positions1d):
                # test for BE 
                self.positions1d = positions1d

