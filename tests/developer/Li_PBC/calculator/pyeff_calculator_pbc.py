from __future__ import division

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, \
    ReadError

from pyeff_energy_force import E_ke
from pyeff_ewald_system import * 
from pyeff_ewald_energy_force import * 
from ase.atom import chemical_symbols
from ase.atoms import Atoms 
import os 

#from ase.units import Ha,Bohr

# global variables 
# numerical zero 
# note: do not make it to small !
zero = 0.0000001 # to avoid divide by zero error 

class pyeff_pbc(FileIOCalculator):
	implemented_properties = ['energy', 'forces']
    	default_parameters = {'p_cfg': None,'fix_nuc': None,'scale':1.0}

        def __init__(self, restart=None, ignore_bad_restart_file=False,
                label=os.curdir, atoms=None, **kwargs):
                """Constructor pyeff calculator.
                   by Sebastian Schwalbe                
                """

                FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                          label, atoms, **kwargs)
                valid_args = ('p_cfg','fix_nuc','scale')
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
		[Lx,Ly,Lz,pyeff_calc,Z,R,r,s,delta] = read_cfg_pbc(self.p_cfg,scale=self.scale)
		# Get box lengths 
		self.Lx = Lx 
		self.Ly = Ly 
		self.Lz = Lz 
		# Get molecular data 
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
		# Set up periodic ewald calculator 
	        ewald_calc = ewald_system(Lx=Lx,Ly=Ly,Lz=Lz,system=pyeff_calc)
	        # Check boundary
                # ewald_calc.check_boundary()
                # Initialize grids 
		rgrid = ewald_calc.make_r_space()
	        kgrid = ewald_calc.make_k_space()
	        # Clear molecular data 
	        for n in range(len(ewald_calc.system.nuc)):
			ewald_calc.system.nuc[n].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*ewald_calc.system.nuc[n].energy,fx=-1*ewald_calc.system.nuc[n].fx,fy=-1*ewald_calc.system.nuc[n].fy,fz=-1*ewald_calc.system.nuc[n].fz,fr=-1*ewald_calc.system.nuc[n].fr)
        	for e in range(len(ewald_calc.system.elec)):
                	ewald_calc.system.elec[e].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*ewald_calc.system.elec[e].energy,fx=-1*ewald_calc.system.elec[e].fx,fy=-1*ewald_calc.system.elec[e].fy,fz=-1*ewald_calc.system.elec[e].fz,fr=-1*ewald_calc.system.elec[e].fr)

		# Periodic energies and forces 
		# E_kin 
		Eke = E_ke(ewald_calc.system,s)
		# E_PauliPeriodic
		EPauli =E_PauliPeriodic(ewald_calc)
		# E_self 
		Eself = E_self(ewald_calc)
		# E_kspace
		Ekspace = E_kspace(ewald_calc)
		# E_rspace
		Erspace = E_rspace(ewald_calc)
		# E_uniform
		Euniform = E_uniform(ewald_calc)
		# E_total 
		ewald_calc.add_ewald_energy_force()

		self.results['energy'] = ewald_calc.system.total_energy() 
		self.results['forces'] = np.array(ewald_calc.system.forces())
		# atoms object 
		bohr = 0.5291772105638411
		# unit conversion in pyeff positons have bohr: bohr2angstroem
		atoms = Atoms(''.join(self.get_chemical_symbols()),positions=np.array(ewald_calc.system.positions()))# *bohr)
		self.atoms = atoms 
		fr = ewald_calc.system.rforces()
                self.fr = fr
		# 1d vectors for pyeff optimization 
		positions1d = ewald_calc.system.positions1d()
                self.positions1d = positions1d
		self.positions = positions1d
		forces1d = ewald_calc.system.forces1d(fix_nuc=self.fix_nuc)
		self.forces1d = forces1d
		self.forces = forces1d 
		# general access of the pyeff object 
		self.ewald_calc = ewald_calc

	def calculate(self, atoms, properties=['energy'],system_changes=['positions']):
		Calculator.calculate(self, atoms, properties, system_changes)

		atoms = self.get_atoms()
		atoms = self.atoms 
		self.write_input(atoms)
		# pyeff_energy_force calculation 
                [Lx,Ly,Lz,pyeff_calc,Z,R,r,s,delta] = read_cfg_pbc('pyeff.cfg')
                # Get box lengths 
                self.Lx = Lx
                self.Ly = Ly
                self.Lz = Lz
                # Get molecular data 
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
                # Set up periodic ewald calculator 
                ewald_calc = ewald_system(Lx=Lx,Ly=Ly,Lz=Lz,system=pyeff_calc)
                # Check boundary
                # ewald_calc.check_boundary() 
                # Initialize grids 
                rgrid = ewald_calc.make_r_space()
                kgrid = ewald_calc.make_k_space()
                # Clear molecular data 
                for n in range(len(ewald_calc.system.nuc)):
                        ewald_calc.system.nuc[n].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*ewald_calc.system.nuc[n].energy,fx=-1*ewald_calc.system.nuc[n].fx,fy=-1*ewald_calc.system.nuc[n].fy,fz=-1*ewald_calc.system.nuc[n].fz,fr=-1*ewald_calc.system.nuc[n].fr)
                for e in range(len(ewald_calc.system.elec)):
                        ewald_calc.system.elec[e].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*ewald_calc.system.elec[e].energy,fx=-1*ewald_calc.system.elec[e].fx,fy=-1*ewald_calc.system.elec[e].fy,fz=-1*ewald_calc.system.elec[e].fz,fr=-1*ewald_calc.system.elec[e].fr)

                # Periodic energies and forces 
                # E_kin 
                Eke = E_ke(ewald_calc.system,s)
                # E_PauliPeriodic
                EPauli =E_PauliPeriodic(ewald_calc)
                # E_self 
                Eself = E_self(ewald_calc)
                # E_kspace
                Ekspace = E_kspace(ewald_calc)
                # E_rspace
                Erspace = E_rspace(ewald_calc)
                # E_uniform
                Euniform = E_uniform(ewald_calc)
                # E_total 
                ewald_calc.add_ewald_energy_force()


		self.results['energy'] = ewald_calc.system.total_energy()
                self.results['forces'] = np.array(ewald_calc.system.forces())
                # atoms object 
                bohr = 0.5291772105638411
                # unit conversion in pyeff positons have bohr: bohr2angstroem
                atoms = Atoms(''.join(self.get_chemical_symbols()),positions=np.array(ewald_calc.system.positions()))# *bohr)
                self.atoms = atoms
                fr = ewald_calc.system.rforces()
                self.fr = fr
                # 1d vectors for pyeff optimization 
                positions1d = ewald_calc.system.positions1d()
                self.positions1d = positions1d
                self.positions = positions1d
                forces1d = ewald_calc.system.forces1d(fix_nuc=self.fix_nuc)
                self.forces1d = forces1d
                self.forces = forces1d
                # general access of the pyeff object 
                self.ewald_calc = ewald_calc

	def write_input(self,atoms):
		def wrap(xi,Li):
                     if xi < 0:
                         xi = xi + Li
                     if xi > Li:
                         xi = xi - Li
                     return xi

		# writes a cfg input file 
		# input variables 
		#positions = atoms.positions 
		self.pos_transform(atoms) 
                positions = self.positions	
                sizes = self.sizes
		# write new file  
		o = open('pyeff.cfg','w')
		o.write('@params\n')
		# x_bound = 0.000000 8.352000
		# y_bound = 0.000000 8.352000
		# z_bound = 0.000000 8.352000
		o.write('x_bound = 0.000000 %0.6f\n' % self.Lx)
		o.write('y_bound = 0.000000 %0.6f\n' % self.Ly)
		o.write('z_bound = 0.000000 %0.6f\n' % self.Lz)
		o.write('calc = minimize\n')
		o.write('@nuclei\n')
		for p in range(len(self.types)):
			dx = positions[p][0]
			dy = positions[p][1]
			dz = positions[p][2]
			#dx = wrap(dx,self.Lx)
			#dy = wrap(dy,self.Ly)
			#dz = wrap(dz,self.Lz)
			if self.types[p] == 'nuclei':
				o.write('%0.5f %0.5f %0.5f %s\n' %(dx,dy,dz,self.chemical_symbols[p]))
				#o.write(str(wrap(positions[p][0],self.Lx))+' '+str(wrap(positions[p][1],self.Ly))+' '+str(wrap(positions[p][2],self.Lz))+' '+str(self.chemical_symbols[p])+'\n')
		o.write('@electrons\n')
		for p in range(len(self.types)):
                        if self.types[p] == 'electron':
				dx = positions[p][0]
		                dy = positions[p][1]
                	        dz = positions[p][2]
                        	#dx = wrap(dx,self.Lx)
                        	#dy = wrap(dy,self.Ly)
                        	#dz = wrap(dz,self.Lz)
				o.write('%0.5f %0.5f %0.5f %i %0.5f\n' %(dx,dy,dz,self.spins[p],sizes[p]))
                                #o.write(str(wrap(positions[p][0],self.Lx))+' '+str(wrap(positions[p][1],self.Ly))+' '+str(wrap(positions[p][2],self.Lz))+' '+str(self.spins[p])+' '+str(sizes[p])+'\n')
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
				positions.append([self.positions1d[idx+0],self.positions1d[idx+1],self.positions1d[idx+2]])
				forces.append([self.forces1d[idx+0],self.forces1d[idx+1],self.forces1d[idx+2]])
				sizes.append(0)
				fr.append(0) 
				idx = idx + 3 
			if self.types[p] == 'electron':
				positions.append([self.positions1d[idx+0],self.positions1d[idx+1],self.positions1d[idx+2]])
				forces.append([self.forces1d[idx+0],self.forces1d[idx+1],self.forces1d[idx+2]])
				sizes.append(np.exp(self.positions1d[idx+3]))
				fr.append(self.forces1d[idx+3]/(self.positions1d[idx+3]+zero))
				idx = idx + 4 
		self.sizes = sizes 
		self.positions = np.array(positions)
		self.forces = forces 
		self.fr = fr 
	
	def show_all(self):
		# print formated output (forces and energy) 
		self.ewald_calc.system.show_all()

	def view(self):
		# use ase.gui for visualization 
		# the import is only needed here 
		from ase.visualize import view
		bohr = 0.5291772105638411
                # unit conversion from bohr to angstroem 
		pos = np.array(self.ewald_calc.system.positions())*bohr 
		sym = self.get_chemical_symbols()
		atoms = Atoms(sym,pos)
		view(atoms)
	
	def write_xyz(self):
		# writes a xyz file 
		# the import is only needed here
		from ase.io import write 
		bohr = 0.5291772105638411
		# unit conversion from bohr to angstroem 
		pos = np.array(self.ewald_calc.system.positions())*bohr 
                sym = self.get_chemical_symbols()
                atoms = Atoms(sym,pos)
		write('pyeff.xyz',atoms,'xyz')

        def write_frmorb(self):
                # writes a cfg input file 
                # input variables 
		positions = np.array(self.ewald_calc.system.positions())
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

        def check_boundary(self):
                self.ewald_calc.check_boundary()
		# 1d vectors for pyeff optimization 
                positions1d = self.ewald_calc.system.positions1d()
                self.positions1d = positions1d
                self.positions = positions1d
                
if __name__ == "__main__":
	import time 
	t0 = time.clock()
	# single point 
	print 'Single Point Calculation'
	p_cfg = 'li_solid_222.cfg'
	#p_cfg = 'Li_222_eff_opt.cfg'
        #p_cfg = 'eff.cfg'
	calc = pyeff_pbc(p_cfg=p_cfg,scale=1.01)
	calc.initialize()
	calc.show_all()
	print(calc.get_energy())
	t = time.clock() - t0
	print('Timing: %0.5f s' % t)

