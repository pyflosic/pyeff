'Author: S. Schwalbe'
'part of pyeff' 

import numpy as np 

class particle:
	'mother class'
        'particle (nuclei or electron)' 
	'with properties: typ, label, position (px,py,pz), spin, energy, force (fx,fy,fz,fr)'
        def __init__(self,typ,label,px,py,pz,pr,spin,energy,fx,fy,fz,fr):
                self.typ = typ
		self.label = label 
		# position
                self.px = px
                self.py = py
                self.pz = pz
		self.pr = pr 
		# spin 
		self.spin = spin 
		# energy 
                self.energy = energy
		# force components 
                self.fx = fx
		self.fy = fy
		self.fz = fz
		self.fr = fr
		
        def update(self,px,py,pz,pr,spin,energy,fx,fy,fz,fr):
		self.px = self.px + px
                self.py = self.py + py
                self.pz = self.pz + pz
                self.pr = self.pr + pr
                self.spin = self.spin + spin 
		self.energy = self.energy + energy
		self.fx = self.fx + fx
		self.fy = self.fy + fy
		self.fz = self.fz + fz		
                self.fr = self.fr + fr

	def show(self):
		print '%10s\t%10.0f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f' %(self.typ,float(self.label),float(self.energy),float(self.fx),float(self.fy),float(self.fz),float(self.fr))

class system(particle):
        'contains all nuclei and electrons of the system'
        def __init__(self,nuc_chemical_symbols,elec_chemical_symbols):
                self.nuc_chemical_symbols = nuc_chemical_symbols
                self.elec_chemical_symbols = elec_chemical_symbols 
		self.nuc = []
                self.elec = []
                for n in range(len(nuc_chemical_symbols)):
                        n_tmp = particle(typ='nuclei',
					label=nuc_chemical_symbols[n],
					px=0,
					py=0,
					pz=0,
					pr=0,
					spin=0,
					energy=0,
					fx=0,
					fy=0,
					fz=0,
					fr=0)
                        self.nuc.append(n_tmp)
                for e in range(len(elec_chemical_symbols)):
                        e_tmp = particle(typ='electron',
					label=elec_chemical_symbols[e],
					energy=0,
	                                px=0,
                                        py=0,
                                        pz=0,
					pr=0,
					spin=0,
					fx=0,
					fy=0,
					fz=0,
					fr=0)
                        self.elec.append(e_tmp)

	def total_energy(self):
		etot = 0
	        for n in self.nuc:
	                etot = etot + float(n.energy)
        	for e in self.elec:
                	etot = etot+ float(e.energy)
        	return etot
	
	def types(self):
                types = []
                for n in self.nuc:
                        types.append(n.typ)
                for e in self.elec:
                        types.append(e.typ)
                return types 


	def chemical_symbols(self):
	        chemical_symbols = []
                for n in self.nuc:
                        chemical_symbols.append(n.label)
                for e in self.elec:
                        chemical_symbols.append(e.label)
                return chemical_symbols

	def positions(self):
		# positions as triple for all particles (nuclei, electrons) 
                positions = []
                for n in self.nuc:
                        positions.append([n.px,n.py,n.pz])
                for e in self.elec:
                        positions.append([e.px,e.py,e.pz])
                return positions

	def forces(self):
		# forces as triple for all particles (nuclei, electrons) 
		forces = []
                for n in self.nuc:
			forces.append([n.fx,n.fy,n.fz])
                for e in self.elec:
			forces.append([e.fx,e.fy,e.fz])
		return forces

	def rforces(self):
		# radial forces as 1d vector 
		rforces = []
                for n in self.nuc:
                        rforces.append(n.fr)
                for e in self.elec:
                        rforces.append(e.fr)
                return rforces

	def forces1d(self,fix_nuc=None):
		# forces as 1d vector for optimization 
	 	forces1d = []
                for n in self.nuc:
			if fix_nuc != None:
				n.fx = 0	
				n.fy = 0 
				n.fz = 0 
                        forces1d.extend([-1*n.fx,-1*n.fy,-1*n.fz])
                for e in self.elec:
                        forces1d.extend([-1*e.fx,-1*e.fy,-1*e.fz,-1*(e.pr)*e.fr])
                return forces1d
	 

	def positions1d(self):
		# forces as 1d vector for optimization 
                positions1d = []
                for n in self.nuc:
                        positions1d.extend([n.px,n.py,n.pz])
                for e in self.elec:
                        positions1d.extend([e.px,e.py,e.pz,np.log(e.pr)])
                return positions1d

	def sizes(self):
		sizes = []
                for n in self.nuc:
                        sizes.append(n.pr)
                for e in self.elec:
                        sizes.append(e.pr)
                return sizes

	def spins(self):
		spins = []
		for n in self.nuc:
                        spins.append(n.spin)
                for e in self.elec:
                        spins.append(e.spin)
                return spins

	def show_all(self):
		print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('----------','----------','----------','----------','----------','----------','----------')
		print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('typ','label','energy','fx','fy','fz','fr')
		print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('----------','----------','----------','----------','----------','----------','----------')
		for n in self.nuc:
                	n.show()
		for e in self.elec:
                        e.show()
		print '\n'
                print '---------------'
                print ' Total Energy  '
                print '---------------'
		print '\n'
                print 'Etot =\t %10s' %(self.total_energy())
		print '\n'


if __name__ == "__main__":
	def main():
		# functionality test 
		nuc = ['C','H','H','H']
		elec = ['He','He','X','X']
		sys = system(nuc,elec)
		print sys.nuc[0].label
		print len(sys.elec)
		(sys.elec[0]).update(px=1,py=0,pz=0,pr=0,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
	main()
