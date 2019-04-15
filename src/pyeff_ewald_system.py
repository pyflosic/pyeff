'Author: S. Schwalbe'
'part of pyeff'

import numpy as np

class ewald_charge:
        'mother class'
        'ewald periodic boundary condition (pbc) class'
	'one ewald descriptor' 
        'with properties: position (px,py,pz),palpha (Gaussian exponent), charge q, energy, force (fx,fy,fz), falpha, force Gaussian exponent'
	def __init__(self,px,py,pz,palpha,q,energy,fx,fy,fz,falpha):
		# position
		self.px = px
		self.py = py
		self.pz = pz
		self.palpha = palpha
		self.q = q 
		# energy 
		self.energy = energy
		# force components 
		self.fx = fx
		self.fy = fy
		self.fz = fz
		self.falpha = falpha

        def update(self,px,py,pz,palpha,q,energy,fx,fy,fz,falpha):
                self.px = self.px + px
                self.py = self.py + py
                self.pz = self.pz + pz
                self.palpha = self.palpha + palpha
		self.q = self.q + q 
                self.energy = self.energy + energy
                self.fx = self.fx + fx
                self.fy = self.fy + fy
                self.fz = self.fz + fz
                self.falpha = self.falpha + falpha

        def show(self):
                print 'Ewald charge \t%10.9f\t%10.5f\t%10.5f\t%10.5f\t%10.5f' %(float(self.energy),float(self.fx),float(self.fy),float(self.fz),float(self.falpha))

# 1st a molecular calculation using system 
# 2nd give system to ewald_system 
# 3rd update system with ewald parameters 

class ewald_system(ewald_charge):
        'contains all ewald descriptors of the system'
	def __init__(self,Lx,Ly,Lz,system):
                # cell 
                # simulation box lengths 
                self.Lx = Lx
                self.Ly = Ly
                self.Lz = Lz
		self.system = system 
		
		#  InitializeEwald(double Lx, double Ly, double Lz, double max_alpha, double max_r, double max_k, double nuc_alpha)
		#  InitializeEwald(params.x_bound[1] - params.x_bound[0],
	        #    params.y_bound[1] - params.y_bound[0],
        	#    params.z_bound[1] - params.z_bound[0],
            	#    2.0 / (params.ewald_re_cutoff * params.ewald_re_cutoff), params.ewald_r_cutoff, params.ewald_k_cutoff, 1.0 / (params.ewald_nuc_r * params.ewald_nuc_r));

		# init default values 
                self.init_defaults()
		# max values 
		self.max_alpha = 2./(self.ewald_re_cutoff**2.)
		# max real space value 
		#self.max_r = self.ewald_r_cutoff
		self.max_r = self.get_ewald_r_cutoff()
		# max k space value 
		#self.max_k = self.ewald_k_cutoff
		self.max_k = self.get_ewald_k_cutoff()
		# collection of all ewald descriptors 
		self.ewald = []
		for n in range(len(system.nuc)):
			ewald_tmp = ewald_charge(px=system.nuc[n].px,
                                          py=system.nuc[n].py,
                                          pz=system.nuc[n].pz,
                        		  # density = |psi|^2
			                  palpha=(2.0/(self.ewald_nuc_r**2.)),
                                          q = system.nuc[n].label,
					  energy=0,
                                          fx=0,
                                          fy=0,
                                          fz=0,
                                          falpha=0)
                        self.ewald.append(ewald_tmp)
		
                for e in range(len(system.elec)):
        		ewald_tmp = ewald_charge(px=system.elec[e].px,
                                        py=system.elec[e].py,
                                        pz=system.elec[e].pz,
                                        # density = |psi|^2
                                        palpha=2.*(1.0/(system.elec[e].pr**2.)),
                                        q=-1.0,
                                        energy=0,
                                        fx=0,
                                        fy=0,
                                        fz=0,
                                        falpha=0)
                        self.ewald.append(ewald_tmp)       
		# Number of Ewald charges 
		self.numcharges = len(self.system.nuc)+len(self.system.elec)
			 
	
        def init_defaults(self):
		# -log10 desired ewald precision 
                self.ewald_log_precision = -6.0
             	# needed for max alpha 
		self.ewald_re_cutoff     = 3.54
		# automatic set r and k space cutoffs?  
                self.ewald_autoset       = 1
		# manual set r-space cutoff 
                self.ewald_r_cutoff      = 7.0
		# manual set k-space cutoff 
                self.ewald_k_cutoff      = 8.0
		# radius of nucleus in ewald 
                self.ewald_nuc_r         = 1e-10
		# maxium electron size ewald can describe 
		self.ewald_max_re        = 4.5
		# variables needed for 
		# get_ewald_r_cutoff
		# and  
		# get_ewald_k_cutoff 
	        self.precision = self.ewald_log_precision
                self.alpha_cutoff = 2./(self.ewald_re_cutoff**2.)
                self.min_alpha = 2./(self.ewald_max_re**2.)

	def get_ewald_r_cutoff(self):
		# calculate r space cutoff 
		# input: precision, alpha_cutoff, min_alpha 
		return np.sqrt((-1.*np.log(10.0) * self.precision + 3.) / (self.alpha_cutoff / (1. + self.alpha_cutoff / self.min_alpha)))
	
	def get_ewald_k_cutoff(self):
		# calculate k space cutoff 
		# input: precision, alpha_cutoff
		return np.sqrt(4. * self.alpha_cutoff * (-np.log(10.0) * self.precision + 5.));

	def get_volume(self):
		# volume of the system 
		# assuming a cubic box 
		return self.Lx*self.Ly*self.Lz

        def get_total_charge(self):
		# total charge of the system 
		total_charge= 0
		for c in range(self.numcharges):
			total_charge = total_charge+self.ewald[c].q
		return total_charge
	
	def add_ewald_energy_force(self):
		# add ewald energy and force to system(elec,nuc)
		idx = 0
		for n in range(len(self.system.nuc)):
			# system:nuc: self,px,py,pz,pr,spin,energy,fx,fy,fz,fr
                        self.system.nuc[n].update(px=0,
                                          py=0,
                                          pz=0,
                                          pr=0,
					  spin=0, 
                                          energy=self.ewald[idx].energy,
                                          fx=self.ewald[idx].fx,
                                          fy=self.ewald[idx].fy,
                                          fz=self.ewald[idx].fz,
                                          fr=0)
			idx = idx +1 
		for e in range(len(self.system.elec)):
			# system:elec: self,px,py,pz,pr,spin,energy,fx,fy,fz,fr
                        self.system.elec[e].update(px=0,
                                        py=0,
                                        pz=0,
                                        pr=0,
					spin=0,
                                        energy=self.ewald[idx].energy,
                                        fx=self.ewald[idx].fx,
                                        fy=self.ewald[idx].fy,
                                        fz=self.ewald[idx].fz,
                                        fr=self.ewald[idx].falpha*(-4.)/(self.system.elec[e].pr**3.))
					#fr=self.ewald[idx].falpha)
					#fr=-1*self.ewald[idx].falpha*1/((-4.)/(self.system.elec[e].pr**3.)))
			idx = idx  + 1

	def make_r_space(self):
		# generate real (r) space grid 
		Lx = self.Lx
		Ly = self.Ly
		Lz = self.Lz
		max_r = self.max_r
		dx = int(max_r/Lx+1)
		dy = int(max_r/Ly+1)
		dz = int(max_r/Lz+1)
		
		rgrid_x = [0]
		rgrid_y = [0]
		rgrid_z = [0]
		for nx in range(-dx,dx+1):
		        for ny in range(-dy,dy+1):
		                for nz in range(-dz,dz+1):
		                        if nx != 0 or ny != 0 or nz != 0:
			                        rgrid_x.append(nx*Lx)
			                        rgrid_y.append(ny*Ly)
			                        rgrid_z.append(nz*Lz)

	        # number of grid points (gps)
                #num_gps =  int((2*dx)*(2*dy)*(2*dz)	
		num_gps = len(rgrid_x)
		self.rgrid_x = rgrid_x
		self.rgrid_y = rgrid_y
		self.rgrid_z = rgrid_z
		self.num_gps = num_gps 

		return [rgrid_x,rgrid_y,rgrid_z]

	def make_k_space(self):
                # generate reciproke k space grid 
                Lkx = 2.*np.pi/self.Lx
                Lky = 2.*np.pi/self.Ly
                Lkz = 2.*np.pi/self.Lz
                max_k = self.max_k
                dkx = int(max_k/Lkx)
                dky = int(max_k/Lky)
                dkz = int(max_k/Lkz)

                kgrid_x = []
                kgrid_y = []
                kgrid_z = []
		kgrid_k2 = []
                for nkx in range(-dkx,dkx+1):
                        for nky in range(-dky,dky+1):
                                for nkz in range(-dkz,dkz+1):
                                        if nkx != 0 or nky != 0 or nkz != 0:
						ksquared = (nkx*Lkx)**2.+(nky*Lky)**2.+(nkz*Lkz)**2.
						#print ksquared
						#print self.max_k**2.
						if ksquared < self.max_k**2.: 
							#print 'yes'
                                                	kgrid_x.append(nkx*Lkx)
                                                	kgrid_y.append(nky*Lky)
                                                	kgrid_z.append(nkz*Lkz)
							kgrid_k2.append(ksquared)

		# number of k points (kps)
                #num_kps =  int((2*dkx)*(2*dky)*(2*dkz))
		num_kps = len(kgrid_x)
                self.kgrid_x = kgrid_x
                self.kgrid_y = kgrid_y
                self.kgrid_z = kgrid_z
		self.kgrid_k2 = kgrid_k2
                self.num_kps = num_kps

                return [kgrid_x,kgrid_y,kgrid_z]

        def check_boundary(self):
                def wrap(xi,Li):
                        shift = 0
                        if xi < 0:
                            shift = + Li 
                        if xi > Li:
                            shift = - Li
                        return shift 
                for n in range(len(self.system.nuc)):
                        shift_x = wrap(self.system.nuc[n].px,self.Lx) 
                        shift_y = wrap(self.system.nuc[n].py,self.Ly)
                        shift_z = wrap(self.system.nuc[n].pz,self.Lz) 
			print('%0.9f %0.9f  %0.9f ' % (shift_x,shift_y,shift_z)) 
                        self.system.nuc[n].update(px=shift_x,
                                          py=shift_y,
                                          pz=shift_z,
                                          pr=0,
                                          spin=0,
                                          energy=0,
                                          fx=0,
                                          fy=0,
                                          fz=0,
                                          fr=0)
                for e in range(len(self.system.elec)):
			shift_x = wrap(self.system.elec[e].px,self.Lx)
			shift_y = wrap(self.system.elec[e].py,self.Ly)
			shift_z = wrap(self.system.elec[e].pz,self.Lz)
                        print('%0.9f %0.9f  %0.9f ' % (shift_x,shift_y,shift_z))
			self.system.elec[e].update(px=shift_x,
                                          py=shift_y,
                                          pz=shift_z,
                                          pr=0,
                                          spin=0,
                                          energy=0,
                                          fx=0,
                                          fy=0,
                                          fz=0,
                                          fr=0)

	def show_all(self):
                print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('----------','----------','----------','----------','----------','----------')
                print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('typ','energy','fx','fy','fz','falpha')
                print '%10s\t%10s\t%10s\t%10s\t%10s\t%10s' %('----------','----------','----------','----------','----------','----------')
                for e in self.ewald:
                        e.show()

	def show_details(self):
		print 'Total volume %0.5f' % self.get_volume()
        	print 'Total charge %i' % self.get_total_charge()
        	print 'Ewald r cutoff %0.9f' % self.get_ewald_r_cutoff()
        	print 'Ewald k cutoff %0.9f' % self.get_ewald_k_cutoff()
        	print 'Number of real space grid points %i' % self.num_gps
        	print 'Number of k space grid points %i' % self.num_kps


if __name__ == "__main__":
	from pyeff_system import * 
        def main():
                # functionality test
		# CH4 example for charge neutrality  
                nuc = [6,1,1,1,1]
                elec = ['He','He','He','He','He','X','X','X','X','X']
                sys = system(nuc,elec)
                print sys.nuc[0].label
                print len(sys.elec)
                (sys.elec[0]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		(sys.elec[1]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		(sys.elec[2]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		(sys.elec[3]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		(sys.elec[4]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		(sys.elec[5]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
                (sys.elec[6]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
                (sys.elec[7]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
                (sys.elec[8]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
                (sys.elec[9]).update(px=1,py=0,pz=0,pr=1,spin=0,energy=0,fx=0,fy=0,fz=0,fr=0)
		sys_ewald = ewald_system(Lx=5,Ly=5,Lz=5,system=sys)
		print sys_ewald.get_volume()
		print sys_ewald.get_total_charge()
		print sys_ewald.get_ewald_r_cutoff()
		print sys_ewald.get_ewald_k_cutoff()
		print sys_ewald.max_k
		print sys_ewald.make_r_space()
		sys_ewald.show_all()
        main()

