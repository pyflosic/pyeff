import numpy as np
from scipy.special import erf
from pyeff_energy_force import system 
from numba import jit

# notes 
#--------------------------------------------------------------------
# erf           ...     takes complex arguments
# deriverative  ...     d/dz erf(z) = 2/ sqrt(Pi) e^(-z^2)

zero = 0.0000001
np.seterr(all='ignore')

# depends on system and ewald_system 

def read_cfg_pbc(p_cfg,print_data=None):
        # 
        # - reads the orginal eff code (cfg file) input for eff 
        # 
        f = open(p_cfg,'r')
        ll = f.readlines()
        # search tags 
        # nuclei positions
        tag_nuc = '@nuclei'
        # electron positions 
        tag_elec = '@electrons'
	# x_bound = 0.000000 8.352000
	# y_bound = 0.000000 8.352000
	# z_bound = 0.000000 8.352000
	# periodic = true
	tag_Lx = 'x_bound'
	tag_Ly = 'y_bound'
	tag_Lz = 'z_bound'
        for l in range(len(ll)):
                if ll[l].find(tag_nuc) != -1:
                        idx_nuc_pos_start = l+1
                if ll[l].find(tag_elec) != -1:
                        idx_elec_pos_start = l+1
		if ll[l].find(tag_Lx) != -1:
			Lx = float(ll[l].split()[-1])
		if ll[l].find(tag_Ly) != -1:
                        Ly = float(ll[l].split()[-1])
		if ll[l].find(tag_Lz) != -1:
                        Lz = float(ll[l].split()[-1])
        Z = []
        R = []
        r = []
        s = []
        delta = []
        for n in range(idx_nuc_pos_start,idx_elec_pos_start-1):
                # format p[0],p[1],p[2],element 
                tmp =  ll[n].split()
                R.append([float(tmp[0]),float(tmp[1]),float(tmp[2])])
                Z.append(float(tmp[3]))
        for e in range(idx_elec_pos_start,len(ll)):
                # format p[0],p[1],p[2],spin,eradius    
                tmp = ll[e].split()
                r.append([float(tmp[0]),float(tmp[1]),float(tmp[2])])
                delta.append(float(tmp[3]))
                s.append(float(tmp[4]))
        f.close()
        if print_data != None:
                print Z
                print R
                print r
                print s
                print delta
        # creates a systems object      
        nuc = Z
        elec = delta
        calc = system(nuc,elec)
        for n in range(len(R)):
                        calc.nuc[n].px = R[n][0]
                        calc.nuc[n].py = R[n][1]
                        calc.nuc[n].pz = R[n][2]
        for e in range(len(r)):
                        calc.elec[e].px = r[e][0]
                        calc.elec[e].py = r[e][1]
                        calc.elec[e].pz = r[e][2]
                        calc.elec[e].pr = s[e]
                        calc.elec[e].spin = delta[e]
        return [Lx,Ly,Lz,calc,Z,R,r,s,delta]



# E_Ewald = E_k-space + E_r-space + E_self + E_uniform 

# E_r-space = \sum_R_{L} \sum_{i\neqj} E_{ij}(r_{ij}-R_{L})
# E_{ij} = 1/r_{ij} erf(\sqrt{\alpha_i*\alpha_j/(\alpha_i+\alpha_j)}r_{ij}
#        - 1/r_{ij} erf(\sqrt{\alpha_i*\alpha_j^{max}/(\alpha_i+\alpha_j^{max})}r_{ij}
# only for interaction between a charge an onther charge with \alpha > \alpha_{max}

@jit(nopython=True, parallel=True)
def E_rspace(es):
	# es ... ewald_system 
	E = []
	E_i = 0
	count = 0
        for i in range(es.numcharges):
		E_i = 0
                for j in range(es.numcharges):
			if (es.ewald[j].palpha > es.max_alpha):
				a1  = np.sqrt(es.ewald[i].palpha*es.ewald[j].palpha/(es.ewald[i].palpha+es.ewald[j].palpha))
                                a2  = np.sqrt(es.ewald[i].palpha*es.max_alpha/(es.ewald[i].palpha+es.max_alpha))
				for gp in range(es.num_gps):
					if i !=j or gp != 0:
						dx = es.ewald[i].px - es.ewald[j].px  - es.rgrid_x[gp]
                                                dy = es.ewald[i].py - es.ewald[j].py  - es.rgrid_y[gp]
                                                dz = es.ewald[i].pz - es.ewald[j].pz  - es.rgrid_z[gp]
                                                r  = np.sqrt(dx**2. + dy**2. + dz**2.)
						if (r < es.max_r):
							# this is important!
							if r == 0:
								r = zero
							# E_rspace = qi*qj*(erf(r*sqrt(alphai*alphaj/(alphai + alphaj)))/r - erf(r*sqrt(alphai*alphamax/(alphai + alphamax)))/r)
							E1 = erf(r*a1)/r
							E2 = -1.*erf(r*a2)/r
							energy = (es.ewald[i].q*es.ewald[j].q)*(E1+E2)
							E_i = E_i + energy
							# Debugging 
							# print('i j gp f1 f2 energy[i] %i %i %i %0.6f  %0.6f %0.6f' % (i,j,gp,E1/a1,E2/a2,E_i))
							
							# dE_rspace/dr = -qi*qj*(dE1+dE2+dE3)
							# dE1 = 2.*sqrt(alphai*alphaj/(alphai + alphaj))*exp(-alphai*alphaj*r**2/(alphai + alphaj))/(sqrt(pi)*r)
							dE1 = 2.*a1*np.exp(-1.*(a1**2.)*(r**2.))/(np.sqrt(np.pi)*r)
							# dE2 = -2.*sqrt(alphai*alphamax/(alphai + alphamax))*exp(-alphai*alphamax*r**2/(alphai + alphamax))/(sqrt(pi)*r)
							dE2 = -2.*a2*np.exp(-1.*(a2**2.)*(r**2.))/(np.sqrt(np.pi)*r)
							# dE3 = -1.*erf(r*sqrt(alphai*alphaj/(alphai + alphaj)))/r**2 + erf(r*sqrt(alphai*alphamax/(alphai + alphamax)))/r**2
							dE3 = -1.*erf(r*a1)/r**2. + erf(r*a2)/r**2.
							dE_rspacedr = es.ewald[i].q*es.ewald[j].q*(dE1+dE2+dE3)
	
							fr = dE_rspacedr/r 
							fx = fr*dx
							fy = fr*dy
							fz = fr*dz

							# falpha_i = -1*qi*qj(dE1dalphai + dE2dalphai)
							# dE1dalphai = -alphaj*sqrt(alphai*alphaj/(alphai + alphaj))*exp(-alphai*alphaj*r**2/(alphai + alphaj))/(sqrt(pi)*alphai*(alphai + alphaj))
							# dE2dalphai = alphamax*sqrt(alphai*alphamax/(alphai + alphamax))*exp(-alphai*alphamax*r**2/(alphai + alphamax))/(sqrt(pi)*alphai*(alphai + alphamax))
							dE1dalphai = -1*es.ewald[j].palpha*a1*np.exp(-1*(a1**2.)*(r**2.))/(np.sqrt(np.pi)*es.ewald[i].palpha*(es.ewald[i].palpha + es.ewald[j].palpha))
							dE2dalphai = es.max_alpha*a2*np.exp(-1*(a2**2)*(r**2))/(np.sqrt(np.pi)*es.ewald[i].palpha*(es.ewald[i].palpha + es.max_alpha))
							falpha_i = -1*es.ewald[i].q*es.ewald[j].q*(dE1dalphai + dE2dalphai)

							# falpha_j = -1*qi*qj(dE1dalphaj + dE2dalphaj)
							# dE1dalphaj = -alphai*sqrt(alphai*alphaj/(alphai + alphaj))*exp(-alphai*alphaj*r**2/(alphai + alphaj))/(sqrt(pi)*alphaj*(alphai + alphaj))
							# dE2dalphaj = 0

							dE1dalphaj = -es.ewald[i].palpha*a1*np.exp(-1.*(a1**2.)*(r**2))/(np.sqrt(np.pi)*es.ewald[j].palpha*(es.ewald[i].palpha + es.ewald[j].palpha))
							dE2dalphaj = 0
							falpha_j =-1* es.ewald[i].q*es.ewald[j].q*(dE1dalphaj)
							
							# update  
							es.ewald[i].update(px=0,py=0,pz=0,palpha=0,q=0,energy=0.5*energy,fx=-0.5*fx,fy=-0.5*fy,fz=-0.5*fz,falpha=-0.5*falpha_i)
                                                        es.ewald[j].update(px=0,py=0,pz=0,palpha=0,q=0,energy=0,fx=(-0.5)*-1*fx,fy=(-0.5)*-1*fy,fz=(-0.5)*-1*fz,falpha=(-0.5)*falpha_j)
							
	return energy
						
# E_k-space = 2Pi/V \sum_k \rho(k)*\rho_max(k)  
# \rho(k) = \sum_i q_i e^{-k^2/(4*\alpha_i)}*e^{-i*k*r_1}
# \rho_max(k) = \rho(k)|_{\alpha_i=\alpha_{max}} for \alpha_i > \alpha_{max}
 
#@jit(nopython=True, parallel=True)
def E_kspace(es):
	# 1st calculate densities 
	ewald_energy = 0
	ewald_S_real = [] 
	ewald_S_imag = [] 
        ewald_S_cap_real = [] 
        ewald_S_cap_imag = [] 
	rho_k_real = 0
	rho_k_imag = 0
	rho_k_cap_real = 0
	rho_k_cap_imag = 0
	V = es.get_volume()
	fac_fr = 2.*np.pi/V
	fac_falpha = -0.5*np.pi/V
	for i in range(es.num_kps):
		ewald_S_real = []
	        ewald_S_imag = [] 
       	 	ewald_S_cap_real = []
        	ewald_S_cap_imag = [] 
		rho_k_real = 0
		rho_k_imag = 0
		rho_k_cap_real = 0
		rho_k_cap_imag  = 0
		for j in range(es.numcharges):
			angle =-1.*(es.kgrid_x[i]*es.ewald[j].px+es.kgrid_y[i]*es.ewald[j].py+es.kgrid_z[i]*es.ewald[j].pz) 
			# print angle 
			mag = es.ewald[j].q*np.exp(-1*es.kgrid_k2[i]/(4.*es.ewald[j].palpha))
			ewald_S_real.append(mag*np.cos(angle))
			ewald_S_imag.append(mag*np.sin(angle))
			if es.ewald[j].palpha > es.max_alpha:
				mag_cap = es.ewald[j].q*np.exp(-1.*es.kgrid_k2[i]/(4.*es.max_alpha))
				ewald_S_cap_real.append(mag_cap*np.cos(angle))
				ewald_S_cap_imag.append(mag_cap*np.sin(angle))
			else:
				mag_cap = mag 
				#ewald_S_cap_real.append(ewald_S_real)
				ewald_S_cap_real.append(mag*np.cos(angle))
				#ewald_S_cap_imag.append(ewald_S_imag)
				ewald_S_cap_imag.append(mag*np.sin(angle))
			rho_k_real = rho_k_real + ewald_S_real[j]
			rho_k_imag = rho_k_imag + ewald_S_imag[j]
			rho_k_cap_real = rho_k_cap_real + ewald_S_cap_real[j]
                        rho_k_cap_imag = rho_k_cap_imag + ewald_S_cap_imag[j]


		# 2nd calculate energy and forces with these densities 	
		j = 0 
		for j in range(es.numcharges):
			E_kspacedr = ewald_S_real[j]*rho_k_cap_imag-ewald_S_imag[j]*rho_k_cap_real+ewald_S_cap_real[j]*rho_k_imag-ewald_S_cap_imag[j]*rho_k_real
			E_kspacedr = E_kspacedr/es.kgrid_k2[i]
			fx = E_kspacedr*es.kgrid_x[i]
			fy = E_kspacedr*es.kgrid_y[i]
			fz = E_kspacedr*es.kgrid_z[i]  
			if es.ewald[j].palpha > es.max_alpha: 
				falpha = (ewald_S_real[j]*rho_k_cap_real+ewald_S_imag[j]*rho_k_cap_imag)/(es.ewald[j].palpha**2.) 

			else: 	
				falpha = (ewald_S_real[j]*rho_k_cap_real+ewald_S_imag[j]*rho_k_cap_imag+ewald_S_cap_real[j]*rho_k_real+ewald_S_cap_imag[j]*rho_k_imag)
			 	falpha = falpha/(es.ewald[j].palpha**2.)

			# update forces 	
			fx = fac_fr*fx 
			fy = fac_fr*fy
			fz = fac_fr*fz
			falpha = fac_falpha*falpha 
			es.ewald[j].update(px=0,py=0,pz=0,palpha=0,q=0,energy=0,fx=fx,fy=fy,fz=fz,falpha=falpha)
		
		ewald_i = (1.0 / es.kgrid_k2[i]) * (rho_k_real * rho_k_cap_real + rho_k_imag * rho_k_cap_imag)
		ewald_energy = ewald_energy + ewald_i
		#print('%0.7f' % ewald_energy)
	ewald_energy = ewald_energy * (2*np.pi/V)/float(es.numcharges)#/es.num_kps
	# print('%0.7f' % ewald_energy)
	i = 0
	# update energy 
	for i in range(es.numcharges):
		es.ewald[i].update(px=0,py=0,pz=0,palpha=0,q=0,energy=ewald_energy,fx=0,fy=0,fz=0,falpha=0)
	return ewald_energy

#
# E_self 
#
# The self energy is the energy of a Gaussian charge with itself - this quantity
# must be subtracted out, since in reality, charges do not repel themselves.
#
#
# E_self = -1/2 \sum_i q_i^2 2/\pi *
# \sqrt{\alpha_i*\alpha_{max}/(\alpha_i+\alpha_max)} for \alpha_i < \alpha_{max}
# \sqrt{\alpha_i*\alpha_i/(\alpha_i+\alpha_i)} otherwise 

#@jit(nopython=True, parallel=True)
def E_self(es):
	# self energy 
	E_self_sum = 0
	for i in range(es.numcharges):
		if es.ewald[i].palpha >= es.max_alpha:
			a2 = np.sqrt(es.ewald[i].palpha*es.max_alpha/(es.ewald[i].palpha + es.max_alpha))
			# sympy 
			# Eself1 = -1.0*qi**2.0*sqrt(alphai*alphamax/(alphai + alphamax))/sqrt(pi)
			# dEself1dalphai = 0.5*alphamax*qi**2.0*sqrt(alphai*alphamax/(alphai + alphamax))/(sqrt(pi)*alphai*(alphai + alphamax))
			Eself = -1.0*(es.ewald[i].q**2.)*a2/np.sqrt(np.pi)
			dEselfdalphai = 0.5*es.max_alpha*(es.ewald[i].q**2.)*a2/(np.sqrt(np.pi)*es.ewald[i].palpha*(es.ewald[i].palpha + es.max_alpha))
			# eff code 
			# dEselfdalphai =(es.ewald[i].q**2.)*(0.5 * a2 * a2 * a2 / (es.ewald[i].palpha**2.)) / np.sqrt(np.pi)
		else: 
			# sympy
			# Eself2 = -0.5*sqrt(2)*sqrt(alphai)*qi**2.0/sqrt(pi)
                        # dEself2dalphai = 0.25*sqrt(2)*qi**2.0/(sqrt(pi)*sqrt(alphai))
			Eself = -0.5*np.sqrt(2)*np.sqrt(es.ewald[i].palpha)*(es.ewald[i].q**2.)/np.sqrt(np.pi)
			dEselfdalphai =	0.25*np.sqrt(2)*(es.ewald[i].q**2.)/(np.sqrt(np.pi)*np.sqrt(es.ewald[i].palpha))
			# eff code 
			# a1 = np.sqrt(es.ewald[i].palpha/2.)
			# dEselfdalphai = (es.ewald[i].q**2.) / (4. * a1) / np.sqrt(np.pi)
		es.ewald[i].update(px=0,py=0,pz=0,palpha=0,q=0,energy=Eself,fx=0,fy=0,fz=0,falpha=dEselfdalphai)
		E_self_sum = E_self_sum + Eself 
	return E_self_sum 
	

# 
# E_uniform 
#
# The uniform energy represents the interaction of each Gaussian charge with the
# uniform neutralizing background charge. For neutral systems (Q=0) this energy has no effect. 
# 
 
# E_{unifrom} = 1/4*Q*\sum_i \pi/V (1/\alpha_{max} - 1/\alpha_{i}) for \alpha_i> \alpha_{max}

#@jit(nopython=True, parallel=True)			
def E_uniform(es):
	
	E_uniform_sum = 0
	V = es.get_volume()
	Q = es.get_total_charge()
	#print 'Q = %i' % Q 
	fac = -0.5*Q*np.pi
	
	for i in range(es.numcharges):
		if es.ewald[i].palpha >= es.max_alpha:
			E_uniform = fac *es.ewald[i].q*(1/es.max_alpha - 1/es.ewald[i].palpha)
			dE_uniformdalphai = -1*fac* es.ewald[i].q / (es.ewald[i].palpha**2.)
		else:
			E_uniform = 0
			dE_uniformdalphai = 0
		E_uniform_sum = E_uniform_sum + E_uniform 
		es.ewald[i].update(px=0,py=0,pz=0,palpha=0,q=0,energy=E_uniform,fx=0,fy=0,fz=0,falpha=dE_uniformdalphai)
	return E_uniform_sum 

def bound(di,Li):
	#if abs(di) < -0.5*abs(Li):
	#	res = 1*(abs(di) + abs(Li))
	#if abs(di) > 0.5*abs(Li):
	#	res = 1*(abs(di) - abs(Li)) 
	#print('di Li %0.9f  %0.9f' % (di,Li))
	if di < -0.5*abs(Li):
		#print(1)
		res = di + abs(Li)
	elif di > 0.5*abs(Li):
		#print(2)
                res =  di - abs(Li)
	else: 
		#print(3)
		res = di  
	return res 

#@jit(nopython=True, parallel=True)
def E_PauliPeriodic(es):
	# this function is similar to the E_Pauli function 
	# but includes minima image coniditions for periodic boundary conditions (pbcs) 

	Pauli_rho = -0.2
	Pauli_r = 1.125
	Pauli_s = 0.9


	def E_up_up(rho,r,si,sj):
                # from sympy calculations
                E_up_up =  (8.0*(-rho + 1.0)*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) + 8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(-8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0))*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)
                # from sympy calculations
                dE_up_updr = -1*1.125*r*(-8.0*si**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(256.0*(si*sj/(si**2 + sj**2))**3.0 - 32.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**6.0*(-256.0*rho + 256.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 256.0*(si*sj/(si**2 + sj**2))**6.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**2*(si**2 + sj**2)**3*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_updr) == True:
                        dE_up_updr = 0

                # from sympy calculations 
                dE_up_upds1 =  -1*0.9*(sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*si**4*(si**2 + sj**2) + si**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(-rho + 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) - 32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*si**2 + 24.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(256.0*r**2*si**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**3*sj**2*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_upds1) == True:
                        dE_up_upds1 = 0

                # from sympy calculations 
                dE_up_upds2 =  -1*0.9*(si**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*sj**4*(si**2 + sj**2) + sj**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) - (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(rho - 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*sj**2 + 24.0*(si**2 - sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(-256.0*r**2*sj**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**3*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_upds2) == True:
                        dE_up_upds2 = 0

                return E_up_up,dE_up_updr,dE_up_upds1,dE_up_upds2
	
	def E_up_down(rho,r,si,sj):
                # from sympy calculation 
                E_up_down = -8.0*rho*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0)

                # from sympy calculations 
                dE_up_downdr = -1*1.125*r*rho*(-64.0*si**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 32.0*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + (si*sj/(si**2 + sj**2))**6.0*(si**2*sj**2*(-1024.0*r**2 + 1536.0*si**2 + 1536.0*sj**2) - 384.0*si**2*(si**2 + sj**2)**2 - 384.0*sj**2*(si**2 + sj**2)**2))/(si**2*sj**2*(si**2 + sj**2)**3*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_downdr) == True:
                        dE_up_downdr = 0

                # from sympy calculations 
                dE_up_downds1 = -1*0.9*rho*(-32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*si**4*(si**2 + sj**2) + si**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*(si*sj/(si**2 + sj**2))**6.0*(32.0*r**2*si**2 + 24.0*(-si**2 + sj**2)*(si**2 + sj**2))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2))/(si**3*sj**2*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_downds1) == True:
                        dE_up_downds1 = 0

                # from sympy calculations 
                dE_up_downds2 = -1.0*0.9*rho*(-32.0*r**2*sj**2.0*(si*sj/(si**2.0 + sj**2.0))**3.0*(8.0*(si*sj/(si**2.0 + sj**2.0))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*si**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*sj**4*(si**2 + sj**2) + sj**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3) - 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*(si*sj/(si**2 + sj**2))**6.0*(32.0*r**2*sj**2 + 24.0*(si**2 - sj**2)*(si**2 + sj**2))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2))/(si**2*sj**3*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
                if np.isnan(dE_up_downds2) == True:
                        dE_up_downds2 = 0


                return E_up_down,dE_up_downdr,dE_up_downds1,dE_up_downds2

	for i in range(len(es.system.elec)):
		x1 = es.system.elec[i].px 
		y1 = es.system.elec[i].py
		z1 = es.system.elec[i].pz
		s_i = es.system.elec[i].pr
		spin_i = es.system.elec[i].spin 
		for j in range(0,i):
			#
			# minima image
			# look if the positions is outside the simulation box 
			# if it is outside wrap the positions back in the box 
			#
			#print(es.Lx)
			#print(es.Ly)
			#print(es.Lz)
			dx = es.system.elec[j].px - x1
			dy = es.system.elec[j].py - y1
			dz = es.system.elec[j].pz - z1
			x_ij = np.sqrt(dx**2.+dy**2.+dz**2.)
			#print('pure dx dy dz r %0.5f %0.5f %0.5f %0.5f' % (dx,dy,dz,x_ij))	
			dx = bound(dx,es.Lx)
			dy = bound(dy,es.Ly)
                        dz = bound(dz,es.Lz)
			x_ij = np.sqrt(dx**2.+dy**2.+dz**2.) 
			#print('minima_image dx dy dz r %0.5f %0.5f %0.5f %0.5f' % (dx,dy,dz,x_ij))
			if x_ij == 0:
                        	x_ij = zero
			s_j = es.system.elec[j].pr
			spin_j = es.system.elec[j].spin 
			if spin_i == spin_j:
				E1, E1dr, E1ds1, dE1ds2 = E_up_up(Pauli_rho,Pauli_r*x_ij,Pauli_s*s_i,Pauli_s*s_j)
				#print E1
				fr = E1dr
				fs1 = E1ds1
				fs2 = dE1ds2
				e = E1
				m_x_ij = x_ij
				fx = np.float64(fr/m_x_ij * dx)
				fx = fr/m_x_ij * dx
				fy = fr/m_x_ij * dy
				fz = fr/m_x_ij * dz
				es.system.elec[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=fs1)
				es.system.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs2)
				#print('%i %i elec[i] E fx dy dz dr %0.5f %0.5f %0.5f %0.5f %0.5f' % (i,j,e/2.,-1*fx,-1*fy,-1*fz,fs1))
				#print('%i %i elec[j] E fx dy dz dr %0.5f %0.5f %0.5f %0.5f %0.5f' % (i,j,e/2.,fx,fy,fz,fs2))
			if spin_i != spin_j:
				E2, E2dr, E2ds1, dE2ds2 = E_up_down(Pauli_rho,Pauli_r*x_ij,Pauli_s*s_i,Pauli_s*s_j)
                                #print E2
                                fr = E2dr
                                fs1 = E2ds1
                                fs2 = dE2ds2
                                e = E2
                                m_x_ij = x_ij
                                fx = np.float64(fr/m_x_ij * dx)
                                fx = fr/m_x_ij * dx
                                fy = fr/m_x_ij * dy
                                fz = fr/m_x_ij * dz
                                es.system.elec[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=fs1)
                                es.system.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs2)
				#print('%i %i elec[i] E fx dy dz dr %0.5f %0.5f %0.5f %0.5f %0.5f' % (i,j,e/2.,-1*fx,-1*fy,-1*fz,fs1))
				#print('%i %i elec[j] E fx dy dz dr %0.5f %0.5f %0.5f %0.5f %0.5f' % (i,j,e/2.,fx,fy,fz,fs2))

	
if __name__ == "__main__":
        from pyeff_system import *
	from pyeff_ewald_system_numba import *
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
		print sys_ewald.make_r_space()
		print sys_ewald.num_gps 
		np.set_printoptions()
		print len(sys_ewald.rgrid_x)
		print sys_ewald.make_k_space()
		print E_rspace(sys_ewald)
		sys_ewald.show_all()
		print E_kspace(sys_ewald)
		sys_ewald.show_all()
		print E_self(sys_ewald)
                sys_ewald.show_all()
        main()

