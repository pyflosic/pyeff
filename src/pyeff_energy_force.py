# pyeff - A birthday present for Prof. Jens Kortus -  
# electron force field (eff) type eff1 in python 
# author: Sebastian Schwalbe
# co-authors: Simon Liebing  
# date:   14.09.2017 

import os
import numpy as np 
from scipy.special import erf
from pyeff_system import * 
#from pyeff_sympy_functions import * 

# global variables 
# numerical zero 
# note: do not make it to small !
zero = 0.0000001 # to avoid divide by zero error 

# variables 
#--------------------------------------------------------------------
# Z 		... nuclei charges 
# R 		... nuclei positions 
# r 		... electron positions
# s 		... size (radius?)  
# delta		... spin value for each electron (+1/-1) 

# notes 
#--------------------------------------------------------------------
# erf           ...     takes complex arguments 
# deriverative  ...	d/dz erf(z) = 2/ sqrt(Pi) e^(-z^2)
# eff-code 	... 	the main physics is in the eff_update.c file 
# lammps eff	...	v_estatic = Enucnuc + Enucelec + Eelecelc  

# TO-DO
#--------------------------------------------------------------------
# (1) implement forces for eff1 
# (2) implement a calculator object for ase, which makes optimization 
#     and molecular dynamics possible 
# (3) implement eff2 and eff3 of the original proposal (eff_code)
# (4) implement eff ecp mode 

def read_cfg(p_cfg,print_data=None):
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
	for l in range(len(ll)):
		if ll[l].find(tag_nuc) != -1:
			idx_nuc_pos_start = l+1 
		if ll[l].find(tag_elec) != -1:
                        idx_elec_pos_start = l+1
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
 	return [calc,Z,R,r,s,delta] 

def read_eff(p_eff,print_data=None):
        # 
        # - reads the orginal eff code (eff file) for comparison with pyeff 
        # - returns a calculation object 
	#
	tag_nuc = '[nuc]'
	tag_elec = '[elec]'
	tag_nuc_pos = '[position_nuc]'
	tag_elec_pos = '[position_elec]'
	tag_nuc_frc = '[energy_force_nuc]'
	tag_elec_frc = '[energy_force_elec]'
	f = open(p_eff,'r')
        ll = f.readlines()
	elec = []
	nuc = [] 
	for l in range(len(ll)):
		if ll[l].find(tag_nuc) != -1:
			tmp_nuc = ll[l].split()[-2]
			nuc.append(tmp_nuc)
		if ll[l].find(tag_elec) != -1:
			tmp_elec = ll[l].split()[-1]
			elec.append(tmp_elec)
        calc = system(nuc,elec)
	# positions 
	n = 0
        e = 0
        for l in range(len(ll)):
                if ll[l].find(tag_nuc_pos) != -1:
                        tmp = ll[l].split()
			calc.nuc[n].px = float(tmp[2])
                        calc.nuc[n].py = float(tmp[3])
                        calc.nuc[n].pz = float(tmp[4])
                        n = n +1
                if ll[l].find(tag_elec_pos) != -1:
                        tmp = ll[l].split()
                        calc.elec[e].px = float(tmp[2])
                        calc.elec[e].py = float(tmp[3])
                        calc.elec[e].pz = float(tmp[4])
                        calc.elec[e].pr = float(tmp[5])
                   	calc.elec[e].spin = elec[e]
			e = e +1
	# energy and forces 
        n = 0 
	e = 0 
	for l in range(len(ll)):
		if ll[l].find(tag_nuc_frc) != -1:
			tmp = ll[l].split()
			calc.nuc[n].energy = float(tmp[2])
			calc.nuc[n].fx = float(tmp[3])
			calc.nuc[n].fy = float(tmp[4])
			calc.nuc[n].fz = float(tmp[5])
			n = n +1 
		if ll[l].find(tag_elec_frc) != -1:
			tmp = ll[l].split()
                        calc.elec[e].energy = float(tmp[2])
                        calc.elec[e].fx = float(tmp[3])
                        calc.elec[e].fy = float(tmp[4])
                        calc.elec[e].fz = float(tmp[5])
			calc.elec[e].fr = float(tmp[6])
                        e = e +1
	f.close() 
	return calc 

def E_ke(calc,s,debug=None):
	'$E_ke = \sum_i 3/2 1/s_i^2$'
	E_ke = 0
	F_ke = []  
	for i in range(len(s)):
		# energy 
		e = (3./2.)*(1./(s[i]**2))
		# force  
		fs = 3 *1./(s[i]**3)
		if debug != None:
			fs_r = sympy_Eke(s[i])
			print 'pyeff fs =\t', fs 
			print 'analytic fs =\t',fs_r[1]
			print 'numerical fs =\t',fs_r[3] 
		E_ke = E_ke + e #+ e2  	
		F_ke.append(fs)
		# update the calculations object (calc) 
		calc.elec[i].update(energy=e,px=0,py=0,pz=0,pr=0,spin=0,fx=0,fy=0,fz=0,fr=fs)
	return E_ke 

def E_nuc_nuc(calc,Z,R,debug=None):
	E_nuc_nuc = 0
	F_nuc_nuc = []  
	for i in range(len(Z)):
		F_nuc = []
		for j in range(len(Z)):
			if i < j:
				# Energy contribution
				dx = R[j][0]-R[i][0]
				dy = R[j][1]-R[i][1]
				dz = R[j][2]-R[i][2]
				R_ij = np.sqrt(abs((dx)**2+(dy)**2+(dz)**2))
				e = Z[i]*Z[j]/R_ij 
				E_nuc_nuc = E_nuc_nuc + e 
				# Force contribution 
				fr =  Z[i]*Z[j]/(R_ij**2)
				if debug != None:
					fr_r = sympy_E_nucnuc(Z[i],Z[j],R_ij)
                			print 'pyeff fr =\t', fr
                			print 'analytic fr =\t',fr_r[1]
                			print 'numerical fr =\t',fr_r[3]
				fx = fr/R_ij * dx
				fy = fr/R_ij * dy
				fz = fr/R_ij * dz
				F_n = [fx,fy,fz] 
				F_nuc_nuc.append(F_n)
				# update the calculations object (calc) 
                		# note: compared to eff-code update for i and j are switched
				calc.nuc[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=0)
				calc.nuc[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=0)

	return E_nuc_nuc 

def E_nuc_elec(calc,Z,R,r,s,debug=None):
	E_nuc_elec = 0
	F_nuc_elec = [] 
	for i in range(len(Z)):
		for j in range(len(r)):
			dx = R[i][0]-r[j][0]
			dy = R[i][1]-r[j][1]
			dz = R[i][2]-r[j][2]
			R_ij = np.sqrt(abs((dx)**2+(dy)**2+(dz)**2)) 
			# obmit divide by zero error 
			if R_ij == 0:
				R_ij = zero
			e = -1*Z[i]/R_ij*erf(np.sqrt(2.)*R_ij/s[j]) # s[i] or s[j] ?
			E_nuc_elec = E_nuc_elec + e
			# from sympy calc 
			# dE/dR 
			# dENEdR = Sum(2*sqrt(2)*Z(i)*exp(-2*R(i, j)**2/s(j)**2)/(sqrt(pi)*R(i, j)*s(j)) - Z(i)*erf(sqrt(2)*R(i, j)/s(j))/R(i, j)**2, (i, 1, n), (j, 1, m))
			fr = 2.*np.sqrt(2.)*Z[i]*np.exp(-2.*R_ij**2/s[j]**2.)/(np.sqrt(np.pi)*R_ij*s[j])-Z[i]*erf(np.sqrt(2.)*R_ij/s[j])/R_ij**2
			fr = -1*fr 
                        if debug != None:
				fr_r = sympy_E_nucelec(Z[i],R_ij,s[j])
                        	print 'pyeff fr =\t', fr
                        	print 'analytic fr =\t',fr_r[1]
                        	print 'numerical fr =\t',fr_r[3]
 			fx = fr/R_ij * dx
			fy = fr/R_ij * dy
			fz = fr/R_ij * dz
			F_ne = [fx,fy,fz]
			F_nuc_elec.append(F_ne)
			# from sympy calc 
			# dE/ds 
			fs = -2*np.sqrt(2)*Z[i]*np.exp(-2*R_ij**2/s[j]**2)/(np.sqrt(np.pi)*s[j]**2)
			# update the calculations object (calc) 
			calc.nuc[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=0)
			# for electron we need dE/ds 
			calc.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs)
	
	return E_nuc_elec

def E_elec_elec(calc,r,s,debug=None):
	E_elec_elec = 0
	F_elec_elec = []
	for i in range(len(r)):
		for j in range(len(r)):
			if j < i:
				dx = r[j][0]-r[i][0]
				dy = r[j][1]-r[i][1]
				dz = r[j][2]-r[i][2]
				x_ij = np.sqrt(abs((dx)**2+(dy)**2+(dz)**2))
				# obmit divide by zero error 
				if x_ij == 0: 
					x_ij = zero
				e = 1./x_ij*erf(np.sqrt(2.)*x_ij/(np.sqrt(s[i]**2+s[j]**2)))
				E_elec_elec = E_elec_elec + e
				# from sympy calc 
  				# dE/dR
				fr = (erf(np.sqrt(2)*x_ij/np.sqrt(s[i]**2 + s[j]**2))/x_ij - 2.*np.sqrt(2)*np.exp(-2.*x_ij**2/(s[i]**2 + s[j]**2))/(np.sqrt(np.pi)*np.sqrt(s[i]**2 + s[j]**2)))/x_ij
				if debug != None:
					fr_r = sympy_E_elecelec(x_ij,s[i],s[j])
                        		print 'pyeff fr =\t', fr
                        		print 'analytic fr =\t',fr_r[1]
                        		print 'numerical fr =\t',fr_r[3]

				fx = fr/x_ij * dx
        	                fy = fr/x_ij * dy
        	                fz = fr/x_ij * dz
				F_ee = [fx,fy,fz]
			        F_elec_elec.append(F_ee)
				# from sympy calc 
        	                # dE/ds1
				fs1 = 2.*np.sqrt(2)*s[i]*np.exp(-2.*x_ij**2/(s[i]**2 + s[j]**2))/(np.sqrt(np.pi)*(s[i]**2 + s[j]**2)**(3./2.))
		                # from sympy calc 
        		        # dE/ds2
				fs2 = 2.*np.sqrt(2)*s[j]*np.exp(-2.*x_ij**2/(s[i]**2 + s[j]**2))/(np.sqrt(np.pi)*(s[i]**2 + s[j]**2)**(3./2.))
				# update the calculations object (calc) 
        	                # for electron we need dE/ds 
				# switch force update compared to eff-code 
				calc.elec[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=fs1)
				calc.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs2)
	

	return E_elec_elec


def E_Pauli(calc,delta,s,r,debug=None,parameters=None):
	if parameters == None:
		Pauli_rho = -0.2
		Pauli_r = 1.125 
		Pauli_s = 0.9 
	if parameters != None:
		Pauli_rho = parameters[0]
		Pauli_r = parameters[1]
		Pauli_s= parameters[2] 
	def Delta_T_ij(r,si,sj):
		# energy 
		Delta_T_ij = (3./2.)*((1./(si**2))+(1./(sj**2)))-2.*(3.*((si**2)+(sj**2))-2.*(r**2))/(((si**2)+(sj**2))**2)
		return Delta_T_ij
		
	def S_ij(r,si,sj):
		# energy 
		S_ij = ((2./((si/sj)+(sj/si)))**(3./2.))*np.exp((-1.*(r**2))/(si**2+sj**2))
		return S_ij
	
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

	E_1 = 0 # E_up_up contribution 
	E_2 = 0 # E_up_down contribution 
	F_1 = []
	F_2 = []
	for i in range(len(delta)):
		for j in range(len(delta)):
			if j < i:
				s_i = s[i]
               			s_j = s[j]
				dx = r[j][0]-r[i][0]
				dy = r[j][1]-r[i][1]
				dz = r[j][2]-r[i][2]
	                        x_ij = np.sqrt(abs((dx)**2+(dy)**2+(dz)**2))
	                        Tij = Delta_T_ij(x_ij,s_i,s_j)
	                        Sij = S_ij(x_ij,s_i,s_j)
				if delta[i] == delta[j]:
					if x_ij == 0:
						x_ij = zero 
					E1, E1dr, E1ds1, dE1ds2 = E_up_up(Pauli_rho,Pauli_r*x_ij,Pauli_s*s_i,Pauli_s*s_j)
					fr = E1dr 
					fs1 = E1ds1
					fs2 = dE1ds2 
					if debug != None:
						# analytic 
						E_r,fr_r,fs1_r,fs2_r = sympy_E_Pauli(x_ij,s_i,s_j,samespin=True,debug=True)
                                	        print 'pyeff S =\t',Sij
                                	        print 'pyeff Delta_T =\t',Tij
                                	        print 'pyeff fr =\t', fr
                                	        print 'analytic fr =\t',fr_r
					E_1 = E_1 + E1 
					e = E1 
					m_x_ij = x_ij
					fx = np.float64(fr/m_x_ij * dx)
					fx = fr/m_x_ij * dx
	                        	fy = fr/m_x_ij * dy
	                        	fz = fr/m_x_ij * dz
					F1 = [fx,fy,fz]
	                        	F_1.append(F1)
					# update the calculations object (calc) 
					calc.elec[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=fs1)
					calc.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs2)
					#print 'e1 = '+str(e)+' fx = '+str(fx)+' fy = '+str(fy)+' fz = '+str(fz)
				if delta[i] != delta[j]:
	                                if x_ij == 0:
	                                        x_ij = zero
					E2, E2dr,E2ds1,E2ds2 = E_up_down(Pauli_rho,Pauli_r*x_ij,Pauli_s*s_i,Pauli_s*s_j)
					fr = E2dr
                                        fs1 = E2ds1
                                        fs2 = E2ds2
					if debug != None:
						# analytic 
						E_r,fr_r,fs1_r,fs2_r = sympy_E_Pauli(x_ij,s_i,s_j,samespin=False,debug=True)
						print 'pyeff S =\t',Sij
						print 'pyeff Delta_T =\t',Tij
						print 'pyeff fr =\t', fr 
                                        	print 'analytic fr =\t',fr_r
                                        fr = E2dr 
					fs1 = E2ds1
					fs2 = E2ds2
					E_2 = E_2 + E2 
					e = E2
					m_x_ij = x_ij
	                                fx = fr/m_x_ij * dx
	                                fy = fr/m_x_ij * dy
	                                fz = fr/m_x_ij * dz
					F2 = [fx,fy,fz]
	                                F_2.append(F2)
					# update the calculations object (calc) 
					calc.elec[i].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=-1*fx,fy=-1*fy,fz=-1*fz,fr=fs1)
					calc.elec[j].update(energy=e/2.,px=0,py=0,pz=0,pr=0,spin=0,fx=fx,fy=fy,fz=fz,fr=fs2)
					#print 'e2 = '+str(e)+' fx = '+str(fx)+' fy = '+str(fy)+' fz = '+str(fz)	
	E_Pauli = (E_1 + E_2)  
	return E_Pauli 

if __name__ == "__main__":

	def E_tot(p_cfg):
		[calc,Z,R,r,s,delta] = read_cfg(p_cfg)
		#
		print '\nEke'
		Eke = E_ke(calc,s)
		calc.show_all()
		#
		print '\nEnucnuc'
		Enucnuc = E_nuc_nuc(calc,Z,R)
		calc.show_all()
		#
		print '\nEnucelec'
		Enucelec = E_nuc_elec(calc,Z,R,r,s)	
		calc.show_all()
		#
		print '\nEelecelec'
		Eelecelec = E_elec_elec(calc,r,s)
		calc.show_all()
		#
		print '\nEPauli'
	        EPauli = E_Pauli(calc,delta,s,r)
	       	calc.show_all()
		#
		print '\n'	
		print '--------------------'
		print 'Energy contributions'
		print '--------------------'
		print 'Eke      =\t', Eke 
		print 'Enucnuc  =\t', Enucnuc
		print 'Enucelec =\t', Enucelec
		print 'Eelecelc =\t',Eelecelec
		print 'EPauli   =\t', EPauli
		print 'Etot    =\t', calc.total_energy()

	def main():
		# 
		# pyeff calculation 
		#
		p_cfg = '../tests/CH4/structs/CH4.cfg'
		E_tot(p_cfg)
		#
		# reference calculation eff-code  
		#
		print '\n'
	        print '--------------------'
	        print ' Reference system   '
	        print '--------------------'
		# reference calculation  
		p_eff = '../tests/CH4/reference/start/CH4.eff'	
		calc_ref = read_eff(p_eff,print_data=None)
		calc = calc_ref
		calc.show_all()
		#
	        print '\n'
	        print '--------------------'
	        print ' Reference Energy   '
	        print '--------------------'
	
		print 'Etot    =\t', calc.total_energy()	
		
	main()
