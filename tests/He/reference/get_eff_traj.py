import io 

def get_eff_traj(p_eff):
	f = open(p_eff,'r')
	o = open('eff_traj.xyz','w')
	ll = f.readlines()
	# tags 
	# chemical symbols list 
	# atoms and electrons  
	atoms_tag = '[nuc]' 
	elec_tag ='[elec]'
	# get start and end of each snapshot 
	start_atoms_pos_tag = '[position_nuc_header]'
	end_atoms_pos_tag = '[position_elec_header]'	
	start_elec_pos_tag = '[position_elec_header]' 
	end_elec_pos_tag = '[velocity_nuc_header]'	
	sym_atoms = {}
	sym_elec = {}
	idx_start_atoms = []
	idx_end_atoms = []
	idx_start_elec = []
	idx_end_elec = []
	for l in range(len(ll)):
		if ll[l].find(atoms_tag) != -1: 
			tmp = ll[l].split()
			sym_atoms[tmp[1]]= tmp[-1]
	        if  ll[l].find(elec_tag) != -1:
                        tmp = ll[l].split()
                        if tmp[-1] == '1':
                                tmp[-1] = 'X'
                        if tmp[-1] == '-1':
                                tmp[-1] = 'He'
                        sym_elec[tmp[1]]= tmp[-1]

		if ll[l].find(start_atoms_pos_tag) != -1:
			idx_start_atoms.append(l+1)
		#if ll[l].find(end_atoms_pos_tag) != -1:
                #        idx_end_atoms.append(l)
		if ll[l].find(start_elec_pos_tag) != -1:
                        idx_start_elec.append(l+1)
                #if ll[l].find(end_elec_pos_tag) != -1:
                #        idx_end_elec.append(l)
	# c  ... configuration 
	# p  ... position atoms
	# pe ... position electrons 
	for c in range(len(idx_start_atoms)):
		o.write(str(len(sym_atoms)+(len(sym_elec)))+'\n')
		o.write('created by hand\n')
		for p in range(idx_start_atoms[c],idx_start_atoms[c]+len(sym_atoms)):		
			tmp =  ll[p].split()
			o.write(str(sym_atoms[tmp[1]])+' '+str(tmp[-3])+' '+str(tmp[-2])+' '+str(tmp[-1])+'\n')
		for pe in range(idx_start_elec[c],idx_start_elec[c]+len(sym_elec)):
                        tmp =  ll[pe].split()
                        o.write(str(sym_elec[tmp[1]])+' '+str(tmp[-4])+' '+str(tmp[-3])+' '+str(tmp[-2])+'\n')
		#o.write('\n')
	f.close()
	o.close()

#p_eff = 'C6H6_MD.eff'
p_eff = 'he.eff'
get_eff_traj(p_eff)
