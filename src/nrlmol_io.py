# author: S. Schwalbe 
# usage: functions to convert nrlmol input to pyeff input and reverse 

# chemical symbols used in nrlmol 
nrlmol_symbols =['HYD','HEL','LIT','BER','BOR','CAR','NIT','OXY',
		'FLU','NEO','SOD','MAG','ALU','SIL','PHO','SUL',
		'CHL','ARG','POT','CAL','SCA','TIT','VAN','CHR',
		'MAN','IRO','COB','NIC','COP','ZIN','GAL','GER',
		'ARS','SEL','BRO','KRY','RUB','STR','YTR','ZIR',
		'NIO','MOL','TEC','RHU','RHO','PAL','SLV','CAD',
		'IND','TIN','ANT','TEL','IOD','XEN','CES','BAR']

def get_element_numbers(e):
	# converts nrlmol.chemical_symbols to element numbers 
	return nrlmol_symbols.index(e)+1

def nrlmol2cfg(p_symbol,p_frmorb,er):
	# needs SYMBOLS and FRMORB file 
	# SYMBOLS/FRMORB untis are atomic units (Bohr)
	# EFF units are also Bohr
	# search tags 
	tag_nuc_start ='1  CALCULATION SETUP BY ISETUP'
	tag_nuc_end = 'ELECTRONS  = '
	f = open(p_symbol,'r')
	ll = f.readlines()
	# symbol reading 
	for l in range(len(ll)):
		if ll[l].find(tag_nuc_start) != -1:
			idx_nuc_start = l + 1 
		if ll[l].find(tag_nuc_end) != -1:
			idx_nuc_end = l 
	f2 = open(p_frmorb,'r')
	ll2 = f2.readlines()
	num_elec = ll2[0].split()
	# writing 
	o = open('nrlmol.cfg','w')
        o.write('@params\n')
        o.write('calc = minimize\n')
        o.write('@nuclei\n')
       	for n in range(idx_nuc_start,idx_nuc_end):
                tmp = ll[n].split()
                sym =  str('%.3s' %(tmp[0].split('-')[1]))
                o.write(str(tmp[2])+' '+str(tmp[3])+' '+str(tmp[4])+' '+str(get_element_numbers(sym))+'\n') 
	o.write('@electrons\n')
	for e in range(1,int(num_elec[0])+1):
		tmp = ll2[e].split()
		# 1st spin channel gets pyeff spin 1 
		o.write(str(tmp[0])+' '+str(tmp[1])+' '+str(tmp[2])+' 1 '+str(er)+'\n')
	for e in range(int(num_elec[0])+1,int(num_elec[0])+int(num_elec[1])+1):
                tmp = ll2[e].split()
		# 2nd spin channel gets pyeff spin -1 
                o.write(str(tmp[0])+' '+str(tmp[1])+' '+str(tmp[2])+' -1 '+str(er)+'\n')
	f.close()
	o.close()


if __name__ == "__main__":

	nrlmol2cfg(p_symbol='./SYMBOL',p_frmorb='./FRMORB',er=0.7)
