import glob 
def get_values(p_file):
	f = open(p_file,'r') 
	ll = f.readlines() 
	f.close() 
	# tags 
	tag_vol = 'Total volume'
	tag_ene = 'Etot' 
	E = []
	V = []
	for l in range(len(ll)):
		if ll[l].find(tag_vol) != -1:
			v = float(ll[l].split()[-1])
			V.append(v) 
		if ll[l].find(tag_ene) != -1:
                	e = float(ll[l].split()[-1])
			E.append(e)
	return V[-1],E[-1] 

def get_all_values():
	o = open('pyeff_vol_ene.dat','w')
	dirs = glob.glob('./vol_*/os_opt/pyeff.log')
	print dirs
	for d in dirs:
		relvol =  float(d.split('/')[1].split('_')[-1])
		[vol,ene] = get_values(p_file=d)
		o.write('%0.9f %0.9f %0.9f\n' % (relvol,vol,ene))
	o.close()

get_all_values()
