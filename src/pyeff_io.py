import numpy as np 

def write_cfg(f_name,types,chemical_symbols,spins,x):
        # write a new cfg file based on the x vector 
        o = open(f_name,'w')
        o.write('@params\n')
        o.write('calc = minimize\n')
        o.write('@nuclei\n')
        idx = 0
        for p in range(len(types)):
                if types[p] == 'nuclei':
                        o.write(str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])+' '+str(chemical_symbols[p])+'\n')
                        idx = idx + 3
        o.write('@electrons\n')
        for p in range(len(types)):
                if types[p] == 'electron':
                        o.write(str(x[idx+0])+' '+str(x[idx+1])+' '+str(x[idx+2])+' '+str(spins[p])+' '+str(np.exp(x[idx+3]))+'\n')
                        idx = idx +4
        o.close()

