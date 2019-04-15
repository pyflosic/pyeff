from ase.io import read
from ase.ga.utilities import get_rdf,get_nnmat
from matplotlib.pyplot import *
from ase.atoms import Atoms
import numpy as np
import glob

l = 5.0
nmax = 100

atoms = read(glob.glob('final.xyz')[0])
atoms.set_cell(np.eye(3)*l)
print atoms
up = Atoms()
up.set_cell(atoms.get_cell())
up_target = 'X'
dn = Atoms()
dn.set_cell(atoms.get_cell())
dn_target = 'He'
nuclei = Atoms()
nuclei.set_cell(atoms.get_cell())


for a in range(len(atoms)):
        if atoms[a].symbol == up_target:
                up.append(atoms[a])
        if atoms[a].symbol == dn_target:
                dn.append(atoms[a])
        if atoms[a].symbol != up_target and atoms[a].symbol != dn_target:
                nuclei.append(atoms[a])

data_up = get_rdf(up,l,nmax)
data_dn = get_rdf(dn,l,nmax)
data_nuclei  = get_rdf(nuclei,l,nmax)
data_all = get_rdf(atoms,l,nmax)
# data_up[0]/max(data_up[0])*len(up)
# data_dn[0]/max(data_dn[0])*len(dn)
bar(data_up[1],data_up[0],label='up',width=0.02)
bar(data_dn[1],data_dn[0],label='dn',width=0.01)
bar(data_nuclei[1],data_nuclei[0],label='nuclei',width=0.02)
#bar(data_all[1],data_all[0],label='all',width=0.05)
legend()
show()
