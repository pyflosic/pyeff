#   Copyright 2019 PyEFF developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
from ase.io import read 
from ase.atoms import atomic_numbers 
try: 
    from ase.atoms import atomic_numbers
except:
    from ase.data import atomic_numbers
from ase.units import Bohr 
import numpy as np
def xyz_to_nuclei_fod(ase_atoms):
	# this function is from the pyflosic flosic.py routine 
        # get the number of atoms + fod
        syssize = len(ase_atoms)
        # Determine, if fod is included or not.
        # A fod is included, if it is no closer then fodeps to another fod.
        fodeps = 0.1
        included = []
        nombernotincluded = 0
        dist = -1.0
        for i in range(0,syssize):
                # Cores are always included.
                if ase_atoms[i].symbol != 'X' and ase_atoms[i].symbol != 'He':
                        included.append(1)
                if ase_atoms[i].symbol == 'X' or ase_atoms[i].symbol == 'He':
                        distold = fodeps
                        for j in range(0,i):
                                dist = np.sqrt((ase_atoms[i].position[0]-ase_atoms[j].position[0])**2+(ase_atoms[i].position[1]-ase_atoms[j].position[1])**2+(ase_atoms[i].position[2]-ase_atoms[j].position[2])**2)
                                if dist < distold and included[j]==1:
                                        distold = dist
                        if distold < fodeps:
                                included.append(0)
                                nombernotincluded += 1
                        else:
                                included.append(1)
        # process fods and cores separetely.
        tmp = []
        nuclei = Atoms([])
        fod1 = Atoms([])
        fod2 = Atoms([])
        nrofnuclei = 0
        for i in range(0,syssize):
                if ase_atoms[i].symbol != 'X':
                        nrofnuclei = nrofnuclei + 1
                if ase_atoms[i].symbol == 'X':
                        break
        for i in range(0,syssize):
                if i < nrofnuclei:
                        tmp.append(ase_atoms[i].symbol+' '+str(ase_atoms[i].position[0])+' '+str(ase_atoms[i].position[1])+' '+str(ase_atoms[i].position[2])+'\n')
                        nuclei.append(ase_atoms[i])
                elif ase_atoms[i].symbol == 'X':
                        fod1.append(ase_atoms[i])
                        if included[i] == 1:
                                tmp.append('ghost1'+' '+str(ase_atoms[i].position[0])+' '+str(ase_atoms[i].position[1])+' '+str(ase_atoms[i].position[2])+'\n')
                elif ase_atoms[i].symbol == 'He':
                        fod2.append(ase_atoms[i])
                        if included[i] == 1:
                                tmp.append('ghost2'+' '+str(ase_atoms[i].position[0])+' '+str(ase_atoms[i].position[1])+' '+str(ase_atoms[i].position[2])+'\n')
        geo = ''.join(tmp)
        return geo,nuclei,fod1,fod2,included


def xyz2cfg(p_xyz,er):
        # needs xyz containing nuclei and fod positions 
	# xyz units are Angstroem 
        # EFF units are Bohr
	# read
	ase_atoms = read(p_xyz)  
	[geo,nuclei,fod1,fod2,included] = xyz_to_nuclei_fod(ase_atoms) 

        # writing 
        o = open('pyflosic.cfg','w')
        o.write('@params\n')
        o.write('calc = minimize\n')
        o.write('@nuclei\n')
	nuc_sym = nuclei.get_chemical_symbols()
	nuc_pos = nuclei.get_positions()/Bohr
	for n in range(len(nuc_sym)):
                o.write(str(nuc_pos[n][0])+' '+str(nuc_pos[n][1])+' '+str(nuc_pos[n][2])+' '+str(atomic_numbers[nuc_sym[n]])+'\n')
        o.write('@electrons\n')
        fod1_sym = fod1.get_chemical_symbols()
        fod1_pos = fod1.get_positions()/Bohr
	for e in range(len(fod1_sym)):
                # 1st spin channel gets pyeff spin 1 
                o.write(str(fod1_pos[e][0])+' '+str(fod1_pos[n][1])+' '+str(fod1_pos[n][2])+' 1 '+str(er)+'\n')
	fod2_sym = fod2.get_chemical_symbols()
        fod2_pos = fod2.get_positions()/Bohr
        for e in range(len(fod2_sym)):
                # 2nd spin channel gets pyeff spin -1 
        	o.write(str(fod2_pos[e][0])+' '+str(fod2_pos[n][1])+' '+str(fod2_pos[n][2])+' -1 '+str(er)+'\n')
        o.close()


if __name__ == "__main__":
	xyz2cfg(p_xyz='B.xyz',er=0.7)

