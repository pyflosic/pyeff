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

