# atoms, molecules, cluster (no pbc)
| file 	| description |  
| ------------- |:-------------:|
| pyeff_system.py | defines the mother class particle and the daugther class system | 
| pyeff_energy_force.py | analytical definitions for energies and forces (EFF1) | 
| pyeff_io.py | write cfg files | 
| pyeff_calculator.py | ase pyeff calculator | 
| pyeff_optimizer.py | ase pyeff optimizer | 
| pyeff_sympy_functions.py | using sympy for checking derivatives | 
| pyeff_wavefunction.py | plotting wavefunctions | 
| nrlmol_io.py | convert NRLMOL input (SYMBOL, FRMORB) files to pyeff input files (cfg) | 
| pyflosic_io.py | convert PyFLOSIC (xyz) files to pyeff input files (cfg) | 

# periodic boundary conditions (pbcs)
| file 	| description |  
| ------------- |:-------------:|
| pyeff_ewald_system.py | similar to pyeff_system.py new mother class ewald_charge and new daugther class ewald_system | 
| pyeff_ewald_energy_force.py | similar to pyeff_energy_force.py but defining energy terms needed for pbcs | 
| pyeff_calculator_pbc.py | ase calculator for pbc calculations | 
| pyeff_optimizer_pbc.py | ase optimizer for pbc calculations | 

# development and unfinished 
| file 	| description |  
| ------------- |:-------------:|
| pyeff_ewald_system_numba.py | try to speed up  pyeff_ewald_system.py with numba| 
| pyeff_ewald_energy_force_numba.py | | try to speed up pyeff_ewald_energy_force.py  with numba | 

