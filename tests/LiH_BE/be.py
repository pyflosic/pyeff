# Calculate binding energy (BE) 
# based on the optimized values for all atoms and the system 

E_atom1_pyeff = -6.08004291199
E_atom2_pyeff = -0.424375627902
E_system_pyeff = -6.59337109926

BE = (E_atom1_pyeff+E_atom2_pyeff)-E_system_pyeff 
# Hartree to kcal/mol 
f_Ha2kcalmol = 630 
print 'Binding energy = '+str(BE*f_Ha2kcalmol)+' kcal/mol'
