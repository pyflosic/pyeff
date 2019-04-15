# Calculate binding energy (BE) 
# based on the optimized values for all atoms and the system 

E_atom1_pyeff = -0.424375627902 
E_system_pyeff = -0.955791165197

BE = (2*E_atom1_pyeff)-E_system_pyeff 

# Hartree to kcal/mol 
f_Ha2kcalmol = 630 
print 'Binding energy = '+str(BE*f_Ha2kcalmol)+' kcal/mol'
