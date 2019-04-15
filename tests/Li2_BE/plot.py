from matplotlib.pyplot import * 

bohr = 0.5291772105638411

pyeff = np.loadtxt('BE.dat')

E_atom_pyeff =  -6.08004291199


fig1 = figure(1) 

plot(pyeff[:,0],pyeff[:,1],'o-',label=r'pyeff Li$_{2}$')
plot([pyeff[1,0],pyeff[-1,0]],[2*E_atom_pyeff,2*E_atom_pyeff],label='pyeff 2 x Li atom')
xlabel('Distance r$_{Li-Li}$ [$\AA$]')
ylabel('Total energy [Hartree]')
legend()
show()

