from matplotlib.pyplot import * 

bohr = 0.5291772105638411

pyeff = np.loadtxt('BE.dat')

E_atom1_pyeff = -6.08004291199
E_atom2_pyeff = -0.424375627902


fig1 = figure(1) 

plot(pyeff[:,0],pyeff[:,1],'o-',label=r'pyeff LiH')
plot([pyeff[0,0],pyeff[-1,0]],[E_atom1_pyeff+E_atom2_pyeff,E_atom1_pyeff+E_atom2_pyeff],label='pyeff Li atom + H atom')
xlabel('Distance r$_{Li-H}$ [$\AA$]')
ylabel('Total energy [Hartree]')
legend()
show()

fig2 = figure(2)

flosic = np.loadtxt('./reference/FLOSIC.dat')
plot(pyeff[:,0],pyeff[:,1]-pyeff[:,1].min(),'o-',label=r'pyeff LiH')
plot(flosic[:,0],flosic[:,1]-flosic[:,1].min(),'o-',label=r'FLOSIC LiH')
xlabel('Distance r$_{Li-H}$ [$\AA$]')
ylabel('Total energy [Hartree]')
legend()
show()

