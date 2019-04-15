from matplotlib.pyplot import * 

bohr = 0.5291772105638411

ref = np.loadtxt('./reference/ref.dat')
pyeff = np.loadtxt('BE.dat')

E_H_pyeff = -0.424343544816

fig1 = figure(1) 

plot(pyeff[:,0],pyeff[:,1],label=r'pyeff H$_{2}$')
plot([pyeff[1,0],pyeff[-1,0]],[2*E_H_pyeff,2*E_H_pyeff],label='pyeff 2 x H atom')
xlabel('Distance r$_{H-H}$ [$\AA$]')
ylabel('Total energy [Hartree]')
legend()
show()

fig2 = figure(2)

plot(ref[:,0]*bohr,ref[:,1]-ref[-1,1],label=r'FCI H$_{2}$')
plot(pyeff[:,0],pyeff[:,1]-pyeff[-1,1],label=r'pyeff H$_{2}$')
xlabel('Distance r$_{H-H}$ [$\AA$]')
ylabel('$E_{tot}-E_{\infty}$ [Hartree]')
legend()
show()

