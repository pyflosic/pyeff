import numpy as np 
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io.cube import * 
from ase.atoms import Atoms

# we need the masses 
# we need a grid 

def psi(x,px,r,s,ps):
	# x 	...	position of the Gaussian wave packet 
	# px 	... 	momenta of the electron 
	# r 	... 	position in space
	# s	...	size of the electron 
	# ps	...	momenta of the size of the electron/ radial momenta
	psi = np.exp(-1*(np.complex((1/s**2),0)-np.complex(0,(2*ps/2)))*(np.linalg.norm(r-x))**2)*np.exp(np.complex(0,np.dot(px,x)))
	return psi 

def momenta(x,s,m_elec=1):
	px = m_elec*x
	ps = 3*m_elec/4*s
	return px,ps

def psi2grid(grid,x,s,cell=[5,5,5]):
	# grid	...	triple grid=[gp_x,gp_y,gp_z]
	# move to the center 
	disp = np.array([grid[0]/2.,grid[1]/2.,grid[2]/2.])
	x_init = x 
	x = x + disp
	print x 
	# grid properties 
	# gp	...	grid points 
	gp_x= grid[0]
	gp_y= grid[1]
	gp_z= grid[2]
	grid = np.zeros([gp_x,gp_y,gp_z])
	# write cube like values 
	for idx_z in range(gp_z):
		for idx_y in range(gp_y):
			for idx_x in range(gp_x):
				r = np.array([idx_x,idx_y,idx_z])
				px,ps = momenta(x,s)
				grid[idx_x,idx_y,idx_z] = psi(x,px,r,s,ps).real
	# gs 	...	grid spacing 
	gs_x = cell[0]/float(gp_x)
	gs_y = cell[1]/float(gp_y)
	gs_z = cell[2]/float(gp_z)
	x = np.array([x[0]*gs_x,x[1]*gs_y,x[2]*gs_z])
	struct = Atoms('H',cell=cell,positions=[x])
	#f = open('test.cube','w')
        #write_cube(fileobj=f,atoms=struct,data=grid)
	#f.close()
	return x,grid 

def get_cube(struct,grid):
	f = open('test.cube','w')
        write_cube(fileobj=f,atoms=struct,data=grid)
	f.close()


def main():
	grid = [40,40,40]
	# Atom1 
	x = np.array([-4,0,0])
	s = 3.2 
	x1,grid1 = psi2grid(grid,x,s)
	# Atom2 
	x = np.array([8,0,0])
        s = 3.2
        x2,grid2 = psi2grid(grid,x,s)
	# 
	struct = Atoms('XX',cell=[5,5,5],positions=[x1,x2])	
	get_cube(struct,grid1+grid2)
main()	
