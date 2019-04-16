from matplotlib.pyplot import * 
import numpy as np

data = np.loadtxt('results.dat',skiprows=1)
fz = 25 
fig1 = figure(1) 
title('Li solid: single point pbc calculations',fontsize=fz)
ax = subplot(111)
plot(data[:,0],data[:,2],'rD',label='eff code',markersize=9)
plot(data[:,0],data[:,1],'go',label='PyEFF')
xlabel('Volume [a.u.]',fontsize=fz)
ylabel('Energy [a.u.]',fontsize=fz)
setp(ax.get_xticklabels(), fontsize=fz)
setp(ax.get_yticklabels(), fontsize=fz)
ax.legend(loc=2,prop={'size':fz-10})
show()
 
