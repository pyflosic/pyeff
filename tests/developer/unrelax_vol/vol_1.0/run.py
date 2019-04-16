from pyeff_system import *
from pyeff_energy_force import * 
from pyeff_ewald_system import *
from pyeff_ewald_energy_force import * 

def E_tot(p_cfg):
        [calc,Z,R,r,s,delta] = read_cfg(p_cfg)
        #
        #print '\nEke'
        Eke = E_ke(calc,s)
        #calc.show_all()
        #
        #print '\nEnucnuc'
        Enucnuc = E_nuc_nuc(calc,Z,R)
        #calc.show_all()
        #
        #print '\nEnucelec'
        Enucelec = E_nuc_elec(calc,Z,R,r,s)
        #calc.show_all()
        #
        #print '\nEelecelec'
        Eelecelec = E_elec_elec(calc,r,s)
        #calc.show_all()
        #
        #print '\nEPauli'
        EPauli = E_Pauli(calc,delta,s,r)
        calc.show_all()
        #
        print '\n'
        print '--------------------'
        print 'Energy contributions'
        print '--------------------'
        print 'Eke      =\t', Eke
        print 'Enucnuc  =\t', Enucnuc
        print 'Enucelec =\t', Enucelec
        print 'Eelecelc =\t',Eelecelec
        print 'EPauli   =\t', EPauli
        print 'Etot    =\t', calc.total_energy()
	return calc 

def main():
        # functionality test
	f_cfg = 'li_solid_111.cfg' 
       	[calc,Z,R,r,s,delta] = read_cfg(f_cfg) 
	a = 8.352000
	calc = E_tot(f_cfg)
	calc.show_all()
        sys_ewald = ewald_system(Lx=a,Ly=a,Lz=a,system=calc)
        print 'Total volume %0.5f' % sys_ewald.get_volume()
        print ' okay '
	print 'Total charge %i' % sys_ewald.get_total_charge()
        print 'Ewald r cutoff %0.9f' % sys_ewald.get_ewald_r_cutoff()
        print ' okay '
	print 'Ewald k cutoff %0.9f' % sys_ewald.get_ewald_k_cutoff()
	print ' okay '
	rgrid = sys_ewald.make_r_space()
	np.savetxt('rgrid.dat',np.transpose(rgrid))
        print 'Number of real space grid points %i' % sys_ewald.num_gps
	print 'Number of real space grid points %i' % len(sys_ewald.rgrid_x)
	kgrid = sys_ewald.make_k_space()
	np.savetxt('kgrid.dat',np.transpose(kgrid))
	print 'Number of k space grid points %i' % sys_ewald.num_kps
	for c in range(sys_ewald.numcharges):
		print(sys_ewald.ewald[c].q)
	#sys_ewald.num_kps = 1 
	#sys_ewald.kgrid_x[0] = -0.25 
	#sys_ewald.kgrid_y[0] = 0
	#sys_ewald.kgrid_z[0] = 0 
	#print sys_ewald.kgrid_x[0],sys_ewald.kgrid_y[0],sys_ewald.kgrid_z[0]	
	print '----------------------'
	print ' Clear molecular data '
	print '----------------------'
        for n in range(len(sys_ewald.system.nuc)):
                sys_ewald.system.nuc[n].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*sys_ewald.system.nuc[n].energy,fx=-1*sys_ewald.system.nuc[n].fx,fy=-1*sys_ewald.system.nuc[n].fy,fz=-1*sys_ewald.system.nuc[n].fz,fr=-1*sys_ewald.system.nuc[n].fr)
        for e in range(len(sys_ewald.system.elec)):
                sys_ewald.system.elec[e].update(px=0,py=0,pz=0,pr=0,spin=0,energy=-1*sys_ewald.system.elec[e].energy,fx=-1*sys_ewald.system.elec[e].fx,fy=-1*sys_ewald.system.elec[e].fy,fz=-1*sys_ewald.system.elec[e].fz,fr=-1*sys_ewald.system.elec[e].fr)
	print '-------'
	print ' E_kin '
	print '-------'
	print ' okay ' 
        Eke = E_ke(sys_ewald.system,s)
	sys_ewald.system.show_all()
	print '---------------'
        print 'E_PauliPeriodic'
        print '---------------'
        print ' okay          '
	E_PauliPeriodic(sys_ewald)
        sys_ewald.system.show_all()
	print '--------'
        print ' E_self '
        print '--------'
        print ' okay   '
        print E_self(sys_ewald)
        sys_ewald.show_all()
	print '--------'
	print 'E_kspace'
	print '--------' 
	print ' okay   '
	print E_kspace(sys_ewald)
	sys_ewald.show_all()
	print '--------'
	print 'E_rspace'
	print '--------'
	print ' okay   '
        print E_rspace(sys_ewald)
        sys_ewald.show_all()
	print '---------'
	print 'E_uniform'
	print '---------'
	print E_uniform(sys_ewald)
	sys_ewald.show_all()
	#sys_ewald.add_ewald_energy_force()
	#sys_ewald.system.show_all()
	#
	#
	#
	sys_ewald.add_ewald_energy_force()
	sys_ewald.system.show_all()
	sys_ewald.system.total_energy()	
main()

