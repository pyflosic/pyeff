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
from sympy import *
from sympy import print_python

# eff1 analytic/numerical force derivations
#
# usage: for debugging errors 
#	 delivers anyltic and numerical derivatives 
#  
# author: S. Schwalbe
#
# variables: 
#		tag ana -- analytic 
#		tag num -- numerical  
#
# NOTE:	may only use H2 test for debugging 
# 	because analytic expression are evaluated 
#	at each iteration step 
# Eke 
#
def sympy_Eke(value_s,debug=None):
	# analytic 
	s = symbols('s') 
	Eke = Function('Eke')(s)
	dEkeds = Function('dEkeds')(s)
	E = Function('E')(s)
	fs = Function('fs')(s)
	Eke = (3./2.)*(1./s**2)
	dEkeds = -1*diff(Eke,s)
	Eke = lambdify(s,Eke,'numpy')
	dEkeds = lambdify(s,dEkeds,'numpy')
	# numerical 
	def fun(s):
		return Eke(s)

	def d_fun(s):
	    	# differential variable 
	    	x=s
	    	h = 1e-5
	    	return (fun(x+h)-fun(x-h))/(2*h)

	ana_Eke = Eke(value_s)	
	ana_dEkeds = dEkeds(value_s)
	num_Eke = fun(value_s)
	num_dEkeds = -1*d_fun(value_s)

	if debug != None:	
		print 'ana_Eke = '+str(ana_Eke)
		print 'ana_dEkeds = '+str(ana_dEkeds)
        	print 'num_Eke = '+str(num_Eke)
        	print 'num_dEkeds = '+str(num_dEkeds)

	return ana_Eke,ana_dEkeds,num_Eke,num_dEkeds

sympy_Eke(2.0)

#
# E_nucnuc
#
def sympy_E_nucnuc(value_Zi,value_Zj,value_R,debug=None): 
	# analytic
	Zi = symbols('Zi')
	Zj = symbols('Zj')
	R  = symbols('R')
	ENN = Function('ENN')(Zi,Zj,R)
	dENNdR = Function('dENNdR')(Zi,Zj,R)
	ENN = Zi*Zj/R
	dENNdR = -1*diff(ENN,R)
	ENN = lambdify([Zi,Zj,R],ENN,'numpy')
	dENNdR = lambdify([Zi,Zj,R],dENNdR,'numpy')
	# numerical
	def fun(Zi,Zj,R):
                return ENN(Zi,Zj,R)

        def d_fun(Zi,Zj,R):
                # differential variable 
                x=R
                h = 1e-5
                return (fun(Zi,Zj,x+h)-fun(Zi,Zj,x-h))/(2*h)
	
	ana_ENN = ENN(value_Zi,value_Zj,value_R)
	ana_dENNdR  = dENNdR(value_Zi,value_Zj,value_R)
	num_ENN = fun(value_Zi,value_Zj,value_R)
	num_dENNdR = -1*d_fun(value_Zi,value_Zj,value_R)

	if debug != None:
	        print 'ana_ENN = '+str(ana_ENN)
	        print 'ana_dENNdR = '+str(ana_dENNdR)
		print 'num_ENN = '+str(num_ENN )
	        print 'num_dENNdR = '+str(num_dENNdR)

	return ana_ENN,ana_dENNdR,num_ENN,num_dENNdR

#
# E_nucelec
# 
def sympy_E_nucelec(value_Zi,value_R,value_sj,debug=None):
	# analytic 
	Zi = symbols('Zi')
	R = symbols('R')
	sj = symbols('sj')
	ENE = Function('ENE')(Zi,R,sj)
	dENEdR = Function('dENEdR')(Zi,R,sj)
	#
	ENE = -1*Zi/R*erf((sqrt(2)*R/sj))
	dENEdR = -1*diff(ENE,R)
	dENEds = -1*diff(ENE,sj)
	#
	ENE = lambdify([Zi,R,sj],ENE,modules=['numpy','sympy'])
        dENEdR = lambdify([Zi,R,sj],dENEdR,modules=['numpy','sympy'])
	# numerical 
	def fun(Zi,R,sj):
                return ENE(Zi,R,sj)

        def d_fun(Zi,R,sj):
                # differential variable 
                x=R
                h = 1e-5
                return (fun(Zi,x+h,sj)-fun(Zi,x-h,sj))/(2*h)
	
	ana_ENE = ENE(value_Zi,value_R,value_sj)
	ana_dENEdR = dENEdR(value_Zi,value_R,value_sj)
	num_ENE = fun(value_Zi,value_R,value_sj)
	num_dENEdR = -1*d_fun(value_Zi,value_R,value_sj)
	if debug != None:
		print 'ana_ENE = '+str(ana_ENE)
        	print 'ana_dENEdR = '+str(ana_dENEdR)
        	print 'num_ENE = '+str(num_ENE)
        	print 'num_dENEdR = '+str(num_dENEdR)
	return ana_ENE,ana_dENEdR,num_ENE,num_dENEdR

# 
# E_elecelec
# 
def sympy_E_elecelec(value_r,value_si,value_sj,debug=None):
        r = symbols('r')
        si = symbols('si')
	sj = symbols('sj')
	EEE = Function('EEE')(r,si,sj)
	dEEEdr = Function('dEEEdr')(r,si,sj)
        EEE = 1/r*erf((sqrt(2)*r/(sqrt(si**2+sj**2))))
        # dE/dR
        dEEEdr = simplify(-1*diff(EEE,r))
        # dE/ds1
        dEEEds1 = simplify(-1*diff(EEE,si))
        # dE/ds2
        dEEEds2 = simplify(-1*diff(EEE,sj))
	# 
	EEE = lambdify([r,si,sj],EEE,modules=['numpy','sympy'])
        dEEEdr = lambdify([r,si,sj],dEEEdr,modules=['numpy','sympy'])
	#
	# numerical 
        def fun(r,si,sj):
                return EEE(r,si,sj)

        def d_fun(r,si,sj):
                # differential variable 
                x=r
                h = 1e-5
                return (fun(x+h,si,sj)-fun(x-h,si,sj))/(2*h)

	ana_EEE = EEE(value_r,value_si,value_sj)
        ana_dEEEdr = dEEEdr(value_r,value_si,value_sj)
        num_EEE = fun(value_r,value_si,value_sj)
        num_dEEEdr = -1*d_fun(value_r,value_si,value_sj)
        if debug != None:
                print 'ana_ENE = '+str(ana_EEE)
                print 'ana_dENEdR = '+str(ana_dEEEdr)
                print 'num_ENE = '+str(num_EEE)
                print 'num_dENEdR = '+str(num_dEEEdr)
        return ana_EEE,ana_dEEEdr,num_EEE,num_dEEEdr

def sympy_E_Pauli(value_r,value_si,value_sj,samespin,debug=None):
	# analytic
	print value_r 
	print value_si 
	print value_sj
	value_rho = -0.2
	value_r = 1.125 * value_r 
	value_si = 0.9 * value_si
	value_sj = 0.9 * value_sj
	# 
	rho = symbols('rho')
	r = symbols('r')
        si = symbols('si')
        sj = symbols('sj')
	# declarations 
 	Delta_T = Function('Delta_T')(r,si,sj)
        S = Function('S')(r,si,sj)
	E_up_up = Function('E_up_up')(rho,S)
	E_up_down = Function('E_up_down')(rho,S)
	# defintions 
	Delta_T = (3./2.)*((1./(si**2))+(1./(sj**2)))-2.*((3.*((si**2)+(sj)**2))-2.*(r**2))/(((si**2)+(sj**2))**2)
	S = ((2./((si/sj)+(sj/si)))**(3./2.))*exp((-1.*(r**2))/(si**2+sj**2))
	# derivatives 
	# note: we use here not forces; positive derivatives 
	dDelta_Tdr = simplify(diff(Delta_T,r))
	dSdr = simplify(diff(S,r))
	if samespin == True:
		# E 
		E_up_up = (S**2/(1.0-(S**2))+(1.-rho)*S**2/(1.0+(S**2)))*Delta_T
		print 'E_up_up = ', E_up_up
		# dE/dr 
		dE_up_updr = simplify(diff(E_up_up,r))
		print 'dE_up_updr = ', dE_up_updr
		dE_up_updr = lambdify([rho,r,si,sj],dE_up_updr,modules=['numpy','sympy'])
                dE_up_updr = dE_up_updr(value_rho,value_r,value_si,value_sj)
		# -1 = force factor, 1.125 = d_overline(r)/dr = 1.125 (iternal derivative) 
		ana_fr = -1*1.125*dE_up_updr
		# dE/ds1 
                dE_up_updsi = simplify(diff(E_up_up,si))
                print 'dE_up_updsi = ',dE_up_updsi
                dE_up_updsi = lambdify([rho,r,si,sj],dE_up_updsi,modules=['numpy','sympy'])
                dE_up_updsi = dE_up_updsi(value_rho,value_r,value_si,value_sj)
                ana_fs1 = -1*0.9*dE_up_updsi
		# dE/ds2
                dE_up_updsj = simplify(diff(E_up_up,sj))
                print 'dE_up_updsj = ',dE_up_updsj
                dE_up_updsj = lambdify([rho,r,si,sj],dE_up_updsj,modules=['numpy','sympy'])
                dE_up_updsj = dE_up_updsj(value_rho,value_r,value_si,value_sj)
                ana_fs2 = -1*0.9*dE_up_updsj
		# eval(E_up_up)
		E_up_up = lambdify([rho,r,si,sj],E_up_up,modules=['numpy','sympy'])
                E_up_up = E_up_up(value_rho,value_r,value_si,value_sj)
                ana_E = E_up_up
		# eval(S) 
                S = lambdify([r,si,sj],S,modules=['numpy'])
                S = S(value_r,value_si,value_sj)
                ana_S = S
                # eval(Delta_T)
                Delta_T = lambdify([r,si,sj],Delta_T,modules=['numpy'])
                Delta_T = Delta_T(value_r,value_si,value_sj)
                ana_Delta_T = Delta_T

	if samespin == False:
		# E  
		E_up_down = -1*rho*S**2/(1.0+(S**2))*Delta_T
		#print E_up_down
		# dE/dr 
		dE_up_downdr = simplify(diff(E_up_down,r))
		#print dE_up_downdr
		dE_up_downdr = lambdify([rho,r,si,sj],dE_up_downdr,modules=['numpy','sympy'])
                dE_up_downdr = dE_up_downdr(value_rho,value_r,value_si,value_sj)
		# -1 = force factor, 1.125 = d_overline(r)/dr = 1.125 (iternal derivative) 
		ana_fr = -1*1.125*dE_up_downdr
		# dE/ds1 
		dE_up_downdsi = simplify(diff(E_up_down,si))
		print dE_up_downdsi 
		dE_up_downdsi = lambdify([rho,r,si,sj],dE_up_downdsi,modules=['numpy','sympy'])
                dE_up_downdsi = dE_up_downdsi(value_rho,value_r,value_si,value_sj)
		ana_fs1 = -1*0.9*dE_up_downdsi
		# dE/ds2 
                dE_up_downdsj = simplify(diff(E_up_down,sj))
                print dE_up_downdsj
                dE_up_downdsj = lambdify([rho,r,si,sj],dE_up_downdsj,modules=['numpy','sympy'])
                dE_up_downdsj = dE_up_downdsj(value_rho,value_r,value_si,value_sj)
                ana_fs2 = -1*0.9*dE_up_downdsj
		#
		E_up_down = lambdify([rho,r,si,sj],E_up_down,modules=['numpy','sympy'])
                E_up_down = E_up_down(value_rho,value_r,value_si,value_sj)
                ana_E = E_up_down
		#
		S = lambdify([r,si,sj],S,modules=['numpy'])
		S = S(value_r,value_si,value_sj)
                ana_S = S
		# 
                Delta_T = lambdify([r,si,sj],Delta_T,modules=['numpy'])
                Delta_T = Delta_T(value_r,value_si,value_sj)
                ana_Delta_T = Delta_T


	if debug != None:
                print 'ana_S = '+str(ana_S)
		print 'ana_Delta_T = '+str(ana_Delta_T)
		print 'ana_E_Pauli = '+str(ana_E)
                print 'ana_d_Paulidr = '+str(ana_fr)
		print 'ana_d_Paulids1 = '+str(ana_fs1)
        return ana_E,ana_fr,ana_fs1,ana_fs2


if __name__ == "__main__":
	def main():
		# functionality test 
		# values from H2 example 
		print sympy_Eke(2.0)	
		print sympy_E_nucnuc(1,1,1.4)
		print sympy_E_nucelec(1,1.4,2.0)
	main()

