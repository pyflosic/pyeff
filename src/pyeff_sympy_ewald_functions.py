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

#
# E_r_space 
#
# E_r-space = \sum_R_{L} \sum_{i\neqj} E_{ij}(r_{ij}-R_{L})
# E_{ij} = 1/r_{ij} erf(\sqrt{\alpha_i*\alpha_j/(\alpha_i+\alpha_j)}r_{ij}
#        - 1/r_{ij} erf(\sqrt{\alpha_i*\alpha_j^{max}/(\alpha_i+\alpha_j^{max})}r_{ij}
# only for interaction between a charge an onther charge with \alpha > \alpha_{max}


def sympy_E_rspace(value_r,value_qi,value_qj,value_alphai,value_alphaj,value_alphamax,debug=None):
        # analytic 
        qi = symbols('qi')
	qj = symbols('qj')
        r = symbols('r')
        alphai = symbols('alphai')
	alphaj = symbols('alphaj')
	alphamax = symbols('alphamax')
	a1 = Function('a1')(alphai,alphaj)
	a2 = Function('a2')(alphai,alphamax)
        Erspace = Function('Erspace')(r,qi,qj,alphai,alphaj,alphamax)
	Erspace_red = Function('Erspace')(r,qi,qj,a1,a2)
        dErspacedr = Function('dENEdR')(r,qi,qj,alphai,alphaj,alphamax)
        #
        Erspace = qi*qj*(1/r*erf((sqrt((alphai*alphaj)/(alphai+alphaj))*r))-1/r*erf((sqrt((alphai*alphamax)/(alphai+alphamax))*r)))
	E1 = qi*qj*(1/r*erf((sqrt((alphai*alphaj)/(alphai+alphaj))*r)))
	E2 = qi*qj*(-1/r*erf((sqrt((alphai*alphamax)/(alphai+alphamax))*r)))
	Erspace_red = qi*qj*(1/r*erf(a1*r)-1/r*erf(a2*r))
	print 'Erspace = '
	print Erspace
        dErspacedr = -1*diff(Erspace,r)
	print 'dErspacedr = '
	print dErspacedr
        dErspacedalphai = -1*diff(Erspace,alphai)
	print 'dErspacedalphai = '
	print dErspacedalphai
	dE1dalphai = simplify(-1*diff(E1,alphai))
	print 'dE1dalphai = '
        print dE1dalphai
	dE2dalphai = simplify(-1*diff(E2,alphai))
        print 'dE2dalphai = '
        print dE2dalphai
	dErspacedalphaj = simplify(-1*diff(Erspace,alphaj))
	print 'dErspacedalphaj = '
	print dErspacedalphaj
	dE1dalphaj = simplify(-1*diff(E1,alphaj))
        print 'dE1dalphaj = '
        print dE1dalphaj
        dE2dalphaj = simplify(-1*diff(E2,alphaj))
        print 'dE2dalphaj = '
        print dE2dalphaj
	#
	#dErspace_reddalphai = simplify(-1*diff(Erspace_red,alphai))
        #print 'dErspace_reddalphai = '
        #print dErspace_reddalphai
        #dErspace_reddalphaj= simplify(-1*diff(Erspace_red,alphaj))
        #print 'dErspacedalphaj = '
        #print dErspace_reddalphaj

        #
        Erspace = lambdify([r,qi,qj,alphai,alphaj,alphamax],Erspace,modules=['numpy','sympy'])
        dErspacedr = lambdify([r,qi,qj,alphai,alphaj,alphamax],dErspacedr,modules=['numpy','sympy'])
	dErspacedalphai = lambdify([r,qi,qj,alphai,alphaj,alphamax],dErspacedalphai,modules=['numpy','sympy'])
	dErspacedalphaj = lambdify([r,qi,qj,alphai,alphaj,alphamax],dErspacedalphaj,modules=['numpy','sympy'])

        # numerical 

        def fun(r,qi,qj,alphai,alphaj,alphamax):
                return Erspace(r,qi,qj,alphai,alphaj,alphamax)

        def dfundr(r,qi,qj,alphai,alphaj,alphamax):
                # differential variable 
                x=r
                h = 1e-5
                return (fun(x+h,qi,qj,alphai,alphaj,alphamax)-fun(x-h,qi,qj,alphai,alphaj,alphamax))/(2*h)

        ana_Erspace = Erspace(value_r,value_qi,value_qj,value_alphai,value_alphaj,value_alphamax)
        ana_dErspacedr = dErspacedr(value_r,value_qi,value_qj,value_alphai,value_alphaj,value_alphamax)
        num_Erspace = fun(value_r,value_qi,value_qj,value_alphai,value_alphaj,value_alphamax)
        num_dErspacedr = -1*dfundr(value_r,value_qi,value_qj,value_alphai,value_alphaj,value_alphamax)
      
	if debug != None:
                print 'ana_Erspace = '+str(ana_Erspace)
                print 'ana_dErspacedr = '+str(ana_dErspacedr)
                print 'num_Erspace = '+str(num_Erspace)
                print 'num_dErspacedr = '+str(num_dErspacedr)
        return ana_Erspace,ana_dErspacedr,num_Erspace,num_dErspacedr

def sympy_E_self(value_qi,value_alphai,debug=None):
        # analytic 
        qi = symbols('qi')
        alphai = symbols('alphai')
        alphamax = symbols('alphamax')
        Eself = Function('Eself')(qi,alphai,alphamax)
        dEselfalphai = Function('dEselfdalphai')(qi,alphai,alphamax)
        #
        Eself1 = -1/2.*qi**2.*(2./sqrt(pi))*(sqrt((alphai*alphamax)/(alphai+alphamax)))
	Eself2 = -1/2.*qi**2.*(2./sqrt(pi))*(sqrt((alphai*alphai)/(alphai+alphai)))

        print 'Eself1 = '
        print Eself1
        dEself1dalphai = simplify(-1*diff(Eself1,alphai))
        print 'dEself1dalphai = '
        print dEself1dalphai
	
	print 'Eself2 = '
	print Eself2
        dEself2dalphai = simplify(-1*diff(Eself2,alphai))
        print 'dEself2dalphai = '
        print dEself2dalphai
	return Eself1,Eself2,dEself1dalphai,dEself2dalphai
if __name__ == "__main__":
        def main():
                # functionality test 
                print sympy_E_rspace(2.0,1,1,1,1,2,debug=True)
		print sympy_E_self(2.0,1,debug=True)

        main()

