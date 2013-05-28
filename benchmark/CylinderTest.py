#!/usr/bin/python
# vim: set expandtab ts=4
from Cylinder import *
import unittest 
from sympy import bspline_basis
from Spherical import *
#from sympy import *
#from sympy import pi
import mpmath as mp
import numpy as np

class CylinderTest(unittest.TestCase):
    def test_Dl_eq104(self):
        startset={5}
        alpha=5.1
        l=1
        etha=0.2
	S=Spherical()
	r,phi,theta=S.U
        C=Cylinder(etha,l,startset,S)
        w=1/sqrt(r)*besselj(-(l+Rational(1,2)),alpha*r)
        expr=C.Dl(w)+alpha**2*w
        #that means that w is an eigenfuction of Dl with the eigenvalue -alpha^2
        S.testnumrad(expr,etha,1,0.01)

    def test_sphericalBesselOde(self):
        # this tests eq. 145
        startset={5}
        alpha=5.1
        for l in range(1):
            nu=l+Rational(1,2)
            etha=0.2
	    S=Spherical()
	    r,phi,theta=S.U
            C=Cylinder(etha,l,startset,S)
            alphas=C.alphafinder(l,etha,startset)
            for alpha in alphas:
                co=l+Rational(1,2)
                c1=besselj(-co,alpha*etha)
                c2=besselj(+co,alpha*etha)
                w=1/sqrt(r)*C.cyl_maker(alpha)(r)
                exp=C.Dl(w)+alpha**2*w
                #that means that w is an eigenfuction of Dl with the eigenvalue -alpha^2
                S.testnumrad(exp,etha,1,0.01)


    def test_cyl_roots(self):
        # We test that the base functions are zero at the boundaries
        startset={9,18,28,38}
        #startset={5}
        l=2
        etha=0.2
	S=Spherical()
	r,phi,theta=S.U
        C=Cylinder(etha,l,startset,S)
        alphas=C.alphafinder(l,etha,startset)
        for alpha in alphas:
            alpha_cyl=C.cyl_maker(alpha)
            self.assertTrue(abs(alpha_cyl(1.0).evalf())<10**(-9))
            self.assertTrue(abs(alpha_cyl(etha).evalf())<10**(-9))
            #def f(x):#
            #    ex=diff(alpha_cyl(r),r)
            #    return(ex.subs(r,x).evalf())
            #mp.plot(f,[etha,1.0])

    def test_cprime(self):
        #  we test both sides of eq. 139
        etha=0.2
        l=1
        startset={5}# up to one (n=1) root only 
	S=Spherical()
	#r,phi,theta=S.U
        C=Cylinder(etha,l,startset,S)
        alpha=C.rootlist[0]
        cp=C.cprime_maker(alpha)
        #print(cp(1))
        ref=(-1)**l*2/(pi*alpha)*besselj(l+1.0/2,alpha*etha)/besselj(l+1.0/2,alpha)
        ref2=C.C_prime_l_one_half(0)
        self.assertTrue(abs(cp(1)-ref)<10**(-9))
        self.assertTrue(abs(cp(1)-ref2)<10**(-9))
        
        ref=(-1)**l*2/(pi*alpha*etha)
        ref2=C.C_prime_l_one_half_etha(0)
        self.assertTrue(abs(cp(etha)-ref)<10**(-9))
        self.assertTrue(abs(cp(etha)-ref2)<10**(-9))

    
        
if  __name__ == '__main__':
     unittest.main()
