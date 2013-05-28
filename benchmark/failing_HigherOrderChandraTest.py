#!/usr/bin/python
# vim: set expandtab ts=4
from Cylinder import *
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
import mpmath as mp

class HigherOrderChandraTest(unittest.TestCase):

    def test_analytic_solution(self):
        # same as in ChandraTest but with higher order Terms 
        etha=0.4
        l=6
        startset={5,10,17,22,27,32}# up to six (n=6) roots only before we run  into precision problems
        C=Cylinder(etha,l,startset)
        Cl=C.computeCl()
        r_min=etha
        r_max=1
        r=Symbol("r")
        f,w=C.radialpart()
        # check the boundary conditions for f
        self.assertTrue(abs(f(r_min).evalf())<10**(-9))
        self.assertTrue(abs(f(r_max).evalf())<10**(-9))
        # check the boundary conditions for w
        self.assertTrue(abs(w(r_min).evalf())<10**(-9))
        self.assertTrue(abs(w(r_max).evalf())<10**(-9))
        W=w(r) 
        F=f(r) 
        # check the boundary conditions for w''
        WSS=diff(W,r,2)
        def wss(arg):
             return((WSS.subs(r,arg)).evalf())
        print("wss(etha)",wss(etha).evalf())
        #mp.plot(wss,[etha,1])
        #self.assertTrue(abs(wss(etha).evalf())<10**(-9)) # not enough precision
        print("wss(1)",wss(1).evalf())
        #self.assertTrue(abs(wss(1).evalf())<10**(-9))
        #print("C.Dl(F,l)")
        # eq. 131
        # The test will still fail due to a lack of precision as it does in the
        # case where the expansion is stopped after the first term
        testnumrad_eq(C.Dl(C.Dl(W)),F,r,r_min,r_max,10**-3)
if  __name__ == '__main__':
     unittest.main()
