#!/usr/bin/python
# vim: set expandtab ts=4
from Cylinder import *
import unittest 
from sympy import bspline_basis
from Spherical import *
import mpmath as mp

class ChandraTest(unittest.TestCase):
    #def test_solenoidal(self):
    #    sp=Spherical()
    #    r,theta, phi = symbols("r theta phi",real=True)
    #    for l in range(1,2):
    #        for m in range(-l,l):
    #            T=r
    #            A=tor_comp(l,m,T)
    #            res=simplify(sp.div(A)) #test if the function is solenoidal
    def test_diff(self):
        #startset={9,18,28,38}
        startset={5}
        etha=0.2
        l=2
	S=Spherical()
	r,phi,theta=S.U
        C=Cylinder(etha,l,startset,S)
        F,W=C.radialpart()
    def test_cyl_roots(self):
        #startset={9,18,28,38}
        startset={5}
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
        #print(diff(ac,r))
    
    def test_chandraCl(self):
        # we reproduce the solutions mentioned in Chandrasekhar p.245 Table XXII
        # we start with the case of two free surfaces
        #startset={49,18,28,38}
        mp.dps =30 

        # set a startvalue for the numeric solver that computes
        # the roots the number of found roots also determines the number of
        # cylinderfunction used in the solution
	S=Spherical()
	r,phi,theta=S.U
        startset={5} 
        Clref={0.2:
               [5.211e3,
               5.708e3,
               8.882e3,
               1.400e4,
               2.121e4]
               ,
               0.3:
               [8.503e3,
                7.113e3,
                9.552e3,
                1.428e4,
                2.131e4,
                3.089e4]
               ,
               0.4:
               [1.682e4,
                1.091e4,
                1.196e4,
                1.585e4,
                2.227e4,
                3.143e4,
                4.365e4]
               ,
               0.5:
               [4.188e4,
                2.181e4,
                1.924e4,
                2.146e4,
                2.673e4,
                3.492e4,
                4.629e4,
                6.125e4,
                8.027e4,
                1.039e5,
                1.325e5,
                1.669e5,
                2.074e5,
                2.545e5,
                3.099e5]
               ,
               0.6:
               [1.403e5,
                6.133e4,
                4.424e4,
                4.076e4,
                4.313e4,
                4.945e4,
                5.933e4,
                7.292e4,
                9.057e4,
                1.128e5,
                1.401e5,
                1.732e5,
                2.126e5,
                2.590e5,
                3.131e5]
               ,
               0.8:
               [7.789e6,
                2.753e6,
                1.500e6,
                1.005e6,
                7.656e5,
                6.368e5,
                5.651e5,
                5.270e5,
                5.109e5,
                5.104e5,
                5.223e5,
                5.448e5,
                5.767e5,
                6.178e5,
                6.678e5]
              }
        for etha in Clref.viewkeys():
            #print(etha)
            for l in range(1,len(Clref[etha])):
                #print("l=",l)
                C=Cylinder(etha,l,startset,S)
                Cl=C.computeCl()
                self.assertTrue((abs(Cl-Clref[etha][l-1])/Clref[etha][l-1])<0.002)



    def test_B2B4(self):
        etha=0.2
        l=1
        #startset={5}# up to one (n=1) root only 
        startset={9,18,28,38}
	S=Spherical()
	r,phi,theta=S.U
        C=Cylinder(etha,l,startset,S)
        #C=Cylinder(etha,l,startset)
        for j in range(len(C.rootlist)):
            alpha=C.rootlist[j]
            B1,B2,B3,B4=C.bs(alpha)

            self.assertTrue(abs((B1+B2+B3+B4).evalf())<10**(-9)) 
            self.assertTrue(abs(B2.evalf()-C.B2(j))<10**(-9))
            self.assertTrue(abs(B4.evalf()-C.B4(j))<10**(-9))

if  __name__ == '__main__':
     unittest.main()


