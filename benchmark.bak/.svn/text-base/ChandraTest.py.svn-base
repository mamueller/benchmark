#!/usr/bin/python
# vim: set expandtab ts=4
from Cylinder import *
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
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
        C=Cylinder(etha,l,startset)
        F,W=C.radialpart()
    def test_cyl_roots(self):
        #startset={9,18,28,38}
        startset={5}
        r=Symbol("r")
        l=2
        etha=0.2
        C=Cylinder(etha,l,startset)
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
                C=Cylinder(etha,l,startset)
                Cl=C.computeCl()
                self.assertTrue((abs(Cl-Clref[etha][l-1])/Clref[etha][l-1])<0.002)



    def test_B2B4(self):
        etha=0.2
        l=1
        #startset={5}# up to one (n=1) root only 
        startset={9,18,28,38}
        C=Cylinder(etha,l,startset)
        for j in range(len(C.rootlist)):
            alpha=C.rootlist[j]
            B1,B2,B3,B4=C.bs(alpha)

            self.assertTrue(abs((B1+B2+B3+B4).evalf())<10**(-9)) 
            self.assertTrue(abs(B2.evalf()-C.B2(j))<10**(-9))
            self.assertTrue(abs(B4.evalf()-C.B4(j))<10**(-9))

class AnalyticalSolutionTest(unittest.TestCase):  
    def setUp(self):
        etha=0.2
        l=1
        startset={5}# up to one (n=1) root only 
        C=Cylinder(etha,l,startset)
        self.C=C
        f,w=C.radialpart()
        self.sol=(f,w)

    def test_analytic_solution(self):
        # - compute as many roots of the bessel functions as needed
        # - construct the cylinderfunctions up to order n
        # - compute the secular matrix (of order n) and find the root(s) Cl of its
        #   determinant. Since the matrix is symmetric the roots of the
        #   characteristic polynom must be real. Note however that the Cl are
        #   not identical with them. The Cl are (up to a constant)            
        #   the Rayleigh numbers for which F and W can  actually be constructed from the cylinder
        #   functions.
        # - compute the eigenvectors (the Ajs) by substituting the Cl back 
        #   in the secular matrix and computing its kernel.
        # - Use the Ajs to build an analytical solution for F
        # - determine the constants B1...B4 for all j <n and construct the
        #   basefunctions for W
        # - combine the W-base funtions with the Ajs to an analytical solution
        #   for W
        #   Test! that this is indeed an analytical solution by inserting it in ( eq. 132 Chandra)
        #   
        #   The test fails because according to the above procedure (and so in
        #   Chandrasekar book the solution is only sought in a week form that
        #   is only for a number of testfunctions the integral over the product
        #   vanishes. In Chandrasekhars book the number of testfunctions is
        #   actually 1!!!! Thus the ODE 132 is not fulfilled very exactly and
        #   this "analytical" solution should not be used as a benchmark to detect errors in a
        #   numerically obtained solution.
        
        #   Note that jmax is arbitrary and can in principle be increased ad
        #   libitum as long as we are able to compute the eigenvectors of the
        #   secular matrix. In chandrasekhar this has been done for order 1 and
        #   2 another possibility to vary the solutions is not to increase the
        #   matrix size but to chose different roots.
        #   
        # - After the solution for F and W are found we can extend them with the
        #   help of spherical harmonics to solution for the whole field
        #   This step is ambigous however since the parameter m of Yml is not
        #   determined by the procedure thus far since the linear problem does
        #   not distinguish these 2l+1 solutions. Every one of them represents a
        #   steady state solution.
        #   Only the parameter l has a physical relevance
        #   since the physically realized solutions will be those with the value
        #   l that leads to the minimal value of Cl (which l will minimize Cl
        #   also depends on the geometry of the sphere (rmin/rmax)

        #   We start with a the case n=1, where the matrix is one dimensional
        #   the eigenvector arbitrary and Cl is known from Chandrasekhars book 
        C=self.C
        f,w=self.sol
        r=Symbol("r")
        # check the boundary conditions for f
        self.assertTrue(abs(f(C.etha).evalf())<10**(-9))
        self.assertTrue(abs(f(1).evalf())<10**(-9))
        # check the boundary conditions for w
        self.assertTrue(abs(w(C.etha).evalf())<10**(-9))
        self.assertTrue(abs(w(1).evalf())<10**(-9))
        W=w(r) 
        F=f(r) 
        # check the boundary conditions for w''
        WSS=diff(W,r,2)
        def wss(arg):
             return(WSS.subs(r,arg))
        self.assertTrue(abs(wss(C.etha).evalf())<10**(-9))
        self.assertTrue(abs(wss(1).evalf())<10**(-9))

    def test_eq_131(self):
        C=self.C
        etha=C.etha
        f,w=self.sol
        r=Symbol("r")
        W=w(r) 
        F=f(r) 
        # This test will fail showing that the accuracy of the expansion is to
        # small if one takes into account only the first part 
        testnumrad((C.Dl(C.Dl(W))-F),etha,1,0.1)


    def test_eq_132(self):
        C=self.C
        l=C.l
        alpha=C.rootlist[0]
        tCl=alpha**6/l*(l+1)
        Cl=C.computeCl()
        print("testCl",tCl,Cl)
        etha=C.etha
        f,w=self.sol
        r=Symbol("r")
        W=w(r) 
        F=f(r) 
        # This test will fail showing that the accuracy of the expansion is to
        # small if one takes into account only the first part 
        testnumrad(1/r*C.Dl(F/r)+l*(l+1)*Cl*W,etha,1,0.01)


        # - Application as a benchmark for a convection code: 
        #   - We first observe that the convection problem in the Stokes
        #     approximation is linear. It follows therefore that the solutions
        #     obtained for the linearized Navier Stokes system
        #     in the stability analysis of ChandraSekhar are actually
        #     solutions not only for infinitesimal amplitudes but for every
        #     amplitude.
        #   - We secondly observe that since there are (for constant viscosity)
        #     no nonlinear terms in the equations that it is not possible to
        #     distinguish the many solutions of the linear problem by means of a
        #     solvebility condition derived from those nonlinear terms as in
        #     Busses fundamental 1975 paper. Consequently we are left with the
        #     following situation: 
        #     for every combination of  ratio etha=rmin/rmax the number and kind
        #     of roots, used for the construction of the radial basis functions, 
        #     we have one value l which leads to the minimal Rayleigh number and
        #     thus 2l+1 solution for steady state convection. That is 2l+1
        #     pairs for temperature and velocity, the pressure field should be
        #     obtainable however by integration of the momentum equation.
        #   - To test a solver one could use the analytic temperature solution
        #     to compute the righthandside of the momentum equation and start
        #     the Stokes solver. The ensuing numeric solutions for v and p can
        #     be tested against the analytic ones.
        #     Together with the given temperature field these v and p should not
        #     change the temperature in the convective solver (that includes T
        #     also) This is however actually only a test for the Stokes solver
        #     with the additional benefit of using a force field stemming from a
        #     temperature distribution.
        #   - Instead it  would be nice to actually have a time dependent solution where
        #     the initial state would produce a predictable time series for the
        #     fields considered. This would be an optimal test for a convective
        #     solver. 
        #     If we only focus on beginning and end a candidate for such a case
        #     would be "disturbed stable solution " where the system moves from
        #     one steady state to another triggered by a disturbance great
        #     enought. This could be achieved for instance by
        #     means of an increased negative temperature gradient starting
        #     convection (wiht the l corresponding to the lowest Raykeigh
        #     number. The convective solver should be able to reproduce the
        #     correct l. The numerical results can then be expanded into
        #     spherical harmonics and the correct order l ensured. Another
        #     option would be to trigger a special pattern Ylm by the form of the
        #     disturbance which could be choosen accordingly.
        #     The task of the solver would then to amplify the correct solution
        #     until the correct amplitude is reached. The correct 
        #     solution must preserve energy while the disturbance must not (to
        #     test anything at all). 
        
if  __name__ == '__main__':
     unittest.main()


