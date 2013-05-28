#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
import mpmath as mp

class OperatorTest(unittest.TestCase):
    def test_LaplaceBeltrami(self):
        sp=Spherical()
        theta=Symbol("theta",real=True)
        phi=Symbol("phi",real=True)
        #l=0 m=0 
        f=Lambda((theta,phi),Ylm(0,0, theta, phi))
        res=sp.laplaceBeltrami(f)
        self.assertEqual(res,Lambda((theta,phi),0))
        # we now test if the spherical harmonics fulfill the eigenvalue problem of 
        # the laplaceBeltrami operator lb on the sphere
        # lb(ylm(l,m))=-(l*(l+1))*ylm(l,m) 
        l,m=1,1
        for l in range(1,2):
            for m in range(-l,l):
                f=Lambda((theta,phi),Ylm(l,m, theta, phi))
                g=Lambda((theta,phi),Zlm(l,m, theta, phi))
                res_f=sp.laplaceBeltrami(f)
                res_g=sp.laplaceBeltrami(g)
                ev=-l*(l+1)
                testval_f=simplify(res_f(phi,theta)-ev*f(phi,theta))
                testval_g=simplify(res_g(phi,theta)-ev*g(phi,theta))
                #print testval_f
                #print testval_g
                # sometimes the ensuing expressions are to complicated for simplify to be able              
                #to detect equality. In that case  we evaluate the functions on certain points
                if testval_f!=Lambda((theta,phi),0):
                    expr_f=expand(res_f(theta,phi)-ev*f(theta,phi))
                    expr_g=expand(res_g(theta,phi)-ev*g(theta,phi))
                    testnumsphere(expr_f,0.001)
                    testnumsphere(expr_g,0.001)
               
        
                 #  print("f=",f)
                 #  print("res=",res)
                 #  print("testval=",testval)
    def test_Divergence(self):
        sp=Spherical()
        r,theta, phi = symbols("r theta phi",real=True)
        A=Matrix(3,1,[r,0,0])
        res=sp.div(A)
        self.assertTrue(res==3)
        A=Matrix(3,1,[1,phi,0])
        res=sp.div(A)
        self.assertTrue(res==2/r + 1/(r*sin(theta)))
        A=Matrix(3,1,[1,0,theta])
        res=sp.div(A)
        self.assertTrue(res==(theta*cos(theta) + sin(theta))/(r*sin(theta)) + 2/r)
        #now we test  solenoidal functions 
        #1 torroidal
        # we form a basis with the help of spherical harmonics T(r) is arbitrary function of radius
        l,m=1,1
        for l in range(1,2):
            for m in range(-l,l):
                T=r
                A=tor_comp(l,m,T)
                res=simplify(sp.div(A)) #test if the function is solenoidal
                #print("A=",A)
                #print("res=",res)
                self.assertTrue(res==0)
                #solenoidal part
                S=r#arbitrary function of radius
                A=pol_comp(l,m,S)
                res=simplify(sp.div(A)) #test if the function is solenoidal
                #print("A=",A)
                #print("res=",res)
                self.assertTrue(res==0)
                T=r
                A=tor_real(l,m,T)
                res=simplify(sp.div(A)) #test if the function is solenoidal
                #print("A=",A)
                #print("res=",res)
                self.assertTrue(res==0)
                #poloidal part
                S=r#arbitrary function of radius
                A=pol_real(l,m,S)
                res=simplify(sp.div(A)) #test if the function is solenoidal
                #print("A=",A)
                #print("res=",res)
	    self.assertTrue(res==0)
    def test_solenoidal(self):
        sp=Spherical()
        r,theta, phi = symbols("r theta phi",real=True)
        for l in range(1,2):
            for m in range(-l,l):
                T=r
                A=tor_comp(l,m,T)
                res=simplify(sp.div(A)) #test if the function is solenoidal

    def t_v(rmin,rmax,jmax):
        etha=rmin/rmax

    def test_rot(self):
        sp=Spherical()
        r,theta, phi = symbols("r theta phi",real=True)
        A=Matrix(3,1,[0,sin(theta)*r,0]) #a vector field symmetric to the z axis
        res=sp.rot(A)
        #print("A=",A)
        #print("res=",res(r,phi,theta))
        self.assertTrue(res(r,phi,theta)==Matrix(3,1,[2*cos(theta),0,-2*sin(theta)]))
        
        A=Matrix(3,1,[0,0,cos(theta)*cos(phi)*r]) #a vector field symmetric to the y axis
        res=sp.rot(A)
        #print("A=",A)
        #print(res(1,0,0))
        self.assertTrue(res(1,0,0)==Matrix(3,1,[0,2,0]))
    def test_bspline_solution3(self):
        r,theta, phi = symbols("r theta phi",real=True)
        sp=Spherical()
        ## We try to build an analytic solution for the stokes system with solenoidal basisfunctions that
        ## furthermore fulfill the boundary conditions for a rigid lid
        ## 
        ## let the basisfunction be vb
        ## first we create the radial part
        ## To this end we use splines that we combine from the b-spline basis
        ## One of the advantages of b-spline basis functions is that they and their derivatives 
        ## vanish at the ends of their support interval
        ## Therefore it is very easy to fulfill the boundary conditions 
        ## 
        r_min=10
        r_max=20
        # cubic splines should be smooth enough for the stokes problem since only 2 derivatives are needed
        # on the other hand the radial part of the soluton composed of the bsplines are differentiated 
        # once before to  build up the 3d Field. So
        degree=3  
        knotnumber=degree+2
        bs=[]
        for i in range(r_min,r_max+2-knotnumber):
            knots=range(i,i+knotnumber)
            b=bspline_basis(degree,knots,0,r)
            bs.append(b)
        #print(bs)
        # we now combine the basisfunctions to a singe radial function
        T=piecewise_fold(bs[0]+bs[-1])
        #print(T)
        
        ### we have to show that the functions are solenoidal
        l=1;m=1
        vb=pol_real(l,m,T)
        res=sp.div(vb)
        testnumshell(res,r_min,r_max,10**-7)
        #print("vb=",vb)
        #print("res=",res)
        
        ### now we test for the boundary conditions
        ### vanishing of T at the inner and outer surface. 
        self.assertTrue(T.subs(r,r_min)==0)
        self.assertTrue(T.subs(r,r_max)==0)
        print(T.subs(r,r_max))
        ### vanishing of T' at the inner and outer surface. 
        dT=diff(T,r)
        self.assertTrue(dT.subs(r,r_min)==0)
        self.assertTrue(dT.subs(r,r_max)==0)
    
if  __name__ == '__main__':
     unittest.main()
