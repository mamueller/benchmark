#!/usr/bin/python
# vim: set expandtab ts=4
from Splines import *
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
import mpmath as mp

class SplineTest(unittest.TestCase):

    def test_H301system(self):
        #etha=0.4
        etha=0.0
        r=1
        h=(r-etha)/3
        print(h)
        x0=etha
        x1=etha+h
        x2=etha+2*h
        ## plot a function and its spline approximation
        v=Symbol("v")
        g=sin(1.0/(r-etha)*4*pi*(v-etha))
        gp=diff(g,v)
        def gf(val):
            return(g.subs(v,val).evalf())
        def gpf(val):
            return(gp.subs(v,val).evalf())
        mref=Matrix(12,12,
            [1,x0,x0**2,  x0**3,0, 0,    0,       0,0, 0,    0,       0,
             0, 1, 2*x0,3*x0**2,0, 0,    0,       0,0, 0,    0,       0,
             0, 0,    0,      0,0, 0,    0,       0,1, r, r**2,    r**3,
             0, 0,    0,      0,0, 0,    0,       0,0, 1,2* r , 3* r**2,
             1,x1,x1**2,  x1**3,0, 0,    0,       0,0, 0,    0,       0,
             0, 0,    0,      0,1,x1,x1**2,   x1**3,0, 0,    0,       0,
             0, 1, 2*x1,3*x1**2,0,-1,-2*x1,-3*x1**2,0, 0,    0,       0,
             0, 0,    2,   6*x1,0, 0,-  -2,   -6*x1,0, 0,    0,       0,
             0, 0,    0,      0,1,x2,x2**2,   x2**3,0, 0,    0,       0,
             0, 0,    0,      0,0, 0,    0,       0,1,x2,x2**2,   x2**3,
             0, 0,    0,      0,0, 1, 2*x2, 3*x2**2,0,-1,-2*x2,-3*x2**2,
             0, 0,    0,      0,0, 0,    2,    6*x2,0, 0,-  -2,   -6*x2]
        )
        rhs_ref=Matrix([gf(etha),gpf(etha),gf(r),gpf(r),gf(x1),gf(x1),0,0,gf(x2),gf(x2),0,0])
        st=Matrix([x0,x1,x2,r])
        vals=map(gf,st)
        dsl=gpf(st[0])
        dsr=gpf(st[-1])
        S=H3Spline01(st,vals,dsl,dsr)
        m,rhs=S.system()
        self.assertTrue(m==mref)
        self.assertTrue(rhs==rhs_ref)
    def test_H301rpoly(self):
        # this time we don't chose equidistant points
        etha=0.0
        h1=0.1
        h2=0.2
        h3=0.3
        h4=0.4
        rmax=etha+h1+h2+h3+h4
        st=Matrix([etha,etha+h1,etha+h1+h2,etha+h1+h2+h3,rmax])
        vals=Matrix([0,1,-1,2,0])
        dsl=0
        dsr=0
        S=H3Spline01(st,vals,dsl,dsr)
        r=Symbol("r")
        pr=S.r_poly(r)
        def prf(v):
            return(pr.subs(r,v).evalf())
        # we check if the values at the given points are actually correct
        pvals=Matrix(map(prf,st))
        dv=(pvals-vals).norm(2)
        self.assertTrue(dv<10**-10)
        print(dv)
        dpr=diff(pr,r)
        def dprf(v):
            return(dpr.subs(r,v).evalf())
        # we check the values of the first derivative the begin and end
        self.assertTrue(dprf(etha)==dsl)
        self.assertTrue(dprf(rmax)==dsr)
        #mp.plot(prf,[etha,rmax])


    def test_H302system(self):
        # same as in ChandraTest but with 2 
        #etha=0.4
        etha=0.0
        r=1
        h=(r-etha)/3
        print(h)
        x0=etha
        x1=etha+h
        x2=etha+2*h
        ## plot a function and its spline approximation
        v=Symbol("v")
        g=sin(1.0/(r-etha)*4*pi*(v-etha))
        gp=diff(g,v)
        gpp=diff(g,v,2)
        def gf(val):
            return(g.subs(v,val).evalf())
        def gpf(val):
            return(gp.subs(v,val).evalf())
        def gppf(val):
            return(gpp.subs(v,val).evalf())
        mref=Matrix(12,12,
            [1,x0,x0**2,  x0**3,0, 0,    0,       0,0, 0,    0,       0,
             0, 0,    2,   6*x0,0, 0,    0,       0,0, 0,    0,       0,
             0, 0,    0,      0,0, 0,    0,       0,1, r, r**2,    r**3,
             0, 0,    0,      0,0, 0,    0,       0,0, 0,    2,     6*r,
             1,x1,x1**2,  x1**3,0, 0,    0,       0,0, 0,    0,       0,
             0, 0,    0,      0,1,x1,x1**2,   x1**3,0, 0,    0,       0,
             0, 1, 2*x1,3*x1**2,0,-1,-2*x1,-3*x1**2,0, 0,    0,       0,
             0, 0,    2,   6*x1,0, 0,-  -2,   -6*x1,0, 0,    0,       0,
             0, 0,    0,      0,1,x2,x2**2,   x2**3,0, 0,    0,       0,
             0, 0,    0,      0,0, 0,    0,       0,1,x2,x2**2,   x2**3,
             0, 0,    0,      0,0, 1, 2*x2, 3*x2**2,0,-1,-2*x2,-3*x2**2,
             0, 0,    0,      0,0, 0,    2,    6*x2,0, 0,-  -2,   -6*x2]
        )
        rhs_ref=Matrix([gf(etha),gppf(etha),gf(r),gppf(r),gf(x1),gf(x1),0,0,gf(x2),gf(x2),0,0])
        st=Matrix([x0,x1,x2,r])
        vals=map(gf,st)
        ddsl=gppf(st[0])
        ddsr=gppf(st[-1])
        S=H3Spline02(st,vals,ddsl,ddsr)
        m,rhs=S.system()
        self.assertTrue(m==mref)
        self.assertTrue(rhs==rhs_ref)

    def test_H302rpoly(self):
        # this time we don't chose equidistant points
        etha=0.0
        h1=0.1
        h2=0.2
        h3=0.3
        h4=0.4
        rmax=etha+h1+h2+h3+h4
        st=Matrix([etha,etha+h1,etha+h1+h2,etha+h1+h2+h3,rmax])
        vals=Matrix([0,1,-1,2,0])
        ddsl=0
        ddsr=0
        S=H3Spline02(st,vals,ddsl,ddsr)
        r=Symbol("r")
        pr=S.r_poly(r)
        def prf(v):
            return(pr.subs(r,v).evalf())
        # we check if the values at the given points are actually correct
        pvals=Matrix(map(prf,st))
        dv=(pvals-vals).norm(2)
        self.assertTrue(dv<10**-10)
        print(dv)
        ddpr=diff(pr,r,2)
        def ddprf(v):
            return(ddpr.subs(r,v).evalf())
        # we check the values of the first derivative the begin and end
        self.assertTrue(ddprf(etha)==ddsl)
        self.assertTrue(ddprf(rmax)==ddsr)
        #mp.plot(prf,[etha,rmax])




if  __name__ == '__main__':
     unittest.main()
