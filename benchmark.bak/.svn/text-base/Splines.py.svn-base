#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
import numpy as np
import mpmath as mp
class H3SplineSuper(object): #(new style class necessary for super to work)
    # This is a Superclass that is not useful as stand alone class
    # since it does not contain boundary conditions which are specified 
    # in the subclasses
    def __init__(self,st,vals):
        self.st=st
        self.vals=vals

    def addInterior(self,m,rhs):
        st=self.st
        vals=self.vals
        imax=st.rows-1
        n=imax*4
        for i in range(1,imax):
            x=st[i]
            m[i*4:(i+1)*4,(i-1)*4:(i+1)*4]  =Matrix((
                [1, x, x**2,   x**3,0 ,0,    0,       0],
                [0, 0,    0,      0,1, x, x**2,    x**3],
                [0, 1,  2*x, 3*x**2,0,-1, -2*x, -3*x**2],
                [0, 0,    2,    6*x,0, 0,-  -2,    -6*x]
                ))
            rhs[i*4:i*4+2,0:1]=Matrix([vals[i],vals[i]])
        return(m,rhs)
    def r_poly(self,r):
        m,rhs=self.system()
        st=self.st
        imax=st.rows-1
        p=m.LUsolve(rhs)
        r=Symbol("r")
        sol=Piecewise((p[0]+r*p[1]+r**2*p[2]+r**3*p[3],r<st[1]))
        for i in range(1,imax):
            i4=i*4
            sol=Piecewise((sol,r<st[i]),(p[i4]+r*p[i4+1]+r**2*p[i4+2]+r**3*p[i4+3],r>=st[i]))
        return(sol)

class H3Spline01(H3SplineSuper):
    # this class implements a special form of hermite spline which satisfies the following conditions:
    # - At the points in st it takes the values given in vals
    # - At st[0] the first derivative will have the value dsl
    # - At st[-1] the first derivative will have the value dsr
    # - The derivatives up to order 2 are continous 

    def __init__(self,st,vals,dsl,dsr):
        pp=super(H3Spline01,self)
        pp.__init__(st,vals)
        self.dsl=dsl
        self.dsr=dsr

    def system(self):
        st  =self.st
        vals=self.vals
        dsl =self.dsl
        dsr =self.dsr

        imax=st.rows-1
        n=imax*4
        m=zeros(n)
        rhs=zeros((n,1))
        xl=st[0]
        m[0:2,0:4]  =Matrix(([1,xl,xl**2,  xl**3],[ 0, 1, 2*xl, 3*xl**2]))
        xr=st[-1]
        m[2:4,n-4:n]=Matrix(([1,xr,xr**2,  xr**3],[ 0, 1, 2*xr, 3*xr**2]))
        rhs[0]=vals[0]
        rhs[1]=dsl
        rhs[2]=vals[-1]
        rhs[3]=dsr
        m,rhs=self.addInterior(m,rhs)
        return(m,rhs)

class H3Spline02(H3SplineSuper):
    # this class implements a special form of hermite spline which satisfies the following conditions:
    # - At the points in st it takes the values given in vals
    # - At st[0] the >>second<< derivative will have the value ddsl
    # - At st[-1] the >>second<< derivative will have the value ddsr
    # - The derivatives up to order 2 are continous 

    def __init__(self,st,vals,ddsl,ddsr):
        pp=super(H3Spline02,self)
        pp.__init__(st,vals)
        self.ddsl=ddsl
        self.ddsr=ddsr

    def system(self):
        st  =self.st
        vals=self.vals
        ddsl =self.ddsl
        ddsr =self.ddsr

        imax=st.rows-1
        n=imax*4
        m=zeros(n)
        rhs=zeros((n,1))
        xl=st[0]
        m[0:2,0:4]  =Matrix(([1,xl,xl**2,  xl**3],[ 0, 0, 2, 6*xl]))
        xr=st[-1]
        m[2:4,n-4:n]=Matrix(([1,xr,xr**2,  xr**3],[ 0, 0, 2, 6*xr]))
        rhs[0]=vals[0]
        rhs[1]=ddsl
        rhs[2]=vals[-1]
        rhs[3]=ddsr
        m,rhs=self.addInterior(m,rhs)
        return(m,rhs)
