#!/usr/bin/python
# vim: set expandtab ts=4
# from sympy import *
#import unittest 
from Tensor import *
import numpy as np
class Spherical(Coords):
    def __init__(self):
        x,y,z=symbols("x,y,z")
        X=x,y,z
        # spherical part
        # Unfortunately there is a lot of confusion about the symbols used, most notably by
        # the opposite role which phi and theta play in mathematics and physics.
        # It is therefore necessary to derive all expressions from a single definition.
        # Here we stick to the convention in Chandrasekhars book, which uses the system most common
        # in physics. 
        
        # We define: 
        # phi         to be the azimuthal angle in the xy-plane from the x-axis with  0<=phi<2pi ,
        # theta       to be the polar (zenithal) angle  from the positive z-axis with 0<=theta<=pi, and 
        # r           to be distance (radius) from a point to the origin. 
        #
        # This is >> not << the convention commonly used in mathematics,where theta and phi are reversed.

        r=symbols("r",positive=True,real=True)
        phi,theta=symbols("phi,theta",real=True)
        # Order:
        # We always use the conventions 
        # A_x=A[0] A_y=A[1] A_z=A[2]
        # A_r=A[0] A_phi=A[1] A_theta=A[2]
        # regardless if A_r is a cellar, roof or physical component.
        # Transformation:
        # Let X and U be  ordered tupels
        #   |x|    |  r|
        # X=|y|  U=|phi|
        #   |z|    |theta|
        # then X(U) is given by:
        x  =  r*cos(phi)*sin(theta)
        y  =  r*sin(phi)*sin(theta)
        z  =           r*cos(theta)
        XofU=[x,y,z]
        U=[r,phi,theta]
        XofU=[x,y,z]
        sc=super(Spherical,self)
        inst=sc.__init__(X,U,XofU)
    
    def scalarSimp(self,exp):	
        r,phi,theta=self.U
        res=trigsimp(exp)
        res=res.subs(sin(phi)**2*sin(theta)**2-cos(phi)**2*cos(theta)**2+cos(phi)**2+cos(theta)**2,1)
	res=simplify(res.subs((sin(phi))**2,1-(cos(phi))**2))
        return(res)

    def __repr__(self):
        return("Spherical()")
   
    def tor_gen(self,fang,T):
           r,phi,theta=self.U
           Tr=0		#r component
           Tphi=-1/r*T*diff(fang,theta)#phi 
           Ttheta=1/(r*sin(theta))*T*diff(fang,phi) #theta component
           A=Matrix(3,1,[Tr,Tphi,Ttheta])
           return(A)
    def pol_gen(self,l,fang,S):    
           r,phi,theta=self.U
           Sr=l*(l+1)/r**2*S*fang
           Sphi=1/(r*sin(theta))*diff(S,r)*diff(fang,phi)
           Stheta=1/r*diff(S,r)*diff(fang,theta)
           A=Matrix(3,1,[Sr,Sphi,Stheta])
           return(A)
    def tor_real(self,l,m,T):
           r,phi,theta=self.U
           fang=Zlm(l,m, theta, phi)
           A=self.tor_gen(fang,T)
           return(A)
    def pol_real(self,l,m,S):    
           r,phi,theta=self.U
           fang=Zlm(l,m, theta, phi)
           A=self.pol_gen(l,fang,S)
           return(A)
    def tor_comp(self,l,m,T):
           r,phi,theta=self.U
           fang=Ylm(l,m, theta, phi)
           A=self.tor_gen(fang,T)
           return(A)
    def pol_comp(self,l,m,S):    
           r,phi,theta=self.U
           fang=Ylm(l,m, theta, phi)
           A=self.pol_gen(l,fang,S)
           return(A)

    def testnumsphere(self,expr,prec):
        r,phi,theta=self.U
        n=8
        increment=pi/n
        for i in range(1,n):
            for j in range(1,n):
                t=float(i*increment)
                p=float(j*increment)
                test=(expr.subs(theta,t)).subs(phi,p)
                if not(Abs(test)<prec):
                    raise NumTestError(test)
    
    def testnumsphere_r(self,expr,r_vals,prec):
        r,phi,theta=self.U
        def f(x):
            return(expr.subs(r,x))
        expressions= map(f,r_vals)
        n=len(r_vals)
        precs=np.ones(n)*prec
        map(self.testnumsphere,expressions,precs)


    def testnumshell(self,expr,r_min,r_max,prec):
        r,phi,theta = self.U
        def f(x):
            return(expr.subs(r,x))
        n=10
        expressions=map(f,np.linspace(r_min,r_max,n))
        precs=np.ones(n)*prec
        map(self.testnumsphere,expressions,precs)
     
    def testnumrad(self,expr1,r_min,r_max,prec):
        r,phi,theta = self.U
        def g(x):
            return(expr1.subs(r,x))
        args=np.linspace(r_min,r_max,10)
        vals=map(g,args)
        maxd=max(map(abs,vals))
        if maxd >=prec:
            raise NumTestError(maxd)
#    def laplaceBeltrami(self,fundef):
#        r,theta, phi = symbols("r theta phi")
#        f=fundef(theta,phi)
#        theta_part=1/(sin(theta))*diff(sin(theta)*diff(f,theta),theta)
#        phi_part=1/((sin(theta))**2)*diff(diff(f,phi),phi)
#        lb=Lambda((theta,phi),theta_part+phi_part)
#        return(lb)
#    def div(self,A):
#        # this function computes the divergence of a vector given in physical
#        # components.
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        div=1/r**2*diff(r**2*Ar,r) +1/(r*sin(theta))*diff(sin(theta)*Atheta,theta) +1/(r*sin(theta))*diff(Aphi,phi)
#        return(div)
#    def vectorLaplace(self,A):
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        res=self.grad(self.div(A))-self.rot(self.rot(A))
#        return(res)
#    def rot(self,A):
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        rr=1/(r*sin(theta))*(diff(Aphi*sin(theta),theta)-diff(Atheta,phi))
#        rphi=1/r*(diff((r*Atheta),r)-diff(Ar,theta))
#        rtheta=1/r*(1/sin(theta)*diff(Ar,phi)-diff(r*Aphi,r))
#        def res(r_val,phi_val,theta_val):
#            subsdict={r:r_val,phi:phi_val,theta:theta_val}
#            mat=Matrix(3,1,[rr.subs(subsdict),rphi.subs(subsdict),rtheta.subs(subsdict)])
#            return(mat)
#        return(res)
class NumTestError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
        return("The difference was: "+str(self.d))


def testnumrad_eq(expr1,expr2,sym,r_min,r_max,prec):
    #print("expr1",expr1)
    #print("expr2",expr2)
    def f1(x):
        return((expr1.subs(sym,x)).evalf())
    def f2(x):
        return((expr2.subs(sym,x)).evalf())
    args=np.linspace(r_min,r_max,50)
    vals1=list2numpy(map(f1,args))
    vals2=list2numpy(map(f2,args))
    test=(vals1-vals2)/(vals1+vals2)
    maxd=max(map(abs,test))
    if maxd >=prec:
        raise NumTestError(maxd)

    


