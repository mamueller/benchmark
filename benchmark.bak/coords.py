from sympy import *
import unittest 
import numpy as np

class Spherical:
   def laplaceBeltrami(self,fundef):
       r,theta, phi = symbols("r theta phi")
       f=fundef(theta,phi)
       theta_part=1/(sin(theta))*diff(sin(theta)*diff(f,theta),theta)
       phi_part=1/((sin(theta))**2)*diff(diff(f,phi),phi)
       lb=Lambda((theta,phi),theta_part+phi_part)
       return(lb)
   def div(self,A):
       r,theta, phi = symbols("r theta phi",real=True)
       Ar=A[0]
       Aphi=A[1]
       Atheta=A[2]
       div=1/r**2*diff(r**2*Ar,r) +1/(r*sin(theta))*diff(sin(theta)*Atheta,theta) +1/(r*sin(theta))*diff(Aphi,phi)
       return(div)
   def vectorLaplace(self,A):
       r,theta, phi = symbols("r theta phi",real=True)
       Ar=A[0]
       Aphi=A[1]
       Atheta=A[2]
       res=self.grad(self.div(A))-self.rot(self.rot(A))
       return(res)
   def rot(self,A):
       r,theta, phi = symbols("r theta phi",real=True)
       Ar=A[0]
       Aphi=A[1]
       Atheta=A[2]
       rr=1/(r*sin(theta))*(diff(Aphi*sin(theta),theta)-diff(Atheta,phi))
       rphi=1/r*(diff((r*Atheta),r)-diff(Ar,theta))
       rtheta=1/r*(1/sin(theta)*diff(Ar,phi)-diff(r*Aphi,r))
       def res(r_val,phi_val,theta_val):
           subsdict={r:r_val,phi:phi_val,theta:theta_val}
           mat=Matrix(3,1,[rr.subs(subsdict),rphi.subs(subsdict),rtheta.subs(subsdict)])
           return(mat)
       return(res)
class NumTestError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
        return("The difference was: "+str(self.d))


def tor_gen(fang,T):
       r,theta, phi = symbols("r theta phi",real=True)
       Tr=0		#r component
       Tphi=-1/r*T*diff(fang,theta)#phi 
       Ttheta=1/(r*sin(theta))*T*diff(fang,phi) #theta component
       A=Matrix(3,1,[Tr,Tphi,Ttheta])
       return(A)
def pol_gen(l,fang,S):    
       r,theta, phi = symbols("r theta phi",real=True)
       Sr=l*(l+1)/r**2*S*fang
       Sphi=1/(r*sin(theta))*diff(S,r)*diff(fang,phi)
       Stheta=1/r*diff(S,r)*diff(fang,theta)
       A=Matrix(3,1,[Sr,Sphi,Stheta])
       return(A)
def tor_real(l,m,T):
       r,theta, phi = symbols("r theta phi",real=True)
       fang=Zlm(l,m, theta, phi)
       A=tor_gen(fang,T)
       return(A)
def pol_real(l,m,S):    
       r,theta, phi = symbols("r theta phi",real=True)
       fang=Zlm(l,m, theta, phi)
       A=pol_gen(l,fang,S)
       return(A)
def tor_comp(l,m,T):
       r,theta, phi = symbols("r theta phi",real=True)
       fang=Ylm(l,m, theta, phi)
       A=tor_gen(fang,T)
       return(A)
def pol_comp(l,m,S):    
       r,theta, phi = symbols("r theta phi",real=True)
       fang=Ylm(l,m, theta, phi)
       A=pol_gen(l,fang,S)
       return(A)
def testnumrad_eq(expr1,expr2,r_min,r_max,prec):
    #print("expr1",expr1)
    #print("expr2",expr2)
    r=Symbol("r")
    def f1(x):
        return((expr1.subs(r,x)).evalf())
    def f2(x):
        return((expr2.subs(r,x)).evalf())
    args=np.linspace(r_min,r_max,50)
    vals1=list2numpy(map(f1,args))
    vals2=list2numpy(map(f2,args))
    test=(vals1-vals2)/(vals1+vals2)
    maxd=max(map(abs,test))
    if maxd >=prec:
        raise NumTestError(maxd)

    

def testnumrad(expr,r_min,r_max,prec):
    r=Symbol("r")
    def f(x):
        return((expr.subs(r,x)).evalf())
    args=np.linspace(r_min,r_max,50)
    vals=map(f,args)
    for val in vals:
        test=val.evalf()
        if not(Abs(test)<prec):
            print(val,test)
            print(args,vals)
            raise NumTestError(test)
def testnumsphere(expr,prec):
    n=16
    increment=pi/n
    for i in range(1,n):
        for j in range(1,n):
            theta=float(i*increment)
            phi=float(j*increment)
            test=simplify(exp)
            if not(Abs(test)<prec):
                raise NumTestError(test)
def testnumshell(expr,r_min,r_max,prec):
    r,theta, phi = symbols("r theta phi",real=True)
    def f(x):
        return(expr.subs(r,x))
    n=10
    expressions=map(f,np.linspace(r_min,r_max,n))
    precs=np.ones(n)*prec
    map(testnumsphere,expressions,precs)
    

