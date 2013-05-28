#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
from sympy import bspline_basis
from Spherical import *
import mpmath as mp
#from ParametrizedTestCase import *
from RigidLidRigidCMBTest import *
from FreeSlipTest  import *
from Solution import *
from Splines import *
Sp=Spherical()
r,phi,theta=Sp.U
###########################################################################################
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
vb=Sp.pol_real(l,m,T)
BsplineSol=Solution(vb,r_min,r_max,Sp)
############################################################################
# this time we don't chose equidistant points
etha=0.2 # take care not to set this to zero because this will cause a NaN
assert(etha!=0)
h1=0.1
h2=0.2
h3=0.3
h4=0.4
rmax=etha+h1+h2+h3+h4
st=Matrix([etha,etha+h1,etha+h1+h2,etha+h1+h2+h3,rmax])
vals=Matrix([0,1,-1,2,0])
S=H3Spline01(st,vals,0,0)
Sr=S.r_poly(r)
l=1;m=1
vb=Sp.pol_real(l,m,Sr)
HsplineSol=Solution(vb,etha,rmax,Sp)
############################################################################
W=H3Spline02(st,vals,0,0)
Sr=r*W.r_poly(r)/(l*(l+1))
vb=Sp.pol_real(l,m,Sr)
HsplineSol2=Solution(vb,etha,rmax,Sp)
############################################################################
r_min=10
r_max=20
Sr=r-r_min
l=1;m=0
#vb=Sp.pol_real(l,m,Sr)
#vb=Matrix(3,1,[Sr,0,0])
v_r,v_phi,v_theta=symbols("u_r,u_phi,u_theta")
vb=Tensor(Sp,["phys"],{(1,):v_r(r,phi,theta),(2,):v_phi(r,phi,theta),(3,):v_theta(r,phi,theta)})
TestSol=Solution(vb,r_min,r_max,Sp)

############################################################################
def suite():
    s = unittest.TestSuite()
    #s.addTest(ParametrizedTestCase.parametrize(RigidLidRigidCMBTest, param=BsplineSol))
    #s.addTest(ParametrizedTestCase.parametrize(RigidLidRigidCMBTest, param=HsplineSol))
    s.addTest(ParametrizedTestCase.parametrize(FreeSlipTest, param=TestSol))
    #s.addTest(ParametrizedTestCase.parametrize(FreeSlipTest, param=HsplineSol2))
    return(s)

if  __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
