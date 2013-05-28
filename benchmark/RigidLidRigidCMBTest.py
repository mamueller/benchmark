#!/usr/bin/python
# vim: set expandtab ts=4
from ParametrizedTestCase import *
from Spherical import *
from RigidLidFreeCMBTest import *
class RigidLidRigidCMBTest(RigidLidFreeCMBTest):
    def test_inner_boundary(self):
        # We have to ensure that the velocityfield has zero velocity in both angular directions
        # at the boundary
        s=self.param
        vb=s.expr
        rmin=s.rmin
        rmax=s.rmax
        sp=s.coords
        # print(vb)
        # we first test that all components of the velocity are zero at the boundary
        sp.testnumsphere_r(vb[0],[rmin],1e-7)
        sp.testnumsphere_r(vb[1],[rmin],1e-7)
        sp.testnumsphere_r(vb[2],[rmin],1e-7)
        r,phi,theta=sp.U
        dr=diff(vb[0],r)
        # now we test that the derivative 
        sp.testnumsphere_r(dr,[rmin],1e-7)

	




