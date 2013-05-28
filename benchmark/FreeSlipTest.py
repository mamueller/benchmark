#!/usr/bin/python
# vim: set expandtab ts=4
from ParametrizedTestCase import *
from Spherical import *
from SolutionTest import *
from Exceptions import *

class FreeSlipTest(SolutionTest):
	
#    def test_inner_boundary(self):
#        s=self.param
#        vb=s.expr
#        rmin=s.rmin
#        sp=s.coords
#        # we first test that the radial component of the velocity is zero at the boundary
#        sp.testnumsphere_r(vb[0],[rmin],1e-7)
#        # now we look at the tangential stresses at the boundary
#        self.freebound(rmin)
#




        #raise(NotImplementedError)
    def test_outer_boundary(self):
        s=self.param
        vb=s.expr
        rmax=s.rmax
        sp=s.coords
        # we first test that the radial component of the velocity is zero at the boundary
#        sp.testnumsphere_r(vb[0],[rmax],1e-7)
        # now we look at the tangential stresses at the boundary
        self.freebound(rmax)
        #raise(NotImplementedError)



