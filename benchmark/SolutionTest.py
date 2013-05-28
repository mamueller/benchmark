#!/usr/bin/python
# vim: set expandtab ts=4
from ParametrizedTestCase import *
from Spherical import *
class SolutionTest(ParametrizedTestCase):
#    def test_solenoidal(self):
#        s=self.param
#        vb=s.expr
#        rmin=s.rmin
#        rmax=s.rmax
#        sp=s.coords
#        res=sp.phys_div(vb)
#        print(res)
#        sp.testnumshell(res,rmin,rmax,1e-7)
    
    def freebound(self,rval): 
        s=self.param
        vb=s.expr
	#print("vb="+str(vb))
        sp=s.coords
        r,phi,theta=sp.U

	# We compute the stresstensor:
        # vb is given by its physical components
        # so we first translate it to roof componets
        #vbr=sp.phys2roof(vb)
        #gt=sp.cellarComponentsOfNablaOnRoofComponents(vbr)#.subs(r,rval) #gt means gradient transposed
        gt=vb.nabla()

	#print("gt="+str(gt))
        g=gt.transpose()
	#print("g="+str(g))
        S_cr=(g+gt)/2
	#print("S_cr="+str(S_cr))
        # the resulting tensor is given by a mixed set of cellar_roof components
        # so it can be applied to a test vector test_c 
        # given by its cellar components by
        # matrix multiplication s_c= S_cr*test_c
        # The resulting tension s_c is given by its cellar components.
        
        # Our test vector is the direction normal to the surface which is given by the 
	# cellar base vector e_r.
	# To apply our stresstensor using the matrix product we express e_r by its cellar components
        # with respect to the roof basis. e^r, e^phi and e^theta.
        np=Matrix(3,1,[1,0,0])
        nc=sp.phys2cellar(np)
	#print("nc="+str(nc))
        snc=S_cr*nc
        #print("snc="+str(snc))
	snp=sp.cellar2phys(snc)
        #print("snp="+str(snp))
        # we now check that the resulting tension has no tangential
        # components.
	prec=10e-7
	s_phi=snc[1]
        print("s_phi="+str(s_phi))
	s_theta=snc[2]
	if s_phi!=0: 
	    sp.testnumsphere(s_phi,prec)
	if s_pheta!=0: 
	    sp.testnumsphere(s_pheta,prec)

