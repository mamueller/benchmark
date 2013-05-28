#!/usr/bin/python
import unittest 
from Spherical import *
from Cartesian import *

class OperatorTest(unittest.TestCase):
    #def setUp(self):
    # for all tests we create 2 instances of a general coordinate system
    # a cartesian and a spherical one
    def test_tensordiv(self):  
        '''This function tests some of the basic properties of the divergence of a second order tensor in general coordinates. '''
        C=Cartesian()
        #########x,y,z=C.X
        ux,uy,uz=C.U
        # we start with a tensor valued function T given by its cartesian
        # components  
        T=Tensor(C,["roof","roof"],{\
        (0,0):ux,(0,1):uy,(0,2):uz,\
        (1,0):ux,(1,1):uy,(1,2):uz,\
        (2,0):ux,(2,1):uy,(2,2):uz\
        })
        refres=Tensor(C,["roof"],{(0,):1,(1,):1,(2,):1})
        res=T.div()
        self.assertTrue(res==refres)
        
        T=Tensor(C,["roof","roof"],{\
        (0,0):ux,(0,1):ux,(0,2):ux,\
        (1,0):uy,(1,1):uy,(1,2):uy,\
        (2,0):uz,(2,1):uz,(2,2):uz\
        })
        refres=Tensor(C,["roof"],{(0,):3,(1,):3,(2,):3})
        res=T.div()
        print(res)
        self.assertTrue(res==refres)
        
        # we now test that the divergence of a tensor is an invariant
        S=Spherical()
        r,phi,theta=S.U
        xu,yu,zu=S.XofU
        x,y,z=S.X
        TX_cartcart=Tensor(S,["cart","cart"],{\
        (0,0):x,(0,1):x,(0,2):x,\
        (1,0):y,(1,1):y,(1,2):y,\
        (2,0):z,(2,1):z,(2,2):z\
        })
        #{x,x,x],[y,y,y],[z,z,z]])
        # we first express the componenst in terms of the new variables
        # Note that they will still refer to the cartesian base vectors.
        TU_cartcart=TX_cartcart.subs({x:xu,y:yu,z:zu})

        print("TU_cartcart="+str(TU_cartcart))
        ## now we transform the componets to roofroof components
        TU_rr=TU_cartcart.transform2(["roof","roof"])    
        print("TU_rr="+str(TU_rr))
        #TU_roofroof=S.cart_cart2rr(TU_cartcart)
        ## print("TU_roofroof=\n"+str(TU_roofroof))
        ## now we can apply the divergence operator to the components
            #res_roof=S.roofroof_div_T(TU_roofroof)
        ## print("roof_div=\n"+str(res_roof))
        ## The result is expressed in roof components
        ## to be able to compare it to the original result we have to transform it to the
        ## original cartesian basis.
        #res_cart=S.roof2cart(res_roof)
        ## print("cart_div=\n"+str(res_cart))
        #self.assertEqual(res_cart,refres)
if  __name__ == '__main__':
     unittest.main()
    # vim:set ff=unix expandtab ts=4 sw=4:
