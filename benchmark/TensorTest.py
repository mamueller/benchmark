#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from Spherical import *

class TensorTest(unittest.TestCase):
###########################################################
    def test_Tensor_initialization(self):
        sp=Spherical()
        # test some allowed use cases
        u=Tensor(sp,["cellar"],{})#empty tensor (or null)
        u=Tensor(sp,["cellar"],{0:1})
        u=Tensor(sp,["cart"],{0:1})
        u=Tensor(sp,["cellar","roof"],{(0,0):1,(1,1):1,(2,2):1}) #identitytensor
        # retrieve values from the tensor
        #ucr=u.components(["cellar","roof"])
        
        #test if the correct exception is raised for wrong input
        with self.assertRaises(UnknownComponentType):
        # typo
            Tensor(sp,["cellarr"],{(1,):1})
        with self.assertRaises(ComponentTypesTupelMismatch):
            # wrong number of component types
             Tensor(sp,["cellar"],{(1,1):1})

###########################################################
    def test_Vector_initialization(self):
        sp=Spherical()
        # test some allowed use cases
        u=Tensor(sp,["cellar"],{0:1})
        u=Tensor(sp,["cart"],{0:1})
        u=Tensor(sp,["cellar","roof"],{(0,0):1}) 
        

###########################################################
    def test_equal(self):
        # we test if a==b behaves like expected
        sp=Spherical()
        sp2=Spherical()
        u=Tensor(sp,["cellar"],{(1,):1})
        v=Tensor(sp2,["cellar"],{(1,):1})
        self.assertTrue(u==v)
        u=Tensor(sp,["cellar"],{(1,):1})
        v=Tensor(sp,["cellar"],{(1,):1})
        self.assertTrue(u==v)
        w=Tensor(sp,["cellar"],{(1,):2})
        self.assertTrue(not(u==w))
        w=Tensor(sp,["cellar"],{(2,):1})
        self.assertTrue(not(u==w))
        w=Tensor(sp,["roof"],{(1,):1})
        self.assertTrue(not(u==w))

        # test if a zero component is treated like a missing one
        u=Tensor(sp,["cellar"],{})
        v=Tensor(sp2,["cellar"],{(0,):0,(1,):0,(2,):0})
        self.assertTrue(u==v)
###########################################################

    def test_raise_and_lower(self):
        sp=Spherical()
        u=Tensor(sp,["cellar"],{(1,):1})
        v=u.raise_first_index().lower_first_index()
        self.assertTrue(u==v)
        v=u.raise_last_index().lower_first_index()
        self.assertTrue(u==v)
        
###########################################################
    def testScalarMul(self):
        sp=Spherical()
        u=Tensor(sp,["roof"],{(1,):1})
        v=u*5
        w=Tensor(sp,["roof"],{(1,):5})
        self.assertEqual(v,w)
        v=5*u
        self.assertEqual(v,w)
    
###########################################################
    def testAdd(self):
        sp=Spherical()
        u=Tensor(sp,["roof"],{(1,):1})
        v=Tensor(sp,["roof"],{(1,):2})
        ref=Tensor(sp,["roof"],{(1,):3})
        res=u+v
        self.assertEqual(res,ref)
        res=v+u
        self.assertEqual(res,ref)
        u=Tensor(sp,["cart"],{(1,):1})
        v=Tensor(sp,["cart"],{(1,):2})
        ref=Tensor(sp,["cart"],{(1,):3})
        res=u+v
        self.assertEqual(res,ref)
        res=v+u
        u=Tensor(sp,["roof"],{(1,):1})
        v=Tensor(sp,["cart"],{(1,):2})
        with self.assertRaises(NotImplementedError):
            u+v
        u=Tensor(sp,["roof","roof"],{(1,1):1})
        v=Tensor(sp,["cart"],{(1,):2})
        with self.assertRaises(ArgumentSizeError):
           u+v
        
        with self.assertRaises(ArgumentTypeError):
           u+3 
        
###########################################################
    def testSub(self):
        sp=Spherical()
        u=Tensor(sp,["roof"],{(1,):1})
        v=Tensor(sp,["roof"],{(1,):2})
        ref=Tensor(sp,["roof"],{(1,):-1})
        res=u-v
        self.assertEqual(res,ref)
    
###########################################################
    def test_extractVector(self):
        sp=Spherical()
        u  =Tensor(sp,["cart","cart"],{(0,0):1,(0,1):2})
        self.assertEqual(u.extractVector((0,"*")),Tensor(sp,["cart"],{(0,):1,(1,):2}))
        self.assertEqual(u.extractVector(("*",0)),Tensor(sp,["cart"],{(0,):1}))
        self.assertEqual(u.extractVector(("*",1)),Tensor(sp,["cart"],{(0,):2}))
        self.assertEqual(u.extractVector(("*",2)),Tensor(sp,["cart"],{}))
        
        u  =Tensor(sp,["roof","cart","cellar"],{(0,0,0):1,(0,1,0):2})
        self.assertEqual(u.extractVector((0,"*",0)),Tensor(sp,["cart"],{(0,):1,(1,):2}))
        self.assertEqual(u.extractVector(("*",0,0)),Tensor(sp,["roof"],{(0,):1}))
        
        with self.assertRaises(ArgumentSizeError):
            u  =Tensor(sp,["cart"],{(0,):1,(1,):2})
            u.extractVector(("*",1))
        #self.assertEqual(u.extractVector(("*",2)),Tensor(sp,["cart"],{}))

###########################################################

    def test_scalarProduct(self):
        sp=Spherical()

        # second order Tensor times vector
        u  =Tensor(sp,["cart","cart"],{(0,0):1})
        w  =Tensor(sp,["cart","cart"],{(1,1):1})
        v  =Tensor(sp,["cart"],{(1,):1})
        
        ref=Tensor(sp,["cart"],{})
        self.assertEqual(u|v,ref)
        self.assertEqual(v|u,ref)
        
        self.assertEqual(u|u,u)
        self.assertEqual(u|w,Tensor(sp,["cart","cart"],{}))
        # roof and cellar components 
        u  =Tensor(sp,["roof","cellar"],{(0,0):1})
        self.assertEqual(u|u,u)
        
        v  =Tensor(sp,["roof"],{(1,):1})
        self.assertEqual(v|u,Tensor(sp,["roof"],{}))
        # mixed with cartesian components 
        # (We allow this although it is very unusual 
        # e.g. e_x e_phi is a valid dyadic product 
        # the components in this form thus refer to roof and cartesian base
        # vectors. Usually the whole component set would be either cartesian
        # or some combination of roof and cellar
        w  =Tensor(sp,["roof","cellar"],{(1,1):1})
        v  =Tensor(sp,["roof"],{(1,):1})



        
###########################################################
    def test_Vector_scalarProduct(self):
        sp=Spherical()
        
        u=Tensor(sp,["cellar"],{(1):1})
        v=Tensor(sp,["roof"],{(1):1})
        self.assertEqual(u|v,1)
        self.assertEqual(v|u,1)
        
        u=Tensor(sp,["cellar"],{(1,):1})
        v=Tensor(sp,["cellar"],{(1,):1})
        self.assertEqual(u|v,sp.g_rr()[1,1])
        
        u=Tensor(sp,["roof"],{(1,):1})
        v=Tensor(sp,["roof"],{(1,):1})
        self.assertEqual(u|v,sp.g_cc()[1,1])
       
        #orthogonal
        u=Tensor(sp,["roof"],{(1,):1})
        v=Tensor(sp,["roof"],{(2,):1})
        self.assertEqual(u|v,0)
        
        ###test transformation invariance ####
        eur=Tensor(sp,["cellar"],{(0,):1}) 
        eor=Tensor(sp,["roof"],{(0,):1})
        self.assertEqual(eur|eor,1)
        eur_cart=eur.transform2(["cart"])
        eor_cart=eor.transform2(["cart"])
        print(eur_cart)
        print(eor_cart)
        sp.testnumshell((eur_cart|eor_cart)-1,0.1,2,1e-7)
        self.assertEqual(eur_cart|eor_cart,1)
###########################################################
    def test_indextupel(self):
        t=indextupel(1)
        self.assertEqual(t,{(0,),(1,),(2,)})
        t=indextupel(2)
        self.assertEqual(t,{(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)})
        t=indextupel(3)
        self.assertEqual(t,{(0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1),(0,1,2),(0,2,0),(0,2,1),(0,2,2),\
	                  (1,0,0),(1,0,1),(1,0,2),(1,1,0),(1,1,1),(1,1,2),(1,2,0),(1,2,1),(1,2,2),\
	                  (2,0,0),(2,0,1),(2,0,2),(2,1,0),(2,1,1),(2,1,2),(2,2,0),(2,2,1),(2,2,2),\
			  })
              
###########################################################
    def test_add_index(self):
        t=indextupel(1)
        res=add_index(t,0)
        self.assertEqual(res,t)
        
        res=add_index(t,1)
        self.assertEqual(res,indextupel(2))
        
        res=add_index(t)
        self.assertEqual(res,indextupel(2))


###########################################################
    def test_del_index(self):
        t=indextupel(1)
        with self.assertRaises(ToShortTupelError):
            res=del_index(t,0)
        t={}
        with self.assertRaises(EmptyIndexSetError):
            res=del_index(t,0)
        
        t=indextupel(2)
        res=del_index(t,0)
        self.assertEqual(res,indextupel(1))
        
        t=indextupel(2)
        res=del_index(t,1)
        self.assertEqual(res,indextupel(1))
        
        t=indextupel(2)
        res=del_index(t,-1)
        self.assertEqual(res,indextupel(1))
        

        t=indextupel(3)
        res=del_index(t,0)
        self.assertEqual(res,indextupel(2))

        t=indextupel(3)
        res=del_index(t,1)
        self.assertEqual(res,indextupel(2))

        t=indextupel(3)
        res=del_index(t,2)
        self.assertEqual(res,indextupel(2))

        t=indextupel(3)
        res=del_index(t,-1)
        self.assertEqual(res,indextupel(2))

        t={(1,2,3)}
        res=del_index(t,1)
        self.assertEqual(res,{(1,3)})


###########################################################
    def test_transform2_vector(self):
        sp=Spherical()
        
        X=Tensor(sp,["roof"],{(0,):1})
        Y=X.transform2(["cart"])
        self.assertEqual(Y,sp.t_gc[0])
        
        X=sp.t_gc[0] #the first cellar base vector in cartesian components
        Y=X.transform2(["roof"])
        zero=Y-Tensor(sp,["roof"],{(0,):1})
        zc=zero.components	
        sp.testnumshell(zc[(0,)],0.1,2,1e-7)
        
            
        X=Tensor(sp,["cellar"],{(0,):1})
        Y=X.transform2(["cart"])
        self.assertEqual(Y,sp.t_gr[0])
        
        X=sp.t_gr[0] #the first roof base vector in cartesian components
        Y=X.transform2(["cellar"]) # the cellar components must be 1,0,0
        zero=Y-Tensor(sp,["cellar"],{(0,):1})
        zc=zero.components	
        sp.testnumshell(zc[(0,)],0.1,2,1e-7)
    
###########################################################
    def test_transform2_secondOrderTensors(self):
        
        ##  construct a two dimensional testcase
        #X=Tensor(sp,["cart","cart"],{(0,0):1})
        #Y=X.transform2(["roof","roof"])
        
        ## we first define a linear mapping P_x that acts as a projector of any vector to the 
        # x axis then we can express the application of P_x to v by the dot
        # product *
        #      P_x =  e  e * v
        #              x  x
        #
        # assume that we have a vector v given by its cartesian components
        #    x    y     z
        # v=v e +v e + v e
        #      x    y     z
           ##                
        #                                        x
        # When we apply P_x to v we get P_x(v) =v e
        #                                          x
        # The cartesian components of P_x = 1*ex ex + 0*ex ey +0 ...
        #
        # Thus if v is expressed in cartesian components 
        # the result P_x(v) expressed in cartesian coordinates is given by the
        # matrix multiplication res=P_x v
        sp=Spherical()
        Pcartcart=Tensor(sp,["cart","cart"],{(0,0):1})
        # we chose the cartesian components of the first roof base vector
        v_cart=sp.t_gr[0]
        res_cart=Pcartcart|v_cart
        
#        # now we want to express the same Projection in terms of 
#        # spherical roof components 
#        Prr=sp.cart_cart2rr(Pcartcart)
#        # If v is expressed in cellar components 
#        v_c=sp.cart2cellar(v_cart)
#        # the application of 
#        # P can again be expressed by a matrix multiplication 
#        #                      i           i
#        # using the fact that g * g  =delta
#        #                          j       j
#        res_r=Prr*v_c
#        # however the result is given in terms of the cellar base vectors 
#        # (in it its roof components) 
#        # We have to transform the result back to cartesian coordinates
#        res_cart_new=sp.roof2cart(res_r)
#        #print("new="+str(res_cart_new))
#        self.assertTrue(res_cart==res_cart_new)
#        
#        
#        # we also test the backtransformation of the componentset of the whole
#        # tensor
#        Pcartcart_new=sp.rr2CartCart(Prr)
#        #print("2Pcartcart_new="+str(Pcartcart_new))
#        self.assertTrue(Pcartcart==Pcartcart_new)
#
#        #print("new="+str(res_cart_new.subs({phi:pi/4,theta:pi/4,r:3})))


        
        
###########################################################
    def test_vec_grad(self):  
        sp=Spherical()
        x,y,z=sp.X
        xu,yu,zu=sp.XofU
        r,phi,theta=sp.U
        # we start with the following vector valued function f given in 
        # cartesian coordinates
        fX=Tensor(sp,["cart"],{(2,):x}) #=x*e
        #                                    z
        # if a vector v is given in cartesian coordinates
        # the vector of change of fX deltafX according to a change of position
        # deltaX can be expressend by the Tensor equation 
        # delta_fX=A|delta_X+O(|delta_X|**2) where A
        # is the second order tensor to represent the 
        # gradient of v and | is the scalar product 
        
        # Althouhg it is simple to derive A in this case we 
        #                                              i
        # compute it with Simmonds eq(4.10): nabla v= g v,i
        # where nabla v is the transposed gradient
        #            x       x            x     x          
        # nabla v= =e  v,x =e  (x*e ),x =e e  =e e  =Tensor(sp,["cart","cart"],{(0,2):1})
        #                          z        z     z
        # since 
        Aref=Tensor(sp,["cart","cart"],{(2,0):1})
        # first we express the (cartesian) components of f as functions of r,phi and theta
        cx=fX.components    
        cu={}
        for k in cx.keys():
            cu[k]=simplify(cx[k].subs({x:xu,y:yu,z:zu}).subs(cos(phi)**2,(1-sin(phi)**2)))
        fUc=Tensor(sp,["cart"],cu)
        # the vector f is still given with respect to the cartesian basis
        # since the vectorgradient is very easy to define if the vector is given
        # in roof components we compute these first
        fUr=fUc.transform2(["roof"])
	cr=fUr.nabla()
        
        #cr=sp.cellarComponentsOfNablaOnRoofComponents(fUr)
        ## transform into roofroof
        #rr=sp.cr2rr(cr)
        ## transform back to cartesian basis
        #cartcart=sp.rr2CartCart(rr)
        #print(cartcart)
###########################################################
    def test_partder(self):
        sp=Spherical()
        r,phi,theta=sp.U
        
        v=Tensor(sp,["roof"],{(0,):1})
        drv=v.partder(0)#with respect to r
        ref=Tensor(sp,["roof"],{(0,):0,(1,):0,(2,):0})
        self.assertEqual(drv,ref)
        
        v=Tensor(sp,["cellar"],{(0,):1})
        # this is the roof base vector e^r
        drv=v.partder(0)#with respect to r
        ref=Tensor(sp,["cellar"],{})
        self.assertEqual(drv,ref)
       
        v=Tensor(sp,["roof"],{(1,):1}) #e_phi
        drv=v.partder(0)#with respect to r
        ref=Tensor(sp,["roof"],{(0,):0,(1,):1/r,(2,):0})
        self.assertEqual(drv,ref)
        
        v=Tensor(sp,["cellar"],{(1,):1}) #e_phi
        drv=v.partder(0)#with respect to r
        ref=Tensor(sp,["cellar"],{(0,):0,(1,):-1/r,(2,):0})
        self.assertEqual(drv,ref)
        print(drv)
        
        
        


if  __name__ == '__main__':
    unittest.main()
