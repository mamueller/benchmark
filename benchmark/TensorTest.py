#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from Spherical import *
from Cartesian import *

class TensorTest(unittest.TestCase):
###########################################################
    def test_Tensor_initialization(self):
        sp=Spherical()
        # test some allowed use cases
        u=Tensor.Tensor(sp,["cellar"],{})#empty tensor (or null)
        u=Tensor.Tensor(sp,["cellar"],{0:1})
        u=Tensor.Tensor(sp,["cart"],{0:1})
        u=Tensor.Tensor(sp,["cellar","roof"],{(0,0):1,(1,1):1,(2,2):1}) #identitytensor
        # retrieve values from the tensor
        #ucr=u.components(["cellar","roof"])
        
        #test if the correct exception is raised for wrong input
        with self.assertRaises(UnknownComponentType):
        # typo
            Tensor.Tensor(sp,["cellarr"],{(1,):1})
        with self.assertRaises(ComponentTypesTupelMismatch):
            # wrong number of component types
             Tensor.Tensor(sp,["cellar"],{(1,1):1})

###########################################################
    def test_Vector_initialization(self):
        sp=Spherical()
        # test some allowed use cases
        u=Tensor.Tensor(sp,["cellar"],{0:1})
        u=Tensor.Tensor(sp,["cart"],{0:1})
        u=Tensor.Tensor(sp,["cellar","roof"],{(0,0):1}) 
        

###########################################################
    def test_equal(self):
        # we test if a==b behaves like expected
        sp=Spherical()
        sp2=Spherical()
        u=Tensor.Tensor(sp,["cellar"],{(1,):1})
        v=Tensor.Tensor(sp2,["cellar"],{(1,):1})
        self.assertTrue(u==v)
        u=Tensor.Tensor(sp,["cellar"],{(1,):1})
        v=Tensor.Tensor(sp,["cellar"],{(1,):1})
        self.assertTrue(u==v)
        w=Tensor.Tensor(sp,["cellar"],{(1,):2})
        self.assertTrue(not(u==w))
        w=Tensor.Tensor(sp,["cellar"],{(2,):1})
        self.assertTrue(not(u==w))
        w=Tensor.Tensor(sp,["roof"],{(1,):1})
        self.assertTrue(not(u==w))

        # test if a zero component is treated like a missing one
        u=Tensor.Tensor(sp,["cellar"],{})
        v=Tensor.Tensor(sp2,["cellar"],{(0,):0,(1,):0,(2,):0})
        self.assertTrue(u==v)
###########################################################

    def test_raise_and_lower(self):
        sp=Spherical()
        u=Tensor.Tensor(sp,["cellar"],{(1,):1})
        v=u.raise_first_index().lower_first_index()
        self.assertTrue(u==v)
        v=u.raise_last_index().lower_first_index()
        self.assertTrue(u==v)
        
###########################################################
    def test_ScalarMul(self):
        sp=Spherical()
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=u*5
        w=Tensor.Tensor(sp,["roof"],{(1,):5})
        self.assertEqual(v,w)
        v=5*u
        self.assertEqual(v,w)
    
###########################################################
    def testAdd(self):
        sp=Spherical()
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=Tensor.Tensor(sp,["roof"],{(1,):2})
        # make a copy for later comparisons
        cu=copy.deepcopy(u)
        cv=copy.deepcopy(v)
        ref=Tensor.Tensor(sp,["roof"],{(1,):3})
        res=u+v
        self.assertEqual(res,ref)
        # test that we do not modyfy the ingredients inadvertatly
        # so u=cu and v=cv should still hold
        self.assertEqual(u,cu)
        self.assertEqual(v,cv)
        # now use __radd__ by changing the order of the parts
        res_reverse=v+u
        self.assertEqual(res_reverse,ref)
        # test that we do not modyfy the ingredients inadvertatly
        # so u=cu and v=cv should still hold
        self.assertEqual(u,cu)
        self.assertEqual(v,cv)
        
        u=Tensor.Tensor(Cartesian(),['cellar', 'roof'],{(0, 2): 1})
        v=Tensor.Tensor(Cartesian(),['cellar', 'roof'],{(1, 2): 1})

        cu=copy.deepcopy(u)
        cv=copy.deepcopy(v)
        ref=Tensor.Tensor(Cartesian(),['cellar', 'roof'],{(1, 2): 1, (0, 2): 1})
        res=u+v
        self.assertEqual(res,ref)
        # test that we do not modyfy the ingredients inadvertatly
        # so u=cu and v=cv should still hold
        self.assertEqual(u,cu)
        self.assertEqual(v,cv)
        # now use __radd__ by changing the order of the parts
        res_reverse=v+u
        self.assertEqual(res_reverse,ref)
        # test that we do not modyfy the ingredients inadvertatly
        # so u=cu and v=cv should still hold
        self.assertEqual(u,cu)
        self.assertEqual(v,cv)


        self.assertEqual(res,ref)
        u=Tensor.Tensor(sp,["cart"],{(1,):1})
        v=Tensor.Tensor(sp,["cart"],{(1,):2})
        ref=Tensor.Tensor(sp,["cart"],{(1,):3})
        res=u+v
        self.assertEqual(res,ref)
        res=v+u
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=Tensor.Tensor(sp,["cart"],{(1,):2})
        with self.assertRaises(NotImplementedError):
            u+v
        u=Tensor.Tensor(sp,["roof","roof"],{(1,1):1})
        v=Tensor.Tensor(sp,["cart"],{(1,):2})
        with self.assertRaises(ArgumentSizeError):
           u+v
        
        with self.assertRaises(ArgumentTypeError):
           u+3 
        
###########################################################
    def testSub(self):
        sp=Spherical()
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=Tensor.Tensor(sp,["roof"],{(1,):2})
        ref=Tensor.Tensor(sp,["roof"],{(1,):-1})
        res=u-v
        self.assertEqual(res,ref)
    
###########################################################
    def test_extractVector(self):
        sp=Spherical()
        u  =Tensor.Tensor(sp,["cart","cart"],{(0,0):1,(0,1):2})
        self.assertEqual(u.extractVector((0,"*")),Tensor.Tensor(sp,["cart"],{(0,):1,(1,):2}))
        self.assertEqual(u.extractVector(("*",0)),Tensor.Tensor(sp,["cart"],{(0,):1}))
        self.assertEqual(u.extractVector(("*",1)),Tensor.Tensor(sp,["cart"],{(0,):2}))
        self.assertEqual(u.extractVector(("*",2)),Tensor.Tensor(sp,["cart"],{}))
        
        u  =Tensor.Tensor(sp,["roof","cart","cellar"],{(0,0,0):1,(0,1,0):2})
        self.assertEqual(u.extractVector((0,"*",0)),Tensor.Tensor(sp,["cart"],{(0,):1,(1,):2}))
        self.assertEqual(u.extractVector(("*",0,0)),Tensor.Tensor(sp,["roof"],{(0,):1}))
        
        with self.assertRaises(ArgumentSizeError):
            u  =Tensor.Tensor(sp,["cart"],{(0,):1,(1,):2})
            u.extractVector(("*",1))
        #self.assertEqual(u.extractVector(("*",2)),Tensor.Tensor(sp,["cart"],{}))

###########################################################
    def test_outerProduct(self):
        sp=Spherical()
        u=Tensor.Tensor(sp,["cellar"],{(1,):1})
        v=Tensor.Tensor(sp,["roof"],{(1,):1})
        res=u*v
        self.assertEqual(res,Tensor.Tensor(sp,["cellar","roof"],{(1,1):1}))

        u=Tensor.Tensor(sp,["cellar","cellar"],{(0,1):3,(1,1):4})
        v=Tensor.Tensor(sp,["roof","cellar"],{(1,0):2})
        self.assertEqual(u*v,Tensor.Tensor(sp,["cellar","cellar","roof","cellar"],{(0,1,1,0):6,(1,1,1,0):8}))
        
        # test that res=a*b does not change when we modify a or b later on 
        # this ensures that we did not forget to copy the igredients instead of referencing them 
        a=Tensor.Tensor(sp,["cellar"],{(1,):1})
        b=Tensor.Tensor(sp,["roof"],{(1,):1})
        res=a*b
        self.assertEqual(res,Tensor.Tensor(sp,["cellar","roof"],{(1,1):1}))

        sp2=Spherical()
        a=Tensor.Tensor(sp2,["roof"],{(2,):2})
        b=Tensor.Tensor(sp2,["cellar"],{(0,):2})
        self.assertEqual(res,Tensor.Tensor(sp,["cellar","roof"],{(1,1):1}))


###########################################################
    def test_innerProduct(self):
        sp=Spherical()

        # second order Tensor times vector
        u  =Tensor.Tensor(sp,["cart","cart"],{(0,0):1})
        w  =Tensor.Tensor(sp,["cart","cart"],{(1,1):1})
        v  =Tensor.Tensor(sp,["cart"],{(1,):1})
        
        ref=Tensor.Tensor(sp,["cart"],{})
        self.assertEqual(u|v,ref)
        self.assertEqual(v|u,ref)
        
        self.assertEqual(u|u,u)
        self.assertEqual(u|w,Tensor.Tensor(sp,["cart","cart"],{}))
        # roof and cellar components 
        u  =Tensor.Tensor(sp,["roof","cellar"],{(0,0):1})
        self.assertEqual(u|u,u)
        
        v  =Tensor.Tensor(sp,["roof"],{(1,):1})
        self.assertEqual(v|u,Tensor.Tensor(sp,["roof"],{}))
        # mixed with cartesian components 
        # (We allow this although it is very unusual 
        # e.g. e_x e_phi is a valid dyadic product 
        # the components in this form thus refer to roof and cartesian base
        # vectors. Usually the whole component set would be either cartesian
        # or some combination of roof and cellar
        w  =Tensor.Tensor(sp,["roof","cellar"],{(1,1):1})
        v  =Tensor.Tensor(sp,["roof"],{(1,):1})
        
        # test that res=a!b does not change when we modify a or b later on 
        # this ensures that we did not forget to copy the igredients instead of referencing them 
        a=Tensor.Tensor(sp,["cellar"],{(1,):2})
        b=Tensor.Tensor(sp,["roof"],{(1,):3})
        res=a|b
        self.assertEqual(res,6)

        sp2=Spherical()
        a=Tensor.Tensor(sp2,["roof"],{(2,):2})
        b=Tensor.Tensor(sp2,["cellar"],{(0,):2})
        self.assertEqual(res,6)



        
###########################################################
    def test_Vector_innerProduct(self):
        sp=Spherical()
        
        u=Tensor.Tensor(sp,["cellar"],{(1):1})
        v=Tensor.Tensor(sp,["roof"],{(1):1})
        self.assertEqual(u|v,1)
        self.assertEqual(v|u,1)
        
        u=Tensor.Tensor(sp,["cellar"],{(1,):1})
        v=Tensor.Tensor(sp,["cellar"],{(1,):1})
        self.assertEqual(u|v,sp.g_rr()[1,1])
        
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=Tensor.Tensor(sp,["roof"],{(1,):1})
        self.assertEqual(u|v,sp.g_cc()[1,1])
       
        #orthogonal
        u=Tensor.Tensor(sp,["roof"],{(1,):1})
        v=Tensor.Tensor(sp,["roof"],{(2,):1})
        self.assertEqual(u|v,0)
        
        ###test transformation invariance ####
        eur=Tensor.Tensor(sp,["cellar"],{(0,):1}) 
        eor=Tensor.Tensor(sp,["roof"],{(0,):1})
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
        
        X=Tensor.Tensor(sp,["roof"],{(0,):1})
        Y=X.transform2(["cart"])
        self.assertEqual(Y,sp.t_gc[0])
        
        X=sp.t_gc[0] #the first cellar base vector in cartesian components
        Y=X.transform2(["roof"])
        zero=Y-Tensor.Tensor(sp,["roof"],{(0,):1})
        zc=zero.components	
        sp.testnumshell(zc[(0,)],0.1,2,1e-7)
        
            
        X=Tensor.Tensor(sp,["cellar"],{(0,):1})
        Y=X.transform2(["cart"])
        self.assertEqual(Y,sp.t_gr[0])
        
        X=sp.t_gr[0] #the first roof base vector in cartesian components
        Y=X.transform2(["cellar"]) # the cellar components must be 1,0,0
        zero=Y-Tensor.Tensor(sp,["cellar"],{(0,):1})
        zc=zero.components	
        sp.testnumshell(zc[(0,)],0.1,2,1e-7)
    
###########################################################
    def test_transform2_secondOrderTensors(self):
        
        ##  construct a two dimensional testcase
        #X=Tensor.Tensor(sp,["cart","cart"],{(0,0):1})
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
        Pcartcart=Tensor.Tensor(sp,["cart","cart"],{(0,0):1})
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
        
    def test_nabla(self):
        # we start with the following vector valued function f given in 
        # cartesian coordinates
        c=Cartesian()
        x,y,z=c.U
        fX=Tensor.Tensor(c,["roof"],{(2,):x}) #=x*e_z
        pp("fX",locals())
        fXn=fX.nabla()


###########################################################
    def test_vec_grad_cart(self):  
        # we first verify that our implementation 
        # works for the comparativly simple case
        # of cartisian coordinates. Here the base vectors are constant
        # so the derivative of a vector has only to take into account
        # the components. Also roof and cellar base vectors are 
        # indistinguishable in this case
        
        # Consider a vector valued function F at positon X.
        # If we change X by delta_X also F will change by delta_F
        # The vector of change delta_F according to a change of position
        # delta_X can be expressend by a linear function A acting on 
        # delta_X and a remainder. We call the linear part dF
        #
        # delta_F = A(delta_X) + O(|delta_X|**2) 
        #
        #         = dF         + O(|delta_X|**2) 
        # 
        # we can express the application of A on delta_X by
        # the scalar product. So we get 
        #
        # d_F=A|delta_X
        #
        # where A is the second order tensor to represent the 
        # gradient of v and | is the scalar product 
        # To test our implementation we choose a linear function
        # F which then coincides with its own derivative.
        # we then have 
        # 
        # delta_F=d_F
        #
        c=Cartesian()
        x,y,z=c.U
        # for F we chose
        #                                         z       
        F=Tensor.Tensor(c,["roof"],{(2,):x}) #=x*e   =x*e   (roof and cellar base vectors are equal in this case) 
        #                                  z
        # for X0 we chose the origin
        x0=0;y0=0;z0=0
        # for X1 we chose 
        x1=1;y1=0;z1=0
        # which means that  delta_X =X1-X1
        # which in cartesian coordinates can be expressed very easily
        delta_X=Tensor.Tensor(c,["roof"],{(0,):x1-x0,(1,):y1-y0,(2,):z1-z0})
        #
        # now we first compute F0=F(X0)=F(x0,y0,z0)
        F0=F.subs({x:x0,y:y0,z:z0})
        F1=F.subs({x:x1,y:y1,z:z1})
        delta_F=F1-F0
        pp("delta_F",locals())

        # now we compute the same result using the grad operator 
        A=F.grad()
        d_F=A|delta_X
        
        self.assertEqual(d_F,delta_F)
        
        
        # we now do this for several examles for F and delta_X
        # where the common property is that F is linear
        # for F we chose
        F=Tensor.Tensor(c,["roof"],{(2,):x+y}) 
        # for X0 we chose the origin
        x0=1;y0=0;z0=0
        # for X1 we chose 
        x1=2;y1=2;z1=1
        # which means that  delta_X =X1-X1
        # which in cartesian coordinates can be expressed very easily
        delta_X=Tensor.Tensor(c,["roof"],{(0,):x1-x0,(1,):y1-y0,(2,):z1-z0})
        #
        # now we first compute F0=F(X0)=F(x0,y0,z0)
        F0=F.subs({x:x0,y:y0,z:z0})
        F1=F.subs({x:x1,y:y1,z:z1})
        delta_F=F1-F0

        # now we compute the same result using the grad operator 
        A=F.grad()
        d_F=A|delta_X
        
        self.assertEqual(d_F,delta_F)

###########################################################
    def test_vec_grad_lin(self):  
        # we now adapt the result of test_vec_grad_cart to the 
        # case of spherical coordinates
        sp=Spherical()
        x,y,z=sp.X
        # for F we chose
        #                                         z       
        F=Tensor.Tensor(sp,["cart"],{(2,):x}) #=x*e   =x*e   
        # we first express the vector function by its cartesian coordinates as 
        # above but this time using the fact that every instance of Tensor can 
        # convert its components to the cartesian coordinate system regardless 
        # of its own coordinate frame. This is posible since all coordinate frame
        # objects  contain information about the cartesian components of 
        # their bases.
        #                                  z
        # for X0 we chose the origin
        x0=0;y0=0;z0=0
        # for X1 we chose 
        x1=1;y1=0;z1=0
        # which means that  delta_X =X1-X1
        # which in cartesian coordinates can be expressed very easily
        delta_X=Tensor.Tensor(sp,["cart"],{(0,):x1-x0,(1,):y1-y0,(2,):z1-z0})
        #
        # now we first compute F0=F(X0)=F(x0,y0,z0)
        F0=F.subs({x:x0,y:y0,z:z0})
        F1=F.subs({x:x1,y:y1,z:z1})
        delta_F=F1-F0
        pp("delta_F",locals())

        # now we compute the same result using the grad operator 
        # but in order to do so we have to convert it to eihter 
        # roof or cellar components. (we will chose roof here) 
        # Note that the conversions will consist of  2 steps
        # 1) express x,y,z by the new variables
        # 2) express ex, ey, ez by the cellar base vectors
        XofU=sp.XofU
        # 1):
        F=F.subs({x:XofU[0],y:XofU[1],z:XofU[2]}) 
        # 2):
        F=F.transform2(["roof"])                   
        pp("F",locals())
        A=F.grad()
        d_F=A|delta_X
        pp("d_F",locals())
        # to be able to compare this result to delta_F
        # we have to express delta_F in the new coordinates too
        delta_F=delta_F.subs({x:XofU[0],y:XofU[1],z:XofU[2]}) 
        delta_F=delta_F.transform2(["roof"])                   
        pp("delta_F",locals())
        self.assertEqual(d_F,delta_F)
        
        
        # we now do this for other examles for F and delta_X
        # where the common property is that F is linear
        # for F we chose
        F=Tensor.Tensor(sp,["cart"],{(2,):x+y}) 
        # for X0 we chose the origin
        x0=1;y0=0;z0=0
        # for X1 we chose 
        x1=2;y1=2;z1=1
        # which means that  delta_X =X1-X0
        # which in cartesian coordinates can be expressed very easily
        delta_X=Tensor.Tensor(sp,["cart"],{(0,):x1-x0,(1,):y1-y0,(2,):z1-z0})
        #
        # now we first compute F0=F(X0)=F(x0,y0,z0)
        F0=F.subs({x:x0,y:y0,z:z0})
        F1=F.subs({x:x1,y:y1,z:z1})
        delta_F=F1-F0
        pp("delta_F",locals())

        pp("delta_X",locals())
        #
        XofU=sp.XofU
        # 1):
        F=F.subs({x:XofU[0],y:XofU[1],z:XofU[2]}) 
        # 2):
        F=F.transform2(["roof"])                   
        pp("F",locals())
        A=F.grad()
        d_F=A|delta_X
        pp("d_F",locals())
        # to be able to compare this result to delta_F
        # we have to express delta_F in the new coordinates too
        delta_F=delta_F.subs({x:XofU[0],y:XofU[1],z:XofU[2]}) 
        delta_F=delta_F.transform2(["roof"])                   
        pp("delta_F",locals())
        self.assertEqual(d_F,delta_F)



###########################################################
    def test_transpose(self):
        c=Cartesian()
        # note that roof and cellar
        # components are equal in this case
        x,y,z=c.U
        A=Tensor.Tensor(c,["roof","cellar"],{(2,0):1})
        a=Tensor.Tensor(c,["cellar"],{(2,):x})
        b=Tensor.Tensor(c,["roof"],{(2,):x})
        res=a|A|b        
        res_trans=b|A.transpose()|a  #definition of the transposed Tensor       
        self.assertEqual(res,res_trans)
        
        sp=Spherical()
        r,phi,theta=sp.U
        
        A=Tensor.Tensor(sp,["roof","cellar"],{(2,0):1})
        a=Tensor.Tensor(sp,["cellar"],{(0,):r,(1,):phi,(2,):r*phi})
        b=Tensor.Tensor(sp,["roof"],{(2,):phi})
        res=a|A|b        
        res_trans=b|A.transpose()|a  #definition of the transposed Tensor       
        self.assertEqual(res,res_trans)

        A=Tensor.Tensor(sp,["cart","cart"],{(2,0):1})
        a=Tensor.Tensor(sp,["cellar"],{(1,):r})
        b=Tensor.Tensor(sp,["roof"],{(2,):phi})
        res=a|A|b        
        res_trans=b|A.transpose()|a  #definition of the transposed Tensor       
        self.assertEqual(res,res_trans)

        A=Tensor.Tensor(sp,["roof","roof"],{(2,0):1})
        a=Tensor.Tensor(sp,["cellar"],{(0,):r})
        b=Tensor.Tensor(sp,["cellar"],{(2,):phi})
        res=a|A|b        
        res_trans=b|A.transpose()|a  #definition of the transposed Tensor       
        self.assertEqual(res,res_trans)
###########################################################
    def test_partder(self):
        sp=Spherical()
        x,y,z=sp.X
        r,phi,theta=sp.U
        
        v=Tensor.Tensor(sp,["roof"],{(0,):1})
        drv=v.partder(0)#with respect to r
        ref=Tensor.Tensor(sp,["roof"],{(0,):0,(1,):0,(2,):0})
        self.assertEqual(drv,ref)
        
        v=Tensor.Tensor(sp,["cellar"],{(0,):1})
        # this is the roof base vector e^r
        drv=v.partder(0)#with respect to r
        ref=Tensor.Tensor(sp,["cellar"],{})
        self.assertEqual(drv,ref)
       
        v=Tensor.Tensor(sp,["roof"],{(1,):1}) #e_phi
        drv=v.partder(0)#with respect to r
        ref=Tensor.Tensor(sp,["roof"],{(0,):0,(1,):1/r,(2,):0})
        self.assertEqual(drv,ref)
        
        v=Tensor.Tensor(sp,["cellar"],{(1,):1}) #e_phi
        drv=v.partder(0)#with respect to r
        ref=Tensor.Tensor(sp,["cellar"],{(0,):0,(1,):-1/r,(2,):0})
        self.assertEqual(drv,ref)
        
        
        


if  __name__ == '__main__':
    unittest.main()
