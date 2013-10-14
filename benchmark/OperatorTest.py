#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
from Spherical import Spherical 
from Cartesian import Cartesian
import Tensor
import copy
import Exceptions
from helperFunctions import pp

class OperatorTest(unittest.TestCase):
    def setUp(self):
        # for all tests we create 2 instances of a general coordinate system
	# a cartesian and a spherical one
        
        
	self.C=Cartesian()
	self.S=Spherical()



###########################################################
    def test_grad_v_sym(self):
        sp=Spherical()
        r,phi,theta=sp.U
        u_r,u_phi,u_theta=symbols("u_r,u_phi,u_theta")
	# we define a symbolic velocity field by its phyical components
        u_phys=Matrix(3,1,[u_r(r,phi,theta),u_phi(r,phi,theta),u_theta(r,phi,theta)])
	# we translate it to roof components
        u_roof=sp.phys2roof(u_phys)
	print("u_roof="+str(u_roof))
	du_roof_dr=sp.part_der_v(0,u_roof) 
	print("u_roof,r="+str(du_roof_dr))

        

    def test_metricTensor(self):
        # we know that cartesian coordinates are orthogonal this 
        # means that the scalar product of a base vector with another 
        # base vector has to be zero
        C=self.C
        # first we attend to the cellar basis
        gcc=C.g_cc()
        for i in range(0,3):
            for j in range(0,3):
                if(i!=j):
                    self.assertTrue(gcc[i,j]==0)

        # now the cellar basis
        grr=C.g_rr()
        for i in range(0,3):
            for j in range(0,3):
                if(i!=j):
                    self.assertTrue(grr[i,j]==0)
        # since g_cc =J^T*J and g_rr Jinv*Jinv^T gcc*grr should be the identety matrix            
        self.assertTrue(grr*gcc==eye(3))
       
        # we know that spherical coordinates are also orthogonal this 
        # means that the scalar product of a base vector with another 
        # base vector has to be zero
        sp=self.S
        # first we attend to the cellar basis
        gcc=sp.g_cc()
        for i in range(0,3):
            for j in range(0,3):
                if(i!=j):
                    self.assertTrue(gcc[i,j]==0)

        # now the cellar basis
        grr=sp.g_rr()
        for i in range(0,3):
            for j in range(0,3):
                if(i!=j):
                    self.assertTrue(grr[i,j]==0)
        # since g_cc =J^T*J and g_rr Jinv*Jinv^T gcc*grr should be the identety matrix            
        self.assertTrue(grr*gcc==eye(3))

    def test_roof_cart_transformations(self):
        sp=self.S
        x,y,z=sp.X
        xu,yu,zu=sp.XofU
        r,phi,theta=sp.U
        # we start with vectors
        #                 i
        # suppose that v=v g
        #                   i
        # is given by its roof components
        vs=[Matrix(3,1,[1,0,0]),Matrix(3,1,[0,1,0]),Matrix(3,1,[0,0,1])]  
        for i in range(len(vs)):
            v_r=vs[i] 
            # we compute the  cartesian components
            v_cart=sp.roof2cart(v_r)
            # which should match the cartesian components of the appropriate
            # cellar base vector g
            #                     r
            self.assertTrue(v_cart==sp.gc[i])
            # now we can backtransform v_cart to roof components
            v_r_new=sp.cart2roof(v_cart)
            self.assertTrue(v_r_new==v_r)
    def test_cellar_cart_transformations(self):
        sp=self.S
        x,y,z=sp.X
        xu,yu,zu=sp.XofU
        r,phi,theta=sp.U
        # we start with vectors
        #                   i
        # suppose that v=v g
        #                 i  
        # is given by its roof components
        # now test the transformations from cellar to cart components and back
        vs=[Matrix(3,1,[1,0,0]),Matrix(3,1,[0,1,0]),Matrix(3,1,[0,0,1])]  
        for i in range(len(vs)):
            v_c=vs[i] 
            # we compute the  cartesian components
            v_cart=sp.cellar2cart(v_c)
            # which should match the cartesian components of the appropriate
            #                    i
            # roof base vector g
            #                     
            self.assertEqual(v_cart,sp.gr[i])
            # now we can backtransform v_cart to roof components
            v_c_new=sp.cart2cellar(v_cart)
            self.assertEqual(v_c_new,v_c)

#    def test_roof_cart_tensor_transformations(self):
#        sp=self.S
#        x,y,z=sp.X
#        xu,yu,zu=sp.XofU
#        r,phi,theta=sp.U
#        ## we first define a linear mapping P_x that acts as a projector of any vector to the 
#        # x axis then we can express the application of P_x to v by the dot
#        # product *
#        #      P_x =  e  e * v
#        #              x  x
#        #
#        # assume that we have a vector v given by its cartesian components
#        #    x    y     z
#        # v=v e +v e + v e
#        #      x    y     z
#           ##                
#        #                                        x
#        # When we apply P_x to v we get P_x(v) =v e
#        #                                          x
#        # The cartesian components of P_x = 1*ex ex + 0*ex ey +0 ...
#        #
#        # Thus if v is expressed in cartesian components 
#        # the result P_x(v) expressed in cartesian coordinates is given by the
#        # matrix multiplication res=P_x v
#        Pcartcart=Matrix([[0,1,0],[0,0,0],[0,0,0]])
#        # we chose the cartesian components of the first roof base vector
#        v_cart=sp.gr[0] 
#        res_cart=Pcartcart*v_cart
#        
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
    def test_roof_cellar_transformations(self):
        sp=self.S
        x,y,z=sp.X
        xu,yu,zu=sp.XofU
        r,phi,theta=sp.U
        # we now attend to the cellar components
        # since (in our case of spherical case) the roof base vectors are parallel to  
        # the first cellar base vectors we expect to see the scale factors
        # suared
        self.assertTrue(sp.roof2cellar(Matrix(3,1,[1,0,0]))==Matrix(3,1,[1,0,0]))
        # 
        self.assertTrue(sp.roof2cellar(Matrix(3,1,[0,1,0]))==Matrix(3,1,[0,r**2*sin(theta)**2,0]))
        
        self.assertTrue(sp.roof2cellar(Matrix(3,1,[0,0,1]))==Matrix(3,1,[0,0,r**2]))
	#

    def test_normal_stress(self):
        # we compute the normal viscous stress at a spherical surface
        # It is related to the velocity by the definition
        # S=mu(nabla v +(nabla v)^t)
        # Assume that we have the components s_cr=s_i,^k  
        # of the stresstensor with respect to the dyadic products
        # g^ig_k. Those are the "cellar roof " components.

        # s_cr[1,1]*g^1g_1 + s_cr[1,2]*g^1g_2 + s_cr[1,3]*g^1g_3
        # s_cr[2,1]*g^2g_1 + s_cr[2,2]*g^2g_2 + s_cr[2,3]*g^2g_3
        # s_cr[3,1]*g^3g_1 + s_cr[3,2]*g^3g_2 + s_cr[3,3]*g^3g_3

        # Assume further that we have the cellar components n_c of the 
        # normal vector n  (given with respect to the roof basis g^k)
        # n_c=[1,0,0] -> n=1*g^1
        
        # using the fact that the scalar products <g_i,g^j>=delta_i^j
        # we can then compute the normal stress  vector sn=S n
        # sn=S n_c[1]g^1+ n_c[2] g^2 + n_c[3] g^3 
        # =(s_cr[1,1]n_c[1]+s_cr[1,2]n_c[2]+s_cr[1,3]n_c[3])g^1
        # +(s_cr[2,1]n_c[1]+s_cr[2,2]n_c[2]+s_cr[2,3]n_c[3])g^2
        # +(s_cr[3,1]n_c[1]+s_cr[3,2]n_c[2]+s_cr[3,3]n_c[3])g^3

        # thus we can express the cellar components sn_c of sn 
        # by a the matrix product sn_c=s_cr*n_c)

        print("bla")
    def test_roofCellar_transpose(self):
        sp=self.S
        # we check the definition of the transposed tensor
        vr,vphi,vtheta,ur,uphi,utheta=symbols("vr,vphi,vtheta,ur,uphi,utheta",real=True)
	T11,T12,T13,T21,T22,T23,T31,T32,T33=symbols("T11,T12,T13,T21,T22,T23,T31,T32,T33",real=True)
	T=Tensor(sp,["cellar","roof"],{(1,1):T11,(1,2):T12,(1,3):T13,(2,1):T21,(2,2):T22,(2,3):T23,(3,1):T31,(3,2):T32,(3,3):T33})
	print("T=",T)
        #u_c=Matrix(3,1,[ur,uphi,utheta])
	u=Tensor(sp,["cellar"],{(1,):ur})
	v=Tensor(sp,["roof"],{(2,):vphi})
        res1=v|(T|u)
        ## the definition of the transposed Tensor
	    ## for all u,v : v | A |u =u |A^t| v
        TT=T.transpose()
        res2=u|(TT|v)
        print("res1-res2=",sp.scalarSimp(res1-res2))



if  __name__ == '__main__':
     unittest.main()
