#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from MarkusIndexed import VIB, OIB, VI ,IncompatibleShapeException,VectorFieldBase,OneFormFieldBase
from Bases import PatchWithMetric

from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
from sympy.tensor import  Idx
from sympy import eye,Symbol, symbols, Lambda, diff , Derivative
from sympy import pi,sin,cos,trigsimp,symbols,Abs,simplify
from sympy.diffgeom import Manifold, CoordSystem


class IndexedTest(unittest.TestCase):
    def test_getSetSingleIntegerIndex(self):
        x=VIB("x")
        x[0]=3
        self.assertEqual(x[0],3)
        
    def test_getSetMultipleIntegerIndeces(self):
        x=VIB("x")
        x[0,0]=3
        self.assertEqual(x[0,0],3)
        #

    def test_getSetSingleSympolicIndex(self):
        x=VIB("x")
        y=VIB("y")
        x[0]=3
        i, j        = map(Idx, ['i', 'j'])
        y[i]=x[j]
        self.assertEqual(y[0],3)
    
    def test_getMultipleSymbolicIndices(self):
        x=VIB("x")
        y=VIB("y")
        x[0,1]=3
        i, j        = map(Idx, ['i', 'j'])
        y[i,j]=x[i,j]
        self.assertEqual(y[0,1],3)
        
        y[i,j]=x[j,i]
        self.assertEqual(y[1,0],3)
    
    def test_getMultipleMixedIntegerAndSympolicIndices(self):
        x=VIB("x")
        z=VIB("z")
        x[0,0,1]=1
        x[1,0,1]=2
        x[0,1,1]=3
        x[1,1,1]=4
        x[0,0,2]=10
        x[1,0,2]=20
        x[0,1,2]=30
        x[1,1,2]=40
        i, j       = map(Idx, ['i', 'j'])
        z[i,j]=x[i,j,1]
        self.assertEqual(z[0,0],1)
        self.assertEqual(z[1,0],2)
        self.assertEqual(z[0,1],3)
        self.assertEqual(z[1,1],4)
    def test_setMultipleMixedIntegerAndSympolicIndices(self):
        x=VIB("x")
        y=VIB("y")
        z=VIB("z")
        x[0,0,1]=1
        x[1,0,1]=2
        x[0,1,1]=3
        x[1,1,1]=4
        i, j       = map(Idx, ['i', 'j'])
        y[1,i,j]=x[j,i,1]
        self.assertEqual(y[1,0,0],1)
        self.assertEqual(y[1,0,1],2)
        self.assertEqual(y[1,1,0],3)
        self.assertEqual(y[1,1,1],4)
    def test_getWithWrongShape(self):
        x=VIB("x")
        x[1]=1
        i, j       = map(Idx, ['i', 'j'])
        with self.assertRaises(IncompatibleShapeException):
            x[i,j]
    def test_setWithWrongShape(self):
        x=VIB("x")
        x[1]=1
        with self.assertRaises(IncompatibleShapeException):
            x[1,1]=1
        
        bc=VectorFieldBase()
        br=OneFormFieldBase(bc)
        x=VIB("x",[bc,br,bc,br])
        with self.assertRaises(IncompatibleShapeException):
            x[0,0]=3
    def test_getWithContraction(self):
        bc=VectorFieldBase()
        br=OneFormFieldBase(bc)
        x=VIB("x",[bc,br])
        x[0,0]=3
        x[1,1]=4
        i=Idx('i')
        self.assertEqual(x[i,i],7)
        
        x=VIB("x",[bc,br,bc,br])
        y=VIB("y")
        x[0,0,0,0]=3
        x[1,1,0,0]=4
        i, j, k       = map(Idx, ['i', 'j','k'])
        y[j,k]=x[i,i,j,k]
        self.assertEqual(y[0,0],7)
        self.assertEqual( x[i,i,j,j],7)
        
        x=VIB("x",[bc,br,bc,br,bc])
        y=VIB("y")
        x[0,0,0,0,0]=3
        x[1,1,0,0,0]=4
        i, j, k, l       = map(Idx, ['i', 'j','k','l'])
        z=x[i,i,j,j,0]
        self.assertEqual(z,7)

    def test_mult(self):
        ## cases with contraction
        bc=VectorFieldBase()
        br=OneFormFieldBase(bc)
        res=VIB("res")
        x=VIB("x",[br])
        A=VIB("A",[bc,br])
        i, j, k,l        = map(Idx, ['i', 'j', 'k','l'])
        x[0]=3
        A[0,0]=2
        res[j]= A[i,j]*x[i]
        self.assertEqual(res[0],6)
        
        with self.assertRaises(ContractionIncompatibleBaseException):
            ## bases are not compatible
            res[j]= A[i,j]*x[j]

        x=VIB("x",[br])
        A=VIB("A",[bc,br])
        x[0]=10
        x[1]=20
        A[0,0]=1
        A[1,0]=2
        A[0,1]=3
        A[1,1]=4
        res[j]= A[i,j]*x[i]
        self.assertEqual(res[0],50)
        self.assertEqual(res[1],110)
        
        ## cases without contraction
        res=VIB("res")
        res[i,j,k]= A[i,j]*x[k]
        self.assertEqual(res[0,0,0],10)
        
        res=VIB("res")
        with self.assertRaises(IncompatibleShapeException):
            res[j]= A[i,j]*x[k]
        
        res=VIB("res")
        x=VIB("x",[br])
        A=VIB("A",[bc,br,bc])
        with self.assertRaises(IncompatibleShapeException):
            res[i,j,l]= A[i,j,k]*x[k]

    def test_del(self):
        #sp=Spherical()
        i, j, k,l        = map(Idx, ['i', 'j', 'k','l'])
        bc=VectorFieldBase()
        br=OneFormFieldBase(bc)
        r,phi,theta=symbols("r phi theta")
        f=Symbol("f")
        x=VIB("x",[bc])
        nabla=OIB("nabla",[br])
        # note that the partial derivative operator is overloaded
        # for VI objects
        nabla[0]=Lambda(f,diff(f,r))
        nabla[1]=Lambda(f,diff(f,phi))
        nabla[2]=Lambda(f,diff(f,theta))
        #y=VIB("y")
        #y[i,j]=nabla[i]*x[j]
    ###def test_div(self):

    ##    y=delop[j]*x[j]
    ##    y[i]=delop[j]*A[j,i]
    #def VID2sumOfdyads(self):
    #    # sum of base vectors
    #    bc=VectorFieldBase()
    #    br=OneFormFieldBase(bc)
    #    x=VIB("x",[br])
    #    r=Symbol("r")
    #    x[0]=r
    #    i=Idx("i")
    #    res=x[i].base_dyads()
    #    d=
    #    ref=r*
        
    #def test_partial_derivative_of_a_baseVector(self):
    #def test_partial_derivative_of_a_Dyad_of_base_vectors(self):
    def test_free_symbols(self):
        n=3
        m = Manifold('M', n)
        patch = PatchWithMetric('P', m)
        cart= CoordSystem('cart', patch)
        patch.setMetrix("cart",eye(n))
        # cellar base
        bc=VectorFieldBase(cart)
        # roof base
        br=OneFormFieldBase(cart)
        x=VIB("x",[bc])
        r,phi=symbols("r,phi")
        x[0]=r**2
        x[1]=phi**2
        i=Idx("i")
        res=x[i].free_symbols
        self.assertEqual(res,set({phi,r}))
        x=VIB("x",[br])
        r,phi=symbols("r,phi")
        x[0]=r**2
        i=Idx("i")
        res=x[i].free_symbols
        self.assertEqual(res,set({r}))
        

    def test_partial_derivative_of_VI_wrt_a_coordinate(self):
        # note that we not even need the partial derivative of the components 
        # but also of the partial derivative of the base vectors which are 
        # expressed by the cristoffel symbols which in turn can be computed from
        # known cristoffel symbols of a given coordinate system and the connection 
        # between the coordinates 
        # 
        # A case of special interest is the cartesian coordinate system whose 
        # coordinate functions are given as the partial derivatives of functions 
        # on the manifold in the directions of the set of base vectors e_x,e_y,ez
        # We therefore start with the cartesian base or equivalently 
        # with the cartesian coordinatesystem defined by this vectorfield 


        #get the catesian cellar base vectors
        dim=3
        manyf = Manifold('M', dim)
        patch = PatchWithMetric('P', manyf)
        cart= CoordSystem('cartesian', patch)
        cc=VectorFieldBase(cart)
        g=VIB("g",[cc,cc])
        for i in range(dim):
            g[i,i]=1
        i=Idx("i")
        j=Idx("j")
        patch.setMetric("cartesian",g[i,j])
        spherical= CoordSystem('spherical', patch)
        # now connect the two coordsystems by a transformation
        x,y,z=symbols("x,y,z")
        r=symbols("r",positive=True,real=True)
        phi,theta=symbols("phi,theta",real=True)
        spherical.connect_to(cart,
            [r,
             phi,
             theta
            ]
            ,
            [r*cos(phi)*sin(theta) ,
             r*sin(phi)*sin(theta),
             r*cos(theta)
            ],
            inverse=False
        )
        # cellar base
        bc=VectorFieldBase(spherical)
        # roof base
        br=OneFormFieldBase(spherical)
        x=VIB("x",[bc])
        r=Symbol("r")
        x[0]=r**2
        i=Idx("i")
        print(type(x[i]))
        #res=Derivative(x[i],r,evaluate=True)
        res=diff(x[i],r)
        print((res))
        

        
if  __name__ == '__main__':
    unittest.main()
