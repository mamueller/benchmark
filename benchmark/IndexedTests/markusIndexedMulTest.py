#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from MarkusIndexed import VIB, OIB, VI ,VectorFieldBase,OneFormFieldBase, PatchWithMetric

from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
from sympy.tensor import  Idx
from sympy import eye,zeros,Symbol, symbols, Lambda, diff , Derivative
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
        
        n=3
        m = Manifold('M', n)
        patch = PatchWithMetric('P', m)
        cart= CoordSystem('cart', patch)
        bc=VectorFieldBase("bc",cart)
        br=OneFormFieldBase(bc)
        x=VIB("x",[bc,br,bc,br])
        with self.assertRaises(IncompatibleShapeException):
            x[0,0]=3
    def test_getWithContraction(self):
        n=3
        m = Manifold('M', n)
        patch = PatchWithMetric('P', m)
        cart= CoordSystem('cart', patch)
        bc=VectorFieldBase("bc",cart)
        br=OneFormFieldBase(bc)
        #full trace
        x=VIB("x",[bc,br])
        x[0,0]=3
        x[1,1]=4
        i=Idx('i')
        self.assertEqual(x[i,i],7)
        # partial trace
        x=VIB("x",[bc,br,bc,br])
        y=VIB("y")
        x[0,0,0,0]=3
        x[1,1,0,0]=4
        i, j, k       = map(Idx, ['i', 'j','k'])
        print("before test")
        y[j,k]=x[i,i,j,k]
        self.assertEqual(y.bases,[bc,br])
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
        n=2
        m = Manifold('M', n)
        patch = PatchWithMetric('P', m)
        cart= CoordSystem('cart', patch)
        ## cases with contraction
        bc=VectorFieldBase("cart_st",cart)
        br=OneFormFieldBase(bc)
        res=VIB("res")
        x=VIB("x",[br])
        A=VIB("A",[bc,br])
        i, j, k,l        = map(Idx, ['i', 'j', 'k','l'])
        x[0]=3
        A[0,0]=2
        print("befor mult")
        res[j]= A[i,j]*x[i]
        print("after mult")
        print("res.bases")
        print(res.bases)
        self.assertEqual(res[0],6)
        self.assertEqual(res.bases,[br])
        
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

 #   def test_del(self):
 #       #sp=Spherical()
 #       i, j, k,l        = map(Idx, ['i', 'j', 'k','l'])
 #       bc=VectorFieldBase(name,cart)
 #       br=OneFormFieldBase(bc)
 #       r,phi,theta=symbols("r phi theta")
 #       f=Symbol("f")
 #       x=VIB("x",[bc])
 #       nabla=OIB("nabla",[br])
 #       # note that the partial derivative operator is overloaded
 #       # for VI objects
 #       nabla[0]=Lambda(f,diff(f,r))
 #       nabla[1]=Lambda(f,diff(f,phi))
 #       nabla[2]=Lambda(f,diff(f,theta))
 #       #y=VIB("y")
 #       #y[i,j]=nabla[i]*x[j]
 #   ###def test_div(self):

 #   ##    y=delop[j]*x[j]
 #   ##    y[i]=delop[j]*A[j,i]
    #def VID2sumOfdyads(self):
    #    # sum of base vectors
    #    bc=VectorFieldBase(name,cart)
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
        # cellar base
        bc=VectorFieldBase("bc",cart)
        # roof base
        br=OneFormFieldBase(bc)
        g=VIB("g",[bc,bc])
        for i in range(n):
            g[i,i]=1
        i=Idx("i")
        j=Idx("j")
        patch.setMetric("cart",g[i,j])
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
        # note that we not only need the partial derivative of the components 
        # but also the partial derivatives of the base vectors expressed as linear combinations
        # of those base vectors which are  expressed by the cristoffel symbols which in turn can be computed from
        # known cristoffel symbols of a given coordinate system and the connection 
        # between the coordinates 
        # 
        # A case of special interest is the cartesian coordinate system whose 
        # coordinate functions are given as the partial derivatives of functions 
        # on the manifold in the directions of the set of base vectors e_x,e_y,ez
        # We therefore start with the cartesian base or equivalently 
        # with the cartesian coordinate system defined by this vectorfield 


        #get the catesian cellar base vectors
        dim=3
        manyf = Manifold('M', dim)
        patch = PatchWithMetric('P', manyf)
        cart= CoordSystem('cart', patch)
        cc=VectorFieldBase("cc",cart)
        g=VIB("g",[cc,cc])
        for i in range(dim):
            g[i,i]=1
        i=Idx("i")
        j=Idx("j")
        patch.setMetric("cart",g[i,j])
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
        bc=VectorFieldBase("bc",spherical)
        # roof base
        br=OneFormFieldBase(bc)
        x=VIB("x",[bc])
        r=Symbol("r")
        x[0]=r**2
        i=Idx("i")
        print(type(x[i]))
        #res=Derivative(x[i],r,evaluate=True)
        res=diff(x[i],r)
        print((res))
        
    def test_change_of_base(self):
        # assume that between two cellar bases {B_j,j=1..n} and {b_i ,i=1..n}
        # the following linear mapping exists
        #                       i 
        # B[j]=A[i,j]b[i] B  = A   b
        #                  j     j  i
        #                                             j 
        # The indices V wrt the new reciprocal base B
        #              j                                                          i
        # of a vector represented by its indices v  wrt the old reciprocal base b 
        #                                         i             
        # can be represented by
        #       i  
        #V   = A  v      ( simmonds)
        # j     j  i
        #                                                                  k
        # We also can represent this by a multiplication of the vector v  b
        #                                                               k
        #                            i     j    
        # with an identity tensor I_A  b  B  
        #                            j  i     
        # 
        # the result is:                             
        #
        #    i  j        k        i  j        k       i      j
        # I_A  B  b  v  b    = I_A  B  v  b  b   = I_A   v  B      
        #    j     i  k           j     k  i          j   i       
        # so the components of I_A and A are the same.
        n=2
        m = Manifold('M', n)
        P = PatchWithMetric('P', m)
        cart= CoordSystem('cart', P)
        # cellar base
        bc=VectorFieldBase("bc",cart)
        # roof base
        br=OneFormFieldBase(bc)
        # metric
        # We define the components of g wrt the roof base
        # therefore the bilinear form is applicable directly to cellar base 
        # vectors and since g induces a scalar product our components
        # should be chosen as the scalar product of cellar base vectors.
        # this is equivalent to the definition in exercise 2.8 on page 39 of Simmonds
        # nomenclature g_cc means cellar cellar components, which means components with respect to roof roof bases)
        g_cc=VIB("IAB",[br,br])
        # to make thing visible we chose some symbols
        g00,g11=symbols("g00,g11")
        g_cc[0,0]=g00
        g_cc[1,1]=g11
        i, j, k        = map(Idx, ['i', 'j','k'])
        P.setMetric("cart",g_cc[i,j])
        # new cellar base
        Bc=VectorFieldBase("Bc")
        # new roof base
        Br=OneFormFieldBase(bc)
        # define the transformation
        # nomeclature: the name rc reflects that the first index is a roof and the second a cellar index
        IA_rc=VIB("IA_bc_Br",[bc,Br])
        a,b=symbols("a,b")
        IA_rc[0,0]=a
        IA_rc[1,1]=b
        x_c=VIB("x",[br])
        x1,x2=symbols("x1,x2")
        x_c[0]=x1
        x_c[1]=x2
        # Apply the Identity Tranformation.
        y_c=VIB("y")
        y_c[j]=IA_rc[i,j]*x_c[i]
        self.assertEqual(y_c.bases,[Br])
        self.assertEqual(y_c[0],a*x1)
        self.assertEqual(y_c[1],b*x2)
        #raise("not finished")
        # a possible alternative would be to register the Transformation
        # Bc.connect_to(bc,IA_bc_Br[i,j])
        # from now on finding the components wrt a registered base could be
        # done x[j].change2bases([Br])
        # self.assertEqual(x[0],a*x1)
        # self.assertEqual(x[2],b*x2)

        # now we want to transform the roof components of a vector (given wrt  the old cellar base)
        # to roof components wrt to the new cellar base
        x_r=VIB("x",[bc])
        x1,x2=symbols("x1,x2")
        x_r[0]=x1
        x_r[1]=x2
        # to be able to apply our Identity transformation we have to lower the first index and raise the second.
        # We achieve this by two other transformations which are called the 
        # musicial isomorphisms # (for raising) and b (for lowering)  which also can be formally described
        # as a tensor multiplications with the metric tensor represented wrt the appropriate bases
        # first we transform the roof components x_r back to x_c which we already can handle 
        # we need the representation of the metric tensor wrt the old cellar base vectors
        #g_rr=P.getMetric([bc,bc])
        # we allready have the metric wrt the roof base vectors
        mat_cc=zeros(n,n)
        dat=g_cc.data
        for tup in dat.keys():
            mat_cc[tup]=dat[tup]
        mat_rr=mat_cc.inverse_LU()   
        g_rr=VIB("IAB",[bc,bc])
        for i in range(n):
            for j in range(n):
                g_rr.data[(i,j)]=mat_rr[i,j]
                
        #g_rr=secondRankTensor(g_cc[i,j].to_matrix().inverseLU(),[bc,bc])
        IA_rr=VIB("IA_rr")
        IA_rr[i,j]=IA_rc[i,k]*g_rr[k,j]
        #y_r[j]=IA_br_Bc[i,j]*x[i]
        #self.assertEqual(y.bases,[Br])
        #self.assertEqual(y[0],1./a*x1)
        #self.assertEqual(y[2],1./b*x2)





        
if  __name__ == '__main__':
    unittest.main()
