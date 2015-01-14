#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from MarkusIndexed import TensorIndexSet, OIB, VI ,VectorFieldBase,OneFormFieldBase, PatchWithMetric
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException,IndexRangeException,MissingBaseException
from sympy.tensor import  Idx
from sympy import eye,zeros,Symbol, symbols, Lambda, diff , Derivative
from sympy import pi,sin,cos,trigsimp,symbols,Abs,simplify,Matrix
from sympy.diffgeom import Manifold, CoordSystem
from concurrencytest import ConcurrentTestSuite, fork_for_tests

class TestChangeOfBaseGivenByCellarVectorTransform(unittest.TestCase):

    def setUp(s):
        # we assume that we know the matrix connecting
        # the new cellar base vectors to the old cellar base vectors
        n=2
        m = Manifold('M', n)
        P = PatchWithMetric('P', m)
        cart= CoordSystem('cart', P)
        # oldcellar base
        s.bc=VectorFieldBase("bc",cart)
        # old roof base
        s.br=OneFormFieldBase(s.bc)
        # new cellar base
        s.Bc=VectorFieldBase("Bc")
        # new roof base
        s.Br=OneFormFieldBase(s.Bc)
        a,b=symbols("a,b")
        s.symbols=[a,b]
        M=Matrix([[a,0],[0,b]]) #columns contain representation of new 
        # base vectors w.r.t. old 
        s.Bc.connect(s.bc,M)

    def test_oldRoofToNewRoof(s):
        # old roof components to new roof components
        # assume that between two cellar bases {B_j,j=1..n} and {b_i ,i=1..n}
        # the following linear mapping exists
        #                       i 
        # B[j]=A[i,j]b[i] B  = A   b
        #                  j     j  i
        #              j                            
        # The indices V  w.r.t. the new  base B
        #                                      j
        #                                         i
        # of a vector represented by its indices v  wrt the old base b 
        #                                                             i             
        # can be represented by
        #  i          i  j
        # V   = (A^-1)  v      ( Simmonds p.37  eq. 2.36)
        #             j  
        #                                                               k
        # We also can represent this by a multiplication of the vector v  b
        #                                                                  k
        #                            i      j 
        # with an identity tensor   T   B  b  
        #                             j  i       
        # 
        # the result is:                             
        #
        # k       i      j  k         i      k  j        i   j        i
        #v  b  = T   B  b  v  b    = T   B  v  b  b   = T   v  B    =V  B 
        #    k     j  i        k       j  i        k      j     i        i
        #
        #                                     i        i
        # so by comparison the components    T  =(A^-1) 
        #                                     j        j
        x=TensorIndexSet("x",[s.bc])
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        y=TensorIndexSet("y")
        i=Idx("i")
        a,b=s.symbols
        y[i]=x[i].rebase2([s.Bc])
        s.assertEqual(y.bases,[s.Bc])
        s.assertEqual(y[0],x1/a)
        s.assertEqual(y[1],x2/b)
         
        # 1 manual alternative ##################################################
        
        # define the (Identety)transformation
        IA=TensorIndexSet("IA",[s.br,s.Bc])
        IA[0,0]=1./a
        IA[1,1]=1./b
        # Apply the Identity Tranformation.
        j=Idx("j")
        y[j]=IA[i,j]*x[i]
        s.assertEqual(y.bases,[s.Bc])
        s.assertEqual(simplify(y[0]-x1/a),0)
        s.assertEqual(simplify(y[1]-x2/b),0)
        
    def test_oldRoofToNewRoofSecondOrder(s):
        # old roof components to new roof components
        # assume that between two cellar bases {B_j,j=1..n} and {b_i ,i=1..n}
        # the following linear mapping exists
        #                           i 
        # Bn[j]=A[i,j]b[i]   Bn  = A   b
        #                      j    j   i
        #               i,j                            
        # The indices Tn     w.r.t. the new  base Bn
        #                                           l
        #                                                       i,j
        # of a second order tensor  represented by its indices t  wrt the old base b 
        #                                                                           k             
        # can be represented by
        #   i,j          i      j  k,p
        # Tn      = (A^-1) (A^-1) t      ( Simmonds p.37  eq. 2.36)
        #                k      p 
        #                                                                   k,l
        # We also can represent this by two multiplications of the tensor t     b  b
        #                                                                        k  l
        #                            i       j 
        # with an identity tensor   I   Bn  b  
        #                             j   i       
        # 
        # the result is:                             
        #
        #   k,p             i       m   j       n    k,p      
        # t     b  b      =I   Bn  b   I   Bn  b   t     b  b    
        #        k  p        m   i      *n   j            k  p
        #
        #
        #                   i       m   j      k,p        
        #                 =I   Bn  b   I     t     Bn  b     
        #                    m   i       p           j  k 
        #                                         i
        #
        #                   i   j      k,p          
        #                 =I   I     t    Bn   Bn      
        #                    k   p          i    j  
        #
        #                   i,j
        #                 =T   Bn   Bn      
        #                        i    j  
        #                                        
        #                                        
        #
        #                                     i        i
        # so by comparison the components    I  =(A^-1) 
        #                                     *j       *j
        #
        t=TensorIndexSet("t",[s.bc,s.bc])
        t_00,t_11=symbols("t_00,t_11")
        t[0,0]=t_00
        t[1,1]=t_11
        inter=TensorIndexSet("inter")
        T=TensorIndexSet("T")
        # 1 manual alternative #########################################
        # define the (Identety)transformation
        IA=TensorIndexSet("IA",[s.Bc,s.br])
        a,b=s.symbols
        IA[0,0]=1./a
        IA[1,1]=1./b
        # Apply the Identity Tranformation.
        i,j,k,p=map(Idx,["i","j","k","p"])
        inter[j,p]=IA[j,k]*t[k,p] #pos 0
        s.assertEqual(inter.bases,[s.Bc,s.bc])
        T[i,j]=IA[i,p]*inter[j,p]
        s.assertEqual(T.bases,[s.Bc,s.Bc])
        #s.assertEqual(T[0,0],t_00/a**2)
        #s.assertEqual(T[1,1],t_11/b**2)

        T[i,j]=t[i,j].rebase2([s.Bc,s.Bc])
        s.assertEqual(T.bases,[s.Bc,s.Bc])
        s.assertEqual(T[0,0],t_00/a**2)
        s.assertEqual(T[1,1],t_11/b**2)
         


    def test_oldCellarToNewCellar(s):
        # old cellar components to new cellar components
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
        #V   = A  v      ( Simmonds p.37  eq. 2.63)
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
        x=TensorIndexSet("x",[s.br])
        print("blub6")
        print(x.bases)
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        y=TensorIndexSet("y")
        i=Idx("i")
        a,b=s.symbols
        y[i]=x[i].rebase2([s.Br])
        s.assertEqual(y.bases,[s.Br])
        s.assertEqual(y[0],a*x1)
        s.assertEqual(y[1],b*x2)
         
        # manual alternative 
        # define the (Identity)transformation
        IA=TensorIndexSet("IA",[s.bc,s.Br])
        IA[0,0]=a
        IA[1,1]=b
        # Apply the identity tranformation.
        y=TensorIndexSet("y")
        j=Idx("j")
        y[j]=IA[i,j]*x[i]
        s.assertEqual(y.bases,[s.Br])
        s.assertEqual(y[0],a*x1)
        s.assertEqual(y[1],b*x2)


class TestChangeOfBaseWithMetric(unittest.TestCase):
    def test_oldRoofToNewCellar(s):
        ###################  manual alternative with metric #########################
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
        # there is already a test demonstrating this:
        n=2
        i, j, k        = map(Idx, ['i', 'j','k'])
        m = Manifold('M', n)
        P = PatchWithMetric('P', m)
        cart= CoordSystem('cart', P)
        # cellar base
        bc=VectorFieldBase("bc",cart)
        # roof base
        br=OneFormFieldBase(bc)
        # new cellar base
        Bc=VectorFieldBase("Bc")
        # new roof base
        Br=OneFormFieldBase(Bc)
        # define the (Identety)transformation
        IA=TensorIndexSet("IA",[bc,Br])
        a,b=symbols("a,b")
        IA[0,0]=a
        IA[1,1]=b
        x=TensorIndexSet("x",[br])
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        # Apply the Identity Tranformation.
        y=TensorIndexSet("y")
        y[j]=IA[i,j]*x[i]
        s.assertEqual(y.bases,[Br])
        s.assertEqual(y[0],a*x1)
        s.assertEqual(y[1],b*x2)

        # now we want to transform the roof components of a vector
        # to new cellar components 
        x=TensorIndexSet("x",[bc])
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        
        # to be able to apply our Identity transformation we have to lower the index of x
        # We achieve this by the musicial isomorphisms b (for lowering)  
        # which also can be formally described
        # as a tensor multiplications with the metric tensor represented wrt the 
        # the old cellar base vectors

        # We define the components of g wrt the roof base
        # therefore the bilinear form is applicable directly to cellar base 
        # vectors and since g induces a scalar product our components
        # should be chosen as the scalar product of cellar base vectors.
        # this is equivalent to the definition in exercise 2.8 on page 39 of Simmonds
        
        g=TensorIndexSet("g",[br,br])
        # to make thing visible we chose some symbols
        g00,g11=symbols("g00,g11")
        g[0,0]=g00
        g[1,1]=g11
        y[j]=IA[i,j]*g[i,k]*x[k]
        s.assertEqual(y.bases,[Br])
        raise(Exception("The result should be checked against something. E.g. a hand computed solution"))
        
     
        
class TestChangeOfBaseGivenByRoofVectorTransform(unittest.TestCase):

    def setUp(s):
        # we assume that we know the matrix connecting
        # the new roof base vectors to the old roof base vectors 
        # a connection between the reciprocal bases
        n=2
        m = Manifold('M', n)
        P = PatchWithMetric('P', m)
        cart= CoordSystem('cart', P)
        # oldcellar base
        s.bc=VectorFieldBase("bc",cart)
        # old roof base
        s.br=OneFormFieldBase(s.bc)
        # new cellar base
        s.Bc=VectorFieldBase("Bc")
        # new roof base
        s.Br=OneFormFieldBase(s.Bc)
        a,b=symbols("a,b")
        s.symbols=[a,b]
        M=Matrix([[a,0],[0,b]]) #columns contain representation of new 
        # base vectors w.r.t. old 
        s.Br.connect(s.br,M)

    def test_oldRoofToNewRoof(s):
        # this is equivalent to the case where the connection between the cellar bases is given an 
        # the cellar cellar transformation is sought
        # so the inverse Matrix will be needed
        x=TensorIndexSet("x",[s.br])
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        y=TensorIndexSet("y")
        i=Idx("i")
        a,b=s.symbols
        y[i]=x[i].rebase2([s.Br])
        s.assertEqual(y.bases,[s.Br])
        s.assertEqual(y[0],x1/a)
        s.assertEqual(y[1],x2/b)
         
        # manual alternative 
        
        # define the (Identety)transformation
        IA=TensorIndexSet("IA",[s.bc,s.Br])
        IA[0,0]=1./a
        IA[1,1]=1./b
        # Apply the Identity Tranformation.
        j=Idx("j")
        y[j]=IA[i,j]*x[i]
        s.assertEqual(y.bases,[s.Br])
        s.assertEqual(simplify(y[0]-x1/a),0)
        s.assertEqual(simplify(y[1]-x2/b),0)
     
    def test_oldCellarToNewCellar(s):
        # this is equivalent to the case where the connection between the cellar bases is given an 
        # the roof roof transformation is sought
        # so the original matrix can be used 
        x=TensorIndexSet("x",[s.bc])
        print("blub6")
        print(x.bases)
        x1,x2=symbols("x1,x2")
        x[0]=x1
        x[1]=x2
        y=TensorIndexSet("y")
        i=Idx("i")
        a,b=s.symbols
        y[i]=x[i].rebase2([s.Bc])
        s.assertEqual(y.bases,[s.Bc])
        s.assertEqual(y[0],a*x1)
        s.assertEqual(y[1],b*x2)
         
        # manual alternative 
        # define the (Identity)transformation
        IA=TensorIndexSet("IA",[s.br,s.Bc])
        IA[0,0]=a
        IA[1,1]=b
        # Apply the identity tranformation.
        y=TensorIndexSet("y")
        j=Idx("j")
        y[j]=IA[i,j]*x[i]
        s.assertEqual(y.bases,[s.Bc])
        s.assertEqual(y[0],a*x1)
        s.assertEqual(y[1],b*x2)

def suite():        
    suite = unittest.TestLoader().loadTestsFromTestCase(TestChangeOfBaseGivenByCellarVectorTransform)
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestChangeOfBaseGivenByRoofVectorTransform))
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestChangeOfBaseWithMetric))
    return(suite)

if  __name__ == '__main__':
    concurrent_suite = ConcurrentTestSuite(suite(), fork_for_tests(16))
    runner = unittest.TextTestRunner()
    runner.run(concurrent_suite)
