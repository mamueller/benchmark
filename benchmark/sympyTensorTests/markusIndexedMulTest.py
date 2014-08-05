#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from MarkusIndexed import VIB, VI ,IncompatibleShapeException,VectorFieldBase,OneFormFieldBase

from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
from sympy.tensor import  Idx


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
        i, j, k        = map(Idx, ['i', 'j', 'k'])
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

        
if  __name__ == '__main__':
    unittest.main()
