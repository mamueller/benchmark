#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
#import SolutionSuite
from OperatorTest import *
#from ChandraTest import *
#from CylinderTest import *
#from SplineTest import *
from TensorTest import *

def suite():
    s=unittest.TestSuite()
    #s.addTest(OperatorTest("test_grad_v_sym"))
    #s.addTest(OperatorTest("test_roofCellar_transpose"))
    #s.addTest(TensorTest("test_transform2_secondOrderTensors"))
    s.addTest(OperatorTest("test_roofCellar_transpose"))
    #s.addTest(TensorTest("test_vec_grad"))
    #s.addTest(TensorTest("test_partder"))
   #  s.addTest(TensorTest("test_Tensor_initialization"))
   # s.addTest(TensorTest("test_scalarProduct"))
    #s.addTest(TensorTest("test_Vector_scalarProduct"))
    #s.addTest(TensorTest("test_transform2_vector"))
    #s.addTest(TensorTest("test_extractVector"))
    #s.addTest(TensorTest("test_toVector"))
    #s.addTest(TensorTest("test_Vector_initialization"))
    #s.addTest(TensorTest("test_del_index"))
    #s.addTest(TensorTest("test_add_index"))
    #s.addTest(TensorTest("test_raise_and_lower"))
    #s.addTest(TensorTest("test_equal"))
    #s.addTest(TensorTest("test_cart"))
    #s.addTest(TensorTest("testScalarMul"))
    #s.addTest(TensorTest("testAdd"))
    #s.addTest(TensorTest("testSub"))
    return(s)

if  __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
