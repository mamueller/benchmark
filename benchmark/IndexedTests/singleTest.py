#!/usr/bin/python3
# vim: set expandtab ts=4
from ChangeOfBaseTests import *

def suite():
    s=unittest.TestSuite()
    #s.addTest(IndexedTest("test_getWithContraction"))
   # s.addTest(IndexedTest("test_mult"))
    s.addTest(TestChangeOfBaseGivenByCellarVectorTransform("test_oldRoofToNewRoofSecondOrder"))
    #s.addTest(IndexedTest("test_getSetItems"))
    #s.addTest(IndexedTest("test_getMultipleMixedIntegerAndSympolicIndices"))
    #s.addTest(IndexedTest("test_partial_derivative_of_VI_wrt_a_coordinate"))
    #s.addTest(IndexedTest("test_free_symbols"))
    return(s)

if  __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
