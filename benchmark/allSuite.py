#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
import SolutionSuite
from OperatorTest import *
from ChandraTest import *
from CylinderTest import *
from SplineTest import *
from TensorTest import *

def suite():
    s=SolutionSuite.suite()
    cases=[ChandraTest,CylinderTest,OperatorTest,SplineTest,TensorTest]
    for c in cases:
        s.addTests(unittest.TestLoader().loadTestsFromTestCase(c))
    return(s)

if  __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
