#!/usr/bin/python
# vim: set expandtab ts=4
from markusIndexedMulTest import *

def suite():
    s=unittest.TestSuite()
    #s.addTest(IndexedTest("test_diff"))
    s.addTest(IndexedTest("test_free_symbols"))
    return(s)

if  __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
