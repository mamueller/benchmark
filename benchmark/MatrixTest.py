#!/usr/bin/python
# vim: set expandtab ts=4
from Cylinder import *
import unittest 
from sympy import bspline_basis
from coords import *
from sympy import *
from sympy import pi
import mpmath as mp

class ChandraTest_2(unittest.TestCase):

    def test_matrix(self):
        nq=Matrix(2,2,[2,1,0,1])
        p= Matrix(2,2,[3,0,0,1])
        T=nq.inv()
        t=T.det()
        pT=p*T
        ev=pT.eigenvects()
        print("ev=",ev)
        c1=ev[0][0]
        v1=Matrix(ev[0][2])
        c2=ev[1][0]
        v2=Matrix(ev[1][2])
        m1=p-c1*nq
        m2=p-c2*nq
        print("determinants",m1.det(),m2.det())
        print("nullspaces",m1.nullspace(),m2.nullspace())


if  __name__ == '__main__':
     unittest.main()
