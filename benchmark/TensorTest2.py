#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import unittest 
from Spherical import Spherical 
from Cartesian import Cartesian
import Tensor2
import copy
import Exceptions
import sympy
from helperFunctions import pp
class TensorTest2(unittest.TestCase):
    def test_ChangeOfBase(self):
        # Note that a change of base for linear mappings 
        # is possible without a scalar product or a metric
        # Only the description 
        # of one base as linear combination of the other is needed
        # Actually the relationship to cartesian coordinates is not 
        # needed to transform the indices of a vector from one base to another 
        # A scalar product is however needed for the definition of a roof and cellar basis
        sp=Spherical()
        # as a first example we look at a vector
        v_src=Tensor2.Tensor2(sp,["src_cellar"],{(0,):1})
        mat=sympy.eye(3)*1./2.
        src2doubleSrc= Tensor2.ChangeOfBase("src_cellar","doubleSrc_cellar",mat)
        # The matrix columns describe the source cellar base vectors in terms of the 
        # (new) target basis (doubleSrc)  cellarbase vectors
        # Since the doubleSrc base vectors are twice as long the src base
        # vectors a src base vector is half the size of a doubleSrc vector
        # hence the 1/2 on the main diagonal
        v_doubleSrc=src2doubleSrc.transform(v_src,0)
        self.assertEqual(v_doubleSrc,Tensor2.Tensor2(sp,["doubleSrc_cellar"],{(0,):1./2.}))
        # 
        # now we look at a projection tensor which in (src) coordinates preserves the first component of a vector. 
        # (to help imagination we could assume the src base to have a fixed length of e.g. one although the mention of length is not necessary for the base change.)
        P_s=Tensor2.Tensor2(sp,["src_cellar","src_roof"],{(0,0):1})
        # first we observe its action on a vector
        v_s=Tensor2.Tensor2(sp,["src_cellar"],{(0,):1,(1,):2,(2,):3})
        v_sp=P_s(i,j)*v_s(j)

        # P_ds_s=src2doubleSrc.transform(P_s_s,0)
        # self.assertEqual(P_ds_s,Tensor2.Tensor2(sp,["doubleSrc","src"],{(0,0):1./2.}))
        # # now we exchage the other base
        # P_s_ds=src2doubleSrc.transform(P_s_s,1)
        # self.assertEqual(P_ds_s,Tensor2.Tensor2(sp,["doubleSrc","src"],{(0,0):1./2.}))
        # P_ds_ds=src2doubleSrc.transform(P_ds_s,1)
        # self.assertEqual(P_ds_ds,Tensor2.Tensor2(sp,["doubleSrc","doubleSrc"],{(0,0):1./4.}))
