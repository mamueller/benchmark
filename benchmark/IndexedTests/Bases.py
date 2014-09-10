#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
from sympy.diffgeom import Manifold, Patch, CoordSystem
class BaseVector(object):
   def __init__(self,coordSystem,ind):
        self.ind=ind
        self.coord=coordSystem
   def _eval_derivative(self):
        return(5)
   
##########################################################
class VectorFieldBase(object):
    """This class stores relationships between bases on the same patch.
    Every base b_i has one and only one dual (or reciprocal) base b^j
    with b^j(b_i)= <b^j,b_i>= delta_i^j . 
    """
    def __init__(self,coordSystem):
        self.coord=coordSystem
        self.vecs=[BaseVector(coordSystem,i) for i in range(coordSystem.dim)]

##########################################################
class OneFormFieldBase(VectorFieldBase):
    def __init__(self,coordSystem):
        self.oneforms=coordSystem.base_oneforms()
