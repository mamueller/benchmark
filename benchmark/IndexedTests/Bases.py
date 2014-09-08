#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
##########################################################
class VectorFieldBase(object):
    """This class stores relationships between bases on the same patch.
    Every base b_i has one and only one dual (or reciprocal) base b^j
    with b^j(b_i)= <b^j,b_i>= delta_i^j . 
    This relationship is expressed by a reference every 
    base object holds its dual base."""
    @classmethod
    def checkDualType(cls,dual):
        if type(dual).__name__!="OneFormFieldBase":
            raise DualBaseExeption()

    def __init__(self,dual=None):
        if dual:
            type(self).checkDualType(dual)
            self.dual=dual
            #register self as the dual of the dual base
            dual.dual=self

##########################################################
class OneFormFieldBase(VectorFieldBase):
    @classmethod
    def checkDualType(cls,dual):
        if not(isinstance(dual,VectorFieldBase)):
            raise DualBaseExeption()

##########################################################
class CartesianVectorFieldBase(VectorFieldBase):
    def firstKindChristoffelSymbols(self):
        pass
