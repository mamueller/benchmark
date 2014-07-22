from __future__ import print_function, division
from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from sympy import symbols, default_sort_key,Expr
from sympy.tensor import IndexedBase, Idx, Indexed
from helperFunctions import pp
from sympy.core import Expr, Tuple, Symbol, sympify, S
from Exceptions import BaseMisMatchExeption,DualBaseExeption
from Tensor import extract_keys
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

class OneFormFieldBase(VectorFieldBase):
    @classmethod
    def checkDualType(cls,dual):
        if type(dual).__name__!="VectorFieldBase":
            raise DualBaseExeption()

class VIB(IndexedBase):
    def __new__(cls, label,baseFieldList,components,shape=None,**kw_args):
        if isinstance(label, string_types):
            label = Symbol(label)
        elif isinstance(label, Symbol):
            pass
        else:
            raise TypeError("Base label should be a string or Symbol.")
        
        obj = Expr.__new__(cls, label, **kw_args)
        obj.baseFieldList=baseFieldList
        obj.components=components
        if is_sequence(shape):
            obj._shape = Tuple(*shape)
        else:
            obj._shape = sympify(shape)

        return obj
    def __setitem__(self,indices,**kwargs):

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch.")
            return VI(self, *indices, **kw_args)
        else:
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch.")
            return VI(self, indices, **kw_args)
class VI(Indexed):
    def __mul__(self,other):
        """multiplication is only allowed if the second factor
        is represented in the dual baseField of the first factor"""
        #extract the indexedBase objects from Base
        sIB=self.base
        sc=sIB.components
        sBFL=sIB.baseFieldList[0]
        
        oIB=other.base
        oc=oIB.components
        oBFL=oIB.baseFieldList[0]

        if sBFL!=oBFL.dual:
            raise BaseMisMatchExeption()
        sDummy=extract_keys(sc.keys(),0,0)
        oDummy=extract_keys(oc.keys(),0,0)
        print(sDummy)
        print(oDummy)
        nonZeroDummyKeys=set(sDummy).intersection(oDummy)
        print(nonZeroDummyKeys)
        res=sum(map(lambda key:sc[key]*oc[key],nonZeroDummyKeys))
        return(res)


def tt():
    bc=VectorFieldBase()
    br=OneFormFieldBase(bc)
    i=Idx('i',2)
    x=VIB("x",[bc],{(0,):2,(1,):2})
    y=VIB("y",[br],{(0,):2,(1,):2})
    res=x[i]*y[i]
    print(res)
    assert(res==8)
    j=Idx('j',2)
    A=VIB("A",[bc,bc],{(0,0):3})
    x
    #res=A[i,j]*y[i]
    #res[j]=A[i,j]*y[i]
    #assert(res.baseFieldList==[bc])
    #assert(res.components=={(0,):6})
tt()
