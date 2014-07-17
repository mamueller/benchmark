from __future__ import print_function, division
from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from sympy import symbols, default_sort_key,Expr
from sympy.tensor import IndexedBase, Idx, Indexed
from helperFunctions import pp
from sympy.core import Expr, Tuple, Symbol, sympify, S
class VIB(IndexedBase):
    def __new__(cls, label,data,shape=None,**kw_args):
        if isinstance(label, string_types):
            label = Symbol(label)
        elif isinstance(label, Symbol):
            pass
        else:
            raise TypeError("Base label should be a string or Symbol.")
        
        obj = Expr.__new__(cls, label, **kw_args)
        obj.data=data
        if is_sequence(shape):
            obj._shape = Tuple(*shape)
        else:
            obj._shape = sympify(shape)

        return obj
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
        ds=self.base.data
        do=other.base.data
        
        return(4)

class Singleton(type):
    def __init__(cls, name, bases, dict):
        super(Singleton, cls).__init__(name, bases, dict)
        cls.instance = None 
        
    def __call__(cls,*args,**kw):
        if cls.instance is None:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
        return cls.instance
        
class VectorFieldBase(object):
    __metaclass__ = Singleton
    def __init__(self,partner=None):
        self.dual=dual
class OneFormFieldBase(VectorFieldBase):

class B1(VectorFieldBase):
    pass

def tt():
    b1=B1()
    i=Idx('i',2)
    x=VIB("x",{(0,):2,(1,):2})
    y=VIB("y",{(0,):2,(1,):2})
    x[i]
    res=x[i]*y[i]
    print(res)
    assert(res==8)
tt()
