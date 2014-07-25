#!/usr/bin/python3
from __future__ import print_function, division
from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from sympy import symbols, default_sort_key,Expr
from sympy.tensor import IndexedBase, Idx, Indexed
from helperFunctions import pp
from sympy.core import Expr, Tuple, Symbol, sympify, S
from Exceptions import BaseMisMatchExeption,DualBaseExeption

class VIB(IndexedBase):
    def __getitem__(self, indices, **kw_args):
            return VI(self, indices, **kw_args)
class VI(Indexed):
    pass
class ViExpr(VI):
    _op_priority = 11.0
    is_commutative = False
class ViMul(ViExpr):
    def __new__(cls, coeff, *args, **kw_args):
        print("this is ViMul")
        print(*args)
        #coeff = sympify(coeff)
        base=VIB("X")
        obj=VI.__new__(cls,base,*args,**kw_args)
        return obj

def tt():
    A=VIB("A")
    B=VIB("B")
    i=Idx('i',2)
    res=A[i]*B[i]
    print(res)
tt()
