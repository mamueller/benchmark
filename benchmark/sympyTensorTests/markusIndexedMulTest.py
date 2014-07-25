from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable
#
class VIB(IndexedBase):
    def __getitem__(self, indices, **kw_args):
        #possible shortcut
        #I=super(VIB,self).__getitem__(indices, **kw_args)
        #VI=VI_from_Indexed(I) ? not yet implemented
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
    def __new__(cls, base, *args, **kw_args):
        from sympy.utilities.misc import filldedent
        if not args:
            raise IndexException("Indexed needs at least one index.")
        if isinstance(base, (string_types, Symbol)):
            base = VIB(base)
        elif not isinstance(base, VIB):
            raise TypeError(filldedent("""
                Indexed expects string, Symbol or VIB as base."""))

        args = list(map(sympify, args))
        return Expr.__new__(cls, base, *args, **kw_args)

    
    def __mul__(self, other):
        print("in mul")
        print(self)
        ## we now create an Indexed object and can then use all the methods 
        expr=super(VI,self).__mul__(other)
        print(type(expr)) 
        print(expr) 
        print(get_contraction_structure(expr))
        print(get_indices(expr))
        # now we can uses this information to implement the the multiplication for our subtype.
        

x, A = map(VIB, ['x', 'A'])
i, j = map(Idx, ['i', 'j'])
print(type(x[i]))
A[i,j]*x[i]
#

