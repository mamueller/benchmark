
from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable

def permuteTuple(t,newPositions):
    nl=[t[p] for p in newPositions]
    return(tuple(nl))



class VIB(IndexedBase):
    def __new__(cls,label, shape=None, **kw_args):
        obj=IndexedBase.__new__(cls, label, shape=None, **kw_args)
        obj.data=dict()
        return(obj)
    def __getitem__(self, indices, **kw_args):
        #possible shortcut
        #I=super(VIB,self).__getitem__(indices, **kw_args)
        #VI=VI_from_Indexed(I) ? not yet implemented
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch.")
            if type(indices[0]) is int:
                return(self.data[indices])
            else:    
                return VI(self, *indices, **kw_args)
        else:
            # only one index 
            index=indices 
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch.")
            if type(index) is int:
                return(self.data[index])
            else:    
                return VI(self, index, **kw_args)
    
    
    
    def __setitem__(self, indices, value,**kw_args):
        if  not(is_sequence(indices)): #only one index
            index=indices
            if type(index) is int:
                self.data[indices]=value
            elif type(index) is Idx:
                self.data=value.base.data
        else: #general case with more than one indices
            #if type(indices[0]) is int:
            if all(type(k) is int for k in indices):
                # all indices are 
                self.data[indices]=value
            
            elif all( type(i) is Idx for i in indices):
                # all indices are symbolic
                vb=value.base
                vbd=vb.data
                if set(indices)==set(value.indices):
                    if indices==value.indices:
                        self.data=vbd
                    else:
                        # a permutation of indices, necessiating a permutation of 
                        newPositions=[indices.index(i) for i in value.indices]
                        vKeys=vbd.keys()
                        for k in vKeys:
                            self.data[permuteTuple(k,newPositions)]=vbd[k]
                            





         

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
        

