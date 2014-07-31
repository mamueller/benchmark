
from sympy.tensor.index_methods import get_contraction_structure, get_indices, _remove_repeated
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from copy import deepcopy

def permuteTuple(t,newPositions):
    nl=[t[p] for p in newPositions]
    return(tuple(nl))

def extractIndices(t,pos):
    nl=[t[p] for p in range(len(t)) if p in pos]
    return(tuple(nl))

def changedTuple(origTuple,newValTuple,newPositions):
    nl=list(origTuple)
    for p in range(len(newValTuple)):
        nl[newPositions[p]]=newValTuple[p]
    return(tuple(nl))


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
        if type(dual).__name__!="VectorFieldBase":
            raise DualBaseExeption()

##########################################################
class DualBaseExeption(BaseException):
    pass
##########################################################
class BaseMisMatchExeption(BaseException):
    pass
##########################################################
class IncompatibleShapeException(BaseException):
    pass

##########################################################
class VIB(IndexedBase):
    def __new__(cls,label,bases=None, **kw_args):
        obj=IndexedBase.__new__(cls, label,  **kw_args)
        obj.data=dict()
        if bases:
            obj.bases=bases
        return(obj)
    def __getitem__(self, indices, **kw_args):
        if not(is_sequence(indices)):
            # only one index 
            index=indices 
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch.")
            if type(index) is int:
                return(self.data[(index,)])
            else:    
                return VI(self, index, **kw_args)

        else:
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch.")
            if all(type(k) is int for k in indices):
                # all indices are integers
                return(self.data[indices])
            elif all( type(i) is Idx for i in indices):
                # all indices are symbolic
                expr=VI(self, *indices, **kw_args)
                cs=get_contraction_structure(expr)
                csk=cs.keys()
                if None in csk:
                    #no contraction
                    return(expr)
                else:
                    # compute the contracted indexed object or number
                    # since we are in __getitem__
                    #we assume that we deal here only with one VI instance not 
                    # with a many term expression
                    dummySuspects=csk[0]
                    for d in dummySuspects:
                        # there could be multiple contractions like x[_i,^i,^i,_i]
                        # which is the same as x[_i,^i,^j,_j] (or x[i,i,j,j] 
                        positions=[ p for p,ind in enumerate(indices) if ind==d ]
                        # find up down pairs
                        if (len(positions)%2!=0):
                            raise(ContractionException)
                        i1=indi


            elif all( type(i) is Idx or type(i) is int  for i in indices):
                # some indices are symbolic some are integers
                fixedIndexPositions=[i  for i in range(len(indices)) if not(type(indices[i]) is Idx) ]
                symbolicIndexPositions=[i  for i in range(len(indices)) if type(indices[i]) is Idx ]
                freeIndices=[i for i in indices if type(i) is Idx]
                #print(fixedIndexPositions)
                freeKeys=[k for k in self.data.keys() if all(k[p]==indices[p] for p in fixedIndexPositions )]
                newBase=VIB("intermediate")
                for k in freeKeys:
                    newBase.data[extractIndices(k,symbolicIndexPositions)]=self.data[k]
                return(VI(newBase,*freeIndices,**kw_args))
                
                
    
    
    
    def __setitem__(self, indices, value,**kw_args):
        l=len(self.data.keys())
        if l>0:
            #if we already have values set we make sure that the shape of the tensor is not changed
            if not(all([len(k)==len(indices)for k in self.data.keys()])):
                raise(IncompatibleShapeException())
        if not(isinstance(value,Indexed)):
            # on the rigth hand side is a "normal" python expression without indices
            if  not(is_sequence(indices)): #only one index
                index=indices
                if type(index) is int:
                    self.data[(index,)]=value
                else:
                    raise(IncompatibleShapeException("If v is a scalar (e.g. number,symbol,expression) \
                    then i in A[i]=v must be an integer and cannot be a symbolic index") )
            else:        
                if all(type(k) is int for k in indices):
                    # all indices are integers
                    self.data[indices]=value
                else:
                    raise(IncompatibleShapeException("If v is a scalar then i,j,..,k in A[i,j,...,k]=v must be integers and can not be symbolic.")) 
        
        else:
            # on the rigth hand side is an expression
            if  not(is_sequence(indices)): #only one index
                index=indices
                if type(index) is int:
                    self.data[indices]=value
                elif type(index) is Idx:
                    self.data=value.base.data
            else: #general case with more than one indices
                #if type(indices[0]) is int:
                if all(type(k) is int for k in indices):
                    # all indices are integers
                    self.data[indices]=value.indices
                
                elif all( type(i) is Idx for i in indices):
                    # all indices are symbolic
                    vb=value.base
                    vbd=vb.data
                    if set(indices)==set(value.indices):
                        if indices==value.indices:
                            self.data=vbd
                        else:
                            # a permutation of indices, necessiating a permutation of the components 
                            newPositions=[indices.index(i) for i in value.indices]
                            vKeys=vbd.keys()
                            for k in vKeys:
                                self.data[permuteTuple(k,newPositions)]=vbd[k]
                elif all( type(i) is Idx or type(i) is int  for i in indices):
                    # some indices are symbolic some are integers
                    newPositions=[indices.index(i) for i in value.indices]
                    vb=value.base
                    vbd=vb.data
                    vKeys=vbd.keys()

                    for k in vKeys:
                        self.data[changedTuple(indices,k,newPositions)]=vbd[k]
                    
         

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
        

