
from sympy.tensor.index_methods import get_contraction_structure, get_indices, _remove_repeated
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption
from TupelHelpers import  permuteTuple, extractIndices, changedTuple, deleteIndices, tupleLen
from Bases import OneFormFieldBase, VectorFieldBase
from copy import deepcopy
##########################################################
class VIB(IndexedBase):
    def __new__(cls,label,bases=None, **kw_args):
        obj=IndexedBase.__new__(cls, label,  **kw_args)
        obj.data=dict()
        if bases:
            obj.bases=bases
        else:
            obj.bases=None
        return(obj)

    def __getitem__(self, indices, **kw_args):
        self.checkShapeCompatibility(indices)
        
        if not(is_sequence(indices)):
            # only one index 
            index=indices 
            if self.shape and len(self.shape) != 1:
                raise IncompatibleShapeException()
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

            elif all( type(i) is Idx or type(i) is int  for i in indices):
                # some indices are symbolic some are integers
                # first handle the integer indices
                fixedIndexPositions=[i  for i in range(len(indices)) if not(type(indices[i]) is Idx) ]
                symbolicIndexPositions=[i  for i in range(len(indices)) if type(indices[i]) is Idx ]
                freeIndices=[i for i in indices if type(i) is Idx]
                #print(fixedIndexPositions)
                freeKeys=[k for k in self.data.keys() if all(k[p]==indices[p] for p in fixedIndexPositions )]
                if self.bases:
                    remainingBases=[b for p,b in enumerate(self.bases) if p in symbolicIndexPositions]
                else:
                    remainingBases=None

                newVIB=VIB("intermediate",remainingBases)
                newVIB.data={}
                for k in freeKeys:
                    newVIB.data[extractIndices(k,symbolicIndexPositions)]=self.data[k]
                print(newVIB.data)
                # now check the remaining part for contraction
                expr=VI(newVIB,*freeIndices,**kw_args)
                cs=get_contraction_structure(expr)
                csk=cs.keys()
                if None in csk:
                    #no contraction
                    return(expr)
                else:
                    # compute the contracted indexed object or number.
                    # Since we are in __getitem__  we assume that we deal here 
                    # only with one VI instance not with a many term expression
                    dummySuspects=csk[0]
                    #print(dummySuspects)
                    #freeIndices=indices
                    for d in dummySuspects:
                        cs=get_contraction_structure(expr)
                        csk=cs.keys()
                        # there could be multiple contractions like x[_i,^i,^i,_i]
                        # which is the same as x[_i,^i,^j,_j] (or x[_i,^j,^i,_j]
                        # but we do not implement this because it is not very clear  
                        dummyPositions=[ p for p,ind in enumerate(freeIndices) if ind==d ]
                        # find pairs
                        n=len(dummyPositions)
                        if n!=2:
                            # exclude things like [i,i,i,j] where one index occures more than 2 times
                            # rendering the contraction ambigous
                            raise(ContractionException(d,n))
                        p0=dummyPositions[0]     
                        p1=dummyPositions[1]     
                        b0=newVIB.bases[p0].dual    
                        b1=newVIB.bases[p1].dual    
                        if b1!=b0.dual:
                            # implement only natural pairing (up and down indices)
                            raise(ContractionIncompatibleBaseException(i,p0,p1,b0,b1))
                        newBases=[b for i,b in enumerate(newVIB.bases) if not(i in dummyPositions)] 
                        nonDummyPositions=[ p for p,ind in enumerate(freeIndices) if ind!=d ]
                        if len(newBases)==0:
                            # the result will be a scalar
                            res=sum([newVIB.data[k] for k in newVIB.data.keys() if k[p0]==k[p1]])
                            return(res)   
                        else:
                            newData={}
                            newKeys=set([deleteIndices(k,dummyPositions) for k in newVIB.data.keys()])
                            for nk in newKeys:
                                newData[nk]=sum([newVIB.data[k] for k in newVIB.data.keys() if k[p0]==k[p1] and deleteIndices(k,dummyPositions)==nk])
                            newVIB.data=newData
                            newVIB.bases=newBases
                            freeIndices=[i for p,i in enumerate(freeIndices) if p in nonDummyPositions]
                            expr=VI(newVIB, *freeIndices, **kw_args)
                        
                    return(expr)
                
                
    
    
    def checkShapeCompatibility(self,indices):
        if self.bases:
            lb=len(self.bases)
        else:
            lb=None

        ld=len(self.data)
        if ld>0:
           lk=tupleLen(self.data.keys()[0])
        else:
           lk=None

        if lk:
            #if we already have set values we make sure that the shape of the 
            # tensor is not changed by the new assignment
            if not(is_sequence(indices)):
                if lk!=1:
                    raise(IncompatibleShapeException())
            else:
                li=len(indices)
                if not(all([len(k)==li for k in self.data.keys()])):
                    raise(IncompatibleShapeException())
        if lb:
            #if we already have set  bases we make sure that the shape of the 
            # tensor is not changed by the new assignment
            if not(is_sequence(indices)):
                if lb!=1:
                    raise(IncompatibleShapeException())
            else:
                li=len(indices)
                if lb!=li: 
                    raise(IncompatibleShapeException())
                   

    def __setitem__(self, indices, value,**kw_args):
        self.checkShapeCompatibility(indices)
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
        #print("in mul")
        #print(self)
        ## we now create an Indexed object and can then use all the methods 
        expr=super(VI,self).__mul__(other)
        #print(type(expr)) 
        #print(expr) 
        #print(get_contraction_structure(expr))
        #print(get_indices(expr))
        # now we can uses this information to implement the the multiplication for our subtype.
        

