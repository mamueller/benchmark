
from sympy.tensor.index_methods import get_contraction_structure, get_indices, _remove_repeated
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
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
                freeKeys=[k for k in self.data.keys() if all(k[p]==indices[p] for p in fixedIndexPositions )]
                if self.bases:
                    remainingBases=[b for p,b in enumerate(self.bases) if p in symbolicIndexPositions]
                else:
                    remainingBases=None

                newVIB=VIB("intermediate",remainingBases)
                newVIB.data={}
                for k in freeKeys:
                    newVIB.data[extractIndices(k,symbolicIndexPositions)]=self.data[k]
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
           lk=tupleLen(list(self.data.keys())[0])
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
        if not(isinstance(value,VI)):
            # on the rigth hand side is a "normal" python expression without indices
            if  not(is_sequence(indices)): #only one index
                index=indices
                if type(index) is int:
                    self.data[(index,)]=value
                else:
                    print(value)
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
                if len(value.indices)!=1:
                    raise(IncompatibleShapeException())
                if type(index) is int:
                    self.data[indices]=value
                elif type(index) is Idx:
                    self.data=value.base.data
            else: #general case with more than one indices
                if len(value.indices)!=len([i for i in indices if type(i) is Idx]):
                    # the value indices are preprocessed by []
                    raise(IncompatibleShapeException())
                # some indices are symbolic some are integers
                newPositions=[indices.index(i) for i in value.indices]
                vb=value.base
                vbd=vb.data
                vKeys=vbd.keys()

                for k in vKeys:
                    self.data[changedTuple(indices,k,newPositions)]=vbd[k]
                    
##########################################################
class OIB(VIB):
    pass
##########################################################

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
       
        sb=self.base
        sbd=sb.data
        si=self.indices
    
        oi=other.indices
        ob=other.base
        obd=ob.data
        if not(sb.bases and ob.bases):
            raise "one argument has no base defined"
        allbases=sb.bases+ob.bases
        
        if type(sb) == VIB:
            # we create an intermediate outer product and 
            # invisibly contract it with the use of the [] operator 
            # and the original indices if indices are repeated 
    
            intermediate=VIB("intermediate",allbases)
            for ks in sb.data.keys():
                for ko in ob.data.keys():
                    intermediate[ks+ko]=sbd[ks]*obd[ko]
            return(intermediate[si+oi])
        elif type(sb)==OIB:
            # this applies to the nabla operator 
            
            # we apply each component of the operator first to the other tensor 
            # and multiply the result to a tensor (not operator) with the same shape
            # as 


            parts={} 
            for ks in sb.data.keys():
                    ib=VI("substituteMultiplikator",sb.bases)
                    ib.data[ks]=1.0,
                    parts[ks]=ib[si]*sbd[ks](other)# this takes care of possible contractions 
                    # e.g in the case of divergene
            res=sum([parts[ks] for ks in parts.keys()  ])
    @property
    def free_symbols(self):
        data=self.base.data
        s=set({})
        for v in data.values():
            s=s.union(v.free_symbols)
        return(s)


    def _eval_derivative(self, x):
        # we treat an indexed object as a sum of dyads
        # The derivative is then the sum of the derivatives of the dyads
        # The derivative of a dyad is computed by the product rule from 
        # the derivative of product of the components and the base dyads
        # the derivative of a basedyad is computed by the product rule
        # from the derivatives of the base vectors
        
        #for the moment just face the dyads and components and go directly to the
        # vector derivative
        IB=self.base
        # bases of the indexed object
        bases=IB.bases
        # extract the first vector fieldbase
        b0=bases[0]
        print(type(b0))
        # extract the first base vector 
        er=b0.vecs[0]
        print(type(er))
        return(er._eval_derivative())
        

