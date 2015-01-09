
from sympy.tensor.index_methods import get_contraction_structure, get_indices, _remove_repeated
from sympy.diffgeom import Manifold, Patch, CoordSystem
from sympy.tensor import Idx
from sympy.printing import pprint
from sympy import symbols, default_sort_key, Matrix
from sympy.tensor import IndexedBase, Indexed, Idx
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException,IndexRangeException,MissingBaseException, MissingConnectionException, NotImplementedException
from TupelHelpers import  permuteTuple, extractIndices,changedTupleAtPos,changedTuple, deleteIndices, tupleLen
from copy import deepcopy
class PatchWithMetric(Patch):
    def __init__(self, name, manifold):
        super(PatchWithMetric, self).__init__(name, manifold)
        # a patch has only one metric but the
        # representation of the metric tensor depends on
        # the coordinate system.
        # Given one representation in one coordinate system A
        # onother representation for coordinate system B
        # can be inferred from the connection of A and B 
        # since the computation can be expensive we cache 
        # the results 
        self.metricRepresentations=dict()
        # same here
        self.ChristoffelCache=dict()

    def setMetric(self,CordSystem,vi):
        # vi is a representation of the metric tensor (a TensorIndexSet)
        # in either roof,roof or cellar,cellar components
        # since the mixed components would be the cronecker delta
        # which lacks the essential metric information.
        bases=vi.bases
        if all([isinstance(b,VectorFieldBase) for b in bases]) or all([isinstance(b,OneFormFieldBase) for b in bases]):
            # create a dictionary indexed by CoordSystems
            self.metricRepresentations[CordSystem]=vi
        else:    
            raise "the metric tensor has to be given in roof roof, or cellar cellar components since the metric information can not be inferred from the cronnecker delta"


    def  christoffelFromMetric(self,metric):
        # compute the christophelsymbols (of the second kind) 
        # we need 2 component sets of the metric tensor
        # roof,roof and cellar,cellar
        i=Idx("i")
        j=Idx("j")
        k=Idx("k")
        l=Idx("p")
        bases=metric.bases
        coords=bases[0].coords
        reciprocalBase=OneFormFieldBase(coords)
        # 1. case metric is given as roof roof
        grr=metric 
        if all([isinstance(b,VectorFieldBase) for b in bases]): 
            gcc=TensorIndexSet("gcc",[reciprocalBase,reciprocalBase])
            mat=Matrix(self.dim,self.dim)
            for i in range(dim):
                for j in range(dim):
                    mat[i,j]=grr[i,j]
            InvMat=mat.inverse_LU()
            for i in range(dim):
                for j in range(dim):
                    gcc[i,j]= InvMat[i,j]#exercise. 2.10 Simmonds
        
        cs=dict() 
        dim=self.dim
        for k in range(dim):
            for i in range(dim):
                for j in range(dim):
                    cs[k,i,j]=sum(
                    [1/2*grr[k,p]
                    *(
                     diff(gcc[i,p],j)
                    +diff(gcc[j,p],i)
                    -diff(gcc[i,j],p)
                    ) for p in range(dim)])
        return(cs)            
                    
    def deriveMetricFromConnection(self,nameOfTargetCoordinateSystem):  
        # check if the patch stores a coordsystem with the right name
        matches=[cs for cs in self.coord_systems if cs.name==nameOfTargetCoordinateSystem]
        if len(matches)<1:
            raise "no coordinate system of the specified name found"
        elif len(matches)>1:
            raise "more than one coord system with the specified name found"
        else:   
            targetCoordinateSystem=matches[0]
        connections=targetCoordinateSystem.transforms
        #connections is a dictionary whose keys are CoordSystems
        # check if there is at least one metric representation in a connected coord system
        connectionsWithMetricRepresentation={key:val for key,val  in connections.items() if key in self.metricRepresentations }
        # in most cases there will be only one match (e.g. the cartesian cs) but it is possible to have more
        # so we chose the first that comes our way
        
        firstCoordSystem=list(connectionsWithMetricRepresentation.keys())[0]
        # now the question is to which of the metric representations of the patch it refers
        
        symbols,equations=connectionsWithMetricRepresentation[firstCoordSystem]
        #compute the jacobian
        J=targetCoordinateSystem.jacobian(firstCoordSystem,symbols)
        pprint(J)
        # the first column of the jacobian contains the components of the 
        # new base vectors in the old base
        # therefor the components of the new metric tensor
        # can be derived from the components of the old
        # by J.transposeLU()*g_old*J          
        g_old=self.metricRepresentations[firstCoordSystem]
        oldBases=g_old.bases
        oB0=oldBases[0] # first base
        if  isinstance(oB0,VectorFieldBase):
            # the second base must be of the same type 
            # in this case also
            newBase=VectorFieldBase("newBase",targetCoordinateSystem)
            newBase.connect(oB0,J)
            print("blub2")
            print(newBase)
        i, j= map(Idx, ['i', 'j'])
        g_new[i,j]=g_old[i,j].rebase2([newBase,newBase])
        




#        raise "continue here to derive the components of the metric tensor in the new base"
    def getMetricRepresentation(self,nameOfCoordinateSystem):   
        if nameOfCoordinateSystem in self.metricRepresentations.keys():
            metric=self.metricRepresentations[nameOfCoordinateSystem]
        else:
            s=self.deriveMetricFromConnection(nameOfCoordinateSystem)
            # update chache
            self.metricRepresentations[nameOfCoordinateSystem]=s
        return(s)
    def christoffelSymbols(self,nameOfCoordinateSystem):
        print("blub1")
        print(nameOfCoordinateSystem)
        if nameOfCoordinateSystem in self.ChristoffelCache.keys():
            s=self.ChristoffelCache[nameOfCoordinateSystem]
        else:
            metric=self.getMetricRepresentation(nameOfCoordinateSystem)
            print("blub3")
            print(metric)
            s=self.christoffelFromMetric(metric)
            # update chache
            self.ChristoffelCache[nameOfCoordinateSystem]=s
        return(s)

##########################################################
class BaseVector(object):
    # represents a cellar base vector
    def __init__(self,coordSystem,ind):
        self.ind=ind
        self.coord=coordSystem
    def _eval_derivative(self):
        # find christoffel symbols
        coord=self.coord
        cname=coord.name
        patch=coord.patch
        cs=patch.christoffelSymbols(cname)
        return(5)
        
   
   
##########################################################
class VectorFieldBase(object):
    """This class stores relationships between bases on the same patch.
    Every base b_i has one and only one dual (or reciprocal) base b^j
    with b^j(b_i)= <b^j,b_i>= delta_i^j but can be connected to other bases
    by linear mappings.
    """
    def __init__(self,name,coordSystem=None):
        self.name=name
        self.baseConnections={}
        if coordSystem:
            self.coord=coordSystem
            self.vecs=[BaseVector(coordSystem,i) for i in range(coordSystem.dim)]


    def getDual(self):
        if 'dual' in self.__dict__:
            return(self.dual)
        else:
            return(OneFormFieldBase(self))

    def connect(self,oldBase,Mat):
        self.baseConnections[oldBase]=Mat

##########################################################
class OneFormFieldBase(object):
    #def __init__(self,vfb):
    #    self.dual=vfb
    #    # register in the base
    #    vfb.dual=self
           
    def __new__(cls,VectorFieldBaseInstance=None):
        # a OneFormFieldBase is by definition the dual base
        # to a VectorFieldBaseInstance
        # So we check if such a dual base has already been created
        # and registered in the VectorFieldBaseInstance
        # If so we copy the object instead of creating a new one
        # this makes sure that there is only one dual base
        if VectorFieldBaseInstance:
            if "dual" in VectorFieldBaseInstance.__dict__:
                obj=VectorFieldBaseInstance.dual
            else:
                obj=object.__new__(cls)
                obj.dual=VectorFieldBaseInstance
                # register in VectorFieldBaseInstance
                VectorFieldBaseInstance.dual=obj
        else:
            obj=object.__new__(cls)
           
        obj.baseConnections={}
        return(obj)
    def getDual(self):
        if 'dual' in self.__dict__:
            return(self.dual)
        else:
            b=VectorFieldBase("")
            b.dual=self
            return(b)

    def connect(self,oldBase,Mat):
        self.baseConnections[oldBase]=Mat

##########################################################
class TensorIndexSet(IndexedBase):
    # this thing represents an indexset of a  tensor w.r.t a fixed base
    # It is not the tensor itself which is a multilinear mapping expressable w.r.t to any base 
    # but one such expression
    def __deepcopy__(self,memo):
        cls = self.__class__
        print(cls)
        result = cls.__new__(cls,"intermediate")
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        # for the bases we do not want copies but references    
        result.bases=self.bases
        return result
    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result


    def __new__(cls,label,bases=None, **kw_args):
        obj=IndexedBase.__new__(cls, label,  **kw_args)
        obj.data=dict()
        if bases:
            obj.bases=bases
        else:
            obj.bases=None
            #raise(MissingBaseException)
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

                newTensorIndexSet=TensorIndexSet("intermediate",remainingBases)
                newTensorIndexSet.data={}
                for k in freeKeys:
                    newTensorIndexSet.data[extractIndices(k,symbolicIndexPositions)]=self.data[k]
                # now check the remaining part for contraction
                expr=VI(newTensorIndexSet,*freeIndices,**kw_args)
                cs=get_contraction_structure(expr)
                csk=list(cs.keys())
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
                        b0=newTensorIndexSet.bases[p0].dual    
                        b1=newTensorIndexSet.bases[p1].dual    
                        if b1!=b0.dual:
                            # implement only natural pairing (up and down indices)
                            raise(ContractionIncompatibleBaseException(d,p0,p1,b0,b1))
                        newBases=[b for i,b in enumerate(newTensorIndexSet.bases) if not(i in dummyPositions)] 
                        nonDummyPositions=[ p for p,ind in enumerate(freeIndices) if ind!=d ]
                        if len(newBases)==0:
                            # the result will be a scalar
                            res=sum([newTensorIndexSet.data[k] for k in newTensorIndexSet.data.keys() if k[p0]==k[p1]])
                            return(res)   
                        else:
                            newData={}
                            newKeys=set([deleteIndices(k,dummyPositions) for k in newTensorIndexSet.data.keys()])
                            for nk in newKeys:
                                newData[nk]=sum([newTensorIndexSet.data[k] for k in newTensorIndexSet.data.keys() if k[p0]==k[p1] and deleteIndices(k,dummyPositions)==nk])
                            newTensorIndexSet.data=newData
                            newTensorIndexSet.bases=newBases
                            freeIndices=[i for p,i in enumerate(freeIndices) if p in nonDummyPositions]
                            expr=VI(newTensorIndexSet, *freeIndices, **kw_args)
                        
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
                if len(value.indices)==1:
                    vIndex=value.indices[0]
                else:
                    raise(IncompatibleShapeException())

                if type(index) is int:
                    self.data[indices]=value
                elif type(index) is Idx:
                    if index.upper==vIndex.upper and index.lower==vIndex.lower:
                        self.data=value.base.data
                        vb=value.base
                        self.bases=vb.bases
                    else:
                        raise(IndexRangeException)
            else: #general case with more than one index
                symbolicIndices=[i for i in indices if type(i) is Idx]
                #integerIndices=[i for i in indices if type(i) is int]
                lsi=len(symbolicIndices)
                if len(value.indices)!=lsi:
                    raise(IncompatibleShapeException())
                    # note that the value indices are preprocessed by []
                else:
                    if any([symbolicIndices[i].lower!=value.indices[i].lower or symbolicIndices[i].upper!=value.indices[i].upper for i in range(lsi)]):
                        raise(IndexRangeException)
                # The lhs expressions can contain indices that are also 
                # present on the rhs, those are matched first
                inter=set(symbolicIndices).intersection(set(value.indices))
                li=len(inter)
                # first handle positions marked by common indices
                newPositions=[indices.index(i) for i in inter]
                oldPositions=[value.indices.index(i) for i in inter]
                # now determine free positions on the lhs
                freeLhs=set(symbolicIndices).difference(inter)
                freeRhs=set(value.indices).difference(inter)
                lf=len(freeLhs)
                freeLhsPositions=[indices.index(i) for i in freeLhs]
                freeRhsPositions=[value.indices.index(i) for i in freeRhs]
                vb=value.base
                self.bases=vb.bases
                vbd=vb.data
                vKeys=vbd.keys()
                newKey=indices 
                # the new key can still contain symbolic indices
                # first we handle those which are present in the rhs
                for k in vKeys:
                    for i in range(li):
                        newKey=changedTupleAtPos(newKey,k[oldPositions[i]],newPositions[i])
                    # now the remaining free symbolic idices are overwritte
                    # by the values of the rhs key
                    for i in range(lf):
                        newKey=changedTupleAtPos(newKey,k[freeRhsPositions[i]],freeLhsPositions[i])
                    # now the remaining free symbolic idices are overwritte
                    self.data[newKey]=vbd[k]
                    
##########################################################
class OIB(TensorIndexSet):
    pass
##########################################################

class VI(Indexed):
    def __new__(cls, base, *args, **kw_args):
        from sympy.utilities.misc import filldedent
        if not args:
            raise IndexException("Indexed needs at least one index.")
        if isinstance(base, (string_types, Symbol)):
            base = TensorIndexSet(base)
        elif not isinstance(base, TensorIndexSet):
            raise TypeError(filldedent("""
                Indexed expects string, Symbol or TensorIndexSet as base."""))

        args = list(map(sympify, args))
        return Expr.__new__(cls, base, *args, **kw_args)



    @staticmethod
    def getBaseTransformationTensor(nb,ob):
        dob=ob.getDual()
        I=TensorIndexSet("I",[dob,nb])
        # Note that a connection of bases is a transitiv relation
        # e.g.  O N E  connection of 
        # C1.) new cellar   to old cellar base  ,      which immediately implies  its inversese 
        # C2.) old cellar   to new cellar base  ,      and also by the definition of the dual bases
        # C3.) new roof     to old roof         ,      and its inverse 
        # C4.) old roof     to new roof         
        
        # so up to now it's bye one, get four

        # if you now add a metric you add more possibilities
        # because you can now raise and lower indices which means
        # that you can connect a base to its dual which gives you
        # C5.) new roof     to new cellar
        # C6.) new cellar   to new roof
        # C7.) old roof     to old cellar
        # C8.) old cellar   to old roof
        # 
        # which you can combine with C1,..,C4
        # to get some new possibilities
        # 
        # C9.)  new roof   to old cellar         
        # C10.) old cellar to new roof 
        # C11.) new cellar to old roof 
        # C12.) old roof   to new cellar
        
        # This was a first example of trasitivity 
        # we can extend it by a base 
        # A1 connected to A2 and B1
        #   which is connected to B2 and C1
        #        which is connected to ....
        
        # This implies a lot more connections 
        # e.g. from A1 to C2.... which we could compute by composition (or matrix multiplication)
        

        # This raises questions how to implement this generaly
        #
        # For now we restrict ourselfes willfully to the already tested cases
        # these are 
        

        if nb.__class__==ob.__class__: # no raising and lowering involved
            dob=ob.getDual()
            dnb=nb.getDual()
            # try if the directly usable connection is already stored
            if dob in dnb.baseConnections.keys():
                TensMat=dnb.baseConnections[dob]
            elif ob in nb.baseConnections.keys():
                TensMat=(nb.baseConnections[ob]).inverse_LU()
            else:
                # we could try to look further and try to infer 
                # the connection but 
                raise(NotImplementedException)
        else:
            # raising and lowering would be involved 
            # we do not do this at the moment
            raise(NotImplementedException)
        # assemble the Transformation Tensor from the matrix        
        dim=TensMat.shape[0]
        for i in range(dim):
            for j in range(dim):
                I[i,j]= TensMat[i,j]
        return(I)


    def rebase2(self,newBases):
        indices=self.indices
        oldBases=self.base.bases
        lo=len(oldBases)
        ln=len(newBases)
        if lo!=ln:
            raise(IncompatibleShapeException())
        newTensorIndexSet=deepcopy(self.base)
        for pos in range(ln):
            nb=newBases[pos]
            ob=oldBases[pos]
            I=self.getBaseTransformationTensor(nb,ob)
            dummy=Idx("dummy")
            orgind=indices[pos]
            newkey=changedTupleAtPos(indices,dummy,pos)
            newTensorIndexSet[indices]=newTensorIndexSet[newkey]*I[dummy,orgind]

        return(newTensorIndexSet[indices])
    
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
        
        if type(sb) == TensorIndexSet:
            # we create an intermediate outer product and 
            # invisibly contract it with the use of the [] operator 
            # and the original indices if indices are repeated 
    
            intermediate=TensorIndexSet("intermediate",allbases)
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
        # extract the first vector field base
        b0=bases[0]
        print(type(b0))
        # extract the first base vector 
        er=b0.vecs[0]
        print(type(er))
        return(er._eval_derivative())
            
            
        

