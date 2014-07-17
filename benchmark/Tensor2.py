#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import copy 
from helperFunctions import pp
import Exceptions 
from sympy.matrices import zeros
from sympy.tensor import IndexedBase, Idx, Indexed

def components2vec(components):
    vec=zeros(3,1)
    for i in range(0,3):
        ind=(i,)
        if (ind in components.keys()): 
            vec[i]=components[ind]
    return(vec)
###########################################################
def vec2components(vec):
    comp={}
    for i in range(0,len(vec)):
        comp[(i,)]=vec[i]
    return(comp)

###########################################################
class Tensor2(Indexed):
    
    def __init__(self,coords,baseNames,components):
        
        self.components={} # initialize empty (Null Tensor)

        indextupels=components.keys()
        
        # if any indextuples are given then perform some compatibility tests
        # otherwise proceed with the initialization of the Null tensor
        if len(indextupels)!=0:
            t0=list(indextupels)[0]
            # for first order Tensors we allow itegers as indeces instead of tuples
            if isinstance(t0,int): 
                wrongtype=lambda y:not(isinstance(y,int))
                # make sure that all other indeces are integers to
                if any(map(wrongtype,indextupels)):
                    raise(Exceptions.IndexTupelError(indextupels))
                # make sure that it really is a first order Tensor
                if len(baseNames)!=1:
                    raise(Exceptions.ComponentTypesTupelMismatch(baseNames,t0))
            # for higher order tensors we expect tuples as indeces    
            # we test if all indextupels have the same length and type
            elif isinstance(t0,tuple):
                l0=len(t0)
                wronglength=lambda y:len(y)!=l0
                if any(map(wronglength,indextupels)):
                    raise(Exceptions.IndexTupelError(indextupels))
                if len(baseNames)!=l0:
                    raise(Exceptions.ComponentTypesTupelMismatch(baseNames,t0))
            else:
               raise(Exceptions.IndexTupelTypeError(indextupels))
            # if the tests were successful copy the components
            for it in indextupels:
               self.components[it]=components[it]
        
        self.coords=coords
        self.baseNames=baseNames
        r=range(0,coords.n)
        self.r=r
        self.purge()
    ##########################################################
    def purge(self):
        # If a component is zero it can be safely erased
        # this is important to  make vectors and tensors
        # with zero components compareable
        c=self.components
        for k in list(c.keys()):
            if c[k]==0:
               del(self.components[k])

    ##########################################################
    def __eq__(self,other):
        self.purge()
        other.purge()
        boolval=\
        type(self.coords)==type(other.coords) and \
        self.baseNames==other.baseNames and\
        self.components==other.components
        return(boolval)
    
###########################################################
class ChangeOfBase(object):
    '''This class represents a change of base vectors'''
    '''It is needed to compute the components of vectors and Tensors with respect to a different base '''
    def __init__(self,source,dest,mat):
        # The matrix columns describe the source base  in terms of the 
        #  target basis 
        self.mat=mat
        self.dest=dest
        self.source=source
        nr,nc=mat.shape
        if nr !=nc:
            raise(Exception,"a coordinate Transformation always has to be quadratic")
        self.mat=mat
        self.n=nr
    ###########################################################
    def transform(self,tens,pos):
    # transform a tensor to another base in the pos th component
        T=copy.deepcopy(tens)
        c=T.coords
        cT=T.baseNames
        if cT[pos]!=self.source:
            raise(Exceptions.TransformError(self.source,self.dest,cT,pos))
        lcT=len(cT)
        pp("lcT",locals())
        Tc=T.components
        mat=self.mat
        tupelLength=len(list(Tc.keys())[0])
        if tupelLength==1:
            vec=components2vec(Tc)
            resvec=mat*vec
            rescomponents=vec2components(resvec)
            T.components=rescomponents
            T.baseNames=[self.dest]
            return(T)
            
        elif tupelLength==2:
            Tck=Tc.keys()
            if pos==0: 
                right_keys=extract_keys(Tck,1,1)
                pp("right_keys",locals())
                newComponents={}
                for k in right_keys:
                    v=extractVectorComponents(Tc,k+("*",))
                    vec=components2vec(v)

                    for l in range(0,self.n):
                        val=(mat.row(l)).dot(vec)
                        if val !=0:
                            newComponents[(k+(l,))]=val
            elif pos==1:
                left_keys=extract_keys(Tck,0,0)
                pp("left_keys",locals())
                newComponents={}
                for k in left_keys:
                    v=extractVectorComponents(Tc,k+("*",))
                    vec=components2vec(v)

                    for l in range(0,self.n):
                        val=(mat.row(l)).dot(vec)
                        if val !=0:
                            newComponents[(k+(l,))]=val
            T.components=newComponents
            T.baseNames[pos]=self.dest
            return(T)        
    
            
#############################################################
    def invTransform(self,TensorComponents,pos):
        #pp("TensorComponents",locals())
        invMat=self.mat.inv()
        tupelLength=len(TensorComponents.keys()[0])
        if tupelLength==1:
            vec=components2vec(TensorComponents)
            resvec=invMat*vec
            rescomponents=vec2components(resvec)
            raise("return a tensor")
            return(rescomponents)
        elif tupelLength==2:
            cck=TensorComponents.keys()
            if pos==0: 
                right_keys=extract_keys(cck,1,1)
                pp("right_keys",locals())
                newComponents={}
                for k in right_keys:
                    v=extractVectorComponents(TensorComponents,k+("*",))
                    vec=components2vec(v)

                    for l in range(0,self.n):
                        val=(invMat.row(l)).dot(vec)
                        if val !=0:
                            newComponents[(k+(l,))]=val
                raise("return a tensor")
                return(newComponents)        
            elif pos==1:
                left_keys=extract_keys(cck,0,0)
                pp("left_keys",locals())
                newComponents={}
                for k in left_keys:
                    v=extractVectorComponents(TensorComponents,k+("*",))
                    vec=components2vec(v)

                    for l in range(0,self.n):
                        val=(invMat.row(l)).dot(vec)
                        if val !=0:
                            newComponents[(k+(l,))]=val
                raise("return a tensor")
                return(newComponents)        
                #pp("left_keys",locals())
