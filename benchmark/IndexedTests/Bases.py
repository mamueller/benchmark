#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
from Exceptions import IncompatibleShapeException, DualBaseExeption, BaseMisMatchExeption,ContractionIncompatibleBaseException
from sympy.diffgeom import Manifold, Patch, CoordSystem
from sympy.tensor import Idx

class PatchWithMetric(Patch):
    def __init__(self, name, manifold):
        super(PatchWithMetric, self).__init__(name, manifold)
        self.metricsCache=dict()
        self.ChristoffelCache=dict()

    def setMetric(self,nameOfCordSystem,mat):
        # mat is considered to be g_ij (cellar components)
        self.metricsCache[nameOfCordSystem]=mat
    def  christoffelFromMetric(self,metric):
        # compute the christophelsymbols (of the second kind) 
        # we need 2 component sets of the metric tensor
        # roof,roof and cellar,cellar
        i=Idx("i")
        j=Idx("j")
        k=Idx("k")
        l=Idx("p")
        IB=metric.base
        bases=IB.bases
        coords=bases[0].coords
        reciprocalBase=OneFormFieldBase(coords)
        # metric is given as roof roof
        #grr=VIB("grr",bases) #maybe we can get rid of the bases here and write 
        grr=VIB("grr")
        grr[k,p]=metric 
        if all([isinstance(b,VectorFieldBase) for b in bases]): 
            gcc=VIB("gcc",[reciprocalBase,reciprocalBase])
            mat=Matrix(self.dim,self.dim)
            for i in range(dim):
                for j in range(dim):
                    mat[i,j]=grr[i,j]
            InvMat=mat.inverseLU()
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
                    
    def deriveMetricFromConnection(self,nameOfCoordinateSystem):  

        coord=[cs for cs in self.coord_systems if cs.name==nameOfCoordinateSystem][0]
        print(coord.transforms)
        raise "continue here to derive the components of the metric tensor in the new base")
    def metrics(self,nameOfCoordinateSystem):   
        if nameOfCoordinateSystem in self.metricsCache.keys():
            metric=self.metricsCache[nameOfCoordinateSystem]
        else:
            s=self.deriveMetricFromConnection(nameOfCoordinateSystem)
            # update chache
            self.metricsCache[nameOfCoordinateSystem]=s
        return(s)
    def christoffelSymbols(self,nameOfCoordinateSystem):
        if nameOfCoordinateSystem in self.ChristoffelCache.keys():
            s=self.ChristoffelCache[nameOfCoordinateSystem]
        else:
            metric=self.metrics(nameOfCoordinateSystem)
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
    with b^j(b_i)= <b^j,b_i>= delta_i^j . 
    """
    def __init__(self,coordSystem):
        self.coord=coordSystem
        self.vecs=[BaseVector(coordSystem,i) for i in range(coordSystem.dim)]

##########################################################
class OneFormFieldBase(VectorFieldBase):
    def __init__(self,coordSystem):
        self.oneforms=coordSystem.base_oneforms()
