class Tensor2(object):
    ##########################################################
    def __init__(self,coords,baseName,componentTypes,components):
        # test if we can handle the component types 
        #if not(set(componentTypes).issubset({"roof","cellar","cart","phys"})):
        #    raise(Exceptions.UnknownComponentType(componentTypes))
        
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
                if len(componentTypes)!=1:
                    raise(Exceptions.ComponentTypesTupelMismatch(componentTypes,t0))
            # for higher order tensors we expect tuples as indeces    
            # we test if all indextupels have the same length and type
            elif isinstance(t0,tuple):
                l0=len(t0)
                wronglength=lambda y:len(y)!=l0
                if any(map(wronglength,indextupels)):
                    raise(Exceptions.IndexTupelError(indextupels))
                if len(componentTypes)!=l0:
                    raise(Exceptions.ComponentTypesTupelMismatch(componentTypes,t0))
            else:
               raise(Exceptions.IndexTupelTypeError(indextupels))
            # if the tests were successful copy the components
            for it in indextupels:
               self.components[it]=components[it]
        
        self.coords=coords
        self.componentTypes=componentTypes
        self.baseName=baseName
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

