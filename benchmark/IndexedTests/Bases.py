#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
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
    def connect(self,oldBase,Mat):
        self.baseConnections[oldBase]=Mat

##########################################################
class OneFormFieldBase(object):
    #def __init__(self,vfb):
    #    self.dual=vfb
    #    # register in the base
    #    vfb.dual=self
    def __new__(cls,vfb=None,**kwargs):
        # check if the there already is a dual base to vfb
        # if so copy the object instead of creating a new one
        # this makes sure that there is only one dual base
        if vfb:
            if "dual" in vfb.__dict__:
                obj=vfb.dual
            else:
                obj=object.__new__(cls)
                obj.dual=vfb
                # register in vfb
                vfb.dual=obj
        else:
            obj=object.__new__(cls)
        
        if "label" in kwargs:
            obj.label=kwargs["label"]
            
           
        return(obj)

bc=VectorFieldBase("bc")
# make sure that the dual base stays the same object
br=OneFormFieldBase(bc)
print(bc.dual)
br2=OneFormFieldBase(bc)
print(bc.dual)
br3=OneFormFieldBase()
br4=OneFormFieldBase(label="br")
print(br4.label)
