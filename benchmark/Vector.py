#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import copy 
from Tensor import *
class Vector(Tensor):
    def __init__(self,coords,componentTypes,components):
        if len(componentTypes)>1:
            raise(ComponentsLengthError(components))
        sc=super(self.__class__,self)
        sc.__init__(coords,componentTypes,components)
        

