#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import copy 
class S:
    def __init__(self,val):
        self.val=val
    
    def __radd__(self,other):
        cs= copy.deepcopy(self)
        if other==0:
            return(cs)
        else:    
            return(cs.val+other.val)

a=S(1)
b=S(2)
print(sum([a,b]))
