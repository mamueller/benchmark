#!/usr/bin/python
# vim: set expandtab ts=4
from sympy import *
class Coords:
    def __init__(self,val,scalarSimp):
        self.val=val
        self.scalarSimp=scalarSimp
        def matSimp(Mat):
            s=Mat.shape
            res=zeros(s)
            for i in range (s[0]):
                for j in range (s[1]): 
                    res[i,j]=self.scalarSimp(Mat[i,j])
            return(res)
        self.matSimp=matSimp
    
    def simplifyMats(self):
    #example for internal use of the manufactured method
        A=eye(3)
        print(self.matSimp(A))




        
# example coordinate frame
# we define a new coordinate frame simply as an istance of class coords
def testsimp(x):
    return(x*2)


sp=Coords(5,testsimp)
print(sp.matSimp(eye(2)))
sp.simplifyMats()

# additional special properties or even methods 
# of the coordinate frame are then 
# added later
sp.val2=5

def f(argthatWillBeBoundToSelf):
    print("this is f")
    print(argthatWillBeBoundToSelf.val)
import types
sp.f= types.MethodType( f, sp )
# call it
sp.f()
