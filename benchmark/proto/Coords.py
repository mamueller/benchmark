#!/usr/bin/python
# vim: set expandtab ts=4
from sympy import *
import unittest 
class CoordTest(unittest.TestCase):
   # def test_frames_as_instances(s):
   # ## we define a concrete coordinate frame simply as an istance of class coords
   #     class Coords(object):
   #         def __init__(self,X,scalarSimp):
   #             self.X=X
   #             self.scalarSimp=scalarSimp
   #             def matSimp(Mat):
   #                 s=Mat.shape
   #                 res=zeros(s)
   #                 for i in range (s[0]):
   #                     for j in range (s[1]): 
   #                         res[i,j]=self.scalarSimp(Mat[i,j])
   #                 return(res)
   #             self.matSimp=matSimp

   #         # example for internal use of the method
   #         def simplifyMats(self):
   #             return(self.matSimp(self.X*eye(3)))

   #     
   #         
   #  ################### begin checks###############################            
   #     x,y,z=symbols("x y z")
   #     def testsimp(exp):
   #         return(exp.subs({x:y}))

   #     sp=Coords(x,testsimp)
   #     s.assertEqual(y*eye(2),sp.matSimp(x*eye(2)))
   #     s.assertEqual(y*eye(3),sp.simplifyMats())
   #     
   #     # additional special properties or even methods 
   #     # of the coordinate frame are then 
   #     # added later
   #     sp.val2=5
   #     
   #     def f(argthatWillBeBoundToSelf):
   #         print("this is f")
   #         print(argthatWillBeBoundToSelf.X)
   #     import types
   #     sp.f=types.MethodType( f, sp )
   #     # call it
   #     sp.f()
   #     
   #     
   #     
   #     
        
        




################################## subclassing ###################################
    # The idea to use factories for different coordinate frames can be seen as 
    # an intermediate way to a subclass (a class is in fact a factory for instances)
    def test_coordframes_as_subclasses(s):
        class Coords(object):
            def scalarSimp(self,exp):
                pass

            def matSimp(self,Mat):
                s=Mat.shape
                res=zeros(s)
                for i in range (s[0]):
                    for j in range (s[1]): 
                        res[i,j]=self.scalarSimp(Mat[i,j])
                return(res)

        class Spherical(Coords):
            def scalarSimp(self,x):
                return(4*x)

     ################### begin checks###############################            

        sp=Spherical()
        s.assertEqual(4*eye(2),sp.matSimp(eye(2)))



if  __name__ == '__main__':
    unittest.main()
