#!/usr/bin/python
# vim: set expandtab ts=4
class O:
   def __init__(self,val):
      self.val=val

   def __eq__(self,other):
       return(self.val==other.val)

   def __mul__(self,rf):
      t=type(rf)
      if t==int or t==float or t==complex:
	 return(O(self.val*rf))
      if t==type(self):
	 return(O(self.val*rf.val))
      else:
	 ret(NotImplemented)

   def __rmul__(self,lf):
      return(self.__mul__(lf))

   def __str__(self):
      return(str(self.val))

a=O(1)
b=O(2)
print(a*b)
print(a*5)
print(6*a)
c=a
print(a)
print(c)
print(a==c)
c=O(1)
print(a)
print(c)
print(a==c)
import unittest 
class OTest(unittest.TestCase):
    def test_equal(self):
        a=O(1)
        c=O(1)
        self.assertEqual(a,c)

if  __name__ == '__main__':
   unittest.main()

