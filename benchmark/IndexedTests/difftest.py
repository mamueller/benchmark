from sympy import diff, Derivative, Expr, Pow, Basic, symbols
from sympy.tensor import IndexedBase, Indexed, Idx
#class MarkusTest(Expr):
class MarkusTest(Indexed):
  def __new__(cls, b):
    obj = Expr.__new__(cls,b)
    return(obj)


  def _eval_derivative(self,v):
    return(5)

r=symbols("r")
a=MarkusTest(r**2)
print(type(a))
#print(diff(a,r))
print(Derivative(a,r,evaluate=True))
#e=r*r
#o=Pow(r,2)
#print(diff(e,r))
#print(Derivative(e,r,evaluate=True))
#print(Derivative(o,r,evaluate=True))

