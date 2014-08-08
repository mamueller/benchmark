from sympy import diff, Derivative, Expr, Pow, Basic, symbols
class MarkusTest(Expr):
  def __new__(cls, b, evaluate=None):
    obj = Expr.__new__(cls,b)
    return(obj)


  def _eval_derivative(self,v):
    return(5)

r=symbols("r")
e=r*r
o=Pow(r,2)
a=MarkusTest(r)
print(type(a))
#print(diff(a,r))
#print(diff(e,r))
print(Derivative(a,r,evaluate=True))
print(Derivative(e,r,evaluate=True))
print(Derivative(o,r,evaluate=True))

