from sympy import diff, Derivative, Expr, Pow, Basic, symbols
#from sympy.tensor import IndexedBase, Indexed, Idx
from MarkusIndexed import *
r=symbols("r")
bc=VectorFieldBase()
br=OneFormFieldBase(bc)
x=VIB("x",[br])
r,phi=symbols("r,phi")
x[0]=r**2
#x[1]=phi**2
i=Idx("i")
res=x[i].free_symbols
print(res)
#print(Derivative(x[i],r,evaluate=True))
#print((x[i])._eval_derivative(r))

