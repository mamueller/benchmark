from sympy import diff, Derivative, Expr, Pow, Basic, symbols
from sympy.tensor import IndexedBase, Indexed, Idx
from MarkusIndexed import VIB, OIB, VI ,IncompatibleShapeException,VectorFieldBase,OneFormFieldBase
r=symbols("r")
bc=VectorFieldBase()
br=OneFormFieldBase(bc)
x=VIB("x",[br])
r=symbols("r")
x[0]=r**2
i=Idx("i")
print(Derivative(x[i],r,evaluate=True))
print((x[i])._eval_derivative(r))

