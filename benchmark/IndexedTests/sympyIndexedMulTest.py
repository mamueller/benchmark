from sympy.tensor.index_methods import get_contraction_structure
from sympy import symbols, default_sort_key
from sympy.tensor import IndexedBase, Idx
x, y = map(IndexedBase, ['x', 'y'])
i, j = map(Idx, ['i', 'j'])
x[i]*y[j]
