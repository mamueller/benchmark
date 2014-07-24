#!/usr/bin/python3
from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
def tt():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
    g = Lorentz.metric
    p, q = tensorhead('p,q', [Lorentz], [[1]])
    p(m1)
    #t = p(m1)*g(m0,m2)
    #t.get_indices()
tt()
