#!/usr/bin/python
from sympy import Symbol, cos,sqrt, series
h = Symbol('h')
c = Symbol('c')
r0 = Symbol('r0')
P=r0/sqrt(r0**2-c*h**2)
print(series(P,h))
DP=diff(P,h)
print(series(P,h))
