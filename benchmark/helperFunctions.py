
# vim:set ff=unix expandtab ts=4 sw=4:
from sympy import *
###########################################################
def pp(string,gl):
    print("\n")	
    print(string,"=",eval(string,gl))
