#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
class A:
    def afunc(self,arg):
        print("I am afunc, and i got "+str(arg))
    def bfunc(self,arg):
        print("I am bfunc, and i got "+str(arg))
       
inst=A()
caseStr="afunc"
call=getattr(inst,caseStr)
call(5)
caseStr="bfunc"
call=getattr(inst,caseStr)
call(6)

