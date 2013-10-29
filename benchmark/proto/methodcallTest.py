#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
class A:
    def afunc(self,arg):
        print("I am afunc, and i got "+str(arg))
    def bfunc(self,arg):
        print("I am bfunc, and i got "+str(arg))
       
inst=A()
call=getattr(inst,"afunc")
call(5)
call=getattr(inst,"bfunc")
call(6)

