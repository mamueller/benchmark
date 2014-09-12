#!/usr/bin/python3
# vim:set ff=unix expandtab ts=4 sw=4:
from string import Template


##########################################################
class DualBaseExeption(BaseException):
    pass
##########################################################
class BaseMisMatchExeption(BaseException):
    pass
##########################################################
class IncompatibleShapeException(BaseException):
    pass
##########################################################
class ContractionException(BaseException):
    def __init__(self,ind,n):
        self.ind=ind
        self.n=n
    def __str__(self):
        t=Template("The index $ind occurs n=$n times in the  expression. but we implement contraction only for n=2")
        return(t.substitute(ind=self.ind,n=self.n))
##########################################################
class ContractionIncompatibleBaseException(Exception):
    def __init__(self,ind,p0,p1,b0,b1):
        self.ind=ind
        self.p0=p0
        self.p1=p1
        self.b0=b0
        self.b1=b1
    def __str__(self):
        t=Template("The index $ind occurs in positions $p0 and $p1 related to the bases $b0 and $b1 which are not dual to each other.")
        #return "blub"
        return(t.substitute( ind=str(self.ind), p0=str(self.p0), p1=str(self.p1), b0=str(self.b0), b1=str(self.b1) ))
        


