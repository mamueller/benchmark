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
class ContractionIncompatibleBaseException(BaseException):
    def __init__(self,ind,p0,p1,b0,b1):
        self.ind=ind
        self.p=p0
        self.p=p1
        self.b=b0
        self.b=b1
    def __str__(self):
        t=Template("The index $ind occurs in positions $p0 and $p1 related to the bases $b0 and $b1 which are not dual to each other.")
        return(t.substitute( ind=selfind, p0=self.p0, p1=self.p1, b0=self.b0, b1=self.b1 ))


