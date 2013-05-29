from string import Template
##########################################################
class ComponentsLengthError(Exception):
    def __init__(self,components):
        self.components=components
    def __str__(self):
        return("The initialize method of class Vector accepts exactly one component type , but got: "+str(components))
##########################################################
class ComponentTypesTupelMismatch(Exception):
    def __init__(self,tupel,types):
        self.tupel=tupel
        self.types=types
    def __str__(self):
       return("The indextupels and list of component types differ in length:\n The types are: "+str(self.types)+"while the first tupel is:"+self.tupel+".\n")

##########################################################
class UnknownComponentType(Exception):
    def __init__(self,t):
        self.t=t
    def __str__(self):
        return("The component type : "+str(self.t)+" is not supported.")

##########################################################
class TransformError(Exception):
    def __init__(self,source,dest,componentTypes,n):
        self.source=source 
        self.dest=dest 
        self.componentTypes=componentTypes 
        self.n=n 
    def __str__(self):
        source		=self.source 
        dest		=self.dest 
        componentTypes	=self.componentTypes 
        n		=self.n 
	
        return("The Transformation cannot be applied because it tranforms from :"+str(self.source)+" to "+str(self.dest)+" while the "+str(self.n)+" th component of "+str(self.componentTypes) +" is \'"+ str(self.componentTypes[n])+"\'.")
#	return("The Transformation cannot be applied because it tranforms from :"+str(source)+" to "+str(dest)) 

##########################################################
class RaiseError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
        return("The index is a roof index already and connot be raised: "+str(self.d))

##########################################################
class EmptyIndexSetError(Exception):
    def str(self):
        return("The indexset is empty")

##########################################################
class ToShortTupelError(Exception):
    def str(self):
        return("The tupel is to small to delete one index from it")


##########################################################
class IndexTupelError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
        return("The indextupels have different length: "+str(self.d))
##########################################################
class IndexTupelTypeError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
       return("The indices must be either integers or tupels of integers. Instead got this:: "+str(self.d))


##########################################################
class LowerError(Exception):
    def __init__(self,d):
        self.d=d
    def __str__(self):
        return("The index is a cellar index already and cannot be lowered: "+str(self.d))

##########################################################
class NotImplementedError(Exception):
    def __str__(self):
        return("The feature is not implemented yet: ")
##########################################################
class ArgumentTypeError(Exception):
    def __init__(self,cl1,cl2,op):
        self.cl1=cl1
        self.cl2=cl2
        self.op=op
    def __str__(self):
        t=Template("The operation $op is not supported for classes $c1 and $c2.")
        return(t.substitute(op=self.op,c1=self.cl1,c2=self.cl2))

##########################################################
class ArgumentSizeError(Exception):
    def __init__(self,s1,s2,op):
        self.s1=s1
        self.s2=s2
        self.op=op
    def __str__(self):
        t=Template("The operation $op is not supported for sizes $s1 and $s2.")
        return(t.substitute(op=self.op,s1=self.s1,s2=self.s2))

