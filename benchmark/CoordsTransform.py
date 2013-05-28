#!/usr/bin/python
# vim: set expandtab ts=4
from sympy import *
import unittest 
import numpy as np
import copy
import Tensor 
from Exceptions import *

class CoordsTransform(object):
    '''This class represents a change of basis'''
    '''It is needed to compute the components of vectors and Tensors in a new
    coordinate frame'''
    def __init__(self,source,dest,mat):
        self.source=source
        self.dest=dest
        self.mat=mat
    def transform(self,tens,n):
        T=copy.deepcopy(tens)
        c=T.coords
        cT=T.componentTypes
        cc=T.components
        if cT[n]!=self.source:
            raise(TransformError(self.source,self.dest,cT,n))
	
        vec=Tensor.components2vec(cc)
        resvec=self.mat*vec
        resvec=c.matSimp(resvec)
        rescomponents=Tensor.vec2components(resvec)
        T.components=rescomponents
        T.componentTypes[n]=self.dest
	return(T)
    def invTransform(self,tens,n):
        T=copy.deepcopy(tens)
        c=T.coords
        cT=T.componentTypes
        cc=T.components
        if cT[n]!=self.dest:
            raise(TransformError(self.dest,self.source,cT,n))
	
        vec=Tensor.components2vec(cc)
        resvec=self.mat.inv()*vec
        resvec=c.matSimp(resvec)
        rescomponents=Tensor.vec2components(resvec)
        T.components=rescomponents
        T.componentTypes[n]=self.source
	return(T)
        

