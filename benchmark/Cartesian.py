#!/usr/bin/python
# vim: set expandtab ts=4
from Coords import *

class Cartesian(Coords):
    def __init__(self):
        # cartesian part
       
        # every general coordinate system is related to the cartesian basis
        x,y,z=symbols("x,y,z")
        X=x,y,z
	# the "new" coordinates are:
        ux,uy,uz=symbols("ux,uy,uz")
        U=ux,uy,uz
	# they are related to the "old" cartesian one by the identity transformation 
        x  =  ux
        y  =  uy
        z  =  uz
        XofU=[x,y,z]
	# every instance of our general coordinate systems has a function attached to it
	# that adds some rules to simplify results. It is called by nearly all the methods
	# before the result is returned. 
	# In the case of the cartesian example this function does nothing 
	sc=super(Cartesian,self)
	sc.__init__(X,U,XofU)

    def __repr__(self):
        return("Cartesian()")
   
