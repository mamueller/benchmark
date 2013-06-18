#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
class Coords(object):
    cV=50
    @staticmethod
    def sm():
        print("This is sm")
        print("I am a static method I can only acces variables of a specific class and not of subclasse") 
        print(Coords.cV)
        #print(????.subclassParameter)
    
    @classmethod
    def cm(clsname):
        print("This is cm")
        #print(__name__)
        print("Iam a class method ")
        print(clsname.cV)
        print(clsname.subclassParameter)

class Spherical(Coords):
    subclassParameter=51



Spherical.sm() 
Spherical.cm() 

