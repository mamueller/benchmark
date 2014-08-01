def permuteTuple(t,newPositions):
    nl=[t[p] for p in newPositions]
    return(tuple(nl))

def extractIndices(t,pos):
    nl=[t[p] for p in range(len(t)) if p in pos]
    return(tuple(nl))

def changedTuple(origTuple,newValTuple,newPositions):
    nl=list(origTuple)
    for p in range(len(newValTuple)):
        nl[newPositions[p]]=newValTuple[p]
    return(tuple(nl))

def deleteIndices(origTuple,pos):
    nl=[origTuple[p] for p in range(len(origTuple)) if not(p in pos)]
    return(tuple(nl))

def tupleLen(t):
        if type(t)==int:
           lk=1
        else:   
           lk=len(t)
        return(lk)
