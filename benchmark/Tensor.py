#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:
import copy 
from sympy import *
from Exceptions import *
from Coords import *
from CoordsTransform import *
###########################################################
def indextupel(dim):
   keyset={(0,),(1,),(2,)}
   for j in range(1,dim):
      keyset=add_index(keyset)
   return(keyset) 

###########################################################
def del_index(keyset,pos):
    keyset=list(keyset)
    if len(keyset)==0:
        raise(EmptyIndexSetError)
    l=len(keyset[0]) #length of first tupol in the set
    if l <2 :
        raise(ToShortTupelError)
    newset=set()
    pos=pos%l
    for key in keyset:
        before=key[0:pos]
        after =key[pos+1:l]
        #print("\n key="+str(key)+'\n before='+ str(before)+"\n after="+ str(after))
        newkey=before+after
        newset.add(newkey)
    return(newset)   

###########################################################
def add_index(keyset,n=1):
    if n==0:
        return(keyset)
    else:
        newset=set()
        r=range(0,3)
        for l in range(0,n):
            for key in keyset:
                for i in r:
                    newkey=key+(i,)
                    newset.add(newkey)
        return(add_index(newset,n-1) )          

###########################################################
def indexTensProd(is1,is2):
    newset=set()
    for k1 in is1:
        for k2 in is2:
            newkey=k1+k2
            newset.add(newkey)
    return(newset)
               
###########################################################
def changedTupel(origTupel,position,newvalue):
    # tuples are immutable, therefor we have to convert to a list
    # first, change what we want to change and then convert back

    l=list(origTupel)
    l[position]=newvalue
    return(tuple(l))

	   
###########################################################
def components2vec(components):
    vec=zeros((3,1))
    for i in range(0,3):
       ind=(i,)
       if (ind in components.keys()): 
	   vec[i]=components[ind]
    print(vec)	   
    return(vec)

###########################################################
def vec2components(vec):
    comp={}
    for i in range(0,len(vec)):
        comp[(i,)]=vec[i]
    return(comp)
   
class Tensor(object):
###########################################################
    def __init__(self,coords,componentTypes,components):
        # test if we can handle the component types 
        if not(set(componentTypes).issubset({"roof","cellar","cart","phys"})):
            raise(UnknownComponentType(componentTypes))
        
        self.components={} # initialize empty (Null Tensor)

        indextupels=components.keys()
        
        # if any indextuples are given then perform some compatibility tests
        # otherwise proceed with the initialization of the Null tensor
        if len(indextupels)!=0:
            t0=indextupels[0]
            # for first order Tensors we allow itegers as indeces instead of tuples
            if isinstance(t0,int): 
                wrongtype=lambda y:not(isinstance(y,int))
                # make sure that all other indeces are integers to
                if any(map(wrongtype,indextupels)):
                    raise(IndexTupelError(indextupels))
                # make sure that it really is a first order Tensor
                if len(componentTypes)!=1:
                    raise(ComponentTypesTupelMismatch(componentTypes,t0))
            # for higher order tensors we expect tuples as indeces    
            # we test if all indextupels have the same length and type
            elif isinstance(t0,tuple):
                l0=len(t0)
                wronglength=lambda y:len(y)!=l0
                if any(map(wronglength,indextupels)):
                    raise(IndexTupelError(indextupels))
                if len(componentTypes)!=l0:
                    raise(ComponentTypesTupelMismatch(componentTypes,t0))
            else:
               raise(IndexTupelTypeError(indextupels))
            # if the tests were successful copy the components
            for it in indextupels:
               self.components[it]=components[it]
        
        self.coords=coords
        self.componentTypes=componentTypes
        r=range(0,coords.n)
        self.r=r
        # initialize coordinate transformations
        self.roof2cart=CoordsTransform("roof","cart",coords.J)
        self.cellar2cart=CoordsTransform("cellar","cart",coords.Jinv.transpose())

###########################################################
    def __str__(self):
        res="class:"+self.__class__.__name__+"\n"+"componentTypes="+str(self.componentTypes)+"\ncomponents="+str(self.components)+"\n"
        return(res)
###########################################################
    def __repr__(self):
        # this function should return an expression that when passed to eval
        # creates an identical copy of the object
        res="Tensor("+repr(self.coords)+","+repr(self.componentTypes)+","+repr(self.components)+")"
        return(res)
        
    
###########################################################
    def raise_index(self,pos):
        # raise index pos of the tensor
        # simmonds P39 2.9
        coords=self.coords
        g_rr=coords.g_rr()
        sc=self.componentTypes
        Res=self
        if sc[pos]!="cellar":
            raise(RaiseError(self))
        # e.g:
        #  j       ij
        # v   = v g
        #        i
        # e.g:
        #  j        ij
        # T   = T  g
        #  *l    il
        # e.g:
        #  *j        kj
        # T   = T   g
        #  l     lk  
        vc=self.components
        for key in vc.keys():
            j=key[pos]
            f=lambda i:changedTupel(key,pos,i)*g_rr[i,j]
            Res.components[key]=sum(map(f,self.r))
           
        Res.componentTypes[pos]="roof"
        return(Res)
    
###########################################################
    def raise_first_index(self):
        # raise index ind of the tensor
        # simmonds P39 2.9
        return(self.raise_index(0))
        
###########################################################
    def raise_last_index(self):
        # raise index ind of the tensor
        # simmonds P39 2.9
        return(self.raise_index(-1))
    
###########################################################
    def lower_index(self,pos):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        coords=self.coords
        g_cc=coords.g_cc()
        sc=self.componentTypes
        Res=self
        if sc[pos]!="roof":
            raise(LowerError(self))
        # e.g:
        #        i  
        # v   = v  g
        #  j        ij
    
        #  l     lk  
        # T   = T   g
        #  *j         kj
        
        #  *k    ik  
        # T   = T   g
        #  l         il
        vc=self.components
        for key in vc.keys():
            j=key[pos]
            f=lambda i:changedTupel(key,pos,i)*g_cc[i,j]
            Res.components[key]=sum(map(f,self.r))
           
        Res.componentTypes[pos]="cellar"
        return(Res)
    
###########################################################
    def lower_first_index(self):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        return(self.lower_index(0))
###########################################################
    def lower_last_index(self):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        return(self.lower_index(-1))
###########################################################
    def __add__(self,other):
        cself=copy.deepcopy(self)
        t=other.__class__.__name__
        #print("t="+str(t))
        if t==self.__class__.__name__:
            cto=other.componentTypes
            ctc=cself.componentTypes
            if len(cto)!=len(ctc):
                raise(ArgumentSizeError("+",(len(cto),len(ctc))))
            if cto==ctc:
                cc=cself.components
                co=other.components
                cok=co.keys()
                for k in cc.keys():
                    if not(k in cok):
                        co[k]=0
                    cc[k]+=co[k]
                return(cself)    
            else:        
                print("adding tensors has up to now only been implemented for "\
                      +"Tensors with the same component types. It could of course "\
                      +"be implemented by a implicit conversion to a common"\
                      +"component Type sequence")
                raise(NotImplementedError)
        else:
            raise(ArgumentTypeError("+",self.__class__.__name__,t))
###########################################################
    def __sub__(self,other):
       return(self+(-1)*other) 
    
                    
    
###########################################################
    def __rmul__(self,lf):
        cself=copy.deepcopy(self)
        # in case that the left factor lf in a product of the form lf*rf
        # does not belong to this class we swap the order
        return(cself.__mul__(lf))
    
###########################################################
    def __mul__(self,other):
        cself=copy.deepcopy(self)
        t=type(other)
        #print("t="+str(t))
        if t==int or t==float or t==complex:
            res=cself
            d=cself.components
            res.components={k:other*d[k] for k in d.keys()}
            return(res)
        scoords=cself.coords
        ocoords=other.coords
        if scoords!=ocoords:
            raise(NotImplementedError) 
        else:
            sc=cself.componentTypes
            oc=other.componentTypes
            print(sc)
            print(oc)
            sl=sc[-1]
            of=oc[0]
            print("sl,of",sl,of)
            if (sl=="roof" and of=="cellar") or (sl=="cellar" and of=="roof"):
                mat=cself.mat*other.mat
                l=sc[0:-1]+oc[1:]
                res=Tensor(scoords,l,mat)
            elif (sl=="cellar" and of=="cellar"):
                # raise first index of ohter tensor
                new=other.raise_first_index()
                res=cself.__mul__(cself,new)
            elif (sl=="roof" and of=="roof"):
                # lower first index of ohter tensor
                new=other.lower_first_index()
                res=cself.__mul__(new)
            else: 
               raise(NotImplementedError) 
            return(res)

###########################################################
    
    def vectorScalarProduct(self,other): #(scalar Product | )
        scT=self.componentTypes
        ocT=other.componentTypes
        sc=self.components
        oc=other.components
        sck=sc.keys()
        ock=oc.keys()
        commonKeys=set(sck).intersection(ock)
        ns=len(scT)
        no=len(ocT)
        if not(ns==1 and no==1):
            print("ns="+str(ns))
            print("no="+str(no))
            raise(NotImplementedError) 
        else:
            ll=scT[-1] #last left 
            fr=ocT[0]  #first right
            if (ll=="roof" and fr=="cellar") or (ll=="cellar" and fr=="roof") or (ll=="cart" and fr=="cart"):
                res=sum(map(lambda key:sc[key]*oc[key],commonKeys))
                res=self.coords.scalarSimp(sympify(res))
            else: 
                scart=self.transform2(["cart"])
                ocart=other.transform2(["cart"])
                res=scart.vectorScalarProduct(ocart)
        return(res)

###########################################################
    def tensorVectorScalarProduct(self,other): #(scalar Product | )
        # self is a tensor of arbitrary order
        # other is a vector (a tensor of first order)
        scT=self.componentTypes
        ocT=other.componentTypes
        sc=self.components
        sck=sc.keys()
        ns=len(scT)
        no=len(ocT)
        if not(ns>1 and no==1):
            print("ns="+str(ns))
            print("no="+str(no))
            raise(NotImplementedError) 
        left_keys=del_index(sck,-1)
        resComponents={}
        for lk in left_keys:
            # extract all tupels of the original lefthand side tensor 
            # that start with lk and build a first order tensor 
            # from it
            lvec=self.extractVector(lk+("*",))
            val=lvec.vectorScalarProduct(other)
            if val!=0:
                resComponents[lk]=val

        resComponentTypes=scT[0:-1]
        res=Tensor(self.coords,resComponentTypes,resComponents)
        return(res)

###########################################################
    def tensorTensorScalarProduct(self,other): #(scalar Product | )
        # self is a tensor of arbitrary order
        # other is a tensor of arbitrary order
        scT=self.componentTypes
        ocT=other.componentTypes
        sc=self.components
        oc=other.components
        sck=sc.keys()
        ock=oc.keys()
        ns=len(scT)
        no=len(ocT)
        if not(ns>1 and no>1):
            print("ns="+str(ns))
            print("no="+str(no))
            raise(NotImplementedError) 
        left_keys=del_index(sck,-1)
        right_keys=del_index(ock,0)
        resComponents={}
        for lk in left_keys:
            # extract all tupels of the original lefthand side tensor 
            # that start with lk and build a first order tensor 
            # from it
            lvec=self.extractVector(lk+("*",))
            
            for rk in right_keys:
                # extract all tupels of the original right hand side tensor 
                # that end with rk and build a first order tensor 
                # from it
                rvec=other.extractVector(("*",)+rk)
            val=lvec.vectorScalarProduct(rvec)
            if val!=0:
                resComponents[lk+rk]=val

            resComponentTypes=scT[0:-1]+ocT[1:]
        res=Tensor(self.coords,resComponentTypes,resComponents)
        return(res)

###########################################################
    def extractVector(self,tup):
        #use in scalar product
        scT=self.componentTypes
        l=len(scT)
        if l <2:
            raise(ArgumentSizeError(s=(l),op="extractVector"))
        sc=self.components
        sck=sc.keys()
        pos=tup.index("*")
        matchfun=lambda testtup: (testtup[0:pos]==tup[0:pos] and testtup[pos+1:]==tup[pos+1:])
        matchingKeys=filter(matchfun,sck)
        vec_components={}
        for m in matchingKeys:
            vec_components[(m[pos],)]=sc[m]
        vec=Tensor(self.coords,[scT[pos]],vec_components)
        return(vec)


###########################################################
    def __or__(self,other): #(scalar Product | )
        scoords=self.coords
        ocoords=other.coords
        if scoords!=ocoords:
          raise(NotImplementedError) 
        else:
            scT=self.componentTypes
            ocT=other.componentTypes
            sc=self.components
            oc=other.components
            sck=sc.keys()
            ock=oc.keys()
            ns=len(scT)
            no=len(ocT)
            #both tensors are vectors
            if (ns==1 and no==1):
                return(self.vectorScalarProduct(other))
                
            ##the left tensor has at least dimensio 2
            elif (ns>1 and no==1):
                return(self.tensorVectorScalarProduct(other))

            ##the right hand side tensor has at least dimensio 2
            elif (ns==1 and no>1):
                return(other.tensorVectorScalarProduct(self))
            
            ##both tensors have at least dimensio 2
            elif (ns>1 and no>1):
                return(self.tensorTensorScalarProduct(other))
            else: 
                raise(NotImplementedError) 

    
###########################################################
    def transform2(self,newComponentTypes):
        cs=copy.deepcopy(self)
        c=cs.coords
        csT=cs.componentTypes
        cc=cs.components
        if (csT==["roof"] and newComponentTypes==["cart"]):
            res=self.roof2cart.transform(self,0)
            return(res)
        if (csT==["cart"] and newComponentTypes==["roof"]):
            res=self.roof2cart.invTransform(self,0)
            return(res)
        if (csT==["cellar"] and newComponentTypes==["cart"]):
            res=self.cellar2cart.transform(self,0)
            return(res)
        if (csT==["cart"] and newComponentTypes==["cellar"]):
            res=self.cellar2cart.invTransform(self,0)
            return(res)
        
        if (csT==["cart","cart"] and newComponentTypes==["roof","roof"]):
            # for conviniece
            cc=self.roof2cart.invTransform(cs,0) #transform the first component, the result of this operation is the representation of mixed components, that is the result is given in the new basis while acting on vectors given in the old basis  
            cc=self.roof2cart.invTransform(cc,1) # transform the second component 
            
        
        
        return(cs)
    
###########################################################
    def div(self):
        sc=copy.deepcopy(self)
        cT=sc.componentTypes
        co=sc.coords
        l=len(cT)
        c=self.components
        n={}
        if (l==2):
            for k in range(0,3):
                f=lambda j:sc.covder(j,(j,k))
                n[(k,)]=sum(map(f,range(0,3)))
            return(Tensor(co,[cT[1]],n))
###########################################################
    def transpose(self):
        new=self
        coords=self.coords
        sc=self.componentTypes
        if len(sc)!=2:
          raise(NotImplementedError) 
        elif sc==["cellar","roof"]:
           new.mat=self.mat.transpose()
           new.componentTypes=["roof","cellar"]
           new=new.lower_first_index().raise_last_index()
        return(new)
###########################################################
    def str(self):
        description="coords="+str(self.coords)\
                +"componentTypes="+str(self.componentTypes)\
                +"components="+str(self.components)
        return(description)
###########################################################
    def __eq__(self,other):
        boolval=\
        type(self.coords)==type(other.coords) and \
        self.componentTypes==other.componentTypes and\
        self.components==other.components
        return(boolval)
    
###########################################################
    def subs(self,*args):
        sc=copy.deepcopy(self)
        c=sc.components
        n={}
        for k in c.keys():
            n[k]=c[k].subs(*args)
        sc.components=n
        return(sc)

###########################################################
    def partder(self,i):
        sct=self.componentTypes
        sco=self.coords
        sC=self.components
        keyset=indextupel(1)
        sm=components2vec(sC)
        rc={}
        for k in range(0,sco.n):
            #rc[(k,)]=sco.cov_der_v(i,k,sm)
            rc[(k,)]=self.covder(i,(k,))
        
        res=Tensor(sco,sct,rc)
        return(res)
###########################################################
    def covder(self,i,k):
        # this function computes the covariant derivative of a tensor
        # is the index of the partial derivative adn k is an index tupel
        # if the tensor is given in roof components v
        # (these are the components of v with respect to the >>cellar<< base vectors since v=vr[i]gc[i])
        # the return value is interpreted as a roof component 
        # of the partial derivative also 
        # if the tensor is given by its cellar components then the 
        # result is to be interpreted as a cellar component of the partial 
        # derivative
        sct=self.componentTypes
        sco=self.coords
        n=sco.n
        sC=self.components
        keys=sC.keys()
        if sct==["roof"]:#p.79 eq.4.17
            if not((k,) in keys):
                s=0
            else:
                s=diff(sC[(k,)],sco.U[i])
            for j in keys:
                    s+=sco.Gamma[k[0],i,j[0]]*sC[j]
        if sct==["cellar"]:#p.79 eq.4.19
            if not((k,) in keys):
                s=0
            else:
                s=diff(sC[(k,)],sco.U[i])
            for j in keys:
                    s=s-sco.Gamma[j[0],i,k[0]]*sC[j]
        if sct==["roof","roof"]:#p.98  eq 4.117
            j= k[0] #first component of the index
            ka=k[1] #second component of the index
            Gamma=sco.Gamma
            if not(k in keys):
                s=0
            else:
                s=diff(sC[k],sco.U[i])
            for p in range(0,n):
                s+=sco.Gamma[j,p,i]*sC[(p,ka)]
            for p in range(0,n):
                if (j,p) in keys:
                    #print(sC[(j,p)])
                    #print(Gamma[ka,p,i])
                    s+=Gamma[ka,p,i]*sC[(j,p)]
            return(simplify(sco.scalarSimp(s)))
        
        return(simplify(sco.scalarSimp(s)))
            


###########################################################
    def nabla(self):
        # This function applies the Nabla operator to the tensor.
        cs=copy.deepcopy(self)
        cst=cs.componentTypes
        #            i
        # nabla v = g * v,i (* = "direct Product")
        keys=cs.components.keys()
        components={}
        t=indextupel(1)
        for i in t:
            comp[i]=cs.partder(i)
        res=Tensor(cs.coords,"cellar",comp)
