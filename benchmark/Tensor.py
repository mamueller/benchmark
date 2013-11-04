#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:

import copy 
from sympy import *
from helperFunctions import *
import Exceptions 

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
        raise(Exceptions.EmptyIndexSetError)
    l=len(keyset[0]) #length of first tupol in the set
    if l <2 :
        raise(Exceptions.ToShortTupelError)
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
def exchangeEnds(origTupel):
    fe=origTupel[0]  #first index
    le=origTupel[-1] #last index
    l=list(origTupel)
    l[0]=le
    l[-1]=fe
    return(tuple(l))
	   
###########################################################
def components2vec(components):
    vec=zeros((3,1))
    for i in range(0,3):
       ind=(i,)
       if (ind in components.keys()): 
	   vec[i]=components[ind]
    return(vec)

###########################################################
def vec2components(vec):
    comp={}
    for i in range(0,len(vec)):
        comp[(i,)]=vec[i]
    return(comp)
###########################################################
###########################################################
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
            raise(Exceptions.TransformError(self.source,self.dest,cT,n))
	
        vec=components2vec(cc)
        resvec=self.mat*vec
        resvec=c.matSimp(resvec)
        rescomponents=vec2components(resvec)
        T.components=rescomponents
        T.componentTypes[n]=self.dest
	return(T)
    def invTransform(self,tens,n):
        T=copy.deepcopy(tens)
        c=T.coords
        cT=T.componentTypes
        cc=T.components
        if cT[n]!=self.dest:
            raise(Exceptions.TransformError(self.dest,self.source,cT,n))
	
        vec=components2vec(cc)
        resvec=self.mat.inv()*vec
        resvec=c.matSimp(resvec)
        rescomponents=vec2components(resvec)
        T.components=rescomponents
        T.componentTypes[n]=self.source
	return(T)
###########################################################
###########################################################
class Coords(object):
    def __init__(self,X,U,XofU):
        self.X=X
        self.U=U
        self.XofU=XofU
        self.n=len(self.X)
        self.g_cc_=False
        self.g_rr_=False
        self.sf_=False
        self.setup()
    def scalarSimp(self,exp):
    	return(exp) #do nothing but return the input

    def matSimp(self,Mat):
        s=Mat.shape
        res=zeros(s)
        for i in range (s[0]):
            for j in range (s[1]): 
                res[i,j]=self.scalarSimp(Mat[i,j])
        return(res)

    def setup(self):
        # Matrix of the derivative dX/dU (Jacobi Matrix)
        def dd(i,j):
            return(diff(self.XofU[i],self.U[j]))
        
        self.J=Matrix(self.n,self.n,dd)
        
        self.Jinv=self.matSimp(self.J.inv())

        #
        # cellar Base vectors:
        self.gc={}
        self.t_gc={}
        for i in range(0,self.n):
            self.gc[i]=self.J[:,i]
            self.t_gc[i]=Tensor(self,["cart"],vec2components(self.gc[i]))

	    #define the cartesian base vectors e_x e_y and e_z
        
        #e_r_c    =self.t_gc[0]
        #e_phi_c  =self.t_gc[1]
        #e_theta_c=self.t_gc[2]
        
        
        # Jacobian (Determinant) of Jacobi Matrix
        detJ=self.J.det()
        # help sympy to use the trigonometric pythagoras
        r,phi,theta=self.U
        self.detJ=self.scalarSimp(detJ)
        
        # the inverse roof basis
        # The inverse base vectors are the columns of the transposed inverse Jacobi Matrix
        self.Jinv[0,0]=self.scalarSimp(self.Jinv[0,0])
        self.gr={}
        self.t_gr={}
        for i in range(0,self.n):
            self.gr[i]=self.Jinv[i,:].transpose()
            self.t_gr[i]=Tensor(self,["cart"],vec2components(self.gr[i]))
        
        #e_r_r    =self.gr[0]
        #e_phi_r  =self.gr[1]
        #e_theta_r=self.gr[2]
        

        # the scale factors (lengts of the basevectors)
        # roof
        self.Hroof=zeros(self.n,self.n)
        for i in range(0,self.n):
            self.Hroof[i,i]=self.gr[i].norm(2)
        #pp("self.Hroof",locals())
        # cellar
        self.Hcellar=zeros(self.n,self.n)
        for i in range(0,self.n):
            self.Hcellar[i,i]=self.gc[i].norm(2)



        # the cristoffel symbols
        # 1.) compute the derivatives of the cellar base vectos
        dgc={}
        # first create a function that differentiates matrices
        def mdiff(Mat,sym):
            s=Mat.shape
            res=zeros(s)
            for i in range(0,s[0]):
                for j in range(0,s[1]):
                    res[i,j]=diff(Mat[i,j],sym)
            return(res)        
        
        for i in range(0,self.n):
            for j in range(0,self.n):
                dgc[i,j]=mdiff(self.gc[i],self.U[j])
        
        # 2.) dgr[i,j]=Gamma[1,ij*self.gc[1]+Gamma[2,ij*self.gc[2]+Gamma[3,ij*self.gc[3]
        # multiply with the roof basic vectors to extract the Gammas
        self.Gamma={}
        for i in range(0,self.n):
            for j in range(0,self.n):
                 for k in range(0,self.n):
                     self.Gamma[k,i,j]=simplify(dgc[i,j].dot(self.gr[k]))
        
        self.roof2cart=CoordsTransform("roof","cart",self.J)
        self.cellar2cart=CoordsTransform("cellar","cart",self.Jinv.transpose())
        # note that for the physical roof components the length of the cellar base is needed and vice versa
        self.roof2roof_phys=CoordsTransform("roof","roof_phys",self.Hcellar)
        self.cellar2cellar_phys=CoordsTransform("cellar","cellar_phys",self.Hroof)

    def part_der_v(self,i,r_tup):
        # this function assumes a tupel r_tup representing the roof_components of a vector v
        # (these are the components of v with respect to the >>cellar<< base vectors since v=vr[i]gc[i])
        # it computes the roof components of the partial derivative of the vector v,sym 
        # with respect to the cellar base vectors
        s=r_tup.shape
        # insist on a vector since 
        if s[1]!=1: 
            raise(Exceptions.BaseException)
        res=zeros(s)
        for k in range(0,self.n):
            res[k]=self.cov_der_v(i,k,r_tup)
        return(res)

    def g_rr(self):
        # roof components of the Identity tensor
        if(self.g_rr_):
            return(self.g_rr_)
        self.g_rr_={}
        r=range(0,self.n)
        for i in r:
            for j in r:
                self.g_rr_[(i,j)]=self.scalarSimp(self.gr[i].dot(self.gr[j]))
        return(self.g_rr_)
    def g_cc(self):
        # cellar components of the metric tensor
        if(self.g_cc_):
            return(self.g_cc_)
        self.g_cc_={}
        r=range(0,self.n)
        for i in r:
            for j in r:
                self.g_cc_[(i,j)]=self.scalarSimp(self.gc[i].dot(self.gc[j]))
        return(self.g_cc_)

    def cc2cr(self,cc):
        # This function returns the mixed components 
        #    *j   
        #  T     
        #    i
        # of a second order tensor 
        #        *j   i
        # T=   T     g   g
        #        i       j
        # which is currently given by its cellar cellar Components
        # in this case the second index has to be raised which is done
        #                                            k 
        # by a dot product with the identity tensor g g  
        #                                              k 
        #                                        i  j        i  j    k
        # in component form this reads T   =T   g  g  = T   g  g  * g  g
        #                                    ij          ij             k
        # 
        #                                                    i  jk
        #                                             = T   g  g      g
        #                                                ij            k
        #
        #                                                *k   i 
        #                                             = T   g   g
        #                                                i*       k
        # 
        #    *k       jk  
        # ->T  = T   g    
        #    i*   ij
        # This is a matrix product cc* self.g_rr
        cr=cc*self.g_rr()
        return(cr)

    def cc2rc(self,cc):
        # This function returns the mixed components 
        #    i   
        #  T     
        #    *j
        # of a second order tensor 
        #        i        j
        # T=   T     g   g
        #        *j   i     
        # which is currently given by its cellar cellar Components
        # in this case the first index has to be raised which is done
        #                                              k 
        # by a dot product with the identity tensor g g  
        #                                            k   
        #                                        i  j       k       i  j   
        # in component form this reads T   =T   g  g  = g  g * T   g  g  
        #                                    ij          k      ij         
        # 
        #                                                        ki    j
        #                                             = T   g   g     g
        #                                                ij  k          
        #
        #                                             !  k       j
        #                                             = T   g   g
        #                                                *j  k    
        # 
        #    k    ki      
        # ->T  = g  T     
        #    *j      ij
        # This is a matrix product cc* self.g_rr
        cr=self.g_rr()*cc
        return(cr)
    def cr2rr(self,cr):
        # This function returns the roof roof components 
        #    ij   
        #  T     
        #    
        # of a second order tensor 
        #        ij   
        # T=   T     g   g
        #             i   j
        # which is currently given by its cellar roof Components cr
        # in this case the first index has to be raised which is done
        #                                               k
        # by a dot product with the identity tensor g  g 
        #                                            k   
        #                                    *j  i          k   *j    i      
        # in component form this reads T   =T   g  g  = g  g  *T     g  g   
        #                                    i      j    k       i        j    
        # 
        #                                                *j     ki
        #                                             = T   g  g      g
        #                                                i   k          j
        #
        #                                                kj    
        #                                             = T   g   g
        #                                                    k   j
        # 
        #    kj   *j  ki     ki  *j
        # ->T  = T   g    = g   T
        #          i             i
        # This is a matrix product self.g_rr*cc 
        cr=self.g_rr()*cr
        return(cr)
    def rr2CartCart(self,rr):
        #  This function computes the cartesian components of second order Tensor given  
        #  by its RoofRoof components
        #          ij 
        #  T=     T   g g
        #              i j
        # The cellar base vectors are the columns of the Jacobian
        # with  we g = J[k,j]e
        #           j         k
        #       ij
        # -> T=T  J[k,i] e  J[l,j]e
        #                 k        l
        #     kl             ij
        # -> T    J[k,i] * T   * J[l,j]
        # expressed as matrix multiplication this is
        #        T
        #  J*rr*J where rr is the component set written as matrix
        res=self.matSimp(self.J*rr*self.J.transpose())
        return(res)
    def cart_cart2rr(self,cartcart):
        # Given the cartesian components of a second order tensor this function
        # computes its spherical roof roof components
        #          ij 
        #  T=     T   e e
        #              i j
        # We can imagine the application of T to a vector
        # given by its cellar components v_c as a composition
        #
        # 1.) we compute the cartesian components v_cart of the vector v 
        #     using the same procedure as function cellar2cart
        # 2.) we apply T to v using the cartesian representations
        #     and matrix multiplication :
        #     res_cart=cartcart*v_cart  
        # 3.) Compute the roof components of the result which is 
        #     up to now given by its cartesian components 
        #     this is done as in function cart2roof
        # 
        # Combining the 3 steps we can write the resulting components as
        # matrix product
        # 
        res=self.Jinv*cartcart*self.Jinv.transpose()
        res=self.matSimp(res)
        # 

        return(res)

    def cellar_nab(self,f):
        # This function computes the cellar components of the gradient of the
        # scalar function f
        g=zeros([self.n,1])
        for i in range(0,self.n):
            g[i]=diff(f,self.U[i])
        return(g)


                


        
    def cellar2cart(self,cellarComponents):
        r,theta, phi = self.U
        # This function computes the cartesian Components of a vector given 
        # the cellar component of a vector (that is the components according to the 
        # roof basis vectors
        res=self.Jinv.transpose()*cellarComponents
        res=self.matSimp(res)
        return(res)
    def cart2cellar(self,cartComponents):
        r,theta, phi = self.U
        # compute the cellar components with respect to the reciproke basis 
	# from the cartesian components
        res=self.J.transpose()*cartComponents
        res=self.matSimp(res)
        return(res)

    def cart2roof(self,cartComponents):
        r,phi,theta=self.U
        #                                             i             i
        # This function computes the roof components v of a vector v g
        #                                                             i
	# from its cartesian components
	# 
        res=self.Jinv*cartComponents
        res=self.matSimp(res)
        return(res)
	



    def roof2cart(self,roofComponents):
        # This function computes the cartesian Components of a vector given 
        # the roof component of a vector (that is the components according to the 
        # cellar basis vectors
        res=self.J*roofComponents
        res=self.matSimp(res)
        return(res)

    def roof2cellar(self,rc):
        # This function converts the roof Components of a vector to its cellar 
        # components. 
        res=self.g_cc()*rc
        res=self.matSimp(res)
        return(res)
    

    def cellar2roof(self,cellarComponents):
        # This function converts the cellar Components of a vector to its roof 
        # components. 
        res=self.g_rr()*cellarComponents
        res=self.matSimp(res)
        return(res)
   
    def cellar2phys(self,cellarComponents):
        rC=self.cellar2roof(cellarComponents)
        pc=self.roof2phys(rC)
        return(pc)


    def roof2phys(self,roofComponents):
        s=self.scale_factors()
        rc=zeros([self.n,1])
        for i in range(0,self.n):
            rc[i]=roofComponents[i]*s[i]
        return(rc)

    def phys2roof(self,physComponents):
        s=self.scale_factors()
        rc=zeros([self.n,1])
        for i in range(0,self.n):
            rc[i]=physComponents[i]/s[i]
        return(rc)
    
    def phys2cellar(self,physComp):
        rC=self.phys2roof(physComp)
	cC=self.roof2cellar(rC)
	return(cC)


    def scale_factors(self):
        if self.sf_:
            return(self.sf_)
        #compute the lengths of the cellar base vectors
        r,theta, phi = self.U
        s=zeros([self.n,1])
        for i in range(0,self.n):
            s[i]=sqrt(self.g_cc()[i,i])
        s=self.matSimp(s)
        for i in range(0,self.n):
            s[i]=s[i].subs({
                sqrt(sin(phi)**2):sin(phi),
                sqrt(sin(theta)**2):sin(theta)
            })
        self.sf_=s    
        return(s)
    def phys2cart(self,physComponents):
        rc=self.phys2roof(physComponents)
        cc=self.roof2cart(rc)
        return(cc)

    def cart2phys(self,Components):
        # compute the physical components
        # (the components with respect to the normalized cellar base vectors)
        # of a vector from the cartesian compontens
        rc=self.cart2roof(Components)
        pc=self.roof2phys(rc)
        return(pc)

    def cellarComponentsOfNablaOnRoofComponents(self,rc):
        #                                                             *j   
        # This function computes the mixed (cellar_roof) components Nv    of the 
        #                                                             i*
        # tensor Nabla v from the roof components rc of the vector v
        # That means that Nabla v is given as linear combination of the direct
        #           i
        # products g  g
        #              j
        # Nv={}
	#print("rc="+str(rc))
        Nv=zeros(3,3)
        for i in range(0,self.n):
            for j in range(0,self.n):
                Nv[i,j]=self.cov_der_v(i,j,rc)
                #print("Nv["+str(i)+","+str(j)+"]="+str(Nv[i,j]))
	i=0;j=0
	#print("cov_der_v("+str(i)+","+str(j)+","+str(rc)+")="+str(self.cov_der_v(i,j,rc)))
        return(Nv)





    def roof_div(self,vr):		
	# compute the divergence of a vector given by its roof_components 
    # (these are the component with respect to the cellar base vectors since v=vr[i]*gc[i]
	res=0
	for i in range(0,self.n):
	    res+=self.cov_der_v(i,i,vr)
	return(simplify(res))    
    
#    def roofroof_div_T(self,Trr):		
#	# compute the divergence of a tensor T given by its roofroof_components Trr
#        # (these are the component with respect to the cellar base vectors since 
#	# T=Trr[i,j] gc[i] gc[j])
#	# (The result are the roof components of a vector 
#	# (with resprect to the cellar basis))
#	res=zeros([3,1])
#	for k in range(0,self.n):
#	    rv=0 # the value of the row of the result
#	    for i in range(0,self.n):
#	        rv+=self.cov_der_T(i,i,k,Trr)
#		rv=simplify(rv)
#	    res[k]=rv	
#	return(res)    
#
    def phys_div(self,vr):		
    # compute the divergence of a vector given by its physical componets 
        # compute the roof components
        rc=self.phys2roof(vr)
        # apply the divergence for roof_components 
        res=self.roof_div(rc)
        return(res)



    

	
	
#    def cov_der_v(self,i,k,vr):
#        # this function computes the covariant derivative of the roof components of a vector v
#        # (these are the components of v with respect to the >>cellar<< base vectors since v=vr[i]gc[i])
#        s=diff(vr[k],self.U[i])
#        for j in range(0,self.n):
#            s+=self.Gamma[k,i,j]*vr[j]
#        return(simplify(s))
#    
#    def cov_der_T(self,i,j,k,Trr):
#        # this function computes the covariant derivative of the roof components Trr of a tensor T
#        # (these are the components of T with respect to the >>cellar<< base vectors since T=Tr[i,j]gc[i]gc[j])
#        s=diff(Trr[j,k],self.U[i])
#        for p in range(0,self.n):
#            s+=self.Gamma[j,p,i]*Trr[p,k]
#        for p in range(0,self.n):
#            s+=self.Gamma[k,p,i]*Trr[j,p]
#        return(simplify(s))
   
##            raise NumTestError(maxd)
#    def laplaceBeltrami(self,fundef):
#        r,theta, phi = symbols("r theta phi")
#        f=fundef(theta,phi)
#        theta_part=1/(sin(theta))*diff(sin(theta)*diff(f,theta),theta)
#        phi_part=1/((sin(theta))**2)*diff(diff(f,phi),phi)
#        lb=Lambda((theta,phi),theta_part+phi_part)
#        return(lb)
#    def div(self,A):
#        # this function computes the divergence of a vector given in physical
#        # components.
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        div=1/r**2*diff(r**2*Ar,r) +1/(r*sin(theta))*diff(sin(theta)*Atheta,theta) +1/(r*sin(theta))*diff(Aphi,phi)
#        return(div)
#    def vectorLaplace(self,A):
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        res=self.grad(self.div(A))-self.rot(self.rot(A))
#        return(res)
#    def rot(self,A):
#        r,theta, phi = symbols("r theta phi",real=True)
#        Ar=A[0]
#        Aphi=A[1]
#        Atheta=A[2]
#        rr=1/(r*sin(theta))*(diff(Aphi*sin(theta),theta)-diff(Atheta,phi))
#        rphi=1/r*(diff((r*Atheta),r)-diff(Ar,theta))
#        rtheta=1/r*(1/sin(theta)*diff(Ar,phi)-diff(r*Aphi,r))
#        def res(r_val,phi_val,theta_val):
#            subsdict={r:r_val,phi:phi_val,theta:theta_val}
#            mat=Matrix(3,1,[rr.subs(subsdict),rphi.subs(subsdict),rtheta.subs(subsdict)])
#            return(mat)
#        return(res)
#
#class NumTestError(Exception):
#    def __init__(self,d):
#        self.d=d
#    def __str__(self):
#        return("The difference was: "+str(self.d))
#
#
#def testnumrad_eq(expr1,expr2,sym,r_min,r_max,prec):
#    #print("expr1",expr1)
#    #print("expr2",expr2)
#    def f1(x):
#        return((expr1.subs(sym,x)).evalf())
#    def f2(x):
#        return((expr2.subs(sym,x)).evalf())
#    args=np.linspace(r_min,r_max,50)
#    vals1=list2numpy(map(f1,args))
#    vals2=list2numpy(map(f2,args))
#    test=(vals1-vals2)/(vals1+vals2)
#    maxd=max(map(abs,test))
#    if maxd >=prec:
#        raise NumTestError(maxd)
#
    
###########################################################
###########################################################
class Scalar(object):
    # zero order tensors are tensors too 
    # Unlike the arguments of an arbitrary function 
    # the arguments of a scalar have a spatial meaning 
    # they describe a >>location<<.
    # If the representation of this location (the coordinates) 
    # is changed then the way the Scalar is computed from the coordinates
    # has to change too to make it only dependent on >>location<<

    def __init__(self,coords,definition):
        self.coords=coords
        self.definition=definition
   ###########################################################
    def nabla(self):
        # This function applies the Nabla operator (always from left) 
        # to the tensor.
        csc=self.coords
        gr=csc.t_gr
        #            i
        # nabla v = g * v,i (* = "direct Product")
        f=lambda i :Tensor(csc,["cellar"],{(i,):1})*diff(self.definition,csc.U[i])
        res=sum(map(f,range(0,csc.n)))
        return(res)
        
###########################################################
###########################################################
class Tensor(object):
    ##########################################################
    def __init__(self,coords,componentTypes,components):
        # test if we can handle the component types 
        if not(set(componentTypes).issubset({"roof","cellar","cart","phys"})):
            raise(Exceptions.UnknownComponentType(componentTypes))
        
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
                    raise(Exceptions.IndexTupelError(indextupels))
                # make sure that it really is a first order Tensor
                if len(componentTypes)!=1:
                    raise(Exceptions.ComponentTypesTupelMismatch(componentTypes,t0))
            # for higher order tensors we expect tuples as indeces    
            # we test if all indextupels have the same length and type
            elif isinstance(t0,tuple):
                l0=len(t0)
                wronglength=lambda y:len(y)!=l0
                if any(map(wronglength,indextupels)):
                    raise(Exceptions.IndexTupelError(indextupels))
                if len(componentTypes)!=l0:
                    raise(Exceptions.ComponentTypesTupelMismatch(componentTypes,t0))
            else:
               raise(Exceptions.IndexTupelTypeError(indextupels))
            # if the tests were successful copy the components
            for it in indextupels:
               self.components[it]=components[it]
        
        self.coords=coords
        self.componentTypes=componentTypes
        r=range(0,coords.n)
        self.r=r
        self.purge()
    ##########################################################
    def purge(self):
        # If a component is zero it can be safely erased
        # this is important to  make vectors and tensors
        # with zero components compareable
        c=self.components
        for k in c.keys():
            if c[k]==0:
               del(self.components[k])

    ##########################################################
    def __str__(self):
        res="class:"+self.__class__.__name__+"\n"+"componentTypes="+str(self.componentTypes)+"\ncomponents="+str(self.components)+"\n"
        return(res)
    ##########################################################
    def __repr__(self):
        # this function should return an expression that when passed to eval
        # creates an identical copy of the object
        res="Tensor("+repr(self.coords)+","+repr(self.componentTypes)+","+repr(self.components)+")"
        return(res)
        
    
    ##########################################################
    def raise_index(self,pos):
        # raise index pos of the tensor
        # simmonds P39 2.9
        coords=self.coords
        g_rr=coords.g_rr()
        sc=self.componentTypes
        Res=self
        if sc[pos]!="cellar":
            raise(Exceptions.RaiseError(self))
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
    
    ##########################################################
    def raise_first_index(self):
        # raise index ind of the tensor
        # simmonds P39 2.9
        return(self.raise_index(0))
        
    ##########################################################
    def raise_last_index(self):
        # raise index ind of the tensor
        # simmonds P39 2.9
        return(self.raise_index(-1))
    
    ##########################################################
    def lower_index(self,pos):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        coords=self.coords
        g_cc=coords.g_cc()
        sc=self.componentTypes
        Res=self
        if sc[pos]!="roof":
            raise(Exceptions.LowerError(self))
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
    
    ##########################################################
    def lower_first_index(self):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        return(self.lower_index(0))
    ##########################################################
    def lower_last_index(self):
        # lower index ind of the tensor self
        # simmonds P39 2.9
        return(self.lower_index(-1))
    ##########################################################
    def __add__(self,other):
        cself=copy.deepcopy(self)
        other=copy.deepcopy(other)
        t=other.__class__.__name__
        if t==self.__class__.__name__:
            cto=other.componentTypes
            ctc=cself.componentTypes
            if len(cto)!=len(ctc):
                raise(Exceptions.ArgumentSizeError("+",(len(cto),len(ctc))))
            if cto==ctc:
                cc=cself.components
                cck=cc.keys()
                co=other.components
                cok=co.keys()
                for k in set(cck).union(cok):
                    #if k in cck and not(k in cok): 
                        #do nothing
                    if (k in cck) and (k in cok): 
                        # sum components
                        cc[k]+=co[k]
                    elif (not(k in cck)) and (k in cok): 
                        # add new component 
                        cc[k]=co[k]
                    #elif k in cck and not(k in cok): 
                        #do nothing because cc[k] does not change
                   
                cself.components=cc
                return(cself)    
            else:        
                str="adding tensors has up to now only been implemented for "\
                      +"Tensors with the same component types. It could of course "\
                      +"be implemented by a implicit conversion to a common"\
                      +"component Type sequence"
                raise(Exceptions.NotImplementedError(str))
        else:
            raise(Exceptions.ArgumentTypeError(self.__class__.__name__,t,"+"))
    ###########################################################
    def __radd__(self,other):
        cs= copy.deepcopy(self)
        # this is necessary to make sum() work
        # it only calls >>add<< but in different oder
        # and it filters the first summand >>0<< that is
        # introduced by >>sum<< which is a proplem since 
        # >>add<< for Tensors does not accept integers of course
        
        if other==0:
            return(cs)
        else:    
            return(cs.__add__(other))
    
    #########################################################
    def __sub__(self,other):
       return(self+(-1)*other) 
    
                    
    
    ##########################################################
    def __rmul__(self,lf):
        cself=copy.deepcopy(self)
        # in case that the left factor lf in a product of the form lf*rf
        # does not belong to this class we swap the order
        return(cself.__mul__(lf))
    
    ##########################################################
    def __mul__(self,other):
        t=type(other)
        #print("t="+str(t))
        cs=copy.deepcopy(self)
        if t==type(self):
            scoords=self.coords
            ocoords=other.coords
            if scoords!=ocoords:
                raise(Exceptions.CoordsMissMatchError( scoords, ocoords))
            
            co=copy.deepcopy(other)
            sc=cs.components
            sck=sc.keys()
            #print("sck",sck)
            sct=cs.componentTypes
            oc=co.components
            ock=oc.keys()
            #print("ock",ock)
            oct=co.componentTypes
            nct=sct+oct #new component Types
            nc={}#new components

            for k1 in sck:
                for k2 in ock:
                    newkey=k1+k2
                    nc[newkey]=sc[k1]*oc[k2]
            res=Tensor(cs.coords,nct,nc)
            return(res)
        else:
            # one of the factor is a number or an expression
            res=cs
            d=self.components
            res.components={k:other*d[k] for k in d.keys()}
            return(res)

    ##########################################################
    
    def vectorScalarProduct(self,other): #(scalar Product | )
        sco=self.coords
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
            raise(Exceptions.NotImplementedError) 
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
        return(sco.scalarSimp(res))

    ##########################################################
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
            raise(Exceptions.NotImplementedError) 
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

    ##########################################################
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
            raise(Exceptions.NotImplementedError) 
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

    ##########################################################
    def extractVector(self,tup):
        #use in scalar product
        scT=self.componentTypes
        l=len(scT)
        if l <2:
            raise(Exceptions.ArgumentSizeError(s=(l),op="extractVector"))
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

    ##########################################################
    def simp(self):
        c=self.components
        co=self.coords
        for k in c.keys():
            c[k]=co.scalarSimp(c[k])
        self.purge()    

    ##########################################################
    def __or__(self,other): #(scalar Product | )
        scoords=self.coords
        ocoords=other.coords
        if type(scoords)!=type(ocoords):
            raise(Exceptions.NotImplementedError) 
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
                res=self.vectorScalarProduct(other)
                
            ##the left tensor has at least dimensio 2
            elif (ns>1 and no==1):
                res=self.tensorVectorScalarProduct(other)

            ##the right hand side tensor has at least dimensio 2
            elif (ns==1 and no>1):
                res=other.tensorVectorScalarProduct(self)
            
            ##both tensors have at least dimensio 2
            elif (ns>1 and no>1):
                res=self.tensorTensorScalarProduct(other)
            else: 
                raise(Exceptions.NotImplementedError) 
            return(res)

    
    ##########################################################
    def transform2(self,newComponentTypes):
        cs=copy.deepcopy(self)
        c=cs.coords
        csT=cs.componentTypes
        cc=cs.components
        # we iterate over the pairs of src and target component types
        # and treat every pair separately
        for i in range(0,len(csT)):
            caseString=csT[i]+"2"+newComponentTypes[i]
            print(caseString)

            if caseString=="roof2cart":
                cs=c.roof2cart.transform(cs,i)
            elif  caseString=="cart2roof":
                cs=c.roof2cart.invTransform(cs,i)
            
            elif caseString=="cellar2cart":
                cs=c.cellar2cart.transform(cs,i)
            elif caseString=="cart2cellar":
                cs=c.cellar2cart.invTransform(cs,i)
            
            elif caseString=="roof2roof_phys":
                cs=c.roof2roof_phys.transform(cs,i)
            elif caseString=="roof_phys2roof":
                cs=c.roof2roof_phys.invTransform(cs,i)
            
            elif caseString=="cellar2cellar_phys":
                cs=c.cellar2cellar_phys.transform(cs,i)
            elif caseString=="cellar_phys2cellar":
                cs=c.cellar2cellar_phys.invTransform(cs,i)
            
        return(cs)
    
    ##########################################################
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
    ##########################################################
    def transpose(self):
        #switch the component types
        ct=copy.deepcopy(self.componentTypes)
        fct=ct[0 ]
        lct=ct[-1]
        ct[0]=lct
        ct[-1]=fct
        #switch the components 
        c=self.components
        keys=c.keys()
        nc={} 
        for k in keys:
            newkey=exchangeEnds(k)
            nc[newkey]=c[k]
        new=Tensor(self.coords,ct,nc) 
        return(new)
    ##########################################################
    def str(self):
        description="coords="+str(self.coords)\
                +"componentTypes="+str(self.componentTypes)\
                +"components="+str(self.components)
        return(description)
    ##########################################################
    def __eq__(self,other):
        self.purge()
        other.purge()
        boolval=\
        type(self.coords)==type(other.coords) and \
        self.componentTypes==other.componentTypes and\
        self.components==other.components
        return(boolval)
    
    ##########################################################
    def subs(self,*args):
        sc=copy.deepcopy(self)
        c=sc.components
        n={}
        for k in c.keys():
            n[k]=c[k].subs(*args)
        sc.components=n
        return(sc)


    ##########################################################
    def partder(self,i):
        sct=self.componentTypes
        sco=self.coords
        sC=self.components
        rc={}
        for k in sC.keys():
            rc[k]=self.covder(i,(k))
        
        res=Tensor(sco,sct,rc)
        return(res)
    ##########################################################
    def covder(self,i,k):
        # this function computes the covariant derivative of a tensor
        # where i is the index of the partial derivative and k is an index tupel
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
        #if sct==["cart"]:
        #    if not(k in keys):
        #        s=0
        #    else:
        #        s=diff(sC[k],sco.X[i])
        #    return(simplify(s))

        #elif sct==["roof"]:#p.79 eq.4.17
        if sct==["roof"]:#p.79 eq.4.17
            if not(k in keys):
                s=0
            else:
                s=diff(sC[k],sco.U[i])
            for j in keys:
                    s+=sco.Gamma[k[0],i,j[0]]*sC[j]
        elif sct==["cellar"]:#p.79 eq.4.19
            if not(k in keys):
                s=0
            else:
                s=diff(sC[k],sco.U[i])
            for j in keys:
                    s=s-sco.Gamma[j[0],i,k[0]]*sC[j]
        elif sct==["roof","roof"]:#p.98  eq 4.117
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
        else:
            raise NotImplementedError("The covariant derivative is implemented for roof and cellar components only, at the moment for tensors of order 2 at most.")
            
        
        return(simplify(sco.scalarSimp(s)))
            
    ##########################################################
    def nabla(self):
        # This function applies the Nabla operator (always from left) 
        # to the tensor.
        csc=self.coords
        gr=csc.t_gr
        #            i
        # nabla v = g * v,i (* = "direct Product")
        f=lambda i :Tensor(csc,["cellar"],{(i,):1})*self.partder(i)
        res=sum(map(f,range(0,csc.n)))
        return(res)

    ##########################################################
    def grad(self):
        cs=copy.deepcopy(self)
        return((cs.nabla()).transpose())
