#!/usr/bin/python
# vim:set ff=unix expandtab ts=4 sw=4:

import copy 
from sympy import *
#from helperFunctions import *
#import Exceptions 
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
        res=zeros(s[0],s[1])
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
            res=zeros(s[0],s[1])
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
        g=zeros(self.n,1)
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
        rc=zeros(self.n,1)
        for i in range(0,self.n):
            rc[i]=roofComponents[i]*s[i]
        return(rc)

    def phys2roof(self,physComponents):
        s=self.scale_factors()
        rc=zeros(self.n,1)
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
        s=zeros(self.n,1)
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
#    # compute the divergence of a tensor T given by its roofroof_components Trr
#        # (these are the component with respect to the cellar base vectors since 
#    # T=Trr[i,j] gc[i] gc[j])
#    # (The result are the roof components of a vector 
#    # (with resprect to the cellar basis))
#    res=zeros(3,1)
#    for k in range(0,self.n):
#        rv=0 # the value of the row of the result
#        for i in range(0,self.n):
#            rv+=self.cov_der_T(i,i,k,Trr)
#        rv=simplify(rv)
#        res[k]=rv    
#    return(res)    
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
