#!/usr/bin/python
# vim: set expandtab ts=4
import unittest 
from sympy import bspline_basis
from Spherical import *
import numpy as np
import mpmath as mp

class Cylinder:    
    def __init__(self,etha,l,startset,coord):
        self.etha=etha
        self.l=l
        self.rootlist=self.alphafinder(l,etha,startset)
    self.coords=coord

    def compute_matrix(self,Cl):    
        rootlist=self.rootlist
        l=self.l
        n=len(rootlist)
        ra=range(n)
        N=zeros(n)
        for j in ra:
            alpha=rootlist[j]
            N[j,j]=self.n_l_one_half_j(j)
        Np=N*alpha**2
        Ncl=N/alpha**4
        Q=zeros(n)
        for k in ra:
            for j in ra:
                Q[k,j]=self.Q(k,j)
        #print("n=",n)
        MQ=l*(l+1)*(Ncl+Q)
        M=Np-Cl*MQ
        #m=zeros(n)
        #for k in ra:
        #    for j in ra:
        #        m[k,j]=self.n_l_one_half_j(k)*(rootlist[k]**2/(l*(l+1)*Cl)-1/rootlist[k]**4)*self.kron(k,j)-self.Q(k,j)
        return(M,Np,MQ)

    def computeCl(self):
        Cl=Symbol("Cl")
        rootlist=self.rootlist
        l=self.l
        m,Np,MQ=self.compute_matrix(Cl)
        T=MQ.inv()
        NpT=Np*T
        ev=NpT.eigenvals()
        #dm=m.det()
        #def cl_root(v):
        #    res=dm.subs(Cl,v)
        #    return(res)
        #Clval=solve(dm,Cl)[0]
        #print("Clval=",Clval)
        def f(val):
            return(re(val)>0 and abs(im(val))<1e-16)
        Clval=re(min(filter(f,ev)))
        # print("filtered ev=",filter(f,ev))
        print("Clval=smallest ev=",Clval.evalf())
        #print("Clval=",Clval)
        #print(mp.log(float(Clval),10))
        #Clvaln=[mp.findroot(cl_root,k) for k in [Clval]] 
        #print("Clvalnum=",Clvaln)
        #print(mp.log(Clvaln[0],10))
        #return(5.211*10**3)#faked val
        # another different approach to compute Cl
        #j=0
        #a2=l*(l+1)*Cl-self.n_l_one_half_j(j)*rootlist[j]**6/(self.n_l_one_half_j(j)+self.Q(j,j)*rootlist[j]**4)
        #Clval2=solve(a2,Cl)[0]
        #print("Clval2=",Clval2)
        return(Clval)#faked val

    # Anohter ingredient to the analytic solution is a function that computes
    # the integral property (eq. (138) chandr)
    def n_l_one_half_j(self,j):
        etha=self.etha
        l=self.l
        rootlist=self.rootlist
        res=2/(pi*rootlist[j])**2\
        *(mp.besselj(l+mp.mpf(1)/mp.mpf(2),rootlist[j]*etha)**2/mp.besselj(l+mp.mpf(1)/mp.mpf(2),rootlist[j])**2 -1)
        return(res)
    
    #define a function to compute the Kronecker delta
    def kron(self,x,y):
        if y==x:
            res=1
        else:
            res=0
        return(res)
    
    # define a function to coumpute the Qkj according to (eq. (166)  chandr)
    def Q(self,k,j):
        rootlist=self.rootlist
        l=self.l
        res=(2/self.rootlist[k]**3)*((2*l+3)*self.g_l_one_half(k)*self.B2(j)-(2*l-1)*self.h_l_one_half(k)*self.B4(j))
        return(res)

    def cprime_maker(self,alpha):
        etha=self.etha
        l=self.l
    S=self.coords
    r,phi,theta=S.U
        alpha_cyl=self.cyl_maker(alpha)
        expr=1.0/alpha*diff(alpha_cyl(r),r)
        def f(arg):
            return(expr.subs(r,arg))
        return(f)           

    def C_prime_l_one_half(self,j):
    # define a function to coumpute the C_prime_l_one_half according to
    # (eq. (139))
        alpha=self.rootlist[j]
        l=self.l
        etha=self.etha
        res=-1**l*2.0/(pi*alpha)*mp.besselj(l+mp.mpf(1)/mp.mpf(2),etha*alpha)/mp.besselj(l+mp.mpf(1)/mp.mpf(2),alpha)
        return(res)
    def C_prime_l_one_half_etha(self,j):
        alpha=self.rootlist[j]
        l=self.l
        etha=self.etha
        res=-1**l*2.0/(pi*alpha*etha)
        return(res)
    
    
    # define a function to coumpute the Glj according to (eq. 159)
    def g_l_one_half(self,k):
        l=self.l
        etha=self.etha
        res=self.C_prime_l_one_half(k)-etha**(l+mp.mpf(3)/mp.mpf(2))*self.C_prime_l_one_half_etha(k)
        return(res)
    
    # define a function to coumpute the Glj according to (eq. 160)
    def h_l_one_half(self,k):
        l=self.l
        etha=self.etha
        res=self.C_prime_l_one_half(k)-etha**(-l+mp.mpf(1)/mp.mpf(2))*self.C_prime_l_one_half_etha(k)
        return(res)
    
    # define the parameter B2 for free surfaces (eq. 164)
    def B2(self,j):
        l=self.l
        etha=self.etha
        #res=mp.mpf(1)/(mp.mpf(2*l+1)*(1-etha**(2*l+3))*self.rootlist[j]**3)*self.g_l_one_half(j)
        res=mp.mpf(1)/(mp.mpf(2*l+1)*(1-etha**(2*l+3))*self.rootlist[j]**3)*self.g_l_one_half(j)
        return(res)
    
    # define the parameter B4 for free surfaces (eq. 165)
    def B4(self,j):
        l=self.l
        etha=self.etha
        res=mp.mpf(1)/(mp.mpf(2*l+1)*(etha**(-2*l+1)-1)*self.rootlist[j]**3)*self.h_l_one_half(j)
        return(res)
    
    # costruct the radial part of the solution
    def bs(self,alpha):
        l=self.l
        et=self.etha
        # the matrix contains eq. 162,163,147,148
        bmat=Matrix(4,4,
            [l*(l-1)            ,(l+2)*(l+1)         ,(l+2)*(l+1)            ,l*(l-1),
             l*(l-1)*et**(l-2)  ,(l+2)*(l+1)*et**l   ,(l+2)*(l+1)*et**-(l+3) ,l*(l-1)*et**-(l+1),
             1                  ,1                   ,1                      ,1,
             et**l              ,et**(l+2)           ,et**-(l+1)             ,et**-(l-1)]
        )
        cp=self.cprime_maker(alpha)
        rhs=Matrix(4,1,
                   [2*cp(1)/alpha**3,
                    2*et**(-3.0/2)*cp(et)/alpha**3,
                    0,
                    0]
        )
        system=bmat.row_join(rhs)
        B1,B2,B3,B4=symbols("B1 B2 B3 B4")
        sol=solve_linear_system(system,B1,B2,B3,B4)
        B1=sol[B1]
        B2=sol[B2]
        B3=sol[B3]
        B4=sol[B4]
        # print("Bs",sol)
        # print("sum of Bs=",(B1+B2+B3+B4).evalf())
        return(B1,B2,B3,B4)

    def radialpart(self):
        # We construct the solutions for F an W according to (eq. 143,146)
        #r=Symbol("r")
    S=self.coords
    r,phi,theta=S.U
        rootlist=self.rootlist
        l=self.l
        etha=self.etha
        n=len(rootlist)
        # print("n=",n)
        # print("rootlist",rootlist)
        # We have to know the eigenvectors of the secular matrix 
        # which are determined up to a constanct factor only
        # and therefore an arbitrary constant only in the one dimensional
        # case
        if n==1:
            Ajs=[1]
        else:
            # we have chosen Cl to make the determinant of m zero
            # that is why we expect one eigenvalue to be zero but since our result
            # is obtained numerically we migth have to correct the matrix a bit 
            # to be actually able to compute this eigenspace. To this end we compute
            # the eigenvalues and choose the one  with the smallest absolute value 
            # and call it "sv"
            # we check that sv is nearly zero. Then we substract and sv*I from m
            # and so correct the small deviation from  now compute the nullspace
            # if one eigenvalue is actually 0 (as expected) nothing gets changed
            # sice sympi seams not able to handle this kind of matrix and refuses to
            # compute the nullspce although it computes the eigenvalues we will try
            # numpy instead which has the advantage that te
            Cl=self.computeCl()
            m,Np,MQ=self.compute_matrix(Cl)
            U,s,Vh=np.linalg.svd(m)
            sv=s[-1] # always take the last one witch should be the one with the
            if sv > 1e-9:
               raise Exception('matrix not singular', 'eggs')
            # print("u,s,Vh,sv",U,s,Vh,sv)
            Ajs=U[-1] 
            # print("Ajs=",Ajs)
            am=Matrix(n,1,Ajs)
            # print("am",am)
            # print("m*am",m*am)
            # We can construct a solution as follows:
        F=0
        W=0
        # To construct the analytic solution for W we have to determine the
        # remaining constants B1(j) and B3(j)
        # We can do so by backsubstituting the solutions for B2(j) and B4(j) into
        # (eq. 147) and (eq. 148) 
        # there is a part of the linear system independent of j 
        for j in range(n):
            alpha=rootlist[j]
            alpha_cyl=self.cyl_maker(alpha)
            Fj=Ajs[j]*1/sqrt(r)*alpha_cyl(r)
            F=F+Fj
            B1,B2,B3,B4=self.bs(alpha)
            Bpart=B1*r**l+B2*r**(l+2)+B3*r**-(l+1)+B4*r**-(l-1)
            # testnumrad((self.Dl(self.Dl(Bpart))),r,etha,1,0.01)
            # #test eq.145
            # testnumrad_eq(self.Dl(alpha_cyl(r)/sqrt(r)),-alpha**2*alpha_cyl(r)/sqrt(r),r,etha,1,0.01)
            # #test eq. 145 applied twice
            # testnumrad_eq(self.Dl(self.Dl(alpha_cyl(r)/sqrt(r))),alpha**4*alpha_cyl(r)/sqrt(r),r,etha,1,0.01)
            Wj=Ajs[j]*1.0/(alpha**4*sqrt(r))*alpha_cyl(r)+Bpart
            #testnumrad((self.Dl(self.Dl(Wj))-Fj),r,etha,1)
            W=W+Wj

        #rescale the solution
        F=F*1.0/sqrt(r)
        def f(arg):
            return(F.subs(r,arg))
        def w(arg):
            return(W.subs(r,arg))
        return(f,w)

    def cyl_maker(self,alpha):
        etha=self.etha
        l=self.l
        #r=Symbol("r")
    S=self.coords
    r,phi,theta=S.U
        co=l+Rational(1,2)
        c1=besselj(-co,alpha*etha)
        c2=besselj(+co,alpha*etha)
        exp=(c1*besselj(co,alpha*r)-c2*besselj(-co,alpha*r))
        # note that the construction in eq. 135 is somewhat misleading due to
        # the  argument (alpha*r) We do not stick to this but rather understand
        # f as a function of r which will accordingly vanish for r=etha and r=1
        def f(arg):
            return(exp.subs(r,arg))
        return(f)
    def Dl(self,w):
        l=self.l
        #r= Symbol("r")
    S=self.coords
    r,phi,theta=S.U
        res=diff(w,r,2)+2/r*diff(w,r)+-l*(l+1)/r**2*w
        return(res)
    def alphafinder(self,l,etha,startset):
        def zero(alpha):    
            #res= mp.besselj(-(l+mp.mpf(1)/mp.mpf(2)),etha*alpha)*mp.besselj(l+mp.mpf(1)/mp.mpf(2),alpha) -mp.besselj(l+mp.mpf(1)/mp.mpf(2),etha*alpha) *mp.besselj(-(l+mp.mpf(1)/mp.mpf(2)),alpha)
            res=\
            +besselj(-(l+1.0/2),etha*alpha)*besselj(l+1.0/2,alpha)\
            -besselj(l+1.0/2,etha*alpha)*besselj(-(l+1.0/2),alpha)
            return(re(res))
        # we find the roots (taking approximate startvalues for the
        # root findig routine
        roots=set({mp.findroot(zero,k) for k in startset})
        #print(sorted(roots))
        #we use only the first two approximatons (sizes of the secular matrix)
        rootlist=list(sorted(roots))
        #mp.plot(zero,[10,40])
        #print(rootlist)
        return(rootlist)
    def cyl_one_half_nu_maker(self,l,nu,alpha,etha):
        # r=Symbol("r")
        S=self.coords
        r,phi,theta=S.U
        co=l+Rational(1,2)
        c1=besselj(-co,alpha*etha)
        c2=besselj(+co,alpha*etha)
        exp=(c1*besselj(nu,alpha*r)-c2*besselj(-nu,alpha*r))
        def f(arg):
            return(exp.subs(r,arg))
            
        return(f)

