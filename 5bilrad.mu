// The task is to provide the projection that maps any arbitrarily chosen gridcell to the reference cell
read("functions.mu");
read("/home/mm/Markus/Aufgaben/texhelper/texlib.mu");
cellplot:=proc(B1,B2,B3,T1,T2,T3,c)
begin
lw:=1:
list:=[
   //the bottom
   G_ab(B1,B2,c,lw)
   ,
   G_ab(B2,B3,c,lw)
   ,
   G_ab(B3,B1,c,lw)
   ,
   //the roof
   G_ab(T1,T2,c,lw)
   ,
   G_ab(T2,T3,c,lw)
   ,
   G_ab(T3,T1,c,lw)
   ,
   //the connection lines
   plot::Polygon3d([B1,T1],LineWidth=lw,Color=c)
   ,
   plot::Polygon3d([B2,T2],LineWidth=lw,Color=c)
   ,
   plot::Polygon3d([B3,T3],LineWidth=lw,Color=c)
]:
return(list):
end_proc:
proj_plot:=proc(ma,mb,P)
local T,Range;
begin
delete t;
   T:=G(ma,mb)[1]:
   Range:=G(ma,mb)[2]:
   C:=plot::Curve3d(P(T(t)),
     t=Range[1]..Range[2],
     Color=RGB::Red,
     LineStyle=Dashed
   ):
   return(C):
end_proc:
rot:=proc(B1,B2,B3,T1,T2,T3)
   // The task is to provide the projection that maps any 
   // arbitrarily chosen gridcell to the reference cell
   local XD,X,XDI,a,b,c,pim,im,spin_axis,spin_angle,Spin_abc,Spin_xyz,InvSpin,Trans,Trans_xyz;
   save x,y,z;
   option escape;
   begin
   /////////////////////////////////////////////////////////////////////////////////////
   // begin rotation
   /////////////////////////////////////////////////////////////////////////////////////
   // The first step towards the reference cell is a rotation of the cell 
   // to a plane parallel to z=0
   // To compute the rotation we define the preimage and the image of the desired rotation
   // The direction of the image vector is the z-axis its
   // length is the same as the length of the preimage
   // The preimage is the vector from the center of the 
   // sphere to the plane spanned by the 3 bottom points of 
   // the spherical triangle orthogonal to that plane
   // first we compute its direction
   a:=B2-B1:
   b:=B3-B1:
   c:=-linalg::crossProduct(a,b): 
   c:=c/norm(c,2):
   print(norm(c,2)):
   //now the preimage is given by
   pim:=c:
   //and the desired image by
   im:=matrix([0,0,1]):
   //now we can construct the rotation 
   //first we compute the spin axis
   //it must be orthogonal to the the image and preimage
   spin_axis:=linalg::crossProduct(pim,im):
   //next comes the spin angle
   spin_angle:=arccos(linalg::scalarProduct(im,pim)):
   // now we can construct a new orthogonal basis 
   //a,b,c from spin_axis and z_new
   c:=spin_axis:
   a:=im:
   b:=linalg::crossProduct(c,a):
   a:=a/norm(a,2):
   b:=b/norm(b,2):
   c:=c/norm(c,2):
   //now we can construct a spin matrix with respect to the 
// basis a,b,c
Spin_abc :=matrix(
  [
    [ cos(spin_angle),-sin(spin_angle),0],
    [ sin(spin_angle), cos(spin_angle),0],
    [              0 ,             0 ,1]
  ]
 ): 
// This spin matrix must be transformed to the old basis x y z  
Trans:=matrix(3,3):
Trans:=linalg::substitute(Trans,a,1,1):
Trans:=linalg::substitute(Trans,b,1,2):
Trans:=linalg::substitute(Trans,c,1,3):

Spin_xyz:=Trans*Spin_abc*linalg::inverseLU(Trans):
delete x,y,z;
X:=matrix(3,1,[x,y,z]):
XD:=Spin_xyz*X:
//we construct a function that maps the vector Y to the vector XD
Spinfunk:=(Y)-->subs(
   XD,
   x=Y[1],
   y=Y[2],
   z=Y[3]
):
InvSpin:=linalg::inverseLU(Spin_xyz):
XDI:=Spin_xyz*X:
InvSpinfunk:=(Y)-->subs(
   XDI,
   x=Y[1],
   y=Y[2],
   z=Y[3]
):
return([Spinfunk,InvSpinfunk]):
end_proc:

proj:=proc()
option escape;
local XI,X,Xb,Intersection,Projection;
save x,y,z;
////////////////////////////////////////////////////////////////////////////////////////////
//begin linearization
//////////////////////////////////////////////////////////////////////////////////////////
//The next step is to map the spherical surfaces to planes 
//so that every point in the spherical prism is mapped to a point in an linear one.
//This is done by a projection that can be described by a ray starting from the center
//of the sphere going throug every point of the sperical surface attached to every r level
//of the original cell and finally hitting the plane z=1
begin
delete x,y,z;
XI:=matrix(
   [
      x/z,
      y/z,
      1
   ]
):
Intersection:=Y-->subs(
  XI,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
X:=XI:
X[3]:=sqrt(x^2+y^2+z^2):

//from curved to straight 
Projection:=Y-->subs(
  X,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
 
X_b:=matrix(
  [
    x/sqrt(x^2+y^2+1)*z,
    y/sqrt(x^2+y^2+1)*z,
    1/sqrt(x^2+y^2+1)*z
  ]
):
//from straight to curved
Back_Projection:=Y-->subs(
  X_b,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
return([Back_Projection,Projection]):
end_proc:

linkoeffs:=proc(Tri,Tri2)
option escape;
local B,k,Bz,kz,a,b,c,d,e,f,g,h,funk;
save X,x,y,z;
// The last step mappes the Prism with a linear Projection to the reference Prism
// Since the triangles at the bottom and the top of the prism are identical
// it suffices to find the transformation of one of them.
// The affine Transformation of a reference triangle T to a mesh triangle T' is described
// by the following Matrix equation 
//  X'=M X + K  
// Where M is a matrix and the rest are 2D vectors
// We arrange the variables in the following order
// (x') =(a b) (x) +(e)
// (y')  (c d) (y) +(f)
// Since we know the coordinates x',y' of the triange as well
// as the coordinates x, y for the of the reference triangle 
// (for all three vertices)  we can compute the unknowns which we rearrange in the 
// following order.
// (x1  y1   0   0  1   0)  (a) = (x1')
// (x2  y2   0   0  1   0)  (b) = (x2')
// (x3  y3   0   0  1   0)  (c) = (x3')
// (0   0    x1  y1 0   1)  (d) = (y1')
// (0   0    x2  y2 0   1)  (e) = (y2')
// (0   0    x3  y3 0   1)  (f) = (y3')
// The affine mapping of the z Koordinates is described by 
// z'=g*z+h
// we can integrate it in the former system
// (x') =(a b 0) (x) +(e)
// (y')  (c d 0) (y) +(f)
// (z')  (0 0 g) (z) +(h)
// 
// (z1 1) (g) = (z1') 
// (z2 1) (h) = (z2') 

begin
x1_s:=Tri[1][1]:
x2_s:=Tri[2][1]:
x3_s:=Tri[3][1]:
y1_s:=Tri[1][2]:
y2_s:=Tri[2][2]:
y3_s:=Tri[3][2]:
z1_s:=Tri[1][3]:

z1_t:=Tri2[1][3]:

TriRef:=refTri(0):
x1:=TriRef[1][1]:
x2:=TriRef[2][1]:
x3:=TriRef[3][1]:
y1:=TriRef[1][2]:
y2:=TriRef[2][2]:
y3:=TriRef[3][2]:
B:=matrix(6,6,[
     [x1, y1,  0,  0, 1, 0],
     [x2, y2,  0,  0, 1, 0],
     [x3, y3,  0,  0, 1, 0],
     [ 0,  0, x1, y1, 0, 1],
     [ 0,  0, x2, y2, 0, 1],
     [ 0,  0, x3, y3, 0, 1]
]):

r:=matrix(6,1,[x1_s,x2_s,x3_s,y1_s,y2_s,y3_s]):
k:=linalg::matlinsolve(B,r):
[a,b,c,d,e,f]:=[op(k)]:
z1:=0:
z2:=1:
Bz:=matrix(2,2,[
   [z1,1],
   [z2,1]
 ]):
rz:=matrix(2,1,[z1_s,z1_t]):
kz:=linalg::matlinsolve(Bz,rz):
[g,h]:=[op(kz)]:

M:=matrix(3,3,[
   [a,b,0], 
   [c,d,0],
   [0,0,g]
]):
K:=matrix(3,1,[e,f,h]):
delete x,y,z:
X:=matrix(3,1,[x,y,z]):  
XS:=M*X+K: //-> M^-1(XS-K)=X
print("XS=",XS):
//from reference to general
funk:=Y->subs(
  XS,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
//from general to reference
Minv:=linalg::inverseLU(M):
YS:=Minv*(X-K):
inv_funk:=Y->subs(
  YS,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
print("f(X)=",funk(X)):
return([funk,inv_funk]): 
end_proc:
refcell:=proc()
local U,B,RefTri;
begin
U:=refTri(1):
B:=refTri(0):
RefTri:=append(U,op(B)):
print(RefTri):
return(RefTri):
end_proc:

refTri:=proc(h)
local TriRef;
begin
TriRef:=[
   matrix(1,3,[0,0,h]),
   matrix(1,3,[1,0,h]),
   matrix(1,3,[0,1,h])
]:
return(TriRef):
end_proc:   


did:=1://diamondid
//n:=9:
n:=5:
//n:=33:
ri:=0.7:
ro:=0.5:
xi:=grdgen(n,[did],ri):
xo:=grdgen(n,[did],ro):
// select a gridcell 
i1:=1:
i2:=n-1:
print(xi[i1  ,i2,did]):
B1:=matrix(xi[i1  ,i2,did]):
B2:=matrix(xi[i1+1,i2,did]):
B3:=matrix(xi[i1+1,i2+1,did]):
T1:=matrix(xo[i1  ,i2,did]):
T2:=matrix(xo[i1+1,i2,did]):
T3:=matrix(xo[i1+1,i2+1,did]):
// draw it
oc:=RGB::Red:
orgcell1:=cellplot(B1,B2,B3,T1,T2,T3,oc):
//plot(op(orgcell),Scaling=Constrained):

// compute the Mapping that maps the grid cell to the reference cell
delete x,y,z;
X:=matrix(3,1,[x,y,z]):
[finv,f]:=rot(B1,B2,B3,T1,T2,T3):
//print("f(X)=",f(X)):
//now we are able to transform the coordinates of all the points 
//to the new system
B1_rot:=f(B1):
B2_rot:=f(B2):
B3_rot:=f(B3):
T1_rot:=f(T1):
T2_rot:=f(T2):
T3_rot:=f(T3):
//print(B1,B1_rot):
//now we can plot the rotatet cell
rc:=RGB::Blue:
rotcell:=cellplot(B1_rot,B2_rot,B3_rot,T1_rot,T2_rot,T3_rot,rc):
//plot::setDefault(plot::Box::LinesVisible = FALSE):
plot(op(orgcell1),op(rotcell),Axes=Origin,Scaling=Constrained):

//we project the cell to a prism 
[ginv,g]:=proj():
//print("g(f(X))=",g(f(X))):
//and draw the projected cell
cP_proj:=[]:
C1:=proj_plot(B1_rot,B2_rot,g):
C2:=proj_plot(B2_rot,B3_rot,g):
C3:=proj_plot(B3_rot,B1_rot,g):
D1:=proj_plot(T1_rot,T2_rot,g):
D2:=proj_plot(T2_rot,T3_rot,g):
D3:=proj_plot(T3_rot,T1_rot,g):
cP_proj:=[C1,C2,C3,D1,D2,D3]:
//plot(op(cP_proj),Scaling=Constrained):

// the last remaining step is to find the linear transformation that 
// projects to the reference cell
B1_p:=g(f(B1)):
B2_p:=g(f(B2)):
B3_p:=g(f(B3)):
T1_p:=g(f(T1)):
T2_p:=g(f(T2)):
T3_p:=g(f(T3)):
B:=[B1_p,B2_p,B3_p]:
T:=[T1_p,T2_p,T3_p]:
delete x,y,z;
X:=matrix(3,1,[x,y,z]):
[hinv,h]:=linkoeffs(T,B):
//print("X=",X);
//print("f(X)=",f(X)):
//print("g(f(X))=",g(f(X))):
print("F(X)=",h(g(f(X)))):
//Now test if the test grid cell really is mapped to the reference prism
B1_r:=h(g(f(B1))):
B2_r:=h(g(f(B2))):
B3_r:=h(g(f(B3))):
T1_r:=h(g(f(T1))):
T2_r:=h(g(f(T2))):
T3_r:=h(g(f(T3))):
print( B1_r, B2_r, B3_r, T1_r, T2_r, T3_r):
//now we start from the reference cell and project it to a real grid cell
rfc:=refcell():
[t1,t2,t3,b1,b2,b3]:=rfc:
XS:=finv(ginv(hinv(X))):
m:=Y->subs(
  XS,
  x=Y[1],
  y=Y[2],
  z=Y[3]
):
print("back=",map(rfc,m)):
print("original=",B1,B2,B3,T1,T2,T3):
dJ:=linalg::det(linalg::jacobian(XS,[x,y,z])):
XGinv:=ginv(X);
TeX_writeVariable("XGinv");
JGinv:=linalg::jacobian(XGinv,[x,y,z]):
TeX_writeVariable("JGinv");
assume(x,Type::Real):
assume(y,Type::Real):
assume(z,Type::Real):
nJGinv:=simplify((norm(JGinv,Frobenius)^2)):
nJGinv:=simplify(
   subs(
      nJGinv,
      abs(x)^2=x^2,
      abs(y)^2=y^2,
      abs(z)^2=z^2,
      abs(x*y*z)^2=(x*y*z)^2,
      abs(y*z)^2=(y*z)^2,
      abs(x*z)^2=(x*z)^2
   )
):
TeX_writeVariable("nJGinv");
//dJGinv:=linalg::det(linalg::jacobian(XGinv,[x,y,z]));
//gjg:=linalg::grad(JGinv,[x,y,z]);

XG:=g(X);

//JS:=linalg::jacobian(XS,[x,y,z]);
//ev:=linalg::eigenvalues(JS);




//c1:=l_plot(b1,b2,m):
//c2:=l_plot(b2,b3,m):
//c3:=l_plot(b3,b1,m):
//d1:=l_plot(t1,t2,m):
//d2:=l_plot(t2,t3,m):
//d3:=l_plot(t3,t1,m):
cP_proj:=[c1,c2,c3,d1,d2,d3]:
//plot(op(cP_proj),Scaling=Constrained):
