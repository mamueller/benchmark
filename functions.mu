scale:=proc(
l ,//:Type::ListOf(Type::Real),
sf:Type::Real)
save l;
local n,k;
begin
n:=nops(l);
for k from 1 to n do
  l[k]:=l[k]*sf:
end_for:
return(l);
end_proc:

xnextend:=proc(
	      n			:Type::PosInt,
	      scalefactor=1	:Type::Real)
	      
  //This function creates an index shifted copy of xntest whith additional ghostno    des
  //while grdgen creates an table with indices from 1 to nt+1 for both i1 and 12 in xnex the idices have values from -1 to nt+1 for i1 and 0 to nt+2 for i2
local   
	xn_org :TableOfEntry(Type::Real)
	,xntest	
	,l
	,idlr
	,idur
	,idul
	,idll
	,id3r
	,did
        ,i1
	,i2
  ;
  begin
    l:=[1,2,3,4,5,6,7,8,9,10]:
    xn_org:=grdgen(n,l,scalefactor):
    //loading the none ghost node coordinates  with shifted i1 in a copy
    //to achieve the same indeces as used in terra
 
      for did from 1 to 10 do
	for i1 from 0 to n do:
	  for i2 from 1 to n+1 do:
	    xntest[i1,i2,did]:=xn_org[i1+1,i2,did]:
	  end_for:
	end_for:
      end_for:	
    //loading ghost node coordinates
    
      for did from 1 to 10 do
         idlr := did + 5:
         idur := modp(did  ,5) + 1:
         idul := modp(did+3,5) + 1:
         idll := idul + 5:
         id3r := modp(did+2,5) + 1:
 
         if(did >= 6) then
            idlr := idur:
            idur := idur + 5:
            idul := idul + 5:
            idll := did   - 5:
            id3r := id3r + 5:
         end_if:
 
            for i2 from 1 to n+1 do:
               for i1 from 0 to n do :
                  xnex[i1,i2,did] := xntest[i1,i2,did]:
               end_for:
            end_for:
	    
            for i from 1 to n+1 do
               xnex[i,0,did] := xntest[1 ,i,idul]:
               xnex[-1,i+1,did] := xntest[i-1,2,idur]:
               xnex[n+1,i,did] := xntest[n+1-i,n ,idll]:
               xnex[i-1,n+2 ,did] := xntest[n-1,n+2-i ,idlr]:
            end_for:
 
            xnex[  0,  0,did] :=  xntest[1,  1,id3r]:
            xnex[ -1,  0,did] :=  xntest[1,  1,id3r]:
            xnex[ -1,  1,did] :=  xntest[1,  1,id3r]:
            xnex[n+1,n+1,did] := xnex[n,n+2,did ]:
            xnex[n+1,n+2,did] := xnex[n,n+2,did ]:
 
       end_for:
 
    return(xnex):
  end_proc;

G:=proc(
  a:Dom::Matrix,
  b:Dom::Matrix
)
  option escape;
  save 
  Kart,
  c,
  bneu,
  T,
  t,
  w_max,
  as,
  bs,
  K,
  w,
  r,
  ra,
  rb,
  Range;
  begin
    //Überprüfung der Eingangsdaten auf gleiche Radien
    ra:=float(norm(a,2)):
    rb:=float(norm(b,2)):
    if abs((ra-rb)) >= abs(0.01*ra) then
      error("Two points of different radii can not be connected by a part
      of a great circle")
    end_if:
    r:=ra:
    //Bestimmung der Drehachse
    c:=linalg::crossProduct(a,b):
     
    
    //Normierung auf Einheitslänge
    a:=a/norm(a,2):
    b:=b/norm(b,2):
    c:=c/norm(c,2):
    
    //Orthogonalisierung der Vektoren a und b
    
    bneu:=linalg::crossProduct(c,a):
    
    //Damit haben wir ein Rechtssystem für den gedrehten Großkreis
    //und damit die Transformation eines Referenzrechtssystems
    T:=matrix(3,3):
    T:=linalg::setCol(T,1,a):
    T:=linalg::setCol(T,2,bneu):
    T:=linalg::setCol(T,3,c):
    Kart:=t-->T*matrix([r*cos(t),r*sin(t),0]):
    
    w_max:=float(arccos(linalg::scalarProduct(a,b))):
    //w_max:=2*PI:
    Range:=[0,w_max]:
    //Kart ist eine Function (Closure) die dem Parameter t die 
    //kartesischen Kurvenpunkte zuordne:t
    return([Kart,Range]);
end_proc:

G_ab:=proc(
  a:Dom::Matrix,
  b:Dom::Matrix,
  color)
save G,T,R,t,C;
begin
  T:=G(a,b)[1]:
  R:=G(a,b)[2]:
  C:=plot::Curve3d(T(t),
    t=R[1]..R[2],
    Color=color,
    LineWidth=2,
    LineStyle=Dashed):
    return(C):
end_proc:

grdgen:=proc(
	    n		:Type::PosInt,
	    list	:Type::ListOf(Type::PosInt),
	    scalefactor=1:Type::Real) 
	    local d;
    save    cosw,  did, fifthpi, i1, i2, j1, j2, k, l, l2, lvt, m, 
	    midpt, nt, phi, scale, sgn, sinw, w, xn;
    begin
    lvt     := floor(1.45*ln(n)):
    nt:=2^lvt:
    //xn:=TableOfEntry(Type::Real);
    fifthpi:=float(PI/5):
    w       := 2.0*arccos(1.0/(2.0*sin(fifthpi))):
    cosw    := cos(w):
    sinw    := sin(w):
    for d from 1 to nops(list) do
      did:=list[d];
      if did >=6.0
        then sgn := -1.0:
	else sgn := 1.0:
      end_if:
      phi := (2*(did mod 5) - 3 + floor((did - 1)/5))*fifthpi:
      xn[   1,   1,did] :=[
    			  0.0,
    			  0.0,
    			  sgn
    		      ]:
      xn[nt+1,   1,did] :=[
    			  sinw*cos(phi),
    			  sinw*sin(phi),
    			  cosw*sgn
    		      ]:
      xn[   1,nt+1,did] :=[
    			sinw*cos(phi + fifthpi + fifthpi),
    			sinw*sin(phi + fifthpi + fifthpi),
    			cosw*sgn
    			]:
      xn[nt+1,nt+1,did] :=[
    			sinw*cos(phi + fifthpi),
    			sinw*sin(phi + fifthpi),
    			-cosw*sgn
    		      ]:
    
      for k from 0 to lvt-1 do
    
         m  := 2^k:
         l  := nt/m:
         l2 := l/2:
     
    //rows of diamond--
    
        for j1 from 1 to m+1 do
          for j2 from 1 to m do
            i1 := (j1-1)*l + 1:
            i2 := (j2-1)*l + l2 + 1:
            xn[i1,i2,did]:=midpt(xn[i1,i2-l2,did],xn[i1,i2+l2,did]):
          end_for:
       end_for:
       
    //columns of diamond--:
        for j1 from 1 to m+1 do
          for j2 from 1 to m do
    	i1 := (j2-1)*l + l2 + 1:
    	i2 := (j1-1)*l + 1:
    	xn[i1,i2,did]:=midpt(xn[i1-l2,i2,did], xn[i1+l2,i2,did]):
          end_for :
        end_for:
     
    //diagonals of diamond--:
    
        for j1 from 1 to m do
          for j2 from 1 to m do
    	  i1 := (j1-1)*l + l2 + 1:
    	  i2 := (j2-1)*l + l2 + 1:
    	  xn[i1,i2,did]:=midpt(xn[i1-l2,i2+l2,did],xn[i1+l2,i2-l2,did]):
          end_for:
        end_for:
      end_for:
    end_for:

    for d from 1 to nops(list)do 
      did:=list[d];
        for i1 from 1 to nt+1 do 
            for i2 from 1 to nt+1 do 
	      xn[i1,i2,did]:=scale(xn[i1,i2,did],scalefactor):
	    end_for
	end_for
    end_for;
    return(xn):
end_proc:

midpt:=proc(
    x1:Type::ListOf(Type::Real),
    x2:Type::ListOf(Type::Real))
    local xnorm;
begin
    x:=[0.0,0.0,0.0]:
    for i from 1 to 3 do
      x[i]:=x1[i]+x2[i]:
    end_for:
    //xnorm:=float(norm(matrix(3,1,x),2)):
    xnorm:=sqrt(_plus(x[i]^2$i=1..3)):
    for i from 1 to 3 do
      x[i]:=x[i]/xnorm:
    end_for:
    //print(norm(matrix(3,1,x),2));
    return(x):
end_proc:






dp:=proc(
	n		:Type::PosInt,
	list		:Type::ListOf(Type::PosInt),
        colors		:Type::ListOf(Type::ListOf(Type::NonNegative)),
	scalefactor=1	:Type::Real
    )
  local lvt,idmax,nt,p,pol1,pol2,v1,v2,v3,v4,c,xn;
  begin
    lvt:=floor(log(2,n)):
    nt:=2^lvt;
    xn:=grdgen(n,list,scalefactor):
    idmax:=nops(list):
    p:=[]:
    a:=Solid:
    for i from 1 to idmax do 
        did:=list[i]:
        c:=colors[(i-1 mod nops(colors))+1]:
        for i1 from 1 to nt do 
            for i2 from 1 to nt do 
    // 	create the verticesvector for 1 triangle
                v1:=xn[i1   ,i2    ,did]:         
                v2:=xn[i1   ,i2+1  ,did]:
                v3:=xn[i1+1 ,i2    ,did]:
                v4:=xn[i1+1 ,i2+1  ,did]:
		print [v4,v2,v3]:
    //      Create  Polygons.
	       	pol1:=plot::Polygon3d([v1,v2,v3],
		  Color=c,
		  Closed=TRUE,
		  Filled=TRUE,
		  LineWidth=10,
		  LineStyle=a
		):
    	    
		pol2:=plot::Polygon3d([v4,v2,v3],
		  Color=c,
		  Closed=TRUE,
		  Filled=TRUE,
		  LineWidth=10,
		  LineStyle=a
		):
                p:=append(p,pol1,pol2):
	  end_for
	end_for
    end_for;

    
  return(p);
end_proc;



pointplot:=proc(
		n		:Type::PosInt,
		list		:Type::ListOf(Type::PosInt)
		)
  local lvt,idmax,nt,p,pol1,pol2,v1,v2,v3,v4;
  begin
    lvt:=floor(log(2,n)):
    nt:=2^lvt;
    xn:=grdgen(n,list):
    idmax:=nops(list):
    p:=[]:
    c:=[1,0,0];
    ps:=FilledCircles:
    pw:=60: 
    for i from 1 to idmax do 
        did:=list[i]:
        for i1 from 1 to nt do 
            for i2 from 1 to nt do 
    // 	create the verticesvector for 1 triangle
                v1:=plot::Point3d(xn[i1   ,i2    ,did],
				Color=c,
				PointStyle=ps,
				PointWidth=pw
				):         
                v2:=plot::Point3d(xn[i1   ,i2+1  ,did],
				Color=c,
				PointStyle=ps,
				PointWidth=pw
				):
                v3:=plot::Point3d(xn[i1+1 ,i2    ,did],
				Color=c,
				PointStyle=ps,
				PointWidth=pw
				):
                v4:=plot::Point3d(xn[i1+1 ,i2+1  ,did],
				Color=c,
				PointStyle=ps,
				PointWidth=pw
				):
		p:=append(p,v1,v2,v3,v4):
	  end_for
       	end_for
      end_for;
    return(p);
end_proc;

