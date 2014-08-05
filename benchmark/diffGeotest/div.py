from sympy import symbols, sin, cos, pi
from sympy.diffgeom import Manifold, Patch, CoordSystem
r, theta = symbols('r, theta')
m = Manifold('M', 2)
patch = Patch('P', m)
rect = CoordSystem('rect', patch)
polar = CoordSystem('polar', patch)
rect in patch.coord_systems
polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])
polar.coord_tuple_transform_to(rect, [0, 2])
polar.coord_tuple_transform_to(rect, [2, pi/2])
rect.coord_tuple_transform_to(polar, [1, 1])
polar.jacobian(rect, [r, theta])
p = polar.point([1, 3*pi/4])
rect.point_to_coords(p)
#Define a basis scalar field (i.e. a coordinate function), that takes a point and returns its coordinates. It is an instance of BaseScalarField.
rect.coord_function(0)(p)
rect.coord_function(1)(p)
#Define a basis vector field (i.e. a unit vector field along the coordinate line). Vectors are also differential operators on scalar fields. It is an instance of BaseVectorField.
v_x = rect.base_vector(0)
x = rect.coord_function(0)
v_x(x)
v_x(v_x(x))
#Define a basis oneform field:
dx = rect.base_oneform(0)
dx(v_x)
#If you provide a list of names the fields will print nicely: - without provided names:
x, v_x, dx
#with provided names
rect = CoordSystem('rect', patch, ['x', 'y'])
rect.coord_function(0), rect.base_vector(0), rect.base_oneform(0)
# Examples
from sympy import symbols, sin, cos, pi
from sympy.diffgeom import (
       Manifold, Patch, CoordSystem, Point)
r, theta = symbols('r, theta')
m = Manifold('M', 2)
p = Patch('P', m)
rect = CoordSystem('rect', p)
polar = CoordSystem('polar', p)
polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])
# Define a point using coordinates from one of the coordinate systems:
p = Point(polar, [r, 3*pi/4])
p.coords()
p.coords(rect)
# Use the predefined R2 manifold, setup some boilerplate.
from sympy import symbols, pi, Function
from sympy.diffgeom.rn import R2, R2_p, R2_r
from sympy.diffgeom import BaseVectorField
from sympy import pprint
x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0')
# Points to be used as arguments for the field:
point_p = R2_p.point([r0, theta0])
point_r = R2_r.point([x0, y0])
#Scalar field to operate on:
g = Function('g')
s_field = g(R2.x, R2.y)
s_field.rcall(point_r)
s_field.rcall(point_p)
#Vector field:
v = BaseVectorField(R2_r, 1)
pprint(v(s_field))
pprint(v(s_field).rcall(point_r).doit())
pprint(v(s_field).rcall(point_p).doit())
# Differential Examples
# scalar 0 forms
from sympy import Function
from sympy.diffgeom.rn import R2
from sympy.diffgeom import Differential
from sympy import pprint
g = Function('g')
s_field = g(R2.x, R2.y)
#Vector fields:
e_x, e_y, = R2.e_x, R2.e_y
#Differentials:
dg = Differential(s_field)
dg
pprint(dg(e_x))
pprint(dg(e_y))
# TensorProduct Examples
from sympy import Function
from sympy.diffgeom.rn import R2
from sympy.diffgeom import TensorProduct
from sympy import pprint
TensorProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y)
TensorProduct(R2.dx, R2.dy)(R2.e_y, R2.e_x)
TensorProduct(R2.dx, R2.x*R2.dy)(R2.x*R2.e_x, R2.e_y)
#You can nest tensor products.
tp1 = TensorProduct(R2.dx, R2.dy)
TensorProduct(tp1, R2.dx)(R2.e_x, R2.e_y, R2.e_x)
#You can make partial contraction for instance when ‘raising an index’. Putting None in the second argument of rcall means that the respective position in the tensor product is left as it is.
TP = TensorProduct
metric = TP(R2.dx, R2.dx) + 3*TP(R2.dy, R2.dy)
metric.rcall(R2.e_y, None)
# WedgeProduct Examples
from sympy import Function
from sympy.diffgeom.rn import R2
from sympy.diffgeom import WedgeProduct
from sympy import pprint

WedgeProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y)
WedgeProduct(R2.dx, R2.dy)(R2.e_y, R2.e_x)
WedgeProduct(R2.dx, R2.x*R2.dy)(R2.x*R2.e_x, R2.e_y)
#You can nest wedge products
wp1 = WedgeProduct(R2.dx, R2.dy)
WedgeProduct(wp1, R2.dx)(R2.e_x, R2.e_y, R2.e_x)
# LieDerivative Examples
from sympy.diffgeom import (LieDerivative, TensorProduct)
from sympy.diffgeom.rn import R2
LieDerivative(R2.e_x, R2.y)
LieDerivative(R2.e_x, R2.x)
LieDerivative(R2.e_x, R2.e_x)
# The Lie derivative of a tensor field by another tensor field is equal to their commutator:
LieDerivative(R2.e_x, R2.e_r)
LieDerivative(R2.e_x + R2.e_y, R2.x)
tp = TensorProduct(R2.dx, R2.dy)
LieDerivative(R2.e_x, tp)
LieDerivative(R2.e_x, tp).doit()
#BaseCovarDerivativeOp
from sympy.diffgeom.rn import R2, R2_r
from sympy.diffgeom import BaseCovarDerivativeOp
from sympy.diffgeom import metric_to_Christoffel_2nd, TensorProduct
TP = TensorProduct
ch = metric_to_Christoffel_2nd(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
ch
cvd = BaseCovarDerivativeOp(R2_r, 0, ch)
cvd(R2.x)
cvd(R2.x*R2.e_x)
#CovarDerivativeOp
from sympy.diffgeom.rn import R2
from sympy.diffgeom import CovarDerivativeOp
from sympy.diffgeom import metric_to_Christoffel_2nd, TensorProduct
TP = TensorProduct
ch = metric_to_Christoffel_2nd(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
ch
cvd = CovarDerivativeOp(R2.x*R2.e_x, ch)
cvd(R2.x)
cvd(R2.x*R2.e_x)
#intcurve_series
from sympy.abc import t, x, y
from sympy.diffgeom.rn import R2, R2_p, R2_r
from sympy.diffgeom import intcurve_series
# Specify a starting point and a vector field:
start_point = R2_r.point([x, y])
vector_field = R2_r.e_x
# Calculate the series:
intcurve_series(vector_field, t, start_point, n=3)
# Or get the elements of the expansion in a list:
series = intcurve_series(vector_field, t, start_point, n=3, coeffs=True)
series[0]
series[1]
series[2]
#intcurve_diffequ
from sympy.abc import t
from sympy.diffgeom.rn import R2, R2_p, R2_r
from sympy.diffgeom import intcurve_diffequ
#Specify a starting point and a vector field:
start_point = R2_r.point([0, 1])
vector_field = -R2.y*R2.e_x + R2.x*R2.e_y
# get the equations, 
equations, init_cond = intcurve_diffequ(vector_field, t, start_point)
equations
init_cond
# The series in the polar coordinate system:
equations, init_cond = intcurve_diffequ(vector_field, t, start_point, R2_p)
equations
init_cond
#vectors_in_basis
#Transform all base vectors in base vectors of a specified coord basis.
#While the new base vectors are in the new coordinate system basis, any coefficients are kept in the old system.
from sympy.diffgeom import vectors_in_basis
from sympy.diffgeom.rn import R2_r, R2_p
vectors_in_basis(R2_r.e_x, R2_p)
vectors_in_basis(R2_p.e_r, R2_r)
#twoform_to_matrix
from sympy.diffgeom.rn import R2
from sympy.diffgeom import twoform_to_matrix, TensorProduct
TP = TensorProduct
twoform_to_matrix(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
twoform_to_matrix(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
twoform_to_matrix(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy) - TP(R2.dx, R2.dy)/2)
# metric_to_Christoffel_1st(expr) Return the nested list of Christoffel symbols for the given metric.  This returns the Christoffel symbol of first kind that represents the Levi-Civita connection for the given metric.
from sympy.diffgeom.rn import R2
from sympy.diffgeom import metric_to_Christoffel_1st, TensorProduct
TP = TensorProduct
metric_to_Christoffel_1st(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
metric_to_Christoffel_1st(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
# metric_to_Riemann_components(expr) Return the components of the Riemann tensor expressed in a given basis.  Given a metric it calculates the components of the Riemann tensor in the canonical basis of the coordinate system in which the metric expression is given.
from sympy import pprint, exp
from sympy.diffgeom.rn import R2
from sympy.diffgeom import metric_to_Riemann_components, TensorProduct
TP = TensorProduct
metric_to_Riemann_components(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))

non_trivial_metric = exp(2*R2.r)*TP(R2.dr, R2.dr) +         R2.r**2*TP(R2.dtheta, R2.dtheta)
non_trivial_metric
riemann = metric_to_Riemann_components(non_trivial_metric)
riemann[0]
riemann[1]
#metric_to_Ricci_components Return the components of the Ricci tensor expressed in a given basis.  Given a metric it calculates the components of the Ricci tensor in the canonical basis of the coordinate system in which the metric expression is given.
from sympy import pprint, exp
from sympy.diffgeom.rn import R2
from sympy.diffgeom import metric_to_Ricci_components, TensorProduct
TP = TensorProduct
metric_to_Ricci_components(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
non_trivial_metric = exp(2*R2.r)*TP(R2.dr, R2.dr) +                              R2.r**2*TP(R2.dtheta, R2.dtheta)
non_trivial_metric
metric_to_Ricci_components(non_trivial_metric) #TODO why is this not simpler

