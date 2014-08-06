
from sympy import symbols, sin, cos, pi, Function
from sympy.diffgeom import Manifold, Patch, CoordSystem
from sympy.diffgeom import BaseVectorField, Differential, TensorProduct, WedgeProduct
from sympy import pprint
r, phi, theta = symbols('r, phi, theta')
m = Manifold('M', 3)
patch = Patch('P', m)
rect = CoordSystem('rect', patch,['x','y','z'])
polar = CoordSystem('polar', patch,['r','phi','theta'])
polar.connect_to(rect, [r, phi, theta], [r*cos(phi)*sin(theta), r*sin(phi)*sin(theta),r*cos(theta)])
# define coordfuntions
m.x=rect.coord_function(0)
m.y=rect.coord_function(1)
m.z=rect.coord_function(2)
m.e_x = BaseVectorField(rect, 0)
m.e_y= BaseVectorField(rect, 1)
m.e_z= BaseVectorField(rect, 2)
m.dx = Differential(m.x)
m.dy = Differential(m.y)
m.dz = Differential(m.z)
omega_x = Function('omega_x')
omega_y = Function('omega_y')
omega_z = Function('omega_z')


m.r=polar.coord_function(0)
m.phi=polar.coord_function(1)
m.theta=polar.coord_function(2)
m.e_r = BaseVectorField(polar, 0)
m.e_phi= BaseVectorField(polar, 1)
m.e_theta = BaseVectorField(polar, 2)
m.dr = Differential(m.r)
m.dphi = Differential(m.phi)
m.dtheta = Differential(m.theta)

omega_r = Function('omega_r')
omega_phi = Function('omega_phi')
omega_theta = Function('omega_theta')

omega=omega_x(m.x,m.y,m.z)*WedgeProduct(m.dy,m.dz)+omega_y(m.x,m.y,m.z)*WedgeProduct(m.dz,m.dx)+omega_z(m.x,m.y,m.z)*WedgeProduct(m.dx,m.dy)
d_omega=Differential(omega)
pprint(d_omega)
