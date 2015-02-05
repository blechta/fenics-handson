#Author: Jaroslav Hron <jaroslav.hron@mff.cuni.cz>

import os
from dolfin import *
import math
import mshr
# get file name
fileName = os.path.splitext(__file__)[0]

# Domain
Center = Point(0.2, 0.2)
Radius = 0.05
L=2.2
W=0.41
geometry = mshr.Rectangle(Point(0.0,0.0), Point(L,W))-mshr.Circle(Center,Radius,10)
info(geometry)

# Mesh
begin('Generating mesh...')
mesh=mshr.generate_mesh(geometry,50)
info(mesh)
end()

bndry = FacetFunction("size_t", mesh)
# Mark the boundaries
bndry.set_all(0)
for f in facets(mesh):
    mp=f.midpoint()
    if near(mp[0],0.0) : bndry[f]=1     # inflow
    if near(mp[0],L)  : bndry[f]=2      # outflow
    if near(mp[1],0.0) or near(mp[1],W) : bndry[f]=3   # walls
    if between( math.hypot(mp[0]-Center[0],mp[1]-Center[1]) , ( 0.0, Radius) ) : bndry[f]=5

ds = Measure("ds")[bndry]

# Define function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
W = MixedFunctionSpace([V, P])

# No-slip boundary condition for velocity on walls and cylinder - boundary id 3
noslip = Constant((0, 0))
bc0 = DirichletBC(W.sub(0), noslip, bndry, 3)
bc_cylinder = DirichletBC(W.sub(0), noslip, bndry, 5)

# Inflow boundary condition for velocity - boundary id 1
v_in = Expression(("0.3 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )","0.0"))
bc_in = DirichletBC(W.sub(0), v_in, bndry, 1)

# Collect boundary conditions
bcs = [bc_cylinder, bc0, bc_in]

# Define unknown and test function(s)
(_v, _p) = TestFunctions(W)

w = Function(W)
(v, p) = split(w)

n = FacetNormal(mesh)

# Define variational form for Stokes
def a(u,v): return inner(grad(u),grad(v))*dx
def b(p,v): return p*div(v)*dx
def L(v): return inner(Constant((0.0,0.0)),v)*dx

ST = a(v,_v) + b(p,_v) + b(_p,v) - L(_v)

J = derivative(ST, w)
ffc_options = {"optimize": True, "eliminate_zeros": True, "quadrature_degree": 4}
problem=NonlinearVariationalProblem(ST, w, bcs, J, form_compiler_parameters=ffc_options)
solver=NonlinearVariationalSolver(problem)

# Solve
solver.solve()

# Extract solutions:
(v, p) = w.split(True)

print "Inflow flux:  %e" % assemble(inner(v,n)*ds(1))
print "Outflow flux: %e" % assemble(inner(v,n)*ds(2))

# Save solution in PVD format
vfile = File("results_%s/v.pvd" % (fileName))
pfile = File("results_%s/p.pvd" % (fileName))
v.rename("v", "velocity")
p.rename("p", "pressure")
vfile << v
pfile << p




