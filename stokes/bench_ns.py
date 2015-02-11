# Copyright (C) 2014, 2015 Jaroslav Hron, Jan Blechta
#
# This file is part of FEniCS tutorial suite.
#
# The suite is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the suite.  If not, see <http://www.gnu.org/licenses/>.

# Begin code
from dolfin import *
import mshr

# Define domain
center = Point(0.2, 0.2)
radius = 0.05
L = 2.2
W = 0.41
geometry = mshr.Rectangle(Point(0.0,0.0), Point(L, W)) \
         - mshr.Circle(center, radius, 10)

# Build mesh
mesh = mshr.generate_mesh(geometry, 50)

# Material parameters
nu = Constant(0.001)

# Construct facet markers
bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
    mp = f.midpoint()
    if near(mp[0], 0.0): # inflow
        bndry[f] = 1
    elif near(mp[0], L): # outflow
        bndry[f] = 2
    elif near(mp[1], 0.0) or near(mp[1], W): # walls
        bndry[f] = 3
    elif mp.distance(center) <= radius: # cylinder
        bndry[f] = 5

ds = Measure("ds")[bndry]


# Define function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh, "CG", 1)
W = MixedFunctionSpace([V, P])

# No-slip boundary condition for velocity on walls and cylinder - boundary id 3
noslip = Constant((0, 0))
bc_walls = DirichletBC(W.sub(0), noslip, bndry, 3)
bc_cylinder = DirichletBC(W.sub(0), noslip, bndry, 5)

# Inflow boundary condition for velocity - boundary id 1
v_in = Expression(("0.3 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"))
bc_in = DirichletBC(W.sub(0), v_in, bndry, 1)

# Collect boundary conditions
bcs = [bc_cylinder, bc_walls, bc_in]

# Define unknown and test function(s)
(_v, _p) = TestFunctions(W)
w = Function(W)
(v, p) = split(w)

n = FacetNormal(mesh)
I = Identity(v.geometric_dimension())    # Identity tensor

D = 0.5*(grad(v)+grad(v).T)
S = 2*nu*D
T = -p*I + S

# Define variational forms
F = inner(T, grad(_v))*dx + _p*div(v)*dx + inner(grad(v)*v, _v)*dx 

J = derivative(F, w)
problem=NonlinearVariationalProblem(F,w,bcs,J)
solver=NonlinearVariationalSolver(problem)

# Solve
solver.solve()

# Extract solutions:
(v, p) = w.split(True)

force = T*n
D=(force[0]/0.002)*ds(5)
L=(force[1]/0.002)*ds(5)

# Assemble functionals
drag = assemble(D)
lift = assemble(L)

print "drag= %e    lift= %e" % (drag , lift)

# Save solution
vfile = File("results/v.xdmf")
pfile = File("results/p.xdmf")
v.rename("v", "velocity")
p.rename("p", "pressure")
vfile << v
pfile << p
