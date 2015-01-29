# Copyright (C) 2015 Jan Blechta
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

# Build mesh and function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, 'Lagrange', 3)

# Right-hand side
f = Expression("sin(6.0*pi*x[0])*sin(2.0*pi*x[1])")

# Define boundary condition
def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, 0.0, boundary)

# Prepare variational formulation
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

# Init function and solve variational problem
u = Function(V)
solve(a == L, u, bc)

# Plot solution
plot(u, interactive=True)


# Now for Neumann
def boundary(x, on_boudary):
    return near(x[1], 0.0) or near(x[1], 1.0)
bc = DirichletBC(V, 0.0, boundary)
solve(a == L, u, bc)
plot(u, interactive=True)


# And now for non-linear elliptic operator
# TODO: This implementation does not stress that linear form (without
#       TrialFunction) with non-linearity in Function is needed. But
#       attendees will probably encounter this when modifying the code.
k = 1e6
F = ( (1.0 + Constant(k)*u*u)*inner(grad(u), grad(v)) - f*v )*dx
solve(F == 0, u, bc)
plot(u, interactive=True)
