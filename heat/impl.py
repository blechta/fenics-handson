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

# Parameters
K = Constant(1.0e-2)
t_end = 10
dt = 0.1

# Create mesh and build function space
mesh = UnitSquareMesh(40, 40, 'crossed')
V = FunctionSpace(mesh, "Lagrange", 1)

# Create boundary markers
left = AutoSubDomain(lambda x: near(x[0], 0.0))
left.mark(boundary_parts, 1)

# Define boundary measure
# TODO: This is deprecated syntax
ds = Measure("ds")[boundary_parts]


# Initial condition and right-hand side
ic = Expression("((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))<0.2*0.2)?(-25*((pow(x[0]-0.25,2)+pow(x[1]-0.25,2))-0.2*0.2)):(0.0)")
f = Expression("((pow(x[0]-0.75,2)+pow(x[1]-0.75,2))<0.2*0.2)?(1.0):(0.0)")

# Define unknown and test function
u = TrialFunction(V)
v = TestFunction(V)
u0 = Function(V)

# Equation coefficients
g = Constant(0.1)
b = Expression( ("-(x[1]-0.5)","(x[0]-0.5)") )
theta = Constant(1.0)

# Define boundary condition
bc = DirichletBC(V, Constant(0.0), "near(x[0], 1.0) || near(x[1], 0.0)")

# Define variational forms
def operator(u, v):
    return ( K*inner(grad(u), grad(v)) - f*v )*dx - g*v*ds(1)
F = (1.0/dt)*inner(u-u0, v)*dx \
  + theta*operator(u, v) \
  + (1.0-theta)*operator(u0, v)

# Create file for storing results
f = File("results/u.xdmf")

# Prepare initial condition
u0.interpolate(ic)

# Prepare solution function and solver
u = Function(V)
problem = LinearVariationalProblem(lhs(F), rhs(F), u, bc)
solver  = LinearVariationalSolver(problem)

# Time-stepping
t = 0.0
while t < t_end:

    # Solve the problem
    solver.solve()

    # Store solution to file and plot
    f << u
    plot(u, title='Solution at t = %g' % t)

    # Move to next time step
    u0.assign(u)
    t += dt

    # Report flux
    flux = assemble(K*grad(u)*ds(1))
    info('t = %g, flux = %g' % (t, flux))
