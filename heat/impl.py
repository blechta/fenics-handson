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

# Create mesh and build function space
mesh = UnitSquareMesh(40, 40, 'crossed')
V = FunctionSpace(mesh, "Lagrange", 1)

# Create boundary markers
boundary_parts = FacetFunction('size_t', mesh)
left   = AutoSubDomain(lambda x: near(x[0], 0.0))
right  = AutoSubDomain(lambda x: near(x[0], 1.0))
bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
left  .mark(boundary_parts, 1)
right .mark(boundary_parts, 2)
bottom.mark(boundary_parts, 2)

# Initial condition and right-hand side
ic = Expression("""pow(x[0] - 0.25, 2) + pow(x[1] - 0.25, 2) < 0.2*0.2
                   ? -25.0 * ((pow(x[0] - 0.25, 2) + pow(x[1] - 0.25, 2)) - 0.2*0.2)
                   : 0.0""")
f = Expression("""pow(x[0] - 0.75, 2) + pow(x[1] - 0.75, 2) < 0.2*0.2
                  ? 1.0
                  : 0.0""")

# Equation coefficients
K = Constant(1e-2) # thermal conductivity
g = Constant(0.1) # Neumann heat flux
b = Expression(("-(x[1] - 0.5)", "x[0] - 0.5")) # convecting velocity

# Define boundary measure on Neumann part of boundary
dsN = Measure("ds", subdomain_id=1, subdomain_data=boundary_parts)

# Define steady part of the equation
def operator(u, v):
    return ( K*inner(grad(u), grad(v)) - f*v )*dx - g*v*dsN

# Define trial and test function and solution at previous time-step
u = TrialFunction(V)
v = TestFunction(V)
u0 = Function(V)

# Time-stepping parameters
t_end = 10
dt = 0.1
theta = Constant(0.5) # Crank-Nicolson scheme

# Define time discretized equation
F = (1.0/dt)*inner(u-u0, v)*dx \
  + theta*operator(u, v) \
  + (1.0-theta)*operator(u0, v)

# Define boundary condition
bc = DirichletBC(V, Constant(0.0), boundary_parts, 2)

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
    n = FacetNormal(mesh)
    flux = assemble(K*dot(grad(u), n)*dsN)
    info('t = %g, flux = %g' % (t, flux))
