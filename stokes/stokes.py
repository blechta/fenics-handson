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
geometry = mshr.Rectangle(Point(0.0, 0.0), Point(L, W)) \
         - mshr.Circle(center, radius, 10)

# Build mesh
mesh = mshr.generate_mesh(geometry, 50)

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

# Build function spaces (Taylor-Hood)
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

# Facet normal, identity tensor and boundary measure
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
ds = Measure("ds", subdomain_data=bndry)
nu = Constant(0.001)


def stokes():

    # Define variational forms for Stokes
    def a(u,v):
        return inner(nu*grad(u), grad(v))*dx
    def b(p,v):
        return p*div(v)*dx
    def L(v):
        return inner(Constant((0.0,0.0)), v)*dx

    # Solve the problem
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    w = Function(W)
    solve(a(u, v) - b(p, v) - b(q, u) == L(v), w, bcs)

    # Report in,out fluxes
    u = w.sub(0, deepcopy=False)
    info("Inflow flux:  %e" % assemble(inner(u, n)*ds(1)))
    info("Outflow flux: %e" % assemble(inner(u, n)*ds(2)))

    return w


def navier_stokes():

    v, q = TestFunctions(W)
    w = Function(W)
    u, p = split(w)

    # Define variational forms
    T = -p*I + 2.0*nu*sym(grad(u))
    F = inner(T, grad(v))*dx - q*div(u)*dx + inner(grad(u)*u, v)*dx

    solve(F == 0, w, bcs)

    # Report drag and lift
    force = dot(T, n)
    D = (force[0]/0.002)*ds(5)
    L = (force[1]/0.002)*ds(5)
    drag = assemble(D)
    lift = assemble(L)
    info("drag= %e    lift= %e" % (drag , lift))

    # Report pressure difference
    a_1 = Point(0.15, 0.2)
    a_2 = Point(0.25, 0.2)
    p = w.sub(1, deepcopy=False)
    p_diff = p(a_1) - p(a_2)
    info("p_diff = %e" % p_diff)

    return w


if __name__ == "__main__":

    for problem in [stokes, navier_stokes]:
        begin("Running '%s'" % problem.__name__)

        # Call solver
        w = problem()

        # Extract solutions
        u, p = w.split()

        # Save solution
        ufile = File("results_%s/u.xdmf" % problem.__name__)
        pfile = File("results_%s/p.xdmf" % problem.__name__)
        u.rename("u", "velocity")
        p.rename("p", "pressure")
        ufile << u
        pfile << p

        plot(u, title='velocity %s' % problem.__name__)
        plot(p, title='pressure %s' % problem.__name__)

        end()

    interactive()
