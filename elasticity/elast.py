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

# Use UFLACS to speed-up assembly
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4


def solve_elasticity(facet_function, E, nu, dt, T_end, output_dir):

    # Get mesh and prepare boundary measure
    mesh = facet_function.mesh()
    ds = Measure("ds", subdomain_data=facet_function)

    # Build function space
    U = VectorFunctionSpace(mesh, "Lagrange", 2)
    P = FunctionSpace(mesh, "Lagrange", 1)
    W = MixedFunctionSpace([U, U, P])

    # Prepare BCs
    bcs = [DirichletBC(W.sub(i), gdim*(0.0,), facet_function, 1)
           for i in [0, 1]]

    # Define constitutive law
    mu = Constant(E/(2.0*(1.0 + nu)))
    def stress(u, p):
        F = I + grad(u)
        J = det(F)
        B = F * F.T
        T = -p*I + mu*(B-I) # Cauchy stress
        if nu == 0.5:
            # Incompressible
            pp = J-1.0
        else:
            # Compressible
            lmbd = Constant(E*nu/((1.0 + nu)*(1.0 - 2.0*nu)))
            pp = 1.0/lmbd*p + (J*J-1.0)
        return T, pp

    # Timestepping theta-method parameters
    q = Constant(0.5)
    dt = Constant(dt)

    # Unknowns, values at previous step and test functions
    w = Function(W)
    (u, v, p) = split(w)
    w0 = Function(W)
    (u0, v0, p0) = split(w0)
    (_u, _v, _p) = TestFunctions(W)

    I = Identity(W.mesh().geometry().dim())

    # Balance of momentum
    T, pp = stress(u, p)
    T0, pp0 = stress(u0, p0)
    F1 = (1.0/dt)*inner(u-u0, _u)*dx - q*inner(v, _u)*dx + (1.0-q)*inner(v0, _u)*dx
    F2a = mu*inner(T, grad(_v))*dx + pp*_p*dx
    F2b = mu*inner(T0, grad(_v))*dx + pp0*_p*dx
    F2 = (1.0/dt)*inner(v-v0, _v)*dx + q*F2a + (1.0-q)*F2b

    # Traction at boundary
    F = I + grad(u)
    bF_magnitude = Constant(0.0)
    bF_direction = {2: Constant((0.0, 1.0)), 3: Constant((0.0, 0.0, 1.0))}[gdim]
    bF = det(F)*dot(inv(F).T, bF_magnitude*bF_direction)
    FF = inner(bF, _v)*ds(2)

    # Whole system and its Jacobian
    F = F1 + F2 + FF
    J = derivative(F, w)

    # Initialize solver
    problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
    solver.parameters['newton_solver']['linear_solver'] = 'mumps'

    # Extract solution components
    (u, v, p) = w.split()
    u.rename("u", "displacement")
    v.rename("v", "velocity")
    p.rename("p", "pressure")

    # Create files for storing solution
    vfile = File("%s/velo.xdmf" % output_dir)
    ufile = File("%s/disp.xdmf" % output_dir)
    pfile = File("%s/pres.xdmf" % output_dir)

    # Prepare plot window
    plt = plot(u, mode="displacement", interactive=False, wireframe=True)

    # Time-stepping loop
    t = 0.0
    while t <= T_end:
        print "Time: %g" % t
        t += float(dt)

        # Increase traction
        bF_magnitude.assign(100.0*t)

        # Prepare to solve and solve
        w0.assign(w)
        solver.solve()

        # Store solution to files and plot
        ufile << (u, t)
        vfile << (v, t)
        pfile << (p, t)
        plt.plot(u)


def geometry_2d():
    n = 4
    x0 = 0.0
    x1 = 20.0
    y0 = 0.0
    y1 = 1.0
    mesh = RectangleMesh(x0, y0, x1, y1, int((x1-x0)*n), int((y1-y0)*n), 'crossed')
    boundary_parts = FacetFunction('size_t', mesh)
    left  = AutoSubDomain(lambda x: near(x[0], x0))
    right = AutoSubDomain(lambda x: near(x[0], x1))
    left .mark(boundary_parts, 1)
    right.mark(boundary_parts, 2)
    return boundary_parts


def geometry_3d():
    mesh = Mesh('lego_beam.xml')
    gdim = mesh.geometry().dim()
    x0 = mesh.coordinates()[:, 0].min()
    x1 = mesh.coordinates()[:, 0].max()
    boundary_parts = FacetFunction('size_t', mesh)
    left  = AutoSubDomain(lambda x: near(x[0], x0))
    right = AutoSubDomain(lambda x: near(x[0], x1))
    left .mark(boundary_parts, 1)
    right.mark(boundary_parts, 2)
    return boundary_parts


if __name__ == '__main__':

    solve_elasticity(geometry_2d(), 1e3, 0.3, 0.1, 20.0)
    solve_elasticity(geometry_2d(), 1e3, 0.5, 0.1, 20.0)
    solve_elasticity(geometry_3d(), 1e3, 0.3, 0.1, 20.0)
