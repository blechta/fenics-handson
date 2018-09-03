# Copyright (C) 2014 Jan Blechta
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
import matplotlib.pyplot as plt

mesh = UnitSquareMesh(40, 40)
V = FunctionSpace(mesh, 'Lagrange', 1)
f = Expression("1.+cos(2.*pi*x[0])*sin(2.*pi*x[1])", degree=2)

def p_laplace(p, eps, u0=None):
    """Solves regularized p-Laplacian with mesh, space and right-hand side
    defined above. Returns solution, regularized energy and energy."""
    p = Constant(p)
    eps = Constant(eps)

    # Initial approximation for Newton
    u = u0.copy(deepcopy=True) if u0 else Function(V)

    # Energy functional
    E = ( 1./p*(eps + dot(grad(u), grad(u)))**(0.5*p) - f*u ) * dx

    # First Gateaux derivative
    F = derivative(E, u)

    # Solve nonlinear problem; Newton is used by default
    bc = DirichletBC(V, 0.0, lambda x,onb: onb)
    solver_parameters = {'newton_solver': {'maximum_iterations': 1000}}
    solve(F == 0, u, bc, solver_parameters=solver_parameters)

    plt.gcf().show()
    plt.clf()
    plot(u, mode="warp", title='p-Laplace, p=%g, eps=%g'%(float(p), float(eps)))
    plt.gcf().canvas.draw()

    # Compute energies
    energy_regularized = assemble(E)
    eps.assign(0.0)
    energy = assemble(E)

    return u, energy_regularized, energy
