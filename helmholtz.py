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

from dolfin import *


def helmholtz_illposed(n):
    """For given mesh division 'n' solves ill-posed problem
        (-Laplace - 5*pi^2) u = x + y   on [0, 1]*[0, 1],
                            u = 0       on boundary,
    and returns space dimension, energy_error (on discrete subspace) and energy."""
    mesh = UnitSquareMesh(n, n)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    bc = DirichletBC(V, 0.0, lambda x, b: b)
    u, v = TrialFunction(V), TestFunction(V)
    a = (inner(grad(u), grad(v)) - Constant(5.0*pi*pi)*u*v)*dx
    L = Expression('x[0] + x[1]')*v*dx
    u = Function(V)
    solve(a == L, u, bc)
    energy_error = assemble(action(Constant(1.0)*action(a, u) - L, u))
    energy       = assemble(action(Constant(0.5)*action(a, u) - L, u))
    return V.dim(), energy_error, energy


def helmholtz_wellposed(n):
    """For given mesh division 'n' solves well-posed problem
        (-Laplace - 5*pi^2) u = f       on [0, 1]*[0, 1],
                            u = 0       on boundary,
    with f orthogonal to eigenspace of 5*pi^2.
    and returns space dimension, energy_error (on discrete subspace) and energy."""
    # Assemble Laplacian
    mesh = UnitSquareMesh(n, n)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    bc = DirichletBC(V, 0.0, lambda x, b: b)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = Constant(0.0)*v*dx
    m = u*v*dx
    A, _ = assemble_system(a, L, bc)
    B = assemble(m)

    # Search for eigenspace for eigenvalue around 5*pi*pi
    eigensolver = SLEPcEigenSolver(as_backend_type(A), as_backend_type(B))
    eigensolver.parameters['problem_type'] = 'gen_hermitian'
    eigensolver.parameters['spectrum'] = 'target real'
    eigensolver.parameters['spectral_shift'] = 5.0*pi*pi
    eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
    eigensolver.parameters['tolerance'] = 1e-6
    #eigensolver.parameters['verbose'] = True
    eigensolver.solve(5)
    print eigensolver.get_iteration_number(), V.dim(), eigensolver.get_number_converged()
    eig = Function(V)
    temp = eig.vector()
    space = []
    for j in range(eigensolver.get_number_converged()):
        r, c, rx, cx = eigensolver.get_eigenpair(j)
        assert near(c/r, 0.0, 1e-6)
        assert near(cx.norm('linf')/rx.norm('linf'), 0.0, 1e-6)
        temp[:] = rx
        if near(r, 5*pi*pi, 0.5*pi*pi):
            space.append(eig.copy(deepcopy=True))
    # Check that we got whole eigenspace - last eigenvalue is different one
    assert not near(r, 5*pi*pi, 0.5*pi*pi)
    print 'Eigenspace for 5*pi^2 has dimension', len(space)

    # Orthogonalize right-hand side to 5*pi^2 eigenspace
    f = Expression('x[0] + x[1]')
    f = project(f, V)
    orthogonalize(space+[f])
    L = f*v*dx

    # Assemble Helmholtz system
    a = (inner(grad(u), grad(v)) - Constant(5.0*pi*pi)*u*v)*dx
    A, b = assemble_system(a, L, bc)

    # Solve well-posed resonant Helmoltz system
    u = Function(V)
    x = u.vector()
    solver = PETScLUSolver('mumps')
    solver.set_operator(A)
    solver.parameters['symmetric'] = True
    solver.parameters['reuse_factorization'] = True
    solver.solve(x, b)

    energy_error = assemble(action(Constant(1.0)*action(a, u) - L, u))
    energy       = assemble(action(Constant(0.5)*action(a, u) - L, u))
    return V.dim(), energy_error, energy


def orthogonalize(A):
    """L^2-orthogonalizes set of Functions A. Stores the result in-place to A.
    Uses classical Gramm-Schmidt algorithm for brevity. For numerical stability
    modified Gramm-Schmidt would be better."""
    assert all(isinstance(a, Function) for a in A)
    if len(A) <= 1:
        return
    orthogonalize(A[:-1])
    f = A[-1]
    for v in A[:-1]:
        # NOTE: L^2 inner product could be preassembled to reduce computation
        #       of r to matvecs.
        r = assemble(inner(f, v)*dx)/assemble(inner(v, v)*dx)
        if f.function_space() == v.function_space():
            f.vector().axpy(-r, v.vector())
        else:
            raise NotImplementedError


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    # Demonstrate that energy of ill-posed Helmholtz goes to minus infinity
    results = np.array(map(helmholtz_illposed, [2**i for i in range(2, 9)]))
    plt.subplot(2, 1, 1)
    plt.plot(results[:, 0], results[:, 2], 'o-')
    plt.title('Ill-posed Helmholtz')
    plt.xlabel('dimension')
    plt.ylabel('energy')
    plt.show(block=False)

    # Demonstrate that energy of well-posed Helmholtz converges
    results = np.array(map(helmholtz_wellposed, [2**i for i in range(2, 9)]))
    plt.subplot(2, 1, 2)
    plt.plot(results[:, 0], results[:, 2], 'o-')
    plt.title('Well-posed Helmholtz')
    plt.xlabel('dimension')
    plt.ylabel('energy')
    plt.tight_layout()
    plt.show(block=True)
