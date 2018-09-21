from dolfin import *
import matplotlib.pyplot as plt


def solve_helmholtz(V, lambd, f):
    """Solve Helmholtz problem

        -\Delta u - lambd u = f  in \Omega
                          u = 0  on \partial\Omega

    and return u.
    """

    bc = DirichletBC(V, 0, lambda x, on_boundary: on_boundary)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx - Constant(lambd)*u*v*dx
    L = f*v*dx
    u = Function(V)
    solve(a == L, u, bc)
    return u


def build_laplacian_eigenspace(V, lambd, maxdim, tol):
    """For given space V finds eigenspace of Laplacian
    (with zero Dirichlet BC) corresponding to eigenvalues
    close to lambd by given tolerance tol. Return list
    with basis functions of the space.
    """

    # Assemble Laplacian A and mass matrix B
    bc = DirichletBC(V, 0, lambda x, on_boundary: on_boundary)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L_dummy = Constant(0)*v*dx
    m = u*v*dx
    A, _ = assemble_system(a, L_dummy, bc)
    B = assemble(m)

    # Prepare eigensolver for
    #
    #    A x = lambda B x
    eigensolver = SLEPcEigenSolver(as_backend_type(A), as_backend_type(B))
    eigensolver.parameters['problem_type'] = 'gen_hermitian'
    eigensolver.parameters['spectrum'] = 'target real'
    eigensolver.parameters['spectral_shift'] = float(lambd)
    eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
    eigensolver.parameters['tolerance'] = 1e-6
    #eigensolver.parameters['verbose'] = True  # for debugging

    # Solve for given number of eigenpairs
    eigensolver.solve(maxdim)

    # Iterate over converged eigenpairs
    space = []
    for j in range(eigensolver.get_number_converged()):

        # Get eigenpair
        r, c, rx, cx = eigensolver.get_eigenpair(j)

        # Check that eigenvalue is real
        assert near(c/r, 0, 1e-6)

        # Consider found eigenvalues close to the target eigenvalue
        if near(r, lambd, tol*lambd):
            print('Found eigenfunction with eigenvalue {} close to target {} '
                  'within tolerance {}'.format(r, lambd, tol))

            # Store the eigenfunction
            eig = Function(V)
            eig.vector()[:] = rx
            space.append(eig)

    # Check that we got whole eigenspace, i.e., last eigenvalue is different one
    assert not near(r, lambd, tol), "Possibly don't have whole eigenspace!"

    # Report
    print('Eigenspace for {} has dimension {}'.format(lambd, len(space)))

    return space


def orthogonalize(A):
    """L^2-orthogonalize a list of Functions living on the same
    function space. Modify the functions in-place.
    Use classical Gramm-Schmidt algorithm for brevity.
    For numerical stability modified Gramm-Schmidt would be better.
    """

    # Set of single function is orthogonal
    if len(A) <= 1:
        return

    # Orthogonalize overything but the last function
    orthogonalize(A[:-1])

    # Orthogonalize the last function to the previous ones
    f = A[-1]
    for v in A[:-1]:
        r = assemble(inner(f, v)*dx) / assemble(inner(v, v)*dx)
        assert f.function_space() == v.function_space()
        f.vector().axpy(-r, v.vector())


def task_1():

    # Problem data
    f = Expression('x[0] + x[1]', degree=1)
    omega2 = 5*pi**2

    # Iterate over refined meshes
    ndofs, energies = [], []
    for n in (2**i for i in range(2, 7)):

        mesh = UnitSquareMesh(n, n)
        V = FunctionSpace(mesh, "Lagrange", 1)
        u = solve_helmholtz(V, omega2, f)

        # Store energy to check convergence
        ndofs.append(u.function_space().dim())
        energies.append(norm(u, norm_type='H10'))

    # Plot energies against number dofs
    plt.plot(ndofs, energies, 'o-')
    plt.xlabel('dimension')
    plt.ylabel('energy')
    plt.show()


def tasks_2_3_4():

    # Problem data
    f = Expression('x[0] + x[1]', degree=1)
    omega2 = 5*pi**2

    # Iterate over refined meshes
    ndofs, energies = [], []
    for n in (2**i for i in range(2, 7)):

        mesh = UnitSquareMesh(n, n)
        V = FunctionSpace(mesh, "Lagrange", 1)

        # Build eingenspace of omega2
        eigenspace = build_laplacian_eigenspace(V, omega2, 10, 0.1)

        # Orthogonalize f to the eigenspace
        f_perp = project(f, V)
        orthogonalize(eigenspace+[f_perp])

        # Find particular solution with orthogonalized rhs
        u = solve_helmholtz(V, omega2, f_perp)

        # Store energy to check convergence
        ndofs.append(u.function_space().dim())
        energies.append(norm(u, norm_type='H10'))

    # Plot energies against number dofs
    plt.plot(ndofs, energies, 'o-')
    plt.xlabel('dimension')
    plt.ylabel('energy')
    plt.show()

    # Create and save w(t, x) for plotting in Paraview
    omega = omega2**0.5
    Pf = project(f - f_perp, V)
    T_per = 2*pi/omega
    create_and_save_w(omega, Pf, u, 20*T_per, 0.1*T_per)


def create_and_save_w(omega, Pf, u, T, dt):
    """Create and save w(t, x) on (0, T) with time
    resolution dt
    """

    # Extract common function space
    V = u.function_space()
    assert V == Pf.function_space()

    w = Function(V)

    Pf_vec = Pf.vector()
    u_vec = u.vector()
    w_vec = w.vector()

    f = XDMFFile('w.xdmf')
    f.parameters['rewrite_function_mesh'] = False
    f.parameters['functions_share_mesh'] = True

    def c1(t):
        return t*sin(omega*t)/(2*omega), -t*cos(omega*t)/(2*omega)

    def c2(t):
        return cos(omega*t), sin(omega*t)

    t = 0
    while t < T:
        c1_real, c1_imag = c1(t)
        c2_real, c2_imag = c2(t)

        # Store real part
        w.vector().zero()
        w.vector().axpy(c1_real, Pf_vec)
        w.vector().axpy(c2_real, u_vec)
        w.rename('w_real', 'w_real')
        f.write(w, t)

        # Store imaginary part
        w.vector().zero()
        w.vector().axpy(c1_imag, Pf_vec)
        w.vector().axpy(c2_imag, u_vec)
        w.rename('w_imag', 'w_imag')
        f.write(w, t)

        t += dt

    f.close()


if __name__ == '__main__':

    task_1()
    tasks_2_3_4()
