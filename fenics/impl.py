from dolfin import *
import matplotlib.pyplot as plt


def solve_task3(V):
    """Return solution of Task 3 on space V"""

    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x, on_boundary):
        return on_boundary and x[0] > 1.0 - DOLFIN_EPS

    # Define boundary condition
    uD = Expression("x[1]", degree=1)
    bc = DirichletBC(V, uD, boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression("x[0]", degree=1)
    g = Expression("sin(5*x[0])*exp(x[1])", degree=3)
    a = inner(grad(u), grad(v))*dx + 6*u*v*dx
    L = f*v*dx + g*v*ds

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    return u


def solve_task5(V):
    """Return solution of Task 5 on space V"""

    # Define variational problem
    u = Function(V)
    v = TestFunction(V)
    f = Expression("x[0]", degree=1)
    g = Expression("sin(5*x[0])*exp(x[1])", degree=3)
    F = inner(grad(u), grad(v))*dx + (u**3 + u)*v*dx - f*v*dx - g*v*ds

    # Compute solution
    solve(F == 0, u)

    return u


def solve_task6(V):
    """Return solution of Task 6 on space V"""

    # Define Dirichlet boundary
    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    # Define boundary condition
    boundary = Boundary()
    uD = Expression("x[1]", degree=1)
    bc = DirichletBC(V, uD, boundary)

    # Define variational problem
    u = Function(V)
    v = TestFunction(V)
    c = Expression("0.1 + 0.5*(x[0]*x[0] + x[1]*x[1])", degree=2)
    f = Expression("100*x[0]", degree=1)
    F = c*inner(grad(u), grad(v))*dx + (10*u**3 + u)*v*dx - f*v*dx

    # Compute solution
    solve(F == 0, u, bc)

    return u


def solve_task7(V):
    """Return solution of Task 7 on space V"""

    # Define variational problem
    u = Function(V)
    f = Expression("x[0]", degree=1)
    g = Expression("sin(5*x[0])*exp(x[1])", degree=3)
    E = ( grad(u)**2/2 + u**4/4 + u**2/2 - f*u )*dx - g*u*ds
    F = derivative(E, u)

    # Compute solution
    solve(F == 0, u)

    return u


def solve_task8(V):
    """Return solution of Task 8 on space V"""

    # Define variational problem
    u = Function(V)
    v = TestFunction(V)
    f = Expression("0.5*(x[0] + x[1])", degree=1)
    A = as_matrix((
         (0.1 + u**2,        0),
         (         0, 1 + u**2),
    ))
    F = inner(A*grad(u), grad(v))*dx + u*v*dx - f*v*dx

    # Compute solution
    solve(F == 0, u)

    return u


if __name__ == '__main__':

    # Create mesh and define function space
    mesh = UnitSquareMesh(32, 32)
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Solve all problems
    u3 = solve_task3(V)
    u5 = solve_task5(V)
    u6 = solve_task6(V)
    u7 = solve_task7(V)
    u8 = solve_task8(V)

    # Compare solution which should be same
    err = ( grad(u7 - u5)**2 + (u7 - u5)**2 )*dx
    err = assemble(err)
    print("||u7 - u5||_H1 =", err)

    # Plot all solutions into separate figures
    plt.figure()
    plot(u3, title='u3', mode="warp")

    plt.figure()
    plot(u5, title='u5', mode="warp")

    plt.figure()
    plot(u6, title='u6', mode="warp")

    plt.figure()
    plot(u7, title='u7', mode="warp")

    plt.figure()
    plot(u8, title='u8', mode="warp")

    # Display all plots
    plt.show()
