from dolfin import *
import matplotlib.pyplot as plt


def step(get_data, dsN, theta, t, dt, u_old, u_new):
    """Perform timestep by theta-scheme. Given function
    get_data(t) returning data (f(t), g(t)); given solution
    u_old at time t; store solution at time t + dt into
    provided Function u_new."""

    f_n, g_n = get_data(t)
    f_np1, g_np1 = get_data(t+dt)

    # Extract function space
    V = u_new.function_space()

    # Prepare weak formulation
    u, v = TrialFunction(V), TestFunction(V)
    idt = Constant(1/dt)
    theta = Constant(theta)
    F = ( idt*(u - u_old)*v*dx
        + inner(grad(theta*u + (1-theta)*u_old), grad(v))*dx
        - (theta*f_np1 + (1-theta)*f_n)*v*dx
        - (theta*g_np1 + (1-theta)*g_n)*v*dsN
    )
    a, L = lhs(F), rhs(F)

    # Push log level
    old_level = get_log_level()
    set_log_level(LogLevel.WARNING)

    # Run the solver
    solve(a == L, u_new)

    # Pop log level
    set_log_level(old_level)


def timestepping(V, dsN, theta, T, dt, u_0, get_data):
    """Perform timestepping using theta-scheme with
    final time T, timestep dt, initial datum u_0 and
    function get_data(t) returning (f(t), g(t))"""

    # Initialize solution function
    u = Function(V)

    # Set initial condition
    u.interpolate(u_0)

    # Open plot window
    fig = init_plot()

    # Perform timestepping
    t = 0
    while t < T:
        print("t =", t, "dt =", dt)

        # Perform time step
        step(get_data, dsN, theta, t, dt, u, u)
        t += dt

        # Update plot
        update_plot(fig, u)


def timestepping_adaptive(V, dsN, theta, T, tol, u_0, get_data):
    """Perform adaptive timestepping using theta-scheme with
    final time T, tolerance tol, initial datum u_0 and
    function get_data(t) returning (f(t), g(t))"""

    # Initialize needed functions
    u_n = Function(V)
    u_np1_low = Function(V)
    u_np1_high = Function(V)

    # Initial time step; the values does not really matter
    dt = T/2

    # Set initial conditions
    u_n.interpolate(u_0)

    # Open plot window
    fig = init_plot()

    # Perform timestepping
    t = 0
    while t < T:
        print("t =", t, "dt =", dt)

        # Compute tentative time steps
        step(get_data, dsN, theta, t, dt, u_n, u_np1_low)
        step(get_data, dsN, theta, t, dt/2, u_n, u_np1_high)
        step(get_data, dsN, theta, t+dt/2, dt/2, u_np1_high, u_np1_high)

        # Compute error estimate and new timestep
        est = compute_est(theta, u_np1_low, u_np1_high)
        dt_new = compute_new_dt(theta, est, tol, dt)

        if est > tol:
            # Tolerance not met; repeat the step with new timestep
            dt = dt_new
            continue

        # Move to next time step
        u_n.vector()[:] = u_np1_high.vector()
        t += dt
        dt = dt_new

        # Update plot
        update_plot(fig, u_n)


def compute_est(theta, u_L, u_H):
    """Return error estimate by Richardson extrapolation"""
    p = 2 if theta == 0.5 else 1
    est = sqrt(assemble((u_L - u_H)**2*dx)) / (2**p - 1)
    return est


def compute_new_dt(theta, est, tol, dt):
    """Return new time step"""
    p = 2 if theta == 0.5 else 1
    rho = 0.9
    dt_new = dt * ( rho * tol / est )**(1/p)
    return dt_new


def init_plot():
    """Open plot window and return its figure object"""
    fig = plt.figure()
    fig.show()
    return fig


def update_plot(fig, u, zlims=(0, 2)):
    """Plot u in 3D warp mode with colorbar into figure fig;
    use zlims as limits on z-axis"""
    fig.clear()
    p = plot(u, mode="warp")
    fig.colorbar(p)
    fig.gca().set_zlim(zlims)
    fig.canvas.draw()


def create_function_space():
    """Return (arbitrary) H^1 conforming function space on
    unit square domain"""
    mesh = UnitSquareMesh(32, 32)
    V = FunctionSpace(mesh, "P", 1)
    return V


def create_surface_measure_left(mesh):
    """Return surface measure on the left boundary of unit
    square"""
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0)
    facets = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    Left().mark(facets, 1)
    ds_left = Measure("ds", mesh, subdomain_data=facets, subdomain_id=1)
    return ds_left


def get_data_1(t):
    """Get data for Problem 1"""
    f = Constant(0)
    g = Constant(1)
    return f, g


def get_data_2(t):
    """Get data for Problem 2"""
    f = Constant(2-t)
    g = Constant(t)
    return f, g


def get_data_3(t):
    """Get data for Problem 3"""
    f = Constant(0)
    g = Constant(max(0, 1-t)/2)
    return f, g


if __name__ == '__main__':
    # Common data
    V = create_function_space()
    ds_left = create_surface_measure_left(V.mesh())
    T = 2
    u_0 = Expression("x[0]", degree=1)

    # Problem 1, implicit Euler
    theta = 1
    dt = 0.1
    timestepping(V, ds_left, theta, T, dt, u_0, get_data_1)

    # Problem 1, explicit Euler
    theta = 0
    dt = 0.1
    timestepping(V, ds_left, theta, T, dt, u_0, get_data_1)

    # Problem 2, implicit Euler
    theta = 1
    dt = 0.1
    timestepping(V, ds_left, theta, T, dt, u_0, get_data_2)

    # Problem 3 with fixed timestep, implicit Euler
    theta = 1
    dt = 0.1
    timestepping(V, ds_left, theta, T, dt, u_0, get_data_3)

    # Problem 3 by adaptive implicit Euler
    theta = 1
    tol = 1e-3
    timestepping_adaptive(V, ds_left, theta, T, tol, u_0, get_data_3)

    # Problem 3 by adaptive Crank-Nicolson
    theta = 0.5
    tol = 1e-3
    timestepping_adaptive(V, ds_left, theta, T, tol, u_0, get_data_3)

    # Hold plots before quitting
    plt.show()
