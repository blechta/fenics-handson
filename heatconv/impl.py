from dolfin import *
import matplotlib.pyplot as plt

# Create mesh and build function space
mesh = UnitSquareMesh(40, 40, 'crossed')
V = FunctionSpace(mesh, "Lagrange", 1)

# Create boundary markers
tdim = mesh.topology().dim()
boundary_parts = MeshFunction('size_t', mesh, tdim-1)
left   = AutoSubDomain(lambda x: near(x[0], 0.0))
right  = AutoSubDomain(lambda x: near(x[0], 1.0))
bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
left  .mark(boundary_parts, 1)
right .mark(boundary_parts, 2)
bottom.mark(boundary_parts, 2)

# Initial condition and right-hand side
ic = Expression("""pow(x[0] - 0.25, 2) + pow(x[1] - 0.25, 2) < 0.2*0.2
                   ? -25.0 * ((pow(x[0] - 0.25, 2) + pow(x[1] - 0.25, 2)) - 0.2*0.2)
                   : 0.0""", degree=1)
f = Expression("""pow(x[0] - 0.75, 2) + pow(x[1] - 0.75, 2) < 0.2*0.2
                  ? 1.0
                  : 0.0""", degree=1)

# Equation coefficients
K = Constant(1e-2) # thermal conductivity
g = Constant(0.01) # Neumann heat flux
b = Expression(("-(x[1] - 0.5)", "x[0] - 0.5"), degree=1) # convecting velocity

# Define boundary measure on Neumann part of boundary
dsN = Measure("ds", subdomain_id=1, subdomain_data=boundary_parts)

# Define steady part of the equation
def operator(u, v):
    return ( K*inner(grad(u), grad(v)) - f*v + dot(b, grad(u))*v )*dx - K*g*v*dsN

# Define trial and test function and solution at previous time-step
u = TrialFunction(V)
v = TestFunction(V)
u0 = Function(V)

# Time-stepping parameters
t_end = 10
dt = 0.1
theta = Constant(0.5) # Crank-Nicolson scheme

# Define time discretized equation
F = (1.0/dt)*inner(u-u0, v)*dx + theta*operator(u, v) + (1.0-theta)*operator(u0, v)

# Define boundary condition
bc = DirichletBC(V, Constant(0.0), boundary_parts, 2)

# Prepare solution function and solver
u = Function(V)
problem = LinearVariationalProblem(lhs(F), rhs(F), u, bc)
solver  = LinearVariationalSolver(problem)

# Prepare initial condition
u0.interpolate(ic)

# Create file for storing results
f = XDMFFile("results/u.xdmf")

# Time-stepping
t = 0.0
u.rename("u", "temperature")
u.interpolate(ic)

# Save initial solution
f.write(u, t)

# Open figure for plots
fig = plt.figure()
plt.show(block=False)

while t < t_end:

    # Solve the problem
    solver.solve()

    # Store solution to file and plot
    f.write(u, t)
    p = plot(u, title='Solution at t = %g' % t)
    if p is not None:
        plt.colorbar(p)
    fig.canvas.draw()
    plt.clf()

    # Move to next time step
    u0.assign(u)
    t += dt

    # Report flux
    n = FacetNormal(mesh)
    flux = assemble(K*dot(grad(u), n)*dsN)
    info('t = %g, flux = %g' % (t, flux))
