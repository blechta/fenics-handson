from dolfin import *
import matplotlib.pyplot as plt

T = 2
dt = Constant(0.1)
theta = Constant(1.0)

u_0 = Expression("x[0]", degree=1)
f_n = Constant(0)
f_np1 = Constant(0)
g_n = Constant(0)
g_np1 = Constant(0)

def update_data_a(t, dt, f_n, f_np1, g_n, g_np1):
    f_n.assign(0)
    f_np1.assign(0)
    g_n.assign(1)
    g_np1.assign(1)

def update_data_b(t, dt, f_n, f_np1, g_n, g_np1):
    tn = t - dt
    tnp1 = t
    f_n.assign(2-tn)
    f_np1.assign(2-tnp1)
    g_n.assign(tn)
    g_np1.assign(tnp1)

def update_data_c(t, dt, f_n, f_np1, g_n, g_np1):
    tn = t - dt
    tnp1 = t
    f_n.assign(0)
    f_np1.assign(0)
    g_n.assign(max(0, 1-tn)/2)
    g_np1.assign(max(0, 1-tnp1)/2)

update_data = update_data_c

mesh = UnitSquareMesh(32, 32)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0)
facets = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
Left().mark(facets, 1)
dsN = Measure("ds", mesh, subdomain_data=facets, subdomain_id=1)

V = FunctionSpace(mesh, "P", 1)
u_np1, v = TrialFunction(V), TestFunction(V)
u_n = Function(V)
F = ( 1/dt*(u_np1 - u_n)*v*dx
    + inner(grad(theta*u_np1 + (1-theta)*u_n), grad(v))*dx
    - (theta*f_np1 + (1-theta)*f_n)*v*dx(mesh)
    - (theta*g_np1 + (1-theta)*g_n)*v*dsN
)

a, L = lhs(F), rhs(F)
u = Function(V)
problem = LinearVariationalProblem(a, L, u)
solver  = LinearVariationalSolver(problem)

u.interpolate(u_0)

fig = plt.figure()
p = plot(u, mode="warp")
plt.colorbar(p)
plt.show(block=True)

fig = plt.figure()
plt.show(block=False)

set_log_level(LogLevel.WARNING)

t = 0
while t < T:
    t += float(dt)
    u_n.vector()[:] = u.vector()

    update_data(t, float(dt), f_n, f_np1, g_n, g_np1)

    solver.solve()

    plt.clf()
    p = plot(u, mode="warp")
    plt.colorbar(p)
    fig.canvas.draw()

plt.show()
