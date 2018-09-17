from dolfin import *
import mshr
import matplotlib.pyplot as plt

# Define domain
center = Point(0.2, 0.2)
radius = 0.05
L = 2.2
W = 0.41
geometry = mshr.Rectangle(Point(0.0, 0.0), Point(L, W)) \
         - mshr.Circle(center, radius, 128)

# Build mesh
mesh = mshr.generate_mesh(geometry, 64)
plot(mesh)
plt.show()

# Construct facet markers
bndry = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
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

# Dump facet markers to file
with XDMFFile('facets.xdmf') as f:
    f.write(bndry)

# Build function spaces (Taylor-Hood)
P2 = VectorElement("P", mesh.ufl_cell(), 2)
P1 = FiniteElement("P", mesh.ufl_cell(), 1)
TH = MixedElement([P2, P1])
W = FunctionSpace(mesh, TH)

# No-slip boundary condition for velocity on walls and cylinder - boundary id 3
noslip = Constant((0, 0))
bc_walls = DirichletBC(W.sub(0), noslip, bndry, 3)
bc_cylinder = DirichletBC(W.sub(0), noslip, bndry, 5)

# Inflow boundary condition for velocity - boundary id 1
v_in = Expression(("0.3 * 4.0 * x[1] * (0.41 - x[1]) / ( 0.41 * 0.41 )", "0.0"),
                  degree=2)
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

    # Report drag and lift
    p = w.sub(1, deepcopy=False)
    T = -p*I + nu*grad(u)
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


def navier_stokes():

    v, q = TestFunctions(W)
    w = Function(W)
    u, p = split(w)

    # Define variational forms
    T = -p*I + 2.0*nu*sym(grad(u))
    F = inner(T, grad(v))*dx - q*div(u)*dx + inner(grad(u)*u, v)*dx
    F += - nu*dot(dot(grad(u), v), n)*ds(2)

    solve(F == 0, w, bcs)

    # Report in,out fluxes
    u = w.sub(0, deepcopy=False)
    info("Inflow flux:  %e" % assemble(inner(u, n)*ds(1)))
    info("Outflow flux: %e" % assemble(inner(u, n)*ds(2)))

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
        print("Running '%s'" % problem.__name__)

        # Call solver
        w = problem()

        # Extract solutions
        u, p = w.split()

        # Save solution
        u.rename("u", "velocity")
        p.rename("p", "pressure")
        with XDMFFile("results_%s/u.xdmf" % problem.__name__) as f:
            f.write(u)
        with XDMFFile("results_%s/p.xdmf" % problem.__name__) as f:
            f.write(p)

        plt.figure()
        plot(u, title='velocity %s' % problem.__name__)
        plt.figure()
        plot(p, title='pressure %s' % problem.__name__)
        plt.show()
