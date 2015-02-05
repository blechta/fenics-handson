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

# Create mesh
n =  4
height = 1.0
length = 20.0
mesh = RectangleMesh(0.0, 0.0, length, height, int(length*n), int(height*n), 'crossed')

# Create boundary markers
boundary_parts = FacetFunction('size_t', mesh)
right    = AutoSubDomain(lambda x: near(x[0], length))
left     = AutoSubDomain(lambda x: near(x[0], 0.0))
left.mark(boundary_parts, 1)
right.mark(boundary_parts, 2)
ds = Measure("ds")[boundary_parts]

U = VectorFunctionSpace(mesh, "CG", 2)
V = VectorFunctionSpace(mesh, "CG", 2)
P = FunctionSpace(mesh,"CG",1)
W = MixedFunctionSpace([U,V,P])


bcs = list()
bcs.append( DirichletBC(W.sub(0), Constant((0.0, 0.0)), boundary_parts, 1) )
bcs.append( DirichletBC(W.sub(1), Constant((0.0, 0.0)), boundary_parts, 1) )

# Material parameters
E=1e3   #Young's modulus
nu=0.3  #Poisson ratio

if nu==0.5 :   #incompressible
    mu = Constant(E/(2*(1 + nu)))
    lmbd = Constant(0.0)
    ilmbd= Constant(0.0)   
    def stress(u,p):
        F = I + grad(u)
        J=det(F)
        B = F * F.T
        C = F.T * F
        E = 0.5*(C-I)
        #W = stored energy...., then P=dW/dE
        T = -p*I+mu*(B-I)  #cauchy stress
        P = J*T*inv(F).T   #1st Piola-Kirchhoff
        pp= J-1.0   # incompressibility constraint or if W=f(E)+p*(J-1) then dW/dp
        return(P,pp)
else :   #compressible
    mu = Constant(E/(2*(1 + nu)))
    lmbd = Constant(E*nu/((1 + nu)*(1 - 2*nu)))
    ilmbd= 1.0/lmbd
    def stress(u,p):
        F = I + grad(u)
        J=det(F)
        B = F * F.T
        C = F.T * F
        E = 0.5*(C-I)
        #W = stored energy...., then P=dW/dE
        T = -p*I+mu*(B-I)   #cauchy stress
        P= J*T*inv(F).T     #1st Piola-Kirchhoff stress
        pp= ilmbd*p + (J*J-1.0)
        return(P,pp)

I=Identity(W.mesh().geometry().dim())

# Timestepping theta-method parameters
q=Constant(0.5)
dt=Constant(0.1)

w=Function(W)
(u,v,p)=split(w)

w0=Function(W)
(u0,v0,p0)=split(w0)

(_u,_v,_p)=TestFunctions(W)


(T,pp)=stress(u,p)
(T0,pp0)=stress(u0,p0)

F1= (1.0/dt)*inner(u-u0,_u)*dx - ( q*inner(v, _u)*dx + (1-q)*inner(v0,_u)*dx )

F2a= mu*inner(T,grad(_v))*dx + (pp)*_p*dx
F2b= mu*inner(T0,grad(_v))*dx + (pp0)*_p*dx
F2= (1.0/dt)*inner(v-v0,_v)*dx + ( q*F2a + (1-q)*F2b )

N = FacetNormal(mesh)
F=I+grad(u)
n = det(F)*dot(N,inv(F).T)
bFmag=Constant(0.0)
bF=bFmag*det(F)*dot(inv(F).T,Constant((0.0,1.0)))
#bF=bFmag*Constant((0.0,1.0))
G = inner(Constant((0.0,1.0)),_v)*dx
FF= inner(bF,_v)*ds(2)

F = F1 + F2 + FF

J = derivative(F, w)

step = 0
nsteps = 200

# Create files for storing solution
vfile = File("results/velo.xdmf")
ufile = File("results/disp.xdmf")
pfile = File("results/pres.xdmf")

parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, "eliminate_zeros": False, "quadrature_degree": 4 }

problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J, form_compiler_parameters=ffc_options)
solver = NonlinearVariationalSolver(problem)
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-6
solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
solver.parameters['newton_solver']['maximum_iterations'] = 20
solver.parameters['newton_solver']['report'] = True

plt=plot(u, mode="displacement", interactive=False, wireframe=True)
t=0.0
while step < nsteps:
    print "Step: " + str(step) + "Time: " + str(t)
    t+=float(dt)
    bFmag.assign(100.0*t)
    try:    
        solver.solve()
    except RuntimeError: 
        print 'diverge'
        quit()
    print 'after solve'
    
    (u,v,p)=w.split()
    u.rename("u","displacement")
    v.rename("v","velocity")
    p.rename("p","pressure")
    ufile << (u,t)
    vfile << (v,t)
    pfile << (p,t)
    plt.plot(u)

    w0.assign(w)
    step += 1
    
print 'end'
