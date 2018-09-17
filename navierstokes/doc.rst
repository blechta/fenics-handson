Navier-Stokes equations
=======================

.. sidebar:: Goals

    Learn how to deal with mixed finite elements.
    Get remembered how fragile can numerical solutions be.
    Reproduce some cool physics -- Kármán vortex street.

.. todo::

    Add a task to analyze discretization error in the lift!

.. todo::

    Add :math:`\mathrm{Re}=100`!


Stokes flow around cylinder
---------------------------

Solve the following linear system of PDEs

.. math::

   - \operatorname{div}(\nabla u) + \nabla p &= f
        &&\quad\text{ in }\Omega,

   \operatorname{div} u &= 0
        &&\quad\text{ in }\Omega,

   u &= 0
        &&\quad\text{ on }\Gamma_\mathrm{D},

   u &= u_\mathrm{IN}
        &&\quad\text{ on }\Gamma_\mathrm{IN},

   \tfrac{\partial u}{\partial\mathbf{n}} &= g
        &&\quad\text{ on }\Gamma_\mathrm{N},

using FE discretization with data

* :math:`\Omega = (0, 2.2)\times(0, 0.41) - B_{0.05}\left((0.2,0.2)\right)`
* :math:`\Gamma_\mathrm{N} = \left\{ x = 2.2 \right\}`
* :math:`\Gamma_\mathrm{IN} = \left\{ x = 0.0 \right\}`
* :math:`\Gamma_\mathrm{D} = \Gamma_\mathrm{W} \cup \Gamma_\mathrm{S}`
* :math:`u_\mathrm{IN} = \left( 0.3 \frac{4}{0.41^2} y (0.41-y) , 0 \right)`

where :math:`B_R(\mathbf{z})` is a disc of radius :math:`R` and center
:math:`\mathbf{z}`

  .. image:: geometry.png
     :align: center
     :width: 70%
     :target: http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html

..


.. admonition:: Task 1

  Write the variational formulation of the problem and
  discretize the equation by a mixed finite element method.


.. admonition:: Task 2

  Build a mesh, prepare a mesh function marking
  :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
  check its correctness.

  .. hint::

      Use the FEniCS meshing tool ``mshr``, see `mshr documentation
      <https://bitbucket.org/benjamik/mshr/wiki/API>`_.

      .. code-block:: python

         import mshr

         # Define domain
         center = Point(0.2, 0.2)
         radius = 0.05
         L = 2.2
         W = 0.41
         geometry =  mshr.Rectangle(Point(0.0,0.0), Point(L, W)) \
                    -mshr.Circle(center, radius, 10)

         # Build mesh
         mesh = mshr.generate_mesh(geometry, 50)


  .. hint::

      Try yet another way to mark the boundaries by direct
      access to the mesh entities by
      `vertices(mesh) <dolfin.cpp.mesh.vertices>`,
      `facets(mesh) <dolfin.cpp.mesh.facets>`,
      `cells(mesh) <dolfin.cpp.mesh.cells>`
      mesh-entity iterators::

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


.. admonition:: Task 3

    Construct the mixed finite element space and the
    bilinear and linear forms together with the `DirichletBC <dolfin.fem.bcs.DirichletBC>` object.

    .. hint::

        Use for example the stable Taylor-Hood finite elements::

            # Build function spaces (Taylor-Hood)
            P2 = VectorElement("P", mesh.ufl_cell(), 2)
            P1 = FiniteElement("P", mesh.ufl_cell(), 1)
            TH = MixedElement([P2, P1])
            W = FunctionSpace(mesh, TH)

    .. hint::

        To define Dirichlet BC on subspace use the
        `W.sub() <dolfin.functions.functionspace.FunctionSpace.sub>` method::

            noslip = Constant((0, 0))
            bc_walls = DirichletBC(W.sub(0), noslip, bndry, 3)

    .. hint::

        To build the forms use::

            # Define trial and test functions
            u, p = TrialFunctions(W)
            v, q = TestFunctions(W)

        Then you can define forms on mixed space using
        ``u``, ``p``, ``v``, ``q`` as usual.


.. admonition:: Task 4

    Now modify the problem to the Navier-Stokes equations
    and compute the `DFG-flow around cylinder benchmark
    <http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html>`_

    .. hint::

        You can use generic `solve <dolfin.fem.solving.solve>` function or
        `NonlinearVariationalProblem <dolfin.fem.solving.NonlinearVariationalProblem>`
        and `NonlinearVariationalSolver <dolfin.cpp.fem.NonlinearVariationalSolver>`
        classes::


            # Define test functions
            v, q = TestFunctions(W)
            w = Function(W)
            u, p = split(w)

            # Facet normal, identity tensor and boundary measure
            n = FacetNormal(mesh)
            I = Identity(mesh.geometry().dim())
            ds = Measure("ds", subdomain_data=bndry)
            nu = Constant(0.001)

            # Define variational forms
            T = -p*I + 2.0*nu*sym(grad(u))
            F = inner(T, grad(v))*dx - q*div(u)*dx + inner(grad(u)*u, v)*dx
            F += - nu*dot(dot(grad(u), v), n)*ds(2)

    .. hint::

        Use `assemble <dolfin.fem.assembling.assemble>` function
        to evaluate the lift and drag functionals::


            # Report drag and lift
            force = dot(T, n)
            D = (force[0]/0.002)*ds(5)
            L = (force[1]/0.002)*ds(5)
            drag = assemble(D)
            lift = assemble(L)
            info("drag= %e    lift= %e" % (drag , lift))


.. only:: priv

    Reference solution
    ------------------
    .. toggle-header::
        :header: **Show/Hide Code**

        .. literalinclude:: stokes.py
