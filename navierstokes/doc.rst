Navier-Stokes equations
=======================

.. sidebar:: Goals

    Learn how to deal with mixed finite elements.
    Remember how fragile can numerical solutions be.
    Reproduce some cool physics -- Kármán vortex street.


Stokes flow around cylinder
---------------------------

Solve the following linear system of PDEs

.. math::
    :label: stokes

    - \nu\Delta\mathbf{u} + \nabla p &= \mathbf{0}
        &&\quad\text{ in }\Omega,

    \operatorname{div}\mathbf{u} &= 0
        &&\quad\text{ in }\Omega,

    \mathbf{u} &= 0
        &&\quad\text{ on }\Gamma_\mathrm{D},

    \mathbf{u} &= \mathbf{u}_\mathrm{IN}
        &&\quad\text{ on }\Gamma_\mathrm{IN},

    \nu\tfrac{\partial\mathbf{u}}{\partial\mathbf{n}} -p\mathbf{n} &= 0
        &&\quad\text{ on }\Gamma_\mathrm{N}

using FE discretization with data

.. math::
    :label: dfg20

    \Omega &= (0, 2.2)\times(0, 0.41) - B_{0.05}\left((0.2,0.2)\right),

    \Gamma_\mathrm{N} &= \left\{ x = 2.2 \right\} = \text{(green)},

    \Gamma_\mathrm{IN} &= \left\{ x = 0.0 \right\} = \text{(red)},

    \Gamma_\mathrm{D} &= \partial\Omega\setminus(\Gamma_\mathrm{N}\cup\Gamma_\mathrm{IN}) = \text{(black)},

    u_\mathrm{IN} &= \left( \frac{4Uy (0.41-y)}{0.41^2}  , 0 \right),

    \nu &= 0.001, \qquad U = 0.3

where :math:`B_R(\mathbf{z})` is a disc of radius :math:`R` and center
:math:`\mathbf{z}`

  .. image:: geometry.png
     :align: center
     :width: 80%
     :target: http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html


.. admonition:: Task 1

    Write the weak formulation of the problem and
    a spatial discretization by a mixed finite element method.


.. admonition:: Task 2

    Build a mesh, prepare a mesh function marking :math:`\Gamma_\mathrm{IN}`,
    :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
    check its correctness.

    .. hint::

        Use the FEniCS meshing tool ``mshr``, see `mshr documentation
        <https://bitbucket.org/benjamik/mshr/wiki/API>`_.

        .. code-block:: python

            import mshr

            # Discretization parameters
            N_circle = 16
            N_bulk = 64

            # Define domain
            center = Point(0.2, 0.2)
            radius = 0.05
            L = 2.2
            W = 0.41
            geometry =  mshr.Rectangle(Point(0.0, 0.0), Point(L, W)) \
                       -mshr.Circle(center, radius, N_circle)

            # Build mesh
            mesh = mshr.generate_mesh(geometry, N_bulk)


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
                if near(mp[0], 0.0):  # inflow
                    bndry[f] = 1
                elif near(mp[0], L):  # outflow
                    bndry[f] = 2
                elif near(mp[1], 0.0) or near(mp[1], W):  # walls
                    bndry[f] = 3
                elif mp.distance(center) <= radius:  # cylinder
                    bndry[f] = 5

            # Dump facet markers to file to plot in Paraview
            with XDMFFile('facets.xdmf') as f:
                f.write(bndry)


.. admonition:: Task 3

    Construct the mixed finite element space and the
    bilinear and linear forms together with appropriate
    `DirichletBC <dolfin.fem.bcs.DirichletBC>` object.

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

            bc_walls = DirichletBC(W.sub(0), (0, 0), bndry, 3)

    .. hint::

        To build the forms use::

            # Define trial and test functions
            u, p = TrialFunctions(W)
            v, q = TestFunctions(W)

        Then you can define forms on mixed space using
        ``u``, ``p``, ``v``, ``q`` as usual.


Steady Navier-Stokes flow
-------------------------

.. admonition:: Task 4

    Modify the problem into the Navier-Stokes equations given by

    .. math::
       :label: navierstokes

       - \nu\Delta\mathbf{u} + \mathbf{u}\cdot\nabla\mathbf{u} + \nabla p = 0
            \quad\text{ in }\Omega

    together with :eq:`stokes`:math:`_2`--:eq:`stokes`:math:`_5`.
    Compute the `DFG-flow around cylinder benchmark 2D-1, laminar case, Re=20
    <http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html>`_
    given by :eq:`navierstokes`,
    :eq:`stokes`:math:`_2`--:eq:`stokes`:math:`_5`, :eq:`dfg20`.

    .. hint::

        As usual get rid of `TrialFunctions` in favour of
        nonlinear dependence on `Function`. You can split
        a ``Function`` on a mixed space into components::

            w = Function(W)
            u, p = split(w)

            F = nu*inner(grad(u), grad(v))*dx + ...


.. admonition:: Task 5

    Add computation of lift and drag coefficients :math:`C_\mathrm{D}`,
    :math:`C_\mathrm{L}` and pressure difference :math:`p_\mathrm{diff}`
    as defined on `the DFG 2D-1 website
    <http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html>`_.

    .. hint::

        Use `assemble <dolfin.fem.assembling.assemble>` function
        to evaluate the lift and drag functionals.

        Use either
        `Function.split() <dolfin.functions.function.Function.split>`
        or `Function.sub() <dolfin.functions.function.Function.sub>`
        to extract pressure ``p`` from solution ``w`` for evaluation.
        Evaluate the pressure ``p`` at point ``a = Point(234, 567)``
        by calling ``p(a)``.


.. admonition:: Task 6

    Check computed pressure difference and lift/drag coefficents
    against the reference. Investigate if/how the lift coefficent
    is sensitive to changes in the discretization parameters --
    conduct a convergence study.


Kármán vortex street
--------------------

.. admonition:: Task 7

    Consider evolutionary Navier-Stokes equations

    .. math::

       u_t - \nu\Delta\mathbf{u} + \mathbf{u}\cdot\nabla\mathbf{u} + \nabla p = 0.

    Prepare temporal discretization and timestepping
    to compute the `DFG-flow around cylinder benchmark 2D-1,
    fixed time interval, Re=100
    <http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark3_re100.html>`_ with datum

    .. math::

        U(t) = 1.5 \sin(\tfrac{\pi}{8} t)

    instead of :eq:`dfg20`:math:`_{6b}`. Plot the transient solution.


.. only:: priv

    Reference solution
    ------------------
    .. toggle-header::
        :header: **Show/Hide Code**

        :download:`Download Code <stokes.py>`

        .. literalinclude:: stokes.py
