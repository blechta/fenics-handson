Stokes equation
===============

Stokes flow around cylinder
---------------------------

Solve the following linear system of PDEs

.. math::
   - \operatorname{div}(\nabla u) + \nabla p &= f
        \quad\text{ in }\Omega, \\
   \operatorname{div} u &= 0
        \quad\text{ in }\Omega, \\
   u &= 0
        \quad\text{ on }\Gamma_\mathrm{D}, \\
   u &= u_\mathrm{IN}
        \quad\text{ on }\Gamma_\mathrm{IN}, \\
   \tfrac{\partial u}{\partial\mathbf{n}} &= g
        \quad\text{ on }\Gamma_\mathrm{N}, \\

using FE discretization with data

* :math:`\Omega = [0, 2.2]x[0, 0.41] - B_{0.05}\left([0.2,0.2]\right)`
* :math:`\Gamma_\mathrm{N} = \left\{ x = 2.2 \right\}`
* :math:`\Gamma_\mathrm{IN} = \left\{ x = 0.0 \right\}` 
* :math:`\Gamma_\mathrm{D} = \Gamma_\mathrm{W} \cup \Gamma_\mathrm{S}` 
* :math:`u_\mathrm{IN} = \left( 0.3 \frac{4}{0.41^2} y (0.41-y) , 0 \right)`

where :math:`B_R(\mathbf{z})` is a circle of radius :math:`R` and center
:math:`\mathbf{z}` 

  .. image:: geometry.png
     :align: center
     :width: 70%
     :target: http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html
      
..

  **Task 1.** Write the variational formulation of the problem and
  discretize the equation by mixed finite element method.

  **Task 2.** Build mesh, prepare facet function marking
  :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
  check its correctness.

      *Hint.* Use the mshr component of fenics - see `mshr documentation
      <https://bitbucket.org/benjamik/mshr/wiki/API>`_
      
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


      *Hint.* Try yet another way to mark the boundaries by direct
      access to the mesh entities by ``facets(mehs)``,
      ``vertices(mesh)``, ``cells(mesh)``
      
      .. code-block:: python

         # Construct facet markers
         bndry = FacetFunction("size_t", mesh)
         for f in facets(mesh):
              mp = f.midpoint()
              if near(mp[0], 0.0): bndry[f] = 1  # inflow
              elif near(mp[0], L): bndry[f] = 2  # outflow
              elif near(mp[1], 0.0) or near(mp[1], W): bndry[f] = 3  # walls
              elif mp.distance(center) <= radius:      bndry[f] = 5  # cylinder
         
         plot(boundary_parts, interactive=True)

                                                          
  **Task 3.** Construct the mixed finite element space and the
  bilinear and linear forms together with the ``DirichletBC`` object.

      *Hint.* Use for example the stable Taylor-Hood finite elements.
      
      .. code-block:: python

         # Build function spaces (Taylor-Hood)
         V = VectorFunctionSpace(mesh, "Lagrange", 2)
         P = FunctionSpace(mesh, "Lagrange", 1)
         W = MixedFunctionSpace([V, P])

      *Hint.* To define Dirichlet BC on subspace use the ``W.sub`` method.

      .. code-block:: python
                      
         noslip = Constant((0, 0))
         bc_walls = DirichletBC( W.sub(0) , noslip , bndry , 3 )

      *Hint.* To build the forms use the ``split`` method or function
      to access the components of the mixed space.

      .. code-block:: python
         
         # Define unknown and test function(s)
         (v_, p_) = TestFunctions(W)
         (v , p)  = TrialFunctions(W)


      Then you can define the forms for example as:

      .. code-block:: python

          def a(u,v): return inner(grad(u), grad(v))*dx
          def b(p,v): return p*div(v)*dx
          def L(v):   return inner(f, v)*dx

          F = a(v,v_) + b(p,v_) + b(p_,v) - L(v_)


      And solve by:

      .. code-block:: python

          w = Function(W)
          solve(lhs(F)==rhs(F), w, bcs)
          (v,p)=w.split(w)
            
                          
  **Task 4.** Now modify the problem to the Navier-Stokes equations
  and compute the `DFG-flow around cylinder benchmark
  <http://www.featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html>`_

    *Hint.* You can use generic ``solve`` function or
    ``NonlinearVariationalProblem`` and ``NonlinearVariationalSolver``
    classes. 

    .. code-block:: python

       (_v, _p) = TestFunctions(W)
       w = Function(W)
       (v, p) = split(w)

       I = Identity(v.geometric_dimension())    # Identity tensor

       D = 0.5*(grad(v)+grad(v).T)  # or D=sym(grad(v))
       T = -p*I + 2*nu*D

       # Define variational forms
       F = inner(T, grad(_v))*dx + _p*div(v)*dx + inner(grad(v)*v,_v)*dx

       
    
    *Hint.* Use ``Assemble`` function to evaluate the lift and drag functionals.

    .. code-block:: python

       (v, p) = w.split(True)

       force = T*n
       D=(force[0]/0.002)*ds(5)
       L=(force[1]/0.002)*ds(5)

       drag = assemble(D)
       lift = assemble(L)

       print "drag= %e    lift= %e" % (drag , lift)

                        
.. only:: solution

   Reference solution - Stokes problem
   -----------------------------------

   .. literalinclude:: stokes.py
      :start-after: # Begin code


   Reference solution - benchmark problem
   --------------------------------------

   .. literalinclude:: bench_ns.py
      :start-after: # Begin code
