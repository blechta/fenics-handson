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
   u &= u_\mathrm{D}
        \quad\text{ in }\Omega_\mathrm{D}, \\
   \tfrac{\partial u}{\partial\mathbf{n}} &= g
        \quad\text{ on }\Gamma_\mathrm{N}, \\

using FE discretization with data

* :math:`\Omega = [0, 2.2]x[0, 0.41] - B_{0.05}\left([0.2,0.2]\right)`
* :math:`\Gamma_\mathrm{N} = \left\{ x = 2.2 \right\}`
* :math:`\Gamma_\mathrm{D} = ` otherwise

where :math:`B_R(\mathbf{z})` is a ball of radius :math:`R` and center
:math:`\mathbf{z}` 

..

  **Task 1.** Discretize the equation in time and write variational formulation
  of the problem.

  **Task 2.** Build mesh, prepare facet function marking
  :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
  check its correctness.

      *Hint.* Use the mshr component of fenics `subdomains-poisson demo
      <http://fenicsproject.org/documentation/dolfin/1.5.0/python/demo/
      documented/subdomains-poisson/python/documentation.html#implementation>`_.
      (Follow a construction of ``boundaries`` object therein.)

      .. code-block:: python

         mesh = UnitSquareMesh(10, 10, 'crossed')

         # Create boundary markers
         boundary_parts = FacetFunction('size_t', mesh)
         left   = AutoSubDomain(lambda x: near(x[0], 0.0))
         right  = AutoSubDomain(lambda x: near(x[0], 1.0))
         bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
         left  .mark(boundary_parts, 1)
         right .mark(boundary_parts, 2)
         bottom.mark(boundary_parts, 2)
         plot(boundary_parts, interactive=True)

                      
  **Task 3.** Define expressions :math:`\mathbf{b}`, :math:`f`, :math:`u_0`
  and plot them.

        *Hint.*
        According to your personal taste, get hint at `Expression class documentation
        <http://fenicsproject.org/documentation/dolfin/1.5.0/python/
        programmers-reference/functions/expression/Expression.html>`_ or any
        `documented demo <http://fenicsproject.org/documentation/dolfin/1.5.0/
        python/demo/index.html>`_. Use any kind of expression you wish (subclassing
        Python ``Expression``, oneline C++, subclassing C++ ``Expression``).

        .. code-block:: python

           # python Expression subclass
           class b_Expression(Expression):
               def eval(self, value, x):
                   vx = x[0] - 0.5
                   vy = x[1] - 0.5
                   value[0] = -vy
                   value[1] =  vx
               def value_shape(self):
                   return (2,)

           b = b_Expression()

        .. code-block:: python

           # oneline C++
           b = Expression(("-(x[1] - 0.5)", "x[0] - 0.5"))
           
                                                        
  **Task 4.** Use facet markers from Task 2 to define ``DirichletBC`` object
  and ``Measure`` for integration along :math:`\Gamma_\mathrm{N}`.

     *Hint* See `UFL class Measure <http://fenicsproject.org/documentation/ufl/1.5.0/ufl.html#ufl.classes.Measure>`_ 

     .. code-block:: python

        dsN = Measure("ds", subdomain_id=1, subdomain_data=boundary_parts)
             
                          
  **Task 5.** Now proceed to variational formulation and time-stepping loop.
  Write bilinear and linear form representing PDE. How is solution at previous
  time-step represented therein?

    *Hint.* Use ``LinearVariationalProblem`` and ``LinearVariationalSolver``
    classes so that ``solve`` method of an instance of the latter is called
    every time-step while nothing else is touched excepted updating value
    of solution from previous time-step figuring in variational form. You
    can use for instance ``Function.assign`` method to do that.


  **Task 6.** Add solution output for external visualisation, like
  Paraview.

     *Hint* See `Poisson demo <http://fenicsproject.org/documentation/dolfin/1.5.0/python/demo/documented/poisson/python/documentation.html#index-0>`_

     .. code-block:: python
                             
        # Create file for storing results
        f = File("results/u.xdmf")

        u.rename("u", "temperature")
        f << u

                        
.. only:: solution

   Reference solution
   ------------------

   .. literalinclude:: impl.py
      :start-after: # Begin code


.. only:: solution

   Reference solution
   ------------------

   .. literalinclude:: stokes.py
      :start-after: # Begin code
