Equation of elasticity
======================

Deformation of elastic material
-------------------------------

Find approximate solution to following non-linear system of PDEs

.. math::
   \vec{u}_t  &= \vec{v}
        \quad\text{ in }\Omega\times[0, T], \\
   \vec{v}_t  &= \operatorname{div} (J \mathbb{T} \mathbb{F}^{-\top})
        \quad\text{ in }\Omega\times[0, T], \\
   J^2 - 1 &= \left\{\begin{array}{ll}
             0         && \text{incompressible case} \\
            -p/\lambda && \text{compressible case}
        \end{array}\right.
        \quad\text{ in }\Omega\times[0, T], \\
   \vec{u} = \vec{v} &= 0
        \quad\text{ on }\Gamma_\mathrm{D}\times[0, T], \\
   \mathbb{T}\vec{n} &= \vec{g}
        \quad\text{ on }\Gamma_\mathrm{N}\times[0, T], \\
   \mathbb{T}\vec{n} &= 0
        \quad\text{ on }\Omega\backslash\Gamma_\mathrm{D}\cup\Gamma_\mathrm{N}\times[0, T], \\
   \vec{u} = \vec{v} &= 0
        \quad\text{ on }\Omega\times{0} \\

where

.. math::
   \mathbb{F} &= \mathbb{I} + \nabla\vec{u}, \\
   J &= \det{\mathbb{F}}, \\
   \mathbb{B} &= \mathbb{F}\,\mathbb{F}^\top, \\
   T &= -p\mathbb{I} + \mu (\mathbb{B-I})

using :math:`\theta`-scheme discretization in time and arbitrary discretization
in space with data

.. math::
   \Omega &= \left\{\begin{array}{ll}
               [0, 20] \times [0, 1]
               && \text{in 2D} \\
               \text{lego brick } 10 \times 2 \times 1H
               && \text{in 3D}
        \end{array}\right. \\
   \Gamma_\mathrm{D} &= \left\{\begin{array}{ll}
               \left\{ x=0 \right\}
               && \text{in 2D} \\
               \left\{ x = \inf_{\vec{x}\in\Omega}{x} \right\}
               && \text{in 3D}
        \end{array}\right. \\
   \Gamma_\mathrm{N} &= \left\{\begin{array}{ll}
               \left\{ x=20 \right\}
               && \text{in 2D} \\
               \left\{ x = \sup_{\vec{x}\in\Omega}{x} \right\}
               && \text{in 3D}
        \end{array}\right. \\
   T &= 5, \\
   \vec{g} &= \left\{\begin{array}{ll}
             J \mathbb{F}^{-\top}
               \Bigl[\negthinspace\begin{smallmatrix}0\\100t\end{smallmatrix}\Bigr]
               && \text{in 2D} \\
             J \mathbb{F}^{-\top}
               \Bigl[\negthinspace\begin{smallmatrix}0\\0\\100t\end{smallmatrix}\Bigr]
               && \text{in 3D}
        \end{array}\right. \\
   \mu &= \frac{E}{2(1+\nu)}, \\
   \lambda &= \left\{\begin{array}{ll}
             \infty && \text{incompressible case} \\
             \frac{E\nu}{(1+\nu)(1-2\nu)} && \text{compressible case}
        \end{array}\right. \\
   E &= 10^5, \\
   \nu &= \left\{\begin{array}{ll}
             1/2 && \text{incompressible case} \\
             0.3 && \text{compressible case}
        \end{array}\right.

Mesh file of lego brick :download:`lego_beam.xml`.


..

   **Task 1.** Discretize the equation in time and write variational formulation
   of the problem.

   **Task 2.** Build mesh, prepare facet function marking
   :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
   check its correctness.

     *Hint.*
     You can get coordinates of :math:`\Gamma_\mathrm{D}` by something like
     ``x0 = mesh.coordinates()[:, 0].min()`` for lego mesh. Analogically
     for :math:`\Gamma_\mathrm{N}`.

   **Task 3.** Define Cauchy stress and variational formulation of the problem.

     *Hint.*
     Get geometric dimension by ``gdim = mesh.geometry().dim()`` to be able
     to write the code independently of the dimension.

   **Task 4.** Prepare a solver and
   write simple time-stepping loop.

     Prepare a solver by

     .. code-block:: python

        problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J)
        solver = NonlinearVariationalSolver(problem)
        solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
        solver.parameters['newton_solver']['linear_solver'] = 'mumps'

     to increase the tolerance reasonably and employ powerful LU solver MUMPS.

     Prepare nice plotting of displacement by

     .. code-block:: python

        plt = plot(u, mode="displacement", interactive=False, wireframe=True)

     and then just update a plot by ``plt.plot(u)`` every time-step.

   **Task 5.** Tune the code for getting a 3D solution in a reasonable time.

     Use a following optimization

     .. code-block:: python

        parameters['form_compiler']['representation'] = 'uflacs'
        parameters['form_compiler']['optimize'] = True
        parameters['form_compiler']['quadrature_degree'] = 4

     and P1/P1/P1 spaces.

     You can also try to run the 3D problem in parallel. You can disable
     plotting from commandline by

     .. code-block:: bash

        DOLFIN_NOPLOT=1 mpirun -n 4 python spam_eggs.py


.. only:: solution

   Reference solution
   ------------------

   .. literalinclude:: elast.py
      :start-after: # Begin code
