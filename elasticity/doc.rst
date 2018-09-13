Equation of elasticity
======================

.. sidebar:: Goals

    Further mixed elements, further cool physics...

.. todo::

    Update this sheet!


Deformation of elastic material
-------------------------------

Find approximate solution to following non-linear system of PDEs

.. math::
   :nowrap:

   \begin{align*}
   \vec{u}_t  &= \vec{v}
        &&\quad\text{ in }\Omega\times(0, T), \\
   \vec{v}_t  &= \operatorname{div} (J \mathbb{T} \mathbb{F}^{-\top})
        &&\quad\text{ in }\Omega\times(0, T), \\
   J^2 - 1 &= \begin{cases}
             0         & \text{incompressible case} \\
            -p/\lambda & \text{compressible case}
        \end{cases}
        &&\quad\text{ in }\Omega\times(0, T), \\
   \vec{u} = \vec{v} &= 0
        &&\quad\text{ on }\Gamma_\mathrm{D}\times(0, T), \\
   \mathbb{T}\vec{n} &= \vec{g}
        &&\quad\text{ on }\Gamma_\mathrm{N}\times(0, T), \\
   \mathbb{T}\vec{n} &= 0
        &&\quad\text{ on }\partial\Omega\backslash(\Gamma_\mathrm{D}\cup\Gamma_\mathrm{N})\times(0, T), \\
   \vec{u} = \vec{v} &= 0
        &&\quad\text{ on }\Omega\times\{0\}
   \end{align*}

where

.. math::
   \mathbb{F} &= \mathbb{I} + \nabla\vec{u}, \\
   J &= \det{\mathbb{F}}, \\
   \mathbb{B} &= \mathbb{F}\,\mathbb{F}^\top, \\
   T &= -p\mathbb{I} + \mu (\mathbb{B-I})

using :math:`\theta`-scheme discretization in time and arbitrary discretization
in space with data

.. math::
   \Omega &=\begin{cases}
               [0, 20] \times [0, 1]
               & \text{in 2D} \\
               \text{lego brick } 10 \times 2 \times 1H
               & \text{in 3D}
        \end{cases} \\
   \Gamma_\mathrm{D} &=\begin{cases}
               \left\{ x=0 \right\}
               & \text{in 2D} \\
               \left\{ x = \inf_{\vec{x}\in\Omega}{x} \right\}
               & \text{in 3D}
        \end{cases} \\
   \Gamma_\mathrm{N} &=\begin{cases}
               \left\{ x=20 \right\}
               & \text{in 2D} \\
               \left\{ x = \sup_{\vec{x}\in\Omega}{x} \right\}
               & \text{in 3D}
        \end{cases} \\
   T &= 5, \\
   \vec{g} &=\begin{cases}
             J \mathbb{F}^{-\top}
               \Bigl[\negthinspace\begin{smallmatrix}0\\100t\end{smallmatrix}\Bigr]
               & \text{in 2D} \\
             J \mathbb{F}^{-\top}
               \Bigl[\negthinspace\begin{smallmatrix}0\\0\\100t\end{smallmatrix}\Bigr]
               & \text{in 3D}
        \end{cases} \\
   \mu &= \frac{E}{2(1+\nu)}, \\
   \lambda &=\begin{cases}
             \infty & \text{incompressible case} \\
             \frac{E\nu}{(1+\nu)(1-2\nu)} & \text{compressible case}
        \end{cases} \\
   E &= 10^5, \\
   \nu &=\begin{cases}
             1/2 & \text{incompressible case} \\
             0.3 & \text{compressible case}
        \end{cases}

Mesh file of lego brick :download:`lego_beam.xml`.


.. admonition:: Task 1

   Discretize the equation in time and write variational formulation
   of the problem.


.. admonition:: Task 2

   Build mesh, prepare facet function marking
   :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
   check its correctness.

   .. hint::

       You can get coordinates of :math:`\Gamma_\mathrm{D}` by something like
       ``x0 = mesh.coordinates()[:, 0].min()`` for lego mesh. Analogically
       for :math:`\Gamma_\mathrm{N}`.


.. admonition:: Task 3

   Define Cauchy stress and variational formulation of the problem.

   .. hint::

       Get geometric dimension by ``gdim = mesh.geometry().dim()`` to be able
       to write the code independently of the dimension.


.. admonition:: Task 4

   Prepare a solver and write simple time-stepping loop.

   Prepare a solver by::

      problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J)
      solver = NonlinearVariationalSolver(problem)
      solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
      solver.parameters['newton_solver']['linear_solver'] = 'mumps'

   to increase the tolerance reasonably and employ
   powerful sparse direct solver MUMPS.

   Prepare nice plotting of displacement by::

      plot(u, mode="displacement")

   Manipulate the plot how shown in
   :ref:`the Matplotlib note <unsteady-matplotlib>`.


.. admonition:: Task 5

    Tune the code for getting a 3D solution in a reasonable time.

    Use a following optimization::

        parameters['form_compiler']['quadrature_degree'] = 4

    and P1/P1/P1 spaces.

    You can also try to run the 3D problem in parallel. You can disable
    plotting from commandline by

    .. code-block:: bash

        DOLFIN_NOPLOT=1 mpirun -n 4 python spam_eggs.py


.. only:: solution

    Reference solution
    ------------------

    .. toggle-header::
        :header: **Show/Hide Code**

        .. literalinclude:: elast.py
