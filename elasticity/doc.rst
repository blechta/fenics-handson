Hyperelasticity
===============

Find approximate solution to following non-linear system of PDEs

.. math::

    \mathbf{u}_t  &= \mathbf{v}
        &&\quad\text{ in }\Omega\times(0, T),

    \mathbf{v}_t  &= \operatorname{div} (J \mathbb{T} \mathbb{F}^{-\top})
        &&\quad\text{ in }\Omega\times(0, T),

    J^2 - 1 &= \begin{cases}
             0         & \text{incompressible case} \newline
            -p/\lambda & \text{compressible case}
        \end{cases}
        &&\quad\text{ in }\Omega\times(0, T),

    \mathbf{u} = \mathbf{v} &= 0
        &&\quad\text{ on }\Gamma_\mathrm{D}\times(0, T),

    J \mathbb{T} \mathbb{F}^{-\top} \mathbf{n} &= \mathbf{g}
        &&\quad\text{ on }\Gamma_\mathrm{N}\times(0, T),

    J \mathbb{T} \mathbb{F}^{-\top} \mathbf{n} &= 0
        &&\quad\text{ on }\partial\Omega\backslash(\Gamma_\mathrm{D}\cup\Gamma_\mathrm{N})\times(0, T),

    \mathbf{u} = \mathbf{v} &= 0
        &&\quad\text{ on }\Omega\times\{0\}

where

.. math::

    \mathbb{F} &= \mathbb{I} + \nabla\mathbf{u},

    J &= \det{\mathbb{F}},

    \mathbb{B} &= \mathbb{F}\,\mathbb{F}^\top,

    \mathbb{T} &= -p\mathbb{I} + \mu (\mathbb{B-I})

using :math:`\theta`-scheme discretization in time and arbitrary discretization
in space with data

.. math::

    \Omega &=
        \begin{cases}
            (0, 20) \times (0, 1)
            & \text{in 2D} \newline
            \text{lego brick } 10 \times 2 \times 1H
            & \text{in 3D}
        \end{cases}

    \Gamma_\mathrm{D} &=
        \begin{cases}
            \left\{ x=0 \right\}
            & \text{in 2D} \newline
            \left\{ x = \inf_{\mathbf{x}\in\Omega}{x} \right\}
            & \text{in 3D}
        \end{cases}

    \Gamma_\mathrm{N} &=
        \begin{cases}
            \left\{ x=20 \right\}
            & \text{in 2D} \newline
            \left\{ x = \sup_{\mathbf{x}\in\Omega}{x} \right\}
            & \text{in 3D}
        \end{cases}

    T &= 5,

    \mathbf{g} &=
        \begin{cases}
            J \mathbb{F}^{-\top}
            \Bigl[\negthinspace\begin{smallmatrix}0\newline100t\end{smallmatrix}\Bigr]
            & \text{in 2D} \newline
            J \mathbb{F}^{-\top}
            \Bigl[\negthinspace\begin{smallmatrix}0\newline0\newline100t\end{smallmatrix}\Bigr]
            & \text{in 3D}
        \end{cases}

    \mu &= \frac{E}{2(1+\nu)},

    \lambda &=
        \begin{cases}
            \infty & \text{incompressible case} \newline
            \frac{E\nu}{(1+\nu)(1-2\nu)} & \text{compressible case}
        \end{cases}

    E &= 10^5,

    \nu &=
        \begin{cases}
           1/2 & \text{incompressible case} \newline
           0.3 & \text{compressible case}
        \end{cases}

Mesh file of lego brick :download:`lego_beam.xml`.
Within shell download by

.. raw:: html

    <div class="highlight-shell notranslate"><div class="highlight">
    <pre><span></span>wget <script>
        var a = document.createElement('a');
        a.href = "../_downloads/lego_beam.xml";
        document.write(a.href);
    </script></pre></div></div>


.. admonition:: Task 1

    Discretize the equation in time using the Crank-Nicolson
    scheme and derive a variational formulation of the problem.
    Consider discretization using P1/P1/P1 mixed element.


.. admonition:: Task 2

    Build 2D mesh::

        mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), 100, 5, 'crossed')

    Prepare facet function marking
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
    Use time step :math:`\Delta t=\tfrac14`.

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


.. admonition:: Task 4

    Solve the compressible 2D problem.

    Solve the incompressible 2D problem.


.. admonition:: Task 5

    Solve the 3D compressible problem.
    Use time step :math:`\Delta t=\tfrac12`.

    Load mesh by::

        mesh = Mesh('lego_beam.xml')

    Use the following optimization::

        # Limit quadrature degree
        dx = dx(degree=4)
        ds = ds(degree=4)

    You can also try to run the 3D problem in parallel:

    .. code-block:: bash

        # Disable plotting
        export MPLBACKEND=template
        export DOLFIN_NOPLOT=1

        # Run the code on <np> processors
        mpirun -n <np> python <script>.py


.. admonition:: Task 6

    Plot computed displacement :math:`u` in Paraview
    using ``Warp by vector`` filter.


.. only:: pub

    Reference solution
    ------------------

    .. toggle-header::
        :header: **Show/Hide Code**

        :download:`Download Code <elast.py>`

        .. literalinclude:: elast.py
