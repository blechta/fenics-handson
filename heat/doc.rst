Heat equation
=============

.. sidebar:: Goals

    Learn how to deal with time-dependent problems.
    Solve heat equation by :math:`\theta`-scheme.
    Plot some nice figures.


We will be interested in solving heat equation:

.. math::

    u_t - \Delta u &= f
        &&\quad\text{ in }\Omega\times(0, T),

    \tfrac{\partial u}{\partial\mathbf{n}} &= g
        &&\quad\text{ on }\partial\Omega\times(0, T),

    u &= u_0
        &&\quad\text{ on }\Omega\times\{0\}

using :math:`\theta`-scheme discretization in time and arbitrary FE discretization
in space with given data :math:`f`, :math:`g`, :math:`u_0`.
:math:`\theta`-scheme time-discrete heat equation reads:

.. math::
    :label: time-discrete

    \frac{1}{\Delta t} \Bigl(u^{n+1} - u^n\Bigr)
    - \theta\Delta u^{n+1} - (1-\theta)\Delta u^n
    &= \theta f(t_{n+1}) + (1-\theta) f(t_n)
        &&\quad\text{ in }\Omega, \; n=0,1,2,\ldots

    \tfrac{\partial u^n}{\partial\mathbf{n}} &= g(t_n)
        &&\quad\text{ on }\partial\Omega, \; n=0,1,2,\ldots

    u^0 &= u_0
        &&\quad\text{ in }\Omega

for a certain sequence :math:`0=t_0 < t_1 < t_2 < ... \leq T`.
Special cases are:

.. list-table::

    * - :math:`\theta=0`
      - explicit Euler scheme,
    * - :math:`\theta=\frac12`
      - Crank-Nicolson scheme,
    * - :math:`\theta=1`
      - implicit Euler scheme.


.. _unsteady-task1:

.. admonition:: Task 1


    Test :eq:`time-discrete` by functions from
    :math:`H^1(\Omega)` and derive a weak formulation
    of :math:`\theta`-scheme for heat equation.


First steps
-----------

Consider data

.. math::
    :label: data0

    \Omega &= (0, 1)^2,

    T &= 2,

    f &= 0,

    g &= 0,

    u_0(x,y) &= x.

.. _unsteady-task2:

.. admonition:: Task 2

    Write FEniCS code implementing problem :eq:`time-discrete`,
    :eq:`data0`, assuming general :math:`\theta`, and arbitrary
    but fixed :math:`\Delta t`. In particular assume::

        from dolfin import *

        mesh = UnitSquareMesh(32, 32)
        V = FunctionSpace(mesh, "Lagrange", 1)

        theta = Constant(0.5)
        dt = Constant(0.1)

    Proceed step-by-step.

        #. **Define all relevant data from** :eq:`data0`.
           Use `Constant <dolfin.functions.constant.Constant>` or
           `Expression <dolfin.functions.expression.Expression>` classes
           to define :math:`f`, :math:`g`, :math:`u_0`.

        #. Define a finite element function for holding
           solution at a particular time step::

               u_n = Function(V)

           and arguments of linear and bilinear forms::

               u, v = TrialFunction(V), TestFunction(V)

        #. **Define bilinear and linear forms describing
           Galerkin descretization of the weak formulation
           derived in**
           :ref:`Task 1 <unsteady-task1>`
           **on the space** ``V``.

           You can conveniently mix bilinear and
           linear terms into a single expression::

               F = 1/dt*(u - u_n)*v*dx + ...

           and separate bilinear and linear part
           using `lhs <ufl.lhs>`, `rhs <ufl.rhs>`::

               a, L = lhs(F), rhs(F)

           .. tip::

               It is good to execute your code every once in a while,
               even when it is not doing anything useful so far,
               e.g., does not have time-stepping yet. You will
               catch the bugs early and fix them easily.

        #. **Prepare for the beggining of time-stepping.**
           Assume ``u0`` is an ``Expression`` or ``Constant``.
           You can use `Function.interpolate()
           <dolfin.cpp.function.Function.interpolate>`
           or `interpolate() <dolfin.fem.interpolation.interpolate>`::

               u_n.interpolate(u0)
               # or
               u_n = interpolate(u0, V)

        #. **Implement time-stepping.** Write a control flow
           statement (for example a :ref:`while <tut-firststeps>` loop) which executes
           the solver for problem ``a == L`` repeatedly while
           updating what needed.

           .. hint::

               Having `Function`\s ``f``, ``g`` on the same space
               you can perform assignment :math:`f := g` by
               ::

                   f.vector()[:] = g.vector()


        #. **Run with different values of**
           :math:`\theta=1,\frac12,0`.

           As a first indicator of correctness of the implementation
           you can drop into the loop lines
           like::

               energy = assemble(u_n*dx)
               print("Energy =", energy)

           Are you observing expected value?


Data IO, plotting
-----------------

There are several possibilities for visualization of data.

.. toggle-header::
    :header: **XDMF output and Paraview**

    One possibility is to use IO facilities of FEniCS and
    visualize using external software, for example Paraview.

    .. Note::

        This approach allows to separate

        * actual computation, which can happen in headless HPC
          environment, for example big parallel clusters of
          thousands of CPU cores,

        * and visualization, which many times needs human
          interaction.

    One can used `XDMFFile <dolfin.cpp.io.XDMFFile>` to store data::

        # Open file for XDMF IO
        f = XDMFFile('solution.xdmf')

        while t < T:

            # Compute time step
            perform_timestep(u_n, t, dt)
            t += dt

            # Save the result to file at time t
            f.write(u_n, t)

    Then you can open Paraview by shell command

    .. code-block:: bash

        paraview &

    and visualize the file ``solution.xdmf``.

.. _unsteady-matplotlib:

.. toggle-header::
    :header: **Matplotlib -- native plotting in Python**


    Another possibility is to use Python plotting library
    `Matplotlib <https://matplotlib.org/>`_.

    .. Note::

        `Matplotlib <https://matplotlib.org/>`_ is Python native
        plotting library, which is programmable and supports

        * interactive use from Python interpreters, including
          popular shells like `Jupyter <https://jupyter.org/>`_,
        * high-quality vector output suitable for scientific
          publishing.

        FEniCS ``plot(obj, **kwargs)`` function implements
        plotting using Matplotlib for several different types
        of ``obj``, for instance ``Function``, ``Expression``,
        ``Mesh``, ``MeshFunction``. As Matplotlib is highly
        programmable and customizable, FEniCS ``plot()`` is
        typically accompanied by some native matplotlib
        commands. Mimimal example of
        interaction of FEniCS and matplotlib::

            from dolfin import *
            import matplotlib.pyplot as plt

            mesh = UnitSquareMesh(64, 64)
            plot(mesh)
            plt.savefig('mesh_64_64.pdf')  # Render to PDF
            plt.show()  # Render into interactive window

    Add something along the lines of::

        import matplotlib.pyplot as plt

        # Open a plot window
        fig = plt.figure()
        fig.show()

        while t < T:

            # Compute time step
            perform_timestep(u_n, t, dt)
            t += dt

            # Update plot to current time step
            fig.clear()
            p = plot(u_n, mode="warp")
            fig.colorbar(p)
            fig.gca().set_zlim((0, 2))
            fig.canvas.draw()

    .. warning::

        Matplotlib's interactive capabalities aparently
        depend on used `Matplotlib backend
        <https://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`_.
        In particular updating the contents of the plot window seems
        to work fine with ``TkAgg`` backend. Issue shell command

        .. code-block:: shell

            export MPLBACKEND=tkagg

        to choose ``TkAgg`` in the current shell session.


.. admonition:: Task 3

    Implement at least one of the aforementioned ways to
    plot your solutions in time. Check that your solution
    of :ref:`Task 2 <unsteady-task2>` looks reasonable.


Nonhomogeneous Neumann BC
-------------------------

Consider :eq:`time-discrete`, :eq:`data0` but now with
nonhomogeneous Neumann data

.. math::
    :label: data1

    g &= 1 \text{ on } \{ x = 0 \},

    g &= 0 \text{ elsewhere}.


.. admonition:: Task 3

    #. Derive weak formulation describing
       :eq:`time-discrete`, :eq:`data0`, :eq:`data1`.

    #. Define surface measure supported on the left
       boundary of the unit square mesh by following
       steps:

        #. subclass `SubDomain <dolfin.cpp.mesh.SubDomain>`,
        #. define `MeshFunction <dolfin.cpp.mesh.MeshFunction>`,
        #. mark the mesh function using
           `SubDomain.mark <dolfin.cpp.mesh.SubDomain.mark>` method,
        #. define integration `Measure <ufl.Measure>`.

    .. hint::
        .. toggle-header::
            :header: **Show/Hide Code**

            ::

                # Define instance of SubDomain class
                class Left(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[0], 0)
                left = Left()

                # Define and mark mesh function on facets
                facets = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
                left.mark(facets, 1)

                # Define exterior facet measure where facets==1
                ds_left = Measure("ds", mesh, subdomain_data=facets, subdomain_id=1)

    3. Using the surface measure, modify the implementation
       from :ref:`Task 2 <unsteady-task2>` to incorporate boundary
       condition :eq:`data1`.

    #. Run the code with :math:`\theta=1` and check that the
       results look as expected.


Time-dependent BC
-----------------

Consider time-dependent data

.. math::
    :label: data2

    f(x, t) &= 2-t,

    g(x, t) &= \left\{\begin{array}{ll}
                   t & x = 0, \newline
                   0 & \text{otherwise}.
               \end{array}\right.

.. admonition:: Task 4

    Modify solution of the previous task to use data :eq:`data2`.

    .. hint::

        You can use `Constant.assign()
        <dolfin.cpp.function.Constant.assign>` or
        ``Expression.<param> = <value>``
        to change existing `Constant` or `Expression`.
        Look for *User defined parameters* in `Expression
        <dolfin.functions.expression.Expression>`
        documentation.


Now consider different time-dependent data

.. math::
    :label: data3

    f(x, t) &= 0,

    g(x, t) &=
        \left\{\begin{array}{ll}
            \max\bigl(0, \tfrac{1-t}{2}\bigr) & x = 0, \newline
            0                                 & \text{otherwise}.
        \end{array}\right.

.. admonition:: Task 5

    Modify solution of the previous task to use data :eq:`data3`.


Adaptive time-stepping
----------------------

Consider solution of *low* precision generated by timestep
:math:`\Delta t`:

.. math::
    :label: steplo

    \frac{1}{\Delta t} \Bigl(u^{n+1}_\mathrm{low} - u^n\Bigr)
    - \theta\Delta u^{n+1}_\mathrm{low} - (1-\theta)\Delta u^n
    = \theta f(t_{n+1}) + (1-\theta) f(t_n)

and solution of *high* precision computed by two timesteps
of a half size:

.. math::
    :label: stephi

    \frac{1}{\Delta t/2} \Bigl(u^{n+1/2}_\mathrm{high} - u^n\Bigr)
    - \theta\Delta u^{n+1/2}_\mathrm{high} - (1-\theta)\Delta u^n
    &= \theta f(t_{n+1/2}) + (1-\theta) f(t_n),

    \frac{1}{\Delta t/2} \Bigl(u^{n+1}_\mathrm{high} - u^{n+1/2}_\mathrm{high}\Bigr)
    - \theta\Delta u^{n+1}_\mathrm{high} - (1-\theta)\Delta u^{n+1/2}_\mathrm{high}
    &= \theta f(t_{n+1}) + (1-\theta) f(t_{n+1/2}).

By `Richardson extrapolation
<https://en.wikipedia.org/wiki/Richardson_extrapolation>`_
one can estimate the error of discretization (in time)
by quantity:

.. math::
    :label: estimator

    \eta :=
    \frac{\|u^{n+1}_\mathrm{high} - u^{n+1}_\mathrm{low}\|_{L^2(\Omega)}}
         {2^p - 1}

where

.. math::
    :label: order

    p = \left\{\begin{array}{ll}
        2 && \theta=\tfrac12, \newline
        1 && \text{otherwise}
        \end{array}\right.

is a theoretical order of accuracy of the :math:`\theta`-scheme.
Given a tolerance :math:`\mathrm{Tol}` set the new timestep to

.. math::
    :label: dtnew

    \Delta t^* :=
    \left( \frac{\rho\,\mathrm{Tol}}{\eta} \right)^\frac1p \Delta t.

Here :math:`0<\rho\leq1` is a chosen safety factor. That
asymptotically ensures that the error (or at least the
estimator) committed with the new time step is
:math:`\rho`-multiple of the tolerance.

Now consider an algorithm:

    #. compute :math:`u^{n+1}_\mathrm{low}`
       and :math:`u^{n+1}_\mathrm{high}`
    #. compute :math:`\eta`
    #. compute :math:`\Delta t^*`
    #. | if :math:`\eta\leq\mathrm{Tol}`:
       |     :math:`u^{n+1}:=u^{n+1}_\mathrm{high}`
       |     :math:`n \mathrel{+}= 1`
    #. update timestep :math:`\Delta t:=\Delta t^*`

.. admonition:: Task 6

    Solve :eq:`time-discrete`, :eq:`data0`:math:`_{1,2,5}`,
    :eq:`data3` using the adaptive strategy described above.


.. only:: priv

    Reference solution
    ------------------
    .. toggle-header::
        :header: **Show/Hide Code**

        .. note::

            The reference solution follows `DRY principle
            <https://en.wikipedia.org/wiki/Don%27t_repeat_yourself>`_.
            Hands-on participants are not expected to write such
            a modularized code during a session.

            More clean design would be achieved by employing classes.
            It is in general a good idea to start just with free
            functions and refactor later into classes when developing
            an object oriented code.

        .. literalinclude:: heat.py
