.. _fenics-intro:

Poisson in a hundred ways
=========================

First touch
-----------

.. toggle-header::
  :header: Login by SSH to ``tyche`` and type: **How to login**

    Open a terminal window by hitting ``CTRL+ALT+T``.
    Use ``ssh`` to connect to a remote system::

        ssh -X -C tyche

    .. note::

        * ``-X`` enable ``X11`` forwarding (allows processes on the remote machine
          opening windows of graphical applications on the local machine)
        * ``-C`` enables compression which is mainly beneficial for access from
          a remote network
        * ``tyche`` stands here for machine ``tyche.mathematik.tu-chemnitz.de``;
          username on the local machine is used by default to login to the remote
          machine; the machine is not accessible from outside the univerity, so
          one would login through a jump host

          .. code-block:: bash

              luigivercotti@local_machine:~$ ssh -X -C user@login.tu-chemnitz.de
              user@login:~$ ssh -X -C tyche
              user@tyche:~$

          Alternatively one can use VPN.

      .. toggle-header::
          :header: **For experts: Kerberos + public key + jump host**

          The most comfortable solution for password-less logins
          outside of the university:

          * add

            .. code-block:: cfg

                Host login
                   Hostname login.tu-chemnitz.de
                   User <user>
                   ForwardAgent yes
                   ForwardX11 yes
                   ForwardX11Trusted yes
                   GSSAPIAuthentication yes
                   GSSAPIDelegateCredentials yes

                Host tyche
                   Hostname tyche.mathematik.tu-chemnitz.de
                   User <user>
                   ForwardAgent yes
                   ForwardX11 yes
                   ForwardX11Trusted yes
                   ProxyCommand ssh login -W %h:%p

            to ``~/.ssh/config``,

          * upload public key to ``login`` and ``tyche``
            by

            .. code-block:: bash

                ssh-copy-id login
                ssh-copy-id tyche

          * `install Kerberos
            <https://www.tu-chemnitz.de/urz/security/ssh.html/>`_

          which allows to get Kerberos ticket by

          .. code-block:: bash

              kinit <user>@TU-CHEMNITZ.DE

          and during its validity password-less login from anywhere

          .. code-block:: bash

              ssh tyche


.. code-block:: bash

    source /LOCAL/opt/fenics-2017.2.0/fenics.conf

to prepare environment for using FEniCS. Now fire up interactive
Python 3 interpreter:

.. code-block:: bash

    python3

You should see something like::

    Python 3.6.5 (default, Apr  1 2018, 05:46:30)
    [GCC 7.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

Now type::

    >>> from dolfin import *
    >>> import matplotlib.pyplot as plt
    >>> 
    >>> mesh = UnitSquareMesh(13, 8)
    >>> plot(mesh)
    [<matplotlib.lines.Line2D object at 0x7fe0003d65c0>, <matplotlib.lines.Line2D object at 0x7fe0003d6748>]
    >>> plt.show()

.. hint::

    Click on ``>>>`` in the right top corner
    of the code snippet to make the code copyable.


A graphical plot of the mesh should appear. If any of the
steps above failed, you're not correctly set up to use FEniCS.
If everything went fine, close the plot window and hit ``^D`` to
quit the interpreter.


Run and modify Poisson demo
---------------------------

.. admonition:: Task 1

    Get the Poisson demo from FEniCS install dir and run it:

    .. code-block:: bash

        mkdir -p work/fenics/poisson
        cd work/fenics/poisson
        cp /LOCAL/opt/fenics-2017.2.0/share/dolfin/demo/documented/poisson/python/demo_poisson.py .
        python3 demo_poisson.py

    You should see some console output and a plot of the solution.


Now login to ``tyche`` from another terminal window and open
the demo file using your favourite editor (if you don't have any
you can use ``gedit``, ``nano``, ...):

.. code-block:: bash

    cd work/fenics/poisson
    <editor> demo_poisson.py


.. admonition:: Task 2

    Now add :ref:`keyword argument <python:tut-keywordargs>`
    ``warp='mode'`` to the `plot <dolfin.common.plotting.plot>` function
    call by applying the following diff:

    .. code-block:: diff

         # Plot solution
         import matplotlib.pyplot as plt
        -plot(u)
        +plot(u, mode='warp')
         plt.show()

    and run the demo again by ``python3 demo_poisson.py``.


.. sidebar:: Hint

    `Constant <dolfin.functions.constant.Constant>`,
    `Expression <dolfin.functions.expression.Expression>`,
    and similar  are clickable links leading to their documentation.

Open :doc:`Poisson demo documentation <demos/poisson/python/demo_poisson.py>`
on the FEniCS website. Notice that the doc page is generated from
the demo file. Go quickly through the docpage while paying attention
to

* definition of weak formulation through forms ``a`` and ``L``,
* usage of `Constant <dolfin.functions.constant.Constant>` and
  `Expression <dolfin.functions.expression.Expression>` classes.


.. admonition:: Task 3

    Modify the code to solve the following problem instead:

    .. math::

        -\Delta u + c u &= f
        &&\text{in } \Omega,

        u &= u_\mathrm{D}
        &&\text{on } \Gamma_\mathrm{D},

        \tfrac{\partial u}{\partial\mathbf{n}} &= g
        &&\text{on } \Gamma_\mathrm{N}

    with

    .. math::
        :nowrap:

        \begin{gather}
            \Omega = (0,1)^2,
            \qquad
            \Gamma_\mathrm{D} = \{(x, y), x=1, 0<y<1\},
            \qquad
            \Gamma_\mathrm{N} = \partial\Omega\setminus\Gamma_\mathrm{D},
        \\
            c = 6,
            \qquad
            f(x, y) = x,
            \qquad
            u_\mathrm{D}(x, y) = y,
            \qquad
            g(x, y) = \sin(5x) \exp(y).
        \end{gather}


Semilinear Poisson equation
---------------------------

.. admonition:: Task 4

    Derive weak formulation for the following semilinear
    Poisson problem:

    .. math::
        :label: nonlinear1

        -\Delta u + u^3 + u &= f
        &&\text{in } \Omega,

        \tfrac{\partial u}{\partial\mathbf{n}} &= g
        &&\text{on } \partial\Omega

    with

    .. math::
        :label: nonlinear2

        \Omega = (0,1)^2,
        \qquad
        f(x, y) = x,
        \qquad
        g(x, y) = \sin(5x) \exp(y).

Notice that the weak formulation has the form

    Find :math:`u\in H^1(\Omega)` such that

.. math::

    F(u; v) = 0
    \qquad \text{for all } v\in H^1(\Omega)

with certain :math:`F` depending on :math:`u` in nonlinear
fashion but being linear in test functions :math:`v`. One
can find the solution iteratively by the Newton method:

    #. Choose :math:`u_0\in H^1(\Omega)`,

    #. For :math:`k=1,2,\ldots` do

        #. Find :math:`\delta u\in H^1(\Omega)` such that

            .. math::
                :label: newton-step

                \frac{\partial F}{\partial u}(u_k; v, \delta u) = -F(u_k; v)
                \qquad \text{for all } v\in H^1(\Omega),


        #. Set :math:`u_{k+1} = u_k + \delta u`.

        #. Check certain convergence criterion and eventually stop iterating.

Here Jacobian :math:`\frac{\partial F}{\partial u}(u; v, \delta u)` is
`Gâteaux derivative <https://en.wikipedia.org/wiki/G%C3%A2teaux_derivative>`_
of :math:`F`. It is generally nonlinear in :math:`u`, but linear in :math:`v`
and :math:`\delta u`. Hence with fixed :math:`u_k\in H^1(\Omega)`
the left-hand side and the right-hand side of :eq:`newton-step`
are a bilinear and linear form respectively and :eq:`newton-step`
is just ordinary linear problem.

.. _fenics-task5:

.. admonition:: Task 5

    Modify the previous code to adapt it to problem
    :eq:`nonlinear1`, :eq:`nonlinear2`.
    Define :math:`F` by filing the gaps in the following code::

        u = Function(V)
        v = TestFunction(V)
        f = Expression(...)
        g = Expression(...)

        F = ...

    If in doubts, peek into :doc:`Nonlinear Poisson demo documentation
    <demos/nonlinear-poisson/python/demo_nonlinear-poisson.py>`.

    Look into documentation of `solve <dolfin.fem.solving.solve>`
    function, read section *Solving nonlinear variational problems*.
    Now you should be able to call the `solve <dolfin.fem.solving.solve>`
    function to obtain the solution.


Nonlinear Dirichlet problem
---------------------------

.. admonition:: Task 6

    Modify the code to solve the following Dirichlet problem:

    .. math::

        -\operatorname{div}(c\nabla u) + 10 u^3 + u &= f
        &&\text{in } \Omega,

        u &= u_\mathrm{D}
        &&\text{on } \partial\Omega

    with

    .. math::

        \Omega = (0,1)^2,
        \qquad
        f(x, y) = 100 x,
        \qquad
        u_\mathrm{D}(x, y) = y,
        \qquad
        c(x, y) = \tfrac{1}{10} + \tfrac12(x^2+y^2).

    .. hint::

        Supply instance of `SubDomain <dolfin.cpp.mesh.SubDomain>`
        class to `DirichletBC <dolfin.fem.bcs.DirichletBC>`.
        How do you tell `SubDomain <dolfin.cpp.mesh.SubDomain>`
        to define :math:`\partial\Omega`? What do you fill in?
        ::

            class Boundary(SubDomain):
                def inside(self, x, on_boundary):
                    return ...

        `on_boundary` argument evaluates to `True` on boundary
        facets, `False` otherwise.


Variational formulation
-----------------------

For :math:`u\in H^1(\Omega)` consider functional

.. math::

    E(u) =
    \int_\Omega \bigl(
        \tfrac12|\nabla u|^2 + \tfrac14u^4 + \tfrac12u^2 - fu
    \bigr) \,\mathrm{d}x
    - \int_{\partial\Omega} gu \,\mathrm{d}s.

Convince yourself that minimization of :math:`F` over :math:`H^1(\Omega)`
is equivalent to problem :eq:`nonlinear1`.

.. admonition:: Task 7

    By filling the following code::

        u = Function(V)
        f = Expression(...)
        g = Expression(...)

        E = ...

    define :math:`E(u)` for data :eq:`nonlinear2`.
    Remember that functionals (zero-forms) do
    not have any test and trial functions.

    Obtain :math:`F(u;v):=\frac{\partial E}{\partial u}(u; v)`
    using `derivative <dolfin.fem.formmanipulations.derivative>`::

        F = derivative(E, u)

    and run the solver like in :ref:`Task 5 <fenics-task5>`.
    Check you get the same solution.


Yet another nonlinearity
------------------------

Consider quasilinear equation in divergence form

.. math::
    :label: nonlinear3

    -\operatorname{div}(\mathcal{A}\nabla u) + u &= f
    &&\text{in } \Omega,

    \tfrac{\partial u}{\partial\mathcal{A}^\top\mathbf{n}} &= 0
    &&\text{on } \partial\Omega,

    \mathcal{A} &= \begin{bmatrix}
        \tfrac{1}{10} + u^2 & 0       \newline
        0                   & 1 + u^2
    \end{bmatrix}
    &&\text{in } \Omega

with data

.. math::
    :label: nonlinear4

    \Omega = (0,1)^2,
    \qquad
    f(x, y) = \tfrac12(x+y).


.. admonition:: Task 8

    Derive weak formulation for the problem :eq:`nonlinear3`.

    Solve the problem :eq:`nonlinear3`, :eq:`nonlinear4`
    using FEniCS. Employ `as_matrix <ufl.tensors.as_matrix>`
    function to define :math:`\mathcal{A}`::

        u = Function(V)
        v = TestFunction(V)

        A = as_matrix((
            (..., ...),
            (..., ...),
        ))
        F = inner(A*grad(u), grad(v))*dx + ...


.. only:: pub

    Reference solution
    ------------------

    .. toggle-header::
        :header: **Show/Hide Code**

        .. literalinclude:: impl.py
