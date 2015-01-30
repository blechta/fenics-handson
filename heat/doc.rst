Heat equation
=============

Heat equation in moving media
-----------------------------

Find approximate solution to following linear PDE

.. math::
   u_t + \mathbf{b}\cdot\nabla{u} - \operatorname{div}(K \nabla u) &= f
        \quad\text{ in }\Omega\times[0, T], \\
   u &= u_\mathrm{D}
        \quad\text{ in }\Omega_\mathrm{D}\times[0, T], \\
   \tfrac{\partial u}{\partial\mathbf{n}} &= g
        \quad\text{ on }\Gamma_\mathrm{N}\times[0, T], \\
   u &= u_0
        \quad\text{ on }\Omega\times{0} \\

using :math:`\theta`-scheme discretization in time and arbitrary discretization
in space with data

* :math:`\Omega = [0, 1]^2`
* :math:`T = 10`
* :math:`\Gamma_\mathrm{N} = \left\{ x = 0 \right\}`
* :math:`\Gamma_\mathrm{D} = \left\{ x = 1 \right\} \cup \left\{ y = 0 \right\}`
* :math:`g = 0.1`
* :math:`K = 0.01`
* :math:`\mathbf{b} = \left( -y-\tfrac{1}{2}, x-\tfrac{1}{2} \right)`
* :math:`f = \chi_{ B_{1/5}\left(\left[\frac{3}{4}, \frac{3}{4}\right]\right) }`
* :math:`u_0(\mathbf{x}) = \left( 1 - 25
  \operatorname{dist}\left(\mathbf{x}, \left[\frac{1}{4}, \frac{1}{4}\right]\right)
  \right)
  \chi_{ B_{1/5}\left(\left[\frac{1}{4}, \frac{1}{4}\right]\right) }`

where :math:`\chi_X` is a characteristic function of set :math:`X`,
:math:`B_R(\mathbf{z})` is a ball of radius :math:`R` and center
:math:`\mathbf{z}` and :math:`\operatorname{dist}(\mathbf{p}, \mathbf{q})`
is Euclidian distance between points :math:`\mathbf{p}`, :math:`\mathbf{q}`.

..

   **Task 1.** Discretize the equation in time and write variational formulation
   of the problem.

   **Task 2.** Build mesh, prepare facet function marking
   :math:`\Gamma_\mathrm{N}` and :math:`\Gamma_\mathrm{D}` and plot it to
   check its correctness.

      *Hint.* You can follow the procedure from `subdomains-poisson demo
      <http://fenicsproject.org/documentation/dolfin/1.5.0/python/demo/
      documented/subdomains-poisson/python/documentation.html#implementation>`_.
      (Follow a construction of ``boundaries`` object therein.)

   **Task 3.** Define expressions :math:`\mathbf{b}`, :math:`f`, :math:`u_0`
   and plot them.

     *Hint.*
     According to your personal taste, get hint at `Expression class documentation
     <http://fenicsproject.org/documentation/dolfin/1.5.0/python/
     programmers-reference/functions/expression/Expression.html>`_ or any
     `documented demo <http://fenicsproject.org/documentation/dolfin/1.5.0/
     python/demo/index.html>`_. Use any kind of expression you wish (subclassing
     Python ``Expression``, oneline C++, subclassing C++ ``Expression``).

  **Task 4.** Use facet markers from Task 2 to define ``DirichletBC`` object
  and ``Measure`` for integration along :math:`\Gamma_\mathrm{N}`.

  **Task 5.** Now proceed to variational formulation and time-stepping loop.
  Write bilinear and linear form representing PDE. How is solution at previous
  time-step represented therein?

    *Hint.* Use ``LinearVariationalProblem`` and ``LinearVariationalSolver``
    classes so that ``solve`` method of an instance of the latter is called
    every time-step while nothing else is touched excepted updating value
    of solution from previous time-step figuring in variational form. You
    can use for instance ``Function.assign`` method to do that.
