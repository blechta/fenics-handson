Time-dependent problems
=======================

Heat equation
-------------

.. sidebar:: Goals

    Learn how to deal with time-dependent problems.
    Solve heat equation by :math:`\theta`-scheme.

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

    \frac{1}{\Delta t} (u^{n+1} - u^n)
    - \theta\Delta u^{n+1} - (1-\theta)\Delta u^n
    &= \theta f(t_{n+1}) + (1-\theta) f(t_n)
        &&\quad\text{ in }\Omega, \; n=0,1,2,\ldots

    \tfrac{\partial u^n}{\partial\mathbf{n}} &= g(t_n)
        &&\quad\text{ on }\partial\Omega, \; n=0,1,2,\ldots

    u^0 &= u_0
        &&\quad\text{ in }\Omega

for a certain sequence :math:`0=t_0 < t_1 < t_2 < ... \leq T`.



.. only:: solution

   Reference solution
   ------------------

   .. toggle-header::
       :header: **Show/Hide Code**

       .. note::

           The reference solution follows `DRY principle
           <https://en.wikipedia.org/wiki/Don%27t_repeat_yourself>`_.
           Hands-on participants are not expected to write such
           a modularized code during a session.

       .. literalinclude:: heat.py
