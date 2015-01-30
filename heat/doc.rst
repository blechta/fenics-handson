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
* :math:`f(x) = \left\{ \begin{array}{ll} 1 && \text{ if }\operatorname{dist}(x, \frac{3}{4}) < \frac{1}{5} \\
                                          0 && \text{ otherwise } \end{array} \right.`
* :math:`u_0 = \ldots`

..

   **Task 1.** Discretize the equation in time and write variational formulation
   of the problem.
