p-Laplace equation
==================

.. todo::

    Update sheet and test reference solution!


Potential for Laplace equation
------------------------------

.. admonition:: Task 1

    Formulate Laplace equation

    .. math::
        -\Delta u &= f &&\text{ in } \Omega,

                u &= 0 &&\text{ on } \partial\Omega,

    as a variational problem (i.e., find potential for the equation) and solve it
    using FEniCS with data

    * :math:`\Omega = [0, 1] \times [0, 1]`,
    * :math:`f = 1 + \cos(2 \pi x) \sin(2 \pi y)`.

    Use `UFL function derivative
    <https://fenics.readthedocs.io/projects/ufl/en/stable/api-doc/ufl.html#ufl.formoperators.derivative>`_
    to automatically obtain Galerkin formulation
    from the potential. Don't assume linearity of the PDE - solve it as nonliner
    (Newton will converge in 1 step).


Potential for p-Laplace equation
--------------------------------

.. admonition:: Task 2

    Replace every occurrence of number :math:`2` in potential for
    Laplace equation by :math:`p`. This is called :math:`p`-Laplacian for
    :math:`1 < p < +\infty`.

    Convince yourself that resulting PDE is non-linear whenever :math:`p \neq 2`.


.. admonition:: Task 3

    Run the algorithm from Task 1 with :math:`p = 1.1` and
    :math:`p = 11`.

    .. hint::

        Use DOLFIN class ``Constant`` to avoid form recompilation by FFC for
        distinct values of :math:`p`.

   Do you know where is the problem? If not, compute second Gateaux derivative
   of the potential (which serves as Jacobian for the Newton-Rhapson algorithm)
   and look at its value for :math:`u = 0`.


.. admonition:: Task 4

    Add some regularization to the potential to make it
    non-singular/non-degenerate. (Prefer regularization without square roots.)
    Find a solution of the original :math:`p`-Laplace problem with
    :math:`p = 1.1` and :math:`p = 11` using careful approximation by
    regularized problem. Report :math:`max_\Omega u_h` for :math:`u_h` being the
    approximate solution.

    .. hint::

        For ``u_h`` Lagrange degree 1 ``Function`` the maximum matches
        maximal nodal value so it is ``u_h.vector().max()`` because nodal basis
        is used. *Warning.* This does not generally hold for higher polynomial
        degrees!


.. only:: priv

   .. note::

      *Lecturer note.* Totally wrong solution may be easily obtained by
      overregularization. Student must find a (heuristic) way to show that
      convergence is established.

      The trick in regularization-convergence for :math:`p = 11` is that Newton
      solver needs to be initialized by previous result computed using higher
      regularization - see ``for``-loop in ``p_large.py``.

      Convergence in discretization parameter is not so important here compared
      to regularization and good solution is obtained on mild meshes. Also exact,
      appropriately high quadrature for large p is not needed here -
      ``parameters['form_compiler']['quadrature_degree'] = 11`` yields the same
      results.

      Solution for :math:`p = 1.1,\,11`  has maximum around ``1.3E-7`` and
      ``0.41`` respectively.

.. only:: priv

    Reference solution
    ------------------

    Reference solution consists of three files -- one module file:

        .. toggle-header::
            :header: **Show/Hide File** ``p_laplace.py``

            .. literalinclude:: p_laplace.py

    and two executable scripts:

        .. toggle-header::
            :header: **Show/Hide File** ``p_small.py``

            .. literalinclude:: p_small.py

        .. toggle-header::
            :header: **Show/Hide File** ``p_large.py``

            .. literalinclude:: p_large.py
