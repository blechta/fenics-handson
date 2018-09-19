Eigenfunctions of Laplacian and Helmholtz equation
==================================================

.. todo::

    Fix sheet and test reference solution!


Wave equation with harmonic forcing
-----------------------------------

Let's have wave equation with special right-hand side

.. math::

   w_{tt} - \Delta w &= f\, e^{i\omega t}
       &&\quad\text{ in }\Omega\times(0,T),

                   w &= 0
       &&\quad\text{ on }\partial\Omega\times(0,T),

with :math:`f \in L^2(\Omega)`. Such a problem has a solution (in some proper
sense; being unique when enriched by initial conditions), see [Evans]_,
chapter 7.2.

.. admonition:: Task 1

    Try seeking for a particular solution of this equation while
    taking advantage of special structure of right-hand side. Assuming ansatz

    .. math::
        w := u\, e^{i\omega t}, \quad u\in H_0^1(\Omega)

    derive non-homogeneous Helmholtz equation for :math:`u` using the Fourier
    method and try solving it using FEniCS with

    * :math:`\Omega = [0,1]\times[0,1]`,
    * :math:`\omega = \sqrt{5}\pi`,
    * :math:`f = x + y`.


.. admonition:: Task 2

    Plot solution energies against number of degrees of freedom.

    .. hint::

        Having list of number of degrees of freedom ``ndofs`` and list of
        energies ``energies`` do::

           import matplotlib.pyplot as plt
           plt.plot(ndofs, energies, 'o-')
           plt.show()

    What does it mean? Is the problem well-posed?


Helmholtz equation and eigenspaces of Laplacian
-----------------------------------------------

Define eigenspace of Laplacian (with zero BC) corresponding to :math:`\omega^2`

.. math::

   E_{\omega^2} = \{ u\in H_0^1(\Omega): -\Delta u = \omega^2 u \}.

If :math:`E_{\omega^2}\neq{0}` then :math:`\omega^2` is eigenvalue. Then by
testing the non-homogeneous Helmholtz equation (derived in previous section) by
non-trivial :math:`v\in E_{\omega^2}` one can see that
:math:`f\perp E_{\omega^2}` is required (**check it!**), otherwise the problem
is ill-posed. Hence the assumed ansatz is generally wrong. In fact, the
condition :math:`f\perp E_{\omega^2}` is sufficient condition for well-posedness
of the problem, see [Evans]_, chapter 6.2.3.

The resolution is to seek for a particular solution for :math:`f^\parallel` and
:math:`f^\perp` (:math:`L^2`-projections of :math:`f` to :math:`E_{\omega^2}`
and :math:`^\perp E_{\omega^2}` respectively) separately. As :math:`E_{\omega^2}`
has finite dimension (due to the Fredholm theory), the former can be obtained by
solving forced harmonic oscillator equation which is easily done in hand (once
:math:`E_{\omega^2}` is known). This is equivalent to picking blowing-up ansatz
:math:`w = u\, t\, e^{i t\omega},\, u\in H_0^1(\Omega)`. So let's focus to the
latter part.

.. todo::

   Make the exposition above simpler.


.. admonition:: Task 3

    Use SLEPc eigensolver to find :math:`E_{\omega^2}`.

    .. hint::

        Having assembled matrices ``A``, ``B``, the eigenvectors solving

        .. math::

            A x = \lambda B x

        with :math:`\lambda` close to target ``lambd`` can be found by::

            eigensolver = SLEPcEigenSolver(as_backend_type(A), as_backend_type(B))
            eigensolver.parameters['problem_type'] = 'gen_hermitian'
            eigensolver.parameters['spectrum'] = 'target real'
            eigensolver.parameters['spectral_shift'] = lambd
            eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
            eigensolver.parameters['tolerance'] = 1e-6
            #eigensolver.parameters['verbose'] = True  # for debugging
            eigensolver.solve(number_of_requested_eigenpairs)

            eig = Function(V)
            eig_vec = eig.vector()
            space = []
            for j in range(eigensolver.get_number_converged()):
                r, c, rx, cx = eigensolver.get_eigenpair(j)
                eig_vec[:] = rx
                plot(eig, title='Eigenvector to eigenvalue %g'%r)
                plt.show()


.. admonition:: Task 4

    Write function which takes a tuple of functions and
    :math:`L^2`-orthogonalizes them using Gramm-Schmidt algorithm.


.. admonition:: Task 5

    Compute :math:`f^\perp` for :math:`f` from Task 1 and solve the
    Helmholtz equation with :math:`f^\perp` on right-hand side. Again, plot
    energies of solutions against number of degrees of freedom.


    .. only:: priv

        .. note::

            *Lecturer note.* Student must not include eigenvectors corresponding
            to other eigenvalues. SLEPc returns these after last targeted one. For
            this case the dimension of :math:`E_{\omega^2}` is 2. Let\'s denote
            this bunch of vectors by ``E``.

            GS orthogonalization is called to tuple ``E+[f]``. This first
            orthogonalizes eigenvectors themself (for sure -- SLEPc doc is not
            conclusive about this) and then orthogonalizes ``f`` to
            :math:`E_{\omega^2}`.


.. only:: priv

    Reference solution
    ------------------

    .. toggle-header::
        :header: **Show/Hide Code**

        :download:`Download Code <impl.py>`

        .. literalinclude:: impl.py


.. [Evans] Lawrence C. Evans. *Partial Differential Equations.* Second edition.
           1998, 2010 AMS, Rhode Island.
