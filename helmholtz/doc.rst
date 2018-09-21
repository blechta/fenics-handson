Eigenfunctions of Laplacian and Helmholtz equation
==================================================

Wave equation with time-harmonic forcing
----------------------------------------

Let's have wave equation with special right-hand side

.. math::
    :label: wavespec

    w_{tt} - \Delta w &= f\, e^{i\omega t}
        &&\quad\text{ in }\Omega\times(0,T),

    w &= 0
        &&\quad\text{ on }\partial\Omega\times(0,T)

with :math:`f \in L^2(\Omega)`. Assuming ansatz

.. math::

    w(t, x) = u(x) e^{i\omega t}

we observe that :math:`u` has to fulfill

.. math::
    :label: helmholtz1

    -\Delta u - \omega^2 u &= f
        &&\quad\text{ in }\Omega,

    u &= 0
        &&\quad\text{ on }\partial\Omega.


.. admonition:: Task 1

    Try solving :eq:`helmholtz1` in FEniCS with data

    .. math::
        :label: helmholtzdata

        \Omega &= (0,1)\times(0,1),

        \omega &= \sqrt{5}\pi,

        f &= x + y

    on series of refined meshes. Observe behavior
    of solution energy :math:`\|\nabla u\|_2` with refinement.
    Is there a convergence or not?


Define eigenspace of Laplacian (with zero BC) corresponding
to :math:`\omega^2` as

.. math::

    E_{\omega^2} := \biggl\{ u\in H_0^1(\Omega): -\Delta u = \omega^2 u \biggr\}.

:math:`E_{\omega^2}\neq\{0\}` if and only if
:math:`\omega^2` is an eigenvalue. Note that
:math:`E_{\omega^2}` is finite-dimensional. Now define
:math:`P_{\omega^2}` as :math:`L^2`-orthogonal projection
onto :math:`E_{\omega^2}`. It is not difficult
to check that the function

.. math::
    :label: wavespecsol

    w(t, x) = \frac{t e^{i\omega t}}{2i\omega} (P_{\omega^2} f)(x)
    + e^{i\omega t} u(x)

solves :eq:`wavespec` provided :math:`u` fulfills

.. math::
    :label: helmholtz2

    -\Delta u - \omega^2 u &= (1-P_{\omega^2}) f
        &&\quad\text{ in }\Omega,

    u &= 0
        &&\quad\text{ on }\partial\Omega.

Note that problem :eq:`helmholtz2` has a solution which is
uniquely determined up to arbitrary function from :math:`E_{\omega^2}`.


.. admonition:: Task 2

    Construct basis of :math:`E_{\omega^2}` by numerically solving
    the corresponding eigenproblem with data :eq:`helmholtzdata`.

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

    Implement projection :math:`P_{\omega^2}`. Use it to solve
    problem :eq:`helmholtz2` with data :eq:`helmholtzdata`.


.. admonition:: Task 5

    Construct the solution :math:`w(t, x)` of the wave
    equations :eq:`wavespec` using formula :eq:`wavespecsol`.
    Plot temporal evolution of its real and imaginary
    part.


Mesh generation by Gmsh
-----------------------

.. admonition:: Task 6

    Modify a Gmsh demo to mesh a half ball

    .. math::

        \{(x,y,z), x^2 + y^2 + z^2 < 1,  y>0\}

    using the following code:

    .. code-block:: shell

        wget https://gitlab.onelab.info/gmsh/gmsh/blob/ad0ab3d5c310e7048ffa6e032ccd4e8f0108aa12/demos/api/boolean.py
        source /LOCAL/opt/gmsh-4.0.0/gmsh.conf
        python3 boolean.py
        meshio-convert -p -o xdmf-xml boolean.msh boolean.xdmf
        paraview boolean.xdmf &

        <edit> boolean.py

        python3 boolean.py
        meshio-convert -p -o xdmf-xml boolean.msh boolean.xdmf

    If in a need peek into
    ::

    >>> import gmsh
    >>> help(gmsh.model.occ.addSphere)


.. admonition:: Task 7

    Find :math:`E_{\omega^2}` with :math:`\omega^2 \approx 70`
    on the half ball. Plot the eigenfunctions in Paraview.

    .. hint::

        Use ``Glyph`` filter, ``Sphere`` glyph type, decrease
        the scale factor to ca. 0.025.

        Use ``Clip`` filter. Drag the clip surface by mouse,
        hit ``Alt+A`` to refresh.


.. only:: priv

    Reference solution
    ------------------

    .. toggle-header::
        :header: **Show/Hide Code**

        :download:`Download Code <impl.py>`

        .. literalinclude:: impl.py
