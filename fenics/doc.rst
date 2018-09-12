Very short introduction to FEniCS
=================================

First touch
-----------

Login by SSH to ``tyche`` and type:

.. code-block:: bash

    source /LOCAL/Software/FEniCS-2018.1/setup_env

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


Run and modify demos
--------------------

.. admonition:: Task 1

    Get the Poisson demo from FEniCS install dir and run it:

    .. code-block:: bash

        mkdir -p work/fenics/poisson
        cd work/fenics/poisson
        cp /LOCAL/Software/FEniCS-2018.1/share/dolfin/demo/python/documented/poisson/demo_poisson.py .
        python3 demo_poisson.py

    You should see some console output and a plot of the solution.


Now login to ``tyche`` from another terminal window and open
the demo file using your favourite editor (if you don't have any
you can use ``gedit``, ``nano``, ...):

.. code-block:: bash

    cd work/fenics/poisson
    <editor> demo_poisson.py


.. admonition:: Task 2

    Now add *keyword argument* ``warp='mode'`` to the ``plot`` function
    call by applying the following diff:

    .. code-block:: diff

         # Plot solution
         import matplotlib.pyplot as plt
        -plot(u)
        +plot(u, mode='warp')
         plt.show()

    and run the demo again by ``python3 demo_poisson.py``.


Open `Poisson demo documentation
<https://fenicsproject.org/docs/dolfin/2018.1.0/python/demos/poisson/demo_poisson.py.html>`_
on the FEniCS website. Notice that the doc page is generated from
the demo file. Go quickly through the docpage while paying attention
to

* definition of weak formulation through forms ``a`` and ``L``,
* usage of ``Constant`` and ``Expression`` classes.


.. admonition:: Task 3

    Now modify the problem to use the following data instead:

        .. math::

            a(u, v) &= \int_\Omega \nabla u\cdot\nabla v\,\mathrm{d}x
                     +         \int_\Omega c\,u v\,\mathrm{d}x

            c       &= 100

            f       &= x

            g       &= \sin(5x) \exp(y)

            u       &= y \qquad \text{on } \Gamma_\mathrm{D}
