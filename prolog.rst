Prolog
======

.. sidebar:: Motto

    "My theory by A. Elk.  Brackets Miss, brackets.  This theory goes
    as follows and begins now.  All brontosauruses are thin at one
    end, much much thicker in the middle and then thin again at the
    far end.  That is my theory, it is mine, and belongs to me and I
    own it, and what it is too."

    -- Anne Elk (Miss)


This document served primarily as task sheets for
`FEniCS hands-on lectures held on Chemnitz University
of Technology in September 2018
<https://www.tu-chemnitz.de/mathematik/part_dgl/teaching/WS2018_FEniCS>`_.
Nevertheless it is not excluded that these sheets could not be
used separately or for any other occassion.


Target
------

This tutorial gives lectures on usage of FEniCS version 2017.2.0
through its Python 3 user interface. It is specifically intended
for newcomers to FEniCS and as such does not assume any knowledge in
Python programming. Rather than taking a Python tutorial first,
the intent is to learn-by-doing. As a consequence first steps
consist of modificating existing FEniCS demos while gradually
taking bigger and bigger tasks in writing original code.


Resources
---------

.. todo::

    Add FEniCS and Python resources


FEniCS installation
-------------------

Obviously we will need a working installation of FEniCS.
`FEniCS can be installed in different ways
<https://fenicsproject.org/download/>`_ which all of them
have some pros and cons. On TU Chemnitz this taken care
of by organizers of the hands-on and participants do not
have to worry about this.

Nevertheless participants might want to install FEnicS
to their laptops, workstation, home computers to practice
or use FEniCS outside of the tutorial classes. The easiest
option for new FEniCS users on Ubuntu is to install using
APT from FEniCS PPA:

.. code-block:: bash

    sudo apt-get install software-properties-common
    sudo add-apt-repository ppa:fenics-packages/fenics
    sudo apt-get update
    sudo apt-get install --no-install-recommends fenics

On the other hand FEniCS images for Docker provide the most portable
solution, with arbitrary FEniCS version choice, for systems where
`Docker CE <https://www.docker.com/community-edition>`_ can be installed
and run; see https://fenicsproject.org/download/.

Resources
---------

* `DOLFIN docs <https://fenics.readthedocs.io/projects/dolfin/en/2017.2.0>`_
* `DOLFIN Python API docs <https://fenicsproject.org/docs/dolfin/2017.2.0/python/index.html>`_
* `Python docs <https://docs.python.org/3>`_
* `FEniCS AllAnswered <https://www.allanswered.com/community/s/fenics-project/>`_
