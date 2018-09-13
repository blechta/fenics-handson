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

* Source code: <https://github.com/blechta/fenics-handson>
* Official build: <https://fenics-handson.readthedocs.io/>


Prerequisities
--------------

Linux shell
^^^^^^^^^^^

This tutorial focuses solely on Python (version 3) user interface
of FEniCS. It is intended for complete newcomers to FEniCS.
Nevertheless it is assumed that participants (readers) have certain
experience in programming, for example in are of numerical algorithms.
Taking this assumption as granted it will not be presented how to use
Bourne Again Shell (Bash). We strongly suggest participants who have
no experience with Linux terminal to approach the closest Linux
machine at their workplace, make a coffee, find some tutorial
lectures (e.g. https://linuxjourney.com/lesson/the-shell) and
spend one or two hours on that.

FEniCS installation
^^^^^^^^^^^^^^^^^^^

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


Version
-------

.. todo:: Clarify FEniCS version used
