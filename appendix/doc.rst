Appendix
========

Mesh import
-----------

.. todo::

    Update: ``pygmsh``, ``meshio``, native ``gmsh`` bindings, ...

GMSH + dolfin-convert [http://geuz.org/gmsh]

1.  Use GMSH to create your domain; mark boundaries either by gmsh gui
    or by preparing gmsh script (see :download:`bench.geo` example).
2.  Generate the mesh by ``gmsh -2 yourdomain.geo`` or
    ``gmsh -3 yourdomain.geo`` for 3D mesh
3.  Run: ``dolfin-convert yourdomain.msh yourdomain.xml``
4.  Then in the FEniCS script you can read the XML mesh:

::

   mesh = Mesh("yourdomain.xml")
   cell_function = MeshFunction("size_t", mesh, "yourdomain_physical_region.xml")
   facet_function = MeshFunction("size_t", mesh, "yourdomain_facet_region.xml")


Local mesh refinement
---------------------

::

   markers = MeshFunction("bool", mesh, mesh.topology().dim())
   markers.set_all(False)
   for c in cells(mesh):
       # Mark cells with facet midpoints near y == 1.0
       for f in facets(c):
           if near(f.midpoint()[1], 1.0):
               markers[c] = True
   mesh = refine(mesh, markers, redistribute=False)

This will create new mesh object.
For adaptivity see `adaptive Poisson demo
<https://fenicsproject.org/docs/dolfin/2018.1.0/python/demos/auto-adaptive-poisson/demo_auto-adaptive-poisson.py.html>`_.


What every computer scientist should know about floating-point arithmetic
-------------------------------------------------------------------------

(Goldberg, 1991) `pdf <http://www.karlin.mff.cuni.cz/~hron/NMMO403/
What_every_computer_scientist_should_know_about_floating-point_
arithmetic-Goldberg-1991.pdf>`_


Periodic table of finite elements
---------------------------------

`femtable.org <http://femtable.org/>`_
