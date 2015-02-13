Appendix
========

Mesh import
-----------

GMSH + dolfin-convert [http://geuz.org/gmsh]

1.  use GMSH to create your domain, mark boundaries either by gmsh gui or by preparing gmsh yourdomain.geo script
2.  generate the mesh by `gmsh -2 yourdomain.geo` or `gmsh -3 yourdomain.geo` for 3D mesh 
3.  run: dolfin-convert yourdomain.msh yourdomain.xml
4.  then in the FEniCS script you can read the xml mesh:

.. code-block:: python

   mesh = Mesh("yourdomain.xml")
   domain= MeshFunction("size_t", mesh, "yourdomain_physical_region.xml")
   bndry= MeshFunction("size_t", mesh, "yourdomain_facet_region.xml")


Local mesh refinement
---------------------

.. code-block:: python

   markers = CellFunction("bool", mesh)
   markers.set_all(False)
   for c in cells(mesh):
       for f in facets(c):
           if near(f.midpoint()[1],1.0) : markers[c] = True
   mesh = refine(mesh, markers, redistribute=False)

This will create new mesh object.
For adaptivity see `adaptive Poisson demo <http://fenicsproject.org/documentation/dolfin/1.5.0/python/demo/documented/auto-adaptive-poisson/python/documentation.html>`_ and `adapt method reference <http://fenicsproject.org/documentation/dolfin/1.5.0/python/programmers-reference/cpp/fem/adapt.html#dolfin.cpp.fem.adapt>`_.

Application build on top of FEniCS
----------------------------------

See `Applications page of FEniCS project <http://fenicsproject.org/applications/>`_.
