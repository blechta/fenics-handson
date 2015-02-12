Apendix
=======

Mesh import
-----------

GMSH + dolfin-convert [http://geuz.org/gmsh]

1.  use GMSH to create your domain, mark boundaries either by gmsh gui or by preparing gmsh yourmesh.geo script
2.  generate the mesh by `gmsh -2 yourmesh.geo` 
3.  run: dolfin-convert yourmesh.msh yourmesh.xml
4.  then in the FEniCS script you can read the xml mesh:

.. code-block:: python

   mesh = Mesh("yourmesh.xml")
   domain= MeshFunction("size_t", mesh, "yourmesh_physical_region.xml")
   bndry= MeshFunction("size_t", mesh, "yourmesh_facet_region.xml")


Local mesh refinement
---------------------

.. code-block:: python

   markers = CellFunction("bool", mesh)
   markers.set_all(False)
   for c in cells(mesh):
       for f in facets(c):
           if near(f.midpoint()[1],1.0) : markers[c] = True
   mesh = refine(mesh, markers, redistribute=False)

More details on ...
