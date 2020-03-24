# ElmerGrid mesh generation and manipulation

This directory includes files for ElmerGrid which is a mesh
generation utility bundled with Elmer project. 

ElmerGrid can create simple structured 2D meshes consisting
of quadrilaterals or triangles. Also limited 3D functionality is provided
in terms of extruded 2D meshes resulting to hexahedrons or prisms.

The program may also be used as a mesh import and export utility. It    
is able to read several different formats (Gmsh, UNIVERSAL file format,
COMSOL mphtxt files,..) and writes mainly Elmer input 
and output formats, and also VTK format.

The meshes may also be given some simple operations such as scaling,
rotation, cloning, extrusion ect. 

ElmerGrid includes also partitioning capabilities of Elmer meshes.
The partitioning may be done using internal geometric division or
Metis library. 

