Elmer/Ice Release Notes for version 8.5
=======================================

Previous release: **8.5**
Period covered: **18 Dec 2018 - 30 Aug 2020**
Number of commits: **~110** (excluding merges of other branches)

These release notes provide information on most essential changes in Elmer/Ice functionalioties. Starting from the uppermost directory of the source tree, you can inquire changes inside the elmerice-directory using
```bash
git log --since="2018-12-18" -- elmerice
```

New Solver/Userfunction Modules
--------------------------------
- `Calving3D_lset`: return calving as a level set function (work in progress).
- `CalvingRemeshMMG`: cut a calving event directly out of a 3D mesh without external gmsh or mesh extrusion. Initial work on allowing calving margins to migrate.

Enhanced Solver/Userfunction Modules
------------------------------------
- `GlaDSCoupledSolver`: modified to work on a secondary hydrology mesh (as opposed to the primary ice mesh) and to discriminate properly between fjord-connected ungrounded areas and isolated ungrounded patches inland. Also should work on the basal boundary of an internally extruded 3D mesh.
- `GlaDSchannelSolver`:  changes to achieve the same outcome as above.

ElmerSolver library functionality
---------------------------------
- Added Zoltan repartitioning capabilities to permit continuous runtime load balancing and to assist with calving remeshing.
- Added support for MMG3D remeshing/mesh adaptation.
