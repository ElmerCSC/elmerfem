Elmer/Ice Release Notes for version 8.5
=======================================

Previous release: **8.5**

Period covered: **18 Dec 2018 - 30 Aug 2020**

Number of commits: **~110** (excluding merges of other branches)

These release notes provide information on most essential changes in Elmer/Ice functionalities. Starting from the uppermost directory of the source tree, you can inquire changes inside the elmerice-directory using
```bash
git log --since="2018-12-18" -- elmerice
```

New Solver/Userfunction Modules
--------------------------------
- `Calving3D_lset.F90`: Return calving as a level set function (work in progress).
- `CalvingRemeshMMG.F90`: Cut a calving event directly out of a 3D mesh without external gmsh or mesh extrusion. Initial work on allowing calving margins to migrate.
- `PlumeSolver.F90` (associated ODEPack library files: `opkda1.F`, `opkda2.F`, `opkdmain.F`): Provides plume melt rates across the calving front of a glacier. Fed by output from GlaDS solvers. Simulates a continuous sheet-style plume across entire front, split up into segments defined by frontal nodes and mesh resolution.

- `CalvingHydroInterp.F90`: Interpolates required variables between 3D ice mesh and 2D hydrology mesh, if using a multi-mesh approach. This is more complicated than it sounds.

- `HydroRestart.F90`: Allows separate 2D hydrology mesh to be restarted in a multi-mesh simulation.
- `USF_SourceCalcCalving.F90`: Uuser function that calculates the source term for GlaDS as a combination of surface melt (provided in some user-specified variable or input file) and basal melt (worked out automatically from the residual of the TemperateIce solver)
- `BasalMelt3D.F90`: Solver that works out basal melt on ungrounded portions of a glacier.
- `GMValid.F90`: Solver that discriminates between ungrounded areas that are connected to the fjord and isolated ungrounded patches inland.

Enhanced Solver/Userfunction Modules
------------------------------------
- `GlaDSCoupledSolver.F90`: modified to work on a secondary hydrology mesh (as opposed to the primary ice mesh) and to discriminate properly between fjord-connected ungrounded areas and isolated ungrounded patches inland. Also should work on the basal boundary of an internally extruded 3D mesh.
- `GlaDSchannelSolver.F90`:  changes to achieve the same outcome as above.
- `CalvingRemesh.F90` and `Calving3D.F90`: changed to avoid interpolating hydrology-specific solvers to the ice mesh after calving. Also changed to allow ice solvers and calving to run at different timestep to hydrology.
- `GroundedSolver.F90: Minor tweak to allow frontal grounded basal nodes to be listed as grounding-line nodes, so that the plume solver knows where to stick plumes


ElmerSolver library functionality
---------------------------------
- Added Zoltan repartitioning capabilities to permit continuous runtime load balancing and to assist with calving remeshing.
- Added support for MMG3D remeshing/mesh adaptation.
