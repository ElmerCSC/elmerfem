Elmer/Ice Release Notes for version 26.1
========================================

Previous release: **9.0**  
Period covered: **Nov 11, 2020 - Jan 19, 2025**  

Trying to exclude the for Elmer/Ice relevant commits is difficult. The list obtained with
```bash
git log --since="2020-11-11" --stat -- elmerice
```
revealed 371 commits since Nov. 11th, 2020.

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from Markus Mützel, Fabien G-C, Iain W.,  Rupert G., Julien B., Samuel C.,  Benjamin R., Mondher C., Olivier G., Joe T.,  Lucas B. to Elmer/Ice related contributions in this release.

### New Versioning Scheme

From this version onward we migrate for Elmer (and hence also Elmer/Ice therein) to Calendar Versioning such that
- First number (major) is the year of the 21st century, e.g. 26
- Second number (minor) is an ordinal number of releases in that year
- Third number (micro) is a growing number which for releases is always 0 and may be omitted.  


I. Overview of changes
---------------------
Just in time to the previous release the Elmer/Ice community introduced the concept of documentation inside this repository, namely,
- `elmerice/Solvers/Documentation` for Solvers
- `elmerice/UserFunctions/Documentation`for user-functions
Since then, these directories have been filled with the documents.



II. New Solver Modules
----------------------
The list obtained with
```bash
for i in $(ls *)
do echo $(git log --diff-filter=A --follow --format=%as -- "$i"), $i
done |sort
```
revealed for new Solver modules:

### IcyMaskSolver.F90

- Solver for creating a mask on whether there is ice-thicknes above a given threshold (H> Hmin) or not (H<Hmin).
- Similar to the 'groundedmask', we have  +1= glaciated (H>Hmin), -1= Ice Free (H<Hmin), 0=contour of the glacier, plus values <-1 for isolated nodes

### Scalar_OUTPUT_Glacier.F90

- This solver is to be used to output some scalar quantities for a glacier configuration (domain without ice characterised by an IcyMask < 0).
- The quantities are:

      - glacier volume
      - glacier area
      - ablation area
      - accumulation area
      - SMB total
      - SMB Ablation
      - SMB Accumulation
      - Front elevation

### UGridDataReader.F90

- This solver reads variables in an unstructured netcdf file (e.g. following the UGRID format) at node and element locations.
- It can be used to e.g. read variables that have been produced with the XIOSOutPutSolver or that have been conservatively interpolated on the mesh using e.g. cdo.
- The input file structure should correspond to the current serial mesh, and variables should be arranged using the node and element ordering.

### Weertman2Coulomb.F90
- Converts linear Weertman coefficient (e.g. from inversion) to "Coulomb"
sliding law parameters
- The following "Conversion mode" options are available (see above document
for more info):
    - "Threshold": A threshold value of the Weertman sliding coefficient is given.  Either side of this threshold one or other of the "Coulomb" coefficients is held constant while the other is derived.
    -  "Smooth": A Weertman equation is used to calculate the As "Coulomb" coefficient.  This is then scaled toward zero for regions where the Coulomb limit is approached  (currently using a tanh function based on effective pressure).  The `C`-coefficient is then derived as a function of `As.
    
### CalvingRemeshparMMG.F90
- Remesh the calving model using MMG3D
- runs in parallel but remeshing is serial!
- Takes a level set which defines a calving event (or multiple calving events). Level
- set is negative inside a calving event, and positive in the remaining domain. This hasn't actually been implemented yet, we use a test function.

New Elmer/Ice user-functions since last release are

### USF_proj.F90
- generic user functions to compute longitude and latitude from projected x,y coordinates and conversely
- relies on generic utilities in the module file projUtils

### USF_GlacierMeshMetric.F90
- Computes the anisotropic target element size based on distance from calving front


III. Enhancement of existing solvers
------------------------------------

