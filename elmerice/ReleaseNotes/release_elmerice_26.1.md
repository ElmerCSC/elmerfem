Elmer/Ice Release Notes for version 26.1
========================================

Previous release: **9.0**  
Period covered: **Nov 11, 2020 - Jan 19, 2026**  

Trying to select the for Elmer/Ice relevant commits is difficult. The list obtained with
```bash
git log --since="2020-11-11" --stat -- elmerice
```
revealed 371 commits since Nov. 11th, 2020. We try to give a summary of the most relevant new developments below.

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from Markus Mützel, Fabien G-C, Iain W.,  Rupert G., Julien B., Samuel C.,  Benjamin R., Mondher C., Olivier G., Joe T.,  Eef v.D, and Lucas B. to Elmer/Ice related code in this release. Many thanks to all contributors supporting the open source philosophy of Elmer/Ice!

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
### CalvingRemeshparMMG.F90
- Similar to CalvingRemeshMMG.F90 except the second remesh (improvement of element quality) occurs in parallel using ParMMG.
- First stage of remeshing (element alignment along zero contour) still occurs in serial
- Takes a level set which defines a calving event (or multiple calving events). Levelset is negative inside a calving event, and positive in the remaining domain.
- ParMMG is currently no longer being developed  so CalvingRemeshMMG is more stable and actively maintained.

### CovarianceUtils.F90:
   - The CovarianceUtils module contains routines and functions for the generic operations involving a covariance matrix $C$ or its inverse
   - In particular an approximation of $B$ for Matérn functions can be obtain in an unstructured mesh by successive applications of a diffusion operator (Guillet et al.,2019)
   - This is used in:
     1. **BackgroundErrorCostSolver.F90**:
        - This solver is mainly intended to be used with the adjoint inverse methods implemented in Elmer/Ice and is an alternative to the regu
larisation solver
        - It allows to parametrise background errors with the Matérn correlation functions     
     2. **GaussianSimulationSolver.F90**:
         - Generate non-conditionnal realizations of a random variable from a gaussian distribution defined by the covariance matrix $C$
         - This can be use for, e.g., ensemble uncertainty quantification
  
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
    

  
New Elmer/Ice user-functions since last release are:

### USF_proj.F90
- generic user functions to compute longitude and latitude from projected x,y coordinates and conversely
- relies on generic utilities in the module file projUtils

### USF_GlacierMeshMetric.F90
- Computes the anisotropic target element size based on distance from calving front


III. Enhancement of existing solvers
------------------------------------
This is a list of selected improvements to existing solvers/methods

### Calving (CalvingRemeshMMG, Calving3DLset and CalvingGlacierAdvance3D)
- Update and production of CalvingRemeshMMG and Calving3DLset as usable solvers. New 3D calving model stable and works for any setup. See Wheel et al., 2024, GMD (https://doi.org/10.5194/gmd-17-5759-2024) for details.
- CalvingGlacierAdvance3D also now a stable solver.
- Stochastic calving CDL  (Crevasse Depth Law) added.
- Add option to remesh full terminus or just area that is calving
- Squashed commit of developments related to 3D calving of ice in St.Andrews between in years 2019-2023. Most of these updates covered in Wheel et al., (2024).
see https://github.com/ElmerCSC/elmerfem/commit/a9564a49d0a0e9d478d894886f5390833f28bd48
  - There are >300 commits over a broad spectrum of features needed in modeling of calving.
  - Most of the commits related to Iain Wheel's iw43@st-andrews.ac.uk work.
  - Also numerous commits by Joe Todd and Eef van Dongen
- This is a huge effort! Thank you!
- A squash was done in order not to mess up with the history of devel branch of Elmer over past four years.

### Inverse Methods
- optimisation of the (non-linear) Weertman friction coefficient
- add Material for the Mass conservation method
- "Regularisation" using a background error (see new **BackgroundErrorCostSolver.F90**)

### Permafrost model
-  dissable mass-matrix computation in steady state simulation of PermaFrost model
- added Interfrost TH1 special case
- introduced possibility to read solvent database into permafrost model (earlier default to pure water)
- Add flux relaxation to Darcy solver
- Added possibility of constant (exported) Temperature in PermafrostPorosityEvolution

### GlaDS
- allowing GlaDS to prescribe Moulins using a mask variable
- Changing default water density name in GlaDS from "water density" to "fresh water density" due to naming conflict with floatation code.
- Added possibility to limit the effective pressure in GlaDS channel evolution to positive values
- GlaDS volume source is now treated as handle (can alos be elment or IP value)
- Enable consistent norm for GladsCoupled. This should be consistent for parallel computation and for cases where there is different number of passive nodes and edges in the problem. Also utilize 'Mesh => Solver % Mesh' as it is used tens of times.

### ComputeDevStressNS.F90
ComputeDevStress is now compatible with both Stokes and Porous. Operations required to compute effective viscosity and compressibility parameter are performed in separated modules  through call of functions EffectiveViscosity PorousEffectiveViscosity

### PorousSolve.F90
Operations required to compute effective viscosity and compressibility parameter are now performed in a separated module PorousMaterialModels.F90 through call of function PorousEffectiveViscosity (replacing  fAandfB_in.F90)

### SSASolver    
- enable CutFEM to be used in assembly
- added linesearch for non-linear iteration into SSABasalSolver
- add support to compute grounding line and calving front fluxes
- add post-processing options:
  1. compute element-average basal stress
  2. compute nodal effective friction coefficient
- Move the friction law in a separate module (SSAMaterialModels) to ease use in other pieces of code
- add "regularised coulomb" friction law
- read sealevel at each visit so that it can change with time
    
