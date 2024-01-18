# Calving Remesh MMG
## General Information
- **Solver Fortran File:** CalvingRemeshMMG.F90
- **Solver/User Function Name:** CalvingRemeshMMG
- **Required Output Variable(s):** {None}
- **Required Input Variable(s):** Calving Lset - currently hardcoded but should update so the user can specify
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None
- **Solver Keywords:**
  Keywords listed in solver:
  - Mesh Hausd(n) = Real 50.0 40.0 ... specify the Hausdorff distance for the first step of remeshing. See algorithm document for info on multiple values.
  - Mesh HMin(n) = Real 50.0 40.0 ... specify the the minimum element edge size for the first step of remeshing.
  - Mesh Hmax = Real 500.0 specify the maximum element edge size for the first step of remeshing.
  - Mesh Hgrad = Real 2.5 specify the gradation value from the first step of remeshing.
  - Mesh Min Quality = Real 0.00001 specify the minimum acceptable element quality for the first step of remeshing. See algorithm document for more details (Acceptable values 0-1).
  - Remeshing Distance = Real 1500.0. The distance from the terminus that will be remeshed.
  - Mesh Rmc Value = Real 1e-15 specify the volume threshold to remove negative parasitic elements. Value is a proportion of the total remeshing volume.
  - Front Advance Solver = String "Front Advance" - give the location of the Front Advance Solver so the fjord rails can be acccounted for when remeshing.
  - Fix Nodes on Rails = Logical True. Remeshing can cause the glacier to 'cut corners' of the fjord geometry. This option enforces the fjord boundaries. See algorithm document for more details.
  - Save MMGLS Meshes = Logical False. Save the Mmg .mesh mesh file prior and after the first step of remeshing. Useful for debugging.
  - Save MMGLS Sol = Logical False. Save the Mmg metric file.
  - If either the Mmg mesh or metric is being saved then the following is required.
  -- Pre MMGLS Mesh Name = File "outputmesh". Output mesh name prior to remeshing. Note no extension as the extension plus timestep will be automatically added.
  -- MMGLS Out Mesh Name = File "postremesh". Output mesh name post remeshing.
  - Mesh Update Variable n = String "variablename". Give the mesh update variables if present so they can be turned on/off by the calving solver.
  - Free Surface Variable n = String "variablename". Give the free surface variables if present so they can be turned on/off by the calving solver.
  - Switch Off Equation n = String "equationname". Give the name of any solver which will alter the mesh e.g., Front advance solver.
  - Pause Solvers Minimum Iceberg Volume = Real 1e7. This is the threshold of calved ice at which the adaptive timestepping will kickin. This is the maximum volume of one iceberg not the cumulative iceberg discharge.
  - Calving Pause Max Steps = Integer 3. Prescribe the maximum number of timesteps calving can pause the other solvers.
  - Calving Stats File Name = File "calvingstats.txt". Name and location of calving stats output file.
  - Save Terminus = Logical True. Option to turn on/off saving the terminus as a txt file.
  - Output Terminus File Name = File "terminusposition.txt". Name and location of terminus output file.
  - Suppress Calving = Logical False. Option to prevent calving. USeful for debugging or if just remeshing with advance is required.
  - Repartition Method = String "Parmetis". Other option "Zoltan". See Zoltan users guide for details on the difference in algorithms and options (./contrib/Zoltan/doc/Zoltan_pdf/ug.pdf in the Elmer repository).
  - Repartition Approach = String "Repartition". Determine time spent on remeshing. Options: 'Partition' - from scratch, 'Repartition' - repartition taking into account current distribution - recommended for dynamic load balancing, and 'Refine' - quickly improve current distribution.
  - Repartition Out Level - verbosity.
  - If ParMetis selected then the following options are required:
  -- Repartition ParMetis Library = String "Kway". Select repartitioning algorithm. See Zoltan users docs for options.
  - If Zoltan selected then the following options are required:
  -- Repartition Zoltan Library = String "graph". Select repartitioning algorithm. See Zoltan users docs for options.
  -- Repartition Zoltan Graph Package = String "phg". Select the graph or hypergraph package. See Zoltan users docs for options.

  Further to the options above that are defined in the solver some options are defined in the Materials. They are related to the second step of remeshing and are listed below.

  - MMG Hausd(n) = Real 50.0 40.0 ... specify the Hausdorff distance for the second step of remeshing. See algorithm document for info on multiple values.
  - MMG HMin(n) = Real 50.0 40.0 ... specify the the minimum element edge size for the second step of remeshing.
  - MMG Multiple Inputs = Logical True allows the mutliple values of Hmin and Hausd to be specified as above. Required for calving.
  - MMG Max Remesh Iterations = Integer 5 specify the remeshing iterations e.g. the size of Hmin and Hausd array.
  - MMG Hmax = Real 500.0 specify the maximum element edge size for the second step of remeshing.
  - MMG Hgrad = Real 2.5 specify the gradation value from the second step of remeshing.
  - MMG Min Quality = Real 0.00001 specify the minimum acceptable element quality for the second step of remeshing. See algorithm document for more details (Acceptable values 0-1).
  - MMG Anisotropic = Logical True. Select whether remeshing is anisotropic.
  - MMG Target Length(3) = ... Define the mesh metric. Three input values if anisotropic, one if isotropic. The user function GlacierMeshMetricAniso can be used for this and shown in the example below.
  - Save MMG Meshes = Logical False. Save the Mmg .mesh mesh file prior and after the first step of remeshing. Useful for debugging.
  - Save MMG Sol = Logical False. Save the Mmg metric file.
  - If either the Mmg mesh or metric is being saved then the following is required.
  -- Pre MMG Mesh Name = File "outputmesh". Output mesh name prior to remeshing. Note no extension as the extension plus timestep will be automatically added.
  -- MMG Out Mesh Name = File "postremesh". Output mesh name post remeshing.
  
## General Description
This solver takes a level set or signed distance variable where the zero contour indicates the new terminus and physically implements it by producing a new mesh. An overview of how the calving 3D solvers fit together is described in Calving.md. Further details of the new calving algorithm can be found [here](https://zenodo.org/records/10182710).

### Requirements {}
Elmer needs to be compiled with Mmg. This can be done by downloading and compiling Mmg (version 5.5.4 or greater) then adding -DMMG_INCLUDE_DIR="/include/file/locations/" and -DMMG_LIBRARY="/library/location/libmmg.so" to your Elmer installation script.

## Known Bugs and Limitations
Sometimes remeshing fails due to complex topography. Regular checking pointing will allow the model restart prior to this error occuring. See further documentation for details.

## SIF Contents
The required keywords in the SIF file for this solver/USF are:

```
Solver n
  Equation = "Remesh"
  Procedure = "ElmerIceSolvers" "CalvingRemeshMMG"
  Exec Solver = "After timestep"
  Solver Timing = Logical True

  Mesh Hausd (10) = Real 15.0 15.0 15.0 15.0 10.0 5.0 4.0 3.0 2.0 1.0
  Mesh Hmin (10) = Real 50.0 40.0 30.0 20.0 10.0 5.0 4.0 3.0 2.0 1.0
  Mesh Hmax = Real 500.0
  Mesh Hgrad = Real 2.5 ! needs to be very low for discretize an isovalue
  ! elem quality is out of 1
  Mesh Min Quality = Real 0.000001
  Remeshing Distance = Real 1500.0
  Mesh Rmc Value = Real 1e-15

  Front Advance Solver = String "Front Advance"
  !enforces lateral margins default = true
  Fix Nodes On Rails = Logical True

  Save MMGLS Meshes = Logical False
  Save MMGLS Sols = Logical False
  ! don't add extensions
  Pre MMGLS Mesh Name = File "mmgadvancefix/prels_"
  MMGLS Output Mesh Name = File "mmgadvancefix/ls_"

  Mesh Update Variable 1 = String "Vertical Mesh Update"
  Mesh Update Variable 2 = String "Longitudinal Mesh Update"
  FreeSurface Variable 1 = String "Zs Top"
  FreeSurface Variable 2 = String "Zs Bottom"

  Switch Off Equation 1 = String "Front Advance"

  !now total calve volume
  Pause Solvers Minimum Iceberg Volume = Real 1.0E7
  Calving Pause Max Steps = Integer 05
  !this is calving pause timestep
  Pseudo SS dt = Real 1e-10

  Calving Stats File Name = File "calvingstats.txt"
  Save Terminus = Logical True
  Output Terminus File Name = File "terminusposition.txt"
  !Suppress Calving = Logical True
   
  Repartition Method = String "ParMetis"
  Repartition Approach = String "Partition"
  Repartition ParMetis Library = String "PartKway"
  Repartition Output Level = Integer 0
  Repartition Zoltan Library = String "graph"
  Repartition Zoltan Graph Package = String "phg"
End

Material 1
  MMG Hmin (5) = Real 50.0 40.0 30.0 25.0 20.0
  MMG Hmax = Real 500.0 !500.0
  MMG Hgrad = Real 2.5
  MMG Hausd (5) = Real 50.0 40.0 30.0 25.0 20.0
  ! elem quality out of 1
  MMG Min Quality = Real 0.0001
  MMG Anisotropic = Logical True
  MMG Target Length(3) = Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "GlacierMeshMetricAniso"
  Save MMG Meshes = Logical False
  Save MMG Sols = Logical False
  ! don't add extensions
  Pre MMG Mesh Name = File "mmgadvancefix/preremesh_"
  MMG Output Mesh Name = File "mmgadvancefix/remesh_"
  MMG Multiple Inputs = Logical True
  MMG Max Remesh Iterations = Integer 5

  GlacierMeshMetric Max Distance = Real 2000.0
  GlacierMeshMetric Min Distance = Real 200.0
  GlacierMeshMetric Max LC = Real 500.0
  GlacierMeshMetric Min LC = Real 50.0
  GlacierMeshMetric Vertical LC = Real 50.0
End


```
## Examples
An example using the full calving algorithm can be found here [ELMER_TRUNK]/elmerice/Tests/Calving3D_lset

## References
Iain Wheel PhD thesis - DOI: https://doi.org/10.17630/sta/611.
Full calving algorithm detailed [here](https://zenodo.org/records/10182710).

