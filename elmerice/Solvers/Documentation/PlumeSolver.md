# Solver PlumeSolver
## General Information
- **Solver/User Function Fortran File:** PlumeSolver.F90
- **Solver/User Function Name:** Plume
- **Required Output Variable(s):** MeltRate
- **Required Input Variable(s):**
  - Elevation (assumed to be called 'elevation')
  - Zb (assumed to be called Zb - defined on the hydrology mesh)
  - GroundedMask (assumed to be called this - defined on the hydrology mesh)
  - Channel Flux (assumed to be called this)
  - Sheet Discharge (assumed to be called this)
  - Sheet Thickness (assumed to be called this)
- **Optional Output Variable(s):**
  - Toe Calving (Exported Variable 1 = String 'Toe Calving')
- **Optional Input Variable(s):** 
  - GMCheck (assumed to be called this - the output of [GMValid.F90](./GMValid.md))
- **Solver Keywords:**
  - Plume Melt Mode = String... (can be 'constant', 'seasonal' or 'off')
  - (Summer/Winter) Salinity Temp Depth Input File = File... (if using constant plume forcing, just specify the one file; if using seasonal, specify a summer and winter file)
  - Force Toe Calving = Logical... (activate or suppress the toe calving routine (see below))
  - Mesh Resolution = Real... (the nominal mesh resolution at the calving front; see below)
  
## General Description
This is a new solver that allows meltwater plumes and their resulting melt rate to be modelled at the calving front. The actual plume model itself is a 1D ODE model, which is an adaptation of Tom Cowton’s MITgcm plume model, which is itself an adaptation of Donald Slater’s MATLAB plume model. It uses the ODEPack library (consisting, here, of odpka1.F, odpka2.F and odpkmain.F), which should be compiled alongside it. ODEPack is written in FORTRAN 77 – I made the minimum changes necessary to get it to compile with the elmerf90 compiler, but if the compiler gets updated at some point, it may find more things inside ODEPack that it doesn’t like. Also note that ODEPack is very particular about the format of its inputs and code that calls it (this is all detailed in odpkmain.F), hence why you’ll notice that outdated things like implicit variables and common blocks appear in the bowels of the PlumeSolver.F90 file. ODEPack is not currently included as part of the Elmer repository, and if you download your own version, it won’t compile. We’re working on setting it up to be included (or to use an alternative library), but, for now, you’re best contacting me (samuel.cook@univ-grenoble-alpes.fr) to get the library files.
PlumeSolver.F90 has three main subroutines:
* Plume: this is the wrapper routine that handles the interaction with the rest of Elmer, chiefly getting the necessary inputs from GlaDS and solver options, and then turning the outputs into a melt rate across the calving front
* PlumeSolver: this is the actual plume model that takes the input discharges and works out a resulting plume profile
* SheetPlume: this defines the system of equations actually solved by ODEPack
The Plume subroutine is a substantially modified and expanded version of a solver written by Joe Todd that takes a provided set of plume profiles and melt rates at fixed locations and applies them to the relevant bits of the calving front. All this functionality still exists and works, if anyone wants to use it (and some of it is hijacked by the new code anyway), but I’m not going to detail it here, because it’s not what I’ve focused on.
The overall strategy of the solver is to dynamically model a continuous line plume along the whole calving front. Each grounding-line node on the hydrology mesh is assigned to the nearest basal frontal node on the ice mesh, with the discharge of the plume at that point being the sum of the sheet and channel discharge across all of the hydrology nodes assigned to that ice node. The plumes are all then modelled to get a set of melt-rate profiles across the calving front, and the melt rate at each node on the calving front is then interpolated from this. This allows melt rates across the calving front to vary as the subglacial drainage evolves over time, and avoids the user having to specify the location or size of any of the plumes.
As things stand, the solver assumes the calving front is a flat, vertical plane. The plume model itself can handle non-vertical profiles, but I haven’t yet got round to working out how to extract this information from Elmer/how far the nature of the internally-extruded meshes in Elmer I’m using even allows non-vertical profiles to exist.

### Salinity Temp Depth Input File Formatting
This is the file that contains the data on the ambient water conditions. This should be a .csv or other text file with three columns – depth (ordered from 0 downwards, so a depth of 50 m should be expressed as -50 in the file), salinity (PSU) and temperature (°C)
Because of the strategy used by the solver to make plumes work properly in parallel, the sequence of depths in the ambient data file is that used to model the plume (it means the solver knows that the size of all the output profiles will be the same, which makes MPI much simpler). Therefore, the ambient data file **must** extend to the maximum depth of the calving front, so that any plume can be emplaced at the correct depth. This may mean you have to just add a few lines on the bottom of the file, copying your deepest data point downwards in increments of a few metres. Essentially, that’s what the toe calving routine (see below) does anyway, so it’s a reasonable assumption to make.

### Toe Calving
The upshot of this routine is to extend melt rates downwards, if your ambient data doesn’t extend as deep as the calving front. This stops unphysical toes forming at the calving front and messing up the mesh. It’s useful to export the toe calving variable anyway, even if you don’t want to activate the routine, as it (or a similar variable) is needed by GroundedSolver.F90 to help it define the grounding line

### Mesh Resolution
The solver needs this to work out the width of each plume segment and how far inland it should look for grounding-line nodes. This should work, even if you have a grounding line a long way inland, because the solver applies a couple of tests to work out if it should ignore a given grounding-line node on the hydrology mesh (e.g. it’s possible to have closed loops of ungrounded areas inland that don’t connect to the front and which the solver should ignore), but you may need to fiddle with this number slightly if you find it’s not doing what it should.

## Known Bugs and Limitations
None.

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "Plume"
  Procedure = File "ElmerIceSolvers" "Plume"
  Exec Solver = "After Timestep"

  Variable = "MeltRate"

  Force Toe Calving = Logical True
  Exported Variable 1 = "Toe Calving"
  
  Mesh Resolution = Real 100.0
   
  Plume Melt Mode = String "seasonal" !seasonal, constant, off
  ! Uncomment the first line if using a constant plume, or the remaining four if using a 
  ! seasonal one
  !! Salinity Temp Depth Input File = File "./Data/AmbientDataSummer.txt"
  !! Summer Salinity Temp Depth Input File = File "./Data/AmbientDataSummer.txt"
  !! Winter Salinity Temp Depth Input File = File "./Data/AmbientDataWinter.txt"
  !! Plume Melt Summer Start = Real $MeltSummerStart
  !! Plume Melt Summer Stop = Real $MeltSummerStop
End
```
## Examples
TO DO
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
None
