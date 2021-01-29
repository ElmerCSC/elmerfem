# Hydro-Calving-Plumes Elmer/Ice Documentation
## Modified Files:
* GlaDSCoupledSolver.F90
* GlaDSchannelSolver.F90
* MeshUtils.F90
* CalvingRemesh.F90
* GroundedSolver.F90
* InterpVarToVar.F90
* ModelDescription.F90
* ElmerSolver.F90
## New Files:
* PlumeSolver.F90 (and associated ODEPack library files: opkda1.F, opkda2.F, opkdmain.F)
* CalvingHydroInterp.F90
* HydroRestart.F90
* USF_SourceCalcCalving.F90
* BasalMelt3D.F90
* GMValid.F90

## Background
This suite of code changes and additions to the Elmer/Ice framework are designed to allow the coupling and concurrent use of the 3D calving solvers, GlaDS hydrology module and a new plume model within a single Elmer/Ice simulation. They are aimed at situations where you may want to model a large, fast-flowing tidewater glacier and were developed using the test case of Store Glacier in Greenland. If your glacier is small, slow or land-terminating, or you’re modelling an entire ice sheet, you probably don’t need to worry about this, or at least not all of it! In intent, all these changes should be backwards-compatible and shouldn’t break anything that doesn’t use them, but may be using some of the components. If you do find that your old simulations don’t work anymore, it’s probably that I’ve forgotten to put an if statement somewhere to only use the new code in the right circumstances, in which case I apologise. Below I list the changes made to existing files, including new SIF keywords, and the purpose and mode-of-use of new ones, as well as some more general considerations for making everything work.

## THE MOST IMPORTANT CHANGE
All of this code will only function properly if you put `Calving = Logical True` in the Simulation section of the sif. Without that, none of it will activate properly (or indeed, at all); if you use it, it will all try to run in any of the modified or new files that are used. As things stand, this isn’t therefore terribly modular – you either turn it all on or none of it. If you just want calving, or just GlaDS, they work fine as standalones, but the plume solver requires the use of GlaDS to function (though it can work without the 3D calving solvers being used – this is explained more below).

## Timestepping
One thing you need to consider with this set-up is the timesteps different parts of the model need to run on. Obviously, liquid water moves and evolves much faster than solid water (i.e. ice), so, in practice, the hydrological solvers want to be running at a much smaller timestep than the ice solvers. However, you can’t run everything at the smallest timestep, unless you’ve got infinite computer time…. So, I found a timestep of 0.1 days for the hydrology worked well, with a timestep of 1 day for the ice (giving a runtime of about 30 hours for one year of simulation). This can be achieved by adding ‘Exec Intervals = 10’ to all the ice solvers (assuming your base timestep is 0.1 days, this means the ice solvers will only run every 10 timesteps, i.e. every day). If you do this, though, you **MUST** also add ‘Timestep Scale = Real 10.0’ to all those solvers, because Elmer is not intelligent enough to recognise, when it comes to time-dependent solvers, that they’re not running at the smallest timestep (essentially, Elmer just looks at the timestep size, which will be the smallest timestep, and assumes everything is running at that, rather than actually checking the time). Otherwise, your ice will evolve 10 times too slowly. Yes, I found this out the hard way.
You should also avoid, if possible, having more than two different timestep sizes for the sake of model stability. So, have your base timestep for the hydrology and some multiple of that for the ice, but try to not start having intermediate ones, because it’s more likely the model will fall over if you do that.

## Restarting runs
If you need to restart from fully coupled model runs, set up the restart machinery as normal, but then include the HydroRestart solver (see below) as a solver that is executed ‘Before Simulation’. This will allow all the ice and hydrology variables to be restarted properly. If calving has been going on, you’ll also need to make sure you use the right ice mesh – this will be the mesh found in the directory you’ve listed under ‘Remesh Move Mesh Dir’ in the Remesh solver. You’ll find one folder created in that directory for every time the model remeshed, so make sure you pick the right one (usually, this will be the one from the latest timestep, but if you’re restarting from a crashed run, you want to pick the one that lines up with the last result output timestep).

## Meshes
The main issue with coupling calving and hydrology within Elmer/Ice is that of meshing. The calving solvers rely on modifying the existing ice mesh, creating a new one, and interpolating all the variables between the two. With GlaDS, this doesn’t work, because the channel variables are defined on the edge elements of the mesh itself so i) modifying the mesh causes problems and ii) the channel variables are not obviously interpolatable. To get round this, this setup makes use of two meshes: the standard 3D, internally extruded ice mesh, on which all the usual ice and calving stuff happens, and a second 2D plane hydrology mesh that the GlaDS solvers work on.
As such, you need to create a secondary mesh that the model can use. Ideally, this should be the same footprint as your ice mesh, except you probably want it to extend a bit past the frontal boundary of your ice mesh. The calving model could see the glacier advance or retreat, so the hydrology mesh needs to have a footprint that covers the entire potential area that the calving front could reach. Resolution-wise, things tend to work best if the two meshes are as similar as possible – large differences in resolution will lead to lots of interpolation artefacts (discussed more fully under CalvingHydroInterp.F90 below) that will most likely mess up your simulation eventually. But, you do want the hydrology mesh to be at a finer resolution than the ice, so try a few things and see what works best.
Once you have got your 2D footprint hydrology mesh, you need to get it into Elmer format in the usual way (whatever that might be for you). When you do this, you’ll want to use the -bcoffset option to increment the boundary condition numbers on the new mesh so that they follow on from the boundaries on the ice mesh (so, if you have 4 boundaries on your ice mesh, plus a surface and basal boundary once it’s extruded, you’ll want to offset the hydrology mesh boundaries by 6 so that they start at 7). Once you’ve done that (and before you partition it), you need to edit the mesh.elements file and ensure the **second** column is exclusively populated with the number 2, rather than 1, which is what it’ll be by default. This tells Elmer that this is a different body to the main ice mesh. I find the awk command the easiest way to perform this replacement (awk ‘$2=”2”’). Once you’ve done that, you can partition the mesh as usual. You also need to copy the new mesh directory to create a renamed version that will be used by solvers that need the hydrology mesh, but that you don’t want to create results files from (this is something that may get fixed, but, as it stands, the results output solver will try to output vtu files from every individual solver mesh – if they’re all using the same mesh, all the vtu files will have the same name and will overwrite each other).
When you’ve done all this, there are a couple of things you need to do in the SIF:
* In the Simulation section of the SIF, put `Need Edges 2D = Logical True` and `Need Edges 3D = Logical False` – this will ensure edges are generated on all the 2D hydrology meshes, but not on the 3D ice mesh (edges are important for channel variables)
* In the solver section for GlaDSCoupledSolver, put `Mesh = “.” “<pathtoyourprimaryhydrologymesh>"` (assuming your mesh directory is some subdirectory of the working directory the SIF lives in)
* In the solver section for the other two GlaDS solvers (the thickness and channel output ones), and for any solvers where you’re reading in a variable that you want to be applied to the hydrology mesh (say, an expanded basal DEM covering the larger footprint area), put the same line, but replace the path with the equivalent for your renamed secondary hydrology mesh directory
* You need to use the Target Bodies feature to differentiate the ice and hydrology meshes. Assuming Body 1 is your main ice body, put the line `Target Bodies(1) = 1` in the Body 1 section of the SIF
* Then, you need to create a new Body 2 for the hydrology mesh (in theory, this can be any number – you could make it Body 7, if you already have Bodies 2-6, say, I think). This needs the line `Target Bodies(1) = 2` in it, though you could replace the ‘2’ with whatever number you substituted into the mesh.elements file earlier. Though I reckon Elmer would be angry if you, say, used Target Bodies(1) = 3 if there isn’t a Target Bodies(1) = 2 somewhere else. This body will have the same material, but a different equation, body force and initial condition to the ice, just as if you were defining a basal boundary body for the hydrology on a standard 3D mesh
* You can then list all the other bodies you’re going to need to define on the boundaries of the ice mesh in the usual way without needing to put in any more Target Bodies statements
* In the Boundary Condition section, you can list BCs for the hydrology mesh as usual, making sure you use the correct BC number (so, in the example above, you’d probably have BCs 7-10 defined on the hydrology mesh, and 1-6 on the ice mesh). But, make sure that all the hydrology-mesh-defined BCs include `Body ID = 2` (I’m not sure if, assuming you’ve numbered the BCs correctly or pointed them at the correct Target Boundary if you’re not setting them by numbering, you actually need this, but it doesn’t hurt)
* For the initial condition and the equation, just list them as usual – the equation will be any solvers you want to execute on the hydrology mesh, typically the three GlaDS solvers, any reader solvers that are reading variables onto the hydrology mesh, the hydrology weights calculator, and the hydrology restart solver (if relevant – all of these are properly explained below)
* In the relevant Body Force section, a few conditional boundary conditions need to be added (turns out Elmer is perfectly happy with boundary conditions in the body force section – they get applied to all nodes on the relevant body, which can be really useful). These work with some of the changes in GlaDSCoupledSolver.F90 (see below) to stop GlaDS trying to do hydrology on ungrounded parts of the mesh – strictly speaking, I should probably either put it all in the code, or put it all in as body force boundary conditions, but this is the currently working set-up:
  * `Hydraulic Potential = Real 0.0`
  * `Hydraulic Potential Condition = Variable GMCheck, NormalStress`
    * `Real MATC “if(tx(0)>0.0){1}else{(tx(1)*(-1))+1E-32}”`
  * The same form of boundary condition should be set for effective pressure, sheet thickness, sheet discharge and sheet storage, as well as the no channel BC (`= Logical True`)
  * Water pressure should instead be set to:
  * `Water Pressure = Variable Zb` (or whatever you’ve called your basal DEM variable)
    * `Real MATC “abs(tx*RhoWS*g)”`
    * Because water pressure just depends on depth if the ice is ungrounded
  * A simpler form of the first BC above is:
  * `Hydraulic Potential = Real 0.0`
  * `Hydraulic Potential Condition = Variable GroundedMask`
    * `Real MATC “(tx*-1)-0.5”`
  * This is fine, provided you don’t have any discrete ungrounded areas inland of the calving front and unconnected to it; if you do, you need the more complicated BC to stop these ungrounded patches turning into unphysical sinks of water

## Modified Files:
### GlaDSCoupledSolver.F90
Modifications here are about making sure everything points at Solver % Mesh, rather than Model % Mesh, and to do with correctly importing and saving the other hydrology variables defined on the other hydrology solvers. There are also several blocks that deal with determining whether to ignore elements based on their groundedness (if it’s ungrounded, there’s no hydrology, though the model will treat isolated ungrounded patches inland as grounded here, because otherwise they become unphysical sinks of water – this discrimination is done based on the GMCheck variable from GMValid.F90), and a few lines that stop channel area from blowing up to stupid proportions, which is a problem on fine-resolution meshes. These are activated by specifying `Max Channel Area = Real…` and `Max Sheet Thickness = Real…` in the solver section of the SIF, should you wish to use them.
In practical terms in the SIF, this doesn’t require anything else new, beyond the listing of the hydrology mesh in the Solver section. Otherwise, all that’s required is a line saying `Need Edges = True` in the solver section, which is picked up by MeshUtils.F90 (see below) and allows the channel variables to be added to the mesh from their (identical) solver mesh.

### GlaDSchannelSolver.F90
Similarly, changes here are to point everything at Solver % Mesh and towards the appropriate location of the channel variables on the primary hydrology mesh, rather than the secondary hydrology mesh this solver will initially be using. None of this requires anything new in the SIF.

### MeshUtils.F90
All the changes in here are pretty much to do with making sure internally-extruded meshes keep edges after extrusion, to make GlaDS work on them, and then making solver-specific meshes able to have edges without the primary solver variable being defined on them. This should all just work – in the current situation, the only line required is the Need Edges line in the GlaDSCoupledSolver section to make that solver meshes have edges, so it can have the channel variables added to it, even though hydraulic potential isn’t an edge variable. If you’re just using GlaDS on the base of a 3D internally extruded mesh, the relevant keyword is instead `Preserve Edges = Logical True` in the Simulation section of the sif.

### CalvingRemesh.F90
Very minor change to stop the solver trying to add the hydrology variables to the post-calving ice mesh. This is activated by a new solver option for the CalvingRemesh solver:
* `Solvers To Ignore(n) = Integer n1 n2 n3…`
  * n1, n2, n3 and so on should be the number of the solvers assigned to any of the hydrology meshes; n should be the total number of solvers concerned
The remesh solver needs to be run every timestep for complicated reasons, but will only actually (possibly) remesh if it’s a timestep the rest of the calving solvers have been run on. This is done by an addition to the file that compares the current time to the time a calving event occurred (which is recorded by an addition to Calving3D.F90). If the time is the same (i.e. the other calving solvers were also run this timestep), remeshing proceeds as normal; if the time is different (i.e. we’re on an intermediary hydrology timestep between ice timesteps), then the CalvingOccurs Boolean is set to false, and no remeshing occurs. The extra computational time required by running the solver in this way every timestep is negligible.

### GroundedSolver.F90
Now should list all grounded nodes on the frontal ice boundary as grounding-line nodes, which is needed for the PlumeSolver, provided a line in the solver section saying `Front Variable = String “<name>”` is added, specifying the name of a variable defined on the calving front (the toe calving variable is useful here – see below). Also added an option to force the entire domain to be grounded, which can be activated by putting a line saying `All Grounded = Logical True` in the solver section. This is useful because if you’re running a hydrology-only simulation where you want to use two meshes (say, for mesh resolution reasons, or because you want to restart a full calving-hydro-plume simulation from it), you’ll need to use the `Calving = True` switch, which will make GlaDSCoupledSolver look for the groundedness of elements to determine whether it needs to bother with them, so you’ll need to include the GroundedSolver in the sif, but you don’t actually want it to do anything beyond marking everything as grounded.

### InterpVarToVar.F90
A couple of modifications to allow it to deal properly with interpolating between 2D and 3D meshes, rather than just between 3D meshes, and for situations where mesh connectivity is not perfect (i.e. it now allocates perms that will actually work).

### ModelDescription.F90
One modification to stop it trying to re-allocate a variable in the case of multiple restarts happening in the same simulation (i.e. if you’ve got more than one mesh). Restarting in the usual manner will load the first mesh, which, unless you’ve set things up strangely, should be the ice mesh (which will also be the mesh defined as `Model % Mesh`). The hydrology mesh can be restarted using the HydroRestart solver described below.

### ElmerSolver.F90
A couple of minor modifications to stop the model trying to restart all the meshes. If you do nothing, the model will only restart the first mesh, which will usually be the ice mesh, if the `Calving=Logical True` switch is set in the sif. You can also specify in the Simulation section a new keyword, `Meshes To Restart(n) = n1, n2, n3…` if you want to modify this default behaviour.

## New Files:
### PlumeSolver.F90
This is the new solver I’ve written that allows meltwater plumes and their resulting melt rate to be modelled at the calving front. The actual plume model itself is a 1D ODE model, which is an adaptation of Tom Cowton’s MITgcm plume model, which is itself an adaptation of Donald Slater’s MATLAB plume model. It uses the ODEPack library (consisting, here, of odpka1.F, odpka2.F and odpkmain.F), which should be compiled alongside it. ODEPack is written in FORTRAN 77 – I made the minimum changes necessary to get it to compile with the elmerf90 compiler, but if the compiler gets updated at some point, it may find more things inside ODEPack that it doesn’t like. Also note that ODEPack is very particular about the format of its inputs and code that calls it (this is all detailed in odpkmain.F), hence why you’ll notice that outdated things like implicit variables and common blocks appear in the bowels of the PlumeSolver.F90 file. ODEPack is not currently included as part of the Elmer repository, and if you download your own version, it won’t compile. We’re working on setting it up to be included (or to use an alternative library), but, for now, you’re best contacting me (samuel.cook@univ-grenoble-alpes.fr) to get the library files.
PlumeSolver.F90 has three main subroutines:
* Plume: this is the wrapper routine that handles the interaction with the rest of Elmer, chiefly getting the necessary inputs from GlaDS and solver options, and then turning the outputs into a melt rate across the calving front
* PlumeSolver: this is the actual plume model that takes the input discharges and works out a resulting plume profile
* SheetPlume: this defines the system of equations actually solved by ODEPack
The Plume subroutine is a substantially modified and expanded version of a solver written by Joe Todd that takes a provided set of plume profiles and melt rates at fixed locations and applies them to the relevant bits of the calving front. All this functionality still exists and works, if anyone wants to use it (and some of it is hijacked by the new code anyway), but I’m not going to detail it here, because it’s not what I’ve focused on.
The overall strategy of the solver is to dynamically model a continuous line plume along the whole calving front. Each grounding-line node on the hydrology mesh is assigned to the nearest basal frontal node on the ice mesh, with the discharge of the plume at that point being the sum of the sheet and channel discharge across all of the hydrology nodes assigned to that ice node. The plumes are all then modelled to get a set of melt-rate profiles across the calving front, and the melt rate at each node on the calving front is then interpolated from this. This allows melt rates across the calving front to vary as the subglacial drainage evolves over time, and avoids the user having to specify the location or size of any of the plumes.
As things stand, the solver assumes the calving front is a flat, vertical plane. The plume model itself can handle non-vertical profiles, but I haven’t yet got round to working out how to extract this information from Elmer/how far the nature of the internally-extruded meshes in Elmer I’m using even allows non-vertical profiles to exist.
As regards solver options and inputs:
* `Plume Melt Mode = String “…”`
  * Options are `constant`, `seasonal` or `off`
  * If `constant`, you need to specify `Salinity Temp Depth Input File` (see below)
  * If `seasonal`, you need to specify `Summer Salinity Temp Depth Input File` and `Winter Salinity Temp Depth Input File` (see below)
  * Additionally, you need to specify `Plume Melt Summer Start = Real…` and `Plume Melt Summer Stop = Real…` in model time to say when your summer conditions start and stop
* `(Summer/Winter) Salinity Temp Depth Input File = File “Name”`
  * This is the file that contains the data on the ambient water conditions. This should be a .csv or other text file with three columns – depth (ordered from 0 downwards, so a depth of 50 m should be expressed as -50 in the file), salinity (PSU) and temperature (°C)
  * Because of the strategy used by the solver to make plumes work properly in parallel, the sequence of depths in the ambient data file is that used to model the plume (it means the solver knows that the size of all the output profiles will be the same, which makes MPI much simpler). Therefore, the ambient data file **must** extend to the maximum depth of the calving front, so that any plume can be emplaced at the correct depth. This may mean you have to just add a few lines on the bottom of the file, copying your deepest data point downwards in increments of a few metres. Essentially, that’s what the toe calving routine (see below) does anyway, so it’s a reasonable assumption to make
* `Force Toe Calving = Logical True/False`
  * This is something Joe implemented, the upshot of which is to extend melt rates downwards, if your ambient data doesn’t extend as deep as the calving front. This stops unphysical toes forming at the calving front and messing up the mesh
  * If using this, specify `Exported Variable 1 = “Name”` for the toe calving variable
  * It’s useful to export the toe calving variable anyway, even if you don’t want to activate the routine, as it (or a similar variable) is needed by GroundedSolver.F90 to help it define the grounding line
* `Mesh Resolution = Real n`
  * This is the nominal mesh resolution at the calving front, which the solver needs to work out the width of each plume segment and how far inland it should look for grounding-line nodes. This should work, even if you have a grounding line a long way inland, because the solver applies a couple of tests to work out if it should ignore a given grounding-line node on the hydrology mesh (e.g. it’s possible to have closed loops of ungrounded areas inland that don’t connect to the front and which the solver should ignore), but you may need to fiddle with this number slightly if you find it’s not doing what it should
Known issues:
* None yet

### CalvingHydroInterp.F90
This represents the other major block of new code written as part of this suite. It handles the interpolation of necessary variables between the ice and hydrology meshes, but also moves read-in variables (using GridDataReader or similar) from their solver-specific secondary hydrology meshes to the primary hydrology mesh associated with GlaDSCoupledSolver. It also corrects interpolation artefacts that will otherwise nix your simulation sooner or later, and ensures conservation of the temperature residual (one of the interpolated variables) to stop the glacier accidentally destroying or creating some energy….
The file contains two main subroutines, imaginatively titled “IceToHydroInterp” and “HydroToIceInterp”. Make sure you get them the right way round. IceToHydroInterp interpolates the ice normal stress, velocity, grounded mask and temperature residual over to the hydrology mesh and then spends a lot of time clearing up artefacts and conserving the temperature residual. HydroToIceInterp is much simpler, as the hydrology mesh is usually finer than the ice mesh, so the interpolation routine doesn’t create anywhere near as many artefacts in problematic locations. Therefore, it pretty much just interpolates the water pressure, effective pressure and sheet discharge onto the ice mesh.
There are also two small subroutines: “HydroWeightsSolver” and “IceWeightsSolver”. These calculate the boundary weights used in the main routines (the reason this happens in a separate solver is complicated – suffice to say it does exist). These need to be called as solvers before the relevant interpolation routines (IceToHydro or HydroToIce) are called, otherwise they’ll crash. IceWeightsSolver needs to run every time the ice mesh is updated (probably every n timesteps); HydroWeightsSolver every time the hydrology mesh is updated (probably never, so it can just run once at the start of the simulation).
Solver options and inputs (for IceToHydroInterp; others have nothing fancy):
* `Load Reader Variables = Logical True/False`
  * Set this option to true if you’ve got variables read onto secondary hydrology meshes (say, a basal DEM or a runoff raster) that need to be transferred to the main hydrology mesh
* `Number of Variables To Read = Integer n`
  * If you are loading read-in variables, say how many there are. Note: the solver is currently only set up to deal with a maximum of 10 variables; if you have more than that, you’ll need to modify the source code
* `Reader Solver 1 = Integer n` and `Reader V1 = String “name”`
  * If loading variables, say which solver number the reader solver is and what the name of the variable is
  * You should provide as many `Reader Solver` and `Reader V` entries as the number you’ve defined under `Number Of Variables To Read`, numbered sequentially
* `Reference Node(3) = Integer x y z` and `Threshold Distance = Real n`
  * These are used as part of the artefact correction for the grounded mask – if specified, all nodes on the hydrology mesh greater than the Threshold Distance from the Reference Node will be automatically set to grounded. This can be useful if dealing with artefacts a long way inland
* `Side = Logical True/False`
  * This isn’t a solver option, but something that should be set in the boundary condition section of the SIF. The interpolation routine tends to create a lot of artefacts on the lateral boundaries of the domain – if you set `Side = Logical True` in the boundary condition sections for the sidewalls of the hydrology mesh, this will force them to be grounded and remove the frequent ungrounded artefacts.
Known issues:
* None, though the artefact correction could probably be improved

### HydroRestart.F90
Largely a direct copy of the Restart() subroutine within the Elmer source code, with a few modifications to disentangle it from all the other restart machinery and to make it pick up the right variables from the right place and send them to the right place.
Solver options and inputs:
* `hp: Restart Variable 1 = String “Name”`
  * All the variables that should be restarted on the primary hydrology mesh (i.e. all those that are normally calculated by or associated with GlaDSCoupledSolver) should be listed like this – just number the entries consecutively
* `channel: Restart Variable 1 = String “Name”`
  * Same as above, but for the channel variables (i.e. channel area and channel flux)
* `sheet: Restart Variable 1 = String “Name”`
  * Same as above, but for the sheet variables (i.e. sheet thickness – it’s the only one associated with the separate sheet thickness solver)
* The changes to GlaDSCoupledSolver.F90 will then ensure all the variables end up associated with the primary hydrology mesh, but it’s best to restart them on their own solver meshes
Known issues:
* None

### USF_SourceCalcCalving.F90
This is a USF that calculates the Hydraulic Potential Volume Source term required in the Body Force section of the SIF for GlaDS. If provided with a surface runoff variable (I usually load it in from a raster) and the temperature residual variable, it will calculate the resulting internal melt and add on the surface melt for each node on the hydrology mesh (or, at the base of your 3D ice mesh, if you’re using GlaDS without all the other bells and whistles), so you can easily vary the source term spatially across your domain.
USF options and inputs:
* These all go in the same Body Force section as where you define the source term
* `Internal Melt = Logical True/False`
  * Switch for whether you want to work out internal melt or not
* `Surface Melt = Logical True/False`
  * Same for surface melt
* `Internal Melt Variable Name = String “Name”`
  * If you are using internal melt, give the name of the variable that you want it worked out from (note: the USF is set up on the assumption that this will be the temperature residual from the TemperateIceSolver. If you want to use something else, you’ll need to change the code in the USF)
* `Surface Melt Variable Name = String “Name”`
  * If using surface melt, the name of the variable that contains it
* Finally, when defining the source term, you call this USF just like any other, and it does not matter what variable you use in the call – the USF will ignore it.

### BasalMelt3D.F90
This is a very simple solver written by Joe Todd that applies a specified basal melt rate to any ungrounded parts of the glacier base.
Solver options and inputs:
* `Basal Melt Stats File = String …`
  * The path to write a file containing some basal melt stats to.
* `GroundedMask Variable = String …`
  * The name of the variable used by the grounded solver
* `Basal Melt Mode = String …`
  * Can be ‘seasonal’ or ‘off’. The latter is fairly self-explanatory; the former requires the following additional options:
* `Basal Melt Summer Rate = Real …`
  * The rate to apply basal melt at in summer.
* `Basal Melt Winter Rate = Real …`
  * The rate to apply basal melt at in winter.
* `Basal Melt Summer Start = Real …`
  * The time in the simulation to begin using summer melt rates (expressed as a number between 0 and 1)
* `Basal Melt Summer Stop = Real …`
  * The time in the simulation to begin using summer melt rates (expressed as a number between 0 and 1)

### GMValid.F90
This is a very simple solver that’s more-or-less just a stripped-down copy of BasalMelt3D.F90 and exists to set up a mask variable for which ungrounded areas are connected to the fjord and which aren’t.
Solver options and inputs:
* None – just the usual lines to define the equation, etc.

## Problems
If you find something that doesn’t work or you can’t easily fix or isn’t listed here, email Samuel Cook (samuel.cook@univ-grenoble-alpes.fr)

## Known Bugs
* There is an unresolved issue somewhere in the calving code that means, sometimes, the coupled model will just crash or hang at the end of the SwitchMesh routine for a reason I haven’t yet pinned down – it’s something to do with particularly complicated or large calving events (I think it's down to calving producing strange meshes). If this happens, rerunning the simulation will usually fix things (because the calving code has an element of randomness, the same event won’t recur in exactly the same way). However, be prepared to give it a few goes (5-6) if it keeps happening. If you’re running long simulations, the model **WILL** crash on this at some point, so the best plan is to just restart from the last full result output timestep the run completed and carry on from there, rather than trying to get a complete run in one go. If the model just keeps crashing at the same point, consider running the previous simulation you’ve restarted from – I found that helps.

## Other Modifications
There have also been other changes to the source code that are not of particular relevance to the end user, but are included here for completeness:
* New MeshTag field defined as a new integer field on the Mesh_t variable type in **Types.F90**. This integer is incremented by **CalvingRemesh.F90** every time remeshing occurs and is then used by all solvers and USFs to decide whether memory allocation needs to be redone, rather than relying on the Mesh % Changed Boolean, which is unreliable
* Line adding CalvingTime as an entry to the Model % Simulation list in **Calving3D.F90**, which is then used by **CalvingRemesh.F90** as described above
