# Ronne-Filchner these case
  
Optimisation of the viscosity and friction for the Ronne-Filchned ice-shelf including ice rises.

To run this experiment some data sets are required (see instructions in the DATA directory):

  -  observations of the surface velocities
  -  Bedmachine for the topography
  

The model domain is defined in the shapefile domain/bassin.shp and contains a contour of the Ronne Filchner ice shelf
that also includes the ice rises (Berkner Island).

## Mesh generation:

Use the python code [Contour2geo.py](https://github.com/ElmerCSC/elmerfem/tree/elmerice/elmerice/Meshers/GIS) 
to generate the Gmsh .geo file from the shapefile.

```bash
python $ELMER_SRC/elmerice/Meshers/GIS/Contour2geo.py -r 5000 -i domain/bassin.shp -o mesh.geo --spline
```

Convert to Elmer format:

```bash
gmsh -2 -format msh2 mesh.geo
ElmerGrid 14 2 mesh.msh -autoclean
# optional partition the mesh
ElmerGrid 2 2 mesh -metis 2
```

## Running the examples

Required configuration files (.sif) are under the SIF directory.

- INIT.sif : initialise the geometry (thickness, bedrock, top and bottom surfaces, groundedMask, and observed velocities)  
Requires the [Scattered2DDataInterpolator solver](http://elmerfem.org/elmerice/wiki/doku.php?id=solvers:scattered). Alternatively the [Griddatareader](http://elmerfem.org/elmerice/wiki/doku.php?id=solvers:griddatareader) can be used.

- OPTIM.sif : run the optimisation of both viscosity and friction  

- Physical_Params.IN: contains the definition of some physical parameters   

- Makefile: compile the USFs_RonneFilchner user function   
