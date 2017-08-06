## make the mesh
ElmerGrid 1 2 rectangle


#UpdateExport.f90: is a solver used to update auxiliary variables
#declared as "Exported Variable" in the solver section; accorded to their definition in the body Force.
#E.G. Here surface Elevation is a variable computed as Zs=Zb+H where H is solution of the thickness evolution equation.

## run the test

## ismip-hom test in 1D
ElmerSolver ismip_SSA_1D.sif

## ismip-hom test in 2D
#ElmerSolver ismip_SSA_2D.sif

## ismip-hom test in 2D but with a 3D mesh
#ElmerSolver ismip_SSA_3D.sif

## Grow an ice cap
#ElmerSolver SSA_IceSheet.sif 

## Test the Weertman friction law with SSA
#ElmerSolver ismip_SSA_2D_Weertman.sif 

## Test the Coulomb friction law with SSA
#ElmerSolver ismip_SSA_2D_Coulomb.sif 
