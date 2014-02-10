/* femfilein.h */
/* Routines for importing existing FEM meshes */


/* These have been removed from the ElmerGUI version
   int LoadElmerInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
   int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info); */

int LoadAbaqusInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadFidapInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadAnsysInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadNastranInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadFieldviewInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadTriangleInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadMeditInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadComsolMesh(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadGidInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadGmshInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadUniversalMesh(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadCGsimMesh(struct FemType *data,char *prefix,int info);
