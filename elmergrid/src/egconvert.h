/* femfilein.h -> egconvert.h */
/* Routines for importing meshes and data from other formats. */

/* These routines are C. ElmerGUI compiles the ElmerGrid core as a library and
   calls into it from C++, so the declarations need C linkage there. */
#ifdef __cplusplus
extern "C" {
#endif

int LoadAbaqusInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadAbaqusOutput(struct FemType *data,char *prefix,int info);
int LoadFidapInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadAnsysInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadNastranInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadFieldviewInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadTriangleInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadMeditInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadComsolMesh(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadGidInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadGmshInput(struct FemType *data,struct BoundaryType *bound,
		  char *prefix,int keeporphans,int dim,int info);
int LoadGeoInput(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadFvcomMesh(struct FemType *data,struct BoundaryType *bound,char *filename,int info);
int LoadUniversalMesh(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadCGsimMesh(struct FemType *data,char *prefix,int info);
int LoadFluxMesh(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);
int LoadFluxMesh3D(struct FemType *data,struct BoundaryType *bound,char *prefix,int info);

#ifdef __cplusplus
}
#endif
