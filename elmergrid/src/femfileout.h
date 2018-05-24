/* femfileout.h */
/* Routines for exporting data to other
   FEM packages, such as Abaqus and Fidap. */

int SaveAbaqusInput(struct FemType *data,char *prefix,int info);
int SaveFidapOutput(struct FemType *data,char *prefix,int info,
		    int vctrs,Real *vect1, ...);
int SaveMeshGmsh(struct FemType *data,struct BoundaryType *bound,
		 int nobound,char *prefix,int decimals,int info);
int SaveMeshVtu(struct FemType *data,struct BoundaryType *bound,
		 int nobound,char *prefix,int dummyzero,int info);
