/* feminfo.h -> egextra.h */
/* Functions providing the user additional information on mesh creation,
   solution and calculation. */

int SaveSolutionDens(struct FemType *data,char *prefix,int info);
int SaveCellInfo(struct GridType *grid,struct CellType *cell,
		 char *prefix,int info);
int SaveBoundary(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info);
int CreateBoundaryChain(struct FemType *data,struct BoundaryType *bound,int info);
int SaveBoundariesChain(struct FemType *data,struct BoundaryType *bound,
			char *prefix,int info);
int SaveBoundaryLine(struct FemType *data,int direction,
		     Real c0,char* prefix,int info);
int SaveBoundaryForm(struct FemType *data,struct CellType *cell, 
		     char* filename,int info);
int SaveSubcellForm(struct FemType *data,struct CellType *cell, 
		    char* filename,int info);
 
int ShowCorners(struct FemType *knot,int variable,Real offset);
int InspectElement(struct FemType *data,int idx);

int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info);
int SaveSolutionElmer(struct FemType *data,struct BoundaryType *bound,
		      int nobound,char *prefix,int decimals,int info);
int SaveSizeInfo(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info);
int SaveElmerInputFemBem(struct FemType *data,struct BoundaryType *bound,
			 char *prefix,int decimals, int info);
int SolutionFromMeshToMesh(struct CellType *cell1, struct GridType *grid1, 
			   struct FemType *data1,
			   struct CellType *cell2, struct GridType *grid2, 
			   struct FemType *data2,
			   int mapgeo,int variable,int info);

