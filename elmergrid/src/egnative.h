/* femelmer.h -> egnative.h */
int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info);
int LoadElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int nonames, int info);
int SaveSolutionElmer(struct FemType *data,struct BoundaryType *bound,
		      int nobound,char *prefix,int decimals,int info);
int SaveSolutionElmerTriangles(struct FemType *data,char *prefix,int info);
int SaveElmerInput(struct FemType *data,struct BoundaryType *bound,
		   char *prefix,int decimals,int nooverwrite, int info);
int SaveSizeInfo(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info);
int SaveElmerInputFemBem(struct FemType *data,struct BoundaryType *bound,
			 char *prefix,int decimals, int info);
