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
 
int InspectElement(struct FemType *data,int idx);

int LoadSolutionElmer(struct FemType *data,int results,char *prefix,int info);
int SaveSolutionElmer(struct FemType *data,struct BoundaryType *bound,
		      int nobound,char *prefix,int decimals,int info);
int SaveSizeInfo(struct FemType *data,struct BoundaryType *bound,
		 char *prefix,int info);
int SaveElmerInputFemBem(struct FemType *data,struct BoundaryType *bound,
			 char *prefix,int decimals, int info);

void InspectVector(Real *vector,int first,int last,Real *min,
		   Real *max,int *mini,int *maxi);
int  Steepest(Real *vector,int first,int last);
Real MeanVector(Real *vector,int first,int last);
Real AbsMeanVector(Real *vector,int first,int last);
Real DifferVector(Real *vector1,Real *vector2,int first,int last);
void ReformVector(Real *vector1,int n1,Real *vector2,int n2);
void AdjustVector(Real max,Real min,Real *vector,int first,int last);
int  ReadRealVector(Real *vector,int first,int last,char *filename);
void SaveRealVector(Real *vector,int first,int last,char *filename);
int  ReadRealMatrix(Real **matrix,int row_first,int row_last,
		    int col_first,int col_last,char *filename);
void SaveRealMatrix(Real **matrix,int row_first,int row_last,
		    int col_first,int col_last,char *filename);
int  ReadIntegerVector(int *vector,int first,int last,char *filename);
void SaveIntegerVector(int *vector,int first,int last,char *filename);
int  ReadIntegerMatrix(int **matrix,int row_first,int row_last,
		       int col_first,int col_last,char *filename);
void SaveIntegerMatrix(int **matrix,int row_first,int row_last,
		       int col_first,int col_last,char *filename);
void SaveNonZeros(Real **matrix,int row_first,int row_last,
		  int col_first,int col_last,char *filename);
int EchoFile(char *filename);
int SetDiscontinuousPoints(struct FemType *data,struct PointType *point,
			   int material);
